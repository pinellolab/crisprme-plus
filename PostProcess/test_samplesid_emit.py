#!/usr/bin/env python3
"""Regression tests for the samplesID SELF-COMPLETENESS fix (v2.3.2).

Problem being fixed: neither ``build-index-only`` nor ``publish_index`` used to
write/upload any samplesID file, so a user who BUILDS + PUBLISHES their OWN merged
variant index was not self-complete -- the combined ``<vcf>.samplesID.txt`` the
search's ``--samplesID`` needs was produced ONLY download-side, unioning per-db
lists that had to already be on HF. These tests pin the two halves of the fix:

  (a) BUILD (crisprme.py._emit_combined_samplesid): a MERGED-panel build emits
      ``samplesIDs/<vcf>.samplesID.txt`` (union of the per-db lists) INTO the
      install, and normalizes the per-db source files to the canonical
      ``<ref>_<db>.samplesID.txt`` names the shared synthesizer resolves by (even
      when the --samplesID inputs were arbitrarily named).
  (b) PUBLISH + DOWNLOAD (crisprme_hf): publish BUNDLES the samplesID lists into
      the main index tarball (arcname samplesIDs/<name>) and sets
      manifest['has_samplesids']; download installs them into <workdir>/samplesIDs/
      BEFORE the 2.3.1 fallback runs (which then no-ops) -> `download --what index`
      alone yields the combined samplesID the search needs, no `--what samples`.
      Never clobbers an already-present samplesID.
  (c) NO-OP: single-dataset index (no combined file needed) + reference-only index
      (no samplesID block at all) publish byte-for-byte as before.

There is ONE union implementation (crisprme_hf.synthesize_combined_samplesid); the
build reuses it. Fully offline: crisprme_hf imports huggingface_hub lazily and the
HF client is stubbed. ``crisprme.py`` is loaded via importlib with ``Bio.Seq``
stubbed (its only unavailable top-level import in the light CI), so the build-side
helper is exercised as-shipped, not re-implemented.
"""
import ast
import errno
import io
import os
import shutil
import sys
import tarfile
import tempfile
import threading
import types
import unittest
from unittest import mock

import crisprme_hf as hf


# --------------------------------------------------------------------------- #
# Load crisprme.py (for _emit_combined_samplesid / _build_db_to_samplesid) with
# Bio.Seq stubbed -- the ONLY module crisprme.py imports at top level that the
# network-free unit-test env lacks (everything else is PostProcess-local or stdlib).
# --------------------------------------------------------------------------- #
def _load_crisprme():
    """Exec crisprme.py's DEFINITIONS (with Bio.Seq stubbed) but NOT its top-level
    CLI dispatch. crisprme.py has no ``if __name__ == '__main__'`` guard -- it runs
    ``if len(sys.argv) < 2 ... elif sys.argv[1] == ...`` at module level -- so we
    truncate the source at that dispatch block and exec only the part above it. That
    gives us the real, as-shipped ``_emit_combined_samplesid`` / ``_build_db_to_samplesid``
    functions without launching the CLI or its side effects (check_crisprme_dirtree)."""
    here = os.path.dirname(os.path.abspath(__file__))
    crisprme_py = os.path.join(os.path.dirname(here), "crisprme.py")
    # crisprme.py inserts its own PostProcess dir on sys.path; make sure this one
    # (where crisprme_hf etc. live) resolves first regardless.
    sys.path.insert(0, here)
    with open(crisprme_py) as fh:
        src = fh.read()
    marker = "\nif len(sys.argv) < 2:"
    cut = src.index(marker)  # KeyError-loud if the dispatch shape ever changes
    defs_only = src[:cut]
    mod = types.ModuleType("crisprme_under_test")
    mod.__file__ = crisprme_py
    code = compile(defs_only, crisprme_py, "exec")
    # crisprme.py does `from Bio.Seq import Seq` -- a SUBMODULE import that a MagicMock
    # parent can't satisfy ("'Bio' is not a package"), so register Bio + Bio.Seq in
    # sys.modules explicitly (Seq only needs to be an importable name). The functions
    # under test don't use Bio.
    if "Bio.Seq" not in sys.modules:
        _bio = types.ModuleType("Bio")
        _bio_seq = types.ModuleType("Bio.Seq")
        _bio_seq.Seq = object
        _bio.Seq = _bio_seq
        sys.modules.setdefault("Bio", _bio)
        sys.modules["Bio.Seq"] = _bio_seq
    # The network-free unit-test env also lacks pandas (imported by assembly_reconcile,
    # which even uses pd.DataFrame in def-time annotations) and possibly other pure-
    # scientific modules. Stub ONLY those, TARGETED -- never a catch-all: a catch-all
    # would also stub optional-import-fallback targets (e.g. requests' `try: import
    # simplejson`), handing them a MagicMock instead of the expected absence and breaking
    # them (MagicMock.JSONDecodeError as an exception base -> metaclass conflict). Use a
    # MagicMock so attribute access (pd.DataFrame) works. requests/numpy are installed on
    # CI; a module that is actually present is imported for real (we stub only absents).
    from unittest.mock import MagicMock
    for _m in ("pandas", "scipy", "sklearn", "matplotlib", "seaborn",
               "statsmodels", "intervaltree", "CRISPR_BULGE_score"):
        if _m in sys.modules:
            continue
        try:
            __import__(_m)
        except Exception:
            sys.modules[_m] = MagicMock()
    exec(code, mod.__dict__)
    return mod


crisprme = _load_crisprme()


def _mkfile(path, content="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(content)


# per-db samplesID DATA rows (header-less, tab-separated: SAMPLE_ID<TAB>...)
_1000G_ROWS = "S_EUR_1\tGBR\tEUR\tmale\nS_EUR_2\tFIN\tEUR\tfemale\nSHARED\tX\tX\tX\n"
_HGDP_ROWS = "S_AFR_1\tYRI\tAFR\tmale\nSHARED\tX\tX\tX\n"  # SHARED is a dup vs 1000G


# --------------------------------------------------------------------------- #
# (a) BUILD-side emission
# --------------------------------------------------------------------------- #
class TestBuildEmitCombinedSamplesID(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="sid_build_")
        self.sdir = os.path.join(self.workdir, "samplesIDs")
        os.makedirs(self.sdir, exist_ok=True)
        self.ref = "hg38"
        self.vcf = "hg38_1000G_HGDP"

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def test_merged_build_emits_combined_and_canonical_perdb(self):
        """A merged build (>=2 datasets) writes the combined union AND normalizes
        the per-db files to canonical <ref>_<db>.samplesID.txt names -- EVEN when
        the --samplesID inputs were arbitrarily named (the load-bearing step that
        lets the shared synthesizer's convention-based lookup find them)."""
        # arbitrarily-named per-db sources (as a build listing may reference)
        src_1000g = os.path.join(self.sdir, "samplesIDs.1000G.txt")
        src_hgdp = os.path.join(self.sdir, "hgdp_panel.samplesID.txt")
        _mkfile(src_1000g, _1000G_ROWS)
        _mkfile(src_hgdp, _HGDP_ROWS)
        # ORDERED {db_label: path} exactly as crisprme._build_db_to_samplesid yields
        db_map = {"1000G": src_1000g, "HGDP": src_hgdp}

        written = crisprme._emit_combined_samplesid(
            self.workdir, self.vcf, self.ref, db_map
        )

        combined = os.path.join(self.sdir, f"{self.vcf}.samplesID.txt")
        self.assertEqual(written, combined, "combined path not returned/written")
        self.assertTrue(os.path.isfile(combined))
        # per-db files copied to CANONICAL names so the synthesizer resolves them
        self.assertTrue(os.path.isfile(os.path.join(self.sdir, "hg38_1000G.samplesID.txt")))
        self.assertTrue(os.path.isfile(os.path.join(self.sdir, "hg38_HGDP.samplesID.txt")))
        # union: header-less, dedup by SAMPLE_ID (SHARED once), stable order (1000G first)
        ids = [ln.split("\t", 1)[0] for ln in open(combined) if ln.strip()]
        self.assertEqual(ids, ["S_EUR_1", "S_EUR_2", "SHARED", "S_AFR_1"],
                         f"combined union wrong / not deduped / wrong order: {ids}")

    def test_single_dataset_build_is_noop(self):
        """(c) A single-dataset map (<2 entries) needs no combined file -> NO-OP:
        the synthesizer would no-op anyway, and no combined file is written."""
        src = os.path.join(self.sdir, "hg38_1000G.samplesID.txt")
        _mkfile(src, _1000G_ROWS)
        written = crisprme._emit_combined_samplesid(
            self.workdir, "hg38_1000G", self.ref, {"1000G": src}
        )
        self.assertIsNone(written, "single-dataset build must not synthesize a combined file")
        self.assertFalse(os.path.isfile(os.path.join(self.sdir, "hg38_1000G_HGDP.samplesID.txt")))

    def test_empty_map_is_noop(self):
        """(c) No --samplesID (empty map, e.g. reference-only reached this path) ->
        NO-OP, nothing written."""
        self.assertIsNone(crisprme._emit_combined_samplesid(self.workdir, self.vcf, self.ref, {}))
        self.assertFalse(os.listdir(self.sdir))

    def test_existing_combined_not_clobbered(self):
        """A shipped combined file is never overwritten (synthesizer's first guard)."""
        _mkfile(os.path.join(self.sdir, "hg38_1000G.samplesID.txt"), _1000G_ROWS)
        _mkfile(os.path.join(self.sdir, "hg38_HGDP.samplesID.txt"), _HGDP_ROWS)
        combined = os.path.join(self.sdir, f"{self.vcf}.samplesID.txt")
        _mkfile(combined, "PRESHIPPED\tX\tX\tX\n")
        written = crisprme._emit_combined_samplesid(
            self.workdir, self.vcf, self.ref,
            {"1000G": os.path.join(self.sdir, "hg38_1000G.samplesID.txt"),
             "HGDP": os.path.join(self.sdir, "hg38_HGDP.samplesID.txt")},
        )
        self.assertIsNone(written, "must not report a write when the file already exists")
        self.assertEqual(open(combined).read(), "PRESHIPPED\tX\tX\tX\n",
                         "shipped combined samplesID was clobbered")


class _FullDisk:
    """A file being written that takes half of what it is given, then fails, as a
    file on a disk that has just filled up does."""

    def __init__(self, fh):
        self.fh = fh

    def write(self, text):
        self.fh.write(text[: len(text) // 2])
        self.fh.flush()
        raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.fh.close()


def _load_ensure_samplesid(workdir):
    """pages/main_page._ensure_samplesid, as shipped, with the web app's working
    directory workdir (the page module needs Dash; this needs only its source)."""
    page = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "pages", "main_page.py")
    src = open(page).read()
    node = next(n for n in ast.parse(src).body
                if isinstance(n, ast.FunctionDef) and n.name == "_ensure_samplesid")
    ns = {"os": os, "threading": threading, "current_working_directory": workdir}
    exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
    return ns["_ensure_samplesid"]


class TestSamplesIDWrittenWholeOrNotAtAll(unittest.TestCase):
    """build-index-only (each per-dataset copy, then the union), the download
    (synthesize_combined_samplesid) and a web search (pages/main_page
    ._ensure_samplesid) wrote each samplesID file in place, and every later run skips
    one that is there: a full disk or a kill as it was written left it cut short, kept
    for good and searched with. Each is now written beside it, then renamed into
    place, and a full disk fails the build."""

    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="sid_whole_")
        self.sdir = os.path.join(self.workdir, "samplesIDs")
        self.src = os.path.join(self.workdir, "lists")
        _mkfile(os.path.join(self.src, "1000G.txt"), _1000G_ROWS)
        _mkfile(os.path.join(self.src, "HGDP.txt"), _HGDP_ROWS)
        self.map = {"1000G": os.path.join(self.src, "1000G.txt"),
                    "HGDP": os.path.join(self.src, "HGDP.txt")}
        self.union = _1000G_ROWS + "S_AFR_1\tYRI\tAFR\tmale\n"
        self.combined = os.path.join(self.sdir, "hg38_1000G_HGDP.samplesID.txt")

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def full_disk(self):
        """Every file opened for writing in samplesIDs/, and every copy into it,
        takes half of its data, then fails as on a full disk."""
        real_open, real_copy, sdir = open, shutil.copy2, self.sdir

        def fake_open(path, mode="r", *args, **kwargs):
            fh = real_open(path, mode, *args, **kwargs)
            inside = os.path.dirname(os.path.abspath(str(path))) == sdir
            return _FullDisk(fh) if inside and "w" in mode else fh

        def fake_copy(src, dst, *args, **kwargs):
            with real_open(src) as fh, real_open(dst, "w") as out:
                _FullDisk(out).write(fh.read())

        return mock.patch("builtins.open", fake_open), mock.patch("shutil.copy2", fake_copy)

    def written(self):
        """What samplesIDs/ holds: {name: text}, hidden files too."""
        out = {}
        for name in sorted(os.listdir(self.sdir)) if os.path.isdir(self.sdir) else ():
            with open(os.path.join(self.sdir, name)) as fh:
                out[name] = fh.read()
        return out

    def test_a_full_disk_leaves_no_file_of_the_build_cut_short(self):
        for which in ("copy", "union"):
            with self.subTest(which=which):
                shutil.rmtree(self.sdir, ignore_errors=True)
                if which == "union":  # (the copies there already)
                    _mkfile(os.path.join(self.sdir, "hg38_1000G.samplesID.txt"), _1000G_ROWS)
                    _mkfile(os.path.join(self.sdir, "hg38_HGDP.samplesID.txt"), _HGDP_ROWS)
                before = self.written()
                opened, copied = self.full_disk()
                with opened, copied, mock.patch("sys.stderr", new_callable=io.StringIO) as err:
                    with self.assertRaises(SystemExit) as end:
                        crisprme._emit_combined_samplesid(self.workdir, "hg38_1000G_HGDP",
                                                          "hg38", self.map)
                self.assertEqual(end.exception.code, 1)
                self.assertTrue(err.getvalue().startswith(
                    f"Error: Writing the samplesID files of hg38_1000G_HGDP into {self.sdir} "
                    f"failed ([Errno {errno.ENOSPC}] {os.strerror(errno.ENOSPC)}"), err.getvalue())
                self.assertTrue(err.getvalue().endswith(
                    "): run build-index-only again once that is fixed, which writes them.\n"))
                self.assertEqual(self.written(), before)  # nothing cut short, nothing left
                # once there is room, the next run writes them whole
                self.assertEqual(crisprme._emit_combined_samplesid(
                    self.workdir, "hg38_1000G_HGDP", "hg38", self.map), self.combined)
                self.assertEqual(self.written(), {"hg38_1000G.samplesID.txt": _1000G_ROWS,
                                                  "hg38_HGDP.samplesID.txt": _HGDP_ROWS,
                                                  "hg38_1000G_HGDP.samplesID.txt": self.union})

    def test_a_full_disk_leaves_no_union_of_the_download_or_the_web_cut_short(self):
        ensure = _load_ensure_samplesid(self.workdir)
        for who, write in (
                ("download", lambda: hf.synthesize_combined_samplesid(
                    self.workdir, "hg38_1000G_HGDP", ref="hg38")),
                ("web search", lambda: ensure("hg38", "hg38_1000G_HGDP"))):
            with self.subTest(who=who):
                shutil.rmtree(self.sdir, ignore_errors=True)
                _mkfile(os.path.join(self.sdir, "hg38_1000G.samplesID.txt"), _1000G_ROWS)
                _mkfile(os.path.join(self.sdir, "hg38_HGDP.samplesID.txt"), _HGDP_ROWS)
                before = self.written()
                opened, _ = self.full_disk()
                with opened:
                    if who == "download":  # (it says so, and goes on)
                        with self.assertRaises(OSError):
                            write()
                    else:  # (best effort: the search then fails for want of it)
                        write()
                self.assertEqual(self.written(), before)
                write()
                self.assertEqual(self.written(), dict(before, **{
                    "hg38_1000G_HGDP.samplesID.txt": self.union}))


# --------------------------------------------------------------------------- #
# (b) PUBLISH + DOWNLOAD round-trip (self-complete `--what index`)
# --------------------------------------------------------------------------- #
class _FakeApi:
    def __init__(self, sink_dir):
        self.sink_dir = sink_dir
        self.uploads = []

    def upload_file(self, path_or_fileobj, path_in_repo, repo_id, repo_type, token):
        dst = os.path.join(self.sink_dir, os.path.basename(path_in_repo))
        shutil.copy(path_or_fileobj, dst)
        self.uploads.append((dst, path_in_repo))


class _FakeHF:
    def __init__(self, api):
        self._api = api

    def HfApi(self):
        return self._api


def _mk_variant_layout(root, index_name, ref, with_samplesids=True, perdb_only=False):
    """Build a minimal variant index + Dictionaries siblings + (optionally) the
    samplesIDs/ that _emit_combined_samplesid would have produced at build time."""
    vcf = index_name.partition("+")[2]
    idx = os.path.join(root, "genome_library", index_name)
    ind = idx + "_INDELS"
    _mkfile(os.path.join(idx, "TSTgenome.NGG.bin"))
    _mkfile(os.path.join(ind, "NGG_2_fakechr22"))
    droot = os.path.join(root, "Dictionaries")
    _mkfile(os.path.join(droot, f"log_indels_{vcf}", "logchr22.txt"), "chr22\n")
    _mkfile(os.path.join(droot, f"dictionaries_{vcf}", "my_dict_chr22.json"), "{}")
    if with_samplesids:
        sdir = os.path.join(root, "samplesIDs")
        _mkfile(os.path.join(sdir, f"{ref}_1000G.samplesID.txt"), _1000G_ROWS)
        if not perdb_only:
            _mkfile(os.path.join(sdir, f"{ref}_HGDP.samplesID.txt"), _HGDP_ROWS)
            # combined (as the build fix emits)
            _mkfile(os.path.join(sdir, f"{vcf}.samplesID.txt"),
                    _1000G_ROWS + "S_AFR_1\tYRI\tAFR\tmale\n")
    return idx


def _publish(index_dir, sink_dir, dictless=False):
    api = _FakeApi(sink_dir)
    orig_require, orig_token = hf._require_hf, hf.resolve_token
    hf._require_hf = lambda: _FakeHF(api)
    hf.resolve_token = lambda t=None: "faketoken"
    try:
        hf.publish_index(index_dir, repo="local/test", token="x", dictless=dictless)
    finally:
        hf._require_hf, hf.resolve_token = orig_require, orig_token
    return api


class TestPublishBundlesSamplesID(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="sid_pub_")
        self.sink = tempfile.mkdtemp(prefix="sid_pubhf_")
        self.ref = "hg38"
        self.name = "NGG_2_hg38+hg38_1000G_HGDP"
        self.vcf = self.name.partition("+")[2]

    def tearDown(self):
        for d in (self.root, self.sink):
            shutil.rmtree(d, ignore_errors=True)

    def _main_tarball(self, api):
        return next(p for p, r in api.uploads
                    if r.endswith(f"{self.name}.tar.gz") and "genotypes_" not in r)

    def test_publish_bundles_combined_and_perdb(self):
        """(b) publish routes the combined + per-db samplesID lists into the main
        tarball under samplesIDs/<name>, and sets manifest['has_samplesids']."""
        idx = _mk_variant_layout(self.root, self.name, self.ref)
        api = _publish(idx, self.sink)
        members = tarfile.open(self._main_tarball(api)).getnames()
        self.assertIn(f"samplesIDs/{self.vcf}.samplesID.txt", members,
                      f"combined samplesID not bundled: {members}")
        self.assertIn(f"samplesIDs/{self.ref}_1000G.samplesID.txt", members)
        self.assertIn(f"samplesIDs/{self.ref}_HGDP.samplesID.txt", members)
        # manifest flag
        import json
        mbytes = tarfile.open(self._main_tarball(api)).extractfile("manifest.json").read()
        self.assertTrue(json.loads(mbytes).get("has_samplesids"))

    def test_publish_dictless_ref_marker_canonicalized(self):
        """(risk) a '-dictless'-marked index name must still bundle per-db files
        named <ref>_<db> (NOT <ref>-dictless_<db>) so download (--ref hg38) finds
        them. manifest['genome'] is 'hg38-dictless' here; the bundling must strip it."""
        marked = "NGG_2_hg38-dictless+hg38_1000G_HGDP"
        idx = _mk_variant_layout(self.root, marked, self.ref)
        api = _publish(idx, self.sink, dictless=True)
        main = next(p for p, r in api.uploads
                    if r.endswith(f"{marked}.tar.gz") and "genotypes_" not in r)
        members = tarfile.open(main).getnames()
        self.assertIn(f"samplesIDs/{self.ref}_1000G.samplesID.txt", members,
                      f"per-db samplesID must use canonicalized ref token: {members}")
        self.assertIn(f"samplesIDs/{self.ref}_HGDP.samplesID.txt", members)
        self.assertFalse(any("hg38-dictless_" in m for m in members),
                         f"per-db files wrongly named with the -dictless marker: {members}")

    def test_reference_only_index_no_samplesid_block(self):
        """(c) a reference-only index (no '+') publishes with NO samplesIDs/ members
        and has_samplesids=false -> byte-for-byte-unchanged behavior."""
        name = "NGG_2_hg38"
        idx = os.path.join(self.root, "genome_library", name)
        _mkfile(os.path.join(idx, "TSTgenome.NGG.bin"))
        api = _publish(idx, self.sink)
        members = tarfile.open(api.uploads[0][0]).getnames()
        self.assertFalse(any(m.startswith("samplesIDs/") for m in members),
                         f"reference-only index must bundle no samplesIDs: {members}")
        import json
        mbytes = tarfile.open(api.uploads[0][0]).extractfile("manifest.json").read()
        self.assertFalse(json.loads(mbytes).get("has_samplesids"))

    def test_legacy_variant_index_no_samplesid_files_bundles_none(self):
        """(c) a legacy variant index built before the emit fix (no samplesIDs/ on
        disk) bundles zero samplesID files -> has_samplesids=false; download-side
        2.3.1 fallback still covers it."""
        idx = _mk_variant_layout(self.root, self.name, self.ref, with_samplesids=False)
        api = _publish(idx, self.sink)
        members = tarfile.open(self._main_tarball(api)).getnames()
        self.assertFalse(any(m.startswith("samplesIDs/") for m in members),
                         f"legacy variant index must bundle no samplesIDs: {members}")


class TestDownloadInstallsBundledSamplesID(unittest.TestCase):
    """Full publish->download round-trip: a fresh `download --what index` yields the
    combined samplesID the search needs, with NO `--what samples` and NO HF per-db
    fetch (the 2.3.1 fallback strictly no-ops)."""

    def setUp(self):
        self.src = tempfile.mkdtemp(prefix="sid_dlsrc_")
        self.sink = tempfile.mkdtemp(prefix="sid_dlhf_")
        self.workdir = tempfile.mkdtemp(prefix="sid_dldst_")
        self.ref = "hg38"
        self.name = "NGG_2_hg38+hg38_1000G_HGDP"
        self.vcf = self.name.partition("+")[2]
        self._orig_snapshot = hf._hf_snapshot

    def tearDown(self):
        hf._hf_snapshot = self._orig_snapshot
        for d in (self.src, self.sink, self.workdir):
            shutil.rmtree(d, ignore_errors=True)

    def _serve_from_sink(self, fetch_log):
        """snapshot stub: serve indexes/* from the sink; record every requested
        pattern so we can PROVE no per-db samplesID HF fetch happened."""
        def fake_snapshot(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
            fetch_log.extend(allow_patterns)
            for pat in allow_patterns:
                base = os.path.basename(pat)
                if pat.startswith("indexes/"):
                    dest = os.path.join(local_dir, "indexes")
                    os.makedirs(dest, exist_ok=True)
                    srcf = os.path.join(self.sink, base)
                    if os.path.isfile(srcf):
                        shutil.copy(srcf, os.path.join(dest, base))
                else:
                    # a samplesIDs/ fetch would land here; we DON'T serve it, to
                    # prove the bundled files made this fetch unnecessary (no-op).
                    pass
            return local_dir
        hf._hf_snapshot = fake_snapshot

    def test_self_complete_download_no_samples_component(self):
        """(b) publish (with bundled samplesIDs) -> download --what index installs the
        combined + per-db lists into <workdir>/samplesIDs/, and the 2.3.1 per-db HF
        fetch never fires (bundled files short-circuit it)."""
        idx = _mk_variant_layout(self.src, self.name, self.ref)
        _publish(idx, self.sink)
        fetch_log = []
        self._serve_from_sink(fetch_log)
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, ref=self.ref, token="x",
                              genotypes=False)
        sdir = os.path.join(self.workdir, "samplesIDs")
        combined = os.path.join(sdir, f"{self.vcf}.samplesID.txt")
        self.assertTrue(os.path.isfile(combined),
                        "combined samplesID not installed from the bundled tarball")
        self.assertTrue(os.path.isfile(os.path.join(sdir, f"{self.ref}_1000G.samplesID.txt")))
        self.assertTrue(os.path.isfile(os.path.join(sdir, f"{self.ref}_HGDP.samplesID.txt")))
        # PROOF the download-side fallback no-op'd: no samplesIDs/ HF pattern fetched
        self.assertFalse(any(p.startswith("samplesIDs/") for p in fetch_log),
                         f"bundled files should make per-db HF fetch unnecessary: {fetch_log}")

    def test_download_never_clobbers_existing_samplesid(self):
        """(KEY INVARIANT) a samplesID already on disk (e.g. `--what all` fetched a
        fresher one) is NEVER overwritten by the bundled install."""
        idx = _mk_variant_layout(self.src, self.name, self.ref)
        _publish(idx, self.sink)
        self._serve_from_sink([])
        sdir = os.path.join(self.workdir, "samplesIDs")
        os.makedirs(sdir, exist_ok=True)
        pre = os.path.join(sdir, f"{self.ref}_1000G.samplesID.txt")
        _mkfile(pre, "USER_FRESH\tX\tX\tX\n")
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, ref=self.ref, token="x",
                              genotypes=False)
        self.assertEqual(open(pre).read(), "USER_FRESH\tX\tX\tX\n",
                         "bundled install clobbered a pre-existing samplesID")

    def test_legacy_index_download_still_works(self):
        """(c) a legacy variant index (no bundled samplesIDs) downloads without error;
        no samplesIDs/ dir is created from the empty bundle."""
        idx = _mk_variant_layout(self.src, self.name, self.ref, with_samplesids=False)
        _publish(idx, self.sink)
        self._serve_from_sink([])
        # must not raise; the index itself still installs
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, ref=self.ref, token="x",
                              genotypes=False)
        self.assertTrue(os.path.isdir(os.path.join(self.workdir, "genome_library", self.name)))


if __name__ == "__main__":
    unittest.main()
