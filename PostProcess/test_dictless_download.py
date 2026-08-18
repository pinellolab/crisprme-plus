#!/usr/bin/env python3
"""Regression tests for the two DOWNLOAD-side out-of-the-box fixes for the
already-published dict-less index (no re-publish; the ~22GB HF artifact stays
valid). Both fixes are additive and strict NO-OPs for a normal (canonical-named,
combined-samplesID-present) dict index, so existing dict-index downloads stay
byte-identical.

GAP 3 — search-resolvable install name. The dict-less index was published as
``NRG_3_hg38-dictless+hg38_1000G_HGDP`` so it extracts as
``genome_library/NRG_3_hg38-dictless+hg38_1000G_HGDP/``, but the search resolves
by the convention ``<pam>_<N>_<ref>+<vcf>`` (ref == genome-folder basename). The
``-dictless`` in the REF segment breaks resolution. Fix: install under the
canonical dir name (strip the marker from the ref segment), preserving the
``+<vcf>`` segment (shared with the genotypes companion).

GAP 2 — combined samplesID. The dict-less install ships per-db samplesIDs
(``hg38_1000G.samplesID.txt`` + ``hg38_HGDP.samplesID.txt``) but NOT the combined
``hg38_1000G_HGDP.samplesID.txt`` the search's ``--samplesID`` listing expects.
Fix: after install, union the per-db DATA rows (header-less, dedup by SAMPLE_ID,
stable dataset order) into the combined file when it is missing.

Fully offline: crisprme_hf imports huggingface_hub lazily and these tests never
touch HF (the snapshot fetch is monkeypatched to a local copy), mirroring
test_index_bundle / test_dictless_publish. The helpers are called directly with
no HF stubbing.
"""
import os
import shutil
import tarfile
import tempfile
import unittest

import crisprme_hf as hf


# --- helpers ----------------------------------------------------------------

def _mkfile(path, content="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(content)


def _build_index_tarball(src_root, index_name, sink_dir, dicts=True):
    """Build a minimal variant-index tarball named after the ORIGINAL (possibly
    '-dictless'-marked) index_name and drop it in sink_dir. Returns the tarball
    path. The archive's internal top-level dir is the ORIGINAL index_name (that
    is how the artifact was really published)."""
    vcf = index_name.partition("+")[2]
    idx = os.path.join(src_root, "genome_library", index_name)
    ind = idx + "_INDELS"
    _mkfile(os.path.join(idx, "TSTgenome.NGG.bin"))
    _mkfile(os.path.join(ind, "NRG_3_fakechr22"))
    droot = os.path.join(src_root, "Dictionaries")
    _mkfile(os.path.join(droot, f"log_indels_{vcf}", "logchr22.txt"), "chr22\n")
    _mkfile(os.path.join(droot, f"registry_{vcf}", "reg_chr22.bin"))
    dict_dirs = [os.path.join(droot, f"log_indels_{vcf}"),
                 os.path.join(droot, f"registry_{vcf}")]
    if dicts:
        _mkfile(os.path.join(droot, f"dictionaries_{vcf}", "d_chr22.json"), "{}")
        dict_dirs.insert(0, os.path.join(droot, f"dictionaries_{vcf}"))
    tarball = os.path.join(sink_dir, f"{index_name}.tar.gz")
    hf._make_index_tarball(tarball, idx, index_name, ind,
                           {"name": index_name}, dict_dirs)
    return tarball


def _serve_snapshot(sink_dir, samples_dir=None, fetch_log=None):
    """Return a _hf_snapshot stub that serves the tarball(s) from sink_dir AND —
    when ``samples_dir`` is given — the per-db samplesID files it holds.

    Each allow_pattern is a repo-relative path like ``indexes/<name>.tar.gz`` or
    ``samplesIDs/<ref>_<db>.samplesID.txt``. The stub reproduces the repo tree under
    ``local_dir`` (as huggingface_hub.snapshot_download does): it routes an
    ``indexes/*`` pattern to ``sink_dir`` and a ``samplesIDs/*`` pattern to
    ``samples_dir``, copying the file to ``local_dir/<prefix>/<basename>`` only when
    the source exists (an unpublished file simply isn't copied — mirroring HF's
    "pattern matched nothing"). Appended to ``fetch_log`` (if given) so a test can
    assert exactly which patterns were requested (and that NONE were on the no-op
    path).
    """
    def fake_snapshot(repo, allow_patterns, local_dir, token=None):
        if fetch_log is not None:
            fetch_log.extend(allow_patterns)
        for pat in allow_patterns:
            prefix = pat.split("/", 1)[0]
            base = os.path.basename(pat)
            if prefix == "samplesIDs" and samples_dir is not None:
                srcf = os.path.join(samples_dir, base)
            else:
                srcf = os.path.join(sink_dir, base)
            if os.path.isfile(srcf):
                dest = os.path.join(local_dir, prefix)
                os.makedirs(dest, exist_ok=True)
                shutil.copy(srcf, os.path.join(dest, base))
        return local_dir
    return fake_snapshot


def _make_remote_samplesids(dbs=("1000G", "HGDP")):
    """Materialize a fake HF ``samplesIDs/`` component in a temp dir: one per-db
    file per requested dataset. Returns the dir path."""
    remote = tempfile.mkdtemp(prefix="dl_samples_")
    bodies = {"1000G": _S_1000G, "HGDP": _S_HGDP}
    for db in dbs:
        _mkfile(os.path.join(remote, f"hg38_{db}.samplesID.txt"),
                bodies.get(db, "#SAMPLE_ID\nX\n"))
    return remote


# per-db samplesID fixtures (header + a couple of data rows each)
_S_1000G = (
    "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
    "HG00096\tGBR\tEUR\tmale\n"
    "HG00097\tGBR\tEUR\tfemale\n"
)
_S_HGDP = (
    "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
    "HGDP00001\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n"
    "HGDP00003\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n"
)


def _write_per_db_samplesids(workdir):
    sdir = os.path.join(workdir, "samplesIDs")
    _mkfile(os.path.join(sdir, "hg38_1000G.samplesID.txt"), _S_1000G)
    _mkfile(os.path.join(sdir, "hg38_HGDP.samplesID.txt"), _S_HGDP)
    return sdir


# --- GAP 3: canonical_index_name (pure) -------------------------------------

class TestCanonicalIndexName(unittest.TestCase):
    def test_strips_dictless_from_ref_segment(self):
        self.assertEqual(
            hf.canonical_index_name("NRG_3_hg38-dictless+hg38_1000G_HGDP"),
            "NRG_3_hg38+hg38_1000G_HGDP",
        )

    def test_reference_only_name_unchanged(self):
        self.assertEqual(hf.canonical_index_name("NRG_3_hg38"), "NRG_3_hg38")

    def test_canonical_variant_name_unchanged(self):
        self.assertEqual(
            hf.canonical_index_name("NRG_3_hg38+hg38_1000G_HGDP"),
            "NRG_3_hg38+hg38_1000G_HGDP",
        )

    def test_marker_in_vcf_segment_preserved(self):
        """Only the base (ref) segment before the FIRST '+' is mutated; a marker
        in the vcf segment (shared w/ the genotypes companion) is preserved."""
        self.assertEqual(
            hf.canonical_index_name("NRG_3_hg38+hg38-dictless"),
            "NRG_3_hg38+hg38-dictless",
        )

    def test_idempotent(self):
        once = hf.canonical_index_name("NRG_3_hg38-dictless+hg38_1000G_HGDP")
        self.assertEqual(hf.canonical_index_name(once), once)


# --- GAP 2: synthesize_combined_samplesid (pure) ----------------------------

class TestSynthesizeCombinedSamplesID(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="sid_")

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def test_generates_from_per_db_with_order_and_headerless(self):
        _write_per_db_samplesids(self.workdir)
        out = hf.synthesize_combined_samplesid(self.workdir, "hg38_1000G_HGDP")
        self.assertEqual(out, os.path.join(self.workdir, "samplesIDs",
                                           "hg38_1000G_HGDP.samplesID.txt"))
        with open(out) as fh:
            body = fh.read()
        # header-less, all four data rows, 1000G before HGDP, trailing newline
        self.assertEqual(
            body,
            "HG00096\tGBR\tEUR\tmale\n"
            "HG00097\tGBR\tEUR\tfemale\n"
            "HGDP00001\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n"
            "HGDP00003\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n",
        )
        self.assertNotIn("#", body)

    def test_dedup_by_sample_id_first_seen_wins(self):
        sdir = os.path.join(self.workdir, "samplesIDs")
        # HG00096 appears in BOTH files with different rows -> first (1000G) wins
        _mkfile(os.path.join(sdir, "hg38_1000G.samplesID.txt"),
                "HG00096\tGBR\tEUR\tmale\n")
        _mkfile(os.path.join(sdir, "hg38_HGDP.samplesID.txt"),
                "HG00096\tOTHER\tOTHER\tfemale\nHGDP00001\tBrahui\tCSA\tmale\n")
        out = hf.synthesize_combined_samplesid(self.workdir, "hg38_1000G_HGDP")
        with open(out) as fh:
            body = fh.read()
        self.assertEqual(
            body,
            "HG00096\tGBR\tEUR\tmale\nHGDP00001\tBrahui\tCSA\tmale\n",
        )

    def test_noop_when_combined_already_present(self):
        sdir = _write_per_db_samplesids(self.workdir)
        target = os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt")
        _mkfile(target, "PRESHIPPED\tX\tY\tZ\n")
        self.assertIsNone(
            hf.synthesize_combined_samplesid(self.workdir, "hg38_1000G_HGDP"))
        with open(target) as fh:  # untouched
            self.assertEqual(fh.read(), "PRESHIPPED\tX\tY\tZ\n")

    def test_noop_single_dataset(self):
        sdir = os.path.join(self.workdir, "samplesIDs")
        _mkfile(os.path.join(sdir, "hg38_1000G.samplesID.txt"), _S_1000G)
        # vcf_name "hg38_1000G" -> dataset "1000G" (no '_') -> NO-OP
        self.assertIsNone(
            hf.synthesize_combined_samplesid(self.workdir, "hg38_1000G"))

    def test_noop_missing_component(self):
        sdir = os.path.join(self.workdir, "samplesIDs")
        _mkfile(os.path.join(sdir, "hg38_1000G.samplesID.txt"), _S_1000G)
        # HGDP component absent -> no half-union written
        self.assertIsNone(
            hf.synthesize_combined_samplesid(self.workdir, "hg38_1000G_HGDP"))
        self.assertFalse(os.path.isfile(
            os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt")))


# --- end-to-end download (offline, stubbed snapshot) ------------------------

class TestDictlessDownloadEndToEnd(unittest.TestCase):
    def setUp(self):
        self.src = tempfile.mkdtemp(prefix="dl_src_")
        self.sink = tempfile.mkdtemp(prefix="dl_sink_")
        self.workdir = tempfile.mkdtemp(prefix="dl_dst_")
        self._orig_snapshot = hf._hf_snapshot

    def tearDown(self):
        hf._hf_snapshot = self._orig_snapshot
        for d in (self.src, self.sink, self.workdir):
            shutil.rmtree(d, ignore_errors=True)

    def test_dictless_index_installs_under_canonical_name_and_syncs_samplesid(self):
        """(a) a '-dictless'-named index installs under the search-resolvable
        canonical dir name; (b) the combined samplesID is generated from the
        per-db files."""
        published = "NRG_3_hg38-dictless+hg38_1000G_HGDP"
        canonical = "NRG_3_hg38+hg38_1000G_HGDP"
        _build_index_tarball(self.src, published, self.sink, dicts=False)
        _write_per_db_samplesids(self.workdir)
        hf._hf_snapshot = _serve_snapshot(self.sink)

        dest = hf.download_component(
            "index", self.workdir, repo="local/test",
            index_name=published, token="x", genotypes=False)

        gl = os.path.join(self.workdir, "genome_library")
        # (a) installed under the CANONICAL name, NOT the published '-dictless' one
        self.assertTrue(os.path.isdir(os.path.join(gl, canonical)),
                        "index not installed under the canonical name")
        self.assertTrue(os.path.isdir(os.path.join(gl, canonical + "_INDELS")),
                        "_INDELS not installed under the canonical name")
        self.assertFalse(os.path.isdir(os.path.join(gl, published)),
                         "the '-dictless' dir must NOT be present after install")
        self.assertFalse(os.path.isdir(os.path.join(gl, published + "_INDELS")))
        # return value points at the canonical install dir
        self.assertEqual(dest, os.path.join(gl, canonical))

        # (b) combined samplesID generated from the per-db lists, named for the
        # canonical vcf segment (hg38_1000G_HGDP), header-less, ordered.
        combined = os.path.join(self.workdir, "samplesIDs",
                                "hg38_1000G_HGDP.samplesID.txt")
        self.assertTrue(os.path.isfile(combined),
                        "combined samplesID not synthesized on install")
        with open(combined) as fh:
            body = fh.read()
        self.assertEqual(
            body,
            "HG00096\tGBR\tEUR\tmale\n"
            "HG00097\tGBR\tEUR\tfemale\n"
            "HGDP00001\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n"
            "HGDP00003\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n",
        )

    def test_canonical_dict_index_is_noop_no_rename_no_regen(self):
        """(c) a CANONICAL index (no marker) with a combined samplesID already
        present downloads with NO rename and NO regeneration — the classic
        dict-index path is byte-identical."""
        canonical = "NRG_3_hg38+hg38_1000G_HGDP"
        _build_index_tarball(self.src, canonical, self.sink, dicts=True)
        # ship the combined samplesID + the per-db lists (classic dict artifact)
        sdir = _write_per_db_samplesids(self.workdir)
        preshipped = os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt")
        _mkfile(preshipped, "PRESHIPPED\tX\tY\tZ\n")
        hf._hf_snapshot = _serve_snapshot(self.sink)

        dest = hf.download_component(
            "index", self.workdir, repo="local/test",
            index_name=canonical, token="x", genotypes=False)

        gl = os.path.join(self.workdir, "genome_library")
        # installs under the SAME name it was published as (no rename)
        self.assertTrue(os.path.isdir(os.path.join(gl, canonical)))
        self.assertTrue(os.path.isdir(os.path.join(gl, canonical + "_INDELS")))
        self.assertEqual(dest, os.path.join(gl, canonical))
        # classic dict index still ships its per-sample SNP dicts
        self.assertTrue(os.path.isfile(os.path.join(
            self.workdir, "Dictionaries", "dictionaries_hg38_1000G_HGDP",
            "d_chr22.json")))
        # the pre-shipped combined samplesID is UNCHANGED (no regeneration)
        with open(preshipped) as fh:
            self.assertEqual(fh.read(), "PRESHIPPED\tX\tY\tZ\n")
        # no stray '-dictless' dir anywhere in genome_library
        self.assertFalse(any("-dictless" in n for n in os.listdir(gl)))


# --- GAP 2b: _derive_perdb_datasets (pure, shared derivation) ----------------

class TestDerivePerdbDatasets(unittest.TestCase):
    """The ensure-step and the synthesizer MUST agree byte-for-byte on which
    per-db files matter; both now derive from this single helper."""

    def test_merged_two_datasets(self):
        self.assertEqual(
            hf._derive_perdb_datasets("hg38_1000G_HGDP", "hg38"),
            ["1000G", "HGDP"],
        )

    def test_merged_three_datasets(self):
        self.assertEqual(
            hf._derive_perdb_datasets("hg38_1000G_HGDP_gnomAD", "hg38"),
            ["1000G", "HGDP", "gnomAD"],
        )

    def test_single_dataset_is_empty(self):
        # "hg38_1000G" -> strip ref -> "1000G" (no '_') -> [] (nothing to ensure)
        self.assertEqual(hf._derive_perdb_datasets("hg38_1000G", "hg38"), [])

    def test_no_ref_prefix(self):
        # a vcf_name that does not start with "<ref>_" is used verbatim
        self.assertEqual(
            hf._derive_perdb_datasets("1000G_HGDP", "hg38"), ["1000G", "HGDP"])


# --- GAP 2b: _ensure_perdb_samplesids (offline, stubbed snapshot) ------------

class TestEnsurePerdbSamplesids(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="ens_")
        self.samples_remote = _make_remote_samplesids()
        self._orig_snapshot = hf._hf_snapshot
        self.fetch_log = []
        hf._hf_snapshot = _serve_snapshot(
            self.workdir, samples_dir=self.samples_remote, fetch_log=self.fetch_log)

    def tearDown(self):
        hf._hf_snapshot = self._orig_snapshot
        for d in (self.workdir, self.samples_remote):
            shutil.rmtree(d, ignore_errors=True)

    def test_fetches_both_missing_per_db(self):
        hf._ensure_perdb_samplesids(
            "local/test", self.workdir, "hg38_1000G_HGDP", "hg38", token="x")
        sdir = os.path.join(self.workdir, "samplesIDs")
        self.assertTrue(os.path.isfile(os.path.join(sdir, "hg38_1000G.samplesID.txt")))
        self.assertTrue(os.path.isfile(os.path.join(sdir, "hg38_HGDP.samplesID.txt")))
        # exactly the two exact-path patterns were requested (no whole-folder glob)
        self.assertEqual(
            sorted(self.fetch_log),
            ["samplesIDs/hg38_1000G.samplesID.txt",
             "samplesIDs/hg38_HGDP.samplesID.txt"],
        )
        # no staging dir left behind
        self.assertFalse(os.path.isdir(os.path.join(self.workdir, ".hf_stage_samplesids")))

    def test_noop_when_all_present(self):
        _write_per_db_samplesids(self.workdir)  # both already on disk
        hf._ensure_perdb_samplesids(
            "local/test", self.workdir, "hg38_1000G_HGDP", "hg38", token="x")
        # NO fetch happened at all (no _hf_snapshot call)
        self.assertEqual(self.fetch_log, [])

    def test_fetches_only_the_missing_one(self):
        sdir = os.path.join(self.workdir, "samplesIDs")
        _mkfile(os.path.join(sdir, "hg38_1000G.samplesID.txt"), _S_1000G)  # present
        hf._ensure_perdb_samplesids(
            "local/test", self.workdir, "hg38_1000G_HGDP", "hg38", token="x")
        # only the absent HGDP file was fetched; the present 1000G was not clobbered
        self.assertEqual(self.fetch_log, ["samplesIDs/hg38_HGDP.samplesID.txt"])
        self.assertTrue(os.path.isfile(os.path.join(sdir, "hg38_HGDP.samplesID.txt")))

    def test_single_dataset_is_noop(self):
        hf._ensure_perdb_samplesids(
            "local/test", self.workdir, "hg38_1000G", "hg38", token="x")
        self.assertEqual(self.fetch_log, [])  # single dataset -> no fetch

    def test_combined_already_present_short_circuits(self):
        sdir = os.path.join(self.workdir, "samplesIDs")
        _mkfile(os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt"), "PRE\n")
        hf._ensure_perdb_samplesids(
            "local/test", self.workdir, "hg38_1000G_HGDP", "hg38", token="x")
        self.assertEqual(self.fetch_log, [])  # combined present -> nothing to ensure


# --- GAP 2b end-to-end: standalone `--what index` (offline, stubbed) ---------

class TestStandaloneIndexSamplesidFetch(unittest.TestCase):
    """A standalone `download --what index` (no prior samples fetch) must still
    end up with the combined samplesID the search needs."""

    def setUp(self):
        self.src = tempfile.mkdtemp(prefix="dl2_src_")
        self.sink = tempfile.mkdtemp(prefix="dl2_sink_")
        self.workdir = tempfile.mkdtemp(prefix="dl2_dst_")
        self.samples_remote = _make_remote_samplesids()
        self._orig_snapshot = hf._hf_snapshot
        self.fetch_log = []

    def tearDown(self):
        hf._hf_snapshot = self._orig_snapshot
        for d in (self.src, self.sink, self.workdir, self.samples_remote):
            shutil.rmtree(d, ignore_errors=True)

    def test_a_standalone_index_fetches_per_db_then_synthesizes(self):
        """(a) STANDALONE `--what index` variant install with NO per-db samplesIDs
        present -> the per-db fetch is triggered (stubbed) and the combined file is
        synthesized."""
        published = "NRG_3_hg38-dictless+hg38_1000G_HGDP"
        _build_index_tarball(self.src, published, self.sink, dicts=False)
        hf._hf_snapshot = _serve_snapshot(
            self.sink, samples_dir=self.samples_remote, fetch_log=self.fetch_log)
        # NOTE: no _write_per_db_samplesids(self.workdir) — the standalone case

        hf.download_component(
            "index", self.workdir, repo="local/test",
            index_name=published, token="x", genotypes=False)

        sdir = os.path.join(self.workdir, "samplesIDs")
        # the two per-db files were fetched into place...
        self.assertTrue(os.path.isfile(os.path.join(sdir, "hg38_1000G.samplesID.txt")))
        self.assertTrue(os.path.isfile(os.path.join(sdir, "hg38_HGDP.samplesID.txt")))
        # ...and the exact-path per-db patterns were among the snapshot requests
        self.assertIn("samplesIDs/hg38_1000G.samplesID.txt", self.fetch_log)
        self.assertIn("samplesIDs/hg38_HGDP.samplesID.txt", self.fetch_log)
        # combined samplesID synthesized from them, named for the canonical vcf seg
        combined = os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt")
        self.assertTrue(os.path.isfile(combined),
                        "combined samplesID not synthesized on standalone install")
        with open(combined) as fh:
            body = fh.read()
        self.assertEqual(
            body,
            "HG00096\tGBR\tEUR\tmale\n"
            "HG00097\tGBR\tEUR\tfemale\n"
            "HGDP00001\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n"
            "HGDP00003\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n",
        )

    def test_b_what_all_path_no_extra_fetch_still_synthesizes(self):
        """(b) NO-REGRESSION for `--what all`: the per-db files are already present
        (the 'samples' component ran first) -> NO extra samplesID fetch, and the
        combined is synthesized exactly as before."""
        published = "NRG_3_hg38-dictless+hg38_1000G_HGDP"
        _build_index_tarball(self.src, published, self.sink, dicts=False)
        _write_per_db_samplesids(self.workdir)  # the --what all preamble
        hf._hf_snapshot = _serve_snapshot(
            self.sink, samples_dir=self.samples_remote, fetch_log=self.fetch_log)

        hf.download_component(
            "index", self.workdir, repo="local/test",
            index_name=published, token="x", genotypes=False)

        sdir = os.path.join(self.workdir, "samplesIDs")
        # NO per-db samplesID fetch happened (only the index tarball was requested)
        self.assertNotIn("samplesIDs/hg38_1000G.samplesID.txt", self.fetch_log)
        self.assertNotIn("samplesIDs/hg38_HGDP.samplesID.txt", self.fetch_log)
        self.assertFalse(any(p.startswith("samplesIDs/") for p in self.fetch_log),
                         "no samplesID fetch may occur on the --what all path")
        # combined synthesized as before
        combined = os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt")
        self.assertTrue(os.path.isfile(combined))
        with open(combined) as fh:
            self.assertEqual(
                fh.read(),
                "HG00096\tGBR\tEUR\tmale\n"
                "HG00097\tGBR\tEUR\tfemale\n"
                "HGDP00001\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n"
                "HGDP00003\tBrahui\tCENTRAL_SOUTH_ASIA\tmale\n",
            )
        # no staging dir left behind
        self.assertFalse(os.path.isdir(
            os.path.join(self.workdir, ".hf_stage_samplesids")))

    def test_c_per_db_absent_from_repo_warns_and_no_crash(self):
        """(c) a per-db file ABSENT from the repo -> warned to STDOUT, no crash, and
        NO half-union combined file written (synthesizer no-ops gracefully)."""
        published = "NRG_3_hg38-dictless+hg38_1000G_HGDP"
        _build_index_tarball(self.src, published, self.sink, dicts=False)
        # remote has ONLY 1000G published; HGDP is absent
        partial_remote = _make_remote_samplesids(dbs=("1000G",))
        self.addCleanup(shutil.rmtree, partial_remote, True)
        hf._hf_snapshot = _serve_snapshot(
            self.sink, samples_dir=partial_remote, fetch_log=self.fetch_log)

        # must not raise
        dest = hf.download_component(
            "index", self.workdir, repo="local/test",
            index_name=published, token="x", genotypes=False)

        gl = os.path.join(self.workdir, "genome_library")
        canonical = "NRG_3_hg38+hg38_1000G_HGDP"
        self.assertEqual(dest, os.path.join(gl, canonical))  # index still installed
        sdir = os.path.join(self.workdir, "samplesIDs")
        # the available per-db (1000G) was fetched; the missing one (HGDP) was not
        self.assertTrue(os.path.isfile(os.path.join(sdir, "hg38_1000G.samplesID.txt")))
        self.assertFalse(os.path.isfile(os.path.join(sdir, "hg38_HGDP.samplesID.txt")))
        # NO combined file (a missing component => graceful no-op, never a half-union)
        self.assertFalse(os.path.isfile(
            os.path.join(sdir, "hg38_1000G_HGDP.samplesID.txt")),
            "must not write a half-union combined samplesID")
        # no staging dir left behind on the warn path
        self.assertFalse(os.path.isdir(
            os.path.join(self.workdir, ".hf_stage_samplesids")))


if __name__ == "__main__":
    unittest.main()
