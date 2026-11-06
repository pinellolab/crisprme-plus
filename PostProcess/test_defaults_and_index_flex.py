#!/usr/bin/env python3
"""The default PAM and index, and which index serves a search (CRISPRme 2.9.0).

SpCas9 NRR (NAA, NAG, NGA, NGG) replaced NRG as the default PAM and index: the PAM file
ships in PAMs/, the web form and its Load Example select it when it is installed, a
download of the PAM files installs it where the data repo has none, ``download --what
all`` fetches its reference index, and the advice of the CLI, complete-test and the
report names the NRR index. NRG and NGG stay supported (TheDefaultIsNRR).

Run with:

    python -m unittest test_defaults_and_index_flex -v   (from PostProcess/)

It needs no web stack: the form's functions are read out of pages/main_page.py (AST) and
run with their data lookups stubbed, as test_samplesid_emit does.
"""

import ast
import io
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent
sys.path.insert(0, str(_HERE))

import crisprme_hf as hf  # noqa: E402
import test_index_resolution as tir  # noqa: E402  (the shipped script's harness)


def _web_functions(*names, page="main_page.py", **namespace):
    """The functions ``names`` of pages/<page>, run in ``namespace`` (their data lookups
    stubbed) with the module-level assignments they read."""
    src = (_REPO / "pages" / page).read_text()
    typing = __import__("typing")
    ns = {**{n: getattr(typing, n) for n in typing.__all__}, **namespace}
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
        elif isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "DEFAULT_PAM" for t in node.targets):
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102
    return ns


def _pams(*names):
    """get_available_PAM's options for the PAM files ``names``."""
    return lambda: [{"label": n, "value": n} for n in names]


class TheDefaultIsNRR(unittest.TestCase):

    def test_the_nrr_pam_file_ships(self):
        line = (_REPO / "PAMs" / "20bp-NRR-SpCas9.txt").read_text()
        # (as the NRG file beside it: the motif line, with its newline)
        self.assertEqual(line, "N" * 20 + "NRR 3\n")
        self.assertEqual((_REPO / "PAMs" / "20bp-NRG-SpCas9.txt").read_text(),
                         "N" * 20 + "NRG 3\n")

    def test_the_form_selects_nrr_then_the_earlier_defaults(self):
        for installed, want in (
                (("20bp-NGG-SpCas9", "20bp-NRG-SpCas9", "20bp-NRR-SpCas9"), "20bp-NRR-SpCas9"),
                (("20bp-NRR-SpCas9", "20bp-NRG-SpCas9"), "20bp-NRR-SpCas9"),
                # an install from before 2.9.0: what it selected then
                (("20bp-NGG-SpCas9", "20bp-NRG-SpCas9"), "20bp-NRG-SpCas9"),
                (("20bp-NAG-SpCas9", "20bp-NGG-SpCas9"), "20bp-NGG-SpCas9"),
                (("20bp-NAG-SpCas9",), "20bp-NAG-SpCas9")):
            with self.subTest(installed=installed):
                ns = _web_functions("_default_pam", get_available_PAM=_pams(*installed))
                self.assertEqual(ns["_default_pam"]("SpCas9"), want)
        ns = _web_functions("_default_pam", get_available_PAM=_pams("23bp-TTTV-Cas12a"))
        self.assertIsNone(ns["_default_pam"]("SpCas9"))
        self.assertEqual(ns["DEFAULT_PAM"], "20bp-NRR-SpCas9")

    def _example(self, pams, genomes):
        ns = _web_functions(
            "load_example_data", "_default_pam", "_default_cas", "_default_genome",
            get_available_PAM=_pams(*pams),
            get_available_CAS=lambda: [{"label": "SpCas9", "value": "SpCas9"}],
            get_available_genomes=lambda: [{"label": g, "value": g} for g in genomes],
            get_variant_dataset_options=lambda genome: [],
            _preferred_variant=lambda options, genome: "ref")
        return ns["load_example_data"](1)

    def test_load_example_selects_the_installed_default_pam_and_genome(self):
        # it computed the installed PAM, then returned "20bp-NRG-SpCas9" and "hg38" anyway
        out = self._example(["20bp-NGG-SpCas9", "20bp-NRR-SpCas9"], ["hg38"])
        self.assertEqual(out[1:3], ["20bp-NRR-SpCas9", "hg38"])
        out = self._example(["20bp-NGG-SpCas9"], ["susScr11"])
        self.assertEqual(out[1:3], ["20bp-NGG-SpCas9", "susScr11"])
        # nothing installed: the default names, which Submit then says are missing
        out = self._example([], [])
        self.assertEqual(out[1:3], ["20bp-NRR-SpCas9", "hg38"])

    def test_download_all_fetches_the_nrr_reference_index(self):
        tree = ast.parse((_REPO / "crisprme.py").read_text())
        value = next(ast.literal_eval(n.value) for n in tree.body
                     if isinstance(n, ast.Assign)
                     and any(getattr(t, "id", "") == "DEFAULT_REFERENCE_INDEX_HG38"
                             for t in n.targets))
        self.assertEqual(value, "NRR_3_hg38")

    def test_the_advice_names_the_nrr_index(self):
        # what a user is told to download: complete-test with no index, a search with no
        # variant index, the setup and download help, the per-sample banner and the report
        for path in ("crisprme.py", "PostProcess/complete_test.py",
                     "PostProcess/generate_report.py", "PostProcess/crisprme_hf.py"):
            text = (_REPO / path).read_text()
            with self.subTest(path=path):
                self.assertIn("--index-name NRR_3_hg38+hg38_HGDP1kGP", " ".join(text.split()))
                for line in text.splitlines():
                    if "--index-name NRG_3_hg38+hg38_HGDP1kGP" in line:
                        self.fail(f"{path} still recommends the NRG index: {line.strip()}")

    def test_the_one_click_installers_download_the_nrr_index(self):
        # the route the README and the website call the easiest
        for path in ("docs/install.sh", "docs/install.ps1"):
            text = (_REPO / path).read_text(encoding="utf-8")
            with self.subTest(path=path):
                self.assertIn("--index-name NRR_3_hg38+hg38_HGDP1kGP --path /DATA", text)
                self.assertNotIn("NRG_3_hg38+hg38_HGDP1kGP", text)

    def test_install_ps1_starts_nothing_after_a_failed_download(self):
        # Its first Start ran the two downloads and the app with ";" between them, so a failed
        # download went on to the next and started the app on half the data (backlog 4.4).
        # (Windows PowerShell 5.1 has no "&&": each step tests the exit status before it.)
        line = next(ln for ln in (_REPO / "docs" / "install.ps1").read_text(
            encoding="utf-8").splitlines() if "download --what all" in ln)
        index_step = line.index("download --what index")
        self.assertIn("; if (`$LASTEXITCODE -eq 0) { docker run", line[:index_step])
        self.assertIn("; if (`$LASTEXITCODE -eq 0) { $run;", line[index_step:])


class APamDownloadInstallsTheDefaultPam(unittest.TestCase):
    """The data repo's pams/ held NGG, NNN and NRG only: the default PAM file comes with
    CRISPRme, and a download of the PAM files installs it beside them."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def download(self, hosted):
        def snapshot(repo, patterns, local_dir, token=None, **kw):
            os.makedirs(os.path.join(local_dir, "pams"), exist_ok=True)
            for name, line in hosted.items():
                Path(local_dir, "pams", name).write_text(line)
            return local_dir

        out = io.StringIO()
        with mock.patch.object(hf, "_hf_snapshot", snapshot), redirect_stdout(out):
            dest = hf.download_component("pams", str(self.tmp), repo="me/data")
        return Path(dest), out.getvalue()

    def test_it_is_installed_when_the_hub_has_none(self):
        dest, said = self.download({"20bp-NRG-SpCas9.txt": "N" * 20 + "NRG 3\n",
                                    "20bp-NGG-SpCas9.txt": "N" * 20 + "NGG 3\n"})
        self.assertEqual(sorted(os.listdir(dest)), ["20bp-NGG-SpCas9.txt",
                                                    "20bp-NRG-SpCas9.txt",
                                                    "20bp-NRR-SpCas9.txt"])
        self.assertEqual((dest / "20bp-NRR-SpCas9.txt").read_text(), "N" * 20 + "NRR 3\n")
        self.assertIn("Installed the default PAM file 20bp-NRR-SpCas9.txt", said)
        # never over a file of that name: the user's own, or the hub's once it hosts one
        (dest / "20bp-NRR-SpCas9.txt").write_text("mine\n")
        dest, said = self.download({})
        self.assertEqual((dest / "20bp-NRR-SpCas9.txt").read_text(), "mine\n")
        self.assertNotIn("Installed the default PAM", said)
        dest, said = self.download({"20bp-NRR-SpCas9.txt": "N" * 20 + "NRR 3\n"})
        self.assertEqual((dest / "20bp-NRR-SpCas9.txt").read_text(), "N" * 20 + "NRR 3\n")
        self.assertNotIn("Installed the default PAM", said)


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class ASearchWithNoBulgesBuildsNoIndex(tir._ShippedScript):
    """A search with no bulges scans the genome's FASTA files, and a variant one the
    enriched genome's when it is here whole: it opens no index. It built the reference
    index <PAM>_1_<genome> all the same (hours on hg38), and a variant one, which nothing
    read; under --index-path it stopped on a library that held no reference index."""

    def log(self, run_number):
        return (self.root / f"run{run_number}" / "out" / "log.txt").read_text()

    def test_a_reference_search_builds_none(self):
        run = self.run_pipeline(bulges=0, variant=False)
        self.assertEqual(run.builds, [])
        self.assertIn("No reference index needed: a search with no bulges scans the genome",
                      run.stdout)
        search = next(c for c in run.calls if c[0] == "search")
        self.assertEqual((Path(search[1]).name, "-r" in search), ("hg38", True))
        # the web's progress reads the step as done from this line
        self.assertIn("Index-genome Reference\tEnd", self.log(self.runs))

    def test_nor_stops_on_an_index_path_library_that_holds_none(self):
        shared = self.root / "shared_library"
        shared.mkdir()
        run = self.run_pipeline(bulges=0, variant=False, index_path=shared)
        self.assertNotIn("no matching reference index", run.stderr)
        self.assertEqual(run.builds, [])
        self.assertEqual(Path(next(c for c in run.calls if c[0] == "search")[1]).name, "hg38")

    def test_a_variant_search_of_a_whole_enriched_genome_builds_neither(self):
        enriched = self.home / "Genomes" / f"hg38+{self.VCF}"
        enriched.mkdir()
        (enriched / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        (self.home / "Genomes" / f"hg38+{self.VCF}_INDELS").mkdir()
        run = self.run_pipeline(bulges=0)
        self.assertEqual([c[1] for c in run.builds if c[1] in ("hg38", f"hg38+{self.VCF}")],
                         [])
        self.assertIn("No variant index needed", run.stdout)
        self.assertIn("Index-genome Variant\tEnd", self.log(self.runs))
        self.assertEqual([Path(c[1]).name for c in run.calls if c[0] == "search"],
                         ["hg38", f"hg38+{self.VCF}"])
        # a search with bulges builds both, as before
        run = self.run_pipeline(bulges=1)
        self.assertEqual([c[1] for c in run.builds if c[1] in ("hg38", f"hg38+{self.VCF}")],
                         ["hg38", f"hg38+{self.VCF}"])


def _lib(root, *names, side=None):
    """A genome_library under ``root`` holding the finished indexes ``names`` (a .bin each);
    with ``side`` (3 or 5), each records in a .pam_build that it was built for a PAM on that
    side of a 20-nt guide, as build-index-only does."""
    lib = Path(root, "genome_library")
    for name in names:
        Path(lib, name).mkdir(parents=True, exist_ok=True)
        Path(lib, name, "chr1.bin").write_bytes(b"")
        if side:
            pam = name.split("_")[0]
            Path(lib, name, ".pam_build").write_text(
                f"x {20 + len(pam)} {len(pam) if side == 3 else -len(pam)}\n")
    return lib


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class ADedicatedIndexServesAPamItHolds(unittest.TestCase):
    """Under --prune-pam a dedicated index serves a PAM its own holds at every position
    (NRR_3 for NGG, NAG, NRG or NGA; NRG_3 for NGG or NAG): the engine prunes it to that PAM
    at leaf time, as it does an all-N index. It was refused with --index-path, and a search
    without it built an index of the PAM. The shipped resolver, and the Python mirrors that
    must not disagree with it: index_fit.pam_serves (crisprme.py, unfinished_marks, the
    web), crisprme.py's _index_the_search_takes, the web's index_max_bulges."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def take(self, lib, pam, prune, pos=3, bmax=2):
        return tir.resolve(lib, pam, "hg38", bmax, prune=prune, pos=pos)

    def test_nrr_serves_the_pams_it_holds_pruned_only(self):
        lib = _lib(self.tmp.name, "NRR_3_hg38")
        for pam in ("NGG", "NAG", "NRG", "NGA", "NAA"):
            with self.subTest(pam=pam):
                self.assertEqual(self.take(lib, pam, prune=1), "NRR_3_hg38")
                self.assertEqual(self.take(lib, pam, prune=0), "")
        for pam in ("NGC", "NYG", "NNG", "NRRT", "NG"):
            with self.subTest(pam=pam):
                self.assertEqual(self.take(lib, pam, prune=1), "")

    def test_nrg_serves_ngg_and_nag_and_not_nga(self):
        lib = _lib(self.tmp.name, "NRG_3_hg38")
        self.assertEqual(self.take(lib, "NGG", 1), "NRG_3_hg38")
        self.assertEqual(self.take(lib, "NAG", 1), "NRG_3_hg38")
        self.assertEqual(self.take(lib, "NGA", 1), "")

    def test_the_pams_own_then_a_dedicated_one_then_a_pam_less_one(self):
        lib = _lib(self.tmp.name, "NNN_3_hg38", "NNNNNNNN_3_hg38", "NRR_3_hg38")
        self.assertEqual(self.take(lib, "NGG", 1), "NRR_3_hg38")
        # (unpruned, as before: the all-N one of its width, which the filter narrows)
        self.assertEqual(self.take(lib, "NGG", 0), "NNN_3_hg38")
        _lib(self.tmp.name, "NGG_3_hg38")
        self.assertEqual(self.take(lib, "NGG", 1), "NGG_3_hg38")
        # a smaller N of the dedicated one wins among those, as it did
        _lib(self.tmp.name, "NRG_3_hg38", "NRR_4_hg38")
        self.assertEqual(self.take(lib, "NAG", 1), "NRG_3_hg38")

    def test_only_for_a_pam_on_its_own_side(self):
        lib = _lib(self.tmp.name, "NRR_3_hg38", side=5)
        self.assertEqual(self.take(lib, "NGG", 1), "")
        # a 5' one, recorded, serves a 5' PAM it holds
        lib = _lib(tempfile.mkdtemp(dir=self.tmp.name), "TTTV_2_hg38", side=5)
        self.assertEqual(self.take(lib, "TTTA", 1, pos=-4, bmax=1), "TTTV_2_hg38")
        self.assertEqual(self.take(lib, "TTTA", 1, pos=4, bmax=1), "")
        # with no record, only at the published widths, as a 3' build, as an all-N one
        lib = _lib(tempfile.mkdtemp(dir=self.tmp.name), "TTTV_2_hg38", "NRR_2_hg38")
        self.assertEqual(self.take(lib, "TTTA", 1, pos=-4, bmax=1), "")
        self.assertEqual(self.take(lib, "TTTA", 1, pos=4, bmax=1), "")
        self.assertEqual(self.take(lib, "NGG", 1, bmax=1), "NRR_2_hg38")

    def test_the_python_mirrors_take_the_same(self):
        import index_fit
        import test_bulge_cap as tbc  # (crisprme.py's definitions and pages_utils, loaded)
        names = ["NGG_2_hg38", "NRG_3_hg38", "NRR_3_hg38", "NNN_3_hg38", "NNNNNNNN_3_hg38",
                 "NRR_2_hg38", "NAG_3_hg38"]
        for k in range(1, 1 << len(names)):
            installed = [n for i, n in enumerate(names) if k >> i & 1]
            root = tempfile.mkdtemp(dir=self.tmp.name)
            lib = _lib(root, *installed)
            for pam in ("NGG", "NAG", "NGA", "NRR"):
                for prune in (0, 1):
                    shell = self.take(lib, pam, prune)
                    taken = tbc.crisprme._index_the_search_takes(
                        [str(lib)], pam, "hg38", 2, bool(prune), True, 20, 1, 1, False)
                    with self.subTest(installed=installed, pam=pam, prune=prune):
                        self.assertEqual(Path(taken).name if taken else "", shell)
                        served = {n for n in installed if index_fit.pam_serves(
                            n.split("_")[0], pam, bool(prune), True) is not None}
                        self.assertEqual(bool(served), bool(shell))


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class ASearchPrunesByItself(tir._ShippedScript):
    """When no installed index serves a search unpruned and one serves it pruned (a wider
    all-N one, a dedicated one of a PAM that holds the search's), the search prunes as if
    --prune-pam was given, and says so first. It stopped for want of a variant index (the
    flag told to the user), or built a reference one. The web launches the script with
    pruning off, so its searches do the same. The run's parameters file names the indexes
    the search took, for the report."""

    def params(self):
        return (self.root / f"run{self.runs}" / "out" / "Params.txt").read_text()

    def run_with_params(self, **kw):
        # (complete-search writes the parameters file before it launches the script)
        run = self.launch(**kw)
        (run.work / "out" / "Params.txt").write_text("Genome_selected\thg38\n")
        return self.finish(run)

    def test_an_nrr_variant_index_serves_an_ngg_search_pruned(self):
        self.install("NRR_3")
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0)
        self.assertIn("NOTE: no installed index serves PAM NGG unpruned, and "
                      f"NRR_3_hg38+{self.VCF} serves it pruned to it at leaf time: searching "
                      "with --prune-pam, as if it was given", run.stdout)
        self.assertEqual(run.searched(variant=True), (f"NRR_3_hg38+{self.VCF}", True))
        self.assertEqual(run.searched(variant=False), ("NRR_3_hg38", True))
        self.assertEqual(run.builds, [])
        # beside an index that serves it unpruned, it is searched as before
        for installed, taken in ((("NGG_3", "NRR_3"), "NGG_3"), (("NNN_3", "NRR_3"), "NNN_3")):
            with self.subTest(installed=installed):
                shutil.rmtree(self.library)
                self.library.mkdir()
                self.install(*installed)
                run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0)
                self.assertNotIn("as if it was given", run.stdout)
                self.assertEqual(run.searched(variant=True), (f"{taken}_hg38+{self.VCF}", False))

    def test_a_wider_pam_less_index_alone_serves_pruned(self):
        self.install("NNNNNNNN_3")
        run = self.run_pipeline(pam="20bp-NRR-SpCas9", prune=0)
        self.assertEqual(run.searched(variant=True), (f"NNNNNNNN_3_hg38+{self.VCF}", True))
        self.assertIn("searching with --prune-pam, as if it was given", run.stdout)

    def test_a_reference_search_takes_one_rather_than_build_one(self):
        self.install("NRR_3", variant=False)
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", variant=False, prune=0)
        self.assertEqual(run.builds, [])
        self.assertEqual(run.searched(variant=False), ("NRR_3_hg38", True))

    def test_the_parameters_file_names_the_indexes_searched(self):
        self.install("NRR_3")
        run = self.run_with_params(pam="20bp-NGG-SpCas9", prune=0)
        params = self.params()
        self.assertIn(f"Reference_index\tNRR_3_hg38 (in {self.library}): built for the PAM NRR "
                      "on the 3' side of a guide of a length it does not record, with up to 2 "
                      "bulge(s) of each type, pruned to NGG at leaf time (--prune-pam)\n", params)
        self.assertIn(f"Prune_pam\tauto\n", params)
        self.assertIn(f"Variant_index:{self.VCF}\tNRR_3_hg38+{self.VCF} (in {self.library}): "
                      "built for the PAM NRR on the 3' side", params)
        self.assertIn(f"; its indels: NRR_3_hg38+{self.VCF}_INDELS\n", params)
        # an index whose trees record the guide length says it
        tir.built_for_guide(self.library / "NRR_3_hg38", 20)
        self.run_with_params(pam="20bp-NGG-SpCas9", prune=0)
        self.assertIn("built for the PAM NRR on the 3' side of a 20-nt guide, with up to 2 "
                      "bulge(s) of each type", self.params())
        # with no bulges, no index: the genome is scanned
        run = self.run_with_params(pam="20bp-NRR-SpCas9", bulges=0, variant=False)
        self.assertIn("Reference_index\tnone: a search with no bulges scans the genome\n",
                      self.params())
        self.assertNotIn("Prune_pam", self.params())


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class AnIndexNameMayCarryItsGuideAndSide(unittest.TestCase):
    """A name of the 2.9.0 form carries the guide (spacer) length and the side of the PAM, in
    the order of the target on the strand: <guide>-<PAM>_<N>_... for a 3' PAM, <PAM>-<guide>_...
    for a 5' one, so two builds that took one name (a 20- and a 23-nt Cas12a guide; a 3' and
    a 5' PAM of one motif) can lie side by side. Every reader takes both forms, and of the same
    index under both names, the newer one. Nothing took the new form before."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def take(self, lib, pam, prune=0, pos=3, bmax=2, any_side=""):
        return tir.resolve(lib, pam, "hg38", bmax, prune=prune, pos=pos, any_side=any_side)

    def test_the_parser(self):
        import index_fit
        self.assertEqual(index_fit.parse_index_name("20-NRR_3_hg38+hg38_HGDP1kGP"),
                         ("NRR", 3, "hg38+hg38_HGDP1kGP", 20, 3))
        self.assertEqual(index_fit.parse_index_name("TTTV-23_2_hg38"), ("TTTV", 2, "hg38", 23, 5))
        self.assertEqual(index_fit.parse_index_name("NRR_3_hg38"), ("NRR", 3, "hg38", None, None))
        for junk in ("NRR", "NRR_x_hg38", "20-NRR_INDELS"):  # no index name
            with self.subTest(name=junk):
                self.assertIsNone(index_fit.parse_index_name(junk))
        # a token that is not of the new form is an older name's PAM, as before
        self.assertEqual(index_fit.parse_index_name("-NRR_3_hg38")[3:], (None, None))

    def test_the_resolver_takes_it_for_its_own_side_alone(self):
        lib = _lib(self.tmp.name, "20-NRR_3_hg38")
        self.assertEqual(self.take(lib, "NRR"), "20-NRR_3_hg38")
        self.assertEqual(self.take(lib, "NRR", pos=-3), "")  # a 5' search: not its side
        self.assertEqual(self.take(lib, "NRR", pos=-3, any_side="any"), "20-NRR_3_hg38")
        self.assertEqual(self.take(lib, "NGG", prune=1), "20-NRR_3_hg38")  # it holds NGG
        lib = _lib(tempfile.mkdtemp(dir=self.tmp.name), "TTTV-23_2_hg38")
        self.assertEqual(self.take(lib, "TTTV", pos=-4, bmax=1), "TTTV-23_2_hg38")
        self.assertEqual(self.take(lib, "TTTV", pos=4, bmax=1), "")

    def test_of_one_index_under_both_names_the_new_one(self):
        lib = _lib(self.tmp.name, "NRR_3_hg38", "20-NRR_3_hg38")
        for locale in ("C", "en_US.UTF-8"):
            with self.subTest(locale=locale):
                self.assertEqual(tir.resolve(lib, "NRR", "hg38", 2, locale=locale),
                                 "20-NRR_3_hg38")
        # the rest of the order comes first: a tighter budget under the older name wins
        lib = _lib(tempfile.mkdtemp(dir=self.tmp.name), "NRR_2_hg38", "20-NRR_3_hg38")
        self.assertEqual(self.take(lib, "NRR"), "NRR_2_hg38")

    def test_the_python_mirrors_and_the_web_take_it_too(self):
        import crisprme_hf
        import index_fit
        import test_bulge_cap as tbc
        lib = _lib(self.tmp.name, "NRR_3_hg38", "20-NRR_3_hg38", "NRR-20_3_hg38")
        taken = tbc.crisprme._index_the_search_takes([str(lib)], "NRR", "hg38", 2, False,
                                                     True, 20, 1, 1, False)
        self.assertEqual(Path(taken).name, self.take(lib, "NRR"))
        self.assertTrue(tbc.crisprme._index_serves_search("20-NRR_3_hg38", "hg38", "NGG", 2,
                                                          True, True))
        self.assertFalse(tbc.crisprme._index_serves_search("NRR-20_3_hg38", "hg38", "NRR", 2,
                                                           False, True))
        # the router reads the flank of a PAM-less one in either form
        self.assertEqual(crisprme_hf.index_flank_width("20-NNNNNNNN_3_hg38+hg38_HGDP1kGP"), 8)
        self.assertEqual(crisprme_hf.resolve_index_repo("20-NNNNNNNN_3_hg38"),
                         crisprme_hf.GENERAL_INDEX_HF_REPO)
        # the web form counts it, on its side
        pu = tbc.pages_utils
        root = Path(tempfile.mkdtemp(dir=self.tmp.name))
        _lib(root, "20-NRR_3_hg38", "NNNN-20_3_hg38+hg38_x")
        (root / "PAMs").mkdir()
        (root / "PAMs" / "20bp-NRR-SpCas9.txt").write_text("N" * 20 + "NRR 3\n")
        with mock.patch.object(pu, "current_working_directory", str(root) + os.sep):
            self.assertEqual(pu.index_max_bulges("hg38", "20bp-NRR-SpCas9"), 2)
            self.assertEqual(pu.index_max_bulges("hg38", "20bp-NRR-SpCas9", "x"), 0)
            self.assertEqual(pu.index_build_pam(str(root / "genome_library" /
                                                    "NNNN-20_3_hg38+hg38_x")), (None, 24, -4))
        # complete-test's smoke and the INDELS view read it
        import complete_test
        self.assertEqual(complete_test._parse_index_name("20-NRR_3_hg38+hg38_x"),
                         ("NRR", "hg38", True))
        self.assertEqual(index_fit.index_pam("TTTV-23_2_hg38"), "TTTV")


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class TheSearchScriptTakesTheNewName(tir._ShippedScript):
    """A variant index and its INDELS index under the 2.9.0 name: the scans take them, the
    searches are handed them, and the indel search is told the name."""

    def test_a_download_of_the_new_name_is_searched(self):
        self.install("NRR_3")
        for old in sorted(self.library.iterdir()):
            old.rename(self.library / ("20-" + old.name))
        run = self.run_pipeline(pam="20bp-NRR-SpCas9")
        self.assertEqual(run.reference, "20-NRR_3_hg38")
        self.assertEqual(run.variant, f"20-NRR_3_hg38+{self.VCF}")
        self.assertEqual(run.indels, f"20-NRR_3_hg38+{self.VCF}_INDELS")
        self.assertEqual(run.searched(variant=True), (f"20-NRR_3_hg38+{self.VCF}", False))
        # and an NGG search of it, pruned by itself
        run = self.run_pipeline(pam="20bp-NGG-SpCas9")
        self.assertEqual(run.searched(variant=True), (f"20-NRR_3_hg38+{self.VCF}", True))


class ThePamPickerOfAPamLessIndex(unittest.TestCase):
    """The web form's PAM picker of a PAM-less index (pages_utils, the logic its callbacks
    run; the callbacks themselves: test_web_bulge_notes in the image): which index it reads,
    which PAMs it offers and refuses, and the PAM file it writes. There was none: a PAM-less
    index served from the web only the PAM files installed, of exactly its width."""

    def setUp(self):
        import test_bulge_cap as tbc
        self.pu = tbc.pages_utils
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        patcher = mock.patch.object(self.pu, "current_working_directory", str(self.root) + os.sep)
        patcher.start()
        self.addCleanup(patcher.stop)
        (self.root / "PAMs").mkdir()

    def test_the_widest_pam_less_index_of_the_genome_or_the_dataset(self):
        self.assertIsNone(self.pu.pamless_index_geometry("hg38"))
        _lib(self.root, "NRR_3_hg38", "NNN_3_hg38", "NNNNNNNN_3_hg38",
             "NNN_3_hg38+hg38_HGDP1kGP")
        self.assertEqual(self.pu.pamless_index_geometry("hg38"),
                         {"name": "NNNNNNNN_3_hg38", "k": 8, "side": 3, "guide_len": 20})
        self.assertEqual(self.pu.pamless_index_geometry("hg38", "HGDP1kGP")["name"],
                         "NNN_3_hg38+hg38_HGDP1kGP")
        self.assertIsNone(self.pu.pamless_index_geometry("hg38", "HPRC"))

    def test_what_it_offers_and_refuses(self):
        k8 = {"name": "NNNNNNNN_3_hg38", "k": 8, "side": 3, "guide_len": 20}
        offered = {o["value"]: o for o in self.pu.pamless_presets(k8)}
        self.assertFalse(offered["SpCas9:NRR"]["disabled"])
        self.assertEqual(list(offered)[0], "SpCas9:NRR")  # the default first
        for ok in ("SpG:NGN", "SpRY:NRN", "SpCas9VRER:NGCG", "NmeCas9:NNNNGATT"):
            self.assertFalse(offered[ok]["disabled"], ok)
        self.assertIn("21-nt guide", offered["SaCas9:NNGRRT"]["label"])
        self.assertTrue(offered["Cas12a:TTTV"]["disabled"])
        n3 = {"name": "NNN_3_hg38", "k": 3, "side": 3, "guide_len": 20}
        self.assertIn("longer than this index's 3-nt flank",
                      self.pu.pamless_pam_problem("NGCG", n3))
        self.assertIsNone(self.pu.pamless_pam_problem("ngn", n3))
        self.assertIn("no IUPAC base", self.pu.pamless_pam_problem("NG;", n3))
        five = {"name": "NNNN_3_hg38", "k": 4, "side": 5, "guide_len": 23}
        self.assertIsNone(self.pu.pamless_pam_problem("TTTV", five))
        self.assertIn("as long as the flank", self.pu.pamless_pam_problem("TTV", five))

    def test_the_pam_file_it_writes(self):
        k8 = {"name": "NNNNNNNN_3_hg38", "k": 8, "side": 3, "guide_len": 20}
        name = self.pu.write_pamless_pam_file("ngn", k8, "SpG")
        self.assertEqual(name, "20bp-NGN-SpG")
        self.assertEqual((self.root / "PAMs" / "20bp-NGN-SpG.txt").read_text(), "N" * 20 + "NGN 3\n")
        five = {"name": "NNNN_3_hg38", "k": 4, "side": 5, "guide_len": 23}
        self.assertEqual(self.pu.write_pamless_pam_file("TTTV", five), "TTTV-23bp-Custom")
        self.assertEqual((self.root / "PAMs" / "TTTV-23bp-Custom.txt").read_text(),
                         "TTTV" + "N" * 23 + " -4\n")
        # never over a file of that name
        (self.root / "PAMs" / "20bp-NGN-SpG.txt").write_text("mine\n")
        self.pu.write_pamless_pam_file("NGN", k8, "SpG")
        self.assertEqual((self.root / "PAMs" / "20bp-NGN-SpG.txt").read_text(), "mine\n")


class SetupFailsWhenItsDownloadFails(unittest.TestCase):
    """`crisprme.py setup` caught the legacy download script's failure, wrote "Legacy database
    setup exited with error code 1" with no end of line, and ended with exit status 0: with
    no network, a requests traceback and a status 0 (backlog 4.4, FR2 cli notes (a))."""

    def test_it_ends_with_the_scripts_status(self):
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        (tmp / "setup_legacy_database.py").write_text("import sys\nsys.exit(3)\n")
        err = io.StringIO()
        with mock.patch.object(cm, "script_path", str(tmp)), \
                mock.patch.object(cm, "input_args", ["crisprme.py", "setup", "--path", str(tmp)]), \
                mock.patch.object(sys, "stderr", err), self.assertRaises(SystemExit) as stop:
            cm.setup_database()
        self.assertEqual(stop.exception.code, 3)
        self.assertEqual(err.getvalue(), "Legacy database setup exited with error code 3\n")

    def test_its_progress_lines_come_as_they_are_printed(self):
        # Into a pipe or a file, as a merged log has them, the script's progress lines came
        # at its end, after its errors (backlog 4.4, FR2 cli notes (c)).
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        (tmp / "setup_legacy_database.py").write_text(
            "import sys\nprint('Downloading hg38')\n"
            "sys.stderr.write('ERROR: no network\\n')\nsys.exit(1)\n")
        env = {k: v for k, v in os.environ.items() if k != "PYTHONUNBUFFERED"}
        with open(tmp / "log", "wb") as log:
            saved = os.dup(1), os.dup(2)
            os.dup2(log.fileno(), 1)
            os.dup2(log.fileno(), 2)
            try:
                with mock.patch.object(cm, "script_path", str(tmp)), \
                        mock.patch.object(cm, "input_args",
                                          ["crisprme.py", "setup", "--path", str(tmp)]), \
                        mock.patch.dict(os.environ, env, clear=True), \
                        self.assertRaises(SystemExit):
                    cm.setup_database()
            finally:
                os.dup2(saved[0], 1)
                os.dup2(saved[1], 2)
                os.close(saved[0])
                os.close(saved[1])
        self.assertEqual((tmp / "log").read_text().splitlines()[:2],
                         ["Downloading hg38", "ERROR: no network"])


class TheTestDataUrlsHaveOneSlash(unittest.TestCase):
    """setup's test-data base URL ended in "/" and each URL added another
    (.../test/data//samplesIDs/...): GitHub answered 307 (backlog 4.4, FR2 cli notes (e))."""

    def test_no_double_slash(self):
        import importlib
        import setup_legacy_database as sld
        for env in ({}, {"CRISPRME_TESTDATA_URL": "https://mirror.example/data/"}):
            with self.subTest(env=env), mock.patch.dict(os.environ, env):
                sld = importlib.reload(sld)
                self.assertFalse(sld.TEST_DATA_BASE_URL.endswith("/"))
                self.assertNotIn("//samplesIDs", f"{sld.TEST_DATA_BASE_URL}/samplesIDs/x.txt")
        importlib.reload(sld)


class GenerateReportSaysWhyThereIsNothingToReport(unittest.TestCase):
    """generate-report on a folder with no integrated results -- a search that found no
    off-target, an assembly-search's combined folder -- ended on a FileNotFoundError
    traceback (OF-12). It says why, in one line, and exits 1."""

    def test_a_zero_hit_run_and_a_combined_folder(self):
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        for name, said in (("zero", "holds no integrated results (*integrated_results.tsv), "
                                    "so there is nothing to report"),
                           ("job_combined", "is an assembly-search's combined folder")):
            with self.subTest(folder=name):
                folder = tmp / name
                folder.mkdir()
                (folder / "log.txt").write_text("No off-targets found\n")
                err = io.StringIO()
                with mock.patch.object(cm, "input_args", ["crisprme.py", "generate-report",
                                                          "--result-dir", str(folder)]), \
                        mock.patch.object(sys, "stderr", err), \
                        self.assertRaises(SystemExit) as stop:
                    cm.generate_report()
                self.assertEqual(stop.exception.code, 1)
                self.assertIn(said, err.getvalue())
                self.assertEqual(len(err.getvalue().splitlines()), 1)


class ThePamLessNucleaseIsNamedAsTheWebNamesIt(unittest.TestCase):
    """A search of 20bp-NNN-NO-PAM recorded the nuclease "NO" from the CLI and "NO-PAM" from
    the web (OF-13): the CLI now takes all of the name after the motif."""

    def test_the_name_after_the_motif(self):
        import test_bulge_cap as tbc
        for pam_file, nuclease in (("PAMs/20bp-NNN-NO-PAM.txt", "NO-PAM"),
                                   ("20bp-NRR-SpCas9.txt", "SpCas9"),
                                   ("/x/23bp-TTTV-Cas12a.txt", "Cas12a"),
                                   ("TTTV-23bp-Cas12a.txt", "Cas12a")):
            with self.subTest(pam_file=pam_file):
                self.assertEqual(tbc.crisprme._nuclease_of(pam_file), nuclease)


class TheSourceEnvironmentHasWhatTheImageHas(unittest.TestCase):
    """environment.yml left out pigz, which the Docker image installs and publish-index packs
    an index with on every core (install P9). Every conda package of the image's main
    environment is in environment.yml."""

    def test_the_images_packages_are_in_environment_yml(self):
        import re
        docker = (_REPO / "Dockerfile").read_text()
        env = (_REPO / "environment.yml").read_text()
        listed = {re.split(r"[=<>\s#]", line.strip()[2:])[0]
                  for line in env.splitlines() if line.strip().startswith("- ")}
        # the image's line of the command-line tools, pigz among them
        line = next(ln for ln in docker.splitlines() if " pigz " in f" {ln} ")
        packages = [t for t in line.split() if t != "\\"]
        self.assertIn("pigz", packages)
        for package in packages:
            with self.subTest(package=package):
                self.assertIn(package, listed)


class TheChecklistsHpcRouteIsANewUsersRun(unittest.TestCase):
    """RELEASE_CHECKLIST's HPC route ran its apptainer commands with the host's environment
    (CONDA_* leaked in on ml008), and steps 2-4b reused step 1's writable sandbox with the
    caches it left (list_282.md, from the 2.8.1 amd64 local gate)."""

    def test_it_says_to_run_clean_and_to_rebuild_the_sandbox(self):
        text = " ".join((_REPO / "docs" / "RELEASE_CHECKLIST.md").read_text().split())
        route = text[text.index("**On an HPC node without Docker**"):]
        self.assertIn("Run the block from a clean host environment, e.g. in `env -i", route)
        self.assertIn("after step 1, build the sandbox again from the SIF", route)


class AGenomeFolderWithASpaceIsRefusedUpFront(unittest.TestCase):
    """A genome folder with a space in its path failed later with a misleading path: the
    search takes it into shell commands unquoted (IDX-19). It is refused in one line."""

    def test_refused_naming_it(self):
        import test_bulge_cap as tbc
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        (tmp / "my genome").mkdir()
        (tmp / "hg38").mkdir()
        err = io.StringIO()
        with mock.patch.object(sys, "stderr", err), self.assertRaises(SystemExit) as stop:
            tbc.crisprme._check_genome(["--genome", str(tmp / "my genome")])
        self.assertEqual(stop.exception.code, 1)
        self.assertIn("has a space in its path", err.getvalue())
        self.assertEqual(len(err.getvalue().splitlines()), 1)
        self.assertEqual(tbc.crisprme._check_genome(["--genome", str(tmp / "hg38")]),
                         str(tmp / "hg38"))


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class TheVcfListIsNotEdited(tir._ShippedScript):
    """The search script appended a newline to the user's --vcf list, and cut its last byte
    at the end: a run that ended early left the newline in it, and a read-only list failed
    the write, on stderr, which fails the stage (RV4-1). Its --samplesID list too, when
    writable. It reads copies."""

    def test_the_list_is_left_as_it_was(self):
        # (the run stops at its variant search, as one that ends early does; a read-only
        # list cannot be tested as root, who writes it anyway)
        self.install("NRR_3")
        run = self.launch(pam="20bp-NRR-SpCas9")
        done = self.finish(run)
        self.assertEqual((run.work / "vcfs.txt").read_text(), self.VCF + "\n")
        self.assertEqual(done.searched(variant=True)[0], f"NRR_3_hg38+{self.VCF}")
        self.assertEqual((run.work / "out" / ".vcf_list_searched.txt").read_text(),
                         self.VCF + "\n\n")

    def test_nor_the_samples_list(self):
        # (a last line with no newline, as an editor may leave it: the run added one, and
        # cut it off again at its end, which a run that ends early never reaches)
        self.install("NRR_3")
        samples = self.root / "my_samples.txt"
        samples.write_text(f"{self.VCF}.samplesID.txt")
        self.finish(self.launch(pam="20bp-NRR-SpCas9", samples=samples))
        self.assertEqual(samples.read_text(), f"{self.VCF}.samplesID.txt")


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class AFlankWiderThanA5PrimePamIsSaid(tir._ShippedScript):
    """A 5' search refused for want of an index named none of the PAM-less indexes of its
    side there whose flank is wider than its PAM: the scan skips them without a word, as a 5'
    PAM is read off a flank of its own width only, and the user was left thinking the index
    was recorded wrong (cli5 D7). The refusal says so, and what would serve."""

    SAID = ("       NNNNNNNN_2_hg38 is a PAM-less index whose flank is 8 nt, and TTTV is a 5' "
            "PAM of 4 nt: a 5' PAM is read off a flank of its own width only, so this search "
            "needs NNNN_2_hg38 or TTTV_2_hg38, built for a 5' PAM\n")

    def test_under_index_path(self):
        self.install("NNNNNNNN_2", variant=False)
        tir.built_for(self.library / "NNNNNNNN_2_hg38", -8)
        for prune in (0, 1):
            with self.subTest(prune=prune):
                run = self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=prune, variant=False,
                                        index_path=self.library)
                self.assertIsNone(run.searched(variant=False))
                self.assertIn("ERROR: no matching reference index under --index-path",
                              run.stderr)
                self.assertIn(self.SAID, run.stderr)
        # a 3' PAM is not told it
        run = self.run_pipeline(pam="20bp-NGCG-SpCas9VRER", prune=0, variant=False,
                                index_path=self.library)
        self.assertNotIn("whose flank is", run.stderr)


class AFlankWiderThanA5PrimePamIsSaidByThePreflight(unittest.TestCase):
    """The same, in complete-search's pre-flight refusal of a variant search."""

    def test_in_the_preflight_of_a_variant_search(self):
        import struct
        import test_bulge_cap as tbc
        root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, root, ignore_errors=True)
        cwd = tbc.crisprme.current_working_directory
        tbc.crisprme.current_working_directory = f"{root}/"
        self.addCleanup(setattr, tbc.crisprme, "current_working_directory", cwd)
        for folder in ("Dictionaries/registry_ds", "Dictionaries/log_indels_ds",
                       "genome_library/NNNNNNNN_2_hg38+ds"):
            (root / folder).mkdir(parents=True)
        (root / "genome_library/NNNNNNNN_2_hg38+ds/X_chr1_1.bin").write_bytes(
            struct.pack("=ii", 7, 20))
        err = io.StringIO()
        with mock.patch.object(sys, "stderr", err), self.assertRaises(SystemExit):
            tbc.crisprme._refuse_no_variant_index([str(root / "VCFs" / "ds")], "hg38", "TTTV",
                                                  1, 1, False, False, 20, "_")
        said = AFlankWiderThanA5PrimePamIsSaid.SAID.strip().replace("_hg38", "_hg38+ds")
        self.assertIn(said + ". ", err.getvalue())


class ASearchTheInputCheckRefusesLeavesItsOutputFolderEmpty(unittest.TestCase):
    """A search the input check refused (a guide without its PAM's N padding, say) left
    .command_line.txt, .version.txt, Params.txt and guides.txt in its output folder, so the
    corrected command, with the same --output, stopped on "Output folder ... not empty!"
    (found running this track's chr22 searches). It leaves the folder as it found it."""

    def test_the_same_command_again_is_refused_for_its_input_again(self):
        import subprocess
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        (tmp / "Genomes" / "g").mkdir(parents=True)
        (tmp / "Genomes" / "g" / "chr1.fa").write_text(">chr1\n" + "ACGTTGCAAC" * 60 + "\n")
        (tmp / "PAMs").mkdir()
        (tmp / "PAMs" / "20bp-NGG-SpCas9.txt").write_text("N" * 21 + "GG 3\n")
        (tmp / "g.txt").write_text("ACTGAAATCTGTAAGCAGGC\n")  # without its NNN
        env = {**os.environ, "CRISPRME_SKIP_SCORER_ENV": "1"}
        try:
            import Bio.Seq  # noqa: F401
        except ImportError:  # (the search reads nothing of it before its input check)
            (tmp / "stub" / "Bio").mkdir(parents=True)
            (tmp / "stub" / "Bio" / "__init__.py").write_text("")
            (tmp / "stub" / "Bio" / "Seq.py").write_text("Seq = str\n")
            env["PYTHONPATH"] = str(tmp / "stub")
        for attempt in (1, 2):
            run = subprocess.run(
                [sys.executable, str(_REPO / "crisprme.py"), "complete-search", "--genome",
                 "Genomes/g", "--guide", "g.txt", "--pam", "PAMs/20bp-NGG-SpCas9.txt", "--mm",
                 "4", "--output", "o", "--thread", "1"],
                cwd=tmp, env=env, capture_output=True, text=True, timeout=300)
            with self.subTest(attempt=attempt):
                self.assertEqual(run.returncode, 1, run.stderr[-2000:])
                self.assertIn("expects a padded length of 23 nt", run.stdout + run.stderr)
                self.assertNotIn("not empty", run.stderr)
                self.assertEqual(os.listdir(tmp / "Results" / "o"), [])


class PoolSearchIndelsKeepsNoDeadHelper(unittest.TestCase):
    """pool_search_indels.py kept _index_pam, which nothing called (emptyidx F8)."""

    def test_every_function_of_it_is_called(self):
        import re
        src = (_HERE / "pool_search_indels.py").read_text()
        for node in ast.parse(src).body:
            if isinstance(node, ast.FunctionDef):
                with self.subTest(function=node.name):
                    self.assertGreater(len(re.findall(rf"\b{node.name}\b", src)), 1)


class TheDesignDocSaysWhatShips(unittest.TestCase):
    """DYNAMIC_PAM_INDEX_COMPAT.md called a search that prunes by itself, the web's pruning
    and the PAM picker *planned*, and its CLI section described what the CLI did not do
    (IDX-22): 2.9.0 does them, and the sections say so, and how it differs."""

    def test_the_shipped_sections_are_not_planned(self):
        text = (_REPO / "docs" / "DYNAMIC_PAM_INDEX_COMPAT.md").read_text()
        for heading in ("## 4. CLI", "## 5. Web", "## 8. Web PAM picker"):
            section = text[text.index(heading):].split("\n## ", 1)[0]
            with self.subTest(heading=heading):
                self.assertNotIn("*planned*", section)
                self.assertIn("2.9.0", section)


class OnePerfectMatchPerGuideIsNoWarning(unittest.TestCase):
    """The web results page put up its red "Multiple perfect matches" banner for a search
    of two guides with one perfect match each: it counted the sites of every guide together
    (RV1-F5). As the report's banner does, it counts them per guide."""

    HEADER = ["Spacer+PAM", "Chromosome", "Start_coordinate_(highest_CFD)",
              "Strand_(highest_CFD)", "Mismatches+bulges_(highest_CFD)"]

    def alert(self, rows, count_guides):
        import types
        tsv = Path(tempfile.mkdtemp()) / "o.integrated_results.tsv"
        self.addCleanup(shutil.rmtree, tsv.parent, ignore_errors=True)
        tsv.write_text("\n".join("\t".join(r) for r in [self.HEADER] + rows) + "\n")
        dash = types.SimpleNamespace(
            Alert=lambda children, color: (color, "".join(map(str, children))),
            Strong=lambda text: text, Code=lambda text: text)
        ns = _web_functions("_perfect_match_sites", "_perfect_match_alert",
                            page="results_page.py", dbc=dash, html=dash)
        return ns["_perfect_match_alert"](str(tsv), count_guides)

    def test_one_perfect_match_per_guide(self):
        rows = [["GAGTCCGAGCAGAAGAAGAANNN", "chr1", "100", "+", "0"],
                ["CTAACAGTTGCTTTTATCACNNN", "chr2", "200", "-", "0"],
                ["CTAACAGTTGCTTTTATCACNNN", "chr2", "900", "-", "2"]]
        self.assertIsNone(self.alert(rows, 2))

    def test_a_guide_with_two_is_named(self):
        rows = [["GAGTCCGAGCAGAAGAAGAANNN", "chr1", "100", "+", "0"],
                ["CTAACAGTTGCTTTTATCACNNN", "chr2", "200", "-", "0"],
                ["CTAACAGTTGCTTTTATCACNNN", "chr5", "300", "+", "0"]]
        color, text = self.alert(rows, 2)
        self.assertEqual(color, "danger")
        self.assertIn("Guide CTAACAGTTGCTTTTATCACNNN matches 2 genomic sites with 0 "
                      "mismatches and 0 bulges (chr2:200, chr5:300).", text)
        self.assertNotIn("GAGTCCGAGCAGAAGAAGAANNN", text)

    def test_one_guide_with_two_as_before(self):
        rows = [["CTAACAGTTGCTTTTATCACNNN", "chr2", "200", "-", "0"],
                ["CTAACAGTTGCTTTTATCACNNN", "chr5", "300", "+", "0"]]
        color, text = self.alert(rows, 1)
        self.assertIn("This guide matches 2 genomic sites with 0 mismatches and 0 bulges "
                      "(chr2:200, chr5:300). Each is a candidate cut site", text)
        self.assertIsNone(self.alert(rows[:1], 1))


if __name__ == "__main__":
    unittest.main()
