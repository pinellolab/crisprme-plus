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


def _web_functions(*names, **namespace):
    """The functions ``names`` of pages/main_page.py, run in ``namespace`` (their data
    lookups stubbed) with the module-level assignments they read."""
    src = (_REPO / "pages" / "main_page.py").read_text()
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
        # with no bulges, no index: the genome is scanned
        run = self.run_with_params(pam="20bp-NRR-SpCas9", bulges=0, variant=False)
        self.assertIn("Reference_index\tnone: a search with no bulges scans the genome\n",
                      self.params())
        self.assertNotIn("Prune_pam", self.params())


if __name__ == "__main__":
    unittest.main()
