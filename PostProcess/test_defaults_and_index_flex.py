#!/usr/bin/env python3
"""The default PAM and index, and which index serves a search (CRISPRme 2.10.0).

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
import re
import shutil
import subprocess
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
                # an install from before 2.10.0: what it selected then
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

    def test_no_surface_calls_the_nrg_index_the_default(self):
        # the build recipes' notes and the paper's table still named NRG the default
        for path, gone in (
                ("seq_script/merge_panels/README.md",
                 "default is the re-phased `NRG_3_hg38+hg38_HGDP1kGP`"),
                ("seq_script/merge_panels/README.md",
                 "The current shipped default (`NRG_3_hg38+hg38_HGDP1kGP`)"),
                ("seq_script/merge_panels/README.md",
                 "`NRG_3_hg38+hg38_HGDP1kGP`** (4,091 samples = 925 HGDP + 3,166 1000 Genomes), "
                 "which is the current **recommended variant default**"),
                ("seq_script/merge_panels/hgdp1kgp_build.sh",
                 'the DEFAULT "real-haplotypes" variant index NRG_3_hg38+hg38_HGDP1kGP'),
                ("seq_script/merge_panels/onekg2021_build.sh",
                 "default is the SHAPEIT5 re-phased NRG_3_hg38+hg38_HGDP1kGP"),
                ("paper/reproducibility/variant-datasets/README.md",
                 "`NRG_3_hg38+hg38_HGDP1kGP` *(default)*"),
                ("test/web/test_variant_dropdown.py",
                 "`NRG_3_hg38+hg38_HGDP1kGP`, the current default")):
            with self.subTest(path=path, gone=gone[:40]):
                # (the text of the lines, without a comment's "#" or a quote's ">")
                text = " ".join(line.lstrip("#> \t")
                                for line in (_REPO / path).read_text().splitlines())
                self.assertNotIn(gone, " ".join(text.split()))

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
        self.assertIn("; if (`$LASTEXITCODE -eq 0) { Remove-Item -Force -LiteralPath `\"$mark`\"; "
                      "$run;", line[index_step:])

    def test_a_start_after_a_failed_download_sets_up_again(self):
        # Start chose between setting up and opening the app by Genomes/ alone, which the
        # download makes in its first minutes: a Start after a failed download opened the app
        # on part of the data, as the failure's own advice to click Start again led to. The
        # setup marks the data folder until both downloads have ended well. (Checked by hand
        # too: the app's AppleScript with osascript, install.ps1 with PowerShell 7, each with a
        # stub docker that fails the index download, then not.)
        ps1 = (_REPO / "docs" / "install.ps1").read_text(encoding="utf-8")
        self.assertIn("function Has-Data { (Test-Path (Join-Path $Data 'Genomes')) -and -not "
                      "(Test-Path (Setup-Mark)) }", ps1)
        start = ps1[ps1.index("$mark = Setup-Mark"):]
        self.assertLess(start.index("New-Item -ItemType File -Force -Path $mark"),
                        start.index("Run-InConsole"))
        sh = (_REPO / "docs" / "install.sh").read_text(encoding="utf-8")
        self.assertIn('" >/dev/null 2>&1 && ! test -e " & quoted form of my setupMark() & '
                      '" && echo yes || echo no") is "yes"\nend hasData', sh)
        start = sh[sh.index("set mark to quoted form of my setupMark()"):]
        self.assertLess(start.index('do shell script "touch " & mark'),
                        start.index("runInTerminal("))
        self.assertIn('runInTerminal("if " & dlRef & " && " & dlVar & "; then rm -f " & mark & '
                      '"; " & runCmd', start)
        self.assertIn("else echo; echo 'The download did not finish: see above, then open "
                      "CRISPRme and click Start again to finish setting up.'; fi", start)


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


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class ASearchLeavesNoViewOfItsIndexes(tir._ShippedScript):
    """The engine searches an index folder that holds .bin files alone, so a search hands it
    a view of an index with a sidecar (.pam_build: a download, and every index a search or
    build-index-only builds): links to its .bin files, in <output>/.binview/. The view stayed
    in the result folder (RV3-4, RV3-6): it goes as the search ends, found what it may or
    failed."""

    def sidecars(self, *names):
        for name in names:
            (self.library / name / ".pam_build").write_text("20bp-NRG-SpCas9 23 3\n")

    def test_a_search_that_goes_on(self):
        self.install("NRG_3", variant=False)
        self.sidecars("NRG_3_hg38")
        run = self.run_pipeline(variant=False, env={"STUB_FINDS_NOTHING": "1"})
        self.assertIn("Search Reference completed", run.stdout)
        searched = next(c for c in run.calls if c[0] == "search")[1]
        self.assertEqual(Path(searched).parent.name, ".binview")  # (it was given the view)
        self.assertFalse((self.root / f"run{self.runs}" / "out" / ".binview").exists())

    def test_searches_that_fail(self):
        self.install("NRG_3")
        self.sidecars("NRG_3_hg38", f"NRG_3_hg38+{self.VCF}")
        for last in ("reference", "variant"):  # (the stub fails the run's last search)
            with self.subTest(last=last):
                run = self.run_pipeline(variant=last == "variant")
                self.assertIn(f"ERROR: Off-targets search {last.title()} failed", run.stderr)
                searched = [c for c in run.calls if c[0] == "search"][-1][1]
                self.assertEqual(Path(searched).parent.name, ".binview")
                self.assertFalse((self.root / f"run{self.runs}" / "out" / ".binview").exists())


def index_fit_pam(name):
    """The PAM of index ``name``, either form of the name (index_fit.index_pam)."""
    import index_fit
    return index_fit.index_pam(name)


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


class TheReportNamesTheIndexesSearched(unittest.TestCase):
    """The report's summary names the indexes the search took, which the search script
    records in .indexes_searched.txt beside the run's parameters, and the report's zip
    bundles that record beside Params.txt, which it leaves as the run wrote it."""

    def test_from_the_record_beside_the_parameters(self):
        import zipfile
        import test_generate_report as tgr
        gr = tgr.gr
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        tsv = tmp / f"{tgr._GUIDE}+NRG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv"
        tsv.write_text("\t".join(tgr._HEADER) + "\n"
                       + "".join("\t".join(row) + "\n" for row in tgr._ROWS))
        params = "Genome_selected\thg38\nGenome_idx\tNGG_2_hg38\nPam\tNGG\n"
        (tmp / ".Params.txt").write_text(params)
        (tmp / ".indexes_searched.txt").write_text(
            "Reference_index\tNRR_3_hg38 (in /lib): built for the PAM NRR on the 3' side of "
            "a 20-nt guide, with up to 2 bulge(s) of each type, pruned to NGG at leaf time "
            "(--prune-pam)\nPrune_pam\tauto\n")
        out_zip = gr.build_report(result_dir=str(tmp), integrated_tsv=str(tsv),
                                  out_zip=str(tmp / "r.zip"),
                                  samplesid_dir=str(tmp / "samplesIDs"))
        with zipfile.ZipFile(out_zip) as zf:
            html = zf.read("report.html").decode()
            self.assertEqual(zf.read("data/Params.txt").decode(), params)
            self.assertEqual(zf.read("data/indexes_searched.txt").decode(),
                             (tmp / ".indexes_searched.txt").read_text())
        self.assertIn("<td class=\"k\">Reference index searched</td><td>NRR_3_hg38 (in /lib): "
                      "built for the PAM NRR", html)
        self.assertIn("<td class=\"k\">PAM pruning</td>", html)


class ThePanelIsTheSearchsOwn(unittest.TestCase):
    """The report's 'Panel: N individuals' counted every list in samplesIDs/, the lists of
    the datasets not searched and those Settings writes too (M2): it counts the samples the
    search ran with, which its script gathers in .sampleID.txt, each under its own dataset."""

    def test_of_the_samples_the_search_ran_with(self):
        import zipfile
        import test_generate_report as tgr
        gr = tgr.gr
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        tsv = tmp / f"{tgr._GUIDE}+NRG_hg38+hg38_1000G_6+3_integrated_results.tsv"
        tsv.write_text("\t".join(tgr._HEADER) + "\n"
                       + "".join("\t".join(row) + "\n" for row in tgr._ROWS))
        sids = tmp / "samplesIDs"
        sids.mkdir()
        head = "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
        (sids / "hg38_1000G.samplesID.txt").write_text(
            head + "HG00096\tGBR\tEUR\tmale\nHG00097\tGBR\tEUR\tfemale\n"
            "NA18525\tCHB\tEAS\tfemale\n")
        (sids / "hg38_HGDP.samplesID.txt").write_text(  # (a dataset not searched)
            head + "HGDP00001\tBrahui\tCSA\tmale\nHGDP00003\tBrahui\tCSA\tmale\n")
        (tmp / ".sampleID.txt").write_text(  # (what the search gathered: 1000G's list)
            head + "HG00096\tGBR\tEUR\tmale\nHG00097\tGBR\tEUR\tfemale\n"
            "NA18525\tCHB\tEAS\tfemale\n")
        with zipfile.ZipFile(gr.build_report(result_dir=str(tmp), integrated_tsv=str(tsv),
                                             out_zip=str(tmp / "r.zip"),
                                             samplesid_dir=str(sids))) as zf:
            html = zf.read("report.html").decode()
        self.assertIn(" Panel: <strong>3</strong> individuals.", html)
        self.assertEqual(gr.run_sample_ids(str(tmp)), {"HG00096", "HG00097", "NA18525"})
        # a report of a table alone, with no result folder, as before
        self.assertIsNone(gr.run_sample_ids(None))


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class ASearchPrunesByItself(tir._ShippedScript):
    """When no installed index serves a search unpruned and one serves it pruned (a wider
    all-N one, a dedicated one of a PAM that holds the search's), the search prunes as if
    --prune-pam was given, and says so first. It stopped for want of a variant index (the
    flag told to the user), or built a reference one. The web launches the script with
    pruning off, so its searches do the same. The run records the indexes the search took,
    for the report, beside its parameters file: the web reuses a job only when its
    .Params.txt is a new request's, byte for byte."""

    def params(self):
        """The record of the indexes the run searched, once its parameters are checked to
        be as written."""
        out = self.root / f"run{self.runs}" / "out"
        self.assertEqual((out / "Params.txt").read_text(), "Genome_selected\thg38\n")
        return (out / ".indexes_searched.txt").read_text()

    def run_with_params(self, **kw):
        # (complete-search writes the parameters file before it launches the script)
        run = self.launch(**kw)
        (run.work / "out" / "Params.txt").write_text("Genome_selected\thg38\n")
        return self.finish(run)

    def test_an_nrr_variant_index_serves_an_ngg_search_pruned(self):
        self.install("NRR_3")
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0)
        self.assertIn(f"NOTE: no installed variant index of hg38+{self.VCF} serves PAM NGG "
                      f"unpruned, and NRR_3_hg38+{self.VCF} serves it pruned to it at leaf "
                      "time: searching with --prune-pam, as if it was given", run.stdout)
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

    def test_the_note_names_the_index_that_needs_it(self):
        # the variant index is the PAM's own, the reference one another PAM's: it said that
        # no installed index served the PAM unpruned, and to build a dedicated one
        self.install("NGG_3", reference=False)
        self.install("NRR_3", variant=False)
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0)
        self.assertIn("NOTE: no installed reference index of hg38 serves PAM NGG unpruned, and "
                      "NRR_3_hg38 serves it pruned to it at leaf time", run.stdout)
        self.assertEqual(run.searched(variant=False), ("NRR_3_hg38", True))
        # complete-search's, the same, from the same library (crisprme.py _prune_on_its_own)
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        vcfs = self.home / "VCFs" / self.VCF
        vcfs.mkdir(parents=True, exist_ok=True)
        with mock.patch.object(cm, "current_working_directory", str(self.home)), \
                mock.patch.object(cm, "_engine_prunes", return_value=True):
            for bulges, which in ((1, "reference index of hg38"), (0, None)):
                index = cm._prune_on_its_own("_", "hg38", [str(vcfs)], "NGG", bulges + 1,
                                             True, 20, bulges, bulges)
                with self.subTest(bulges=bulges):
                    if which is None:  # (no bulges: the reference index is not opened)
                        self.assertIsNone(index)
                        continue
                    self.assertEqual(Path(index).name, "NRR_3_hg38")
                    self.assertTrue(cm._prune_note(index, "NGG").startswith(
                        f"NOTE [complete-search]: no installed {which} serves PAM NGG "
                        f"unpruned, and NRR_3_hg38 ({self.library}) serves it pruned"))
            shutil.rmtree(self.library / f"NGG_3_hg38+{self.VCF}")
            self.install("NRR_3", reference=False)
            index = cm._prune_on_its_own("_", "hg38", [str(vcfs)], "NGG", 2, True, 20, 1, 1)
            self.assertIn(f"no installed variant index of hg38+{self.VCF} serves PAM NGG",
                          cm._prune_note(index, "NGG"))

    def test_a_wider_pam_less_index_alone_serves_pruned(self):
        self.install("NNNNNNNN_3")
        run = self.run_pipeline(pam="20bp-NRR-SpCas9", prune=0)
        self.assertEqual(run.searched(variant=True), (f"NNNNNNNN_3_hg38+{self.VCF}", True))
        self.assertIn("searching with --prune-pam, as if it was given", run.stdout)

    def test_it_says_what_a_flank_wider_than_the_pam_costs_where_it_decides(self):
        """complete-search says, of an all-N index wider than the PAM that a search reads
        pruned, what such a flank cannot index, and warns of one built before CRISPRitz 2.9.1
        (crisprme.py _wider_flank_notes), on its console. A web search, which prunes by
        itself since 2.10.0, said neither: the script says both in its log, where it decides.
        Given "auto", complete-search has said them."""
        import test_prune_flank as tpf
        (self.home / "Genomes" / "hg38" / "chr1.fa").write_text(
            ">chr1\n" + "ACGT" * 20 + "N" * 30 + "ACGT" * 20 + "\n")
        self.install("NNNNNNNN_3")
        # (the trees of chr1, whose sequence holds an N: none in a flank, as an index built
        # before the fix; the variant index's INDELS index tells nothing)
        for name in ("NNNNNNNN_3_hg38", f"NNNNNNNN_3_hg38+{self.VCF}"):
            tpf.tree(str(self.library / name / "NNNNNNNN_chr1_1.bin"), 8,
                     ["ACGTACGT", "TTTTACGT"])
        flank = ("NOTE: --prune-pam serves the 3-nt PAM NGG from {}, whose stored flank is 8 "
                 "nt: a + strand PAM with 1 to 5 nt after it to the end of a contig, and a - "
                 "strand PAM in the first 5 nt of a contig, are not indexed")
        old = ("WARNING: {} stores no N in the flank of any target where its genome has one: "
               "it was built before CRISPRitz 2.9.1 (CRISPRme 2.8.0 or earlier), which "
               "dropped such targets, and searched for the shorter PAM NGG it misses the "
               "sites whose flank holds an N. Put a dedicated NGG index beside it, which a "
               "search takes before it (build one with crisprme.py build-index-only or in "
               "the web app's Settings, or download one if it is published), or remove it and "
               "build it again with build-index-only (which keeps an index it finds in "
               "place).")
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0)
        # (nothing on stderr, which fails a stage, but the stub's stop at the last search)
        self.assertEqual(run.stderr, "ERROR: Off-targets search Variant failed\n")
        self.assertIn("searching with --prune-pam, as if it was given", run.stdout)
        # (the flank once, of the first index; the warning of each)
        self.assertEqual(run.stdout.count("whose stored flank is 8 nt"), 1)
        self.assertIn(flank.format("NNNNNNNN_3_hg38"), run.stdout)
        for name in ("NNNNNNNN_3_hg38", f"NNNNNNNN_3_hg38+{self.VCF}"):
            self.assertIn(old.format(name), run.stdout)
        # a reference search reads the reference index alone
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0, variant=False)
        self.assertIn(old.format("NNNNNNNN_3_hg38"), run.stdout)
        self.assertNotIn(f"+{self.VCF} stores no N", run.stdout)
        # given the flag, or "auto", complete-search has said them on its console
        for passed in ("auto", 1):
            with self.subTest(prune=passed):
                run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=passed, index_path="_")
                self.assertEqual(run.searched(variant=True),
                                 (f"NNNNNNNN_3_hg38+{self.VCF}", True))
                self.assertNotIn("stored flank", run.stdout)
                self.assertNotIn("stores no N", run.stdout)
        # an index of the PAM's width serves it unpruned: nothing to say
        shutil.rmtree(self.library)
        self.install("NNN_3")
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=0)
        self.assertNotIn("as if it was given", run.stdout)
        self.assertNotIn("stored flank", run.stdout)

    def test_a_reference_search_takes_one_rather_than_build_one(self):
        self.install("NRR_3", variant=False)
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", variant=False, prune=0)
        self.assertEqual(run.builds, [])
        self.assertEqual(run.searched(variant=False), ("NRR_3_hg38", True))

    def test_the_run_records_the_indexes_searched(self):
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

    def test_a_search_with_no_bulges_records_the_variant_genome_it_scans(self):
        """Its variant genome here whole, a search with no bulges scans it (STEP 3) and
        opens no variant index: the record named one installed of its PAM all the same, with
        the bulges it holds and, where the search prunes by itself, pruned at leaf time, and
        the report's summary named it as the index searched. Only its INDELS index is."""
        enriched = self.home / "Genomes" / f"hg38+{self.VCF}"
        enriched.mkdir(parents=True)
        (enriched / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        self.install("NRR_3")
        for pam, pruned in (("20bp-NRR-SpCas9", False), ("20bp-NGG-SpCas9", True)):
            with self.subTest(pam=pam):
                run = self.run_with_params(pam=pam, bulges=0)
                variant_search = [c for c in run.calls if c[0] == "search"][-1]
                # (the enriched genome, scanned: -r, no -index)
                self.assertEqual((Path(variant_search[1]).name, "-r" in variant_search),
                                 (f"hg38+{self.VCF}", True))
                params = self.params()
                self.assertIn(f"Variant_index:{self.VCF}\tnone: a search with no bulges scans "
                              f"the variant genome; its indels: NRR_3_hg38+{self.VCF}_INDELS\n",
                              params)
                self.assertNotIn(f"NRR_3_hg38+{self.VCF} (in", params)
                # the INDELS index, which the indel search reads, is pruned as the search is
                self.assertEqual("Prune_pam\tauto\n" in params, pruned)
        # without that genome, the variant index answers the search, and the record names it
        shutil.rmtree(enriched)
        run = self.run_with_params(pam="20bp-NRR-SpCas9", bulges=0)
        self.assertEqual(Path([c for c in run.calls if c[0] == "search"][-1][1]).name,
                         f"NRR_3_hg38+{self.VCF}")
        self.assertIn(f"Variant_index:{self.VCF}\tNRR_3_hg38+{self.VCF} (in {self.library}): "
                      "built for the PAM NRR on the 3' side", self.params())

    def test_the_record_says_a_search_pruned_by_itself_wherever_it_was_decided(self):
        """complete-search decides the pruning for a search of its own and passes the flag,
        so the script's own decision never ran and the record held no Prune_pam line: a
        command-line search and a web search of the same library left different records
        (METHODS section 9 lists the line). complete-search says it took the flag itself
        ("auto"), and the record says so for either. Given "auto", the script takes the
        flag without deciding again, and without saying it again: complete-search has said
        so on its console."""
        self.install("NRR_3", variant=False)
        # (0 is a web search, which the script decides; with no --index-path, as from the
        # command line)
        for passed, recorded, said in ((0, True, True), ("auto", True, False),
                                       (1, False, False)):
            with self.subTest(prune=passed):
                run = self.run_with_params(pam="20bp-NGG-SpCas9", variant=False, prune=passed,
                                           index_path="_")
                self.assertEqual(run.searched(variant=False), ("NRR_3_hg38", True))
                self.assertEqual("Prune_pam\tauto\n" in self.params(), recorded)
                self.assertEqual("as if it was given" in run.stdout, said)


class CompleteSearchSaysWhoTookThePruningFlag(unittest.TestCase):
    """complete-search decides the pruning for its own search (_prune_on_its_own) and hands
    the pipeline --prune-pam, so the script's own decision never ran and its record of the
    indexes searched held no `Prune_pam auto` line, which METHODS section 9 lists and a web
    search of the same library wrote. It tells the script which it was: "auto" where it took
    the flag itself, "1" where it was given, "0" where the search prunes nothing."""

    def launch(self, *extra, installed=("NRR_3_g",)):
        """(the last argument of the pipeline command line complete-search builds for a
        reference search of NGG with one bulge of each type, what it said on stdout), with
        ``installed`` in the working directory's library."""
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        (tmp / "Genomes" / "g").mkdir(parents=True)
        (tmp / "Genomes" / "g" / "chr1.fa").write_text(">chr1\n" + "ACGTTGCAAC" * 60 + "\n")
        (tmp / "PAMs").mkdir()
        shutil.copy(_REPO / "PAMs" / "20bp-NGG-SpCas9.txt", tmp / "PAMs")
        (tmp / "g.txt").write_text("ACTGAAATCTGTAAGCAGGC" + "N" * 3 + "\n")
        _lib(tmp, *installed)
        argv = ["crisprme.py", "complete-search", "--genome", "Genomes/g", "--guide", "g.txt",
                "--pam", "PAMs/20bp-NGG-SpCas9.txt", "--mm", "4", "--bDNA", "1", "--bRNA", "1",
                "--output", "o", "--thread", "1", *extra]
        said, launched = io.StringIO(), []
        here = os.getcwd()
        self.addCleanup(os.chdir, here)
        os.chdir(tmp)
        with mock.patch.object(cm, "input_args", argv), \
                mock.patch.object(cm, "current_working_directory", f"{tmp}/"), \
                mock.patch.object(cm, "warn_low_memory", lambda: None), \
                mock.patch.object(cm, "_engine_prunes", return_value=True), \
                mock.patch.object(cm.subprocess, "call",
                                  side_effect=lambda cmd, **kw: launched.append(cmd) or 0), \
                mock.patch.dict(os.environ, {"CRISPRME_SKIP_SCORER_ENV": "1",
                                             "CRISPRME_SKIP_MEM_CHECK": "1"}), \
                redirect_stdout(said):
            cm.complete_search()
        self.assertEqual(len(launched), 1, said.getvalue())
        self.params = (Path(tmp) / "Results" / "o" / ".Params.txt").read_text()
        return launched[0].split()[-1], said.getvalue()

    def test_auto_where_it_took_the_flag_itself(self):
        passed, said = self.launch()
        self.assertEqual(passed, "auto")
        self.assertIn("searching with --prune-pam, as if it was given", said)
        self.assertIn("Prune_pam\tauto\n", self.params)

    def test_one_where_it_was_given_and_none_where_nothing_prunes(self):
        passed, said = self.launch("--prune-pam")
        self.assertEqual(passed, "1")
        self.assertNotIn("as if it was given", said)
        self.assertIn("Prune_pam\tgiven\n", self.params)
        # the PAM's own index serves it unpruned: no flag, and nothing said
        passed, said = self.launch(installed=("NGG_3_g", "NRR_3_g"))
        self.assertEqual(passed, "0")
        self.assertNotIn("as if it was given", said)
        self.assertNotIn("Prune_pam", self.params)


@unittest.skipUnless(shutil.which("python") and tir._gnu_realpath(),
                     "runs the shipped search script as the image does (GNU realpath)")
class AnIndexNameMayCarryItsGuideAndSide(unittest.TestCase):
    """A name of the 2.10.0 form carries the guide (spacer) length and the side of the PAM, in
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

    def test_the_side_a_name_says_is_its_record_for_every_reader(self):
        """A name that says the side records it as a .pam_build does, ahead of one: the
        search script refused, for want of a .pam_build, a PAM-less index named for a 5'
        PAM (NNNN-20_2_hg38) for a 5' search, and one named for a 3' PAM of a width but 3
        and 8 (20-NNNN_2_hg38) for a 3' search, which crisprme.py's pre-flight and the web
        offered. Old and new names, 3' and 5', recorded, unrecorded and recorded otherwise:
        the search script, crisprme.py, the web and complete-test take the same."""
        import test_bulge_cap as tbc
        import complete_test
        cm, pu = tbc.crisprme, tbc.pages_utils
        pams = {"20bp-TTTV-Cas12a": ("TTTV" + "N" * 20, -4),
                "20bp-NNNN-NO-PAM5": ("N" * 24, -4),
                "20bp-NGG-SpCas9": ("N" * 20 + "NGG", 3),
                "20bp-NRRT-x": ("N" * 20 + "NRRT", 4)}
        indexes = ("NNNN-20_2_hg38", "20-NNNN_2_hg38", "TTTV-20_2_hg38", "20-TTTV_2_hg38",
                   "20-NNNNNNNN_2_hg38", "NNNN_2_hg38", "TTTV_2_hg38", "NNNNNNNN_2_hg38")
        taken = {}
        for name in indexes:
            k = len(index_fit_pam(name))
            for record in (None, 3, 5):
                root = Path(tempfile.mkdtemp(dir=self.tmp.name))
                lib = _lib(root, name)
                if record:
                    (lib / name / ".pam_build").write_text(
                        f"x {20 + k} {k if record == 3 else -k}\n")
                (root / "PAMs").mkdir()
                for pam, (seq, pos) in pams.items():
                    (root / "PAMs" / f"{pam}.txt").write_text(f"{seq} {pos}\n")
                for pam, (seq, pos) in pams.items():
                    motif = seq[:-pos] if pos < 0 else seq[-pos:]
                    for prune in (0, 1):
                        shell = self.take(lib, motif, prune=prune, pos=pos)
                        with self.subTest(index=name, record=record, pam=pam, prune=prune):
                            mirror = cm._index_the_search_takes(
                                [str(lib)], motif, "hg38", 2, bool(prune), pos > 0, 20, 1, 1,
                                False)
                            self.assertEqual(Path(mirror).name if mirror else "", shell)
                            if prune:  # (the web counts what the search prunes by itself)
                                with mock.patch.object(pu, "current_working_directory",
                                                       str(root) + os.sep):
                                    self.assertEqual(pu.index_max_bulges("hg38", pam),
                                                     1 if shell else 0)
                        taken[name, record, pam, prune] = shell
                with self.subTest(index=name, record=record, reader="complete-test"):
                    side = 5 if name.split("_")[0].endswith("-20") else (
                        3 if name.startswith("20-") else record)
                    with mock.patch.object(os, "getcwd", return_value=str(root)):
                        five = complete_test._smoke_five_prime(name, index_fit_pam(name))
                    if side:
                        self.assertEqual(five, side == 5)
        # what the reviewers' runs found refused, served: a 5' PAM-less index named for its
        # side, with no .pam_build, for a TTTV search; a 3' one of width 4 for an NGG search
        self.assertEqual(taken["NNNN-20_2_hg38", None, "20bp-TTTV-Cas12a", 0], "NNNN-20_2_hg38")
        self.assertEqual(taken["20-NNNN_2_hg38", None, "20bp-NGG-SpCas9", 1], "20-NNNN_2_hg38")
        # ... the name ahead of a record that says the other side
        self.assertEqual(taken["NNNN-20_2_hg38", 3, "20bp-NNNN-NO-PAM5", 0], "NNNN-20_2_hg38")
        self.assertEqual(taken["20-NNNN_2_hg38", 5, "20bp-NGG-SpCas9", 1], "20-NNNN_2_hg38")
        # ... and never for the other side
        self.assertEqual(taken["NNNN-20_2_hg38", 5, "20bp-NGG-SpCas9", 1], "")
        self.assertEqual(taken["20-TTTV_2_hg38", None, "20bp-TTTV-Cas12a", 1], "")
        # an older name, unrecorded, as before: a PAM-less one of width 4 is not taken for a
        # 3' PAM, a recorded one is
        self.assertEqual(taken["NNNN_2_hg38", None, "20bp-NGG-SpCas9", 1], "")
        self.assertEqual(taken["NNNN_2_hg38", 3, "20bp-NGG-SpCas9", 1], "NNNN_2_hg38")
        self.assertEqual(taken["NNNN_2_hg38", 5, "20bp-TTTV-Cas12a", 0], "NNNN_2_hg38")

    def test_index_side_reads_the_name_first(self):
        for name, side in (("NNNN-20_2_hg38", "5"), ("20-NNNN_2_hg38", "3"),
                           ("NNNN_2_hg38", "3")):
            with self.subTest(name=name):
                lib = _lib(tempfile.mkdtemp(dir=self.tmp.name), name)
                said = subprocess.run(
                    ["bash", "-c", tir.RESOLVER + '\n_index_side "$1"; echo "$?"', "bash",
                     str(lib / name)], capture_output=True, text=True)
                # (an older name with no record is taken for a 3' build, status 1)
                self.assertEqual(said.stdout.split(), [side, "1" if "-" not in name else "0"])


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
    """A variant index and its INDELS index under the 2.10.0 name: the scans take them, the
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

    def test_a_name_that_says_its_side_needs_no_record(self):
        # a PAM-less index named for a 5' PAM, with no .pam_build, for a 5' search under
        # --index-path: it was refused, as "not taken for a 3' PAM", and it records its side
        _lib(self.home, "NNNN-20_2_hg38")
        tttv = "TTTV" + "N" * 20 + " -4\n"
        run = self.run_pipeline(pam="20bp-TTTV-Cas12a", pam_text=tttv, variant=False,
                                index_path=self.library)
        self.assertNotIn("index", run.stderr)  # (the stub stops the run at its search)
        self.assertEqual(run.searched(variant=False), ("NNNN-20_2_hg38", False))
        self.assertIn("NNNN-20_2_hg38 (in", self.params())
        self.assertIn("built for the PAM NNNN on the 5' side", self.params())
        # one named for the other side stays refused, and the refusal says why by its name
        shutil.rmtree(self.library / "NNNN-20_2_hg38")
        _lib(self.home, "20-NNNN_2_hg38")
        run = self.run_pipeline(pam="20bp-TTTV-Cas12a", pam_text=tttv, variant=False,
                                index_path=self.library)
        self.assertEqual(run.searched(variant=False), None)
        self.assertIn("20-NNNN_2_hg38 is a PAM-less index for a PAM on the 3' side of the "
                      "guide (as its name says), and 20bp-TTTV-Cas12a.txt puts the PAM on the "
                      "5' side", run.stderr)
        self.assertNotIn("records no side", run.stderr)
        self.assertNotIn("record that by running", run.stderr)

    def test_an_unrecorded_pam_less_index_is_named_and_how_to_use_it(self):
        # an all-N index of a width but 3 and 8, with no .pam_build, really built for a 3'
        # PAM, under --index-path: a search without --prune-pam named it not at all, and
        # with it said only that it is not taken for a 3' PAM; build-index-only records no
        # side over it, so it could only be rebuilt (R6B RV1-6). Its name can say the side
        shutil.rmtree(self.library)
        _lib(self.home, "NNNN_3_hg38")
        with open(self.library / "NNNN_3_hg38" / "chr1.bin", "wb") as fh:  # (its header)
            fh.write((7).to_bytes(4, sys.byteorder) + (20).to_bytes(4, sys.byteorder))
        for prune in (0, 1):
            with self.subTest(prune=prune):
                run = self.run_pipeline(pam="20bp-NGG-SpCas9", prune=prune, variant=False,
                                        index_path=self.library)
                self.assertEqual(run.searched(variant=False), None)
                self.assertIn("NNNN_3_hg38 is a PAM-less index that records no side of the "
                              "guide in a .pam_build", run.stderr)
                lib = Path(os.path.realpath(self.library))
                self.assertIn(f"if it was built from a 3' PAM file, name it so: mv "
                              f"'{lib / 'NNNN_3_hg38'}' '{lib / '20-NNNN_3_hg38'}', which "
                              f"records that, for a guide of 20 nt", run.stderr)
        # as it says: the search takes it, pruned by itself
        (self.library / "NNNN_3_hg38").rename(self.library / "20-NNNN_3_hg38")
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", variant=False, index_path=self.library)
        self.assertEqual(run.searched(variant=False), ("20-NNNN_3_hg38", True))

    def params(self):
        return (self.root / f"run{self.runs}" / "out" / ".indexes_searched.txt").read_text()


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


class ASearchWithNoTemporaryFolderIsRefused(unittest.TestCase):
    """A search with no temporary folder it can write (docker run --read-only without
    --tmpfs /tmp) ended with exit 0, CRISPR-Bulge scores of -1, no report and a WARNING in
    its error log alone (final confirmation of 2.8.1): complete-search stops first."""

    def test_refused_naming_the_folders_and_the_fix(self):
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        err = io.StringIO()
        with mock.patch.object(sys, "stderr", err), self.assertRaises(SystemExit) as stop:
            cm._refuse_no_temp_folder(["/nonexistent/tmp", "/nonexistent/var/tmp"])
        self.assertEqual(stop.exception.code, 1)
        self.assertEqual(err.getvalue(), (
            "Error: No temporary folder can be written here (/nonexistent/tmp, "
            "/nonexistent/var/tmp): the CRISPR-Bulge scorer and the report write their "
            "temporary files there, and without one a search reports no CRISPR-Bulge score "
            "and writes no report. Set CRISPRME_TMPDIR to a writable folder (e.g. one on the "
            "data volume, -e CRISPRME_TMPDIR=/DATA/crisprme_tmp), or give the container a "
            "/tmp (docker run --tmpfs /tmp)\n"))
        # one that can be written: nothing said, nothing left in it
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        self.assertIsNone(cm._refuse_no_temp_folder(["/nonexistent/tmp", tmp]))
        self.assertEqual(os.listdir(tmp), [])
        self.assertIn("_refuse_no_temp_folder()", (_REPO / "crisprme.py").read_text())

    def test_a_temporary_folder_set_to_the_working_one_is_taken(self):
        # docker run --read-only -v <data>:/DATA -w /DATA -e CRISPRME_TMPDIR=/DATA (or -e
        # TMPDIR=/DATA): $TMPDIR is the working folder. Its search was refused, as the
        # current folder was left out of tempfile's list wherever it stood, not only as
        # tempfile's last resort; 2.8.1 ran it to its report.
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        data = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, data, ignore_errors=True)
        gone = os.path.join(data, "nonexistent")
        real = tempfile._candidate_tempdir_list

        def read_only_container():
            # tempfile's own list, in a container whose /tmp, /var/tmp and /usr/tmp it
            # cannot write
            return [gone + d if d in ("/tmp", "/var/tmp", "/usr/tmp") else d for d in real()]

        here = os.getcwd()
        self.addCleanup(os.chdir, here)
        os.chdir(data)
        env = {k: v for k, v in os.environ.items() if k not in ("TMPDIR", "TEMP", "TMP")}
        with mock.patch.object(tempfile, "_candidate_tempdir_list", read_only_container):
            for tmpdir in (os.getcwd(), data):
                with self.subTest(TMPDIR=tmpdir), \
                        mock.patch.dict(os.environ, dict(env, TMPDIR=tmpdir), clear=True):
                    self.assertIsNone(cm._refuse_no_temp_folder())
            self.assertEqual(os.listdir(data), [])
            # with none set, the working folder is tempfile's last resort alone, which the
            # stages that write temporary files do not run in: refused, as before
            err = io.StringIO()
            with mock.patch.dict(os.environ, env, clear=True), \
                    mock.patch.object(sys, "stderr", err), \
                    self.assertRaises(SystemExit) as stop:
                cm._refuse_no_temp_folder()
        self.assertEqual(stop.exception.code, 1)
        self.assertIn(f"written here ({gone}/tmp, {gone}/var/tmp, {gone}/usr/tmp): ",
                      err.getvalue())


class TheImageIsBuiltOnAPinnedBase(unittest.TestCase):
    """The Dockerfile built on `FROM mambaorg/micromamba`, whatever image was latest that day
    (C20): it names the tag and the digest of the one 2.8.1's image was built on."""

    def test_the_base_has_a_tag_and_a_digest(self):
        import re
        froms = [line for line in (_REPO / "Dockerfile").read_text().splitlines()
                 if line.startswith("FROM ")]
        self.assertEqual(len(froms), 1)
        self.assertRegex(froms[0], r"^FROM mambaorg/micromamba:\d+\.\d+\.\d+@sha256:[0-9a-f]{64}$")


class ASelfTestOfASourceInstallSaysNoNoise(unittest.TestCase):
    """Messages a source or Bioconda install printed that the image does not (install P10):
    the Hub's "You are sending unauthenticated requests to the HF Hub" on every anonymous
    download (a server warning huggingface_hub logs, which the image's HF_HUB_VERBOSITY=error
    hides), and git's detached-HEAD advice as install_from_source.sh checked out the CRISPRitz
    tag."""

    def test_the_hub_logs_at_error_level_as_in_the_image(self):
        self.assertIn("ENV HF_HUB_VERBOSITY=error", (_REPO / "Dockerfile").read_text())
        for given, want in ((None, "error"), ("info", "info")):  # (one set stays)
            env = {k: v for k, v in os.environ.items() if k != "HF_HUB_VERBOSITY"}
            if given:
                env["HF_HUB_VERBOSITY"] = given
            with self.subTest(given=given):
                # (loading crisprme.py's definitions runs its module code)
                run = subprocess.run(
                    [sys.executable, "-c", "import os, test_bulge_cap; "
                     "print(os.environ.get('HF_HUB_VERBOSITY'))"],
                    cwd=_HERE, env=env, capture_output=True, text=True, timeout=300)
                self.assertEqual(run.stdout.split()[-1:], [want], run.stderr[-2000:])

    def test_the_crispritz_tag_is_checked_out_without_advice(self):
        lines = (_REPO / "install_from_source.sh").read_text().splitlines()
        at = next(i for i, line in enumerate(lines)
                  if line.startswith("git ") and " clone " in line)
        clone = (lines[at] + " " + lines[at + 1]).replace("\\", " ")  # (the command's 2 lines)
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        git = ["git", "-c", "user.email=t@t", "-c", "user.name=t"]
        subprocess.run([*git, "init", "-q", str(tmp / "repo")], check=True)
        subprocess.run([*git, "-C", str(tmp / "repo"), "commit", "-q", "--allow-empty", "-m",
                        "x"], check=True)
        subprocess.run([*git, "-C", str(tmp / "repo"), "tag", "v2.9.1"], check=True)
        command = (clone.replace('"${CRISPRITZ_REF}"', "v2.9.1")
                   .replace("https://github.com/pinellolab/CRISPRitz.git", f"file://{tmp}/repo")
                   .replace('"${SRC}"', str(tmp / "src")))
        run = subprocess.run(command, shell=True, capture_output=True, text=True, cwd=tmp)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertNotIn("detached HEAD", run.stderr)
        self.assertTrue((tmp / "src" / ".git").is_dir())


class TheChecklistsHpcRouteIsANewUsersRun(unittest.TestCase):
    """RELEASE_CHECKLIST's HPC route ran its apptainer commands with the host's environment
    (CONDA_* leaked in on ml008), and steps 2-4b reused step 1's writable sandbox with the
    caches it left (list_282.md, from the 2.8.1 amd64 local gate)."""

    def test_it_says_to_run_clean_and_to_rebuild_the_sandbox(self):
        text = " ".join((_REPO / "docs" / "RELEASE_CHECKLIST.md").read_text().split())
        route = text[text.index("**On an HPC node without Docker**"):]
        self.assertIn("Run the block from a clean host environment, e.g. in `env -i", route)
        self.assertIn("after step 1, build the sandbox again from the SIF", route)


class TheArchiveGateSaysWhereAnEmptyIndelsIndexComesFrom(unittest.TestCase):
    """RELEASE_CHECKLIST's gate on the published archives (T-40) said that a present-but-empty
    `_INDELS` is always a broken build. A dataset whose every indel is on a contig with an
    underscore in its name, which no build gives a sub-index (pool_index_indels.py says so),
    gets one from a build that ends well (R6B emptyidx RV1-2, RV2-2, R3-1); download refuses
    it as it refuses a broken one, so the gate's verdict, not to publish it, holds."""

    def test_it_names_both(self):
        text = " ".join((_REPO / "docs" / "RELEASE_CHECKLIST.md").read_text().split())
        gate = text[text.index("**Every published variant archive can find indels"):]
        gate = gate[:gate.index("- [ ]")]
        self.assertNotIn("always a broken build", gate)
        self.assertIn("A present-but-empty `_INDELS` is a broken build (a dataset with no indels "
                      "gets no `_INDELS` at all), or that of a dataset whose every indel is on a "
                      "contig with an underscore in its name, which no build gives a sub-index: "
                      "publish neither, as download refuses both.", gate)


class FutureWorkSaysTheShareOfIndelsInCisAsMethodsDoes(unittest.TestCase):
    """FUTURE_WORK §3 said that two or more indels in cis are common, "~17-19% of real indel
    off-targets": the raw count, mostly repeats, of which METHODS §8 gives ~1-2% of indel loci
    after repeat-masking, and ~0.1-0.2% of indel off-targets missed (plan 4.1)."""

    def test_the_raw_share_is_said_raw_beside_the_masked_one(self):
        text = " ".join((_REPO / "docs" / "FUTURE_WORK.md").read_text().split())
        self.assertIn("common in raw counts (~17–19% of real indel off-targets), most of it in "
                      "STR/VNTR repeats: ~1–2% of indel loci after repeat-masking, and ~0.1–0.2% "
                      "of indel off-targets genuinely missed (METHODS §8)", text)
        methods = " ".join((_REPO / "METHODS.md").read_text().split())
        self.assertIn("falling to **~1–2%** of indel loci after repeat-masking", methods)
        self.assertIn("**~0.1–0.2%** of indel off-targets", methods)


class AGenomeFolderWithASpaceIsRefusedUpFront(unittest.TestCase):
    """A genome folder with a space in its path failed later with a misleading path: the
    search takes it into shell commands unquoted (IDX-19), and the working directory's too.
    complete-search and build-index-only refuse either in one line, the working directory
    first, and the advice for the genome folder works: a link made under a working directory
    with no space in its path."""

    def refusal(self, command, cwd, *args):
        """What ``command`` (complete-search or build-index-only) says, run from ``cwd`` with
        ``args``, as it stops: (status, stderr)."""
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        err, here = io.StringIO(), os.getcwd()
        self.addCleanup(os.chdir, here)
        os.chdir(cwd)
        argv = ["crisprme.py", command, *args]
        with mock.patch.object(cm, "input_args", argv), \
                mock.patch.object(cm, "current_working_directory", f"{os.getcwd()}/"), \
                mock.patch.object(cm, "warn_low_memory", lambda: None), \
                mock.patch.object(sys, "stderr", err), self.assertRaises(SystemExit) as stop:
            {"complete-search": cm.complete_search, "build-index-only": cm.build_index_only}[
                command]()
        return stop.exception.code, err.getvalue()

    def test_refused_naming_it(self):
        import test_bulge_cap as tbc
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        pam = str(_REPO / "PAMs" / "20bp-NRR-SpCas9.txt")
        for cwd in (tmp / "w", tmp / "w 2"):
            for folder in ("Genomes/gx", "Genomes/g x"):
                (cwd / folder).mkdir(parents=True)
                (cwd / folder / "chr1.fa").write_text(">chr1\nACGT\n")
        search = ("--pam", pam, "--guide", pam, "--mm", "2", "--bDNA", "1", "--bRNA", "1",
                  "--output", "o")
        build = ("--pam", pam, "--bDNA", "1", "--bRNA", "1")
        for command, extra in (("complete-search", search), ("build-index-only", build)):
            with self.subTest(command=command, where="the working directory"):
                status, said = self.refusal(command, tmp / "w 2", "--genome", "Genomes/gx", *extra)
                self.assertEqual(status, 1)
                self.assertEqual(said.splitlines(), [
                    f"Error: The working directory '{tmp / 'w 2'}' has a space in its path, which "
                    f"{command} cannot take into its commands: run it from a folder whose whole "
                    "path has no space"
                    + (" (or give a --path that has none)" if command == "build-index-only"
                       else "")])
            with self.subTest(command=command, where="the genome folder"):
                status, said = self.refusal(command, tmp / "w", "--genome", "Genomes/g x", *extra)
                self.assertEqual(status, 1)
                self.assertEqual(len(said.splitlines()), 1)
                self.assertIn(f"The --genome folder '{tmp / 'w' / 'Genomes' / 'g x'}' has a space "
                              f"in its path, which {command} cannot take", said)
                self.assertIn(f"(ln -s '{tmp / 'w' / 'Genomes' / 'g x'}' "
                              f"{tmp / 'w' / 'Genomes' / '<a name without spaces>'})", said)
        # ... and the link it advises is taken
        (tmp / "w" / "Genomes" / "gy").symlink_to(tmp / "w" / "Genomes" / "g x")
        os.chdir(tmp / "w")
        self.addCleanup(os.chdir, _HERE)
        self.assertEqual(tbc.crisprme._check_genome(["--genome", "Genomes/gy"], str(tmp / "w")),
                         str(tmp / "w" / "Genomes" / "gy"))


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


class EachMergeTakesItsOwnSortingCriteria(unittest.TestCase):
    """merge_close_targets_cfd.sh chose the criteria of the CFD and CRISPR-Bulge merges by a
    variable nothing set, so every merge took --sorting-criteria and --sorting-criteria-scoring
    did nothing (wset N1). The score merges take it now, and its default is what they took,
    --sorting-criteria's, so a run with neither flag merges as before."""

    def merge(self, pivot):
        """The criteria the merge hands remove_contiguous_samples.py (a stand-in here) for
        a merge of pivot ``pivot``, with scoring criteria 'bulges' and others 'mm'."""
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        (tmp / "remove_contiguous_samples.py").write_text(
            "import sys\nopen('argv.txt', 'w').write(' '.join(sys.argv[1:]))\n")
        subprocess.run(["bash", str(_HERE / "merge_close_targets_cfd.sh"), "in.txt", "out.txt",
                        "3", pivot, "bulges", "mm"], cwd=tmp, check=True, capture_output=True)
        return (tmp / "argv.txt").read_text().split()[-2:]

    def test_the_score_merges_take_the_scoring_criteria(self):
        self.assertEqual(self.merge("score"), ["score", "bulges"])
        self.assertEqual(self.merge("total"), ["total", "mm"])

    def test_with_no_flag_they_merge_as_before(self):
        import test_bulge_cap as tbc
        cm = tbc.crisprme
        self.assertEqual(cm._check_sorting_criteria_scoring([], False),
                         cm._check_sorting_criteria([], False))
        self.assertEqual(cm._check_sorting_criteria([], False), "mm+bulges,mm")
        page = (_REPO / "pages" / "main_page.py").read_text()
        self.assertIn('sorting_criteria_scoring = "mm+bulges,mm"', page)
        self.assertIn('sorting_criteria = "mm+bulges,mm"', page)


class TheSamplePageComputesNoPathItDoesNotRead(unittest.TestCase):
    """The results page's sample tab computed the zip of the sample's personal targets, which
    its download link takes from the name in div-info-sumbysample-targets, and two paths of
    files it no longer reads, and used none of them (wsearch N8)."""

    def test_every_local_it_writes_is_read(self):
        src = (_REPO / "pages" / "results_page.py").read_text()
        page = next(n for n in ast.parse(src).body
                    if isinstance(n, ast.FunctionDef) and n.name == "sample_page")
        written = {n.id for n in ast.walk(page)
                   if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)}
        read = {n.id for n in ast.walk(page)
                if isinstance(n, ast.Name) and not isinstance(n.ctx, ast.Store)}
        # (genome_type: its job-type check of the parameters, kept as it was)
        self.assertEqual(written - read - {"genome_type"}, set())


class TheMergeHasOneDefinitionOfAnIndelToken(unittest.TestCase):
    """remove_contiguous_samples.py read a variant token as an indel's in two places, the
    representative's tie-break and the carrier resolution, each with its own code (RV3-D):
    the tie-break asks the carrier resolution's _holds_indel_token."""

    def test_the_tie_break_asks_the_one_helper(self):
        src = (_HERE / "remove_contiguous_samples.py").read_text()
        functions = {n.name: n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef)}
        tie_break = functions["_rep_tiebreak"]
        calls = {n.func.id for n in ast.walk(tie_break)
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
        self.assertIn("_holds_indel_token", calls)
        self.assertNotIn("rsplit", ast.get_source_segment(src, tie_break))


class PoolSearchIndelsKeepsNoDeadHelper(unittest.TestCase):
    """pool_search_indels.py kept _index_pam, which nothing called (emptyidx F8)."""

    def test_every_function_of_it_is_called(self):
        import re
        src = (_HERE / "pool_search_indels.py").read_text()
        for node in ast.parse(src).body:
            if isinstance(node, ast.FunctionDef):
                with self.subTest(function=node.name):
                    self.assertGreater(len(re.findall(rf"\b{node.name}\b", src)), 1)


class TheReadOnlySmokeCanFail(unittest.TestCase):
    """validate-benchmarks' read-only smoke made the code folder read-only with chmod -R a-w
    in a container that runs as root, which writes through permission bits with
    CAP_DAC_OVERRIDE: a write into the code folder passed, and the step could not fail
    (CI-1). Its container drops the capability, and works in a folder of its own, which
    root owns. (Run by hand on pinellolab/crisprme:v2.8.1 with this tree: a search script
    planted to write into its own folder failed the step, rc 1, "Permission denied"; without
    the --cap-drop it passed, as on GitHub.)"""

    def test_its_container_drops_the_capability(self):
        text = (_REPO / ".github" / "workflows" / "validate-benchmarks.yml").read_text()
        step = text[text.index("- name: read-only-container smoke"):]
        step = step[:step.index("\n      - name:")]
        self.assertIn('docker run --rm --cap-drop DAC_OVERRIDE $BRANCH_MOUNT "$IMG" bash -c', step)
        self.assertNotIn("rowork:/DATA", step)  # (a folder of the runner's user)
        self.assertLess(step.index("mkdir -p /tmp/rowork && cd /tmp/rowork"),
                        step.index("chmod -R a-w /opt/conda/opt/crisprme"))
        checklist = (_REPO / "docs" / "RELEASE_CHECKLIST.md").read_text()
        row = next(ln for ln in checklist.splitlines() if ln.startswith("| T-38 |"))
        self.assertIn("--cap-drop DAC_OVERRIDE", row)


class TheFunctionBasedModulesRunUnderUnittest(unittest.TestCase):
    """Six test modules whose tests are plain functions, run by their own __main__, reported
    "Ran 0 tests" and "OK" under python -m unittest, a false pass (suite P4): CI runs them as
    scripts, but one run by hand or by a loop of `python -m unittest` tested nothing. Each
    has a load_tests hook that collects every test_ function."""

    def test_each_test_function_is_collected(self):
        import importlib
        import test_ci_workflows as tcw
        modules = tcw.function_based_modules()
        self.assertEqual(len(modules), 6, modules)
        for name in modules:
            with self.subTest(module=name):
                module = importlib.import_module(name)
                functions = [k for k, v in vars(module).items()
                             if k.startswith("test_") and callable(v)]
                suite = unittest.defaultTestLoader.loadTestsFromModule(module)
                self.assertEqual(suite.countTestCases(), len(functions))
                self.assertGreater(len(functions), 1)


class AWebTestModuleWithoutDashSkipsItsClasses(unittest.TestCase):
    """Four web test modules raised SkipTest at import where Dash is not installed, which
    `python -m unittest <module>` reports as a traceback and an error, not as skips (skips
    D7). Each class skips now. Dash is made unimportable for the run, as in the light CI."""

    def test_each_skips_and_none_errs(self):
        for name in ("test_web_search_reuse", "test_web_bulge_notes", "test_web_results_links",
                     "test_results_summary_and_status"):
            with self.subTest(module=name):
                run = subprocess.run(
                    [sys.executable, "-c", "import sys; sys.modules['dash'] = None; "
                     f"import unittest; unittest.main(module=None, argv=['x', '{name}'])"],
                    cwd=_HERE, capture_output=True, text=True, timeout=300)
                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertNotIn("Traceback", run.stderr)
                self.assertRegex(run.stderr, r"\nOK \(skipped=\d+\)\n")


class TheCrisprBulgeBenchmarkDoesNothingOnImport(unittest.TestCase):
    """bench_cbulge.py, a developer's benchmark of the CRISPR-Bulge scorer, changed the
    working folder and the module path at import, to a developer's path, or failed to
    import where that path is not (RV4-3). main() takes the repository now."""

    def test_import_changes_neither_the_folder_nor_the_path(self):
        import importlib
        here, path = os.getcwd(), list(sys.path)
        sys.modules.pop("bench_cbulge", None)
        with mock.patch.dict(os.environ, {"CBULGE_REPO": ""}):
            bench = importlib.import_module("bench_cbulge")
        self.addCleanup(sys.modules.pop, "bench_cbulge", None)
        self.assertEqual((os.getcwd(), sys.path), (here, path))
        # run with no repository given, it says what to give, as a usage error
        err = io.StringIO()
        with mock.patch.object(sys, "argv", ["bench_cbulge.py"]), \
                mock.patch.object(sys, "stderr", err), self.assertRaises(SystemExit) as stop:
            bench.main()
        self.assertEqual(stop.exception.code, 2)
        self.assertIn("give the CRISPR-Bulge repository: --repo <folder>, or CBULGE_REPO",
                      err.getvalue())
        self.assertEqual(os.getcwd(), here)


class TheDesignDocSaysWhatShips(unittest.TestCase):
    """DYNAMIC_PAM_INDEX_COMPAT.md called a search that prunes by itself, the web's pruning
    and the PAM picker *planned*, and its CLI section described what the CLI did not do
    (IDX-22): 2.10.0 does them, and the sections say so, and how it differs."""

    def test_the_shipped_sections_are_not_planned(self):
        text = (_REPO / "docs" / "DYNAMIC_PAM_INDEX_COMPAT.md").read_text()
        for heading in ("## 4. CLI", "## 5. Web", "## 8. Web PAM picker"):
            section = text[text.index(heading):].split("\n## ", 1)[0]
            with self.subTest(heading=heading):
                self.assertNotIn("*planned*", section)
                self.assertIn("2.10.0", section)

    def test_no_sentence_says_what_2_10_0_no_longer_does(self):
        # one index for one PAM (an NRG index serves NGG pruned), no all-N index of another
        # width without the flag (a search prunes by itself), every PAM no index serves
        # refused (a search builds one where it can); and what 2.10.0 does not do: the picker
        # reads no tree (its geometry is index_build_pam's), an all-N index as wide as the PAM
        # serves unpruned, with no note, and the PAM list is what offers every PAM
        for doc, gone in (("METHODS.md", "One index therefore answers one PAM"),
                          ("METHODS.md", "does not take an all-`N` index of another"),
                          ("docs/DYNAMIC_PAM_INDEX_COMPAT.md",
                           "A PAM no installed index serves is refused"),
                          ("docs/DYNAMIC_PAM_INDEX_COMPAT.md",
                           "`.pam_build` and the trees (`.bin`)"),
                          ("README.md", "which a search does by itself when no index of its "
                                        "own PAM is installed"),
                          ("README.md", "from an NRR or a PAM-less index: pruned to NGG"),
                          ("docs/DOCKER_QUICKSTART.md", "a search prunes an index to its PAM "
                                                        "by itself when no index of that PAM "
                                                        "is installed"),
                          ("METHODS.md", "(Its PAM picker prefers a PAM whose geometry")):
            with self.subTest(doc=doc, gone=gone):
                self.assertNotIn(gone, " ".join((_REPO / doc).read_text().split()))


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


def _fenced_blocks(text):
    """The fenced code blocks of the Markdown ``text``, each as its text, a line that goes
    on past a backslash joined to the next."""
    blocks, block = [], None
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            if block is not None:
                blocks.append("\n".join(block).replace("\\\n", " "))
            block = None if block is not None else []
        elif block is not None:
            block.append(line)
    return blocks


class TheDocsNameTheFolderABuildWrites(unittest.TestCase):
    """build-index-only names an index for its guide length and side since 2.10.0 (E2:
    20-NRR_3_hg38+hg38_HGDP1kGP), and publish-index publishes the folder it is given.
    PRECOMPUTED_INDEXES' dict-less flow said that its build wrote NRR_3_hg38+hg38_HGDP1kGP,
    so its publish-index of that folder stopped on "Index directory ... does not exist"; the
    release checklist's build said NNN_3_hg38+hg38_HGDP1kGP; and the panel recipes printed
    the older name after their build, whichever CRISPRme they run."""

    DOCS = ("README.md", "docs/PRECOMPUTED_INDEXES.md", "docs/crisprme_data_setup_051826.md",
            "docs/RELEASE_CHECKLIST.md", "docs/SETTINGS_DATA_MANAGER.md")

    @staticmethod
    def built(command):
        """The folder that ``command``, a build-index-only of a PAM file of PAMs/, writes, as
        crisprme.py names it (index_fit.index_name)."""
        import index_fit
        words = command.split()
        flag = dict(zip(words, words[1:]))
        pam_file = _REPO / flag["--pam"]
        motif, length = pam_file.read_text().split()[:2]
        length = int(length)
        rest = Path(flag["--genome"]).name
        if "--vcf" in flag:
            rest += "+" + Path(flag["--vcf"]).name
        n = max(int(flag.get("--bDNA", 0)), int(flag.get("--bRNA", 0))) + 1
        return index_fit.index_name(motif[-length:] if length > 0 else motif[:-length], n,
                                    rest, index_fit.pam_file_spacer(pam_file), length > 0)

    def test_a_documented_build_names_the_folder_it_writes(self):
        checked = 0
        for doc in self.DOCS:
            for block in _fenced_blocks((_REPO / doc).read_text()):
                builds = [line for line in block.splitlines()
                          if "crisprme.py build-index-only" in line]
                said = re.findall(r"->\s*genome_library/([^\s/{]+)", block)
                if len(builds) != 1 or not said:
                    continue
                want = self.built(builds[0])
                for name in said:
                    with self.subTest(doc=doc, name=name):
                        self.assertIn(name, (want, want + "_INDELS"))
                    checked += 1
        self.assertGreaterEqual(checked, 3)  # PRECOMPUTED_INDEXES' two, the checklist's

    def test_a_documented_publish_names_a_folder_a_build_writes(self):
        import index_fit
        checked = 0
        for doc in self.DOCS:
            for line in (_REPO / doc).read_text().splitlines():
                for name in re.findall(r"publish-index --index genome_library/([^\s`<]+)", line):
                    with self.subTest(doc=doc, name=name):
                        # (a placeholder, genome_library/<name>, is no name)
                        self.assertIsNotNone(index_fit.parse_index_name(name)[3], line)
                    checked += 1
        self.assertGreaterEqual(checked, 5)

    def test_the_panel_recipes_print_the_folder_their_build_wrote(self):
        # as a CRISPRme before 2.10.0 names it (the as-run SIFs and overlay), and as 2.10.0
        # does, for the PAM the recipe takes by default or another one (PAM=...)
        for recipe, ds in (("hgdp1kgp_build.sh", "hg38_HGDP1kGP"),
                           ("hprc_build.sh", "hg38_HPRC"),
                           ("onekg2021_build.sh", "hg38_1000G2021")):
            script = (_REPO / "seq_script" / "merge_panels" / recipe).read_text()
            after_build = script[script.index('echo "BUILD_EXIT='):]
            for name in (f"NRG_3_hg38+{ds}", f"20-NRG_3_hg38+{ds}", f"20-NRR_3_hg38+{ds}"):
                with self.subTest(recipe=recipe, name=name), \
                        tempfile.TemporaryDirectory() as build:
                    for folder in (name, name + "_INDELS"):
                        os.makedirs(os.path.join(build, "genome_library", folder))
                    run = subprocess.run(["bash", "-c", after_build], capture_output=True,
                                         text=True, env={**os.environ, "B": build})
                    self.assertEqual(run.stdout.splitlines()[1:],
                                     [f"index: {build}/genome_library/{name} (+ _INDELS)"])
                    self.assertEqual(run.stderr, "")

    def test_the_mega_indels_recipe_says_which_name_to_set(self):
        # its guard stops on an INDELS index of another name than IDXNAME says, and named
        # the older form alone, <MOTIF>_<N>, which a 2.10.0 checkout does not write
        script = (_REPO / "seq_script" / "merge_panels" / "mega_build_indels.sh").read_text()
        guard = script[script.index("if [ ! -d \"$W/genome_library/$IDXDIR\" ]"):]
        guard = guard[:guard.index("continue")]
        self.assertIn("<guide>-<MOTIF>_<N>", guard)
        self.assertIn("e.g. 20-NRG_3", guard)


class TheDocsSayWhatARunWrites(unittest.TestCase):
    """Texts of the backlog (plan 4.8) that said what a run does not write: the web user
    guide called every ALT row "created by a variant allele" (RV1-WEB2-3), METHODS put
    targets-scoring's CFD-only banner "in its log" (D3R3-3), the README and the help of
    complete-search said that "the CRISPR-Bulge columns read -1" in a CFD-only run, whose
    CRISPR-Bulge risk score reads 0 (D3R3-2), and no text said that a line of the SNP+indel
    companion names 50 carriers at most (LP RV4-cooc-carrier-cap)."""

    @staticmethod
    def text(path):
        return " ".join((_REPO / path).read_text().split())

    def says(self, text, phrase, where):
        self.assertTrue(phrase in text, f"{where} does not say {phrase!r}")

    def says_not(self, text, phrase, where):
        self.assertFalse(phrase in text, f"{where} still says {phrase!r}")

    def test_an_alt_row_is_a_variant_alignment(self):
        # as the results page's caption says: an ALT row is a variant alignment, of a site
        # a variant creates (Not_found_in_REF y) or of a reference site it changes
        self.says(self.text("pages/results_page.py"),
                  "a variant either creates the site or changes a reference", "results_page")
        guide = self.text("docs/crisprme_web_interface_user_guide.md")
        self.says_not(guide, "ALT = created by a variant allele", "the user guide")
        self.says(guide, "a variant either creates the site (`Not_found_in_REF` `y`) or "
                         "changes a reference site", "the user guide")

    def test_targets_scoring_says_cfd_only_on_stdout(self):
        # scorer_runner prints the banner on stdout: a search's log_verbose.txt holds it,
        # and targets-scoring writes no log
        self.says(self.text("PostProcess/scorer_runner.py"), "print(banner, flush=True)",
                  "scorer_runner")
        methods = self.text("METHODS.md")
        self.says_not(methods, "a banner saying so in its log", "METHODS")
        self.says(methods, "a search in its `log_verbose.txt` and `targets-scoring` on stdout",
                  "METHODS")

    def test_a_cfd_only_run_reads_0_for_the_crispr_bulge_risk_score(self):
        # add_risk_score.py writes ALT less REF: -1 less -1
        tree = ast.parse((_REPO / "crisprme.py").read_text())
        helper = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                      and node.name == "print_help_complete_search")
        help_text = " ".join(" ".join(node.value for node in ast.walk(helper)
                                      if isinstance(node, ast.Constant)
                                      and isinstance(node.value, str)).split())
        for where, text, said in (
                ("README.md", self.text("README.md"),
                 "(the CRISPR-Bulge scores read `-1`, and its risk score `0`)"),
                ("complete-search's help", help_text,
                 "(the CRISPR-Bulge scores read -1, and its risk score 0)")):
            with self.subTest(where=where):
                self.says_not(text, "CRISPR-Bulge columns read", where)
                self.says(text, said, where)

    def test_the_cooc_companion_says_it_names_50_carriers(self):
        # n_cis counts every carrier, cis_samples names the first 50 (analisi_indels_NNN.py)
        self.says(self.text("PostProcess/analisi_indels_NNN.py"),
                  'str(len(samples)) + "\\t" + ",".join(sorted(samples)[:50])',
                  "analisi_indels_NNN.py")
        for doc in ("README.md", "METHODS.md"):
            with self.subTest(doc=doc):
                self.says(self.text(doc), "the first 50 by name", doc)


if __name__ == "__main__":
    unittest.main()
