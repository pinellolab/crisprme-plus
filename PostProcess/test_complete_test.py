"""Fast, hermetic unit tests for ``complete_test.run_crisprme_test`` dispatch.

REGRESSION GUARD (the exact situation this was written for): the v2.6.0 complete-test
rewrite made ``run_crisprme_test`` *always* require a prebuilt index and run a single
example-guide smoke, which silently dropped the ``validate-benchmarks`` CI's download +
per-benchmark on-demand-build flow -> ``validate-test`` had no predictions to compare
against the brute-force ground truth and the gate went red (unnoticed because the heavy
benchmark job only runs on PR/weekly-cron, not push-to-main).

These tests pin the two-mode contract WITHOUT any download / index build / subprocess, so
they run in the fast unit suite (unit-tests.yml, which DOES run on push-to-main):

  * CI / brute-force gate  (env CRISPRME_ALLOW_ONDEMAND_BUILD set): download the test data
    and run EVERY registered benchmark via complete-search with ON-DEMAND build (no
    --index-path), one per-benchmark output dir (crisprme-test-out_<name>).
  * default USER path (flag unset): require a prebuilt index, run a single smoke against it
    via --index-path, and NEVER download/build.

The smoke honours --chrom. It used not to receive it, so `complete-test --chrom chr22`
searched every chromosome of the index's genome: on a two-chromosome genome it returned
10 rows on chr22 and 10 on chrA, and with NRG_3_hg38+hg38_HGDP1kGP it was a genome-wide
variant search. Now a chromosome other than "all" runs the same complete-search from a
per-run view of the working directory in which it is the only chromosome, removed on
every exit, and one the genome or the index lacks is refused in one line before anything
is searched. A chromosome holding no off-target of the example guide passes with an
empty result, as complete-search ends it. A SIGTERM or SIGHUP is passed on to the
search, which is waited for before the view goes, and one ignored at entry (nohup's
SIGHUP) stays ignored. ChromView checks the view in-process on a small working
directory, with complete-search stubbed; RealSignals sends real signals to
complete_test.py run as a process, with a stub crisprme.py on PATH that starts a few
processes; StopSearch checks the kill of a search that outlasts its wait; the CLI test
checks that crisprme.py passes the refusal on without a traceback.
"""

import contextlib
import io
import os
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import time
import types
import unittest
from unittest import mock

import complete_test as ct

HERE = os.path.dirname(os.path.abspath(__file__))


_REGISTRY = {
    "thresholds": {"mm": 4, "bDNA": 1, "bRNA": 1},
    "benchmarks": [
        {"name": "cas9_a", "nuclease": "SpCas9", "pam_name": "pA.txt",
         "pam_content": "NNNNNNNNNNNNNNNNNNNNNGG 3", "guide_file": "gA.txt",
         "guide_crisprme": "A" * 20 + "NGG"},
        {"name": "cas12a_b", "nuclease": "Cas12a", "pam_name": "pB.txt",
         "pam_content": "TTTVNNNNNNNNNNNNNNNNNNNNNNNN 3", "guide_file": "gB.txt",
         "guide_crisprme": "TTTV" + "C" * 23},
    ],
}


class DispatchContract(unittest.TestCase):
    def _run(self, allow_ondemand, chrom="chr22", skip_heavy=False, registry=_REGISTRY):
        """Run run_crisprme_test with every side-effecting dependency mocked; return the
        list of complete-search command strings + a flag of whether require_prebuilt_index
        was consulted. The keyword arguments of each call are kept in self.kwargs."""
        calls = []
        self.kwargs = []
        pam_file = mock.patch.object(ct, "_smoke_pam_file", side_effect=lambda _line, name: name)
        pam_file.start()
        self.addCleanup(pam_file.stop)
        env = dict(os.environ)
        env.pop("CRISPRME_ALLOW_ONDEMAND_BUILD", None)
        env.pop("CRISPRME_SKIP_HEAVY", None)
        if allow_ondemand:
            env["CRISPRME_ALLOW_ONDEMAND_BUILD"] = "1"
        if skip_heavy:
            env["CRISPRME_SKIP_HEAVY"] = "1"
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(ct, "check_crisprme_directory_tree", lambda *_a, **_k: None), \
             mock.patch.object(ct, "check_output", lambda *_a, **_k: None), \
             mock.patch.object(ct, "subprocess") as msub, \
             mock.patch.object(ct, "require_prebuilt_index") as mreq, \
             mock.patch.object(ct, "download_genome_data", return_value="Genomes/hg38_chr22") as mdl, \
             mock.patch.object(ct, "download_vcf_data", lambda *_a, **_k: None), \
             mock.patch.object(ct, "write_vcf_config", return_value="vcf.config.txt"), \
             mock.patch.object(ct, "download_samples_ids_data", lambda *_a, **_k: None), \
             mock.patch.object(ct, "write_samplesids_config", return_value="samples.config.txt"), \
             mock.patch.object(ct, "download_annotation_data", return_value=("gencode.bed", "encode.bed")), \
             mock.patch.object(ct, "load_benchmarks", return_value=registry), \
             mock.patch.object(ct, "write_pamfile", side_effect=lambda name, *_a, **_k: name), \
             mock.patch.object(ct, "write_guidefile", side_effect=lambda name, *_a, **_k: name), \
             mock.patch.object(ct, "write_sg1617_guidefile", return_value="sg1617.guide.txt"), \
             mock.patch.object(ct, "_parse_index_name", return_value=("NRG", "hg38_HGDP1kGP", True)), \
             mock.patch.object(ct, "glob", return_value=["x.integrated_results.tsv"]), \
             mock.patch("os.path.isdir", return_value=True), \
             mock.patch("os.path.isfile", return_value=True), \
             contextlib.redirect_stderr(io.StringIO()):
            mreq.return_value = "NRG_3_hg38+hg38_HGDP1kGP"
            msub.call.side_effect = lambda cmd, **k: calls.append(cmd) or self.kwargs.append(k) or 0
            ct.run_crisprme_test(chrom, "1000G", 4, False)
        return calls, mreq.called, mdl.called

    def test_ci_mode_downloads_and_builds_per_benchmark_on_demand(self):
        """CRISPRME_ALLOW_ONDEMAND_BUILD set -> download + one on-demand-build complete-search
        per registered benchmark; NEVER require_prebuilt_index. This is the exact flow the
        brute-force benchmark gate needs; its loss is what this guard catches."""
        calls, req_called, dl_called = self._run(allow_ondemand=True)
        self.assertTrue(dl_called, "CI mode must download the test genome")
        self.assertFalse(req_called, "CI mode must NOT gate on a prebuilt index")
        self.assertEqual(len(calls), len(_REGISTRY["benchmarks"]),
                         "CI mode must run one complete-search per registered benchmark")
        for cmd in calls:
            self.assertIn("complete-search", cmd)
            self.assertNotIn("--index-path", cmd,
                             "CI mode must BUILD on demand (no --index-path), not reuse a prebuilt index")
        # each benchmark gets its own crisprme-test-out_<name> dir (validate-test reads these)
        for bench in _REGISTRY["benchmarks"]:
            self.assertTrue(any(f"{ct.COMPLETETESTRESDIR}_{bench['name']}" in c for c in calls),
                            f"missing per-benchmark output dir for {bench['name']}")

    def test_a_heavy_benchmark_runs_unless_the_opt_out_is_set(self):
        """A "heavy" (bulge-2) benchmark is part of the self-test: the whole set fits a
        16 GiB, 4-CPU runner. CRISPRME_SKIP_HEAVY=1 leaves it out on a very small machine,
        and only it."""
        registry = dict(_REGISTRY, benchmarks=[
            _REGISTRY["benchmarks"][0], dict(_REGISTRY["benchmarks"][1], heavy=True)])
        heavy_dir = f"{ct.COMPLETETESTRESDIR}_{registry['benchmarks'][1]['name']}"
        calls, _, _ = self._run(allow_ondemand=True, registry=registry)
        self.assertEqual(len(calls), 2)
        self.assertTrue(any(heavy_dir in c for c in calls))
        calls, _, _ = self._run(allow_ondemand=True, registry=registry, skip_heavy=True)
        self.assertEqual(len(calls), 1)
        self.assertFalse(any(heavy_dir in c for c in calls))

    def test_user_mode_requires_prebuilt_index_and_never_builds(self):
        """Flag unset -> require a prebuilt index + a single smoke via --index-path; NEVER
        download/build (Luca's rule: CRISPRme does not build an index for the user)."""
        calls, req_called, dl_called = self._run(allow_ondemand=False, chrom="all")
        self.assertTrue(req_called, "user mode must require a prebuilt index")
        self.assertFalse(dl_called, "user mode must NOT download/build test data")
        self.assertEqual(len(calls), 1, "user mode runs exactly one example smoke")
        self.assertIn("--index-path genome_library", calls[0],
                      "user-mode smoke must reuse the prebuilt index (--index-path)")

    def test_user_mode_whole_genome_command_is_unchanged(self):
        """Without --chrom (or with all) the smoke runs the command it always ran, from
        the working directory."""
        calls, _, _ = self._run(allow_ondemand=False, chrom="all")
        genome = os.path.join(os.getcwd(), "Genomes", "hg38_HGDP1kGP")
        self.assertEqual(calls, [
            f"crisprme.py complete-search --genome {genome} --mm 2 --bDNA 1 --bRNA 1 "
            "--max-total-edits 4 --merge 3 --pam 20bp-NRG-smoke.txt --guide "
            "sg1617.guide.txt --vcf list_vcf.txt --samplesID list_samplesID.txt "
            "--index-path genome_library --output crisprme-test-out --thread 4  --ci-cd-test"
        ])
        self.assertTrue(self.kwargs[0].get("shell"))
        self.assertIsNone(self.kwargs[0].get("cwd"))


class AnEarlierRunInTheWayIsNamedWhole(unittest.TestCase):
    """complete-test refuses to run over the results of an earlier run. It named one folder
    (of the benchmark mode's six, crisprme-test-out_cas12a_hbg: a user who removed it was
    stopped five more times), refused the prebuilt-index smoke on those six too, which it
    does not write, and exited 0, which a script or a CI step reads as a pass though
    nothing ran. It names every folder of this run's mode there, in one line (the benchmark
    mode's crisprme-test-out_<name> of any benchmark, the smoke's crisprme-test-out), and
    exits 2 before it downloads or searches anything; crisprme.py passes the 2 on."""

    SCRIPT = os.path.join(HERE, "complete_test.py")

    def setUp(self):
        self.wd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.wd)
        cwd = os.getcwd()
        os.chdir(self.wd)
        self.addCleanup(os.chdir, cwd)
        self.results = os.path.join(os.getcwd(), "Results")
        self.benchmarks = sorted(f"{ct.COMPLETETESTRESDIR}_{bench['name']}"
                                 for bench in ct.load_benchmarks()["benchmarks"])

    @staticmethod
    def env(ondemand, skip_heavy=None):
        env = {k: v for k, v in os.environ.items()
               if k not in ("CRISPRME_ALLOW_ONDEMAND_BUILD", "CRISPRME_SKIP_HEAVY")}
        if ondemand:
            env["CRISPRME_ALLOW_ONDEMAND_BUILD"] = "1"
        if skip_heavy:
            env["CRISPRME_SKIP_HEAVY"] = skip_heavy
        return env

    def run_test(self, folders, ondemand, skip_heavy=None):
        """run_crisprme_test with ``folders`` in Results/, the benchmark mode for
        ``ondemand``: its exit status, its stderr, and the mode it ran (none if refused)."""
        for name in folders:
            os.makedirs(os.path.join(self.results, name), exist_ok=True)
        ran, err = [], io.StringIO()
        with mock.patch.dict(os.environ, self.env(ondemand, skip_heavy), clear=True), \
             mock.patch.object(ct, "_run_benchmarks_ondemand", lambda *a: ran.append("bench")), \
             mock.patch.object(ct, "_run_prebuilt_index_smoke", lambda *a: ran.append("smoke")), \
             contextlib.redirect_stderr(err):
            try:
                ct.run_crisprme_test("chr22", "1000G", 4, False)
                code = 0
            except SystemExit as stop:
                code = stop.code
        return code, err.getvalue(), ran

    def said(self, names):
        return ("Complete-test already run once. Please delete the complete-test results "
                f"folder{'s' if len(names) > 1 else ''} before running it again: "
                + " ".join(os.path.join(self.results, n) for n in names) + "\n")

    def test_the_benchmark_modes_folders_are_named_in_one_line(self):
        self.assertEqual(len(self.benchmarks), 6)
        # a heavy one's too under CRISPRME_SKIP_HEAVY=1 (a validate-test without it reads
        # it), and, as before, one of a benchmark an earlier release ran
        named = sorted(self.benchmarks + [f"{ct.COMPLETETESTRESDIR}_cas9"])
        for skip_heavy in (None, "1"):
            with self.subTest(CRISPRME_SKIP_HEAVY=skip_heavy):
                code, err, ran = self.run_test(named + [ct.COMPLETETESTRESDIR], True, skip_heavy)
                self.assertEqual((code, ran), (2, []))
                self.assertEqual(err, self.said(named))

    def test_only_those_that_are_there(self):
        code, err, ran = self.run_test(self.benchmarks[3:5], True)
        self.assertEqual((code, err, ran), (2, self.said(self.benchmarks[3:5]), []))

    def test_the_smoke_is_refused_on_its_own_folder_only(self):
        code, err, ran = self.run_test(self.benchmarks, False)
        self.assertEqual((code, err, ran), (0, "", ["smoke"]))
        code, err, ran = self.run_test([ct.COMPLETETESTRESDIR], False)
        self.assertEqual((code, err, ran), (2, self.said([ct.COMPLETETESTRESDIR]), []))

    def test_the_benchmark_mode_is_not_refused_on_the_smokes(self):
        code, err, ran = self.run_test([ct.COMPLETETESTRESDIR], True)
        self.assertEqual((code, err, ran), (0, "", ["bench"]))

    def test_the_script_exits_2_having_downloaded_nothing(self):
        os.makedirs(os.path.join(self.results, self.benchmarks[0]))
        res = subprocess.run([sys.executable, self.SCRIPT, "chr22", "1000G", "4", "False"],
                             env=self.env(True), capture_output=True, text=True, timeout=120)
        self.assertEqual((res.returncode, res.stdout, res.stderr),
                         (2, "", self.said(self.benchmarks[:1])))
        self.assertEqual(os.listdir("Genomes"), [])


class NothingToSearchIsNoPass(unittest.TestCase):
    """The prebuilt-index smoke (complete-test without CRISPRME_ALLOW_ONDEMAND_BUILD) in a
    folder with no prebuilt index, or with one whose genome is not there, says how to get
    it and exited 0, which a script or a CI step reads as a pass though nothing ran, as the
    re-run guard beside it says (AnEarlierRunInTheWayIsNamedWhole). It exits 2 now, as
    that guard does, having said the same and searched and downloaded nothing; crisprme.py
    passes the 2 on (TheCliPassesTheRefusalOn)."""

    SCRIPT = os.path.join(HERE, "complete_test.py")

    def setUp(self):
        self.wd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.wd)

    def run_script(self):
        env = {k: v for k, v in os.environ.items() if k != "CRISPRME_ALLOW_ONDEMAND_BUILD"}
        return subprocess.run([sys.executable, self.SCRIPT, "chr22", "1000G", "4", "False"],
                              cwd=self.wd, env=env, capture_output=True, text=True,
                              timeout=120)

    def assertNothingRan(self, res, says):
        self.assertEqual(res.returncode, 2, res.stderr)
        self.assertIn(says, res.stderr)
        self.assertNotIn("Traceback", res.stderr)
        for folder in ("Genomes", "Results", "VCFs"):
            self.assertEqual(os.listdir(os.path.join(self.wd, folder)), [], folder)

    def test_no_prebuilt_index(self):
        self.assertNothingRan(
            self.run_script(),
            "No prebuilt index found for the complete-test (looked in "
            f"{os.path.join(os.path.realpath(self.wd), 'genome_library')}).\n"
            "CRISPRme does not build one automatically.\n\nDownload an index first, e.g.:")

    def test_an_index_whose_genome_is_not_there(self):
        index = os.path.join(self.wd, "genome_library", "NRG_3_hg38")
        os.makedirs(index)
        with open(os.path.join(index, "X_chr22_1.bin"), "wb") as fh:
            fh.write(b"t")
        self.assertNothingRan(
            self.run_script(),
            "The prebuilt index 'NRG_3_hg38' is present, but its reference genome "
            "(Genomes/hg38) is not.")


class WhatTheDataRepoDoesNotHostIsNoFailure(unittest.TestCase):
    """The documented self-test (CRISPRME_ALLOW_ONDEMAND_BUILD=1 complete-test --chrom chr22
    --vcf_dataset 1000G) always said 'HuggingFace fetch failed for
    ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz (404 Client
    Error ... Entry Not Found ...); falling back to original source', then downloaded it
    from EBI and passed: the data repo hosts no raw 1000 Genomes or HGDP VCF, nor
    samplesIDs.HGDP.txt, nor the unpacked FASTA of a chromosome but chr22. A file the Hub
    answers it does not hold is said in one line on stdout; a failure to fetch one (no
    connection, a checksum that does not match) is still said on stderr, as before.

    The Hub's errors are huggingface_hub's own where it is installed; the light CI env has
    none, and stand-ins of their hierarchy take its place for the test."""

    VCF = "ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

    def setUp(self):
        try:
            from huggingface_hub.utils import EntryNotFoundError, LocalEntryNotFoundError
        except ImportError:
            class EntryNotFoundError(Exception):
                pass

            class LocalEntryNotFoundError(FileNotFoundError, EntryNotFoundError):
                pass

            hub = types.ModuleType("huggingface_hub")
            hub.utils = types.ModuleType("huggingface_hub.utils")
            hub.utils.EntryNotFoundError = EntryNotFoundError
            hub.utils.LocalEntryNotFoundError = LocalEntryNotFoundError
            patcher = mock.patch.dict(sys.modules, {"huggingface_hub": hub,
                                                    "huggingface_hub.utils": hub.utils})
            patcher.start()
            self.addCleanup(patcher.stop)
        # what hf_hub_download raises for a file the repo does not hold (named
        # RemoteEntryNotFoundError from huggingface_hub 1.0 on), and with no connection
        # and the file not cached
        self.missing = type("RemoteEntryNotFoundError", (EntryNotFoundError,), {})(
            "404 Client Error. Entry Not Found for url: https://huggingface.co/datasets/...")
        self.offline = LocalEntryNotFoundError(
            "An error happened while trying to locate the file on the Hub and we cannot find "
            "the requested files in the local cache.")
        self.wd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.wd)
        cwd = os.getcwd()
        os.chdir(self.wd)
        self.addCleanup(os.chdir, cwd)

    def fetch(self, call, hf_error, md5s=()):
        """``call`` with hf_fetch raising ``hf_error`` (None: fetching the file), the
        other downloads faked, and compute_md5 answering ``md5s``, then each file's MD5
        from the tables: its stdout, its stderr and the URLs it downloaded from."""
        tables = dict(ct.MD51000G, **ct.MD5HGDP, **ct.MD5SAMPLES)
        md5s, fetched = list(md5s), []

        def hf_fetch(_name, dest):
            if hf_error is not None:
                raise hf_error
            open(dest, "w").close()

        def download(dest, fname=None, ftp_conn=False, ftp_server=None, ftp_path=None,
                     http_url=None):
            fetched.append(http_url or f"ftp://{ftp_server}{ftp_path}")
            path = os.path.join(dest, os.path.basename(fetched[-1]))
            open(path, "w").close()
            return path

        def gunzip(_gz, out):
            open(out, "w").close()
            return out

        def compute_md5(path):
            return md5s.pop(0) if md5s else tables.get(os.path.basename(path))

        out, err = io.StringIO(), io.StringIO()
        with mock.patch.object(ct, "hf_fetch", hf_fetch), \
             mock.patch.object(ct, "download", download), \
             mock.patch.object(ct, "gunzip", gunzip), \
             mock.patch.object(ct, "compute_md5", compute_md5), \
             contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            call()
        return out.getvalue(), err.getvalue(), fetched

    def test_the_raw_vcfs_of_the_self_test(self):
        os.makedirs("VCFs")
        for ds, fname, source in (
                ("1000G", self.VCF, "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
                 "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"),
                ("HGDP", "hgdp_wgs.20190516.full.chr22.vcf.gz",
                 "ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/")):
            with self.subTest(dataset=ds):
                out, err, fetched = self.fetch(
                    lambda: ct.download_vcf_data("chr22", "VCFs", ds), self.missing)
                self.assertEqual(out, f"{fname} is not on HuggingFace ({ct.HF_DATA_REPO}); "
                                      "downloading it from its original source\n")
                self.assertEqual(err, "Downloading VCF data for chromsome(s) chr22\n")
                self.assertEqual(fetched, [source + fname])

    def test_the_hgdp_sample_list(self):
        out, err, fetched = self.fetch(lambda: ct.download_samples_ids_data("HGDP"),
                                       self.missing)
        self.assertEqual(out, f"samplesIDs.HGDP.txt is not on HuggingFace ({ct.HF_DATA_REPO}); "
                              "downloading it from its original source\n")
        self.assertEqual(err, "Downloading sample ids for dataset(s) HGDP\n")
        self.assertEqual(fetched, [f"{ct.TESTDATAURL.rstrip('/')}/samplesIDs/samplesIDs.HGDP.txt"])
        self.assertTrue(os.path.isfile(os.path.join("samplesIDs", "hg38_HGDP.samplesID.txt")))

    def test_the_fasta_of_another_chromosome(self):
        os.makedirs("Genomes")
        out, err, fetched = self.fetch(lambda: ct.download_genome_data("chr21", "Genomes"),
                                       self.missing)
        self.assertEqual(out, f"chr21.fa is not on HuggingFace ({ct.HF_DATA_REPO}); "
                              "downloading it from UCSC\n")
        self.assertEqual(err, "Downloading fasta file for chromosome(s) chr21\n")
        self.assertEqual(fetched, [f"{ct.HG38URL}/chromosomes/chr21.fa.gz"])

    def test_a_failure_to_fetch_is_still_said(self):
        os.makedirs("VCFs")
        for hf_error, md5s, why in ((self.offline, (), self.offline),
                                    (None, ("0" * 32,), "MD5 mismatch after HuggingFace fetch")):
            with self.subTest(why=str(why)[:40]):
                out, err, fetched = self.fetch(
                    lambda: ct.download_vcf_data("chr22", "VCFs", "1000G"), hf_error, md5s)
                self.assertEqual(out, "")
                self.assertEqual(err, "Downloading VCF data for chromsome(s) chr22\n"
                                      f"HuggingFace fetch failed for {self.VCF} ({why}); "
                                      "falling back to original source\n")
                self.assertEqual(len(fetched), 1)
                os.remove(os.path.join("VCFs", "hg38_1000G", self.VCF))
        os.makedirs("Genomes")
        out, err, _ = self.fetch(lambda: ct.download_genome_data("chr22", "Genomes"),
                                 self.offline)
        self.assertEqual((out, err), ("", "Downloading fasta file for chromosome(s) chr22\n"
                                          f"HuggingFace fetch failed for chr22.fa "
                                          f"({self.offline}); falling back to UCSC\n"))

    def test_what_tells_one_from_the_other(self):
        self.assertTrue(ct.hf_missing(self.missing))
        for error in (self.offline, ValueError("MD5 mismatch after HuggingFace fetch"),
                      ImportError("No module named 'huggingface_hub'"), OSError("timed out")):
            self.assertFalse(ct.hf_missing(error), error)
        with mock.patch.dict(sys.modules, {"huggingface_hub": None,
                                           "huggingface_hub.utils": None}):
            self.assertFalse(ct.hf_missing(self.missing))  # (no huggingface_hub to ask)


# --------------------------------------------------------------------------- #
# A small working directory, as `download --what index` or build-index-only leave one:
# genome toy (chr2, chr22, a chr22_ random contig, and chrY no index has), a reference
# and a variant index with its INDELS companion and dictionaries, plus what a view must
# leave out: another genome's index, a locked one, a hidden one, another dataset's.
# --------------------------------------------------------------------------- #
_FILES = {
    "Genomes/toy/chr2.fa": ">chr2\nACGT\n",
    "Genomes/toy/chr22.fa": ">chr22\nACGT\n",
    "Genomes/toy/chr22_KI270731v1_random.fa": ">chr22_KI270731v1_random\nACGT\n",
    "Genomes/toy/chrY.fa": ">chrY\nACGT\n",
    "Genomes/toy+toy_panel/chr22.enriched.fa": ">chr22\nACGT\n",
    "genome_library/NGG_2_toy/.pam_build": "NGG 3\n",
    "genome_library/NGG_2_toy/NGG_chr2_1.bin": "t",
    "genome_library/NGG_2_toy/NGG_chr22_1.bin": "t",
    "genome_library/NGG_2_toy/NGG_chr22_2.bin": "t",
    "genome_library/NGG_2_toy/NGG_chr22_KI270731v1_random_1.bin": "t",
    "genome_library/NGG_2_toy+toy_panel/.pam_build": "NGG 3\n",
    "genome_library/NGG_2_toy+toy_panel/.display_label": "toy panel\n",
    "genome_library/NGG_2_toy+toy_panel/NGG_chr2_1.bin": "t",
    "genome_library/NGG_2_toy+toy_panel/NGG_chr22_1.bin": "t",
    "genome_library/NGG_2_toy+toy_panel/NGG_chr22_KI270731v1_random_1.bin": "t",
    "genome_library/NGG_2_toy+toy_panel/NGG_chrM_1.bin": "t",
    "genome_library/NGG_2_toy+toy_panel_INDELS/NGG_2_fakechr2/NGG_fakechr2_1.bin": "t",
    "genome_library/NGG_2_toy+toy_panel_INDELS/NGG_2_fakechr22/NGG_fakechr22_1.bin": "t",
    "genome_library/NRG_3_toy/NRG_chr22_1.bin": "t",
    "genome_library/NGG_4_toy/NGG_chr22_1.bin": "t",
    "genome_library/NGG_4_toy.lock/owner": "1\n",
    "genome_library/.NGG_2_toy/NGG_chr22_1.bin": "t",
    "genome_library/NGG_3_other/NGG_chr22_1.bin": "t",
    "genome_library/NGG_2_toy+zz_panel/NGG_chr22_1.bin": "t",
    "genome_library/NGG_2_toy+zz_panel_INDELS/NGG_2_fakechr22/NGG_fakechr22_1.bin": "t",
    "Dictionaries/dictionaries_toy_panel/my_dict_chr22.json.gz": "d",
    "Dictionaries/registry_toy_panel/reg_chr22.idx": "{}",
    "Dictionaries/genotypes_toy_panel/gt_chr22.bin": "g",
    "Dictionaries/log_indels_toy_panel/logchr2.txt.gz": "l",
    "Dictionaries/log_indels_toy_panel/logchr22.txt.gz": "l",
    "Dictionaries/log_indels_zz_panel/logchr22.txt.gz": "l",
    "VCFs/toy_panel/chr2.vcf.gz": "v",
    "VCFs/toy_panel/chr22.vcf.gz": "v",
    "samplesIDs/toy.samplesID.txt": "#SAMPLE_ID\n",
    "Annotations/genes.bed": "",
    "PAMs/.keep": "",
    "PAMs/20bp-NGG-SpCas9.txt": "NNNNNNNNNNNNNNNNNNNNNGG 3\n",
    "Results/.keep": "",
    "list_vcf.txt": "toy_panel\n",
    "list_samplesID.txt": "toy.samplesID.txt\n",
    "notes.txt": "a file of the user's own\n",
}

# the view a --chrom chr22 smoke runs from, every path in it: None for a folder, else
# "link", a symlink to the same path in the working directory
_VIEW_CHR22 = {
    "Annotations": "link", "PAMs": "link", "Results": "link", "samplesIDs": "link",
    "list_vcf.txt": "link", "list_samplesID.txt": "link", "notes.txt": "link",
    "sg1617_test_guide.txt": "link",
    "Genomes": None, "Genomes/toy": None, "Genomes/toy/chr22.fa": "link",
    "VCFs": None,
    "genome_library": None,
    "genome_library/NGG_2_toy": None,
    "genome_library/NGG_2_toy/NGG_chr22_1.bin": "link",
    "genome_library/NGG_2_toy/NGG_chr22_2.bin": "link",
    "genome_library/NGG_2_toy+toy_panel": None,
    "genome_library/NGG_2_toy+toy_panel/NGG_chr22_1.bin": "link",
    "genome_library/NGG_2_toy+toy_panel_INDELS": None,
    "genome_library/NGG_2_toy+toy_panel_INDELS/NGG_2_fakechr22": "link",
    "genome_library/NRG_3_toy": None,
    "genome_library/NRG_3_toy/NRG_chr22_1.bin": "link",
    "Dictionaries": None,
    "Dictionaries/dictionaries_toy_panel": "link",
    "Dictionaries/registry_toy_panel": "link",
    "Dictionaries/genotypes_toy_panel": "link",
    "Dictionaries/log_indels_zz_panel": "link",
    "Dictionaries/log_indels_toy_panel": None,
    "Dictionaries/log_indels_toy_panel/logchr22.txt.gz": "link",
}


def _make_workdir(root):
    for rel, text in _FILES.items():
        path = os.path.join(root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)


def _tree(root):
    """{relative path: None (folder) | ("link", target) | ("file", content)}, not
    following symlinks."""
    out = {}

    def walk(d):
        for e in os.scandir(d):
            rel = os.path.relpath(e.path, root)
            if e.is_symlink():
                out[rel] = ("link", os.readlink(e.path))
            elif e.is_dir():
                out[rel] = None
                walk(e.path)
            else:
                with open(e.path) as fh:
                    out[rel] = ("file", fh.read())
    walk(root)
    return out


def _views(wd):
    return [n for n in os.listdir(wd) if n.startswith(ct.CHROMVIEWPREFIX)]


class _SmokeDir(unittest.TestCase):
    """A working directory with the shipped PAM files, and indexes by name."""

    def setUp(self):
        self._cwd = os.getcwd()
        self.wd = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.wd, ignore_errors=True)
        os.makedirs(os.path.join(self.wd, "genome_library"))
        shutil.copytree(os.path.join(os.path.dirname(HERE), "PAMs"), os.path.join(self.wd, "PAMs"))
        os.chdir(self.wd)
        self.addCleanup(os.chdir, self._cwd)
        self.addCleanup(ct._smoke_pam_dir.clear)

    def index(self, name, record=None, guide=None):
        folder = os.path.join(self.wd, "genome_library", name)
        os.makedirs(folder)
        with open(os.path.join(folder, "X_chr1_1.bin"), "wb") as fh:
            fh.write(struct.pack("=ii", 7, guide) if guide else b"t")
        if record:
            with open(os.path.join(folder, ".pam_build"), "w") as fh:
                fh.write(record)


class WhichIndexTheSmokeTakes(_SmokeDir):
    """It took the first name in sorted order: NNNNNNNN_3 and NNGRRT_2 sort before NRG_3,
    so that installing the published K=8 index beside the default broke the smoke."""

    def taken(self, *names):
        for name in names:
            self.index(name)
        return ct.require_prebuilt_index()

    def test_a_dedicated_3_nt_index_before_a_wider_pam_less_one(self):
        self.assertEqual(self.taken("NNNNNNNN_3_hg38", "NRG_3_hg38"), "NRG_3_hg38")

    def test_before_a_dedicated_wider_one(self):
        self.assertEqual(self.taken("NNGRRT_2_hg38", "NRG_3_hg38"), "NRG_3_hg38")

    def test_before_a_3_nt_pam_less_one(self):
        self.assertEqual(self.taken("NNN_3_hg38", "NRG_3_hg38"), "NRG_3_hg38")

    def test_a_3_nt_pam_less_one_before_a_wider_one(self):
        self.assertEqual(self.taken("NNNNNNNN_3_hg38", "NNN_3_hg38"), "NNN_3_hg38")

    def test_a_variant_index_comes_first_still(self):
        self.assertEqual(self.taken("NRG_3_hg38", "NNNNNNNN_3_hg38+ds"), "NNNNNNNN_3_hg38+ds")
        shutil.rmtree(os.path.join(self.wd, "genome_library", "NRG_3_hg38"))
        self.index("NRG_3_hg38+ds")
        self.assertEqual(ct.require_prebuilt_index(), "NRG_3_hg38+ds")

    def test_a_5_prime_3_nt_index_after_a_3_prime_one(self):
        self.index("ATT_3_hg38", "23bp-ATT-X 26 -3\n")
        self.index("NRG_3_hg38")
        self.assertEqual(ct.require_prebuilt_index(), "NRG_3_hg38")

    def test_one_with_no_bulges_last(self):
        # the first web search's: the smoke searches with one bulge of each type
        self.assertEqual(self.taken("NGG_1_hg38", "NNNNNNNN_3_hg38"), "NNNNNNNN_3_hg38")

    def test_the_only_one_is_taken_whatever_it_is(self):
        self.assertEqual(self.taken("NNNNNNNN_3_hg38"), "NNNNNNNN_3_hg38")


class WhatTheSmokeSearchesWith(_SmokeDir):
    """The smoke wrote 'N'*20 + the name's PAM on the 3' side, and the 23-nt guide, whatever
    the index: 'guide is 23 nt but the PAM expects 28', a fix the user cannot apply, on any
    PAM wider than 3 nt, and a search on the wrong side of a 5' index, which finds nothing
    and passed under --chrom. The failed run also left that PAM file in PAMs/."""

    def inputs(self, name, record=None, guide=None):
        self.index(name, record, guide)
        pam, guide_file, search, prune, five = ct._smoke_inputs(name, name.split("_")[0])
        with open(pam) as fh:
            line = fh.read().strip()
        with open(guide_file) as fh:
            seq = fh.read().strip()
        return line, seq, search, prune, five, pam

    def assertNoPamWritten(self):
        shipped = os.listdir(os.path.join(os.path.dirname(HERE), "PAMs"))
        self.assertEqual(sorted(os.listdir(os.path.join(self.wd, "PAMs"))), sorted(shipped))

    def test_the_published_wider_pam_less_index_is_searched_for_nrr_with_prune(self):
        # it has no .pam_build, as the published archives have none; NRR is the default
        # PAM (it was NRG before 2.9.0)
        line, seq, search, prune, five, pam = self.inputs("NNNNNNNN_3_hg38")
        self.assertEqual(line, "N" * 20 + "NRR 3")
        self.assertEqual(seq, "CTAACAGTTGCTTTTATCACNNN")
        self.assertEqual((search, prune, five), ("NRR", " --prune-pam", False))
        self.assertEqual(pam, os.path.join(self.wd, "PAMs", "20bp-NRR-SpCas9.txt"))
        self.assertNoPamWritten()

    def test_a_3_nt_index_keeps_the_shipped_file_and_the_guide(self):
        for name in ("NRG_3_hg38", "NGG_3_hg38"):
            line, seq, search, prune, five, pam = self.inputs(name)
            token = name[:3]
            self.assertEqual(line, "N" * 20 + f"{token} 3")
            self.assertEqual(seq, "CTAACAGTTGCTTTTATCACNNN")
            self.assertEqual((search, prune, five), (token, "", False))
            self.assertEqual(pam, os.path.join(self.wd, "PAMs", f"20bp-{token}-SpCas9.txt"))
        self.assertNoPamWritten()

    def test_a_3_nt_pam_less_index_is_searched_as_it_was(self):
        line, seq, search, prune, five, pam = self.inputs("NNN_3_hg38")
        self.assertEqual((line, seq, prune), ("N" * 23 + " 3", "CTAACAGTTGCTTTTATCACNNN", ""))
        self.assertEqual(pam, os.path.join(self.wd, "PAMs", "20bp-NNN-NO-PAM.txt"))

    def test_a_wider_3_prime_pam_pads_the_guide_to_it(self):
        line, seq, search, prune, five, pam = self.inputs("NNGRRT_3_hg38")
        self.assertEqual(line, "N" * 20 + "NNGRRT 6")
        self.assertEqual(seq, "CTAACAGTTGCTTTTATCACNNNNNN")
        self.assertEqual((search, prune, five), ("NNGRRT", "", False))
        self.assertTrue(os.path.basename(pam).endswith("-smoke.txt"))
        self.assertNoPamWritten()

    def test_a_2_nt_pam_pads_two(self):
        line, seq, *_ = self.inputs("NG_2_hg38")
        self.assertEqual((line, seq), ("N" * 20 + "NG 2", "CTAACAGTTGCTTTTATCACNN"))

    def test_a_5_prime_index_is_searched_on_its_own_side(self):
        line, seq, search, prune, five, pam = self.inputs("TTTV_2_hg38", "23bp-TTTV-Cas12a 27 -4\n",
                                                          guide=23)
        self.assertEqual(line, "TTTV" + "N" * 20 + " -4")
        self.assertEqual(seq, "NNNNCTAACAGTTGCTTTTATCAC")
        self.assertEqual((search, prune, five), ("TTTV", "", True))
        self.assertNoPamWritten()

    def test_a_short_5_prime_index_too(self):
        line, seq, _, _, five, _ = self.inputs("TTN_2_hg38", "23bp-TTN-X 26 -3\n")
        self.assertEqual((line, seq, five), ("TTN" + "N" * 20 + " -3", "NNNCTAACAGTTGCTTTTATCAC", True))

    def test_a_5_prime_pam_less_index_is_not_pruned(self):
        line, seq, search, prune, five, _ = self.inputs("NNNNNNNN_3_hg38", "28bp-N-5p 28 -8\n")
        self.assertEqual(line, "N" * 28 + " -8")
        self.assertEqual(seq, "NNNNNNNNCTAACAGTTGCTTTTATCAC")
        self.assertEqual((search, prune, five), ("NNNNNNNN", "", True))

    def test_a_dedicated_index_with_no_record_takes_its_side_from_the_pam_files(self):
        # a search builds a dedicated index and records nothing on it
        line, _, _, _, five, _ = self.inputs("TTTV_2_hg38")
        self.assertEqual((line, five), ("TTTV" + "N" * 20 + " -4", True))

    def test_one_a_failed_run_left_a_3_prime_file_for_does_not_turn_it_round(self):
        with open(os.path.join(self.wd, "PAMs", "20bp-TTTV-SpCas9.txt"), "w") as fh:
            fh.write("N" * 20 + "TTTV 4\n")
        _, _, _, _, five, _ = self.inputs("TTTV_2_hg38")
        self.assertTrue(five)

    def test_a_wider_dedicated_index_of_unknown_side_is_refused_in_one_line(self):
        self.index("GGTT_2_hg38")
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit) as stop:
            ct._smoke_inputs("GGTT_2_hg38", "GGTT")
        self.assertEqual(stop.exception.code, 2)
        self.assertIn("cannot tell which side of the guide the PAM of index 'GGTT_2_hg38'",
                      stderr.getvalue())
        self.assertEqual(len(stderr.getvalue().splitlines()), 1)

    def test_a_shorter_3_nt_index_with_no_record_is_taken_for_3_prime(self):
        # NAA_3_hg38 has none, and no shipped PAM file holds NAA (NRR stood in here until
        # its file shipped, with 2.9.0)
        line, _, _, _, five, _ = self.inputs("NAA_3_hg38")
        self.assertEqual((line, five), ("N" * 20 + "NAA 3", False))
        self.assertNoPamWritten()

    def test_an_index_built_for_a_shorter_guide_gets_the_guide_it_fits(self):
        line, seq, *_ = self.inputs("NGG_3_hg38", guide=18)
        self.assertEqual((line, seq), ("N" * 18 + "NGG 3", "AACAGTTGCTTTTATCACNNN"))

    def test_one_built_for_a_longer_guide_gets_the_20_nt_example(self):
        line, seq, *_ = self.inputs("NGG_2_hg38", guide=23)
        self.assertEqual((line, seq), ("N" * 20 + "NGG 3", "CTAACAGTTGCTTTTATCACNNN"))

    def test_a_pam_file_of_the_users_is_used_not_overwritten(self):
        mine = os.path.join(self.wd, "PAMs", "20bp-NRG-Mine.txt")
        os.remove(os.path.join(self.wd, "PAMs", "20bp-NRG-SpCas9.txt"))
        with open(mine, "w") as fh:
            fh.write("N" * 20 + "NRG 3\n")
        *_, pam = self.inputs("NRG_3_hg38")
        self.assertEqual(pam, mine)

    def test_a_file_it_writes_is_gone_when_the_smoke_ends(self):
        # as complete-test runs it, through the whole smoke, with the search mocked
        self.index("NNGRRT_3_hg38")
        os.makedirs(os.path.join(self.wd, "Genomes", "hg38"))
        seen = []
        with mock.patch.object(ct, "_search_in_smoke",
                               side_effect=lambda *a, **k: seen.append(a[5])), \
                contextlib.redirect_stderr(io.StringIO()):
            ct._run_prebuilt_index_smoke(2, False)
        self.assertTrue(seen[0].endswith("20bp-NNGRRT-smoke.txt"), seen)
        self.assertFalse(os.path.exists(seen[0]))
        self.assertEqual([n for n in os.listdir(self.wd) if n.startswith(ct.SMOKEPAMPREFIX)], [])


class ZeroHitWithTheIndexsOwnGeometry(unittest.TestCase):
    """An SpCas9 guide can have no off-target in a whole genome with another nuclease's PAM,
    on its index's own geometry: an empty result, and a pass. With an SpCas9-class PAM it is
    still a failure of the smoke, which always finds the example's off-targets."""

    def smoke(self, *inputs, chrom="all"):
        wd = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, wd, ignore_errors=True)
        cwd = os.getcwd()
        os.chdir(wd)
        self.addCleanup(os.chdir, cwd)
        res = os.path.join(wd, "Results", ct.COMPLETETESTRESDIR)
        os.makedirs(res)
        with open(os.path.join(res, "log.txt"), "w") as fh:
            fh.write("No off-targets found\n")
        stderr = io.StringIO()
        with mock.patch.object(ct.subprocess, "call", return_value=0), \
                contextlib.redirect_stderr(stderr):
            ct._search_in_smoke(*inputs, 2, False, chrom)
        return stderr.getvalue()

    def inputs(self, name, token, search, prune, five):
        return (name, token, "hg38", False, "Genomes/hg38", "p.txt", "g.txt", search, prune, five)

    def test_a_5_prime_index_passes(self):
        out = self.smoke(*self.inputs("TTTV_2_hg38", "TTTV", "TTTV", "", True))
        self.assertIn("complete-test PASSED", out)
        self.assertIn("holds none of its off-targets with PAM TTTV", out)

    def test_a_wider_3_prime_pam_passes(self):
        self.assertIn("complete-test PASSED",
                      self.smoke(*self.inputs("NNGRRT_2_hg38", "NNGRRT", "NNGRRT", "", False)))

    def test_nrg_on_a_whole_genome_does_not(self):
        with self.assertRaises(SystemExit) as stop:
            self.smoke(*self.inputs("NRG_3_hg38", "NRG", "NRG", "", False))
        self.assertEqual(stop.exception.code, 1)

    def test_nor_nrg_under_prune(self):
        with self.assertRaises(SystemExit) as stop:
            self.smoke(*self.inputs("NNNNNNNN_3_hg38", "NNNNNNNN", "NRG", " --prune-pam", False))
        self.assertEqual(stop.exception.code, 1)


class ChromView(unittest.TestCase):
    def setUp(self):
        self._cwd = os.getcwd()
        self.wd = os.path.realpath(tempfile.mkdtemp())
        _make_workdir(self.wd)
        os.chdir(self.wd)
        env = dict(os.environ)
        env.pop("CRISPRME_ALLOW_ONDEMAND_BUILD", None)
        self._env = mock.patch.dict(os.environ, env, clear=True)
        self._env.start()

    def tearDown(self):
        self._env.stop()
        os.chdir(self._cwd)
        shutil.rmtree(self.wd)

    def smoke(self, chrom, search=None):
        """Run the smoke for chrom with complete-search stubbed by search(cmd, cwd), which
        by default writes an integrated_results file and returns 0: the whole-genome
        search's subprocess.call, and _search_in_view, which runs a one-chromosome one
        from its view. Returns the (cmd, cwd, tree of cwd) of each call, what went to
        stderr, and the SystemExit code (None when it returned)."""
        calls = []

        def run(cmd, cwd):
            calls.append((cmd, cwd, _tree(cwd) if cwd else None))
            return (search or self._result)(cmd, cwd)

        def call(cmd, shell=False, cwd=None):
            self.assertTrue(shell)
            return run(cmd, cwd)

        err = io.StringIO()
        code = None
        # create=True: the tests also run against a commit before _search_in_view, which
        # ran the one-chromosome search with subprocess.call too
        with mock.patch.object(ct.subprocess, "call", side_effect=call), \
             mock.patch.object(ct, "_search_in_view", side_effect=run, create=True), \
             contextlib.redirect_stderr(err):
            try:
                ct.run_crisprme_test(chrom, "1000G", 2, False)
            except SystemExit as e:
                code = e.code
        return calls, err.getvalue(), code

    def _result(self, cmd, cwd):
        """complete-search with off-targets: an integrated_results table."""
        res = os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR)
        os.makedirs(res, exist_ok=True)
        open(os.path.join(res, "x_integrated_results.tsv"), "w").close()
        return 0

    def _empty(self, cmd, cwd):
        """complete-search's zero-hit end, as submit_job has it: exit 0, no
        integrated_results table, and the marker the web status page reads in log.txt."""
        res = os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR)
        os.makedirs(res, exist_ok=True)
        with open(os.path.join(res, "log.txt"), "w") as fh:
            fh.write("Post-analysis\tEnd\tSat\nNo off-targets found\nJob\tEnd\tSat\n")
        return 0

    def test_the_search_runs_from_a_view_holding_only_the_chromosome(self):
        """Its genome folder holds chr22.fa alone; each index the search could resolve,
        chr22's trees alone (not chr2's, nor the chr22_ contig's, nor the sidecars); the
        INDELS companion, fakechr22's sub-index; the indel log folder, chr22's log; VCFs/
        nothing; everything else is linked as it is. Another genome's index, a locked or
        hidden one and another dataset's are left out."""
        calls, _, code = self.smoke("chr22")
        self.assertIsNone(code)
        self.assertEqual(len(calls), 1)
        cmd, cwd, tree = calls[0]
        self.assertIsNotNone(cwd, "complete-search must run from the view")
        self.assertEqual(os.path.dirname(cwd), self.wd)
        self.assertTrue(os.path.basename(cwd).startswith(ct.CHROMVIEWPREFIX), cwd)
        self.assertEqual(set(tree), set(_VIEW_CHR22))
        for rel, kind in _VIEW_CHR22.items():
            if kind is None:
                self.assertIsNone(tree[rel], rel)
            else:  # a symlink to the same path in the working directory: nothing copied
                self.assertEqual(tree[rel], ("link", os.path.join(self.wd, rel)), rel)
        self.assertIn(f" --genome {os.path.join(cwd, 'Genomes', 'toy')} ", cmd)
        self.assertIn(f" --output {os.path.join(self.wd, 'Results', 'crisprme-test-out')} ", cmd)
        self.assertIn(" --index-path genome_library ", cmd)
        self.assertIn(" --vcf list_vcf.txt --samplesID list_samplesID.txt ", cmd)

    def test_it_is_the_whole_genome_search_but_for_the_genome_and_output_paths(self):
        whole, _, _ = self.smoke("all")
        self.assertIsNone(whole[0][1], "the whole-genome smoke runs from the working dir")
        shutil.rmtree(os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR))
        one, _, _ = self.smoke("chr22")
        cmd, view, _ = one[0]
        self.assertIsNotNone(view, "complete-search must run from the view")
        self.assertEqual(
            cmd.replace(os.path.join(view, "Genomes", "toy"),
                        os.path.join(self.wd, "Genomes", "toy"))
               .replace(os.path.join(self.wd, "Results", "crisprme-test-out"),
                        "crisprme-test-out"),
            whole[0][0],
        )

    def test_a_chromosome_without_an_indel_log_gets_an_empty_one(self):
        """The contig is in both indexes but has no indel log: the view gives it an empty
        one, which the indel search reads as no indel there, as the whole-genome run
        does by leaving it out of its list. Without it the indel search would find no
        chromosome at all to search, and fail."""
        calls, _, code = self.smoke("chr22_KI270731v1_random")
        self.assertIsNone(code)
        _, cwd, tree = calls[0]
        self.assertIsNotNone(cwd, "complete-search must run from the view")
        logs = {r: v for r, v in tree.items()
                if r.startswith("Dictionaries/log_indels_toy_panel/")}
        self.assertEqual(logs, {
            "Dictionaries/log_indels_toy_panel/logchr22_KI270731v1_random.txt": ("file", ""),
        })
        self.assertEqual(
            sorted(r for r in tree if r.startswith("genome_library/")),
            ["genome_library/NGG_2_toy",
             "genome_library/NGG_2_toy+toy_panel",
             "genome_library/NGG_2_toy+toy_panel/NGG_chr22_KI270731v1_random_1.bin",
             "genome_library/NGG_2_toy+toy_panel_INDELS",
             "genome_library/NGG_2_toy/NGG_chr22_KI270731v1_random_1.bin",
             "genome_library/NRG_3_toy"],
        )

    def test_the_view_is_removed_and_the_data_left_as_it_was(self):
        before = _tree(self.wd)
        calls, _, code = self.smoke("chr22")
        self.assertIsNone(code)
        self.assertIsNotNone(calls[0][1], "complete-search must run from the view")
        self.assertFalse(os.path.lexists(calls[0][1]))
        self.assertEqual(_views(self.wd), [])
        after = _tree(self.wd)
        written = {"sg1617_test_guide.txt",
                   "Results/crisprme-test-out",
                   "Results/crisprme-test-out/x_integrated_results.tsv"}
        self.assertEqual(set(after) - set(before), written)
        self.assertEqual({r: after[r] for r in before}, before)

    def test_a_5_prime_pam_less_index_keeps_its_record_in_the_view(self):
        # The search takes an all-N index with no record for a 3' build, and refuses it for
        # a 5' PAM: an index built for one failed under --chrom ("Genome enrichment failed"),
        # and passed on a 3-nt PAM with the search on the wrong side. Its .binview, which the
        # record brings, goes with the view.
        for name in ("NGG_2_toy", "NGG_2_toy+toy_panel", "NRG_3_toy"):
            shutil.rmtree(name.join(("genome_library/", "")))
        shutil.rmtree("genome_library/NGG_2_toy+toy_panel_INDELS")
        os.makedirs("genome_library/NNN_2_toy")
        with open("genome_library/NNN_2_toy/NNN_chr22_1.bin", "w") as fh:
            fh.write("t")
        with open("genome_library/NNN_2_toy/.pam_build", "w") as fh:
            fh.write("23bp-NNN-5p 23 -3\n")
        os.makedirs("genome_library/NNN_3_toy")
        with open("genome_library/NNN_3_toy/NNN_chr22_1.bin", "w") as fh:
            fh.write("t")
        with open("genome_library/NNN_3_toy/.pam_build", "w") as fh:
            fh.write("20bp-NNN-NO-PAM 23 3\n")

        def search(cmd, cwd):
            os.makedirs(os.path.join(self.wd, "Results", "crisprme-test-out", ".binview"))
            return self._result(cmd, cwd)

        calls, _, code = self.smoke("chr22", search)
        self.assertIsNone(code)
        tree = calls[0][2]
        self.assertIn("genome_library/NNN_2_toy/.pam_build", tree)
        self.assertNotIn("genome_library/NNN_3_toy/.pam_build", tree)  # 3': as without it
        self.assertFalse(os.path.exists(os.path.join("Results", "crisprme-test-out", ".binview")))

    def _assert_removed_on(self, search, want_code=None, want_exc=None):
        seen = []

        def run(cmd, cwd):
            self.assertIsNotNone(cwd, "complete-search must run from the view")
            self.assertTrue(os.path.isdir(cwd))
            seen.append(cwd)
            return search(cmd, cwd)

        if want_exc:
            with self.assertRaises(want_exc):
                self.smoke("chr22", run)
        else:
            _, _, code = self.smoke("chr22", run)
            self.assertEqual(code, want_code)
        self.assertEqual(len(seen), 1)
        self.assertTrue(os.path.basename(seen[0]).startswith(ct.CHROMVIEWPREFIX))
        self.assertFalse(os.path.lexists(seen[0]))
        self.assertEqual(_views(self.wd), [])

    def test_the_view_is_removed_when_the_search_fails(self):
        self._assert_removed_on(lambda cmd, cwd: 1, want_code=1)

    def test_the_view_is_removed_on_ctrl_c(self):
        def interrupted(cmd, cwd):
            raise KeyboardInterrupt
        self._assert_removed_on(interrupted, want_exc=KeyboardInterrupt)

    def test_the_view_is_removed_on_sigterm_and_sighup(self):
        for sig in (signal.SIGTERM, signal.SIGHUP):
            with self.subTest(signal=sig.name):
                # at its default, whatever this test run inherited: an ignored one stays so
                before = signal.signal(sig, signal.SIG_DFL)
                try:
                    def stopped(cmd, cwd):
                        os.kill(os.getpid(), sig)
                        time.sleep(5)  # the handler runs first
                        return 0
                    self._assert_removed_on(stopped, want_code=128 + sig)
                    self.assertIs(signal.getsignal(sig), signal.SIG_DFL)
                finally:
                    signal.signal(sig, before)

    def test_a_signal_ignored_at_entry_stays_ignored(self):
        """nohup ignores SIGHUP: complete_test.py must leave it so, and the search it
        starts then inherits it (exec keeps an ignored signal, and resets a handled one
        to its default), so a hangup stops neither, as without --chrom. The other
        signal is still handled."""
        for sig, other in ((signal.SIGHUP, signal.SIGTERM), (signal.SIGTERM, signal.SIGHUP)):
            with self.subTest(signal=sig.name):
                before = signal.signal(sig, signal.SIG_IGN)
                before_other = signal.signal(other, signal.SIG_DFL)
                try:
                    seen = []

                    def hung_up(cmd, cwd):
                        seen.append((signal.getsignal(sig), signal.getsignal(other)))
                        os.kill(os.getpid(), sig)
                        time.sleep(0.2)
                        return self._result(cmd, cwd)

                    calls, err, code = self.smoke("chr22", hung_up)
                    self.assertIsNone(code, err)
                    self.assertIsNotNone(calls[0][1], "complete-search must run from the view")
                    self.assertEqual(seen, [(signal.SIG_IGN, ct._exit_on_signal)])
                    self.assertIs(signal.getsignal(sig), signal.SIG_IGN)
                    self.assertIn("complete-test PASSED", err)
                    self.assertEqual(_views(self.wd), [])
                finally:
                    signal.signal(sig, before)
                    signal.signal(other, before_other)
                    shutil.rmtree(os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR),
                                  ignore_errors=True)

    def test_a_chromosome_without_an_off_target_passes_with_an_empty_result(self):
        """The example guide has no off-target there (chrM, most unplaced contigs):
        complete-search ends at its zero-hit short-circuit, exit 0 and no
        integrated_results table, as the whole-genome test then reports none there
        either. The test passes, saying so in one line, rather than failing on the
        missing table, with crisprme.py's traceback pointing at an empty error log."""
        calls, err, code = self.smoke("chr22", self._empty)
        self.assertIsNone(code, err)
        self.assertIsNotNone(calls[0][1], "complete-search must run from the view")
        last = err.splitlines()[-1]
        self.assertTrue(last.startswith(
            "complete-test PASSED: searched the example guide against "
            "'NGG_2_toy+toy_panel' on chr22, which holds none of its off-targets"), last)
        self.assertNotIn("ERROR", err)
        self.assertEqual(_views(self.wd), [])

    def test_an_empty_result_still_fails_without_the_zero_hit_end(self):
        """No table, and no zero-hit marker in log.txt: the search did not end as an
        empty one, so the test still fails."""
        _, err, code = self.smoke("chr22", lambda cmd, cwd: 0)
        self.assertEqual(code, 1)
        self.assertIn("ERROR: complete-test produced no integrated_results.tsv", err)

    def test_an_empty_whole_genome_result_still_fails(self):
        """Without --chrom the smoke is unchanged: the example guide has off-targets in
        the genome, and a run without a table fails, zero-hit marker or not."""
        calls, err, code = self.smoke("all", self._empty)
        self.assertIsNone(calls[0][1])
        self.assertEqual(code, 1)
        self.assertIn("ERROR: complete-test produced no integrated_results.tsv", err)

    def test_a_view_left_half_built_is_removed(self):
        real = os.symlink
        made = []

        def failing(src, dst):
            if len(made) == 5:
                raise OSError("no space left on device")
            made.append(dst)
            real(src, dst)

        with mock.patch.object(ct.os, "symlink", side_effect=failing), \
             self.assertRaises(OSError):
            self.smoke("chr22")
        self.assertEqual(len(made), 5)
        self.assertEqual(_views(self.wd), [])

    def _refused(self, chrom):
        calls, err, code = self.smoke(chrom)
        self.assertEqual(calls, [], "nothing may be searched")
        self.assertEqual(code, 2)
        lines = err.splitlines()
        self.assertEqual(len(lines), 1, err)
        self.assertTrue(lines[0].startswith(f"ERROR: complete-test --chrom {chrom}: "), lines[0])
        self.assertEqual(_views(self.wd), [])
        for written in ("sg1617_test_guide.txt",):
            self.assertFalse(os.path.exists(written), written)
        return lines[0]

    def test_a_chromosome_neither_has_is_refused_naming_those_available(self):
        line = self._refused("chrQ")
        self.assertIn("Genomes/toy and the index 'NGG_2_toy+toy_panel' have no chrQ", line)
        self.assertTrue(line.endswith("; available: chr2, chr22, chr22_KI270731v1_random"), line)

    def test_a_chromosome_the_index_lacks_is_refused(self):
        line = self._refused("chrY")
        self.assertIn(": the index 'NGG_2_toy+toy_panel' has no chrY; ", line)

    def test_a_chromosome_the_genome_lacks_is_refused(self):
        line = self._refused("chrM")
        self.assertIn(": Genomes/toy has no chrM; ", line)

    def test_a_long_list_names_the_primary_chromosomes(self):
        lib = os.path.join(self.wd, "genome_library", "NGG_2_toy+toy_panel")
        primary = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]
        for c in primary + [f"chr1_KI2707{n:02d}v1_random" for n in range(20)]:
            open(os.path.join(self.wd, "Genomes", "toy", f"{c}.fa"), "w").close()
            open(os.path.join(lib, f"NGG_{c}_1.bin"), "w").close()
        line = self._refused("chrQ")
        self.assertTrue(line.endswith(f"; available: {', '.join(primary)} and 21 more"), line)

    def _fake_indel_genome(self, *chroms):
        """Genomes/toy+toy_panel_INDELS as add-variants writes it, a header-only
        fake-indel FASTA for each of chroms (no indel there to index)."""
        for chrom in chroms:
            path = os.path.join(self.wd, "Genomes", "toy+toy_panel_INDELS",
                                f"fake_toy_panel_{chrom}", f"fake{chrom}.fa")
            os.makedirs(os.path.dirname(path))
            with open(path, "w") as fh:
                fh.write(f">fake{chrom}\n")

    def test_a_variant_index_without_its_indels_index_gets_its_fake_indel_genome(self):
        """build-index-only builds no INDELS index for an indel-free dataset, and a build
        stopped before its INDELS index leaves the variant index without it: the search
        builds it from Genomes/<ref>+<vcf>_INDELS, or finds no indel there to index, as
        the whole-genome test does. The view had no such folder, so the one-chromosome
        smoke failed, saying to re-download an index that does not exist. It now holds
        it, with the chromosome's fake-indel genome alone."""
        shutil.rmtree(os.path.join(self.wd, "genome_library", "NGG_2_toy+toy_panel_INDELS"))
        self._fake_indel_genome("chr2", "chr22")
        calls, err, code = self.smoke("chr22")
        self.assertIsNone(code, err)
        _, cwd, tree = calls[0]
        self.assertIsNotNone(cwd, "complete-search must run from the view")
        rel = "Genomes/toy+toy_panel_INDELS/fake_toy_panel_chr22"
        self.assertEqual({r: v for r, v in tree.items() if r.startswith("Genomes/toy+")}, {
            "Genomes/toy+toy_panel_INDELS": None,
            rel: ("link", os.path.join(self.wd, rel)),
        })
        self.assertEqual(sorted(r for r in tree if r.startswith("genome_library/NGG_2_toy+")),
                         ["genome_library/NGG_2_toy+toy_panel",
                          "genome_library/NGG_2_toy+toy_panel/NGG_chr22_1.bin"])
        self.assertEqual(_views(self.wd), [])

    def test_the_fake_indel_genomes_mark_goes_into_the_view_with_it(self):
        """Its INDELS index in place, the search does not read the fake-indel genome,
        and its mark does not stop the smoke; were the search to build an INDELS index
        from it, the view's copy of the folder carries the mark, which the search then
        refuses, as it does outside."""
        self._fake_indel_genome("chr22")
        mark = os.path.join(self.wd, "Genomes", "toy+toy_panel_INDELS", ".building")
        os.makedirs(mark)
        calls, err, code = self.smoke("chr22")
        self.assertIsNone(code, err)
        _, cwd, tree = calls[0]
        self.assertIsNotNone(cwd, "complete-search must run from the view")
        self.assertEqual(tree["Genomes/toy+toy_panel_INDELS/.building"], ("link", mark))

    def _unfinished(self, marked, chrom="chr22"):
        """The smoke for chrom with each of the working directory's folders marked
        (relative paths): the ERROR lines on stderr and the exit code, having checked
        that nothing was searched, and that neither a view nor a results folder was
        left."""
        for rel in marked:
            os.makedirs(os.path.join(self.wd, rel, ".building"))
        calls, err, code = self.smoke(chrom)
        self.assertEqual(calls, [], "nothing may be searched")
        self.assertEqual(_views(self.wd), [])
        self.assertFalse(os.path.lexists(os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR)))
        return [l for l in err.splitlines() if l.startswith("ERROR")], code

    def test_a_part_marked_unfinished_is_refused_before_the_view_is_made(self):
        """complete-search refuses a search with a dataset that build-index-only has
        left a part of unfinished, but looks for the marks where it runs from: from the
        view, which holds no Genomes/<ref>+<vcf>, the one-chromosome smoke passed where
        every real search of the dataset is refused, and for the Dictionaries parts it
        named their paths in the view, gone once the run ended, which its rm -rf then
        pointed to. Each part is now looked at in the working directory first, and
        named by its real path, before anything is searched or written."""
        shutil.rmtree(os.path.join(self.wd, "genome_library", "NGG_2_toy+toy_panel_INDELS"))
        self._fake_indel_genome("chr22")
        rerun = "run build-index-only again with toy_panel"
        cases = (
            ("Genomes/toy+toy_panel", "the enrichment of toy_panel", f"{rerun}, which enriches "
             "it again"),
            ("Genomes/toy+toy_panel_INDELS", "the SNP codes of the toy_panel fake-indel genome",
             f"{rerun}, which finishes them"),
            ("Dictionaries/registry_toy_panel", "the dict-less tiers of toy_panel",
             f"{rerun} and its --samplesID, which finishes them, or remove them to search "
             f"without them (rm -rf '{self.wd}/Dictionaries/registry_toy_panel' "
             f"'{self.wd}/Dictionaries/genotypes_toy_panel')"),
            ("Dictionaries/indel_genotypes_toy_panel", "the phased indel genotype store of "
             "toy_panel", f"{rerun}, which finishes it, or remove it to search without it "
             f"(rm -rf '{self.wd}/Dictionaries/indel_genotypes_toy_panel')"),
        )
        for rel, what, then in cases:
            with self.subTest(part=rel):
                errors, code = self._unfinished([rel])
                self.assertEqual(code, 1)
                self.assertEqual(errors, [
                    f"ERROR: {os.path.join(self.wd, rel)} is unfinished: a build-index-only run "
                    f"is writing {what} there, or stopped before it was done, and the search "
                    f"would use what is there as if it were whole. Let that run finish, or "
                    f"{then}, then run complete-test again."])
                shutil.rmtree(os.path.join(self.wd, rel, ".building"))
        with self.subTest(part="all of them"):
            errors, code = self._unfinished([rel for rel, _, _ in cases])
            self.assertEqual((code, len(errors)), (1, 4), errors)
        for rel, _, _ in cases:
            shutil.rmtree(os.path.join(self.wd, rel, ".building"))
        # a Dictionaries entry that is a symlink (as the view links each): its real path
        real = os.path.join(self.wd, "elsewhere", "registry_toy_panel")
        os.makedirs(os.path.dirname(real))
        os.rename(os.path.join(self.wd, "Dictionaries", "registry_toy_panel"), real)
        os.symlink(real, os.path.join(self.wd, "Dictionaries", "registry_toy_panel"))
        os.makedirs(os.path.join(real, ".building"))
        with self.subTest(part="a linked registry"):
            calls, err, code = self.smoke("chr22")
            self.assertEqual((calls, code), ([], 1))
            self.assertIn(f"ERROR: {real} is unfinished: ", err)
            self.assertIn(f"(rm -rf '{real}' '{self.wd}/Dictionaries/genotypes_toy_panel')", err)

    def test_what_a_search_does_not_read_is_not_refused(self):
        """With the SNP+indel co-occurrence off, a search does not read the indel
        genotype store, and a reference smoke reads no part of a dataset: their marks do
        not stop it, as they do not stop complete-search. (Nor does the fake-indel
        genome's, with its INDELS index in place: test_the_fake_indel_genomes_mark_goes_
        into_the_view_with_it.)"""
        store = "Dictionaries/indel_genotypes_toy_panel"
        for value in ("0", ""):
            with self.subTest(part=store, CRISPRME_INDEL_SNP=value):
                os.makedirs(os.path.join(self.wd, store, ".building"), exist_ok=True)
                with mock.patch.dict(os.environ, {"CRISPRME_INDEL_SNP": value}):
                    calls, err, code = self.smoke("chr22")
                self.assertIsNone(code, err)
                self.assertEqual(len(calls), 1)
                shutil.rmtree(os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR))
        with self.subTest(part="a reference smoke"):
            os.remove(os.path.join(self.wd, "list_vcf.txt"))
            calls, err, code = self.smoke("chr22")
            self.assertIsNone(code, err)
            self.assertEqual(len(calls), 1)

    def test_a_dataset_the_list_names_twice_gets_one_folder_in_the_view(self):
        """A --vcf list may name a dataset twice, or as toy_panel, toy_panel/ and
        VCFs/toy_panel, which the search reads as one. The view made the dataset's
        fake-indel genome folder once a line, and the second one failed with a
        FileExistsError traceback wherever Genomes/<ref>+<vcf>_INDELS is (any dataset
        built here), and a part marked unfinished was refused once a line. The view,
        and the refusal, are now those of the dataset named once."""
        self._fake_indel_genome("chr22")
        for case in ("its INDELS index in place", "no INDELS index"):
            with self.subTest(case=case):
                if case == "no INDELS index":
                    shutil.rmtree(os.path.join(self.wd, "genome_library",
                                               "NGG_2_toy+toy_panel_INDELS"))
                trees = []
                for listed in ("toy_panel\n", "toy_panel\ntoy_panel/\nVCFs/toy_panel\n"):
                    with open(os.path.join(self.wd, "list_vcf.txt"), "w") as fh:
                        fh.write(listed)
                    calls, err, code = self.smoke("chr22")
                    self.assertIsNone(code, err)
                    self.assertEqual(len(calls), 1)
                    trees.append(calls[0][2])
                    shutil.rmtree(os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR))
                self.assertEqual(trees[1], trees[0])
                self.assertIn("Genomes/toy+toy_panel_INDELS/fake_toy_panel_chr22", trees[1])
                self.assertEqual(_views(self.wd), [])
        with self.subTest(case="a part marked unfinished"):
            errors, code = self._unfinished(["Genomes/toy+toy_panel"])
            self.assertEqual(code, 1)
            self.assertEqual(errors, [
                f"ERROR: {self.wd}/Genomes/toy+toy_panel is unfinished: a build-index-only run "
                "is writing the enrichment of toy_panel there, or stopped before it was done, "
                "and the search would use what is there as if it were whole. Let that run "
                "finish, or run build-index-only again with toy_panel, which enriches it again, "
                "then run complete-test again."])


# complete-search's stand-in for RealSignals: it and a child shell, as the real search
# has a tree of processes, note their PIDs, and any SIGTERM or SIGHUP they get, with
# whether the view they run from is still there; released, it writes the table
# complete-search would in its --output folder, and ends.
_STUB = r"""#!/bin/sh
b=@BIN@
echo $$ >> "$b/pids"
trap 'echo "$$ TERM $([ -d "$PWD" ] && echo in-view || echo view-gone)" >> "$b/got"; exit 143' TERM
trap 'echo "$$ HUP $([ -d "$PWD" ] && echo in-view || echo view-gone)" >> "$b/got"; exit 129' HUP
if [ "$1" = child ]; then
	while [ ! -e "$b/release" ]; do sleep 0.05; done
	exit 0
fi
"$0" child &
while [ "$(grep -c '' "$b/pids")" -lt 2 ]; do sleep 0.05; done
pwd -P > "$b/started.tmp" && mv "$b/started.tmp" "$b/started"
while [ ! -e "$b/release" ]; do sleep 0.05; done
wait
for a; do [ "$prev" = --output ] && out=$a; prev=$a; done
mkdir -p "$out" && : > "$out/x_integrated_results.tsv"
"""

_HAS_PROC = os.path.isdir("/proc/self")


def _running(pid):
    """Whether process pid still runs, from /proc (Linux): a zombie, ended but not yet
    reaped, does not."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


class RealSignals(unittest.TestCase):
    """complete_test.py run as a process, as crisprme.py runs it, in its own session with
    SIGINT at its default (as under a terminal), with crisprme.py on PATH the stub
    above. The view must go however the run is stopped: Ctrl-C (SIGINT to the whole
    process group), or SIGTERM or SIGHUP to complete_test.py alone, which it passes on
    to the search, waiting for it to end before the view goes. Under nohup a hangup
    stops neither."""

    def setUp(self):
        self.wd = os.path.realpath(tempfile.mkdtemp())
        self.bin = os.path.realpath(tempfile.mkdtemp())
        _make_workdir(self.wd)
        self.mark = os.path.join(self.bin, "started")
        stub = os.path.join(self.bin, "crisprme.py")
        with open(stub, "w") as fh:
            fh.write(_STUB.replace("@BIN@", self.bin))
        os.chmod(stub, 0o755)

    def tearDown(self):
        shutil.rmtree(self.wd)
        shutil.rmtree(self.bin)

    def _lines(self, name):
        try:
            with open(os.path.join(self.bin, name)) as fh:
                return fh.read().split("\n")[:-1]
        except FileNotFoundError:
            return []

    def _run(self, stop, nohup=False):
        """Start complete_test.py --chrom chr22 and, once the stub runs in the view,
        call stop(proc); return its exit code, its stderr, the stub's PIDs, those of
        them still running as it ended, and what the stub noted of the signals."""
        env = {k: v for k, v in os.environ.items() if k != "CRISPRME_ALLOW_ONDEMAND_BUILD"}
        env["PATH"] = self.bin + os.pathsep + env.get("PATH", "")

        def preexec():
            # as under a terminal, or nohup, whatever this test run inherited
            for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
                signal.signal(s, signal.SIG_DFL)
            if nohup:
                signal.signal(signal.SIGHUP, signal.SIG_IGN)

        with tempfile.TemporaryFile() as err:
            proc = subprocess.Popen(
                [sys.executable, os.path.join(HERE, "complete_test.py"), "chr22", "1000G",
                 "1", "False"],
                cwd=self.wd, env=env, stdout=subprocess.DEVNULL, stderr=err,
                start_new_session=True, preexec_fn=preexec,
            )
            try:
                deadline = time.time() + 30
                while not os.path.exists(self.mark) and proc.poll() is None \
                        and time.time() < deadline:
                    time.sleep(0.05)
                err.seek(0)
                self.assertTrue(os.path.exists(self.mark), err.read().decode())
                with open(self.mark) as fh:
                    view = fh.read().strip()
                self.assertEqual(os.path.dirname(view), self.wd)
                self.assertTrue(os.path.basename(view).startswith(ct.CHROMVIEWPREFIX), view)
                self.assertTrue(os.path.isdir(view))
                stop(proc)
                proc.wait(timeout=30)
                pids = [int(p) for p in self._lines("pids")]
                left = [p for p in pids if _running(p)]
            finally:
                with contextlib.suppress(ProcessLookupError, PermissionError):
                    os.killpg(proc.pid, signal.SIGKILL)  # the stub, if still there
                proc.wait()
            err.seek(0)
            stderr = err.read().decode()
        self.assertFalse(os.path.lexists(view))
        self.assertEqual(_views(self.wd), [])
        return proc.returncode, stderr, pids, left, self._lines("got")

    def test_ctrl_c(self):
        rc, *_ = self._run(lambda proc: os.killpg(proc.pid, signal.SIGINT))
        self.assertNotEqual(rc, 0)

    def _alone(self, sig):
        """sig to complete_test.py alone, as `kill <its PID>` sends it: complete-search
        and the processes under it, in the same process group, get none, so it passes
        sig on to them and waits for them to end before it removes the view, lest they
        go on without it."""
        rc, stderr, pids, left, got = self._run(lambda proc: proc.send_signal(sig))
        self.assertEqual(rc, 128 + sig, stderr)
        self.assertEqual(len(pids), 2)
        if _HAS_PROC:  # where complete_test.py finds the search's processes
            self.assertEqual(left, [], "the search outlived complete-test")
            self.assertEqual(sorted(got), sorted(f"{p} {sig.name[3:]} in-view" for p in pids))

    def test_sigterm_to_complete_test_alone(self):
        self._alone(signal.SIGTERM)

    def test_sighup_to_complete_test_alone(self):
        self._alone(signal.SIGHUP)

    def test_a_hangup_under_nohup(self):
        """nohup ignores SIGHUP, and so must complete_test.py and the search it starts:
        a hangup (SIGHUP to the whole process group) stops neither, and the test
        passes once the search ends, as without --chrom."""
        def hang_up(proc):
            os.killpg(proc.pid, signal.SIGHUP)
            time.sleep(0.5)
            self.assertIsNone(proc.poll(), "complete-test ended on an ignored hangup")
            open(os.path.join(self.bin, "release"), "w").close()

        rc, stderr, pids, left, got = self._run(hang_up, nohup=True)
        self.assertEqual(rc, 0, stderr)
        self.assertIn("complete-test PASSED", stderr)
        self.assertEqual(got, [], "the search was sent the ignored SIGHUP")
        self.assertEqual(left, [])


class StopSearch(unittest.TestCase):
    """_stop_search on a real search that ignores the signal passed on: it is killed,
    with the processes under it, once SEARCHSTOPWAIT seconds are over."""

    @unittest.skipUnless(_HAS_PROC, "complete_test.py finds a search's processes in /proc")
    def test_a_search_outlasting_the_wait_is_killed(self):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp)
        pids = os.path.join(tmp, "pids")
        proc = subprocess.Popen(
            # ignored here, the signal is ignored in every process under it too
            f"trap '' TERM; sh -c 'echo $$ >> {pids}; while :; do sleep 0.05; done' & "
            f"echo $$ >> {pids}; while :; do sleep 0.05; done",
            shell=True,
        )
        self.addCleanup(lambda: proc.poll() is None and proc.kill())

        def noted():
            with contextlib.suppress(FileNotFoundError), open(pids) as fh:
                return [int(p) for p in fh.read().split()]
            return []

        deadline = time.time() + 10
        while time.time() < deadline and len(noted()) < 2:
            time.sleep(0.05)
        tree = noted()
        self.assertEqual(len(tree), 2)
        t0 = time.monotonic()
        with mock.patch.object(ct, "SEARCHSTOPWAIT", 0.5):
            ct._stop_search(proc, signal.SIGTERM)
        took = time.monotonic() - t0
        self.assertGreaterEqual(took, 0.5)
        self.assertLess(took, 5)
        self.assertEqual(proc.returncode, -signal.SIGKILL)
        # _stop_search waits for proc alone: a process under it that was just sent
        # SIGKILL can still show as running for a moment, under load
        deadline = time.time() + 5
        while time.time() < deadline and any(_running(p) for p in tree):
            time.sleep(0.05)
        self.assertEqual([p for p in tree if _running(p)], [])


# --------------------------------------------------------------------------- #
# Load crisprme.py DEFINITIONS only (no CLI dispatch), Bio.Seq + scientific
# absents stubbed -- identical to test_validate_test._load_crisprme.
# --------------------------------------------------------------------------- #
def _load_crisprme():
    crisprme_py = os.path.join(os.path.dirname(HERE), "crisprme.py")
    sys.path.insert(0, HERE)
    with open(crisprme_py) as fh:
        src = fh.read()
    marker = "\nif len(sys.argv) < 2:"
    cut = src.index(marker)  # KeyError-loud if the dispatch shape ever changes
    mod = types.ModuleType("crisprme_under_test_complete")
    mod.__file__ = crisprme_py
    code = compile(src[:cut], crisprme_py, "exec")
    # Stubs exist ONLY so crisprme.py's top-level imports resolve while we exec its
    # definitions; they are removed again in the finally, so they cannot shadow the real
    # modules for a test collected later.
    _inserted = []
    if "Bio.Seq" not in sys.modules:
        _bio = types.ModuleType("Bio")
        _bio_seq = types.ModuleType("Bio.Seq")
        _bio_seq.Seq = object
        _bio.Seq = _bio_seq
        if "Bio" not in sys.modules:
            sys.modules["Bio"] = _bio
            _inserted.append("Bio")
        sys.modules["Bio.Seq"] = _bio_seq
        _inserted.append("Bio.Seq")
    for _m in ("pandas", "scipy", "sklearn", "matplotlib", "seaborn",
               "statsmodels", "intervaltree", "CRISPR_BULGE_score"):
        if _m in sys.modules:
            continue
        try:
            __import__(_m)
        except Exception:
            sys.modules[_m] = mock.MagicMock()
            _inserted.append(_m)
    try:
        exec(code, mod.__dict__)
    finally:
        for _m in _inserted:
            sys.modules.pop(_m, None)
    return mod


class TheCliPassesTheRefusalOn(unittest.TestCase):
    """complete_test.py refuses a --chrom the genome or the index lacks, or to run over an
    earlier run's results, with exit 2, having said why in one line; crisprme.py exits 2
    too, adding nothing: neither the OSError traceback (whose 'See
    Results/crisprme-test-out/log_error.txt' names a file no search wrote) nor the line of
    a failure, which would name the earlier run's error logs. Any other failure ends with
    exit 1 and one line, which names the search's error log only where one was written."""

    ARGV = ["crisprme.py", "complete-test", "--chrom", "chrQ"]

    @classmethod
    def setUpClass(cls):
        cls.crisprme = _load_crisprme()

    def _cli(self, code, argv=None, passed="chrQ 1000G 4 False"):
        with mock.patch.object(self.crisprme, "input_args", argv or self.ARGV), \
             mock.patch.object(self.crisprme, "warn_low_memory", lambda *a, **k: None), \
             mock.patch.object(self.crisprme, "subprocess") as msub:
            msub.call.return_value = code
            try:
                self.crisprme.complete_test_crisprme()
            finally:
                (cmd,), _ = msub.call.call_args
                self.assertTrue(cmd.endswith(f"complete_test.py {passed}"), cmd)

    def test_a_refused_chrom_exits_2_without_a_traceback(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            self._cli(2)
        self.assertEqual((cm.exception.code, err.getvalue()), (2, ""))

    def test_a_failed_test_ends_in_one_line_naming_the_log_it_left(self):
        """It raised an OSError, a traceback ending in 'See
        Results/crisprme-test-out/log_error.txt', which a search refused before it
        starts never writes, and which the benchmark mode writes under another name. The
        logs of the other mode are not this run's: complete_test.py runs beside them."""
        smoke = "Results/crisprme-test-out/log_error.txt"
        bench = "Results/crisprme-test-out_cas12a_hbg/log_error.txt"
        cases = (
            (False, {}, "see the messages above"),
            (False, {smoke: ""}, "see the messages above"),
            (False, {smoke: "ERROR: x\n"}, f"see {smoke} for details"),
            (True, {"Results/crisprme-test-out_cas9_sg1617/log_error_no_check.txt": "",
                    bench: "ERROR: y\n"}, f"see {bench} for details"),
            (False, {bench: "ERROR: y\n"}, "see the messages above"),
            (True, {smoke: "ERROR: x\n", bench: "ERROR: y\n"}, f"see {bench} for details"),
        )
        for benchmarks, logs, says in cases:
            with self.subTest(benchmarks=benchmarks, logs=sorted(logs)):
                wd = tempfile.mkdtemp()
                self.addCleanup(shutil.rmtree, wd)
                for rel, text in logs.items():
                    os.makedirs(os.path.dirname(os.path.join(wd, rel)), exist_ok=True)
                    with open(os.path.join(wd, rel), "w") as fh:
                        fh.write(text)
                env = {k: v for k, v in os.environ.items()
                       if k != "CRISPRME_ALLOW_ONDEMAND_BUILD"}
                if benchmarks:
                    env["CRISPRME_ALLOW_ONDEMAND_BUILD"] = "1"
                cwd = os.getcwd()
                os.chdir(wd)
                err = io.StringIO()
                try:
                    with mock.patch.dict(os.environ, env, clear=True), \
                         contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as cm:
                        self._cli(1)
                finally:
                    os.chdir(cwd)
                self.assertEqual(cm.exception.code, 1)
                self.assertEqual(err.getvalue(),
                                 f"\nCRISPRme complete test encountered an error: {says}\n")

    def test_vcf_dataset_without_the_benchmark_switch_is_said_on_stdout(self):
        """--vcf_dataset is the benchmark mode's, which CRISPRME_ALLOW_ONDEMAND_BUILD
        alone selects: without it, complete-test ran the prebuilt-index smoke with no
        word of it, as it does in Docker when the variable is set in front of `docker
        run`, on the host. One line on stdout now says so; the run, and its exit status,
        are as before."""
        argv = ["crisprme.py", "complete-test", "--chrom", "chr22", "--vcf_dataset", "1000G"]
        note = ("NOTE [complete-test]: --vcf_dataset is for the benchmark mode, whose output "
                "validate-test checks, and that mode runs only with "
                "CRISPRME_ALLOW_ONDEMAND_BUILD=1 set (Docker: pass it with -e, as in docker run "
                "-e CRISPRME_ALLOW_ONDEMAND_BUILD=1 ...). It is not "
                "set here, so complete-test runs the example guide against a prebuilt index "
                "instead.")
        for switch, given, said in ((None, argv, [note]), ("1", argv, []),
                                    (None, argv[:4], [])):
            with self.subTest(CRISPRME_ALLOW_ONDEMAND_BUILD=switch, argv=given[2:]):
                env = {k: v for k, v in os.environ.items()
                       if k != "CRISPRME_ALLOW_ONDEMAND_BUILD"}
                if switch:
                    env["CRISPRME_ALLOW_ONDEMAND_BUILD"] = switch
                out = io.StringIO()
                with mock.patch.dict(os.environ, env, clear=True), \
                     contextlib.redirect_stdout(out):
                    self._cli(0, given, "chr22 1000G 4 False")
                self.assertEqual(out.getvalue().splitlines(), said)


if __name__ == "__main__":
    unittest.main()
