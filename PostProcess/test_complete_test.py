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
is searched. ChromView checks the view in-process on a small working directory, with
complete-search stubbed; RealSignals sends real signals to complete_test.py run as a
process, with a stub crisprme.py on PATH; the CLI test checks that crisprme.py passes the
refusal on without a traceback.
"""

import contextlib
import io
import os
import shutil
import signal
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
    def _run(self, allow_ondemand, chrom="chr22"):
        """Run run_crisprme_test with every side-effecting dependency mocked; return the
        list of complete-search command strings + a flag of whether require_prebuilt_index
        was consulted. The keyword arguments of each call are kept in self.kwargs."""
        calls = []
        self.kwargs = []
        env = dict(os.environ)
        env.pop("CRISPRME_ALLOW_ONDEMAND_BUILD", None)
        env.pop("CRISPRME_SKIP_HEAVY", None)
        if allow_ondemand:
            env["CRISPRME_ALLOW_ONDEMAND_BUILD"] = "1"
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
             mock.patch.object(ct, "load_benchmarks", return_value=_REGISTRY), \
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
            "--max-total-edits 4 --merge 3 --pam 20bp-NRG-SpCas9.txt --guide "
            "sg1617.guide.txt --vcf list_vcf.txt --samplesID list_samplesID.txt "
            "--index-path genome_library --output crisprme-test-out --thread 4  --ci-cd-test"
        ])
        self.assertTrue(self.kwargs[0].get("shell"))
        self.assertIsNone(self.kwargs[0].get("cwd"))


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
        by default writes an integrated_results file and returns 0. Returns the (cmd,
        cwd, tree of cwd) of each call, what went to stderr, and the SystemExit code
        (None when it returned)."""
        calls = []

        def call(cmd, shell=False, cwd=None):
            self.assertTrue(shell)
            calls.append((cmd, cwd, _tree(cwd) if cwd else None))
            if search:
                return search(cmd, cwd)
            res = os.path.join(self.wd, "Results", ct.COMPLETETESTRESDIR)
            os.makedirs(res, exist_ok=True)
            open(os.path.join(res, "x_integrated_results.tsv"), "w").close()
            return 0

        err = io.StringIO()
        code = None
        with mock.patch.object(ct.subprocess, "call", side_effect=call), \
             contextlib.redirect_stderr(err):
            try:
                ct.run_crisprme_test(chrom, "1000G", 2, False)
            except SystemExit as e:
                code = e.code
        return calls, err.getvalue(), code

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
        written = {"sg1617_test_guide.txt", "PAMs/20bp-NGG-SpCas9.txt",
                   "Results/crisprme-test-out",
                   "Results/crisprme-test-out/x_integrated_results.tsv"}
        self.assertEqual(set(after) - set(before), written)
        self.assertEqual({r: after[r] for r in before}, before)

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
                before = signal.getsignal(sig)

                def stopped(cmd, cwd):
                    os.kill(os.getpid(), sig)
                    time.sleep(5)  # the handler runs first
                    return 0
                self._assert_removed_on(stopped, want_code=128 + sig)
                self.assertIs(signal.getsignal(sig), before)

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
        for written in ("sg1617_test_guide.txt", os.path.join("PAMs", "20bp-NGG-SpCas9.txt")):
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


class RealSignals(unittest.TestCase):
    """complete_test.py run as a process, as crisprme.py runs it, in its own session with
    SIGINT at its default (as under a terminal), with crisprme.py on PATH a stub that
    notes where complete-search was started and waits. The view must go however the
    run is stopped: Ctrl-C (SIGINT to the whole process group), or SIGTERM or SIGHUP
    to complete_test.py alone."""

    def setUp(self):
        self.wd = os.path.realpath(tempfile.mkdtemp())
        self.bin = tempfile.mkdtemp()
        _make_workdir(self.wd)
        self.mark = os.path.join(self.bin, "started")
        stub = os.path.join(self.bin, "crisprme.py")
        with open(stub, "w") as fh:
            fh.write(f'#!/bin/sh\npwd -P > "{self.mark}.tmp" && mv "{self.mark}.tmp" "{self.mark}"\n'
                     "exec sleep 60\n")
        os.chmod(stub, 0o755)

    def tearDown(self):
        shutil.rmtree(self.wd)
        shutil.rmtree(self.bin)

    def _stop(self, sig, group):
        env = {k: v for k, v in os.environ.items() if k != "CRISPRME_ALLOW_ONDEMAND_BUILD"}
        env["PATH"] = self.bin + os.pathsep + env.get("PATH", "")
        with tempfile.TemporaryFile() as err:
            proc = subprocess.Popen(
                [sys.executable, os.path.join(HERE, "complete_test.py"), "chr22", "1000G",
                 "1", "False"],
                cwd=self.wd, env=env, stdout=subprocess.DEVNULL, stderr=err,
                start_new_session=True,
                preexec_fn=lambda: signal.signal(signal.SIGINT, signal.SIG_DFL),
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
                if group:
                    os.killpg(proc.pid, sig)
                else:
                    proc.send_signal(sig)
                proc.wait(timeout=30)
            finally:
                with contextlib.suppress(ProcessLookupError, PermissionError):
                    os.killpg(proc.pid, signal.SIGKILL)  # the stub's sleep, if still there
                proc.wait()
        self.assertFalse(os.path.lexists(view))
        self.assertEqual(_views(self.wd), [])
        return proc.returncode

    def test_ctrl_c(self):
        self.assertNotEqual(self._stop(signal.SIGINT, group=True), 0)

    def test_sigterm_to_complete_test_alone(self):
        self.assertEqual(self._stop(signal.SIGTERM, group=False), 128 + signal.SIGTERM)

    def test_sighup_to_complete_test_alone(self):
        self.assertEqual(self._stop(signal.SIGHUP, group=False), 128 + signal.SIGHUP)


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
    """complete_test.py refuses a --chrom the genome or the index lacks with exit 2,
    having said why in one line; crisprme.py exits 2 too, without the OSError traceback
    (whose 'See Results/crisprme-test-out/log_error.txt' names a file no search wrote).
    Any other failure still raises it."""

    ARGV = ["crisprme.py", "complete-test", "--chrom", "chrQ"]

    @classmethod
    def setUpClass(cls):
        cls.crisprme = _load_crisprme()

    def _cli(self, code):
        with mock.patch.object(self.crisprme, "input_args", self.ARGV), \
             mock.patch.object(self.crisprme, "warn_low_memory", lambda *a, **k: None), \
             mock.patch.object(self.crisprme, "subprocess") as msub:
            msub.call.return_value = code
            try:
                self.crisprme.complete_test_crisprme()
            finally:
                (cmd,), _ = msub.call.call_args
                self.assertTrue(cmd.endswith("complete_test.py chrQ 1000G 4 False"), cmd)

    def test_a_refused_chrom_exits_2_without_a_traceback(self):
        with self.assertRaises(SystemExit) as cm:
            self._cli(2)
        self.assertEqual(cm.exception.code, 2)

    def test_a_failed_test_still_raises(self):
        with self.assertRaises(OSError):
            self._cli(1)


if __name__ == "__main__":
    unittest.main()
