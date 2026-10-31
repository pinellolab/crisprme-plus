#!/usr/bin/env python3
"""validate-test where there is no complete-test benchmark output to check, and on a
chromosome that a benchmark's brute-force reference does not cover.

validate-test compares the off-targets that complete-test's BENCHMARK mode wrote
(CRISPRME_ALLOW_ONDEMAND_BUILD=1 ... --vcf_dataset 1000G, as validate-benchmarks.yml
runs it) with the committed brute-force ground truth, and starts by reading the
vcf.config.test.txt that mode writes. Run anywhere else -- a fresh folder, or after
the default prebuilt-index smoke, which writes no such file -- it ended in two Python
tracebacks: validate.py's FileNotFoundError, then crisprme.py's OSError. It now says
in one line what to run first, and exits 1.

That line names the pair for the chromosome asked for, but the Cas12a reference and
both single-"n edits" ones cover chr22 only ("scope" in benchmarks.json): compared on
chr21, or on a whole genome, each site CRISPRme found off chr22 counted as extra and
the benchmark FAILED. Such a benchmark is now skipped there.

validate.py is run as a subprocess, as crisprme.py runs it, and its
run_test_validation in-process with each comparison stubbed. crisprme.py's
validate_test is exec'd without the CLI dispatch (the truncate-and-exec loader of
test_bulge_cap), with its subprocess call stubbed.
"""
import contextlib
import io
import os
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATE = os.path.join(HERE, "validate.py")


# --------------------------------------------------------------------------- #
# Load crisprme.py DEFINITIONS only (no CLI dispatch), Bio.Seq + scientific
# absents stubbed -- identical to test_bulge_cap._load_crisprme.
# --------------------------------------------------------------------------- #
def _load_crisprme():
    crisprme_py = os.path.join(os.path.dirname(HERE), "crisprme.py")
    sys.path.insert(0, HERE)
    with open(crisprme_py) as fh:
        src = fh.read()
    marker = "\nif len(sys.argv) < 2:"
    cut = src.index(marker)  # KeyError-loud if the dispatch shape ever changes
    mod = types.ModuleType("crisprme_under_test_validate")
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


crisprme = _load_crisprme()
import validate  # noqa: E402 -- PostProcess is on sys.path now (_load_crisprme)


class ValidatePyWithoutBenchmarkOutput(unittest.TestCase):
    def test_one_line_on_stderr_and_exit_1(self):
        with tempfile.TemporaryDirectory() as d:  # a folder complete-test never ran in
            out = subprocess.run([sys.executable, VALIDATE, "chr22"], cwd=d,
                                 capture_output=True, text=True)
        self.assertEqual(out.returncode, 1, out.stderr)
        self.assertEqual(out.stdout, "")
        self.assertNotIn("Traceback", out.stderr)
        lines = out.stderr.splitlines()
        self.assertEqual(len(lines), 1, out.stderr)
        # it says what is missing and what to run first: the documented benchmark mode
        for part in ("vcf.config.test.txt", "CRISPRME_ALLOW_ONDEMAND_BUILD=1",
                     "complete-test --chrom chr22 --vcf_dataset 1000G"):
            self.assertIn(part, lines[0])

    def hint(self, chrom):
        """validate.py's one line, run as validate-test runs it for `chrom`."""
        with tempfile.TemporaryDirectory() as d:
            out = subprocess.run([sys.executable, VALIDATE, chrom], cwd=d,
                                 capture_output=True, text=True)
        self.assertEqual(out.returncode, 1, out.stderr)
        lines = out.stderr.splitlines()
        self.assertEqual(len(lines), 1, out.stderr)
        return lines[0]

    def test_the_commands_are_for_the_chromosome_asked_for(self):
        # validate-test --chrom chr21 reads hg38_chr21 reports, which a chr22 run never writes
        line = self.hint("chr21")
        self.assertIn("complete-test --chrom chr21 --vcf_dataset 1000G", line)
        self.assertIn("then crisprme.py validate-test --chrom chr21", line)
        self.assertNotIn("chr22", line)

    def test_without_chrom_they_are_the_documented_chr22_pair(self):
        # validate-test without --chrom ("all") checks a whole-genome run: rather than
        # that, the chr22 self-test, and the validate-test that then has something to check
        line = self.hint("all")
        self.assertIn("complete-test --chrom chr22 --vcf_dataset 1000G", line)
        self.assertIn("then crisprme.py validate-test --chrom chr22", line)


class ABenchmarkOffItsReferenceIsSkipped(unittest.TestCase):
    """After the `complete-test --chrom chr21 --vcf_dataset 1000G` that the line above
    names for chr21, `validate-test --chrom chr21` reported '[cas12a_hbg] FAILED:
    mismatch. CRISPRme sites: 223  Brute-force sites: 0': its reference covers chr22
    only, so every chr21 site counted as extra. A benchmark is now compared only on
    the chromosome its reference covers, and skipped elsewhere."""

    # what the benchmark mode runs with CRISPRME_SKIP_HEAVY=1, as that line has it
    LIGHT = ["cas9_sg1617", "cas12a_hbg"]

    def run_validation(self, chrom, ran):
        """run_test_validation(chrom) on the shipped registry, with CRISPRME_SKIP_HEAVY
        unset as in the validate-test that line names, where the benchmarks in `ran`
        left their output, each comparison stubbed to a match. Returns (exit code,
        stderr lines, the benchmarks compared)."""
        compared = []
        err = io.StringIO()
        env = dict(os.environ)
        env.pop("CRISPRME_SKIP_HEAVY", None)

        def found(name, *_a):
            return (f"{name}.ref.txt", f"{name}.alt.txt") if name in ran else None

        def match(name, *_a):
            compared.append(name)
            return True

        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(validate, "check_variant_dataset"), \
             mock.patch.object(validate, "find_crispritz_targets", side_effect=found), \
             mock.patch.object(validate, "download_brute_force_targets", return_value="bf.tsv"), \
             mock.patch.object(validate, "load_targets", return_value=(None, None)), \
             mock.patch.object(validate, "validate", side_effect=match), \
             contextlib.redirect_stderr(err):
            with self.assertRaises(SystemExit) as cm:
                validate.run_test_validation(chrom)
        return cm.exception.code, err.getvalue().splitlines(), compared

    def test_chr21_compares_only_what_the_references_cover(self):
        code, lines, compared = self.run_validation("chr21", self.LIGHT)
        self.assertEqual(compared, ["cas9_sg1617"])
        self.assertIn("[cas12a_hbg] SKIPPED: its brute-force reference covers chr22 only, "
                      "not chr21.", lines)
        self.assertIn("Validation summary: 1 passed, 0 failed, 3 skipped.", lines)
        self.assertEqual(code, 0)

    def test_so_does_a_whole_genome_run(self):
        # validate-test without --chrom
        code, lines, compared = self.run_validation("all", self.LIGHT)
        self.assertEqual(compared, ["cas9_sg1617"])
        self.assertIn("[cas12a_hbg] SKIPPED: its brute-force reference covers chr22 only, "
                      "not the whole genome.", lines)
        self.assertIn("Validation summary: 1 passed, 0 failed, 3 skipped.", lines)
        self.assertEqual(code, 0)

    def test_chr22_still_compares_every_benchmark_that_ran(self):
        # the documented self-test ("expect: 2 passed, 0 failed, 2 skipped")
        code, lines, compared = self.run_validation("chr22", self.LIGHT)
        self.assertEqual(compared, self.LIGHT)
        self.assertIn("Validation summary: 2 passed, 0 failed, 2 skipped.", lines)
        self.assertEqual(code, 0)
        # and a run of all four, the heavy ones included
        every = [b["name"] for b in validate.load_benchmarks()["benchmarks"]]
        code, lines, compared = self.run_validation("chr22", every)
        self.assertEqual(compared, every)
        self.assertEqual(code, 0)

    def test_validate_py_no_longer_fails_chr21_sites_against_the_chr22_reference(self):
        # validate.py itself and the committed Cas12a reference, on chr21 output of
        # cas12a_hbg's: before, the real comparison counted each of these sites as extra
        registry = validate.load_benchmarks()
        bench = next(b for b in registry["benchmarks"] if b["name"] == "cas12a_hbg")
        th = registry["thresholds"]
        site = "X\tTTTVCCTTGTCAAGGCTATTGGTC\tTTTACCTTGTCAAGGCTATTGGTC\tchr21\t{}\t{}\t+\t0\t0\t0\n"
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "vcf.config.test.txt"), "w") as fout:
                fout.write("hg38_1000G\n")
            tdir = os.path.join(d, validate.targets_dir("cas12a_hbg"))
            os.makedirs(tdir)
            names = validate.target_report_names(bench["pam_name"], th["mm"], th["bDNA"],
                                                 th["bRNA"], "chr21")
            for fname, positions in zip(names, ((5010000, 5020000), (5030000,))):
                with open(os.path.join(tdir, fname), "w") as fout:
                    fout.writelines(site.format(p, p) for p in positions)
            out = subprocess.run([sys.executable, VALIDATE, "chr21"], cwd=d,
                                 capture_output=True, text=True)
        self.assertIn("[cas12a_hbg] SKIPPED: its brute-force reference covers chr22 only, "
                      "not chr21.", out.stderr.splitlines(), out.stderr)
        self.assertNotIn("FAILED", out.stderr)


class TheCliCallerAddsNoTraceback(unittest.TestCase):
    """validate.py has said why on stderr; crisprme.py only passes the failure on."""

    ARGV = ["crisprme.py", "validate-test", "--chrom", "chr22"]

    def test_a_failed_validation_exits_1(self):
        with mock.patch.object(crisprme, "input_args", self.ARGV), \
             mock.patch.object(crisprme, "subprocess") as msub:
            msub.call.return_value = 1
            with self.assertRaises(SystemExit) as cm:
                crisprme.validate_test()
        self.assertEqual(cm.exception.code, 1)
        (cmd,), _ = msub.call.call_args
        self.assertTrue(cmd.endswith("validate.py chr22"), cmd)

    def test_a_passed_validation_returns(self):
        with mock.patch.object(crisprme, "input_args", self.ARGV), \
             mock.patch.object(crisprme, "subprocess") as msub:
            msub.call.return_value = 0
            self.assertIsNone(crisprme.validate_test())


if __name__ == "__main__":
    unittest.main(verbosity=2)
