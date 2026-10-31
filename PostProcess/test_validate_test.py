#!/usr/bin/env python3
"""validate-test where there is no complete-test benchmark output to check.

validate-test compares the off-targets that complete-test's BENCHMARK mode wrote
(CRISPRME_ALLOW_ONDEMAND_BUILD=1 ... --vcf_dataset 1000G, as validate-benchmarks.yml
runs it) with the committed brute-force ground truth, and starts by reading the
vcf.config.test.txt that mode writes. Run anywhere else -- a fresh folder, or after
the default prebuilt-index smoke, which writes no such file -- it ended in two Python
tracebacks: validate.py's FileNotFoundError, then crisprme.py's OSError. It now says
in one line what to run first, and exits 1.

validate.py is run as a subprocess, as crisprme.py runs it. crisprme.py's
validate_test is exec'd without the CLI dispatch (the truncate-and-exec loader of
test_bulge_cap), with its subprocess call stubbed.
"""
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
