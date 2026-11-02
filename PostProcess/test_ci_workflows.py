#!/usr/bin/env python3
"""What the CI workflows and the documented self-test run, read from their files.

A job that is cancelled before its end, never starts on the change at hand, or skips
what it is there to test, passes for nothing, and the release rule is all-green CI.

Read as text: the CI env has no YAML parser (and YAML 1.1 reads the `on` key as true).
"""
import json
import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORKFLOWS = REPO / ".github" / "workflows"


def workflow(name):
    return (WORKFLOWS / name).read_text()


@unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
class TheUnitTestsJobHasTheTimeToFinish(unittest.TestCase):
    """Its limit was 10 minutes. The hosted job took 523 s of those 600 at aa9c54b, and
    side by side under one load the steps of the trial integration of the fixes after it
    took 413 s against 334 s, which predicts about 640 s: the job would be cancelled, and
    a tag needs it green."""

    def test_its_limit_is_at_least_twenty_minutes(self):
        limit = re.search(r"^\s*timeout-minutes:\s*(\d+)\s*$", workflow("unit-tests.yml"), re.M)
        self.assertGreaterEqual(int(limit.group(1)), 20)


class TheDocumentedSelfTestRunsEveryBenchmark(unittest.TestCase):
    """Every documented self-test and the CI job set CRISPRME_SKIP_HEAVY=1, on the claim
    that the bulge-2 index builds need more than 16 GB, so each left out what it could have
    checked. Measured: 3.66 GB peak RSS, and a pass in 16 GiB, 4 CPUs, no swap. The variable
    stays as an opt-out for very small machines, so it may be named, not put in a command."""

    SURFACES = ("README.md", "docs/DOCKER_QUICKSTART.md", "crisprme.py", "install_from_source.sh",
                "PostProcess/validate.py", "test/benchmark/README.md")
    IN_A_COMMAND = re.compile(r"CRISPRME_ALLOW_ONDEMAND_BUILD=1[\s\\>]+CRISPRME_SKIP_HEAVY"
                              r"|-e CRISPRME_SKIP_HEAVY")
    OLD_CLAIM = re.compile(r"need more than[\s>]+16 GB|OOMs? a (?:standard )?16 ?GB"
                           r"|too large for this runner")

    def text(self, name):
        return (REPO / name).read_text()

    def benchmarks(self):
        return json.loads(self.text("test/benchmark/benchmarks.json"))["benchmarks"]

    def test_no_documented_command_sets_the_opt_out(self):
        for name in self.SURFACES:
            with self.subTest(surface=name):
                self.assertIsNone(self.IN_A_COMMAND.search(self.text(name)))

    @unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
    def test_nor_does_the_ci_job(self):
        self.assertNotIn("CRISPRME_SKIP_HEAVY", workflow("validate-benchmarks.yml"))

    def test_the_expected_output_is_every_benchmark_passed(self):
        expected = f"{len(self.benchmarks())} passed, 0 failed, 0 skipped"
        for name, count in (("README.md", 3), ("docs/DOCKER_QUICKSTART.md", 1)):
            with self.subTest(surface=name):
                said = re.findall(r"# expect: (.+)", self.text(name))
                self.assertEqual(said, [expected] * count)

    def test_the_more_than_16_gb_claim_is_gone(self):
        names = self.SURFACES + ("PostProcess/complete_test.py", "test/benchmark/benchmarks.json")
        if WORKFLOWS.is_dir():
            names += (".github/workflows/validate-benchmarks.yml",)
        for name in names:
            with self.subTest(surface=name):
                self.assertIsNone(self.OLD_CLAIM.search(self.text(name)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
