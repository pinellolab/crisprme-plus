#!/usr/bin/env python3
"""What the CI workflows run, read from their files.

A job that is cancelled before its end, never starts on the change at hand, or skips
what it is there to test, passes for nothing, and the release rule is all-green CI.

Read as text: the CI env has no YAML parser (and YAML 1.1 reads the `on` key as true).
"""
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
