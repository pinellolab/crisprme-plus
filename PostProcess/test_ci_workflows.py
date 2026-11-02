#!/usr/bin/env python3
"""What the CI workflows and the documented self-test run, read from their files.

A job that is cancelled before its end, never starts on the change at hand, or skips
what it is there to test, passes for nothing, and the release rule is all-green CI.

Read as text: the CI env has no YAML parser (and YAML 1.1 reads the `on` key as true).
"""
import ast
import fnmatch
import json
import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORKFLOWS = REPO / ".github" / "workflows"


def workflow(name):
    return (WORKFLOWS / name).read_text()


def pull_request_paths(name):
    """The `on: pull_request: paths:` entries of a workflow."""
    lines = workflow(name).splitlines()
    i = lines.index("    paths:", lines.index("  pull_request:", lines.index("on:")))
    paths = []
    for line in lines[i + 1:]:
        entry = line.strip()
        if not entry or entry.startswith("#"):
            continue
        if not entry.startswith("- "):
            break
        paths.append(entry[2:].strip("\"'"))
    return paths


def function_based_modules():
    """The test modules whose tests are plain functions run by their own __main__."""
    found = []
    for path in sorted((REPO / "PostProcess").glob("test_*.py")):
        body = ast.parse(path.read_text()).body
        if any(isinstance(n, ast.FunctionDef) and n.name.startswith("test_") for n in body):
            found.append(path.stem)
    return found


@unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
class TheUnitTestsJobHasTheTimeToFinish(unittest.TestCase):
    """Its limit was 10 minutes. The hosted job took 523 s of those 600 at aa9c54b, and
    side by side under one load the steps of the trial integration of the fixes after it
    took 413 s against 334 s, which predicts about 640 s: the job would be cancelled, and
    a tag needs it green."""

    def test_its_limit_is_at_least_twenty_minutes(self):
        limit = re.search(r"^\s*timeout-minutes:\s*(\d+)\s*$", workflow("unit-tests.yml"), re.M)
        self.assertGreaterEqual(int(limit.group(1)), 20)


@unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
class APullRequestToWhatAJobReadsStartsIt(unittest.TestCase):
    """A pull request starts a workflow only for the paths it lists. web-e2e.yml did not list
    PostProcess/, though it runs test_web_* and test_results_* from there, on code the pages
    import. unit-tests.yml listed two pages/ files, though its tests read main_page.py and
    the whole folder, and missed scripts/, test/benchmark/ and the docs its tests scan. The
    release PR touches everything, so this bit later PRs: green with the job never started."""

    def test_web_e2e_runs_on_postprocess(self):
        self.assertIn("PostProcess/**", pull_request_paths("web-e2e.yml"))

    def test_unit_tests_runs_on_what_its_tests_read(self):
        paths = pull_request_paths("unit-tests.yml")
        for path in ("pages/**", "scripts/**", "docs/**", "README.md", "METHODS.md",
                     "test/benchmark/**", "PAMs/**", "CHANGELOG.md", ".github/workflows/**"):
            with self.subTest(path=path):
                self.assertIn(path, paths)

    def test_the_two_pages_entries_it_replaces_are_gone(self):
        paths = pull_request_paths("unit-tests.yml")
        self.assertNotIn("pages/history_page.py", paths)
        self.assertNotIn("pages/settings_page.py", paths)

    def test_a_path_this_module_reads_is_listed_for_it(self):
        # this module reads them, so it must run on a change to any
        paths = pull_request_paths("unit-tests.yml")
        for read in ("docs/RELEASE_CHECKLIST.md", "README.md", "docs/DOCKER_QUICKSTART.md",
                     "install_from_source.sh", "test/benchmark/benchmarks.json",
                     ".github/workflows/web-e2e.yml", "crisprme.py", "PostProcess/validate.py"):
            with self.subTest(path=read):
                self.assertTrue(any(fnmatch.fnmatch(read, p.replace("**", "*")) for p in paths))


@unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
class TheTestsThatSkipInTheLightJobRunInCI(unittest.TestCase):
    """unit-tests.yml's job has no BEDOPS, minimap2, impg, Dash or Biopython, so its step for
    post_process.sh ran one test and skipped it, the impg truth tests and the flag-gate
    behaviour tests skipped too, and test_assembly_annotate was in no workflow. All pass in the
    image, where no job ran them."""

    MODULES = ("test_assembly_annotate", "test_assembly_reconcile", "test_complete_search_flags",
               "test_post_process_tempdir", "test_index_resolution")

    def in_image_step(self):
        text = workflow("web-e2e.yml")
        start = text.index("""- name: "Tests that need the image's tools""")
        return text[start:text.index("- name:", start + 1)]

    def test_the_image_runs_them(self):
        step = self.in_image_step()
        self.assertIn("crisprme:web", step)
        self.assertIn("-w /opt/conda/opt/crisprme/PostProcess", step)
        modules = re.search(r"python -m unittest ((?:[\w.]+\s*\\?\s*)+?) -v", step)
        self.assertEqual(sorted(modules.group(1).replace("\\", " ").split()), sorted(self.MODULES))

    def test_a_skip_there_fails_the_step(self):
        # bash -e alone would not see the tests fail through the tee, and a test the image
        # no longer has the tool for would skip, and pass, as they did in the light job
        step = self.in_image_step()
        self.assertIn("set -o pipefail", step)
        self.assertRegex(step, r"""if grep "skipped '" \S+ \\\n\s+\| grep -vE "en_US[^"]*\|no ps here[^"]*"; then""")
        self.assertRegex(step, r"FAIL: .*\n?.*exit 1")

    def test_every_test_module_is_named_in_a_workflow(self):
        text = "\n".join(p.read_text() for p in WORKFLOWS.glob("*.yml"))
        unnamed = [p.stem for p in sorted((REPO / "PostProcess").glob("test_*.py"))
                   if not re.search(rf"\b{p.stem}\b", text)]
        self.assertEqual(unnamed, [])

    def test_a_function_based_module_runs_as_a_script(self):
        # `python -m unittest` finds no test in these, and says OK: 0 tests, a false pass
        text = "\n".join(p.read_text() for p in WORKFLOWS.glob("*.yml"))
        modules = function_based_modules()
        self.assertEqual(len(modules), 6, modules)
        for name in modules:
            with self.subTest(module=name):
                self.assertRegex(text, rf"python {name}\.py\b")
                self.assertNotRegex(text, rf"-m unittest(?:\s+\w+)*\s+{name}\b")

    def test_the_real_hugging_face_download_runs_in_the_smoke_step(self):
        # the job already downloads from there; the unit-tests job has no huggingface_hub
        text = workflow("validate-benchmarks.yml")
        start = text.index("- name: Smoke")
        step = text[start:text.index("- name:", start + 1)]
        self.assertIn("-e CRISPRME_RUN_HF_TESTS=1", step)
        self.assertRegex(step, r"python -m unittest test_crisprme_hf test_hf_downloads\n")


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
