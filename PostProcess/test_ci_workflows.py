#!/usr/bin/env python3
"""What the CI workflows and the documented self-test run, read from their files.

A job that is cancelled before its end, never starts on the change at hand, or skips
what it is there to test, passes for nothing, and the release rule is all-green CI.

Read as text: the CI env has no YAML parser (and YAML 1.1 reads the `on` key as true).
"""
import ast
import json
import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORKFLOWS = REPO / ".github" / "workflows"


def read(path):
    return Path(path).read_text(encoding="utf-8")


def workflow(name):
    return read(WORKFLOWS / name)


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


def listed(path, patterns):
    """Whether a path filter matches `path` as GitHub matches it: * stops at a slash, ** does
    not, and **/ may match nothing (fnmatch's * crosses slashes). The filters here use no
    other wildcard."""
    for pattern in patterns:
        rx = "".join({"**/": "(?:.*/)?", "**": ".*", "*": "[^/]*"}.get(part, re.escape(part))
                     for part in re.split(r"(\*\*/|\*\*|\*)", pattern))
        if re.fullmatch(rx, path):
            return True
    return False


def function_based_modules():
    """The test modules whose tests are plain functions run by their own __main__."""
    found = []
    for path in sorted((REPO / "PostProcess").glob("test_*.py")):
        body = ast.parse(read(path)).body
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
    release PR touches everything, so this bit later PRs: green with the job never started.
    A replay of both workflows' steps, logging each file read, found more: the report tests
    read assets/ and test/data/samplesIDs/, test_complete_search_flags reads every script in
    the repository, seq_script/ and the paper's included, for a complete-search caller, and
    the web pages import seq_script/ and download_hprc_assembly.py."""

    # read in that replay, by a test the workflow runs, and matched by none of its filters before
    SEEN = {"unit-tests.yml": ("assets/crisprme-logo.svg", "assets/crisprme_bg_report.jpg",
                               "test/data/samplesIDs/samplesIDs.1000G.txt", "environment.yml",
                               "seq_script/merge_panels/build_combined_panel.sh",
                               "download_hprc_assembly.py", "plot_generation_paper/plot_generation.sh",
                               "paper/reproducibility/variant-datasets/hg38_AoU/csv_to_vcf.py"),
            "web-e2e.yml": ("seq_script/convert_pam.py", "seq_script/extract_seq.py",
                            "download_hprc_assembly.py", "PAMs/20bp-NRG-SpCas9.txt")}

    def test_web_e2e_runs_on_postprocess(self):
        self.assertIn("PostProcess/**", pull_request_paths("web-e2e.yml"))

    def test_unit_tests_runs_on_what_its_tests_read(self):
        paths = pull_request_paths("unit-tests.yml")
        for path in ("pages/**", "scripts/**", "docs/**", "README.md", "METHODS.md",
                     "test/benchmark/**", "PAMs/**", "CHANGELOG.md", ".github/workflows/**",
                     "assets/**", "test/data/**", "**.py", "**.sh", "**.yml", "**.yaml"):
            with self.subTest(path=path):
                self.assertIn(path, paths)

    def test_the_two_pages_entries_it_replaces_are_gone(self):
        paths = pull_request_paths("unit-tests.yml")
        self.assertNotIn("pages/history_page.py", paths)
        self.assertNotIn("pages/settings_page.py", paths)

    def test_the_match_is_githubs(self):
        # with fnmatch's, "*.py" would pass for a script in a folder
        self.assertTrue(listed("seq_script/merge_panels/hprc_build.sh", ["**.sh"]))
        self.assertFalse(listed("seq_script/convert_pam.py", ["*.py"]))
        self.assertTrue(listed("docs/hf/card.md", ["docs/**"]))
        self.assertFalse(listed("docs/hf/card.md", ["docs/*"]))
        self.assertFalse(listed("READMEXmd", ["README.md"]))
        self.assertTrue(listed("README.md", ["**/README.md"]))

    def test_a_path_this_module_reads_is_listed_for_it(self):
        # this module reads them, so it must run on a change to any
        paths = pull_request_paths("unit-tests.yml")
        for read in ("docs/RELEASE_CHECKLIST.md", "README.md", "docs/DOCKER_QUICKSTART.md",
                     "install_from_source.sh", "test/benchmark/benchmarks.json",
                     ".github/workflows/web-e2e.yml", "crisprme.py", "PostProcess/validate.py"):
            with self.subTest(path=read):
                self.assertTrue(listed(read, paths))

    def test_a_file_its_steps_were_seen_to_read_is_listed(self):
        for name, files in self.SEEN.items():
            paths = pull_request_paths(name)
            for read in files:
                with self.subTest(workflow=name, path=read):
                    self.assertTrue(listed(read, paths))

    def test_every_file_the_flag_gate_scans_is_listed(self):
        # a pull request adding a complete-search call with a dead flag anywhere must start it
        import test_complete_search_flags as flag_gate
        paths = pull_request_paths("unit-tests.yml")
        scanned = [os.path.relpath(p, flag_gate.REPO) for p in flag_gate.scanned_files()]
        self.assertIn("PostProcess/complete_test.py", scanned)
        self.assertEqual([p for p in scanned if not listed(p, paths)], [])

    def test_web_e2e_runs_on_every_repository_module_the_pages_import(self):
        # the web unit tests and the browser tests load the pages, and so all they import
        paths = pull_request_paths("web-e2e.yml")
        files = set()
        for src in sorted(REPO.glob("pages/*.py")) + [REPO / "app.py", REPO / "index.py"]:
            for node in ast.walk(ast.parse(read(src))):
                if isinstance(node, ast.Import):
                    names = [a.name for a in node.names]
                elif isinstance(node, ast.ImportFrom) and not node.level and node.module:
                    names = [node.module] + [f"{node.module}.{a.name}" for a in node.names]
                else:
                    continue
                for name in names:
                    stem = name.replace(".", "/")
                    files |= {f for f in (f"{stem}.py", f"{stem}/__init__.py") if (REPO / f).is_file()}
        self.assertIn("seq_script/convert_pam.py", files)
        self.assertEqual(sorted(f for f in files if not listed(f, paths)), [])


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

    def guard(self, log):
        """The exit status of the step's check of its log, run as the step runs it."""
        step = self.in_image_step()
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "in-image-tests.log").write_text(log)
            return subprocess.run(
                ["bash", "-e", "-c", step[step.index("if grep"):]], cwd=d, capture_output=True,
                text=True).returncode

    def test_a_failed_test_fails_the_step(self):
        # bash -e alone would not see it through the tee
        self.assertIn("set -o pipefail", self.in_image_step())

    def test_a_skip_there_fails_the_step(self):
        # a test whose tool has gone from the image skips, and would pass as in the light job.
        # unittest gives a reason holding an apostrophe in double quotes.
        for line in ("test_a (m.C.test_a) ... skipped 'bedops (sort-bed/closest-features) not on PATH'",
                     'test_b (m.C.test_b) ... skipped "needs the checkout\'s .github/workflows"',
                     "setUpClass (m.C) ... skipped 'crisprme.py dependencies not installed'"):
            with self.subTest(line=line):
                self.assertEqual(self.guard(f"test_c (m.C.test_c) ... ok\n{line}\n"), 1)

    def test_but_the_two_the_image_cannot_help(self):
        log = ("test_c (m.C.test_c) ... ok\n"
               "test_d (m.C.test_d) ... skipped 'en_US.UTF-8 is not installed here, or collates like C'\n"
               "test_e (m.C.test_e) ... skipped 'no ps here to compare with, as in the image'\n")
        self.assertEqual(self.guard(log), 0)

    def test_every_test_module_is_named_in_a_workflow(self):
        text = "\n".join(read(p) for p in WORKFLOWS.glob("*.yml"))
        unnamed = [p.stem for p in sorted((REPO / "PostProcess").glob("test_*.py"))
                   if not re.search(rf"\b{p.stem}\b", text)]
        self.assertEqual(unnamed, [])

    def test_a_function_based_module_runs_as_a_script(self):
        # `python -m unittest` finds no test in these, and says OK: 0 tests, a false pass
        text = "\n".join(read(p) for p in WORKFLOWS.glob("*.yml"))
        modules = function_based_modules()
        self.assertTrue(modules)
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
        return read(REPO / name)

    def benchmarks(self):
        return json.loads(self.text("test/benchmark/benchmarks.json"))["benchmarks"]

    def test_no_documented_command_sets_the_opt_out(self):
        for name in self.SURFACES:
            with self.subTest(surface=name):
                self.assertIsNone(self.IN_A_COMMAND.search(self.text(name)))

    @unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
    def test_nor_does_the_ci_job(self):
        self.assertNotIn("CRISPRME_SKIP_HEAVY", workflow("validate-benchmarks.yml"))

    @unittest.skipUnless(WORKFLOWS.is_dir(), "needs the checkout's .github/workflows")
    def test_a_benchmark_that_did_not_run_fails_the_ci_job(self):
        # validate-test exits 0 with one skipped: complete-test left it no output, say
        text = workflow("validate-benchmarks.yml")
        start = text.index("- name: validate-test (chr22)")
        step = text[start:text.index("- name:", start + 1)]
        self.assertIn("set -o pipefail", step)
        check = step[step.index("if ! grep"):]
        for summary, status in (("6 passed, 0 failed, 0 skipped.", 0), ("5 passed, 0 failed, 1 skipped.", 1),
                                ("2 passed, 0 failed, 4 skipped.", 1)):
            with self.subTest(summary=summary), tempfile.TemporaryDirectory() as d:
                (Path(d) / "validate-test.log").write_text(f"\nValidation summary: {summary}\nEnjoy CRISPRme!\n")
                got = subprocess.run(["bash", "-e", "-c", check], cwd=d, capture_output=True, text=True)
                self.assertEqual(got.returncode, status, got.stdout)

    def test_the_expected_output_is_every_benchmark_passed(self):
        expected = f"{len(self.benchmarks())} passed, 0 failed, 0 skipped"
        for name in ("README.md", "docs/DOCKER_QUICKSTART.md"):
            with self.subTest(surface=name):
                said = re.findall(r"# expect: (\d+ passed, \d+ failed, \d+ skipped)", self.text(name))
                self.assertTrue(said)
                self.assertEqual(set(said), {expected})

    def test_the_more_than_16_gb_claim_is_gone(self):
        names = self.SURFACES + ("PostProcess/complete_test.py", "test/benchmark/benchmarks.json")
        if WORKFLOWS.is_dir():
            names += (".github/workflows/validate-benchmarks.yml",)
        for name in names:
            with self.subTest(surface=name):
                self.assertIsNone(self.OLD_CLAIM.search(self.text(name)))


class TheReleaseGateLeavesNothingSkipped(unittest.TestCase):
    """docs/RELEASE_CHECKLIST.md's local gate: the commands that run on the release image
    what CI cannot, with the only skips it allows named. Its checks are greps, so what
    they look for has to stay what the tests and validate-test print."""

    def checklist(self):
        return read(REPO / "docs" / "RELEASE_CHECKLIST.md")

    def test_it_expects_every_benchmark(self):
        count = len(json.loads(read(REPO / "test/benchmark/benchmarks.json"))["benchmarks"])
        self.assertIn(f"Validation summary: {count} passed, 0 failed, 0 skipped", self.checklist())

    def test_it_sets_no_opt_out(self):
        self.assertNotRegex(self.checklist(), r"-e CRISPRME_SKIP_HEAVY|--env CRISPRME_SKIP_HEAVY")

    def test_the_skips_it_allows_are_the_skips_the_tests_print(self):
        text = self.checklist()
        allowed = re.search(r"""grep -hE "skipped \['\\"\]" \S+ \| grep -vE "([^"]+)"$""", text,
                            re.M).group(1).split("|")
        self.assertEqual(len(allowed), 3, allowed)
        sources = {"no TensorFlow here": "test_crispr_bulge_fidelity.py",
                   "en_US.UTF-8 is not installed": "test_index_resolution.py",
                   "no ps here to compare": "test_index_resolution.py"}
        self.assertEqual(sorted(allowed), sorted(sources))
        for reason, module in sources.items():
            with self.subTest(skip=reason):
                self.assertIn(reason, read(REPO / "PostProcess" / module))
                self.assertIn(f"`{reason}`", text)
        # the one more an Apptainer run allows
        self.assertIn("no group to give genome_library",
                      read(REPO / "PostProcess/test_build_index_only_atomic.py"))
        self.assertIn("`no group to give genome_library`", text)

    def test_its_grep_catches_any_other_skip(self):
        line = next(ln for ln in self.checklist().splitlines() if ln.startswith('grep -hE "skipped'))
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "test_x.log"
            log.write_text("test_a (m.C.test_a) ... skipped 'no ps here to compare with, as in the image'\n"
                           "test_b (m.C.test_b) ... skipped 'no TensorFlow here: run it in the cbulge env'\n"
                           "test_c (m.C.test_c) ... skipped 'en_US.UTF-8 is not installed here, or...'\n")
            said = subprocess.run(["bash", "-c", line], env={"GATE": d, "PATH": os.environ["PATH"]},
                                  capture_output=True, text=True).stdout
            self.assertEqual(said, "")
            log.write_text('test_d (m.C.test_d) ... skipped "no group to give it but this user\'s own"\n'
                           "test_e (m.C.test_e) ... skipped 'bedops not on PATH'\n")
            said = subprocess.run(["bash", "-c", line], env={"GATE": d, "PATH": os.environ["PATH"]},
                                  capture_output=True, text=True).stdout
            self.assertEqual(len(said.splitlines()), 2, said)

    def test_it_runs_the_function_based_modules_as_scripts(self):
        text = self.checklist()
        for name in function_based_modules():
            with self.subTest(module=name):
                self.assertIn(f"`{name}`", text)
        self.assertIn('if grep -q "^Ran 0 tests"', text)
        self.assertIn('python "$f" >', text)

    def test_it_takes_the_commit_the_image_provisions(self):
        self.assertIn("print(scorer_env.CBULGE_PIN)", self.checklist())

    def test_it_names_every_host_its_commands_reach(self):
        # a maintainer who opens only the hosts it names must be able to pass it: without
        # the Debian mirror step 2's locale tests skip, without PyPI step 5 has no playwright
        text = self.checklist()
        gate = text[text.index("#### The local gate"):text.index("**On an HPC node")]
        needs = re.search(r"It needs Docker, git-lfs, and network to ([^.]+)\.", gate).group(1)
        for command, host in (("IMG=pinellolab/crisprme:", "Docker Hub"),
                              ("-e CRISPRME_RUN_HF_TESTS=1", "Hugging Face"),
                              ("complete-test --chrom chr22 --vcf_dataset 1000G", "EBI"),
                              ("git clone https://github.com/", "GitHub"),
                              ("mcr.microsoft.com/", "MCR"),
                              ("pip install", "PyPI"),
                              ("apt-get install", "Debian package mirror")):
            with self.subTest(host=host):
                self.assertIn(command, gate)
                self.assertIn(host, needs)


if __name__ == "__main__":
    unittest.main(verbosity=2)
