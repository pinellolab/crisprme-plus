#!/usr/bin/env python3
"""install_from_source.sh refuses macOS before it does anything else.

The source build is Linux-only: environment.yml has packages with no macOS build, and
CRISPRitz's buildTST needs GCC's parallel STL, which no macOS toolchain provides. The
script used to start regardless and fail later, on a solver or compiler error. It now
stops at once on Darwin, with one line pointing to Docker.

Run with a stub `uname` on PATH and no conda env, so that without the guard the script
stops at its next check (the activated env), never at the build; a stub `git` would
fail loudly if it ever got that far.

What the script leaves in the checkout is ignored by git, and a pull request that
changes only the script or the .gitignore still runs this module in CI.
"""
import ast
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "install_from_source.sh"
WORKFLOW = SCRIPT.parent / ".github" / "workflows" / "unit-tests.yml"


def search_data_folders():
    """crisprme.py's CRISPRMEDIRS, the data folders a search creates in the folder it
    runs in. Read from the source: importing crisprme.py runs its command line."""
    tree = ast.parse((SCRIPT.parent / "crisprme.py").read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id == "CRISPRMEDIRS":
                    return ast.literal_eval(node.value)
    raise AssertionError("CRISPRMEDIRS not found as a module-level literal in crisprme.py")


def pull_request_paths():
    """The `on: pull_request: paths:` entries of unit-tests.yml, read as text: the CI env
    has no YAML parser (and YAML 1.1 reads the `on` key as true)."""
    lines = WORKFLOW.read_text().splitlines()
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


class TheSourceBuildRefusesMacOS(unittest.TestCase):
    def run_script(self, system):
        with tempfile.TemporaryDirectory() as d:
            for name, body in (("uname", f"echo {system}"), ("git", "echo git >&2; exit 97")):
                stub = Path(d) / name
                stub.write_text(f"#!/bin/sh\n{body}\n")
                stub.chmod(0o755)
            env = {k: v for k, v in os.environ.items() if k != "CONDA_PREFIX"}
            env["PATH"] = f"{d}{os.pathsep}{env.get('PATH', '')}"
            return subprocess.run(["bash", str(SCRIPT)], env=env, capture_output=True,
                                  text=True)

    def test_darwin_stops_with_one_line_pointing_to_docker(self):
        out = self.run_script("Darwin")
        self.assertEqual(out.returncode, 1, out.stderr)
        self.assertEqual(out.stdout, "")
        lines = out.stderr.splitlines()
        self.assertEqual(len(lines), 1, out.stderr)
        self.assertIn("Docker", lines[0])
        self.assertIn("README 1.2", lines[0])

    def test_linux_goes_on_to_the_conda_env_check(self):
        out = self.run_script("Linux")
        self.assertNotEqual(out.returncode, 0)
        self.assertIn("CONDA_PREFIX", out.stderr)
        self.assertNotIn("Docker", out.stderr)


@unittest.skipUnless(shutil.which("git") and (SCRIPT.parent / ".gitignore").is_file(),
                     "needs git, and the checkout's .gitignore beside the script")
class WhatTheInstallerLeavesInTheCheckoutIsIgnored(unittest.TestCase):
    """The script runs `crisprme.py scorer-env create` and `check` in the checkout, and
    they save their state to <checkout>/Annotations/.scorer_env.json. Annotations/ is
    one of the data folders a search creates in the folder it runs in (crisprme.py's
    CRISPRMEDIRS), and the only one of them the checkout's .gitignore left out, so
    `git status` showed it untracked after every install. Checked against the repo's
    .gitignore, in a scratch repository."""

    def ignored(self, *paths):
        """Which of `paths` git ignores, created in a repository with that .gitignore."""
        with tempfile.TemporaryDirectory() as d:
            subprocess.run(["git", "init", "-q", d], check=True)
            shutil.copy(SCRIPT.parent / ".gitignore", d)
            for path in paths:
                (Path(d) / path).parent.mkdir(parents=True, exist_ok=True)
                (Path(d) / path).write_text("")
            out = subprocess.run(["git", "-C", d, "check-ignore", *paths],
                                 capture_output=True, text=True)
        return out.stdout.split()

    def test_the_scorer_env_state_is_ignored(self):
        self.assertEqual(self.ignored("Annotations/.scorer_env.json"),
                         ["Annotations/.scorer_env.json"])

    def test_a_new_file_in_any_data_folder_a_search_creates_is_ignored(self):
        new = [f"{d}/new.txt" for d in search_data_folders()]
        self.assertEqual(sorted(self.ignored(*new)), sorted(new))

    def test_a_new_file_beside_the_tracked_test_annotations_is_not(self):
        self.assertEqual(self.ignored("test/data/Annotations/new.bed.tar.gz"), [])


@unittest.skipUnless(WORKFLOW.is_file(), "needs the checkout's unit-tests.yml beside the script")
class APullRequestToWhatThisChecksRunsIt(unittest.TestCase):
    """unit-tests.yml runs on a pull request only for the paths it lists, and none of
    them covered the two files this module checks outside PostProcess/: a pull request
    that changed only install_from_source.sh or the .gitignore (dropping /Annotations/,
    say) did not run it before the merge."""

    def test_the_installer_and_the_gitignore_are_pull_request_paths(self):
        paths = pull_request_paths()
        self.assertIn("install_from_source.sh", paths)
        self.assertIn(".gitignore", paths)


if __name__ == "__main__":
    unittest.main(verbosity=2)
