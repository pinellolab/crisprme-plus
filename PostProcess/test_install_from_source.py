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
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SCRIPT = Path(__file__).resolve().parent.parent / "install_from_source.sh"
REPO = SCRIPT.parent
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
        """Which of `paths` git ignores, created in a repository with that .gitignore.

        git runs with a home of its own, no system configuration and none of the
        caller's GIT_* variables, so that only that .gitignore decides: a global
        excludes file of the developer's (core.excludesFile, or git's default
        ~/.config/git/ignore) holding, say, *.tar.gz would change the answer."""
        with tempfile.TemporaryDirectory() as d:
            env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
            env.update(HOME=d, XDG_CONFIG_HOME=d, GIT_CONFIG_NOSYSTEM="1")
            repo = Path(d) / "checkout"
            subprocess.run(["git", "init", "-q", str(repo)], check=True, env=env)
            shutil.copy(SCRIPT.parent / ".gitignore", repo)
            for path in paths:
                (repo / path).parent.mkdir(parents=True, exist_ok=True)
                (repo / path).write_text("")
            out = subprocess.run(["git", "-C", str(repo), "check-ignore", *paths],
                                 capture_output=True, text=True, env=env)
        # 1 is "none of them": anything else but 0 is an error, not an answer
        self.assertIn(out.returncode, (0, 1), out.stderr)
        return out.stdout.split()

    def test_the_scorer_env_state_is_ignored(self):
        self.assertEqual(self.ignored("Annotations/.scorer_env.json"),
                         ["Annotations/.scorer_env.json"])

    def test_a_new_file_in_any_data_folder_a_search_creates_is_ignored(self):
        new = [f"{d}/new.txt" for d in search_data_folders()]
        self.assertEqual(sorted(self.ignored(*new)), sorted(new))

    def test_a_new_file_beside_the_tracked_test_annotations_is_not(self):
        self.assertEqual(self.ignored("test/data/Annotations/new.bed.tar.gz"), [])

    def test_the_answers_do_not_depend_on_the_developers_own_git_configuration(self):
        # git configured, the ways a developer's machine can be, to ignore every *.tar.gz:
        # core.excludesFile in ~/.gitconfig, git's default global excludes file, and the
        # system configuration
        with tempfile.TemporaryDirectory() as home:
            excludes = Path(home, "excludes")
            excludes.write_text("*.tar.gz\n")
            config = f"[core]\n\texcludesFile = {excludes}\n"
            Path(home, ".gitconfig").write_text(config)
            Path(home, "gitconfig_system").write_text(config)
            Path(home, ".config", "git").mkdir(parents=True)
            Path(home, ".config", "git", "ignore").write_text("*.tar.gz\n")
            with mock.patch.dict(os.environ, {
                    "HOME": home, "XDG_CONFIG_HOME": str(Path(home, ".config")),
                    "GIT_CONFIG_GLOBAL": str(Path(home, ".gitconfig")),
                    "GIT_CONFIG_SYSTEM": str(Path(home, "gitconfig_system"))}):
                self.assertEqual(self.ignored("test/data/Annotations/new.bed.tar.gz"), [])
                self.assertEqual(self.ignored("Annotations/.scorer_env.json"),
                                 ["Annotations/.scorer_env.json"])


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


class TheDocsNameTheCRISPRitzThisScriptBuilds(unittest.TestCase):
    """A bare `install_from_source.sh` in the docs is this script, which builds the
    CRISPRitz tag CRISPRITZ_REF defaults to. The 2.8.1 notes and INDEX_PROVENANCE said
    2.9.0 is what "the 2.8.0 image and `install_from_source.sh`" build: 2.8.0's script
    does, this one builds v2.9.1. The scan joins the wrapped lines of README, METHODS,
    docs/ and the CHANGELOG down to its newest release into sentences. One in which the
    script builds something may name another CRISPRitz of the same line only when it
    says whose script that is ("its", "2.8.0's")."""

    BUILDS = re.compile(r"`install_from_source\.sh`(?: script)? (?:build|builds|built|compiles)\b")
    OWNED = re.compile(r"(?:\bits|\b\d+\.\d+\.\d+['’]s|\bthe \d+\.\d+\.\d+)"
                       r" `install_from_source\.sh`")

    @staticmethod
    def _sentences():
        """Yield (path, sentence) for every sentence that names the script."""
        paths = [REPO / "README.md", REPO / "METHODS.md", *sorted((REPO / "docs").rglob("*.md"))]
        texts = [(p, p.read_text(encoding="utf-8")) for p in paths]
        log = (REPO / "CHANGELOG.md").read_text(encoding="utf-8")
        releases = [m.start() for m in re.finditer(r"^## \[\d", log, flags=re.M)]
        texts.append((REPO / "CHANGELOG.md", log[:releases[1]]))
        for path, text in texts:
            for block in re.split(r"\n\s*\n|\n(?=\s*(?:[-*|#]|\d+\.)\s)", text):
                for sentence in re.split(r"(?<=[.;!?])\s+(?=[A-Z`*(])", " ".join(block.split())):
                    if "install_from_source.sh" in sentence:
                        yield str(path.relative_to(REPO)), sentence

    def test_scan_finds_the_sentence_on_what_2_8_0_builds(self):
        """Guard the scan itself: if it finds nothing it proves nothing."""
        paths = {path for path, sentence in self._sentences() if self.BUILDS.search(sentence)}
        self.assertIn("docs/INDEX_PROVENANCE.md", paths, f"scan found only: {sorted(paths)}")

    def test_a_sentence_on_what_the_script_builds_names_its_default(self):
        default = re.search(r'^CRISPRITZ_REF="\$\{CRISPRITZ_REF:-v(\d+\.\d+\.\d+)\}"',
                            SCRIPT.read_text(), flags=re.M).group(1)
        same_line = re.compile(rf"\bv?({re.escape(default.rsplit('.', 1)[0])}\.\d+)\b")
        offenders = []
        for path, sentence in self._sentences():
            bare = self.OWNED.sub("", sentence)
            if self.BUILDS.search(bare) and set(same_line.findall(bare)) - {default}:
                offenders.append(f"{path}: {sentence[:80]}")
        self.assertEqual(offenders, [], f"this script builds CRISPRitz {default}; a sentence "
                                        "that has it build another must say whose script "
                                        "that is: " + "; ".join(offenders))


class TheVersionSweepComparesLikeWithLike(unittest.TestCase):
    """docs/RELEASE_CHECKLIST.md's version sweep listed the Dockerfile's CRISPRitz ARG among
    the values that must read CRISPRme's own X.Y.Z: v2.9.1 against 2.8.1, a rule no correct
    release could pass, and one a maintainer could 'fix' by building CRISPRitz v2.8.1, a tag
    that exists. The ARG is the CRISPRitz pin: the tag this script builds by default."""

    def test_the_dockerfile_builds_the_crispritz_this_script_builds(self):
        default = re.search(r'^CRISPRITZ_REF="\$\{CRISPRITZ_REF:-(v\d+\.\d+\.\d+)\}"',
                            SCRIPT.read_text(), flags=re.M).group(1)
        ref = re.search(r"^ARG crispritz_ref=(v\d+\.\d+\.\d+)\s*$",
                        (REPO / "Dockerfile").read_text(), flags=re.M).group(1)
        self.assertEqual(ref, default)

    def test_the_sweep_asks_x_y_z_of_crisprme_and_the_crispritz_pin_of_the_dockerfile(self):
        sweep = next(ln for ln in (REPO / "docs" / "RELEASE_CHECKLIST.md").read_text(
            encoding="utf-8").splitlines() if ln.startswith("- [ ] **Version sweep**"))
        self.assertNotIn("`Dockerfile` (CRISPRitz ARG)", sweep)
        # it pins none, but its header comment names the release (the next test)
        self.assertNotIn("carries no CRISPRme version", sweep)
        self.assertIn("The `Dockerfile` pins no CRISPRme version (its header comment names the "
                      "release: bump it with the rest", sweep)
        self.assertIn("its `ARG crispritz_ref` is the CRISPRitz pin", sweep)
        self.assertIn("must equal `install_from_source.sh`'s `CRISPRITZ_REF` default", sweep)

    def test_the_dockerfile_header_names_the_release_crisprme_py_is(self):
        """The Dockerfile pins no CRISPRme version, but its header comment names the release,
        bumped at every release since 2.7.9 (it named 2.4.0 at every tag from v2.5.0 to
        v2.7.8). scripts/prepare_release.py bumps crisprme.py alone, so the release docs bump
        the comment by hand, and a release that forgets it fails here."""
        version = re.search(r'^version = "(\d+\.\d+\.\d+)"',
                            (REPO / "crisprme.py").read_text(encoding="utf-8"), flags=re.M)
        header = re.search(r"^# CRISPRme\+ (\d+\.\d+\.\d+) — ",
                           (REPO / "Dockerfile").read_text(encoding="utf-8"), flags=re.M)
        self.assertIsNotNone(version, "crisprme.py has no version line")
        self.assertIsNotNone(header, "the Dockerfile's header comment names no release")
        self.assertEqual(header.group(1), version.group(1),
                         "the Dockerfile's header comment names another release than crisprme.py")

    def test_the_release_docs_bump_the_header_and_nothing_else_of_the_dockerfile(self):
        """RELEASING.md and the release skill said the Dockerfile carries no CRISPRme version
        (the skill: 'Dockerfile -> nothing', and a bump whose diff 'must show ONLY
        crisprme.py'), so a release that followed them left the header comment behind. Each now
        gives the command that bumps it: run here, with X.Y.Z filled in, on a copy of the
        Dockerfile, it changes that one line and leaves no file behind."""
        dockerfile = (REPO / "Dockerfile").read_text(encoding="utf-8")
        for doc in (REPO / "docs" / "RELEASING.md",
                    REPO / ".claude" / "skills" / "release-crisprme" / "SKILL.md"):
            text = doc.read_text(encoding="utf-8")
            with self.subTest(doc=doc.name):
                for wrong in ("carries no CRISPRme version", "carried no CRISPRme version",
                              "`Dockerfile` → nothing", "must show ONLY crisprme.py\n"):
                    self.assertFalse(wrong in text, f"{doc.name} still says {wrong!r}")
                bump = [ln for ln in text.splitlines()
                        if ln.startswith("sed ") and "# CRISPRme+ X.Y.Z" in ln]
                self.assertEqual(len(bump), 1, "no one command bumps the Dockerfile's header")
                with tempfile.TemporaryDirectory() as tmp:
                    Path(tmp, "Dockerfile").write_text(dockerfile, encoding="utf-8")
                    subprocess.run(["bash", "-c", bump[0].replace("X.Y.Z", "9.10.11")],
                                   cwd=tmp, check=True)
                    self.assertEqual(os.listdir(tmp), ["Dockerfile"])
                    after = Path(tmp, "Dockerfile").read_text(encoding="utf-8")
                changed = [(a, b) for a, b in zip(dockerfile.split("\n"), after.split("\n"))
                           if a != b]
                self.assertEqual(len(after.split("\n")), len(dockerfile.split("\n")))
                self.assertEqual(len(changed), 1, changed)
                self.assertRegex(changed[0][0], r"^# CRISPRme\+ \d+\.\d+\.\d+ — ")
                self.assertEqual(changed[0][1],
                                 re.sub(r"\d+\.\d+\.\d+", "9.10.11", changed[0][0], count=1))


if __name__ == "__main__":
    unittest.main(verbosity=2)
