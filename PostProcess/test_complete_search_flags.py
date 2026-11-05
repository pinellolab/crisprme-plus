"""Twin test for complete-search's flag allowlist (``crisprme.py``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_complete_search_flags.py' -v

or directly:

    <env>/bin/python3 PostProcess/test_complete_search_flags.py -v

complete-search parses ``sys.argv`` by hand rather than with argparse, so a flag
no code reads is simply ignored. That makes a typo silent AND expensive: with
``--prune-pams`` instead of ``--prune-pam`` the search runs genome-wide with
pruning off, finishes without complaint, and is only detectable by diffing the
results against a dedicated index. ``_check_unknown_flags`` exists to reject
those, and the allowlist it checks against is maintained by hand -- so the real
risk is drift: a flag documented in the help, or honoured by the code, that the
allowlist forgets, which would make a legitimate command fail.

The static checks below pin that invariant in both directions by parsing the
source (no import, so the suite needs none of CRISPRme's dependencies). The
subprocess checks confirm the behaviour end to end; they need the module to be
importable and are skipped otherwise.
"""

import ast
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CRISPRME = os.path.join(REPO, "crisprme.py")
FLAG = re.compile(r"--[A-Za-z0-9_-]+")
FLAG_LITERAL = re.compile(r'"(--[A-Za-z0-9_-]+)"')


def _source_lines():
    with open(CRISPRME, encoding="utf-8") as fh:
        return fh.read().split("\n")


def _block(lines, header):
    """Source lines of the top-level `def header...` block."""
    start = next(n for n, l in enumerate(lines) if l.startswith(header))
    end = next(
        (n for n, l in enumerate(lines[start + 1 :], start + 1) if l.startswith("def ")),
        len(lines),
    )
    return lines[start:end]


def documented_flags(lines):
    return set(FLAG.findall("\n".join(_block(lines, "def print_help_complete_search"))))


def allowlist(lines):
    src = "\n".join(lines)
    start = src.index("_COMPLETE_SEARCH_FLAGS = frozenset(")
    return set(FLAG_LITERAL.findall(src[start : src.index("\n)\n", start)]))


def global_flags(lines):
    """Flags read at MODULE scope off the full argv, so valid on every subcommand.

    Scoped to the module-level statements on purpose, a function's body left out:
    `input_args` is also how OTHER subcommands read their own flags, and those
    are not valid for complete-search. (The lines above the first `def` were
    taken for them until crisprme.py defined, above its `--ci-cd-test`, the
    functions it runs before its imports.)
    """
    src = "\n".join(lines)
    return {
        f
        for node in ast.parse(src).body
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        for line in lines[node.lineno - 1:node.end_lineno]
        if "input_args" in line
        for f in FLAG_LITERAL.findall(line)
    }


def honoured_flags(lines):
    """Flags complete-search actually reads: its own body, the helpers it calls,
    and the ones read at module scope from the full argv."""
    body = _block(lines, "def complete_search")
    text = "\n".join(body)
    found = set(FLAG_LITERAL.findall(text)) | global_flags(lines)
    for helper in set(
        re.findall(r"\b(_check_[a-z_]+|_registry_data_type_for_search)\s*\(", text)
    ):
        try:
            found |= set(FLAG_LITERAL.findall("\n".join(_block(lines, f"def {helper}"))))
        except StopIteration:
            pass
    return found


class TestAllowlistCoversEverything(unittest.TestCase):
    """The allowlist must not reject anything the help advertises or the code reads."""

    def setUp(self):
        self.lines = _source_lines()

    def test_every_documented_flag_is_accepted(self):
        missing = documented_flags(self.lines) - allowlist(self.lines)
        self.assertEqual(
            missing,
            set(),
            f"documented in --help but rejected by _check_unknown_flags: {sorted(missing)}",
        )

    def test_every_honoured_flag_is_accepted(self):
        missing = honoured_flags(self.lines) - allowlist(self.lines)
        self.assertEqual(
            missing,
            set(),
            f"read by complete-search but rejected by _check_unknown_flags: {sorted(missing)}",
        )

    def test_allowlist_has_no_dead_entries(self):
        """Every allowed flag is either documented or read, so the list cannot rot."""
        known = documented_flags(self.lines) | honoured_flags(self.lines)
        dead = allowlist(self.lines) - known
        self.assertEqual(
            dead, set(), f"allowed but neither documented nor read: {sorted(dead)}"
        )

    def test_prune_pam_is_documented_and_accepted(self):
        self.assertIn("--prune-pam", documented_flags(self.lines))
        self.assertIn("--prune-pam", allowlist(self.lines))


def _importable():
    return (
        subprocess.run(
            [sys.executable, "-c", "import Bio.Seq"], capture_output=True
        ).returncode
        == 0
    )


def scanned_files():
    """Yield every file the caller scan below reads: each .py, .sh, .yml and .yaml
    file in the checkout, wherever it is, but the tests and crisprme.py.

    test_ci_workflows checks that a pull request changing any of them starts the
    job that runs this scan.
    """
    for path in sorted(pathlib.Path(REPO).rglob("*")):
        if path.suffix not in (".py", ".sh", ".yml", ".yaml") or not path.is_file():
            continue
        if (path.name.startswith("test_") or ".git" in path.parts
                or path.resolve() == pathlib.Path(CRISPRME).resolve()):
            continue
        yield path


class TestCallersUseAcceptedFlags(unittest.TestCase):
    """Every complete-search flag the project itself passes must be accepted.

    The allowlist was checked against the help text and against the flags the
    parser reads, but not against the commands this repo builds. It shipped
    rejecting ``--bmax``, which ``complete_test.py`` had been passing for the
    brute-force benchmarks -- complete-search never had such a flag and silently
    ignored it, so nothing failed until the unknown-flag check arrived and broke
    the benchmark run. A static scan is cheap; finding this by running the
    benchmarks is not.
    """

    # flags that are not complete-search's own: a subcommand's or a shell's
    IGNORE = frozenset({"--help"})

    def _invocations(self):
        """Yield (path, flags) for every complete-search command this repo builds.

        Deliberately narrow, because a loose scan reports noise and gets muted:
        the trigger is a non-comment line literally building `crisprme.py
        complete-search`, and the command is then followed only through adjacent
        string pieces. crisprme.py itself is excluded -- it defines the flags, and
        its help text is covered by test_every_documented_flag_is_accepted.
        """
        for path in scanned_files():
            try:
                lines = path.read_text(encoding="utf-8").split("\n")
            except (OSError, UnicodeDecodeError):
                continue
            for i, line in enumerate(lines):
                if "crisprme.py complete-search" not in line:
                    continue
                if line.lstrip().startswith("#"):
                    continue
                block = [line]
                for nxt in lines[i + 1:i + 30]:
                    if nxt.lstrip().startswith("#"):
                        break
                    if '"' in nxt or "'" in nxt or block[-1].rstrip().endswith("\\"):
                        block.append(nxt)
                    else:
                        break
                flags = set(FLAG.findall(" ".join(block))) - self.IGNORE
                if flags:
                    yield path.relative_to(REPO), flags

    def test_scan_finds_the_benchmark_invocation(self):
        """Guard the scan itself: if it finds nothing it proves nothing."""
        paths = {str(p) for p, _ in self._invocations()}
        self.assertIn("PostProcess/complete_test.py", paths,
                      f"scan found only: {sorted(paths)}")

    def test_no_caller_passes_a_flag_the_allowlist_rejects(self):
        allowed = allowlist(_source_lines())
        offenders = []
        for path, flags in self._invocations():
            for flag in sorted(flags - allowed):
                offenders.append(f"{path}: {flag}")
        self.assertEqual(offenders, [], "complete-search flags no code accepts: "
                                        + ", ".join(offenders))

    def test_bmax_is_gone(self):
        """The specific regression: --bmax was never a complete-search flag."""
        self.assertNotIn("--bmax", allowlist(_source_lines()))
        for path, flags in self._invocations():
            self.assertNotIn("--bmax", flags, f"{path} still passes --bmax")


class TestDocumentedCommandsPassTheFlagChecks(unittest.TestCase):
    """Every complete-search command the docs give must get past complete-search's own
    checks of its flags, or a user who copies it gets an error and nothing else.

    The general-index Hugging Face card's only example ran ``--prune-pam`` without
    ``--index-path``, which complete-search refuses up front ("--prune-pam requires
    --index-path"). The scan reads the fenced code blocks of the README, METHODS and
    docs/ (the Hugging Face cards included) and follows backslash continuations.
    """

    IGNORE = frozenset({"--help"})

    def _commands(self):
        """Yield (path, flags) for every `crisprme.py complete-search` in a code block."""
        paths = [pathlib.Path(REPO, "README.md"), pathlib.Path(REPO, "METHODS.md")]
        paths += sorted(pathlib.Path(REPO, "docs").rglob("*.md"))
        for path in paths:
            lines = path.read_text(encoding="utf-8").split("\n")
            fenced = False
            for i, line in enumerate(lines):
                if line.lstrip().startswith("```"):
                    fenced = not fenced
                    continue
                if (not fenced or "crisprme.py complete-search" not in line
                        or line.lstrip().startswith("#")):
                    continue
                cmd = [line.split(" #", 1)[0]]
                while cmd[-1].rstrip().endswith("\\") and i + len(cmd) < len(lines):
                    cmd.append(lines[i + len(cmd)].split(" #", 1)[0])
                args = " ".join(cmd).split("crisprme.py complete-search", 1)[1]
                yield path.relative_to(REPO), set(FLAG.findall(args)) - self.IGNORE

    def test_scan_finds_the_general_index_example(self):
        """Guard the scan itself: if it finds nothing it proves nothing."""
        paths = {str(p) for p, _ in self._commands()}
        self.assertIn("docs/hf/crisprme-data-general-index_README.md", paths,
                      f"scan found only: {sorted(paths)}")

    def test_no_documented_command_passes_a_flag_the_allowlist_rejects(self):
        allowed = allowlist(_source_lines())
        offenders = [f"{path}: {flag}" for path, flags in self._commands()
                     for flag in sorted(flags - allowed)]
        self.assertEqual(offenders, [], "documented complete-search flags no code "
                                        "accepts: " + ", ".join(offenders))

    def test_prune_pam_is_always_given_with_index_path(self):
        offenders = sorted({str(path) for path, flags in self._commands()
                            if "--prune-pam" in flags and "--index-path" not in flags})
        self.assertEqual(offenders, [], "--prune-pam without --index-path, which "
                                        "complete-search refuses: " + ", ".join(offenders))


class TestNoDocumentedCommandIsABareCrisprme(unittest.TestCase):
    """A bare ``crisprme.py`` prints its list of commands and exits 1, and first
    creates the CRISPRme folders (``Genomes/``, ``Results/`` and five more) in the
    directory it runs in. README §1.1.2 gave it to check an install, and §3.1 to
    show the help; ``--version`` and ``--help`` create nothing."""

    def test_no_code_block_runs_crisprme_py_alone(self):
        paths = [pathlib.Path(REPO, "README.md"), pathlib.Path(REPO, "METHODS.md")]
        paths += sorted(pathlib.Path(REPO, "docs").rglob("*.md"))
        offenders = []
        for path in paths:
            fenced = False
            for n, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
                if line.lstrip().startswith("```"):
                    fenced = not fenced
                elif fenced and line.split(" #", 1)[0].strip().lstrip("$ ") == "crisprme.py":
                    offenders.append(f"{path.relative_to(REPO)}:{n}")
        self.assertEqual(offenders, [], "a bare crisprme.py, which exits 1 and creates "
                                        "seven folders: " + ", ".join(offenders))


class TestThePrunePamHelpSaysWhichSideItServes(unittest.TestCase):
    """The published PAM-less indexes store their flank on the 3' side, and a 5'
    search of one found nothing, exit 0. The --prune-pam help said only that a 5'
    PAM "must equal" the flank width, and that the result is identical to a
    dedicated index's."""

    def _entry(self):
        block = "\n".join(_block(_source_lines(), "def print_help_complete_search"))
        text = "".join(re.findall(r'^\s*"(.*)"\s*$', block, flags=re.M))
        start = text.index("\\t--prune-pam,")
        return text[start:text.index("\\n", start)]

    def test_it_says_a_5_prime_pam_needs_an_index_built_for_one(self):
        entry = self._entry()
        self.assertIn("A 5' PAM (Cas12a) needs an index built for a 5' PAM", entry)
        self.assertIn("published PAM-less indexes store the flank on the 3' side", entry)
        self.assertNotIn("a 5' PAM must equal it", entry)

    def test_it_says_what_a_flank_wider_than_the_pam_cannot_index(self):
        # (it said "The result is identical to a dedicated index built for that PAM", which
        # does not hold for a PAM at a contig's very end or start: idx audit IDX-4)
        entry = self._entry()
        self.assertIn("From a flank as wide as the PAM the result is identical to a dedicated "
                      "index built for that PAM", entry)
        self.assertIn("are not indexed (the search says so)", entry)
        self.assertIn("built by CRISPRme 2.8.1 or later", entry)
        self.assertNotIn("The result is identical to a dedicated index", entry)


class TestTheComputeBackendHelpSaysWhatAMissingScorerEnvDoes(unittest.TestCase):
    """The --compute-backend help said "if the cbulge scorer env is absent the run completes
    on CFD only". complete-search refuses to start unless CRISPRME_SKIP_SCORER_ENV=1 (since
    2.7.3), so a user who read the help got exit 1. The help now says so; the refusal itself
    is pinned by test_annotation_tmpdirs, which runs complete-search with an env not ready."""

    def _entry(self):
        block = "\n".join(_block(_source_lines(), "def print_help_complete_search"))
        text = "".join(re.findall(r'^\s*"(.*)"\s*$', block, flags=re.M))
        start = text.index("\\t--compute-backend,")
        return text[start:text.index("\\n", start)]

    def test_it_says_complete_search_refuses_and_names_the_fix_and_the_opt_out(self):
        entry = self._entry()
        self.assertIn("complete-search refuses to start while the cbulge scorer env or its "
                      "weights are missing or broken", entry)
        self.assertIn("`crisprme.py scorer-env create`", entry)
        self.assertIn("`crisprme.py scorer-env doctor`", entry)
        self.assertIn("set CRISPRME_SKIP_SCORER_ENV=1 to run CFD-only", entry)
        self.assertNotIn("completes on CFD only", entry)

    def test_the_opt_out_it_names_is_the_one_complete_search_reads(self):
        body = "\n".join(_block(_source_lines(), "def complete_search"))
        self.assertIn('os.environ.get("CRISPRME_SKIP_SCORER_ENV", "0") == "1"', body)


# Needs crisprme.py to be importable: these run it as a subprocess, and it imports
# Bio.Seq at module scope. The static checks above parse the source instead and must
# stay unguarded -- gating them on a dependency they do not use is how the --bmax
# regression guard silently stopped running in CI.
@unittest.skipUnless(_importable(), "crisprme.py dependencies not installed")
class TestRejectionBehaviour(unittest.TestCase):
    """End to end: the gate runs before any other validation, so these need no data."""

    def _run(self, *flags):
        # From a folder of its own: past the gate crisprme.py makes its data folders in the
        # current one, which on a read-only code directory (a SIF) is a traceback before
        # "--genome is required", and on a writable one leaves them beside the code.
        with tempfile.TemporaryDirectory() as tmp:
            return subprocess.run(
                [sys.executable, CRISPRME, "complete-search", *flags],
                capture_output=True,
                text=True,
                cwd=tmp,
            )

    def test_it_leaves_no_data_folder_beside_the_code(self):
        """What a run of this class from the code directory did to an install, and what
        stopped it on a read-only one."""
        before = set(os.listdir(REPO))
        self._run("--prune-pam")
        self.assertEqual(set(os.listdir(REPO)) - before, set())

    def test_typo_is_rejected_with_a_suggestion(self):
        r = self._run("--prune-pams")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("unrecognized complete-search flag", r.stderr)
        self.assertIn("did you mean --prune-pam?", r.stderr)

    def test_real_flag_passes_the_gate(self):
        """--prune-pam must get past the flag gate and fail later, on missing inputs."""
        r = self._run("--prune-pam")
        self.assertNotEqual(r.returncode, 0)
        self.assertNotIn("unrecognized complete-search flag", r.stderr)
        self.assertIn("--genome is required", r.stderr)

    def test_retired_flag_still_warns_rather_than_failing_as_unknown(self):
        r = self._run("--full")
        self.assertNotIn("unrecognized complete-search flag", r.stderr)

    def test_unknown_flag_with_no_close_match_is_still_rejected(self):
        r = self._run("--zzzzzzzz")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("unrecognized complete-search flag", r.stderr)
        self.assertIn("--zzzzzzzz", r.stderr)


@unittest.skipUnless(_importable(), "crisprme.py dependencies not installed")
class TestAnInputErrorIsOneLineBeforeAnyWork(unittest.TestCase):
    """A PAM file named against the <length>-<motif>-<CasName>.txt convention, one holding
    two motif lines or a first line that is not '<motif> <length>', and --index-path or
    --vcf-filter-pass-values with no value ended complete-search in a Python traceback (a
    raw ValueError), all but the two motif lines' only once Results/<output> was made, which
    stayed, empty, and the misnamed file only once the annotation was sorted (~30 s with the
    shipped one). Each is one 'Error: ...' line and exit 1 now, as the other input errors,
    checked before the annotation is read or anything is made."""

    GUIDE = "CTAACAGTTGCTTTTATCACNNN"

    def setUp(self):
        self.wd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.wd, ignore_errors=True)
        os.makedirs(os.path.join(self.wd, "Genomes", "g"))
        with open(os.path.join(self.wd, "Genomes", "g", "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 30 + "\n")
        with open(os.path.join(self.wd, "guide.txt"), "w") as fh:
            fh.write(self.GUIDE + "\n")
        os.makedirs(os.path.join(self.wd, "PAMs"))

    def pam(self, name, text="NNNNNNNNNNNNNNNNNNNNNGG 3\n"):
        with open(os.path.join(self.wd, "PAMs", name), "w") as fh:
            fh.write(text)
        return os.path.join("PAMs", name)

    def search(self, pam, *flags):
        """complete-search of the guide on Genomes/g into Results/out (no scorer env needed:
        what is checked comes before it)."""
        env = dict(os.environ, CRISPRME_SKIP_SCORER_ENV="1")
        return subprocess.run(
            [sys.executable, CRISPRME, "complete-search", "--genome", "Genomes/g", "--guide",
             "guide.txt", "--pam", pam, "--mm", "4", "--output", "out", *flags],
            capture_output=True, text=True, cwd=self.wd, env=env)

    def assertRefused(self, run, says):
        """Exit 1, with ``says`` the last line on stderr (a low-memory host's warning may
        come first), no traceback, and no Results/out."""
        self.assertEqual(run.returncode, 1, run.stderr)
        self.assertNotIn("Traceback", run.stderr)
        self.assertEqual(run.stderr.splitlines()[-1:], [f"Error: {says}"], run.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.wd, "Results", "out")))

    def test_a_pam_file_named_against_the_convention(self):
        self.assertRefused(
            self.search(self.pam("mypam.txt")),
            "Invalid PAM filename 'mypam.txt': expected the '<length>-<motif>-<CasName>.txt' "
            "convention (e.g. 20bp-NGG-SpCas9.txt).")

    def test_before_the_annotation_is_read(self):
        # the check of an annotation that is not there would refuse it first
        self.assertRefused(
            self.search(self.pam("mypam.txt"), "--annotation", "no_such_annotation.bed"),
            "Invalid PAM filename 'mypam.txt': expected the '<length>-<motif>-<CasName>.txt' "
            "convention (e.g. 20bp-NGG-SpCas9.txt).")

    def test_a_pam_file_of_two_motif_lines(self):
        pam = self.pam("20bp-NRG-SpCas9.txt",
                       "NNNNNNNNNNNNNNNNNNNNNGG 3\nNNNNNNNNNNNNNNNNNNNNNAG 3\n")
        self.assertRefused(
            self.search(pam),
            f"Only one PAM motif per file is supported; the PAM file "
            f"'{os.path.join(os.path.realpath(self.wd), pam)}' contains more than one "
            f"non-empty line. Combine the motifs into a single IUPAC motif "
            f"(e.g. NAG + NGG -> NRG).")

    def test_a_pam_file_whose_line_is_not_a_motif_and_a_length(self):
        pam = self.pam("20bp-NGG-SpCas9.txt", "NNNNNNNNNNNNNNNNNNNNNGG\n")
        self.assertRefused(
            self.search(pam),
            f"The PAM file '{os.path.join(os.path.realpath(self.wd), pam)}' must hold the "
            f"motif and the PAM's length on its first line, as in "
            f"'NNNNNNNNNNNNNNNNNNNNNGG 3' (a negative length for a PAM at the 5' end of the "
            f"guide).")

    def test_a_flag_with_no_value(self):
        pam = self.pam("20bp-NGG-SpCas9.txt")
        for flag in ("--index-path", "--vcf-filter-pass-values"):
            with self.subTest(flag=flag, given="last"):
                self.assertRefused(self.search(pam, flag), f"Missing input for {flag}")
            with self.subTest(flag=flag, given="before another flag"):
                self.assertRefused(self.search(pam, flag, "--thread", "2"),
                                   f"Please input a value for flag {flag}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
