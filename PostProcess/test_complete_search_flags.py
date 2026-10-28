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

import os
import re
import subprocess
import sys
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

    Scoped to the lines above the first `def` on purpose: further down,
    `input_args` is also how OTHER subcommands read their own flags, and those
    are not valid for complete-search.
    """
    first_def = next((n for n, l in enumerate(lines) if l.startswith("def ")), len(lines))
    return {
        f
        for line in lines[:first_def]
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


@unittest.skipUnless(_importable(), "crisprme.py dependencies not installed")
class TestRejectionBehaviour(unittest.TestCase):
    """End to end: the gate runs before any other validation, so these need no data."""

    def _run(self, *flags):
        return subprocess.run(
            [sys.executable, CRISPRME, "complete-search", *flags],
            capture_output=True,
            text=True,
            cwd=REPO,
        )

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


if __name__ == "__main__":
    unittest.main(verbosity=2)
