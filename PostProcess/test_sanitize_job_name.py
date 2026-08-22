"""Regression test for the job-name sanitizer in pages/main_page.py.

A user-supplied job name becomes a results-directory component and is
interpolated into the shell search pipeline. A name like ``SBDS(T>C)`` used to
reach ``mkdir <path>`` run through ``/bin/sh`` and crash with
``Syntax error: "(" unexpected``. ``_sanitize_job_name`` restricts the name to a
filesystem/shell-safe charset. This test extracts the real function from source
(so it needs neither Dash nor a configured app) and pins its contract.
"""

import ast
import os
import re
import unittest

_MAIN_PAGE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pages", "main_page.py"
)


def _load_sanitizer():
    src = open(_MAIN_PAGE).read()
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "_sanitize_job_name":
            ns = {"re": re}
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
            return ns["_sanitize_job_name"]
    raise AssertionError("_sanitize_job_name not found in pages/main_page.py")


_SHELL_META = set(" \t()<>|&;$`'\"\\/*?![]{}#~")


class TestSanitizeJobName(unittest.TestCase):
    sanitize = staticmethod(_load_sanitizer())

    def test_reported_crash_name_is_safe(self):
        # the exact name from the bug report
        self.assertEqual(self.sanitize("SBDS(T>C)"), "SBDS_T_C")
        self.assertEqual(self.sanitize("DAN_SBDSP1/SBDS(T>C) "), "DAN_SBDSP1_SBDS_T_C")

    def test_no_shell_or_path_metacharacters_survive(self):
        for name in [
            "SBDS(T>C)",
            "a; rm -rf /",
            "guide & echo hi",
            "name|pipe",
            "sub/dir/name",
            "$(whoami)",
            "back`tick`",
            'quote"and\'quote',
            "star*and?glob",
        ]:
            out = self.sanitize(name)
            self.assertFalse(
                any(c in _SHELL_META for c in out),
                f"unsafe char left in {out!r} from {name!r}",
            )

    def test_normal_names_are_preserved(self):
        for name in ["normal_name-v2.1", "run42", "TRAC_NRG", "a.b_c-d"]:
            self.assertEqual(self.sanitize(name), name)

    def test_all_unsafe_reduces_to_empty(self):
        # caller falls back to the random-only job id when this is empty
        for name in ["((()))", "   ", "///", "><>"]:
            self.assertEqual(self.sanitize(name), "")

    def test_collapses_runs_and_trims(self):
        self.assertEqual(self.sanitize("  a   b  "), "a_b")
        self.assertEqual(self.sanitize("__lead_trail__"), "lead_trail")


if __name__ == "__main__":
    unittest.main()
