"""Unit tests for the 2.8.1 web Settings download-progress + failure-relay
helpers in pages/settings_page.py.

The Settings page imports a live Dash ``app`` (and dash/dbc), which the
network-free CI does NOT install, so — like test_history_page_assembly.py /
test_sanitize_job_name.py — the pure helper functions are EXTRACTED from source
via ast and exec'd with a tiny namespace (stdlib + light stubs for html/dbc).
This keeps the progress-JSON parsing, ETA formatting, full-log assembly, the
genotypes-degraded detector, and the retry-command round-trip covered in CI
without a Dash runtime.
"""

import ast
import json
import os
import tempfile
import time
import unittest
from typing import List, Optional

_SETTINGS_PAGE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "pages", "settings_page.py",
)

_NAMES = [
    "_read_progress",
    "_fmt_eta",
    "_render_progress_bar",
    "_full_job_log",
    "_genotypes_degraded_warning",
    "_fmt_size",
]


class _Stub:
    """A minimal stand-in for an html/dbc element that records its args so tests
    can assert on what was rendered without a real Dash runtime."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __repr__(self):
        return f"_Stub({self.args!r}, {self.kwargs!r})"


class _StubModule:
    """Any attribute access returns the _Stub factory (html.Div, dbc.Progress...)."""

    def __getattr__(self, _name):
        return _Stub


def _load_functions():
    src = open(_SETTINGS_PAGE).read()
    tree = ast.parse(src)
    ns = {
        "os": os,
        "time": time,
        "json": json,
        "List": List,
        "Optional": Optional,
        "html": _StubModule(),
        "dbc": _StubModule(),
        # module-level constants the helpers reference
        "SETTINGS_DIR": "Settings",
        "_PROGRESS_STALE_SECONDS": 120,
        "_BENIGN_ERR_PATTERNS": (),
        "current_working_directory": "",
    }
    found = set()
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in _NAMES:
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
            found.add(node.name)
    missing = set(_NAMES) - found
    if missing:
        raise AssertionError(f"Not found in pages/settings_page.py: {sorted(missing)}")
    return ns


_NS = _load_functions()


class TestReadProgress(unittest.TestCase):
    read_progress = staticmethod(_NS["_read_progress"])

    def _jobdir(self):
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(d, ignore_errors=True))
        return d

    def _write(self, jobdir, payload):
        with open(os.path.join(jobdir, "progress.json"), "w") as fh:
            json.dump(payload, fh)

    def test_absent_is_none(self):
        self.assertIsNone(self.read_progress(self._jobdir()))

    def test_fresh_payload_parsed(self):
        d = self._jobdir()
        self._write(d, {"phase": "Downloading index", "bytes_done": 5,
                        "bytes_total": 10, "pct": 50.0, "updated": time.time()})
        got = self.read_progress(d)
        self.assertIsNotNone(got)
        self.assertEqual(got["pct"], 50.0)

    def test_stale_is_none(self):
        d = self._jobdir()
        self._write(d, {"pct": 10.0, "updated": time.time() - 10_000})
        self.assertIsNone(self.read_progress(d))  # frozen -> fall back to spinner

    def test_malformed_is_none(self):
        d = self._jobdir()
        with open(os.path.join(d, "progress.json"), "w") as fh:
            fh.write("{not json")
        self.assertIsNone(self.read_progress(d))

    def test_non_dict_is_none(self):
        d = self._jobdir()
        with open(os.path.join(d, "progress.json"), "w") as fh:
            json.dump([1, 2, 3], fh)
        self.assertIsNone(self.read_progress(d))


class TestFmtEta(unittest.TestCase):
    fmt_eta = staticmethod(_NS["_fmt_eta"])

    def test_none_and_zero(self):
        self.assertEqual(self.fmt_eta(None), "")
        self.assertEqual(self.fmt_eta(0), "")
        self.assertEqual(self.fmt_eta(-5), "")

    def test_seconds(self):
        self.assertEqual(self.fmt_eta(30), "~30 sec left")

    def test_minutes(self):
        self.assertEqual(self.fmt_eta(180), "~3 min left")


class TestRenderProgressBar(unittest.TestCase):
    render = staticmethod(_NS["_render_progress_bar"])

    def test_none_when_no_bytes_no_total(self):
        self.assertIsNone(self.render("stage", {"bytes_done": 0, "bytes_total": 0,
                                                 "pct": 0}))

    def test_renders_with_total(self):
        out = self.render("Downloading index",
                          {"bytes_done": 500, "bytes_total": 1000, "pct": 50.0,
                           "phase": "Downloading index 1/3", "detail": "NRG",
                           "eta_seconds": 120})
        self.assertIsNotNone(out)  # a stubbed html.Div

    def test_renders_with_unknown_total(self):
        # size API offline: no total but bytes are moving -> still render a bar
        out = self.render("dl", {"bytes_done": 123, "bytes_total": 0, "pct": 0})
        self.assertIsNotNone(out)


class TestFullJobLog(unittest.TestCase):
    full_log = staticmethod(_NS["_full_job_log"])

    def _jobdir(self):
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(d, ignore_errors=True))
        return d

    def test_no_logs(self):
        self.assertIn("No log files", self.full_log(self._jobdir()))

    def test_includes_error_and_verbose(self):
        d = self._jobdir()
        with open(os.path.join(d, "log_error.txt"), "w") as fh:
            fh.write("Traceback: boom\n")
        with open(os.path.join(d, "log_verbose.txt"), "w") as fh:
            fh.write("line1\nline2\n")
        out = self.full_log(d)
        self.assertIn("log_error.txt", out)
        self.assertIn("boom", out)
        self.assertIn("log_verbose.txt", out)
        self.assertIn("line2", out)

    def test_verbose_tail_truncation(self):
        d = self._jobdir()
        with open(os.path.join(d, "log_verbose.txt"), "w") as fh:
            for i in range(1000):
                fh.write(f"row{i}\n")
        out = self.full_log(d, verbose_tail=50)
        self.assertIn("showing last 50 lines", out)
        self.assertIn("row999", out)
        self.assertNotIn("row0\n", out)  # head dropped


class TestGenotypesDegradedWarning(unittest.TestCase):
    warn = staticmethod(_NS["_genotypes_degraded_warning"])

    def _jobdir(self):
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(d, ignore_errors=True))
        return d

    def test_none_without_verbose(self):
        self.assertIsNone(self.warn(self._jobdir()))

    def test_none_when_installed(self):
        d = self._jobdir()
        with open(os.path.join(d, "log_verbose.txt"), "w") as fh:
            fh.write("Installed Tier-1 genotype store genotypes_hg38_X/\n")
        self.assertIsNone(self.warn(d))

    def test_warns_on_degraded(self):
        d = self._jobdir()
        with open(os.path.join(d, "log_verbose.txt"), "w") as fh:
            fh.write("NOTE: could not fetch Tier-1 genotype store "
                     "'genotypes_hg38_X.tar.gz' (timeout); continuing without "
                     "it — off-target detection works, Samples degraded.\n")
        self.assertIsNotNone(self.warn(d))  # a stubbed dbc.Alert

    def test_no_false_positive_on_plain_success(self):
        d = self._jobdir()
        with open(os.path.join(d, "log_verbose.txt"), "w") as fh:
            fh.write("Downloaded index NRG_3_hg38 (built 2026-01-01)\n")
        self.assertIsNone(self.warn(d))


if __name__ == "__main__":
    unittest.main()
