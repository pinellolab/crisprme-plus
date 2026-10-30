"""Unit tests for the 2.8.1 web Settings download-progress + failure-relay
helpers in pages/settings_page.py.

The Settings page imports a live Dash ``app`` (and dash/dbc), which the
network-free CI does NOT install, so — like test_history_page_assembly.py /
test_sanitize_job_name.py — the pure helper functions are EXTRACTED from source
via ast and exec'd with a tiny namespace (stdlib + light stubs for html/dbc).
This keeps the progress-JSON parsing, ETA formatting, full-log assembly, the
genotypes-degraded detector, and the retry-command round-trip covered in CI
without a Dash runtime. The job-id guard (_settings_jobdir) is covered the same
way, through the undecorated Retry / Show-full-log / progress callbacks.
"""

import ast
import json
import os
import random
import re
import shutil
import string
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


def _load_functions(names=_NAMES, **overrides):
    """Exec the named functions (and module-level assignments) from source into
    a fresh namespace; ``overrides`` adds or replaces names in it (stubs)."""
    with open(_SETTINGS_PAGE) as fh:
        src = fh.read()
    tree = ast.parse(src)
    ns = {
        "os": os,
        "time": time,
        "json": json,
        "re": re,
        "random": random,
        "string": string,
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
    ns.update(overrides)
    found = set()
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
            found.add(node.name)
        elif isinstance(node, ast.Assign):
            targets = {t.id for t in node.targets if isinstance(t, ast.Name)}
            if targets & set(names):
                exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
                found |= targets
    missing = set(names) - found
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


# --- job-id guard --------------------------------------------------------------
# A Settings job id comes back from the browser as dcc.Store data, which a client
# can set to anything. It must never become a path outside Settings/: Retry runs
# <jobdir>/cmd.json with shell=True, and the log readers return file contents.
_GOOD_ID = "set_ABCDEFGH12"  # the shape _new_job_id makes
_LINK_ID = "set_LINKOUT000"  # job-shaped, but a symlink out of Settings/


class _PreventUpdate(Exception):
    """Stand-in for dash.exceptions.PreventUpdate."""


_NO_UPDATE = object()  # stand-in for dash.no_update


class _JobDirCase(unittest.TestCase):
    """A temp data dir with Settings/ and, beside it, an ``outside`` dir holding
    a job-shaped cmd.json + logs that no callback may ever read or run."""

    names = ["_JOB_ID_RE", "_settings_jobdir"]
    stubs = {}

    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)
        self.settings = os.path.join(self.cwd, "Settings")
        os.makedirs(os.path.join(self.settings, "uploads"))
        self.outside = os.path.join(self.cwd, "outside")
        os.makedirs(self.outside)
        with open(os.path.join(self.outside, "cmd.json"), "w") as fh:
            json.dump({"cmd": "echo SHOULD-NOT-RUN", "stage": "Outside"}, fh)
        for fname in ("log_error.txt", "log_verbose.txt"):
            with open(os.path.join(self.outside, fname), "w") as fh:
                fh.write("OUTSIDE-SECRET\n")
        with open(os.path.join(self.outside, "log.txt"), "w") as fh:
            fh.write("Outside\tQueued\nOutside\tStart\nOutside\tFAILED\n")
        os.symlink(self.outside, os.path.join(self.settings, _LINK_ID))
        self.launched = []  # (cmd, stage) pairs handed to _launch_job_cmd
        self.opened = []  # real path of every file the extracted code opens

        def _launch(cmd, stage):  # records instead of running anything
            self.launched.append((cmd, stage))
            return "set_NEWJOB0001"

        def _open(path, *args, **kwargs):
            self.opened.append(os.path.realpath(path))
            return open(path, *args, **kwargs)

        self.ns = _load_functions(
            self.names,
            current_working_directory=self.cwd,
            _launch_job_cmd=_launch,
            open=_open,
            ONLINE=False,
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
            **self.stubs,
        )

    def bad_ids(self):
        return [
            "../outside",
            "../../outside",
            _GOOD_ID + "/../../outside",
            "set_../../outside",
            self.outside,  # absolute path
            "/" + _GOOD_ID,  # absolute path with a job-shaped tail
            _LINK_ID,  # symlink that escapes Settings/
            "..", ".", "", "uploads",  # Settings/'s parent, itself, its uploads dir
            "set_abcdefgh12",  # wrong alphabet
            "set_ABCDEFGH1", "set_ABCDEFGH123",  # wrong length
            _GOOD_ID + "\n", _GOOD_ID + "/", " " + _GOOD_ID,
            "pub_ABCDEFGH12",  # a prefix _new_job_id is never called with
            None, 5, [_GOOD_ID], {"id": _GOOD_ID},  # not a string
        ]

    def outside_reads(self):
        root = os.path.realpath(self.outside) + os.sep
        return [p for p in self.opened if p.startswith(root)]

    def make_job(self, job_id, files):
        jobdir = os.path.join(self.settings, job_id)
        os.makedirs(jobdir)
        for fname, text in files.items():
            with open(os.path.join(jobdir, fname), "w") as fh:
                fh.write(text)
        return jobdir


class TestSettingsJobdir(_JobDirCase):
    names = ["_JOB_ID_RE", "_settings_jobdir", "_new_job_id"]

    def test_valid_id_maps_to_settings_subdir(self):
        want = os.path.join(self.cwd, "Settings", _GOOD_ID)
        self.assertEqual(self.ns["_settings_jobdir"](_GOOD_ID), want)  # not yet created
        os.makedirs(want)
        self.assertEqual(self.ns["_settings_jobdir"](_GOOD_ID), want)

    def test_app_cwd_with_trailing_slash(self):
        # app.py sets current_working_directory = f"{os.getcwd()}/"
        self.ns["current_working_directory"] = self.cwd + "/"
        self.assertEqual(self.ns["_settings_jobdir"](_GOOD_ID),
                         os.path.join(self.cwd + "/", "Settings", _GOOD_ID))
        self.assertIsNone(self.ns["_settings_jobdir"]("../outside"))
        self.assertIsNone(self.ns["_settings_jobdir"](_LINK_ID))

    def test_ids_from_new_job_id_are_accepted(self):
        # the guard must accept everything the generator makes, or Retry and the
        # progress panel would silently stop working for real jobs
        for _ in range(200):
            jid = self.ns["_new_job_id"]()
            self.assertEqual(self.ns["_settings_jobdir"](jid),
                             os.path.join(self.cwd, "Settings", jid))

    def test_traversal_absolute_and_malformed_ids_rejected(self):
        for bad in self.bad_ids():
            with self.subTest(job_id=bad):
                self.assertIsNone(self.ns["_settings_jobdir"](bad))

    def test_symlink_escaping_settings_rejected(self):
        self.assertTrue(os.path.isfile(os.path.join(self.settings, _LINK_ID, "cmd.json")))
        self.assertIsNone(self.ns["_settings_jobdir"](_LINK_ID))

    def test_symlinked_settings_root_still_accepted(self):
        # a data dir whose Settings/ is itself a symlink (e.g. onto a volume) is a
        # legitimate layout: job dirs inside it must still resolve
        volume = os.path.join(self.cwd, "volume", "Settings")
        os.makedirs(os.path.join(volume, _GOOD_ID))
        data = os.path.join(self.cwd, "data")
        os.makedirs(data)
        os.symlink(volume, os.path.join(data, "Settings"))
        self.ns["current_working_directory"] = data
        self.assertEqual(self.ns["_settings_jobdir"](_GOOD_ID),
                         os.path.join(data, "Settings", _GOOD_ID))


class TestJobCallbacksRefuseBadIds(_JobDirCase):
    """Retry, Show-full-log and the progress poll, undecorated: a bad id is
    refused before anything is read or run; a real job id works as before."""

    names = [
        "_JOB_ID_RE",
        "_settings_jobdir",
        "_relaunch_job",
        "_full_job_log",
        "_retry_settings_job",
        "_toggle_settings_full_log",
        "refresh_settings_job",
    ]
    stubs = {
        "_render_failure": lambda stage, jobdir, job_id: ("failure", stage, jobdir, job_id),
        "_genotypes_degraded_warning": lambda jobdir: None,
        "_render_all_tables": lambda: "tables",
        "_render_storage": lambda: "storage",
        "_read_progress": lambda jobdir: None,
        "_render_progress_bar": lambda stage, prog: None,
        "_JOB_STALE_SECONDS": 1800,
    }

    def test_retry_refuses_bad_ids_without_reading_or_running(self):
        for bad in self.bad_ids():
            with self.subTest(job_id=bad):
                try:
                    out = self.ns["_retry_settings_job"](1, bad)
                except _PreventUpdate:
                    out = None
                if out is not None:  # the neutral "could not retry" message
                    self.assertIs(out[0], _NO_UPDATE)
                    self.assertIs(out[1], _NO_UPDATE)
        self.assertEqual(self.launched, [])
        self.assertEqual(self.outside_reads(), [])

    def test_retry_valid_id_relaunches_recorded_command(self):
        cmd = "crisprme.py download --what index --index-name NRG_3_hg38 --path /data"
        self.make_job(_GOOD_ID, {"cmd.json": json.dumps({"cmd": cmd, "stage": "Download index"})})
        out = self.ns["_retry_settings_job"](1, _GOOD_ID)
        self.assertEqual(out[0], "set_NEWJOB0001")
        self.assertIs(out[1], False)
        self.assertEqual(self.launched, [(cmd, "Download index")])

    def test_show_full_log_refuses_bad_ids_without_reading(self):
        for bad in self.bad_ids():
            with self.subTest(job_id=bad):
                self.assertEqual(self.ns["_toggle_settings_full_log"](1, True, bad),
                                 (False, "No job to show."))
        self.assertEqual(self.outside_reads(), [])

    def test_show_full_log_valid_id(self):
        self.make_job(_GOOD_ID, {"log_error.txt": "Traceback: boom\n"})
        hidden, text = self.ns["_toggle_settings_full_log"](1, True, _GOOD_ID)
        self.assertIs(hidden, False)
        self.assertIn("boom", text)

    def test_progress_poll_refuses_bad_ids_without_reading(self):
        for bad in self.bad_ids():
            with self.subTest(job_id=bad):
                with self.assertRaises(_PreventUpdate):
                    self.ns["refresh_settings_job"](1, bad)
        self.assertEqual(self.opened, [])

    def test_progress_poll_valid_id_unchanged(self):
        jobdir = self.make_job(_GOOD_ID, {
            "log.txt": "Download index\tQueued\nDownload index\tStart\nDownload index\tFAILED\n",
        })
        out = self.ns["refresh_settings_job"](1, _GOOD_ID)
        self.assertEqual(out[0], ("failure", "Download index", jobdir, _GOOD_ID))
        self.assertIs(out[1], True)
        self.assertIsNone(out[4])


if __name__ == "__main__":
    unittest.main()
