"""Unit tests for the 2.8.1 web Settings download-progress + failure-relay
helpers in pages/settings_page.py.

The Settings page imports a live Dash ``app`` (and dash/dbc), which the
network-free CI does NOT install, so — like test_history_page_assembly.py /
test_sanitize_job_name.py — the pure helper functions are EXTRACTED from source
via ast and exec'd with a tiny namespace (stdlib + light stubs for html/dbc).
This keeps the progress-JSON parsing, ETA formatting, full-log assembly, the
genotypes-degraded detector, and the retry-command round-trip covered in CI
without a Dash runtime. The job-id guard (_settings_jobdir) is covered the same
way, through the undecorated Retry / Show-full-log / progress callbacks, and so
is the genome .tar.gz upload's member vetting (_finalize_upload), run on
archives the tests build.
"""

import ast
import gzip
import hashlib
import io
import json
import os
import random
import re
import shutil
import stat
import string
import tarfile
import tempfile
import threading
import time
import unittest
from typing import List, Optional

import genome_in_use  # (PostProcess/: the rule _genome_upload_problem reads)

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


def _publish_stage():
    """The stage publish_index_cb launches its job with. That callback is only
    defined under ``if MAINTAINER_MODE:``, so it is read from the AST."""
    with open(_SETTINGS_PAGE) as fh:
        tree = ast.parse(fh.read())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "publish_index_cb":
            for call in ast.walk(node):
                if (isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                        and call.func.id == "launch_settings_job"):
                    return call.args[1].value
    raise AssertionError("publish_index_cb's launch_settings_job(...) call not found")


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
            MAINTAINER_MODE=False,
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

    def test_symlink_to_settings_itself_rejected(self):
        # a job-shaped entry that resolves to Settings/ itself is not a job dir:
        # its log.txt / cmd.json would be Settings/'s own
        for job_id, target in (("set_ROOTLINK00", self.settings), ("set_DOTLINK000", ".")):
            with self.subTest(target=target):
                link = os.path.join(self.settings, job_id)
                os.symlink(target, link)
                self.assertEqual(os.path.realpath(link), os.path.realpath(self.settings))
                self.assertIsNone(self.ns["_settings_jobdir"](job_id))

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
    refused before anything is read or run; a real job id works as before.
    Retry re-runs only a job the FAILED banner is shown for, and never a
    'Publish index' job outside MAINTAINER_MODE."""

    names = [
        "_JOB_ID_RE",
        "_settings_jobdir",
        "_job_failed",
        "_relaunch_job",
        "_full_job_log",
        "_retry_settings_job",
        "_toggle_settings_full_log",
        "refresh_settings_job",
    ]
    stubs = {
        "_render_failure": lambda stage, jobdir, job_id: ("failure", stage, jobdir, job_id),
        "_genotypes_degraded_warning": lambda jobdir: None,
        "_build_dirs_left_warning": lambda jobdir: None,
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

    def test_retry_failed_job_relaunches_recorded_command(self):
        cmd = "crisprme.py download --what index --index-name NRG_3_hg38 --path /data"
        self.make_job(_GOOD_ID, {
            "log.txt": "Download index\tQueued\nDownload index\tStart\nDownload index\tFAILED\n",
            "cmd.json": json.dumps({"cmd": cmd, "stage": "Download index"}),
        })
        out = self.ns["_retry_settings_job"](1, _GOOD_ID)
        self.assertEqual(out[0], "set_NEWJOB0001")
        self.assertIs(out[1], False)
        self.assertEqual(self.launched, [(cmd, "Download index")])

    def test_retry_only_reruns_a_job_the_banner_shows_failed(self):
        # Retry is offered only on the FAILED banner, but the job id is client
        # data: a job in any other state must not be re-run, and the decision
        # must match what the progress poll renders for the same log.txt
        cases = [
            ("failed", "S\tQueued\nS\tStart\nS\tFAILED\n", True),
            ("failed, no Queued marker", "S\tStart\nS\tFAILED\n", True),
            ("no log.txt", None, False),
            ("empty log.txt", "", False),
            ("queued", "S\tQueued\n", False),
            ("running", "S\tQueued\nS\tStart\n", False),
            ("done", "S\tQueued\nS\tStart\nS\tEnd\n", False),
            ("End beside FAILED (shown as done)", "S\tStart\nS\tFAILED\nS\tEnd\n", False),
            ("FAILED for another stage", "S\tQueued\nS\tStart\nT\tFAILED\n", False),
        ]
        for i, (label, log, retry) in enumerate(cases):
            with self.subTest(label):
                job_id = f"set_STATE{i:05d}"
                files = {"cmd.json": json.dumps({"cmd": f"echo job{i}", "stage": "S"})}
                if log is not None:
                    files["log.txt"] = log
                jobdir = self.make_job(job_id, files)
                try:
                    banner = self.ns["refresh_settings_job"](1, job_id)[0]
                except _PreventUpdate:
                    banner = None
                self.assertEqual(banner == ("failure", "S", jobdir, job_id), retry)
                del self.launched[:]
                out = self.ns["_retry_settings_job"](1, job_id)
                if retry:
                    self.assertEqual(self.launched, [(f"echo job{i}", "S")])
                    self.assertEqual(out[0], "set_NEWJOB0001")
                else:
                    self.assertEqual(self.launched, [])
                    self.assertIs(out[0], _NO_UPDATE)  # the neutral "could not retry"

    def test_retry_refuses_publish_index_outside_maintainer_mode(self):
        # publish_index_cb only exists in MAINTAINER_MODE: a publish recorded
        # while it was on must not be re-run through Retry once it is off
        stage = _publish_stage()
        cmd = "crisprme.py publish-index --index /data/genome_library/NRG_3_hg38 --path /data"
        self.make_job(_GOOD_ID, {
            "log.txt": f"{stage}\tQueued\n{stage}\tStart\n{stage}\tFAILED\n",
            "cmd.json": json.dumps({"cmd": cmd, "stage": stage}),
        })
        out = self.ns["_retry_settings_job"](1, _GOOD_ID)
        self.assertEqual(self.launched, [])
        self.assertIs(out[0], _NO_UPDATE)
        self.ns["MAINTAINER_MODE"] = True  # a maintainer may retry a failed publish
        out = self.ns["_retry_settings_job"](1, _GOOD_ID)
        self.assertEqual(self.launched, [(cmd, stage)])
        self.assertEqual(out[0], "set_NEWJOB0001")

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


# --- genome .tar.gz upload -----------------------------------------------------
# _finalize_upload unpacks a browser-uploaded genome archive into
# Genomes/<assembly>/. Its member names are client data: one that is absolute,
# has a '..' component, is a link, device or FIFO, or resolves outside that
# folder must refuse the WHOLE archive before anything is written (unvetted, a
# member could plant a job-shaped Settings/<id>/cmd.json for Retry to run).
_FA1 = b">chr1\nACGTACGTNNNN\n"
_FA2 = b">chr2\nTTTTGGGG\n"
# the start of a ._<name> member macOS tar writes (an AppleDouble file: not UTF-8)
_APPLEDOUBLE = b"\x00\x05\x16\x07\x00\x02\x00\x00Mac OS X        \x00\x02\x00\x00\x00\t\x8e\x00"


def _write_tar_gz(path, entries):
    """Write ``entries``, (name, kind, payload) tuples, as a .tar.gz at ``path``.
    kind: "file" (payload = bytes), "dir", "symlink" / "hardlink" (payload = the
    link target), "chrdev", "blkdev" or "fifo"."""
    with tarfile.open(path, "w:gz") as tar:
        for name, kind, payload in entries:
            info = tarfile.TarInfo(name)
            info.mtime = 1_700_000_000
            data = None
            if kind == "file":
                info.size, data = len(payload), io.BytesIO(payload)
            elif kind == "dir":
                info.type, info.mode = tarfile.DIRTYPE, 0o755
            elif kind in ("symlink", "hardlink"):
                info.type = tarfile.SYMTYPE if kind == "symlink" else tarfile.LNKTYPE
                info.linkname = payload
            elif kind in ("chrdev", "blkdev"):
                info.type = tarfile.CHRTYPE if kind == "chrdev" else tarfile.BLKTYPE
                info.devmajor, info.devminor = 1, 3
            else:
                info.type = tarfile.FIFOTYPE
            tar.addfile(info, data)


def _tree(root, meta=False):
    """{path relative to root: what is there} for everything under ``root``,
    links not followed; ``meta`` adds permission bits and a file's mtime."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        for name in dirnames + filenames:
            path = os.path.join(dirpath, name)
            st = os.lstat(path)
            if stat.S_ISLNK(st.st_mode):
                what = ("symlink", os.readlink(path))
            elif stat.S_ISDIR(st.st_mode):
                what = ("dir",)
            elif stat.S_ISREG(st.st_mode):
                with open(path, "rb") as fh:
                    what = ("file", fh.read(), st.st_nlink)
            else:  # a device or FIFO: never opened
                what = ("special", stat.S_IFMT(st.st_mode))
            if meta:
                what += (stat.S_IMODE(st.st_mode),)
                if stat.S_ISREG(st.st_mode):
                    what += (int(st.st_mtime),)
            out[os.path.relpath(path, root)] = what
    return out


class TestGenomeTarUpload(unittest.TestCase):
    """The genome .tar.gz branch of _finalize_upload, run on archives built here:
    a bad member refuses the whole archive with nothing written anywhere, and a
    normal genome archive extracts exactly as tarfile's own extractall does
    (which is all the code did before)."""

    # Each bad case comes after a benign chroms/chr1.fa, so an extractor that
    # vets member by member would already have written something.
    BAD = [
        ("'..' component (../x)", lambda out, dest: [("../x", "file", b"x")]),
        ("'..' out to outside/",
         lambda out, dest: [("chroms/../../../../outside/x", "file", b"x")]),
        ("'..' that stays inside", lambda out, dest: [("chroms/../chr2.fa", "file", _FA2)]),
        ("absolute path", lambda out, dest: [(os.path.join(out, "abs.txt"), "file", b"x")]),
        ("absolute path inside dest",
         lambda out, dest: [(os.path.join(dest, "chr2.fa"), "file", _FA2)]),
        ("symlink", lambda out, dest: [("lnk", "symlink", out), ("lnk/x", "file", b"x")]),
        ("relative symlink",
         lambda out, dest: [("lnk", "symlink", "../../../outside"), ("lnk/x", "file", b"x")]),
        ("hardlink", lambda out, dest: [
            ("hl", "hardlink", os.path.join(out, "secret.txt")), ("hl", "file", b"PWNED\n")]),
        ("character device", lambda out, dest: [("null", "chrdev", None)]),
        ("block device", lambda out, dest: [("blk", "blkdev", None)]),
        ("fifo", lambda out, dest: [("pipe", "fifo", None)]),
    ]
    GOOD = [
        ("chroms/ folder, like hg38.chromFa.tar.gz", [
            ("chroms", "dir", None),
            ("chroms/chr1.fa", "file", _FA1),
            ("chroms/chr2.fa", "file", _FA2),
        ]),
        ("files at the top, like hg19 chromFa.tar.gz", [
            ("chr1.fa", "file", _FA1),
            ("chr2.fa", "file", _FA2),
        ]),
        ("./ names, like tar -C dir .", [
            (".", "dir", None),
            ("./chr1.fa", "file", _FA1),
            ("./alt", "dir", None),
            ("./alt/chr2.fa", "file", _FA2),
        ]),
    ]

    def sandbox(self):
        """A fresh root holding the app data dir (data/, with Settings/uploads/)
        and, beside it, outside/ with a file a hardlink member could target.
        Returns (root, cwd, outside, ns)."""
        root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, root, ignore_errors=True)
        cwd = os.path.join(root, "data")
        outside = os.path.join(root, "outside")
        os.makedirs(os.path.join(cwd, "Settings", "uploads"))
        os.makedirs(outside)
        with open(os.path.join(outside, "secret.txt"), "w") as fh:
            fh.write("ORIGINAL\n")
        ns = _load_functions(
            ["_validate_name", "_unusable_in_shell", "_unsafe_tar_member",
             "_genome_upload_assembly", "_built_on_genome", "_genome_upload_problem",
             "_install_uploaded_genome", "_finalize_upload"],
            current_working_directory=cwd,
            gzip=gzip,
            tarfile=tarfile,
            shutil=shutil,
            tempfile=tempfile,
            genome_in_use=genome_in_use,
        )
        return root, cwd, outside, ns

    def upload(self, cwd, name, entries):
        """The finished .part file the chunked-upload route hands over."""
        part = os.path.join(cwd, "Settings", "uploads", name + ".part")
        _write_tar_gz(part, entries)
        return part

    def test_unsafe_member_refuses_whole_archive_writing_nothing(self):
        for label, bad in self.BAD:
            with self.subTest(label):
                root, cwd, outside, ns = self.sandbox()
                dest = os.path.join(cwd, "Genomes", "evil")
                part = self.upload(cwd, "evil.tar.gz",
                                   [("chroms/chr1.fa", "file", _FA1)] + bad(outside, dest))
                before = _tree(root)
                with self.assertRaises(ValueError):
                    ns["_finalize_upload"]("genome", part, "evil.tar.gz", "", "")
                self.assertEqual(_tree(root), before)  # nothing written, anywhere
                self.assertFalse(os.path.exists(os.path.join(cwd, "Genomes")))

    def test_member_through_a_symlink_already_in_dest_refused(self):
        # a plain relative member: only the realpath check sees that
        # Genomes/<assembly>/chroms already points out of the assembly folder
        root, cwd, outside, ns = self.sandbox()
        dest = os.path.join(cwd, "Genomes", "evil")
        os.makedirs(dest)
        os.symlink(outside, os.path.join(dest, "chroms"))
        part = self.upload(cwd, "evil.tar.gz", [("chroms/chr1.fa", "file", _FA1)])
        before = _tree(root)
        with self.assertRaises(ValueError):
            ns["_finalize_upload"]("genome", part, "evil.tar.gz", "", "")
        self.assertEqual(_tree(root), before)

    def test_normal_genome_archive_extracts_flat(self):
        # every FASTA lands at the top of Genomes/<assembly>/, as index-genome and the
        # search take only those (a UCSC chromFa.tar.gz holds them in chroms/): the
        # files tarfile's own extractall writes, byte for byte, mode and mtime too
        for label, members in self.GOOD:
            with self.subTest(label):
                root, cwd, _outside, ns = self.sandbox()
                part = self.upload(cwd, "myasm.chromFa.tar.gz", members)
                ref = os.path.join(root, "reference")
                os.makedirs(ref)
                with gzip.open(part, "rb") as fin, tarfile.open(fileobj=fin, mode="r") as tar:
                    tar.extractall(ref)
                want = {os.path.basename(path): what
                        for path, what in _tree(ref, meta=True).items() if what[0] == "file"}
                got = ns["_finalize_upload"]("genome", part, "myasm.chromFa.tar.gz", "", "")
                self.assertEqual(got, "myasm")
                self.assertFalse(os.path.exists(part))
                dest = os.path.join(cwd, "Genomes", "myasm")
                self.assertEqual(sorted(want), ["chr1.fa", "chr2.fa"])
                self.assertEqual(_tree(dest, meta=True), want)
                self.assertEqual(os.listdir(os.path.join(cwd, "Genomes")), ["myasm"])

    def test_symlinked_genomes_root_still_accepted(self):
        # Genomes/ on another volume through a symlink is a legitimate layout:
        # containment is checked against the real path of Genomes/<assembly>/
        root, cwd, _outside, ns = self.sandbox()
        volume = os.path.join(root, "volume")
        os.makedirs(volume)
        os.symlink(volume, os.path.join(cwd, "Genomes"))
        part = self.upload(cwd, "myasm.tar.gz", self.GOOD[0][1])
        self.assertEqual(ns["_finalize_upload"]("genome", part, "myasm.tar.gz", "", ""), "myasm")
        with open(os.path.join(volume, "myasm", "chr1.fa"), "rb") as fh:
            self.assertEqual(fh.read(), _FA1)

    # --- the layout every search needs: one FASTA a sequence, at the top -------------
    # (a Direct URL download gives it: utils.download_reference_genome)
    def plain_upload(self, cwd, name, data):
        part = os.path.join(cwd, "Settings", "uploads", name + ".part")
        with open(part, "wb") as fh:
            fh.write(data)
        return part

    def assert_one_file_a_sequence(self, cwd, assembly, want):
        dest = os.path.join(cwd, "Genomes", assembly)
        got = {}
        for f in sorted(os.listdir(dest)):
            with open(os.path.join(dest, f), "rb") as fh:
                got[f] = fh.read()
        self.assertEqual(got, want)
        # nothing else in Genomes/: no staging folder left, no other genome
        self.assertEqual(os.listdir(os.path.join(cwd, "Genomes")), [assembly])

    def test_multi_record_fasta_upload_is_split(self):
        both = _FA1 + _FA2
        for name, data in (("walkdel.fa", both), ("walkdel.fasta", both),
                           ("walkdel.fa.gz", gzip.compress(both))):
            with self.subTest(name):
                root, cwd, _outside, ns = self.sandbox()
                part = self.plain_upload(cwd, name, data)
                self.assertEqual(ns["_finalize_upload"]("genome", part, name, "", ""), "walkdel")
                self.assertFalse(os.path.exists(part))
                self.assert_one_file_a_sequence(cwd, "walkdel", {"chr1.fa": _FA1, "chr2.fa": _FA2})

    def test_single_record_fasta_upload_is_kept_as_is(self):
        root, cwd, _outside, ns = self.sandbox()
        part = self.plain_upload(cwd, "chrM.fa", b">chrM\nACGT\n")
        self.assertEqual(ns["_finalize_upload"]("genome", part, "chrM.fa", "", ""), "chrM")
        self.assert_one_file_a_sequence(cwd, "chrM", {"chrM.fa": b">chrM\nACGT\n"})

    def test_single_sequence_upload_is_named_after_its_sequence(self):
        # a one-chromosome subset or a bacterial genome, named otherwise: it went in as
        # walkone.fa, and every search of it found its off-targets and then dropped them
        # all, a search taking each file for the sequence it is named after (round 3
        # review, F2); the build of its index said done
        cases = [
            ("walkone.fa", self.plain_upload, _FA1, {"chr1.fa": _FA1}),
            ("walkone.fa.gz", self.plain_upload, gzip.compress(_FA1), {"chr1.fa": _FA1}),
            ("walkone.tar.gz", self.upload, [("walkone/hg38_chr1.fa", "file", _FA1),
                                             ("walkone/chr2.fa", "file", _FA2)],
             {"chr1.fa": _FA1, "chr2.fa": _FA2}),
        ]
        for name, upload, data, want in cases:
            with self.subTest(name):
                root, cwd, _outside, ns = self.sandbox()
                part = upload(cwd, name, data)
                self.assertEqual(ns["_finalize_upload"]("genome", part, name, "", ""), "walkone")
                self.assert_one_file_a_sequence(cwd, "walkone", want)
        # a sequence id no file can be named after is refused, as a split refuses one
        root, cwd, _outside, ns = self.sandbox()
        part = self.plain_upload(cwd, "walkone.fa", b">chr1/a\nACGT\n")
        with self.assertRaises(ValueError) as cm:
            ns["_finalize_upload"]("genome", part, "walkone.fa", "", "")
        self.assertIn("sequence id 'chr1/a' cannot be used as a file name", str(cm.exception))
        self.assertEqual(os.listdir(os.path.join(cwd, "Genomes")), [])

    def test_header_of_more_than_its_id_upload_goes_in_as_the_id(self):
        # an NCBI, GENCODE or Ensembl FASTA, or CRLF line ends, a file already named after
        # its sequence among them: it went in as it was, with an index built as done, and
        # every search of it lost results, often all of them (round 5 review, F1)
        for name, data, want in (
                ("ecoli.fa", b">NC_000913.3 Escherichia coli K-12\nACGT\nAC\n",
                 {"NC_000913.3.fa": b">NC_000913.3\nACGT\nAC\n"}),
                ("gencode.fa.gz", gzip.compress(b">chr1 1\nACGT\n>chr2 2\nGG\n"),
                 {"chr1.fa": b">chr1\nACGT\n", "chr2.fa": b">chr2\nGG\n"}),
                ("chr1.fa", b">chr1 1\nACGT\n", {"chr1.fa": b">chr1\nACGT\n"}),
                ("chr1.fa", b">chr1\r\nACGT\r\nAC\r\n", {"chr1.fa": b">chr1\nACGT\nAC\n"}),
                # a CR that only ends the file: a build reads it, or a variant build does
                # (integration round 2 review, F1)
                ("chr1.fa", b">chr1\nACGT\nAC\r\n", {"chr1.fa": b">chr1\nACGT\nAC\n"}),
                ("chr1.fa", b">chr1\nACGT\nAC\n\r", {"chr1.fa": b">chr1\nACGT\nAC\n"}),
                ("chr1.fa", b">chr1\nACGT\nAC\r\n\r\n", {"chr1.fa": b">chr1\nACGT\nAC\n\n"})):
            with self.subTest(name, data=data):
                root, cwd, _outside, ns = self.sandbox()
                part = self.plain_upload(cwd, name, data)
                assembly = name.split(".")[0]
                self.assertEqual(ns["_finalize_upload"]("genome", part, name, "", ""), assembly)
                self.assert_one_file_a_sequence(cwd, assembly, want)
        # in an archive, a file already named after its sequence too
        root, cwd, _outside, ns = self.sandbox()
        part = self.upload(cwd, "asm.tar.gz", [("chroms/chr1.fa", "file", b">chr1 1\nACGT\n"),
                                               ("chroms/chr2.fa", "file", _FA2)])
        self.assertEqual(ns["_finalize_upload"]("genome", part, "asm.tar.gz", "", ""), "asm")
        self.assert_one_file_a_sequence(cwd, "asm", {"chr1.fa": b">chr1\nACGT\n", "chr2.fa": _FA2})
        # an id the pipeline cannot take is refused, with nothing installed
        root, cwd, _outside, ns = self.sandbox()
        part = self.plain_upload(cwd, "ecoli.fa", b">gi|49175990|ref|NC_000913.2| E. coli\nACGT\n")
        with self.assertRaises(ValueError) as cm:
            ns["_finalize_upload"]("genome", part, "ecoli.fa", "", "")
        self.assertIn("sequence id 'gi|49175990|ref|NC_000913.2|' holds '|'", str(cm.exception))
        self.assertEqual(os.listdir(os.path.join(cwd, "Genomes")), [])

    def test_two_files_of_one_sequence_refused_installing_nothing(self):
        # named after its sequence, one file would take the other's name
        for label, entries in (
                ("another file of its name", [("chr1.fa", "file", _FA1),
                                              ("extra.fa", "file", b">chr1\nGGGG\n")]),
                ("a split piece of its name", [("asm.fa", "file", _FA1 + _FA2),
                                               ("more.fa", "file", b">chr2\nGGGG\n")])):
            with self.subTest(label):
                root, cwd, _outside, ns = self.sandbox()
                part = self.upload(cwd, "asm.tar.gz", entries)
                with self.assertRaises(ValueError) as cm:
                    ns["_finalize_upload"]("genome", part, "asm.tar.gz", "", "")
                self.assertRegex(str(cm.exception), r"^(extra|more)\.fa: cannot name it after "
                                 r"its sequence -- chr[12]\.fa already exists$")
                self.assertEqual(os.listdir(os.path.join(cwd, "Genomes")), [])

    def test_archive_member_with_several_sequences_is_split(self):
        root, cwd, _outside, ns = self.sandbox()
        part = self.upload(cwd, "asm.tar.gz", [("asm/asm.fa", "file", _FA1 + _FA2),
                                               ("asm/README", "file", b"not a genome\n")])
        self.assertEqual(ns["_finalize_upload"]("genome", part, "asm.tar.gz", "", ""), "asm")
        self.assert_one_file_a_sequence(cwd, "asm", {"chr1.fa": _FA1, "chr2.fa": _FA2})

    def test_archive_of_fasta_files_installs_them_as_fa(self):
        # index-genome indexes the .fa files only, and a search lists a genome's
        # chromosomes as its *.fa: kept as chr1.fasta, an archive's genome was given an
        # index with no .bin, and then a build that refused it said to upload it again,
        # which laid it out the same (round 2 review, R2-1)
        cases = [
            ("files at the top", [("chr1.fasta", "file", _FA1), ("chr2.fasta", "file", _FA2)]),
            ("chroms/ folder of .fasta and .fa", [("chroms", "dir", None),
                                                  ("chroms/chr1.fasta", "file", _FA1),
                                                  ("chroms/chr2.fa", "file", _FA2)]),
            ("a .fasta of several sequences", [("asm/asm.fasta", "file", _FA1 + _FA2)]),
        ]
        for label, entries in cases:
            with self.subTest(label):
                root, cwd, _outside, ns = self.sandbox()
                part = self.upload(cwd, "walkfs2.tar.gz", entries)
                self.assertEqual(
                    ns["_finalize_upload"]("genome", part, "walkfs2.tar.gz", "", ""), "walkfs2")
                self.assert_one_file_a_sequence(cwd, "walkfs2",
                                                {"chr1.fa": _FA1, "chr2.fa": _FA2})

    def test_macos_tar_attribute_files_are_left_out(self):
        # macOS tar writes each file's extended attributes as a ._<name> member beside
        # it; installed, ._chr1.fa was read as a chromosome by every search, which
        # failed (UnicodeDecodeError), though its build said done
        cases = [
            ("files at the top", [("._chr1.fa", "file", _APPLEDOUBLE), ("chr1.fa", "file", _FA1),
                                  ("._chr2.fa", "file", _APPLEDOUBLE), ("chr2.fa", "file", _FA2)]),
            ("chroms/ folder", [("._chroms", "file", _APPLEDOUBLE), ("chroms", "dir", None),
                                ("chroms/._chr1.fa", "file", _APPLEDOUBLE),
                                ("chroms/chr1.fa", "file", _FA1),
                                ("chroms/._chr2.fasta", "file", _APPLEDOUBLE),
                                ("chroms/chr2.fasta", "file", _FA2)]),
        ]
        for label, entries in cases:
            with self.subTest(label):
                root, cwd, _outside, ns = self.sandbox()
                part = self.upload(cwd, "walkad.tar.gz", entries)
                self.assertEqual(
                    ns["_finalize_upload"]("genome", part, "walkad.tar.gz", "", ""), "walkad")
                self.assert_one_file_a_sequence(cwd, "walkad",
                                                {"chr1.fa": _FA1, "chr2.fa": _FA2})

    def test_unusable_upload_refused_installing_nothing(self):
        cases = [
            ("no FASTA in the archive", "asm.tar.gz", [("asm/README", "file", b"x\n")], None),
            ("hidden FASTA files only", "asm.tar.gz", [("._chr1.fa", "file", _APPLEDOUBLE)],
             "no .fa or .fasta file in asm.tar.gz"),
            ("two FASTA files of one name", "asm.tar.gz",
             [("a/chr1.fa", "file", _FA1), ("b/chr1.fa", "file", _FA1)],
             "would both be installed as chr1.fa (a/chr1.fa, b/chr1.fa)"),
            ("a .fa and a .fasta of one name", "asm.tar.gz",
             [("chr1.fa", "file", _FA1), ("chroms/chr1.fasta", "file", _FA1)],
             "would both be installed as chr1.fa (chr1.fa, chroms/chr1.fasta)"),
            ("two sequences of one id", "asm.tar.gz", [("asm.fa", "file", _FA1 + _FA1)], None),
        ]
        for label, name, entries, want in cases:
            with self.subTest(label):
                root, cwd, _outside, ns = self.sandbox()
                part = self.upload(cwd, name, entries)
                with self.assertRaises(ValueError) as cm:
                    ns["_finalize_upload"]("genome", part, name, "", "")
                if want:
                    self.assertIn(want, str(cm.exception))
                self.assertEqual(os.listdir(os.path.join(cwd, "Genomes")), [])

    def test_upload_replaces_a_genome_of_the_same_name_whole(self):
        root, cwd, _outside, ns = self.sandbox()
        old = os.path.join(cwd, "Genomes", "asm")
        os.makedirs(old)
        with open(os.path.join(old, "chrOld.fa"), "wb") as fh:
            fh.write(b">chrOld\nAAAA\n")
        part = self.plain_upload(cwd, "asm.fa", _FA1 + _FA2)
        ns["_finalize_upload"]("genome", part, "asm.fa", "", "")
        self.assert_one_file_a_sequence(cwd, "asm", {"chr1.fa": _FA1, "chr2.fa": _FA2})

    # --- never under the name of a genome something installed was built on -----------
    # (final review round 1, F2): an index or a variant genome holds that genome's
    # sequences and positions, and searches and builds of the new one took it as it was
    BUILT_ON = [
        ("its index", ["genome_library/NRG_2_asm"]),
        ("a variant index of it", ["genome_library/NRG_3_asm+walk",
                                   "genome_library/NRG_3_asm+walk_INDELS"]),
        ("a variant genome of it", ["Genomes/asm+walk", "Genomes/asm+walk_INDELS"]),
    ]

    def test_upload_refused_while_something_is_built_on_a_genome_of_its_name(self):
        for label, built in self.BUILT_ON:
            # the genome there (replaced), or removed before as the 2.8.1 notes say to
            # do with one uploaded before (its index of the one FASTA it was is still
            # there, and gave chr2's sites at chr1)
            for there in (True, False):
                with self.subTest(label, genome_there=there):
                    root, cwd, _outside, ns = self.sandbox()
                    if there:
                        os.makedirs(os.path.join(cwd, "Genomes", "asm"))
                        with open(os.path.join(cwd, "Genomes", "asm", "asm.fa"), "wb") as fh:
                            fh.write(_FA1 + _FA2)
                    for d in built:
                        os.makedirs(os.path.join(cwd, d))
                        with open(os.path.join(cwd, d, "NRG_chr1_1.bin"), "wb") as fh:
                            fh.write(b"x")
                    part = self.plain_upload(cwd, "asm.fa", _FA2 + _FA1)
                    before = _tree(cwd)
                    with self.assertRaises(ValueError) as cm:
                        ns["_finalize_upload"]("genome", part, "asm.fa", "", "")
                    self.assertTrue(str(cm.exception).startswith(
                        "Indexes or variant genomes built on a genome named 'asm' are "
                        f"installed ({', '.join(built)}): "), cm.exception)
                    self.assertIn("Remove them first (an index under 'Remove installed "
                                  "data'; a Genomes/asm+<dataset> folder from your data "
                                  "folder), or upload this genome under another name.",
                                  str(cm.exception))
                    self.assertEqual(_tree(cwd), before)  # nothing changed

    def test_a_direct_url_download_is_refused_alike_with_its_own_way_out(self):
        # (add_genome says so at once; crisprme.py download refuses it with the same rule)
        root, cwd, _outside, ns = self.sandbox()
        os.makedirs(os.path.join(cwd, "genome_library", "NRG_2_asm"))
        said = ns["_genome_upload_problem"]("asm", "download")
        self.assertTrue(said.startswith(
            "Indexes or variant genomes built on a genome named 'asm' are installed "
            "(genome_library/NRG_2_asm): searches of a genome downloaded under that name "
            "would go on using them"), said)
        self.assertTrue(said.endswith("Remove them first (an index under 'Remove installed "
                                      "data'; a Genomes/asm+<dataset> folder from your data "
                                      "folder), or download this genome under another name."),
                        said)
        self.assertIsNone(ns["_genome_upload_problem"]("asm2", "download"))

    def test_what_is_built_on_other_genomes_does_not_refuse_it(self):
        root, cwd, _outside, ns = self.sandbox()
        for d in ("genome_library/NRG_2_asm2", "genome_library/NRG_2_hg38+asm",
                  "genome_library/NRG_2_myasm", "genome_library/.build-index-only-NRG_2_asm",
                  "Genomes/asm2+walk", "Genomes/hg38+asm"):
            os.makedirs(os.path.join(cwd, d))
        self.assertIsNone(ns["_genome_upload_problem"]("asm"))
        part = self.plain_upload(cwd, "asm.fa", _FA1 + _FA2)
        self.assertEqual(ns["_finalize_upload"]("genome", part, "asm.fa", "", ""), "asm")
        self.assertEqual(sorted(os.listdir(os.path.join(cwd, "Genomes", "asm"))),
                         ["chr1.fa", "chr2.fa"])

    def test_genome_linked_to_a_shared_copy_is_not_written_through(self):
        root, cwd, _outside, ns = self.sandbox()
        shared = os.path.join(root, "shared")
        os.makedirs(shared)
        with open(os.path.join(shared, "chrS.fa"), "wb") as fh:
            fh.write(b">chrS\nCCCC\n")
        os.makedirs(os.path.join(cwd, "Genomes"))
        os.symlink(shared, os.path.join(cwd, "Genomes", "asm"))
        part = self.plain_upload(cwd, "asm.fa", _FA1 + _FA2)
        ns["_finalize_upload"]("genome", part, "asm.fa", "", "")
        self.assertFalse(os.path.islink(os.path.join(cwd, "Genomes", "asm")))
        self.assert_one_file_a_sequence(cwd, "asm", {"chr1.fa": _FA1, "chr2.fa": _FA2})
        self.assertEqual(os.listdir(shared), ["chrS.fa"])


# --- the sample list of a VCF dataset added in Settings (Gate B walk P3) -------------
# A web variant search reads samplesIDs/<dataset>.samplesID.txt, and failed at its end
# for a dataset added only in Settings: nothing wrote one, and the Settings build never
# passed --samplesID, so the tiers per-sample analysis reads were never built either.
_SID_NAMES = ["SAMPLESID_HEADER", "_UNKNOWN_SUPERPOPULATION", "_SAMPLESID_LOCK",
              "_vcf_header_samples", "_own_lists_record", "_own_list", "_record_own_lists",
              "_installed_sample_metadata", "_listed_samples", "_write_dataset_samplesid",
              "_samplesid_listing"]
_HEADER = "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
# the list the shipped hg38_mega index installs: one pseudo-sample a source dataset. Its
# VCFs are sites-only, and it is built without --samplesID (docs/DESIGN_mega_index.md):
# its registry is made from their INFO allele frequencies (build_mega_registry.py)
_MEGA_LIST = _HEADER + "".join(f"{d}\t{d}\t{d}\tNA\n"
                               for d in ("1000G2021", "HGDP", "gnomAD", "TOPMed", "AoU"))


def _write_vcf(path, samples, records=("chr1\t10\t.\tA\tG\t.\tPASS\t.",)):
    """A gzipped VCF whose #CHROM header has these sample columns (none: sites-only)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        cols += ["FORMAT"] + list(samples)
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n" + "\t".join(cols) + "\n")
        for rec in records:
            fh.write(rec + ("\tGT" + "\t0|1" * len(samples) if samples else "") + "\n")


class _SamplesIDCase(unittest.TestCase):
    names = _SID_NAMES
    stubs = {}

    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)
        for d in ("VCFs", "samplesIDs", "Genomes/hg38", "Settings/uploads"):
            os.makedirs(os.path.join(self.cwd, d))
        self.ns = _load_functions(self.names, current_working_directory=self.cwd + "/",
                                  gzip=gzip, shutil=shutil, threading=threading, hashlib=hashlib,
                                  **self.stubs)

    def vcf(self, dataset, fname, samples, **kw):
        _write_vcf(os.path.join(self.cwd, "VCFs", dataset, fname), samples, **kw)

    def sid(self, dataset):
        return os.path.join(self.cwd, "samplesIDs", f"{dataset}.samplesID.txt")

    def read(self, path):
        with open(path) as fh:
            return fh.read()


class TestDatasetSamplesID(_SamplesIDCase):
    def test_list_written_from_the_vcf_headers(self):
        # the union of the files' samples, each once, in the order first seen (files
        # by name), in the shipped 4-column layout; hprc_build.sh's values for a
        # sample with no metadata: the dataset, unknown, unknown
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2", "S3"])
        self.vcf("walk", "walk.chr2.vcf.gz", ["S1", "S2", "S3"])
        self.vcf("walk", "walk.chrY.vcf.gz", ["S1", "S4"])
        got = self.ns["_write_dataset_samplesid"]("walk")
        self.assertEqual(got, self.sid("walk"))
        self.assertEqual(self.read(got), _HEADER + "".join(
            f"{s}\twalk\tunknown\tunknown\n" for s in ("S1", "S2", "S3", "S4")))

    def test_never_more_samples_than_the_vcfs_genotype(self):
        # a listed sample no VCF genotypes counts as hom-ref in every AN (#46)
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        self.ns["_write_dataset_samplesid"]("walk")
        rows = [ln.split("\t")[0] for ln in self.read(self.sid("walk")).splitlines()[1:]]
        self.assertEqual(rows, ["S1", "S2"])

    def test_sites_only_dataset_gets_no_list(self):
        self.vcf("agg", "agg.vcf.gz", [])
        self.assertIsNone(self.ns["_write_dataset_samplesid"]("agg"))
        self.assertFalse(os.path.exists(self.sid("agg")))

    def test_unreadable_or_empty_dataset_gets_no_list(self):
        os.makedirs(os.path.join(self.cwd, "VCFs", "empty"))
        bad = os.path.join(self.cwd, "VCFs", "bad", "bad.vcf.gz")
        os.makedirs(os.path.dirname(bad))
        with open(bad, "wb") as fh:
            fh.write(b"not gzip at all")
        for dataset in ("empty", "bad", "absent"):
            with self.subTest(dataset):
                self.assertIsNone(self.ns["_write_dataset_samplesid"](dataset))
                self.assertFalse(os.path.exists(self.sid(dataset)))

    def test_a_list_from_elsewhere_is_never_changed(self):
        # a list the user wrote, or a download installed, carries real populations
        mine = _HEADER + "S1\tGBR\tEUR\tmale\n"
        with open(self.sid("walk"), "w") as fh:
            fh.write(mine)
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        self.assertEqual(self.ns["_write_dataset_samplesid"]("walk"), self.sid("walk"))
        self.assertEqual(self.read(self.sid("walk")), mine)

    def test_sites_only_dataset_leaves_a_list_from_elsewhere_alone(self):
        # the mega index's list names datasets, not samples its VCFs genotype
        with open(self.sid("hg38_mega"), "w") as fh:
            fh.write(_MEGA_LIST)
        self.vcf("hg38_mega", "mega.chr1.vcf.gz", [])
        self.assertIsNone(self.ns["_write_dataset_samplesid"]("hg38_mega"))
        self.assertEqual(self.read(self.sid("hg38_mega")), _MEGA_LIST)

    def test_a_list_naming_no_sample_the_vcfs_genotype_is_not_passed(self):
        # tiers emitted from it would hold no sample, and so no allele frequency
        other = _HEADER + "X1\tGBR\tEUR\tmale\n"
        with open(self.sid("walk"), "w") as fh:
            fh.write(other)
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        self.assertIsNone(self.ns["_write_dataset_samplesid"]("walk"))
        self.assertEqual(self.read(self.sid("walk")), other)

    def test_a_list_of_its_header_alone_is_never_changed(self):
        # the page never writes one: someone put it there (as the Submit refusal of a
        # search said to), for a sites-only dataset or a genotyped one, and it stays
        # (final review round 1, F3); it names no sample, so no build passes it
        for dataset, samples in (("agg", []), ("walk", ["S1", "S2"])):
            with self.subTest(dataset=dataset):
                with open(self.sid(dataset), "w") as fh:
                    fh.write(_HEADER)
                self.vcf(dataset, f"{dataset}.chr1.vcf.gz", samples)
                self.assertIsNone(self.ns["_write_dataset_samplesid"](dataset))
                self.assertEqual(self.read(self.sid(dataset)), _HEADER)

    def test_own_list_follows_the_dataset(self):
        # a file of the dataset added after the list was written (the uploader takes
        # one file at a time): the list this page wrote is written again
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1"])
        self.ns["_write_dataset_samplesid"]("walk")
        self.vcf("walk", "walk.chr2.vcf.gz", ["S1", "S2"])
        self.ns["_write_dataset_samplesid"]("walk")
        self.assertEqual(self.read(self.sid("walk")),
                         _HEADER + "S1\twalk\tunknown\tunknown\nS2\twalk\tunknown\tunknown\n")
        # unchanged samples: not written again
        before = os.stat(self.sid("walk")).st_mtime_ns
        time.sleep(0.01)
        self.ns["_write_dataset_samplesid"]("walk")
        self.assertEqual(os.stat(self.sid("walk")).st_mtime_ns, before)
        # replaced by a sites-only file: its own list would name samples it lacks
        shutil.rmtree(os.path.join(self.cwd, "VCFs", "walk"))
        self.vcf("walk", "walk.vcf.gz", [])
        self.assertIsNone(self.ns["_write_dataset_samplesid"]("walk"))
        self.assertFalse(os.path.exists(self.sid("walk")))

    def test_registered_server_folder(self):
        # 'Register a server folder' links VCFs/<name> to the folder
        server = os.path.join(self.cwd, "server")
        _write_vcf(os.path.join(server, "cohort.chr1.vcf.gz"), ["A", "B"])
        os.symlink(server, os.path.join(self.cwd, "VCFs", "cohort"))
        self.assertEqual(self.ns["_write_dataset_samplesid"]("cohort"), self.sid("cohort"))
        self.assertIn("B\tcohort\tunknown\tunknown\n", self.read(self.sid("cohort")))

    def test_dataset_names_that_are_not_a_folder_of_vcfs_are_refused(self):
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1"])
        for bad in ("../walk", "walk/..", ".hidden", "", None):
            with self.subTest(bad=bad):
                self.assertIsNone(self.ns["_write_dataset_samplesid"](bad))
        self.assertEqual(os.listdir(os.path.join(self.cwd, "samplesIDs")), [])

    def test_listing_names_the_list_as_build_index_only_reads_one(self):
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1"])
        self.ns["_write_dataset_samplesid"]("walk")
        listing = self.ns["_samplesid_listing"]("walk")
        self.assertEqual(listing, os.path.join(self.cwd + "/", "samplesIDs",
                                               ".walk.samplesID.list"))
        self.assertEqual(self.read(listing), "walk.samplesID.txt\n")

    def test_two_threads_writing_the_list_at_once_write_it_one_after_the_other(self):
        # the web server runs each request in a thread of one process, so two uploads
        # into one dataset that finish together write its list at once: through one
        # temporary file, one truncated the other's and its rename found nothing there
        # (round 2 review, R2-2), and the record of the list this page wrote could name
        # the other's. One writes it at a time; the other then finds it written.
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        inside, overlaps, errors = [], [], []
        read_headers = self.ns["_vcf_header_samples"]

        def slow_read(vcf_dir):  # (called by each call, while it decides what to write)
            inside.append(1)
            if len(inside) > 1:
                overlaps.append(vcf_dir)
            time.sleep(0.3)
            try:
                return read_headers(vcf_dir)
            finally:
                inside.pop()

        def write():
            try:
                self.ns["_write_dataset_samplesid"]("walk")
            except Exception as exc:  # noqa: BLE001 - reported below
                errors.append(exc)

        self.ns["_vcf_header_samples"] = slow_read
        threads = [threading.Thread(target=write) for _ in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(30)
        self.assertEqual(errors, [])
        self.assertEqual(overlaps, [])
        self.assertEqual(self.read(self.sid("walk")), _HEADER
                         + "S1\twalk\tunknown\tunknown\nS2\twalk\tunknown\tunknown\n")
        self.assertEqual(self.ns["_own_list"]("walk"), self.read(self.sid("walk")))
        self.assertEqual(sorted(os.listdir(os.path.join(self.cwd, "samplesIDs"))),
                         [".walk.samplesID.sha256", "walk.samplesID.txt"])


# --- each sample with the metadata an installed list gives it (round 3 review, F1) ----
# The shareable report reads every list of samplesIDs/, the last (by name) naming a
# sample giving its superpopulation: a list written "<sample> <dataset> unknown unknown"
# for every sample turned a panel's samples unknown in the report of every search of
# that panel (S1..S6 SUPX/SUPY of hg38_walk, read before walkz).
_PANEL = _HEADER + ("S1\tPOPA\tSUPX\tmale\nS2\tPOPA\tSUPX\tfemale\nS3\tPOPB\tSUPY\tmale\n")


class TestSampleListTakesInstalledMetadata(_SamplesIDCase):
    def put(self, name, text):
        with open(os.path.join(self.cwd, "samplesIDs", name), "w") as fh:
            fh.write(text)

    def test_each_sample_takes_what_an_installed_list_gives_it(self):
        self.put("hg38_walk.samplesID.txt", _PANEL)
        # read after it: S3's row replaces the panel's, as the report reads them; S2's,
        # which names no superpopulation, gives nothing
        self.put("hg38_zz.samplesID.txt", "S3\tPOPZ\tSUPZ\tfemale\nS2\tzz\tunknown\tunknown\n")
        self.put("samplesIDs.config.txt", "S1\tBAD\tBAD\tBAD\n")  # a config, not a list
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1", "S2", "S3", "S9"])
        self.assertEqual(self.ns["_write_dataset_samplesid"]("walkz"), self.sid("walkz"))
        self.assertEqual(self.read(self.sid("walkz")), _HEADER
                         + "S1\tPOPA\tSUPX\tmale\nS2\tPOPA\tSUPX\tfemale\n"
                         + "S3\tPOPZ\tSUPZ\tfemale\nS9\twalkz\tunknown\tunknown\n")

    def test_the_report_keeps_the_superpopulation_a_panel_gives(self):
        # the report's own reader, on what the page wrote before the panel was installed
        # (its list then gives unknown) and once it writes it again (a build)
        try:
            import generate_report
        except ImportError as exc:  # (pandas, matplotlib)
            self.skipTest(f"generate_report not importable: {exc}")
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1", "S2", "S3", "S9"])
        self.ns["_write_dataset_samplesid"]("walkz")
        self.put("hg38_walk.samplesID.txt", _PANEL)  # read before walkz.samplesID.txt
        sdir = os.path.join(self.cwd, "samplesIDs")
        want = {"S1": "SUPX", "S2": "SUPX", "S3": "SUPY", "S9": "unknown"}
        self.assertEqual(generate_report.load_sample_superpop(sdir), want)
        self.ns["_write_dataset_samplesid"]("walkz")
        self.assertEqual(self.read(self.sid("walkz")), _HEADER
                         + "S1\tPOPA\tSUPX\tmale\nS2\tPOPA\tSUPX\tfemale\n"
                         + "S3\tPOPB\tSUPY\tmale\nS9\twalkz\tunknown\tunknown\n")
        self.assertEqual(generate_report.load_sample_superpop(sdir), want)

    def test_own_list_follows_the_installed_lists(self):
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1", "S2"])
        self.ns["_write_dataset_samplesid"]("walkz")
        placeholders = _HEADER + "S1\twalkz\tunknown\tunknown\nS2\twalkz\tunknown\tunknown\n"
        self.assertEqual(self.read(self.sid("walkz")), placeholders)
        self.put("hg38_walk.samplesID.txt", _PANEL)  # a panel installed since
        self.ns["_write_dataset_samplesid"]("walkz")
        self.assertEqual(self.read(self.sid("walkz")),
                         _HEADER + "S1\tPOPA\tSUPX\tmale\nS2\tPOPA\tSUPX\tfemale\n")
        os.remove(os.path.join(self.cwd, "samplesIDs", "hg38_walk.samplesID.txt"))
        self.ns["_write_dataset_samplesid"]("walkz")  # and removed
        self.assertEqual(self.read(self.sid("walkz")), placeholders)

    def test_a_list_edited_since_it_was_written_is_never_changed(self):
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1", "S2"])
        self.ns["_write_dataset_samplesid"]("walkz")
        mine = _HEADER + "S1\tGBR\tEUR\tmale\nS2\tYRI\tAFR\tfemale\n"
        with open(self.sid("walkz"), "w") as fh:
            fh.write(mine)  # the populations, as the Settings card says to give them
        self.assertIsNone(self.ns["_own_list"]("walkz"))
        self.put("hg38_walk.samplesID.txt", _PANEL)
        self.vcf("walkz", "walkz.chr2.vcf.gz", ["S1", "S2", "S3"])
        self.assertEqual(self.ns["_write_dataset_samplesid"]("walkz"), self.sid("walkz"))
        self.assertEqual(self.read(self.sid("walkz")), mine)

    def test_a_list_of_this_page_is_one_whatever_stops_its_writing(self):
        # recorded as this page's before and after it goes in: stopped between the two,
        # the list there, old or new, is still one this page wrote and may write again
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1"])
        self.ns["_write_dataset_samplesid"]("walkz")
        old = self.read(self.sid("walkz"))
        self.vcf("walkz", "walkz.chr2.vcf.gz", ["S1", "S2"])
        replace, stops = os.replace, []

        def stopped_replace(src, dst):
            if dst == self.sid("walkz"):
                stops.append(dst)
                raise OSError("stopped")
            return replace(src, dst)

        self.ns["os"] = type("os_", (), {"__getattr__": lambda _s, n: getattr(os, n),
                                         "replace": staticmethod(stopped_replace)})()
        with self.assertRaises(OSError):
            self.ns["_write_dataset_samplesid"]("walkz")
        self.ns["os"] = os
        self.assertEqual(stops, [self.sid("walkz")])
        self.assertEqual(self.read(self.sid("walkz")), old)
        self.assertEqual(self.ns["_own_list"]("walkz"), old)
        self.ns["_write_dataset_samplesid"]("walkz")
        self.assertEqual(self.read(self.sid("walkz")), _HEADER
                         + "S1\twalkz\tunknown\tunknown\nS2\twalkz\tunknown\tunknown\n")


class TestRemovingADatasetRemovesItsOwnList(_SamplesIDCase):
    """_delete_targets of 'vcf:<name>': the list this page wrote goes with the dataset,
    which no longer genotypes its samples (it gave other datasets' reports their
    superpopulations, round 3 review F1), with its record and a build's listing."""

    names = _SID_NAMES + ["_delete_targets"]

    def test_own_list_goes_with_the_dataset(self):
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1"])
        self.ns["_write_dataset_samplesid"]("walkz")
        self.ns["_samplesid_listing"]("walkz")
        sdir = os.path.join(self.cwd + "/", "samplesIDs")
        self.assertEqual(self.ns["_delete_targets"]("vcf", "walkz"), [
            os.path.join(self.cwd + "/", "VCFs", "walkz"),
            os.path.join(self.cwd + "/", "VCFs", ".walkz.refgenome"),
            os.path.join(sdir, "walkz.samplesID.txt"),
            os.path.join(sdir, ".walkz.samplesID.sha256"),
            os.path.join(sdir, ".walkz.samplesID.list"),
        ])

    def test_a_list_from_elsewhere_stays(self):
        self.vcf("walkz", "walkz.chr1.vcf.gz", ["S1"])
        self.ns["_write_dataset_samplesid"]("walkz")
        with open(self.sid("walkz"), "w") as fh:
            fh.write(_HEADER + "S1\tGBR\tEUR\tmale\n")  # edited since: yours
        targets = self.ns["_delete_targets"]("vcf", "walkz")
        self.assertNotIn(os.path.join(self.cwd + "/", "samplesIDs", "walkz.samplesID.txt"),
                         targets)
        self.assertIn(os.path.join(self.cwd + "/", "VCFs", "walkz"), targets)


class TestVcfUploadWritesItsSampleList(_SamplesIDCase):
    names = _SID_NAMES + ["_validate_name", "_unusable_in_shell", "_write_vcf_marker",
                          "_vcf_upload_problem", "_finalize_upload"]

    def test_uploaded_vcf_gets_its_marker_and_list(self):
        part = os.path.join(self.cwd, "Settings", "uploads", "walk.chr1.vcf.gz.part")
        _write_vcf(part, ["S1", "S2"])
        got = self.ns["_finalize_upload"]("vcf", part, "walk.chr1.vcf.gz", "walkvcf", "hg38")
        self.assertEqual(got, "walkvcf")
        self.assertEqual(self.read(os.path.join(self.cwd, "VCFs", ".walkvcf.refgenome")), "hg38")
        self.assertEqual(self.read(self.sid("walkvcf")), _HEADER
                         + "S1\twalkvcf\tunknown\tunknown\nS2\twalkvcf\tunknown\tunknown\n")


    def test_upload_keeps_a_list_of_its_header_alone(self):
        with open(self.sid("walkagg"), "w") as fh:
            fh.write(_HEADER)
        part = os.path.join(self.cwd, "Settings", "uploads", "walkagg.vcf.gz.part")
        _write_vcf(part, [])
        self.ns["_finalize_upload"]("vcf", part, "walkagg.vcf.gz", "walkagg", "hg38")
        self.assertEqual(self.read(self.sid("walkagg")), _HEADER)


class TestSettingsBuildPassesSampleList(_SamplesIDCase):
    """build_index, undecorated: the command a Settings build launches."""

    names = _SID_NAMES + ["build_index"]

    def setUp(self):
        self.launched = []
        self.stubs = {
            "launch_settings_job":
                lambda argv, stage: self.launched.append((argv, stage)) or "job",
            "_start": lambda job: (job, False, "started"),
            "ONLINE": False,
            "PreventUpdate": _PreventUpdate,
            "no_update": _NO_UPDATE,
        }
        super().setUp()

    def argv(self, vcf, name=None):
        self.ns["build_index"](1, "hg38", "20bp-NRG-SpCas9", 1, 1, vcf, name)
        return self.launched[-1]

    def base_argv(self):
        return ["build-index-only", "--genome", os.path.join(self.cwd + "/", "Genomes", "hg38"),
                "--pam", os.path.join(self.cwd + "/", "PAMs", "20bp-NRG-SpCas9.txt"),
                "--bDNA", "1", "--bRNA", "1"]

    def test_variant_build_passes_the_datasets_sample_list(self):
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        argv, stage = self.argv("walk")
        listing = os.path.join(self.cwd + "/", "samplesIDs", ".walk.samplesID.list")
        self.assertEqual(stage, "Build variant index")
        self.assertEqual(argv, self.base_argv() + [
            "--vcf", os.path.join(self.cwd + "/", "VCFs", "walk"), "--samplesID", listing])
        self.assertEqual(self.read(listing), "walk.samplesID.txt\n")
        self.assertTrue(os.path.isfile(self.sid("walk")))

    def test_a_list_from_elsewhere_is_passed_unchanged(self):
        mine = _HEADER + "S1\tGBR\tEUR\tmale\nS2\tYRI\tAFR\tfemale\n"
        with open(self.sid("walk"), "w") as fh:
            fh.write(mine)
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        argv, _stage = self.argv("walk")
        self.assertIn("--samplesID", argv)
        self.assertEqual(self.read(self.sid("walk")), mine)

    def test_sites_only_dataset_built_without_one(self):
        self.vcf("agg", "agg.vcf.gz", [])
        argv, _stage = self.argv("agg")
        self.assertEqual(argv, self.base_argv() + [
            "--vcf", os.path.join(self.cwd + "/", "VCFs", "agg")])

    def test_sites_only_dataset_with_a_list_from_elsewhere_built_without_it(self):
        # a Settings build of the mega index (final review round 1, F1): passed its
        # list, the build emitted tiers from it in place of the registry made from the
        # VCFs' INFO allele frequencies, and every variant off-target of a search was gone
        with open(self.sid("hg38_mega"), "w") as fh:
            fh.write(_MEGA_LIST)
        self.vcf("hg38_mega", "mega.chr1.vcf.gz", [])
        self.vcf("hg38_mega", "mega.chr2.vcf.gz", [])
        argv, _stage = self.argv("hg38_mega")
        self.assertEqual(argv, self.base_argv() + [
            "--vcf", os.path.join(self.cwd + "/", "VCFs", "hg38_mega")])
        self.assertEqual(self.read(self.sid("hg38_mega")), _MEGA_LIST)
        self.assertEqual(os.listdir(os.path.join(self.cwd, "samplesIDs")),
                         ["hg38_mega.samplesID.txt"])

    def test_sites_only_dataset_keeps_a_list_of_its_header_alone(self):
        with open(self.sid("agg"), "w") as fh:
            fh.write(_HEADER)
        self.vcf("agg", "agg.vcf.gz", [])
        argv, _stage = self.argv("agg")
        self.assertEqual(argv, self.base_argv() + [
            "--vcf", os.path.join(self.cwd + "/", "VCFs", "agg")])
        self.assertEqual(self.read(self.sid("agg")), _HEADER)

    def test_reference_build_command_unchanged(self):
        self.assertEqual(self.argv("", "my index"),
                         (self.base_argv() + ["--name", "my index"], "Build index"))


# --- the false 'Samples degraded' alert after a build (final review FR#8) ------------
class TestDegradedAlertOnlyForIndexDownloads(unittest.TestCase):
    """The detector matches crisprme_hf's own notes only, and refresh_settings_job runs
    it for a 'Download index' job only: every Settings variant build printed
    '[indel-snp] phased indel genotype store: ...' and got the alert."""

    BUILD = ("[indel-snp] phased indel genotype store: 2 chromosome(s) -> "
             "indel_genotypes_x\nVariant index ready: NRG_2_hg38+x (+ _INDELS).\n")
    PUBLISH = "Published Tier-1 genotype store to lucapinello/crisprme-data:indexes/x\n"
    FETCH_FAILED = ("NOTE: could not fetch Tier-1 genotype store 'genotypes_x.tar.gz' "
                    "(timeout); continuing without it — off-target detection works, "
                    "Samples degraded.\n")
    ABSENT = ("NOTE: index 'NRG_3_hg38+x' has no genotypes_x.tar.gz companion in "
              "repo (registry-only/detection publish); off-target detection works, "
              "per-sample Samples degraded.\n")
    SKIPPED = ("Skipping Tier-1 genotype store for 'NRG_3_hg38+x' (--no-genotypes): "
               "off-target detection works, but per-sample Samples will be degraded "
               "until genotypes_x/ is present.\n")

    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)
        self.ns = _load_functions(
            ["_JOB_ID_RE", "_settings_jobdir", "_genotypes_degraded_warning",
             "refresh_settings_job"],
            current_working_directory=self.cwd,
            _build_dirs_left_warning=lambda jobdir: None,
            _render_all_tables=lambda: "tables",
            _render_storage=lambda: "storage",
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
        )

    def job(self, stage, verbose):
        jobdir = os.path.join(self.cwd, "Settings", _GOOD_ID)
        os.makedirs(jobdir, exist_ok=True)
        with open(os.path.join(jobdir, "log.txt"), "w") as fh:
            fh.write(f"{stage}\tQueued\n{stage}\tStart\n{stage}\tEnd\n")
        with open(os.path.join(jobdir, "log_verbose.txt"), "w") as fh:
            fh.write(verbose)
        return jobdir

    def test_detector_matches_crisprme_hf_notes_only(self):
        for label, text, alert in (
                ("build-index-only's store line", self.BUILD, False),
                ("publish-index's line", self.PUBLISH, False),
                ("companion fetch failed", self.FETCH_FAILED, True),
                ("companion absent", self.ABSENT, True),
                ("--no-genotypes", self.SKIPPED, True),
                ("installed after a failed try",
                 self.FETCH_FAILED + "Installed Tier-1 genotype store genotypes_x/ into D\n",
                 False)):
            with self.subTest(label):
                got = self.ns["_genotypes_degraded_warning"](self.job("Download index", text))
                self.assertIs(got is not None, alert)

    def test_only_an_index_download_is_checked(self):
        for stage, alert in (("Download index", True), ("Build variant index", False),
                             ("Publish index", False), ("Download VCF", False)):
            with self.subTest(stage):
                self.job(stage, self.FETCH_FAILED)
                children = self.ns["refresh_settings_job"](1, _GOOD_ID)[0]
                # done alone, or done with the alert
                self.assertIs(children.args[0] != f"{stage}: done.", alert, children)

    def test_variant_build_shows_a_plain_done(self):
        self.job("Build variant index", self.BUILD)
        children = self.ns["refresh_settings_job"](1, _GOOD_ID)[0]
        self.assertEqual(children.args, ("Build variant index: done.",))


# --- build directories a stopped build left (final review FR#7) ----------------------
class TestBuildDirsLeftWarning(unittest.TestCase):
    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)
        self.ns = _load_functions(
            ["_JOB_ID_RE", "_settings_jobdir", "_dir_size", "_fmt_size",
             "_build_dirs_left_warning", "refresh_settings_job"],
            current_working_directory=self.cwd + "/",
            _genotypes_degraded_warning=lambda jobdir: None,
            _render_failure=lambda stage, jobdir, job_id: ("failure", stage),
            _render_all_tables=lambda: "tables",
            _render_storage=lambda: "storage",
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
        )
        self.left = os.path.join(self.cwd, "genome_library", ".build-index-only-0ln00w08")
        os.makedirs(self.left)
        with open(os.path.join(self.left, "NRG_chr1_1.bin"), "wb") as fh:
            fh.write(b"\0" * 1500)
        self.jobdir = os.path.join(self.cwd, "Settings", _GOOD_ID)
        os.makedirs(self.jobdir)

    def verbose(self, text):
        with open(os.path.join(self.jobdir, "log_verbose.txt"), "w") as fh:
            fh.write(text)

    def leaving(self, path):
        # crisprme.py _reclaim_build_dirs, for an owner on another host/PID namespace
        return (f"Leaving {path}: the build-index-only run that made it (pid 111 on "
                f"c8d74a49b10d, started 2026-10-31T01:02:03Z) ran in another PID "
                f"namespace (another container, or this one before a restart) or on "
                f"another host, so this run cannot tell whether it still runs. If it "
                f"does not, remove it: rm -rf '{path}'\n")

    def test_names_each_directory_left_with_its_size(self):
        self.verbose(self.leaving(self.left) + "Building reference index NRG_2_toyg...\n")
        alert = self.ns["_build_dirs_left_warning"](self.jobdir)
        self.assertIsNotNone(alert)
        code = alert.args[0][1]
        self.assertEqual(code.args, ("genome_library/.build-index-only-0ln00w08 (1.5 KB)",))

    def test_nothing_to_say(self):
        for label, text in (("no Leaving line", "Index built: x\n"),
                            ("reclaimed since", self.leaving(self.left + "_gone"))):
            with self.subTest(label):
                self.verbose(text)
                self.assertIsNone(self.ns["_build_dirs_left_warning"](self.jobdir))
        os.remove(os.path.join(self.jobdir, "log_verbose.txt"))
        self.assertIsNone(self.ns["_build_dirs_left_warning"](self.jobdir))

    def test_done_banner_shows_it(self):
        with open(os.path.join(self.jobdir, "log.txt"), "w") as fh:
            fh.write("Build index\tQueued\nBuild index\tStart\nBuild index\tEnd\n")
        self.verbose(self.leaving(self.left))
        children = self.ns["refresh_settings_job"](1, _GOOD_ID)[0]
        done, alert = children.args[0]
        self.assertEqual(done.args, ("Build index: done.",))
        self.assertIn("unfinished build directories", alert.args[0][0])


class TestFailureBannerShowsBuildDirsLeft(unittest.TestCase):
    def test_failed_banner_includes_the_warning(self):
        ns = _load_functions(["_render_failure"], _error_tail=lambda path: "boom",
                             _build_dirs_left_warning=lambda jobdir: "LEFT", dcc=_StubModule())
        banner = ns["_render_failure"]("Build index", "/x", _GOOD_ID)
        self.assertIn("LEFT", banner.args[0])


# --- a custom annotation BED the page accepts can always be enabled (walk P11) --------
class TestAnnotationChoices(unittest.TestCase):
    def choices(self, genome, beds, genomes=("hg38", "susScr11")):
        ns = _load_functions(
            ["_annotation_for_genome", "_annotation_choices"],
            get_custom_annotations=lambda: [{"label": b, "value": b} for b in beds],
            get_available_genomes=lambda: [{"label": g, "value": g} for g in genomes],
            resolve_builtin_annotation=lambda: None,
            read_enabled_annotations=lambda g: [],
            BUILTIN_INTOGEN_HG38="intogen_drivers.hg38.bed.gz",
        )
        return [o["value"] for o in ns["_annotation_choices"](genome)[0]]

    def test_bed_naming_no_genome_offered_for_every_genome(self):
        beds = ["walk_regions.bed", "track.hg38.bed", "pig.susScr11.bed"]
        self.assertEqual(self.choices("hg38", beds), ["walk_regions.bed", "track.hg38.bed"])
        self.assertEqual(self.choices("susScr11", beds), ["walk_regions.bed", "pig.susScr11.bed"])

    def test_genome_named_bed_still_offered_for_its_genome_only(self):
        self.assertEqual(self.choices("hg38", ["pig.susScr11.bed"]), [])
        # a genome with spaces in its dropdown value (get_available_genomes)
        self.assertEqual(self.choices("hg38 chr22", ["x.hg38_chr22.bed"],
                                      genomes=("hg38 chr22",)), ["x.hg38_chr22.bed"])


class TestAnnotationUploadResetsItsContents(unittest.TestCase):
    """dcc.Upload fires on a change of its contents: the callback empties them, so
    a next file with the same bytes is added too."""

    def test_contents_emptied_after_each_upload(self):
        import base64

        cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, cwd, ignore_errors=True)
        ns = _load_functions(
            ["_validate_name", "_unusable_in_shell", "_fb_refresh", "_add_annotation",
             "add_annotation"],
            current_working_directory=cwd,
            base64=base64,
            validate_annotation_bed=lambda path: None,
            get_available_genomes=lambda: [{"label": "hg38", "value": "hg38"}],
            _render_all_tables=lambda: "tables",
            _render_storage=lambda: "storage",
            ONLINE=False,
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
        )
        bed = "data:text/plain;base64," + base64.b64encode(b"chr1\t1\t2\tx\n").decode()
        for name in ("walk_regions.bed", "walk_regions_copy.bed"):
            out = ns["add_annotation"](bed, name)
            self.assertEqual(len(out), 4)
            self.assertIsNone(out[3])
            self.assertIn(name, out[0].args[0])
            self.assertTrue(os.path.isfile(os.path.join(cwd, "Annotations", name)))
        self.assertIn("so every genome offers it", out[0].args[0])
        out = ns["add_annotation"](bed, "track.hg38.bed")
        self.assertIn("'Manage annotations' for hg38", out[0].args[0])


# --- the 'Started' message (walk P15 (2)) ---------------------------------------------
class TestStartedFeedback(unittest.TestCase):
    def setUp(self):
        self.ns = _load_functions(
            ["_STARTED", "_start", "_shows_started", "_LOCKABLE_BUTTON_IDS",
             "_JOB_FEEDBACK_IDS", "_lock_buttons_while_busy"],
            no_update=_NO_UPDATE,
        )

    def test_started_in_the_neutral_colour(self):
        job, enabled, feedback = self.ns["_start"]("set_ABCDEFGH12")
        self.assertEqual((job, enabled), ("set_ABCDEFGH12", False))
        self.assertEqual(feedback.args, (self.ns["_STARTED"],))
        self.assertEqual(feedback.kwargs["style"]["color"], "#555")

    def test_cleared_once_the_job_has_ended(self):
        started = {"props": {"children": self.ns["_STARTED"], "style": {"color": "#555"}},
                   "type": "Span", "namespace": "dash_html_components"}
        error = "Select an installed genome and PAM."
        n = len(self.ns["_LOCKABLE_BUTTON_IDS"])
        out = self.ns["_lock_buttons_while_busy"]("set_ABCDEFGH12", started, error, None, started)
        self.assertEqual(out, [True] * n + [_NO_UPDATE] * 4)  # running: kept
        out = self.ns["_lock_buttons_while_busy"](None, started, error, None, self.ns["_STARTED"])
        self.assertEqual(out, [False] * n + ["", _NO_UPDATE, _NO_UPDATE, ""])


# --- names CRISPRme takes into shell commands -----------------------------------------
# Settings took the name of a genome (typed, or an upload's file name), a VCF dataset, a
# nuclease and an HPRC sample with any character but a space, '+', '/', '\' and '..', an
# uploaded file's name with any but '..', and a genome download's URL as typed. CRISPRme
# takes such names into shell commands as they are: a search's post-analysis the genome's,
# the dataset's and the PAM file's (named after the nuclease), an index build the genome's
# and the dataset's, a genome download its own and the name its URL ends in. Each name
# below holds what a shell reads there: a command in it would have run, a pattern stood for
# other names, a blank split it, and a '-' or '.' first made it an option or a hidden file.
_SHELL_NAMES = [
    "x$(touch${IFS}MARK)", "x;touch${IFS}MARK;", "x`touch${IFS}MARK`", "x|touch${IFS}MARK",
    "x&&touch${IFS}MARK", "x>MARK", "x<MARK", "x'MARK'", 'x"MARK"', "x#MARK", "x*", "x?",
    "x[ab]", "x\ttouch", "x\ntouch${IFS}MARK", "-x", ".x",
]
# names the catalogs, HPRC and the docs use: they go in as before
_USUAL_NAMES = ["hg38", "susScr11", "T2T-CHM13v2.0", "GRCh38.p14", "hg38_HGDP1kGP",
                "my_pig_vcf", "HG01255", "SpCas9", "e2epanel"]


class _Request:
    """A stand-in for flask.request: one chunk of an upload, as the uploader sends it."""

    def __init__(self, target, file_name, body, dataset="", genome=""):
        self.headers = {"X-Target": target, "X-File-Name": file_name, "X-Dataset": dataset,
                        "X-Genome": genome, "X-Chunk-Index": "0", "X-Total-Chunks": "1"}
        self.body = body

    def get_data(self):
        return self.body


class TestSettingsRefusesNamesAShellReads(unittest.TestCase):
    """The callbacks that add a genome, a VCF dataset or a nuclease, fetch an HPRC sample
    or take an upload, undecorated, on a data folder holding genome hg38: a name a shell
    reads is refused, with nothing written or launched, and the usual names go in."""

    def setUp(self):
        import personal_assembly  # (PostProcess/: a bundle upload's import)

        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)
        for d in ("Genomes/hg38", "VCFs", "PAMs", "samplesIDs", "Settings/uploads"):
            os.makedirs(os.path.join(self.cwd, d))
        self.launched = []
        self.ns = _load_functions(
            ["_validate_name", "_unusable_in_shell", "_genome_upload_problem", "add_genome",
             "_write_vcf_marker", "add_vcf", "_fb_refresh", "add_pam", "_shlex_quote",
             "fetch_from_hprc", "_uploads_dir", "_unsafe_tar_member", "_genome_upload_assembly",
             "_install_uploaded_genome", "_vcf_upload_problem", "_finalize_upload",
             "_MAX_UPLOAD_BYTES", "_fmt_size", "_settings_upload_chunk"],
            current_working_directory=self.cwd,
            __file__=_SETTINGS_PAGE,
            sys=__import__("sys"),
            gzip=gzip,
            tarfile=tarfile,
            shutil=shutil,
            tempfile=tempfile,
            genome_in_use=genome_in_use,
            personal_assembly=personal_assembly,
            ONLINE=False,
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
            launch_settings_job=lambda argv, stage: self.launched.append(argv) or "job",
            launch_settings_job_raw=lambda cmd, stage: self.launched.append(cmd) or "job",
            _start=lambda job: (job, False, "started"),
            _preflight_disk=lambda need: None,
            _catalog_size=lambda component, name: 0,
            _write_dataset_samplesid=lambda dataset: None,
            _render_all_tables=lambda: "tables",
            _render_storage=lambda: "storage",
            validate_annotation_bed=lambda path: None,
            LIFTOVER_DIR="LiftoverFiles",
        )

    def refused(self, call, said):
        """``call()`` refused a name, saying why, with nothing written or launched."""
        before = _tree(self.cwd)
        said(call())
        self.assertEqual(_tree(self.cwd), before)
        self.assertEqual(self.launched, [])

    def assert_says_shell(self, message):
        self.assertIn("CRISPRme takes it into shell commands", str(message))

    def test_a_name_a_shell_reads_is_refused_and_the_usual_go_in(self):
        for name in _SHELL_NAMES:
            with self.subTest(name=name):
                self.assert_says_shell(self.ns["_validate_name"](name))
        for name in _USUAL_NAMES + [" hg38 "]:  # (the blanks around a name typed go)
            with self.subTest(name=name):
                self.assertIsNone(self.ns["_validate_name"](name))

    def test_a_genome_download(self):
        url = "https://example.org/x.fa.gz"
        for name in _SHELL_NAMES:
            with self.subTest(name=name):
                self.refused(lambda: self.ns["add_genome"](1, name, "url", url, None),
                             lambda out: self.assert_says_shell(out[2]))
                # the download is saved under the name its URL ends in
                self.refused(lambda: self.ns["add_genome"](
                    1, "susScr11", "url", f"https://example.org/{name}.fa.gz", None),
                    lambda out: self.assert_says_shell(out[2]))
        self.assertEqual(self.ns["add_genome"](1, "T2T-CHM13v2.0", "url", url, None),
                         ("job", False, "started"))
        self.assertEqual(self.launched, [["download", "--what", "genome", "--ref",
                                          "T2T-CHM13v2.0", "--source", "url", "--url", url]])

    def test_a_vcf_dataset_registered_from_a_folder(self):
        folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, folder, ignore_errors=True)
        for name in _SHELL_NAMES:
            with self.subTest(name=name):
                self.refused(lambda: self.ns["add_vcf"](1, name, "path", folder, None, "hg38"),
                             lambda out: self.assert_says_shell(out[2]))
        out = self.ns["add_vcf"](1, "hg38_HGDP1kGP", "path", folder, None, "hg38")
        self.assertEqual(out[2], "Registered VCF dataset 'hg38_HGDP1kGP' for hg38.")
        self.assertEqual(os.readlink(os.path.join(self.cwd, "VCFs", "hg38_HGDP1kGP")), folder)

    def test_a_nuclease(self):
        for name in _SHELL_NAMES:
            with self.subTest(name=name):
                self.refused(lambda: self.ns["add_pam"](1, name, 20, "NGG", "3"),
                             lambda out: self.assert_says_shell(out[0]))
        self.assertEqual(self.ns["add_pam"](1, "SpCas9", 20, "NGG", "3")[1:],
                         ("tables", "storage"))
        self.assertEqual(os.listdir(os.path.join(self.cwd, "PAMs")), ["20bp-NGG-SpCas9.txt"])

    def test_an_hprc_sample(self):
        for name in _SHELL_NAMES:
            with self.subTest(name=name):
                self.refused(lambda: self.ns["fetch_from_hprc"](1, name),
                             lambda out: self.assert_says_shell(out[2]))
        self.assertEqual(self.ns["fetch_from_hprc"](1, "HG01255"), ("job", False, "started"))
        self.assertEqual(len(self.launched), 1)
        self.assertIn(" HG01255 --data-dir ", self.launched[0])

    def upload(self, *args, **kwargs):
        self.ns["request"] = _Request(*args, **kwargs)
        return self.ns["_settings_upload_chunk"]()

    def test_an_upload(self):
        vcf = io.BytesIO()
        with gzip.GzipFile(fileobj=vcf, mode="wb", mtime=0) as fh:
            fh.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for name in _SHELL_NAMES:
            with self.subTest(name=name):
                said = lambda out: (self.assertEqual(out[1], 400), self.assert_says_shell(out[0]))
                # a genome, under its file's name; a VCF dataset; a VCF file, the part of
                # whose name that holds 'chr' a variant search names its indels after
                self.refused(lambda: self.upload("genome", f"{name}.fa", _FA1), said)
                self.refused(lambda: self.upload("vcf", "walk.chr1.vcf.gz", vcf.getvalue(),
                                                 dataset=name, genome="hg38"), said)
                self.refused(lambda: self.upload("vcf", f"{name}chr1.vcf.gz", vcf.getvalue(),
                                                 dataset="walk", genome="hg38"), said)
        self.assertEqual(self.upload("genome", "susScr11.fa", _FA1 + _FA2),
                         ("complete:susScr11", 200))
        self.assertEqual(sorted(os.listdir(os.path.join(self.cwd, "Genomes", "susScr11"))),
                         ["chr1.fa", "chr2.fa"])
        self.assertEqual(self.upload("vcf", "walk.chr1.vcf.gz", vcf.getvalue(), dataset="walk",
                                     genome="hg38"), ("complete:walk", 200))
        self.assertEqual(os.listdir(os.path.join(self.cwd, "VCFs", "walk")), ["walk.chr1.vcf.gz"])
        # a personal-assembly bundle is named after its archive's top folder: the file's
        # name only says how it is packed, and goes in with a blank too (as a second
        # download of it is named)
        meta = {"format": "crisprme-personal-assembly", "format_version": 1,
                "individual": "HG01255", "haplotypes": {
                    hap: {"genome_dir": f"{hap}/genome", "chain": f"{hap}/x.chain",
                          "chromalias": f"{hap}/x.chromAlias.txt"}
                    for hap in ("paternal", "maternal")}}
        part = os.path.join(self.cwd, "bundle.tar.gz")
        entries = [("HG01255/metadata.json", "file", json.dumps(meta).encode())]
        for hap in ("paternal", "maternal"):
            entries += [(f"HG01255/{hap}/genome/chr1.fa", "file", _FA1),
                        (f"HG01255/{hap}/x.chain", "file", b"x\n"),
                        (f"HG01255/{hap}/x.chromAlias.txt", "file", b"x\n")]
        _write_tar_gz(part, entries)
        with open(part, "rb") as fh:
            body = fh.read()
        os.remove(part)
        self.assertEqual(self.upload("assembly-bundle", "HG01255 (1).tar.gz", body),
                         ("complete:HG01255", 200))
        self.assertTrue(os.path.isfile(os.path.join(self.cwd, "Assemblies", "HG01255",
                                                    "metadata.json")))


# --- what a Settings dropdown offers --------------------------------------------------
# Dash does not check a dropdown's value against its options, and a request can send a
# callback any value, which names the folders and files the callback removes, writes,
# archives or reads: outside the data folder too, unless it takes only what it offers.
def _listing(root, folder, keep=os.path.isdir, value=lambda name: name):
    """A stand-in for a pages_utils listing (get_available_genomes, ...) of the entries of
    <root>/<folder> that ``keep`` takes, hidden ones left out, as its dropdown offers them."""
    def listed():
        top = os.path.join(root, folder)
        return [{"label": value(n), "value": value(n)} for n in sorted(os.listdir(top))
                if not n.startswith(".") and keep(os.path.join(top, n))]
    return listed


class TestRemoveTakesAnItemTheCardLists(unittest.TestCase):
    """ask_delete and do_delete, undecorated, on a data folder (data/) beside a folder of
    someone else's (outside/). The item's name was only basename-guarded, which leaves
    '..', '.' and '' as they are: '<kind>:..' removed the whole data folder, every job's
    results with it, and, through a folder of it that is a link, the folder its target is
    in; '<kind>:.' and '<kind>:' removed every item of the kind. Only an item the card lists
    is measured and removed; of another, nothing is, and Delete says it is not found."""

    KINDS = ("genome", "index", "vcf", "annotation", "pam", "assembly")

    def setUp(self):
        import assembly_reconcile
        import personal_assembly

        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.cwd = os.path.join(self.root, "data")
        for path, text in (("Genomes/hg38/chr1.fa", ">chr1\nACGTACGT\n"),
                           ("Genomes/delme/chr1.fa", ">chr1\nACGT\n"),
                           ("genome_library/NRG_2_hg38/NRG_chr1_1.bin", "x"),
                           ("VCFs/walk/walk.chr1.vcf.gz", "x"),
                           ("Annotations/a.bed", "chr1\t1\t2\tx\n"),
                           ("PAMs/20bp-NRG-SpCas9.txt", "NNNNNNNNNNNNNNNNNNNNNRG 3\n"),
                           ("Results/JOB1/results.txt", "a job's results\n"),
                           ("../outside/keep.txt", "someone else's\n")):
            path = os.path.join(self.cwd, path)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write(text)
        os.makedirs(os.path.join(self.cwd, "samplesIDs"))
        self.ns = _load_functions(
            ["_deletable_options", "_deletable", "_delete_targets", "_delete_dependents",
             "_delete_summary", "_dir_size", "_fmt_size", "_own_lists_record", "_own_list",
             "ask_delete", "do_delete"],
            current_working_directory=self.cwd + "/",
            shutil=shutil,
            hashlib=hashlib,
            personal_assembly=personal_assembly,
            assembly_reconcile=assembly_reconcile,
            get_available_genomes=_listing(self.cwd, "Genomes"),
            get_available_indexes=_listing(
                self.cwd, "genome_library",
                keep=lambda p: os.path.isdir(p) and not p.endswith("_INDELS")),
            get_all_vcf_datasets=_listing(self.cwd, "VCFs"),
            get_custom_annotations=_listing(self.cwd, "Annotations", keep=os.path.isfile),
            get_available_PAM=_listing(self.cwd, "PAMs", keep=os.path.isfile,
                                       value=lambda name: name[:-len(".txt")]),
            installed_assemblies=lambda: [],
            index_build_pam=lambda index_dir: (None, 0, 0),
            _render_all_tables=lambda: "tables",
            _render_storage=lambda: "storage",
            GENOMES_DIR="Genomes",
            LIFTOVER_DIR="LiftoverFiles",
            ONLINE=False,
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
        )
        self.measured = []  # what ask_delete measures, as the size the dialog says it frees
        dir_size = self.ns["_dir_size"]
        self.ns["_dir_size"] = lambda path: self.measured.append(path) or dir_size(path)

    def crafted(self):
        outside = os.path.join(self.root, "outside")
        return ([f"{kind}:{name}" for kind in self.KINDS for name in ("..", ".", "")]
                + ["genome:../../outside", f"genome:{outside}", "genome:hg38/../..",
                   "index:NRG_2_hg38/..", "vcf:walk/../../../outside", "nosuchkind:.."])

    def test_an_item_the_card_does_not_list_is_neither_measured_nor_removed(self):
        before = _tree(self.root)
        for item in self.crafted():
            with self.subTest(item=item):
                self.assertEqual(self.ns["ask_delete"](1, item),
                                 (True, "Delete this item? This cannot be undone."))
                said = self.ns["do_delete"](1, item)[0]
                self.assertEqual(said.args, (f"{item} not found (already removed?).",))
                self.assertEqual(_tree(self.root), before)
        self.assertEqual(self.measured, [])

    def test_through_a_folder_of_the_data_folder_that_is_a_link(self):
        # genome_library kept on other storage, beside someone else's data: '<link>/..' is
        # the folder the link's target is in
        shared = os.path.join(self.root, "shared")
        os.makedirs(os.path.join(shared, "other_lab_data"))
        with open(os.path.join(shared, "other_lab_data", "keep.txt"), "w") as fh:
            fh.write("someone else's\n")
        os.rename(os.path.join(self.cwd, "genome_library"),
                  os.path.join(shared, "genome_library"))
        os.symlink(os.path.join(shared, "genome_library"),
                   os.path.join(self.cwd, "genome_library"))
        before = _tree(self.root)
        for item in self.crafted():
            with self.subTest(item=item):
                self.ns["do_delete"](1, item)
                self.assertEqual(_tree(self.root), before)
        # an index of it goes as before, and only it
        said = self.ns["do_delete"](1, "index:NRG_2_hg38")[0]
        self.assertEqual(said.args, ("Deleted index:NRG_2_hg38 (freed 1 B). You can download "
                                     "or rebuild it later.",))
        self.assertEqual(os.listdir(os.path.join(shared, "genome_library")), [])
        self.assertTrue(os.path.isfile(os.path.join(shared, "other_lab_data", "keep.txt")))

    def test_an_item_the_card_lists_is_removed_as_before(self):
        self.assertTrue(self.ns["ask_delete"](1, "genome:delme")[1].startswith(
            "Delete reference genome “delme”?\n\nThis frees about 11 B"))
        for item, path in (("genome:delme", "Genomes/delme"),
                           ("index:NRG_2_hg38", "genome_library/NRG_2_hg38"),
                           ("vcf:walk", "VCFs/walk"),
                           ("annotation:a.bed", "Annotations/a.bed"),
                           ("pam:20bp-NRG-SpCas9", "PAMs/20bp-NRG-SpCas9.txt")):
            with self.subTest(item=item):
                before = _tree(self.root)
                said = self.ns["do_delete"](1, item)[0]
                self.assertTrue(said.args[0].startswith(f"Deleted {item} (freed "), said)
                gone = {p for p in before if p == os.path.join("data", path)
                        or p.startswith(os.path.join("data", path) + os.sep)}
                self.assertTrue(gone)
                self.assertEqual(set(_tree(self.root)), set(before) - gone)
        # deleted since the card was drawn: not found
        said = self.ns["do_delete"](1, "genome:delme")[0]
        self.assertEqual(said.args, ("genome:delme not found (already removed?).",))

    def test_a_name_that_is_not_an_entry_of_its_folder_is_refused(self):
        for kind in self.KINDS + ("nosuchkind",):
            for name in ("", ".", "..", "a/b", "../x", "/x", "hg38/.."):
                with self.subTest(kind=kind, name=name), self.assertRaises(ValueError):
                    self.ns["_delete_targets"](kind, name)


class TestVcfUploadIntoARegisteredFolder(unittest.TestCase):
    """The upload route, add_vcf (which registers a folder of the server as a VCF dataset,
    linking VCFs/<name> to it) and _finalize_upload, undecorated. An upload into a
    registered dataset went through the link into that folder, under any file name, new
    or replacing a file there: a request could write any file into any folder the server
    can write. A registered folder is read where it is, and nothing is uploaded into it,
    whenever it was registered; an upload takes a VCF dataset's own files only."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.cwd = os.path.join(self.root, "data")
        for d in ("Genomes/hg38", "VCFs", "samplesIDs", "Settings/uploads"):
            os.makedirs(os.path.join(self.cwd, d))
        self.server = os.path.join(self.root, "server_folder")  # a folder of the server's
        os.makedirs(self.server)
        for name in ("existing.vcf.gz", "notes.txt"):
            with open(os.path.join(self.server, name), "w") as fh:
                fh.write("ORIGINAL\n")
        self.ns = _load_functions(
            ["_validate_name", "_unusable_in_shell", "_write_vcf_marker", "add_vcf",
             "_uploads_dir", "_vcf_upload_problem", "_finalize_upload", "_MAX_UPLOAD_BYTES",
             "_fmt_size", "_settings_upload_chunk"],
            current_working_directory=self.cwd,
            shutil=shutil,
            ONLINE=False,
            PreventUpdate=_PreventUpdate,
            no_update=_NO_UPDATE,
            _write_dataset_samplesid=lambda dataset: None,
            _render_all_tables=lambda: "tables",
            _render_storage=lambda: "storage",
        )

    def upload(self, name, body, dataset, index=0, total=1):
        """One chunk of an upload, as the uploader sends it: the route's (text, status)."""
        self.ns["request"] = _Request("vcf", name, body, dataset=dataset, genome="hg38")
        self.ns["request"].headers.update({"X-Chunk-Index": str(index),
                                           "X-Total-Chunks": str(total)})
        return self.ns["_settings_upload_chunk"]()

    def register(self, dataset):
        said = self.ns["add_vcf"](1, dataset, "path", self.server, None, "hg38")[2]
        self.assertEqual(said, f"Registered VCF dataset '{dataset}' for hg38.")

    def assert_refused(self, answer, status=400):
        self.assertEqual(answer[1], status, answer)
        self.assertIn(" was registered from a folder on the server: its files are read "
                      "there, and nothing is uploaded into it.", answer[0])

    def test_nothing_is_uploaded_into_a_registered_folder(self):
        self.register("regdir")
        before = _tree(self.server, meta=True)
        for name in ("new.chr1.vcf.gz", "existing.vcf.gz", "existing.vcf.gz.tbi"):
            with self.subTest(name=name):
                self.assert_refused(self.upload(name, b"UPLOADED\n", "regdir"))
        self.assertEqual(_tree(self.server, meta=True), before)
        self.assertEqual(os.listdir(os.path.join(self.cwd, "Settings", "uploads")), [])
        # registered while an upload into a dataset of its name was on its way, and named
        # by a later chunk of another upload (each chunk sends its own names)
        self.assertEqual(self.upload("x.chr1.vcf.gz", b"PART1\n", "later", 0, 2),
                         ("chunk received", 200))
        self.register("later")
        self.assertEqual(self.upload("y.chr1.vcf.gz", b"PART1\n", "mine", 0, 2),
                         ("chunk received", 200))
        before = _tree(self.server, meta=True)
        for name, dataset in (("x.chr1.vcf.gz", "later"), ("y.chr1.vcf.gz", "regdir")):
            with self.subTest(dataset=dataset):
                text, status = self.upload(name, b"PART2\n", dataset, 1, 2)
                self.assert_refused((text, status), 500)
        self.assertEqual(_tree(self.server, meta=True), before)
        self.assertEqual(os.listdir(os.path.join(self.cwd, "Settings", "uploads")), [])
        self.assertFalse(os.path.exists(os.path.join(self.cwd, "VCFs", "mine")))

    def test_only_a_vcf_datasets_own_files_are_uploaded(self):
        for name in ("sitecustomize.py", "notes.txt", "walk.chr1.vcf", "walk.chr1.vcf.gz.bak",
                     "walk.tbi"):
            with self.subTest(name=name):
                text, status = self.upload(name, b"x\n", "walk")
                self.assertEqual(status, 400)
                self.assertIn("a VCF dataset takes bgzip-compressed VCF files (.vcf.gz) and "
                              "their index (.vcf.gz.tbi or .vcf.gz.csi)", text)
        self.assertEqual(os.listdir(os.path.join(self.cwd, "VCFs")), [])
        # a dataset's own files go in as before
        for name in ("walk.chr1.vcf.gz", "walk.chr1.vcf.gz.tbi", "walk.chr2.vcf.gz",
                     "walk.chr2.vcf.gz.csi"):
            with self.subTest(name=name):
                self.assertEqual(self.upload(name, b"x\n", "walk"), ("complete:walk", 200))
        self.assertEqual(sorted(os.listdir(os.path.join(self.cwd, "VCFs", "walk"))),
                         ["walk.chr1.vcf.gz", "walk.chr1.vcf.gz.tbi", "walk.chr2.vcf.gz",
                          "walk.chr2.vcf.gz.csi"])


if __name__ == "__main__":
    unittest.main()
