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
            ["_validate_name", "_unsafe_tar_member", "_genome_upload_assembly",
             "_built_on_genome", "_genome_upload_problem", "_install_uploaded_genome",
             "_finalize_upload"],
            current_working_directory=cwd,
            gzip=gzip,
            tarfile=tarfile,
            shutil=shutil,
            tempfile=tempfile,
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

    def test_unusable_upload_refused_installing_nothing(self):
        cases = [
            ("no FASTA in the archive", "asm.tar.gz", [("asm/README", "file", b"x\n")], None),
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
_SID_NAMES = ["SAMPLESID_HEADER", "_vcf_header_samples", "_placeholder_samples",
              "_listed_samples", "_write_dataset_samplesid", "_samplesid_listing"]
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
                                  gzip=gzip, shutil=shutil, threading=threading,
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

    def test_two_threads_writing_the_list_at_once_each_use_their_own_file(self):
        # the web server runs each request in a thread of one process, so two uploads
        # into one dataset that finish together write its list at once: through one
        # temporary file, one truncated the other's and its rename found nothing there
        # (round 2 review, R2-2). Each thread holds its temporary file open until the
        # other has opened its own.
        self.vcf("walk", "walk.chr1.vcf.gz", ["S1", "S2"])
        both_open = threading.Barrier(2, timeout=10)
        opened, errors = [], []

        def open_(path, mode="r", *args, **kwargs):
            fh = open(path, mode, *args, **kwargs)
            if "w" in mode and os.path.basename(path).startswith(".walk.samplesID.txt."):
                opened.append(path)
                both_open.wait()
            return fh

        def write():
            try:
                self.ns["_write_dataset_samplesid"]("walk")
            except Exception as exc:  # noqa: BLE001 - reported below
                errors.append(exc)

        self.ns["open"] = open_
        threads = [threading.Thread(target=write) for _ in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(30)
        self.assertEqual(errors, [])
        self.assertEqual(len(set(opened)), 2, opened)
        self.assertEqual(self.read(self.sid("walk")), _HEADER
                         + "S1\twalk\tunknown\tunknown\nS2\twalk\tunknown\tunknown\n")
        self.assertEqual(os.listdir(os.path.join(self.cwd, "samplesIDs")),
                         ["walk.samplesID.txt"])


class TestVcfUploadWritesItsSampleList(_SamplesIDCase):
    names = _SID_NAMES + ["_validate_name", "_write_vcf_marker", "_finalize_upload"]

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
            ["_validate_name", "_fb_refresh", "_add_annotation", "add_annotation"],
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


if __name__ == "__main__":
    unittest.main()
