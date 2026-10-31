"""Tests for the configurable concurrent-search-jobs setting (app_settings.py).

Default stays a conservative 2; env var overrides a persisted value; values are
clamped to [1, MAX_JOBS_CAP] so the pool is never oversubscribed. app.py, loaded
the way the web app is really started, must size its search pool from it, and
the Settings page must say when the env var overrides the value saved there."""

import ast
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock

import app_settings as a

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Loads app.py the way the web app really loads it. `crisprme.py web-interface`
# runs `subprocess.run([sys.executable, index.py])`: a fresh interpreter whose
# sys.path holds only the install root (index.py's directory, no PostProcess/),
# whose cwd is the user's working directory, and whose index.py imports app.py
# first, before pages/pages_utils.py puts PostProcess/ on sys.path. -I drops
# PYTHONPATH, user site-packages and the cwd entry; the root is then put back as
# sys.path[0]. The web stack (not in the light CI env) plays no part in the
# import-path logic, so it is faked through sys.modules.
_APP_IMPORT_DRIVER = textwrap.dedent('''
    import json, sys, types

    class _Stub:
        """Stands in for any dash/flask object: accepts any call, attribute or item."""
        def __init__(self, *args, **kwargs):
            pass
        def __call__(self, *args, **kwargs):
            return _Stub()
        def __getattr__(self, name):
            if name.startswith("__"):
                raise AttributeError(name)
            value = _Stub()
            setattr(self, name, value)
            return value
        def __setitem__(self, key, value):
            pass

    def _module_getattr(name):
        if name.startswith("__"):
            raise AttributeError(name)
        return _Stub()

    for _name in ("flask_caching", "dash_bootstrap_components", "dash", "flask"):
        _mod = types.ModuleType(_name)
        _mod.__getattr__ = _module_getattr
        sys.modules[_name] = _mod

    sys.path.insert(0, sys.argv[1])
    import app

    print(json.dumps({"max_workers": app.pool_executor._max_workers}))
''')


class TestMaxConcurrentJobs(unittest.TestCase):
    def setUp(self):
        os.environ.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)

    def tearDown(self):
        os.environ.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)

    def test_default_is_two(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(a.get_max_concurrent_jobs(d), 2)

    def test_set_and_read_back(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(a.set_max_concurrent_jobs(4, d), 4)
            self.assertEqual(a.get_max_concurrent_jobs(d), 4)

    def test_clamped_to_cap_and_floor(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(a.set_max_concurrent_jobs(999, d), a.MAX_JOBS_CAP)
            self.assertEqual(a.set_max_concurrent_jobs(0, d), 1)

    def test_env_overrides_file(self):
        with tempfile.TemporaryDirectory() as d:
            a.set_max_concurrent_jobs(3, d)
            os.environ["CRISPRME_MAX_CONCURRENT_JOBS"] = "8"
            self.assertEqual(a.get_max_concurrent_jobs(d), 8)

    def test_corrupt_config_falls_back_to_default(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, a.CONFIG_FILE), "w") as f:
                f.write("{not valid json")
            self.assertEqual(a.get_max_concurrent_jobs(d), 2)

    def test_env_override_only_for_a_whole_number(self):
        # ASCII digits only: int() rejects some other str.isdigit() characters,
        # such as superscript two and circled one, and env_override() raised
        # on them.
        self.assertIsNone(a.env_override())
        cases = (
            ("8", 8), ("5", 5), (" 3 ", 3), (" 5 ", 5), ("99", 16), ("0", 1),
            ('"5"', None), ("", None), ("\u00b2", None), ("\u2460", None),
            ("\u0661\u0662", None),  # Arabic-Indic 12: int() takes it, not ASCII
        )
        for raw, pinned in cases:
            with self.subTest(raw=raw):
                os.environ["CRISPRME_MAX_CONCURRENT_JOBS"] = raw
                self.assertEqual(a.env_override(), pinned)

    def test_env_override_never_raises_on_too_many_digits(self):
        # int() converts at most sys.get_int_max_str_digits() (4300) digits.
        os.environ["CRISPRME_MAX_CONCURRENT_JOBS"] = "9" * 5000
        self.assertIsNone(a.env_override())

    def test_ignored_env_value_falls_back_to_the_saved_value(self):
        with tempfile.TemporaryDirectory() as d:
            a.set_max_concurrent_jobs(4, d)
            for raw in ('"5"', "\u00b2", "\u2460"):
                with self.subTest(raw=raw):
                    os.environ["CRISPRME_MAX_CONCURRENT_JOBS"] = raw
                    self.assertEqual(a.get_max_concurrent_jobs(d), 4)


class TestAppSearchPoolSize(unittest.TestCase):
    """app.py used to import app_settings with only the install root on sys.path,
    so under `crisprme.py web-interface` the import always failed and a bare
    except silently pinned the search pool to 2, whatever the env var or the
    Settings page said."""

    def _import_app(
        self, root=_REPO_ROOT, env_jobs=None, saved_jobs=None, io_encoding=None
    ):
        """Import <root>/app.py like the web app does; return (search pool
        max_workers, everything app.py printed on stdout). ``io_encoding`` is the
        child's PYTHONIOENCODING, which -I ignores: the child then runs with -s -P
        (the rest of -I) and without PYTHONPATH."""
        env = dict(os.environ)
        env.pop("CRISPRME_CONFIG_DIR", None)
        env.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)
        if env_jobs is not None:
            env["CRISPRME_MAX_CONCURRENT_JOBS"] = str(env_jobs)
        flags = ["-I"]
        if io_encoding is not None:
            env["PYTHONIOENCODING"] = io_encoding
            env.pop("PYTHONPATH", None)
            flags = ["-s", "-P"]
        with tempfile.TemporaryDirectory() as workdir:  # the user's cwd, not the repo
            if saved_jobs is not None:  # what the Settings page persists
                a.set_max_concurrent_jobs(saved_jobs, workdir)
            res = subprocess.run(
                [sys.executable, *flags, "-c", _APP_IMPORT_DRIVER, root],
                cwd=workdir, env=env, capture_output=True, text=True, timeout=120,
            )
        self.assertEqual(res.returncode, 0, res.stderr)
        lines = res.stdout.strip().splitlines()
        return json.loads(lines[-1])["max_workers"], "\n".join(lines[:-1])

    def test_env_var_sizes_the_search_pool(self):
        workers, out = self._import_app(env_jobs=5)
        self.assertEqual(workers, 5)
        self.assertEqual(out, "")

    def test_settings_page_value_sizes_the_search_pool(self):
        workers, out = self._import_app(saved_jobs=4)
        self.assertEqual(workers, 4)
        self.assertEqual(out, "")

    def test_default_stays_two(self):
        workers, out = self._import_app()
        self.assertEqual(workers, 2)
        self.assertEqual(out, "")

    def test_missing_app_settings_falls_back_to_two_and_says_so(self):
        with tempfile.TemporaryDirectory() as root:  # an install with no PostProcess/
            shutil.copy(os.path.join(_REPO_ROOT, "app.py"), root)
            workers, out = self._import_app(root=root, env_jobs=5)
        self.assertEqual(workers, 2)
        self.assertEqual(len(out.splitlines()), 1)
        self.assertIn("app_settings", out)
        self.assertRegex(out, r"\b2\b")

    # An env value app_settings cannot use as given used to change the pool size
    # without a word: ignored if it is not a whole number, clamped into [1, 16].

    def test_env_value_not_a_whole_number_is_ignored_and_says_so(self):
        # `docker run --env-file` keeps quotes: a CRISPRME_MAX_CONCURRENT_JOBS="5"
        # line arrives as the value "5" with its quotes.
        for raw in ('"5"', "4.0", "-1", "five"):
            with self.subTest(raw=raw):
                workers, out = self._import_app(env_jobs=raw)
                self.assertEqual(workers, 2)
                self.assertEqual(len(out.splitlines()), 1, out)
                self.assertIn(f"CRISPRME_MAX_CONCURRENT_JOBS={raw!r}", out)
                self.assertRegex(out, r"\busing 2\b")

    def test_ignored_env_value_names_the_saved_value_used_instead(self):
        # Superscript two and circled one pass str.isdigit() but not int():
        # app_settings raised on them, so the pool fell back to 2 and the saved
        # value went unused.
        for raw in ('"5"', "\u00b2", "\u2460"):
            with self.subTest(raw=raw):
                workers, out = self._import_app(env_jobs=raw, saved_jobs=3)
                self.assertEqual(workers, 3)
                self.assertEqual(len(out.splitlines()), 1, out)
                self.assertRegex(out, r"\busing 3\b")

    def test_env_value_too_long_for_int_is_ignored_and_says_so(self):
        workers, out = self._import_app(env_jobs="9" * 5000)
        self.assertEqual(workers, 2)
        self.assertEqual(len(out.splitlines()), 1, out[:200])
        self.assertRegex(out, r"\busing 2\b")

    def test_out_of_range_env_value_is_clamped_and_says_so(self):
        for raw, used in (("99", 16), ("0", 1)):
            with self.subTest(raw=raw):
                workers, out = self._import_app(env_jobs=raw)
                self.assertEqual(workers, used)
                self.assertEqual(len(out.splitlines()), 1, out)
                self.assertIn(f"CRISPRME_MAX_CONCURRENT_JOBS={raw!r}", out)
                self.assertRegex(out, rf"\busing {used}\b")

    # The warnings run at import: one that stdout cannot encode used to raise
    # there and stop the web server, e.g. a latin-1 locale and a value pasted
    # with curly quotes. They now print the value and the error ASCII-escaped.

    def test_non_ascii_env_value_is_reported_on_an_ascii_stdout(self):
        for raw in ("\u201c5\u201d", "\u00b2", "\u2460"):
            with self.subTest(raw=raw):
                workers, out = self._import_app(env_jobs=raw, io_encoding="ascii")
                self.assertEqual(workers, 2)
                self.assertEqual(len(out.splitlines()), 1, out)
                self.assertIn(f"CRISPRME_MAX_CONCURRENT_JOBS={raw!a}", out)
                self.assertRegex(out, r"\busing 2\b")

    def test_unreadable_setting_is_reported_on_an_ascii_stdout(self):
        with tempfile.TemporaryDirectory() as root:  # app_settings fails, non-ASCII
            shutil.copy(os.path.join(_REPO_ROOT, "app.py"), root)
            os.mkdir(os.path.join(root, "PostProcess"))
            with open(os.path.join(root, "PostProcess", "app_settings.py"), "w") as fh:
                fh.write('raise OSError("cannot read /srv/caf\\u00e9")\n')
            workers, out = self._import_app(root=root, io_encoding="ascii")
        self.assertEqual(workers, 2)
        self.assertEqual(len(out.splitlines()), 1, out)
        self.assertIn("app_settings", out)
        self.assertIn("caf\\xe9", out)

    def test_env_value_used_as_given_prints_nothing(self):
        # Whole numbers 1-16 are used as given, blanks and a leading zero
        # included; an empty value counts as unset.
        for raw, used in (("1", 1), ("16", 16), (" 7 ", 7), ("05", 5), ("", 2)):
            with self.subTest(raw=raw):
                workers, out = self._import_app(env_jobs=raw)
                self.assertEqual(workers, used)
                self.assertEqual(out, "")


_SETTINGS_PAGE = os.path.join(_REPO_ROOT, "pages", "settings_page.py")


class _Component:
    """Stands in for a dash html/dcc component: keeps its children and props."""

    def __init__(self, *args, **kwargs):
        self.args, self.kwargs = args, kwargs


class _Components:
    """Stands in for dash's html / dcc modules: every name is a _Component."""

    def __getattr__(self, _name):
        return _Component


def _load_concurrency_ui():
    """The Settings page's "Concurrent search jobs" card and its Save callback,
    extracted from pages/settings_page.py (which needs a live Dash app to import).
    Returns (render_card, save): render_card() runs the card's statements from
    settings_page() and returns the (title, blurb, body) they pass to _add_card;
    save is _max_jobs_save itself."""
    src = open(_SETTINGS_PAGE).read()
    funcs = {n.name: n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef)}
    page = funcs["settings_page"].body
    card = next(
        i for i, node in enumerate(page)
        if isinstance(node, ast.Assign)
        and [ast.unparse(t) for t in node.targets] == ["concurrency_card"]
    )
    if not isinstance(page[card - 1], ast.Try):  # the try: import app_settings ...
        raise AssertionError("settings_page(): concurrency_card lost its try block")
    save = funcs["_max_jobs_save"]
    save.decorator_list = []  # @app.callback needs the live Dash app
    ns = {
        "html": _Components(), "dcc": _Components(), "ONLINE": False,
        "PreventUpdate": type("PreventUpdate", (Exception,), {}),
        "_add_card": lambda title, blurb, body: (title, blurb, body),
    }
    save_stmts = ast.Module(body=[save], type_ignores=[])
    exec(compile(save_stmts, _SETTINGS_PAGE, "exec"), ns)  # noqa: S102 - own source
    card_stmts = ast.Module(body=page[card - 1:card + 1], type_ignores=[])
    card_code = compile(card_stmts, _SETTINGS_PAGE, "exec")

    def render_card():
        exec(card_code, ns)  # noqa: S102 - trusted own source
        return ns["concurrency_card"]

    return render_card, ns["_max_jobs_save"]


class TestSettingsPageEnvOverride(unittest.TestCase):
    """CRISPRME_MAX_CONCURRENT_JOBS overrides the value saved from the Settings
    page, which used to say nothing about it: with the env var at 5, saving 3
    answered "Restart the app to apply." though every restart kept 5."""

    @classmethod
    def setUpClass(cls):
        render_card, save = _load_concurrency_ui()
        cls.render_card, cls.save = staticmethod(render_card), staticmethod(save)

    def _settings(self, env_jobs, save):
        """With the env var at ``env_jobs`` (None: unset), render the card, then
        type ``save`` and press Save; return (card blurb, Save feedback, the value
        now saved)."""
        with tempfile.TemporaryDirectory() as cfg, \
                mock.patch.dict(os.environ, {"CRISPRME_CONFIG_DIR": cfg}):
            os.environ.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)
            if env_jobs is not None:
                os.environ["CRISPRME_MAX_CONCURRENT_JOBS"] = env_jobs
            blurb = self.render_card()[1]
            feedback = self.save(1, save).args[0]
            with open(os.path.join(cfg, a.CONFIG_FILE)) as fh:
                saved = json.load(fh)["max_concurrent_jobs"]
        return blurb, feedback, saved

    def test_override_and_value_in_effect_are_named(self):
        blurb, feedback, saved = self._settings("5", save=3)
        self.assertEqual(saved, 3)  # still saved, for when the env var is dropped
        for text in (blurb, feedback):
            self.assertIn("CRISPRME_MAX_CONCURRENT_JOBS", text)
            self.assertIn("overrides", text)
            self.assertRegex(text, r"\b5 concurrent search job\(s\) in effect")
        self.assertNotIn("Restart the app to apply.", feedback)

    def test_clamped_override_names_the_clamped_value(self):
        blurb, feedback, _ = self._settings("99", save=3)
        for text in (blurb, feedback):
            self.assertRegex(text, r"\b16 concurrent search job\(s\) in effect")

    def test_no_override_keeps_the_plain_texts(self):
        # Unset, or not a whole number: then the env var is ignored (app.py warns
        # about it) and the saved value does apply after a restart. Superscript
        # two and circled one pass str.isdigit() but not int(): Save raised on
        # them (HTTP 500, no feedback, though the value was written).
        plain_save = (
            "Saved: up to 3 concurrent search job(s). Restart the app to apply."
        )
        for env_jobs in (None, "", '"5"', "\u00b2", "\u2460"):
            with self.subTest(env_jobs=env_jobs):
                blurb, feedback, saved = self._settings(env_jobs, save=3)
                self.assertTrue(
                    blurb.endswith("Takes effect after the app restarts."), blurb
                )
                self.assertEqual(feedback, plain_save)
                self.assertEqual(saved, 3)

    def test_save_reports_what_it_cannot_persist(self):
        # A config dir it cannot write, or a value that is not a number (only a
        # hand-made request sends one), used to raise out of the callback.
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(
            os.environ, {"CRISPRME_CONFIG_DIR": os.path.join(tmp, "missing")}
        ):
            os.environ.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)
            for value in (3, "three"):
                with self.subTest(value=value):
                    feedback = self.save(1, value).args[0]
                    self.assertTrue(
                        feedback.startswith("Could not persist the setting"), feedback
                    )


if __name__ == "__main__":
    unittest.main()
