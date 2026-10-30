"""Tests for the configurable concurrent-search-jobs setting (app_settings.py).

Default stays a conservative 2; env var overrides a persisted value; values are
clamped to [1, MAX_JOBS_CAP] so the pool is never oversubscribed. app.py, loaded
the way the web app is really started, must size its search pool from it."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest

import app_settings as a

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Loads app.py the way the web app really loads it. `crisprme.py web-interface`
# runs `subprocess.run([sys.executable, index.py])`: a fresh interpreter whose
# sys.path holds only the install root (index.py's directory, no PostProcess/),
# whose cwd is the user's working directory, and whose index.py imports app.py
# first, before pages/pages_utils.py puts PostProcess/ on sys.path. -I drops this
# test's own cwd (PostProcess/) and any PYTHONPATH; the root is then put back as
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


class TestAppSearchPoolSize(unittest.TestCase):
    """app.py used to import app_settings with only the install root on sys.path,
    so under `crisprme.py web-interface` the import always failed and a bare
    except silently pinned the search pool to 2, whatever the env var or the
    Settings page said."""

    def _import_app(self, root=_REPO_ROOT, env_jobs=None, saved_jobs=None):
        """Import <root>/app.py like the web app does; return (search pool
        max_workers, everything app.py printed on stdout)."""
        env = dict(os.environ)
        env.pop("CRISPRME_CONFIG_DIR", None)
        env.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)
        if env_jobs is not None:
            env["CRISPRME_MAX_CONCURRENT_JOBS"] = str(env_jobs)
        with tempfile.TemporaryDirectory() as workdir:  # the user's cwd, not the repo
            if saved_jobs is not None:  # what the Settings page persists
                a.set_max_concurrent_jobs(saved_jobs, workdir)
            res = subprocess.run(
                [sys.executable, "-I", "-c", _APP_IMPORT_DRIVER, root],
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


if __name__ == "__main__":
    unittest.main()
