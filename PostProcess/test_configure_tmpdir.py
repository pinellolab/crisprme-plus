"""Unit test for crisprme.py::_configure_tmpdir (the $CRISPRME_TMPDIR knob) and
crisprme.py::_configure_mplconfigdir (a writable folder for matplotlib).

crisprme.py runs its command dispatch at import time, so it cannot be imported in a
unit test; following test_indel_pam_creation, the function's source is extracted and
exec'd in a controlled namespace instead.
"""

import ast
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
CRISPRME = os.path.join(HERE, "..", "crisprme.py")


def _extract_func(name):
    """The source of crisprme.py's top-level function ``name`` (not the module-level code
    after it: the ones run before the imports are followed by the command's setup)."""
    with open(CRISPRME) as fh:
        src = fh.read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(src, node)
    raise AssertionError(f"could not extract {name} from crisprme.py")


class TestConfigureTmpdir(unittest.TestCase):
    def setUp(self):
        ns = {"os": os, "sys": sys}
        exec(_extract_func("_configure_tmpdir"), ns)
        self.fn = ns["_configure_tmpdir"]
        self._saved = (os.environ.get("CRISPRME_TMPDIR"),
                       os.environ.get("TMPDIR"), tempfile.tempdir)

    def tearDown(self):
        c, t, td = self._saved
        for k, v in (("CRISPRME_TMPDIR", c), ("TMPDIR", t)):
            os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)
        tempfile.tempdir = td

    def test_sets_tmpdir_and_creates_missing_dir(self):
        with tempfile.TemporaryDirectory() as d:
            sub = os.path.join(d, "crisprme_tmp")  # does not exist yet -> created
            os.environ["CRISPRME_TMPDIR"] = sub
            os.environ.pop("TMPDIR", None)
            self.fn()
            self.assertTrue(os.path.isdir(sub))
            self.assertEqual(os.environ.get("TMPDIR"), sub)
            self.assertEqual(tempfile.tempdir, sub)

    def test_noop_when_unset(self):
        os.environ.pop("CRISPRME_TMPDIR", None)
        os.environ.pop("TMPDIR", None)
        tempfile.tempdir = None
        self.fn()
        self.assertIsNone(os.environ.get("TMPDIR"))
        self.assertIsNone(tempfile.tempdir)

    def test_unwritable_dir_warns_and_falls_back(self):
        os.environ["CRISPRME_TMPDIR"] = "/proc/nonexistent/cannot_create"
        os.environ.pop("TMPDIR", None)
        tempfile.tempdir = None
        self.fn()  # must not raise; must not set TMPDIR to an unusable dir
        self.assertIsNone(os.environ.get("TMPDIR"))


# what matplotlib says on stderr when it has to make a temporary folder of its own
_MPL_SAYS = "Matplotlib created a temporary cache directory"


class TestMatplotlibGetsAWritableFolder(unittest.TestCase):
    """With a home that cannot hold matplotlib's folders (apptainer --no-home, or --home on
    a read-only folder; docker --read-only), every command printed on stderr that matplotlib
    had made a temporary folder ('... is not a writable directory', 'Matplotlib created a
    temporary cache directory ...'), as crisprme.py imports it, and the fc-list its font-cache
    build runs that fontconfig had no cache folder it could write ($XDG_CACHE_HOME/fontconfig).
    A pipeline step's stderr is its error log, and the step fails when that is not empty; the
    steps stayed silent only as matplotlib passed them its folder. crisprme.py now gives
    matplotlib a temporary folder then ($MPLCONFIGDIR, removed at exit), and fontconfig too
    ($XDG_CACHE_HOME) when the home's ~/.cache cannot be written, before any command runs;
    with a writable home, or $MPLCONFIGDIR set, nothing changes. A home under a file stands
    for one that cannot be written in, as root too (the image's tests run as root, who writes
    in a folder of any mode)."""

    def setUp(self):
        ns = {"os": os, "sys": sys, "shutil": shutil}
        for name in ("_can_write_folder", "_configure_mplconfigdir"):
            exec(_extract_func(name), ns)
        self.fn = ns["_configure_mplconfigdir"]
        self.top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.top, ignore_errors=True)
        with open(os.path.join(self.top, "a_file"), "w"):
            pass
        self.unusable = os.path.join(self.top, "a_file", "home")
        self.env = {k: v for k, v in os.environ.items()
                    if k not in ("MPLCONFIGDIR", "XDG_CONFIG_HOME", "XDG_CACHE_HOME")}

    def configure(self, home, **env):
        """_configure_mplconfigdir with HOME ``home`` (and ``env``): the environment it
        leaves, and the functions it has run at exit."""
        at_exit = []
        with mock.patch.dict(os.environ, dict(self.env, HOME=home, **env), clear=True), \
                mock.patch("atexit.register", lambda f, *a, **k: at_exit.append((f, a, k))):
            self.fn()
            return dict(os.environ), at_exit

    def matplotlib_says(self, env):
        """What matplotlib writes on stderr in ``env`` when it finds its folders (as an
        import of it does; not its font cache, which it builds on its own time)."""
        return subprocess.run(
            [sys.executable, "-c", "import matplotlib; matplotlib.get_cachedir()"], env=env,
            capture_output=True, text=True, timeout=300).stderr

    @unittest.skipUnless(subprocess.run([sys.executable, "-c", "import matplotlib"],
                                        capture_output=True).returncode == 0,
                         "no matplotlib here")
    def test_a_home_that_cannot_hold_them(self):
        # what it did: matplotlib's own temporary folder, said on stderr
        self.assertIn(_MPL_SAYS, self.matplotlib_says(dict(self.env, HOME=self.unusable)))
        env, at_exit = self.configure(self.unusable)
        folder = env.get("MPLCONFIGDIR")
        if folder:
            self.addCleanup(shutil.rmtree, folder, ignore_errors=True)
        self.assertTrue(folder and os.path.isdir(folder) and os.access(folder, os.W_OK), env)
        self.assertEqual(os.path.dirname(folder), tempfile.gettempdir())
        self.assertEqual(env.get("XDG_CACHE_HOME"), folder)  # (~/.cache cannot be written)
        self.assertEqual(self.matplotlib_says(env), "")
        # removed at exit
        self.assertEqual(at_exit, [(shutil.rmtree, (folder,), {"ignore_errors": True})])
        for f, a, k in at_exit:
            f(*a, **k)
        self.assertFalse(os.path.exists(folder))

    def test_a_home_that_holds_them_is_left_as_it_is(self):
        home = os.path.join(self.top, "home")
        os.makedirs(home)
        for given in ({}, {"XDG_CACHE_HOME": os.path.join(self.top, "cache")}):
            with self.subTest(**given):
                env, at_exit = self.configure(home, **given)
                self.assertNotIn("MPLCONFIGDIR", env)
                self.assertEqual(env.get("XDG_CACHE_HOME"), given.get("XDG_CACHE_HOME"))
                self.assertEqual(at_exit, [])
        self.assertEqual(os.listdir(home), [])  # (nothing made there)
        # nor one that is not there yet, in a folder that can hold it
        env, at_exit = self.configure(os.path.join(self.top, "new_home"))
        self.assertNotIn("MPLCONFIGDIR", env)

    def test_a_folder_given_is_kept(self):
        env, at_exit = self.configure(self.unusable, MPLCONFIGDIR=self.top)
        self.assertEqual((env["MPLCONFIGDIR"], env.get("XDG_CACHE_HOME"), at_exit),
                         (self.top, None, []))
        # and a cache folder given, when matplotlib's folders cannot be made
        cache = os.path.join(self.top, "cache")
        env, at_exit = self.configure(self.unusable, XDG_CACHE_HOME=cache)
        self.addCleanup(shutil.rmtree, env["MPLCONFIGDIR"], ignore_errors=True)
        self.assertEqual(env["XDG_CACHE_HOME"], cache)

    def test_a_cache_folder_that_cannot_be_made(self):
        # the config folder can be, the cache folder ($XDG_CACHE_HOME, which matplotlib
        # reads on Linux; on macOS it keeps both in ~/.matplotlib) cannot
        home = os.path.join(self.top, "home")
        os.makedirs(home)
        env, at_exit = self.configure(home, XDG_CACHE_HOME=self.unusable)
        if env.get("MPLCONFIGDIR"):
            self.addCleanup(shutil.rmtree, env["MPLCONFIGDIR"], ignore_errors=True)
        self.assertEqual("MPLCONFIGDIR" in env, sys.platform.startswith(("linux", "freebsd")))
        self.assertEqual(env["XDG_CACHE_HOME"], self.unusable)  # (given: kept)

    @unittest.skipUnless(
        subprocess.run([sys.executable, "-c", "import Bio.Seq, pandas, matplotlib"],
                       capture_output=True).returncode == 0,
        "crisprme.py's and generate-report's dependencies not installed")
    def test_generate_report_says_nothing_of_it(self):
        # the command, end to end: nothing on stderr, and the folder gone when it ends
        import test_generate_report as fixture

        res_dir = os.path.join(self.top, "Results", "job")
        os.makedirs(res_dir)
        tsv = os.path.join(res_dir, f"{fixture._GUIDE}+NRG_hg38+hg38_1000G_HGDP_6+3"
                                    "_integrated_results.tsv")
        with open(tsv, "w") as fh:
            fh.write("\t".join(fixture._HEADER) + "\n")
            fh.writelines("\t".join(row) + "\n" for row in fixture._ROWS)
        tmp = os.path.join(self.top, "tmp")
        os.makedirs(tmp)
        run = subprocess.run(
            [sys.executable, CRISPRME, "generate-report", "--integrated-results", tsv,
             "--output", os.path.join(self.top, "report.zip")],
            env=dict(self.env, HOME=self.unusable, TMPDIR=tmp), cwd=self.top,
            capture_output=True, text=True, timeout=600)
        # (matplotlib notes that it builds its font cache when that takes over 5 s, as on a
        # host with many fonts: into the new folder, by crisprme.py itself, which imports
        # generate_report before any command runs, so once a run and before any step)
        said = [line for line in run.stderr.splitlines()
                if line != "Matplotlib is building the font cache; this may take a moment."]
        self.assertEqual((run.returncode, said), (0, []), run.stderr)
        self.assertTrue(os.path.isfile(os.path.join(self.top, "report.zip")))
        self.assertEqual(os.listdir(tmp), [])

    def test_every_command_gets_it(self):
        # at module level, after $CRISPRME_TMPDIR is honoured, and before the PostProcess
        # modules are imported (generate_report imports matplotlib) and any command runs
        with open(CRISPRME) as fh:
            src = fh.read()
        calls = [m.start() for m in re.finditer(r"^_configure_mplconfigdir\(\)", src, re.M)]
        self.assertEqual(len(calls), 1)
        self.assertLess(src.index("\n_configure_tmpdir()"), calls[0])
        for later in ("\nsys.path.insert(0, script_path)", "\nfrom generate_report import",
                      "\nif len(sys.argv) < 2:"):
            with self.subTest(before=later.strip()):
                self.assertLess(calls[0], src.index(later))


if __name__ == "__main__":
    unittest.main()
