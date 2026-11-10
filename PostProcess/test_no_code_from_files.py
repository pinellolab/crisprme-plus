#!/usr/bin/env python3
"""No CRISPRme code loads a file in a way that can run code, nor takes a temporary name that
another user can make first (2.10.0).

A search writes its files in its output folder and reads some of them back, and the folders
it works in can be shared: an output or data folder under a umask of 002 (its group can write
it), a lab's data folder, /tmp. The post-analysis in pieces once handed its pieces a pickle of
the chromosome's dictionary keys in the output folder, and loading a pickle runs whatever it
names, so a member of that folder's group could run code as the user who ran the search (the
review of lane S, round 3, S-R3-1). The copy is JSON now (test_postproc_pieces: a planted
pickle runs nothing); this module keeps the class out of every module.

Over the parse tree of every Python file of the repository but the tests (test_*.py and
test/), there is:
  * no import of a module that loads objects which can run code (pickle, cPickle, _pickle,
    marshal, shelve, dill, cloudpickle, joblib), nor such an import by name; the data a run
    writes and reads back is JSON or text;
  * no call of the builtins eval or exec;
  * no PyYAML load other than safe_load, no pandas read_pickle, no allow_pickle=True;
  * no tempfile.mktemp: it gives a name, not a file, which another user of a shared folder can
    make first (a link to a file of the user's, which the run then writes through);
    tempfile.mkstemp, mkdtemp and TemporaryDirectory make their own.
And in every shell script (but comments): no eval, no mktemp -u / --dry-run (a name only),
and no name made by adding to the one mktemp made ("$(mktemp ...).ext": the file mktemp made
is another).

Run with:
    cd PostProcess && python3 -m unittest test_no_code_from_files -v
"""

import ast
import os
import re
import unittest
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# modules whose loading of a file or of bytes builds objects that can run code
CODE_LOADERS = {"pickle", "cPickle", "_pickle", "marshal", "shelve", "dill", "cloudpickle",
                "joblib"}
_YAML_LOADS = {"load", "unsafe_load", "full_load", "load_all", "unsafe_load_all",
               "full_load_all"}
_SKIP_DIRS = {".git", "__pycache__", "node_modules"}


def _shipped(suffix):
    """The repository's files ending in ``suffix``, but the tests (test_*.py, test/)."""
    out = []
    for dirpath, dirnames, names in os.walk(REPO):
        rel = os.path.relpath(dirpath, REPO)
        dirnames[:] = sorted(d for d in dirnames if d not in _SKIP_DIRS
                             and not (rel == "." and d == "test"))
        for name in sorted(names):
            if name.endswith(suffix) and not name.startswith("test_"):
                out.append(os.path.join(dirpath, name))
    return out


def _const_str(node):
    return node.value if isinstance(node, ast.Constant) and isinstance(node.value, str) else None


def python_findings(source, name="<source>"):
    """(line, what) for each use in ``source`` of a way to load code or to take a temporary
    name, as the module's docstring lists them."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # an old invalid escape is no finding
        tree = ast.parse(source, name)
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] in CODE_LOADERS:
                    out.append((node.lineno, "import " + alias.name))
        elif isinstance(node, ast.ImportFrom):
            root = (node.module or "").split(".")[0]
            if root in CODE_LOADERS:
                out.append((node.lineno, "from %s import ..." % node.module))
            if node.module == "tempfile" and any(a.name == "mktemp" for a in node.names):
                out.append((node.lineno, "from tempfile import mktemp"))
        elif isinstance(node, ast.Call):
            func = node.func
            owner = func.value.id if (isinstance(func, ast.Attribute)
                                      and isinstance(func.value, ast.Name)) else None
            attr = func.attr if isinstance(func, ast.Attribute) else None
            if isinstance(func, ast.Name) and func.id in ("eval", "exec"):
                out.append((node.lineno, func.id + "()"))
            if (isinstance(func, ast.Name) and func.id == "__import__") or attr == "import_module":
                target = _const_str(node.args[0]) if node.args else None
                if target and target.split(".")[0] in CODE_LOADERS:
                    out.append((node.lineno, "import of %s by name" % target))
            if owner == "tempfile" and attr == "mktemp":
                out.append((node.lineno, "tempfile.mktemp()"))
            if owner == "yaml" and attr in _YAML_LOADS:
                out.append((node.lineno, "yaml.%s()" % attr))
            if attr == "read_pickle":
                out.append((node.lineno, "read_pickle()"))
            for kw in node.keywords:
                if kw.arg == "allow_pickle" and not (isinstance(kw.value, ast.Constant)
                                                     and kw.value.value is False):
                    out.append((node.lineno, "allow_pickle=%s" % ast.dump(kw.value)))
    return sorted(out)


_SH_EVAL = re.compile(r"(^|[;&|({`]|\bthen|\bdo|\belse)\s*eval\b")
_SH_DRY = re.compile(r"\bmktemp\b[^;&|)\n]*\s(-u\b|--dry-run\b|-[a-zA-Z]*u[a-zA-Z]*\b)")
_SH_BESIDE = re.compile(r"\$\(\s*mktemp\b[^)]*\)[^\s\"';)|&]")


def shell_findings(text):
    """(line, what) for each eval, dry-run mktemp, or name made beside mktemp's, but in
    comments."""
    out = []
    for i, line in enumerate(text.splitlines(), 1):
        code = line.split("#", 1)[0] if not line.lstrip().startswith("#!") else ""
        if not code.strip():
            continue
        if _SH_EVAL.search(code):
            out.append((i, "eval"))
        if _SH_DRY.search(code):
            out.append((i, "mktemp -u (a name only)"))
        if _SH_BESIDE.search(code):
            out.append((i, "a name beside the file mktemp made"))
    return out


class NoCodeIsLoadedFromAFile(unittest.TestCase):
    def test_no_python_module_does(self):
        files = _shipped(".py")
        self.assertIn(os.path.join(HERE, "new_simple_analysis.py"), files)
        self.assertIn(os.path.join(REPO, "crisprme.py"), files)
        found = []
        for path in files:
            with open(path, encoding="utf-8") as fh:
                source = fh.read()
            found += ["%s:%d: %s" % (os.path.relpath(path, REPO), line, what)
                      for line, what in python_findings(source, path)]
        self.assertEqual(found, [])

    def test_no_shell_script_does(self):
        files = _shipped(".sh")
        self.assertIn(os.path.join(REPO, "docs", "install.sh"), files)
        found = []
        for path in files:
            with open(path, encoding="utf-8", errors="replace") as fh:
                found += ["%s:%d: %s" % (os.path.relpath(path, REPO), line, what)
                          for line, what in shell_findings(fh.read())]
        self.assertEqual(found, [])


class TheScanSeesEachForm(unittest.TestCase):
    """The scan itself: each form it is there for is found, and what is safe is not."""

    def test_python_forms(self):
        for source in ("import pickle", "import cPickle as p", "from pickle import loads",
                       "import marshal", "import shelve", "import dill",
                       "from joblib import load", "import _pickle",
                       "__import__('pickle')", "importlib.import_module('shelve')",
                       "eval(x)", "exec(code)", "tempfile.mktemp(suffix='.fa')",
                       "from tempfile import mktemp", "yaml.load(s)",
                       "yaml.unsafe_load(s)", "pd.read_pickle(p)",
                       "np.load(p, allow_pickle=True)"):
            with self.subTest(source=source):
                self.assertEqual(len(python_findings(source)), 1, source)
        for source in ("import json\njson.load(fh)", "np.load(p)",
                       "np.load(p, allow_pickle=False)", "yaml.safe_load(s)",
                       "tempfile.mkstemp()", "tempfile.mkdtemp(prefix='x')",
                       "re.compile('eval')", "model.eval()", "obj.exec()",
                       "'import pickle'"):
            with self.subTest(source=source):
                self.assertEqual(python_findings(source), [], source)

    def test_shell_forms(self):
        for text in ('eval "$cmd"', 'x=1; eval "$y"', 'tmp=$(mktemp -u)',
                     'tmp=$(mktemp --dry-run /tmp/x.XXXX)', 'SRC="$(mktemp -t app).applescript"',
                     'if true; then eval "$x"; fi'):
            with self.subTest(text=text):
                self.assertEqual(len(shell_findings(text)), 1, text)
        for text in ('# eval is not run here', 'TMP="$(mktemp -d)"',
                     'SRC_DIR="$(mktemp -d -t app)"', 'WORK=$(mktemp -d -p "$OUT")',
                     'echo "evaluate"', 'x="$(mktemp)"; y="$x.done"'):
            with self.subTest(text=text):
                self.assertEqual(shell_findings(text), [], text)


if __name__ == "__main__":
    unittest.main()
