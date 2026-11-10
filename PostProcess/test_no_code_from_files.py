#!/usr/bin/env python3
"""The ways the source shows of running code that a file or a folder holds, and of temporary
names another user could make first, kept out of CRISPRme (2.10.0).

A search writes its files in its output folder and reads some of them back, and the folders
it works in can be shared: an output or data folder under a umask of 002 (its group can write
it), a lab's data folder, /tmp. The post-analysis in pieces once handed its pieces a pickle of
the chromosome's dictionary keys in the output folder, and loading a pickle runs whatever it
names, so a member of that folder's group could run code as the user who ran the search (the
review of lane S, round 3, S-R3-1). The copy is JSON now (test_postproc_pieces: a planted
pickle runs nothing). Round 4 found two more of the class, older: every complete-search ran
the scorer's health check as "python -c" from the data folder, which put that folder first on
Python's path (S-R4-1), and the web app's Settings Retry ran through a shell the command a
job's cmd.json held (S-R4-2; test_settings_progress); and the sweep that followed, values of
a job's files that went into shell commands as the shell's own text.

Over the parse tree of every Python file of the repository but the tests (test_*.py and
test/), there is:
  * no import of a module that loads objects which can run code (pickle, cPickle, _pickle,
    marshal, shelve, dill, cloudpickle, joblib), nor such an import by name; the data a run
    writes and reads back is JSON or text;
  * no call of the builtins eval or exec;
  * no PyYAML load other than safe_load, no pandas read_pickle, no allow_pickle=True;
  * no tempfile.mktemp: it gives a name, not a file, which another user of a shared folder can
    make first (a link to a file of the user's, which the run then writes through);
    tempfile.mkstemp, mkdtemp and TemporaryDirectory make their own;
  * no Python run on -c, -m or - (in an argument list or a shell command) that is not
    isolated (-I): those put the working directory first on its path.
And in every shell script (but comments): no eval, no mktemp -u / --dry-run (a name only),
no name made by adding to the one mktemp made ("$(mktemp ...).ext": the file mktemp made
is another), and no Python run on -c, -m or - that is not isolated. The Dockerfile's build
step and the CI workflows, which run in folders the build and the checkout make, are not
scanned.

WhatASearchReadsBackRunsNothing runs the scripts the sweep fixed on a line or a name another
user could put in a search's folder, and TheDebugFlagRunsNoModuleOfTheWorkingDirectory
crisprme.py --debug from a folder holding a PostProcess/ of someone else's. These are the forms a scan of the source and those runs
show: a folder others can write still lets them swap a folder a search made there, put links
it writes through and change what older parts of the pipeline read, so keep data, output and
temporary folders writable by you alone (docs/TROUBLESHOOTING.md section 7).

Run with:
    cd PostProcess && python3 -m unittest test_no_code_from_files -v
"""

import ast
import os
import re
import shutil
import subprocess
import sys
import tempfile
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


# Python run on a program in its argument (-c), on a module (-m) or on its standard input (-)
# puts the working directory first on sys.path, so a module another user put there runs on
# its first import; isolated (-I), the working directory, the PYTHON* variables and the
# user's site-packages stay out of its path. A search runs from the data folder, which a
# umask of 002 leaves to its group, and the scorer's health check ran there as "python -c"
# (2.10.0, the review of lane S, round 4, S-R4-1).
_PY_TOKEN = re.compile(r"(?:^|/)python[0-9.]*$")
_PY_VAR = re.compile(r"^\$\{?[A-Za-z_]*(?:PYTHON|python|PY|py)[A-Za-z0-9_]*\}?$")
_PY_NAMES = {"py", "python", "python3", "pyexe", "py_exe", "python_exe", "python_bin",
             "interpreter", "env_py", "env_python"}
_RUNS_FROM_CWD = re.compile(r"^-[A-Za-z]*[cm]$")  # -c, -m, or a group ending in one (-Bc)
_FLAGS = re.compile(r"^-[A-Za-z]+$")


def _not_isolated(tokens):
    """Whether ``tokens``, what follows a python, run it on -c, -m or - with no -I before."""
    flags = []
    for tok in tokens:
        if tok is None or not tok.startswith("-"):
            return False
        flags.append(tok)
        if tok == "-" or _RUNS_FROM_CWD.match(tok):
            return not any("I" in f[1:] for f in flags)
        if not _FLAGS.match(tok):
            return False
    return False


def _is_python(node):
    """Whether the AST ``node`` names a Python interpreter: a string "python", "python3",
    ".../python3.10"; sys.executable; a name such as py or env_python; env_python(...)."""
    text = _const_str(node)
    if text is not None:
        return bool(_PY_TOKEN.search(text))
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        if node.value.id == "sys" and node.attr == "executable":
            return True
    if isinstance(node, ast.Name):
        return node.id in _PY_NAMES
    if isinstance(node, ast.Attribute):
        return node.attr in _PY_NAMES
    if isinstance(node, ast.Call):
        func = node.func
        name = func.id if isinstance(func, ast.Name) else getattr(func, "attr", "")
        return name == "env_python"
    return False


def _command_findings(text):
    """(offset, what) for each python in the shell command ``text`` run on -c, -m or -
    without -I."""
    out = []
    words = [(m.start(), m.group().strip("'\"")) for m in re.finditer(r"[^\s;&|()`]+", text)]
    for i, (at, word) in enumerate(words):
        if _PY_TOKEN.search(word) or _PY_VAR.match(word):
            if _not_isolated([w for _, w in words[i + 1:]]):
                out.append((at, "python -c, -m or - without -I"))
    return out


def python_isolation_findings(source, name="<source>"):
    """(line, what) for each python ``source`` runs on -c, -m or - without -I: in an argument
    list ([sys.executable, "-c", code]) or in a shell command (a string, an f-string with
    {sys.executable} in it), docstrings aside."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = ast.parse(source, name)
    skip = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and _const_str(body[0].value) is not None:
                skip.add(id(body[0].value))
        if isinstance(node, ast.JoinedStr):
            skip.update(id(v) for v in node.values)
    out = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.List, ast.Tuple)):
            elts = node.elts
            for j, elt in enumerate(elts):
                if _is_python(elt) and _not_isolated([_const_str(e) for e in elts[j + 1:]]):
                    out.append((node.lineno, "python -c, -m or - without -I"))
        elif isinstance(node, ast.JoinedStr):
            text = "".join(_const_str(v) if _const_str(v) is not None
                           else ("python" if _is_python(v.value) else "X")
                           for v in node.values)
            out += [(node.lineno, what) for _, what in _command_findings(text)]
        elif _const_str(node) is not None and id(node) not in skip:
            out += [(node.lineno, what) for _, what in _command_findings(node.value)]
    return sorted(out)


def shell_isolation_findings(text):
    """(line, what) for each python a shell script runs on -c, -m or - without -I, but in
    comments."""
    out = []
    for i, line in enumerate(text.splitlines(), 1):
        code = line.split("#", 1)[0] if not line.lstrip().startswith("#!") else ""
        out += [(i, what) for _, what in _command_findings(code)]
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

    def test_every_python_run_on_a_program_or_a_module_runs_isolated(self):
        found = []
        for path in _shipped(".py"):
            with open(path, encoding="utf-8") as fh:
                found += ["%s:%d: %s" % (os.path.relpath(path, REPO), line, what)
                          for line, what in python_isolation_findings(fh.read(), path)]
        for path in _shipped(".sh"):
            with open(path, encoding="utf-8", errors="replace") as fh:
                found += ["%s:%d: %s" % (os.path.relpath(path, REPO), line, what)
                          for line, what in shell_isolation_findings(fh.read())]
        self.assertEqual(found, [])


class WhatASearchReadsBackRunsNothing(unittest.TestCase):
    """A search reads back files of its output folder, which others may write (its group,
    under a umask of 002), and the names in it; up to 2.10.0 some of what it read went into
    shell commands as the shell's own text (the review of lane S, round 4): each runs, from
    the folder a search runs it in, nothing a line or a name holds."""

    GUIDE = "CTAACAGTTGCTTTTATCACNNN"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def markers(self):
        return sorted(os.path.relpath(os.path.join(root, f), self.tmp)
                      for root, _dirs, files in os.walk(self.tmp)
                      for f in files if f.startswith("MARK"))

    def test_the_summaries_of_a_guide_list_with_a_line_put_in_it(self):
        # process_summaries.py looks each guide of the search's guide file up with grep,
        # sort, mv and rm, through a shell
        row = ["X", self.GUIDE, "CTAACAGTTGCTTTTATCACTGG", "NA", "chr1", "100", "100", "+",
               "0", "0", "0", "NA", "NA", "NA", "NA", self.GUIDE, "NA", "NA", "NA", "1",
               "1.0", "1.0"]
        best = os.path.join(self.tmp, "best.txt")
        with open(best, "w") as fh:
            fh.write("\t".join(row) + "\n")
        guides = os.path.join(self.tmp, "guides.txt")
        planted = "CTAACAG;touch${IFS}MARK_semicolon;#"
        with open(guides, "w") as fh:
            fh.write(f"{self.GUIDE}\n{planted}\nCTAACAG$(touch${{IFS}}MARK_substitution)\n")
        out = os.path.join(self.tmp, "job")
        os.makedirs(out)
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "process_summaries.py"), best, guides, "_",
             "4", "1", "1", out, "ref", "CFD"],
            cwd=self.tmp, stdin=subprocess.DEVNULL, capture_output=True, text=True)
        self.assertEqual(self.markers(), [])
        self.assertEqual(done.returncode, 0, done.stderr)
        with open(os.path.join(out, f"job.summary_by_guide.{self.GUIDE}_CFD.txt")) as fh:
            self.assertEqual(fh.read().splitlines()[1], "X\t0\t0\t1\t0\t1\t0")  # its target

    def test_the_cfd_graphs_of_a_folder_with_a_name_put_in_it(self):
        # assemble_cfd_graphs.py (post-analysis) sums the folder's CFD graph files and
        # removes each, through a shell
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        graph = "ref\tvar\n" + "1\t2\n" * 101
        names = ("chr1.snps.CFDGraph.txt", "x;touch${IFS}MARK_name;.snps.CFDGraph.txt",
                 "fake$(touch${IFS}MARK_fake).CFDGraph.txt")
        for name in names:
            with open(os.path.join(out, name), "w") as fh:
                fh.write(graph)
        done = subprocess.run([sys.executable, os.path.join(HERE, "assemble_cfd_graphs.py"), out],
                              cwd=self.tmp, stdin=subprocess.DEVNULL, capture_output=True,
                              text=True)
        self.assertEqual(self.markers(), [])
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(sorted(os.listdir(out)), ["indels.CFDGraph.txt", "snps.CFDGraph.txt"])
        with open(os.path.join(out, "snps.CFDGraph.txt")) as fh:
            self.assertEqual(fh.read().splitlines()[1], "2\t4")


try:  # crisprme.py imports Biopython at once, which the light unit-test CI does not have
    import Bio  # noqa: F401
except ImportError:  # pragma: no cover
    Bio = None


@unittest.skipUnless(Bio, "Biopython not installed: crisprme.py does not start")
class TheDebugFlagRunsNoModuleOfTheWorkingDirectory(unittest.TestCase):
    """crisprme.py --debug took the working directory's PostProcess/ for its own whenever there
    was one, and loaded its memory_check.py at start-up: run from a data folder its group can
    write (a umask of 002), a PostProcess/ another member put there ran as the user, whatever
    the subcommand (2.10.0, the review of lane S, round 4). It keeps the PostProcess/ beside
    the crisprme.py that runs."""

    def test_a_planted_postprocess_runs_nothing(self):
        data = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, data, ignore_errors=True)
        os.makedirs(os.path.join(data, "PostProcess"))
        for mod in ("memory_check", "validate_inputs", "genome_layout"):
            with open(os.path.join(data, "PostProcess", mod + ".py"), "w") as fh:
                fh.write("import os\nopen(os.path.join(%r, 'MARK_%s'), 'w').close()\n"
                         % (data, mod))
        done = subprocess.run([sys.executable, os.path.join(REPO, "crisprme.py"), "--version",
                               "--debug"], cwd=data, stdin=subprocess.DEVNULL,
                              capture_output=True, text=True, timeout=600)
        self.assertEqual(sorted(f for f in os.listdir(data) if f.startswith("MARK")), [])
        self.assertIn("DEBUG MODE", done.stdout)


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

    def test_python_run_from_the_working_directory(self):
        for source in ('subprocess.run([py, "-c", probe])',
                       'subprocess.run([sys.executable, "-m", "pip", "install", x])',
                       '_run([exe, "run", "-n", name, "python", "-c", "import sys"])',
                       'subprocess.run(["python3", "-u", "-c", code])',
                       'subprocess.run([env_python(n), "-Bc", code])',
                       'subprocess.run(["/opt/conda/bin/python3.10", "-", "x"], input=code)',
                       'os.system("python -c \'import annotation\' " + f)',
                       'os.system(f"{sys.executable} -m json.tool {path}")',
                       'subprocess.call(f"cd {d} && python3 -E -s -c \'print(1)\'", shell=True)'):
            with self.subTest(source=source):
                self.assertEqual(len(python_isolation_findings(source)), 1, source)
        for source in ('subprocess.run([py, "-I", "-c", probe])',
                       'subprocess.run([sys.executable, "-Im", "pip"])',
                       'subprocess.run([py, "-u", worker])',
                       'subprocess.run(["gunzip", "-c", path])',
                       'args += ["-c", channel]',
                       'subprocess.run(["bash", "-c", script])',
                       'subprocess.call([sys.executable, script, "-c", "x"])',
                       'os.system(f"{sys.executable} {script} -m 4")',
                       'os.system("python -I -c \'import annotation\' " + f)',
                       'def f():\n    """Run as ``python -m PostProcess.generate_report``."""'):
            with self.subTest(source=source):
                self.assertEqual(python_isolation_findings(source), [], source)
        for text in ('python -c "import annotation,sys" "$f"',
                     'speedup=$(python3 -c "print(1)")', 'x=1; python3 -m pip install y',
                     '"$PYTHON" -c "$code"', 'python3 - <<EOF', '${python_bin} -u -c x'):
            with self.subTest(text=text):
                self.assertEqual(len(shell_isolation_findings(text)), 1, text)
        for text in ('python -I -c "import annotation,sys" "$f"', 'python3 -Ic "print(1)"',
                     'python annotation.py "$f"', './process_summaries.py -c 3',
                     '# python -c "x"', 'rm -f lib/python*/site-packages/pyproject.toml',
                     'python $starting_dir/populations_distribution.py "$f" $total'):
            with self.subTest(text=text):
                self.assertEqual(shell_isolation_findings(text), [], text)


if __name__ == "__main__":
    unittest.main()
