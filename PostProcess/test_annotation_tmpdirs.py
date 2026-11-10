#!/usr/bin/env python3
"""complete-search and assembly-search remove the annotation temp folders they make.

_sort_annotation (crisprme.py) sorts each --annotation / --gene_annotation into a new
crisprme_annot_* folder in the working directory ($CRISPRME_TMPDIR if set), and
_process_personal_annotation builds --personal_annotation in a crisprme_pannot_* one.
The pipeline reads those files, and gunzips back into them, until it ends. Nothing
removed them, so every annotated search left two behind, ~1.7 GB with the shipped
annotations.

--personal_annotation with --annotation is merged into one sorted file in the
crisprme_pannot_* folder; the merge looked for the plain copy of the sorted --annotation,
which bgzip had removed, and stopped the search. And the personal file was read as text,
so a bgzipped one, which --help asked for, stopped it too: it is now decompressed first.

A search whose CRISPR-Bulge scorer env is not ready stops before the pipeline, after the
annotations are sorted, as --help now says: it leaves none of their folders either.

These tests run complete-search the way the CLI does -- the dispatch statement is taken
verbatim from crisprme.py, whose definitions are exec'd without the dispatch (the
truncate-and-exec loader of test_bulge_cap) -- on a tiny synthetic genome, with the
pipeline launch stubbed to succeed or fail, and sort-bed / bgzip replaced by Python
stand-ins (the light CI env has neither). Every other step of complete-search runs.

A pipeline that may still be running keeps them: after a Ctrl-C, subprocess.call has
killed only its /bin/sh (dash does not exec the pipeline, which runs on; bash does, and
the step it was running runs on), and a signal that killed the pipeline's bash need not
have killed its workers. One test sends SIGINT to a real crisprme.py process alone while
the worker of a stub pipeline, which goes on under either shell, is still to read them.
"""
import contextlib
import gzip
import io
import os
import re
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import textwrap
import time
import types
import unittest
import zlib
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
CRISPRME_PY = os.path.join(os.path.dirname(HERE), "crisprme.py")
_MARKER = "\nif len(sys.argv) < 2:"


# --------------------------------------------------------------------------- #
# Load crisprme.py DEFINITIONS only (no CLI dispatch), Bio.Seq + scientific
# absents stubbed -- identical to test_bulge_cap._load_crisprme.
# --------------------------------------------------------------------------- #
def _load_crisprme():
    sys.path.insert(0, HERE)
    with open(CRISPRME_PY) as fh:
        src = fh.read()
    cut = src.index(_MARKER)  # KeyError-loud if the dispatch shape ever changes
    mod = types.ModuleType("crisprme_under_test_annot_tmpdirs")
    mod.__file__ = CRISPRME_PY
    code = compile(src[:cut], CRISPRME_PY, "exec")
    # Stubs exist ONLY so crisprme.py's top-level imports resolve while we exec its
    # definitions; they are removed again in the finally, so they cannot shadow the real
    # modules for a test collected later.
    _inserted = []
    if "Bio.Seq" not in sys.modules:
        _bio = types.ModuleType("Bio")
        _bio_seq = types.ModuleType("Bio.Seq")
        _bio_seq.Seq = object
        _bio.Seq = _bio_seq
        if "Bio" not in sys.modules:
            sys.modules["Bio"] = _bio
            _inserted.append("Bio")
        sys.modules["Bio.Seq"] = _bio_seq
        _inserted.append("Bio.Seq")
    for _m in ("pandas", "scipy", "sklearn", "matplotlib", "seaborn",
               "statsmodels", "intervaltree", "CRISPR_BULGE_score"):
        if _m in sys.modules:
            continue
        try:
            __import__(_m)
        except Exception:
            sys.modules[_m] = mock.MagicMock()
            _inserted.append(_m)
    try:
        exec(code, mod.__dict__)
    finally:
        for _m in _inserted:
            sys.modules.pop(_m, None)
    return mod


def _dispatch(subcommand):
    """The statement crisprme.py's CLI dispatch runs for `subcommand`, verbatim."""
    with open(CRISPRME_PY) as fh:
        tail = fh.read().split(_MARKER, 1)[1]
    m = re.search(rf'^elif sys\.argv\[1\] == "{re.escape(subcommand)}":[^\n]*\n'
                  r'((?:    [^\n]*\n)+)', tail, re.M)
    return textwrap.dedent(m.group(1))


crisprme = _load_crisprme()


def _sort_bed(fname, outfname):
    """sort-bed stand-in: the same lines, sorted."""
    with open(fname) as src, open(outfname, "w") as dst:
        dst.writelines(sorted(src))
    return outfname


def _compress_file(fname):
    """`bgzip -f` stand-in: gzip in place, removing the plain file as bgzip does."""
    with open(fname, "rb") as src, gzip.open(f"{fname}.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(fname)
    return f"{fname}.gz"


def _bgzf(data, block=65280):
    """The bytes `bgzip` writes for ``data`` (the light CI env has no bgzip): gzip members of
    at most 64 KiB, each with the BC extra field that holds the member's size less one, then
    bgzip's empty end-of-file member."""
    out = b""
    for i in range(0, len(data), block):
        chunk = data[i:i + block]
        deflate = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = deflate.compress(chunk) + deflate.flush()
        out += (b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6)
                + b"BC" + struct.pack("<HH", 2, len(body) + 25) + body
                + struct.pack("<II", zlib.crc32(chunk) & 0xFFFFFFFF, len(chunk)))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


# The pipeline, as far as the annotation folders go: it starts a worker in the
# background, as the real one starts its annotation.py workers, and that worker outlives
# it and its /bin/sh. Once told to (a "go" file), the worker reads both sorted
# annotations ($5 --annotation, $16 --gene_annotation) and says whether they were there.
_STUB_PIPELINE = r"""#!/bin/bash
d=$(dirname "$0")
(
	for _ in $(seq 1200); do [ -e "$d/go" ] && break; sleep 0.1; done
	if cat "$5" "${16}" >/dev/null 2>&1; then echo read; else echo gone; fi >"$d/worker.tmp"
	mv "$d/worker.tmp" "$d/worker.result"
) &
touch "$d/started"
wait
"""


def _search_until_interrupted(work, stub):
    """Run in a child process by the SIGINT test: complete-search as the CLI dispatch runs
    it, in `work`, sort-bed / bgzip stood in for, and the real subprocess.call launching
    `stub` in place of the pipeline script."""
    signal.signal(signal.SIGINT, signal.default_int_handler)  # whatever it inherited
    pipeline = os.path.join(crisprme.script_path, "submit_job_automated_new_multiple_vcfs.sh")
    real_call = subprocess.call

    def call(cmd, *a, **k):
        if isinstance(cmd, str):
            cmd = cmd.replace(pipeline, stub)
        return real_call(cmd, *a, **k)

    argv = ["crisprme.py", "complete-search", "--genome", "Genomes/synth",
            "--pam", "PAMs/20bp-NGG-SpCas9.txt", "--guide", "guide.txt",
            "--mm", "1", "--bDNA", "0", "--bRNA", "0", "--output", "int",
            "--thread", "1", "--annotation", "myann/func.bed",
            "--gene_annotation", "myann/genes.bed"]
    with mock.patch.object(crisprme, "input_args", argv), \
         mock.patch.object(crisprme, "current_working_directory", work + os.sep), \
         mock.patch.object(crisprme, "_sort_bed", _sort_bed), \
         mock.patch.object(crisprme, "_compress_file", _compress_file), \
         mock.patch.object(crisprme.subprocess, "call", side_effect=call):
        exec(_dispatch("complete-search"), crisprme.__dict__)


class AnnotatedSearchRemovesItsTempFolders(unittest.TestCase):
    def setUp(self):
        # realpath: the folders are made under os.getcwd(), which resolves symlinks
        self.work = os.path.realpath(tempfile.mkdtemp(prefix="annot_tmpdirs_"))
        self.addCleanup(shutil.rmtree, self.work, True)
        cwd = os.getcwd()
        os.chdir(self.work)
        self.addCleanup(os.chdir, cwd)
        os.makedirs("Genomes/synth")
        with open("Genomes/synth/chrS.fa", "w") as fh:
            fh.write(">chrS\n" + "ACGTACGTACGTACGTACGTAGG" * 200 + "\n")
        os.makedirs("PAMs")
        with open("PAMs/20bp-NGG-SpCas9.txt", "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")
        with open("guide.txt", "w") as fh:
            fh.write("ACGTACGTACGTACGTACGTNNN\n")
        os.makedirs("myann")
        with open("myann/func.bed", "w") as fh:
            fh.write("chrS\t2000\t4600\tpromoter_test\nchrS\t0\t2000\texon_test\n")
        with open("myann/genes.bed", "w") as fh:
            fh.write("chrS\t100\t3000\texon:ENSTX.1:1\t.\t+\tENSEMBL\texon\t.\t"
                     "ID=exon:ENSTX.1:1;gene_id=ENSGX.1;gene_name=TESTG\n")
        with open("personal.bed", "w") as fh:
            fh.write("chrS\t10\t200\tmy_site\n")
        if hasattr(crisprme, "_ANNOTATION_TMPDIRS"):
            crisprme._ANNOTATION_TMPDIRS[:] = []
        env = {k: v for k, v in os.environ.items()
               if k not in ("CRISPRME_TMPDIR", "CRISPRME_KEEP_INTERMEDIATES")}
        env.update(CRISPRME_SKIP_SCORER_ENV="1", CRISPRME_SKIP_MEM_CHECK="1")
        for patcher in (mock.patch.dict(os.environ, env, clear=True),
                        mock.patch.object(crisprme, "current_working_directory",
                                          self.work + os.sep),
                        mock.patch.object(crisprme, "_sort_bed", _sort_bed),
                        mock.patch.object(crisprme, "_compress_file", _compress_file)):
            patcher.start()
            self.addCleanup(patcher.stop)

    def leftovers(self):
        return sorted(d for d in os.listdir(self.work)
                      if d.startswith(("crisprme_annot_", "crisprme_pannot_")))

    def search(self, *extra, pipeline_rc=0, output="ann"):
        """Run `crisprme.py complete-search ...` as the CLI dispatch does, the pipeline
        stubbed; return the annotation folders that existed while it ran. What it wrote to
        stderr is left in self.stderr, and whether it launched the pipeline in
        self.launched."""
        argv = ["crisprme.py", "complete-search", "--genome", "Genomes/synth",
                "--pam", "PAMs/20bp-NGG-SpCas9.txt", "--guide", "guide.txt",
                "--mm", "1", "--bDNA", "0", "--bRNA", "0", "--output", output,
                "--thread", "1", *extra]
        during = []
        real_call = subprocess.call
        self.launched = False

        def call(cmd, *a, **k):
            if isinstance(cmd, str) and "submit_job_automated_new_multiple_vcfs.sh" in cmd:
                self.launched = True
                # the pipeline is handed the sorted files in folders that still exist
                during.extend(d for d in self.leftovers()
                              if os.path.join(self.work, d) in cmd)
                if pipeline_rc is KeyboardInterrupt:  # Ctrl-C while it runs
                    raise KeyboardInterrupt
                return pipeline_rc
            return real_call(cmd, *a, **k)

        err = io.StringIO()
        try:
            with mock.patch.object(crisprme, "input_args", argv), \
                 mock.patch.object(crisprme.subprocess, "call", side_effect=call), \
                 contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(err):
                exec(_dispatch("complete-search"), crisprme.__dict__)
        finally:
            self.stderr = err.getvalue()
        return during

    ANNOTATED = ("--annotation", "myann/func.bed", "--gene_annotation", "myann/genes.bed")

    def test_a_completed_search_leaves_none(self):
        during = self.search(*self.ANNOTATED)
        self.assertEqual(len(during), 2, "the pipeline ran on both sorted annotations")
        self.assertTrue(all(d.startswith("crisprme_annot_") for d in during), during)
        self.assertEqual(self.leftovers(), [])

    def test_a_failed_pipeline_leaves_none(self):
        with self.assertRaises(SystemExit) as stop:  # complete-search: "CRISPRme run failed."
            self.search(*self.ANNOTATED, pipeline_rc=1)
        self.assertEqual(stop.exception.code, 1)
        self.assertEqual(self.leftovers(), [])

    def test_an_input_error_after_sorting_leaves_none(self):
        # both annotations are sorted before --max-total-edits is checked
        with self.assertRaises(SystemExit):
            self.search(*self.ANNOTATED, "--max-total-edits", "-1")
        self.assertEqual(self.leftovers(), [])

    def test_a_personal_annotation_folder_is_removed_too(self):
        during = self.search("--personal_annotation", "personal.bed")
        self.assertEqual([d[:16] for d in during], ["crisprme_pannot_"])
        self.assertEqual(self.leftovers(), [])

    def test_a_personal_annotation_is_added_to_the_annotation(self):
        # with --annotation too, the pipeline gets the sorted --annotation and the
        # personal intervals, tagged _personal, in one sorted file (kept here to read it)
        during = self.search("--annotation", "myann/func.bed",
                             "--personal_annotation", "personal.bed", "--keep-intermediates")
        self.assertEqual([d[:16] for d in during], ["crisprme_pannot_"])
        merged = os.path.join(self.work, during[0], "annotation+personal.bed.gz")
        with gzip.open(merged, "rt") as fh:
            self.assertEqual(fh.read().splitlines(), sorted([
                "chrS\t2000\t4600\tpromoter_test", "chrS\t0\t2000\texon_test",
                "chrS\t10\t200\tmy_site_personal"]))

    def merged_personal(self, personal, output):
        """The one sorted file the pipeline gets for --annotation with --personal_annotation
        ``personal``, as lines (the folder kept to read it, then removed)."""
        during = self.search("--annotation", "myann/func.bed", "--personal_annotation",
                             personal, "--keep-intermediates", output=output)
        self.assertEqual([d[:16] for d in during], ["crisprme_pannot_"])
        with gzip.open(os.path.join(self.work, during[0], "annotation+personal.bed.gz"),
                       "rt") as fh:
            lines = fh.read().splitlines()
        shutil.rmtree(os.path.join(self.work, during[0]))
        return lines

    def test_a_bgzipped_or_gzipped_personal_annotation_is_read_as_the_plain_one(self):
        # --help said the file "must be compressed using bgzip", but the tagging read it as
        # text: a bgzipped one reached sort-bed garbled ("No tabs/spaces found at line 1")
        # and the search stopped. A compressed file is told by its content, not its name.
        with open("personal.bed", "rb") as fh:
            plain = fh.read()
        expected = self.merged_personal("personal.bed", "plain")
        self.assertIn("chrS\t10\t200\tmy_site_personal", expected)
        for i, (name, data) in enumerate((("personal.bed.gz", _bgzf(plain)),
                                          ("personal.gzip.bed.gz", gzip.compress(plain)),
                                          ("personal_bgzipped_unsuffixed.bed", _bgzf(plain)))):
            with self.subTest(file=name):
                with open(name, "wb") as fh:
                    fh.write(data)
                self.assertEqual(self.merged_personal(name, f"compressed{i}"), expected)

    def test_a_damaged_compressed_personal_annotation_stops_the_search(self):
        # cut short, as an interrupted copy leaves it: one line that says so, no search
        with open("personal.bed", "rb") as fh:
            data = gzip.compress(fh.read())
        with open("personal.bed.gz", "wb") as fh:
            fh.write(data[:-12])
        with self.assertRaises(SystemExit):
            self.search("--personal_annotation", "personal.bed.gz")
        self.assertIn("Error: Decompressing personal annotation file failed", self.stderr)
        self.assertFalse(self.launched)
        self.assertEqual(self.leftovers(), [])

    def test_a_search_whose_scorer_env_is_not_ready_stops_before_the_pipeline(self):
        # --help, the README and TROUBLESHOOTING said such a run completes on CFD only;
        # complete-search stops, naming the fix and the opt-out (CRISPRME_SKIP_SCORER_ENV=1,
        # which every other test here sets). An env of a name no install has, so that the
        # release image, whose cbulge env is ready, gives the same answer as a bare runner.
        import scorer_env
        env = {k: v for k, v in os.environ.items() if k != "CRISPRME_SKIP_SCORER_ENV"}
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(scorer_env, "DEFAULT_ENV", "crisprme_test_no_such_env"):
            with self.assertRaises(SystemExit):
                self.search(*self.ANNOTATED)
        self.assertIn("'crisprme_test_no_such_env' CRISPR-Bulge scorer env is not ready",
                      self.stderr)
        self.assertIn("crisprme.py scorer-env create", self.stderr)
        self.assertIn("set CRISPRME_SKIP_SCORER_ENV=1", self.stderr)
        self.assertFalse(self.launched)
        self.assertEqual(self.leftovers(), [])

    def test_a_personal_annotation_with_an_annotation_leaves_none(self):
        self.search("--annotation", "myann/func.bed", "--personal_annotation", "personal.bed")
        self.assertEqual(self.leftovers(), [])

    def test_another_jobs_folder_is_left_alone(self):
        other = os.path.join(self.work, "crisprme_annot_otherjob")
        os.makedirs(other)
        with open(os.path.join(other, "annotation.sorted.bed.gz"), "wb") as fh:
            fh.write(b"in use")
        self.search(*self.ANNOTATED)
        self.assertEqual(self.leftovers(), ["crisprme_annot_otherjob"])

    def test_keep_intermediates_keeps_them(self):
        during = self.search(*self.ANNOTATED, "--keep-intermediates")
        self.assertEqual(len(during), 2)
        self.assertEqual(self.leftovers(), sorted(during))

    def test_a_ctrl_c_while_the_pipeline_runs_keeps_them(self):
        # subprocess.call has killed only the pipeline's /bin/sh: the pipeline may go on
        with self.assertRaises(KeyboardInterrupt):
            self.search(*self.ANNOTATED, pipeline_rc=KeyboardInterrupt)
        self.assertEqual(len(self.leftovers()), 2)

    def test_a_pipeline_ended_by_a_signal_keeps_them(self):
        # 128 + N: its bash was killed by signal N, its workers need not have been
        # (an INT, TERM or HUP trap ends it at once in a `wait`); -N: its /bin/sh was
        for i, rc in enumerate((130, 137, 143, -15)):
            with self.subTest(rc=rc):
                for d in self.leftovers():
                    shutil.rmtree(os.path.join(self.work, d))
                with self.assertRaises(SystemExit) as stop:  # "CRISPRme run failed."
                    self.search(*self.ANNOTATED, pipeline_rc=rc, output=f"sig{i}")
                self.assertEqual(stop.exception.code, 1)
                self.assertEqual(len(self.leftovers()), 2)

    def test_a_sigint_to_crisprme_py_alone_leaves_the_running_pipeline_its_folders(self):
        # kill -INT <crisprme.py pid>, not to its process group: the pipeline, here a
        # worker of it that outlives its /bin/sh, goes on after crisprme.py has exited,
        # and still finds both sorted annotations
        sdir = os.path.join(self.work, "stub")
        os.makedirs(sdir)
        stub = os.path.join(sdir, "pipeline.sh")
        with open(stub, "w") as fh:
            fh.write(_STUB_PIPELINE)
        os.chmod(stub, 0o755)
        go, started, result, log = (os.path.join(sdir, f) for f in
                                    ("go", "started", "worker.result", "child.log"))
        with open(log, "wb") as fh:
            child = subprocess.Popen(
                [sys.executable, "-c", "import sys; sys.path.insert(0, sys.argv[1]); "
                 "import test_annotation_tmpdirs as t; t._search_until_interrupted("
                 "*sys.argv[2:])", HERE, self.work, stub],
                cwd=self.work, env=dict(os.environ), stdout=fh, stderr=subprocess.STDOUT)

        def child_log():
            with open(log, errors="replace") as fh:
                return fh.read()

        def until(done, what):
            deadline = time.monotonic() + 120
            while not done():
                if time.monotonic() > deadline:
                    self.fail(f"no {what} after 120 s:\n{child_log()}")
                time.sleep(0.05)

        def stop():  # never leave the child or the worker running
            open(go, "w").close()
            if child.poll() is None:
                child.kill()
            child.wait()
            deadline = time.monotonic() + 30
            while (os.path.exists(started) and not os.path.exists(result)
                   and time.monotonic() < deadline):
                time.sleep(0.05)
        self.addCleanup(stop)

        until(lambda: os.path.exists(started) or child.poll() is not None,
              "pipeline started")
        self.assertIsNone(child.poll(), child_log())
        os.kill(child.pid, signal.SIGINT)
        child.wait(timeout=120)
        self.assertIn("KeyboardInterrupt", child_log())
        self.assertEqual(len(self.leftovers()), 2, "the running pipeline's folders")
        open(go, "w").close()
        until(lambda: os.path.exists(result), "worker result")
        with open(result) as fh:
            self.assertEqual(fh.read().strip(), "read")


class TheCosmicFreeCopyOfAnAnnotationItsFolderCannotHold(unittest.TestCase):
    """When the folder of the annotation cannot hold its COSMIC-free copy (a read-only or a
    shared install), _apply_cosmic_license strips it elsewhere: into $CRISPRME_TMPDIR, the
    working directory, or tempfile's folder (/tmp). It wrote there under a fixed name,
    .nocosmic.<name>, the strip beside it under one made of its process id, and in a shared
    /tmp any other user of the machine can make such names first: a link to a file of the
    user's, which the strip overwrote with the annotation's lines (2.10.0's sweep after the
    review of lane S, round 3). It strips into a folder of its own there now (mkdtemp: 0700, a
    new name), removed with the annotation temp folders."""

    def setUp(self):
        self.work = os.path.realpath(tempfile.mkdtemp(prefix="nocosmic_"))
        self.addCleanup(shutil.rmtree, self.work, True)
        os.makedirs(os.path.join(self.work, "Annotations"))  # no licence attested there
        self.ann_dir = os.path.join(self.work, "install")
        os.makedirs(self.ann_dir)
        self.ann = os.path.join(self.ann_dir, "func.bed")
        with open(self.ann, "w") as fh:
            fh.write("chrS\t0\t10\tTP53_COSMIC\nchrS\t20\t30\tpromoter_test\n")
        self.shared = os.path.join(self.work, "shared")  # a /tmp other users write
        os.makedirs(self.shared)
        crisprme._ANNOTATION_TMPDIRS[:] = []
        self.addCleanup(crisprme._ANNOTATION_TMPDIRS.clear)
        import cosmic_license
        real, ann_dir = cosmic_license.strip_cosmic, self.ann_dir

        def strip(src, out):
            # the annotation's folder cannot be written, whoever runs the test (root too)
            if os.path.dirname(out) == ann_dir:
                raise PermissionError(13, "Permission denied", out)
            return real(src, out)
        env = {k: v for k, v in os.environ.items() if k != "CRISPRME_TMPDIR"}
        env["CRISPRME_TMPDIR"] = self.shared
        for patcher in (mock.patch.dict(os.environ, env, clear=True),
                        mock.patch.object(crisprme, "current_working_directory",
                                          self.work + os.sep),
                        mock.patch.object(cosmic_license, "strip_cosmic", strip)):
            patcher.start()
            self.addCleanup(patcher.stop)

    def test_it_strips_into_a_folder_of_its_own(self):
        victim = os.path.join(self.work, "victim.txt")
        with open(victim, "w") as fh:
            fh.write("mine\n")
        fixed = os.path.join(self.shared, ".nocosmic.func.bed")
        link = fixed + ".%d.tmp" % os.getpid()
        os.symlink(victim, link)  # another user's, made first
        with open(fixed, "w") as fh:
            fh.write("chrS\t0\t99\tplanted_row\n")
        with contextlib.redirect_stdout(io.StringIO()):
            out = crisprme._apply_cosmic_license([], self.ann)
        with open(victim) as fh:
            self.assertEqual(fh.read(), "mine\n")
        folder = os.path.dirname(out)
        self.assertEqual(os.path.dirname(folder), self.shared)
        self.assertTrue(os.path.basename(folder).startswith("crisprme_nocosmic_"), out)
        self.assertEqual(os.stat(folder).st_mode & 0o777, 0o700)
        self.assertEqual(crisprme._ANNOTATION_TMPDIRS, [folder])
        with open(out) as fh:
            self.assertEqual(fh.read(), "chrS\t20\t30\tpromoter_test\n")
        # what the other user made is left as it was
        self.assertTrue(os.path.islink(link))
        with open(fixed) as fh:
            self.assertEqual(fh.read(), "chrS\t0\t99\tplanted_row\n")


class TheDispatchRemovesThemForEveryCallerThatSortsAnnotations(unittest.TestCase):
    """complete-search and assembly-search are the two subcommands that call
    _check_annotation / _check_gene_annotation (the web's assembly-search job runs
    `crisprme.py assembly-search` too; its complete-search job prepares its annotation
    in the job folder, not in a crisprme_annot_* one)."""

    def test_both_run_through_the_remover(self):
        self.assertEqual(_dispatch("complete-search").strip(),
                         "_run_and_remove_annotation_tmpdirs(complete_search)")
        self.assertEqual(_dispatch("assembly-search").strip(),
                         "_run_and_remove_annotation_tmpdirs(assembly_search)")

    @staticmethod
    def _callers(call_pattern):
        """Top-level functions of crisprme.py whose code calls `call_pattern`."""
        with open(CRISPRME_PY) as fh:
            src = fh.read().split(_MARKER, 1)[0]
        callers = set()
        for m in re.finditer(r"^def (\w+)\(", src, re.M):
            body = src[m.end():].split("\ndef ", 1)[0]
            code = "\n".join(ln.split("#", 1)[0] for ln in body.splitlines())
            if re.search(call_pattern, code):
                callers.add(m.group(1))
        return callers

    def test_no_other_subcommand_makes_them(self):
        self.assertEqual(self._callers(r"\b_sort_annotation\("),
                         {"_check_annotation", "_check_gene_annotation"})
        self.assertEqual(self._callers(r"\b_process_personal_annotation\("),
                         {"_check_personal_annotation"})
        self.assertEqual(
            self._callers(r"\b_check_(gene_|personal_)?annotation\("),
            {"complete_search", "assembly_search"})


class TheHelpSaysWhatKeepIntermediatesKeeps(unittest.TestCase):
    """complete-search --help is where a user looks for what --keep-intermediates keeps.
    It listed the pipeline's working tables only, so nothing but the CHANGELOG said that
    the flag also keeps the annotation temp folders, and it named the bestMerge table,
    which the cleanup never removes, among the files removed."""

    @staticmethod
    def _entry():
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            try:
                crisprme.print_help_complete_search()
            except SystemExit:
                pass
        return next(ln for ln in err.getvalue().split("\n")
                    if ln.lstrip().startswith("--keep-intermediates,"))

    def test_it_names_both_folders_and_where_they_are(self):
        entry = self._entry()
        for name in ("crisprme_annot_*", "crisprme_pannot_*", "$CRISPRME_TMPDIR"):
            self.assertIn(name, entry)

    def test_it_lists_the_best_merge_table_among_the_kept_outputs(self):
        removed, kept = self._entry().split("kept either way")[0].split("canonical outputs")
        self.assertNotIn("bestMerge", removed)
        self.assertIn("bestMerge", kept)


class TheHelpSaysWhichPersonalAnnotationFilesWork(unittest.TestCase):
    """--help is the only place that says what --personal_annotation reads, and it said the
    file "must be compressed using bgzip": since 2.2.0 only a plain one had worked, and a
    bgzipped one stopped the search. Both work now, and gzip too."""

    def test_it_says_plain_bgzip_and_gzip_work(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            try:
                crisprme.print_help_complete_search()
            except SystemExit:
                pass
        entry = next(ln for ln in err.getvalue().split("\n")
                     if ln.lstrip().startswith("--personal_annotation,"))
        self.assertIn("may be plain text or compressed with bgzip or gzip", entry)
        self.assertNotIn("must be compressed", entry)


if __name__ == "__main__":
    unittest.main(verbosity=2)
