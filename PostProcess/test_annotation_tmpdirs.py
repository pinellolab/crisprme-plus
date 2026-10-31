#!/usr/bin/env python3
"""complete-search and assembly-search remove the annotation temp folders they make.

_sort_annotation (crisprme.py) sorts each --annotation / --gene_annotation into a new
crisprme_annot_* folder in the working directory ($CRISPRME_TMPDIR if set), and
_process_personal_annotation builds --personal_annotation in a crisprme_pannot_* one.
The pipeline reads those files, and gunzips back into them, until it ends. Nothing
removed them, so every annotated search left two behind, ~1.7 GB with the shipped
annotations.

These tests run complete-search the way the CLI does -- the dispatch statement is taken
verbatim from crisprme.py, whose definitions are exec'd without the dispatch (the
truncate-and-exec loader of test_bulge_cap) -- on a tiny synthetic genome, with the
pipeline launch stubbed to succeed or fail, and sort-bed / bgzip replaced by Python
stand-ins (the light CI env has neither). Every other step of complete-search runs.
"""
import contextlib
import gzip
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import types
import unittest
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
        stubbed; return the annotation folders that existed while it ran."""
        argv = ["crisprme.py", "complete-search", "--genome", "Genomes/synth",
                "--pam", "PAMs/20bp-NGG-SpCas9.txt", "--guide", "guide.txt",
                "--mm", "1", "--bDNA", "0", "--bRNA", "0", "--output", output,
                "--thread", "1", *extra]
        during = []
        real_call = subprocess.call

        def call(cmd, *a, **k):
            if isinstance(cmd, str) and "submit_job_automated_new_multiple_vcfs.sh" in cmd:
                # the pipeline is handed the sorted files in folders that still exist
                during.extend(d for d in self.leftovers()
                              if os.path.join(self.work, d) in cmd)
                return pipeline_rc
            return real_call(cmd, *a, **k)

        with mock.patch.object(crisprme, "input_args", argv), \
             mock.patch.object(crisprme.subprocess, "call", side_effect=call), \
             contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            exec(_dispatch("complete-search"), crisprme.__dict__)
        return during

    ANNOTATED = ("--annotation", "myann/func.bed", "--gene_annotation", "myann/genes.bed")

    def test_a_completed_search_leaves_none(self):
        during = self.search(*self.ANNOTATED)
        self.assertEqual(len(during), 2, "the pipeline ran on both sorted annotations")
        self.assertTrue(all(d.startswith("crisprme_annot_") for d in during), during)
        self.assertEqual(self.leftovers(), [])

    def test_a_failed_pipeline_leaves_none(self):
        with self.assertRaises(OSError):  # complete-search: "CRISPRme run failed!"
            self.search(*self.ANNOTATED, pipeline_rc=1)
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
