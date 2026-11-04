"""Text a request sends to the web app runs no command, and '../' in it reaches nothing
outside the data folder.

Up to 2.8.0 the web app took such text into shell commands as it was:

* the genomic-sequence input with BED coordinates (seq_script/extract_seq.py, run by the
  Submit check, the base-editing window list and the search launch): the sequence's '>'
  header named the BED it wrote in the data folder, and the BED line's chromosome the
  FASTA bedtools read, in a command run through a shell; a header holding '../' wrote the
  BED outside the data folder, and the index bedtools wrote beside the FASTA was then
  removed, one there before too.

Each case sends a value holding a command that would create a marker file, and a '../'
path, and checks that no marker and no file outside the data folder appear, and that the
normal value still works. The extraction's checks run everywhere (those that read a
sequence need bedtools); the pages' need the web stack (Dash), and run in the image.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_shell_safety.py' -v
"""

import os
import shutil
import sys
import tempfile
import unittest
import warnings
from unittest.mock import MagicMock, patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from seq_script import extract_seq  # noqa: E402  (no Dash, no Biopython)

try:  # the pages need the web stack (Dash), which the light unit-test CI does not have
    import dash  # noqa: F401
except ImportError:  # pragma: no cover
    dash = None
if dash is not None:
    sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object
    import index  # noqa: E402  (check_directories: the data folder the web app starts with)
    import pages.load_page as lp  # noqa: E402
    import pages.main_page as mp  # noqa: E402
    import pages.pages_utils as pu  # noqa: E402
    import pages.results_page as rp  # noqa: E402

BEDTOOLS = shutil.which("bedtools")
NEEDS_DASH = unittest.skipUnless(dash, "dash not installed: web page tests skipped")
NEEDS_BEDTOOLS = unittest.skipUnless(BEDTOOLS, "bedtools not on PATH")
PAM = "20bp-NRG-SpCas9"
# chr1 of the test genome: one 20-nt guide and its NGG, between T runs
SEQ = "T" * 10 + "CTAACAGTTGCTTTTATCACAGG" + "T" * 10


class _Folder(unittest.TestCase):
    """The web app's data folder (data/, its working directory) with genome hg38 (chr1)
    and the NRG PAM, beside a folder outside it (outside/) that '../' reaches, holding a
    genome and a PAM of its own."""

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)  # files the page code leaves open
        self.top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.top, ignore_errors=True)
        self.data = os.path.join(self.top, "data")
        self.outside = os.path.join(self.top, "outside")
        for d in (os.path.join(self.data, "Genomes"), os.path.join(self.data, "PAMs"),
                  self.outside):
            os.makedirs(d, exist_ok=True)
        for folder in (os.path.join(self.data, "Genomes", "hg38"),
                       os.path.join(self.outside, "hg38")):
            os.makedirs(folder)
            with open(os.path.join(folder, "chr1.fa"), "w") as fh:
                fh.write(">chr1\n" + SEQ + "\n")
        for folder in (os.path.join(self.data, "PAMs"), self.outside):
            with open(os.path.join(folder, f"{PAM}.txt"), "w") as fh:
                fh.write("N" * 20 + "NRG 3\n")
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.data)  # where the web app runs, and a command it starts
        self.before = self.outside_files()

    def markers(self):
        """The marker files a command in the values sent would have created."""
        return sorted(os.path.relpath(os.path.join(root, f), self.top)
                      for root, _dirs, files in os.walk(self.top)
                      for f in files if f.startswith("MARK"))

    def outside_files(self):
        """Every file and folder of the test's folder outside the data folder."""
        found = set()
        for root, dirs, files in os.walk(self.top):
            dirs[:] = [d for d in dirs if os.path.join(root, d) != self.data]
            found |= {os.path.relpath(os.path.join(root, n), self.top) for n in dirs + files}
        return found

    def assertNothingRun(self):
        self.assertEqual(self.markers(), [])
        self.assertEqual(self.outside_files(), self.before)


class TestExtractSequence(_Folder):
    """seq_script/extract_seq.py, which the genomic-sequence input's BED coordinates go
    through."""

    def test_a_header_holding_a_command_runs_nothing(self):
        # the header named the BED, in the shell command: it ran with no FASTA too
        for name in ("s;touch${IFS}MARK_a;#", "s$(touch${IFS}MARK_b)", "s`touch${IFS}MARK_c`"):
            with self.subTest(name=name):
                with self.assertRaises(FileNotFoundError):
                    extract_seq.extractSequence(name, "chrZ:0-43", "hg38")
                self.assertNothingRun()
                self.assertEqual(sorted(os.listdir(self.data)), ["Genomes", "PAMs"])

    def test_a_header_holding_dot_dot_writes_nothing_outside(self):
        # its BED, written as ../outside/x.bed, stayed there when bedtools failed
        with self.assertRaises(FileNotFoundError):
            extract_seq.extractSequence("../outside/x", "chrZ:0-43", "hg38")
        self.assertNothingRun()

    def test_only_a_fasta_of_an_installed_genome_is_read(self):
        for chrom, genome in (("chr1$(touch${IFS}MARK_d)", "hg38"),
                              ("chr1;touch${IFS}MARK_e;#", "hg38"),
                              ("../../../outside/hg38/chr1", "hg38"),
                              ("chr1", "../../outside/hg38"),
                              ("chr1", ".")):
            with self.subTest(chrom=chrom, genome=genome):
                with self.assertRaises(FileNotFoundError):
                    extract_seq.extractSequence("s", f"{chrom}:0-43", genome)
                self.assertNothingRun()  # bedtools indexed ../outside's FASTA, there

    @NEEDS_BEDTOOLS
    def test_the_sequence_is_read_as_before_and_no_index_is_left_or_removed(self):
        genome = os.path.join(self.data, "Genomes", "hg38")
        self.assertEqual(extract_seq.extractSequence("s", "chr1:0-43", "hg38"), SEQ)
        # the header is not a file name now: one with a command reads the same
        self.assertEqual(
            extract_seq.extractSequence("s;touch${IFS}MARK_f;#", "chr1:10-33", "hg38"),
            SEQ[10:33])
        self.assertNothingRun()
        self.assertEqual(sorted(os.listdir(genome)), ["chr1.fa"])
        self.assertEqual(sorted(os.listdir(self.data)), ["Genomes", "PAMs"])
        # an index that was there stays (it was removed), and a .fasta's is not left
        with open(os.path.join(genome, "chr1.fa.fai"), "w") as fh:
            fh.write(f"chr1\t{len(SEQ)}\t6\t{len(SEQ)}\t{len(SEQ) + 1}\n")
        self.assertEqual(extract_seq.extractSequence("s", "chr1:0-43", "hg38"), SEQ)
        self.assertTrue(os.path.isfile(os.path.join(genome, "chr1.fa.fai")))
        os.makedirs(os.path.join(self.data, "Genomes", "mm"))
        with open(os.path.join(self.data, "Genomes", "mm", "chr2.fasta"), "w") as fh:
            fh.write(">chr2\n" + SEQ + "\n")
        self.assertEqual(extract_seq.extractSequence("s", "chr2:0-43", "mm"), SEQ)
        self.assertEqual(os.listdir(os.path.join(self.data, "Genomes", "mm")), ["chr2.fasta"])


class _WebFolder(_Folder):
    """The data folder as the web app makes it (index.check_directories), the pages
    reading it, and stand-ins for the scripts the pages start (in an app folder's
    PostProcess/)."""

    def setUp(self):
        super().setUp()
        index.check_directories(self.data)
        self.app = tempfile.mkdtemp()  # (the app's folder: not the data folder's side)
        self.addCleanup(shutil.rmtree, self.app, ignore_errors=True)
        os.makedirs(os.path.join(self.app, "PostProcess"))
        for name, text in (
            # the search, which records its arguments
            ("submit_job_automated_new_multiple_vcfs.sh",
             '#!/bin/bash\nprintf "%s\\n" "$@" > "$(dirname "$0")/args.txt"\n'),
        ):
            path = os.path.join(self.app, "PostProcess", name)
            with open(path, "w") as fh:
                fh.write(text)
            os.chmod(path, 0o755)
        for mod in (mp, pu, lp, rp):
            p = patch.object(mod, "current_working_directory", self.data + os.sep)
            p.start()
            self.addCleanup(p.stop)
        for mod in (mp, rp):
            p = patch.object(mod, "app_directory", self.app)
            p.start()
            self.addCleanup(p.stop)
        self.before = self.outside_files()


@NEEDS_DASH
class TestSubmitCheckAndBaseEditing(_WebFolder):
    """The genomic-sequence input with BED coordinates, as Submit checks it and as the
    base-editing window list reads it at each change."""

    def submit(self, text):
        out = mp.check_input(1, None, "hg38", PAM, "GS", text, 4, 0, 0, "ref", False)
        return out[0], [li.children for li in out[-1].children[1].children]

    @NEEDS_BEDTOOLS
    def test_submit_check(self):
        self.assertEqual(self.submit(">s\nchr1 0 43"), (1, []))  # read, its guide searched
        self.assertEqual(self.submit(">s;touch${IFS}MARK_a;#\nchr1 0 43"), (1, []))
        for text in (">s\nchr1$(touch${IFS}MARK_b) 0 43", ">../outside/x\nchrZ 0 43"):
            with self.subTest(text=text):
                submit, lines = self.submit(text)
                self.assertIsNone(submit)
                self.assertTrue(any(line.startswith("Cannot read chr") for line in lines))
        self.assertNothingRun()

    @NEEDS_BEDTOOLS
    def test_base_editing_window_list(self):
        for text in (">s\nchr1 0 43", ">s;touch${IFS}MARK_c;#\nchr1 0 43"):
            with self.subTest(text=text):
                options, _ = mp.update_base_editing_dropdown(text, "GS", "hg38")
                self.assertEqual([o["value"] for o in options], list(range(1, len(SEQ) + 1)))
        for text in (">s\nchr1$(touch${IFS}MARK_d) 0 43", ">../outside/x\nchrZ 0 43"):
            with self.subTest(text=text), self.assertRaises(FileNotFoundError):
                mp.update_base_editing_dropdown(text, "GS", "hg38")
        self.assertNothingRun()


@NEEDS_DASH
class TestSearchLaunch(_WebFolder):
    """Submit's launch of the search (change_url), which a request can call with any
    values: the job pool is a stand-in, whose command runs as the pool runs it."""

    def setUp(self):
        super().setUp()
        self.pool = MagicMock()
        p = patch.object(mp, "pool_executor", self.pool)
        p.start()
        self.addCleanup(p.stop)

    def submit(self, genome="hg38", pam=PAM, variant="ref", guides="CTAACAGTTGCTTTTATCAC",
               guide_type="IP", mms=4, email=None, name=None):
        """Submit, as the form sends it (Advanced options open), and the search it ran."""
        mp.change_url(1, "http://127.0.0.1:8080/index", genome, variant, pam, guide_type,
                      guides, mms, 1, 1, "N", None, None, None, ["email"] if email else [],
                      email, name, 6, True, "population-level", [])
        self.assertEqual(self.pool.submit.call_count, 1)
        fn, *args = self.pool.submit.call_args[0]
        cmd = args[0]
        fn(*args, **self.pool.submit.call_args[1])
        self.pool.reset_mock()
        with open(os.path.join(self.app, "PostProcess", "args.txt")) as fh:
            return cmd, fh.read().split("\n")[:-1]

    def jobs(self):
        return sorted(os.listdir(os.path.join(self.data, "Results")))

    @NEEDS_BEDTOOLS
    def test_a_genomic_sequence_header_runs_nothing(self):
        self.submit(guides=">s;touch${IFS}MARK_b;#\nchr1 0 43", guide_type="GS")
        job = os.path.join(self.data, "Results", self.jobs()[0])
        with open(os.path.join(job, ".guides.txt")) as fh:
            self.assertIn("CTAACAGTTGCTTTTATCACNNN", fh.read().split("\n"))
        self.assertNothingRun()


if __name__ == "__main__":
    unittest.main()
