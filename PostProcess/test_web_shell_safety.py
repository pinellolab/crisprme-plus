"""Text a request sends to the web app runs no command, and '../' in it reaches nothing
outside the data folder.

Up to 2.8.0 the web app took such text into shell commands as it was:

* the genomic-sequence input with BED coordinates (seq_script/extract_seq.py, run by the
  Submit check, the base-editing window list and the search launch): the sequence's '>'
  header named the BED it wrote in the data folder, and the BED line's chromosome the
  FASTA bedtools read, in a command run through a shell; a header holding '../' wrote the
  BED outside the data folder, and the index bedtools wrote beside the FASTA was then
  removed, one there before too;
* the search launch (main_page.change_url): the e-mail address, and the mismatches and
  genome, which a request can set past the Submit check;
* the Personal assembly tab's search launch (main_page.submit_assembly_search_job): the
  guides, which it wrote to the search's guide file as sent, and the search passes each
  guide of that file to a shell command (process_summaries.py); and the PAM and the six
  folders and files of the individual, which it searched as sent, '../' too;
* the active annotation of a genome (pages_utils.build_active_annotation), whose file is
  named after the genome, a name Settings took as typed (test_settings_progress.py checks
  that it no longer takes one a shell reads; a genome folder can be put there by hand);
* Remove result on a job's status page (load_page.remove_result): the job id of the URL;
* the results pages: the guide, sample, chromosome and position of a results link, the
  guide and mismatches the tables send back, a sample's name (a VCF's), and the results
  file's name, which holds the genome's.

pages/sequence_utils.py, which no module imports, ran bedtools getfasta through a shell
the same way; it is checked too (with Biopython and bedtools), for a page that would.

And up to 2.8.0 the pages of a job took the job id of their URL (/load?job=,
/result?job=) into paths as it was, with no shell: '../' in it reached a folder outside
Results/, outside the data folder too. Remove result of a personal-assembly search
removed that folder (and the haplotype folders its .Params.txt names, joined the same
way), Show full log and the status of a personal-assembly search read its log files,
the results pages wrote the table filter they keep (.dropdown.json) into it, and their
tables, lists, download links and the routes to a job read its files. The callbacks that
fill in and check the search form (its variant datasets, analysis modes, bulge options and
the note under them, the Submit check) joined the genome, PAM and dataset a request sends
into the paths they list and read the same way, and what they found outside the data
folder changed what they answered.

Each case sends a value holding a command that would create a marker file, and a '../'
path, and checks that no marker and no file outside the data folder appear, and that the
normal value still works; a job id holding '../' reads, writes and removes nothing
outside Results/. The extraction's checks run everywhere (those that read a sequence
need bedtools); the pages' need the web stack (Dash), and run in the image.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_shell_safety.py' -v
"""

import builtins
import contextlib
import gzip
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import unittest
import warnings
from unittest.mock import MagicMock, patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from seq_script import extract_seq  # noqa: E402  (no Dash, no Biopython)

try:  # pages/sequence_utils.py needs Biopython, which the light unit-test CI does not have
    from pages import sequence_utils  # noqa: E402
except ImportError:  # pragma: no cover
    sequence_utils = None
try:  # the pages need the web stack (Dash), which the light unit-test CI does not have
    import dash  # noqa: F401
except ImportError:  # pragma: no cover
    dash = None
if dash is not None:
    from dash.exceptions import PreventUpdate
    import pandas as pd

    sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object
    import index  # noqa: E402  (check_directories: the data folder the web app starts with)
    import pages.load_page as lp  # noqa: E402
    import pages.main_page as mp  # noqa: E402
    import pages.pages_utils as pu  # noqa: E402
    import pages.results_page as rp  # noqa: E402

BEDTOOLS = shutil.which("bedtools")
NEEDS_DASH = unittest.skipUnless(dash, "dash not installed: web page tests skipped")
NEEDS_BEDTOOLS = unittest.skipUnless(BEDTOOLS, "bedtools not on PATH")
NEEDS_BIOPYTHON = unittest.skipUnless(sequence_utils, "Biopython not installed")
PAM = "20bp-NRG-SpCas9"
# chr1 of the test genome: one 20-nt guide and its NGG, between T runs
SEQ = "T" * 10 + "CTAACAGTTGCTTTTATCACAGG" + "T" * 10
GUIDE = "CTAACAGTTGCTTTTATCACNNN"
CRITERIA = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


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


@NEEDS_BIOPYTHON
@NEEDS_BEDTOOLS
class TestSequenceUtils(_Folder):
    """pages/sequence_utils.py, which no module imports: its bedtools getfasta ran
    through a shell on a BED named after the header and the FASTA of the genome's
    folder, so a page that imported it would have run what either holds."""

    def test_the_sequence_is_read_as_before_and_runs_nothing(self):
        genome = os.path.join(self.data, "Genomes", "hg38")
        self.assertEqual(sequence_utils.extract_sequence("s", "chr1:0-43", "hg38"), SEQ)
        for name in ("s;touch${IFS}MARK_a;#", "s$(touch${IFS}MARK_b)", "s`touch${IFS}MARK_c`"):
            with self.subTest(name=name):
                self.assertEqual(
                    sequence_utils.extract_sequence(name, "chr1:10-33", "hg38"), SEQ[10:33])
        os.rename(genome, os.path.join(self.data, "Genomes", "g;touch${IFS}MARK_d;#"))
        self.assertEqual(
            sequence_utils.extract_sequence("s", "chr1:0-43", "g;touch${IFS}MARK_d;#"), SEQ)
        self.assertNothingRun()
        # as before: the BED and the FASTA's index are removed
        self.assertEqual(sorted(os.listdir(self.data)), ["Genomes", "PAMs"])
        self.assertEqual(
            os.listdir(os.path.join(self.data, "Genomes", "g;touch${IFS}MARK_d;#")), ["chr1.fa"])


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
            # the Graphical Reports radar chart, which records its arguments
            ("generate_img_radar_chart.py",
             "import json, os, sys\nwith open(os.path.join(os.path.dirname(__file__), "
             "'radar_args.json'), 'w') as fh:\n    json.dump(sys.argv[1:], fh)\n"),
            # the Personal Risk Card plots, which draws empty ones
            ("CRISPRme_plots_personal.py",
             "import os, sys\nfor kind in ('fewest', 'CFD', 'CRISPR_Bulge'):\n"
             "    open(os.path.join(sys.argv[2], 'CRISPRme_%s_top_1000_log_for_main_text_%s"
             ".png' % (kind, sys.argv[3])), 'wb').close()\n"),
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

    def test_a_search_runs_as_before_and_the_email_address_is_one_argument(self):
        cmd, args = self.submit(name="run1")
        # nothing is quoted: the command is the one of 2.8.0
        self.assertNotIn("'", cmd)
        job = os.path.join(self.data, "Results", self.jobs()[0])
        self.assertEqual(args[:3], [os.path.join(self.data, "Genomes", "hg38"),
                                    os.path.join(job, ".list_vcfs.txt"),
                                    os.path.join(job, ".guides.txt")])
        self.assertEqual((len(args), args[7], args[16], args[25]), (26, "4", "_", "6"))
        email = "a@b;touch${IFS}MARK_a;#"
        _cmd, args = self.submit(email=email, name="run2")
        self.assertEqual(len(args), 26)
        self.assertEqual(args[16], email)  # the search's 17th argument
        self.assertNothingRun()

    @NEEDS_BEDTOOLS
    def test_a_genomic_sequence_header_runs_nothing(self):
        self.submit(guides=">s;touch${IFS}MARK_b;#\nchr1 0 43", guide_type="GS")
        job = os.path.join(self.data, "Results", self.jobs()[0])
        with open(os.path.join(job, ".guides.txt")) as fh:
            self.assertIn("CTAACAGTTGCTTTTATCACNNN", fh.read().split("\n"))
        self.assertNothingRun()

    def test_what_the_form_does_not_offer_is_refused(self):
        for kwargs in (dict(mms="4;touch${IFS}MARK_c;#"), dict(mms="x"),
                       dict(genome="hg38;touch${IFS}MARK_d;#"),
                       dict(genome="../../outside/hg38"), dict(pam=f"../../outside/{PAM}"),
                       dict(variant="1000G$(touch${IFS}MARK_e)"), dict(genome=None)):
            with self.subTest(**kwargs):
                with self.assertRaises(PreventUpdate):
                    self.submit(**kwargs)
                self.assertEqual(self.pool.submit.call_count, 0)
                self.assertEqual(self.jobs(), [])
        self.assertNothingRun()

    @unittest.skipUnless(shutil.which("sort-bed") and shutil.which("bgzip"),
                         "BEDOPS sort-bed or bgzip not on PATH")
    def test_a_genome_named_with_a_command_runs_nothing(self):
        # Settings took a genome's name as typed, and a genome folder can be put in Genomes/
        # by hand; its active annotation is named after it
        genome = "x$(touch${IFS}MARK_f)"
        os.makedirs(os.path.join(self.data, "Genomes", genome))
        shutil.copy(os.path.join(self.data, "Genomes", "hg38", "chr1.fa"),
                    os.path.join(self.data, "Genomes", genome))
        for name, row in (("a.bed", "chr1\t10\t500\tpromoter\n"),
                          ("b.bed", "chr1\t2\t9\tenhancer\n")):
            with open(os.path.join(self.data, "Annotations", name), "w") as fh:
                fh.write(row)
        pu.write_enabled_annotations(genome, ["a.bed", "b.bed"])
        # (the form offers a genome folder with its underscores as spaces)
        _cmd, args = self.submit(genome=genome.replace("_", " "))
        self.assertEqual(args[0], os.path.join(self.data, "Genomes", genome))
        self.assertTrue(args[4].endswith(f".active.{genome}.bed.gz"))
        self.assertTrue(os.path.isfile(args[4]))
        self.assertNothingRun()


@contextlib.contextmanager
def _paths_read():
    """A list of the real path of every path listed, tested or opened under it
    (os.listdir, os.scandir, os.path.isdir / isfile / exists, open, gzip.open)."""
    seen = []

    def recording(func):
        def recorded(*args, **kwargs):
            if args and isinstance(args[0], (str, bytes, os.PathLike)):
                seen.append(os.fsdecode(os.path.realpath(args[0])))
            return func(*args, **kwargs)
        return recorded

    with contextlib.ExitStack() as stack:
        for owner, name in ((os, "listdir"), (os, "scandir"), (os.path, "isdir"),
                            (os.path, "isfile"), (os.path, "exists"), (builtins, "open"),
                            (gzip, "open")):
            stack.enter_context(patch.object(owner, name, recording(getattr(owner, name))))
        yield seen


def _vcf(path, samples):
    """A gzipped VCF of a header only: with a sample column, or none (sites-only)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                 + ("\tFORMAT\tS1" if samples else "") + "\n")


@NEEDS_DASH
class TestSearchFormReadsItsOwnValues(_WebFolder):
    """The callbacks that fill in and check the search form, which a request can call with
    any genome, PAM and variant dataset: they joined them into the paths of the folders and
    files they list and read, so an absolute path or '../' reached outside the data folder,
    and what was there changed what they answered (a bit about any folder or file of the
    server a request). Values the form does not offer are answered as none is, and nothing
    outside is read; those it offers are answered as before."""

    def setUp(self):
        super().setUp()
        # an installed dataset: a VCF with samples, its sample list and genotype store
        _vcf(os.path.join(self.data, "VCFs", "hg38_walk", "walk.chr1.vcf.gz"), True)
        with open(os.path.join(self.data, "samplesIDs", "hg38_walk.samplesID.txt"), "w") as fh:
            fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                     "S1\twalk\tunknown\tunknown\n")
        os.makedirs(os.path.join(self.data, "Dictionaries", "genotypes_hg38_walk"))
        open(os.path.join(self.data, "Dictionaries", "genotypes_hg38_walk", "chr1"), "w").close()
        # outside the data folder: VCF folders with and without samples, a genome folder,
        # a full folder and an empty one, and two PAM files
        for folder, samples in (("cohort_geno", True), ("cohort_sites", False),
                                ("a/cohort_1000G", True), ("b/cohort_1000G", False)):
            _vcf(os.path.join(self.outside, folder, "x.chr1.vcf.gz"), samples)
        os.makedirs(os.path.join(self.outside, "genome_x"))
        shutil.copy(os.path.join(self.outside, "hg38", "chr1.fa"),
                    os.path.join(self.outside, "genome_x"))
        os.makedirs(os.path.join(self.outside, "store"))
        open(os.path.join(self.outside, "store", "chr1"), "w").close()
        os.makedirs(os.path.join(self.outside, "empty"))
        for name in ("pam", "x-20bp-y"):  # a PAM motif the guide's sequence does not hold
            with open(os.path.join(self.outside, f"{name}.txt"), "w") as fh:
                fh.write("N" * 20 + "CCC 3\n")
        self.before = self.outside_files()
        self.cohort = os.path.join(self.outside, "cohort")  # '<genome>_<dataset>' reaches
        self.through = "walk/../../../outside/"  # VCFs/hg38_<this>... is outside/...

    def assertNothingReadOutside(self, seen):
        outside = os.path.realpath(self.outside) + os.sep
        self.assertEqual([p for p in seen if p.startswith(outside)], [])
        self.assertNothingRun()

    def test_the_variant_datasets_offered(self):
        none = mp.change_variant_dataset_options(None)
        with _paths_read() as seen:
            for folder in ("a", "b", "c"):  # samples, sites-only, no folder
                with self.subTest(folder=folder):
                    self.assertEqual(mp.change_variant_dataset_options(
                        os.path.join(self.outside, folder, "cohort")), none)
        self.assertNothingReadOutside(seen)
        options, value = mp.change_variant_dataset_options("hg38")
        self.assertEqual(([o["value"] for o in options], value), (["ref", "walk"], "walk"))

    def test_the_analysis_modes(self):
        none = mp.update_search_mode_availability(None, None, "population-level")
        self.assertTrue(none[0][1]["disabled"])  # per-sample: no genotypes to resolve
        with _paths_read() as seen:
            for dataset in (self.through + "store", self.through + "empty"):
                with self.subTest(dataset=dataset):
                    self.assertEqual(mp.update_search_mode_availability(
                        dataset, "hg38", "population-level"), none)
            self.assertEqual(mp.update_search_mode_availability(
                "walk", self.cohort, "population-level"), none)
        self.assertNothingReadOutside(seen)
        options, _value = mp.update_search_mode_availability("walk", "hg38", "population-level")
        self.assertFalse(options[1]["disabled"])  # the genotype store of walk

    def test_the_bulge_options(self):
        none = mp.limit_bulges_to_index(None, PAM, "ref", 1, 1, 4, 6, {})
        ref = mp.limit_bulges_to_index("hg38", PAM, "ref", 1, 1, 4, 6, {})
        with _paths_read() as seen:
            for genome, pam, variant, want in (
                    (os.path.join(self.outside, "genome_x"), PAM, "ref", none),
                    (os.path.join(self.outside, "nothing_here"), PAM, "ref", none),
                    ("hg38", "../../outside/pam", "ref", none),
                    ("hg38", "../../outside/nothing_here", "ref", none),
                    ("hg38", PAM, self.through + "cohort_geno", ref),
                    ("hg38", PAM, self.through + "nothing_here", ref)):
                with self.subTest(genome=genome, pam=pam, variant=variant):
                    self.assertEqual(
                        mp.limit_bulges_to_index(genome, pam, variant, 1, 1, 4, 6, {}), want)
        self.assertNothingReadOutside(seen)
        # walk, which has no variant index: 0 bulges, as before
        opts, _, dna, rna, _, _ = mp.limit_bulges_to_index("hg38", PAM, "walk", 1, 1, 4, 6, {})
        self.assertEqual(([o["value"] for o in opts], dna, rna), ([0], 0, 0))

    def test_the_bulge_note(self):
        ref = mp.update_bulge_note("hg38", PAM, "ref", 1, 1)
        with _paths_read() as seen:
            for genome, pam, variant, want in (
                    (self.cohort, PAM, "sites", ""), (self.cohort, PAM, "geno", ""),
                    (self.cohort, PAM, "absent", ""),
                    (os.path.join(self.outside, "genome_x"), PAM, "ref", ""),
                    ("hg38", "../../outside/pam", "ref", ""),
                    ("hg38", "../../outside/nothing_here", "ref", ""),
                    ("hg38", PAM, self.through + "cohort_sites", ref),
                    ("hg38", PAM, self.through + "cohort_geno", ref),
                    ("hg38", PAM, self.through + "nothing_here", ref)):
                with self.subTest(genome=genome, pam=pam, variant=variant):
                    self.assertEqual(mp.update_bulge_note(genome, pam, variant, 1, 1), want)
        self.assertNothingReadOutside(seen)
        self.assertTrue(mp.update_bulge_note("hg38", PAM, "walk", 1, 1).startswith(
            "No walk variant index for this PAM is installed yet, so only a 0-bulge search "
            "is offered"))

    def submit(self, pam=PAM, variant="ref", guide_type="IP", guides="CTAACAGTTGCTTTTATCAC"):
        """check_input, as the Submit button runs it: (runs, the dialog's lines)."""
        out = mp.check_input(1, None, "hg38", pam, guide_type, guides, 4, 0, 0, variant, False)
        return out[0] == 1, [li.children for li in out[-1].children[1].children]

    def test_the_submit_check(self):
        with _paths_read() as seen:
            # a PAM file outside, and one not there: the same lines, but for its name
            answers = []
            for pam in ("../../outside/x-20bp-y", "../../outside/z-20bp-y"):
                ran, lines = self.submit(pam=pam, guide_type="GS", guides=">s\n" + SEQ)
                answers.append((ran, [line.replace(pam, "<pam>") for line in lines]))
            self.assertEqual(answers[0], answers[1])
            self.assertEqual(answers[0], (False, [
                "PAM: <pam> is not installed. Select an installed PAM, or add it in Settings "
                "(Define a nuclease / PAM)."]))
            for variant in (self.through + "cohort_sites", self.through + "cohort_geno",
                            self.through + "nothing_here", "walk+" + self.through + "empty"):
                with self.subTest(variant=variant):
                    self.assertEqual(self.submit(variant=variant), (False, [
                        f"Variants: {variant} is not installed for hg38. Select an installed "
                        "variant dataset, or add it in Settings (Add a VCF dataset)."]))
        self.assertNothingReadOutside(seen)
        self.assertEqual(self.submit(variant="walk"), (True, []))
        self.assertEqual(self.submit(guide_type="GS", guides=">s\n" + SEQ), (True, []))


@NEEDS_DASH
class TestAssemblySearchLaunch(_WebFolder):
    """The Personal assembly tab's Submit (submit_assembly_search_job), which a request can
    call with any guides, PAM and six fields, and the search's look-up of each guide of the
    guide file it writes (process_summaries.py, which runs a shell command holding the
    guide). The job pool is a stand-in, whose command runs as the pool runs it, and so is
    assembly-search, which records its arguments."""

    IND = "HG1"

    def setUp(self):
        super().setUp()
        self.pool = MagicMock()
        p = patch.object(mp, "pool_executor", self.pool)
        p.start()
        self.addCleanup(p.stop)
        # a personal assembly as Settings installs one (Assemblies/<individual>/)
        root = os.path.join(self.data, "Assemblies", self.IND)
        haplotypes = {}
        for hap in ("paternal", "maternal"):
            os.makedirs(os.path.join(root, hap, "genome"))
            shutil.copy(os.path.join(self.data, "Genomes", "hg38", "chr1.fa"),
                        os.path.join(root, hap, "genome"))
            for name in (f"{self.IND}_{hap}_vs_GRCh38.chain", f"{self.IND}_{hap}.chromAlias.txt"):
                with open(os.path.join(root, hap, name), "w") as fh:
                    fh.write("x\n")
            haplotypes[hap] = {"assembly_name": f"{self.IND}_{hap}", "genome_dir": f"{hap}/genome",
                               "chain": f"{hap}/{self.IND}_{hap}_vs_GRCh38.chain",
                               "chromalias": f"{hap}/{self.IND}_{hap}.chromAlias.txt"}
        with open(os.path.join(root, "metadata.json"), "w") as fh:
            json.dump({"format": "crisprme-personal-assembly", "format_version": 1,
                       "individual": self.IND, "haplotypes": haplotypes}, fh)
        with open(os.path.join(self.app, "crisprme.py"), "w") as fh:
            fh.write("import json, os, sys\nwith open(os.path.join(os.path.dirname(__file__), "
                     "'assembly_args.json'), 'w') as fh:\n    json.dump(sys.argv[1:], fh)\n")
        self.before = self.outside_files()

    def submit(self, guides, pam=PAM, fields=None, individual=IND):
        """Submit, as the form sends it with the individual picked (and the six fields the
        picker fills in for it, unless ``fields``): the guide file of the search it ran and
        None, or None and the lines of the warning it gave instead."""
        fields = fields or mp.populate_assembly_fields_from_individual(self.IND)
        result = mp.submit_assembly_search_job(
            1, *fields, pam, guides, 4, 1, 1, None, individual, 6, False)
        if not self.pool.submit.call_count:
            return None, [li.children for li in result[3].children]
        fn, *args = self.pool.submit.call_args[0]
        fn(*args)
        self.pool.reset_mock()
        with open(os.path.join(self.app, "assembly_args.json")) as fh:
            argv = json.load(fh)
        return argv[argv.index("--guide") + 1], None

    def look_up(self, guide_file):
        """process_summaries.py on guide_file, as the search runs it (here on a search that
        found nothing), in the data folder, where a file a command leaves is found."""
        best = os.path.join(self.data, "Results", "best.txt")
        open(best, "w").close()
        return subprocess.run(
            [sys.executable, os.path.join(HERE, "process_summaries.py"), best, guide_file, "_",
             "4", "1", "1", os.path.dirname(guide_file), "ref", "CFD"],
            cwd=self.data, stdin=subprocess.DEVNULL, capture_output=True, text=True)

    def test_a_guide_holding_a_command_or_a_path_runs_nothing(self):
        # the guides went into the guide file as sent, upper-cased: the look-up of one ran
        # the command in it, and wrote where a path in it led
        for guide in ("CTAACAG;>MARK_A;ATCA", "CTAACAG$(>MARK_B)TCA", "CTAACAG;>../MARK_C;A"):
            with self.subTest(guide=guide):
                guide_file, _ = self.submit(f"{guide}\n{GUIDE[:-3].lower()}")
                done = self.look_up(guide_file)
                self.assertNothingRun()
                self.assertEqual(done.returncode, 0, done.stderr)
                with open(guide_file) as fh:
                    self.assertEqual(fh.read(), f"{GUIDE}\n")
                job = os.path.dirname(guide_file)
                self.assertTrue(os.path.isfile(os.path.join(
                    job, f"{os.path.basename(job)}.summary_by_guide.{GUIDE}_CFD.txt")))
        # alone, it leaves no guide to search
        self.assertEqual(self.submit("CTAACAG;>MARK_A;ATCA"),
                         (None, ["at least one guide sequence"]))
        self.assertNothingRun()

    def test_a_guide_is_searched_in_the_letters_a_reference_search_keeps(self):
        # an N, a U and a space went into the guide file too
        guide_file, _ = self.submit("CTAACAGTTGCTTTTATCAN\nCTAACAGUUGCTTTTATCAC\n"
                                    "CTAACAG TGCTTTTATCAC\nctaacagttgcttttatcac\n"
                                    "RYSWKMBDHVACGTACGTAC")
        with open(guide_file) as fh:
            self.assertEqual(fh.read(), f"{GUIDE}\nRYSWKMBDHVACGTACGTACNNN\n")

    def test_only_an_installed_pam_and_the_files_of_the_individual_picked_are_searched(self):
        # any PAM file '../' led to, and any folders and files as the six fields the picker
        # fills in, were searched as sent: the search takes the PAM file's name, and that of
        # a genome folder of Genomes/, into its shell commands
        with open(os.path.join(self.outside, "NRG-20bp.txt"), "w") as fh:
            fh.write("N" * 20 + "NRG 3\n")  # (a name whose guide length the submit reads)
        self.before = self.outside_files()
        fields = list(mp.populate_assembly_fields_from_individual(self.IND))
        outside_genome = os.path.join("..", "outside", "hg38")
        outside_file = os.path.join(self.outside, f"{PAM}.txt")
        pick = ["select an individual on the Personal assembly tab"]
        for kwargs, said in (
                (dict(pam="../../outside/NRG-20bp"),
                 ["PAM (file unreadable)", "at least one guide sequence"]),
                (dict(fields=[outside_genome] + fields[1:]), pick),
                (dict(fields=fields[:3] + [outside_file] + fields[4:]), pick),
                (dict(individual="HG2"), pick)):
            with self.subTest(**kwargs):
                self.assertEqual(self.submit(GUIDE[:-3], **kwargs), (None, said))
                self.assertEqual(os.listdir(os.path.join(self.data, "Results")), [])
        guide_file, _ = self.submit(GUIDE[:-3])  # (as the form sends it)
        with open(guide_file) as fh:
            self.assertEqual(fh.read(), f"{GUIDE}\n")
        self.assertNothingRun()


@NEEDS_DASH
@unittest.skipUnless(shutil.which("sort-bed") and shutil.which("bgzip"),
                     "BEDOPS sort-bed or bgzip not on PATH")
class TestActiveAnnotation(_WebFolder):
    """The annotation a search applies, assembled from those enabled for its genome in
    Settings, in a file named after the genome."""

    def test_it_is_built_whatever_the_genome_is_named(self):
        for name, row in (("a.bed", "chr1\t10\t500\tpromoter\n"),
                          ("b.bed", "chr1\t2\t9\tenhancer\n")):
            with open(os.path.join(self.data, "Annotations", name), "w") as fh:
                fh.write(row)
        for genome in ("hg38", "x$(touch${IFS}MARK_a)", "y;touch${IFS}MARK_b;#"):
            with self.subTest(genome=genome):
                pu.write_enabled_annotations(genome, ["a.bed", "b.bed"])
                active, _gencode = pu.build_active_annotation(genome)
                self.assertEqual(active, f".active.{genome}.bed.gz")
                path = os.path.join(self.data, "Annotations", active)
                self.assertEqual(pd.read_csv(path, sep="\t", header=None)[3].tolist(),
                                 ["enhancer", "promoter"])  # both, sorted
        self.assertNothingRun()


@NEEDS_DASH
class TestRemoveResult(_WebFolder):
    """Remove result on a job's status page, which removes the job of its URL (and, for a
    personal-assembly search, the haplotype searches its .Params.txt names)."""

    def test_only_a_job_folder_is_removed(self):
        results = os.path.join(self.data, "Results")
        os.makedirs(os.path.join(results, "job1_ABC"))
        os.makedirs(os.path.join(self.outside, "victim"))
        self.before = self.outside_files()
        for search in ("?job=X;touch${IFS}MARK_a;#", "?job=../../outside/victim",
                       "?job=..", "?job=", "?job=nosuchjob"):
            with self.subTest(search=search), self.assertRaises(PreventUpdate):
                lp.remove_result(1, search)
        self.assertNothingRun()
        self.assertEqual(lp.remove_result(1, "?job=job1_ABC").children, "Results deleted")
        self.assertEqual(os.listdir(results), [])

    def test_a_link_leaving_results_is_not_a_job_folder(self):
        os.makedirs(os.path.join(self.outside, "victim"))
        os.symlink(os.path.join(self.outside, "victim"),
                   os.path.join(self.data, "Results", "linked"))
        self.assertIsNone(pu.results_job_dir("linked"))
        os.makedirs(os.path.join(self.data, "Results", "job_2.a-b"))
        self.assertEqual(pu.results_job_dir("job_2.a-b"),
                         os.path.join(self.data, "Results", "job_2.a-b"))

    def test_only_the_folders_of_an_assembly_search_are_removed(self):
        results = os.path.join(self.data, "Results")
        # asm1's .Params.txt as the Personal assembly tab writes it; asm2's names a
        # haplotype folder outside Results/
        for job, paternal in (("asm1_ABC", "asm1_ABC_paternal"),
                              ("asm2_ABC", "../../outside/victim")):
            os.makedirs(os.path.join(results, f"{job}_combined"))
            with open(os.path.join(results, f"{job}_combined", ".Params.txt"), "w") as fh:
                fh.write(f"1\tGenome_type\tassembly\n2\tPaternal_dir\t{paternal}\n"
                         f"3\tMaternal_dir\t{job}_maternal\n")
        for name in ("asm1_ABC_paternal", "asm1_ABC_maternal"):
            os.makedirs(os.path.join(results, name))
        os.makedirs(os.path.join(self.outside, "victim"))
        open(os.path.join(self.outside, "victim", "keep.txt"), "w").close()
        self.before = self.outside_files()
        for search in ("?job=../../outside/victim", "?job=X;touch${IFS}MARK_a;#", "?job=..",
                       "?job=", "?job=nosuchjob_combined"):
            with self.subTest(search=search), self.assertRaises(PreventUpdate):
                lp.remove_assembly_result(1, search)
        self.assertNothingRun()
        self.assertEqual(lp.remove_assembly_result(1, "?job=asm2_ABC_combined").children,
                         "Results deleted")
        self.assertNothingRun()  # its job folder only
        self.assertEqual(lp.remove_assembly_result(1, "?job=asm1_ABC_combined").children,
                         "Results deleted")
        self.assertEqual(os.listdir(results), [])  # the three folders of the search


@NEEDS_DASH
class TestStatusPage(_WebFolder):
    """A job's status page and the routes to a job's pages, which read the job of their
    URL: a search's and a personal-assembly search's job folders, and a folder beside the
    data folder that '../' reaches, holding the files they read."""

    OUT = "../../outside/logs"  # Results/<this> is the folder beside the data folder

    def setUp(self):
        super().setUp()
        results = os.path.join(self.data, "Results")
        for folder, files in (
            (os.path.join(results, "job1_ABC"),
             {"log_error.txt": "the job's error\n", "log_verbose.txt": "the job's step\n"}),
            (os.path.join(results, "asm1_ABC_combined"),
             {".Params.txt": "1\tGenome_type\tassembly\n",
              "log.txt": "Running paternal haplotype search\n"}),
            (os.path.join(self.outside, "logs"),
             {".Params.txt": "1\tGenome_type\tassembly\n", "log_error.txt": "OUTSIDE\n",
              "log_verbose.txt": "OUTSIDE\n",
              "log.txt": "OUTSIDE\nReconciliation complete. Wrote x\n"}),
        ):
            os.makedirs(folder)
            for name, text in files.items():
                with open(os.path.join(folder, name), "w") as fh:
                    fh.write(text)
        p = patch.object(index, "current_working_directory", self.data + os.sep)
        p.start()
        self.addCleanup(p.stop)
        self.before = self.outside_files()

    def test_the_logs_of_a_job_folder_only_are_shown(self):
        self.assertEqual(
            lp.show_full_log(1, "?job=job1_ABC", True),
            (False, "===== log_error.txt =====\nthe job's error\n\n"
                    "===== log_verbose.txt (tail) =====\nthe job's step"))
        self.assertEqual(lp.show_assembly_full_log(1, "?job=asm1_ABC_combined", True),
                         (False, "Running paternal haplotype search"))
        for search in (f"?job={self.OUT}", "?job=..", "?job=nosuchjob"):
            with self.subTest(search=search):
                self.assertEqual(lp.show_full_log(1, search, True),
                                 (False, "No log files found for this job yet."))
                self.assertEqual(lp.show_assembly_full_log(1, search, True),
                                 (False, "No log file found for this job yet."))
        self.assertNothingRun()

    def test_the_status_of_an_assembly_search_is_read_from_its_folder_only(self):
        status = lp.refresh_assembly_search(1, "?job=asm1_ABC_combined")
        self.assertEqual(status[1].children, "Searching paternal haplotype...")
        self.assertEqual(status[5], "")
        for search in (f"?job={self.OUT}", "?job=nosuchjob"):
            with self.subTest(search=search):
                status = lp.refresh_assembly_search(1, search)
                self.assertEqual(status[5].children, "The selected result does not exist")
        self.assertNothingRun()

    def test_the_routes_to_a_job_take_a_job_folder_only(self):
        def route(path, search):
            return index.change_page(f"http://127.0.0.1:8080{path}{search}", path, search, "")

        layout, link = route("/load", "?job=asm1_ABC_combined")
        self.assertTrue(link.endswith("/load?job=asm1_ABC_combined"))
        self.assertIn("div-assembly-status-report", repr(layout))  # its status page
        with patch.object(rp, "result_page_assembly", side_effect=lambda job: f"page of {job}"):
            self.assertEqual(route("/result", "?job=asm1_ABC_combined")[0],
                             "page of asm1_ABC_combined")
            for search in (f"?job={self.OUT}", "?job=..", "?job="):
                with self.subTest(search=search):
                    self.assertEqual(route("/load", search)[1], "/index")  # no job to show
                    page, _link = route("/result", search)
                    self.assertEqual(page.children.children,
                                     "The selected result does not exist")
        self.assertNothingRun()


SAMPLE_CMD = "S3;touch${IFS}MARK_s;"  # a VCF's sample name; Settings uploads one as it is


@NEEDS_DASH
@unittest.skipUnless(shutil.which("zip"), "zip not on PATH")
class TestResultsPages(_WebFolder):
    """The results pages of a finished search, and the same search beside the data
    folder, which '../' in a job id reaches."""

    JOB = "job_ABC"
    OUT = "../../outside/job_ABC"  # Results/<this> is the copy beside the data folder

    def setUp(self):
        super().setUp()
        self.job = self.finished(os.path.join(self.data, "Results", self.JOB))
        self.finished(os.path.join(self.outside, self.JOB))
        self.before = self.outside_files()

    def finished(self, d, integrated=f"{GUIDE}+NRG_hg38+__4+2_integrated_results.tsv"):
        """What the results pages read of a finished reference search (a minimal copy)."""
        job = os.path.basename(d)
        os.makedirs(os.path.join(d, "imgs"))
        files = {
            ".Params.txt": "Genome_selected\thg38\nGenome_ref\thg38\nGenome_idx\tNRG_2_hg38\n"
                           "Pam\tNRG\nMax_bulges\t2\nMismatches\t4\nDNA\t1\nRNA\t1\n"
                           "Max_total_edits\t6\nThreshold_mode\tadvanced\nNuclease\tSpCas9\n"
                           "Ref_comp\tFalse\n",
            ".guides.txt": f"{GUIDE}\n",
            "log.txt": "Job\tStart\tFri\nJob\tDone\tFri\n",
            ".dropdown.json": "CFD",
            f".{job}.acfd_CFD.txt": f"{GUIDE}\t0.95\tNA\tNA\n",
            f".{job}.general_target_count.{GUIDE}_CFD.txt":
                "0\t1\t2\t3\t4\n1\t3\t0\t0\t1\n3\t1\t2\t0\t2\n0\t0\t1\t0\t3\n0\t0\t0\t0\t0\n"
                "0\t0\t0\t0\t0\n0\t0\t0\t0\t0\n",
            f"{GUIDE}+NRG_hg38+__4+2_all_results_with_alternative_alignments.tsv": "",
        }
        for crit in ("CFD", "CRISPR-Bulge", "fewest"):
            files[f"{job}.summary_by_samples.{GUIDE}_{crit}.txt"] = (
                f"{GUIDE}\nSample\tSex\tPopulation\tSuper Population\tTargets in Sample\t"
                "Targets in Population\tTargets in Super Population\tPAM Creation\n"
                f"S1\tmale\tPOPA\tSUPX\t2\t2\t3\t0\nS2\tfemale\tPOPB\tSUPX\t2\t2\t3\t1\n"
                f"{SAMPLE_CMD}\tmale\tPOPA\tSUPX\t1\t1\t1\t0\n")
        for name, text in files.items():
            with open(os.path.join(d, name), "w") as fh:
                fh.write(text)
        rows = []
        for i, (mmb, samples) in enumerate(((0, "S1,S2"), (0, "S1"), (3, "S2"), (4, "NA"),
                                            (2, SAMPLE_CMD))):
            row = {"Spacer+PAM": GUIDE, "Chromosome": "chr1"}
            for crit in CRITERIA:
                row[f"Start_coordinate_({crit})"] = 1000 * (i + 1)
                row[f"Mismatches_({crit})"] = mmb
                row[f"Bulges_({crit})"] = 0
                row[f"Bulge_type_({crit})"] = "X"
                row[f"Mismatches+bulges_({crit})"] = mmb
                row[f"Variant_samples_({crit})"] = samples
                row[f"Variant_info_genome_({crit})"] = "NA"
            row["CFD_score_(highest_CFD)"] = 1.0 - 0.2 * i
            row["CRISPR_BULGE_score_(highest_CRISPR_BULGE)"] = 1.0 - 0.2 * i
            rows.append(row)
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(d, integrated), sep="\t", index=False)
        with sqlite3.connect(os.path.join(d, f".{job}.db")) as conn:
            df.to_sql("final_table", conn, index=False)
        return d

    def wait_for(self, path):
        """A file zip writes in the background."""
        for _ in range(100):
            if os.path.isfile(path):
                return True
            time.sleep(0.1)
        return False

    def test_a_region_link(self):
        # the region's file, as an earlier visit leaves it (the lookup is the shell's)
        open(os.path.join(self.job, f"{self.JOB}.chr1_1000.{GUIDE}.txt"), "w").close()
        page = rp.cluster_page(self.JOB, f"{GUIDE}-Pos-chr1-1000")
        self.assertEqual(page.children[0].children, "Selected Position: chr1 - 1000")
        for hash_term in (f"{GUIDE};touch${{IFS}}MARK_a;#-Pos-chr1-1000",
                          f"{GUIDE}-Pos-chr1$(touch${{IFS}}MARK_b)-1000",
                          f"{GUIDE}-Pos-chr1-1000;touch${{IFS}}MARK_c;#",
                          f"{GUIDE}-Pos-chr1'-1000", "OTHERGUIDE-Pos-chr1-1000"):
            with self.subTest(hash=hash_term):
                page = rp.cluster_page(self.JOB, hash_term)
                self.assertEqual(page.children.children, "The selected position does not exist")
        page = rp.cluster_page(self.OUT, f"{GUIDE}-Pos-chr1-1000")
        self.assertEqual(page.children.children, "The selected result does not exist")
        self.assertNothingRun()

    def test_a_sample_link(self):
        data, _cols = rp.update_table_sample(0, 10, [], None, f"?job={self.JOB}",
                                             f"#{GUIDE}-Sample-S1")
        self.assertEqual(len(data), 2)
        self.assertTrue(self.wait_for(
            os.path.join(self.job, f"{self.JOB}.S1.{GUIDE}.personal_targets.zip")))
        data, _cols = rp.update_table_sample(0, 10, [], None, f"?job={self.JOB}",
                                             f"#{GUIDE}-Sample-{SAMPLE_CMD}")
        self.assertEqual(len(data), 1)
        self.assertTrue(self.wait_for(
            os.path.join(self.job, f"{self.JOB}.{SAMPLE_CMD}.{GUIDE}.personal_targets.zip")))
        with self.assertRaises(PreventUpdate):
            rp.update_table_sample(0, 10, [], None, f"?job={self.OUT}", f"#{GUIDE}-Sample-S1")
        self.assertNothingRun()

    def test_a_show_targets_link(self):
        data, _cols = rp.update_table_subset(0, 10, [], None, [], f"?job={self.JOB}",
                                             f"#{GUIDE}newX00")
        self.assertEqual(len(data), 2)
        self.assertTrue(self.wait_for(
            os.path.join(self.job, f"{self.JOB}.X.0.0.{GUIDE}.targets.zip")))
        guide = "G;touch${IFS}MARK_d;#"
        data, _cols = rp.update_table_subset(0, 10, [], None, [], f"?job={self.JOB}",
                                             f"#{guide}newX00")
        self.assertEqual(data, [])
        self.assertTrue(self.wait_for(
            os.path.join(self.job, f"{self.JOB}.X.0.0.{guide}.targets.zip")))
        with self.assertRaises(PreventUpdate):
            rp.update_table_subset(0, 10, [], None, [], f"?job={self.OUT}", f"#{GUIDE}newX00")
        self.assertNothingRun()

    def test_the_result_summary(self):
        data, _cells = rp.update_table_general_profile(0, 10, [], "", "CFD",
                                                       f"?job={self.JOB}")
        self.assertEqual([row["Guide"] for row in data], [GUIDE])
        self.assertTrue(self.wait_for(os.path.join(
            self.job, f"{GUIDE}+NRG_hg38+__4+2_integrated_results.zip")))
        # a search of a genome named in Settings with a command: its results file holds
        # the name
        named = "named_ABC"
        integrated = f"{GUIDE}+NRG_x$(touch${{IFS}}MARK_e)+__4+2_integrated_results.tsv"
        d = self.finished(os.path.join(self.data, "Results", named), integrated=integrated)
        rp.update_table_general_profile(0, 10, [], "", "CFD", f"?job={named}")
        self.assertTrue(self.wait_for(os.path.join(d, integrated[:-3] + "zip")))
        with self.assertRaises(PreventUpdate):
            rp.update_table_general_profile(0, 10, [], "", "CFD", f"?job={self.OUT}")
        self.assertNothingRun()

    def test_the_graphical_reports(self):
        recorded = os.path.join(self.app, "PostProcess", "radar_args.json")
        rp.update_images_tabs("2", [{"row": 0}], "CFD", f"?job={self.JOB}",
                              [{"Guide": GUIDE}])
        with open(recorded) as fh:
            # the words of the command 2.8.0 ran
            self.assertEqual(json.load(fh), [
                GUIDE, f"{self.job}/.guide_dict_{GUIDE}_CFD.json",
                f"{self.job}/.motif_dict_{GUIDE}_CFD.json", "2", "0", "TOTAL_CFD",
                f"{self.job}/imgs/"])
        for mm, guide in (("2", "G;touch${IFS}MARK_f;#"), ("2;touch${IFS}MARK_g;#", GUIDE)):
            with self.subTest(mm=mm, guide=guide):
                rp.update_images_tabs(mm, [{"row": 0}], "CFD", f"?job={self.JOB}",
                                      [{"Guide": guide}])
                with open(recorded) as fh:
                    self.assertEqual(json.load(fh)[:4:3], [guide, mm])
        with self.assertRaises(PreventUpdate):
            rp.update_images_tabs("2", [{"row": 0}], "CFD", f"?job={self.OUT}",
                                  [{"Guide": GUIDE}])
        self.assertNothingRun()

    def test_the_personal_risk_card(self):
        for sample, private_sites in (("S1", 1), (SAMPLE_CMD, 1)):
            with self.subTest(sample=sample):
                rp.generate_sample_card(1, "CFD", sample, [{"row": 0}], [{"Guide": GUIDE}],
                                        f"?job={self.JOB}")
                private = os.path.join(self.job, f"{self.JOB}.{sample}.{GUIDE}.private_targets")
                self.assertEqual(len(pd.read_csv(private + ".tsv", sep="\t")), private_sites)
                self.assertTrue(os.path.isfile(private + ".zip"))
                self.assertTrue(os.path.isfile(os.path.join(
                    self.job, "imgs",
                    f"CRISPRme_CFD_top_1000_log_for_main_text_{GUIDE}.{sample}.personal.png")))
        with self.assertRaises(PreventUpdate):
            rp.generate_sample_card(1, "CFD", "S1", [{"row": 0}], [{"Guide": GUIDE}],
                                    f"?job={self.OUT}")
        self.assertNothingRun()

    def test_the_table_filter_is_kept_in_the_job_folder_only(self):
        # the criterion of the results pages' dropdown, which their tables read back
        self.assertIsNone(rp.sendto_write_json("fewest", f"?job={self.JOB}"))
        self.assertEqual(pu.read_json(self.JOB), "fewest")
        with self.assertRaises(PreventUpdate):
            rp.sendto_write_json("fewest", f"?job={self.OUT}")
        for call in (lambda: pu.write_json("fewest", self.OUT), lambda: pu.read_json(self.OUT)):
            with self.assertRaises(FileNotFoundError):
                call()
        with open(os.path.join(self.outside, self.JOB, ".dropdown.json")) as fh:
            self.assertEqual(fh.read(), "CFD")  # the copy's, as it was
        self.assertNothingRun()

    def test_the_tables_lists_and_links_read_the_job_folder_only(self):
        """Each of the other callbacks that reads the job of the URL: the job's files, and
        no file of its copy beside the data folder."""
        var, asm = "var_ABC", "asm_ABC_combined"  # a variant search's, an assembly one's
        targets = ["#Bulge type", "crRNA", "DNA", "Chromosome", "Position",
                   "Cluster Position", "Direction", "Mismatches", "Bulge Size", "Total"]
        targets += [f"col{i}" for i in range(len(targets), 38)]
        for top in (os.path.join(self.data, "Results"), self.outside):
            files = {
                ".sampleID.txt": "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                                 "S1\tPOPA\tSUPX\tmale\nS2\tPOPB\tSUPX\tfemale\n",
                f"{self.JOB}.chr1_1000.{GUIDE}.txt": "Chromosome\tPosition\tSamples\n"
                                                     "chr1\t1000\tn\n",
                f"{self.JOB}.targets.txt": "\t".join(targets) + "\n" + "\t".join(
                    ["X", GUIDE, "CTAACAGTTGCTTTTATCACAGG", "chr1", "1000", "1000", "+"]
                    + ["0"] * (len(targets) - 7)) + "\n",
                f"{self.JOB}_report.zip": "",
                f"{self.JOB}.S1.{GUIDE}.personal_targets.zip": "",
                f"{self.JOB}.X.0.0.{GUIDE}.targets.zip": "",
            }
            for name, text in files.items():
                with open(os.path.join(top, self.JOB, name), "w") as fh:
                    fh.write(text)
            d = self.finished(os.path.join(top, var))
            with open(os.path.join(d, ".Params.txt")) as fh:
                params = fh.read().replace("Ref_comp\tFalse", "Ref_comp\tTrue")
            with open(os.path.join(d, ".Params.txt"), "w") as fh:
                fh.write(params)
            with open(os.path.join(d, f"{var}.chr1_1000.{GUIDE}.scomposition.txt"), "w") as fh:
                fh.write("Chromosome\tPosition\tSamples\nchr1\t1000\tS1\n")
            os.makedirs(os.path.join(top, asm))
            with open(os.path.join(top, asm, ".Params.txt"), "w") as fh:
                fh.write("1\tGenome_type\tassembly\n")
            pd.DataFrame([{"origin": "both", "hg38_chr": "chr1", "hg38_start": 1000,
                           "Spacer+PAM_paternal": GUIDE, "Spacer+PAM_maternal": GUIDE,
                           "CFD_score_(fewest_mm+b)_paternal": 0.9}]).to_csv(
                os.path.join(top, asm, f"{asm}_hg38.tsv"), sep="\t", index=False)
        self.before = self.outside_files()
        # the tables' loaders cache what they read in Cache/ of the folder the app runs in
        os.makedirs(os.path.join(self.data, "Cache"), exist_ok=True)
        sel, guides, region = [{"row": 0}], [{"Guide": GUIDE}], f"#{GUIDE}-Pos-chr1-1000"
        for job, read in (
            (asm, lambda s: rp.update_assembly_site_set_controls("mappable", s)),
            (asm, lambda s: rp.update_assembly_region_chrom_options("hg38", "mappable", s)),
            (asm, lambda s: rp.update_assembly_results_table("mappable", "hg38_start", "asc",
                                                             None, s)),
            (self.JOB, lambda s: rp.download_link_sample(
                1, f"{self.JOB}.summary_by_samples.{GUIDE}_CFD", s)),
            (self.JOB, lambda s: rp.download_report(1, f"{self.JOB}_report.zip", s)),
            (self.JOB, lambda s: rp.download_link_sumbysample(
                1, f"{self.JOB}.S1.{GUIDE}.personal_targets", s)),
            (self.JOB, lambda s: rp.downloadLinkGuide(1, f"{self.JOB}.X.0.0.{GUIDE}.targets", s)),
            (var, lambda s: rp.update_iupac_decomposition_table_cluster(0, 10, "", s, region)),
            (self.JOB, lambda s: rp.update_table_cluster(0, 10, [], "", "", s, region)),
            (self.JOB, lambda s: rp.load_distribution_populations(sel, guides, s)),
            (self.JOB, lambda s: rp.filter_position_table("chr1,0,10000", 1, "CFD", s, sel,
                                                          guides, "1/1")),
            (self.JOB, lambda s: rp.filter_sample_table(None, None, "None,None,None", "CFD", 1,
                                                        s, sel, guides, "1/1")),
            (self.JOB, lambda s: rp.update_sample_drop("POPA", s)),
            (self.JOB, lambda s: rp.update_population_drop("SUPX", s)),
            (self.JOB, lambda s: rp.update_content_tab("tab-summary-by-sample", sel, "CFD",
                                                       guides, s, "both")),
            (self.JOB, lambda s: rp.update_table(0, 10, [], "", s, f"#{GUIDE}")),
            (self.JOB, lambda s: rp.update_output(1, 0, "CFD", 10, sel, "CFD", "Mismatches",
                                                  guides, None, None, "ASC", None, s, False)),
        ):
            with self.subTest(callback=read.__code__.co_names[1]):  # rp.<the callback>(...)
                self.assertTrue(read(f"?job={job}"))  # the job's own
                with self.assertRaises(PreventUpdate):
                    read(f"?job=../../outside/{job}")  # its copy beside the data folder
        self.assertNothingRun()


if __name__ == "__main__":
    unittest.main()
