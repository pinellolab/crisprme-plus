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
  guide of that file to a shell command (process_summaries.py);
* the active annotation of a genome (pages_utils.build_active_annotation), whose file is
  named after the genome, a name Settings takes as typed;
* Remove result on a job's status page (load_page.remove_result): the job id of the URL;
* the results pages: the guide, sample, chromosome and position of a results link, the
  guide and mismatches the tables send back, a sample's name (a VCF's), and the results
  file's name, which holds the genome's.

Each case sends a value holding a command that would create a marker file, and a '../'
path, and checks that no marker and no file outside the data folder appear, and that the
normal value still works. The extraction's checks run everywhere (those that read a
sequence need bedtools); the pages' need the web stack (Dash), and run in the image.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_shell_safety.py' -v
"""

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
    def test_a_genome_named_with_a_command_in_settings_runs_nothing(self):
        # Settings takes a genome's name as typed; its active annotation is named after it
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


@NEEDS_DASH
class TestAssemblySearchLaunch(_WebFolder):
    """The Personal assembly tab's Submit (submit_assembly_search_job), which a request can
    call with any guides, and the search's look-up of each guide of the guide file it writes
    (process_summaries.py, which runs a shell command holding the guide). The job pool is a
    stand-in, whose command runs as the pool runs it, and so is assembly-search, which
    records its arguments."""

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

    def submit(self, guides):
        """Submit, as the form sends it with the individual picked: the guide file of the
        search it ran and None, or None and the lines of the warning it gave instead."""
        result = mp.submit_assembly_search_job(
            1, *mp.populate_assembly_fields_from_individual(self.IND), PAM, guides, 4, 1, 1,
            None, self.IND, 6, False)
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
    """Remove result on a job's status page, which removes the job of its URL."""

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


if __name__ == "__main__":
    unittest.main()
