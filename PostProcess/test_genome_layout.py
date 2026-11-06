"""Tests for the one-FASTA-per-sequence genome layout (``genome_layout.py``) and
the places that enforce it (``utils.download_reference_genome``,
``validate_inputs.check_genome_fasta``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_genome_layout.py' -v

Background: CRISPRme keys its chromosome list by FASTA file name, so a genome
delivered as ONE multi-sequence file (UCSC's ``susScr11.fa.gz``) indexes and
searches fine and then silently returns an empty result set.
"""

import errno
import gzip
import io
import os
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import unittest
import warnings
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import genome_layout as gl  # noqa: E402
import safe_tar  # noqa: E402
import utils  # noqa: E402
import validate_inputs as vi  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "genome_layout.py")

MULTI = (
    ">chr1 first sequence\nACGTACGT\nACGT\n"
    ">chr2\nGGGGCCCC\n"
    ">chrUn_x_random description here\nTTTT\nAA\n"
)


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def _read(path):
    with open(path) as f:
        return f.read()


class TestRecordIds(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_ids_are_first_token_of_each_header(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, MULTI)
        self.assertEqual(gl.fasta_record_ids(p), ["chr1", "chr2", "chrUn_x_random"])

    def test_limit_stops_early(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, MULTI)
        self.assertEqual(gl.fasta_record_ids(p, limit=1), ["chr1"])

    def test_single_record_and_empty_and_headerless(self):
        one = os.path.join(self.tmp, "chr1.fa")
        _write(one, ">chr1\nACGT\n")
        self.assertEqual(gl.fasta_record_ids(one), ["chr1"])
        empty = os.path.join(self.tmp, "e.fa")
        _write(empty, "")
        self.assertEqual(gl.fasta_record_ids(empty), [])
        junk = os.path.join(self.tmp, "j.fa")
        _write(junk, "ACGT\nACGT\n")
        self.assertEqual(gl.fasta_record_ids(junk), [])

    def test_last_header_without_trailing_newline(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, ">a\nAC\n>b")
        self.assertEqual(gl.fasta_record_ids(p), ["a", "b"])

    def test_empty_sequence_between_headers(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, ">a\n>b\nAC\n")
        self.assertEqual(gl.fasta_record_ids(p), ["a", "b"])


class TestFindProblems(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_normal_per_chromosome_genome_is_clean(self):
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nACGT\n")
        _write(os.path.join(self.tmp, "chr2.fa"), ">chr2\nACGT\n")
        self.assertEqual(gl.find_multi_record_fastas(self.tmp), {})
        self.assertEqual(gl.find_name_mismatched_fastas(self.tmp), {})

    def test_multi_sequence_file_named_after_assembly_is_found(self):
        _write(os.path.join(self.tmp, "susScr11.fa"), MULTI)
        found = gl.find_multi_record_fastas(self.tmp)
        self.assertEqual(list(found), ["susScr11.fa"])
        self.assertEqual(len(found["susScr11.fa"]), 3)

    def test_fast_mode_skips_files_named_after_their_first_sequence(self):
        # chr1.fa holding chr1 + extra sequences: only the deep scan catches it
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nAC\n>extra\nGG\n")
        self.assertEqual(gl.find_multi_record_fastas(self.tmp), {})
        self.assertEqual(list(gl.find_multi_record_fastas(self.tmp, deep=True)), ["chr1.fa"])

    def test_single_sequence_with_wrong_name_is_a_name_mismatch_only(self):
        _write(os.path.join(self.tmp, "NC_000001.fa"), ">chr1\nACGT\n")
        self.assertEqual(gl.find_multi_record_fastas(self.tmp), {})
        self.assertEqual(gl.find_name_mismatched_fastas(self.tmp), {"NC_000001.fa": "chr1"})

    def test_non_fasta_files_ignored(self):
        _write(os.path.join(self.tmp, "notes.txt"), ">a\n>b\n")
        self.assertEqual(gl.find_multi_record_fastas(self.tmp, deep=True), {})


class TestSplit(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_split_is_lossless_and_named_after_records(self):
        p = os.path.join(self.tmp, "susScr11.fa")
        _write(p, MULTI)
        pieces = gl.split_multi_fasta(p)
        self.assertFalse(os.path.exists(p))  # original removed
        self.assertEqual(
            sorted(os.path.basename(x) for x in pieces),
            ["chr1.fa", "chr2.fa", "chrUn_x_random.fa"],
        )
        joined = "".join(_read(os.path.join(self.tmp, n)) for n in ("chr1.fa", "chr2.fa", "chrUn_x_random.fa"))
        # each header as >id, the sequence byte for byte
        self.assertEqual(joined, ">chr1\nACGTACGT\nACGT\n>chr2\nGGGGCCCC\n>chrUn_x_random\nTTTT\nAA\n")
        self.assertEqual(gl.find_multi_record_fastas(self.tmp, deep=True), {})
        self.assertEqual(gl.find_name_mismatched_fastas(self.tmp), {})  # names now match
        self.assertEqual([f for f in os.listdir(self.tmp) if f.endswith(".part")], [])

    def test_duplicate_ids_refused_and_nothing_written(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, ">a\nAC\n>a\nGG\n")
        with self.assertRaises(ValueError):
            gl.split_multi_fasta(p)
        self.assertEqual(os.listdir(self.tmp), ["g.fa"])

    def test_unusable_id_refused(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, ">ok\nAC\n>bad/id\nGG\n")
        with self.assertRaises(ValueError):
            gl.split_multi_fasta(p)
        self.assertEqual(os.listdir(self.tmp), ["g.fa"])

    def test_refuses_to_overwrite_an_existing_file(self):
        p = os.path.join(self.tmp, "g.fa")
        _write(p, ">a\nAC\n>b\nGG\n")
        _write(os.path.join(self.tmp, "b.fa"), ">b\nKEEP\n")
        with self.assertRaises(FileExistsError):
            gl.split_multi_fasta(p)
        self.assertEqual(_read(os.path.join(self.tmp, "b.fa")), ">b\nKEEP\n")
        self.assertTrue(os.path.exists(p))

    def test_first_record_may_reuse_the_original_file_name(self):
        p = os.path.join(self.tmp, "chr1.fa")
        _write(p, ">chr1\nAC\n>chr2\nGG\n")
        gl.split_multi_fasta(p)
        self.assertEqual(_read(p), ">chr1\nAC\n")
        self.assertEqual(_read(os.path.join(self.tmp, "chr2.fa")), ">chr2\nGG\n")


class TestNameAfterSequence(unittest.TestCase):
    """A FASTA of one sequence named otherwise (ecoli.fa holding >NC_000913.3, or a
    one-chromosome subset hg38_chr22.fa): CRISPRme takes each file for the sequence it
    is named after, so every search of it found its off-targets and then dropped them
    all. split_multi_record_fastas, which every genome download and Settings upload
    runs, names it after its sequence (round 3 review, F2)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_renamed_beside_a_split_and_a_named_file(self):
        _write(os.path.join(self.tmp, "ecoli.fa"), ">NC_000913.3 E. coli K-12\nACGT\nAC\n")
        _write(os.path.join(self.tmp, "asm.fa"), ">chr1\nAC\n>chr2\nGG\n")
        _write(os.path.join(self.tmp, "chrM.fa"), ">chrM\nTT\n")
        done = gl.split_multi_record_fastas(self.tmp)
        self.assertEqual(sorted(done), ["asm.fa", "ecoli.fa"])
        self.assertEqual(done["ecoli.fa"], [os.path.join(self.tmp, "NC_000913.3.fa")])
        self.assertEqual(sorted(os.listdir(self.tmp)),
                         ["NC_000913.3.fa", "chr1.fa", "chr2.fa", "chrM.fa"])
        self.assertEqual(_read(os.path.join(self.tmp, "NC_000913.3.fa")), ">NC_000913.3\nACGT\nAC\n")
        self.assertEqual(gl.find_name_mismatched_fastas(self.tmp), {})

    def test_refuses_a_name_another_file_has(self):
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nKEEP\n")
        _write(os.path.join(self.tmp, "extra.fa"), ">chr1\nGG\n")
        with self.assertRaises(FileExistsError) as cm:
            gl.split_multi_record_fastas(self.tmp)
        self.assertIn("extra.fa: cannot name it after its sequence", str(cm.exception))
        self.assertEqual(_read(os.path.join(self.tmp, "chr1.fa")), ">chr1\nKEEP\n")
        self.assertEqual(_read(os.path.join(self.tmp, "extra.fa")), ">chr1\nGG\n")

    def test_a_name_that_differs_in_case_only_is_the_file_itself(self):
        # on a case-insensitive file system (a macOS data folder), chr1.fa is CHR1.fa;
        # on another, a second link to the file stands for that
        _write(os.path.join(self.tmp, "CHR1.fa"), ">chr1\nAC\n")
        if not os.path.exists(os.path.join(self.tmp, "chr1.fa")):
            os.link(os.path.join(self.tmp, "CHR1.fa"), os.path.join(self.tmp, "chr1.fa"))
        gl.split_multi_record_fastas(self.tmp)
        self.assertEqual(os.listdir(self.tmp), ["chr1.fa"])
        self.assertEqual(_read(os.path.join(self.tmp, "chr1.fa")), ">chr1\nAC\n")


class TestHeaderIsExactlyTheId(unittest.TestCase):
    """CRISPRitz names a chromosome after the whole header line, CR and all, and the
    pipeline splits on blanks: an NCBI, GENCODE or Ensembl header, or CRLF line ends,
    indexed as done, and every search of it lost results, often all of them with no
    error (round 5 review, F1)."""

    # file name, what the file holds, its id
    STYLES = {
        "NCBI": ("NC_000913.3.fa", ">NC_000913.3 Escherichia coli str. K-12, complete genome\nACGT\nAC\n", "NC_000913.3"),
        "GENCODE": ("chr1.fa", ">chr1 1\nACGT\nAC\n", "chr1"),
        "Ensembl": ("1.fa", ">1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF\nACGT\nAC\n", "1"),
        "GATK": ("chr1.fa", ">chr1  AC:CM000663.2  gi:568336023  LN:248956422\nACGT\nAC\n", "chr1"),
        "tab": ("chr1.fa", ">chr1\tdescription\nACGT\nAC\n", "chr1"),
        "a blank after >": ("chr1.fa", "> chr1\nACGT\nAC\n", "chr1"),
        "a blank at the end": ("chr1.fa", ">chr1 \nACGT\nAC\n", "chr1"),
        "CRLF": ("chr1.fa", ">chr1\r\nACGT\r\nAC\r\n", "chr1"),
        "CRLF, a description": ("chr1.fa", ">chr1 1\r\nACGT\r\nAC\r\n", "chr1"),
        "LF after the header, CRLF after it": ("chr1.fa", ">chr1\nACGT\r\nAC\r\n", "chr1"),
    }

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _fresh(self):
        shutil.rmtree(self.tmp, ignore_errors=True)
        os.makedirs(self.tmp)

    def _writeb(self, name, text):
        with open(os.path.join(self.tmp, name), "wb") as f:
            f.write(text.encode())

    def _readb(self, name):
        with open(os.path.join(self.tmp, name), "rb") as f:
            return f.read().decode()

    def test_each_style_is_found_refused_and_written_as_the_id(self):
        for style, (name, text, rid) in self.STYLES.items():
            with self.subTest(style):
                self._fresh()
                self._writeb(name, text)
                self.assertEqual(list(gl.find_unclean_fastas(self.tmp)), [name])
                self.assertEqual(list(gl.find_unclean_fastas(self.tmp, deep=True)), [name])
                problems = gl.describe_header_problems(self.tmp)
                self.assertEqual(len(problems), 1)
                self.assertIn(f"genome_layout.py --split {self.tmp}", problems[0])
                self.assertIn("remove each index built from the folder", problems[0])
                issues, _ = vi.check_genome_fasta(self.tmp)
                self.assertEqual([i.severity for i in issues], [vi.ERROR])
                self.assertEqual(issues[0].message, problems[0])
                done = gl.split_multi_record_fastas(self.tmp)
                self.assertEqual(done, {name: [os.path.join(self.tmp, name)]})
                self.assertEqual(os.listdir(self.tmp), [name])
                self.assertEqual(self._readb(name), f">{rid}\nACGT\nAC\n")
                self.assertEqual(gl.find_unclean_fastas(self.tmp, deep=True), {})
                self.assertEqual(vi.check_genome_fasta(self.tmp)[0], [])

    def test_a_description_in_another_encoding(self):
        with open(os.path.join(self.tmp, "chr1.fa"), "wb") as f:
            f.write(b">chr1 caf\xe9\nACGT\n")
        issues, _ = vi.check_genome_fasta(self.tmp)  # (it read the file as UTF-8, and crashed)
        self.assertEqual([i.severity for i in issues], [vi.ERROR])
        self.assertIn("its header line reads '>chr1 caf", issues[0].message)
        gl.split_multi_record_fastas(self.tmp)
        self.assertEqual(self._readb("chr1.fa"), ">chr1\nACGT\n")

    def test_a_clean_genome_is_not_touched(self):
        for name in ("chr1.fa", "chr2.fa"):
            self._writeb(name, f">{name[:-3]}\nACGT\nAC")  # (no last line end, as it was)
        before = {n: os.stat(os.path.join(self.tmp, n)) for n in os.listdir(self.tmp)}
        self.assertEqual(gl.find_unclean_fastas(self.tmp, deep=True), {})
        self.assertEqual(gl.describe_header_problems(self.tmp), [])
        self.assertEqual(gl.split_multi_record_fastas(self.tmp), {})
        for name, stat in before.items():
            after = os.stat(os.path.join(self.tmp, name))
            self.assertEqual((after.st_ino, after.st_mtime_ns), (stat.st_ino, stat.st_mtime_ns))

    def test_the_named_wrongly_the_several_and_the_rest_of_the_sequence(self):
        self._writeb("ecoli.fa", ">NC_000913.3 E. coli\r\nACGT\r\nAC\r\n")
        self._writeb("asm.fa", ">chr1 first\r\nACGT\r\n>chr2\tsecond\r\nGG\r\nCC")
        gl.split_multi_record_fastas(self.tmp)
        self.assertEqual(sorted(os.listdir(self.tmp)), ["NC_000913.3.fa", "chr1.fa", "chr2.fa"])
        self.assertEqual(self._readb("NC_000913.3.fa"), ">NC_000913.3\nACGT\nAC\n")
        self.assertEqual(self._readb("chr1.fa"), ">chr1\nACGT\n")
        self.assertEqual(self._readb("chr2.fa"), ">chr2\nGG\nCC")  # no line end, as it had none

    def test_a_cr_further_on_is_for_the_deep_check_and_the_fix_to_find(self):
        self._writeb("chr1.fa", ">chr1\nACGT\nAC\r\nGG\n")
        self.assertEqual(gl.find_unclean_fastas(self.tmp), {})  # two lines a file, per search
        self.assertEqual(gl.find_unclean_fastas(self.tmp, deep=True), {"chr1.fa": "its lines end in CR"})
        gl.split_multi_record_fastas(self.tmp)
        self.assertEqual(self._readb("chr1.fa"), ">chr1\nACGT\nAC\nGG\n")

    def test_a_fasta_extension_goes_to_fa_and_only_fa_files_are_a_searchs_concern(self):
        self._writeb("chr1.fa", ">chr1\nACGT\n")
        self._writeb("notes.fasta", ">notes a description\nACGT\n")
        self.assertEqual(list(gl.find_unclean_fastas(self.tmp)), ["notes.fasta"])
        self.assertEqual(gl.describe_header_problems(self.tmp), [])  # index-genome reads .fa only
        self._fresh()
        self._writeb("chr1.fasta", ">chr1 description\nACGT\n")
        gl.split_multi_record_fastas(self.tmp)
        self.assertEqual(os.listdir(self.tmp), ["chr1.fa"])
        self.assertEqual(self._readb("chr1.fa"), ">chr1\nACGT\n")

    def test_lines_that_end_in_cr_alone_are_refused_with_nothing_lost(self):
        # read as one line, whose id alone would have been kept
        text = ">chr1\rACGT\rAC\r"
        self._writeb("chr1.fa", text)
        self.assertEqual(gl.find_unclean_fastas(self.tmp), {"chr1.fa": "its lines end in CR"})
        with self.assertRaises(ValueError) as cm:
            gl.split_multi_record_fastas(self.tmp)
        self.assertIn("its lines end in CR alone", str(cm.exception))
        self.assertEqual(os.listdir(self.tmp), ["chr1.fa"])
        self.assertEqual(self._readb("chr1.fa"), text)

    def test_an_id_the_pipeline_cannot_take_is_refused(self):
        for char in "|;'\"$()&<>`#":
            with self.subTest(char):
                self._fresh()
                rid = f"chr1{char}x"
                # the id is fine as a file name: as it is, the search stops, or drops all
                self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
                problems = gl.describe_header_problems(self.tmp)
                self.assertEqual(len(problems), 1)
                self.assertIn(f"sequence id {rid!r} holds {char!r}", problems[0])
                self.assertIn("Rename the sequence in the FASTA header", problems[0])
                self.assertEqual([i.severity for i in vi.check_genome_fasta(self.tmp)[0]], [vi.ERROR])
                self._fresh()
                self._writeb("asm.fa", f">chr2\nGG\n>{rid} a description\nACGT\n")
                with self.assertRaises(ValueError) as cm:
                    gl.split_multi_record_fastas(self.tmp)
                self.assertIn(f"sequence id {rid!r} holds {char!r}", str(cm.exception))
                self.assertEqual(os.listdir(self.tmp), ["asm.fa"])

    def test_ids_the_pipeline_takes_are_not_refused(self):
        # measured so through the real CLI, with a build and a search, alone in the folder
        # or beside one other; one that holds '*', '?' or '[' as long as it matches no
        # other name there, nor a file of the PostProcess folder, where the post-analysis
        # runs ('a*' matches adjust_cols.py, and a search of it failed: this check does not
        # see it); and '~' at the start as long as no user of the machine has the name
        # that follows ('~root' made a search fail, '~nouser' worked)
        for rid in ("NC_000913.3", "chr1_KI270706v1_random", "chr1:1-50000", "1", "HLA-A*01:01",
                    "seq=1,2", "X@y%z^~!+", "c[12]", "seq?1", "x-n", "a~b", "_x", ":x", "+x",
                    "~nouser"):
            with self.subTest(rid):
                self._fresh()
                self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
                self.assertEqual(gl.describe_header_problems(self.tmp), [])
                self.assertEqual(vi.check_genome_fasta(self.tmp)[0], [])

    def test_an_id_that_starts_with_a_dash_or_a_dot_is_refused(self):
        # measured through the real CLI: a search of -x (or -n) or .x beside b ended with
        # exit status 0 and only b's off-targets, none of the other's, its on-target among
        # them (the post-analysis grep took -x for an option; the search script lists the
        # chromosomes as the folder's *.fa, which leaves out .x.fa)
        for rid in ("-x", "-n", "--help", ".x"):
            with self.subTest(rid):
                self._fresh()
                self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
                self._writeb("b.fa", ">b\nACGT\n")
                problems = gl.describe_header_problems(self.tmp)
                self.assertEqual(len(problems), 1, problems)
                other = "." if rid[0] == "-" else "-"
                self.assertTrue(problems[0].startswith(
                    f"{rid}.fa: sequence id {rid!r} starts with {rid[0]!r}, which the pipeline "
                    f"cannot take (nor one that starts with {other!r}): a search of it finds "
                    "its off-targets and then drops them all, with no error. Rename the "
                    "sequence in the FASTA header and the file after it, then remove each "
                    "index built from the folder"), problems[0])
                issues, _ = vi.check_genome_fasta(self.tmp)
                self.assertEqual([(i.severity, i.message) for i in issues], [(vi.ERROR, problems[0])])
                # renamed, the folder is as a search needs it
                os.remove(os.path.join(self.tmp, f"{rid}.fa"))
                self._writeb("x.fa", ">x\nACGT\n")
                self.assertEqual(gl.describe_header_problems(self.tmp), [])
                # a split (an upload, a download) refuses it before anything is written
                self._fresh()
                self._writeb("asm.fa", f">chr2\nGG\n>{rid} a description\nACGT\n")
                with self.assertRaises(ValueError) as cm:
                    gl.split_multi_record_fastas(self.tmp)
                self.assertIn(f"sequence id {rid!r} starts with '-', which the pipeline cannot "
                              "take" if rid[0] == "-" else
                              f"sequence id {rid!r} cannot be used as a file name", str(cm.exception))
                self.assertEqual(os.listdir(self.tmp), ["asm.fa"])
        # one message for each kind, whatever the number of files
        self._fresh()
        for rid in ("-x", "-y", "gi|1", "gi|2"):
            self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
        problems = gl.describe_header_problems(self.tmp)
        self.assertEqual([p.split(":", 1)[0] for p in problems], ["gi|1.fa", "-x.fa"])

    def test_an_id_that_matches_another_name_as_a_shell_pattern_is_refused(self):
        # the search script and its post-analysis take each name into shell commands
        # unquoted: measured through the real CLI, each of these searches failed, or with
        # --thread 1 ended with exit status 0 and none of the first one's off-targets (all
        # but the second pair; the third is the case of HLA-C*02:11 and HLA-C*01:02:11, and
        # of HLA-DQA1*03:01:01 and HLA-DQA1*01:03:01:01, in the GATK GRCh38 FASTA)
        for ids in (("c[12]", "c1"), ("seq?1", "seqA1"), ("HLA-X*03:01:01", "HLA-X*01:03:01:01"),
                    ("seq?1", "seq01"), ("HLA-A*02:10", "HLA-A*24:02:10"), ("x[^a]", "xb")):
            with self.subTest(ids):
                self._fresh()
                for rid in ids:
                    self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
                self.assertEqual(gl.find_pattern_clashes(self.tmp), {f"{ids[0]}.fa": [f"{ids[1]}.fa"]})
                problems = gl.describe_header_problems(self.tmp)
                self.assertEqual(len(problems), 1)
                self.assertTrue(problems[0].startswith(
                    f"{ids[0]}.fa: its sequence id {ids[0]!r}, read as a shell pattern, matches "
                    f"{ids[1]}.fa. "), problems[0])
                self.assertIn("a search of this folder fails, or leaves out the off-targets of "
                              "such a sequence with no error", problems[0])
                self.assertIn("then remove each index built from the folder", problems[0])
                issues, _ = vi.check_genome_fasta(self.tmp)
                self.assertEqual([(i.severity, i.message) for i in issues], [(vi.ERROR, problems[0])])
                # renamed with none of those characters, the folder is as a search needs it
                os.remove(os.path.join(self.tmp, f"{ids[0]}.fa"))
                self._writeb("renamed.fa", ">renamed\nACGT\n")
                self.assertEqual(gl.describe_header_problems(self.tmp), [])

    def test_an_id_that_matches_no_other_name_is_taken(self):
        # (c[12] matches neither cA nor itself; a '[' with no ']' is a character; [!a]
        # does not match xa; a hidden file is never a match)
        for ids in (("c[12]", "cA"), ("HLA-A*01:01", "chr1"), ("HLA-A*01:01:01:01", "HLA-A*02:01:01:01"),
                    ("c[12", "c1"), ("x[!a]", "xa")):
            with self.subTest(ids):
                self._fresh()
                for rid in ids:
                    self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
                self.assertEqual(gl.find_pattern_clashes(self.tmp), {})
                self.assertEqual(gl.describe_header_problems(self.tmp), [])
                self.assertEqual(vi.check_genome_fasta(self.tmp)[0], [])
        self._fresh()
        self._writeb("a*.fa", ">a*\nACGT\n")
        self._writeb("._ab.fa", "\0\5\26\7")  # (macOS's file attributes of ab.fa)
        self.assertEqual(gl.find_pattern_clashes(self.tmp), {})
        self._writeb("ab.fasta", ">ab\nACGT\n")  # (which no search reads)
        self.assertEqual(gl.find_pattern_clashes(self.tmp), {})

    def test_each_file_of_such_an_id_is_named(self):
        # HLA ids of the GATK GRCh38 FASTA, 10 of whose 525 match another's name so
        for rid in ("HLA-A*01:02", "HLA-A*29:02:01:02", "HLA-A*30:02:01:02", "HLA-A*31:01:02",
                    "HLA-A*68:01:01:02", "HLA-A*02:10", "HLA-A*24:02:10", "HLA-A*01:01:01:01"):
            self._writeb(f"{rid}.fa", f">{rid}\nACGT\n")
        self.assertEqual(gl.find_pattern_clashes(self.tmp), {
            "HLA-A*01:02.fa": ["HLA-A*29:02:01:02.fa", "HLA-A*30:02:01:02.fa", "HLA-A*31:01:02.fa",
                               "HLA-A*68:01:01:02.fa"],
            "HLA-A*02:10.fa": ["HLA-A*24:02:10.fa"]})
        problems = gl.describe_header_problems(self.tmp)
        self.assertEqual(len(problems), 1)
        self.assertIn("HLA-A*01:02.fa: its sequence id 'HLA-A*01:02', read as a shell pattern, "
                      "matches HLA-A*29:02:01:02.fa, HLA-A*30:02:01:02.fa, HLA-A*31:01:02.fa, ...; "
                      "2 files here hold such an id (HLA-A*01:02.fa, HLA-A*02:10.fa). ", problems[0])

    def test_a_file_of_several_sequences_is_left_to_its_split(self):
        # its split names the files after its sequences, and it is reported as such
        self._writeb("a*.fa", ">chr1\nACGT\n>chr2\nGG\n")
        self._writeb("ab.fa", ">ab\nACGT\n")
        self.assertEqual(list(gl.find_pattern_clashes(self.tmp)), ["a*.fa"])
        self.assertEqual(gl.describe_header_problems(self.tmp, reported={"a*.fa": ["chr1", "chr2"]}), [])
        issues, _ = vi.check_genome_fasta(self.tmp)
        self.assertEqual([i.severity for i in issues], [vi.ERROR])
        self.assertIn("a*.fa: holds 2 sequences", issues[0].message)

    def test_an_index_built_from_such_a_header_is_found(self):
        # it holds those names once the FASTA reads >id: every search of it would still
        # lose its results, and a build took it as present
        lib = os.path.join(self.tmp, "genome_library")
        for name, bins in (("NRG_2_ecoli", ["NRG_NC_000913.3 Escherichia coli_1.bin"]),
                           ("NRG_2_ecoli+ds", ["NRG_chr1\r_1.bin"]),
                           ("NRG_2_clean", ["NRG_NC_000913.3 Escherichia coli_1.bin"]),
                           ("NRG_2_my_ecoli", ["NRG_NC_000913.3 Escherichia coli_1.bin"]),
                           ("NRG_2_ecoli2", ["NRG_chr1_1.bin", ".pam_build"])):
            os.makedirs(os.path.join(lib, name))
            for b in bins:
                _write(os.path.join(lib, name, b), "")
        problems = gl.describe_index_problems(lib, "ecoli")
        self.assertEqual(len(problems), 2)
        self.assertIn(f"The index {lib}/NRG_2_ecoli holds 'NRG_NC_000913.3 Escherichia coli_1.bin'",
                      problems[0])
        self.assertIn(f"rm -rf '{lib}/NRG_2_ecoli'", problems[0])
        self.assertIn(f"rm -rf '{lib}/NRG_2_ecoli+ds'", problems[1])
        self.assertIn("and the variant genome Genomes/ecoli+ds and run build-index-only", problems[1])
        self.assertNotIn("variant genome", problems[0])
        self.assertEqual(gl.describe_index_problems(lib, "ecoli2"), [])
        self.assertEqual(gl.describe_index_problems(lib + "_none", "ecoli"), [])


class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _run(self, *args):
        return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)

    def test_clean_genome_exits_zero(self):
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nAC\n")
        self.assertEqual(self._run(self.tmp).returncode, 0)

    def test_multi_sequence_genome_exits_one_and_says_how_to_fix(self):
        _write(os.path.join(self.tmp, "susScr11.fa"), MULTI)
        r = self._run(self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertIn("ONE FASTA file per sequence", r.stderr)
        self.assertIn("--split", r.stderr)

    def test_split_flag_fixes_in_place(self):
        _write(os.path.join(self.tmp, "susScr11.fa"), MULTI)
        r = self._run("--split", self.tmp)
        self.assertEqual(r.returncode, 0)
        self.assertEqual(self._run(self.tmp).returncode, 0)
        self.assertEqual(sorted(os.listdir(self.tmp)), ["chr1.fa", "chr2.fa", "chrUn_x_random.fa"])

    def test_split_flag_names_a_file_of_one_sequence_after_it(self):
        _write(os.path.join(self.tmp, "ecoli.fa"), ">NC_000913.3\nAC\n")
        # the check a search runs first, and which a search that dropped all its rows says
        # to run, passed it: a search of it finds that sequence's off-targets and drops
        # them all (M1). It is refused, saying how to fix it
        r = self._run(self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(r.stderr.splitlines()[:2], [
            "ERROR: ecoli.fa: holds the sequence NC_000913.3 and is named otherwise. CRISPRme "
            "reads a sequence from the file named after it (NC_000913.3.fa), so a search of "
            "this folder finds that sequence's off-targets and then drops them all, with no "
            "error. Name each file after its sequence with:",
            f"    python PostProcess/genome_layout.py --split {self.tmp}"])
        r = self._run("--split", self.tmp)
        self.assertEqual(r.returncode, 0)
        self.assertEqual(r.stdout, f"Renamed ecoli.fa to NC_000913.3.fa in {self.tmp}\n")
        self.assertEqual(os.listdir(self.tmp), ["NC_000913.3.fa"])
        self.assertEqual(self._run(self.tmp).returncode, 0)

    def test_a_header_of_more_than_the_id_exits_one_and_split_writes_it(self):
        # what a web search runs (the search script), which has no validator
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1 1\nAC\n")
        r = self._run(self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertIn("ERROR: chr1.fa: its header line reads '>chr1 1'", r.stderr)
        self.assertIn(f"genome_layout.py --split {self.tmp}", r.stderr)
        r = self._run("--split", self.tmp)
        self.assertEqual(r.returncode, 0)
        self.assertIn("Wrote chr1.fa's header as >id and its lines with no CR", r.stdout)
        self.assertIn("any CRISPRitz index built from the old layout must be rebuilt", r.stdout)
        self.assertEqual(_read(os.path.join(self.tmp, "chr1.fa")), ">chr1\nAC\n")
        self.assertEqual(self._run(self.tmp).returncode, 0)
        self.assertEqual(self._run("--split", self.tmp).stdout, "")

    def test_deep_check_finds_a_cr_past_the_first_lines(self):
        with open(os.path.join(self.tmp, "chr1.fa"), "wb") as f:
            f.write(b">chr1\nAC\nGT\r\nAC\n")
        self.assertEqual(self._run(self.tmp).returncode, 0)
        r = self._run("--deep", self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertIn("ERROR: chr1.fa: its lines end in CR", r.stderr)

    def test_a_cr_that_only_ends_the_file_is_refused_too(self):
        # the CRLF of its last line, or CRs after its last LF: the index build takes a CR on
        # a line that ends in LF for a base, and the variant enricher keeps every CR, so with
        # a PAM in the last bases of the sequence the build wrote warnings on stderr and a
        # search of the index crashed (SIGSEGV) or reported other rows, through the real
        # engine (integration round 2 review, F1). The two-line check refuses a CR in its two
        # lines, --deep one anywhere, and --split writes the file with none
        path = os.path.join(self.tmp, "chr1.fa")
        # (the file, refused by the two-line check, as --split writes it)
        for text, fast, fixed in ((b">chr1\nAC\nGT\r\n", False, b">chr1\nAC\nGT\n"),
                                  (b">chr1\nAC\nGT\n\r", False, b">chr1\nAC\nGT\n"),
                                  (b">chr1\nAC\nGT\n\r\r\r", False, b">chr1\nAC\nGT\n"),
                                  (b">chr1\nAC\nGT\r\n\r\n", False, b">chr1\nAC\nGT\n\n"),
                                  (b">chr1\nAC\nGT\r\r\n", False, b">chr1\nAC\nGT\n"),
                                  (b">chr1\nAC\nGT\r\n\n", False, b">chr1\nAC\nGT\n\n"),
                                  (b">chr1\nAC\nGT\n\r\n", False, b">chr1\nAC\nGT\n\n"),
                                  (b">chr1\nACGT\r\n", True, b">chr1\nACGT\n"),
                                  (b">chr1\nACGT\r", True, b">chr1\nACGT"),
                                  (b">chr1\nACGT\r\nAC\n", True, b">chr1\nACGT\nAC\n"),
                                  (b">chr1\nAC\nGT\r\n\nAC\n", False, b">chr1\nAC\nGT\n\nAC\n")):
            with open(path, "wb") as f:
                f.write(text)
            for options, refused in (((), fast), (("--deep",), True)):
                with self.subTest(text=text, options=options):
                    r = self._run(*options, self.tmp)
                    if refused:
                        self.assertEqual(r.returncode, 1, r.stderr)
                        self.assertIn("ERROR: chr1.fa: its lines end in CR", r.stderr)
                    else:
                        self.assertEqual((r.returncode, r.stderr), (0, ""))
            with self.subTest(text=text, check="complete-search's"):
                issues, _ = vi.check_genome_fasta(self.tmp)
                self.assertEqual([i.severity for i in issues], [vi.ERROR] if fast else [])
            with self.subTest(text=text, options="--split"):
                r = self._run("--split", self.tmp)
                self.assertEqual((r.returncode, r.stderr), (0, ""))
                self.assertIn("Wrote chr1.fa's header as >id and its lines with no CR", r.stdout)
                with open(path, "rb") as f:
                    self.assertEqual(f.read(), fixed)
                self.assertEqual(self._run("--deep", self.tmp).returncode, 0)

    def test_the_check_of_a_genome_finds_an_index_left_from_its_old_headers(self):
        # the web search's guard: a data folder, its Genomes/<genome>, its genome_library
        genome = os.path.join(self.tmp, "Genomes", "ecoli")
        os.makedirs(genome)
        _write(os.path.join(genome, "NC_000913.3.fa"), ">NC_000913.3\nAC\n")
        self.assertEqual(self._run(genome).returncode, 0)
        index = os.path.join(self.tmp, "genome_library", "NRG_2_ecoli")
        os.makedirs(index)
        _write(os.path.join(index, "NRG_NC_000913.3 Escherichia coli_1.bin"), "")
        r = self._run(genome)
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"ERROR: The index {index} holds", r.stderr)

    def test_split_names_the_indexes_left_from_the_old_headers(self):
        genome = os.path.join(self.tmp, "Genomes", "ecoli")
        os.makedirs(genome)
        _write(os.path.join(genome, "NC_000913.3.fa"), ">NC_000913.3 E. coli\nAC\n")
        library = os.path.join(self.tmp, "genome_library")
        for name in ("NRG_2_ecoli", "NRG_2_ecoli+ds", "NRG_2_other"):
            os.makedirs(os.path.join(library, name))
        r = self._run("--split", genome)
        self.assertEqual(r.returncode, 0)
        self.assertIn(f"NOTE: in {library}, built on this genome: NRG_2_ecoli, NRG_2_ecoli+ds\n",
                      r.stdout)

    def test_split_of_an_id_it_cannot_take_says_so_and_changes_nothing(self):
        _write(os.path.join(self.tmp, "asm.fa"), ">gi|49175990|ref|NC_000913.2| E. coli\nAC\n")
        r = self._run("--split", self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(r.stderr, "ERROR: asm.fa: sequence id 'gi|49175990|ref|NC_000913.2|' holds "
                         "'|', which the pipeline cannot take; rename the sequence in the FASTA "
                         "header and retry\n")
        self.assertEqual(os.listdir(self.tmp), ["asm.fa"])

    def test_an_id_that_matches_another_name_as_a_shell_pattern_exits_one(self):
        # what a web search runs (the search script), which has no validator
        for rid in ("HLA-A*02:10", "HLA-A*24:02:10"):
            _write(os.path.join(self.tmp, f"{rid}.fa"), f">{rid}\nAC\n")
        r = self._run(self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertTrue(r.stderr.startswith(
            "ERROR: HLA-A*02:10.fa: its sequence id 'HLA-A*02:10', read as a shell pattern, "
            "matches HLA-A*24:02:10.fa. "), r.stderr)
        os.rename(os.path.join(self.tmp, "HLA-A*24:02:10.fa"), os.path.join(self.tmp, "chr1.fa"))
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nAC\n")
        r = self._run(self.tmp)
        self.assertEqual((r.returncode, r.stderr), (0, ""))

    def test_bad_usage(self):
        self.assertEqual(self._run().returncode, 2)
        self.assertEqual(self._run("/no/such/dir").returncode, 2)


class TestDownloadSplits(unittest.TestCase):
    """download_reference_genome must never leave a multi-sequence file behind."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.genomes = os.path.join(self.tmp, "Genomes")

    def _fake_download(self, payload: bytes):
        def fake(staging, http_url=None, **_kw):
            gz = os.path.join(staging, "asm.fa.gz")
            with gzip.open(gz, "wb") as f:
                f.write(payload)
            return gz
        return fake

    def test_single_multi_fasta_url_download_is_split(self):
        with mock.patch.object(utils, "download", self._fake_download(MULTI.encode())):
            dest = utils.download_reference_genome(
                "susScr11", self.genomes, source="url", url="http://x/asm.fa.gz"
            )
        self.assertEqual(
            sorted(os.listdir(dest)), ["chr1.fa", "chr2.fa", "chrUn_x_random.fa"]
        )
        self.assertEqual(gl.find_multi_record_fastas(dest, deep=True), {})

    def test_url_download_of_a_header_with_a_description_is_written_as_the_id(self):
        for payload, name in ((b">NC_000913.3 E. coli K-12\nACGT\n", "NC_000913.3.fa"),
                              (b">chr1 1\r\nACGT\r\n", "chr1.fa")):
            with self.subTest(name):
                with mock.patch.object(utils, "download", self._fake_download(payload)):
                    dest = utils.download_reference_genome(
                        "asm", self.genomes, source="url", url="http://x/asm.fa.gz"
                    )
                self.assertEqual(os.listdir(dest), [name])
                self.assertEqual(_read(os.path.join(dest, name)), f">{name[:-3]}\nACGT\n")

    def test_single_sequence_url_download_is_named_after_it(self):
        # it went in as <assembly>.fa, and every search of it dropped all its results
        with mock.patch.object(utils, "download", self._fake_download(b">NC_000913.3\nACGT\n")):
            dest = utils.download_reference_genome(
                "ecoli", self.genomes, source="url", url="http://x/asm.fa.gz"
            )
        self.assertEqual(os.listdir(dest), ["NC_000913.3.fa"])

    def test_failed_split_leaves_no_half_built_genome_folder(self):
        dup = b">a\nAC\n>a\nGG\n"
        with mock.patch.object(utils, "download", self._fake_download(dup)):
            with self.assertRaises(ValueError):
                utils.download_reference_genome(
                    "badgenome", self.genomes, source="url", url="http://x/asm.fa.gz"
                )
        self.assertFalse(os.path.exists(os.path.join(self.genomes, "badgenome")))


class TestDownloadRefusesUnsafeArchive(unittest.TestCase):
    """--source url extracts an archive from a URL the web client chose: a member
    that is absolute, has a '..' component, or is a link/device/FIFO refuses the
    whole archive before anything is extracted; a chromFa archive still works, its
    files with safe modes and none of the archive's owners."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.data = os.path.join(self.tmp, "data")
        self.genomes = os.path.join(self.data, "Genomes")
        os.makedirs(self.genomes)  # as in a real data dir (the function makes it anyway)
        os.makedirs(os.path.join(self.data, "Settings"))

    def _fake_download(self, members):
        def fake(staging, http_url=None, **_kw):
            path = os.path.join(staging, "asm.tar.gz")
            with tarfile.open(path, "w:gz") as tar:
                for info, payload in members:
                    tar.addfile(info, io.BytesIO(payload) if payload else None)
            return path
        return fake

    @staticmethod
    def _file(name, data=b">chr1\nACGT\n"):
        info = tarfile.TarInfo(name)
        info.size = len(data)
        return info, data

    @staticmethod
    def _special(name, kind, linkname=""):
        info = tarfile.TarInfo(name)
        info.type, info.linkname = kind, linkname
        return info, None

    def _snapshot(self):
        return sorted(
            os.path.relpath(os.path.join(d, n), self.tmp)
            for d, dirs, files in os.walk(self.tmp) for n in dirs + files
        )

    def test_unsafe_member_refused_before_anything_is_extracted(self):
        outside = os.path.join(self.tmp, "outside.txt")
        secret = os.path.join(self.tmp, "secret.txt")  # a hardlink target in the sandbox
        with open(secret, "w") as fh:
            fh.write("ORIGINAL\n")
        cases = {
            "'..'": self._file("../../Settings/set_EVILJOB000/cmd.json", b"{}"),
            "absolute": self._file(outside, b"x"),
            "symlink": self._special("lnk", tarfile.SYMTYPE, self.tmp),
            "hardlink": self._special("hl", tarfile.LNKTYPE, secret),
            "fifo": self._special("pipe", tarfile.FIFOTYPE),
        }
        for label, bad in cases.items():
            with self.subTest(label):
                before = self._snapshot()
                fake = self._fake_download([self._file("chr1.fa"), bad])
                with mock.patch.object(utils, "download", fake):
                    with self.assertRaises(ValueError):
                        utils.download_reference_genome(
                            "evil", self.genomes, source="url", url="http://x/asm.tar.gz")
                self.assertEqual(self._snapshot(), before)

    def test_chromfa_archive_still_installs(self):
        members = [self._file("chroms/chr1.fa"), self._file("chroms/chr2.fa", b">chr2\nGG\n")]
        with mock.patch.object(utils, "download", self._fake_download(members)):
            dest = utils.download_reference_genome(
                "myasm", self.genomes, source="url", url="http://x/asm.tar.gz")
        self.assertEqual(sorted(os.listdir(dest)), ["chr1.fa", "chr2.fa"])

    def test_files_go_in_with_safe_modes_and_none_of_the_archive_s_owners(self):
        # tarfile gave each file the mode its member has, setuid and setgid too, and, run
        # as root (as the Docker image runs the web app, whose data folder is mounted from
        # the host), the owner it names
        def owned(name, mode, data=None):
            info, data = (self._file(name, data) if data else
                          self._special(name, tarfile.DIRTYPE))
            info.mode, info.uid, info.gid = mode, 1234, 1234
            return info, data

        members = [owned("chroms", 0o3777), owned("chroms/chr1.fa", 0o4755, b">chr1\nACGT\n"),
                   owned("chroms/chr2.fa", 0o2666, b">chr2\nGG\n")]
        for fallback in (False, True):  # (True: Python before 3.11.4, with no 'data' filter)
            with self.subTest(fallback=fallback):
                owners = []
                record = lambda path, uid, gid, *a, **kw: owners.append((path, uid, gid))  # noqa: E731
                with mock.patch.object(utils, "download", self._fake_download(members)), \
                        mock.patch("os.geteuid", return_value=0), \
                        mock.patch("os.chown", record), mock.patch("os.lchown", record), \
                        mock.patch.object(safe_tar, "_DATA_FILTER",
                                          None if fallback else safe_tar._DATA_FILTER), \
                        warnings.catch_warnings():
                    warnings.simplefilter("ignore", DeprecationWarning)
                    dest = utils.download_reference_genome(
                        "modes", self.genomes, source="url", url="http://x/asm.tar.gz")
                modes = {name: stat.S_IMODE(os.stat(os.path.join(dest, name)).st_mode)
                         for name in os.listdir(dest)}
                self.assertEqual(modes, {"chr1.fa": 0o755, "chr2.fa": 0o644})
                self.assertEqual([o for o in owners if o[1:] != (-1, -1)], [])
                self.assertEqual(_read(os.path.join(dest, "chr2.fa")), ">chr2\nGG\n")

    def test_an_archive_goes_in_as_the_settings_guide_says(self):
        # docs/SETTINGS_DATA_MANAGER.md, Direct URL: each header written as >id; unlike an
        # upload, a .fasta file that needs nothing stays so, and a hidden file goes in
        # (the guide said once it was laid out as an upload is; change both together)
        members = [self._file("asm/chr1.fasta"), self._file("asm/._chr2.fa", b"\0\5\26\7"),
                   self._file("asm/chr2.fa", b">chr2\nGG\n"),
                   self._file("asm/chr3.fa", b">chr3 3\nCC\n"),
                   self._file("asm/hg38_chr4.fa", b">chr4\nTT\n")]
        with mock.patch.object(utils, "download", self._fake_download(members)):
            dest = utils.download_reference_genome(
                "uasm", self.genomes, source="url", url="http://x/asm.tar.gz")
        self.assertEqual(sorted(os.listdir(dest)),
                         ["._chr2.fa", "chr1.fasta", "chr2.fa", "chr3.fa", "chr4.fa"])
        self.assertEqual(_read(os.path.join(dest, "chr3.fa")), ">chr3\nCC\n")
        self.assertEqual(_read(os.path.join(dest, "chr4.fa")), ">chr4\nTT\n")

    def test_a_fasta_whose_only_cr_ends_it_goes_in_with_none(self):
        # files named after their sequence, with a plain header: only their CR is to write.
        # A build reads such a CR, or a variant build does, and a search of the index
        # crashed or reported other rows (integration round 2 review, F1)
        want = {"chr1.fa": (b">chr1\nAC\nGT\r\n", b">chr1\nAC\nGT\n"),
                "chr2.fa": (b">chr2\nGG\n\r", b">chr2\nGG\n"),
                "chr3.fa": (b">chr3\nCC\r\n\r\n", b">chr3\nCC\n\n")}
        members = [self._file(f"asm/{name}", data) for name, (data, _) in want.items()]
        with mock.patch.object(utils, "download", self._fake_download(members)):
            dest = utils.download_reference_genome(
                "crasm", self.genomes, source="url", url="http://x/asm.tar.gz")
        self.assertEqual(sorted(os.listdir(dest)), sorted(want))
        for name, (_, written) in want.items():
            with open(os.path.join(dest, name), "rb") as fh:
                self.assertEqual(fh.read(), written, name)

    def test_an_archive_with_two_fasta_files_of_one_name_is_refused(self):
        # each went in under its file name alone, the second over the first, with exit
        # status 0 and no message: an assembly's two haplotypes went in mixed (finB round
        # 3 review). Refused, as an upload is, and a genome installed under the name is kept
        installed = os.path.join(self.genomes, "udup")
        os.makedirs(installed)
        _write(os.path.join(installed, "chr9.fa"), ">chr9\nACGT\n")
        members = [self._file("asm/a/chr1.fa", b">chr1\nAAAA\n"),
                   self._file("asm/b/chr1.fa", b">chr1\nCCCC\n"),
                   self._file("asm/chr3.fa", b">chr3\nGG\n")]
        with mock.patch.object(utils, "download", self._fake_download(members)):
            with self.assertRaises(RuntimeError) as cm:
                utils.download_reference_genome(
                    "udup", self.genomes, source="url", url="http://x/asm.tar.gz")
        message = str(cm.exception)
        self.assertTrue(message.startswith(
            "The download holds two FASTA files that would both go in as chr1.fa ("), message)
        self.assertIn("asm/a/chr1.fa", message)
        self.assertIn("asm/b/chr1.fa", message)
        self.assertIn("), one replacing the other, so nothing was installed. If they belong "
                      "to two genomes (the two haplotypes of an assembly, say), install each "
                      "under a name of its own.", message)
        self.assertEqual(os.listdir(self.genomes), ["udup"])  # (no staging folder left)
        self.assertEqual(os.listdir(installed), ["chr9.fa"])
        # one of each name goes in as before
        members[1] = self._file("asm/b/chr2.fa", b">chr2\nCCCC\n")
        with mock.patch.object(utils, "download", self._fake_download(members)):
            dest = utils.download_reference_genome(
                "udup", self.genomes, source="url", url="http://x/asm.tar.gz")
        self.assertEqual(sorted(os.listdir(dest)), ["chr1.fa", "chr2.fa", "chr3.fa"])


class TestARefusedDownloadKeepsTheInstalledGenome(unittest.TestCase):
    """A download that its layout refuses (an id the pipeline cannot take, two sequences
    of one id, a sequence that would go in as another file) leaves a genome installed
    under its name as it was. The download removed that genome before it laid out its
    own files there, so a refused one left no genome under the name; at a5379e0 it put
    the bad genome in its place (finB round 3 review). It is laid out in the staging
    folder now, as an upload is. By Direct URL (a .fa.gz, a .tar.gz) and from UCSC (a
    chromFa archive, or the one FASTA of an assembly that has none)."""

    SOURCES = ("url .fa.gz", "url .tar.gz", "ucsc chromFa", "ucsc .fa.gz")
    INSTALLED = ">chr9\nACGT\n"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.genomes = os.path.join(self.tmp, "Genomes")
        self.installed = os.path.join(self.genomes, "g")

    def install(self):
        shutil.rmtree(self.genomes, ignore_errors=True)
        os.makedirs(self.installed)
        _write(os.path.join(self.installed, "chr9.fa"), self.INSTALLED)

    def download(self, source, files):
        """download_reference_genome of g from source, whose archive holds files (path in
        the archive -> content; a .fa.gz, the content of its one FASTA)."""
        def fake(staging, http_url=None, **_kw):
            if source.endswith(".fa.gz"):
                path = os.path.join(staging, "g.fa.gz")
                (content,) = files.values()
                with gzip.open(path, "wb") as fh:
                    fh.write(content)
                return path
            path = os.path.join(staging, "g.tar.gz")
            with tarfile.open(path, "w:gz") as tar:
                for name, content in files.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(content)
                    tar.addfile(info, io.BytesIO(content))
            return path

        def exists(url):  # (UCSC: its chromFa archive, or the one FASTA of one with none)
            return url.endswith(".chromFa.tar.gz" if source == "ucsc chromFa" else "/g.fa.gz")

        kind, ext = source.split()
        with mock.patch.object(utils, "download", fake), \
                mock.patch.object(utils, "_url_exists", exists):
            return utils.download_reference_genome(
                "g", self.genomes, source=kind, url="http://x/g" + ext if kind == "url" else None)

    def assert_kept(self):
        self.assertEqual(os.listdir(self.genomes), ["g"])  # (no staging folder left)
        self.assertEqual(os.listdir(self.installed), ["chr9.fa"])
        self.assertEqual(_read(os.path.join(self.installed, "chr9.fa")), self.INSTALLED)

    def test_a_refused_download_keeps_the_installed_genome(self):
        one = {  # (a .fa.gz: one FASTA, named g.fa)
            "an id that starts with '-'": (
                b">-x\nACGT\n>b\nGG\n",
                "g.fa: sequence id '-x' starts with '-', which the pipeline cannot take"),
            "two sequences of one id": (b">a\nAC\n>a\nGG\n", "g.fa: duplicate sequence id 'a'"),
        }
        archive = {
            "an id that starts with '-'": (
                {"asm/-x.fa": b">-x a description\nACGT\n", "asm/b.fa": b">b\nGG\n"},
                "-x.fa: sequence id '-x' starts with '-', which the pipeline cannot take"),
            "two sequences of one id": (
                {"asm/all.fa": b">a\nAC\n>a\nGG\n"}, "all.fa: duplicate sequence id 'a'"),
            "a sequence that would go in as another file": (
                {"asm/chr1.fa": b">chr1\nAC\n", "asm/all.fa": b">chr1\nGG\n>chr2\nTT\n"},
                "all.fa: cannot split -- chr1.fa already exists"),  # (no staging folder)
        }
        for source in self.SOURCES:
            cases = one if source.endswith(".fa.gz") else archive
            for label, (content, says) in cases.items():
                with self.subTest(source=source, case=label):
                    self.install()
                    files = {"g.fa": content} if source.endswith(".fa.gz") else content
                    with self.assertRaises(ValueError) as cm:
                        self.download(source, files)
                    self.assertTrue(str(cm.exception).startswith(says), str(cm.exception))
                    self.assert_kept()
            with self.subTest(source=source, case="a download that goes in"):
                # replaces the installed genome, as before
                self.install()
                files = ({"g.fa": b">chr1 one\nAC\n>chr2\nGG\n"} if source.endswith(".fa.gz")
                         else {"asm/chr1.fa": b">chr1 one\nAC\n", "asm/chr2.fa": b">chr2\nGG\n"})
                self.assertEqual(self.download(source, files), self.installed)
                self.assertEqual(os.listdir(self.genomes), ["g"])
                self.assertEqual(sorted(os.listdir(self.installed)), ["chr1.fa", "chr2.fa"])
                self.assertEqual(_read(os.path.join(self.installed, "chr1.fa")), ">chr1\nAC\n")

    def test_a_download_that_cannot_be_put_in_place_puts_the_installed_genome_back(self):
        real = os.rename

        def rename(src, dst, **kw):  # (the laid-out genome's move into place fails)
            if dst == self.installed and os.path.basename(src) == "new":
                raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC))
            return real(src, dst, **kw)

        self.install()
        with mock.patch.object(os, "rename", rename):
            with self.assertRaises(OSError):
                self.download("url .fa.gz", {"g.fa": b">chr1\nAC\n"})
        self.assert_kept()

    def test_an_installed_genome_that_is_a_link_to_a_shared_copy(self):
        # the link alone is replaced, and kept when the download is refused; the copy it
        # points to is never written to (it was left as it was, the download failing)
        shared = os.path.join(self.tmp, "shared", "g")
        os.makedirs(shared)
        _write(os.path.join(shared, "chr9.fa"), self.INSTALLED)
        os.makedirs(self.genomes)
        os.symlink(shared, self.installed)
        with self.assertRaises(ValueError):
            self.download("url .fa.gz", {"g.fa": b">-x\nACGT\n"})
        self.assertEqual(os.readlink(self.installed), shared)
        self.assert_kept()
        self.assertEqual(self.download("url .fa.gz", {"g.fa": b">chr1\nAC\n"}), self.installed)
        self.assertFalse(os.path.islink(self.installed))
        self.assertEqual(os.listdir(self.installed), ["chr1.fa"])
        self.assertEqual(os.listdir(shared), ["chr9.fa"])
        self.assertEqual(_read(os.path.join(shared, "chr9.fa")), self.INSTALLED)


class TestValidator(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_multi_sequence_file_is_an_error_with_the_fix(self):
        _write(os.path.join(self.tmp, "susScr11.fa"), MULTI)
        issues, _ = vi.check_genome_fasta(self.tmp)
        errors = [i for i in issues if i.severity == vi.ERROR]
        self.assertEqual(len(errors), 1)
        self.assertIn("3 sequences", errors[0].message)
        self.assertIn("genome_layout.py --split", errors[0].message)

    def test_wrong_name_single_sequence_is_refused(self):
        # a search of it drops all chr1's off-targets, and its own check (genome_layout.py)
        # refuses it now, as here: this was only a warning (M1). A .fasta file, which no
        # search reads, is still only named
        _write(os.path.join(self.tmp, "NC_000001.fa"), ">chr1\nACGT\n")
        issues, _ = vi.check_genome_fasta(self.tmp)
        self.assertEqual([i.severity for i in issues], [vi.ERROR])
        self.assertTrue(issues[0].message.startswith(
            "NC_000001.fa: holds the sequence chr1 and is named otherwise."))
        os.rename(os.path.join(self.tmp, "NC_000001.fa"), os.path.join(self.tmp, "NC_000001.fasta"))
        _write(os.path.join(self.tmp, "chr2.fa"), ">chr2\nACGT\n")
        issues, _ = vi.check_genome_fasta(self.tmp)
        self.assertEqual([i.severity for i in issues], [vi.WARN])
        self.assertIn("chr1", issues[0].message)

    def test_the_search_check_reports_an_index_left_from_old_headers(self):
        genome = os.path.join(self.tmp, "Genomes", "ecoli")
        os.makedirs(genome)
        _write(os.path.join(genome, "chr1.fa"), ">chr1\nACGT\n")
        index = os.path.join(self.tmp, "genome_library", "NRG_2_ecoli")
        os.makedirs(index)
        _write(os.path.join(index, "NRG_chr1 1_1.bin"), "")
        guide = os.path.join(self.tmp, "guide.txt")
        _write(guide, "CTAACAGTTGCTTTTATCACNNN\n")
        pam = os.path.join(self.tmp, "pam.txt")
        _write(pam, "NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        report = vi.run_lightweight(genome, [], guide, pam, current_working_directory=self.tmp)
        self.assertTrue(report.has_errors())
        self.assertIn(f"[ERROR] The index {index} holds 'NRG_chr1 1_1.bin'", report.render())

    def test_standard_genome_has_no_issues(self):
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nACGT\n")
        _write(os.path.join(self.tmp, "chr2.fa"), ">chr2\nACGT\n")
        issues, chroms = vi.check_genome_fasta(self.tmp)
        self.assertEqual(issues, [])
        self.assertEqual(sorted(chroms), ["chr1", "chr2"])


if __name__ == "__main__":
    unittest.main()
