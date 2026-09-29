"""Tests for the one-FASTA-per-sequence genome layout (``genome_layout.py``) and
the places that enforce it (``utils.download_reference_genome``,
``validate_inputs.check_genome_fasta``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_genome_layout.py' -v

Background: CRISPRme keys its chromosome list by FASTA file name, so a genome
delivered as ONE multi-sequence file (UCSC's ``susScr11.fa.gz``) indexes and
searches fine and then silently returns an empty result set.
"""

import gzip
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import genome_layout as gl  # noqa: E402
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
        self.assertEqual(joined, MULTI)  # byte-for-byte
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

    def test_failed_split_leaves_no_half_built_genome_folder(self):
        dup = b">a\nAC\n>a\nGG\n"
        with mock.patch.object(utils, "download", self._fake_download(dup)):
            with self.assertRaises(ValueError):
                utils.download_reference_genome(
                    "badgenome", self.genomes, source="url", url="http://x/asm.fa.gz"
                )
        self.assertFalse(os.path.exists(os.path.join(self.genomes, "badgenome")))


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

    def test_wrong_name_single_sequence_is_only_a_warning(self):
        _write(os.path.join(self.tmp, "NC_000001.fa"), ">chr1\nACGT\n")
        issues, _ = vi.check_genome_fasta(self.tmp)
        self.assertEqual([i.severity for i in issues], [vi.WARN])
        self.assertIn("chr1", issues[0].message)

    def test_standard_genome_has_no_issues(self):
        _write(os.path.join(self.tmp, "chr1.fa"), ">chr1\nACGT\n")
        _write(os.path.join(self.tmp, "chr2.fa"), ">chr2\nACGT\n")
        issues, chroms = vi.check_genome_fasta(self.tmp)
        self.assertEqual(issues, [])
        self.assertEqual(sorted(chroms), ["chr1", "chr2"])


if __name__ == "__main__":
    unittest.main()
