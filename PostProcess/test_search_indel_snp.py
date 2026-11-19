#!/usr/bin/env python3
"""What a search that enriches a dataset itself does for SNP+indel targets
(search_indel_snp.py): the SNP overlay of the fake-indel genome, and the phased indel
genotype store. Run through its command line, as the search script runs it: its exit
status and stderr are what a stage reads."""

import gzip
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import search_indel_snp  # noqa: E402

_TOOL = Path(__file__).resolve().parent / "search_indel_snp.py"
SEQ = "ACGTTGCAAG" * 6
PLAIN = f">fakechr1\n{SEQ[:26]}T{SEQ[26:]}\n"
CODED = f">fakechr1\n{SEQ[:10]}R{SEQ[11:26]}T{SEQ[26:]}\n"
LOG = ("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
       f"chr1_0-60_1\tS2\trs1\t0.25\tchr1_26_{SEQ[25]}_{SEQ[25]}T\t0,61\t{SEQ}\n")
VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"


def tool(*args, env=None):
    return subprocess.run([sys.executable, str(_TOOL), *map(str, args)], capture_output=True,
                          text=True, env=dict(os.environ, **(env or {})))


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)


class TheOverlay(_Tmp):
    def setUp(self):
        super().setUp()
        self.fake = self.root / "indels" / "fake_ds_chr1" / "fakechr1.fa"
        self.fake.parent.mkdir(parents=True)
        self.enriched = self.root / "enriched"
        self.enriched.mkdir()
        (self.enriched / "chr1.enriched.fa").write_text(f">chr1\n{SEQ[:10]}R{SEQ[11:]}\n")
        self.logs = self.root / "logs"
        self.logs.mkdir()
        (self.logs / "logchr1.txt").write_text(LOG)

    def overlay(self):
        return tool("overlay", self.root / "indels", self.enriched, self.logs)

    def test_it_puts_the_codes_on_and_exits_0_however_many(self):
        # (the module's own command line exits with the number of codes written)
        self.fake.write_text(PLAIN)
        out = self.overlay()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertEqual(self.fake.read_text().split("\n")[1], CODED.split("\n")[1])

    def test_a_genome_coded_already_is_not_written(self):
        self.fake.write_text(CODED)
        before = self.fake.stat().st_mtime_ns
        out = self.overlay()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertEqual(self.fake.stat().st_mtime_ns, before)

    def test_so_a_read_only_one_coded_already_is_no_failure(self):
        self.fake.write_text(CODED)
        self.fake.parent.chmod(0o555)
        self.addCleanup(self.fake.parent.chmod, 0o755)
        out = self.overlay()
        self.assertEqual((out.returncode, out.stderr), (0, ""))

    def test_a_chromosome_without_its_inputs_is_skipped(self):
        # complete-test --chrom's view holds a chromosome's fake-indel genome, and no
        # enriched genome
        self.fake.write_text(PLAIN)
        (self.enriched / "chr1.enriched.fa").unlink()
        out = self.overlay()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertIn("skip chr1 (missing inputs)", out.stdout)
        self.assertEqual(self.fake.read_text(), PLAIN)

    def test_a_failure_says_so_on_stderr(self):
        out = tool("overlay", self.root / "no_such_folder", self.enriched, self.logs)
        self.assertNotEqual(out.returncode, 0)
        self.assertTrue(out.stderr)

    def test_a_genome_it_cannot_write_stops_in_one_line(self):
        # (the codes are not on it, so an index built from it would be blind.) Patched, as
        # a root run writes through a read-only mode.
        err = io.StringIO()
        with mock.patch.object(search_indel_snp.overlay_indel_snps, "main",
                               side_effect=PermissionError(13, "Permission denied", "x")), \
                mock.patch("sys.stderr", err):
            rc = search_indel_snp.main(["search_indel_snp", "overlay", "gen", "enr", "log"])
        self.assertEqual(rc, 1)
        self.assertEqual(len(err.getvalue().strip().splitlines()), 1, err.getvalue())
        self.assertIn("ERROR: the SNP codes cannot be put on the fake-indel genome in gen",
                      err.getvalue())

    def test_the_opt_out_is_a_no_op(self):
        self.fake.write_text(PLAIN)
        out = tool("overlay", self.root / "indels", self.enriched, self.logs,
                   env={"CRISPRME_INDEL_SNP": "0"})
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertEqual(self.fake.read_text(), PLAIN)


class TheStore(_Tmp):
    def setUp(self):
        super().setUp()
        self.vcfs = self.root / "vcfs"
        self.vcfs.mkdir()
        self.store = self.root / "indel_genotypes_ds"
        self.scratch = self.root / "scratch"

    def vcf(self, name, chrom):
        with gzip.open(self.vcfs / name, "wt") as fh:
            fh.write(VCF_HEAD + f"{chrom}\t26\trs1\tA\tAT\t.\tPASS\tAF=0.25\tGT\t0|0\t1|0\n")

    def write_store(self, threads=2):
        return tool("store", self.vcfs, self.store, self.scratch, threads)

    def test_each_chromosome_gets_its_store(self):
        self.vcf("x.chr1.vcf.gz", "chr1")
        self.vcf("x.chr2.vcf.gz", "2")  # a VCF without the chr prefix, as 1000G's
        out = self.write_store()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        for chrom in ("chr1", "chr2"):
            with gzip.open(self.store / f"gt_indel_{chrom}.tsv.gz", "rt") as fh:
                self.assertIn("S2", fh.read())
            self.assertTrue((self.store / f"gt_indel_{chrom}.tsv.gz.done").exists())
        self.assertEqual(sorted(p.name for p in self.scratch.iterdir()), [])

    def test_one_thread_writes_the_same(self):
        self.vcf("x.chr1.vcf.gz", "chr1")
        out = self.write_store(threads=1)
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertTrue((self.store / "gt_indel_chr1.tsv.gz").exists())

    def test_a_complete_store_is_kept_and_a_truncated_one_written_again(self):
        self.vcf("x.chr1.vcf.gz", "chr1")
        self.vcf("x.chr2.vcf.gz", "chr2")
        self.store.mkdir()
        (self.store / "gt_indel_chr1.tsv.gz").write_bytes(gzip.compress(b"kept\n"))
        # (a store of 2.10.0: its marker says it holds the records of PASS and '.' alone)
        (self.store / "gt_indel_chr1.tsv.gz.done").write_text(
            f"1\n{search_indel_snp.build_indel_genotypes.FILTER_MARK}\n")
        (self.store / "gt_indel_chr2.tsv.gz").write_bytes(gzip.compress(b"cut\n" * 99)[:12])
        out = self.write_store()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        with gzip.open(self.store / "gt_indel_chr1.tsv.gz", "rt") as fh:
            self.assertEqual(fh.read(), "kept\n")
        with gzip.open(self.store / "gt_indel_chr2.tsv.gz", "rt") as fh:
            self.assertIn("S2", fh.read())

    def test_a_store_written_before_2_10_0_is_written_again(self):
        """Its marker holds the record count alone: it may hold the indels of records whose
        FILTER the search skips (neither PASS nor '.'), which no fake-indel window holds, so
        a carrier would be taken off a row for one. Written again, it holds none of them."""
        with gzip.open(self.vcfs / "x.chr1.vcf.gz", "wt") as fh:
            fh.write(VCF_HEAD + "chr1\t26\trs1\tA\tAT\t.\tPASS\tAF=0.25\tGT\t0|0\t1|0\n"
                     "chr1\t40\trs2\tCA\tC\t.\tLowQual\tAF=0.25\tGT\t1|0\t0|0\n")
        self.store.mkdir()
        old = gzip.compress(b"26_A_AT\tS2:1|0\n40_CA_C\tS1:1|0\n")
        (self.store / "gt_indel_chr1.tsv.gz").write_bytes(old)
        (self.store / "gt_indel_chr1.tsv.gz.done").write_text("2\n")
        out = self.write_store()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        with gzip.open(self.store / "gt_indel_chr1.tsv.gz", "rt") as fh:
            self.assertEqual(fh.read(), "26_A_AT\tS2:1|0\n")
        self.assertEqual((self.store / "gt_indel_chr1.tsv.gz.done").read_text(),
                         f"1\n{search_indel_snp.build_indel_genotypes.FILTER_MARK}\n")

    def test_a_vcf_with_no_record_writes_nothing(self):
        with gzip.open(self.vcfs / "x.chr1.vcf.gz", "wt") as fh:
            fh.write(VCF_HEAD)
        out = self.write_store()
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertFalse(self.store.exists())


class TheCommandLine(_Tmp):
    def test_a_wrong_call_says_how_on_stderr(self):
        out = tool("nonsense")
        self.assertEqual(out.returncode, 2)
        self.assertIn("search_indel_snp.py overlay", out.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
