#!/usr/bin/env python3
"""Tests for build_mega_indel_af -- the sites-only INDEL allele-frequency sidecar.

Covers: only INDELs/MNVs are stored (SNPs excluded, complementing the SNP
registry's SNP-only gate); per-dataset AF + AF_max round-trip through the reader;
comma-ALT + no-AF rows skipped; atomic write leaves no .tmp; and the records read as
CRISPRitz 2.10.0's enricher reads them (FILTER PASS or '.', the alleles it windows, one
REF>ALT on two records the more frequent one's, none it gives no window). STDLIB only.
"""
import gzip
import os
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import build_mega_indel_af as bi  # noqa: E402


def _write_vcf(path, rows):
    """rows: (chrom, pos, id, ref, alt, info[, FILTER]) (FILTER PASS by default)."""
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            c, p, i, r, a, info = row[:6]
            filt = row[6] if len(row) > 6 else "PASS"
            fh.write("\t".join([c, str(p), i, r, a, ".", filt, info]) + "\n")


class TestIndelClassify(unittest.TestCase):
    def test_is_indel(self):
        self.assertFalse(bi._is_indel("A", "G"))   # SNP
        self.assertFalse(bi._is_indel("a", "g"))   # lowercase SNP
        self.assertTrue(bi._is_indel("T", "TA"))   # insertion
        self.assertTrue(bi._is_indel("TA", "T"))   # deletion
        self.assertTrue(bi._is_indel("AT", "GC"))  # MNV (equal-length, not a SNP)
        # a spanning deletion's '*' over a one-base REF: CRISPRitz's enricher windows
        # nothing for it (and codes nothing), so no row can name it
        self.assertFalse(bi._is_indel("A", "*"))
        self.assertFalse(bi._is_indel("A", "<DEL>"))  # symbolic


class TestBuildAndRead(unittest.TestCase):
    def test_store_roundtrip(self):
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr22.afmax.vcf.gz")
        _write_vcf(vcf, [
            ("chr22", 100, "rs1", "A", "G", "AF_HGDP=0.03;AF_max=0.03"),        # SNP -> excluded
            ("chr22", 200, "rs2", "T", "TA", "AF_gnomAD=0.02;AF_TOPMed=0.05;AF_max=0.05"),  # ins
            ("chr22", 300, "rs3", "GAT", "G", "AF_1000G2021=0.1;AF_max=0.1"),   # del
            ("chr22", 400, ".", "C", "T,TA", "AF_HGDP=0.4"),                    # comma-ALT -> skip
            ("chr22", 500, ".", "AT", "GC", "DP=10"),                          # MNV, no AF -> skip
        ])
        n, out = bi.build_indel_af_store(vcf, "chr22", d)
        self.assertEqual(n, 2)  # only the insertion + deletion
        self.assertFalse(os.path.exists(out + ".tmp"))  # atomic: no leftover tmp

        r = bi.IndelAfReader(out)
        self.assertEqual(len(r), 2)
        ins = r.lookup(200, "T", "TA")
        self.assertAlmostEqual(ins["gnomAD"], 0.02)
        self.assertAlmostEqual(ins["TOPMed"], 0.05)
        self.assertAlmostEqual(ins["AF_max"], 0.05)
        self.assertNotIn("HGDP", ins)  # absent dataset not stored
        dele = r.lookup(300, "gat", "g")  # case-insensitive
        self.assertAlmostEqual(dele["1000G2021"], 0.1)
        self.assertIsNone(r.lookup(100, "A", "G"))   # the SNP is not in the indel store
        self.assertIsNone(r.lookup(999, "A", "AA"))  # absent


class TestParseAfMax(unittest.TestCase):
    def test_guards(self):
        self.assertAlmostEqual(bi._parse_af_max("AF_HGDP=0.1;AF_max=0.3"), 0.3)
        self.assertIsNone(bi._parse_af_max("AF_HGDP=0.1"))          # absent
        self.assertIsNone(bi._parse_af_max("AF_max=."))            # missing
        self.assertIsNone(bi._parse_af_max("AF_max=inf"))         # non-finite
        self.assertIsNone(bi._parse_af_max("AF_max=1e400"))       # overflow
        self.assertIsNone(bi._parse_af_max("AF_max=nan"))         # NaN
        self.assertIsNone(bi._parse_af_max("AF_max=1.5"))         # out of range
        self.assertIsNone(bi._parse_af_max("AF_max=0.0"))         # non-positive
        self.assertAlmostEqual(bi._parse_af_max("AF_max=1.0"), 1.0)  # exactly 1 ok


class TestDedupAndFilters(unittest.TestCase):
    def test_duplicate_keys_merged_into_the_more_frequent(self):
        # two identical (pos,ref,alt) indels (as bcftools norm can emit in VNTRs) with
        # different per-dataset AFs -> one store row, the more frequent record's (its
        # INFO/AF, which the build input sets to AF_max), as CRISPRitz's enricher merges
        # one allele's entries; up to CRISPRme 2.8.1 the store took each source's
        # highest, a set of frequencies no record has
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr1.afmax.vcf.gz")
        _write_vcf(vcf, [
            ("chr1", 100, ".", "T", "TA", "AF_gnomAD=0.02;AF_HGDP=0.01;AF_max=0.02;AF=0.02"),
            ("chr1", 100, ".", "T", "TA", "AF_gnomAD=0.05;AF_TOPMed=0.03;AF_max=0.05;AF=0.05"),
            ("chr1", 200, ".", "G", "GC", "AF_HGDP=0.04;AF_max=0.04"),   # no INFO/AF:
            ("chr1", 200, ".", "G", "GC", "AF_TOPMed=0.4;AF_max=0.4"),   # a tie, the earlier
        ])
        n, out = bi.build_indel_af_store(vcf, "chr1", d)
        self.assertEqual(n, 2)  # deduped to one unique key each
        r = bi.IndelAfReader(out)
        self.assertEqual(len(r), 2)
        self.assertEqual(r.lookup(100, "T", "TA"),
                         {"gnomAD": 0.05, "TOPMed": 0.03, "AF_max": 0.05})
        self.assertEqual(r.lookup(200, "G", "GC"), {"HGDP": 0.04, "AF_max": 0.04})

    def test_afmax_absent_roundtrip(self):
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr1.afmax.vcf.gz")
        _write_vcf(vcf, [("chr1", 200, ".", "C", "CT", "AF_HGDP=0.02")])  # no AF_max
        _, out = bi.build_indel_af_store(vcf, "chr1", d)
        rec = bi.IndelAfReader(out).lookup(200, "C", "CT")
        self.assertAlmostEqual(rec["HGDP"], 0.02)
        self.assertNotIn("AF_max", rec)  # absent AF_max stays absent (not fabricated)

    def test_chrom_filter_excludes_other_contigs(self):
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "multi.vcf.gz")
        _write_vcf(vcf, [
            ("chr1", 300, ".", "A", "AT", "AF_HGDP=0.02;AF_max=0.02"),
            ("chr2", 300, ".", "A", "AT", "AF_HGDP=0.09;AF_max=0.09"),  # same POS, other chrom
        ])
        _, out = bi.build_indel_af_store(vcf, "chr1", d)
        r = bi.IndelAfReader(out)
        self.assertEqual(len(r), 1)                       # chr2 excluded
        self.assertAlmostEqual(r.lookup(300, "A", "AT")["HGDP"], 0.02)  # chr1's value, not chr2's


class TestPartition(unittest.TestCase):
    def test_registry_and_indel_store_partition_the_alleles_the_enricher_reads(self):
        # every allele CRISPRitz's enricher reads lands in EXACTLY ONE store (SNP
        # registry XOR indel sidecar) -- no double-count, no gap -- and one it does not
        # read (a record outside PASS and '.', an ALT it does not code) in neither. REF
        # and ALT are read in upper case, as CRISPRitz 2.10.0 reads them: an SNV written
        # a>t is the registry's A>T, an indel written ga>g the sidecar's GA>G
        import build_mega_registry as bmr
        from tier0_registry import RegistryReader
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr1.afmax.vcf.gz")
        rows = [
            ("chr1", 100, ".", "A", "G", "AF_HGDP=0.03;AF_max=0.03"),   # SNP
            ("chr1", 200, ".", "a", "t", "AF_gnomAD=0.02;AF_max=0.02"), # lower-case SNP
            ("chr1", 300, ".", "T", "TA", "AF_HGDP=0.01;AF_max=0.01"),  # ins
            ("chr1", 400, ".", "ga", "g", "AF_AoU=0.02;AF_max=0.02"),   # del, lower case
            ("chr1", 500, ".", "AT", "GC", "AF_HGDP=0.02;AF_max=0.02"), # MNV -> indel side
            ("chr1", 600, ".", "C", "G", "AF_HGDP=0.02;AF_max=0.02", "LOW_VQSLOD"),
            ("chr1", 700, ".", "CA", "C", "AF_HGDP=0.02;AF_max=0.02", "ExcHet"),
            ("chr1", 800, ".", "G", "n", "AF_HGDP=0.02;AF_max=0.02"),   # not coded
        ]
        _write_vcf(vcf, rows)
        _, rbin, ridx, stats = bmr.build(vcf, "chr1", os.path.join(d, "reg"))
        n_ind, iout = bi.build_indel_af_store(vcf, "chr1", os.path.join(d, "ind"))
        reg = RegistryReader(rbin, ridx)
        ind = bi.IndelAfReader(iout)
        self.assertEqual(stats["written"], 2)   # A>G, a>t
        self.assertEqual(n_ind, 3)              # ins + del + MNV
        self.assertEqual(stats["not_coded"], 1)  # the n
        self.assertEqual(stats["written"] + n_ind + stats["not_coded"] + stats["filtered"],
                         len(rows))
        # each variant in exactly one store
        self.assertIsNotNone(reg.lookup(100, "G"));  self.assertIsNone(ind.lookup(100, "A", "G"))
        self.assertIsNotNone(reg.lookup(200, "T"));  self.assertIsNone(ind.lookup(200, "A", "T"))
        self.assertIsNotNone(ind.lookup(300, "T", "TA")); self.assertIsNone(reg.lookup(300, "TA"))
        self.assertIsNotNone(ind.lookup(400, "GA", "G")); self.assertIsNone(reg.lookup(400, "G"))
        self.assertIsNotNone(ind.lookup(500, "AT", "GC")); self.assertIsNone(reg.lookup(500, "GC"))
        for pos, ref, alt in ((600, "C", "G"), (700, "CA", "C"), (800, "G", "N")):
            self.assertIsNone(reg.lookup(pos, alt))
            self.assertIsNone(ind.lookup(pos, ref, alt))


class TestReadAsTheEnricher(unittest.TestCase):
    """The indel frequencies of the alleles CRISPRitz 2.10.0's enricher windows, read as
    it reads the VCF."""

    def build(self, rows, chrom="chr22"):
        d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, d, True)
        vcf = os.path.join(d, "mega.%s.build.vcf.gz" % chrom)
        _write_vcf(vcf, rows)
        stats = {}
        _n, out = bi.build_indel_af_store(vcf, chrom, d, stats=stats)
        return bi.IndelAfReader(out), stats

    def test_records_outside_pass_are_left_out(self):
        r, stats = self.build([
            ("chr22", 100, ".", "TA", "T", "AF_HGDP=0.1;AF_max=0.1", "LOW_VQSLOD"),
            ("chr22", 200, ".", "C", "CT", "AF_HGDP=0.2;AF_max=0.2", "."),
        ])
        self.assertIsNone(r.lookup(100, "TA", "T"))
        self.assertIsNotNone(r.lookup(200, "C", "CT"))
        self.assertEqual(stats["filtered"], 1)

    def test_every_indel_of_a_line_of_several_alts(self):
        r, _stats = self.build([
            ("chr22", 100, ".", "CA", "C,CAA,*", "AF_HGDP=0.1,0.2,0.3;AF_max=0.1,0.2,0.3"),
        ])
        self.assertEqual(r.lookup(100, "CA", "C"), {"HGDP": 0.1, "AF_max": 0.1})
        self.assertEqual(r.lookup(100, "CA", "CAA"), {"HGDP": 0.2, "AF_max": 0.2})
        # a '*' over a two-base REF is windowed by the enricher (another length): kept
        self.assertEqual(r.lookup(100, "CA", "*"), {"HGDP": 0.3, "AF_max": 0.3})

    def test_no_window_within_26_bases_of_the_start(self):
        # indel_to_fasta writes no window for an indel at POS < 26 (its window would start
        # before the chromosome), so no row names it
        r, stats = self.build([
            ("chr22", 25, ".", "TA", "T", "AF_HGDP=0.1;AF_max=0.1"),
            ("chr22", 26, ".", "GA", "G", "AF_HGDP=0.1;AF_max=0.1"),
        ])
        self.assertIsNone(r.lookup(25, "TA", "T"))
        self.assertIsNotNone(r.lookup(26, "GA", "G"))
        self.assertEqual(stats["no_window"], 1)


    def test_the_window_is_the_one_the_indel_log_shows(self):
        """The left flank of the enricher's windows is read from the indel log the build
        wrote (--log-dir): with CRISPRitz 2.9.1's 25/27 window an indel at POS 26 has a
        window, with a wider one (2.10.0's 30/30, lane W's, makes the flank 31 here) it has
        none, and no row can name it. Up to review r1 of lane M the 26 was fixed in the code, which
        a wider window would have left wrong at the start of every chromosome."""
        d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, d, True)
        rows = [("chr22", 26, ".", "GA", "G", "AF_HGDP=0.1;AF_max=0.1"),
                ("chr22", 31, ".", "TC", "T", "AF_HGDP=0.1;AF_max=0.1")]
        vcf = os.path.join(d, "mega.chr22.build.vcf.gz")
        _write_vcf(vcf, rows)
        logs = os.path.join(d, "log_indels_hg38_mega")
        os.makedirs(logs)
        with gzip.open(os.path.join(logs, "logchr22.txt.gz"), "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
                     "chr22_0-62_1\t\t.\t0.1\tchr22_31_TC_T\t0,61\tN\n")
        self.assertEqual(bi.window_left_of_log(logs, "chr22"), 31)
        self.assertIsNone(bi.window_left_of_log(logs, "chr21"))
        out = os.path.join(d, "store")
        import subprocess
        done = subprocess.run([sys.executable, os.path.join(_HERE, "build_mega_indel_af.py"),
                               "--vcf", vcf, "--chrom", "chr22", "--out-dir", out,
                               "--log-dir", logs], capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("1 within 31 bases of the start", done.stderr)
        r = bi.IndelAfReader(os.path.join(out, "indel_af_chr22.tsv.gz"))
        self.assertIsNone(r.lookup(26, "GA", "G"))
        self.assertIsNotNone(r.lookup(31, "TC", "T"))
        # without a log, the window of CRISPRitz 2.10.0's enricher
        done = subprocess.run([sys.executable, os.path.join(_HERE, "build_mega_indel_af.py"),
                               "--vcf", vcf, "--chrom", "chr22", "--out-dir", out],
                              capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIsNotNone(bi.IndelAfReader(os.path.join(out, "indel_af_chr22.tsv.gz"))
                             .lookup(26, "GA", "G"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
