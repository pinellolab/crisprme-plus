#!/usr/bin/env python3
"""Tests for build_mega_indel_af -- the sites-only INDEL allele-frequency sidecar.

Covers: only INDELs/MNVs are stored (SNPs excluded, complementing the SNP
registry's SNP-only gate); per-dataset AF + AF_max round-trip through the reader;
comma-ALT + no-AF rows skipped; atomic write leaves no .tmp. STDLIB only.
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
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, i, r, a, info in rows:
            fh.write("\t".join([c, str(p), i, r, a, ".", "PASS", info]) + "\n")


class TestIndelClassify(unittest.TestCase):
    def test_is_indel(self):
        self.assertFalse(bi._is_indel("A", "G"))   # SNP
        self.assertFalse(bi._is_indel("a", "g"))   # lowercase SNP
        self.assertTrue(bi._is_indel("T", "TA"))   # insertion
        self.assertTrue(bi._is_indel("TA", "T"))   # deletion
        self.assertTrue(bi._is_indel("AT", "GC"))  # MNV (equal-length, not a SNP)
        self.assertTrue(bi._is_indel("A", "*"))    # spanning deletion


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
    def test_duplicate_keys_merged_by_max(self):
        # two identical (pos,ref,alt) indels (as bcftools norm can emit in VNTRs)
        # with DIFFERENT per-dataset AFs -> one store row, max per dataset, no drop.
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr1.afmax.vcf.gz")
        _write_vcf(vcf, [
            ("chr1", 100, ".", "T", "TA", "AF_gnomAD=0.02;AF_HGDP=0.01;AF_max=0.02"),
            ("chr1", 100, ".", "T", "TA", "AF_gnomAD=0.05;AF_TOPMed=0.03;AF_max=0.05"),
        ])
        n, out = bi.build_indel_af_store(vcf, "chr1", d)
        self.assertEqual(n, 1)  # deduped to one unique key
        r = bi.IndelAfReader(out)
        self.assertEqual(len(r), 1)
        rec = r.lookup(100, "T", "TA")
        self.assertAlmostEqual(rec["gnomAD"], 0.05)  # max(0.02, 0.05)
        self.assertAlmostEqual(rec["HGDP"], 0.01)    # only in record 1
        self.assertAlmostEqual(rec["TOPMed"], 0.03)  # only in record 2
        self.assertAlmostEqual(rec["AF_max"], 0.05)

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
    def test_registry_and_indel_store_partition_all_variants(self):
        # every VCF variant lands in EXACTLY ONE store (SNP registry XOR indel
        # sidecar) -- no double-count, no gap.
        import build_mega_registry as bmr
        from tier0_registry import RegistryReader
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr1.afmax.vcf.gz")
        rows = [
            ("chr1", 100, ".", "A", "G", "AF_HGDP=0.03;AF_max=0.03"),   # SNP
            ("chr1", 200, ".", "a", "t", "AF_gnomAD=0.02;AF_max=0.02"), # lowercase SNP
            ("chr1", 300, ".", "T", "TA", "AF_HGDP=0.01;AF_max=0.01"),  # ins
            ("chr1", 400, ".", "GA", "G", "AF_AoU=0.02;AF_max=0.02"),   # del
            ("chr1", 500, ".", "AT", "GC", "AF_HGDP=0.02;AF_max=0.02"), # MNV -> indel side
        ]
        _write_vcf(vcf, rows)
        _, rbin, ridx, stats = bmr.build(vcf, "chr1", os.path.join(d, "reg"))
        n_ind, iout = bi.build_indel_af_store(vcf, "chr1", os.path.join(d, "ind"))
        reg = RegistryReader(rbin, ridx)
        ind = bi.IndelAfReader(iout)
        self.assertEqual(stats["snps"], 2)   # A>G + a>t
        self.assertEqual(n_ind, 3)           # ins + del + MNV
        self.assertEqual(stats["snps"] + n_ind, len(rows))  # partition covers all
        # each variant in exactly one store
        self.assertIsNotNone(reg.lookup(100, "G"));  self.assertIsNone(ind.lookup(100, "A", "G"))
        self.assertIsNotNone(ind.lookup(300, "T", "TA")); self.assertIsNone(reg.lookup(300, "TA"))
        self.assertIsNotNone(ind.lookup(500, "AT", "GC")); self.assertIsNone(reg.lookup(500, "GC"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
