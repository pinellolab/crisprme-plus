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


if __name__ == "__main__":
    unittest.main(verbosity=2)
