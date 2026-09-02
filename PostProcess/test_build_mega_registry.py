#!/usr/bin/env python3
"""Tests for build_mega_registry -- the sites-only merged-VCF -> registry driver.

Covers the INFO AF parser (per-dataset AF_<ds>, AF_max skipped, multiallelic-comma
guard, non-positive AF drop) and build() end to end (synthetic gzip VCF -> registry
-> RegistryReader round-trip). STDLIB only.
"""
import gzip
import os
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import build_mega_registry as bmr  # noqa: E402
from tier0_registry import (  # noqa: E402
    RegistryReader, GLOBAL_GROUP_ID, db_group_id,
)

_ALL = set(bmr.DEFAULT_DATASET_META)


class TestParseInfoAfs(unittest.TestCase):
    def test_reads_per_dataset_skips_afmax(self):
        info = "AF_HGDP=0.0305;AF_TOPMed=0.0956;AF_max=0.0956"
        self.assertEqual(bmr._parse_info_afs(info, _ALL),
                         {"HGDP": 0.0305, "TOPMed": 0.0956})

    def test_drops_nonpositive_and_unknown(self):
        info = "AF_HGDP=0.0;AF_gnomAD=0.01;AF_bogus=0.5;DP=10"
        self.assertEqual(bmr._parse_info_afs(info, _ALL), {"gnomAD": 0.01})

    def test_empty_when_no_af_fields(self):
        self.assertEqual(bmr._parse_info_afs("DP=10;NS=3", _ALL), {})

    def test_rejects_nonfinite_and_out_of_range(self):
        # inf/1e400 must be rejected (else round(inf*AN) OverflowError aborts the
        # whole chromosome build); nan and AF>1 rejected too.
        self.assertEqual(bmr._parse_info_afs("AF_gnomAD=inf", _ALL), {})
        self.assertEqual(bmr._parse_info_afs("AF_gnomAD=1e400", _ALL), {})
        self.assertEqual(bmr._parse_info_afs("AF_gnomAD=nan", _ALL), {})
        self.assertEqual(bmr._parse_info_afs("AF_gnomAD=1.5", _ALL), {})
        self.assertEqual(bmr._parse_info_afs("AF_gnomAD=1.0;AF_HGDP=0.3", _ALL),
                         {"gnomAD": 1.0, "HGDP": 0.3})  # exactly 1.0 is allowed


def _write_vcf(path, rows):
    """rows: list of (chrom, pos, id, ref, alt, info)."""
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, i, r, a, info in rows:
            fh.write("\t".join([c, str(p), i, r, a, ".", "PASS", info]) + "\n")


class TestBuildEndToEnd(unittest.TestCase):
    def test_build_and_read_back(self):
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr22.afmax.vcf.gz")
        _write_vcf(vcf, [
            ("chr22", 10516173, "rs1", "A", "G",
             "AF_HGDP=0.030541;AF_TOPMed=0.0956251;AF_max=0.0956251"),
            ("chr22", 500, "rs3", "C", "T",
             "AF_1000G2021=0.10;AF_gnomAD=0.02;AF_max=0.10"),
            ("chr22", 500, "rs4", "C", "A",  # multiallelic split -> distinct alt
             "AF_AoU=0.004;AF_max=0.004"),
            ("chr22", 999, ".", "G", "A",     # no rsID; comma-AF -> float() drops
             "AF_HGDP=0.01,0.02"),
            ("chr22", 1200, ".", "G", "A,T",  # comma-ALT guard -> whole row skipped
             "AF_HGDP=0.5"),
            ("chr22", 2000, "rs9", "T", "TA",  # insertion -> counted as indel, not
             "AF_gnomAD=0.03;AF_max=0.03"),    #   in the SNP registry
            ("chr22", 3000, "rs10", "c", "t",  # lowercase SNP -> uppercased, counted
             "AF_gnomAD=0.02;AF_max=0.02"),    #   as a SNP (not misfiled as indel)
        ])

        manifest, binp, idxp, stats = bmr.build(vcf, "chr22", d)
        self.assertEqual(manifest["aggregation"], "info_af")
        # full accounting: 4 SNPs (3 upper + 1 lowercase), 1 indel, 1 multiallelic
        # (the A,T ALT), 1 no-AF (the comma-AF at 999). Nothing silently dropped.
        self.assertEqual(stats["snps"], 4)
        self.assertEqual(stats["indels"], 1)
        self.assertEqual(stats["multiallelic"], 1)
        self.assertEqual(stats["no_af"], 1)
        reader = RegistryReader(binp, idxp)

        # the lowercase SNP (c>t) is in the registry, uppercased to (C>T)
        self.assertIn(db_group_id("gnomAD"), reader.lookup(3000, "T"))

        g1 = reader.lookup(10516173, "G")
        self.assertEqual(set(g1), {db_group_id("HGDP"), db_group_id("TOPMed"),
                                   GLOBAL_GROUP_ID})
        self.assertLessEqual(
            abs(g1[db_group_id("HGDP")].allele_freq() - 0.030541), 1e-3)

        # multiallelic at pos 500: two distinct alts, distinct dataset sets
        gt = reader.lookup(500, "T")
        ga = reader.lookup(500, "A")
        self.assertEqual(set(gt), {db_group_id("1000G2021"), db_group_id("gnomAD"),
                                   GLOBAL_GROUP_ID})
        self.assertEqual(set(ga), {db_group_id("AoU"), GLOBAL_GROUP_ID})

        # pos 999 has a comma-AF (malformed for biallelic) -> the float() drops it,
        # leaving no AF -> record absent.
        self.assertIsNone(reader.lookup(999, "A"))
        # pos 1200 ALT "A,T" -> comma guard skips the whole record.
        self.assertIsNone(reader.lookup(1200, "A"))
        # pos 2000 T>TA is an indel -> not in the SNP registry.
        self.assertIsNone(reader.lookup(2000, "TA"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
