#!/usr/bin/env python3
"""Tests for build_mega_registry -- the sites-only merged-VCF -> registry driver.

Covers the INFO AF parser (per-dataset AF_<ds>, AF_max skipped, a comma list read per
ALT, non-positive AF drop) and build() end to end (synthetic gzip VCF -> registry
-> RegistryReader round-trip; a record listing several ALTs read one ALT at a time), and
the records it reads, as CRISPRitz 2.10.0's enricher reads them: FILTER PASS or '.', the
SNP ALTs it codes, one REF>ALT on two records merged into the more frequent. STDLIB only.
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
    """rows: list of (chrom, pos, id, ref, alt, info[, FILTER]) (FILTER PASS by default)."""
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            c, p, i, r, a, info = row[:6]
            filt = row[6] if len(row) > 6 else "PASS"
            fh.write("\t".join([c, str(p), i, r, a, ".", filt, info]) + "\n")


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
            ("chr22", 1200, ".", "G", "A,T",  # two ALTs, one AF -> neither has its own
             "AF_HGDP=0.5"),
            ("chr22", 1300, "rs7;rs8", "G", "A,T",  # two ALTs: each its AF and its ID
             "AF_HGDP=0.5,0.25;AF_gnomAD=.,0.01;AF_max=0.5,0.25"),
            ("chr22", 1400, "rs11", "G", "C,GT",  # a SNP and an indel on one line
             "AF_HGDP=0.2,0.3"),
            ("chr22", 2000, "rs9", "T", "TA",  # insertion -> counted as indel, not
             "AF_gnomAD=0.03;AF_max=0.03"),    #   in the SNP registry
            ("chr22", 3000, "rs10", "c", "t",  # a lower-case ALT: CRISPRitz codes none,
             "AF_gnomAD=0.02;AF_max=0.02"),    #   so the registry holds none
        ])

        manifest, binp, idxp, stats = bmr.build(vcf, "chr22", d)
        self.assertEqual(manifest["aggregation"], "info_af")
        # full accounting, per ALT: 9 SNP alleles read (3 upper, the comma-AF at 999, the
        # two of 1200 and of 1300, the C of 1400), 6 of them written (3 no-AF: the
        # comma-AF at 999, both ALTs of 1200), 2 indels (2000, the GT of 1400), 1 ALT
        # CRISPRitz does not code (the lower-case t); 3 records listed several ALTs.
        # Nothing silently dropped.
        self.assertEqual(stats["snps"], 9)
        self.assertEqual(stats["written"], 6)
        self.assertEqual(stats["indels"], 2)
        self.assertEqual(stats["not_coded"], 1)
        self.assertEqual(stats["multiallelic"], 3)
        self.assertEqual(stats["no_af"], 3)
        reader = RegistryReader(binp, idxp)

        # a record listing several ALTs: each ALT its own AF (the k-th of each list)
        # and its own ID (one per ALT in the ID column); it was skipped whole
        self.assertEqual(reader.alts_at(1300), ["A", "T"])
        ga, gt = reader.lookup(1300, "A"), reader.lookup(1300, "T")
        self.assertLessEqual(abs(ga[db_group_id("HGDP")].allele_freq() - 0.5), 1e-3)
        self.assertNotIn(db_group_id("gnomAD"), ga)
        self.assertLessEqual(abs(gt[db_group_id("HGDP")].allele_freq() - 0.25), 1e-3)
        self.assertLessEqual(abs(gt[db_group_id("gnomAD")].allele_freq() - 0.01), 1e-3)
        self.assertEqual((reader.rsid(1300, "A"), reader.rsid(1300, "T")), ("rs7", "rs8"))
        self.assertEqual(reader.alts_at(1400), ["C"])
        self.assertLessEqual(
            abs(reader.lookup(1400, "C")[db_group_id("HGDP")].allele_freq() - 0.2), 1e-3)
        self.assertEqual(reader.rsid(1400, "C"), "rs11")

        # the lower-case SNP (c>t) is not: CRISPRitz's enricher codes no lower-case ALT
        # (up to CRISPRme 2.8.1 it was upper-cased into the registry, a C>T the genome
        # never coded)
        self.assertIsNone(reader.lookup(3000, "T"))

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
        # pos 1200 ALT "A,T" with one AF value: no ALT has an AF of its own.
        self.assertIsNone(reader.lookup(1200, "A"))
        self.assertIsNone(reader.lookup(1200, "T"))
        # pos 2000 T>TA is an indel -> not in the SNP registry.
        self.assertIsNone(reader.lookup(2000, "TA"))


class TestChromFilterAndGlobal(unittest.TestCase):
    def test_chrom_filter_excludes_other_contigs(self):
        # a multi-contig VCF must not leak other chroms into the registry (the
        # registry is keyed by (pos,alt) with NO chrom -> same-POS would collide).
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "multi.vcf.gz")
        _write_vcf(vcf, [
            ("chr1", 500, "rsA", "C", "T", "AF_HGDP=0.10;AF_max=0.10"),
            ("chr2", 500, "rsB", "C", "T", "AF_HGDP=0.90;AF_max=0.90"),  # same POS, other chrom
        ])
        manifest, binp, idxp, stats = bmr.build(vcf, "chr1", d)
        self.assertEqual(stats["snps"], 1)  # chr2 excluded
        reg = RegistryReader(binp, idxp)
        g = reg.lookup(500, "T")
        self.assertLessEqual(abs(g[db_group_id("HGDP")].allele_freq() - 0.10), 1e-3)  # chr1's, not 0.90

    def test_global_matches_vcf_afmax_token(self):
        # the GLOBAL group's AF must equal the VCF's ACTUAL INFO/AF_max token, not
        # just a value re-derived inside the test.
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "mega.chr1.afmax.vcf.gz")
        _write_vcf(vcf, [("chr1", 700, "rsC", "C", "T",
                          "AF_1000G2021=0.0377889;AF_HGDP=0.030541;AF_max=0.0377889")])
        _, binp, idxp, _ = bmr.build(vcf, "chr1", d)
        gg = RegistryReader(binp, idxp).lookup(700, "T")[GLOBAL_GROUP_ID]
        self.assertLessEqual(abs(gg.allele_freq() - 0.0377889), 1e-4)  # == the AF_max token


class TestReadAsTheEnricher(unittest.TestCase):
    """The registry holds the SNP alleles CRISPRitz 2.10.0's enricher codes, read as it
    reads the VCF: the hg38_mega VCFs carry HGDP-2019's FILTER (LOW_VQSLOD, ExcHet) onto
    merged records, which the enricher skips, while the registry of CRISPRme 2.8.1 held
    them and its SNP rows named them (2,442,989 SNVs genome-wide)."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, self.d, True)
        self.vcf = os.path.join(self.d, "mega.chr22.build.vcf.gz")

    def build(self, rows):
        _write_vcf(self.vcf, rows)
        manifest, binp, idxp, stats = bmr.build(self.vcf, "chr22", self.d)
        return manifest, RegistryReader(binp, idxp), stats

    def test_records_outside_pass_are_left_out(self):
        _m, reg, stats = self.build([
            ("chr22", 100, "rs1", "A", "G", "AF_HGDP=0.1;AF_max=0.1", "LOW_VQSLOD"),
            ("chr22", 200, "rs2", "C", "T", "AF_HGDP=0.2;AF_max=0.2", "ExcHet;LOW_VQSLOD"),
            ("chr22", 300, "rs3", "G", "A", "AF_HGDP=0.3;AF_max=0.3", "."),
            ("chr22", 400, "rs4", "T", "C", "AF_HGDP=0.4;AF_max=0.4", "PASS"),
            ("chr22", 400, "rs5", "T", "G", "AF_gnomAD=0.01;AF_max=0.01", "pass"),
        ])
        self.assertIsNone(reg.lookup(100, "G"))
        self.assertIsNone(reg.lookup(200, "T"))
        self.assertIsNone(reg.lookup(400, "G"))  # "pass" is not PASS
        self.assertIsNotNone(reg.lookup(300, "A"))
        self.assertIsNotNone(reg.lookup(400, "C"))
        self.assertEqual(reg.alts_at(400), ["C"])
        self.assertEqual((stats["filtered"], stats["written"]), (3, 2))

    def test_alts_the_enricher_does_not_code(self):
        # the '*' of a spanning deletion, an N, a lower-case base: no code, no record;
        # every ALT it codes of the same line is read
        _m, reg, stats = self.build([
            ("chr22", 100, ".", "G", "A,*", "AF_HGDP=0.1,0.05;AF_max=0.1,0.05"),
            ("chr22", 200, ".", "C", "N", "AF_HGDP=0.2;AF_max=0.2"),
            ("chr22", 300, ".", "A", "g", "AF_HGDP=0.3;AF_max=0.3"),
        ])
        self.assertEqual(reg.alts_at(100), ["A"])
        self.assertEqual(reg.alts_at(200), [])
        self.assertEqual(reg.alts_at(300), [])
        self.assertEqual((stats["written"], stats["not_coded"]), (1, 3))

    def test_one_allele_on_two_records_is_the_more_frequent(self):
        # the enricher merges one REF>ALT of two records into one entry, which takes the
        # AF and ID of the more frequent record (the higher INFO/AF, found by its key;
        # on a tie, the earlier)
        _m, reg, stats = self.build([
            ("chr22", 100, "rsA", "A", "G", "AF_HGDP=0.01;AF_max=0.01;AF=0.01"),
            ("chr22", 100, "rsB", "A", "G", "AF_gnomAD=0.03;AF_max=0.03;AF=0.03"),
            ("chr22", 200, "rsC", "C", "T", "AF_HGDP=0.2;AF_max=0.2;AF=0.2"),
            ("chr22", 200, "rsD", "C", "T", "AF_TOPMed=0.1;AF_max=0.1;AF=0.1"),
            ("chr22", 300, "rsE", "G", "C", "AF_HGDP=0.05;AF_max=0.05"),  # no INFO/AF:
            ("chr22", 300, "rsF", "G", "C", "AF_AoU=0.5;AF_max=0.5"),     # a tie
        ])
        self.assertEqual(len(reg), 3)
        g = reg.lookup(100, "G")
        self.assertEqual(set(g), {db_group_id("gnomAD"), GLOBAL_GROUP_ID})
        self.assertEqual(reg.rsid(100, "G"), "rsB")
        g = reg.lookup(200, "T")
        self.assertEqual(set(g), {db_group_id("HGDP"), GLOBAL_GROUP_ID})
        self.assertEqual(reg.rsid(200, "T"), "rsC")
        self.assertEqual(reg.rsid(300, "C"), "rsE")
        self.assertEqual(stats["merged"], 3)

    def test_the_manifest_says_which_records(self):
        manifest, _reg, _stats = self.build([("chr22", 100, ".", "A", "G", "AF_HGDP=0.1")])
        self.assertIn("PASS", manifest["vcf_records"])
        self.assertIn("CRISPRitz 2.10.0", manifest["vcf_records"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
