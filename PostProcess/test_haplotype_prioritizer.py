#!/usr/bin/env python
"""The carrier gate (haplotype_prioritizer.py): the tiers of a row's variants against phased
stores the repository's own writers build (tier1_genotypes, build_indel_genotypes), the
panels a search reads, and the command the search pipeline runs, on an integrated results
table it must leave as it is.

Run: cd PostProcess && python -m unittest test_haplotype_prioritizer -v
"""

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import haplotype_prioritizer as hp  # noqa: E402
import tier1_genotypes as t1  # noqa: E402
import build_indel_genotypes as big  # noqa: E402
import test_build_index_only_atomic as bio  # noqa: E402  (the real crisprme.py, a stub engine)
import test_dictless_publish as tdp  # noqa: E402  (publish-index with the HF client stubbed)

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule

CHROM = "chr22"
PANEL = "hg38_TESTPANEL"
# 20 samples: 40 haplotypes on an autosome
SAMPLES = [f"S{i:02d}" for i in range(20)]
SAMPLE_META = {s: ("TEST", "AFR" if i < 12 else "EUR", "female") for i, s in enumerate(SAMPLES)}


def build_panel(dicts, snps, indels=(), name=PANEL):
    """A phased store of ``name`` under ``dicts``: ``snps`` [(pos, REF, ALT, {sample: gt})]
    through tier1_genotypes, ``indels`` [(pos, REF, ALT, {sample: gt})] through a VCF and
    build_indel_genotypes, as build-index-only writes them."""
    gdir = os.path.join(dicts, f"genotypes_{name}")
    os.makedirs(gdir, exist_ok=True)
    t1.compile_genotypes([(p, a, dict(g), r) for p, r, a, g in snps],
                         t1.build_sample_axis(SAMPLE_META),
                         os.path.join(gdir, f"gt_{CHROM}.bin"),
                         os.path.join(gdir, f"gt_{CHROM}.idx"))
    if indels:
        idir = os.path.join(dicts, f"indel_genotypes_{name}")
        os.makedirs(idir, exist_ok=True)
        vcf = os.path.join(dicts, f"{name}.indels.vcf")
        with open(vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(SAMPLES) + "\n")
            for p, r, a, gts in indels:
                fh.write(f"{CHROM}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t"
                         + "\t".join(gts.get(s, "0|0") for s in SAMPLES) + "\n")
        big.compile_indel_genotypes(vcf, os.path.join(idir, f"gt_indel_{CHROM}.tsv.gz"))
    return gdir


def v(pos, ref, alt):
    return (CHROM, pos, ref, alt)


class TestParseAndSupport(unittest.TestCase):
    def test_parse_variants(self):
        self.assertEqual(hp.parse_variants("NA"), [])
        self.assertEqual(hp.parse_variants("n"), [])
        self.assertEqual(
            hp.parse_variants("chr1_KI270706v1_random_100_A_G,chr2_5_AT_A,chr2_5_AT_A"),
            [("chr1_KI270706v1_random", 100, "A", "G"), ("chr2", 5, "AT", "A")])

    def test_cis_trans_unphased_haploid(self):
        a = {"S1": "1|0", "S2": "1|0", "S3": "1/0", "S4": "1", "S5": "1|1"}
        b = {"S1": "1|0", "S2": "0|1", "S3": "0/1", "S4": "1", "S5": "0|1"}
        haps, cis, unknown, carrying = hp.support([a, b])
        # S1: cis on slot 0; S2: trans; S3: unphased; S4: haploid, both on its one slot;
        # S5: the second slot holds both
        self.assertEqual((haps, cis, unknown, carrying), (3, {"S1", "S4", "S5"}, {"S3"}, 3))

    def test_breaker_on_the_slot_or_the_other(self):
        snp = {"S1": "1|0", "S2": "1|0", "S3": "1|0"}
        indel = {"S1": "1|0", "S2": "0|1", "S3": "1/0"}
        haps, cis, unknown, _c = hp.support([snp], breakers=[indel])
        # S1 loses its slot to the indel; S2's indel is on the other slot; S3's indel call
        # is unphased, so its slot may or may not hold the target
        self.assertEqual((haps, cis, unknown), (1, {"S2"}, {"S3"}))

    def test_an_indel_that_leaves_the_targets_bases(self):
        # GGGG + 10 CA units + TTTT; a target over 6 units inside the repeat
        ref = "GGGG" + "CA" * 10 + "TTTT"

        def seq_at(a, b):
            return ref[max(a, 0):max(b, 0)]
        # a unit lost or gained inside it: the repeat goes on past the target, so the target's
        # bases are there, read from its first base on
        self.assertTrue(hp.holds_target(seq_at, 8, 19, {}, (12, "ACA", "A")))
        self.assertTrue(hp.holds_target(seq_at, 8, 19, {}, (12, "A", "ACA")))
        # with the row's SNP allele kept on the haplotype
        self.assertTrue(hp.holds_target(seq_at, 8, 19, {10: "G"}, (14, "ACA", "A")))
        # over the whole repeat and past its end: one unit fewer is another sequence
        self.assertFalse(hp.holds_target(seq_at, 2, 25, {}, (12, "ACA", "A")))
        # a deletion in a stretch with no repeat
        plain = "ACGTTGCAAGTCCATGGTCA"
        self.assertFalse(hp.holds_target(lambda a, b: plain[max(a, 0):b], 2, 17, {},
                                         (8, "CAA", "C")))

    def test_two_alts_on_one_slot_is_unknown(self):
        alt_g = {"S1": "1|0", "S2": "1|0"}
        alt_t = {"S1": "1|0"}  # S1's first slot holds G and T at one position
        haps, cis, unknown, _c = hp.support([alt_g], rivals=[alt_t])
        self.assertEqual((haps, cis, unknown), (1, {"S2"}, {"S1"}))


class TestClassify(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def panel(self, snps, indels=()):
        build_panel(self.dicts, snps, indels)
        p = hp.PanelChrom(PANEL, self.dicts, CHROM)
        self.addCleanup(p.close)
        return p

    def test_a_two_variants_on_one_haplotype(self):
        p = self.panel([
            (100, "A", "G", {"S00": "1|0", "S01": "1|0", "S02": "0|1"}),
            (105, "C", "T", {"S00": "1|0", "S01": "1|0", "S03": "1|0"}),
        ])
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p], span=(95, 117))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("A", "", 2))
        self.assertEqual(ver.samples, ("S00", "S01"))
        self.assertEqual(ver.panel_haplotypes, 40)
        self.assertEqual(ver.label(), f"together on 2 of 40 {PANEL} haplotypes")

    def test_a_one_haplotype(self):
        p = self.panel([
            (100, "A", "G", {"S00": "1|0", "S02": "0|1"}),
            (105, "C", "T", {"S00": "1|0", "S02": "1|0"}),
        ])
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p])
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("A", hp.ONE_HAPLOTYPE, 1))

    def test_c1_never_together_with_and_without_power(self):
        # common variants, always in trans: 40 x 0.5 x 0.5 = 10 haplotypes expected
        common_a = {s: "1|0" for s in SAMPLES}
        common_b = {s: "0|1" for s in SAMPLES}
        p = self.panel([(100, "A", "G", common_a), (105, "C", "T", common_b),
                        (110, "G", "A", {"S00": "1|0"}), (112, "T", "C", {"S01": "1|0"})])
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p])
        self.assertEqual((ver.tier, ver.state), ("C1", ""))
        self.assertAlmostEqual(ver.expected, 10.0)
        self.assertEqual(ver.label(), f"not together on any of 40 {PANEL} haplotypes "
                                      "(expected 10)")
        # two singletons: 40 x 1/40 x 1/40 = 0.025 expected, a zero says nothing
        rare = hp.classify([v(110, "G", "A"), v(112, "T", "C")], [p])
        self.assertEqual((rare.tier, rare.state), ("C1", hp.TOO_RARE))
        self.assertIn("too rare to tell (expected 0.025)", rare.label())
        # the threshold is a parameter: 0 turns the state off
        off = hp.classify([v(110, "G", "A"), v(112, "T", "C")], [p], min_expected=0)
        self.assertEqual((off.tier, off.state), ("C1", ""))

    def test_c2_and_b(self):
        p = self.panel([(100, "A", "G", {"S00": "1|0"})])
        c2 = hp.classify([v(100, "A", "G"), v(300, "C", "T")], [p],
                         af_of=lambda var: {100: 0.02, 300: 0.0001}[var[1]])
        self.assertEqual((c2.tier, c2.missing, c2.n_variants), ("C2", 1, 2))
        self.assertEqual(c2.label(), f"not checked: 1 of 2 variants not in {PANEL}; frequency "
                                     "product 2.0e-06, rarest 1.0e-04 (below the 1e-06 / "
                                     "0.001 bound)")
        self.assertEqual(hp.classify([v(100, "A", "G")], [p]).tier, "B")
        absent = hp.classify([v(300, "C", "T")], [p])
        self.assertEqual((absent.tier, absent.state), ("B", hp.NOT_IN_PANEL))
        # a REF the panel does not have is not that variant
        self.assertEqual(hp.classify([v(100, "C", "G"), v(300, "C", "T")], [p]).missing, 2)

    def test_indel_in_the_target_takes_the_slot(self):
        # the SNP at 100 on slot 0 of S00 and S01; a phased deletion at 104 on the same
        # slots of both, inside a target over 90..112: no haplotype has the row's sequence
        snps = [(100, "A", "G", {"S00": "1|0", "S01": "1|0"})]
        indels = [(104, "CTT", "C", {"S00": "1|0", "S01": "1|0", "S05": "0|1"})]
        p = self.panel(snps, indels)
        ver = hp.classify([v(100, "A", "G")], [p], span=(90, 112))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("C1", hp.INDEL_IN_TARGET, 2))
        self.assertEqual(ver.label(), f"all 2 {PANEL} haplotypes carrying it have an indel "
                                      "in the target")
        # a target that ends before the deletion's bases keeps them
        self.assertEqual(hp.classify([v(100, "A", "G")], [p], span=(90, 103)).tier, "B")

    def test_a_repeat_unit_lost_in_the_target_keeps_the_haplotype(self):
        # a CA repeat on chr22 bases 100..119 (0-based); the SNP C>G at 104 and a CA unit
        # lost after 109 on the same slots: given the reference, the target over the inner
        # units keeps its bases, so the carriers keep it; without it, positions alone
        ref = "G" * 100 + "CA" * 10 + "T" * 20
        snps = [(105, "C", "G", {"S00": "1|0", "S01": "1|0"})]
        indels = [(110, "ACA", "A", {"S00": "1|0", "S01": "1|0"})]
        self.assertEqual(ref[104] + ref[109:112], "CACA")
        p = self.panel(snps, indels)
        span = (102, 113)
        p.seq_at = lambda a, b: ref[max(a, 0):max(b, 0)]
        self.assertEqual(hp.classify([v(105, "C", "G")], [p], span=span).tier, "B")
        p.seq_at = None
        self.assertEqual(hp.classify([v(105, "C", "G")], [p], span=span).state,
                         hp.INDEL_IN_TARGET)

    def test_two_snps_together_only_beside_an_indel_in_the_target(self):
        snps = [(100, "A", "G", {"S00": "1|0", "S01": "1|0"}),
                (102, "C", "T", {"S00": "1|0", "S01": "1|0"})]
        indels = [(104, "CTT", "C", {"S00": "1|0", "S01": "1|0"})]
        p = self.panel(snps, indels)
        ver = hp.classify([v(100, "A", "G"), v(102, "C", "T")], [p], span=(90, 112))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("C1", hp.INDEL_IN_TARGET, 2))
        self.assertEqual(ver.label(), f"all 2 {PANEL} haplotypes carrying them together have "
                                      "an indel in the target")

    def test_indel_on_the_other_slot_unphased_or_at_the_snp(self):
        snps = [(100, "A", "G", {"S00": "1|0", "S01": "1|0", "S02": "1|0"})]
        indels = [
            (104, "CTT", "C", {"S00": "0|1", "S01": "1/0"}),  # other slot; unphased
            (100, "AT", "A", {"S02": "1|0"}),  # at the SNP's own position: an artefact
        ]
        p = self.panel(snps, indels)
        self.assertEqual(hp.classify([v(100, "A", "G")], [p], span=(90, 112)).tier, "B")
        only_s02 = self.panel([(100, "A", "G", {"S02": "1|0"})], indels)
        self.assertEqual(hp.classify([v(100, "A", "G")], [only_s02], span=(90, 112)).tier, "B")

    def test_an_indel_over_the_snps_takes_the_slot(self):
        """F1: a deletion whose REF holds the row's SNPs' bases at other positions (TRAC's
        chr22:46177784 AGGC>A under 46177786 G>C and 46177787 C>T) puts two records of each
        base on a slot that holds them all: its bases cannot be read, and the slot does not
        hold the row, as a per-sample row reads its carriers (METHODS section 4). S00 holds
        the SNPs and the deletion on its first slot, S01 the SNPs alone; S02's slot holds a
        deletion at the first SNP's own position, a phasing artefact."""
        snps = [(106, "G", "C", {"S00": "1|0", "S01": "1|0", "S02": "1|0"}),
                (107, "C", "T", {"S00": "1|0", "S01": "1|0"})]
        indels = [(104, "AGGC", "A", {"S00": "1|0"}), (106, "GC", "G", {"S02": "1|0"})]
        ref = "ACGT" * 26 + "AGGCTTCA" + "ACGT" * 8
        p = self.panel(snps, indels)
        for seq_at in (None, lambda a, b: ref[max(a, 0):max(b, 0)]):
            p.seq_at = seq_at
            pair = hp.classify([v(106, "G", "C"), v(107, "C", "T")], [p], span=(100, 122))
            self.assertEqual((pair.tier, pair.state, pair.haplotypes, pair.samples),
                             ("A", hp.ONE_HAPLOTYPE, 1, ("S01",)))
        alone = self.panel([(106, "G", "C", {"S00": "1|0"})], indels[:1])
        ver = hp.classify([v(106, "G", "C")], [alone], span=(100, 122))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("C1", hp.INDEL_IN_TARGET, 1))
        changing = hp.changing_indels(p.indels, 99, 121)
        self.assertEqual(hp.breakers_of(changing, {106: "C", 107: "T"}, 99, 121),
                         [{"S00": "1|0"}])
        # the artefact at a SNP's own position takes nothing off, and S02's slot holds two
        # alleles of 106: phase not known
        self.assertEqual(hp.breakers_of(changing, {106: "C"}, 99, 121), [{"S00": "1|0"}])
        only_s02 = self.panel([(106, "G", "C", {"S02": "1|0"})], indels[1:])
        self.assertEqual(hp.classify([v(106, "G", "C")], [only_s02], span=(100, 122)).tier,
                         "B")

    def test_multi_alt_genotypes_are_read_per_alt(self):
        # one VCF line with two ALTs: S00 is G on slot 0 and T on slot 1 (1|2), listed under
        # both; with 105 C>T on slot 1, T and 105 are together, G and 105 are not
        p = self.panel([
            (100, "A", "G", {"S00": "1|2", "S01": "1|0"}),
            (100, "A", "T", {"S00": "1|2", "S02": "0|2"}),
            (105, "C", "T", {"S00": "0|1", "S02": "0|1"}),
        ])
        t = hp.classify([v(100, "A", "T"), v(105, "C", "T")], [p])
        self.assertEqual((t.tier, t.haplotypes, t.samples), ("A", 2, ("S00", "S02")))
        g = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p])
        self.assertEqual(g.tier, "C1")

    def test_two_panels_the_one_that_holds_them_decides(self):
        build_panel(self.dicts, [(100, "A", "G", {"S00": "1|0"}),
                                 (105, "C", "T", {"S00": "1|0"})], name="hg38_P1")
        build_panel(self.dicts, [(100, "A", "G", {"S03": "1|0"})], name="hg38_P2")
        p1 = hp.PanelChrom("hg38_P1", self.dicts, CHROM)
        p2 = hp.PanelChrom("hg38_P2", self.dicts, CHROM)
        self.addCleanup(p1.close)
        self.addCleanup(p2.close)
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p2, p1])
        self.assertEqual((ver.tier, ver.panel), ("A", "hg38_P1"))


class TestAVariantTheScoreDoesNotNeed(unittest.TestCase):
    """A row's worst case takes every allele that matches the guide best, and some of them
    leave its CFD, its mismatches and its PAM as they are (an allele at a PAM base the
    search leaves free, a mismatch for a mismatch the CFD table scores alike). A haplotype
    carrying the rest alone has the row's off-target at the row's score, so the row is
    tier A, not C1 (RV2-1); one whose allele the row's PAM needs, or whose CFD is another
    at the three decimals the results report, is not (RV3-1)."""

    # a 20-nt spacer + NGG, the alignment CFD is defined for; the ALT protospacer takes a
    # mismatch at spacer position 3 (103) and the PAM's N base (121), which CFD never reads
    GUIDE = "CTAACAGTTGCTTTTATCACNNN"
    REF = "CTAACAGTTGCTTTTATCACAGG"
    ALT = "CTTACAGTTGCTTTTATCACGGG"
    # a second guide whose position 7 scores dC and dT alike: both alleles leave the CFD
    GUIDE_TIE = "CTAACAATTGCTTTTATCACNNN"
    ALT_TIE = "CTAACACTTGCTTTTATCACGGG"  # REF has T at 107, the ALT C, scored the same
    # the search's PAM: the first base of an NRG PAM is free, its second is A or G
    PAM = hp.Pam("NRG", True)

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_score_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.dicts = os.path.join(self.tmp, "Dictionaries")

    def panel(self, snps, indels=()):
        build_panel(self.dicts, snps, indels)
        p = hp.PanelChrom(PANEL, self.dicts, CHROM)
        self.addCleanup(p.close)
        return p

    def alignment(self, guide=GUIDE, aligned=ALT, strand="+", start=100, cfd=None,
                  pam=PAM, block=hp.BLOCKS[0]):
        return hp.Alignment(1, block, guide, start, strand, guide, aligned, 3,
                            (), "NA", cfd, pam)

    def test_the_columns_of_a_variant_on_either_strand(self):
        plus = hp.variant_columns(100, self.ALT, "+",
                                  [v(103, "A", "T"), v(121, "A", "G")])
        self.assertEqual(plus, [(2, "A"), (20, "A")])
        # a minus-strand alignment is the reverse complement, its first column the last base
        minus = hp.variant_columns(100, self.ALT, "-",
                                   [v(121, "T", "A"), v(103, "T", "C")])
        self.assertEqual(minus, [(2, "A"), (20, "A")])
        # an allele the alignment does not take, and one outside it
        self.assertIsNone(hp.variant_columns(100, self.ALT, "+", [v(103, "A", "C")]))
        self.assertIsNone(hp.variant_columns(100, self.ALT, "+", [v(130, "A", "G")]))

    def test_which_variants_the_cfd_needs(self):
        placed = hp.variant_columns(100, self.ALT, "+", [v(103, "A", "T"), v(121, "A", "G")])
        self.assertEqual(hp.scoring_variants(self.GUIDE, self.ALT, placed, pam=self.PAM), [0])
        # the row's own CFD must be the one the gate reads, or it says nothing
        self.assertEqual(hp.scoring_variants(self.GUIDE, self.ALT, placed, shown=0.706,
                                             pam=self.PAM), [0])
        self.assertIsNone(hp.scoring_variants(self.GUIDE, self.ALT, placed, shown=0.5,
                                              pam=self.PAM))
        # another geometry (a 4-base PAM) is not CFD's
        self.assertIsNone(hp.alignment_cfd(self.GUIDE, self.ALT, pam_len=4))

    # a row of three mismatches whose CFD, 0.000, the results report as 0.000 whatever
    # happens at its PAM: its PAM is GAG, an NRG, and the reference has GCG, which is not
    ALT_PAM_MADE = "CTAACAGTTGCTGTGAGCACGAG"

    def test_an_allele_the_pam_needs_is_needed(self):
        """The guide's side of an alignment writes its PAM as N, which every base
        matches, and the CFD of a row scoring 0.000 moves within the decimals it reports
        whatever the PAM does: the allele that makes the PAM is read against the search's
        own PAM instead. Without it the sequence has no NRG PAM, so there is no off-target
        there at all and the row stays C1 (RV3-1; 461 rows of the genome-wide run,
        chr9:92157298(+) among them, were labelled "same score" for it)."""
        alt = self.ALT_PAM_MADE
        placed = hp.variant_columns(100, alt, "+", [v(113, "T", "G"), v(122, "C", "A")])
        self.assertEqual(placed, [(12, "T"), (21, "C")])
        # the CFD cannot see it: with the PAM's R base back at the reference both the row
        # and the sequence a haplotype without that allele has are reported as 0.000
        self.assertTrue(hp.same_cfd(hp.alignment_cfd(self.GUIDE, alt[:21] + "C" + alt[22:]),
                                    hp.alignment_cfd(self.GUIDE, alt)))
        self.assertEqual(hp.scoring_variants(self.GUIDE, alt, placed, shown=0.0,
                                             pam=self.PAM), [0, 1])
        p = self.panel([(113, "T", "G", {"S00": "1|0"}), (122, "C", "A", {"S02": "1|0"})])
        ver = hp.classify([v(113, "T", "G"), v(122, "C", "A")], [p],
                          alignment=self.alignment(aligned=alt, cfd=0.0))
        self.assertEqual(ver.tier, "C1")
        self.assertNotIn(ver.state, hp.SAME_SCORE_STATES)

    def test_another_allele_the_pam_needs_takes_the_slot(self):
        """The same of an allele of the panel the row does not take: the row's PAM is GAG,
        and the one slot carrying the variant its score needs also carries an allele at the
        PAM's R base that leaves GCG there, which is no NRG. That slot has no off-target at
        all, so the row stays C1 (score_breakers reads the PAM region too); the row's own
        allele at the PAM's free N base needs no haplotype."""
        alt = self.ALT_PAM_MADE  # the variant at column 20 is its N base, A -> G
        variants = [v(113, "T", "G"), v(121, "A", "G")]
        kept = self.panel([(113, "T", "G", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"})])
        ver = hp.classify(variants, [kept],
                          alignment=self.alignment(aligned=alt, cfd=0.0))
        self.assertEqual((ver.tier, ver.state, ver.kept),
                         ("A", hp.SAME_SCORE, (f"{CHROM}_113_T_G",)))
        # the same panel (a second store of its own) with an A -> C at the PAM's R base,
        # on the slot of S00, the sample carrying the variant the score needs
        self.dicts = os.path.join(self.tmp, "Dictionaries2")
        broken = self.panel([(113, "T", "G", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"}),
                             (122, "A", "C", {"S00": "1|0"})])
        # the CFD reports 0.000 with it as without it
        self.assertTrue(hp.same_cfd(hp.alignment_cfd(self.GUIDE, alt[:21] + "C" + alt[22:]),
                                    hp.alignment_cfd(self.GUIDE, alt)))
        ver = hp.classify(variants, [broken],
                          alignment=self.alignment(aligned=alt, cfd=0.0))
        self.assertEqual(ver.tier, "C1")

    def test_the_crispr_bulge_alignment_has_no_same_score_state(self):
        """The exemption is the CFD's. The gate reads no CRISPR-Bulge score (the model is
        not run again), and an allele the CFD never reads -- the PAM's N base -- does
        change that score, so the alignment the report's second ranking takes is C1 where
        no haplotype carries its combination (RV3-2: chr6:133354312(-) was 4th by
        CRISPR-Bulge at 0.0702, the haplotype its label named scoring 0.0124)."""
        p = self.panel([(103, "A", "T", {"S00": "1|0", "S01": "0|1"}),
                        (121, "A", "G", {"S02": "1|0"})])
        variants = [v(103, "A", "T"), v(121, "A", "G")]
        by_cfd = hp.classify(variants, [p], alignment=self.alignment(block=hp.BLOCKS[0]))
        self.assertEqual((by_cfd.tier, by_cfd.state), ("A", hp.SAME_SCORE))
        by_cb = hp.classify(variants, [p], alignment=self.alignment(block=hp.CB_BLOCK))
        self.assertEqual(by_cb.tier, "C1")
        self.assertNotIn(by_cb.state, hp.SAME_SCORE_STATES)

    def test_without_the_searchs_pam_every_allele_in_a_pam_region_is_needed(self):
        """The gate reads an alignment's PAM region against the search's PAM file. Run
        without it (--pam), it cannot say which bases make a PAM there and takes every
        allele of either end as one the row needs."""
        placed = hp.variant_columns(100, self.ALT, "+", [v(103, "A", "T"), v(121, "A", "G")])
        self.assertEqual(hp.scoring_variants(self.GUIDE, self.ALT, placed, pam=None), [0, 1])
        p = self.panel([(103, "A", "T", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"})])
        ver = hp.classify([v(103, "A", "T"), v(121, "A", "G")], [p],
                          alignment=self.alignment(pam=None))
        self.assertEqual(ver.tier, "C1")

    def test_a_cfd_that_differs_at_the_decimals_shown_is_a_variant_the_score_needs(self):
        """"The same CFD" is the CFD the integrated results report, of three decimals: an
        allele that moves it from 0.008 to 0.007 changes what the row says (RV3-1; the
        gate took any difference up to 5e-4 for none)."""
        guide = "CTAACAGTTGCTTTTATCA-CNNN"
        alt = "CTAATAAT-GATTTAATCAGCTGG"
        placed = [(4, "G"), (21, "C")]
        self.assertEqual(
            [round(hp.alignment_cfd(guide, alt), 6),
             round(hp.alignment_cfd(guide, alt[:4] + "G" + alt[5:]), 6)],
            [0.007654, 0.007216])
        # the allele at column 4 moves the CFD the results report; the other, at the PAM's
        # free first base, leaves the CFD and gives an NRG PAM all the same
        self.assertEqual(hp.scoring_variants(guide, alt, placed, shown=0.008, pam=self.PAM),
                         [0])

    def test_a_haplotype_with_the_scoring_variant_alone_is_tier_a(self):
        p = self.panel([(103, "A", "T", {"S00": "1|0", "S01": "0|1"}),
                        (121, "A", "G", {"S02": "1|0"})])
        variants = [v(103, "A", "T"), v(121, "A", "G")]
        # without the alignment the gate sees only that no haplotype carries both
        self.assertEqual(hp.classify(variants, [p]).tier, "C1")
        ver = hp.classify(variants, [p], alignment=self.alignment())
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("A", hp.SAME_SCORE, 2))
        self.assertEqual(ver.samples, ("S00", "S01"))
        self.assertEqual(ver.kept, (f"{CHROM}_103_A_T",))
        self.assertEqual(ver.label(), "the same CFD on 2 of 40 hg38_TESTPANEL haplotypes, "
                                      f"which carry {CHROM}_103_A_T: the other variant does "
                                      "not change it")

    def test_a_scoring_variant_nobody_carries_with_the_others_stays_c1(self):
        # the variant the score needs is in the panel, but on no haplotype that holds the
        # target: an indel of the store takes the target's bases off the one slot carrying it
        p = self.panel([(103, "A", "T", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"})],
                       indels=[(110, "GTT", "G", {"S00": "1|0"})])
        variants = [v(103, "A", "T"), v(121, "A", "G")]
        ver = hp.classify(variants, [p], span=(100, 122), alignment=self.alignment())
        self.assertEqual(ver.tier, "C1")
        self.assertNotIn(ver.state, hp.SAME_SCORE_STATES)
        # the same panel without that indel: the slot holds the target, so the row is A
        no_indel = hp.classify(variants, [p], alignment=self.alignment())
        self.assertEqual((no_indel.tier, no_indel.state), ("A", hp.SAME_SCORE))

    def test_every_variant_scored_alike_reads_as_the_reference(self):
        p = self.panel([(107, "T", "C", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"})])
        ver = hp.classify([v(107, "T", "C"), v(121, "A", "G")], [p],
                          alignment=self.alignment(self.GUIDE_TIE, self.ALT_TIE))
        self.assertEqual((ver.tier, ver.state), ("A", hp.SAME_SCORE_REF))
        self.assertEqual(ver.label(), "the same CFD as the reference here: none of the 2 "
                                      "variants changes it")

    def test_an_allele_that_makes_a_column_match_the_guide_is_needed(self):
        """An allele can leave the CFD as it is and still be needed: at position 1 the CFD
        table scores a mismatch 1.0, but without the allele the alignment has one mismatch
        more than the row reports, which the search's budget need not allow -- the row's
        off-target at the row's score is the one with it."""
        guide = "ATAACAGTTGCTTTTATCACNNN"
        alt = "ATAACAGTTGCTTTTATCACGGG"  # column 1 matches the guide; the reference has T
        placed = hp.variant_columns(100, alt, "+", [v(101, "T", "A"), v(121, "A", "G")])
        self.assertEqual(hp.alignment_cfd(guide, alt),
                         hp.alignment_cfd(guide, "T" + alt[1:]))
        self.assertEqual(hp.scoring_variants(guide, alt, placed, pam=self.PAM), [0])
        p = self.panel([(101, "T", "A", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"})])
        ver = hp.classify([v(101, "T", "A"), v(121, "A", "G")], [p],
                          alignment=self.alignment(guide, alt))
        self.assertEqual((ver.tier, ver.state, ver.kept),
                         ("A", hp.SAME_SCORE, (f"{CHROM}_101_T_A",)))

    def test_another_allele_in_the_window_that_changes_the_score_takes_the_slot(self):
        """A haplotype has the row's off-target at the row's score only if it holds nothing
        else inside the alignment that changes it: the slot carrying the variant the score
        needs also carries another allele of the panel there, so the row stays C1."""
        p = self.panel([(103, "A", "T", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"}),
                        (110, "T", "A", {"S00": "1|0"})])
        ver = hp.classify([v(103, "A", "T"), v(121, "A", "G")], [p],
                          alignment=self.alignment())
        self.assertEqual(ver.tier, "C1")
        # on the other slot of that sample it takes nothing off
        other = self.panel([(103, "A", "T", {"S00": "1|0"}), (121, "A", "G", {"S02": "1|0"}),
                            (110, "T", "A", {"S00": "0|1"})])
        kept = hp.classify([v(103, "A", "T"), v(121, "A", "G")], [other],
                           alignment=self.alignment())
        self.assertEqual((kept.tier, kept.state), ("A", hp.SAME_SCORE))

    def test_a_combination_every_variant_of_which_scores_stays_c1(self):
        # both alleles are in the spacer and change the CFD: no haplotype has the row's score
        p = self.panel([(103, "A", "T", {"S00": "1|0"}), (105, "T", "G", {"S02": "1|0"})])
        alt = "CTTAGAGTTGCTTTTATCACAGG"
        ver = hp.classify([v(103, "A", "T"), v(105, "T", "G")], [p],
                          alignment=self.alignment(aligned=alt))
        self.assertEqual(ver.tier, "C1")

    def test_the_search_writes_the_tier_and_the_scoring_variants(self):
        """The command the pipeline runs: the companion says which variants the CFD needs,
        and the sidecar counts the rows a haplotype reaches the CFD of. The search gives
        it its PAM file (--pam); run without it, no allele of a PAM region is read as one
        the row can do without, and the same row is C1."""
        build_panel(self.dicts, [(103, "A", "T", {"S00": "1|0", "S01": "0|1"}),
                                 (105, "T", "G", {"S02": "1|0"}),
                                 (121, "A", "G", {"S03": "1|0"})])
        out = os.path.join(self.tmp, "Results", "job")
        os.makedirs(out)
        vcfs = os.path.join(self.tmp, "list_vcf.txt")
        with open(vcfs, "w") as fh:
            fh.write(f"VCFs/{PANEL}\n")
        results = os.path.join(out, "G+NRG_hg38+" + PANEL + "_4+2_integrated_results.tsv")
        same = dict(crrna=self.GUIDE, pam="GGG")
        write_results(results, [
            # 1: the PAM's N base does not change the CFD, so the haplotypes carrying the
            # other allele alone have the row's off-target at its score
            {"Spacer+PAM": self.GUIDE, "Chromosome": CHROM,
             **block("highest_CFD", 100, "+", self.ALT, "alt",
                     f"{CHROM}_103_A_T,{CHROM}_121_A_G", "0.05,0.05", **same),
             **block("fewest_mm+b", 100, "+", self.ALT, "ref", "NA", **same),
             **block("highest_CRISPR_BULGE", 100, "+", self.ALT, "ref", "NA", **same)},
            # 2: both alleles change it, and no haplotype carries them together
            {"Spacer+PAM": self.GUIDE, "Chromosome": CHROM,
             **block("highest_CFD", 100, "+", "CTTAGAGTTGCTTTTATCACAGG", "alt",
                     f"{CHROM}_103_A_T,{CHROM}_105_T_G", "0.05,0.05",
                     crrna=self.GUIDE, pam="AGG"),
             **block("fewest_mm+b", 100, "+", self.ALT, "ref", "NA", **same),
             **block("highest_CRISPR_BULGE", 100, "+", self.ALT, "ref", "NA", **same)},
        ])
        pam_file = os.path.join(self.tmp, "20bp-NRG-SpCas9.txt")
        with open(pam_file, "w") as fh:
            fh.write("N" * 20 + "NRG 3\n")

        def run(*extra):
            done = subprocess.run(
                [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"),
                 "--results", results, "--dictionaries", self.dicts, "--vcf-list", vcfs,
                 "--mode", "population-level", *extra], capture_output=True, text=True)
            self.assertEqual((done.returncode, done.stderr), (0, ""))
            with open(hp.companion_path(results)) as fh:
                body = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#")]
            with open(os.path.join(out, hp.META_NAME)) as fh:
                return done.stdout, body, json.load(fh)["counts"]

        stdout, body, counts = run("--pam", pam_file)
        self.assertIn("A 1 (1 of them on a haplotype that reaches the row's CFD without "
                      "every variant)", stdout)
        tier, state = (hp.COMPANION_COLUMNS.index(c) for c in ("Tier", "State"))
        scoring = hp.COMPANION_COLUMNS.index("Scoring_variants")
        self.assertEqual([(r[0], r[tier], r[state], r[scoring]) for r in body],
                         [("1", "A", hp.SAME_SCORE, f"{CHROM}_103_A_T"),
                          ("2", "C1", hp.TOO_RARE, "-")])
        self.assertEqual(counts, {"A": 1, "C1": 1, hp.SAME_SCORE_COUNT: 1})
        _stdout, body, counts = run()
        self.assertEqual([(r[0], r[tier], r[state]) for r in body],
                         [("1", "C1", hp.TOO_RARE), ("2", "C1", hp.TOO_RARE)])
        self.assertEqual(counts, {"C1": 2})


def write_results(path, rows):
    """An integrated results table with the columns the gate reads (and a filler)."""
    cols = ["Spacer+PAM", "Chromosome"]
    for block in hp.BLOCKS:
        cols += [f"{c}_({block})" for c in hp._BLOCK_COLUMNS]
    cols.append("Annotation_GENCODE")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "NA") for c in cols) + "\n")


def block(name, start, strand, aligned, origin, variants, maf="NA", crrna=None, pam=None):
    """One alignment's columns; ``crrna`` and ``pam`` (the aligned spacer+PAM and the PAM)
    where the test wants the gate to read the alignment's CFD."""
    cols = {f"Start_coordinate_({name})": str(start), f"Strand_({name})": strand,
            f"Aligned_protospacer+PAM_ALT_({name})": aligned,
            f"REF/ALT_origin_({name})": origin,
            f"Variant_info_genome_({name})": variants, f"Variant_MAF_({name})": maf}
    if crrna is not None:
        cols[f"Aligned_spacer+PAM_({name})"] = crrna
    if pam is not None:
        cols[f"PAM_({name})"] = pam
    return cols


class TestGateCommand(unittest.TestCase):
    """The command the pipeline runs: the companion, the sidecar, and an integrated results
    table left byte for byte as it was."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_cli_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        build_panel(self.dicts, [
            (100, "A", "G", {"S00": "1|0", "S01": "1|0"}),
            (105, "C", "T", {"S00": "1|0", "S02": "0|1"}),
            (110, "G", "A", {"S03": "0|1"}),
        ])
        self.out = os.path.join(self.tmp, "Results", "job")
        os.makedirs(self.out)
        self.vcfs = os.path.join(self.tmp, "list_vcf.txt")
        with open(self.vcfs, "w") as fh:
            fh.write(f"VCFs/{PANEL}\n")
        self.results = os.path.join(
            self.out, "GUIDE+NRG_hg38+" + PANEL + "_4+2_integrated_results.tsv")
        aln = "ACGTACGTACGTACGTACGTAGG"
        two = f"{CHROM}_100_A_G,{CHROM}_105_C_T"
        trans = f"{CHROM}_105_C_T,{CHROM}_110_G_A"
        rows = [
            # 1: reference row
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 90, "+", aln, "ref", "NA"),
             **block("fewest_mm+b", 90, "+", aln, "ref", "NA"),
             **block("highest_CRISPR_BULGE", 90, "+", aln, "ref", "NA")},
            # 2: CFD alignment A, CRISPR-Bulge alignment C1
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 90, "+", aln, "alt", two, "0.05,0.05"),
             **block("fewest_mm+b", 90, "+", aln, "alt", f"{CHROM}_100_A_G", "0.05"),
             **block("highest_CRISPR_BULGE", 90, "-", aln, "alt", trans, "0.05,0.025")},
            # 3: a variant no panel holds
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 200, "+", aln, "alt", f"{CHROM}_100_A_G,{CHROM}_210_T_C"),
             **block("fewest_mm+b", 200, "+", aln, "ref", "NA"),
             **block("highest_CRISPR_BULGE", 200, "+", aln, "ref", "NA")},
        ]
        write_results(self.results, rows)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_gate(self, *extra):
        cmd = [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"),
               "--results", self.results, "--dictionaries", self.dicts,
               "--vcf-list", self.vcfs, *extra]
        return subprocess.run(cmd, capture_output=True, text=True)

    def md5(self):
        with open(self.results, "rb") as fh:
            return hashlib.md5(fh.read()).hexdigest()

    def test_population_level_writes_the_companion_and_leaves_the_rows(self):
        before = self.md5()
        done = self.run_gate("--mode", "population-level")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(done.stderr, "")
        self.assertIn("Carrier gate: the rows' variant combinations checked", done.stdout)
        self.assertEqual(self.md5(), before)  # the off-target set is the same file
        comp = hp.companion_path(self.results)
        with open(comp) as fh:
            lines = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("# ")]
        header, body = lines[0], lines[1:]
        self.assertEqual(header[0], "#Row")
        got = [(r[0], r[1], r[7], r[8]) for r in body]
        self.assertEqual(got, [
            ("2", "highest_CFD", "A", hp.ONE_HAPLOTYPE),
            ("2", "fewest_mm+b", "B", "-"),
            ("2", "highest_CRISPR_BULGE", "C1", hp.TOO_RARE),
            ("3", "highest_CFD", "C2", "-"),
        ])
        with open(os.path.join(self.out, hp.META_NAME)) as fh:
            meta = json.load(fh)
        self.assertEqual(meta["status"], "run")
        self.assertEqual(meta["counts"], {"A": 1, "C2": 1})
        self.assertEqual(meta["panels"][0]["name"], PANEL)
        self.assertEqual(meta["panels"][0]["samples"], 20)
        self.assertEqual(meta["panels"][0]["superpopulations"], {"AFR": 12, "EUR": 8})
        self.assertEqual(meta["companion"], os.path.basename(comp))

    def test_off_and_per_sample_write_no_companion(self):
        for extra, status in ((("--mode", "population-level", "--prioritize-observed", "off"),
                               "off"), (("--mode", "per-sample"), "per-sample")):
            done = self.run_gate(*extra)
            self.assertEqual((done.returncode, done.stderr), (0, ""))
            with open(os.path.join(self.out, hp.META_NAME)) as fh:
                self.assertEqual(json.load(fh)["status"], status)
            self.assertFalse(os.path.exists(hp.companion_path(self.results)))


class TestTheSearchGivesTheGateItsPam(unittest.TestCase):
    """The pipeline's carrier-gate step passes the search's PAM file: without it the gate
    reads no PAM region of an alignment, and takes every allele in one for an allele the
    row's score needs (RV3-1)."""

    def test_the_gate_step_passes_the_pam_file(self):
        with open(os.path.join(HERE, "submit_job_automated_new_multiple_vcfs.sh")) as fh:
            lines = fh.read().splitlines()
        i = next(n for n, ln in enumerate(lines)
                 if "haplotype_prioritizer.py" in ln and "python " in ln)
        step = ""
        while True:
            step += lines[i]
            if not lines[i].rstrip().endswith("\\"):
                break
            i += 1
        self.assertIn('--pam "$pam_file"', step)


class TestTheFrequencyBoundOfASitesOnlyRow(unittest.TestCase):
    """A sites-only dataset's row of a SNP and an indel the genotype oracle does not hold
    (C2) shows the frequency bound of both: the SNP's frequency from the dataset's registry,
    the indel's from its indel frequencies, as the mega index's builders write them."""

    def setUp(self):
        import build_mega_indel_af as bia
        import build_mega_registry as bmr
        self.tmp = tempfile.mkdtemp(prefix="gate_bound_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        dicts = os.path.join(self.tmp, "Dictionaries")
        build_panel(dicts, [(300, "A", "G", {"S00": "1|0"})], name=hp.ORACLE_PANEL)
        vcf = os.path.join(self.tmp, "mega.chr22.vcf.gz")
        with gzip.open(vcf, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write(f"{CHROM}\t100\trs1\tA\tG\t.\tPASS\tAF_gnomAD=0.05;AF_max=0.05\n")
            fh.write(f"{CHROM}\t104\trs2\tCTT\tC\t.\tPASS\tAF_gnomAD=0.02;AF_max=0.02\n")
        bmr.build(vcf, CHROM, os.path.join(dicts, "registry_hg38_mega"))
        bia.build_indel_af_store(vcf, CHROM, os.path.join(dicts, "indel_af_hg38_mega"))
        out = os.path.join(self.tmp, "Results", "job")
        os.makedirs(out)
        vcfs = os.path.join(self.tmp, "list_vcf.txt")
        with open(vcfs, "w") as fh:
            fh.write("VCFs/hg38_mega\n")
        self.results = os.path.join(out, "GUIDE+NRG_hg38+hg38_mega_4+2_integrated_results.tsv")
        aln = "ACGTACGTACGTACGTACGTAGG"
        # the row's MAF gives no indel frequency (a sites-only index's indel log has none)
        write_results(self.results, [
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 90, "+", aln, "alt", f"{CHROM}_100_A_G,{CHROM}_104_CTT_C",
                     "NA,0.05"),
             **block("fewest_mm+b", 90, "+", aln, "ref", "NA"),
             **block("highest_CRISPR_BULGE", 90, "+", aln, "ref", "NA")}])
        self.cmd = [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"),
                    "--results", self.results, "--dictionaries", dicts, "--vcf-list", vcfs,
                    "--mode", "population-level", "--prioritize-observed", "on"]

    def test_the_snp_and_the_indel_frequency(self):
        done = subprocess.run(self.cmd, capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        with open(hp.companion_path(self.results)) as fh:
            body = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#")]
        self.assertEqual([(r[1], r[7]) for r in body], [("highest_CFD", "C2")])
        label = hp.COMPANION_COLUMNS.index("Label")
        self.assertEqual(body[0][label], "not checked: 2 of 2 variants not in hg38_HGDP1kGP; "
                                      "frequency product 1.0e-03, rarest 2.0e-02 "
                                      "(within the 1e-06 / 0.001 bound)")


class TestReportRanksByTheGate(unittest.TestCase):
    """generate_report on a population-level run the gate checked: the C1 row leaves the
    top-1000, the panel and the risk card for its own table, the A and C2 rows keep their
    place with their tier as their Observed cell, and the rows themselves are the search's."""

    GUIDE = "CTCTCAGCTGGTACACGGCANNN"
    ROWS = [
        # chrom pos  mm  cfd   cb     origin variants (CFD and CRISPR-Bulge alignment)
        ("1000", "0", "1.0", "1.0", "ref", "NA"),
        ("500", "3", "0.30", "0.10", "ref", "NA"),
        ("1000", "3", "0.50", "0.40", "alt", f"{CHROM}_1010_A_G,{CHROM}_1015_C_T"),  # A
        ("2000", "3", "0.95", "0.90", "alt", f"{CHROM}_2010_G_A,{CHROM}_2015_T_C"),  # C1
        ("3000", "3", "0.40", "0.30", "alt", f"{CHROM}_3010_A_C,{CHROM}_3020_G_T"),  # C2
        ("4000", "3", "0.20", "0.20", "alt", f"{CHROM}_4010_C_G"),                   # B
    ]

    def setUp(self):
        import generate_report  # noqa: F401  (pandas + matplotlib: installed in CI)
        self.tmp = tempfile.mkdtemp(prefix="gate_report_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        build_panel(self.dicts, [
            (1010, "A", "G", {"S00": "1|0", "S01": "1|0"}),
            (1015, "C", "T", {"S00": "1|0", "S01": "1|0", "S02": "0|1"}),
            (2010, "G", "A", {s: "1|0" for s in SAMPLES}),
            (2015, "T", "C", {s: "0|1" for s in SAMPLES}),
            (3010, "A", "C", {"S04": "1|0"}),
            (4010, "C", "G", {"S05": "0|1"}),
            # the pair of the CRISPR-Bulge test below: on no one haplotype
            (5011, "A", "G", {"S00": "1|0"}),
            (5021, "A", "T", {"S02": "1|0"}),
        ])
        self.rd = os.path.join(self.tmp, "Results", "job")
        os.makedirs(self.rd)
        with open(os.path.join(self.rd, ".search_mode"), "w") as fh:
            fh.write("population-level")
        with open(os.path.join(self.rd, ".Params.txt"), "w") as fh:
            fh.write("Mismatches\t4\nDNA\t1\nRNA\t1\nMax_bulges\t2\nPam\tNGG\n")
        self.vcfs = os.path.join(self.tmp, "list_vcf.txt")
        with open(self.vcfs, "w") as fh:
            fh.write(f"{PANEL}\n")
        self.results = os.path.join(
            self.rd, f"{self.GUIDE}+NGG_hg38+{PANEL}_4+2_integrated_results.tsv")
        self.cols = ["Spacer+PAM", "Chromosome", "Not_found_in_REF", "Annotation_GENCODE"]
        per_block = ["Start_coordinate", "Strand", "Aligned_spacer+PAM",
                     "Aligned_protospacer+PAM_REF", "Aligned_protospacer+PAM_ALT", "PAM",
                     "Mismatches", "Bulges", "Mismatches+bulges", "REF/ALT_origin",
                     "PAM_creation", "Variant_info_genome", "Variant_MAF", "Variant_rsID",
                     "Variant_samples"]
        for b in hp.BLOCKS:
            self.cols += [f"{c}_({b})" for c in per_block]
        self.cols += ["CFD_score_(highest_CFD)", "CFD_score_REF_(highest_CFD)",
                      "CFD_score_ALT_(highest_CFD)",
                      "CRISPR_BULGE_score_(highest_CRISPR_BULGE)",
                      "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)",
                      "CRISPR_BULGE_score_ALT_(highest_CRISPR_BULGE)"]
        self.pam_file = os.path.join(self.tmp, "20bp-NGG-SpCas9.txt")
        with open(self.pam_file, "w") as fh:
            fh.write("N" * 20 + "NGG 3\n")
        self.write_results(self.ROWS)

    def write_results(self, rows):
        with open(self.results, "w") as fh:
            fh.write("\t".join(self.cols) + "\n")
            for row in rows:
                self.write_row(fh, *row)

    def gate(self):
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"), "--results",
             self.results, "--dictionaries", self.dicts, "--vcf-list", self.vcfs,
             "--mode", "population-level", "--pam", self.pam_file],
            capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))

    def write_row(self, fh, pos, mm, cfd, cb, origin, variants, cb_variants=None,
                  aln="CTCTCAGCTGGTACACGGCATGG"):
        """One row of the integrated results: ``variants`` are those of each alignment but
        the highest-CRISPR-Bulge one when ``cb_variants`` is given; ``aln`` its aligned
        protospacer+PAM."""
        r = {"Spacer+PAM": self.GUIDE, "Chromosome": CHROM,
             "Not_found_in_REF": "y" if origin == "alt" else "n",
             "Annotation_GENCODE": "intergenic"}
        for b in hp.BLOCKS:
            r.update({f"Start_coordinate_({b})": pos, f"Strand_({b})": "+",
                      f"Aligned_spacer+PAM_({b})": self.GUIDE,
                      f"Aligned_protospacer+PAM_REF_({b})": aln,
                      f"Aligned_protospacer+PAM_ALT_({b})": aln if origin == "alt" else "NA",
                      f"PAM_({b})": "TGG", f"Mismatches_({b})": mm, f"Bulges_({b})": "0",
                      f"Mismatches+bulges_({b})": mm, f"REF/ALT_origin_({b})": origin,
                      f"PAM_creation_({b})": "NA",
                      f"Variant_info_genome_({b})": cb_variants
                      if cb_variants and b == "highest_CRISPR_BULGE" else variants,
                      f"Variant_MAF_({b})": "0.05" if origin == "alt" else "NA",
                      f"Variant_rsID_({b})": ".", f"Variant_samples_({b})":
                          "S00,S01,S02" if origin == "alt" else "NA"})
        r.update({"CFD_score_(highest_CFD)": cfd, "CFD_score_REF_(highest_CFD)": "0.1",
                  "CFD_score_ALT_(highest_CFD)": cfd,
                  "CRISPR_BULGE_score_(highest_CRISPR_BULGE)": cb,
                  "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": "0.05",
                  "CRISPR_BULGE_score_ALT_(highest_CRISPR_BULGE)": cb})
        fh.write("\t".join(r.get(c, "NA") for c in self.cols) + "\n")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def report(self):
        import generate_report as gr
        import zipfile
        out = gr.build_report(result_dir=self.rd, integrated_tsv=self.results,
                              report_full=True)
        files = {}
        with zipfile.ZipFile(out) as zf:
            for name in zf.namelist():
                data = zf.read(name)
                files[name] = gzip.decompress(data) if name.endswith(".gz") else data
        return {k: v.decode() for k, v in files.items()}

    def positions(self, tsv):
        lines = tsv.splitlines()
        col = lines[0].split("\t").index("Position")
        return [ln.split("\t")[col] for ln in lines[1:]]

    def test_c1_leaves_the_rankings_for_its_own_table(self):
        with open(self.results, "rb") as fh:
            source = fh.read()
        # unchecked (no gate), the rows rank by CFD alone
        files = self.report()
        self.assertEqual(self.positions(files["data/top1000.tsv"]),
                         ["2000", "1000", "3000", "500", "4000"])
        self.assertNotIn("data/unobserved_worst_case.tsv", files)
        self.assertNotIn("Carrier gate", files["report.html"])
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"), "--results",
             self.results, "--dictionaries", self.dicts, "--vcf-list", self.vcfs,
             "--mode", "population-level"], capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        files = self.report()
        top = self.positions(files["data/top1000.tsv"])
        self.assertEqual(top, ["1000", "3000", "500", "4000"])  # by CFD, without 2000
        self.assertNotIn("2000", self.positions(files["data/panel_top100.tsv"]))
        self.assertNotIn("2000", self.positions(files["data/top1000_crispr_bulge.tsv"]))
        self.assertEqual(self.positions(files["data/unobserved_worst_case.tsv"]), ["2000"])
        html = files["report.html"]
        self.assertIn("Carrier gate", html)
        self.assertIn("variant combinations no haplotype of the panel carries", html)
        self.assertIn(f"together on 2 of 40 {PANEL} haplotypes", html)
        self.assertIn(f"not together on any of 40 {PANEL} haplotypes (expected 10)", html)
        self.assertIn(f"not checked: 1 of 2 variants not in {PANEL}", html)
        self.assertIn("AFR 12, EUR 8", html)
        self.assertNotIn("chimera", html.split("Carrier gate", 1)[1].split("</tr>", 1)[0])
        manifest = json.loads(files["data/run_manifest.json"])
        self.assertEqual(manifest["carrier_gate"]["tiers"], {"A": 1, "B": 1, "C1": 1, "C2": 1})
        # the off-target set is the search's: the bundled table is the source file
        self.assertEqual(files["data/integrated_results.tsv.gz"].encode(), source)
        self.assertIn("data/carrier_gate.tsv.gz", files)

    def test_a_one_variant_row_the_gate_demoted_reads_its_label_and_the_legend_has_it(self):
        """A row of one variant every carrying haplotype of which has an indel in the
        target is C1: its Observed cell gives the gate's label, in place of its carriers,
        and the report's legend and the inputs card say so (RV2-3)."""
        build_panel(self.dicts, [
            (1010, "A", "G", {"S00": "1|0", "S01": "1|0"}),
            (1015, "C", "T", {"S00": "1|0", "S01": "1|0", "S02": "0|1"}),
            (2010, "G", "A", {s: "1|0" for s in SAMPLES}),
            (2015, "T", "C", {s: "0|1" for s in SAMPLES}),
            (3010, "A", "C", {"S04": "1|0"}),
            (4010, "C", "G", {"S05": "0|1"}),
            (7010, "C", "T", {"S06": "1|0"}),
        ], indels=[(7015, "GTT", "G", {"S06": "1|0"})])
        with open(self.results, "a") as fh:
            self.write_row(fh, "7000", "3", "0.45", "0.35", "alt", f"{CHROM}_7010_C_T")
        self.gate()
        files = self.report()
        demoted = files["data/unobserved_worst_case.tsv"]
        head = demoted.splitlines()[0].split("\t")
        pos, obs = head.index("Position"), head.index("Observed")
        rows = {ln.split("\t")[pos]: ln.split("\t")[obs] for ln in demoted.splitlines()[1:]}
        self.assertEqual(rows["7000"],
                         f"the 1 {PANEL} haplotype carrying it has an indel in the target")
        html = files["report.html"]
        self.assertIn("A row of one variant reads its carriers, unless every haplotype "
                      "carrying it has an indel that takes the target's bases off it",
                      html)
        self.assertIn("the tier of a row the gate checked &mdash; a row of two or more "
                      "variants, and a row of one whose every carrying haplotype has an "
                      "indel in the target", html)

    def test_the_plots_of_section_2_rank_what_their_tables_do(self):
        """Panels (c) and (d) of Section 2 rank by CRISPR-Bulge score, as the CRISPR-Bulge
        table does, so the plot and the table show the same sites: a row whose highest-CFD
        alignment is C1, and whose highest-CRISPR-Bulge alignment is not, stays in both."""
        import generate_report as gr
        from unittest import mock
        with open(self.results, "a") as fh:
            self.write_row(fh, "6000", "3", "0.70", "0.70", "alt",
                           f"{CHROM}_2010_G_A,{CHROM}_2015_T_C",
                           cb_variants=f"{CHROM}_4010_C_G")
        self.gate()
        seen = []

        def record(sub, *args, **kwargs):
            seen.append(list(sub["Start_coordinate_(highest_CFD)"]))
            return "data:,"

        with mock.patch.object(gr, "_cfd_style_scatter", record):
            files = self.report()
        self.assertEqual(len(seen), 4)  # (a) CFD, (b) CFD delta, (c) and (d) CRISPR-Bulge
        by_cfd, _delta, by_cb, by_cb_delta = seen
        self.assertNotIn("6000", by_cfd)
        self.assertIn("6000", by_cb)
        self.assertIn("6000", self.positions(files["data/top1000_crispr_bulge.tsv"]))
        self.assertNotIn("6000", self.positions(files["data/top1000.tsv"]))
        # each plot shows the sites of the table ranked by the same score
        self.assertEqual(sorted(by_cfd), sorted(self.positions(files["data/top1000.tsv"])))
        for plot in (by_cb, by_cb_delta):
            self.assertEqual(sorted(plot),
                             sorted(self.positions(files["data/top1000_crispr_bulge.tsv"])))

    def section4(self, files):
        """The validation panel's totals, its threshold files and run_manifest.json's
        counts."""
        import re
        html = files["report.html"]
        bridge = re.search(r"Off-targets \(([\d,]+) total sites &minus; ([\d,]+) perfect "
                           r"match(?:es)? = ([\d,]+)\)", html).groups()
        created = re.search(r"Variant-created off-targets \(Not_found_in_REF\)</td>"
                            r"<td class='num'>([\d,]+)", html).group(1)
        tiers = {name[len("data/"):]: sorted(self.positions(tsv)) for name, tsv in files.items()
                 if name.startswith(("data/cfd_ge_", "data/crispr_bulge_ge_",
                                     "data/variant_created"))}
        counts = json.loads(files["data/run_manifest.json"])["counts"]
        return bridge, created, tiers, counts

    def test_the_totals_count_every_row(self):
        """The validation panel's off-targets, perfect matches and variant-created
        off-targets, and run_manifest.json's counts, are the search's: the rows the gate
        ranks after the others are among them, as in the Section-1 matrix, and a row says
        how many there are. Its thresholds leave them out, as the risk card does: the CFD
        ones the rows whose highest-CFD alignment is C1, the CRISPR-Bulge ones those whose
        highest-CRISPR-Bulge alignment is (here 5000's, whose CFD alignment is one
        variant)."""
        with open(self.results, "a") as fh:
            self.write_row(fh, "5000", "2", "0.60", "0.80", "alt", f"{CHROM}_4010_C_G",
                           cb_variants=f"{CHROM}_2010_G_A,{CHROM}_2015_T_C")
        bridge, created, tiers, counts = self.section4(self.report())
        self.assertEqual(bridge, ("7", "1", "6"))
        self.assertEqual(created, "5")
        self.assertEqual(tiers["cfd_ge_0.50.tsv"], ["1000", "2000", "5000"])
        self.assertEqual(tiers["crispr_bulge_ge_0.50.tsv"], ["2000", "5000"])
        self.assertEqual(tiers["variant_created.tsv"], ["1000", "2000", "3000", "4000", "5000"])
        self.gate()
        files = self.report()
        gated = self.section4(files)
        self.assertEqual(gated[:2], (bridge, created))
        self.assertEqual({k: v for k, v in gated[3].items() if k != "recommended_panel_size"},
                         {k: v for k, v in counts.items() if k != "recommended_panel_size"})
        self.assertEqual(gated[3]["variant_created_off_targets"], 5)
        self.assertEqual(gated[2], {"cfd_ge_0.50.tsv": ["1000", "5000"],
                                    "cfd_ge_0.20.tsv": ["1000", "3000", "4000", "500", "5000"],
                                    "cfd_ge_0.05.tsv": ["1000", "3000", "4000", "500", "5000"],
                                    "crispr_bulge_ge_0.20.tsv": ["1000", "3000", "4000"],
                                    "variant_created.tsv": tiers["variant_created.tsv"]})
        html = files["report.html"]
        self.assertIn("Ranked after every other row by the carrier gate (tier C1: no haplotype "
                      "of the panel has the sequence of the row's highest-CFD alignment; their "
                      "table follows the top-1000, in Section 6): counted above, left out of the "
                      "CFD and edit-distance thresholds, the panel and the risk card&rsquo;s CFD "
                      "bands</td><td class='num'>1</td>", html)
        self.assertIn("&hellip; at the highest-CRISPR-Bulge alignment: left out of the "
                      "CRISPR-Bulge thresholds and the risk card&rsquo;s CRISPR-Bulge "
                      "bands</td><td class='num'>2</td>", html)
        self.assertIn("5000", self.positions(files["data/top1000.tsv"]))
        self.assertNotIn("5000", self.positions(files["data/top1000_crispr_bulge.tsv"]))

    def test_the_crispr_bulge_views_keep_no_row_by_a_same_score_verdict(self):
        """A row whose two alleles no haplotype carries together, one of them at the PAM's
        N base: a haplotype carrying the other alone has the row's CFD, so the row keeps
        its place in the views ranked by CFD; the CRISPR-Bulge score of such a haplotype
        is another, which the gate does not read, so the row leaves the views ranked by
        that score (RV3-2)."""
        # three mismatches with the guide, CFD 0.183; its variants are at the guide's own
        # base in column 10 (the score needs it: without it the alignment has a mismatch
        # more) and at the PAM's N base
        self.write_results(self.ROWS + [
            ("5000", "3", "0.183", "0.90", "alt", f"{CHROM}_5011_A_G,{CHROM}_5021_A_T",
             None, "CACACCGCTGGTACACGGCATGG")])
        self.gate()
        with open(hp.companion_path(self.results)) as fh:
            body = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#")]
        tier, state = (hp.COMPANION_COLUMNS.index(c) for c in ("Tier", "State"))
        mine = [(r[1], r[tier], r[state]) for r in body if r[4] == "5000"]
        self.assertEqual(mine, [("highest_CFD", "A", hp.SAME_SCORE),
                                ("fewest_mm+b", "A", hp.SAME_SCORE),
                                ("highest_CRISPR_BULGE", "C1", hp.TOO_RARE)])
        files = self.report()
        self.assertIn("5000", self.positions(files["data/top1000.tsv"]))
        self.assertNotIn("5000", self.positions(files["data/top1000_crispr_bulge.tsv"]))

    def test_a_c1_row_whose_crispr_bulge_alignment_is_carried(self):
        """2000's highest-CRISPR-Bulge alignment one variant (B): it leaves the CFD
        thresholds but stays in the CRISPR-Bulge ones, which then leave no row out."""
        self.write_results([r + (f"{CHROM}_4010_C_G",) if r[0] == "2000" else r
                            for r in self.ROWS])
        self.gate()
        files = self.report()
        _bridge, _created, tiers, _counts = self.section4(files)
        self.assertEqual(tiers["cfd_ge_0.50.tsv"], ["1000"])
        self.assertEqual(tiers["crispr_bulge_ge_0.50.tsv"], ["2000"])
        self.assertIn("2000", self.positions(files["data/top1000_crispr_bulge.tsv"]))
        self.assertNotIn("at the highest-CRISPR-Bulge alignment: left out", files["report.html"])


class _PostAnalysis(unittest.TestCase):
    """The real new_simple_analysis.py on one variant window of chrT, with the stores the
    repository's writers build: the Tier-0 registry + Tier-1 genotypes (the dict-less
    tiers), or a classic per-sample dictionary, and the phased indel genotypes."""

    GUIDE = "CTAACAGTTGCTTTTATCAC"
    FLANK = "ACGTGCATTC" * 10

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_nsa_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.at = len(self.FLANK)
        self.chrom = "chrT"

    def world(self, site, snps, indels, tiers, right=None, chrom="chrT", meta=SAMPLE_META):
        """``snps`` [(column, REF, ALT, {sample: gt})] and ``indels`` [(column, REF, ALT,
        {sample: gt})], columns of ``site``; ``tiers``: "tiers" or "dict"; ``right``: the
        flank after the site (by default the one before it); ``chrom`` and ``meta``: the
        chromosome and the panel's samples (on chrX a male's calls are haploid)."""
        import tier0_compile as t0c
        import tier0_registry as t0
        self.chrom = chrom
        self.fasta = os.path.join(self.tmp, f"{chrom}.fa")
        with open(self.fasta, "w") as fh:
            fh.write(f">{chrom}\n" + self.FLANK + site + (right or self.FLANK) + "\n")
        d = os.path.join(self.tmp, "Dictionaries")
        recs = [(self.at + c + 1, r, a, f"rs{c}", dict(g)) for c, r, a, g in snps]
        if tiers == "tiers":
            rd = os.path.join(d, "registry_ds")
            gd = os.path.join(d, "genotypes_ds")
            os.makedirs(rd)
            os.makedirs(gd)
            t0.compile_registry_panel(recs, meta, None, t0c.ploidy_of_for_chrom(chrom),
                                      os.path.join(rd, f"reg_{chrom}.bin"),
                                      os.path.join(rd, f"reg_{chrom}.idx"))
            t1.compile_genotypes([(p, a, g, r) for p, r, a, _rs, g in recs],
                                 t1.build_sample_axis(meta),
                                 os.path.join(gd, f"gt_{chrom}.bin"),
                                 os.path.join(gd, f"gt_{chrom}.idx"))
        else:
            dd = os.path.join(d, "dictionaries_ds")
            os.makedirs(dd)
            entries = {}
            for p, r, a, rs, g in recs:
                entry = "%s;%s,%s;%s;%s" % (
                    ",".join("%s:%s" % kv for kv in sorted(g.items())), r, a, rs, "0.05")
                key = f"{chrom},{p}"
                entries[key] = entries[key] + "$" + entry if key in entries else entry
            with open(os.path.join(dd, f"my_dict_{chrom}.json"), "w") as fh:
                json.dump(entries, fh)
        if indels:
            igd = os.path.join(d, "indel_genotypes_ds")
            os.makedirs(igd)
            with gzip.open(os.path.join(igd, f"gt_indel_{chrom}.tsv.gz"), "wt") as fh:
                for c, r, a, g in indels:
                    fh.write("%d_%s_%s\t%s\n" % (self.at + c + 1, r, a, ",".join(
                        "%s:%s" % kv for kv in sorted(g.items()))))
        self.dict_path = os.path.join(d, "dictionaries_ds", f"my_dict_{chrom}.json")
        self.pam = os.path.join(self.tmp, "pam.txt")
        with open(self.pam, "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")

    def analyse(self, window, mode, budget, extra_env=None):
        """The rows of new_simple_analysis.py for the variant search's row ``window`` (the
        site with an IUPAC code at each SNP column): [(DNA, mismatches, samples, SNPs)]."""
        targets = os.path.join(self.tmp, f"targets_{mode}.txt")
        row = ["X", self.GUIDE + "NNN", window, self.chrom, str(self.at), str(self.at), "+",
               "0", "0", "0", "NA", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
        with open(targets, "w") as fh:
            fh.write("#header\n" + "\t".join(row) + "\n")
        out = os.path.join(self.tmp, f"out_{mode}")
        env = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
        env.update(CRISPRME_FAST_MODE=mode, CRISPRME_SKIP_SCORER_ENV="1", **(extra_env or {}))
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "new_simple_analysis.py"), self.fasta,
             self.dict_path, targets, self.pam, out, str(budget)],
            cwd=HERE, env=env, capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""), done.stdout)
        with open(out + ".bestCFD.txt") as fh:
            lines = [ln.rstrip("\n").split("\t") for ln in fh][1:]
        return sorted((f[2].upper(), int(f[7]), f[12], f[17]) for f in lines)


class TestPerSampleIndelInTheTarget(_PostAnalysis):
    """F2: a per-sample SNP row named a carrier whose haplotype holding the SNP also holds
    an indel inside the target (its sequence is not the row's). S00 holds the SNP and a
    deletion on its first slot; S02 has the deletion on the other slot; S03's deletion call
    is unphased; S01 has no indel."""

    # the guide with a mismatch at column 4 (A for C), which the SNP restores
    SITE = "CTAAAAGTTGCTTTTATCACTGG"
    SNPS = [(4, "A", "C", {"S00": "1|0", "S01": "1|0", "S02": "0|1", "S03": "1|0"})]
    # a deletion of the T at column 13, inside the target
    INDELS = [(12, "TT", "T", {"S00": "1|0", "S02": "1|0", "S03": "1/0"})]
    WINDOW = SITE[:4] + "M" + SITE[5:]

    def check(self, tiers):
        self.world(self.SITE, self.SNPS, self.INDELS, tiers)
        per_sample = self.analyse(self.WINDOW, "0", 0)
        # (the classic dictionary's path gives a row per haplotype slot, which the merge
        # step joins)
        self.assertEqual({r[:2] + (r[3],) for r in per_sample},
                         {(self.GUIDE + "TGG", 0, f"chrT_{self.at + 5}_A_C")})
        self.assertEqual(sorted(s for r in per_sample for s in r[2].split(",")),
                         ["S01", "S02", "S03"])
        # a population-level row lists every carrier of its alleles; the carrier gate
        # reads the indel genotypes for it
        population = self.analyse(self.WINDOW, "1", 0)
        self.assertEqual(sorted(s for r in population for s in r[2].split(",")),
                         ["S00", "S01", "S02", "S03"])

    def test_dict_less_tiers(self):
        self.check("tiers")

    def test_classic_dictionary(self):
        self.check("dict")


class TestPerSampleIndelInTheTargetOfAMale(_PostAnalysis):
    """F2 on a haploid call: a male's chrX genotypes are his one haplotype's, so an indel
    of his that changes the target takes him off the row, and his SNPs there are never a
    haplotype of their own. S00 is male and carries the row's SNP, a second SNP of the
    window and a deletion inside the target; S01 (female) carries the row's SNP alone
    (RV3-3: the fix took only phased calls off, leaving 778 haploid carrier listings in
    the per-sample genome-wide run)."""

    META = {s: ("TEST", "AFR", "male" if s == "S00" else "female") for s in SAMPLES}
    SITE = "CTAAAAGTTGCTTTTATCACTGG"
    # the row's SNP at column 4, a second at column 17, a deletion of the T at column 13
    SNPS = [(4, "A", "C", {"S00": "1", "S01": "1|0"}), (17, "T", "A", {"S00": "1"})]
    INDELS = [(12, "TT", "T", {"S00": "1"})]
    WINDOW = SITE[:4] + "M" + SITE[5:17] + "W" + SITE[18:]

    def check(self, tiers):
        self.world(self.SITE, self.SNPS, self.INDELS, tiers, chrom="chrX", meta=self.META)
        rows = self.analyse(self.WINDOW, "0", 2)
        named = [(r[3], r[2]) for r in rows if r[3] not in ("n", "NA")]
        self.assertEqual(named, [(f"chrX_{self.at + 5}_A_C", "S01")])

    def test_dict_less_tiers(self):
        self.check("tiers")

    def test_classic_dictionary(self):
        self.check("dict")

    def test_his_haplotype_is_his_calls_together(self):
        """Without an indel to take it off, the male's row is the haplotype he has -- both
        his SNPs -- and not every subset of them (which the enumerator gave an unphased
        call, and read a haploid one as)."""
        self.world(self.SITE, self.SNPS, (), "tiers", chrom="chrX", meta=self.META)
        rows = self.analyse(self.WINDOW, "0", 2)
        named = sorted((r[3], r[2]) for r in rows if r[3] not in ("n", "NA"))
        self.assertEqual(named, [
            (f"chrX_{self.at + 5}_A_C", "S01"),
            (f"chrX_{self.at + 5}_A_C,chrX_{self.at + 18}_T_A", "S00")])


class TestPerSampleIndelThatLeavesTheTarget(_PostAnalysis):
    """An indel inside the target that leaves its bases where they were takes no carrier
    off: the PAM's GG runs on into GGGG past the site, so a G lost there leaves the target's
    bases, read from its first base on (S00 keeps the row). A T lost in the TTTT the site
    holds does not (S01 loses it); S02 has the SNP alone."""

    SITE = "CTAAAAGTTGCTTTTATCACTGG"
    RIGHT = "GGGG" + "ACGTGCATTC" * 10
    SNPS = [(4, "A", "C", {"S00": "1|0", "S01": "1|0", "S02": "1|0"})]
    INDELS = [(21, "GG", "G", {"S00": "1|0"}), (12, "TT", "T", {"S01": "1|0"})]
    WINDOW = SITE[:4] + "M" + SITE[5:]

    def check(self, tiers):
        self.world(self.SITE, self.SNPS, self.INDELS, tiers, right=self.RIGHT)
        per_sample = self.analyse(self.WINDOW, "0", 0)
        self.assertEqual({r[:2] for r in per_sample}, {(self.GUIDE + "TGG", 0)})
        self.assertEqual(sorted(s for r in per_sample for s in r[2].split(",")),
                         ["S00", "S02"])

    def test_dict_less_tiers(self):
        self.check("tiers")

    def test_classic_dictionary(self):
        self.check("dict")


class TestPerSampleDenseWindowOnTheDictionary(_PostAnalysis):
    """A per-sample window denser than CRISPRME_IUPAC_CAP on an index without the dict-less
    tiers was reported by its greedy worst-case alignment, which took every allele that
    matches the guide and listed the carriers of any of them: here all three SNPs, 0
    mismatches, for S00, S01 and S02, none of whom carries more than two. Now its rows are
    the haplotypes its carriers have, as on the tiers."""

    # three mismatches with the guide (columns 2, 8 and 15), each restored by a SNP
    SNPS = [(2, "C", "A", {"S00": "1|0"}), (8, "A", "T", {"S00": "1|0", "S02": "0|1"}),
            (15, "T", "A", {"S01": "1|0"})]

    def test_rows_are_the_carriers_haplotypes(self):
        site = list(self.GUIDE + "TGG")
        for col, ref, _alt, _g in self.SNPS:
            site[col] = ref
        site = "".join(site)
        window = list(site)
        for col, ref, alt, _g in self.SNPS:
            window[col] = {frozenset("AC"): "M", frozenset("AT"): "W"}[frozenset(ref + alt)]
        self.world(site, self.SNPS, (), "dict")
        rows = self.analyse("".join(window), "0", 3, {"CRISPRME_IUPAC_CAP": "2"})
        snps = {r[3]: (r[1], r[2]) for r in rows if r[3] not in ("n", "NA")}
        p = lambda c: f"chrT_{self.at + c + 1}"  # noqa: E731
        self.assertEqual(snps, {
            f"{p(2)}_C_A,{p(8)}_A_T": (1, "S00"),
            f"{p(8)}_A_T": (2, "S02"),
            f"{p(15)}_T_A": (2, "S01"),
        })
        self.assertNotIn(0, [r[1] for r in rows])  # no row nobody carries


class TestGenotypeOracleDownload(unittest.TestCase):
    """download --what genotype-oracle: the phased HGDP + 1000 Genomes genotypes and indel
    genotypes into Dictionaries/, where the gate finds them, from the archives publish-index
    uploads of the panel's variant index (the HuggingFace client stubbed: the repository is a
    folder)."""

    def setUp(self):
        import crisprme_hf as hf
        self.hf = hf
        self.tmp = tempfile.mkdtemp(prefix="gate_oracle_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.hub = os.path.join(self.tmp, "hub")
        os.makedirs(self.hub)
        orig = hf._hf_snapshot
        self.addCleanup(setattr, hf, "_hf_snapshot", orig)
        self.fetched = []

        def snapshot(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
            for pat in allow_patterns:
                self.fetched.append(pat)
                f = os.path.join(self.hub, os.path.basename(pat))
                if os.path.isfile(f):
                    os.makedirs(os.path.join(local_dir, os.path.dirname(pat)), exist_ok=True)
                    shutil.copy(f, os.path.join(local_dir, pat))
            return local_dir
        hf._hf_snapshot = snapshot
        self.work = os.path.join(self.tmp, "work")
        os.makedirs(self.work)

    def publish(self, indels=True):
        """publish-index of the panel's variant index, built with its phased genotypes and
        (``indels``) its phased indel genotypes: the uploads are the repository's files.
        Returns their paths in the repository."""
        root = os.path.join(self.tmp, "src")
        idx = tdp._mklayout(root, f"NRG_3_hg38+{hp.ORACLE_PANEL}", registry=False,
                            genotypes=False)
        build_panel(os.path.join(root, "Dictionaries"), [(100, "A", "G", {"S00": "1|0"})],
                    [(104, "CTT", "C", {"S00": "1|0"})] if indels else (),
                    name=hp.ORACLE_PANEL)
        return [remote for _copy, remote in tdp._publish(idx, self.hub).uploads]

    def test_installs_both_stores_where_the_gate_reads_them(self):
        self.assertIn(f"indexes/indel_genotypes_{hp.ORACLE_PANEL}.tar.gz", self.publish())
        dests = self.hf.download_genotype_oracle(self.work, repo="local/test")
        dicts = os.path.join(self.work, "Dictionaries")
        self.assertEqual(dests, [os.path.join(dicts, f"{s}_{hp.ORACLE_PANEL}")
                                 for s in ("genotypes", "indel_genotypes")])
        self.assertEqual(hp.find_oracle(dicts), (hp.ORACLE_PANEL, dicts))
        p = hp.PanelChrom(hp.ORACLE_PANEL, dicts, CHROM, "oracle")
        self.addCleanup(p.close)
        self.assertEqual(p.carriers(v(100, "A", "G")), {"S00": "1|0"})
        self.assertEqual(p.carriers(v(104, "CTT", "C")), {"S00": "1|0"})
        # installed: kept, not fetched again
        self.fetched.clear()
        self.hf.download_genotype_oracle(self.work, repo="local/test")
        self.assertEqual(self.fetched, [])

    def test_without_indel_genotypes_on_the_hub_it_says_so(self):
        self.assertNotIn(f"indexes/indel_genotypes_{hp.ORACLE_PANEL}.tar.gz",
                         self.publish(indels=False))
        import contextlib
        import io
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            dests = self.hf.download_genotype_oracle(self.work, repo="local/test")
        self.assertEqual(len(dests), 1)
        self.assertIn("has no phased indel genotypes of hg38_HGDP1kGP", out.getvalue())

    def test_no_genotypes_on_the_hub_is_an_error(self):
        with self.assertRaises(ValueError):
            self.hf.download_genotype_oracle(self.work, repo="local/test")


class TestTheReportSaysWhyTheGateDidNotRun(unittest.TestCase):
    """The inputs card's Carrier gate row for each way the gate did not check the rows."""

    def test_each_status(self):
        import generate_report as gr
        said = {status: gr.render_gate_row({"status": status})
                for status in ("no-oracle", "opt-in", "off", "untiered", "per-sample",
                               "no-panel")}
        self.assertIn("crisprme.py download --what genotype-oracle", said["no-oracle"][1])
        self.assertIn("only on request", said["opt-in"][1])
        self.assertIn("Turned off", said["off"][1])
        self.assertIn("build-index-only --samplesID", said["untiered"][1])
        # a per-sample search's rows name their cis carriers; no phased store, nothing to say
        self.assertIsNone(said["per-sample"])
        self.assertIsNone(said["no-panel"])
        self.assertIsNone(gr.render_gate_row(None))
        # a companion the report could not read
        self.assertIn("not read", gr.render_gate_row({"status": "run", "used": False})[1])


class TheSearchSaysWhatTheGateDoes(bio._Fixture):
    """complete-search as a user runs it (a stub engine): which panel the gate reads, or why
    it does not run, and its flags checked before anything starts."""

    def test_it_names_the_panel_or_says_why_not(self):
        rc, out, err = self.run_bio()  # ds, with its phased genotype store
        self.assertEqual(rc, 0, err)
        rc, out, err = self.run_search("--population-level")
        self.assertIn("Launching job ", out, out + err)
        self.assertIn("[complete-search] Carrier gate: each row's variant combination is "
                      "checked against the phased genotypes of ds; the report ranks a "
                      "combination no haplotype carries after every other row", out)
        for options in (("--population-level", "--no-prioritize-observed"), ("--per-sample",)):
            with self.subTest(options=options):
                rc, out, err = self.run_search(*options)
                self.assertIn("Launching job ", out, out + err)
                self.assertNotIn("Carrier gate", out)
        # the same dataset with the classic dictionaries alone: said, not run
        shutil.rmtree(self.path("Dictionaries", "genotypes_ds"))
        shutil.rmtree(self.path("Dictionaries", "registry_ds"))
        rc, out, err = self.run_search("--population-level")
        self.assertIn("[complete-search] Carrier gate: not run, as the genotyped datasets of "
                      "this search have the classic per-sample dictionaries alone", out)

    def test_its_flags_are_checked_first(self):
        for options, said in (
                (("--prioritize-observed", "--no-prioritize-observed"),
                 "--prioritize-observed and --no-prioritize-observed are mutually exclusive"),
                (("--gate-min-expected", "-1"), "--gate-min-expected needs a number >= 0"),
                (("--gate-min-expected", "three"), "--gate-min-expected needs a number >= 0")):
            with self.subTest(options=options):
                rc, out, err = self.run_search(*options)
                self.assertEqual(rc, 1)
                self.assertIn(said, err)
                self.assertNotIn("Launching job", out)


class TestPlanPanels(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_plan_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        os.makedirs(os.path.join(self.dicts, "genotypes_hg38_G"))
        reg = os.path.join(self.dicts, "registry_hg38_mega")
        os.makedirs(reg)
        with open(os.path.join(reg, "reg_chr22.idx"), "w") as fh:
            json.dump({"aggregation": "info_af", "data_type": "sites-only"}, fh)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_modes(self):
        d = self.dicts
        self.assertEqual(hp.plan_panels(d, ["hg38_G"], "population-level", "auto"),
                         ("run", [("hg38_G", d, "own")]))
        self.assertEqual(hp.plan_panels(d, ["hg38_G"], "per-sample", "auto"),
                         ("per-sample", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_G"], "population-level", "off"), ("off", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_mega"], "population-level", "auto"),
                         ("opt-in", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_mega"], "population-level", "on"),
                         ("no-oracle", []))
        # a genotyped dataset with the classic dictionaries alone: said, not run
        os.makedirs(os.path.join(d, "dictionaries_hg38_1000G"))
        self.assertEqual(hp.plan_panels(d, ["hg38_1000G"], "population-level", "auto"),
                         ("untiered", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_1000G"], "per-sample", "auto"),
                         ("no-panel", []))

    def test_the_oracle_conventional_or_given(self):
        os.makedirs(os.path.join(self.dicts, f"genotypes_{hp.ORACLE_PANEL}"))
        self.assertEqual(hp.plan_panels(self.dicts, ["hg38_mega"], "population-level", "on"),
                         ("run", [(hp.ORACLE_PANEL, self.dicts, "oracle")]))
        other = os.path.join(self.tmp, "elsewhere")
        os.makedirs(os.path.join(other, "genotypes_hg38_HGDP1kGP_v2"))
        given = os.path.join(other, "genotypes_hg38_HGDP1kGP_v2")
        self.assertEqual(hp.find_oracle(self.dicts, given), ("hg38_HGDP1kGP_v2", other))
        # a genotyped dataset searched beside the mega is read once, as its own panel
        self.assertEqual(
            hp.plan_panels(self.dicts, [hp.ORACLE_PANEL, "hg38_mega"], "population-level", "on"),
            ("run", [(hp.ORACLE_PANEL, self.dicts, "own")]))


if __name__ == "__main__":
    unittest.main()
