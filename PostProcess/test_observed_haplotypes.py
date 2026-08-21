"""Unit tests for the observed-haplotype enumerator (dict-less multi-variant path).

Covers the LOCKED-DESIGN test matrix at the PURE-FUNCTION level (no genome / dict /
arg files, no numpy / pysam): phased cis, phased trans, unphased, mixed panel, phase-
set block, chimera exclusion, dense region, multiallelic alt-index resolution, and the
carrier/phase-state dedup policy.

STDLIB ONLY. ``observed_haplotypes`` takes ``parse_haplotypes`` as an argument, so we
import the real ``population_summary.parse_haplotypes`` (itself stdlib-only via
tier0_registry) as the single source of allele-slot truth -- no MagicMock needed for
the slot math. A trivial ``ploidy_of`` stub (everyone diploid) stands in for the
chromosome-ploidy model.

Run with:
    cd PostProcess && python3 -m unittest test_observed_haplotypes -v
"""

import os
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

import observed_haplotypes as oh
from population_summary import parse_haplotypes


def _ploidy_diploid(sid, sex):
    return 2


def _col(pos_c, alts, carrier_gts, alt_index=None, info=None, ps=None):
    """Build one column descriptor. ``carrier_gts`` is a list (per alt) of
    {sample: gt}; ``alt_index`` defaults to "1" per alt; ``info`` defaults to a
    placeholder [rsID, AF, snp_info] per alt."""
    n = len(alts)
    if alt_index is None:
        alt_index = ["1"] * n
    if info is None:
        info = [["rs%d_%d" % (pos_c, i), "0.10", "chrT_%d_C_%s" % (pos_c + 1, alts[i])]
                for i in range(n)]
    return {
        "pos_c": pos_c,
        "alts": list(alts),
        "carrier_gts": list(carrier_gts),
        "alt_index": list(alt_index),
        "info": list(info),
        "ps": ps,
    }


REF = "AAAAA"  # 5-nt reference window; ambiguity columns applied per test


def _enumerate(positions, max_putative=-1):
    return oh.enumerate_observed_haplotypes(
        positions, REF, parse_haplotypes, _ploidy_diploid, "chrT",
        max_putative=max_putative,
    )


def _by_set(haps):
    """Map frozenset[(pos,alt)] -> ObservedHaplotype for easy assertions."""
    return {h.variant_set: h for h in haps}


class PhasedCis(unittest.TestCase):
    def test_phased_cis_single_confirmed_haplotype(self):
        # One sample S carries alt on hap0 at columns 0,1,2 (all "1|0") -> ONE
        # CONFIRMED haplotype carrying all three (pos,alt).
        positions = [
            _col(0, ["G"], [{"S": "1|0"}]),
            _col(1, ["G"], [{"S": "1|0"}]),
            _col(2, ["G"], [{"S": "1|0"}]),
        ]
        haps = _enumerate(positions)
        self.assertEqual(len(haps), 1)
        h = haps[0]
        self.assertEqual(h.variant_set,
                         frozenset({(0, "G"), (1, "G"), (2, "G")}))
        self.assertEqual(h.carriers, {"S"})
        self.assertEqual(h.phase_state, oh.CONFIRMED)
        self.assertEqual("".join(h.seq), "GGGAA")


class PhasedTrans(unittest.TestCase):
    def test_phased_trans_two_confirmed_singletons(self):
        # Alts on ALTERNATING slots (col0 on hap0, col1 on hap1) -> TWO CONFIRMED
        # single-variant haplotypes, NOT a combined cis haplotype (no chimera).
        positions = [
            _col(0, ["G"], [{"S": "1|0"}]),
            _col(1, ["G"], [{"S": "0|1"}]),
        ]
        haps = _by_set(_enumerate(positions))
        self.assertEqual(len(haps), 2)
        self.assertIn(frozenset({(0, "G")}), haps)
        self.assertIn(frozenset({(1, "G")}), haps)
        self.assertNotIn(frozenset({(0, "G"), (1, "G")}), haps)
        for h in haps.values():
            self.assertEqual(h.phase_state, oh.CONFIRMED)
            self.assertEqual(h.carriers, {"S"})


class Unphased(unittest.TestCase):
    def test_unphased_enumerates_subcombinations(self):
        # 0/1 at two positions, UNPHASED -> the cis/trans phasing is unknown, so emit
        # EVERY non-empty subset as a candidate cis haplotype (a variant that would BREAK
        # a target -- e.g. disrupt the PAM -- must be droppable, so the maximal union
        # cannot be the only set). All PUTATIVE, sample present in each.
        positions = [
            _col(0, ["G"], [{"S": "0/1"}]),
            _col(1, ["G"], [{"S": "0/1"}]),
        ]
        haps = _by_set(_enumerate(positions))
        self.assertEqual(
            set(haps),
            {
                frozenset({(0, "G")}),
                frozenset({(1, "G")}),
                frozenset({(0, "G"), (1, "G")}),
            },
        )
        for h in haps.values():
            self.assertEqual(h.phase_state, oh.PUTATIVE)
            self.assertEqual(h.carriers, {"S"})


class MixedPanel(unittest.TestCase):
    def test_mixed_phased_unphased_trans(self):
        # A phased-cis (both cols hap0), B unphased, C phased-trans.
        positions = [
            _col(0, ["G"], [{"A": "1|0", "B": "0/1", "C": "1|0"}]),
            _col(1, ["G"], [{"A": "1|0", "B": "0/1", "C": "0|1"}]),
        ]
        haps = _by_set(_enumerate(positions))
        # A cis set {(0,G),(1,G)} CONFIRMED; B (unphased) now contributes EVERY subset
        # {0},{1},{0,1} as PUTATIVE. The full set coincides with A's cis -> downgraded to
        # PUTATIVE, carriers {A,B}.
        cis = haps[frozenset({(0, "G"), (1, "G")})]
        self.assertEqual(cis.carriers, {"A", "B"})
        self.assertEqual(cis.phase_state, oh.PUTATIVE)  # downgraded by B's union subset
        # C contributes CONFIRMED singletons {0} and {1}; B contributes the SAME sets as
        # PUTATIVE -> shared -> downgraded to PUTATIVE with carriers {B,C}.
        s0 = haps[frozenset({(0, "G")})]
        s1 = haps[frozenset({(1, "G")})]
        self.assertEqual(s0.carriers, {"B", "C"})
        self.assertEqual(s0.phase_state, oh.PUTATIVE)
        self.assertEqual(s1.carriers, {"B", "C"})
        self.assertEqual(s1.phase_state, oh.PUTATIVE)


class PhaseSetBlock(unittest.TestCase):
    def test_different_ps_ids_enumerates_subcombinations(self):
        # Two phased calls in DIFFERENT phase sets -> cannot prove cis -> treated as
        # unphased -> enumerate every non-empty subset, all PUTATIVE.
        positions = [
            _col(0, ["G"], [{"S": "1|0"}], ps=[{"S": 100}]),
            _col(1, ["G"], [{"S": "1|0"}], ps=[{"S": 200}]),
        ]
        haps = _by_set(_enumerate(positions))
        self.assertEqual(
            set(haps),
            {
                frozenset({(0, "G")}),
                frozenset({(1, "G")}),
                frozenset({(0, "G"), (1, "G")}),
            },
        )
        for h in haps.values():
            self.assertEqual(h.phase_state, oh.PUTATIVE)

    def test_unphased_dense_window_falls_back_to_union(self):
        # GUARD: a sample carrying more than max_putative variants in the window falls
        # back to the union only (the 2^k blow-up is confined to that sample).
        positions = [_col(i, ["G"], [{"S": "0/1"}]) for i in range(4)]
        haps = _enumerate(positions, max_putative=2)  # 4 > 2 -> union only
        self.assertEqual(len(haps), 1)
        self.assertEqual(
            haps[0].variant_set,
            frozenset({(0, "G"), (1, "G"), (2, "G"), (3, "G")}),
        )

    def test_same_ps_id_is_confirmed(self):
        positions = [
            _col(0, ["G"], [{"S": "1|0"}], ps=[{"S": 100}]),
            _col(1, ["G"], [{"S": "1|0"}], ps=[{"S": 100}]),
        ]
        haps = _enumerate(positions)
        self.assertEqual(len(haps), 1)
        self.assertEqual(haps[0].variant_set, frozenset({(0, "G"), (1, "G")}))
        self.assertEqual(haps[0].phase_state, oh.CONFIRMED)

    def test_no_ps_tag_single_block_confirmed(self):
        # PS absent (ps=None) => single whole-chromosome block => cis => CONFIRMED.
        positions = [
            _col(0, ["G"], [{"S": "1|0"}]),
            _col(1, ["G"], [{"S": "1|0"}]),
        ]
        haps = _enumerate(positions)
        self.assertEqual(haps[0].phase_state, oh.CONFIRMED)


class ChimeraExclusion(unittest.TestCase):
    def test_no_cross_individual_chimera(self):
        # A carries only col0, B carries only col1; NO single sample carries both ->
        # the combined {(0,G),(1,G)} haplotype must NEVER be emitted.
        positions = [
            _col(0, ["G"], [{"A": "1|0"}]),
            _col(1, ["G"], [{"B": "1|0"}]),
        ]
        haps = _by_set(_enumerate(positions))
        self.assertIn(frozenset({(0, "G")}), haps)
        self.assertIn(frozenset({(1, "G")}), haps)
        self.assertNotIn(frozenset({(0, "G"), (1, "G")}), haps)
        self.assertEqual(haps[frozenset({(0, "G")})].carriers, {"A"})
        self.assertEqual(haps[frozenset({(1, "G")})].carriers, {"B"})


class DenseRegion(unittest.TestCase):
    def test_dense_window_only_real_haplotypes(self):
        # k=5 columns but only 3 distinct real haplotypes exist -> exactly 3 rows
        # (NOT 2^5=32, NOT one greedy row). Build 5 columns; three samples each with a
        # distinct cis haplotype.
        # S1: cols 0,1,2 on hap0 ; S2: cols 2,3 on hap0 ; S3: col 4 on hap0.
        cols = {i: {} for i in range(5)}
        for i in (0, 1, 2):
            cols[i]["S1"] = "1|0"
        for i in (2, 3):
            cols[i]["S2"] = "1|0"
        cols[4]["S3"] = "1|0"
        positions = [_col(i, ["G"], [cols[i]]) for i in range(5)]
        haps = _by_set(_enumerate(positions))
        self.assertEqual(len(haps), 3)
        self.assertIn(frozenset({(0, "G"), (1, "G"), (2, "G")}), haps)
        self.assertIn(frozenset({(2, "G"), (3, "G")}), haps)
        self.assertIn(frozenset({(4, "G")}), haps)


class MultiallelicAltIndex(unittest.TestCase):
    def test_altB_foreign_index_not_dropped(self):
        # A genuine multiallelic column: altA (token "1") and altB (token "2") at the
        # SAME position. Carrier X of altB has genotype "1|2" (altA on hap0, altB on
        # hap1). Without per-record alt-index resolution, altB (tested against "1")
        # would be DROPPED. With it, both alts are enumerated on the right slots.
        col = _col(
            0,
            ["G", "T"],                     # altA=G (index 1), altB=T (index 2)
            [{"X": "1|2"}, {"X": "1|2"}],   # X is a carrier of BOTH records
            alt_index=["1", "2"],
        )
        positions = [col]
        haps = _by_set(_enumerate(positions))
        # X: hap0 carries altA(G) (token 1), hap1 carries altB(T) (token 2) -> two
        # CONFIRMED singletons, altB NOT dropped.
        self.assertIn(frozenset({(0, "G")}), haps)
        self.assertIn(frozenset({(0, "T")}), haps)
        self.assertEqual(haps[frozenset({(0, "G")})].carriers, {"X"})
        self.assertEqual(haps[frozenset({(0, "T")})].carriers, {"X"})

    def test_alt_index_resolution_helper(self):
        # The resolver picks the shared non-ref token of a record's carriers.
        # Biallelic default -> "1".
        self.assertEqual(
            oh._alt_index_for_record(["1|0", "0|1", "0/1"], parse_haplotypes, 2), "1")
        # A record whose carriers all show "2" (altB of a multiallelic) -> "2".
        self.assertEqual(
            oh._alt_index_for_record(["1|2", "2|0", "0/2"], parse_haplotypes, 2), "2")


class HomozygousAndBoundary(unittest.TestCase):
    def test_phased_hom_single_confirmed_set(self):
        # "1|1" at two cols => both slots identical => ONE CONFIRMED set (not two).
        positions = [
            _col(0, ["G"], [{"S": "1|1"}]),
            _col(1, ["G"], [{"S": "1|1"}]),
        ]
        haps = _enumerate(positions)
        self.assertEqual(len(haps), 1)
        self.assertEqual(haps[0].variant_set, frozenset({(0, "G"), (1, "G")}))
        self.assertEqual(haps[0].phase_state, oh.CONFIRMED)

    def test_boundary_pos_c_past_refseq_no_indexerror(self):
        # An ambiguity column index past len(REF) must be applied defensively (seq
        # untouched) without IndexError. REF len is 5; use pos_c=9.
        positions = [_col(9, ["G"], [{"S": "1|0"}])]
        haps = _enumerate(positions)  # must not raise
        self.assertEqual(len(haps), 1)
        self.assertEqual("".join(haps[0].seq), REF)  # unchanged (out of range)


if __name__ == "__main__":
    unittest.main()
