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
import site_alleles
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


def _enumerate(positions, max_putative=-1, **caps):
    return oh.enumerate_observed_haplotypes(
        positions, REF, parse_haplotypes, _ploidy_diploid, "chrT",
        max_putative=max_putative, **caps
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

    def test_each_alt_read_against_its_own_index_among_the_position_s_alts(self):
        # The resolver new_simple_analysis hands the enumerator (site_alleles): an alt's
        # index is the token its carriers share, narrowed by what the position's other
        # alts need. altA's two carriers both hold 1 and 2 (1|2, 2|1): the most frequent
        # token of altA's own carriers (_alt_index_for_record, removed) ties and took 2,
        # putting altA on altB's slots; altB's 0|2 carrier says 2 is altB's.
        site = [("G", {"X": "1|2", "Y": "2|1"}),
                ("T", {"X": "1|2", "Y": "2|1", "Z": "0|2"})]
        self.assertEqual(site_alleles.alt_tokens(site), ["1", "2"])
        got = dict(site_alleles.per_alt_entries(site))
        self.assertEqual(got["G"], {"X": "1|0", "Y": "0|1"})
        self.assertEqual(got["T"], {"X": "0|1", "Y": "1|0", "Z": "0|1"})
        positions = [_col(0, ["G", "T"], [got["G"], got["T"]])]
        haps = _by_set(_enumerate(positions))
        self.assertEqual(haps[frozenset({(0, "G")})].carriers, {"X", "Y"})
        self.assertEqual(haps[frozenset({(0, "T")})].carriers, {"X", "Y", "Z"})
        for h in haps.values():
            self.assertEqual(h.phase_state, oh.CONFIRMED)


class TwoAltsOfOnePosition(unittest.TestCase):
    """CRISPRitz 2.10.0 codes every ALT of a split site, and the panel's split records
    were phased one by one: NA18869 is 0|1 under G and 0|1 under T at chr22:10,661,629
    (672 such sample-positions on chr22 of HGDP + 1000 Genomes). Up to 2.8.1 the slot
    became one CONFIRMED set holding both alts: a row naming two alleles at one position,
    as confirmed cis, and the sample a carrier of neither single-alt set."""

    def _row(self, g_gt, t_gt, extra=None):
        cols = [_col(1, ["G", "T"], [{"S": g_gt, "O": "1|0"}, {"S": t_gt}])]
        if extra is not None:
            cols.append(_col(3, ["C"], [extra]))
        return _by_set(_enumerate(cols))

    def _no_set_with_two_alts_of_a_column(self, haps):
        for vset in haps:
            cols = [pos_c for pos_c, _alt in vset]
            self.assertEqual(len(cols), len(set(cols)), sorted(vset))

    def test_two_alts_on_one_slot_make_the_sample_putative_under_each(self):
        for g_gt, t_gt in (("0|1", "0|1"), ("1|0", "1|0"), ("1|1", "0|1")):
            with self.subTest(g=g_gt, t=t_gt):
                haps = self._row(g_gt, t_gt)
                self._no_set_with_two_alts_of_a_column(haps)
                self.assertIn("S", haps[frozenset({(1, "G")})].carriers)
                self.assertIn("S", haps[frozenset({(1, "T")})].carriers)
                # {G} is shared with O (1|0, CONFIRMED): the label is per set
                self.assertEqual(haps[frozenset({(1, "G")})].phase_state, oh.PUTATIVE)
                self.assertEqual(haps[frozenset({(1, "T")})].phase_state, oh.PUTATIVE)

    def test_three_alts_on_one_slot(self):
        cols = [_col(1, ["C", "G", "T"], [{"S": "1|0"}, {"S": "1|0"}, {"S": "1|0"}])]
        haps = _by_set(_enumerate(cols))
        self.assertEqual(set(haps), {frozenset({(1, a)}) for a in "CGT"})
        for h in haps.values():
            self.assertEqual((h.carriers, h.phase_state), ({"S"}, oh.PUTATIVE))

    def test_a_second_column_on_the_conflicting_slot(self):
        haps = self._row("0|1", "0|1", extra={"S": "0|1"})
        self._no_set_with_two_alts_of_a_column(haps)
        for vset in ({(1, "G"), (3, "C")}, {(1, "T"), (3, "C")}, {(3, "C")}):
            self.assertIn("S", haps[frozenset(vset)].carriers)
            self.assertEqual(haps[frozenset(vset)].phase_state, oh.PUTATIVE)

    def test_alts_on_different_slots_stay_confirmed(self):
        haps = self._row("1|0", "0|1", extra={"S": "0|1"})
        self.assertEqual(set(haps), {frozenset({(1, "G")}),
                                     frozenset({(1, "T"), (3, "C")})})
        for h in haps.values():
            self.assertEqual(h.phase_state, oh.CONFIRMED)
        self.assertEqual(haps[frozenset({(1, "G")})].carriers, {"S", "O"})

    def test_an_unphased_heterozygote_of_two_alts_takes_one_alt_a_set(self):
        # a 1/2 heterozygote, split into 1/0 under A and 0/1 under C
        cols = [_col(3, ["A", "C"], [{"U": "1/0"}, {"U": "0/1"}]),
                _col(4, ["G"], [{"U": "0/1"}])]
        haps = _by_set(_enumerate(cols))
        self.assertEqual(set(haps), {
            frozenset({(3, "A")}), frozenset({(3, "C")}), frozenset({(4, "G")}),
            frozenset({(3, "A"), (4, "G")}), frozenset({(3, "C"), (4, "G")})})

    def test_the_union_above_the_cap_is_one_per_choice_of_alt(self):
        cols = [_col(1, ["A", "C"], [{"U": "1/0"}, {"U": "0/1"}]),
                _col(2, ["G"], [{"U": "0/1"}]),
                _col(3, ["G"], [{"U": "0/1"}])]
        haps = _by_set(_enumerate(cols, max_putative=2))
        self.assertEqual(set(haps), {
            frozenset({(1, "A"), (2, "G"), (3, "G")}),
            frozenset({(1, "C"), (2, "G"), (3, "G")})})

    def test_the_unions_above_the_cap_are_at_most_max_unions(self):
        """A sample PUTATIVE in a window (here a phased slot holding two alts at column
        0) that carries two alts at k columns: one union per choice of alt made 2^k of
        them (NA19042 carries two or three at 16 positions of one 26-bp IGL window of
        chr22: 8.5 million sets and 90 GB). They stay while there are at most max_unions (16);
        past that there is one, with each column's first alt, or the one ``choose``
        picks (new_simple_analysis: the alt nearest the guide)."""
        def cols(k):
            return [_col(c, ["A", "T"], [{"S": "0|1"}, {"S": "1|0" if c else "0|1"}])
                    for c in range(k)]
        haps = _enumerate(cols(4), max_putative=2)
        self.assertEqual(len(haps), 16)
        self.assertEqual({h.phase_state for h in haps}, {oh.PUTATIVE})
        haps = _enumerate(cols(14), max_putative=2)
        self.assertEqual([h.variant_set for h in haps],
                         [frozenset((c, "A") for c in range(14))])
        self.assertEqual((haps[0].carriers, haps[0].phase_state), ({"S"}, oh.PUTATIVE))
        haps = _enumerate(cols(14), max_putative=2, choose=lambda c, alts: alts[c % 2])
        self.assertEqual([h.variant_set for h in haps],
                         [frozenset((c, "AT"[c % 2]) for c in range(14))])
        haps = _enumerate(cols(5), max_putative=2, max_unions=32)
        self.assertEqual(len(haps), 32)


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
