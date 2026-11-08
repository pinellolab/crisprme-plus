#!/usr/bin/env python
"""Per-carrier-haplotype resolution of SNP-overlaid indel targets (2.8.1).

The SNP overlay of the fake-indel genome writes each enriched-genome IUPAC code onto
the flanks of a fake-indel window, and the engine matches a code when ANY of its
alleles fits, so one raw indel target stands for every combination of the window's
SNP alleles. Up to 2.8.0 the indel post-analysis turned each code into one base
column by column (the guide's base when it was an allele, else the code's first
allele) without any genotype, which produced:

  * chimeras: one sample's indel reported with SNP alleles only other samples carry
    (chr18:12386282(-), CFD 0.708, "1 cis carrier" NA19072, who is 0|0 at the SNP
    whose only carrier is NA20541), and
  * PAMs the search never matched: an R at NRG's last base always read A.

``indel_snp_cis.resolve_target`` replaces that: every indel carrier slot gets its
own allele at each code column. These tests pin

  (a) a property over every IUPAC code at every PAM position, for NRG, NGG, NRR and
      a 5' TTTV PAM: whenever the raw PAM can satisfy the motif (pam_filter keeps
      the target), every resolved row's PAM satisfies it, with and without
      genotypes, and is concrete;
  (b) table-driven per-haplotype cases on a small phased fixture: a SNP in cis / in
      trans / absent / unphased / with a missing allele / with no genotype source,
      the chr18 geometry, a minus-strand PAM SNP, a bulge next to the PAM, split
      multi-allelic records and several alts on one VCF line, overlapping indels,
      haploid chrX calls, the edit budget, the worst case against a brute force, the
      PAM-window SNPs of the co-occurrence decode, and which targets of a fake-indel
      window hold its indel at all;
  (c) that analisi_indels_NNN.py routes every overlaid target through it, with the
      search budget plumbed in.

Each fails on 2.8.0 / aa9c54b (no resolve_target; the column-by-column resolver).
Hermetic, stdlib only (+ the repo's CFD tables).

Run: python -m unittest test_indel_haplotype_resolution -v
"""

import itertools
import os
import random
import unittest

import indel_snp_cis as isc
import pam_filter
from bulge_cfd_penalty import cfd_score_bulge_aware

HERE = os.path.dirname(os.path.abspath(__file__))
GUIDE = "ACTGAAATCTGTAAGCAGGC"
CODES = "ACGTRYSWKMBDHVN"  # the 15 IUPAC characters
MOTIF_BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
               "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
               "H": "ACT", "V": "ACG", "N": "ACGT"}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _cfd(guide, seq):
    return cfd_score_bulge_aware(guide, seq)


def _col(col, ref_fwd, alts=None, code=None, strand="+", known=True, real=None):
    """A code column: forward ref, {forward alt: {sample: gt}}, the code shown."""
    if code is None:
        shown = {ref_fwd} | set(alts or {})
        if strand == "-":
            shown = {COMP[b] for b in shown}
        code = next(c for c, b in MOTIF_BASES.items()
                    if set(b) == shown and c != "N") if len(shown) > 1 else shown.pop()
    return isc.SnpColumn(col, 1000 + col if real is None else real, ref_fwd, alts,
                         code, strand, known=known)


def _resolve(guide, target, cols, gts, motif="NRG", pam_at_start=False, max_mm=4,
             max_total=-1, bulges=0, phase_known=True, others=(), scored=True):
    motif_at = isc.pam_columns(len(target), motif, pam_at_start)
    n = len(motif)
    spacer = range(n, len(target)) if pam_at_start else range(len(target) - n)
    if scored and not pam_at_start and len(motif) == 3:
        def score(bases):
            return _cfd(guide, isc.realize(target, guide, cols, bases))
    else:
        score = None
    return isc.resolve_target(guide, target, cols, gts, motif_at, spacer, max_mm,
                              max_total, bulges, phase_known=phase_known,
                              score_fn=score, others=others)


# --------------------------------------------------------------------------- #
# (a) property: a resolved PAM always satisfies the motif when the raw PAM can
# --------------------------------------------------------------------------- #
class ResolvedPamSatisfiesTheMotif(unittest.TestCase):
    MOTIFS = (("NRG", False), ("NGG", False), ("NRR", False), ("TTTV", True))

    def _targets(self, motif, pam_at_start):
        """Every raw PAM over the 15 IUPAC characters at every position (both cases)
        that pam_filter would keep, on a perfect-match spacer."""
        n = len(motif)
        for pam in itertools.product(CODES, repeat=n):
            pam = "".join(pam)
            if not pam_filter._pam_can_match(pam, motif):
                continue
            for case in (str.upper, str.lower):
                p = "".join(case(c) if c in "RYSWKMBDHV" else c for c in pam)
                if pam_at_start:
                    yield "N" * n + GUIDE, p + GUIDE, list(range(n))
                else:
                    yield GUIDE + "N" * n, GUIDE + p, list(range(len(GUIDE), len(GUIDE) + n))

    def _check_rows(self, rows, motif, pam_at_start, target):
        n = len(motif)
        for r in rows:
            pam = r["seq"][:n] if pam_at_start else r["seq"][-n:]
            self.assertTrue(set(pam) <= set("ACGT"), (target, r["seq"]))
            self.assertTrue(pam_filter._pam_can_match(pam, motif), (target, r["seq"]))

    def test_without_genotypes_the_worst_case_keeps_the_pam(self):
        """A sites-only panel: one PUTATIVE row, the worst case over the codes."""
        for motif, start in self.MOTIFS:
            n_checked = 0
            for guide, target, pcols in self._targets(motif, start):
                amb = isc.ambiguous_columns(target)
                if not amb:
                    continue
                cols = []
                for j in amb:
                    bases = MOTIF_BASES[target[j].upper()]
                    cols.append(_col(j, bases[0], None, target[j], known=False))
                rows = _resolve(guide, target, cols, None, motif, start)
                self.assertEqual(len(rows), 1, (motif, target))
                self.assertEqual(rows[0]["phase"], isc.PUTATIVE)
                self._check_rows(rows, motif, start, target)
                n_checked += 1
            self.assertGreater(n_checked, 100, motif)

    def test_each_carrier_allele_either_keeps_the_pam_or_drops_the_row(self):
        """A genotyped panel: every carrier gets its own allele; a row survives only
        with a PAM the motif allows, and every allele that allows it gives one."""
        for motif, start in self.MOTIFS:
            for guide, target, pcols in self._targets(motif, start):
                amb = isc.ambiguous_columns(target)
                if len(amb) != 1:
                    continue  # one code column: each (ref, alt) split is exhaustive
                j = amb[0]
                bases = MOTIF_BASES[target[j].upper()]
                motif_char = motif[pcols.index(j)]
                for ref, alt in itertools.permutations(bases, 2):
                    cols = [_col(j, ref, {alt: {"ALT": "1|0"}}, target[j])]
                    rows = _resolve(guide, target, cols,
                                    {"ALT": "1|0", "REF": "1|0"}, motif, start)
                    self._check_rows(rows, motif, start, target)
                    got = {r["seq"][j]: r["samples"] for r in rows}
                    want = {b: s for b, s in ((ref, ["REF"]), (alt, ["ALT"]))
                            if b in MOTIF_BASES[motif_char]}
                    self.assertEqual(got, want, (motif, target, ref, alt))

    def test_the_observed_nrg_case(self):
        """The audit's 1,088 rows: an R at NRG's last base (ref G, alt A). The base
        2.8.0 showed (A) is a PAM the search never matched."""
        guide, target = GUIDE + "NNN", GUIDE + "TGR"
        cols = [_col(22, "G", {"A": {"S2": "1|0", "S4": "1|0"}}, "R")]
        rows = _resolve(guide, target, cols, {"S1": "1|0", "S2": "1|0", "S3": "0|1"})
        self.assertEqual([(r["seq"], r["samples"]) for r in rows],
                         [(GUIDE + "TGG", ["S1", "S3"])])
        self.assertAlmostEqual(_cfd(guide, rows[0]["seq"]), 1.0)


# --------------------------------------------------------------------------- #
# (b) table-driven per-haplotype cases
# --------------------------------------------------------------------------- #
class PerHaplotypeCases(unittest.TestCase):
    # chr18:12386282(-): raw hit A-TaAAAcaWGcAAGCAGGCAGG, 4mm+1b; the W is the T>A
    # SNP chr18:12386296 (forward), shown A>T on the minus strand; the guide wants T.
    CHR18_GUIDE = "ACTGAAATCTGTAAGCAGGCNNN"
    CHR18_DNA = "A-TaAAAcaWGcAAGCAGGCAGG"

    def _chr18(self, carriers):
        return [_col(9, "T", {"A": carriers}, "W", strand="-", real=12386295)]

    def test_chr18_the_trans_carrier_chimera_is_dropped(self):
        rows = _resolve(self.CHR18_GUIDE, self.CHR18_DNA,
                        self._chr18({"NA20541": "0|1"}), {"NA19072": "0|1"}, bulges=1)
        self.assertEqual(rows, [], "NA19072 is 0|0 at the SNP: 5mm+1b, past --mm 4")

    def test_chr18_the_cis_haplotype_is_kept_and_named(self):
        cols = self._chr18({"NA19072": "0|1"})
        rows = _resolve(self.CHR18_GUIDE, self.CHR18_DNA, cols, {"NA19072": "0|1"},
                        bulges=1)
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual((r["seq"], r["mm"], r["samples"], r["phase"], r["slots"]),
                         ("A-TaAAAcaTGcAAGCAGGCAGG", 4, ["NA19072"], isc.CONFIRMED, 1))
        self.assertAlmostEqual(_cfd(self.CHR18_GUIDE, r["seq"]), 0.708, places=3)
        self.assertEqual(r["enabling"], [0])  # the SNP the target requires

    def test_chr18_same_sample_trans_is_dropped(self):
        rows = _resolve(self.CHR18_GUIDE, self.CHR18_DNA,
                        self._chr18({"NA19072": "1|0"}), {"NA19072": "0|1"}, bulges=1)
        self.assertEqual(rows, [])

    def test_cis_trans_absent_split_into_their_own_rows(self):
        """One spacer SNP the guide matches only with the alt (C>A at column 5): the
        cis slot has a perfect target, the trans and absent ones a 1-mismatch one."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"CIS": "1|0", "TRANS": "0|1", "OTHER": "1|1"}}, "M")]
        rows = _resolve(guide, target, cols,
                        {"CIS": "1|0", "TRANS": "1|0", "ABSENT": "0|1"})
        got = sorted((r["mm"], r["seq"], tuple(r["samples"]), r["slots"]) for r in rows)
        self.assertEqual(got, [
            (0, GUIDE + "TGG", ("CIS",), 1),
            (1, GUIDE[:5] + "c" + GUIDE[6:] + "TGG", ("ABSENT", "TRANS"), 2),
        ])
        self.assertTrue(all(r["phase"] == isc.CONFIRMED for r in rows))

    def test_a_homozygous_indel_gives_one_slot_per_haplotype(self):
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"HOM": "0|1"}}, "M")]
        rows = _resolve(guide, target, cols, {"HOM": "1|1"})
        self.assertEqual(sorted((r["mm"], r["samples"], r["slots"]) for r in rows),
                         [(0, ["HOM"], 1), (1, ["HOM"], 1)])

    def test_unphased_is_putative_and_lists_only_co_carriers(self):
        """With '/' genotypes the slot is unknown: the worst case its OWN alleles
        allow, PUTATIVE; a sample that lacks the alt never gets the alt's row."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"HET": "0/1", "HOMALT": "1/1"}}, "M")]
        rows = _resolve(guide, target, cols,
                        {"HET": "0/1", "HOMALT": "0/1", "NOSNP": "0/1"})
        got = sorted((r["mm"], tuple(r["samples"]), r["phase"]) for r in rows)
        self.assertEqual(got, [
            (0, ("HET", "HOMALT"), isc.PUTATIVE),   # the alt's row: only its carriers
            (1, ("NOSNP",), isc.CONFIRMED),         # no SNP: the same on any slot
        ])

    def test_an_unphased_homozygous_alt_has_no_reference_allele(self):
        """1/1 at a SNP that spoils the target: no slot of this carrier holds the
        reference base, so it is never listed for the perfect target. The reference
        is its candidate only with fewer alt copies than its ploidy."""
        guide = GUIDE + "NNN"
        target = GUIDE[:3] + "K" + GUIDE[4:] + "TGG"  # ref G (the guide's), alt T
        cols = [_col(3, "G", {"T": {"HOM": "1/1", "HET": "0/1"}}, "K")]
        rows = _resolve(guide, target, cols, {"HOM": "0/1"})
        self.assertEqual([(r["seq"][3], r["mm"], r["samples"], r["phase"]) for r in rows],
                         [("t", 1, ["HOM"], isc.PUTATIVE)])
        rows = _resolve(guide, target, cols, {"HET": "0/1"})
        self.assertEqual([(r["seq"][3], r["mm"], r["samples"], r["phase"]) for r in rows],
                         [("G", 0, ["HET"], isc.PUTATIVE)])

    def test_a_missing_allele_on_the_indel_slot_is_putative(self):
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"S": ".|1"}}, "M")]
        rows = _resolve(guide, target, cols, {"S": "1|0"})
        self.assertEqual([(r["mm"], r["samples"], r["phase"]) for r in rows],
                         [(0, ["S"], isc.PUTATIVE)])

    def test_no_phased_indel_genotype_uses_the_log_carriers_putatively(self):
        """The indel store lacks the indel: its log carriers, phase unknown."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"A1": "1|0"}}, "M")]
        carriers = isc.IndelCarriers({"A1": "0/1", "B1": "0/1"}, phase_known=False)
        rows = _resolve(guide, target, cols, carriers)
        got = sorted((r["mm"], tuple(r["samples"]), r["phase"]) for r in rows)
        self.assertEqual(got, [(0, ("A1",), isc.PUTATIVE), (1, ("B1",), isc.CONFIRMED)])

    def test_no_genotype_source_keeps_the_log_carriers(self):
        """Sites-only: one PUTATIVE worst-case row; samples None (the log's)."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGR"
        cols = [_col(5, "C", None, "M", known=False),
                _col(22, "G", None, "R", known=False)]
        rows = _resolve(guide, target, cols, None)
        self.assertEqual([(r["seq"], r["samples"], r["phase"]) for r in rows],
                         [(GUIDE + "TGG", None, isc.PUTATIVE)])

    def test_reference_only_reading(self):
        """CRISPRME_INDEL_SNP=0: each code reads as the reference base."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGR"
        cols = [_col(5, "C", None, "C", known=False), _col(22, "A", None, "A", known=False)]
        self.assertEqual(_resolve(guide, target, cols, None), [])  # ref PAM TGA
        cols = [_col(5, "C", None, "C", known=False), _col(22, "G", None, "G", known=False)]
        rows = _resolve(guide, target, cols, None)
        self.assertEqual([(r["seq"], r["mm"]) for r in rows],
                         [(GUIDE[:5] + "c" + GUIDE[6:] + "TGG", 1)])

    def test_minus_strand_pam_snp(self):
        """Site B of the synthetic oracle: forward C>T at the PAM's last base reads
        G>A on the minus strand; the cis carrier loses its PAM."""
        guide = GUIDE + "NNN"
        target = GUIDE + "AGR"
        cols = [_col(22, "C", {"T": {"S06": "1|0", "S07": "0|1", "S08": "1|0"}}, "R",
                     strand="-")]
        rows = _resolve(guide, target, cols, {"S05": "1|0", "S06": "1|0", "S07": "1|0"})
        self.assertEqual([(r["seq"], r["samples"]) for r in rows],
                         [(GUIDE + "AGG", ["S05", "S07"])])

    def test_a_bulge_next_to_the_pam(self):
        """A DNA bulge (gap in the crRNA) right before the PAM, its base a code: the
        bulge column is no mismatch; an RNA bulge before an R in the PAM."""
        guide = GUIDE[:-1] + "-" + GUIDE[-1] + "NNN"
        target = GUIDE[:-1] + "Y" + GUIDE[-1] + "TGG"
        cols = [_col(19, "C", {"T": {"S1": "1|0"}}, "Y")]
        rows = _resolve(guide, target, cols, {"S1": "1|0", "S2": "1|0"}, bulges=1)
        got = sorted((r["seq"], r["mm"], tuple(r["samples"])) for r in rows)
        self.assertEqual(got, [(GUIDE[:-1] + "C" + GUIDE[-1] + "TGG", 0, ("S2",)),
                               (GUIDE[:-1] + "T" + GUIDE[-1] + "TGG", 0, ("S1",))])
        guide = GUIDE + "NNN"
        target = GUIDE[:-1] + "-TGR"
        cols = [_col(22, "G", {"A": {"S1": "0|1"}}, "R")]
        rows = _resolve(guide, target, cols, {"S1": "1|1"}, bulges=1)
        self.assertEqual([(r["seq"], r["samples"], r["slots"]) for r in rows],
                         [(GUIDE[:-1] + "-TGG", ["S1"], 1)])

    def test_split_multiallelic_records(self):
        """C>A and C>T at one position (two VCF records): the enricher's code holds
        one alt (Y = C/T), but each slot gets the allele it carries, A included."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "y" + GUIDE[6:] + "TGG"  # guide A; code Y has no A
        cols = [_col(5, "C", {"A": {"S1": "1|0"}, "T": {"S1": "0|1", "S2": "1|0"}}, "Y")]
        rows = _resolve(guide, target, cols, {"S1": "1|1", "S2": "1|0", "S3": "1|0"})
        got = sorted((r["seq"][5], r["mm"], tuple(r["samples"])) for r in rows)
        self.assertEqual(got, [("A", 0, ("S1",)), ("c", 1, ("S3",)), ("t", 1, ("S1", "S2"))])

    def test_two_alts_on_one_slot_is_not_a_haplotype(self):
        """Two ALTs of one position on the indel's slot: PUTATIVE, among those ALTs -- its
        two ALT copies fill a diploid's slots, so the reference is no option (METHODS
        section 4); the indel on the other slot, which holds the reference, is exact."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "y" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"S1": "1|0"}, "T": {"S1": "1|0"}}, "Y")]
        rows = _resolve(guide, target, cols, {"S1": "1|0"})
        self.assertEqual([r["phase"] for r in rows], [isc.PUTATIVE])
        cols = [_col(5, "C", {"A": {"S1": "1|0", "S2": "1|0"},
                              "T": {"S1": "1|0", "S2": "1|0"}}, "Y")]
        exact, putative = isc.resolve_carriers(cols, {"S1": "1|0", "S2": "0|1"})
        self.assertEqual(putative, {(frozenset("AT"),): {"S1"}})
        self.assertEqual(exact, {("C",): [{"S2"}, 1]})

    def test_changes_span(self):
        cs = isc.changes_span
        # deletion TAC>T at 101: deletes 0-based 101, 102
        self.assertTrue(cs(101, "TAC", "T", 95, 101))
        self.assertTrue(cs(101, "TAC", "T", 102, 120))
        self.assertFalse(cs(101, "TAC", "T", 103, 120))   # after the deleted bases
        self.assertFalse(cs(101, "TAC", "T", 80, 100))    # the anchor (0-based 100) only
        # insertion T>TGG at 101: between 0-based 100 and 101
        self.assertTrue(cs(101, "T", "TGG", 90, 110))
        self.assertFalse(cs(101, "T", "TGG", 101, 110))   # starts after the gap
        self.assertFalse(cs(101, "T", "TGG", 90, 100))    # ends at the anchor
        self.assertFalse(cs(101, "A", "A", 0, 1000))      # no change
        # bases the window's own indel has deleted already: 0-based 100..120
        self.assertFalse(cs(106, "CGT", "C", 90, 130, (100, 120)))   # all inside
        self.assertTrue(cs(119, "CGTA", "C", 90, 130, (100, 120)))   # 121 is outside
        self.assertTrue(cs(99, "CGTA", "C", 90, 130, (100, 120)))    # 99 is outside
        self.assertTrue(cs(106, "C", "CGT", 90, 130, (100, 120)))    # bases are added
        self.assertEqual(isc.kept_bases("TAC", "T"), 1)
        self.assertEqual(isc.kept_bases("T", "TGG"), 1)
        self.assertEqual(isc.kept_bases("AT", "CG"), 0)

    def test_target_ref_span_takes_in_the_indel_place(self):
        """The ALT bases, anchor included, map to no reference position, so a target
        that ends inside them would stop one base short of the indel; the span takes
        in the indel's REF. The genome-wide case: chr2:148517281(-) on the window of
        A>ACTT at 148517302, whose only carrier HG00344 also holds T>TC at 148517301
        on the same slot (a C inserted next to the anchor, inside the target)."""
        dna = "AgTaAGAAgMaGTAAG-AGGCTRG"
        real_of = isc.build_offset_to_real(dna, 48764595, "-", 48764590, 148517276,
                                           "A", "ACTT")
        reals = [r for r in real_of if r is not None]
        self.assertEqual((min(reals), max(reals)), (148517281, 148517300))
        self.assertFalse(isc.changes_span(148517301, "T", "TC", 148517281, 148517300))
        span = isc.target_ref_span(dna, real_of, 148517302, "A")
        self.assertEqual(span, (148517281, 148517301))
        self.assertTrue(isc.changes_span(148517301, "T", "TC", *span))
        # a deletion's span takes in the bases its ALT replaces
        self.assertEqual(isc.target_ref_span("ACGT", [10, 11, None, 15], 13, "ATT"),
                         (10, 15))
        # a target with no ALT base (a flank, or a bulge column) is its own bases
        self.assertEqual(isc.target_ref_span("AC-T", [10, 11, None, 12], 40, "A"), (10, 12))
        self.assertIsNone(isc.target_ref_span("--", [None, None], 40, "A"))

    def test_overlapping_indels(self):
        """A second indel inside the target on the SAME slot: that slot does not
        have the window's sequence. On the other slot it does not matter; unphased,
        the sample is PUTATIVE."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"BOTH": "1|0", "OTHERSLOT": "1|0", "UNPH": "1|0"}}, "M")]
        gts = {"BOTH": "1|0", "OTHERSLOT": "1|0", "UNPH": "1|0", "PLAIN": "1|0"}
        second = {"BOTH": "1|0", "OTHERSLOT": "0|1", "UNPH": "0/1"}
        rows = _resolve(guide, target, cols, gts, others=[second])
        got = sorted((r["mm"], tuple(r["samples"]), r["phase"]) for r in rows)
        self.assertEqual(got, [(0, ("OTHERSLOT", "UNPH"), isc.PUTATIVE),
                               (1, ("PLAIN",), isc.CONFIRMED)])
        # and the carriers of an indel that never shares the window are never listed
        rows = _resolve(guide, target, cols, {"PLAIN": "1|0"})
        self.assertEqual([r["samples"] for r in rows], [["PLAIN"]])

    def test_haploid_chrx_calls(self):
        """A haploid call has one haplotype, so cis is certain: exact. A haploid
        indel against a diploid SNP call cannot be put on a slot: PUTATIVE. Carriers
        of one sequence share its row, PUTATIVE when any of them is."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"MALE1": "1", "MALE3": "0|1"}}, "M")]
        gts = {"MALE1": "1", "MALE2": "1", "MALE3": "1"}
        exact, putative = isc.resolve_carriers(cols, gts)
        self.assertEqual(exact, {("A",): [{"MALE1"}, 1], ("C",): [{"MALE2"}, 1]})
        self.assertEqual(putative, {(frozenset("AC"),): {"MALE3"}})
        rows = _resolve(guide, target, cols, gts)
        got = sorted((r["mm"], tuple(r["samples"]), r["phase"]) for r in rows)
        self.assertEqual(got, [(0, ("MALE1", "MALE3"), isc.PUTATIVE),
                               (1, ("MALE2",), isc.CONFIRMED)])

    def test_the_edit_budget(self):
        guide = GUIDE + "NNN"
        # 4 engine mismatches + a code the guide matches only through the alt
        target = "tgac" + GUIDE[4:5] + "M" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"ALT": "1|0"}}, "M")]
        gts = {"ALT": "1|0", "REF": "1|0"}
        rows = _resolve(guide, target, cols, gts)
        self.assertEqual([(r["mm"], r["samples"]) for r in rows], [(4, ["ALT"])])
        rows = _resolve(guide, target, cols, gts, max_mm=5)
        self.assertEqual(sorted((r["mm"], tuple(r["samples"])) for r in rows),
                         [(4, ("ALT",)), (5, ("REF",))])
        rows = _resolve(guide, target, cols, gts, max_mm=5, max_total=4)
        self.assertEqual([(r["mm"], r["samples"]) for r in rows], [(4, ["ALT"])])
        self.assertEqual(_resolve(guide, target, cols, gts, max_total=5, bulges=2), [])

    def test_the_worst_case_is_the_cfd_maximum(self):
        """Brute force: over random codes, the PUTATIVE worst case is the highest
        CFD any combination of the candidate alleles reaches with a valid PAM."""
        rng = random.Random(281)
        guide = GUIDE + "NNN"
        amb_codes = "RYSWKMBDHV"
        for _ in range(300):
            target = list(GUIDE + rng.choice("ACGT") + "GG")
            cols_at = sorted(rng.sample(range(23), rng.randint(1, 4)))
            for j in cols_at:
                target[j] = rng.choice(amb_codes)
            target = "".join(target)
            if not pam_filter._pam_can_match(target[-3:], "NRG"):
                continue
            cols = [_col(j, MOTIF_BASES[target[j]][0], None, target[j], known=False)
                    for j in cols_at]
            rows = _resolve(guide, target, cols, None, max_mm=20)
            best = None
            for combo in itertools.product(*[MOTIF_BASES[target[j]] for j in cols_at]):
                seq = isc.realize(target, guide, cols, combo)
                if not isc.pam_satisfied(seq, isc.pam_columns(23, "NRG")):
                    continue
                s = _cfd(guide, seq)
                best = s if best is None else max(best, s)
            self.assertEqual(len(rows), 1, target)
            self.assertAlmostEqual(_cfd(guide, rows[0]["seq"]), best, places=12,
                                   msg=target)

    def test_pam_window_snps_count_in_the_cooccurrence_decode(self):
        """indel_snp_cis.used_snps_for_target: a SNP that CREATES the PAM is used
        (X3); without the motif map a PAM column is skipped, as before."""
        guide, target = GUIDE + "NNN", GUIDE + "TGR"
        snps = {122: ("A", "G", "rsPAM", {"S": "1|0"})}
        motif_at = isc.pam_columns(len(target), "NRG")
        used = isc.used_snps_for_target(guide, target, lambda j: 100 + j, snps.get,
                                        pam_motif_at=motif_at)
        self.assertEqual(used, [(122, "rsPAM", {"S": "1|0"})])
        self.assertEqual(isc.used_snps_for_target(guide, target, lambda j: 100 + j,
                                                  snps.get), [])
        # an alt that does not create the PAM (ref already allows it) is incidental
        snps = {122: ("G", "A", "rsX", {"S": "1|0"})}
        self.assertEqual(isc.used_snps_for_target(guide, target, lambda j: 100 + j,
                                                  snps.get, pam_motif_at=motif_at), [])

    def test_enabling_columns(self):
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "M" + GUIDE[6:] + "TGR"
        cols = [_col(5, "C", {"A": {"S": "1|0"}}, "M"), _col(22, "A", {"G": {"S": "1|0"}}, "R")]
        rows = _resolve(guide, target, cols, {"S": "1|0"})
        self.assertEqual([(r["seq"], r["enabling"]) for r in rows],
                         [(GUIDE + "TGG", [0, 1])])

    def test_classic_targets_have_no_code_column(self):
        self.assertEqual(isc.ambiguous_columns("ACGTACGTacgt-ACGTN"), [])
        self.assertEqual(isc.ambiguous_columns("ACrTyACG"), [2, 4])


class ASecondIndelIsReadOnTheBases(unittest.TestCase):
    """indel_snp_cis.local_haplotype / target_pairs / holds_target_bases: a carrier whose
    haplotype holds a second indel inside an indel target lacks the window's sequence only
    when that indel changes the target's bases. Decided on the position alone
    (changes_span), as at 9c6b4046, in a tandem repeat, where a unit gained or lost at the
    target's edge leaves the bases as they are, carriers that do hold the row's sequence
    were taken off, and a row left with none vanished. The three real shapes below are
    chr22 of the HGDP + 1000 Genomes panel (hg38 bases, 0-based offsets); dataset C is
    the review's synthetic case (a 1-bp deletion in the PAM's G run, and a spacer deletion
    that does change the target); dataset A its record left unnormalized inside the
    target's 5' repeat. Each "holds" fails on the position rule."""

    # chr22:48348960-48349080 and 48349080-48349180 (0-based, end exclusive), hg38
    C22A = ("CTCACACACACTCTCACACAGACTCACACTCACACACACTCTCACACAGACACGCTCACACACATTCACATACTTACA"
            "CTCTCACACAGACTCACACTCACACATACACTGTCACACAGT")
    C22B = ("CACATGCACTCATATACATACTTACACACACTCTCTCACACACACTCAGCCAGTCACACACACACCCACAGAGACTTT"
            "CACACACAGACTCACACACACA")
    # chr22:48634700-48634800
    C22C = ("TCATTCATTCACCCTCTCATTCATCTACTCATTCACTCTCTCACTCATTCATTCACTCAGTCACTCAGCCACTCATTTAG"
            "TCACTCACTCATTCACTGTC")

    @staticmethod
    def _bases(seq, at):
        def value(src):
            return src if isinstance(src, str) else seq[src - at]
        return value

    def _holds(self, seq, at, window, others, first, n, value=None):
        """Whether a slot carrying ``window`` and ``others`` (VCF POS, REF, ALT) holds the
        window's n bases from the place ``first`` on; ``seq`` holds the reference from
        0-based ``at``."""
        lo, hi = at, at + len(seq) - 1
        on_window = isc.local_haplotype(lo, hi, window)
        keys = [k for k, _s in on_window]
        last = keys[keys.index(first) + n - 1]
        hap = isc.local_haplotype(lo, hi, window, others)
        ways = isc.target_pairs(on_window, hap, first, last)
        return isc.holds_target_bases(ways, value or self._bases(seq, at))

    def _target(self, seq, at, window, first, n):
        on_window = isc.local_haplotype(at, at + len(seq) - 1, window)
        keys = [k for k, _s in on_window]
        i = keys.index(first)
        value = self._bases(seq, at)
        return "".join(value(s) for _k, s in on_window[i:i + n])

    def test_vegfa3_a_deletion_takes_the_last_base_and_the_repeat_puts_it_back(self):
        """chr22:48349010(-), the VEGFA3 guide: A>ACT makes the target; on 1,180 of its
        1,199 carriers the same slot holds CAT>C at 48349029, whose AT deletion takes the
        target's last forward base, and the repeat puts an A back. HG00100 slot 0."""
        window, second = (48349010, "A", "ACT"), (48349029, "CAT", "C")
        first = isc.window_place(25 + 1, 48349010, 48349010 - 26, "A", "ACT")
        self.assertEqual(first, (48349008, 2, 1))  # the first inserted base
        self.assertEqual(self._target(self.C22A, 48348960, window, first, 22),
                         "CTCACGCTCACACACATTCACA")
        self.assertTrue(isc.changes_span(*second, 48349009, 48349029))  # the position rule
        self.assertTrue(self._holds(self.C22A, 48348960, window, [second], first, 22))
        # each way is read on its own (_breaking_slots), so there are always two: a target
        # the window written out does not hold has neither
        on_window = isc.local_haplotype(48348960, 48348960 + len(self.C22A) - 1, window)
        self.assertEqual(isc.target_pairs(on_window, on_window, (48340000, 0), first),
                         [None, None])

    def test_a_repeat_unit_inserted_before_the_last_base(self):
        """chr22:48349112(-) (AGT>A at 48349132; NA19327's slot also holds A>ACAGT at
        48349136 and ACACT>A at 48349108, just before the target) and chr22:48634743(-)
        (CTCAT>C at 48634744; NA18907 holds C>CCACT at 48634769 on both slots): the
        inserted unit's first base is the one it displaces. Each was the only carrier: the
        row vanished."""
        window = (48349132, "AGT", "A")
        first = (48349112, 0)
        self.assertEqual(self._target(self.C22B, 48349080, window, first, 23),
                         "CTCTCACACACACTCAGCCACAC")
        others = [(48349108, "ACACT", "A"), (48349136, "A", "ACAGT")]
        self.assertTrue(self._holds(self.C22B, 48349080, window, others, first, 23))
        window = (48634744, "CTCAT", "C")
        first = isc.window_place(25, 48634744, 48634744 - 26, "CTCAT", "C")
        self.assertEqual(first, (48634742, 2, 0))  # the anchor: the window writes the ALT
        self.assertEqual(self._target(self.C22C, 48634700, window, first, 23),
                         "CTCATTCACTCAGTCACTCAGCC")
        self.assertTrue(self._holds(self.C22C, 48634700, window, [(48634769, "C", "CCACT")],
                                    first, 23))
        # an insertion there whose first base is not the one it displaces changes it
        self.assertFalse(self._holds(self.C22C, 48634700, window,
                                     [(48634769, "C", "CTACTCACT")], first, 23))

    # dataset C of the review (make_cases.py): a T, the guide CACACACTTGCT with an extra G
    # the window's TG>T deletes, TTTATCAC, then the PAM TGG with one more G and an A
    SITE_C = "T" + "CACACACTTGCT" + "G" + "TTTATCAC" + "TGGG" + "A" + "CCCC"

    def test_a_deletion_in_the_pams_g_run_and_a_spacer_deletion(self):
        seq, at = self.SITE_C, 1999
        window = (2012, "TG", "T")  # 0-based 2011 is the T before the extra G
        first = (2000, 0)
        self.assertEqual(self._target(seq, at, window, first, 23),
                         "CACACACTTGCTTTTATCACTGG")
        pam_run = (2022, "TG", "T")   # S1: the PAM keeps its G from the run
        spacer = (2016, "TA", "T")    # S2: the spacer loses a base
        self.assertTrue(isc.changes_span(*pam_run, 2000, 2023))
        self.assertTrue(self._holds(seq, at, window, [pam_run], first, 23))
        self.assertFalse(self._holds(seq, at, window, [spacer], first, 23))
        self.assertFalse(self._holds(seq, at, window, [pam_run, spacer], first, 23))

    def test_a_unit_inserted_in_the_5prime_repeat_by_a_record_left_unnormalized(self):
        """Dataset A: the CA insertion written at a POS inside the target's 5' CA repeat
        (bcftools norm moves it upstream, out of the target). The slot holds the 23-mer
        at the PAM, up to the target's last base, not from its first."""
        guide = "CACACACTTGCTTTTATCAC"
        seq = "GTA" + "CACACACACACA" + guide[:12] + "G" + guide[12:] + "TGG" + "ATC"
        at = 1985
        g0 = at + 3 + 12
        window = (g0 + 12, "TG", "T")
        first = (g0, 0)
        self.assertEqual(self._target(seq, at, window, first, 23), guide + "TGG")
        self.assertTrue(self._holds(seq, at, window, [(g0 + 2, "A", "ACA")], first, 23))
        # the same unit, left-aligned (bcftools norm), is outside the target: no change
        self.assertFalse(isc.changes_span(at + 3, "A", "ACA", g0, g0 + 23))

    def test_a_snp_the_shift_brings_into_the_target(self):
        """The bases a deletion brings in from past the target are the slot's own: a SNP
        alt there, on the same slot, changes the target; on the other slot it does not."""
        seq, at = self.SITE_C, 1999
        window, pam_run = (2012, "TG", "T"), (2022, "TG", "T")
        first = (2000, 0)
        brought = 2024  # the run's last G: the PAM's last base on the slot

        def on_slot(alt_here):
            def value(src):
                if isinstance(src, str):
                    return src
                return alt_here if src == brought and alt_here else seq[src - at]
            return value

        self.assertTrue(self._holds(seq, at, window, [pam_run], first, 23, on_slot(None)))
        self.assertFalse(self._holds(seq, at, window, [pam_run], first, 23, on_slot("A")))
        # an allele the genotypes cannot put on the slot is not the same base
        self.assertFalse(self._holds(seq, at, window, [pam_run], first, 23,
                                     lambda s: s if isinstance(s, str) else (
                                         None if s == brought else seq[s - at])))

    # chr22:33466660-33466760
    C22D = ("AGCCCTGGCCTCATTCACCCCCTCTCTCTCTCTACACACACACACACACACACATACACACACACACTACACACACACACC"
            "TTAAGAGTTTGTTATCACC")

    def test_the_bases_must_come_through_the_indel(self):
        """chr22:33466691(-) on the window of CACACAT>C at 33466709: HGDP01094's slot also
        holds T>TACACAC at 33466693, inside the target. From the target's first base its
        haplotype has the target's bases, but the three CA units gained upstream push the
        deletion's junction out of that stretch: it holds the bases through its insertion,
        not through the deletion the row names, and is not listed for it (the audit's
        carrier check calls such a listing a failure)."""
        window, ins = (33466709, "CACACAT", "C"), (33466693, "T", "TACACAC")
        first = (33466691, 0)
        self.assertEqual(self._target(self.C22D, 33466660, window, first, 23),
                         "CTACACACACACACACACACACA")
        lo, hi = 33466660, 33466660 + len(self.C22D) - 1
        hap = isc.local_haplotype(lo, hi, window, [ins])
        hkeys = [k for k, _s in hap]
        stretch = hap[hkeys.index(first):hkeys.index(first) + 23]
        value = self._bases(self.C22D, 33466660)
        self.assertEqual("".join(value(src) for _k, src in stretch), "CTACACACACACACACACACACA")
        self.assertFalse(self._holds(self.C22D, 33466660, window, [ins], first, 23))

    def test_variants_that_overlap_are_no_haplotype(self):
        self.assertIsNone(isc.local_haplotype(100, 140, (120, "TCA", "T"), [(121, "CA", "C")]))
        self.assertIsNone(isc.local_haplotype(100, 140, (120, "T", "TC"),
                                              [(110, "A", "AG"), (110, "A", "AT")]))
        # Two second indels whose REFs share a base, in either order, though the bases each
        # changes do not meet: at one anchor a tandem duplication and the deletion of the
        # same bases (the review's X1, C>CTA and CTA>C), which f1b75ba added together into
        # the reference when the insertion came first, as _breaking_slots sorts them; an
        # insertion anchored on a base another record deletes; a record left unnormalized
        # (TCA>TA) whose shared last base anchors an insertion.
        window = (130, "T", "TC")
        for pair in ([(110, "A", "ACA"), (110, "ACA", "A")],
                     [(110, "ACA", "A"), (112, "A", "AG")],
                     [(110, "TCA", "TA"), (112, "A", "AG")]):
            for others in (pair, pair[::-1]):
                self.assertIsNone(isc.local_haplotype(105, 116, window, others), others)
        # an indel inserted right before the window's anchor comes first
        hap = isc.local_haplotype(118, 122, (121, "G", "GTT"), [(120, "C", "CA")])
        self.assertEqual([s for _k, s in hap], [118, 119, "A", "G", "T", "T", 121, 122])
        # an insertion where the bases a record that keeps none deletes start (TA>A at the
        # next base) comes first, whatever the order: f1b75ba gave None with the deletion first
        pair = [(110, "A", "ACA"), (111, "TA", "A")]
        for others in (pair, pair[::-1]):
            hap = isc.local_haplotype(108, 131, window, others)
            self.assertEqual([s for _k, s in hap],
                             [108, 109, "C", "A"] + list(range(111, 129)) + ["T", "C", 130, 131])

    def test_one_answer_whatever_the_order_of_the_variants(self):
        """Every set of records, in every order, gives one answer: None when two REFs share
        a base, else the bases of the haplotype that holds them all (each REF replaced by
        its ALT, by position), or None for two insertions at one point, whose bases have no
        one order. At f1b75ba the answer of an insertion and a deletion at one anchor, or
        of an insertion where the bases a record keeping none deletes start, depended on
        which came first."""
        rng = random.Random(2810)
        ref = "".join(rng.choice("ACGT") for _ in range(60))

        def record(p1):
            r0 = ref[p1 - 1]
            new = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3)))
            kind = rng.randrange(5)
            if kind == 0:
                return (p1, r0, r0 + new)                       # an insertion
            if kind == 1:
                return (p1, ref[p1 - 1:p1 + rng.randint(1, 3)], r0)  # a deletion
            if kind == 2 and new[0] != r0:
                return (p1, r0, new + r0)                       # inserted before REF
            if kind == 3 and ref[p1] != r0:
                return (p1, ref[p1 - 1:p1 + 1], ref[p1])        # deletes its first base
            return (p1, ref[p1 - 1:p1 + 2], r0 + ref[p1 + 1])   # left unnormalized

        def holding_all(records):
            out, at = [], 0
            for p1, r, a in sorted(records):
                out += [ref[at:p1 - 1], a]
                at = p1 - 1 + len(r)
            return "".join(out) + ref[at:]

        seen = {"overlap": 0, "bases": 0, "one point": 0}
        for _ in range(3000):
            window, *others = [record(rng.randint(20, 32)) for _ in range(rng.randint(2, 4))]
            answers = set()
            for order in itertools.permutations(others):
                hap = isc.local_haplotype(0, len(ref) - 1, window, list(order))
                answers.add(None if hap is None else tuple(hap))
            self.assertEqual(len(answers), 1, (window, others))
            hap = answers.pop()
            refs = sorted((p1 - 1, p1 - 1 + len(r)) for p1, r, _a in [window] + others)
            if any(b[0] < a[1] for a, b in zip(refs, refs[1:])):
                self.assertIsNone(hap, (window, others))
                seen["overlap"] += 1
            elif hap is not None:
                self.assertEqual("".join(s if isinstance(s, str) else ref[s] for _k, s in hap),
                                 holding_all([window] + others), (window, others))
                seen["bases"] += 1
            else:
                self.assertGreaterEqual(sum(len(a) > len(r) for _p, r, a in others), 2)
                seen["one point"] += 1
        self.assertTrue(all(seen.values()), seen)


class TargetsThatHoldTheirIndel(unittest.TestCase):
    """indel_snp_cis.holds_indel: a fake-indel window is 25 bp of plain flank, the ALT,
    then the downstream flank. A target wholly in a flank is reference (or SNP) sequence
    that needs no indel; its row showed the indel's token at a coordinate the indel
    shifts, and could take the real row's place in its cluster."""

    POS1 = 1026   # VCF POS: the window starts at 1000, so the REF's first base is offset 25
    START = 1000

    def _holds(self, first, last, ref, alt):
        return isc.holds_indel(first, last, self.POS1, self.START, ref, alt)

    def test_a_deletion_is_held_across_its_junction(self):
        ref, alt = "TC", "T"
        for first, last, want, why in (
                (3, 25, False, "ends at the anchor: reference sequence"),
                (3, 26, True, "the last base is past the junction"),
                (25, 47, True, "starts at the anchor"),
                (26, 48, False, "starts past the junction: the downstream flank"),
                (40, 62, False, "well inside the downstream flank"),
                (0, 22, False, "inside the left flank")):
            self.assertEqual(self._holds(first, last, ref, alt), want, why)

    def test_a_long_deletion_is_the_same_on_the_window(self):
        for ref in ("TCAG", "TCAGGTCAGGTCAGGTCAGGTCAGGTCAGG"):
            self.assertFalse(self._holds(3, 25, ref, "T"))
            self.assertTrue(self._holds(3, 26, ref, "T"))
            self.assertFalse(self._holds(26, 48, ref, "T"))

    def test_an_insertion_is_held_through_an_inserted_base(self):
        ref, alt = "G", "GCAT"  # ALT at 25..28, inserted bases 26..28
        for first, last, want, why in (
                (3, 25, False, "holds the anchor alone: its sequence is the reference's"),
                (4, 26, True, "ends at the first inserted base"),
                (26, 48, True, "starts at the first inserted base"),
                (28, 50, True, "starts at the last inserted base"),
                (29, 51, False, "starts past the insertion: the downstream flank")):
            self.assertEqual(self._holds(first, last, ref, alt), want, why)
        # a target wholly inside a long insertion spans no reference base, and needs it
        self.assertTrue(self._holds(30, 52, "G", "G" + "CAT" * 20))

    def test_a_replacement_is_held_through_any_replacing_base(self):
        ref, alt = "CT", "GAA"  # nothing kept: ALT at 25..27
        for first, last, want in ((3, 24, False), (3, 25, True), (27, 49, True),
                                  (28, 50, False)):
            self.assertEqual(self._holds(first, last, ref, alt), want, (first, last))

    def test_more_than_one_kept_base(self):
        """An ALT sharing two bases with REF: the junction is after both."""
        ref, alt = "ACG", "AC"  # ALT at 25..26, the junction between 26 and 27
        for first, last, want in ((3, 26, False), (3, 27, True), (26, 48, True),
                                  (27, 49, False)):
            self.assertEqual(self._holds(first, last, ref, alt), want, (first, last))

    def test_a_record_left_unnormalized(self):
        """TCA>TA keeps the A: the record deletes one C (bases either side: 25, 26), so a
        target starting at that A is reference sequence. TA>TCA adds one C at offset 26."""
        for first, last, want in ((3, 25, False), (3, 26, True), (26, 48, False)):
            self.assertEqual(self._holds(first, last, "TCA", "TA"), want, (first, last))
        for first, last, want in ((3, 25, False), (4, 26, True), (26, 48, True),
                                  (27, 49, False)):
            self.assertEqual(self._holds(first, last, "TA", "TCA"), want, (first, last))
        self.assertEqual(isc.shared_tail("TCA", "TA"), 1)
        self.assertEqual(isc.shared_tail("TA", "TCA"), 1)
        self.assertEqual(isc.shared_tail("TC", "T"), 0)
        self.assertEqual(isc.shared_tail("AA", "A"), 0)  # no base is counted twice
        # the A a second variant keeps changes nothing either: only 0-based 101 is lost
        self.assertTrue(isc.changes_span(101, "TCA", "TA", 95, 101))
        self.assertFalse(isc.changes_span(101, "TCA", "TA", 102, 120))
        self.assertFalse(isc.changes_span(101, "TA", "TCA", 101, 120))  # past the insertion
        self.assertTrue(isc.changes_span(101, "TA", "TCA", 95, 110))

    def test_a_shorter_left_flank(self):
        """The window of an indel near the chromosome's start: the REF's first base is
        found from the window's own start, not at a fixed offset."""
        self.assertTrue(isc.holds_indel(0, 22, 10, 0, "TC", "T"))
        self.assertFalse(isc.holds_indel(10, 32, 10, 0, "TC", "T"))

    def test_window_offsets_count_the_bases_only(self):
        self.assertEqual(isc.window_offsets("ACGT", 110, 100), (10, 13))
        self.assertEqual(isc.window_offsets("A-GT", 110, 100), (10, 12))  # an RNA bulge
        self.assertEqual(isc.window_offsets("A-Gt--A", 0, 0), (0, 3))


class TheBadTargetFilter(unittest.TestCase):
    """remove_bad_indel_targets.py reads a row's indel from its variant list."""

    def _filter(self, rows):
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "best.txt")
            with open(path, "w") as fh:
                fh.write("header\n")
                for start, var in rows:
                    cols = ["n"] * 22
                    cols[1] = GUIDE + "NNN"
                    cols[5] = str(start)
                    cols[18] = var
                    fh.write("\t".join(cols) + "\n")
            subprocess.run([sys.executable, os.path.join(HERE, "remove_bad_indel_targets.py"),
                            path], check=True, cwd=d)
            with open(path) as fh:
                return [line.split("\t")[18] for line in fh.readlines()[1:]]

    def test_a_row_is_judged_on_its_indel_alone(self):
        """'chr1_100_GA_G,chr1_120_A_G'.split('_') used to read the ALT as 'G,chr1', a
        range six bases wide."""
        one = "chr1_100_GA_G"
        two = one + ",chr1_120_A_G"
        self.assertEqual(self._filter([(90, one), (90, two)]), [one, two])
        self.assertEqual(self._filter([(102, one), (102, two)]), [],
                         "past the indel's own range, whatever else the row names")

    def test_a_contig_name_with_underscores(self):
        tok = "chr1_KI270706v1_random_100_GA_G"
        self.assertEqual(self._filter([(90, tok)]), [tok])
        self.assertEqual(self._filter([(102, tok)]), [])

    def test_reference_rows_pass(self):
        self.assertEqual(self._filter([(5, "n")]), ["n"])


class MultiAllelicRecordsOnOneLine(unittest.TestCase):
    """indel_snp_cis.biallelic_site: a record with several alts on one VCF line keeps
    its raw genotype under every alt its sample carries (Tier-1 and the per-sample
    dictionary alike), and the resolution reads "1" as "this slot holds this alt"."""

    def test_the_table(self):
        site = {"A": {"S26": "0|1", "S27": "2|1"},
                "T": {"S25": "2|2", "S27": "2|1", "S28": "0|2", "S29": "1|2", "S30": "2"}}
        site["A"]["S29"] = "1|2"
        got = isc.biallelic_site(site)
        self.assertEqual(got["A"], {"S26": "0|1", "S27": "0|1", "S29": "1|0"})
        self.assertEqual(got["T"], {"S25": "1|1", "S27": "1|0", "S28": "0|1", "S29": "0|1",
                                    "S30": "1"})

    def test_unphased_and_missing_alleles_stay_so(self):
        got = isc.biallelic_site({"A": {"U": "1/2", "M": ".|1"}, "T": {"U": "1/2", "V": "0/2",
                                                                       "W": ".|2"}})
        self.assertEqual(got["A"], {"U": "1/0", "M": ".|1"})
        self.assertEqual(got["T"], {"U": "0/1", "V": "0/1", "W": ".|1"})

    def test_a_single_carrier_needs_the_other_alt(self):
        """2|1, the only carrier of A: its genotype alone holds both tokens; T's other
        carrier is 2|2, so A is 1. The most frequent token (the SNP path's rule) ties."""
        got = isc.biallelic_site({"A": {"S": "2|1"}, "T": {"S": "2|1", "X": "2|2"}})
        self.assertEqual(got["A"], {"S": "0|1"})
        self.assertEqual(got["T"], {"S": "1|0", "X": "1|1"})

    def test_a_lone_foreign_index(self):
        got = isc.biallelic_site({"T": {"S": "2|0", "X": "2|2"}})
        self.assertEqual(got["T"], {"S": "1|0", "X": "1|1"})

    def test_genotypes_that_cannot_tell_lose_their_phase(self):
        """Every carrier of both alts is 1|2: which slot holds which alt is not known."""
        got = isc.biallelic_site({"A": {"S": "1|2"}, "T": {"S": "1|2"}})
        self.assertEqual(got, {"A": {"S": "1/0"}, "T": {"S": "1/0"}})
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "H" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", got, "H")]
        rows = _resolve(guide, target, cols, {"S": "1|0"})
        self.assertEqual([(r["seq"][5], r["samples"], r["phase"]) for r in rows],
                         [("A", ["S"], isc.PUTATIVE)])

    def test_split_and_biallelic_records_come_back_as_they_are(self):
        for site in ({"A": {"S": "1|0"}}, {"A": {"S": "1|0"}, "T": {"S": "0|1", "X": "1|1"}},
                     {"A": {"S": "1/0", "H": "1", "M": ".|1"}}, {}):
            self.assertIs(isc.biallelic_site(site), site)

    def test_each_slot_gets_its_own_allele(self):
        """The reviewer's site: C>A,T on one line, A the guide's base, a deletion on one
        slot. S26's deletion slot holds A, S27's holds T (2|1), S25 is 2|2 (no C at all),
        S28 has the reference: all exact, none PUTATIVE."""
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "H" + GUIDE[6:] + "TGG"
        site = isc.biallelic_site({"A": {"S26": "0|1", "S27": "2|1"},
                                   "T": {"S25": "2|2", "S27": "2|1"}})
        cols = [_col(5, "C", site, "H")]
        rows = _resolve(guide, target, cols,
                        {"S25": "1|0", "S26": "0|1", "S27": "1|0", "S28": "1|0"})
        got = sorted((r["seq"][5], tuple(r["samples"]), r["phase"]) for r in rows)
        self.assertEqual(got, [("A", ("S26",), isc.CONFIRMED),
                               ("c", ("S28",), isc.CONFIRMED),
                               ("t", ("S25", "S27"), isc.CONFIRMED)])

    def test_unphased_homozygous_foreign_index_has_no_reference_allele(self):
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "H" + GUIDE[6:] + "TGG"
        site = isc.biallelic_site({"A": {"X": "1/0"}, "T": {"S": "2/2", "H": "0/2"}})
        cols = [_col(5, "C", site, "H")]
        rows = _resolve(guide, target, cols, {"S": "1|0"})
        self.assertEqual([(r["seq"][5], r["mm"], r["phase"]) for r in rows],
                         [("t", 1, isc.PUTATIVE)])
        # a heterozygote does have one: its worst case over C and T
        rows = _resolve(guide, target, cols, {"H": "1|0"})
        self.assertEqual([r["phase"] for r in rows], [isc.PUTATIVE])
        self.assertIn(rows[0]["seq"][5], "ct")


# --------------------------------------------------------------------------- #
# (c) the analyzer routes every overlaid target through the resolver
# --------------------------------------------------------------------------- #
class TheAnalyzerUsesIt(unittest.TestCase):
    def _src(self, name):
        with open(os.path.join(HERE, name)) as fh:
            return fh.read()

    def test_overlaid_targets_go_through_resolve_target(self):
        src = self._src("analisi_indels_NNN.py")
        self.assertIn("_isc.ambiguous_columns(line[2])", src)
        self.assertIn("_isc.resolve_target(", src)
        self.assertIn("_isc.pam_columns(len(line[2]), pam, pam_at_beginning)", src)
        self.assertIn("allowed_mms, max_total_edits", src)
        self.assertNotIn("_resolve_overlaid_iupac", src)

    def test_the_total_edits_cap_reaches_the_analyzer(self):
        self.assertIn('max_total_edits = int(sys.argv[11]) if len(sys.argv) > 11 else -1',
                      self._src("analisi_indels_NNN.py"))
        self.assertIn('max_total_edits=${13:--1}', self._src("analisi_indels_NNN.sh"))
        self.assertIn('max_total_edits=${14:--1}', self._src("post_analisi_indel.sh"))
        self.assertIn("{chrom} {max_total_edits}", self._src("pool_post_analisi_indel.py"))
        sub = self._src("submit_job_automated_new_multiple_vcfs.sh")
        self.assertIn('$final_res $final_res_alt $ncpus "$max_total_edits"', sub)


if __name__ == "__main__":
    unittest.main(verbosity=2)
