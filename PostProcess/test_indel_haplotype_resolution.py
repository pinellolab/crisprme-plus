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
      multi-allelic records, overlapping indels, haploid chrX calls, the edit
      budget, the worst case against a brute force, and the PAM-window SNPs of the
      co-occurrence decode;
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
        guide = GUIDE + "NNN"
        target = GUIDE[:5] + "y" + GUIDE[6:] + "TGG"
        cols = [_col(5, "C", {"A": {"S1": "1|0"}, "T": {"S1": "1|0"}}, "Y")]
        rows = _resolve(guide, target, cols, {"S1": "1|0"})
        self.assertEqual([r["phase"] for r in rows], [isc.PUTATIVE])

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
