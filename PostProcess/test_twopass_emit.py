#!/usr/bin/env python3
"""Tests for twopass_emit (fast-mode worst-possible representative).

The load-bearing property: the greedy per-column representative achieves the TRUE
minimum mismatch over all 2^k allele combinations (so fast mode never under-reports a
window's worst-possible off-target). Brute-forced against every combination on random
windows, including bulges (both strands) and multiallelic columns. STDLIB only.
"""
import os
import random
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import twopass_emit as te  # noqa: E402

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "-": "-"}


def _revcom(s):
    return "".join(_COMP.get(b, "N") for b in reversed(s))


class TestGreedyWorstCase(unittest.TestCase):
    def test_matches_bruteforce_random(self):
        random.seed(11)
        for _ in range(2000):
            glen = random.randint(6, 12)
            guide = "".join(random.choice("ACGT") for _ in range(glen))
            ref = "".join(random.choice("ACGT") for _ in range(glen))
            revert = random.random() < 0.5
            # a plain (ungapped) window: real_target has no '-'; pos_beg=0, pos_end=None
            real_target = ref
            # pick a few variant columns, each with 1-3 alt bases
            ncols = random.randint(1, min(5, glen))
            cols = []
            for pc in sorted(random.sample(range(glen), ncols)):
                alts = random.sample([b for b in "ACGT" if b != ref[pc]],
                                     random.randint(1, 3))
                cols.append({
                    "pos_c": pc,
                    "candidates": [{"alt": a, "carriers": {a}, "info": [a, "0.1", "s"]}
                                   for a in alts],
                })
            got = te.greedy_worst_case(cols, ref, real_target, guide, revert,
                                       0, None, _revcom)
            got_mm = te.aligned_mismatches(got["seq"], real_target, guide, revert,
                                           0, None, _revcom)
            brute = te.brute_force_min_mismatch(cols, ref, real_target, guide, revert,
                                                0, None, _revcom)
            self.assertEqual(got_mm, brute,
                             f"guide={guide} ref={ref} cols={cols} revert={revert}")

    def test_matches_bruteforce_with_bulges(self):
        random.seed(23)
        for _ in range(1000):
            glen = random.randint(6, 10)
            guide = "".join(random.choice("ACGT") for _ in range(glen))
            ref = "".join(random.choice("ACGT") for _ in range(glen))
            # inject 1-2 DNA-bulge gaps into the aligned target (extra window bases)
            real_target = list(ref)
            for _ in range(random.randint(0, 2)):
                real_target.insert(random.randint(0, len(real_target)), "-")
            real_target = "".join(real_target)
            revert = random.random() < 0.5
            ncols = random.randint(1, min(4, glen))
            cols = []
            for pc in sorted(random.sample(range(glen), ncols)):
                alts = random.sample([b for b in "ACGT" if b != ref[pc]], random.randint(1, 2))
                cols.append({
                    "pos_c": pc,
                    "candidates": [{"alt": a, "carriers": set(), "info": [a, "0.1", "s"]}
                                   for a in alts],
                })
            got = te.greedy_worst_case(cols, ref, real_target, guide, revert, 0, None, _revcom)
            got_mm = te.aligned_mismatches(got["seq"], real_target, guide, revert, 0, None, _revcom)
            brute = te.brute_force_min_mismatch(cols, ref, real_target, guide, revert, 0, None, _revcom)
            self.assertEqual(got_mm, brute)

    def test_carriers_and_info_unioned_for_chosen_alts(self):
        # guide == ref except two columns where an alt restores the match; both alts
        # must be chosen and their carriers/info unioned.
        guide = "ACGTACGT"
        ref = "AXGTAYGT".replace("X", "C").replace("Y", "C")  # ref = ACGTACGT? no:
        # make ref differ at col1 (C vs guide C -> same) -- construct explicitly:
        ref = "AAGTAAGT"           # differs from guide at pos1 (A vs C) and pos5 (A vs C)
        cols = [
            {"pos_c": 1, "candidates": [{"alt": "C", "carriers": {"s1"}, "info": ["rs1", "0.2", "i1"]}]},
            {"pos_c": 5, "candidates": [{"alt": "C", "carriers": {"s2"}, "info": ["rs2", "0.3", "i2"]}]},
        ]
        got = te.greedy_worst_case(cols, ref, ref, guide, False, 0, None, _revcom)
        self.assertEqual("".join(got["seq"]), "ACGTACGT")          # both alts applied
        self.assertEqual(got["carriers"], {"s1", "s2"})
        self.assertEqual(got["info"], [["rs1", "0.2", "i1"], ["rs2", "0.3", "i2"]])
        self.assertEqual(te.aligned_mismatches(got["seq"], ref, guide, False, 0, None, _revcom), 0)

    def test_pam_aware_picks_valid_pam_allele_over_lex(self):
        # A PAM-region variant column (mismatch-neutral: outside the protospacer window
        # ref[0:-1]) with alts A and G, where only G yields a valid PAM. Lexicographic
        # order prefers A (invalid). Without pam_valid_fn the greedy picks A (the rep is
        # then PAM-invalid and the off-target would be dropped downstream -- the real
        # chr22 141-miss bug). With pam_valid_fn the greedy must pick G.
        guide = "ACGT"
        ref = "ACGTC"          # protospacer ref[0:-1]="ACGT" matches guide (mm 0); last = PAM
        cols = [{"pos_c": 4, "candidates": [
            {"alt": "A", "carriers": {"sA"}, "info": ["rsA", "0.1", "iA"]},
            {"alt": "G", "carriers": {"sG"}, "info": ["rsG", "0.2", "iG"]}]}]
        # PAM valid iff the last base is G
        pam_ok = lambda seq: seq[4] == "G"

        # PAM-blind: lex order -> A (invalid PAM)
        blind = te.greedy_worst_case(cols, ref, ref, guide, False, 0, -1, _revcom)
        self.assertEqual(blind["seq"][4], "A")
        self.assertFalse(pam_ok(blind["seq"]))

        # PAM-aware: mismatch-neutral -> prefer the valid-PAM allele G
        aware = te.greedy_worst_case(cols, ref, ref, guide, False, 0, -1, _revcom,
                                     pam_valid_fn=pam_ok)
        self.assertEqual(aware["seq"][4], "G")
        self.assertTrue(pam_ok(aware["seq"]))
        self.assertEqual(aware["carriers"], {"sG"})
        self.assertEqual(aware["info"], [["rsG", "0.2", "iG"]])

    def test_joint_pam_two_columns_resolves_valid(self):
        # A valid PAM that needs a specific allele at TWO variant columns -- here the
        # REFERENCE base at col2 AND the alt at col3 (GMR-style). A single left-to-right
        # pass picks the alt at col2 (lex beats the ref on an invalid tie) and gets stuck;
        # the multi-pass / brute-force fallback must still yield a valid-PAM rep.
        guide = "AC"
        ref = "ACAT"                 # protospacer ref[0:-2]="AC" matches guide; PAM cols 2,3
        cols = [
            {"pos_c": 2, "candidates": [{"alt": "G", "carriers": {"s2"}, "info": ["r2", "0.1", "i"]}]},
            {"pos_c": 3, "candidates": [{"alt": "G", "carriers": {"s3"}, "info": ["r3", "0.2", "i"]}]},
        ]
        pam_ok = lambda seq: seq[2] == "A" and seq[3] == "G"   # ref@col2, alt@col3
        r = te.greedy_worst_case(cols, ref, ref, guide, False, 0, -2, _revcom, pam_valid_fn=pam_ok)
        self.assertTrue(pam_ok(r["seq"]), "".join(r["seq"]))
        self.assertEqual("".join(r["seq"]), "ACAG")
        self.assertEqual(r["carriers"], {"s3"})           # col2 kept reference
        self.assertEqual(r["info"], [["r3", "0.2", "i"]])

    def test_pam_repair_helper_finds_min_mm_valid(self):
        guide, ref = "AC", "ACAT"
        cols = [
            {"pos_c": 2, "candidates": [{"alt": "G", "carriers": set(), "info": ["r", "0.1", "i"]}]},
            {"pos_c": 3, "candidates": [{"alt": "G", "carriers": set(), "info": ["r", "0.2", "i"]}]},
        ]
        pam_ok = lambda seq: seq[2] == "A" and seq[3] == "G"
        got = te._pam_valid_min_mm(cols, ref, ref, guide, False, 0, -2, _revcom, pam_ok, 4096)
        self.assertIsNotNone(got)
        self.assertEqual("".join(got["seq"]), "ACAG")
        self.assertIsNone(te._pam_valid_min_mm(cols, ref, ref, guide, False, 0, -2,
                                               _revcom, lambda s: False, 4096))
        self.assertIsNone(te._pam_valid_min_mm(cols, ref, ref, guide, False, 0, -2,
                                               _revcom, pam_ok, 1))   # over cap -> None

    def test_min_mm_chosen_when_it_is_pam_valid(self):
        # When the min-mismatch allele is ALSO PAM-valid it is chosen (mismatch minimized,
        # no PAM override) -- the common case.
        guide = "ACGT"
        ref = "AAGT"           # pos1 ref A != guide C (mm 1); alt C fixes it (mm 0)
        cols = [{"pos_c": 1, "candidates": [
            {"alt": "C", "carriers": {"sC"}, "info": ["rC", "0.2", "s"]}]}]
        r = te.greedy_worst_case(cols, ref, ref, guide, False, 0, None, _revcom,
                                 pam_valid_fn=lambda seq: True)
        self.assertEqual("".join(r["seq"]), "ACGT")   # C chosen (mm 0)
        self.assertEqual(te.aligned_mismatches(r["seq"], ref, guide, False, 0, None, _revcom), 0)

    def test_pam_valid_off_target_preferred_over_invalid_lower_mm(self):
        # The worst-possible off-target must have a VALID PAM. When the LOWER-mismatch
        # allele is PAM-invalid (not a real off-target -> would be dropped by pam_ok,
        # the 141-miss bug) and a higher-mismatch allele is PAM-valid, the greedy must
        # return the VALID one. Mismatch is minimized AMONG the valid-PAM combos.
        guide = "ACGT"
        ref = "AAGT"
        cols = [{"pos_c": 1, "candidates": [
            {"alt": "A", "carriers": {"sA"}, "info": ["rA", "0.1", "s"]},   # mm1, PAM-valid
            {"alt": "C", "carriers": {"sC"}, "info": ["rC", "0.2", "s"]}]}]  # mm0, PAM-invalid
        pam_ok = lambda seq: seq[1] == "A"
        r = te.greedy_worst_case(cols, ref, ref, guide, False, 0, None, _revcom,
                                 pam_valid_fn=pam_ok)
        self.assertTrue(pam_ok(r["seq"]))
        self.assertEqual("".join(r["seq"]), "AAGT")   # valid-PAM (mm1) beats invalid-PAM (mm0)

    def test_unproductive_column_left_at_reference(self):
        # an alt that only ADDS a mismatch is not chosen; its carriers do not leak in.
        guide = "ACGTACGT"
        ref = "ACGTACGT"           # already a perfect match
        cols = [{"pos_c": 3, "candidates": [{"alt": "A", "carriers": {"sX"}, "info": ["r", "0.1", "s"]}]}]
        got = te.greedy_worst_case(cols, ref, ref, guide, False, 0, None, _revcom)
        self.assertEqual("".join(got["seq"]), "ACGTACGT")
        self.assertEqual(got["carriers"], set())
        self.assertEqual(got["info"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
