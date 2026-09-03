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
