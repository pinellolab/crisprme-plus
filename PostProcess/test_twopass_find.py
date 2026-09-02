#!/usr/bin/env python3
"""Tests for twopass_find (Pass-1 collapsed-IUPAC min-edit `D`).

The load-bearing guarantee (design doc section 5): `D` is a LOWER BOUND on the
true edit distance of EVERY concrete haplotype consistent with the window, so
`D <= k` never misses an in-budget carried off-target (for SNP + bulge windows).
This is brute-forced here: for random guides + random allele-set windows, for
EVERY concrete haplotype H in the window, assert `D <= d(guide, H)`.

`d(guide, H)` is computed by the SAME DP with singleton allele sets, so a single
implementation defines both the bound and the ground truth. STDLIB only.
"""
import itertools
import os
import random
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import twopass_find as tp  # noqa: E402


def _singletons(seq):
    return [frozenset(c) for c in seq]


class TestBasics(unittest.TestCase):
    def test_perfect_match(self):
        g = "ACGTACGTAC"
        self.assertEqual(tp.collapsed_min_edit(g, _singletons(g), 0, 0)[0], 0)

    def test_one_mismatch(self):
        g = "ACGTACGTAC"
        w = _singletons("ACGTACGTAT")  # last base differs (C->T)
        self.assertEqual(tp.collapsed_min_edit(g, w, 0, 0)[0], 1)

    def test_fully_ambiguous_window_is_zero(self):
        # the "NNNN is EASY" case: a perfect protospacer is formable -> D=0.
        g = "ACGTACGTACGTACGTACGT"
        self.assertEqual(tp.window_min_edit(g, "N" * 20, 0, 0)[0], 0)

    def test_dna_bulge(self):
        # window one base longer than the guide -> exactly one DNA bulge, D=1.
        g = "ACGTACGT"
        w = _singletons("ACGTAACGT")  # inserted 'A' after pos 5
        D, mm, db, rb = tp.collapsed_min_edit(g, w, max_bdna=1, max_brna=0)
        self.assertEqual((D, db, rb), (1, 1, 0))
        self.assertEqual(mm, 0)

    def test_rna_bulge(self):
        # window one base shorter than the guide -> one RNA bulge, D=1.
        g = "ACGTAACGT"
        w = _singletons("ACGTACGT")
        D, mm, db, rb = tp.collapsed_min_edit(g, w, max_bdna=0, max_brna=1)
        self.assertEqual((D, rb), (1, 1))

    def test_bulge_caps_exceeded(self):
        g = "ACGTACGT"
        w = _singletons("ACGTACGTAAA")  # needs 3 DNA bulges
        self.assertIsNone(tp.collapsed_min_edit(g, w, max_bdna=2, max_brna=0))

    def test_iupac_decode(self):
        # R = A/G : a guide 'A' or 'G' at an R position matches (0), 'C'/'T' don't.
        self.assertEqual(tp.window_min_edit("A", "R", 0, 0)[0], 0)
        self.assertEqual(tp.window_min_edit("G", "R", 0, 0)[0], 0)
        self.assertEqual(tp.window_min_edit("C", "R", 0, 0)[0], 1)

    def test_min_over_windows_picks_best(self):
        g = "ACGTACGT"
        # window 0: 2 mismatches; window 1 (a fake-indel contig): perfect
        r = tp.min_edit_over_windows(g, ["ACGTAAAA", "ACGTACGT"], 0, 0)
        self.assertEqual((r[0], r[4]), (0, 1))  # D=0 from window index 1


class TestLowerBoundBruteForce(unittest.TestCase):
    """D must be <= the concrete edit of EVERY haplotype in the window."""

    def _check(self, max_bdna, max_brna, dlen):
        rnd = random.Random((max_bdna, max_brna, dlen).__hash__())
        for _ in range(1500):
            L = rnd.randint(4, 7)
            guide = "".join(rnd.choice("ACGT") for _ in range(L))
            # window length may differ from guide by <= the bulge budget
            W = max(1, L + rnd.randint(-max_brna, max_bdna))
            sets = [frozenset(rnd.sample("ACGT", rnd.randint(1, 3))) for _ in range(W)]
            D = tp.collapsed_min_edit(guide, sets, max_bdna, max_brna)
            if D is None:
                # if the collapsed (most permissive) alignment can't fit the caps,
                # no concrete haplotype can either -> nothing to check.
                for combo in itertools.product(*sets):
                    self.assertIsNone(
                        tp.collapsed_min_edit(guide, _singletons(combo),
                                              max_bdna, max_brna))
                continue
            for combo in itertools.product(*sets):
                dH = tp.collapsed_min_edit(guide, _singletons(combo),
                                           max_bdna, max_brna)
                self.assertIsNotNone(dH)  # H is a concrete subset -> also alignable
                self.assertLessEqual(
                    D[0], dH[0],
                    f"LOWER-BOUND VIOLATION D={D[0]} > d_H={dH[0]} "
                    f"guide={guide} sets={[''.join(sorted(s)) for s in sets]} H={''.join(combo)}")

    def test_no_bulges(self):
        self._check(0, 0, 0)

    def test_with_bulges(self):
        self._check(1, 1, 1)

    def test_dna_bulge_only(self):
        self._check(2, 0, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
