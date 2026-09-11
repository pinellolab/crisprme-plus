#!/usr/bin/env python3
"""Tests for twopass_score (Pass-2 worst-possible per-window scoring).

Verifies the integration of Pass-1 min-edit + exact worst-case CFD, and the
CRISTA best-effort shortlist (ordering toward low-edit haplotypes, truncation
flag, worst = max over the evaluated set). CFD-vs-brute-force is re-checked as an
oracle here too. STDLIB only (a mock crista_fn stands in for the RandomForest).
"""
import itertools
import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import twopass_score as ts  # noqa: E402
import twopass_cfd_exact as cfd  # noqa: E402

MM, PAM = cfd.load_scores()


class TestScoreWindow(unittest.TestCase):
    def test_min_edit_and_cfd_present(self):
        guide = "ACGTACGTACGTACGTACGT"
        sets = [frozenset(b) for b in guide]  # perfect window
        r = ts.score_window(guide, sets, ["GG"], MM, PAM)
        self.assertEqual(r["min_edit"][0], 0)              # perfect -> D=0
        self.assertEqual(r["cfd_worst"], PAM.get("GG", 0.0))  # perfect proto -> max PAM
        self.assertIsNone(r["crista_worst"])              # no crista_fn -> skipped

    def test_cfd_worst_matches_bruteforce(self):
        # short guide so itertools.product is tractable; the scorer's cfd_worst must
        # equal the brute-force max CFD over concrete haplotypes x PAMs.
        import random
        random.seed(5)
        for _ in range(300):
            L = random.randint(4, 7)
            guide = "".join(random.choice("ACGT") for _ in range(L))
            sets = [sorted(set(random.choice("ACGT-") for _ in range(random.randint(1, 3))))
                    for _ in range(L)]
            pam_set = ["GG", "AG", "TG"]
            got = ts.score_window(guide, sets, pam_set, MM, PAM)["cfd_worst"]
            best = 0.0
            for combo in itertools.product(*sets):
                for p in pam_set:
                    v = cfd.cfd_concrete(guide, "".join(combo), p, MM, PAM)
                    best = max(best, v)
            self.assertEqual(got, best, f"guide={guide} sets={sets}")


class TestCristaShortlist(unittest.TestCase):
    def test_low_edit_ordering_puts_guide_base_first(self):
        guide = "ACGT"
        sets = [frozenset("CA"), frozenset("C"), frozenset("GA"), frozenset("T")]
        haps, trunc = ts._iter_low_edit_haplotypes(guide, sets, cap=100)
        self.assertFalse(trunc)
        # the very first haplotype should match the guide at every ambiguous pos it can
        self.assertEqual(haps[0], "ACGT")

    def test_truncation_flag(self):
        guide = "AAAAAA"
        sets = [frozenset("ACGT")] * 6  # 4^6 = 4096 > cap
        haps, trunc = ts._iter_low_edit_haplotypes(guide, sets, cap=10)
        self.assertTrue(trunc)
        self.assertEqual(len(haps), 10)
        self.assertEqual(haps[0], "AAAAAA")  # guide-matching first

    def test_crista_worst_is_max_over_shortlist_and_flags_approx(self):
        guide = "ACGTAC"
        sets = [frozenset("ACGT")] * 6  # huge product -> truncated
        # mock CRISTA: score = number of positions equal to the guide (max at the
        # guide itself, which the low-edit ordering surfaces first).
        def mock_crista(triples):
            return [sum(a == b for a, b in zip(g, d)) for (g, d, p) in triples]
        r = ts.score_window(guide, sets, ["GG"], MM, PAM,
                            crista_fn=mock_crista, crista_cap=20)
        self.assertTrue(r["crista_approx"])                 # truncated -> approximate
        self.assertEqual(r["crista_worst"], len(guide))     # guide hap scores full length
        self.assertEqual(r["crista_n_evaluated"], 20)       # cap x 1 pam

    def test_crista_exhaustive_not_flagged(self):
        guide = "ACG"
        sets = [frozenset("A"), frozenset("C"), frozenset("G")]  # single haplotype
        r = ts.score_window(guide, sets, ["GG"], MM, PAM,
                            crista_fn=lambda tr: [1.0] * len(tr))
        self.assertFalse(r["crista_approx"])
        self.assertEqual(r["crista_n_evaluated"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
