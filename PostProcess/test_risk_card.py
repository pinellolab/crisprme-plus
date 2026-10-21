#!/usr/bin/env python3
"""Tests for the off-target risk card + CSS-only modal (risk_card.py)."""
import os
import unittest

import numpy as np

import risk_card
import offtarget_calibration as oc

_CALDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration")


def _scores(n=5000, a=0.3, b=8, seed=0):
    rng = np.random.default_rng(seed)
    return list(rng.beta(a, b, n))


class TestRiskCard(unittest.TestCase):
    def setUp(self):
        oc._CACHE.clear()
        oc.set_artifact_dir(_CALDIR)

    def tearDown(self):
        oc._CACHE.clear()
        oc.set_artifact_dir(_CALDIR)

    def test_card_structure_and_self_contained(self):
        card, modal = risk_card.render_risk_card(
            {"guide_display": "ACTGAAATCTGTAAGCAGGC", "nuclease": "SpCas9", "pam": "NRG"},
            _scores(seed=1), _scores(a=0.3, b=12, seed=2), n_perfect=2)
        self.assertIn('class="risk-card"', card)
        self.assertIn("OFF-TARGET RISK", card)
        self.assertIn('href="#risk-card-modal"', card)       # CSS-:target modal link
        self.assertIn("not unique", card)                    # on-target warning (n_perfect>=2)
        self.assertIn('id="risk-card-modal"', modal)
        # self-contained: no JavaScript anywhere in card or modal
        self.assertNotIn("<script", (card + modal).lower())
        # plain header rows only (never <thead>/<tbody>/<th>, which would hijack the
        # report's top-1000 table parser)
        for frag in (card, modal):
            self.assertNotIn("<thead", frag.lower())
            self.assertNotIn("<tbody", frag.lower())
            self.assertNotIn("<th>", frag.lower())

    def test_no_warning_when_on_target_unique(self):
        card, _ = risk_card.render_risk_card({"guide_display": "G"}, _scores(), None, n_perfect=1)
        self.assertNotIn("not unique", card)

    def test_prevalence_note_qualitative_no_arbitrary_ratios(self):
        # the prior-shift caveat is present but carries NO arbitrary ratios/counts
        # (we do not know a deployed guide's own proneness, so we don't assert x1/3/x3)
        card, _ = risk_card.render_risk_card({"guide_display": "G"}, _scores(), None, n_perfect=0)
        self.assertIn("off-target-prone", card.lower())
        self.assertIn("shift this count substantially", card)
        for banned in ("one-third as prone", "three times as prone",
                       "one-third", "three times", "&rarr; ~"):
            self.assertNotIn(banned, card)

    def test_graceful_degradation_missing_artifact(self):
        oc.set_artifact_dir("/nonexistent/calibration/dir")
        oc._CACHE.clear()
        card, modal = risk_card.render_risk_card({"guide_display": "G"}, _scores(), None)
        self.assertIn("OFF-TARGET RISK", card)   # still renders (counts-only)
        self.assertIn("&mdash;", card)           # expected-true column falls back to em-dash
        self.assertNotIn("<script", (card + modal).lower())

    def test_artifact_exposes_ind_quantities(self):
        for scorer in ("cfd", "crispr_bulge"):
            s = oc.artifact_summary(scorer)
            self.assertIsNotNone(s["aupr"]["point"])
            self.assertEqual(len(s["aupr"]["ci"]), 2)
            cov = s["loo_total_count_coverage_nominal90"]["ge0.2"]
            self.assertGreaterEqual(cov, 0.85)   # single-guide coverage ~>=90%
            self.assertIn("ge0.2", s["phi_bands"])

    def test_aggregate_uses_own_band_not_summed_subbands(self):
        # ge0.2 is computed on the pooled band, so its interval need NOT equal the
        # sum of the sub-band intervals (the floor can be below the sub-band floors).
        sc = _scores(n=20000, seed=5)
        res = {r["band"]: r for r in oc.expected_counts(sc, "cfd", n_boot=400, seed=7)}
        self.assertIn("ge0.2", res)
        self.assertIn("0.2-0.5", res)
        self.assertIn("0.5-1", res)
        # point of ge0.2 ~= point(0.2-0.5) + point(0.5-1)
        self.assertAlmostEqual(
            res["ge0.2"]["point"], res["0.2-0.5"]["point"] + res["0.5-1"]["point"], delta=0.5)


if __name__ == "__main__":
    unittest.main()
