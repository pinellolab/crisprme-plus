"""Scorer-aware report column: the second-score column's display label + tier
thresholds follow the active scorer (CRISTA | CRISPR-Bulge), read from .Params.txt.
Physical columns stay '*_(highest_CRISTA)'; only display switches."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import generate_report as g  # imports pandas + matplotlib (present in unit CI)
    HAVE = True
except Exception:
    HAVE = False


@unittest.skipUnless(HAVE, "generate_report deps (pandas/matplotlib) absent")
class TestScorerLabel(unittest.TestCase):
    def tearDown(self):
        g._set_active_scorer("crista")  # restore default for other tests

    def test_labels_and_thresholds_by_scorer(self):
        g._set_active_scorer("crista")
        self.assertEqual(g.scorer_label(), "CRISTA")
        self.assertEqual(g.scorer_thresholds(), (0.6, 0.4, 0.2))
        g._set_active_scorer("crispr-bulge")
        self.assertEqual(g.scorer_label(), "CRISPR-Bulge")
        self.assertEqual(g.scorer_thresholds(), (0.5, 0.2, 0.1))

    def test_unknown_scorer_falls_back_to_crista(self):
        g._set_active_scorer("nonsense")
        self.assertEqual(g.scorer_label(), "CRISTA")
        self.assertEqual(g.scorer_thresholds(), (0.6, 0.4, 0.2))
        g._set_active_scorer(None)
        self.assertEqual(g.scorer_label(), "CRISTA")

    def test_curated_header_relabeled(self):
        g._set_active_scorer("crispr-bulge")
        h = g.curated_headers(True)
        self.assertIn("CRISPR-Bulge", h)
        self.assertNotIn("CRISTA", h)
        g._set_active_scorer("crista")
        self.assertIn("CRISTA", g.curated_headers(True))
        self.assertNotIn("CRISPR-Bulge", g.curated_headers(True))

    def test_ml_column_dropped_when_not_computed(self):
        # when the ML score wasn't computed, neither label appears (column dropped)
        g._set_active_scorer("crispr-bulge")
        h = g.curated_headers(False)
        self.assertNotIn("CRISPR-Bulge", h)
        self.assertNotIn("CRISTA", h)


if __name__ == "__main__":
    unittest.main()
