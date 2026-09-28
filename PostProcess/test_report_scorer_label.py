"""Second-score (ML) report column: CRISPR-Bulge is the sole ML scorer, so its
display label and tier thresholds are fixed. The physical column stays
'*_(highest_CRISPR_BULGE)' and is shown only when the ML score was computed."""
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
    def test_label_and_thresholds(self):
        self.assertEqual(g.scorer_label(), "CRISPR-Bulge")
        self.assertEqual(g.scorer_thresholds(), (0.5, 0.2, 0.1))

    def test_set_active_scorer_is_noop(self):
        # retained for call compatibility; any argument keeps CRISPR-Bulge
        g._set_active_scorer("anything")
        self.assertEqual(g.scorer_label(), "CRISPR-Bulge")
        g._set_active_scorer(None)
        self.assertEqual(g.scorer_label(), "CRISPR-Bulge")

    def test_curated_header_present_when_computed(self):
        h = g.curated_headers(True)
        self.assertIn("CRISPR-Bulge", h)

    def test_ml_column_dropped_when_not_computed(self):
        # when the ML score wasn't computed, the label does not appear (column dropped)
        h = g.curated_headers(False)
        self.assertNotIn("CRISPR-Bulge", h)


if __name__ == "__main__":
    unittest.main()
