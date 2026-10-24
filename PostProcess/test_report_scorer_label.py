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
        self.assertEqual(g.scorer_thresholds(), (0.5, 0.2))
        # Section-4 tiers add the sensitive operating point (benchmark 95%-recall default,
        # cb_operating_point) after the >= cuts; file key stays crispr_bulge_gt0 (stable id).
        import cb_operating_point as cbop
        specs = g._crispr_bulge_tier_specs()
        self.assertEqual([s[0] for s in specs],
                         ["&ge; 0.5", "&ge; 0.2", cbop.sensitive_label_html()])
        self.assertEqual(specs[-1][1], "crispr_bulge_gt0")
        self.assertEqual(g._tier_filename("crispr_bulge_gt0"), "crispr_bulge_gt_0.tsv")
        # the sensitive tier predicate is ">= threshold", not "> 0"
        thr = cbop.sensitive_threshold()
        self.assertTrue(specs[-1][2](thr))            # at threshold -> included
        self.assertFalse(specs[-1][2](thr / 2.0))     # below -> excluded

    def test_sensitive_tier_file_contents_predicate(self):
        """The crispr_bulge_gt0 tier predicate -- which selects the rows written to the
        sensitive-tier download crispr_bulge_gt_0.tsv -- keeps EXACTLY the sites with
        CRISPR-Bulge >= the 95%-recall sensitive threshold, NOT every score > 0. A
        sub-threshold noise-floor score (and NaN / the -1 N/A sentinel) must be EXCLUDED
        from the tier file. Pins the file CONTENTS, not just the scalar boundary."""
        import pandas as pd
        import cb_operating_point as cbop
        thr = cbop.sensitive_threshold()
        label, key, pred = g._crispr_bulge_tier_specs()[-1]
        self.assertEqual(key, "crispr_bulge_gt0")           # stable download file key
        s = pd.Series([1.0, 0.5, thr, thr / 2.0, 1e-9, 0.0, -1.0, float("nan")])
        self.assertEqual(list(pred(s)),
                         [True, True, True, False, False, False, False, False])

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
