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

    def test_cb_sensitive_band_counts_down_to_threshold_not_just_ge02(self):
        """The CRISPR-Bulge "sensitive" band counts every site AT/ABOVE the benchmark
        95%-recall operating point (cb_operating_point default 1.9e-4) -- not just >=0.2,
        and not every CB>0 (scores in the sub-threshold noise floor are excluded).
        Regression for the bug where render_risk_card pre-filtered CB to >=0.2 before
        counting, so the band collapsed to the >=0.2 count. CFD stays >=0.2."""
        import re, html
        import cb_operating_point as cbop
        sens = f"CRISPR-Bulge {cbop.sensitive_label_text()}, sensitive"
        # 2 at >=0.5, 2 in [0.2,0.5), 5 in [threshold,0.2): sensitive=9, ge0.2=4, 0.5-1=2.
        # Plus one value BELOW the 95%-recall threshold -> excluded from the sensitive band.
        cb = [0.87, 0.55, 0.33, 0.21, 0.13, 0.08, 0.04, 0.02, 0.005,
              cbop.sensitive_threshold() / 2.0]
        card, _ = risk_card.render_risk_card(
            {"guide_display": "G", "nuclease": "SpCas9", "pam": "NRG"},
            cfd_scores=[0.6, 0.3], cb_scores=cb, n_perfect=0)
        counts = {}
        for row in re.findall(r"<tr[^>]*>(.*?)</tr>", card, re.S):
            cells = [html.unescape(re.sub("<[^>]+>", "", c)).strip()
                     for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)]
            cells = [c for c in cells if c != ""]
            if cells and cells[0].startswith("CRISPR-Bulge") and len(cells) > 1:
                counts[cells[0]] = int(cells[1].replace(",", ""))
        self.assertEqual(counts.get(sens), 9)                 # >= threshold (sub-threshold one excluded)
        self.assertEqual(counts.get("CRISPR-Bulge ≥ 0.2"), 4)  # the >=0.2 subset
        self.assertEqual(counts.get("CRISPR-Bulge ≥ 0.5"), 2)
        # the bug made the sensitive band == ge0.2; assert they differ when sub-0.2 scores exist
        self.assertGreater(counts[sens], counts["CRISPR-Bulge ≥ 0.2"])

    def test_pr_curve_section_caption_names_operating_points(self):
        """The 'how to read this' modal embeds the PR-curve image AND captions the marked
        operating points, so a reader can interpret the curve: CFD >=0.2/>=0.5 and the
        CRISPR-Bulge SENSITIVE (95%-recall) default via cb_operating_point. Regression for
        the previously un-captioned PR image (the reader couldn't tell what the markers were)."""
        import cb_operating_point as cbop
        _, modal = risk_card.render_risk_card(
            {"guide_display": "G", "nuclease": "SpCas9", "pam": "NRG"},
            cfd_scores=[0.6, 0.3], cb_scores=[0.8, 0.3, 0.1], n_perfect=0)
        self.assertIn("rc-pr", modal)                       # PR image embedded (ships in-repo)
        self.assertIn("operating points", modal)            # caption present
        self.assertIn(cbop.sensitive_label_html(), modal)   # names the sensitive default
        self.assertIn("95%-recall threshold", modal)        # explains what the sensitive point is


class TestModalSaysWhatWasMeasured(unittest.TestCase):
    """The explanation page must describe the benchmark it was measured on.

    Every claim below was contradicted by a re-measurement of the raw benchmark
    (2,940,588 rows / 359 positives / 46 guides; CRISPR-Bulge held-out 2,703,869 /
    287 / 40) and is pinned here against the shipped artifacts so it cannot come back.
    """

    def setUp(self):
        oc._CACHE.clear()
        oc.set_artifact_dir(_CALDIR)

    def tearDown(self):
        oc._CACHE.clear()
        oc.set_artifact_dir(_CALDIR)

    def _modal(self):
        _card, modal = risk_card.render_risk_card(
            {"guide_display": "G", "nuclease": "SpCas9", "pam": "NRG"},
            cfd_scores=[0.6, 0.3, 0.25], cb_scores=[0.6, 0.3, 0.05])
        return modal

    def test_modal_discloses_the_on_target_positives(self):
        """The benchmark counts each guide's own on-target as a positive; the card's score
        arrays do not (generate_report builds them with `_mmb > 0`, "exclude perfect"). The
        page must say so where the numbers are, with the deployment-matched values beside
        them: measured on the benchmark with the 46/40 perfect-match rows dropped, CFD's
        trapezoid AUPR is 0.025 (from 0.168), CRISPR-Bulge's tie-grouped full-precision AP
        0.521 (from 0.622), and CFD's >=0.5 precision 0.048 (from 0.085)."""
        modal = self._modal()
        self.assertIn("on-target as a positive", modal)
        self.assertIn("46 of CFD’s 359", modal)
        self.assertIn("40 of CRISPR-Bulge’s 287", modal)
        for matched in ("0.025", "0.521", "0.511", "0.048"):
            self.assertIn(matched, modal)
        # and it must not be buried: the disclosure sits with the AUPR table, before the
        # figure whose panel titles print the same numbers.
        self.assertLess(modal.index("on-target as a positive"), modal.index("rc-pr"))

    def test_modal_does_not_call_the_cfd_column_held_out(self):
        """CFD has no training set: its support is all 46 guides, the 6 CRISPR-Bulge
        trained on included, so nothing is held out for it. The table header, the note and
        the figure caption said "held-out" of both scorers."""
        modal = self._modal()
        self.assertNotIn("held-out AUPR", modal)
        self.assertIn("external-benchmark AUPR", modal)
        self.assertIn("nothing is held out for it", modal)
        self.assertNotIn("on held-out validated off-target sites", modal)

    def test_modal_counts_positives_as_positives_not_sites(self):
        """"CFD 359 sites / 46 guides" said the AUPR rested on 359 rows. 359 and 287 are
        positive counts; the sites are 2,940,588 and 2,703,869, and both come from the
        artifacts that the table's own numbers come from."""
        modal = self._modal()
        self.assertNotIn("359 sites", modal)
        for scorer, pos in (("cfd", 359), ("crispr_bulge", 287)):
            s = oc.artifact_summary(scorer)
            self.assertEqual(int(s["n_positives"]), pos)
            self.assertIn(f"{pos} validated off-targets" if scorer == "cfd" else f"{pos} among",
                          modal)
            self.assertIn(f"{int(s['n_sites']):,}", modal)

    def test_modal_pr_caption_has_the_axes_the_figures_plot(self):
        """The figures plot recall on x, so higher recall moves RIGHT and down. The caption
        said "left-and-down", which reads the trade-off backwards."""
        modal = self._modal()
        self.assertNotIn("left-and-down", modal)
        self.assertIn("right-and-down", modal)
        self.assertIn("recall is the x-axis", modal)

    def test_modal_shows_each_bands_benchmark_support(self):
        """band_counts is in both artifacts and returned by artifact_summary, and nothing
        displayed it: the card gave an expected count and a 90% range for a band without
        saying how many positives, or how many guides, it rests on. Every number in the
        new table must come from the artifact."""
        modal = self._modal()
        self.assertIn("benchmark positives in it", modal)
        self.assertIn("guides carrying one", modal)
        for scorer, spec in (("cfd", risk_card._CFD_ROWS), ("crispr_bulge", risk_card._CB_ROWS)):
            counts = oc.artifact_summary(scorer)["band_counts"]
            for label, disp in spec:
                c = counts[label]
                row = (f'<tr><td>{disp}</td><td>{c["positives"]:,}</td>'
                       f'<td>{c["guides_with_positive"]}</td>'
                       f'<td>{c["guides_contributing"]}</td></tr>')
                with self.subTest(band=disp):
                    self.assertIn(row, modal)
        # the band the support table exists for: 41 positives from 12 of 40 guides
        cb = oc.artifact_summary("crispr_bulge")["band_counts"]["0.2-0.5"]
        self.assertEqual((cb["positives"], cb["guides_with_positive"]), (41, 12))

    def test_modal_states_the_measured_value_of_the_two_cfd_rescue_arms(self):
        """On the 40 held-out guides the two CFD arms add 9,660 and 614 sites to the
        6,076-site sensitive tier and recover 1 of the 14 validated off-targets it misses;
        10 of those 14 carry a bulge and the CFD bulge arm recovers none of them. The page
        asserted the opposite ("rescues the high-CFD sites CRISPR-Bulge misses")."""
        modal = self._modal()
        self.assertNotIn("the high-CFD sites CRISPR-Bulge misses", modal)
        for n in ("9,660", "614", "6,076", "2.7&times;"):
            self.assertIn(n, modal)
        self.assertIn("recover <b>1</b> of the 14", modal)
        self.assertIn("none of the 10 missed sites that carry a bulge", modal)

    def test_modal_recall_growth_is_the_measured_sixteen_fold(self):
        """95% -> 99% recall on the held-out benchmark is 6,076 -> 99,477 sites = 16x at the
        shipped thresholds (16.4x at the exact recall quantiles), not the "roughly 6x" the
        page carried, which came from a worked search."""
        modal = self._modal()
        self.assertNotIn("roughly 6&times; the sites", modal)
        self.assertIn("99%-recall list is <b>16&times;</b>", modal)
        self.assertIn("99,477", modal)

    def test_modal_coverage_claim_matches_the_shipped_artifacts(self):
        """The page claimed "~90% coverage for a single novel guide" while METHODS said the
        interval cannot be validated for a novel clinical guide. What the artifacts measure
        is per-held-out-BENCHMARK-guide coverage, 89-100% over the displayed bands, lowest
        CFD 0.2-0.5. The modal must quote that range, name the band that falls short, and
        keep METHODS' limit on a clinically novel guide."""
        modal = self._modal()
        cov = []
        for scorer, spec in (("cfd", risk_card._CFD_ROWS), ("crispr_bulge", risk_card._CB_ROWS)):
            table = oc.artifact_summary(scorer)["loo_total_count_coverage_nominal90"]
            cov += [table[label] for label, _disp in spec if table.get(label) is not None]
        self.assertEqual(f"{int(round(min(cov) * 100))}–{int(round(max(cov) * 100))}%", "89–100%")
        self.assertIn("89–100% of folds against a nominal 90%", modal)
        self.assertIn("CFD 0.2–0.5, which measures 89%", modal)
        self.assertNotIn("reaches ~90% coverage for a single novel guide", modal)
        self.assertIn("cannot be, validated for a clinically novel guide", modal)

    def test_modal_does_not_claim_the_shipped_map_is_leave_one_guide_out(self):
        """p_point is the ALL-guide isotonic fit (it closes on the benchmark by
        construction); leave-one-guide-out enters only the band errors, phi and the
        coverage check. "built so no guide ever set its own probabilities" described the
        validation, not the map that ships."""
        modal = self._modal()
        self.assertNotIn("no guide ever set its own probabilities", modal)
        self.assertIn("fit on all of those guides at once", modal)
        self.assertIn("not the shipped map itself", modal)


if __name__ == "__main__":
    unittest.main()
