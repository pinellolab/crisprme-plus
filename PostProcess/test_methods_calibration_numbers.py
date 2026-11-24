#!/usr/bin/env python3
"""docs/offtarget_risk_summary_METHODS.md must quote the artifacts that SHIP.

The methods document travels with the off-target risk report. Between 2.7.11 and
2.10.1 the CFD score column was refitted and the document was not: its φ values, its
decision-band φ, its per-held-out-guide coverage table, its ≥0.2 aggregate φ and its
global LOGO residual were all still the superseded ``cfd_calibration.npz.pre2711``'s,
while the card beside them rendered the new artifact. Nothing caught it, because no
test read the prose against the ``.npz``.

These tests do. Every number below is re-read from
``PostProcess/calibration/{cfd,crispr_bulge}_calibration.npz`` at run time and matched
against the sentence that quotes it, so a refit that does not update the document
fails here rather than in a reviewer's hands. The handful of quantities that are not
artifact fields -- the on-target-excluded AUPRs, the dominant-guide φ move, the
rescue-arm counts -- are pinned as literals with their provenance, and cross-checked
between the document and ``risk_card`` so the two cannot drift apart.

Run as CI runs it::

    cd PostProcess && python3 test_methods_calibration_numbers.py
"""
import json
import os
import re
import unittest

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_CALDIR = os.path.join(_HERE, "calibration")
_METHODS = os.path.join(_REPO, "docs", "offtarget_risk_summary_METHODS.md")


def _meta(scorer):
    path = os.path.join(_CALDIR, f"{scorer}_calibration.npz")
    return json.loads(str(np.load(path, allow_pickle=False)["meta"]))


def _text():
    with open(_METHODS, encoding="utf-8") as fh:
        return fh.read()


def _one(pattern, text, flags=0):
    """The single match of ``pattern``; fails loudly if the anchor moved."""
    found = re.findall(pattern, text, flags)
    if len(found) != 1:
        raise AssertionError(f"expected exactly one match of {pattern!r}, got {found!r}")
    return found[0]


class TestMethodsQuotesTheShippedArtifacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = _text()
        cls.cfd = _meta("cfd")
        cls.cb = _meta("crispr_bulge")

    def test_section5_phi_values_are_the_shipped_ones(self):
        """§5's low/action-band φ and both decision-band φ come from phi_bands."""
        low = _one(r"CFD φ ≈ ([\d.]+) / ([\d.]+) / ([\d.]+)", self.txt)
        for printed, label in zip(low, ("0-0.1", "0.1-0.2", "0.2-0.5")):
            with self.subTest(band=label):
                self.assertAlmostEqual(float(printed), self.cfd["phi_bands"][label],
                                       delta=0.005)
        decision = float(_one(r"\*\*([\d.]+) in the decision band\*\*", self.txt))
        self.assertAlmostEqual(decision, self.cfd["phi_bands"]["0.5-1"], delta=0.005)
        cb_decision = float(_one(r"CRISPR-Bulge's decision\s+band is \*\*([\d.]+)\*\*",
                                 self.txt))
        self.assertAlmostEqual(cb_decision, self.cb["phi_bands"]["0.5-1"], delta=0.005)

    def test_the_ge02_aggregate_phi_is_the_shipped_one(self):
        """§5 and §7 both quote φ for the pooled ≥0.2 band; both must be the artifact's."""
        agg = self.cfd["phi_bands"]["ge0.2"]
        s5 = float(_one(r"\*\*≥0\.2 aggregate\*\* \(\*\*([\d.]+)\*\*, estimated on all", self.txt))
        s7 = float(_one(r"overdispersion \(shipped φ ≈ \*\*([\d.]+)\*\*", self.txt))
        self.assertAlmostEqual(s5, agg, delta=0.0005)
        self.assertAlmostEqual(s7, agg, delta=0.0005)
        # and §7's dominant-guide move starts from the same value
        lo, hi = _one(r"\*\*φ_≥0\.2 ([\d.]+) → ([\d.]+)\*\*", self.txt)
        self.assertAlmostEqual(float(lo), agg, delta=0.0005)
        # 0.737: φ(ge0.2) re-estimated with the 79-positive dominant guide dropped,
        # measured on the benchmark (0.8255 -> 0.7367). Not an artifact field.
        self.assertAlmostEqual(float(hi), 0.737, delta=0.0005)

    def test_the_coverage_table_row_is_the_shipped_coverage(self):
        """The shipped row of §5's coverage table is
        loo_total_count_coverage_nominal90, in the table's own column order."""
        row = _one(r"\| \*\*two-level \+ overdispersion \(shipped\)\*\* \|(.+)\|\n", self.txt)
        printed = [float(c.strip().strip("*").rstrip("%")) for c in row.split("|") if c.strip()]
        expected = [self.cfd["loo_total_count_coverage_nominal90"]["0.2-0.5"],
                    self.cfd["loo_total_count_coverage_nominal90"]["0.5-1"],
                    self.cfd["loo_total_count_coverage_nominal90"]["ge0.2"],
                    self.cb["loo_total_count_coverage_nominal90"]["ge0.2"]]
        self.assertEqual(len(printed), 4)
        for got, want in zip(printed, expected):
            self.assertAlmostEqual(got, want * 100, delta=0.05)

    def test_the_band_that_falls_below_nominal_is_named_and_not_hidden(self):
        """"All shipped bands reach ≥90%" was true of the pre-2.7.11 artifact only.
        CFD 0.2-0.5 ships at 89.1%; the document must name it, and must not repeat the
        blanket claim anywhere (§5 and §8 both carried it)."""
        cov = dict(self.cfd["loo_total_count_coverage_nominal90"])
        self.assertLess(cov["0.2-0.5"], 0.90)                     # the premise still holds
        want = f'{cov["0.2-0.5"] * 100:.1f}%'
        self.assertTrue(want in self.txt, f"METHODS does not print {want!r}")
        for banned in ("All shipped bands reach ≥90%", "≥90% in all shipped bands"):
            self.assertFalse(banned in self.txt, f"METHODS still claims {banned!r}")

    def test_decision_band_support_is_the_artifacts_band_counts(self):
        """"the decision band rests on only 1-3 guides" matched no artifact. The counts
        the document gives must be band_counts', for both scorers."""
        for meta in (self.cfd, self.cb):
            bc = meta["band_counts"]["0.5-1"]
            want = (f'{bc["positives"]} positives from all '
                    f'{bc["guides_with_positive"]} guides')
            with self.subTest(scorer=meta["scorer"]):
                self.assertTrue(want in self.txt, f"METHODS does not say {want!r}")
        for banned in ("rests on only **1–3 guides", "rests on **1–3 guides"):
            self.assertFalse(banned in self.txt, f"METHODS still claims {banned!r}")
        # the band whose support really is thin, from the artifact
        thin = self.cb["band_counts"]["0.2-0.5"]
        want = (f'{thin["positives"]} positives from {thin["guides_with_positive"]} '
                f'of {self.cb["n_guides"]} guides')
        self.assertTrue(want in self.txt, f"METHODS does not say {want!r}")

    def test_support_counts_and_prevalences_are_the_artifacts(self):
        for meta in (self.cfd, self.cb):
            with self.subTest(scorer=meta["scorer"]):
                for want in (f'{meta["n_sites"]:,}', str(meta["n_positives"])):
                    self.assertTrue(want in self.txt, f"METHODS does not print {want!r}")

    def test_no_superseded_value_is_stated_as_current(self):
        """The pre-2.7.11 numbers may appear only inside a sentence that says they are
        superseded -- they are kept so a reader of the earlier report can see what
        changed, and they must never read as a live claim."""
        superseded = ("1.83 / 1.92 / 0.99", "0.806", "+3.8%", "95 / 93 / 95",
                      "≈6×", "+0.66", "1–3 guides", "358.1")
        marks = ("earlier draft", "pre-2.7.11", "superseded", "not recomputed",
                 "matched no artifact", "no definition we can find")
        # markdown paragraphs: the disclaimer and the value it disowns belong to the
        # same block of prose, which is wrapped over several source lines.
        for para in re.split(r"\n\s*\n", self.txt):
            low = para.lower()
            for bad in superseded:
                if bad in para:
                    with self.subTest(value=bad, para=para[:60].replace("\n", " ")):
                        self.assertTrue(
                            any(m in low for m in marks),
                            f"{bad!r} is stated without saying it is superseded")


class TestMethodsAndModalAgree(unittest.TestCase):
    """The two documents a reviewer reads must carry the same corrected numbers."""

    @classmethod
    def setUpClass(cls):
        cls.txt = _text()

    def test_the_on_target_disclosure_is_in_both_and_matches(self):
        """The benchmark's positives include each guide's own on-target; the report's
        ranked list excludes them. Deployment-matched values, re-measured on the raw
        benchmark with the 46 / 40 perfect-match rows dropped: CFD trapezoid AUPR
        0.025146, CRISPR-Bulge tie-grouped full-precision AP 0.520929, CRISPR-Bulge
        trapezoid on the 3-decimal column 0.510504, CFD t>=0.5 precision 0.047576."""
        import risk_card
        note = risk_card._ONTARGET_NOTE
        for value in ("0.025", "0.521", "0.511", "0.048"):
            with self.subTest(value=value):
                self.assertTrue(value in note, f"the modal drops {value!r}")
                self.assertTrue(value in self.txt, f"METHODS drops {value!r}")
        for name, doc in (("modal", note), ("METHODS", self.txt)):
            self.assertTrue("46 of CFD" in doc, f"{name} does not count the on-targets")

    def test_the_rescue_arms_carry_the_measured_counts_in_both(self):
        """Measured on the 40 held-out guides with CRISPR-Bulge at full precision:
        base tier 6,076 sites / 273 of 287; +9,660 sites and +1 off-target from
        CFD>=0.2 on 0-bulge; +614 and +0 from CFD>=0.5 on 1-bulge; 10 of the 14 missed
        positives carry a bulge and the bulge arm recovers none of them."""
        import risk_card
        _card, modal = risk_card.render_risk_card(
            {"guide_display": "G"}, [0.6, 0.3], [0.6, 0.3, 0.05])
        for value in ("9,660", "614", "6,076"):
            with self.subTest(value=value):
                self.assertTrue(value in modal, f"the modal drops {value!r}")
                self.assertTrue(value in self.txt, f"METHODS drops {value!r}")
        self.assertTrue("recovers none of the 10" in self.txt,
                        "METHODS does not say the bulge arm recovers none of the 10")
        self.assertIn("none of the 10 missed sites that carry a bulge", modal)

    def test_the_recall_growth_ratio_is_the_measured_one_in_both(self):
        """95% -> 99% recall on the held-out benchmark: 6,076 -> 99,477 sites at the
        shipped thresholds (16.4x at the exact recall quantiles, 6,048 -> 99,406)."""
        import risk_card
        _card, modal = risk_card.render_risk_card(
            {"guide_display": "G"}, [0.6, 0.3], [0.6, 0.3, 0.05])
        self.assertIn("99,477", modal)
        self.assertTrue("99,477" in self.txt, "METHODS drops the 99%-recall list size")
        self.assertNotIn("roughly 6&times; the sites", modal)

    def test_the_aupr_gap_is_attributed_to_the_stored_columns_defect_in_both(self):
        """0.613 (table) vs 0.62 (figure) is 0.0090, and 54% of it is the stored CRISPR-Bulge
        column's PAM-internal-bulge scoring defect, 26% the figure's full-precision rescoring
        and 20% its tie-grouped AP rule (trapezoid 0.6130872 / 0.6179856 / 0.6182664,
        tie-grouped 0.6148478 / 0.6197042 / 0.6220800). Both documents said "different score
        precision and AP rule", which names the two minority terms; the same defect is
        described elsewhere in this document, so the two sections contradicted each other on
        a discrepancy the reader can see on the page."""
        import risk_card
        self.assertNotIn("different score", self.txt)
        self.assertNotIn("precision and AP rule", self.txt)
        self.assertNotIn("only the score precision and the AP rule differ",
                         risk_card._CB_AUPR_FIGURE_NOTE)
        for name, doc in (("METHODS", self.txt), ("modal", risk_card._CB_AUPR_FIGURE_NOTE)):
            with self.subTest(doc=name):
                self.assertIn("54%", doc)
                self.assertIn("26%", doc)
                self.assertIn("20%", doc)
                self.assertTrue("inside the 0.56–0.70 CI" in doc or "inside the 90% CI" in doc,
                                f"{name} drops the CI remark")
        self.assertIn("0.6131 → 0.6180", self.txt)      # the trapezoid anchor, matching §2

    def test_the_bulge_ap_is_not_called_the_prevalence_floor_in_either(self):
        """CFD bulge=1 AP 0.009175 against a 5.75e-5 bulge=1 floor is ~160x no-skill, while
        the page calls CFD's 0.168 good at ~1,379x the overall floor: "essentially the
        prevalence floor" / "stays at the floor" contradicted that yardstick. §3.5's
        comparison against CRISPR-Bulge's 0.285 is the defensible framing and is now used
        in all three places."""
        import risk_card
        _card, modal = risk_card.render_risk_card(
            {"guide_display": "G"}, [0.6, 0.3], [0.6, 0.3, 0.05])
        for name, doc in (("METHODS", self.txt), ("modal", modal)):
            with self.subTest(doc=name):
                self.assertNotIn("essentially the prevalence floor", doc)
                self.assertNotIn("at the floor", doc)
                self.assertIn("0.285", doc)
        self.assertEqual(self.txt.count("bulge=1 AP ≈ 0.01 vs CRISPR-Bulge 0.285"), 2)

    def test_the_on_target_refit_moves_the_decision_band_down_not_up(self):
        """The published map prints x1.082 what a map refit without the on-target rows
        prints in CFD ≥0.5 and x1.021 over ≥0.2, so the deployment-matched count is LOWER
        (x0.924, x0.979), not higher. The magnitudes were right and the direction was
        inverted -- in the band a reviewer acts on."""
        self.assertAlmostEqual(1 / 1.082, 0.924, places=3)
        self.assertAlmostEqual(1 / 1.021, 0.979, places=3)
        para = _one(r"(?s)(effect on the card's \*expected counts\*.*?top bins\.)", self.txt)
        self.assertIn("lowers the CFD ≥0.5 expected total by about 8% (×0.92)", para)
        self.assertIn("≥0.2 total by 2% (×0.98)", para)
        self.assertIn("×1.082", para)
        self.assertIn("×1.021", para)
        self.assertNotIn("would print moves ×1.08", para)

    def test_the_benchmark_scope_and_assay_floor_reached_the_page(self):
        """METHODS §2 has both; the modal is what the collaborator reads and had neither."""
        import risk_card
        _card, modal = risk_card.render_risk_card(
            {"guide_display": "G"}, [0.6, 0.3], [0.6, 0.3, 0.05])
        for phrase in ("GUIDE-seq", "assay-detected", "reference-genome",
                       "zero ≥2-gap rows", "22%", "27.5%"):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, modal)
        # the same dominant-guide share the document reports (METHODS wraps the sentence)
        self.assertRegex(self.txt, r"79 of the 359 positives\s+\(79/359 = 22\.0%\)")
        self.assertIn("79 of CFD’s 359 positives (22%)", modal)


if __name__ == "__main__":
    unittest.main(verbosity=2)
