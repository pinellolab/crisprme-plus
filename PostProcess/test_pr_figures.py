#!/usr/bin/env python3
"""The PR figures the risk card embeds are reproducible, and say what they are.

Both PNGs under ``PostProcess/calibration/`` used to be produced by scripts kept
outside the repository. That is how two defects reached a collaborator-facing
report: ``pr_curves.png`` drew its CRISPR-Bulge precision from a 300,000-negative
sample of 2.7 M negatives, and ``pr_bulge_stratified.png`` called its CRISPR-Bulge
panel "(in-domain)" while filtering on a benchmark column that is true for every
row, so the panel silently included the 6 guides CRISPR-Bulge was trained on.

``scripts/render_pr_figures.py`` now generates both from a committed curve-data
file. These tests check the chain end to end:

  * every number each figure prints is recomputed here from the committed curve
    data (no plotting needed), so data and figure can never disagree;
  * the CRISPR-Bulge panels carry held-out counts, the CFD panels carry the full
    benchmark, and the titles say which;
  * ``pr_figures.sha256`` ties both PNGs *and* the curve data together, so
    changing either without re-rendering fails;
  * when the local matplotlib/freetype matches the one that rendered the
    committed PNGs, re-rendering reproduces them byte for byte;
  * ``risk_card`` still finds and embeds both PNGs.
"""
import base64
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_CALDIR = os.path.join(_HERE, "calibration")
_REPO = os.path.dirname(_HERE)
_SCRIPT = os.path.join(_REPO, "scripts", "render_pr_figures.py")
if os.path.join(_REPO, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(_REPO, "scripts"))

import render_pr_figures as rpf  # noqa: E402

_DATA = os.path.join(_CALDIR, "pr_figure_data.npz")
_PNGS = ("pr_curves.png", "pr_bulge_stratified.png")


def _sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _manifest():
    out, head = {}, ""
    with open(os.path.join(_CALDIR, rpf.MANIFEST)) as fh:
        for line in fh:
            if line.startswith("#"):
                head = line.strip()
                continue
            digest, name = line.split()
            out[name] = digest
    return out, head


class TestCurveData(unittest.TestCase):
    """Every number the two figures print, recomputed from the committed data."""

    @classmethod
    def setUpClass(cls):
        cls.d = np.load(_DATA, allow_pickle=False)
        cls.meta = json.loads(str(cls.d["meta"]))

    def test_average_precisions_are_what_the_polylines_say(self):
        # The stored AP is the number printed in each panel title/legend; it must
        # be the average precision of the stored polyline, not a leftover value.
        for key in ("cfd", "cb", "scfd_all", "scfd_b0", "scfd_b1",
                    "scb_all", "scb_b0", "scb_b1"):
            with self.subTest(curve=key):
                ap = rpf.average_precision(self.d[f"{key}_rec"], self.d[f"{key}_prec"])
                self.assertAlmostEqual(ap, float(self.d[f"{key}_ap"]), places=12)

    def test_curves_are_monotone_in_recall_and_bounded(self):
        for key in ("cfd", "cb", "scfd_all", "scfd_b0", "scfd_b1",
                    "scb_all", "scb_b0", "scb_b1"):
            with self.subTest(curve=key):
                rec, prec = self.d[f"{key}_rec"], self.d[f"{key}_prec"]
                self.assertTrue(np.all(np.diff(rec) >= 0), "recall must not go back")
                self.assertTrue(np.all((rec >= 0) & (rec <= 1)))
                self.assertTrue(np.all((prec >= 0) & (prec <= 1)))
                self.assertAlmostEqual(float(rec[-1]), 1.0, places=9)

    def test_crispr_bulge_panels_are_held_out_and_cfd_panels_are_not(self):
        # The defect this file exists to prevent: a CRISPR-Bulge panel drawn over
        # the 6 training guides. CFD has no training set, so its panels keep the
        # whole benchmark -- the support the shipped cfd_calibration.npz uses.
        self.assertEqual(int(self.d["cb_pos"]), 287)
        self.assertEqual(int(self.d["scb_all_pos"]), 287)
        self.assertEqual(int(self.d["scb_guides"]), 40)
        self.assertEqual(int(self.d["cfd_pos"]), 359)
        self.assertEqual(int(self.d["scfd_all_pos"]), 359)
        self.assertEqual(int(self.d["scfd_guides"]), 46)
        self.assertLess(int(self.d["cb_n"]), int(self.d["cfd_n"]))
        self.assertIn("cb_heldout", self.meta["cb_support"])

    def test_bulge_strata_partition_their_panel(self):
        for panel in ("scfd", "scb"):
            with self.subTest(panel=panel):
                self.assertEqual(int(self.d[f"{panel}_b0_pos"]) + int(self.d[f"{panel}_b1_pos"]),
                                 int(self.d[f"{panel}_all_pos"]))
                self.assertEqual(int(self.d[f"{panel}_b0_n"]) + int(self.d[f"{panel}_b1_n"]),
                                 int(self.d[f"{panel}_all_n"]))

    def test_no_sampling_and_no_rng(self):
        # Defect A was a 300k-negative weighted sample; precision at CB >= 0.5 rested
        # on 53 negatives in the whole space and moved by 0.25 between samples.
        self.assertEqual(self.meta["sampling"], "none -- every negative is used")
        with open(_SCRIPT) as fh:
            src = fh.read()
        self.assertNotIn("import random", src)
        self.assertNotIn("NEG_SAMPLE", src)

    def test_marked_points_are_the_reports_own_operating_points(self):
        import cb_operating_point as cbop
        marks = json.loads(str(self.d["marks"]))
        thresholds = {(panel, round(t, 12)) for panel, t, _l, _p, _r, _tp, _fp in marks}
        self.assertEqual(thresholds, {("cfd", 0.2), ("cfd", 0.5), ("cb", 0.2), ("cb", 0.5),
                                      ("cb", round(cbop.sensitive_threshold(), 12))})
        for panel, t, _lab, prec, rec, tp, fp in marks:
            with self.subTest(panel=panel, t=t):
                n_pos = int(self.d[f"{panel}_pos"])
                self.assertAlmostEqual(prec, tp / (tp + fp), places=12)
                self.assertAlmostEqual(rec, tp / n_pos, places=12)
                # The dot must sit ON the drawn polyline. A threshold that falls
                # between two observed scores lands inside a vertical run (constant
                # recall, falling precision), whose interior points the stored curve
                # drops as collinear -- so check the segment, not a vertex.
                rc, pc = self.d[f"{panel}_rec"], self.d[f"{panel}_prec"]
                on = np.isclose(rc, rec, rtol=0, atol=1e-12)
                self.assertTrue(on.any(), "marked recall is not on the curve")
                self.assertLessEqual(prec, float(pc[on].max()) + 1e-12)
                self.assertGreaterEqual(prec, float(pc[on].min()) - 1e-12)

    def test_sensitive_point_still_recovers_95_percent(self):
        import cb_operating_point as cbop
        marks = json.loads(str(self.d["marks"]))
        sens = [m for m in marks if m[0] == "cb"
                and abs(m[1] - cbop.sensitive_threshold()) < 1e-18][0]
        self.assertGreaterEqual(sens[4], 0.94)       # recall
        self.assertLessEqual(sens[4], 0.96)
        self.assertIn("recall", sens[2])

    def test_input_provenance_is_recorded(self):
        for key in ("benchmark_md5", "scored_md5", "generator", "ap_recipe"):
            self.assertTrue(self.meta.get(key), f"meta.{key} missing")
        self.assertEqual(len(self.meta["benchmark_md5"]), 32)
        self.assertEqual(self.meta["benchmark_rows"], 2940588)


class TestCurveMaths(unittest.TestCase):
    """The pure functions the committed numbers come out of."""

    def test_average_precision_matches_a_hand_worked_example(self):
        rec, prec, ap, n_pos, n = rpf.pr_curve([0.9, 0.8, 0.7, 0.6], [1, 0, 1, 0])
        self.assertEqual((n_pos, n), (2, 4))
        self.assertAlmostEqual(ap, 0.5 * 1.0 + 0.5 * (2 / 3), places=12)

    def test_ties_do_not_depend_on_row_order(self):
        a = rpf.pr_curve([0.5, 0.5, 0.5, 0.1], [1, 0, 0, 1])[2]
        b = rpf.pr_curve([0.5, 0.5, 0.5, 0.1], [0, 0, 1, 1])[2]
        self.assertAlmostEqual(a, b, places=12)

    def test_decimation_preserves_the_polyline_and_the_average_precision(self):
        rng = np.random.default_rng(0)
        s = rng.random(5000)
        y = (rng.random(5000) < 0.02).astype(int)
        y[0] = 1
        rec, prec, ap, _, _ = rpf.pr_curve(s, y)
        drec, dprec = rpf.decimate(rec, prec)
        self.assertLess(drec.size, rec.size)
        self.assertAlmostEqual(rpf.average_precision(drec, dprec), ap, places=12)
        # dropped points lie exactly on the kept segments (same recall, so the
        # rendered line is identical)
        self.assertEqual(set(np.round(drec, 12)), set(np.round(rec, 12)))

    def test_point_at_counts_the_list_a_threshold_would_produce(self):
        prec, rec, tp, fp = rpf.point_at([0.9, 0.4, 0.2, 0.1], [1, 1, 0, 1], 0.2)
        self.assertEqual((tp, fp), (2, 1))
        self.assertAlmostEqual(prec, 2 / 3, places=12)
        self.assertAlmostEqual(rec, 2 / 3, places=12)

    def test_compact_sci_prints_the_threshold_the_report_quotes(self):
        self.assertEqual(rpf.compact_sci(1.9e-4), "1.9e-4")
        self.assertEqual(rpf.compact_sci(2.7e-3), "2.7e-3")


class TestCommittedFigures(unittest.TestCase):
    """The committed PNGs are the ones the committed data produces."""

    def test_manifest_covers_the_data_and_both_pngs(self):
        digests, _head = _manifest()
        self.assertEqual(set(digests), {os.path.basename(_DATA)} | set(_PNGS))
        for name, digest in digests.items():
            with self.subTest(name=name):
                self.assertEqual(_sha256(os.path.join(_CALDIR, name)), digest,
                                 f"{name} changed without re-running "
                                 "scripts/render_pr_figures.py render")

    def test_rerender_reproduces_the_committed_pngs(self):
        _digests, head = _manifest()
        if rpf.toolchain() not in head:
            self.skipTest(f"committed PNGs were rendered with '{head}'; "
                          f"this machine has '{rpf.toolchain()}'. The numeric checks "
                          "above still cover the content.")
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([sys.executable, _SCRIPT, "render", "--data", _DATA,
                            "--outdir", tmp], check=True, capture_output=True)
            for name in _PNGS:
                with self.subTest(png=name):
                    self.assertEqual(_sha256(os.path.join(tmp, name)),
                                     _sha256(os.path.join(_CALDIR, name)))

    def test_risk_card_embeds_both_figures(self):
        import risk_card
        for uri, name in ((risk_card._pr_curve_uri(), "pr_curves.png"),
                          (risk_card._png_uri("pr_bulge_stratified.png"),
                           "pr_bulge_stratified.png")):
            with self.subTest(png=name):
                self.assertTrue(uri and uri.startswith("data:image/png;base64,"))
                blob = base64.b64decode(uri.split(",", 1)[1])
                with open(os.path.join(_CALDIR, name), "rb") as fh:
                    self.assertEqual(blob, fh.read())

    def test_panel_titles_name_their_support_not_in_domain(self):
        # "(in-domain)" was false: in_domain is n_bulge <= 1, true for every row of
        # a <=1-gap benchmark, so the filter dropped nothing and the panel included
        # the 6 training guides. The title must state the support it really has.
        d = np.load(_DATA, allow_pickle=False)
        cfd_title, cb_title = rpf.stratified_titles(d)
        self.assertNotIn("in-domain", cb_title)
        self.assertIn("40 held-out guides", cb_title)
        self.assertIn("all 46 guides", cfd_title)
        self.assertIn("held-out", rpf.overall_titles(d)[1])

    def test_legend_counts_are_the_held_out_counts(self):
        d = np.load(_DATA, allow_pickle=False)
        legend = rpf.stratified_legend(d, "scb")
        self.assertIn("pos 287/", legend[0])
        self.assertNotIn("pos 359/", " ".join(legend))
        self.assertIn("pos 359/", rpf.stratified_legend(d, "scfd")[0])


_BENCH_HEAD_OLD = ("guide\tofftarget_seq\tchrom\tpos\tstrand\tn_mismatch\tn_bulge\t"
                   "bulge_type\tlabel\tin_domain\tsource\tcb_heldout\tcrrna_aligned\n")
_BENCH_ROWS_OLD = ("G1\tAAC\tchr1\t10\t+\t1\t0\tX0\t1\t1\tTsai_2015\t1\tAAA\n"
                   "G2\tTTC\tchr2\t20\t-\t2\t1\tRNA\t0\t1\tChen_2017\t0\tTT-A\n")
_SCORED = ("guide\tofftarget_seq\tn_mismatch\tn_bulge\tin_domain\tlabel\tcfd_score\t"
           "crispr_bulge_score\nG1\tAAC\t1\t0\t1\t1\t0.5\t0.9\n"
           "G2\tTTC\t2\t1\t1\t0\t0.1\t0.2\n")


def _drop_in_domain(text):
    """The layout scripts/build_offtarget_benchmark.py now writes."""
    return "\n".join("\t".join(c for i, c in enumerate(line.split("\t")) if i != 9)
                     for line in text.splitlines()) + "\n"


class TestCalibrationLoaderSurvivesTheColumnRemoval(unittest.TestCase):
    """fit_calibration.py used to read cb_heldout by index, one past in_domain."""

    @staticmethod
    def _load(tmp, bench_text):
        import fit_calibration
        b = os.path.join(tmp, "benchmark.tsv")
        s = os.path.join(tmp, "scored.tsv")
        with open(b, "w") as fh:
            fh.write(bench_text)
        with open(s, "w") as fh:
            fh.write(_SCORED)
        return fit_calibration.load(s, b)

    def test_same_values_with_and_without_the_dropped_column(self):
        old_text = _BENCH_HEAD_OLD + _BENCH_ROWS_OLD
        with tempfile.TemporaryDirectory() as tmp:
            old = self._load(tmp, old_text)
        with tempfile.TemporaryDirectory() as tmp:
            new = self._load(tmp, _drop_in_domain(old_text))
        for a, b in zip(old, new):
            np.testing.assert_array_equal(a, b)
        # and it is still the leakage flag that selects the CRISPR-Bulge fit
        np.testing.assert_array_equal(new[4], np.array([True, False]))

    def test_misaligned_inputs_are_refused_not_silently_zipped(self):
        swapped = _BENCH_HEAD_OLD + "".join(reversed(_BENCH_ROWS_OLD.splitlines(True)))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit):
                self._load(tmp, swapped)


if __name__ == "__main__":
    unittest.main(verbosity=2)
