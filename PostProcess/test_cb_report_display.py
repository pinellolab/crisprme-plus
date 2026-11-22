#!/usr/bin/env python3
"""How the report SHOWS a CRISPR-Bulge score: its text and its axis.

A collaborator reviewing a pre-IND off-target assessment could not use the report's
CRISPR-Bulge columns or panels:

  * the curated tables (the in-report top-1000 and every download) printed the score
    with 4 decimals, so the sensitive operating point (1.9e-4, cb_operating_point) fell
    between "0.0001" and "0.0002" and nearly every row of a top-1000 BY CRISPR-Bulge read
    "0.0000" -- the column could not tell a site above the threshold from one below it;
  * the score panels put the score on a LINEAR 0-1 y-axis, where every point of a
    genome-wide run lies on the bottom rule and the threshold cannot be drawn at all.

Now both the curated cells and the two plot scripts take their text and their axis from
cb_operating_point (score_text / log_axis / draw_log_score_axis), so the display and the
threshold they display cannot drift apart. The zeros -- the model's float32 output
underflows, a few percent of the rows of a genome-wide run -- are kept: a log axis cannot
show 0, so they are clipped onto a marked floor one decade below the smallest non-zero
score and named in the legend, never dropped in silence.

CFD is deliberately untouched: it is calibrated, uses the whole 0-1 range, and keeps its
4 decimals and its linear axis.

The SEARCH RESULT TABLES are not touched by any of this: the post-analyses already write
CRISPR-Bulge with 6 significant digits, and this module only formats for display.

Run with:

    cd PostProcess && python3 test_cb_report_display.py
"""

import os
import sys
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    import cb_operating_point as cbop
    import generate_report as gr

    _HAVE_DEPS = True
    _SKIP_REASON = ""
except Exception as exc:  # pragma: no cover - dependency-less CI leg
    _HAVE_DEPS = False
    _SKIP_REASON = f"matplotlib/pandas/report deps unavailable: {exc}"


_SPACER = "CTCTCAGCTGGTACACGGCANNN"
_REF_ALN = "CTCTCAGCTGGTACACGGCAtGG"
_ALT_ALN = "CTCTCAGCTGGTACACGGCAAGG"


def _frame(rows):
    """A minimal highest-CFD/highest-CRISPR_BULGE frame, as the report resolves one."""
    header = {
        "Spacer+PAM": _SPACER,
        "Chromosome": "chr5",
        "Start_coordinate_(highest_CFD)": "400",
        "Strand_(highest_CFD)": "+",
        "Aligned_spacer+PAM_(highest_CFD)": _SPACER,
        "Aligned_protospacer+PAM_REF_(highest_CFD)": _REF_ALN,
        "Aligned_protospacer+PAM_ALT_(highest_CFD)": _ALT_ALN,
        "Mismatches_(highest_CFD)": "2",
        "Bulges_(highest_CFD)": "0",
        "Mismatches+bulges_(highest_CFD)": "2",
        "REF/ALT_origin_(highest_CFD)": "alt",
        "CFD_score_(highest_CFD)": "0.8",
        "CFD_score_REF_(highest_CFD)": "0.2",
        "CRISPR_BULGE_score_(highest_CRISPR_BULGE)": "0.5",
        "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": "0.5",
    }
    return pd.DataFrame([{**header, **row} for row in rows])


def _curated(rows, has_crispr_bulge=True):
    df = _frame(rows)
    cols = gr._resolve(df.columns, list(gr._COLS.keys()))
    return gr.build_curated_frame(df, cols, has_crispr_bulge)


def _axes_of(**kwargs):
    """Run _cfd_style_scatter and return the Axes it drew on (never written to a URI)."""
    seen = {}

    def _capture(fig, dpi=120):
        seen["fig"] = fig
        seen["ax"] = fig.axes[0]
        return "data:,"

    with mock.patch.object(gr, "_fig_to_data_uri", _capture):
        gr._cfd_style_scatter(**kwargs)
    return seen["ax"], seen["fig"]


def _hlines(ax):
    """y of every full-width horizontal line drawn on ax (axhline)."""
    return [line.get_ydata()[0] for line in ax.lines
            if len(set(line.get_ydata())) == 1 and len(line.get_ydata()) == 2]


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestScoreText(unittest.TestCase):
    """cb_operating_point.score_text: the ONE text of a CRISPR-Bulge score."""

    def test_the_threshold_is_distinguishable_from_just_above_and_just_below(self):
        """What the collaborator asked for: at 4 decimals the operating point and the
        scores either side of it are the same string, so the column cannot say which
        side of the threshold a site is on."""
        thr = cbop.sensitive_threshold()
        below, above = thr * 0.9, thr * 1.1
        texts = [cbop.score_text(v) for v in (below, thr, above)]
        self.assertEqual(len(set(texts)), 3, texts)
        # and each still reads back as its own number (6 significant digits)
        for value, text in zip((below, thr, above), texts):
            self.assertAlmostEqual(float(text), value, delta=abs(value) * 1e-6)
        # the order is legible, not only distinguishable
        self.assertEqual([float(t) for t in texts], sorted(float(t) for t in texts))
        # the format it replaces could not do this (all three are "0.0002")
        self.assertEqual(len({f"{v:.4f}" for v in (below, thr, above)}), 1)

    def test_the_threshold_text_is_the_threshold_label_text(self):
        """score_text and the tier/panel labels must print the operating point the same
        way, or a reader compares a cell with a label that does not match it."""
        self.assertEqual(cbop.score_text(cbop.sensitive_threshold()), cbop._thr_str())
        self.assertIn(cbop.score_text(cbop.sensitive_threshold()),
                      cbop.sensitive_label_text())

    def test_zero_one_and_a_denormal_score(self):
        self.assertEqual(cbop.score_text(0.0), "0")        # a true 0 still reads as 0
        self.assertEqual(cbop.score_text(1.0), "1")
        # the smallest non-zero score of a 2.10.0 genome-wide run
        self.assertEqual(cbop.score_text(2.35e-39), "2.35e-39")
        # a mid-range score keeps 6 significant digits, as the post-analyses write it
        self.assertEqual(cbop.score_text(0.03521684), "0.0352168")

    def test_the_format_is_the_one_the_post_analyses_write(self):
        """cbulge_rows / targets_scoring have written "{:.6g}" since 2.8.0: the curated
        cell re-formats the value it reads rather than rounding it a second time."""
        self.assertEqual(cbop.SCORE_FORMAT, "{:.6g}")
        for raw in ("1.06094e-21", "0.319816", "2.35e-39", "0"):
            self.assertEqual(cbop.score_text(float(raw)), raw)


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestCuratedScoreCells(unittest.TestCase):
    """The curated tables (the in-report table and every download share them)."""

    def test_crispr_bulge_cells_place_a_row_against_the_threshold(self):
        thr = cbop.sensitive_threshold()
        frame = _curated([
            {"CRISPR_BULGE_score_(highest_CRISPR_BULGE)": repr(thr * 1.1),
             "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": repr(thr)},
            {"CRISPR_BULGE_score_(highest_CRISPR_BULGE)": repr(thr * 0.9),
             "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": "2.35e-39"},
        ])
        label = gr.scorer_label()
        self.assertEqual(list(frame[label]), ["0.000209", "0.000171"])
        self.assertEqual(list(frame[label + "_ref"]), ["0.00019", "2.35e-39"])
        # every cell is on the right side of the threshold, and says so
        for cell, expected in zip(frame[label], (True, False)):
            self.assertEqual(float(cell) >= thr, expected, cell)

    def test_a_true_zero_reads_zero_and_an_unscored_row_keeps_the_placeholder(self):
        frame = _curated([
            {"CRISPR_BULGE_score_(highest_CRISPR_BULGE)": "0",
             "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": "0.0"},
            # -1 is the pipeline's "this row was not scored" sentinel, NA a missing value
            {"CRISPR_BULGE_score_(highest_CRISPR_BULGE)": "-1",
             "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": "NA"},
        ])
        label = gr.scorer_label()
        self.assertEqual(list(frame[label]), ["0", gr.CURATED_MISSING])
        self.assertEqual(list(frame[label + "_ref"]), ["0", gr.CURATED_MISSING])

    def test_cfd_is_untouched(self):
        """CFD is calibrated and uses the whole 0-1 range: it keeps its 4 decimals."""
        frame = _curated([{"CFD_score_(highest_CFD)": "0.8",
                           "CFD_score_REF_(highest_CFD)": "0.123456"},
                          {"CFD_score_(highest_CFD)": "-1",
                           "CFD_score_REF_(highest_CFD)": "0"}])
        self.assertEqual(list(frame["CFD"]), ["0.8000", gr.CURATED_MISSING])
        self.assertEqual(list(frame["CFD_ref"]), ["0.1235", "0.0000"])


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestRefProtospacerColumn(unittest.TestCase):
    """The REF protospacer alignment, its own column beside the ALT one, so that the
    three lines (spacer, REF protospacer, ALT protospacer) read as an alignment and a
    reader can see which base the variant changes."""

    def setUp(self):
        gr._HAS_VARIANTS = None
        gr._DROP_KINDS = None
        self.addCleanup(setattr, gr, "_HAS_VARIANTS", None)
        self.addCleanup(setattr, gr, "_DROP_KINDS", None)

    def test_the_three_alignment_columns_are_adjacent_and_in_reading_order(self):
        headers = gr.curated_headers(has_crispr_bulge=True)
        i = headers.index("Aligned_spacer+PAM")
        self.assertEqual(headers[i:i + 3],
                         ["Aligned_spacer+PAM", "Aligned_protospacer+PAM_REF",
                          "Aligned_protospacer+PAM"])

    def test_both_alignments_of_a_variant_row(self):
        frame = _curated([{}])
        self.assertEqual(frame["Aligned_protospacer+PAM_REF"][0], _REF_ALN)
        self.assertEqual(frame["Aligned_protospacer+PAM"][0], _ALT_ALN)
        self.assertEqual(frame["Aligned_spacer+PAM"][0], _SPACER)

    def test_a_row_with_no_reference_alignment(self):
        """A variant-created site absent from the reference (Not_found_in_REF): there is
        no reference protospacer, so the new column is "-" -- and the ALT column keeps
        today's value, the sequence that does exist."""
        frame = _curated([{"Aligned_protospacer+PAM_REF_(highest_CFD)": "NA"}])
        self.assertEqual(frame["Aligned_protospacer+PAM_REF"][0], gr.CURATED_MISSING)
        self.assertEqual(frame["Aligned_protospacer+PAM"][0], _ALT_ALN)

    def test_the_alt_column_still_falls_back_to_the_reference_alignment(self):
        """A reference row of a variant run has no ALT alignment: unchanged behaviour,
        the ALT column shows the REF one -- and now the REF column shows it too, which
        is the truth about such a row (one alignment, no variant)."""
        frame = _curated([{"Aligned_protospacer+PAM_ALT_(highest_CFD)": "NA"}])
        self.assertEqual(frame["Aligned_protospacer+PAM"][0], _REF_ALN)
        self.assertEqual(frame["Aligned_protospacer+PAM_REF"][0], _REF_ALN)

    def test_a_reference_only_run_grows_no_duplicate_column(self):
        """No variant off-target anywhere: every row IS the reference, so a REF column
        would repeat "Aligned_protospacer+PAM" on every row. Dropped, as Observed is."""
        gr._HAS_VARIANTS = False
        headers = gr.curated_headers(has_crispr_bulge=True)
        self.assertNotIn("Aligned_protospacer+PAM_REF", headers)
        self.assertIn("Aligned_protospacer+PAM", headers)
        frame = _curated([{"Aligned_protospacer+PAM_ALT_(highest_CFD)": "NA"}])
        self.assertNotIn("Aligned_protospacer+PAM_REF", frame.columns)
        self.assertEqual(frame["Aligned_protospacer+PAM"][0], _REF_ALN)

    def test_an_assembly_run_grows_no_duplicate_column(self):
        """assembly-search never runs with a --vcf (each haplotype IS the genome), so it
        drops the REF/ALT columns -- the new one with them."""
        gr._DROP_KINDS = gr._NO_VARIANT_KINDS
        self.assertNotIn("Aligned_protospacer+PAM_REF",
                         gr.curated_headers(has_crispr_bulge=True))
        self.assertIn("Aligned_protospacer+PAM", gr.curated_headers(has_crispr_bulge=True))


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestLogAxis(unittest.TestCase):
    """cb_operating_point.log_axis: where the floor goes, and that the operating point
    is always on the axis."""

    def test_the_floor_is_one_decade_below_the_smallest_non_zero_score(self):
        axis = cbop.log_axis([0.0, 2.35e-39, 1e-20, 0.4])
        self.assertEqual(axis.floor, 1e-40)
        self.assertLess(axis.bottom, axis.floor)
        self.assertGreater(axis.band_top, axis.floor)
        self.assertEqual(axis.top, 1.0)
        # the smallest real score stays clear of the "= 0" band
        self.assertGreater(2.35e-39, axis.band_top)

    def test_the_threshold_is_always_inside_the_axis(self):
        """Every score above the operating point: the axis must still reach down to it,
        or the panel's line -- the whole point of the panel -- would be off-screen."""
        axis = cbop.log_axis([0.4, 0.9, 1.0])
        self.assertLess(axis.floor, cbop.sensitive_threshold())

    def test_scores_that_are_not_scores_set_no_bound(self):
        """-1 is the not-scored sentinel, not a score of 0; NaN and text are not scores."""
        plain = cbop.log_axis([])
        for values in ([-1.0, -1.0], [float("nan")], ["NA", ""], [2.0]):
            self.assertEqual(cbop.log_axis(values), plain, values)
        self.assertLess(plain.floor, cbop.sensitive_threshold())

    def test_an_all_zero_panel_still_has_an_axis(self):
        axis = cbop.log_axis([0.0, 0.0])
        self.assertGreater(axis.floor, 0.0)
        self.assertLess(axis.floor, cbop.sensitive_threshold())


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestScorePanels(unittest.TestCase):
    """generate_report's Section-2 score panels."""

    @staticmethod
    def _realistic(n=400):
        """A run as they come: most sites far below 1e-20, a few near 1, some exactly 0."""
        rng = np.random.default_rng(11)
        cb = 10.0 ** rng.uniform(-39.0, -18.0, n)
        cb[:4] = [0.97, 0.31, cbop.sensitive_threshold() * 2, 1e-5]
        cb[-int(0.04 * n):] = 0.0          # ~4% exact zeros, as in a real run
        return pd.DataFrame({
            "cb": cb, "cb_ref": cb, "cb_alt": np.clip(cb * 10.0, 0.0, 1.0),
            "maf": ["0.01"] * n, "samples": ["HG1,HG2"] * n, "rsid": ["rs1"] * n,
        })

    def _panel(self, sub, **kwargs):
        return _axes_of(
            sub=sub, score_key="cb", ref_key="cb_ref", alt_key="cb_alt",
            xlabel="x", title="t", score_name="CRISPR-Bulge",
            maf_col="maf", samp_col="samples", rsid_col="rsid", **kwargs,
        )

    def test_a_realistic_crispr_bulge_panel_is_log_and_shows_the_threshold(self):
        sub = self._realistic()
        ax, fig = self._panel(sub, log_y=True)
        self.addCleanup(plt.close, fig)
        self.assertEqual(ax.get_yscale(), "log")
        # the sensitive operating point is drawn, at ITS value (never hard-coded here)
        self.assertIn(cbop.sensitive_threshold(), _hlines(ax))
        # ...and named in the legend, with what the bottom band means
        labels = [t.get_text() for leg in ax.findobj(match=plt.matplotlib.legend.Legend)
                  for t in leg.get_texts()]
        self.assertIn(cbop.sensitive_plot_label(), labels)
        self.assertTrue([l for l in labels if l.startswith("score 0")], labels)
        # the axis reaches below the smallest non-zero score, and above it to 1
        floor = cbop.log_axis(sub["cb"].tolist() + sub["cb_alt"].tolist()).floor
        bottom, top = ax.get_ylim()
        self.assertLess(bottom, floor)
        self.assertGreaterEqual(top, 1.0)

    def test_no_point_is_dropped_and_the_zeros_sit_on_the_floor(self):
        """A log axis cannot show 0: the zeros are clipped onto the floor, inside the
        marked band, not silently left out of the panel."""
        sub = self._realistic()
        ax, fig = self._panel(sub, log_y=True)
        self.addCleanup(plt.close, fig)
        axis = cbop.log_axis(sub["cb"].tolist() + sub["cb_alt"].tolist())
        drawn = ax.collections[0].get_offsets()
        self.assertEqual(len(drawn), len(sub))                  # every site is there
        n_zero = int((sub["cb"] == 0).sum())
        self.assertGreater(n_zero, 0)
        on_floor = [p for p in drawn.tolist() if p[1] == axis.floor]
        self.assertEqual(len(on_floor), n_zero)
        self.assertTrue(all(p[1] >= axis.floor for p in drawn.tolist()))

    def test_the_not_scored_sentinel_is_not_drawn_as_a_zero(self):
        """-1 means "the model did not score this row". On the linear axis it fell off
        the 0-1 window; on the log one it must not land on the zero floor."""
        sub = self._realistic(50)
        sub.loc[sub.index[:10], ["cb", "cb_ref", "cb_alt"]] = -1.0
        ax, fig = self._panel(sub, log_y=True)
        self.addCleanup(plt.close, fig)
        floor = cbop.log_axis(sub["cb"].tolist()).floor
        # a masked offset (NaN y) is a point matplotlib draws nowhere: None in tolist()
        ys = [p[1] for p in ax.collections[0].get_offsets().tolist()]
        self.assertEqual(sum(1 for y in ys if y is None), 10)   # the 10 sentinel rows
        drawn = [y for y in ys if y is not None]
        self.assertEqual(sum(1 for y in drawn if y == floor), int((sub["cb"] == 0).sum()))
        self.assertEqual(len(drawn), len(sub) - 10)

    def test_the_ref_to_alt_arrows_survive_the_log_axis(self):
        """A pair that goes from 1e-30 to 1e-20 has moved ten decades -- most of the
        panel -- though ALT minus REF is 1e-20. The arrow is drawn by the RATIO of the
        two, or every arrow of a genome-wide run would be dropped as "no movement"."""
        ref = np.array([10.0 ** (-2.0 - i) for i in range(25)])
        sub = pd.DataFrame({"cb": ref, "cb_ref": ref,
                            "cb_alt": np.clip(ref * 1e6, 0.0, 1.0)})
        ax, fig = self._panel(sub, log_y=True)
        self.addCleanup(plt.close, fig)
        arrows = [t for t in ax.texts if getattr(t, "arrow_patch", None) is not None]
        self.assertEqual(len(arrows), len(sub))
        for arrow, yr, ya in zip(arrows, sub["cb_ref"], sub["cb_alt"]):
            self.assertAlmostEqual(arrow.xyann[1], yr)   # tail on the REF point
            self.assertAlmostEqual(arrow.xy[1], ya)      # head on the ALT point
        # a site whose two scores are the same still gets no arrow
        same = pd.DataFrame({"cb": ref, "cb_ref": ref, "cb_alt": ref})
        ax2, fig2 = self._panel(same, log_y=True)
        self.addCleanup(plt.close, fig2)
        self.assertEqual([t for t in ax2.texts
                          if getattr(t, "arrow_patch", None) is not None], [])

    def test_the_allele_frequency_sizing_is_the_same_on_both_axes(self):
        """The AF encoding (marker area) is a property of the row, not of the y-axis."""
        sub = self._realistic(60)
        sub["maf"] = [f"{(i % 10) / 10 + 0.001:.3f}" for i in range(len(sub))]
        sizes = {}
        for log_y in (True, False):
            ax, fig = self._panel(sub, log_y=log_y)
            self.addCleanup(plt.close, fig)
            sizes[log_y] = [list(c.get_sizes()) for c in ax.collections]
        self.assertEqual(sizes[True], sizes[False])
        self.assertGreater(len(set(sizes[True][1])), 1)  # AF really is encoded

    def test_the_cfd_panel_is_unchanged(self):
        sub = self._realistic()
        ax, fig = self._panel(sub, log_y=False)
        self.addCleanup(plt.close, fig)
        self.assertEqual(ax.get_yscale(), "linear")
        self.assertEqual(ax.get_ylim(), (0.0, 1.0))
        self.assertNotIn(cbop.sensitive_threshold(), _hlines(ax))
        labels = [t.get_text() for leg in ax.findobj(match=plt.matplotlib.legend.Legend)
                  for t in leg.get_texts()]
        self.assertNotIn(cbop.sensitive_plot_label(), labels)

    def test_a_panel_with_no_valid_score_still_raises_for_the_placeholder(self):
        """The caller renders its labelled placeholder on ValueError: a run whose scores
        are all the -1 sentinel must not come back as an empty log plot."""
        sub = self._realistic(20)
        sub[["cb", "cb_ref", "cb_alt"]] = -1.0
        for log_y in (True, False):
            with self.subTest(log_y=log_y), self.assertRaises(ValueError):
                self._panel(sub, log_y=log_y)

    def test_the_panels_say_in_their_caption_what_the_axis_does(self):
        caption = gr._cb_log_caption("CRISPR-Bulge")
        self.assertIn("log10", caption)
        self.assertIn(cbop.sensitive_label_html(), caption)
        for word in ("0", "shaded band", "smallest non-zero"):
            self.assertIn(word, caption)


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestTopPlotScripts(unittest.TestCase):
    """The results page's top-1000 plots and the Personal Risk Card's, which draw the
    same two scores: the CRISPR-Bulge one gets the same log axis as the report's."""

    CRITERIA = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")

    @classmethod
    def _sites(cls, n=60):
        rng = np.random.default_rng(5)
        cb = 10.0 ** rng.uniform(-38.0, -9.0, n)
        cb[0], cb[1] = 0.84, cbop.sensitive_threshold() * 3
        cb[-3:] = 0.0
        rows = []
        for i in range(n):
            row = {"Spacer+PAM": _SPACER, "Chromosome": "chr1",
                   "Start_coordinate": 1000 * (i + 1)}
            for crit in cls.CRITERIA:
                row[f"Mismatches+bulges_({crit})"] = 0 if i == 0 else 3
                row[f"Bulges_({crit})"] = 0
                row[f"REF/ALT_origin_({crit})"] = "alt" if i % 2 else "ref"
                row[f"Variant_MAF_({crit})"] = 0.2 if i % 2 else None
                row[f"Variant_samples_({crit})"] = "S1,S2" if i % 2 else None
                row[f"Aligned_protospacer+PAM_REF_({crit})"] = _REF_ALN
            row["CFD_score_(highest_CFD)"] = 0.3
            row["CFD_score_REF_(highest_CFD)"] = 0.3
            row["CFD_score_ALT_(highest_CFD)"] = 0.3
            row["CRISPR_BULGE_score_(highest_CRISPR_BULGE)"] = cb[i]
            row["CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)"] = cb[i]
            row["CRISPR_BULGE_score_ALT_(highest_CRISPR_BULGE)"] = (
                min(1.0, cb[i] * 1e3) if i % 2 else cb[i])
            rows.append(row)
        return pd.DataFrame(rows)

    @staticmethod
    def _drawn_axis(plot, df, *args):
        """Run a plot script's function and return its Axes, captured at savefig."""
        seen = {}

        def _capture(*_a, **_k):
            seen.setdefault("ax", plt.gca())
            seen.setdefault("points", [c.get_offsets().tolist()
                                       for c in plt.gca().collections])
            seen.setdefault("ylim", plt.gca().get_ylim())
            seen.setdefault("yscale", plt.gca().get_yscale())
            seen.setdefault("hlines", _hlines(plt.gca()))
            legend = plt.gca().get_legend()
            seen.setdefault("legend",
                            [t.get_text() for t in legend.get_texts()] if legend else [])

        import tempfile

        with tempfile.TemporaryDirectory() as out, \
                mock.patch.object(plt, "savefig", _capture):
            plot(df, out + os.sep, *args)
        return seen

    def test_both_scripts_draw_crispr_bulge_on_the_log_axis_with_the_threshold(self):
        import CRISPRme_plots as plots
        import CRISPRme_plots_personal as personal

        df = self._sites()
        axis = cbop.log_axis(df["CRISPR_BULGE_score_(highest_CRISPR_BULGE)"].tolist())
        for name, plot, table in (
            ("CRISPRme_plots", plots.plot_with_CRISPR_BULGE_score,
             plots.filter_table(df, "crispr_bulge").copy()),
            ("CRISPRme_plots_personal", personal.plot_with_CRISPR_BULGE_score, df.copy()),
        ):
            with self.subTest(script=name):
                seen = self._drawn_axis(plot, table, "GUIDE")
                self.assertEqual(seen["yscale"], "log")
                self.assertIn(cbop.sensitive_threshold(), seen["hlines"])
                self.assertIn(cbop.sensitive_plot_label(), seen["legend"])
                self.assertTrue([l for l in seen["legend"] if l.startswith("score 0")],
                                seen["legend"])
                self.assertLess(seen["ylim"][0], axis.floor)   # room for the "= 0" band
                # the zeros are drawn, on the floor
                ref_points = seen["points"][0]
                self.assertEqual(len(ref_points), len(table))
                self.assertEqual(sum(1 for p in ref_points if p[1] == axis.floor), 3)

    def test_the_cfd_plots_of_both_scripts_are_unchanged(self):
        import CRISPRme_plots as plots
        import CRISPRme_plots_personal as personal

        df = self._sites()
        for name, plot, table in (
            ("CRISPRme_plots", plots.plot_with_CFD_score,
             plots.filter_table(df, "cfd").copy()),
            ("CRISPRme_plots_personal", personal.plot_with_CFD_score, df.copy()),
        ):
            with self.subTest(script=name):
                seen = self._drawn_axis(plot, table, "GUIDE")
                self.assertEqual(seen["yscale"], "linear")
                self.assertEqual(seen["ylim"], (0.0, 1.0))
                self.assertNotIn(cbop.sensitive_threshold(), seen["hlines"])
                self.assertNotIn(cbop.sensitive_plot_label(), seen["legend"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
