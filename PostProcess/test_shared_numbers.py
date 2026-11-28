#!/usr/bin/env python3
"""One quantity, one source: every number stated on more than one surface.

The off-target work is spread over surfaces that are each reviewed on their own --
``docs/offtarget_risk_summary_METHODS.md``, ``METHODS.md``, ``README.md``, the
risk-card modal, two PNG figures, ``PostProcess/cb_operating_point.py``,
``PostProcess/offtarget_calibration.py``, ``scripts/fit_calibration.py``,
``scripts/render_pr_figures.py`` and the changelog. Reviewing them one at a time is
exactly how the same quantity ends up stated three different ways. It has happened
repeatedly:

  * the 95%->99% validation-list ratio read 6x in the modal, ~6x in
    ``cb_operating_point``'s docstring and 7.1x in the methods document, where the
    benchmark's own figure is 16.4x;
  * the CRISPR-Bulge bulge=1 average precision read 0.25 in one place and 0.285 on
    the figure that measures it;
  * the held-out CRISPR-Bulge AUPR was printed twice in one modal without saying why
    the two differed;
  * ``287`` and ``359`` -- counts of validated POSITIVES -- were printed as counts of
    SITES in three places;
  * a whole methods block quoted a calibration fit that had been superseded;
  * the sites-listed-per-validated-site column read 9 / 21 / 354 / 790 where the
    committed curve data gives 9.05 / 22.15 / 348.79 / 786.08.

Each of those was a single surface read in isolation. This file reads them together.
Every number below is recomputed from the one artifact that owns it --

    PostProcess/calibration/pr_figure_data.npz      the curves, their average
                                                    precisions, their supports and
                                                    the report's marked operating
                                                    points
    PostProcess/calibration/*_calibration.npz       AUPR + CI, phi per band,
                                                    band_counts, per-held-out-guide
                                                    coverage, prevalence
    PostProcess/cb_operating_point.py               the three recall targets, their
                                                    thresholds, CIs and label text

-- and then matched against every copy of it, prose included, so a surface that stops
stating the source's value fails here rather than in a collaborator's hands.

Be exact about what that buys, because the first version of this paragraph was not. A
presence check sees a value DISAPPEAR from a surface -- which is what a drift does when
the surface states the quantity once, and it is what every case of the mutation harness
does. It does not see one copy of a twice-stated quantity drift while the other copy
keeps it present. The two helpers that count are what cover that: assertNumberOnEvery
accepts exactly one rounding per surface, and assertOnEvery(..., times=n) pins how many
times. Where a surface legitimately repeats a number, the second copy is held by those
and not by the plain assertion.

Nothing here re-measures the benchmark: the benchmark is not in the repository. What
is in the repository is the curve data the figures are drawn from, and every list
size, ratio and average precision this file checks is derived from it at run time.

Run as CI runs it::

    cd PostProcess && python3 test_shared_numbers.py
"""
import json
import math
import os
import re
import sys
import unittest

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_CALDIR = os.path.join(_HERE, "calibration")
_SCRIPTS = os.path.join(_REPO, "scripts")
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

import cb_operating_point as cbop            # noqa: E402
import render_pr_figures as rpf              # noqa: E402  (compact_sci: one formatter)

_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

# Every surface in this file describes the SHIPPED operating point, and the committed
# curve data marks only that one. Under CRISPRME_CB_RECALL_TARGET=90 or 99 the module
# returns a threshold no document quotes and no mark carries -- which is not drift, it is
# the variable doing its job. The checks that key on it say so and step aside; everything
# that does not (the AUPRs, the supports, the bulge APs, the calibration claims) still runs.
_only_default = unittest.skipUnless(
    cbop.recall_target() == cbop.DEFAULT_TARGET,
    f"CRISPRME_CB_RECALL_TARGET={cbop.recall_target()} selects an operating point the "
    f"documents and the committed curve data do not describe (they carry the shipped "
    f"{cbop.DEFAULT_TARGET}%-recall default)")


def _read(*parts):
    with open(os.path.join(_REPO, *parts), encoding="utf-8") as fh:
        return fh.read()


def _curve_data():
    return np.load(os.path.join(_CALDIR, "pr_figure_data.npz"), allow_pickle=False)


def _artifact(scorer):
    path = os.path.join(_CALDIR, f"{scorer}_calibration.npz")
    return json.loads(str(np.load(path, allow_pickle=False)["meta"]))


def _modal():
    import risk_card
    _card, modal = risk_card.render_risk_card(
        {"guide_display": "G"}, [0.6, 0.3, 0.05], [0.6, 0.3, 0.05])
    return modal


_PAYLOAD = re.compile(r"base64,[A-Za-z0-9+/=]+")


def _prose(html):
    """The modal with its embedded image payloads taken out.

    The rendered modal is 261k characters and 249k of them are the base64 of two PNGs
    -- random-looking text in which "16" occurs 146 times, "46" 75 and "40" 41. Every
    assertion below about a SHORT token was therefore satisfiable by image noise rather
    than by the sentence it is about: `assertIn("16", modal)` could not fail. The
    figures are checked as figures (by sha256, in TestTheShippedExamplePage); what the
    prose says is checked here, against the prose.
    """
    return _PAYLOAD.sub("base64,<payload>", html)


def _newest_changelog_section():
    """The section the entries this suite holds are written in.

    While a release is being written the entries are in ``[Unreleased]``. At the release
    they move into the dated ``[X.Y.Z]`` section below it (``docs/RELEASING.md`` step 3).
    So BOTH are returned, always.

    This used to return whichever of the two was non-empty, newest first. That guessed
    that anything in ``[Unreleased]`` must be these entries -- and the first unrelated
    entry added there (a companion file, a new column) hid the released section and made
    every assertion below fail at once, reporting a dozen numbers as missing from a
    changelog that states all of them. Reading both keeps the protection the single-
    section read was for -- the numbers must be stated in a CURRENT entry, not in some
    ancient one -- without assuming what else the open section contains.
    """
    txt = _read("CHANGELOG.md")
    start = txt.index("## [Unreleased]")
    end = txt.index("\n## [", start + 1)
    nxt = end + 1
    return txt[start:end] + "\n" + txt[nxt:txt.index("\n## [", nxt + 1)]


def _block(text, opening, closing):
    """One section of a long document.

    METHODS.md and README.md are thousands of lines; a bare substring test on "0.56"
    or "46" would pass against some unrelated table and so would not notice the risk
    card's own numbers disappearing. Every assertion about those two files is made
    against the risk-card section alone, which is what this cuts out.
    """
    start = text.index(opening)
    return text[start:text.index(closing, start + len(opening))]


def _times_ten(x, digits=1):
    """0.00019 -> '1.9×10⁻⁴', the way the markdown documents write a threshold."""
    mant, exp = f"{x:.{digits}e}".split("e")
    return f"{mant}×10{str(int(exp)).translate(_SUP)}"


def _list_at_recall(d, panel, k):
    """(min list, max list) a threshold that recovers exactly ``k`` positives produces.

    Recall is k/n_pos; the curve may hold a whole vertical run at that recall (constant
    recall, falling precision), whose interior points the stored polyline drops as
    collinear. The tightest list -- the one the k-th best positive's own score gives --
    is the run's highest precision; the loosest is its lowest.
    """
    n_pos = int(d[f"{panel}_pos"])
    rec, prec = d[f"{panel}_rec"], d[f"{panel}_prec"]
    on = np.isclose(rec, k / n_pos, rtol=0, atol=1e-12)
    if not on.any():
        raise AssertionError(f"recall {k}/{n_pos} is not a vertex of the {panel} curve")
    return k / float(prec[on].max()), k / float(prec[on].min())


def _tightest_list(d, panel, k):
    return int(round(_list_at_recall(d, panel, k)[0]))


class _Surfaces(unittest.TestCase):
    """The documents every subclass reads, loaded once."""

    @classmethod
    def setUpClass(cls):
        cls.d = _curve_data()
        cls.cfd = _artifact("cfd")
        cls.cb = _artifact("crispr_bulge")
        cls.doc = _read("docs", "offtarget_risk_summary_METHODS.md")
        cls.methods = _read("METHODS.md")
        cls.readme = _read("README.md")
        cls.methods_card = _block(
            cls.methods, "**Top-of-report off-target risk summary", "\n\n**")
        cls.readme_card = _block(
            cls.readme, "- **Off-target risk summary (top of the report).**", "\n- **")
        cls.changelog = _newest_changelog_section()
        cls.cbop_doc = cbop.__doc__
        cls.calib_src = _read("PostProcess", "offtarget_calibration.py")
        cls.fit_src = _read("scripts", "fit_calibration.py")
        cls.modal_render = _modal()
        # the SURFACE is what the modal says; the image bytes are checked as images
        cls.modal = _prose(cls.modal_render)
        assert len(cls.modal) < len(cls.modal_render), (
            "the modal no longer embeds a base64 payload: _prose() is now a no-op and "
            "the short-token assertions below have lost the guard it gives them")

    def surfaces(self, *names):
        return [(n, getattr(self, n)) for n in names]

    def assertOnEvery(self, value, *names, times=None):
        """``value`` is on each named surface -- and, with ``times``, how often.

        Counting rather than ``assertIn`` for the same reason assertNumberOnEvery
        counts: a quantity written twice on one surface can drift in one copy and keep
        the test green on the other, so where a surface should state something exactly
        once the test says so instead of accepting "at least once".
        """
        for name, text in self.surfaces(*names):
            found = text.count(value)
            with self.subTest(surface=name, value=value):
                self.assertTrue(found, f"{name} does not state {value!r}")
                if times is not None:
                    self.assertEqual(found, times,
                                     f"{name} states {value!r} {found} times, "
                                     f"expected {times}")

    def assertPatternOnEvery(self, pattern, *names):
        """A quantity in the shape of the sentence it belongs to, not as a bare token."""
        for name, text in self.surfaces(*names):
            with self.subTest(surface=name, pattern=pattern):
                self.assertRegex(text, pattern, f"{name} does not state {pattern!r}")

    def assertOnNone(self, value, *names):
        for name, text in self.surfaces(*names):
            with self.subTest(surface=name, value=value):
                self.assertNotIn(value, text, f"{name} still states {value!r}")

    def assertNumberOnEvery(self, value, *names, places=(1, 2)):
        """The same NUMBER on each surface, whatever precision that surface writes it to.

        A docstring table that prints 22.3 and a document that prints 22.26 agree; a
        document that prints 21 does not. Holding every surface to one rounding would
        make the test enforce typography, and it would be the thing that drifts.
        """
        wanted = [f"{value:.{p}f}" for p in places]
        for name, text in self.surfaces(*names):
            present = [w for w in wanted if w in text]
            with self.subTest(surface=name, value=wanted):
                self.assertTrue(present, f"{name} states none of {wanted}")
                # and only one of them: a surface that writes the same number two ways
                # is where the next divergence starts, and it defeats this check -- one
                # spelling can then drift while the other keeps the test green
                self.assertEqual(len(present), 1,
                                 f"{name} writes the same number two ways: {present}")


# --------------------------------------------------------------------------- #
# 1. The operating points and their labels
# --------------------------------------------------------------------------- #
@_only_default
class TestOperatingPoints(_Surfaces):
    """cb_operating_point.CB_RECALL_POINTS is the one source; nothing may restate it
    differently, and the curve data's marked point must be the same number."""

    def test_the_label_every_reader_sees_is_built_from_the_threshold(self):
        thr, target = cbop.sensitive_threshold(), cbop.recall_target()
        self.assertEqual(cbop.sensitive_label_text(), f"≥ {thr:.2g} ({target}% recall)")
        self.assertEqual(cbop.sensitive_label_html(), f"&ge; {thr:.2g} ({target}% recall)")
        self.assertEqual(cbop.sensitive_plot_label(),
                         f"sensitive threshold {thr:.2g} ({target}% recall)")
        # and the label, not a typed copy of it, is what the card and the modal print
        self.assertIn(cbop.sensitive_label_text(), self.modal)
        self.assertIn(cbop.sensitive_label_html(), self.modal)

    def test_the_marked_point_on_the_figure_is_the_modules_threshold(self):
        meta = json.loads(str(self.d["meta"]))
        self.assertEqual(meta["sensitive_threshold"], cbop.sensitive_threshold())
        self.assertEqual(meta["recall_target"], cbop.recall_target())
        marks = json.loads(str(self.d["marks"]))
        sens = [m for m in marks
                if m[0] == "cb" and m[1] == cbop.sensitive_threshold()]
        self.assertEqual(len(sens), 1, "the sensitive point is not marked exactly once")

    def test_every_document_quotes_the_threshold_the_module_returns(self):
        thr = cbop.sensitive_threshold()
        # the three spellings the surfaces use, all derived from the one value
        self.assertOnEvery(rpf.compact_sci(thr), "cbop_doc", "readme", "changelog")
        self.assertOnEvery(_times_ten(thr), "doc", "methods")
        self.assertOnEvery(f"{thr:.2g}", "readme", "modal")

    def test_the_95_percent_interval_is_the_modules_on_every_surface(self):
        _thr, lo, hi = cbop.CB_RECALL_POINTS[95]
        self.assertOnEvery(f"[{rpf.compact_sci(lo)}, {rpf.compact_sci(hi)}]", "cbop_doc")
        self.assertOnEvery(f"[{_times_ten(lo)}, {_times_ten(hi)}]", "doc", "methods")

    def test_the_other_two_recall_points_are_stated_once_and_completely(self):
        """90% and 99% live only in cb_operating_point. Its docstring table must print
        each threshold with its own CI, so the table cannot drift from the dict."""
        for target in (90, 95, 99):
            thr, lo, hi = cbop.CB_RECALL_POINTS[target]
            row = (rf"{target}%\s+{re.escape(rpf.compact_sci(thr))}\s+"
                   rf"\[{re.escape(rpf.compact_sci(lo))}, "
                   rf"{re.escape(rpf.compact_sci(hi))}\]")
            with self.subTest(target=target):
                self.assertRegex(self.cbop_doc, row)
        self.assertEqual(cbop.DEFAULT_TARGET, 95)
        self.assertEqual(sorted(cbop.CB_RECALL_POINTS), [90, 95, 99])


# --------------------------------------------------------------------------- #
# 2. The two AUPRs and their intervals
# --------------------------------------------------------------------------- #
class TestAuprsAndIntervals(_Surfaces):

    def test_each_aupr_and_ci_is_the_artifacts_wherever_it_is_quoted(self):
        for meta in (self.cfd, self.cb):
            point = meta["aupr"]["point"]
            lo, hi = meta["aupr"]["ci"]
            self.assertOnEvery(f"{point:.3f}", "doc", "methods_card", "readme_card")
            for bound in (lo, hi):
                self.assertOnEvery(f"{bound:.2f}", "doc", "methods_card", "readme_card")
        # the modal renders both rows straight out of the artifacts
        for meta in (self.cfd, self.cb):
            self.assertIn(f'{meta["aupr"]["point"]:.3f}', self.modal)
            lo, hi = meta["aupr"]["ci"]
            self.assertIn(f"[{lo:.2f}–{hi:.2f}]", self.modal)

    def test_the_figures_aupr_is_the_curve_datas_wherever_it_is_quoted(self):
        ap = float(self.d["cb_ap"])
        self.assertEqual(rpf.overall_titles(self.d)[1].count(f"{ap:.2f}"), 1)
        self.assertOnEvery(f"{ap:.2f}", "doc", "modal", "changelog")
        self.assertOnEvery(f"{ap:.3f}", "doc")
        self.assertOnEvery(f"{ap:.6f}", "changelog")

    def test_the_two_crispr_bulge_auprs_are_a_gap_both_pages_explain_the_same_way(self):
        """The artifact's 0.613 and the figure's 0.622 are the same 287 positives. The
        gap, its three terms and the CI remark must read the same in the document and
        in the modal; §2 is the one place the defect itself is measured."""
        gap = float(self.d["cb_ap"]) - self.cb["aupr"]["point"]
        self.assertAlmostEqual(gap, 0.0089928, places=6)
        lo, hi = self.cb["aupr"]["ci"]
        self.assertLess(lo, float(self.d["cb_ap"]))
        self.assertLess(float(self.d["cb_ap"]), hi)
        import risk_card
        for name, text in (("METHODS", self.doc), ("modal", risk_card._CB_AUPR_FIGURE_NOTE)):
            for term in ("54%", "26%", "20%"):
                with self.subTest(surface=name, term=term):
                    self.assertIn(term, text)
        self.assertIn("0.6131 → 0.6180", self.doc)

    def test_no_live_surface_calls_the_cfd_column_held_out(self):
        """CFD has no training set, so nothing is held out for it: its support is all
        46 benchmark guides. "held-out AUPR" over a CFD number was corrected on the
        figure caption, the modal and §3 and was still on METHODS.md and README.md."""
        for name in ("doc", "methods", "readme", "modal"):
            with self.subTest(surface=name):
                self.assertIsNone(re.search(r"[Hh]eld-out AUPR", getattr(self, name)),
                                  f"{name} calls an AUPR held-out without naming the scorer")
        self.assertIn("CFD: all 46, nothing held out", rpf.OVERALL_SUPTITLE)
        for name, text in self.surfaces("doc", "methods_card", "readme_card", "modal"):
            with self.subTest(surface=name):
                self.assertTrue(
                    "no training set" in text or "nothing is held out" in text,
                    f"{name} does not say why CFD's column is not held out")

    def test_the_on_target_tie_group_closes_arithmetically_on_both_pages(self):
        """Why CFD moves most when the on-targets are dropped: its 46 on-targets tie at
        exactly 1.0 at the head of the ranking, and that one tie group is worth most of
        the headline. The document and the modal both spell the multiplication out, so
        the three numbers have to multiply and the share has to be a share of the
        artifact's own AUPR -- otherwise a reader can check the arithmetic and find it
        does not close."""
        import risk_card
        precision, recall, share = 0.902, 0.128, 0.116
        self.assertAlmostEqual(precision * recall, share, places=2)
        self.assertEqual(f"{share / self.cfd['aupr']['point'] * 100:.0f}%", "69%")
        card_src = _read("PostProcess", "risk_card.py")
        for name, text in (("METHODS", self.doc), ("modal", risk_card._ONTARGET_NOTE)):
            with self.subTest(surface=name):
                self.assertIn("69%", text)
        for value in ("0.902", "0.128", "0.116"):
            with self.subTest(value=value):
                self.assertIn(value, self.doc)
                self.assertIn(value, card_src)

    def test_the_deployment_matched_values_agree_across_the_three_places(self):
        """0.025 / 0.521 are measured, not derivable here, so the three surfaces that
        print them are held to each other -- and to being BELOW the published ones."""
        import risk_card
        for value in ("0.025", "0.521"):
            self.assertOnEvery(value, "doc", "modal")
            self.assertIn(value, rpf.OVERALL_SUPTITLE)
        self.assertIn("0.025", risk_card._ONTARGET_NOTE)
        self.assertLess(0.025, self.cfd["aupr"]["point"])
        self.assertLess(0.521, float(self.d["cb_ap"]))


# --------------------------------------------------------------------------- #
# 3. The support counts
# --------------------------------------------------------------------------- #
class TestSupportCounts(_Surfaces):

    def test_the_artifacts_and_the_curve_data_describe_the_same_two_supports(self):
        self.assertEqual((self.cfd["n_positives"], self.cfd["n_sites"], self.cfd["n_guides"]),
                         (int(self.d["cfd_pos"]), int(self.d["cfd_n"]),
                          int(self.d["scfd_guides"])))
        self.assertEqual((self.cb["n_positives"], self.cb["n_sites"], self.cb["n_guides"]),
                         (int(self.d["cb_pos"]), int(self.d["cb_n"]),
                          int(self.d["scb_guides"])))

    def test_every_surface_quotes_those_counts_and_not_others(self):
        """Each count in the phrase that makes it that quantity, not as a bare token.

        "46" and "40" are two characters. Asserted bare they were satisfied by any two
        digits anywhere on the surface -- which, before _prose(), included 75 and 41
        hits inside the modal's base64 images. The phrase is what a reader reads and it
        is what has to survive: a support count that loses its noun has still drifted.
        """
        for meta in (self.cfd, self.cb):
            pos, sites = meta["n_positives"], meta["n_sites"]
            with self.subTest(scorer=meta["scorer"]):
                self.assertOnEvery(f"{pos} positives", "doc", "changelog")
                # the modal names the two together ("359 validated off-targets among
                # 2,940,588 candidate sites", "287 among 2,703,869 sites")
                self.assertPatternOnEvery(rf"\b{pos}\b[^<>]*among {sites:,}", "modal")
                self.assertOnEvery(f"{sites:,}", "doc", "modal", "changelog")
                self.assertOnEvery(f'{meta["n_guides"]} guides', "doc", "modal")
        # the 359/46 pair also carries METHODS.md's and README.md's risk-card sections
        self.assertOnEvery(f'{self.cfd["n_positives"]} experimentally confirmed '
                           f'off-targets', "methods_card", "readme_card")
        self.assertOnEvery(f'{self.cfd["n_guides"]} guides', "methods_card", "readme_card")
        # and README's section names CRISPR-Bulge's own support beside it
        self.assertOnEvery(f'{self.cb["n_positives"]} positives', "readme_card")
        self.assertOnEvery(f'{self.cb["n_guides"]} guides', "readme_card")

    def test_a_positive_is_never_printed_as_a_site(self):
        """"CFD 359 sites / 46 guides" named the positive count as the site count, four
        orders of magnitude out. 2,940,588 sites carry those 359."""
        for pos in (self.cfd["n_positives"], self.cb["n_positives"]):
            for bad in (f"{pos} sites", f"{pos} held-out sites", f"{pos} validated sites"):
                self.assertOnNone(bad, "doc", "methods", "readme", "modal")

    def test_the_bulge_strata_partition_each_support(self):
        for panel, meta in (("scfd", self.cfd), ("scb", self.cb)):
            with self.subTest(panel=panel):
                b0, b1 = int(self.d[f"{panel}_b0_pos"]), int(self.d[f"{panel}_b1_pos"])
                self.assertEqual(b0 + b1, meta["n_positives"])
                self.assertEqual(int(self.d[f"{panel}_b0_n"]) + int(self.d[f"{panel}_b1_n"]),
                                 meta["n_sites"])
        # §2 prints both splits in words
        self.assertIn(f'{int(self.d["scb_b0_pos"])} + {int(self.d["scb_b1_pos"])}', self.doc)
        self.assertIn(f'{int(self.d["scfd_b0_pos"])} with 0 bulges, '
                      f'{int(self.d["scfd_b1_pos"])} with 1', self.doc)

    def test_the_band_support_the_document_prints_is_the_artifacts(self):
        for meta in (self.cfd, self.cb):
            bc = meta["band_counts"]["0.5-1"]
            with self.subTest(scorer=meta["scorer"]):
                self.assertIn(f'{bc["positives"]} positives from all '
                              f'{bc["guides_with_positive"]} guides', self.doc)
        thin = self.cb["band_counts"]["0.2-0.5"]
        self.assertIn(f'{thin["positives"]} positives from '
                      f'{thin["guides_with_positive"]} guides', self.modal)

    def test_the_dominant_guides_share_is_the_same_fraction_on_both_pages(self):
        """One guide supplies 79 positives. Its share of each support is arithmetic, so
        the two percentages the page prints must be that arithmetic."""
        dominant = 79
        self.assertIn(f"{dominant}/{self.cfd['n_positives']} = "
                      f"{dominant / self.cfd['n_positives'] * 100:.1f}%", self.doc)
        self.assertEqual(f"{dominant / self.cb['n_positives'] * 100:.1f}%", "27.5%")
        self.assertOnEvery("27.5%", "doc", "modal")


# --------------------------------------------------------------------------- #
# 4. The validation list: sizes, ratios and the recall trade
# --------------------------------------------------------------------------- #
class TestValidationListSizes(_Surfaces):
    """Re-derived from the committed curve data, which is the only in-tree record of
    what each threshold lists. These are the numbers that read 6x, 6.2x and 7.1x in
    three files while the benchmark's own figure was 16.4x."""

    def quantile_lists(self):
        n_pos = int(self.d["cb_pos"])
        out = {}
        for target in (0.90, 0.95, 0.99, 1.00):
            k = math.ceil(target * n_pos)
            out[target] = (k, _tightest_list(self.d, "cb", k))
        return out

    def test_the_docstring_table_is_the_curve_datas_own_list_column(self):
        lists = self.quantile_lists()
        for target in (0.90, 0.95, 0.99):
            k, size = lists[target]
            with self.subTest(target=target):
                self.assertIn(f"{size:,}", self.cbop_doc)
                self.assertIn(f"{size / k:.1f}", self.cbop_doc)

    def test_the_methods_document_prints_the_same_column_to_two_decimals(self):
        lists = self.quantile_lists()
        for target in (0.90, 0.95, 0.99, 1.00):
            k, size = lists[target]
            with self.subTest(target=target):
                self.assertIn(f"{size / k:.2f}", self.doc,
                              f"§3.5 does not give {target:.0%} as {size / k:.2f}")

    @_only_default
    def test_the_shipped_threshold_lists_what_the_marked_point_says(self):
        marks = json.loads(str(self.d["marks"]))
        _p, thr, _lab, _prec, _rec, tp, fp = [
            m for m in marks if m[0] == "cb" and m[1] == cbop.sensitive_threshold()][0]
        listed = tp + fp
        self.assertEqual((tp, fp, listed), (273, 5803, 6076))
        self.assertOnEvery(f"{listed:,}", "doc", "modal", "changelog")
        self.assertNumberOnEvery(listed / tp, "doc", "cbop_doc")
        self.assertIn(f"{tp} of {int(self.d['cb_pos'])}", self.doc)

    def test_the_ninety_nine_percent_ratio_is_sixteen_fold_on_every_surface(self):
        lists = self.quantile_lists()
        _k95, l95 = lists[0.95]
        _k99, l99 = lists[0.99]
        ratio = l99 / l95
        self.assertAlmostEqual(ratio, 16.44, places=2)
        self.assertIn(f"{ratio:.1f}×", self.doc)
        self.assertIn(f"{l95:,} → {l99:,}", self.doc)
        self.assertIn(f"{l95:,} -> {l99:,}", self.cbop_doc)
        # the ROUNDED form, in the spelling each surface actually writes. This line
        # used to be assertIn(f"{round(ratio)}", surface) -- the two characters "16",
        # which the modal's own "0.168" contains three times over. It could not fail,
        # so it was not a check on the ratio at all.
        rounded = round(ratio)
        self.assertOnEvery(f"≈{rounded}×", "doc", times=1)
        self.assertOnEvery(f"{rounded}&times;", "modal", times=1)
        self.assertOnEvery(f"~{rounded}x", "cbop_doc", times=1)

    def test_the_shipped_99_percent_list_is_inside_the_curves_own_run(self):
        """99,477 is the list at the ROUNDED 1.5e-8, which sits below the exact recall
        quantile, so it is an interior point of a vertical run the stored polyline
        decimates. It cannot be read off the curve -- but it can be bounded by it, and
        its ratio against the shipped 95% list must still round to the 16x both the
        document and the modal print."""
        n_pos = int(self.d["cb_pos"])
        k = math.ceil(0.99 * n_pos)
        tight, loose = _list_at_recall(self.d, "cb", k)
        shipped_99, shipped_95 = 99477, 6076
        self.assertLessEqual(tight, shipped_99)
        self.assertLessEqual(shipped_99, loose)
        self.assertEqual(round(shipped_99 / shipped_95), 16)
        self.assertOnEvery(f"{shipped_99:,}", "doc", "modal")

    def test_no_live_surface_still_says_the_list_only_sixfolds(self):
        for bad in ("~6x", "roughly 6&times;", "≈6× the 95%"):
            self.assertOnNone(bad, "cbop_doc", "modal", "methods", "readme")
        # the document may say it only while disowning it
        for para in re.split(r"\n\s*\n", self.doc):
            if "≈6×" in para:
                with self.subTest(para=para[:50].replace("\n", " ")):
                    self.assertIn("earlier draft", para.lower())

    def test_the_two_cfd_arms_close_arithmetically_on_both_pages(self):
        base, arm0, arm1, union = 6076, 9660, 614, 16350
        self.assertEqual(base + arm0 + arm1, union)
        self.assertEqual(f"{union / base:.1f}", "2.7")
        self.assertEqual(arm0 + arm1, 10274)
        for value in (f"{arm0:,}", f"{arm1:,}", f"{base:,}"):
            self.assertOnEvery(value, "doc", "modal")
        self.assertOnEvery(f"{union / base:.1f}×", "doc")
        self.assertOnEvery(f"{arm0 + arm1:,}", "doc", "modal")

    def test_the_rounding_cost_is_split_the_same_way_wherever_it_is_stated(self):
        """A 5e-4 display floor leaves the sensitive tier holding 267 of 287 -- a 7.0%
        total miss. The part the ROUNDING causes is 273 -> 267: 6 sites, 2.1 points.
        "dropped ~7% of validated off-targets" fused the two."""
        cbulge = _read("PostProcess", "cbulge_rows.py")
        full, floored, n_pos = 273, 267, int(self.d["cb_pos"])
        self.assertEqual(f"{(full - floored) / n_pos * 100:.1f}", "2.1")
        self.assertEqual(f"{(n_pos - floored) / n_pos * 100:.1f}", "7.0")
        for name, text in (("METHODS", self.doc), ("cbulge_rows", cbulge)):
            with self.subTest(surface=name):
                self.assertIn(f"{floored} of the {n_pos}", text)
                self.assertNotIn("dropped ~7% of validated off-targets", text)


# --------------------------------------------------------------------------- #
# 5. The bulge-stratified average precisions
# --------------------------------------------------------------------------- #
class TestBulgeAveragePrecisions(_Surfaces):

    def test_each_stratum_ap_is_the_curve_datas_wherever_it_is_quoted(self):
        cb1 = float(self.d["scb_b1_ap"])
        cfd1 = float(self.d["scfd_b1_ap"])
        cfd0 = float(self.d["scfd_b0_ap"])
        cb0 = float(self.d["scb_b0_ap"])
        self.assertOnEvery(f"{cb1:.3f}", "doc", "modal", "changelog")
        self.assertOnEvery(f"{cfd1:.2f}", "doc", "modal")
        self.assertOnEvery(f"{cfd0:.3f}", "doc")
        self.assertOnEvery(f"{cfd0:.2f}", "modal")
        self.assertOnEvery(f"{cb0:.3f}", "changelog")
        # the figure legend computes all of them, so it cannot be the one that drifts
        legend = " ".join(rpf.stratified_legend(self.d, "scb")
                          + rpf.stratified_legend(self.d, "scfd"))
        for ap in (cb1, cfd1, cfd0, cb0):
            self.assertIn(f"AP={ap:.3f}", legend)

    def test_the_superseded_bulge_ap_is_on_no_live_surface(self):
        """0.25 was the leaky stratified panel's bulge=1 AP (0.252 on a support that
        included the 6 training guides). The held-out figure measures 0.285."""
        cb1 = float(self.d["scb_b1_ap"])
        self.assertGreater(cb1, 0.28)
        for bad in ("bulge=1 AP ≈ 0.25", "AP ≈ 0.25)", "CRISPR-Bulge’s 0.25 "):
            self.assertOnNone(bad, "doc", "methods", "readme", "modal")

    def test_the_bulge_ap_is_set_against_crispr_bulges_not_against_the_floor(self):
        """CFD's 0.009 sits ~160x its own bulge=1 prevalence floor, while the same page
        calls CFD's 0.168 "far above" a floor it beats by ~1,379x. "essentially the
        prevalence floor" contradicted that yardstick, so the comparison is against
        CRISPR-Bulge's 0.285 instead -- in the document, the modal and the changelog."""
        floor = int(self.d["scfd_b1_pos"]) / int(self.d["scfd_b1_n"])
        self.assertEqual(f"{floor:.2e}".replace("e-05", "e-5"), "5.75e-5")
        over_floor = float(self.d["scfd_b1_ap"]) / floor
        self.assertTrue(155 < over_floor < 165, f"{over_floor} is not ~160x")
        headline = self.cfd["aupr"]["point"] / self.cfd["prevalence"]
        self.assertEqual(f"{headline:,.0f}", "1,379")
        self.assertIn("160x", self.changelog)
        self.assertIn("1,379x", self.changelog)
        for bad in ("essentially the prevalence floor", "at the floor"):
            self.assertOnNone(bad, "doc", "methods", "readme", "modal")

    def test_the_no_skill_floor_each_page_quotes_is_the_artifacts_prevalence(self):
        self.assertOnEvery(_times_ten(self.cfd["prevalence"], 2), "doc", "modal")
        self.assertIn(_times_ten(self.cb["prevalence"], 2), self.modal)


# --------------------------------------------------------------------------- #
# 6. What the calibration IS -- the claims that travel with every number above
# --------------------------------------------------------------------------- #
class TestCalibrationMethodClaims(_Surfaces):

    def test_no_surface_calls_the_shipped_map_a_leave_one_guide_out_fit(self):
        """p_point is the all-guide isotonic fit; leave-one-guide-out produces the band
        errors, phi and the coverage check. METHODS.md and README.md presented the
        in-sample closure as an out-of-sample guarantee."""
        for name in ("methods", "readme", "modal"):
            with self.subTest(surface=name):
                self.assertNotIn("fit leave-one-guide-out", getattr(self, name))
        self.assertIn("fit on all benchmark guides at once", self.doc)
        for name, text in self.surfaces("methods_card", "readme_card"):
            with self.subTest(surface=name):
                self.assertRegex(text, r"fit on all (46 )?benchmark guides at\s+once")
        # §4's one mention of the phrase is the sentence that disowns it
        for para in re.split(r"\n\s*\n", self.doc):
            if "fit leave-one-guide-out" in para:
                with self.subTest(para=para[:50].replace("\n", " ")):
                    self.assertTrue("Saying the shipped map" in para
                                    or "still called the shipped map" in para)

    def test_the_coverage_range_is_the_artifacts_on_every_surface_that_gives_one(self):
        """89.1% to 100% over the bands the card shows, CFD 0.2-0.5 the only one below
        nominal. "about 90% coverage" and "~90% coverage" rounded the shortfall away."""
        import risk_card
        cov = []
        for scorer, spec in (("cfd", risk_card._CFD_ROWS),
                             ("crispr_bulge", risk_card._CB_ROWS)):
            table = _artifact(scorer)["loo_total_count_coverage_nominal90"]
            cov += [table[label] for label, _d in spec if table.get(label) is not None]
        lo, hi = min(cov) * 100, max(cov) * 100
        self.assertEqual((f"{lo:.1f}", f"{hi:.0f}"), ("89.1", "100"))
        self.assertOnEvery(f"{lo:.1f}%", "doc", "methods_card", "readme_card")
        self.assertOnEvery(f"{lo:.1f}% to {hi:.0f}%", "methods_card", "readme_card",
                           "calib_src")
        self.assertIn(f"{round(lo)}–{round(hi)}% of folds", self.modal)
        for name in ("methods", "readme", "calib_src"):
            with self.subTest(surface=name):
                self.assertNotIn("to about 90%", getattr(self, name))
                self.assertNotIn("~90% coverage", getattr(self, name))

    def test_a_held_out_benchmark_guide_is_never_equated_with_a_novel_one(self):
        """§5 is explicit that the per-held-out-guide check is not the novel-guide case.
        Two runtime modules equated them in a parenthesis.

        The first version of this guard banned the two exact spellings it had just
        deleted, "(= one novel guide)" and "(=one novel guide)". It therefore missed a
        third, "(one novel guide)" with no "=", which survived in fit_calibration.py 118
        lines below the docstring the same commit corrected -- the same shape as the
        bare "+0.66" the test below exists for. So the rule is now about the CLAIM, not
        a spelling: in these two modules every mention of a novel guide has to be on a
        line that denies it, which is the only thing §5 lets them say about one."""
        mentions = re.compile(r"novel[- ]guide|novel one|novel\s+guide")
        denial = re.compile(r"\bNOT\b|\bnot\b|cannot|never")
        for name in ("calib_src", "fit_src"):
            text = getattr(self, name)
            with self.subTest(surface=name):
                self.assertNotIn("(= one novel guide)", text)
                self.assertNotIn("(=one novel guide)", text)
                self.assertIn("novel", text)   # it still says what it is NOT
            lines = [l for l in text.splitlines() if mentions.search(l)]
            self.assertTrue(lines, f"{name} no longer says what the check is NOT")
            for line in lines:
                with self.subTest(surface=name, line=line.strip()[:70]):
                    self.assertRegex(
                        line.strip(), denial,
                        f"{name} mentions a novel guide without denying it: "
                        f"{line.strip()!r}")
        self.assertIn("cannot be, validated for a clinically novel guide", self.modal)

    def test_the_cross_band_correlation_is_stated_once_and_is_the_measured_one(self):
        """"~+0.66" reproduces under no definition anyone could find; the measured
        correlation is +0.91 Pearson / +0.73 Spearman / +0.79 on the LOGO residuals,
        and §7 is the one place it is measured. It outlived its retirement in two
        docstrings, which is exactly the drift this file exists to catch."""
        self.assertIn("+0.91 Pearson", self.doc)
        for name in ("calib_src", "fit_src", "doc"):
            text = getattr(self, name)
            for value in ("+0.91", "+0.73", "+0.79"):
                with self.subTest(surface=name, value=value):
                    self.assertIn(value, text,
                                  f"{name} drops the measured correlation {value}")
            # the retired figure may be NAMED -- in quotes, beside the sentence that
            # retires it -- but never stated bare, which is how it survived two
            # retirements and two documents
            for m in re.finditer(r"\+0\.66", text):
                with self.subTest(surface=name, at=text[max(0, m.start() - 30):m.end()]):
                    self.assertIn('"', text[max(0, m.start() - 6):m.start()],
                                  f"{name} states +0.66 as a live value, unquoted")
            if "+0.66" in text:
                with self.subTest(surface=name):
                    self.assertTrue(
                        any(w in text for w in ("retired", "no definition we can find",
                                                "reproduced under no")),
                        f"{name} names +0.66 without saying it was retired")

    def test_the_phi_and_band_values_the_document_quotes_are_the_artifacts(self):
        """Already covered per-sentence by test_methods_calibration_numbers; here the
        two runtime modules are held to the same artifacts, since they are what the
        report actually applies."""
        for scorer in ("cfd", "crispr_bulge"):
            meta = _artifact(scorer)
            self.assertIn("ge0.2", meta["phi_bands"])
            self.assertEqual(set(meta["band_counts"]), set(meta["phi_bands"]))
        import offtarget_calibration as oc
        summary = oc.artifact_summary("cfd")
        self.assertEqual(summary["n_positives"], self.cfd["n_positives"])
        self.assertEqual(summary["n_sites"], self.cfd["n_sites"])
        self.assertEqual(summary["aupr"], self.cfd["aupr"])


# --------------------------------------------------------------------------- #
# 7. The changelog entry that tells a reader what moved
# --------------------------------------------------------------------------- #
class TestTheShippedExamplePage(_Surfaces):
    """`docs/example-report.html` is the report a reader meets before they run anything,
    and it is GENERATED -- which is how it kept the pre-fix modal for a release and a
    half after every other surface was corrected. It is rebuilt from this tree now, so
    the one check that keeps it honest is that the modal it carries is still the modal
    this tree renders: a later edit to `risk_card` that is not followed by a rebuild
    fails here instead of shipping a page that contradicts the code behind it."""

    @classmethod
    def page(cls):
        return _read("docs", "example-report.html")

    @_only_default
    def test_the_example_page_embeds_the_modal_this_tree_renders(self):
        import risk_card
        import offtarget_calibration as oc
        self.assertIn(risk_card._render_modal(oc), self.page(),
                      "docs/example-report.html no longer carries the modal this tree "
                      "renders; rebuild it (crisprme.py generate-report) -- see the "
                      "changelog entry for how, and backlog item 42")

    def test_the_example_page_embeds_the_committed_figures(self):
        import risk_card
        page = self.page()
        for name in ("pr_curves.png", "pr_bulge_stratified.png"):
            with self.subTest(png=name):
                self.assertIn(risk_card._png_uri(name), page,
                              f"the example page's {name} is not the committed one")

    def test_the_example_page_carries_no_superseded_claim(self):
        page = self.page()
        for bad in ("essentially the prevalence floor", "held-out AUPR",
                    f'{self.cfd["n_positives"]} sites', f'{self.cb["n_positives"]} sites',
                    "roughly 6&times;", "left-and-down", "AP ≈ 0.25"):
            with self.subTest(claim=bad):
                self.assertNotIn(bad, page)


class TestChangelogQuotesTheArtifacts(_Surfaces):

    @_only_default
    def test_the_figure_entrys_after_numbers_are_the_committed_curve_datas(self):
        """Every "now reads" value in the PR-figure entry is re-derived here, so the
        entry cannot describe a figure the repository does not contain."""
        marks = {m[1]: m for m in json.loads(str(self.d["marks"])) if m[0] == "cb"}
        for thr in (cbop.sensitive_threshold(), 0.2, 0.5):
            prec = marks[thr][3]
            with self.subTest(threshold=thr):
                self.assertIn(f"{prec:.6f}", self.changelog)
        self.assertIn(f'{float(self.d["cb_ap"]):.6f}', self.changelog)
        self.assertIn(f'{int(self.d["cb_pos"])} positives of '
                      f'{int(self.d["cb_n"]):,} sites', self.changelog)

    def test_the_stratified_entrys_three_pairs_end_on_the_curve_datas_values(self):
        for key in ("scb_all_ap", "scb_b0_ap", "scb_b1_ap"):
            with self.subTest(curve=key):
                self.assertIn(f'{float(self.d[key]):.3f}', self.changelog)
        self.assertIn(f'{float(self.d["scfd_all_ap"]):.3f}', self.changelog)


class TestTheChangelogReadItself(unittest.TestCase):
    """The helper every assertion in this file depends on.

    It used to return whichever of [Unreleased] / the newest dated section was non-empty,
    newest first. Adding an unrelated entry to [Unreleased] -- a new companion file --
    therefore hid the released section and failed a dozen assertions at once, each saying
    the changelog does not state a number the changelog states.
    """

    def test_it_returns_both_sections(self):
        got = _newest_changelog_section()
        self.assertIn("## [Unreleased]", got)
        self.assertRegex(got, r"## \[\d+\.\d+\.\d+\]")

    def test_an_unrelated_unreleased_entry_does_not_hide_the_released_one(self):
        """The exact break. A number in the dated section must still be found when
        [Unreleased] has content of its own."""
        got = _newest_changelog_section()
        self.assertIn("### Added", got)
        # the released section's own text survives alongside whatever is open
        txt = _read("CHANGELOG.md")
        start = txt.index("## [Unreleased]")
        nxt = txt.index("\n## [", start + 1) + 1
        dated = txt[nxt:txt.index("\n## [", nxt + 1)]
        probe = next((l for l in dated.splitlines() if len(l.strip()) > 40), None)
        self.assertIsNotNone(probe, "the dated section has no prose to probe with")
        self.assertIn(probe, got)

    def test_it_stops_before_the_section_below(self):
        """Both sections, not the whole file: a number stated only in an ancient entry
        must not satisfy an assertion about a current one."""
        txt = _read("CHANGELOG.md")
        self.assertLess(len(_newest_changelog_section()), len(txt) / 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
