"""CRISPR-Bulge sensitive operating point — single source of truth.

The report's Section-4 panel selection, the Section-4 "sensitive" tier, the top
risk-card band, the PR-curve marked point, and the methods text ALL key on the
value returned here, so they can never drift apart.

Default = the **95%-recall** threshold on the 40-guide held-out external benchmark
(Tsai 2015 + Refined-TrueOT + Chen 2017 + Listgarten 2018; 287 validated
positives, CRISPR-Bulge held-out; see WHICH POPULATION below -- 40 of those 287 are
on-targets). Point estimates per recall target, each with
the guide-cluster bootstrap 95% CI:

    recall  CB threshold   guide-bootstrap 95% CI         list   sites/true-OT
    90%     2.7e-3         [3.4e-4, 1.4e-2]              2,343      9.0
    95%     1.9e-4         [5.8e-6, 3.1e-3]   <-DEFAULT  6,048     22.2
    99%     1.5e-8         [7.1e-10, 3.2e-4]            99,406    348.8

The last two columns are the list each threshold produces on that benchmark and how
many sites it lists per validated off-target it finds, both at the EXACT recall
quantiles. At the rounded thresholds printed here the 95% list is 6,076 (273 of the
287, 22.3 listed per true one) and the 99% list 99,477. Every number here is
re-derived from the committed curve data (``calibration/pr_figure_data.npz``) by
``PostProcess/test_shared_numbers.py``, which also holds the methods document and the
risk-card modal to the same values, so this table cannot drift from them -- except
99,477, which that data BOUNDS rather than gives (it is an interior point of a run the
stored polyline decimates).

WHICH POPULATION THIS RECALL IS AGAINST. The 287 are the held-out validated POSITIVES,
and they include each guide's own ON-TARGET: 40 perfect (0-mismatch, 0-gap) rows, one
per guide. The deployed ranking never contains one -- the report builds every
CRISPR-Bulge array as "off-targets (exclude perfect)" -- so the default above is a
quantile over a population 40 of whose members it is never applied to. Measured both
ways on the same benchmark, and under the best-alignment view the report's own ranked
tables use:

    population / alignment              recall of this default   its own 95% point
    287 positives, assay alignment      273/287 = 95.1%          1.9e-4 (the default)
    247 off-targets, assay alignment    233/247 = 94.3%          3.8e-5 (5.1x lower)
    247 off-targets, best-CB alignment  241/247 = 97.6%          1.1e-3 (5.7x higher)

(guide-cluster 95% CI on that 94.3%: 0.902 to 0.994.) The 99% row of the table above is
the SAME number on both populations: ceil(0.99*287) and ceil(0.99*247) each permit 2
misses and the on-targets occupy the top of the ranking, so the index shifts by 40 and
lands on the same score. At the 95% target the on-targets buy TWO permitted misses, not
forty -- 14 against 12. The default STAYS at 1.9e-4: the correction has no determinate
sign (5.1x down on the assay alignment, 5.7x up on the best-alignment view, and higher
again if each guide is weighted equally instead of each site), and each candidate point
lies well inside the other's bootstrap interval -- P(resampled 287-point <= 3.8e-5) =
0.375, P(resampled 247-point >= 1.9e-4) = 0.447, the two being 0.705 decades apart
inside an interval 2.72 decades wide. Written out for a reader in
``docs/offtarget_risk_summary_METHODS.md`` section 3.5.

The tail quantile is intrinsically noisy (wide CI), so this is the documented
default operating point, not a hard physical constant. Higher recall buys recall
only by inflating the validation list: on this benchmark 95% -> 99% is ~16x more
sites to assay (6,048 -> 99,406, i.e. 16.4x). A single search's own ratio is a
different quantity and need not match it; the worked example lives in exactly one
place, ``docs/offtarget_risk_summary_METHODS.md`` section 3.5, and is not restated
here.
Select a different target with the ``CRISPRME_CB_RECALL_TARGET={90|95|99}``
environment variable.

NOTE: the default threshold (1.9e-4) is BELOW the legacy 3-decimal display floor
(5e-4), so CRISPR-Bulge scores must be written at full precision ("{:.6g}") for
this operating point to be representable (see new_simple_analysis.preprocess_
CRISPR_BULGE_score / analisi_indels_NNN) -- ``score_text`` below -- and SHOWN on a
log10 axis -- ``log_axis`` / ``draw_log_score_axis`` below. That axis is also why a score
panel's legends go underneath it instead of in it -- ``legend_below`` below.
"""
import math
import os
from collections import namedtuple

# recall target -> (CB threshold, bootstrap 95% CI lo, hi)
CB_RECALL_POINTS = {
    90: (2.7e-3, 3.4e-4, 1.4e-2),
    95: (1.9e-4, 5.8e-6, 3.1e-3),
    99: (1.5e-8, 7.1e-10, 3.2e-4),
}
DEFAULT_TARGET = 95


def recall_target() -> int:
    """Selected recall target (default 95); overridable via CRISPRME_CB_RECALL_TARGET."""
    try:
        t = int(os.environ.get("CRISPRME_CB_RECALL_TARGET", DEFAULT_TARGET))
    except (TypeError, ValueError):
        return DEFAULT_TARGET
    return t if t in CB_RECALL_POINTS else DEFAULT_TARGET


def sensitive_threshold() -> float:
    """CRISPR-Bulge score at/above which a site is hard-included (the sensitive default)."""
    return CB_RECALL_POINTS[recall_target()][0]


def _thr_str() -> str:
    return f"{sensitive_threshold():.2g}"


def sensitive_label_html() -> str:
    """e.g. '&ge; 0.00019 (95% recall)' for report tier labels / panel note."""
    return f"&ge; {_thr_str()} ({recall_target()}% recall)"


def sensitive_label_text() -> str:
    """e.g. '≥ 0.00019 (95% recall)' for the risk-card band + plain-text copy."""
    return f"≥ {_thr_str()} ({recall_target()}% recall)"


# --------------------------------------------------------------------------- #
# How a CRISPR-Bulge score is WRITTEN
# --------------------------------------------------------------------------- #
# The one text of every place that writes a CRISPR-Bulge score for a reader. The
# post-analyses already write the search's own tables this way (cbulge_rows.score_rows,
# targets_scoring: "{:.6g}" / "%.6g" since 2.8.0), so the report's curated tables
# re-formatting the value they read with SCORE_FORMAT reproduces the source text rather
# than rounding it a second time.
SCORE_FORMAT = "{:.6g}"


def score_text(score) -> str:
    """A CRISPR-Bulge score as its writers print it: 6 significant digits.

    Six SIGNIFICANT digits, not fixed decimals: the operating points (1.9e-4 at the
    default 95% recall, 1.5e-8 at 99%) and nearly every score sit below the floor of a
    fixed-decimal format, which collapses them all onto one string -- 4 decimals print
    the same "0.0002" on both sides of 1.9e-4, and "0.0000" for everything below 5e-5.
    "{:.6g}" keeps a plain decimal where one is readable (0.00019, 0.0352168) and turns
    to scientific notation only below 1e-4 (2.35e-39), where no reader was going to count
    zeros anyway; 6 digits are inside the precision of the model's own float32 output
    (~7 significant decimal digits), so no digit printed is numerical noise.

    A true zero reads "0" (not "0.0000"), and the not-scored sentinel is not this
    function's business: the caller keeps its own placeholder for an unscored row.
    """
    return SCORE_FORMAT.format(float(score))


# --------------------------------------------------------------------------- #
# How a CRISPR-Bulge score is DRAWN
# --------------------------------------------------------------------------- #
# A linear 0-1 y-axis cannot show this score family: the sensitive operating point is
# 1.9e-4 and whole genome-wide runs sit below it, so every point collapses onto the
# bottom rule and the threshold is indistinguishable from zero. The score panels put
# CRISPR-Bulge on a log10 axis instead, with the operating point drawn on it.
#
# A log axis cannot show 0, and real runs have plenty of exact zeros (the model's float32
# output underflows: in a 2.10.0 genome-wide run about 4% of rows scored exactly 0 and the
# smallest non-zero score was 2.35e-39). Those rows are NOT dropped: the axis gets a FLOOR
# one full decade below the smallest score it has to show, the zeros are drawn on it, and
# a shaded band plus a legend entry say that the band means "0", so a reader never mistakes
# a floored point for a real small score and never silently loses one.
LOG_FLOOR_DECADES = 1
# The "= 0" lane: from a decade below the floor (the axis' bottom) to half a decade above
# it. The lowest real score is a full decade above the floor, so it stays out of the lane,
# and the lane is tall enough to read as a lane on an axis tens of decades tall.
_BAND_HALF_WIDTH = 10.0 ** 0.5
# A guard, never reached by a real score (float32's smallest subnormal is 1.4e-45): a
# floor below float64's own smallest subnormal would be 0, which no log axis can take.
_MIN_FLOOR = 1e-300
THRESHOLD_COLOR = "#b7791f"
ZERO_BAND_COLOR = "#6b7280"

# floor: where a 0 is drawn; bottom/top: the axis limits; band_top: the top of the "= 0" band
LogAxis = namedtuple("LogAxis", "floor bottom top band_top")


def _finite_scores(values):
    """The values in [0,1] (the score's range): the -1 not-scored sentinel, NaN and
    anything unparseable are not scores and must not set an axis bound."""
    out = []
    for v in values:
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if math.isfinite(f) and 0.0 <= f <= 1.0:
            out.append(f)
    return out


def log_axis(values) -> LogAxis:
    """The log10 y-axis for these CRISPR-Bulge scores: floor, limits and zero band.

    The floor is one full decade below the smallest value the axis must show -- the
    smallest NON-ZERO score, or the sensitive threshold when every score is above it, so
    the operating-point line is always inside the axis. Zeros (and anything at or below
    the floor) are drawn ON the floor, inside the band from ``bottom`` to ``band_top``.
    """
    scores = [v for v in _finite_scores(values) if v > 0.0]
    low = min(scores) if scores else sensitive_threshold()
    low = min(low, sensitive_threshold())
    floor = 10.0 ** (math.floor(math.log10(low)) - LOG_FLOOR_DECADES)
    floor = max(floor, _MIN_FLOOR)
    return LogAxis(floor=floor, bottom=floor / 10.0 ** LOG_FLOOR_DECADES, top=1.0,
                   band_top=floor * _BAND_HALF_WIDTH)


def sensitive_plot_label() -> str:
    """e.g. 'sensitive threshold 0.00019 (95% recall)' for a panel's legend."""
    return f"sensitive threshold {_thr_str()} ({recall_target()}% recall)"


def zero_band_label(floor) -> str:
    """e.g. 'score 0 (underflow), drawn at 1e-40' -- what the shaded bottom band means."""
    return f"score 0 (underflow), drawn at {floor:.0e}"


def draw_log_score_axis(ax, values, zeros=False):
    """Put ``ax`` on the log10 CRISPR-Bulge y-axis and draw the operating point on it.

    Shared by the report's Section-2 panels (generate_report._cfd_style_scatter) and the
    results page's / Personal Risk Card's top plots (CRISPRme_plots,
    CRISPRme_plots_personal), so the three can never drift apart. ``values`` is every
    score the panel draws; ``zeros`` says whether any of them is a 0 that will be drawn on
    the floor -- only then is the band shaded and named, so a panel without zeros gets no
    band to explain.

    Returns ``(LogAxis, legend_handles)``; the caller plots its points clipped to
    ``axis.floor`` and adds the handles to its own legend.
    """
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches

    axis = log_axis(values)
    ax.set_yscale("log")
    ax.set_ylim(axis.bottom if zeros else axis.floor, axis.top)
    handles = [mlines.Line2D([], [], color=THRESHOLD_COLOR, linestyle="--", linewidth=1.0,
                             label=sensitive_plot_label())]
    ax.axhline(sensitive_threshold(), color=THRESHOLD_COLOR, linestyle="--",
               linewidth=1.0, zorder=0)
    if zeros:
        ax.axhspan(axis.bottom, axis.band_top, color=ZERO_BAND_COLOR, alpha=0.22,
                   linewidth=0, zorder=0)
        # the band is a thin strip of an axis tens of decades tall and the markers drawn
        # in it are taller than it is, so its top edge is also a line: THAT is what reads
        # as "everything in this lane is a zero" on a printed page
        ax.axhline(axis.band_top, color=ZERO_BAND_COLOR, linewidth=0.8, alpha=0.6,
                   zorder=0)
        handles.append(mpatches.Patch(color=ZERO_BAND_COLOR, alpha=0.22,
                                      label=zero_band_label(axis.floor)))
    return axis, handles


# --------------------------------------------------------------------------- #
# Where a panel's legends GO
# --------------------------------------------------------------------------- #
# Both legends of a score panel used to sit INSIDE the axes ("upper center" for the
# allele-frequency size key, "upper right" for the colours), which worked while the score
# was CFD on a linear 0-1 axis: a top-1000 CFD curve falls away from the top-right corner,
# so the boxes landed on empty paper. On the CRISPR-Bulge log10 axis it is the other way
# round -- the axis reaches down tens of decades to the floor, so every well-scoring site
# is pressed into the top band of the panel, exactly under those two boxes. Measured on a
# real genome-wide run's top 1000 (axis 1e-31..1): the colour legend alone covered 1481 of
# the 2000 drawn points, scores 1.5e-9 to 1.1e-3 -- the whole neighbourhood of the
# sensitive operating point, the one thing the panel exists to show. Adding the threshold
# and "= 0" entries to that box had made it taller still.
#
# So the legends go BELOW the axes, in one row: the size key on the left (it needs its own
# title and its markers read as a scale), the colours/lines on the right. Below, not
# beside: the panels are wide and short (8.0 x 3.2 in in the report, 7.5 x 2.25 in in the
# plot scripts) and the personal report puts two of them side by side in one row, so
# horizontal room is what the data is shortest of -- a right-hand column cost 32% of the
# data area's width (7.08 -> 4.79 in), against 0.37 in of figure height for a row
# underneath, and the longest entry ("sensitive threshold ... (95% recall)") -- one of the
# two the log axis added -- is what sets that column's width.
# The row is placed under the x-label and the FIGURE IS GROWN to fit it, so the data area
# keeps the size it had: nothing is traded for the move except page height.
PANEL_LEGEND_GAP = 0.12     # inches between the x-label and the legend row
PANEL_LEGEND_MARGIN = 0.04  # inches of figure left below the legend row


def _renderer(fig):
    """A renderer to measure artists with, without forcing a full raster draw."""
    try:
        return fig.canvas.get_renderer()
    except AttributeError:  # pragma: no cover - non-Agg canvas
        fig.canvas.draw()
        return fig.canvas.get_renderer()


def _axes_box_inches(ax, fig):
    """The axes' box as (x0, y0, w, h) in inches from the figure's bottom-left."""
    pos = ax.get_position()
    w, h = fig.get_size_inches()
    return pos.x0 * w, pos.y0 * h, pos.width * w, pos.height * h


def legend_below(ax, size_handles=(), series_handles=(), size_title="Allele frequency",
                 fontsize=None, series_ncol=2):
    """Draw this score panel's legends BELOW the axes, out of the data.

    The ONE placement of every score panel the report and the results page draw
    (generate_report._cfd_style_scatter, CRISPRme_plots and CRISPRme_plots_personal's
    plot_with_* functions), so the panels of one report cannot drift apart and no future
    panel re-learns why "upper right" does not work on a log axis.

    ``size_handles`` is the allele-frequency marker-size key (it keeps its own box and its
    ``size_title``: its markers are a scale, not a series); ``series_handles`` is
    everything else the panel names -- Reference, Alternative (variant), any perfect-match
    star, and on the log axis the sensitive-threshold line and the "= 0" underflow band
    (draw_log_score_axis's handles). Either may be empty.

    Lays the axes out first (``tight_layout``, so the caller must NOT call it again),
    places the legends just under the x-label, then GROWS the figure by what they need and
    moves the axes up by the same amount, so the data area keeps the size it had and the
    legends cannot be clipped -- including by a ``savefig`` without ``bbox_inches="tight"``.
    That last step moves ``ax`` alone, so the figure is expected to hold just this panel,
    which is what every score panel is.

    Returns the Legend artists, in the order drawn.
    """
    fig = ax.figure
    fig.tight_layout()  # the axes' final box, with its title/labels and NO legend in it
    rend = _renderer(fig)
    box = _axes_box_inches(ax, fig)
    # top of the legend row, in axes coordinates: under everything the axes already draw
    # below themselves (tick labels, x-label), which is where the axes' tight box ends
    axbb = ax.get_window_extent(rend)
    below = (ax.get_tightbbox(rend).y0 - axbb.y0) / axbb.height
    y = below - PANEL_LEGEND_GAP * fig.dpi / axbb.height
    common = dict(fontsize=fontsize, frameon=False, borderaxespad=0.0)
    legends = []
    if len(size_handles):
        legends.append(ax.legend(handles=list(size_handles), title=size_title,
                                 ncol=len(size_handles), loc="upper left",
                                 bbox_to_anchor=(0.0, y), **common))
        # add_artist is how an axes keeps a SECOND legend -- and it clips what it adds to
        # the axes' own patch, which makes matplotlib skip the artist in every layout
        # calculation: the figure's tight bounding box then stops at the axes and
        # ``savefig(bbox_inches="tight")`` cuts this legend off the bottom of the image
        # (it did: the CFD panel's size key lost its markers, which hung lower than the
        # two-entry colour legend). It is deliberately outside the axes, so it is
        # deliberately not clipped to them.
        ax.add_artist(legends[-1])
        legends[-1].set_clip_on(False)
        legends[-1].set_clip_path(None)
    if len(series_handles):
        # right of the size key when there is one, centred when the panel has only this
        legends.append(ax.legend(
            handles=list(series_handles), ncol=series_ncol,
            loc="upper right" if len(size_handles) else "upper center",
            bbox_to_anchor=(1.0 if len(size_handles) else 0.5, y), **common))
    if not legends:
        return legends
    # room for the row: grow the figure DOWNWARD and keep the axes' box in inches, so the
    # data area is not squeezed to pay for the legends (the panels are only ~3 in tall)
    low = min(leg.get_window_extent(rend).y0 for leg in legends) / fig.dpi
    deficit = PANEL_LEGEND_MARGIN - low
    if deficit > 0:
        w, h = fig.get_size_inches()
        new_h = h + deficit
        fig.set_size_inches(w, new_h, forward=False)
        x0, y0, bw, bh = box
        ax.set_position([x0 / w, (y0 + deficit) / new_h, bw / w, bh / new_h])
    return legends
