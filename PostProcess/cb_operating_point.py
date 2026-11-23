"""CRISPR-Bulge sensitive operating point — single source of truth.

The report's Section-4 panel selection, the Section-4 "sensitive" tier, the top
risk-card band, the PR-curve marked point, and the methods text ALL key on the
value returned here, so they can never drift apart.

Default = the **95%-recall** threshold on the 40-guide held-out external benchmark
(Tsai 2015 + Refined-TrueOT + Chen 2017 + Listgarten 2018; 287 validated
off-targets, CRISPR-Bulge held-out). Point estimates per recall target, each with
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
log10 axis -- ``log_axis`` / ``draw_log_score_axis`` below.
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
