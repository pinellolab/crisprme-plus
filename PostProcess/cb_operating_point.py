"""CRISPR-Bulge sensitive operating point — single source of truth.

The report's Section-4 panel selection, the Section-4 "sensitive" tier, the top
risk-card band, the PR-curve marked point, and the methods text ALL key on the
value returned here, so they can never drift apart.

Default = the **95%-recall** threshold on the 40-guide held-out external benchmark
(Tsai 2015 + Refined-TrueOT + Chen 2017 + Listgarten 2018; 287 validated
off-targets, CRISPR-Bulge held-out). Point estimates per recall target, each with
the guide-cluster bootstrap 95% CI:

    recall  CB threshold   guide-bootstrap 95% CI       sites/true-OT
    90%     2.7e-3         [3.4e-4, 1.4e-2]             ~9
    95%     1.9e-4         [5.8e-6, 3.1e-3]   <-DEFAULT ~21
    99%     1.5e-8         [7.1e-10, 3.2e-4]            ~354

The tail quantile is intrinsically noisy (wide CI), so this is the documented
default operating point, not a hard physical constant. Higher recall buys recall
only by inflating the validation list: 95% -> 99% is ~6x more sites to assay
(ACTGA/NRG worked example: 206 -> 1,286 loci). Select a different target with
the ``CRISPRME_CB_RECALL_TARGET={90|95|99}`` environment variable.

NOTE: the default threshold (1.9e-4) is BELOW the legacy 3-decimal display floor
(5e-4), so CRISPR-Bulge scores must be written at full precision ("{:.6g}") for
this operating point to be representable (see new_simple_analysis.preprocess_
CRISPR_BULGE_score / analisi_indels_NNN).
"""
import os

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
