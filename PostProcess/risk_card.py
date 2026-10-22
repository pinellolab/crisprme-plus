#!/usr/bin/env python3
"""Top-of-report OFF-TARGET RISK card + "How to read this" modal.

Replaces the legacy ``100/(100+Σ)`` aggregate specificity rows. For each scorer
(CFD, CRISPR-Bulge -- shown side by side, never blended) it reports the count of
<=1-gap off-target sites per score band and a calibrated EXPECTED number of true
off-targets with a coverage-valid 90% range, plus a qualitative note that a
deployed guide's own prior shift can move the count further. The numbers come
from ``offtarget_calibration`` (isotonic point map + two-level
guide-cluster bootstrap + Gamma overdispersion; see METHODS). Degrades
gracefully to counts-only if the calibration artifact/module is unavailable.

Kept out of ``generate_report.py`` to keep that file's blast radius small; the
caller extracts the <=1-gap off-target score arrays (it owns the df/cols idioms)
and passes them here. The modal reuses the observed-haplotypes modal CSS/JS
(``oh-open`` link + ``oh-modal`` div), so no new script is needed.
"""
from __future__ import annotations
import html as _html
import os
import sys
from typing import List, Optional, Sequence

MODAL_ID = "risk-card-modal"

# report bands the card surfaces (label in artifact -> display). Order high->low.
_CFD_ROWS = [("0.5-1", "CFD ≥ 0.5"), ("0.2-0.5", "CFD 0.2–0.5"), ("ge0.2", "CFD ≥ 0.2")]
_CB_ROWS = [("0.5-1", "CRISPR-Bulge ≥ 0.5"), ("0.2-0.5", "CRISPR-Bulge 0.2–0.5"),
            ("ge0.2", "CRISPR-Bulge ≥ 0.2")]
_MIN_BAND = 0.2  # lowest score edge the card displays; scores below are pre-filtered out


def _esc(x) -> str:
    return _html.escape("" if x is None else str(x))


def _calib():
    """Import the calibration module, or None if unavailable."""
    try:
        import offtarget_calibration as cal  # same package dir
        return cal
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"risk-card: calibration unavailable ({exc}); counts-only\n")
        return None


def _band_map(rows):
    return {r["band"]: r for r in rows}


def _fmt_ci(lo, hi):
    return f"[{lo:.0f}, {hi:.0f}]"


def _scorer_block(cal, scorer, scores, row_spec):
    """Return (html_rows, ge02_point) for one scorer. Falls back to counts-only (no
    expected-true column) if the calibration module OR its artifact is unavailable."""
    import numpy as np
    scores = [s for s in scores if s is not None]
    res = None
    if cal is not None and scores:
        try:
            res = _band_map(cal.expected_counts(scores, scorer))
        except Exception as exc:  # noqa: BLE001 - missing artifact etc. -> counts-only
            sys.stderr.write(f"risk-card: calibration failed for {scorer} ({exc}); counts-only\n")
            res = None
    if res is None:
        a = np.asarray(scores, dtype=float)
        rows = []
        for label, disp in row_spec:
            if label == "0.5-1":
                n = int((a >= 0.5).sum()) if a.size else 0
            elif label == "0.2-0.5":
                n = int(((a >= 0.2) & (a < 0.5)).sum()) if a.size else 0
            else:
                n = int((a >= 0.2).sum()) if a.size else 0
            rows.append(f'<tr><td>{_esc(disp)}</td><td>{n:,}</td><td>&mdash;</td></tr>')
        return "".join(rows), None
    rows = []
    ge02_point = None
    for label, disp in row_spec:
        r = res.get(label)
        if not r:
            continue
        n = r["n_sites"]
        pt = r["point"]
        ci = _fmt_ci(r["lo90"], r["hi90"])
        rows.append(
            f'<tr><td>{_esc(disp)}</td><td>{n:,}</td>'
            f'<td>≈{pt:.0f} &nbsp;<span style="color:#64748b">{_esc(ci)}</span></td></tr>')
        if label == "ge0.2":
            ge02_point = pt
    return "".join(rows), ge02_point


def _pr_curve_uri() -> Optional[str]:
    """Return a data: URI for the bundled PR-curve PNG, or None if absent."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration", "pr_curves.png")
    if not os.path.isfile(path):
        return None
    import base64
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


# clinician-reviewed modal copy (see METHODS; signed statistician/clinician/genome-editor)
_ISOTONIC = (
    "Isotonic calibration is a simple translation step: it takes a site's raw CFD or "
    "CRISPR-Bulge score, which on its own is only a ranking number, and converts it into an "
    "estimated chance the site is a true off-target, with one rule: a higher score never maps "
    "to a lower probability. Each scorer's curve was learned from experimentally confirmed "
    "off-targets across four published studies (359 for CFD, 287 for CRISPR-Bulge; see "
    "METHODS §2), built so no guide ever set its own probabilities. The “expected true” "
    "count is the sum of these per-site probabilities; "
    "its 90% range is a guide-level bootstrap that reaches ~90% coverage for a single novel "
    "guide. Both scorers are sequence-based, so the calibration applies to any ≤1-gap "
    "alignment (reference or variant allele); ≥2-gap sites and per-sample carrier incidence "
    "are out of scope.")
_READ = (
    "Two independent scorers, shown separately, never blended. <b>CFD</b> is a broad "
    "screen (flags almost every real off-target, most flags false). <b>CRISPR-Bulge</b> is a "
    "prioritizer (high scores are usually real, but it confirms only about half of "
    "assay-detected off-targets and is defined only on ≤1-gap alignments). A low "
    "CRISPR-Bulge count is <b>not</b> clearance. The expected count is an estimate, not a "
    "ceiling or a shortlist: every listed site remains an assay candidate. The expected count "
    "assumes this guide is about as off-target-prone as a typical benchmark guide; a guide "
    "that is genuinely more or less prone at its own locus can shift the count substantially "
    "in either direction, and that prior shift is usually a larger uncertainty than the "
    "statistical range shown (it is not captured by that range).")


def render_risk_card(meta: dict, cfd_scores: Sequence[float],
                     cb_scores: Optional[Sequence[float]], n_perfect: int = 0):
    """Return (card_html, modal_html). ``*_scores`` are the <=1-gap OFF-target raw
    scores (on-targets excluded); ``cb_scores`` may be None if CRISPR-Bulge did not run."""
    cal = _calib()
    guide = meta.get("guide_display") or meta.get("guide") or ""
    header = f"{_esc(guide)} &middot; {_esc(meta.get('nuclease',''))} &middot; {_esc(meta.get('pam',''))}"

    # NB: risk-card tables use plain header ROWS (no <thead>/<tbody>/<th>), so the card never
    # adds a table the report's top-1000 parser could mistake for the ottable. (That parser is
    # anchored on class="ottable" and bounded at the following </table>; observed-haplotype
    # modals also emit <thead>/<tbody> tables, so the card deliberately stays out of that markup.)
    _hdr = ('<tr class="rc-hdr"><td>band (≤1-gap)</td><td>sites</td>'
            '<td>expected true (90%)</td></tr>')
    # the card shows only >=0.2 bands; pre-filter so the calibrator does not bootstrap the
    # huge near-zero mass (performance) -- all displayed bands start at 0.2.
    _cfd = [s for s in (cfd_scores or []) if s is not None and s >= _MIN_BAND]
    cfd_rows, cfd_pt = _scorer_block(cal, "cfd", _cfd, _CFD_ROWS)
    blocks = [f'<div class="rc-scorer"><table class="rc-tbl">{_hdr}{cfd_rows}</table></div>']
    if cb_scores is not None:
        _cb = [s for s in cb_scores if s is not None and s >= _MIN_BAND]
        cb_rows, _ = _scorer_block(cal, "crispr_bulge", _cb, _CB_ROWS)
        blocks.append(f'<div class="rc-scorer"><table class="rc-tbl">{_hdr}{cb_rows}</table></div>')

    prev = ""
    if cal is not None and cfd_pt is not None:
        prev = ('<div class="rc-prev">The expected count assumes this guide is about as '
                'off-target-prone as a typical benchmark guide. A guide that is genuinely '
                'more or less prone at its own locus can shift this count substantially in '
                'either direction &mdash; often the largest uncertainty, and not something '
                'the range above captures.</div>')

    warn = ""
    if n_perfect and n_perfect >= 2:
        warn = (f'<div class="rc-warn">⚠ On-target is not unique: {int(n_perfect)} '
                f'perfect (0-mismatch, 0-gap) matches.</div>')

    note = ('<div class="rc-note">'
            '<b>Expected true</b> estimates how many of the listed sites are real off-targets '
            '&mdash; it is not a shortlist; every listed site stays an assay candidate (see the '
            'validation panel). Sites needing &ge;2 bulges show &ldquo;&ndash;&rdquo; (not scored, '
            'not no-risk) and are kept but demoted there. The two scorers are different views of '
            'the same sites, not additive. '
            f'<a class="rc-open" href="#{MODAL_ID}">How to read this →</a></div>')

    card = (
        f'<style>{CARD_CSS}</style>'
        '<div class="risk-card">'
        f'<div class="rc-head">OFF-TARGET RISK &middot; {header}</div>'
        f'{warn}'
        f'<div class="rc-scorers">{"".join(blocks)}</div>'
        f'{prev}{note}</div>')

    modal = _render_modal(cal)
    return card, modal


def _render_modal(cal) -> str:
    summary_html = ""
    if cal is not None:
        def _fmt(scorer):
            # per-scorer + None-guarded: a missing/corrupt artifact for ONE scorer must
            # drop only that row, not both, and a None coverage/CI must not raise.
            try:
                a = cal.artifact_summary(scorer)
                ap = a.get("aupr") or {}
                aupr = ap.get("point")
                if aupr is None:
                    return ""
                ci = ap.get("ci") or [None, None]
                cov = (a.get("loo_total_count_coverage_nominal90") or {}).get("ge0.2")
                ci_txt = (f'[{ci[0]:.2f}–{ci[1]:.2f}]'
                          if ci[0] is not None and ci[1] is not None else "")
                cov_txt = f"{int(cov * 100)}%" if cov is not None else "—"
                return (f'<tr><td>{_esc(a.get("scorer"))}</td>'
                        f'<td>{aupr:.3f} <span style="color:#64748b">{ci_txt}</span></td>'
                        f'<td>{cov_txt}</td></tr>')
            except Exception:  # noqa: BLE001 - one scorer's artifact issue must not drop the other
                return ""
        rows = _fmt("cfd") + _fmt("crispr_bulge")
        if rows:
            summary_html = (
                '<table class="rc-tbl"><tr class="rc-hdr"><td>scorer</td>'
                '<td>held-out AUPR (90% CI)</td><td>≥0.2 single-guide coverage</td></tr>'
                f'{rows}</table>')
    aupr_note = (
        '<p style="font-size:0.85em;color:#64748b">No-skill AUPR baseline = the off-target '
        'prevalence ≈ 1.22×10⁻⁴; both scorers are far above it. The two AUPRs sit on different '
        'held-out benchmark sets (CFD 359 sites / 46 guides, CRISPR-Bulge 287 / 40) and are '
        'not directly comparable across scorers.</p>'
        if summary_html else "")
    pr = _pr_curve_uri()
    pr_html = (f'<div class="rc-pr"><img alt="Precision-recall of each scorer" src="{pr}" '
               f'style="max-width:100%"></div>') if pr else ""
    # CSS-only (:target) modal -- no JavaScript, so the report stays self-contained.
    return (
        f'<div class="rc-modal" id="{MODAL_ID}"><div class="rc-modal-box">'
        '<a href="#" class="rc-modal-close" aria-label="Close">&times;</a>'
        '<div class="rc-modal-head">How to read the off-target risk</div>'
        '<div class="rc-modal-body">'
        f'<p>{_READ}</p>'
        f'{summary_html}{aupr_note}'
        f'{pr_html}'
        f'<p><b>What “expected true” means.</b> {_ISOTONIC}</p>'
        '<p style="color:#64748b">Full method: METHODS § off-target risk summary.</p>'
        '</div></div></div>')


# minimal CSS; injected once by the caller alongside the existing report styles
CARD_CSS = """
.risk-card{border:1px solid #e2e8f0;border-radius:8px;padding:14px 16px;margin:10px 0;background:#fff}
.risk-card .rc-head{font-weight:700;margin-bottom:8px}
.risk-card .rc-warn{color:#b91c1c;font-weight:600;margin:6px 0}
.risk-card .rc-scorers{display:flex;gap:24px;flex-wrap:wrap}
.risk-card .rc-scorer{flex:1;min-width:280px}
.risk-card table.rc-tbl{border-collapse:collapse;width:100%;font-size:0.85em}
.risk-card table.rc-tbl td{text-align:left;padding:3px 8px;border-bottom:1px solid #f1f5f9}
table.rc-tbl .rc-hdr td{font-weight:600;color:#334155}
.risk-card .rc-prev{margin-top:8px;font-size:0.85em}
.risk-card .rc-note{margin-top:6px;font-size:0.85em;color:#475569}
.rc-open{cursor:pointer;color:#2c5282;text-decoration:underline}
.rc-modal{display:none;position:fixed;inset:0;z-index:60;background:rgba(20,30,50,.45);align-items:flex-start;justify-content:center;padding:4vh 12px;overflow:auto}
.rc-modal:target{display:flex}
.rc-modal-box{background:#fff;border-radius:10px;max-width:900px;width:100%;padding:18px 22px;position:relative;box-shadow:0 8px 40px rgba(20,40,80,.3)}
.rc-modal-close{position:absolute;top:8px;right:14px;font-size:1.5em;text-decoration:none;color:#64748b}
.rc-modal-head{font-weight:700;font-size:1.1em;margin-bottom:10px}
.rc-modal-body p{margin:.6em 0;font-size:0.85em;line-height:1.5}
.rc-modal-body table.rc-tbl{margin:10px 0}
"""
