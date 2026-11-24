#!/usr/bin/env python3
"""Render the two precision-recall PNGs the off-target risk card embeds.

``PostProcess/risk_card.py`` embeds ``calibration/pr_curves.png`` and
``calibration/pr_bulge_stratified.png`` as data URIs. Both used to be produced by
one-off scripts kept outside the repository, which is how they drifted from the
analysis they illustrate. This script is the in-tree generator, and it runs in two
stages so the second one can be re-run (and tested) by anyone with a checkout:

  stage 1 ``build``   needs the external benchmark (hundreds of MB) and the
                      CRISPR-Bulge weights. It re-scores every held-out site at
                      full precision, computes the precision-recall curves over
                      the WHOLE negative set, and writes a small curve-data file
                      (``pr_figure_data.npz``, a few hundred points per curve)
                      that carries every number the figures print plus the md5 of
                      each input.

  stage 2 ``render``  needs numpy + matplotlib only. It draws both PNGs from the
                      committed curve data. ``PostProcess/test_pr_figures.py``
                      runs this stage and compares the result with the committed
                      PNGs, so the figures can never again drift from the numbers
                      they claim to show.

Curve data is stored losslessly for plotting *and* for the average precision:
points strictly inside a vertical run (same recall, falling precision) are
collinear with the run's endpoints, so dropping them changes neither the rendered
polyline nor the average precision, which is a sum of ``delta_recall * precision``.

Stage 1, as used for the committed PNGs (nice -n 19, ~11 min on CPU)::

    CBULGE_REPO=<crispr-bulge-checkout> \\
    python3 scripts/render_pr_figures.py build \\
        --benchmark benchmark.tsv --scored scored.tsv \\
        --out PostProcess/calibration/pr_figure_data.npz \\
        --diagnostics pr_figure_diagnostics.json

Stage 2 (deterministic, no RNG anywhere -- the old 300k-negative sample is gone)::

    python3 scripts/render_pr_figures.py render \\
        --data PostProcess/calibration/pr_figure_data.npz \\
        --outdir PostProcess/calibration

Supports of the four panels, which are NOT the same and must not be read as if
they were:

  * CFD is not trained on anything, so every benchmark guide is fair game for it:
    both CFD panels use all sites / all guides, which is the support the shipped
    ``cfd_calibration.npz`` AUPR is computed on.
  * CRISPR-Bulge was trained on 6 of the benchmark's guides, so both CRISPR-Bulge
    panels use ``cb_heldout == 1`` only. (The benchmark's ``in_domain`` column is
    *not* a leakage flag -- it is ``n_bulge <= 1``, true for every row of a
    <=1-gap benchmark -- so it is ignored here; ``cb_heldout`` is the flag that
    means what "held out" means.)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_POSTPROCESS = os.path.join(os.path.dirname(_HERE), "PostProcess")
if _POSTPROCESS not in sys.path:
    sys.path.insert(0, _POSTPROCESS)
import cb_operating_point as cbop  # noqa: E402  (single source of truth for the marked point)

CFD_COLOR = "#2b6cb0"
CB_COLOR = "#c05621"
# (label, n_bulge or None for "all", colour, linestyle, linewidth)
STRATA = (("all ≤1-gap", None, "#2b6cb0", "-", 2.0),
          ("bulge=0", 0, "#38a169", "--", 1.6),
          ("bulge=1", 1, "#e53e3e", "-", 2.4))
CFD_MARKS = (0.2, 0.5)

# Caption of pr_curves.png. Two things it used to get wrong, both of which a reviewer
# reads straight off the image without opening the modal:
#   * it called BOTH panels "held-out". Nothing is held out for CFD -- CFD has no
#     training set, so its panel is all 46 benchmark guides, the 6 CRISPR-Bulge trained
#     on included. "External" is the word that is true of both panels; the CRISPR-Bulge
#     panel's own title still says "held-out", because that one is.
#   * it did not say that the benchmark's positives include each guide's own ON-TARGET
#     (46 perfect matches of CFD's 359, 40 of CRISPR-Bulge's 287), while the report's
#     ranked list excludes them. The AUPRs in the panel titles are therefore measured on
#     a population the report does not rank. Dropping those rows and re-running the same
#     estimators gives CFD 0.025 (trapezoid, from 0.168) and CRISPR-Bulge 0.521
#     (tie-grouped full precision, from 0.622) -- measured on the benchmark, see
#     risk_card._ONTARGET_NOTE and METHODS section 3.
# Caption text only: no curve, AP, count or marked point moves (pr_figure_data.npz is
# untouched, and test_pr_figures re-derives every printed number from it).
OVERALL_SUPTITLE = (
    "Precision-recall on externally validated off-target sites "
    "(CRISPR-Bulge: 40 held-out guides; CFD: all 46, nothing held out)\n"
    "positives include each guide's own on-target (46 of 359, 40 of 287); "
    "excluding them: CFD 0.025, CRISPR-Bulge 0.521\n"
    "marked points = the report's own thresholds (CB sensitive default = 95%-recall)")


# ---------------------------------------------------------------- curve maths
def pr_curve(scores, labels):
    """Tie-grouped precision-recall curve and step average precision.

    One point per distinct score, precision taken at the END of each tie group,
    so the result does not depend on the order of equally scored rows (the old
    stratified-figure recipe sorted ties by benchmark row order, which made its
    average precision a property of the input file rather than of the data).
    ``ap`` equals ``sklearn.metrics.average_precision_score``.
    """
    s = np.asarray(scores, dtype=np.float64)
    y = np.asarray(labels).astype(np.int64)
    if s.size != y.size:
        raise ValueError("scores and labels differ in length")
    n_pos = int(y.sum())
    if n_pos == 0 or s.size == 0:
        raise ValueError("a precision-recall curve needs at least one positive")
    order = np.argsort(-s, kind="stable")
    s = s[order]
    y = y[order]
    group_end = np.append(s[1:] != s[:-1], True)
    tp = np.cumsum(y)[group_end].astype(np.float64)
    fp = np.cumsum(1 - y)[group_end].astype(np.float64)
    rec = tp / n_pos
    prec = tp / (tp + fp)
    return rec, prec, average_precision(rec, prec), n_pos, int(s.size)


def average_precision(rec, prec):
    """Step average precision of a precision-recall polyline.

    ``sum(delta_recall * precision)``. Points that do not advance recall
    contribute nothing, which is why the stored (decimated) polyline gives the
    same number as the full one.
    """
    rec = np.asarray(rec, dtype=np.float64)
    prec = np.asarray(prec, dtype=np.float64)
    return float(np.sum(np.diff(np.concatenate(([0.0], rec))) * prec))


def decimate(rec, prec):
    """Drop points strictly inside a vertical run; the drawn polyline is identical.

    Within a run of equal recall, precision falls monotonically (true positives
    fixed, false positives growing), so every interior point lies exactly on the
    segment joining the run's endpoints.
    """
    rec = np.asarray(rec, dtype=np.float64)
    prec = np.asarray(prec, dtype=np.float64)
    keep = np.ones(rec.size, dtype=bool)
    if rec.size > 2:
        interior = (rec[1:-1] == rec[:-2]) & (rec[1:-1] == rec[2:])
        keep[1:-1] = ~interior
    return rec[keep], prec[keep]


def point_at(scores, labels, threshold):
    """(precision, recall, n_tp, n_fp) of the ``score >= threshold`` list."""
    s = np.asarray(scores, dtype=np.float64)
    y = np.asarray(labels).astype(bool)
    sel = s >= threshold
    tp = int(np.count_nonzero(sel & y))
    fp = int(np.count_nonzero(sel & ~y))
    n_pos = int(np.count_nonzero(y))
    prec = tp / (tp + fp) if tp + fp else float("nan")
    return prec, (tp / n_pos if n_pos else float("nan")), tp, fp


def compact_sci(x):
    """1.9e-4 -> '1.9e-4' (what the figure has always printed)."""
    mant, exp = f"{x:.1e}".split("e")
    return f"{mant}e{int(exp)}"


# ------------------------------------------------------------------- stage 1
BENCH_COLS = ("guide", "offtarget_seq", "n_mismatch", "n_bulge", "label",
              "cb_heldout", "crrna_aligned")
SCORED_COLS = ("guide", "offtarget_seq", "n_bulge", "label", "cfd_score",
               "crispr_bulge_score")


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_inputs(benchmark, scored):
    """Read the two row-aligned TSVs, checking that they really are row-aligned."""
    out = {k: [] for k in ("guide", "crrna", "off", "label", "n_bulge",
                           "cb_heldout", "in_domain", "cfd", "cb3")}
    csv.field_size_limit(10 ** 7)
    with open(benchmark) as bf, open(scored) as sf:
        rb, rs = csv.reader(bf, delimiter="\t"), csv.reader(sf, delimiter="\t")
        hb = {k: i for i, k in enumerate(next(rb))}
        hs = {k: i for i, k in enumerate(next(rs))}
        for col in BENCH_COLS:
            if col not in hb:
                raise SystemExit(f"{benchmark}: missing column {col!r}")
        for col in SCORED_COLS:
            if col not in hs:
                raise SystemExit(f"{scored}: missing column {col!r}")
        for i, (b, s) in enumerate(zip(rb, rs)):
            if b[hb["guide"]] != s[hs["guide"]] or b[hb["offtarget_seq"]] != s[hs["offtarget_seq"]]:
                raise SystemExit(f"benchmark/scored rows diverge at data row {i}")
            out["guide"].append(b[hb["guide"]])
            out["crrna"].append(b[hb["crrna_aligned"]])
            out["off"].append(b[hb["offtarget_seq"]])
            out["label"].append(b[hb["label"]] == "1")
            out["n_bulge"].append(int(b[hb["n_bulge"]]))
            out["cb_heldout"].append(b[hb["cb_heldout"]] == "1")
            # `in_domain` is a historical column (it is `n_bulge <= 1`, constant 1 on
            # a <=1-gap benchmark). NOTHING filters on it -- filtering on it was the
            # bug. It is read, when present, only so the diagnostics can show that.
            if "in_domain" in hb:
                out["in_domain"].append(b[hb["in_domain"]] == "1")
            out["cfd"].append(float(s[hs["cfd_score"]]))
            out["cb3"].append(float(s[hs["crispr_bulge_score"]]))
        if next(rb, None) is not None or next(rs, None) is not None:
            raise SystemExit("benchmark and scored have different row counts")
    return out


def _score_crispr_bulge(crrna, off, chunk, device):
    import crispr_bulge_score as cbs
    out = np.empty(len(crrna), dtype=np.float64)
    for i in range(0, len(crrna), chunk):
        j = min(i + chunk, len(crrna))
        out[i:j] = np.asarray(cbs.CRISPR_BULGE_predict_list(crrna[i:j], off[i:j],
                                                            device=device), dtype=np.float64)
        sys.stderr.write(f"  CRISPR-Bulge {j:,}/{len(crrna):,}\n")
        sys.stderr.flush()
    return out


def _curve_record(prefix, scores, labels):
    rec, prec, ap, n_pos, n = pr_curve(scores, labels)
    drec, dprec = decimate(rec, prec)
    if abs(average_precision(drec, dprec) - ap) > 1e-12:
        raise AssertionError("decimation changed the average precision")
    return {f"{prefix}_rec": drec, f"{prefix}_prec": dprec,
            f"{prefix}_ap": np.float64(ap), f"{prefix}_pos": np.int64(n_pos),
            f"{prefix}_n": np.int64(n)}


def build(args):
    data = _read_inputs(args.benchmark, args.scored)
    n = len(data["label"])
    label = np.array(data["label"])
    n_bulge = np.array(data["n_bulge"], dtype=np.int8)
    heldout = np.array(data["cb_heldout"])
    cfd = np.array(data["cfd"], dtype=np.float64)
    cb3 = np.array(data["cb3"], dtype=np.float64)
    sys.stderr.write(f"rows {n:,}  positives {int(label.sum()):,}  "
                     f"held-out rows {int(heldout.sum()):,}  "
                     f"held-out positives {int(label[heldout].sum()):,}\n")

    cb = _score_crispr_bulge(data["crrna"], data["off"], args.chunk, args.device)
    if args.save_scores:
        np.save(args.save_scores, cb)
        sys.stderr.write(f"wrote {args.save_scores}\n")
    failed = cb < 0
    if failed.any():
        sys.stderr.write(f"  WARNING: {int(failed.sum()):,} rows failed to score; dropped\n")
    ok = ~failed

    store = {}
    # pr_curves.png -- CFD over every site, CRISPR-Bulge over every held-out site.
    store.update(_curve_record("cfd", cfd, label))
    hm = heldout & ok
    store.update(_curve_record("cb", cb[hm], label[hm]))

    marks = []
    for t in CFD_MARKS:
        p, r, tp, fp = point_at(cfd, label, t)
        marks.append(("cfd", float(t), f"≥{t:g}", p, r, tp, fp))
    sens = cbop.sensitive_threshold()
    for t, lab in ((sens, f"≥{compact_sci(sens)}\n({cbop.recall_target()}% recall)"),
                   (0.2, "≥0.2"), (0.5, "≥0.5")):
        p, r, tp, fp = point_at(cb[hm], label[hm], t)
        marks.append(("cb", float(t), lab, p, r, tp, fp))
    store["marks"] = np.array(json.dumps(marks))

    # pr_bulge_stratified.png -- same supports and scores, split by bulge class.
    for panel, scores, mask in (("scfd", cfd, np.ones(n, dtype=bool)), ("scb", cb, hm)):
        for key, nb in (("all", None), ("b0", 0), ("b1", 1)):
            sel = mask if nb is None else (mask & (n_bulge == nb))
            store.update(_curve_record(f"{panel}_{key}", scores[sel], label[sel]))
    store["scfd_guides"] = np.int64(len(set(data["guide"])))
    store["scb_guides"] = np.int64(len({g for g, h in zip(data["guide"], hm) if h}))

    store["meta"] = np.array(json.dumps({
        "generator": "scripts/render_pr_figures.py build",
        "benchmark_md5": _md5(args.benchmark), "scored_md5": _md5(args.scored),
        "benchmark_rows": n,
        "cb_support": "cb_heldout == 1 (CRISPR-Bulge training guides excluded)",
        "cfd_support": "all rows (CFD has no training set)",
        "cb_precision": "full precision, re-scored by PostProcess/crispr_bulge_score.py "
                        "(the 3-decimal column in scored.tsv floors at 5e-4, above the "
                        "report's own 1.9e-4 operating point, and ties most of the "
                        "benchmark at 0.000)",
        "ap_recipe": "tie-grouped step average precision (sklearn average_precision_score)",
        "sampling": "none -- every negative is used",
        "sensitive_threshold": sens, "recall_target": cbop.recall_target(),
    }))
    np.savez_compressed(args.out, **store)
    sys.stderr.write(f"wrote {args.out}\n")

    if args.diagnostics:
        _write_diagnostics(args.diagnostics, label, n_bulge, heldout, ok, cfd, cb, cb3,
                           data["in_domain"], data["guide"])


def _roworder_ap(scores, labels):
    """The old stratified-figure recipe: ties broken by input row order."""
    s = np.asarray(scores, dtype=np.float64)
    y = np.asarray(labels).astype(np.int64)[np.argsort(-s, kind="mergesort")]
    n_pos = int(y.sum())
    tp = np.cumsum(y)
    prec = tp / np.arange(1, y.size + 1)
    return float((prec * y).sum() / n_pos) if n_pos else float("nan")


def _write_diagnostics(path, label, n_bulge, heldout, ok, cfd, cb, cb3, in_domain, guide):
    """Every before/after number the figure change moves, for the review record."""
    indom = np.array(in_domain, dtype=bool)
    out = {
        "rows": int(label.size), "positives": int(label.sum()),
        "guides": len(set(guide)),
        # the defect: a filter on in_domain drops nothing, so "(in-domain)" was not
        # a held-out panel. Absent entirely once the benchmark is rebuilt in-tree.
        "in_domain_present": bool(indom.size),
        "in_domain_true": int(indom.sum()),
        "in_domain_equals_n_bulge_le_1": bool(indom.size and np.array_equal(indom, n_bulge <= 1)),
        "in_domain_drops": int((~indom).sum()),
        "heldout_rows": int(heldout.sum()), "heldout_positives": int(label[heldout].sum()),
        "heldout_guides": len({g for g, h in zip(guide, heldout) if h}),
        "training_rows": int((~heldout).sum()),
        "training_positives": int(label[~heldout].sum()),
        "training_guides": len({g for g, h in zip(guide, heldout) if not h}),
        "strata": {},
    }
    hm = heldout & ok
    for name, scores, mask in (("cb_full_heldout", cb, hm),
                               ("cb_full_all", cb, ok),
                               ("cb_full_training", cb, (~heldout) & ok),
                               ("cb3_heldout", cb3, heldout),
                               ("cb3_all", cb3, np.ones(label.size, dtype=bool)),
                               ("cb3_training", cb3, ~heldout),
                               ("cfd_all", cfd, np.ones(label.size, dtype=bool)),
                               ("cfd_heldout", cfd, heldout)):
        out["strata"][name] = {}
        for key, nb in (("all", None), ("b0", 0), ("b1", 1)):
            sel = mask if nb is None else (mask & (n_bulge == nb))
            if not label[sel].any():
                continue
            _, _, ap, n_pos, nn = pr_curve(scores[sel], label[sel])
            out["strata"][name][key] = {
                "ap_step": round(ap, 6),
                "ap_roworder_legacy": round(_roworder_ap(scores[sel], label[sel]), 6),
                "positives": n_pos, "sites": nn}
    out["marked_points_exact"] = {}
    for t in CFD_MARKS:
        p, r, tp, fp = point_at(cfd, label, t)
        out["marked_points_exact"][f"cfd>={t:g}"] = {"precision": round(p, 6),
                                                     "recall": round(r, 6), "tp": tp, "fp": fp}
    for t in (cbop.sensitive_threshold(), 0.2, 0.5):
        p, r, tp, fp = point_at(cb[hm], label[hm], t)
        out["marked_points_exact"][f"cb>={t:g}"] = {"precision": round(p, 6),
                                                    "recall": round(r, 6), "tp": tp, "fp": fp}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    sys.stderr.write(f"wrote {path}\n")


# ------------------------------------------------------------------- stage 2
def overall_titles(d):
    """Panel titles of pr_curves.png, in drawing order."""
    return (f"CFD (Doench 2016)  (AUPR={float(d['cfd_ap']):.2f})",
            "CRISPR-Bulge (held-out, full precision)  "
            f"(AUPR={float(d['cb_ap']):.2f})")


def stratified_titles(d):
    """Panel titles of pr_bulge_stratified.png, in drawing order.

    Each names its own support: the CRISPR-Bulge panel is restricted to guides the
    model was not trained on, the CFD panel is not (CFD has no training set). The
    panel this replaced said "(in-domain)", which described no filter at all.
    """
    return (f"CFD  —  PR by bulge class (all {int(d['scfd_guides'])} guides)",
            "CRISPR-Bulge  —  PR by bulge class "
            f"({int(d['scb_guides'])} held-out guides, full precision)")


def stratified_legend(d, panel):
    """Legend entries of one stratified panel, in drawing order."""
    out = []
    for (lbl, nb, _colour, _style, _lw) in STRATA:
        p = f"{panel}_{ {None: 'all', 0: 'b0', 1: 'b1'}[nb] }"
        out.append(f"{lbl}: AP={float(d[f'{p}_ap']):.3f} "
                   f"(pos {int(d[f'{p}_pos'])}/{int(d[f'{p}_n'])})")
    return out


def _axes(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=.25)


def render_overall(d, path, dpi=130):
    """pr_curves.png -- one panel per scorer, with the report's operating points."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    marks = json.loads(str(d["marks"]))
    titles = overall_titles(d)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(d["cfd_rec"], d["cfd_prec"], color=CFD_COLOR, lw=2)
    ax[0].set_title(titles[0], fontsize=10)
    ax[1].plot(d["cb_rec"], d["cb_prec"], color=CB_COLOR, lw=2)
    ax[1].set_title(titles[1], fontsize=10)
    for panel, _t, lab, prec, rec, _tp, _fp in marks:
        a, colour = (ax[0], CFD_COLOR) if panel == "cfd" else (ax[1], CB_COLOR)
        a.plot(rec, prec, "o", color=colour, ms=6)
        a.annotate(lab, (rec, prec), textcoords="offset points", xytext=(6, 4),
                   fontsize=8, fontweight="bold")
    for a in ax:
        a.set_xlabel("Recall")
        a.set_ylabel("Precision")
        _axes(a)
    fig.suptitle(OVERALL_SUPTITLE, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def render_stratified(d, path, dpi=145):
    """pr_bulge_stratified.png -- precision-recall split by bulge class."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.7))
    for a, panel, title in zip(ax, ("scfd", "scb"), stratified_titles(d)):
        for (_lbl, nb, colour, style, lw), label in zip(STRATA, stratified_legend(d, panel)):
            p = f"{panel}_{ {None: 'all', 0: 'b0', 1: 'b1'}[nb] }"
            a.plot(d[f"{p}_rec"], d[f"{p}_prec"], style, color=colour, lw=lw, label=label)
        a.set_title(title, fontsize=11)
        a.set_xlabel("recall")
        a.set_ylabel("precision")
        a.legend(fontsize=8, loc="upper right")
        _axes(a)
    fig.suptitle("Off-target score precision/recall, stratified by bulge count\n"
                 "(benchmark: Tsai2015 + Refined-TrueOT + Chen2017 + Listgarten2018)",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


MANIFEST = "pr_figures.sha256"


def toolchain():
    """'matplotlib X.Y.Z freetype A.B.C' -- PNG bytes depend on both."""
    import matplotlib
    from matplotlib import ft2font
    return f"matplotlib {matplotlib.__version__} freetype {ft2font.__freetype_version__}"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(outdir, data_path):
    """Chain the PNGs to the curve data they were drawn from.

    The manifest covers the data file as well as the PNGs, so changing the curve
    data without re-rendering is as loud a test failure as editing a PNG.
    """
    lines = [f"# rendered by scripts/render_pr_figures.py with {toolchain()}\n"]
    for p in (data_path, os.path.join(outdir, "pr_curves.png"),
              os.path.join(outdir, "pr_bulge_stratified.png")):
        lines.append(f"{_sha256(p)}  {os.path.basename(p)}\n")
    path = os.path.join(outdir, MANIFEST)
    with open(path, "w") as fh:
        fh.writelines(lines)
    return path


def render(args):
    d = np.load(args.data, allow_pickle=False)
    os.makedirs(args.outdir, exist_ok=True)
    overall = os.path.join(args.outdir, "pr_curves.png")
    strat = os.path.join(args.outdir, "pr_bulge_stratified.png")
    render_overall(d, overall)
    render_stratified(d, strat)
    manifest = write_manifest(args.outdir, args.data)
    for p in (overall, strat):
        print(f"{_md5(p)}  {p}")
    with open(manifest) as fh:
        sys.stdout.write(fh.read())


# ----------------------------------------------------------------------- cli
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="score the benchmark and write the curve data")
    b.add_argument("--benchmark", required=True)
    b.add_argument("--scored", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--diagnostics")
    b.add_argument("--save-scores", help="also dump the raw CRISPR-Bulge scores (.npy)")
    b.add_argument("--chunk", type=int, default=50000)
    b.add_argument("--device", default="cpu")
    b.set_defaults(func=build)
    r = sub.add_parser("render", help="draw both PNGs from the curve data")
    r.add_argument("--data", required=True)
    r.add_argument("--outdir", required=True)
    r.set_defaults(func=render)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
