#!/usr/bin/env python
"""Self-contained, shareable CRISPRme off-target report (v2 -- IND briefing book).

Given a CRISPRme result folder (or a bare ``*integrated_results.tsv``), this
module produces a single, easily-transferable ZIP::

    <jobid>_report.zip
      report.html                 # self-contained: base64 PNG plots, inline
                                   #   top-1000 table, inline CSS, opens offline
      integrated_results.tsv.gz   # the full results, the HTML links to it
                                   #   with a RELATIVE href (resolves post-unzip)
      top1000.tsv                 # the top-1000-by-CFD rows shown in the table,
                                   #   its own RELATIVE download link

The report is a portable *digest* of the full interactive CRISPRme website
result (personal risk cards, etc. stay in the website). It is meant for a
collaborator preparing an IND briefing book -- e.g. designing a targeted-NGS
rhAMP-Seq confirmation panel from the predicted off-targets (Saha lab, TRAC
guide, SpCas9 NRG).

Report structure (top -> bottom)
--------------------------------
1. HEADER + GLOBAL SUMMARY: mirrors the web result-page top table. Left card:
   gRNA (spacer+PAM), nuclease, Aggregated Specificity Score (0-100, from
   ``.<jobid>.acfd_CFD.txt`` if present, else "CFD score not available"). Right:
   the "Off-targets by Mismatch (MM) and Bulge (B)" matrix, grouped REFERENCE vs
   VARIANT, then by bulge count, columns Total + 0MM..<mm>MM.
2. KEY GRAPHICAL REPORT: the CRISPRme-paper CFD scatter (site index log-x vs
   score, red REF / blue ALT points sized by allele frequency, ref->alt arrows,
   top site rsID-annotated), in up to three panels: by CFD, by DELTA=ALT-REF, and
   by CRISTA (only when CRISTA was computed).
3. SIMPLIFIED reference-vs-population plot (v1 single view).
4. RECOMMENDED VALIDATION PANEL: candidate counts at CFD / edit-distance
   thresholds, a suggested tiered panel (union of top-100-by-CFD + all mm+b<=2 +
   all variant-created with CFD>=0.05), and an IND rationale.
5. DOWNLOADS: full integrated_results.tsv.gz + top1000.tsv (both bundled).
6. SCROLLABLE TOP-1000 TABLE (by CFD desc, mm+b<=1 excluded), with a PAM-creation
   column and CRISTA when computed.
7. FOOTER: CRISPRme version + provenance stamp + fixed research-only disclaimer.

Design goals / robustness posture
---------------------------------
* Pure stdlib + matplotlib + pandas + numpy (the pipeline conda env). No Dash,
  no Jinja, no network. The HTML has no <script> and no external <link>, so it
  opens with ``file://`` on any machine.
* Columns are selected BY NAME from the header (never fixed indices), because
  the dict-less 85-col schema and the dict-based schema differ in column set /
  order. Missing columns (CRISTA, per-dataset, annotation on dict-less) are
  optional -- the report degrades, it never crashes.
* ONE canonical partition everywhere: variant-created := Not_found_in_REF=="y";
  reference := the rest; on-target := mm+b==0 (reported separately).
  variant + reference + on-target == total, exactly.
* Every plot/table is wrapped in try/except; a failure inlines a small
  placeholder rather than failing the report (issue-#143 "preserve results").

Runnable both as ``python -m PostProcess.generate_report --result-dir ...`` and
via ``crisprme.py generate-report --result-dir ...``.
"""

from __future__ import annotations

import argparse
import base64
import glob
import gzip
import html
import io
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime

import matplotlib

# never touch an X server (pipeline convention, CRISPRme_plots.py:21)
matplotlib.use("Agg")

import matplotlib.lines as mlines  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# report-generator version -- bumped here, stamped in the footer provenance line.
REPORT_GENERATOR_VERSION = "2.0"

# --------------------------------------------------------------------------- #
# Recommended-validation-panel thresholds (module-level constants, section 4)
# --------------------------------------------------------------------------- #
CFD_THRESHOLDS = (0.5, 0.2, 0.1, 0.05)
MMB_THRESHOLDS = (1, 2, 3, 4)
# tiered-panel construction knobs
PANEL_TOP_N_BY_CFD = 100
PANEL_MMB_MAX = 2
PANEL_VARIANT_CFD_MIN = 0.05

# --------------------------------------------------------------------------- #
# Column-name resolution helpers
# --------------------------------------------------------------------------- #
# The "highest_CFD" projection is the one the report uses (CFD is the primary
# ranking the IND reviewer cares about). We resolve every column by its base
# name against whichever suffixes a given schema uses.
_PROJ = "(highest_CFD)"
_CRISTA_PROJ = "(highest_CRISTA)"

# base-name -> preferred header suffix list (first match wins)
_COLS = {
    "guide": ["Spacer+PAM"],
    "chrom": ["Chromosome"],
    "pos": [f"Start_coordinate_{_PROJ}", "Start_coordinate"],
    "strand": [f"Strand_{_PROJ}", "Strand"],
    "aln_ref": [f"Aligned_protospacer+PAM_REF_{_PROJ}", "Aligned_protospacer+PAM_REF"],
    "aln_alt": [f"Aligned_protospacer+PAM_ALT_{_PROJ}", "Aligned_protospacer+PAM_ALT"],
    "pam": [f"PAM_{_PROJ}", "PAM"],
    "mm": [f"Mismatches_{_PROJ}", "Mismatches"],
    "bulges": [f"Bulges_{_PROJ}", "Bulges"],
    "mmb": [f"Mismatches+bulges_{_PROJ}", "Mismatches+bulges"],
    "origin": [f"REF/ALT_origin_{_PROJ}", "REF/ALT_origin"],
    "pam_creation": [f"PAM_creation_{_PROJ}", "PAM_creation"],
    "cfd": [f"CFD_score_{_PROJ}", "CFD_score"],
    "cfd_ref": [f"CFD_score_REF_{_PROJ}"],
    "cfd_alt": [f"CFD_score_ALT_{_PROJ}"],
    "var_genome": [f"Variant_info_genome_{_PROJ}", "Variant_info_genome"],
    "maf": [f"Variant_MAF_{_PROJ}", "Variant_MAF"],
    "rsid": [f"Variant_rsID_{_PROJ}", "Variant_rsID"],
    "samples": [f"Variant_samples_{_PROJ}", "Variant_samples"],
    "not_in_ref": ["Not_found_in_REF"],
    "gene_name": ["Annotation_closest_gene_name"],
    "gene_dist": ["Annotation_closest_gene_distance_(kb)"],
    # CRISTA projection (present only when CRISTA was computed this run)
    "crista": [f"CRISTA_score_{_CRISTA_PROJ}", "CRISTA_score"],
    "crista_ref": [f"CRISTA_score_REF_{_CRISTA_PROJ}"],
    "crista_alt": [f"CRISTA_score_ALT_{_CRISTA_PROJ}"],
    "crista_mmb": [f"Mismatches+bulges_{_CRISTA_PROJ}"],
    "crista_maf": [f"Variant_MAF_{_CRISTA_PROJ}"],
    "crista_samples": [f"Variant_samples_{_CRISTA_PROJ}"],
    "crista_rsid": [f"Variant_rsID_{_CRISTA_PROJ}"],
}

_NA_TOKENS = {"", "na", "n", ".", "nan", "none", "-1"}


def _resolve(header, keys):
    """Map logical column keys -> actual header names present in the TSV."""
    present = set(header)
    resolved = {}
    for key in keys:
        for candidate in _COLS[key]:
            if candidate in present:
                resolved[key] = candidate
                break
    return resolved


def _is_na(value) -> bool:
    if value is None:
        return True
    try:
        if isinstance(value, float) and np.isnan(value):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() in _NA_TOKENS


def _min_maf(raw):
    """Min AF across a comma-joined multi-SNP haplotype (CRISPRme_plots.py:58-60)."""
    if _is_na(raw):
        return None
    vals = []
    for tok in str(raw).split(","):
        tok = tok.strip()
        if _is_na(tok):
            continue
        try:
            vals.append(float(tok))
        except ValueError:
            continue
    return min(vals) if vals else None


def _first_non_na(raw):
    """First non-NA token from a comma-joined field (e.g. rsID list)."""
    if _is_na(raw):
        return None
    for tok in str(raw).split(","):
        tok = tok.strip()
        if not _is_na(tok):
            return tok
    return None


def _to_int_series(series):
    return pd.to_numeric(series, errors="coerce").fillna(-1).astype(int)


def _to_float_series(series):
    return pd.to_numeric(series, errors="coerce")


def crista_computed(df, cols):
    """True when the CRISTA_score projection column exists AND has >=1 real value.

    The stable schema always ships the CRISTA columns; a dict-less / CRISTA-off
    run either omits them or leaves them blank. We include the CRISTA scatter
    only when a real (numeric, non-NA) CRISTA score is present in this run.
    """
    if "crista" not in cols or cols["crista"] not in df.columns:
        return False
    vals = _to_float_series(df[cols["crista"]])
    return bool(vals.notna().any())


# --------------------------------------------------------------------------- #
# Canonical partition (variant / reference / on-target) -- single source
# --------------------------------------------------------------------------- #
def partition_masks(df, cols):
    """Return (variant_mask, reference_mask, ontarget_mask) over df.

    variant-created := Not_found_in_REF == "y" (fallback: origin == alt)
    on-target       := mm+b == 0
    reference       := everything else

    By construction the three masks are disjoint and cover every row, so
    variant.sum() + reference.sum() + ontarget.sum() == len(df) EXACTLY.
    """
    n = len(df)
    if "mmb" in cols:
        ontarget = _to_int_series(df[cols["mmb"]]) == 0
    else:
        ontarget = pd.Series([False] * n, index=df.index)

    if "not_in_ref" in cols:
        variant_raw = (
            df[cols["not_in_ref"]].astype(str).str.strip().str.lower().eq("y")
        )
    elif "origin" in cols:
        variant_raw = df[cols["origin"]].astype(str).str.strip().str.lower().eq("alt")
    else:
        variant_raw = pd.Series([False] * n, index=df.index)

    # on-target takes precedence, so the three are mutually exclusive
    variant = variant_raw & (~ontarget)
    reference = (~variant) & (~ontarget)
    return variant, reference, ontarget


# --------------------------------------------------------------------------- #
# Aggregated Specificity Score (.<jobid>.acfd_CFD.txt)
# --------------------------------------------------------------------------- #
def read_specificity_score(result_dir, job_id, guides):
    """Read the aggregated CFD specificity score (0-100) from the acfd sidecar.

    File format (process_summaries.py:211-215): ``<guide>\\t<score>\\tNA\\tNA``
    where <score> = 100/(100+sum_cfds) in [0,1). The web converts it to 0-100 via
    x*100 (results_page.py:2881-2888) and shows "CFD score not available" when
    x>=1 (i.e. no scoreable off-targets). Returns a display string.
    """
    if not result_dir:
        return "CFD score not available"
    acfd_file = os.path.join(result_dir, f".{job_id}.acfd_CFD.txt")
    if not os.path.isfile(acfd_file):
        return "CFD score not available"
    try:
        with open(acfd_file) as handle:
            lines = [ln for ln in handle.read().strip().split("\n") if ln.strip()]
    except OSError:
        return "CFD score not available"
    # map guide -> score; if a single guide, take the only/first row
    by_guide = {}
    for ln in lines:
        parts = ln.split("\t")
        if len(parts) >= 2:
            try:
                by_guide[parts[0].strip()] = float(parts[1])
            except ValueError:
                continue
    if not by_guide:
        return "CFD score not available"

    def _fmt(x):
        if 0 <= x < 1:
            return f"{x * 100:.3f}"
        return "CFD score not available"

    if len(by_guide) == 1:
        return _fmt(next(iter(by_guide.values())))
    # multi-guide: label per guide (in the order they appear in the table)
    ordered = guides if guides else list(by_guide.keys())
    parts = []
    for g in ordered:
        if g in by_guide:
            parts.append(f"{g}: {_fmt(by_guide[g])}")
    return "; ".join(parts) if parts else "CFD score not available"


# --------------------------------------------------------------------------- #
# Summary metadata: .Params.txt / .version.txt / filename fallback
# --------------------------------------------------------------------------- #
def _read_kv_sidecar(path):
    """Read a tab-separated ``key\\tvalue`` sidecar (Params.txt/.version.txt)."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    try:
        with open(path) as handle:
            for line in handle:
                line = line.rstrip("\n")
                if "\t" not in line:
                    continue
                key, _, val = line.partition("\t")
                out[key.strip()] = val.strip()
    except OSError:
        pass
    return out


def _parse_results_filename(tsv_path):
    """Fallback summary from ``<guide>+<pam>_<ref>+<vcf>_<mm>+<bmax>_integrated_results.tsv``.

    e.g. CTCTCAGCTGGTACACGGCANNN+NRG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv
      -> guide=CTCTCAGCTGGTACACGGCANNN, pam=NRG, genome=hg38,
         datasets=1000G+HGDP, mm=6, bMax=3
    """
    info = {}
    base = os.path.basename(tsv_path)
    for suffix in (
        "_integrated_results.tsv.gz",
        "_integrated_results.tsv",
        ".bestMerge.txt.integrated_results.tsv",
    ):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    # <guide>+<pam>_<ref>+<vcf>_<mm>+<bmax>
    try:
        guide_pam, rest = base.split("_", 1)
        if "+" in guide_pam:
            info["guide"], info["pam"] = guide_pam.split("+", 1)
        else:
            info["guide"] = guide_pam
        # rest = <ref>+<vcf>_<mm>+<bmax>  (vcf may itself contain '_')
        if "_" in rest:
            genome_vcf, mm_b = rest.rsplit("_", 1)
        else:
            genome_vcf, mm_b = rest, ""
        if "+" in genome_vcf:
            info["genome"], vcf = genome_vcf.split("+", 1)
            # vcf like "hg38_1000G_HGDP" -> datasets after the ref token
            vcf_tokens = vcf.split("_")
            if vcf_tokens and vcf_tokens[0] == info["genome"]:
                vcf_tokens = vcf_tokens[1:]
            info["datasets"] = "+".join(t for t in vcf_tokens if t) or vcf
        else:
            info["genome"] = genome_vcf
        if "+" in mm_b:
            info["mm"], info["bmax"] = mm_b.split("+", 1)
    except ValueError:
        pass
    return info


def build_summary_meta(result_dir, tsv_path, df, cols, params_override=None):
    """Collect summary metadata used by the header card and matrix.

    Returns a dict with resolved guide list, nuclease, pam, genome, datasets, mm,
    bdna, brna, bmax, max_edits, date, version, and the canonical partition
    counts (n_total, n_variant, n_reference, n_ontarget, n_offtarget).
    """
    params = {}
    version = None
    date = None
    if result_dir:
        params = _read_kv_sidecar(os.path.join(result_dir, ".Params.txt"))
        if not params:
            params = _read_kv_sidecar(os.path.join(result_dir, "Params.txt"))
        ver = _read_kv_sidecar(os.path.join(result_dir, ".version.txt"))
        version = ver.get("crisprme_version")
        for probe in (".Params.txt", "Params.txt", os.path.basename(tsv_path)):
            candidate = os.path.join(result_dir, probe)
            if os.path.isfile(candidate):
                date = datetime.fromtimestamp(
                    os.path.getmtime(candidate)
                ).strftime("%Y-%m-%d")
                break
    if params_override:
        params = {**params, **params_override}

    fn = _parse_results_filename(tsv_path)

    guides = []
    if "guide" in cols and cols["guide"] in df.columns:
        guides = [g for g in df[cols["guide"]].dropna().unique().tolist()]
    if not guides and fn.get("guide"):
        guides = [fn["guide"]]

    pam = params.get("Pam") or fn.get("pam") or "n/a"
    nuclease = params.get("Nuclease") or "n/a"
    genome = (
        params.get("Genome_selected")
        or params.get("Genome_ref")
        or fn.get("genome")
        or "n/a"
    )
    datasets = fn.get("datasets")
    if not datasets:
        idx = params.get("Genome_idx", "")
        ref_comp = params.get("Ref_comp", "")
        if "True" in str(ref_comp):
            datasets = "reference + variants"
        elif idx and idx != "None":
            datasets = idx
        else:
            datasets = "reference only"

    def _int_or_none(x):
        try:
            return int(str(x).strip())
        except (TypeError, ValueError):
            return None

    mm = params.get("Mismatches") or fn.get("mm") or "n/a"
    mm_int = _int_or_none(mm)
    bdna = params.get("DNA")
    brna = params.get("RNA")
    bmax = params.get("Max_bulges") or fn.get("bmax") or "n/a"
    bmax_int = _int_or_none(bmax)
    max_edits = params.get("Max_total_edits") or "n/a"

    variant, reference, ontarget = partition_masks(df, cols)
    n_total = len(df)
    n_variant = int(variant.sum())
    n_reference = int(reference.sum())
    n_ontarget = int(ontarget.sum())
    n_offtarget = n_total - n_ontarget

    return {
        "guides": guides,
        "guide_display": ", ".join(guides) if guides else "n/a",
        "nuclease": nuclease,
        "pam": pam,
        "genome": genome,
        "datasets": datasets,
        "mm": str(mm),
        "mm_int": mm_int,
        "bdna": bdna,
        "brna": brna,
        "bmax": str(bmax),
        "bmax_int": bmax_int,
        "max_edits": str(max_edits),
        "date": date or "n/a",
        "version": version,
        "n_total": n_total,
        "n_variant": n_variant,
        "n_reference": n_reference,
        "n_ontarget": n_ontarget,
        "n_offtarget": n_offtarget,
    }


# --------------------------------------------------------------------------- #
# SECTION 1: global off-target-by-MM-and-B matrix (REFERENCE vs VARIANT)
# --------------------------------------------------------------------------- #
def build_mmb_matrix(df, cols, meta):
    """Off-targets binned by origin(REFERENCE/VARIANT) x bulge x mismatch.

    Mirrors the web result-page top matrix: columns Total, 0MM..<mm>MM; rows
    grouped REFERENCE then VARIANT, within each a row per bulge count 0..maxbulges.
    Uses the canonical partition (variant := Not_found_in_REF=="y"; on-target
    rows mm+b==0 are excluded from the off-target matrix). Returns a dict with
    ``mm_cols`` (0..mm) and ``groups`` = list of (label, [ (bulge, total, [per-mm
    counts]) ... ]).
    """
    variant, reference, _ontarget = partition_masks(df, cols)

    mm_series = _to_int_series(df[cols["mm"]]) if "mm" in cols else None
    b_series = _to_int_series(df[cols["bulges"]]) if "bulges" in cols else None
    if mm_series is None or b_series is None:
        return None

    # column extent: prefer the run's configured mm; else observed max
    max_mm = meta.get("mm_int")
    if max_mm is None:
        max_mm = int(mm_series[mm_series >= 0].max()) if (mm_series >= 0).any() else 0
    # row extent: prefer configured max bulges; else observed max
    max_b = meta.get("bmax_int")
    if max_b is None:
        max_b = int(b_series[b_series >= 0].max()) if (b_series >= 0).any() else 0

    mm_cols = list(range(0, max_mm + 1))

    def _rows_for(mask):
        rows = []
        for b in range(0, max_b + 1):
            b_mask = mask & (b_series == b)
            per_mm = [int((b_mask & (mm_series == m)).sum()) for m in mm_cols]
            total = int(b_mask.sum())
            rows.append((b, total, per_mm))
        return rows

    groups = [
        ("REFERENCE", _rows_for(reference)),
        ("VARIANT", _rows_for(variant)),
    ]
    return {"mm_cols": mm_cols, "groups": groups}


def render_summary_and_matrix(meta, spec_score, matrix):
    """Section 1 HTML: header card (left) + MM/B matrix (right)."""
    bdna = meta["bdna"]
    brna = meta["brna"]
    bulge_disp = (
        f"{bdna if bdna is not None else 'n/a'} / "
        f"{brna if brna is not None else 'n/a'}"
        + (f" (max {meta['bmax']})" if meta["bmax"] != "n/a" else "")
    )

    left_rows = [
        ("gRNA (spacer+PAM)", meta["guide_display"]),
        ("Nuclease", meta["nuclease"]),
        ("PAM", meta["pam"]),
        ("Genome", meta["genome"]),
        ("Variant dataset(s)", meta["datasets"]),
        ("Mismatches", meta["mm"]),
        ("Bulges (DNA / RNA)", bulge_disp),
        ("Max total edits", meta["max_edits"]),
        ("Aggregated Specificity Score (0-100)", spec_score),
    ]
    left_html = "".join(
        f'<tr><td class="k">{_esc(k)}</td><td>{_esc(v)}</td></tr>'
        for k, v in left_rows
    )

    # right: the MM/B matrix
    if matrix is None:
        matrix_html = "<p>Off-target matrix unavailable (missing MM/bulge columns).</p>"
    else:
        mm_cols = matrix["mm_cols"]
        head_cells = ['<th class="grp">Origin</th>', "<th>Bulges (B)</th>", "<th>Total</th>"]
        head_cells += [f"<th>{m}MM</th>" for m in mm_cols]
        head_html = "".join(head_cells)

        body = []
        for label, rows in matrix["groups"]:
            grp_total = sum(r[1] for r in rows)
            first = True
            for b, total, per_mm in rows:
                cells = []
                if first:
                    cells.append(
                        f'<td class="grp" rowspan="{len(rows)}">{_esc(label)}'
                        f'<br><span class="grp-total">({grp_total:,})</span></td>'
                    )
                    first = False
                cells.append(f"<td>{b}B</td>")
                cells.append(f"<td class='tot'>{total:,}</td>")
                cells += [f"<td>{c:,}</td>" for c in per_mm]
                body.append("<tr>" + "".join(cells) + "</tr>")
        matrix_html = (
            '<div class="matrix-wrap"><table class="matrix">'
            f"<thead><tr>{head_html}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>"
        )

    return f"""
<div class="summary-grid">
  <div class="summary-card">
    <table class="summary-table"><tbody>{left_html}</tbody></table>
  </div>
  <div class="matrix-card">
    <div class="matrix-title">Off-targets by Mismatch (MM) and Bulge (B)</div>
    {matrix_html}
    <p class="caption">Off-target counts (on-target row mm+b=0 excluded), grouped
    by REFERENCE vs VARIANT origin (variant-created := Not_found_in_REF), then by
    bulge count. Total column = row sum across mismatches.</p>
  </div>
</div>
"""


# --------------------------------------------------------------------------- #
# samplesID -> superpopulation mapping (for the simplified population plot)
# --------------------------------------------------------------------------- #
def load_sample_superpop(samplesid_dir, datasets_hint=""):
    """sample_id -> SUPERPOPULATION_ID, loaded from any samplesID files present.

    Mirrors process_summaries.py:38-63 (header
    ``#SAMPLE_ID\\tPOPULATION_ID\\tSUPERPOPULATION_ID\\tSEX``). Returns an empty
    dict when no samplesID files are resolvable; the caller then falls back to
    per-dataset (1000G vs HGDP) provenance.
    """
    mapping = {}
    if not samplesid_dir or not os.path.isdir(samplesid_dir):
        return mapping
    for path in sorted(glob.glob(os.path.join(samplesid_dir, "samplesID*.txt"))):
        try:
            with open(path) as handle:
                handle.readline()  # header
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3 and parts[0]:
                        mapping[parts[0]] = parts[2]
        except OSError:
            continue
    return mapping


def _dataset_of(sample_id):
    """Infer dataset provenance from an ID prefix (preserve-provenance rule)."""
    sid = sample_id.strip()
    if sid.startswith("HGDP"):
        return "HGDP"
    if sid.startswith("HG") or sid.startswith("NA"):
        return "1000G"
    return "other"


# --------------------------------------------------------------------------- #
# Plots (matplotlib -> in-memory PNG -> base64 data URI)
# --------------------------------------------------------------------------- #
def _fig_to_data_uri(fig, dpi=120):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode()
    return f"data:image/png;base64,{encoded}"


def _placeholder_uri(message):
    """A tiny figure that says a plot could not be drawn (never fails the run)."""
    try:
        fig, ax = plt.subplots(figsize=(6, 1.6))
        ax.axis("off")
        ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11, wrap=True)
        return _fig_to_data_uri(fig, dpi=90)
    except Exception:  # pragma: no cover - graphics totally unavailable
        # 1x1 transparent PNG
        return (
            "data:image/png;base64,"
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk"
            "+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=="
        )


def _cfd_style_scatter(sub, cols, score_key, ref_key, alt_key, xlabel, title):
    """The CRISPRme-paper CFD/CRISTA scatter, adapted from
    CRISPRme_plots.py:plot_with_CFD_score (~L347).

    ``sub`` is an ALREADY-ORDERED (top-first) top-N frame. Red REF points and
    blue ALT points, both sized by allele frequency; gray arrows connect the
    ref->alt pair for the same site; the top-ranked off-target is annotated with
    its rsID. Returns a base64 data URI.
    """
    work = sub.reset_index(drop=True).copy()
    n = len(work)
    if n == 0:
        raise ValueError("no rows to plot")
    work["_index"] = np.arange(1, n + 1)

    # min-AF across a multi-SNP haplotype; -1 => AF unknown (CRISPRme_plots.py:54-60)
    maf_col = cols.get("_maf_for_" + score_key)
    if maf_col is None or maf_col not in work.columns:
        maf_col = cols.get("maf") if score_key == "cfd" else cols.get("crista_maf")
    if maf_col and maf_col in work.columns:
        af = work[maf_col].astype(str).str.split(",").apply(
            lambda toks: min(
                [float(t) for t in toks if _num_ok(t)] or [-1.0]
            )
        )
    else:
        af = pd.Series([-1.0] * n)
    work["AF"] = pd.to_numeric(af, errors="coerce").fillna(-1.0)

    samp_col = (
        cols.get("samples") if score_key == "cfd" else cols.get("crista_samples")
    )
    if samp_col and samp_col in work.columns:
        _samp = work[samp_col].astype(str)
        has_var = _samp.str.len().gt(1) & ~_samp.str.lower().isin(
            ["nan", "na", "n", "."]
        )
    else:
        has_var = pd.Series([False] * n)
    work["_has_variant"] = has_var.values

    # point sizes (CRISPRme_plots.py:70-80)
    work["plot_AF"] = np.where(
        work["AF"] >= 0,
        np.sqrt(np.clip(work["AF"], 0, None) * 1000 + 0.001 * 1000) + 6.0,
        np.where(work["_has_variant"], 12.0, np.nan),
    )
    work["ref_AF"] = np.sqrt(np.clip(1 - work["AF"], 0, None) * 1000)

    y_ref = _to_float_series(work[ref_key]) if ref_key in work.columns else None
    y_alt = _to_float_series(work[alt_key]) if alt_key in work.columns else None
    # fall back to the combined score column if REF/ALT variants are missing
    y_combined = _to_float_series(work[score_key]) if score_key in work.columns else None
    if y_ref is None:
        y_ref = y_combined
    if y_alt is None:
        y_alt = y_combined

    transparent_red = "#e5323280"
    transparent_blue = "#2b6cb080"

    plt.rcParams.update({"font.size": 8})
    fig, ax = plt.subplots(figsize=(8.0, 3.2))

    ax.scatter(
        work["_index"], y_ref, s=work["ref_AF"], c=transparent_red, zorder=1,
        label="_ref",
    )
    ax.scatter(
        work["_index"], y_alt, s=work["plot_AF"], c=transparent_blue, zorder=2,
        edgecolors="black", linewidths=0.4, label="_alt",
    )
    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(score_label(score_key) + " score")
    ax.set_title(title)
    ax.set_xlim(0.9, max(1000, n))
    ax.set_ylim(0, 1)

    # ref -> alt arrows (CRISPRme_plots.py:160-177)
    for x, yr, ya in zip(work["_index"], y_ref, y_alt):
        if pd.isna(yr) or pd.isna(ya):
            continue
        z = ya - yr
        if abs(z) < 1e-9:
            continue
        try:
            ax.arrow(
                x, yr + 0.02, 0, z - 0.04,
                color="gray",
                head_width=(x * (10 ** 0.005 - 10 ** (-0.005))),
                head_length=0.02, length_includes_head=True, zorder=0, alpha=0.5,
            )
        except (ValueError, FloatingPointError):
            continue

    # annotate top-ranked off-target with its rsID
    rsid_col = cols.get("rsid") if score_key == "cfd" else cols.get("crista_rsid")
    if rsid_col and rsid_col in work.columns:
        top_rsid = _first_non_na(work.iloc[0][rsid_col])
        top_y = y_alt.iloc[0] if pd.notna(y_alt.iloc[0]) else y_ref.iloc[0]
        if top_rsid and pd.notna(top_y):
            ax.annotate(
                f"top site: {top_rsid}",
                xy=(1, top_y),
                xytext=(2.2, min(0.96, top_y + 0.12)),
                fontsize=8,
                arrowprops=dict(arrowstyle="->", color="#333", lw=0.7),
            )

    # allele-frequency size legend (CRISPRme_plots.py:88-114)
    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", color="black",
                      markersize=math.sqrt(math.sqrt((v + 0.001) * 1000)), label=lab)
        for v, lab in ((1.0, "1"), (0.1, "0.1"), (0.01, "0.01"))
    ]
    leg1 = ax.legend(handles=handles, title="Allele frequency", ncol=3, loc="upper center", fontsize=7)
    ax.add_artist(leg1)
    ax.legend(
        handles=[
            mpatches.Patch(color=transparent_red, label="Reference"),
            mpatches.Patch(color=transparent_blue, label="Alternative (variant)"),
        ],
        loc="upper right", fontsize=7,
    )
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def _num_ok(tok):
    tok = str(tok).strip()
    if tok.lower() in _NA_TOKENS:
        return False
    try:
        float(tok)
        return True
    except ValueError:
        return False


def score_label(score_key):
    return "CRISTA" if score_key == "crista" else "CFD"


def plot_scatter_panels(df, cols, n=1000, include_crista=False):
    """Produce the key graphical-report scatter panels (section 2).

    Returns a list of (title, caption, data_uri). Panels:
      (a) top-N by CFD score
      (b) top-N by DELTA = ALT CFD - REF CFD (largest variant increase first)
      (c) top-N by CRISTA score  -- only when include_crista is True
    Every panel is guarded; a failing panel yields a placeholder URI (never
    aborts).
    """
    panels = []

    # candidate frame: drop on-/near-on-target rows (mm+b<=1), same as the table
    base = df
    if "mmb" in cols:
        base = base[_to_int_series(base[cols["mmb"]]) > 1]

    # (a) by CFD score
    def _cfd_sorted():
        cfd = _to_float_series(base[cols["cfd"]]).fillna(-1.0)
        return base.assign(_s=cfd).sort_values("_s", ascending=False).head(n)

    # (b) by DELTA = ALT - REF (descending)
    def _delta_sorted():
        alt = _to_float_series(base[cols["cfd_alt"]])
        ref = _to_float_series(base[cols["cfd_ref"]])
        delta = (alt - ref)
        return (
            base.assign(_delta=delta)
            .sort_values("_delta", ascending=False, na_position="last")
            .head(n)
        )

    # (c) by CRISTA
    def _crista_sorted():
        cr = _to_float_series(base[cols["crista"]]).fillna(-1.0)
        return base.assign(_s=cr).sort_values("_s", ascending=False).head(n)

    ref_key = cols.get("cfd_ref", "")
    alt_key = cols.get("cfd_alt", "")
    score_key = cols.get("cfd", "")

    # (a)
    try:
        uri = _cfd_style_scatter(
            _cfd_sorted(), cols, score_key, ref_key, alt_key,
            xlabel="Candidate off-target site (ranked by CFD)",
            title=f"Top {min(n, len(base))} candidates by CFD score",
        )
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: CFD scatter unavailable: {exc}\n")
        uri = _placeholder_uri("CFD scatter unavailable")
    panels.append((
        "By CFD score",
        "Candidate off-targets ranked by CFD (log x). Red = reference, blue = "
        "variant-created (sized by allele frequency); gray arrows connect the "
        "reference->variant pair for the same site; the top site is labeled with "
        "its rsID.",
        uri,
    ))

    # (b) only when the REF/ALT CFD columns exist
    if "cfd_ref" in cols and "cfd_alt" in cols:
        try:
            uri = _cfd_style_scatter(
                _delta_sorted(), cols, score_key, ref_key, alt_key,
                xlabel="Candidate off-target site (ranked by variant effect ALT-REF)",
                title=f"Top {min(n, len(base))} by variant-induced CFD increase (ALT-REF)",
            )
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"generate-report: delta scatter unavailable: {exc}\n")
            uri = _placeholder_uri("Delta scatter unavailable")
        panels.append((
            "By variant effect (ALT - REF CFD)",
            "Same axes, re-ranked by the variant-induced CFD change (ALT-REF, "
            "descending): the variants that most raise the CFD score come first, "
            "foregrounding the risk-relevant variant-created sites.",
            uri,
        ))

    # (c) CRISTA, only when computed
    if include_crista:
        cr_score = cols.get("crista", "")
        cr_ref = cols.get("crista_ref", "")
        cr_alt = cols.get("crista_alt", "")
        try:
            uri = _cfd_style_scatter(
                _crista_sorted(), cols, cr_score, cr_ref, cr_alt,
                xlabel="Candidate off-target site (ranked by CRISTA)",
                title=f"Top {min(n, len(base))} candidates by CRISTA score",
            )
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"generate-report: CRISTA scatter unavailable: {exc}\n")
            uri = _placeholder_uri("CRISTA scatter unavailable")
        panels.append((
            "By CRISTA score",
            "Independent CRISTA scoring model, ranked by CRISTA score. Included "
            "because CRISTA scores were computed for this run.",
            uri,
        ))

    return panels


def plot_population(df, cols, sample_superpop):
    """SECTION 3: SIMPLIFIED population view at the run's SELECTED parameters.

    Two panels over the FULL result set (no per-total-mm faceting):
      left  : reference vs variant-created off-target counts
      right : among variant-created, one bar per superpopulation (or per
              dataset when superpop mapping is unavailable).
    """
    variant, _reference, _ontarget = partition_masks(df, cols)
    variant_mask = variant
    n_variant = int(variant_mask.sum())
    n_reference = len(df) - n_variant  # (reference + on-target) as the non-variant bar

    group_counts = {}
    use_superpop = bool(sample_superpop)
    samples_col = cols.get("samples")
    if samples_col and samples_col in df.columns:
        for raw in df.loc[variant_mask, samples_col].astype(str):
            if _is_na(raw):
                continue
            seen = set()
            for sid in raw.split(","):
                sid = sid.strip()
                if _is_na(sid):
                    continue
                if use_superpop:
                    grp = sample_superpop.get(sid)
                    if grp is None:
                        grp = _dataset_of(sid)
                else:
                    grp = _dataset_of(sid)
                seen.add(grp)
            for grp in seen:
                group_counts[grp] = group_counts.get(grp, 0) + 1

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8))

    ax1.bar(
        ["reference", "variant-created"],
        [n_reference, n_variant],
        color=["#a0aec0", "#dd6b20"],
    )
    ax1.set_ylabel("number of off-targets")
    ax1.set_title("Reference vs variant-created")
    for i, val in enumerate([n_reference, n_variant]):
        ax1.text(i, val, f"{val:,}", ha="center", va="bottom", fontsize=9)

    if group_counts:
        labels = sorted(group_counts, key=lambda k: (-group_counts[k], k))
        values = [group_counts[k] for k in labels]
        ax2.bar(labels, values, color="#2b6cb0")
        ax2.set_ylabel("variant-created off-targets")
        title = "By superpopulation" if use_superpop else "By dataset (provenance)"
        ax2.set_title(title)
        if len(labels) > 4:
            ax2.tick_params(axis="x", rotation=45)
    else:
        ax2.axis("off")
        ax2.text(
            0.5, 0.5,
            "No per-population breakdown\n(no sample IDs / mapping available)",
            ha="center", va="center", fontsize=10,
        )
    fig.tight_layout()
    return _fig_to_data_uri(fig)


# --------------------------------------------------------------------------- #
# SECTION 4: recommended validation panel
# --------------------------------------------------------------------------- #
def build_validation_panel(df, cols):
    """Compute candidate counts at CFD / edit-distance thresholds and a suggested
    tiered panel. Returns a dict for rendering.

    Counting is over the OFF-TARGET set (on-target row mm+b==0 excluded), so
    thresholds report actionable candidates for a confirmation panel.
    """
    variant, _reference, ontarget = partition_masks(df, cols)
    offt = df[~ontarget]  # off-targets only

    cfd = _to_float_series(offt[cols["cfd"]]) if "cfd" in cols else pd.Series([], dtype=float)
    mmb = _to_int_series(offt[cols["mmb"]]) if "mmb" in cols else pd.Series([], dtype=int)
    var_off = variant[~ontarget]

    cfd_counts = [(t, int((cfd >= t).sum())) for t in CFD_THRESHOLDS]
    mmb_counts = [(t, int((mmb <= t).sum())) for t in MMB_THRESHOLDS]

    n_variant = int(var_off.sum())
    n_variant_cfd = int((var_off & (cfd >= PANEL_VARIANT_CFD_MIN)).sum())

    # suggested tiered panel = union of index sets over the OFF-TARGET frame
    idx = offt.index
    top_by_cfd = set(
        offt.assign(_s=cfd.fillna(-1.0))
        .sort_values("_s", ascending=False)
        .head(PANEL_TOP_N_BY_CFD)
        .index
    )
    all_mmb = set(idx[(mmb <= PANEL_MMB_MAX).values])
    all_var_cfd = set(idx[(var_off & (cfd >= PANEL_VARIANT_CFD_MIN)).values])
    panel = top_by_cfd | all_mmb | all_var_cfd

    return {
        "cfd_counts": cfd_counts,
        "mmb_counts": mmb_counts,
        "n_variant": n_variant,
        "n_variant_cfd": n_variant_cfd,
        "panel_size": len(panel),
        "tier_sizes": {
            "top_cfd": len(top_by_cfd),
            "mmb": len(all_mmb),
            "variant_cfd": len(all_var_cfd),
        },
    }


def render_validation_panel(vp):
    """Section 4 HTML: threshold tables + suggested panel + IND rationale."""
    cfd_rows = "".join(
        f"<tr><td>CFD &ge; {t}</td><td class='num'>{c:,}</td></tr>"
        for t, c in vp["cfd_counts"]
    )
    mmb_rows = "".join(
        f"<tr><td>mismatches + bulges &le; {t}</td><td class='num'>{c:,}</td></tr>"
        for t, c in vp["mmb_counts"]
    )
    ts = vp["tier_sizes"]
    rationale = (
        "IND rationale: no single metric captures off-target risk, so the "
        "suggested confirmation panel combines a CFD-score cutoff with an "
        "edit-distance (mismatch + bulge) cutoff. It is the union of the top "
        f"{PANEL_TOP_N_BY_CFD} candidates by CFD, every site within "
        f"{PANEL_MMB_MAX} mismatches + bulges of the protospacer, and every "
        f"variant-created site with CFD &ge; {PANEL_VARIANT_CFD_MIN}. This "
        "captures both the highest-scoring predicted cleavage sites and the "
        "near-cognate sequences that scoring models can under-weight, while "
        "specifically retaining the population-variant-created sites that a "
        "reference-only analysis would miss. The full integrated results table "
        "is provided below (and as a bundled download) so the panel can be "
        "freely re-subset for a given assay budget."
    )
    return f"""
<div class="panel-grid">
  <div>
    <table class="thr-table"><thead><tr><th>CFD threshold</th><th>Candidates</th></tr></thead>
    <tbody>{cfd_rows}</tbody></table>
  </div>
  <div>
    <table class="thr-table"><thead><tr><th>Edit-distance threshold</th><th>Candidates</th></tr></thead>
    <tbody>{mmb_rows}</tbody></table>
  </div>
</div>
<table class="thr-table" style="max-width:640px">
  <tbody>
    <tr><td>Variant-created off-targets (Not_found_in_REF)</td><td class="num">{vp['n_variant']:,}</td></tr>
    <tr><td>&hellip; of those with CFD &ge; {PANEL_VARIANT_CFD_MIN}</td><td class="num">{vp['n_variant_cfd']:,}</td></tr>
    <tr class="panel-hi"><td><strong>Suggested tiered panel &mdash; total sites</strong><br>
      <span class="caption">union of top-{PANEL_TOP_N_BY_CFD}-by-CFD ({ts['top_cfd']:,})
      &cup; all mm+b&le;{PANEL_MMB_MAX} ({ts['mmb']:,})
      &cup; variant-created CFD&ge;{PANEL_VARIANT_CFD_MIN} ({ts['variant_cfd']:,})</span></td>
      <td class="num"><strong>{vp['panel_size']:,}</strong></td></tr>
  </tbody>
</table>
<p class="caption">{rationale}</p>
"""


# --------------------------------------------------------------------------- #
# Top-1000 selection + HTML table + top1000.tsv
# --------------------------------------------------------------------------- #
def select_top(df, cols, n=1000):
    """Top-N off-targets by CFD desc, dropping on-/near-on-target rows.

    Filter mm+b > 1 to drop the on-target and near-on-target rows, matching
    CRISPRme_plots.py:527-530 filter_table.
    """
    work = df.copy()
    if "mmb" in cols:
        work = work[_to_int_series(work[cols["mmb"]]) > 1]
    if "cfd" in cols:
        work = work.assign(
            _cfd=pd.to_numeric(work[cols["cfd"]], errors="coerce").fillna(-1.0)
        ).sort_values("_cfd", ascending=False).drop(columns=["_cfd"])
    return work.head(n)


def _esc(value):
    if _is_na(value):
        return ""
    return html.escape(str(value))


def build_table_html(top_df, cols, has_crista):
    """Scrollable inline top-N table, sorted by CFD desc; no JS (opens offline).

    Columns: rank, chr, position, strand, aligned protospacer+PAM (ALT/REF
    fallback), MM, bulges, mm+b, CFD, [CRISTA if computed], REF/ALT origin, PAM
    creation, variant (rsID | genomic key), MAF (em-dash + footnote when blank),
    gene+distance.
    """
    headers = [
        "Rank", "Chr", "Position", "Strand", "Aligned protospacer+PAM (ALT)",
        "MM", "Bulges", "MM+B", "CFD",
    ]
    if has_crista:
        headers.append("CRISTA")
    headers += [
        "Origin", "PAM creation", "Variant (rsID | genomic)", "MAF", "Gene",
    ]
    head_html = "".join(f"<th>{h}</th>" for h in headers)

    maf_missing_seen = False
    body_rows = []
    for rank, (_, row) in enumerate(top_df.iterrows(), start=1):
        aligned = ""
        if "aln_alt" in cols:
            aligned = row.get(cols["aln_alt"], "")
        if _is_na(aligned) and "aln_ref" in cols:
            aligned = row.get(cols["aln_ref"], "")

        variant = _first_non_na(row.get(cols["rsid"])) if "rsid" in cols else None
        if variant is None and "var_genome" in cols:
            variant = _first_non_na(row.get(cols["var_genome"]))

        maf = _min_maf(row.get(cols["maf"])) if "maf" in cols else None
        if isinstance(maf, float):
            maf_txt = f"{maf:.2e}"
        else:
            maf_txt = "&mdash;"
            maf_missing_seen = True

        gene = row.get(cols["gene_name"]) if "gene_name" in cols else None
        gene_txt = "" if _is_na(gene) else _esc(gene)
        if gene_txt and "gene_dist" in cols:
            dist = row.get(cols["gene_dist"])
            if not _is_na(dist):
                gene_txt += f" ({_esc(dist)} kb)"

        cfd = pd.to_numeric(row.get(cols["cfd"]), errors="coerce") if "cfd" in cols else None
        cfd_txt = f"{cfd:.4f}" if pd.notna(cfd) else ""

        pam_creation = row.get(cols["pam_creation"]) if "pam_creation" in cols else None
        pam_txt = "" if _is_na(pam_creation) else _esc(pam_creation)

        cells = [
            str(rank),
            _esc(row.get(cols["chrom"])) if "chrom" in cols else "",
            _esc(row.get(cols["pos"])) if "pos" in cols else "",
            _esc(row.get(cols["strand"])) if "strand" in cols else "",
            f"<code>{_esc(aligned)}</code>",
            _esc(row.get(cols["mm"])) if "mm" in cols else "",
            _esc(row.get(cols["bulges"])) if "bulges" in cols else "",
            _esc(row.get(cols["mmb"])) if "mmb" in cols else "",
            cfd_txt,
        ]
        if has_crista:
            crista = (
                pd.to_numeric(row.get(cols["crista"]), errors="coerce")
                if "crista" in cols else None
            )
            cells.append(f"{crista:.4f}" if pd.notna(crista) else "")
        cells += [
            _esc(row.get(cols["origin"])).upper() if "origin" in cols else "",
            pam_txt,
            _esc(variant),
            maf_txt,
            gene_txt,
        ]
        body_rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")

    footnote = ""
    if maf_missing_seen:
        footnote = (
            '<p class="caption">MAF shown as &mdash; when no allele frequency is '
            "annotated for the variant (the variant still exists in the callset; "
            "only its population frequency is unavailable from the AF source).</p>"
        )

    table = (
        '<div class="ottable-wrap"><table class="ottable">'
        f"<thead><tr>{head_html}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table></div>"
    )
    return table + footnote


def write_top1000_tsv(top_df, path):
    """Write the exact top-N rows (all original columns) to a TSV."""
    top_df.to_csv(path, sep="\t", index=False)


# --------------------------------------------------------------------------- #
# HTML assembly
# --------------------------------------------------------------------------- #
_CSS = """
:root { color-scheme: light; }
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       max-width: 1180px; margin: 24px auto; padding: 0 20px; color: #1a202c;
       line-height: 1.45; }
h1 { font-size: 1.6em; margin-bottom: 0.1em; }
h2 { font-size: 1.25em; margin-top: 1.8em; border-bottom: 1px solid #e2e8f0;
     padding-bottom: 0.2em; }
.subtitle { color: #4a5568; margin-top: 0; }
.caption { color: #718096; font-size: 0.86em; margin: 0.3em 0 1.2em; }
.summary-grid { display: flex; flex-wrap: wrap; gap: 24px; align-items: flex-start; }
.summary-card { flex: 0 0 auto; }
.matrix-card { flex: 1 1 460px; min-width: 420px; }
.matrix-title { font-weight: 600; margin-bottom: 0.4em; }
table.summary-table { border-collapse: collapse; margin: 0.2em 0; width: 100%;
                      max-width: 520px; }
table.summary-table td { border: 1px solid #e2e8f0; padding: 6px 10px;
                         vertical-align: top; }
table.summary-table td.k { background: #f7fafc; font-weight: 600; width: 250px; }
.matrix-wrap { overflow-x: auto; }
table.matrix { border-collapse: collapse; font-size: 0.82em; }
table.matrix th, table.matrix td { border: 1px solid #e2e8f0; padding: 4px 9px;
                                   text-align: right; }
table.matrix th { background: #f7fafc; }
table.matrix td.grp { background: #edf2f7; font-weight: 600; text-align: center;
                      vertical-align: middle; }
table.matrix td.tot { font-weight: 600; background: #f7fafc; }
.grp-total { color: #718096; font-weight: 400; font-size: 0.9em; }
.plot { margin: 0.6em 0 0.6em; }
.plot img { max-width: 100%; height: auto; border: 1px solid #edf2f7; border-radius: 4px; }
a.download { display: inline-block; background: #2b6cb0; color: #fff;
             padding: 8px 16px; border-radius: 4px; text-decoration: none;
             margin: 4px 8px 4px 0; }
a.download:hover { background: #2c5282; }
.panel-grid { display: flex; flex-wrap: wrap; gap: 24px; }
table.thr-table { border-collapse: collapse; margin: 0.4em 0; }
table.thr-table th, table.thr-table td { border: 1px solid #e2e8f0; padding: 5px 12px; }
table.thr-table th { background: #f7fafc; text-align: left; }
table.thr-table td.num { text-align: right; font-variant-numeric: tabular-nums; }
tr.panel-hi td { background: #ebf8ff; }
.ottable-wrap { max-height: 560px; overflow: auto; border: 1px solid #ccc;
                border-radius: 4px; }
table.ottable { border-collapse: collapse; width: 100%; font-size: 0.8em; }
table.ottable th, table.ottable td { padding: 4px 8px; text-align: left;
                                     white-space: nowrap; }
table.ottable thead th { position: sticky; top: 0; background: #fff;
                         border-bottom: 2px solid #cbd5e0; z-index: 1; }
table.ottable tbody tr:nth-child(even) { background: #f6f6f6; }
table.ottable code { font-family: SFMono-Regular, Menlo, Consolas, monospace;
                     font-size: 0.95em; }
footer { margin-top: 2.5em; color: #4a5568; font-size: 0.82em;
         border-top: 1px solid #e2e8f0; padding-top: 0.9em; }
footer .disclaimer { color: #742a2a; background: #fff5f5; border: 1px solid #fed7d7;
                     border-radius: 4px; padding: 10px 14px; margin-top: 0.8em; }
""".strip()

DISCLAIMER = (
    "This report is provided for research purposes only. CRISPRme off-target "
    "predictions are computational, may contain false positives and false "
    "negatives, and are NOT a substitute for experimental validation. The "
    "authors, contributors, and their institutions make no warranties, express "
    "or implied, and accept no liability for any decision, result, clinical or "
    "regulatory application, or other use arising from this report."
)


def render_html(
    job_id, summary_matrix_html, scatter_panels, pop_uri, validation_html,
    table_html, tsv_gz_name, top1000_name, footer_html,
):
    scatter_html = []
    for title, caption, uri in scatter_panels:
        scatter_html.append(
            f'<h3 style="margin:0.8em 0 0.2em">{_esc(title)}</h3>'
            f'<div class="plot"><img alt="{_esc(title)}" src="{uri}"></div>'
            f'<p class="caption">{caption}</p>'
        )
    scatter_block = "\n".join(scatter_html)

    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(job_id)} CRISPRme off-target report</title>
<style>{_CSS}</style>
</head><body>

<h1>CRISPRme off-target prediction report</h1>
<p class="subtitle">IND briefing-book digest &mdash; targeted-NGS (rhAMP-Seq)
confirmation panel design</p>

<h2>1. Summary</h2>
{summary_matrix_html}

<h2>2. Key graphical report</h2>
<p class="caption">Reference vs variant off-target scores across the top-ranked
candidates (CRISPRme-paper style). The same scatter is shown under multiple
rankings so the highest-scoring sites and the highest-variant-effect sites are
both foregrounded.</p>
{scatter_block}

<h2>3. Reference vs population origin</h2>
<div class="plot"><img alt="Reference vs population origin" src="{pop_uri}"></div>
<p class="caption">Left: reference vs variant-created off-target counts at the
run's selected parameters. Right: variant-created sites broken down by
superpopulation (or dataset provenance when superpopulation mapping is
unavailable). A site is counted once per group with at least one carrier.</p>

<h2>4. Recommended validation panel</h2>
{validation_html}

<h2>5. Downloads</h2>
<p>
  <a class="download" href="{_esc(tsv_gz_name)}" download>Complete integrated results (TSV, gzip)</a>
  <a class="download" href="{_esc(top1000_name)}" download>Top-1000 off-targets (TSV)</a>
</p>
<p class="caption">Both files are bundled alongside this HTML in the same ZIP;
the links resolve after unzipping on any machine. The top-1000 TSV contains
exactly the rows shown in the table below.</p>

<h2>6. Top 1000 off-targets (by CFD)</h2>
{table_html}

{footer_html}

</body></html>
"""


def build_footer(meta, version, tsv_basename):
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    crisprme_version = version or meta.get("version") or "n/a"
    return f"""<footer>
<p>CRISPRme version: {_esc(crisprme_version)}
&nbsp;&middot;&nbsp; report generator v{_esc(REPORT_GENERATOR_VERSION)}
&nbsp;&middot;&nbsp; generated {_esc(stamp)}
&nbsp;&middot;&nbsp; source: {_esc(tsv_basename)}</p>
<p class="disclaimer">{_esc(DISCLAIMER)}</p>
</footer>"""


# --------------------------------------------------------------------------- #
# TSV resolution
# --------------------------------------------------------------------------- #
def resolve_integrated_tsv(result_dir):
    """Find the finalized integrated_results TSV in a result folder.

    Matches results_page.py:2298-2300 (glob ``*integrated_results.tsv``). Prefers
    the renamed ``<guide>+..._integrated_results.tsv`` over the intermediate
    ``.bestMerge.txt.integrated_results.tsv``.
    """
    candidates = sorted(
        glob.glob(os.path.join(result_dir, "*integrated_results.tsv"))
    )
    candidates += sorted(
        glob.glob(os.path.join(result_dir, "*integrated_results.tsv.gz"))
    )
    if not candidates:
        return None
    candidates.sort(key=lambda p: ("bestMerge" in os.path.basename(p), p))
    return candidates[0]


def _job_id_from(result_dir, tsv_path):
    if result_dir:
        jid = os.path.basename(os.path.normpath(result_dir))
        if jid:
            return jid
    base = os.path.basename(tsv_path)
    for suffix in ("_integrated_results.tsv.gz", "_integrated_results.tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def _default_samplesid_dir(result_dir):
    """Best-effort samplesID directory near the install (never required)."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = []
    if result_dir:
        candidates.append(os.path.join(result_dir, "samplesIDs"))
        candidates.append(os.path.join(result_dir, "..", "..", "samplesIDs"))
    candidates.append(os.path.join(here, "..", "samplesIDs"))
    candidates.append(os.path.join(here, "..", "test", "data", "samplesIDs"))
    candidates.append(os.path.join(os.getcwd(), "samplesIDs"))
    for cand in candidates:
        if cand and os.path.isdir(cand):
            return cand
    return None


def _package_version():
    """Best-effort CRISPRme package version (footer fallback)."""
    try:
        from importlib.metadata import PackageNotFoundError, version as _v

        try:
            return _v("CRISPRme")
        except PackageNotFoundError:
            return None
    except Exception:  # noqa: BLE001
        return None


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def build_report(
    result_dir=None,
    integrated_tsv=None,
    out_zip=None,
    samplesid_dir=None,
    params_override=None,
    top_n=1000,
):
    """Build ``<jobid>_report.zip`` (report.html + integrated_results.tsv.gz +
    top1000.tsv).

    Parameters
    ----------
    result_dir : str, optional
        A CRISPRme result folder. The TSV, .Params.txt, .version.txt, the acfd
        specificity sidecar and the output ZIP location are all resolved from
        here when given.
    integrated_tsv : str, optional
        An explicit integrated_results TSV (``.tsv`` or ``.tsv.gz``). Required
        when ``result_dir`` has no discoverable TSV.
    out_zip : str, optional
        Output ZIP path. Defaults to ``<result_dir or cwd>/<jobid>_report.zip``.
    samplesid_dir : str, optional
        Directory of samplesID files for the superpopulation plot. Auto-detected
        near the install when omitted; the plot degrades to per-dataset bars if
        unresolved.
    params_override : dict, optional
        Override / supply summary fields when no .Params.txt exists (e.g. the
        synthesized summary for a bare sample TSV).
    top_n : int
        Table / top-N TSV size (default 1000).

    Returns
    -------
    str
        Absolute path to the written ZIP.
    """
    if integrated_tsv is None:
        if not result_dir:
            raise ValueError("Provide result_dir or integrated_tsv")
        integrated_tsv = resolve_integrated_tsv(result_dir)
        if integrated_tsv is None:
            raise FileNotFoundError(
                f"No *integrated_results.tsv found in {result_dir}"
            )
    integrated_tsv = os.path.abspath(integrated_tsv)
    if not os.path.isfile(integrated_tsv):
        raise FileNotFoundError(integrated_tsv)

    job_id = _job_id_from(result_dir, integrated_tsv)

    if out_zip is None:
        base_dir = result_dir if result_dir else os.path.dirname(integrated_tsv)
        out_zip = os.path.join(os.path.abspath(base_dir), f"{job_id}_report.zip")
    out_zip = os.path.abspath(out_zip)

    if samplesid_dir is None:
        samplesid_dir = _default_samplesid_dir(result_dir)

    # ---- load data (all columns as str; resolve by name) -------------------
    compression = "gzip" if integrated_tsv.endswith(".gz") else None
    df = pd.read_csv(
        integrated_tsv,
        sep="\t",
        dtype=str,
        compression=compression,
        na_filter=False,
        low_memory=False,
    )
    cols = _resolve(df.columns, list(_COLS.keys()))

    meta = build_summary_meta(
        result_dir, integrated_tsv, df, cols, params_override=params_override
    )
    version = meta.get("version") or _package_version()
    has_crista = crista_computed(df, cols)

    spec_score = read_specificity_score(result_dir, job_id, meta["guides"])

    fn = _parse_results_filename(integrated_tsv)
    sample_superpop = load_sample_superpop(samplesid_dir, fn.get("datasets", ""))

    top_df = select_top(df, cols, n=top_n)

    # ---- SECTION 1: summary + matrix ---------------------------------------
    try:
        matrix = build_mmb_matrix(df, cols, meta)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: matrix unavailable: {exc}\n")
        matrix = None
    try:
        summary_matrix_html = render_summary_and_matrix(meta, spec_score, matrix)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: summary section unavailable: {exc}\n")
        summary_matrix_html = "<p>Summary unavailable.</p>"

    # ---- SECTION 2: scatter panels -----------------------------------------
    try:
        scatter_panels = plot_scatter_panels(
            df, cols, n=top_n, include_crista=has_crista
        )
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: scatter panels unavailable: {exc}\n")
        scatter_panels = [
            ("Graphical report", "Plot unavailable.", _placeholder_uri("Plot unavailable"))
        ]

    # ---- SECTION 3: population plot -----------------------------------------
    try:
        pop_uri = plot_population(df, cols, sample_superpop)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: population plot unavailable: {exc}\n")
        pop_uri = _placeholder_uri("Population plot unavailable")

    # ---- SECTION 4: validation panel ---------------------------------------
    try:
        vp = build_validation_panel(df, cols)
        validation_html = render_validation_panel(vp)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: validation panel unavailable: {exc}\n")
        validation_html = "<p>Validation panel unavailable.</p>"

    # ---- SECTION 6: top-1000 table -----------------------------------------
    try:
        table_html = build_table_html(top_df, cols, has_crista)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: table unavailable: {exc}\n")
        table_html = "<p>Top-1000 table unavailable.</p>"

    # ---- SECTION 7: footer --------------------------------------------------
    footer_html = build_footer(meta, version, os.path.basename(integrated_tsv))

    tsv_gz_name = "integrated_results.tsv.gz"
    top1000_name = "top1000.tsv"
    html_doc = render_html(
        job_id, summary_matrix_html, scatter_panels, pop_uri, validation_html,
        table_html, tsv_gz_name, top1000_name, footer_html,
    )

    # ---- stage the 3 files and zip -j (flat) -------------------------------
    staging = tempfile.mkdtemp(prefix="crisprme_report_")
    try:
        html_path = os.path.join(staging, "report.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(html_doc)

        top1000_path = os.path.join(staging, top1000_name)
        try:
            write_top1000_tsv(top_df, top1000_path)
        except Exception as exc:  # noqa: BLE001 - never abort on the bundled TSV
            sys.stderr.write(f"generate-report: top1000.tsv unavailable: {exc}\n")
            with open(top1000_path, "w") as handle:
                handle.write("\t".join(df.columns) + "\n")

        gz_path = os.path.join(staging, tsv_gz_name)
        if integrated_tsv.endswith(".gz"):
            shutil.copyfile(integrated_tsv, gz_path)
        else:
            with open(integrated_tsv, "rb") as src, gzip.open(gz_path, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)

        if os.path.exists(out_zip):
            os.remove(out_zip)
        os.makedirs(os.path.dirname(out_zip), exist_ok=True)
        # `zip -j` matches submit_job_automated_new_multiple_vcfs.sh:1187 and
        # flat-decompresses to exactly the 3 files. Fall back to Python zipfile
        # if the `zip` binary is unavailable.
        try:
            subprocess.run(
                ["zip", "-j", "-q", out_zip, html_path, gz_path, top1000_path],
                check=True,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            import zipfile

            with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(html_path, "report.html")
                zf.write(gz_path, tsv_gz_name)
                zf.write(top1000_path, top1000_name)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    return out_zip


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="generate-report",
        description="Build a self-contained, shareable CRISPRme off-target report ZIP.",
    )
    parser.add_argument(
        "--result-dir",
        help="CRISPRme result folder (Results/<jobid>). The TSV, .Params.txt, "
        ".version.txt and default output location are resolved from here.",
    )
    parser.add_argument(
        "--integrated-results",
        dest="integrated_results",
        help="Explicit integrated_results TSV (.tsv or .tsv.gz). Required if "
        "--result-dir has no discoverable TSV.",
    )
    parser.add_argument(
        "--samplesID-dir",
        dest="samplesid_dir",
        help="Directory of samplesID files for the superpopulation plot "
        "(auto-detected when omitted).",
    )
    parser.add_argument(
        "--output",
        "--out",
        dest="output",
        help="Output ZIP path (default: <result-dir>/<jobid>_report.zip).",
    )
    return parser


def main(argv=None):
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    if not args.result_dir and not args.integrated_results:
        parser.error("provide --result-dir or --integrated-results")
    out = build_report(
        result_dir=args.result_dir,
        integrated_tsv=args.integrated_results,
        out_zip=args.output,
        samplesid_dir=args.samplesid_dir,
    )
    sys.stdout.write(f"Report written: {out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
