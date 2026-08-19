#!/usr/bin/env python
"""Self-contained, shareable CRISPRme off-target report.

Given a CRISPRme result folder (or a bare ``*integrated_results.tsv``), this
module produces a single, easily-transferable ZIP::

    <jobid>_report.zip
      report.html                 # self-contained: base64 PNG plots, inline
                                   #   top-1000 table, inline CSS, opens offline
      integrated_results.tsv.gz   # the full results, the HTML links to it
                                   #   with a RELATIVE href (resolves post-unzip)

The report is a portable *digest* of the full interactive CRISPRme website
result (personal risk cards, etc. stay in the website). It is meant for a
collaborator preparing an IND briefing book -- e.g. designing a targeted-NGS
rhAMP-Seq confirmation panel from the predicted off-targets.

Design goals / robustness posture
---------------------------------
* Pure stdlib + matplotlib + pandas + numpy (the pipeline conda env). No Dash,
  no Jinja, no network. The HTML has no <script> and no external <link>, so it
  opens with ``file://`` on any machine.
* Columns are selected BY NAME from the header (never fixed indices), because
  the dict-less 85-col schema and the dict-based schema differ in column set /
  order. Missing columns (CRISTA, per-dataset, annotation on dict-less) are
  optional -- the report degrades, it never crashes.
* Every plot is wrapped in try/except; a failing plot inlines a small
  "plot unavailable" placeholder rather than failing the report (mirrors the
  issue-#143 "preserve results" posture used in the pipeline).

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
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime

import matplotlib

# never touch an X server (pipeline convention, CRISPRme_plots.py:21)
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# --------------------------------------------------------------------------- #
# Column-name resolution helpers
# --------------------------------------------------------------------------- #
# The "highest_CFD" projection is the one the report uses (CFD is the primary
# ranking the IND reviewer cares about). We resolve every column by its base
# name against whichever suffixes a given schema uses.
_PROJ = "(highest_CFD)"

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


def build_summary(result_dir, tsv_path, df, cols, params_override=None):
    """Assemble the top-of-page summary rows (mirrors the web result summary)."""
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

    # guide(s): prefer the actual TSV values (handles multi-guide runs)
    guides = []
    if "guide" in cols and cols["guide"] in df.columns:
        guides = [g for g in df[cols["guide"]].dropna().unique().tolist()]
    if not guides and fn.get("guide"):
        guides = [fn["guide"]]
    guide_display = ", ".join(guides) if guides else "n/a"

    pam = params.get("Pam") or fn.get("pam") or "n/a"
    nuclease = params.get("Nuclease") or "n/a"
    genome = params.get("Genome_selected") or params.get("Genome_ref") or fn.get(
        "genome"
    ) or "n/a"
    # variant datasets: prefer filename decode (explicit), else Ref_comp hint
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
    mm = params.get("Mismatches") or fn.get("mm") or "n/a"
    bdna = params.get("DNA")
    brna = params.get("RNA")
    bmax = params.get("Max_bulges") or fn.get("bmax") or "n/a"
    max_edits = params.get("Max_total_edits") or "n/a"

    # counts straight from the data (single source of truth)
    n_rows = len(df)
    n_ontarget = 0
    n_variant = 0
    n_ref = 0
    if "mmb" in cols:
        n_ontarget = int((_to_int_series(df[cols["mmb"]]) == 0).sum())
    if "not_in_ref" in cols:
        n_variant = int(
            df[cols["not_in_ref"]].astype(str).str.strip().str.lower().eq("y").sum()
        )
    if "origin" in cols:
        n_ref = int(
            df[cols["origin"]].astype(str).str.strip().str.lower().eq("ref").sum()
        )
    n_offtarget = n_rows - n_ontarget

    rows = [
        ("Guide RNA (spacer+PAM)", guide_display),
        ("Nuclease", nuclease),
        ("PAM", pam),
        ("Genome", genome),
        ("Variant dataset(s)", datasets),
        ("Mismatches", str(mm)),
        (
            "Bulges (DNA / RNA)",
            f"{bdna if bdna is not None else 'n/a'} / "
            f"{brna if brna is not None else 'n/a'}"
            + (f" (max bulges {bmax})" if bmax != "n/a" else ""),
        ),
        ("Max total edits", str(max_edits)),
        ("Total off-targets found", f"{n_offtarget:,}"),
        ("On-target site(s)", f"{n_ontarget:,}"),
        ("Variant-created off-targets", f"{n_variant:,}"),
        ("Reference off-targets", f"{n_ref:,}"),
        ("Date", date or "n/a"),
        ("CRISPRme version", version or "n/a"),
    ]
    return rows


def _to_int_series(series):
    return pd.to_numeric(series, errors="coerce").fillna(-1).astype(int)


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


def plot_top1000_by_cfd(top_df, cols):
    """PLOT 2a-i: distribution of the top-1000 CFD scores (histogram)."""
    fig, ax = plt.subplots(figsize=(7, 3.6))
    scores = pd.to_numeric(top_df[cols["cfd"]], errors="coerce").dropna()
    scores = scores[(scores >= 0) & (scores <= 1)]
    ax.hist(scores, bins=20, range=(0, 1), color="#2b6cb0", edgecolor="white")
    ax.set_xlabel("CFD score")
    ax.set_ylabel("number of off-targets")
    ax.set_title(f"Top {len(top_df)} candidates by CFD score")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def plot_top1000_by_mmb(top_df, cols):
    """PLOT 2a-ii: top-1000 grouped by mismatches+bulges, stacked ref vs variant."""
    mmb = _to_int_series(top_df[cols["mmb"]])
    if "origin" in cols:
        origin = top_df[cols["origin"]].astype(str).str.strip().str.lower()
    else:
        origin = pd.Series(["ref"] * len(top_df), index=top_df.index)
    is_alt = origin.eq("alt")

    max_mmb = int(mmb.max()) if len(mmb) and mmb.max() >= 0 else 0
    buckets = list(range(0, max_mmb + 1))
    ref_counts = [int(((mmb == b) & (~is_alt)).sum()) for b in buckets]
    alt_counts = [int(((mmb == b) & (is_alt)).sum()) for b in buckets]

    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.bar(buckets, ref_counts, color="#a0aec0", label="reference")
    ax.bar(
        buckets,
        alt_counts,
        bottom=ref_counts,
        color="#dd6b20",
        label="variant-created",
    )
    ax.set_xlabel("mismatches + bulges")
    ax.set_ylabel("number of off-targets")
    ax.set_title(f"Top {len(top_df)} candidates by edit distance")
    ax.set_xticks(buckets)
    ax.legend(frameon=False)
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def plot_population(df, cols, sample_superpop):
    """PLOT 2b: SIMPLIFIED population view at the run's SELECTED parameters.

    Two panels over the FULL result set (no per-total-mm faceting):
      left  : reference vs variant-created off-target counts
      right : among variant-created, one bar per superpopulation (or per
              dataset when superpop mapping is unavailable).
    """
    # left panel: reference vs variant-created
    if "not_in_ref" in cols:
        variant_mask = df[cols["not_in_ref"]].astype(str).str.strip().str.lower().eq("y")
    elif "origin" in cols:
        variant_mask = df[cols["origin"]].astype(str).str.strip().str.lower().eq("alt")
    else:
        variant_mask = pd.Series([False] * len(df), index=df.index)
    n_variant = int(variant_mask.sum())
    n_reference = len(df) - n_variant

    # right panel: per-superpopulation (or per-dataset fallback) among variants.
    # A site is counted ONCE per superpop that has >=1 carrier (seen-set
    # semantics, process_summaries.py:154/163-164).
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
                        grp = _dataset_of(sid)  # keep provenance for unknowns
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
        title = (
            "By superpopulation" if use_superpop else "By dataset (provenance)"
        )
        ax2.set_title(title)
        if len(labels) > 4:
            ax2.tick_params(axis="x", rotation=45)
    else:
        ax2.axis("off")
        ax2.text(
            0.5,
            0.5,
            "No per-population breakdown\n(no sample IDs / mapping available)",
            ha="center",
            va="center",
            fontsize=10,
        )
    fig.tight_layout()
    return _fig_to_data_uri(fig)


# --------------------------------------------------------------------------- #
# Top-1000 selection + HTML table
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
        ).sort_values("_cfd", ascending=False)
    return work.head(n)


def _esc(value):
    if _is_na(value):
        return ""
    return html.escape(str(value))


def build_table_html(top_df, cols):
    """Scrollable inline top-N table, sorted by CFD desc; no JS (opens offline)."""
    headers = [
        "Rank",
        "Chr",
        "Position",
        "Strand",
        "Aligned protospacer+PAM (ALT)",
        "MM",
        "Bulges",
        "CFD",
        "Origin",
        "Variant (rsID | genomic)",
        "MAF",
        "Gene",
    ]
    head_html = "".join(f"<th>{h}</th>" for h in headers)

    body_rows = []
    for rank, (_, row) in enumerate(top_df.iterrows(), start=1):
        # aligned ALT with fall back to REF when ALT is NA
        aligned = ""
        if "aln_alt" in cols:
            aligned = row.get(cols["aln_alt"], "")
        if _is_na(aligned) and "aln_ref" in cols:
            aligned = row.get(cols["aln_ref"], "")

        # variant cell: rsID if any non-NA, else genomic
        variant = _first_non_na(row.get(cols["rsid"])) if "rsid" in cols else None
        if variant is None and "var_genome" in cols:
            variant = _first_non_na(row.get(cols["var_genome"]))

        maf = _min_maf(row.get(cols["maf"])) if "maf" in cols else None
        maf_txt = f"{maf:.2e}" if isinstance(maf, float) else ""

        gene = row.get(cols["gene_name"]) if "gene_name" in cols else None
        gene_txt = "" if _is_na(gene) else _esc(gene)
        if gene_txt and "gene_dist" in cols:
            dist = row.get(cols["gene_dist"])
            if not _is_na(dist):
                gene_txt += f" ({_esc(dist)} kb)"

        cfd = pd.to_numeric(row.get(cols["cfd"]), errors="coerce") if "cfd" in cols else None
        cfd_txt = f"{cfd:.4f}" if pd.notna(cfd) else ""

        cells = [
            str(rank),
            _esc(row.get(cols["chrom"])) if "chrom" in cols else "",
            _esc(row.get(cols["pos"])) if "pos" in cols else "",
            _esc(row.get(cols["strand"])) if "strand" in cols else "",
            f"<code>{_esc(aligned)}</code>",
            _esc(row.get(cols["mm"])) if "mm" in cols else "",
            _esc(row.get(cols["bulges"])) if "bulges" in cols else "",
            cfd_txt,
            _esc(row.get(cols["origin"])).upper() if "origin" in cols else "",
            _esc(variant),
            maf_txt,
            gene_txt,
        ]
        body_rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")

    return (
        '<div class="ottable-wrap"><table class="ottable">'
        f"<thead><tr>{head_html}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table></div>"
    )


# --------------------------------------------------------------------------- #
# HTML assembly
# --------------------------------------------------------------------------- #
_CSS = """
:root { color-scheme: light; }
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       max-width: 1100px; margin: 24px auto; padding: 0 20px; color: #1a202c;
       line-height: 1.45; }
h1 { font-size: 1.6em; margin-bottom: 0.1em; }
h2 { font-size: 1.25em; margin-top: 1.6em; border-bottom: 1px solid #e2e8f0;
     padding-bottom: 0.2em; }
.subtitle { color: #4a5568; margin-top: 0; }
.caption { color: #718096; font-size: 0.88em; margin: 0.2em 0 1.2em; }
table.summary-table { border-collapse: collapse; margin: 1em 0; width: 100%;
                      max-width: 640px; }
table.summary-table td { border: 1px solid #e2e8f0; padding: 6px 10px;
                         vertical-align: top; }
table.summary-table td.k { background: #f7fafc; font-weight: 600; width: 240px; }
.plot { margin: 0.6em 0 1.4em; }
.plot img { max-width: 100%; height: auto; border: 1px solid #edf2f7; border-radius: 4px; }
a.download { display: inline-block; background: #2b6cb0; color: #fff;
             padding: 8px 16px; border-radius: 4px; text-decoration: none; }
a.download:hover { background: #2c5282; }
.ottable-wrap { max-height: 520px; overflow-y: auto; border: 1px solid #ccc;
                border-radius: 4px; }
table.ottable { border-collapse: collapse; width: 100%; font-size: 0.82em; }
table.ottable th, table.ottable td { padding: 4px 8px; text-align: left;
                                     white-space: nowrap; }
table.ottable thead th { position: sticky; top: 0; background: #fff;
                         border-bottom: 2px solid #cbd5e0; z-index: 1; }
table.ottable tbody tr:nth-child(even) { background: #f6f6f6; }
table.ottable code { font-family: SFMono-Regular, Menlo, Consolas, monospace;
                     font-size: 0.95em; }
footer { margin-top: 2.5em; color: #a0aec0; font-size: 0.8em;
         border-top: 1px solid #e2e8f0; padding-top: 0.8em; }
""".strip()


def render_html(job_id, summary_rows, plots, table_html, tsv_gz_name):
    summary_html = "".join(
        f'<tr><td class="k">{_esc(k)}</td><td>{_esc(v)}</td></tr>'
        for k, v in summary_rows
    )
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(job_id)} CRISPRme off-target report</title>
<style>{_CSS}</style>
</head><body>

<h1>CRISPRme off-target prediction report</h1>
<p class="subtitle">Prepared for targeted-NGS (rhAMP-Seq) confirmation panel design</p>

<table class="summary-table"><tbody>{summary_html}</tbody></table>
<p class="caption">This portable digest uses the highest-CFD projection of each
off-target cluster. It is a shareable summary of the full interactive CRISPRme
result (which additionally provides per-sample risk cards and interactive
browsing).</p>

<h2>Graphical report (top-1000 candidates)</h2>
<div class="plot"><img alt="Top-1000 by CFD score" src="{plots['cfd']}"></div>
<p class="caption">Distribution of CFD scores across the top-1000 off-target
candidates (ranked by CFD, on-/near-on-target rows excluded).</p>
<div class="plot"><img alt="Top-1000 by edit distance" src="{plots['mmb']}"></div>
<p class="caption">Top-1000 candidates grouped by total mismatches + bulges,
stacked by reference vs variant-created origin.</p>

<h2>Reference vs population origin</h2>
<div class="plot"><img alt="Reference vs population origin" src="{plots['pop']}"></div>
<p class="caption">Left: how many candidate sites exist in the reference vs.
those created only by a variant, at the run's selected parameters. Right:
variant-created sites broken down by superpopulation (or dataset provenance when
superpopulation mapping is unavailable). A site is counted once per group with at
least one carrier.</p>

<h2>Full results</h2>
<p><a class="download" href="{_esc(tsv_gz_name)}" download>Download the complete
integrated results (TSV, gzip)</a></p>
<p class="caption">Bundled alongside this HTML in the same ZIP; the link resolves
after unzipping on any machine.</p>

<h2>Top 1000 off-targets (by CFD)</h2>
{table_html}

<footer>Generated by CRISPRme generate-report. Off-target predictions are
computational candidates; experimental confirmation (e.g. targeted NGS /
rhAMP-Seq) is required before any clinical interpretation.</footer>

</body></html>
"""


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
    # deprioritise intermediate bestMerge names
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
    """Build ``<jobid>_report.zip`` (report.html + integrated_results.tsv.gz).

    Parameters
    ----------
    result_dir : str, optional
        A CRISPRme result folder. The TSV, .Params.txt, .version.txt and the
        output ZIP location are all resolved from here when given.
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

    summary_rows = build_summary(
        result_dir, integrated_tsv, df, cols, params_override=params_override
    )

    fn = _parse_results_filename(integrated_tsv)
    sample_superpop = load_sample_superpop(samplesid_dir, fn.get("datasets", ""))

    top_df = select_top(df, cols, n=top_n)

    # ---- plots (each guarded so a graphics error never fails the run) ------
    plots = {}
    for key, fn_plot in (
        ("cfd", lambda: plot_top1000_by_cfd(top_df, cols)),
        ("mmb", lambda: plot_top1000_by_mmb(top_df, cols)),
        ("pop", lambda: plot_population(df, cols, sample_superpop)),
    ):
        try:
            plots[key] = fn_plot()
        except Exception as exc:  # noqa: BLE001 - report must never crash on plots
            sys.stderr.write(f"generate-report: plot '{key}' unavailable: {exc}\n")
            plots[key] = _placeholder_uri(f"Plot unavailable ({key})")

    try:
        table_html = build_table_html(top_df, cols)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: table unavailable: {exc}\n")
        table_html = "<p>Top-1000 table unavailable.</p>"

    tsv_gz_name = "integrated_results.tsv.gz"
    html_doc = render_html(job_id, summary_rows, plots, table_html, tsv_gz_name)

    # ---- stage the 2 files and zip -j (flat) -------------------------------
    staging = tempfile.mkdtemp(prefix="crisprme_report_")
    try:
        html_path = os.path.join(staging, "report.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(html_doc)

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
        # flat-decompresses to exactly the 2 files. Fall back to Python zipfile
        # if the `zip` binary is unavailable.
        try:
            subprocess.run(
                ["zip", "-j", "-q", out_zip, html_path, gz_path],
                check=True,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            import zipfile

            with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(html_path, "report.html")
                zf.write(gz_path, tsv_gz_name)
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
