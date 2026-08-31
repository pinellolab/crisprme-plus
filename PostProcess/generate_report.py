#!/usr/bin/env python
"""Self-contained, shareable CRISPRme off-target report (v2.4 -- IND briefing book).

Given a CRISPRme result folder (or a bare ``*integrated_results.tsv``), this
module produces a single, easily-transferable ZIP::

    <jobid>_report.zip
      report.html                 # self-contained: base64 PNG plots, inline
                                   #   top-1000 table, inline CSS, opens offline
      integrated_results.tsv.gz   # the full RAW results (all 85 columns), the
                                   #   HTML links to it with a RELATIVE href
      top1000.tsv                 # the top-1000-by-CFD rows shown in the table,
                                   #   CURATED columns, its own RELATIVE link
      panel_top100.tsv            # the hybrid worst-case top-100 validation
                                   #   panel (section 4), CURATED columns
      cfd_ge_0.50.tsv ...         # per-tier subsets (CFD>=0.5/0.2/0.05,
      mmb_le_1.tsv ...            #   mm+b<=1/2/3/4, variant_created), each in
      variant_created.tsv         #   CURATED columns, only when non-empty

Curated columns (report v2.4)
-----------------------------
ONE curated, readable, Excel-ready column set (``CURATED_COLUMNS``) is shared by
BOTH the in-report top-1000 table AND every exported download file (top1000.tsv,
panel_top100.tsv, and every per-tier TSV). It carries the ranking columns PLUS
the annotation columns (gene, distance, GENCODE, ENCODE, DHS). Columns are
resolved BY NAME from the highest_CFD projection; a missing source degrades to
``-``. The complete raw 85-column dump stays as integrated_results.tsv.gz.

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
2. KEY GRAPHICAL REPORT: the CRISPRme-paper ref/alt scatter (site index log-x vs
   score, red REF / blue ALT points sized by allele frequency, ref->alt arrows,
   top site rsID-annotated), in up to FOUR panels: by CFD, by variant effect
   (CFD ALT-REF delta), by CRISTA, and by variant effect (CRISTA ALT-REF delta).
   The two CRISTA panels appear only when CRISTA was computed -> 4 panels with
   CRISTA, 2 without.
3. SIMPLIFIED reference-vs-population plot (v1 single view).
4. RECOMMENDED VALIDATION PANEL: the full threshold table (candidate counts at
   CFD>= {0.5,0.2,0.05}, mm+b<= {1,2,3,4}, variant-created counts) PLUS the
   HYBRID worst-case top-100 panel: HARD-INCLUDE every site with mm+b<=2 OR
   CFD>=0.5, then FILL to 100 by worst-case severity across CFD desc / CRISTA
   desc (if computed) / mm+b asc (no variant quota). An explicit in-report
   methods note (plain-language, real constants) sits under the panel. The panel
   and each non-empty threshold tier are exported (curated columns) and linked.
5. DOWNLOADS: full RAW integrated_results.tsv.gz + curated top1000.tsv +
   panel_top100.tsv + every non-empty per-tier curated TSV.
6. SCROLLABLE TOP-1000 TABLE (by CFD desc, mm+b<=1 excluded) in the CURATED
   columns, including the annotation columns (gene, distance, GENCODE, ENCODE,
   DHS) and CRISTA when computed.
7. ANNOTATION LEGEND: plain-language meaning of every annotation column value
   (GENCODE / DHS / ENCODE SCREEN v4 cCREs / COSMIC Cancer Gene Census).
FOOTER (unnumbered): CRISPRme version + provenance stamp + fixed research-only
   disclaimer.

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
import json
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
REPORT_GENERATOR_VERSION = "2.4"

# --------------------------------------------------------------------------- #
# Recommended-validation-panel thresholds (module-level constants, section 4)
# --------------------------------------------------------------------------- #
CFD_THRESHOLDS = (0.5, 0.2, 0.05)
# CRISTA is on a DIFFERENT scale than CFD -- reusing CFD's cut points made CRISTA>=0.05
# select ~98% of off-targets (a meaningless "shortlist"). These CRISTA-appropriate,
# higher cut points keep the tiers graduated + model-relative. CRISTA remains the
# score to lean on for bulge/gapped sites (where CFD is out of its training domain).
CRISTA_THRESHOLDS = (0.6, 0.4, 0.2)
MMB_THRESHOLDS = (1, 2, 3, 4)
# threshold-table variant-created CFD floor (kept for the full threshold table)
PANEL_VARIANT_CFD_MIN = 0.05

# --------------------------------------------------------------------------- #
# HYBRID worst-case top-100 panel (section 4).
# --------------------------------------------------------------------------- #
# Over the OFF-TARGET set (on-target mm+b==0 excluded), the panel is built in two
# stages:
#   1. HARD-INCLUDE every site that is close by sequence OR high-scoring, i.e.
#      mm+bulges <= PANEL_FLOOR_MMB  OR  CFD >= PANEL_FLOOR_CFD. These are always
#      in the panel even if they exceed the cap (a low-edit-distance or high-CFD
#      site is never dropped from the confirmation panel).
#   2. FILL the remaining slots up to PANEL_CAP by worst-case severity: each site
#      is ranked independently by CFD (desc), CRISTA (desc; only when computed),
#      and mm+bulges (asc, fewer = closer = worse); a site's SEVERITY is the BEST
#      (minimum) rank across its available metrics, so a site that is worst by ANY
#      single metric floats up. Ties: CFD desc -> CRISTA desc -> mm+b asc.
# There is NO variant-created quota: variant-created sites enter through the same
# floors/ranks as reference sites.
PANEL_CAP = 100
PANEL_FLOOR_MMB = 2  # hard-include every off-target with mm+bulges <= this
PANEL_FLOOR_CFD = 0.5  # hard-include every off-target with CFD >= this
# metrics contributing to the worst-case severity (logical key, direction).
# direction "desc" => higher is worse; "asc" => lower is worse (closer sequence).
PANEL_WORSTCASE_METRICS = (
    ("cfd", "desc"),
    ("crista", "desc"),
    ("mmb", "asc"),
)
# bundled worst-case-panel filename (extra download alongside top1000.tsv)
PANEL_TOP100_NAME = "panel_top100.tsv"

# per-tier download files (section 4/5) are plain .tsv unless the curated TSV
# would exceed this size, in which case they are gzipped (.tsv.gz).
TIER_GZIP_BYTES = 2 * 1024 * 1024  # ~2 MB

# --------------------------------------------------------------------------- #
# ONE curated column set, shared by the in-report top-1000 table AND every
# exported download file (top1000.tsv, panel_top100.tsv, per-tier tsvs).
# --------------------------------------------------------------------------- #
# Each entry is (display_header, kind) where kind selects the value builder in
# ``curated_row`` / ``build_curated_frame``. Columns are resolved BY NAME from
# the integrated_results header (highest_CFD projection); a column whose source
# is missing degrades to "-" rather than being dropped, so every download and the
# table always carry the same, readable, Excel-ready schema. "rank" and "crista"
# are handled specially (rank is 1-based row order; CRISTA appears only when
# crista_computed()). The full 85-column raw dump stays as integrated_results.tsv.gz.
CURATED_COLUMNS = (
    ("rank", "rank"),
    ("Chromosome", "chrom"),
    ("Position", "pos"),
    ("Strand", "strand"),
    ("Aligned_protospacer+PAM", "aligned"),  # ALT with REF fallback
    ("Mismatches", "mm"),
    ("Bulges", "bulges"),
    ("Mismatches+bulges", "mmb"),
    ("Perfect_match", "perfect_match"),  # "Yes" when mm+b==0: a perfect genomic match
    ("CFD", "cfd"),
    ("CRISTA", "crista"),  # emitted only when crista_computed()
    ("REF/ALT_origin", "origin"),
    ("PAM_creation", "pam_creation"),
    ("Variant", "variant"),  # rsID | genomic key when rsID absent
    ("MAF", "maf"),  # em-dash when blank
    ("Gene", "gene_name"),
    ("Gene_distance_kb", "gene_dist"),
    ("GENCODE", "gencode"),
    ("ENCODE", "encode"),
    ("DHS", "dhs"),
    ("COSMIC_cancer_gene", "cosmic"),  # Cancer Gene Census tier/role; "-" when none
    ("High_complexity_region", "complex_region"),  # dense window: greedy shown, more exist
)
# value used when a curated column's source is missing / blank
CURATED_MISSING = "-"

# When True (set via build_report(drop_maf=True) / the --no-maf CLI flag), the MAF
# column is omitted ENTIRELY -- from the curated header, the in-report table, and
# every download file. Use this when a run's allele frequencies are not yet
# finalized (e.g. the Tier-0 AF denominator rebuild is still pending) and showing
# a MAF column would be misleading. Nothing else changes; the raw 85-column dump
# still carries whatever MAF the source TSV had.
_DROP_MAF = False

# Annotation curated kinds that are dropped when their source column is absent from
# the run's TSV, so a run without a given annotation (legacy no-COSMIC bundle, a
# dict-less no-annotation run, a non-human genome) does not show an all-"-" column
# that falsely implies the screen was performed. None => keep all (backward-compat);
# build_report sets it to {kind : kind present in cols} for the run.
_ANNOTATION_KINDS = frozenset(
    {"gencode", "encode", "dhs", "cosmic", "gene_name", "gene_dist"}
)
_PRESENT_ANN_KINDS = None


def _active_columns():
    """``CURATED_COLUMNS`` minus MAF (when ``_DROP_MAF``) and minus any annotation
    column whose source is absent from the run (when ``_PRESENT_ANN_KINDS`` is set)."""
    cols = CURATED_COLUMNS
    if _DROP_MAF:
        cols = tuple(c for c in cols if c[1] != "maf")
    if _PRESENT_ANN_KINDS is not None:
        cols = tuple(
            c for c in cols
            if c[1] not in _ANNOTATION_KINDS or c[1] in _PRESENT_ANN_KINDS
        )
    return cols

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
    "gencode": ["Annotation_GENCODE"],
    "encode": ["Annotation_ENCODE"],
    "dhs": ["Annotation_DHS"],
    "cosmic": ["Annotation_COSMIC"],
    "complex_region": ["High_variant_density_region"],
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


def _scored(num):
    """CFD/CRISTA are defined on [0,1]. Anything else -- the ``do_scores=False``
    sentinel ``-1.000`` that the pipeline writes for non-SpCas9 / non-3nt-PAM /
    5'-PAM nucleases (Cas12a etc.), plus NaN/None -- is 'not scored' -> None."""
    if num is None:
        return None
    try:
        if pd.isna(num):
            return None
        f = float(num)
    except (TypeError, ValueError):
        return None
    return f if 0.0 <= f <= 1.0 else None


def _in_range(series):
    """A float Series masked to [0,1]; out-of-range / unparseable -> NaN, so the
    -1.000 not-scored sentinel never counts, sorts, or thresholds as a real score."""
    s = _to_float_series(series)
    return s.where((s >= 0.0) & (s <= 1.0))


def cfd_computed(df, cols):
    """True when the CFD column exists AND has >=1 in-range [0,1] value (i.e. real
    scores, not the all-(-1) sentinel emitted for out-of-CFD-regime nucleases)."""
    if "cfd" not in cols or cols["cfd"] not in df.columns:
        return False
    return bool(_in_range(df[cols["cfd"]]).notna().any())


def crista_computed(df, cols):
    """True when the CRISTA_score projection column exists AND has >=1 in-range
    [0,1] value. The stable schema always ships the CRISTA columns; a dict-less /
    CRISTA-off / non-SpCas9 run either omits them or fills them with the -1
    sentinel. We include the CRISTA path only when a real score is present.
    """
    if "crista" not in cols or cols["crista"] not in df.columns:
        return False
    return bool(_in_range(df[cols["crista"]]).notna().any())


# --------------------------------------------------------------------------- #
# Curated column projection (ONE set shared by the table + every download file)
# --------------------------------------------------------------------------- #
def _curated_cell(kind, row, cols):
    """Compute the display value for one curated column of one row.

    Values are resolved BY NAME (via ``cols``) from the highest_CFD projection;
    anything missing / blank degrades to ``CURATED_MISSING`` ("-"). ``rank`` and
    ``crista`` are handled by the caller (rank is positional; CRISTA is dropped
    entirely when not computed). Returns a plain string, Excel-ready.
    """
    def _get(key):
        return row.get(cols[key]) if key in cols else None

    if kind == "chrom":
        v = _get("chrom")
    elif kind == "pos":
        v = _get("pos")
    elif kind == "strand":
        v = _get("strand")
    elif kind == "aligned":
        # ALT with REF fallback
        v = row.get(cols["aln_alt"]) if "aln_alt" in cols else None
        if _is_na(v):
            v = row.get(cols["aln_ref"]) if "aln_ref" in cols else None
    elif kind == "mm":
        v = _get("mm")
    elif kind == "bulges":
        v = _get("bulges")
    elif kind == "mmb":
        v = _get("mmb")
    elif kind == "cfd":
        num = _scored(pd.to_numeric(_get("cfd"), errors="coerce")) if "cfd" in cols else None
        return f"{num:.4f}" if num is not None else CURATED_MISSING
    elif kind == "crista":
        num = _scored(pd.to_numeric(_get("crista"), errors="coerce")) if "crista" in cols else None
        return f"{num:.4f}" if num is not None else CURATED_MISSING
    elif kind == "origin":
        v = _get("origin")
        if not _is_na(v):
            return str(v).upper()
        return CURATED_MISSING
    elif kind == "pam_creation":
        v = _get("pam_creation")
    elif kind == "variant":
        # rsID | genomic key when rsID absent
        v = _first_non_na(_get("rsid")) if "rsid" in cols else None
        if v is None and "var_genome" in cols:
            v = _first_non_na(_get("var_genome"))
    elif kind == "maf":
        maf = _min_maf(_get("maf")) if "maf" in cols else None
        # em-dash when blank (MAF footnote explains the blanks)
        return f"{maf:.2e}" if isinstance(maf, float) else CURATED_MISSING
    elif kind == "gene_name":
        v = _get("gene_name")
    elif kind == "gene_dist":
        v = _get("gene_dist")
    elif kind == "gencode":
        v = _get("gencode")
    elif kind == "encode":
        v = _get("encode")
    elif kind == "dhs":
        v = _get("dhs")
    elif kind == "cosmic":
        v = _get("cosmic")
    elif kind == "complex_region":
        # high-variant-density flag: compact "Yes (N var)" for the table; the full
        # note + IUPAC live in integrated_results.tsv and the bundled regions BED.
        v = _get("complex_region")
        if _is_na(v):
            return CURATED_MISSING
        s = str(v)
        if "(" in s and " variants" in s:
            try:
                n = s.split("(", 1)[1].split(" variants")[0].strip()
                return f"Yes ({n} var)"
            except Exception:
                pass
        return "Yes"
    elif kind == "perfect_match":
        # "Yes" for a perfect genomic match (0 mismatches + 0 bulges): a candidate
        # cut site with no a-priori on/off-target distinction. Blank otherwise.
        raw = _get("mmb")
        num = pd.to_numeric(raw, errors="coerce") if raw is not None else None
        return "Yes" if (num is not None and pd.notna(num) and int(num) == 0) else CURATED_MISSING
    else:
        v = None

    if _is_na(v):
        return CURATED_MISSING
    return str(v)


def curated_headers(has_crista):
    """The curated display headers, dropping CRISTA when not computed (and MAF
    when ``_DROP_MAF`` is set)."""
    return [h for h, kind in _active_columns() if kind != "crista" or has_crista]


def build_curated_frame(sub_df, cols, has_crista, start_rank=1):
    """Project a sub-frame onto the ONE curated column set (rows in input order).

    The result is a plain string DataFrame with the curated display headers as
    columns (``rank`` first, CRISTA only when computed), used BOTH to write the
    exported TSVs (top1000/panel/per-tier) and, via ``build_table_html``, the
    in-report table -- so the table and every download share exactly the same
    columns, in the same order, resolved by name (missing -> "-").
    """
    headers = curated_headers(has_crista)
    kinds = [kind for _h, kind in _active_columns() if kind != "crista" or has_crista]
    data = {h: [] for h in headers}
    for offset, (_idx, row) in enumerate(sub_df.iterrows()):
        for h, kind in zip(headers, kinds):
            if kind == "rank":
                data[h].append(str(start_rank + offset))
            else:
                data[h].append(_curated_cell(kind, row, cols))
    return pd.DataFrame(data, columns=headers)


def write_curated_tsv(sub_df, cols, has_crista, path, start_rank=1):
    """Write a sub-frame as a curated-column TSV (shared by every download)."""
    build_curated_frame(sub_df, cols, has_crista, start_rank=start_rank).to_csv(
        path, sep="\t", index=False
    )


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


def dedupe_reference_rows(df, cols):
    """Drop duplicate REFERENCE rows, keeping one per (chrom,pos,strand,target).

    The locus-completeness rule (METHODS.md sec 4) emits the reference off-target
    for every candidate window, so the same reference site can appear more than
    once (once per co-located variant haplotype / per score projection). The
    reference off-target is variant-independent, so those repeats are the SAME
    site and would over-count it in the matrix / plots / table / panel.

    Only REFERENCE rows are de-duplicated, keyed on their site identity
    (chromosome + position + strand + aligned REF protospacer, using whichever of
    those columns are resolvable). VARIANT rows are left untouched -- distinct
    haplotypes at one locus are genuinely different variant off-targets -- and so
    are on-target rows. Distinct reference loci (different chrom/pos/strand/target)
    are preserved. The raw integrated_results dump bundled in the ZIP is a direct
    copy of the source file and is NOT affected. No-op when no identity column is
    resolvable or there are no duplicates. Returns a (possibly) shorter DataFrame.
    """
    _variant, reference, _ontarget = partition_masks(df, cols)
    if not bool(reference.any()):
        return df
    key_keys = [k for k in ("chrom", "pos", "strand", "aln_ref") if k in cols]
    if not key_keys:
        return df  # cannot form a site identity -> leave rows as-is
    key = pd.Series([""] * len(df), index=df.index)
    for k in key_keys:
        key = key.str.cat(df[cols[k]].astype(str), sep="\x1f")
    # duplicate reference rows: reference AND a repeated identity key, keeping first
    dup_ref = reference & key.duplicated(keep="first")
    if not bool(dup_ref.any()):
        return df
    return df[~dup_ref]


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


def _normalize_genome(genome):
    """Display form of a genome token: strip a trailing ``_ref``/``_reference``
    marker some web/CLI configs append to the reference-genome directory name
    (e.g. ``hg38_ref`` -> ``hg38``), so the report header shows the assembly, not
    an internal folder suffix. Genome-agnostic (works for mm10, custom, ...)."""
    g = str(genome or "").strip()
    for suf in ("_reference", "_ref"):
        if g.lower().endswith(suf) and len(g) > len(suf):
            return g[: -len(suf)]
    return g


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
            info["genome"] = _normalize_genome(info["genome"])
            # vcf like "hg38_1000G_HGDP" -> datasets after the ref token
            vcf_tokens = vcf.split("_")
            if vcf_tokens and vcf_tokens[0] == info["genome"]:
                vcf_tokens = vcf_tokens[1:]
            # drop a bare "ref"/"reference" marker (reference-only web/CLI runs
            # name the vcf token "<genome>_ref" -> it must NOT become a dataset)
            vcf_tokens = [
                t for t in vcf_tokens if t and t.lower() not in ("ref", "reference")
            ]
            if vcf_tokens:
                info["datasets"] = "+".join(vcf_tokens)
            # else: leave datasets unset -> build_summary_meta -> "reference only"
        else:
            info["genome"] = _normalize_genome(genome_vcf)
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
    genome = _normalize_genome(
        params.get("Genome_selected")
        or params.get("Genome_ref")
        or fn.get("genome")
        or ""
    ) or "n/a"
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

    # observed max mismatches+bulges on the FEWEST-mm+b (score-neutral) alignment --
    # the basis the off-target matrix uses. A site can still exceed the search cap
    # here only when its minimal alignment against the shown allele does (it is then
    # within budget against the OTHER allele). Falls back to the default mmb column.
    obs_max_mmb = None
    _obs_col = "Mismatches+bulges_(fewest_mm+b)"
    if _obs_col not in df.columns:
        _obs_col = cols["mmb"] if ("mmb" in cols and cols["mmb"] in df.columns) else None
    if _obs_col is not None:
        _mm = _to_int_series(df[_obs_col])
        _mm = _mm[_mm >= 0]
        if len(_mm):
            obs_max_mmb = int(_mm.max())

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
        "obs_max_mmb": obs_max_mmb,
    }


# --------------------------------------------------------------------------- #
# SECTION 1: global off-target-by-MM-and-B matrix (REFERENCE vs VARIANT)
# --------------------------------------------------------------------------- #
def build_mmb_matrix(df, cols, meta):
    """Off-targets binned by origin(REFERENCE/VARIANT) x bulge x mismatch.

    Mirrors the web result-page top matrix: columns Total, 0MM..<mm>MM; rows
    grouped REFERENCE then VARIANT, within each a row per bulge count 0..maxbulges.
    Split by origin only (variant := Not_found_in_REF=="y"); PERFECT matches
    (mm+b==0) are INCLUDED as putative off-targets in the 0MM/0B cell, because a
    guide's intended on-target cannot be told from a perfect-match off-target a
    priori. So REFERENCE + VARIANT matrix totals == the grand total of all sites.

    Extent is chosen so EVERY off-target lands in a cell:
      * bulge rows span 0..(bDNA + bRNA)  -- NOT Max_bulges, which under-counts
        when both a DNA and an RNA bulge co-occur (bDNA=2, bRNA=2 -> rows 0..4);
      * mm columns span 0..Mismatches (0..6).
    If the OBSERVED max ever exceeds the configured span (defensive), the extent
    grows to cover it so no off-target is silently dropped -- guaranteeing the
    matrix REFERENCE + VARIANT totals equal the canonical partition off-target
    totals (REFERENCE + VARIANT + on-target == grand total). Returns a dict with
    ``mm_cols`` (0..mm) and ``groups`` = list of (label, [ (bulge, total, [per-mm
    counts]) ... ]).
    """
    # Origin split ONLY (perfect matches mm+b==0 are kept, not excluded), so the
    # 0MM/0B cell reports the guide's perfect genomic match(es) as putative
    # off-targets alongside everything else.
    if "not_in_ref" in cols:
        variant = df[cols["not_in_ref"]].astype(str).str.strip().str.lower().eq("y")
    elif "origin" in cols:
        variant = df[cols["origin"]].astype(str).str.strip().str.lower().eq("alt")
    else:
        variant = pd.Series([False] * len(df), index=df.index)
    reference = ~variant

    # Place each site by its FEWEST-mismatch+bulge alignment -- the score-NEUTRAL
    # view (it minimizes edits regardless of CFD vs CRISTA), so the matrix does not
    # privilege the higher-edit highest-CFD alignment. Every off-target was found
    # because its minimal alignment is within budget, so almost all sites land in
    # budget here; the few remaining beyond-budget cells are sites within budget
    # against their OTHER allele. Falls back to the default columns if absent.
    _mm_col = "Mismatches_(fewest_mm+b)"
    _b_col = "Bulges_(fewest_mm+b)"
    if _mm_col not in df.columns:
        _mm_col = cols["mm"] if "mm" in cols else None
    if _b_col not in df.columns:
        _b_col = cols["bulges"] if "bulges" in cols else None
    mm_series = _to_int_series(df[_mm_col]) if _mm_col is not None else None
    b_series = _to_int_series(df[_b_col]) if _b_col is not None else None
    if mm_series is None or b_series is None:
        return None

    obs_max_mm = int(mm_series[mm_series >= 0].max()) if (mm_series >= 0).any() else 0
    obs_max_b = int(b_series[b_series >= 0].max()) if (b_series >= 0).any() else 0

    # column extent: run's configured Mismatches, grown to cover observed
    max_mm = meta.get("mm_int")
    max_mm = obs_max_mm if max_mm is None else max(max_mm, obs_max_mm)

    # row extent: bDNA + bRNA (total simultaneous bulges), grown to cover observed.
    # Fall back to Max_bulges, then observed, when DNA/RNA counts are unavailable.
    def _int_or_none(x):
        try:
            return int(str(x).strip())
        except (TypeError, ValueError):
            return None

    bdna = _int_or_none(meta.get("bdna"))
    brna = _int_or_none(meta.get("brna"))
    if bdna is not None and brna is not None:
        max_b = bdna + brna
    elif meta.get("bmax_int") is not None:
        max_b = meta["bmax_int"]
    else:
        max_b = obs_max_b
    max_b = max(max_b, obs_max_b)

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


def render_inputs_criteria(meta, variant_created_name=None, dataset_counts=None,
                           variant_count=None):
    """'Analysis inputs & criteria' box (FDA off-target guidance VII.F.i / F.ii).

    States, in one place, the variant database(s), the variant inclusion policy
    (all common + rare, genic + intergenic, no CRISPRme-applied AF threshold), the
    allele-frequency basis, and the off-target search criteria (variant homology
    gain via fewer mismatches / lowered gaps + variant-created PAM detection).
    Everything is read from the run's meta (no hardcoded database names).
    ``variant_created_name`` is the ACTUAL bundled filename (``.tsv`` or ``.tsv.gz``)
    so the link never breaks; when absent the name is shown as plain text.
    ``dataset_counts`` (``{dataset: n_individuals}`` from the samplesID files) and
    ``variant_count`` (``{'n_records', 'databases'}`` from the Tier-0 registry) add
    the panel size + database SNP-variant count to the inclusion row when known
    (registry sample counts preferred -- they match the AC/AN denominator).
    """
    # bare <code>name</code>; _linkify_bundled_filenames makes every inline mention
    # of a bundled file (this one included) clickable in one uniform pass.
    vc_html = (
        f"<code>{_esc(variant_created_name)}</code>"
        if variant_created_name else "<code>variant_created.tsv</code>"
    )
    ds = meta.get("datasets", "n/a")
    mm = meta.get("mm", "n/a")
    bdna, brna = meta.get("bdna"), meta.get("brna")
    bulges = f"{bdna if bdna is not None else 'n/a'} / {brna if brna is not None else 'n/a'}"
    max_edits = meta.get("max_edits", "n/a")
    rows = [
        ("Variant database(s)", _esc(ds)),
        ("Variants included",
         "All variants present in the database(s) &mdash; common and rare, genic "
         "and intergenic, SNPs and insertions/deletions (indels)."
         + panel_and_variants_note(dataset_counts, variant_count)),
        ("Allele-frequency basis",
         f"Combined-panel minor/alternate allele frequency over the merged "
         f"{_esc(ds)} panel: for genotyped databases this is AC/AN across the "
         f"panel; a frequency-only database contributes its reported population "
         f"allele frequency directly. Per-dataset and per-superpopulation "
         f"frequencies are available in the interactive website."),
        ("Off-target search criteria",
         f"Sites where a variant increases homology to the guide by reducing "
         f"mismatches (up to {_esc(mm)}) and/or lowering gaps (DNA/RNA bulges up to "
         f"{_esc(bulges)}; max total edits {_esc(max_edits)}). Variant-created PAM "
         f"sequences are detected and flagged (<code>PAM_creation</code>). The "
         f"edit cap bounds the search against the variant-enriched (IUPAC-expanded) "
         f"genome; individual reported alignments &mdash; especially "
         f"variant-expanded ones &mdash; may show more total edits."),
        ("Variant-contributed sites",
         "Flagged by REF/ALT origin and <code>PAM_creation</code>; the full list is "
         f"bundled as {vc_html}."),
    ]
    body = "".join(f'<tr><td class="k">{k}</td><td>{v}</td></tr>' for k, v in rows)
    return (
        '<div class="inputs-card" style="margin-top:1em;padding:14px 18px;'
        'background:#f7fafc;border:1px solid #cbd5e0;border-radius:10px">'
        '<div class="matrix-title">Analysis inputs &amp; criteria</div>'
        f'<table class="summary-table"><tbody>{body}</tbody></table>'
        '<p class="caption">Provided to support an off-target analysis accounting '
        "for human genetic variation. The scientific justification for the database "
        "choice, any allele-frequency threshold, and any population stratification "
        "is the sponsor&rsquo;s to provide.</p></div>"
    )


def render_summary_and_matrix(meta, spec_score, matrix):
    """Section 1 HTML: header card (left) + MM/B matrix (right)."""
    bdna = meta["bdna"]
    brna = meta["brna"]

    def _int_or_none(x):
        try:
            return int(str(x).strip())
        except (TypeError, ValueError):
            return None

    _bd, _br = _int_or_none(bdna), _int_or_none(brna)
    if _bd is not None and _br is not None:
        # bDNA / bRNA are PER-TYPE caps; a single alignment may carry both, so the
        # total bulge count reaches bDNA+bRNA (the matrix spans 0..bDNA+bRNA). The
        # old "(max {Max_bulges})" = max(bDNA,bRNA) read as a total cap and
        # contradicted the populated 3B/4B rows.
        bulge_disp = (
            f"{_bd} / {_br} (up to {_bd + _br} total; a DNA and an RNA bulge "
            "may co-occur in one alignment)"
        )
    else:
        bulge_disp = (
            f"{bdna if bdna is not None else 'n/a'} / "
            f"{brna if brna is not None else 'n/a'}"
        )

    # Each SITE is one row carrying both its reference and variant-carrier alignments
    # as columns; the matrix places it ONCE, by its FEWEST-mismatch+bulge alignment
    # (score-neutral). A few land beyond the max-total-edits budget: the reported
    # carrier allele exceeds it while the same row's reference allele stays in scope
    # (the raw search caps mismatches and bulges independently, so an alignment's total
    # can exceed the budget) -- distinct sites, not duplicates, not extra risk; greyed.
    # The full explanation is rendered as a caption BELOW the matrix (_greyed_note).
    _me = meta["max_edits"]
    _obs = meta.get("obs_max_mmb")
    _me_display, _greyed_note = _me, ""
    try:
        if _obs is not None and int(_obs) > int(_me):
            _me_display = f"{_me} (search cap)"
            _cap_bits = []
            if _bd is not None and _br is not None:
                _cap_bits.append(
                    f"up to {meta['mm']} mismatches and {_bd + _br} bulges "
                    f"({_bd} DNA + {_br} RNA)"
                )
            _cap_phrase = (
                f" The raw search caps mismatches and bulges independently &mdash; "
                f"{_cap_bits[0]} &mdash; so an alignment&rsquo;s total can exceed the budget."
                if _cap_bits else ""
            )
            _greyed_note = (
                f" A few sites reach up to {int(_obs)} total edits &mdash; the "
                f"<strong>greyed</strong> cells beyond the max-total-edits budget of "
                f"{_me}. There the reported (variant-carrier) allele exceeds the budget "
                f"while the <em>same row&rsquo;s</em> reference allele stays within scope."
                f"{_cap_phrase} These are the same distinct sites, kept for full locus "
                "coverage &mdash; not duplicates and not additional off-target risk."
            )
    except (TypeError, ValueError):
        pass

    left_rows = [
        ("gRNA (spacer+PAM)", meta["guide_display"]),
        ("Nuclease", meta["nuclease"]),
        ("PAM", meta["pam"]),
        ("Genome", meta["genome"]),
        ("Variant dataset(s)", meta["datasets"]),
        ("Mismatches", meta["mm"]),
        ("Bulges (DNA / RNA)", bulge_disp),
        ("Max total edits", _me_display),
        ("Aggregated Specificity Score (0-100; higher = more specific)", spec_score),
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

        # beyond-budget boundary: cells whose MM+B exceeds the search budget hold
        # reference/variant-reconstructed alignments of within-budget sites (kept for
        # coverage, not extra risk) -- gray them so the in-budget region reads clearly.
        try:
            _budget = int(_me)
        except (TypeError, ValueError):
            _budget = None
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
                for mi, c in enumerate(per_mm):
                    # highlight the 0 MM / 0 B cell -- the perfect match(es)
                    if b == 0 and mm_cols[mi] == 0 and c > 0:
                        cells.append(
                            f"<td style='font-weight:700;background:#fef2f2'>{c:,}</td>"
                        )
                    elif _budget is not None and (b + mm_cols[mi]) > _budget:
                        cells.append(
                            f"<td style='background:#f1f5f9;color:#94a3b8' "
                            f"title='{b + mm_cols[mi]} edits &gt; search budget {_budget}: "
                            f"reference/variant reconstruction, not extra off-target risk'>{c:,}</td>"
                        )
                    else:
                        cells.append(f"<td>{c:,}</td>")
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
    <div class="matrix-title">On-target(s) and Putative Off-targets by Mismatch (MM) and Bulge (B)</div>
    {matrix_html}
    <p class="caption">Each cell counts <strong>distinct putative off-target sites</strong>.
    The highlighted <strong>0&nbsp;MM / 0&nbsp;B</strong> cell holds the guide&rsquo;s
    <strong>perfect genomic match(es)</strong> &mdash; the candidate on-target(s), forced
    to the top of the validation panel below (the intended on-target cannot be told from a
    perfect-match off-target by sequence alone). Rows are split by origin &mdash;
    <strong>REFERENCE</strong> (the target exists in the reference genome, even if a
    variant also alters it) vs <strong>VARIANT</strong> (the target exists only because a
    variant creates it; <code>Not_found_in_REF</code>) &mdash; then by bulge count;
    Total&nbsp;= row sum across mismatches.</p>
    <p class="caption">Every site is <strong>one row</strong> that carries <em>both</em> its
    reference-genome alignment and its variant-carrier alignment as side-by-side columns,
    and is placed here <strong>once</strong> &mdash; by its <strong>fewest-mismatch+bulge</strong>
    alignment, the score-neutral view (it does not prefer CFD over CRISTA). A site is never
    double-counted across cells.{_greyed_note}</p>
  </div>
</div>
"""


# --------------------------------------------------------------------------- #
# samplesID -> superpopulation mapping (for the simplified population plot)
# --------------------------------------------------------------------------- #
def _iter_samplesid_files(samplesid_dir):
    """All samplesID files under EITHER naming convention, deduped and sorted:

      * the finalized install layout ``<genome>_<db...>.samplesID.txt`` (e.g.
        ``hg38_1000G.samplesID.txt``, ``hg38_1000G_HGDP.samplesID.txt``) -- this
        is what a real install actually carries (setup_legacy_database.py renames
        the download to this form); AND
      * the classic/download name ``samplesID*.txt`` (e.g. ``samplesIDs.1000G.txt``).

    Config sidecars (``*.config.*``) are skipped. The old superpop loader only
    matched the second form, so on a real ``hg38_*.samplesID.txt`` install the
    superpopulation plot silently went blank -- this locator fixes that.
    """
    if not samplesid_dir or not os.path.isdir(samplesid_dir):
        return []
    seen, out = set(), []
    for pat in ("*.samplesID.txt", "samplesID*.txt"):
        for path in glob.glob(os.path.join(samplesid_dir, pat)):
            if ".config." in os.path.basename(path):
                continue
            if path not in seen:
                seen.add(path)
                out.append(path)
    return sorted(out)


def _samplesid_dataset_label(basename):
    """``(label, n_tokens)`` for a samplesID filename, or ``None`` if not a per-db
    sample list. ``n_tokens`` is the dataset-token count -- LOWER = more specific
    (a native per-db file), so the caller keeps the most-specific label per sample.
    Handles both ``<genome>_<db...>.samplesID.txt`` and classic ``samplesIDs.<db>.txt``.
    """
    if basename.endswith(".samplesID.txt"):
        base = basename[: -len(".samplesID.txt")]
        parts = base.split("_")
        if len(parts) < 2:
            return None  # need <genome>_<db...>
        return "_".join(parts[1:]), len(parts) - 1  # drop the genome token
    if basename.startswith("samplesID") and basename.endswith(".txt"):
        # classic samplesIDs.<db>.txt -> a single native db
        segs = basename[: -len(".txt")].split(".")
        if len(segs) < 2 or not segs[1]:
            return None
        return ".".join(segs[1:]), 1
    return None


def load_sample_superpop(samplesid_dir, datasets_hint=""):
    """sample_id -> SUPERPOPULATION_ID, loaded from any samplesID files present.

    Mirrors process_summaries.py:38-63 (header
    ``#SAMPLE_ID\\tPOPULATION_ID\\tSUPERPOPULATION_ID\\tSEX``). Returns an empty
    dict when no samplesID files are resolvable; the caller then falls back to
    per-dataset (1000G vs HGDP) provenance.
    """
    mapping = {}
    for path in _iter_samplesid_files(samplesid_dir):
        try:
            with open(path) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue  # header (either #-prefixed or a blank first line)
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3 and parts[0] and parts[0].upper() != "SAMPLE_ID":
                        mapping[parts[0]] = parts[2]
        except OSError:
            continue
    return mapping


def load_sample_dataset(samplesid_dir):
    """sample_id -> native dataset label, READ FROM the per-db samplesID files.

    A merged panel ships a combined ``<genome>_<db1>_<db2>...samplesID.txt`` plus
    per-db ``<genome>_<db>.samplesID.txt`` files. Each sample is assigned the label
    from the MOST SPECIFIC file it appears in (fewest dataset tokens) -- i.e. its
    native per-db provenance. Fully dynamic for ANY database(s); nothing hardcoded.
    The genome token (first ``_`` segment) is stripped to form the label. Returns
    ``{}`` when no samplesID files are resolvable (caller falls back to a heuristic).
    """
    best = {}  # sample_id -> (n_tokens, label)
    for path in _iter_samplesid_files(samplesid_dir):
        parsed = _samplesid_dataset_label(os.path.basename(path))
        if parsed is None:
            continue
        label, n_tokens = parsed
        try:
            with open(path) as handle:
                for line in handle:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    sid = line.split("\t", 1)[0].strip()
                    if not sid or sid.upper() == "SAMPLE_ID":
                        continue
                    prev = best.get(sid)
                    if prev is None or n_tokens < prev[0]:
                        best[sid] = (n_tokens, label)
        except OSError:
            continue
    return {sid: lbl for sid, (_n, lbl) in best.items()}


def _dataset_of(sample_id):
    """Last-resort dataset-provenance heuristic when the per-db samplesID files are
    unavailable (see :func:`load_sample_dataset`, which is preferred)."""
    sid = sample_id.strip()
    if sid.startswith("HGDP"):
        return "HGDP"
    if sid.startswith("HG") or sid.startswith("NA"):
        return "1000G"
    return "other"


def dataset_individual_counts(samplesid_dir):
    """``{dataset_label: n_individuals}`` from the samplesID files -- CHEAP (just
    counts sample IDs; no VCF/genotype parsing). Each sample is counted under its
    MOST-SPECIFIC native dataset (see :func:`load_sample_dataset`), so 1000G and
    HGDP individuals are separated and the combined-panel superset does not double
    count. Returns ``{}`` when no samplesID files are resolvable."""
    ds = load_sample_dataset(samplesid_dir)  # sid -> native dataset label
    counts = {}
    for label in ds.values():
        counts[label] = counts.get(label, 0) + 1
    return counts


def _registry_variant_count(result_dir, meta):
    """Cheap database variant-count + CORRECTED genotyped panel sizes from the
    Tier-0 registry manifests -- index-only (no VCF access, no bcftools, no full
    scan). Returns ``{'n_records': int, 'databases': {label: sample_count}}`` or
    ``None``. Prefers a build-time ``variant_count.json`` sidecar; else sums the
    ``n_records`` key across the per-chrom ``reg_*.idx`` JSON headers (each a small
    manifest). Never raises, never scans VCFs.

    The registry ``sample_count`` is the GENOTYPED panel behind AC/AN (e.g. 3,477
    for 1000G+HGDP), which is authoritative over the samplesID-derived count (that
    can over-list individuals relative to the AC/AN denominator). ``n_records`` is
    the search-usable SNP-variant count (indels / no-carrier records excluded)."""
    try:
        if not result_dir:
            return None
        dict_dir = None
        for cand in (
            os.path.join(result_dir, "..", "..", "Dictionaries"),
            os.path.join(result_dir, "..", "Dictionaries"),
            os.path.join(os.getcwd(), "Dictionaries"),
        ):
            if os.path.isdir(cand):
                dict_dir = cand
                break
        if dict_dir is None:
            return None
        reg_dirs = [d for d in sorted(glob.glob(os.path.join(dict_dir, "registry_*")))
                    if os.path.isdir(d)]
        if not reg_dirs:
            return None
        # pick the registry whose name carries the run's dataset tokens
        tokens = [t for t in str(meta.get("datasets", "")).replace("+", "_").split("_") if t]
        chosen = None
        for d in reg_dirs:
            name = os.path.basename(d).lower()
            if tokens and all(t.lower() in name for t in tokens):
                chosen = d
                break
        if chosen is None:
            chosen = reg_dirs[0] if len(reg_dirs) == 1 else None
        if chosen is None:
            return None

        def _dbs(manifest):
            return {
                k: v.get("sample_count")
                for k, v in (manifest.get("databases") or {}).items()
                if isinstance(v, dict) and v.get("sample_count")
            }

        sidecar = os.path.join(chosen, "variant_count.json")
        if os.path.isfile(sidecar):
            with open(sidecar) as fh:
                m = json.load(fh)
            n = int(m.get("n_records", 0) or 0)
            if n > 0:
                _ni = m.get("n_indels")
                return {
                    "n_records": n,
                    "n_indels": int(_ni) if _ni else None,
                    "databases": _dbs(m),
                }
        # fallback: sum n_records across the per-chrom reg_*.idx headers. The idx
        # files carry no indel count, so n_indels is only available from the
        # build-time sidecar above.
        total, dbs = 0, {}
        for p in sorted(glob.glob(os.path.join(chosen, "reg_*.idx"))):
            try:
                with open(p) as fh:
                    m = json.load(fh)
            except (OSError, ValueError):
                continue
            total += int(m.get("n_records", 0) or 0)
            if not dbs:
                dbs = _dbs(m)
        return ({"n_records": total, "n_indels": None, "databases": dbs}
                if total > 0 else None)
    except Exception:  # noqa: BLE001 - a count is optional; never break the report
        return None


def panel_and_variants_note(dataset_counts, variant_count=None):
    """Panel size + database variant count for the 'Variants included' row.

    Prefers the registry's GENOTYPED ``sample_count`` (``variant_count['databases']``)
    over the samplesID-derived ``dataset_counts`` -- the latter can over-list
    individuals relative to the AC/AN denominator. Appends the search-usable SNP
    variant count when known. Returns "" when nothing is available."""
    counts = None
    if variant_count and variant_count.get("databases"):
        counts = {k: v for k, v in variant_count["databases"].items() if v}
    if not counts:
        counts = dataset_counts or {}
    if not counts:
        return ""
    total = sum(counts.values())
    if total <= 0:
        return ""
    parts = ", ".join(f"{lbl} n={n:,}" for lbl, n in sorted(counts.items()))
    per_ds = f" ({parts})" if len(counts) > 1 else ""
    lead = f" Panel: <strong>{total:,}</strong> individuals{per_ds}."
    n_snp = variant_count.get("n_records") if variant_count else None
    n_indel = variant_count.get("n_indels") if variant_count else None
    if n_snp:
        # SNPs come from the registry; indels from the separate indel index. Show
        # BOTH counts as searched when the indel count is KNOWN (build-time
        # manifest) -- including a genuine 0 for a SNP-only database; when the
        # indel count is UNKNOWN (n_indel is None: legacy install / idx fallback)
        # state indels are also searched without inventing a number. Fully
        # database-agnostic -- nothing here assumes a particular panel.
        if n_indel is not None:
            lead += (
                f" The database contributes <strong>{n_snp:,}</strong> SNPs and "
                f"<strong>{n_indel:,}</strong> indels, all searched."
            )
        else:
            lead += (
                f" The database contributes <strong>{n_snp:,}</strong> SNPs; "
                "insertions/deletions from the same panel are also searched."
            )
    return lead


# --------------------------------------------------------------------------- #
# Plots (matplotlib -> in-memory figure -> base64 data URI)
# --------------------------------------------------------------------------- #
# SVG (vector) keeps every figure crisp at any zoom / print resolution -- IND-grade,
# publication-quality -- and is the SAME matplotlib drawing as the PNG, just not
# rasterized. Text is emitted as vector OUTLINES (svg.fonttype='path'), i.e. the
# glyphs are embedded as paths taken from the exact font matplotlib used, so the
# figure renders byte-identically in any browser/viewer with NO font dependency
# (no fallback-font substitution). Flip _FIG_FORMAT to "png" to fall back to raster.
_FIG_FORMAT = "svg"
matplotlib.rcParams["svg.fonttype"] = "path"  # embed text as vector outlines


def _fig_to_data_uri(fig, dpi=120):
    buf = io.BytesIO()
    if _FIG_FORMAT == "svg":
        fig.savefig(buf, format="svg", bbox_inches="tight")
        plt.close(fig)
        encoded = base64.b64encode(buf.getvalue()).decode()
        return f"data:image/svg+xml;base64,{encoded}"
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


def _cfd_style_scatter(
    sub, score_key, ref_key, alt_key, xlabel, title,
    score_name="CFD", maf_col=None, samp_col=None, rsid_col=None,
):
    """The CRISPRme-paper CFD/CRISTA scatter, adapted from
    CRISPRme_plots.py:plot_with_CFD_score (~L347).

    Generalized over the score family: ``score_key`` / ``ref_key`` / ``alt_key``
    are the resolved combined/REF/ALT score-column names (CFD or CRISTA), and
    ``maf_col`` / ``samp_col`` / ``rsid_col`` are the resolved allele-frequency /
    samples / rsID column names for that same projection. ``score_name`` is the
    y-axis / annotation label ("CFD" or "CRISTA").

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
    if maf_col and maf_col in work.columns:
        af = work[maf_col].astype(str).str.split(",").apply(
            lambda toks: min(
                [float(t) for t in toks if _num_ok(t)] or [-1.0]
            )
        )
    else:
        af = pd.Series([-1.0] * n)
    work["AF"] = pd.to_numeric(af, errors="coerce").fillna(-1.0)

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
    # marker sizes must be finite: a reference off-target (no AF, no variant) gets
    # NaN above, which would crash fig.savefig (reference-only runs). NaN -> a small
    # visible marker.
    work["plot_AF"] = np.nan_to_num(np.asarray(work["plot_AF"], dtype=float), nan=6.0)
    work["ref_AF"] = np.sqrt(np.clip(1 - work["AF"], 0, None) * 1000)

    y_ref = _to_float_series(work[ref_key]) if ref_key in work.columns else None
    y_alt = _to_float_series(work[alt_key]) if alt_key in work.columns else None
    # fall back to the combined score column if REF/ALT variants are missing
    y_combined = _to_float_series(work[score_key]) if score_key in work.columns else None
    if y_ref is None:
        y_ref = y_combined
    if y_alt is None:
        y_alt = y_combined
    # nothing scoreable in [0,1] (e.g. the -1 sentinel for a non-CFD-regime nuclease)
    # -> raise so the caller renders a labeled placeholder, not an empty off-axis plot.
    _yr = y_ref if y_ref is not None else pd.Series([], dtype=float)
    _ya = y_alt if y_alt is not None else pd.Series([], dtype=float)
    if not (_in_range(_yr).notna().any() or _in_range(_ya).notna().any()):
        raise ValueError(f"no valid {score_name} score to plot")

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
    ax.set_ylabel(f"{score_name} score")
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


def plot_scatter_panels(df, cols, n=1000, include_crista=False):
    """Produce the key graphical-report scatter panels (section 2).

    Returns a list of (title, caption, data_uri). Panels:
      (a) top-N by CFD score
      (b) top-N by CFD DELTA = ALT CFD - REF CFD (largest variant increase first)
      (c) top-N by CRISTA score               -- only when include_crista is True
      (d) top-N by CRISTA DELTA = ALT - REF    -- only when include_crista is True
    So: CRISTA computed -> 4 panels; CRISTA absent -> 2 panels. Every panel is
    guarded; a failing panel yields a placeholder URI (never aborts).
    """
    panels = []

    # candidate frame: drop on-/near-on-target rows (mm+b<=1), same as the table
    base = df
    if "mmb" in cols:
        base = base[_to_int_series(base[cols["mmb"]]) > 1]
    n_shown = min(n, len(base))

    def _score_sorted(score_col):
        s = _to_float_series(base[score_col]).fillna(-1.0)
        return base.assign(_s=s).sort_values("_s", ascending=False).head(n)

    def _delta_sorted(alt_col, ref_col):
        delta = _to_float_series(base[alt_col]) - _to_float_series(base[ref_col])
        return (
            base.assign(_delta=delta)
            .sort_values("_delta", ascending=False, na_position="last")
            .head(n)
        )

    # (a) by CFD score
    cfd_score = cols.get("cfd", "")
    cfd_ref = cols.get("cfd_ref", "")
    cfd_alt = cols.get("cfd_alt", "")
    try:
        uri = _cfd_style_scatter(
            _score_sorted(cfd_score), cfd_score, cfd_ref, cfd_alt,
            xlabel="Candidate off-target site (ranked by CFD)",
            title=f"Top {n_shown} candidates by CFD score",
            score_name="CFD",
            maf_col=cols.get("maf"), samp_col=cols.get("samples"),
            rsid_col=cols.get("rsid"),
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

    # (b) by CFD DELTA -- only when the REF/ALT CFD columns exist
    if "cfd_ref" in cols and "cfd_alt" in cols:
        try:
            uri = _cfd_style_scatter(
                _delta_sorted(cfd_alt, cfd_ref), cfd_score, cfd_ref, cfd_alt,
                xlabel="Candidate off-target site (ranked by variant effect ALT-REF)",
                title=f"Top {n_shown} by variant-induced CFD increase (ALT-REF)",
                score_name="CFD",
                maf_col=cols.get("maf"), samp_col=cols.get("samples"),
                rsid_col=cols.get("rsid"),
            )
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"generate-report: CFD delta scatter unavailable: {exc}\n")
            uri = _placeholder_uri("CFD delta scatter unavailable")
        panels.append((
            "By variant effect (CFD ALT - REF delta)",
            "Same axes, re-ranked by the variant-induced CFD change (ALT-REF, "
            "descending): the variants that most raise the CFD score come first, "
            "foregrounding the risk-relevant variant-created sites.",
            uri,
        ))

    # (c) + (d) CRISTA, only when computed
    if include_crista:
        cr_score = cols.get("crista", "")
        cr_ref = cols.get("crista_ref", "")
        cr_alt = cols.get("crista_alt", "")
        cr_maf = cols.get("crista_maf")
        cr_samp = cols.get("crista_samples")
        cr_rsid = cols.get("crista_rsid")

        # (c) by CRISTA score
        try:
            uri = _cfd_style_scatter(
                _score_sorted(cr_score), cr_score, cr_ref, cr_alt,
                xlabel="Candidate off-target site (ranked by CRISTA)",
                title=f"Top {n_shown} candidates by CRISTA score",
                score_name="CRISTA",
                maf_col=cr_maf, samp_col=cr_samp, rsid_col=cr_rsid,
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

        # (d) by CRISTA DELTA -- only when the REF/ALT CRISTA columns exist
        if "crista_ref" in cols and "crista_alt" in cols:
            try:
                uri = _cfd_style_scatter(
                    _delta_sorted(cr_alt, cr_ref), cr_score, cr_ref, cr_alt,
                    xlabel="Candidate off-target site (ranked by variant effect ALT-REF)",
                    title=f"Top {n_shown} by variant-induced CRISTA increase (ALT-REF)",
                    score_name="CRISTA",
                    maf_col=cr_maf, samp_col=cr_samp, rsid_col=cr_rsid,
                )
            except Exception as exc:  # noqa: BLE001
                sys.stderr.write(
                    f"generate-report: CRISTA delta scatter unavailable: {exc}\n"
                )
                uri = _placeholder_uri("CRISTA delta scatter unavailable")
            panels.append((
                "By variant effect (CRISTA)",
                "The same ref/alt scatter as for the CFD score above, now with "
                "CRISTA on the y-axis, "
                "re-ranked by the variant-induced CRISTA change (ALT-REF, "
                "descending): the population variants that most raise the "
                "independent CRISTA cleavage score come first. Included because "
                "CRISTA scores were computed for this run.",
                uri,
            ))

    return panels


def plot_population(df, cols, sample_superpop, sample_dataset=None):
    """SECTION 3: SIMPLIFIED population view at the run's SELECTED parameters.

    Two panels over the FULL result set (no per-total-mm faceting):
      left  : reference vs variant-created off-target counts
      right : among variant-created, one bar per superpopulation (or per
              dataset when superpop mapping is unavailable).
    """
    variant, reference, _ontarget = partition_masks(df, cols)
    variant_mask = variant
    n_variant = int(variant_mask.sum())
    # OFF-TARGETS only: exclude the on-target (mm+b==0) from the reference bar --
    # the chart counts off-targets, and the on-target is not one (canonical
    # partition: reference + variant + on-target == total).
    n_reference = int(reference.sum())

    group_counts = {}
    use_superpop = bool(sample_superpop)
    _ds = sample_dataset or {}  # native per-db provenance (dynamic; see load_sample_dataset)
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
                    # keep the superpopulation axis pure: a sample with no
                    # superpop mapping goes to "unknown", NOT a dataset name (which
                    # would read as if a dataset were a superpopulation). With the
                    # samplesID-glob fix this fallback is rarely hit.
                    grp = sample_superpop.get(sid) or "unknown"
                else:
                    grp = _ds.get(sid) or _dataset_of(sid)  # native provenance
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
def select_worstcase_panel(df, cols, cap=PANEL_CAP):
    """Return the rows for the HYBRID worst-case top-100 validation panel (sec 4).

    Every PERFECT match (mm+b==0) is placed at the TOP of the panel first -- a
    candidate cut site with no a-priori on/off-target distinction, never dropped.
    The OFF-TARGET rows are then selected in two stages (over the off-target set):

    1. HARD-INCLUDE every site that is close by sequence OR high-scoring:
       ``mm+bulges <= PANEL_FLOOR_MMB (2)`` OR ``CFD >= PANEL_FLOOR_CFD (0.5)``.
       These are always kept; if the hard-includes already exceed ``cap`` we keep
       them all (a low-edit-distance / high-CFD site is never dropped).
    2. FILL the remaining slots up to ``cap`` by worst-case severity. Each site
       is ranked independently by every available metric in
       ``PANEL_WORSTCASE_METRICS``: CFD (desc), CRISTA (desc; only when computed)
       and mm+bulges (asc, fewer = closer sequence = worse); rank 1 == worst. A
       site's SEVERITY is the BEST (minimum) rank across its available metrics,
       so a site that is worst by ANY single metric floats up. Fill sites are
       taken by ascending severity; ties are broken by CFD desc -> CRISTA desc ->
       mm+b asc.

    There is NO variant-created quota: variant-created sites qualify through the
    same floors/ranks as reference sites. Returns the selected sub-frame (original
    columns): hard-includes first (severity-ordered), then the fill (also
    severity-ordered).
    """
    _variant, _reference, ontarget = partition_masks(df, cols)
    # Perfect matches (0 mismatch + 0 bulge) are candidate cut sites with no
    # a-priori on/off-target distinction -- ALWAYS include every one, at the TOP of
    # the panel, never subject to the cap (a perfect-match off-target cuts as
    # efficiently as the intended site and must be validated).
    perfect = df[ontarget].copy()
    offt = df[~ontarget].copy()
    if len(offt) == 0:
        return perfect

    cfd = (
        _to_float_series(offt[cols["cfd"]]).fillna(-1.0)
        if "cfd" in cols else pd.Series(-1.0, index=offt.index)
    )
    crista = (
        _to_float_series(offt[cols["crista"]])
        if "crista" in cols else pd.Series(np.nan, index=offt.index)
    )
    mmb = (
        _to_int_series(offt[cols["mmb"]])
        if "mmb" in cols else pd.Series(10 ** 6, index=offt.index)
    )
    # unparseable mmb -> -1 (from _to_int_series); clamp to a large sentinel so it is
    # neither hard-included (mmb <= floor) nor ranked most-severe (mmb asc rank=1)
    mmb = mmb.where(mmb >= 0, 10 ** 6)
    has_crista = ((crista >= 0) & (crista <= 1)).any()

    # per-metric ranks (rank 1 == worst). ascending flag flips per direction:
    #   desc metric -> higher is worse -> rank ascending=False
    #   asc  metric -> lower  is worse -> rank ascending=True
    series_by_key = {"cfd": cfd, "crista": crista, "mmb": mmb}
    rank_frames = []
    for key, direction in PANEL_WORSTCASE_METRICS:
        if key == "crista" and not has_crista:
            continue  # CRISTA absent -> this metric does not contribute
        s = series_by_key.get(key)
        if s is None:
            continue
        ascending = direction == "asc"
        rank_frames.append(s.rank(method="min", ascending=ascending))

    if rank_frames:
        # severity = BEST (minimum) available rank across the contributing metrics
        severity = pd.concat(rank_frames, axis=1).min(axis=1)
    else:
        severity = pd.Series(1.0, index=offt.index)

    # STAGE 1: hard-includes (mm+b <= floor OR CFD >= floor)
    hard_mask = (mmb <= PANEL_FLOOR_MMB) | (cfd >= PANEL_FLOOR_CFD)

    ordered = offt.assign(
        _severity=severity, _cfd=cfd, _crista=crista.fillna(-1.0), _mmb=mmb,
        _hard=hard_mask,
    ).sort_values(
        # hard-includes first, then by worst-case severity; ties CFD/CRISTA/mm+b
        ["_hard", "_severity", "_cfd", "_crista", "_mmb"],
        ascending=[False, True, False, False, True],
    )

    n_hard = int(hard_mask.sum())
    # if the hard-includes already exceed the cap keep them ALL; otherwise fill
    keep = max(cap, n_hard)
    drop = ["_severity", "_cfd", "_crista", "_mmb", "_hard"]
    off_panel = ordered.head(keep).drop(columns=drop)
    # prepend the perfect matches (disjoint from off_panel by construction)
    if len(perfect):
        return pd.concat([perfect, off_panel])
    return off_panel


def build_validation_panel(df, cols):
    """Compute candidate counts at CFD / edit-distance thresholds and the
    worst-case suggested panel. Returns a dict for rendering.

    Counting is over the OFF-TARGET set (on-target row mm+b==0 excluded), so
    thresholds report actionable candidates for a confirmation panel. The FULL
    threshold table is kept (CFD>= {0.5,0.2,0.05}, mm+b<= {1,2,3,4}, and
    variant-created counts); the suggested tier is the worst-case top-100.
    """
    variant, _reference, ontarget = partition_masks(df, cols)
    offt = df[~ontarget]  # off-targets only

    # Perfect matches (0 mm + 0 bulge): distinct candidate cut sites. With >1 there
    # is no a-priori on-target, so these drive the red warning banner + panel.
    perfect_sites = []
    perfect_by_guide = {}
    if bool(ontarget.any()):
        seen = set()
        gcol = cols.get("guide")
        for _, r in df[ontarget].iterrows():
            c = str(r[cols["chrom"]]) if "chrom" in cols else "?"
            p = str(r[cols["pos"]]) if "pos" in cols else "?"
            s = str(r[cols["strand"]]) if "strand" in cols else ""
            g = str(r[gcol]) if gcol and gcol in df.columns else ""
            # origin (ref vs variant-created) + MAF distinguish a UNIVERSAL reference
            # perfect match from a RARE variant-created one in the banner.
            o = (str(r[cols["origin"]]).strip().lower()
                 if "origin" in cols and cols["origin"] in df.columns else "")
            m = (_min_maf(r[cols["maf"]])
                 if "maf" in cols and cols["maf"] in df.columns else None)
            # dedup per (guide, locus): in a MULTI-guide run each guide has its own
            # on-target -- ambiguity is per-guide, not across guides.
            if (g, c, p, s) not in seen:
                seen.add((g, c, p, s))
                perfect_sites.append(
                    {"guide": g, "chrom": c, "pos": p, "strand": s,
                     "origin": o, "maf": m}
                )
                perfect_by_guide[g] = perfect_by_guide.get(g, 0) + 1

    # INDEX-ALIGNED sentinel fallbacks (mirror select_worstcase_panel): a length-0
    # Series here would make (cfd>=t)/(mmb<=t) length-0 and raise "wrong length"
    # (ValueError) that silently wipes the entire Section-4 validation panel for a
    # CRISTA-only / mm+b-only / CFD-skipped run. Full-length sentinels -> empty tiers.
    cfd = _to_float_series(offt[cols["cfd"]]) if "cfd" in cols else pd.Series(-1.0, index=offt.index)
    mmb = _to_int_series(offt[cols["mmb"]]) if "mmb" in cols else pd.Series(10 ** 6, index=offt.index)
    crista = (
        _to_float_series(offt[cols["crista"]]) if "crista" in cols
        else pd.Series(np.nan, index=offt.index)
    )
    var_off = variant[~ontarget]

    cfd_counts = [(t, int((cfd >= t).sum())) for t in CFD_THRESHOLDS]
    mmb_counts = [(t, int(((mmb >= 0) & (mmb <= t)).sum())) for t in MMB_THRESHOLDS]
    crista_counts = (
        [(t, int((crista >= t).sum())) for t in CRISTA_THRESHOLDS]
        if "crista" in cols and ((crista >= 0) & (crista <= 1)).any() else []
    )

    n_variant = int(var_off.sum())
    n_variant_cfd = int((var_off & (cfd >= PANEL_VARIANT_CFD_MIN)).sum())

    has_crista = crista_computed(df, cols)

    # hybrid worst-case top-100 panel
    panel_df = select_worstcase_panel(df, cols, cap=PANEL_CAP)
    panel_size = len(panel_df)
    # how many of the selected panel are variant-created
    if "not_in_ref" in cols and cols["not_in_ref"] in panel_df.columns:
        panel_variant = int(
            panel_df[cols["not_in_ref"]].astype(str).str.strip().str.lower().eq("y").sum()
        )
    else:
        panel_variant = 0

    # per-tier off-target subsets (for the bundled curated downloads + links).
    # Each entry: (logical tier key, display label, sub-frame). Only non-empty
    # tiers are bundled/linked (decided by the caller).
    tiers = build_tier_frames(df, cols, offt, variant, ontarget, cfd, mmb, crista)

    return {
        "cfd_counts": cfd_counts,
        "crista_counts": crista_counts,
        "mmb_counts": mmb_counts,
        "n_variant": n_variant,
        "n_variant_cfd": n_variant_cfd,
        "n_offtarget": int((~ontarget).sum()),
        "has_crista": has_crista,
        "panel_size": panel_size,
        "panel_variant": panel_variant,
        "panel_df": panel_df,
        "tiers": tiers,
        "n_perfect": len(perfect_sites),
        "perfect_sites": perfect_sites,
        # per-guide multiplicity: the red "no unambiguous on-target" banner fires
        # only when a SINGLE guide has >= 2 perfect matches (a multi-guide run with
        # one perfect match each is normal, not ambiguous).
        "max_perfect_per_guide": (
            max(perfect_by_guide.values()) if perfect_by_guide else 0
        ),
        "n_guides_with_perfect": len(perfect_by_guide),
    }


def _tier_filename(key):
    """Canonical bundled filename for a per-tier subset (curated TSV)."""
    if key == "panel":
        return PANEL_TOP100_NAME
    if key.startswith("cfd_"):
        return f"cfd_ge_{key.split('_', 1)[1]}.tsv"
    if key.startswith("crista_"):
        return f"crista_ge_{key.split('_', 1)[1]}.tsv"
    if key.startswith("mmb_"):
        return f"mmb_le_{key.split('_', 1)[1]}.tsv"
    return f"{key}.tsv"


def build_tier_frames(df, cols, offt, variant, ontarget, cfd, mmb, crista=None):
    """Off-target subsets for each threshold tier (section-4 downloads).

    Returns a list of dicts ``{key, label, filename, df}`` for the CFD>= and
    mm+b<= tiers and the variant-created tier, over the OFF-TARGET set (on-target
    mm+b==0 excluded). ``cfd`` / ``mmb`` are the OFF-TARGET-aligned series already
    computed by the caller. Empty tiers are still returned; the caller skips
    bundling/linking empties. Sub-frames carry the ORIGINAL columns so the
    curated projection is applied uniformly at write time.
    """
    tiers = []
    var_off = variant[~ontarget]
    for t in CFD_THRESHOLDS:
        sub = offt[(cfd >= t).values]
        tiers.append({
            "key": f"cfd_{t:.2f}",
            "label": f"CFD &ge; {t}",
            "filename": _tier_filename(f"cfd_{t:.2f}"),
            "df": sub,
        })
    if crista is not None and "crista" in cols and ((crista >= 0) & (crista <= 1)).any():
        for t in CRISTA_THRESHOLDS:
            sub = offt[(crista >= t).values]
            tiers.append({
                "key": f"crista_{t:.2f}",
                "label": f"CRISTA &ge; {t}",
                "filename": _tier_filename(f"crista_{t:.2f}"),
                "df": sub,
            })
    for t in MMB_THRESHOLDS:
        # mask invalid mmb (_to_int_series maps unparseable -> -1; -1 <= t would
        # falsely pull a garbage site into EVERY low-edit tier + the panel)
        sub = offt[((mmb >= 0) & (mmb <= t)).values]
        tiers.append({
            "key": f"mmb_{t}",
            "label": f"mismatches + bulges &le; {t}",
            "filename": _tier_filename(f"mmb_{t}"),
            "df": sub,
        })
    tiers.append({
        "key": "variant_created",
        "label": "variant-contributed off-target sites (variant-created)",
        "filename": "variant_created.tsv",
        "df": offt[var_off.values],
    })
    return tiers


def render_validation_panel(vp, panel_tsv_name=None, tier_links=None):
    """Section 4 HTML: full threshold table + hybrid panel + methods note + the
    per-tier download table.

    ``panel_tsv_name`` (when given) is the bundled hybrid-panel filename.
    ``tier_links`` maps a tier key ("panel", "cfd_0.50", ..., "mmb_1", ...,
    "variant_created") to its bundled filename, so a Download column / link
    appears only for the tiers that were actually bundled (non-empty ones).
    """
    tier_links = tier_links or {}

    def _tier_count_link(key, count):
        fname = tier_links.get(key)
        if not fname:
            # no bundled file: an empty tier (0) reads as a blank/missing value if
            # shown as a bare "0" -- label it as a genuine count of zero sites
            return "0 <span class='caption'>(none)</span>" if count == 0 else f"{count:,}"
        return (
            f"{count:,} &nbsp;<a class=\"tier-dl\" href=\"{_dl_href(fname)}\" "
            f"download>{_esc(fname)}</a>"
        )

    cfd_rows = "".join(
        f"<tr><td>CFD &ge; {t}</td><td class='num'>"
        f"{_tier_count_link(f'cfd_{t:.2f}', c)}</td></tr>"
        for t, c in vp["cfd_counts"]
    )
    mmb_rows = "".join(
        f"<tr><td>mismatches + bulges &le; {t}</td><td class='num'>"
        f"{_tier_count_link(f'mmb_{t}', c)}</td></tr>"
        for t, c in vp["mmb_counts"]
    )
    # CRISTA thresholds, symmetric with CFD (only when CRISTA was computed)
    crista_table_block = ""
    if vp.get("crista_counts"):
        crista_rows = "".join(
            f"<tr><td>CRISTA &ge; {t}</td><td class='num'>"
            f"{_tier_count_link(f'crista_{t:.2f}', c)}</td></tr>"
            for t, c in vp["crista_counts"]
        )
        crista_table_block = (
            '<div><table class="thr-table"><thead><tr>'
            "<th>CRISTA threshold</th><th>Candidates (download)</th></tr></thead>"
            f"<tbody>{crista_rows}</tbody></table></div>"
        )
    metric_names = ["CFD (desc)"]
    if vp.get("has_crista"):
        metric_names.append("CRISTA (desc)")
    metric_names.append("mismatches+bulges (asc)")
    metric_list = ", ".join(metric_names)

    # C) EXPLICIT IN-REPORT METHODS NOTE (plain-language, using the real
    #    constants). Two hard-include floors, then fill by worst-case severity.
    metric_or = (
        "CFD, CRISTA, or mm+b" if vp.get("has_crista") else "CFD or mm+b"
    )
    note = (
        f"How the panel was chosen (hybrid, ~{PANEL_CAP} sites &mdash; may be more "
        f"when many sites are hard-included). "
        f"First, every off-target that is CLOSE by sequence OR HIGH-scoring is "
        f"hard-included &mdash; specifically every site with mismatches+bulges "
        f"&le; {PANEL_FLOOR_MMB} OR CFD &ge; {PANEL_FLOOR_CFD}. These are always "
        f"kept (if the hard-included sites already exceed {PANEL_CAP}, they are "
        f"all kept). The remaining slots up to {PANEL_CAP} are then filled by "
        f"worst-case severity: each site is ranked independently by "
        f"{metric_list}, and a site is prioritized if it is worst by ANY single "
        f"one of those metrics ({metric_or}) &mdash; so the highest-scoring "
        f"predicted cleavage sites AND the near-cognate low-edit-distance "
        f"sequences that scoring models can under-weight both surface. There is "
        f"NO category quota: variant-created sites qualify through the same "
        f"floors and ranks as reference sites. The full threshold table above, "
        f"the per-threshold subset files, and the complete raw integrated results "
        f"(bundled downloads) let the panel be reviewed or expanded for any assay "
        f"budget."
    )
    _panel_n = vp.get("panel_size", 0)
    _panel_dl_n = f"~{PANEL_CAP} sites" + (
        f", {_panel_n:,} here" if _panel_n > PANEL_CAP else ""
    )
    panel_dl = ""
    if panel_tsv_name:
        panel_dl = (
            f'<p><a class="download" href="{_dl_href(panel_tsv_name)}" download>'
            f"Recommended hybrid worst-case panel (TSV, {_panel_dl_n})"
            f"</a></p>"
        )
    # off-target-count bridge: show the total-sites minus perfect-match subtraction
    # so summing the Section-1 matrix (which includes perfect matches) reconciles.
    _np = int(vp.get("n_perfect", 0) or 0)
    _off_label = "Off-targets (on-target mm+b=0 excluded)"
    if _np > 0:
        _total_sites = vp["n_offtarget"] + _np
        _off_label = (
            f"Off-targets ({_total_sites:,} total sites &minus; {_np:,} perfect "
            f"match{'es' if _np != 1 else ''} = {vp['n_offtarget']:,})"
        )
    # panel can EXCEED PANEL_CAP because all hard-included sites are always kept
    _panel_hdr = f"Recommended panel &mdash; hybrid worst-case (~{PANEL_CAP} sites"
    if vp.get("panel_size", 0) > PANEL_CAP:
        _panel_hdr += f"; {vp['panel_size']:,} here &mdash; all hard-included sites kept"
    _panel_hdr += ")"
    return f"""
<div class="panel-grid">
  <div>
    <table class="thr-table"><thead><tr><th>CFD threshold</th><th>Candidates (download)</th></tr></thead>
    <tbody>{cfd_rows}</tbody></table>
  </div>
  {crista_table_block}
  <div>
    <table class="thr-table"><thead><tr><th>Edit-distance threshold</th><th>Candidates (download)</th></tr></thead>
    <tbody>{mmb_rows}</tbody></table>
  </div>
</div>
<table class="thr-table" style="max-width:720px">
  <tbody>
    <tr><td>{_off_label}</td><td class="num">{vp['n_offtarget']:,}</td></tr>
    <tr><td>Variant-created off-targets (Not_found_in_REF)</td><td class="num">{_tier_count_link('variant_created', vp['n_variant'])}</td></tr>
    <tr><td>&hellip; of those with CFD &ge; {PANEL_VARIANT_CFD_MIN}</td><td class="num">{vp['n_variant_cfd']:,}</td></tr>
    <tr class="panel-hi"><td><strong>{_panel_hdr}</strong><br>
      <span class="caption">hard-include (mm+b &le; {PANEL_FLOOR_MMB} OR CFD &ge;
      {PANEL_FLOOR_CFD}), then fill up to {PANEL_CAP} by worst-case severity across
      {metric_list}; of the selected sites, <strong>{vp['panel_variant']:,}</strong>
      are variant-created.</span></td>
      <td class="num"><strong>{_tier_count_link('panel', vp['panel_size'])}</strong></td></tr>
  </tbody>
</table>
{panel_dl}
<p class="caption"><strong>Methods.</strong> {note}</p>
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


def select_top_crista(df, cols, n=1000):
    """Top-N off-targets by CRISTA desc (mm+b > 1), when CRISTA is computed.

    Mirrors :func:`select_top` (same on-/near-on-target filter, same curated
    columns) but ranks by the CRISTA score instead of CFD. Returns an EMPTY frame
    when CRISTA is absent, so the caller renders the CRISTA table only when it
    exists.
    """
    if "crista" not in cols:
        return df.iloc[0:0]
    work = df.copy()
    if "mmb" in cols:
        work = work[_to_int_series(work[cols["mmb"]]) > 1]
    work = work.assign(
        _cr=pd.to_numeric(work[cols["crista"]], errors="coerce").fillna(-1.0)
    ).sort_values("_cr", ascending=False).drop(columns=["_cr"])
    return work.head(n)


def _esc(value):
    if _is_na(value):
        return ""
    return html.escape(str(value))


# All bundled downloads live under this subfolder in the report ZIP, so the top
# level contains ONLY report.html -- it is then self-evident that the reader opens
# report.html. The in-HTML links point here; the zip stores files under it.
DATA_SUBDIR = "data"


def _dl_href(name):
    """href to a bundled download (which lives under ``data/``). The visible link
    TEXT stays the bare filename; only the href is prefixed."""
    return f"{DATA_SUBDIR}/{_esc(name)}"


def _linkify_bundled_filenames(html, names):
    """Make EVERY inline ``<code>FILE</code>`` mention of a bundled file clickable
    (a download link into ``data/``), so any reference to a download in the report
    prose -- not just the dedicated download buttons -- is itself a link. ``names``
    is the set of files actually bundled, so a mention is only linked when the file
    exists. No inline ``<code>`` filename is pre-linked (see ``vc_html``), so a plain
    substring replace is safe and idempotent; a filename that is a substring of
    another never collides because the whole ``<code>...</code>`` token must match."""
    for name in names:
        if not name or name == "report.html":
            continue
        token = f"<code>{_esc(name)}</code>"
        if token in html:
            html = html.replace(
                token, f'<a href="{_dl_href(name)}" download>{token}</a>'
            )
    return html


# MAF footnote (report v2.4): explains every blank / em-dash MAF cell in the table
# and the curated download files. The panel name is READ FROM THE RUN (meta
# datasets), never hardcoded, so it stays correct for any variant database(s).
def maf_footnote(datasets=None):
    ds = datasets or "the variant panel"
    return (
        "MAF = combined-panel minor/alternate allele frequency. "
        "MAF blank (&mdash;) = reference off-target (no variant), an indel-derived "
        "variant (the frequency registry is SNP-only), or a SNP not in the frequency "
        "panel; for SNP variant off-targets the frequency is AC/AN over the genotyped "
        f"panel &mdash; here the merged {ds} union panel (the combined global "
        "AF), not a single ancestry. A variant present in the panel but whose "
        "source allele frequency is 0 (e.g. a secondary allele of a multiallelic "
        "site) is shown at a display floor of 1&times;10<sup>&minus;5</sup> so it "
        "renders on log-scale plots &mdash; read it as &ldquo;present, frequency "
        "effectively 0&rdquo;, not as a measured frequency of 10<sup>&minus;5</sup>."
    )


def build_table_html(top_df, cols, has_crista, datasets=""):
    """Scrollable inline top-N table, sorted by CFD desc; no JS (opens offline).

    Renders EXACTLY the ONE curated column set (``CURATED_COLUMNS``) that every
    download file uses -- so the table and top1000.tsv / panel_top100.tsv / the
    per-tier TSVs all show the same columns, including the annotation columns
    (Gene, Gene_distance_kb, GENCODE, ENCODE, DHS) and CRISTA when computed. Cells
    come from ``build_curated_frame`` (values resolved by name; missing -> "-").
    The aligned protospacer+PAM cell is wrapped in <code>; the MAF em-dash keeps
    its footnote. Extra columns don't break the scroll box / sticky header.
    """
    curated = build_curated_frame(top_df, cols, has_crista, start_rank=1)
    headers = list(curated.columns)
    head_html = "".join(f"<th>{_esc(h)}</th>" for h in headers)

    aligned_hdr = "Aligned_protospacer+PAM"
    maf_hdr = "MAF"
    maf_missing_seen = False
    body_rows = []
    for _idx, row in curated.iterrows():
        cells = []
        for h in headers:
            val = row[h]
            if h == aligned_hdr:
                cells.append(f"<code>{_esc(val)}</code>")
            elif h == maf_hdr:
                if val == CURATED_MISSING:
                    maf_missing_seen = True
                    cells.append("&mdash;")
                else:
                    cells.append(_esc(val))
            else:
                cells.append(_esc(val))
        body_rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")

    footnote = ""
    if maf_missing_seen:
        footnote = f'<p class="caption">{maf_footnote(datasets)}</p>'

    table = (
        '<div class="ottable-wrap"><table class="ottable">'
        f"<thead><tr>{head_html}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table></div>"
    )
    return table + footnote


def write_top1000_tsv(top_df, cols, has_crista, path):
    """Write the top-N rows as a CURATED-column TSV (shared curated schema)."""
    write_curated_tsv(top_df, cols, has_crista, path, start_rank=1)


# --------------------------------------------------------------------------- #
# HTML assembly
# --------------------------------------------------------------------------- #
_CSS = """
:root { color-scheme: light; }
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       margin: 0; padding: 26px 14px; color: #1a202c; line-height: 1.45;
       background-color: #eef2f7; }
.page { max-width: 1180px; margin: 0 auto; background: #ffffff; padding: 26px 34px 34px;
        border-radius: 10px; box-shadow: 0 2px 16px rgba(20,40,80,0.14); }
.report-header { display: flex; align-items: center; gap: 18px;
                 border-bottom: 2px solid #edf2f7; padding-bottom: 14px;
                 margin-bottom: 1.2em; }
.report-header img.logo { height: 92px; width: auto; }
.report-header .titles { flex: 1 1 auto; }
.report-header h1 { margin: 0 0 0.05em; }
.legend { display: flex; flex-direction: column; gap: 10px; margin: 0.5em 0 1em; }
.legend-item { display: flex; gap: 16px; align-items: baseline; border: 1px solid #e6edf5;
               border-radius: 6px; padding: 10px 14px; background: #fbfcfe; }
.legend-term { flex: 0 0 215px; font-weight: 600; color: #2c5282; }
.legend-def { flex: 1 1 auto; color: #2d3748; font-size: 0.92em; }
.legend-def ul { margin: 0.4em 0 0; padding-left: 1.2em; }
.legend-def li { margin: 0.15em 0; }
.legend code { font-family: SFMono-Regular, Menlo, Consolas, monospace; font-size: 0.9em;
               background: #edf2f7; padding: 0 3px; border-radius: 3px; }
h1 { font-size: 1.6em; margin-bottom: 0.1em; }
h2 { font-size: 1.25em; margin-top: 1.8em; border-bottom: 1px solid #e2e8f0;
     padding-bottom: 0.2em; }
.subtitle { color: #4a5568; margin-top: 0; }
.caption { color: #718096; font-size: 0.86em; margin: 0.3em 0 1.2em; }
.summary-grid { display: flex; flex-wrap: wrap; gap: 24px; align-items: flex-start; }
.summary-card { flex: 0 1 auto; max-width: 460px; min-width: 0; }
.summary-card .caption { overflow-wrap: break-word; }
.matrix-card { flex: 1 1 460px; min-width: 420px; }
.matrix-title { font-weight: 600; margin-bottom: 0.4em; }
table.summary-table { border-collapse: collapse; margin: 0.2em 0; width: 100%;
                      max-width: 520px; }
table.summary-table td { border: 1px solid #e2e8f0; padding: 6px 10px;
                         vertical-align: top; }
table.summary-table td.k { background: #f7fafc; font-weight: 600; width: 250px; }
/* The standalone 'Analysis inputs & criteria' card is full page-width, so let its
   table fill the card (matching the caption below it) instead of the 520px cap used
   for the narrow flex-column summary-card table. */
.inputs-card table.summary-table { max-width: none; }
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
.tier-downloads a.download { background: #4a5568; font-size: 0.86em;
                             padding: 6px 12px; }
.tier-downloads a.download:hover { background: #2d3748; }
.panel-grid { display: flex; flex-wrap: wrap; gap: 24px; }
table.thr-table { border-collapse: collapse; margin: 0.4em 0; }
table.thr-table th, table.thr-table td { border: 1px solid #e2e8f0; padding: 5px 12px; }
table.thr-table th { background: #f7fafc; text-align: left; }
table.thr-table td.num { text-align: right; font-variant-numeric: tabular-nums; }
a.tier-dl { color: #2b6cb0; text-decoration: none; font-size: 0.9em;
            font-variant-numeric: normal; }
a.tier-dl:hover { text-decoration: underline; }
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
    "This report is provided for research purposes only and on an \"AS IS\" basis, "
    "without warranty of any kind, express or implied, including without limitation "
    "any warranty of merchantability, fitness for a particular purpose, accuracy, "
    "completeness, or non-infringement. CRISPRme+ off-target predictions are "
    "computational and may contain false positives and false negatives; they are "
    "NOT a substitute for experimental validation and must not be the sole basis "
    "for any clinical, diagnostic, therapeutic, or regulatory decision. Results "
    "depend on the software version, algorithms, reference genome, PAM, search "
    "parameters, and variant datasets used, and MAY CHANGE as CRISPRme+, its methods, "
    "or the underlying data are updated or improved — a report reflects only the "
    "inputs and version stated above and is not a fixed or guaranteed output. To the "
    "maximum extent permitted by law, the authors, contributors, and their "
    "institutions make no warranties and accept no liability for any loss, injury, "
    "damage, cost, or claim, or for any decision, result, or clinical, diagnostic, "
    "therapeutic, regulatory, commercial, or other use of or reliance on this report "
    "or the software."
)


def _asset_data_uri(name):
    """Return a base64 data URI for a bundled asset (logo/background), or "".

    Resolves ``assets/<name>`` relative to the repo (this file lives in
    PostProcess/, assets/ is its sibling). Missing asset -> "" so the report
    degrades gracefully (no logo / plain background) instead of failing.
    """
    import base64
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "..", "assets", name),
                 os.path.join(here, "assets", name)):
        cand = os.path.abspath(cand)
        if os.path.isfile(cand):
            ext = name.rsplit(".", 1)[-1].lower()
            mime = {"svg": "image/svg+xml", "jpg": "image/jpeg",
                    "jpeg": "image/jpeg"}.get(ext, "image/" + ext)
            with open(cand, "rb") as fh:
                b64 = base64.b64encode(fh.read()).decode("ascii")
            return "data:%s;base64,%s" % (mime, b64)
    return ""


# Annotation legend (Section 7): plain-language meaning of every annotation
# column value, so a reviewer never has to guess what "dELS" or "Tier1_TSG" is.
_ANNOTATION_LEGEND = [
    ("gene_name", "Gene", "Symbol of the nearest gene to the off-target site (from "
     "the supplied gene annotation); use with <code>Gene_distance_kb</code> to see "
     "whether a site falls in or near a gene of interest."),
    ("gene_dist", "Gene_distance_kb", "Signed distance in kilobases to the nearest "
     "gene: <code>0.0</code> means the site lies WITHIN a gene (the GENCODE column "
     "is then a genic feature, not <code>intergenic</code>); a non-zero value&rsquo;s "
     "magnitude is the distance to the nearest gene boundary and its sign indicates "
     "the side (upstream vs downstream of the gene)."),
    ("gencode", "GENCODE", "Gene-model context of the site as labeled by the "
     "supplied GENCODE annotation: commonly <code>exon</code>, <code>CDS</code> "
     "(protein-coding sequence), <code>UTR</code>, <code>transcript</code>, "
     "<code>start_codon</code>/<code>stop_codon</code>, or <code>intergenic</code> "
     "(outside annotated genes); other feature terms (e.g. <code>protein_coding</code>, "
     "<code>intron</code>, <code>lincRNA</code>) may appear depending on the bundle."),
    ("dhs", "DHS", "DNase I Hypersensitive Site &mdash; a region of open, accessible "
     "chromatin (often regulatory), labeled by the tissue / organ system in "
     "which it is active (e.g. <code>Lymphoid</code>, <code>Neural</code>, "
     "<code>Tissue_invariant</code>)."),
    ("encode", "ENCODE / functional region", "Functional-region annotation "
     "(the built-in bundle uses ENCODE SCREEN v4 candidate <i>cis</i>-regulatory "
     "elements; a custom functional-region BED is reported in this column too). "
     "SCREEN v4 cCRE classes:<ul>"
     "<li><b>PLS</b> &mdash; promoter-like signature (at a transcription start "
     "site).</li>"
     "<li><b>pELS</b> &mdash; proximal enhancer-like signature (near a TSS, "
     "&le;2 kb).</li>"
     "<li><b>dELS</b> &mdash; distal enhancer-like signature (far from a TSS).</li>"
     "<li><b>CA-CTCF</b> &mdash; chromatin-accessible + CTCF binding (often "
     "insulator/architectural).</li>"
     "<li><b>CA-H3K4me3</b> &mdash; chromatin-accessible + H3K4me3, away from a "
     "TSS.</li>"
     "<li><b>CA-TF</b> &mdash; chromatin-accessible + transcription-factor "
     "binding.</li>"
     "<li><b>CA</b> &mdash; chromatin-accessible only.</li>"
     "<li><b>TF</b> &mdash; transcription-factor binding only (no high "
     "accessibility).</li></ul>"),
    ("cosmic", "COSMIC (Cancer Gene Census)", "The off-target lies in a gene catalogued in "
     "the COSMIC Cancer Gene Census &mdash; a curated list of genes causally "
     "implicated in cancer. Labels combine a confidence <b>tier</b> and the "
     "gene's documented <b>role(s)</b>:<ul>"
     "<li><b>Tier 1</b> &mdash; extensive, curated evidence of a direct, causal "
     "role in cancer.</li>"
     "<li><b>Tier 2</b> &mdash; strong but less-extensively-curated evidence.</li>"
     "<li><b>oncogene</b> &mdash; drives cancer when activated/over-active.</li>"
     "<li><b>TSG</b> &mdash; tumor-suppressor gene (drives cancer when "
     "lost/inactivated).</li>"
     "<li><b>fusion</b> &mdash; recurrently involved in cancer gene fusions.</li>"
     "</ul>A blank cell (&ndash;) means the site is not in a Cancer Gene Census "
     "gene. This flag is context for prioritization, not evidence of risk on its "
     "own."),
]


# Scores & columns legend (Section 7, always shown): what every score/column in the
# tables means, with a primary citation AND a plain-language rule of thumb -- so a
# bench scientist can judge HOW risky a site is and a reviewer can trace provenance.
_SCORE_LEGEND = [
    ("CFD",
     "<b>Cutting Frequency Determination</b> score (Doench <i>et al.</i>, "
     "<i>Nat. Biotechnol.</i> 2016) &mdash; a 0&ndash;1 estimate of how efficiently "
     "SpCas9 cuts a mismatched site relative to a perfect match (1 = as efficient as "
     "the on-target). Higher = more likely to be cut. <b>Rule of thumb:</b> treat "
     "CFD&nbsp;&ge;&nbsp;0.2 as worth validating (the recommended panel already "
     "hard-includes CFD&nbsp;&ge;&nbsp;0.5). <b>Caveat:</b> CFD was trained on "
     "single-base mismatches; CFD values for sites containing DNA/RNA bulges "
     "(insertions/deletions) are an extrapolation beyond the model&rsquo;s training "
     "domain &mdash; weigh CRISTA there."),
    ("CRISTA",
     "CRISTA score (Abadi <i>et al.</i>, <i>PLoS Comput. Biol.</i> 2017) &mdash; an "
     "<b>independent</b> machine-learning 0&ndash;1 estimate of cleavage propensity "
     "that also models bulges/indels. Higher = more likely to be cut. Reported "
     "alongside CFD because the two models can disagree; <b>a site scored high by "
     "EITHER model warrants validation</b>. CRISTA&rsquo;s scale is not directly "
     "comparable to CFD&rsquo;s &mdash; the same number means different things in "
     "each, so its threshold tiers are model-relative."),
    ("Mismatches",
     "Number of base substitutions between the guide and the genomic site (fewer = "
     "closer sequence match = generally higher off-target risk)."),
    ("Bulges",
     "DNA/RNA bulges (gaps) in the alignment &mdash; a bulge lets the guide pair with "
     "a genomic site of slightly different length (a small insertion/deletion)."),
    ("Mismatches+bulges",
     "Total edit distance (mismatches + bulges). The smallest values are the "
     "near-cognate sites that scoring models can under-weight, so the panel includes "
     "them regardless of score."),
    ("Perfect_match",
     "<b>Yes</b> when the site matches with <b>0 mismatches AND 0 bulges</b> &mdash; "
     "an exact genomic match, i.e. a candidate cut site indistinguishable from an "
     "intended on-target by sequence alone."),
    ("REF/ALT_origin",
     "Whether the site exists in the <b>reference</b> genome (<code>ref</code>) or is "
     "created/altered by a population <b>variant</b> (<code>alt</code>). Variant "
     "(alt) sites only exist in carriers of that variant."),
    ("PAM_creation",
     "Flagged when a variant creates a new PAM absent from the reference, enabling an "
     "off-target that reference-only tools miss."),
    ("MAF",
     "Minor-allele frequency of the contributing variant over the genotyped panel "
     "(blank for reference sites). A value of <b>1&times;10<sup>&minus;5</sup></b> is "
     "a display floor meaning &ldquo;present, frequency effectively&nbsp;0&rdquo;, "
     "not a measured frequency."),
    ("Aligned protospacer+PAM notation",
     "In the aligned-sequence column: <b>UPPERCASE</b> = a base matching the guide, "
     "<b>lowercase</b> = a mismatch, and a dash <code>-</code> = a bulge/gap."),
]


def build_score_legend_html():
    """Render the scores-&-columns legend (Section 7, always present)."""
    items = "".join(
        '<div class="legend-item"><div class="legend-term">%s</div>'
        '<div class="legend-def">%s</div></div>' % (term, definition)
        for term, definition in _SCORE_LEGEND
    )
    return (
        '<p class="caption">What the score and key columns in the tables above (and '
        "in the downloads) mean, so a site can be judged not just by <i>where</i> it "
        "is but <i>how</i> likely it is to be cut.</p>"
        '<div class="legend">' + items + "</div>"
    )


def build_annotation_legend_html(cols=None, df=None):
    """Render the annotation legend (Section 7). Only include an entry when its
    backing annotation column is actually present in this run (and, when ``df`` is
    given, carries at least one real value), so a run WITHOUT a given annotation is
    never falsely documented (e.g. a legacy/no-COSMIC or non-human/custom-genome
    run must not print a COSMIC Cancer Gene Census glossary). Returns "" when no
    entry qualifies, so the caller omits the whole section. ``cols=None`` keeps the
    legacy full legend (backward-compatible)."""
    def _active(key):
        if cols is None:
            return True
        if key not in cols:
            return False
        if df is not None and cols[key] in df.columns:
            return bool(df[cols[key]].map(lambda v: not _is_na(v)).any())
        return True

    entries = [(t, d) for (k, t, d) in _ANNOTATION_LEGEND if _active(k)]
    if not entries:
        return ""
    items = "".join(
        '<div class="legend-item"><div class="legend-term">%s</div>'
        '<div class="legend-def">%s</div></div>' % (term, definition)
        for term, definition in entries
    )
    return (
        '<p class="caption">What the annotation columns in the table above (and '
        "in the downloads) mean. Annotations describe the genomic context of "
        "each off-target site; they help prioritize which sites to examine "
        "first.</p>"
        '<div class="legend">' + items + "</div>"
    )


def render_next_steps_box(vp):
    """Plain-language 'What to do next' box -- turns the recommended panel into an
    actionable bench-validation workflow for a wet-lab reader."""
    panel_n = vp.get("panel_size", 0)
    n_perfect = int(vp.get("n_perfect", 0) or 0)
    max_pg = int(vp.get("max_perfect_per_guide", n_perfect) or 0)
    steps = []
    if max_pg >= 2:
        steps.append(
            "<b>Resolve the on-target first.</b> This guide has more than one perfect "
            "(0-mismatch / 0-bulge) match, so the intended on-target cannot be chosen "
            "from sequence alone &mdash; confirm it experimentally, or redesign the "
            "guide, before interpreting the rest."
        )
    steps.append(
        f"<b>Validate the recommended panel.</b> Design amplicon / rhAMP-Seq assays "
        f"for the ~{panel_n} sites in <code>panel_top100.tsv</code> "
        "(Section&nbsp;4) &mdash; the worst-case shortlist across all metrics."
    )
    steps.append(
        "<b>Prioritize within the panel</b> sites that are (i) high CFD or CRISTA, "
        "(ii) low edit-distance (few mismatches/bulges), and (iii) inside a gene "
        "&mdash; especially a COSMIC cancer gene."
    )
    steps.append(
        "<b>Read variant-created sites in context.</b> A site flagged "
        "<code>alt</code> only exists in carriers of that variant; its relevance "
        "depends on the allele frequency in your target population."
    )
    items = "".join(f'<li style="margin:0.35em 0">{s}</li>' for s in steps)
    return (
        '<div style="border:1px solid #9ae6b4;background:#f0fff4;border-radius:10px;'
        'padding:14px 18px;margin:1.2em 0">'
        '<div class="matrix-title">What to do next</div>'
        f'<ol style="margin:0.4em 0 0;padding-left:1.3em">{items}</ol></div>'
    )


def render_perfect_match_banner(vp):
    """Prominent banner for perfect (0 mismatch + 0 bulge) matches.

    A guide with several perfect genomic matches has no a-priori on-target -- each
    is an equally-efficient candidate cut site, and a perfect-match OFF-target is
    the highest-risk class. The red "no unambiguous on-target" warning fires only
    when a SINGLE guide has ``>= 2`` perfect matches (``max_perfect_per_guide``);
    a multi-guide run with exactly one perfect match per guide is NORMAL and gets
    the amber "presumed on-target(s)" note. ``0`` perfect matches -> empty string.
    """
    n = int(vp.get("n_perfect", 0) or 0)
    sites = vp.get("perfect_sites", []) or []
    # per-guide multiplicity drives ambiguity; default to n for legacy vp dicts
    max_pg = int(vp.get("max_perfect_per_guide", n) or 0)
    n_guides = int(vp.get("n_guides_with_perfect", 1) or 1)
    if n <= 0:
        return ""

    multi_guide = n_guides > 1

    def _origin_tag(s):
        # distinguish a UNIVERSAL reference match from a RARE variant-created one, so
        # a rare-allele perfect match is not presented as co-equal to the on-target.
        o = str(s.get("origin", "")).lower()
        m = s.get("maf")
        if o == "alt":
            if m is not None and m > 0:
                return f" <em>[variant-created, MAF&nbsp;{m:.2g} &mdash; only in carriers]</em>"
            return " <em>[variant-created &mdash; only in carriers]</em>"
        if o == "ref":
            return " <em>[reference &mdash; present in all genomes]</em>"
        return ""

    def _fmt(s):
        strand = f" ({_esc(str(s.get('strand', '')))})" if s.get("strand") else ""
        loc = f"{_esc(str(s.get('chrom', '?')))}:{_esc(str(s.get('pos', '?')))}{strand}"
        # label the guide when several guides contribute perfect matches
        if multi_guide and s.get("guide"):
            loc = f"{_esc(str(s.get('guide')))} &rarr; {loc}"
        return loc + _origin_tag(s)

    site_list = ", ".join(_fmt(s) for s in sites[:20])
    if len(sites) > 20:
        site_list += ", &hellip;"

    if max_pg >= 2:
        lead = (
            "One or more guides match" if multi_guide else "This guide matches"
        )
        return (
            '<div style="border:2px solid #dc2626;background:#fef2f2;'
            'padding:0.9em 1.1em;margin:1em 0;border-radius:6px">'
            '<p style="margin:0 0 0.4em;color:#991b1b;font-size:1.05em">'
            "<strong>&#9888; Multiple perfect matches &mdash; no unambiguous "
            "on-target.</strong></p>"
            f'<p style="margin:0">{lead} <strong>multiple</strong> genomic '
            "sites with <strong>0 mismatches and 0 bulges</strong>. The intended "
            "on-target cannot be determined from sequence alone &mdash; "
            "<strong>every one is a candidate cut site</strong> (a perfect-match "
            "off-target cuts as efficiently as the on-target). All are placed at "
            "the top of the validation panel (Section&nbsp;4) and flagged "
            "<code>Perfect_match&nbsp;=&nbsp;Yes</code> in the tables. "
            f"Sites: {site_list}.</p></div>"
        )
    # each guide has exactly one perfect match (one guide, or one-per-guide)
    if n == 1:
        lead_txt = (
            "one genomic site matches with 0 mismatches / 0 bulges "
            f"({site_list}) &mdash; the presumed on-target"
        )
    else:
        lead_txt = (
            f"{n} genomic sites match with 0 mismatches / 0 bulges, one per guide "
            f"({site_list}) &mdash; the presumed on-target of each guide"
        )
    return (
        '<div style="border-left:4px solid #d97706;background:#fffbeb;'
        'padding:0.6em 0.9em;margin:1em 0">'
        f"<p style=\"margin:0\"><strong>Perfect match:</strong> {lead_txt}, "
        "included in the validation panel and flagged "
        "<code>Perfect_match&nbsp;=&nbsp;Yes</code>.</p></div>"
    )


def render_html(
    job_id, summary_matrix_html, scatter_panels, pop_uri, validation_html,
    table_html, tsv_gz_name, top1000_name, footer_html,
    panel_top100_name=None, tier_downloads=None,
    hvdr_bundle_name=None, hvdr_n_regions=0, perfect_banner="",
    cooc_bundle_name=None, cooc_n_rows=0, cooc_n_cis=0,
    table_crista_html="", inputs_criteria_html="", legend_html=None,
    next_steps_html="",
):
    scatter_html = []
    for title, caption, uri in scatter_panels:
        scatter_html.append(
            f'<h3 style="margin:0.8em 0 0.2em">{_esc(title)}</h3>'
            f'<div class="plot"><img alt="{_esc(title)}" src="{uri}"></div>'
            f'<p class="caption">{caption}</p>'
        )
    scatter_block = "\n".join(scatter_html)

    panel_download = ""
    panel_caption = ""
    if panel_top100_name:
        panel_download = (
            f'\n  <a class="download" href="{_dl_href(panel_top100_name)}" download>'
            f"Recommended hybrid panel (TSV)</a>"
        )
        panel_caption = (
            f" The recommended hybrid worst-case top-100 panel is also bundled as "
            f"<code>{_esc(panel_top100_name)}</code>."
        )

    # highly-complex (high-variant-density) regions: download link + a prominent
    # callout so a report reader is fully aware these windows exist and where to look.
    hvdr_download = ""
    hvdr_callout = ""
    if hvdr_bundle_name:
        hvdr_download = (
            f'\n  <a class="download" href="{_dl_href(hvdr_bundle_name)}" download>'
            f"Highly complex regions near top sites ({hvdr_n_regions:,}) &mdash; BED</a>"
        )
        hvdr_callout = (
            f'<p class="caption" style="border-left:4px solid #d97706;'
            f'padding-left:0.7em;background:#fffbeb">'
            f"<strong>Highly complex (high-variant-density) regions:</strong> "
            f"{hvdr_n_regions:,} window(s) overlapping the top-ranked reported off-targets "
            f"carry so many overlapping variants that a single greedy worst-case alignment "
            f"is reported for each (flagged in the <code>High_complexity_region</code> "
            f"column of the table) &mdash; ADDITIONAL haplotype alignments may exist there. "
            f"Those regions (span, variant count, carriers, full IUPAC protospacer) are "
            f"bundled as <code>{_esc(hvdr_bundle_name)}</code>. The <strong>complete "
            f"genome-wide</strong> flag is in the <code>High_variant_density_region</code> "
            f"column of the integrated results (every site).</p>"
        )

    # SNP+indel cis co-occurrence companion: download link + a callout naming the
    # CONFIRMED-cis site count (a phased indel + nearby SNP on the same haplotype).
    cooc_download = ""
    cooc_callout = ""
    if cooc_bundle_name:
        cooc_download = (
            f'\n  <a class="download" href="{_dl_href(cooc_bundle_name)}" download>'
            f"SNP + indel cis co-occurrences ({cooc_n_cis:,} confirmed-cis) &mdash; TSV</a>"
        )
        cooc_callout = (
            f'<p class="caption" style="border-left:4px solid #2563eb;'
            f'padding-left:0.7em;background:#eff6ff">'
            f"<strong>SNP + indel cis co-occurrences:</strong> "
            f"{cooc_n_cis:,} confirmed-cis site(s) (of {cooc_n_rows:,} candidate "
            f"co-occurrences) where an indel and a nearby SNP fall on the <em>same "
            f"haplotype</em> (phased), so both edits are carried together by the same "
            f"individuals &mdash; a joint off-target that neither variant produces "
            f"alone. Each row lists the indel, the cis SNP (rsID), the phase, the "
            f"joint allele frequency, and the carrier sample(s); the full list is "
            f"bundled as <code>{_esc(cooc_bundle_name)}</code>.</p>"
        )

    # per-tier curated downloads (Section 5). ``tier_downloads`` is a list of
    # (label, filename) for every non-empty threshold tier that was bundled.
    tier_download_html = ""
    tier_caption = ""
    if tier_downloads:
        links = "".join(
            f'\n  <a class="download" href="{_dl_href(fname)}" download>'
            f"{_esc(label)}</a>"
            for label, fname in tier_downloads
        )
        tier_download_html = (
            f'\n<p class="tier-downloads">{links}\n</p>'
        )
        tier_caption = (
            " Per-threshold subsets (CFD and mismatch+bulge tiers, plus the "
            "variant-created subset) are bundled as curated-column TSVs so the "
            "panel can be expanded to any tier for review."
        )

    logo_uri = _asset_data_uri("crisprme-logo.svg") or _asset_data_uri("crisprme-logo.png")
    # High-resolution seamless tile (1254px, ~0.85 MB JPEG) shown at 640px via
    # background-size, so the pattern keeps its density but stays crisp on
    # high-DPI/retina screens (the old 640px PNG was upscaled ~2x and looked soft).
    # Fall back to the full web-UI PNG tile if the report JPEG isn't bundled.
    bg_uri = (_asset_data_uri("crisprme_bg_report.jpg")
              or _asset_data_uri("crisprme_bg_tile.png"))
    logo_html = (f'<img class="logo" src="{logo_uri}" alt="CRISPRme+ logo">'
                 if logo_uri else "")
    bg_style = (f"<style>body {{ background-image: url('{bg_uri}');"
                f" background-repeat: repeat; background-size: 640px 640px; }}</style>"
                if bg_uri else "")
    if legend_html is None:
        legend_html = build_score_legend_html()
    # Section 7 always present (the scores/columns legend is always relevant)
    legend_section = (
        f"\n<h2>7. Legend &mdash; scores, columns &amp; annotations</h2>\n{legend_html}"
        if legend_html else ""
    )
    crista_block = ""
    if table_crista_html:
        crista_block = (
            '<h3 style="margin:1.4em 0 0.3em 0">Ranked by CRISTA score</h3>\n'
            + table_crista_html
        )

    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(job_id)} CRISPRme+ off-target assessment report</title>
<style>{_CSS}</style>
{bg_style}
</head><body>
<div class="page">

<div class="report-header">
{logo_html}
<div class="titles">
<h1>Off-target assessment report</h1>
<p class="subtitle">CRISPRme+ &mdash; genome-wide off-target prediction accounting
for human genetic variation</p>
</div>
</div>

<h2>1. Summary</h2>
{summary_matrix_html}
{inputs_criteria_html}
{perfect_banner}

<h2>2. Key graphical report</h2>
<p class="caption">Reference vs variant off-target scores across the top-ranked
candidates. The same scatter is shown under multiple
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
{next_steps_html}

<h2>5. Downloads</h2>
<p>
  <a class="download" href="{_dl_href(tsv_gz_name)}" download>Complete raw integrated results (all columns, TSV gzip)</a>
  <a class="download" href="{_dl_href(top1000_name)}" download>Top-1000 off-targets (curated TSV)</a>{panel_download}{hvdr_download}{cooc_download}
</p>{tier_download_html}
{hvdr_callout}
{cooc_callout}
<p class="caption">All files are bundled alongside this HTML in the same ZIP;
the links resolve after unzipping on any machine. The top-1000 TSV, the panel,
and the per-tier subsets share the SAME curated, readable columns as the table
below; the complete integrated results (all columns) stays as the raw
<code>{_esc(tsv_gz_name)}</code>.{panel_caption}{tier_caption}</p>

<p class="caption" style="border-left:4px solid #2b6cb0;padding-left:0.7em;background:#eef6ff">
<strong>Need more detail?</strong> The interactive CRISPRme+ web interface goes beyond
this static report: it can generate a <strong>personalized off-target report for any
individual included in the variant panel</strong> (using that person's specific
genotypes), let you explore per-sample and per-population results interactively, and
run analyses against a <strong>custom personal genome assembly</strong>. Launch it
locally with <code>crisprme.py web-interface</code>.</p>

<h2>6. Top 1000 putative off-targets</h2>
<h3 style="margin:0.6em 0 0.3em 0">Ranked by CFD score</h3>
{table_html}
{crista_block}

{legend_section}

{footer_html}

</div>
</body></html>
"""


def build_footer(meta, version, tsv_basename):
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    crisprme_version = version or meta.get("version") or "n/a"
    return f"""<footer>
<p class="license"><strong>License.</strong> CRISPRme+ is <strong>free for academic
and non-profit research use</strong> (AGPL-3.0), for the user's own non-commercial
research and teaching. <strong>Any commercial or for-profit use &mdash; of CRISPRme+
or of any result or report it produces, including in a clinical trial, product, or
development program, and regardless of who ran the software &mdash; requires a
commercial license.</strong> To obtain a license, please contact
<a href="mailto:lpinello@mgh.harvard.edu">Luca Pinello (lpinello@mgh.harvard.edu)</a>.</p>
<p class="disclaimer">{_esc(DISCLAIMER)}</p>
<p class="feedback">Feedback or a bug to report? Please open an issue at
<a href="https://github.com/pinellolab/crisprme-plus/issues">github.com/pinellolab/crisprme-plus/issues</a>.</p>
<p>CRISPRme+ version: {_esc(crisprme_version)}
&nbsp;&middot;&nbsp; report generator v{_esc(REPORT_GENERATOR_VERSION)}
&nbsp;&middot;&nbsp; generated {_esc(stamp)}
&nbsp;&middot;&nbsp; source: {_esc(tsv_basename)}</p>
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
    """Best-effort CRISPRme version for the footer.

    Order: (1) installed package metadata (pip/conda), then (2) the canonical
    ``version = "X"`` assignment in the sibling ``crisprme.py`` -- the source /
    Docker install case, where CRISPRme is laid down by ``COPY .`` and has NO
    package metadata, so importlib.metadata fails and the footer used to read
    "n/a". Returns None only if neither is available.
    """
    try:
        from importlib.metadata import PackageNotFoundError, version as _v

        try:
            return _v("CRISPRme")
        except PackageNotFoundError:
            pass
    except Exception:  # noqa: BLE001
        pass
    # source / Docker install: parse the canonical version out of crisprme.py.
    # PostProcess/ sits under the source root (../crisprme.py) and, in the
    # Bioconda/Docker bin layout, crisprme.py is also copied alongside -- try both.
    try:
        import re

        here = os.path.dirname(os.path.abspath(__file__))
        for cand in (
            os.path.join(here, os.pardir, "crisprme.py"),
            os.path.join(here, "crisprme.py"),
        ):
            if os.path.isfile(cand):
                with open(cand, "r", errors="replace") as fh:
                    for line in fh:
                        m = re.match(r'version\s*=\s*["\']([^"\']+)["\']', line)
                        if m:
                            return m.group(1)
    except Exception:  # noqa: BLE001
        pass
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
    drop_maf=False,
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
    global _DROP_MAF
    _DROP_MAF = bool(drop_maf)
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
    # drop annotation curated columns whose source is absent from THIS run (no
    # all-"-" COSMIC/ENCODE/... column implying a screen that was never performed)
    global _PRESENT_ANN_KINDS
    _PRESENT_ANN_KINDS = {k for k in _ANNOTATION_KINDS if k in cols}

    # De-duplicate REFERENCE off-target rows (locus-completeness can emit the same
    # variant-independent reference site once per co-located haplotype). Applied to
    # the report's working frame only; the raw integrated_results.tsv.gz bundled in
    # the ZIP is copied from the source file and stays complete.
    df = dedupe_reference_rows(df, cols)

    meta = build_summary_meta(
        result_dir, integrated_tsv, df, cols, params_override=params_override
    )
    version = meta.get("version") or _package_version()
    has_crista = crista_computed(df, cols)

    spec_score = read_specificity_score(result_dir, job_id, meta["guides"])

    fn = _parse_results_filename(integrated_tsv)
    sample_superpop = load_sample_superpop(samplesid_dir, fn.get("datasets", ""))
    sample_dataset = load_sample_dataset(samplesid_dir)  # native per-db provenance (dynamic)

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
        pop_uri = plot_population(df, cols, sample_superpop, sample_dataset)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: population plot unavailable: {exc}\n")
        pop_uri = _placeholder_uri("Population plot unavailable")

    # ---- SECTION 4: validation panel (hybrid worst-case top-100) ------------
    # Compute the panel + the per-tier subsets FIRST, then write each non-empty
    # tier as a CURATED-column TSV (gzip when > ~2 MB) so the section-4 table can
    # link the exact bundled filenames. We stage the files into ``staging`` and
    # record: panel_top100_name, tier_links (key->filename for the sec-4 links),
    # tier_downloads (label,filename for section 5), and staged_tier_paths.
    staging = tempfile.mkdtemp(prefix="crisprme_report_")
    tsv_gz_name = "integrated_results.tsv.gz"
    top1000_name = "top1000.tsv"
    panel_top100_name = None
    tier_links = {}
    tier_downloads = []
    staged_tier_paths = []

    def _stage_curated(sub_df, base_name):
        """Write ``sub_df`` as a curated TSV under ``staging``; gzip if > ~2 MB.

        Returns the bundled filename (basename), or None on failure. Writes a
        plain .tsv first, then re-packs to .tsv.gz when it exceeds TIER_GZIP_BYTES
        (plain .tsv otherwise, per spec)."""
        plain = os.path.join(staging, base_name)
        try:
            write_curated_tsv(sub_df, cols, has_crista, plain, start_rank=1)
        except Exception as exc:  # noqa: BLE001 - never abort on a bundled TSV
            sys.stderr.write(f"generate-report: {base_name} unavailable: {exc}\n")
            return None
        if os.path.getsize(plain) > TIER_GZIP_BYTES:
            gz_name = base_name + ".gz"
            gz = os.path.join(staging, gz_name)
            with open(plain, "rb") as src, gzip.open(gz, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            os.remove(plain)
            staged_tier_paths.append(gz)
            return gz_name
        staged_tier_paths.append(plain)
        return base_name

    perfect_banner = ""
    try:
        vp = build_validation_panel(df, cols)
        perfect_banner = render_perfect_match_banner(vp)
        panel_top100_df = vp.get("panel_df")

        # panel_top100.tsv (curated) -- always a hard bundle when non-empty.
        # Linked in section 4 (tier_links) and section 5 (panel_download in
        # render_html), so it is NOT added to tier_downloads (avoid duplicate).
        if panel_top100_df is not None and len(panel_top100_df) > 0:
            fname = _stage_curated(panel_top100_df, PANEL_TOP100_NAME)
            if fname:
                panel_top100_name = fname
                tier_links["panel"] = fname

        # per-threshold tiers -- only the non-empty ones
        for tier in vp.get("tiers", []):
            sub = tier["df"]
            if sub is None or len(sub) == 0:
                continue
            fname = _stage_curated(sub, tier["filename"])
            if not fname:
                continue
            tier_links[tier["key"]] = fname
            label = f"{tier['label'].replace('&ge;', '>=').replace('&le;', '<=')} ({len(sub):,})"
            tier_downloads.append((label, fname))

        validation_html = render_validation_panel(
            vp, panel_tsv_name=panel_top100_name, tier_links=tier_links
        )
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: validation panel unavailable: {exc}\n")
        validation_html = "<p>Validation panel unavailable.</p>"

    # ---- SECTION 6: top-1000 table (curated columns incl. annotations) ------
    try:
        table_html = build_table_html(top_df, cols, has_crista, datasets=meta.get("datasets", ""))
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"generate-report: table unavailable: {exc}\n")
        table_html = "<p>Top-1000 table unavailable.</p>"
    # Second top-1000 table ranked by CRISTA (only when CRISTA was computed), with
    # a bundled curated TSV alongside top1000.tsv.
    table_crista_html = ""
    top_crista_df = None
    if has_crista:
        try:
            top_crista_df = select_top_crista(df, cols, n=top_n)
            if len(top_crista_df):
                table_crista_html = build_table_html(
                    top_crista_df, cols, has_crista, datasets=meta.get("datasets", "")
                )
                tier_downloads.append(
                    ("Top-1000 by CRISTA (curated TSV)", "top1000_crista.tsv")
                )
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"generate-report: CRISTA table unavailable: {exc}\n")
            table_crista_html, top_crista_df = "", None

    # ---- high-variant-density ("highly complex") regions BED (#144) ----------
    # Merge the per-chromosome beds the search wrote into ONE bundled file (single
    # header) so report readers can see + dig into every dense window where a greedy
    # worst-case alignment is reported and additional alignments may exist. Computed
    # BEFORE render_html so the downloads section can link it + show the count.
    hvdr_bundle_name = None
    hvdr_n_regions = 0
    try:
        import bisect

        _hvdr_src_dir = result_dir if result_dir else os.path.dirname(integrated_tsv)
        _hvdr_files = sorted(
            glob.glob(os.path.join(_hvdr_src_dir, "*high_variant_density_regions.bed"))
        )
        # Scope the report's HVDR list to the TOP-N reported off-targets. Genome-wide,
        # a permissive search flags tens of thousands of dense windows, most of them
        # low-relevance -- so keep only regions overlapping a top-N site here; the FULL
        # per-site flag stays in the High_variant_density_region column of the
        # integrated results. Fall back to the full list if positions can't be read.
        _top_pos = {}
        try:
            if "chrom" in cols and "pos" in cols and len(top_df):
                _cs = top_df[cols["chrom"]].astype(str).tolist()
                _ps = pd.to_numeric(top_df[cols["pos"]], errors="coerce").tolist()
                for _c, _p in zip(_cs, _ps):
                    if pd.notna(_p):
                        _top_pos.setdefault(_c, []).append(int(_p))
                for _c in _top_pos:
                    _top_pos[_c].sort()
        except Exception:  # noqa: BLE001
            _top_pos = {}
        _filter_hvdr = bool(_top_pos)
        _TOL = 25  # guide+PAM slack for matching a site to its dense window

        def _hits_top(chrom, start, end):
            arr = _top_pos.get(chrom)
            if not arr:
                return False
            i = bisect.bisect_left(arr, start - _TOL)
            return i < len(arr) and arr[i] <= end + _TOL

        if _hvdr_files:
            _hvdr_path = os.path.join(staging, "high_variant_density_regions.bed")
            _hvdr_seen = set()  # dedup: per-chrom beds carry duplicate region rows
            with open(_hvdr_path, "w") as _out:
                _out.write(
                    "#chrom\tstart\tend\tguide\tn_variants\t"
                    "samples_with_alt\tiupac_protospacer\n"
                )
                for _bf in _hvdr_files:
                    for _ln in open(_bf):
                        if _ln.startswith("#") or not _ln.strip():
                            continue
                        if _filter_hvdr:
                            _p = _ln.split("\t")
                            try:
                                if not _hits_top(_p[0], int(_p[1]), int(_p[2])):
                                    continue
                            except (IndexError, ValueError):
                                continue
                        # count + write each unique region ONCE (the per-chrom writer
                        # emits a row per flagged window per alignment/pass, so the
                        # same window recurs; key on the full normalized row).
                        _key = _ln.rstrip("\n")
                        if _key in _hvdr_seen:
                            continue
                        _hvdr_seen.add(_key)
                        _out.write(_ln if _ln.endswith("\n") else _ln + "\n")
                        hvdr_n_regions += 1
            if hvdr_n_regions:
                hvdr_bundle_name = "high_variant_density_regions.bed"
            else:
                os.remove(_hvdr_path)
    except Exception as exc:  # noqa: BLE001 - never abort on the sidecar bundle
        sys.stderr.write(f"generate-report: HVDR bed bundle unavailable: {exc}\n")
        hvdr_bundle_name = None

    # ---- SNP+indel cis co-occurrence companion (*.indel_snp_cooc.tsv) --------
    # The indel post-analysis writes ONE 12-column *.indel_snp_cooc.tsv per
    # chromosome (header at analisi_indels_NNN.py); nothing merges them, so the
    # flagship SNP+indel cis co-occurrence output is otherwise invisible in the
    # report. Concat them here into a single bundled file (keep the FIRST header
    # only) and count the CONFIRMED-cis rows (phase field contains "cis") so the
    # Downloads section can link the file + show the count. Mirrors the HVDR merge
    # above; never aborts the report on a malformed sidecar.
    cooc_bundle_name = None
    cooc_n_rows = 0
    cooc_n_cis = 0
    try:
        _cooc_src_dir = result_dir if result_dir else os.path.dirname(integrated_tsv)
        _cooc_files = sorted(
            glob.glob(os.path.join(_cooc_src_dir, "*indel_snp_cooc.tsv"))
        )
        if _cooc_files:
            _cooc_path = os.path.join(staging, "indel_snp_cooc.tsv")
            _cooc_header_written = False
            _cooc_seen = set()  # dedup: the writer emits a cooc row per alignment
            #                     pass for the same off-target (~10x duplication on
            #                     real data), so count + bundle each UNIQUE row once.
            with open(_cooc_path, "w") as _out:
                for _cf in _cooc_files:
                    with open(_cf) as _src:
                        for _i, _ln in enumerate(_src):
                            if _i == 0:  # each per-chrom file has its own header
                                if not _cooc_header_written:
                                    _out.write(_ln if _ln.endswith("\n") else _ln + "\n")
                                    _cooc_header_written = True
                                continue
                            if not _ln.strip():
                                continue
                            _key = _ln.rstrip("\n")
                            if _key in _cooc_seen:
                                continue
                            _cooc_seen.add(_key)
                            _out.write(_ln if _ln.endswith("\n") else _ln + "\n")
                            cooc_n_rows += 1
                            # phase (col index 8, 0-based) is CONFIRMED or PUTATIVE
                            # -- EVERY row is already a cis co-occurrence; CONFIRMED
                            # means the same-haplotype phasing is proven (all carriers
                            # phased). Count the CONFIRMED subset (indel_snp_cis.py).
                            _parts = _key.split("\t")
                            if len(_parts) > 8 and _parts[8].strip().upper() == "CONFIRMED":
                                cooc_n_cis += 1
            if cooc_n_rows:
                cooc_bundle_name = "indel_snp_cooc.tsv"
            else:
                os.remove(_cooc_path)
    except Exception as exc:  # noqa: BLE001 - never abort on the sidecar bundle
        sys.stderr.write(f"generate-report: indel_snp_cooc bundle unavailable: {exc}\n")
        cooc_bundle_name = None

    # ---- FOOTER (unnumbered; section 7 is the annotation legend, built in
    #      render_html via build_annotation_legend_html) -----------------------
    footer_html = build_footer(meta, version, os.path.basename(integrated_tsv))

    # inputs/criteria box built HERE (after tiers staged) so the variant_created
    # link uses the ACTUAL bundled filename (.tsv or .tsv.gz), never a broken guess.
    # per-dataset panel sizes for the 'Variants included' row (reuse the already
    # loaded native-provenance map; cheap -- just tallies sample IDs)
    _ds_counts = {}
    for _lbl in sample_dataset.values():
        _ds_counts[_lbl] = _ds_counts.get(_lbl, 0) + 1
    _reg_vc = _registry_variant_count(result_dir, meta)
    inputs_criteria_html = render_inputs_criteria(
        meta, tier_links.get("variant_created"),
        dataset_counts=_ds_counts,
        variant_count=_reg_vc,
    )
    html_doc = render_html(
        job_id, summary_matrix_html, scatter_panels, pop_uri, validation_html,
        table_html, tsv_gz_name, top1000_name, footer_html,
        panel_top100_name=panel_top100_name,
        tier_downloads=tier_downloads,
        hvdr_bundle_name=hvdr_bundle_name,
        hvdr_n_regions=hvdr_n_regions,
        cooc_bundle_name=cooc_bundle_name,
        cooc_n_rows=cooc_n_rows,
        cooc_n_cis=cooc_n_cis,
        perfect_banner=perfect_banner,
        next_steps_html=render_next_steps_box(vp),
        table_crista_html=table_crista_html,
        inputs_criteria_html=inputs_criteria_html,
        legend_html=(
            '<h3 style="margin:0.6em 0 0.3em">Scores &amp; columns</h3>'
            + build_score_legend_html()
            + (
                '<h3 style="margin:1.2em 0 0.3em">Annotations</h3>' + _ann_legend
                if (_ann_legend := build_annotation_legend_html(cols, df)) else ""
            )
        ),
    )

    # make every inline <code>FILE</code> reference to a bundled download clickable
    _linkable = {tsv_gz_name, top1000_name}
    if panel_top100_name:
        _linkable.add(panel_top100_name)
    if hvdr_bundle_name:
        _linkable.add(hvdr_bundle_name)
    if cooc_bundle_name:
        _linkable.add(cooc_bundle_name)
    if top_crista_df is not None and len(top_crista_df):
        _linkable.add("top1000_crista.tsv")
    _linkable.update(v for v in tier_links.values() if v)
    html_doc = _linkify_bundled_filenames(html_doc, _linkable)

    # ---- stage the remaining files and zip -j (flat) -----------------------
    try:
        html_path = os.path.join(staging, "report.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(html_doc)

        # top1000.tsv -- curated columns (same schema as the table + downloads)
        top1000_path = os.path.join(staging, top1000_name)
        try:
            write_top1000_tsv(top_df, cols, has_crista, top1000_path)
        except Exception as exc:  # noqa: BLE001 - never abort on the bundled TSV
            sys.stderr.write(f"generate-report: top1000.tsv unavailable: {exc}\n")
            with open(top1000_path, "w") as handle:
                handle.write("\t".join(curated_headers(has_crista)) + "\n")

        # top1000_crista.tsv -- the CRISTA-ranked companion (same curated schema)
        if top_crista_df is not None and len(top_crista_df):
            try:
                write_top1000_tsv(
                    top_crista_df, cols, has_crista,
                    os.path.join(staging, "top1000_crista.tsv"),
                )
            except Exception as exc:  # noqa: BLE001 - never abort on the bundled TSV
                sys.stderr.write(f"generate-report: top1000_crista.tsv unavailable: {exc}\n")

        # the complete RAW results (all 85 columns) stay as the gzip
        gz_path = os.path.join(staging, tsv_gz_name)
        if integrated_tsv.endswith(".gz"):
            shutil.copyfile(integrated_tsv, gz_path)
        else:
            with open(integrated_tsv, "rb") as src, gzip.open(gz_path, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)

        bundle = [html_path, gz_path, top1000_path]
        _crista_path = os.path.join(staging, "top1000_crista.tsv")
        if os.path.isfile(_crista_path):
            bundle.append(_crista_path)
        bundle += staged_tier_paths
        if hvdr_bundle_name:
            bundle.append(os.path.join(staging, hvdr_bundle_name))
        if cooc_bundle_name:
            bundle.append(os.path.join(staging, cooc_bundle_name))

        # machine-readable run manifest (IND traceability) + the raw .Params.txt, so
        # the ZIP is self-sufficient for re-execution. Never break the report on it.
        try:
            _manifest = {
                "crisprme_version": meta.get("version") or "2.4.0",
                "report_generator": "v2.4",
                "generated_date": meta.get("date"),
                "guides": meta.get("guides"),
                "nuclease": meta.get("nuclease"),
                "pam": meta.get("pam"),
                "genome": meta.get("genome"),
                "variant_datasets": meta.get("datasets"),
                "search": {
                    "mismatches": meta.get("mm"),
                    "bulges_dna": meta.get("bdna"),
                    "bulges_rna": meta.get("brna"),
                    "max_total_edits": meta.get("max_edits"),
                    "observed_max_mismatches_plus_bulges": meta.get("obs_max_mmb"),
                },
                "counts": {
                    "total_sites": meta.get("n_total"),
                    "off_targets": meta.get("n_offtarget"),
                    "on_target_perfect_matches": (
                        vp.get("n_perfect") if isinstance(vp, dict) else None),
                    "variant_created_off_targets": (
                        vp.get("n_variant") if isinstance(vp, dict) else None),
                    "recommended_panel_size": (
                        vp.get("panel_size") if isinstance(vp, dict) else None),
                },
                "database": _reg_vc,  # {n_records SNPs, n_indels, databases} or None
                "source_integrated_results": os.path.basename(integrated_tsv),
                "crispritz_note": (
                    "search-engine (CRISPRitz) version is recorded in the build; the "
                    "pinellolab/crisprme:v2.4.0 image builds CRISPRitz v2.8.2. See "
                    "Params.txt for the exact run parameters."
                ),
            }
            _man_path = os.path.join(staging, "run_manifest.json")
            with open(_man_path, "w") as _mf:
                json.dump(_manifest, _mf, indent=2, sort_keys=True)
            bundle.append(_man_path)
            if result_dir:
                for _pn in (".Params.txt", "Params.txt"):
                    _pp = os.path.join(result_dir, _pn)
                    if os.path.isfile(_pp):
                        _dst = os.path.join(staging, "Params.txt")
                        shutil.copyfile(_pp, _dst)
                        bundle.append(_dst)
                        break
        except Exception as exc:  # noqa: BLE001 - manifest is optional
            sys.stderr.write(f"generate-report: run_manifest unavailable: {exc}\n")

        if os.path.exists(out_zip):
            os.remove(out_zip)
        os.makedirs(os.path.dirname(out_zip), exist_ok=True)
        # Layout: report.html at the TOP level, every other bundled file under
        # data/ -- so unzipping shows a single obvious report.html plus a data/
        # folder. The in-HTML links point at data/<name> (see _dl_href), so the
        # report stays fully self-contained and openable in place.
        import zipfile

        def _arcname(path):
            base = os.path.basename(path)
            return base if base == "report.html" else f"{DATA_SUBDIR}/{base}"

        with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            for p in bundle:
                zf.write(p, _arcname(p))
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
    parser.add_argument(
        "--no-maf",
        dest="no_maf",
        action="store_true",
        help="Omit the MAF column entirely (table + every download) -- use when "
        "allele frequencies are not yet finalized so a MAF column would mislead.",
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
        drop_maf=args.no_maf,
    )
    sys.stdout.write(f"Report written: {out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
