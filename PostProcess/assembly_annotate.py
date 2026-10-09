"""hg38 functional annotation for assembly-search's reconciled results.

Assembly-search runs each haplotype through `complete-search` with no
`--annotation`, so nothing in the reconciled `combined_hg38.tsv` says what
genomic feature (gene, regulatory element, ...) a site falls in. Every
mappable row (`both`, `paternal_only`, `maternal_only`) has an hg38
coordinate from liftOver, so it can be annotated against the same hg38
annotation bundle `complete-search` uses -- once, here, after reconciliation,
not once per haplotype in each haplotype's own coordinate space (an hg38
annotation file has no meaning there).

Sites with no hg38 coordinate (one-sided non-mappable, and the
`both_haplotype_private` sites) get no annotation: they aren't in hg38 space
at all, so there is nothing to look up.

The annotation is the feature list at the *lifted* hg38 locus. In a region
where the haplotype diverges from hg38 it describes the corresponding
reference locus, not the haplotype's own sequence.

Reuses `annotation.load_annotation_bed()` / `annotation.annotate_target()`
(the exact lookup `complete-search` uses) rather than reimplementing the
tabix query.

Closest-gene annotation (`--gene_annotation`) is a separate, OPTIONAL step
handled by `add_closest_gene_columns()`. It is a nearest-feature query rather
than an overlap query, so it uses BEDOPS `closest-features` exactly as
`complete-search` does in `post_process.sh`, against the same per-transcript
GENCODE BED whose attributes column carries `gene_id=`/`gene_name=`.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Optional, Tuple

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from annotation import (  # noqa: E402
    ANNOTATION_SPLIT_COLS,
    annotate_target,
    load_annotation_bed,
    split_annotation_labels,
)
from utils import gene_region_class  # noqa: E402

ANNOTATION_COL = "Annotation"
NO_ANNOTATION = "n"  # same "nothing here" marker complete-search's annotation uses
GENE_REGION_COL = "Annotation_gene_region"
# Per-kind columns, in the order complete-search's own results table lists them.
# `Annotation` (the whole comma-joined feature list) is KEPT alongside these: it is
# what the web results table shows, and it is the only place a label lives if a
# future bundle introduces a kind this split doesn't know about.
_SPLIT_COL_ORDER = (
    "Annotation_GENCODE",
    GENE_REGION_COL,
    "Annotation_ENCODE",
    "Annotation_DHS",
    "Annotation_COSMIC",
    "Annotation_INTOGEN",
    "Annotation_personal",
)

# Closest-gene columns, named exactly as complete-search's own
# `resultIntegrator.py` names them, so `generate_report._COLS` resolves them
# with no report-side change (`gene_name` -> Gene, `gene_dist` -> Gene_distance_kb).
CLOSEST_GENE_NAME_COL = "Annotation_closest_gene_name"
CLOSEST_GENE_ID_COL = "Annotation_closest_gene_ID"
CLOSEST_GENE_DIST_COL = "Annotation_closest_gene_distance_(kb)"
CLOSEST_GENE_COLS = (CLOSEST_GENE_NAME_COL, CLOSEST_GENE_ID_COL, CLOSEST_GENE_DIST_COL)

# complete-search's fallback: when a site overlaps NO gencode feature but a
# closest gene exists at a non-zero distance, its GENCODE cell reads
# "intergenic" (resultIntegrator.py sets it in the closest-gene block, then
# overwrites it from the overlap labels only `if len(gencode_annotations)` --
# so it is a fallback, never a clobber). Mirrored here with the same precedence.
GENCODE_SPLIT_COL = "Annotation_GENCODE"
INTERGENIC = "intergenic"

# Genomic-side aligned target (gaps = RNA bulges shorter than the guide,
# DNA bulges add bases) -- its ungapped length is the span the site really
# covers on the genome, the same span `annotation.compute_target_coords()`
# uses for complete-search.
_TARGET_COLS = (
    "Aligned_protospacer+PAM_REF_(fewest_mm+b)_paternal",
    "Aligned_protospacer+PAM_REF_(fewest_mm+b)_maternal",
)
# Column the new one is placed after, so it sits with the hg38 coordinates
# in the results table; appended at the end if that column isn't there.
_INSERT_AFTER = "hg38_end_maternal"


def _target_length(row: pd.Series) -> int:
    for col in _TARGET_COLS:
        seq = row.get(col)
        if isinstance(seq, str) and seq:
            return len(seq.replace("-", ""))
    return 1  # no aligned target on either side: annotate the lifted point


def annotate_combined(combined: pd.DataFrame, annotation_path: str) -> pd.Series:
    """Feature list (comma-separated, sorted) at each row's hg38 locus.

    Rows without an hg38 coordinate get NaN (not "n"): "n" means "looked, and
    nothing overlaps", which is a different statement from "not in hg38".
    """
    annotation = load_annotation_bed(annotation_path)
    try:
        out = pd.Series(index=combined.index, dtype=object)
        has_hg38 = combined["hg38_chr"].notna() & combined["hg38_start"].notna()
        for idx in combined.index[has_hg38]:
            row = combined.loc[idx]
            start = int(row["hg38_start"])
            stop = start + _target_length(row)
            found = annotate_target(str(row["hg38_chr"]), start, stop, annotation)
            out.at[idx] = found if found else NO_ANNOTATION
        return out
    finally:
        annotation.close()


def split_annotation_frame(values: pd.Series) -> pd.DataFrame:
    """Per-kind annotation columns derived from the comma-joined feature lists.

    complete-search's results table carries one column PER annotation kind
    (`Annotation_GENCODE`, `Annotation_DHS`, ...), and the report's curated view
    reads those columns by name. Assembly-search used to emit only the combined
    `Annotation` string, so every one of those curated cells rendered "-" even
    though the annotation had been performed -- the report then showed six screens
    as "not done" when they had been. Splitting here, with the same classification
    complete-search uses, is what makes those cells populate; the report itself
    needs no assembly-specific knowledge.

    A kind is only given a column when at least one row actually has a value for
    it. An all-blank column would reintroduce exactly the "declared but never
    populated" problem this split exists to remove -- and the report's own
    `_PRESENT_ANN_KINDS` filter keys off column presence, so omitting the column is
    also what makes it correctly drop the matching curated column.

    Args:
        values: Comma-joined feature lists, as `annotate_combined` returns them
            (NaN where the row has no hg38 coordinate).

    Returns:
        A frame indexed like `values`, holding only the non-empty kind columns.
    """
    split = {col: [] for col in ANNOTATION_SPLIT_COLS}
    regions = []
    for value in values:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            # no hg38 coordinate -> nothing was looked up (distinct from "looked,
            # found nothing"), so every kind stays blank for this row
            for col in ANNOTATION_SPLIT_COLS:
                split[col].append("")
            regions.append("")
            continue
        buckets = split_annotation_labels(value)
        for col in ANNOTATION_SPLIT_COLS:
            split[col].append(buckets[col])
        # same helper, and the same `missing`, resultIntegrator uses -- no
        # closest-gene distance exists for an assembly run, so it is not passed
        regions.append(gene_region_class(buckets["Annotation_GENCODE"], missing="NA"))
    split[GENE_REGION_COL] = regions
    frame = pd.DataFrame(split, index=values.index)
    populated = [c for c in _SPLIT_COL_ORDER if c in frame.columns and frame[c].any()]
    return frame[populated]


def add_annotation_column(
    combined: pd.DataFrame, annotation_path: Optional[str]
) -> pd.DataFrame:
    """Returns `combined` with the `Annotation` column and the per-kind
    `Annotation_*` columns added (or unchanged if no annotation file was given)."""
    if not annotation_path:
        return combined
    combined = combined.copy()
    values = annotate_combined(combined, annotation_path)
    split = split_annotation_frame(values)
    stale = [ANNOTATION_COL, *ANNOTATION_SPLIT_COLS, GENE_REGION_COL]
    combined = combined.drop(columns=[c for c in stale if c in combined.columns])
    if _INSERT_AFTER in combined.columns:
        pos = list(combined.columns).index(_INSERT_AFTER) + 1
    else:
        pos = len(combined.columns)
    combined.insert(pos, ANNOTATION_COL, values)
    for offset, col in enumerate(split.columns, start=1):
        combined.insert(pos + offset, col, split[col])
    return combined


def _require_bedops() -> None:
    """Fails with a clear message if BEDOPS isn't installed.

    `closest-features`/`sort-bed` are the same binaries complete-search's
    `post_process.sh` already depends on, so a working install has them; this
    only turns a cryptic FileNotFoundError into an actionable one.
    """
    missing = [t for t in ("sort-bed", "closest-features") if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"{' and '.join(missing)} not found on PATH -- required for "
            "--gene_annotation (BEDOPS). Install with `mamba install -c bioconda bedops`."
        )


def _parse_gene_attributes(attrs: str) -> Tuple[str, str]:
    """`(gene_name, gene_id)` from a GFF attributes field, or `("", "")`.

    Same parse complete-search uses (`resultIntegrator.py`): split the
    attribute string on ';' and read the `gene_name=`/`gene_id=` entries. The
    GENCODE BED's records are per-transcript features (exon, five_prime_UTR,
    ...), so these attributes are what turns "nearest feature" into "nearest
    gene".
    """
    name = gid = ""
    for field in attrs.split(";"):
        field = field.strip()
        if field.startswith("gene_name="):
            name = field.split("=", 1)[1]
        elif field.startswith("gene_id="):
            gid = field.split("=", 1)[1]
    return name, gid


def _sorted_plain_bed(path: str, tmpdir: str) -> str:
    """Decompresses (if needed) and `sort-bed`s an annotation file.

    Sorted unconditionally, as `post_process.sh` does: `closest-features`
    requires BEDOPS sort order and silently misbehaves without it, and the
    caller may hand us either a plain .bed or a .bed.gz.
    """
    plain = path
    if path.endswith(".gz"):
        plain = os.path.join(tmpdir, "gene_annotation.bed")
        with open(plain, "wb") as dst:
            if subprocess.call(["gunzip", "-c", path], stdout=dst) != 0:
                raise RuntimeError(f"failed decompressing {path}")
    out = os.path.join(tmpdir, "gene_annotation.sorted.bed")
    with open(out, "wb") as dst:
        if subprocess.call(["sort-bed", plain], stdout=dst) != 0:
            raise RuntimeError(f"sort-bed failed on {plain}")
    return out


def closest_gene_frame(combined: pd.DataFrame, gene_annotation_path: str) -> pd.DataFrame:
    """Closest-gene name/ID/distance for every row with an hg38 coordinate.

    Rows with no hg38 coordinate (one-sided non-mappable and
    `both_haplotype_private`) are not queried and come back NaN: they are not
    in hg38 space, so "nearest hg38 gene" is not a statement we can make.

    The queried span is `hg38_start` to `hg38_start + ungapped target length`,
    the same span `annotate_combined()` uses and the same one complete-search
    builds (`$7 + length($3)` in `post_process.sh`).
    """
    _require_bedops()
    out = pd.DataFrame(
        {c: pd.Series(index=combined.index, dtype=object) for c in CLOSEST_GENE_COLS}
    )
    has_hg38 = combined["hg38_chr"].notna() & combined["hg38_start"].notna()
    if not has_hg38.any():
        return out
    tmpdir = tempfile.mkdtemp(prefix="crisprme_closest_gene_")
    try:
        # 4-column BED keyed by POSITION in the index, not by the index label:
        # the key round-trips through closest-features as text, and a positional
        # integer is unambiguous where an arbitrary index label need not be.
        keys = list(combined.index[has_hg38])
        targets = os.path.join(tmpdir, "targets.bed")
        with open(targets, "w") as fh:
            for key_pos, idx in enumerate(keys):
                row = combined.loc[idx]
                start = int(row["hg38_start"])
                fh.write(
                    f"{row['hg38_chr']}\t{start}\t{start + _target_length(row)}\t{key_pos}\n"
                )
        targets_sorted = os.path.join(tmpdir, "targets.sorted.bed")
        with open(targets_sorted, "wb") as dst:
            if subprocess.call(["sort-bed", targets], stdout=dst) != 0:
                raise RuntimeError("sort-bed failed on the target BED")
        genes_sorted = _sorted_plain_bed(gene_annotation_path, tmpdir)
        found = os.path.join(tmpdir, "found.bed")
        with open(found, "wb") as dst:
            code = subprocess.call(
                ["closest-features", "--closest", "--delim", "\t", "--dist",
                 targets_sorted, genes_sorted],
                stdout=dst,
            )
        if code != 0:
            raise RuntimeError("closest-features failed")
        with open(found) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                try:
                    idx = keys[int(fields[3])]
                except (ValueError, IndexError):
                    continue
                # Guard the distance parse the same way resultIntegrator does:
                # a "NA" (no closest feature) or an empty trailing field once
                # crashed the whole integration through float('').
                try:
                    dist = float(fields[-1].strip())
                except ValueError:
                    continue
                attrs = next((f for f in fields if "gene_id=" in f), "")
                name, gid = _parse_gene_attributes(attrs)
                if not (name or gid):
                    continue
                out.at[idx, CLOSEST_GENE_NAME_COL] = name
                out.at[idx, CLOSEST_GENE_ID_COL] = gid
                out.at[idx, CLOSEST_GENE_DIST_COL] = str(dist / 1000)
        return out
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def add_closest_gene_columns(
    combined: pd.DataFrame, gene_annotation_path: Optional[str]
) -> pd.DataFrame:
    """Returns `combined` with the three closest-gene columns added.

    Unchanged when no gene annotation is given -- the columns are then simply
    absent, and the report's existing `_PRESENT_ANN_KINDS` filter keeps
    dropping Gene/Gene_distance_kb exactly as it does today.

    Also applies complete-search's `intergenic` fallback, but only where a
    split `Annotation_GENCODE` column exists AND is empty for that row, so it
    can never overwrite a real overlap label. That column is produced by the
    per-kind annotation split; where it is absent (no --annotation, or a build
    without the split) this is a no-op. "intergenic" is deliberately NOT added
    to the combined `Annotation` string: that string is a list of OVERLAPPING
    features, and a suffix-less token there would be mis-bucketed by the
    per-kind splitter.
    """
    if not gene_annotation_path:
        return combined
    combined = combined.copy()
    frame = closest_gene_frame(combined, gene_annotation_path)
    for col in CLOSEST_GENE_COLS:
        if col in combined.columns:
            combined = combined.drop(columns=[col])
    for col in CLOSEST_GENE_COLS:
        combined[col] = frame[col]
    if GENCODE_SPLIT_COL in combined.columns:
        gencode = combined[GENCODE_SPLIT_COL]
        blank = gencode.isna() | (gencode.astype(str).str.strip() == "")
        dist = pd.to_numeric(combined[CLOSEST_GENE_DIST_COL], errors="coerce")
        combined.loc[blank & dist.notna() & (dist != 0), GENCODE_SPLIT_COL] = INTERGENIC
    return combined
