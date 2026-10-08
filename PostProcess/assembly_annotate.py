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
"""

import os
import sys
from typing import Optional

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
