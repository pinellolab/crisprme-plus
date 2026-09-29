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

from annotation import annotate_target, load_annotation_bed  # noqa: E402

ANNOTATION_COL = "Annotation"
NO_ANNOTATION = "n"  # same "nothing here" marker complete-search's annotation uses

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


def add_annotation_column(
    combined: pd.DataFrame, annotation_path: Optional[str]
) -> pd.DataFrame:
    """Returns `combined` with an `Annotation` column added (or unchanged if
    no annotation file was given)."""
    if not annotation_path:
        return combined
    combined = combined.copy()
    values = annotate_combined(combined, annotation_path)
    if ANNOTATION_COL in combined.columns:
        combined = combined.drop(columns=[ANNOTATION_COL])
    if _INSERT_AFTER in combined.columns:
        pos = list(combined.columns).index(_INSERT_AFTER) + 1
    else:
        pos = len(combined.columns)
    combined.insert(pos, ANNOTATION_COL, values)
    return combined
