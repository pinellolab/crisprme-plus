"""ADDITIVE phase-confirmation companion writer (CRISPRme+ dict-less redesign).

Surfaces the per-off-target CONFIRMED-vs-PUTATIVE distinction the observed-haplotype
enumerator computes, the SAFE ADDITIVE way -- exactly like
``population_summary_companion``: a SEPARATE companion TSV next to the output stem,
one row per dict-less VARIANT off-target, keyed by the SAME identity columns the
bestMerge uses (Chromosome, Position, Direction, crRNA, DNA) plus the creating SNP
column, so a downstream/web join is a straight key match.

WHY A COMPANION FILE (not a new bestMerge column): the finalized ``final_line``
addresses its Reference / ref-score-sentinel / tmp_pos_mms tail by NEGATIVE index in
the CFD (target[-3]/target[-4]) and CRISTA (target[-2]/target[-3]) scorers, so
appending a trailing column silently corrupts EVERY variant row's score. And the
downstream positional consumers -- ``adjust_cols.py`` (positional slice) and
``change_headers_bestMerge.py`` (``chunk[new_order]``, a fixed NAME list that SILENTLY
DROPS unknown columns) -- would either desync or drop the flag. A separate joinable
TSV avoids all of that while still surfacing the phase confidence.

GATE: the caller only invokes this on the dict-less branch (``mygt is not None``) with
a non-empty row list, so on a legacy / dict-only install NOTHING is written (byte-
identical legacy behavior). Any error is caught by the caller (log + skip) so the
companion can NEVER break a run.

PHASE SEMANTICS (documented in the header comment of the emitted file):
  * CONFIRMED -- every carrier reached this haplotype via a PHASED, same-phase-set cis
    path (a real cis haplotype).
  * PUTATIVE  -- >= 1 carrier reached it unphased, cross-phase-set, or (PS absent) the
    single-whole-chromosome-block assumption could not confirm cis; the haplotype is
    NEVER dropped, only flagged phase-unconfirmed.

STDLIB ONLY.
"""

from __future__ import annotations

# Identity columns first (join key), then the flag. Mirrors the population-summary
# companion's identity-first layout.
COMPANION_HEADER = (
    "#Chromosome\tPosition\tDirection\tcrRNA\tDNA\tSNP\tPhase_Confirmation"
)


def write_companion(out_path, rows):
    """Write the phase-confirmation companion TSV.

    Args:
      out_path: destination path (``<outputFile>.phase_confirmation.tsv``).
      rows: list of dicts, each with keys Chromosome, Position, Direction, crRNA, DNA,
        SNP, Phase_Confirmation (as accumulated in new_simple_analysis).

    Returns the number of rows written (0 => no file content beyond the header). The
    caller gates on ``rows`` being non-empty, but we still handle the empty case
    gracefully (write just the documented header) so a partial run is self-describing.
    """
    with open(out_path, "w") as fh:
        fh.write(
            "# CRISPRme+ dict-less phase-confirmation companion. One row per variant\n"
            "# off-target (deduped by identity). Join FROM bestMerge/bestCFD/bestCRISTA\n"
            "# by (Chromosome, Position, Direction, crRNA, DNA[, SNP]).\n"
            "#   CONFIRMED = real phased same-phase-set cis haplotype.\n"
            "#   PUTATIVE  = >=1 carrier unphased / cross-phase-set / PS-absent single-\n"
            "#               block; haplotype flagged phase-unconfirmed, never dropped.\n"
        )
        fh.write(COMPANION_HEADER + "\n")
        n = 0
        for rec in rows:
            fh.write(
                "\t".join((
                    str(rec.get("Chromosome", ".")),
                    str(rec.get("Position", ".")),
                    str(rec.get("Direction", ".")),
                    str(rec.get("crRNA", ".")),
                    str(rec.get("DNA", ".")),
                    str(rec.get("SNP", ".")),
                    str(rec.get("Phase_Confirmation", ".")),
                )) + "\n"
            )
            n += 1
    return n
