"""
This module reconciles CRISPRme off-target predictions across two personal
genome haplotype assemblies (e.g. paternal/maternal), producing one combined
report in hg38 coordinates.

Each haplotype's predicted off-targets are lifted independently to hg38 via
`liftOver` (using that haplotype's own chain file), then joined directly on
hg38 coordinates: found on both haplotypes is homozygous-equivalent, found on
only one is heterozygous-equivalent, and predictions that don't lift at all
are haplotype-non-mappable -- invisible to any reference-genome-based search.

This logic was developed and validated interactively against real HG01255
(HPRC) data before being added here.

Only `*_integrated_results.tsv` is read (2026-09-10 -- previously this
module also read `*_all_results_with_alternative_alignments.tsv` and
outer-merged it in, on the theory that `integrated_results.tsv` alone was
missing real distinct sites; that theory was wrong, see
`load_crisprme_predictions()`'s own docstring for the full derivation).
`integrated_results.tsv` already has exactly one row per real physical
cluster by construction (`merge_contiguous_targets.py`'s own
`retrieve_best_target()` already picks one representative per cluster,
sorted by fewest mismatches+bulges first for the `_(fewest_mm+b)` column
family this module uses throughout).

One real, independent bug is still guarded against here: ~5.7% of loci land
1-3bp off their "true" anchor when a bulge shifts the alignment's
registration -- the *same* physical site reported at a shifted position, a
real artifact even within one already-deduplicated `integrated_results.tsv`
file. `cluster_collapse()` guards against this by reusing the exact greedy
chained-gap clustering algorithm CRISPRme's own `--merge` step uses
(`merge_contiguous_targets.py:531-596`), rather than exact-coordinate
matching, sorted by fewest mismatches+bulges (matching
`retrieve_best_target()`'s own primary criterion, not CFD -- picking a
DIFFERENT criterion than the file's own already-correct choice was this
module's real historical bug, now fixed at both call sites). The `merge_bp`
passed to every function in this module should match the `--merge` value
used for the underlying `complete-search` runs, since this is reusing
CRISPRme's own definition of "one site."
"""

from typing import Dict, List, Optional, Tuple

import glob
import hashlib
import os
import shutil
import subprocess
import time

from datetime import datetime

import pandas as pd

PRED_COLS = [
    "Spacer+PAM", "Chromosome", "Start_coordinate_(fewest_mm+b)",
    "Strand_(fewest_mm+b)", "Aligned_spacer+PAM_(fewest_mm+b)",
    "Aligned_protospacer+PAM_REF_(fewest_mm+b)", "Aligned_protospacer+PAM_ALT_(fewest_mm+b)",
    "Mismatches_(fewest_mm+b)", "Bulges_(fewest_mm+b)", "CFD_score_(fewest_mm+b)",
]
# if more than this fraction of a haplotype's predictions have a chromosome
# missing from the chromAlias file, that's almost certainly a genome/alias
# naming mismatch, not a handful of legitimately-unmappable decoy contigs --
# raise instead of silently folding all of them into "non-mappable"
CHROM_ALIAS_MISMATCH_ERROR_RATIO = 0.5
# same principle, applied to liftOver's own rejection rate (distinct failure
# class from the chromAlias-coverage check above -- e.g. a wrong/swapped
# chain file, not an unrecognized chromosome name)
LIFTOVER_FAILURE_ERROR_RATIO = 0.5


def cluster_collapse(
    df: pd.DataFrame, chrom_col: str, strand_col: str, pos_col: str,
    score_col: str, merge_bp: int, ascending: bool = False,
) -> pd.DataFrame:
    """Collapses rows to one per proximity cluster.

    Uses the same greedy chained-gap algorithm CRISPRme's own `--merge` step
    uses: sort by position, start a new cluster whenever the gap from the
    previous point (same chrom+strand) exceeds `merge_bp` -- clusters chain,
    this is not a fixed window -- then keep the best-scoring row per cluster.

    Args:
        df: Rows to collapse.
        chrom_col: Column name holding the chromosome.
        strand_col: Column name holding the strand.
        pos_col: Column name holding the position to cluster on.
        score_col: Column name to pick the best row per cluster.
        merge_bp: Maximum gap (bp) between consecutive points to stay in the
            same cluster. Should match the `--merge` value used to generate
            the underlying CRISPRme results.
        ascending: False (default) picks the HIGHEST `score_col` per cluster
            (e.g. CFD -- higher is better). True picks the LOWEST (e.g.
            mismatches+bulges -- fewer is better, matching
            `merge_contiguous_targets.py`'s own `_(fewest_mm+b)` criterion).
            Either way NaN sorts last, never picked over a real value.

    Returns:
        One row per cluster, the best-`score_col` row in each.
    """
    df = df.sort_values([chrom_col, strand_col, pos_col]).reset_index(drop=True)
    cluster_ids: List[int] = []
    cluster_id = 0
    prev_chrom = prev_strand = prev_pos = None
    for chrom, strand, pos in zip(df[chrom_col], df[strand_col], df[pos_col]):
        if chrom != prev_chrom or strand != prev_strand or (pos - prev_pos) > merge_bp:
            cluster_id += 1
        cluster_ids.append(cluster_id)
        prev_chrom, prev_strand, prev_pos = chrom, strand, pos
    df = df.assign(_cluster_id=cluster_ids)
    df = df.sort_values(score_col, ascending=ascending, na_position="last")
    df = df.drop_duplicates(subset=["_cluster_id"], keep="first")
    return df.drop(columns=["_cluster_id"]).reset_index(drop=True)


def find_results_prefix(results_dir: str) -> str:
    """Derives a haplotype's CRISPRme results filename prefix by locating its
    `*_integrated_results.tsv` file, rather than reconstructing CRISPRme's
    internal naming convention (guide+PAM+genome+mm+bMax) by hand.

    Args:
        results_dir: A `complete-search` output folder.

    Returns:
        The filename prefix shared by `<prefix>_integrated_results.tsv` and
        `<prefix>_all_results_with_alternative_alignments.tsv`.

    Raises:
        FileNotFoundError: If zero or more than one `*_integrated_results.tsv`
            is found (e.g. a multi-guide run, not yet supported here).
    """
    suffix = "_integrated_results.tsv"
    matches = glob.glob(os.path.join(results_dir, f"*{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one *{suffix} file in {results_dir}, found "
            f"{len(matches)}: {matches}"
        )
    return os.path.basename(matches[0])[: -len(suffix)]


# submit_job_automated_new_multiple_vcfs.sh renames log_error.txt to this,
# unconditionally, as its very last meaningful step (right after echoing
# "JOB END") -- so its presence is the one signal that a run reached its
# true end, not just that some result files happened to get written along
# the way.
LOG_ERROR_NO_CHECK_FILENAME = "log_error_no_check.txt"
# Same filenames pages/pages_utils.py's own PARAMS_FILE/GUIDES_FILE constants
# name -- kept as separate constants here (not a shared import) since
# crisprme.py has no dependency on the pages/ web layer and shouldn't gain
# one just for this.
PARAMS_FILE = ".Params.txt"
GUIDES_FILE = ".guides.txt"


def haplotype_search_complete(results_dir: str) -> bool:
    """Checks whether a haplotype's `complete-search` output represents a
    genuinely finished run, not just one that produced some result files
    before failing partway through.

    `*_integrated_results.tsv` (what `find_results_prefix` looks for) is
    written well before the pipeline's actual end -- a crash during the
    zip/db_creation/mail-send/cleanup steps that follow it would still leave
    that file behind. `log_error.txt` only gets renamed to
    `log_error_no_check.txt` as the pipeline's unconditional last step, so
    requiring both together is the accurate "did this really finish" check,
    used to decide whether `assembly_search` can skip re-running an
    already-completed haplotype search.

    2026-09-10: no longer also requires
    `*_all_results_with_alternative_alignments.tsv` to exist --
    `load_crisprme_predictions` doesn't read that file anymore (see its own
    docstring), so requiring it here would reject a genuinely complete run
    over a file this pipeline no longer needs.

    Args:
        results_dir: A `complete-search` output folder.

    Returns:
        True if the integrated results file and the post-completion log
        rename are both present.
    """
    if not os.path.isdir(results_dir):
        return False
    try:
        find_results_prefix(results_dir)
    except FileNotFoundError:
        return False
    return os.path.isfile(os.path.join(results_dir, LOG_ERROR_NO_CHECK_FILENAME))


def clean_incomplete_haplotype_output(
    output_dir: str, reason: str = "an incomplete previous attempt"
) -> None:
    """Removes a haplotype's stale output directory left behind by a
    previous incomplete/crashed `complete-search` attempt (or one built with
    different search parameters, see `haplotype_params_match`), so a retry
    can actually proceed.

    `complete_search()`'s own `_check_output()` guard refuses to run into a
    non-empty `--output` folder -- exactly right for protecting against
    clobbering a genuinely different run, but it also blocks retrying the
    *same* haplotype search after a crash, since a crash always leaves
    partial files behind (its own error message says as much: "please
    delete ... before running a new CRISPRme search"). Only call this once
    the caller has already confirmed `output_dir` does NOT represent a
    usable, current result -- at that point nothing in it is worth keeping,
    so removing it to unblock the retry is safe.

    No-op if the directory doesn't exist yet (a fresh, first-time run has
    nothing to clean up).

    Args:
        output_dir: The specific haplotype's own output directory (e.g.
            `Results/<name>_paternal`) -- never call this with anything
            broader, since it recursively deletes everything under it.
        reason: Human-readable phrase describing why this directory is
            being removed, logged so the user isn't left wondering (e.g.
            override with "results built with different search parameters"
            when that's why, rather than the default incomplete-run wording).
    """
    if os.path.isdir(output_dir):
        print(f"Detected {reason} at {output_dir}, removing before retry")
        shutil.rmtree(output_dir)


_HAPLOTYPE_PARAM_FLAGS = [
    "--genome", "--guide", "--pam", "--mm", "--bDNA", "--bRNA", "--merge",
    "--max-total-edits",
]
COMMAND_LINE_FILENAME = ".command_line.txt"


def _read_recorded_params(results_dir: str) -> Optional[Dict[str, str]]:
    """Reads the recorded `--genome`/`--guide`/`--pam`/`--mm`/`--bDNA`/
    `--bRNA`/`--merge` values from a previous `complete-search` run's
    `.command_line.txt` -- written unconditionally by every `complete-search`
    invocation (including each haplotype's own, via `_run_haplotype_search`),
    so this needs no new files or directory-naming changes to read.

    Returns:
        A dict of flag -> value for whichever of the flags above were found,
        or `None` if the file is missing or doesn't parse (e.g. a genuinely
        incomplete/crashed run that never got far enough to write it).
    """
    command_line_path = os.path.join(results_dir, COMMAND_LINE_FILENAME)
    if not os.path.isfile(command_line_path):
        return None
    with open(command_line_path) as f:
        content = f.read()
    prefix = "input_command\t"
    if not content.startswith(prefix):
        return None
    tokens = content[len(prefix):].strip().split()
    recorded = {}
    for flag in _HAPLOTYPE_PARAM_FLAGS:
        if flag in tokens:
            idx = tokens.index(flag)
            if idx + 1 < len(tokens):
                recorded[flag] = tokens[idx + 1]
    return recorded


def haplotype_params_match(
    results_dir: str, genome_dir: str, guide_file: str, pam_file: str,
    mm: int, bDNA: int, bRNA: int, merge_bp: int, max_total_edits: int,
) -> bool:
    """Checks whether a haplotype's existing output directory was produced
    with the same search parameters as the current invocation.

    Without this, `assembly_search` would silently reuse a stale result if
    rerun with the same `--output` name but different `--genome-*`/`--guide`/
    `--pam`/`--mm`/`--bDNA`/`--bRNA`/`--merge`/`--max-total-edits` values --
    `haplotype_search_complete` only checks that *a* complete result exists,
    not that it's the *current* one. Reads the existing `.command_line.txt`
    `complete_search()` already writes; no new files, no directory-naming
    changes.

    Args:
        results_dir: The haplotype's output directory to check.
        genome_dir, guide_file, pam_file: Absolute paths, as resolved by
            `assembly_search()`'s own arg-checking (must match exactly).
        mm, bDNA, bRNA, merge_bp: The current invocation's values.
        max_total_edits: The current invocation's resolved total-edits cap
            (an explicit `--max-total-edits` override if the user gave one,
            else `assembly_search()`'s own `mm+bDNA+bRNA` default) -- without
            this, two runs with identical mm/bDNA/bRNA but different
            `--max-total-edits` overrides would be wrongly treated as
            reusable, silently reusing a search built with the wrong cap.

    Returns:
        True only if every recorded parameter matches exactly; False if
        anything differs, or if the recorded parameters can't be read at all
        (treated the same as "don't reuse" -- consistent with how
        `haplotype_search_complete` already treats an unreadable result as
        not-yet-complete).
    """
    recorded = _read_recorded_params(results_dir)
    if recorded is None:
        return False
    current = {
        "--genome": genome_dir, "--guide": guide_file, "--pam": pam_file,
        "--mm": str(mm), "--bDNA": str(bDNA), "--bRNA": str(bRNA),
        "--merge": str(merge_bp), "--max-total-edits": str(max_total_edits),
    }
    return recorded == current


def _read_haplotype_pam_label(haplotype_results_dir: str) -> str:
    """Reads the `Pam` value complete_search() already recorded in a
    haplotype's own `.Params.txt` (2-column `key\\tvalue`, e.g. `Pam\\tNGG`)
    -- reused here instead of re-deriving the PAM motif from the PAM file
    ourselves, which would duplicate complete_search()'s own file-format
    parsing (multi-line rejection, signed offset, etc.) for a value that's
    purely for display. Falls back to "?" if unreadable -- consistent with
    how the web page already shows "?" for a missing Pam key."""
    params_path = os.path.join(haplotype_results_dir, PARAMS_FILE)
    if not os.path.isfile(params_path):
        return "?"
    try:
        with open(params_path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0] == "Pam":
                    return fields[1]
    except OSError:
        pass
    return "?"


def write_combined_params_file(
    combined_dir: str, output_base: str,
    genome_paternal: str, genome_maternal: str,
    chain_paternal: str, chain_maternal: str,
    chrom_alias_paternal: str, chrom_alias_maternal: str,
    paternal_results_dir: str, mm: int, bDNA: int, bRNA: int,
) -> None:
    """Writes the combined job's own `.Params.txt` sidecar.

    `complete_search()` unconditionally writes a `.Params.txt` for every run
    it produces, CLI or web -- that's what lets the web UI recognize and
    render any completed job regardless of how it was launched
    (`index.py`'s `_is_assembly_job()`, `history_page.py`'s
    `retrieve_resultsdirs()`). `assembly_search()` had no equivalent for its
    own `{output_base}_combined` output: that file was only ever written by
    `submit_assembly_search_job()` (pages/main_page.py), the web form's own
    submit handler, before it shells out to `assembly-search` as a
    subprocess. A user running the documented `assembly-search` CLI
    subcommand directly got a combined job directory the web UI couldn't
    render at all -- this closes that gap the same way `complete_search()`
    already handles it.

    Same 3-column (index, key, value) format and key set
    `submit_assembly_search_job()` writes, plus one addition (`Job_start`,
    absent from the web layer's own version): `history_page.py`'s job list
    needs a real timestamp, and assembly-search's own `log.txt` (a raw
    merged stdout/stderr dump, unlike complete_search()'s structured
    per-stage log) has no equivalent line to read one back from.

    Args:
        combined_dir: The `{output_base}_combined` output directory.
        output_base: The `--output` value the run was invoked with.
        genome_paternal, genome_maternal: Absolute paths to the two haplotype
            genome folders.
        chain_paternal, chain_maternal, chrom_alias_paternal,
            chrom_alias_maternal: Absolute paths to the liftOver/chromAlias
            files.
        paternal_results_dir: The paternal haplotype's own (already-complete)
            `complete-search` output dir, read back for its recorded `Pam`
            label -- see `_read_haplotype_pam_label`.
        mm, bDNA, bRNA: The current invocation's search parameters.
    """
    params_lines = [
        ("Genome_type", "assembly"),
        ("Genome_paternal", os.path.basename(genome_paternal.rstrip(os.sep))),
        ("Genome_maternal", os.path.basename(genome_maternal.rstrip(os.sep))),
        ("Chain_paternal", chain_paternal),
        ("Chain_maternal", chain_maternal),
        ("ChromAlias_paternal", chrom_alias_paternal),
        ("ChromAlias_maternal", chrom_alias_maternal),
        ("Pam", _read_haplotype_pam_label(paternal_results_dir)),
        ("Mismatches", str(mm)),
        ("DNA", str(bDNA)),
        ("RNA", str(bRNA)),
        ("Output_base", output_base),
        ("Paternal_dir", f"{output_base}_paternal"),
        ("Maternal_dir", f"{output_base}_maternal"),
        ("Combined_dir", f"{output_base}_combined"),
        ("Job_start", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
    ]
    with open(os.path.join(combined_dir, PARAMS_FILE), "w") as pf:
        for i, (key, value) in enumerate(params_lines, start=1):
            pf.write(f"{i}\t{key}\t{value}\n")


def write_combined_guides_file(combined_dir: str, guide_file: str) -> None:
    """Copies the run's guide file into the combined job's own `.guides.txt`.

    Same gap as `write_combined_params_file()`, confirmed the same way: a
    real (not hand-built) `assembly-search` CLI run leaves the combined
    directory with no `.guides.txt` at all, because -- exactly like
    `.Params.txt` -- `complete_search()` unconditionally copies the guide
    file into every output dir it produces (`crisprme.py`: copies to
    `guides.txt`, then renames to `.guides.txt`), but `assembly_search()`
    had no equivalent for its own combined output; only
    `submit_assembly_search_job()` (the web form's submit handler) ever
    wrote one. `history_page.py`'s `count_guides()` (and any future reader
    of a combined job's guide list) needs this file to exist regardless of
    whether the job was launched from the CLI or the web.

    Args:
        combined_dir: The `{output_base}_combined` output directory.
        guide_file: The `--guide` file path the run was invoked with.
    """
    shutil.copyfile(guide_file, os.path.join(combined_dir, GUIDES_FILE))


def load_chrom_alias(chrom_alias_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Loads a chromAlias file into assembly->ucsc and ucsc->genbank mappings.

    Args:
        chrom_alias_file: Path to an HPRC-style `*.chromAlias.txt` file
            (tab-separated, columns `# assembly`, `ucsc`, `genbank`).

    Returns:
        `(assembly_to_ucsc, ucsc_to_genbank)`.
    """
    chrom_alias_df = pd.read_csv(chrom_alias_file, sep="\t")
    assembly_to_ucsc = dict(zip(chrom_alias_df["# assembly"], chrom_alias_df["ucsc"]))
    ucsc_to_genbank = dict(zip(chrom_alias_df["ucsc"], chrom_alias_df["genbank"]))
    return assembly_to_ucsc, ucsc_to_genbank


def load_crisprme_predictions(
    results_dir: str, prefix: str, merge_bp: int, cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """Loads one haplotype's CRISPRme output into one deduplicated,
    one-row-per-genomic-locus table.

    Reads ONLY `*_integrated_results.tsv` -- CRISPRme's own merge step
    (`merge_contiguous_targets.py`'s `retrieve_best_target()`) already picks
    exactly one representative row per real physical cluster for the
    `_(fewest_mm+b)` column family (sorted by fewest mismatches+bulges
    FIRST, CFD isn't even in that primary key), so `integrated_results.tsv`
    already has exactly one row per real site by construction.

    2026-09-10, decided after real re-derivation (Manuel agreed; Luca
    confirmed 2026-09-10): this function used to ALSO read
    `*_all_results_with_alternative_alignments.tsv` and outer-merge it in,
    on the theory that `integrated_results.tsv` alone was missing real
    distinct sites. That theory was wrong -- the union's own "duplicates"
    (up to 233 rows at one locus) were alt-file-only artifacts (multiple
    alternative alignments for ONE already-counted site), not additional
    real off-targets. Re-clustering the union by CFD (this function's own
    old default) doesn't recover anything real; it just occasionally swaps
    in a less-consistent representative than `integrated_results.tsv`'s own
    already-correct fewest-mm+b choice.

    `cluster_collapse()` is still applied (not redundant): it's real,
    independent protection against a bulge-registration-drift artifact
    documented at the top of this module (~5.7% of loci land 1-3bp off
    their "true" anchor when a bulge shifts the alignment's registration) --
    a real phenomenon `integrated_results.tsv` on its own doesn't rule out,
    unrelated to the alt-file question above. Sorts by fewest
    mismatches+bulges (ascending), matching the SAME `_(fewest_mm+b)`
    criterion `retrieve_best_target()` already used to build this file in
    the first place -- not CFD, which was this function's own prior (now
    fixed) inconsistency, not `retrieve_best_target()`'s.

    Args:
        results_dir: A `complete-search` output folder.
        prefix: Filename prefix (from `find_results_prefix`).
        merge_bp: Passed to `cluster_collapse` -- should match the `--merge`
            value used for this haplotype's `complete-search` run.
        cols: Result columns to keep. Defaults to `PRED_COLS`.

    Returns:
        One row per genomic locus, with a stable `off_target_id` column.
    """
    if cols is None:
        cols = PRED_COLS

    integrated_file = os.path.join(results_dir, f"{prefix}_integrated_results.tsv")
    combined = pd.read_csv(integrated_file, sep="\t", usecols=cols)

    mmb = (
        pd.to_numeric(combined["Mismatches_(fewest_mm+b)"], errors="coerce")
        + pd.to_numeric(combined["Bulges_(fewest_mm+b)"], errors="coerce")
    )
    combined = cluster_collapse(
        combined.assign(_mmb=mmb), "Chromosome", "Strand_(fewest_mm+b)",
        "Start_coordinate_(fewest_mm+b)", "_mmb", merge_bp, ascending=True,
    ).drop(columns=["_mmb"])
    combined["off_target_id"] = combined.index.astype(str)
    return combined


def check_chrom_alias_coverage(
    preds: pd.DataFrame, ucsc_to_genbank: Dict[str, str], name: str
) -> None:
    """Sanity-checks that the chromAlias file actually names this haplotype's
    chromosomes, before treating any mismatch as ordinary haplotype-private
    biology.

    A handful of predictions on contigs the chromAlias file doesn't cover
    (e.g. small unplaced/decoy scaffolds) is expected and gets folded into
    the non-mappable count downstream. But if a *large* fraction of
    predictions have an unrecognized chromosome, that's a strong signal the
    genome FASTA's chromosome naming doesn't actually match the supplied
    chromAlias file (wrong file, wrong assembly, or a naming convention the
    file doesn't cover) -- worth failing clearly and immediately rather than
    silently reporting a mostly-meaningless "almost everything is
    haplotype-private" result.

    Args:
        preds: This haplotype's predictions (must have a `Chromosome` column).
        ucsc_to_genbank: This haplotype's chromAlias `ucsc`->`genbank` mapping.
        name: Haplotype name, for the error message.

    Raises:
        ValueError: If the unmatched-row fraction exceeds
            `CHROM_ALIAS_MISMATCH_ERROR_RATIO`.
    """
    total = len(preds)
    if total == 0:
        return
    known = set(ucsc_to_genbank)
    unmatched = preds.loc[~preds["Chromosome"].isin(known)]
    ratio = len(unmatched) / total
    if ratio > CHROM_ALIAS_MISMATCH_ERROR_RATIO:
        examples = sorted(unmatched["Chromosome"].unique())[:5]
        raise ValueError(
            f"{name}: {len(unmatched)}/{total} predictions ({ratio:.0%}) have a "
            f"chromosome not found in the {name} chromAlias file's 'ucsc' column "
            f"(e.g. {examples}) -- this usually means the genome FASTA's "
            f"chromosome naming doesn't match the --chrom-alias-{name} file "
            "provided, not that these are genuinely unmappable contigs. Both "
            "haplotype searches already completed successfully; only "
            "reconciliation is affected by this. Fix the chromAlias file or "
            "genome naming, then re-run assembly-search -- already-completed "
            "haplotype searches will be reused, not re-run."
        )


def check_liftover_failure_rate(attempted_count: int, unlifted_count: int, name: str) -> None:
    """Sanity-checks that `liftOver` itself isn't rejecting a suspiciously
    high fraction of inputs, before treating a large non-mappable count as
    ordinary haplotype-private biology.

    Distinct from `check_chrom_alias_coverage`, which only catches
    predictions whose chromosome isn't even recognized by the chromAlias
    file, checked *before* liftOver ever runs. This catches a different
    failure class: liftOver runs and genuinely rejects most of its input,
    which can happen with a wrong chain file (e.g. the paternal/maternal
    `--chain-*` files swapped, or a stale/mismatched chain file) -- same
    "fail clearly instead of silently reporting a mostly-meaningless result"
    principle as the chromAlias check, applied to the post-liftover signal.

    Args:
        attempted_count: Predictions actually handed to liftOver (i.e.
            excluding ones already dropped for lacking a chromAlias entry --
            those are a distinct, already-guarded failure mode).
        unlifted_count: How many of those liftOver itself rejected.
        name: Haplotype name, for the error message.

    Raises:
        ValueError: If the rejection rate exceeds `LIFTOVER_FAILURE_ERROR_RATIO`.
    """
    if attempted_count == 0:
        return
    ratio = unlifted_count / attempted_count
    if ratio > LIFTOVER_FAILURE_ERROR_RATIO:
        raise ValueError(
            f"{name}: liftOver rejected {unlifted_count}/{attempted_count} "
            f"predictions ({ratio:.0%}) -- this usually means the wrong chain "
            f"file was supplied for --chain-{name} (e.g. paternal/maternal "
            "chain files swapped, or a stale/mismatched chain file), not that "
            "these are genuinely haplotype-private sites. Both haplotype "
            "searches already completed successfully; only reconciliation is "
            "affected by this. Fix the chain file, then re-run assembly-search "
            "-- already-completed haplotype searches will be reused, not "
            "re-run."
        )


def build_offtarget_bed(
    preds: pd.DataFrame, ucsc_to_genbank: Dict[str, str], bed_path: str
) -> Tuple[str, set]:
    """Writes a BED file of one haplotype's predicted off-target coordinates,
    keyed by each row's `off_target_id`, in the chain file's chromosome
    naming (GenBank accessions, per the HPRC chromAlias convention).

    Returns:
        `(bed_path, dropped_ids)` -- `dropped_ids` are the `off_target_id`s
        excluded because their chromosome has no chromAlias mapping. These
        are just as unmappable-to-hg38 as a genuine liftOver failure (there's
        no way to look them up in the chain file at all), so callers should
        fold them into the same non-mappable accounting as
        `load_unlifted_ids`, not just discard them.
    """
    bed = preds[["Chromosome", "Start_coordinate_(fewest_mm+b)", "off_target_id"]].copy()
    bed = bed.rename(columns={"Start_coordinate_(fewest_mm+b)": "chromEnd"})
    bed["chromStart"] = bed["chromEnd"] - 1
    bed["chrom"] = bed["Chromosome"].map(ucsc_to_genbank)
    dropped_ids = set(bed.loc[bed["chrom"].isna(), "off_target_id"])
    bed = bed.dropna(subset=["chrom"])
    bed = bed[["chrom", "chromStart", "chromEnd", "off_target_id"]]
    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    return bed_path, dropped_ids


def check_liftover_available() -> None:
    """Verifies the `liftOver` binary is on PATH before any expensive work
    starts.

    `liftOver` (UCSC tool, bioconda package `ucsc-liftover`) is a required
    dependency for `assembly-search` -- bundled into the Docker image
    alongside crispritz/crisprme, but not yet part of CRISPRme's own bioconda
    recipe for the plain `mamba install crisprme` path (a separate,
    not-yet-done follow-up). Checking this before launching two multi-hour
    haplotype searches means a missing dependency fails immediately and
    clearly instead of only surfacing at the very last reconciliation step.

    Raises:
        RuntimeError: If `liftOver` isn't found on PATH.
    """
    if shutil.which("liftOver") is None:
        raise RuntimeError(
            "liftOver (UCSC tool) not found on PATH -- required for "
            "assembly-search. Install it with `mamba install -c bioconda "
            "ucsc-liftover` (or `conda install`), or use a CRISPRme Docker "
            "image built after this dependency was added."
        )


def run_liftover(bed_path: str, chain_file: str, mapped_path: str, unmapped_path: str) -> Tuple[str, str]:
    """Runs UCSC `liftOver` on a BED file. Assumes `liftOver` is on PATH --
    see `check_liftover_available`.

    Raises:
        RuntimeError: If `liftOver` exits non-zero.
    """
    cmd = f"liftOver {bed_path} {chain_file} {mapped_path} {unmapped_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"liftOver failed:\n{result.stderr}")
    return mapped_path, unmapped_path


def load_lifted_bed(mapped_path: str) -> pd.DataFrame:
    return pd.read_csv(
        mapped_path, sep="\t", header=None,
        names=["hg38_chr", "hg38_start", "hg38_end", "off_target_id"],
        dtype={"off_target_id": str},
    )


def load_unlifted_ids(unmapped_path: str) -> set:
    """liftOver writes failed features back out in BED format, with
    '#'-prefixed comment lines explaining why each one failed to map."""
    ids = set()
    with open(unmapped_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.strip().split("\t")[3])
    return ids


# =============================================================================
# Direct haplotype-vs-haplotype reconciliation (impg/minimap2)
#
# Resolves a subset of each haplotype's `_non_mappable` (no hg38 equivalent)
# predictions by aligning the two haplotype genomes directly to each other,
# bypassing hg38 entirely. Validated against real HG01255 (HPRC) data before
# being added here; the findings that shaped the design are recorded inline
# below: (1) impg's `-r` query silently breaks if both sides of a PAF share
# identical sequence names, fixed by prefixing each side's sequence names
# before indexing; (2) a locus with more than one candidate alignment hit is
# the direct, checkable signature of a repeat-family region -- collapsing on
# such a hit risks pairing the wrong copy (confirmed with a real example, a
# chr9 pericentromeric repeat family with 3 copies on one haplotype and 6 on
# the other), so `resolve_haplotype_private` requires exactly one hit before
# treating a pair as confirmed; (3) that single-hit check is directional
# (queries FROM one haplotype INTO the other) and genuinely asymmetric --
# a real genome-wide measurement found a locus can pass cleanly in one
# direction while the reverse query for the same physical site returns
# multiple candidates, almost always with tied identity scores (no
# per-hit signal to break the tie). `reconcile_haplotypes` therefore calls
# `resolve_haplotype_private_bidirectional`, not `resolve_haplotype_private`
# directly; (4) minimap2 itself is not symmetric between query and target,
# so even a single alignment's own choice of which genome to index is an
# arbitrary decision that measurably changes the result (confirmed
# 2026-09-02: swapping the fixed orientation changed the reconciled count
# by 65 pairs on real data, with no principled way to prefer one fixed
# choice over the other) -- `reconcile_haplotypes` therefore builds BOTH
# orientations (`build_or_reuse_haplotype_alignments_both_orientations`)
# and routes each query direction to its own matched alignment, rather
# than picking one orientation for both. See
# `resolve_haplotype_private_bidirectional`'s docstring for the real
# numbers and the comparison against a symmetric wfmash alignment this
# choice is based on.
# =============================================================================

HAPLOTYPE_ALIGNMENTS_DIRNAME = "haplotype_alignments"
ALIGNMENT_PARAMS_FILENAME = ".alignment_genomes.txt"


def _alignment_paf_filename(name_target: str) -> str:
    """Filename for the alignment PAF that has `name_target` as minimap2's
    indexed target side, within one (order-independent) genome-pair cache
    directory.

    Both orientations of a genome pair are real, different alignments --
    minimap2 is not symmetric between query and target (confirmed directly,
    2026-09-02: swapping which genome is indexed changes which candidate
    hits exist at all, especially in repeat regions -- see
    `resolve_haplotype_private_bidirectional`'s docstring). One cache
    directory therefore holds up to two PAF files, one per target choice,
    distinguished by this filename rather than by a separate cache
    directory per orientation -- keeps `_alignment_cache_key` itself
    order-independent (one cache entry per unordered genome pair, matching
    its own docstring) while still caching both real alignments."""
    return f"alignment.target_{name_target}.paf"


def check_impg_available() -> None:
    """Verifies `minimap2` and `impg` are both on PATH before any expensive
    work starts -- mirrors `check_liftover_available` exactly.

    Both are required unconditionally for assembly-search's direct
    haplotype-vs-haplotype reconciliation (unlike `liftOver`, which is only
    needed when a chain file is actually supplied).

    Raises:
        RuntimeError: If either `minimap2` or `impg` isn't found on PATH.
    """
    missing = [tool for tool in ("minimap2", "impg") if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(
            f"{' and '.join(missing)} not found on PATH -- required for "
            "assembly-search's direct haplotype-vs-haplotype reconciliation. "
            "Install with `mamba install -c bioconda minimap2 impg=0.5.0` "
            "(pinned: this module depends on impg 0.5.x-specific CLI/output "
            "behavior), or use a CRISPRme Docker image built after this "
            "dependency was added."
        )


def _alignment_cache_key(genome_a: str, genome_b: str) -> str:
    """Short, order-independent identity for a genome pair, used as the
    cache directory name. Genome folder paths can't be embedded directly in
    a directory name the way `genome_library`'s short `<pam>_<budget>_<ref>`
    identifiers are, so a hash stands in for the name; the literal paths are
    still recorded in a sidecar file (`_read_alignment_cache_params`) so a
    human can audit which genomes a cache entry came from."""
    resolved = sorted([os.path.abspath(genome_a), os.path.abspath(genome_b)])
    return hashlib.sha256("\n".join(resolved).encode()).hexdigest()[:16]


def _read_alignment_cache_params(cache_dir: str) -> Optional[List[str]]:
    """Reads the two literal genome paths recorded for a cache entry, or
    `None` if the sidecar file is missing/malformed."""
    path = os.path.join(cache_dir, ALIGNMENT_PARAMS_FILENAME)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    if len(lines) != 2:
        return None
    return lines


def _haplotype_alignment_cache_valid(
    cache_dir: str, genome_a: str, genome_b: str, name_target: str,
) -> bool:
    """Whether `cache_dir` already holds a usable, current alignment for
    exactly this genome pair, with `name_target` as the indexed target --
    checks the recorded genome paths actually match (not just that *some*
    cache entry exists at this hash, avoiding the directory-existence-only
    weakness `genome_library`'s own reuse check was flagged for) and that
    the PAF + impg index for this specific orientation are both present.
    """
    recorded = _read_alignment_cache_params(cache_dir)
    if recorded is None:
        return False
    current = sorted([os.path.abspath(genome_a), os.path.abspath(genome_b)])
    if recorded != current:
        return False
    paf_path = os.path.join(cache_dir, _alignment_paf_filename(name_target))
    return os.path.isfile(paf_path) and os.path.isfile(paf_path + ".impg")


def _concat_genome_fasta(genome_dir: str, out_path: str) -> str:
    """Concatenates every `*.fa` file in a haplotype genome folder (one file
    per chromosome/scaffold, CRISPRme's own `--genome-*` convention) into a
    single whole-genome FASTA for `minimap2`."""
    fasta_files = sorted(glob.glob(os.path.join(genome_dir, "*.fa")))
    if not fasta_files:
        raise RuntimeError(f"No .fa files found in genome directory: {genome_dir}")
    with open(out_path, "wb") as out:
        for fasta_file in fasta_files:
            with open(fasta_file, "rb") as f:
                shutil.copyfileobj(f, out)
    return out_path


def _rename_paf_sequences(raw_paf: str, out_paf: str, name_a: str, name_b: str) -> str:
    """Prefixes each side's sequence names (PAF columns 1 and 6) with its
    haplotype name -- fixes a real impg identical-sequence-name query bug:
    both haplotype genomes use bare chromosome names (e.g. `chr19` on both
    sides), which silently breaks impg's `-r` name-based lookup unless the
    two sides are made unique first."""
    with open(raw_paf) as src, open(out_paf, "w") as dst:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            fields[0] = f"{name_a}_{fields[0]}"
            fields[5] = f"{name_b}_{fields[5]}"
            dst.write("\t".join(fields) + "\n")
    return out_paf


def build_or_reuse_haplotype_alignment(
    genome_a: str, genome_b: str, name_a: str, name_b: str,
    cache_root: str, threads: int = 4,
) -> str:
    """Returns the path to a ready-to-query, impg-indexed PAF aligning
    `genome_a` (query) against `genome_b` (target), building it on a cache
    miss (`minimap2 -cx asm5 --cs` + sequence-name prefixing + `impg
    index`) and reusing it on a cache hit.

    `name_a`/`name_b` (e.g. `"paternal"`/`"maternal"`) prefix each side's
    sequence names in the renamed PAF -- callers must use the same names
    again when querying (see `query_haplotype_alignment`).

    Concurrency: a lockfile (`<cache_dir>.lock`, atomic `os.mkdir`) guards
    against two processes building the same pair simultaneously, mirroring
    `genome_library`'s own lockfile pattern in
    `submit_job_automated_new_multiple_vcfs.sh`.

    Raises:
        RuntimeError: If `minimap2` or `impg index` exits non-zero.
    """
    os.makedirs(cache_root, exist_ok=True)
    cache_dir = os.path.join(cache_root, _alignment_cache_key(genome_a, genome_b))
    paf_path = os.path.join(cache_dir, _alignment_paf_filename(name_b))

    if _haplotype_alignment_cache_valid(cache_dir, genome_a, genome_b, name_b):
        return paf_path

    # lock is keyed on (cache_dir, target name) -- building the OTHER
    # orientation for the same genome pair (see
    # `build_or_reuse_haplotype_alignments_both_orientations`) must not
    # block on this one's lock, since they're independent minimap2 runs
    # writing different files in the same shared cache_dir.
    lock_path = f"{cache_dir}.target_{name_b}.lock"
    while True:
        try:
            os.mkdir(lock_path)
            break
        except FileExistsError:
            time.sleep(5)
    try:
        # re-check now that the lock is held: another process may have
        # finished building this exact pair while we were waiting
        if _haplotype_alignment_cache_valid(cache_dir, genome_a, genome_b, name_b):
            return paf_path

        os.makedirs(cache_dir, exist_ok=True)
        fasta_a = _concat_genome_fasta(genome_a, os.path.join(cache_dir, f"genome_a.{name_a}.fa"))
        fasta_b = _concat_genome_fasta(genome_b, os.path.join(cache_dir, f"genome_b.{name_b}.fa"))

        raw_paf = os.path.join(cache_dir, f"alignment.target_{name_b}.raw.paf")
        cmd = f"minimap2 -cx asm5 --cs -t {threads} {fasta_b} {fasta_a} > {raw_paf}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed:\n{result.stderr}")

        _rename_paf_sequences(raw_paf, paf_path, name_a, name_b)
        os.remove(raw_paf)
        os.remove(fasta_a)
        os.remove(fasta_b)

        result = subprocess.run(
            ["impg", "index", "-a", paf_path, "-t", str(threads)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"impg index failed:\n{result.stderr}")

        with open(os.path.join(cache_dir, ALIGNMENT_PARAMS_FILENAME), "w") as f:
            for path in sorted([os.path.abspath(genome_a), os.path.abspath(genome_b)]):
                f.write(path + "\n")
    finally:
        os.rmdir(lock_path)
    return paf_path


def build_or_reuse_haplotype_alignments_both_orientations(
    genome_a: str, genome_b: str, name_a: str, name_b: str,
    cache_root: str, threads: int = 4,
) -> Tuple[str, str]:
    """Builds (or reuses) BOTH minimap2 orientations for one genome pair --
    `genome_a` as target / `genome_b` as target -- and returns the two PAF
    paths as `(paf_a_target, paf_b_target)`.

    Why both: minimap2 is not symmetric between query and target for
    assembly-to-assembly alignment (confirmed directly, 2026-09-02, on real
    HG01255 data -- swapping which genome is indexed changed which
    candidate hits exist at all in repeat regions, and measurably changed
    the final reconciled count in both directions). A single orientation
    means one query direction gets the alignment best-suited to detect
    ambiguity in the OTHER genome, and the other direction doesn't. Real,
    measured comparison against multiple independent estimates (a second
    minimap2 orientation, and a symmetric wfmash alignment) found this
    "matched" approach -- pairing each query direction with the alignment
    where the genome being searched INTO was the indexed target -- captures
    close to as many real matches as the larger single-orientation
    alternative while disagreeing with the best available estimate
    noticeably less often (Jaccard 0.505 vs 0.466-0.470 for either single
    orientation alone). See `resolve_haplotype_private_bidirectional` for
    how the two returned PAFs get routed to their matched query direction.

    Cost: builds minimap2 twice (once per orientation) instead of once --
    roughly doubles the one-time, cached-per-genome-pair alignment build
    cost. Each orientation is independently lockfile-guarded so a second
    process building the same pair can reuse whichever orientation is
    already done rather than waiting on both.
    """
    paf_a_target = build_or_reuse_haplotype_alignment(
        genome_b, genome_a, name_b, name_a, cache_root, threads=threads,
    )
    paf_b_target = build_or_reuse_haplotype_alignment(
        genome_a, genome_b, name_a, name_b, cache_root, threads=threads,
    )
    return paf_a_target, paf_b_target


def query_haplotype_alignment(
    paf_path: str, from_name: str, chrom: str, start: int, end: int,
) -> List[Tuple[str, int, int]]:
    """Queries a built alignment (see `build_or_reuse_haplotype_alignment`)
    for the coordinate(s) on the *other* haplotype corresponding to a native
    `(chrom, start, end)` locus on the `from_name` haplotype.

    Returns:
        A list of `(chrom, start, end)` hits on the other haplotype, with
        the haplotype-name prefix stripped back off. Empty if no hit.
        Zero, one, or more than one hit is meaningful to the caller: more
        than one is the signature of a repeat-family region (see this
        section's module-level note) -- deliberately returned as-is
        (not collapsed to "the best hit") so `resolve_haplotype_private`
        can apply its own single-hit requirement.
    """
    # impg 0.5.0 rejects query ranges below 101bp ("below minimum ... Lower
    # --min-transitive-len or use a longer range") even for a plain,
    # non-transitive `-r` query -- confirmed directly against a real index.
    # Off-target windows here are ~1bp; pad symmetrically around the
    # midpoint to clear the minimum rather than pass a fragile CLI override.
    # The projected hit's own midpoint (what callers actually use) is
    # unaffected as long as the padded window stays inside one contiguous
    # alignment block, true at this scale for anything but a query sitting
    # exactly on a block boundary.
    MIN_QUERY_LEN = 101
    midpoint = (start + end) // 2
    if end - start < MIN_QUERY_LEN:
        half = MIN_QUERY_LEN // 2 + 1
        # clamp to 0: a locus within `half` bp of a contig's own start would
        # otherwise produce a negative coordinate, which impg rejects outright
        # (a real query, not a hypothetical -- short contigs/scaffolds are
        # common in personal assemblies).
        start, end = max(0, midpoint - half), midpoint + half

    seq_name = f"{from_name}_{chrom}"
    cmd = [
        "impg", "query", "-a", paf_path,
        "-r", f"{seq_name}:{start}-{end}",
        "-d", "0", "-o", "bedpe",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # A haplotype-A sequence with no alignment data at all (e.g. a short
        # unplaced scaffold minimap2 never emitted any PAF record for
        # against haplotype B) isn't a genuine query failure -- it's a real,
        # valid "no counterpart found" outcome, not a bug. Confirmed against
        # a real genome-wide query (a real HG01255 `chrUn_*` scaffold
        # triggered exactly this). Anything else stays a hard failure.
        if "not found in index" in result.stderr:
            return []
        raise RuntimeError(f"impg query failed:\n{result.stderr}")

    hits = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        target_chrom, target_start, target_end = fields[0], int(fields[1]), int(fields[2])
        stripped_chrom = target_chrom
        if stripped_chrom.startswith(f"{from_name}_"):
            # bidirectional impg queries can echo the query's own side back;
            # only real cross-haplotype hits are meaningful here
            continue
        idx = stripped_chrom.find("_")
        if idx != -1:
            stripped_chrom = stripped_chrom[idx + 1:]
        hits.append((stripped_chrom, target_start, target_end))
    return hits


def resolve_haplotype_private(
    unlifted_ids: Dict[str, set],
    predictions: Dict[str, pd.DataFrame],
    paf_path: str,
    names: List[str],
    merge_bp: int,
) -> Tuple[pd.DataFrame, Dict[str, set]]:
    """Resolves a subset of each haplotype's non-mappable (haplotype-private)
    predictions by direct alignment, per the module-level note above.

    For each haplotype-A-private off-target, looks up its native coordinate
    and queries the direct alignment for haplotype B's corresponding
    coordinate. Collapses a pair into `both_haplotype_private` only when
    BOTH hold:

    1. Reciprocal confirmation: the projected coordinate lands within
       `merge_bp` of an off-target haplotype B's own search *also*
       independently called, and that off-target is *also* in
       `unlifted_ids[B]` -- the same "both haplotypes independently found
       this, neither maps to hg38" shape as the `both` category, for the
       private bucket.
    2. Single-hit only: the direct-alignment query for the haplotype-A
       locus returns exactly one candidate hit. More than one hit is the
       signature of a repeat-family region (see module-level note) --
       left unresolved (stays in `<name>_non_mappable`) rather than
       collapsed on a guess.

    Args:
        unlifted_ids: Per-haplotype sets of haplotype-private `off_target_id`s
            (from `reconcile_haplotypes`'s own accounting).
        predictions: Per-haplotype prediction DataFrames (same ones already
            loaded in `reconcile_haplotypes`), used to look up each
            haplotype-private id's native (non-hg38) coordinate.
        paf_path: Path returned by `build_or_reuse_haplotype_alignment`.
        names: The two haplotype names, in the same order used to build the
            alignment (`name_a`, `name_b`).
        merge_bp: Same locus-clustering tolerance used everywhere else in
            this module.

    Returns:
        `(resolved_rows, absorbed_ids)` -- `resolved_rows` has one row per
        confirmed pair (native + haplotype id for both sides), and
        `absorbed_ids` is `{name: set of ids removed from that haplotype's
        non_mappable count}` (every id in a confirmed pair, both sides).
    """
    a, b = names[0], names[1]
    preds_by_id = {name: predictions[name].set_index("off_target_id") for name in names}

    # index haplotype B's own private predictions by (chrom, start) for a
    # fast "did B's own search independently call something near here" check
    b_private = preds_by_id[b].loc[preds_by_id[b].index.isin(unlifted_ids[b])]
    b_private_by_chrom: Dict[str, List[Tuple[int, str]]] = {}
    for b_id, row in b_private.iterrows():
        b_private_by_chrom.setdefault(row["Chromosome"], []).append(
            (row["Start_coordinate_(fewest_mm+b)"], b_id)
        )

    rows = []
    absorbed = {a: set(), b: set()}
    for a_id in sorted(unlifted_ids[a]):
        a_row = preds_by_id[a].loc[a_id]
        a_chrom = a_row["Chromosome"]
        a_pos = int(a_row["Start_coordinate_(fewest_mm+b)"])

        hits = query_haplotype_alignment(paf_path, a, a_chrom, a_pos - 1, a_pos)
        if len(hits) != 1:
            continue  # zero hits (no counterpart) or multiple (ambiguous/repeat-family): leave unresolved

        hit_chrom, hit_start, hit_end = hits[0]
        hit_mid = (hit_start + hit_end) // 2

        candidates = b_private_by_chrom.get(hit_chrom, [])
        # pick the NEAREST candidate within merge_bp, not the first one found in
        # iteration order -- iteration order is incidental (whatever order
        # b_private.iterrows() produced), not meaningful, so breaking on the
        # first in-range candidate picked an arbitrary one whenever more than
        # one fell within tolerance of the same hit. Measured directly against
        # real HG01255 data (2026-09-02): this exact scenario never actually
        # triggered on that dataset (0 of 1,330 single-hit queries had more
        # than one in-range candidate), so it wasn't the cause of any
        # disagreement observed there -- this is a defensive correctness fix
        # for other genomes/guides with denser repeat clusters, not a fix for
        # a bug already shown to explain real data.
        in_range = [(abs(b_pos - hit_mid), b_id) for b_pos, b_id in candidates
                    if abs(b_pos - hit_mid) <= merge_bp]
        if not in_range:
            continue  # no reciprocal confirmation: leave unresolved
        match_id = min(in_range)[1]

        rows.append({
            f"off_target_id_{a}": a_id, f"off_target_id_{b}": match_id,
            "Chromosome": a_chrom, "Start_coordinate": a_pos,
            "origin": "both_haplotype_private",
        })
        absorbed[a].add(a_id)
        absorbed[b].add(match_id)

    return pd.DataFrame(rows), absorbed


def resolve_haplotype_private_bidirectional(
    unlifted_ids: Dict[str, set],
    predictions: Dict[str, pd.DataFrame],
    paf_path_target_a: str,
    paf_path_target_b: str,
    names: List[str],
    merge_bp: int,
) -> Tuple[pd.DataFrame, Dict[str, set]]:
    """Runs `resolve_haplotype_private` in BOTH directions, each against its
    own MATCHED alignment, and keeps only the pairs both directions
    independently confirm -- the intersection, not the union.

    `resolve_haplotype_private` only ever queries FROM the first name in
    its own `names` INTO the second, so a pair can pass in one direction
    (single hit, reciprocally confirmed) while the reverse query for the
    exact same physical locus is genuinely ambiguous (multiple candidate
    hits) -- real, measured directional asymmetry, not a hypothetical edge
    case.

    Why TWO alignments, not one: minimap2 is not symmetric between query
    and target for assembly-to-assembly alignment -- confirmed directly
    (2026-09-02, real HG01255 data): building the SAME genome pair with
    target and query swapped changes which candidate hits exist at all,
    especially in repeat regions, and measurably changes the final
    reconciled count (352 with one fixed orientation vs. up to 417 with
    the other -- neither more "correct" than the other, just two arbitrary
    choices that disagree by 65 pairs with no principled way to prefer
    one). A query FROM haplotype X INTO haplotype Y is best served by the
    alignment where Y was minimap2's indexed TARGET, since that's what
    determines whether Y has look-alike copies at all -- so this function
    takes two alignments, each already built with the matched target, and
    routes each query direction to its own:
      - querying a->b (into b) uses `paf_path_target_b` (b was the target)
      - querying b->a (into a) uses `paf_path_target_a` (a was the target)
    Measured against a symmetric wfmash alignment (the closest thing to a
    trustworthy independent estimate available), this matched-alignment
    approach agrees more closely (Jaccard 0.505) than either single fixed
    orientation alone (0.466 / 0.470) -- see
    `build_or_reuse_haplotype_alignments_both_orientations` for the full
    real-data comparison this is based on.

    Why intersection, not union, within a single alignment's own two
    directions: characterizing real disagreement pairs found every one is
    a genuine multi-hit case in the direction that failed (never a plain
    zero-hit or an out-of-tolerance near-miss), and most of those are a
    dead tie -- every candidate hit scores identical alignment identity,
    so there is no available signal to prefer one candidate over another.
    Trusting whichever direction happened to query cleanly, in a genuine
    multi-copy repeat, is exactly the false-merge risk the single-hit
    filter in `resolve_haplotype_private` exists to prevent in the first
    place -- requiring both directions to agree closes that gap. This does
    leave some real pairs unresolved rather than merged, which is the
    intended, conservative tradeoff: favor never merging two genuinely
    distinct sites over recovering every possible match.

    Args:
        paf_path_target_a: Alignment PAF with `names[0]` as minimap2's
            indexed target (used for the `names[1] -> names[0]` query).
        paf_path_target_b: Alignment PAF with `names[1]` as minimap2's
            indexed target (used for the `names[0] -> names[1]` query).
        Other args: same as `resolve_haplotype_private`.

    Returns:
        Same shape as `resolve_haplotype_private`: `(resolved_rows,
        absorbed_ids)`, restricted to bidirectionally-confirmed pairs only.
    """
    a, b = names[0], names[1]
    resolved_fwd, _ = resolve_haplotype_private(unlifted_ids, predictions, paf_path_target_b, [a, b], merge_bp)
    resolved_rev, _ = resolve_haplotype_private(unlifted_ids, predictions, paf_path_target_a, [b, a], merge_bp)

    # both calls use the real haplotype-name strings (not positional "a"/
    # "b") for their column names, so both DataFrames share the same
    # column names regardless of query order -- no reordering needed.
    col_a, col_b = f"off_target_id_{a}", f"off_target_id_{b}"

    def _pairset(df: pd.DataFrame) -> set:
        if len(df) == 0:
            return set()
        return set(zip(df[col_a], df[col_b]))

    confirmed = _pairset(resolved_fwd) & _pairset(resolved_rev)

    # Per-pair native-locus lookup from each direction's OWN result.
    # `resolved_fwd` was queried as [a, b], so its internal "Chromosome"/
    # "Start_coordinate" (see `resolve_haplotype_private`) is haplotype `a`'s
    # own native locus; `resolved_rev` was queried as [b, a], so ITS
    # "Chromosome"/"Start_coordinate" is really haplotype `b`'s native locus.
    # Neither haplotype has an hg38 coordinate here by definition (that's why
    # the site is haplotype-private), so each side's own native coordinate is
    # the most specific location this row can carry -- same convention used
    # for non-mappable sites elsewhere in this pipeline.
    #
    # Real bug, fixed 2026-09-10: this function used to build each row as
    # bare `{col_a: a_id, col_b: b_id, "origin": ...}`, discarding this locus
    # data even though `resolve_haplotype_private` computes it for both
    # directions. Confirmed against real production output: every one of 376
    # `both_haplotype_private` rows had ONLY the two internal ids + origin --
    # no chromosome, no position, nothing a user could locate. Never caught
    # earlier because every prior validation checked row COUNTS, not content.
    fwd_loc = (
        resolved_fwd.set_index([col_a, col_b])[["Chromosome", "Start_coordinate"]]
        if len(resolved_fwd) else None
    )
    rev_loc = (
        resolved_rev.set_index([col_a, col_b])[["Chromosome", "Start_coordinate"]]
        if len(resolved_rev) else None
    )

    rows = []
    absorbed = {a: set(), b: set()}
    for a_id, b_id in sorted(confirmed):
        a_chrom, a_pos = fwd_loc.loc[(a_id, b_id)]
        b_chrom, b_pos = rev_loc.loc[(a_id, b_id)]
        rows.append({
            col_a: a_id, col_b: b_id,
            f"Chromosome_{a}": a_chrom, f"Start_coordinate_{a}": a_pos,
            f"Chromosome_{b}": b_chrom, f"Start_coordinate_{b}": b_pos,
            "origin": "both_haplotype_private",
        })
        absorbed[a].add(a_id)
        absorbed[b].add(b_id)

    return pd.DataFrame(rows), absorbed


LOG_VERBOSE_FILENAME = "log_verbose.txt"
LOG_ERROR_FILENAME = "log_error.txt"

_SHARED_KEY_COLS = ["hg38_chr", "Strand_(fewest_mm+b)", "hg38_start"]

# `pandas.merge_asof` requires the "on" column to be sorted globally across
# the whole frame, not just within each `by` group (verified directly against
# the pinned dev-env's pandas 1.2.5). `_SHARED_KEY_COLS`' order -- (chrom,
# strand) primary, position last -- is correct for `cluster_collapse`'s
# chained clustering, but leaves `hg38_start` non-monotonic overall once more
# than one (chrom, strand) group is present (e.g. any real result set with
# both + and - strand hits), which raises "left/right keys must be sorted".
# Putting `hg38_start` first fixes this; `merge_asof`'s `by=` still correctly
# restricts matches to within the same (chrom, strand) group regardless of
# how rows from different groups are interleaved.
_MERGE_ASOF_SORT_COLS = ["hg38_start", "hg38_chr", "Strand_(fewest_mm+b)"]


def _combine_across_haplotypes(
    hg38_predictions: Dict[str, pd.DataFrame], names: List[str], merge_bp: int
) -> pd.DataFrame:
    """Combines two haplotypes' already within-haplotype-collapsed hg38
    predictions into one report, matching "same site" by proximity (within
    `merge_bp`) rather than exact-coordinate equality -- absorbs the 1-3bp
    bulge-registration drift documented in `cluster_collapse`'s docstring,
    which independently affects each haplotype's own lift-over anchor
    position.

    Deliberately does NOT reuse `cluster_collapse`'s transitive chaining:
    chaining is only correct within one haplotype's own results, where
    nearby rows are known artifacts of one real alignment event (bulge
    bookkeeping shifting the counted anchor position). Across two
    independent haplotypes, a bridging point could just as easily be an
    unrelated, genuinely distinct real site rather than a shared artifact --
    chaining across them risks merging real, distinct loci into one. Instead
    this does a direct, non-chaining nearest-within-tolerance match
    (`pandas.merge_asof`, one hop per row, no transitive linking).

    `merge_asof` matches each haplotype-`a` row independently, so the same
    haplotype-`b` row can end up the nearest match for more than one
    haplotype-`a` row (possible when two haplotype-`a` sites are each within
    `merge_bp` of the same haplotype-`b` site but more than `merge_bp` apart
    from each other -- otherwise they'd already have been collapsed within
    haplotype `a`'s own `cluster_collapse` pass). A single real haplotype-`b`
    site can't genuinely confirm two different haplotype-`a` sites as the
    same physical locus, so only the closest claim is kept; the other
    haplotype-`a` row is demoted back to `"<a>_only"`.

    Args:
        hg38_predictions: Each haplotype's already-`cluster_collapse`d
            (one row per real site, within that haplotype) predictions in
            hg38 coordinates.
        names: The two haplotype names, in the order defining which side's
            columns get which suffix.
        merge_bp: Same tolerance used for `cluster_collapse` -- should match
            the `--merge` value used for both haplotypes' `complete-search`
            runs.

    Returns:
        One row per reconciled site, with an `origin` column
        (`"<a>_only"`, `"<b>_only"`, or `"both"`).
    """
    a, b = names
    left = hg38_predictions[a].sort_values(_MERGE_ASOF_SORT_COLS).reset_index(drop=True)
    right = hg38_predictions[b].sort_values(_MERGE_ASOF_SORT_COLS).reset_index(drop=True).copy()

    # an entirely non-mappable haplotype leaves an empty frame here (e.g.
    # every prediction failed to lift over) -- pandas infers an `object`
    # dtype for an empty column, which merge_asof rejects when the other
    # side is a real int64 `hg38_start` ("incompatible merge keys" error).
    # Short-circuit rather than coerce dtypes into a merge that has nothing
    # to actually match against either way.
    if left.empty and right.empty:
        return pd.DataFrame(columns=list(left.columns) + ["origin"])
    if right.empty:
        out = left.rename(columns={
            col: f"{col}_{a}" for col in left.columns if col not in _SHARED_KEY_COLS
        })
        out["origin"] = f"{a}_only"
        return out
    if left.empty:
        out = right.rename(columns={
            col: f"{col}_{b}" for col in right.columns if col not in _SHARED_KEY_COLS
        })
        out["origin"] = f"{b}_only"
        return out

    # merge_asof collapses the "on" column to a single (left's) value, so
    # keep a copy of the right side's own position to compute the true gap
    # for tie-breaking below
    right["_right_hg38_start"] = right["hg38_start"]

    matched = pd.merge_asof(
        left, right, on="hg38_start", by=["hg38_chr", "Strand_(fewest_mm+b)"],
        tolerance=merge_bp, direction="nearest", suffixes=(f"_{a}", f"_{b}"),
    )
    b_id_col = f"off_target_id_{b}"
    matched["_gap"] = (matched["hg38_start"] - matched["_right_hg38_start"]).abs()

    is_matched = matched[b_id_col].notna()
    closest_claim_idx = matched[is_matched].groupby(b_id_col)["_gap"].idxmin()
    demote = is_matched & ~matched.index.isin(closest_claim_idx)
    b_cols = [c for c in matched.columns if c.endswith(f"_{b}")] + ["_right_hg38_start"]
    matched.loc[demote, b_cols] = pd.NA

    matched["origin"] = matched[b_id_col].notna().map({True: "both", False: f"{a}_only"})
    matched = matched.drop(columns=["_gap", "_right_hg38_start"])

    matched_b_ids = set(matched.loc[matched["origin"] == "both", b_id_col])
    b_only = right[~right["off_target_id"].isin(matched_b_ids)].drop(columns=["_right_hg38_start"]).copy()
    rename_map = {col: f"{col}_{b}" for col in b_only.columns if col not in _SHARED_KEY_COLS}
    b_only = b_only.rename(columns=rename_map)
    b_only["origin"] = f"{b}_only"

    return pd.concat([matched, b_only], ignore_index=True)


def reconcile_haplotypes(
    haplotypes: Dict[str, dict], workdir: str, merge_bp: int = 3,
    alignment_cache_root: Optional[str] = None, threads: int = 4,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Runs the full reconciliation pipeline for exactly two haplotypes.

    Writes its own progress to `<workdir>/log_verbose.txt` and, on failure,
    the error to `<workdir>/log_error.txt` -- the same filenames every other
    CRISPRme run (`complete-search`, `complete-test`) already writes, opened
    the same way `complete_search()` opens them: truncated fresh at the
    start of every call, not accumulated (`crisprme.py`'s own
    `open(..., "w")` on both files, right before launching the pipeline).
    `complete_search()` never actually re-enters an existing output folder
    (an earlier check blocks that), but `reconcile_haplotypes` deliberately
    does support being retried into the same `workdir` (see
    `haplotype_search_complete` in this module and its use in
    `crisprme.py`'s `assembly_search`), so each retry's log reflects only
    that attempt.

    Args:
        haplotypes: Exactly two entries, keyed by haplotype name (e.g.
            "paternal", "maternal"). Each value must have:
            `chrom_alias_file`, `chain_file`, `results_dir`, `genome_dir`
            (the haplotype's own genome folder, used for direct
            haplotype-vs-haplotype alignment). `results_prefix` is derived
            automatically via `find_results_prefix` if not given.
        workdir: Directory for intermediate BED files and the log files
            described above (created if missing).
        merge_bp: Locus-clustering threshold -- should match the `--merge`
            value used for both haplotypes' `complete-search` runs.
        alignment_cache_root: Directory for the cached, reusable direct
            haplotype-vs-haplotype alignment (see
            `build_or_reuse_haplotype_alignment`) -- a top-level cache,
            sibling of `genome_library/`, shared across runs against the
            same two genomes, not nested under `workdir`. `crisprme.py`'s
            `assembly_search()` always passes this (the direct-alignment
            step is on by default for real runs); left optional here, and
            skipped entirely if omitted, so this function stays usable in
            isolation (e.g. existing tests / callers that only care about
            the hg38-mapping path) without requiring haplotype genome
            folders or the `minimap2`/`impg` dependency.
        threads: Passed to the `minimap2`/`impg index` calls on a cache
            miss.

    Returns:
        `(combined, summary)` where `combined` is one row per reconciled
        locus with an `origin` column (`"<a>_only"`, `"<b>_only"`, `"both"`,
        or `"both_haplotype_private"` -- see `resolve_haplotype_private`),
        and `summary` has counts per category plus each haplotype's
        remaining non-mappable count (after direct-alignment resolution).

    Raises:
        ValueError: If `haplotypes` doesn't have exactly two entries.
    """
    if len(haplotypes) != 2:
        raise ValueError(
            f"reconcile_haplotypes requires exactly 2 haplotypes, got {len(haplotypes)} "
            "(N-way ploidy is not yet supported)"
        )
    names = list(haplotypes.keys())
    a, b = names[0], names[1]

    os.makedirs(workdir, exist_ok=True)
    log_verbose_path = os.path.join(workdir, LOG_VERBOSE_FILENAME)
    log_error_path = os.path.join(workdir, LOG_ERROR_FILENAME)
    open(log_verbose_path, "w").close()
    open(log_error_path, "w").close()

    def log(msg: str) -> None:
        with open(log_verbose_path, "a") as f:
            f.write(msg + "\n")

    try:
        predictions: Dict[str, pd.DataFrame] = {}
        lifted: Dict[str, pd.DataFrame] = {}
        unlifted_ids: Dict[str, set] = {}

        for name, cfg in haplotypes.items():
            log(f"Loading {name} predictions...")
            try:
                prefix = cfg.get("results_prefix") or find_results_prefix(cfg["results_dir"])
            except FileNotFoundError:
                # find_results_prefix raises FileNotFoundError for BOTH zero and
                # >1 *_integrated_results.tsv matches. Only the ZERO case is a
                # benign result to absorb: a haplotype whose search produced no
                # off-targets at all (a guide with no hits on that assembly) is
                # legitimate, not an error. The >1 case is an unsupported
                # multi-guide run and must still surface -- so re-raise it.
                n_matches = len(
                    glob.glob(os.path.join(cfg["results_dir"], "*_integrated_results.tsv"))
                )
                if n_matches != 0:
                    raise
                # Empty haplotype: substitute a well-formed empty prediction set
                # so reconciliation still runs -- every locus on the OTHER
                # haplotype then becomes <other>_only (never "both"), and this
                # haplotype contributes 0 non-mappable. Column schemas match the
                # real loaders (load_crisprme_predictions -> PRED_COLS +
                # off_target_id; load_lifted_bed -> hg38 cols + off_target_id)
                # so the merge and cluster_collapse below get a valid empty frame.
                log(f"{name}: no off-targets found -- treating as an empty result")
                predictions[name] = pd.DataFrame(columns=PRED_COLS + ["off_target_id"])
                lifted[name] = pd.DataFrame(
                    columns=["hg38_chr", "hg38_start", "hg38_end", "off_target_id"]
                )
                unlifted_ids[name] = set()
                continue
            _, ucsc_to_genbank = load_chrom_alias(cfg["chrom_alias_file"])

            predictions[name] = load_crisprme_predictions(cfg["results_dir"], prefix, merge_bp)
            check_chrom_alias_coverage(predictions[name], ucsc_to_genbank, name)
            log(f"{name}: {len(predictions[name])} unique predicted loci")

            bed_path = os.path.join(workdir, f"{name}_offtargets.bed")
            _, no_chrom_alias_ids = build_offtarget_bed(predictions[name], ucsc_to_genbank, bed_path)

            log(f"Lifting {name} predictions to hg38...")
            mapped_path = os.path.join(workdir, f"{name}_offtargets_lifted.bed")
            unmapped_path = os.path.join(workdir, f"{name}_offtargets_not_lifted.bed")
            run_liftover(bed_path, cfg["chain_file"], mapped_path, unmapped_path)

            lifted[name] = load_lifted_bed(mapped_path)
            liftover_rejected_ids = load_unlifted_ids(unmapped_path)
            check_liftover_failure_rate(
                len(predictions[name]) - len(no_chrom_alias_ids),
                len(liftover_rejected_ids),
                name,
            )
            # a missing chromAlias entry is just as unmappable-to-hg38 as a
            # genuine liftOver failure -- there's no chain-file name to look
            # it up under at all -- so it belongs in the same non-mappable
            # count, not silently dropped (see build_offtarget_bed's docstring)
            unlifted_ids[name] = liftover_rejected_ids | no_chrom_alias_ids
            log(f"{name}: {len(unlifted_ids[name])} non-mappable (no hg38 equivalent)")

        resolved_private = pd.DataFrame()
        if alignment_cache_root:
            log("Aligning haplotypes directly (bypassing hg38) to resolve haplotype-private predictions "
                "(both orientations, see resolve_haplotype_private_bidirectional)...")
            paf_target_a, paf_target_b = build_or_reuse_haplotype_alignments_both_orientations(
                haplotypes[a]["genome_dir"], haplotypes[b]["genome_dir"], a, b,
                alignment_cache_root, threads=threads,
            )
            resolved_private, absorbed = resolve_haplotype_private_bidirectional(
                unlifted_ids, predictions, paf_target_a, paf_target_b, [a, b], merge_bp,
            )
            log(f"Direct alignment resolved {len(resolved_private)} bidirectionally-confirmed haplotype-private pairs")
            unlifted_ids[a] -= absorbed[a]
            unlifted_ids[b] -= absorbed[b]

        hg38_predictions: Dict[str, pd.DataFrame] = {}
        for name in haplotypes:
            preds = predictions[name].copy()
            preds["off_target_id"] = preds["off_target_id"].astype(str)
            merged = preds.merge(lifted[name], on="off_target_id", how="inner")
            # Same fewest-mm+b criterion as load_crisprme_predictions()'s own
            # cluster_collapse() call, applied here post-liftover: two
            # native-coordinate loci can land on the same (or adjacent) hg38
            # destination after liftOver -- a real liftOver-collision case,
            # independent of the alt-file question load_crisprme_predictions()
            # handles, but with the identical CFD-vs-fewest-mm+b criterion bug
            # (both call sites fixed together -- this one is not a new
            # finding, just the same criterion bug applied post-liftover).
            mmb = (
                pd.to_numeric(merged["Mismatches_(fewest_mm+b)"], errors="coerce")
                + pd.to_numeric(merged["Bulges_(fewest_mm+b)"], errors="coerce")
            )
            merged = cluster_collapse(
                merged.assign(_mmb=mmb), "hg38_chr", "Strand_(fewest_mm+b)", "hg38_start",
                "_mmb", merge_bp, ascending=True,
            ).drop(columns=["_mmb"])
            hg38_predictions[name] = merged

        log("Combining lifted predictions across haplotypes...")
        combined = _combine_across_haplotypes(hg38_predictions, [a, b], merge_bp)
        if len(resolved_private):
            combined = pd.concat([combined, resolved_private], ignore_index=True)

        origin_counts = combined["origin"].value_counts()
        summary = {
            "both": int(origin_counts.get("both", 0)),
            f"{a}_only": int(origin_counts.get(f"{a}_only", 0)),
            f"{b}_only": int(origin_counts.get(f"{b}_only", 0)),
            "both_haplotype_private": int(origin_counts.get("both_haplotype_private", 0)),
            f"{a}_non_mappable": len(unlifted_ids[a]),
            f"{b}_non_mappable": len(unlifted_ids[b]),
        }
        log("Reconciliation complete:")
        for category, count in summary.items():
            log(f"  {category}: {count}")
    except Exception as e:
        with open(log_error_path, "a") as f:
            f.write(f"{e}\n")
        raise
    return combined, summary
