#!/usr/bin/env python
"""The notes a chromosome's post-analysis ends with (stdout, log_verbose.txt).

One process writes them at its end. When a chromosome's post-analysis runs in pieces
(2.10.0, option P: postproc_pool.py), a piece writes none of them: it leaves its counts in
<piece>.piece.json, and merge_pieces.py writes the chromosome's notes from their sums and
from the merged files, with these same functions, so the lines are those of one process.

Stdout only: the search treats any stderr write as a failure. Stdlib only, apart from
window_density, so merge_pieces.py imports it without the analysis scripts' module globals.
"""

import window_density as _wd


def snp_counts():
    """The counts of a SNP post-analysis its notes are written from: the pre-check and the
    K-skip, --max-degeneracy, and the SNP rows' carriers an indel takes off."""
    return {
        "precheck": {"seen": 0, "skipped": 0, "kseen": 0, "kskipped": 0},
        "degeneracy": {"cut": 0, "withheld": 0},
        "snp_indel": {"targets": 0, "listings": 0, "gone": 0, "putative": 0, "kept": 0},
    }


def indel_counts():
    """The counts of an indel post-analysis its notes are written from: how its overlaid
    targets were resolved, and its rows of another contig."""
    return {
        "overlaid": {"targets": 0, "rows": 0, "dropped": 0, "putative": 0, "unmapped": 0,
                     "second_indel": 0, "flank": 0, "prechecked": 0},
        "foreign_rows": 0,
    }


def add_counts(total, more):
    """Add ``more`` (nested dicts of integers) into ``total``; returns ``total``."""
    for key, value in more.items():
        if isinstance(value, dict):
            add_counts(total.setdefault(key, {}), value)
        else:
            total[key] = total.get(key, 0) + value
    return total


def dense_window_counts(letters, iupac_cap):
    """(over the cap, hyperdense): how many of the windows whose letters are ``letters``
    spell more than 2^``iupac_cap`` sequences, and 2^HYPERDENSE_BITS or more."""
    n_over = n_hyper = 0
    for window in letters:
        n_seq = _wd.sequences(window)
        n_over += _wd.over_cap(n_seq, iupac_cap)
        n_hyper += _wd.hyperdense(n_seq)
    return n_over, n_hyper


def density_lines(chrom, out_name, n_four_allele, n_over, n_hyper, iupac_cap):
    """The notes on the four-allele positions and the windows of the high-variant-density
    BED (``out_name``: the output's file name, without its folder)."""
    lines = []
    if n_four_allele:
        lines.append(f"Four-allele positions (Q: every base is an observed allele there) in "
                     f"the windows of {chrom}: {n_four_allele}")
    if n_over or n_hyper:
        over = (f"{n_over} spell more than 2^{iupac_cap} sequences "
                f"(dense: over CRISPRME_IUPAC_CAP), " if iupac_cap >= 0 else "")
        lines.append(f"High-variant-density windows of {chrom}: {over}"
                     f"{n_hyper} spell 2^{_wd.HYPERDENSE_BITS} or more "
                     f"(hyperdense); listed in {out_name}.high_variant_density_regions.bed")
    return lines


def snp_lines(chrom, counts, precheck_on, fast_mode, max_degeneracy):
    """The notes on the pre-check or the K-skip, --max-degeneracy, and the SNP rows'
    carriers an indel takes off, from ``counts`` (snp_counts)."""
    pre, deg, si = counts["precheck"], counts["degeneracy"], counts["snp_indel"]
    lines = []
    if pre["seen"]:
        lines.append(f"Haplotype pre-check on {chrom}: {pre['skipped']} of {pre['seen']} raw "
                     f"target(s) with a variant letter skipped (no sample's genotypes can make "
                     f"them a row within the budget and the PAM; their reference rows are "
                     f"written)")
    elif not precheck_on and (not fast_mode or max_degeneracy >= 0):
        lines.append(f"Haplotype pre-check on {chrom}: off (CRISPRME_HAPLOTYPE_PRECHECK=0)")
    if pre["kseen"]:
        lines.append(f"K-skip on {chrom}: {pre['kskipped']} of {pre['kseen']} raw target(s) "
                     f"with a variant letter skipped (every row within the budget and the PAM "
                     f"would rest on more than --max-degeneracy {max_degeneracy} ambiguous "
                     f"column(s))")
    if max_degeneracy >= 0:
        lines.append(f"--max-degeneracy {max_degeneracy} on {chrom}: {deg['cut']} worst-case "
                     f"alignment(s) cut back to rest on at most {max_degeneracy} ambiguous "
                     f"column(s); {deg['withheld']} row(s) resting on more, which no phased "
                     f"carrier confirms, withheld")
    if si["targets"] or si["putative"] or si["kept"]:
        lines.append("snp-indel: %d SNP target(s) lost %d carrier listing(s) whose haplotype "
                     "holding the SNPs also holds an indel that changes the target's bases; %d "
                     "haplotype(s) left with no carrier are not reported; %d carrier "
                     "haplotype(s) kept, PUTATIVE, beside an unphased indel; %d kept, the "
                     "indel leaving the target's bases as they are"
                     % (si["targets"], si["listings"], si["gone"], si["putative"], si["kept"]))
    return lines


def indel_lines(chrom, counts):
    """The notes of an indel post-analysis, from ``counts`` (indel_counts)."""
    ovl = counts["overlaid"]
    lines = []
    if counts["foreign_rows"]:
        lines.append("WARNING: %d indel target row(s) in the %s subset carried another "
                     "contig's label and were dropped; they would otherwise have been "
                     "reported against %s indels they do not belong to."
                     % (counts["foreign_rows"], chrom, chrom))
    if ovl["flank"]:
        lines.append("indel-snp: %d target(s) of the fake-indel windows lie in a flank and "
                     "hold no base of the indel: not indel off-targets, dropped" % ovl["flank"])
    if ovl["second_indel"]:
        lines.append("indel-snp: %d target(s) without a SNP code lost the carriers whose "
                     "haplotype holds a second indel that changes their bases"
                     % ovl["second_indel"])
    if not ovl["targets"]:
        return lines
    lines.append("indel-snp: %s overlaid target(s) resolved per carrier haplotype into %d "
                 "row(s) (%d PUTATIVE); %d had no carrier with a target there"
                 % (ovl["targets"], ovl["rows"], ovl["putative"], ovl["dropped"]))
    if ovl["prechecked"]:
        lines.append("indel-snp: %d of them told by the carriers-first pre-check, before the "
                     "second indels of the phased store were read" % ovl["prechecked"])
    if ovl["unmapped"]:
        lines.append("WARNING [indel-snp]: %d IUPAC column(s) of %s did not map onto the "
                     "reference base they code; they were resolved without genotypes"
                     % (ovl["unmapped"], chrom))
    return lines


def write(lines):
    """Print ``lines`` on stdout, flushed (each its own line)."""
    for line in lines:
        print(line, flush=True)
