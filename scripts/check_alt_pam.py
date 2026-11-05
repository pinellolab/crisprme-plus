#!/usr/bin/env python3
"""Release gate: no ALT-origin off-target may display a PAM the search did not match, and
no reported target may hold an N.

usage: check_alt_pam.py [--motif NRG] <integrated_results.tsv | results dir> ...

ALT PAM. For every row of an integrated results table, and for each of its three blocks
(highest CFD, fewest mismatches+bulges, highest CRISPR-Bulge), an ALT-origin block's
``Aligned_protospacer+PAM_ALT`` must hold a PAM that satisfies the searched motif:
every PAM base a concrete A/C/G/T that the motif allows at its position. The PAM is
the trailing motif-wide window when the aligned spacer ends with N's (a 3' PAM, e.g.
SpCas9 NRG), the leading one when it starts with them (a 5' PAM, e.g. Cas12a TTTV).
The motif comes from ``--motif`` or from the ``Pam`` line of the run's
``Params.txt`` / ``.Params.txt``.

Up to 2.8.0 the indel post-analysis read an IUPAC code at a PAM column of a SNP-overlaid
fake-indel target as the code's first allele, so an R (A/G) at NRG's last base was
shown as A: a genome-wide HGDP + 1000 Genomes search reported 1,088 such rows in its
highest-CFD block. A run of a release must report none.

BLIND SPOT: the check sees only a base the motif forbids. Under a motif that allows the
base the misreading produced, e.g. NRR (A allowed at the last position) or an N
position, the same misreading passes unseen; the per-haplotype unit tests
(PostProcess/test_indel_haplotype_resolution.py) and the synthetic oracle
(test/indel_snp_oracle) cover those.

NO N. In every block of every row, REF- and ALT-origin alike, neither
``Aligned_protospacer+PAM_REF`` nor ``Aligned_protospacer+PAM_ALT`` may hold an N (or an n,
an N the guide does not match), PAM included. An N of the genome is an assembly gap (4.88%
of hg38, 150,630,719 bases in 1,001 runs), which matches no guide and no PAM: the engine
drops a target with an N in its guide region and matches an N to no PAM base. A reported
target with an N means that guard gave way, and the genome's gaps could come out as
off-targets. An empty or NA cell (the ALT of a REF row) holds no target.

Exit status 0 when every block passes both checks, 1 on any violation (each listed on
stdout), 2 on a usage error.
"""

import csv
import glob
import os
import sys

_BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
          "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
          "V": "ACG", "N": "ACGT"}
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _motif_from_params(result_dir):
    for name in ("Params.txt", ".Params.txt"):
        path = os.path.join(result_dir, name)
        if os.path.isfile(path):
            with open(path) as fh:
                for line in fh:
                    key, _, value = line.rstrip("\n").partition("\t")
                    if key == "Pam" and value:
                        return value.strip().upper()
    return None


def pam_ok(pam, motif):
    """True iff every base of ``pam`` is a concrete base the motif allows there."""
    return len(pam) == len(motif) and all(
        b in "ACGT" and b in _BASES.get(m, "") for b, m in zip(pam.upper(), motif.upper()))


def check(path, motif):
    """(rows, ALT blocks checked, [violations]) for one integrated results table."""
    csv.field_size_limit(1 << 30)
    n = len(motif)
    rows = checked = 0
    bad = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        head = next(reader)
        col = {name: i for i, name in enumerate(head)}
        blocks = []
        for b in BLOCKS:
            try:
                blocks.append((b, col["REF/ALT_origin_(%s)" % b],
                               col["Aligned_spacer+PAM_(%s)" % b],
                               col["Aligned_protospacer+PAM_ALT_(%s)" % b],
                               col["Start_coordinate_(%s)" % b], col["Strand_(%s)" % b]))
            except KeyError:
                continue
        if not blocks:
            raise ValueError("%s: no REF/ALT origin columns -- not an integrated "
                             "results table" % path)
        for row in reader:
            rows += 1
            for b, i_or, i_sp, i_alt, i_start, i_strand in blocks:
                if row[i_or].strip().lower() != "alt" or row[i_alt] in ("", "NA"):
                    continue
                spacer, alt = row[i_sp], row[i_alt]
                five = spacer[:n].upper() == "N" * n and spacer[-n:].upper() != "N" * n
                pam = alt[:n] if five else alt[-n:]
                checked += 1
                if not pam_ok(pam, motif):
                    bad.append((b, row[col.get("Chromosome", 1)], row[i_start], row[i_strand],
                                alt, pam))
    return rows, checked, bad


def check_no_n(path):
    """(targets checked, [violations]) for one integrated results table: every
    ``Aligned_protospacer+PAM_REF`` and ``_ALT`` cell of each block that holds a target
    must hold no N or n."""
    csv.field_size_limit(1 << 30)
    checked = 0
    bad = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        head = next(reader)
        col = {name: i for i, name in enumerate(head)}
        cells = []
        for b in BLOCKS:
            for side in ("REF", "ALT"):
                try:
                    cells.append((b, side, col["Aligned_protospacer+PAM_%s_(%s)" % (side, b)],
                                  col["Start_coordinate_(%s)" % b], col["Strand_(%s)" % b]))
                except KeyError:
                    continue
        if not cells:
            raise ValueError("%s: no Aligned_protospacer+PAM columns -- not an integrated "
                             "results table" % path)
        for row in reader:
            for b, side, i_seq, i_start, i_strand in cells:
                seq = row[i_seq]
                if seq in ("", "NA"):
                    continue
                checked += 1
                if "N" in seq.upper():
                    bad.append((b, side, row[col.get("Chromosome", 1)], row[i_start],
                                row[i_strand], seq))
    return checked, bad


def main(argv):
    motif = None
    paths = []
    it = iter(argv[1:])
    for arg in it:
        if arg == "--motif":
            motif = next(it, None)
        elif arg in ("-h", "--help"):
            print(__doc__)
            return 0
        else:
            paths.append(arg)
    if not paths:
        sys.stdout.write(__doc__)
        return 2
    tables = []
    for p in paths:
        if os.path.isdir(p):
            found = sorted(glob.glob(os.path.join(p, "*integrated_results.tsv")))
            if not found:
                print("check_alt_pam: no *integrated_results.tsv in %s" % p)
                return 2
            tables.extend(found)
        else:
            tables.append(p)
    failed = 0
    for t in tables:
        m = motif or _motif_from_params(os.path.dirname(os.path.abspath(t)))
        if not m:
            print("check_alt_pam: no --motif and no Pam line in %s's Params.txt" % t)
            return 2
        rows, checked, bad = check(t, m)
        status = "FAIL" if bad else "PASS"
        print("%s  %s: %d rows, %d ALT blocks, %d with a PAM that violates %s"
              % (status, t, rows, checked, len(bad), m))
        for b in bad[:50]:
            print("    %s %s:%s(%s) %s PAM %s" % b)
        if len(bad) > 50:
            print("    ... %d more" % (len(bad) - 50))
        targets, with_n = check_no_n(t)
        print("%s  %s: %d reported targets (REF and ALT, every block), %d with an N"
              % ("FAIL" if with_n else "PASS", t, targets, len(with_n)))
        for b in with_n[:50]:
            print("    %s %s %s:%s(%s) %s" % b)
        if len(with_n) > 50:
            print("    ... %d more" % (len(with_n) - 50))
        failed += bool(bad) or bool(with_n)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
