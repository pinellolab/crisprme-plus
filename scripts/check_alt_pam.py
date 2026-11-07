#!/usr/bin/env python3
"""Release gate: no ALT-origin off-target may display a PAM the search did not match, and
no reported target may hold an N.

usage: check_alt_pam.py [--motif NRG] [--genome DIR|FASTA]
                        <integrated_results.tsv | results dir> ...

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

One exception (2.10.0): the REF text of an indel row (one whose Variant_info_genome names a
variant whose REF and ALT differ in length) is the reference at the row's coordinates, as
many bases as its target holds, and next to an assembly gap it can run into the gap while
the target, on the haplotype carrying the indel, holds none (on HPRC's chr22, the README
guide's rows at chr22:11978143 (-) and 11980225 (+)). Such an N passes where the reference
genome itself holds an N at that position, read from the row's start and strand (the - strand
reverse-complemented) and the text's gaps; any other N in it, and every N of an ALT text,
still fails. The genome is ``--genome`` (a folder of <chromosome>.fa files, as Genomes/hg38,
or one FASTA file), else the run's own: Genomes/<Genome_ref of the run's Params.txt> beside
its Results folder. Without one, such an N fails, as before.

Exit status 0 when every block passes both checks, 1 on any violation (each listed on
stdout), 2 on a usage error.
"""

import csv
import glob
import gzip
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


def _params(result_dir, key):
    for name in ("Params.txt", ".Params.txt"):
        path = os.path.join(result_dir, name)
        if os.path.isfile(path):
            with open(path) as fh:
                for line in fh:
                    k, _, value = line.rstrip("\n").partition("\t")
                    if k == key and value:
                        return value.strip()
    return None


class Genome(object):
    """The reference genome an indel row's REF text is read from, one chromosome at a time:
    a folder of <chromosome>.fa (.fasta, .fa.gz) files or one FASTA file."""

    def __init__(self, path):
        self.path = path
        self.name = None
        self.seq = None

    def _files(self, chrom):
        if os.path.isdir(self.path):
            for ext in (".fa", ".fasta", ".fa.gz", ".fasta.gz"):
                f = os.path.join(self.path, chrom + ext)
                if os.path.isfile(f):
                    return [f]
            return []
        return [self.path]

    def chromosome(self, chrom):
        """The chromosome's sequence (upper case), or None when the genome lacks it."""
        if chrom == self.name:
            return self.seq
        self.name, self.seq = chrom, None
        whole = os.path.isdir(self.path)  # a file per chromosome: all of it
        for f in self._files(chrom):
            op = gzip.open if f.endswith(".gz") else open
            parts, inside = [], whole
            with op(f, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if whole:
                            continue
                        if inside:
                            break  # the next chromosome of one FASTA file
                        inside = line[1:].split()[0] == chrom
                        continue
                    if inside:
                        parts.append(line.strip())
            if parts:
                self.seq = "".join(parts).upper()
        return self.seq


def run_genome(table):
    """The genome of the run a table belongs to: Genomes/<Genome_ref> beside its Results
    folder (<data>/Results/<job>/table), or None."""
    result_dir = os.path.dirname(os.path.abspath(table))
    name = _params(result_dir, "Genome_ref") or _params(result_dir, "Genome_selected")
    if not name:
        return None
    path = os.path.join(result_dir, "..", "..", "Genomes", name)
    return Genome(path) if os.path.isdir(path) else None


def _names_indel(variants):
    for tok in variants.split(","):
        p = tok.rsplit("_", 3)
        if len(p) == 4 and len(p[2]) != len(p[3]):
            return True
    return False


def n_where_the_reference_has_n(seq, start, strand, chrom, genome):
    """True iff every N (or n) of an indel row's REF text ``seq`` sits where the genome
    holds an N: the text's bases, gaps aside, are the reference from ``start`` on, in the
    target's orientation (reverse-complemented on the - strand)."""
    if genome is None:
        return False
    ref = genome.chromosome(chrom)
    if ref is None:
        return False
    bases = [c for c in seq if c != "-"]
    nb = len(bases)
    for j, c in enumerate(bases):
        if c.upper() != "N":
            continue
        q = start + j if strand == "+" else start + nb - 1 - j
        if not 0 <= q < len(ref) or ref[q] != "N":
            return False
    return True


def check_no_n(path, genome=None):
    """(targets checked, [violations], REF texts of indel rows that hold an N only where the
    reference does) for one integrated results table: every ``Aligned_protospacer+PAM_REF``
    and ``_ALT`` cell of each block that holds a target must hold no N or n, but the REF
    text of an indel row may where the reference genome itself holds an N
    (n_where_the_reference_has_n)."""
    csv.field_size_limit(1 << 30)
    checked = accepted = 0
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
                                  col["Start_coordinate_(%s)" % b], col["Strand_(%s)" % b],
                                  col.get("Variant_info_genome_(%s)" % b)))
                except KeyError:
                    continue
        if not cells:
            raise ValueError("%s: no Aligned_protospacer+PAM columns -- not an integrated "
                             "results table" % path)
        for row in reader:
            for b, side, i_seq, i_start, i_strand, i_var in cells:
                seq = row[i_seq]
                if seq in ("", "NA"):
                    continue
                checked += 1
                if "N" not in seq.upper():
                    continue
                chrom = row[col.get("Chromosome", 1)]
                if (side == "REF" and i_var is not None and _names_indel(row[i_var])
                        and n_where_the_reference_has_n(seq, int(row[i_start]), row[i_strand],
                                                        chrom, genome)):
                    accepted += 1
                    continue
                bad.append((b, side, chrom, row[i_start], row[i_strand], seq))
    return checked, bad, accepted


def main(argv):
    motif = None
    genome_path = None
    paths = []
    it = iter(argv[1:])
    for arg in it:
        if arg == "--motif":
            motif = next(it, None)
        elif arg == "--genome":
            genome_path = next(it, None)
            if not genome_path or not os.path.exists(genome_path):
                print("check_alt_pam: --genome %s: no such file or folder" % genome_path)
                return 2
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
    given = Genome(genome_path) if genome_path else None
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
        genome = given or run_genome(t)
        targets, with_n, accepted = check_no_n(t, genome)
        print("%s  %s: %d reported targets (REF and ALT, every block), %d with an N%s"
              % ("FAIL" if with_n else "PASS", t, targets, len(with_n),
                 "; %d indel rows' REF texts hold an N only where the reference genome does"
                 % accepted if accepted else ""))
        for b in with_n[:50]:
            print("    %s %s %s:%s(%s) %s" % b)
        if len(with_n) > 50:
            print("    ... %d more" % (len(with_n) - 50))
        if with_n and genome is None and any(side == "REF" for _b, side, *_r in with_n):
            print("    (no reference genome to read an indel row's REF text against: pass "
                  "--genome, or run this on a results folder inside its data folder)")
        failed += bool(bad) or bool(with_n)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
