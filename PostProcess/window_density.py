#!/usr/bin/env python3
"""The ambiguity budget of a variant search: how permissive a window of the variant genome
is, and how much of a reported alignment rests on its ambiguous columns.

A variant genome writes each SNP position as the IUPAC letter of the bases observed there (R
for A or G, B for C, G or T, ...), and a four-allele position, where the reference base and
three ALTs cover A, C, G and T, as Q. N stays an assembly gap: it matches nothing, so it adds
no sequence to a window.

F, a window's bits: log2 of the number of sequences its letters spell, the product of the
bases each letter stands for (1 for A, C, G, T and N, 2 for R, Y, S, W, K and M, 3 for B, D,
H and V, 4 for Q). A window of more than 2^CRISPRME_IUPAC_CAP sequences is dense, as the
report and the docs have called it since 2.8.1: a per-sample search of an index without the
dict-less tiers reports it by the haplotypes its carriers have, read from the dictionaries as
the tiers' path reads them, instead of every combination of its alleles (a population-level
search reports every window by its greedy worst-case alignment; a per-sample one on the tiers
any window by its carriers' own haplotypes). One of 2^HYPERDENSE_BITS sequences or more is
hyperdense: it matches nearly any guide somewhere, and is flagged whatever
CRISPRME_IUPAC_CAP is.

D, the ambiguous columns an alignment rests on: those where the reference base would not fit
(a mismatch against the guide, or a base the PAM does not admit) and the alignment's base
does. complete-search's --max-degeneracy caps it.

Stdlib only and free of new_simple_analysis's module globals (which open sys.argv[1] at
import), so it is unit-tested alone, as twopass_emit is.
"""

import math

# A window spelling 2^HYPERDENSE_BITS sequences or more is flagged hyperdense. The one place
# to tune it: new_simple_analysis.py and generate_report.py read it here.
HYPERDENSE_BITS = 20

# The sites-only density guard (CRISPRme 2.10.0): a candidate of a sites-only dataset (the
# mega) whose window spells 2^SITES_ONLY_GUARD_BITS sequences or more is withheld from the
# SNP and the indel post-analysis, its window listed instead (sites_only_guard.py). A
# sites-only dataset has no genotypes to tell which of a window's combinations anyone
# carries, and every allele coded makes the immunoglobulin loci spell 2^40 and more: such a
# window matches nearly any guide, and its rows are combinations of 10 to 17 variants.
# complete-search --sites-only-density-guard sets it (CRISPRME_SITES_ONLY_GUARD), and
# --no-sites-only-density-guard turns the guard off.
SITES_ONLY_GUARD_BITS = HYPERDENSE_BITS

# CRISPRME_IUPAC_CAP's default, in bits: above 2^10 sequences the dictionary path reports a
# window by the haplotypes its carriers have, not every combination of its alleles. That is ten
# two-allele columns, which the cap allowed when it counted ambiguous positions, up to 2.8.1.
DEFAULT_IUPAC_CAP = 10

# The letter of a four-allele position in a variant genome (CRISPRitz 2.10.0 and later).
FOUR_ALLELE = "Q"

# The bases each letter of a variant genome stands for (an unknown letter counts one).
_WIDTH = {}
for _letters, _width in (("ACGTN", 1), ("RYSWKM", 2), ("BDHV", 3), (FOUR_ALLELE, 4)):
    for _c in _letters:
        _WIDTH[_c] = _WIDTH[_c.lower()] = _width

# The bases each letter of a PAM admits (a PAM's N is any base).
_PAM_ADMITS = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def sequences(window):
    """How many sequences ``window`` spells: the product, over its letters, of the bases each
    stands for. A bulge gap ('-') is no base of the genome and counts one."""
    n = 1
    for c in window:
        n *= _WIDTH.get(c, 1)
    return n


# the letters of more than one base, by how many they stand for (guard_withholds)
_KEEP_AMBIGUOUS = {i: None for i in range(256)}
for _c in "RYSWKMBDHV" + FOUR_ALLELE:
    _KEEP_AMBIGUOUS.pop(ord(_c), None)
    _KEEP_AMBIGUOUS.pop(ord(_c.lower()), None)


def guard_withholds(window, guard_bits):
    """True when ``window`` spells 2^guard_bits sequences or more, as sequences() counts
    them (a bulge gap and a letter of one base count one): the sites-only density guard
    withholds such a candidate. Fast, as it runs over every raw target of a sites-only
    dataset: the letters standing for more than one base are kept by one translate, and
    fewer than guard_bits / 2 of them cannot reach the bound. A negative bound is none."""
    if guard_bits < 0:
        return False
    amb = window.translate(_KEEP_AMBIGUOUS)
    if 2 * len(amb) < guard_bits:
        return False
    n4 = amb.count(FOUR_ALLELE) + amb.count(FOUR_ALLELE.lower())
    n3 = sum(amb.count(c) for c in "BDHVbdhv")
    n2 = len(amb) - n3 - n4
    return (1 << (n2 + 2 * n4)) * 3 ** n3 >= 1 << guard_bits


def bits(n_sequences):
    """F: log2 of a window's number of sequences."""
    return math.log2(n_sequences)


def over_cap(n_sequences, cap_bits):
    """True when a window spells more than 2^cap_bits sequences; a negative cap is none."""
    return cap_bits >= 0 and n_sequences > (1 << cap_bits)


def hyperdense(n_sequences):
    """True when a window spells 2^HYPERDENSE_BITS sequences or more."""
    return n_sequences >= (1 << HYPERDENSE_BITS)


def variant_columns(window):
    """How many columns of ``window`` are variant positions: letters standing for two bases
    or more, a four-allele position (Q) among them."""
    return sum(1 for c in window if _WIDTH.get(c, 1) > 1)


def four_allele_offsets(window):
    """The offsets of the four-allele positions (Q) of ``window``."""
    return [i for i, c in enumerate(window) if c in "Qq"]


def degeneracy(row, ref, guide, pam, pos_beg, pos_end, pam_begin, pam_end):
    """D of one alignment: its columns where the reference base would not fit and the row's
    base does. ``row`` and ``ref`` are the aligned target and its reference, bulge gaps in
    place, in either case (new_simple_analysis lower-cases a mismatch); ``guide`` is the
    aligned guide without its PAM (a gap at a DNA bulge); ``pam`` the PAM's letters, and the
    four slice bounds new_simple_analysis's (the protospacer is ``[pos_beg:pos_end]``, the
    PAM ``[pam_begin:pam_end]``). A column with a gap in either sequence or in the guide is
    a bulge, which no base fits."""
    d = 0
    for i, (r, a) in enumerate(zip(ref[pos_beg:pos_end], row[pos_beg:pos_end])):
        g = guide[i] if i < len(guide) else "-"
        if "-" in (g, r, a):
            continue
        if r.upper() != g.upper() and a.upper() == g.upper():
            d += 1
    for j, (r, a) in enumerate(zip(ref[pam_begin:pam_end], row[pam_begin:pam_end])):
        admits = _PAM_ADMITS.get(pam[j].upper(), "") if j < len(pam) else ""
        if r.upper() not in admits and a.upper() in admits and a != "-":
            d += 1
    return d
