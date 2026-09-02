#!/usr/bin/env python3
"""Pass-1 of the 2.5.1 two-pass fast-mode search: score-free region find.

Design: docs/DESIGN_2.5.1_two_pass_fast_mode.md sections 2-3, 5-6. For every
candidate window this computes the COLLAPSED-IUPAC MIN EDIT `D` -- the minimum
mismatches+bulges of the guide against a window where each position admits an
ALLELE SET (the variant-enriched IUPAC), with NO haplotype enumeration:

    D = min over alignments of  ( #positions where guide_i not-in S_j ) + ( #bulges )

`D` is a LOWER BOUND on the true edit distance of any real haplotype carried at
that window (proof, doc section 5: per-position set-containment for SNPs +
alignment-reuse for bulges -- adversarially verified, 0 failures over ~1.3M
cases). So thresholding `D <= k` is LOSSLESS FOR DETECTION, for the cases the
bound covers. It BREAKS only for >=2 cis indels in one protospacer and
equal-length MNVs (doc section 6); those windows are flagged, not silently
dropped. This module is the enumeration-free replacement for the intractable 2^k
haplotype-lattice post-analysis (49h+ on the 4x dense panel).

STDLIB only. Pass-2 scoring (exact CFD via twopass_cfd_exact, min-edit, CRISTA
shortlist) runs only on the windows this pass keeps.
"""
import functools

# IUPAC ambiguity codes -> the concrete base set they denote (the enriched
# reference / fake-indel genome stores overlaid SNPs as IUPAC; a plain ACGT base
# is the singleton set). '-'/'.' are NOT window content (gaps come from bulges).
IUPAC_TO_SET = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "U": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("GC"),
    "W": frozenset("AT"), "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"), "H": frozenset("ACT"),
    "V": frozenset("ACG"), "N": frozenset("ACGT"),
}

_INF = (1 << 30, 0, 0, 0)  # (D, mm, dna_bulges, rna_bulges) sentinel = unreachable


def iupac_window_to_sets(window):
    """Decode an IUPAC window string (enriched genome region) into per-position
    allele sets. An unknown character maps to the empty set (nothing matches ->
    forces a mismatch there), never a crash."""
    return [IUPAC_TO_SET.get(b.upper(), frozenset()) for b in window]


def collapsed_min_edit(guide, allele_sets, max_bdna=2, max_brna=2):
    """The collapsed-IUPAC min edit `D` of ``guide`` GLOBALLY aligned to a window
    whose position j admits any base in ``allele_sets[j]``.

    A position "matches" (0 mismatch) iff ``guide[i] in allele_sets[j]`` -- the
    set-membership relaxation that makes this a lower bound without enumerating
    which concrete allele is present. Bulges: a DNA bulge consumes a WINDOW base
    with no guide base (target insertion); an RNA bulge consumes a GUIDE base with
    no window base (guide insertion); each costs 1. Both ends of the guide AND the
    window are fully consumed (the window is the exact candidate protospacer, so
    its length is G - rna_bulges .. G + dna_bulges).

    Returns ``(D, mm, dna_bulges, rna_bulges)`` for the minimal-`D` alignment
    (ties broken toward fewer mismatches), or ``None`` if no alignment fits the
    bulge caps. `D == mm + dna_bulges + rna_bulges`.
    """
    G, W = len(guide), len(allele_sets)
    g = guide.upper()

    @functools.lru_cache(maxsize=None)
    def rec(i, j, dbl, rbl):
        if i == G and j == W:
            return (0, 0, 0, 0)
        best = _INF
        if i < G and j < W:  # consume a guide base against a window base
            sub = rec(i + 1, j + 1, dbl, rbl)
            mmc = 0 if g[i] in allele_sets[j] else 1
            cand = (sub[0] + mmc, sub[1] + mmc, sub[2], sub[3])
            if cand < best:
                best = cand
        if j < W and dbl > 0:  # DNA bulge: extra window base, no guide base
            sub = rec(i, j + 1, dbl - 1, rbl)
            cand = (sub[0] + 1, sub[1], sub[2] + 1, sub[3])
            if cand < best:
                best = cand
        if i < G and rbl > 0:  # RNA bulge: extra guide base, no window base
            sub = rec(i + 1, j, dbl, rbl - 1)
            cand = (sub[0] + 1, sub[1], sub[2], sub[3] + 1)
            if cand < best:
                best = cand
        return best

    res = rec(0, 0, max_bdna, max_brna)
    rec.cache_clear()
    return None if res[0] >= _INF[0] else res


def window_min_edit(guide, iupac_window, max_bdna=2, max_brna=2):
    """collapsed_min_edit on an IUPAC window STRING (decodes to allele sets first).
    Returns the same ``(D, mm, dna_bulges, rna_bulges)`` tuple, or None."""
    return collapsed_min_edit(guide, iupac_window_to_sets(iupac_window),
                              max_bdna=max_bdna, max_brna=max_brna)


def min_edit_over_windows(guide, windows, max_bdna=2, max_brna=2):
    """`D = min` over a set of candidate windows -- the enriched-reference window
    UNION each single-indel fake-contig window (doc section 3). Each ``windows``
    item is an IUPAC string. Returns ``(D, mm, dna_bulges, rna_bulges, which)``
    where ``which`` is the index of the winning window, or None if none align."""
    best, best_i = None, -1
    for idx, w in enumerate(windows):
        r = window_min_edit(guide, w, max_bdna=max_bdna, max_brna=max_brna)
        if r is not None and (best is None or r[0] < best[0]):
            best, best_i = r, idx
    return None if best is None else (best[0], best[1], best[2], best[3], best_i)


def detects(guide, iupac_window, max_mm, max_bdna, max_brna):
    """Pass-1 detection predicate: True iff the window's collapsed min edit is
    within the search budget (`D <= max_mm + max_bdna + max_brna` and the bulge
    components fit). Because `D` lower-bounds every real haplotype's edit, a False
    here PROVABLY has no in-budget carried off-target at this window (for the cases
    the bound covers -- multi-indel/MNV windows must be flagged separately)."""
    r = window_min_edit(guide, iupac_window, max_bdna=max_bdna, max_brna=max_brna)
    if r is None:
        return False
    D, mm, db, rb = r
    return mm <= max_mm and db <= max_bdna and rb <= max_brna
