#!/usr/bin/env python3
"""Exact factorized worst-case CFD — 2.5.1 two-pass fast mode (Pass-2 scoring).

Design: docs/DESIGN_2.5.1_two_pass_fast_mode.md §7. `calc_cfd`
(new_simple_analysis.py) is a pure product of per-position factors × one PAM
factor, each factor in [0,1] and position-independent, so the MAX CFD over a
per-position allele SET factorizes EXACTLY:

    max_CFD = ( ∏_i  max_{a ∈ S_i} f(guide_i, a, i) )  ×  max_{p ∈ PamSet} pam(p)

computed in O(len guide) with NO haplotype enumeration. For a CONCRETE target it
is bit-identical to `new_simple_analysis.calc_cfd(..., do_scores=True)` (the panel
verified 0 ULP over 30K + ~945K cases). This is a fast-mode CONSTRUCT — validated
against `calc_cfd` (see test_twopass_cfd_exact.py), not yet an authoritative
replacement. A future C++ port must pass the same differential/oracle/shadow gates
(design §7) before being trusted. CRISPR-Bulge is NOT factorizable → out of scope here.

Caveats replicated from `calc_cfd` exactly: T→U on guide AND target before compare
and key; match → 1.0; mismatch → mm_scores["r{guide}:d{revcom(dna)},{pos}"] with a
missing key → 1.0 (the position-1 gap is free); an untranslatable IUPAC base →
0.0 (calc_cfd zeroes the whole score). Only meaningful in the canonical SpCas9
regime (20-nt guide, per-position factors); the caller gates that.

STDLIB + the shipped mismatch_score.pkl / PAM_scores.pkl.
"""
import os
import pickle

_HERE = os.path.dirname(os.path.realpath(__file__))
_BASECOMP = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "-": "-"}


def load_scores():
    """Load the shipped CFD pickles (mismatch_score.pkl, PAM_scores.pkl)."""
    with open(os.path.join(_HERE, "mismatch_score.pkl"), "rb") as fh:
        mm = pickle.load(fh)
    with open(os.path.join(_HERE, "PAM_scores.pkl"), "rb") as fh:
        pam = pickle.load(fh)
    return mm, pam


def _revcom(s):
    """EXACT `new_simple_analysis.revcom` semantics: reverse-complement over
    A/C/G/T/U/-, returning None if any base is untranslatable (IUPAC/N)."""
    try:
        return "".join(_BASECOMP[b] for b in reversed(s))
    except KeyError:
        return None


# Position-aware bulge scoring (2.7.11): the per-column factor + the guide-frame position
# rule live in bulge_cfd_penalty so the matrix scorer and the live scorer
# (cfd_score_bulge_aware) share ONE definition. Both are STDLIB; a pass-through fallback
# keeps the module importable if the helper is ever absent (mismatch-only stays exact).
try:
    from bulge_cfd_penalty import column_factor as _column_factor
except Exception:  # pragma: no cover - shipped alongside
    _column_factor = None


def position_factor(guide_u, dna_u, pos1, mm_scores):
    """Per-position CFD factor. For a match/mismatch this replicates calc_cfd's inner loop
    exactly; for a bulge column ('-') it applies the 2.7.11 position-aware factor via the
    shared primitive. guide_u / dna_u are single bases in T→U space; pos1 is the 1-based
    GUIDE-FRAME position (DNA-bulge columns deleted). Untranslatable base → 0.0."""
    if _column_factor is not None:
        return _column_factor(guide_u, dna_u, pos1, mm_scores)
    if guide_u == dna_u:
        return 1.0
    rc = _revcom(dna_u)
    if rc is None:
        return 0.0
    return mm_scores.get("r" + guide_u + ":d" + rc + "," + str(pos1), 1.0)


def cfd_concrete(guide, dna, pam, mm_scores, pam_scores):
    """Position-aware bulge-correct CFD for one concrete alignment (2.7.11); bit-identical
    to cfd_score_bulge_aware on the same alignment, and to the old calc_cfd for mismatch-only
    rows. ``guide``/``dna`` are the column-aligned spacer (may carry '-' bulges); ``pam`` is
    the 2-base PAM. DNA-bulge columns (guide '-') are deleted from the guide frame so the
    remaining columns keep their true Cas9 positions."""
    g = guide.replace("T", "U")
    d = dna.replace("T", "U")
    score = 1.0
    dna_gaps = 0  # guide '-' columns deleted so far -> shift subsequent guide-frame positions
    for i, (gb, db) in enumerate(zip(g, d)):
        score *= position_factor(gb, db, (i + 1) - dna_gaps, mm_scores)
        if gb == "-":
            dna_gaps += 1
    return score * pam_scores.get(pam, 0.0)


def cfd_worst_case(guide, allele_sets, pam_set, mm_scores, pam_scores):
    """EXACT max CFD over the per-position DNA allele SETS + the PAM set, via
    factorization — O(len guide), NO enumeration. Equals
    max over itertools.product(*allele_sets) × max over pam_set of cfd_concrete.

    allele_sets[i] is an iterable of DNA bases possible at protospacer position i
    (e.g. the variant-enriched allele set); pam_set is an iterable of possible
    PAM strings. Empty position set → 0.0 (no allele = no match).

    The factorization still holds with bulges: the DNA-bulge deletion pattern is fixed by
    the (allele-independent) ``guide``, so every column has a fixed guide-frame position and
    the worst case is the per-column max of the bulge-aware ``position_factor``."""
    g = guide.replace("T", "U")
    score = 1.0
    dna_gaps = 0
    for i, gb in enumerate(g):
        pos = (i + 1) - dna_gaps
        best = 0.0
        for a in allele_sets[i]:
            f = position_factor(gb, a.replace("T", "U"), pos, mm_scores)
            if f > best:
                best = f
        score *= best
        if gb == "-":
            dna_gaps += 1
    best_pam = max((pam_scores.get(p, 0.0) for p in pam_set), default=0.0)
    return score * best_pam
