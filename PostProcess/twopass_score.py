#!/usr/bin/env python3
"""Pass-2 of the 2.5.1 two-pass fast mode: per-score worst-POSSIBLE annotation.

Design: docs/DESIGN_2.5.1_two_pass_fast_mode.md sections 2, 4, 7. Runs ONLY on the
windows Pass-1 (twopass_find) kept, and computes, for each requested score, the
WORST-POSSIBLE value over the window's allele sets -- enumeration-free where the
score factorizes:

  * min-edit : `D` from Pass-1 directly (already the worst-case = smallest edit any
               concrete haplotype can achieve; a lower bound on the true edit).
  * CFD      : EXACT worst-case via factorization (twopass_cfd_exact.cfd_worst_case),
               O(guide) per window, bit-identical to calc_cfd, NO enumeration.
  * CRISTA   : a RandomForest -> NOT factorizable. Best-effort: evaluate on a bounded
               candidate SET of low-edit concrete haplotypes and FLAG it approximate
               when the set is truncated. Never claimed exhaustive (see the report).

"Worst-POSSIBLE" (not worst-observed) is deliberate: unphased / aggregate datasets
have no per-sample haplotypes, so this bounds risk over ALL formable alleles; the
phased tiers are only an ANNOTATION ("is this worst case actually carried?").

STDLIB + twopass_cfd_exact + twopass_find (all in-repo).
"""
import itertools
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import twopass_cfd_exact as _cfd  # noqa: E402
import twopass_find as _find  # noqa: E402

# Bound on how many concrete haplotypes CRISTA is evaluated over before the result
# is flagged approximate (the factorization does not apply to a RandomForest).
CRISTA_CANDIDATE_CAP = 256


def _iter_low_edit_haplotypes(guide, allele_sets, cap):
    """Yield up to ``cap`` concrete haplotypes from the window's allele sets,
    ordered so the LOWEST-edit (closest-to-guide) ones come first: at each position
    prefer the guide's own base when it is in the set. Deterministic. Returns
    (haplotypes, truncated) where ``truncated`` is True if the full product exceeded
    ``cap`` (=> CRISTA is approximate over the shortlist, not exhaustive)."""
    ordered = []
    total = 1
    for i, s in enumerate(allele_sets):
        bases = sorted(s)
        gb = guide[i] if i < len(guide) else None
        if gb in bases:  # put the guide-matching base first (drives the low-edit end)
            bases = [gb] + [b for b in bases if b != gb]
        ordered.append(bases)
        total *= max(1, len(bases))
    truncated = total > cap
    out = []
    for combo in itertools.product(*ordered):
        out.append("".join(combo))
        if len(out) >= cap:
            break
    return out, truncated


def score_window(guide, allele_sets, pam_set, mm_scores, pam_scores,
                 crista_fn=None, crista_cap=CRISTA_CANDIDATE_CAP,
                 max_bdna=2, max_brna=2):
    """Worst-possible scores for ONE Pass-1 candidate window.

    Args:
      guide: the protospacer guide (ACGT, length G).
      allele_sets: per-position iterable of admissible DNA bases (Pass-1 decoded the
        enriched-IUPAC window; length ~G).
      pam_set: iterable of possible PAM strings at the window.
      mm_scores/pam_scores: the shipped CFD pickles (twopass_cfd_exact.load_scores()).
      crista_fn: optional callable(list_of_(guide,dna,pam)) -> list_of_scores. When
        None, CRISTA is skipped (min-edit + CFD still returned).

    Returns a dict:
      min_edit         : (D, mm, dna_bulges, rna_bulges) from Pass-1, or None.
      cfd_worst        : exact worst-case CFD over the allele sets x pam_set.
      crista_worst     : max CRISTA over the evaluated shortlist (None if no crista_fn).
      crista_approx    : True if the shortlist was truncated (CRISTA not exhaustive).
      crista_n_evaluated: shortlist size actually scored.
    """
    sets = [frozenset(s) for s in allele_sets]
    out = {
        "min_edit": _find.collapsed_min_edit(guide, sets, max_bdna, max_brna),
        "cfd_worst": _cfd.cfd_worst_case(guide, sets, list(pam_set),
                                         mm_scores, pam_scores),
        "crista_worst": None,
        "crista_approx": False,
        "crista_n_evaluated": 0,
    }
    if crista_fn is not None:
        haps, truncated = _iter_low_edit_haplotypes(guide, sets, crista_cap)
        pams = list(pam_set)
        triples = [(guide, h, p) for h in haps for p in pams]
        scores = crista_fn(triples) if triples else []
        out["crista_worst"] = max(scores) if scores else None
        out["crista_approx"] = truncated
        out["crista_n_evaluated"] = len(triples)
    return out
