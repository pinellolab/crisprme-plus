#!/usr/bin/env python3
"""Pass-2 emission: the single worst-POSSIBLE off-target row per IUPAC window.

Design: docs/DESIGN_2.5.1_two_pass_fast_mode.md sections 2, 4. In fast mode the
post-analysis does NOT enumerate the 2^k IUPAC haplotype lattice (nor the observed
per-sample haplotypes). Instead, for each candidate window it emits ONE concrete
representative haplotype -- the greedy per-column MIN-MISMATCH allele choice, which
(mismatch being additive per column under the window's fixed bulge alignment) is the
EXACT argmin edit over all 2^k allele combinations, and therefore the maximum-CFD /
worst-possible off-target the window can form. This is the enumeration-free core that
replaces the intractable dense-panel post-analysis (49h+ on the 4x panel).

This is the SAME greedy representative the shipped code already builds for dense
(`CRISPRME_IUPAC_CAP`) windows (validated == brute-force argmin over all allele
combinations); fast mode simply applies it to EVERY window and skips the enumeration.

Kept STDLIB-only and free of ``new_simple_analysis`` module globals (which open
sys.argv[1] at import) so it is unit-testable in isolation, mirroring twopass_find /
twopass_score / twopass_cfd_exact. The alignment parameters (pos_beg / pos_end / the
reverse-complement function) are passed in explicitly by the caller.
"""


def aligned_mismatches(seq_prerevert, real_target, guide_no_pam, revert,
                       pos_beg, pos_end, revcom_fn):
    """Mismatch count of an ungapped, pre-revert candidate sequence against the guide.

    Byte-for-byte the semantics of ``new_simple_analysis._aligned_mm``: reverse-
    complement on the minus strand, re-insert bulge gaps at ``real_target``'s '-'
    positions, then compare the protospacer window (``[pos_beg:pos_end]``) to
    ``guide_no_pam`` (a '-' window base is never a mismatch)."""
    seq = revcom_fn("".join(seq_prerevert)) if revert else "".join(seq_prerevert)
    t = list(seq)
    for pos, ch in enumerate(real_target):
        if ch == "-":
            t.insert(pos, "-")
    window = t[pos_beg:pos_end]
    mm = 0
    for i, ch in enumerate(window):
        if i < len(guide_no_pam) and ch != "-" and ch.upper() != guide_no_pam[i]:
            mm += 1
    return mm


_PAM_INVALID = 1   # worse than _PAM_OK in the selection key (prefer a valid PAM on ties)
_PAM_OK = 0
_REF_SENTINEL = "\xff"   # sorts after every ACGT base -> an alt wins a full (mm,pam) tie


def greedy_worst_case(columns, ref_seq, real_target, guide_no_pam, revert,
                      pos_beg, pos_end, revcom_fn, pam_valid_fn=None):
    """Build the greedy worst-possible representative haplotype for one window.

    Args:
      columns: the variant columns of this window, each a dict:
        ``{"pos_c": int, "candidates": [{"alt": base, "carriers": set, "info": [...]}]}``
        (``pos_c`` indexes into ``ref_seq``; ``candidates`` are the admissible alt
        alleles at that column, each with its carrier sample set + [rsID, AF, snp_info]).
      ref_seq: the pre-revert reference sequence spanning the window.
      real_target / guide_no_pam / revert / pos_beg / pos_end / revcom_fn: the fixed
        bulge-alignment context (see ``aligned_mismatches``).
      pam_valid_fn: OPTIONAL ``callable(seq_prerevert_list) -> bool`` reporting whether the
        candidate produces a VALID PAM (mirrors the finalizer's ``pam_ok`` gate). When
        given, a mismatch-neutral column at a PAM position picks the PAM-*valid* allele
        instead of the lexicographic one -- WITHOUT this, a PAM-region variant (e.g. an
        ``R`` = A/G column where only G yields a valid NRG PAM) is decided by lex order,
        so ~half the time the greedy rep is PAM-invalid and the off-target is silently
        DROPPED (verified on real chr22 data: 141 PAM-creation off-targets missed). Never
        raises mismatch to gain a PAM (mismatch is the primary key), so it is lossless.

    Returns ``{"seq": <list of bases>, "carriers": <set>, "info": <list of [rsID,AF,snp]>}``.
    At each column pick the alt minimizing the selection key ``(mismatch, pam_invalid,
    alt)``: fewest mismatches first (additive per column under the fixed alignment, so the
    greedy equals the argmin edit over all 2^k combinations), then -- among mismatch ties
    -- a valid PAM, then the lexicographically-smaller alt (deterministic). An alt is
    preferred over the reference base on a full tie (keeps PAM-creating / present variants).
    Columns whose every alt only raises the mismatch count are left at the reference base."""
    greedy_seq = list(ref_seq)
    carriers, info = set(), []

    def _pam_flag(seq):
        return _PAM_INVALID if (pam_valid_fn is not None and not pam_valid_fn(seq)) else _PAM_OK

    for col in columns:
        pc = col["pos_c"]
        base_mm = aligned_mismatches(greedy_seq, real_target, guide_no_pam, revert,
                                     pos_beg, pos_end, revcom_fn)
        # key of KEEPING the reference base here (no alt); the sentinel lex makes any
        # alt win a full (mm, pam) tie, matching the legacy "prefer an alt on ties".
        best_key = (base_mm, _pam_flag(greedy_seq), _REF_SENTINEL)
        best = None  # (alt, carriers, info)
        for cand in col["candidates"]:
            trial = list(greedy_seq)
            trial[pc] = cand["alt"]
            m = aligned_mismatches(trial, real_target, guide_no_pam, revert,
                                   pos_beg, pos_end, revcom_fn)
            key = (m, _pam_flag(trial), cand["alt"])
            if key < best_key:
                best_key, best = key, (cand["alt"], cand["carriers"], cand["info"])
        if best is not None:
            greedy_seq[pc] = best[0]
            carriers |= set(best[1])
            info.append(best[2])
    return {"seq": greedy_seq, "carriers": carriers, "info": info}


def brute_force_min_mismatch(columns, ref_seq, real_target, guide_no_pam, revert,
                             pos_beg, pos_end, revcom_fn):
    """Reference oracle for the tests: the true minimum mismatch count over EVERY
    combination of one admissible allele (or the reference base) per column. Exponential
    in the number of columns -- test-only."""
    import itertools
    choices = []
    for col in columns:
        opts = [ref_seq[col["pos_c"]]] + [c["alt"] for c in col["candidates"]]
        choices.append([(col["pos_c"], o) for o in opts])
    best = None
    for combo in itertools.product(*choices):
        seq = list(ref_seq)
        for pc, base in combo:
            seq[pc] = base
        m = aligned_mismatches(seq, real_target, guide_no_pam, revert,
                               pos_beg, pos_end, revcom_fn)
        best = m if best is None else min(best, m)
    return best
