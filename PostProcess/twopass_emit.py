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
# Bound on the alt-combo brute-force used ONLY to repair the rare joint-PAM residual the
# greedy leaves PAM-invalid. 2^12: a window with <=12 variant columns is fully repaired;
# denser windows (already flagged as high-variant-density) skip the fallback.
_PAM_REPAIR_CAP = 1 << 12


def _pam_valid_min_mm(columns, ref_seq, real_target, guide_no_pam, revert,
                      pos_beg, pos_end, revcom_fn, pam_valid_fn, cap):
    """Exact min-mismatch VALID-PAM assignment over the window's alt combos, or None if
    no valid-PAM combo exists (or the product exceeds ``cap``). Each column may take its
    reference base or any candidate alt. Used as the fallback when the greedy + multi-pass
    leaves the rep PAM-invalid but slow (enumeration) still forms a valid-PAM off-target
    (e.g. a valid PAM needing the REFERENCE allele at a variant PAM column). Returns
    ``{"seq": [...], "chosen": {pos_c: candidate}}``."""
    import itertools
    opts, total = [], 1
    for col in columns:
        col_opts = [None] + list(col["candidates"])   # None = keep the reference base
        total *= len(col_opts)
        if total > cap:
            return None
        opts.append([(col["pos_c"], o) for o in col_opts])
    best = None   # (mm, seq, chosen)
    for combo in itertools.product(*opts):
        seq = list(ref_seq)
        chosen = {}
        for pc, o in combo:
            if o is not None:
                seq[pc] = o["alt"]
                chosen[pc] = o
        if not pam_valid_fn(seq):
            continue
        mm = aligned_mismatches(seq, real_target, guide_no_pam, revert,
                                pos_beg, pos_end, revcom_fn)
        if best is None or mm < best[0]:
            best = (mm, seq, chosen)
    return None if best is None else {"seq": best[1], "chosen": best[2]}


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

    def _pam_flag(seq):
        return _PAM_INVALID if (pam_valid_fn is not None and not pam_valid_fn(seq)) else _PAM_OK

    def _select(col):
        """Pick the best allele for ``col`` given the CURRENT greedy_seq (every other
        column fixed at its current value). Mutates greedy_seq[pc] to the choice. Returns
        the chosen candidate dict, or None when the reference base wins."""
        pc = col["pos_c"]
        greedy_seq[pc] = ref_seq[pc]                      # reference baseline for this col
        base_mm = aligned_mismatches(greedy_seq, real_target, guide_no_pam, revert,
                                     pos_beg, pos_end, revcom_fn)
        # sentinel lex on the reference makes any alt win a full (mm, pam) tie.
        best_key, best = (base_mm, _pam_flag(greedy_seq), _REF_SENTINEL), None
        for cand in col["candidates"]:
            greedy_seq[pc] = cand["alt"]
            m = aligned_mismatches(greedy_seq, real_target, guide_no_pam, revert,
                                   pos_beg, pos_end, revcom_fn)
            key = (m, _pam_flag(greedy_seq), cand["alt"])
            if key < best_key:
                best_key, best = key, cand
        greedy_seq[pc] = best["alt"] if best is not None else ref_seq[pc]
        return best

    chosen = {col["pos_c"]: _select(col) for col in columns}

    # A PAM needing a specific allele at 2+ variant columns (e.g. GMR = M:A AND R:G) may
    # not resolve in one left-to-right pass (one column is decided before the other is
    # set). Re-select each column given the others' current choices until the PAM is valid
    # or a fixpoint; each pass fixes >=1 PAM column so it converges in <= (#columns) passes.
    # Skipped when not PAM-aware (idempotent there) so pam_valid_fn=None stays single-pass.
    # Mismatch is the primary key, so a re-selection never raises the mismatch count
    # (position-independent) -- losslessness is preserved.
    if pam_valid_fn is not None:
        for _ in range(len(columns)):
            if pam_valid_fn(greedy_seq):
                break
            changed = False
            for col in columns:
                prev = chosen[col["pos_c"]]
                cur = _select(col)
                if (cur["alt"] if cur else None) != (prev["alt"] if prev else None):
                    changed = True
                chosen[col["pos_c"]] = cur
            if not changed:
                break

        # Fallback for the rare joint-PAM residual the coordinate-ascent can't resolve (a
        # valid PAM that needs the REFERENCE allele at a variant PAM column, or 2+
        # simultaneous non-lex flips): bounded brute-force for the exact min-mismatch
        # VALID-PAM assignment. Runs ONLY when the rep is still PAM-invalid, so it is off
        # the common path. A win here is a real off-target the greedy would otherwise drop.
        if not pam_valid_fn(greedy_seq):
            repaired = _pam_valid_min_mm(columns, ref_seq, real_target, guide_no_pam,
                                         revert, pos_beg, pos_end, revcom_fn,
                                         pam_valid_fn, _PAM_REPAIR_CAP)
            if repaired is not None:
                greedy_seq = repaired["seq"]
                chosen = {col["pos_c"]: repaired["chosen"].get(col["pos_c"])
                          for col in columns}

    carriers, info = set(), []
    for col in columns:                                  # rebuild in column order
        ch = chosen[col["pos_c"]]
        if ch is not None:
            carriers |= set(ch["carriers"])
            info.append(ch["info"])
    return {"seq": greedy_seq, "carriers": carriers, "info": info}


def _max_score_brute(columns, ref_seq, real_target, guide_no_pam, revert,
                     pos_beg, pos_end, revcom_fn, score_fn, cap, valid_fn=None):
    """EXACT max-``score_fn`` assignment over the window's alt combos (each column: its
    reference base or any candidate alt), or None if the product exceeds ``cap`` (or no
    combo satisfies ``valid_fn``). Used to make the worst-case CFD EXACT: CFD does NOT
    fully factorize (the PAM factor is a joint lookup over the PAM bases), so a per-column
    greedy is not guaranteed optimal. ``valid_fn(seq, chosen) -> bool`` restricts the max
    to EMITTABLE combos (mirroring the finalizer: a variant off-target must have >=1
    carrier + be in the mismatch budget + have a valid PAM) -- WITHOUT it the argmax can be
    a carrier-less combo the finalizer drops, so a weaker emittable combo (the true
    worst-possible carried off-target) is reported instead. Returns ``{"seq": list,
    "chosen": {pos_c: candidate}}``."""
    import itertools
    opts, total = [], 1
    for col in columns:
        col_opts = [None] + list(col["candidates"])   # None = keep the reference base
        total *= len(col_opts)
        if total > cap:
            return None
        opts.append([(col["pos_c"], o) for o in col_opts])
    best = None   # (score, seq, chosen)
    for combo in itertools.product(*opts):
        seq = list(ref_seq)
        chosen = {}
        for pc, o in combo:
            if o is not None:
                seq[pc] = o["alt"]
                chosen[pc] = o
        if valid_fn is not None and not valid_fn(seq, chosen):
            continue
        s = score_fn(seq)
        if best is None or s > best[0]:
            best = (s, seq, chosen)
    return None if best is None else {"seq": best[1], "chosen": best[2]}


def greedy_max_score(columns, ref_seq, real_target, guide_no_pam, revert,
                     pos_beg, pos_end, revcom_fn, score_fn, cap=_PAM_REPAIR_CAP,
                     valid_fn=None):
    """The strongest (worst) off-target the window can form under a score where HIGHER =
    worse (e.g. worst-case CFD) -- a min-MISMATCH representative does NOT maximize CFD (CFD
    is position-weighted, so a seed mismatch outweighs several distal ones). EXACT: over a
    small window (allele-combo product <= ``cap``) it brute-forces the true argmax, because
    CFD does NOT fully factorize -- the PAM factor is a JOINT lookup over the PAM bases, so
    a per-column greedy can stall at a local max and UNDER-state the worst case (measured on
    real chr22 double-bulge windows). Over a large (dense, separately-flagged) window it
    falls back to a greedy per-column argmax + bounded multi-pass. ``score_fn(seq_prerevert
    _list) -> float`` MUST already fold in the alignment + PAM. Returns ``{"seq": list,
    "carriers": set, "info": [[rsID,AF,snp], ...]}``."""
    exact = _max_score_brute(columns, ref_seq, real_target, guide_no_pam, revert,
                             pos_beg, pos_end, revcom_fn, score_fn, cap, valid_fn=valid_fn)
    if exact is not None:
        carriers, info = set(), []
        for col in columns:
            ch = exact["chosen"].get(col["pos_c"])
            if ch is not None:
                carriers |= set(ch["carriers"])
                info.append(ch["info"])
        return {"seq": exact["seq"], "carriers": carriers, "info": info}

    # dense window (> cap combos): greedy per-column argmax + bounded multi-pass (best-effort).
    greedy_seq = list(ref_seq)

    def _select(col):
        pc = col["pos_c"]
        greedy_seq[pc] = ref_seq[pc]                       # reference baseline for this col
        # maximize score: minimize the key (-score, is_ref, alt) -> higher score first,
        # an alt preferred over the reference on a full tie, then lex (deterministic).
        best_key, best = (-score_fn(greedy_seq), 1, _REF_SENTINEL), None
        for cand in col["candidates"]:
            greedy_seq[pc] = cand["alt"]
            key = (-score_fn(greedy_seq), 0, cand["alt"])
            if key < best_key:
                best_key, best = key, cand
        greedy_seq[pc] = best["alt"] if best is not None else ref_seq[pc]
        return best

    chosen = {col["pos_c"]: _select(col) for col in columns}
    prev = score_fn(greedy_seq)
    for _ in range(len(columns)):
        changed = False
        for col in columns:
            p = chosen[col["pos_c"]]
            c = _select(col)
            if (c["alt"] if c else None) != (p["alt"] if p else None):
                changed = True
            chosen[col["pos_c"]] = c
        cur = score_fn(greedy_seq)
        if not changed or cur <= prev:            # fixpoint or no improvement
            break
        prev = cur

    carriers, info = set(), []
    for col in columns:
        ch = chosen[col["pos_c"]]
        if ch is not None:
            carriers |= set(ch["carriers"])
            info.append(ch["info"])
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
