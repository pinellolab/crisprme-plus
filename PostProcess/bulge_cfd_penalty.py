"""Position-aware CFD bulge penalty (CRISPRme+ 2.7.11).

Background
----------
CFD (Doench 2016) ships a 152-entry gap matrix: 76 RNA-bulge keys (``rX:d-,p``) and
76 DNA-bulge keys (``r-:dX,p``), position-weighted (seed ≈ 0, PAM-distal ≈ 0.7-0.97).
But in CRISPRme's scoring path only the RNA-bulge half is ever applied: the DNA-bulge
caller trims the looped-out base off the 5' end and realigns to 20 nt, which collapses
the bulge to column 1 where ``calc_cfd`` finds no key and takes its ``except: pass``
free-pass (factor 1.0). Net effect: **every DNA bulge — and any RNA bulge at the
undefined terminal columns 1/21 — is unpenalized**, regardless of how close it sits to
the PAM. A PAM-proximal bulge (maximally disruptive) scores the same as none.

This module recovers the bulge's TRUE position from the display alignment (where the
gap is in place) and returns the multiplicative correction that applies Doench's own
``r-:dX`` / ``rX:d-`` factor at that position — activating the dead DNA-bulge half and
filling the undefined terminal columns by extrapolation (5' tip → lenient, at/into the
PAM → abolish). Mismatch scoring and in-range RNA bulges (already handled by calc_cfd)
are left byte-identical; only the previously free-passed bulges change.

The correction is a *multiplier* on the existing CFD: for the bulges calc_cfd already
free-passed (contributed 1.0), multiply by the real position factor; for in-range RNA
bulges calc_cfd already scored, contribute 1.0 (no double-count).
"""

SPACER_LEN = 20  # SpCas9 protospacer length; PAM is the trailing columns beyond it


def _terminal_factor(pos):
    """Doench's gap matrix is defined only for spacer positions 2..20. Extrapolate the
    two undefined terminal columns a bulge can land on:
      * 5' PAM-distal tip (pos <= 1): bulges there are well tolerated (pos-2 is already
        0.69-0.97) -> lenient.
      * at / into the PAM (pos >= SPACER_LEN+1): a bulge at the PAM-spacer junction
        abolishes recognition -> 0.
    Returns the factor, or None when the position is in-range (caller uses the matrix)."""
    if pos <= 1:
        return 0.9
    if pos >= SPACER_LEN + 1:
        return 0.0
    return None


def _revcom_base(b):
    return {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "-": "-"}.get(b.upper())


def bulge_cfd_correction(aligned_crrna, aligned_protospacer, mm_scores, pam_len=3):
    """Multiplicative CFD correction applied ON TOP of ``cfd_score_bulge_aware``'s base.

    ``aligned_crrna`` / ``aligned_protospacer`` are the DISPLAY alignment (same length,
    gaps '-' in place, PAM included as the trailing ``pam_len`` columns). A gap in the
    crRNA = DNA bulge (extra genomic base); a gap in the protospacer = RNA bulge (extra
    guide base). Returns a float in [0, 1]: 1.0 when nothing needs correcting.

    Convention (must match ``cfd_score_bulge_aware``):
      * DNA bulge (gap in crRNA): the base score DELETES that column, so calc_cfd applies
        NO DNA-bulge factor -> here we supply the FULL position-aware ``r-:dX,pos`` factor
        (with the undefined terminal columns extrapolated).
      * RNA bulge (gap in protospacer): the base score KEEPS that column, so calc_cfd
        already applied ``rX:d-,pos`` for the in-range columns (2..20) -> here we correct
        ONLY the undefined terminal columns (no double-count in range).
    """
    g = aligned_crrna.upper()
    d = aligned_protospacer.upper()
    n = min(len(g), len(d))
    spacer_cols = n - pam_len  # the PAM occupies the last pam_len columns
    corr = 1.0
    for c in range(spacer_cols):
        pos = c + 1  # 1-based from the 5' end (Doench frame)
        gap_in_crrna = g[c] == "-"
        gap_in_proto = d[c] == "-"
        if not (gap_in_crrna or gap_in_proto):
            continue
        if gap_in_crrna:
            # DNA bulge: supply the full position factor (base score deleted this column).
            base = _revcom_base(d[c])
            key = "r-:d" + str(base) + "," + str(pos)
            tf = _terminal_factor(pos)
            corr *= tf if tf is not None else mm_scores.get(key, _terminal_factor(SPACER_LEN + 1) or 0.0)
        else:
            # RNA bulge: calc_cfd ALREADY applied rX:d-,pos for in-range columns (2..20);
            # only the undefined terminal columns were free-passed -> correct just those.
            tf = _terminal_factor(pos)
            if tf is not None:
                corr *= tf
    return corr


def cfd_score_bulge_aware(aligned_crrna, aligned_protospacer, mm_scores, pam_scores,
                          calc_cfd_fn, pam_len=3, do_scores=True):
    """Position-aware, byte-stable replacement for the per-caller ``calc_cfd`` bulge paths.

    Takes the IN-PLACE display alignment (``aligned_crrna`` / ``aligned_protospacer`` same
    length, gaps '-' in place, PAM = trailing ``pam_len`` cols) and returns the CFD in
    [0, 1], scoring the alignment that is actually displayed/searched:

      1. DNA-bulge columns (gap in crRNA) are DELETED from BOTH strands, so calc_cfd scores
         the 20 guide bases at their true Cas9 positions (the historical front-trim dropped
         the 5' guide base and mis-paired every downstream column; the in-place '-' shifted
         them). calc_cfd therefore applies no DNA-bulge factor.
      2. RNA-bulge columns (gap in protospacer) stay in place, so calc_cfd applies the
         in-range ``rX:d-,pos`` factor exactly as before.
      3. The result is multiplied by ``bulge_cfd_correction`` which supplies the full
         position-aware DNA-bulge factor and fixes the undefined terminal (1 / PAM-proximal)
         columns that calc_cfd free-passed.

    Mismatch-only alignments (no gaps) are byte-identical to the previous
    ``calc_cfd(crRNA, proto[:-pam], proto[-2:])``; only bulge rows change.

    ``calc_cfd_fn`` is injected by the caller (``calc_cfd`` in each module) to avoid a
    circular import.
    """
    cr = aligned_crrna
    pr = aligned_protospacer.upper()
    if "-" in cr:  # DNA bulge(s): drop the looped-out genomic base(s) from both strands
        keep = [i for i, ch in enumerate(cr) if ch != "-"]
        gcr = "".join(cr[i] for i in keep)
        gpr = "".join(pr[i] for i in keep)
    else:
        gcr, gpr = cr, pr
    base = calc_cfd_fn(gcr, gpr[:-pam_len], gpr[-2:], mm_scores, pam_scores, do_scores)
    if base == 0 or base == -1:
        return base  # untranslatable IUPAC (0) or no-score sentinel (-1) short-circuit
    return base * bulge_cfd_correction(aligned_crrna, aligned_protospacer, mm_scores, pam_len)
