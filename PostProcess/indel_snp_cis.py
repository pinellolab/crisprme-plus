#!/usr/bin/env python
"""Cis co-occurrence of an indel + overlapping SNP alt(s) (feature/indel-snp).

Given an indel off-target whose aligned window overlaps one or more SNP alt
alleles (the alleles the target actually USES, determined by the caller from the
IUPAC-decoded alignment), decide -- per sample -- whether that sample carries the
indel AND all those SNP alts on the SAME haplotype (cis), and whether that cis
call is CONFIRMED (all contributing genotypes phased) or PUTATIVE.

Semantics mirror observed_haplotypes (CONFIRMED = every cis carrier reached the
set via a phased same-slot path; PUTATIVE if any carrier is unphased):
  * A sample is a candidate iff it carries the indel AND every used SNP alt.
  * Phased candidate: it is a cis carrier iff some haplotype slot h carries the
    indel alt AND every used SNP alt; cis-copy count = number of such slots.
  * Unphased candidate (any contributing GT has '/' or is haploid/single-token):
    slot alignment can't be proven -> counted as ONE PUTATIVE cis copy and forces
    the overall phase_state to PUTATIVE.

Genotypes are the biallelic normalized strings produced by
build_indel_genotypes.normalize_gt_for_alt / the SNP tier ("1|0","0|1","1|1",
"0/1", "1"). STDLIB ONLY. Pure + unit-tested (test_indel_snp_cis.py).

Since 2.8.1 this module also turns an overlaid indel target (IUPAC codes at the SNP
columns of the fake-indel window) into the rows its carriers really have:
resolve_target gives every carrier slot its own allele at each code column and keeps
the sequences that still hold the PAM and fit the search budget (see the section
"Resolving an overlaid indel target per carrier haplotype" below;
test_indel_haplotype_resolution.py).
"""

import bisect
import itertools

import site_alleles
from overlay_indel_snps import map_fake_offset_to_real

CONFIRMED = "CONFIRMED"
PUTATIVE = "PUTATIVE"


def _slots(gt):
    """Return (slot_list, phased) for a normalized biallelic GT.

    A '|' is phased; '/' or a single haploid token is treated as UNPHASED (we
    cannot prove same-slot cis for it).
    """
    if "|" in gt:
        return gt.split("|"), True
    if "/" in gt:
        return gt.split("/"), False
    return [gt], False  # haploid / single token -> unconfirmed phase


def cis_cooccurrence(indel_gt_by_sample, snp_gts_by_sample):
    """Compute the cis co-occurrence of the indel + all used SNP alts.

    Args:
      indel_gt_by_sample: {sample_id: normalized_gt} for the indel.
      snp_gts_by_sample: list of {sample_id: normalized_gt}, ONE dict per used SNP
        alt (empty list => just the indel, trivially "cis").

    Returns (cis_carriers:set[str], phase_state:str, ac:int):
      * cis_carriers -- samples carrying the indel + all SNP alts in cis (phased
        same-slot) or putatively-cis (unphased but carries all).
      * phase_state -- CONFIRMED iff every cis carrier is phased-same-slot, else
        PUTATIVE.
      * ac -- cis haplotype copy count (allele count) for joint AF = ac / AN.
    """
    cis_carriers = set()
    ac = 0
    any_putative = False

    for sample, igt in indel_gt_by_sample.items():
        # must carry every used SNP alt too
        snp_slot_data = []
        carries_all = True
        for snp in snp_gts_by_sample:
            g = snp.get(sample)
            if g is None:
                carries_all = False
                break
            snp_slot_data.append(_slots(g))
        if not carries_all:
            continue

        islots, iphased = _slots(igt)
        all_phased = iphased and all(ph for _, ph in snp_slot_data)

        if all_phased:
            # count slots where the indel alt AND every SNP alt co-occur
            copies = 0
            for h in range(len(islots)):
                if islots[h] != "1":
                    continue
                if all(
                    (sl[h] if h < len(sl) else "0") == "1" for sl, _ in snp_slot_data
                ):
                    copies += 1
            if copies:
                cis_carriers.add(sample)
                ac += copies
        else:
            # unphased: carries the indel + all SNP alts, but slot alignment
            # unprovable -> ONE putative cis copy, forces PUTATIVE.
            cis_carriers.add(sample)
            ac += 1
            any_putative = True

    phase_state = PUTATIVE if any_putative else CONFIRMED
    return cis_carriers, phase_state, ac


def joint_af(ac, an):
    """Joint indel+SNP cis allele frequency (ac / an), or 0.0 if an falsy."""
    return (ac / an) if an else 0.0


def putative_joint_af(marginal_afs):
    """PUTATIVE joint allele frequency when genotypes are UNAVAILABLE (sites-only /
    aggregate panels such as the mega): with no per-sample genotypes we cannot count
    cis carriers, but the cis co-occurrence frequency of a set of variants can never
    exceed the frequency of the RAREST participating allele -- so report min(marginal
    AFs) as a conservative upper bound (matches the PUTATIVE = upper-bound semantics).
    Returns 0.0 if the set is empty or any marginal AF is missing (None), so we never
    emit a misleadingly-high joint frequency from incomplete information."""
    afs = list(marginal_afs)
    if not afs or any(a is None for a in afs):
        return 0.0
    return min(afs)


# The IUPAC ambiguity chars an overlaid SNP position can carry in the target seq, and Q,
# the code of a position whose reference base and SNP ALTs hold all four bases (N is a
# gap, never a variant).
_IUPAC = set("RYSWKMBDHVQryswkmbdhvq")

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(base):
    """Single-base complement (identity for non-ACGT, e.g. IUPAC/N)."""
    return base.translate(_COMP)


def build_offset_to_real(dna_aligned, fake_start, strand, indel_fake_start,
                         indel_start, ref, alt):
    """Map each aligned column of an indel off-target's DNA to a real 0-based ref
    position (or None for a bulge/inserted column).

    Args:
      dna_aligned: the aligned DNA target (search col 2; '-' at RNA-bulge columns).
      fake_start:  the target's leftmost fake position (search col 4, forward coord).
      strand:      '+' or '-' (search col 6).
      indel_fake_start: the indel's fake-region start (log FAKEPOS start).
      indel_start, ref, alt: the indel's real 0-based start + REF/ALT (for the
        piecewise fake->real map).

    Returns a list ``real_of`` aligned to ``dna_aligned`` (real 0-based pos or None).

    Each non-'-' DNA base consumes one forward fake position: for '+' strand the
    leftmost base is at ``fake_start`` and positions increase left->right; for '-'
    strand the displayed DNA is the reverse complement, so the leftmost base is at
    the HIGHEST forward position and they decrease left->right. The forward fake
    offset is then run through the validated piecewise indel-boundary map.

    VALIDATED: 2189/2189 per-base reference matches (1582 '+' / 607 '-') against
    105 real chr22 indel off-targets -- see the commit's validation note.
    """
    base_cols = [j for j, ch in enumerate(dna_aligned) if ch != "-"]
    n = len(base_cols)
    real_of = [None] * len(dna_aligned)
    for i, j in enumerate(base_cols):
        fpos = fake_start + i if strand == "+" else fake_start + (n - 1 - i)
        real_of[j] = map_fake_offset_to_real(
            fpos - indel_fake_start, indel_start, len(ref), len(alt)
        )
    return real_of


def used_snps_for_target(guide_seq, target_seq, offset_to_real, snp_at,
                         pam_motif_at=None):
    """Identify the SNP alt alleles an indel off-target actually USES.

    Args:
      guide_seq, target_seq: the aligned spacer and DNA target (equal length; may
        contain '-' gaps for bulges; target_seq carries IUPAC codes at overlaid
        SNP positions).
      offset_to_real(j): map aligned column j -> real reference position, or None
        for a bulge/inserted column with no reference position. (The caller builds
        this to handle bulge gaps, strand, and the piecewise indel-boundary shift.)
      snp_at(real_pos): -> (ref_base, alt_base, rsID, {sample: normalized_gt}) or
        None if no SNP at that position.
      pam_motif_at: {column: PAM motif character} for the PAM columns of the target
        (pam_columns). Without it a PAM column is skipped, as before 2.8.1.

    Returns a list of (real_pos, rsID, alt_gt_by_sample) for every column where the
    target carries an IUPAC ambiguity AND the off-target REQUIRES the SNP's alt
    allele (a genuine SNP+indel co-occurrence), not merely tolerates the reference:
      * a spacer column uses the alt iff the guide base equals the alt, not the REF;
      * a PAM column (the guide holds 'N' there) uses the alt iff the alt satisfies
        the PAM motif there and the REF does not. The engine kept the target only
        because some allele of the code satisfies the motif, so a PAM the reference
        cannot form was created by that SNP.
    Ref-satisfied ambiguities are skipped (the SNP is incidental, not enabling).
    """
    used = []
    for j, (g, t) in enumerate(zip(guide_seq, target_seq)):
        if t not in _IUPAC:
            continue
        gu = g.upper()
        motif = pam_motif_at.get(j) if pam_motif_at else None
        if gu == "-" or (gu == "N" and motif is None):
            continue
        real = offset_to_real(j)
        if real is None:
            continue
        snp = snp_at(real)
        if snp is None:
            continue
        ref_base, alt_base, rsid, gt = snp
        if motif is not None:
            allowed = _CODE.get(motif.upper(), "")
            if alt_base.upper() in allowed and ref_base.upper() not in allowed:
                used.append((real, rsid, gt))
            continue
        # USES the alt iff the guide matches the alt but not the ref
        if gu == alt_base.upper() and gu != ref_base.upper():
            used.append((real, rsid, gt))
    return used


# ---------------------------------------------------------------------------
# Resolving an overlaid indel target per carrier haplotype (2.8.1)
# ---------------------------------------------------------------------------
# The SNP overlay writes each enriched-genome IUPAC code onto the flanks of a
# fake-indel window, and the engine matches a code when ANY of its alleles fits the
# guide (or the PAM). A raw indel target therefore stands for every combination of
# the window's SNP alleles at once, and only the genotypes say which combinations
# exist. Up to 2.8.0 the post-analysis turned each code into one base column by
# column -- the guide's base when it was an allele, else the code's first allele --
# without looking at any genotype: one sample's indel was reported with SNP alleles
# that only other samples carry (often scored high although no individual carries
# the sequence), and a code at a PAM column always read as its first allele (an R
# at NRG's last base read A, a PAM the search never matched).
#
# Here every carrier slot of the indel gets its OWN allele at each code column
# (the SNP alt the slot carries, else the reference base); slots that end up with
# the same sequence form one row, and a sequence whose PAM fails the motif or whose
# recounted edits exceed the search budget is dropped. A carrier whose phase cannot
# be read off its genotypes (unphased, haploid against diploid, a missing allele, no
# phased indel genotype) is resolved to the strongest target its own alleles allow
# and marked PUTATIVE: it is never listed for an allele it does not carry. With no
# genotypes at all (a sites-only panel) every carrier gets that worst case over the
# code's alleles, PUTATIVE, as the co-occurrence fallback of 2.5.2.

# bases each IUPAC character stands for (upper case; a lowercase code is the same code,
# its case only marks a mismatch in the engine's DNA column)
_CODE = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT", "Q": "ACGT",
}
_BASE_ORDER = "ACGT"


def ambiguous_columns(target_seq):
    """The columns of an aligned target that hold an IUPAC ambiguity code (either
    case): the SNPs the overlay put on the fake-indel window there."""
    return [j for j, ch in enumerate(target_seq) if ch in _IUPAC]


def pam_columns(target_len, motif, pam_at_start=False):
    """{column: motif character} for the PAM window of an aligned target of
    ``target_len`` columns: the last ``len(motif)`` columns for a 3' PAM, the first
    ones for a 5' PAM. The engine never puts a bulge inside the PAM, so the window is
    exactly the motif's width."""
    n = len(motif)
    start = 0 if pam_at_start else target_len - n
    return {start + k: m for k, m in enumerate(motif.upper())}


class SnpColumn(object):
    """One code column of an indel target, in the target's own orientation.

      col    the column in the aligned target
      real   its 0-based reference position (forward strand)
      pos1   the 1-based position the registry / genotype tiers key it by
      ref    the reference base as the target shows it (complemented on '-')
      ref_fwd the forward-strand reference base
      alts   {base as the target shows it: (forward alt, {sample: genotype})}: every
             alt allele the genotype source records there, with its carriers
      code   the bases of the IUPAC code the engine matched there
      known  False when no genotype source covers the position: every carrier is
             then unresolved there (a sites-only panel, a record a tier lacks)
    """

    __slots__ = ("col", "real", "pos1", "ref", "ref_fwd", "alts", "code", "known")

    def __init__(self, col, real, ref_fwd, alts, code, strand, known=True):
        self.col = col
        self.real = real
        self.pos1 = real + 1
        self.ref_fwd = ref_fwd.upper()
        self.ref = complement(self.ref_fwd) if strand == "-" else self.ref_fwd
        self.alts = {}
        for alt_fwd, carriers in (alts or {}).items():
            a = alt_fwd.upper()
            shown = complement(a) if strand == "-" else a
            self.alts[shown] = (a, carriers)
        self.code = frozenset(_CODE.get(code.upper(), ""))
        self.known = known


# (allele tokens, kind) of a genotype: kind is "phased" ('|'), "unphased" ('/') or
# "haploid" (one token, e.g. a male chrX call)
_gt_tokens = site_alleles.gt_tokens


def kept_bases(ref, alt):
    """How many bases REF and ALT share at the start (the VCF anchor base): the REF
    bases a haplotype carrying the variant keeps. The rest of REF is deleted or
    replaced."""
    c = 0
    while c < min(len(ref), len(alt)) and ref[c].upper() == alt[c].upper():
        c += 1
    return c


def shared_tail(ref, alt):
    """How many bases REF and ALT share at the end, past the ones they share at the
    start: a record left unnormalized (TCA>TA, from a multi-allelic line) keeps them,
    and a haplotype carrying it holds them unchanged. bcftools norm leaves none."""
    c = kept_bases(ref, alt)
    t = 0
    while t < min(len(ref), len(alt)) - c and ref[-1 - t].upper() == alt[-1 - t].upper():
        t += 1
    return t


def changes_span(pos1, ref, alt, lo, hi, gone=None):
    """True iff a variant (VCF POS / REF / ALT) changes the reference bases lo..hi
    (0-based, inclusive) of a haplotype that carries it: a deleted or replaced base
    inside, or bases inserted between two bases of it. The bases REF and ALT share
    at the start (the VCF anchor base) or at the end are kept, so they change nothing.
    ``gone`` (first, last; 0-based, inclusive): bases the haplotype has lost already,
    to the indel of the window being read: losing them again changes nothing."""
    a = pos1 - 1
    c = kept_bases(ref, alt)
    t = shared_tail(ref, alt)
    if len(ref) - t > c:  # bases a+c .. a+len(ref)-t-1 are deleted or replaced
        first, last = a + c, a + len(ref) - t - 1
        parts = [(first, last)] if gone is None else [
            (first, min(last, gone[0] - 1)), (max(first, gone[1] + 1), last)]
        if any(f <= l and f <= hi and l >= lo for f, l in parts):
            return True
    if len(alt) - t > c:  # bases inserted between a+c-1 and a+c
        return lo <= a + c - 1 and a + c <= hi
    return False


def target_ref_span(dna_aligned, real_of, indel_pos1, ref):
    """(lo, hi): the reference bases (0-based, inclusive) an indel off-target covers
    on a haplotype carrying its indel, or None for a target with no base.

    ``real_of`` (build_offset_to_real) maps the ALT bases to no reference position,
    the anchor base included, so a target that ends inside them would stop one
    base short of the indel. A target holding any ALT base also holds the indel's
    REF place (the anchor; for a deletion, the bases the ALT replaces), so the span
    takes in REF: a second variant there, or inserted next to the anchor, is inside
    the target."""
    reals = [r for r in real_of if r is not None]
    if any(r is None and c != "-" for c, r in zip(dna_aligned, real_of)):
        reals += [indel_pos1 - 1, indel_pos1 - 2 + len(ref)]
    if not reals:
        return None
    return min(reals), max(reals)


def window_offsets(dna_aligned, fake_start, indel_fake_start):
    """(first, last): the offsets (inclusive) of a target's bases in its fake-indel
    window. The aligned DNA's bases are contiguous on the window, a gap being an RNA
    bulge, and ``fake_start`` is the leftmost of them on either strand."""
    first = fake_start - indel_fake_start
    return first, first + len(dna_aligned) - dna_aligned.count("-") - 1


def holds_indel(first, last, indel_pos1, indel_start, ref, alt):
    """True iff a target at the window offsets first..last (window_offsets) holds the
    window's indel: a base the ALT adds or replaces, or the two bases on either side
    of the deleted ones. Anything else lies wholly in a flank: plain reference (or SNP)
    sequence that does not need the indel, at a coordinate the indel shifts.
    ``indel_start`` is the real 0-based start of the window (the log's), so the REF's
    first base is at offset indel_pos1 - 1 - indel_start (overlay_indel_snps.LEFT_FLANK
    on a full window)."""
    alt_at = indel_pos1 - 1 - indel_start
    added_from = alt_at + kept_bases(ref, alt)        # first base the ALT adds or replaces
    added_to = alt_at + len(alt) - shared_tail(ref, alt)  # and the base past the last
    if added_from < added_to:
        return first < added_to and last >= added_from
    return first < added_from and last >= added_from  # the bases either side of a deletion


def reference_start(offset, indel_pos1, indel_start, ref, alt):
    """The 0-based reference coordinate an indel row starts at when its first aligned
    base (the leftmost, on either strand) is at ``offset`` of a fake-indel window whose
    real 0-based start is ``indel_start``.

    The window is the reference up to the indel, its ALT, then the reference past its
    REF, so ``indel_start + offset`` is the base's place on the haplotype that carries
    the indel. That is its reference coordinate for a flank base before the indel, for
    a base REF and ALT share at the start (the anchor) and for an ALT base standing in
    a REF base's place (an MNV's, or a base an indel replaces). A base an insertion
    adds has no reference coordinate: there, the haplotype's coordinate is a base of
    the reference 1 to len(ALT) - len(REF) bases past the insertion's anchor, which the
    target does not hold. The row starts instead at the reference base the added bases
    follow (a normalized insertion's anchor: its VCF POS, 1-based), the place
    target_ref_span already counts as the target's. A base past the change (a shared
    tail, the downstream flank) is the reference base it is, len(REF) - len(ALT) on.

    Up to 2.8.1 every row was placed at ``indel_start + offset``, so a row whose first
    base is one an insertion adds started past the insertion, its REF text was read
    there, and its Variant_info_spacer+PAM put the insertion, before its start, at
    position 1 (the - strand: the last)."""
    alt_at = indel_pos1 - 1 - indel_start
    k = offset - alt_at
    c = kept_bases(ref, alt)
    t = shared_tail(ref, alt)
    if k < 0 or len(ref) == len(alt) or k < c:
        return indel_start + offset  # a flank base, an MNV's, or a base REF and ALT share
    if k >= len(alt) - t:
        return indel_start + offset - len(alt) + len(ref)  # past the change
    j = k - c  # a base of the ALT's changed middle
    if j < len(ref) - c - t:
        return indel_start + offset  # it stands in the place of the REF base it replaces
    # a base the ALT adds: the place of the last REF base before it
    return indel_pos1 - 2 + len(ref) - t


# ---------------------------------------------------------------------------
# Counting along an indel row's target (2.10.0)
# ---------------------------------------------------------------------------
# A SNP row's target is reference bases, and Variant_info_spacer+PAM numbers each of
# its variants by its distance from the row's start along the reference
# (resultIntegrator). An indel row's target is bases of the haplotype carrying its
# indel: an insertion's bases are in it, a deletion's are not. Its variants are
# numbered the same way along that haplotype: from the target's first base on it to
# each variant's base on it. Up to 2.8.1 they were numbered along the reference, so a
# variant past the indel was off by len(ALT) - len(REF), also on a row starting on a
# base an insertion adds, which 2.8.1 placed at that base's place on the haplotype;
# counted along the reference from the insertion's anchor, where 2.10.0 places such a
# row, it would be off by the inserted bases the target holds, less one.


def haplotype_place(pos1, indel_pos1, ref, alt):
    """The place (1-based), on a haplotype carrying the indel (VCF POS ``indel_pos1``,
    ``ref``, ``alt``), of the reference base at ``pos1`` (1-based): its own before the
    bases the indel changes and on the bases REF and ALT share at the start,
    len(ALT) - len(REF) on past them (the bases they share at the end included). A base
    the indel deletes or replaces is on no haplotype carrying it: None."""
    a = indel_pos1 - 1
    q = pos1 - 1
    if q < a + kept_bases(ref, alt):
        return pos1
    if q >= a + len(ref) - shared_tail(ref, alt):
        return pos1 + len(alt) - len(ref)
    return None


def haplotype_base_kind(h, indel_pos1, ref, alt):
    """("ref", p) for the base at place ``h`` (0-based) of a haplotype carrying the indel
    when it is reference position p (a flank base, a base REF and ALT share at either
    end); ("alt", p) for an ALT base standing in the place of the REF base p it replaces
    (an MNV's, a complex record's); ("added", None) for a base the ALT adds, which has no
    reference position."""
    a = indel_pos1 - 1
    c = kept_bases(ref, alt)
    t = shared_tail(ref, alt)
    if h < a + c:
        return "ref", h
    if h >= a + len(alt) - t:
        return "ref", h - len(alt) + len(ref)
    if h - a - c < len(ref) - c - t:
        return "alt", h
    return "added", None


def row_first_base(strand, start, alt_aligned, ref_aligned, indel, snps=()):
    """The place (0-based) of an indel row's first aligned base (the leftmost, on either
    strand) on the haplotype carrying its indel ``indel`` (VCF POS, REF, ALT) and its SNP
    alts ``snps`` [(POS, REF, ALT)], from what the row shows: its Start_coordinate
    ``start``, its aligned ALT (the target) and its aligned REF, which holds the reference
    bases from ``start`` on, as many as the target holds (alignRefFromVar).

    That place is ``start`` unless the first base is one the indel adds: such a row
    starts at the base the added bases follow (reference_start), and its first base is
    the k-th added base, k >= 1, at ``start`` + k. The row does not say k, and the bases
    do: the target is the haplotype's bases from its first base on, and the haplotype
    from ``start`` on is the REF text's first base when it is a reference one, the ALT's
    bases, then the REF text's bases past the indel, with the SNP alts the row names. The
    smallest k at which they are the target's is taken: two places can hold the same
    target only in a repeat the target lies wholly in, where they hold the same bases (the
    smaller is the base itself on a row starting on a reference base). None of them,
    which a row the post-analysis wrote never gives: ``start``.

    Cost: the haplotype is cut from the REF text and the ALT in slices and the target
    found in it with ``str.find``, so no Python step is taken per base of the ALT. Built
    base by base (haplotype_base_kind for each of its len(ALT) - len(REF) + 23 bases, then
    a slice compared at each k), a row at the anchor of one of the 402 insertions of 53.9
    to 54.4 kb that HPRC holds at chr22:10690626 took 83 ms, and the merge
    (remove_contiguous_samples) places each of the 16,058 rows of the README guide's chr22
    search that start inside them."""
    ipos, iref, ialt = indel
    if len(ialt) <= len(iref):
        return start  # no base added
    t = shared_tail(iref, ialt)
    if start != ipos - 2 + len(iref) - t:
        return start  # the row starts before the added bases' anchor
    target = alt_aligned.replace("-", "").upper()
    refb = ref_aligned.replace("-", "").upper()
    if not target or len(refb) != len(target):
        return start
    if strand == "-":
        target, refb = _revcomp(target), _revcomp(refb)
    n, more = len(target), len(ialt) - len(iref)
    # The haplotype's n + more bases from ``start`` on, by haplotype_base_kind: ``start``
    # holds the ALT's base number ``at`` (from 0; len(REF) - t - 1). It is the REF text's
    # first base when it is one REF and ALT share at the start (at < c, a plain insertion),
    # and an ALT base standing in a REF base's place otherwise; then come the ALT's bases
    # up to its shared tail, then the REF text's bases from its second on (the tail
    # included).
    c = kept_bases(iref, ialt)
    at = start - ipos + 1
    hap = (refb[0] if at < c else "") + ialt[max(c, at):len(ialt) - t].upper() + refb[1:]
    if snps:
        hap = list(hap)
        for pos1, _r, alt in snps:
            place = haplotype_place(pos1, ipos, iref, ialt)
            if place is not None and 0 <= place - 1 - start < len(hap):
                hap[place - 1 - start] = alt.upper()
        hap = "".join(hap)
    k = hap.find(target)
    return start + k if 0 <= k <= more else start


def variant_places(strand, start, alt_aligned, ref_aligned, tokens):
    """(place of the row's first base, [place of each variant's POS base]) on the
    haplotype carrying the row's indel, for the Variant_info_spacer+PAM numbering:
    ``tokens`` are the row's variants [(POS, REF, ALT)], forward strand, as
    Variant_info_genome names them, the indel's first. A row whose variants have no indel
    (a SNP row, an MNV row: no base added or lost) has the reference's: (``start``, their
    POS). An indel's own POS base, its anchor, is before the bases it changes, so it keeps
    its place; a SNP past the indel is len(ALT) - len(REF) on (haplotype_place); and the
    first base is row_first_base's."""
    indel = next((tk for tk in tokens if len(tk[1]) != len(tk[2])), None)
    if indel is None:
        return start, [tk[0] for tk in tokens]
    snps = [tk for tk in tokens if tk is not indel]
    places = [indel[0] if tk is indel else
              (haplotype_place(tk[0], *indel) or tk[0]) for tk in tokens]
    first = row_first_base(strand, start, alt_aligned, ref_aligned, indel, snps) if snps \
        else start
    return first, places


def _revcomp(seq):
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Where an indel row's PAM is judged on the reference (2.10.0)
# ---------------------------------------------------------------------------
# PAM_creation asks whether the reference holds a valid PAM where the row's PAM is. On an
# indel row that place depends on which of the target's bases are the indel's, and in a
# repeat that is not fixed: the VCF writes an insertion at the left end of its repeat
# (bcftools norm left-aligns it), but the same haplotype holds the same bases with the
# inserted sequence rotated through the repeat and written further right, and a deletion
# likewise. Up to 2.8.1 the PAM was judged where the VCF writes the indel, so the answer
# depended on that choice: an insertion of a repeat unit next to a reference PAM, written
# left of it, puts inserted bases in a PAM the reference holds (created, wrongly), and a
# row that starts on a reference base read its REF alignment's window, which is not the
# reference PAM beside the protospacer when the indel lies between them (on the - strand
# with a 3' PAM, the + strand with a 5' one, where the row starts at the PAM's far end: an
# insertion's bases, or a deletion's, shift that window off the PAM's place).
#
# The rule reads every placement of the indel on the row's haplotype -- the reference
# with the row's own variants: its indel and the SNP alts it names -- and judges the PAM
# where the reference holds it, in that order:
#   1. a placement that puts the whole PAM on reference bases (no inserted base, no
#      deletion's junction inside it): the reference bases it sits on there. Without a
#      SNP alt of the row in the PAM they are the PAM's own bases, so not created;
#   2. else a placement that puts the protospacer's base that touches the PAM on a
#      reference base: the reference PAM beside it, on the PAM's side (the PAM's place
#      beside that protospacer). A deletion's PAM always has one, its bases all being
#      reference bases;
#   3. else (an insertion that holds that base in every placement, the protospacer's end
#      and the PAM's start or the whole PAM): no reference place, created.
# The placements of one kind read the same reference bases (but for a SNP alt of the row in
# a PAM that a repeat spans: placed before the indel and past it, the PAM sits on two
# places, and it is not created if either holds the motif), so the answer depends on the
# haplotype alone, not on where in its repeat the VCF writes the indel.
#
# A complex record replaces REF bases and adds bases or loses some (A>T plus inserted
# bases: HPRC writes many insertions so). It is read as haplotype_base_kind and
# reference_start read it: its first min(len(REF), len(ALT)) changed bases stand in the
# places of the REF bases they replace, as an MNV's, and its other ALT bases are added
# after them (or its other REF bases deleted). The replaced bases are read as a SNP alt of
# the row is, the reference base each replaces, and the added (deleted) bases are placed
# as an insertion's (a deletion's), so a PAM made of a complex record's added bases is
# judged as one made of an insertion's. An MNV moves no base: its PAM is read on the
# reference under it, by the caller.


def _with_alts(genome, lo, hi, alts):
    """genome[lo:hi], upper case, with the SNP alts ``alts`` {0-based position: base} in
    it (the reference as the row's haplotype holds it outside its indel); None past
    either end of the genome."""
    if lo < 0 or hi > len(genome):
        return None
    seq = genome[lo:hi].upper()
    hits = [q for q in alts if lo <= q < hi]
    if not hits:
        return seq
    seq = list(seq)
    for q in hits:
        seq[q - lo] = alts[q]
    return "".join(seq)


def pam_reference_place(genome, indel, snps, first, n, k, pam_right):
    """Where an indel row's PAM is judged on the reference (the rule above): ("ref", los),
    the k reference bases genome[lo:lo + k] (forward strand) of each place ``lo`` the rule
    finds -- the PAM is not created if one of them holds the motif; ("created", ()), no
    placement of the indel gives the PAM a place on the reference, so the indel created
    it; None for a record that adds no base and loses none (an MNV), which the caller
    judges on the reference under the PAM, as 2.8.1.

    A PAM that placements put on the reference has one place there or two -- before the
    indel and past it, which hold the same bases --, and the first is the VCF's placement's
    when it has one. A SNP alt of the row in the PAM, or a base a complex record replaces,
    is read as the reference base it replaces, and then the two can differ: both are given
    (a row without SNP alts or replaced bases gets the first only). The PAM beside the
    protospacer has one place.

    ``indel`` is the row's indel (VCF POS, REF, ALT); ``snps`` [(POS, REF, ALT)] the SNP
    alts the row names next to it; ``first`` the target's first base (the leftmost,
    forward) on the haplotype carrying the indel, in the fake-indel window's coordinates
    carried to the chromosome (the reference's before the indel, its ALT written whole
    from POS, then the reference's shifted by len(ALT) - len(REF)); ``n`` its bases;
    ``k`` the PAM's; ``pam_right``: the PAM is the target's right end, forward (the +
    strand with a 3' PAM, the - strand with a 5' one).

    The indel's changed bases start at POS - 1 + the bases REF and ALT share at the start
    (the anchor); a complex record's first min(gone, added) of them stand in REF bases'
    places, and e is the base past those. An insertion's m added bases lie at [e, e + m)
    of the haplotype and can sit s bases further right while the haplotype's next s bases
    repeat them (H[e + i] = H[e + m + i]), or further left likewise; a deletion's d bases,
    deleted from the reference at [e, e + d), leave the haplotype's junction before its
    base e, which can move while the reference repeats them (R[e + i] = R[e + d + i]).
    Those comparisons read the reference with the row's SNP alts and the bases a complex
    record replaces, so a placement keeps the same haplotype with the same variants. With
    the PAM at [p0, p0 + k) of the haplotype:
      * an insertion placed at x leaves the PAM on the reference, as genome[p0:] when
        x >= p0 + k and as genome[p0 - m:] when x + m <= p0; otherwise the protospacer's
        base by the PAM is on the reference when x >= p0 (PAM at the right) or
        x + m <= p0 + k (at the left), and the PAM beside it is genome[p0:p0 + k] or
        genome[p0 - m:p0 - m + k];
      * a deletion's junction at x leaves the PAM on the reference, as genome[p0:] when
        x >= p0 + k and as genome[p0 + d:] when x <= p0; otherwise the PAM beside the
        protospacer is genome[p0:p0 + k] (PAM at the right) or genome[p0 + d:p0 + d + k].
    Only the shifts those tests need are read, each one slice comparison, so a long
    insertion, or a long run of replaced bases, adds no Python step per base."""
    ipos, iref, ialt = indel
    c = kept_bases(iref, ialt)
    t = shared_tail(iref, ialt)
    gone, added = len(iref) - c - t, len(ialt) - c - t
    if gone == added:
        return None  # an MNV: its ALT bases stand in REF bases' places, none added or lost
    # a complex record's bases that stand in REF bases' places, at [rep_at, e): read as the
    # row's SNP alts are; its other bases are an insertion's (a deletion's) at e
    r = min(gone, added)
    rep_at, rep = ipos - 1 + c, ialt[c:c + r].upper()
    e = rep_at + r
    gone, added = gone - r, added - r
    p0 = first + n - k if pam_right else first
    alts = {}
    for p, _r, a in snps or ():
        alts[p - 1] = a.upper()

    def hap_ref(lo, hi):  # the reference as the row's haplotype holds it outside the indel
        seq = _with_alts(genome, lo, hi, alts)
        if seq is None or not r or hi <= rep_at or lo >= e:
            return seq
        a, b = max(lo, rep_at), min(hi, e)
        return seq[:a - lo] + rep[a - rep_at:b - rep_at] + seq[b - lo:]

    if added:
        m = added
        ins = ialt[c + r:c + r + m].upper()

        def moves(s):  # the inserted bases can sit s bases further right (left if s < 0)
            if s > 0:
                if hap_ref(e, e + 1) != ins[0]:
                    return False  # not a repeat on that side (no slice read)
                x = hap_ref(e, e + s)
                return x is not None and (ins + x)[:s] == x
            if hap_ref(e - 1, e) != ins[-1]:
                return False
            y = hap_ref(e + s, e)
            return y is not None and (y + ins)[s:] == y

        # the PAM before the insertion (x >= p0 + k: genome[p0:]) or past it
        # (x + m <= p0: genome[p0 - m:]), where the VCF has it, then by moving it
        los = [p0] if e >= p0 + k else [p0 - m] if e + m <= p0 else []
        if not los or alts or r:
            if p0 not in los and moves(p0 + k - e):
                los.append(p0)
            if p0 - m not in los and moves(p0 - m - e):
                los.append(p0 - m)
        if los:
            return "ref", tuple(los)
        if pam_right:
            if e >= p0 or moves(p0 - e):
                return "ref", (p0,)
        elif e + m <= p0 + k or moves(p0 + k - e - m):
            return "ref", (p0 - m,)
        return "created", ()
    d = gone

    def slides(lo, s):  # the deleted bases can sit at lo .. lo + d - 1 instead
        x = hap_ref(lo, lo + s)
        z = hap_ref(lo + d, lo + d + s)
        return x is not None and z is not None and x == z

    los = [p0] if e >= p0 + k else [p0 + d] if e <= p0 else []
    if not los or alts or r:
        if p0 not in los and e < p0 + k and slides(e, p0 + k - e):
            los.append(p0)  # the junction moved past the PAM
        if p0 + d not in los and e > p0 and slides(p0, e - p0):
            los.append(p0 + d)  # before it
    if los:
        return "ref", tuple(los)
    return "ref", ((p0,) if pam_right else (p0 + d,))


# ---------------------------------------------------------------------------
# A second indel inside a target, read on the sequence (2.8.1)
# ---------------------------------------------------------------------------
# The fake-indel window holds one indel, so a carrier whose haplotype also holds a
# second indel inside the target may not have the window's sequence there. Whether it
# has is a question of bases, not of positions: in a tandem repeat a unit gained or lost
# at the target's edge leaves the target's bases as they are (the repeat puts the same
# bases back), and that carrier does have the row's sequence. local_haplotype writes out
# a haplotype's bases around the target, each with the place it holds, and
# holds_target_bases asks whether the window's target is on it at the target's own
# place: from its first base on, or up to its last base.


def window_place(offset, indel_pos1, indel_start, ref, alt):
    """The place (local_haplotype) of the base at ``offset`` of a fake-indel window
    whose real 0-based start is ``indel_start``: a reference base before the indel, one
    of the ALT bases (its anchor included: the window writes the ALT as it is), or a
    reference base past the REF."""
    alt_at = indel_pos1 - 1 - indel_start
    if offset < alt_at:
        return (indel_start + offset, 0)
    if offset < alt_at + len(alt):
        return (indel_pos1 - 2, 2, offset - alt_at)
    return (indel_start + offset - len(alt) + len(ref), 0)


def local_haplotype(lo, hi, window, others=()):
    """The bases of a haplotype over the reference positions lo..hi (0-based,
    inclusive), as [(place, source)] in the order the haplotype holds them.

    ``window`` (VCF POS, REF, ALT) is the indel of the fake-indel window, written as the
    window writes it: its REF out, its ALT in, the anchor base included; None for a
    haplotype read around a SNP row's target, which holds no window indel. ``others`` are
    the haplotype's other indels; the bases one shares with its REF at either end
    (kept_bases, shared_tail) stay reference bases. A source is a reference position
    (int), whose base is the haplotype's own allele there, or a base an ALT puts in
    (str). A place sorts as the haplotype runs: (p, 0) for reference position p, and
    (p, rank, k) for the k-th base an ALT puts in after reference position p, another
    indel's (rank 1) before the window's (rank 2), which follows the bases an indel
    inserted before its anchor. None where two of the variants overlap, or insert at
    one point: no one haplotype holds both (a VCF writes such a haplotype as one
    record). Two variants overlap when their REFs share a base, whatever bases each
    keeps: two records at one POS always do (an insertion and a deletion there are two
    alleles of one site, though the bases each changes do not meet). The bases do not
    depend on the order of ``others``."""
    records = ([window] if window is not None else []) + list(others)
    refs = sorted((p1 - 1, p1 - 1 + len(ref)) for p1, ref, _alt in records)
    if any(b[0] < a[1] for a, b in zip(refs, refs[1:])):
        return None  # two REFs share a base
    edits = []
    if window is not None:
        p1, ref, alt = window
        s = p1 - 1
        edits.append((s, s + len(ref), alt, s - 1, 2))
    for p1, ref, alt in others:
        s = p1 - 1
        c = kept_bases(ref, alt)
        t = shared_tail(ref, alt)
        edits.append((s + c, s + len(ref) - t, alt[c:len(alt) - t], s + c - 1, 1))
    # with no two REFs sharing a base, only an insertion can start where another record's
    # changed bases start, and it comes first: its bases go in before them
    edits.sort(key=lambda e: (e[0], e[1], e[4]))
    out = []
    at = lo
    last = None
    for start, end, put, after, rank in edits:
        if last is not None and start == last[0] == last[1] == end and rank == last[2]:
            return None  # two insertions at one point: their bases have no one order
        out.extend(((p, 0), p) for p in range(at, min(start, hi + 1)))
        out.extend(((after, rank, k), b.upper()) for k, b in enumerate(put))
        at = max(at, end)
        last = (start, end, rank)
    out.extend(((p, 0), p) for p in range(at, hi + 1))
    return out


def _same_base(x, y, value):
    """Whether sources x and y (local_haplotype) are the same base on a haplotype: the
    same reference position is, whatever its allele; otherwise their bases are
    compared, and a base ``value`` cannot tell (None) is not the same."""
    if x == y:
        return True
    bx = x if isinstance(x, str) else value(x)
    by = y if isinstance(y, str) else value(y)
    return bx is not None and bx == by


def target_pairs(window_hap, hap, first, last, through=()):
    """The ways ``hap`` can hold the bases ``window_hap`` holds at places first..last
    (the target's first and last base on the window; local_haplotype), at the target's
    own place: one from the place of its first base on, one up to the place of its last
    base. Each is the list of (window source, haplotype source) pairs that are not the
    same reference position, or None where it cannot hold them there: an ALT base
    against another, a stretch that runs past the haplotype written out, or one that
    leaves out what the target holds of the window's indel (its ALT bases on the window,
    and the base past them: for a deletion, the junction). The haplotype must have the
    bases through that indel, as the row names it, not only in a stretch that another of
    its variants has lengthened (a repeat unit gained upstream can push the indel out of
    a stretch whose bases are those of the window). ``through``: more places the stretch
    must hold, the bases of the SNPs a SNP row names (snp_target_held). Always two, in
    that order: a target the window written out does not hold has neither."""
    keys = [k for k, _s in window_hap]
    try:
        i, j = keys.index(first), keys.index(last)
    except ValueError:
        return [None, None]
    want = [s for _k, s in window_hap[i:j + 1]]
    n = len(want)
    lit = [x for x in range(i, j + 1) if len(keys[x]) == 3 and keys[x][1] == 2]
    core = {keys[x] for x in lit} | set(through)
    if lit and lit[-1] < j:
        core.add(keys[lit[-1] + 1])  # the base past the window's ALT: a deletion's junction
    hkeys = [k for k, _s in hap]
    ways = []
    for a in (bisect.bisect_left(hkeys, first), bisect.bisect_right(hkeys, last) - n):
        if a < 0 or a + n > len(hap) or not core <= set(hkeys[a:a + n]):
            ways.append(None)  # past the bases written out, or not through the indel
            continue
        pairs = []
        for x, (_k, y) in zip(want, hap[a:a + n]):
            if x == y:
                continue
            if isinstance(x, str) and isinstance(y, str):
                pairs = None  # two ALT bases that differ
                break
            pairs.append((x, y))
        ways.append(pairs)
    return ways


def holds_target_bases(ways, value):
    """True iff one of ``ways`` (target_pairs) holds: every pair is the same base on the
    haplotype whose alleles ``value`` gives (source -> base, None where not known)."""
    return any(pairs is not None and all(_same_base(x, y, value) for x, y in pairs)
               for pairs in ways)


# ---------------------------------------------------------------------------
# A SNP row's carrier whose haplotype holds an indel inside the target (2.10.0)
# ---------------------------------------------------------------------------
# A SNP row is a target on the reference with SNP alleles: the enriched genome holds
# no indel, so its carriers are read from SNP genotypes. A carrier whose haplotype
# holding the row's SNPs also holds an indel that changes the target's bases does not
# have the row's sequence there. As for an indel row's second indel (2.8.1), that is
# read on the haplotype's bases: in a tandem repeat a unit gained or lost at the
# target's edge leaves them as they are.


def snp_target_span(start, target_len):
    """(lo, hi): the reference bases (0-based, inclusive) a SNP row's target covers: its
    start and as many bases as the target holds, the PAM and any DNA-bulge base
    included (an RNA bulge is a gap of the guide, no base of the target)."""
    return start, start + target_len - 1


def snp_target_frame(lo, hi, indels):
    """(first, last): the reference positions (0-based, inclusive) snp_target_held
    writes a haplotype over: the target's span and, on either side, its length and the
    bases the haplotype's indels ``indels`` [(POS, REF, ALT)] add or delete."""
    margin = hi - lo + 1 + sum(abs(len(r) - len(a)) for _p, r, a in indels)
    return lo - margin, hi + margin


def snp_target_held(lo, hi, snps, held, value, reference):
    """True iff a haplotype holding a SNP row's SNPs and some indels has the row's
    target at its place: the bases of the reference lo..hi (0-based, inclusive) with the
    haplotype's own alleles there, read from the target's first base on or up to its
    last base, the stretch holding the base of every SNP the row names (``snps``, their
    0-based positions), not one another variant pushed out of it.

    ``held`` gives, for each way (from the first base on, up to the last base), the
    haplotype's indels [(POS, REF, ALT)] that way reads, or None where it cannot be read
    (an indel there that cannot be put on a slot). ``value`` (source -> base, None where
    not known) gives the haplotype's alleles, the row's at a position of lo..hi;
    ``reference`` is the chromosome's sequence. A way that reads two records of one base
    (two indels whose REFs share a base, or a SNP alt of the haplotype on a base of an
    indel's REF) cannot be read."""
    through = [(p, 0) for p in snps]
    for side in (0, 1):
        indels = held[side]
        if indels is None:
            continue
        if any(not 0 <= q < len(reference) or value(q) != reference[q]
               for p, r, _a in indels for q in range(p - 1, p - 1 + len(r))):
            continue  # a SNP alt on a base of the indel's REF: two records of one base
        first, last = snp_target_frame(lo, hi, indels)
        hap = local_haplotype(first, last, None, indels)
        if hap is None:
            continue  # two of its indels overlap
        ways = target_pairs(local_haplotype(first, last, None), hap, (lo, 0), (hi, 0),
                            through)
        if holds_target_bases([ways[side]], value):
            return True
    return False


class IndelCarriers(object):
    """The carriers of one indel, parsed once: {sample: (genotype tokens, kind,
    indel-carrying slots)} plus the slot total. ``phase_known`` is False when they
    come without phase (the indel log, no phased indel genotype store): no slot can
    then be tied to a SNP."""

    __slots__ = ("parsed", "total_slots", "phase_known")

    def __init__(self, indel_gts, phase_known=True):
        self.parsed = {}
        self.total_slots = 0
        self.phase_known = phase_known
        for sample, gt in indel_gts.items():
            toks, kind = _gt_tokens(gt)
            islots = [h for h, t in enumerate(toks) if t == "1"]
            if islots:
                self.parsed[sample] = (toks, kind, islots)
                self.total_slots += len(islots)


# A genotype as "this alt or not": the alt's allele index becomes "1", any other allele
# "0", a missing allele stays "."
_for_alt = site_alleles.for_alt


def biallelic_site(rec):
    """{alt: {sample: genotype}} at one position, each genotype read against its own alt.

    The resolution reads "1" as "this slot holds this alt". That is what a split VCF
    record holds, but a record with several alts on one line (C>A,T) keeps its raw
    genotype under every alt its sample carries: S:2|1 is listed under A and under T.
    The allele index of an alt is the token all its carriers share, and when that
    leaves more than one (a single carrier, 2|1), the one no other alt of the position
    needs; the record's own alt then reads "1", any other allele "0": 2|1 is 0|1 under
    A, 1|0 under T. Where the genotypes cannot say which index is the alt (every carrier
    of two alts is 1|2) no slot can be tied to it: "1/0", one copy on an unknown slot,
    so the carrier is resolved PUTATIVE on the alleles it does carry. A position with
    no foreign index (split or biallelic records) comes back as it is. The rule is
    every reader's (site_alleles.per_alt)."""
    return site_alleles.per_alt(rec)


def resolve_carriers(columns, indel_gts, phase_known=True, others=()):
    """Split the indel's carriers into exactly resolved slots and PUTATIVE samples.

    Args:
      columns: the target's code columns (SnpColumn), in column order.
      indel_gts: {sample: normalized indel genotype} ("1|0", "0|1", "1|1", "0/1",
        haploid "1"), or an IndelCarriers built from it.
      phase_known: False when the indel's carriers come without phase (the indel
        log, no phased indel genotype store): no slot can then be tied to a SNP.
        (An IndelCarriers carries its own.)
      others: one {sample: normalized genotype} per OTHER variant that leaves a slot
        without the window's sequence there: a second indel that changes the target's
        bases (the caller reads that on the slot's sequence: local_haplotype), or a
        complex allele split into overlapping records. The fake-indel window holds one
        indel, so a slot marked "1" there is no carrier of any row of this target. A
        sample whose second variant cannot be put on a slot (unphased) is PUTATIVE.

    Returns ``(exact, putative)``:
      exact: {bases: [samples, n_slots]} -- ``bases`` is a tuple with the base each
        indel-carrying slot has at every code column (target orientation).
      putative: {candidates: samples} -- ``candidates`` is a tuple holding, per code
        column, the bases the sample may carry there (any slot): its own alleles
        where it is genotyped, the code's alleles where it is not.

    A slot is exact when the indel and every SNP the sample carries at the columns
    are phased with the same ploidy (or all haploid) and each allele the slot holds
    is known. A sample that carries no SNP alt at any column has the reference base
    on every indel slot whatever its phase. Any other sample is PUTATIVE. Only the
    samples that carry a SNP alt at the columns are looked at one by one, so a
    common indel costs what its SNP co-carriers cost.
    """
    if not isinstance(indel_gts, IndelCarriers):
        indel_gts = IndelCarriers(indel_gts, phase_known)
    phase_known = indel_gts.phase_known
    parsed = indel_gts.parsed
    ref_bases = tuple(c.ref for c in columns)
    any_unknown = any(not c.known for c in columns)
    exact = {}
    putative = {}
    def _pairs(carriers):
        # (sample, genotype) of the carriers that also carry the indel
        if len(carriers) <= len(parsed):
            return [(s, g) for s, g in carriers.items() if s in parsed]
        return [(s, carriers[s]) for s in parsed if s in carriers]

    seen_by = {}  # sample -> [(column index, base shown, genotype tokens, kind)]
    for k, c in enumerate(columns):
        for base, (_alt, carriers) in c.alts.items():
            for s, g in _pairs(carriers):
                toks, kind = _gt_tokens(g)
                seen_by.setdefault(s, []).append((k, base, toks, kind))
    blocked = {}  # sample -> slots with a second indel in the span (None: unknown)
    for jgts in others:
        for s, g in _pairs(jgts):
            itoks, ikind, _islots = parsed[s]
            jtoks, jkind = _gt_tokens(g)
            if (phase_known and ikind != "unphased" and jkind == ikind
                    and len(jtoks) == len(itoks)):
                if blocked.get(s, ()) is not None:
                    blocked[s] = set(blocked.get(s, ())) | {
                        h for h, t in enumerate(jtoks) if t == "1"}
            else:
                blocked[s] = None
    if any_unknown:
        todo = sorted(parsed)
    else:
        todo = sorted(set(seen_by) | set(blocked))
        rest = len(parsed) - len(todo)
        if rest:
            ent = exact.setdefault(ref_bases, [set(), 0])
            ent[0].update(s for s in parsed if s not in seen_by and s not in blocked)
            ent[1] += indel_gts.total_slots - sum(len(parsed[s][2]) for s in todo)
    for sample in todo:
        itoks, ikind, islots = parsed[sample]
        seen = seen_by.get(sample, ())
        unsure = sample in blocked and blocked[sample] is None
        if blocked.get(sample):
            islots = [h for h in islots if h not in blocked[sample]]
            if not islots:
                continue  # every indel slot also holds a variant that breaks the window
        if not seen and not any_unknown and not unsure:
            ent = exact.setdefault(ref_bases, [set(), 0])
            ent[0].add(sample)
            ent[1] += len(islots)
            continue
        slots = None
        if phase_known and not any_unknown and not unsure and ikind != "unphased":
            slots = []
            for h in islots:
                bases = list(ref_bases)
                for k, base, toks, kind in seen:
                    if kind != ikind or len(toks) != len(itoks):
                        slots = None
                        break
                    if toks[h] == "1":
                        if bases[k] != ref_bases[k]:
                            slots = None  # two alts on one slot: not a haplotype
                            break
                        bases[k] = base
                    elif toks[h] != "0":
                        slots = None  # a missing allele on the indel's slot
                        break
                if slots is None:
                    break
                slots.append(tuple(bases))
        if slots is not None:
            for bases in slots:
                ent = exact.setdefault(bases, [set(), 0])
                ent[0].add(sample)
                ent[1] += 1
            continue
        cands = []
        for k, c in enumerate(columns):
            if not c.known:
                cands.append(c.code | {c.ref})
                continue
            carried = set()
            copies = 0
            ploidy = len(itoks)
            unknown = False
            for kk, base, toks, _kind in seen:
                if kk != k:
                    continue
                carried.add(base)
                copies += sum(1 for t in toks if t == "1")
                ploidy = max(ploidy, len(toks))
                unknown = unknown or any(t not in ("0", "1") for t in toks)
            if not carried or unknown or copies < ploidy:
                carried.add(c.ref)
            cands.append(frozenset(carried))
        putative.setdefault(tuple(cands), set()).add(sample)
    return exact, putative


def worst_case_bases(columns, candidates, guide_seq, motif_at, score_fn=None):
    """The strongest target a PUTATIVE carrier's own alleles allow: the bases to put
    at each code column, or None when no allele it may carry satisfies the PAM.

    At a spacer column the guide's own base wins when it is among the candidates;
    otherwise, and at the PAM columns (jointly, among the bases that satisfy the
    motif), the bases that maximize ``score_fn`` (the CFD of the whole alignment,
    which factors per column) win. Ties, and every choice without a ``score_fn``,
    keep the reference base first and then A, C, G, T, so no SNP is named for a
    column it does not change the score of.
    """
    def order(c, cand):
        return sorted(cand, key=lambda b: (b != c.ref, _BASE_ORDER.find(b)))

    picks = []
    free = []    # (column index, candidates) a score decides
    pam_ks = []  # (column index, motif-satisfying candidates)
    for k, c in enumerate(columns):
        cand = order(c, candidates[k])
        motif = motif_at.get(c.col)
        if motif is not None:
            ok = [b for b in cand if b in _CODE.get(motif, "")]
            if not ok:
                return None
            picks.append(ok[0])
            if len(ok) > 1:
                pam_ks.append((k, ok))
            continue
        g = guide_seq[c.col].upper()
        if g in cand:
            picks.append(g)
        else:
            picks.append(cand[0])
            if len(cand) > 1:
                free.append((k, cand))
    if score_fn is None or (not free and not pam_ks):
        return tuple(picks)
    # CFD is a product of one factor per spacer column times a joint factor of the
    # PAM's last two bases, so the PAM columns are chosen together and each spacer
    # column alone. A factor of 0 elsewhere (a PAM CFD scores 0, as NRR's AA) hides
    # every choice, so the passes repeat until neither changes anything.
    for _ in range(len(free) + 2):
        changed = False
        if pam_ks:
            best = None
            for combo in itertools.product(*[ok for _k, ok in pam_ks]):
                trial = list(picks)
                for (k, _ok), b in zip(pam_ks, combo):
                    trial[k] = b
                s = score_fn(tuple(trial))
                if best is None or s > best[0]:
                    best = (s, combo)
            for (k, _ok), b in zip(pam_ks, best[1]):
                if picks[k] != b:
                    picks[k] = b
                    changed = True
        for k, cand in free:
            best = None
            for b in cand:
                trial = list(picks)
                trial[k] = b
                s = score_fn(tuple(trial))
                if best is None or s > best[0]:
                    best = (s, b)
            if best[1] != picks[k]:
                picks[k] = best[1]
                changed = True
        if not changed:
            break
    return tuple(picks)


def realize(target_seq, guide_seq, columns, bases):
    """The concrete aligned target with ``bases`` at the code columns: upper case
    where the base matches the guide (or the guide has 'N' / a bulge gap there),
    lower case where it is a mismatch, as the engine marks its own mismatches."""
    out = list(target_seq)
    for c, b in zip(columns, bases):
        g = guide_seq[c.col].upper()
        out[c.col] = b if (g in ("N", "-") or g == b) else b.lower()
    return "".join(out)


def count_mismatches(guide_seq, target_seq, spacer_cols):
    """Mismatches of a concrete target over its protospacer columns: a column where
    neither side is a bulge gap and the bases differ."""
    n = 0
    for j in spacer_cols:
        g, t = guide_seq[j], target_seq[j]
        if g != "-" and t != "-" and g.upper() != t.upper():
            n += 1
    return n


def pam_satisfied(target_seq, motif_at):
    """True iff every PAM column of a concrete target holds a base the motif allows
    there (a genomic N, or a gap, satisfies nothing)."""
    for j, m in motif_at.items():
        b = target_seq[j].upper()
        if b not in _BASE_ORDER or b not in _CODE.get(m, ""):
            return False
    return True


def enabling_columns(guide_seq, columns, bases, motif_at):
    """Indexes of the code columns whose alt allele the target REQUIRES (the SNPs a
    SNP+indel co-occurrence is made of): the base is an alt, and it is the guide's
    base at a spacer column, or it satisfies the motif where the reference does not
    at a PAM column."""
    out = []
    for k, (c, b) in enumerate(zip(columns, bases)):
        if b == c.ref:
            continue
        motif = motif_at.get(c.col)
        if motif is not None:
            allowed = _CODE.get(motif, "")
            if b in allowed and c.ref not in allowed:
                out.append(k)
        elif guide_seq[c.col].upper() == b:
            out.append(k)
    return out


def parse_dict_value(value):
    """Parse one per-sample SNP dictionary value (``my_dict_<chrom>.json``): entries
    ``<sample:gt,...>;<ref>,<alt>;<rsID>;<AF>`` joined by '$', one per alt. Returns
    ``[(ref, alt, rsid, af, {sample: gt})]`` for the single-base ACGT alts (the SNPs
    the overlay codes), one per REF>ALT, each genotype read against its own alt's
    allele index (site_alleles.read_records) among every alt of one character the value
    lists, a '*' too (the Tier-1 store's compile, tier0_compile.iter_dict_records, reads
    them so: a carrier of the '*' alone tells its index); malformed entries are
    skipped."""
    out = []
    for entry in str(value).split("$"):
        parts = entry.split(";")
        if len(parts) < 3:
            continue
        refalt = parts[1].strip().split(",")
        if len(refalt) < 2:
            continue
        ref, alt = refalt[0].strip().upper(), refalt[1].strip().upper()
        if len(ref) != 1 or len(alt) != 1:
            continue
        carriers = {}
        for tok in parts[0].split(","):
            sid, sep, gt = tok.strip().partition(":")
            if sid and sep:
                gt = gt.strip()
                carriers[sid] = (site_alleles.merge_gt(carriers[sid], gt)
                                 if sid in carriers else gt)
        rsid = parts[2].strip() or "."
        af = parts[3].strip() if len(parts) > 3 else ""
        out.append((ref, alt, rsid, af, carriers))
    # each genotype read against its own alt's allele index, two records of one REF>ALT
    # merged (a slot carries the alt if either puts it there): site_alleles.read_records;
    # then the bases
    return [r for r in site_alleles.read_records(out)
            if r[0] in _BASE_ORDER and r[1] in _BASE_ORDER]


def load_dict_snps(dict_path, chrom, positions):
    """The per-sample SNP dictionary entries at the 1-based ``positions`` of one
    chromosome: {pos1: parse_dict_value(...)}. A classic dictionary install (no
    Tier-1 genotype store, e.g. an index built without --samplesID) holds the
    genotypes the per-haplotype resolution needs only there. Prefers the gzipped
    dictionary when present; streams it with ijson when installed, so memory grows
    with the positions asked for, not the dictionary (else one json.load)."""
    import gzip
    import json
    import os

    if os.path.exists(dict_path + ".gz"):
        dict_path = dict_path + ".gz"
    if not os.path.exists(dict_path) or not positions:
        return {}
    keys = {"%s,%d" % (chrom, p): p for p in positions}
    out = {}

    def _keep(key, value):
        p = keys.get(key)
        if p is not None and value:
            recs = parse_dict_value(value)
            if recs:
                out[p] = recs

    is_gz = dict_path.endswith(".gz")
    try:
        import ijson
    except ImportError:
        ijson = None
    if ijson is not None:
        with (gzip.open(dict_path, "rb") if is_gz else open(dict_path, "rb")) as fh:
            for key, value in ijson.kvitems(fh, ""):
                _keep(key, value)
    else:
        with (gzip.open(dict_path, "rt") if is_gz else open(dict_path)) as fh:
            for key, value in json.load(fh).items():
                _keep(key, value)
    return out


def resolve_target(guide_seq, target_seq, columns, indel_gts, motif_at, spacer_cols,
                   max_mm, max_total, bulge_size, phase_known=True, score_fn=None,
                   max_dna_bulges=None, max_rna_bulges=None, others=()):
    """Resolve one overlaid indel target into the rows its carriers really have.

    Args:
      guide_seq, target_seq: the engine's aligned crRNA and DNA (codes included).
      columns: the target's code columns (SnpColumn), in column order.
      indel_gts, phase_known: as for ``resolve_carriers``. ``indel_gts=None`` means
        there is nothing to resolve the carriers with (no SNP genotypes at all, a
        sites-only panel): one PUTATIVE row then holds the worst case over the
        code's alleles at every column, its ``samples`` None (the caller keeps the
        carriers the indel log lists, as the 2.5.2 fallback did).
      motif_at: {column: motif character} for the PAM columns (``pam_columns``).
      spacer_cols: the protospacer columns (mismatches are counted there).
      max_mm, max_total, bulge_size: the search's --mm, its total-edits cap (< 0:
        none) and the target's bulge count, which resolution never changes.
      score_fn: bases tuple -> score of the realized target (the CFD), used to pick
        a PUTATIVE carrier's worst case; None picks by fewest mismatches.
      max_dna_bulges, max_rna_bulges: the search's bulge limits (the engine kept
        the target within them; checked again for completeness).
      others: carriers of the other indels that change the target's span, as for
        ``resolve_carriers``.

    Returns a list of dicts, one per surviving sequence, in a deterministic order:
      seq      the concrete aligned target (case-marked)
      mm       its mismatches
      bases    the base at each code column (target orientation)
      samples  the sorted carrier samples whose haplotype it is
      slots    exact carrier slots (haplotype copies) + one per PUTATIVE sample
      phase    CONFIRMED when every carrier is resolved exactly, else PUTATIVE
      enabling indexes of the columns whose alt the target requires
    A sequence no carrier has is never returned; nor one whose PAM fails the motif,
    or whose recounted edits exceed the budget.
    """
    groups = {}  # bases -> [samples, slots, putative?]
    if indel_gts is None:
        cands = tuple(c.code | {c.ref} for c in columns)
        bases = worst_case_bases(columns, cands, guide_seq, motif_at, score_fn)
        if bases is not None:
            groups[bases] = [None, 0, True]
    else:
        exact, putative = resolve_carriers(columns, indel_gts, phase_known, others)
        for bases, (samples, slots) in exact.items():
            groups[bases] = [set(samples), slots, False]
        picked = {}
        for cands, samples in sorted(putative.items(), key=lambda kv: sorted(kv[1])):
            if cands not in picked:
                picked[cands] = worst_case_bases(columns, cands, guide_seq, motif_at,
                                                 score_fn)
            bases = picked[cands]
            if bases is None:
                continue
            ent = groups.setdefault(bases, [set(), 0, False])
            ent[0] |= samples
            ent[1] += len(samples)
            ent[2] = True
    dna_bulges = guide_seq.count("-")
    rna_bulges = target_seq.count("-")
    if (max_dna_bulges is not None and dna_bulges > max_dna_bulges) or (
            max_rna_bulges is not None and rna_bulges > max_rna_bulges):
        return []
    rows = []
    for bases in sorted(groups):
        samples, slots, is_putative = groups[bases]
        seq = realize(target_seq, guide_seq, columns, bases)
        if not pam_satisfied(seq, motif_at):
            continue  # this haplotype has no PAM here
        mm = count_mismatches(guide_seq, seq, spacer_cols)
        if mm > max_mm or (max_total is not None and max_total >= 0
                           and mm + bulge_size > max_total):
            continue  # this haplotype is not a target within the search budget
        rows.append({
            "seq": seq, "mm": mm, "bases": bases,
            "samples": None if samples is None else sorted(samples),
            "slots": slots, "phase": PUTATIVE if is_putative else CONFIRMED,
            "enabling": enabling_columns(guide_seq, columns, bases, motif_at),
        })
    return rows
