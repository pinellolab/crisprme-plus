#!/usr/bin/env python3
"""The exact haplotype pre-check of a per-sample search, and the population-level K-skip
(2.10.0).

A raw target of the variant search is a window of the variant genome, its variant positions
written as IUPAC letters, that the guide matches within the budget when each letter may be
any base it stands for. The post-analysis then asks which sequences real people carry there
(new_simple_analysis.py for the SNP windows, analisi_indels_NNN.py for the windows of the
fake-indel genome). In a dense region (the immunoglobulin loci) almost none of these raw
targets is a sequence anybody carries, and expanding each one sample by sample is most of the
cost of a search.

THE PRE-CHECK (per sample). Every row the per-sample post-analysis can write for a raw target
puts, at each variant column, the reference base or an ALT that one of its carriers carries:
a phased slot's alleles, a PUTATIVE sample's own alleles (its subsets, or the unions past
CRISPRME_IUPAC_CAP), the classic dictionaries' lattice of carriers shared by every allele of a
combination. The slot filter of an indel (it takes a carrier off, or keeps it PUTATIVE beside
an unphased indel, on the row its haplotype already gives), --max-degeneracy and the #139
column admission only take alleles or carriers away. So the target can give a variant row only
if one sample, free to
put at each column the reference base or any ALT it carries there (its "union", with no phase),
reaches the budget and the PAM. Where no sample can, the target gives no variant row: it is
skipped, and only its reference row is written (the post-analysis writes that one whatever the
genotypes). The test is a bound, never a guess: it may keep a target that gives no row (one
sample's union is wider than any of its haplotypes), never skip one that gives a row.

The finalizer counts a row's mismatches over the protospacer columns of its aligned target
(bulge gaps included, the bulges then taken off) and needs every PAM column to hold a base the
PAM admits. Mismatches and the PAM are per column, so the test is per column too:

* a protospacer column where the reference base mismatches the guide is FIXED by an ALT that
  is the guide's base there; ``need`` = the reference alignment's mismatches minus the budget
  is how many columns a row must fix;
* a PAM column whose reference base the PAM does not admit needs an ALT it admits.

A sample can produce the target iff it carries an admitted ALT at every such PAM column and a
fixing ALT at ``need`` columns or more. Samples are bits of Python integers, so "some sample
is in ``need`` of these sets" is a few integer operations a column (``at_least``).

THE K-SKIP (population-level, --max-degeneracy K). Every row a capped population-level search
writes rests on at most K ambiguous columns (window_density.degeneracy), and every row within
the budget and the PAM rests on at least D_min = max(0, need) + the PAM columns that need an
ALT, whatever ALT it takes. A raw target whose D_min exceeds K therefore gives no variant row
and is skipped the same way, before its genotypes are read.

Stdlib only and free of the post-analysis scripts' module globals, so it is unit-tested alone
(test_haplotype_precheck.py, against a brute-force enumeration of every haplotype).
"""

# The values that turn off a switch of the 2.10.0 speed steps (CRISPRME_HAPLOTYPE_PRECHECK,
# CRISPRME_POSTPROC_PIECES), in any case, spaces around them aside. Any other value, or none,
# leaves it on.
SWITCH_OFF = ("0", "false", "no", "off")


def switch_on(environ, name):
    """True unless the variable ``name`` of ``environ`` says off (SWITCH_OFF)."""
    return str(environ.get(name, "1")).strip().lower() not in SWITCH_OFF


# The bases each IUPAC letter of a PAM file admits (a PAM's N is any base).
PAM_ADMITS = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(base):
    """The complementary base (a letter that is no base stays as it is)."""
    return base.translate(_COMPLEMENT)


class SampleBits(object):
    """Sample ids as bits of Python integers: one bit per sample, in first-seen order."""

    __slots__ = ("_bit", "names")

    def __init__(self):
        self._bit = {}
        self.names = []

    def bit(self, sample):
        b = self._bit.get(sample)
        if b is None:
            b = self._bit[sample] = len(self.names)
            self.names.append(sample)
        return b

    def mask(self, samples):
        """The integer whose bits are ``samples``."""
        bits = [self.bit(s) for s in samples]
        if not bits:
            return 0
        buf = bytearray(max(bits) // 8 + 1)
        for b in bits:
            buf[b >> 3] |= 1 << (b & 7)
        return int.from_bytes(bytes(buf), "little")

    def samples(self, mask):
        """The sample ids of ``mask``'s bits."""
        bits = bin(mask)[:1:-1]  # bit 0 first
        return [self.names[i] for i, ch in enumerate(bits) if ch == "1"]


def at_least(masks, k):
    """The samples (bits) present in at least ``k`` of ``masks`` (k >= 1)."""
    if k <= 0:
        raise ValueError("k must be at least 1")
    acc = [0] * (k + 1)  # acc[j]: in at least j of the masks seen so far
    for m in masks:
        if not m:
            continue
        for j in range(k, 1, -1):
            acc[j] |= acc[j - 1] & m
        acc[1] |= m
    return acc[k]


def aligned_places(real_target, n_bases, revert):
    """For each base of a raw target's pre-revert sequence (the reference strand's window,
    ``n_bases`` long), its column in the aligned target the finalizer builds: reverse
    complemented on the minus strand, with ``real_target``'s bulge gaps ('-') put back. None
    when the window and the alignment do not have the same number of bases (a window cut at a
    chromosome end), which no test here reads."""
    places = [j for j, ch in enumerate(real_target) if ch != "-"]
    if len(places) != n_bases:
        return None
    if revert:
        return places[::-1]
    return places


def snp_columns(real_target, ref_prerevert, revert, guide_no_pam, pam, pos_beg, pos_end,
                pam_begin, pam_end, variant_columns, admits):
    """What each variant column of a SNP raw target can do for its rows, read with the rules of
    new_simple_analysis's finalizer.

    Args:
      real_target: the raw target's aligned DNA (split[2]: its bulge gaps in place).
      ref_prerevert: the reference window on the forward strand (upper case), as long as the
        raw target's bases.
      revert: True on the minus strand.
      guide_no_pam: the aligned guide without its PAM (split[1][pos_beg:pos_end]).
      pam, pos_beg, pos_end, pam_begin, pam_end: the PAM letters and the four slice bounds of
        the finalizer (the protospacer is [pos_beg:pos_end] of the aligned target, the PAM
        [pam_begin:pam_end]).
      variant_columns: [(pos_c, alts)]: each variant column's pre-revert index and the forward
        ALTs that may stand there (the bases the post-analysis may put there).
      admits: callable(pam_letter, base) -> True when the PAM letter admits the base (the
        finalizer's iupac_code_set test).

    Returns None when the test cannot be read (the window is cut short), else a tuple
    ``(base_mm, fixes, pam_needs, blocked)``:
      base_mm: the reference alignment's mismatches as the finalizer counts them;
      fixes: [(pos_c, [alt, ...])]: each protospacer column whose reference base mismatches
        the guide, with the ALTs that are the guide's base there;
      pam_needs: [(pos_c, [alt, ...])]: each PAM column whose reference base the PAM does not
        admit, with the ALTs it admits;
      blocked: True when a PAM column with no variant holds a base the PAM does not admit (no
        row of the window can then have its PAM).
    """
    places = aligned_places(real_target, len(ref_prerevert), revert)
    if places is None:
        return None
    ref_fwd = ref_prerevert
    n = len(real_target)
    aligned = ["-"] * n
    for i, j in enumerate(places):
        b = ref_fwd[i]
        aligned[j] = complement(b) if revert else b
    ps = list(range(n))[pos_beg:pos_end]
    pm = list(range(n))[pam_begin:pam_end]
    if len(ps) != len(guide_no_pam) or len(pm) != len(pam):
        return None
    ps_at = {j: t for t, j in enumerate(ps)}
    pm_at = {j: t for t, j in enumerate(pm)}
    base_mm = sum(1 for t, j in enumerate(ps) if aligned[j].upper() != guide_no_pam[t])
    variant_place = {}
    fixes, pam_needs = [], []
    for pos_c, alts in variant_columns:
        if not 0 <= pos_c < len(places):
            continue
        j = places[pos_c]
        variant_place[j] = pos_c
        shown = [(a, complement(a) if revert else a) for a in alts]
        ref_here = aligned[j].upper()
        if j in ps_at:
            g = guide_no_pam[ps_at[j]]
            if ref_here != g:
                fixes.append((pos_c, [a for a, s in shown if s.upper() == g]))
        elif j in pm_at:
            letter = pam[pm_at[j]]
            if not admits(letter, ref_here):
                pam_needs.append((pos_c, [a for a, s in shown if admits(letter, s.upper())]))
    blocked = any(j not in variant_place and not admits(pam[t], aligned[j].upper())
                  for t, j in enumerate(pm))
    return base_mm, fixes, pam_needs, blocked


def sample_can_produce(need, fix_masks, pam_masks, window_mask):
    """True when one sample of ``window_mask`` (the samples carrying an ALT of the window)
    carries an ALT the PAM admits at every PAM column of ``pam_masks`` and an ALT that fixes a
    mismatch at ``need`` or more of the columns of ``fix_masks`` (masks of samples, one a
    column)."""
    cand = window_mask
    for m in pam_masks:
        cand &= m
        if not cand:
            return False
    if not cand:
        return False
    if need <= 0:
        return True
    usable = [m & cand for m in fix_masks if m & cand]
    if len(usable) < need:
        return False
    return at_least(usable, need) != 0


def degeneracy_floor(need, fixes, pam_needs, blocked):
    """D_min of a raw target: the fewest ambiguous columns (window_density.degeneracy) any of
    its rows within the budget and the PAM rests on, whatever ALTs it takes; None when no row
    can be within them (a PAM no ALT satisfies, or fewer fixable columns than ``need``)."""
    if blocked or any(not alts for _pc, alts in pam_needs):
        return None
    fixable = sum(1 for _pc, alts in fixes if alts)
    if need > fixable:
        return None
    return max(0, need) + len(pam_needs)


def indel_target_possible(guide_seq, target_seq, motif_at, spacer_cols, max_mm, max_total,
                          bulge_size, columns, carriers, everyone, max_dna_bulges=None,
                          max_rna_bulges=None, code_bases=None):
    """The same pre-check for an overlaid target of the fake-indel genome
    (indel_snp_cis.resolve_target's gates): False when no carrier of the window's indel can
    have one of its rows, True otherwise (a bound: it may keep a target that gives no row).

    Every row resolve_target writes holds, at each code column, the reference base or a base
    one of its carriers carries there (an exact slot's alleles; a PUTATIVE carrier's own
    alleles, the reference among them; the code's alleles at a column no genotype covers). A
    second indel that breaks a slot takes that slot's carrier away, or, where it cannot be put
    on a slot (unphased, or of another ploidy), makes the carrier PUTATIVE
    (indel_snp_cis.resolve_carriers): that can add a row, one the carrier's slots do not hold
    (a base of one slot at one column, of the other at another), but a PUTATIVE row is still
    a choice of the carrier's own alleles at each code column. So a carrier that may put at
    each code column the reference base or any base it carries there (any base of the code
    where nobody's genotype is known) is the widest row a carrier can have, with or without
    second indels.

    Args:
      guide_seq, target_seq: the engine's aligned crRNA and DNA (codes included).
      motif_at: {column: PAM letter} (indel_snp_cis.pam_columns).
      spacer_cols: the protospacer columns (mismatches are counted there).
      max_mm, max_total, bulge_size: the search's --mm, its total-edits cap (< 0 or None:
        none) and the target's bulge count.
      columns: [(col, ref, known, code, {base: mask})]: each code column, its reference base
        and the bases the genotypes know there with the samples that carry each, all in the
        target's orientation; ``known`` False where no genotype covers it (``code``: the bases
        of its IUPAC code).
      carriers: the mask of the indel's carriers.
      everyone: the mask a column no genotype covers grants (the carriers).
      max_dna_bulges, max_rna_bulges: the search's bulge limits.
      code_bases: callable(letter) -> the bases a PAM letter admits (default PAM_ADMITS).
    """
    admits_of = code_bases or (lambda letter: PAM_ADMITS.get(letter.upper(), ""))
    if (max_dna_bulges is not None and guide_seq.count("-") > max_dna_bulges) or (
            max_rna_bulges is not None and target_seq.count("-") > max_rna_bulges):
        return False
    if not carriers:
        return False
    by_col = {c[0]: c for c in columns}
    base_mm = 0  # the mismatches with the reference base at every code column
    fix_masks = []
    for j in spacer_cols:
        g = guide_seq[j]
        if j not in by_col:
            t = target_seq[j]
            if g != "-" and t != "-" and g.upper() != t.upper():
                base_mm += 1
            continue
        if g == "-":
            continue  # a DNA bulge: no base is a mismatch there
        _col, ref, known, code, base_masks = by_col[j]
        gu = g.upper()
        if ref == gu:
            continue
        base_mm += 1
        if not known:
            fix_masks.append(everyone if gu in code else 0)
        else:
            fix_masks.append(base_masks.get(gu, 0))
    pam_masks = []
    for j, letter in motif_at.items():
        allowed = admits_of(letter)
        if j not in by_col:
            t = target_seq[j].upper()
            if t not in "ACGT" or t not in allowed:
                return False
            continue
        _col, ref, known, code, base_masks = by_col[j]
        if ref in allowed:
            continue
        if not known:
            if not (set(code) | {ref}) & set(allowed):
                return False
            continue
        m = 0
        for b, bm in base_masks.items():
            if b in allowed:
                m |= bm
        pam_masks.append(m)
    if max_total is not None and max_total >= 0:
        limit = min(max_mm, max_total - bulge_size)
    else:
        limit = max_mm
    need = base_mm - limit
    return sample_can_produce(need, fix_masks, pam_masks, carriers)
