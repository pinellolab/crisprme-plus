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

import itertools

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


def forward_snp_tokens(reg, chrom, used, global_group_id="global"):
    """Build FORWARD-STRAND variant tokens / rsIDs / MAFs for the cis SNP(s).

    Used by the INDEL post-analysis CONFIRMED-cis branch to enrich the indel row
    with the full real-haplotype variant set (indel + co-occurring SNPs) so the
    report's Variant column + observed_haplotypes companion list ALL the variants
    the off-target actually uses, not the indel alone.

    Args:
      reg: a Tier-0 RegistryReader (duck-typed: ``alts_at(pos1)`` -> list[str],
        ``ref(pos1, alt)`` -> forward-strand ref base, ``rsid(pos1, alt)`` ->
        rsID or None, ``lookup(pos1, alt)`` -> {group_id: Counts} with
        ``.allele_freq()``). All reads are FORWARD-strand by construction.
      chrom: chromosome string for the token prefix.
      used: the ``used_snps_for_target`` output -- a list of ``(real0, rsid, gt)``
        tuples where ``real0`` is a 0-based reference position (registry key =
        real0 + 1).
      global_group_id: the registry's GLOBAL group id (default "global").

    CRITICAL: the ``rsid``/``gt`` carried in ``used`` came from the caller's
    ``snp_at`` closure, whose REF/ALT are strand-COMPLEMENTED to the target strand
    (for the guide-vs-alt test only). We therefore re-read REF/ALT/rsID/AF straight
    from the registry here so every token is the forward-strand
    ``chrom_pos1_REF_ALT`` that resultIntegrator expects -- it splits the variant
    column on ',' and strand-corrects each token itself. The token / rsID / MAF
    format mirrors the SNP path EXACTLY
    (simple_analysis_registry.retrieve_5tuple):
      * token : ``chrom + "_" + str(pos1) + "_" + ref + "_" + alt``
      * rsID  : ``reg.rsid(pos1, alt)`` or "." when missing
      * MAF   : ``"%.6g" % global_allele_freq`` (same formatting as the SNP path)

    Returns three POSITION-ALIGNED lists ``(tokens, rsids, mafs)`` so the caller can
    append them 1:1 to the row's variant / rsID / MAF columns. A SNP whose registry
    record is absent (defensive -- it was just read in this scope) is skipped from
    all three lists together, so the 1:1 alignment is never broken.
    """
    tokens, rsids, mafs = [], [], []
    for (real0, _rs, _g) in used:
        pos1 = real0 + 1
        alts = reg.alts_at(pos1)
        if not alts:
            continue
        alt = alts[0]  # biallelic SNP position (same choice as the decode's snp_at)
        ref = reg.ref(pos1, alt) or "N"
        rsid = reg.rsid(pos1, alt)
        if rsid is None:
            rsid = "."
        groups = reg.lookup(pos1, alt)
        gcnt = groups.get(global_group_id) if groups else None
        af = gcnt.allele_freq() if gcnt is not None else 0.0
        tokens.append(chrom + "_" + str(pos1) + "_" + ref + "_" + alt)
        rsids.append(rsid)
        mafs.append("%.6g" % af)
    return tokens, rsids, mafs


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


# The IUPAC ambiguity chars an overlaid SNP position can carry in the target seq.
_IUPAC = set("RYSWKMBDHVryswkmbdhv")

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
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
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


def _gt_tokens(gt):
    """(allele tokens, kind) of a genotype: kind is "phased" ('|'), "unphased" ('/')
    or "haploid" (one token, e.g. a male chrX call)."""
    if "|" in gt:
        return gt.split("|"), "phased"
    if "/" in gt:
        return gt.split("/"), "unphased"
    return [gt], "haploid"


def kept_bases(ref, alt):
    """How many bases REF and ALT share at the start (the VCF anchor base): the REF
    bases a haplotype carrying the variant keeps. The rest of REF is deleted or
    replaced."""
    c = 0
    while c < min(len(ref), len(alt)) and ref[c].upper() == alt[c].upper():
        c += 1
    return c


def changes_span(pos1, ref, alt, lo, hi, gone=None):
    """True iff a variant (VCF POS / REF / ALT) changes the reference bases lo..hi
    (0-based, inclusive) of a haplotype that carries it: a deleted or replaced base
    inside, or bases inserted between two bases of it. The bases REF and ALT share
    at the start (the VCF anchor base) are kept, so they change nothing. ``gone``
    (first, last; 0-based, inclusive): bases the haplotype has lost already, to the
    indel of the window being read: losing them again changes nothing."""
    a = pos1 - 1
    c = kept_bases(ref, alt)
    if len(ref) > c:  # bases a+c .. a+len(ref)-1 are deleted or replaced
        first, last = a + c, a + len(ref) - 1
        parts = [(first, last)] if gone is None else [
            (first, min(last, gone[0] - 1)), (max(first, gone[1] + 1), last)]
        if any(f <= l and f <= hi and l >= lo for f, l in parts):
            return True
    if len(alt) > c:  # bases inserted between a+c-1 and a+c
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


def resolve_carriers(columns, indel_gts, phase_known=True, others=()):
    """Split the indel's carriers into exactly resolved slots and PUTATIVE samples.

    Args:
      columns: the target's code columns (SnpColumn), in column order.
      indel_gts: {sample: normalized indel genotype} ("1|0", "0|1", "1|1", "0/1",
        haploid "1"), or an IndelCarriers built from it.
      phase_known: False when the indel's carriers come without phase (the indel
        log, no phased indel genotype store): no slot can then be tied to a SNP.
        (An IndelCarriers carries its own.)
      others: one {sample: normalized genotype} per OTHER indel that changes the
        target's reference span (changes_span). The fake-indel window holds one
        indel, so a slot that carries a second one there does not have the window's
        sequence: it is no carrier of any row of this target. A sample whose second
        indel cannot be put on a slot (unphased) is PUTATIVE.

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
                continue  # every indel slot also holds a second indel here
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
    the overlay codes); malformed entries are skipped."""
    out = []
    for entry in str(value).split("$"):
        parts = entry.split(";")
        if len(parts) < 3:
            continue
        refalt = parts[1].strip().split(",")
        if len(refalt) < 2:
            continue
        ref, alt = refalt[0].strip().upper(), refalt[1].strip().upper()
        if len(ref) != 1 or len(alt) != 1 or ref not in _BASE_ORDER or alt not in _BASE_ORDER:
            continue
        carriers = {}
        for tok in parts[0].split(","):
            sid, sep, gt = tok.strip().partition(":")
            if sid and sep:
                carriers[sid] = gt.strip()
        rsid = parts[2].strip() or "."
        af = parts[3].strip() if len(parts) > 3 else ""
        out.append((ref, alt, rsid, af, carriers))
    return out


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
