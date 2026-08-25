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
"""

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


def used_snps_for_target(guide_seq, target_seq, offset_to_real, snp_at):
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

    Returns a list of (real_pos, rsID, alt_gt_by_sample) for every column where the
    target carries an IUPAC ambiguity AND the guide base equals the SNP ALT (not the
    REF) -- i.e. the off-target REQUIRES the alt allele (a genuine SNP+indel
    co-occurrence), not merely tolerates the reference. Ref-satisfied ambiguities
    are skipped (the SNP is incidental, not enabling).
    """
    used = []
    for j, (g, t) in enumerate(zip(guide_seq, target_seq)):
        if t not in _IUPAC:
            continue
        gu = g.upper()
        if gu in ("-", "N"):
            continue
        real = offset_to_real(j)
        if real is None:
            continue
        snp = snp_at(real)
        if snp is None:
            continue
        ref_base, alt_base, rsid, gt = snp
        # USES the alt iff the guide matches the alt but not the ref
        if gu == alt_base.upper() and gu != ref_base.upper():
            used.append((real, rsid, gt))
    return used
