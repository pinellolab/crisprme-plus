"""Pure retrieveFromDict machinery for new_simple_analysis (CRISPRme+ dictless
redesign, Phase 1 wiring).

This module holds the SIDE-EFFECT-FREE core of ``new_simple_analysis.retrieveFromDict``
so it is importable and unit-testable in isolation (new_simple_analysis.py itself
opens ``sys.argv`` files at import time). ``new_simple_analysis`` imports these and
its ``retrieveFromDict`` is a thin wrapper that supplies the module-level registry
reader / mydict entry / current_chr.

The 5-tuple contract (see new_simple_analysis.py ~158) is:
  (snp_list, sample_list, rsID_list, AF_list, snp_info_list)
  * snp_list      : alt allele per '$'-entry / per registry alt
  * sample_list   : per-alt carrier list ("sampleID:genotype" tokens) or [] each
  * rsID_list     : variant id
  * AF_list       : allele-frequency STRING
  * snp_info_list : "<chrom>_<pos1based>_<ref>_<alt>"

STDLIB ONLY. ``tier0_registry`` is passed in (never imported here) so this module
has no hard dependency and cannot break the legacy path.
"""

from __future__ import annotations


def legacy_5tuple_from_entry(entry, chrom, chr_pos):
    """The EXACT legacy dict-entry decode: turn one raw dict value ('$'-joined
    per-ALT entries) into the 5-tuple. Verbatim lift of the original
    ``retrieveFromDict`` hit-branch body -- byte-identical so the no-registry path
    is unchanged (including its IndexError on a malformed entry)."""
    multi_entry = entry.split("$")
    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    for entry in multi_entry:
        split_entry = entry.split(";")
        samples = split_entry[0].strip().split(",")
        if samples[0] == "":
            samples = []
        sample_list.append(samples)
        snp_list.append(split_entry[1].strip().split(",")[1])
        rsID_list.append(split_entry[2].strip())
        AF_list.append(split_entry[3].strip())
        snp_info_list.append(
            chrom
            + "_"
            + str(chr_pos + 1)
            + "_"
            + split_entry[1].split(",")[0]
            + "_"
            + split_entry[1].split(",")[1]
        )
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list


def no_entry_5tuple(chrom, chr_pos):
    """The EXACT legacy no-entry fallback: one fake C>G SNP with no samples.

    Verbatim lift of the original ``retrieveFromDict`` except-branch body."""
    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    sample_list.append([])  # no samples
    snp_list.append("C")  # fake snp
    rsID_list.append(".")  # no rsid
    AF_list.append("0")  # fake AF
    snp_info_list.append(
        chrom + "_" + str(chr_pos + 1) + "_" + "C" + "_" + "G"
    )  # fake snp info list
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list


def dict_alt_to_samples(entry):
    """Map alt allele -> its carrier list from a raw dict value, plus first-seen
    alt order. Parses the '$'-joined per-ALT entries the same way the legacy
    decode does; used by the registry path to reattach the dict's carrier lists
    (preserving per-'$'-alt alignment) to registry-ordered alts. Returns ({}, [])
    for a falsy entry."""
    alt_to_samples = {}
    alt_order = []
    if not entry:
        return alt_to_samples, alt_order
    for sub in entry.split("$"):
        split_entry = sub.split(";")
        if len(split_entry) < 2:
            continue
        refalt = split_entry[1].strip().split(",")
        if len(refalt) < 2:
            continue
        alt = refalt[1]
        samples = split_entry[0].strip().split(",")
        if samples[0] == "":
            samples = []
        if alt not in alt_to_samples:
            alt_to_samples[alt] = samples
            alt_order.append(alt)
    return alt_to_samples, alt_order


def retrieve_5tuple(reader, entry, chrom, chr_pos, global_group_id="global",
                    gtreader=None):
    """Selection-matrix-driven builder of the retrieveFromDict 5-tuple. PURE.

    Args:
      reader: a ``tier0_registry.RegistryReader`` for ``chrom`` or None.
      entry:  the raw ``mydict`` value ('$'-joined per-ALT string) for this
              position, or None if the dict has no entry for it.
      chrom, chr_pos: chromosome name + 0-based position (dict keys are 1-based,
              i.e. ``chr_pos + 1``).
      global_group_id: the registry's GLOBAL group id (``tier0_registry.GLOBAL_GROUP_ID``).
      gtreader: a ``tier1_genotypes.GenotypeReader`` for ``chrom`` or None. Present
              only on a true DICTLESS install (Tier-0 registry + genotype tier, NO
              per-sample dict); supplies the per-alt "sampleID:genotype" carrier
              tokens the dict would have supplied. When None (any dict install, any
              old deploy), the Samples column comes from the dict (or is degraded
              to [] when neither a dict nor a genotype tier is present) -- so the
              legacy and registry+dict paths are BYTE-IDENTICAL to before.

    Selection matrix (WIRING spec):
      registry PRESENT + dict PRESENT                -> UNCHANGED: legacy dict
                                                        decode + registry AF
                                                        override. (gtreader unused:
                                                        the dict already carries the
                                                        exact carrier lists + alt
                                                        order the downstream
                                                        haplotype decomposition
                                                        depends on.)
      registry PRESENT + dict ABSENT + gtreader      -> DICTLESS FULL: iterate the
                                                        registry's alts; metadata/AF
                                                        from the registry, and
                                                        sample_list[i] =
                                                        gtreader.carrier_tokens(pos1,
                                                        alt) (or [] if None).
      registry PRESENT + dict ABSENT + NO gtreader   -> registry-only: metadata/AF
                                                        from registry, sample_list =
                                                        [] per alt (degraded
                                                        Samples), UNCHANGED.
      registry ABSENT                                -> LEGACY path EXACTLY as
                                                        today: the dict hit decode,
                                                        or the no-entry fake-SNP
                                                        fallback.

    ALT ORDER (dictless case): there is no dict to match, so alts are emitted in
    ``reader.alts_at(pos1)`` (sorted) order. The downstream haplotype decomposition
    is alt-order sensitive, so a dictless run MAY differ from a dict run at a few
    boundary bulge targets even though the carrier SETS match (validated separately
    in step 2d).

    5-tuple contract + element ordering + the no-entry fallback are preserved.
    """
    # --- registry ABSENT: byte-identical legacy path --------------------------
    if reader is None:
        if entry is None:
            return no_entry_5tuple(chrom, chr_pos)
        return legacy_5tuple_from_entry(entry, chrom, chr_pos)

    pos1 = chr_pos + 1  # registry is keyed by 1-based genomic position

    # --- registry PRESENT + dict PRESENT: AUGMENT the legacy decode -----------
    # Keep the dict's alt SET, alt ORDER, carrier samples, rsID and snp_info
    # BYTE-IDENTICAL to the legacy path (so the downstream haplotype decomposition,
    # clustering and the Samples column are unchanged), and OVERRIDE ONLY the AF
    # column with the registry's corrected global allele frequency (AC/AN over the
    # FULL panel) for each alt the registry knows. This makes registry+dict output
    # reproduce legacy EXACTLY except the (intentionally) corrected AF -- the legacy
    # dict AF field is empty for ~95% of variants and mis-polarized for some.
    # (Ordering the registry's alts independently would perturb the order-sensitive
    #  decomposition and can add/drop a boundary bulge target, so we do NOT.)
    if entry is not None:
        snp_list, sample_list, rsID_list, AF_list, snp_info_list = \
            legacy_5tuple_from_entry(entry, chrom, chr_pos)
        for i, alt in enumerate(snp_list):
            groups = reader.lookup(pos1, alt)
            if groups:
                gcnt = groups.get(global_group_id)
                if gcnt is not None:
                    AF_list[i] = "%.6g" % gcnt.allele_freq()
        return snp_list, sample_list, rsID_list, AF_list, snp_info_list

    # --- registry PRESENT + dict ABSENT: dictless / registry-only --------------
    # Metadata/AF from the registry, one 5-tuple element per registry alt (in
    # ``reader.alts_at`` sorted order -- there is no dict to match; see the alt-order
    # note in the docstring). For the Samples column:
    #   * gtreader PRESENT (true dictless: Tier-0 registry + genotype tier, NO
    #     per-sample dict) -> DICTLESS FULL: sample_list[i] =
    #     gtreader.carrier_tokens(pos1, alt) (the EXACT legacy "sampleID:genotype"
    #     token list per alt, ascending sample-index order), or [] when the tier has
    #     no record for this (pos, alt).
    #   * gtreader ABSENT (Tier-0-only install) -> registry-only, UNCHANGED:
    #     sample_list = [] per alt (degraded Samples). With an empty carrier list the
    #     downstream finalization guard in new_simple_analysis.py (`if len(samples) >
    #     0`) DROPS the target row, so a Tier-0-only (no dict, no genotype tier)
    #     install still emits NO variant rows -- standalone dictless variant output
    #     needs the genotype tier.
    # Per-alt alignment (snp/sample/rsID/AF/snp_info share the index order) is
    # preserved, and carrier lists are looked up PER ALT so there is no cross-alt
    # leak at a multiallelic site.
    alts = reader.alts_at(pos1)
    if not alts:
        return no_entry_5tuple(chrom, chr_pos)
    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    for alt in alts:
        groups = reader.lookup(pos1, alt)
        gcnt = groups.get(global_group_id) if groups else None
        af = gcnt.allele_freq() if gcnt is not None else 0.0
        rsid = reader.rsid(pos1, alt)
        if rsid is None:
            rsid = "."
        ref = reader.ref(pos1, alt) or "N"
        snp_list.append(alt)
        rsID_list.append(rsid)
        AF_list.append("%.6g" % af)
        snp_info_list.append(chrom + "_" + str(pos1) + "_" + ref + "_" + alt)
        if gtreader is not None:
            tokens = gtreader.carrier_tokens(pos1, alt)
            sample_list.append(tokens if tokens is not None else [])
        else:
            sample_list.append([])
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list
