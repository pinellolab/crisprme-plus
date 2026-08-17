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


def retrieve_5tuple(reader, entry, chrom, chr_pos, global_group_id="global"):
    """Selection-matrix-driven builder of the retrieveFromDict 5-tuple. PURE.

    Args:
      reader: a ``tier0_registry.RegistryReader`` for ``chrom`` or None.
      entry:  the raw ``mydict`` value ('$'-joined per-ALT string) for this
              position, or None if the dict has no entry for it.
      chrom, chr_pos: chromosome name + 0-based position (dict keys are 1-based,
              i.e. ``chr_pos + 1``).
      global_group_id: the registry's GLOBAL group id (``tier0_registry.GLOBAL_GROUP_ID``).

    Selection matrix (WIRING spec):
      registry PRESENT + dict PRESENT -> metadata/AF from the registry, carrier
                                         samples from the dict (per-alt aligned).
      registry PRESENT + dict ABSENT  -> registry-only: metadata/AF from registry,
                                         sample_list = [] per alt (degraded Samples;
                                         the genotype tier is a later phase).
      registry ABSENT                 -> LEGACY path EXACTLY as today: the dict hit
                                         decode, or the no-entry fake-SNP fallback.

    5-tuple contract + element ordering + the no-entry fallback are preserved.
    """
    # --- registry ABSENT: byte-identical legacy path --------------------------
    if reader is None:
        if entry is None:
            return no_entry_5tuple(chrom, chr_pos)
        return legacy_5tuple_from_entry(entry, chrom, chr_pos)

    # --- registry PRESENT -----------------------------------------------------
    pos1 = chr_pos + 1  # registry is keyed by 1-based genomic position
    alts = reader.alts_at(pos1)
    if not alts:
        # Registry has NO record at this position. Fall back to the legacy dict
        # behavior so we never lose an off-target the dict knew about.
        if entry is None:
            return no_entry_5tuple(chrom, chr_pos)
        return legacy_5tuple_from_entry(entry, chrom, chr_pos)

    # Dict carriers keyed by alt (empty on a registry-only / Tier-0-only install).
    # One 5-tuple element per registry alt, in the STABLE registry alt order so
    # downstream alignment against sample_list holds.
    alt_to_samples, _ = dict_alt_to_samples(entry)

    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    for alt in alts:
        groups = reader.lookup(pos1, alt)
        gcnt = groups.get(global_group_id) if groups else None
        # AF from the registry is the CORRECTED global allele frequency AC/AN
        # (re-derived over the FULL panel). This INTENTIONALLY differs from the
        # legacy dict AF field: that field is empty for ~95% of variants and
        # mis-polarized for some. See the Tier-0 registry design notes.
        af = gcnt.allele_freq() if gcnt is not None else 0.0
        rsid = reader.rsid(pos1, alt)
        if rsid is None:
            rsid = "."
        snp_list.append(alt)
        rsID_list.append(rsid)
        # %.6g keeps the AF column compact (like the legacy short strings) rather
        # than a 17-sig-fig float repr; the value is the corrected AC/AN.
        AF_list.append("%.6g" % af)
        snp_info_list.append(chrom + "_" + str(pos1) + "_" + "N" + "_" + alt)
        if entry is not None and alt in alt_to_samples:
            # dict present: reuse its carrier list for this alt (per-'$'-alt aligned)
            sample_list.append(alt_to_samples[alt])
        else:
            # registry-only / dict has no carriers for this alt: degraded Samples.
            # SCOPE NOTE (Phase 1): with an empty carrier list the downstream
            # finalization guard in new_simple_analysis.py (`if len(samples) > 0`)
            # DROPS the target row entirely -- so a true dictless (Tier-0-only, no
            # per-sample dict) install emits NO variant rows yet. Phase 1 therefore
            # AUGMENTS the dict path (corrected AF/rsID from the registry); standalone
            # dictless variant output requires the per-sample GENOTYPE TIER (Phase 2/3).
            sample_list.append([])
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list
