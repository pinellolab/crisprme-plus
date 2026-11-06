"""Observed-haplotype enumerator for the dict-less multi-variant off-target path
(CRISPRme+ dict-less redesign, sensitivity-first phased/unphased handling).

WHY THIS MODULE EXISTS
----------------------
The legacy ``new_simple_analysis.iupac_decomposition`` decomposes an IUPAC (variant)
off-target window into concrete alt-allele haplotypes by (a) splitting every
genotype on ``"|"`` (which mis-places an UNPHASED ``0/1`` onto BOTH haplotype seeds),
(b) growing a full 2^k combination lattice, (c) taking a SINGLE global phased-vs-
unphased decision (``haplotype_check``) rather than a PER-SAMPLE one, and (d) above
``IUPAC_CAP`` falling back to one greedy representative that UNDER-REPORTS risk.

When a Tier-1 genotype store is present (``mygt is not None``, the true dict-less
install) we can do far better: read the genotypes ONCE per window and enumerate
exactly the DISTINCT haplotypes / variant-sets that occur in at least one REAL
individual, each carrying its EXACT carriers. This module is that enumerator, kept
PURE + SIDE-EFFECT-FREE (no ``sys.argv`` / file I/O at import) so it is unit-testable
in isolation, exactly like ``simple_analysis_registry``. ``new_simple_analysis``
imports it and its dict-less branch drives ``enumerate_observed_haplotypes``.

SEMANTICS (LOCKED DESIGN, sensitivity-first -- never miss a putative region)
----------------------------------------------------------------------------
For each candidate window with k variant positions, for each REAL individual that
carries >= 1 alt in the window:
  * PHASED (its GT strings use ``"|"``) AND all covered positions share a phase set
    (PS id equal, or PS absent => single whole-chromosome block) -> group the alts by
    HAPLOTYPE SLOT; each non-empty slot is one CONFIRMED cis variant-set (a real cis
    haplotype the individual carries).
  * UNPHASED (its GT strings use ``"/"``), OR its alleles span DIFFERENT phase sets,
    OR its genotypes mix phased/unphased across the window -> collapse the
    individual's ENTIRE diploid alt collection into ONE PUTATIVE variant-set (phase
    unconfirmed, but NEVER dropped, NEVER split onto two seeds).
Then dedup across individuals into DISTINCT observed variant-sets, union carriers per
set, and mark a set CONFIRMED iff EVERY contributing (individual, path) reached it via
a phased same-PS cis path -- else PUTATIVE (a set reachable both ways is PUTATIVE, the
conservative choice). Cross-individual chimeras (allele combinations no single
individual carries) are EXCLUDED by construction: a variant-set key only exists if a
real sample produced exactly it. The output size is therefore bounded by the number of
DISTINCT observed sets -- at most ~ploidy * n_carriers -- NOT 2^k.

MULTIALLELIC POSITIONS (CRISPRitz 2.10.0 codes every ALT of a split site)
-------------------------------------------------------------------------
The tiers store one (pos, alt) record per alt, and a position can hold two or three of
them. The caller hands each column's genotypes already read against each alt's own
allele index (``site_alleles.per_alt_entries``: "1" holds that alt, also for a record
that listed several ALTs on one VCF line, whose raw "1|2" becomes "1|0" under the first
alt and "0|1" under the second), so ``alt_index`` is "1" there. Two rules keep a set to
one alt per column:
  * a phased sample whose slot holds two alts of one column (split records phased one
    by one: 672 sample-positions on chr22 of HGDP + 1000 Genomes, e.g. NA18869 0|1 under
    G and 0|1 under T at chr22:10,661,629) is no haplotype there: it takes the PUTATIVE
    branch, as an unphased sample does, and is a carrier of {G} and of {T};
  * the PUTATIVE branch's combinations, the union fallback above ``max_putative``
    included, take at most one alt of a column.
The phase label is per set, so a set such a sample shares reads PUTATIVE for all its
carriers.

STDLIB ONLY. ``population_summary`` is imported for ``parse_haplotypes`` /
``_token_is_alt`` ONLY (one source of truth for allele-slot math); that import is
guarded at the single call site in ``new_simple_analysis`` so an old deploy without it
cannot break the legacy path (the dict-less branch is simply not taken).
"""

from __future__ import annotations

import itertools

CONFIRMED = "CONFIRMED"
PUTATIVE = "PUTATIVE"


class ObservedHaplotype(object):
    """One distinct variant-set observed in >= 1 real individual within a window.

    Fields:
      variant_set : frozenset[(pos_c:int, alt:str)] -- the ambiguity-column index
                    ``pos_c`` (0-based within the ungapped replaceTarget) paired with
                    the alt base applied there. Keys the dedup, so it is order- and
                    duplicate-independent.
      carriers    : set[str] -- the EXACT sample ids that carry this set (fixes the
                    Samples column and, downstream, the AF).
      phase_state : CONFIRMED (every carrier reached it via a phased same-PS cis path)
                    or PUTATIVE (>= 1 carrier reached it unphased / cross-PS).
      seq         : list[char] -- refSeq (pre-revert) with each (pos_c, alt) applied,
                    mirroring the legacy ``listReplaceTarget`` construction so it drops
                    straight into the existing finalization block.
      info        : list[[rsID, AF, snp_info]] per applied (pos_c, alt), in ascending
                    pos_c order (the legacy per-position info-row shape).
    """

    __slots__ = ("variant_set", "carriers", "phase_state", "seq", "info")

    def __init__(self, variant_set, carriers, phase_state, seq, info):
        self.variant_set = variant_set
        self.carriers = carriers
        self.phase_state = phase_state
        self.seq = seq
        self.info = info

    def __repr__(self):
        return ("ObservedHaplotype(vset=%r, carriers=%r, phase=%s)"
                % (sorted(self.variant_set), sorted(self.carriers),
                   self.phase_state))


def _gt_is_phased(gt):
    """True iff a genotype string is PHASED (uses '|'), False if unphased ('/') or
    haploid / single-token (no separator -> treat as unconfirmed-phase, PUTATIVE)."""
    if gt is None:
        return False
    return "|" in gt


def _per_sample_variant_sets(positions, parse_haplotypes, ploidy_of, chrom,
                             max_putative=-1):
    """Build each real individual's per-phase-set cis variant-set(s) for the window.

    Args:
      positions: list of column descriptors, one per IUPAC/variant column, each a dict
        with keys:
          "pos_c"       : 0-based ambiguity-column index within replaceTarget.
          "alts"        : list of the column's alt bases (alignment index i).
          "carrier_gts" : list (index i) of {sample_id: gt_string} for that alt.
          "alt_index"   : list (index i) of the resolved alt-index token per alt.
          "ps"          : list (index i) of {sample_id: ps_id} phase-set maps, or
                          None when PS is not stored (single-block assumption).
      parse_haplotypes: ``population_summary.parse_haplotypes`` (allele-slot split).
      ploidy_of: callable(sample_id, sex) -> ploidy (0/1/2). ``sex`` is unknown here
        (the enumerator has no per-sample sex axis at this layer) so we pass None; the
        autosomal model returns 2 regardless, and chrX/Y refinements are a documented
        follow-up (they only shrink slots, never invent carriers).
      chrom: chromosome name (unused today; accepted for a future PS/sex axis).

    Returns dict[sample_id -> list of (frozenset[(pos_c, alt)], phase_state)].
    Each sample contributes:
      * one CONFIRMED set per non-empty haplotype slot when ALL its covered columns
        are phased AND share a phase set AND no slot holds two alts of one column;
      * else PUTATIVE sets: every combination of the variants it carries that takes
        at most one alt of a column (above ``max_putative`` variants, the union of
        everything it carries, one per choice of alt at a column holding two).
    """
    # 1) Gather, per sample, its (pos_c, alt, slot-tokens, phased?, ps) at each column
    #    it carries. We key alts by (pos_c, alt) so a multiallelic column with two alts
    #    contributes two entries for a "1|2" carrier.
    per_sample = {}  # sid -> list of dicts: {pos_c, alt, toks, phased, ps}
    for col in positions:
        pos_c = col["pos_c"]
        for i, alt in enumerate(col["alts"]):
            alt_index = col["alt_index"][i]
            ps_map = col["ps"][i] if col["ps"] is not None else None
            for sid, gt in col["carrier_gts"][i].items():
                ploidy = ploidy_of(sid, None)
                toks = parse_haplotypes(gt, ploidy)
                phased = _gt_is_phased(gt)
                ps_id = ps_map.get(sid) if ps_map is not None else None
                per_sample.setdefault(sid, []).append({
                    "pos_c": pos_c,
                    "alt": alt,
                    "alt_index": alt_index,
                    "toks": toks,
                    "phased": phased,
                    "ps": ps_id,
                })

    out = {}
    for sid, entries in per_sample.items():
        # PUTATIVE trigger: any column unphased, OR the covered phased columns do not
        # all share ONE phase set. PS absent => a single whole-chromosome block => all
        # phased columns are cis (documented assumption).
        all_phased = all(e["phased"] for e in entries)
        ps_ids = set(e["ps"] for e in entries if e["phased"])
        # Ignore a lone None (no PS tag = single block). A mix of {None, some real PS}
        # or {two real PS} means we cannot prove cis across the window -> PUTATIVE.
        real_ps = set(p for p in ps_ids if p is not None)
        mixed_ps = (len(real_ps) > 1) or (len(real_ps) >= 1 and None in ps_ids)

        slots = None
        if all_phased and not mixed_ps:
            # CONFIRMED: group alts by haplotype slot. A slot h carries (pos_c, alt)
            # iff this sample's token on slot h at that column equals the alt's index.
            slots = {}  # slot -> set[(pos_c, alt)]
            for e in entries:
                toks = e["toks"]
                for h, tok in enumerate(toks):
                    if tok == e["alt_index"]:
                        slots.setdefault(h, set()).add((e["pos_c"], e["alt"]))
            # A slot holding two alts of one column is no haplotype (its split records
            # were phased one by one): the sample is PUTATIVE in this window, as the
            # indel path's resolve_carriers reads such a slot.
            for vset in slots.values():
                if len(set(pos_c for pos_c, _alt in vset)) < len(vset):
                    slots = None
                    break
        if slots is not None:
            sets = []
            seen_keys = set()
            for h, vset in slots.items():
                if not vset:
                    continue
                key = frozenset(vset)
                if key in seen_keys:
                    continue  # two identical slots (hom) => one CONFIRMED set
                seen_keys.add(key)
                sets.append((key, CONFIRMED))
            if sets:
                out[sid] = sets
            # A phased sample with no slot actually carrying an alt (all ref/missing)
            # contributes nothing -- it is not a carrier of any set.
        else:
            # PUTATIVE: the sample is unphased (or its phased columns span >1 phase set),
            # so the cis/trans phasing of its carried variants is UNKNOWN. Emit EVERY
            # non-empty SUBSET of the variants it carries as a candidate cis haplotype --
            # any subset could be the sample's true haplotype, and a variant that BREAKS a
            # target (e.g. one that disrupts the PAM) MUST be droppable, so the maximal
            # union alone is NOT sufficient (it hides a valid sub-combination off-target
            # that needs a reference allele at a het column). The downstream finalize
            # mm/PAM gates prune out-of-budget / PAM-invalid subsets, and
            # enumerate_observed_haplotypes dedups subsets shared across samples. These are
            # the sample's OWN variants only -- never cross-individual chimeras.
            # GUARD: a sample carrying more than ``max_putative`` variants in the window
            # falls back to the union (the 2^k blow-up is confined to that one sample; a
            # dense window is already surfaced to the high-variant-density BED upstream).
            # A haplotype holds one allele of a position, so no candidate takes two alts
            # of one column (a 1/2 heterozygote, or a phased slot holding two): the
            # combinations choose, per column, nothing or one of the alts carried there,
            # and the union is one per choice of alt at a column holding two.
            carried = sorted(set(
                (e["pos_c"], e["alt"]) for e in entries
                if any(tok == e["alt_index"] for tok in e["toks"])
            ))
            if carried:
                by_col = {}
                for pos_c, alt in carried:
                    by_col.setdefault(pos_c, []).append(alt)
                cols = sorted(by_col)
                if max_putative >= 0 and len(carried) > max_putative:
                    out[sid] = [
                        (frozenset(zip(cols, choice)), PUTATIVE)
                        for choice in itertools.product(*(by_col[c] for c in cols))
                    ]
                else:
                    out[sid] = [
                        (frozenset((c, a) for c, a in zip(cols, choice)
                                   if a is not None), PUTATIVE)
                        for choice in itertools.product(
                            *([None] + by_col[c] for c in cols))
                        if any(a is not None for a in choice)
                    ]
    return out


def enumerate_observed_haplotypes(positions, refSeq, parse_haplotypes, ploidy_of,
                                  chrom, max_putative=-1):
    """Enumerate the DISTINCT observed haplotypes for one candidate window.

    Args:
      positions: see ``_per_sample_variant_sets``. Each column also carries per-alt
        ``info`` = [rsID, AF, snp_info] used to build the emitted info rows.
      refSeq: the pre-revert reference window (str), len == len(replaceTarget). Each
        ObservedHaplotype.seq is a list(refSeq) copy with its (pos_c, alt) applied.
      parse_haplotypes, ploidy_of, chrom: as in ``_per_sample_variant_sets``.

    Returns list[ObservedHaplotype], one per distinct observed variant-set, carriers
    unioned across samples, phase_state downgraded to PUTATIVE if ANY contributor was
    unphased / cross-PS. EXCLUDES cross-individual chimeras (a set only exists if a
    real sample produced exactly it). Deterministic order: sorted by the set's sorted
    (pos_c, alt) tuple so the emitted rows are reproducible.
    """
    # Per-(pos_c, alt) info lookup so each haplotype's info rows are built in ascending
    # pos_c order regardless of which sample produced the set.
    info_of = {}
    for col in positions:
        pos_c = col["pos_c"]
        for i, alt in enumerate(col["alts"]):
            info_of[(pos_c, alt)] = col["info"][i]

    per_sample = _per_sample_variant_sets(positions, parse_haplotypes, ploidy_of,
                                          chrom, max_putative=max_putative)

    # Dedup: variant_set -> {carriers:set, confirmed_all:bool}. A set is CONFIRMED only
    # if EVERY contributing (sample, path) was CONFIRMED cis.
    dedup = {}
    for sid, sets in per_sample.items():
        for vset, phase_state in sets:
            slot = dedup.get(vset)
            if slot is None:
                slot = {"carriers": set(), "confirmed_all": True}
                dedup[vset] = slot
            slot["carriers"].add(sid)
            if phase_state != CONFIRMED:
                slot["confirmed_all"] = False

    out = []
    for vset in sorted(dedup, key=lambda s: sorted(s)):
        slot = dedup[vset]
        seq = list(refSeq)
        info = []
        for (pos_c, alt) in sorted(vset):
            if pos_c < len(seq):
                seq[pos_c] = alt
            info.append(list(info_of.get((pos_c, alt), [".", "0", "."])))
        phase_state = CONFIRMED if slot["confirmed_all"] else PUTATIVE
        out.append(ObservedHaplotype(
            variant_set=vset,
            carriers=set(slot["carriers"]),
            phase_state=phase_state,
            seq=seq,
            info=info,
        ))
    return out
