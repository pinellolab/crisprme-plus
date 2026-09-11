"""Combination-aware population summaries (CRISPRme+ dictless redesign, Phase 3).

This module is the PURE computation heart of the "population summary" the
dictless post-analysis attaches to every off-target. It has NO pipeline I/O: it
consumes ONLY the two already-built readers (Tier-0 registry + Tier-1 genotype
store) and a sample axis, and returns a ``PopSummary`` value object.

WHY combination-aware (the design ground truth, see the design memory
``crisprme-plus-dictless-population-design`` #3):

    An off-target's FREQUENCY is a HAPLOTYPE (co-occurrence) property, NOT a
    per-SNP allele frequency. For a multi-SNP off-target the frequency is the
    frequency of the SPECIFIC allele COMBINATION that CREATES the off-target
    (the SNPs may be in LD; the product of per-SNP AFs is wrong). The caller
    passes the (pos, alt) pairs whose ALT creates / strengthens the off-target.

Two paths, keyed on k == len(pos_alts):

  * k == 1  (SINGLE variant): the exact per-group Counts already live in the
    Tier-0 registry, computed full-panel and ploidy/sex-aware at build time. We
    read ``tier0_reader.lookup(pos, alt)`` and derive the views + record-level
    stats via ``derive_record_stats`` -- byte-for-byte the already-built path.

  * k >= 2  (COMBINATION): Tier-0 stores only PER-SNP aggregates, which cannot
    express co-occurrence. We consult the Tier-1 genotype store for the CARRIERS
    of each (pi, ai) and intersect them, honoring phase:

      phased=True  -- an individual carries the combination iff SOME haplotype
                      h in {0,1} carries ALL ai IN CIS: for that h, the
                      individual's allele on haplotype h equals ai at EVERY pi.
                      A "1|1" hom carries on both haplotypes; a haploid "1"
                      contributes exactly one haplotype. allele (haplotype)
                      frequency is well-defined: (# cis-haplotypes carrying all
                      ai) / AN. hom(combination) = carries on BOTH haplotypes.

      phased=False -- phase is unknown, so we CANNOT prove cis. We report the
                      ASSUME-CIS UPPER BOUND on carriers: an individual is a
                      (bound) carrier iff it is het-or-hom for ai at EVERY pi.
                      allele frequency is UNDEFINED (reported as None). We ALSO
                      report a LOWER BOUND (see ``_unphased_carrier_lower_bound``).
                      hom(combination) = hom for ai at EVERY pi (an upper bound
                      on true cis-homozygosity).

DENOMINATORS (both k paths, per group G):
    n_called = number of individuals in G on this chromosome (ploidy > 0). The
               per-sample dict has NO missing genotypes, so every panel member is
               "called".
    AN       = sum over G of ploidy_of(sample_id, sex)  (ploidy/sex-aware;
               chrX-nonPAR male = 1, chrY female = 0/absent, etc.). Computed ONCE
               from the axis.
    GLOBAL   = the DEDUPLICATED union of canonical sample_ids (a sample present in
               two databases is counted once globally).

NEVER conflate across databases: per-db and per-(db x subpop) groups use that
db's native labels. max_subpop_af = argmax allele_freq over (db x subpop) groups
with a deterministic tie-break (alphabetically-first group_id). observed =
global n_carrier >= 1.

STDLIB ONLY. PURE (no file I/O beyond the passed readers). Not wired into
new_simple_analysis.py here (that is step 3c).
"""

from __future__ import annotations

import math

from tier0_registry import (
    Counts,
    derive_record_stats,
    GLOBAL_GROUP_ID,
    SEP,
    db_group_id,
    db_subpop_group_id,
)

# NaN sentinel for an UNDEFINED allele frequency (unphased combination). We use a
# real float('nan') AND expose the ``is_undefined`` predicate on GroupSummary so
# callers never accidentally treat NaN as 0.0.
AF_UNDEFINED = float("nan")


# --------------------------------------------------------------------------- #
# Genotype / haplotype parsing (reuses the tier0/tier1 token conventions).
# --------------------------------------------------------------------------- #
def parse_haplotypes(gt, ploidy):
    """Split a genotype string into its per-HAPLOTYPE allele tokens.

    Honors the tier0/tier1 conventions:
      * "a0|a1" (phased) / "a0/a1" (unphased) / "a0" (haploid) all split into
        their raw tokens ("0", "1", ".", or a foreign alt index like "2").
      * ploidy caps the number of haplotype slots we honor. We NEVER pad a short
        GT by duplicating a token (a haploid "1" declared diploid yields ONE
        slot), and we NEVER count more than ``ploidy`` slots (a diploid-style
        "1|1" declared haploid yields ONE slot). This mirrors
        ``tier0_registry._sample_alleles`` exactly so AN/AC stay consistent
        with the registry's full-panel counts.

    Returns a list of raw allele-token strings, length <= ploidy.
    """
    if gt is None:
        toks = []
    else:
        s = gt.strip()
        if s == "":
            toks = []
        elif "|" in s:
            toks = s.split("|")
        elif "/" in s:
            toks = s.split("/")
        else:
            toks = [s]
    if ploidy <= 1:
        return toks[:1]
    return toks[:ploidy]


def _token_is_alt(token, alt_index="1"):
    """True iff a raw allele token IS this record's alt.

    Mirrors ``tier0_registry._classify``: "." is missing (not alt), the
    ``alt_index`` token (default "1" under the bcftools norm -m- biallelic
    convention) is THIS record's alt, and ANY OTHER token -- ref "0" or a foreign
    alt index like "2" -- is NOT this record's alt.
    """
    if token == ".":
        return False
    return token == alt_index


def individual_carries_at(gt, ploidy, alt_index="1"):
    """Unphased carrier test at ONE (pos, alt): het-or-hom for the alt.

    True iff at least one honored haplotype slot carries the record's alt.
    """
    for tok in parse_haplotypes(gt, ploidy):
        if _token_is_alt(tok, alt_index):
            return True
    return False


def individual_hom_at(gt, ploidy, alt_index="1"):
    """Unphased hom test at ONE (pos, alt): fully-called AND every slot is alt.

    Matches ``tier0_registry._sample_counts``'s ``is_hom``: fully called (no
    missing slot among ``ploidy`` slots) and every called slot is this alt. A
    haploid "1" is hom (hemizygous carrier).
    """
    toks = parse_haplotypes(gt, ploidy)
    if len(toks) != ploidy or ploidy == 0:
        return False
    for tok in toks:
        if tok == ".":
            return False
        if not _token_is_alt(tok, alt_index):
            return False
    return True


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #
class GroupSummary(object):
    """Per-group population summary for ONE off-target (single or combination).

    Fields:
      allele_count (AC) : alt HAPLOTYPES carrying the combination in this group.
                          For k==1 this is the registry AC. For a phased k>=2
                          combination it is the number of cis-haplotypes carrying
                          ALL ai. For an unphased k>=2 combination it is None
                          (allele-level counting is UNDEFINED without phase).
      allele_number (AN): full-panel called alleles in the group (sum of ploidy).
      n_carrier         : individuals carrying the combination (>=1 haplotype in
                          cis when phased; het-or-hom at every pi -- the
                          ASSUME-CIS UPPER BOUND -- when unphased).
      n_carrier_lower   : unphased LOWER bound on n_carrier (None when phased,
                          where n_carrier is exact).
      n_hom             : individuals homozygous for the combination.
      n_called          : full-panel called individuals in the group.

    Derived (properties): allele_freq (AC/AN, or NaN sentinel if AC is None),
    carrier_freq (n_carrier/n_called), hom_freq (n_hom/n_called).
    """

    __slots__ = ("allele_count", "allele_number", "n_carrier",
                 "n_carrier_lower", "n_hom", "n_called")

    def __init__(self, allele_count, allele_number, n_carrier, n_hom, n_called,
                 n_carrier_lower=None):
        self.allele_count = allele_count
        self.allele_number = allele_number
        self.n_carrier = n_carrier
        self.n_carrier_lower = n_carrier_lower
        self.n_hom = n_hom
        self.n_called = n_called

    @property
    def allele_freq_defined(self):
        """True iff allele frequency is well-defined (AC is not None)."""
        return self.allele_count is not None

    @property
    def allele_freq(self):
        """AC / AN, or the NaN sentinel when undefined (unphased combination)."""
        if self.allele_count is None:
            return AF_UNDEFINED
        return (self.allele_count / self.allele_number) if self.allele_number else 0.0

    @property
    def carrier_freq(self):
        return (self.n_carrier / self.n_called) if self.n_called else 0.0

    @property
    def carrier_freq_lower(self):
        if self.n_carrier_lower is None:
            return None
        return (self.n_carrier_lower / self.n_called) if self.n_called else 0.0

    @property
    def hom_freq(self):
        return (self.n_hom / self.n_called) if self.n_called else 0.0

    def as_tuple(self):
        return (self.allele_count, self.allele_number, self.n_carrier,
                self.n_carrier_lower, self.n_hom, self.n_called)

    def __eq__(self, other):
        return isinstance(other, GroupSummary) and self.as_tuple() == other.as_tuple()

    def __repr__(self):
        return ("GroupSummary(AC={}, AN={}, n_carrier={}, n_carrier_lower={}, "
                "n_hom={}, n_called={})").format(*self.as_tuple())


class PopSummary(object):
    """Combination-aware population summary for ONE off-target.

    Fields:
      pos_alts      : the (pos, alt) pairs that jointly create the off-target.
      k             : len(pos_alts) (1 == single variant, >=2 == combination).
      phased        : the phasing mode used (bool).
      groups        : dict[group_id] -> GroupSummary for GLOBAL, each <db>, and
                      each <db>::<subpop> group that HAS >= 1 carrier (SPARSE,
                      mirroring the Tier-0 sparsity contract). A group with zero
                      carriers is omitted (its freqs are 0 by construction).
      global_summary: the GLOBAL GroupSummary, or None if no global carrier.
      allele_freq_defined : False for an unphased combination (k>=2), else True.

    Record-level (mirrors ``derive_record_stats``):
      global_af             : GLOBAL allele_freq (NaN sentinel if undefined).
      global_carrier_freq   : GLOBAL n_carrier / n_called.
      global_hom_freq       : GLOBAL n_hom / n_called.
      global_carrier_n      : ABSOLUTE global carrier count (deduped union).
      max_subpop_af         : max allele_freq over (db x subpop) groups.
      max_subpop_af_label   : the group_id achieving that max (None if none /
                              undefined).
      observed              : True iff GLOBAL n_carrier >= 1.
    """

    __slots__ = ("pos_alts", "k", "phased", "groups", "global_summary",
                 "allele_freq_defined", "global_af", "global_carrier_freq",
                 "global_hom_freq", "global_carrier_n", "max_subpop_af",
                 "max_subpop_af_label", "observed")

    def __init__(self, pos_alts, phased, groups):
        self.pos_alts = list(pos_alts)
        self.k = len(self.pos_alts)
        self.phased = bool(phased)
        self.groups = groups

        g = groups.get(GLOBAL_GROUP_ID)
        self.global_summary = g
        self.allele_freq_defined = (g.allele_freq_defined if g is not None
                                    else (self.k == 1 or self.phased))
        self.global_af = g.allele_freq if g is not None else (
            0.0 if self.allele_freq_defined else AF_UNDEFINED)
        self.global_carrier_freq = g.carrier_freq if g is not None else 0.0
        self.global_hom_freq = g.hom_freq if g is not None else 0.0
        self.global_carrier_n = g.n_carrier if g is not None else 0
        self.observed = bool(g is not None and g.n_carrier >= 1)

        self.max_subpop_af, self.max_subpop_af_label = _max_subpop_af(groups)

    def group_ids(self):
        return sorted(self.groups)

    def __repr__(self):
        return ("PopSummary(k={}, phased={}, observed={}, global_carrier_n={}, "
                "global_af={}, max_subpop_af={} @ {})").format(
                    self.k, self.phased, self.observed, self.global_carrier_n,
                    self.global_af, self.max_subpop_af, self.max_subpop_af_label)


def _max_subpop_af(groups):
    """argmax allele_freq over (db x subpop) groups; deterministic tie-break.

    Iterate group ids in SORTED order so ties resolve to the alphabetically-first
    (db x subpop) label. Groups whose allele freq is UNDEFINED (unphased
    combination) are skipped -- an undefined AF cannot be a max. Returns
    (max_af, label) with label None if no eligible (db x subpop) group.
    """
    max_af = 0.0
    max_label = None
    for gid in sorted(groups):
        if gid == GLOBAL_GROUP_ID:
            continue
        if SEP not in gid:  # a db group, not a (db x subpop) cell
            continue
        gs = groups[gid]
        if not gs.allele_freq_defined:
            continue
        af = gs.allele_freq
        if max_label is None or af > max_af:
            max_af = af
            max_label = gid
    return max_af, max_label


# --------------------------------------------------------------------------- #
# Panel (denominators) -- computed ONCE from the axis, ploidy/sex-aware.
# --------------------------------------------------------------------------- #
class Panel(object):
    """Per-group full-panel denominators (AN + n_called), ploidy/sex-aware.

    Built ONCE from the sample axis. For each group -- every (db x subpop) cell,
    every db, and the GLOBAL group over the DEDUPED union of canonical sample_ids
    -- it stores (AN, n_called) where AN = sum of ploidy over the group's members
    and n_called = number of members present on this chromosome (ploidy > 0).

    ploidy 0 (e.g. a female on chrY) => the sample is ABSENT and contributes NO
    alleles to AN and is not a called individual in any group (mirrors
    ``tier0_registry.PanelIndex``).

    Also exposes per-index membership so a carrier (identified by sample index on
    the axis) can be bucketed into its (db, subpop) and counted once globally.
    """

    __slots__ = ("subpop", "db", "global_", "_index_group", "_global_members")

    def __init__(self, axis, ploidy_of):
        # group_id -> [AN, n_called]
        self.subpop = {}
        self.db = {}
        self.global_ = [0, 0]
        # sample index -> (db, subpop, ploidy) for members present (ploidy>0)
        self._index_group = {}
        # canonical sample_id -> ploidy (deduped global membership)
        self._global_members = {}

        n = len(axis)
        for idx in range(n):
            sample_id = axis.sample_id(idx)
            database, subpopulation, sex = axis.sample_meta(idx)
            ploidy = ploidy_of(sample_id, sex)
            if ploidy <= 0:
                continue  # absent on this chromosome
            self._index_group[idx] = (database, subpopulation, ploidy)

            sp_key = db_subpop_group_id(database, subpopulation)
            b = self.subpop.get(sp_key)
            if b is None:
                b = [0, 0]
                self.subpop[sp_key] = b
            b[0] += ploidy
            b[1] += 1

            db_key = db_group_id(database)
            d = self.db.get(db_key)
            if d is None:
                d = [0, 0]
                self.db[db_key] = d
            d[0] += ploidy
            d[1] += 1

            # GLOBAL: dedup by canonical sample_id (same id in two dbs counts once;
            # ploidy taken from the first membership -- identical for one individual).
            if sample_id not in self._global_members:
                self._global_members[sample_id] = ploidy
                self.global_[0] += ploidy
                self.global_[1] += 1

    def index_group(self, idx):
        """(db, subpop, ploidy) for a sample index, or None if absent/ungrouped."""
        return self._index_group.get(idx)

    def is_global_member(self, sample_id):
        return sample_id in self._global_members


# --------------------------------------------------------------------------- #
# Unphased LOWER bound on combination carriers
# --------------------------------------------------------------------------- #
def _unphased_carrier_lower_bound(single_carrier_counts, n_called):
    """A guaranteed LOWER bound on the number of individuals carrying ALL k pairs.

    We know, per group, how many individuals carry EACH pair i alone
    (``single_carrier_counts[i]``) and the group size (``n_called``). The exact
    intersection is unknown without joint genotypes at the individual level, but
    two hard constraints bound it from below:

      (1) The intersection cannot exceed the SMALLEST single set:
              inter <= min_i(single_i).
          (This is the UPPER-bound witness; the assume-cis carrier count is a
          separate, tighter upper bound because it uses the JOINT genotypes.)

      (2) INCLUSION-EXCLUSION / Bonferroni floor. For sets A_1..A_k inside a
          universe of size N, the fraction NOT in a given set is
          (N - |A_i|)/N, and at most sum_i (N - |A_i|) individuals are missing
          from AT LEAST ONE set, so
              |A_1 ∩ ... ∩ A_k| >= N - sum_i (N - |A_i|)
                                 = sum_i|A_i| - (k-1)*N.
          Clamped at 0 (the bound is vacuous when the sets are small).

    The tightest guaranteed lower bound is thus:
        max(0, sum_i|A_i| - (k-1)*N)
    which also never exceeds min_i|A_i| (so it is a valid lower bound below the
    intersection). We return that value.

    NOTE: this is computed PER GROUP using that group's own single-pair carrier
    counts and its own n_called (never mixing groups).
    """
    if not single_carrier_counts:
        return 0
    k = len(single_carrier_counts)
    total = sum(single_carrier_counts)
    floor = total - (k - 1) * n_called
    if floor < 0:
        floor = 0
    # Also never exceed the smallest single set (defensive; the Bonferroni floor
    # already satisfies this, but clamp to be provably valid under rounding).
    smallest = min(single_carrier_counts)
    return min(floor, smallest)


# --------------------------------------------------------------------------- #
# k == 1 : exact single-variant summary straight from the Tier-0 registry.
# --------------------------------------------------------------------------- #
def _summarize_single(pos, alt, tier0_reader):
    """SINGLE-variant summary: exact per-group Counts from the Tier-0 registry.

    Returns a dict[group_id] -> GroupSummary. The registry already computed these
    full-panel and ploidy/sex-aware at build time; we simply wrap each Counts.
    """
    counts = tier0_reader.lookup(pos, alt)
    if counts is None:
        return {}
    groups = {}
    for gid, c in counts.items():
        groups[gid] = GroupSummary(
            allele_count=c.AC,
            allele_number=c.AN,
            n_carrier=c.n_carrier_indiv,
            n_hom=c.n_hom_indiv,
            n_called=c.n_called_indiv,
            n_carrier_lower=None,  # exact -- no bound needed
        )
    return groups


# --------------------------------------------------------------------------- #
# k >= 2 : combination summary from the Tier-1 genotype store + panel.
# --------------------------------------------------------------------------- #
def summarize(pos_alts, tier0_reader, gt_reader, axis, ploidy_of, phased,
              alt_index="1", panel=None):
    """Compute the combination-aware ``PopSummary`` for one off-target.

    Args:
      pos_alts   : list of (pos1based:int, alt:str) pairs that JOINTLY create the
                   off-target (k >= 1). For k==1 the single pair; for k>=2 the
                   full allele combination.
      tier0_reader : a ``tier0_registry.RegistryReader`` (used for k==1).
      gt_reader  : a ``tier1_genotypes.GenotypeReader`` (used for k>=2).
      axis       : a ``tier1_genotypes.SampleAxis`` (index<->id + per-index meta).
                   Typically ``gt_reader``'s own axis.
      ploidy_of  : callable(sample_id, sex) -> ploidy (0/1/2), ploidy/sex-aware.
      phased     : bool. Controls the k>=2 semantics (cis vs assume-cis bound).
      alt_index  : the genotype token denoting the record's alt (default "1", the
                   bcftools norm -m- biallelic convention).
      panel      : optional prebuilt ``Panel`` (denominators). Built from ``axis``
                   + ``ploidy_of`` if not supplied. Reuse across many off-targets.

    Returns a ``PopSummary``.
    """
    pos_alts = list(pos_alts)
    if len(pos_alts) < 1:
        raise ValueError("summarize() needs at least one (pos, alt) pair")

    # ---- k == 1 : exact, straight from the registry (full-panel, ploidy-aware). ---- #
    if len(pos_alts) == 1:
        pos, alt = pos_alts[0]
        groups = _summarize_single(pos, alt, tier0_reader)
        return PopSummary(pos_alts, phased=phased, groups=groups)

    # ---- k >= 2 : combination via the genotype tier + panel denominators. ---- #
    # GRACEFUL DEGRADATION (registry-only build): a k>=2 combination frequency is a
    # HAPLOTYPE (co-occurrence) property that CANNOT be computed without the Tier-1
    # genotype store (carriers) and the sample axis (denominators). When either is
    # absent (a --no-genotypes / registry-only install, gt_reader/axis None) we must
    # NOT build Panel(axis=None) -- that raises len(None) (TypeError) and the caller
    # would SKIP the whole off-target. Instead return an empty-groups PopSummary: the
    # None-group path renders it as a 'requires genotype tier' row (freqs NA for the
    # conservative unphased default, observed False) so the row is EMITTED, not
    # silently dropped -- the contract population_summary_companion already promises.
    if gt_reader is None or axis is None:
        return PopSummary(pos_alts, phased=phased, groups={})

    if panel is None:
        panel = Panel(axis, ploidy_of)

    k = len(pos_alts)

    # Per-pair carrier maps: sid -> gt string (alt-carriers only). A pair with a
    # missing record (no carriers on disk) has an EMPTY map -> zero combination
    # carriers (an off-target no one carries).
    pair_carriers = []
    for (pos, alt) in pos_alts:
        recon = gt_reader.carriers(pos, alt)
        pair_carriers.append(dict(recon) if recon else {})

    # The combination's candidate carriers are the INTERSECTION (by sample_id) of
    # all pairs' carrier sets: an individual must be an alt-carrier at EVERY pair.
    # Start from the smallest set for efficiency.
    smallest = min(range(k), key=lambda i: len(pair_carriers[i]))
    candidates = [sid for sid in pair_carriers[smallest]
                  if all(sid in pair_carriers[i] for i in range(k))]

    # Per-group accumulators. group_id -> dict with AC/haps, n_carrier, n_hom.
    acc = {}  # gid -> [AC_haps_or_None, n_carrier, n_hom]

    def _bump(gid, ac_haps, is_carrier, is_hom):
        a = acc.get(gid)
        if a is None:
            a = [0, 0, 0]
            acc[gid] = a
        if ac_haps is not None:
            a[0] += ac_haps
        a[1] += 1 if is_carrier else 0
        a[2] += 1 if is_hom else 0

    # Per-group single-pair carrier tallies (for the unphased lower bound). We
    # count, per group, how many of THAT group's members carry each pair alone.
    # gid -> list of length k of counts.
    single_tally = {}

    def _single_bump(gid, i):
        t = single_tally.get(gid)
        if t is None:
            t = [0] * k
            single_tally[gid] = t
        t[i] += 1

    # For the unphased lower bound we need per-group per-pair single-carrier
    # counts over the GROUP PANEL. Walk each pair's carriers once, bucket by group.
    if not phased:
        # Track globally-deduped single membership so GLOBAL single counts a
        # canonical id once per pair.
        global_single_seen = [set() for _ in range(k)]
        for i in range(k):
            for sid in pair_carriers[i]:
                idx = axis.get_index(sid)
                if idx is None:
                    continue
                grp = panel.index_group(idx)
                if grp is None:
                    continue
                database, subpopulation, _ploidy = grp
                _single_bump(db_subpop_group_id(database, subpopulation), i)
                _single_bump(db_group_id(database), i)
                if panel.is_global_member(sid) and sid not in global_single_seen[i]:
                    global_single_seen[i].add(sid)
                    _single_bump(GLOBAL_GROUP_ID, i)

    # Walk the (small) candidate set and decide carriage / homozygosity + haplotype
    # count, honoring phase. GLOBAL dedups the canonical sample_id.
    global_seen = set()
    for sid in candidates:
        idx = axis.get_index(sid)
        if idx is None:
            continue  # carrier not on the axis (defensive)
        grp = panel.index_group(idx)
        if grp is None:
            continue  # absent on this chromosome (ploidy 0)
        database, subpopulation, ploidy = grp

        gts = [pair_carriers[i][sid] for i in range(k)]

        if phased:
            ac_haps, is_carrier, is_hom = _phased_decide(gts, ploidy, alt_index)
        else:
            ac_haps, is_carrier, is_hom = _unphased_decide(gts, ploidy, alt_index)

        if not is_carrier:
            continue  # not a combination carrier -> contributes nothing

        sp_key = db_subpop_group_id(database, subpopulation)
        db_key = db_group_id(database)
        _bump(sp_key, ac_haps, is_carrier, is_hom)
        _bump(db_key, ac_haps, is_carrier, is_hom)
        if panel.is_global_member(sid) and sid not in global_seen:
            global_seen.add(sid)
            _bump(GLOBAL_GROUP_ID, ac_haps, is_carrier, is_hom)

    # Assemble SPARSE per-group GroupSummary (only groups with >= 1 carrier), each
    # with FULL-PANEL denominators from the panel.
    groups = {}
    for gid, (ac_haps, n_carrier, n_hom) in acc.items():
        if n_carrier < 1:
            continue
        AN, n_called = _panel_denoms(panel, gid)
        if phased:
            allele_count = ac_haps
            n_carrier_lower = None  # exact
        else:
            allele_count = None  # UNDEFINED without phase
            singles = single_tally.get(gid, [0] * k)
            n_carrier_lower = _unphased_carrier_lower_bound(singles, n_called)
            # The lower bound can never exceed the (assume-cis) upper bound we
            # actually counted; clamp defensively.
            if n_carrier_lower > n_carrier:
                n_carrier_lower = n_carrier
        groups[gid] = GroupSummary(
            allele_count=allele_count,
            allele_number=AN,
            n_carrier=n_carrier,
            n_hom=n_hom,
            n_called=n_called,
            n_carrier_lower=n_carrier_lower,
        )

    return PopSummary(pos_alts, phased=phased, groups=groups)


def _panel_denoms(panel, gid):
    """(AN, n_called) for a group id from the panel."""
    if gid == GLOBAL_GROUP_ID:
        return panel.global_[0], panel.global_[1]
    if SEP in gid:
        b = panel.subpop.get(gid)
        return (b[0], b[1]) if b else (0, 0)
    d = panel.db.get(gid)
    return (d[0], d[1]) if d else (0, 0)


def _phased_decide(gts, ploidy, alt_index="1"):
    """Phased combination decision for ONE individual across the k pairs.

    ``gts`` is the list of k genotype strings (this individual's gt at each pi).
    Returns (ac_haps, is_carrier, is_hom):
      * Parse each pair's gt into per-haplotype tokens (honoring ploidy). All
        pairs are on the SAME individual, so they share the same haplotype axis
        (hap 0, hap 1, ...). We require the SAME number of honored slots across
        pairs; if a gt is short we align on the min length (defensive; real
        phased panels have consistent ploidy per sample).
      * A haplotype h CARRIES the combination iff at EVERY pair the token on hap h
        equals the alt. ac_haps = number of such haplotypes.
      * is_carrier = ac_haps >= 1 (at least one cis haplotype).
      * is_hom = ac_haps == ploidy AND ploidy >= 1 (carries on ALL its
        haplotypes -- both haps of a diploid, or the single hap of a haploid).
    """
    hap_lists = [parse_haplotypes(gt, ploidy) for gt in gts]
    if not hap_lists:
        return 0, False, False
    n_hap = min(len(h) for h in hap_lists)
    if n_hap <= 0:
        return 0, False, False
    ac_haps = 0
    for h in range(n_hap):
        if all(_token_is_alt(hap_lists[i][h], alt_index) for i in range(len(gts))):
            ac_haps += 1
    is_carrier = ac_haps >= 1
    is_hom = (ac_haps == ploidy) and ploidy >= 1
    return ac_haps, is_carrier, is_hom


def _unphased_decide(gts, ploidy, alt_index="1"):
    """Unphased combination decision (ASSUME-CIS UPPER BOUND) for ONE individual.

    ``gts`` is the list of k genotype strings. Returns (ac_haps, is_carrier,
    is_hom):
      * ac_haps is None -- allele/haplotype counting is UNDEFINED without phase.
      * is_carrier = het-or-hom for the alt at EVERY pair (the individual COULD
        carry the combination in cis; this is the upper bound).
      * is_hom = hom for the alt at EVERY pair (an upper bound on true cis
        homozygosity: if hom at every site, both haplotypes carry every alt, so
        the combination IS present on both haps regardless of phase).
    """
    is_carrier = all(individual_carries_at(gt, ploidy, alt_index) for gt in gts)
    is_hom = is_carrier and all(
        individual_hom_at(gt, ploidy, alt_index) for gt in gts)
    return None, is_carrier, is_hom


# --------------------------------------------------------------------------- #
# Convenience: derive_record_stats-compatible dict for a PopSummary.
# --------------------------------------------------------------------------- #
def record_stats(summary):
    """Return a ``derive_record_stats``-shaped dict from a ``PopSummary``.

    For k==1 this equals ``derive_record_stats`` over the registry Counts (the
    exact path). For k>=2 it reports the combination-level equivalents, with
    ``global_af``/``max_subpop_af`` NaN/None-aware for the unphased case.
    """
    return {
        "global_af": summary.global_af,
        "global_carrier_freq": summary.global_carrier_freq,
        "global_hom_freq": summary.global_hom_freq,
        "global_carrier_n": summary.global_carrier_n,
        "max_subpop_af": summary.max_subpop_af,
        "max_subpop_af_label": summary.max_subpop_af_label,
        "observed": summary.observed,
        "allele_freq_defined": summary.allele_freq_defined,
    }


def is_nan(x):
    """True iff ``x`` is the NaN AF sentinel (undefined allele freq)."""
    return isinstance(x, float) and math.isnan(x)
