"""Tier-0 compact variant registry (CRISPRme+ dictless redesign, Phase 1).

This module is the correctness HEART of the "Tier 0" registry that will replace
the 152 GB per-sample SNP dicts. It has two halves:

  (A) A ploidy-aware, pure-function AGGREGATOR (``aggregate_record``) that turns
      one (position, specific-alt) record's per-sample genotypes into sparse
      per-group allele/carrier Counts. It honors per-sample ploidy so a haploid
      (chrX-nonPAR / chrY) male contributes AT MOST 1 to AC and exactly 1 to AN
      -- it NEVER duplicates a haploid GT to "treat it as homozygous" (the legacy
      bug in dictionary_creation_indels.py that doubles AC).

  (B) A little-endian, mmap-friendly BINARY FORMAT + COMPILER + READER. The reader
      does O(log n) binary search over a fixed-width sorted (pos, alt) key array
      via ``bisect`` on an ``mmap`` -- it NEVER parses the whole file into a Python
      dict/list.

Design notes / ground truth (see PostProcess/new_simple_analysis.py 158-290):
  * The existing SNP dict entry is per-ALT: "<samples>;<ref,alt>;<rsID>;<AF>" and
    ONLY lists alt-carriers. We keep the per-ALT record as first class here:
    multiallelic sites become multiple records at the same pos with distinct alt.
  * A genotype token that is "0" is ref, "." is missing, and any other token is
    THIS record's alt (dict entries are already split per-alt, so a non-"0"
    non-"." token means "carries the alt of this record").

STDLIB ONLY (struct, mmap, array, bisect, json) so it runs in the lightweight
unit-tests CI (no numpy / no pysam).
"""

from __future__ import annotations

import bisect
import json
import mmap
import struct

# --------------------------------------------------------------------------- #
# Group id conventions (also mirrored in the JSON manifest taxonomy).
# --------------------------------------------------------------------------- #
GLOBAL_GROUP_ID = "global"
SEP = "::"  # group_id separator, e.g. "1000G::EUR" ; "1000G" ; "global"


def db_group_id(database):
    """Group id for a whole database (aggregate of its subpopulations)."""
    return str(database)


def db_subpop_group_id(database, subpopulation):
    """Group id for one (database x subpopulation) cell."""
    return "{}{}{}".format(database, SEP, subpopulation)


# --------------------------------------------------------------------------- #
# (A) Aggregator core
# --------------------------------------------------------------------------- #
class Counts(object):
    """Per-group allele/carrier counts for ONE (pos, alt) record.

    All five fields are stored as INTEGERS (persisted as u32); AF and the freq
    helpers are DERIVED, never stored.

      AC             number of alt alleles actually called (missing excluded)
      AN             number of called alleles == sum of per-sample called ploidy
                     (missing allele slots excluded) -> AF = AC / AN, NOT AC / 2N
      n_carrier_indiv  individuals with >= 1 alt allele
      n_hom_indiv    individuals whose ALL called alleles are alt (a haploid male
                     with 1 alt copy is hom/hemizygous-carrier)
      n_called_indiv individuals with >= 1 called (non-missing) allele
    """

    __slots__ = ("AC", "AN", "n_carrier_indiv", "n_hom_indiv", "n_called_indiv")

    def __init__(self, AC=0, AN=0, n_carrier_indiv=0, n_hom_indiv=0,
                 n_called_indiv=0):
        self.AC = AC
        self.AN = AN
        self.n_carrier_indiv = n_carrier_indiv
        self.n_hom_indiv = n_hom_indiv
        self.n_called_indiv = n_called_indiv

    def add(self, other):
        self.AC += other.AC
        self.AN += other.AN
        self.n_carrier_indiv += other.n_carrier_indiv
        self.n_hom_indiv += other.n_hom_indiv
        self.n_called_indiv += other.n_called_indiv

    def as_tuple(self):
        return (self.AC, self.AN, self.n_carrier_indiv, self.n_hom_indiv,
                self.n_called_indiv)

    def __eq__(self, other):
        return isinstance(other, Counts) and self.as_tuple() == other.as_tuple()

    def __hash__(self):
        return hash(self.as_tuple())

    def __repr__(self):
        return ("Counts(AC={}, AN={}, n_carrier_indiv={}, n_hom_indiv={}, "
                "n_called_indiv={})").format(*self.as_tuple())

    # ---- derived helpers (never stored) ---- #
    def allele_freq(self):
        return (self.AC / self.AN) if self.AN else 0.0

    def carrier_freq(self):
        return (self.n_carrier_indiv / self.n_called_indiv) if self.n_called_indiv else 0.0

    def hom_freq(self):
        return (self.n_hom_indiv / self.n_called_indiv) if self.n_called_indiv else 0.0


def _parse_genotype(gt):
    """Split a genotype string into its allele tokens.

    Accepts "a0|a1" (phased), "a0/a1" (unphased) and "a0" (haploid). Returns the
    list of raw allele tokens exactly as written (no ploidy assumption here; the
    caller-supplied ``ploidy_of`` decides how many slots count).
    """
    if gt is None:
        return []
    gt = gt.strip()
    if gt == "":
        return []
    # phased '|' or unphased '/'
    if "|" in gt:
        return gt.split("|")
    if "/" in gt:
        return gt.split("/")
    return [gt]


def _sample_alleles(gt, ploidy):
    """Return exactly ``ploidy`` allele tokens for a sample.

    This is where the anti-double-count invariant lives:
      * ploidy 1 (haploid male on chrX-nonPAR / chrY): take a SINGLE token. We do
        NOT re-count it a second time. If the GT was written diploid-style
        ("1|1") for a region the caller has declared haploid, only the first slot
        is honored -> AC/AN contribution is at most 1.
      * ploidy 2: take up to two tokens. A truly haploid GT ("1") declared diploid
        yields ONE slot (the other is simply absent, i.e. treated as an
        uncalled/short GT) rather than being duplicated.

    We never pad a short GT by duplicating a token.
    """
    toks = _parse_genotype(gt)
    if ploidy <= 1:
        return toks[:1]
    return toks[:ploidy]


def _classify(token, alt_index="1"):
    """Classify one allele token against THIS record's SPECIFIC alt.

    Returns 'alt', 'ref', or 'missing'. Under the default pipeline
    (``bcftools norm -m-``, see convert_gnomAD_vcfs.py) multiallelic sites are
    split into biallelic records whose carried allele is recoded to '1', so
    ``alt_index`` defaults to "1". A token equal to ``alt_index`` is THIS
    record's alt; "." is missing; ANY OTHER token -- ref "0" OR a foreign alt
    index like "2" on an unsplit multiallelic genotype -- is NOT this record's
    alt and is treated as 'ref' (a called, non-carrier allele). This is what
    prevents a '1|2' genotype from being miscounted as two copies of alt-1.
    """
    if token == ".":
        return "missing"
    if token == alt_index:
        return "alt"
    return "ref"


def _sample_counts(gt, ploidy, alt_index="1"):
    """Compute (ac, an, is_carrier, is_hom, is_called) for one sample.

    ac  : alt alleles called for this sample (0..ploidy)
    an  : called (non-missing) alleles for this sample (0..ploidy)
    is_carrier : an alt allele was called
    is_hom     : sample is FULLY called (an == ploidy) AND every allele is alt.
                 A half-missing diploid ('1|.') is NOT hom -- the unobserved
                 allele may be ref -- matching bcftools/plink. A haploid male
                 with a single alt copy IS hom/hemizygous (an == ploidy == 1).
    is_called  : >=1 non-missing allele
    """
    alleles = _sample_alleles(gt, ploidy)
    ac = 0
    an = 0
    for tok in alleles:
        kind = _classify(tok, alt_index)
        if kind == "missing":
            continue
        an += 1
        if kind == "alt":
            ac += 1
    is_called = an > 0
    is_carrier = ac > 0
    is_hom = is_called and ac == an and an == ploidy  # fully-called, all alt
    return ac, an, is_carrier, is_hom, is_called


def aggregate_record(alt_genotypes, sample_meta, ploidy_of, alt_index="1"):
    """Aggregate ONE (pos, specific-alt) record into sparse per-group Counts.

    Args:
      alt_genotypes: dict sample_id -> genotype string ("a0|a1"/"a0/a1"/"a0").
        Tokens in {"0","1",...,"."}. Only samples relevant to this record need
        appear; but any sample_id present in ``sample_meta`` and NOT in
        ``alt_genotypes`` simply does not contribute (matches the legacy dict,
        which lists only alt-carriers). Samples in ``alt_genotypes`` but not in
        ``sample_meta`` are ignored (no taxonomy => cannot be grouped).
      sample_meta: sample_id -> (database, subpopulation, sex).
      ploidy_of: callable(sample_id, sex) -> 1 or 2.
      alt_index: the genotype token that denotes THIS record's alt (default "1",
        the bcftools norm -m- biallelic convention). See ``_classify``.

    Returns:
      dict[group_id] -> Counts. SPARSE: a (db x subpop) or db group is emitted
      only if it has >= 1 carrier individual. The GLOBAL group is emitted only if
      the deduplicated union has >= 1 carrier.

    Note on GLOBAL dedup: ``alt_genotypes`` is keyed by sample_id, so within a
    single call each sample contributes once by construction. Deduplicating the
    SAME canonical individual that appears across MULTIPLE databases is the
    CALLER's (merge-layer's) responsibility: it must union the databases into one
    ``alt_genotypes``/``sample_meta`` under a canonical sample_id before calling
    here (the build-time global sample axis, per the design). This function then
    counts each canonical id exactly once globally.
    """
    # Per (db, subpop) accumulation.
    subpop_counts = {}   # (db, subpop) -> Counts
    subpop_carrier = {}  # (db, subpop) -> bool (>=1 carrier seen)
    db_counts = {}       # db -> Counts
    db_carrier = {}      # db -> bool

    # Global dedup: canonical sample_id counted once. We record its per-sample
    # contribution the FIRST time we see it (a sample with the same canonical id
    # in two databases carries identical genotype semantics for the site).
    global_seen = {}     # sample_id -> (ac, an, is_carrier, is_hom, is_called)

    for sample_id, gt in alt_genotypes.items():
        meta = sample_meta.get(sample_id)
        if meta is None:
            continue  # ungroupable
        database, subpopulation, sex = meta
        ploidy = ploidy_of(sample_id, sex)
        ac, an, is_carrier, is_hom, is_called = _sample_counts(gt, ploidy, alt_index)

        # A sample with no called alleles (fully missing) still does not create a
        # group, but it DOES count toward n_called_indiv only if called. Since a
        # fully-missing sample has an==0 and is_called False, it contributes
        # nothing; but we must remember: groups are keyed by carriers, and
        # n_called_indiv counts called individuals within any group that exists.
        cnt = Counts(
            AC=ac,
            AN=an,
            n_carrier_indiv=1 if is_carrier else 0,
            n_hom_indiv=1 if is_hom else 0,
            n_called_indiv=1 if is_called else 0,
        )

        key = (database, subpopulation)
        if key not in subpop_counts:
            subpop_counts[key] = Counts()
            subpop_carrier[key] = False
        subpop_counts[key].add(cnt)
        if is_carrier:
            subpop_carrier[key] = True

        if database not in db_counts:
            db_counts[database] = Counts()
            db_carrier[database] = False
        db_counts[database].add(cnt)
        if is_carrier:
            db_carrier[database] = True

        if sample_id not in global_seen:
            global_seen[sample_id] = (ac, an, is_carrier, is_hom, is_called)

    result = {}

    # Sparse (db x subpop) groups: only those with >=1 carrier.
    for key, cnt in subpop_counts.items():
        if subpop_carrier[key]:
            result[db_subpop_group_id(key[0], key[1])] = cnt

    # Sparse db groups: only those with >=1 carrier.
    for database, cnt in db_counts.items():
        if db_carrier[database]:
            result[db_group_id(database)] = cnt

    # Global group over deduplicated union.
    g = Counts()
    g_carrier = False
    for sample_id, (ac, an, is_carrier, is_hom, is_called) in global_seen.items():
        g.AC += ac
        g.AN += an
        g.n_carrier_indiv += 1 if is_carrier else 0
        g.n_hom_indiv += 1 if is_hom else 0
        g.n_called_indiv += 1 if is_called else 0
        if is_carrier:
            g_carrier = True
    if g_carrier:
        result[GLOBAL_GROUP_ID] = g

    return result


# --------------------------------------------------------------------------- #
# (A') Panel-aware aggregation (FULL-PANEL AN via hom-ref baseline + delta)
# --------------------------------------------------------------------------- #
#
# The legacy SNP dict lists ONLY alt-carriers per record. To re-derive a correct
# population allele frequency AF = AC / AN, every panel sample that is NOT listed
# for a record must be counted as hom-ref: it contributes AC 0, AN == its ploidy,
# is_called (AN slots all observed as ref), NOT a carrier, NOT hom. So AN is the
# FULL-PANEL called-allele count, not the alleles among carriers only.
#
# Naively this is O(#panel) per record (~60M records x ~4400 samples). Instead we
# precompute ONCE, per group, the baseline totals assuming EVERY member is
# hom-ref (PanelIndex). For a record we COPY that baseline and, for each of the
# (few) carriers, add the DELTA between the carrier's TRUE per-sample counts and
# the hom-ref assumption the baseline already made for it (ac 0, an ploidy,
# called 1, carrier 0, hom 0). Cost is O(#carriers) per record.
#
# Sparsity: a group is WRITTEN only if it has >= 1 carrier for the record (a group
# with 0 carriers carries no off-target signal), but the counts we emit for a
# written group are the FULL-PANEL counts (so AF/carrier_freq denominators are
# right), not the carriers-only counts.


class PanelIndex(object):
    """Precomputed per-group hom-ref baselines for a fixed panel.

    Built ONCE per compile. For each group -- every (db x subpop) cell, every db,
    and the GLOBAL group over the DEDUPED union of canonical sample ids -- it
    stores the baseline Counts assuming ALL members are hom-ref:

        AC = 0
        AN = sum(ploidy over the group's members)
        n_carrier_indiv = 0
        n_hom_indiv = 0
        n_called_indiv = number of individuals in the group

    It also records, per sample, the (db, subpop) it belongs to and its ploidy, so
    ``aggregate_record_panel`` can apply carrier deltas without re-deriving group
    membership. The GLOBAL axis dedups by canonical sample_id: a sample present in
    two databases is counted once globally (its ploidy taken from the first
    membership seen; identical across databases for the same individual).
    """

    __slots__ = ("subpop_baseline", "db_baseline", "global_baseline",
                 "_sample_group", "_sample_ploidy_global", "n_global")

    def __init__(self, sample_meta, ploidy_of):
        # Per-group baseline Counts (hom-ref for every member).
        self.subpop_baseline = {}   # (db, subpop) -> Counts
        self.db_baseline = {}       # db -> Counts
        self.global_baseline = Counts()

        # Per-sample resolved (db, subpop, ploidy) for delta application.
        self._sample_group = {}     # sample_id -> (db, subpop)
        # Global ploidy per DEDUPED canonical id (counted once globally).
        self._sample_ploidy_global = {}  # sample_id -> ploidy
        self.n_global = 0

        for sample_id, (database, subpopulation, sex) in sample_meta.items():
            ploidy = ploidy_of(sample_id, sex)
            # ploidy 0 => the sample is ABSENT on this chromosome (e.g. a female on
            # chrY): exclude it from the panel entirely so it contributes NO phantom
            # alleles to AN and is not counted as a called individual in any group.
            if ploidy <= 0:
                continue
            self._sample_group[sample_id] = (database, subpopulation)

            key = (database, subpopulation)
            b = self.subpop_baseline.get(key)
            if b is None:
                b = Counts()
                self.subpop_baseline[key] = b
            b.AN += ploidy
            b.n_called_indiv += 1

            d = self.db_baseline.get(database)
            if d is None:
                d = Counts()
                self.db_baseline[database] = d
            d.AN += ploidy
            d.n_called_indiv += 1

            # GLOBAL: dedup by canonical sample_id.
            if sample_id not in self._sample_ploidy_global:
                self._sample_ploidy_global[sample_id] = ploidy
                self.global_baseline.AN += ploidy
                self.global_baseline.n_called_indiv += 1
                self.n_global += 1

    def sample_group(self, sample_id):
        """(db, subpop) for a sample, or None if not in the panel."""
        return self._sample_group.get(sample_id)

    def is_global_member(self, sample_id):
        return sample_id in self._sample_ploidy_global


def _apply_delta(dst, ac, an, is_carrier, is_hom, is_called):
    """Add the (true - hom-ref-assumed) per-sample delta into a group's Counts.

    The baseline already assumed for this sample: ac 0, an ploidy, called 1,
    carrier 0, hom 0. ``an`` here is the sample's TRUE called-allele count, which
    for a listed carrier equals its ploidy MINUS any missing slots. The baseline
    assumed ``ploidy`` called alleles, but a carrier is by definition a listed
    sample whose true AN is its observed non-missing slots; the AN delta is
    (true_an - ploidy). We do NOT have ``ploidy`` here directly, so the caller
    passes the DELTA-ready values: see ``aggregate_record_panel`` which computes
    (true - assumed) explicitly.
    """
    dst.AC += ac
    dst.AN += an
    dst.n_carrier_indiv += is_carrier
    dst.n_hom_indiv += is_hom
    dst.n_called_indiv += is_called


def aggregate_record_panel(carrier_genotypes, panel_index, sample_meta,
                           ploidy_of, alt_index="1"):
    """Aggregate ONE (pos, specific-alt) record with FULL-PANEL AN.

    Args:
      carrier_genotypes: dict sample_id -> genotype string, listing ONLY the
        alt-carriers for this record (exactly what the legacy dict stores). Any
        panel sample NOT listed is treated as hom-ref (via the baseline).
      panel_index: a ``PanelIndex`` built from the SAME sample_meta/ploidy_of.
      sample_meta: sample_id -> (database, subpopulation, sex).
      ploidy_of: callable(sample_id, sex) -> 1 or 2 (same as the PanelIndex).
      alt_index: genotype token denoting THIS record's alt (default "1").

    Returns:
      dict[group_id] -> Counts. SPARSE on "has a carrier": a group is emitted only
      if >= 1 of its members carries the alt for THIS record. BUT the emitted
      counts are FULL-PANEL (baseline + carrier deltas), so AC/AN/n_called reflect
      the WHOLE group, giving correct AF = AC / AN denominators.

    Cost: O(#carriers). We copy only the baselines of the (few) touched groups.
    """
    # Touched groups (those with >= 1 carrier). We lazily copy each touched
    # group's baseline the first time a carrier hits it, then apply deltas.
    subpop_out = {}   # (db, subpop) -> Counts (copy of baseline + deltas)
    db_out = {}       # db -> Counts
    global_out = None  # Counts (copy of global baseline) once a global carrier hits

    # Global dedup: apply each canonical id's delta once.
    global_seen = set()

    for sample_id, gt in carrier_genotypes.items():
        meta = sample_meta.get(sample_id)
        if meta is None:
            continue  # ungroupable / not in panel
        grp = panel_index.sample_group(sample_id)
        if grp is None:
            continue  # in meta but not in the panel baseline (defensive)
        database, subpopulation = grp
        _, _, sex = meta
        ploidy = ploidy_of(sample_id, sex)

        true_ac, true_an, is_carrier, is_hom, is_called = _sample_counts(
            gt, ploidy, alt_index)

        # Delta vs the hom-ref assumption the baseline already made:
        #   assumed: ac=0, an=ploidy, called=1, carrier=0, hom=0
        d_ac = true_ac - 0
        d_an = true_an - ploidy
        d_car = (1 if is_carrier else 0) - 0
        d_hom = (1 if is_hom else 0) - 0
        d_call = (1 if is_called else 0) - 1

        # (db x subpop) group
        key = (database, subpopulation)
        out = subpop_out.get(key)
        if out is None:
            base = panel_index.subpop_baseline.get(key)
            # base must exist: the sample was registered in the PanelIndex.
            out = Counts(*base.as_tuple())
            subpop_out[key] = out
        _apply_delta(out, d_ac, d_an, d_car, d_hom, d_call)

        # db group
        dout = db_out.get(database)
        if dout is None:
            dbase = panel_index.db_baseline.get(database)
            dout = Counts(*dbase.as_tuple())
            db_out[database] = dout
        _apply_delta(dout, d_ac, d_an, d_car, d_hom, d_call)

        # GLOBAL group (dedup by canonical id)
        if sample_id not in global_seen and panel_index.is_global_member(sample_id):
            global_seen.add(sample_id)
            if global_out is None:
                global_out = Counts(*panel_index.global_baseline.as_tuple())
            _apply_delta(global_out, d_ac, d_an, d_car, d_hom, d_call)

    result = {}
    # Only WRITE groups that have >= 1 carrier for this record; but their counts
    # are already full-panel (baseline + deltas).
    for key, cnt in subpop_out.items():
        if cnt.n_carrier_indiv >= 1:
            result[db_subpop_group_id(key[0], key[1])] = cnt
    for database, cnt in db_out.items():
        if cnt.n_carrier_indiv >= 1:
            result[db_group_id(database)] = cnt
    if global_out is not None and global_out.n_carrier_indiv >= 1:
        result[GLOBAL_GROUP_ID] = global_out
    return result


def derive_record_stats(groups):
    """Derive record-level summary stats from a record's per-group Counts.

    Args:
      groups: dict[group_id] -> Counts (the output of ``aggregate_record`` or a
        reader ``lookup``).

    Returns dict with:
      global_af, global_carrier_freq, global_hom_freq : from the GLOBAL group
      max_subpop_af : max allele_freq over (db x subpop) groups (0.0 if none)
      max_subpop_af_label : the group_id achieving that max (None if none)
      observed : True iff the global group exists with n_carrier_indiv >= 1
    """
    g = groups.get(GLOBAL_GROUP_ID)
    global_af = g.allele_freq() if g else 0.0
    global_carrier_freq = g.carrier_freq() if g else 0.0
    global_hom_freq = g.hom_freq() if g else 0.0
    observed = bool(g and g.n_carrier_indiv >= 1)

    max_af = 0.0
    max_label = None
    # Iterate in sorted group-id order so ties resolve deterministically to the
    # alphabetically-first (db x subpop) label (any argmax is valid, but the
    # reported label must be stable across insertion orders).
    for gid in sorted(groups):
        if gid == GLOBAL_GROUP_ID:
            continue
        # (db x subpop) groups contain the SEP; db groups do not. Only consider
        # (db x subpop) groups for max-subpop.
        if SEP not in gid:
            continue
        af = groups[gid].allele_freq()
        if max_label is None or af > max_af:
            max_af = af
            max_label = gid
    return {
        "global_af": global_af,
        "global_carrier_freq": global_carrier_freq,
        "global_hom_freq": global_hom_freq,
        "max_subpop_af": max_af,
        "max_subpop_af_label": max_label,
        "observed": observed,
    }


# --------------------------------------------------------------------------- #
# (B) Binary format + compiler + mmap reader
# --------------------------------------------------------------------------- #
#
# File layout (all little-endian):
#
#   <out_bin>  primary binary
#     [1] MAGIC  b"T0RG"  (4 bytes)
#     [2] VERSION u32 (=1)
#     [3] n_records u32
#     [4] alt_field_width u8  (== 1 for SNPs; alt strings padded/truncated)
#                             (kept generic so callers see the constant)
#     [5] pad (3 bytes, zero) -> 16-byte header
#     [6] RECORD ARRAY : n_records * RECORD_STRUCT, sorted by (pos, alt_bytes).
#           each record (fixed width, MMAP-FRIENDLY, no Python parsing needed to
#           binary-search):
#             q   pos            (int64)
#             c   alt            (1 byte, ascii)
#             c   ref            (1 byte, ascii)
#             xx  pad (2 bytes)
#             I   rsid_off       (u32, byte offset into STRING POOL; the pool
#                                 stores a u16 length prefix then utf-8 bytes)
#             I   n_groups       (u32, number of group entries for this record)
#             Q   groups_off     (u64, byte offset into GROUP BLOB)
#         => struct format "<qccxxIIQ" = 8+1+1+2+4+4+8 = 28 bytes.
#     [7] GROUP BLOB : concatenated group entries. Each entry:
#             I   gid_off        (u32 offset into STRING POOL for the group_id)
#             I   AC
#             I   AN
#             I   n_carrier_indiv
#             I   n_hom_indiv
#             I   n_called_indiv
#         => struct format "<IIIIII" = 24 bytes. A record's n_groups entries are
#            stored contiguously starting at groups_off.
#     [8] STRING POOL : deduplicated strings. Each string = u16 length + utf-8.
#
#   <out_idx>  json sidecar / manifest, capturing section offsets + the taxonomy
#     so the reader can mmap sections without re-deriving them.
#
# The record array is a fixed-width sorted array => bisect on an mmap gives the
# record without materializing the array. The idx also stores the record-array
# byte offset + record stride so the reader can seek directly.

MAGIC = b"T0RG"
VERSION = 1
_HEADER_STRUCT = struct.Struct("<4sIIBxxx")   # magic, version, n_records, altw, pad
_RECORD_STRUCT = struct.Struct("<qccxxIIQ")   # pos, alt, ref, pad, rsid_off, n_groups, groups_off
_GROUP_STRUCT = struct.Struct("<IIIIII")      # gid_off, AC, AN, ncar, nhom, ncall
HEADER_SIZE = _HEADER_STRUCT.size             # 16
RECORD_SIZE = _RECORD_STRUCT.size             # 28
GROUP_SIZE = _GROUP_STRUCT.size               # 24


class _StringPool(object):
    """Builds a deduplicated string pool; returns byte offsets."""

    def __init__(self):
        self._buf = bytearray()
        self._offsets = {}

    def intern(self, s):
        if s is None:
            s = ""
        if s in self._offsets:
            return self._offsets[s]
        off = len(self._buf)
        data = s.encode("utf-8")
        if len(data) > 0xFFFF:
            raise ValueError("string too long for u16 length prefix: %r" % (s[:32],))
        self._buf += struct.pack("<H", len(data))
        self._buf += data
        self._offsets[s] = off
        return off

    def getbuffer(self):
        return bytes(self._buf)


def _pack_alt(alt):
    """Pack an alt string into the fixed 1-byte record slot.

    SNP registries use single-character alts. We enforce that here; the format
    could grow ``alt_field_width`` but Phase 1 is SNP-first (per-alt records).
    """
    if not alt or len(alt) != 1:
        raise ValueError(
            "tier0_registry Phase-1 supports single-char alt only (per-alt SNP "
            "records); got alt=%r" % (alt,)
        )
    return alt.encode("ascii")


def _pack_ref(ref):
    r = (ref or ".")[:1]
    return r.encode("ascii")


def compile_registry(records, sample_meta, taxonomy, ploidy_of, out_bin, out_idx,
                     alt_index="1"):
    """Compile records into the mmap-friendly binary + json manifest.

    CARRIERS-ONLY aggregation (``aggregate_record``): AN counts only alleles among
    the samples actually listed for each record. For full-panel AN (correct AF
    denominators over the WHOLE panel) use ``compile_registry_panel``.

    Args:
      records: iterable of (pos:int, ref:str1, alt:str1, rsid:str,
               alt_genotypes:dict). Multiallelic => multiple records at the same
               pos with different alt. Need not be pre-sorted; we sort by
               (pos, alt).
      sample_meta: sample_id -> (database, subpopulation, sex).
      taxonomy: dict describing databases -> {sample_count, phased_placeholder,
                subpopulations:[...]}, or None to auto-derive from sample_meta.
      ploidy_of: callable(sample_id, sex) -> 1 or 2.
      out_bin: path for the binary.
      out_idx: path for the json manifest.

    Returns the manifest dict (also written to out_idx).
    """
    recs = list(records)
    recs.sort(key=lambda r: (int(r[0]), r[2]))

    def agg(alt_genotypes, ai):
        return aggregate_record(alt_genotypes, sample_meta, ploidy_of, ai)

    return _write_registry(recs, agg, sample_meta, taxonomy, out_bin, out_idx,
                           alt_index=alt_index, aggregation="carriers")


def compile_registry_panel(records, sample_meta, taxonomy, ploidy_of, out_bin,
                           out_idx, alt_index="1", panel_index=None):
    """Compile records with FULL-PANEL AN (``aggregate_record_panel``).

    Identical binary format + manifest to ``compile_registry``, but every emitted
    group's AC/AN/n_called reflect the WHOLE panel of that group (carriers plus
    all unlisted samples counted as hom-ref), so AF = AC / AN uses the full-panel
    denominator. Groups remain SPARSE on "has a carrier".

    Args are the same as ``compile_registry``; ``records`` here carry
    CARRIER-ONLY genotype dicts (the legacy dict semantics). A ``PanelIndex`` is
    built once from ``sample_meta``/``ploidy_of`` (or reuse a prebuilt one via
    ``panel_index``).
    """
    recs = list(records)
    recs.sort(key=lambda r: (int(r[0]), r[2]))

    if panel_index is None:
        panel_index = PanelIndex(sample_meta, ploidy_of)

    def agg(carrier_genotypes, ai):
        return aggregate_record_panel(carrier_genotypes, panel_index, sample_meta,
                                      ploidy_of, ai)

    return _write_registry(recs, agg, sample_meta, taxonomy, out_bin, out_idx,
                           alt_index=alt_index, aggregation="panel")


def _write_registry(recs, aggregate_fn, sample_meta, taxonomy, out_bin, out_idx,
                    alt_index="1", aggregation="carriers"):
    """Shared binary writer. ``aggregate_fn(alt_genotypes, alt_index)`` returns
    dict[group_id]->Counts for one (pos, alt) record (either the carriers-only
    ``aggregate_record`` path or the full-panel ``aggregate_record_panel`` path).
    """
    pool = _StringPool()
    group_blob = bytearray()
    record_rows = []  # (pos, alt_byte, ref_byte, rsid_off, n_groups, groups_off)

    # Track which group_ids actually appear (for the manifest taxonomy).
    seen_group_ids = {}  # group_id -> (kind, database, subpopulation)

    for (pos, ref, alt, rsid, alt_genotypes) in recs:
        groups = aggregate_fn(alt_genotypes, alt_index)
        # Stable, deterministic group ordering: global last, then db, then
        # (db x subpop) sorted -- but ordering does not affect lookup semantics.
        gids = sorted(groups.keys())
        groups_off = len(group_blob)
        for gid in gids:
            cnt = groups[gid]
            gid_off = pool.intern(gid)
            group_blob += _GROUP_STRUCT.pack(
                gid_off, cnt.AC, cnt.AN,
                cnt.n_carrier_indiv, cnt.n_hom_indiv, cnt.n_called_indiv,
            )
            if gid not in seen_group_ids:
                if gid == GLOBAL_GROUP_ID:
                    seen_group_ids[gid] = ("global", None, None)
                elif SEP in gid:
                    db, sp = gid.split(SEP, 1)
                    seen_group_ids[gid] = ("db_subpop", db, sp)
                else:
                    seen_group_ids[gid] = ("db", gid, None)
        rsid_off = pool.intern(rsid if rsid else ".")
        record_rows.append(
            (int(pos), _pack_alt(alt), _pack_ref(ref), rsid_off, len(gids), groups_off)
        )

    n_records = len(record_rows)

    # Compute section offsets.
    record_array_off = HEADER_SIZE
    group_blob_off = record_array_off + n_records * RECORD_SIZE
    string_pool_off = group_blob_off + len(group_blob)

    # The stored groups_off / rsid_off are RELATIVE to their section start; make
    # them absolute file offsets so the reader can seek with a single addend.
    with open(out_bin, "wb") as fh:
        fh.write(_HEADER_STRUCT.pack(MAGIC, VERSION, n_records, 1))
        for (pos, altb, refb, rsid_off, n_groups, groups_off) in record_rows:
            fh.write(_RECORD_STRUCT.pack(
                pos, altb, refb,
                string_pool_off + rsid_off,
                n_groups,
                group_blob_off + groups_off,
            ))
        # rewrite gid_off references in the blob to absolute pool offsets
        # (they were interned as pool-relative). Do it in a copy.
        blob = bytearray(group_blob)
        for i in range(0, len(blob), GROUP_SIZE):
            gid_off, AC, AN, ncar, nhom, ncall = _GROUP_STRUCT.unpack_from(blob, i)
            _GROUP_STRUCT.pack_into(blob, i, string_pool_off + gid_off,
                                    AC, AN, ncar, nhom, ncall)
        fh.write(bytes(blob))
        fh.write(pool.getbuffer())

    # ---- manifest / taxonomy ---- #
    if taxonomy is None:
        taxonomy = _derive_taxonomy(sample_meta)

    group_taxonomy = {}
    for gid, (kind, db, sp) in seen_group_ids.items():
        group_taxonomy[gid] = {"kind": kind, "database": db, "subpopulation": sp}

    manifest = {
        "magic": MAGIC.decode("ascii"),
        "version": VERSION,
        "n_records": n_records,
        "record_array_off": record_array_off,
        "record_size": RECORD_SIZE,
        "group_blob_off": group_blob_off,
        "group_size": GROUP_SIZE,
        "string_pool_off": string_pool_off,
        "alt_field_width": 1,
        "aggregation": aggregation,  # "carriers" (AN over listed only) or
                                     # "panel" (AN over the full panel)
        "global_group_id": GLOBAL_GROUP_ID,
        "group_sep": SEP,
        "databases": taxonomy,
        "groups": group_taxonomy,
    }
    with open(out_idx, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)

    return manifest


def _derive_taxonomy(sample_meta):
    """Auto-derive per-database {sample_count, phased_placeholder, subpops}."""
    dbs = {}
    for sample_id, (database, subpopulation, sex) in sample_meta.items():
        d = dbs.setdefault(database, {
            "sample_count": 0,
            "phased_placeholder": True,
            "subpopulations": set(),
        })
        d["sample_count"] += 1
        d["subpopulations"].add(subpopulation)
    for d in dbs.values():
        d["subpopulations"] = sorted(d["subpopulations"])
    return dbs


class RegistryReader(object):
    """mmap-backed reader. Binary-searches a sorted fixed-width record array.

    Does NOT load the record array into a Python list/dict; ``lookup`` uses
    ``bisect`` over an on-mmap key accessor.
    """

    def __init__(self, bin_path, idx_path):
        with open(idx_path, "r") as fh:
            self.manifest = json.load(fh)
        self._n = self.manifest["n_records"]
        self._rec_off = self.manifest["record_array_off"]
        self._rec_size = self.manifest["record_size"]
        self._fh = open(bin_path, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        # sanity-check header
        magic, version, n_records, altw = _HEADER_STRUCT.unpack_from(self._mm, 0)
        if magic != MAGIC:
            raise ValueError("bad magic in %s: %r" % (bin_path, magic))
        if n_records != self._n:
            raise ValueError("n_records mismatch between bin and idx")

        # bisect helper: a virtual sorted sequence of record keys, computed on
        # demand from the mmap so we never build a Python list of all records.
        reader = self

        class _KeyView(object):
            def __len__(_self):
                return reader._n

            def __getitem__(_self, i):
                return reader._key_at(i)

        self._keys = _KeyView()

    # ---- low-level record access (no full parse) ---- #
    def _record_base(self, i):
        return self._rec_off + i * self._rec_size

    def _key_at(self, i):
        base = self._record_base(i)
        pos = struct.unpack_from("<q", self._mm, base)[0]
        altb = self._mm[base + 8:base + 9]  # single byte after the int64 pos
        return (pos, altb)

    def _record_at(self, i):
        base = self._record_base(i)
        (pos, altb, refb, rsid_off, n_groups, groups_off) = \
            _RECORD_STRUCT.unpack_from(self._mm, base)
        return pos, altb, refb, rsid_off, n_groups, groups_off

    def _read_string(self, off):
        (n,) = struct.unpack_from("<H", self._mm, off)
        return self._mm[off + 2:off + 2 + n].decode("utf-8")

    def _find_index(self, pos, alt):
        key = (int(pos), alt.encode("ascii"))
        i = bisect.bisect_left(self._keys, key)
        if i < self._n and self._key_at(i) == key:
            return i
        return None

    # ---- public API ---- #
    def __len__(self):
        return self._n

    def lookup(self, pos, alt):
        """Return dict[group_id]->Counts for (pos, alt), or None if absent."""
        i = self._find_index(pos, alt)
        if i is None:
            return None
        _, _, _, _, n_groups, groups_off = self._record_at(i)
        groups = {}
        for k in range(n_groups):
            base = groups_off + k * self.manifest["group_size"]
            gid_off, AC, AN, ncar, nhom, ncall = _GROUP_STRUCT.unpack_from(self._mm, base)
            gid = self._read_string(gid_off)
            groups[gid] = Counts(AC, AN, ncar, nhom, ncall)
        return groups

    def rsid(self, pos, alt):
        """Return the rsID string for (pos, alt), or None if the record absent."""
        i = self._find_index(pos, alt)
        if i is None:
            return None
        _, _, _, rsid_off, _, _ = self._record_at(i)
        return self._read_string(rsid_off)

    def record_key_at(self, i):
        """(pos, alt_str) at index i -- for iteration / first/last tests."""
        pos, altb = self._key_at(i)
        return pos, altb.decode("ascii")

    def derived(self, pos, alt):
        """Return derive_record_stats for a looked-up record (None if absent)."""
        groups = self.lookup(pos, alt)
        if groups is None:
            return None
        return derive_record_stats(groups)

    def close(self):
        try:
            self._mm.close()
        finally:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


# --------------------------------------------------------------------------- #
# Ploidy helpers (caller-supplied ploidy_of factory)
# --------------------------------------------------------------------------- #
def autosomal_ploidy(sample_id, sex):
    """Every sample is diploid on an autosome/PAR."""
    return 2


def make_chr_ploidy(haploid_male=True, haploid_female=False, absent_female=False):
    """Factory for a ploidy_of returning 2 (diploid), 1 (haploid), or 0 (absent).

    * chrX-nonPAR: ``make_chr_ploidy(haploid_male=True)`` -> males haploid (1),
      females diploid (2).
    * chrY: ``make_chr_ploidy(haploid_male=True, absent_female=True)`` -> males
      haploid (1), females ABSENT (0). Females carry no Y chromosome, so a female
      samplesID row must contribute NO alleles to the chrY AN denominator;
      ``PanelIndex`` drops ploidy-0 samples from the panel entirely. (Using haploid
      females here instead would inflate chrY AN by the female count and bias every
      chrY AF low.)
    """
    def _ploidy(sample_id, sex):
        s = (sex or "").strip().lower()
        if s == "female" and absent_female:
            return 0
        if s == "male" and haploid_male:
            return 1
        if s == "female" and haploid_female:
            return 1
        return 2
    return _ploidy
