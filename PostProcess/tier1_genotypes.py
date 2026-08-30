"""Tier-1 compact per-sample genotype store (CRISPRme+ dictless redesign, Phase 2).

The Tier-0 registry (``tier0_registry``) collapsed the 152 GB per-sample SNP dicts
into per-group AGGREGATE counts (AC/AN/carrier/hom/called). Tier-0 answers "how
common is this alt in EUR / in 1000G / globally" but DELIBERATELY forgets WHICH
individuals carry each (pos, alt): it cannot rebuild the legacy CRISPRme "Samples"
column ("HG00096:1|0,NA12878:0|1,...").

Tier-1 (THIS module) restores exactly that -- and only that -- compactly. It is a
position-indexed store of the per-record CARRIER LIST: for each (pos, specific-alt)
it records the set of alt-carrier individuals and their genotype string, so a later
wiring step can reconstruct the "sampleID:genotype" token list per alt WITHOUT the
dict. See new_simple_analysis.py 158-199 (retrieveFromDict) + 261-290: the
reconstruction target is precisely the per-alt list of "sampleID:genotype" tokens
that ``sampleSet[i]`` iterates and ``sample.split(":")`` decodes.

Design (MIRRORS tier0_registry's binary discipline):
  * A GLOBAL SAMPLE AXIS: a deterministic index<->sample_id map over the sorted,
    deduplicated union of all sample_ids, with per-index (database, subpopulation,
    sex). Persisted in the manifest. Indices fit u16 for panel scale (~4400); we
    auto-widen to u32 if the axis exceeds 65535, mirroring tier0's count_width guard.
  * A GT VOCABULARY: the distinct genotype strings ("1|0","0|1","1|1","0/1","1",...),
    each mapped to a small u8 code, stored ONCE in the manifest -- NOT per carrier.
  * A fixed-width RECORD ARRAY sorted by (pos, alt_byte), one row per (pos, alt),
    each pointing at a slice of a CARRIER BLOB. Same bisect-on-mmap lookup as tier0:
    the reader NEVER materializes a Python list of all records.
  * A CARRIER BLOB where each record's carriers are delta-varint + u8 gt_code
    encoded (ascending sample index): ~2-3 bytes/carrier vs ~12+ in the JSON dict.

NOTE (documented, NOT implemented here): small-block compression of the carrier
blob and per-DATABASE separation of the store are Phase-2 refinements. This core
step 2a is a SINGLE COMBINED, UNCOMPRESSED store over the one global sample axis.

STDLIB ONLY (struct, mmap, bisect, json, gzip) so it runs in the lightweight
unit-tests CI (no numpy / no pysam / no zstd).
"""

from __future__ import annotations

import bisect
import json
import mmap
import os
import struct

import tier0_compile as t0c

# --------------------------------------------------------------------------- #
# Binary format constants (all little-endian, mirroring tier0_registry).
# --------------------------------------------------------------------------- #
MAGIC = b"T1GT"     # Tier-1 GenoType store (distinct from tier0's b"T0RG")
VERSION = 1
# Reader accepts only VERSION; any other version is rejected with a clear message.

# Header (16 bytes, mirrors tier0's 16-byte fixed header discipline):
#   [1] MAGIC        b"T1GT"   (4 bytes)
#   [2] VERSION      u16 (=1)
#   [3] index_width  u8   (2 => u16 sample indices, 4 => u32)  -- axis-size guard
#   [4] alt_field_w  u8   (== 1 for SNPs)
#   [5] n_records    u32
#   [6] pad          (4 bytes, zero) -> 16-byte header
_HEADER_STRUCT = struct.Struct("<4sHBBIxxxx")  # magic, ver, idx_w, altw, n_rec, pad
HEADER_SIZE = _HEADER_STRUCT.size              # 16

# Record (16 bytes, fixed-width, sorted by (pos, alt_byte)):
#   I  pos            (u32; human coords < 2^32)
#   c  alt            (1 byte ascii)
#   c  ref            (1 byte ascii)
#   H  n_carriers     (u16 count of carriers for this (pos, alt))
#   I  carriers_off   (u32 absolute byte offset into the CARRIER BLOB)
#   I  carriers_len   (u32 byte length of this record's carrier slice)
# => "<IccHII" = 4+1+1+2+4+4 = 16 bytes.
_RECORD_STRUCT = struct.Struct("<IccHII")
RECORD_SIZE = _RECORD_STRUCT.size              # 16

_MAX_POS = 0xFFFFFFFF          # u32 position ceiling (all human coords fit)
_U16_MAX = 0xFFFF
_U32_MAX = 0xFFFFFFFF

# CARRIER BLOB, per record (contiguous starting at carriers_off, carriers_len bytes):
#   n_carriers                       : LEB128 unsigned varint
#   then, for each carrier, ascending sample-index order:
#       delta_sample_index           : LEB128 unsigned varint (index - prev_index;
#                                      first delta is the index itself, prev=0)
#       gt_code                      : u8 (index into the manifest gt vocab)
# Delta coding keeps the per-carrier index bytes tiny because carriers are a small,
# ascending subset of the axis. Typical cost: 1 varint byte (delta) + 1 gt byte =
# ~2 bytes/carrier; a large gap or a >65535-sample axis may push a delta to 2-3
# varint bytes -> ~3-4 bytes/carrier. Contrast the JSON dict's "HG00096:1|0,"
# (>=12 bytes/carrier).


def _write_uvarint(buf, value):
    """Append an unsigned LEB128 varint of ``value`` to bytearray ``buf``."""
    if value < 0:
        raise ValueError("uvarint cannot encode negative %d" % value)
    while True:
        byte = value & 0x7F
        value >>= 7
        if value:
            buf.append(byte | 0x80)
        else:
            buf.append(byte)
            return


def _read_uvarint(mm, off):
    """Read an unsigned LEB128 varint from ``mm`` at ``off``.

    Returns (value, next_off). Works on an mmap or a bytes/bytearray.
    """
    result = 0
    shift = 0
    while True:
        b = mm[off]
        off += 1
        result |= (b & 0x7F) << shift
        if not (b & 0x80):
            return result, off
        shift += 7


# --------------------------------------------------------------------------- #
# (1) GLOBAL SAMPLE AXIS
# --------------------------------------------------------------------------- #
def _choose_index_width(n_samples):
    """Pick the smallest safe sample-index width (2=u16, 4=u32).

    GUARD (mirrors tier0's ``_choose_count_width``): NEVER silently overflow. The
    axis holds ``n_samples`` samples, so the largest index is n_samples-1. u16 while
    that fits in 16 bits, u32 while it fits in 32 bits, else a clear error.
    """
    if n_samples <= 0:
        return 2
    max_index = n_samples - 1
    if max_index <= _U16_MAX:
        return 2
    if max_index <= _U32_MAX:
        return 4
    raise ValueError(
        "tier1_genotypes: sample axis of %d samples exceeds the u32 index ceiling "
        "(%d); implausibly large panel -- refusing to overflow"
        % (n_samples, _U32_MAX + 1)
    )


class SampleAxis(object):
    """Deterministic index<->sample_id map over the deduplicated sample union.

    ``sample_ids`` is the sorted, deduplicated list; index i <-> sample_ids[i].
    ``meta[i]`` == (database, subpopulation, sex) for sample_ids[i]. ``index_width``
    is 2 (u16) for panel scale, auto-widened to 4 (u32) past 65535 samples.
    """

    __slots__ = ("sample_ids", "meta", "_id_to_index", "index_width")

    def __init__(self, sample_ids, meta):
        self.sample_ids = list(sample_ids)
        self.meta = list(meta)
        self._id_to_index = {sid: i for i, sid in enumerate(self.sample_ids)}
        self.index_width = _choose_index_width(len(self.sample_ids))

    def __len__(self):
        return len(self.sample_ids)

    def index_of(self, sample_id):
        return self._id_to_index[sample_id]

    def get_index(self, sample_id, default=None):
        return self._id_to_index.get(sample_id, default)

    def sample_id(self, index):
        return self.sample_ids[index]

    def sample_meta(self, index):
        return self.meta[index]

    def to_manifest(self):
        """JSON-serializable axis payload for the manifest sidecar."""
        return {
            "n_samples": len(self.sample_ids),
            "index_width": self.index_width,
            # index == position in these parallel lists (the on-disk sample index).
            "sample_ids": self.sample_ids,
            # meta stored as [db, subpop, sex] per index (tuples aren't JSON).
            "sample_meta": [list(m) for m in self.meta],
        }

    @classmethod
    def from_manifest(cls, payload):
        sample_ids = list(payload["sample_ids"])
        meta = [tuple(m) for m in payload["sample_meta"]]
        axis = cls(sample_ids, meta)
        # Trust the persisted width (it was chosen at build time); it must agree
        # with what we'd recompute -- guard against a hand-edited manifest.
        persisted = int(payload["index_width"])
        if persisted != axis.index_width:
            raise ValueError(
                "tier1_genotypes: manifest index_width %d disagrees with the width "
                "implied by %d samples (%d)"
                % (persisted, len(sample_ids), axis.index_width))
        return axis


def build_sample_axis(sample_meta):
    """Build the deterministic global SampleAxis from a sample_meta mapping.

    Args:
      sample_meta: dict sample_id -> (database, subpopulation, sex), e.g. the output
        of ``tier0_compile.build_sample_meta``.

    Returns a ``SampleAxis`` over the SORTED, DEDUPLICATED union of sample_ids (a
    dict is already 1:1, so this is just a stable sort for reproducible indices).
    Per-index (database, subpopulation, sex) is carried so downstream can aggregate
    by db/subpop later (Phase-2 refinement) without re-reading samplesID files.
    """
    sample_ids = sorted(sample_meta.keys())
    meta = [tuple(sample_meta[sid]) for sid in sample_ids]
    return SampleAxis(sample_ids, meta)


# --------------------------------------------------------------------------- #
# (2) COMPACT CARRIER STORE : compiler
# --------------------------------------------------------------------------- #
def _pack_alt(alt):
    """Pack an alt string into the fixed 1-byte record slot (SNP-first)."""
    if not alt or len(alt) != 1:
        raise ValueError(
            "tier1_genotypes step 2a supports single-char alt only (per-alt SNP "
            "records); got alt=%r" % (alt,))
    return alt.encode("ascii")


def _pack_ref(ref):
    # SNP-first: a multi-char ref is an indel/MNV, out of scope here. RAISE rather
    # than silently truncate 'AT' -> 'A' (mirrors _pack_alt); the dict adapter
    # already pre-filters non-single-char refs, so this only guards direct callers.
    if ref is not None and ref != "" and len(ref) != 1:
        raise ValueError(
            "tier1_genotypes step 2a supports single-char ref only (SNP records); "
            "got ref=%r" % (ref,))
    r = (ref or ".")[:1]
    return r.encode("ascii")


def _is_snp_alt(alt):
    return bool(alt) and len(alt) == 1 and alt.upper() in ("A", "C", "G", "T")


def compile_genotypes(records, sample_axis, out_bin, out_idx):
    """Compile per-record carrier lists into the mmap-friendly binary + manifest.

    Args:
      records: iterable of (pos:int, alt:str1, carriers:dict{sample_id: gt_string}).
        Multiallelic sites => multiple records at the same pos with different alt
        (SNP-first: single-char alt). ``carriers`` lists ONLY the alt-carriers, each
        with the exact genotype string to round-trip. Need not be pre-sorted; we
        sort by (pos, alt_byte). An optional 4th element ``ref`` per tuple is honored
        for the record's stored ref base (defaults to ".").
      sample_axis: a ``SampleAxis`` covering every carrier sample_id.
      out_bin: path for the binary store.
      out_idx: path for the JSON manifest sidecar.

    Every carrier sample_id MUST be resolvable via the axis (build the axis from the
    SAME sample_meta). A carrier not on the axis raises (a silent drop would corrupt
    the reconstruction). Returns the manifest dict (also written to out_idx).
    """
    # Materialize + sort by (pos, alt_byte). ``records`` is one row per (pos, alt),
    # so this is a records-sized list at BUILD time only (the READER never does).
    prepared = []
    for rec in records:
        if len(rec) == 4:
            pos, alt, carriers, ref = rec
        else:
            pos, alt, carriers = rec
            ref = "."
        prepared.append((int(pos), alt, ref, carriers))
    prepared.sort(key=lambda r: (r[0], r[1]))
    # (pos, alt) MUST be unique: the reader bisects to exactly one row per key, so a
    # duplicate would make the second row's carriers silently unreachable. The
    # bcftools norm -m- biallelic dict never repeats a (pos, alt); this guards a
    # future non-dict caller from a silent carrier drop.
    for _i in range(1, len(prepared)):
        if (prepared[_i][0], prepared[_i][1]) == (prepared[_i - 1][0], prepared[_i - 1][1]):
            raise ValueError(
                "tier1_genotypes: duplicate (pos, alt) record (%s, %r) -- each "
                "(pos, alt) must be unique" % (prepared[_i][0], prepared[_i][1]))

    # GT vocabulary: intern distinct genotype strings -> small u8 codes, in
    # first-seen order (deterministic given the sorted record order).
    gt_code = {}       # gt_string -> code
    gt_vocab = []      # code -> gt_string (the manifest table)

    def _intern_gt(gt):
        code = gt_code.get(gt)
        if code is None:
            code = len(gt_vocab)
            if code > 0xFF:
                raise ValueError(
                    "tier1_genotypes: GT vocabulary exceeds 256 distinct genotype "
                    "strings (u8 gt_code ceiling); saw %r" % (gt,))
            gt_code[gt] = code
            gt_vocab.append(gt)
        return code

    idx_width = sample_axis.index_width

    carrier_blob = bytearray()
    record_rows = []   # (pos, altb, refb, n_carriers, carriers_off, carriers_len)

    for (pos, alt, ref, carriers) in prepared:
        if pos < 0 or pos > _MAX_POS:
            raise ValueError(
                "tier1_genotypes: position %d does not fit in u32 (0..%d)"
                % (pos, _MAX_POS))
        altb = _pack_alt(alt)
        refb = _pack_ref(ref)

        # Resolve carriers to (index, gt_code), then sort ASCENDING by index so the
        # delta coding + reconstruction order are deterministic.
        resolved = []
        for sid, gt in carriers.items():
            index = sample_axis.get_index(sid)
            if index is None:
                raise ValueError(
                    "tier1_genotypes: carrier %r at pos %d alt %s is not on the "
                    "sample axis (build the axis from the same sample_meta)"
                    % (sid, pos, alt))
            resolved.append((index, _intern_gt(gt)))
        resolved.sort(key=lambda t: t[0])

        n_carriers = len(resolved)
        if n_carriers > _U16_MAX:
            raise ValueError(
                "tier1_genotypes: record at pos %d alt %s has %d carriers "
                "(> u16 max %d)" % (pos, alt, n_carriers, _U16_MAX))

        carriers_off = len(carrier_blob)
        _write_uvarint(carrier_blob, n_carriers)
        prev = 0
        for index, code in resolved:
            delta = index - prev  # ascending -> non-negative
            _write_uvarint(carrier_blob, delta)
            carrier_blob.append(code & 0xFF)
            prev = index
        carriers_len = len(carrier_blob) - carriers_off

        record_rows.append((pos, altb, refb, n_carriers, carriers_off, carriers_len))

    n_records = len(record_rows)

    # Section offsets. carriers_off in each row is relative to the blob start; we
    # store ABSOLUTE file offsets so the reader seeks with one addend.
    record_array_off = HEADER_SIZE
    carrier_blob_off = record_array_off + n_records * RECORD_SIZE

    with open(out_bin, "wb") as fh:
        fh.write(_HEADER_STRUCT.pack(MAGIC, VERSION, idx_width, 1, n_records))
        for (pos, altb, refb, n_carriers, coff, clen) in record_rows:
            fh.write(_RECORD_STRUCT.pack(
                pos, altb, refb, n_carriers, carrier_blob_off + coff, clen))
        fh.write(bytes(carrier_blob))

    manifest = {
        "magic": MAGIC.decode("ascii"),
        "version": VERSION,
        "n_records": n_records,
        "record_array_off": record_array_off,
        "record_size": RECORD_SIZE,
        "carrier_blob_off": carrier_blob_off,
        "carrier_blob_len": len(carrier_blob),
        "index_width": idx_width,
        "alt_field_width": 1,
        # GT vocabulary: index == on-disk gt_code -> genotype string. Stored ONCE.
        "gt_vocab": gt_vocab,
        # Global sample axis (index<->sample_id + per-index db/subpop/sex).
        "sample_axis": sample_axis.to_manifest(),
        # Documented Phase-2 refinements NOT applied to this core store.
        "compression": "none",
        "store_layout": "combined-global-axis",
    }
    # .idx written LAST and ATOMICALLY: a present .idx marks this chromosome's store
    # complete (the .bin is closed above), so a resume keys on .idx and never skips a
    # kill-truncated .bin. tmp+replace keeps a partial .idx off the final path.
    _idx_tmp = out_idx + ".tmp"
    with open(_idx_tmp, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    os.replace(_idx_tmp, out_idx)

    return manifest


def compile_genotypes_from_dict(dict_path, db_to_samplesid, chrom, out_bin, out_idx,
                                *, subpop_field="superpopulation"):
    """Compile a legacy SNP dict into a Tier-1 genotype store.

    Builds the global sample axis via ``tier0_compile.build_sample_meta`` (reading
    the per-database samplesID files), streams the dict with
    ``tier0_compile.iter_dict_records``, and compiles each (pos, alt) carrier list.

    Non-SNP alts (multi-char alt, multi-char ref, or a non-ACGT placeholder such as
    the '*' spanning-deletion marker) are SKIPPED with a counter (never crash),
    mirroring ``tier0_compile.compile_from_dict``. Records with no listed carriers
    are also skipped (nothing to reconstruct).

    Returns a dict of build stats:
      {"manifest", "n_written", "n_skipped_indel", "n_skipped_empty",
       "n_positions", "overlaps", "n_samples"}.
    """
    sample_meta, overlaps = t0c.build_sample_meta(
        db_to_samplesid, subpop_field=subpop_field)
    sample_axis = build_sample_axis(sample_meta)

    stats = {
        "n_written": 0,
        "n_skipped_indel": 0,
        "n_skipped_empty": 0,
        "n_positions": 0,
    }
    seen_positions = set()

    def record_stream():
        for (pos, ref, alt, rsid, carrier_gts) in t0c.iter_dict_records(dict_path, chrom):
            seen_positions.add(pos)
            # SNP-first: single-base substitution only. Skip indels / non-ACGT.
            if not _is_snp_alt(alt) or (ref and len(ref) != 1):
                stats["n_skipped_indel"] += 1
                continue
            # Only keep carriers actually on the axis (all should be, but a dict
            # sample absent from every samplesID cannot be placed and would corrupt
            # reconstruction -> drop it defensively rather than raise here).
            filtered = {sid: gt for sid, gt in carrier_gts.items()
                        if sample_axis.get_index(sid) is not None}
            if not filtered:
                stats["n_skipped_empty"] += 1
                continue
            stats["n_written"] += 1
            # (pos, alt, carriers, ref) 4-tuple form so the ref base is preserved.
            yield (pos, alt, filtered, ref)

    manifest = compile_genotypes(record_stream(), sample_axis, out_bin, out_idx)

    stats["n_positions"] = len(seen_positions)
    stats["manifest"] = manifest
    stats["overlaps"] = overlaps
    stats["n_samples"] = len(sample_axis)
    return stats


# --------------------------------------------------------------------------- #
# (3) mmap + bisect READER
# --------------------------------------------------------------------------- #
class GenotypeReader(object):
    """mmap-backed reader. Binary-searches a sorted fixed-width record array.

    Does NOT load the record array into a Python list/dict; ``carriers`` uses
    ``bisect`` over an on-mmap key accessor (mirrors tier0_registry.RegistryReader).
    """

    def __init__(self, bin_path, idx_path):
        with open(idx_path, "r") as fh:
            self.manifest = json.load(fh)
        self._n = self.manifest["n_records"]
        self._rec_off = self.manifest["record_array_off"]
        self._rec_size = self.manifest["record_size"]
        self._carrier_blob_off = self.manifest["carrier_blob_off"]
        self._index_width = self.manifest["index_width"]
        # GT vocab (code -> gt_string) and the sample axis (index -> id/meta) are
        # small (vocab <= 256; axis ~ n_samples, NOT n_records) so caching them does
        # not violate the "no per-record materialization" contract. We rebuild the
        # axis object once for id/meta accessors.
        self._gt_vocab = self.manifest["gt_vocab"]
        self._axis = SampleAxis.from_manifest(self.manifest["sample_axis"])

        self._fh = open(bin_path, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        # Header sanity-check (also rejects unknown / wrong-version files clearly).
        # On failure, close the just-opened mmap + handle before re-raising so a
        # rejected file does not leak a descriptor.
        try:
            magic, version, idx_w, altw, n_records = \
                _HEADER_STRUCT.unpack_from(self._mm, 0)
            if magic != MAGIC:
                raise ValueError("bad magic in %s: %r" % (bin_path, magic))
            if version != VERSION:
                raise ValueError(
                    "tier1_genotypes: unsupported format version %d in %s (this "
                    "reader supports v%d only; rebuild the store)"
                    % (version, bin_path, VERSION))
            if n_records != self._n:
                raise ValueError("n_records mismatch between bin and idx")
            if idx_w != self._index_width:
                raise ValueError("header/manifest index_width mismatch in %s"
                                 % bin_path)
        except Exception:
            try:
                self._mm.close()
            finally:
                self._fh.close()
            raise

        # bisect helper: a virtual sorted sequence of record (pos, alt_byte) keys,
        # computed on demand from the mmap -- no Python list of all records.
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
        pos = struct.unpack_from("<I", self._mm, base)[0]  # u32 pos
        altb = self._mm[base + 4:base + 5]                 # single byte after pos
        return (pos, altb)

    def _record_at(self, i):
        base = self._record_base(i)
        (pos, altb, refb, n_carriers, carriers_off, carriers_len) = \
            _RECORD_STRUCT.unpack_from(self._mm, base)
        return pos, altb, refb, n_carriers, carriers_off, carriers_len

    def _find_index(self, pos, alt):
        key = (int(pos), alt.encode("ascii"))
        i = bisect.bisect_left(self._keys, key)
        if i < self._n and self._key_at(i) == key:
            return i
        return None

    def _decode_carriers(self, n_carriers, carriers_off):
        """Decode a record's carrier slice into [(sample_index, gt_code), ...] in
        ascending sample-index order, directly off the mmap (no whole-blob parse)."""
        off = carriers_off
        stored_n, off = _read_uvarint(self._mm, off)
        # stored_n (in the varint) and n_carriers (in the record) agree by
        # construction; trust the record's u16 for the loop count and cross-check.
        if stored_n != n_carriers:
            raise ValueError(
                "tier1_genotypes: carrier-count mismatch (record=%d blob=%d) -- "
                "corrupt store" % (n_carriers, stored_n))
        out = []
        prev = 0
        for _ in range(n_carriers):
            delta, off = _read_uvarint(self._mm, off)
            index = prev + delta
            code = self._mm[off]
            off += 1
            out.append((index, code))
            prev = index
        return out

    # ---- public API ---- #
    def __len__(self):
        return self._n

    def carriers(self, pos, alt):
        """Return the reconstructed carrier list for (pos, alt), or None if absent.

        Result: list of (sample_id:str, gt_string:str) in ASCENDING sample-index
        order (deterministic). sample_id comes from the axis; gt_string from the GT
        vocab -- an EXACT round-trip of what was compiled.
        """
        i = self._find_index(pos, alt)
        if i is None:
            return None
        _, _, _, n_carriers, carriers_off, _ = self._record_at(i)
        decoded = self._decode_carriers(n_carriers, carriers_off)
        return [(self._axis.sample_id(index), self._gt_vocab[code])
                for (index, code) in decoded]

    def carrier_tokens(self, pos, alt):
        """["sample_id:gt", ...] for (pos, alt) (legacy per-alt token list), or None.

        This is the exact form ``sampleSet[i]`` iterates in new_simple_analysis.py
        (each token later split on ':'), so a wiring step drops it straight into
        ``sample_list``. Ascending sample-index order (deterministic)."""
        recon = self.carriers(pos, alt)
        if recon is None:
            return None
        return ["%s:%s" % (sid, gt) for (sid, gt) in recon]

    def alts_at(self, pos):
        """Alt alleles present at ``pos`` in on-disk (sorted) order (bisect + walk).

        Multiallelic sites are consecutive fixed-width records with the same pos and
        distinct alt. O(log n + #alts), no full parse."""
        pos = int(pos)
        lo = bisect.bisect_left(self._keys, (pos, b""))
        alts = []
        prev = None
        i = lo
        while i < self._n:
            p, altb = self._key_at(i)
            if p != pos:
                break
            a = altb.decode("ascii")
            if a != prev:
                alts.append(a)
                prev = a
            i += 1
        return alts

    def ref(self, pos, alt):
        """Reference base for (pos, alt) as a str, or None if the record is absent."""
        i = self._find_index(pos, alt)
        if i is None:
            return None
        _, _, refb, _, _, _ = self._record_at(i)
        return refb.decode("ascii")

    def record_key_at(self, i):
        """(pos, alt_str) at index i -- for iteration / first/last tests."""
        pos, altb = self._key_at(i)
        return pos, altb.decode("ascii")

    # ---- sample-axis accessors ---- #
    def n_samples(self):
        return len(self._axis)

    def sample_id(self, index):
        return self._axis.sample_id(index)

    def sample_meta(self, index):
        """(database, subpopulation, sex) for a sample index."""
        return self._axis.sample_meta(index)

    def gt_vocab(self):
        """The distinct genotype strings; list index == on-disk gt_code."""
        return list(self._gt_vocab)

    def index_width_used(self):
        """On-disk sample-index width in bytes (2 => u16, 4 => u32)."""
        return self._index_width

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
