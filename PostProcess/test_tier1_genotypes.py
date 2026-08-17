"""Unit tests for the Tier-1 compact per-sample genotype store (Phase 2 step 2a).

STDLIB ONLY (unittest + tempfile + json). Covers:
  (a) synthetic carrier round-trip incl. a multiallelic '$'-style pos (two alts with
      different carrier sets), a haploid-style "1" gt, phased "1|0"/"0|1", and hom
      "1|1" -- assert carrier_tokens EXACTLY equal the input "sample_id:gt" lists in
      ascending sample-index order;
  (b) absent (pos, alt) -> None; first / last record keys;
  (c) a dict fixture: write a tiny my_dict_<chrom>.json + two per-db samplesID files,
      compile via compile_genotypes_from_dict, assert carrier_tokens match the dict's
      per-alt carrier lists; a non-SNP alt is skipped;
  (d) mmap/bisect structural check (reader holds no per-record Python sequence);
  (e) sample-index auto-widen past u16 (a >65535-sample synthetic axis with a small
      record set) plus an over-u32 guard.
"""

from __future__ import annotations

import json
import mmap as _mmap
import os
import tempfile
import unittest

import tier1_genotypes as t1


# --------------------------------------------------------------------------- #
# Fixture writers (match tier0_compile's samplesID + dict formats exactly).
# --------------------------------------------------------------------------- #
def write_samplesid(path, rows):
    """rows: list of (sample_id, population, superpopulation, sex)."""
    with open(path, "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for sid, pop, spop, sex in rows:
            fh.write("\t".join([sid, pop, spop, sex]) + "\n")


def write_dict(path, mapping):
    """mapping: {"<chrom>,<pos>": "<entry>$<entry>..."} -> json file."""
    with open(path, "w") as fh:
        json.dump(mapping, fh)


def _tokens(carrier_dict, axis):
    """The EXPECTED carrier_tokens for a carrier dict: "sid:gt" in ascending
    sample-index order (the store's deterministic order)."""
    ordered = sorted(carrier_dict.items(), key=lambda kv: axis.index_of(kv[0]))
    return ["%s:%s" % (sid, gt) for sid, gt in ordered]


# --------------------------------------------------------------------------- #
# (a)(b)(d) Synthetic round-trip via compile_genotypes.
# --------------------------------------------------------------------------- #
class TestSyntheticRoundTrip(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.binp = os.path.join(self.d, "gt.bin")
        self.idxp = os.path.join(self.d, "gt.idx.json")
        # 6-sample panel across two dbs / three subpops / both sexes.
        self.sample_meta = {
            "S_EUR_M1": ("1000G", "EUR", "male"),
            "S_EUR_F1": ("1000G", "EUR", "female"),
            "S_AFR_M1": ("1000G", "AFR", "male"),
            "S_AFR_F1": ("1000G", "AFR", "female"),
            "H_EAS_M1": ("HGDP", "EAS", "male"),
            "H_EAS_F1": ("HGDP", "EAS", "female"),
        }
        self.axis = t1.build_sample_axis(self.sample_meta)

        # Records. Includes: phased "1|0"/"0|1", hom "1|1", haploid "1"; a
        # multiallelic pos 300 with two alts (G and T) that have DIFFERENT carriers.
        self.recs_by_key = {
            # pos 100, alt A: phased + hom mix
            (100, "A"): {"S_EUR_M1": "1|0", "S_EUR_F1": "0|1", "S_AFR_F1": "1|1"},
            # pos 100, alt C: a second alt at the SAME position (multiallelic)
            (100, "C"): {"H_EAS_F1": "0|1"},
            # pos 300, alt G: one carrier set
            (300, "G"): {"S_EUR_M1": "1|0", "H_EAS_M1": "1"},   # haploid "1"
            # pos 300, alt T: a DIFFERENT carrier set at the same pos
            (300, "T"): {"S_AFR_M1": "1", "S_EUR_F1": "1|1", "H_EAS_F1": "0/1"},
            # pos 500, alt A: unphased "0/1"
            (500, "A"): {"S_AFR_F1": "0/1"},
        }
        # ref bases per record (round-tripped too).
        refs = {(100, "A"): "G", (100, "C"): "G", (300, "G"): "C",
                (300, "T"): "C", (500, "A"): "T"}
        records = [(pos, alt, carriers, refs[(pos, alt)])
                   for (pos, alt), carriers in self.recs_by_key.items()]
        # Deliberately pass UNSORTED input; the compiler must sort by (pos, alt).
        records = list(reversed(records))
        t1.compile_genotypes(records, self.axis, self.binp, self.idxp)

    def test_multichar_ref_raises(self):
        # SNP-first: a multi-char ref (indel) must be REJECTED, not silently
        # truncated 'AT' -> 'A'. (The dict adapter pre-filters these; guard direct callers.)
        with self.assertRaises(ValueError):
            t1.compile_genotypes([(700, "A", {"S_EUR_M1": "1|0"}, "AT")],
                                 self.axis, self.binp + ".x", self.idxp + ".x")

    def test_duplicate_pos_alt_raises(self):
        # Two records with the SAME (pos, alt) would make the second unreachable
        # via bisect -> reject at compile so carriers are never silently dropped.
        with self.assertRaises(ValueError):
            t1.compile_genotypes(
                [(800, "A", {"S_EUR_M1": "1|0"}), (800, "A", {"S_AFR_F1": "0|1"})],
                self.axis, self.binp + ".y", self.idxp + ".y")

    def test_exact_carrier_token_roundtrip(self):
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            for (pos, alt), carriers in self.recs_by_key.items():
                expected = _tokens(carriers, self.axis)
                got = r.carrier_tokens(pos, alt)
                self.assertEqual(got, expected,
                                 "carrier_tokens mismatch at %d %s" % (pos, alt))
                # carriers() tuples must match too (id, gt) in the same order.
                recon = r.carriers(pos, alt)
                self.assertEqual(
                    recon,
                    [(sid, carriers[sid]) for sid in
                     sorted(carriers, key=self.axis.index_of)])
        finally:
            r.close()

    def test_multiallelic_alts_are_independent(self):
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            # Same pos 300, two alts with different carriers.
            self.assertEqual(r.alts_at(300), ["G", "T"])
            self.assertEqual(r.alts_at(100), ["A", "C"])
            self.assertNotEqual(r.carrier_tokens(300, "G"),
                                r.carrier_tokens(300, "T"))
            # A haploid "1" gt round-trips verbatim (no forced diploidization).
            self.assertIn("H_EAS_M1:1", r.carrier_tokens(300, "G"))
            self.assertIn("S_AFR_M1:1", r.carrier_tokens(300, "T"))
        finally:
            r.close()

    def test_ascending_index_order(self):
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            recon = r.carriers(300, "T")
            indices = [self.axis.index_of(sid) for sid, _ in recon]
            self.assertEqual(indices, sorted(indices))
        finally:
            r.close()

    def test_absent_lookups_return_none(self):
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            self.assertIsNone(r.carrier_tokens(999, "A"))   # absent pos
            self.assertIsNone(r.carriers(999, "A"))
            self.assertIsNone(r.carrier_tokens(100, "T"))   # present pos, absent alt
            self.assertIsNone(r.carrier_tokens(300, "A"))   # present pos, absent alt
            self.assertIsNone(r.ref(999, "A"))
        finally:
            r.close()

    def test_first_and_last_record(self):
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            n = len(r)
            self.assertEqual(n, 5)
            # Sorted by (pos, alt_byte): first == (100, "A"), last == (500, "A").
            self.assertEqual(r.record_key_at(0), (100, "A"))
            self.assertEqual(r.record_key_at(n - 1), (500, "A"))
            self.assertIsNotNone(r.carrier_tokens(100, "A"))
            self.assertIsNotNone(r.carrier_tokens(500, "A"))
        finally:
            r.close()

    def test_ref_and_meta_accessors(self):
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            self.assertEqual(r.ref(100, "A"), "G")
            self.assertEqual(r.ref(300, "T"), "C")
            self.assertEqual(r.n_samples(), 6)
            # index<->id<->meta round-trip through the reader's axis.
            idx = self.axis.index_of("H_EAS_F1")
            self.assertEqual(r.sample_id(idx), "H_EAS_F1")
            self.assertEqual(r.sample_meta(idx), ("HGDP", "EAS", "female"))
        finally:
            r.close()

    def test_gt_vocab_stored_once(self):
        # The distinct gt strings are interned into a small vocab (not per carrier).
        with open(self.idxp) as fh:
            manifest = json.load(fh)
        vocab = manifest["gt_vocab"]
        self.assertEqual(sorted(vocab),
                         sorted({"1|0", "0|1", "1|1", "1", "0/1"}))
        # Vocab is far smaller than the total carrier count.
        total_carriers = sum(len(c) for c in self.recs_by_key.values())
        self.assertLess(len(vocab), total_carriers)

    def test_context_manager(self):
        with t1.GenotypeReader(self.binp, self.idxp) as r:
            self.assertEqual(len(r), 5)

    def test_reader_does_not_materialize_records(self):
        # Structural guarantee (mirrors tier0): the reader keeps no Python
        # list/dict of all records; only the mmap + manifest + small axis/vocab.
        r = t1.GenotypeReader(self.binp, self.idxp)
        try:
            n = len(r)
            self.assertTrue(
                any(isinstance(v, _mmap.mmap) for v in vars(r).values()),
                "reader must hold an mmap")
            for name, val in vars(r).items():
                if name in ("manifest", "_gt_vocab", "_axis"):
                    continue  # small n_samples/vocab-sized state, NOT per-record
                if isinstance(val, (list, tuple)):
                    self.fail("reader holds a Python sequence attr %r "
                              "(records must stay on the mmap)" % name)
                if isinstance(val, dict) and n > 1 and len(val) >= n:
                    self.fail("reader materialized a record-sized dict: %s" % name)
            # The key view is a lazy accessor, not a real list.
            self.assertFalse(isinstance(r._keys, (list, tuple)))
            self.assertEqual(len(r._keys), n)
            self.assertIsNotNone(r.carrier_tokens(100, "A"))  # bisect still works
        finally:
            r.close()

    def test_version_reject(self):
        # Corrupt the on-disk VERSION and confirm the reader rejects it clearly.
        bad = os.path.join(self.d, "bad.bin")
        with open(self.binp, "rb") as fh:
            data = bytearray(fh.read())
        # VERSION is a u16 right after the 4-byte magic.
        data[4:6] = (t1.VERSION + 7).to_bytes(2, "little")
        with open(bad, "wb") as fh:
            fh.write(bytes(data))
        with self.assertRaises(ValueError) as ctx:
            t1.GenotypeReader(bad, self.idxp)
        self.assertIn("version", str(ctx.exception).lower())

    def test_bad_magic_reject(self):
        bad = os.path.join(self.d, "badmagic.bin")
        with open(self.binp, "rb") as fh:
            data = bytearray(fh.read())
        data[0:4] = b"XXXX"
        with open(bad, "wb") as fh:
            fh.write(bytes(data))
        with self.assertRaises(ValueError) as ctx:
            t1.GenotypeReader(bad, self.idxp)
        self.assertIn("magic", str(ctx.exception).lower())


# --------------------------------------------------------------------------- #
# (c) Dict-fixture round-trip via compile_genotypes_from_dict.
# --------------------------------------------------------------------------- #
G1000_ROWS = [
    ("S_EUR_M1", "GBR", "EUR", "male"),
    ("S_EUR_F1", "GBR", "EUR", "female"),
    ("S_AFR_M1", "YRI", "AFR", "male"),
    ("S_AFR_F1", "YRI", "AFR", "female"),
]
HGDP_ROWS = [
    ("H_EAS_M1", "Han", "EAS", "male"),
    ("H_EAS_F1", "Han", "EAS", "female"),
]


class TestCompileFromDict(unittest.TestCase):
    CHROM = "chr1"

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.g1000_sid = os.path.join(self.d, "hg38_1000G.samplesID.txt")
        self.hgdp_sid = os.path.join(self.d, "hg38_HGDP.samplesID.txt")
        write_samplesid(self.g1000_sid, G1000_ROWS)
        write_samplesid(self.hgdp_sid, HGDP_ROWS)
        self.db_to_sid = {"1000G": self.g1000_sid, "HGDP": self.hgdp_sid}
        # Hand-built axis for computing expected ordered tokens (independent of
        # the module's build, though it must match).
        meta = {}
        for sid, pop, spop, sex in G1000_ROWS:
            meta[sid] = ("1000G", spop, sex)
        for sid, pop, spop, sex in HGDP_ROWS:
            meta[sid] = ("HGDP", spop, sex)
        self.axis = t1.build_sample_axis(meta)

        # Dict entry format:  "<samples>;<ref,alt>;<rsID>;<AF>"
        #   <samples> = comma-joined "sampleID:genotype" (alt-carriers only).
        # pos 100: single SNP (A>G).
        e100 = "S_EUR_M1:1|0,S_AFR_F1:1|1;A,G;rs100;0.5"
        # pos 300: MULTIALLELIC '$'-joined -> alt C and alt T, DIFFERENT carriers.
        e300_C = "S_EUR_F1:0|1,H_EAS_M1:1;A,C;rs300a;0.3"
        e300_T = "S_AFR_M1:1,H_EAS_F1:0/1;A,T;rs300b;0.2"
        # pos 500: a NON-SNP alt (indel) -> must be SKIPPED.
        e500_indel = "S_EUR_M1:1|0;AT,A;rs500;0.1"
        # pos 700: an entry whose carriers are ALL off-axis -> skipped as empty.
        e700_offaxis = "GHOST1:1|0,GHOST2:0|1;C,G;rs700;0.01"

        self.mapping = {
            "%s,100" % self.CHROM: e100,
            "%s,300" % self.CHROM: e300_C + "$" + e300_T,
            "%s,500" % self.CHROM: e500_indel,
            "%s,700" % self.CHROM: e700_offaxis,
        }
        # Expected carrier dicts per (pos, alt) for the SNP records.
        self.expected = {
            (100, "G"): {"S_EUR_M1": "1|0", "S_AFR_F1": "1|1"},
            (300, "C"): {"S_EUR_F1": "0|1", "H_EAS_M1": "1"},
            (300, "T"): {"S_AFR_M1": "1", "H_EAS_F1": "0/1"},
        }

    def _compile(self):
        dpath = os.path.join(self.d, "my_dict_%s.json" % self.CHROM)
        write_dict(dpath, self.mapping)
        binp = os.path.join(self.d, "gt.bin")
        idxp = os.path.join(self.d, "gt.idx.json")
        stats = t1.compile_genotypes_from_dict(
            dpath, self.db_to_sid, self.CHROM, binp, idxp)
        return binp, idxp, stats

    def test_dict_carrier_tokens_match(self):
        binp, idxp, stats = self._compile()
        r = t1.GenotypeReader(binp, idxp)
        try:
            for (pos, alt), carriers in self.expected.items():
                expected = _tokens(carriers, self.axis)
                self.assertEqual(r.carrier_tokens(pos, alt), expected,
                                 "dict token mismatch at %d %s" % (pos, alt))
            # The two 300 alts stay independent.
            self.assertEqual(r.alts_at(300), ["C", "T"])
            # pos 500 (indel) produced NO record.
            self.assertEqual(r.alts_at(500), [])
            self.assertIsNone(r.carrier_tokens(500, "A"))
            # pos 700 (all off-axis) produced NO record.
            self.assertEqual(r.alts_at(700), [])
        finally:
            r.close()

    def test_dict_stats(self):
        binp, idxp, stats = self._compile()
        # 3 SNP records written (100/G, 300/C, 300/T).
        self.assertEqual(stats["n_written"], 3)
        # 1 indel skipped (500) + 1 all-off-axis skipped (700).
        self.assertEqual(stats["n_skipped_indel"], 1)
        self.assertEqual(stats["n_skipped_empty"], 1)
        self.assertEqual(stats["n_samples"], 6)
        self.assertEqual(stats["manifest"]["n_records"], 3)
        # 4 dict positions seen.
        self.assertEqual(stats["n_positions"], 4)

    def test_dict_gt_vocab(self):
        binp, idxp, stats = self._compile()
        r = t1.GenotypeReader(binp, idxp)
        try:
            self.assertEqual(sorted(r.gt_vocab()),
                             sorted({"1|0", "1|1", "0|1", "1", "0/1"}))
        finally:
            r.close()


# --------------------------------------------------------------------------- #
# (e) Sample-index width guard / auto-widen.
# --------------------------------------------------------------------------- #
class TestIndexWidth(unittest.TestCase):
    def test_u16_default(self):
        meta = {"S%05d" % i: ("db", "sp", "male") for i in range(100)}
        axis = t1.build_sample_axis(meta)
        self.assertEqual(axis.index_width, 2)

    def test_auto_widen_past_u16(self):
        # A >65535-sample axis forces u32 indices, with a SMALL record set so the
        # test stays fast. We reach past the u16 ceiling and store a record whose
        # carrier is a high index (> 65535) to exercise the wide delta path.
        self.d = tempfile.mkdtemp()
        n = 70000  # > 65535 => index_width must widen to 4
        meta = {"S%06d" % i: ("db", "sp", "male") for i in range(n)}
        axis = t1.build_sample_axis(meta)
        self.assertEqual(axis.index_width, 4)
        self.assertEqual(len(axis), n)

        high_id = "S%06d" % (n - 1)     # index 69999 (> 65535)
        low_id = "S%06d" % 3
        carriers = {high_id: "1|1", low_id: "0|1"}
        binp = os.path.join(self.d, "wide.bin")
        idxp = os.path.join(self.d, "wide.idx.json")
        manifest = t1.compile_genotypes(
            [(42, "A", carriers)], axis, binp, idxp)
        self.assertEqual(manifest["index_width"], 4)

        r = t1.GenotypeReader(binp, idxp)
        try:
            self.assertEqual(r.index_width_used(), 4)
            expected = _tokens(carriers, axis)
            self.assertEqual(r.carrier_tokens(42, "A"), expected)
            # The high-index carrier round-trips exactly.
            self.assertIn("%s:1|1" % high_id, r.carrier_tokens(42, "A"))
        finally:
            r.close()

    def test_over_u32_guard(self):
        # A guard test: an axis past the u32 index ceiling must RAISE, not overflow.
        # We do NOT allocate 4.3e9 samples; we call the width chooser directly.
        with self.assertRaises(ValueError) as ctx:
            t1._choose_index_width(0xFFFFFFFF + 2)
        self.assertIn("u32", str(ctx.exception).lower())
        # And just under the ceiling is fine (returns 4).
        self.assertEqual(t1._choose_index_width(0xFFFFFFFF + 1), 4)

    def test_gt_vocab_over_u8_guard(self):
        # More than 256 distinct gt strings would overflow the u8 gt_code; guard it.
        meta = {"S%05d" % i: ("db", "sp", "male") for i in range(300)}
        axis = t1.build_sample_axis(meta)
        # One record whose 300 carriers each have a UNIQUE gt string.
        carriers = {"S%05d" % i: "g%d|0" % i for i in range(300)}
        d = tempfile.mkdtemp()
        binp = os.path.join(d, "v.bin")
        idxp = os.path.join(d, "v.idx.json")
        with self.assertRaises(ValueError) as ctx:
            t1.compile_genotypes([(1, "A", carriers)], axis, binp, idxp)
        self.assertIn("256", str(ctx.exception))


if __name__ == "__main__":
    unittest.main(verbosity=2)
