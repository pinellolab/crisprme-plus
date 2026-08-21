"""Tests for the Tier-0 registry v3 block-compression (transcode + reader).

STDLIB ONLY (unittest + tempfile). We build a small/medium v2 registry with the
existing panel-aware compile helper, transcode it to a v3 (zlib block-compressed)
registry, and prove the v3 reader is byte-for-byte semantically identical to the
v2 reader for EVERY record: keys, per-alt lookup (allele_freq via each group's
Counts), rsid, ref, plus misses. We also prove the v3 .bin shrinks, that a single
process can read a v2 and a v3 file at once, and that small block_records forces
multi-block data so cross-block lookups / alts_at straddling a block boundary are
exercised.

Reuses registry_fix_an.registries_equal as the primary equality oracle (it walks
record keys + per-group Counts + rsid + ref), and adds direct assertions on the
allele_freq path + block layout on top.
"""

import os
import tempfile
import unittest

import tier0_registry as t0
from tier0_registry import (
    RegistryReader,
    compile_registry_panel,
    transcode_registry,
    autosomal_ploidy,
    GLOBAL_GROUP_ID,
    VERSION,
    VERSION_COMPRESSED,
    CODEC_ZLIB,
    CODEC_RAW,
)
import registry_fix_an as rfa


# --------------------------------------------------------------------------- #
# Fixtures: a panel-aware registry with GLOBAL + db + db x subpop groups.
# (Same shape as test_tier0_registry.py's panel fixtures.)
# --------------------------------------------------------------------------- #
SAMPLE_META = {
    "S1": ("1000G", "EUR", "male"),
    "S2": ("1000G", "EUR", "female"),
    "S3": ("1000G", "AFR", "male"),
    "S4": ("1000G", "AFR", "female"),
    "H1": ("HGDP", "EAS", "female"),
    "H2": ("HGDP", "EAS", "male"),
}


def _small_records():
    """A handful of records incl. a multiallelic site (same pos, distinct alt)."""
    return [
        (100, "A", "G", "rs100", {"S1": "1|0", "S2": "1|1", "H1": "0|1"}),
        (100, "A", "T", "rs100b", {"S3": "1|0", "H2": "1|1"}),   # multiallelic
        (250, "C", "T", "rs250", {"S2": "0|1", "S4": "1|1"}),
        (900, "T", "G", ".", {"S1": "1|1"}),
        (900, "T", "A", "rs900a", {"H1": "1|0", "H2": "0|1"}),   # multiallelic
        (1500, "G", "C", "rs1500", {"S3": "1|0"}),
    ]


def _many_records(n=4000):
    """A few thousand records so zlib clearly wins on the record array + pool.

    Distinct rsIDs (mostly unique) + varied carriers keep it realistic while
    staying compressible (fixed-width records + repetitive small counts).
    """
    recs = []
    samples = list(SAMPLE_META.keys())
    for i in range(n):
        pos = 1000 + i  # strictly increasing -> unique (pos, alt) keys
        carriers = {}
        # a deterministic but varied carrier set per record
        for j, sid in enumerate(samples):
            if (i + j) % 3 == 0:
                carriers[sid] = "1|1" if (i + j) % 2 == 0 else "1|0"
        if not carriers:
            carriers[samples[i % len(samples)]] = "1|0"
        recs.append((pos, "A", "G", "rs%d" % (1000000 + i), carriers))
    return recs


def _compile_v2(records, tmpdir, name="reg"):
    b = os.path.join(tmpdir, name + ".bin")
    i = os.path.join(tmpdir, name + ".idx")
    compile_registry_panel(records, SAMPLE_META, None, autosomal_ploidy, b, i)
    return b, i


class _ParityMixin(object):
    """Assert a v3 reader reproduces a v2 reader exactly (keys/lookup/af/rsid/ref
    + misses). Shared by the several transcode scenarios below."""

    def assert_readers_equal(self, b2, i2, b3, i3, extra_miss=()):
        # (1) The shared oracle: record keys + per-group Counts + rsid + ref.
        ok, diffs = rfa.registries_equal(b2, i2, b3, i3)
        self.assertTrue(ok, "registries_equal diffs: %r" % (diffs[:10],))

        r2 = RegistryReader(b2, i2)
        r3 = RegistryReader(b3, i3)
        try:
            self.assertEqual(len(r2), len(r3))
            for idx in range(len(r2)):
                k2 = r2.record_key_at(idx)
                k3 = r3.record_key_at(idx)
                self.assertEqual(k2, k3, "key[%d]" % idx)
                pos, _ = k2
                # alts_at parity (order-sensitive; must match on-disk order).
                self.assertEqual(r2.alts_at(pos), r3.alts_at(pos),
                                 "alts_at(%d)" % pos)
                for alt in r2.alts_at(pos):
                    g2 = r2.lookup(pos, alt)
                    g3 = r3.lookup(pos, alt)
                    self.assertIsNotNone(g3, "v3 lookup(%d,%s) None" % (pos, alt))
                    self.assertEqual(set(g2), set(g3),
                                     "group set @(%d,%s)" % (pos, alt))
                    for gid in g2:
                        self.assertEqual(g2[gid].as_tuple(), g3[gid].as_tuple(),
                                         "counts %s @(%d,%s)" % (gid, pos, alt))
                        # exercise the allele_freq DERIVED path explicitly.
                        self.assertAlmostEqual(g2[gid].allele_freq(),
                                               g3[gid].allele_freq(),
                                               msg="af %s @(%d,%s)" % (gid, pos, alt))
                    self.assertEqual(r2.rsid(pos, alt), r3.rsid(pos, alt),
                                     "rsid @(%d,%s)" % (pos, alt))
                    self.assertEqual(r2.ref(pos, alt), r3.ref(pos, alt),
                                     "ref @(%d,%s)" % (pos, alt))
            # (2) misses: a position/alt not present must be None on BOTH readers.
            misses = [(1, "A"), (999999, "G")] + list(extra_miss)
            for (pos, alt) in misses:
                self.assertEqual(r2.lookup(pos, alt), r3.lookup(pos, alt))
                self.assertIsNone(r3.lookup(pos, alt), "expected miss (%d,%s)"
                                  % (pos, alt))
                self.assertEqual(r2.rsid(pos, alt), r3.rsid(pos, alt))
                self.assertEqual(r2.alts_at(pos), r3.alts_at(pos))
        finally:
            r2.close()
            r3.close()


class TestTranscodeParitySmall(_ParityMixin, unittest.TestCase):
    """Small panel-aware registry, transcoded with the default block size."""

    def test_v3_equals_v2_everywhere(self):
        d = tempfile.mkdtemp()
        b2, i2 = _compile_v2(_small_records(), d)
        b3 = os.path.join(d, "reg.v3.bin")
        i3 = os.path.join(d, "reg.v3.idx")
        manifest = transcode_registry(b2, i2, b3, i3)
        # v3 manifest is well-formed.
        self.assertEqual(manifest["version"], VERSION_COMPRESSED)
        self.assertEqual(manifest["format"], VERSION_COMPRESSED)
        self.assertEqual(manifest["codec"], CODEC_ZLIB)
        self.assertIn("record_blocks", manifest)
        self.assertIn("group_blocks", manifest)
        self.assertIn("pool_blocks", manifest)
        # a present-position/absent-alt miss + a present-pos miss for coverage.
        self.assert_readers_equal(b2, i2, b3, i3,
                                  extra_miss=[(100, "C"), (250, "G")])


class TestTranscodeShrinks(unittest.TestCase):
    """The v3 .bin must be SMALLER than the v2 .bin on a few-thousand-record set."""

    def test_v3_smaller_than_v2(self):
        d = tempfile.mkdtemp()
        b2, i2 = _compile_v2(_many_records(4000), d, name="big")
        b3 = os.path.join(d, "big.v3.bin")
        i3 = os.path.join(d, "big.v3.idx")
        transcode_registry(b2, i2, b3, i3)
        s2 = os.path.getsize(b2)
        s3 = os.path.getsize(b3)
        self.assertLess(s3, s2,
                        "v3 .bin (%d) not smaller than v2 (%d)" % (s3, s2))


class TestCoexistence(_ParityMixin, unittest.TestCase):
    """One process reads BOTH a v2 and a v3 registry correctly (different files)."""

    def test_v2_and_v3_in_same_process(self):
        d = tempfile.mkdtemp()
        b2, i2 = _compile_v2(_small_records(), d, name="coex")
        b3 = os.path.join(d, "coex.v3.bin")
        i3 = os.path.join(d, "coex.v3.idx")
        transcode_registry(b2, i2, b3, i3)
        r2 = RegistryReader(b2, i2)   # v2 path
        r3 = RegistryReader(b3, i3)   # v3 path
        try:
            # both open simultaneously; assert their internal framing differs...
            self.assertFalse(r2._compressed)
            self.assertTrue(r3._compressed)
            # ...yet they return identical data for a shared record.
            g2 = r2.lookup(100, "G")
            g3 = r3.lookup(100, "G")
            self.assertEqual(set(g2), set(g3))
            for gid in g2:
                self.assertEqual(g2[gid].as_tuple(), g3[gid].as_tuple())
            self.assertEqual(r2.rsid(100, "G"), r3.rsid(100, "G"))
            self.assertEqual(r2.ref(100, "T"), r3.ref(100, "T"))
        finally:
            r2.close()
            r3.close()
        # and the full parity oracle holds across both.
        self.assert_readers_equal(b2, i2, b3, i3)


class TestBlockBoundary(_ParityMixin, unittest.TestCase):
    """A small block_records forces the data to span multiple blocks, exercising
    cross-block record lookups + alts_at + group/pool reads that straddle a block
    boundary."""

    def test_multiblock_cross_boundary(self):
        d = tempfile.mkdtemp()
        records = _many_records(500)
        b2, i2 = _compile_v2(records, d, name="mb")
        b3 = os.path.join(d, "mb.v3.bin")
        i3 = os.path.join(d, "mb.v3.idx")
        # 64 records/block -> ~8 record blocks; group blob + pool chunked at
        # 64*RECORD_SIZE bytes so many records' group entries / rsIDs straddle a
        # block edge.
        manifest = transcode_registry(b2, i2, b3, i3, block_records=64)
        self.assertEqual(manifest["block_records"], 64)
        self.assertGreater(len(manifest["record_blocks"]), 1)
        self.assertGreater(len(manifest["group_blocks"]), 1)
        self.assert_readers_equal(b2, i2, b3, i3)

    def test_alts_at_straddles_block(self):
        # Put a multiallelic site whose alts land on either side of a record-block
        # boundary: with block_records=1 every record is its own block, so a
        # 2-alt site is guaranteed to straddle two record blocks.
        d = tempfile.mkdtemp()
        records = [
            (500, "A", "C", "rsC", {"S1": "1|1"}),
            (500, "A", "T", "rsT", {"S2": "0|1", "S3": "1|0"}),
            (500, "A", "G", "rsG", {"H1": "1|1"}),
        ]
        b2, i2 = _compile_v2(records, d, name="straddle")
        b3 = os.path.join(d, "straddle.v3.bin")
        i3 = os.path.join(d, "straddle.v3.idx")
        transcode_registry(b2, i2, b3, i3, block_records=1)
        r3 = RegistryReader(b3, i3)
        try:
            # alts_at must walk ACROSS the per-record blocks and return all 3.
            self.assertEqual(r3.alts_at(500), ["C", "G", "T"])
            for alt in ("C", "G", "T"):
                self.assertIsNotNone(r3.lookup(500, alt))
        finally:
            r3.close()
        self.assert_readers_equal(b2, i2, b3, i3)


class TestRawCodecFraming(_ParityMixin, unittest.TestCase):
    """codec=0 (v3 framing, uncompressed blocks) is a legal v3 file the reader
    treats via the raw-block path -- still block-indexed, no zlib."""

    def test_codec_raw_roundtrip(self):
        d = tempfile.mkdtemp()
        b2, i2 = _compile_v2(_small_records(), d, name="raw")
        b3 = os.path.join(d, "raw.v3.bin")
        i3 = os.path.join(d, "raw.v3.idx")
        manifest = transcode_registry(b2, i2, b3, i3, codec=CODEC_RAW,
                                      block_records=2)
        self.assertEqual(manifest["codec"], CODEC_RAW)
        self.assertEqual(manifest["version"], VERSION_COMPRESSED)
        r3 = RegistryReader(b3, i3)
        try:
            self.assertTrue(r3._compressed)   # v3 framing (block-indexed)
            self.assertEqual(r3._codec, CODEC_RAW)
        finally:
            r3.close()
        self.assert_readers_equal(b2, i2, b3, i3)


class TestDirectCompressedBuild(_ParityMixin, unittest.TestCase):
    """compile_registry_panel(compress=True) writes v3 directly; it must equal a
    v2 build of the same records."""

    def test_direct_build_equals_transcode(self):
        d = tempfile.mkdtemp()
        records = _small_records()
        b2, i2 = _compile_v2(records, d, name="direct2")
        b3 = os.path.join(d, "direct3.bin")
        i3 = os.path.join(d, "direct3.idx")
        compile_registry_panel(records, SAMPLE_META, None, autosomal_ploidy,
                               b3, i3, compress=True, block_records=2)
        self.assert_readers_equal(b2, i2, b3, i3)


class TestTranscodeGuards(unittest.TestCase):
    def test_rejects_already_v3(self):
        d = tempfile.mkdtemp()
        b2, i2 = _compile_v2(_small_records(), d, name="guard")
        b3 = os.path.join(d, "guard.v3.bin")
        i3 = os.path.join(d, "guard.v3.idx")
        transcode_registry(b2, i2, b3, i3)
        # transcoding a v3 source must be refused (it is not a raw v2 file).
        b3b = os.path.join(d, "guard.v3b.bin")
        i3b = os.path.join(d, "guard.v3b.idx")
        with self.assertRaises(ValueError):
            transcode_registry(b3, i3, b3b, i3b)

    def test_rejects_unknown_codec(self):
        d = tempfile.mkdtemp()
        b2, i2 = _compile_v2(_small_records(), d, name="codecguard")
        with self.assertRaises(ValueError):
            transcode_registry(b2, i2, os.path.join(d, "x.bin"),
                               os.path.join(d, "x.idx"), codec=7)


if __name__ == "__main__":
    unittest.main()
