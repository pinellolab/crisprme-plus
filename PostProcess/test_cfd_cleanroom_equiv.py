"""Bit-identity gate for the CFD clean-room (2.8.0).

The live CFD path used to inject a per-module ``calc_cfd`` (an inherited transliteration)
into ``bulge_cfd_penalty.cfd_score_bulge_aware`` as the mismatch-only base scorer. 2.8.0
replaces that base with the own-implementation ``cfd_score.cfd`` (Doench et al.,
Nat. Biotechnol. 2016), scored from the committed JSON tables, and the inherited
``calc_cfd`` / ``revcom`` / ``get_mm_pam_scores`` / ``.pkl`` are removed from the tree.

This test pins the clean-room path to the historical behaviour via a GOLDEN DATA fixture
(``cfd_golden.tsv.gz``): ``(aligned_crRNA, aligned_protospacer) -> expected CFD`` rows,
each expectation the raw IEEE-754 bit pattern produced by the former
``cfd_score_bulge_aware`` (legacy base) before the swap. The fixture is pure numbers --
no inherited code -- so it keeps proving equivalence without re-introducing the lineage.
It spans mismatch-only, DNA/RNA/terminal/PAM-proximal bulges, every gap key, IUPAC/N and
non-canonical PAM. The suite also checks clean-room internal consistency (``cfd_score.cfd``
vs ``twopass_cfd_exact.cfd_concrete``) and the published Doench anchor values.

Run: ``<env>/bin/python3 -m unittest test_cfd_cleanroom_equiv -v``
"""

import gzip
import os
import random
import struct
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bulge_cfd_penalty as bp  # noqa: E402
import twopass_cfd_exact as tp  # noqa: E402
from cfd_score import cfd as cleanroom_cfd, cfd_tables  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_GOLDEN = os.path.join(_HERE, "cfd_golden.tsv.gz")
_BASES = "ACGT"


def _bits(x):
    return struct.pack("<d", float(x))


class CfdGoldenBitIdentity(unittest.TestCase):
    """The production clean-room path reproduces the historical CFD bit-for-bit."""

    def test_golden_fixture_reproduced_exactly(self):
        self.assertTrue(os.path.exists(_GOLDEN), "missing golden fixture %s" % _GOLDEN)
        n = 0
        with gzip.open(_GOLDEN, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cr, pr, want_hex = line.rstrip("\n").split("\t")
                got = bp.cfd_score_bulge_aware(cr, pr, pam_len=3)  # clean-room default
                self.assertEqual(_bits(got).hex(), want_hex,
                                 "golden mismatch cr=%r pr=%r got=%r" % (cr, pr, got))
                n += 1
        self.assertGreater(n, 6000, "golden fixture unexpectedly small (%d rows)" % n)


class CfdCleanroomInternalConsistency(unittest.TestCase):
    """The independent clean-room scorers agree bit-for-bit, so none silently drifts."""

    def test_cfd_score_equals_cfd_concrete_mismatch_only(self):
        tables = cfd_tables()
        mm, pam = tables
        rng = random.Random(99)
        for _ in range(30000):
            guide = "".join(rng.choice(_BASES) for _ in range(20))
            proto = list(guide)
            for p in rng.sample(range(20), rng.randint(0, 6)):
                proto[p] = rng.choice(_BASES)
            proto = "".join(proto)
            pamk = rng.choice(_BASES) + rng.choice(_BASES)
            a = cleanroom_cfd(guide, proto, pamk, tables)
            b = tp.cfd_concrete(guide, proto, pamk, mm, pam)
            self.assertEqual(_bits(a), _bits(b),
                             "cfd_score.cfd != cfd_concrete on %r/%r/%r" % (guide, proto, pamk))

    def test_bulge_aware_equals_cfd_concrete(self):
        """cfd_score_bulge_aware (production) == cfd_concrete on the same alignment
        including bulges -- both clean-room, one shared column_factor."""
        mm, pam = cfd_tables()
        gb = "ACGTACGTACGTACGTACGT"
        for pos in range(1, 20):
            # DNA bulge
            cr = gb[:pos] + "-" + gb[pos:] + "NGG"
            pr = gb[:pos] + "A" + gb[pos:] + "AGG"
            want = bp.cfd_score_bulge_aware(cr, pr, pam_len=3)
            got = tp.cfd_concrete(cr, pr[:-3], pr[-2:], mm, pam)
            self.assertEqual(_bits(got), _bits(want), "DNA-bulge col%d" % (pos + 1))


class CfdDoenchAnchors(unittest.TestCase):
    """Published Doench 2016 values the clean-room must carry exactly."""

    def test_table_anchor_values(self):
        mm, pam = cfd_tables()
        self.assertEqual(mm["rC:d-,2"], 0.96875)
        self.assertEqual(mm["r-:dA,2"], 0.692307692308)
        self.assertEqual(mm["r-:dA,20"], 0.6)
        self.assertEqual(pam["GG"], 1.0)
        self.assertEqual(len(mm), 392)   # 240 single-mismatch + 152 gap
        self.assertEqual(len(pam), 16)

    def test_perfect_match_is_pam_factor(self):
        g = "ACGTACGTACGTACGTACGT"
        self.assertEqual(cleanroom_cfd(g, g, "GG", cfd_tables()), 1.0)

    def test_noncanonical_pam_zeroes(self):
        g = "ACGTACGTACGTACGTACGT"
        self.assertEqual(cleanroom_cfd(g, g, "ZZ", cfd_tables()), 0.0)

    def test_iupac_protospacer_zeroes(self):
        g = "ACGTACGTACGTACGTACGT"
        self.assertEqual(cleanroom_cfd(g, g[:7] + "N" + g[8:], "GG", cfd_tables()), 0.0)


if __name__ == "__main__":
    unittest.main()
