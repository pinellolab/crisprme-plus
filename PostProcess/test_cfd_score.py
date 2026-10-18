"""Bit-identity of the clean-room CFD scorer (`cfd_score.cfd`) against the
reference implementation that produced every CRISPRme+ result to date
(`targets_scoring.calc_cfd`).

Equality is checked on the raw IEEE-754 bytes of the two doubles, so -0.0, a
one-ULP divergence or a denormal all fail loudly. The reference is used purely
as a black-box oracle here.

Run:
    cd PostProcess && python3 -m unittest test_cfd_score -v
"""
import itertools
import os
import random
import struct
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cfd_score
import targets_scoring as _oracle_mod

_ORACLE_MM, _ORACLE_PAM = _oracle_mod.get_mm_pam_scores()
_TABLES = cfd_score.cfd_tables()
_BASES = "ACGT"


def _oracle(spacer, protospacer, pam):
    return _oracle_mod.calc_cfd(spacer, protospacer, pam, _ORACLE_MM, _ORACLE_PAM, True)


def _bits(x):
    return struct.pack("<d", float(x))


class TestCfdBitIdentity(unittest.TestCase):
    def _assert_same(self, spacer, protospacer, pam):
        got = cfd_score.cfd(spacer, protospacer, pam, _TABLES)
        want = _oracle(spacer, protospacer, pam)
        self.assertEqual(
            _bits(got), _bits(want),
            "CFD differs for spacer=%s proto=%s pam=%s: %r vs oracle %r"
            % (spacer, protospacer, pam, got, want),
        )

    def test_exhaustive_single_mismatch_all_positions(self):
        # every guide base x every position (1..20) x every substitution, PAM=GG
        for gbase in _BASES:
            spacer = gbase * 20
            for pos in range(20):
                for tbase in _BASES:
                    proto = list(spacer)
                    proto[pos] = tbase
                    self._assert_same(spacer, "".join(proto), "GG")

    def test_all_pam_keys_including_non_canonical(self):
        spacer = "ACGTACGTACGTACGTACGT"
        for a, b in itertools.product(_BASES, repeat=2):
            # perfect protospacer so only the PAM factor varies
            self._assert_same(spacer, spacer, a + b)

    def test_random_multi_mismatch(self):
        rng = random.Random(20240617)
        for _ in range(20000):
            spacer = "".join(rng.choice(_BASES) for _ in range(20))
            proto = list(spacer)
            for pos in rng.sample(range(20), rng.randint(1, 6)):
                proto[pos] = rng.choice(_BASES)
            pam = rng.choice(_BASES) + rng.choice(_BASES)
            self._assert_same(spacer, "".join(proto), pam)

    def test_gap_keys_dna_and_rna_bulge(self):
        spacer = "ACGTACGTACGTACGTACGT"
        # DNA gap (d-): a '-' in the protospacer at each position
        for pos in range(20):
            proto = list(spacer)
            proto[pos] = "-"
            self._assert_same(spacer, "".join(proto), "GG")
        # RNA gap (r-): a '-' in the spacer at each position (protospacer keeps a base)
        for pos in range(20):
            sp = list(spacer)
            sp[pos] = "-"
            self._assert_same("".join(sp), spacer, "GG")

    def test_ambiguous_protospacer_base_zeroes_score(self):
        spacer = "ACGTACGTACGTACGTACGT"
        for pos in (0, 7, 19):
            proto = list(spacer)
            proto[pos] = "N"
            self._assert_same(spacer, "".join(proto), "GG")

    def test_perfect_match_is_pam_factor(self):
        spacer = "ACGTACGTACGTACGTACGT"
        self._assert_same(spacer, spacer, "GG")  # == pam['GG'] == 1.0


if __name__ == "__main__":
    unittest.main()
