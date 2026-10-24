"""Bit-identity of the clean-room CFD scorer (`cfd_score.cfd`) against an
independent clean-room witness.

The primary oracle is `twopass_cfd_exact.cfd_concrete`, the exact-factorization
scorer used by two-pass fast mode: a *separate* clean-room implementation of the
same Doench 2016 definition (one a direct per-position product, the other the
concrete evaluation of the factorized worst-case scorer). Agreement across the
mismatch / PAM / IUPAC / perfect-match case classes is a genuine
cross-implementation bit-identity check, not a tautology. Both scorers were
proven 0-ULP over 50k+ mismatch-only cases (see
`test_cfd_cleanroom_equiv.py`).

The one exception is the DNA/RNA gap (bulge) case class. `cfd_score.cfd` is
bulge-AGNOSTIC: a gap column is scored by a flat per-position gap-key lookup
(`r<spacer>:d<complement>,<pos>`, a missing key -> factor 1.0) -- exactly the
flat-product semantics of the removed legacy `calc_cfd`. `cfd_concrete` is
bulge-AWARE (it applies the 2.8.0 position-aware bulge factor and shifts the
guide frame across DNA-bulge columns), so by design it DISAGREES on gap columns
and is not a valid witness for them. That case class therefore keeps the
bulge-agnostic flat gap-key-lookup oracle; see `_assert_same_flat_gap`.

Equality is checked on the raw IEEE-754 bytes of the two doubles, so -0.0, a
one-ULP divergence or a denormal all fail loudly. The witness is used purely as
a black-box oracle here. Both oracles load the committed JSON tables (never a
pickle).

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
import twopass_cfd_exact as _oracle_mod

_TABLES = cfd_score.cfd_tables()
_ORACLE_MM, _ORACLE_PAM = _oracle_mod.load_scores()
_BASES = "ACGT"


def _oracle(spacer, protospacer, pam):
    return _oracle_mod.cfd_concrete(spacer, protospacer, pam, _ORACLE_MM, _ORACLE_PAM)


# Watson-Crick complement for the "d<base>" half of a gap/mismatch key; a gap maps
# to a gap, U collapses to A (RNA read). Written locally so the flat-gap oracle below
# carries no inherited scorer body.
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "-": "-"}


def _oracle_flat_gap(spacer, protospacer, pam):
    """Independent bulge-AGNOSTIC CFD oracle: the flat per-position gap-key lookup
    (the product of ``mm["r<spacer>:d<complement>,<pos>"]``, a missing key -> 1.0,
    times the PAM factor). Computed straight from the committed JSON tables with no
    inherited scorer body, so the gap/RNA-gap rows keep a real second witness
    instead of comparing ``cfd_score.cfd`` to itself. (``cfd_concrete`` is
    bulge-AWARE and intentionally diverges on gap columns, so it cannot serve here.)"""
    sp = spacer.replace("T", "U")
    pr = protospacer.replace("T", "U")
    product = 1.0
    for pos, (sb, pb) in enumerate(zip(sp, pr), start=1):
        if sb == pb:
            continue
        comp = _COMPLEMENT.get(pb)
        if comp is None:
            return 0.0
        product *= _ORACLE_MM.get("r%s:d%s,%d" % (sb, comp, pos), 1.0)
    return product * _ORACLE_PAM.get(pam, 0.0)


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

    def _assert_same_flat_gap(self, spacer, protospacer, pam):
        got = cfd_score.cfd(spacer, protospacer, pam, _TABLES)
        want = _oracle_flat_gap(spacer, protospacer, pam)
        self.assertEqual(
            _bits(got), _bits(want),
            "CFD (bulge-agnostic gap) differs for spacer=%s proto=%s pam=%s: %r vs oracle %r"
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
        # cfd_score.cfd is bulge-AGNOSTIC here (flat gap-key lookup, legacy semantics),
        # so the oracle is the independent flat-gap reference, not the bulge-aware witness.
        spacer = "ACGTACGTACGTACGTACGT"
        # DNA gap (d-): a '-' in the protospacer at each position
        for pos in range(20):
            proto = list(spacer)
            proto[pos] = "-"
            self._assert_same_flat_gap(spacer, "".join(proto), "GG")
        # RNA gap (r-): a '-' in the spacer at each position (protospacer keeps a base)
        for pos in range(20):
            sp = list(spacer)
            sp[pos] = "-"
            self._assert_same_flat_gap("".join(sp), spacer, "GG")

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
