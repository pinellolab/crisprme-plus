"""Fast-mode (2.5.1 two-pass) end-to-end equivalence against the slow enumeration.

Drives the REAL ``new_simple_analysis.iupac_decomposition`` through the AST-load harness
(``test_phased_haplotype._load_pure_functions``) both ways -- default (enumerate) and
``_FAST_MODE=True`` (single worst-possible representative) -- on the SAME fixtures, and
asserts the two load-bearing properties from docs/DESIGN_2.5.1_two_pass_fast_mode.md:

  1. LOSSLESS DETECTION: fast emits a variant off-target at EVERY locus the slow path
     does (a superset -- worst-possible may also surface a formable-but-unobserved
     combination), and always emits the candidate's reference off-target.
  2. WORST-POSSIBLE: the fast representative's mismatch count is <= the MINIMUM mismatch
     over all slow (observed / enumerated) haplotypes at that locus -- i.e. fast never
     under-states the worst-case off-target (it is the argmin edit over all allele combos).

Both the dict-less OBSERVED path (mygt present) and the legacy 2^k path (mygt None) are
exercised. STDLIB only (+ the harness's numpy/pandas/CRISTA stubs).

Run: cd PostProcess && python3 -m unittest test_twopass_fast_equivalence -v
"""
import io
import os
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN
from test_observed_haplotype_integration import (
    _build, _GtSentinel, _is_reference_row, _variant_rows, _reference_rows,
)


def _run(n_snps, fast, legacy=False, **kw):
    """Drive iupac_decomposition once. ``fast`` toggles _FAST_MODE; ``legacy`` drops the
    genotype tier (mygt None) to take the legacy 2^k path instead of the observed one."""
    split, overrides = _build(n_snps, **kw)
    if legacy:
        overrides["mygt"] = None            # -> legacy 2^k / greedy-cap path
    if fast:
        overrides["_FAST_MODE"] = True      # applied AFTER exec, so it wins over the env
    ns = _load_pure_functions(overrides)
    ns["_phase_confirmation_rows"].clear()
    ns["_phase_confirmation_keys"].clear()
    cluster = []
    ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
    return cluster


def _min_variant_mm(rows):
    var = _variant_rows(rows)
    return min((int(r[7]) for r in var), default=None)


class _EquivalenceMixin(object):
    legacy = False

    def _assert_equiv(self, n_snps, **kw):
        slow = _run(n_snps, fast=False, legacy=self.legacy, **kw)
        fast = _run(n_snps, fast=True, legacy=self.legacy, **kw)

        slow_var = _variant_rows(slow)
        fast_var = _variant_rows(fast)

        # (1) fast collapses the enumeration to ONE representative (never MORE rows).
        self.assertLessEqual(len(fast_var), max(1, len(slow_var)) if slow_var else 1)

        # (2) lossless DETECTION: if the slow path found any in-budget variant off-target
        # at this locus, fast must also emit one.
        if slow_var:
            self.assertTrue(fast_var,
                            "fast dropped a locus the slow path detected (n=%d)" % n_snps)

        # (3) worst-POSSIBLE: fast's mismatch <= the min over all slow haplotypes.
        if slow_var and fast_var:
            self.assertLessEqual(
                _min_variant_mm(fast), _min_variant_mm(slow),
                "fast representative understates the worst case (n=%d): fast=%r slow=%r"
                % (n_snps, _min_variant_mm(fast), _min_variant_mm(slow)))

        # (4) the reference off-target is emitted by BOTH (locus coverage preserved).
        self.assertEqual(len(_reference_rows(fast)), len(_reference_rows(slow)))
        return slow, fast


class ObservedPathEquivalence(_EquivalenceMixin, unittest.TestCase):
    legacy = False

    def test_phased_cis(self):
        for n in (1, 2, 3, 4, 6):
            self._assert_equiv(n, phased=True)

    def test_phased_trans(self):
        # trans: no individual carries the full cis combo, so slow emits sub-haplotypes
        # only. Fast (worst-possible) may surface the formable full combo -- still a
        # superset with <= mismatch. Detection + worst-case bound must hold.
        for n in (2, 4, 6):
            self._assert_equiv(n, phased=True, trans=True)

    def test_unphased(self):
        for n in (1, 2, 3, 4):
            self._assert_equiv(n, phased=False)

    def test_unphased_n4_collapses_15_to_1(self):
        slow, fast = self._assert_equiv(4, phased=False)
        # slow enumerates all 2^4-1 subsets; fast emits exactly one 0-mm representative.
        self.assertEqual(len(_variant_rows(slow)), 15)
        self.assertEqual(len(_variant_rows(fast)), 1)
        rep = _variant_rows(fast)[0]
        self.assertEqual(int(rep[7]), 0)          # the worst case is a perfect match
        self.assertEqual(rep[12], "S")            # its carrier is surfaced


class LegacyPathEquivalence(_EquivalenceMixin, unittest.TestCase):
    legacy = True

    def test_phased_cis(self):
        for n in (1, 2, 3, 4, 6):
            self._assert_equiv(n, phased=True)

    def test_unphased(self):
        for n in (1, 2, 3, 4):
            self._assert_equiv(n, phased=False)

    def test_legacy_fast_collapses_enumeration(self):
        # dict-only legacy path, N=5 phased-cis: slow builds the multi-SNP lattice; fast
        # forces the single greedy representative (0-mm, carrier S).
        slow, fast = self._assert_equiv(5, phased=True)
        self.assertGreaterEqual(len(_variant_rows(slow)), len(_variant_rows(fast)))
        self.assertTrue(any(int(r[7]) == 0 for r in _variant_rows(fast)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
