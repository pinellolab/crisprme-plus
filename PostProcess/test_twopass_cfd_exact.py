#!/usr/bin/env python3
"""Equivalence tests for twopass_cfd_exact (2.5.1 fast-mode exact CFD).

Two guarantees (design docs/DESIGN_2.5.1_two_pass_fast_mode.md §7):
  (1) cfd_concrete is BIT-identical to the shipped new_simple_analysis.calc_cfd
      (raw-double ==, not abs<eps) over random + edge inputs.
  (2) cfd_worst_case (factorized max over per-position allele sets) EQUALS the
      brute-force max over itertools.product of concrete haplotypes × the PAM set.

The real calc_cfd / revcom / get_mm_pam_scores are loaded by exec'ing ONLY their
`def` blocks from new_simple_analysis.py (that module runs analysis at import, so
it cannot be imported directly). STDLIB only.
"""
import itertools
import os
import pickle
import random
import re
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import twopass_cfd_exact as tp  # noqa: E402


def _extract_defs(src, names):
    """Return the source of the named top-level `def` blocks, concatenated."""
    lines = src.splitlines(keepends=True)
    out, i = [], 0
    while i < len(lines):
        m = re.match(r"def (\w+)\(", lines[i])
        if m and m.group(1) in names:
            j = i + 1
            while j < len(lines) and (lines[j][:1] in (" ", "\t") or not lines[j].strip()):
                j += 1
            out.append("".join(lines[i:j]))
            i = j
        else:
            i += 1
    return "".join(out)


_NSA = os.path.join(_HERE, "new_simple_analysis.py")
_ref = {"os": os, "pickle": pickle, "__file__": _NSA}
exec(_extract_defs(open(_NSA).read(), {"revcom", "calc_cfd", "get_mm_pam_scores"}), _ref)
_calc_cfd = _ref["calc_cfd"]
MM, PAM = _ref["get_mm_pam_scores"]()

_BASES = "ACGT"
_DNA_POOL = "ACGTACGT-N"          # weighted toward ACGT, incl. a bulge gap + an N
_PAMS = [a + b for a in "ACGTN" for b in "ACGTN"]


class TestCfdConcreteBitIdentical(unittest.TestCase):
    def test_random_bit_identical(self):
        random.seed(1)
        for _ in range(20000):
            guide = "".join(random.choice(_BASES) for _ in range(20))
            dna = "".join(random.choice(_DNA_POOL) for _ in range(20))
            pam = random.choice(_PAMS)
            ours = tp.cfd_concrete(guide, dna, pam, MM, PAM)
            ref = _calc_cfd(guide, dna, pam, MM, PAM, True)
            self.assertEqual(ours, ref,
                             f"guide={guide} dna={dna} pam={pam}: {ours!r} != {ref!r}")

    def test_edge_cases(self):
        g = "ACGTACGTACGTACGTACGT"
        cases = [
            (g, g, "GG"),                              # all match
            (g, "-" + g[1:], "GG"),                    # position-1 gap (free)
            (g, g[:8] + "-" + g[9:], "GG"),            # mid gap
            (g, "N" + g[1:], "GG"),                    # untranslatable base -> 0
            (g, g, "NG"),                              # N in PAM
            (g, g, "ZZ"),                              # non-canonical PAM -> 0
            ("TTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAA", "GG"),  # T->U, all mm
        ]
        for guide, dna, pam in cases:
            self.assertEqual(tp.cfd_concrete(guide, dna, pam, MM, PAM),
                             _calc_cfd(guide, dna, pam, MM, PAM, True),
                             f"edge {guide}/{dna}/{pam}")


class TestWorstCaseFactorizationOracle(unittest.TestCase):
    def test_matches_bruteforce(self):
        # short guides so itertools.product is tractable; factorization is
        # length-independent, so this fully exercises the max = product-of-maxes law.
        random.seed(2)
        for _ in range(2000):
            L = random.randint(4, 8)
            guide = "".join(random.choice(_BASES) for _ in range(L))
            allele_sets = [sorted(set(random.choice("ACGT-")
                                      for _ in range(random.randint(1, 3))))
                           for _ in range(L)]
            pam_set = sorted(set(random.choice(_PAMS)
                                 for _ in range(random.randint(1, 3))))
            ours = tp.cfd_worst_case(guide, allele_sets, pam_set, MM, PAM)
            best = 0.0
            for combo in itertools.product(*allele_sets):
                dna = "".join(combo)
                for p in pam_set:
                    v = tp.cfd_concrete(guide, dna, p, MM, PAM)
                    if v > best:
                        best = v
            self.assertEqual(ours, best,
                             f"worst-case {guide} sets={allele_sets} pam={pam_set}: "
                             f"{ours!r} != brute {best!r}")

    def test_full_ambiguity_perfect_match(self):
        # a fully-ambiguous window (every base possible) -> the guide is formable
        # exactly -> worst-case CFD = max PAM (the design's "NNNN is easy" case).
        guide = "ACGTACGTACGTACGTACGT"
        allele_sets = [list("ACGT") for _ in range(20)]
        got = tp.cfd_worst_case(guide, allele_sets, ["GG"], MM, PAM)
        self.assertEqual(got, PAM.get("GG", 0.0))


if __name__ == "__main__":
    unittest.main(verbosity=2)
