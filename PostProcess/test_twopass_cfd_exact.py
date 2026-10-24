#!/usr/bin/env python3
"""Equivalence tests for twopass_cfd_exact (2.5.1 fast-mode exact CFD).

Two guarantees (design docs/DESIGN_2.5.1_two_pass_fast_mode.md §7):
  (1) cfd_concrete equals the clean-room CFD oracle: on mismatch-only alignments it is
      BIT-identical to cfd_score.cfd (raw-double ==, not abs<eps), and on bulge-bearing
      alignments it matches the live bulge-aware scorer (cfd_score_bulge_aware).
  (2) cfd_worst_case (factorized max over per-position allele sets) EQUALS the
      brute-force max over itertools.product of concrete haplotypes × the PAM set.

The CFD tables and oracle are the clean-room scorer (cfd_score / twopass_cfd_exact),
loaded from the committed JSON tables. STDLIB only.
"""
import itertools
import os
import random
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import twopass_cfd_exact as tp  # noqa: E402
import twopass_emit as te  # noqa: E402
import bulge_cfd_penalty as bp  # noqa: E402  (the live bulge-aware scorer, 2.8.0)
import cfd_score as cs  # noqa: E402  (the clean-room CFD oracle)


# Clean-room CFD tables (JSON), shared by every scorer under test.
MM, PAM = tp.load_scores()

_BASES = "ACGT"
_DNA_POOL = "ACGTACGT-N"          # weighted toward ACGT, incl. a bulge gap + an N
_PAMS = [a + b for a in "ACGTN" for b in "ACGTN"]


def _live(guide, dna, pam2):
    """The live bulge-aware scorer (cfd_score_bulge_aware) on the same alignment, so the
    matrix cfd_concrete can be pinned to it. Builds the full display alignment: a dummy
    3-nt PAM whose last 2 bases == the matrix pam, and a dummy guide PAM (ignored)."""
    return bp.cfd_score_bulge_aware(guide + "NNN", dna + "N" + pam2, MM, PAM, pam_len=3)


class TestCfdConcreteBitIdentical(unittest.TestCase):
    def test_random_matches_live_scorer(self):
        # cfd_concrete must equal the live bulge-aware scorer (cfd_score_bulge_aware) on EVERY
        # alignment, and additionally be BIT-identical to the clean-room cfd_score.cfd on the
        # mismatch-only subset (no gaps anywhere). The DNA pool includes a bulge gap + an N.
        random.seed(1)
        for _ in range(20000):
            guide = "".join(random.choice(_BASES) for _ in range(20))
            # inject occasional DNA bulges (guide '-') too, not just RNA bulges
            if random.random() < 0.15:
                k = random.randrange(20)
                guide = guide[:k] + "-" + guide[k + 1:]
            dna = "".join(random.choice(_DNA_POOL) for _ in range(20))
            pam = random.choice(_PAMS)
            ours = tp.cfd_concrete(guide, dna, pam, MM, PAM)
            self.assertAlmostEqual(ours, _live(guide, dna, pam), places=12,
                                   msg=f"guide={guide} dna={dna} pam={pam}")
            if "-" not in guide and "-" not in dna:   # mismatch-only -> raw-double identical
                self.assertEqual(ours, cs.cfd(guide, dna, pam, (MM, PAM)),
                                 f"mm-only {guide}/{dna}/{pam}")

    def test_edge_cases(self):
        g = "ACGTACGTACGTACGTACGT"
        # mismatch-only / PAM edges -> bit-identical to the clean-room cfd_score.cfd
        mm_only = [
            (g, g, "GG"),                              # all match
            (g, "N" + g[1:], "GG"),                    # untranslatable base -> 0
            (g, g, "NG"),                              # N in PAM
            (g, g, "ZZ"),                              # non-canonical PAM -> 0
            ("TTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAA", "GG"),  # T->U, all mm
        ]
        for guide, dna, pam in mm_only:
            self.assertEqual(tp.cfd_concrete(guide, dna, pam, MM, PAM),
                             cs.cfd(guide, dna, pam, (MM, PAM)),
                             f"edge {guide}/{dna}/{pam}")
        # bulge rows -> position-aware; pinned to the live scorer (NOT the old free-pass)
        bulge = [
            (g, "-" + g[1:], "GG"),                    # RNA bulge pos1 (terminal -> 0.9, not free)
            (g, g[:8] + "-" + g[9:], "GG"),            # RNA bulge mid (in-range)
            ("-" + g[1:], g, "GG"),                    # DNA bulge pos1 (terminal 5' tip)
            (g[:19] + "-", g, "GG"),                   # DNA bulge pos20 (seed, heavy)
        ]
        for guide, dna, pam in bulge:
            self.assertAlmostEqual(tp.cfd_concrete(guide, dna, pam, MM, PAM),
                                   _live(guide, dna, pam), places=12,
                                   msg=f"bulge edge {guide}/{dna}/{pam}")


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


class TestGreedyMaxScoreVsFactorizedOracle(unittest.TestCase):
    """2.5.1 item (c): the PRODUCTION worst-case-CFD path
    (``twopass_emit.greedy_max_score`` -> bounded brute-force with the clean-room
    cfd_concrete as score_fn) must agree with the INDEPENDENT factorized oracle
    (``twopass_cfd_exact.cfd_worst_case``). Two separately-written EXACT implementations
    returning the same worst-case CFD -- including the JOINT-PAM case that breaks a
    per-column greedy (the PAM factor is a 2-base joint lookup). Windows are kept small so
    the combo product stays <= greedy_max_score's brute cap (exact path, not the dense
    greedy fallback). PAM sets are built as a per-column PRODUCT so the same allele space
    is representable both as greedy_max_score PAM COLUMNS and as the oracle's ``pam_set``."""

    def _one(self, guide, proto_sets, pam_c1, pam_c2):
        L = len(guide)
        ref = [sorted(s)[0] for s in proto_sets] + [sorted(pam_c1)[0], sorted(pam_c2)[0]]
        ref_str = "".join(ref)
        columns = []
        # protospacer variant columns
        for i, s in enumerate(proto_sets):
            alts = [b for b in sorted(s) if b != ref[i]]
            if alts:
                columns.append({"pos_c": i, "candidates":
                                [{"alt": b, "carriers": {"S"}, "info": ["rs", "0.1", "snp"]}
                                 for b in alts]})
        # the two CFD-relevant PAM columns (positions L, L+1 = the last 2 PAM bases)
        for j, s in ((0, pam_c1), (1, pam_c2)):
            pc = L + j
            alts = [b for b in sorted(s) if b != ref[pc]]
            if alts:
                columns.append({"pos_c": pc, "candidates":
                                [{"alt": b, "carriers": {"S"}, "info": ["rs", "0.1", "snp"]}
                                 for b in alts]})

        def score_fn(seq_list):
            s = "".join(seq_list)
            return tp.cfd_concrete(guide, s[:L], s[L:L + 2], MM, PAM)

        rep = te.greedy_max_score(columns, ref_str, ref_str, guide, False, 0, L,
                                  tp._revcom, score_fn)
        rs = "".join(rep["seq"])
        got = tp.cfd_concrete(guide, rs[:L], rs[L:L + 2], MM, PAM)
        pam_set = sorted({a + b for a in pam_c1 for b in pam_c2})
        oracle = tp.cfd_worst_case(guide, [sorted(s) for s in proto_sets], pam_set, MM, PAM)
        # combos stayed under the brute cap (else greedy_max_score used the inexact
        # fallback and the comparison would not be meaningful).
        total = 1
        for c in columns:
            total *= (1 + len(c["candidates"]))
        return got, oracle, total

    def test_random_matches_oracle(self):
        random.seed(11)
        cap = 1 << 12
        for _ in range(4000):
            L = random.randint(4, 6)                      # 3^6 * 2 * 2 = 2916 <= 4096 cap
            guide = "".join(random.choice("ACGT") for _ in range(L))
            proto_sets = [set(random.choice("ACGT") for _ in range(random.randint(1, 3)))
                          for _ in range(L)]
            pam_c1 = set(random.choice("ACGT") for _ in range(random.randint(1, 2)))
            pam_c2 = set(random.choice("ACGT") for _ in range(random.randint(1, 2)))
            got, oracle, total = self._one(guide, proto_sets, pam_c1, pam_c2)
            self.assertLessEqual(total, cap, "window exceeded brute cap (test setup bug)")
            self.assertEqual(got, oracle,
                             "greedy_max_score vs factorized oracle: %r != %r "
                             "(guide=%s proto=%s pam=%s/%s)"
                             % (got, oracle, guide, proto_sets, pam_c1, pam_c2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
