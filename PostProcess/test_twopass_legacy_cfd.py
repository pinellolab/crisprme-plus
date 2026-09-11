"""Legacy (dict / mega) fast-mode WORST-CASE CFD equivalence -- 2.5.1 item (b).

The legacy 2^k / dict path (``new_simple_analysis.iupac_decomposition`` with ``mygt``
None) in fast mode used to emit ONLY the min-mismatch greedy representative. But CFD is a
position-WEIGHTED product, so at a MULTIALLELIC forced-mismatch column the min-mismatch
(lex-tiebreak) allele can score LOWER CFD than another carried allele -> the fast row
UNDER-STATES the worst-case CFD a dict/mega search should report. The fix emits a second
(max-CFD) representative on this path too (mirrors the observed path + the ml007 fast-vs-
slow CFD benchmark).

This test drives the REAL legacy ``iupac_decomposition`` both ways (slow enumerate / fast)
with REAL CFD matrices injected, on a fixture ADAPTIVELY chosen so the min-mm allele is NOT
the max-CFD allele, and asserts fast's worst-case CFD >= slow's (== the true worst case).
Without the fix this fails (fast would carry only the lower-CFD lex allele). STDLIB +
the harness's numpy/pandas/CRISTA stubs + the real mismatch_score/PAM_scores pickles.

Run: cd PostProcess && python3 -m unittest test_twopass_legacy_cfd -v
"""
import io
import os
import pickle
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN
from test_observed_haplotype_integration import _variant_rows

_MM = pickle.load(open(os.path.join(PP, "mismatch_score.pkl"), "rb"))
_PAM = pickle.load(open(os.path.join(PP, "PAM_scores.pkl"), "rb"))

# calc_cfd is a pure function of its args; borrow the real one from the module to score
# rows exactly as preprocess_CFD_score does (guide-with-PAM, protospacer, last-2 PAM).
_NS0 = _load_pure_functions({})
_calc_cfd = _NS0["calc_cfd"]


def _score_row_cfd(final_line):
    """CFD of a finalized variant row, exactly as preprocess_CFD_score computes it."""
    t = final_line[2].upper()
    return _calc_cfd(GUIDE, t[:-3], t[-2:], _MM, _PAM, True)


def _score_single_mm(pos, alt):
    """CFD of a single A>alt mismatch at protospacer ``pos`` (NGG PAM = ...AGG)."""
    proto = list("A" * 20)
    proto[pos] = alt
    target23 = "".join(proto) + "AGG"
    return _calc_cfd(GUIDE, target23[:-3], target23[-2:], _MM, _PAM, True)


class LegacyFastWorstCaseCFD(unittest.TestCase):
    def _pick_divergent(self):
        """First (pos, a1, a2) with a1<a2 (lex) but CFD(a1)<CFD(a2): the min-mm tiebreak
        picks a1 while the true worst case is a2. Guaranteed to exist for real CFD weights
        (they are not reverse-lex at every position)."""
        for pos in range(20):
            cfd = {a: _score_single_mm(pos, a) for a in "CGT"}
            for a1 in "CGT":
                for a2 in "CGT":
                    if a1 < a2 and cfd[a1] < cfd[a2] - 1e-9:
                        return pos, a1, a2
        self.fail("no divergent position found (unexpected for real CFD matrices)")

    def _run(self, pos, a1, a2, fast):
        ref = ["A"] * TARGET_LEN
        ref[pos] = "C"                      # SNP column ref = mismatch vs 'A' guide
        ref[20], ref[21], ref[22] = "A", "G", "G"
        genome = "".join(ref)
        dna = list(genome)
        dna[pos] = "B"                      # IUPAC {C,G,T} -> variant-decomposition branch
        dna = "".join(dna)
        split = [
            "X", GUIDE, dna, "chrT", "0", "0", "+",
            "1", "0", "1", "NGG", "y", "NA", "NA", "NA", "NA", "NA", "NA",
        ]
        # MULTIALLELIC dict entry: both alts are forced mismatches, carried by distinct
        # samples so BOTH are emittable (a1 by S1, a2 by S2).
        mydict = {
            "chrT," + str(pos + 1):
                "S1:1|0;C,%s;rs1;0.10$S2:1|0;C,%s;rs2;0.10" % (a1, a2)
        }
        overrides = dict(
            genomeStr=genome, current_chr="chrT", mydict=mydict,
            myreg=None, mygt=None,          # mygt None -> LEGACY 2^k / dict path
            haplotype_check=True, IUPAC_CAP=10, hvdr_bed=io.StringIO(),
            pam="NGG", pos_beg=0, pos_end=-3, pam_begin=-3, pam_end=None,
            allowed_mms=6, bulge_pos=8,     # split[8] = Bulge_Size (module global)
            do_scores=True, mm_scores=_MM, pam_scores=_PAM,   # ACTIVATE CFD scoring
        )
        if fast:
            overrides["_FAST_MODE"] = True
        ns = _load_pure_functions(overrides)
        ns["_phase_confirmation_rows"].clear()
        ns["_phase_confirmation_keys"].clear()
        cluster = []
        ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
        return cluster

    def test_fast_reports_worst_case_cfd(self):
        pos, a1, a2 = self._pick_divergent()
        slow = _variant_rows(self._run(pos, a1, a2, fast=False))
        fast = _variant_rows(self._run(pos, a1, a2, fast=True))

        self.assertTrue(slow, "slow path produced no variant off-target")
        self.assertTrue(fast, "fast path produced no variant off-target")

        slow_max = max(_score_row_cfd(r) for r in slow)
        fast_max = max(_score_row_cfd(r) for r in fast)
        cfd_a1, cfd_a2 = _score_single_mm(pos, a1), _score_single_mm(pos, a2)

        # sanity: the fixture actually exercises the gap (min-mm allele under-scores).
        self.assertLess(cfd_a1, cfd_a2,
                        "fixture not divergent: min-mm allele already max-CFD")
        # slow (full enumeration) surfaces the true worst case = the higher-CFD allele.
        self.assertAlmostEqual(slow_max, cfd_a2, places=6)
        # THE PROPERTY (b) guarantees: fast never under-states the worst-case CFD.
        # Without the max-CFD rep, fast would carry only the lex/min-mm allele a1 and
        # fast_max would equal cfd_a1 < cfd_a2 -> this assertion would FAIL.
        self.assertGreaterEqual(fast_max, slow_max - 1e-9,
                                "fast UNDER-reports worst-case CFD (pos=%d %s/%s): "
                                "fast=%.4f slow=%.4f" % (pos, a1, a2, fast_max, slow_max))

    def test_fast_still_collapses_to_at_most_two_reps(self):
        # worst-possible = at most min-mm + max-CFD (never the full 2^k lattice).
        pos, a1, a2 = self._pick_divergent()
        fast = _variant_rows(self._run(pos, a1, a2, fast=True))
        self.assertLessEqual(len(fast), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
