#!/usr/bin/env python3
"""Regression test: the two calc_cfd twins must stay byte-for-byte equivalent.

CRISPRme+ carries two copies of the CFD scorer — the canonical SNP/variant path
new_simple_analysis.calc_cfd and the INDEL path analisi_indels_NNN.calc_cfd. A
correctness audit (docs/analyses/cfd_correctness_audit.md) found the INDEL copy
had drifted with two INDEL-only branches:

  * `if "N" == sl: score *= 1` in the mismatch loop — dead code (revcom("N") is
    None, so the key build TypeErrors and the score is zeroed+broken first), but
    confusing; and
  * `if "N" in pam: score *= 1 else: ...get(pam, 0.0)` in the PAM factor — a REAL
    divergence: an N in the scored 2bp PAM kept a non-zero CFD on the INDEL path
    while the SNP path zeroed it.

The INDEL copy was reconciled to the SNP rule. This test pins them together so
they cannot silently diverge again: it exec's ONLY the `revcom` + `calc_cfd` def
blocks from each module (both modules run analysis at import, so neither can be
imported directly) and asserts identical raw-double output over random + edge
inputs — mismatch-only, bulge gaps, IUPAC/N in the off-target base, and crucially
N and other non-canonical bases in the PAM. STDLIB only.
"""
import os
import pickle
import random
import re
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SNP = os.path.join(_HERE, "new_simple_analysis.py")
_INDEL = os.path.join(_HERE, "analisi_indels_NNN.py")


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


def _load_calc_cfd(path):
    ns = {"os": os, "pickle": pickle, "__file__": path}
    exec(_extract_defs(open(path).read(), {"revcom", "calc_cfd"}), ns)
    return ns["calc_cfd"]


_SNP_CFD = _load_calc_cfd(_SNP)
_INDEL_CFD = _load_calc_cfd(_INDEL)

MM = pickle.load(open(os.path.join(_HERE, "mismatch_score.pkl"), "rb"))
PAM = pickle.load(open(os.path.join(_HERE, "PAM_scores.pkl"), "rb"))

# off-target pools weighted toward real bases but including a bulge gap, an N, and
# an IUPAC code so the ambiguous-base + gap paths are exercised on both twins.
_GUIDE_POOL = "ACGT"
_OFF_POOL = "ACGTACGT-NR"
_PAM_BASES = "ACGTN"


class TestCalcCfdTwinsAgree(unittest.TestCase):
    def _assert_equal(self, guide, off, pam):
        a = _SNP_CFD(guide, off, pam, MM, PAM, True)
        b = _INDEL_CFD(guide, off, pam, MM, PAM, True)
        # raw-double equality, not abs<eps: the two implementations must be identical
        self.assertEqual(
            a, b, f"twin CFD mismatch: SNP={a!r} INDEL={b!r} "
                  f"(guide={guide!r} off={off!r} pam={pam!r})"
        )

    def test_random_agreement(self):
        random.seed(7)
        for _ in range(50000):
            guide = "".join(random.choice(_GUIDE_POOL) for _ in range(20))
            off = "".join(random.choice(_OFF_POOL) for _ in range(20))
            pam = "".join(random.choice(_PAM_BASES) for _ in range(2))
            self._assert_equal(guide, off, pam)

    def test_n_in_pam_all_positions(self):
        # the former real divergence: N anywhere in the scored 2bp PAM
        guide = "ACGTACGTACGTACGTACGT"
        off = "ACGTACGTACGTACGTACGT"  # perfect spacer match -> isolates the PAM factor
        for pam in ("NG", "GN", "NN", "NA", "AN", "NC", "CN", "NT", "TN"):
            self._assert_equal(guide, off, pam)
            # and it must zero, matching the SNP/issue-#94 rule
            self.assertEqual(_INDEL_CFD(guide, off, pam, MM, PAM, True), 0.0)

    def test_n_in_offtarget_base(self):
        # the former dead-code branch: N (and other ambiguous bases) in the spacer
        guide = "ACGTACGTACGTACGTACGT"
        for i in range(20):
            for amb in ("N", "R", "Y", "n"):
                off = guide[:i] + amb + guide[i + 1:]
                self._assert_equal(guide, off, "GG")
                # ambiguous off-target base -> CFD collapses to 0 on both twins
                self.assertEqual(_INDEL_CFD(guide, off, "GG", MM, PAM, True), 0.0)

    def test_bulge_gaps_each_position(self):
        guide = "ACGTACGTACGTACGTACGT"
        for i in range(20):
            off = guide[:i] + "-" + guide[i + 1:]
            self._assert_equal(guide, off, "GG")

    def test_canonical_and_noncanonical_pam(self):
        guide = "ACGTACGTACGTACGTACGT"
        off = "ACGTACGTACGTACGTACGT"
        for pam in ("GG", "AG", "TG", "CG", "GA", "AA", "XY", "--", "G-"):
            self._assert_equal(guide, off, pam)


if __name__ == "__main__":
    unittest.main(verbosity=2)
