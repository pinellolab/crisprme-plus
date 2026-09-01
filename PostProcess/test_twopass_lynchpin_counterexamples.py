#!/usr/bin/env python3
"""Design-spec regression fixtures for the #175 two-pass fast-mode LOWER-BOUND lynchpin.

This is NOT testing shipped code — it is an executable specification of the two
adversarially-VERIFIED counterexamples where the collapsed-IUPAC min-edit D
(min over {reference} + {one fake contig per single indel}) is NOT a lower bound
on the true edit distance d(guide, H) of a real multi-indel haplotype H, and a
provably-safe correction that closes the gap. It exists so that whenever the fast
mode is built, these cases are pinned as regressions (and so the gap is documented
in-repo, not just in a workflow transcript).

BACKGROUND (adversarial panel, verifier-confirmed):
  The min-edit lower bound HOLDS for SNPs (per-position set-containment) and for
  bulges (reuse H's own alignment), but BREAKS when >=2 germline indels co-occur
  in cis within one protospacer: the build materializes ONE indel per fake contig,
  so a haplotype carrying two indels lives on NO single contig, and D over-counts
  by up to the total bulge cost of the un-modeled indels -> a FALSE NEGATIVE, even
  at k=0 (the safety-critical perfect-match case).

  A third verified case (equal-length MNV / block substitution silently dropped by
  the enricher) is documented at the bottom; its fix is `bcftools norm -a`
  atomization, a different mechanism, so it is described but not asserted here.

STDLIB ONLY. Run: python3 PostProcess/test_twopass_lynchpin_counterexamples.py
"""

import unittest


def levenshtein(a, b):
    """Unit-cost edit distance (a lower bound on CRISPRme mm+bulge cost)."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def collapsed_min_edit(guide, templates):
    """The fast-mode Step-1 value D = min edit over the enumerated representations
    (the reference + one fake contig per SINGLE indel). This is what the current
    build can express: never a contig with two indels applied at once."""
    return min(levenshtein(guide, t) for t in templates)


# Each counterexample gives PRE-VERIFIED sequences (indices independently checked):
#   guide, the enumerated single-indel templates {ref, contig_indelA, contig_indelB},
#   the REAL cis haplotype H (both indels applied), and the # of cis indels in window.
COUNTEREXAMPLES = [
    {
        "name": "length-cancelling cis indel pair (+1 ins / -1 del) -> perfect off-target",
        "guide": "ATGCGGTATTCGACAGAGCA",
        "templates": {
            "reference":        "ATGCGGTATTGACAGAGCTA",   # no indel applied  (d=2)
            "contig_ins_only":  "ATGCGGTATTCGACAGAGCTA",   # +C only, 21nt     (d=1)
            "contig_del_only":  "ATGCGGTATTGACAGAGCA",     # -T only, 19nt     (d=1)
        },
        "H": "ATGCGGTATTCGACAGAGCA",   # both -> length-cancels -> == guide (d=0)
        "n_cis_indels": 2,
        "indel_bulges": [1, 1],        # |len(ref)-len(alt)| per indel
    },
    {
        "name": "two 2bp cis deletions -> perfect off-target",
        "guide": "ACGTACGTACGTACGTACGT",
        "templates": {
            "reference":        "ACGTGGACGTACGTTTACGTACGT",   # GG+TT present    (d=4)
            "contig_delGG":     "ACGTACGTACGTTTACGTACGT",     # -GG only         (d=2)
            "contig_delTT":     "ACGTGGACGTACGTACGTACGT",     # -TT only         (d=2)
        },
        "H": "ACGTACGTACGTACGTACGT",   # both deletions -> == guide (d=0)
        "n_cis_indels": 2,
        "indel_bulges": [2, 2],
    },
]


class TestMinEditLowerBoundBreaksOnMultiIndel(unittest.TestCase):
    def test_naive_collapsed_D_is_not_a_lower_bound(self):
        """The core break: D > d(guide,H) — the fast mode would MISS a real
        perfect (0-edit) off-target at k=0."""
        for ce in COUNTEREXAMPLES:
            g = ce["guide"]
            d_true = levenshtein(g, ce["H"])
            D = collapsed_min_edit(g, list(ce["templates"].values()))
            self.assertEqual(d_true, 0, ce["name"] + ": H must be a perfect match")
            self.assertGreater(
                D, d_true,
                ce["name"] + ": expected the miss D > d (lower-bound violation)")
            # concretely: thresholding D <= 0 (perfect-match detection) drops a
            # site whose real edit distance IS 0 -> a false negative.
            self.assertFalse(D <= 0, ce["name"] + ": naive threshold wrongly excludes")

    def test_count_only_slack_is_INSUFFICIENT_for_multibp_indels(self):
        """The panel's first-cut slack (subtract n_indels-1) closes the 1bp case
        but NOT the 2bp case — a real refinement this fixture pins."""
        closed = {}
        for ce in COUNTEREXAMPLES:
            g = ce["guide"]
            D = collapsed_min_edit(g, list(ce["templates"].values()))
            D_slack = max(0, D - (ce["n_cis_indels"] - 1))
            closed[ce["name"]] = (D_slack <= levenshtein(g, ce["H"]))
        self.assertTrue(closed[COUNTEREXAMPLES[0]["name"]],
                        "count-slack should close the 1bp cancelling pair")
        self.assertFalse(closed[COUNTEREXAMPLES[1]["name"]],
                         "count-slack must FAIL to close the 2bp-deletion pair "
                         "(over-count scales with indel SIZE, not count)")

    def test_size_aware_slack_restores_the_lower_bound(self):
        """A size-aware slack (subtract the worst un-modeled indel bulge budget =
        sum of bulges minus the smallest, since a single contig models one indel)
        makes D' <= d for BOTH cases -> lossless again."""
        for ce in COUNTEREXAMPLES:
            g = ce["guide"]
            D = collapsed_min_edit(g, list(ce["templates"].values()))
            b = sorted(ce["indel_bulges"])
            unmodeled = sum(b) - b[0]        # worst case: chosen contig modeled the smallest
            D_prime = max(0, D - unmodeled)
            self.assertLessEqual(
                D_prime, levenshtein(g, ce["H"]),
                ce["name"] + ": size-aware slack must not miss")

    def test_flag_all_multiindel_windows_is_provably_lossless(self):
        """The maximally-conservative fix — flag any window with >=2 cis indels
        regardless of D — trivially never misses. Empirically cheap because such
        windows are ~entirely STR/low-complexity loci already flagged as high-
        variant-density regions."""
        k = 0
        for ce in COUNTEREXAMPLES:
            g = ce["guide"]
            D = collapsed_min_edit(g, list(ce["templates"].values()))
            naive_flag = (D <= k)
            safe_flag = (D <= k) or (ce["n_cis_indels"] >= 2)
            self.assertFalse(naive_flag, ce["name"] + ": naive misses")
            self.assertTrue(safe_flag, ce["name"] + ": flag-all-multiindel catches it")


# ---------------------------------------------------------------------------
# THIRD verified counterexample (documented, not asserted — different fix path):
#   Equal-length MNV / block substitution (len(REF)==len(ALT)>1) co-occurring in
#   cis with an indel. The enricher (SNPsProcess) silently DROPS multi-base ALTs,
#   and _is_indel excludes len(ref)==len(alt), so the MNV alt bases are in NO
#   allele-set S_i and on NO fake contig -> membership `S_i superset {H_i}` is
#   violated at the MNV columns -> D over-counts -> verified D=2 vs real d=0.
#   Reachability: the shipped 1000G chr22 VCF is already normalized (0 MNV
#   records), so it does not trigger as-shipped; the risk is UN-NORMALIZED merges.
#   Fix: `bcftools norm -a` atomization (decompose MNVs to per-column SNPs) before
#   enrichment, so their alt bases enter the IUPAC sets.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
