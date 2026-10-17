"""Test the dict-less HIGH-VARIANT-DENSITY behavior in --per-sample mode (2.7.10).

Historically a capped dense window emitted ONE cross-sample greedy min-mismatch
"worst-case" representative so the region always surfaced >=1 off-target row. That rep
unioned carriers across DIFFERENT samples' best-per-column alts -- a chimera NO single
individual carries -- which over-stated the carrier roster and the risk (the real-world
chr14:106770868 15-SNP / 47-"carrier" / CFD-0.57 over-call). 2.7.10 PER-SAMPLE CONTRACT:
an off-target is reported ONLY if >=1 real individual carries it, so this synthetic
cross-sample rep is NO LONGER emitted in --per-sample mode. The dense region is still
surfaced via the high-variant-density BED (and the dense-region flag column) + the
reference off-target row -- we simply never fabricate a carrier-less worst case.

Fixture: 4 UNPHASED variant columns for one sample S, allowed_mms=1. The observed union
{all 4} = 2 mismatches FAILS the budget, the reference {none} = 2 mismatches FAILS, and
capped (IUPAC_CAP=2 < 4 cols) the sub-combination enumeration collapses to the union, so
NO observed haplotype survives. Under the fix the region therefore surfaces NO variant row
(only the .bed entry). Uncapped (IUPAC_CAP=10) the 0mm sub-combination IS enumerated and
surfaces with its real carrier -- the control, unchanged. (A single-sample in-budget
sub-combination lost to the UNPHASED cap is a pre-existing, documented limitation, rare on
the phased default panel where the CONFIRMED per-slot path is uncapped.)

Run with:
    cd PostProcess && python3 -m unittest test_hvdr_greedy_cap -v
"""

import io
import os
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN


class _GtSentinel(object):
    def __len__(self):
        return 1


def _build(iupac_cap):
    # reference: mismatches at cols 0,1 (C vs the all-'A' guide); matches at 2,3; NGG PAM.
    ref = ["A"] * TARGET_LEN
    ref[0] = "C"
    ref[1] = "C"
    ref[21], ref[22] = "G", "G"  # PAM cols 20,21,22 == A,G,G -> matches NGG
    genome = "".join(ref)

    dna = list(genome)
    for i in range(4):
        dna[i] = "N"  # ambiguity marker at the 4 variant columns
    dna = "".join(dna)

    split = [
        "X", GUIDE, dna, "chrT", "0", "0", "+",
        "4", "0", "4", "NGG", "y", "NA",
        "NA", "NA", "NA", "NA", "NA",
    ]
    # UNPHASED (0/1): every non-empty subset is a candidate cis haplotype.
    mydict = {
        "chrT,1": "S:0/1;C,A;rs0;0.10",  # C->A FIXES a mismatch
        "chrT,2": "S:0/1;C,A;rs1;0.10",  # C->A FIXES a mismatch
        "chrT,3": "S:0/1;A,T;rs2;0.10",  # A->T BREAKS a match
        "chrT,4": "S:0/1;A,T;rs3;0.10",  # A->T BREAKS a match
    }
    overrides = dict(
        genomeStr=genome,
        current_chr="chrT",
        mydict=mydict,
        myreg=None,
        mygt=_GtSentinel(),          # selects the observed (dict-less) branch
        haplotype_check=False,
        IUPAC_CAP=iupac_cap,
        hvdr_bed=io.StringIO(),
        pam="NGG",
        pos_beg=0,
        pos_end=-3,
        pam_begin=-3,
        pam_end=None,
        allowed_mms=1,               # union (2mm) fails; greedy (0mm) passes
    )
    return split, overrides


def _run(iupac_cap):
    split, overrides = _build(iupac_cap)
    ns = _load_pure_functions(overrides)
    ns["_phase_confirmation_rows"].clear()
    ns["_phase_confirmation_keys"].clear()
    cluster = []
    ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
    return cluster, ns


def _is_reference_row(row):
    return row[12] == "n" and row[-3] == "n" and row[-2] == 55


def _variant_rows(rows):
    return [r for r in rows if not _is_reference_row(r)]


class TestHvdrGreedyCap(unittest.TestCase):
    def test_capped_region_emits_no_synthetic_chimera_row(self):
        """CAPPED (IUPAC_CAP=2 < 4 cols), every observed haplotype gated by the mismatch
        budget: --per-sample emits NO synthetic cross-sample worst-case row (2.7.10 "no
        observed carrier -> no report"). The region is still surfaced by the reference row
        + the high-variant-density BED, never a fabricated carrier-less chimera."""
        cluster, ns = _run(iupac_cap=2)
        var = _variant_rows(cluster)
        self.assertEqual(
            var, [],
            "capped dense region must NOT fabricate a carrier-less worst-case row; got %r"
            % [(r[7], r[12]) for r in var],
        )
        # the region is STILL surfaced in the high-variant-density BED (visibility kept)
        self.assertIn("chrT", ns["hvdr_bed"].getvalue())

    def test_uncapped_also_finds_it_via_subcombinations(self):
        """CONTROL: uncapped (IUPAC_CAP=10) the sub-combination enumeration already
        finds the same 0mm subset, so the region is never missed there either."""
        cluster, ns = _run(iupac_cap=10)
        var = _variant_rows(cluster)
        self.assertTrue(any(int(r[7]) == 0 for r in var),
                        "uncapped enumeration must also surface the 0mm subset")


if __name__ == "__main__":
    unittest.main()
