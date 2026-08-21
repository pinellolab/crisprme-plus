"""Test the dict-less HIGH-VARIANT-DENSITY greedy-cap fix: in a capped dense window
the observed per-individual union may exceed the mismatch budget, so historically the
region appeared ONLY in the .high_variant_density_regions.bed with NO off-target row in
the results. The fix emits ONE greedy min-mismatch representative per capped region
(parity with the legacy dict cap), so every dense region ALWAYS surfaces >=1 off-target
row -- we never miss a region.

Fixture: 4 unphased variant columns for one sample S. Two alts FIX a reference mismatch
(C->A, matching the all-A guide), two alts BREAK a match (A->T). With allowed_mms=1:
  - the observed UNION {all 4} = 2 mismatches -> FAILS the budget (dropped),
  - the reference {none}          = 2 mismatches -> FAILS,
  - the greedy best {fix the two mismatches, skip the two breaks} = 0 mismatches -> PASSES.
Capped (IUPAC_CAP=2 < 4 columns) the union path emits nothing that survives, so the ONLY
surviving variant row is the greedy representative the fix adds.

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
    def test_capped_region_surfaces_greedy_best(self):
        """CAPPED (IUPAC_CAP=2 < 4 cols): the union (2mm) is dropped, but the greedy
        representative (0mm) is emitted so the dense region still surfaces >=1 row."""
        cluster, ns = _run(iupac_cap=2)
        var = _variant_rows(cluster)
        self.assertTrue(var, "capped dense region must surface >=1 off-target row")
        zero = [r for r in var if int(r[7]) == 0]
        self.assertEqual(
            len(zero), 1,
            "expected exactly one 0-mismatch greedy representative, got %r"
            % [(r[7], r[12]) for r in var],
        )
        g = zero[0]
        self.assertEqual(g[12], "S", "greedy carriers = union of the chosen alts' samples")
        # the .bed logged the region too
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
