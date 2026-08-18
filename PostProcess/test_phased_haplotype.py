"""Regression test for the phased multi-SNP haplotype under-report.

BUG (fix/phased-haplotype-underreport)
--------------------------------------
``new_simple_analysis.iupac_decomposition`` decomposes an IUPAC (variant)
off-target into every concrete alt-allele haplotype. ``totalDict[count][level]``
holds the combination lattice: ``count`` is the phase (hap0 / hap1), ``level`` is
the number-of-SNPs-minus-one, and level ``L+1`` grows by crossing each level-``L``
combo against the LEVEL-0 single-SNP seeds only.

For PHASED data the old code, the moment it created a level-(L+1) combo,
SUBTRACTED that combo's sample set from both its parent AND the level-0 seed it
had just consumed (a dedup device: report a cis haplotype once, on its maximal
combo). But because deeper levels grow ONLY against those very level-0 seeds,
emptying them mid-loop STARVED levels >= 2: a sample carrying k>=4 cis alts got
"used up" into disjoint 2-SNP pairs, no level-0 seed still contained it, and its
maximal k-SNP haplotype could never form. Finalization then dropped the surviving
lossy sub-combos whose residual mismatches exceeded ``allowed_mms`` -> the true
0-mismatch haplotype was silently under-reported (N=4 -> min 2mm, N=6 -> 4mm,
N=8 -> 6mm, N=12 full-enum -> ZERO rows). The UNPHASED path (no subtraction) and
the >cap GREEDY path (single level-0 entry, no combination) were always correct,
so the pipeline was self-inconsistent across the IUPAC_CAP boundary.

FIX
---
Deferred subtraction. Phase A: the growth loop runs mutation-free, so the full
2^N-1 lattice forms (identical to the already-correct unphased enumeration).
Phase B (phased, non-capped only): AFTER the loop, peel each sample off every
shorter (strict (pos,elem)-subset) ancestor combo, so it survives only on its
maximal cis combo -- reproducing the dedup the old subtraction intended, now that
the maximal combo actually exists.

WHAT THIS TEST ASSERTS
----------------------
Driving ``iupac_decomposition`` directly with synthetic fixtures:
  * PHASED cis N in {3,4,6,8,12}: the fully-substituted 0-mismatch haplotype IS
    emitted, carrying ALL N SNPs and the correct sample set (the pre-fix bug is
    gone). N=12 is exercised in full-enumeration (cap raised) -- the pre-fix
    zero-rows case.
  * UNPHASED N>=4: unchanged -- the full lattice (2^N-1 rows) including the 0mm
    maximal combo, every sample on every sub-combo (cis unverifiable).
  * PHASED TRANS (alts on alternating haplotypes of one sample): variants on
    different haplotypes do NOT fully combine (no false cis attribution).
  * >cap GREEDY path: unaffected -- exactly one 0mm representative per haplotype.

The harness AST-loads only the pure functions of ``new_simple_analysis.py``
(everything before the ``inFasta = open(...)`` module body) and injects the
runtime globals the function reads, so no genome / dict / arg files are needed.

Run with:
    cd PostProcess && python3 -m unittest test_phased_haplotype -v
"""

import ast
import io
import os
import sys
import types
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
SRC_PATH = os.path.join(PP, "new_simple_analysis.py")


def _load_pure_functions(global_overrides):
    """Exec only the import / constant / function-definition prologue of
    new_simple_analysis.py (everything BEFORE the side-effecting module body that
    opens ``sys.argv`` files, marked by ``inFasta = open(...)``), then inject the
    runtime globals ``iupac_decomposition`` reads. Returns the module namespace."""
    with open(SRC_PATH) as fh:
        tree = ast.parse(fh.read())
    keep = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            getattr(t, "id", None) == "inFasta" for t in node.targets
        ):
            break  # first side-effecting statement of the module body -> stop
        keep.append(node)
    module = ast.Module(body=keep, type_ignores=[])
    # The kept prologue imports numpy / pandas / CRISTA_score at module top, but
    # iupac_decomposition uses NONE of them at runtime. The light unit-tests CI env
    # has no numpy/pandas/scientific stack (the other PostProcess tests are
    # stdlib-only by design), so exec would die on `import numpy`. Inject harmless
    # stubs ONLY when the real package is unavailable — a full local env keeps using
    # the real ones. (This is exactly why the suite passed locally but failed in CI.)
    for _name, _attrs in (("numpy", ()), ("pandas", ()),
                          ("CRISTA_score", ("CRISTA_predict_list",))):
        try:
            __import__(_name)
        except Exception:
            _stub = types.ModuleType(_name)
            for _a in _attrs:
                setattr(_stub, _a, lambda *a, **k: None)
            sys.modules[_name] = _stub
    namespace = {}
    exec(compile(module, SRC_PATH, "exec"), namespace)
    namespace.update(global_overrides)
    return namespace


# Standard SpCas9 geometry: 20-nt protospacer + 3-nt NGG PAM at the 3' end.
GUIDE_NO_PAM = "A" * 20
GUIDE = GUIDE_NO_PAM + "NGG"
TARGET_LEN = 23


def _build_fixture(n_snps, phased=True, trans=False, cap=10, allowed_mms=6):
    """Construct a synthetic single-target fixture with ``n_snps`` co-occurring
    ambiguity positions in the protospacer.

    Reference genome window = 'C' at each SNP position (a mismatch vs the all-'A'
    guide) and 'A' elsewhere, with an ``NGG`` PAM. Every SNP's ALT allele is 'A'
    (matches the guide), so:
      * the reference target has ``n_snps`` mismatches, and
      * the FULLY-substituted haplotype has ZERO mismatches.
    One sample ``S`` carries every alt. Phasing:
      * phased cis  -> "S:1|0" at every SNP (all alts on hap0),
      * phased trans -> alternate "1|0"/"0|1" (alts split across haplotypes),
      * unphased    -> "S:0/1".
    """
    ref = ["A"] * TARGET_LEN
    for i in range(n_snps):
        ref[i] = "C"
    ref[20], ref[21], ref[22] = "A", "G", "G"  # N G G PAM
    genome = "".join(ref)

    dna = list(genome)
    for i in range(n_snps):
        dna[i] = "W"  # ambiguity code -> triggers the variant decomposition branch
    dna = "".join(dna)

    # split columns: 0 Bulge_type 1 crRNA 2 DNA 3 Chromosome 4 Position
    # 5 Cluster_Position 6 Direction 7 Mismatches 8 Bulge_Size 9 Total 10 PAM_gen
    # 11 Var_uniq 12 Samples 13.. (annotation/rsID/AF/SNP placeholders)
    split = [
        "X", GUIDE, dna, "chrT", "0", "0", "+",
        str(n_snps), "0", str(n_snps), "NGG", "y", "NA",
        "NA", "NA", "NA", "NA", "NA",
    ]

    mydict = {}
    for i in range(n_snps):
        if not phased:
            gt = "S:0/1"
        elif trans:
            gt = "S:1|0" if i % 2 == 0 else "S:0|1"
        else:
            gt = "S:1|0"
        # entry format parsed by retrieveFromDict: "samples;ref,alt;rsID;AF"
        mydict["chrT," + str(i + 1)] = gt + ";C,A;rs%d;0.10" % i

    overrides = dict(
        genomeStr=genome,
        current_chr="chrT",
        mydict=mydict,
        haplotype_check=phased,
        IUPAC_CAP=cap,
        hvdr_bed=io.StringIO(),  # capped path writes a BED line; swallow it
        pam="NGG",
        pos_beg=0,
        pos_end=-3,
        pam_begin=-3,
        pam_end=None,
        allowed_mms=allowed_mms,
    )
    return split, overrides


def _decompose(n_snps, phased=True, trans=False, cap=10, allowed_mms=6):
    """Run iupac_decomposition on a fixture; return the list of emitted rows.

    Each row is the ``final_line`` list. Useful columns: [7]=Mismatches (relative
    to ref), [12]=Samples, [17]=SNP info (comma-joined per SNP)."""
    split, overrides = _build_fixture(n_snps, phased, trans, cap, allowed_mms)
    ns = _load_pure_functions(overrides)
    cluster = []
    ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
    return cluster


def _n_snps_in_row(row):
    """Number of SNPs folded into a row (SNP-info column 17 is comma-joined)."""
    return row[17].count(",") + 1 if row[17] not in ("", "NA") else 0


class PhasedFullHaplotypeFound(unittest.TestCase):
    """The fully-substituted 0-mismatch phased haplotype must be emitted for every
    N of co-occurring cis SNPs -- this is the bug the fix targets."""

    def _assert_full_haplotype(self, n_snps, cap=10, allowed_mms=6):
        rows = _decompose(n_snps, phased=True, cap=cap, allowed_mms=allowed_mms)
        self.assertTrue(rows, f"N={n_snps}: no rows emitted (haplotype dropped)")
        zero_mm = [r for r in rows if int(r[7]) == 0]
        self.assertEqual(
            len(zero_mm), 1,
            f"N={n_snps}: expected exactly one 0-mismatch full haplotype, "
            f"got {[(r[7], r[12], _n_snps_in_row(r)) for r in rows]}",
        )
        row = zero_mm[0]
        self.assertEqual(
            _n_snps_in_row(row), n_snps,
            f"N={n_snps}: the 0mm target must fold in all {n_snps} SNPs, "
            f"got {_n_snps_in_row(row)}",
        )
        self.assertEqual(
            row[12], "S",
            f"N={n_snps}: the 0mm haplotype must be attributed to sample S, "
            f"got {row[12]!r}",
        )

    def test_N3(self):
        # N=3 was already correct pre-fix; must stay correct.
        self._assert_full_haplotype(3)

    def test_N4(self):
        # Primary bug: pre-fix gave min 2mm (true 0mm dropped).
        self._assert_full_haplotype(4)

    def test_N6(self):
        # Pre-fix -> min 4mm.
        self._assert_full_haplotype(6)

    def test_N8(self):
        # Pre-fix -> min 6mm.
        self._assert_full_haplotype(8)

    def test_N12_full_enumeration(self):
        # Pre-fix, in full-enumeration (cap raised above N so it does NOT take the
        # correct greedy path), every surviving sub-combo carried >6 residual
        # mismatches -> ALL dropped -> ZERO rows. The fix recovers the single 0mm
        # 12-SNP haplotype.
        self._assert_full_haplotype(12, cap=20, allowed_mms=6)


class UnphasedUnchanged(unittest.TestCase):
    """Unphased data cannot prove cis, so the code keeps every sample on every
    single/sub-combo AND the maximal combo (the conservative superset). The fix
    (Phase B is phased-only) must leave this behavior untouched."""

    def _assert_full_lattice(self, n_snps):
        rows = _decompose(n_snps, phased=False)
        # full lattice for a single sample = all non-empty subsets = 2^N - 1 rows
        self.assertEqual(
            len(rows), (2 ** n_snps) - 1,
            f"unphased N={n_snps}: expected the full 2^N-1 lattice, got {len(rows)}",
        )
        zero_mm = [r for r in rows if int(r[7]) == 0]
        self.assertEqual(
            len(zero_mm), 1, f"unphased N={n_snps}: one 0mm maximal combo expected"
        )
        self.assertEqual(_n_snps_in_row(zero_mm[0]), n_snps)
        self.assertEqual(zero_mm[0][12], "S")

    def test_unphased_N4(self):
        self._assert_full_lattice(4)

    def test_unphased_N6(self):
        self._assert_full_lattice(6)


class PhasedTransDoesNotCombine(unittest.TestCase):
    """Cis-semantics guard: alts on DIFFERENT haplotypes of one sample must NOT
    fully combine into a single 0mm target (no trans/false attribution)."""

    def test_trans_N4_no_full_combo(self):
        rows = _decompose(4, phased=True, trans=True)
        self.assertTrue(rows)
        # No row may fold all 4 alts together, because they are not all in cis.
        self.assertFalse(
            any(_n_snps_in_row(r) == 4 for r in rows),
            "trans alts wrongly combined into a full 4-SNP haplotype",
        )
        # And there must be NO 0-mismatch full target for sample S.
        self.assertFalse(
            any(int(r[7]) == 0 and _n_snps_in_row(r) == 4 for r in rows),
            "trans variants produced a false 0mm cis haplotype",
        )


class GreedyCapPathUnaffected(unittest.TestCase):
    """The >IUPAC_CAP greedy high-variant-density path (single level-0
    representative per haplotype, no combination, no peel) must be unaffected: it
    already reports one 0mm target and the fix's Phase B is guarded by
    ``not capped``."""

    def test_greedy_over_cap(self):
        # N=4 with cap=2 -> countIUPAC (4) > cap (2) -> greedy branch.
        rows = _decompose(4, phased=True, cap=2)
        self.assertEqual(len(rows), 1, "greedy path must emit exactly one row")
        row = rows[0]
        self.assertEqual(int(row[7]), 0, "greedy representative must be the 0mm haplotype")
        self.assertEqual(_n_snps_in_row(row), 4)
        self.assertEqual(row[12], "S")


if __name__ == "__main__":
    unittest.main()
