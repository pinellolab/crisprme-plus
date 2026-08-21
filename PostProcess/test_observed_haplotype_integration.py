"""Integration tests for the DICT-LESS observed-haplotype branch of
``new_simple_analysis.iupac_decomposition``.

These drive the REAL ``iupac_decomposition`` through the same AST-load harness the
legacy ``test_phased_haplotype`` uses (exec only the pre-``inFasta`` prologue, inject
runtime globals), but with ``mygt`` PRESENT (a truthy sentinel) so the new observed-
haplotype branch is taken. Carrier genotypes are supplied via the legacy ``mydict``
entries (``retrieveFromDict`` returns the "sample:gt" tokens the enumerator splits),
with ``myreg`` left None -- the observed branch only needs the 5-tuple's carrier
tokens, which the dict decode provides.

We assert the finalized rows carry the right Samples column, mismatch count, SNP
folding, and (critically) that the Reference / ref-score sentinel tail is intact at the
expected NEGATIVE indices (target[-3]==33, target[-2]==tmp_pos_mms, target[-4]==
refSeq) so CFD/CRISTA scoring is NOT corrupted by the new path. We also assert the
phase-confirmation accumulator captures CONFIRMED vs PUTATIVE correctly.

STDLIB ONLY (+ the AST harness's numpy/pandas/CRISTA stubs). No genome / dict / arg
files.

Run with:
    cd PostProcess && python3 -m unittest test_observed_haplotype_integration -v
"""

import io
import os
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN


class _GtSentinel(object):
    """Truthy stand-in for a Tier-1 GenotypeReader. The observed branch only checks
    ``mygt is not None`` (it reads carriers via retrieveFromDict from the dict here),
    so a bare truthy object is enough to select the branch."""

    def __len__(self):
        return 1


def _build(n_snps, phased=True, trans=False, allowed_mms=6, multiallelic=False):
    """Build a fixture + overrides that SELECT the observed branch (mygt present).

    Reference window: 'C' at each SNP column (mismatch vs all-'A' guide), 'A'
    elsewhere, NGG PAM. Every SNP's ALT is 'A' (matches guide) so the fully-substituted
    haplotype is 0-mismatch. One sample S carries every alt.
    """
    ref = ["A"] * TARGET_LEN
    for i in range(n_snps):
        ref[i] = "C"
    ref[20], ref[21], ref[22] = "A", "G", "G"
    genome = "".join(ref)

    dna = list(genome)
    for i in range(n_snps):
        dna[i] = "W"  # ambiguity code -> variant decomposition branch
    dna = "".join(dna)

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
        mydict["chrT," + str(i + 1)] = gt + ";C,A;rs%d;0.10" % i

    overrides = dict(
        genomeStr=genome,
        current_chr="chrT",
        mydict=mydict,
        myreg=None,
        mygt=_GtSentinel(),          # SELECTS the observed branch
        haplotype_check=phased,      # unused on the observed branch
        IUPAC_CAP=10,
        hvdr_bed=io.StringIO(),
        pam="NGG",
        pos_beg=0,
        pos_end=-3,
        pam_begin=-3,
        pam_end=None,
        allowed_mms=allowed_mms,
    )
    return split, overrides


def _run(n_snps, **kw):
    split, overrides = _build(n_snps, **kw)
    ns = _load_pure_functions(overrides)
    # reset the module-level phase-confirmation accumulator between runs
    ns["_phase_confirmation_rows"].clear()
    ns["_phase_confirmation_keys"].clear()
    cluster = []
    ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
    return cluster, ns


def _n_snps_in_row(row):
    return row[17].count(",") + 1 if row[17] not in ("", "NA", "n") else 0


def _is_reference_row(row):
    """A candidate's REFERENCE off-target row (LOCUS-COVERAGE FIX): pure reference
    alignment, no carriers (Samples == legacy literal 'n'), and the legacy 'n' /
    sentinel-55 tail (Reference col == 'n', target[-2] == 55) that bins it origin=ref
    downstream. The row is byte-identical to the legacy non-IUPAC else-branch reference
    row: its Samples/rsID/AF/SNP markers are the literal 'n' (NOT 'NA'); remove_n_and_
    dots rewrites those to the literal 'NA' string that resultIntegrator's target[13]/
    target[18] guards depend on."""
    return row[12] == "n" and row[-3] == "n" and row[-2] == 55


def _variant_rows(rows):
    """Only the observed VARIANT haplotype rows (drop the reference off-target row)."""
    return [r for r in rows if not _is_reference_row(r)]


def _reference_rows(rows):
    return [r for r in rows if _is_reference_row(r)]


class ObservedBranchSelected(unittest.TestCase):
    def test_phased_cis_full_haplotype_emitted(self):
        # PHASED cis N=4: ONE 0-mismatch haplotype folding all 4 SNPs, sample S.
        rows, ns = _run(4, phased=True)
        zero = [r for r in rows if int(r[7]) == 0]
        self.assertEqual(len(zero), 1,
                         "expected exactly one 0mm full haplotype, got %r"
                         % [(r[7], r[12], _n_snps_in_row(r)) for r in rows])
        row = zero[0]
        self.assertEqual(_n_snps_in_row(row), 4)
        self.assertEqual(row[12], "S")

    def test_phased_trans_no_full_combo(self):
        # PHASED trans N=4: alts alternate slots -> NO single 4-SNP cis haplotype.
        rows, ns = _run(4, phased=True, trans=True)
        self.assertTrue(rows)
        self.assertFalse(any(_n_snps_in_row(r) == 4 for r in rows),
                         "trans alts wrongly combined into a full cis haplotype")

    def test_unphased_enumerates_subcombinations(self):
        # UNPHASED N=4: the phasing is unknown, so EVERY non-empty subset of the sample's
        # 4 variants is a candidate cis haplotype = 2^4 - 1 = 15 PUTATIVE variant rows
        # (the finalize mm/PAM gates prune any out-of-budget / PAM-invalid subset; here
        # all pass), sample present in each -- PLUS the candidate's single REFERENCE
        # off-target row (LOCUS-COVERAGE FIX). This is the fix for the union-only miss:
        # a valid sub-combination off-target (needing a reference allele at a het column,
        # e.g. PAM creation) is no longer hidden by the maximal union.
        rows, ns = _run(4, phased=False)
        var = _variant_rows(rows)
        self.assertEqual(len(var), 15,
                         "unphased must emit all 2^4-1 subset rows, got %d" % len(var))
        # the maximal union (all 4) is among them, and every row carries sample S.
        self.assertTrue(any(_n_snps_in_row(r) == 4 for r in var))
        for r in var:
            self.assertEqual(r[12], "S")
        # exactly ONE reference off-target row for the candidate.
        self.assertEqual(len(_reference_rows(rows)), 1)


class ScoringSentinelsIntact(unittest.TestCase):
    """The new path must NOT append the phase flag to final_line: the CFD/CRISTA
    scorers read the Reference / ref-score / tmp_pos tail by NEGATIVE index."""

    def test_negative_index_tail_shape(self):
        rows, ns = _run(3, phased=True)
        # Assert the tail shape on a VARIANT row (the finalizer's alt-row tail):
        # [..., refSeq_with_bulges, 33, tmp_pos_mms].
        var = _variant_rows(rows)
        self.assertTrue(var)
        row = var[0]
        self.assertEqual(row[-2], 33, "ref-score sentinel must be at target[-2] "
                         "pre-pop (becomes target[-3] after the CFD append)")
        self.assertIsInstance(row[-1], int, "tmp_pos_mms (int) must be the last token")
        self.assertIsInstance(row[-3], str, "refSeq_with_bulges must be at target[-3]")
        # refSeq_with_bulges is the all-ref window with NGG PAM, length 23.
        self.assertEqual(len(row[-3]), TARGET_LEN)

    def test_reference_row_tail_shape(self):
        # The candidate's REFERENCE off-target row carries the legacy non-IUPAC
        # else-branch tail: [..., "n", 55, tmp_pos_mms]. The "n" + 55 sentinel bin it
        # origin=ref downstream (no separate ref-CFD recompute); the DNA is REF.
        rows, ns = _run(3, phased=True)
        ref = _reference_rows(rows)
        self.assertEqual(len(ref), 1, "exactly one reference off-target expected")
        row = ref[0]
        self.assertEqual(row[-3], "n", "Reference column must be 'n' for a ref row")
        self.assertEqual(row[-2], 55, "ref sentinel 55 (DNA is REF) must be at [-2]")
        self.assertIsInstance(row[-1], int, "tmp_pos_mms (int) must be the last token")
        self.assertEqual(row[12], "n",
                         "reference row Samples must be the legacy literal 'n' "
                         "(round-trips to 'NA' via remove_n_and_dots -> origin=ref)")
        self.assertEqual(row[17], "n",
                         "reference row SNP/snp_info must be the legacy literal 'n' "
                         "(the load-bearing col18 marker downstream)")

    def test_no_extra_column_on_final_line(self):
        # BOTH the variant AND the reference row must have exactly the legacy 18 split
        # columns + 3 appended (refSeq/"n", 33/55, tmp_pos) = 21 tokens; NO
        # phase-confirmation column spliced in.
        rows, ns = _run(1, phased=True)
        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(len(row), 18 + 3)


class PhaseConfirmationAccumulator(unittest.TestCase):
    def test_confirmed_flag_for_phased_cis(self):
        rows, ns = _run(3, phased=True)
        recs = ns["_phase_confirmation_rows"]
        self.assertTrue(recs)
        self.assertTrue(all(r["Phase_Confirmation"] == "CONFIRMED" for r in recs),
                        "phased cis must be CONFIRMED, got %r"
                        % [r["Phase_Confirmation"] for r in recs])

    def test_putative_flag_for_unphased(self):
        rows, ns = _run(3, phased=False)
        recs = ns["_phase_confirmation_rows"]
        self.assertTrue(recs)
        self.assertTrue(all(r["Phase_Confirmation"] == "PUTATIVE" for r in recs),
                        "unphased must be PUTATIVE, got %r"
                        % [r["Phase_Confirmation"] for r in recs])

    def test_companion_writer_round_trip(self):
        import phase_confirmation_companion as pcc
        import tempfile
        rows, ns = _run(2, phased=True)
        recs = ns["_phase_confirmation_rows"]
        with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tf:
            path = tf.name
        try:
            n = pcc.write_companion(path, recs)
            self.assertEqual(n, len(recs))
            with open(path) as fh:
                body = fh.read()
            self.assertIn("Phase_Confirmation", body)
            self.assertIn("CONFIRMED", body)
        finally:
            os.remove(path)


class ChimeraExcludedEndToEnd(unittest.TestCase):
    def test_two_samples_no_chimera_row(self):
        # A carries SNP col0, B carries SNP col1; no sample carries both -> no 2-SNP
        # row. Build a custom fixture.
        n_snps = 2
        ref = ["A"] * TARGET_LEN
        ref[0] = ref[1] = "C"
        ref[20], ref[21], ref[22] = "A", "G", "G"
        genome = "".join(ref)
        dna = list(genome)
        dna[0] = dna[1] = "W"
        dna = "".join(dna)
        split = ["X", GUIDE, dna, "chrT", "0", "0", "+", "2", "0", "2", "NGG",
                 "y", "NA", "NA", "NA", "NA", "NA", "NA"]
        mydict = {
            "chrT,1": "A:1|0;C,A;rs0;0.10",   # only A carries col0
            "chrT,2": "B:1|0;C,A;rs1;0.10",   # only B carries col1
        }
        overrides = dict(
            genomeStr=genome, current_chr="chrT", mydict=mydict, myreg=None,
            mygt=_GtSentinel(), haplotype_check=True, IUPAC_CAP=10,
            hvdr_bed=io.StringIO(), pam="NGG", pos_beg=0, pos_end=-3,
            pam_begin=-3, pam_end=None, allowed_mms=6,
        )
        ns = _load_pure_functions(overrides)
        ns["_phase_confirmation_rows"].clear()
        ns["_phase_confirmation_keys"].clear()
        cluster = []
        ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
        self.assertFalse(any(_n_snps_in_row(r) == 2 for r in cluster),
                         "cross-individual chimera 2-SNP row must NOT be emitted")
        # Each single-SNP VARIANT row carries exactly its one real carrier (the
        # reference off-target row is excluded -- it has no carriers).
        for r in _variant_rows(cluster):
            self.assertIn(r[12], ("A", "B"))
        # exactly one reference off-target row (LOCUS-COVERAGE FIX), no carriers.
        self.assertEqual(len(_reference_rows(cluster)), 1)


class ReferenceOffTargetCoverage(unittest.TestCase):
    """LOCUS-COVERAGE FIX regression: the observed path must emit each IUPAC
    candidate's REFERENCE off-target (matching the legacy embedded-refSeq/33
    representation), so its locus is NEVER dropped -- even when NO observed variant
    haplotype survives -- while still excluding phantom multi-variant combinations and
    keeping exact per-individual carriers.

    Pre-fix, ``_iupac_decomposition_observed`` early-returned on an empty
    ``haplotypes`` and emitted the reference ONLY inside surviving alt rows, so a
    candidate with no productive observed variant lost its entire locus (the ~58k
    dropped loci found by the ml007 genome-wide e2e).
    """

    def _run_gts(self, gt_by_col, n_snps=None, allowed_mms=6):
        """Drive the observed path with an explicit per-column genotype token so we can
        model a no-productive-variant candidate (all-ref / missing genotypes)."""
        if n_snps is None:
            n_snps = len(gt_by_col)
        ref = ["A"] * TARGET_LEN
        for i in range(n_snps):
            ref[i] = "C"
        ref[20], ref[21], ref[22] = "A", "G", "G"
        genome = "".join(ref)
        dna = list(genome)
        for i in range(n_snps):
            dna[i] = "W"
        dna = "".join(dna)
        split = ["X", GUIDE, dna, "chrT", "0", "0", "+",
                 str(n_snps), "0", str(n_snps), "NGG", "y", "NA",
                 "NA", "NA", "NA", "NA", "NA"]
        mydict = {}
        for i in range(n_snps):
            mydict["chrT," + str(i + 1)] = gt_by_col[i] + ";C,A;rs%d;0.10" % i
        overrides = dict(
            genomeStr=genome, current_chr="chrT", mydict=mydict, myreg=None,
            mygt=_GtSentinel(), haplotype_check=True, IUPAC_CAP=10,
            hvdr_bed=io.StringIO(), pam="NGG", pos_beg=0, pos_end=-3,
            pam_begin=-3, pam_end=None, allowed_mms=allowed_mms,
        )
        ns = _load_pure_functions(overrides)
        ns["_phase_confirmation_rows"].clear()
        ns["_phase_confirmation_keys"].clear()
        cluster = []
        ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
        return cluster, ns

    def test_no_productive_variant_still_emits_reference(self):
        # (a) Every sample is ref/ref at the ambiguity columns -> the observed
        # enumeration is EMPTY -> pre-fix this candidate emitted NOTHING and its locus
        # was dropped. Now it emits EXACTLY ONE reference off-target row.
        cluster, ns = self._run_gts(["S:0|0", "S:0|0", "S:0|0"])
        self.assertEqual(len(cluster), 1,
                         "no-productive-variant candidate must still emit its "
                         "reference off-target (got %d rows)" % len(cluster))
        row = cluster[0]
        self.assertTrue(_is_reference_row(row),
                        "the single emitted row must be the reference off-target")
        # Sentinel tail intact + reference sentinels.
        self.assertEqual(row[-3], "n")
        self.assertEqual(row[-2], 55)
        self.assertEqual(row[12], "n")      # no carriers (legacy 'n' marker)
        self.assertEqual(row[17], "n")      # snp_info sentinel (legacy 'n' marker)
        # exactly the legacy 18 + 3 tail columns (no spliced-in phase column).
        self.assertEqual(len(row), 18 + 3)

    def test_reference_alongside_variant_no_duplication(self):
        # (b) A productive phased-cis N=4 candidate now emits the 1 full-haplotype
        # variant row PLUS exactly 1 reference row (locus reference coverage restored),
        # with no duplicate reference row and no cross-individual chimera.
        cluster, ns = _run(4, phased=True)
        refs = _reference_rows(cluster)
        var = _variant_rows(cluster)
        self.assertEqual(len(refs), 1, "exactly one reference off-target row expected")
        # the productive 0mm full haplotype survives among the variant rows.
        zero = [r for r in var if int(r[7]) == 0]
        self.assertEqual(len(zero), 1)
        self.assertEqual(_n_snps_in_row(zero[0]), 4)
        self.assertEqual(zero[0][12], "S")
        # no cross-individual 4-SNP chimera beyond the single real cis haplotype and
        # the reference row.
        self.assertEqual(len(cluster), len(var) + 1)

    def test_reference_row_creates_no_phase_confirmation(self):
        # (c) The reference row must NOT create a phase-confirmation record: those are
        # for VARIANT off-targets only. The accumulator length == number of VARIANT
        # rows (deduped), never inflated by the reference row.
        cluster, ns = _run(3, phased=True)
        recs = ns["_phase_confirmation_rows"]
        var = _variant_rows(cluster)
        self.assertTrue(recs)
        # every recorded phase-confirmation identity must map to a VARIANT row's SNP
        # column (never the reference "NA").
        self.assertTrue(all(r["SNP"] != "NA" for r in recs),
                        "reference row (SNP == NA) must not appear in the "
                        "phase-confirmation accumulator")
        var_snp_keys = {(r[3], r[4], r[6], r[1], r[2], r[17]) for r in var}
        rec_keys = {(r["Chromosome"], r["Position"], r["Direction"], r["crRNA"],
                     r["DNA"], r["SNP"]) for r in recs}
        self.assertTrue(rec_keys <= var_snp_keys,
                        "phase-confirmation records must correspond only to variant rows")

    def test_over_budget_reference_dropped(self):
        # The reference row obeys the SAME budget gate as the finalizer
        # (mm_new_t - bulges > allowed_mms -> dropped), matching stable. A candidate
        # whose reference window has 5 mismatches is dropped at allowed_mms=4.
        cluster, ns = self._run_gts(["S:0|0"] * 5, allowed_mms=4)
        self.assertEqual(len(_reference_rows(cluster)), 0,
                         "an over-budget reference (5 mm > allowed 4) must be dropped, "
                         "matching the legacy finalizer's mm-threshold gate")
        # in-budget (4 mm at allowed 4) is KEPT.
        cluster, ns = self._run_gts(["S:0|0"] * 4, allowed_mms=4)
        self.assertEqual(len(_reference_rows(cluster)), 1,
                         "an in-budget reference (4 mm <= allowed 4) must be kept")


if __name__ == "__main__":
    unittest.main()
