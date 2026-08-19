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
    return row[17].count(",") + 1 if row[17] not in ("", "NA") else 0


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

    def test_unphased_single_putative_union(self):
        # UNPHASED N=4: ONE PUTATIVE union haplotype (all 4), sample present once.
        rows, ns = _run(4, phased=False)
        # Only the full union set exists (one sample, one putative set).
        self.assertEqual(len(rows), 1,
                         "unphased must emit ONE union row, got %d" % len(rows))
        self.assertEqual(_n_snps_in_row(rows[0]), 4)
        self.assertEqual(rows[0][12], "S")


class ScoringSentinelsIntact(unittest.TestCase):
    """The new path must NOT append the phase flag to final_line: the CFD/CRISTA
    scorers read the Reference / ref-score / tmp_pos tail by NEGATIVE index."""

    def test_negative_index_tail_shape(self):
        rows, ns = _run(3, phased=True)
        self.assertTrue(rows)
        row = rows[0]
        # tail appended by the finalizer: [..., refSeq_with_bulges, 33, tmp_pos_mms]
        self.assertEqual(row[-2], 33, "ref-score sentinel must be at target[-2] "
                         "pre-pop (becomes target[-3] after the CFD append)")
        self.assertIsInstance(row[-1], int, "tmp_pos_mms (int) must be the last token")
        self.assertIsInstance(row[-3], str, "refSeq_with_bulges must be at target[-3]")
        # refSeq_with_bulges is the all-ref window with NGG PAM, length 23.
        self.assertEqual(len(row[-3]), TARGET_LEN)

    def test_no_extra_column_on_final_line(self):
        # final_line must have exactly the legacy 18 split columns + 3 appended
        # (refSeq, 33, tmp_pos) = 21 tokens; NO phase-confirmation column spliced in.
        rows, ns = _run(1, phased=True)
        self.assertTrue(rows)
        self.assertEqual(len(rows[0]), 18 + 3)


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
        # Each single-SNP row carries exactly its one real carrier.
        for r in cluster:
            self.assertIn(r[12], ("A", "B"))


if __name__ == "__main__":
    unittest.main()
