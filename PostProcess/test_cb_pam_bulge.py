#!/usr/bin/env python3
"""A DNA bulge inside the PAM must survive the guide-PAM normalization.

The CRISPR-Bulge ensemble is handed ``(aligned guide, aligned off-target)``, with the
guide's own PAM normalized to NGG first. For years that normalization was a fixed
slice, ``crrna[:-pam_len] + "NGG"``. An in-place DNA-bulge alignment is one column
LONGER than the guide, so when the bulge's gap column falls inside the PAM the guide's
three PAM bases occupy four columns (``...NG-G``) and the slice silently ate the gap:

  * the ``total_gaps > 1`` null rule, which runs on the masked string, then counted 0
    gaps on a row that has one, and
  * the model received ``...NNGG`` -- a gapless 24-mer carrying a phantom base, which
    is not an alignment of this site, and scores like a bulge-free off-target.

Measured on the external CRISPR-Bulge benchmark (2,940,588 rows): 30,131 rows have the
gap in those columns, the fabricated score exceeds the true one on 28,704 of them, and
it lifts 167 rows over the sensitive 1.9e-4 default, 22 over 0.2 and 1 over 0.5 -- never
the other way. A CRISPRme *search* cannot produce such a row (it anchors the PAM
exactly), so this reaches results only through externally supplied lists
(``crisprme.py targets-scoring``), which is exactly where the benchmark came from.

These tests pin both halves: the gap survives, and on every row WITHOUT a gap in the
PAM columns the normalizer is byte-identical to the slice it replaced, so no result a
CRISPRme search ever produced can move.

STDLIB only. Hermetic: the scorer is stubbed, no cbulge env or TensorFlow needed.
"""
import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import cbulge_rows as cbr  # noqa: E402
import targets_scoring as ts  # noqa: E402

SPACER = "ACGTACGTACGTACGTACGT"          # 20 nt, synthetic
GUIDE_NGG = SPACER + "NGG"                # how an upstream aligned guide is written
GUIDE_NNN = SPACER + "NNN"                # how a CRISPRme NRR/NNN guide is written


def _old_mask(crrna, pam_len=3):
    """The slice this replaced, kept here as the oracle for the no-gap case."""
    return crrna[: len(crrna) - pam_len] + "NGG"


class TestMaskKeepsTheBulge(unittest.TestCase):
    def test_pam_internal_dna_bulge_survives(self):
        # guide PAM NGG spread over four columns by a 1-bp DNA bulge
        aligned = SPACER + "NG-G"
        self.assertEqual(cbr.mask_guide_pam(aligned), aligned)
        self.assertEqual(cbr.mask_guide_pam(aligned).count("-"), 1)
        # the bug: the slice deletes the gap and invents a base
        self.assertEqual(_old_mask(aligned), SPACER + "NNGG")
        self.assertEqual(_old_mask(aligned).count("-"), 0)

    def test_pam_internal_bulge_on_an_iupac_pam_is_still_normalized(self):
        # an NRR/NNN guide with the bulge inside its PAM: normalized AND gap kept
        self.assertEqual(cbr.mask_guide_pam(SPACER + "NN-N"), SPACER + "NG-G")
        self.assertEqual(cbr.mask_guide_pam(SPACER + "N-NN"), SPACER + "N-GG")
        self.assertEqual(cbr.mask_guide_pam(SPACER + "NNN-"), SPACER + "NGG-")

    def test_spacer_bulge_is_untouched_and_pam_normalized(self):
        # the ordinary case a search produces: gap in the spacer, PAM at the end
        aligned = SPACER[:10] + "-" + SPACER[10:] + "NNN"
        self.assertEqual(cbr.mask_guide_pam(aligned), SPACER[:10] + "-" + SPACER[10:] + "NGG")

    def test_no_gap_in_pam_is_byte_identical_to_the_slice_it_replaced(self):
        cases = [
            GUIDE_NGG, GUIDE_NNN, SPACER + "NRG", SPACER + "AGG",
            SPACER[:3] + "-" + SPACER[3:] + "NGG",      # 5' spacer bulge
            SPACER[:19] + "-" + SPACER[19:] + "NNN",    # PAM-proximal spacer bulge
            "ACGT" + "NGG",                              # short guide
        ]
        for c in cases:
            with self.subTest(c=c):
                self.assertEqual(cbr.mask_guide_pam(c), _old_mask(c))

    def test_other_pam_widths(self):
        # a mask whose width matches pam_len: gap-preserving, same rule
        self.assertEqual(cbr.mask_guide_pam("A" * 20 + "TT-TN", pam_len=4, pam="TTTN"),
                         "A" * 20 + "TT-TN")
        self.assertEqual(cbr.mask_guide_pam("A" * 20 + "NN", pam_len=2, pam="GG"),
                         "A" * 20 + "GG")
        # a pam_len that does NOT match the mask keeps the old slice verbatim, so the
        # caller's length rule still nulls a guide outside the 20+3 shape
        self.assertEqual(cbr.mask_guide_pam("A" * 20 + "TTTN", pam_len=4),
                         _old_mask("A" * 20 + "TTTN", 4))
        self.assertEqual(len(cbr.mask_guide_pam("A" * 20 + "TTTN", pam_len=4)), 23)


class TestTargetsScoringPairing(unittest.TestCase):
    """The CLI's _cb_pair must hand the model the alignment, and judge it as one."""

    def _pair(self, crrna, dna, pam_len=3):
        site = ts.Site(crrna=crrna, dna=dna, pam_len=pam_len, bulge_type="DNA",
                       bulge_size=1, mismatches=0)
        null = []
        sg, off = ts._cb_pair(site, 0, null)
        return sg, off, null

    def test_pam_internal_bulge_reaches_the_model_as_a_bulge(self):
        crrna = SPACER + "NG-G"
        dna = "ACGTACGTACGTACGTACGTAGAG"
        sg, off, null = self._pair(crrna, dna)
        self.assertEqual(null, [], "a single 1-bp bulge is in the model's domain")
        self.assertEqual(sg, crrna)
        self.assertEqual(off, dna)
        self.assertEqual(sg.count("-") + off.count("-"), 1)
        self.assertEqual(len(sg), len(off))

    def test_two_bulges_one_of_them_in_the_pam_is_still_nulled(self):
        # gap in the spacer AND a gap in the PAM = 2 gap columns -> out of domain.
        # The old slice hid the PAM one and let this row through with 1 gap.
        crrna = SPACER[:10] + "-" + SPACER[10:] + "NG-G"
        dna = "A" * 21 + "GGAG"
        self.assertEqual(len(crrna), len(dna))
        sg, off, null = self._pair(crrna, dna)
        self.assertEqual(null, [0])
        self.assertEqual(sg, ts._CB_SG_DUMMY)
        self.assertEqual(_old_mask(crrna).count("-") + dna.count("-"), 1,
                         "the slice saw only one gap here -- the defect this pins")

    def test_mismatch_only_row_is_unchanged(self):
        crrna = GUIDE_NNN
        dna = "ACGTACGTACGTACGTACCTAGG"
        sg, off, null = self._pair(crrna, dna)
        self.assertEqual(null, [])
        self.assertEqual(sg, SPACER + "NGG")
        self.assertEqual(off, dna)


class TestPipelinePairing(unittest.TestCase):
    """cbulge_rows.score_rows (the SNP/indel post-analysis path) uses the same mask."""

    def test_score_rows_sends_the_aligned_pair(self):
        seen = {}

        class _Stub(object):
            @staticmethod
            def CRISPR_BULGE_predict_list(sgs, offs):
                seen.setdefault("calls", []).append((list(sgs), list(offs)))
                return [0.25] * len(sgs)

        crrna = SPACER + "NG-G"
        dna = "ACGTACGTACGTACGTACGTAGAG"
        # target layout the post-analysis uses: [.., crRNA, DNA, .., ref, marker, ..]
        target = ["x", crrna, dna, "chr1", "1", "+", dna, 55, "tail"]
        sys.modules["scorer_runner"] = _Stub
        try:
            out = cbr.score_rows([target], True)
        finally:
            del sys.modules["scorer_runner"]
        self.assertEqual(len(out), 1)
        sgs, offs = seen["calls"][0]
        self.assertEqual(sgs, [crrna], "the gap must reach the model")
        self.assertEqual(offs, [dna])
        self.assertNotIn(ts._CB_SG_DUMMY, sgs, "a 1-bp PAM bulge is not out of domain")


if __name__ == "__main__":
    unittest.main(verbosity=2)
