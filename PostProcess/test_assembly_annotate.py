"""Tests for assembly-search's hg38 annotation step
(``PostProcess/assembly_annotate.py``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_assembly_annotate.py' -v
"""

import os
import sys
import tempfile
import unittest

import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import assembly_annotate as aa  # noqa: E402

REF_COL = "Aligned_protospacer+PAM_REF_(fewest_mm+b)"


def _make_bed(tmpdir):
    bed = os.path.join(tmpdir, "ann.bed")
    with open(bed, "w") as f:
        f.write("chr1\t100\t200\texon_gencode\n")
        f.write("chr1\t150\t300\tpELS_encode\n")
        f.write("chr1\t150\t300\tpELS_encode\n")  # duplicate feature name
        f.write("chr2\t1000\t1100\tCTCF_encode\n")
    pysam.tabix_index(bed, force=True, preset="bed")  # writes ann.bed.gz + .tbi
    return bed + ".gz"


def _row(origin, chrom, start, pat_target=None, mat_target=None):
    return {
        "origin": origin, "hg38_chr": chrom, "hg38_start": start,
        f"{REF_COL}_paternal": pat_target, f"{REF_COL}_maternal": mat_target,
        "hg38_end_maternal": None,
    }


class TestAnnotateCombined(unittest.TestCase):
    def test_overlapping_features_sorted_and_deduplicated(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([_row("both", "chr1", 160, "A" * 23)])
            out = aa.add_annotation_column(df, ann)
            self.assertEqual(out.loc[0, "Annotation"], "exon_gencode,pELS_encode")

    def test_no_overlap_is_n_not_blank(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([_row("both", "chr1", 5000, "A" * 23)])
            self.assertEqual(aa.add_annotation_column(df, ann).loc[0, "Annotation"], "n")

    def test_contig_absent_from_annotation_is_n(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([_row("both", "chrUn_x", 10, "A" * 23)])
            self.assertEqual(aa.add_annotation_column(df, ann).loc[0, "Annotation"], "n")

    def test_rows_without_hg38_coordinate_stay_blank(self):
        # "not in hg38" (NaN) must stay distinguishable from "n" (looked, nothing there)
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([
                _row("paternal_non_mappable", None, None, "A" * 23),
                _row("both", "chr1", 160, "A" * 23),
            ])
            out = aa.add_annotation_column(df, ann)
            self.assertTrue(pd.isna(out.loc[0, "Annotation"]))
            self.assertNotEqual(out.loc[1, "Annotation"], "n")

    def test_window_uses_ungapped_target_length(self):
        # feature starts at 195; a 23-bp target at 172 spans 172-195 (miss),
        # but with 4 RNA-bulge gaps removed the same 27-char string is 23 bp
        with tempfile.TemporaryDirectory() as tmp:
            bed = os.path.join(tmp, "w.bed")
            with open(bed, "w") as f:
                f.write("chr1\t195\t210\tedge_feature\n")
            pysam.tabix_index(bed, force=True, preset="bed")
            ann = bed + ".gz"
            gapped = "A" * 12 + "----" + "A" * 11  # 27 chars, 23 bases
            df = pd.DataFrame([
                _row("both", "chr1", 172, gapped),        # covers 172..195 -> miss
                _row("both", "chr1", 173, "A" * 23),      # covers 173..196 -> hit
            ])
            out = aa.add_annotation_column(df, ann)
            self.assertEqual(out.loc[0, "Annotation"], "n")
            self.assertEqual(out.loc[1, "Annotation"], "edge_feature")

    def test_falls_back_to_maternal_target_for_maternal_only_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([_row("maternal_only", "chr2", 1050, None, "A" * 23)])
            self.assertEqual(aa.add_annotation_column(df, ann).loc[0, "Annotation"], "CTCF_encode")

    def test_column_placed_after_hg38_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([_row("both", "chr1", 160, "A" * 23)])
            cols = list(aa.add_annotation_column(df, ann).columns)
            self.assertEqual(cols[cols.index("hg38_end_maternal") + 1], "Annotation")

    def test_no_annotation_file_leaves_frame_unchanged(self):
        df = pd.DataFrame([_row("both", "chr1", 160, "A" * 23)])
        out = aa.add_annotation_column(df, None)
        self.assertNotIn("Annotation", out.columns)

    def test_input_frame_not_mutated(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = _make_bed(tmp)
            df = pd.DataFrame([_row("both", "chr1", 160, "A" * 23)])
            aa.add_annotation_column(df, ann)
            self.assertNotIn("Annotation", df.columns)


if __name__ == "__main__":
    unittest.main()
