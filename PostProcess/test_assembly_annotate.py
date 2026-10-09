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


def _gene_bed(tmpdir):
    """A per-transcript GENCODE-shaped BED: the gene identity lives in the
    attributes column, exactly as the shipped gencode.protein_coding.bed.gz does."""
    bed = os.path.join(tmpdir, "genes.bed")
    attrs = ("ID=UTR5:ENST1;Parent=ENST1;gene_id={gid};transcript_id=ENST1;"
             "gene_type=protein_coding;gene_name={name};level=2")
    with open(bed, "w") as f:
        # feature the chr1:100-200 query sits INSIDE -> distance 0
        f.write("chr1\t100\t200\tUTR5:ENST1\t.\t+\tHAVANA\tfive_prime_UTR\t.\t"
                + attrs.format(gid="ENSG0001.1", name="AAA") + "\n")
        # feature 1kb downstream of a chr2:1000 query -> non-zero distance
        f.write("chr2\t5000\t5100\texon:ENST2\t.\t-\tHAVANA\texon\t.\t"
                + attrs.format(gid="ENSG0002.2", name="BBB") + "\n")
    return bed


class TestClosestGene(unittest.TestCase):
    """`--gene_annotation`: a NEAREST-feature query, unlike the overlap query above."""

    def test_mappable_row_gets_name_id_and_distance(self):
        with tempfile.TemporaryDirectory() as tmp:
            df = pd.DataFrame([_row("both", "chr1", 100, pat_target="A" * 23)])
            out = aa.add_closest_gene_columns(df, _gene_bed(tmp))
            self.assertEqual(out.at[0, aa.CLOSEST_GENE_NAME_COL], "AAA")
            self.assertEqual(out.at[0, aa.CLOSEST_GENE_ID_COL], "ENSG0001.1")
            # inside the feature -> zero distance
            self.assertEqual(float(out.at[0, aa.CLOSEST_GENE_DIST_COL]), 0.0)

    def test_distance_is_reported_in_kb(self):
        with tempfile.TemporaryDirectory() as tmp:
            df = pd.DataFrame([_row("both", "chr2", 1000, pat_target="A" * 23)])
            out = aa.add_closest_gene_columns(df, _gene_bed(tmp))
            self.assertEqual(out.at[0, aa.CLOSEST_GENE_NAME_COL], "BBB")
            # ~3977bp away; the column is kb, so well under 10 and not raw bp
            self.assertLess(abs(float(out.at[0, aa.CLOSEST_GENE_DIST_COL])), 10)

    def test_row_without_hg38_coordinate_stays_blank(self):
        """A haplotype-private site has no hg38 equivalent, so "nearest hg38
        gene" is not a statement that can be made about it."""
        with tempfile.TemporaryDirectory() as tmp:
            df = pd.DataFrame([
                _row("both", "chr1", 100, pat_target="A" * 23),
                _row("both_haplotype_private", None, None, pat_target="A" * 23),
            ])
            out = aa.add_closest_gene_columns(df, _gene_bed(tmp))
            self.assertEqual(out.at[0, aa.CLOSEST_GENE_NAME_COL], "AAA")
            for col in aa.CLOSEST_GENE_COLS:
                self.assertTrue(pd.isna(out.at[1, col]), f"{col} must stay blank")

    def test_no_gene_annotation_leaves_frame_unchanged(self):
        df = pd.DataFrame([_row("both", "chr1", 100, pat_target="A" * 23)])
        out = aa.add_closest_gene_columns(df, None)
        self.assertEqual(list(out.columns), list(df.columns))
        for col in aa.CLOSEST_GENE_COLS:
            self.assertNotIn(col, out.columns)

    def test_input_frame_not_mutated(self):
        with tempfile.TemporaryDirectory() as tmp:
            df = pd.DataFrame([_row("both", "chr1", 100, pat_target="A" * 23)])
            before = list(df.columns)
            aa.add_closest_gene_columns(df, _gene_bed(tmp))
            self.assertEqual(list(df.columns), before)

    def test_malformed_distance_is_treated_as_no_annotation(self):
        """`float('')` once crashed a whole integration in the complete-search
        path; a non-numeric trailing field must simply yield no annotation."""
        with tempfile.TemporaryDirectory() as tmp:
            df = pd.DataFrame([_row("both", "chr1", 100, pat_target="A" * 23)])
            empty = os.path.join(tmp, "empty.bed")
            open(empty, "w").close()  # no features at all -> "NA" distance
            out = aa.add_closest_gene_columns(df, empty)
            for col in aa.CLOSEST_GENE_COLS:
                self.assertTrue(pd.isna(out.at[0, col]), f"{col} should be blank")

    def test_attributes_parse_pulls_gene_not_transcript(self):
        name, gid = aa._parse_gene_attributes(
            "ID=UTR5:ENST00000642843.2;Parent=ENST00000642843.2;"
            "gene_id=ENSG00000281518.7;transcript_id=ENST00000642843.2;"
            "gene_type=protein_coding;gene_name=FOXO6"
        )
        self.assertEqual((name, gid), ("FOXO6", "ENSG00000281518.7"))

    def test_intergenic_fills_a_blank_gencode_cell(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = _row("both", "chr2", 1000, pat_target="A" * 23)
            row[aa.GENCODE_SPLIT_COL] = ""  # looked, nothing overlapped
            out = aa.add_closest_gene_columns(pd.DataFrame([row]), _gene_bed(tmp))
            self.assertEqual(out.at[0, aa.GENCODE_SPLIT_COL], aa.INTERGENIC)

    def test_intergenic_never_overwrites_a_real_gencode_label(self):
        """complete-search's ordering makes "intergenic" a FALLBACK (it sets it,
        then overwrites from the overlap labels only when there are any). Same
        precedence here: a populated cell is left alone."""
        with tempfile.TemporaryDirectory() as tmp:
            row = _row("both", "chr2", 1000, pat_target="A" * 23)
            row[aa.GENCODE_SPLIT_COL] = "exon,gene"
            out = aa.add_closest_gene_columns(pd.DataFrame([row]), _gene_bed(tmp))
            self.assertEqual(out.at[0, aa.GENCODE_SPLIT_COL], "exon,gene")

    def test_intergenic_not_applied_at_zero_distance(self):
        """Distance 0 means the site is INSIDE the gene -- not intergenic."""
        with tempfile.TemporaryDirectory() as tmp:
            row = _row("both", "chr1", 100, pat_target="A" * 23)
            row[aa.GENCODE_SPLIT_COL] = ""
            out = aa.add_closest_gene_columns(pd.DataFrame([row]), _gene_bed(tmp))
            self.assertEqual(out.at[0, aa.GENCODE_SPLIT_COL], "")

    def test_closest_gene_columns_match_complete_searchs_names(self):
        """generate_report._COLS resolves Gene/Gene_distance_kb by these exact
        names, so a rename here silently empties those curated columns."""
        self.assertEqual(aa.CLOSEST_GENE_NAME_COL, "Annotation_closest_gene_name")
        self.assertEqual(aa.CLOSEST_GENE_ID_COL, "Annotation_closest_gene_ID")
        self.assertEqual(aa.CLOSEST_GENE_DIST_COL,
                         "Annotation_closest_gene_distance_(kb)")


class TestAnnotationLogNote(unittest.TestCase):
    """The web job runner writes its "no annotation" note to log.txt at job
    start (not submit time), before the subprocess's own output."""

    def _run(self, note):
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        sys.argv = ["x"]
        import pages.main_page as mp
        with tempfile.TemporaryDirectory() as tmp:
            open(os.path.join(tmp, mp.QUEUE_FILE), "w").close()
            mp._run_assembly_search_job("echo from-subprocess", tmp, note)
            with open(os.path.join(tmp, mp.LOG_FILE)) as f:
                return f.read(), os.path.isfile(os.path.join(tmp, mp.QUEUE_FILE))

    def test_note_written_before_subprocess_output_and_queue_cleared(self):
        log, queued = self._run("[web] No hg38 annotation is enabled")
        self.assertEqual(log.splitlines(), ["[web] No hg38 annotation is enabled", "from-subprocess"])
        self.assertFalse(queued)

    def test_no_note_leaves_log_unchanged(self):
        log, _ = self._run("")
        self.assertEqual(log.splitlines(), ["from-subprocess"])


if __name__ == "__main__":
    unittest.main()


class TestAnnotationSplit(unittest.TestCase):
    """The per-kind `Annotation_*` columns assembly-search now emits.

    The report's curated view reads those columns by name (`_COLS`), so emitting
    only the combined `Annotation` string made six annotation cells render "-"
    even when --annotation had run.
    """

    def _oracle(self, annotation_str):
        """resultIntegrator.py's own classification, transcribed.

        Pinned here on purpose: if that block is ever changed, these tests fail
        and whoever changes it has to change both, instead of the two surfaces
        silently describing the same bundle differently.
        """
        personal = set(); encode = set(); gencode = set()
        dhs = set(); cosmic = set(); intogen = set()
        for elem in annotation_str.split(","):
            if "_personal" in elem:
                personal.add(elem.replace("_personal", ""))
            elif "_gencode" in elem:
                gencode.add(elem.replace("_gencode", ""))
            elif "_DHS" in elem:
                dhs.add(elem.replace("_DHS", ""))
            elif "_COSMIC" in elem:
                cosmic.add(elem.replace("_COSMIC", ""))
            elif "_INTOGEN" in elem:
                intogen.add(elem.replace("_INTOGEN", ""))
            else:
                encode.add(elem)
        return {
            "Annotation_personal": ",".join(sorted(personal)),
            "Annotation_GENCODE": ",".join(sorted(gencode)),
            "Annotation_DHS": ",".join(sorted(dhs)),
            "Annotation_COSMIC": ",".join(sorted(cosmic)),
            "Annotation_INTOGEN": ",".join(sorted(intogen)),
            "Annotation_ENCODE": ",".join(sorted(encode)),
        }

    def test_classification_matches_complete_searchs_own(self):
        from annotation import split_annotation_labels
        for case in (
            "dELS,gene_gencode,transcript(+)_gencode",
            "CA-CTCF,Lymphoid_DHS,gene_gencode,transcript(-)_gencode",
            "Tier1_TSG_COSMIC,Stromal_B_DHS,dELS,gene_gencode",
            "SDHA_INTOGEN,gene_gencode,transcript(+)_gencode",
            "exon_gencode,CDS_gencode,gene_gencode,UTR_gencode",
            "Neural_DHS,Primitive_/_embryonic_DHS",
            "PLS,pELS",
        ):
            got, want = split_annotation_labels(case), self._oracle(case)
            for col, value in want.items():
                self.assertEqual(got[col], value, f"{case} -> {col}")

    def test_cosmic_and_intogen_are_not_swallowed_by_the_encode_catchall(self):
        # ENCODE is the else-branch, so these two MUST be tested before it
        from annotation import split_annotation_labels
        got = split_annotation_labels("Tier1_TSG_COSMIC,SDHA_INTOGEN,dELS")
        self.assertEqual(got["Annotation_COSMIC"], "Tier1_TSG")
        self.assertEqual(got["Annotation_INTOGEN"], "SDHA")
        self.assertEqual(got["Annotation_ENCODE"], "dELS")

    def test_unknown_label_is_reported_not_dropped(self):
        # a user-supplied bundle may use labels this split has never seen; they
        # land in the ENCODE catch-all rather than vanishing
        from annotation import split_annotation_labels
        got = split_annotation_labels("SomeCustomLabel,AnotherOne")
        self.assertEqual(got["Annotation_ENCODE"], "AnotherOne,SomeCustomLabel")

    def test_nothing_overlaps_markers_give_every_kind_blank(self):
        from annotation import ANNOTATION_SPLIT_COLS, split_annotation_labels
        for marker in ("n", "", "NA", "nan", None):
            got = split_annotation_labels(marker)
            self.assertTrue(
                all(got[c] == "" for c in ANNOTATION_SPLIT_COLS), f"marker {marker!r}"
            )

    def test_only_kinds_with_a_value_get_a_column(self):
        """An all-blank column would reintroduce the very bug this split fixes,
        and leaving the kind unmapped is what makes the report drop its curated
        column instead of rendering a column of "-"."""
        values = pd.Series(["dELS,gene_gencode", "Neural_DHS", "n"])
        frame = aa.split_annotation_frame(values)
        self.assertIn("Annotation_GENCODE", frame.columns)
        self.assertIn("Annotation_ENCODE", frame.columns)
        self.assertIn("Annotation_DHS", frame.columns)
        # no COSMIC/IntOGen/personal label anywhere -> no such column
        self.assertNotIn("Annotation_COSMIC", frame.columns)
        self.assertNotIn("Annotation_INTOGEN", frame.columns)
        self.assertNotIn("Annotation_personal", frame.columns)

    def test_row_without_an_hg38_coordinate_stays_blank(self):
        # NaN means "not in hg38 space, nothing was looked up" -- distinct from
        # "looked and found nothing", and must not be reported as a feature
        values = pd.Series(["gene_gencode", float("nan")])
        frame = aa.split_annotation_frame(values)
        self.assertEqual(frame["Annotation_GENCODE"].tolist(), ["gene", ""])

    def test_gene_region_is_derived_with_the_shared_helper(self):
        values = pd.Series(["CDS_gencode,gene_gencode", "gene_gencode", "n"])
        frame = aa.split_annotation_frame(values)
        self.assertEqual(
            frame[aa.GENE_REGION_COL].tolist(), ["CDS", "intron", "NA"]
        )

    def test_combined_annotation_column_is_kept_alongside_the_split(self):
        tmp = tempfile.mkdtemp()
        bed = _make_bed(tmp)
        combined = pd.DataFrame({
            "hg38_chr": ["chr1"], "hg38_start": [100],
            "hg38_end_maternal": [123],
            REF_COL + "_paternal": ["A" * 23],
        })
        out = aa.add_annotation_column(combined, bed)
        self.assertIn(aa.ANNOTATION_COL, out.columns)
        # and the split columns sit immediately after it
        cols = list(out.columns)
        self.assertEqual(cols.index(aa.ANNOTATION_COL) + 1,
                         min(cols.index(c) for c in cols
                             if c.startswith("Annotation_")))
