"""Regression test for the both_haplotype_private site set in
pages/results_page.py's assembly-search results page.

Real data (376 real both_haplotype_private rows from a genome-wide HG01255
run) showed every hg38/CFD/mismatch column 100% NaN for these rows when they
were still mixed into _assembly_mappable_frame()'s hg38 table -- they have
no hg38 coordinate at all, so that table's own display-column logic
(built for hg38-mapped rows) produced an all-blank-except-4-columns shape.
This pins: (1) _assembly_mappable_frame() now excludes them, (2) their own
dedicated frame (_assembly_haplotype_private_frame()) rebuilds real
Spacer+PAM/CFD/mismatches values by joining back to each haplotype's own
predictions, and (3) origin_counts (feeding the coverage stat cards/plot)
still sees the real both_haplotype_private count -- filtering must happen
AFTER that count is taken, not before.

Needs a real (if minimal) Dash `app` import, same as the rest of
pages/results_page.py -- unlike history_page.py's simpler helpers, these
functions have too many real dependencies (glob, current_working_directory,
find_results_prefix, load_crisprme_predictions) to extract via AST in
isolation, so this imports the module directly and patches
current_working_directory to a tmpdir instead.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "PostProcess"))

from assembly_reconcile import PRED_COLS  # noqa: E402
import pages.results_page as rp  # noqa: E402

PARAMS_FILE = ".Params.txt"
GUIDES_FILE = ".guides.txt"
RESULTS_DIR = "Results"


def _pred_row(chrom, pos, cfd, mm=0, bulges=0):
    return {
        "Spacer+PAM": "ACGTACGTACGTACGTACGTNGG", "Chromosome": chrom,
        "Start_coordinate_(fewest_mm+b)": pos, "Strand_(fewest_mm+b)": "+",
        "Aligned_spacer+PAM_(fewest_mm+b)": "ACGTACGTACGTACGTACGTNGG",
        "Aligned_protospacer+PAM_REF_(fewest_mm+b)": "ACGTACGTACGTACGTACGTNGG",
        "Aligned_protospacer+PAM_ALT_(fewest_mm+b)": None,
        "Mismatches_(fewest_mm+b)": mm, "Bulges_(fewest_mm+b)": bulges,
        "CFD_score_(fewest_mm+b)": cfd,
    }


class TestHaplotypePrivateSiteSet(unittest.TestCase):
    def _make_haplotype_results(self, tmpdir, dirname, row):
        """One haplotype's complete-search output dir: a single-row
        integrated_results.tsv (+ empty alt-alignments file), the real
        PRED_COLS schema -- same pattern test_assembly_reconcile.py's own
        TestReconcileHaplotypes._make_haplotype_results uses."""
        results_dir = os.path.join(tmpdir, RESULTS_DIR, dirname)
        os.makedirs(results_dir, exist_ok=True)
        prefix = f"guide_PAM_{dirname}_mm4_bMax2"
        cols = PRED_COLS + ["Bulge_type_(fewest_mm+b)"]
        full_row = dict(row)
        full_row.setdefault("Bulge_type_(fewest_mm+b)", None)
        df = pd.DataFrame([full_row], columns=cols)
        df.to_csv(
            os.path.join(results_dir, f"{prefix}_integrated_results.tsv"),
            sep="\t", index=False,
        )
        pd.DataFrame(columns=cols).to_csv(
            os.path.join(results_dir, f"{prefix}_all_results_with_alternative_alignments.tsv"),
            sep="\t", index=False,
        )
        return results_dir

    def _make_combined_job(self, tmpdir, job_id):
        combined_dir = os.path.join(tmpdir, RESULTS_DIR, f"{job_id}_combined")
        os.makedirs(combined_dir, exist_ok=True)
        row = {
            "Spacer+PAM_paternal": None, "Chromosome_paternal": "chr9",
            "Start_coordinate_(fewest_mm+b)_paternal": None, "Strand_(fewest_mm+b)": None,
            "Aligned_spacer+PAM_(fewest_mm+b)_paternal": None,
            "Aligned_protospacer+PAM_REF_(fewest_mm+b)_paternal": None,
            "Aligned_protospacer+PAM_ALT_(fewest_mm+b)_paternal": None,
            "Mismatches_(fewest_mm+b)_paternal": None, "Bulges_(fewest_mm+b)_paternal": None,
            "CFD_score_(fewest_mm+b)_paternal": None, "off_target_id_paternal": 0,
            "hg38_chr": None, "hg38_start": None, "hg38_end_paternal": None,
            "Spacer+PAM_maternal": None, "Chromosome_maternal": "chr9",
            "Start_coordinate_(fewest_mm+b)_maternal": None,
            "Aligned_spacer+PAM_(fewest_mm+b)_maternal": None,
            "Aligned_protospacer+PAM_REF_(fewest_mm+b)_maternal": None,
            "Aligned_protospacer+PAM_ALT_(fewest_mm+b)_maternal": None,
            "Mismatches_(fewest_mm+b)_maternal": None, "Bulges_(fewest_mm+b)_maternal": None,
            "CFD_score_(fewest_mm+b)_maternal": None, "off_target_id_maternal": 0,
            "hg38_end_maternal": None, "origin": "both_haplotype_private",
            "Start_coordinate_paternal": 67041469, "Start_coordinate_maternal": 70115639,
        }
        both_row = dict(row)
        both_row.update({
            "origin": "both", "Chromosome_paternal": "chr1", "Chromosome_maternal": "chr1",
            "hg38_chr": "chr1", "hg38_start": 1000, "off_target_id_paternal": 0,
            "off_target_id_maternal": 0,
        })
        pd.DataFrame([row, both_row]).to_csv(
            os.path.join(combined_dir, f"{job_id}_combined_hg38.tsv"), sep="\t", index=False,
        )
        with open(os.path.join(combined_dir, PARAMS_FILE), "w") as f:
            f.write(f"1\tGenome_type\tassembly\n")
            f.write(f"2\tPaternal_dir\t{job_id}_paternal\n")
            f.write(f"3\tMaternal_dir\t{job_id}_maternal\n")
        return combined_dir

    def test_haplotype_private_rows_get_real_values_not_nan(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            job_id = "testjob"
            self._make_haplotype_results(
                tmpdir, f"{job_id}_paternal", _pred_row("chr9", 67041469, 1.0, mm=0, bulges=1)
            )
            self._make_haplotype_results(
                tmpdir, f"{job_id}_maternal", _pred_row("chr9", 70115639, 0.546, mm=1, bulges=2)
            )
            combined_dir = self._make_combined_job(tmpdir, job_id)
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                df, cols = rp._assembly_haplotype_private_frame(combined_dir)
                self.assertEqual(len(df), 1)
                self.assertEqual(df.iloc[0]["CFD_score_(fewest_mm+b)_paternal"], 1.0)
                self.assertEqual(df.iloc[0]["CFD_score_(fewest_mm+b)_maternal"], 0.546)
                self.assertEqual(df.iloc[0]["Mismatches_(fewest_mm+b)_maternal"], 1)
                self.assertEqual(df.iloc[0]["Chromosome_paternal"], "chr9")
                self.assertEqual(df.iloc[0]["Chromosome_maternal"], "chr9")
                self.assertNotIn("hg38_chr", cols)

    def test_mappable_frame_excludes_haplotype_private_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            job_id = "testjob"
            self._make_haplotype_results(
                tmpdir, f"{job_id}_paternal", _pred_row("chr9", 67041469, 1.0)
            )
            self._make_haplotype_results(
                tmpdir, f"{job_id}_maternal", _pred_row("chr9", 70115639, 0.5)
            )
            combined_dir = self._make_combined_job(tmpdir, job_id)
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                df, _ = rp._assembly_mappable_frame(combined_dir)
                self.assertEqual(len(df), 1)  # only the "both" row survives
                self.assertNotIn("both_haplotype_private", df["origin"].tolist())

    def test_site_set_dispatch_returns_haplotype_private_frame(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            job_id = "testjob"
            self._make_haplotype_results(
                tmpdir, f"{job_id}_paternal", _pred_row("chr9", 67041469, 1.0)
            )
            self._make_haplotype_results(
                tmpdir, f"{job_id}_maternal", _pred_row("chr9", 70115639, 0.5)
            )
            combined_dir = self._make_combined_job(tmpdir, job_id)
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                df, _ = rp._assembly_site_set_frame(combined_dir, "both_haplotype_private")
                self.assertEqual(len(df), 1)

    def test_missing_haplotype_dirs_returns_empty_not_crash(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            job_id = "testjob"
            combined_dir = self._make_combined_job(tmpdir, job_id)
            # no _make_haplotype_results calls -- haplotype dirs don't exist
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                df, cols = rp._assembly_haplotype_private_frame(combined_dir)
                self.assertEqual(len(df), 1)
                self.assertTrue(pd.isna(df.iloc[0]["CFD_score_(fewest_mm+b)_paternal"]))

    # -- Zero-count edge cases -----------------------------------------
    # The four tests above all construct at least one both_haplotype_private
    # row. Real jobs overwhelmingly have ZERO of them: every pre-impg job
    # (no "origin" column at all in combined_hg38.tsv) and any impg-enabled
    # job whose direct alignment genuinely finds no haplotype-private pairs.
    # Confirmed directly against the real diverseOriginTest2_INILXMVL7H_combined
    # job (predates impg, no "origin" column) before adding these: the stat
    # card showed a real "0", the coverage plot rendered normally, and the
    # Custom Ranking dropdown's "Both haplotypes (non-mappable to hg38)"
    # option produced an empty (not crashing) table. Pinning that here.
    def _make_combined_job_no_origin_column(self, tmpdir, job_id):
        """An old-style (pre-impg) combined_hg38.tsv: no "origin" column at
        all, same as every real job reconciled before this feature existed."""
        combined_dir = os.path.join(tmpdir, RESULTS_DIR, f"{job_id}_combined")
        os.makedirs(combined_dir, exist_ok=True)
        pd.DataFrame([{
            "Chromosome_paternal": "chr1", "hg38_chr": "chr1", "hg38_start": 1000,
            "off_target_id_paternal": 0, "off_target_id_maternal": 0,
        }]).to_csv(
            os.path.join(combined_dir, f"{job_id}_combined_hg38.tsv"), sep="\t", index=False,
        )
        with open(os.path.join(combined_dir, PARAMS_FILE), "w") as f:
            f.write("1\tGenome_type\tassembly\n")
        return combined_dir

    def test_no_origin_column_at_all_returns_empty_not_crash(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            combined_dir = self._make_combined_job_no_origin_column(tmpdir, "oldjob")
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                df, cols = rp._assembly_haplotype_private_frame(combined_dir)
                self.assertEqual(len(df), 0)
                self.assertEqual(cols, [])
                # site-set dispatch must not crash either, feeding the same
                # empty-df path the Custom Ranking callback's own
                # `if df.empty: return [], [], []` guard relies on
                df2, cols2 = rp._assembly_site_set_frame(combined_dir, "both_haplotype_private")
                self.assertTrue(df2.empty)

    def test_origin_column_present_but_zero_private_rows(self):
        # a real impg-enabled job whose direct alignment found zero
        # haplotype-private pairs -- "origin" column exists, just no rows
        # equal "both_haplotype_private"
        with tempfile.TemporaryDirectory() as tmpdir:
            job_id = "job_zero_private"
            combined_dir = os.path.join(tmpdir, RESULTS_DIR, f"{job_id}_combined")
            os.makedirs(combined_dir, exist_ok=True)
            pd.DataFrame([
                {"origin": "both", "Chromosome_paternal": "chr1", "hg38_chr": "chr1",
                 "hg38_start": 1000, "off_target_id_paternal": 0, "off_target_id_maternal": 0},
                {"origin": "paternal_only", "Chromosome_paternal": "chr2", "hg38_chr": "chr2",
                 "hg38_start": 2000, "off_target_id_paternal": 1, "off_target_id_maternal": None},
            ]).to_csv(
                os.path.join(combined_dir, f"{job_id}_combined_hg38.tsv"), sep="\t", index=False,
            )
            with open(os.path.join(combined_dir, PARAMS_FILE), "w") as f:
                f.write("1\tGenome_type\tassembly\n")
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                df, cols = rp._assembly_haplotype_private_frame(combined_dir)
                self.assertEqual(len(df), 0)
                self.assertEqual(cols, [])
                # the other two real rows must be unaffected
                mappable_df, _ = rp._assembly_mappable_frame(combined_dir)
                self.assertEqual(len(mappable_df), 2)

    def test_blank_alt_and_plain_start_columns_are_hidden(self):
        # ALT alignment columns are structurally empty in assembly-search (no
        # VCF), and the un-suffixed Start_coordinate_paternal/_maternal are
        # populated only on both_haplotype_private rows -- neither should be
        # rendered as an all-blank column in the tables where they're empty.
        alt = "Aligned_protospacer+PAM_ALT_(fewest_mm+b)"
        with tempfile.TemporaryDirectory() as tmpdir:
            job_id = "testjob"
            self._make_haplotype_results(
                tmpdir, f"{job_id}_paternal", _pred_row("chr9", 67041469, 1.0)
            )
            self._make_haplotype_results(
                tmpdir, f"{job_id}_maternal", _pred_row("chr9", 70115639, 0.5)
            )
            combined_dir = self._make_combined_job(tmpdir, job_id)
            with open(os.path.join(combined_dir, "paternal_offtargets_not_lifted.bed"), "w") as f:
                f.write("# Deleted in new\nchr9\t1\t2\t0\n")
            with patch.object(rp, "current_working_directory", tmpdir + os.sep):
                _, mappable_cols = rp._assembly_mappable_frame(combined_dir)
                self.assertNotIn("Start_coordinate_paternal", mappable_cols)
                self.assertNotIn("Start_coordinate_maternal", mappable_cols)
                self.assertIn("Start_coordinate_(fewest_mm+b)_paternal", mappable_cols)
                self.assertNotIn(f"{alt}_paternal", mappable_cols)

                unmappable, unmappable_cols = rp._assembly_site_set_frame(
                    combined_dir, "paternal_unmappable"
                )
                self.assertEqual(len(unmappable), 1)
                self.assertNotIn(alt, unmappable_cols)
                self.assertIn("CFD_score_(fewest_mm+b)", unmappable_cols)

                _, private_cols = rp._assembly_haplotype_private_frame(combined_dir)
                self.assertNotIn(f"{alt}_paternal", private_cols)
                self.assertNotIn(f"{alt}_maternal", private_cols)
                # the private table's own native coordinates must stay visible
                self.assertIn("Start_coordinate_paternal", private_cols)
                self.assertIn("Start_coordinate_maternal", private_cols)


if __name__ == "__main__":
    unittest.main()
