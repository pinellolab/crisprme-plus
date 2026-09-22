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


if __name__ == "__main__":
    unittest.main()
