"""Regression tests for PostProcess/change_headers_bestMerge.py.

The script converts the 3-scorer `.bestMerge.txt` (CFD + MMBLG + CRISPR_BULGE
blocks, emitted by the submit_job pr-join) into the per-cluster web display /
download format. Its fixed `new_order` projection predates the CRISPR-Bulge
scorer and would SILENTLY DROP the 24 CRISPR_BULGE_* columns (a pandas column
selection, not a KeyError) -- the per-cluster download lost every CRISPR-Bulge
column with no error. These tests lock in that the block is preserved, and that
a CFD-only input (older runs, no CRISPR_BULGE columns) is unaffected.
"""
import os
import subprocess
import sys
import tempfile
import unittest

try:
    import pandas  # noqa: F401
    _HAVE_PANDAS = True
except Exception:
    _HAVE_PANDAS = False

_SCRIPT = os.path.join(os.path.dirname(__file__), "change_headers_bestMerge.py")

# the exact 48-column new_order the script selects (CFD + MMBLG blocks)
_NEW_ORDER = [
    "Real_Guide", "Chromosome", "Position", "Direction", "Cluster_Position", "crRNA",
    "Reference", "DNA", "Mismatches", "Bulge_Size", "Total", "#Bulge_type", "PAM_gen",
    "CFD", "CFD_ref", "Highest_CFD_Risk_Score", "Highest_CFD_Absolute_Risk_Score",
    "SNP", "AF", "rsID", "Samples", "Var_uniq", "#Seq_in_cluster",
    "MMBLG_Real_Guide", "MMBLG_Chromosome", "MMBLG_Position", "MMBLG_Cluster_Position",
    "MMBLG_Direction", "MMBLG_crRNA", "MMBLG_Reference", "MMBLG_DNA", "MMBLG_Mismatches",
    "MMBLG_Bulge_Size", "MMBLG_Total", "MMBLG_#Bulge_type", "MMBLG_PAM_gen", "MMBLG_CFD",
    "MMBLG_CFD_ref", "MMBLG_CFD_Risk_Score", "MMBLG_CFD_Absolute_Risk_Score",
    "MMBLG_Var_uniq", "MMBLG_SNP", "MMBLG_AF", "MMBLG_rsID", "MMBLG_Samples",
    "MMBLG_#Seq_in_cluster", "MMBLG_Annotation_Type", "Annotation_Type",
]
# the CRISPR-Bulge representative-alignment block appended by the 3-way bestMerge join
_CB_COLS = ["CRISPR_BULGE_CFD", "CRISPR_BULGE_Samples", "CRISPR_BULGE_rsID", "CRISPR_BULGE_AF"]


@unittest.skipUnless(_HAVE_PANDAS, "pandas required")
class TestChangeHeadersBestMerge(unittest.TestCase):
    def _run(self, header, row):
        d = tempfile.mkdtemp(prefix="chm_test_")
        self.addCleanup(lambda: __import__("shutil").rmtree(d, ignore_errors=True))
        src = os.path.join(d, "in.txt")
        out = os.path.join(d, "out.txt")
        with open(src, "w") as h:
            h.write("\t".join(header) + "\n")
            h.write("\t".join(row) + "\n")
        res = subprocess.run([sys.executable, _SCRIPT, src, out], capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        with open(out) as h:
            out_header = h.readline().rstrip("\n").split("\t")
            out_row = h.readline().rstrip("\n").split("\t")
        return out_header, out_row

    def test_crispr_bulge_block_is_preserved(self):
        header = _NEW_ORDER + _CB_COLS
        # distinct sentinel per column so we can prove the CB values survive intact
        row = [f"v_{c}" for c in header]
        # the REF/ALT swap keys on == "NA"; keep them non-NA so it's a passthrough
        row[header.index("Reference")] = "ACGT"
        row[header.index("DNA")] = "ACGA"
        row[header.index("MMBLG_Reference")] = "ACGT"
        row[header.index("MMBLG_DNA")] = "ACGA"
        out_header, out_row = self._run(header, row)
        # every CRISPR_BULGE column must still be present (the bug dropped all of them)
        for c in _CB_COLS:
            self.assertIn(c, out_header, f"{c} was silently dropped from the per-cluster output")
        # and carry its original value (not blanked / shuffled)
        self.assertEqual(out_row[out_header.index("CRISPR_BULGE_CFD")], "v_CRISPR_BULGE_CFD")
        # the CFD display transform still happened (sanity)
        self.assertIn("CFD_score_(highest_CFD)", out_header)

    def test_cfd_only_input_is_unaffected(self):
        # an older run with no CRISPR_BULGE columns -> no CB columns in output, no crash
        header = list(_NEW_ORDER)
        row = [f"v_{c}" for c in header]
        row[header.index("Reference")] = "ACGT"
        row[header.index("DNA")] = "ACGA"
        row[header.index("MMBLG_Reference")] = "ACGT"
        row[header.index("MMBLG_DNA")] = "ACGA"
        out_header, _ = self._run(header, row)
        self.assertFalse([c for c in out_header if c.startswith("CRISPR_BULGE")])
        self.assertIn("CFD_score_(highest_CFD)", out_header)


if __name__ == "__main__":
    unittest.main()
