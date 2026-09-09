"""Unit test for synth_sites_gt: sites-only VCF -> synthetic presence-GT VCF that
the CRISPRitz enricher will accept for fake-indel genome construction. stdlib only."""
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_sites_gt import synth  # noqa: E402


class TestSynthSitesGT(unittest.TestCase):
    def test_adds_presence_gt_and_preserves_sites(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "sites.vcf.gz")
            out = os.path.join(d, "gt.vcf.gz")
            with gzip.open(inp, "wt") as w:
                w.write("##fileformat=VCFv4.2\n")
                w.write('##INFO=<ID=AF_max,Number=A,Type=Float,Description="x">\n')
                w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                w.write("chr22\t10510343\t.\tCTAA\tC\t.\tPASS\tAF_max=0.015;AF=0.015\n")
                w.write("chr22\t10510355\t.\tAT\tA\t.\tPASS\tAF_max=0.012\n")
            n = synth(inp, out, "TEST")
            self.assertEqual(n, 2)
            body = gzip.open(out, "rt").read()
            self.assertIn("##FORMAT=<ID=GT", body)          # FORMAT header injected
            self.assertIn("\tFORMAT\tTEST", body)            # one synthetic sample col
            self.assertEqual(body.count("\tGT\t0/1"), 2)     # presence GT on every record
            self.assertIn("AF_max=0.015", body)              # INFO preserved verbatim
            row = [l for l in body.splitlines()
                   if l.startswith("chr22\t10510343")][0].split("\t")
            # the 8 fixed columns (incl. the indel REF/ALT) are untouched
            self.assertEqual(row[3], "CTAA")
            self.assertEqual(row[4], "C")
            self.assertEqual(row[8], "GT")
            self.assertEqual(row[9], "0/1")


if __name__ == "__main__":
    unittest.main()
