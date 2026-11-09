"""Unit test for synth_sites_gt: sites-only VCF -> synthetic presence-GT VCF that
the CRISPRitz enricher will accept for fake-indel genome construction. stdlib only."""
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_sites_gt import synth  # noqa: E402

try:
    import pysam  # noqa: F401  (the BGZF writer)
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False


@unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF) not installed")
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
            # the .vcf.gz MUST be BGZF (blocked gzip) so the enricher recipe can
            # tabix-index it -- a plain gzip.open() stream is rejected by tabix
            # ("not BGZF"). BGZF sets the gzip FEXTRA flag (byte 3 == 0x04) and
            # carries a "BC" subfield; plain gzip has byte 3 == 0x00.
            with open(out, "rb") as fh:
                head = fh.read(18)
            self.assertEqual(head[:4], b"\x1f\x8b\x08\x04", "output is not BGZF")
            self.assertIn(b"BC", head, "BGZF 'BC' subfield missing")
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




class TestEveryAltOfALine(unittest.TestCase):
    def test_every_alt_of_a_line_is_carried_and_the_header_says_so(self):
        # CRISPRitz's enricher windows an indel a sample carries: a presence genotype of
        # 0/1 on a line listing several ALTs carried its first alone, so the others got
        # no fake-indel window and were never searched. And the VCF's one sample is no
        # individual: the header line tells build-index-only so (crisprme.py reads the
        # dataset as sites-only).
        import synth_sites_gt as sg
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "sites.vcf")
            out = os.path.join(d, "gt.vcf")
            with open(inp, "w") as w:
                w.write("##fileformat=VCFv4.2\n")
                w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                w.write("chr22\t100\t.\tCA\tC,CAA\t.\tPASS\tAF_max=0.1,0.2\n")
                w.write("chr22\t200\t.\tG\tA,C,T\t.\tPASS\tAF_max=0.1,0.2,0.3\n")
                w.write("chr22\t300\t.\tAT\tA\t.\tPASS\tAF_max=0.1\n")
            self.assertEqual(synth(inp, out, "MEGA"), 3)
            with open(out) as fh:
                body = fh.read().splitlines()
            self.assertIn(sg.SITES_ONLY_MARK, body)
            self.assertLess(body.index(sg.SITES_ONLY_MARK),
                            next(i for i, l in enumerate(body) if l.startswith("#CHROM")))
            gts = [l.split("\t")[9] for l in body if not l.startswith("#")]
            self.assertEqual(gts, ["1/2", "1/2/3", "0/1"])


if __name__ == "__main__":
    unittest.main()
