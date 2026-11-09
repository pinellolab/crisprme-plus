#!/usr/bin/env python3
"""prepare_mega_build_vcf.py: the build input of the sites-only mega index, from its merged
VCF: the records CRISPRitz's enricher reads (FILTER PASS or '.'), each with INFO/AF its
AF_max, every other byte as it was. STDLIB only, but for the BGZF writer (pysam).

Run: cd PostProcess && python -m unittest test_prepare_mega_build_vcf -v
"""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prepare_mega_build_vcf as pm  # noqa: E402

try:
    import pysam  # noqa: F401
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False

HEADER = ("##fileformat=VCFv4.2\n"
          '##INFO=<ID=AF_max,Number=A,Type=Float,Description="max">\n'
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
ROWS = [
    "chr22\t100\trs1\tA\tG\t.\tPASS\tAF_HGDP=0.1;AF_max=0.1",
    "chr22\t200\trs2\tC\tT\t.\tLOW_VQSLOD\tAF_HGDP=0.2;AF_max=0.2",
    "chr22\t300\t.\tCTT\tC\t.\t.\tAF_gnomAD=0.03;AF_max=0.03",
    "chr22\t400\t.\tG\tA,GT\t.\tPASS\tAF=0.9;AF_HGDP=0.1,0.2;AF_max=0.1,0.2",
    "chr22\t500\t.\tT\tC\t.\tExcHet\tAF_HGDP=0.5;AF_max=0.5",
    "chr22\t600\t.\tT\tTA\t.\tPASS\tAF_HGDP=.;AF_max=.",
]


class TestPrepare(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, self.d, True)
        self.inp = os.path.join(self.d, "mega.chr22.afmax.vcf.gz")
        with gzip.open(self.inp, "wt") as fh:
            fh.write(HEADER + "\n".join(ROWS) + "\n")

    def read(self, path):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            return fh.read().splitlines()

    def test_pass_records_with_af_max_as_af(self):
        out = os.path.join(self.d, "build.vcf")
        stats = pm.prepare(self.inp, out)
        lines = self.read(out)
        body = [l.split("\t") for l in lines if not l.startswith("#")]
        self.assertEqual([r[1] for r in body], ["100", "300", "400", "600"])
        self.assertEqual([r[7] for r in body], [
            "AF_HGDP=0.1;AF_max=0.1;AF=0.1",
            "AF_gnomAD=0.03;AF_max=0.03;AF=0.03",
            "AF_HGDP=0.1,0.2;AF_max=0.1,0.2;AF=0.1,0.2",   # the AF it had is replaced
            "AF_HGDP=.;AF_max=.",                          # no AF_max: no AF
        ])
        self.assertEqual([r[:7] for r in body], [ln.split("\t")[:7] for ln in
                                                 (ROWS[0], ROWS[2], ROWS[3], ROWS[5])])
        self.assertEqual(stats, {"records": 6, "kept": 4, "filtered": 2, "af": 3,
                                 "no_af_max": 1})
        header = [l for l in lines if l.startswith("##")]
        self.assertIn(pm.AF_HEADER, header)
        self.assertIn(pm.BUILD_MARK, header)
        self.assertTrue(lines[len(header)].startswith("#CHROM"))

    def test_an_af_header_in_place_is_kept_once(self):
        with gzip.open(self.inp, "wt") as fh:
            fh.write(HEADER.replace("#CHROM", '##INFO=<ID=AF,Number=A,Type=Float,'
                                              'Description="x">\n#CHROM')
                     + ROWS[0] + "\n")
        out = os.path.join(self.d, "build.vcf")
        pm.prepare(self.inp, out)
        self.assertEqual(sum(1 for l in self.read(out) if l.startswith("##INFO=<ID=AF,")), 1)

    @unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF) not installed")
    def test_command_writes_bgzf(self):
        out = os.path.join(self.d, "build.vcf.gz")
        done = subprocess.run([sys.executable, os.path.join(HERE, "prepare_mega_build_vcf.py"),
                               self.inp, out], capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("6 records -> 4 kept", done.stderr)
        with open(out, "rb") as fh:
            self.assertEqual(fh.read(4), b"\x1f\x8b\x08\x04")  # BGZF, for tabix
        self.assertEqual(sum(1 for l in self.read(out) if not l.startswith("#")), 4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
