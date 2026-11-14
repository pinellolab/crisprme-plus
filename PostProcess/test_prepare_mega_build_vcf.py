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
        self.assertEqual(stats, {"records": 6, "kept": 4, "filtered": 2, "collapsed": 0,
                                 "alts_collapsed": 0, "af": 3, "no_af_max": 1,
                                 "af_max_derived": 0})
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

    def test_one_record_of_each_ref_alt(self):
        """Two PASS records of one REF>ALT are one allele: the build input keeps the more
        frequent (INFO/AF, set to AF_max; on a tie the earlier), as the enricher merges two
        SNP records of one allele. The enricher gives each indel record a fake-indel
        window and an indel-log row of its own, while the stores hold an allele once: with
        both records the recipe's checks failed (chr21: "1 allele(s) on more than one row",
        review r1 of lane M, M-R1-3). A record of another FILTER is out first; a line of
        several ALTs loses only the ALT another record has, with its AF values."""
        deletion = "AGGTGTAATTCCAGTGCTGG"
        rows = [
            # the as-run merge's chr21 case: one deletion on two lines, one AF_max
            "chr21\t45918284\t.\t%s\tA\t32075\tPASS\tAF_HGDP=0.00271739;AF_max=0.00271739"
            % deletion,
            "chr21\t45918284\t.\t%s\tA\t9066.62\tPASS\tAF_HGDP=0.00135135;AF_max=0.00271739"
            % deletion,
            # an SNV on two lines, the second more frequent; a third of another FILTER
            "chr21\t45918300\trs1\tC\tT\t.\tPASS\tAF_HGDP=0.1;AF_max=0.1",
            "chr21\t45918300\trs2\tC\tT\t.\t.\tAF_gnomAD=0.3;AF_max=0.3",
            "chr21\t45918300\trs3\tC\tT\t.\tLOW_VQSLOD\tAF_HGDP=0.9;AF_max=0.9",
            # a line of two ALTs, one of them on a more frequent line of its own
            "chr21\t45918400\t.\tG\tA,GTT\t.\tPASS\tAF_HGDP=0.2,0.01;AF_max=0.2,0.01",
            "chr21\t45918400\t.\tG\tGTT\t.\tPASS\tAF_HGDP=0.05;AF_max=0.05",
            # one REF>ALT at two positions: two alleles
            "chr21\t45918500\t.\tT\tTA\t.\tPASS\tAF_HGDP=0.1;AF_max=0.1",
            "chr21\t45918501\t.\tT\tTA\t.\tPASS\tAF_HGDP=0.1;AF_max=0.1",
        ]
        with gzip.open(self.inp, "wt") as fh:
            fh.write(HEADER.replace("#CHROM", '##INFO=<ID=AF_HGDP,Number=A,Type=Float,'
                                              'Description="x">\n#CHROM')
                     + "\n".join(rows) + "\n")
        out = os.path.join(self.d, "build.vcf")
        stats = pm.prepare(self.inp, out)
        body = [l.split("\t") for l in self.read(out) if not l.startswith("#")]
        self.assertEqual([(r[1], r[2], r[3][:3], r[4], r[5], r[7]) for r in body], [
            ("45918284", ".", "AGG", "A", "32075",
             "AF_HGDP=0.00271739;AF_max=0.00271739;AF=0.00271739"),
            ("45918300", "rs2", "C", "T", ".", "AF_gnomAD=0.3;AF_max=0.3;AF=0.3"),
            ("45918400", ".", "G", "A", ".", "AF_HGDP=0.2;AF_max=0.2;AF=0.2"),
            ("45918400", ".", "G", "GTT", ".", "AF_HGDP=0.05;AF_max=0.05;AF=0.05"),
            ("45918500", ".", "T", "TA", ".", "AF_HGDP=0.1;AF_max=0.1;AF=0.1"),
            ("45918501", ".", "T", "TA", ".", "AF_HGDP=0.1;AF_max=0.1;AF=0.1"),
        ])
        self.assertEqual(stats, {"records": 9, "kept": 6, "filtered": 1, "collapsed": 2,
                                 "alts_collapsed": 1, "af": 6, "no_af_max": 0,
                                 "af_max_derived": 0})
        # every store and check then reads one record of each allele
        import build_mega_registry as bmr
        import build_mega_indel_af as bmi
        st = {}
        self.assertEqual(len(list(bmr.read_as_enricher(out, "indel", {"HGDP", "gnomAD"},
                                                       stats=st))), 4)
        self.assertEqual(st["merged"], 0)
        self.assertEqual([r[:3] for r in bmi.iter_indel_af_records(out, {"HGDP"})],
                         [("45918284", deletion, "A"), ("45918400", "G", "GTT"),
                          ("45918500", "T", "TA"), ("45918501", "T", "TA")])

    def test_a_record_differing_from_another_only_by_case_is_one_allele(self):
        """CRISPRitz 2.10.0 reads REF and every ALT in upper case, as VCF 4.x does (its
        docs/MULTIALLELIC_ENRICHMENT.md sections 2.1 and 2.3): a deletion written CTT>C on
        one line and ctt>c on another is one allele, which the enricher would window twice
        (two log rows of one allele: the recipe's check fails), and an SNV written C>G and
        c>g is one dictionary entry. The build input keeps one record of each, the more
        frequent, written as the input has it (the enricher upper-cases it); a line of
        several ALTs loses only the one another line has in another case. Keyed as written,
        both records of each stayed."""
        rows = [
            "chr21	1000	.	CTT	C	.	PASS	AF_gnomAD=0.01;AF_max=0.01",
            "chr21	1000	.	ctt	c	.	PASS	AF_gnomAD=0.02;AF_max=0.02",
            "chr21	1100	rs1	C	G	.	PASS	AF_HGDP=0.3;AF_max=0.3",
            "chr21	1100	rs2	c	g	.	.	AF_HGDP=0.1;AF_max=0.1",
            "chr21	1200	.	A	g,AT	.	PASS	AF_HGDP=0.05,0.2;AF_max=0.05,0.2",
            "chr21	1200	.	A	G	.	PASS	AF_HGDP=0.4;AF_max=0.4",
            # a symbolic allele is an ID, read as written: two alleles
            "chr21	1300	.	A	<DEL>	.	PASS	AF_HGDP=0.1;AF_max=0.1",
            "chr21	1300	.	a	<del>	.	PASS	AF_HGDP=0.1;AF_max=0.1",
        ]
        with gzip.open(self.inp, "wt") as fh:
            fh.write(HEADER.replace("#CHROM", '##INFO=<ID=AF_HGDP,Number=A,Type=Float,'
                                              'Description="x">\n#CHROM')
                     + "\n".join(rows) + "\n")
        out = os.path.join(self.d, "build.vcf")
        stats = pm.prepare(self.inp, out)
        body = [l.split("\t") for l in self.read(out) if not l.startswith("#")]
        self.assertEqual([(r[1], r[2], r[3], r[4], r[7]) for r in body], [
            ("1000", ".", "ctt", "c", "AF_gnomAD=0.02;AF_max=0.02;AF=0.02"),
            ("1100", "rs1", "C", "G", "AF_HGDP=0.3;AF_max=0.3;AF=0.3"),
            ("1200", ".", "A", "AT", "AF_HGDP=0.2;AF_max=0.2;AF=0.2"),
            ("1200", ".", "A", "G", "AF_HGDP=0.4;AF_max=0.4;AF=0.4"),
            ("1300", ".", "A", "<DEL>", "AF_HGDP=0.1;AF_max=0.1;AF=0.1"),
            ("1300", ".", "a", "<del>", "AF_HGDP=0.1;AF_max=0.1;AF=0.1"),
        ])
        self.assertEqual((stats["collapsed"], stats["alts_collapsed"], stats["kept"]),
                         (2, 1, 6))
        # every store then reads one record of each allele, in upper case
        import build_mega_registry as bmr
        st = {}
        self.assertEqual([r[:3] for r in bmr.read_as_enricher(out, "indel", {"gnomAD"},
                                                              stats=st)],
                         [(1000, "CTT", "C")])
        self.assertEqual(st["merged"], 0)
        st = {}
        self.assertEqual(sorted(r[:3] for r in bmr.read_as_enricher(out, "snp", {"HGDP"},
                                                                    stats=st)),
                         [(1100, "C", "G"), (1200, "A", "G")])
        self.assertEqual(st["merged"], 0)

    def test_a_record_without_af_max_gets_the_one_its_sources_give(self):
        """The as-run merge declares AF_max and writes it on none of chrY's records (its
        177,968 PASS records; every other chromosome has it on all): with no AF_max the build
        input had no AF, so no indel row of chrY showed a frequency and the frequency bound
        (D6) had none to read. Such a record gets the largest of its sources' frequencies,
        per ALT and as the source wrote it, as AF_max and AF both (found by the census of
        review r1's fix of lane M); one whose sources give none gets neither."""
        rows = [
            "chrY\t2781613\t.\tTCAGC\tT\t1232.18\tPASS\tAF_HGDP=0.007788",
            "chrY\t2781653\t.\tG\tT\t301.73\tPASS\tAF_HGDP=0.001575;AF_AoU=0.001277",
            "chrY\t2781700\t.\tG\tT,GA\t.\tPASS\tAF_HGDP=0.002,.;AF_AoU=0.001,0.01",
            "chrY\t2781716\t.\tT\tC\t728.3\tPASS\tAF_HGDP=.",
        ]
        with gzip.open(self.inp, "wt") as fh:
            fh.write(HEADER + "\n".join(rows) + "\n")
        out = os.path.join(self.d, "build.vcf")
        stats = pm.prepare(self.inp, out)
        body = [l.split("\t") for l in self.read(out) if not l.startswith("#")]
        self.assertEqual([r[7] for r in body], [
            "AF_HGDP=0.007788;AF_max=0.007788;AF=0.007788",
            "AF_HGDP=0.001575;AF_AoU=0.001277;AF_max=0.001575;AF=0.001575",
            "AF_HGDP=0.002,.;AF_AoU=0.001,0.01;AF_max=0.002,0.01;AF=0.002,0.01",
            "AF_HGDP=.",
        ])
        self.assertEqual((stats["af"], stats["af_max_derived"], stats["no_af_max"]), (3, 3, 1))
        # the indel frequencies and the checks read it (verify_mega_chrom derives it alike)
        import build_mega_indel_af as bmi
        recs = list(bmi.iter_indel_af_records(out, {"HGDP", "AoU"}))
        self.assertEqual([(r[0], r[4]) for r in recs], [("2781613", 0.007788), ("2781700", 0.01)])

    def test_a_duplicate_on_a_line_with_samples_stops(self):
        with gzip.open(self.inp, "wt") as fh:
            fh.write(HEADER.replace("\tINFO\n", "\tINFO\tFORMAT\tS1\n")
                     + "chr21\t10\t.\tA\tG,T\t.\tPASS\tAF_max=0.1,0.2\tGT\t1|2\n"
                     + "chr21\t10\t.\tA\tT\t.\tPASS\tAF_max=0.3\tGT\t0|1\n")
        with self.assertRaises(SystemExit) as cm:
            pm.prepare(self.inp, os.path.join(self.d, "build.vcf"))
        self.assertIn("sample columns", str(cm.exception))

    @unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF) not installed")
    def test_command_writes_bgzf(self):
        out = os.path.join(self.d, "build.vcf.gz")
        done = subprocess.run([sys.executable, os.path.join(HERE, "prepare_mega_build_vcf.py"),
                               self.inp, out], capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("6 records -> 4 kept (FILTER PASS or '.'), 2 of another FILTER left "
                      "out, 0 left out as a second record of one REF>ALT", done.stderr)
        with open(out, "rb") as fh:
            self.assertEqual(fh.read(4), b"\x1f\x8b\x08\x04")  # BGZF, for tabix
        self.assertEqual(sum(1 for l in self.read(out) if not l.startswith("#")), 4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
