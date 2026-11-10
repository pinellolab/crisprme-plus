#!/usr/bin/env python3
"""The sites-only density guard (sites_only_guard.py, CRISPRme 2.10.0): of a sites-only
dataset, the raw candidates whose window spells 2^bits sequences or more leave the
targets before the SNP and indel post-analysis read them, and their regions are listed;
a genotyped dataset's targets, and any with the guard off, are left as they are.

Run: cd PostProcess && python -m unittest test_sites_only_guard -v
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import sites_only_guard as sog  # noqa: E402
import window_density as wd  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCACNNN"


def target(dna, pos, chrom="chr22", crrna=GUIDE, strand="+"):
    """A raw CRISPRitz target line: Bulge_type crRNA DNA Chromosome Position
    Cluster_Position Direction Mismatches Bulge_Size Total."""
    return "\t".join(["X", crrna, dna, chrom, str(pos), str(pos), strand, "3", "0", "3"]) + "\n"


# windows of 23 letters spelling 2^20 sequences (ten four-allele letters, Q), 2^24 (twelve),
# 2^9 (nine two-allele letters), and a plain one
DENSE = "QQQQQQQQQQACGTACGTACGTA"          # 4^10 = 2^20: withheld at 20 bits
DENSER = "QQQQQQQQQQQQACGTACGTACG"         # 4^12 = 2^24
SPARSE = "RYRYRYRYRACGTACGTACGTAC"          # 2^9
PLAIN = "CTAACAGTTGCTTTTATCACAGG"


class TestTheDensityTest(unittest.TestCase):
    def test_it_is_sequences_against_the_bound(self):
        for w in (DENSE, DENSER, SPARSE, PLAIN, "BDHV" * 5 + "AGG", "Q-Q-Q" * 4,
                  "ryswkmryswkm" + "ACGTACGTACG"):
            for bits in (0, 1, 10, 19, 20, 21, 24, 40):
                self.assertEqual(wd.guard_withholds(w, bits),
                                 wd.sequences(w) >= (1 << bits), (w, bits))
        self.assertFalse(wd.guard_withholds(DENSER, -1))  # off
        self.assertEqual(wd.SITES_ONLY_GUARD_BITS, wd.HYPERDENSE_BITS)


class _Dicts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="guard_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        self.registry("hg38_mega", {"aggregation": "info_af", "data_type": "sites-only"})
        self.registry("hg38_G", {"aggregation": "carriers", "data_type": "genotyped-phased"})
        os.environ.pop(sog.ENV, None)
        self.addCleanup(os.environ.pop, sog.ENV, None)

    def registry(self, name, manifest):
        d = os.path.join(self.dicts, f"registry_{name}")
        os.makedirs(d)
        with open(os.path.join(d, "reg_chr22.idx"), "w") as fh:
            json.dump(manifest, fh)

    def targets(self, lines, name="t.txt"):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as fh:
            fh.writelines(lines)
        return path

    def read(self, path):
        with open(path) as fh:
            return fh.readlines()


class TestTheSnpSide(_Dicts):
    def test_a_sites_only_datasets_dense_candidates_are_withheld(self):
        lines = ["#Bulge_type\tcrRNA\tDNA\tChromosome\n", target(DENSE, 100),
                 target(PLAIN, 300), target(DENSER, 104, strand="-"), target(SPARSE, 500),
                 target(DENSE, 2000, crrna="ACTGAAATCTGTAAGCAGGCNNN")]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "x.snp" + sog.BED_SUFFIX)
        got = sog.run("snp", path, self.dicts, "hg38_mega", "chr22", bed)
        self.assertEqual(got, (5, 3))
        self.assertEqual(self.read(path), [lines[0], lines[2], lines[4]])
        body = [l.rstrip("\n").split("\t") for l in self.read(bed) if not l.startswith("#")]
        # one region per guide of overlapping windows (100..123 and 104..127 are one)
        self.assertEqual(body, [
            ["chr22", "100", "127", GUIDE, "2", "24.00"],
            ["chr22", "2000", "2023", "ACTGAAATCTGTAAGCAGGCNNN", "1", "20.00"]])
        head = sog.read_bed_counts(bed)
        self.assertEqual(head, {"bits": 20, "dataset": "hg38_mega", "chrom": "chr22",
                                "kind": "snp", "candidates": 5, "withheld": 3})

    def test_a_genotyped_dataset_and_the_guard_off_are_left_alone(self):
        lines = [target(DENSE, 100), target(PLAIN, 300)]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "g.bed")
        self.assertIsNone(sog.run("snp", path, self.dicts, "hg38_G", "chr22", bed))
        self.assertIsNone(sog.run("snp", path, self.dicts, "hg38_unknown", "chr22", bed))
        os.environ[sog.ENV] = "-1"
        self.assertIsNone(sog.run("snp", path, self.dicts, "hg38_mega", "chr22", bed))
        self.assertEqual(self.read(path), lines)
        self.assertFalse(os.path.exists(bed))

    def test_the_bound_is_the_environments(self):
        lines = [target(DENSE, 100), target(SPARSE, 300)]
        path = self.targets(lines)
        os.environ[sog.ENV] = "9"
        got = sog.run("snp", path, self.dicts, "hg38_mega", "chr22",
                      os.path.join(self.tmp, "b.bed"))
        self.assertEqual(got, (2, 2))
        self.assertEqual(self.read(path), [])


class TestTheIndelSide(_Dicts):
    def test_regions_are_placed_on_the_reference_by_the_log(self):
        # two fake-indel contigs: chr22_1000-1053_1 at fake 0..53, chr22_5000-5054_2 at 54..
        log = os.path.join(self.dicts, "log_indels_hg38_mega")
        os.makedirs(log)
        with gzip.open(os.path.join(log, "logchr22.txt.gz"), "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            fh.write("chr22_1000-1053_1\t\t.\t0.1\tchr22_1026_A_AT\t0,53\tN\n")
            fh.write("chr22_5000-5054_2\t\t.\t0.1\tchr22_5026_AC_A\t54,107\tN\n")
        lines = [target(DENSE, 10, chrom="fakechr22"), target(PLAIN, 12, chrom="fakechr22"),
                 target(DENSER, 60, chrom="fakechr22")]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "x.indel" + sog.BED_SUFFIX)
        self.assertEqual(sog.run("indel", path, self.dicts, "hg38_mega", "chr22", bed), (3, 2))
        self.assertEqual(self.read(path), [lines[1]])
        body = [l.rstrip("\n").split("\t")[:3] for l in self.read(bed) if not l.startswith("#")]
        self.assertEqual(body, [["chr22", "1010", "1033"], ["chr22", "5006", "5029"]])

    def test_a_window_it_cannot_place_is_analysed_as_it_is(self):
        # every candidate it withholds is listed: one before the log's first contig (fake 20)
        # has no place on the reference, so the guard leaves it to the post-analysis
        log = os.path.join(self.dicts, "log_indels_hg38_mega")
        os.makedirs(log)
        with gzip.open(os.path.join(log, "logchr22.txt.gz"), "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            fh.write("chr22_1000-1053_1\t\t.\t0.1\tchr22_1026_A_AT\t20,73\tN\n")
        lines = [target(DENSE, 5, chrom="fakechr22"), target(DENSER, 30, chrom="fakechr22")]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "x.indel" + sog.BED_SUFFIX)
        self.assertEqual(sog.run("indel", path, self.dicts, "hg38_mega", "chr22", bed), (2, 1))
        self.assertEqual(self.read(path), [lines[0]])
        body = [l.rstrip("\n").split("\t")[:3] for l in self.read(bed) if not l.startswith("#")]
        self.assertEqual(body, [["chr22", "1010", "1033"]])

    def test_without_the_indel_log_the_indel_side_is_left_unguarded(self):
        lines = [target(DENSE, 10, chrom="fakechr22"), target(PLAIN, 12, chrom="fakechr22")]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "x.indel" + sog.BED_SUFFIX)
        done = subprocess.run([sys.executable, os.path.join(HERE, "sites_only_guard.py"),
                               "indel", path, self.dicts, "hg38_mega", "chr22", bed],
                              capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        self.assertIn("WARNING [sites-only guard]: hg38_mega chr22 indel: not applied", done.stdout)
        self.assertIn("no indel log", done.stdout)
        self.assertEqual(self.read(path), lines)
        self.assertFalse(os.path.exists(bed))


class TestTheWithheldWindowsAreInTheDenseRegionList(_Dicts):
    """The windows the guard withholds go to the dense-region list, the file every row is
    noted by (resultIntegrator's High_variant_density_region) and the report bundles, in
    the columns of the SNP post-analysis's own (new_simple_analysis), their greedy column
    ``withheld``: the SNP side one line a window, as the post-analysis lists a dense one,
    so a row the guard keeps there (a reference row, a genotyped dataset's) keeps its note,
    and the four-allele positions of those windows; the indel side one line a region. Up to
    review r1 of lane M (M-R1-4) they were in the guard's own BED only: on the rebuilt mega's
    chr22 (five guides, NRG 4 + 1 + 1) the dense-region list went from 137,427 lines to
    30,389, 19 rows the guard kept lost their note and 89 got another."""

    HEADER = ("#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\tiupac_protospacer"
              "\tbits\tfour_allele\tgreedy\n")

    def body(self, path):
        return [l.rstrip("\n").split("\t") for l in self.read(path) if not l.startswith("#")]

    def test_the_snp_side_lists_each_window_as_the_post_analysis_does(self):
        # one window read by two raw targets (a mismatch's lower case), one on the minus
        # strand (listed on the forward strand), and a plain target
        lines = [target(DENSE, 100), target(DENSE[:12] + "c" + DENSE[13:], 100),
                 target(PLAIN, 300), target(DENSER, 104, strand="-")]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "x_chr22.snp" + sog.BED_SUFFIX)
        self.assertEqual(sog.run("snp", path, self.dicts, "hg38_mega", "chr22", bed), (4, 3))
        hvdr = os.path.join(self.tmp, "x_chr22.snp.sites_only_density_guard"
                                      ".high_variant_density_regions.bed")
        quads = os.path.join(self.tmp, "x_chr22.snp.sites_only_density_guard"
                                       ".four_allele_positions.bed")
        self.assertEqual(sog.companion_paths(bed), (hvdr, quads))
        self.assertEqual(self.read(hvdr)[0], self.HEADER)
        forward = DENSER.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        self.assertEqual(self.body(hvdr), [
            # of the two lines of one window the one resultIntegrator reads first: the
            # smaller letters (an upper-case base before a lower-case mismatch)
            ["chr22", "100", "123", GUIDE, "10", ".", DENSE, "20.00", "10", "withheld"],
            ["chr22", "104", "127", GUIDE, "12", ".", forward, "24.00", "12", "withheld"]])
        # the four-allele positions of the windows, as the post-analysis lists every one
        # it reads (0-based, the minus-strand window read backwards)
        self.assertEqual(self.body(quads),
                         [["chr22", str(p), str(p + 1)] for p in range(100, 127)
                          if p < 110 or 104 + 23 - 12 <= p])

    def test_the_indel_side_lists_its_regions(self):
        log = os.path.join(self.dicts, "log_indels_hg38_mega")
        os.makedirs(log)
        with gzip.open(os.path.join(log, "logchr22.txt.gz"), "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            fh.write("chr22_1000-1053_1\t\t.\t0.1\tchr22_1026_A_AT\t0,53\tN\n")
        lines = [target(DENSE, 10, chrom="fakechr22"), target(DENSER, 14, chrom="fakechr22"),
                 target(PLAIN, 30, chrom="fakechr22")]
        path = self.targets(lines)
        bed = os.path.join(self.tmp, "x_chr22.indel" + sog.BED_SUFFIX)
        self.assertEqual(sog.run("indel", path, self.dicts, "hg38_mega", "chr22", bed), (3, 2))
        hvdr, quads = sog.companion_paths(bed)
        # one region, 1010..1037, its densest window's letters
        self.assertEqual(self.body(hvdr), [["chr22", "1010", "1037", GUIDE, "12", ".", DENSER,
                                            "24.00", "12", "withheld"]])
        self.assertFalse(os.path.exists(quads))

    def test_a_genotyped_dataset_or_the_guard_off_lists_nothing(self):
        path = self.targets([target(DENSE, 100)])
        bed = os.path.join(self.tmp, "g_chr22.snp" + sog.BED_SUFFIX)
        self.assertIsNone(sog.run("snp", path, self.dicts, "hg38_G", "chr22", bed))
        self.assertFalse(any(os.path.exists(p) for p in sog.companion_paths(bed)))


def _gnu_cut():
    """resultIntegrator drops its unused columns with GNU ``cut --complement``."""
    try:
        res = subprocess.run(["cut", "--complement", "-f1"], input=b"a\tb\n",
                             capture_output=True)
    except OSError:
        return False
    return res.returncode == 0


_INTERVALTREE_STUB = """import collections

Interval = collections.namedtuple("Interval", "begin end data")


class IntervalTree:
    def __init__(self):
        self.items = []

    def __setitem__(self, key, value):
        self.items.append(Interval(key.start, key.stop, value))

    def __getitem__(self, key):
        if isinstance(key, slice):
            return {i for i in self.items if i.begin < key.stop and key.start < i.end}
        return {i for i in self.items if i.begin <= key < i.end}

    def __len__(self):
        return len(self.items)
"""


@unittest.skipUnless(_gnu_cut(), "resultIntegrator.py needs GNU cut --complement")
class TestARowTheGuardKeepsKeepsItsNote(_Dicts):
    """The real resultIntegrator.py on the rows a search writes in a window the guard
    withheld (a reference row of chr22:18,201,798 (-) of the README guide lost its
    "high_variant_density (16 variants, 2^22.51 sequences)" note with the guard on): the
    note names the window, says that the guard withheld a sites-only dataset's alignments
    of it, and its four-allele positions."""

    def test_the_note(self):
        from test_seed_counts_dna_bulge import _block
        with open(os.path.join(self.tmp, "intervaltree.py"), "w") as fh:
            fh.write(_INTERVALTREE_STUB)
        res_dir = os.path.join(self.tmp, "res")
        os.makedirs(res_dir)
        path = self.targets([target(DENSE, 100, chrom="chrT"), target(PLAIN, 3000, "chrT")])
        sog.run("snp", path, self.dicts, "hg38_mega", "chrT",
                os.path.join(res_dir, "job_chrT.snp" + sog.BED_SUFFIX))
        aln = "CTAACAGTTGCTTTTATCACTGG"
        best = os.path.join(res_dir, "job.bestMerge.txt")
        with open(best, "w") as fh:
            for start in (100, 3000):
                block = _block("X", GUIDE, GUIDE, aln, "n", 0, 0, "NA", "NA")
                block[5] = block[6] = str(start)
                fh.write("\t".join(block * 3) + "\n")
        found = os.path.join(res_dir, "job.found.bed")
        with open(found, "w") as fh:
            fh.write("NA\nNA\n")
        empirical = os.path.join(self.tmp, "dummy.txt")
        open(empirical, "w").close()
        guides = os.path.join(self.tmp, "guides.txt")
        with open(guides, "w") as fh:
            fh.write(GUIDE + "\n")
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "resultIntegrator.py"), best, empirical,
             found, guides, res_dir + "/", "true", "hg38", "vuota", "1", "0", "none"],
            cwd=HERE, env=dict(os.environ, PYTHONPATH=self.tmp), capture_output=True,
            text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        import csv
        with open(best + ".integrated_results.tsv", newline="") as fh:
            notes = {int(r["Start_coordinate_(highest_CFD)"]): r["High_variant_density_region"]
                     for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(
            notes[100],
            "high_variant_density (10 variants, 2^20.00 sequences): the sites-only density "
            "guard withheld a sites-only dataset's alignments of this window; other haplotype "
            "alignments may exist; full_IUPAC=%s; four_allele_positions=%s"
            % (DENSE, ",".join("chrT:%d" % p for p in range(101, 111))))
        self.assertEqual(notes[3000], "NA")


class TestTheCommandNeverFailsTheRun(_Dicts):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, os.path.join(HERE, "sites_only_guard.py"),
                               *args], capture_output=True, text=True)

    def test_a_bad_targets_file_is_left_as_it_was(self):
        path = os.path.join(self.tmp, "missing.txt")
        done = self.run_cli("snp", path, self.dicts, "hg38_mega", "chr22",
                            os.path.join(self.tmp, "o.bed"))
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        self.assertIn("WARNING [sites-only guard]", done.stdout)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "o.bed")))

    def test_it_says_what_it_withheld(self):
        path = self.targets([target(DENSE, 100), target(PLAIN, 200)])
        done = self.run_cli("snp", path, self.dicts, "hg38_mega", "chr22",
                            os.path.join(self.tmp, "o.bed"))
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        self.assertIn("Sites-only density guard (hg38_mega, chr22, SNP): 1 of 2 candidate(s) "
                      "in windows of 2^20 sequences or more withheld, in 1 region(s)",
                      done.stdout)


class TestThePostAnalysisRunsIt(unittest.TestCase):
    """post_analisi_snp.sh and post_analisi_indel.sh run the guard on the chromosome's
    variant targets before the analysis reads them, with the dataset's Dictionaries."""

    def step(self, script, kind):
        with open(os.path.join(HERE, script)) as fh:
            text = fh.read()
        call = text.index("python sites_only_guard.py " + kind)
        analysis = text.index("./scriptAnalisiNNN_v3.sh" if kind == "snp"
                              else "./analisi_indels_NNN.sh")
        self.assertLess(text.index('targets_tsv_alt_chrom"\n', 0), call)
        self.assertLess(call, analysis)
        return text[call:text.index("\n\n", call)]

    def test_both_sides(self):
        snp = self.step("post_analisi_snp.sh", "snp")
        self.assertIn('"$targets_tsv_alt_chrom" "$(dirname "${dict_folder%/}")"', snp)
        self.assertIn(".snp.sites_only_density_guard.bed", snp)
        indel = self.step("post_analisi_indel.sh", "indel")
        self.assertIn('"$targets_tsv_alt_chrom" "$dict_folder" "$vcf_name" "$chrom"', indel)
        self.assertIn(".indel.sites_only_density_guard.bed", indel)

    def test_an_indel_subset_the_guard_empties_is_not_taken_for_issue_172(self):
        # the warning of an empty per-chromosome subset (no target matched the fake
        # chromosome) is said before the guard runs: a subset it withholds whole is no such
        with open(os.path.join(HERE, "post_analisi_indel.sh")) as fh:
            text = fh.read()
        self.assertLess(text.index("is EMPTY (see GitHub issue #172"),
                        text.index("python sites_only_guard.py indel"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
