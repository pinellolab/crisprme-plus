#!/usr/bin/env python
"""The carrier gate (haplotype_prioritizer.py): the tiers of a row's variants against phased
stores the repository's own writers build (tier1_genotypes, build_indel_genotypes), the
panels a search reads, and the command the search pipeline runs, on an integrated results
table it must leave as it is.

Run: cd PostProcess && python -m unittest test_haplotype_prioritizer -v
"""

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import haplotype_prioritizer as hp  # noqa: E402
import tier1_genotypes as t1  # noqa: E402
import build_indel_genotypes as big  # noqa: E402
import test_build_index_only_atomic as bio  # noqa: E402  (the real crisprme.py, a stub engine)

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule

CHROM = "chr22"
PANEL = "hg38_TESTPANEL"
# 20 samples: 40 haplotypes on an autosome
SAMPLES = [f"S{i:02d}" for i in range(20)]
SAMPLE_META = {s: ("TEST", "AFR" if i < 12 else "EUR", "female") for i, s in enumerate(SAMPLES)}


def build_panel(dicts, snps, indels=(), name=PANEL):
    """A phased store of ``name`` under ``dicts``: ``snps`` [(pos, REF, ALT, {sample: gt})]
    through tier1_genotypes, ``indels`` [(pos, REF, ALT, {sample: gt})] through a VCF and
    build_indel_genotypes, as build-index-only writes them."""
    gdir = os.path.join(dicts, f"genotypes_{name}")
    os.makedirs(gdir, exist_ok=True)
    t1.compile_genotypes([(p, a, dict(g), r) for p, r, a, g in snps],
                         t1.build_sample_axis(SAMPLE_META),
                         os.path.join(gdir, f"gt_{CHROM}.bin"),
                         os.path.join(gdir, f"gt_{CHROM}.idx"))
    if indels:
        idir = os.path.join(dicts, f"indel_genotypes_{name}")
        os.makedirs(idir, exist_ok=True)
        vcf = os.path.join(dicts, f"{name}.indels.vcf")
        with open(vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(SAMPLES) + "\n")
            for p, r, a, gts in indels:
                fh.write(f"{CHROM}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t"
                         + "\t".join(gts.get(s, "0|0") for s in SAMPLES) + "\n")
        big.compile_indel_genotypes(vcf, os.path.join(idir, f"gt_indel_{CHROM}.tsv.gz"))
    return gdir


def v(pos, ref, alt):
    return (CHROM, pos, ref, alt)


class TestParseAndSupport(unittest.TestCase):
    def test_parse_variants(self):
        self.assertEqual(hp.parse_variants("NA"), [])
        self.assertEqual(hp.parse_variants("n"), [])
        self.assertEqual(
            hp.parse_variants("chr1_KI270706v1_random_100_A_G,chr2_5_AT_A,chr2_5_AT_A"),
            [("chr1_KI270706v1_random", 100, "A", "G"), ("chr2", 5, "AT", "A")])

    def test_cis_trans_unphased_haploid(self):
        a = {"S1": "1|0", "S2": "1|0", "S3": "1/0", "S4": "1", "S5": "1|1"}
        b = {"S1": "1|0", "S2": "0|1", "S3": "0/1", "S4": "1", "S5": "0|1"}
        haps, cis, unknown, carrying = hp.support([a, b])
        # S1: cis on slot 0; S2: trans; S3: unphased; S4: haploid, both on its one slot;
        # S5: the second slot holds both
        self.assertEqual((haps, cis, unknown, carrying), (3, {"S1", "S4", "S5"}, {"S3"}, 3))

    def test_breaker_on_the_slot_or_the_other(self):
        snp = {"S1": "1|0", "S2": "1|0", "S3": "1|0"}
        indel = {"S1": "1|0", "S2": "0|1", "S3": "1/0"}
        haps, cis, unknown, _c = hp.support([snp], breakers=[indel])
        # S1 loses its slot to the indel; S2's indel is on the other slot; S3's indel call
        # is unphased, so its slot may or may not hold the target
        self.assertEqual((haps, cis, unknown), (1, {"S2"}, {"S3"}))

    def test_two_alts_on_one_slot_is_unknown(self):
        alt_g = {"S1": "1|0", "S2": "1|0"}
        alt_t = {"S1": "1|0"}  # S1's first slot holds G and T at one position
        haps, cis, unknown, _c = hp.support([alt_g], rivals=[alt_t])
        self.assertEqual((haps, cis, unknown), (1, {"S2"}, {"S1"}))


class TestClassify(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def panel(self, snps, indels=()):
        build_panel(self.dicts, snps, indels)
        p = hp.PanelChrom(PANEL, self.dicts, CHROM)
        self.addCleanup(p.close)
        return p

    def test_a_two_variants_on_one_haplotype(self):
        p = self.panel([
            (100, "A", "G", {"S00": "1|0", "S01": "1|0", "S02": "0|1"}),
            (105, "C", "T", {"S00": "1|0", "S01": "1|0", "S03": "1|0"}),
        ])
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p], span=(95, 117))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("A", "", 2))
        self.assertEqual(ver.samples, ("S00", "S01"))
        self.assertEqual(ver.panel_haplotypes, 40)
        self.assertEqual(ver.label(), f"together on 2 of 40 {PANEL} haplotypes")

    def test_a_one_haplotype(self):
        p = self.panel([
            (100, "A", "G", {"S00": "1|0", "S02": "0|1"}),
            (105, "C", "T", {"S00": "1|0", "S02": "1|0"}),
        ])
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p])
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("A", hp.ONE_HAPLOTYPE, 1))

    def test_c1_never_together_with_and_without_power(self):
        # common variants, always in trans: 40 x 0.5 x 0.5 = 10 haplotypes expected
        common_a = {s: "1|0" for s in SAMPLES}
        common_b = {s: "0|1" for s in SAMPLES}
        p = self.panel([(100, "A", "G", common_a), (105, "C", "T", common_b),
                        (110, "G", "A", {"S00": "1|0"}), (112, "T", "C", {"S01": "1|0"})])
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p])
        self.assertEqual((ver.tier, ver.state), ("C1", ""))
        self.assertAlmostEqual(ver.expected, 10.0)
        self.assertEqual(ver.label(), f"not together on any of 40 {PANEL} haplotypes "
                                      "(expected 10)")
        # two singletons: 40 x 1/40 x 1/40 = 0.025 expected, a zero says nothing
        rare = hp.classify([v(110, "G", "A"), v(112, "T", "C")], [p])
        self.assertEqual((rare.tier, rare.state), ("C1", hp.TOO_RARE))
        self.assertIn("too rare to tell (expected 0.025)", rare.label())
        # the threshold is a parameter: 0 turns the state off
        off = hp.classify([v(110, "G", "A"), v(112, "T", "C")], [p], min_expected=0)
        self.assertEqual((off.tier, off.state), ("C1", ""))

    def test_c2_and_b(self):
        p = self.panel([(100, "A", "G", {"S00": "1|0"})])
        c2 = hp.classify([v(100, "A", "G"), v(300, "C", "T")], [p],
                         af_of=lambda var: {100: 0.02, 300: 0.0001}[var[1]])
        self.assertEqual((c2.tier, c2.missing, c2.n_variants), ("C2", 1, 2))
        self.assertEqual(c2.label(), f"not checked: 1 of 2 variants not in {PANEL}; frequency "
                                     "product 2.0e-06, rarest 1.0e-04 (below the 1e-06 / "
                                     "0.001 bound)")
        self.assertEqual(hp.classify([v(100, "A", "G")], [p]).tier, "B")
        absent = hp.classify([v(300, "C", "T")], [p])
        self.assertEqual((absent.tier, absent.state), ("B", hp.NOT_IN_PANEL))
        # a REF the panel does not have is not that variant
        self.assertEqual(hp.classify([v(100, "C", "G"), v(300, "C", "T")], [p]).missing, 2)

    def test_indel_in_the_target_takes_the_slot(self):
        # the SNP at 100 on slot 0 of S00 and S01; a phased deletion at 104 on the same
        # slots of both, inside a target over 90..112: no haplotype has the row's sequence
        snps = [(100, "A", "G", {"S00": "1|0", "S01": "1|0"})]
        indels = [(104, "CTT", "C", {"S00": "1|0", "S01": "1|0", "S05": "0|1"})]
        p = self.panel(snps, indels)
        ver = hp.classify([v(100, "A", "G")], [p], span=(90, 112))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("C1", hp.INDEL_IN_TARGET, 2))
        self.assertEqual(ver.label(), f"all 2 {PANEL} haplotypes carrying it have an indel "
                                      "in the target")
        # a target that ends before the deletion's bases keeps them
        self.assertEqual(hp.classify([v(100, "A", "G")], [p], span=(90, 103)).tier, "B")

    def test_two_snps_together_only_beside_an_indel_in_the_target(self):
        snps = [(100, "A", "G", {"S00": "1|0", "S01": "1|0"}),
                (102, "C", "T", {"S00": "1|0", "S01": "1|0"})]
        indels = [(104, "CTT", "C", {"S00": "1|0", "S01": "1|0"})]
        p = self.panel(snps, indels)
        ver = hp.classify([v(100, "A", "G"), v(102, "C", "T")], [p], span=(90, 112))
        self.assertEqual((ver.tier, ver.state, ver.haplotypes), ("C1", hp.INDEL_IN_TARGET, 2))
        self.assertEqual(ver.label(), f"all 2 {PANEL} haplotypes carrying them together have "
                                      "an indel in the target")

    def test_indel_on_the_other_slot_unphased_or_at_the_snp(self):
        snps = [(100, "A", "G", {"S00": "1|0", "S01": "1|0", "S02": "1|0"})]
        indels = [
            (104, "CTT", "C", {"S00": "0|1", "S01": "1/0"}),  # other slot; unphased
            (100, "AT", "A", {"S02": "1|0"}),  # at the SNP's own position: an artefact
        ]
        p = self.panel(snps, indels)
        self.assertEqual(hp.classify([v(100, "A", "G")], [p], span=(90, 112)).tier, "B")
        only_s02 = self.panel([(100, "A", "G", {"S02": "1|0"})], indels)
        self.assertEqual(hp.classify([v(100, "A", "G")], [only_s02], span=(90, 112)).tier, "B")

    def test_multi_alt_genotypes_are_read_per_alt(self):
        # one VCF line with two ALTs: S00 is G on slot 0 and T on slot 1 (1|2), listed under
        # both; with 105 C>T on slot 1, T and 105 are together, G and 105 are not
        p = self.panel([
            (100, "A", "G", {"S00": "1|2", "S01": "1|0"}),
            (100, "A", "T", {"S00": "1|2", "S02": "0|2"}),
            (105, "C", "T", {"S00": "0|1", "S02": "0|1"}),
        ])
        t = hp.classify([v(100, "A", "T"), v(105, "C", "T")], [p])
        self.assertEqual((t.tier, t.haplotypes, t.samples), ("A", 2, ("S00", "S02")))
        g = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p])
        self.assertEqual(g.tier, "C1")

    def test_two_panels_the_one_that_holds_them_decides(self):
        build_panel(self.dicts, [(100, "A", "G", {"S00": "1|0"}),
                                 (105, "C", "T", {"S00": "1|0"})], name="hg38_P1")
        build_panel(self.dicts, [(100, "A", "G", {"S03": "1|0"})], name="hg38_P2")
        p1 = hp.PanelChrom("hg38_P1", self.dicts, CHROM)
        p2 = hp.PanelChrom("hg38_P2", self.dicts, CHROM)
        self.addCleanup(p1.close)
        self.addCleanup(p2.close)
        ver = hp.classify([v(100, "A", "G"), v(105, "C", "T")], [p2, p1])
        self.assertEqual((ver.tier, ver.panel), ("A", "hg38_P1"))


def write_results(path, rows):
    """An integrated results table with the columns the gate reads (and a filler)."""
    cols = ["Spacer+PAM", "Chromosome"]
    for block in hp.BLOCKS:
        cols += [f"{c}_({block})" for c in hp._BLOCK_COLUMNS]
    cols.append("Annotation_GENCODE")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "NA") for c in cols) + "\n")


def block(name, start, strand, aligned, origin, variants, maf="NA"):
    return {f"Start_coordinate_({name})": str(start), f"Strand_({name})": strand,
            f"Aligned_protospacer+PAM_ALT_({name})": aligned,
            f"REF/ALT_origin_({name})": origin,
            f"Variant_info_genome_({name})": variants, f"Variant_MAF_({name})": maf}


class TestGateCommand(unittest.TestCase):
    """The command the pipeline runs: the companion, the sidecar, and an integrated results
    table left byte for byte as it was."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_cli_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        build_panel(self.dicts, [
            (100, "A", "G", {"S00": "1|0", "S01": "1|0"}),
            (105, "C", "T", {"S00": "1|0", "S02": "0|1"}),
            (110, "G", "A", {"S03": "0|1"}),
        ])
        self.out = os.path.join(self.tmp, "Results", "job")
        os.makedirs(self.out)
        self.vcfs = os.path.join(self.tmp, "list_vcf.txt")
        with open(self.vcfs, "w") as fh:
            fh.write(f"VCFs/{PANEL}\n")
        self.results = os.path.join(
            self.out, "GUIDE+NRG_hg38+" + PANEL + "_4+2_integrated_results.tsv")
        aln = "ACGTACGTACGTACGTACGTAGG"
        two = f"{CHROM}_100_A_G,{CHROM}_105_C_T"
        trans = f"{CHROM}_105_C_T,{CHROM}_110_G_A"
        rows = [
            # 1: reference row
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 90, "+", aln, "ref", "NA"),
             **block("fewest_mm+b", 90, "+", aln, "ref", "NA"),
             **block("highest_CRISPR_BULGE", 90, "+", aln, "ref", "NA")},
            # 2: CFD alignment A, CRISPR-Bulge alignment C1
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 90, "+", aln, "alt", two, "0.05,0.05"),
             **block("fewest_mm+b", 90, "+", aln, "alt", f"{CHROM}_100_A_G", "0.05"),
             **block("highest_CRISPR_BULGE", 90, "-", aln, "alt", trans, "0.05,0.025")},
            # 3: a variant no panel holds
            {"Spacer+PAM": "G", "Chromosome": CHROM,
             **block("highest_CFD", 200, "+", aln, "alt", f"{CHROM}_100_A_G,{CHROM}_210_T_C"),
             **block("fewest_mm+b", 200, "+", aln, "ref", "NA"),
             **block("highest_CRISPR_BULGE", 200, "+", aln, "ref", "NA")},
        ]
        write_results(self.results, rows)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_gate(self, *extra):
        cmd = [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"),
               "--results", self.results, "--dictionaries", self.dicts,
               "--vcf-list", self.vcfs, *extra]
        return subprocess.run(cmd, capture_output=True, text=True)

    def md5(self):
        with open(self.results, "rb") as fh:
            return hashlib.md5(fh.read()).hexdigest()

    def test_population_level_writes_the_companion_and_leaves_the_rows(self):
        before = self.md5()
        done = self.run_gate("--mode", "population-level")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(done.stderr, "")
        self.assertIn("Carrier gate: the rows' variant combinations checked", done.stdout)
        self.assertEqual(self.md5(), before)  # the off-target set is the same file
        comp = hp.companion_path(self.results)
        with open(comp) as fh:
            lines = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("# ")]
        header, body = lines[0], lines[1:]
        self.assertEqual(header[0], "#Row")
        got = [(r[0], r[1], r[7], r[8]) for r in body]
        self.assertEqual(got, [
            ("2", "highest_CFD", "A", hp.ONE_HAPLOTYPE),
            ("2", "fewest_mm+b", "B", "-"),
            ("2", "highest_CRISPR_BULGE", "C1", hp.TOO_RARE),
            ("3", "highest_CFD", "C2", "-"),
        ])
        with open(os.path.join(self.out, hp.META_NAME)) as fh:
            meta = json.load(fh)
        self.assertEqual(meta["status"], "run")
        self.assertEqual(meta["counts"], {"A": 1, "C2": 1})
        self.assertEqual(meta["panels"][0]["name"], PANEL)
        self.assertEqual(meta["panels"][0]["samples"], 20)
        self.assertEqual(meta["panels"][0]["superpopulations"], {"AFR": 12, "EUR": 8})
        self.assertEqual(meta["companion"], os.path.basename(comp))

    def test_off_and_per_sample_write_no_companion(self):
        for extra, status in ((("--mode", "population-level", "--prioritize-observed", "off"),
                               "off"), (("--mode", "per-sample"), "per-sample")):
            done = self.run_gate(*extra)
            self.assertEqual((done.returncode, done.stderr), (0, ""))
            with open(os.path.join(self.out, hp.META_NAME)) as fh:
                self.assertEqual(json.load(fh)["status"], status)
            self.assertFalse(os.path.exists(hp.companion_path(self.results)))


class TestReportRanksByTheGate(unittest.TestCase):
    """generate_report on a population-level run the gate checked: the C1 row leaves the
    top-1000, the panel and the risk card for its own table, the A and C2 rows keep their
    place with their tier as their Observed cell, and the rows themselves are the search's."""

    GUIDE = "CTCTCAGCTGGTACACGGCANNN"

    def setUp(self):
        import generate_report  # noqa: F401  (pandas + matplotlib: installed in CI)
        self.tmp = tempfile.mkdtemp(prefix="gate_report_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        build_panel(self.dicts, [
            (1010, "A", "G", {"S00": "1|0", "S01": "1|0"}),
            (1015, "C", "T", {"S00": "1|0", "S01": "1|0", "S02": "0|1"}),
            (2010, "G", "A", {s: "1|0" for s in SAMPLES}),
            (2015, "T", "C", {s: "0|1" for s in SAMPLES}),
            (3010, "A", "C", {"S04": "1|0"}),
            (4010, "C", "G", {"S05": "0|1"}),
        ])
        self.rd = os.path.join(self.tmp, "Results", "job")
        os.makedirs(self.rd)
        with open(os.path.join(self.rd, ".search_mode"), "w") as fh:
            fh.write("population-level")
        with open(os.path.join(self.rd, ".Params.txt"), "w") as fh:
            fh.write("Mismatches\t4\nDNA\t1\nRNA\t1\nMax_bulges\t2\nPam\tNGG\n")
        self.vcfs = os.path.join(self.tmp, "list_vcf.txt")
        with open(self.vcfs, "w") as fh:
            fh.write(f"{PANEL}\n")
        self.results = os.path.join(
            self.rd, f"{self.GUIDE}+NGG_hg38+{PANEL}_4+2_integrated_results.tsv")
        aln = "CTCTCAGCTGGTACACGGCATGG"
        rows = [
            # chrom pos  mm  cfd   cb     origin variants (CFD and CRISPR-Bulge alignment)
            ("1000", "0", "1.0", "1.0", "ref", "NA"),
            ("500", "3", "0.30", "0.10", "ref", "NA"),
            ("1000", "3", "0.50", "0.40", "alt", f"{CHROM}_1010_A_G,{CHROM}_1015_C_T"),  # A
            ("2000", "3", "0.95", "0.90", "alt", f"{CHROM}_2010_G_A,{CHROM}_2015_T_C"),  # C1
            ("3000", "3", "0.40", "0.30", "alt", f"{CHROM}_3010_A_C,{CHROM}_3020_G_T"),  # C2
            ("4000", "3", "0.20", "0.20", "alt", f"{CHROM}_4010_C_G"),                   # B
        ]
        self.cols = ["Spacer+PAM", "Chromosome", "Not_found_in_REF", "Annotation_GENCODE"]
        per_block = ["Start_coordinate", "Strand", "Aligned_spacer+PAM",
                     "Aligned_protospacer+PAM_REF", "Aligned_protospacer+PAM_ALT", "PAM",
                     "Mismatches", "Bulges", "Mismatches+bulges", "REF/ALT_origin",
                     "PAM_creation", "Variant_info_genome", "Variant_MAF", "Variant_rsID",
                     "Variant_samples"]
        for b in hp.BLOCKS:
            self.cols += [f"{c}_({b})" for c in per_block]
        self.cols += ["CFD_score_(highest_CFD)", "CFD_score_REF_(highest_CFD)",
                      "CFD_score_ALT_(highest_CFD)",
                      "CRISPR_BULGE_score_(highest_CRISPR_BULGE)",
                      "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)",
                      "CRISPR_BULGE_score_ALT_(highest_CRISPR_BULGE)"]
        with open(self.results, "w") as fh:
            fh.write("\t".join(self.cols) + "\n")
            for pos, mm, cfd, cb, origin, variants in rows:
                r = {"Spacer+PAM": self.GUIDE, "Chromosome": CHROM,
                     "Not_found_in_REF": "y" if origin == "alt" else "n",
                     "Annotation_GENCODE": "intergenic"}
                for b in hp.BLOCKS:
                    r.update({f"Start_coordinate_({b})": pos, f"Strand_({b})": "+",
                              f"Aligned_spacer+PAM_({b})": self.GUIDE,
                              f"Aligned_protospacer+PAM_REF_({b})": aln,
                              f"Aligned_protospacer+PAM_ALT_({b})": aln if origin == "alt" else "NA",
                              f"PAM_({b})": "TGG", f"Mismatches_({b})": mm, f"Bulges_({b})": "0",
                              f"Mismatches+bulges_({b})": mm, f"REF/ALT_origin_({b})": origin,
                              f"PAM_creation_({b})": "NA", f"Variant_info_genome_({b})": variants,
                              f"Variant_MAF_({b})": "0.05" if origin == "alt" else "NA",
                              f"Variant_rsID_({b})": ".", f"Variant_samples_({b})":
                                  "S00,S01,S02" if origin == "alt" else "NA"})
                r.update({"CFD_score_(highest_CFD)": cfd, "CFD_score_REF_(highest_CFD)": "0.1",
                          "CFD_score_ALT_(highest_CFD)": cfd,
                          "CRISPR_BULGE_score_(highest_CRISPR_BULGE)": cb,
                          "CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)": "0.05",
                          "CRISPR_BULGE_score_ALT_(highest_CRISPR_BULGE)": cb})
                fh.write("\t".join(r.get(c, "NA") for c in self.cols) + "\n")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def report(self):
        import generate_report as gr
        import zipfile
        out = gr.build_report(result_dir=self.rd, integrated_tsv=self.results,
                              report_full=True)
        files = {}
        with zipfile.ZipFile(out) as zf:
            for name in zf.namelist():
                data = zf.read(name)
                files[name] = gzip.decompress(data) if name.endswith(".gz") else data
        return {k: v.decode() for k, v in files.items()}

    def positions(self, tsv):
        lines = tsv.splitlines()
        col = lines[0].split("\t").index("Position")
        return [ln.split("\t")[col] for ln in lines[1:]]

    def test_c1_leaves_the_rankings_for_its_own_table(self):
        with open(self.results, "rb") as fh:
            source = fh.read()
        # unchecked (no gate), the rows rank by CFD alone
        files = self.report()
        self.assertEqual(self.positions(files["data/top1000.tsv"]),
                         ["2000", "1000", "3000", "500", "4000"])
        self.assertNotIn("data/unobserved_worst_case.tsv", files)
        self.assertNotIn("Carrier gate", files["report.html"])
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "haplotype_prioritizer.py"), "--results",
             self.results, "--dictionaries", self.dicts, "--vcf-list", self.vcfs,
             "--mode", "population-level"], capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        files = self.report()
        top = self.positions(files["data/top1000.tsv"])
        self.assertEqual(top, ["1000", "3000", "500", "4000"])  # by CFD, without 2000
        self.assertNotIn("2000", self.positions(files["data/panel_top100.tsv"]))
        self.assertNotIn("2000", self.positions(files["data/top1000_crispr_bulge.tsv"]))
        self.assertEqual(self.positions(files["data/unobserved_worst_case.tsv"]), ["2000"])
        html = files["report.html"]
        self.assertIn("Carrier gate", html)
        self.assertIn("variant combinations no haplotype of the panel carries", html)
        self.assertIn(f"together on 2 of 40 {PANEL} haplotypes", html)
        self.assertIn(f"not together on any of 40 {PANEL} haplotypes (expected 10)", html)
        self.assertIn(f"not checked: 1 of 2 variants not in {PANEL}", html)
        self.assertIn("AFR 12, EUR 8", html)
        self.assertNotIn("chimera", html.split("Carrier gate", 1)[1].split("</tr>", 1)[0])
        manifest = json.loads(files["data/run_manifest.json"])
        self.assertEqual(manifest["carrier_gate"]["tiers"], {"A": 1, "B": 1, "C1": 1, "C2": 1})
        # the off-target set is the search's: the bundled table is the source file
        self.assertEqual(files["data/integrated_results.tsv.gz"].encode(), source)
        self.assertIn("data/carrier_gate.tsv.gz", files)


class _PostAnalysis(unittest.TestCase):
    """The real new_simple_analysis.py on one variant window of chrT, with the stores the
    repository's writers build: the Tier-0 registry + Tier-1 genotypes (the dict-less
    tiers), or a classic per-sample dictionary, and the phased indel genotypes."""

    GUIDE = "CTAACAGTTGCTTTTATCAC"
    FLANK = "ACGTGCATTC" * 10

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_nsa_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.at = len(self.FLANK)

    def world(self, site, snps, indels, tiers):
        """``snps`` [(column, REF, ALT, {sample: gt})] and ``indels`` [(column, REF, ALT,
        {sample: gt})], columns of ``site``; ``tiers``: "tiers" or "dict"."""
        import tier0_registry as t0
        chrom = "chrT"
        self.fasta = os.path.join(self.tmp, "chrT.fa")
        with open(self.fasta, "w") as fh:
            fh.write(">chrT\n" + self.FLANK + site + self.FLANK + "\n")
        d = os.path.join(self.tmp, "Dictionaries")
        recs = [(self.at + c + 1, r, a, f"rs{c}", dict(g)) for c, r, a, g in snps]
        if tiers == "tiers":
            rd = os.path.join(d, "registry_ds")
            gd = os.path.join(d, "genotypes_ds")
            os.makedirs(rd)
            os.makedirs(gd)
            t0.compile_registry_panel(recs, SAMPLE_META, None, t0.autosomal_ploidy,
                                      os.path.join(rd, f"reg_{chrom}.bin"),
                                      os.path.join(rd, f"reg_{chrom}.idx"))
            t1.compile_genotypes([(p, a, g, r) for p, r, a, _rs, g in recs],
                                 t1.build_sample_axis(SAMPLE_META),
                                 os.path.join(gd, f"gt_{chrom}.bin"),
                                 os.path.join(gd, f"gt_{chrom}.idx"))
        else:
            dd = os.path.join(d, "dictionaries_ds")
            os.makedirs(dd)
            entries = {}
            for p, r, a, rs, g in recs:
                entry = "%s;%s,%s;%s;%s" % (
                    ",".join("%s:%s" % kv for kv in sorted(g.items())), r, a, rs, "0.05")
                key = f"{chrom},{p}"
                entries[key] = entries[key] + "$" + entry if key in entries else entry
            with open(os.path.join(dd, f"my_dict_{chrom}.json"), "w") as fh:
                json.dump(entries, fh)
        if indels:
            igd = os.path.join(d, "indel_genotypes_ds")
            os.makedirs(igd)
            with gzip.open(os.path.join(igd, f"gt_indel_{chrom}.tsv.gz"), "wt") as fh:
                for c, r, a, g in indels:
                    fh.write("%d_%s_%s\t%s\n" % (self.at + c + 1, r, a, ",".join(
                        "%s:%s" % kv for kv in sorted(g.items()))))
        self.dict_path = os.path.join(d, "dictionaries_ds", f"my_dict_{chrom}.json")
        self.pam = os.path.join(self.tmp, "pam.txt")
        with open(self.pam, "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")

    def analyse(self, window, mode, budget, extra_env=None):
        """The rows of new_simple_analysis.py for the variant search's row ``window`` (the
        site with an IUPAC code at each SNP column): [(DNA, mismatches, samples, SNPs)]."""
        targets = os.path.join(self.tmp, f"targets_{mode}.txt")
        row = ["X", self.GUIDE + "NNN", window, "chrT", str(self.at), str(self.at), "+",
               "0", "0", "0", "NA", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
        with open(targets, "w") as fh:
            fh.write("#header\n" + "\t".join(row) + "\n")
        out = os.path.join(self.tmp, f"out_{mode}")
        env = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
        env.update(CRISPRME_FAST_MODE=mode, CRISPRME_SKIP_SCORER_ENV="1", **(extra_env or {}))
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "new_simple_analysis.py"), self.fasta,
             self.dict_path, targets, self.pam, out, str(budget)],
            cwd=HERE, env=env, capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""), done.stdout)
        with open(out + ".bestCFD.txt") as fh:
            lines = [ln.rstrip("\n").split("\t") for ln in fh][1:]
        return sorted((f[2].upper(), int(f[7]), f[12], f[17]) for f in lines)


class TestPerSampleIndelInTheTarget(_PostAnalysis):
    """F2: a per-sample SNP row named a carrier whose haplotype holding the SNP also holds
    an indel inside the target (its sequence is not the row's). S00 holds the SNP and a
    deletion on its first slot; S02 has the deletion on the other slot; S03's deletion call
    is unphased; S01 has no indel."""

    # the guide with a mismatch at column 4 (A for C), which the SNP restores
    SITE = "CTAAAAGTTGCTTTTATCACTGG"
    SNPS = [(4, "A", "C", {"S00": "1|0", "S01": "1|0", "S02": "0|1", "S03": "1|0"})]
    # a deletion of the T at column 13, inside the target
    INDELS = [(12, "TT", "T", {"S00": "1|0", "S02": "1|0", "S03": "1/0"})]
    WINDOW = SITE[:4] + "M" + SITE[5:]

    def check(self, tiers):
        self.world(self.SITE, self.SNPS, self.INDELS, tiers)
        per_sample = self.analyse(self.WINDOW, "0", 0)
        # (the classic dictionary's path gives a row per haplotype slot, which the merge
        # step joins)
        self.assertEqual({r[:2] + (r[3],) for r in per_sample},
                         {(self.GUIDE + "TGG", 0, f"chrT_{self.at + 5}_A_C")})
        self.assertEqual(sorted(s for r in per_sample for s in r[2].split(",")),
                         ["S01", "S02", "S03"])
        # a population-level row lists every carrier of its alleles; the carrier gate
        # reads the indel genotypes for it
        population = self.analyse(self.WINDOW, "1", 0)
        self.assertEqual(sorted(s for r in population for s in r[2].split(",")),
                         ["S00", "S01", "S02", "S03"])

    def test_dict_less_tiers(self):
        self.check("tiers")

    def test_classic_dictionary(self):
        self.check("dict")


class TestPerSampleDenseWindowOnTheDictionary(_PostAnalysis):
    """A per-sample window denser than CRISPRME_IUPAC_CAP on an index without the dict-less
    tiers was reported by its greedy worst-case alignment, which took every allele that
    matches the guide and listed the carriers of any of them: here all three SNPs, 0
    mismatches, for S00, S01 and S02, none of whom carries more than two. Now its rows are
    the haplotypes its carriers have, as on the tiers."""

    # three mismatches with the guide (columns 2, 8 and 15), each restored by a SNP
    SNPS = [(2, "C", "A", {"S00": "1|0"}), (8, "A", "T", {"S00": "1|0", "S02": "0|1"}),
            (15, "T", "A", {"S01": "1|0"})]

    def test_rows_are_the_carriers_haplotypes(self):
        site = list(self.GUIDE + "TGG")
        for col, ref, _alt, _g in self.SNPS:
            site[col] = ref
        site = "".join(site)
        window = list(site)
        for col, ref, alt, _g in self.SNPS:
            window[col] = {frozenset("AC"): "M", frozenset("AT"): "W"}[frozenset(ref + alt)]
        self.world(site, self.SNPS, (), "dict")
        rows = self.analyse("".join(window), "0", 3, {"CRISPRME_IUPAC_CAP": "2"})
        snps = {r[3]: (r[1], r[2]) for r in rows if r[3] not in ("n", "NA")}
        p = lambda c: f"chrT_{self.at + c + 1}"  # noqa: E731
        self.assertEqual(snps, {
            f"{p(2)}_C_A,{p(8)}_A_T": (1, "S00"),
            f"{p(8)}_A_T": (2, "S02"),
            f"{p(15)}_T_A": (2, "S01"),
        })
        self.assertNotIn(0, [r[1] for r in rows])  # no row nobody carries


class TestGenotypeOracleDownload(unittest.TestCase):
    """download --what genotype-oracle: the phased HGDP + 1000 Genomes genotypes and indel
    genotypes into Dictionaries/, where the gate finds them, from archives publish-index
    writes (the HuggingFace client stubbed: the repository is a folder)."""

    def setUp(self):
        import crisprme_hf as hf
        self.hf = hf
        self.tmp = tempfile.mkdtemp(prefix="gate_oracle_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.hub = os.path.join(self.tmp, "hub")
        os.makedirs(self.hub)
        src = os.path.join(self.tmp, "src", "Dictionaries")
        build_panel(src, [(100, "A", "G", {"S00": "1|0"})],
                    [(104, "CTT", "C", {"S00": "1|0"})], name=hp.ORACLE_PANEL)
        self.src = src
        orig = hf._hf_snapshot
        self.addCleanup(setattr, hf, "_hf_snapshot", orig)
        self.fetched = []

        def snapshot(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
            for pat in allow_patterns:
                self.fetched.append(pat)
                f = os.path.join(self.hub, os.path.basename(pat))
                if os.path.isfile(f):
                    os.makedirs(os.path.join(local_dir, os.path.dirname(pat)), exist_ok=True)
                    shutil.copy(f, os.path.join(local_dir, pat))
            return local_dir
        hf._hf_snapshot = snapshot
        self.work = os.path.join(self.tmp, "work")
        os.makedirs(self.work)

    def publish(self, store):
        self.hf._make_genotypes_tarball(
            os.path.join(self.hub, f"{store}_{hp.ORACLE_PANEL}.tar.gz"),
            os.path.join(self.src, f"{store}_{hp.ORACLE_PANEL}"))

    def test_installs_both_stores_where_the_gate_reads_them(self):
        self.publish("genotypes")
        self.publish("indel_genotypes")
        dests = self.hf.download_genotype_oracle(self.work, repo="local/test")
        dicts = os.path.join(self.work, "Dictionaries")
        self.assertEqual(dests, [os.path.join(dicts, f"{s}_{hp.ORACLE_PANEL}")
                                 for s in ("genotypes", "indel_genotypes")])
        self.assertEqual(hp.find_oracle(dicts), (hp.ORACLE_PANEL, dicts))
        p = hp.PanelChrom(hp.ORACLE_PANEL, dicts, CHROM, "oracle")
        self.addCleanup(p.close)
        self.assertEqual(p.carriers(v(100, "A", "G")), {"S00": "1|0"})
        self.assertEqual(p.carriers(v(104, "CTT", "C")), {"S00": "1|0"})
        # installed: kept, not fetched again
        self.fetched.clear()
        self.hf.download_genotype_oracle(self.work, repo="local/test")
        self.assertEqual(self.fetched, [])

    def test_without_indel_genotypes_on_the_hub_it_says_so(self):
        self.publish("genotypes")
        import contextlib
        import io
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            dests = self.hf.download_genotype_oracle(self.work, repo="local/test")
        self.assertEqual(len(dests), 1)
        self.assertIn("has no phased indel genotypes of hg38_HGDP1kGP", out.getvalue())

    def test_no_genotypes_on_the_hub_is_an_error(self):
        with self.assertRaises(ValueError):
            self.hf.download_genotype_oracle(self.work, repo="local/test")


class TestTheReportSaysWhyTheGateDidNotRun(unittest.TestCase):
    """The inputs card's Carrier gate row for each way the gate did not check the rows."""

    def test_each_status(self):
        import generate_report as gr
        said = {status: gr.render_gate_row({"status": status})
                for status in ("no-oracle", "opt-in", "off", "untiered", "per-sample",
                               "no-panel")}
        self.assertIn("crisprme.py download --what genotype-oracle", said["no-oracle"][1])
        self.assertIn("only on request", said["opt-in"][1])
        self.assertIn("Turned off", said["off"][1])
        self.assertIn("build-index-only --samplesID", said["untiered"][1])
        # a per-sample search's rows name their cis carriers; no phased store, nothing to say
        self.assertIsNone(said["per-sample"])
        self.assertIsNone(said["no-panel"])
        self.assertIsNone(gr.render_gate_row(None))
        # a companion the report could not read
        self.assertIn("not read", gr.render_gate_row({"status": "run", "used": False})[1])


class TheSearchSaysWhatTheGateDoes(bio._Fixture):
    """complete-search as a user runs it (a stub engine): which panel the gate reads, or why
    it does not run, and its flags checked before anything starts."""

    def test_it_names_the_panel_or_says_why_not(self):
        rc, out, err = self.run_bio()  # ds, with its phased genotype store
        self.assertEqual(rc, 0, err)
        rc, out, err = self.run_search("--population-level")
        self.assertIn("Launching job ", out, out + err)
        self.assertIn("[complete-search] Carrier gate: each row's variant combination is "
                      "checked against the phased genotypes of ds; the report ranks a "
                      "combination no haplotype carries after every other row", out)
        for options in (("--population-level", "--no-prioritize-observed"), ("--per-sample",)):
            with self.subTest(options=options):
                rc, out, err = self.run_search(*options)
                self.assertIn("Launching job ", out, out + err)
                self.assertNotIn("Carrier gate", out)
        # the same dataset with the classic dictionaries alone: said, not run
        shutil.rmtree(self.path("Dictionaries", "genotypes_ds"))
        shutil.rmtree(self.path("Dictionaries", "registry_ds"))
        rc, out, err = self.run_search("--population-level")
        self.assertIn("[complete-search] Carrier gate: not run, as the genotyped datasets of "
                      "this search have the classic per-sample dictionaries alone", out)

    def test_its_flags_are_checked_first(self):
        for options, said in (
                (("--prioritize-observed", "--no-prioritize-observed"),
                 "--prioritize-observed and --no-prioritize-observed are mutually exclusive"),
                (("--gate-min-expected", "-1"), "--gate-min-expected needs a number >= 0"),
                (("--gate-min-expected", "three"), "--gate-min-expected needs a number >= 0")):
            with self.subTest(options=options):
                rc, out, err = self.run_search(*options)
                self.assertEqual(rc, 1)
                self.assertIn(said, err)
                self.assertNotIn("Launching job", out)


class TestPlanPanels(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_plan_")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        os.makedirs(os.path.join(self.dicts, "genotypes_hg38_G"))
        reg = os.path.join(self.dicts, "registry_hg38_mega")
        os.makedirs(reg)
        with open(os.path.join(reg, "reg_chr22.idx"), "w") as fh:
            json.dump({"aggregation": "info_af", "data_type": "sites-only"}, fh)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_modes(self):
        d = self.dicts
        self.assertEqual(hp.plan_panels(d, ["hg38_G"], "population-level", "auto"),
                         ("run", [("hg38_G", d, "own")]))
        self.assertEqual(hp.plan_panels(d, ["hg38_G"], "per-sample", "auto"),
                         ("per-sample", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_G"], "population-level", "off"), ("off", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_mega"], "population-level", "auto"),
                         ("opt-in", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_mega"], "population-level", "on"),
                         ("no-oracle", []))
        # a genotyped dataset with the classic dictionaries alone: said, not run
        os.makedirs(os.path.join(d, "dictionaries_hg38_1000G"))
        self.assertEqual(hp.plan_panels(d, ["hg38_1000G"], "population-level", "auto"),
                         ("untiered", []))
        self.assertEqual(hp.plan_panels(d, ["hg38_1000G"], "per-sample", "auto"),
                         ("no-panel", []))

    def test_the_oracle_conventional_or_given(self):
        os.makedirs(os.path.join(self.dicts, f"genotypes_{hp.ORACLE_PANEL}"))
        self.assertEqual(hp.plan_panels(self.dicts, ["hg38_mega"], "population-level", "on"),
                         ("run", [(hp.ORACLE_PANEL, self.dicts, "oracle")]))
        other = os.path.join(self.tmp, "elsewhere")
        os.makedirs(os.path.join(other, "genotypes_hg38_HGDP1kGP_v2"))
        given = os.path.join(other, "genotypes_hg38_HGDP1kGP_v2")
        self.assertEqual(hp.find_oracle(self.dicts, given), ("hg38_HGDP1kGP_v2", other))
        # a genotyped dataset searched beside the mega is read once, as its own panel
        self.assertEqual(
            hp.plan_panels(self.dicts, [hp.ORACLE_PANEL, "hg38_mega"], "population-level", "on"),
            ("run", [(hp.ORACLE_PANEL, self.dicts, "own")]))


if __name__ == "__main__":
    unittest.main()
