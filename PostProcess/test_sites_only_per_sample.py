#!/usr/bin/env python3
"""A sites-only dataset's post-analysis is the population-level one, whatever the mode.

A sites-only dataset (the mega panel) has no genotypes for a per-sample search to resolve.
Its rows went through the per-sample path's enumeration all the same, which joins two
variants only through carriers they share, and it has none: a site that takes two variants
together lost its row, and a site scored lower. The command line said the population-level
analysis runs instead, and a search of the default dataset beside hg38_mega runs per-sample
by default. Each dataset of a search has a post-analysis of its own, so now, per dataset,
from its registry's data_type, a sites-only one goes the population-level way
(new_simple_analysis.py), and complete-search says so (crisprme.py).

The first class runs the real new_simple_analysis.py on a sites-only registry that the real
build_mega_registry.py builds; the second, the real crisprme.py complete-search against the
stub engine of test_build_index_only_atomic, up to the search it launches.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import build_mega_registry
import test_build_index_only_atomic as bio

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule

HERE = os.path.dirname(os.path.abspath(__file__))
GUIDE = "CTAACAGTTGCTTTTATCAC"
CRRNA = GUIDE + "NNN"
# 2 mismatches from the guide (at its columns 4 and 10), with a variant at each that restores
# it: within a budget of no mismatch only with both variants on the one row
SITE = "CTAAAAGTTGGTTTTATCACTGG"
FLANK = "ACGTGCATTC" * 10
GENOME = ">chrT\n" + FLANK + SITE + FLANK + "\n"
VCF_HEADER = ("##fileformat=VCFv4.2\n"
              '##INFO=<ID=AF_1000G2021,Number=A,Type=Float,Description="AF">\n'
              '##INFO=<ID=AF_HGDP,Number=A,Type=Float,Description="AF">\n'
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


class ThePostAnalysisOfASitesOnlyDataset(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="sites_only_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.fasta = os.path.join(self.tmp, "chrT.fa")
        with open(self.fasta, "w") as fh:
            fh.write(GENOME)
        vcf = os.path.join(self.tmp, "mega.chrT.vcf")
        at = len(FLANK)  # 0-based start of the site
        with open(vcf, "w") as fh:
            fh.write(VCF_HEADER)
            # the two variants of the site, each making one mismatch a match
            for col, ref, alt in ((4, "A", "C"), (10, "G", "C")):
                fh.write(f"chrT\t{at + col + 1}\trsT{col}\t{ref}\t{alt}\t.\tPASS\t"
                         "AF_1000G2021=0.1;AF_HGDP=0.2\n")
        self.dicts = os.path.join(self.tmp, "Dictionaries")
        build_mega_registry.build(vcf, "chrT", os.path.join(self.dicts, "registry_mega"))
        self.pam = os.path.join(self.tmp, "pam.txt")
        with open(self.pam, "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        # the variant search's row: the window, an IUPAC code at each variant column
        self.targets = os.path.join(self.tmp, "chrT.targets.txt")
        window = SITE[:4] + "M" + SITE[5:10] + "S" + SITE[11:]
        row = ["X", CRRNA, window, "chrT", str(at), str(at), "+", "0", "0", "0", "NA", "y",
               "NA", "NA", "NA", "NA", "NA", "NA"]
        with open(self.targets, "w") as fh:
            fh.write("#header\n" + "\t".join(row) + "\n")

    def analyse(self, mode):
        out = os.path.join(self.tmp, f"out_{mode}")
        env = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
        env.update(CRISPRME_FAST_MODE=mode, CRISPRME_SKIP_SCORER_ENV="1")
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "new_simple_analysis.py"), self.fasta,
             os.path.join(self.dicts, "dictionaries_mega", "my_dict_chrT.json"), self.targets,
             self.pam, out, "0"],
            cwd=HERE, env=env, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        rows = {}
        for ext in ("bestCFD", "bestmmblg", "bestCRISPR_BULGE"):
            with open(f"{out}.{ext}.txt") as fh:
                rows[ext] = fh.read()
        return rows, done.stdout

    def test_it_is_the_same_in_either_mode_and_keeps_the_site_that_takes_both_variants(self):
        population, _ = self.analyse("1")
        per_sample, said = self.analyse("0")
        for rows in (population, per_sample):
            for ext, text in rows.items():
                lines = text.splitlines()[1:]
                with self.subTest(ext):
                    # both variants on the one row: the guide with no mismatch
                    both = [l.split("\t") for l in lines
                            if l.split("\t")[2].upper() == GUIDE + "TGG"]
                    self.assertEqual(len(both), 1, lines)
                    self.assertEqual(both[0][7], "0")
                    self.assertEqual(sorted(both[0][15].split(",")), ["rsT10", "rsT4"])
        self.assertEqual(per_sample, population)
        self.assertIn("Sites-only dataset on chrT: population-level post-analysis.", said)


class TheModeOfASearchWithASitesOnlyDataset(bio._Fixture):
    """complete-search as a user runs it next to the indexes: ds, a phased dataset that
    build-index-only built (its tiers), and a sites-only one copied from it."""

    def add_sites_only(self, name, manifest=None):
        for tier in ("dictionaries", "log_indels"):
            shutil.copytree(self.path("Dictionaries", f"{tier}_ds"),
                            self.path("Dictionaries", f"{tier}_{name}"))
        registry = self.path("Dictionaries", f"registry_{name}")
        os.makedirs(registry)
        for chrom in ("chr1", "chr2"):
            with open(os.path.join(registry, f"reg_{chrom}.idx"), "w") as fh:
                json.dump(manifest or {"data_type": "sites-only", "aggregation": "info_af"}, fh)
            open(os.path.join(registry, f"reg_{chrom}.bin"), "w").close()
        for suffix in ("", "_INDELS"):
            shutil.copytree(self.lib(bio.VAR + suffix), self.lib(f"NRG_2_hg38+{name}{suffix}"))
        for suffix in ("", "_INDELS"):
            shutil.copytree(self.path("Genomes", f"hg38+ds{suffix}"),
                            self.path("Genomes", f"hg38+{name}{suffix}"))

    def mode_of(self, vcfs, *options):
        rc, out, err = self.run_search(*options, vcfs=vcfs)
        self.assertIn("Launching job ", out, out + err)
        with open(self.path("Results", "search", ".search_mode")) as fh:
            mode = fh.read()
        with open(self.path("Results", "search", ".alt_mode")) as fh:
            alt = fh.read()
        return mode, alt, out

    def test_phased_dataset_runs_per_sample_beside_it_whichever_name_sorts_first(self):
        self.run_bio()
        for name in ("aaa_sites", "zzz_sites"):
            with self.subTest(name):
                self.add_sites_only(name)
                mode, alt, out = self.mode_of(f"ds\n{name}\n")
                self.assertEqual((mode, alt), ("per-sample", "alt-alignments"))
                self.assertIn("genotyped + PHASED index detected", out)
                self.assertIn(f"[complete-search] Sites-only, with no genotypes to resolve: {name}. "
                              "A sites-only dataset's rows are the population-level analysis's", out)
                self.assertNotIn("WARNING [complete-search]: --per-sample", out)

    def test_per_sample_on_sites_only_datasets_alone_is_the_population_level_analysis(self):
        self.run_bio()
        self.add_sites_only("mega")
        mode, alt, out = self.mode_of("mega\n", "--per-sample")
        self.assertEqual((mode, alt), ("population-level", "no-alt-alignments"))
        self.assertIn("WARNING [complete-search]: --per-sample has nothing to resolve here, as what it "
                      "searches is sites-only, with no genotypes (mega): the population-level analysis "
                      "runs instead.", out)
        self.assertIn("[complete-search] POPULATION-LEVEL analysis", out)
        self.assertNotIn("PER-SAMPLE analysis", out)
        # with no mode flag, as it was
        mode, alt, out = self.mode_of("mega\n")
        self.assertEqual((mode, alt), ("population-level", "no-alt-alignments"))
        self.assertNotIn("sites-only", out.replace("sites-only / unphased", ""))

    def test_a_registry_with_no_data_type_is_sites_only_by_its_aggregation(self):
        # a mega index built before the manifest had one, and not stamped since
        self.run_bio()
        self.add_sites_only("mega", {"aggregation": "info_af"})
        mode, alt, out = self.mode_of("mega\n", "--per-sample")
        self.assertEqual((mode, alt), ("population-level", "no-alt-alignments"))
        self.assertIn("WARNING [complete-search]: --per-sample has nothing to resolve here", out)
        mode, alt, out = self.mode_of("ds\nmega\n")
        self.assertEqual((mode, alt), ("per-sample", "alt-alignments"))
        self.assertIn("[complete-search] Sites-only, with no genotypes to resolve: mega.", out)

    def test_per_sample_on_the_two_says_which_goes_which_way(self):
        self.run_bio()
        self.add_sites_only("mega")
        mode, alt, out = self.mode_of("ds\nmega\n", "--per-sample")
        self.assertEqual((mode, alt), ("per-sample", "alt-alignments"))
        self.assertIn("[complete-search] Sites-only, with no genotypes to resolve: mega.", out)
        self.assertNotIn("WARNING [complete-search]: --per-sample", out)
        mode, alt, out = self.mode_of("ds\nmega\n", "--population-level")
        self.assertEqual((mode, alt), ("population-level", "no-alt-alignments"))
        self.assertNotIn("[complete-search] Sites-only", out)

    def test_no_sites_only_dataset_no_word_of_it(self):
        self.run_bio()
        mode, alt, out = self.mode_of("ds\n")
        self.assertEqual((mode, alt), ("per-sample", "alt-alignments"))
        self.assertNotIn("Sites-only, with no genotypes", out)
        self.assertNotIn("WARNING [complete-search]: --per-sample", out)


if __name__ == "__main__":
    unittest.main()
