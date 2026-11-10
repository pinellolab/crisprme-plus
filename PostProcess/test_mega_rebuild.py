#!/usr/bin/env python3
"""The mega rebuild of CRISPRme 2.10.0: build-index-only on a sites-only dataset, and the
recipe end to end (seq_script/merge_panels/mega_rebuild.sh), with the real crisprme.py and
the stub engine of test_build_index_only_atomic.

  * build-index-only of a sites-only dataset -- VCFs with no sample columns, or with the one
    presence genotype synth_sites_gt.py writes -- says on stdout what builds its registry
    (build_mega_registry.py, from the INFO AF fields), where it printed a WARNING on stderr
    to run it again WITH --samplesID, which a sites-only dataset has none of (the judges of
    D2 saw it on every mega build), and builds no phased indel store of the presence
    genotype (it would name that pseudo-sample every indel's carrier). A genotyped dataset
    without --samplesID still gets the WARNING.
  * mega_rebuild.sh on a two-chromosome install: the build input (PASS and '.' records,
    INFO/AF = AF_max), one build-index-only of the presence-genotype copy, the presence
    sample out of the indel log, the per-sample dictionaries removed, the registry and the
    indel frequencies from the build input, the build input left as the dataset's VCFs, and
    verify_mega_chrom.py's checks of every store against the source VCF passing. The stub
    add-variants codes chr:11 A>G and windows the insertion at 26 whatever the VCF: the
    source VCF holds those as its PASS records, and records of other FILTERs, which the
    real enricher skips, so every store must leave them out.

Run: cd PostProcess && python -m unittest test_mega_rebuild -v
"""

import gzip
import os
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_build_index_only_atomic as bio  # noqa: E402  (the real crisprme.py, a stub engine)
import synth_sites_gt as sg  # noqa: E402

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule

REPO = os.path.dirname(HERE)
RECIPE = os.path.join(REPO, "seq_script", "merge_panels", "mega_rebuild.sh")

try:
    import pysam  # noqa: F401  (BGZF for the recipe's VCFs; tabix through it here)
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False


def write_vcf(path, rows, samples=(), header=()):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if samples:
        cols += "\tFORMAT\t" + "\t".join(samples)
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n" + "".join(h + "\n" for h in header) + cols + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


class ABuildOfASitesOnlyDataset(bio._Fixture):
    """build-index-only of a sites-only dataset says what builds its registry, on stdout."""

    def sites_only(self, presence=False):
        """VCFs/ds as a sites-only dataset: its records with no sample columns, or with the
        presence genotype synth_sites_gt.py writes."""
        for chrom, seq in bio.GENOME.items():
            path = self.path("VCFs", "ds", f"ds.{chrom}.vcf.gz")
            rows = [[chrom, "11", ".", "A", "G", ".", "PASS", "AF_gnomAD=0.25;AF_max=0.25"],
                    [chrom, "26", ".", seq[25], seq[25] + "T", ".", "PASS",
                     "AF_gnomAD=0.25;AF_max=0.25"]]
            if presence:
                write_vcf(path, [r + ["GT", "0/1"] for r in rows], samples=("MEGA",),
                          header=(sg.SITES_ONLY_MARK,))
            else:
                write_vcf(path, rows)

    def build(self):
        command = ["build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam",
                   bio.PAM, "--bDNA", "1", "--bRNA", "1", "--thread", "2", "--path", self.work,
                   "--vcf", self.path("VCFs", "ds")]
        return self.finish(self.start(command=command))

    def test_vcfs_with_no_samples(self):
        self.sites_only()
        rc, out, err = self.build()
        self.assertEqual(rc, 0, out + err)
        self.assertNotIn("INCOMPLETE VARIANT INDEX", err)
        self.assertNotIn("WITH --samplesID", err)
        self.assertEqual(err, "")
        self.assertIn("NOTE [build-index-only]: ds is a sites-only dataset (its VCFs have no "
                      "sample columns)", out)
        self.assertIn("PostProcess/build_mega_registry.py", out)
        self.assertIn("Its indels are searchable only from a VCF carrying the presence "
                      "genotype", out)
        self.assertFalse(os.path.exists(self.path("Dictionaries", "indel_genotypes_ds")))
        self.assertTrue(os.path.isdir(self.lib(bio.VAR)))

    def test_the_presence_genotype(self):
        self.sites_only(presence=True)
        rc, out, err = self.build()
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(err, "")
        self.assertIn("NOTE [build-index-only]: ds is a sites-only dataset (its one sample is "
                      "the presence genotype PostProcess/synth_sites_gt.py writes)", out)
        self.assertIn("Strip the presence sample from Dictionaries/log_indels_ds/", out)
        # no phased indel store naming the pseudo-sample every indel's carrier
        self.assertFalse(os.path.exists(self.path("Dictionaries", "indel_genotypes_ds")))
        self.assertTrue(os.path.isdir(self.lib(bio.INDELS)))

    def test_a_genotyped_dataset_without_samples_is_warned_as_before(self):
        rc, out, err = self.build()
        self.assertEqual(rc, 0, out + err)
        self.assertIn("WARNING: INCOMPLETE VARIANT INDEX", err)
        self.assertNotIn("is a sites-only dataset", out)


TABIX_STUB = r'''#!PYTHON
"""tabix -f -p vcf FILE, through pysam (no htslib on this host)."""
import sys
import pysam
pysam.tabix_index(sys.argv[-1], preset="vcf", force=True, keep_original=True)
'''

CRISPRME_WRAPPER = r'''#!/bin/sh
exec "PYTHON" "CRISPRME" "$@"
'''


@unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF and tabix) not installed")
class TheRecipe(bio._Fixture):
    """mega_rebuild.sh end to end on two chromosomes, the stub engine under it."""

    def setUp(self):
        super().setUp()
        shutil.rmtree(self.path("VCFs", "ds"))
        self.src = os.path.join(self.root, "merged")
        for chrom, seq in bio.GENOME.items():
            write_vcf(os.path.join(self.src, f"mega.{chrom}.afmax.vcf.gz"), [
                # what the stub enricher codes and windows: the PASS records
                [chrom, "11", "rs1", "A", "G", ".", "PASS",
                 "AF_gnomAD=0.2;AF_HGDP=0.25;AF_max=0.25"],
                [chrom, "26", "rs2", seq[25], seq[25] + "T", ".", "PASS",
                 "AF_TOPMed=0.25;AF_max=0.25"],
                # records the enricher skips, so no store may hold
                [chrom, "31", "rs3", seq[30], "T" if seq[30] != "T" else "A", ".",
                 "LOW_VQSLOD", "AF_HGDP=0.5;AF_max=0.5"],
                [chrom, "41", "rs4", seq[40] + seq[41], seq[40], ".", "ExcHet;LOW_VQSLOD",
                 "AF_HGDP=0.5;AF_max=0.5"],
            ], header=('##INFO=<ID=AF_max,Number=A,Type=Float,Description="max">',))
        os.makedirs(self.path("PAMs"), exist_ok=True)
        shutil.copy(bio.PAM, self.path("PAMs", os.path.basename(bio.PAM)))
        for name, text in (("tabix", TABIX_STUB), ("crisprme.py", CRISPRME_WRAPPER)):
            p = os.path.join(self.bin, name)
            with open(p, "w") as fh:
                fh.write(text.replace("#!PYTHON", "#!" + sys.executable)
                         .replace("PYTHON", sys.executable).replace("CRISPRME", bio.CRISPRME))
            os.chmod(p, 0o755)

    def rebuild(self):
        env = {k: v for k, v in os.environ.items() if not k.startswith(("STUB_", "CRISPRME_"))}
        env.update(PATH=os.pathsep.join([self.bin, os.environ.get("PATH", "")]),
                   PYTHONPATH=os.path.join(bio.SHARED, "py"),
                   MPLCONFIGDIR=os.path.join(bio.SHARED, "mpl"), STUB_LOG=self.log,
                   STUB_MARK=self.mark, STUB_RELEASE=self.release,
                   PAM="PAMs/" + os.path.basename(bio.PAM), BDNA="1", BRNA="1", THREADS="2",
                   PYTHON=sys.executable, PP=HERE)
        return subprocess.run(["bash", RECIPE, self.work, self.src], env=env,
                              capture_output=True, text=True, timeout=600)

    def test_every_store_reads_the_records_the_enricher_reads(self):
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("SUMMARY: 2 chroms, 2 SNP alleles + 2 indel-path alleles, ALL PASS",
                      done.stdout)
        self.assertNotIn("INCOMPLETE VARIANT INDEX", done.stderr)
        self.assertIn("NOTE [build-index-only]: hg38_mega is a sites-only dataset (its one "
                      "sample is the presence genotype", done.stdout)
        d = self.path("Dictionaries")
        self.assertEqual(sorted(os.listdir(d)),
                         ["indel_af_hg38_mega", "log_indels_hg38_mega", "registry_hg38_mega"])
        # the presence sample out of the log (the stub names S2 there)
        for chrom in bio.GENOME:
            with gzip.open(os.path.join(d, "log_indels_hg38_mega", f"log{chrom}.txt.gz"),
                           "rt") as fh:
                rows = [ln.rstrip("\n").split("\t") for ln in fh][1:]
            self.assertEqual([r[1] for r in rows], [""])
        # the build input is the dataset's VCFs: its PASS records, INFO/AF = AF_max
        vdir = self.path("VCFs", "hg38_mega")
        self.assertEqual(sorted(f for f in os.listdir(vdir) if f.endswith(".vcf.gz")),
                         ["mega.chr1.build.vcf.gz", "mega.chr2.build.vcf.gz"])
        with gzip.open(os.path.join(vdir, "mega.chr1.build.vcf.gz"), "rt") as fh:
            body = [ln.split("\t") for ln in fh if not ln.startswith("#")]
        self.assertEqual([(r[1], r[6], r[7].rstrip("\n").rsplit(";", 1)[1]) for r in body],
                         [("11", "PASS", "AF=0.25"), ("26", "PASS", "AF=0.25")])
        self.assertFalse(os.path.exists(self.path(".mega_rebuild_hg38_mega")))
        self.assertTrue(os.path.isdir(self.lib("20-NRG_2_hg38+hg38_mega")))
        self.assertTrue(os.path.isdir(self.lib("20-NRG_2_hg38+hg38_mega_INDELS")))

    # the samplesID file the shipped mega install holds: its five sources as pseudo-samples
    SAMPLES_ID = ("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                  "1000G2021\t1000G2021\t1000G2021\tNA\nHGDP\tHGDP\tHGDP\tNA\n"
                  "gnomAD\tgnomAD\tgnomAD\tNA\nTOPMed\tTOPMed\tTOPMed\tNA\nAoU\tAoU\tAoU\tNA\n")

    def test_the_install_has_the_samples_id_file_a_search_reads(self):
        # a search names the dataset's samplesID file (--samplesID), and the post-analysis
        # reads it (a search of the gate failed at hg38_mega without it): the recipe writes
        # the shipped one, the registry's sources as pseudo-samples, and keeps one there
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        sid = self.path("samplesIDs", "hg38_mega.samplesID.txt")
        with open(sid) as fh:
            self.assertEqual(fh.read(), self.SAMPLES_ID)
        import build_mega_registry as bmr
        self.assertEqual([ln.split("\t")[0] for ln in self.SAMPLES_ID.splitlines()[1:]],
                         list(bmr.DEFAULT_DATASET_META))

    def test_a_samples_id_file_already_there_is_kept(self):
        os.makedirs(self.path("samplesIDs"), exist_ok=True)
        sid = self.path("samplesIDs", "hg38_mega.samplesID.txt")
        with open(sid, "w") as fh:
            fh.write(self.SAMPLES_ID.replace("\tNA\n", "\tunknown\n"))
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        with open(sid) as fh:
            self.assertEqual(fh.read(), self.SAMPLES_ID.replace("\tNA\n", "\tunknown\n"))

    def test_a_second_record_of_one_allele_is_left_out_of_every_part(self):
        """One REF>ALT on two PASS records at a position (the as-run merge wrote 8 indels
        so, genome-wide, long deletions with one AF_max on both lines): the build input
        keeps the more frequent record (INFO/AF = AF_max; on a tie the earlier), so the
        enricher, which gives each indel record a window and a log row of its own, the
        indel frequencies and the registry read one record each, and the checks pass. With
        both records the checks of chr21 failed ("1 allele(s) on more than one row";
        review r1 of lane M, M-R1-3). (The stub enricher windows one insertion whatever the
        VCF: the build input is what the real one reads.)"""
        for chrom, seq in bio.GENOME.items():
            path = os.path.join(self.src, f"mega.{chrom}.afmax.vcf.gz")
            with gzip.open(path, "rt") as fh:
                lines = fh.read().splitlines()
            lines += [
                # the insertion at 26 again, its AF_max the first one's (a tie)
                "\t".join([chrom, "26", ".", seq[25], seq[25] + "T", ".", "PASS",
                           "AF_HGDP=0.125;AF_max=0.25"]),
                # the SNV at 11 again, more frequent: its frequencies are the allele's
                "\t".join([chrom, "11", "rs9", "A", "G", ".", ".",
                           "AF_HGDP=0.5;AF_max=0.5"]),
            ]
            head = [ln for ln in lines if ln.startswith("#")]
            body = sorted((ln for ln in lines if not ln.startswith("#")),
                          key=lambda ln: int(ln.split("\t")[1]))
            write_vcf(path, [ln.split("\t") for ln in body],
                      header=[h for h in head if h.startswith("##") and
                              not h.startswith("##fileformat")])
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("SUMMARY: 2 chroms, 2 SNP alleles + 2 indel-path alleles, ALL PASS",
                      done.stdout)
        self.assertIn("second_records_of_one_allele=2", done.stdout)
        self.assertIn("2 left out as a second record of one REF>ALT", done.stderr)
        vdir = self.path("VCFs", "hg38_mega")
        with gzip.open(os.path.join(vdir, "mega.chr1.build.vcf.gz"), "rt") as fh:
            body = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#")]
        self.assertEqual([(r[1], r[2], r[7]) for r in body], [
            ("11", "rs9", "AF_HGDP=0.5;AF_max=0.5;AF=0.5"),
            ("26", "rs2", "AF_TOPMed=0.25;AF_max=0.25;AF=0.25")])
        import build_mega_indel_af as bmi
        store = bmi.IndelAfReader(self.path("Dictionaries", "indel_af_hg38_mega",
                                            "indel_af_chr1.tsv.gz"))
        seq = bio.GENOME["chr1"]
        self.assertEqual(store.lookup(26, seq[25], seq[25] + "T"),
                         {"TOPMed": 0.25, "AF_max": 0.25})
        from tier0_registry import RegistryReader, GLOBAL_GROUP_ID
        reg = RegistryReader(self.path("Dictionaries", "registry_hg38_mega", "reg_chr1.bin"),
                             self.path("Dictionaries", "registry_hg38_mega", "reg_chr1.idx"))
        self.assertAlmostEqual(reg.lookup(11, "G")[GLOBAL_GROUP_ID].allele_freq(), 0.5,
                               places=3)
        reg.close()

    def test_another_pam_on_the_finished_install_names_no_samples_id(self):
        """The recipe's way to another PAM is build-index-only on the finished install
        (its header): that build keeps the registry and the indel frequencies in place and
        says where a sites-only dataset's come from, not to run it again with --samplesID,
        which such a dataset has none of (review r1 of lane M, M-R1-N2)."""
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        nrr = os.path.join(REPO, "PAMs", "20bp-NRR-SpCas9.txt")
        shutil.copy(nrr, self.path("PAMs", os.path.basename(nrr)))
        rc, out, err = self.finish(self.start(command=[
            "build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam",
            self.path("PAMs", os.path.basename(nrr)), "--bDNA", "1", "--bRNA", "1",
            "--thread", "2", "--path", self.work, "--vcf", self.path("VCFs", "hg38_mega")]))
        self.assertEqual((rc, err), (0, ""), out)
        self.assertIn("NOTE [build-index-only]: hg38_mega is a sites-only dataset: its "
                      "registry", out)
        self.assertIn("not from --samplesID", out)
        self.assertNotIn("run this again with --samplesID", out)
        self.assertTrue(os.path.isdir(self.lib("20-NRR_2_hg38+hg38_mega")))

    def test_the_checks_fail_a_registry_that_reads_every_record(self):
        # the shipped mega's registry read every record whatever its FILTER: the checks
        # of the recipe catch it (the LOW_VQSLOD SNV at 31 is no allele the genome codes)
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        import build_mega_registry as bmr
        import verify_mega_chrom as vm
        reg = self.path("Dictionaries", "registry_hg38_mega")
        everything = os.path.join(self.root, "every_record.vcf.gz")
        with gzip.open(os.path.join(self.src, "mega.chr1.afmax.vcf.gz"), "rt") as fh, \
                gzip.open(everything, "wt") as out:
            for line in fh:
                f = line.split("\t")
                out.write(line if line.startswith("#") else "\t".join(f[:6] + ["PASS"] + f[7:]))
        bmr.build(everything, "chr1", reg)
        ok, _ns, _ni = vm.verify(self.src, reg, self.path("Dictionaries", "indel_af_hg38_mega"),
                                 "chr1")
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main(verbosity=2)
