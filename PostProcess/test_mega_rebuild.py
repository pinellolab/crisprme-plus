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
import re
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_build_index_only_atomic as bio  # noqa: E402  (the real crisprme.py, a stub engine)
import synth_sites_gt as sg  # noqa: E402
from test_indel_window import pasted, shell_errors  # noqa: E402  (bash and zsh, as pasted)

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


class _Recipe(bio._Fixture):
    """The recipe's fixture: the two-chromosome genome, the merged VCFs of each chromosome and
    the stub engine (TheRecipe, TheRecipeInAFolderOfTheDataset)."""

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

    def rebuild(self, install=None):
        """The recipe, run in the working directory, or in ``install``."""
        env = {k: v for k, v in os.environ.items() if not k.startswith(("STUB_", "CRISPRME_"))}
        env.update(PATH=os.pathsep.join([self.bin, os.environ.get("PATH", "")]),
                   PYTHONPATH=os.path.join(bio.SHARED, "py"),
                   MPLCONFIGDIR=os.path.join(bio.SHARED, "mpl"), STUB_LOG=self.log,
                   STUB_MARK=self.mark, STUB_RELEASE=self.release,
                   PAM="PAMs/" + os.path.basename(bio.PAM), BDNA="1", BRNA="1", THREADS="2",
                   PYTHON=sys.executable, PP=HERE)
        return subprocess.run(["bash", RECIPE, install or self.work, self.src], env=env,
                              capture_output=True, text=True, timeout=600)


@unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF and tabix) not installed")
class TheRecipe(_Recipe):
    """mega_rebuild.sh end to end on two chromosomes, the stub engine under it."""

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


    def test_the_checks_fail_a_global_group_counted_over_another_allele_number(self):
        """The GLOBAL group is the largest of the sources' frequencies as a count over the
        summed nominal AN of the sources that carry the allele, which the carrier gate's
        frequency bound reads (D6): the checks hold it to that count exactly. A registry
        whose allele numbers are twice the nominal ones keeps every source's frequency
        within 0.5/AN, and its GLOBAL frequency within 1e-3 of AF_max, which the checks let
        through up to review r1 of lane M (M-R1-N3); they fail it now."""
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        import build_mega_registry as bmr
        import verify_mega_chrom as vm
        reg = self.path("Dictionaries", "registry_hg38_mega")
        meta = {ds: dict(m, an_nominal=2 * m["an_nominal"])
                for ds, m in bmr.DEFAULT_DATASET_META.items()}
        build_vcf = self.path("VCFs", "hg38_mega", "mega.chr1.build.vcf.gz")
        bmr.build(build_vcf, "chr1", reg, dataset_meta=meta)
        ok, _ns, _ni = vm.verify(self.src, reg, self.path("Dictionaries", "indel_af_hg38_mega"),
                                 "chr1")
        self.assertFalse(ok)
        # and the recipe's own registry passes
        bmr.build(build_vcf, "chr1", reg)
        ok, _ns, _ni = vm.verify(self.src, reg, self.path("Dictionaries", "indel_af_hg38_mega"),
                                 "chr1")
        self.assertTrue(ok)


class TheChecksReadTheVcfInUpperCase(unittest.TestCase):
    """CRISPRitz 2.10.0 reads a VCF's REF and every ALT in upper case, as VCF 4.x does (its
    docs/MULTIALLELIC_ENRICHMENT.md sections 2.1 and 2.9): it codes a lower-case SNV ALT and
    names a record written in lower case in upper case in its indel log, and the recipe's
    checks (verify_mega_chrom.py) read the source VCF so too. The mega's merged VCFs hold 13
    deletions written in mixed case (chr3, 4, 5, 10, 12, 18 and 21): read as written, the log
    check failed each of those chromosomes ("indel log alleles ...: N more, N missing"), and
    the check of the build input let a REF>ALT on two records in two cases through, which the
    enricher would window twice. The stores here are written by the recipe's own builders
    (build_mega_registry.py, build_mega_indel_af.py, prepare_mega_build_vcf.py), the indel log
    and the enriched genome as CRISPRitz 2.10.0 writes them."""

    SEQ = ("GGCATTACGGATACCATGCAATTGACCGGTATCAGCTAGCCTTAGGCATCGATGCTAGCCTTACTGATCG"
           "TACGTAGCTAGTACGTTGCAATTGGCCAATTGCATTGACCGGTA")

    def setUp(self):
        import tempfile
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        seq = self.SEQ
        assert (seq[10], seq[39], seq[59:62]) == ("A", "C", "CTT"), (seq[10], seq[39],
                                                                    seq[59:62])
        self.gw = os.path.join(self.d, "gw")
        write_vcf(os.path.join(self.gw, "mega.chr1.afmax.vcf.gz"), [
            ["chr1", "11", "rs1", "A", "T", ".", "PASS", "AF_HGDP=0.1;AF_max=0.1"],
            # an SNV written in lower case: CRISPRitz 2.10.0 codes it (S at 40)
            ["chr1", "40", "rs2", "c", "g", ".", "PASS", "AF_gnomAD=0.02;AF_max=0.02"],
            # a deletion in mixed case, as the mega's 13 are: its log row is upper case
            ["chr1", "60", "rs3", "cTT", "c", ".", "PASS", "AF_gnomAD=0.01;AF_max=0.01"],
        ], header=('##INFO=<ID=AF_max,Number=A,Type=Float,Description="max">',))
        for name, text in (("ref", seq), ("enr", seq[:10] + "W" + seq[11:39] + "S" + seq[40:])):
            os.makedirs(os.path.join(self.d, name))
            with open(os.path.join(self.d, name, "chr1.fa" if name == "ref"
                                   else "chr1.enriched.fa"), "w") as fh:
                fh.write(">chr1\n%s\n" % text)

    def log(self, indel):
        """An indel log of the chromosome (the presence sample out), one row naming the
        deletion as ``indel``, its window 30/30 (CRISPRitz 2.10.0)."""
        import tempfile
        d = tempfile.mkdtemp(prefix="log_", dir=self.d)
        with open(os.path.join(d, "logchr1.txt"), "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
                     "chr1_29-92_1\t\trs3\t0.01\tchr1_60_%s\t0,1\tN\n" % indel)
        return d

    def stores(self, source):
        """The registry, the indel frequencies and the build input of ``source`` (a merged
        VCF), built as the recipe builds them, from the build input."""
        import build_mega_registry as bmr
        import build_mega_indel_af as bmi
        import prepare_mega_build_vcf as pm
        build = os.path.join(self.d, "build")
        os.makedirs(build, exist_ok=True)
        plain = os.path.join(build, "mega.chr1.build.vcf")
        pm.prepare(source, plain)
        with open(plain, "rb") as fh, gzip.open(plain + ".gz", "wb") as out:
            out.write(fh.read())
        os.remove(plain)
        reg, ind = os.path.join(self.d, "reg"), os.path.join(self.d, "ind")
        bmr.build(plain + ".gz", "chr1", reg)
        bmi.build_indel_af_store(plain + ".gz", "chr1", ind, window_left=31)
        return reg, ind, build

    def verify(self, log, build=None, reg=None, ind=None):
        import contextlib
        import io
        import verify_mega_chrom as vm
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            ok, ns, ni = vm.verify(self.gw, reg, ind, "chr1", log=log,
                                   enriched=os.path.join(self.d, "enr"),
                                   reference=os.path.join(self.d, "ref"), build_dir=build)
        return ok, ns, ni, out.getvalue()

    def test_records_written_in_lower_case_match_the_upper_case_stores(self):
        reg, ind, build = self.stores(os.path.join(self.gw, "mega.chr1.afmax.vcf.gz"))
        ok, ns, ni, out = self.verify(self.log("CTT_C"), build, reg, ind)
        self.assertTrue(ok, out)
        self.assertEqual((ns, ni), (2, 1), out)
        self.assertIn("chr1: PASS", out)
        self.assertIn("coded_positions=2", out)
        # a log that names the deletion as the VCF writes it (an enricher before CRISPRitz
        # 2.10.0's case rule, whose genome codes no lower-case SNV either) is no log of
        # these records
        ok, _ns, _ni, out = self.verify(self.log("cTT_c"), build, reg, ind)
        self.assertFalse(ok, out)
        self.assertIn("indel log alleles 1, expected 1: 1 more ([(60, 'cTT', 'c')]), "
                      "1 missing ([(60, 'CTT', 'C')])", out)

    def test_a_ref_alt_on_two_records_in_two_cases_is_caught_in_the_build_input(self):
        """verify_mega_chrom.check_build_vcf keys each ALT of the build input as VCF 4.x
        reads it: C>G and c>g on two records is one REF>ALT on more than one record, as two
        records written alike are (prepare_mega_build_vcf.py keeps one record of each). Keyed
        as written, the check let it through."""
        import verify_mega_chrom as vm
        path = os.path.join(self.d, "twin.build.vcf.gz")
        write_vcf(path, [
            ["chr1", "40", "rs2", "C", "G", ".", "PASS", "AF_max=0.02;AF=0.02"],
            ["chr1", "40", "rs4", "c", "g", ".", "PASS", "AF_max=0.01;AF=0.01"],
            ["chr1", "60", "rs3", "CTT", "C", ".", "PASS", "AF_max=0.01;AF=0.01"],
            ["chr1", "60", "rs5", "ctt", "c,CT", ".", "PASS", "AF_max=0.01,0.1;AF=0.01,0.1"],
            # a symbolic allele is an ID, read as written: no repeat
            ["chr1", "70", ".", "A", "<DEL>", ".", "PASS", "AF_max=0.01;AF=0.01"],
            ["chr1", "70", ".", "A", "<del>", ".", "PASS", "AF_max=0.01;AF=0.01"],
        ])
        errs = []
        vm.check_build_vcf(path, 6, errs)
        self.assertEqual(errs, ["build input: 2 REF>ALT(s) on more than one record"])
        # the same twins in the source VCF: the build input keeps one record of each, and
        # every check passes
        write_vcf(os.path.join(self.gw, "mega.chr1.afmax.vcf.gz"), [
            ["chr1", "11", "rs1", "A", "T", ".", "PASS", "AF_HGDP=0.1;AF_max=0.1"],
            ["chr1", "40", "rs2", "C", "G", ".", "PASS", "AF_gnomAD=0.02;AF_max=0.02"],
            ["chr1", "40", "rs4", "c", "g", ".", "PASS", "AF_gnomAD=0.01;AF_max=0.01"],
            ["chr1", "60", "rs3", "cTT", "c", ".", "PASS", "AF_gnomAD=0.01;AF_max=0.01"],
            ["chr1", "60", "rs6", "CTT", "C", ".", "PASS", "AF_gnomAD=0.005;AF_max=0.005"],
        ], header=('##INFO=<ID=AF_max,Number=A,Type=Float,Description="max">',))
        reg, ind, build = self.stores(os.path.join(self.gw, "mega.chr1.afmax.vcf.gz"))
        with gzip.open(os.path.join(build, "mega.chr1.build.vcf.gz"), "rt") as fh:
            body = [ln.split("\t")[:5] for ln in fh if not ln.startswith("#")]
        self.assertEqual(body, [["chr1", "11", "rs1", "A", "T"], ["chr1", "40", "rs2", "C", "G"],
                                ["chr1", "60", "rs3", "cTT", "c"]])
        ok, ns, ni, out = self.verify(self.log("CTT_C"), build, reg, ind)
        self.assertTrue(ok, out)
        self.assertIn("second_records_of_one_allele=2", out)


def _listing(root):
    """Every file, folder and link below ``root``, as relative paths (a link's own name)."""
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        for n in dirnames + filenames:
            out.append(os.path.relpath(os.path.join(dirpath, n), root))
    return sorted(out)


def _content(path):
    """A file's bytes, those it holds uncompressed where it is gzip or BGZF (whose header
    holds the time it was written)."""
    with open(path, "rb") as fh:
        data = fh.read()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


def _section_12_commands():
    """The commands of docs/TROUBLESHOOTING.md section 12, "Building the mega's kit again":
    the first bash block after that paragraph's title."""
    with open(os.path.join(REPO, "docs", "TROUBLESHOOTING.md")) as fh:
        text = fh.read()
    start = text.index("**Building the mega's kit again.**", text.index("## 12. "))
    block = text.index("```bash\n", start) + len("```bash\n")
    return text[block:text.index("```", block)]


@unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF and tabix) not installed")
class TheRecipeInAFolderOfTheDataset(_Recipe):
    """IR1-F2, the integration review r1 of 2.10.0: the advice on a sites-only dataset named
    the recipe after a removal that leaves the dataset's VCFs, registry and indel frequencies in
    the data folder (they come from its VCFs alone), and the recipe, run there with the merged
    VCF as its source, refused the folder for its VCFs, after it had made its staging folder
    (.mega_rebuild_hg38_mega/prepared) and Dictionaries there. It refuses a working directory
    that holds the dataset's VCFs, or any part it builds, before it makes anything; and the flow
    the advice gives, the recipe in a new working directory, then its parts moved in place of
    the dataset's (docs/TROUBLESHOOTING.md section 12, its commands run here as written, the
    chromosomes and the PAM those of this genome), puts a kit back that a build reads. W1C-B1,
    the review of wave 1c-a: so it does from the state a search's note of windows too short
    leaves, the kit in place, whichever removal of the note's is run first."""

    PARTS = [("Genomes", "hg38+hg38_mega"), ("Genomes", "hg38+hg38_mega_INDELS"),
             ("Dictionaries", "dictionaries_hg38_mega"), ("Dictionaries", "log_indels_hg38_mega"),
             ("Dictionaries", "registry_hg38_mega"), ("Dictionaries", "indel_af_hg38_mega"),
             ("Dictionaries", "genotypes_hg38_mega"), ("Dictionaries", "indel_genotypes_hg38_mega"),
             ("genome_library", "20-NRG_2_hg38+hg38_mega"),
             ("genome_library", "20-NRG_2_hg38+hg38_mega_INDELS")]

    def the_kit(self, root):
        """The parts of hg38_mega the recipe builds, in ``root``, as {relative path: bytes}."""
        out = {}
        for top in ("Genomes", "Dictionaries", "genome_library", "VCFs", "samplesIDs"):
            for rel in _listing(os.path.join(root, top)):
                path = os.path.join(root, top, rel)
                if "hg38_mega" in rel and os.path.isfile(path):
                    out[os.path.join(top, rel)] = _content(path)
        return out

    def advice_rm(self):
        """The removal the advice names (indel_window.enrichment_stores): all of the dataset
        that comes from an enrichment; its VCFs, registry and indel frequencies stay."""
        import indel_window as iw
        for p in iw.enrichment_stores(self.work, "hg38", "hg38_mega"):
            shutil.rmtree(p)
        self.assertTrue(os.path.isdir(self.path("VCFs", "hg38_mega")))
        self.assertTrue(os.path.isdir(self.path("Dictionaries", "registry_hg38_mega")))
        self.assertTrue(os.path.isdir(self.path("Dictionaries", "indel_af_hg38_mega")))

    def refused(self, *named):
        """The recipe run in the working directory: refused, each of ``named`` said, nothing
        made there and no engine called."""
        before, calls = _listing(self.work), len(self.calls())
        done = self.rebuild()
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertEqual([p for p in _listing(self.work) if p not in before], [], done.stderr)
        self.assertEqual(_listing(self.work), before)
        self.assertEqual(len(self.calls()), calls)
        for n in named:
            self.assertIn(n, done.stderr)
        self.assertIn("and nothing was made there: run this recipe in a new working directory",
                      done.stderr)
        self.assertIn("docs/TROUBLESHOOTING.md section 12", done.stderr)

    def test_a_folder_that_holds_the_datasets_vcfs_is_refused_with_nothing_made(self):
        # the data folder as the advice's removal leaves it: the merged VCFs, the registry and
        # the indel frequencies of the dataset (advice_flow.sh of the integration review)
        vdir = self.path("VCFs", "hg38_mega")
        os.makedirs(vdir)
        shutil.copy(os.path.join(self.src, "mega.chr1.afmax.vcf.gz"), vdir)
        for part in ("registry_hg38_mega", "indel_af_hg38_mega"):
            os.makedirs(self.path("Dictionaries", part))
        self.refused("%s already holds VCFs: this recipe writes the dataset's VCFs itself" % vdir,
                     self.path("Dictionaries", "registry_hg38_mega") + " is there",
                     self.path("Dictionaries", "indel_af_hg38_mega") + " is there")
        self.assertFalse(os.path.exists(self.path(".mega_rebuild_hg38_mega")))

    def test_a_folder_that_holds_a_part_it_builds_is_refused_with_nothing_made(self):
        # a kit, a log, a store or an index of the dataset in place: build-index-only would
        # reuse the kit (an earlier one's windows), or a store would be left beside this one's
        for part in self.PARTS:
            with self.subTest(part=os.path.join(*part)):
                p = self.path(*part)
                os.makedirs(p)
                try:
                    self.refused("[error] %s is there: this recipe builds it" % p)
                finally:
                    shutil.rmtree(p)
                    for top in ("Genomes", "Dictionaries", "genome_library"):
                        if os.path.isdir(self.path(top)) and not os.listdir(self.path(top)) \
                                and top != "Genomes":
                            os.rmdir(self.path(top))
        # and with none of it there, it builds as before
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    # the shell section 12's commands and a note's removals are pasted into (the class's,
    # bash; TheRecipeInAFolderOfTheDatasetPastedIntoZsh's, zsh)
    SHELL = "bash"

    def paste(self, command, cwd=None):
        """``command``, as a message gives it, pasted into the class's shell (pasted) from
        ``cwd``, the fixture's root by default: it exits 0, and the shell says no error."""
        rc, out, err = pasted(self, self.SHELL, command, cwd or self.root)
        self.assertEqual((rc, shell_errors(err)), (0, []), (command, out, err))

    def run_section_12(self, src=None, chroms=None, **stub):
        """The commands of docs/TROUBLESHOOTING.md section 12, pasted, as a user pastes them,
        into the class's shell, interactive and with no option of its own (pasted: bash, or
        zsh, which reads a '#' as a word and stops a command at a pattern that matches
        nothing), the chromosomes and the PAM those of this genome (or the chromosomes
        ``chroms``), DATA the working directory, NEW a new folder beside it, SRC the merged
        VCFs (or ``src``), REPO the repository; ``stub``: more of the stub engine's
        environment. Returns the run (returncode, stdout, stderr; returncode None where the
        shell never got past the commands). Up to the review of wave 1c-a this test ran them
        in bash -eu, which a user's shell does not add, and so it missed a removal that ran
        after the recipe had failed (W1C-B2R2-1); up to wave 1c-b of 2.10.0, in bash -c alone,
        and so it missed that zsh, the default shell of macOS, ran them itself, and stopped
        at their first comment."""
        commands = _section_12_commands()
        self.assertIn('bash "$REPO/seq_script/merge_panels/mega_rebuild.sh" "$NEW" "$SRC"',
                      commands)
        commands = (commands.replace("chr{1..22} chrX chrY", chroms or " ".join(bio.GENOME))
                    .replace("20bp-NRR-SpCas9.txt", os.path.basename(bio.PAM)))
        env = {k: v for k, v in os.environ.items() if not k.startswith(("STUB_", "CRISPRME_"))}
        env.update(PATH=os.pathsep.join([self.bin, os.environ.get("PATH", "")]),
                   PYTHONPATH=os.path.join(bio.SHARED, "py"),
                   MPLCONFIGDIR=os.path.join(bio.SHARED, "mpl"), STUB_LOG=self.log,
                   STUB_MARK=self.mark, STUB_RELEASE=self.release, BDNA="1", BRNA="1",
                   THREADS="2", PYTHON=sys.executable, PP=HERE, DATA=self.work,
                   NEW=os.path.join(self.root, "new"), SRC=src or self.src, REPO=REPO, **stub)
        rc, out, err = pasted(self, self.SHELL, commands, self.root, env=env)
        return subprocess.CompletedProcess(self.SHELL, rc, out, err)

    def section_12(self):
        """section 12's commands (run_section_12), and what each says checked: the recipe's
        checks pass, nothing of the dataset is left in the new folder but its samplesID file
        (the data folder has its own, which stays), and no staging folder in the data folder.
        Returns the run."""
        new = os.path.join(self.root, "new")
        ran = self.run_section_12()
        self.assertEqual(ran.returncode, 0, ran.stdout + ran.stderr)
        self.assertEqual(shell_errors(ran.stderr), [], ran.stderr)
        self.assertIn("SUMMARY: 2 chroms, 2 SNP alleles + 2 indel-path alleles, ALL PASS",
                      ran.stdout)
        self.assertEqual(sorted(self.the_kit(new)), ["samplesIDs/hg38_mega.samplesID.txt"])
        self.assertFalse(os.path.exists(self.path(".mega_rebuild_hg38_mega")))
        return ran

    def kit_in_place(self, kit):
        """The data folder holds the kit ``kit`` (the_kit's), file for file and no other part
        of the dataset, and a build of another PAM on it reads it (the recipe's header: another
        PAM is build-index-only on the finished install, which reuses its enrichment), with no
        warning."""
        got = self.the_kit(self.work)
        self.assertEqual(sorted(got), sorted(kit))
        for name in kit:
            self.assertEqual(got[name], kit[name], name)
        nrr = os.path.join(REPO, "PAMs", "20bp-NRR-SpCas9.txt")
        shutil.copy(nrr, self.path("PAMs", os.path.basename(nrr)))
        rc, out, err = self.finish(self.start(command=[
            "build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam",
            self.path("PAMs", os.path.basename(nrr)), "--bDNA", "1", "--bRNA", "1",
            "--thread", "2", "--path", self.work, "--vcf", self.path("VCFs", "hg38_mega")]))
        self.assertEqual((rc, err), (0, ""), out)
        self.assertNotIn("WARNING [build-index-only]", out)
        self.assertTrue(os.path.isdir(self.lib("20-NRR_2_hg38+hg38_mega")))
        self.assertTrue(os.path.isdir(self.lib("20-NRR_2_hg38+hg38_mega_INDELS")))

    def test_the_commands_of_section_12_put_a_rebuilt_kit_in_place(self):
        done = self.rebuild()   # the dataset as an install holds it, the data folder its own
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        kit = self.the_kit(self.work)
        self.assertTrue(any(k.startswith("Genomes/hg38+hg38_mega_INDELS/") for k in kit), kit)
        self.advice_rm()
        self.refused("already holds VCFs")
        self.section_12()
        # the same kit, back in the data folder, which a build of another PAM reads
        self.kit_in_place(kit)

    def test_section_12_in_a_data_folder_that_crisprme_2_8_1_set_up(self):
        """W1C-B2R2-1, the review of wave 1c-a: section 12 copied the recipe's PAM file from
        the data folder's PAMs/, where CRISPRme 2.8.1 put none of NRR (2.10.0's default): run
        there, the copy and then the recipe failed, and the removal that came next took the
        whole dataset (on chr22, 16 entries of hg38_mega down to its samplesID file). The PAM
        file comes from the repository, and the commands run as written in a data folder
        whose PAMs/ lacks it."""
        kit = self.an_earlier_kit_in_place()
        os.remove(self.path("PAMs", os.path.basename(bio.PAM)))
        self.section_12()
        self.kit_in_place(kit)

    def test_section_12_leaves_the_data_folder_as_it_was_when_the_recipe_fails(self):
        """W1C-B2R2-1, the review of wave 1c-a: section 12's removal and moves ran whatever the
        recipe did, and a pasted block has no set -e: a recipe that failed, for a mistyped SRC,
        a FASTA of no chromosome it builds in its genome folder, or its build failing (a full
        disk, an engine that stops), left a data folder with no hg38_mega at all, its merged
        VCFs given as SRC among what was removed. The commands stop at the first command that
        fails now, and the data folder is as it was: every file, folder and link, and each
        file's bytes."""
        self.an_earlier_kit_in_place()
        before = bio.tree(self.work)
        new = os.path.join(self.root, "new")
        for why, kw in (("a mistyped SRC", dict(src=os.path.join(self.root, "merged_typo"))),
                        ("a FASTA of no chromosome it builds",
                         dict(chroms=" ".join(list(bio.GENOME) + ["chr3"]))),
                        ("its build failing", dict(STUB_FAIL="add-variants"))):
            with self.subTest(why=why):
                shutil.rmtree(new, ignore_errors=True)
                ran = self.run_section_12(**kw)
                # (the shell got past the commands, and they stopped with an error of theirs)
                self.assertIsNotNone(ran.returncode, ran.stdout + ran.stderr)
                self.assertNotEqual(ran.returncode, 0, ran.stdout + ran.stderr)
                self.assertNotIn("ALL PASS", ran.stdout)
                self.assertEqual(bio.tree(self.work), before, ran.stdout + ran.stderr)
                if why == "its build failing":   # (the recipe got as far as its build)
                    self.assertTrue(os.path.isdir(os.path.join(new, ".mega_rebuild_hg38_mega")),
                                    ran.stderr)
                shutil.rmtree(new, ignore_errors=True)
        # and with nothing in the way, as written, they put the recipe's kit in place
        self.section_12()

    SACAS9 = os.path.join(REPO, "PAMs", "21bp-NNGRRT-SaCas9.txt")
    EARLIER = "an_earlier_kit.txt"

    def an_earlier_kit_in_place(self):
        """The data folder with the dataset's kit in place, as an install holds it, and what
        an earlier kit leaves beside the parts the recipe builds (a file of its own in each of
        its genomes and its indel logs, and the indexes of another PAM, as 2.8.1's mega has
        NGG's): the kit the recipe builds, which must replace it (the_kit's, as built here)."""
        done = self.rebuild()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        kit = self.the_kit(self.work)
        for part in (("Genomes", "hg38+hg38_mega"), ("Genomes", "hg38+hg38_mega_INDELS"),
                     ("Dictionaries", "log_indels_hg38_mega"),
                     ("genome_library", "NGG_2_hg38+hg38_mega"),
                     ("genome_library", "NGG_2_hg38+hg38_mega_INDELS")):
            bio._write(self.path(*part, self.EARLIER), "an earlier kit's\n")
        return kit

    def a_searchs_note(self, bulges=0):
        """complete-search beside the kit in place, as a user runs it (the real CLI, the stub
        engine, whose windows of 25/34 hold up to 26 bases at an MNV): a SaCas9 site, 21 nt and
        NNGRRT, needs 27 with no bulge (28 with ``bulges`` of each type, 1), so it prints its
        note on hg38_mega, the one line returned. With bulges, the search then builds its
        variant and INDELS indexes of hg38_mega from the kit in place (and stops at the stub's
        search, which it has none of)."""
        guide = os.path.join(self.root, "guide_sa.txt")
        for path, text in ((guide, bio.GENOME["chr1"][:21] + "NNNNNN\n"),
                           (os.path.join(self.root, "vcfs_mega.txt"), "hg38_mega\n"),
                           (os.path.join(self.root, "samples_mega.txt"),
                            "hg38_mega.samplesID.txt\n")):
            bio._write(path, text)
        rc, out, err = self.finish(self.start(command=[
            "complete-search", "--genome", self.path("Genomes", "hg38"),
            "--vcf", os.path.join(self.root, "vcfs_mega.txt"),
            "--samplesID", os.path.join(self.root, "samples_mega.txt"), "--guide", guide,
            "--pam", self.SACAS9, "--mm", "0", "--bDNA", str(bulges), "--bRNA", str(bulges),
            "--output", "note", "--thread", "2"], CRISPRME_SKIP_SCORER_ENV="1"), timeout=300)
        notes = [ln for ln in out.splitlines() if ln.startswith(
            "WARNING [complete-search]: the fake-indel windows of hg38_mega hold 25 reference "
            "bases before each indel and 34 after it")]
        self.assertEqual(len(notes), 1, out + err)
        self.assertIn("this search's targets need up to %d (a 27-nt guide and PAM%s)" % (
            27 + bulges, " with 1 DNA bulge" if bulges else ""), notes[0])
        self.assertIn("docs/TROUBLESHOOTING.md section 12 has the commands", notes[0])
        return notes[0]

    @staticmethod
    def removals(said):
        """The command of each removal a message gives, as written, in its order: '(rm -rf
        ...)', or the kit's, '(rm -rf ...; find ...; ...; and any copy of those indexes ...)',
        its find of each index name pattern among it (indel_window.kit_removal)."""
        return [m.group(1) for m in re.finditer(
            r"\((rm -rf .*?)(?:\)|; and any copy of those indexes)", said)]

    def section_12_after_the_note(self, every):
        """W1C-B1, the review of wave 1c-a: a search's note of windows too short on a
        sites-only dataset named only the removal of the dataset's VCFs, registry and indel
        frequencies (relative to the working directory), and none of the kit in place. Run
        from the state it left, the kit in place, section 12's commands, which removed only
        those three too, failed to move the kit the recipe built in ('Directory not empty'),
        and the data folder was left with the old kit's genomes and indel logs beside the new
        indexes, registry and indel frequencies, which build-index-only and every search then
        refused. The note names the removal of all that came from the kit now, by absolute
        path, and section 12 removes every part of the dataset the recipe builds by name.
        Here the note's first removal alone, the three, as the review's flow ran it, or
        (``every``) every removal it names, each from the data folder; then section 12 as
        written: the kit the recipe builds is in place, the earlier kit's parts and its
        indexes of another PAM gone, and a build reads it."""
        kit = self.an_earlier_kit_in_place()
        note = self.a_searchs_note()
        said = self.removals(note)
        self.assertTrue(said, note)
        for command in (said if every else said[:1]):
            self.paste(command, self.work)
        self.section_12()
        self.kit_in_place(kit)
        for name in ("NGG_2_hg38+hg38_mega", "NGG_2_hg38+hg38_mega_INDELS"):
            self.assertFalse(os.path.lexists(self.lib(name)), name)
        return said

    def test_section_12_after_the_first_removal_a_searchs_note_names(self):
        said = self.section_12_after_the_note(every=False)
        # (the three, by absolute path, and then that of the kit in place: its parts by
        # absolute path, the dataset's indexes by a find of their name patterns,
        # indel_window.index_removals)
        import indel_window as iw
        self.assertEqual(len(said), 2, said)
        self.assertEqual(said[0], "rm -rf " + " ".join("'%s'" % self.path(*p) for p in (
            ("VCFs", "hg38_mega"), ("Dictionaries", "registry_hg38_mega"),
            ("Dictionaries", "indel_af_hg38_mega"))))
        for part in (("Genomes", "hg38+hg38_mega"), ("Genomes", "hg38+hg38_mega_INDELS"),
                     ("Dictionaries", "log_indels_hg38_mega")):
            self.assertIn("'%s'" % self.path(*part), said[1])
        for find in iw.index_removals(self.path("genome_library"), "hg38", "hg38_mega"):
            self.assertIn(find, said[1].split("; "))
        self.assertNotIn("registry_hg38_mega", said[1])

    def test_section_12_after_every_removal_a_searchs_note_names(self):
        self.section_12_after_the_note(every=True)

    def test_the_flow_a_note_gives_leaves_no_index_its_search_built_from_the_old_kit(self):
        """W1C-B2R2-2, the review of wave 1c-a: a search's note of windows too short named the
        removal of the indexes of the dataset in place when it was printed, and the search
        then built its own from the kit in place (a SaCas9 search with bulges beside 2.8.1's
        mega: 21-NNGRRT_2_hg38+hg38_mega and its _INDELS, from the kit of 25/27). The note's
        own flow (its removal of the kit, the recipe in a new working directory, its removal of
        the three, the moves it names) put the recipe's kit in place and left those two
        indexes: the next SaCas9 search, refused for the INDELS index, read that variant index
        of the old kit once its own advice had built the INDELS index again (16 targets lost
        on chr22). The note names the dataset's indexes by their patterns now, which take
        those its search builds after it: its flow, followed as it says, leaves none of them,
        and a build reads the kit with no warning."""
        kit = self.an_earlier_kit_in_place()
        built = ("21-NNGRRT_2_hg38+hg38_mega", "21-NNGRRT_2_hg38+hg38_mega_INDELS")
        for name in built:
            self.assertFalse(os.path.lexists(self.lib(name)), name)
        note = self.a_searchs_note(bulges=1)
        for name in built:   # (built after the note, from the kit in place: the old one)
            self.assertTrue(os.path.isdir(self.lib(name)), name)
            bio._write(os.path.join(self.lib(name), self.EARLIER), "an earlier kit's\n")
        three, gone = self.removals(note)
        # the note's removal of the kit in place and all that was built from it, first
        self.paste(gone)
        # the recipe in a new working directory of its own, whose Genomes/hg38 holds the
        # FASTAs of the chromosomes to build and no other, beside the PAM file
        new = os.path.join(self.root, "new")
        os.makedirs(os.path.join(new, "Genomes", "hg38"))
        for chrom in bio.GENOME:
            os.symlink(self.path("Genomes", "hg38", chrom + ".fa"),
                       os.path.join(new, "Genomes", "hg38", chrom + ".fa"))
        os.makedirs(os.path.join(new, "PAMs"))
        shutil.copy(bio.PAM, os.path.join(new, "PAMs"))
        done = self.rebuild(new)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        # then the note's removal of the dataset's VCFs, registry and indel frequencies, and
        # the moves it names
        self.paste(three)
        self.assertIn("move in Genomes/hg38+hg38_mega and Genomes/hg38+hg38_mega_INDELS, the "
                      "indexes of hg38_mega of its genome_library (not its reference index), "
                      "Dictionaries/log_indels_hg38_mega, registry_hg38_mega and "
                      "indel_af_hg38_mega, and VCFs/hg38_mega", note)
        for part in (("Genomes", "hg38+hg38_mega"), ("Genomes", "hg38+hg38_mega_INDELS"),
                     ("Dictionaries", "log_indels_hg38_mega"),
                     ("Dictionaries", "registry_hg38_mega"),
                     ("Dictionaries", "indel_af_hg38_mega"), ("VCFs", "hg38_mega")):
            os.rename(os.path.join(new, *part), self.path(*part))
        for name in os.listdir(os.path.join(new, "genome_library")):
            if "+hg38_mega" in name:
                os.rename(os.path.join(new, "genome_library", name), self.lib(name))
        # no index the search built from the old kit, nor any other of it, is left
        for name in built + ("NGG_2_hg38+hg38_mega", "NGG_2_hg38+hg38_mega_INDELS"):
            self.assertFalse(os.path.lexists(self.lib(name)), name)
        self.kit_in_place(kit)

    def test_section_12_beside_the_datasets_variant_indexes_alone(self):
        """Section 12 removes the dataset's indexes by a find of each name pattern, which may
        match both kinds, none or one: with an index of each kind in the data folder (an
        earlier kit's, test_section_12_in_a_data_folder_that_crisprme_2_8_1_set_up), with none
        (after the advice's removal, test_the_commands_of_section_12_put_a_rebuilt_kit_in_place)
        and here with the dataset's variant indexes alone, its INDELS indexes removed, the
        commands put the recipe's kit in place, pasted into bash or into zsh. Up to wave 1c-b
        of 2.10.0 they gave the two patterns to rm -rf, which zsh, the default shell of macOS,
        expands itself, and so stops the whole command where one matches nothing."""
        kit = self.an_earlier_kit_in_place()
        gone = [n for n in os.listdir(self.path("genome_library"))
                if n.endswith("_hg38+hg38_mega_INDELS")]
        self.assertEqual(sorted(gone), ["20-NRG_2_hg38+hg38_mega_INDELS",
                                        "NGG_2_hg38+hg38_mega_INDELS"])
        for name in gone:
            shutil.rmtree(self.lib(name))
        self.section_12()
        self.kit_in_place(kit)
        self.assertFalse(os.path.lexists(self.lib("NGG_2_hg38+hg38_mega")))


@unittest.skipUnless(HAVE_PYSAM, "pysam (BGZF and tabix) not installed")
@unittest.skipUnless(shutil.which("zsh"), "zsh is not installed")
class TheRecipeInAFolderOfTheDatasetPastedIntoZsh(TheRecipeInAFolderOfTheDataset):
    """Every test above that pastes section 12's commands or a note's removal, pasted into
    zsh, the default shell of macOS (wave 1c-b of 2.10.0): zsh ran section 12's subshell,
    ( set -eu ... ), itself, and stopped at its first comment, a word in an interactive zsh
    (a command not found), and expanded the index patterns that a note gave rm -rf, stopping
    the whole command where one matched nothing. The two tests of the recipe's refusals paste
    nothing."""

    SHELL = "zsh"
    test_a_folder_that_holds_the_datasets_vcfs_is_refused_with_nothing_made = None
    test_a_folder_that_holds_a_part_it_builds_is_refused_with_nothing_made = None


if __name__ == "__main__":
    unittest.main(verbosity=2)
