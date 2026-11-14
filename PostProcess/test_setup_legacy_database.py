"""`crisprme.py setup` (setup_legacy_database.py) tries the Hugging Face data repo first for
the chromosome's FASTA, the 1000 Genomes and HGDP VCFs and the two sample lists, as
complete-test does, and said 'HuggingFace fetch failed for <file> (404 Client Error ... Entry
Not Found ...); falling back to original source' for each one the repo does not hold: every
VCF, samplesIDs.HGDP.txt, and the unpacked FASTA of any chromosome but chr22. `setup --chrom
chr22` said it three times, then downloaded the three files from their sources as it should.
A file the Hub answers it does not hold is said in one line on stdout, as complete-test says
it; a failure to fetch one (no connection, a checksum that does not match) is still said on
stderr, as before.

The Hub's errors are huggingface_hub's own where it is installed; the light CI env has none,
and stand-ins of their hierarchy take its place for the test.
"""

import contextlib
import io
import os
import shutil
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import setup_legacy_database as sld


class WhatTheDataRepoDoesNotHostIsNoFailure(unittest.TestCase):

    VCF_1000G = "ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    VCF_HGDP = "hgdp_wgs.20190516.full.chr22.vcf.gz"
    # what lucapinello/crisprme-data holds of the files `setup --chrom chr22` asks it for
    HOSTED = ("genomes/hg38/chr22.fa", "samplesIDs/samplesIDs.1000G.txt")

    def setUp(self):
        try:
            from huggingface_hub.utils import EntryNotFoundError, LocalEntryNotFoundError
        except ImportError:
            class EntryNotFoundError(Exception):
                pass

            class LocalEntryNotFoundError(FileNotFoundError, EntryNotFoundError):
                pass

            hub = types.ModuleType("huggingface_hub")
            hub.utils = types.ModuleType("huggingface_hub.utils")
            hub.utils.EntryNotFoundError = EntryNotFoundError
            hub.utils.LocalEntryNotFoundError = LocalEntryNotFoundError
            patcher = mock.patch.dict(sys.modules, {"huggingface_hub": hub,
                                                    "huggingface_hub.utils": hub.utils})
            patcher.start()
            self.addCleanup(patcher.stop)
        # what hf_hub_download raises for a file the repo does not hold (named
        # RemoteEntryNotFoundError from huggingface_hub 1.0 on), and with no connection
        # and the file not cached
        self.missing = type("RemoteEntryNotFoundError", (EntryNotFoundError,), {})(
            "404 Client Error. Entry Not Found for url: https://huggingface.co/datasets/...")
        self.offline = LocalEntryNotFoundError(
            "An error happened while trying to locate the file on the Hub and we cannot find "
            "the requested files in the local cache.")
        self.wd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.wd)
        cwd = os.getcwd()
        os.chdir(self.wd)
        self.addCleanup(os.chdir, cwd)

    def fetch(self, call, hosted=HOSTED, hf_error=None, md5s=()):
        """``call`` with hf_fetch fetching the files in ``hosted`` and raising ``hf_error``
        (the Hub's answer that it holds no such file, by default) for any other, the other
        downloads faked, and compute_md5 answering ``md5s``, then each file's MD5 from the
        tables: its stdout, its stderr and the URLs it downloaded from."""
        tables = dict(sld.MD51000G, **sld.MD5HGDP, **sld.MD5SAMPLES, **sld.MD5ANNOTATION)
        md5s, fetched = list(md5s), []

        def hf_fetch(name, dest):
            if name not in hosted:
                raise hf_error or self.missing
            open(dest, "w").close()

        def download(dest, fname=None, ftp_conn=False, ftp_server=None, ftp_path=None,
                     http_url=None):
            fetched.append(http_url or f"ftp://{ftp_server}{ftp_path}")
            path = os.path.join(dest, os.path.basename(fetched[-1]))
            open(path, "w").close()
            return path

        def gunzip(_gz, out):
            open(out, "w").close()
            return out

        def compute_md5(path):
            return md5s.pop(0) if md5s else tables.get(os.path.basename(path))

        out, err = io.StringIO(), io.StringIO()
        with mock.patch.object(sld, "hf_fetch", hf_fetch), \
             mock.patch.object(sld, "download", download), \
             mock.patch.object(sld, "gunzip", gunzip), \
             mock.patch.object(sld, "untar", lambda *_args: None), \
             mock.patch.object(sld, "compute_md5", compute_md5), \
             contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            call()
        return out.getvalue(), err.getvalue(), fetched

    def test_the_raw_vcfs(self):
        os.makedirs("VCFs")
        out, err, fetched = self.fetch(
            lambda: sld._download_vcf_data("chr22", Path("VCFs"), False))
        self.assertEqual(out, "Downloading genetic variants data for dataset 1000G\n"
                              f"{self.VCF_1000G} is not on HuggingFace ({sld.HF_DATA_REPO}); "
                              "downloading it from its original source\n"
                              "Downloading genetic variants data for dataset HGDP\n"
                              f"{self.VCF_HGDP} is not on HuggingFace ({sld.HF_DATA_REPO}); "
                              "downloading it from its original source\n")
        self.assertEqual(err, "Asset missing, will download: "
                              f"{os.path.join('VCFs', 'hg38_1000G', self.VCF_1000G)}\n"
                              "Asset missing, will download: "
                              f"{os.path.join('VCFs', 'hg38_HGDP', self.VCF_HGDP)}\n")
        self.assertEqual(fetched, ["ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
                                   "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
                                   + self.VCF_1000G,
                                   "ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/"
                                   + self.VCF_HGDP])

    def test_the_hgdp_sample_list(self):
        out, err, fetched = self.fetch(lambda: sld._download_samples_ids_data(Path("."), False))
        self.assertEqual(out, "Downloading samples ids data for dataset 1000G\n"
                              "Downloading samples ids data for dataset HGDP\n"
                              f"samplesIDs.HGDP.txt is not on HuggingFace ({sld.HF_DATA_REPO}); "
                              "downloading it from its original source\n")
        self.assertEqual(err, "Asset missing, will download: "
                              f"{os.path.join('samplesIDs', 'hg38_1000G.samplesID.txt')}\n"
                              "Fetched samplesIDs.1000G.txt from HuggingFace "
                              f"({sld.HF_DATA_REPO})\n"
                              "Asset missing, will download: "
                              f"{os.path.join('samplesIDs', 'hg38_HGDP.samplesID.txt')}\n")
        self.assertEqual(fetched, [f"{sld.TEST_DATA_BASE_URL}/samplesIDs/samplesIDs.HGDP.txt"])
        for ds in ("1000G", "HGDP"):
            self.assertTrue(os.path.isfile(os.path.join("samplesIDs", f"hg38_{ds}.samplesID.txt")))

    def test_the_fasta_of_another_chromosome(self):
        os.makedirs("Genomes")
        out, err, fetched = self.fetch(
            lambda: sld._download_genome_data("chr21", Path("Genomes"), False))
        self.assertEqual(out, "Downloading genome data for chromosomes: chr21\n"
                              f"chr21.fa is not on HuggingFace ({sld.HF_DATA_REPO}); "
                              "downloading it from UCSC\n")
        self.assertEqual(err, "Asset missing, will download: "
                              f"{os.path.join('Genomes', 'hg38_chr21', 'chr21.fa')}\n")
        self.assertEqual(fetched, [f"{sld.HG38_BASE_URL}/chromosomes/chr21.fa.gz"])
        self.assertTrue(os.path.isfile(os.path.join("Genomes", "hg38_chr21", "chr21.fa")))

    def test_a_failure_to_fetch_is_still_said(self):
        with self.subTest(why="a file it holds does not arrive whole"):
            expected = sld.MD5SAMPLES["samplesIDs.1000G.txt"]
            out, err, fetched = self.fetch(
                lambda: sld._download_samples_ids_data(Path("."), False), md5s=("0" * 32,))
            self.assertEqual(out, "Downloading samples ids data for dataset 1000G\n"
                                  "Downloading samples ids data for dataset HGDP\n"
                                  "samplesIDs.HGDP.txt is not on HuggingFace "
                                  f"({sld.HF_DATA_REPO}); downloading it from its original "
                                  "source\n")
            self.assertEqual(err, "Asset missing, will download: "
                                  f"{os.path.join('samplesIDs', 'hg38_1000G.samplesID.txt')}\n"
                                  "HuggingFace fetch failed for samplesIDs.1000G.txt (MD5 "
                                  f"mismatch for samplesIDs.1000G.txt: expected '{expected}', "
                                  f"got '{'0' * 32}'); falling back to original source\n"
                                  "Asset missing, will download: "
                                  f"{os.path.join('samplesIDs', 'hg38_HGDP.samplesID.txt')}\n")
            self.assertEqual(len(fetched), 2)

        def every_download():
            sld._download_genome_data("chr22", Path("Genomes"), False)
            sld._download_vcf_data("chr22", Path("VCFs"), False)
            sld._download_samples_ids_data(Path("."), False)

        with self.subTest(why="no connection"):
            shutil.rmtree("samplesIDs", ignore_errors=True)
            os.makedirs("Genomes")
            os.makedirs("VCFs")
            out, err, fetched = self.fetch(every_download, hosted=(), hf_error=self.offline)
            self.assertNotIn("HuggingFace", out)
            self.assertEqual(
                [ln for ln in err.splitlines() if "HuggingFace" in ln],
                [f"HuggingFace fetch failed for chr22.fa ({self.offline}); falling back to UCSC"]
                + [f"HuggingFace fetch failed for {f} ({self.offline}); falling back to "
                   "original source" for f in (self.VCF_1000G, self.VCF_HGDP,
                                               "samplesIDs.1000G.txt", "samplesIDs.HGDP.txt")])
            self.assertEqual(len(fetched), 5)

    def test_setup_of_chr22_says_the_three_files_and_no_failure(self):
        out, err, fetched = self.fetch(
            lambda: sld.run_legacy_database_setup("chr22", Path(self.wd), False))
        self.assertEqual(
            [ln for ln in out.splitlines() if "HuggingFace" in ln],
            [f"{f} is not on HuggingFace ({sld.HF_DATA_REPO}); downloading it from its original "
             "source" for f in (self.VCF_1000G, self.VCF_HGDP, "samplesIDs.HGDP.txt")])
        self.assertEqual([ln for ln in err.splitlines() if "HuggingFace" in ln],
                         [f"Fetched chr22.fa from HuggingFace ({sld.HF_DATA_REPO})",
                          f"Fetched samplesIDs.1000G.txt from HuggingFace ({sld.HF_DATA_REPO})"])
        self.assertEqual([os.path.basename(url) for url in fetched],
                         [self.VCF_1000G, self.VCF_HGDP, "samplesIDs.HGDP.txt",
                          "gencode.protein_coding.bed.tar.gz",
                          "dhs+encode+gencode.hg38.bed.tar.gz"])


class ThePamFilesSetupWritesAreTheShippedOnes(unittest.TestCase):
    """`crisprme.py setup` writes PAM files into the working directory's PAMs/, over any of the
    same name. Its 22bp-NNGRRN-SaCas9.txt held a 23-nt guide (25 N, then GRRN), while the file
    of that name CRISPRme ships, from which the index set of 2.10.0 builds 22-NNGRRN_3_hg38
    (seq_script/rebuild/index_set_2.10.0.tsv), holds a 22-nt one: in a working directory set up
    so, a search of a 22-nt guide with that PAM file, on the published index, read a 23-nt
    PAM file (found in the sweep of the stage-2b review of lane D). Every PAM file setup writes
    under the name of one CRISPRme ships holds what that one holds."""

    def test_same_name_same_content(self):
        shipped_dir = Path(__file__).resolve().parent.parent / "PAMs"
        with tempfile.TemporaryDirectory() as tmp, \
                contextlib.redirect_stdout(io.StringIO()):
            written = sld._write_all_pam_files(Path(tmp))
            shared = sorted(n for n in written if (shipped_dir / n).is_file())
            self.assertIn("22bp-NNGRRN-SaCas9.txt", shared)
            for name in shared:
                self.assertEqual(Path(written[name]).read_text().split(),
                                 (shipped_dir / name).read_text().split(), name)


if __name__ == "__main__":
    unittest.main()
