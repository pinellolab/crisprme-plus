"""
This module provides functionality to execute the CRISPRme test workflow, including
downloading genomic and variant data, preparing input files, and running the CRISPRme
command-line tool. It includes functions for managing directories, downloading data
from various sources, and configuring test parameters.

Key functions include:
- `_assign_genome_directory_name`: Resolves the genome subdirectory name for the 
    requested chromosome selection.
- `download_genome_data`: Downloads genome data for a specified chromosome to the
    destination directory.
- `ensure_vcf_dataset_directory`: Ensures the existence of a directory for a specific
    VCF dataset.
- `download_vcf_data`: Downloads VCF data for a specific chromosome and variant dataset.
- `ensure_samplesids_directory`: Ensures the existence of the 'samplesIDs' directory.
- `download_samples_ids_data`: Downloads samples IDs data for a specific variant dataset.
- `ensure_annotation_directory`: Ensures the existence of the 'annotation' directory.
- `download_annotation_data`: Downloads gencode and encode annotation data to the
    'annotation' directory.
- `write_ngg_pamfile`: Writes a test PAM file containing the NGG sequence.
- `write_sg1617_guidefile`: Writes a test guide file containing the sg1617 guide
    sequence.
- `write_vcf_config`: Writes a test VCF list file for a specific variant dataset.
- `write_samplesids_config`: Writes a test samples ID list file for a specific
    variant dataset.
- `run_crisprme_test`: Executes the CRISPRme test workflow for a specified chromosome
    and dataset.
- `main`: The entry point of the module that orchestrates the test execution.

This module is designed to facilitate the testing and validation of the CRISPRme
tool, ensuring that all necessary data and configurations are correctly handled
before running the analysis.
"""

from utils import (
    check_crisprme_directory_tree,
    download,
    hf_fetch,
    hf_missing,
    HF_DATA_REPO,
    gunzip,
    untar,
    rename,
    compute_md5,
    CHROMS,
    CRISPRME_DIRS,
    MD5GENOME,
    MD51000G,
    MD5HGDP,
    MD5SAMPLES,
    MD5ANNOTATION,
)

from typing import Optional, Tuple

import index_fit
import unfinished_marks

import subprocess
import contextlib
import tempfile
import shutil
import signal
import json
import time
import sys
import os
import re
from fnmatch import fnmatchcase
from glob import glob

# benchmark registry (test/benchmark/benchmarks.json)
BENCHMARKS_JSON = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, "test", "benchmark",
    "benchmarks.json",
)

# define genome data url
HG38URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"

# define 1000G server and url
VCF1000GSERVER = "ftp.1000genomes.ebi.ac.uk"
VCF1000GURL = (
    "/vol1/ftp/data_collections/1000_genomes_project/release/"
    "20190312_biallelic_SNV_and_INDEL/"
    "ALL.{}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)

# define hgdp server and url
VCFHGDPSERVER = "ngs.sanger.ac.uk"
VCFHGDPURL = "/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{}.vcf.gz"

# define test data url
TESTDATAURL = os.environ.get(
    "CRISPRME_TESTDATA_URL",
    "https://raw.githubusercontent.com/pinellolab/crisprme-plus/refs/heads/main/test/data/",
)

# define complete-test results folder name
COMPLETETESTRESDIR = "crisprme-test-out"

# name prefix of the per-run view a one-chromosome smoke searches from (_chrom_view)
CHROMVIEWPREFIX = ".crisprme-test-view-"

# name prefix of the folder a smoke's PAM file goes in, when PAMs/ has none (_smoke_pam_file)
SMOKEPAMPREFIX = ".crisprme-test-pam-"

# seconds a one-chromosome search has to end on the SIGTERM or SIGHUP passed on to it
# before it is killed (_stop_search)
SEARCHSTOPWAIT = 10


def require_prebuilt_index() -> str:
    """complete-test runs a variant search, which needs a prebuilt index.

    CRISPRme does NOT build an index automatically -- the user decides to DOWNLOAD one
    (or build it). If no usable index (a genome_library/<name>/ dir with .bin shards) is
    present in the working directory, print how to obtain one and stop with exit 2, as the
    re-run guard (check_output) does, since nothing ran; do not download test data or build
    anything.

    Returns:
        The name of the first usable prebuilt index directory found.
    """
    lib = os.path.join(os.getcwd(), "genome_library")
    present = []
    if os.path.isdir(lib):
        for name in sorted(os.listdir(lib)):
            if name.endswith("_INDELS"):
                continue  # the fake-indel companion, not the primary search index
            d = os.path.join(lib, name)
            if os.path.isdir(d) and any(f.endswith(".bin") for f in os.listdir(d)):
                present.append(name)
    if not present:
        sys.stderr.write(
            f"\nNo prebuilt index found for the complete-test (looked in {lib}).\n"
            "CRISPRme does not build one automatically.\n\n"
            "Download an index first, e.g.:\n"
            "  crisprme.py download --what index --index-name NRR_3_hg38+hg38_HGDP1kGP\n"
            "  (the recommended default: SpCas9 NRR, phased HGDP + 1000 Genomes, 4,091\n"
            "   samples; or --index-name NRR_3_hg38 for a reference-only scan)\n\n"
            "then re-run:  crisprme.py complete-test --chrom chr22 --thread 4\n"
            "  (chromosome 22 only; without --chrom it searches every chromosome)\n"
        )
        # exit 2, not 0, which a script or a CI step read as a pass though nothing ran (as
        # check_output's); crisprme.py passes it on
        sys.exit(2)
    # prefer a variant index (name contains "+") so the smoke exercises the variant path,
    # then the common path (_smoke_rank): sorted order alone took NNNNNNNN_3 before NRG_3
    variant = [n for n in present if "+" in n]
    return min(variant or present, key=_smoke_rank)


def _recorded_five_prime(index_dir: str) -> Optional[bool]:
    """Whether an index records a 5' PAM: its name, which says the side since 2.9.0
    (index_fit.parse_index_name: NNNN-20_2_hg38 for a 5' PAM), or else its .pam_build
    ("<name> <seq_len> <pam_pos>", a negative pam_pos), as the search script's _index_side
    reads them. None without either."""
    parsed = index_fit.parse_index_name(os.path.basename(os.path.normpath(index_dir)))
    if parsed and parsed[4]:
        return parsed[4] == 5
    try:
        with open(os.path.join(index_dir, ".pam_build")) as fh:
            return int(fh.readline().split()[-1]) < 0
    except (OSError, ValueError, IndexError):
        return None


def _smoke_rank(index_name: str) -> tuple:
    """Which of several indexes the smoke takes first, sorted order within each: a 3'
    dedicated index of a 3-nt PAM, the common path; an all-N one of 3 nt; another PAM or
    a 5' one; and last an all-N index wider than 3 nt, which takes the whole genome's
    flanks to search. An index of N 1 holds no bulge, and the smoke searches with them."""
    parsed = _parse_index_name(index_name)
    if parsed is None:
        return (2, 0, index_name)
    token, bmax = parsed[0], index_fit.parse_index_name(index_name)[1]
    all_n = set(token) == {"N"}
    if all_n and len(token) > 3:
        rank = 3
    elif len(token) != 3 or _recorded_five_prime(
            os.path.join(os.getcwd(), "genome_library", index_name)):
        rank = 2
    else:
        rank = 1 if all_n else 0
    return (int(bmax < 2), rank, index_name)


def check_output(benchmarks: bool) -> None:
    """
    Refuse to run complete-test over the results of an earlier run of the same mode.

    The benchmark mode (``benchmarks``) writes a crisprme-test-out_<name> folder in
    Results/ for each benchmark, which validate-test reads, and the prebuilt-index smoke
    writes crisprme-test-out. complete-search refuses an output folder that is not empty,
    and validate-test would read an earlier run's results as this run's: every folder of
    this mode there (crisprme-test-out_<name> of any benchmark, a heavy one or one an
    earlier release ran too) is named, all in one message, before anything is downloaded
    or searched. A folder of the other mode is in neither's way, nor is an empty one.

    Args:
        benchmarks (bool): Whether this run is the benchmark mode's.

    Raises:
        SystemExit: With status 2 if any is there.
        AssertionError: If the main results directory does not exist.
    """
    results_dir = os.path.abspath(os.path.join(os.getcwd(), CRISPRME_DIRS[1]))
    assert os.path.isdir(results_dir)
    # (an empty one is no earlier run's results: complete-search makes its output folder,
    # then a refusal of its pre-flight leaves it empty, and that stopped every run after
    # it, until the user removed it, R2-F2; complete-search takes an empty folder)
    found = [os.path.join(results_dir, d) for d in sorted(os.listdir(results_dir))
             if (d.startswith(f"{COMPLETETESTRESDIR}_") if benchmarks
                 else d == COMPLETETESTRESDIR)
             and os.path.isdir(os.path.join(results_dir, d))
             and os.listdir(os.path.join(results_dir, d))]
    if found:
        # all in one message (named one at a time, the benchmark mode's six stopped a user
        # five more times), and exit 2, not 0, which a script or a CI step read as a pass
        # though nothing ran; crisprme.py passes it on
        sys.stderr.write(
            "Complete-test already run once. Please delete the complete-test results "
            f"folder{'s' if len(found) > 1 else ''} before running it again: "
            f"{' '.join(found)}\n"
        )
        sys.exit(2)


def _assign_genome_directory_name(chrom: str) -> str:
    """Return the genome subdirectory name for a chromosome selection.

    The whole-genome build (``chrom == "all"``) lives in ``hg38``; a
    single-chromosome build lives in ``hg38_<chrom>`` so per-chromosome test
    assets stay isolated from a full build.

    Args:
        chrom (str): Chromosome in UCSC format (e.g. "chr22"), or "all".

    Returns:
        str: Genome subdirectory name (e.g. "hg38" or "hg38_chr22").
    """
    return "hg38" if chrom == "all" else f"hg38_{chrom}"


def download_genome_data(chrom: str, dest: str) -> str:
    """
    Download genome data for the requested chromosome selection and return the
    directory containing the resulting FASTA file(s).

    The whole-genome build (``chrom == "all"``) is stored in ``<dest>/hg38``.
    A single-chromosome build is stored in ``<dest>/hg38_<chrom>``.

    Note on coordinates: this function only handles file placement; the FASTA
    payload is UCSC hg38 (1-based, closed genomic coordinates) and is not
    modified here.

    Args:
        chrom (str): Chromosome in UCSC format (e.g. "chr22"), or "all".
        dest (str): The CRISPRme "Genomes" directory in which to create the
            genome directory.

    Returns:
        str: Absolute path to the genome directory to pass to
            ``crisprme.py complete-search --genome``.

    Raises:
        ValueError: If the chromosome is invalid or an MD5 check fails.
        FileExistsError: If the destination directory does not exist.
    """
    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    sys.stderr.write(f"Downloading fasta file for chromosome(s) {chrom}\n")
    genome_dir = os.path.join(dest, _assign_genome_directory_name(chrom))
    if chrom == "all":
        chromstar = download(dest, http_url=f"{HG38URL}/bigZips/hg38.chromFa.tar.gz")
        # check genome md5
        if MD5GENOME[os.path.basename(chromstar)] != compute_md5(chromstar):
            raise ValueError(f"Download for {os.path.basename(chromstar)} failed")
        chromsdir = untar(chromstar, dest, "chroms")  # decompress archive
        # rename extracted chroms dir -> Genomes/hg38
        genome_dir = rename(chromsdir, genome_dir)
        assert os.path.isdir(genome_dir)  # was asserting the pre-rename path
    else:
        os.makedirs(genome_dir, exist_ok=True)  # create Genomes/hg38_{chrom}
        chromfa = os.path.join(genome_dir, f"{chrom}.fa")
        # fast path: fetch the extracted FASTA from HuggingFace; on any failure
        # (repo/file missing, no network, huggingface_hub absent) fall back to
        # the original UCSC source transparently. A file the repo does not hold
        # is no failure: one line on stdout says where it comes from instead
        try:
            hf_fetch(f"genomes/hg38/{chrom}.fa", chromfa)
            sys.stderr.write(f"Fetched {chrom}.fa from HuggingFace ({HF_DATA_REPO})\n")
            assert os.path.isfile(chromfa)
            return os.path.abspath(genome_dir)
        except Exception as e:
            if hf_missing(e):
                print(f"{chrom}.fa is not on HuggingFace ({HF_DATA_REPO}); downloading it "
                      "from UCSC", flush=True)
            else:
                sys.stderr.write(
                    f"HuggingFace fetch failed for {chrom}.fa ({e}); falling back to UCSC\n"
                )
        chromgz = download(dest, http_url=f"{HG38URL}/chromosomes/{chrom}.fa.gz")
        chromfa = gunzip(
            chromgz,
            os.path.join(genome_dir, os.path.splitext(os.path.basename(chromgz))[0]),
        )  # decompress chrom FASTA
        assert os.path.isfile(chromfa)
    return os.path.abspath(genome_dir)


def ensure_vcf_dataset_directory(dest: str, dataset: str) -> str:
    """
    Ensure the existence of a directory for a specific VCF dataset within the
    specified destination directory.

    Args:
        dest (str): The destination directory where the VCF dataset directory
            should be created.
        dataset (str): The name or identifier of the VCF dataset.

    Returns:
        str: The path to the VCF dataset directory.
    """

    vcf_dataset_dir = os.path.join(dest, f"hg38_{dataset}")
    if not os.path.exists(vcf_dataset_dir):
        os.mkdir(vcf_dataset_dir)
    return vcf_dataset_dir


def download_vcf_data(chrom: str, dest: str, dataset: str) -> None:
    """
    Download VCF data for a specific chromosome and variant dataset to the destination
    directory.

    Args:
        chrom (str): The chromosome identifier in UCSC format.
        dest (str): The destination directory to save the downloaded VCF data.
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP", or 1000G+HGDP).

    Returns:
        None

    Raises:
        ValueError: If the input chromosome or variant dataset is invalid.
        FileExistsError: If the destination directory does not exist.
    """

    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # create VCF dataset directory within VCFs folder
    sys.stderr.write(f"Downloading VCF data for chromsome(s) {chrom}\n")
    for ds in dataset.split("+"):
        vcf_dataset_dir = ensure_vcf_dataset_directory(dest, ds)
        ftp_server = VCF1000GSERVER if ds == "1000G" else VCFHGDPSERVER
        vcf_url = VCF1000GURL if ds == "1000G" else VCFHGDPURL
        chroms = CHROMS if chrom == "all" else [chrom]
        md5data = MD51000G if ds == "1000G" else MD5HGDP
        for c in chroms:
            expected_fname = os.path.basename(vcf_url.format(c))
            dest_path = os.path.join(vcf_dataset_dir, expected_fname)
            # Resume: skip files already downloaded and checksum-verified in a
            # previous (interrupted) run, so a partially completed multi-hour
            # download can be restarted cheaply instead of starting over.
            if os.path.isfile(dest_path) and md5data.get(
                expected_fname
            ) == compute_md5(dest_path):
                sys.stderr.write(
                    f"{expected_fname} already present and verified; skipping\n"
                )
                continue
            # fast path: fetch the bgzipped VCF from HuggingFace and MD5-verify it;
            # on any failure fall through to the original FTP/HTTPS source below. The
            # default repo hosts none of these raw VCFs (its vcfs/ holds the inputs of
            # the shipped panels), which is no failure: one line on stdout says so
            try:
                hf_fetch(f"vcfs/{ds}/{expected_fname}", dest_path)
                if md5data.get(expected_fname) == compute_md5(dest_path):
                    sys.stderr.write(
                        f"Fetched {expected_fname} from HuggingFace ({HF_DATA_REPO})\n"
                    )
                    continue
                raise ValueError("MD5 mismatch after HuggingFace fetch")
            except Exception as e:
                if hf_missing(e):
                    print(f"{expected_fname} is not on HuggingFace ({HF_DATA_REPO}); "
                          "downloading it from its original source", flush=True)
                else:
                    sys.stderr.write(
                        f"HuggingFace fetch failed for {expected_fname} ({e}); "
                        "falling back to original source\n"
                    )
            # Retry on checksum mismatch (e.g. a truncated transfer from a slow or
            # flaky host) instead of aborting the whole run on the first failure.
            attempts = 3
            for attempt in range(1, attempts + 1):
                if ds == "1000G":
                    # the EBI 1000G host also serves HTTPS; prefer it, since FTP is
                    # frequently blocked on CI runners and institutional/cloud networks
                    vcf = download(
                        vcf_dataset_dir,
                        http_url=f"https://{ftp_server}{vcf_url.format(c)}",
                    )
                else:  # HGDP (Sanger) via FTP
                    vcf = download(
                        vcf_dataset_dir,
                        ftp_conn=True,
                        ftp_server=ftp_server,
                        ftp_path=vcf_url.format(c),
                    )
                if md5data[os.path.basename(vcf)] == compute_md5(vcf):
                    break  # download verified
                sys.stderr.write(
                    f"Checksum mismatch for {os.path.basename(vcf)} "
                    f"(attempt {attempt}/{attempts}); re-downloading\n"
                )
                if attempt == attempts:
                    raise ValueError(
                        f"Download for {os.path.basename(vcf)} failed after "
                        f"{attempts} attempts"
                    )


def ensure_samplesids_directory(dest: str) -> str:
    """
    Ensure the existence of the 'samplesIDs' directory within the specified
    destination directory.

    Args:
        dest (str): The destination directory where the 'samplesIDs' directory
            should be created.

    Returns:
        str: The path to the 'samplesIDs' directory.
    """

    samplesids_dir = os.path.join(dest, CRISPRME_DIRS[6])
    if not os.path.exists(samplesids_dir):
        os.mkdir(samplesids_dir)
    return samplesids_dir


def download_samples_ids_data(dataset: str) -> None:
    """
    Download samples IDs data for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g.,
            "1000G", "HGDP").

    Returns:
        None

    Raises:
        ValueError: If the variant dataset is unknown.
    """

    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # samples ids folder must be located within current directory
    # -- see check_crisprme_directory_tree() for details
    sys.stderr.write(f"Downloading sample ids for dataset(s) {dataset}\n")
    samplesids_dir = ensure_samplesids_directory(os.getcwd())
    for ds in dataset.split("+"):
        samplesid_fname = (
            "samplesIDs.1000G.txt" if ds == "1000G" else "samplesIDs.HGDP.txt"
        )
        # Fetch to the name the config + pipeline expect (write_samplesids_config
        # references hg38_<ds>.samplesID.txt). BOTH the HuggingFace fast path and
        # the original-source fallback must land the file at this final name.
        # Previously only the fallback renamed, so once the file existed on HF the
        # success path left it as "samplesIDs.<ds>.txt" and the run failed with
        # "hg38_<ds>.samplesID.txt: No such file or directory".
        final_path = os.path.join(samplesids_dir, f"hg38_{ds}.samplesID.txt")
        # fast path: HuggingFace, then fall back to the original source (a file the
        # repo does not hold is no failure: one line on stdout says so)
        try:
            hf_fetch(f"samplesIDs/{samplesid_fname}", final_path)
            if MD5SAMPLES.get(samplesid_fname) == compute_md5(final_path):
                sys.stderr.write(
                    f"Fetched {samplesid_fname} from HuggingFace ({HF_DATA_REPO}) "
                    f"-> {os.path.basename(final_path)}\n"
                )
                continue
            raise ValueError("MD5 mismatch after HuggingFace fetch")
        except Exception as e:
            if hf_missing(e):
                print(f"{samplesid_fname} is not on HuggingFace ({HF_DATA_REPO}); "
                      "downloading it from its original source", flush=True)
            else:
                sys.stderr.write(
                    f"HuggingFace fetch failed for {samplesid_fname} ({e}); "
                    "falling back to original source\n"
                )
        samplesids = download(
            samplesids_dir,
            http_url=f"{TESTDATAURL.rstrip('/')}/samplesIDs/{samplesid_fname}",
        )
        if MD5SAMPLES[os.path.basename(samplesids)] != compute_md5(samplesids):
            raise ValueError(f"Download for {os.path.basename(samplesids)} failed")
        rename(samplesids, final_path)

def ensure_annotation_directory(dest: str) -> str:
    """
    Ensure the existence of the 'annotation' directory within the specified
    destination directory.

    Args:
        dest (str): The destination directory where the 'annotation' directory
            should be created.

    Returns:
        str: The path to the 'annotation' directory.
    """

    annotation_dir = os.path.join(dest, CRISPRME_DIRS[4])
    if not os.path.exists(annotation_dir):
        os.mkdir(annotation_dir)
    return annotation_dir


def download_annotation_data() -> Tuple[str, str]:
    """
    Download gencode and encode annotation data to the 'annotation' directory
    within the current working directory.

    Returns:
        Tuple[str, str]: Paths to the downloaded gencode and encode annotation
            files.
    """

    sys.stderr.write("Downloading ENCODE and GENCODE annotation data\n")
    annotation_dir = ensure_annotation_directory(os.getcwd())
    gencode = _retrieve_ann_data(
        annotation_dir,
        "Annotations/gencode.protein_coding.bed.tar.gz",
        "gencode.protein_coding.bed",
    )
    encode = _retrieve_ann_data(
        annotation_dir,
        "Annotations/dhs+encode+gencode.hg38.bed.tar.gz",
        "dhs+encode+gencode.hg38.bed",
    )
    return gencode, encode


def _retrieve_ann_data(annotation_dir: str, url: str, fname: str) -> str:
    """
    Download and extract an annotation file, then compress it with bgzip.

    This function downloads an annotation archive, verifies its integrity, extracts
    the specified file, and compresses it using bgzip. It returns the path to the
    compressed annotation file.

    Args:
        annotation_dir (str): The directory to store the annotation data.
        url (str): The URL of the annotation archive to download.
        fname (str): The name of the file to extract and compress.

    Returns:
        str: The path to the compressed annotation file.

    Raises:
        ValueError: If the downloaded file fails the MD5 check.
    """

    # download gencode annotation
    annfile_tar = download(annotation_dir, http_url=os.path.join(TESTDATAURL, url))
    if MD5ANNOTATION[os.path.basename(annfile_tar)] != compute_md5(annfile_tar):
        raise ValueError(f"Download for {os.path.basename(annfile_tar)} failed")
    return os.path.join(untar(annfile_tar, annotation_dir), fname)


def ensure_pams_directory(dest: str) -> str:
    """
    Ensure the existence of the 'PAMs' directory within the specified destination
    directory.

    Args:
        dest (str): The destination directory where the 'PAMs' directory should
            be created.

    Returns:
        str: The path to the 'PAMs' directory.
    """

    pams_dir = os.path.join(dest, CRISPRME_DIRS[5])
    if not os.path.exists(pams_dir):
        os.mkdir(pams_dir)
    return pams_dir


def write_ngg_pamfile() -> str:
    """
    Write a test PAM file containing the NGG sequence to the 'PAMs' directory
    within the current working directory.

    Returns:
        str: The path to the created test PAM file.
    """

    sys.stderr.write("Creating PAM file\n")
    pams_dir = ensure_pams_directory(
        os.getcwd()
    )  # PAMs directory must be in current working dir
    pamfile = os.path.join(pams_dir, "20bp-NGG-SpCas9.txt")
    try:
        with open(pamfile, mode="w") as outfile:
            outfile.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")  # 20 + 3 bp (NGG)
    except IOError as e:
        raise IOError("An error occurred while writing the test PAM file") from e
    return pamfile


def write_sg1617_guidefile(guide_seq: str = "CTAACAGTTGCTTTTATCACNNN") -> str:
    """
    Write a test guide file containing the sg1617 guide sequence, padded for the PAM it
    is searched with.

    Returns:
        str: The path to the created test guide file.
    """

    sys.stderr.write("Creating guide file\n")
    guidefile = "sg1617_test_guide.txt"
    try:
        with open(guidefile, mode="w") as outfile:
            outfile.write(guide_seq + "\n")  # sg1617 guide
    except IOError as e:
        raise IOError("An error occerred while writing the test guide file") from e
    return guidefile


def load_benchmarks() -> dict:
    """Load the benchmark registry (falls back to the canonical Cas9 case)."""
    try:
        with open(BENCHMARKS_JSON) as fin:
            return json.load(fin)
    except (OSError, ValueError):
        return {
            "thresholds": {"mm": 4, "bDNA": 1, "bRNA": 1},
            "benchmarks": [{
                "name": "cas9_sg1617", "nuclease": "SpCas9",
                "pam_name": "20bp-NGG-SpCas9.txt",
                "pam_content": "NNNNNNNNNNNNNNNNNNNNNGG 3",
                "guide_file": "sg1617_test_guide.txt",
                "guide_crisprme": "CTAACAGTTGCTTTTATCACNNN",
            }],
        }


def write_pamfile(pam_name: str, pam_content: str) -> str:
    """Write a benchmark PAM file into the working-dir PAMs directory."""
    sys.stderr.write(f"Creating PAM file {pam_name}\n")
    pamfile = os.path.join(ensure_pams_directory(os.getcwd()), pam_name)
    with open(pamfile, mode="w") as outfile:
        outfile.write(pam_content.rstrip("\n") + "\n")
    return pamfile


def write_guidefile(guide_file: str, guide_seq: str) -> str:
    """Write a benchmark guide file into the working directory."""
    sys.stderr.write(f"Creating guide file {guide_file}\n")
    with open(guide_file, mode="w") as outfile:
        outfile.write(guide_seq.rstrip("\n") + "\n")
    return guide_file


def write_vcf_config(dataset: str) -> str:
    """
    Write a test VCF list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test VCF list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # config vcf list file
    sys.stderr.write(f"Creating VCF config file for dataset(s) {dataset}\n")
    vcf_config = "vcf.config.test.txt"
    try:
        with open(vcf_config, mode="w") as outfile:
            for ds in dataset.split("+"):
                outfile.write(f"hg38_{ds}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return vcf_config


def write_samplesids_config(dataset: str) -> str:
    """
    Write a test samples ID list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test samples ID list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # configure sample ids list
    sys.stderr.write(f"Creating samples config file for dataset(s) {dataset}\n")
    samples_config = "samplesIDs.config.test.txt"
    fname_map = {"1000G": "hg38_1000G.samplesID.txt", "HGDP": "hg38_HGDP.samplesID.txt"}
    try:
        with open(samples_config, mode="w") as outfile:
            for ds in dataset.split("+"):
                outfile.write(f"{fname_map[ds]}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return samples_config


def _parse_index_name(index_name: str):
    """Parse a CRISPRme index directory name ``<PAM>_<bMax>_<ref>[+<variant>]``, or the
    form that carries the guide length and the side since 2.9.0
    (index_fit.parse_index_name).

    Returns (pam_token, ref_genome, is_variant) or None if it does not parse.
    e.g. ``NRG_3_hg38`` -> ('NRG', 'hg38', False);
         ``NRG_3_hg38+hg38_1000G2021`` -> ('NRG', 'hg38', True);
         ``20-NRR_3_hg38`` -> ('NRR', 'hg38', False).
    """
    parsed = index_fit.parse_index_name(index_name)
    if not parsed:
        return None
    pam_token, _bmax, rest = parsed[:3]
    ref = rest.split("+")[0]
    return pam_token, ref, ("+" in rest)


def _natural_order(name: str) -> list:
    """Sort key that puts chr2 before chr10."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


def _index_chroms(index_dir: str) -> set:
    """Chromosomes an index holds: each is one or more ``<PAM>_<chrom>_<n>.bin`` trees."""
    found = set()
    for f in os.listdir(index_dir):
        tree = re.fullmatch(r"[^_]+_(.+)_\d+\.bin", f)
        if tree:
            found.add(tree.group(1))
    return found


def _require_chrom(chrom: str, genome_dir: str, ref: str, index_name: str) -> None:
    """Refuse a --chrom that Genomes/<ref> or the smoke's index does not have, in one
    line naming the chromosomes both have, before anything is written or searched.
    Exits 2, which crisprme.py passes on without a traceback."""
    in_genome = {f[:-3] for f in os.listdir(genome_dir) if f.endswith(".fa")}
    in_index = _index_chroms(os.path.join(os.getcwd(), "genome_library", index_name))
    lacking = [
        what
        for what, chroms in (
            (f"Genomes/{ref}", in_genome),
            (f"the index '{index_name}'", in_index),
        )
        if chrom not in chroms
    ]
    if not lacking:
        return
    both = sorted(in_genome & in_index, key=_natural_order)
    # hg38 has 455 sequences: past a screenful, name the primary ones
    shown = both if len(both) <= 30 else [c for c in both if "_" not in c][:30] or both[:30]
    more = f" and {len(both) - len(shown)} more" if len(both) > len(shown) else ""
    sys.stderr.write(
        f"ERROR: complete-test --chrom {chrom}: {' and '.join(lacking)} "
        f"{'have' if len(lacking) > 1 else 'has'} no {chrom}; available: "
        f"{', '.join(shown) + more if both else 'none, the two share no chromosome'}\n"
    )
    sys.exit(2)


def _exit_on_signal(signum, _frame) -> None:
    """SIGTERM / SIGHUP: leave through _search_in_view, which stops the search first, and
    the finally that removes the per-run view."""
    sys.exit(128 + signum)


def _proc_stat(pid: int):
    """(state, parent PID, start time) of process ``pid`` from /proc/<pid>/stat (Linux),
    or None once it is gone. Its name, in parentheses, may hold spaces and parentheses,
    so the fields are counted from the last ')'."""
    try:
        with open(f"/proc/{pid}/stat") as fin:
            fields = fin.read().rsplit(")", 1)[1].split()
        return fields[0], int(fields[1]), fields[19]
    except (OSError, IndexError, ValueError):
        return None


def _stop_search(proc, signum: int) -> None:
    """Pass the SIGTERM or SIGHUP ``signum`` complete_test.py got on to the search
    ``proc`` and every process under it, wait for all of them to end, and SIGKILL those
    still running after SEARCHSTOPWAIT seconds.

    A ``kill`` of complete_test.py's PID reaches it alone: the shell under it,
    complete-search, the search script and its children, which share its process
    group, would go on, and fail on the view's paths once it is removed. They are found
    in /proc before any of them ends, as the children of one that ends go to init;
    without /proc (macOS) only ``proc`` is sent the signal. A signal sent to the whole
    group (a hangup, ``kill -- -<pgid>``, ``timeout``) has reached them already, and
    this second one changes nothing: they end on the first. A further SIGTERM or SIGHUP
    (a hangup can bring two) does not cut the wait short; _chrom_view's finally puts
    the handlers back.
    """
    for s in (signal.SIGTERM, signal.SIGHUP):
        if signal.getsignal(s) is _exit_on_signal:
            signal.signal(s, signal.SIG_IGN)
    found = {}  # PID -> start time, of each process found under proc

    def running(pid):
        stat = _proc_stat(pid)
        return stat is not None and stat[0] != "Z" and stat[2] == found[pid]

    def send(sig):
        # find the processes under proc, or under one found before that still runs (its
        # own parent may have ended), then send sig to those still running
        children = {}
        with contextlib.suppress(OSError):
            for name in os.listdir("/proc"):
                stat = _proc_stat(int(name)) if name.isdigit() else None
                if stat is not None:
                    children.setdefault(stat[1], []).append((int(name), stat[2]))
        todo = [pid for pid in found if running(pid)]
        if proc.returncode is None:  # not reaped yet, so its PID is still its own
            todo.append(proc.pid)
        while todo:
            for pid, start in children.get(todo.pop(), ()):
                if pid not in found:
                    found[pid] = start
                    todo.append(pid)
        proc.send_signal(sig)
        for pid in found:
            if running(pid):
                with contextlib.suppress(OSError):
                    os.kill(pid, sig)

    send(signum)
    deadline = time.monotonic() + SEARCHSTOPWAIT
    while time.monotonic() < deadline and (
        proc.poll() is None or any(running(pid) for pid in found)
    ):
        time.sleep(0.05)
    send(signal.SIGKILL)
    proc.wait()


def _search_in_view(cmd: str, view: str) -> int:
    """Run the one-chromosome search, the shell command ``cmd``, from ``view`` and return
    its exit code, as subprocess.call does, but on a SIGTERM or SIGHUP stop the search
    (_stop_search) before the view is removed. On a Ctrl-C, which the terminal sends
    the search too, kill the shell and leave at once, as subprocess.call does: where
    /bin/sh execs complete_test.py, crisprme.py kills it 0.25 s after a Ctrl-C."""
    with subprocess.Popen(cmd, shell=True, cwd=view) as proc:
        try:
            return proc.wait()
        except SystemExit as stop:  # _exit_on_signal's: 128 + the signal's number
            _stop_search(proc, stop.code - 128)
            raise
        except BaseException:
            proc.kill()
            raise


def _log_chrom(fname: str):
    """Chromosome of an indel log ``log<chrom>.txt[.gz]``, chr-prefixed as
    pool_search_indels.py reads it; None for any other file."""
    stem = fname[:-3] if fname.endswith(".gz") else fname
    if not (stem.startswith("log") and stem.endswith(".txt")):
        return None
    chrom = stem[len("log"):-len(".txt")]
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _vcf_names(vcf_list: Optional[str]) -> list:
    """The dataset folder names (<vcf>) of the smoke's --vcf list, as the search reads
    them, each once (a list may name one twice, or as rv, rv/ and VCFs/rv: the view
    gets one folder for it, and a refusal one line); none for a reference smoke (no
    list)."""
    if not vcf_list:
        return []
    with open(vcf_list) as fin:
        return list(dict.fromkeys(os.path.basename(l.strip().rstrip("/")) for l in fin
                                  if l.strip()))


def _refuse_unfinished(ref: str, vcf_names: list, pam: str, prune: bool = False,
                       three_prime: bool = True) -> None:
    """Refuse the one-chromosome smoke while build-index-only has marked unfinished a
    part of a dataset it searches (unfinished_marks.dataset_parts, the parts
    complete-search refuses), before its view is made: the view holds no
    Genomes/<ref>+<vcf>, so the search run from it would not see that one's mark, and
    it would name each part it refuses by its path in the view, gone once the run ends.
    One line on stderr for each, then exit 1, with nothing searched and no results
    folder made. The smoke searches with PAM ``pam``, one bulge of each type (an index
    N of 2), and --prune-pam only on a wider all-N index (``prune``)."""
    cwd = os.getcwd()
    lib = os.path.join(cwd, "genome_library")
    refused = []
    for vcf in vcf_names:
        in_place = unfinished_marks.indels_index_in_place(lib, f"{ref}+{vcf}", pam, 2, prune,
                                                          three_prime)
        for part, what, rerun in unfinished_marks.dataset_parts(cwd, ref, vcf, in_place):
            if unfinished_marks.is_marked(part):
                refused.append(
                    f"ERROR: {part} is unfinished: a build-index-only run is writing {what} "
                    "there, or stopped before it was done, and the search would use what is "
                    f"there as if it were whole. Let that run finish, or {rerun[0].lower()}"
                    f"{rerun[1:]}, then run complete-test again.\n"
                )
    if refused:
        sys.stderr.write("".join(refused))
        sys.exit(1)


@contextlib.contextmanager
def _chrom_view(chrom: str, ref: str, vcf_list: Optional[str]):
    """Yield a per-run view of the working directory in which ``chrom`` is the only
    chromosome, for complete-search to run from, and remove it when done. ``vcf_list``
    is the smoke's --vcf list of dataset folder names (<vcf>), None for a reference one.

    A hidden folder in the working directory that mirrors it with symlinks:
    Genomes/<ref>/ holds only <chrom>.fa; every index the search could resolve
    (genome_library/*_*_<ref> and *_*_<ref>+<vcf>, locked ones aside, as its
    _resolve_index_dir has them) only its <PAM>_<chrom>_<n>.bin trees, and its _INDELS
    companion only the <PAM>_<N>_fake<chrom> sub-index; Genomes/<ref>+<vcf>_INDELS/,
    from which the search builds an INDELS index it does not find, only
    fake_<vcf>_<chrom>/ (and the folder's unfinished mark, if it has one);
    Dictionaries/log_indels_<vcf>/, from which the indel search and post-analysis take
    their chromosomes, only <chrom>'s log (an empty one if the panel has none: no indel
    to search there).
    Everything else in the working directory is linked as it is, but VCFs/, which
    stays empty: a search of a prebuilt index never reads the source VCFs, which would
    list the chromosomes in the logs' place. The index sidecars (.pam_build,
    .display_label) are left out: only the web listings read them, and in their absence
    the search reads these index folders as they are, rather than through a .binview in
    the results whose links would outlive the view. But for the .pam_build of an all-N
    index built for a 5' PAM, which the search needs: without it it takes that index for
    a 3' build, and refuses it. _run_prebuilt_index_smoke removes that .binview.

    The search's results are the whole-genome search's on <chrom> but for one file: for
    a chromosome with an indel log and no indel off-target, the whole-genome run, whose
    indel post-analysis runs for the other chromosomes' indel off-targets, writes a
    header-only fake<chrom>_*.indel_snp_cooc.tsv, which a run with no indel off-target
    at all skips with its post-analysis. It holds no row, and the report leaves
    header-only ones out.

    Nothing is copied or written outside the view, and no listing or later run looks
    for an index in it (they read genome_library/ and Genomes/ only). It is removed on
    success, failure, Ctrl-C, SIGTERM and SIGHUP, the last two passed on to the search
    and waited for first (_search_in_view); a signal ignored at entry, as nohup ignores
    SIGHUP, stays ignored, here and in the search. A SIGKILL leaves the view behind,
    which is harmless: a folder of symlinks that nothing reads. A SIGINT sent to
    complete_test.py alone, not the Ctrl-C a terminal sends the whole process group,
    removes it but leaves the search running, to fail without it.
    """
    cwd = os.getcwd()
    vcf_names = _vcf_names(vcf_list)
    # a signal ignored at entry stays so, and the search inherits it, for exec keeps an
    # ignored signal and resets a handled one to its default: nohup ignores SIGHUP
    saved = {
        s: signal.signal(s, _exit_on_signal)
        for s in (signal.SIGTERM, signal.SIGHUP)
        if signal.getsignal(s) is not signal.SIG_IGN
    }
    view = None
    try:
        view = tempfile.mkdtemp(prefix=CHROMVIEWPREFIX, dir=cwd)

        def link(*rel):
            os.symlink(os.path.join(cwd, *rel), os.path.join(view, *rel))

        genomes, dicts = CRISPRME_DIRS[0], CRISPRME_DIRS[2]
        restricted = (genomes, "genome_library", dicts, CRISPRME_DIRS[3])
        for name in os.listdir(cwd):
            if name not in restricted and not name.startswith((CHROMVIEWPREFIX, SMOKEPAMPREFIX)):
                link(name)
        for name in restricted:
            os.mkdir(os.path.join(view, name))
        os.mkdir(os.path.join(view, genomes, ref))
        link(genomes, ref, f"{chrom}.fa")
        lib = os.path.join(cwd, "genome_library")
        tree = re.compile(r"[^_]+_" + re.escape(chrom) + r"_\d+\.bin")
        suffixes = [ref] + [f"{ref}+{vcf}" for vcf in vcf_names]
        for name in os.listdir(lib):
            idx = os.path.join(lib, name)
            if (
                name.startswith(".")
                or not os.path.isdir(idx)
                or os.path.isdir(f"{idx}.lock")
                or not any(fnmatchcase(name, f"*_*_{s}") for s in suffixes)
            ):
                continue
            os.mkdir(os.path.join(view, "genome_library", name))
            five_prime_all_n = (set(index_fit.index_pam(name)) == {"N"}
                                and _recorded_five_prime(idx))
            for f in os.listdir(idx):
                if tree.fullmatch(f) or (f == ".pam_build" and five_prime_all_n):
                    link("genome_library", name, f)
            indels = f"{name}_INDELS"
            if os.path.isdir(os.path.join(lib, indels)) and not os.path.isdir(
                os.path.join(lib, f"{indels}.lock")
            ):
                os.mkdir(os.path.join(view, "genome_library", indels))
                # <PAM>_<N>_fake<chrom>: the engine's name, whichever form the index's takes
                _parts = index_fit.parse_index_name(name)
                fake = f"{_parts[0]}_{_parts[1]}_fake{chrom}"
                if os.path.isdir(os.path.join(lib, indels, fake)):
                    link("genome_library", indels, fake)
        # A variant index with no _INDELS companion (an indel-free dataset, for which
        # build-index-only builds none; a build stopped before its INDELS index): the
        # search builds it from Genomes/<ref>+<vcf>_INDELS, or finds no indel to index
        # there. It gets <chrom>'s fake-indel genome alone, and the folder's mark if it
        # has one, which the search then refuses as it does outside.
        for vcf in vcf_names:
            fakes = f"{ref}+{vcf}_INDELS"
            if not os.path.isdir(os.path.join(cwd, genomes, fakes)):
                continue
            os.mkdir(os.path.join(view, genomes, fakes))
            for f in (f"fake_{vcf}_{chrom}", ".building"):
                if os.path.exists(os.path.join(cwd, genomes, fakes, f)):
                    link(genomes, fakes, f)
        logdirs = {f"log_indels_{vcf}" for vcf in vcf_names}
        for name in os.listdir(os.path.join(cwd, dicts)):
            if name not in logdirs:
                link(dicts, name)
                continue
            os.mkdir(os.path.join(view, dicts, name))
            logs = [
                f
                for f in os.listdir(os.path.join(cwd, dicts, name))
                if _log_chrom(f) == chrom
            ]
            for f in logs:
                link(dicts, name, f)
            if not logs:
                open(os.path.join(view, dicts, name, f"log{chrom}.txt"), "w").close()
        yield view
    finally:
        if view is not None:
            shutil.rmtree(view, ignore_errors=True)
        for s, handler in saved.items():
            signal.signal(s, handler)


def _zero_hit(res_dir: str) -> bool:
    """Whether the search into ``res_dir`` ended at complete-search's zero-hit
    short-circuit (submit_job): exit 0 with no off-target, so no integrated_results
    table, and its log.txt marked as the web status page reads it."""
    try:
        with open(os.path.join(res_dir, "log.txt")) as fin:
            return "No off-targets found" in fin.read()
    except OSError:
        return False


# the smoke's example guide (sg1617), without its PAM's padding
_EXAMPLE_SPACER = "CTAACAGTTGCTTTTATCAC"


def _pam_files(pams: str) -> list:
    """The names in the PAMs folder ``pams``, sorted, none where it cannot be listed."""
    try:
        return sorted(os.listdir(pams))
    except OSError:
        return []


def _smoke_five_prime(index_name: str, token: str) -> Optional[bool]:
    """Whether the index was built for a PAM on the 5' side of the guide: what its name or
    its .pam_build records (_recorded_five_prime). Without either, an all-N index is a 3'
    build, as the search takes it (_index_side in submit_job_automated_new_multiple_vcfs.sh)
    and every published one is; and a dedicated one, which a search builds with no record,
    is what the PAM files of its motif in PAMs/ say, or else 3' for a PAM of up to 3 nt, as
    every shipped one that short is. None for one wider that no PAM file tells of: its side
    is not known."""
    recorded = _recorded_five_prime(os.path.join(os.getcwd(), "genome_library", index_name))
    if recorded is not None:
        return recorded
    if set(token) == {"N"}:
        return False
    sides = set()
    pams = os.path.join(os.getcwd(), CRISPRME_DIRS[5])
    for name in _pam_files(pams):
        try:
            with open(os.path.join(pams, name)) as fh:
                seq, pos = fh.readline().split()[:2]
            pos = int(pos)
        except (OSError, ValueError, UnicodeDecodeError):
            continue
        if (seq[:-pos] if pos < 0 else seq[-pos:]) == token:
            sides.add(pos < 0)
    if sides:
        return True in sides  # (an earlier complete-test left a 3' file of a 5' motif)
    return False if len(token) <= 3 else None


def _smoke_pam_file(line: str, name: str) -> str:
    """A PAM file that holds ``line``: the one in PAMs/ that does, or else one written
    under the name ``name`` in a hidden folder of the working directory, for the run
    (_smoke_pam_dir): a file written into PAMs/ stayed, and the web offered it as a PAM to
    search with, for a nuclease it is not."""
    pams = os.path.join(os.getcwd(), CRISPRME_DIRS[5])
    for found in _pam_files(pams):
        try:
            with open(os.path.join(pams, found)) as fh:
                if fh.readline().split() == line.split():
                    return os.path.join(pams, found)
        except (OSError, UnicodeDecodeError):
            continue
    folder = tempfile.mkdtemp(prefix=SMOKEPAMPREFIX, dir=os.getcwd())
    _smoke_pam_dir.append(folder)
    path = os.path.join(folder, name)
    with open(path, "w") as fh:
        fh.write(line + "\n")
    return path


# the hidden folders _smoke_pam_file wrote into, to remove when the smoke ends
_smoke_pam_dir = []


def _smoke_inputs(index_name: str, token: str):
    """What the smoke searches an index with, in the index's geometry: (PAM file, guide
    file, the PAM it searches, the search's extra flag, whether that PAM is on the 5' side).

    An all-N index wider than 3 nt on the 3' side (the published NNNNNNNN_3) is searched
    for NRR, the default PAM, with --prune-pam, its documented use, which scores the
    example. Any other is
    searched for its own PAM, on its own side: a 20-nt spacer (fewer on an index built for
    a shorter guide), which any longer build serves, with the example guide padded for
    that PAM. Exits 2, in one line, for an index whose side cannot be told (a search of
    the wrong side finds nothing at all, which would pass under --chrom)."""
    five_prime = _smoke_five_prime(index_name, token)
    if five_prime is None:
        sys.stderr.write(
            f"ERROR: complete-test cannot tell which side of the guide the PAM of index "
            f"'{index_name}' is on: it records no .pam_build, and PAMs/ holds no PAM file of "
            f"{token}. Put its PAM file in PAMs/, or run complete-search on the index with it.\n"
        )
        sys.exit(2)
    if not five_prime and set(token) == {"N"} and len(token) > 3:
        pam = _smoke_pam_file("N" * 20 + "NRR 3", "20bp-NRR-SpCas9.txt")
        return pam, write_sg1617_guidefile(), "NRR", " --prune-pam", False
    built = index_fit.built_guide_len(os.path.join(os.getcwd(), "genome_library", index_name))
    spacer = min(20, built or 20)
    width = len(token)
    if five_prime:
        line = f"{token}{'N' * spacer} -{width}"
        seq = "N" * width + _EXAMPLE_SPACER[:spacer]
    else:
        line = f"{'N' * spacer}{token} {width}"
        seq = _EXAMPLE_SPACER[20 - spacer:] + "N" * width
    pam = _smoke_pam_file(line, f"{spacer}bp-{token}-smoke.txt")
    return pam, write_sg1617_guidefile(seq), token, "", five_prime


def _run_prebuilt_index_smoke(threads: int, debug: bool, chrom: str = "all") -> None:
    """Default (user) complete-test: run the built-in example guide against an index the
    user has DOWNLOADED (or built). CRISPRme never builds an index automatically here; if
    none is present, ``require_prebuilt_index`` prints how to download one and exits.

    Args:
        threads (int): Number of threads for the search.
        debug (bool): Run complete-search in debug mode.
        chrom (str): "all" searches every chromosome of the index's genome; any other
            value only that chromosome, from a per-run view (``_chrom_view``), and one
            that holds no off-target of the example guide passes with an empty result.
    """
    index_name = require_prebuilt_index()
    parsed = _parse_index_name(index_name)
    if parsed is None:
        sys.stderr.write(
            f"Found index '{index_name}' but could not parse it "
            "(expected <PAM>_<bMax>_<ref>[+<variant>]). Aborting.\n"
        )
        sys.exit(1)
    pam_token, ref, is_variant = parsed
    genome_dir = os.path.join(os.getcwd(), CRISPRME_DIRS[0], ref)  # Genomes/<ref>
    if not os.path.isdir(genome_dir):
        sys.stderr.write(
            f"\nThe prebuilt index '{index_name}' is present, but its reference genome "
            f"(Genomes/{ref}) is not.\nCRISPRme does not build data automatically -- "
            "download the genome, e.g.:\n"
            f"  crisprme.py download --what genome --ref {ref}   (or --what all)\n\n"
            "then re-run:  crisprme.py complete-test --chrom chr22 --thread 4\n"
            "  (chromosome 22 only; without --chrom it searches every chromosome)\n"
        )
        sys.exit(2)  # nothing ran (require_prebuilt_index)
    if chrom != "all":
        _require_chrom(chrom, genome_dir, ref, index_name)
    # the built-in example guide, in the index's geometry (_smoke_inputs)
    pam, guide, search_pam, prune_arg, five_prime = _smoke_inputs(index_name, pam_token)
    try:
        _search_in_smoke(index_name, pam_token, ref, is_variant, genome_dir, pam, guide,
                         search_pam, prune_arg, five_prime, threads, debug, chrom)
    finally:
        for folder in _smoke_pam_dir:
            shutil.rmtree(folder, ignore_errors=True)
        del _smoke_pam_dir[:]


def _search_in_smoke(index_name, pam_token, ref, is_variant, genome_dir, pam, guide, search_pam,
                     prune_arg, five_prime, threads, debug, chrom) -> None:
    """The smoke's complete-search, and what it left, for _run_prebuilt_index_smoke."""
    # A variant index searches its variants, in complete-search's context-aware default
    # mode: per-sample on a phased genotyped index such as NRR_3_hg38+hg38_HGDP1kGP,
    # population-level otherwise.
    # `download --what index` writes list_vcf.txt / list_samplesID.txt; a build/DATA dir
    # may instead use vcf_list.txt / samples_list.txt -- accept either.
    vcf_arg = sid_arg = ""
    vcf_list = None
    if is_variant:
        vcf_list = next((f for f in ("list_vcf.txt", "vcf_list.txt") if os.path.isfile(f)), None)
        if vcf_list:
            vcf_arg = f"--vcf {vcf_list}"
            sid_list = next((f for f in ("list_samplesID.txt", "samples_list.txt") if os.path.isfile(f)), None)
            if sid_list:
                sid_arg = f"--samplesID {sid_list}"
        else:
            sys.stderr.write(
                f"NOTE: variant index '{index_name}' present but no list_vcf.txt found; "
                "running a REFERENCE-only smoke. (A `download --what index` creates it.)\n"
            )
    debug_arg = "--debug" if debug else ""
    output_dir = COMPLETETESTRESDIR
    sys.stderr.write(
        f"Running complete-test smoke against prebuilt index '{index_name}' "
        f"(PAM {search_pam}{prune_arg}, genome {ref}, "
        f"{'variant' if (is_variant and vcf_arg) else 'reference'}"
        f"{'' if chrom == 'all' else f', {chrom} only'}); "
        f"mm=2 bDNA=1 bRNA=1 -> Results/{output_dir}\n"
    )

    def complete_search(genome: str, output: str, view=None) -> int:
        # complete-search will NOT build an index.
        cmd = (
            f"crisprme.py complete-search --genome {genome} "
            f"--mm 2 --bDNA 1 --bRNA 1 --max-total-edits 4 "
            f"--merge 3 --pam {pam} --guide {guide} {vcf_arg} {sid_arg} "
            f"--index-path genome_library{prune_arg} --output {output} "
            f"--thread {threads} {debug_arg} --ci-cd-test"
        )
        if view is None:
            return subprocess.call(cmd, shell=True)
        return _search_in_view(cmd, view)

    if chrom == "all":
        returncode = complete_search(genome_dir, output_dir)
    else:
        # the same search, from a view of the working directory where chrom is alone;
        # the --output path (the usual Results folder) given whole, so that the paths
        # complete-search prints and logs are not the view's, gone once it is removed.
        # The parts of its datasets it refuses while they are marked unfinished are
        # looked at here first, as the view does not hold all of them.
        _refuse_unfinished(ref, _vcf_names(vcf_list), search_pam, bool(prune_arg), not five_prime)
        with _chrom_view(chrom, ref, vcf_list) as view:
            returncode = complete_search(
                os.path.join(view, CRISPRME_DIRS[0], ref),
                os.path.join(os.getcwd(), CRISPRME_DIRS[1], output_dir),
                view,
            )
        # (the search's view of an index with a sidecar, whose links went with the view)
        shutil.rmtree(os.path.join(os.getcwd(), CRISPRME_DIRS[1], output_dir, ".binview"),
                      ignore_errors=True)
    if returncode != 0:
        sys.stderr.write(
            f"ERROR: complete-test search failed (exit code {returncode}). See the log above.\n"
        )
        sys.exit(returncode)
    # confirm a report / integrated_results was produced
    res_dir = os.path.join(os.getcwd(), CRISPRME_DIRS[1], output_dir)
    produced = glob(os.path.join(res_dir, "*integrated_results.tsv"))
    if not produced and chrom != "all" and _zero_hit(res_dir):
        # the example guide has no off-target on chrom (chrM, most unplaced contigs), so
        # the whole-genome test reports none there either: an empty result
        sys.stderr.write(
            f"\ncomplete-test PASSED: searched the example guide against '{index_name}' "
            f"on {chrom}, which holds none of its off-targets (an empty result, so no "
            "integrated_results.tsv).\n"
        )
        return
    if not produced and _zero_hit(res_dir) and (five_prime or (len(search_pam) > 3 and not prune_arg)):
        # an SpCas9 guide, searched with another nuclease's PAM, can have no off-target
        # in a whole genome, on its index's own geometry: an empty result too
        sys.stderr.write(
            f"\ncomplete-test PASSED: searched the example guide against '{index_name}', "
            f"which holds none of its off-targets with PAM {search_pam} (an empty result, so "
            "no integrated_results.tsv).\n"
        )
        return
    if not produced:
        sys.stderr.write(
            f"ERROR: complete-test produced no integrated_results.tsv in {res_dir}.\n"
        )
        sys.exit(1)
    sys.stderr.write(
        f"\ncomplete-test PASSED: searched the example guide against '{index_name}'"
        f"{'' if chrom == 'all' else f' on {chrom}'} and "
        f"produced {os.path.basename(produced[0])} in Results/{output_dir}.\n"
    )


def _run_benchmarks_ondemand(chrom: str, dataset: str, threads: int, debug: bool) -> None:
    """CI / brute-force-gate path (enabled by CRISPRME_ALLOW_ONDEMAND_BUILD): download the
    test genome + VCF + samples + annotation for ``chrom``/``dataset``, then run EVERY
    registered benchmark via complete-search, building each benchmark's index on demand
    (permitted only because the caller set CRISPRME_ALLOW_ONDEMAND_BUILD). Produces the
    per-benchmark output dirs (``crisprme-test-out_<name>``) that ``validate-test`` compares
    against the committed brute-force ground truth. This is NOT the default user path -- users
    get ``_run_prebuilt_index_smoke`` (never auto-builds).
    """
    genome_dir = download_genome_data(chrom, CRISPRME_DIRS[0])  # download genome data
    download_vcf_data(chrom, CRISPRME_DIRS[3], dataset)  # download vcf data
    vcf = write_vcf_config(dataset)  # write test vcf list
    download_samples_ids_data(dataset)  # download vcf dataset samples ids
    samplesids = write_samplesids_config(dataset)  # write test samples ids list
    gencode, encode = download_annotation_data()  # gencode + encode annotation
    debug_arg = "--debug" if debug else ""
    # Run one complete-search per registered benchmark. complete-search refuses to run into
    # a non-empty output folder, so each benchmark gets its OWN output dir; validate-test
    # looks in each.
    registry = load_benchmarks()
    global_th = registry.get("thresholds", {"mm": 4, "bDNA": 1, "bRNA": 1})
    # "heavy" benchmarks build a bulge-2 variant index. They run by default: the whole set
    # fits a 16 GiB, 4-CPU hosted runner. CRISPRME_SKIP_HEAVY=1 leaves them out on a very
    # small machine.
    skip_heavy = os.environ.get("CRISPRME_SKIP_HEAVY") == "1"
    for bench in registry["benchmarks"]:
        if skip_heavy and bench.get("heavy"):
            sys.stderr.write(
                f"Skipping heavy benchmark '{bench['name']}' (CRISPRME_SKIP_HEAVY=1; "
                "unset it to run this bulge-2 benchmark too)\n"
            )
            continue
        th = dict(global_th)
        th.update(bench.get("thresholds", {}))
        # bMax (= bDNA+bRNA) is NOT passed: complete-search has never had a --bmax flag,
        # it recomputes bMax from the per-type budgets. It used to be passed anyway and was
        # silently ignored until the unknown-flag check started rejecting it. The binding
        # knob is --max-total-edits: default mm+bDNA+bRNA (non-binding -> per-type mode)
        # unless the case pins a smaller single-n value. The brute-force ground truth must
        # be generated with the SAME per-type budgets AND --max-total-edits (see
        # generate_references.py).
        max_total_edits = th.get("max_total_edits", th["mm"] + th["bDNA"] + th["bRNA"])
        output_dir = f"{COMPLETETESTRESDIR}_{bench['name']}"
        pam = write_pamfile(bench["pam_name"], bench["pam_content"])
        guide = write_guidefile(bench["guide_file"], bench["guide_crisprme"])
        sys.stderr.write(
            f"Running complete-search for benchmark '{bench['name']}' "
            f"({bench.get('nuclease', '')}) mm={th['mm']} bDNA={th['bDNA']} "
            f"bRNA={th['bRNA']} max-total-edits={max_total_edits} -> {output_dir}\n"
        )
        # NO --index-path: complete-search builds this benchmark's index on demand (allowed
        # because run_crisprme_test gated this branch on CRISPRME_ALLOW_ONDEMAND_BUILD).
        crisprme_cmd = (
            f"crisprme.py complete-search --genome {genome_dir} "
            f"--mm {th['mm']} --bDNA {th['bDNA']} --bRNA {th['bRNA']} "
            f"--max-total-edits {max_total_edits} "
            f"--merge 3 --pam {pam} --guide {guide} --vcf {vcf} "
            f"--samplesID {samplesids} --annotation {encode} "
            f"--gene_annotation {gencode} --output {output_dir} "
            f"--thread {threads} {debug_arg} --ci-cd-test"
        )
        returncode = subprocess.call(crisprme_cmd, shell=True)
        if returncode != 0:
            sys.stderr.write(
                "ERROR: complete-test failed during complete-search for benchmark "
                f"'{bench['name']}' (exit code {returncode}). See the log above.\n"
            )
            sys.exit(returncode)


def run_crisprme_test(chrom: str, dataset: str, threads: int, debug: bool) -> None:
    """complete-test entry point, with two modes:

    * DEFAULT (user): run the built-in example guide against an index the user DOWNLOADED
      (or built), on ``chrom`` alone, or on every chromosome for "all". CRISPRme never
      builds an index automatically here.
    * CI / brute-force gate (env ``CRISPRME_ALLOW_ONDEMAND_BUILD`` set): download the test
      data for ``chrom``/``dataset`` and run every registered benchmark with on-demand index
      build, so ``validate-test`` can compare predictions to the committed brute-force ground
      truth. The ``validate-benchmarks`` CI sets this flag.
    """
    check_crisprme_directory_tree(os.getcwd())  # check crisprme directory tree
    benchmarks = bool(os.environ.get("CRISPRME_ALLOW_ONDEMAND_BUILD"))
    check_output(benchmarks)  # refuse to overwrite an earlier run's results
    if benchmarks:
        _run_benchmarks_ondemand(chrom, dataset, threads, debug)
    else:
        _run_prebuilt_index_smoke(threads, debug, chrom)


def main():
    chrom, dataset, threads, debug = sys.argv[1:]  # read commandline args
    debug = debug == "True"
    run_crisprme_test(chrom, dataset, int(threads), debug)  # run crisprme test


if __name__ == "__main__":
    main()
