#!/usr/bin/env python

from typing import List, NoReturn, Tuple
from Bio.Seq import Seq

import subprocess
import itertools
import pathlib
import shutil
import sys
import os
import re


version = "2.5.1-dev"  # CRISPRme version; drop -dev when tagging v2.5.1
__version__ = version

script_path = os.path.dirname(os.path.abspath(__file__))
origin_path = os.path.dirname(os.path.abspath(__file__))
# path where this file is located
# origin_path = os.path.dirname(os.path.realpath(__file__))
# conda path
conda_path = "opt/crisprme/PostProcess/"
# path corrected to use with conda
corrected_origin_path = script_path[:-3] + conda_path
corrected_web_path = f"{origin_path[:-3]}/opt/crisprme/"
# corrected_web_path = os.getcwd()

script_path = corrected_origin_path
current_working_directory = f"{os.getcwd()}/"
# script_path = corrected_web_path+"/PostProcess/"

input_args = sys.argv

if "--debug" in input_args:
    print("DEBUG MODE")
    script_path = current_working_directory + "PostProcess/"
    corrected_web_path = current_working_directory

# Load the non-fatal low-memory warning from the PostProcess directory by
# explicit file path: crisprme.py is installed in bin/, so PostProcess is not
# importable as a package, and crisprme.py otherwise shells out rather than
# importing PostProcess modules. Never let this optional check break startup.
import importlib.util as _ilu

try:
    _mc_spec = _ilu.spec_from_file_location(
        "memory_check", os.path.join(script_path, "memory_check.py")
    )
    _mc_mod = _ilu.module_from_spec(_mc_spec)
    _mc_spec.loader.exec_module(_mc_mod)
    warn_low_memory = _mc_mod.warn_low_memory
except Exception:  # pragma: no cover - defensive; memory check is optional
    def warn_low_memory(*args, **kwargs):
        return None

sys.path.insert(0, script_path)
from validate_inputs import run_lightweight, run_full, resolve_vcf_dataset_dirs  # noqa: E402
from crisprme_hf import (  # noqa: E402  (huggingface_hub imported lazily inside)
    download_component,
    publish_index,
    resolve_repo,
    synthesize_combined_samplesid,
    DEFAULT_HF_REPO,
)
from utils import download_reference_genome  # noqa: E402
from assembly_reconcile import reconcile_haplotypes, check_liftover_available, haplotype_search_complete, clean_incomplete_haplotype_output, haplotype_params_match  # noqa: E402

cicd_test = False
if "--ci-cd-test" in input_args:
    cicd_test = True

VALID_CHARS = {
    "a",
    "A",
    "t",
    "T",
    "c",
    "C",
    "g",
    "G",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "r",
    "y",
    "s",
    "w",
    "k",
    "m",
    "b",
    "d",
    "h",
    "v",
}

CRISPRMEDIRS = [
    "Genomes", "Results", "Dictionaries", "VCFs", "Annotations", "PAMs", "samplesIDs",
]

def is_folder_empty(folder: str) -> bool:
    return any(os.scandir(folder))


# Input chr1:11,130,540-11,130,751
def extractSequence(name, input_range, genome_selected):
    name = "_".join(name.split())
    current_working_directory = os.getcwd() + "/"
    chrom = input_range.split(":")[0]
    start_position = (
        input_range.split(":")[1]
        .split("-")[0]
        .replace(",", "")
        .replace(".", "")
        .replace(" ", "")
    )
    end_position = (
        input_range.split(":")[1]
        .split("-")[1]
        .replace(",", "")
        .replace(".", "")
        .replace(" ", "")
    )

    list_chr = [
        f
        for f in os.listdir(current_working_directory + "Genomes/" + genome_selected)
        if os.path.isfile(
            os.path.join(current_working_directory + "Genomes/" + genome_selected, f)
        )
        and not f.endswith(".fai")
    ]
    add_ext = ".fa"
    if ".fasta" in list_chr[0]:
        add_ext = ".fasta"
    with open(current_working_directory + name + ".bed", "w") as b:
        b.write(chrom + "\t" + start_position + "\t" + end_position)

    output_extract = subprocess.check_output(
        [
            "bedtools getfasta -fi "
            + current_working_directory
            + "Genomes/"
            + genome_selected
            + "/"
            + chrom
            + add_ext
            + " -bed "
            + current_working_directory
            + name
            + ".bed"
        ],
        shell=True,
    ).decode("utf-8")
    try:
        os.remove(
            current_working_directory
            + "Genomes/"
            + genome_selected
            + "/"
            + chrom
            + ".fa.fai"
        )
    except:
        pass
    try:
        os.remove(current_working_directory + name + ".bed")
    except:
        pass
    ret_string = output_extract.split("\n")[1].strip()
    return ret_string


def getGuides(extracted_seq, pam, len_guide, pam_begin):
    len_pam = len(pam)
    # dict
    len_guide = int(len_guide)
    pam_dict = {
        "A": "ARWMDHV",
        "C": "CYSMBHV",
        "G": "GRSKBDV",
        "T": "TYWKBDH",
        "R": "ARWMDHVSKBG",
        "Y": "CYSMBHVWKDT",
        "S": "CYSMBHVKDRG",
        "W": "ARWMDHVYKBT",
        "K": "GRSKBDVYWHT",
        "M": "ARWMDHVYSBC",
        "B": "CYSMBHVRKDGWT",
        "D": "ARWMDHVSKBGYT",
        "H": "ARWMDHVYSBCKT",
        "V": "ARWMDHVYSBCKG",
        "N": "ACGTRYSWKMBDHV",
    }
    list_prod = []
    for char in pam:
        list_prod.append(pam_dict[char])

    iupac_pam = []  # NNNNNNN NGG
    for element in itertools.product(*list_prod):
        iupac_pam.append("".join(element))

    rev_pam = str(Seq(pam).reverse_complement())
    list_prod = []
    for char in rev_pam:
        list_prod.append(pam_dict[char])

    # CCN NNNNNNN  -> results found with this pam must be reverse complemented
    iupac_pam_reverse = []
    for element in itertools.product(*list_prod):
        iupac_pam_reverse.append("".join(element))

    extracted_seq = extracted_seq.upper()
    len_sequence = len(extracted_seq)
    guides = []
    for pam in iupac_pam:
        pos = [m.start() for m in re.finditer("(?=" + pam + ")", extracted_seq)]
        if pos:
            for i in pos:
                if pam_begin:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    guides.append(extracted_seq[i + len_pam : i + len_pam + len_guide])
                else:
                    if i < len_guide:
                        continue
                    # guides.append(extracted_seq[i-len_guide:i+len_pam])           # i is position where first char of pam is found, eg the N char in NNNNNN NGG
                    # print('1 for:' , extracted_seq[i-len_guide:i])
                    guides.append(extracted_seq[i - len_guide : i])
    for pam in iupac_pam_reverse:  # Negative strand
        pos = [m.start() for m in re.finditer("(?=" + pam + ")", extracted_seq)]
        if pos:
            for i in pos:
                if pam_begin:
                    if i < len_guide:
                        continue
                    guides.append(
                        str(Seq(extracted_seq[i - len_guide : i]).reverse_complement())
                    )
                else:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    # guides.append(str(Seq(extracted_seq[i:i+len_pam+len_guide]).reverse_complement()))         # i is position where first char of pam is found, eg the first C char in CCN NNNNNN
                    # print('2 for:', str(Seq(extracted_seq[i + len_pam : i + len_guide + len_pam]).reverse_complement()))
                    guides.append(
                        str(
                            Seq(
                                extracted_seq[i + len_pam : i + len_guide + len_pam]
                            ).reverse_complement()
                        )
                    )
    return guides
    # return guides for when adding to app.py


def check_crisprme_dirtree() -> None:
    """Ensures that the working directory contains the required CRISPRme folder 
    structure.
    
    Checks for the existence of each expected directory and creates any that are 
    missing.
    """
    # check that working folder respect crisprme's directory tree structure
    for directory in CRISPRMEDIRS:
        if not os.path.exists(os.path.join(current_working_directory, directory)):
            # expected folder not found, create it
            crisprmedir = os.path.join(current_working_directory, directory)
            os.makedirs(crisprmedir)

def print_help_complete_search() -> None:
    """Prints detailed help information for the complete-search functionality.

    Outputs a description of the pipeline and lists all available command-line 
    options to stderr, then exits the program.
    """
    # functionality description
    sys.stderr.write(
        "The complete-search functionality is an end-to-end automated pipeline "
        "that takes raw input files and performs the full workflow up to "
        "post-analysis. Starting from the user-provided genome, variants, guides, "
        "PAM and annotation files, it identifies potential CRISPR off-targets "
        "incorporating variant and haplotype information, and scores each "
        "candidate guide. The pipeline performs genome-wide searches, integrates "
        "annotation data, and generates comprehensive reports."
    )
    # options
    sys.stderr.write(
        "Options:\n"
        "\t--genome, specify the reference genome folder [REQUIRED]\n"
        "\t--vcf, specify a file listing VCF folders (one per line) [OPTIONAL]\n"
        "\t--guide, specify a file containing guide RNAs [REQUIRED if --sequence "
        "not provided]\n"
        "\t--sequence, specify a file with DNA sequences or BED coordinates to "
        "extract guides [REQUIRED if --guide not provided]\n"
        "\t--pam, specify a file containing the PAM sequence [REQUIRED]\n"
        "\t--be-window, specify the window to search for base editor "
        "susceptibility (e.g., --be-window 4,8) [OPTIONAL]\n"
        "\t--be-base, the base(s) for the chosen base editor (e.g., --be-base "
        "A,C) [OPTIONAL]\n"
        "\t--annotation, specify BED files with genome annotations (e.g., "
        "regulatory elements, enhancers). The fourth column must contain the "
        "annotation name. The input BED files must be compressed using bgzip "
        "[OPTIONAL]\n"
        "\t--personal_annotation, specify BED files with personal genomic "
        "annotations. The fourth column must contain the annotation name. The "
        "input BED files must be compressed using bgzip [OPTIONAL]\n"
        "\t--samplesID, specify a file listing sample files (one per line) "
        "present in samplesIDs folder [OPTIONAL]\n"
        "\t--gene_annotation, specify gene annotation (e.g., GENCODE) to find "
        "nearest gene for each target (must be bgzip-compressed) [OPTIONAL]\n"
        "\t--mm, number of mismatches allowed in the search [REQUIRED]\n"
        "\t--bDNA, number of DNA bulges allowed in the search [OPTIONAL]\n"
        "\t--bRNA, number of RNA bulges allowed in the search [OPTIONAL]\n"
        "\t--merge, window size (nucleotides) to merge candidate off-targets "
        "using the highest scoring as pivot [default: 3]\n"
        "\t--sorting-criteria-scoring, comma-separated list to sort targets by "
        "scoring criteria: 'mm', 'bulges', or 'mm+bulges' [default: 'mm+bulges']\n"
        "\t--sorting-criteria, comma-separated list to sort targets by 'mm', "
        "'bulges', or 'mm+bulges' [default: 'mm+bulges,mm']\n"
        "\t--output, specify the output folder name; results will be saved in "
        "Results/<name> [REQUIRED]\n"
        "\t--thread, set number of threads to use [default: 8]\n"
        "\t--index-path, use a prebuilt reference-index library at this path "
        "(e.g. one made with 'build-index-only' or downloaded ahead of time) "
        "instead of building the index under the working directory; a missing "
        "matching index is a hard error [OPTIONAL]\n"
        "\t--max-total-edits, cap the TOTAL edits (mismatches + bulges) per "
        "reported alignment, pruned INSIDE the TST search so excess alignments "
        "are never generated (much faster + smaller intermediates). E.g. with "
        "--mm 6 --bDNA 2 --bRNA 2 --max-total-edits 6, a 4mm+1+1 alignment is "
        "kept but a 6mm+2+2 (=10) one is skipped. Default 4; set it >= "
        "mm+bDNA+bRNA to effectively disable. NOTE: the cap is on the alignment "
        "against the searched (possibly variant-enriched) genome; a variant that "
        "matches the guide lowers the searched edit count, so a VARIANT off-target's "
        "REF-based mismatches+bulges reported in the results can exceed the cap "
        "(reference off-targets are always <= the cap) [OPTIONAL]\n"
        "\t--full_input_validate, also run a full per-VCF-record scan (chromosome "
        "coverage, AF/FILTER consistency, POS bounds, multiallelic/breakend/"
        "duplicate/phasing survey) before launching the search; slower than the "
        "default lightweight checks, so opt-in [OPTIONAL]\n")
    sys.exit(1)


def error(msg: str) -> NoReturn:
    """Prints an error message to stderr and exits the program.

    This function is used to display error messages and terminate execution with 
    a non-zero status.
    
    Args:
        msg: The error message to display.
    """
    sys.stderr.write(f"Error: {msg}\n")
    sys.exit(1)


def summarize_pipeline_failure(outputfolder: str) -> str:
    """Build a concise, actionable summary of a failed pipeline run.

    The search pipeline writes per-stage Start/End markers to ``log.txt`` and
    all subprocess stderr to ``log_error.txt``. On failure the user otherwise
    just sees "run failed"; this reads those two files to report WHICH stage
    failed (the last stage that started but never ended) and WHAT the actual
    error was (the tail of the error log), so a non-expert doesn't have to open
    and interpret the log files by hand.
    """
    lines = []
    log_txt = os.path.join(outputfolder, "log.txt")
    log_err = os.path.join(outputfolder, "log_error.txt")
    # the failing stage = the last one that started but has no matching End
    failing_stage = None
    try:
        started, ended = [], set()
        with open(log_txt) as fh:
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[1] == "Start":
                    started.append(parts[0])
                elif len(parts) >= 2 and parts[1] == "End":
                    ended.add(parts[0])
        for stage in reversed(started):
            if stage not in ended:
                failing_stage = stage
                break
    except OSError:
        pass
    if failing_stage:
        lines.append(f"  Failed during stage: {failing_stage}")
    try:
        with open(log_err) as fh:
            errlines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        if errlines:
            lines.append("  Last lines of the error log:")
            lines.extend(f"      {ln}" for ln in errlines[-15:])
    except OSError:
        pass
    lines.append(f"  Full error log: {log_err}")
    return "\n".join(lines)


def _check_mandatory_args(args: List[str]) -> None:
    if "--genome" not in args:
        error("--genome is required")
    if "--guide" not in args and "--sequence" not in args:
        error("No input guide. One between --guide and --sequence must be specified")
    if "--pam" not in args:
        error("--pam is required")
    if "--mm" not in args:
        error("--mm is required")
    if "--output" not in args:
        error("--output is required")

def _check_genome(args: List[str]) -> str:
    """Retrieves and validates the reference genome directory from command-line 
    arguments.

    Ensures the --genome argument is provided and points to an existing directory. 
    Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the reference genome directory.

    Raises:
        SystemExit: If the --genome argument is missing or the specified directory 
            does not exist.
    """
    try:  # read genome input genome folder path
        genomedir = os.path.abspath(args[args.index("--genome") + 1])
    except IndexError:  # no argument for --genome
        error("Missing input for --genome. Reference genome folder must be specified")
    if not os.path.isdir(genomedir):
        error("The folder specified for --genome does not exist")
    return genomedir

def _check_vcf(args: List[str], variant: bool) -> str:
    """Retrieves and validates the VCF configuration file from command-line 
    arguments.

    Ensures the --vcf argument is provided and points to an existing file if 
    variant-aware search is enabled. Raises an error if the argument is missing 
    or invalid.

    Args:
        args: List of command-line arguments.
        variant: Boolean indicating if variant-aware search is enabled.

    Returns:
        The absolute path to the VCF configuration file, or "_" if variant is False.

    Raises:
        SystemExit: If the --vcf argument is missing or the specified file does 
            not exist.
    """
    if not variant:
        return "_"
    try:
        vcfdir = os.path.realpath(args[args.index("--vcf") + 1])
    except IndexError:
        error("Missing input for --vcf. VCF config file must be specified")
    if not os.path.isfile(vcfdir):
        error("The config file specified for --vcf does not exist")
    return vcfdir

def _check_guide(args: List[str], guide: bool) -> str:
    """Retrieves and validates the guide file from command-line arguments.

    Ensures the --guide argument is provided and points to an existing file, and 
    checks for conflicting input flags. Raises an error if the argument is missing,
    invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the guide file.

    Raises:
        SystemExit: If the --guide argument is missing, the specified file does 
            not exist, or there is a conflict with --sequence.
    """
    if not guide:
        return ""
    if "--guide" in args and "--sequence" in args:
        error("Error: Conflicting flags --guide and --sequence. Use only one")
    try:
        guidefile = os.path.abspath(args[args.index("--guide") + 1])
    except IndexError:
        error("Missing input for --guide. Guide file file must be specified")
    if not os.path.isfile(guidefile):
        error("The file specified for --guide does not exist")
    return guidefile

def _check_sequence(args: List[str], sequence: bool) -> str:
    """Retrieves and validates the sequence file from command-line arguments.

    Ensures the --sequence argument is provided and points to an existing file, 
    and checks for conflicting input flags. Raises an error if the argument is 
    missing, invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the sequence file.

    Raises:
        SystemExit: If the --sequence argument is missing, the specified file does 
            not exist, or there is a conflict with --guide.
    """
    if not sequence:
        return ""
    if "--guide" in args and "--sequence" in args:
        error("Error: Conflicting flags --guide and --sequence. Use only one")
    try:
        sequence_file = os.path.abspath(args[args.index("--sequence") + 1])
    except IndexError:
        error("Missing input for --sequence. Guide file file must be specified")
    if not os.path.isfile(sequence_file):
        error("The file specified for --sequence does not exist")
    return sequence_file

def _check_pam(args: List[str]) -> str:
    """Retrieves and validates the PAM file from command-line arguments.

    Ensures the --pam argument is provided and points to an existing file. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the PAM file.

    Raises:
        SystemExit: If the --pam argument is missing or the specified file does 
            not exist.
    """
    try:
        pamfile = os.path.abspath(args[args.index("--pam") + 1])
    except IndexError:
        error("Missing input for --pam. PAM file file must be specified")
        exit(1)
    if not os.path.isfile(pamfile):
        error("The file specified for --pam does not exist")
    return pamfile

def _check_be_window(args: List[str], be_window: bool) -> Tuple[int, int]:
    """Retrieves and validates the base editing window from command-line arguments.

    Ensures the --be-window argument is provided and correctly formatted, and 
    checks for required dependencies. Raises an error if the argument is missing, 
    invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        A tuple containing the start and end positions of the base editing window.

    Raises:
        SystemExit: If the --be-window argument is missing, incorrectly formatted, 
            or if --be-base is not provided when required.
    """
    if not be_window:
        return 1, 0
    if "--be-window" in args and "--be-base" not in args:
        error(
            "Missing --be-base argument. Please input the base editor to check "
            "in the specified window"
        )
    try:
        base_window = args[args.index("--be-window") + 1]
    except IndexError:
        error("Missing input for --be-window. Base editing window file must be specified")
    try:
        base_start, base_end = base_window.strip().split(",")
        base_start, base_end = int(base_start), int(base_end)
    except Exception:
        error("Invalid base editing window specified")
    if base_end < base_start:
        error("Invalid base editing window specified")
    return base_start, base_end

def _check_be_base(args: List[str], be_base: bool) -> str:
    """Retrieves and validates the base editor set from command-line arguments.

    Ensures the --be-base argument is provided and contains only valid nucleotide 
    characters, and checks for required dependencies. Raises an error if the argument 
    is missing, invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        The base editor set as a string.

    Raises:
        SystemExit: If the --be-base argument is missing, contains invalid characters, 
            or --be-window is not provided when required.
    """
    if not be_base:
        return "none"
    if "--be-base" in args and "--be-window" not in args:
        error(
            "Missing --be-window argument. Please input the base editing window "
            "to check for the specified editor"
        )
    try:
        base_set = args[args.index("--be-base") + 1]
    except IndexError:
        error("Missing input for --be-base. Base editor file must be specified")
    if any(nt not in VALID_CHARS for nt in base_set.strip().split(",")):
        error("Invalid editor specified")
    return base_set 

def _decompress_file(fname: str, outfname: str) -> str:
    """Decompresses a bgzipped file to a specified output file.

    Uses gunzip to decompress the input file and writes the result to the output 
    file. Raises an error if decompression fails.

    Args:
        fname: Path to the gzipped input file.
        outfname: Path where the decompressed file will be written.

    Returns:
        The path to the decompressed output file.

    Raises:
        SystemExit: If decompression fails.
    """
    code = subprocess.call(f"gunzip -k -c {fname} > {outfname}", shell=True)
    if code != 0:
        error("Decompressing file failed")
    assert os.path.isfile(outfname)
    return outfname

def _compress_file(fname: str) -> str:
    """Compresses a file using bgzip and returns the path to the compressed file.

    Uses bgzip to compress the specified file. Raises an error if compression fails.

    Args:
        fname: Path to the file to be compressed.

    Returns:
        The path to the compressed file with a .gz extension.

    Raises:
        SystemExit: If compression fails.
    """
    code = subprocess.call(f"bgzip -f {fname}", shell=True)
    if code != 0:
        error("Compressing and indexing file failed")
    assert os.path.isfile(f"{fname}.gz")
    return f"{fname}.gz"

def _sort_bed(fname: str, outfname: str) -> str:
    """Sorts a BED file and writes the sorted output to a new file.

    Uses sort-bed to sort the input BED file and saves the result to the specified 
    output file. Raises an error if sorting fails.

    Args:
        fname: Path to the input BED file.
        outfname: Path where the sorted BED file will be written.

    Returns:
        The path to the sorted BED file.

    Raises:
        SystemExit: If sorting fails.
    """
    code = subprocess.call(f"sort-bed {fname} > {outfname}", shell=True)
    if code != 0:
        error("Sorting BED file failed")
    assert os.path.isfile(outfname)
    return outfname

def _cat_files(fname1: str, fname2: str, outfname: str) -> str:
    """Concatenates two files and writes the result to a new file.

    Uses the cat command to combine the contents of two files into a single output 
    file. Raises an error if concatenation fails.

    Args:
        fname1: Path to the first input file.
        fname2: Path to the second input file.
        outfname: Path where the concatenated file will be written.

    Returns:
        The path to the concatenated output file.

    Raises:
        SystemExit: If concatenation fails.
    """
    code = subprocess.call(f"cat {fname1} {fname2} > {outfname}", shell=True)
    if code != 0:
        error("Concatenating files failed")
    assert os.path.isfile(outfname)
    return outfname

def _mv_file(fname: str, outfname: str) -> str:
    """Renames or moves a file to a new location.

    Uses the mv command to move or rename the specified file. Raises an error if 
    the operation fails.

    Args:
        fname: Path to the source file.
        outfname: Path to the destination file.

    Returns:
        The path to the moved or renamed file.

    Raises:
        SystemExit: If the move or rename operation fails.
    """
    code = subprocess.call(f"mv {fname} {outfname}", shell=True)
    if code != 0:
        error("Renaming file failed")
    assert os.path.isfile(outfname)
    return outfname

def _rm_files(fnames: List[str]) -> None:
    """Removes a list of files from the filesystem.

    Iterates over the provided list of file paths and deletes each file. Raises 
    an error if any file cannot be removed.

    Args:
        fnames: List of file paths to remove.

    Raises:
        SystemExit: If removing any file fails.
    """
    for fname in fnames:
        if os.path.isfile(fname):
            code = subprocess.call(f"rm {fname}", shell=True)
            if code != 0:
                error("Failed removing file")


def _writable_tmp_base() -> "Optional[str]":
    """Pick a roomy, writable base dir for annotation intermediates so a small or
    full ``/tmp`` (common in HPC/containers) does not break annotation prep (#138).

    Prefers an explicit ``$CRISPRME_TMPDIR``, then the current working directory
    (the job/output dir, which lives on the roomy data disk where the search
    already writes its results), and finally ``None`` -- letting ``tempfile`` use
    its default (``$TMPDIR`` or ``/tmp``). Never returns a read-only dir, so the
    #97 read-only-install guarantee is preserved.
    """
    for cand in (os.environ.get("CRISPRME_TMPDIR"), os.getcwd()):
        if cand and os.path.isdir(cand) and os.access(cand, os.W_OK):
            return cand
    return None


def _sort_annotation(annotationfile: str) -> str:
    """Sorts, compresses, and replaces a BED annotation file for downstream
    analysis.

    Sorts the input annotation file, compresses it with bgzip, 
    and replaces the original file. Raises an error if any step fails.

    Args:
        annotationfile: Path to the input annotation file.

    Returns:
        The path to the sorted and compressed annotation file.

    Raises:
        SystemExit: If decompression, sorting, compression, or renaming fails.
    """
    # sort-bed needs an uncompressed BED. Accept either a plain .bed or a
    # bgzipped .bed.gz (setup/complete-test download the latter). Read-only-safe:
    # decompress/sort/compress through a per-invocation temp dir and return the
    # temp .gz, so the (possibly shared or read-only) install annotation dir is
    # never written to and concurrent jobs cannot race on it (#97).
    import tempfile

    tmpdir = tempfile.mkdtemp(prefix="crisprme_annot_", dir=_writable_tmp_base())
    if annotationfile.endswith(".gz"):
        plain_path = os.path.join(tmpdir, "annotation.bed")
        _decompress_file(annotationfile, plain_path)  # gz -> temp plain
    elif os.path.isfile(annotationfile):
        plain_path = annotationfile  # plain input; sort-bed only reads it
    elif os.path.isfile(f"{annotationfile}.gz"):
        plain_path = os.path.join(tmpdir, "annotation.bed")
        _decompress_file(f"{annotationfile}.gz", plain_path)
    else:
        error(f"Annotation file not found: {annotationfile}")
    annotationfile_sorted = _sort_bed(
        plain_path, os.path.join(tmpdir, "annotation.sorted.bed")
    )
    return _compress_file(annotationfile_sorted)  # temp annotation.sorted.bed.gz


def _check_annotation(args: List[str], annotation: bool) -> str:
    """Retrieves and validates the annotation file from command-line arguments.

    Ensures the --annotation argument is provided and points to an existing file, 
    or returns a mock file if not specified. Raises an error if the argument is 
    missing or invalid.

    Args:
        args: List of command-line arguments.
        annotation: Boolean indicating if annotation is required.

    Returns:
        The absolute path to the annotation file, or a mock file if annotation is 
            not required.

    Raises:
        SystemExit: If the --annotation argument is required but missing, or the 
            specified file does not exist.
    """
    if not annotation:
        return os.path.join(script_path, "vuoto.txt")  # mock annotation file
    try:
        annotationfile = os.path.abspath(args[args.index("--annotation") + 1])
    except IndexError:
        error("Missing input for --annotation. Annotation file must be specified")
    if not os.path.isfile(annotationfile):
        error("The file specified for --annotation does not exist")
    return _sort_annotation(annotationfile)  # sort input annotation file


def _tag_personal_annotation(fname: str, outfname) -> str:
    """Tags the fourth column of a BED file as personal and writes the result to 
    a new file.

    Modifies the annotation name in the fourth column by appending '_personal', 
    replaces spaces with tabs, and updates commas in the annotation name. Raises 
    an error if tagging fails.

    Args:
        fname: Path to the input BED file.
        outfname: Path where the tagged BED file will be written.

    Returns:
        The path to the tagged BED file.

    Raises:
        SystemExit: If tagging fails.
    """
    code = subprocess.call(
        f"awk '$4 = $4\"_personal\"' {fname} | sed \"s/ /\t/g\" | sed "
        f"\"s/,/_personal,/g\" > {outfname}", 
        shell=True,
    )
    if code != 0:
        error("Tagging personal annotation file failed")
    assert os.path.isfile(outfname)
    return outfname

def _process_personal_annotation(personal_annotationfile: str, annotationfile: str) -> str:
    """Integrates personal and reference annotation files into a single sorted 
    and compressed BED file.

    Decompresses and tags the personal annotation file, merges it with the reference 
    annotation file if present, sorts the combined file, compresses it, and returns 
    the path to the final file. Raises an error if any step fails.

    Args:
        personal_annotationfile: Path to the personal annotation BED file.
        annotationfile: Path to the reference annotation BED file.

    Returns:
        The path to the sorted and compressed combined annotation file.

    Raises:
        SystemExit: If decompression, tagging, concatenation, sorting, or 
            compression fails.
    """
    # Read-only-safe: build the merged personal+reference annotation entirely in a
    # per-invocation temp dir, never next to the (possibly read-only) inputs (#97).
    import tempfile

    tmpdir = tempfile.mkdtemp(prefix="crisprme_pannot_", dir=_writable_tmp_base())
    pannotation_tag = _tag_personal_annotation(
        personal_annotationfile, os.path.join(tmpdir, "personal.tag.bed")
    )
    concat_annotationfile = os.path.join(tmpdir, "annotation+personal.bed")
    if annotationfile == os.path.join(script_path, "vuoto.txt"):
        concat_annotationfile = _mv_file(pannotation_tag, concat_annotationfile)
    else:  # concatenate personal and annotation file
        concat_annotationfile = _cat_files(
            annotationfile, pannotation_tag, concat_annotationfile
        )
    # sort the concatenated annotation, then bgzip to annotation+personal.bed.gz
    concat_sorted = _sort_bed(
        concat_annotationfile, os.path.join(tmpdir, "annotation+personal.sorted.bed")
    )
    concat_annotationfile = _mv_file(concat_sorted, concat_annotationfile)
    result = _compress_file(concat_annotationfile)  # temp annotation+personal.bed.gz
    assert os.path.isfile(result)
    return result
    

def _check_personal_annotation(args: List[str], annotationfile: str, personal_annotation: bool) -> str:
    """Retrieves and processes the personal annotation file if specified.

    Checks for the --personal_annotation argument, validates the file, and integrates 
    it with the reference annotation file. Returns the processed annotation file path, 
    or the original annotation file if no personal annotation is provided.

    Args:
        args: List of command-line arguments.
        annotationfile: Path to the reference annotation file.
        personal_annotation: Boolean indicating if personal annotation is required.

    Returns:
        The path to the processed annotation file.

    Raises:
        SystemExit: If the --personal_annotation argument is required but missing, 
            or the specified file does not exist.
    """
    if not personal_annotation:
        return annotationfile
    try:
        personal_annotationfile = os.path.abspath(args[args.index("--personal_annotation") + 1])
    except IndexError:
        error("Missing input for --personal_annotation. Annotation file must be specified")
    if not os.path.isfile(personal_annotationfile):
        error("The file specified for --personal_annotation does not exist")
    if annotationfile.endswith(".gz"):
        annotationfile = annotationfile[:-3]
    return _process_personal_annotation(personal_annotationfile, annotationfile)

def _check_samples_ids(args: List[str], variant: bool) -> str:
    """Retrieves and validates the samples ID file from command-line arguments.

    Ensures the --samplesID argument is provided and points to an existing file 
    if variant-aware search is enabled. Returns a mock file if variant is not used. 
    Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        variant: Boolean indicating if variant-aware search is enabled.

    Returns:
        The absolute path to the samples ID file, or a mock file if variant is False.

    Raises:
        SystemExit: If the --samplesID argument is missing or the specified file 
            does not exist.
    """
    if variant and "--samplesID" not in args:
        error("Missing --samplesID argument for variant-aware offtargets search")
    if not variant and "--samplesID" in args:
        error("Missing --samplesID selected, but missing --vcf argument")
    if not variant:  # use mock file for samples if variant not used
        return os.path.join(script_path, "vuoto.txt")
    try:
        samplefile = os.path.abspath(args[args.index("--samplesID") + 1])
    except IndexError:
        error("Missing input for --samplesID. Samples file must be specified")
        exit(1)
    if not os.path.isfile(samplefile):
        error("The file specified for --samplesID does not exist")
    return samplefile

def _check_gene_annotation(args: List[str], geneann: bool) -> str:
    """Retrieves and validates the gene annotation file from command-line arguments.

    Ensures the --gene_annotation argument is provided and points to an existing file, 
    or returns a mock file if not specified. Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        geneann: Boolean indicating if gene annotation is required.

    Returns:
        The absolute path to the gene annotation file, or a mock file if gene 
            annotation is not required.

    Raises:
        SystemExit: If the --gene_annotation argument is required but missing, 
            or the specified file does not exist.
    """
    if not geneann:
        return os.path.join(script_path, "vuoto.txt")
    try:
        gene_annotation = os.path.abspath(args[args.index("--gene_annotation") + 1])
    except IndexError:
        error("Missing input for --gene_annotation. Gene annotation file must be specified")
    if not os.path.isfile(gene_annotation):
        error("The file specified for --gene_annotation does not exist")
    # Sort (and bgzip/tabix) the gene annotation the same way --annotation is
    # handled. Previously this only compressed the file, so an unsorted BED
    # broke downstream processing (FDA item 4/5). _sort_annotation accepts both
    # plain .bed and bgzipped .bed.gz inputs.
    return _sort_annotation(gene_annotation)  # sort input gene annotation file

def _check_mm(args: List[str]) -> int:
    """Retrieves and validates the number of mismatches from command-line arguments.

    Ensures the --mm argument is provided and is a non-negative integer. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.

    Returns:
        The number of mismatches as an integer.

    Raises:
        SystemExit: If the --mm argument is missing or the value is negative.
    """
    try:
        mm = int(args[args.index("--mm") + 1])
    except IndexError:
        error("Missing input for --mm. Mismatches number must be specified")
    if mm < 0:
        error("Invalid number of mismatches specified")
    return mm

def _check_bdna(args: List[str], bdna: bool) -> int:
    """Retrieves and validates the number of DNA bulges from command-line arguments.

    Ensures the --bDNA argument is provided and is a non-negative integer. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        bdna: Boolean indicating if DNA bulges are required.

    Returns:
        The number of DNA bulges as an integer.

    Raises:
        SystemExit: If the --bDNA argument is missing or the value is negative.
    """
    if not bdna:
        return 0
    try:
        bDNA = int(args[args.index("--bDNA") + 1])
    except IndexError:
        error("Missing input for --bDNA. DNA bulges number must be specified")
    if bDNA < 0:
        error("Invalid number of DNA bulges specified")
    return bDNA

def _check_brna(args: List[str], brna: bool) -> int:
    """Retrieves and validates the number of RNA bulges from command-line arguments.

    Ensures the --bRNA argument is provided and is a non-negative integer. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        brna: Boolean indicating if RNA bulges are required.

    Returns:
        The number of RNA bulges as an integer.

    Raises:
        SystemExit: If the --bRNA argument is missing or the value is negative.
    """
    if not brna:
        return 0
    try:
        bRNA = int(args[args.index("--bRNA") + 1])
    except IndexError:
        error("Missing input for --bRNA. RNA bulges number must be specified")
    if bRNA < 0:
        error("Invalid number of RNA bulges specified")
    return bRNA

def _check_merge(args: List[str], merge: bool) -> int:
    """Retrieves and validates the merge threshold from command-line arguments.

    Ensures the --merge argument is provided and is a non-negative integer. 
    Returns a default value if not specified. Raises an error if the argument is 
    missing or invalid.

    Args:
        args: List of command-line arguments.
        merge: Boolean indicating if merge threshold is required.

    Returns:
        The merge threshold as an integer.

    Raises:
        SystemExit: If the --merge argument is missing or the value is negative.
    """
    if not merge:
        return 3
    try:
        merge_t = int(args[args.index("--merge") + 1])
    except IndexError:
        error("Missing input for --merge. Merge threshold must be specified")
    if merge_t < 0:
        error("Invalid merge threshold specified")
    return merge_t

def _check_sorting_criteria_scoring(args: List[str], sorting_criteria: bool) -> str:
    """Retrieves and validates the sorting criteria for scoring from command-line 
    arguments.

    Ensures the --sorting-criteria-scoring argument is provided, contains valid 
    and non-repeated criteria, and does not exceed the allowed number of criteria. 
    Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        sorting_criteria: Boolean indicating if sorting criteria for scoring is 
            required.

    Returns:
        The sorting criteria for scoring as a comma-separated string.

    Raises:
        SystemExit: If the --sorting-criteria-scoring argument is missing, contains 
            forbidden or repeated criteria, or exceeds the allowed number of criteria.
    """
    if not sorting_criteria:
        return "mm+bulges"
    try:
        sorting_criteria_scoring = args[args.index("--sorting-criteria-scoring") + 1]
    except IndexError:
        error(
            "Missing input for --sorting-criteria-scoring. Sorting criteria "
            "(scoring) must be specified"
        )
    if len(sorting_criteria_scoring.split(",")) > len(
        set(sorting_criteria_scoring.split(","))
    ):
        error("Repeated sorting criteria (scoring)\n")
    if len(sorting_criteria_scoring.split(",")) > 3:
        error("Forbidden or repeated sorting criteria (scoring)\n")
    if any(
        c not in ["mm+bulges", "mm", "bulges"]
        for c in sorting_criteria_scoring.split(",")
    ):
        error("Forbidden sorting criteria (scoring) selected\n")
    return sorting_criteria_scoring

def _check_sorting_criteria(args: List[str], sorting_criteria: bool) -> str:
    """Retrieves and validates the sorting criteria from command-line arguments.

    Ensures the --sorting-criteria argument is provided, contains valid and non-repeated 
    criteria, and does not exceed the allowed number of criteria. Raises an error if the 
    argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        sorting_criteria: Boolean indicating if sorting criteria is required.

    Returns:
        The sorting criteria as a comma-separated string.

    Raises:
        SystemExit: If the --sorting-criteria argument is missing, contains forbidden
            or repeated criteria, or exceeds the allowed number of criteria.
    """
    if not sorting_criteria:
        # Match the web default and this flag's documented default ('mm+bulges,mm'):
        # the fewest/total tie-break uses mm+bulges then mm, so CLI and web rank tied
        # targets identically. Previously this returned only 'mm+bulges', disagreeing
        # with both the web (main_page) and this command's own --help.
        return "mm+bulges,mm"
    try:
        sorting_criteria_fewest = args[args.index("--sorting-criteria") + 1]
    except IndexError:
        error(
            "Missing input for --sorting-criteria. Sorting criteria "
            "must be specified"
        )
    if len(sorting_criteria_fewest.split(",")) > len(
        set(sorting_criteria_fewest.split(","))
    ):
        error("Repeated sorting criteria\n")
    if len(sorting_criteria_fewest.split(",")) > 3:
        error("Forbidden or repeated sorting criteria\n")
    if any(
        c not in ["mm+bulges", "mm", "bulges"]
        for c in sorting_criteria_fewest.split(",")
    ):
        error("Forbidden sorting criteria selected\n")
    return sorting_criteria_fewest

def _check_output(args: List[str]) -> str:
    """Retrieves and validates the output folder from command-line arguments.

    Ensures the --output argument is provided and points to a valid directory. 
    Raises an error if the argument is missing, the folder does not exist, or is 
    not empty.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the output folder.

    Raises:
        SystemExit: If the --output argument is missing, the folder does not exist, 
            or is not empty.
    """
    try:
        outputfolder = os.path.join(
            current_working_directory, CRISPRMEDIRS[1], args[args.index("--output") + 1]
        )
    except IndexError:
        error("Missing input for --output. Output folder must be specified")
    if os.path.isdir(outputfolder):  # check whether the folder is present or not
        if is_folder_empty(outputfolder):  # if present check if not empty
            error(
                f"Output folder {outputfolder} not empty!Select another "
                "output folder for the current CRISPRme run.If the previous "
                "run using the following folder threw an error, please delete "
                f"{outputfolder} before running a new CRISPRme search."
            )
    else:  # old folder doesn't exist, create it
        os.makedirs(outputfolder)    
    if not os.path.isdir(outputfolder):
        error("The folder specified for --output does not exist")
    return outputfolder

def _check_threads(args: List[str], threads: bool) -> int:
    """Retrieves and validates the number of threads from command-line arguments.

    Ensures the --thread argument is provided and is a positive integer. Returns 
    a default value if not specified. Raises an error if the argument is missing 
    or invalid.

    Args:
        args: List of command-line arguments.
        threads: Boolean indicating if the number of threads is specified.

    Returns:
        The number of threads as an integer.

    Raises:
        SystemExit: If the --thread argument is missing or the value is less than 1.
    """
    if not threads:
        return 8  # default use 8 threads
    try:
        thread = int(args[args.index("--thread") + 1])
    except IndexError:
        error("Missing input for --thread. Number of threads must be specified")
    if thread < 1:
        error("Invalid number of threads specified")
    return thread

# App-wide bulge ceiling (usable bulges = MAX_BULGES - 1), matching pages_utils.MAX_BULGES.
# A reference index up to this depth is buildable on demand from the shipped raw genome.
MAX_BULGES = 3


def _installed_index_bulge_cap(motif: str, genome: str, is_variant: bool) -> int:
    """Bulge depth a search can REACH for this PAM motif + genome (mirrors the web's
    pages_utils.reference_bulge_capacity/index_max_bulges cap without web dependencies).

    Index folders are named ``<motif>_<N>_<genome>`` (reference) or
    ``<motif>_<N>_<genome>+<enriched>`` (variant); the usable bulge count is N-1 (the +1
    covers alignments that start with a gap). A pamless all-N index of the same length
    serves any PAM.

    The REFERENCE and VARIANT terms are ASYMMETRIC:
      * REFERENCE: the larger of the installed reference index depth and -- when the raw
        genome ``Genomes/<genome>/`` is present (has FASTA chromosome files) -- the app
        ceiling ``MAX_BULGES - 1``, because a reference index up to that depth is buildable
        on demand from the shipped raw genome (the search shell already builds it under an
        mkdir lock). This is the fix for the dict-less bug where a variant tarball ships no
        reference index, so a fresh-install variant complete-search silently derived 0
        bulges even though the raw genome IS shipped and a bmax=2 reference is buildable.
      * VARIANT: STRICTLY the installed variant index (index_max_bulges(...,dataset)); a
        variant/indel index can NOT be built dict-less (its source VCFs are not shipped),
        so this term is never relaxed to a buildable ceiling.

    Returns 0 when no index can supply bulges and no raw genome is present, so the caller
    leaves the search bulge-free rather than forcing an impossible build. Best-effort: any
    error -> 0.
    """
    try:
        lib = os.path.join(current_working_directory, "genome_library")
        if not motif or not genome:
            return 0
        pamless = "N" * len(motif)
        prefixes = (f"{motif}_", f"{pamless}_")
        suffix = f"_{genome}"

        def _scan(want_plus: bool) -> int:
            if not os.path.isdir(lib):
                return 0
            best = 0
            for d in os.listdir(lib):
                if not os.path.isdir(os.path.join(lib, d)) or d.endswith("_INDELS"):
                    continue
                base, plus, _enriched = d.partition("+")
                if want_plus != bool(plus):
                    continue
                pfx = next((p for p in prefixes if base.startswith(p)), None)
                if pfx is None or not base.endswith(suffix):
                    continue
                n_str = base[len(pfx) : -len(suffix)]
                if n_str.isdigit():
                    best = max(best, int(n_str))
            return max(0, best - 1)

        # REFERENCE term: installed reference index OR buildable-from-raw-genome ceiling.
        ref_installed = _scan(want_plus=False)
        raw = os.path.join(current_working_directory, "Genomes", genome)
        ref_buildable = 0
        try:
            if os.path.isdir(raw) and any(
                f.endswith((".fa", ".fasta")) for f in os.listdir(raw)
            ):
                ref_buildable = MAX_BULGES - 1
        except OSError:
            ref_buildable = 0
        ref_cap = max(ref_installed, ref_buildable)
        if not is_variant:
            return ref_cap
        # VARIANT search also needs the (SHIPPED, non-buildable) variant index.
        return min(ref_cap, _scan(want_plus=True))
    except Exception:
        return 0


def complete_search() -> None:
    warn_low_memory()  # non-fatal low-memory warning (Docker Desktop, etc.)
    args = input_args[2:]  # retrieve complete-search input arguments
    if "--help" in args or not args:  # print help
        print_help_complete_search()
    check_crisprme_dirtree()  # check crisprme directory tree structure
    _check_mandatory_args(args)  # check mandatory arguments
    genomedir = _check_genome(args)  # input genome folder
    variant = "--vcf" in args  # variant-aware search?
    vcfdir = _check_vcf(args, variant)  # input variants dataset
    guidefile = _check_guide(args, "--guide" in args)  # guide file
    sequence_file = _check_sequence(args, "--sequence" in args)  # sequence
    sequence_use = bool(sequence_file)  # use sequence file for guide
    assert sum([bool(guidefile), bool(sequence_file)]) == 1
    pamfile = _check_pam(args)  # pam file
    base_start, base_end = _check_be_window(args, "--be-window" in args)  # base editing window
    base_set = _check_be_base(args, "--be-base" in args)  # base editing bases
    annotationfile = _check_annotation(args, "--annotation" in args)  # annotation file
    annotationfile = _check_personal_annotation(args, annotationfile, "--personal_annotation" in args) # personal annotation file
    samplefile = _check_samples_ids(args, variant)  # samples ids file
    gene_annotation = _check_gene_annotation(args, "--gene_annotation" in args)  # gene annotation file
    mm = _check_mm(args)  # mismatches
    _bdna_given, _brna_given = "--bDNA" in args, "--bRNA" in args
    bDNA, bRNA = _check_bdna(args, _bdna_given), _check_brna(args, _brna_given)  # bulges
    bMax = max(bDNA, bRNA)  # maximum number of bulges
    merge_t = _check_merge(args, "--merge" in args)  # merge threshold
    sorting_criteria_scoring = _check_sorting_criteria_scoring(args, "--sorting-criteria-scoring" in args)  # sorting criteria score columns
    sorting_criteria = _check_sorting_criteria(args, "--sorting-criteria" in args)  # sorting criteria (mm+bulges) columns
    outputfolder = _check_output(args)  # output folder
    thread = _check_threads(args, "--thread" in args)  # number of threads
    # VCF FILTER values to treat as passing (default: "PASS,." per VCF spec §1.6.1)
    vcf_filter_pass_values = "PASS,."
    if "--vcf-filter-pass-values" in args:
        try:
            vcf_filter_pass_values = args[args.index("--vcf-filter-pass-values") + 1]
            if vcf_filter_pass_values.startswith("--"):
                raise ValueError("Please input a value for flag --vcf-filter-pass-values")
        except IndexError as e:
            raise ValueError("Missing input for --vcf-filter-pass-values") from e
    full_input_validate = "--full_input_validate" in args

    # optional prebuilt/staged reference-index library (--index-path). When
    # given, the reference index is looked up here (e.g. an index made with
    # build-index-only, or one downloaded ahead of time) rather than built under
    # the working directory; a missing index becomes a hard error downstream.
    index_path = "_"
    if "--index-path" in args:
        try:
            index_path = args[args.index("--index-path") + 1]
            if index_path.startswith("--"):
                raise ValueError("Please input a value for flag --index-path")
        except IndexError as e:
            raise ValueError("Missing input for --index-path") from e
        index_path = os.path.abspath(index_path)
        if not os.path.isdir(index_path):
            error(f"The index path {index_path} does not exist or is not a directory")

    # cap on TOTAL edits per alignment (--max-total-edits, issue #107): prevents
    # combined-edit alignments (e.g. 6mm+2+2 bulges = 10) from bloating the
    # intermediate files, scoring and post-analysis. Enforced INSIDE the TST
    # search (pruned before generation, --max-edits) with a post-search awk drop
    # as a backstop for the -r/brute-force path. Default 4 (a real off-target
    # rarely stacks many mismatches AND several bulges); -1 disables it.
    max_total_edits = 4
    if "--max-total-edits" in args:
        try:
            max_total_edits = int(args[args.index("--max-total-edits") + 1])
        except (IndexError, ValueError):
            error("Please provide a non-negative integer for --max-total-edits")
        if max_total_edits < 0:
            error("--max-total-edits must be a non-negative integer")

    # extract pam seq from file
    pam_len = 0
    total_pam_len = 0
    with open(pamfile, "r") as pam_file:
        pam_char = pam_file.readline()
        # Only the first line is used, so a multi-line PAM file would silently
        # ignore the extra motifs. Reject it with a clear error instead (FDA
        # item): a single IUPAC motif per file is supported.
        if any(line.strip() for line in pam_file):
            raise ValueError(
                "Only one PAM motif per file is supported; the PAM file "
                f"'{pamfile}' contains more than one non-empty line. Combine "
                "the motifs into a single IUPAC motif (e.g. NAG + NGG -> NRG)."
            )
        total_pam_len = len(pam_char.split(" ")[0])
        index_pam_value = pam_char.split(" ")[-1]
        if int(pam_char.split(" ")[-1]) < 0:
            end_idx = int(pam_char.split(" ")[-1]) * (-1)
            pam_char = pam_char.split(" ")[0][0:end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(" ")[-1])
            pam_char = pam_char.split(" ")[0][end_idx * (-1) :]
            pam_len = end_idx
            pam_begin = False

    # Issue #105: a partially-degenerate PAM motif (IUPAC codes W/R/Y/S/K/M/B/D/H/V)
    # combined with bulges can crash the underlying CRISPRitz search engine with a
    # heap-corruption error ("free(): invalid pointer"). The observed trigger is an
    # odd-length degenerate motif (e.g. WTN); even-length ones such as TTTV (Cas12a)
    # are stable, so this is a scoped, NON-FATAL warning (never blocks a valid run).
    # The root cause is fixed in CRISPRitz 2.7.1; CRISPRme will pin to it in a later
    # release. Until then, surface a clear message instead of a cryptic C++ abort.
    if (
        bMax > 0
        and pam_len % 2 == 1
        and any(c in "WRYSKMBDHV" for c in pam_char.upper())
    ):
        sys.stderr.write(
            f"WARNING: PAM motif '{pam_char}' is a partially-degenerate, odd-length "
            "motif used with bulges. This combination can crash the underlying "
            "search engine (issue #105, fixed in CRISPRitz 2.7.1). If the search "
            "aborts with a memory error, retry without bulges or use a "
            "fully-specified PAM.\n"
        )

    genome_ref = os.path.basename(genomedir)
    annotation_name = os.path.basename(annotationfile)
    # The nuclease name is parsed from the PAM filename, which must follow the
    # <length>-<motif>-<CasName>.txt convention (e.g. 20bp-NGG-SpCas9.txt).
    # Validate before indexing to avoid a cryptic IndexError (FDA item).
    pam_name_fields = os.path.basename(pamfile).split(".")[0].split("-")
    if len(pam_name_fields) < 3:
        raise ValueError(
            f"Invalid PAM filename '{os.path.basename(pamfile)}': expected the "
            "'<length>-<motif>-<CasName>.txt' convention "
            "(e.g. 20bp-NGG-SpCas9.txt)."
        )
    nuclease = pam_name_fields[2]
    # Treat --max-total-edits as the single "max edits" knob (mirroring the web slider):
    # when the user did NOT specify per-type bulges, derive bDNA=bRNA from the max-edits
    # budget, BOUNDED by the bulge depth a search can REACH here -- so e.g.
    # `complete-search --vcf ... --max-total-edits 4` searches up to 2 bulges of each type
    # (finding 2mm+1bulge / 2mm+2bulge patterns) with no bulge flags. The reference term of
    # the cap is buildable-aware (the reference index is built on demand from the shipped
    # raw genome, as the search shell does), so a fresh/dict-less install no longer silently
    # derives 0 bulges; the variant term stays strictly installed-index-based (a variant
    # index can't be built dict-less). Explicit --bDNA/--bRNA always win; if no index can
    # supply bulges and no raw genome exists the search stays bulge-free (fast, safe).
    if not _bdna_given and not _brna_given and max_total_edits > 0:
        _idx_cap = _installed_index_bulge_cap(pam_char, genome_ref, variant)
        _derived = min(max_total_edits, _idx_cap)
        if _derived > 0:
            bDNA = bRNA = _derived
            bMax = max(bDNA, bRNA)
            print(
                f"[complete-search] no --bDNA/--bRNA given: deriving up to {_derived} "
                f"bulge(s) of each type from --max-total-edits {max_total_edits} "
                f"(reachable index bulge depth {_idx_cap})."
            )
    # [max-total-edits] Surface the silent combined-edit prune (issue #107): when the
    # requested mm + bulges exceed the cap, any alignment stacking more than
    # max_total_edits edits is PRUNED inside the TST search -- the exact reason a deep
    # off-target (e.g. chr14:63727708 = 5mm + 2 bulges = 7 edits) is missed on a default
    # run. Keep the conservative default (advanced users raise it), but never silent.
    if 0 <= max_total_edits < mm + bDNA + bRNA:
        print(
            f"WARNING [complete-search]: --max-total-edits {max_total_edits} is below the "
            f"requested {mm}mm + {bDNA} DNA + {bRNA} RNA bulges = {mm + bDNA + bRNA} total "
            f"edits. Alignments needing more than {max_total_edits} COMBINED edits are "
            f"PRUNED (e.g. a 5mm+2-bulge = 7-edit off-target is dropped). Raise "
            f"--max-total-edits to {mm + bDNA + bRNA} to keep such deep off-targets.",
            flush=True,
        )
    if bMax != 0:
        search_index = True
    else:
        search_index = False
    if variant:
        genome_idx_list = []
        with open(vcfdir, "r") as vcfs:
            for line in vcfs:
                if line.strip():
                    if line[-2] == "/":
                        line = line[:-2]
                    base_vcf = os.path.basename(line)
                    genome_idx_list.append(
                        pam_char
                        + "_"
                        + str(bMax + 1)  # for alignments starting with gaps
                        + "_"
                        + genome_ref
                        + "+"
                        + base_vcf.strip()
                    )
        genome_idx = ",".join(genome_idx_list)
        ref_comparison = True
    else:
        genome_idx = pam_char + "_" + str(bMax + 1) + "_" + genome_ref
        ref_comparison = False
    # os.chdir(script_path)
    # write crisprme version to file
    with open(outputfolder + "/.command_line.txt", "w") as p:
        p.write("input_command\t" + " ".join(sys.argv[:]))
        p.write("\n")
        p.close()
    with open(outputfolder + "/.version.txt", "w") as p:
        p.write("crisprme_version\t" + __version__)
        p.write("\n")
        p.close()
    # write parameters to file
    with open(outputfolder + "/Params.txt", "w") as p:
        p.write("Genome_selected\t" + genome_ref.replace(" ", "_") + "\n")
        p.write("Genome_ref\t" + genome_ref + "\n")
        if search_index:
            p.write("Genome_idx\t" + genome_idx + "\n")
        else:
            p.write("Genome_idx\t" + "None\n")
        p.write("Pam\t" + pam_char + "\n")
        p.write("Max_bulges\t" + str(bMax) + "\n")
        p.write("Mismatches\t" + str(mm) + "\n")
        p.write("DNA\t" + str(bDNA) + "\n")
        p.write("RNA\t" + str(bRNA) + "\n")
        # the binding total-edits cap the search actually used (the "Max edits" the
        # user set): the report + web read this; without it they showed "n/a".
        p.write("Max_total_edits\t" + str(max_total_edits) + "\n")
        # Persist the silent-prune WARN into the run sidecar so it reaches
        # .Params.txt -> report.zip. The parent-stdout WARN (printed earlier in
        # complete_search) is NOT captured: log_verbose.txt only redirects the child
        # job. generate_report reads .Params.txt as a kv sidecar AND bundles it
        # verbatim. Single physical line (no embedded newline) so the reader parses
        # it cleanly. Guarded by the SAME condition as the stdout WARN.
        if 0 <= max_total_edits < mm + bDNA + bRNA:
            p.write(
                "Pruning_note\t"
                + (
                    f"--max-total-edits {max_total_edits} is below the requested "
                    f"{mm}mm + {bDNA} DNA + {bRNA} RNA bulges = {mm + bDNA + bRNA} total "
                    f"edits; alignments needing more than {max_total_edits} COMBINED edits "
                    f"were PRUNED (e.g. a 5mm+2-bulge = 7-edit off-target is dropped). "
                    f"Raise --max-total-edits to {mm + bDNA + bRNA} to keep such deep "
                    f"off-targets."
                )
                + "\n"
            )
        p.write("Annotation\t" + str(annotation_name) + "\n")
        p.write("Nuclease\t" + str(nuclease) + "\n")
        # p.write('Gecko\t' + str(gecko_comp) + '\n')
        p.write("Ref_comp\t" + str(ref_comparison) + "\n")
        p.close()
    len_guide_sequence = total_pam_len - pam_len
    if sequence_use:
        guides = list()
        text_sequence = str()
        for line in open(sequence_file, "r"):
            text_sequence += line
        for name_and_seq in text_sequence.split(">"):
            if "" == name_and_seq:
                continue
            name = name_and_seq[: name_and_seq.find("\n")]
            seq = name_and_seq[name_and_seq.find("\n") :]
            # seq = seq.strip().split()
            # seq = ''.join(seq)
            seq = seq.strip()
            # name, seq = name_and_seq.strip().split('\n')
            if "chr" in seq:
                # extracted_seq = extract_seq.extractSequence(
                #         name, seq, genome_ref.replace(' ', '_'))
                for single_row in seq.split("\n"):
                    if "" == single_row:
                        continue
                    pieces_of_row = single_row.strip().split()
                    seq_to_extract = (
                        pieces_of_row[0]
                        + ":"
                        + pieces_of_row[1]
                        + "-"
                        + pieces_of_row[2]
                    )
                    extracted_seq = extractSequence(
                        name, seq_to_extract, genome_ref.replace(" ", "_")
                    )
                    guides.extend(
                        getGuides(
                            extracted_seq, pam_char, len_guide_sequence, pam_begin
                        )
                    )
            else:
                seq = seq.split()
                seq = "".join(seq)
                extracted_seq = seq.strip()
                guides.extend(
                    getGuides(extracted_seq, pam_char, len_guide_sequence, pam_begin)
                )
        temp_guides = list()
        for guide in guides:
            addN = "N" * pam_len
            if pam_begin:
                temp_guides.append(addN + guide)
            else:
                temp_guides.append(guide + addN)
        if len(temp_guides) > 1000000000:
            temp_guides = temp_guides[:1000000000]
        guides = temp_guides
        extracted_guides_file = open(outputfolder + "/guides.txt", "w")
        for guide in guides:
            extracted_guides_file.write(guide + "\n")
        extracted_guides_file.close()
    # print(guides)
    # exit(0)
    void_mail = "_"
    if sequence_use == False:
        # shutil (no shell) so a guide/output path containing spaces does not break
        shutil.copyfile(guidefile, os.path.join(outputfolder, "guides.txt"))

    # pre-flight input validation (lightweight tier, always on): catches
    # misconfigurations that would otherwise only surface deep into the run
    vcf_dataset_dirs = (
        resolve_vcf_dataset_dirs(vcfdir, current_working_directory) if variant else []
    )
    validation_report = run_lightweight(
        genomedir,
        vcf_dataset_dirs,
        os.path.join(outputfolder, "guides.txt"),
        pamfile,
        gene_annotation,
        samplefile,
        current_working_directory,
    )
    validation_report.write()

    # opt-in full-file scan (--full_input_validate): slower, so only runs on
    # request; still checked before the pipeline subprocess launches
    full_validation_report = None
    if full_input_validate and variant:
        full_validation_report = run_full(genomedir, vcf_dataset_dirs)
        full_validation_report.write()

    if validation_report.has_errors() or (
        full_validation_report is not None and full_validation_report.has_errors()
    ):
        sys.exit(1)

    print(
        f"Launching job {outputfolder}. The stdout is redirected in log_verbose.txt and stderr is redirected in log_error.txt"
    )
    # start search with set parameters
    with open(f"{outputfolder}/log_verbose.txt", "w") as log_verbose:
        with open(f"{outputfolder}/log_error.txt", "w") as log_error:
            crisprme_run = (
                f"{os.path.join(script_path, 'submit_job_automated_new_multiple_vcfs.sh')} "
                f"{genomedir} {vcfdir} {os.path.join(outputfolder, 'guides.txt')} "
                f"{pamfile} {annotationfile} {samplefile} {bMax + 1} {mm} {bDNA} {bRNA} "
                f"{merge_t} {outputfolder} {script_path} {thread} {current_working_directory} "
                f"{gene_annotation} {void_mail} {base_start} {base_end} {base_set} "
                f"{sorting_criteria_scoring} {sorting_criteria} {cicd_test} "
                f"{vcf_filter_pass_values} {index_path} {max_total_edits}"
            )
            code = subprocess.call(
                crisprme_run, shell=True, stderr=log_error, stdout=log_verbose
            )
            if code != 0:
                # surface WHERE it failed + the actual error, not just "failed"
                sys.stderr.write(
                    "\nCRISPRme run failed.\n"
                    + summarize_pipeline_failure(outputfolder)
                    + "\n"
                )
                raise OSError(
                    f"\nCRISPRme run failed! See {os.path.join(outputfolder, 'log_error.txt')} for details\n"
                )
            # subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), str(vcfdir), str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(
            #     samplefile), str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
    # else:
    #     with open(f"{outputfolder}/log_verbose.txt", 'w') as log_verbose:
    #         with open(f"{outputfolder}/log_error.txt", 'w') as log_error:
    #             subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), '_', str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(script_path+'vuoto.txt'),
    #                             str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
    # change name of guide and param files to hidden
    os.system(f"mv {outputfolder}/guides.txt {outputfolder}/.guides.txt")
    os.system(f"mv {outputfolder}/Params.txt {outputfolder}/.Params.txt")


def print_help_build_index() -> None:
    """Prints detailed help information for the build-index-only functionality.

    Outputs a description of the index-build step and lists all available
    command-line options to stderr, then exits the program.
    """
    sys.stderr.write(
        "The build-index-only functionality pre-builds the CRISPRitz reference "
        "genome index that bulge-enabled searches rely on, WITHOUT running a "
        "search. The index is written under <path>/genome_library/ and is then "
        "reused automatically by any later 'complete-search' run (launched from "
        "the same working directory, or pointed at it with --index-path) that "
        "uses the same genome, PAM and bulge settings. Building the index once "
        "and reusing it avoids repeating the single most expensive step across "
        "many searches, and lets an index be staged (or downloaded) ahead of "
        "time on a large machine.\n"
    )
    sys.stderr.write(
        "Options:\n"
        "\t--genome, specify the reference genome folder [REQUIRED]\n"
        "\t--pam, specify a file containing the PAM sequence [REQUIRED]\n"
        "\t--bDNA, number of DNA bulges the index must support [OPTIONAL, "
        "default 0]\n"
        "\t--bRNA, number of RNA bulges the index must support [OPTIONAL, "
        "default 0]\n"
        "\t--thread, set number of threads to use [default: 8]\n"
        "\t--vcf, a VCF dataset directory (e.g. VCFs/1000G). When given, also "
        "pre-builds the variant-aware index: enriches the genome with the VCF "
        "and indexes the enriched (SNP) and indels genomes, so the first "
        "variant-aware search does not pay the enrichment/indexing cost "
        "[OPTIONAL]\n"
        "\t--samplesID, a listing file (one samplesID filename per line, under "
        "samplesIDs/; a combined panel lists both 1000G and HGDP). When given "
        "with --vcf, the build ALSO emits the additive dictless tiers (Tier-0 "
        "registry + Tier-1 genotype store) per chromosome alongside the dicts, so "
        "the built index carries the fast-post-analysis tiers; omit to build dicts "
        "only [OPTIONAL]\n"
        "\t--path, working directory under which genome_library/ is created "
        "[OPTIONAL, default: current directory]\n"
        "\t--name, human-friendly label for the finished index, written to a "
        "'.display_label' sidecar so the web index list / search form show it "
        "instead of the auto <motif>_<N>_<genome> name (falls back to the "
        "convention when omitted) [OPTIONAL]\n"
    )
    sys.exit(1)


def _write_display_label(index_dir: str, label: "str | None") -> None:
    """Persists an optional human-friendly label for a built index.

    Writes ``<index_dir>/.display_label`` when a non-empty ``label`` is given so
    the web index list / search form (pages_utils.get_available_indexes) shows it
    instead of the auto <motif>_<N>_<genome> convention. A no-op for a blank/None
    label, so the convention-based fallback stays the default. Publishing bundles
    the sidecar into the HF tarball (crisprme_hf.publish_index).
    """
    if not label or not label.strip():
        return
    if not os.path.isdir(index_dir):
        return
    try:
        with open(os.path.join(index_dir, ".display_label"), "w") as fd:
            fd.write(label.strip() + "\n")
    except OSError as exc:  # non-fatal: the index is built, only the label failed
        sys.stderr.write(f"Warning: could not write index display label: {exc}\n")


def _db_name_from_samplesid(samples_filename: str) -> str:
    """Derive a clean per-database label from a samplesID filename.

    The label is used ONLY as the dataset-provenance name carried in the Tier-0 /
    Tier-1 per-sample meta (``read_samplesid(path, database)``), so it just needs to
    be stable and human-readable. samplesID files ship under a few conventions
    (``hg38_1000G.samplesID.txt``, ``samplesIDs.HGDP.txt``, ``HGDP.samplesID.txt``),
    so strip the common ``samplesID(s)`` token and the reference-genome prefix and
    keep the meaningful middle (e.g. "1000G", "HGDP", "gnomad.v41"). Never fails:
    falls back to the bare basename.
    """
    base = os.path.basename(samples_filename).strip()
    # drop a trailing .txt / .txt.gz
    for suf in (".txt.gz", ".txt", ".gz"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    # drop the samplesID token wherever it sits (prefix "samplesIDs.", infix
    # ".samplesID", suffix ".samplesID")
    for tok in ("samplesIDs.", "samplesIDs_", "samplesID."):
        if base.startswith(tok):
            base = base[len(tok):]
    for tok in (".samplesID", "_samplesID", ".samplesIDs", "_samplesIDs"):
        base = base.replace(tok, "")
    # drop a leading reference-genome token like "hg38_"
    for ref_tok in ("hg38_", "hg19_", "GRCh38_", "GRCh37_", "T2T_", "mm10_", "mm39_"):
        if base.startswith(ref_tok):
            base = base[len(ref_tok):]
            break
    base = base.strip("._-")
    return base or os.path.basename(samples_filename)


def _build_db_to_samplesid(samples_listing: str, workdir: str):
    """Build the ORDERED {db_name: samplesID_path} map for tier emission.

    ``samples_listing`` is the ``--samplesID`` LISTING file (the SAME format
    complete-search uses: one samplesID FILENAME per line, each resolved under
    ``<workdir>/samplesIDs/``). For the combined 1000G+HGDP panel this lists BOTH
    files, so the returned map has BOTH entries (order preserved). Blank lines and
    ``#`` comment lines are skipped. Returns {} if nothing usable is found (the
    caller then skips tier emission -- additive, never fatal). Missing listed files
    are noted to STDOUT and skipped (a partial map still emits usable tiers).
    """
    import collections

    db_map = collections.OrderedDict()
    if not samples_listing or not os.path.isfile(samples_listing):
        return db_map
    samples_dir = os.path.join(workdir, "samplesIDs")
    with open(samples_listing) as fh:
        for line in fh:
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            # accept either a bare filename (resolved under samplesIDs/) or an
            # already-absolute/relative path that exists as given.
            candidate = name
            if not os.path.isabs(candidate) or not os.path.isfile(candidate):
                under_dir = os.path.join(samples_dir, name)
                if os.path.isfile(under_dir):
                    candidate = under_dir
            if not os.path.isfile(candidate):
                print(
                    "build-index-only: NOTE samplesID entry %r not found under "
                    "%s (skipping it for tier emission)" % (name, samples_dir),
                    flush=True,
                )
                continue
            db_name = _db_name_from_samplesid(name)
            if db_name in db_map:
                # Two samplesID files reduced to the SAME dataset label. Overwriting
                # silently would drop the earlier dataset's per-sample provenance
                # (violating "never conflate datasets") -- make it VISIBLE on stdout.
                print(
                    "build-index-only: NOTE two samplesID files map to the same "
                    "dataset label %r (%r overwrites the earlier entry); per-sample "
                    "provenance for the earlier dataset would be lost -- give the "
                    "files distinct names to keep datasets separate." % (db_name, name),
                    flush=True,
                )
            db_map[db_name] = os.path.abspath(candidate)
    return db_map


def _emit_combined_samplesid(workdir, vcf_name, genome_ref, db_to_samplesid):
    """Emit the combined ``samplesIDs/<vcf_name>.samplesID.txt`` for a MERGED
    variant index into the install, reusing the ONE union implementation
    (:func:`crisprme_hf.synthesize_combined_samplesid`).

    This is the build-side complement of the download-side synthesizer: a user
    who BUILDS + PUBLISHES their own merged index must ship the combined
    samplesID so ``download --what index`` + search is self-complete (the
    search's ``--samplesID`` listing expects that single combined file, and the
    download-side generation only unions per-db files that are already on HF).

    ``db_to_samplesid`` is the ORDERED ``{db_label: samplesID_path}`` map already
    computed by :func:`_build_db_to_samplesid` (reused, never recomputed). The map
    KEYS are clean db labels (``1000G``, ``HGDP``) but the map VALUES may point at
    arbitrarily-named source files (e.g. ``samplesIDs.HGDP.txt``); the shared
    synthesizer keys off the STRICT ``<ref>_<db>.samplesID.txt`` convention, so we
    first copy each per-db source into ``samplesIDs/`` under that CANONICAL name
    (only when absent — never clobber a shipped file), then call the synthesizer.

    NO-OP (returns ``None``) when the map has < 2 usable datasets (single dataset
    => no combined file is needed; the per-db list IS the search's samplesID) — the
    exact guard the synthesizer no-ops on. Fully guarded: any error is reported to
    STDOUT (stderr is fatal in the post-analysis context) and swallowed so a build
    never aborts over this additive convenience file. Returns the written combined
    path, or ``None`` (no-op / already present / a component missing).
    """
    try:
        # Gate: only a MERGED panel (>= 2 datasets) needs a combined file. A single
        # entry (or two files collapsing to one label) => single dataset => NO-OP.
        if not db_to_samplesid or len(db_to_samplesid) < 2:
            return None
        import shutil as _shutil

        sdir = os.path.join(workdir, "samplesIDs")
        os.makedirs(sdir, exist_ok=True)
        # Normalize the per-db source files to the strict canonical names the
        # shared synthesizer resolves by (<ref>_<db>.samplesID.txt). Copy in only
        # when the canonical target is absent, so a shipped file is never touched.
        for db_label, src in db_to_samplesid.items():
            canonical = os.path.join(sdir, f"{genome_ref}_{db_label}.samplesID.txt")
            if not os.path.isfile(canonical) and os.path.isfile(src):
                _shutil.copy2(src, canonical)
        # Single implementation of the union (header-less, dedup by SAMPLE_ID,
        # stable dataset order) — mirrors pages/main_page._ensure_samplesid.
        written = synthesize_combined_samplesid(workdir, vcf_name, ref=genome_ref)
        if written:
            print(
                f"build-index-only: wrote combined samplesID "
                f"{os.path.basename(written)} into {sdir} (union of "
                f"{', '.join(db_to_samplesid)}); a published index is now "
                f"self-complete for search.",
                flush=True,
            )
        else:
            print(
                f"build-index-only: combined samplesID for {vcf_name} not "
                f"synthesized (already present, single dataset, or a per-db "
                f"component missing) -> no-op.",
                flush=True,
            )
        return written
    except Exception as _sid_err:  # additive convenience -> never abort a build
        print(
            f"WARNING [build-index-only]: could not emit combined samplesID for "
            f"{vcf_name} ({_sid_err}); the index/dicts are still built. Provide "
            f"samplesIDs/{vcf_name}.samplesID.txt manually if a merged search "
            f"needs it.",
            flush=True,
        )
        return None


def _write_pam_build(index_dir: str, pamfile: str) -> None:
    """Records the PAM an index was built with, into ``<index_dir>/.pam_build``.

    A CRISPRitz index is built over k-mers of an exact length + PAM orientation, so
    it can only search PAMs of matching geometry. Persisting the source PAM's name
    and geometry (``<pam_name> <seq_len> <pam_pos>``) lets the web PAM selector
    (pages_utils.get_pam_options / index_build_pam) offer only compatible PAMs —
    critical for pamless (NNN) indexes, which otherwise look like they serve every
    PAM. The sidecar lives inside the index dir, so it travels in the publish
    tarball and download automatically. Non-fatal on any error.
    """
    if not os.path.isdir(index_dir):
        return
    try:
        with open(pamfile) as fh:
            parts = fh.readline().split()
        seq, pos = parts[0], int(parts[1])
        name = os.path.splitext(os.path.basename(pamfile))[0]
        with open(os.path.join(index_dir, ".pam_build"), "w") as fd:
            fd.write(f"{name} {len(seq)} {pos}\n")
    except (OSError, ValueError, IndexError) as exc:
        sys.stderr.write(f"Warning: could not write index PAM provenance: {exc}\n")


def build_index_only() -> None:
    """Pre-builds the CRISPRitz reference index for a genome+PAM+bulge combo.

    Exposes the index-build step that ``complete-search`` otherwise performs
    implicitly, so an index can be constructed once (e.g. on a large machine)
    and reused across many searches. The produced directory
    (``genome_library/<PAM>_<bMax+1>_<genome>``) is byte-for-byte the same one
    ``complete-search`` looks for, so no extra bookkeeping is required: a later
    search with matching --genome/--pam/--bDNA/--bRNA finds and reuses it.
    """
    args = input_args[2:]  # retrieve build-index-only input arguments
    if "--help" in args:
        print_help_build_index()
    genomedir = _check_genome(args)  # reference genome folder
    pamfile = _check_pam(args)  # PAM file
    thread = _check_threads(args, "--thread" in args)  # number of threads
    bDNA = _check_bdna(args, "--bDNA" in args)  # DNA bulges
    bRNA = _check_brna(args, "--bRNA" in args)  # RNA bulges
    bMax = max(bDNA, bRNA)  # maximum number of bulges
    # optional human-friendly label for the finished index; written to a
    # .display_label sidecar so the web index list / search form shows it
    # instead of the auto <motif>_<N>_<genome> convention (falls back to the
    # convention when absent — see pages_utils._friendly_index_label).
    display_name = None
    if "--name" in args:
        display_name = args[args.index("--name") + 1].strip()
    if bMax == 0:
        sys.stderr.write(
            "Nothing to do: a genome index is only required for bulge-enabled "
            "searches (--bDNA/--bRNA > 0). Zero-bulge searches run directly on "
            "the FASTA with no index.\n"
        )
        sys.exit(0)
    # working directory under which genome_library/ is created (mirrors the
    # current_working_directory complete-search uses to locate genome_library/)
    workdir = os.getcwd()
    if "--path" in args:
        workdir = os.path.abspath(args[args.index("--path") + 1])
        if not os.path.isdir(workdir):
            error(f"The working directory {workdir} does not exist")
    # derive the true PAM string exactly as complete-search does (see the PAM
    # parsing block in complete_search), so the index folder name matches the
    # one a later search will look for
    with open(pamfile, "r") as pf:
        pam_char = pf.readline()
        idx_val = int(pam_char.split(" ")[-1])
        end_idx = abs(idx_val)
        if idx_val < 0:  # 5' PAM (e.g. Cas12a)
            pam_char = pam_char.split(" ")[0][0:end_idx]
        else:  # 3' PAM (e.g. SpCas9)
            pam_char = pam_char.split(" ")[0][end_idx * (-1):]
    genome_ref = os.path.basename(genomedir)
    # complete-search indexes with bMax+1 ("for alignments starting with gaps")
    index_name = f"{pam_char}_{bMax + 1}_{genome_ref}"
    idx_folder = os.path.join(workdir, "genome_library", index_name)
    vcf_given = "--vcf" in args
    # ---- reference index ------------------------------------------------------
    if os.path.isdir(idx_folder):
        print(f"Reference index already present: {idx_folder}", flush=True)
    else:
        os.makedirs(os.path.join(workdir, "genome_library"), exist_ok=True)
        print(
            f"Building reference index {index_name} (bMax {bMax + 1}, "
            f"{thread} thread(s))...",
            flush=True,
        )
        # exactly the invocation complete-search uses, run from workdir so the
        # index lands in <workdir>/genome_library/ under the expected name
        index_cmd = (
            f"crispritz.py index-genome {genome_ref} {genomedir}/ {pamfile} "
            f"-bMax {bMax + 1} -th {thread}"
        )
        code = subprocess.call(index_cmd, shell=True, cwd=workdir)
        if code != 0 or not os.path.isdir(idx_folder):
            error(f"Reference genome indexing failed (expected {idx_folder})")
        print(f"Index built: {idx_folder}", flush=True)
    _write_pam_build(idx_folder, pamfile)
    if not vcf_given:
        _write_display_label(idx_folder, display_name)
        print(
            "It will be reused automatically by complete-search runs launched from "
            "this working directory (or pointed here with --index-path) that use the "
            "same genome, PAM and bulge settings.",
            flush=True,
        )
        return
    # ---- variant (enriched SNP + indels) index --------------------------------
    # Pre-build the variant-aware index so the first variant search does not pay
    # the (slow) enrichment + indexing cost. This mirrors STEP 1-2 of the search
    # pipeline (submit_job_automated_new_multiple_vcfs.sh) exactly, producing the
    # same folder names a later search looks for.
    import shutil
    import tempfile
    import gzip
    from glob import glob as _glob

    vcfdir = os.path.abspath(args[args.index("--vcf") + 1])
    if not os.path.isdir(vcfdir):
        error(f"The VCF dataset directory {vcfdir} does not exist")
    vcf_name = os.path.basename(vcfdir.rstrip("/"))
    # OPTIONAL --samplesID: a listing file (one samplesID filename per line under
    # samplesIDs/, the same format complete-search uses; combined panels list both
    # 1000G and HGDP). When given, the build ALSO emits the ADDITIVE dictless tiers
    # (Tier-0 registry + Tier-1 genotype store) per chromosome, alongside the dicts,
    # so a freshly-built variant index carries the fast-post-analysis tiers. Absent
    # => dicts only (legacy behavior, unchanged).
    samples_listing = None
    if "--samplesID" in args:
        try:
            samples_listing = os.path.abspath(args[args.index("--samplesID") + 1])
        except IndexError:
            error("Missing input for --samplesID. A samples listing file must be "
                  "specified")
        if not os.path.isfile(samples_listing):
            error(f"The file specified for --samplesID does not exist: "
                  f"{samples_listing}")
    enriched = os.path.join(workdir, "Genomes", f"{genome_ref}+{vcf_name}")
    indels_out = os.path.join(workdir, "Genomes", f"{genome_ref}+{vcf_name}_INDELS")
    dict_folder = os.path.join(workdir, "Dictionaries", f"dictionaries_{vcf_name}")
    indel_dict = os.path.join(workdir, "Dictionaries", f"log_indels_{vcf_name}")
    snp_idx = os.path.join(
        workdir, "genome_library", f"{pam_char}_{bMax + 1}_{genome_ref}+{vcf_name}"
    )
    indels_idx = snp_idx + "_INDELS"
    # STEP 1: enrich the genome with the VCF (crispritz add-variants). add-variants
    # writes a fixed-name variants_genome/ in its cwd, so run it in a throwaway
    # temp dir to avoid collisions, then move the outputs into place.
    if not os.path.isdir(enriched):
        print(f"Enriching {genome_ref} with {vcf_name} (add-variants)...", flush=True)
        tmp = tempfile.mkdtemp(prefix="run_", dir=os.path.join(workdir, "Genomes"))
        variants_tmp = os.path.join(tmp, "variants_genome")
        code = subprocess.call(
            f"crispritz.py add-variants {vcfdir}/ {genomedir}/ true",
            shell=True,
            cwd=tmp,
        )
        if code != 0:
            shutil.rmtree(tmp, ignore_errors=True)
            error(f"Genome enrichment (add-variants) failed on {vcf_name}")
        # This block only runs on a FRESH enrichment (enriched was absent), so any
        # pre-existing dataset-specific output dirs are stale/partial from a prior
        # aborted build. Clear them first, else shutil.move onto an existing
        # fake_chrN / dict raises "Destination path already exists" (issue #54).
        for _d in (indels_out, dict_folder, indel_dict):
            if os.path.isdir(_d):
                shutil.rmtree(_d)
        os.makedirs(indels_out, exist_ok=True)
        shutil.move(
            os.path.join(variants_tmp, "SNPs_genome", f"{genome_ref}_enriched"),
            enriched,
        )
        for f in _glob(os.path.join(variants_tmp, "fake*")):
            shutil.move(f, indels_out)
        os.makedirs(dict_folder, exist_ok=True)
        os.makedirs(indel_dict, exist_ok=True)
        for f in _glob(os.path.join(variants_tmp, "SNPs_genome", "*.json")):
            shutil.move(f, dict_folder)
        for f in _glob(os.path.join(variants_tmp, "SNPs_genome", "log*.txt")):
            shutil.move(f, indel_dict)
        # Compress the per-chromosome dicts in place (my_dict_*.json -> .json.gz,
        # log*.txt -> .txt.gz) so the published index ships them ~3.5x smaller
        # (~40-50GB not ~152GB for 1000G+HGDP) and the variant post-analysis reads
        # them gzipped on the fly (no 150GB decompress). pigz (parallel) when
        # available, else Python gzip. Mirrors _make_index_tarball's pigz preference.
        _pigz = shutil.which("pigz")
        for _d, _pat in ((dict_folder, "*.json"), (indel_dict, "log*.txt")):
            _files = _glob(os.path.join(_d, _pat))
            if not _files:
                continue
            if _pigz:
                subprocess.call([_pigz, "-f", "-p", str(thread), *_files])
            else:
                for _f in _files:
                    with open(_f, "rb") as _src, gzip.open(_f + ".gz", "wb") as _dst:
                        shutil.copyfileobj(_src, _dst)
                    os.remove(_f)
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print(f"Enriched genome already present: {enriched}", flush=True)
    # STEP 1b (ADDITIVE): emit the dictless tiers (Tier-0 registry + Tier-1 genotype
    # store) per chromosome, FROM the per-sample dicts just produced -- WITHOUT
    # touching the dicts. The search side (new_simple_analysis.py) already CONSUMES
    # these as sibling registry_<vcf>/ + genotypes_<vcf>/ dirs; this is the build
    # side that PRODUCES them. GATED on a usable --samplesID map; GUARDED so any
    # failure logs a WARNING to STDOUT (stderr is fatal here) and never aborts the
    # dict/index build. Path derivation MIRRORS the search resolvers exactly (the
    # helper swaps the dictionaries_ prefix), so a built install actually uses them.
    # NOTE (AF denominator, #117/#121): the panel denominator (AN) is built by
    # tier0_compile.build_sample_meta over EXACTLY the per-db samplesID files listed
    # here, so those files MUST be VCF-FILTERED to the genotyped samples (the merge
    # script build_combined_panel.sh writes filtered per-db lists + a listing) or the
    # panel would over-count phantom hom-ref individuals and inflate AN.
    db_to_samplesid = _build_db_to_samplesid(samples_listing, workdir)
    if not db_to_samplesid:
        if samples_listing:
            print(
                "build-index-only: no usable samplesID entries -> skipping dictless "
                "tier emission (dicts unaffected; search falls back to the dict "
                "path).",
                flush=True,
            )
        else:
            print(
                "build-index-only: no --samplesID given -> skipping dictless tier "
                "emission (dicts still built; pass --samplesID to also emit the "
                "fast-post-analysis tiers).",
                flush=True,
            )
    else:
        try:
            import build_dictless_tiers as _bdt
        except Exception as _bdt_err:  # tier modules absent in this deploy
            _bdt = None
            print(
                "WARNING [build-index-only]: dictless tier modules unavailable "
                f"({_bdt_err}); building dicts only. Search falls back to the dict "
                "path.",
                flush=True,
            )
        if _bdt is not None:
            # Enumerate the chromosomes from the dict files just written. A chromosome
            # normally has exactly ONE of my_dict_<chrom>.json / .json.gz (the gzip
            # step removes the plain form); if a partial/interrupted gzip left BOTH,
            # dedupe by <chrom> (preferring the plain form) so its tiers emit once.
            def _chrom_of(_p):
                _s = os.path.basename(_p)[len("my_dict_"):]
                for _suf in (".json.gz", ".json"):
                    if _s.endswith(_suf):
                        return _s[: -len(_suf)]
                return _s
            _by_chrom = {}
            for _p in _glob(os.path.join(dict_folder, "my_dict_*.json.gz")):
                _by_chrom[_chrom_of(_p)] = _p
            for _p in _glob(os.path.join(dict_folder, "my_dict_*.json")):
                _by_chrom[_chrom_of(_p)] = _p  # plain overwrites gz -> plain preferred
            _dict_files = [_by_chrom[_c] for _c in sorted(_by_chrom)]
            if not _dict_files:
                print(
                    "build-index-only: no my_dict_*.json[.gz] found in "
                    f"{dict_folder} -> nothing to emit dictless tiers from.",
                    flush=True,
                )
            print(
                f"Emitting dictless tiers for {len(_dict_files)} chromosome(s) from "
                f"{os.path.basename(dict_folder)} (databases: "
                f"{', '.join(db_to_samplesid)})...",
                flush=True,
            )
            # [tier-parallel] Tier emission is embarrassingly parallel per chromosome
            # (each call writes its own reg_<chrom>/gt_<chrom>). At high-coverage scale
            # the per-chrom dict is multi-GB and peaks at ~150-180GB RSS, so bound
            # concurrency by RAM: CRISPRME_TIER_WORKERS (default 4 -> ~700GB peak).
            # Biggest dicts first so peak-RAM chromosomes don't pile up at the tail.
            # emit_dictless_tiers_guarded never raises and returns a picklable metadata
            # dict (or None), so it is safe to fan out over a fork pool.
            def _tier_chrom(_dfile):
                _stem = os.path.basename(_dfile)[len("my_dict_"):]
                for _suf in (".json.gz", ".json"):
                    if _stem.endswith(_suf):
                        return _stem[: -len(_suf)]
                return _stem
            # Skip chromosomes already emitted so a restart resumes instead of
            # re-reading multi-GB dicts. Key on the .idx sidecars: each writer emits
            # the .bin, then writes its .idx LAST and ATOMICALLY (tmp+os.replace), so a
            # present .idx is the completion marker. The .bin alone is NOT -- a kill
            # mid-write leaves a truncated .bin that a presence check would wrongly
            # accept as done (silently shipping a short tier store).
            _reg_dir = os.path.join(os.path.dirname(dict_folder), f"registry_{vcf_name}")
            _gt_dir = os.path.join(os.path.dirname(dict_folder), f"genotypes_{vcf_name}")
            _tier_jobs = []
            _tier_present = 0
            for _dfile in sorted(_dict_files, key=os.path.getsize, reverse=True):
                _c = _tier_chrom(_dfile)
                if (os.path.isfile(os.path.join(_reg_dir, f"reg_{_c}.idx"))
                        and os.path.isfile(os.path.join(_gt_dir, f"gt_{_c}.idx"))):
                    _tier_present += 1
                    continue
                _tier_jobs.append((_dfile, db_to_samplesid, _c, dict_folder))
            if _tier_present:
                print(f"  {_tier_present} chromosome(s) already emitted; resuming "
                      f"{len(_tier_jobs)} remaining", flush=True)
            try:
                _tier_workers = int(os.environ.get("CRISPRME_TIER_WORKERS", "4"))
            except ValueError:
                _tier_workers = 4
            _tier_workers = max(1, min(_tier_workers, len(_tier_jobs) or 1))
            _emitted = _tier_present
            if _tier_workers > 1:
                print(f"  emitting tiers with {_tier_workers} parallel worker(s) "
                      "(RAM-bounded; set CRISPRME_TIER_WORKERS to tune)", flush=True)
                import multiprocessing as _mp
                try:
                    with _mp.get_context("fork").Pool(_tier_workers) as _tpool:
                        for _r in _tpool.starmap(
                            _bdt.emit_dictless_tiers_guarded, _tier_jobs
                        ):
                            if _r is not None:
                                _emitted += 1
                except Exception as _tperr:  # noqa: BLE001 - degrade to sequential
                    print(f"WARNING [tier-parallel]: pool failed ({_tperr}); "
                          "re-emitting sequentially.", flush=True)
                    _emitted = 0
                    for _job in _tier_jobs:
                        if _bdt.emit_dictless_tiers_guarded(*_job) is not None:
                            _emitted += 1
            else:
                for _job in _tier_jobs:
                    if _bdt.emit_dictless_tiers_guarded(*_job) is not None:
                        _emitted += 1
            print(
                f"Dictless tier emission complete: {_emitted}/{len(_dict_files)} "
                "chromosome(s) emitted (registry_<vcf>/ + genotypes_<vcf>/ siblings "
                "of the dicts).",
                flush=True,
            )
            # Persist a tiny variant-count manifest so reports can state the
            # database size (SNP variants + genotyped panel sizes) with ZERO VCF
            # access -- summed from the per-chrom registry .idx headers just
            # written. Purely additive; a failure here never affects the build.
            if _emitted:
                try:
                    import json as _json
                    _reg_dir = os.path.join(
                        os.path.dirname(dict_folder), f"registry_{vcf_name}"
                    )
                    _n_records, _dbs = 0, {}
                    for _ix in sorted(_glob(os.path.join(_reg_dir, "reg_*.idx"))):
                        try:
                            with open(_ix) as _fh:
                                _m = _json.load(_fh)
                        except (OSError, ValueError):
                            continue
                        _n_records += int(_m.get("n_records", 0) or 0)
                        if not _dbs and _m.get("databases"):
                            _dbs = _m["databases"]
                    # exact indel count = data lines (minus 1-line header) across the
                    # log_indels log*.txt[.gz] files (one line per indel allele). The
                    # SNP registry does NOT hold indels, so this is the only place the
                    # indel total is captured for the report.
                    _n_indels = None
                    try:
                        import gzip as _gzip
                        _logs = (sorted(_glob(os.path.join(indel_dict, "log*.txt.gz")))
                                 or sorted(_glob(os.path.join(indel_dict, "log*.txt"))))
                        if _logs:
                            _ni = 0
                            for _lg in _logs:
                                _op = (_gzip.open(_lg, "rb") if _lg.endswith(".gz")
                                       else open(_lg, "rb"))
                                _c = 0
                                with _op as _lfh:
                                    while True:
                                        _b = _lfh.read(1 << 20)
                                        if not _b:
                                            break
                                        _c += _b.count(b"\n")
                                _ni += max(0, _c - 1)  # minus header line
                            _n_indels = _ni
                    except Exception:  # noqa: BLE001 - indel count is optional
                        _n_indels = None
                    if _n_records > 0:
                        _payload = {"n_records": _n_records, "databases": _dbs}
                        if _n_indels is not None:
                            _payload["n_indels"] = _n_indels
                        with open(
                            os.path.join(_reg_dir, "variant_count.json"), "w"
                        ) as _out:
                            _json.dump(_payload, _out, indent=2, sort_keys=True)
                        print(
                            f"Wrote variant-count manifest: {_n_records:,} SNPs"
                            + (f" + {_n_indels:,} indels." if _n_indels is not None
                               else " (indel count unavailable)."),
                            flush=True,
                        )
                except Exception as _vc_err:  # noqa: BLE001 - manifest is optional
                    print(
                        f"WARNING [build-index-only]: variant-count manifest not "
                        f"written ({_vc_err}); reports fall back to the .idx headers.",
                        flush=True,
                    )
    # [indel-snp] STEP 1c (ADDITIVE, gated): emit the PHASED indel genotype store so
    # the indel post-analysis can do CONFIRMED-cis SNP+indel co-occurrence. One
    # gt_indel_<chrom>.tsv.gz per chromosome under indel_genotypes_<vcf>/ (sibling of
    # the SNP tiers). GUARDED (stdout warning, never abort); absent -> the post-
    # analysis degrades to PUTATIVE via the log's unphased carriers. Enabled by
    # default (opt-out): set CRISPRME_INDEL_SNP=0 to disable.
    if os.environ.get("CRISPRME_INDEL_SNP", "1") in ("1", "true", "True", "yes"):
        try:
            import build_indel_genotypes as _big
            _igt_dir = os.path.join(workdir, "Dictionaries", f"indel_genotypes_{vcf_name}")
            os.makedirs(_igt_dir, exist_ok=True)
            _keep = None
            if db_to_samplesid:  # union of the per-db panels (match the SNP tiers)
                _keep = set()
                for _sid in db_to_samplesid.values():
                    if _sid and os.path.isfile(_sid):
                        _keep |= _big._load_panel(_sid)
                _keep = _keep or None
            _n_ig = 0
            # Build (vcf, out, keep) jobs; detect chrom from the first data line (the
            # filename may not carry it). Skip chromosomes already built so a restart
            # resumes instead of re-streaming 20GB VCFs.
            _ig_jobs = []
            for _vf in sorted(_glob(os.path.join(vcfdir, "*.vcf.gz"))
                              + _glob(os.path.join(vcfdir, "*.vcf"))):
                _chrom = None
                _op = gzip.open(_vf, "rt") if _vf.endswith(".gz") else open(_vf)
                with _op as _fh:
                    for _ln in _fh:
                        if _ln.startswith("#"):
                            continue
                        _cc = _ln.split("\t", 1)[0]
                        _chrom = _cc if _cc.startswith("chr") else "chr" + _cc
                        break
                if _chrom is None:
                    continue
                _igt_out = os.path.join(_igt_dir, f"gt_indel_{_chrom}.tsv.gz")
                # Resume only if the store is a COMPLETE, non-truncated gzip -- a
                # size>0 check would skip a file left partial by an OOM/SIGKILL
                # mid-write and silently ship a short store (the gt_indel_chr1 bug).
                if os.path.isfile(_igt_out) and _big.store_is_complete(_igt_out):
                    _n_ig += 1  # already built (resume)
                    continue
                _ig_jobs.append((_vf, _igt_out, _keep))
            # compile_indel_genotypes streams the VCF (low RAM, I/O + gzip bound), so
            # fan out wide -- CRISPRME_INDELGT_WORKERS (default 8). Guarded worker so a
            # single bad chromosome cannot sink the store.
            if _ig_jobs:
                try:
                    _igw = int(os.environ.get("CRISPRME_INDELGT_WORKERS", "8"))
                except ValueError:
                    _igw = 8
                _igw = max(1, min(_igw, len(_ig_jobs)))
                print(f"[indel-snp] building phased indel GT for {len(_ig_jobs)} "
                      f"chromosome(s) with {_igw} parallel worker(s)", flush=True)
                if _igw > 1:
                    import multiprocessing as _mp
                    with _mp.get_context("fork").Pool(_igw) as _igpool:
                        _ig_res = _igpool.starmap(
                            _big.compile_indel_genotypes_safe, _ig_jobs)
                    _n_ig += sum(1 for _r in _ig_res if _r is not None and _r >= 0)
                else:
                    for _job in _ig_jobs:
                        if _big.compile_indel_genotypes_safe(*_job) >= 0:
                            _n_ig += 1
            print(f"[indel-snp] phased indel genotype store: {_n_ig} chromosome(s) "
                  f"-> {os.path.basename(_igt_dir)}", flush=True)
        except Exception as _ig_err:  # noqa: BLE001 - store is optional (-> PUTATIVE)
            print(f"WARNING [indel-snp]: phased indel genotype store not built "
                  f"({_ig_err}); indel post-analysis falls back to PUTATIVE.", flush=True)
    # STEP 2: index the enriched (SNP) genome
    if not os.path.isdir(snp_idx):
        print(f"Building variant index {os.path.basename(snp_idx)}...", flush=True)
        code = subprocess.call(
            f"crispritz.py index-genome {genome_ref}+{vcf_name} {enriched}/ "
            f"{pamfile} -bMax {bMax + 1} -th {thread}",
            shell=True,
            cwd=workdir,
        )
        if code != 0 or not os.path.isdir(snp_idx):
            error(f"Variant genome indexing failed (expected {snp_idx})")
    else:
        print(f"Variant index already present: {snp_idx}", flush=True)
    # STEP 3: index the indels genome (pool_index_indels.py, as the pipeline does)
    if not os.path.isdir(indels_idx):
        print(f"Building indels index {os.path.basename(indels_idx)}...", flush=True)
        # [indel-snp] Overlay SNP IUPAC codes onto the fake-indel flanks BEFORE
        # indexing so this _INDELS index also finds SNP+indel co-occurring
        # off-targets (searchTST matches IUPAC codes for free -- no -var needed).
        # Enabled by default (opt-out): set CRISPRME_INDEL_SNP=0 to disable. GUARDED
        # (stdout warning, never abort -- stderr is fatal here), mirroring the tiers above.
        if os.environ.get("CRISPRME_INDEL_SNP", "1") in ("1", "true", "True", "yes"):
            try:
                import overlay_indel_snps as _ois
                _n = _ois.main(["overlay_indel_snps", indels_out, enriched, indel_dict])
                print(
                    f"[indel-snp] SNP-overlaid the fake-indel genome "
                    f"({_n} flank bases IUPAC-coded).",
                    flush=True,
                )
            except Exception as _ois_err:
                print(
                    f"WARNING [indel-snp]: fake-indel SNP overlay failed "
                    f"({_ois_err}); this _INDELS index will be SNP-blind.",
                    flush=True,
                )
        pool = os.path.join(script_path, "pool_index_indels.py")
        code = subprocess.call(
            [
                sys.executable,
                pool,
                f"{indels_out}/",
                pamfile,
                pam_char,
                genome_ref,
                vcf_name,
                str(bMax + 1),
                str(thread),
            ],
            cwd=workdir,
        )
        if code != 0 or not os.path.isdir(indels_idx):
            error(f"Indels indexing failed (expected {indels_idx})")
    else:
        print(f"Indels index already present: {indels_idx}", flush=True)
    _write_pam_build(snp_idx, pamfile)
    _write_display_label(snp_idx, display_name)
    # ADDITIVE: emit the combined samplesIDs/<vcf_name>.samplesID.txt for a MERGED
    # panel so a user who BUILDS + PUBLISHES their own merged index is self-complete
    # (the search's --samplesID listing wants that single combined file). Reuses the
    # already-computed db_to_samplesid map and the ONE union implementation; a strict
    # NO-OP for a single-dataset panel / when --samplesID was absent (empty map).
    # Fully guarded (STDOUT-only; never aborts) — the index/dicts are already built.
    _emit_combined_samplesid(workdir, vcf_name, genome_ref, db_to_samplesid)
    print(
        f"Variant index ready: {os.path.basename(snp_idx)} (+ _INDELS). A later "
        f"variant-aware search on {genome_ref} with {vcf_name} reuses it.",
        flush=True,
    )


def print_help_download() -> None:
    """Prints detailed help information for the download functionality."""
    sys.stderr.write(
        "The download functionality fetches CRISPRme reference data from a "
        "HuggingFace dataset repository over HF's CDN — typically much faster and "
        "more reliable than the original UCSC/FTP sources — and places it in the "
        "canonical CRISPRme directory layout (Genomes/, Annotations/, PAMs/, "
        "samplesIDs/, VCFs/, genome_library/). FASTA is decompressed after "
        f"download; VCFs are kept bgzipped. Default repo: {DEFAULT_HF_REPO} "
        "(override with --hf-repo or the CRISPRME_HF_REPO environment "
        "variable).\n"
    )
    sys.stderr.write(
        "Options:\n"
        "\t--what, component to fetch: genome | annotations | pams | samples | "
        "vcf | index | all (all = genome+annotations+pams+samples, plus the "
        "default hg38 reference index so a reference scan works out of the box) "
        "[REQUIRED]\n"
        "\t--ref, reference genome name for --what genome/index [default: hg38]\n"
        "\t--dataset, variant dataset name for --what vcf (e.g. 1000G, HGDP) "
        "[REQUIRED for --what vcf]\n"
        "\t--index-name, precomputed index directory name for --what index "
        "(e.g. the default NRG_3_hg38, or the variant-aware NRG_3_hg38+hg38_1000G_HGDP) "
        "[REQUIRED for --what index]\n"
        "\t--hf-repo, HuggingFace dataset repo id to fetch from [OPTIONAL]\n"
        "\t--source, for --what genome: hf (HuggingFace, default) | ucsc "
        "(UCSC goldenPath by assembly name, e.g. --ref susScr11) | url "
        "(explicit --url link) [OPTIONAL]\n"
        "\t--url, explicit genome download URL when --source url [OPTIONAL]\n"
        "\t--no-genotypes, for --what index: skip fetching the separate Tier-1 "
        "genotype store (genotypes_<vcf>.tar.gz). Off-target detection still "
        "works via the bundled Tier-0 registry, but per-sample Samples are "
        "degraded until the store is present [OPTIONAL, default: fetch it]\n"
        "\t--path, working directory the CRISPRme dir-tree lives under "
        "[OPTIONAL, default: current directory]\n"
    )
    sys.exit(1)


# Default REFERENCE index published on HF for hg38 (the default NRG PAM, bMax 3);
# fetched by ``download --what all`` so a simple reference scan works out of the box.
DEFAULT_REFERENCE_INDEX_HG38 = "NRG_3_hg38"


def download_data() -> None:
    """Fetches CRISPRme reference components from a HuggingFace dataset repo.

    A fast-CDN alternative to the FTP/UCSC downloads performed by ``setup``:
    genome, annotations, PAMs, sample-ID files, variant VCFs and precomputed
    indexes are pulled from a HuggingFace dataset into the canonical CRISPRme
    layout under the working directory.
    """
    args = input_args[2:]  # retrieve download input arguments
    if "--help" in args or "--what" not in args:
        print_help_download()
    what = args[args.index("--what") + 1]
    repo = args[args.index("--hf-repo") + 1] if "--hf-repo" in args else None
    ref = args[args.index("--ref") + 1] if "--ref" in args else "hg38"
    dataset = args[args.index("--dataset") + 1] if "--dataset" in args else None
    index_name = args[args.index("--index-name") + 1] if "--index-name" in args else None
    # genome source: 'hf' (HuggingFace, default), 'ucsc' (goldenPath by assembly),
    # or 'url' (explicit download link). Only meaningful for --what genome.
    source = args[args.index("--source") + 1] if "--source" in args else "hf"
    url = args[args.index("--url") + 1] if "--url" in args else None
    # ADDITIVE: for --what index, --no-genotypes skips the big Tier-1 genotype
    # store (detection-only; degraded Samples). Default = fetch it.
    genotypes = "--no-genotypes" not in args
    workdir = os.getcwd()
    if "--path" in args:
        workdir = os.path.abspath(args[args.index("--path") + 1])
        if not os.path.isdir(workdir):
            error(f"The working directory {workdir} does not exist")
    # "all" fetches the always-needed reference bundle (variant VCFs are dataset
    # specific, so they are requested explicitly). On hg38 it ALSO fetches the
    # default-PAM REFERENCE index (below) so a simple reference-genome scan works
    # out of the box.
    if what == "all":
        components = ["genome", "annotations", "pams", "samples"]
    else:
        components = [what]
    for comp in components:
        try:
            if comp == "genome" and source in ("ucsc", "url"):
                # non-HuggingFace reference genome (e.g. a UCSC assembly such as
                # the pig susScr11) -> download straight into Genomes/<ref>/
                dest = download_reference_genome(
                    ref, os.path.join(workdir, "Genomes"), source=source, url=url
                )
            else:
                dest = download_component(
                    comp,
                    workdir,
                    repo=repo,
                    ref=ref,
                    dataset=dataset,
                    index_name=index_name,
                    genotypes=genotypes,
                )
        except (ValueError, ImportError, RuntimeError) as e:
            error(str(e))
        print(f"Downloaded {comp} -> {dest}", flush=True)

    # For --what all on hg38, also fetch the default-PAM REFERENCE index so a
    # simple reference-genome scan works immediately. Without it, the first
    # reference search must build the index on demand from the raw genome
    # (correct, but slow). Non-fatal: on any hiccup the genome is still present
    # and the index is built on demand at search time.
    if what == "all" and ref == "hg38":
        ref_index = index_name or DEFAULT_REFERENCE_INDEX_HG38
        try:
            dest = download_component(
                "index", workdir, repo=repo, ref=ref,
                index_name=ref_index, genotypes=False,
            )
            print(f"Downloaded reference index {ref_index} -> {dest}", flush=True)
        except (ValueError, ImportError, RuntimeError) as e:
            sys.stderr.write(
                f"WARNING: could not fetch the default reference index "
                f"{ref_index} ({e}). The reference genome is present; a reference "
                f"index will be built on demand at the first search.\n"
            )


def print_help_publish_index() -> None:
    """Prints detailed help information for the publish-index functionality."""
    sys.stderr.write(
        "The publish-index functionality uploads a locally built CRISPRitz "
        "reference index (a genome_library/<name> directory, e.g. from "
        "build-index-only) to a HuggingFace dataset repository as a single "
        "compressed archive, so it can later be fetched with "
        "'crisprme.py download --what index'. An HF write token is required "
        "(provide --token or set HF_TOKEN; never commit it).\n"
    )
    sys.stderr.write(
        "Options:\n"
        "\t--index, path to the genome_library/<name> index directory to publish "
        "[REQUIRED]\n"
        "\t--hf-repo, HuggingFace dataset repo id to upload to [OPTIONAL]\n"
        "\t--token, HuggingFace write token [OPTIONAL if HF_TOKEN is set]\n"
        "\t--name, optional human-friendly display name for the index (else a "
        "clear convention label is derived from the folder name) [OPTIONAL]\n"
        "\t--dictless, EXCLUDE the per-sample SNP dictionaries "
        "(dictionaries_<vcf>/) from the main tarball — the additive Tier-0 "
        "registry + Tier-1 genotype tiers replace them — while keeping the indel "
        "logs. Without this flag the classic dicts are bundled exactly as before "
        "[OPTIONAL, default: off]. In BOTH modes the small registry_<vcf>/ is "
        "added when present, and a separate genotypes_<vcf>.tar.gz companion is "
        "uploaded when a genotype store exists.\n"
    )
    sys.exit(1)


def publish_index_cmd() -> None:
    """Uploads a locally built reference index to a HuggingFace dataset repo."""
    args = input_args[2:]  # retrieve publish-index input arguments
    if "--help" in args or "--index" not in args:
        print_help_publish_index()
    index_dir = os.path.abspath(args[args.index("--index") + 1])
    repo = args[args.index("--hf-repo") + 1] if "--hf-repo" in args else None
    token = args[args.index("--token") + 1] if "--token" in args else None
    # optional human-friendly display name (else the UI parses a convention label)
    name = args[args.index("--name") + 1] if "--name" in args else None
    # ADDITIVE: --dictless drops the 152GB per-sample SNP dicts from the main
    # tarball (tiers replace them); default off => byte-for-byte-unchanged publish.
    dictless = "--dictless" in args
    try:
        remote_path = publish_index(
            index_dir, repo=repo, token=token, display_name=name, dictless=dictless
        )
    except (ValueError, ImportError) as e:
        error(str(e))
    print(f"Published index to {resolve_repo(repo)}:{remote_path}", flush=True)


def print_help_assembly_search() -> None:
    """Prints detailed help information for the assembly-search functionality.

    Outputs a description of the pipeline and lists all available command-line
    options to stderr, then exits the program.
    """
    sys.stderr.write(
        "The assembly-search functionality searches a fully assembled personal "
        "diploid genome directly -- two haplotype assemblies (e.g. paternal and "
        "maternal) -- instead of inferring variants from population data via a "
        "reference genome + VCF. No --vcf is used: each haplotype assembly IS the "
        "individual's genome. Each haplotype is searched independently with the "
        "same underlying complete-search pipeline, then predictions are lifted "
        "to hg38 (via the supplied chain files) and reconciled: found on both "
        "haplotypes is homozygous-equivalent, found on only one is "
        "heterozygous-equivalent, and predictions with no hg38 equivalent at all "
        "are haplotype-non-mappable -- invisible to any reference-based search.\n"
    )
    sys.stderr.write(
        "Options:\n"
        "\t--genome-paternal, --genome-maternal, one haplotype's per-chromosome "
        "FASTA folder each [REQUIRED]\n"
        "\t--chain-paternal, --chain-maternal, one haplotype's liftOver chain "
        "file vs. GRCh38 each [REQUIRED]\n"
        "\t--chrom-alias-paternal, --chrom-alias-maternal, one haplotype's "
        "chromAlias file each (tab-separated, columns '# assembly', 'ucsc', "
        "'genbank' -- HPRC-style) [REQUIRED]\n"
        "\t--guide, specify a file containing guide RNAs [REQUIRED]\n"
        "\t--pam, specify a file containing the PAM sequence [REQUIRED]\n"
        "\t--mm, number of mismatches allowed in the search [REQUIRED]\n"
        "\t--bDNA, number of DNA bulges allowed in the search [OPTIONAL]\n"
        "\t--bRNA, number of RNA bulges allowed in the search [OPTIONAL]\n"
        "\t--merge, window size (nucleotides) to merge candidate off-targets "
        "using the highest scoring as pivot [default: 3] -- also used as the "
        "locus-clustering threshold when reconciling the two haplotypes, so it "
        "must describe both runs consistently\n"
        "\t--output, base output name; each haplotype's results are saved in "
        "Results/<name>_paternal and Results/<name>_maternal, and the "
        "reconciled combined report in Results/<name>_combined [REQUIRED]\n"
        "\t--thread, set number of threads to use [default: 8]\n"
        "\t--debug, debug mode (passed through to each haplotype's search)\n"
    )
    sys.exit(1)


def _check_named_dir(args: List[str], flag: str, description: str) -> str:
    """Generic version of `_check_genome` for a caller-specified flag name."""
    try:
        d = os.path.abspath(args[args.index(flag) + 1])
    except IndexError:
        error(f"Missing input for {flag}. {description} must be specified")
    if not os.path.isdir(d):
        error(f"The folder specified for {flag} does not exist")
    return d


def _check_named_file(args: List[str], flag: str, description: str) -> str:
    """Generic version of the existing file-checking helpers for a
    caller-specified flag name."""
    try:
        f = os.path.abspath(args[args.index(flag) + 1])
    except IndexError:
        error(f"Missing input for {flag}. {description} must be specified")
    if not os.path.isfile(f):
        error(f"The file specified for {flag} does not exist")
    return f


def _check_mandatory_args_assembly_search(args: List[str]) -> None:
    required = [
        "--genome-paternal", "--genome-maternal",
        "--chain-paternal", "--chain-maternal",
        "--chrom-alias-paternal", "--chrom-alias-maternal",
        "--guide", "--pam", "--mm", "--output",
    ]
    for flag in required:
        if flag not in args:
            error(f"{flag} is required")


def _run_haplotype_search(
    genomedir: str, guidefile: str, pamfile: str, mm: int, bDNA: int, bRNA: int,
    merge_t: int, output_name: str, thread: int, debug: bool,
) -> str:
    """Runs `complete-search` for one haplotype as a subprocess, mirroring the
    pattern already established by `complete_test_crisprme()` -- a subcommand
    invoking `complete-search` as a fresh process rather than calling
    `complete_search()` in-process, which is necessary here since
    `complete_search()` derives its arguments from the module-level
    `input_args = sys.argv`, not from parameters, so it can't cleanly be
    called twice in-process with different `--genome` values. Uses the exact
    interpreter and script path running right now (rather than relying on a
    bare `crisprme.py` being on PATH, as `complete_test.py` does) since this
    is a new, not-yet-installed subcommand. Output streams straight to the
    terminal rather than being captured, also matching
    `complete_test_crisprme()` -- the real per-step pipeline logs are
    written to `Results/<output_name>/log_verbose.txt`/`log_error.txt`
    regardless of how this outer process's own stdout/stderr are handled, so
    capturing them here bought no real diagnostic value while hiding all
    progress output during what can be a multi-hour run.

    Returns:
        The absolute path to the haplotype's `complete-search` output folder.
    """
    python_exe = sys.executable
    crisprme_script = os.path.abspath(__file__)
    debug_flag = "--debug" if debug else ""
    cmd = (
        f"{python_exe} {crisprme_script} complete-search --genome {genomedir} "
        f"--guide {guidefile} --pam {pamfile} --mm {mm} --bDNA {bDNA} --bRNA {bRNA} "
        f"--merge {merge_t} --output {output_name} --thread {thread} {debug_flag}"
    )
    output_folder = os.path.join(current_working_directory, CRISPRMEDIRS[1], output_name)
    code = subprocess.call(cmd, shell=True, cwd=current_working_directory)
    if code != 0:
        raise OSError(
            f"\nHaplotype search failed for --output {output_name}! See "
            f"{os.path.join(output_folder, 'log_error.txt')} for details\n"
        )
    return output_folder


def assembly_search() -> None:
    """Searches a personal diploid genome assembly (two haplotypes) for
    off-targets and reconciles predictions across both, mapped to hg38.

    No `--vcf` is used or accepted: each haplotype assembly already is the
    individual's genome, so there's nothing to infer from population data.
    Each haplotype is searched independently via `complete-search`
    (`_run_haplotype_search`), then predictions are reconciled via
    `assembly_reconcile.reconcile_haplotypes` -- see that module for the
    reconciliation algorithm and the real-data validation behind it.
    """
    args = input_args[2:]
    if "--help" in args or not args:
        print_help_assembly_search()
    check_crisprme_dirtree()
    _check_mandatory_args_assembly_search(args)
    check_liftover_available()

    genome_paternal = _check_named_dir(args, "--genome-paternal", "Paternal genome folder")
    genome_maternal = _check_named_dir(args, "--genome-maternal", "Maternal genome folder")
    chain_paternal = _check_named_file(args, "--chain-paternal", "Paternal liftOver chain file")
    chain_maternal = _check_named_file(args, "--chain-maternal", "Maternal liftOver chain file")
    chrom_alias_paternal = _check_named_file(args, "--chrom-alias-paternal", "Paternal chromAlias file")
    chrom_alias_maternal = _check_named_file(args, "--chrom-alias-maternal", "Maternal chromAlias file")
    guidefile = _check_guide(args, True)
    pamfile = _check_pam(args)
    mm = _check_mm(args)
    bDNA = _check_bdna(args, "--bDNA" in args)
    bRNA = _check_brna(args, "--bRNA" in args)
    merge_t = _check_merge(args, "--merge" in args)
    thread = _check_threads(args, "--thread" in args)
    debug = "--debug" in args

    try:
        output_base = args[args.index("--output") + 1]
    except IndexError:
        error("Missing input for --output. Output base name must be specified")

    combined_output = os.path.join(current_working_directory, CRISPRMEDIRS[1], f"{output_base}_combined")
    combined_tsv = os.path.join(combined_output, f"{output_base}_combined_hg38.tsv")
    # only guard against clobbering a previously *completed* run -- leftover
    # reconciliation intermediates (BED files) from a prior failed attempt
    # are fine to overwrite and are regenerated fresh below regardless
    if os.path.isfile(combined_tsv):
        error(
            f"{combined_tsv} already exists from a previously completed "
            "assembly-search run! Select another --output name, or delete "
            "it to re-run reconciliation."
        )
    os.makedirs(combined_output, exist_ok=True)

    paternal_output_name = f"{output_base}_paternal"
    maternal_output_name = f"{output_base}_maternal"
    paternal_output_dir = os.path.join(current_working_directory, CRISPRMEDIRS[1], paternal_output_name)
    maternal_output_dir = os.path.join(current_working_directory, CRISPRMEDIRS[1], maternal_output_name)

    # avoid re-running a multi-hour haplotype search if it already completed
    # successfully -- e.g. a retry after a failure that only affected
    # reconciliation (a chromAlias/liftOver issue, say), not the searches
    paternal_complete = haplotype_search_complete(paternal_output_dir)
    paternal_reusable = paternal_complete and haplotype_params_match(
        paternal_output_dir, genome_paternal, guidefile, pamfile, mm, bDNA, bRNA, merge_t
    )
    if paternal_reusable:
        print(f"Found existing completed paternal search results at Results/{paternal_output_name}, reusing (not re-running)")
        paternal_results = paternal_output_dir
    else:
        if paternal_complete:
            clean_incomplete_haplotype_output(
                paternal_output_dir, reason="results built with different search parameters"
            )
        else:
            clean_incomplete_haplotype_output(paternal_output_dir)
        print(f"Running paternal haplotype search -> Results/{paternal_output_name}")
        paternal_results = _run_haplotype_search(
            genome_paternal, guidefile, pamfile, mm, bDNA, bRNA, merge_t,
            paternal_output_name, thread, debug,
        )

    maternal_complete = haplotype_search_complete(maternal_output_dir)
    maternal_reusable = maternal_complete and haplotype_params_match(
        maternal_output_dir, genome_maternal, guidefile, pamfile, mm, bDNA, bRNA, merge_t
    )
    if maternal_reusable:
        print(f"Found existing completed maternal search results at Results/{maternal_output_name}, reusing (not re-running)")
        maternal_results = maternal_output_dir
    else:
        if maternal_complete:
            clean_incomplete_haplotype_output(
                maternal_output_dir, reason="results built with different search parameters"
            )
        else:
            clean_incomplete_haplotype_output(maternal_output_dir)
        print(f"Running maternal haplotype search -> Results/{maternal_output_name}")
        maternal_results = _run_haplotype_search(
            genome_maternal, guidefile, pamfile, mm, bDNA, bRNA, merge_t,
            maternal_output_name, thread, debug,
        )

    print("Reconciling paternal and maternal predictions against hg38...")
    haplotypes = {
        "paternal": {
            "chrom_alias_file": chrom_alias_paternal,
            "chain_file": chain_paternal,
            "results_dir": paternal_results,
        },
        "maternal": {
            "chrom_alias_file": chrom_alias_maternal,
            "chain_file": chain_maternal,
            "results_dir": maternal_results,
        },
    }
    combined, summary = reconcile_haplotypes(haplotypes, combined_output, merge_bp=merge_t)
    combined.to_csv(combined_tsv, sep="\t", index=False)

    print(f"Reconciliation complete. Wrote {combined_tsv}")
    for category, count in summary.items():
        print(f"  {category}: {count}")


def target_integration():
    if "--help" in input_args:
        print(
            "This is the automated integration process that process the final result file to generate a usable target panel."
        )
        print("These are the flags that must be used in order to run this function:")
        print(
            "\t--targets, used to specify the final result file to use in the panel creation process"
        )
        print(
            "\t--empirical_data, used to specify the file that contains empirical data provided by the user to assess in-silico targets"
        )
        print("\t--output, used to specify the output folder for the results")
        exit(0)

    if "--targets" not in input_args:
        print("--targets must be contained in the input")
        exit(1)
    else:
        try:
            target_file = os.path.abspath(input_args[input_args.index("--targets") + 1])
        except IndexError:
            print("Please input some parameter for flag --targets")
            exit(1)
        if not os.path.isfile(target_file):
            print("The file specified for --target_file does not exist")
            exit(1)

    # if "--vcf_dir" not in input_args:
    #     print("--vcf_dir non in input, multi-variant haplotype will not be calculated")
    #     vcf_dir = script_path+'vuota/'
    #     # exit(1)
    # else:
    #     try:
    #         vcf_dir = os.path.abspath(
    #             input_args[input_args.index("--vcf_dir")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --vcf_dir")
    #         exit(1)
    #     if not os.path.isdir(vcf_dir):
    #         print("The folder specified for --vcf_dir does not exist")
    #         exit(1)

    # if "--genome_version" not in input_args:
    #     print("--genome_version must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         genome_version = input_args[input_args.index(
    #             "--genome_version")+1]
    #     except IndexError:
    #         print("Please input some parameter for flag --genome")
    #         exit(1)

    # if "--guide" not in input_args:
    #     guidefile = script_path+'vuoto.txt'
    #     # print("--guide must be contained in the input")
    #     # exit(1)
    # else:
    #     try:
    #         guidefile = os.path.abspath(
    #             input_args[input_args.index("--guide")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --guide")
    #         exit(1)
    #     if not os.path.isfile(guidefile):
    #         print("The file specified for --guide does not exist")
    #         exit(1)

    if "--empirical_data" not in input_args:
        print("--empirical_data not in input, proceeding without empirical data")
        empiricalfile = script_path + "vuoto.txt"
        # exit(1)
    else:
        try:
            empiricalfile = os.path.abspath(
                input_args[input_args.index("--empirical_data") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --empirical_data")
            exit(1)
        if not os.path.isfile(empiricalfile):
            print("The file specified for --empirical_data does not exist")
            exit(1)

    # if "--gencode" not in input_args:
    #     print("--gencode must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         gencode_file = os.path.abspath(
    #             input_args[input_args.index("--gencode")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --gencode")
    #         exit(1)
    #     if not os.path.isfile(gencode_file):
    #         print("The file specified for --gencode does not exist")
    #         exit(1)

    if "--output" not in input_args:
        print("--output must be contained in the input")
        exit(1)
    else:
        try:
            outputfolder = os.path.abspath(input_args[input_args.index("--output") + 1])
        except IndexError:
            print("Please input some parameter for flag --output")
            exit(1)
        if not os.path.isdir(outputfolder):
            print("The folder specified for --output does not exist")
            exit(1)

    os.system(
        f"{script_path}./empirical_integrator.py {target_file} {empiricalfile} {outputfolder}"
    )


def print_help_gnomad_converter():
    """
    Prints the help information for the gnomAD converter functionality, providing
    details on the conversion process from gnomAD VCFs to VCFs compatible with
    CRISPRme. It outlines the options available for specifying directories,
    sample IDs, variant filtering, multiallelic site handling, and thread usage
    during the conversion process.

    Raises:
        SystemExit: If the help information is displayed to guide users on using
        the gnomAD converter functionality.
    """

    # functionality description
    sys.stderr.write(
        "The gnomAD converter functionality simplifies the conversion process "
        "of gnomAD VCFs (versions 3.1 and 4.0) into VCFs supported by CRISPRme. "
        "It ensures a seamless transition while maintaining compatibility with "
        "CRISPRme's requirements, focusing on the structure and content of "
        "precomputed sample IDs file \n\n"
    )
    # options
    sys.stderr.write(
        "Options:\n"
        "\t--gnomAD_VCFdir, specifies the directory containing gnomAD VCFs. "
        "Files must have the BGZ extension\n"
        "\t--samplesID, specifies the precomputed sample IDs file necessary "
        "for incorporating population-specific information into the output "
        "VCFs\n"
        "\t--joint, optional flag to specify the input GnomAD VCF contain joint "
        "allele frequencies\n"
        "\t--keep, optional flag to retain all variants, regardless of their "
        "filter flag. By default, variants with a filter flag different from "
        "PASS are discarded\n"
        "\t--multiallelic, optional flag to merge variants mapped to the "
        "same position, creating multiallelic sites in the output VCFs. By "
        "default, each site remains biallelic\n"
        "\t--threads, used to set the number of threads used in the conversion "
        "process [default 8]\n"
    )
    sys.exit(1)


def gnomAD_converter():
    """
    Runs the gnomAD converter functionality based on specified arguments, converting
    gnomAD VCF files into formats compatible with CRISPRme.

    Raises:
        ValueError: If mandatory arguments are missing or have incorrect values.
        FileExistsError: If the specified gnomAD VCF directory cannot be located.
        FileNotFoundError: If the specified sample IDs file cannot be found.
        subprocess.SubprocessError: If an error occurs during the gnomAD VCF
            conversion process.
    """

    args = input_args[2:]  # recover gnomAD converter args
    if "--help" in args or not args:  # print help
        print_help_gnomad_converter()
        sys.exit(1)
    if "--gnomAD_VCFdir" not in args:
        raise ValueError(
            "--gnomAD_VCFdir is a mandatory argument required for the conversion "
            "process. Please specify the directory containing gnomAD VCFs using "
            "this option\n"
        )
    if "--samplesID" not in args:
        raise ValueError(
            "--samplesID is a mandatory argument required for the conversion "
            "process. Please specify the sample IDs file this option\n"
        )
    # read gnomAD directory arg
    try:
        gnomad_dir = args[args.index("--gnomAD_VCFdir") + 1]
        if gnomad_dir.startswith("--"):
            raise ValueError("Please input some parameter for flag --gnomAD_VCFdir\n")
        gnomad_dir = os.path.abspath(gnomad_dir)  # first sanity check passed
        if not os.path.isdir(gnomad_dir):
            raise FileExistsError(f"Unable to locate {gnomad_dir}")
    except IndexError as e:
        raise ValueError(
            "Please input some parameter for flag --gnomAD_VCFdir\n"
        ) from e
    # read samples ids arg
    try:
        samples_ids = args[args.index("--samplesID") + 1]
        if samples_ids.startswith("--"):
            raise ValueError("Please input some parameter for flag --samplesID")
        samples_ids = os.path.abspath(samples_ids)  # first sanity check passed
        if not os.path.isfile(samples_ids):
            raise FileNotFoundError(f"Unable to locate {samples_ids}")
    except IndexError as e:
        raise ValueError("Please input some parameter for flag --samplesID") from e
    # read joint gnomad vcf files
    joint = "--joint" in args
    # read keep arg
    keep = "--keep" in args  # keep all variants regardless of filter label
    # read multiallelic arg
    multiallelic = "--multiallelic" in args  # merge variants in multiallelic sites
    # read threads arg
    threads = 8
    if "--threads" in args:
        try:
            threads = int(args[args.index("--threads") + 1])
            if threads <= 0:
                raise ValueError(f"Forbidden number of threads ({threads})")
        except IndexError as e:
            raise ValueError("Missing or forbidden threads value") from e
    # read vcf filter pass values arg
    vcf_filter_pass_values = "PASS,."
    if "--vcf-filter-pass-values" in args:
        try:
            vcf_filter_pass_values = args[args.index("--vcf-filter-pass-values") + 1]
            if vcf_filter_pass_values.startswith("--"):
                raise ValueError("Please input a value for flag --vcf-filter-pass-values")
        except IndexError as e:
            raise ValueError("Missing input for --vcf-filter-pass-values") from e
    # run gnom AD converter
    gnomad_converter_script = os.path.join(script_path, "convert_gnomAD_vcfs.py")
    cmd = (
        f"python {gnomad_converter_script} {gnomad_dir} {samples_ids} {joint} "
        f"{keep} {multiallelic} {threads} {vcf_filter_pass_values}"
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise subprocess.SubprocessError(
            f"An error occurred while converting gnomAD VCFs in {gnomad_dir}"
        )


def personal_card():
    if "--help" in input_args:
        print(
            "This is the personal card generator that creates a files with all the private targets for the input sample"
        )
        print("These are the flags that must be used in order to run this function:")
        print(
            "\t--result_dir, directory containing the result from which extract the targets to generate the card"
        )
        print(
            "\t--guide_seq, sequence of the guide to use in order to exctract the targets"
        )
        print("\t--sample_id, ID of the sample to use in order to generate the card")
        exit(0)

    if "--result_dir" not in input_args:
        print("--result_dir not in input, please input a result directory")
        exit(1)
    else:
        try:
            result_dir = os.path.abspath(
                input_args[input_args.index("--result_dir") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --result_dir")
            exit(1)
        if not os.path.isdir(result_dir):
            print("The folder specified for --result_dir does not exist")
            exit(1)

    if "--guide_seq" not in input_args:
        print(
            "--guide_seq must be contained in the input, e.g. CTAACAGTTGCTTTTATCACNNN"
        )
        exit(1)
    else:
        try:
            guide = input_args[input_args.index("--guide_seq") + 1]
        except IndexError:
            print("Please input some parameter for flag --guide_seq")
            exit(1)
    if "--sample_id" not in input_args:
        print("--sample_id must be contained in the input, e.g. HG00001")
        exit(1)
    else:
        try:
            sample_id = input_args[input_args.index("--sample_id") + 1]
        except IndexError:
            print("Please input some parameter for flag --sample_id")
            exit(1)

    os.system(
        script_path
        + "./generate_sample_card.py "
        + result_dir
        + " "
        + guide
        + " "
        + sample_id
        + " "
        + script_path
    )

def print_help_web_interface():
    # functionality description
    sys.stderr.write(
        "This function starts a local server to use the web interface.\n"
        "Open your browser at http://127.0.0.1:8080\n"
    )
    # options
    sys.stderr.write(
        "Options:\n"
        "\t--debug, debug mode\n"
    )
    sys.exit(1)


def web_interface():
    args = input_args[2:]
    if "--help" in args or len(input_args) < 2:  # print help and exit
        print_help_web_interface()
    # resolve index.py relative to this script's location
    # regardless of conda/source install layout
    index_script = os.path.join(corrected_web_path, "index.py")
    if not os.path.isfile(index_script):
        sys.stderr.write(
            f"Error: Cannot find index.py at {index_script}\n"
            "The web interface requires index.py to be co-located with crisprme.py.\n"
        )
        sys.exit(1)
    try: 
        subprocess.run([sys.executable, index_script], check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            f"\nWeb interface exited with error code {e.returncode}.\n"
            "Check above for traceback. Common causes:\n"
            "  - Missing dependencies: pip install dash flask flask-caching "
            "dash-bootstrap-components\n"
            "  - Port 8080 already in use: lsof -i :8080\n"
        )
        sys.exit(e.returncode)
    except FileNotFoundError:
        sys.stderr.write(f"Error: Python interpreter not found at {sys.executable}\n")
        sys.exit(1)
    except KeyboardInterrupt:
        sys.stderr.write("\nWeb interface stopped\n")
        sys.exit(0)


def crisprme_version():
    if len(input_args) != 2:
        sys.stderr.write("Wrong number of arguments for crisprme.py version\n")
        sys.exit(1)
    sys.stdout.write(f"v{__version__}\n")


def print_help_generate_report() -> None:
    """Prints detailed help for the generate-report functionality."""
    sys.stderr.write(
        "The generate-report functionality builds a SELF-CONTAINED, easily "
        "shareable report for a completed CRISPRme run. It produces a single "
        "ZIP (<jobid>_report.zip) bundling report.html (offline, no external "
        "dependencies -- plots embedded as base64 PNG, the top-1000 off-target "
        "table inline, CSS inline), integrated_results.tsv.gz (the full results "
        "the HTML links to with a relative href), the top-1000 and top-100 "
        "validation-panel tables, and the non-empty per-tier curated TSVs. It is "
        "a portable digest of the full interactive "
        "website result, aimed at sharing off-target predictions (e.g. to "
        "design a targeted-NGS / rhAMP-Seq confirmation panel).\n"
    )
    sys.stderr.write(
        "Options:\n"
        "\t--result-dir, a CRISPRme result folder (Results/<jobid>). The "
        "integrated_results TSV, .Params.txt, .version.txt and the default "
        "output location are resolved from here [REQUIRED unless "
        "--integrated-results is given]\n"
        "\t--integrated-results, an explicit integrated_results TSV (.tsv or "
        ".tsv.gz) to report on [OPTIONAL]\n"
        "\t--samplesID-dir, directory of samplesID files used for the "
        "superpopulation breakdown plot; auto-detected near the install when "
        "omitted [OPTIONAL]\n"
        "\t--output, output ZIP path [default: "
        "<result-dir>/<jobid>_report.zip]\n"
    )
    sys.exit(1)


def generate_report() -> None:
    """Standalone entry point for the shareable-report generator.

    Delegates to PostProcess/generate_report.build_report so the same command
    (re)generates the shareable ZIP for both CLI and website runs (their
    Results/<jobid> folders are identical).
    """
    args = input_args[2:]
    if "--help" in args:
        print_help_generate_report()

    def _opt(*names):
        for name in names:
            if name in args:
                idx = args.index(name)
                if idx + 1 < len(args):
                    return args[idx + 1]
                sys.stderr.write(f"ERROR! Missing value for {name}\n")
                print_help_generate_report()
        return None

    result_dir = _opt("--result-dir")
    integrated_results = _opt("--integrated-results")
    samplesid_dir = _opt("--samplesID-dir")
    output = _opt("--output", "--out")

    if not result_dir and not integrated_results:
        sys.stderr.write(
            "ERROR! generate-report requires --result-dir or "
            "--integrated-results\n\n"
        )
        print_help_generate_report()

    # import lazily so the rest of the CLI never pays the matplotlib import cost
    sys.path.insert(0, os.path.join(script_path, "PostProcess"))
    import generate_report as _report_module

    out = _report_module.build_report(
        result_dir=result_dir,
        integrated_tsv=integrated_results,
        out_zip=output,
        samplesid_dir=samplesid_dir,
    )
    sys.stdout.write(f"Report written: {out}\n")


def print_help_complete_test():
    """
    Prints the help information for executing comprehensive testing of the
    complete-search functionality provided by CRISPRme.

    Raises:
        SystemExit: If the help information is displayed to guide users on
            executing comprehensive testing.
    """

    # write intro message to stdout
    sys.stderr.write(
        "Execute comprehensive testing for complete-search functionality "
        "provided by CRISPRme\n\n"
    )
    # list functionality options
    sys.stderr.write(
        "Options:\n"
        "\t--chrom, test the complete-search functionality on the specified "
        "chromosome (e.g., chr22). By default, the test is conducted on all "
        "chromosomes\n"
        "\t--vcf_dataset, VCFs dataset to be used during CRISPRme testing. "
        "Available options include 1000 Genomes (1000G) and Human Genome "
        "Diversity Project (HGDP). To use the combined dataset type '1000G+HGDP' "
        "The default dataset is 1000 Genomes.\n"
        "\t--thread, number of threads.\n"
        "\t--debug, debug mode.\n"
    )
    sys.exit(1)


def complete_test_crisprme():
    """
    Executes comprehensive testing for the complete-search functionality provided
    by CRISPRme based on specified arguments.

    Raises:
        OSError: If the CRISPRme test fails, indicating an issue with the testing
            process.
    """

    warn_low_memory()  # non-fatal low-memory warning (Docker Desktop, etc.)
    if "--help" in input_args or len(input_args) < 3:
        print_help_complete_test()
        sys.exit(1)
    chrom = "all"
    if "--chrom" in input_args:  # individual chrom to test
        try:
            chrom = input_args[input_args.index("--chrom") + 1]
            if chrom.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --chrom\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --chrom\n")
            sys.exit(1)
    vcf_dataset = "1000G"
    if "--vcf_dataset" in input_args:  # specified variant dataset
        try:
            vcf_dataset = input_args[input_args.index("--vcf_dataset") + 1]
            if vcf_dataset.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --vcf_dataset\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --vcf_dataset\n")
            sys.exit(1)
    threads = 4
    if "--thread" in input_args:  # number of threads to use during test
        try:
            threads = input_args[input_args.index("--thread") + 1]
            if threads.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --thread\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some value for flag --thread\n")
            sys.exit(1)
    debug = "--debug" in input_args  # run local or via conda/Docker
    # begin crisprme test
    script_test = os.path.join(script_path, "complete_test.py")
    code = subprocess.call(
        f"python {script_test} {chrom} {vcf_dataset} {threads} {debug}", shell=True
    )
    if code != 0:
        raise OSError(
            "\nCRISPRme complete test encountered an Error! See Results/crisprme-test-out/log_error.txt for details\n"
        )
    
def print_help_validate_test() -> None:
    """
    Prints the help information for validating off-target sites produced by the
    complete-test functionality in CRISPRme.

    Explains the validation requirements, describes how brute-force search and
    alignment results are compared to CRISPRme predictions, and lists the
    available command-line options.

    Raises:
        SystemExit: If the help information is displayed to guide users on
            executing the validation procedure.
    """
    # write intro message to stderr
    sys.stderr.write(
        "Validate off-target sites generated by the complete-test functionality. "
        "This functionality compares the off-targets identified by CRISPRme "
        "against results obtained via brute-force search and alignment.\n" 
        "Requirements:\n"
        "\t- The `complete-test` functionality must have been executed beforehand\n"
        "\t- Validation is supported only when complete-test was run using "
        "1000 Genomes variant data\n\n"
    )
    # list functionality options
    sys.stderr.write(
        "Options:\n"
        "\t--chrom, Validate off-target sites on a specific chromosome "
        "(e.g., chr22). If not provided, validation is performed across "
        "all available chromosomes\n"
        "\t--debug, debug mode\n"
    )
    sys.exit(1)


def validate_test():
    """
    Executes validation of off-target sites produced by the complete-test
    functionality in CRISPRme.

    Runs an external validation script to compare CRISPRme-predicted off-targets
    against reference results on a specified chromosome or across all chromosomes.

    Raises:
        OSError: If the validation script returns a non-zero exit code,
            indicating that an error occurred during off-target validation.
    """
    if "--help" in input_args or len(input_args) < 2:
        print_help_validate_test()
        sys.exit(1)
    chrom = "all"
    if "--chrom" in input_args:  # individual chrom to test
        try:
            chrom = input_args[input_args.index("--chrom") + 1]
            if chrom.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --chrom\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --chrom\n")
            sys.exit(1)
    # begin crisprme test
    script_validation = os.path.join(script_path, "validate.py")
    code = subprocess.call(
        f"{sys.executable} {script_validation} {chrom}", shell=True
    )
    if code != 0:
        raise OSError("CRISPRme off-target sites validation encountered an Error!")
    

def print_help_setup_database_test() -> None:
    # write intro message to stderr
    sys.stderr.write(
        "This command initializes the CRISPRme legacy database by downloading "
        "all reference genomes, variant datasets, PAM definition files, and "
        "associated resources originally distributed through the CRISPRme web "
        "server.\n\n"
        "The downloaded resources can then be reused across analyses without "
        "requiring additional downloads\n\n"
    )
    # list functionality options
    sys.stderr.write(
        "Options:\n"
        "\t--path, Path to the directory where the legacy database and "
        "associated resources will be installed "
        "[default: current working directory]\n"
        "\t--chrom, download data for the specified chromsome only "
        "(e.g., chr22) [default: all]\n"
        "\t--force, force data download even if the database is already present "
        "[default: do not force]\n"
        "\t--debug, debug mode\n"
    )
    sys.exit(1)
    

def setup_database():
    """
    Setup CRISPRme legacy dataset downloading all genome, variant datasets, PAM
    files originally available in the CRISPRme website.

    Runs an external download script to retrieve the data and build the legacy 
    database.

    Raises:
        OSError: If the download script returns a non-zero exit code,
            indicating that an error occurred while setting up the legacy
            database.
    """
    if "--help" in input_args or len(input_args) < 2:
        print_help_setup_database_test()
        sys.exit(1)
    working_dir = os.path.abspath(os.getcwd())
    if "--path" in input_args:
        try:
            working_dir = os.path.abspath(input_args[input_args.index("--path") + 1])
            if working_dir.startswith("--"):
                raise ValueError
        except (IndexError, ValueError):
            sys.stderr.write(
                "Please provide a value for --path (e.g., /path/to/my/folder)"
            )
            sys.exit(1)
    chrom = "all"
    if "--chrom" in input_args:  # individual chrom to test
        try:
            chrom = input_args[input_args.index("--chrom") + 1]
            if chrom.startswith("--"):
                raise ValueError
        except (IndexError, ValueError):
            sys.stderr.write("Please provide a value for --chrom (e.g. chr22 or all)\n")
            sys.exit(1)
    force = "--force" in input_args
    # begin crisprme test
    script_setup = os.path.join(script_path, "setup_legacy_database.py")
    try:
        subprocess.run(
            [sys.executable, script_setup, chrom, working_dir, str(force)], check=True
        )
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
                f"Legacy database setup exited with error code {e.returncode}"
            )


# HELP FUNCTION
def crisprme_help() -> None:
    """
    Prints the general help information for CRISPRme, describing each available 
    functionality.

    Outputs usage instructions, requirements for input files, and a summary of 
    all main commands to stderr, then exits the program.
    """
    # print crisprme help; describe each functionality
    sys.stderr.write(
        "Help:\n\n"
        "- ALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND SEPARATED BY CHROMOSOME\n"
        "- ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED (WITH BGZIP) AND SEPARATED BY CHROMOSOME\n\n"
        "Functionalities:\n\n"
        "crisprme.py complete-search\n"
        "\tPerforms genome-wide off-targets search (reference and variant, if "
        "specified), including CFD and CRISTA analysis, and target selection\n\n"
        "crisprme.py complete-test\n"
        "\tTest the complete CRISPRme pipeline on single chromosomes or complete "
        "genomes\n\n"
        "crisprme.py build-index-only\n"
        "\tPre-builds the reusable CRISPRitz reference index for bulge-enabled "
        "searches (genome + PAM + bulges) without running a search\n\n"
        "crisprme.py download\n"
        "\tFast-downloads reference data (genome, annotations, PAMs, samples, "
        "VCFs, precomputed indexes) from a HuggingFace dataset repository\n\n"
        "crisprme.py publish-index\n"
        "\tUploads a locally built reference index to a HuggingFace dataset "
        "repository for later reuse via 'download --what index'\n\n"
        "crisprme.py assembly-search\n"
        "\tSearches a personal diploid genome assembly (two haplotypes, no VCF) "
        "and reconciles off-target predictions across both, mapped to hg38\n\n"
        "crisprme.py validate-test\n"
        "\tValidate targets obtained from complete-test by comparing them against "
        "brute-force search and alignment results\n\n"
        "crisprme.py targets-integration\n"
        "\tIntegrates in-silico targets with empirical data to generate a usable "
        "panel\n\n"
        "crisprme.py gnomAD-converter\n"
        "\tConverts gnomAD VCF files into CRISPRme compatible VCFs (supports "
        "gnomAD >= v3.1)\n\n"
        "crisprme.py generate-personal-card\n"
        "\tGenerates a personal card for specific samples by extracting all "
        "private targets\n\n"
        "crisprme.py generate-report\n"
        "\tBuilds a self-contained, shareable report (<jobid>_report.zip: "
        "offline report.html + integrated_results.tsv.gz) for a completed "
        "run\n\n"
        "crisprme.py setup\n"
        "\tInitializes the legacy database by downloading all reference "
        "genomes, variant datasets, PAM definition files, and associated "
        "resources\n\n"
        "crisprme.py web-interface\n"
        "\tActivates CRISPRme's web interface for local browser use\n\n"
        "crisprme.py --version\n"
        "\tPrints CRISPRme version to stdout and exit\n\n"
        "For additional information on each CRISPRme functionality type <function> "
        "--help (e.g. 'crisprme.py complete-search --help')\n"
    )
    sys.exit(1)  # stop execution


if len(sys.argv) < 2:
    check_crisprme_dirtree()  # check crisprme directory tree structure
    crisprme_help()  # no arg? print help
elif sys.argv[1] == "complete-search":  # run complete search
    complete_search()
elif sys.argv[1] == "complete-test":  # run complete test
    complete_test_crisprme()
elif sys.argv[1] == "assembly-search":  # run diploid assembly search
    assembly_search()
elif sys.argv[1] == "build-index-only":  # pre-build reusable reference index
    build_index_only()
elif sys.argv[1] == "download":  # fast HuggingFace download of reference data
    download_data()
elif sys.argv[1] == "publish-index":  # upload a prebuilt index to HuggingFace
    publish_index_cmd()
elif sys.argv[1] == "validate-test":  # run validate complete-test
    validate_test()
elif sys.argv[1] == "targets-integration":  # run targets integration
    target_integration()
elif sys.argv[1] == "gnomAD-converter":  # run gnomad converter
    gnomAD_converter()
elif sys.argv[1] == "generate-personal-card":  # run create personal card
    personal_card()
elif sys.argv[1] == "generate-report":  # build shareable self-contained report
    generate_report()
elif sys.argv[1] == "setup":  # run legacy database setup
    setup_database()
elif sys.argv[1] == "web-interface":  # run web interface
    web_interface()
elif sys.argv[1] == "--version":  # print version
    crisprme_version()
elif sys.argv[1] in ("--help", "-h", "help"):  # explicit help request (no error)
    crisprme_help()
else:
    sys.stderr.write(f"ERROR! {sys.argv[1]} is not an allowed command!\n\n")
    crisprme_help()  # print help if invalid command is given
