""" """

from pysam import TabixFile
from pysam.utils import SamtoolsError
from typing import Dict, List, Tuple, Optional
from time import time

import fcntl
import pysam
import sys
import os

TBI = "tbi"

# A valid bgzipped file ends with this 28-byte empty-BGZF-block "EOF marker";
# its absence is htslib's definition of a truncated file.
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _bgzf_eof_ok(path: str) -> bool:
    """True if ``path`` ends with the BGZF EOF marker (i.e. not truncated)."""
    try:
        if os.path.getsize(path) < len(_BGZF_EOF):
            return False
        with open(path, "rb") as fh:
            fh.seek(-len(_BGZF_EOF), os.SEEK_END)
            return fh.read(len(_BGZF_EOF)) == _BGZF_EOF
    except OSError:
        return False


def _reindex(annotation_fname: str) -> None:
    """(Re)build the tabix index for a bgzipped BED, dropping any stale/corrupt one."""
    try:
        os.remove(annotation_fname + ".tbi")
    except OSError:
        pass
    pysam.tabix_index(annotation_fname, force=True, preset="bed")


def parse_commandline(args: List[str]) -> Tuple[str, str, str]:
    """Parses and validates command line arguments for the annotation script.

    Ensures the correct number of arguments are provided and that input files exist. 
    Returns the file paths for the offtargets, annotation, and output files.

    Args:
        args: List of command line arguments.

    Returns:
        Tuple containing the offtargets file name, annotation file name, and output 
            file name.

    Raises:
        ValueError: If the number of arguments is incorrect.
        FileNotFoundError: If any of the specified files do not exist.
    """
    if len(args) != 3:  # check input arguments consistency
        raise ValueError("Too many/few input arguments for annotation script")
    offtargets_fname = args[0]  # offtargets table report
    if not os.path.isfile(offtargets_fname):
        raise FileNotFoundError(f"Cannot find off-targets file {offtargets_fname}")
    annotation_fname = args[1]
    if os.path.basename(annotation_fname) != "vuoto.txt" and not os.path.isfile(annotation_fname):
        raise FileNotFoundError(
            f"Cannot find annotation files {annotation_fname}"
        )
    offtargets_out_fname = args[2]  # annotated offtargets table report
    assert offtargets_out_fname
    return offtargets_fname, annotation_fname, offtargets_out_fname


def load_annotation_bed(annotation_fname: str) -> TabixFile:
    """Loads a BED annotation file and ensures a Tabix index is present.

    Checks for the existence of a Tabix index for the given annotation BED file 
    and creates one if necessary. Returns a TabixFile object for querying the 
    annotation data.

    Args:
        annotation_fname: Path to the annotation BED file.

    Returns:
        TabixFile object for the annotation BED file.

    Raises:
        SamtoolsError: If the annotation BED file cannot be loaded or indexed.
    """
    # Guard against a truncated/corrupt annotation download BEFORE htslib does,
    # so the failure is an actionable message instead of the cryptic
    # "EOF marker is absent ... could not open index" traceback. A bgzipped BED
    # (.gz) must end with the BGZF EOF marker; its absence means it is truncated.
    if annotation_fname.endswith(".gz") and not _bgzf_eof_ok(annotation_fname):
        raise SamtoolsError(
            f"Annotation file '{annotation_fname}' is truncated or corrupt "
            "(incomplete download). Re-download it and retry:\n"
            "  crisprme.py download --what annotations --path <DATA_DIR>\n"
            "The download now verifies integrity and re-fetches only bad files. "
            "(Or delete this file and its .tbi, then re-run the download.)"
        )
    # check that tabix index is available and up-to-date; use a lock file so that
    # parallel annotation processes don't write the .tbi simultaneously. A
    # missing / empty / stale .tbi is rebuilt -- a shipped index whose mtime
    # predates its .bed.gz, or a half-written one, would otherwise fail the load.
    tbi_path = annotation_fname + ".tbi"
    lock_path = annotation_fname + ".tbi.lock"
    with open(lock_path, "w") as lock_fh:
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
        if (not os.path.isfile(tbi_path)
                or os.path.getsize(tbi_path) == 0
                or os.path.getmtime(tbi_path) < os.path.getmtime(annotation_fname)):
            _reindex(annotation_fname)
    try:  # return tabix indexes for each annotation bed
        return pysam.TabixFile(annotation_fname)
    except Exception:
        # self-heal once: a corrupt .tbi that slipped past the checks above
        # (e.g. same mtime as the bed) -- rebuild it and retry before giving up.
        with open(lock_path, "w") as lock_fh:
            fcntl.flock(lock_fh, fcntl.LOCK_EX)
            _reindex(annotation_fname)
        try:
            return pysam.TabixFile(annotation_fname)
        except Exception as e:
            raise SamtoolsError(
                f"Could not load annotation file '{annotation_fname}' even after "
                "rebuilding its index -- the file may be corrupt. Re-download it: "
                "crisprme.py download --what annotations --path <DATA_DIR>"
            ) from e


def compute_target_coords(fields: List[str]) -> Tuple[str, int, int]:
    """Computes the genomic coordinates for a target from a list of fields.

    Returns the chromosome, start, and stop positions for the target based on 
    the provided fields.

    Args:
        fields: List of string fields containing target information.

    Returns:
        Tuple containing the chromosome (str), start (int), and stop (int) 
            positions.
    """
    start = int(fields[5])  # retrieve target start position
    stop = start + len(fields[1].replace("-", ""))  # compute target stop
    return fields[4], start, stop


def annotate_target(chrom: str, start: int, stop: int, annotation: TabixFile) -> str:
    """Annotates a genomic target using a Tabix-indexed annotation file.

    Retrieves and returns a comma-separated list of annotation features overlapping 
    the specified genomic region.

    Args:
        chrom: Chromosome name.
        start: Start position of the target.
        stop: Stop position of the target.
        annotation: TabixFile object for annotation data.

    Returns:
        Comma-separated string of annotation features for the target.
    """
    if chrom not in annotation.contigs:
        return "n"  # contig not in input annotation file 
    target_anns = {
        feature.split()[3]
        for feature in annotation.fetch(chrom, start, stop)
    }  # retrieve annotations for current target
    return ",".join(sorted(target_anns))  # report as comma-separated list


# Annotation-label buckets, and the column each bucket is written to.
#
# Deliberately identical to the classification resultIntegrator.py already
# performs on its own raw annotation field, so that a complete-search report and
# an assembly-search report describe the same annotation bundle with the same
# column names and the same values. Two properties of that classification are
# load-bearing and must not be "tidied":
#   * the suffixed kinds are tested BEFORE the ENCODE catch-all -- a COSMIC or
#     IntOGen label would otherwise be bucketed as ENCODE and silently mislabelled;
#   * the suffix is matched as a SUBSTRING and stripped with str.replace, which is
#     what turns "gene_gencode" into the "gene" that complete-search reports.
# ENCODE is the catch-all, so an unrecognised label from a user-supplied bundle is
# never dropped -- it is reported there rather than lost.
ANNOTATION_SUFFIX_COLS = (
    ("_personal", "Annotation_personal"),
    ("_gencode", "Annotation_GENCODE"),
    ("_DHS", "Annotation_DHS"),
    ("_COSMIC", "Annotation_COSMIC"),
    ("_INTOGEN", "Annotation_INTOGEN"),
)
ANNOTATION_CATCHALL_COL = "Annotation_ENCODE"
ANNOTATION_SPLIT_COLS = tuple(col for _, col in ANNOTATION_SUFFIX_COLS) + (
    ANNOTATION_CATCHALL_COL,
)
# "nothing overlaps here" markers: annotate_target()'s own "n", plus the usual
# blank spellings. A kind with no labels is reported blank, matching
# complete-search (whose per-kind columns are "NA"/empty when the kind is absent).
_NO_ANNOTATION_TOKENS = frozenset({"", "n", "na", "nan", "none", "."})


def split_annotation_labels(annotation: Optional[str]) -> Dict[str, str]:
    """Buckets one comma-separated annotation feature list into per-kind columns.

    Args:
        annotation: A feature list as `annotate_target` returns it (e.g.
            "dELS,gene_gencode,Neural_DHS"), or a blank/"n" marker.

    Returns:
        `{column_name: comma-joined sorted value}` for every column in
        `ANNOTATION_SPLIT_COLS`; a kind with no labels maps to "".
    """
    buckets: Dict[str, set] = {col: set() for col in ANNOTATION_SPLIT_COLS}
    text = "" if annotation is None else str(annotation).strip()
    if text.lower() in _NO_ANNOTATION_TOKENS:
        return {col: "" for col in ANNOTATION_SPLIT_COLS}
    for elem in text.split(","):
        elem = elem.strip()
        if not elem:
            continue
        for suffix, col in ANNOTATION_SUFFIX_COLS:
            if suffix in elem:
                buckets[col].add(elem.replace(suffix, ""))
                break
        else:
            buckets[ANNOTATION_CATCHALL_COL].add(elem)
    # sorted() so column order is deterministic across processes (bare set
    # iteration is PYTHONHASHSEED-dependent) -- same reason resultIntegrator sorts
    return {col: ",".join(sorted(vals)) for col, vals in buckets.items()}


def annotate_offtargets(
    offtargets_fname: str, annotation: Optional[TabixFile], offtargets_out_fname: str,
) -> None:
    """Annotates a genomic target using a Tabix-indexed annotation file.

    Retrieves and returns a comma-separated list of annotation features overlapping 
    the specified genomic region.

    Args:
        chrom: Chromosome name.
        start: Start position of the target.
        stop: Stop position of the target.
        annotation: TabixFile object for annotation data.

    Returns:
        Comma-separated string of annotation features for the target.
    """
    try:
        with open(offtargets_fname, mode="r") as infile, open(
            offtargets_out_fname, mode="w"
        ) as outfile:
            header = infile.readline().strip()  # preserve header in annotated report
            outfile.write(f"{header}\n")
            for offtarget in infile:  # read offtargets
                # explicit tab-split (not the default whitespace-collapsing
                # .split()): a genuinely empty column (e.g. a blank AF value)
                # sits between two tabs and must survive as its own "" field,
                # not be silently dropped, which would shift every later
                # column left by one and break the fixed-index readers
                # downstream (add_risk_score.py, resultIntegrator.py).
                fields = offtarget.rstrip("\n").split("\t")
                # annotate current off-target using input annotation datasets
                fields[14] = (
                    "n" if annotation is None else annotate_target(
                        *compute_target_coords(fields), annotation
                    )
                )
                offtarget_annotated = "\t".join(fields)
                outfile.write(f"{offtarget_annotated}\n")  # write annotated offtarget
    except (IOError, Exception) as e:
        raise OSError(f"Annotation failed on off-targets in {offtargets_fname}") from e


def main() -> None:
    """Annotates a list of off-targets using an optional annotation dataset.

    Reads off-targets from a file, annotates each entry with features from the 
    provided annotation dataset, and writes the results to an output file.

    Args:
        offtargets_fname: Path to the input file containing off-targets.
        annotation: TabixFile object for annotation data, or None if no annotation.
        offtargets_out_fname: Path to the output file for annotated off-targets.

    Raises:
        OSError: If annotation fails due to file I/O or processing errors.
    """
    # parse input command line arguments
    offtargets_fname, annotation, offtargets_out_fname = parse_commandline(sys.argv[1:])
    start = time()  # track annotation time
    empty = os.path.basename(annotation) == "vuoto.txt"
    # load annotation bed files
    annotation = None if empty else load_annotation_bed(annotation)
    # annotate offtargets using input bed files
    annotate_offtargets(offtargets_fname, annotation, offtargets_out_fname)
    sys.stdout.write(f"Annotation completed in {time() - start:.2f}s\n")


if __name__ == "__main__":
    main()