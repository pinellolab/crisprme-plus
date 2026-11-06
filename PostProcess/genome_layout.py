"""Genome folder layout: one FASTA file per sequence.

CRISPRme identifies a genome's chromosomes by FASTA *file name* (``chr1.fa``
-> ``chr1``) -- `submit_job_automated_new_multiple_vcfs.sh` builds its
chromosome list from ``$ref_folder/*.fa``, and the post-analysis, result
concatenation and variant enrichment are all keyed by that name. CRISPRitz
itself happily indexes and searches a single multi-record FASTA, so a genome
delivered as one file (e.g. UCSC's ``susScr11.fa.gz`` for the pig, which has no
per-chromosome ``chromFa`` archive) passes the index and search steps and then
silently yields an EMPTY result: the whole genome is treated as one
"chromosome" named after the file, and no hit matches it.

This module detects that layout and fixes it by splitting each multi-record
file into one file per sequence (named after the record id, the way UCSC's own
per-chromosome archives are). A file of one sequence not named after it
(``hg38_chr22.fa`` holding ``>chr22``, ``ecoli.fa`` holding ``>NC_000913.3``)
fails the same way, and gets the same fix: it is named after its sequence.

A header line must also be ``>id`` and nothing else, on lines that end in LF alone.
CRISPRitz names a chromosome after the whole header line, carriage return and all, and
the pipeline splits on blanks and takes each name into shell commands: an NCBI, GENCODE
or Ensembl header (``>NC_000913.3 Escherichia coli ...``, ``>chr1 1``, ``>1 dna:chromosome
...``), or CRLF line ends, indexed as done, and every search of it then lost results,
often all of them with no error. The same fix writes such a file as ``>id`` with no
carriage return. An id the pipeline cannot take (``gi|49175990|ref|...``, a PanSN
``HG01255#1#CM086702.1``, one that starts with ``-`` or ``.``) is refused: it is for the
user to rename. So is a folder where an id that holds ``*``, ``?`` or ``[`` matches, read
as a shell pattern, the name of another of its files (``HLA-A*02:10`` and
``HLA-A*24:02:10`` of the GATK GRCh38 FASTA): a search of it failed, or lost the
off-targets of one of them with no error.

This module is used by:

* ``utils.download_reference_genome`` -- split right after a UCSC/URL download;
* ``pages/settings_page.py`` -- the same, for a genome uploaded in the web Settings
  page;
* ``validate_inputs.check_genome_fasta`` -- fail fast (with the fix) before a
  ``complete-search`` starts;
* ``submit_job_automated_new_multiple_vcfs.sh`` -- the same guard for the
  web-launched path, which does not run the validator.

Command line::

    python genome_layout.py <genome_dir>            # fast check; exit 1 if a file has >1 sequence,
                                                    # a header that is not >id, CR line ends, or
                                                    # an id a search cannot take
    python genome_layout.py --deep <genome_dir>     # same, reading every file completely
    python genome_layout.py --split <genome_dir>    # split such files in place (implies --deep),
                                                    # name each file of one sequence after it, and
                                                    # write every header as >id, with no CR
"""

import fnmatch
import mmap
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

FASTA_EXTENSIONS = (".fa", ".fasta")
# In a sequence id: each makes a shell command that holds the chromosome name fail (or,
# '#', makes the tables that hold it lose its rows). Measured through the real CLI, with
# ":", "." and "_" among those that work inside an id.
UNUSABLE_IN_ID = "|;'\"$()&<>`#"
# At the start of a sequence id: the post-analysis greps for the name, unquoted, so grep
# took '-x' for an option, and the search script lists the chromosomes as the folder's
# *.fa, which leaves out '.x.fa'. A search of either ended with exit status 0 and none of
# its off-targets (measured through the real CLI).
UNUSABLE_ID_START = "-."
# In a sequence id, what makes the shell read it as a pattern: the search script and its
# post-analysis take each name into commands unquoted, where it stands for every file
# name it matches. Measured through the real CLI, an id that matches another name of its
# folder made a search fail, or drop its off-targets with no error, and one that matches
# a file of the PostProcess folder, where the post-analysis runs ('a*' matches
# adjust_cols.py), made it fail; find_pattern_clashes looks for the first only. Ids that
# matched neither worked (HLA-A*01:01, KIR2DL1*001, c[12] beside cA).
PATTERN_CHARS = "*?["
HEADER_LIMIT = 1 << 16  # (the most of a line read to judge it; no id comes near this)


def fasta_record_ids(path: str, limit: Optional[int] = None) -> List[str]:
    """Returns the record ids (first whitespace-delimited token of each '>'
    header line) of a FASTA file, in order.

    Scans the file as a memory map (a multi-GB genome is a matter of seconds)
    instead of iterating lines. ``limit`` stops after that many headers.
    """
    ids: List[str] = []
    if os.path.getsize(path) == 0:
        return ids
    with open(path, "rb") as handle:
        with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            pos = 0 if mm[:1] == b">" else mm.find(b"\n>") + 1
            if pos <= 0 and mm[:1] != b">":
                return ids  # no header at all
            while pos >= 0:
                end = mm.find(b"\n", pos)
                header = mm[pos + 1 : end if end != -1 else len(mm)]
                token = header.split(None, 1)[0] if header.strip() else b""
                ids.append(token.decode("utf-8", errors="replace"))
                if limit is not None and len(ids) >= limit:
                    break
                nxt = mm.find(b"\n>", end if end != -1 else len(mm))
                pos = nxt + 1 if nxt != -1 else -1
    return ids


def _stem(fname: str) -> str:
    return fname.rsplit(".", 1)[0]


def find_multi_record_fastas(genomedir: str, deep: bool = False) -> Dict[str, List[str]]:
    """Maps each FASTA file in ``genomedir`` that holds more than one sequence
    to its record ids.

    Default (fast) mode reads only each file's FIRST header: a file whose first
    record id equals its file name (``chr1.fa`` starting ``>chr1``) is taken to
    be a normal one-sequence-per-file chromosome and is not scanned further; any
    other file is scanned in full. That costs one line per file on a standard
    genome instead of reading every byte (~30 s for hg38 on network storage), and
    catches the real cases (a genome shipped as one multi-sequence file is
    named after the assembly, not its first sequence). ``deep=True`` scans every
    file completely, which also finds e.g. a ``chr1.fa`` with extra sequences
    appended.
    """
    found: Dict[str, List[str]] = {}
    for fname in sorted(os.listdir(genomedir)):
        if not fname.endswith(FASTA_EXTENSIONS):
            continue
        path = os.path.join(genomedir, fname)
        if not deep:
            first = fasta_record_ids(path, limit=1)
            if first and first[0] == _stem(fname):
                continue
        ids = fasta_record_ids(path)
        if len(ids) > 1:
            found[fname] = ids
    return found


def find_name_mismatched_fastas(genomedir: str) -> Dict[str, str]:
    """Maps single-sequence FASTA files whose file name differs from their
    sequence id to that id (``NC_000001.fa`` holding ``>chr1``). Cheap: one
    header line per file."""
    out: Dict[str, str] = {}
    for fname in sorted(os.listdir(genomedir)):
        if not fname.endswith(FASTA_EXTENSIONS):
            continue
        first = fasta_record_ids(os.path.join(genomedir, fname), limit=1)
        if first and first[0] != _stem(fname):
            out[fname] = first[0]
    return out


def _unclean_reason(path: str, deep: bool = False) -> Optional[str]:
    """Why a FASTA of one sequence is not as the pipeline needs it, or None: its first
    line is more than ``>id``, or its lines end in CR (the first two lines tell; ``deep``
    looks for a CR in the whole file).

    A CR that only ends the file counts too. The index build takes each CR on a line
    that ends in LF (the CRLF of the last line) for a base of the sequence, and the
    variant enricher keeps every CR, one after the last LF too, in the one line of
    sequence it writes, which the build of the variant index then reads. Measured
    through the real engine, a PAM in the last bases of the sequence then gave a target
    window that held the CR: the build wrote warnings on stderr at exit status 0, and a
    search of the index crashed (SIGSEGV) or reported other rows than the file with no
    CR."""
    with open(path, "rb") as handle:
        head = handle.readline(HEADER_LIMIT)
        first = handle.readline(HEADER_LIMIT)
    if not head.startswith(b">"):
        return None  # no header there: check_genome_fasta reports it
    line = head[1:].rstrip(b"\n")
    if b"\r" in line or b"\r" in first:
        return "its lines end in CR"
    if line.split() != [line]:
        return f"its header line reads '>{line.decode('utf-8', errors='replace')[:70]}'"
    if deep and os.path.getsize(path):
        with open(path, "rb") as handle:
            with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                if mm.find(b"\r") != -1:
                    return "its lines end in CR"
    return None


def find_unclean_fastas(genomedir: str, deep: bool = False) -> Dict[str, str]:
    """Maps each FASTA file in ``genomedir`` whose header is not exactly ``>id``, or
    whose lines end in CR, to what is wrong with it (``_unclean_reason``). Cheap unless
    ``deep``: two lines per file. A file of several sequences is only looked at for
    its first header: a split writes each as it should be."""
    found: Dict[str, str] = {}
    for fname in sorted(os.listdir(genomedir)):
        if fname.endswith(FASTA_EXTENSIONS):
            reason = _unclean_reason(os.path.join(genomedir, fname), deep)
            if reason:
                found[fname] = reason
    return found


def find_unterminated_fastas(genomedir: str) -> List[str]:
    """The FASTA files of ``genomedir`` whose last line has no line end (LF). The index
    build drops such a line with no error: a search lost the off-targets of its bases,
    all of a sequence written as one line (LP unterminated-last-line-dropped). Cheap: one
    byte a file."""
    found: List[str] = []
    for fname in sorted(os.listdir(genomedir)):
        path = os.path.join(genomedir, fname)
        if not fname.endswith(FASTA_EXTENSIONS) or not os.path.isfile(path):
            continue
        with open(path, "rb") as handle:
            if handle.read(1) != b">":  # (no FASTA: check_genome_fasta says so)
                continue
            handle.seek(-1, os.SEEK_END)
            if handle.read(1) != b"\n":
                found.append(fname)
    return found


def find_pattern_clashes(genomedir: str) -> Dict[str, List[str]]:
    """Maps each ``.fa`` file of ``genomedir`` whose name, read as a shell pattern (it
    holds ``*``, ``?`` or ``[``), matches the names of other ``.fa`` files there to those
    names: ``HLA-A*02:10.fa`` to ``['HLA-A*24:02:10.fa']``. The search script lists the
    chromosomes as ``basename $file``, and its post-analysis names each chromosome's
    files so: with ``c[12].fa`` beside ``c1.fa`` it took c1 twice and c[12] never, so
    a search with ``--thread 1`` ended with no c[12] off-target and no error, and with
    more threads the chromosomes' steps removed or read each other's files and it
    failed. Hidden files are left out, as the shell leaves them out of a match."""
    names = sorted(f for f in os.listdir(genomedir)
                   if f.split(".")[-1] == "fa" and not f.startswith("."))
    found: Dict[str, List[str]] = {}
    for fname in names:
        if not any(c in fname for c in PATTERN_CHARS):
            continue
        # (the shell's [^...] is fnmatch's [!...]; a class such as [[:digit:]], fnmatch
        # does not read as the shell does)
        pattern = fname.replace("[^", "[!")
        others = [o for o in names if o != fname and fnmatch.fnmatchcase(o, pattern)]
        if others:
            found[fname] = others
    return found


def _unusable_char(rid: str) -> Optional[str]:
    return next((c for c in UNUSABLE_IN_ID if c in rid), None)


def _unusable_start(rid: str) -> Optional[str]:
    return rid[0] if rid and rid[0] in UNUSABLE_ID_START else None


def _validate_ids(fname: str, ids: List[str]) -> None:
    seen = set()
    for rid in ids:
        if not rid or any(c in rid for c in ("/", "\\", "\0")) or rid.startswith("."):
            raise ValueError(
                f"{fname}: sequence id {rid!r} cannot be used as a file name; "
                "rename the sequence in the FASTA header and retry"
            )
        if _unusable_char(rid):
            raise ValueError(
                f"{fname}: sequence id {rid!r} holds {_unusable_char(rid)!r}, which the "
                "pipeline cannot take; rename the sequence in the FASTA header and retry"
            )
        if _unusable_start(rid):
            raise ValueError(
                f"{fname}: sequence id {rid!r} starts with {_unusable_start(rid)!r}, which "
                "the pipeline cannot take; rename the sequence in the FASTA header and retry"
            )
        if rid in seen:
            raise ValueError(f"{fname}: duplicate sequence id {rid!r}")
        seen.add(rid)


def split_multi_fasta(path: str, dest_dir: Optional[str] = None) -> List[str]:
    """Splits a multi-record FASTA into ``<record id>.fa`` files in ``dest_dir``
    (default: the file's own directory) and removes the original. A FASTA of one
    record goes the same way: it is renamed after its record id.

    Each header is written as ``>id``, and every line with no CR (a CR that ends a line
    alone ends it with LF), and each piece ends with LF; the rest is copied verbatim, so
    the sequence of the pieces is the original's. (The index build drops a last line
    with no LF: lone CRs deleted, not turned into LF, left a line of tens of kb with none.) Refuses
    (before anything is put in place) if two records share an id, an id is not a
    usable file name or holds what the pipeline cannot take, or the file's lines end in
    CR alone, and never overwrites an existing file other than the original itself (the
    same file under a name that differs in case only, on a case-insensitive file
    system, is the original).
    """
    dest_dir = dest_dir or os.path.dirname(os.path.abspath(path))
    fname = os.path.basename(path)
    ids = fasta_record_ids(path)
    _validate_ids(fname, ids)
    original = os.path.abspath(path)
    targets = {rid: os.path.join(dest_dir, f"{rid}.fa") for rid in ids}
    for rid, target in targets.items():
        if (os.path.exists(target) and os.path.abspath(target) != original
                and not os.path.samefile(target, original)):
            raise FileExistsError(
                f"{fname}: cannot split -- {target} already exists" if len(ids) > 1
                else f"{fname}: cannot name it after its sequence -- {target} already exists"
            )
    written: List[str] = []
    out = None
    tmp_names: List[Tuple[str, str]] = []
    try:
        with open(path, "rb") as fin:
            ended = True  # (whether what the piece holds so far ends in LF)
            for line in fin:
                if line.startswith(b">"):
                    if b"\r" in line.rstrip(b"\n")[:-1]:
                        # (a line holds the file: its header's id alone would be kept)
                        raise ValueError(
                            f"{fname}: its lines end in CR alone, so it reads as one line; "
                            "convert them to LF (tr '\\r' '\\n') and retry"
                        )
                    if out is not None:
                        if not ended:
                            out.write(b"\n")
                        out.close()
                    rid_bytes = line[1:].split(None, 1)[0]
                    rid = rid_bytes.decode("utf-8", errors="replace")
                    tmp = targets[rid] + ".part"
                    out = open(tmp, "wb")
                    tmp_names.append((tmp, targets[rid]))
                    out.write(b">" + rid_bytes + b"\n")
                    ended = True
                elif out is not None:
                    # a CR between bases ends a line (a file of lone-CR line ends); any
                    # other is dropped, as with CRLF line ends
                    line = re.sub(rb"(?<=[^\r\n])\r+(?=[^\r\n])", b"\n",
                                  line.replace(b"\r\n", b"\n")).replace(b"\r", b"")
                    if line:
                        out.write(line)
                        ended = line.endswith(b"\n")
        if out is not None:
            if not ended:
                out.write(b"\n")
            out.close()
            out = None
        os.remove(path)  # before the renames: a piece may reuse the original's name
        for tmp, final in tmp_names:
            os.replace(tmp, final)
            written.append(final)
    except BaseException:
        if out is not None:
            out.close()
        for tmp, _final in tmp_names:
            if os.path.exists(tmp):
                os.remove(tmp)
        raise
    return written


def split_multi_record_fastas(genomedir: str, deep: bool = True) -> Dict[str, List[str]]:
    """Splits every multi-record FASTA in ``genomedir`` in place, then names each
    FASTA of one sequence that is named otherwise after it (``walkone.fa`` holding
    ``>chr1`` becomes ``chr1.fa``): CRISPRme takes a file for the sequence it is
    named after, so every search of such a folder found its off-targets and then
    dropped them all. Each file it writes, and each other file whose header is not
    exactly ``>id``, whose lines end in CR, or whose last line has no LF, then reads
    ``>id`` with no CR, and ends with LF. A file that needs none of this is not touched. Returns ``{original file name: [new file
    paths]}`` (empty if nothing needed it). Deep by default: this is a one-off fix, not
    a per-search check."""
    result: Dict[str, List[str]] = {}
    for fname in find_multi_record_fastas(genomedir, deep=deep):
        result[fname] = split_multi_fasta(os.path.join(genomedir, fname), genomedir)
    for fname in find_name_mismatched_fastas(genomedir):
        result[fname] = split_multi_fasta(os.path.join(genomedir, fname), genomedir)
    for fname in find_unclean_fastas(genomedir, deep=deep):
        result[fname] = split_multi_fasta(os.path.join(genomedir, fname), genomedir)
    for fname in find_unterminated_fastas(genomedir):  # (each piece ends with LF)
        result[fname] = split_multi_fasta(os.path.join(genomedir, fname), genomedir)
    return result


def describe_problem(genomedir: str, fname: str, ids: List[str]) -> str:
    """Human-readable error for a multi-record FASTA, including the fix."""
    shown = ", ".join(ids[:3]) + (", ..." if len(ids) > 3 else "")
    return (
        f"{fname}: holds {len(ids)} sequences ({shown}). CRISPRme needs ONE FASTA "
        f"file per sequence, named after it (e.g. chr1.fa); with a multi-sequence "
        f"file the search runs but every result is silently dropped. Split it with:\n"
        f"    python PostProcess/genome_layout.py --split {genomedir}"
    )


def describe_header_problems(genomedir: str, deep: bool = False, reported=()) -> List[str]:
    """Human-readable errors, each with its fix, for the ``.fa`` files of ``genomedir``
    (those a search reads) that the search cannot use though they hold one sequence named
    after it: a header that is not exactly ``>id``, CR line ends, a last line with no LF,
    an id the pipeline cannot take, an id that matches, read as a shell pattern, another
    of their names. One message per kind, whatever the number of files; none for the
    files in ``reported``, as holding several sequences, whose fix is the same."""
    genome = os.path.basename(os.path.normpath(genomedir))
    errors: List[str] = []
    unclean = {f: why for f, why in find_unclean_fastas(genomedir, deep).items()
               if f.split(".")[-1] == "fa" and f not in reported}
    if unclean:
        fname, why = next(iter(unclean.items()))
        more = f" ({len(unclean) - 1} more file(s) here too)" if len(unclean) > 1 else ""
        errors.append(
            f"{fname}: {why}{more}. A header must be >id and nothing else, on lines that "
            f"end in LF: CRISPRitz names a chromosome after the whole line, so every "
            f"search of this folder loses results, all of them on many such headers, with "
            f"no error. Write them so with:\n"
            f"    python PostProcess/genome_layout.py --split {genomedir}\n"
            f"Then remove each index built from the folder as it was "
            f"(genome_library/*_{genome}, genome_library/*_{genome}+*, and the variant "
            f"genomes Genomes/{genome}+*), which keeps those names, and run build-index-only "
            f"again. (A genome uploaded from the web Settings page is written this way: "
            f"remove this one there, with the indexes built on it, and upload it again.)"
        )
    unended = [f for f in find_unterminated_fastas(genomedir)
               if f.split(".")[-1] == "fa" and f not in reported and f not in unclean]
    if unended:
        more = f" ({len(unended) - 1} more file(s) here too)" if len(unended) > 1 else ""
        errors.append(
            f"{unended[0]}: its last line has no line end (LF){more}. The index build drops "
            f"a last line with none, so every search of this folder loses the off-targets "
            f"of its bases (all of a sequence written on one line), with no error. End it "
            f"with one with:\n"
            f"    python PostProcess/genome_layout.py --split {genomedir}\n"
            f"Then remove each index built from the folder as it was "
            f"(genome_library/*_{genome}, genome_library/*_{genome}+*, and the variant "
            f"genomes Genomes/{genome}+*), and run build-index-only again."
        )
    held = started = None  # (the first file of each kind)
    for fname in sorted(os.listdir(genomedir)):
        if fname.split(".")[-1] != "fa" or fname in reported:
            continue
        ids = fasta_record_ids(os.path.join(genomedir, fname), limit=1)
        if ids and held is None and _unusable_char(ids[0]):
            held = (
                f"{fname}: sequence id {ids[0]!r} holds {_unusable_char(ids[0])!r}, which "
                f"the pipeline cannot take (nor any of {' '.join(UNUSABLE_IN_ID)}): a search "
                f"of it fails, or drops all its off-targets. Rename the sequence in the "
                f"FASTA header and the file after it, then remove each index built from the "
                f"folder and run build-index-only again."
            )
        if ids and started is None and _unusable_start(ids[0]):
            start = _unusable_start(ids[0])
            others = " or ".join(repr(c) for c in UNUSABLE_ID_START if c != start)
            started = (
                f"{fname}: sequence id {ids[0]!r} starts with {start!r}, which the pipeline "
                f"cannot take (nor one that starts with {others}): a search of it finds its "
                f"off-targets and then drops them all, with no error. Rename the sequence in "
                f"the FASTA header and the file after it, then remove each index built from "
                f"the folder and run build-index-only again."
            )
        if held and started:
            break
    errors += [problem for problem in (held, started) if problem]
    clashes = {f: o for f, o in find_pattern_clashes(genomedir).items() if f not in reported}
    if clashes:
        fname, others = next(iter(clashes.items()))
        shown = ", ".join(others[:3]) + (", ..." if len(others) > 3 else "")
        every = list(clashes)
        also = (f"; {len(every)} files here hold such an id ({', '.join(every[:10])}"
                f"{', ...' if len(every) > 10 else ''})" if len(every) > 1 else "")
        errors.append(
            f"{fname}: its sequence id {_stem(fname)!r}, read as a shell pattern, matches "
            f"{shown}{also}. A search takes each sequence's name into shell commands as it "
            f"is, where '*', '?' and '[' make it such a pattern: a search of this folder "
            f"fails, or leaves out the off-targets of such a sequence with no error. Rename "
            f"each such sequence in its FASTA header and its file after it, with none of "
            f"those characters (or remove its file), then remove each index built from the "
            f"folder and run build-index-only again. (A genome uploaded from the web "
            f"Settings page: remove it there, rename the sequences in your FASTA and upload "
            f"it again.)"
        )
    return errors


def describe_name_problems(genomedir: str, reported=()) -> List[str]:
    """The error, with its fix, for the ``.fa`` files of ``genomedir`` that hold one
    sequence and are named otherwise (walkone.fa holding >chr1), but those in
    ``reported``: a search takes each file for the sequence it is named after, so it finds
    that sequence's off-targets and drops them all. One message, whatever the number of
    files. The search, which runs this module's check first, ran its whole course for none,
    then said to run the check, which passed such a file (M1)."""
    misnamed = {f: rid for f, rid in find_name_mismatched_fastas(genomedir).items()
                if f.split(".")[-1] == "fa" and f not in reported}
    if not misnamed:
        return []
    genome = os.path.basename(os.path.normpath(genomedir))
    fname, rid = next(iter(misnamed.items()))
    more = f" ({len(misnamed) - 1} more file(s) here too)" if len(misnamed) > 1 else ""
    return [
        f"{fname}: holds the sequence {rid} and is named otherwise{more}. CRISPRme reads a "
        f"sequence from the file named after it ({rid}.fa), so a search of this folder finds "
        f"that sequence's off-targets and then drops them all, with no error. Name each file "
        f"after its sequence with:\n"
        f"    python PostProcess/genome_layout.py --split {genomedir}\n"
        f"An index built from the folder as it was is named after the sequence, and serves "
        f"once the files are renamed. (A genome uploaded from the web Settings page is named "
        f"this way: remove this one there, with the indexes built on it, and upload it again.)"
    ]


def _indexes_of(library: str, genome: str) -> List[str]:
    """The folders of the genome library named as an index of ``genome``:
    ``<PAM>_<N>_<genome>``, or ``<PAM>_<N>_<genome>+<dataset>`` of a variant one."""
    try:
        names = sorted(os.listdir(library))
    except OSError:
        return []
    named = re.compile(r"[A-Za-z]+_[0-9]+_" + re.escape(genome) + r"(\+.*)?\Z")
    return [name for name in names if named.match(name)]


def describe_index_problems(library: str, genome: str) -> List[str]:
    """Human-readable errors, each with its fix, for the indexes of ``genome`` in the
    genome library ``library`` that were built from a header of more than ``>id``, or
    CR line ends: each .bin is named after the whole line, and every search of the
    index loses results even once the FASTA is written as it should be."""
    errors: List[str] = []
    for name in _indexes_of(library, genome):
        index = os.path.join(library, name)
        try:
            bins = sorted(os.listdir(index))
        except OSError:
            continue
        unclean = [b for b in bins if b.endswith(".bin") and len(b.split()) != 1]
        if unclean:
            # (a variant index: the enriched genome it was built from holds the same headers)
            variant = name.split("_", 2)[2]
            also = f" and the variant genome Genomes/{variant}" if "+" in variant else ""
            errors.append(
                f"The index {index} holds {unclean[0]!r}: it was built from a FASTA header "
                f"of more than >id, or with CR line ends, and every search of it loses "
                f"results. Remove it (rm -rf '{index}', or remove the index in the web "
                f"Settings page){also} and run build-index-only again."
            )
    return errors


def main(argv: List[str]) -> int:
    split = "--split" in argv
    deep = "--deep" in argv or split
    args = [a for a in argv if a not in ("--split", "--deep")]
    if len(args) != 1 or not os.path.isdir(args[0]):
        sys.stderr.write("usage: genome_layout.py [--split|--deep] <genome_dir>\n")
        return 2
    genomedir = args[0]
    # (the library of the data folder that holds Genomes/<genome>, as for a web search)
    library = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(genomedir))),
                           "genome_library")
    genome = os.path.basename(os.path.normpath(genomedir))
    problems = find_multi_record_fastas(genomedir, deep=deep)
    if split:
        unclean = find_unclean_fastas(genomedir, deep=True)
        unended = find_unterminated_fastas(genomedir)
        try:
            done = split_multi_record_fastas(genomedir)
        except (ValueError, OSError) as err:
            sys.stderr.write(f"ERROR: {err}\n")
            return 1
        for fname, pieces in done.items():
            if fname in problems:
                print(f"Split {fname} into {len(pieces)} files in {genomedir}")
            elif os.path.basename(pieces[0]) != fname:
                # one sequence, named otherwise: an index built from it names it right
                print(f"Renamed {fname} to {os.path.basename(pieces[0])} in {genomedir}")
            if fname in unclean:
                print(f"Wrote {fname}'s header as >id and its lines with no CR ({unclean[fname]})")
            if fname in unended and fname not in unclean:
                print(f"Ended {fname}'s last line with a line end (LF)")
        if problems or unclean or unended:
            print(
                "NOTE: any CRISPRitz index built from the old layout must be rebuilt "
                "(delete it and re-run build-index-only)."
            )
            built = _indexes_of(library, genome)
            if built:
                print(f"NOTE: in {library}, built on this genome: {', '.join(built)}")
        return 0
    header_problems = describe_header_problems(genomedir, deep, problems)
    header_problems += describe_name_problems(genomedir, problems)
    header_problems += describe_index_problems(library, genome)
    if not problems and not header_problems:
        return 0
    for fname, ids in problems.items():
        sys.stderr.write("ERROR: " + describe_problem(genomedir, fname, ids) + "\n")
    for problem in header_problems:
        sys.stderr.write("ERROR: " + problem + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
