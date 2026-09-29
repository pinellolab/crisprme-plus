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
per-chromosome archives are). It is used by:

* ``utils.download_reference_genome`` -- split right after a UCSC/URL download;
* ``validate_inputs.check_genome_fasta`` -- fail fast (with the fix) before a
  ``complete-search`` starts;
* ``submit_job_automated_new_multiple_vcfs.sh`` -- the same guard for the
  web-launched path, which does not run the validator.

Command line::

    python genome_layout.py <genome_dir>            # fast check; exit 1 if a file has >1 sequence
    python genome_layout.py --deep <genome_dir>     # same, reading every file completely
    python genome_layout.py --split <genome_dir>    # split such files in place (implies --deep)
"""

import mmap
import os
import sys
from typing import Dict, List, Optional, Tuple

FASTA_EXTENSIONS = (".fa", ".fasta")


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


def _validate_ids(fname: str, ids: List[str]) -> None:
    seen = set()
    for rid in ids:
        if not rid or any(c in rid for c in ("/", "\\", "\0")) or rid.startswith("."):
            raise ValueError(
                f"{fname}: sequence id {rid!r} cannot be used as a file name; "
                "rename the sequence in the FASTA header and retry"
            )
        if rid in seen:
            raise ValueError(f"{fname}: duplicate sequence id {rid!r}")
        seen.add(rid)


def split_multi_fasta(path: str, dest_dir: Optional[str] = None) -> List[str]:
    """Splits a multi-record FASTA into ``<record id>.fa`` files in ``dest_dir``
    (default: the file's own directory) and removes the original.

    Lines are copied verbatim, so the concatenation of the pieces is byte-for-
    byte the original. Refuses (before writing anything) if two records share
    an id or an id is not a usable file name, and never overwrites an existing
    file other than the original itself.
    """
    dest_dir = dest_dir or os.path.dirname(os.path.abspath(path))
    fname = os.path.basename(path)
    ids = fasta_record_ids(path)
    _validate_ids(fname, ids)
    original = os.path.abspath(path)
    targets = {rid: os.path.join(dest_dir, f"{rid}.fa") for rid in ids}
    for rid, target in targets.items():
        if os.path.exists(target) and os.path.abspath(target) != original:
            raise FileExistsError(
                f"{fname}: cannot split -- {target} already exists"
            )
    written: List[str] = []
    out = None
    tmp_names: List[Tuple[str, str]] = []
    try:
        with open(path, "rb") as fin:
            for line in fin:
                if line.startswith(b">"):
                    if out is not None:
                        out.close()
                    rid = line[1:].split(None, 1)[0].decode("utf-8", errors="replace")
                    tmp = targets[rid] + ".part"
                    out = open(tmp, "wb")
                    tmp_names.append((tmp, targets[rid]))
                if out is not None:
                    out.write(line)
        if out is not None:
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
    """Splits every multi-record FASTA in ``genomedir`` in place. Returns
    ``{original file name: [new file paths]}`` (empty if nothing needed
    splitting). Deep by default: this is a one-off fix, not a per-search check."""
    result: Dict[str, List[str]] = {}
    for fname in find_multi_record_fastas(genomedir, deep=deep):
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


def main(argv: List[str]) -> int:
    split = "--split" in argv
    deep = "--deep" in argv or split
    args = [a for a in argv if a not in ("--split", "--deep")]
    if len(args) != 1 or not os.path.isdir(args[0]):
        sys.stderr.write("usage: genome_layout.py [--split|--deep] <genome_dir>\n")
        return 2
    genomedir = args[0]
    problems = find_multi_record_fastas(genomedir, deep=deep)
    if not problems:
        return 0
    if split:
        done = split_multi_record_fastas(genomedir)
        for fname, pieces in done.items():
            print(f"Split {fname} into {len(pieces)} files in {genomedir}")
        print(
            "NOTE: any CRISPRitz index built from the old single-file layout must be "
            "rebuilt (delete it and re-run build-index-only)."
        )
        return 0
    for fname, ids in problems.items():
        sys.stderr.write("ERROR: " + describe_problem(genomedir, fname, ids) + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
