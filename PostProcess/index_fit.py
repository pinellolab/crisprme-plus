"""Whether an index, whoever built it, can answer a search of a guide of a given length.

An index folder is named ``<PAM>_<N>_<genome>``, which says nothing of the guide it was
built for. The engine reads that length from each tree (``len_guide_used``, the second
int32 of a .bin, written by buildTST since 2019, whatever built the index), but only to
offset coordinates: it searches any guide it is given, and past what the tree stores it
reports wrong alignments or loses targets, at exit 0.

The rule, measured through the real CLI against an index dedicated to the guide (246 cells
of built length, N, guide length and bulges, dedicated, all-N and all-N under --prune-pam,
and 70 more of variant searches and a 5' PAM; submit_job_automated_new_multiple_vcfs.sh
_index_fits_guide says the same, and the two must not disagree):

    built + N >= guide + bDNA + (1 if bRNA else 0)

Inside it the results are byte-identical to a dedicated index; outside it, in every cell,
they differ. So a 21-nt guide is served at up to one bulge of each type by any index built
for 20 nt with N=3, and a 23-nt guide is served by none of the published ones.
"""

import os
import sys

# The first four bytes of an index file that CRISPRitz 2.9.2 and later write (index format
# 2: CRISPRitz's docs/INDEX_FORMAT.md): the magic, then the format version as an int32.
# A format-1 file (CRISPRitz up to 2.9.1, the only format this CRISPRme reads) opens with
# its leaf count, a positive int32, where the magic reads as a negative one.
CZI_MAGIC = b"CZI\xfe"


def bin_format(path):
    """The index format of the CRISPRitz .bin ``path``: 1 for a file with no header (up
    to CRISPRitz 2.9.1), else the version that follows the magic (2 from CRISPRitz 2.9.2).
    None where it cannot tell: a file too short to say, or that cannot be read."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError:
        return None
    if head[:4] != CZI_MAGIC:
        return 1 if len(head) == 8 else None
    if len(head) < 8:
        return None
    return int.from_bytes(head[4:], sys.byteorder, signed=True)


def newer_bin(index_dir):
    """(path, format) of the first .bin of ``index_dir``, or of a folder in it (an INDELS
    index holds one a chromosome), written in an index format newer than 1, which this
    CRISPRme does not read; None where there is none."""
    try:
        for top, dirs, files in os.walk(index_dir):
            dirs.sort()
            for name in sorted(files):
                if name.endswith(".bin"):
                    fmt = bin_format(os.path.join(top, name))
                    if fmt is not None and fmt != 1:
                        return os.path.join(top, name), fmt
    except OSError:
        pass
    return None


def newer_format_error(index_dir):
    """The ERROR for an index ``index_dir`` built by a newer CRISPRitz than this CRISPRme
    reads (newer_bin), or None. Its version would otherwise be read as a guide length (a
    format-2 index as one built for a 2-nt guide), and its leaf count as no tree."""
    found = newer_bin(index_dir)
    if found is None:
        return None
    path, fmt = found
    return (f"{index_dir} was built by a newer CRISPRitz ({os.path.relpath(path, index_dir)} "
            f"is in index format {fmt}, which CRISPRitz 2.9.2 and later write): this "
            f"CRISPRme reads index format 1 only (CRISPRitz up to 2.9.1), and needs a newer "
            f"release to search it. Use a CRISPRme that reads index format {fmt}, or install "
            f"CRISPRitz 2.9.1 and build the index again, with {index_dir} moved out of the "
            f"way first")


def built_guide_len(index_dir):
    """The guide length ``index_dir`` was built for, from the header of a .bin of it: the
    second int32, in the byte order of the machine that reads it, as the engine reads it.
    None where it cannot tell: no .bin, one too short, or one in an index format newer
    than 1 (newer_format_error), whose second int32 is its format version."""
    try:
        for name in sorted(os.listdir(index_dir)):
            if name.endswith(".bin"):
                with open(os.path.join(index_dir, name), "rb") as fh:
                    head = fh.read(8)
                if len(head) == 8 and head[:4] != CZI_MAGIC:
                    built = int.from_bytes(head[4:], sys.byteorder, signed=True)
                    return built if built > 0 else None
                return None
    except OSError:
        pass
    return None


def pam_file_spacer(pam_file):
    """The guide length a PAM file asks for (its sequence less its PAM), or None where the
    file cannot be read as one."""
    try:
        with open(pam_file) as fh:
            parts = fh.readline().split()
        return len(parts[0]) - abs(int(parts[1]))
    except (OSError, ValueError, IndexError):
        return None


def fits(built, n, guide, bdna, brna):
    """Whether an index built for a guide of ``built`` nt with N = ``n`` answers a guide of
    ``guide`` nt with these bulges, by the rule above. An index or a search whose length
    is not known (None or 0: a short-form PAM file has no spacer) fits, as before."""
    if not built or not guide:
        return True
    return built + n >= guide + bdna + (1 if brna else 0)


def bulge_cap(built, n, guide):
    """The most bulges of each type an index built for a guide of ``built`` nt with N =
    ``n`` serves for a guide of ``guide`` nt: N-1 (the name's budget), or fewer where the
    guide is longer than the index stores; 0 where it does not fit with one of each."""
    cap = n - 1
    if built and guide:
        cap = min(cap, built + n - guide - 1)
    return max(cap, 0)


# The bases each IUPAC code stands for (pam_filter.py's table).
_IUPAC_BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def pam_holds(index_pam, motif):
    """Whether PAM ``index_pam``, of a dedicated index, holds PAM ``motif``: the two are as
    long, and at each position every base ``motif`` stands for ``index_pam`` stands for too
    (NRR holds NGG, NAG, NRG and NGA; NRG holds NGG and NAG). Such an index stores every
    site of ``motif``, and under --prune-pam the engine writes out only those, at leaf time,
    as it does an all-N index's (the search script's _pam_holds, the same rule)."""
    if not motif or len(index_pam) != len(motif):
        return False
    try:
        return all(set(_IUPAC_BASES[m]) <= set(_IUPAC_BASES[b])
                   for m, b in zip(motif, index_pam))
    except KeyError:  # a letter that is no IUPAC code
        return False


def pam_serves(index_pam, motif, prune, three_prime):
    """How an index built for PAM ``index_pam`` serves a search for PAM ``motif``, as the
    search script's _resolve_index_dir ranks it (its _exact), with which this must not
    disagree: 2 as its own PAM; 1, under --prune-pam, as a dedicated index of a PAM that
    holds it (pam_holds); 0 as an all-N one of its width, or, under --prune-pam and for a
    3' PAM, wider (a front pad would land inside a 5' PAM); None not at all. The side of
    the guide the index was built for is checked apart."""
    if index_pam == motif:
        return 2
    if not index_pam:
        return None
    if set(index_pam) != {"N"}:
        return 1 if prune and pam_holds(index_pam, motif) else None
    if len(index_pam) == len(motif) or (prune and three_prime and len(index_pam) > len(motif)):
        return 0
    return None


def pam_breadth(pam):
    """How many sequences PAM ``pam`` stands for (NRR 4, NRG 2, NGG 1): of two dedicated
    indexes of PAMs that hold a search's, the search takes the narrower, which stores fewer
    sites, and of two as broad the one whose PAM sorts first (the search script's
    _pam_narrower)."""
    breadth = 1
    for code in pam:
        breadth *= len(_IUPAC_BASES.get(code, ""))
    return breadth
