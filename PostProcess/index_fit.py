"""Whether an index, whoever built it, can answer a search: in an index format this CRISPRme
and the installed CRISPRitz read, and of a guide of a given length.

An index folder named ``<PAM>_<N>_<genome>``, as every index was before CRISPRme 2.10.0,
says nothing of the guide it was built for (the name of 2.10.0 does: parse_index_name). The
engine reads that length from each tree (``len_guide_used``, the second int32 of a format-1
.bin, written by buildTST since 2019, whatever built the index; the header's guide_len in
format 2), but only to offset coordinates: it searches any guide it is given, and past what
the tree stores it reports wrong alignments or loses targets, at exit 0.

The rule, measured through the real CLI against an index dedicated to the guide (246 cells
of built length, N, guide length and bulges, dedicated, all-N and all-N under --prune-pam,
and 70 more of variant searches and a 5' PAM; submit_job_automated_new_multiple_vcfs.sh
_index_fits_guide says the same, and the two must not disagree):

    built + N >= guide + bDNA + (1 if bRNA else 0)

Inside it the results are byte-identical to a dedicated index; outside it, in every cell,
they differ. So a 21-nt guide is served at up to one bulge of each type by any index built
for 20 nt with N=3, and a 23-nt guide is served by none of the published ones. (On an index
of format 2 the rows inside the rule are those of a dedicated index but at contig ends and
N runs, where an index built for the shorter guide holds sites one built for the longer
guide cannot: CRISPRitz's docs/INDEX_FORMAT.md, section 5.1.)

Index formats (CRISPRitz's docs/INDEX_FORMAT.md): CRISPRitz up to 2.9.1 writes format 1, a
.bin with no header; 2.10.0 and later write format 2, a .bin that opens and ends with a
magic and holds in its header the parameters it was built with, and keeps the sites at
contig ends and N runs that format 1 cannot hold. A CRISPRitz reads every format up to the
one it writes (its crispritz.py's INDEX_FORMAT; none up to 2.9.1), so an index is searched
only where both this CRISPRme (READ_FORMATS) and the installed CRISPRitz read its format.
"""

import json
import os
import re
import shutil
import struct
import sys

# The first four bytes, and the last four, of an index file that CRISPRitz 2.10.0 and later
# write (index format 2: CRISPRitz's docs/INDEX_FORMAT.md): the magic, then the format
# version as an int32. A format-1 file (CRISPRitz up to 2.9.1) opens with its leaf count, a
# positive int32, where the magic reads as a negative one.
CZI_MAGIC = b"CZI\xfe"
# The index formats this CRISPRme reads
READ_FORMATS = (1, 2)
# The header flags of format 2 (CRISPRitz's include/indexFormat.h): the PAM at the 5' end of
# the guide, an all-N build PAM, and the four-allele symbol Q in the tree or the flanks. A
# reader refuses a file with another bit set, as searchTST does.
FLAG_PAM_AT_START = 1
FLAG_PAMLESS = 2
FLAG_FOUR_ALLELE = 4
_FLAGS_KNOWN = FLAG_PAM_AT_START | FLAG_PAMLESS | FLAG_FOUR_ALLELE
# The record of what an index holds, in its own folder (write_index_meta). A dot file:
# crispritz.py search refuses a folder that holds any other file but a .bin, with an ERROR
# on stdout and no targets file, at exit 0 (CRISPRitz 2.9.1 and 2.10.0 alike).
INDEX_META = ".index_meta.json"


def read_index_header(path):
    """The index format and build parameters of the CRISPRitz .bin ``path``, read from
    its header, the two ints that open its body and its last 4 bytes, without loading
    the tree (docs/INDEX_FORMAT.md section 7). An "error" entry says why searchTST
    refuses the file: every check searchTST makes on a file's own bytes and name
    before it reads its leaves is made here, in the same order. A format-1 file has no
    end marker, so one cut short after its first 8 bytes is not seen here.

    CRISPRitz's own reader (crispritz.py read_index_header, which index-info runs), as
    CRISPRitz 2.10.0 has it, with the header flags it knows (_FLAGS_KNOWN): the test of
    this module, and test/index_format2/check_real_engine.py against the installed
    CRISPRitz, check the two return the same for the same files."""

    def refused(info, why):
        info["error"] = why
        return info

    def damaged(info, header_len, width, size):
        # a leaf record takes at least 4 + ceil(K/2) + 1 bytes (window start, flank
        # nibbles, marker or flags byte), so a leaf count the file cannot hold is damage
        need = header_len + 12 + info["num_leaves"] * (4 + (width + 1) // 2 + 1)
        if need > size:
            return refused(info, "damaged: its %d leaves need at least %d bytes, the file holds %d"
                           % (info["num_leaves"], need, size))
        return info

    # the PAM the file is named for, <PAM>_<contig>_<part>.bin: searchTST reads every
    # .bin of a directory at the width of that PAM
    named = os.fsencode(os.path.basename(path)).split(b"_", 1)[0]
    try:
        fh = open(path, "rb")
    except OSError as exc:
        return refused({"format": None}, "cannot be opened (%s)" % exc.strerror)
    with fh:
        size = os.fstat(fh.fileno()).st_size
        head = fh.read(28)
        if head[:4] != CZI_MAGIC:
            # format 1 (CRISPRitz <= 2.9.1): no header, the file opens with numLeaves
            # and len_guide
            if len(head) < 8:
                return refused({"format": None}, "too short to be a CRISPRitz index (.bin) file")
            leaves, guide = struct.unpack("=ii", head[:8])
            info = {"format": 1, "num_leaves": leaves, "guide_len": guide}
            if leaves <= 0:
                return refused(info, "not a valid CRISPRitz index (.bin) file (numLeaves=%d)" % leaves)
            return damaged(info, 0, len(named), size)
        if len(head) < 28:
            return refused({"format": None}, "cut-short index header")
        version, header_len, flags, guide, pam_len, mb = struct.unpack("=6i", head[4:28])
        info = {
            "format": version if version >= 2 else None,
            "pam_at_start": bool(flags & FLAG_PAM_AT_START),
            "pamless": bool(flags & FLAG_PAMLESS),
            "four_allele": bool(flags & FLAG_FOUR_ALLELE),
            "guide_len": guide,
            "pam_len": pam_len,
            "max_bulges": mb,
        }
        if version > 2:
            return refused(info, "index format %d, newer than this CRISPRitz reads (formats 1 and "
                                 "2): update CRISPRitz, or rebuild the index with this one" % version)
        # max_bulges fits pad_valid's 5 bits
        if (version != 2 or flags & ~_FLAGS_KNOWN or not 1 <= pam_len <= 1000
                or header_len != 28 + pam_len or guide < 1 or not 0 <= mb <= 31):
            return refused(info, "malformed index header (version %d, flags %d, PAM length %d, "
                                 "header length %d, guide length %d, max bulges %d)"
                           % (version, flags, pam_len, header_len, guide, mb))
        pam = fh.read(pam_len)
        info["pam"] = pam.decode("ascii", "replace")
        if len(pam) < pam_len:
            return refused(info, "cut-short index header")
        fh.seek(-4, os.SEEK_END)
        info["complete"] = fh.read(4) == CZI_MAGIC
        if not info["complete"]:
            return refused(info, "cut short: it does not end in the end-of-index marker")
        fh.seek(header_len)
        body = fh.read(8)
    if len(body) < 8:
        return refused(info, "cut short: no leaf count and guide length after its header")
    leaves, body_guide = struct.unpack("=ii", body)
    info["num_leaves"] = leaves
    if body_guide != guide:
        return refused(info, "malformed index: its header gives guide length %d, its body %d"
                       % (guide, body_guide))
    # searchTST reads a format-2 file only if its header PAM is its name's prefix
    if pam != named:
        return refused(info, "holds the PAM/flank region '%s' but is named for '%s'"
                       % (info["pam"], os.fsdecode(named)))
    if leaves <= 0:
        return refused(info, "not a valid CRISPRitz index (.bin) file (numLeaves=%d)" % leaves)
    return damaged(info, header_len, pam_len, size)


def bin_format(path):
    """The index format of the CRISPRitz .bin ``path``: 1 for a file with no header (up
    to CRISPRitz 2.9.1), else the version that follows the magic (2 from CRISPRitz 2.10.0).
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


def is_tree(name):
    """Whether ``name``, a file of an index folder, is named as one of its trees,
    ``<PAM>_<contig>_<part>.bin``, as the search script's links to them are (_bin_only_view):
    a dot file is none. macOS tar and cp leave an AppleDouble file ``._<name>`` beside each
    file they copy onto a volume without extended attributes; ``._NGG_chr22_1.bin`` opens
    with no magic, so it reads as a format-1 tree built for a 512-nt guide (its second
    int32), and searchTST stops on a folder that holds one (CRISPRitz's docs/INDEX_FORMAT.md,
    section 5)."""
    return name.endswith(".bin") and not name.startswith(".")


_engine_cache = {}


def _engine_value(name):
    """The value of the ``<name> = ...`` line of the crispritz.py on PATH, the engine this
    install runs, read from the file rather than run, so that nothing reaches stderr (fatal
    to a pipeline stage), as generate_report reads its VERSION; "" where that file has no
    such line, None where there is no crispritz.py on PATH."""
    exe = shutil.which("crispritz.py")
    if not exe:
        return None
    try:
        key = (exe, os.stat(exe).st_mtime_ns)
        if key not in _engine_cache:
            values = {}
            with open(exe, errors="replace") as fh:
                for line in fh:
                    m = re.match(r"""(VERSION|INDEX_FORMAT|INDEL_LEFT|INDEL_RIGHT)\s*=\s*"""
                                 r"""["']?([^"'\s#]+)""", line)
                    if m and m.group(1) not in values:
                        values[m.group(1)] = m.group(2)
            _engine_cache[key] = values
        return _engine_cache[key].get(name, "")
    except OSError:
        return None


def engine_version():
    """The VERSION of the installed CRISPRitz, or None where it cannot be told."""
    return _engine_value("VERSION") or None


def engine_index_format():
    """The newest index format the installed CRISPRitz reads, and the one its buildTST writes:
    the INDEX_FORMAT of its crispritz.py, which CRISPRitz states from 2.10.0; 1 for one that
    states none (up to 2.9.1, which read and write format 1 only); None with no crispritz.py
    on PATH."""
    value = _engine_value("INDEX_FORMAT")
    if value is None:
        return None
    return int(value) if value.isdigit() else 1


def readable_formats():
    """The index formats this install searches: those this CRISPRme reads that the installed
    CRISPRitz reads too, all of READ_FORMATS where no crispritz.py is on PATH to tell."""
    engine = engine_index_format()
    return tuple(f for f in READ_FORMATS if engine is None or f <= engine)


def unreadable_bin(index_dir, every=False):
    """(path, format) of the first tree of ``index_dir`` (is_tree), or of a folder in it (an
    INDELS index holds one a chromosome), in an index format this install does not search
    (readable_formats); None where there is none. The first tree of each folder is read, as
    the search script's _unreadable_bin reads it: one build writes a folder whole (a folder of
    two builds is refused as such, index_coverage); with ``every``, every tree, as a build
    just made is checked."""
    readable = readable_formats()
    try:
        for top, dirs, files in os.walk(index_dir):
            dirs.sort()
            bins = sorted(name for name in files if is_tree(name))
            for name in bins if every else bins[:1]:
                fmt = bin_format(os.path.join(top, name))
                if fmt is not None and fmt not in readable:
                    return os.path.join(top, name), fmt
    except OSError:
        pass
    return None


def unreadable_format_error(index_dir):
    """The ERROR for an index ``index_dir`` in an index format this CRISPRme, or the installed
    CRISPRitz, does not read (unreadable_bin), naming its format and what reads it, or None.
    Searched, its format version would be read as a guide length (a format-2 index by
    CRISPRme 2.8.1 as one built for a 2-nt guide), or the engine would refuse it, or abort."""
    found = unreadable_bin(index_dir)
    if found is None:
        return None
    path, fmt = found
    held = (f"{index_dir} was built by a newer CRISPRitz ({os.path.relpath(path, index_dir)} "
            f"is in index format {fmt}")
    if fmt not in READ_FORMATS:
        return (f"{held}): this CRISPRme reads index formats 1 and 2 (which CRISPRitz writes up "
                f"to 2.9.1 and from 2.10.0), and needs a newer release to search it. Use a "
                f"CRISPRme that reads index format {fmt}, or build the index again with the "
                f"CRISPRitz installed here, with {index_dir} moved out of the way first")
    engine = engine_version()
    return (f"{held}, which CRISPRitz 2.10.0 and later write) than the one installed here"
            f"{f', CRISPRitz {engine}' if engine else ''}, which reads index format "
            f"{engine_index_format()} only. Install the CRISPRitz this CRISPRme installs "
            f"with (2.10.0 or later; the Docker image, install_from_source.sh and the Bioconda "
            f"package carry it), or build the index again with the CRISPRitz installed here, "
            f"with {index_dir} moved out of the way first")


def _first_bin_header(index_dir):
    """read_index_header of the first tree of ``index_dir`` (is_tree), or None where it has
    none."""
    try:
        for name in sorted(os.listdir(index_dir)):
            if is_tree(name):
                return read_index_header(os.path.join(index_dir, name))
    except OSError:
        pass
    return None


# What read_index_header says of a format-2 file whose header it could not read; any other
# refusal leaves the header's build parameters read (searchTST then refuses the file, as it
# does a format-1 file whose leaves are cut short, whose guide length is read all the same)
_HEADER_UNREAD = ("cut-short index header", "malformed index header", "index format ")


def _readable_header(info):
    """Whether ``info`` (read_index_header) holds the build parameters of a header read."""
    return (info is not None and info.get("format") in READ_FORMATS
            and not info.get("error", "").startswith(_HEADER_UNREAD))


def built_guide_len(index_dir):
    """The guide length ``index_dir`` was built for, from the header of a .bin of it: in
    format 1 the second int32, in the byte order of the machine that reads it, as the engine
    reads it; in format 2 the header's guide_len. None where it cannot tell: no .bin, one
    too short, or one in an index format this CRISPRme does not read (unreadable_bin), whose
    second int32 is its format version."""
    info = _first_bin_header(index_dir)
    if not _readable_header(info):
        return None
    built = info.get("guide_len")
    return built if built and built > 0 else None


def built_five_prime(index_dir):
    """Whether ``index_dir`` was built for a PAM at the 5' end of the guide, as the header of
    its .bin files records it (format 2); None where it records no side (format 1 holds none:
    the index's .pam_build may) or cannot be read."""
    info = _first_bin_header(index_dir)
    if not _readable_header(info) or info["format"] != 2:
        return None
    return info["pam_at_start"]


def bin_formats(index_dir):
    """{format: number of .bin files} of the trees of ``index_dir`` itself (is_tree), from
    their first eight bytes (bin_format; None for one too short to tell)."""
    found = {}
    try:
        names = sorted(os.listdir(index_dir))
    except OSError:
        return found
    for name in names:
        if is_tree(name):
            fmt = bin_format(os.path.join(index_dir, name))
            found[fmt] = found.get(fmt, 0) + 1
    return found


def read_index_meta(index_dir):
    """The record write_index_meta keeps in ``index_dir`` (a dict), or None."""
    try:
        with open(os.path.join(index_dir, INDEX_META)) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return None
    return meta if isinstance(meta, dict) else None


def recorded_format(index_dir):
    """The index format ``index_dir`` records it was built in (read_index_meta), or None."""
    fmt = (read_index_meta(index_dir) or {}).get("index_format")
    return fmt if isinstance(fmt, int) and not isinstance(fmt, bool) else None


def _crisprme_version():
    """This CRISPRme's version, as crisprme.py beside its PostProcess folder states it."""
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        with open(os.path.join(os.path.dirname(here), "crisprme.py"), errors="replace") as fh:
            for line in fh:
                m = re.match(r"""version\s*=\s*["']([^"']+)["']""", line)
                if m:
                    return m.group(1)
    except OSError:
        pass
    return None


def write_index_meta(index_dir, pam_file=None, built=False, indel_window=None):
    """Records in ``index_dir`` (INDEX_META, JSON) what the index holds: its index format,
    build PAM, the side of the guide that PAM is on, the flank width (the build PAM's
    length), the guide length and max bulges it was built for, whether its build PAM is all
    N and whether it holds the four-allele symbol Q, and, for the sub-index of an INDELS
    index, ``indel_window``: the (left, right) reference bases of the fake-indel windows it
    was built from, which its indel log tells (indel_window.META_KEY; a search refuses the
    sub-index beside a log of other windows). Each is read from the headers of its
    .bin files where they hold it (format 2 holds them all; format 1 the guide length, its
    file names the PAM); else from ``pam_file``, the PAM file of the build, and the N of the
    folder's name, and the format, for an index with no .bin, from the installed CRISPRitz
    (the one that just built it) -- only with ``built``, for an index this run has just
    built: of one found in place nothing is taken that the index does not hold itself (a
    side its .pam_build records excepted, which buildTST format 1 does not store). With
    ``built``, the versions of the CRISPRitz and CRISPRme that built it are recorded too;
    without it, an index that has a record keeps it. Returns None, or why the record could
    not be written (an OSError's text); nothing is written for an index whose .bin files
    differ in format, or that has none and was not built here (no format to record)."""
    path = os.path.join(index_dir, INDEX_META)
    if not built and os.path.isfile(path):
        return None  # (the build's own record, with what built it)
    meta = {}
    formats = bin_formats(index_dir)
    if len(formats) > 1:
        return None
    flags = set()
    if formats:
        fmt = next(iter(formats))
        if fmt not in READ_FORMATS:
            return None
        meta["index_format"] = fmt
        for name in sorted(os.listdir(index_dir)):
            if not is_tree(name):
                continue
            info = read_index_header(os.path.join(index_dir, name))
            if fmt == 2 and "error" not in info:
                if "build_pam" not in meta:
                    meta["build_pam"] = info["pam"]
                    meta["side"] = "5'" if info["pam_at_start"] else "3'"
                    meta["guide_len"] = info["guide_len"]
                    meta["max_bulges"] = info["max_bulges"]
                flags.add(bool(info.get("four_allele")))
            elif fmt == 1 and "build_pam" not in meta:
                meta["build_pam"] = name.split("_", 1)[0]
                if info.get("guide_len", 0) > 0:
                    meta["guide_len"] = info["guide_len"]
            if fmt == 1 and "build_pam" in meta:
                break
    elif built:
        meta["index_format"] = engine_index_format() or 1
    else:
        return None
    seq, pos = None, None
    if pam_file:
        try:
            with open(pam_file) as fh:
                parts = fh.readline().split()
            seq, pos = parts[0].upper(), int(parts[1])
        except (OSError, ValueError, IndexError):
            seq = None
    if built and seq and pos:
        meta.setdefault("build_pam", seq[:-pos] if pos < 0 else seq[-pos:])
        meta.setdefault("side", "5'" if pos < 0 else "3'")
        meta.setdefault("guide_len", len(seq) - abs(pos))
    if "side" not in meta:
        try:
            with open(os.path.join(index_dir, ".pam_build")) as fh:
                recorded = int(fh.readline().split()[-1])
            if recorded:
                meta["side"] = "5'" if recorded < 0 else "3'"
        except (OSError, ValueError, IndexError):
            pass
    name = os.path.basename(os.path.normpath(index_dir))
    n = name.split("_")[1] if name.count("_") >= 2 else ""
    if "max_bulges" not in meta and n.isdigit():
        meta["max_bulges"] = int(n)
    if "build_pam" in meta:
        meta["flank_len"] = len(meta["build_pam"])
        meta["pamless"] = set(meta["build_pam"]) == {"N"}
    meta["four_allele"] = True in flags
    if indel_window is not None:
        meta["indel_window"] = [int(indel_window[0]), int(indel_window[1])]
    if built:
        meta["built_by"] = {"crispritz": engine_version(), "crisprme": _crisprme_version()}
    try:
        with open(path + ".part", "w") as fh:
            json.dump(meta, fh, indent=1, sort_keys=True)
            fh.write("\n")
        os.replace(path + ".part", path)
    except OSError as exc:
        try:
            os.remove(path + ".part")
        except OSError:
            pass
        return f"could not write {path} ({exc})"
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


def parse_index_name(name):
    """The parts of an index folder's (or archive's) name, ``(pam, n, rest, guide, side)``,
    or None where the name is not an index's.

    An index is named ``<PAM>_<N>_<genome>[+<dataset>]`` (``NRR_3_hg38+hg38_HGDP1kGP``),
    which says nothing of the guide (spacer) length it was built for, nor of the side of the
    guide its PAM is on: indexes of one PAM for two guide lengths, or for a 3' and a 5' PAM
    (``N...NTTTV 4`` and the Cas12a ``TTTV`` file both build ``TTTV_<N>_<genome>``), take one
    name. Since CRISPRme 2.10.0 a name carries both, as builds write it (index_name):
    ``<guide>-<PAM>_<N>_...`` for a 3' PAM (``20-NRR_3_hg38``), ``<PAM>-<guide>_<N>_...`` for
    a 5' one (``TTTV-23_3_hg38``), in the order of the target on the strand. ``guide`` and
    ``side`` (3 or 5) are None for an older name. Every reader takes both forms, and of two
    indexes that serve a search alike, the newer form (the search script's
    _index_name_parts, the same rule). ``rest`` is the name
    after ``<N>_``: the genome, and ``+<dataset>`` for a variant index."""
    token, sep, rest = os.path.basename(str(name)).partition("_")
    n, sep2, rest = rest.partition("_")
    if not sep or not sep2 or not token or not n.isdigit():
        return None
    pam, guide, side = token, None, None
    head, dash, tail = token.partition("-")
    if dash and head.isdigit() and tail and "-" not in tail:
        pam, guide, side = tail, int(head), 3
    elif dash and tail.isdigit() and head:
        pam, guide, side = head, int(tail), 5
    return pam, int(n), rest, guide, side


def index_pam(name):
    """The PAM of the index named ``name`` (either form: parse_index_name), or the name's
    first field where it is no index name."""
    parsed = parse_index_name(name)
    return parsed[0] if parsed else os.path.basename(str(name)).split("_")[0]


def index_name(pam, n, rest, guide, three_prime):
    """The name of an index of PAM ``pam`` built with -bMax ``n`` for ``rest`` (the genome,
    and ``+<dataset>`` for a variant index), for a guide of ``guide`` nt with the PAM on its
    3' side (``three_prime``) or its 5' side: the form of CRISPRme 2.10.0, which carries both
    (parse_index_name), under which build-index-only and a search build an index, so that
    the indexes of one PAM for the two sides, or for two guide lengths, take a name each;
    the older ``<PAM>_<N>_<rest>`` where the guide length is not known (0 or None: a PAM
    file of its PAM alone). The search script's _build_name, the same rule."""
    if not guide:
        return f"{pam}_{n}_{rest}"
    return f"{guide}-{pam}_{n}_{rest}" if three_prime else f"{pam}-{guide}_{n}_{rest}"


def main(argv):
    """For the search script, which records each index it builds as build-index-only does:

        python index_fit.py meta <index dir> <PAM file of the build>

    writes the index's record (write_index_meta, as built here); a record it cannot write
    is said on stdout (stderr is fatal to the stage), and the status is 0 either way."""
    if len(argv) == 3 and argv[0] == "meta":
        why = write_index_meta(argv[1], argv[2], built=True)
        if why:
            print(f"Note: {why}")
        return 0
    sys.stderr.write("usage: index_fit.py meta <index dir> <PAM file>\n")
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
