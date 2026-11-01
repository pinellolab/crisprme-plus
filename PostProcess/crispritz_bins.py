#!/usr/bin/env python
"""Whether the .bin files of a CRISPRitz index are whole.

buildTST (CRISPRitz 2.9.1 and before) saves each tree with unchecked writes, so a
full disk leaves its .bin empty or cut short while the engine still exits 0, and
every search then fails on that index (``numLeaves=``, ``numNodes=0``, or a crash
in searchTST). build-index-only and pool_index_indels.py check each .bin the engine
wrote before they take an index for built.

A .bin is one tree, as saveTST (CRISPRitz's mainParallel.cpp) writes it:

* int numLeaves (at least 1) and int len_guide;
* numLeaves targets: int guideIndex, the target's PAM in ceil(P/2) bytes (two bases
  a byte; P is the length of the PAM that starts the file name,
  ``<PAM>_<chromosome>_<part>.bin``), then ``0``, or ``_`` and the int -(i + 1) of
  the target i with the same sequence inserted before it;
* int numNodes (at least 1);
* the tree, node by node (serialize): each node's base, then for each missing kid
  ``0`` (lo, hi) or ``_`` (eq: the node ends a sequence) followed by the int
  -(i + 1) of the last target i inserted with that sequence; four bits each, two a
  byte, a ``_`` ending its byte.

So with D distinct sequences (the targets whose flag is ``0``) the tree takes
(3 * numNodes + 1 + odd) / 2 bytes and 4 * D more, ``odd`` (the runs between two
``_`` of odd length) being at most D, and it ends with the eq end of the first
sequence inserted, the middle target's (insall), whose chain of earlier targets
ends at that middle target. A file cut short anywhere fails one of these. Ints are
in the byte order of the machine that wrote them, which is this one.

buildTST does not check that a .bin opened either: with the disk's inode table or
a file-count quota full, it writes nothing for the contigs left and exits 0, and
every .bin there is whole. A contig with no target gets no .bin, so missing_bins
tells which contigs without one hold a target, by buildTST's own rule for one.
"""

import array
import os
import re
import struct


def bin_problem(path):
    """What is wrong with the CRISPRitz index file ``path``, or None when it is
    whole."""
    pam = len(os.path.basename(path).split("_", 1)[0])
    width = 4 + (pam + 1) // 2  # a target, up to its flag
    size = os.path.getsize(path)
    if size < 8:
        return "empty" if not size else f"cut short ({size} bytes)"
    with open(path, "rb") as fh:
        leaves = struct.unpack("=i", fh.read(4))[0]
        if leaves <= 0:
            return f"not a tree (numLeaves={leaves})"
        if leaves > (size - 8) // (width + 1):
            return f"cut short in its {leaves} targets"
        fh.seek(8)
        block = fh.read(min(size - 8, leaves * (width + 5) + 4))
        fh.seek(size - 5)
        tail = fh.read(5)
    nexts = array.array("i", bytes(4 * leaves))
    distinct, pos = 0, 0
    try:
        for i in range(leaves):
            flag = block[pos + width]
            if flag == 0x30:  # "0"
                distinct += 1
                pos += width + 1
            elif flag == 0x5F:  # "_"
                nexts[i] = struct.unpack_from("=i", block, pos + width + 1)[0]
                pos += width + 5
            else:
                return f"corrupt (target {i} of {leaves})"
        nodes = struct.unpack_from("=i", block, pos)[0]
    except (IndexError, struct.error):
        return f"cut short in its {leaves} targets"
    if nodes <= 0:
        return f"not a tree (numNodes={nodes})"
    odd = 2 * (size - 8 - pos - 4 - 4 * distinct) - (3 * nodes + 1)
    if not 0 <= odd <= distinct:
        return f"cut short in its tree of {nodes} nodes"
    last = struct.unpack_from("=i", tail, 1)[0]
    i, steps = -last - 1, 0
    while 0 <= i < leaves and nexts[i] and steps < leaves:
        i, steps = -nexts[i] - 1, steps + 1
    if not (tail[0] == 0xF0 or tail[0] & 0x0F == 0x0F) or i != (leaves - 1) // 2 \
            or nexts[i]:
        return f"cut short at the end of its tree of {nodes} nodes"
    return None


def incomplete_bins(directory):
    """``[(path, problem)]`` for each .bin under ``directory``, an index (an INDELS
    index has one directory a chromosome), that is not whole."""
    found = []
    for top, dirs, files in os.walk(directory):
        dirs.sort()
        for name in sorted(files):
            if name.endswith(".bin"):
                path = os.path.join(top, name)
                problem = bin_problem(path)
                if problem:
                    found.append((path, problem))
    return found


# The four bits buildTST gives a base of the genome (genome_bit_conversion, in
# CRISPRitz's pamGenerator.cpp; N 0000) and of the PAM (pam_bit_conversion; N 1111):
# a genome base matches a PAM base when their bits meet, so a genome N matches none.
_BITS = dict(zip("ACGTRYSWKMBDHV", (1, 2, 4, 8, 5, 10, 6, 9, 12, 3, 14, 13, 11, 7)))
_PAM_BITS = dict(_BITS, N=15)
_REVERSE = str.maketrans("ACGTRYMKHDBV", "TGCAYRKMDHVB")  # reversenuc, less the reversal


def holds_target(seq, pam, pam_len, at_start, max_bulges):
    """Whether buildTST (CRISPRitz 2.9.1) indexes a target in contig ``seq`` (upper
    case, each base one of ACGTN and the IUPAC codes), for ``pam`` (the PAM itself,
    at the start of the pam_len bases of a target with ``at_start``, as a 5' PAM),
    with ``max_bulges``: as searchPAMonGenome finds its PAM on either strand,
    guideRegionHasN drops a window of pam_len + max_bulges bases with an N outside its
    PAM, and main takes each kept window (a window at 0 as on the minus strand)."""
    width, size, n = pam_len + max_bulges, len(pam), len(seq)
    if n < width:
        return False
    scanned = n - (pam_len if at_start else size)  # where searchPAMonGenome looks

    def clean(a, b):  # no N in seq[a:b]
        return seq.find("N", a, b) < 0

    def plus(nt):  # a PAM on the plus strand at nt
        if at_start:
            return nt <= n - width and (
                clean(nt + size, nt + width) if nt else clean(0, width - size))
        start = nt + size - width
        return start >= 0 and (clean(start, nt) if start else clean(size, width))

    def minus(nt):  # one on the minus strand
        if at_start:
            return nt + size - width >= 0 and clean(nt + size - width, nt)
        return nt <= n - width and clean(nt + size, nt + width)

    if all(_PAM_BITS[base] == 15 for base in pam):
        # Every position is a PAM on both strands: try those whose window can be
        # clean, from each N-free run of its length, and 0 and the window at 0 (whose
        # N test differs)
        def tried():
            yield 0
            yield width - size
            for m in re.finditer("(?=[^N]{%d})" % (width - size), seq):
                yield m.start() + width - size
                yield m.start() - size
        return any(0 <= nt < scanned and (plus(nt) or minus(nt)) for nt in tried())
    for strand, motif in ((plus, pam), (minus, pam.translate(_REVERSE)[::-1])):
        regex = "(?=%s)" % "".join(
            "[%s]" % "".join(b for b in _BITS if _BITS[b] & _PAM_BITS[base])
            for base in motif)
        for m in re.finditer(regex, seq):
            if m.start() >= scanned:
                break
            if strand(m.start()):
                return True
    return False


def _read_pam(pam_file):
    """(the PAM, pam_len, at_start) as buildTST reads them from the first line of
    ``pam_file``, or None where this module would read it otherwise."""
    try:
        with open(pam_file) as fh:
            line = fh.readline().rstrip("\n").upper()
        motif, _, limit = line.partition(" ")
        limit = int(limit.strip())
    except (OSError, ValueError):
        return None
    if not motif or not 0 < abs(limit) <= len(motif) or set(motif) - set(_PAM_BITS):
        return None
    return (motif[:-limit] if limit < 0 else motif[-limit:]), len(motif), limit < 0


def _read_contig(fasta):
    """(header, sequence) of a FASTA as buildTST reads it: its first line less one
    character, then every whole line after it, upper case; or None where this module
    would read it otherwise (a byte other than an IUPAC base, CRLF line ends, a
    first line that is no header)."""
    try:
        with open(fasta, "rb") as fh:
            lines = fh.read().split(b"\n")
    except OSError:
        return None
    # (getline's last line, with no \n after it, is not taken; nor is the empty
    # string after a last \n)
    head, body = lines[0], b"".join(lines[1:-1]).upper()
    if not head.startswith(b">") or body.translate(None, b"ACGTNRYSWKMBDHV"):
        return None
    return head[1:].decode("utf-8", "surrogateescape"), body.decode("ascii")


def missing_bins(directory, genome_dir, pam_file, max_bulges):
    """``([(FASTA, contig)], [FASTA])``, of the FASTAs of ``genome_dir`` that
    index-genome indexed into ``directory`` (every ``*.fa``) with ``pam_file`` and
    ``max_bulges`` and that have no ``<PAM>_<contig>_<part>.bin`` there: those whose
    contig holds a target (holds_target), and those this module would read otherwise
    than buildTST, which may (a byte other than an IUPAC base, CRLF line ends, a
    first line that is no header, a PAM file it cannot read)."""
    pam = _read_pam(pam_file)
    try:
        bins = os.listdir(directory)
    except OSError:
        bins = []
    named = re.compile(r"%s_(.*)_[0-9]+\.bin\Z" % re.escape(pam[0]) if pam else r"\Z.")
    have = {m.group(1) for m in map(named.match, bins) if m}  # the contigs with a .bin
    found, unread = [], []
    for name in sorted(os.listdir(genome_dir)):
        path = os.path.join(genome_dir, name)
        if name.strip().split(".")[-1] != "fa" or not os.path.isfile(path):
            continue
        try:
            with open(path, "rb") as fh:
                head = fh.readline().rstrip(b"\n")[1:]
            head = head.decode("utf-8", "surrogateescape")
        except OSError:
            unread.append(name)
            continue
        if head in have:
            continue
        contig = _read_contig(path) if pam is not None else None
        if contig is None:
            unread.append(name)
        elif holds_target(contig[1], *pam, max_bulges):
            found.append((name, contig[0]))
    return found, unread
