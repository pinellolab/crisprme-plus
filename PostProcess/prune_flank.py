#!/usr/bin/env python
"""Whether an all-N index was built by an engine that dropped a target whose flank holds an N.

Under --prune-pam a search for a PAM of L nt reads an all-N index whose stored flank is K nt.
CRISPRitz before 2.9.1 (CRISPRme 2.8.0 and earlier) dropped every target whose PAM/flank region
held an N, so such an index of K > L misses the sites whose flank holds an N past the L bases
the PAM covers, which an index of the PAM's own width, or a dedicated one, finds: 3 of 36 rows
on the synthetic genome, 2.6% of the indel rows on chr22 (idx audit IDX-5). The index records
nothing of the engine that built it (.pam_build has the same format), but its trees tell: the
engine since 2.9.1 stores such a target with an N in its flank, a nibble of 0 (an N matches no
base), which no earlier build holds. So a tree of a contig whose sequence holds an N beside its
bases holds a 0 nibble in a flank when the index was built since, and none when before. The
fake-indel genome of an INDELS index is N-separated, one window an indel, so each of its trees
that holds more targets than one window gives does too.

A tree (.bin) is int numLeaves and int len_guide_used, then each target: int guideIndex, its
flank in ceil(K/2) bytes, two bases a byte, high nibble first (for odd K the last byte's high
nibble is padding, always 0, and its low nibble the last base), then '0', or '_' and an int;
then the tree itself (crispritz_bins has the whole layout).
"""

import os
import re
import struct

# More targets than one fake-indel window gives (both strands of its 2 x ~25 bases and the
# indel's, short in every published dataset): a tree that holds more spans several windows, and
# so an N between them.
ONE_WINDOW = 50000


def flank_holds_n(bin_path, k):
    """Whether a target of the tree ``bin_path`` (an index of flank width ``k``) stores an N in
    its flank; None where the file cannot be read as a tree."""
    width = (k + 1) // 2
    padding = 2 * width - 2 if k % 2 else None  # (the nibble that holds no base)
    try:
        with open(bin_path, "rb") as fh:
            head = fh.read(8)
            if len(head) < 8:
                return None
            leaves = struct.unpack("=i", head[:4])[0]
            if leaves <= 0:
                return None
            for _ in range(leaves):
                target = fh.read(4 + width + 1)
                if len(target) < 4 + width + 1:
                    return None
                for i, byte in enumerate(target[4:4 + width]):
                    if (byte >> 4 == 0 and 2 * i != padding) or byte & 15 == 0:
                        return True
                flag = target[-1:]
                if flag == b"_":
                    fh.seek(4, os.SEEK_CUR)
                elif flag != b"0":
                    return None
    except OSError:
        return None
    return False


def _bins(index_dir, pam):
    """[(size, path, contig, part)] of the trees of ``index_dir``, smallest first."""
    named = re.compile(r"%s_(.*)_([0-9]+)\.bin\Z" % re.escape(pam))
    found = []
    for f in os.listdir(index_dir):
        m = named.match(f)
        if m:
            path = os.path.join(index_dir, f)
            found.append((os.path.getsize(path), path, m.group(1), int(m.group(2))))
    return sorted(found)


def _holds_n_and_bases(fasta, limit=1 << 26):
    """Whether FASTA ``fasta`` holds an N and a base beside it, read up to ``limit`` bytes of
    sequence; None where it cannot be read."""
    seen_n = seen_base = False
    try:
        with open(fasta) as fh:
            fh.readline()
            read = 0
            for line in fh:
                line = line.strip().upper()
                read += len(line)
                if "N" in line:
                    seen_n = True
                if line.strip("N"):
                    seen_base = True
                if (seen_n and seen_base) or read > limit:
                    break
    except (OSError, UnicodeDecodeError):
        return None
    return seen_n and seen_base


def built_before_n_fix(index_dir, genome_dir):
    """True when all-N reference index ``index_dir`` holds no N in any flank of the smallest
    one-tree contig of ``genome_dir`` whose sequence holds an N beside its bases; False when it
    holds one there; None where no such contig tells."""
    pam = os.path.basename(os.path.normpath(index_dir)).split("_")[0]
    trees = _bins(index_dir, pam)
    parts = {}
    for _size, _path, contig, part in trees:
        parts[contig] = max(parts.get(contig, 0), part)
    for _size, path, contig, _part in trees:
        if parts[contig] != 1:
            continue  # (a contig of several trees shares its targets out among them)
        fasta = os.path.join(genome_dir, f"{contig}.fa")
        if not os.path.isfile(fasta) or not _holds_n_and_bases(fasta):
            continue
        holds = flank_holds_n(path, len(pam))
        if holds is not None:
            return not holds
    return None


def indels_built_before_n_fix(indels_dir):
    """The same of an INDELS index ``indels_dir`` (one sub-index a chromosome), by the smallest
    tree of a fake-indel contig that spans more than one window."""
    trees = []
    for sub in os.listdir(indels_dir):
        folder = os.path.join(indels_dir, sub)
        if os.path.isdir(folder):
            trees += _bins(folder, sub.split("_")[0])
    for _size, path, _contig, _part in sorted(trees):
        try:
            with open(path, "rb") as fh:
                leaves = struct.unpack("=i", fh.read(4))[0]
        except (OSError, struct.error):
            continue
        if leaves <= ONE_WINDOW:
            continue
        holds = flank_holds_n(path, len(os.path.basename(path).split("_")[0]))
        if holds is not None:
            return not holds
    return None
