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
then the tree itself (crispritz_bins has the whole layout). That is index format 1; an index
of format 2 (CRISPRitz 2.10.0 and later, whose trees open with a header) was built since the
fix, whatever its flanks hold.

Command line, for the search script, of a search it prunes by itself (one the web launched,
which never passes --prune-pam): what wider_flank_notes says of the indexes it reads, on
stdout, one line each (nothing on stderr, which fails the stage)::

    python prune_flank.py <PAM> <genome dir> <index dir> ...
"""

import os
import re
import struct
import sys

import index_fit

# More targets than one fake-indel window gives (both strands of its 2 x 30 reference bases,
# 25 and 27 up to CRISPRitz 2.9.1, and the indel's, short in every published dataset): a tree
# that holds more spans several windows, and so an N between them.
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
    # (its trees are named for the PAM, whichever form the index's name takes)
    pam = index_fit.index_pam(os.path.normpath(index_dir))
    trees = _bins(index_dir, pam)
    if trees and index_fit.bin_format(trees[0][1]) == 2:
        return False
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
    if trees and index_fit.bin_format(sorted(trees)[0][1]) == 2:
        return False
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


def wider_flank_notes(indexes, motif, genome_dir, who=None):
    """What a search under --prune-pam, given or taken by itself, should say of the all-N
    indexes ``indexes`` it reads whose stored flank is wider than the PAM ``motif`` (the
    others: nothing), one line each: the sites such a flank cannot index in index format 1
    (METHODS section 9; format 2, which CRISPRitz writes from 2.10.0, holds them), with, under
    a CRISPRitz that reads format 1 alone (2.9.1), the one its search reports on the wrong
    strand (2.10.0's reads it right in either format), and, for an index built before the
    engine kept a target whose flank holds an N, the sites it misses and what serves the
    search whole: that index built again, or a dedicated index of the PAM beside it, which
    the search takes first. A search without --prune-pam
    and --index-path built a dedicated index before 2.10.0; it now prunes from that same index
    by itself. ``who``: "complete-search", which says them on its console for a search of its
    own; None for the search script, which says them in its log for a search it prunes by
    itself, one the web launched, and so names the web app's route beside the command."""
    tag = f" [{who}]" if who else ""
    notes, seen = [], set()
    for index in indexes:
        pam = index_fit.index_pam(index) if index else ""
        if not pam or set(pam) != {"N"} or len(pam) <= len(motif) or index in seen:
            continue
        seen.add(index)
        # (of format 1, or of no format its trees tell, as of an index not there; one of
        # format 2 holds those sites)
        formats = index_fit.bin_formats(index)
        if not notes and (1 in formats or not formats):
            d = len(pam) - len(motif)
            # (a - strand PAM d nt from the start: its window is at index 0, which CRISPRitz
            # 2.9.1's search read as a + strand one)
            misread = (f", and one {d} nt from the start is reported on the + strand at "
                       f"another position" if index_fit.engine_index_format() == 1 else "")
            notes.append(
                f"NOTE{tag}: --prune-pam serves the {len(motif)}-nt PAM {motif} from "
                f"{os.path.basename(index)}, whose stored flank is {len(pam)} nt: a + strand PAM "
                f"with 1 to {d} nt after it to the end of a contig, and a - strand PAM in the "
                f"first {d} nt of a contig, are not indexed{misread} (METHODS section 9). An "
                f"index whose flank is the PAM's width, or a dedicated one, has no such limit, "
                f"nor has one built by CRISPRitz 2.10.0 or later (index format 2).")
        try:
            if os.path.basename(index).count("+"):
                before = indels_built_before_n_fix(index + "_INDELS") \
                    if os.path.isdir(index + "_INDELS") else None
                if before is None:
                    before = built_before_n_fix(index, genome_dir)
            else:
                before = built_before_n_fix(index, genome_dir)
        except OSError:
            before = None
        if before:
            if who:
                serves = (f"Remove it and build it again with build-index-only (which keeps an "
                          f"index it finds in place), or put a dedicated {motif} index beside "
                          f"it, which a search takes before it (build one with build-index-only, "
                          f"or download one if it is published).")
            else:
                serves = (f"Put a dedicated {motif} index beside it, which a search takes before "
                          f"it (build one with crisprme.py build-index-only or in the web app's "
                          f"Settings, or download one if it is published), or remove it and "
                          f"build it again with build-index-only (which keeps an index it finds "
                          f"in place).")
            notes.append(
                f"WARNING{tag}: {os.path.basename(index)} stores no N in the flank of any target "
                f"where its genome has one: it was built before CRISPRitz 2.9.1 (CRISPRme 2.8.0 "
                f"or earlier), which dropped such targets, and searched for the shorter PAM "
                f"{motif} it misses the sites whose flank holds an N. {serves}")
    return notes


if __name__ == "__main__":
    for _note in wider_flank_notes(sys.argv[3:], sys.argv[1], sys.argv[2]):
        print(_note)
