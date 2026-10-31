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
"""

import array
import os
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
