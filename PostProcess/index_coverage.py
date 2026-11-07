#!/usr/bin/env python
"""Whether an index holds the trees of every chromosome it should.

A search takes an index that holds one .bin as a whole one. One that lacks a chromosome's
trees -- a partial copy, a build killed between chromosomes by an earlier release, or an
index of a genome folder that gained a FASTA since -- loses every site of that chromosome
with exit 0 and no message, from the web and the CLI alike, wherever the search opens the
index (with bulges; a variant one also without, where no whole enriched genome is there to
scan). So the search script checks an index before it searches it, and build-index-only one
it finds "already present".

buildTST writes no tree for a contig that holds no site of the PAM (all N, shorter than one
window, no PAM), so a FASTA without a tree is missing only where it holds a site, by
buildTST's own rule (crispritz_bins.holds_target, by the rule of the index format the index
was built in, as its trees' headers say); a contig with more sites than one tree
takes several, ``<PAM>_<contig>_<1..n>.bin``, every one but the last full, so a gap in
those numbers, or a full last tree of a contig that holds more, is a lost tree too
(crispritz_bins.missing_bins reads a build the same way). The losing of a contig's last
trees when the one before is not full cannot be told from the files.

A reference index answers for the FASTAs of its genome folder. A variant index holds trees
for the chromosomes of its dataset only, so it answers for those: the chromosomes the search
takes from the dataset's VCFs or indel logs, each read in the reference FASTA of its name,
whose every site the variant genome also holds (a variant base widens what it matches, and
adds no N).

Command line, for the search script (nothing on stdout or stderr when the index is whole;
the ERROR on stderr and exit 1 when it is not; a note on stdout and exit 0 where it cannot
tell)::

    python index_coverage.py <index dir> <PAM position of the search's PAM file> <genome dir> \
        [chromosome ...]
"""

import os
import re
import sys

import crispritz_bins
import index_fit


# How problems says that an index holds trees of more than one index format: index-genome
# removes a folder before it builds into it, so no build leaves one, and docs/INDEX_FORMAT.md
# (section 8) has such a folder rebuilt.
_MIXED = "holds trees of index formats"


def _five_prime(index_dir, five_prime_search):
    """Whether ``index_dir`` was built for a 5' PAM: as the header of its trees records
    (index format 2), else as its .pam_build records ("<name> <seq_len> <pam_pos>"), else as
    the search's PAM file (``five_prime_search``) for a dedicated index, and 3' for an all-N
    one (as the search script takes them)."""
    side = index_fit.built_five_prime(index_dir)
    if side is not None:
        return side
    try:
        with open(os.path.join(index_dir, ".pam_build")) as fh:
            return int(fh.readline().split()[-1]) < 0
    except (OSError, ValueError, IndexError):
        pam = os.path.basename(os.path.normpath(index_dir)).split("_")[0]
        return five_prime_search and set(pam) != {"N"}


def _trees(index_dir, pam):
    """{contig: {part numbers}} of the trees of ``index_dir`` named for PAM ``pam``, and the
    PAMs its other trees are named for (index_fit.is_tree: a macOS ._<tree>.bin is none)."""
    named = re.compile(r"%s_(.*)_([0-9]+)\.bin\Z" % re.escape(pam))
    tree = re.compile(r"([^_]+)_.+_[0-9]+\.bin\Z")
    parts, others = {}, set()
    for f in filter(index_fit.is_tree, os.listdir(index_dir)):
        m = named.match(f)
        if m:
            parts.setdefault(m.group(1), set()).add(int(m.group(2)))
        elif tree.match(f):
            others.add(tree.match(f).group(1))
    return parts, others


def problems(index_dir, genome_dir, five_prime_search=False, chroms=None):
    """``[(FASTA, what)]``: the FASTAs of ``genome_dir`` (every ``*.fa``, or those of the
    chromosomes ``chroms`` named ``<chromosome>.fa``) whose trees ``index_dir`` lacks: "no
    tree", or which trees. A FASTA this module cannot read as buildTST reads it, or an index
    whose guide length is not known, is left out (no answer is no refusal)."""
    name = os.path.basename(os.path.normpath(index_dir))
    pam, _, rest = name.partition("_")
    n = rest.partition("_")[0]
    guide = index_fit.built_guide_len(index_dir)
    if not pam or not n.isdigit() or not guide:
        return []
    parts, others = _trees(index_dir, pam)
    if not parts and others:
        return [("", f"its trees are named for {', '.join(sorted(others))}, not {pam}: the "
                     f"folder was renamed")]
    formats = sorted(f for f in index_fit.bin_formats(index_dir) if f is not None)
    if len(formats) > 1:
        return [("", f"{_MIXED} {' and '.join(map(str, formats))}, which no build writes "
                     f"together: it was put together from two builds")]
    fmt = crispritz_bins.index_format(index_dir)
    if chroms is None:
        fastas = sorted(f for f in os.listdir(genome_dir) if f.strip().split(".")[-1] == "fa")
    else:
        fastas = [f"{c}.fa" for c in chroms]
    at_start = _five_prime(index_dir, five_prime_search)
    found = []
    for fasta in fastas:
        path = os.path.join(genome_dir, fasta)
        try:
            with open(path, "rb") as fh:
                head = fh.readline().rstrip(b"\n")[1:].decode("utf-8", "surrogateescape")
        except OSError:
            continue  # (a chromosome of the dataset the genome folder has no FASTA of)
        have = parts.get(head)
        if have and len(have) < max(have):
            found.append((fasta, "trees " + ", ".join(
                str(i) for i in range(1, max(have)) if i not in have) + f" of {max(have)}"))
            continue
        last = os.path.join(index_dir, f"{pam}_{head}_{max(have)}.bin") if have else None
        if have and crispritz_bins._leaves(last) != crispritz_bins.TARGETS_A_TREE:
            continue
        contig = crispritz_bins._read_contig(path)
        if contig is None:
            continue
        if not have:
            if crispritz_bins.holds_target(contig[1], pam, guide + len(pam), at_start, int(n),
                                           fmt):
                found.append((fasta, "no tree"))
            continue
        count = crispritz_bins.count_targets(contig[1], pam, guide + len(pam), at_start, int(n),
                                             fmt)
        if count is not None and count > len(have) * crispritz_bins.TARGETS_A_TREE:
            found.append((fasta, f"trees after {max(have)} (its last tree is full, and it holds "
                                 f"{count} sites)"))
    return found


def vcf_chroms(vcf_dir):
    """The chromosomes of a dataset as the search script takes them from its VCFs: of each
    ``*.vcf.gz`` file name, the last dot-separated part that holds "chr"."""
    chroms = []
    for f in sorted(os.listdir(vcf_dir)):
        if f.endswith(".vcf.gz"):
            named = [part for part in f.split(".") if "chr" in part]
            if named and named[-1] not in chroms:
                chroms.append(named[-1])
    return chroms


def message(index_dir, found, rerun):
    """The ERROR for index ``index_dir`` lacking ``found`` (``problems``), saying how to get
    it whole; ``rerun`` says what to run again once it is removed."""
    name = os.path.basename(os.path.normpath(index_dir))
    if found and not found[0][0] and found[0][1].startswith(_MIXED):
        return (f"index '{name}' ({index_dir}) {found[0][1]}: a tree of format 1 lacks the "
                f"sites at contig ends and N runs that one of format 2 holds, and either build "
                f"may lack trees of the other's FASTAs. Remove it (rm -rf '{index_dir}') and "
                f"{rerun}")
    if found and not found[0][0]:
        what = found[0][1]
    else:
        what = "lacks trees of FASTAs that hold sites of its PAM: " + "; ".join(
            f"{f} ({why})" for f, why in found[:6]) + (
            f"; and {len(found) - 6} more" if len(found) > 6 else "")
    return (f"index '{name}' ({index_dir}) {what}: a "
            f"search of it would lose them with no error. It is a partial copy, or a build "
            f"stopped partway, or the genome folder has changed since it was built. Remove it "
            f"(rm -rf '{index_dir}') and {rerun}")


def main(argv):
    index_dir, pos, genome_dir, chroms = argv[0], int(argv[1]), argv[2], argv[3:] or None
    try:
        found = problems(index_dir, genome_dir, pos < 0, chroms)
    except (OSError, ValueError) as exc:  # (no answer is no refusal)
        print(f"Note: could not check that {index_dir} holds every chromosome's trees ({exc})")
        return 0
    if not found:
        return 0
    sys.stderr.write("ERROR: " + message(
        index_dir, found, "run this search again, which builds it anew where it can, or "
        "download it again or build it with crisprme.py build-index-only") + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
