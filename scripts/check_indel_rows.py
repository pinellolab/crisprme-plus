#!/usr/bin/env python3
"""An independent check of CRISPRme's indel (and MNV) off-target rows on real data.

    check_indel_rows.py --results <results dir> --reference <genome dir or FASTA>
                        --vcf <bgzipped VCF> [--vcf ...] | --vcf-dir <dir>
                        --pam <PAM file> --guide <guide or guide file>
                        --mm 4 --bDNA 1 --bRNA 1 [--max-total-edits 6]
                        [--chrom chr22 | --region chr22:1-500000]
                        [--mode soundness|completeness|both] [--snp-rows report|off]
                        [--jobs N] [--known known.tsv] [--out failures.tsv]

WHAT IT IS. A second opinion on the indel path, built only from what a user has: the
reference FASTA, the VCF(s) the search read, the guide, the PAM file, the search budget
and the result folder. It imports nothing but Python's standard library and no code of
CRISPRme's variant or indel machinery (``PostProcess/indel_snp_cis.py``,
``analisi_indels_NNN.py``, ``overlay_indel_snps*``, ``resultIntegrator``,
``remove_contiguous_samples``, the CRISPRitz enricher, the dictionaries, the indel log or
the fake contigs), and none of it is copied here. It reads the VCF with its own BGZF +
tabix reader, writes each sample's haplotypes itself from the VCF 4.x definition of a
record, keeps its own map from haplotype bases to reference coordinates and aligns the
guide with its own exhaustive aligner. The CFD tables CRISPRme ships
(``PostProcess/cfd_*_scores.json``) are read as data.

``crisprme.py validate-test`` compares the engine's SITES (chromosome, start, strand) on
the reference and on the SNV-only 1000 Genomes chr22 genome, so an indel row's
coordinate, its aligned strings, its REF/ALT text, its variants, its carriers and the
alternative-alignment table are never checked against an independent truth on real data.
This is that check.

THE TRUTH IS THE VCF. A sample's haplotype is the reference with every ALT allele its
genotype puts on that slot, each in place of its REF bases (VCF 4.x). Two alleles of one
slot whose REF spans meet (a SNP at an insertion's own position, a SNP on a base a
deletion removes, a complex allele split into overlapping records, two ALTs at one
position phased onto one haplotype) describe one stretch twice: the VCF gives no sequence
there, and the slot is counted under its own class
(limit-overlapping-records) wherever a target touches it. A sample's alleles are on known
slots when its calls are phased ('|', or haploid); where a call that touches the target
is unphased or has a missing allele, which of its alleles share the record's haplotype is
unknown, and the truth is then CRISPRme's documented PUTATIVE rule (METHODS section 4):
the sample is listed for the strongest target its own alleles allow -- the guide's base
wherever it carries it; elsewhere, and at the PAM jointly, the alleles of highest CFD
among those the motif allows -- and never for an allele it does not carry.

THE DOCUMENTED VIEW ONLY CLASSIFIES. Where CRISPRme's output and the truth part, the
rules CRISPRme documents for reading a carrier (METHODS section 4, CHANGELOG 2.8.1 Known
limits), re-implemented here from that text, say whether the difference is a documented
limit (counted and listed under its own class, non-gating) or not (gating):
  * the fake-indel window holds one indel-path record (an indel or an MNV) and the SNPs of
    the enriched genome, which keeps only the last of several SNP records at one position
    (limit-last-alt-wins);
  * a second phased indel inside the target is read on the slot's bases, from the
    target's first base on or up to its last base, through the record; an unreadable or
    different reading takes the slot off, and a target that needs two indel-path records
    is never reported (limit-two-indels);
  * another allele at the record's own position is a phasing artefact, a variant only on
    bases the record deletes changes nothing, and an overlapping indel at another position
    takes the slot off (limit-overlapping-records; a listed carrier the rule takes off is
    carrier-off-by-rule, gating);
  * the carriers of an MNV, of an alt of a multi-ALT SNP line whose number the genotypes
    cannot tell, and of a second indel with an unphased call are resolved PUTATIVE
    (limit-putative-carrier where that differs from their known haplotype).
  A disagreement no documented rule explains gates (unexplained-disagreement: a defect of
  CRISPRme or of this checker).

TWO MODES.

  SOUNDNESS (every reported indel/MNV row is real). Each row block of the integrated
  results (highest-CFD, fewest mm+b, highest CRISPR-Bulge) and of the alternative
  alignments whose variant list names an indel or an MNV:
    * is labelled REF/ALT_origin alt, names exactly one indel-path record (an indel or an
      MNV) and SNPs, every one a record of the VCF, forward strand, on the row's
      chromosome;
    * its aligned spacer, gaps removed, is the guide (--guide, never the table's own
      column) followed (or, for a 5' PAM, preceded) by the PAM; its mismatches, bulge
      types and sizes and their total recount exactly, within the budget; no gap falls in
      the PAM; the PAM matches the searched motif; no IUPAC code is left;
    * its target (gaps removed, PAM included) is the reference with exactly the named
      variants applied, read from the reported start under the convention below, and
      holds a base of its indel; every named SNP changes a base of the target;
    * its REF text is the reference read from the start, gapped as the target;
    * every listed carrier holds it, and no carrier is missing (CARRIERS).

  COMPLETENESS (no real indel off-target is missed -- the real-data brute force). For
  every indel and MNV record of the VCF in the checked region: each carrier's haplotypes
  around it are written (its phased slots; a call of unknown phase is tried both ways),
  identical sequences are aligned once, every in-budget alignment on both strands with the
  PAM that holds a base the record changes (the junction rule) is enumerated and mapped to
  the reference, and each (sample, target) must be reported for that sample at that site
  -- as an integrated row or an alternative alignment -- or merged into a row of the same
  site by the documented merge. A target that touches a call of unknown phase is demanded
  as the strongest target the sample's alleles allow there (the PUTATIVE rule). A run that
  writes no alternative alignments keeps one row per merge cluster -- a population-level
  run, or one made with --ci-cd-test, as the documented chr22 self-test (.alt_mode =
  no-alt-alignments) -- so there each target's cluster must be reported with a
  fewest-mm+b block no worse than it and a highest-CFD block scoring at least its CFD.

DEFINITIONS (fix/282-indel-alignment; METHODS section 4, "SNP+indel off-targets", and
the merge of ``remove_contiguous_samples.py``).

  Coordinate. A row starts at the 0-based reference coordinate of its target's FIRST
  ALIGNED BASE, the leftmost on either strand. A base an insertion adds has none: such a
  row starts at the reference base the added bases follow, the insertion's ANCHOR (its
  VCF POS - 1, 0-based). An ALT allele's bases are laid on the reference after the bases
  REF and ALT share at either end: the j-th remaining ALT base stands in the place of the
  j-th remaining REF base while there is one, and the rest are added after the last of
  them (after the shared base before them when none is replaced).

  REF text. ``Aligned_protospacer+PAM_REF`` is the reference read from the start for as
  many bases as the target holds, in the target's orientation, with the target's gap
  columns put back at the same indices (letter case is not compared).

  Bulges. A gap in the aligned spacer is a DNA bulge, a gap in the aligned target an RNA
  bulge. None falls in the PAM columns; the spacer's guide end away from the PAM is never
  a gap (a leading RNA bulge, and a DNA bulge right next to the PAM, are allowed).

  Junction rule. A target holds its record when it holds a base the ALT adds or
  replaces; for a pure deletion, when it holds the bases on both sides of the deleted
  ones, side by side. A target that holds none lies in a flank: reference or SNP sequence,
  which the reference and SNP searches report; it is no indel row, and completeness does
  not ask for it.

  CARRIERS (soundness). The pool is the listed samples and every carrier of the row's
  record. A listed sample whose haplotypes do not hold the target at its start with every
  named variant, or that is listed for a weaker target than its strongest (PUTATIVE), fails
  unless the documented view lists it (then a limit, with its cause); a sample whose
  haplotype holds the target at the start with the named variants, or whose strongest
  PUTATIVE target it is, fails when the row omits it, unless the documented view leaves it
  off (a limit). The merge never moves a carrier off its sequence's row (same position,
  variants and aligned target merge carriers; equal carrier sets fold to the best
  alignment), so the row itself must list it.

  Merge (completeness). The pipeline clusters rows by guide and chromosome with Cluster
  Positions within --merge bp (3) of the previous, chained, whatever the strand. Inside a
  cluster, variant rows of one (Cluster Position, variant list) -- and for an indel/MNV
  row one aligned target -- merge their carriers, and of the rows of one (Cluster
  Position, variant list, carrier set) only the best is kept. CRISPRitz's Cluster Position
  is the start moved to the PAM's fixed end (less the RNA bulges and plus the DNA bulges on
  the + strand with a 3' PAM, the - strand with a 5' one), which every alignment of one
  PAM site shares. So a (sample, target) is satisfied by a row reporting exactly it with
  its variants, or by a row of the same Cluster Position naming the same variants whose
  carrier set is the target's own as CRISPRme writes it (every carrier the documented view
  lists for it, and nobody else but carriers it resolves PUTATIVE). In a run with one row
  per cluster, the cluster of a target is the integrated row with a block within --merge
  bp of its Cluster Position, else the nearest row on either side; its fewest-mm+b block
  must have no more edits than the target (the first sorting criterion, mm+bulges by
  default) and its highest-CFD block at least the target's CFD.

LIMITS (non-gating, counted per class, one TSV line per row block or per (record,
target)): limit-last-alt-wins, limit-two-indels, limit-overlapping-records,
limit-putative-carrier (all above), limit-contig-end (a record within one target of a
contig's end, or a row whose REF text runs past it), limit-n-bases (a stretch of N in a
record's window: no target is aligned across it), limit-putative-cap (a carrier with more
than 10 positions of unknown phase within a target's reach of a record: its targets that
touch one are not enumerated), and in the SNP-row report limit-snp-row-indel and
limit-snp-row-unphased.
Rows on a chromosome no VCF covers, or with no reference sequence, and table lines shorter
than their header are reported (gating): nothing is skipped silently.

OUTPUT. A TSV of every finding with its reason class, gating flag and evidence, a
per-class summary on stdout, exit status 1 on any gating finding (2 on a usage error).
``--known`` lists findings already reported upstream: they turn non-gating
(known-defect), and a listed one that no longer shows is reported (known-not-seen).
SNP-only rows of a per-sample run can be checked in report-only mode
(``--snp-rows report``); population-level SNP rows are a documented worst-case union and
are out of scope.

VCF ACCESS. A .tbi beside a VCF makes every read a seek. Without one, each VCF is read
once into an in-memory block index (nothing is written), which costs one decompression of
the file; a canonical install (one VCF per chromosome, 28 GB) is then feasible with
--jobs. ``--jobs N`` splits the work by region over N processes.
"""

import argparse
import bisect
import csv
import glob
import gzip
import io
import itertools
import json
import os
import re
import struct
import sys
import zipfile
import zlib
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
PUTATIVE_CAP = 10  # positions of unknown phase enumerated per carrier and record


def rc(s):
    return s.translate(COMP)[::-1]


def comp(b):
    return b.translate(COMP)


def iupac_match(motif_base, base):
    return base.upper() in IUPAC.get(motif_base.upper(), "")


# the target bases each spacer character (a guide base or an IUPAC PAM code, upper case)
# matches, for the aligner's inner loop
_MATCHES = {m: frozenset(v) for m, v in IUPAC.items()}


# --------------------------------------------------------------------------- #
# BGZF + tabix: read a VCF region without pysam (an index is optional)
# --------------------------------------------------------------------------- #
class Bgzf(object):
    """Just enough BGZF to start reading at a virtual offset. A BGZF file is a series of
    gzip members, each carrying its size in a BC extra field."""

    def __init__(self, path):
        self.path = path
        self.fh = open(path, "rb")
        head = self.fh.read(18)
        self.fh.seek(0)
        self.bgzf = (len(head) >= 18 and head[:2] == b"\x1f\x8b" and (head[3] & 4)
                     and head[12:14] == b"BC")

    def close(self):
        self.fh.close()

    def lines(self, voffset=0):
        """Decoded lines from a virtual offset (coffset << 16 | uoffset). A plain gzip
        file has no blocks to seek to, so it is streamed from its start."""
        if not self.bgzf:
            with gzip.open(self.path, "rt") as fh:
                for line in fh:
                    yield line.rstrip("\n")
            return
        self.fh.seek(voffset >> 16)
        skip = voffset & 0xFFFF
        tail = b""
        while True:
            block = self._member()
            if block is None:
                break
            if skip:
                block = block[skip:]
                skip = 0
            buf = tail + block
            parts = buf.split(b"\n")
            tail = parts.pop()
            for p in parts:
                yield p.decode("utf-8", "replace")
        if tail:
            yield tail.decode("utf-8", "replace")

    def blocks(self):
        """(compressed offset, decompressed bytes) of every block, in order."""
        self.fh.seek(0)
        while True:
            at = self.fh.tell()
            block = self._member()
            if block is None:
                return
            yield at, block

    def _member(self):
        head = self.fh.read(12)
        if len(head) < 12:
            return None
        if head[:2] != b"\x1f\x8b" or not (head[3] & 4):
            raise ValueError("not a BGZF member in %s" % self.path)
        xlen = struct.unpack("<H", head[10:12])[0]
        extra = self.fh.read(xlen)
        bsize = None
        i = 0
        while i + 4 <= xlen:
            slen = struct.unpack("<H", extra[i + 2:i + 4])[0]
            if extra[i] == 66 and extra[i + 1] == 67 and slen == 2:
                bsize = struct.unpack("<H", extra[i + 4:i + 6])[0] + 1
            i += 4 + slen
        if bsize is None:
            raise ValueError("a gzip member without BGZF's BC field in %s" % self.path)
        payload = self.fh.read(bsize - 12 - xlen)
        return zlib.decompressobj(-15).decompress(payload[:-8])


def _reg2bins(beg, end):
    """The tabix bins that can hold a feature overlapping [beg, end) (0-based)."""
    out = [0]
    end -= 1
    for sh, start in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        out.extend(range(start + (beg >> sh), start + (end >> sh) + 1))
    return out


class Tabix(object):
    """A tabix (.tbi) index, read here rather than through pysam."""

    def __init__(self, path):
        with gzip.open(path, "rb") as fh:
            buf = fh.read()
        if buf[:4] != b"TBI\x01":
            raise ValueError("not a tabix index: %s" % path)
        n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack("<8i", buf[4:36])
        names = buf[36:36 + l_nm].split(b"\x00")
        self.refs = {}
        off = 36 + l_nm
        for k in range(n_ref):
            n_bin = struct.unpack("<i", buf[off:off + 4])[0]
            off += 4
            bins = {}
            for _ in range(n_bin):
                b, n_chunk = struct.unpack("<Ii", buf[off:off + 8])
                off += 8
                chunks = []
                for _c in range(n_chunk):
                    cb, ce = struct.unpack("<QQ", buf[off:off + 16])
                    off += 16
                    chunks.append((cb, ce))
                bins[b] = chunks
            n_intv = struct.unpack("<i", buf[off:off + 4])[0]
            off += 4
            intv = struct.unpack("<%dQ" % n_intv, buf[off:off + 8 * n_intv]) if n_intv else ()
            off += 8 * n_intv
            self.refs[names[k].decode()] = (bins, intv)

    def chunks(self, chrom, beg, end):
        """Merged (begin, end) virtual-offset chunks for the region, 0-based [beg, end)."""
        if chrom not in self.refs:
            return []
        bins, intv = self.refs[chrom]
        cand = []
        for b in _reg2bins(beg, end):
            cand.extend(bins.get(b, ()))
        if not cand:
            return []
        lo = intv[min(beg >> 14, len(intv) - 1)] if intv else 0
        cand = sorted((cb, ce) for cb, ce in cand if ce > lo)
        out = []
        for cb, ce in cand:
            cb = max(cb, lo)
            if out and cb <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], ce))
            else:
                out.append((cb, ce))
        return out


class Vcf(object):
    """A bgzipped VCF read by region, with its genotypes. Our own reader."""

    LONGEST = 20000  # how far before a region an unindexed read starts (long deletions)

    def __init__(self, path):
        self.path = path
        self.samples = []
        self.contigs = set()
        self.tbi = None
        self.first_chrom = None
        self.blockidx = None  # {chrom: ([first POS of a block], [block offset])}
        if os.path.exists(path + ".tbi"):
            try:
                self.tbi = Tabix(path + ".tbi")
            except Exception:  # noqa: BLE001 - an unreadable index: read without one
                self.tbi = None
        bg = Bgzf(path)
        self.bgzf = bg.bgzf
        try:
            for line in bg.lines(0):
                if line.startswith("##"):
                    if line.startswith("##contig=<ID="):
                        self.contigs.add(line[13:].split(",")[0].split(">")[0])
                    continue
                if line.startswith("#CHROM"):
                    self.samples = line.rstrip("\n").split("\t")[9:]
                    continue
                if line:
                    self.first_chrom = line.split("\t", 1)[0]
                break
        finally:
            bg.close()

    def build_blockidx(self):
        """One pass over the file, keeping the CHROM and POS of the first whole line of
        each block: an index in memory (nothing is written next to a shared file)."""
        if self.blockidx is not None or self.tbi is not None or not self.bgzf:
            return
        idx = {}
        bg = Bgzf(self.path)
        try:
            for at, block in bg.blocks():
                nl = block.find(b"\n")
                while nl >= 0 and nl + 1 < len(block) and block[nl + 1:nl + 2] == b"#":
                    nl = block.find(b"\n", nl + 1)
                if nl < 0 or nl + 1 >= len(block):
                    continue
                end = block.find(b"\t", nl + 1)
                end2 = block.find(b"\t", end + 1) if end >= 0 else -1
                if end < 0 or end2 < 0:
                    continue
                chrom = block[nl + 1:end].decode()
                try:
                    pos = int(block[end + 1:end2])
                except ValueError:
                    continue
                ent = idx.setdefault(chrom, ([], []))
                ent[0].append(pos)
                ent[1].append(at)
                # a block that starts on one chromosome and holds the next: the next one's
                # first line is only found by scanning the block, which the reader does
                last = block.rfind(b"\n", 0, len(block) - 1)
                if last > nl:
                    e = block.find(b"\t", last + 1)
                    if e > 0:
                        lc = block[last + 1:e].decode()
                        if lc != chrom:
                            idx.setdefault(lc, ([], []))
                            idx[lc][0].append(0)
                            idx[lc][1].append(at)
        finally:
            bg.close()
        self.blockidx = idx

    def contig(self, chrom):
        """What this VCF calls ``chrom`` (chr22 or 22)."""
        names = [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]
        known = (set(self.tbi.refs) if self.tbi is not None
                 else set(self.blockidx) if self.blockidx is not None
                 else self.contigs | ({self.first_chrom} if self.first_chrom else set()))
        for n in names:
            if n in known:
                return n
        return chrom

    def chroms(self):
        """The chromosomes this VCF holds (as it names them)."""
        if self.tbi is not None:
            return list(self.tbi.refs)
        if self.bgzf:
            self.build_blockidx()
            return list(self.blockidx)
        return [self.first_chrom] if self.first_chrom else []

    def covers(self, chrom):
        """True / False when this VCF can tell; None for a plain-gzip file that has to be
        read to know."""
        name = self.contig(chrom)
        if self.tbi is not None:
            return name in self.tbi.refs
        if self.bgzf:
            self.build_blockidx()
            return name in self.blockidx
        return None if self.first_chrom != name else True

    def records(self, chrom, beg, end):
        """Every record whose REF span overlaps [beg, end) (0-based), as a dict."""
        chrom = self.contig(chrom)
        bg = Bgzf(self.path)
        try:
            if self.tbi is not None:
                for cb, ce in self.tbi.chunks(chrom, beg, end):
                    for line in bg.lines(cb):
                        if bg.fh.tell() > (ce >> 16) + (1 << 17):
                            break
                        r = self._parse(line, chrom, beg, end)
                        if r == "past":
                            break
                        if r:
                            yield r
                return
            start = 0
            if self.bgzf:
                self.build_blockidx()
                ent = self.blockidx.get(chrom)
                if not ent:
                    return
                i = bisect.bisect_right(ent[0], beg - self.LONGEST) - 1
                start = ent[1][max(i, 0)] << 16
            seen = False
            for line in bg.lines(start):
                if not line or line[0] == "#":
                    continue
                if not line.startswith(chrom + "\t"):
                    if seen:
                        break  # past the chromosome (a VCF is sorted)
                    continue
                seen = True
                r = self._parse(line, chrom, beg, end)
                if r == "past":
                    break
                if r:
                    yield r
        finally:
            bg.close()

    def _parse(self, line, chrom, beg, end):
        """A record as a dict; the genotype columns stay unsplit (only the lines a window
        keeps are split)."""
        if not line or line.startswith("#"):
            return None
        head, sep, tail = line.rstrip("\n").partition("\t")
        if head != chrom or not sep:
            return None
        f = tail.split("\t", 8)
        if len(f) < 8:
            return None
        pos1 = int(f[0])
        if pos1 - 1 >= end:
            return "past"
        ref = f[2].upper()
        if pos1 - 1 + len(ref) <= beg:
            return None
        return {"chrom": head, "pos1": pos1, "id": f[1], "ref": ref,
                "alts": [a.upper() for a in f[3].split(",")], "info": f[6],
                "format": f[7] if len(f) > 7 else "", "tail": f[8] if len(f) > 8 else ""}


# --------------------------------------------------------------------------- #
# reference FASTA, read without writing an .fai
# --------------------------------------------------------------------------- #
class Reference(object):
    """A per-chromosome FASTA folder or one FASTA of any number of contigs. Only the
    named contig is read (chr22 or 22), and nothing is written beside it."""

    def __init__(self, path):
        self.path = path
        self.cache = {}

    def _names(self, chrom):
        return [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]

    def _file(self, chrom):
        if os.path.isdir(self.path):
            for n in self._names(chrom):
                for name in ("%s.fa" % n, "%s.fasta" % n, "%s.fa.gz" % n, "%s.fasta.gz" % n):
                    p = os.path.join(self.path, name)
                    if os.path.exists(p):
                        return p
            raise FileNotFoundError("no FASTA for %s under %s" % (chrom, self.path))
        return self.path

    def has(self, chrom):
        """Whether the reference holds ``chrom`` (a file of the folder, or a header of
        the one FASTA, read once)."""
        if os.path.isdir(self.path):
            try:
                self._file(chrom)
                return True
            except FileNotFoundError:
                return False
        if not hasattr(self, "_names_in"):
            op = gzip.open if self.path.endswith(".gz") else open
            with op(self.path, "rt") as fh:
                self._names_in = {line[1:].split()[0] for line in fh
                                  if line.startswith(">") and line[1:].split()}
        return bool(set(self._names(chrom)) & self._names_in)

    def seq(self, chrom):
        if chrom in self.cache:
            return self.cache[chrom]
        p = self._file(chrom)
        names = set(self._names(chrom))
        op = gzip.open if p.endswith(".gz") else open
        parts = []
        found = False
        with op(p, "rt") as fh:
            want = False
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0] if line[1:].split() else ""
                    if want:
                        break  # the contig is read whole
                    want = name in names
                    found = found or want
                    continue
                if want:
                    parts.append(line.strip())
        if not found:
            raise KeyError("no contig %s in %s" % (chrom, p))
        s = "".join(parts).upper()
        if not getattr(self, "keep_all", False):
            self.cache = {}
        self.cache[chrom] = s
        return s


# --------------------------------------------------------------------------- #
# ALT alleles and genotypes, as VCF 4.x defines them
# --------------------------------------------------------------------------- #
class Allele(object):
    """One ALT allele of a VCF line. On a haplotype that carries it, ALT stands in place
    of the REF bases [lo, hi) (0-based). The bases REF and ALT share at either end stay
    where they are, so the allele's edit is what lies between them: it replaces the
    reference bases [edit_lo, edit_hi) with ``body``."""

    __slots__ = ("key", "lo", "hi", "line", "number", "edit_lo", "edit_hi", "body",
                 "kind", "_bases")

    def __init__(self, pos1, ref, alt, line, number):
        self.key = (pos1, ref, alt)
        self.lo = pos1 - 1
        self.hi = self.lo + len(ref)
        self.line = line        # the line's id in its Window
        self.number = number    # the allele's number in the GT field (1 = the first ALT)
        head = len(os.path.commonprefix([ref, alt]))
        back = len(os.path.commonprefix([ref[head:][::-1], alt[head:][::-1]]))
        self.edit_lo = self.lo + head
        self.edit_hi = self.hi - back
        self.body = alt[head:len(alt) - back]
        if len(ref) == 1 and len(alt) == 1:
            self.kind = "snv"
        elif len(ref) == len(alt):
            self.kind = "mnv"
        else:
            self.kind = "indel"
        self._bases = None

    def token(self, chrom):
        return "%s_%d_%s_%s" % (chrom, self.key[0], self.key[1], self.key[2])

    def bases(self):
        """[(base, (coordinate, k))] of the edit's ALT bases: the j-th stands at reference
        coordinate edit_lo + j while a replaced base is left to stand in (k = 0); the rest
        are the 1st, 2nd... base added after the last replaced base, or after the shared
        base before the edit when it replaces none (an insertion's anchor)."""
        return list(zip(self.body, self.places()))

    def places(self):
        """The places of the edit's ALT bases (see bases())."""
        if self._bases is None:
            width = self.edit_hi - self.edit_lo
            stand = min(width, len(self.body))
            out = [(self.edit_lo + j, 0) for j in range(stand)]
            anchor = self.edit_lo + width - 1
            out.extend((anchor, k) for k in range(1, len(self.body) - stand + 1))
            self._bases = out
        return self._bases

    def window(self, lo, hi):
        """(first, last) indexes of the edit's ALT bases whose place lies in [lo, hi):
        the bases standing there, and the added ones when their anchor does."""
        width = self.edit_hi - self.edit_lo
        stand = min(width, len(self.body))
        a = max(0, lo - self.edit_lo)
        b = max(a, min(stand, hi - self.edit_lo))
        if len(self.body) > stand and lo <= self.edit_lo + width - 1 < hi:
            b = len(self.body)  # the anchor is in the window, and so are its added bases
        return a, b

    def adds(self):
        return len(self.body) > self.edit_hi - self.edit_lo


def touches(al, a, b):
    """Whether an allele changes what a haplotype holds between reference coordinates a
    and b (inclusive): it replaces or deletes one of them, or adds bases between two."""
    if al.edit_hi > al.edit_lo and al.edit_lo <= b and al.edit_hi - 1 >= a:
        return True
    if al.adds():
        anchor = al.edit_hi - 1 if al.edit_hi > al.edit_lo else al.edit_lo - 1
        return a <= anchor < b
    return False


def changes_target(al, t):
    """Whether an allele changes target ``t``: the target holds a base it wrote (the
    target was read on a haplotype that carries it -- one ending inside an insertion
    anchored at its last coordinate included), or writing it would change the reference
    stretch the target covers (touches())."""
    return any(x is al for x in t.src) or touches(al, t.span[0], t.span[1])


def parse_gt(gt):
    """(allele numbers, phased) of a GT value: None for a missing allele; '|' phased,
    '/' unphased, a single allele (haploid) phased."""
    if "/" in gt:
        parts, phased = gt.replace("|", "/").split("/"), False
    elif "|" in gt:
        parts, phased = gt.split("|"), True
    else:
        parts, phased = [gt], True
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            out.append(None)
    return tuple(out), phased


_HOMREF = frozenset(("0|0", "0/0", "0", "0|0|0", "0/0/0"))


def parse_calls(fmt, tail, names):
    """{sample: (allele numbers, phased)} of every sample whose call is not all reference
    (GT must be the first FORMAT key, as VCF 4.x requires when it is present)."""
    out = {}
    if not tail or not fmt.startswith("GT"):
        return out
    fields = tail.split("\t")
    if fmt != "GT":
        fields = [g.partition(":")[0] for g in fields]
    homref = _HOMREF
    for i in [i for i, g in enumerate(fields[:len(names)]) if g not in homref]:
        alle, phased = parse_gt(fields[i])
        if any(x != 0 for x in alle):
            out[names[i]] = (alle, phased)
    return out


class Window(object):
    """The VCF lines of one region: their ALT alleles and the samples' calls.

    ``shadowed`` and ``untied`` belong to the documented view only: the SNV alleles
    CRISPRme's enriched genome and SNP stores never see (an earlier line of a split
    multi-allelic site: only the last is kept, CHANGELOG 2.8.1 Known limits), and the
    alts of a multi-ALT SNV line whose allele number the genotypes cannot tell CRISPRme
    (METHODS section 4: their carriers are resolved as unphased ones)."""

    def __init__(self, refseq, chrom, beg, end, vcfs, needed=None):
        self.chrom = chrom
        self.refseq = refseq
        self.beg = max(0, beg)
        self.end = end
        self.needed = needed
        self.samples = []
        self.alleles = []
        self.by_key = defaultdict(list)
        self.line_lo = []       # line id -> 0-based start of its REF
        self.line_hi = []
        self.line_alleles = []  # line id -> [Allele or None per ALT]
        self.calls = []         # line id -> {sample: (allele numbers, phased)}
        self.shadowed = set()
        self.untied = set()
        self.snv_at = defaultdict(list)  # reference coordinate -> SNV alleles there
        self._order = None
        seen = set()
        for v in vcfs:
            for s in v.samples:
                if s not in seen:
                    seen.add(s)
                    self.samples.append(s)
            names = v.samples
            snv_lines = defaultdict(list)  # POS -> line ids of SNV-only lines, file order
            for r in v.records(chrom, self.beg, self.end):
                lo = r["pos1"] - 1
                hi = lo + len(r["ref"])
                if not self._wanted(lo, hi):
                    continue
                lid = len(self.line_lo)
                als = []
                for k, alt in enumerate(r["alts"]):
                    if (not alt or alt in (".", "*") or alt.startswith("<")
                            or set(alt) - set("ACGTN")):
                        als.append(None)
                        continue
                    als.append(Allele(r["pos1"], r["ref"], alt, lid, k + 1))
                if not any(a is not None for a in als):
                    continue
                self.line_lo.append(lo)
                self.line_hi.append(hi)
                self.line_alleles.append(als)
                self.calls.append(parse_calls(r["format"], r["tail"], names))
                for a in als:
                    if a is None:
                        continue
                    self.alleles.append(a)
                    self.by_key[a.key].append(a)
                    if a.kind == "snv":
                        self.snv_at[a.lo].append(a)
                snvs = [a for a in als if a is not None]
                if len(r["ref"]) == 1 and all(a.kind == "snv" for a in snvs):
                    snv_lines[r["pos1"]].append(lid)
                    if len(snvs) > 1:
                        self._untie(lid)
            for lids in snv_lines.values():
                if len(lids) < 2:
                    continue
                last = {a.key for a in self.line_alleles[lids[-1]] if a is not None}
                for lid in lids[:-1]:
                    for a in self.line_alleles[lid]:
                        if a is not None and a.key not in last:
                            self.shadowed.add(a)

    def _wanted(self, lo, hi):
        if self.needed is None:
            return True
        i = bisect.bisect_right(self.needed, (lo, float("inf")))
        for j in (i - 1, i):
            if 0 <= j < len(self.needed):
                a, b = self.needed[j]
                if lo < b and hi > a:
                    return True
        return False

    def _untie(self, lid):
        """Which alts of a multi-ALT SNV line CRISPRme can read against their own allele
        number (documented view): those some carrier holds alone on the line, and the last
        one left once every other alt of the line is told."""
        als = [a for a in self.line_alleles[lid] if a is not None]
        calls = self.calls[lid]
        told = set()
        for a in als:
            for alle, _ph in calls.values():
                if {x for x in alle if x} == {a.number}:
                    told.add(a)
                    break
        rest = [a for a in als if a not in told]
        if len(rest) == 1:
            told.add(rest[0])
        for a in als:
            if a not in told:
                self.untied.add(a)

    def lines_in(self, lo, hi):
        """The line ids whose REF span meets [lo, hi), in position order."""
        if self._order is None:
            self._order = sorted(range(len(self.line_lo)), key=lambda i: self.line_lo[i])
            self._starts = [self.line_lo[i] for i in self._order]
            self._widest = max([h - l for l, h in zip(self.line_lo, self.line_hi)] or [1])
        i = bisect.bisect_left(self._starts, lo - self._widest)
        j = bisect.bisect_left(self._starts, hi)
        return [k for k in self._order[i:j] if self.line_hi[k] > lo and self.line_lo[k] < hi]

    def call(self, allele, sample):
        return self.calls[allele.line].get(sample)

    def carriers(self, allele):
        """The samples whose call holds the allele on some slot."""
        return [s for s, (alle, _p) in self.calls[allele.line].items()
                if allele.number in alle]

    def signature(self, sample, lids):
        return tuple((lid, self.calls[lid][sample]) for lid in lids
                     if sample in self.calls[lid])


# --------------------------------------------------------------------------- #
# haplotypes: the reference with a slot's alleles in place
# --------------------------------------------------------------------------- #
class Hap(object):
    """A stretch of one haplotype: its bases, each base's place -- (reference coordinate,
    0) for a base standing there, (coordinate, k) for the k-th base an allele adds after
    it -- the allele that wrote it (None: the reference), and the spans where the alleles
    given overlap (no sequence: ``clash``)."""

    __slots__ = ("seq", "place", "src", "clash")

    def __init__(self, seq, place, src, clash):
        self.seq = seq
        self.place = place
        self.src = src
        self.clash = clash

    def span_at(self, coord):
        """(first, end) indexes of the bases whose place has reference coordinate
        ``coord``: the places run in order along a haplotype, so a bisection finds them."""
        i = bisect.bisect_left(self.place, (coord, 0))
        j = bisect.bisect_left(self.place, (coord + 1, 0), i)
        return i, j

    def at(self, coord):
        """Indexes of the bases whose place has reference coordinate ``coord``."""
        i, j = self.span_at(coord)
        return range(i, j)

    def index(self, place):
        i = bisect.bisect_left(self.place, place)
        return i if i < len(self.place) and self.place[i] == place else None

    def clashes(self, a, b):
        """The overlap zones (lo, hi, alleles) that meet reference coordinates a..b."""
        return [z for z in self.clash if z[0] <= b and z[1] - 1 >= a]


def make_hap(refseq, alleles, lo, hi):
    """[lo, hi) of the haplotype that holds ``alleles``. Alleles whose REF spans meet
    describe one stretch twice: none of them is written and the stretch is in ``clash``.
    Only the window is written: an allele reaching into it from either side puts in the
    bases of its edit that stand (or are added after a base) inside it."""
    hi = min(hi, len(refseq))
    use = sorted((a for a in alleles if a.hi > lo and a.lo < hi), key=lambda a: (a.lo, a.hi))
    runs = []
    for a in use:
        if runs and a.lo < runs[-1][1]:
            runs[-1][0].append(a)
            runs[-1][1] = max(runs[-1][1], a.hi)
        else:
            runs.append([[a], a.hi])
    seq, place, src, clash = [], [], [], []
    p = lo

    def ref_upto(stop):
        stop = min(stop, hi)
        if stop > p:
            seq.append(refseq[p:stop])
            place.extend([(q, 0) for q in range(p, stop)])
            src.extend([None] * (stop - p))
        return max(p, stop)

    for members, end in runs:
        if len(members) > 1:
            clash.append((members[0].lo, end, tuple(members)))
            continue
        a = members[0]
        p = ref_upto(a.edit_lo)
        i, j = a.window(lo, hi)
        if j > i:
            seq.append(a.body[i:j])
            place.extend(a.places()[i:j])
            src.extend([a] * (j - i))
        p = max(p, a.edit_hi)
    ref_upto(hi)
    return Hap("".join(seq), place, src, clash)


def holds_edit(hap, i, n, al):
    """The junction rule on hap[i:i + n]: it holds a base the allele adds or replaces, or,
    for a pure deletion, the bases either side of the deleted ones, side by side."""
    if al.body:
        return any(hap.src[k] is al for k in range(i, i + n))
    before, after = (al.edit_lo - 1, 0), (al.edit_hi, 0)
    for k in range(i, i + n - 1):
        if hap.place[k] == before and hap.place[k + 1] == after:
            return True
    return False


def reads_at(hap, coord, probe, al=None):
    """The indexes where a base of reference coordinate ``coord`` starts ``probe`` (and,
    with ``al``, the stretch holds that allele by the junction rule)."""
    n = len(probe)
    i, j = hap.span_at(coord)
    out = []
    k = hap.seq.find(probe, i, j - 1 + n)
    while 0 <= k < j:
        if al is None or holds_edit(hap, k, n, al):
            out.append(k)
        k = hap.seq.find(probe, k + 1, j - 1 + n)
    return out


def span_of(places):
    cs = [pl[0] for pl in places]
    return min(cs), max(cs)


# --------------------------------------------------------------------------- #
# CFD, from the published tables CRISPRme ships as data
# --------------------------------------------------------------------------- #
class Cfd(object):
    """The CFD score (Doench et al. 2016) as METHODS section 8 documents it: one factor
    per spacer column at its guide position (a DNA bulge's extra base takes no guide
    position), a bulge's gap factor at its position, the two undefined terminal columns
    extrapolated -- lenient at the PAM-distal 5' tip (0.9, the value CRISPRme's scorer
    uses), abolished at the PAM junction -- times the factor of the PAM's last two bases.
    Defined only for a 20-nt guide with a 3-nt 3' PAM, as CRISPRme computes it."""

    TIP, JUNCTION = 0.9, 0.0

    def __init__(self, folder):
        with open(os.path.join(folder, "cfd_mismatch_scores.json")) as fh:
            self.mm = json.load(fh)
        with open(os.path.join(folder, "cfd_pam_scores.json")) as fh:
            self.pam = json.load(fh)

    @staticmethod
    def applies(guide_len, pam_len, pam_5prime):
        return guide_len == 20 and pam_len == 3 and not pam_5prime

    def factor(self, g, d, pos):
        """One spacer column: guide base (or '-') g against target base (or '-') d."""
        g = g.upper().replace("T", "U")
        d = d.upper()
        if g == "-":
            if pos <= 1:
                return self.TIP
            if pos > 20:
                return self.JUNCTION
            return self.mm.get("r-:d%s,%d" % (comp(d), pos), self.JUNCTION)
        if d == "-":
            if pos <= 1:
                return self.TIP
            if pos > 20:
                return self.JUNCTION
            return self.mm.get("r%s:d-,%d" % (g, pos), 1.0)
        if d not in "ACGT":
            return 0.0
        if g == d.replace("T", "U"):
            return 1.0
        return self.mm.get("r%s:d%s,%d" % (g, comp(d), pos), 1.0)

    @staticmethod
    def positions(crrna, n_spacer):
        """The guide position of each spacer column."""
        out = []
        dna = 0
        for c in range(n_spacer):
            out.append(c + 1 - dna)
            if crrna[c] == "-":
                dna += 1
        return out

    def score(self, crrna, target, pam_len=3):
        n = len(crrna) - pam_len
        pos = self.positions(crrna, n)
        s = 1.0
        for c in range(n):
            s *= self.factor(crrna[c], target[c], pos[c])
        return s * self.pam.get(target[-2:].upper(), 0.0)


# --------------------------------------------------------------------------- #
# our own exhaustive aligner
# --------------------------------------------------------------------------- #
_ALIGN_CACHE = {}


def align_pairs(spacer, target, max_mm, max_dna, max_rna, max_total, pam_len,
                pam_5prime):
    """Every alignment of ``spacer`` (guide + PAM motif, IUPAC) to ``target`` (plain
    bases) within the budget. A gap in the spacer is a DNA bulge, a gap in the target an
    RNA bulge; the PAM holds no gap, the spacer's first and last columns never one.
    Yields (aligned spacer, aligned target, mismatches, dna bulges, rna bulges)."""
    ck = (spacer, target, max_mm, max_dna, max_rna, max_total, pam_len, pam_5prime)
    hit = _ALIGN_CACHE.get(ck)
    if hit is not None:
        return hit
    n, m = len(spacer), len(target)
    target = target.upper()
    matches = [_MATCHES.get(c.upper(), frozenset()) for c in spacer]
    if pam_5prime:
        pam_first, pam_last = 0, pam_len - 1
    else:
        pam_first, pam_last = n - pam_len, n - 1
    out = []
    if not _feasible(matches, target, n, m, max_mm, max_dna, max_rna, max_total,
                     pam_first, pam_last, pam_5prime):
        _ALIGN_CACHE[ck] = out
        return out

    def walk(i, j, mm, dn, rn, sa, ta):
        if mm > max_mm or dn > max_dna or rn > max_rna or mm + dn + rn > max_total:
            return
        if n - i > m - j + max_rna or m - j > n - i + max_dna:
            return
        if i == n:
            if j == m and sa and sa[0] != "-" and sa[-1] != "-":
                out.append(("".join(sa), "".join(ta), mm, dn, rn))
            return
        in_pam = pam_first <= i <= pam_last
        if j < m:
            ok = target[j] in matches[i]
            if ok or not in_pam:
                sa.append(spacer[i])
                ta.append(target[j])
                walk(i + 1, j + 1, mm + (0 if ok else 1), dn, rn, sa, ta)
                sa.pop()
                ta.pop()
        dna_ok = (i <= pam_first) if not pam_5prime else (i > pam_last)
        if dna_ok and dn < max_dna and j < m and len(sa) > 0:
            sa.append("-")
            ta.append(target[j])
            walk(i, j + 1, mm, dn + 1, rn, sa, ta)
            sa.pop()
            ta.pop()
        if not in_pam and rn < max_rna:
            sa.append(spacer[i])
            ta.append("-")
            walk(i + 1, j, mm, dn, rn + 1, sa, ta)
            sa.pop()
            ta.pop()

    walk(0, 0, 0, 0, 0, [], [])
    if len(_ALIGN_CACHE) > 200000:
        _ALIGN_CACHE.clear()
    _ALIGN_CACHE[ck] = out
    return out


def _feasible(matches, target, n, m, max_mm, max_dna, max_rna, max_total, pam_first,
              pam_last, pam_5prime):
    """Whether any alignment of the spacer to ``target`` fits the budget: the fewest
    mismatches for each count of DNA and RNA bulges, by dynamic programming over the
    spacer's positions (a DNA bulge takes a target base, an RNA bulge a spacer base, so the
    target position follows). It only rules out pieces with no in-budget alignment, so the
    exhaustive walk that follows finds exactly what it found without it."""
    inf = max_mm + 1
    # cost[(dn, rn)] at spacer position i: fewest mismatches so far
    cost = {(0, 0): 0}
    for i in range(n + 1):
        # DNA bulges at this spacer position (target bases with no spacer base)
        dna_ok = (i <= pam_first) if not pam_5prime else (i > pam_last)
        if dna_ok and i > 0 and max_dna:
            for dn in range(max_dna):
                for rn in range(max_rna + 1):
                    c = cost.get((dn, rn))
                    j = i + dn - rn
                    if c is None or j >= m:
                        continue
                    if c < cost.get((dn + 1, rn), inf):
                        cost[(dn + 1, rn)] = c
        if i == n:
            break
        nxt = {}
        in_pam = pam_first <= i <= pam_last
        for (dn, rn), c in cost.items():
            j = i + dn - rn
            if j < m:
                ok = target[j] in matches[i]
                if ok or not in_pam:
                    c2 = c + (0 if ok else 1)
                    if c2 <= max_mm and c2 < nxt.get((dn, rn), inf):
                        nxt[(dn, rn)] = c2
            if not in_pam and rn < max_rna and c < nxt.get((dn, rn + 1), inf):
                nxt[(dn, rn + 1)] = c
        cost = nxt
        if not cost:
            return False
    for (dn, rn), c in cost.items():
        if i + dn - rn == m and c <= max_mm and c + dn + rn <= max_total:
            return True
    return False


def pam_sites(motif, seq, pam_len):
    out = []
    for q in range(0, len(seq) - pam_len + 1):
        if all(iupac_match(m, b) for m, b in zip(motif, seq[q:q + pam_len])):
            out.append(q)
    return out


def scan_sequence(spacer, seq, max_mm, max_dna, max_rna, max_total, pam_len, pam_5prime,
                  skipped=None):
    """Every in-budget alignment on ``seq``, both strands, anchored on the PAM. Yields
    (strand, offset of the leftmost aligned base on the + reading, aligned spacer,
    aligned target, mm, dna bulges, rna bulges). A piece holding an N is not aligned (no
    base pairs with N); ``skipped`` (a list) gets one entry per such piece."""
    n = len(spacer)
    motif = spacer[:pam_len] if pam_5prime else spacer[n - pam_len:]
    proto = n - pam_len
    for strand, s in (("+", seq), ("-", rc(seq))):
        for q in pam_sites(motif, s, pam_len):
            for width in range(proto - max_rna, proto + max_dna + 1):
                beg = (q - width) if not pam_5prime else q
                end = (q + pam_len) if not pam_5prime else (q + pam_len + width)
                if beg < 0 or end > len(s):
                    continue
                piece = s[beg:end]
                if "N" in piece:
                    if skipped is not None:
                        skipped.append((strand, beg))
                    continue
                for sa, ta, mm, dn, rn in align_pairs(
                        spacer, piece, max_mm, max_dna, max_rna, max_total, pam_len,
                        pam_5prime):
                    off = beg if strand == "+" else len(s) - end
                    yield strand, off, sa, ta, mm, dn, rn


# --------------------------------------------------------------------------- #
# results tables
# --------------------------------------------------------------------------- #
def open_table(path):
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.endswith(".tsv")]
            if not names:
                return io.StringIO("")
            return io.StringIO(z.read(names[0]).decode())
    return open(path, newline="")


def parse_token(tok):
    p = tok.rsplit("_", 3)
    return p[0], int(p[1]), p[2].upper(), p[3].upper()


def is_indel_token(tok):
    p = tok.rsplit("_", 3)
    return len(p) == 4 and (len(p[2]) > 1 or len(p[3]) > 1)


def cluster_of(start, strand, probe_len, full_len, pam_5prime):
    """CRISPRitz's Cluster Position of an alignment: its start, moved to the PAM's fixed
    end -- less the RNA bulges and plus the DNA bulges on the + strand with a 3' PAM (and
    the - strand with a 5' PAM), the start itself otherwise. Every alignment of one PAM
    site shares it, which is what the merge folds on."""
    if (strand == "+") != pam_5prime:
        return start + probe_len - full_len
    return start


def read_rows(path, table, chrom=None, lo=None, hi=None, margin=0, pam_5prime=False,
              loci=None, bad=None):
    """Every block (dict) of a results table that names a variant. Rows whose start is
    within [lo - margin, hi + margin) are read; those within [lo, hi) are marked for the
    soundness check ('own'). ``loci`` (a dict) gets, per chromosome, every integrated
    line's (Cluster Position of each block, fewest-mm+b total, highest-CFD score) for the
    per-cluster completeness check; ``bad`` (a list) gets the lines shorter than the
    header."""
    csv.field_size_limit(1 << 30)
    out = []
    with open_table(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        try:
            head = next(r)
        except StopIteration:
            return out
        col = {n: i for i, n in enumerate(head)}
        if "Chromosome" not in col:
            return out
        for line_no, row in enumerate(r, start=2):
            if not row:
                continue
            if len(row) < len(head):
                if bad is not None:
                    bad.append((table, line_no, len(row), len(head)))
                continue
            ch = row[col["Chromosome"]]
            if chrom and ch != chrom:
                continue
            full_len = len(row[col["Spacer+PAM"]]) if "Spacer+PAM" in col else 23
            line_blocks = []
            for b in BLOCKS:
                try:
                    start = int(row[col["Start_coordinate_(%s)" % b]])
                except (KeyError, ValueError):
                    continue
                if lo is not None and not (lo - margin <= start < hi + margin):
                    continue
                strand = row[col["Strand_(%s)" % b]]
                alt = row[col["Aligned_protospacer+PAM_ALT_(%s)" % b]]
                ref = row[col["Aligned_protospacer+PAM_REF_(%s)" % b]]
                shown = alt if alt not in ("", "NA", "n") else ref
                width = len(shown.replace("-", ""))
                line_blocks.append((b, cluster_of(start, strand, width, full_len,
                                                  pam_5prime)))
                origin = row[col["REF/ALT_origin_(%s)" % b]]
                var = row[col["Variant_info_genome_(%s)" % b]]
                if var in ("", "NA", "n"):
                    continue
                out.append({
                    "table": table, "line": line_no, "block": b, "chrom": ch,
                    "spacer_col": row[col["Spacer+PAM"]],
                    "start": start,
                    "own": lo is None or lo <= start < hi,
                    "strand": strand,
                    "crrna": row[col["Aligned_spacer+PAM_(%s)" % b]],
                    "ref": ref,
                    "alt": alt,
                    "pam": row[col["PAM_(%s)" % b]],
                    "mm": row[col["Mismatches_(%s)" % b]],
                    "bulges": row[col["Bulges_(%s)" % b]],
                    "total": row[col["Mismatches+bulges_(%s)" % b]],
                    "btype": row[col["Bulge_type_(%s)" % b]],
                    "origin": origin,
                    "var": var,
                    "samples": row[col["Variant_samples_(%s)" % b]],
                })
            if loci is not None and line_blocks:
                try:
                    total = int(row[col["Mismatches+bulges_(fewest_mm+b)"]])
                except (KeyError, ValueError):
                    total = None
                try:
                    cfd = float(row[col["CFD_score_(highest_CFD)"]])
                except (KeyError, ValueError):
                    cfd = None
                loci.setdefault(ch, []).append(([c for _b, c in line_blocks], total, cfd,
                                                line_no))
    return out


def result_tables(resdir):
    """(path, label) of the integrated results and the alternative alignments."""
    out = []
    for p in sorted(glob.glob(os.path.join(resdir, "*integrated_results.tsv"))):
        out.append((p, "integrated"))
    if not out:
        for p in sorted(glob.glob(os.path.join(resdir, "*integrated_results.zip"))):
            out.append((p, "integrated"))
    stem = os.path.join(resdir, "*all_results_with_alternative_alignments")
    alt = sorted(glob.glob(stem + ".tsv")) or sorted(glob.glob(stem + ".zip"))
    for p in alt:
        out.append((p, "alternative"))
    return out


def run_kind(resdir):
    """(analysis mode, alternative alignments written, merge window, first sorting
    criterion of the fewest block) of a run, from its own markers. A population-level
    run writes no alternative alignments by design, nor does a run made with --ci-cd-test
    (the documented chr22 self-test): both keep one row per cluster (.alt_mode =
    no-alt-alignments)."""
    mode = alt = None
    merge, criterion = 3, "mm+bulges"
    p = os.path.join(resdir, ".alt_mode")
    if os.path.exists(p):
        with open(p) as fh:
            v = fh.read().strip()
        alt = {"no-alt-alignments": False, "alt-alignments": True}.get(v)
    p = os.path.join(resdir, ".command_line.txt")
    if os.path.exists(p):
        with open(p) as fh:
            v = fh.read()
        if "--population-level" in v:
            mode = "population-level"
        elif "--per-sample" in v or "--ci-cd-test" in v:
            mode = "per-sample"
        m = re.search(r"--merge\s+(\d+)", v)
        if m:
            merge = int(m.group(1))
        m = re.search(r"--sorting-criteria\s+(\S+)", v)
        if m:
            criterion = m.group(1).split(",")[0]
    if mode == "population-level" and alt is None:
        alt = False
    return mode, alt, merge, criterion


def samples_of(row):
    return [s for s in row["samples"].split(",") if s and s not in ("NA", "n", ".")]


# --------------------------------------------------------------------------- #
# findings
# --------------------------------------------------------------------------- #
COLS = ["gating", "class", "table", "line", "block", "chrom", "start", "strand",
        "spacer", "target", "variants", "samples", "detail"]


class Report(object):
    def __init__(self):
        self.rows = []

    def add(self, gating, cls, row, detail, samples=None):
        self.rows.append({
            "gating": "yes" if gating else "no", "class": cls,
            "table": row.get("table", ""), "line": row.get("line", ""),
            "block": row.get("block", ""), "chrom": row.get("chrom", ""),
            "start": row.get("start", ""), "strand": row.get("strand", ""),
            "spacer": row.get("crrna", ""), "target": row.get("alt", ""),
            "variants": row.get("var", ""),
            "samples": samples if samples is not None else row.get("samples", ""),
            "detail": detail,
        })

    def gating(self):
        return sum(1 for r in self.rows if r["gating"] == "yes")

    def counts(self):
        out = defaultdict(int)
        for r in self.rows:
            out[r["class"]] += 1
        return out

    def write(self, path):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            for r in self.rows:
                w.writerow(r)


def recount(crrna, target, pam_len, pam_5prime):
    """(mismatches, dna bulges, rna bulges, problem) of an aligned pair as CRISPRme
    writes it: the spacer's gaps are DNA bulges, the target's RNA bulges, case ignored,
    no gap in the PAM, and the spacer's guide end away from the PAM never a gap."""
    if len(crrna) != len(target):
        return 0, 0, 0, "aligned spacer and target differ in length"
    n = len(crrna)
    pam_cols = set(range(0, pam_len)) if pam_5prime else set(range(n - pam_len, n))
    far = n - 1 if pam_5prime else 0
    if crrna[far] == "-":
        return 0, 0, 0, "a gap at the spacer's guide end"
    mm = dn = rn = 0
    for i, (g, d) in enumerate(zip(crrna, target)):
        if g == "-" and d == "-":
            return 0, 0, 0, "a column with a gap in both strings"
        if i in pam_cols and (g == "-" or d == "-"):
            return 0, 0, 0, "a gap in the PAM"
        if g == "-":
            dn += 1
            continue
        if d == "-":
            rn += 1
            continue
        if i in pam_cols:
            continue
        if not iupac_match(g, d):
            mm += 1
    return mm, dn, rn, None


def base_columns(alt_aligned, strand):
    """For each base of the target in + order, its column in the aligned strings."""
    cols = [j for j, ch in enumerate(alt_aligned) if ch != "-"]
    return cols if strand == "+" else cols[::-1]


class Ctx(object):
    """What every check needs: the budget, the PAM, the guides, the CFD tables."""

    def __init__(self, args, pam_motif, pam_len, pam_5prime, guides, cfd=None):
        self.args = args
        self.mm = args.mm
        self.bdna = args.bdna
        self.brna = args.brna
        self.cap = args.max_total_edits if args.max_total_edits is not None else 99
        self.pam_motif = pam_motif
        self.pam_len = pam_len
        self.pam_5prime = pam_5prime
        self.guides = guides  # [guide without PAM]
        self.cfd = cfd if cfd is not None and guides and Cfd.applies(
            len(guides[0]), pam_len, pam_5prime) else None

    def full(self, guide):
        return (self.pam_motif + guide) if self.pam_5prime else (guide + self.pam_motif)

    def pam_cols(self, n):
        return range(0, self.pam_len) if self.pam_5prime else range(n - self.pam_len, n)


# --------------------------------------------------------------------------- #
# a row's own strings
# --------------------------------------------------------------------------- #
def alignment_checks(row, ctx, rep, gating):
    """The checks of a row's own strings: its guide, its aligned spacer, the recount, the
    bulge type, the budget, the PAM. False when the row cannot be read further."""
    crrna, alt = row["crrna"], row["alt"]
    stripped = alt.replace("-", "").upper()
    if set(stripped) - set("ACGT"):
        rep.add(gating, "iupac-left", row, "the target holds %s"
                % "".join(sorted(set(stripped) - set("ACGT"))))
        return False
    pam_len, p5 = ctx.pam_len, ctx.pam_5prime
    spc = row["spacer_col"].upper()
    sguide = spc[pam_len:] if p5 else spc[:len(spc) - pam_len]
    guide = next((g for g in ctx.guides if g == sguide), None)
    if guide is None:
        rep.add(gating, "guide-unknown", row,
                "the row's Spacer+PAM %s is none of the given guides" % row["spacer_col"])
        return False
    nogap = crrna.replace("-", "").upper()
    g_part = nogap[pam_len:] if p5 else nogap[:len(nogap) - pam_len]
    p_part = nogap[:pam_len] if p5 else nogap[len(nogap) - pam_len:]
    if g_part != guide or not (p_part == "N" * pam_len or p_part == ctx.pam_motif):
        rep.add(gating, "spacer-not-guide", row,
                "the aligned spacer, gaps removed, is %s, not the guide %s with its PAM"
                % (nogap, guide))
        return False
    mm, dn, rn, problem = recount(crrna, alt, pam_len, p5)
    if problem:
        rep.add(gating, "alignment-invalid", row, problem)
        return False
    try:
        r_mm, r_bul, r_tot = int(row["mm"]), int(row["bulges"]), int(row["total"])
    except ValueError:
        rep.add(gating, "alignment-invalid", row, "unreadable counts")
        return False
    if (mm, dn + rn) != (r_mm, r_bul) or r_tot != r_mm + r_bul:
        rep.add(gating, "counts-mismatch", row,
                "reported %d mm + %d bulges = %d; the alignment has %d mm, %d DNA and "
                "%d RNA bulges" % (r_mm, r_bul, r_tot, mm, dn, rn))
    want = ",".join([t for t, c in (("RNA", rn), ("DNA", dn)) if c]) or "X"
    if row["btype"] != want and not (want == "X" and row["btype"] in ("", "NA")):
        rep.add(gating, "bulge-type-mismatch", row,
                "reported %s, the alignment has %d DNA and %d RNA bulges"
                % (row["btype"], dn, rn))
    if mm > ctx.mm or dn > ctx.bdna or rn > ctx.brna or mm + dn + rn > ctx.cap:
        rep.add(gating, "budget-exceeded", row,
                "%d mm, %d DNA, %d RNA bulges against %d/%d/%d, cap %s"
                % (mm, dn, rn, ctx.mm, ctx.bdna, ctx.brna, ctx.args.max_total_edits))
    pam_seq = stripped[:pam_len] if p5 else stripped[-pam_len:]
    if not all(iupac_match(m, b) for m, b in zip(ctx.pam_motif, pam_seq)):
        rep.add(gating, "pam-mismatch", row, "PAM %s against the motif %s"
                % (pam_seq, ctx.pam_motif))
    if row["pam"] and row["pam"] != "NA" and row["pam"].upper() != pam_seq:
        rep.add(gating, "pam-column-mismatch", row,
                "the PAM column says %s, the target has %s" % (row["pam"], pam_seq))
    return True


def ref_text_check(row, win, rep, gating):
    stripped = row["alt"].replace("-", "").upper()
    span = win.refseq[row["start"]:row["start"] + len(stripped)]
    if len(span) != len(stripped):
        rep.add(False, "limit-contig-end", row,
                "the REF text cannot be read: the target runs past the contig's end")
        return
    text = list(span if row["strand"] == "+" else rc(span))
    for i, c in enumerate(row["alt"]):
        if c == "-":
            text.insert(i, "-")
    if row["ref"].upper() != "".join(text).upper():
        rep.add(gating, "ref-text-mismatch", row,
                "REF %s, the reference at %d gapped as the target is %s"
                % (row["ref"], row["start"], "".join(text)))


def named_alleles(row, win, rep, gating):
    """The row's variants as Allele objects (the first line holding each), or None after
    a finding."""
    named = []
    for t in [x for x in row["var"].split(",") if x and x != "NA"]:
        try:
            c, p, r, a = parse_token(t)
        except (ValueError, IndexError):
            rep.add(gating, "variant-unparseable", row, t)
            return None
        found = win.by_key.get((p, r, a)) if c == row["chrom"] else None
        if not found:
            rep.add(gating, "variant-not-in-vcf", row,
                    "%s is not a record of the VCF (forward strand)" % t)
            return None
        if any(x in win.shadowed for x in found) and not any(
                x not in win.shadowed for x in found):
            rep.add(gating, "variant-shadowed", row,
                    "%s is an earlier SNP record of a split site, which CRISPRme's stores "
                    "never hold" % t)
        named.append(found[0])
    return named


# --------------------------------------------------------------------------- #
# a row's target: the reference with its named variants
# --------------------------------------------------------------------------- #
class Target(object):
    """A target holding ``main`` (the row's indel-path record): its + reading ``probe``,
    the places of its bases (on the sequence that makes it), the reference coordinates it
    spans, its SNV code columns (base index, coordinate), the window the checks write
    haplotypes over and the row's aligned strings (None in completeness when unknown)."""

    __slots__ = ("main", "named", "probe", "strand", "start", "place", "span", "codes",
                 "lo", "hi", "crrna", "aligned", "n", "src")

    def __init__(self, win, main, named, probe, strand, place, src, crrna, aligned):
        self.main = main
        self.named = named
        self.probe = probe
        self.strand = strand
        self.place = place
        self.src = src
        self.n = len(probe)
        self.start = place[0][0]
        self.span = span_of(place)
        self.codes = [(i, pl[0]) for i, (pl, s) in enumerate(zip(place, src))
                      if pl[1] == 0 and s is not main and win.snv_at.get(pl[0])]
        self.lo = max(0, min(self.span[0], main.lo) - 2 * self.n - 60)
        self.hi = min(len(win.refseq), max(self.span[1], main.hi) + 2 * self.n + 60)
        self.crrna = crrna
        self.aligned = aligned


def row_target(win, row, named, main):
    """The row's target on the reference with exactly its named variants: (Target,
    problem class, detail)."""
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    # a row of the record starts within one target of its REF span (inserted bases are
    # written after their anchor whatever the window's size)
    width = (main.hi - main.lo) if main else 1
    anchor = main.lo if main else row["start"]
    lo = max(0, min(row["start"], anchor) - n - width - 60)
    hi = min(len(win.refseq), max(row["start"], anchor + width) + 2 * n + width + 60)
    hap = make_hap(win.refseq, named, lo, hi)
    if hap.clash:
        return None, "variants-overlap", (
            "the named variants cannot all be on one haplotype: two of them describe the "
            "reference bases %s" % ",".join("%d..%d" % (z[0], z[1] - 1) for z in hap.clash))
    here = reads_at(hap, row["start"], probe)
    for i in here:
        if main is None or holds_edit(hap, i, n, main):
            return Target(win, main, named, probe, row["strand"], hap.place[i:i + n],
                          hap.src[i:i + n], row["crrna"], row["alt"]), None, None
    through = []
    k = hap.seq.find(probe)
    while k >= 0:
        if main is None or holds_edit(hap, k, n, main):
            through.append(hap.place[k][0])
        k = hap.seq.find(probe, k + 1)
    if through:
        return None, "start-not-on-variant-sequence", (
            "the target is the reference with its named variants%s starting at %s, not at "
            "%d" % (", through its record," if main else "",
                    ",".join(map(str, sorted(set(through))[:6])), row["start"]))
    if here:
        return None, "target-holds-no-indel-base", (
            "the target at %d holds no base its record %s changes (a flank target)"
            % (row["start"], main.token(win.chrom)))
    return None, "target-not-on-variant-sequence", (
        "the target is not the reference with its named variants near %d: a base no "
        "named variant explains" % row["start"])


# --------------------------------------------------------------------------- #
# the truth: the sample's haplotypes, as its calls give them
# --------------------------------------------------------------------------- #
def main_call(win, key, sample):
    """(the line's allele of ``key``, the sample's call) where the sample carries it."""
    for a in win.by_key.get(key, ()):
        c = win.call(a, sample)
        if c is not None and a.number in c[0]:
            return a, c
    return None, None


def record_slots(call, number):
    """(the slots holding allele ``number``, whether the call says which slots they
    are). An unphased heterozygous call holds it on one slot of unknown place (None)."""
    alle, phased = call
    hs = [h for h, x in enumerate(alle) if x == number]
    if not hs:
        return [], True
    if phased or len(alle) == 1 or all(x == number for x in alle):
        return hs, True
    return [None], False


def on_slot(call_main, h, call, number):
    """Whether slot h of the record's call also holds allele ``number`` of another line
    ('yes' / 'no' / '?'): known when both calls are phased alike, or the other call is
    homozygous or lacks the allele; unknown otherwise (an unphased call, a missing allele,
    a haploid call against a diploid one, or the record's own slot unknown)."""
    if call is None:
        return "no"
    alle, phased = call
    if number not in alle and None not in alle:
        return "no"
    ploidy = len(call_main[0])
    if len(alle) == ploidy and all(x == number for x in alle):
        return "yes"
    if h is None or not phased or not call_main[1] or len(alle) != ploidy:
        return "?"
    x = alle[h]
    if x is None:
        return "?"
    return "yes" if x == number else "no"


def slot_alleles(win, sample, main_al, call_main, h, lids):
    """(alleles the slot holds, alleles it may hold) among the lines ``lids``."""
    present, unknown = [], []
    for lid in lids:
        call = win.calls[lid].get(sample)
        if call is None:
            continue
        for a in win.line_alleles[lid]:
            if a is None or a is main_al:
                continue
            if lid == main_al.line:
                continue  # one allele of a line per slot: the slot's is the record's
            r = on_slot(call_main, h, call, a.number)
            if r == "yes":
                present.append(a)
            elif r == "?":
                unknown.append(a)
    return present, unknown


def truth(win, tgt, sample, lids):
    """What the sample's own haplotypes say of the row's target:
      ('none',)                       it carries no allele of the row's record
      ('exact', holds, info)          its alleles where the target reads are on known slots
      ('clash', info)                 a slot of the record holds alleles that overlap there
      ('unknown', info)               which alleles share the record's slot there is unknown
    ``holds``: some slot of the record reads the target from the row's start, through the
    record, and holds every named variant."""
    main_al, call = main_call(win, tgt.main.key, sample)
    if main_al is None:
        return ("none",)
    slots, _known = record_slots(call, main_al.number)
    named = [a for a in tgt.named if a.key != tgt.main.key]
    info = []
    verdicts = []
    for h in slots:
        present, unknown = slot_alleles(win, sample, main_al, call, h, lids)
        hap = make_hap(win.refseq, [main_al] + present, tgt.lo, tgt.hi)
        hits = reads_at(hap, tgt.start, tgt.probe, main_al)
        span = span_of(hap.place[hits[0]:hits[0] + tgt.n]) if hits else tgt.span
        lo_, hi_ = min(span[0], tgt.span[0]), max(span[1], tgt.span[1])
        unk = [a for a in unknown if touches(a, lo_, hi_)]
        clash = hap.clashes(lo_ - 1, hi_ + 1)
        keys = {a.key for a in present}
        has_named = all(a.key in keys for a in named)
        extra = [a for a in present if touches(a, tgt.span[0], tgt.span[1])
                 and a.key not in {x.key for x in named}]
        holds = bool(hits) and has_named
        info.append({"slot": h, "holds": holds, "named": has_named, "extra": extra,
                     "unknown": unk, "clash": clash})
        verdicts.append((holds, bool(unk), bool(clash)))
    if any(hd and not u and not c for hd, u, c in verdicts):
        return ("exact", True, info)
    if any(u for _hd, u, _c in verdicts):
        return ("unknown", info)
    if any(c for _hd, _u, c in verdicts):
        return ("clash", info)
    return ("exact", False, info)


def options(win, sample, q, visible_only=False):
    """The bases the sample's own alleles allow at reference coordinate q: the ALT of every
    SNV allele there it carries, and the reference base unless each of its slots holds an
    ALT there (a missing allele may be the reference). ``visible_only``: as CRISPRme's
    stores see them (an earlier line of a split site reads as the reference)."""
    alts = set()
    held = set()
    ploidy = 0
    missing = False
    for a in win.snv_at.get(q, ()):
        call = win.call(a, sample)
        if call is None:
            continue
        alle = call[0]
        ploidy = max(ploidy, len(alle))
        if visible_only and a in win.shadowed:
            continue
        for h, x in enumerate(alle):
            if x == a.number:
                alts.add(a.body)
                held.add(h)
            elif x is None:
                missing = True
    if missing or len(held) < max(ploidy, 1):
        alts.add(win.refseq[q])
    return alts


def strongest(ctx, tgt, opts_by_code):
    """The documented PUTATIVE choice (METHODS section 4) over the target's code columns:
    {base index: set of allowed bases (+ orientation)} for the spacer columns, and the set
    of allowed tuples for the PAM's code columns (jointly), plus whether every choice is
    determined. The guide's base wins wherever it is an option; elsewhere, and at the PAM,
    the options of highest CFD (every option when CFD does not apply)."""
    if tgt.crrna is None or tgt.aligned is None:
        return None
    n_cols = len(tgt.crrna)
    cols = base_columns(tgt.aligned, tgt.strand)
    pam = set(ctx.pam_cols(n_cols))
    positions = Cfd.positions(tgt.crrna, n_cols - ctx.pam_len)
    ori = (lambda b: b) if tgt.strand == "+" else comp
    allowed = {}
    pam_codes = []
    determined = True
    for i, _q in tgt.codes:
        opts = opts_by_code[i]
        col = cols[i]
        if col in pam:
            pam_codes.append(i)
            continue
        g = tgt.crrna[col].upper()
        o_opts = {ori(b): b for b in opts}
        if g != "-" and g in o_opts:
            allowed[i] = {o_opts[g]}
            continue
        if ctx.cfd is None or len(opts) == 1:
            allowed[i] = set(opts)
            determined = determined and len(opts) == 1
            continue
        score = {b: ctx.cfd.factor(g, ori(b), positions[col]) for b in opts}
        top = max(score.values())
        allowed[i] = {b for b, v in score.items() if v >= top - 1e-12}
        determined = determined and len(allowed[i]) == 1
    pam_ok = None
    if pam_codes:
        target_bases = list(tgt.aligned.replace("-", "").upper())
        # target_bases is in the target's orientation: base index i (+ order) sits at
        # position i on the + strand, n-1-i on the - strand
        def at(i):
            return i if tgt.strand == "+" else len(target_bases) - 1 - i
        combos = []
        for combo in itertools.product(*[sorted(opts_by_code[i]) for i in pam_codes]):
            tb = list(target_bases)
            for i, b in zip(pam_codes, combo):
                tb[at(i)] = ori(b)
            pseq = "".join(tb[:ctx.pam_len] if ctx.pam_5prime else tb[-ctx.pam_len:])
            if not all(iupac_match(m, x) for m, x in zip(ctx.pam_motif, pseq)):
                continue
            score = ctx.cfd.pam.get(pseq[-2:], 0.0) if ctx.cfd is not None else 0.0
            combos.append((score, combo))
        if combos:
            top = max(s for s, _c in combos)
            pam_ok = {c for s, c in combos if s >= top - 1e-12}
        else:
            pam_ok = set()
        determined = determined and len(pam_ok) == 1
    return allowed, pam_codes, pam_ok, determined


def putative_check(win, ctx, tgt, sample, visible_only=False):
    """(allowed, strongest, determined) of the row's bases for a sample read by the
    PUTATIVE rule: every base at a code column is an allele of the sample; the bases are
    its strongest choice; that choice is one only (no tie the rule leaves open)."""
    opts = {i: options(win, sample, q, visible_only) for i, q in tgt.codes}
    row = {i: tgt.probe[i] for i, _q in tgt.codes}
    allowed = all(row[i] in opts[i] for i in opts)
    st = strongest(ctx, tgt, opts)
    if st is None:
        return allowed, allowed, False
    per, pam_codes, pam_ok, determined = st
    strong = all(row[i] in per[i] for i in per)
    if pam_codes:
        strong = strong and tuple(row[i] for i in pam_codes) in pam_ok
    return allowed, allowed and strong, determined


# --------------------------------------------------------------------------- #
# CRISPRme's documented view of a carrier (classification only)
# --------------------------------------------------------------------------- #
def documented(win, ctx, tgt, sample, lids):
    """What CRISPRme's documented rules make of a carrier of the row's record for the
    row's target (METHODS section 4; CHANGELOG 2.8.1 Known limits), re-implemented from
    that text: (verdict, causes), verdict 'lists', 'off' or 'maybe' (PUTATIVE, and the
    rule leaves the choice open), causes a set naming what decided it."""
    main_al, call = main_call(win, tgt.main.key, sample)
    if main_al is None:
        return "off", {"no-record"}
    causes = set()
    if main_al.kind == "mnv":
        causes.add("putative-mnv")  # MNVs are not in the phased indel store
        return _documented_putative(win, ctx, tgt, sample, causes)
    slots, known = record_slots(call, main_al.number)
    if not known:
        causes.add("putative")
        return _documented_putative(win, ctx, tgt, sample, causes)
    s0, s1 = tgt.span
    verdicts = []
    for h in slots:
        snvs, seconds, readable = [], [], True
        off = None
        for lid in lids:
            c = win.calls[lid].get(sample)
            if c is None or lid == main_al.line:
                continue
            for a in win.line_alleles[lid]:
                if a is None:
                    continue
                r = on_slot(call, h, c, a.number)
                if r == "no":
                    continue
                meets = a.lo < main_al.hi and a.hi > main_al.lo
                if meets:
                    if a.key[0] == main_al.key[0]:
                        causes.add("artefact")  # another allele at its own position
                        continue
                    if (not main_al.body and not a.adds() and a.lo >= main_al.edit_lo
                            and a.hi <= main_al.edit_hi):
                        causes.add("deleted")  # only on bases the record deletes
                        continue
                    if r == "?":
                        causes.add("artefact")  # a multi-allelic heterozygote
                        continue
                    if a.kind == "snv":
                        off = off or "undocumented"
                    else:
                        off = "complex"
                    continue
                if a.kind == "snv":
                    if not touches(a, s0, s1):
                        continue
                    if r == "?":
                        causes.add("putative")
                        return _documented_putative(win, ctx, tgt, sample, causes)
                    if a in win.untied:
                        causes.add("putative-untied")
                        return _documented_putative(win, ctx, tgt, sample, causes)
                    if a in win.shadowed:
                        causes.add("last-alt-wins")
                        continue
                    snvs.append(a)
                elif a.kind == "indel":
                    if r == "?":
                        if changes_target(a, tgt):
                            causes.add("putative")
                            return _documented_putative(win, ctx, tgt, sample, causes)
                        continue
                    seconds.append(a)
                else:  # an MNV: not in the phased indel store, not written
                    if changes_target(a, tgt):
                        causes.add("second-indel")
        if off:
            causes.add(off)
            verdicts.append(False)
            continue
        inside = [a for a in seconds if changes_target(a, tgt)]
        window = make_hap(win.refseq, [main_al] + snvs, tgt.lo, tgt.hi)
        if not inside:
            verdicts.append(bool(reads_at(window, tgt.start, tgt.probe, main_al)))
            continue
        causes.add("second-indel")
        # the slot's bases, with its phased store indels and its own SNP alleles around
        snv_all = []
        unread = []
        for lid in lids:
            c = win.calls[lid].get(sample)
            if c is None or lid == main_al.line:
                continue
            for a in win.line_alleles[lid]:
                if a is None or a.kind != "snv" or a in win.shadowed:
                    continue
                if a.lo < main_al.hi and a.hi > main_al.lo:
                    continue
                r = on_slot(call, h, c, a.number)
                if r == "yes":
                    snv_all.append(a)
                elif r == "?":
                    unread.append(a)
        hap = make_hap(win.refseq, [main_al] + snv_all + seconds, tgt.lo, tgt.hi)
        ok = False
        if not hap.clash or not hap.clashes(s0 - tgt.n, s1 + tgt.n):
            n = tgt.n
            i = hap.index(tgt.place[0])
            if i is not None and hap.seq[i:i + n] == tgt.probe and holds_edit(
                    hap, i, n, main_al):
                sp = span_of(hap.place[i:i + n])
                ok = not any(touches(a, sp[0], sp[1]) for a in unread)
            j = hap.index(tgt.place[-1])
            if not ok and j is not None and j - n + 1 >= 0 and \
                    hap.seq[j - n + 1:j + 1] == tgt.probe and holds_edit(
                        hap, j - n + 1, n, main_al):
                sp = span_of(hap.place[j - n + 1:j + 1])
                ok = not any(touches(a, sp[0], sp[1]) for a in unread)
        verdicts.append(ok)
    return ("lists" if any(verdicts) else "off"), causes


def _documented_putative(win, ctx, tgt, sample, causes):
    allowed, strong, determined = putative_check(win, ctx, tgt, sample, visible_only=True)
    if any(a in win.shadowed and win.call(a, sample) is not None
           for _i, q in tgt.codes for a in win.snv_at.get(q, ())):
        causes.add("last-alt-wins")
    if not (allowed and strong):
        return "off", causes
    return ("lists" if determined else "maybe"), causes


LIMIT_OF = (("last-alt-wins", "limit-last-alt-wins"),
            ("second-indel", "limit-two-indels"),
            ("putative-mnv", "limit-putative-carrier"),
            ("putative-untied", "limit-putative-carrier"),
            ("putative", "limit-putative-carrier"),
            ("complex", "limit-overlapping-records"),
            ("undocumented", "limit-overlapping-records"),
            ("artefact", "limit-overlapping-records"),
            ("deleted", "limit-overlapping-records"))


def limit_of(causes):
    for c, cls in LIMIT_OF:
        if c in causes:
            return cls
    return None


# --------------------------------------------------------------------------- #
# soundness
# --------------------------------------------------------------------------- #
CARRIER_TEXT = {
    "carrier-without-indel": "listed without carrying %s",
    "carrier-without-target": "listed, and no haplotype of theirs holds the target at its "
                              "start with %s",
    "carrier-without-variants": "listed, and no slot of theirs holding %s holds every "
                                "named variant",
    "variant-missing": "listed, and their slot of %s holds a variant inside the target "
                       "the row does not name",
    "carrier-without-allele": "listed (phase unknown there) for an allele they do not "
                              "carry, with %s",
    "carrier-not-strongest": "listed (phase unknown there) for a weaker target than the "
                             "strongest their alleles allow, with %s (METHODS section 4: "
                             "the guide's base wherever they carry it, else the highest "
                             "CFD)",
    "carrier-off-by-rule": "listed although their slot of %s holds a record overlapping it "
                           "at another position (CRISPRme documents taking such a slot "
                           "off)",
    "carrier-missing": "hold the target at its start with %s (or it is their strongest), "
                       "are listed by CRISPRme's documented rules, and are not on the row",
    "unexplained-disagreement": "their haplotypes and CRISPRme's documented rules part on "
                                "this row of %s for no documented reason",
    "limit-last-alt-wins": "differ from the row by an earlier SNP record of a split site, "
                           "which CRISPRme does not read, with %s",
    "limit-two-indels": "need a second indel-path record on the slot of %s, or have one "
                        "read by the documented rule",
    "limit-overlapping-records": "hold alleles that overlap %s or each other where the "
                                 "target reads (no one sequence in the VCF there)",
    "limit-putative-carrier": "are carriers of %s that CRISPRme resolves PUTATIVE (an MNV, "
                              "an untied multi-ALT alt, an unphased call) while their "
                              "haplotype says otherwise",
}


def carrier_verdict(win, ctx, tgt, sample, listed, lids):
    """(gating, class) for one sample of a row, or None when it agrees."""
    tv = truth(win, tgt, sample, lids)
    kind = tv[0]
    if kind == "none":
        return (True, "carrier-without-indel") if listed else None
    if kind == "exact":
        holds, info = tv[1], tv[2]
        if holds == listed:
            return None
        verdict, causes = documented(win, ctx, tgt, sample, lids)
        if holds:  # on its haplotype, not on the row
            if verdict == "lists":
                return True, "carrier-missing"
            if verdict == "maybe":
                return None
            return unexplained_or(causes)
        if verdict in ("lists", "maybe"):
            return unexplained_or(causes)
        if "complex" in causes:
            return True, "carrier-off-by-rule"
        return True, diagnose(info)
    if kind == "clash":
        verdict, causes = documented(win, ctx, tgt, sample, lids)
        if listed and verdict == "off" and "complex" in causes:
            return True, "carrier-off-by-rule"
        zones = [z for x in tv[1] for z in x["clash"]]
        return False, "limit-overlapping-records", "; ".join(sorted(clash_text(win, zones)))
    # unknown phase where the target reads: the documented PUTATIVE rule is the truth
    allowed, strong, determined = putative_check(win, ctx, tgt, sample)
    if listed:
        if strong:
            return None
        verdict, causes = documented(win, ctx, tgt, sample, lids)
        if verdict in ("lists", "maybe"):
            return unexplained_or(causes)
        return True, ("carrier-without-allele" if not allowed else "carrier-not-strongest")
    if not (strong and determined):
        return None
    verdict, causes = documented(win, ctx, tgt, sample, lids)
    if verdict == "lists":
        return True, "carrier-missing"
    if verdict == "maybe":
        return None
    return unexplained_or(causes)


def unexplained_or(causes):
    """The limit a documented cause makes of a disagreement; with none, the disagreement
    is unexplained (gating: a defect of CRISPRme or of this checker)."""
    cls = limit_of(causes)
    if cls:
        return False, cls, "documented: %s" % ",".join(sorted(causes))
    return True, "unexplained-disagreement"


def diagnose(info):
    """The class of a listed carrier whose known slots of the record lack the target."""
    if any(x["named"] for x in info):
        if any(x["named"] and x["extra"] for x in info):
            return "variant-missing"
        return "carrier-without-target"
    return "carrier-without-variants"


def local_lines(win, tgt):
    return win.lines_in(tgt.lo, tgt.hi)


def check_row(row, win, ctx, rep, gating=True):
    """Every soundness check of one indel/MNV row block."""
    if row["origin"] != "alt":
        rep.add(gating, "origin-not-alt", row,
                "REF/ALT_origin is %s on a row naming an indel-path record" % row["origin"])
    if not alignment_checks(row, ctx, rep, gating):
        return
    ref_text_check(row, win, rep, gating)
    named = named_alleles(row, win, rep, gating)
    if named is None:
        return
    mains = [a for a in named if a.kind != "snv"]
    if len(mains) != 1:
        rep.add(gating, "variants-not-one-indel", row,
                "the row names %d indel/MNV records; a fake-indel window holds one"
                % len(mains))
        return
    main = mains[0]
    samples = samples_of(row)
    if not samples:
        rep.add(gating, "no-carrier-listed", row, "the row lists no carrier")
        return
    tgt, cls, detail = row_target(win, row, named, main)
    if tgt is None:
        rep.add(gating, cls, row, detail)
        return
    for a in named:
        if a is not main and not any(s is a for s in tgt.src):
            rep.add(gating, "variant-outside-target", row,
                    "%s changes no base of the target (%d..%d)"
                    % (a.token(win.chrom), tgt.span[0], tgt.span[1]))
    if not win.samples:
        return  # a sites-only VCF: no genotype to judge a carrier by (noted once by main)
    lids = local_lines(win, tgt)
    listed = set(samples)
    pool = set(listed)
    for a in win.by_key.get(main.key, ()):
        pool.update(win.carriers(a))
    out = defaultdict(list)
    why = defaultdict(set)
    cache = {}
    for s in sorted(pool):
        sig = (s in listed, win.signature(s, lids))
        v = cache.get(sig)
        if v is None:
            v = carrier_verdict(win, ctx, tgt, s, s in listed, lids) or ("ok",)
            cache[sig] = v
        if v != ("ok",):
            out[v[:2]].append(s)
            if len(v) > 2 and v[2]:
                why[v[:2]].add(v[2])
    tok = main.token(win.chrom)
    for (g, cls), who in sorted(out.items(), key=lambda kv: kv[0][1]):
        extra = sorted(why.get((g, cls), ()))
        detail = "%d sample(s) %s%s: %s" % (
            len(who), CARRIER_TEXT.get(cls, "%s") % tok,
            " (%s%s)" % ("; ".join(extra[:4]), "; ..." if len(extra) > 4 else "")
            if extra else "", ",".join(who[:12]))
        rep.add(g and gating, cls, row, detail, samples=",".join(who))


def check_snp_row(row, win, ctx, rep):
    """The same checks on a SNP-only row of a per-sample run, report-only (CHANGELOG
    2.8.1 Known limits: a SNP row lists carriers from SNP genotypes only)."""
    if not alignment_checks(row, ctx, rep, False):
        return
    ref_text_check(row, win, rep, False)
    named = named_alleles(row, win, rep, False)
    if named is None:
        return
    hap = None
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    lo = max(0, row["start"] - 2 * n - 60)
    hi = min(len(win.refseq), row["start"] + 3 * n + 60)
    hap = make_hap(win.refseq, named, lo, hi)
    hits = reads_at(hap, row["start"], probe)
    if not hits:
        rep.add(False, "target-not-on-variant-sequence", row,
                "the target is not the reference with its named variants at %d"
                % row["start"])
        return
    span = span_of(hap.place[hits[0]:hits[0] + n])
    lids = win.lines_in(lo, hi)
    out = defaultdict(list)
    for s in samples_of(row):
        calls = {lid: win.calls[lid][s] for lid in lids if s in win.calls[lid]}
        if not calls:
            out["carrier-without-target"].append(s)
            continue
        ploidies = {len(c[0]) for c in calls.values()}
        touching = [lid for lid in calls if win.line_lo[lid] <= span[1]
                    and win.line_hi[lid] - 1 >= span[0]]
        if len(ploidies) > 1 or any(
                (not calls[lid][1] and len(set(calls[lid][0])) > 1)
                or None in calls[lid][0] for lid in touching):
            out["limit-snp-row-unphased"].append(s)
            continue
        ploidy = ploidies.pop()
        ok, why = False, None
        for h in range(ploidy):
            keys = [a for lid, c in calls.items() for a in win.line_alleles[lid]
                    if a is not None and c[0][h] == a.number]
            if not all(any(a is k for k in keys) or any(a.key == k.key for k in keys)
                       for a in named):
                continue
            hh = make_hap(win.refseq, keys, lo, hi)
            if reads_at(hh, row["start"], probe):
                ok = True
                break
            inside = [a for a in keys if touches(a, span[0], span[1])
                      and all(a.key != x.key for x in named)]
            if any(a.kind != "snv" for a in inside):
                why = why or "limit-snp-row-indel"
            elif any(a in win.shadowed for a in inside):
                why = why or "limit-last-alt-wins"
            elif hh.clashes(span[0] - 1, span[1] + 1):
                why = why or "limit-overlapping-records"
        if not ok:
            out[why or "carrier-without-target"].append(s)
    text = {"limit-snp-row-indel": "hold an indel (or MNV) inside the target on the slot "
                                   "of the row's SNPs: a SNP row's carriers come from SNP "
                                   "genotypes only (CHANGELOG 2.8.1, Known limits)",
            "limit-last-alt-wins": "hold an earlier SNP record of a split site inside the "
                                   "target (CRISPRme reads the last record)",
            "limit-overlapping-records": "hold overlapping alleles inside the target",
            "limit-snp-row-unphased": "have a call of unknown phase inside the target: a "
                                      "SNP row lists every subset a carrier may hold "
                                      "(METHODS section 4), not checked here",
            "carrier-without-target": "hold no haplotype with this target and every named "
                                      "variant"}
    for cls, who in sorted(out.items()):
        rep.add(False, cls, row, "%d listed sample(s) %s: %s"
                % (len(who), text[cls], ",".join(who[:12])), samples=",".join(who))


# --------------------------------------------------------------------------- #
# completeness
# --------------------------------------------------------------------------- #
class Demand(object):
    """A target a carrier must be reported for: exact (its known haplotype) or the
    strongest of its PUTATIVE choices (any one of ``alts`` when the rule leaves a tie)."""

    __slots__ = ("kind", "alts", "total", "cfd")

    def __init__(self, kind, alts, total, cfd):
        self.kind = kind
        self.alts = alts  # [(Target, variant tokens frozenset)]
        self.total = total
        self.cfd = cfd


def completeness(win, ctx, rep, reported, region, popmode):
    """The real-data brute force over the indel and MNV records of ``region``."""
    lo_r, hi_r = region
    for guide in ctx.guides:
        full = ctx.full(guide)
        one = len(full) + max(ctx.bdna, ctx.brna) + 2   # one target and its bulges
        reach = 2 * one + 20
        for al in list(win.alleles):
            if al.kind == "snv" or not (lo_r <= al.lo < hi_r):
                continue
            if win.by_key[al.key][0] is not al:
                continue  # one record per key (a second line of it is the same record)
            carriers = set()
            for a in win.by_key[al.key]:
                carriers.update(win.carriers(a))
            if not carriers:
                continue
            lo = max(0, al.lo - reach)
            hi = min(len(win.refseq), al.hi + reach)
            if al.lo - one < 0 or al.hi + one > len(win.refseq):
                rep.add(False, "limit-contig-end", {
                    "table": "completeness", "chrom": win.chrom, "start": al.lo,
                    "var": al.token(win.chrom)},
                    "the record is within one target of the contig's end")
                continue
            comp_record(win, ctx, rep, reported, al, carriers, guide, full, lo, hi,
                        popmode, (al.edit_lo - one, al.edit_hi + one))


def scan_cached(cache, ctx, full, seq, skipped):
    hits = cache.get(seq)
    if hits is None:
        hits = list(scan_sequence(full, seq, ctx.mm, ctx.bdna, ctx.brna, ctx.cap,
                                  ctx.pam_len, ctx.pam_5prime, skipped))
        cache[seq] = hits
    return hits


def demand_tokens(win, al, place, src, present):
    """The variants CRISPRme names for a target of the record: the record and the SNV
    alleles whose bases it holds; and the other indel-path alleles that change it (one
    that wrote a base of it, or whose edit lies inside it)."""
    toks = {al.token(win.chrom)}
    others = []
    s0, s1 = span_of(place)
    for s in src:
        if s is not None and s is not al and s.kind == "snv":
            toks.add(s.token(win.chrom))
    for a in present:
        if a is not al and a.kind != "snv" and (any(x is a for x in src)
                                                or touches(a, s0, s1)):
            others.append(a)
    return frozenset(toks), others


def comp_record(win, ctx, rep, reported, al, carriers, guide, full, lo, hi, popmode,
                reach=None):
    """The demands of one record. ``reach``: the reference coordinates a target holding
    the record's edit can cover; an unknown call outside it cannot change such a target,
    so only the calls inside it are tried both ways."""
    reach = reach or (lo, hi)
    chrom = win.chrom
    lids = win.lines_in(lo, hi)
    groups = defaultdict(list)
    for s in carriers:
        _a, call = main_call(win, al.key, s)
        groups[(call, win.signature(s, lids))].append(s)
    cache = {}
    skipped = []
    demands = defaultdict(set)   # (kind, key of the demand) -> samples
    demand_of = {}
    limits = defaultdict(set)    # (class, strand, start, probe, why) -> samples
    clash_who = set()
    clash_how = set()
    capped = set()
    tok = al.token(chrom)
    for (call, _sig), members in groups.items():
        s = members[0]
        main_al, call = main_call(win, al.key, s)
        slots, known = record_slots(call, main_al.number)
        need_putative = not known
        for h in (slots if known else ()):
            present, unknown = slot_alleles(win, s, main_al, call, h, lids)
            by_line = defaultdict(list)
            for a in unknown:
                if touches(a, reach[0], reach[1]) or (a.lo <= reach[1] and a.hi > reach[0]):
                    by_line[a.line].append(a)
            choices = [[None] + v for v in by_line.values()]
            if len(choices) > PUTATIVE_CAP:
                capped.update(members)
                choices = [[None]]
            for combo in itertools.product(*choices):
                chosen = [a for a in combo if a is not None]
                hap = make_hap(win.refseq, [main_al] + present + chosen, lo, hi)
                own = hap.clashes(main_al.lo, main_al.hi - 1)
                if own:
                    clash_who.update(members)  # the record itself overlaps another
                    clash_how.update(clash_text(win, own))
                    continue
                for strand, off, sa, ta, mm, dn, rn in scan_cached(
                        cache, ctx, full, hap.seq, skipped):
                    width = len(ta.replace("-", ""))
                    if off + width > len(hap.seq) or not holds_edit(hap, off, width,
                                                                    main_al):
                        continue
                    place = hap.place[off:off + width]
                    src = hap.src[off:off + width]
                    s0, s1 = span_of(place)
                    if any(touches(a, s0, s1) for a in unknown):
                        need_putative = True
                        continue
                    near = hap.clashes(s0 - 1, s1 + 1)
                    if near:
                        clash_who.update(members)
                        clash_how.update(clash_text(win, near))
                        continue
                    probe = hap.seq[off:off + width]
                    toks, others = demand_tokens(win, main_al, place, src,
                                                 present + chosen)
                    t = Target(win, main_al, [main_al], probe, strand, place, src, sa, ta)
                    # a target that needs a second indel-path record on the slot is
                    # demanded too: a row reporting it (the second one leaves the bases
                    # read from the target's first base) satisfies it, and a miss is the
                    # documented limit (one record per fake-indel window)
                    key = ("second" if others else "exact", strand, t.start, probe, toks)
                    demands[key].update(members)
                    keep_demand(demand_of, key, t, toks, mm + dn + rn,
                                target_cfd(ctx, sa, ta))
        if need_putative:
            putative_targets(win, ctx, rep, al, main_al, s, members, lids, lo, hi, full,
                             known, cache, skipped, demands, demand_of, capped, reach)
    resolve(win, ctx, rep, reported, al, carriers, lids, demands, demand_of, popmode)
    for (cls, strand, start, probe, why), who in sorted(limits.items(),
                                                        key=lambda kv: kv[0][:4]):
        who = sorted(who)
        rep.add(False, cls, {"table": "completeness", "chrom": chrom, "start": start,
                             "strand": strand, "alt": probe, "var": tok},
                "%d sample(s): %s: %s" % (len(who), why, ",".join(who[:12])),
                samples=",".join(who))
    if clash_who:
        who = sorted(clash_who)
        how = sorted(clash_how)
        rep.add(False, "limit-overlapping-records", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d sample(s) hold %s on a slot whose alleles overlap where some of its targets "
            "read, so the VCF gives no one sequence there and they are not demanded "
            "(overlapping on one slot: %s%s): %s"
            % (len(who), tok, "; ".join(how[:4]), "; ..." if len(how) > 4 else "",
               ",".join(who[:12])), samples=",".join(who))
    if capped:
        who = sorted(capped)
        rep.add(False, "limit-putative-cap", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d sample(s) have more than %d positions of unknown phase within a target's "
            "reach of %s: their targets there that touch one are not enumerated: %s"
            % (len(who), PUTATIVE_CAP, tok, ",".join(who[:12])), samples=",".join(who))
    if skipped:
        rep.add(False, "limit-n-bases", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d window(s) around %s hold an N and are not aligned (no base pairs with N)"
            % (len(skipped), tok))


def keep_demand(demand_of, key, t, toks, total, cfd):
    """Record one alignment of a demanded target: the target is one sequence at one
    start, and its many alignments differ in edits and CFD, so the demand keeps the best
    of each (what the cluster's row must match in a run with one row per cluster)."""
    d = demand_of.get(key)
    if d is None:
        demand_of[key] = Demand(key[0], [(t, toks)], total, cfd)
        return
    d.total = min(d.total, total)
    if cfd is not None and (d.cfd is None or cfd > d.cfd):
        d.cfd = cfd


def clash_text(win, zones):
    return {"+".join(a.token(win.chrom) for a in z[2]) for z in zones}


def target_cfd(ctx, sa, ta):
    if ctx.cfd is None:
        return None
    return ctx.cfd.score(sa, ta.upper(), ctx.pam_len)


def putative_targets(win, ctx, rep, al, main_al, s, members, lids, lo, hi, full, known,
                     cache, skipped, demands, demand_of, capped, reach=None):
    """The PUTATIVE demands of one carrier group: every in-budget target its own alleles
    allow with the record (any choice at each SNV position: the documented rule ignores
    phase), grouped by alignment; the strongest choice of each is demanded. With the
    record's slot known, only alignments that touch a call of unknown phase are."""
    main_call_ = main_call(win, al.key, s)[1]
    unknown_all = []
    if known:
        for h in record_slots(main_call_, main_al.number)[0]:
            unknown_all.extend(slot_alleles(win, s, main_al, main_call_, h, lids)[1])
    r_lo, r_hi = reach or (lo, hi)
    positions = sorted({a.lo for a in win.alleles if a.kind == "snv"
                        and max(lo, r_lo) <= a.lo <= min(hi - 1, r_hi)})
    choice = []
    for q in positions:
        opts = options(win, s, q)
        if len(opts) > 1:
            choice.append((q, sorted(opts)))
    if len(choice) > PUTATIVE_CAP:
        capped.update(members)
        return
    fixed = []
    for q in positions:
        opts = options(win, s, q)
        if len(opts) == 1:
            b = next(iter(opts))
            if b != win.refseq[q]:
                fixed.append(next(a for a in win.snv_at[q] if a.body == b))
    groups = defaultdict(dict)  # alignment -> {probe: (Target, mm+dn+rn)}
    for combo in itertools.product(*[opts for _q, opts in choice]):
        chosen = list(fixed)
        for (q, _opts), b in zip(choice, combo):
            if b != win.refseq[q]:
                chosen.append(next(a for a in win.snv_at[q] if a.body == b))
        hap = make_hap(win.refseq, [main_al] + chosen, lo, hi)
        for strand, off, sa, ta, mm, dn, rn in scan_cached(cache, ctx, full, hap.seq,
                                                           skipped):
            width = len(ta.replace("-", ""))
            if off + width > len(hap.seq) or not holds_edit(hap, off, width, main_al):
                continue
            place = hap.place[off:off + width]
            s0, s1 = span_of(place)
            if known and not any(touches(a, s0, s1) for a in unknown_all):
                continue
            probe = hap.seq[off:off + width]
            gap = "".join("-" if c == "-" else "x" for c in ta)
            g = groups[(strand, tuple(place), sa, gap)]
            if probe not in g:
                t = Target(win, main_al, [main_al], probe, strand, place,
                           hap.src[off:off + width], sa, ta)
                g[probe] = (t, mm + dn + rn)
    for (strand, _place, sa, _gap), members_t in groups.items():
        best = []
        for probe, (t, total) in members_t.items():
            if putative_check(win, ctx, t, s)[1]:
                best.append((t, total))
        if not best:
            continue
        alts = []
        for t, _total in best:
            toks = {al.token(win.chrom)}
            for x in t.src:
                if x is not None and x is not main_al and x.kind == "snv":
                    toks.add(x.token(win.chrom))
            alts.append((t, frozenset(toks)))
        t0 = alts[0][0]
        key = ("putative", strand, t0.start, tuple(sorted(t.probe for t, _x in alts)),
               tuple(sorted(tuple(sorted(x)) for _t, x in alts)))
        demands[key].update(members)
        best_total = min(tt for _t, tt in best)
        best_cfd = (max((target_cfd(ctx, t.crrna, t.aligned) or 0.0) for t, _x in alts)
                    if ctx.cfd else None)
        if key not in demand_of:
            demand_of[key] = Demand("putative", alts, best_total, best_cfd)
        else:
            keep_demand(demand_of, key, alts[0][0], alts[0][1], best_total, best_cfd)


def resolve(win, ctx, rep, reported, al, carriers, lids, demands, demand_of, popmode):
    """Each demand against the tables; a miss is classified by the documented view."""
    chrom = win.chrom
    full_len = len(ctx.full(ctx.guides[0]))
    fold_sets = {}
    for key in sorted(demands, key=lambda k: (k[0], k[1], k[2], str(k[3]))):
        d = demand_of[key]
        who = sorted(demands[key])
        t0, toks0 = d.alts[0]
        cl = cluster_of(t0.start, t0.strand, t0.n, full_len, ctx.pam_5prime)
        row = {"table": "completeness", "chrom": chrom, "start": t0.start,
               "strand": t0.strand, "crrna": t0.crrna, "alt": t0.probe,
               "var": ",".join(sorted(toks0))}
        if popmode:
            # one row per cluster: the cluster must hold a row as good as the target, for
            # every carrier CRISPRme's documented rules report it for
            if d.kind == "second":
                ok, why = reported.locus(chrom, cl, d.total, d.cfd)
                if not ok:
                    rep.add(False, "limit-two-indels", row,
                            "%d sample(s) hold this target only with a second "
                            "indel-path record on the slot (one record per fake-indel "
                            "window), and its cluster %s: %s"
                            % (len(who), why, ",".join(who[:12])), samples=",".join(who))
                continue
            ok, why = reported.locus(chrom, cl, d.total, d.cfd)
            if ok:
                continue
            # CRISPRme must report the target when its documented rules list any one of
            # the carriers for it (one call signature is enough to ask)
            causes, reach = set(), None
            seen = set()
            for s in who:
                sig = win.signature(s, lids)
                if sig in seen:
                    continue
                seen.add(sig)
                verdict, cs = documented_for(win, ctx, d, s, lids)
                causes |= cs
                if verdict != "off":
                    reach = s
                    break
            if reach is None:
                g, cls = unexplained_or(causes)[:2]
                rep.add(g, cls, row,
                        "%d sample(s): their cluster %s, and CRISPRme's documented rules "
                        "do not report this target for them (documented: %s): %s"
                        % (len(who), why, ",".join(sorted(causes)), ",".join(who[:12])),
                        samples=",".join(who))
                continue
            rep.add(True, "missing-locus", row,
                    "%d sample(s) hold this target (%d edits%s) at Cluster Position %d "
                    "(%s among them, listed by CRISPRme's documented rules), and %s"
                    % (len(who), d.total, "" if d.cfd is None else ", CFD %.4f" % d.cfd,
                       cl, reach, why), samples=",".join(who))
            continue
        missing = []
        for s in who:
            if any(reported.has(s, t.strand, t.start, t.probe, toks) for t, toks in d.alts):
                continue
            hit = False
            for t, toks in d.alts:
                fk = (t.strand, t.start, t.probe, toks)
                if fk not in fold_sets:
                    fold_sets[fk] = documented_carriers(win, ctx, t, carriers, lids)
                cm, pm = fold_sets[fk]
                if reported.folded(s, cl, toks, cm | set(who), pm):
                    hit = True
                    break
            if not hit:
                missing.append(s)
        if not missing:
            continue
        out = defaultdict(list)
        if d.kind == "second":
            rep.add(False, "limit-two-indels", row,
                    "%d sample(s) hold this target only with a second indel-path record "
                    "on the slot, and no row reports it for them (one record per "
                    "fake-indel window): %s" % (len(missing), ",".join(missing[:12])),
                    samples=",".join(missing))
            continue
        for s in missing:
            verdict, causes = documented_for(win, ctx, d, s, lids)
            if verdict == "lists":
                out[(True, "missing-row" if d.kind == "exact" else "missing-putative-row",
                     "")].append(s)
            elif verdict == "maybe":
                if reported.listed_at(s, cl, al.token(chrom)):
                    continue
                out[(True, "missing-putative-row", "")].append(s)
            else:
                g, cls = unexplained_or(causes)[:2]
                out[(g, cls, ",".join(sorted(causes)))].append(s)
        for (g, cls, why), ss in sorted(out.items(), key=lambda kv: kv[0][1]):
            text = ("hold%s this target at %d (+ reading %s)%s and no row reports it for "
                    "them, nor one at that Cluster Position naming the same variants with "
                    "the target's carriers%s"
                    % ("s" if g else "", t0.start, t0.probe,
                       " as the strongest PUTATIVE target" if d.kind == "putative" else "",
                       "" if g else " (documented: %s)" % why))
            if g:  # a failure is one line per (sample, target)
                for x in ss:
                    rep.add(True, cls, row, "%s %s" % (x, text), samples=x)
            else:
                rep.add(False, cls, row, "%d sample(s) %s: %s"
                        % (len(ss), text, ",".join(ss[:12])), samples=",".join(ss))


def documented_for(win, ctx, d, s, lids):
    """The documented verdict on a demand: 'lists' if CRISPRme's rules report any of its
    alternatives for the sample, 'maybe' if they leave the choice open."""
    seen = set()
    best = "off"
    for t, _toks in d.alts:
        t.named = [t.main] + [x for x in t.src if x is not None and x is not t.main
                              and x.kind == "snv"]
        verdict, causes = documented(win, ctx, t, s, win.lines_in(t.lo, t.hi))
        seen |= causes
        if verdict == "lists":
            return "lists", seen
        if verdict == "maybe":
            best = "maybe"
    return best, seen


def documented_carriers(win, ctx, t, carriers, lids):
    """(listed, PUTATIVE-open): the carriers of the record CRISPRme's documented rules
    list for the target, and those they resolve PUTATIVE without fixing the choice."""
    listed, maybe = set(), set()
    lids = win.lines_in(t.lo, t.hi)
    cache = {}
    for s in carriers:
        sig = win.signature(s, lids)
        v = cache.get(sig)
        if v is None:
            v = documented(win, ctx, t, s, lids)[0]
            cache[sig] = v
        if v == "lists":
            listed.add(s)
        elif v == "maybe":
            maybe.add(s)
    return listed, maybe


class Reported(object):
    """What the tables report, for completeness."""

    def __init__(self, rows, ctx, loci=None, merge=3, criterion="mm+bulges"):
        self.exact = set()
        self.site = defaultdict(list)    # (cluster, variants) -> [carrier sets]
        self.at = set()                  # (sample, cluster, record token)
        self.ctx = ctx
        self.merge = merge
        self.criterion = criterion
        self.loci = {}
        for ch, lines in (loci or {}).items():
            pos = []
            for k, (clusters, _total, _cfd, _ln) in enumerate(lines):
                for c in clusters:
                    pos.append((c, k))
            pos.sort()
            self.loci[ch] = (pos, lines)
        for r in rows:
            t = r["alt"].replace("-", "").upper()
            probe = t if r["strand"] == "+" else rc(t)
            carriers = frozenset(samples_of(r))
            var = frozenset(x for x in r["var"].split(",") if x)
            cl = cluster_of(r["start"], r["strand"], len(probe), len(r["spacer_col"]),
                            ctx.pam_5prime)
            for s in carriers:
                self.exact.add((s, r["strand"], r["start"], probe, var))
                for v in var:
                    if is_indel_token(v):
                        self.at.add((s, cl, v))
            self.site[(cl, var)].append(carriers)

    def has(self, s, strand, start, probe, var):
        return (s, strand, start, probe, var) in self.exact

    def folded(self, sample, cluster, var, exact_who, putative_who):
        """A row of the same Cluster Position naming the same variants whose carriers are
        the target's own: the merge keeps the best of those (_fold_same_carriers)."""
        allowed = set(exact_who) | set(putative_who)
        for carriers in self.site.get((cluster, var), ()):
            if sample in carriers and exact_who <= carriers and carriers <= allowed:
                return True
        return False

    def listed_at(self, sample, cluster, tok):
        return (sample, cluster, tok) in self.at

    LOCUS_REACH = 200

    def locus(self, chrom, c, total, cfd):
        """(ok, why) for a target of Cluster Position c in a run with one row per merge
        cluster: the cluster's row (one with a block within --merge bp of c; else the
        nearest row on either side, one of which is the cluster's if c is in one) has a
        fewest-mm+b block of at most the target's edits and a highest-CFD block of at
        least its CFD (printed to three decimals)."""
        pos, lines = self.loci.get(chrom, ([], []))
        i = bisect.bisect_left(pos, (c - self.merge, -1))
        near = set()
        while i < len(pos) and pos[i][0] <= c + self.merge:
            near.add(pos[i][1])
            i += 1
        if not near:
            j = bisect.bisect_left(pos, (c - self.merge, -1)) - 1
            if j >= 0 and c - pos[j][0] <= self.LOCUS_REACH:
                near.add(pos[j][1])
            k = bisect.bisect_right(pos, (c + self.merge, float("inf")))
            if k < len(pos) and pos[k][0] - c <= self.LOCUS_REACH:
                near.add(pos[k][1])
        if not near:
            return False, "no integrated row lies within %d bp of it" % self.LOCUS_REACH
        why = []
        for k in sorted(near):
            _clusters, r_total, r_cfd, ln = lines[k]
            if self.criterion == "mm+bulges" and r_total is not None and r_total > total:
                why.append("the row at line %d has a fewest-mm+b block of %d edits"
                           % (ln, r_total))
                continue
            if cfd is not None and r_cfd is not None and r_cfd < cfd - 0.0006:
                why.append("the row at line %d has a highest-CFD block of %.3f"
                           % (ln, r_cfd))
                continue
            return True, ""
        return False, "; ".join(why)


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def read_pam(path):
    """(motif, length, 5-prime?) from a CRISPRme PAM file, e.g. 'NNN...NRG 3'."""
    with open(path) as fh:
        text = fh.read().split()
    motif, k = text[0].upper(), int(text[1])
    if k < 0:
        return motif[:-k], -k, True
    return motif[-k:], k, False


def read_guides(values, pam_len, pam_5prime):
    """The guides (without PAM) from --guide values: sequences or files of one guide per
    line, as CRISPRme reads them (a guide may carry N's over the PAM)."""
    out = []
    for v in values:
        items = [v]
        if os.path.isfile(v):
            with open(v) as fh:
                items = [line.strip() for line in fh if line.strip()]
        for g in items:
            g = g.upper()
            if pam_5prime and g.startswith("N" * pam_len):
                g = g[pam_len:]
            elif not pam_5prime and g.endswith("N" * pam_len):
                g = g[:len(g) - pam_len]
            if set(g) - set("ACGTU"):
                raise ValueError("not a guide: %s" % g)
            out.append(g.replace("U", "T"))
    return out


def load_known(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "class":
                continue
            f += [""] * (7 - len(f))
            out.append({"class": f[0], "chrom": f[1], "start": f[2], "strand": f[3],
                        "variants": f[4], "samples": f[5], "reason": f[6]})
    return out


def apply_known(rep, known, chroms, modes):
    used = set()
    for r in rep.rows:
        for i, k in enumerate(known):
            if (r["class"] == k["class"] and r["chrom"] == k["chrom"]
                    and str(r["start"]) == k["start"] and r["strand"] == k["strand"]
                    and r["variants"] == k["variants"]
                    and (not k["samples"] or r["samples"] == k["samples"])):
                used.add(i)
                r["detail"] = "%s: %s | known: %s" % (r["class"], r["detail"], k["reason"])
                r["class"] = "known-defect"
                r["gating"] = "no"
                break
    for i, k in enumerate(known):
        if i in used or k["chrom"] not in chroms:
            continue
        mode = "completeness" if k["class"].startswith("missing") else "soundness"
        if mode not in modes:
            continue
        rep.rows.append({"gating": "no", "class": "known-not-seen", "table": "",
                         "line": "", "block": "", "chrom": k["chrom"], "start": k["start"],
                         "strand": k["strand"], "spacer": "", "target": "",
                         "variants": k["variants"], "samples": k["samples"],
                         "detail": "listed in --known (%s) and not found: fixed, or "
                                   "outside this run?" % k["reason"]})


G = {}  # what the work units share (set before the processes fork)


def _blockidx_of(path):
    v = Vcf(path)
    v.build_blockidx()
    return v.blockidx


def run_unit(unit):
    """One work unit: ('sound' | 'comp', chrom, lo, hi, rows). Returns the findings."""
    kind, ch, tlo, thi, rows = unit
    ctx = G["ctx"]
    rep = Report()
    refseq = G["ref"].seq(ch)
    pad = 2 * (len(ctx.guides[0]) + ctx.pam_len) + 1000
    vcfs = G["vcfs_for"][ch]
    if kind == "sound":
        spans = sorted((max(0, r["start"] - pad), r["start"] + pad) for r in rows)
        needed = []
        for a, b in spans:
            if needed and a <= needed[-1][1]:
                needed[-1] = (needed[-1][0], max(needed[-1][1], b))
            else:
                needed.append((a, b))
        win = Window(refseq, ch, max(0, tlo - pad), thi + pad, vcfs, needed=needed)
        for r in sorted(rows, key=lambda x: x["start"]):
            if any(is_indel_token(x) for x in r["var"].split(",")):
                check_row(r, win, ctx, rep)
            else:
                check_snp_row(r, win, ctx, rep)
    else:
        win = Window(refseq, ch, max(0, tlo - pad), thi + pad, vcfs)
        completeness(win, ctx, rep, G["reported"][ch], (tlo, thi), G["popmode"])
    return rep.rows


def default_cfd_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "PostProcess")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--results", required=True, help="a CRISPRme result folder")
    p.add_argument("--reference", required=True,
                   help="the reference genome: a folder of per-chromosome FASTAs or one "
                        "FASTA (any number of contigs)")
    p.add_argument("--vcf", action="append", default=[], help="a bgzipped VCF (repeatable)")
    p.add_argument("--vcf-dir", action="append", default=[], help="a folder of VCFs")
    p.add_argument("--pam", required=True, help="the PAM file the search used")
    p.add_argument("--guide", action="append", required=True,
                   help="the guide (or the guide file the search read); repeatable")
    p.add_argument("--mm", type=int, required=True)
    p.add_argument("--bDNA", dest="bdna", type=int, required=True)
    p.add_argument("--bRNA", dest="brna", type=int, required=True)
    p.add_argument("--max-total-edits", type=int, default=6)
    p.add_argument("--chrom", default=None)
    p.add_argument("--region", default=None, help="chrom:start-end (1-based, inclusive)")
    p.add_argument("--mode", default="both", choices=("soundness", "completeness", "both"))
    p.add_argument("--analysis", default="auto",
                   choices=("auto", "per-sample", "population-level"),
                   help="the run's analysis mode (auto: read from the result folder)")
    p.add_argument("--snp-rows", default="off", choices=("off", "report"),
                   help="also check SNP-only rows of a per-sample run, report-only")
    p.add_argument("--known", default=None,
                   help="a TSV of findings already reported upstream (columns: class, "
                        "chrom, start, strand, variants, samples (empty: any), reason): "
                        "a matching finding turns non-gating as known-defect, and an "
                        "entry no longer found is reported as known-not-seen")
    p.add_argument("--cfd-tables", default=None,
                   help="the folder of cfd_mismatch_scores.json and cfd_pam_scores.json "
                        "(default: the repository's PostProcess)")
    p.add_argument("--jobs", type=int, default=1, help="processes to split the work over")
    p.add_argument("--out", default=None, help="the TSV of findings")
    p.add_argument("--tile", type=int, default=50000,
                   help="how many bases of a chromosome one work unit covers")
    p.add_argument("--progress", type=int, default=0, help="a line every N work units")
    args = p.parse_args(argv)

    if args.max_total_edits is not None and args.max_total_edits < 0:
        args.max_total_edits = None
    pam_motif, pam_len, pam_5prime = read_pam(args.pam)
    try:
        guides = read_guides(args.guide, pam_len, pam_5prime)
    except ValueError as e:
        sys.stderr.write("ERROR: %s\n" % e)
        return 2
    vcf_paths = list(args.vcf)
    for d in args.vcf_dir:
        vcf_paths.extend(sorted(glob.glob(os.path.join(d, "*.vcf.gz"))))
    if not vcf_paths:
        sys.stderr.write("ERROR: no VCF given (--vcf / --vcf-dir)\n")
        return 2
    chrom = args.chrom
    lo = hi = None
    if args.region:
        m = re.match(r"^([^:]+):(\d+)-(\d+)$", args.region)
        if not m:
            sys.stderr.write("ERROR: --region wants chrom:start-end\n")
            return 2
        chrom, lo, hi = m.group(1), int(m.group(2)) - 1, int(m.group(3))
    tables = result_tables(args.results)
    if not tables:
        sys.stderr.write("ERROR: no results table under %s\n" % args.results)
        return 2
    mode_run, alt_written, merge, criterion = run_kind(args.results)
    if args.analysis != "auto":
        mode_run = args.analysis
        alt_written = False if mode_run == "population-level" else alt_written
    # one row per cluster (no alternative alignments): completeness is checked per cluster
    popmode = alt_written is False
    cfd = None
    try:
        cfd = Cfd(args.cfd_tables or default_cfd_dir())
    except (OSError, ValueError) as e:
        sys.stderr.write("check_indel_rows: no CFD tables (%s): the PUTATIVE choices CFD "
                         "decides are left open and the per-cluster check compares edits "
                         "only\n" % e)
    ctx = Ctx(args, pam_motif, pam_len, pam_5prime, guides, cfd)
    # a target of a record inside the region starts at most one target (and its bulges)
    # before the record and may hold a long REF: rows are read with that margin, so that
    # completeness sees the rows of the region's records; soundness checks the region's own
    full_len = len(ctx.full(guides[0]))
    margin = full_len + max(ctx.bdna, ctx.brna) + Vcf.LONGEST
    rows = []
    loci = {}
    bad = []
    for path, label in tables:
        rows.extend(read_rows(path, label, chrom, lo, hi, margin, pam_5prime,
                              loci=loci if label == "integrated" else None, bad=bad))
    indel_rows = [r for r in rows if any(is_indel_token(t) for t in r["var"].split(","))]
    snp_rows = [r for r in rows if r["origin"] == "alt"
                and not any(is_indel_token(t) for t in r["var"].split(","))]
    sys.stderr.write("check_indel_rows: %d indel/MNV row block(s), %d SNP-only, %s run%s, "
                     "from %s\n" % (sum(1 for r in indel_rows if r["own"]),
                                    sum(1 for r in snp_rows if r["own"]),
                                    mode_run or "unknown",
                                    " without alternative alignments" if popmode else "",
                                    ", ".join(os.path.basename(p) for p, _l in tables)))
    if mode_run == "population-level" and args.snp_rows == "report":
        sys.stderr.write("check_indel_rows: a population-level run's SNP rows are a "
                         "worst-case union; --snp-rows report is ignored\n")
        args.snp_rows = "off"
    ref = Reference(args.reference)
    ref.keep_all = args.jobs > 1  # loaded once in this process, shared with the workers
    vcfs_all = [Vcf(p) for p in vcf_paths]
    if not any(v.samples for v in vcfs_all):
        sys.stderr.write("check_indel_rows: the VCF(s) hold no genotypes (a sites-only "
                         "panel): carriers are not checked, only the rows themselves, and "
                         "completeness has no haplotype to enumerate\n")
    todo = [v for v in vcfs_all if v.tbi is None and v.bgzf]
    if todo:
        sys.stderr.write("check_indel_rows: %d VCF(s) without a .tbi: reading each once "
                         "into an in-memory block index\n" % len(todo))
        if args.jobs > 1 and len(todo) > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(min(args.jobs, len(todo))) as pool:
                for v, idx in zip(todo, pool.map(_blockidx_of, [v.path for v in todo])):
                    v.blockidx = idx
        else:
            for v in todo:
                v.build_blockidx()

    rep = Report()
    for table, line_no, nf, nh in bad:
        rep.add(True, "row-unreadable", {"table": table, "line": line_no},
                "the line has %d fields, the header %d: not read" % (nf, nh))
    check_rows = [r for r in indel_rows + (snp_rows if args.snp_rows == "report" else [])
                  if r["own"]]
    by_chrom = defaultdict(list)
    for r in check_rows:
        by_chrom[r["chrom"]].append(r)
    chroms = set(by_chrom) if args.mode in ("soundness", "both") else set()
    if args.mode in ("completeness", "both"):
        if chrom:
            chroms.add(chrom)
        else:
            # every chromosome a VCF holds (and the reference has), not only those with
            # rows: a chromosome whose indel rows were all lost is checked too
            prefix = any(r["chrom"].startswith("chr") for r in rows) or not rows
            for v in vcfs_all:
                for n in v.chroms():
                    name = n if (n.startswith("chr") or not prefix) else "chr" + n
                    if ref.has(name):
                        chroms.add(name)
            chroms |= set(by_chrom)
            if not chroms:
                sys.stderr.write("ERROR: completeness found no chromosome to check; give "
                                 "--chrom or --region\n")
                return 2
    vcfs_for = {}
    lengths = {}
    for ch in sorted(chroms):
        cov = [v for v in vcfs_all if v.covers(ch) is not False]
        vcfs_for[ch] = cov
        why = None
        if not cov:
            why = "no VCF given covers %s" % ch
        else:
            try:
                lengths[ch] = len(ref.seq(ch))
            except (KeyError, FileNotFoundError) as e:
                why = "no reference sequence for %s (%s)" % (ch, e)
        if why is None:
            continue
        n = 0
        for r in by_chrom.get(ch, ()):
            n += 1
            rep.add(True, "row-not-checked", r, why + ": this row is not checked")
        if args.mode in ("completeness", "both") and chrom == ch:
            rep.add(True, "row-not-checked", {"chrom": ch, "table": "completeness"},
                    why + ": completeness not run")
        sys.stderr.write("check_indel_rows: %s (%d row block(s) not checked)\n" % (why, n))
    if not ref.keep_all:
        ref.cache = {}
    reported = {}
    for ch in chroms:
        here = [r for r in indel_rows if r["chrom"] == ch]
        reported[ch] = Reported(here, ctx, loci={ch: loci.get(ch, [])}, merge=merge,
                                criterion=criterion)
    if popmode and criterion != "mm+bulges" and args.mode != "soundness":
        sys.stderr.write("check_indel_rows: the run sorts its fewest block by %s first: "
                         "the per-cluster check compares the CFD only\n" % criterion)

    units = []
    tile = args.tile
    for ch in sorted(chroms):
        if ch not in lengths:
            continue
        if args.mode in ("soundness", "both"):
            tiles = defaultdict(list)
            for r in by_chrom.get(ch, ()):
                tiles[r["start"] // tile].append(r)
            for t in sorted(tiles):
                units.append(("sound", ch, t * tile, (t + 1) * tile, tiles[t]))
        if args.mode in ("completeness", "both"):
            clo = lo if lo is not None else 0
            chi = min(hi, lengths[ch]) if hi is not None else lengths[ch]
            t = clo
            while t < chi:
                units.append(("comp", ch, t, min(chi, t + tile), []))
                t += tile
    G.update(ctx=ctx, ref=ref, vcfs_for=vcfs_for, reported=reported, popmode=popmode)
    done = 0
    # each unit's findings are appended to <out>.partial as it ends, so a long run that is
    # stopped keeps what it found; the sorted TSV replaces it at the end
    partial = None
    if args.out:
        partial = open(args.out + ".partial", "w", newline="")
        pw = csv.DictWriter(partial, fieldnames=COLS, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
        pw.writeheader()

    def keep(found):
        rep.rows.extend(found)
        if partial is not None:
            for r in found:
                pw.writerow(r)
            partial.flush()
    if args.jobs > 1 and len(units) > 1:
        import multiprocessing as mp
        ctx_mp = mp.get_context("fork")
        with ctx_mp.Pool(args.jobs) as pool:
            for found in pool.imap_unordered(run_unit, units, chunksize=1):
                keep(found)
                done += 1
                if args.progress and done % args.progress == 0:
                    sys.stderr.write("  %d/%d unit(s), %d finding(s)\n"
                                     % (done, len(units), len(rep.rows)))
                    sys.stderr.flush()
    else:
        for u in sorted(units, key=lambda x: (x[1], x[2])):
            keep(run_unit(u))
            done += 1
            if args.progress and done % args.progress == 0:
                sys.stderr.write("  %s %d/%d unit(s), %d finding(s)\n"
                                 % (u[1], done, len(units), len(rep.rows)))
                sys.stderr.flush()
    if args.known:
        modes = {"soundness", "completeness"} if args.mode == "both" else {args.mode}
        apply_known(rep, load_known(args.known), chroms, modes)
    rep.rows.sort(key=lambda r: (r["chrom"], str(r["start"]).zfill(12), r["class"],
                                 r["table"], str(r["line"]), r["block"]))
    if args.out:
        rep.write(args.out)
        partial.close()
        os.remove(args.out + ".partial")
    counts = rep.counts()
    print("check_indel_rows: %d finding(s), %d gating" % (len(rep.rows), rep.gating()))
    for cls in sorted(counts):
        g = sum(1 for r in rep.rows if r["class"] == cls and r["gating"] == "yes")
        print("  %-34s %7d%s" % (cls, counts[cls], "" if g else "  (non-gating)"))
    print("check_indel_rows: %s" % ("FAIL" if rep.gating() else "PASS"))
    return 1 if rep.gating() else 0


if __name__ == "__main__":
    sys.exit(main())
