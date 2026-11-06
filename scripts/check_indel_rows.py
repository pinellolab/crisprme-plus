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
and the result folder. It shares no code with CRISPRme's variant or indel machinery
(``PostProcess/indel_snp_cis.py``, ``analisi_indels_NNN.py``, ``overlay_indel_snps*``,
``resultIntegrator``, ``remove_contiguous_samples``, the CRISPRitz enricher, the
dictionaries, the indel log or the fake contigs) and imports nothing but Python's
standard library: it reads the VCF with its own BGZF + tabix reader, builds every
sample's two haplotypes itself, keeps its own map from haplotype offsets to reference
coordinates, and aligns the guide with its own exhaustive dynamic-programming aligner.
The rules CRISPRme documents (METHODS sections 4 and 8, CHANGELOG 2.8.1 Known limits)
are re-implemented here from that text, so a disagreement is either a CRISPRme defect,
a documented limit (counted and listed under its own class), or a bug in one of the two.

``crisprme.py validate-test`` compares the engine's SITES (chromosome, start, strand) on
the reference and on the SNV-only 1000 Genomes chr22 genome, so an indel row's
coordinate, its aligned strings, its REF/ALT text, its variants, its carriers and the
alternative-alignment table are never checked against an independent truth on real data.
This is that check.

TWO MODES.

  SOUNDNESS (every reported indel/MNV row is real). Each row block of the integrated
  results (highest-CFD, fewest mm+b, highest CRISPR-Bulge) and of the alternative
  alignments whose variant list names an indel or an MNV:
    * names exactly one indel-path record (an indel or an MNV) and SNPs, every one a
      record of the VCF, forward strand, on the row's chromosome;
    * its aligned spacer, gaps removed, is the guide (--guide, never the table's own
      column) followed (or, for a 5' PAM, preceded) by the PAM; its mismatches, bulge
      types and sizes and their total recount exactly, within the budget; no gap falls in
      the PAM; the PAM matches the searched motif; no IUPAC code is left;
    * its target (gaps removed, PAM included) is the reference with exactly the named
      variants applied, read from the reported start under the convention below, and
      holds a base of its indel; every named SNP changes a base of the target;
    * its REF text is the reference read from the start, gapped as the target;
    * every listed carrier really has it, and no carrier is missing (see CARRIERS).

  COMPLETENESS (no real indel off-target is missed -- the real-data brute force). For
  every indel and MNV record of the VCF in the checked region: build each carrier's
  haplotype window around it, deduplicate identical windows, enumerate every in-budget
  alignment on both strands with the PAM that holds a base the record changes (the
  junction rule), map it to the reference, and require that each (sample, target) is
  reported for that sample at that site -- as an integrated row or an alternative
  alignment -- or is merged into a better row by the documented merge key. A carrier
  CRISPRme resolves PUTATIVE must be listed at that site on a row naming the record.
  A run that writes no alternative alignments keeps one row per cluster -- a
  population-level run, or one made with --ci-cd-test, as the documented chr22 self-test
  (.alt_mode = no-alt-alignments) -- so there completeness is checked per locus: some
  integrated row must lie within one target of the site.

DEFINITIONS (fix/282-indel-alignment; METHODS section 4, "SNP+indel off-targets", and
the merge of ``remove_contiguous_samples.py``).

  Coordinate. A row starts at the 0-based reference coordinate of its target's FIRST
  ALIGNED BASE, the leftmost on either strand. A base an insertion adds has none: such a
  row starts at the reference base the added bases follow, the insertion's ANCHOR (its
  VCF POS - 1, 0-based). Up to 2.8.1 it started at that base's place on the haplotype.

  REF text. ``Aligned_protospacer+PAM_REF`` is the reference read from the start for as
  many bases as the target holds, in the target's orientation, with the target's gap
  columns put back at the same indices (letter case is not compared).

  Bulges. A gap in the aligned spacer is a DNA bulge, a gap in the aligned target an RNA
  bulge. None falls in the PAM columns; the spacer's guide end away from the PAM is never
  a gap (a leading RNA bulge, and a DNA bulge right next to the PAM, are allowed).

  Junction rule. A target holds its record when it holds a base the ALT adds or
  replaces; for a pure deletion, when it holds the bases on both sides of the deleted
  ones. A target that holds none lies in a flank: reference or SNP sequence, which the
  reference and SNP searches report; it is no indel row, and completeness does not ask
  for it.

  Readings. A carrier's haplotype holds a target at the target's place when the target's
  bases are on it from the place of its first base on, or up to the place of its last
  base, through the record (its ALT bases and the base past them). That is how CRISPRme
  reads a carrier with a second indel inside the target (a repeat unit gained or lost at
  the target's edge leaves the bases, and the carrier keeps its row).

  CRISPRme's view of a carrier slot (the rules this re-implements, METHODS section 4).
  The fake-indel window holds one indel-path record (an indel or an MNV) and the SNPs of
  the enriched genome. For a carrier of the row's record:
    * another allele at the record's own POS (an indel keeping the same anchor, or a SNP)
      on the same slot is a phasing artefact and is ignored; a variant only on bases the
      record deletes changes nothing;
    * a phased store indel overlapping the record at another position (one complex
      allele split into overlapping records) takes that slot off the target's rows, as
      does a SNP on a base an unnormalized record keeps past its POS;
    * a phased second store indel that changes the target's span keeps the slot only when
      the slot's bases still read as the target (a reading above); unreadable bases
      (an unphased call, two overlapping variants, a SNP on the second indel's REF) take
      it off;
    * a slot is EXACT -- listed for the bases its own SNP alleles give -- when the record
      and every SNP it carries in the target are phased alike; a sample that carries no
      SNP alt in the target is exact whatever its phase. Anything else (an unphased call,
      a missing allele, a haploid call against a diploid one, an unphased second indel,
      an MNV, whose phase the indel store never keeps) is PUTATIVE: listed for the
      strongest target its own alleles allow (the guide's base wherever it carries it),
      never for an allele it does not carry.
    * the enriched genome and the SNP stores keep only the last of several SNP records at
      one position (CHANGELOG 2.8.1, Known limits), so an earlier record's alt is read as
      the reference base there; MNVs are not in the phased indel store.

  CARRIERS (soundness). Each listed carrier and each carrier of the row's record is
  judged twice: on its TRUE haplotypes (every VCF record) and on CRISPRme's documented
  view above. A listed carrier whose true haplotype lacks the target is a gating failure
  unless the documented view explains it (then it is a limit, non-gating, with its
  cause); a carrier whose true haplotype holds the target with the row's variants and
  that the documented view lists is a gating failure when it is missing (a PUTATIVE
  carrier when its strongest target is determined without CFD and is the row's). The
  merge never moves a carrier off its sequence's row (same position, variants and
  aligned target merge carriers; equal carrier sets fold to the best alignment), so the
  row itself must list it.

  Merge (completeness). The pipeline clusters rows by guide and chromosome with Cluster
  Positions within --merge bp of the previous, chained. Inside a cluster, variant rows of
  one (Cluster Position, variant list) -- and for an indel/MNV row one aligned target --
  merge their carriers, and of the rows of one (Cluster Position, variant list, carrier
  set) only the best is kept (``_fold_same_carriers``). CRISPRitz's Cluster Position is
  the start moved to the PAM's fixed end (less the RNA bulges and plus the DNA bulges on
  the + strand with a 3' PAM, the - strand with a 5' one), which every alignment of one PAM
  site shares: CRISPRitz reports dozens of alignments of one site, and the fold keeps one.
  So a (sample, target) is satisfied by a row reporting exactly it with its variants, or
  by a row of the same Cluster Position naming the same variants whose carrier set is the
  target's own as CRISPRme writes it: every carrier the documented view lists for the
  target, and nobody but those and PUTATIVE ones (a carrier listed for a sequence its true
  haplotype lacks -- an earlier split-site allele CRISPRme cannot see -- counts, as it
  counts in CRISPRme's fold).

LIMITS (non-gating, counted per class, one TSV line per row block or per (sample,
target)): limit-last-alt-wins (a split multi-allelic SNP site: only the last record is
searched), limit-two-indels (a target needing a second indel-path record -- an indel or
an MNV -- on the same slot: one record per fake-indel window), limit-complex-allele (a
slot holding a record that overlaps the row's at another position), limit-putative-
carrier (a carrier CRISPRme cannot place on a slot, listed for the strongest target its
alleles allow, which its known haplotype does not have), limit-contig-end (a record
within one target of a contig's end), and in the SNP-row report limit-snp-row-indel (a
SNP row's carrier whose haplotype holds an indel inside the target). Rows on a chromosome
no VCF covers, or with no reference sequence, are reported (row-not-checked, gating):
nothing is skipped silently. A finding of the documented view and the haplotypes parting
for no documented reason is model-gap (gating: a defect of one of the two).

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
ORDER = "ACGT"


def rc(s):
    return s.translate(COMP)[::-1]


def comp(b):
    return b.translate(COMP)


def iupac_match(motif_base, base):
    return base.upper() in IUPAC.get(motif_base.upper(), "")


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
        """A record as a dict; the genotype columns stay unsplit (Window splits only the
        records a window needs)."""
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
                "tail": f[8] if len(f) > 8 else ""}


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
# records and genotypes
# --------------------------------------------------------------------------- #
def kept_bases(ref, alt):
    """The bases REF and ALT share at the start (the VCF anchor)."""
    c = 0
    while c < min(len(ref), len(alt)) and ref[c] == alt[c]:
        c += 1
    return c


def shared_tail(ref, alt):
    """The bases REF and ALT share at the end, past the shared start (an unnormalized
    record, TCA>TA)."""
    c = kept_bases(ref, alt)
    t = 0
    while t < min(len(ref), len(alt)) - c and ref[-1 - t] == alt[-1 - t]:
        t += 1
    return t


def changed(key):
    """(cs, ce, n_alt): the reference bases [cs, ce) the record deletes or replaces and
    how many ALT bases stand between its shared start and end."""
    pos1, ref, alt = key
    c = kept_bases(ref, alt)
    t = shared_tail(ref, alt)
    return pos1 - 1 + c, pos1 - 1 + len(ref) - t, len(alt) - t - c


def gone_of(key):
    """(first, last) of the reference bases the record deletes or replaces, or None."""
    cs, ce, _na = changed(key)
    return (cs, ce - 1) if ce > cs else None


def is_snp(key):
    return len(key[1]) == 1 and len(key[2]) == 1


def is_mnv(key):
    return len(key[1]) == len(key[2]) and len(key[1]) > 1


def is_store_indel(key):
    """An indel: what the phased indel genotype store keeps (not an MNV)."""
    return len(key[1]) != len(key[2])


def changes_span(key, lo, hi, gone=None):
    """True iff the record changes the reference bases lo..hi (inclusive) of a haplotype
    carrying it: a deleted or replaced base inside (not one ``gone`` already), or bases
    inserted between two bases of it."""
    pos1, ref, alt = key
    a = pos1 - 1
    c = kept_bases(ref, alt)
    t = shared_tail(ref, alt)
    if len(ref) - t > c:
        first, last = a + c, a + len(ref) - t - 1
        parts = [(first, last)] if gone is None else [
            (first, min(last, gone[0] - 1)), (max(first, gone[1] + 1), last)]
        if any(f <= l and f <= hi and l >= lo for f, l in parts):
            return True
    if len(alt) - t > c:
        return lo <= a + c - 1 and a + c <= hi
    return False


def gt_tokens(gt):
    """(tokens, kind) of a genotype: kind 'phased' ('|'), 'unphased' ('/') or
    'haploid' (one token)."""
    g = gt.split(":", 1)[0]
    if "|" in g:
        return g.split("|"), "phased"
    if "/" in g:
        return g.split("/"), "unphased"
    return [g], "haploid"


def record_bases(key, rank):
    """(cs, ce, [(base, place)]) of a record's changed part. A base standing in a REF
    base's place keeps that place, (p, 0, 0); the k-th base the ALT adds past them is
    (p, rank, k), p being the last reference base before it (an insertion's anchor)."""
    pos1, ref, alt = key
    c = kept_bases(ref, alt)
    t = shared_tail(ref, alt)
    rmid = ref[c:len(ref) - t]
    amid = alt[c:len(alt) - t]
    cs = pos1 - 1 + c
    out = []
    inplace = min(len(rmid), len(amid))
    for j in range(inplace):
        out.append((amid[j], (cs + j, 0, 0)))
    anchor = cs + len(rmid) - 1
    for j in range(len(rmid), len(amid)):
        out.append((amid[j], (anchor, rank, j - len(rmid))))
    return cs, cs + len(rmid), out


def apply_keys(refseq, keys, lo, hi, main=None):
    """(sequence, places, conflicts) of the reference [lo, hi) with ``keys`` applied.
    ``places`` (sorted) gives each base's place: (p, 0, 0) at reference position p,
    (p, rank, k) for the k-th base an ALT adds after p (rank 2 for ``main``, the row's
    record, 1 for any other). A key whose changed bases meet one already applied, or two
    insertions at one point, is left out and returned in ``conflicts``."""
    edits = []
    for key in keys:
        cs, ce, bases = record_bases(key, 2 if key == main else 1)
        if ce < lo or cs > hi:
            continue
        edits.append((cs, ce, 0 if ce == cs else 1, bases, key))
    edits.sort(key=lambda e: (e[0], e[2], e[1]))
    seq = []
    places = []
    i = lo
    conflicts = []
    last_point = None
    for cs, ce, _ins, bases, key in edits:
        if cs < i:
            if i > lo or ce <= lo:
                conflicts.append(key)
                continue
        if ce == cs:
            if last_point == cs:
                conflicts.append(key)  # two insertions at one point
                continue
            last_point = cs
        if cs > i:
            seq.append(refseq[i:cs])
            places.extend((p, 0, 0) for p in range(i, cs))
        for b, pl in bases:
            if lo <= pl[0] < hi or (pl[1] and lo - 1 <= pl[0] < hi):
                if pl[1] and pl[0] < lo:
                    continue
                seq.append(b)
                places.append(pl)
        i = max(i, ce)
    if i < hi:
        seq.append(refseq[i:hi])
        places.extend((p, 0, 0) for p in range(i, hi))
    return "".join(seq), places, conflicts


def main_places(key, places):
    """The places of a record's ALT bases among ``places`` and the place of the base past
    its REF: the core a reading must hold to be through the record."""
    pos1, ref, _alt = key
    p0 = pos1 - 1
    core = set()
    for pl in places:
        if pl[1] == 2 or (pl[1] == 0 and p0 <= pl[0] < p0 + len(ref)):
            core.add(pl)
    past = (p0 + len(ref), 0, 0)
    return core, past


def holds_main(occ, key):
    """The junction rule on a target's places (a haplotype carrying ``key``): it holds a
    base the ALT adds or replaces, or -- a pure deletion -- the bases either side."""
    cs, ce, na = changed(key)
    nr = ce - cs
    if na > 0:
        for pl in occ:
            if pl[1] == 2:
                return True
            if pl[1] == 0 and cs <= pl[0] < cs + min(nr, na):
                return True
        return False
    for x, y in zip(occ, occ[1:]):
        if x[1] == 0 and y[1] == 0 and x[0] == cs - 1 and y[0] == ce:
            return True
    return False


def reading(seq, places, n, side, first_place, last_place):
    """(begin index, bases, places) of the n bases read from ``first_place`` on (side 0)
    or up to ``last_place`` (side 1), or None past the ends."""
    if side == 0:
        b = bisect.bisect_left(places, first_place)
    else:
        b = bisect.bisect_right(places, last_place) - n
    if b < 0 or b + n > len(seq):
        return None
    return b, seq[b:b + n], places[b:b + n]


def reads_as(seq, places, want, first_place, last_place, core):
    """True iff the haplotype (seq, places) holds the bases ``want`` at the target's
    place, from its first base on or up to its last, through the record (``core``)."""
    n = len(want)
    for side in (0, 1):
        r = reading(seq, places, n, side, first_place, last_place)
        if r is not None and r[1] == want and core <= set(r[2]):
            return True
    return False


class Window(object):
    """The VCF records of one region and their genotypes.

    ``carry[key]`` -> {sample: (tokens, kind)} for the samples carrying the alt;
    ``tok[key]`` the alt's allele index; ``shadowed`` the SNP records CRISPRme's enriched
    genome and stores never see (a later SNP record at the same position of the same
    file); ``untieable`` the alts of a multi-ALT line whose allele index the genotypes
    cannot tell CRISPRme (METHODS section 4: every carrier of two alts is 1|2)."""

    HOMREF = ("0|0", "0/0", "0", "", ".", ".|.", "./.")

    def __init__(self, refseq, chrom, beg, end, vcfs, needed=None):
        self.chrom = chrom
        self.beg = max(0, beg)
        self.end = end
        self.refseq = refseq
        self.keys = []
        self.carry = {}
        self.tok = {}
        self.shadowed = set()
        self.untieable = set()
        self.samples = []
        self.needed = needed
        self._sorted = None
        seen = set()
        for fi, v in enumerate(vcfs):
            for s in v.samples:
                if s not in seen:
                    seen.add(s)
                    self.samples.append(s)
            names = v.samples
            snp_lines = defaultdict(list)  # POS -> [keys of one line] in file order
            for r in v.records(chrom, self.beg, self.end):
                if not self._wanted(r["pos1"] - 1, r["pos1"] - 1 + len(r["ref"])):
                    continue
                carriers = None
                line_keys = []
                alts = r["alts"]
                for k, alt in enumerate(alts):
                    if alt in (".", "*", "<NON_REF>") or not alt or alt.startswith("<"):
                        continue
                    key = (r["pos1"], r["ref"], alt)
                    if carriers is None:
                        carriers = self._carriers(r["tail"], names)
                    tok = str(k + 1)
                    who = {s: tp for s, tp in carriers.items() if tok in tp[0]}
                    if key not in self.carry:
                        self.keys.append(key)
                        self.carry[key] = {}
                        self.tok[key] = tok
                    self.carry[key].update(who)
                    line_keys.append(key)
                if len(r["ref"]) == 1 and line_keys and all(is_snp(k) for k in line_keys):
                    snp_lines[r["pos1"]].append(line_keys)
                if len(line_keys) > 1:
                    self._tie(line_keys)
            for _pos, lines in snp_lines.items():
                last = set(lines[-1])
                for keys in lines[:-1]:
                    self.shadowed.update(k for k in keys if k not in last)
        # reference position -> the SNP records there (visible or not)
        self.snp_at = defaultdict(list)
        for k in self.keys:
            if is_snp(k):
                self.snp_at[k[0] - 1].append(k)

    def _tie(self, keys):
        """Which alts of one multi-ALT line CRISPRme can tie to their allele index: the
        token every carrier of the alt shares, settled across the line's alts."""
        shared = {}
        for k in keys:
            common = None
            for toks, _kind in self.carry[k].values():
                t = {x for x in toks if x not in ("0", ".")}
                common = t if common is None else common & t
            shared[k] = common or set()
        changed_ = True
        while changed_:
            changed_ = False
            taken = {next(iter(c)) for c in shared.values() if len(c) == 1}
            for k, c in shared.items():
                if len(c) > 1 and c & taken:
                    shared[k] = c - taken
                    changed_ = True
        for k, c in shared.items():
            if len(c) != 1:
                self.untieable.add(k)

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

    @classmethod
    def _carriers(cls, tail, names):
        out = {}
        if not tail:
            return out
        fields = tail.split("\t")
        homref = cls.HOMREF
        for i, g in enumerate(fields):
            if i >= len(names):
                break
            if g in homref:
                continue
            toks, kind = gt_tokens(g)
            if any(t not in ("0", ".") for t in toks):
                out[names[i]] = (toks, kind)
        return out

    def is_record(self, key):
        return key in self.carry

    def local_keys(self, lo, hi):
        """The records whose REF span overlaps [lo, hi), in position order."""
        if self._sorted is None:
            self._sorted = sorted(self.keys, key=lambda k: (k[0], len(k[1]), k[2]))
            self._starts = [k[0] - 1 for k in self._sorted]
            self._widest = max([len(k[1]) for k in self._sorted] or [1])
        i = bisect.bisect_left(self._starts, lo - self._widest)
        j = bisect.bisect_left(self._starts, hi)
        return [k for k in self._sorted[i:j] if k[0] - 1 + len(k[1]) > lo and k[0] - 1 < hi]

    def norm(self, key, sample):
        """(tokens as '1' this alt / '0' another allele / '.' missing, kind) of a sample's
        call at ``key``, or None for a sample without the alt."""
        tp = self.carry[key].get(sample)
        if tp is None:
            return None
        toks, kind = tp
        t = self.tok[key]
        return ["1" if x == t else ("." if x == "." else "0") for x in toks], kind


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
    if pam_5prime:
        pam_first, pam_last = 0, pam_len - 1
    else:
        pam_first, pam_last = n - pam_len, n - 1
    out = []

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
            ok = iupac_match(spacer[i], target[j])
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


def pam_sites(motif, seq, pam_len):
    out = []
    for q in range(0, len(seq) - pam_len + 1):
        if all(iupac_match(m, b) for m, b in zip(motif, seq[q:q + pam_len])):
            out.append(q)
    return out


def scan_sequence(spacer, seq, max_mm, max_dna, max_rna, max_total, pam_len, pam_5prime):
    """Every in-budget alignment on ``seq``, both strands, anchored on the PAM. Yields
    (strand, offset of the leftmost aligned base on the + reading, aligned spacer,
    aligned target, mm, dna bulges, rna bulges)."""
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


def read_rows(path, table, chrom=None, lo=None, hi=None, starts=None):
    """Every variant (block, fields) of a results table. ``starts``: a dict that gets
    (chrom, block) -> [start] of every row, any origin (the population-level locus
    check)."""
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
            if len(row) < len(head):
                continue
            ch = row[col["Chromosome"]]
            if chrom and ch != chrom:
                continue
            for b in BLOCKS:
                try:
                    start = int(row[col["Start_coordinate_(%s)" % b]])
                except (KeyError, ValueError):
                    continue
                if lo is not None and not (lo <= start < hi):
                    continue
                if starts is not None:
                    starts.setdefault((ch, b), []).append(start)
                origin = row[col["REF/ALT_origin_(%s)" % b]]
                var = row[col["Variant_info_genome_(%s)" % b]]
                if origin != "alt" or var in ("", "NA", "n"):
                    continue
                out.append({
                    "table": table, "line": line_no, "block": b, "chrom": ch,
                    "spacer_col": row[col["Spacer+PAM"]],
                    "start": start,
                    "strand": row[col["Strand_(%s)" % b]],
                    "crrna": row[col["Aligned_spacer+PAM_(%s)" % b]],
                    "ref": row[col["Aligned_protospacer+PAM_REF_(%s)" % b]],
                    "alt": row[col["Aligned_protospacer+PAM_ALT_(%s)" % b]],
                    "pam": row[col["PAM_(%s)" % b]],
                    "mm": row[col["Mismatches_(%s)" % b]],
                    "bulges": row[col["Bulges_(%s)" % b]],
                    "total": row[col["Mismatches+bulges_(%s)" % b]],
                    "btype": row[col["Bulge_type_(%s)" % b]],
                    "origin": origin,
                    "var": var,
                    "samples": row[col["Variant_samples_(%s)" % b]],
                })
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
    """(analysis mode, alternative alignments written) of a run, from its own markers:
    'population-level' / 'per-sample' / None, and True / False / None. A
    population-level run writes no alternative alignments by design, nor does a run
    made with --ci-cd-test (the documented chr22 self-test): both keep one row per
    cluster (.alt_mode = no-alt-alignments)."""
    mode = alt = None
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
    if mode == "population-level" and alt is None:
        alt = False
    return mode, alt


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
    """What every check needs: the budget, the PAM, the guides."""

    def __init__(self, args, pam_motif, pam_len, pam_5prime, guides):
        self.args = args
        self.mm = args.mm
        self.bdna = args.bdna
        self.brna = args.brna
        self.cap = args.max_total_edits if args.max_total_edits is not None else 99
        self.pam_motif = pam_motif
        self.pam_len = pam_len
        self.pam_5prime = pam_5prime
        self.guides = guides  # [guide without PAM]

    def full(self, guide):
        return (self.pam_motif + guide) if self.pam_5prime else (guide + self.pam_motif)


# --------------------------------------------------------------------------- #
# CRISPRme's documented view of a carrier, and the truth
# --------------------------------------------------------------------------- #
class Target(object):
    """A target on a window that carries ``main`` (the row's indel-path record)."""

    __slots__ = ("probe", "strand", "start", "occ", "first", "last", "core", "span",
                 "code", "cols", "crrna", "aligned", "main", "n")

    def __init__(self, main, probe, strand, occ, crrna, aligned, win):
        self.main = main
        self.probe = probe
        self.strand = strand
        self.occ = occ
        self.n = len(probe)
        self.first = occ[0]
        self.last = occ[-1]
        self.start = occ[0][0]
        core, past = main_places(main, occ)
        self.core = core | ({past} if past in occ else set())
        p0 = main[0] - 1
        iend = p0 + len(main[1])
        reals = [pl[0] for pl in occ if pl[1] == 0 and not (p0 <= pl[0] < iend)]
        if any(pl[1] == 2 or (pl[1] == 0 and p0 <= pl[0] < iend) for pl in occ):
            reals += [p0, iend - 1]
        # a base another record adds sits between two reference bases: both are the
        # target's (so that record changes the target's span)
        for pl in occ:
            if pl[1] == 1:
                reals += [pl[0], pl[0] + 1]
        self.span = (min(reals), max(reals))
        # code columns: reference positions of the target outside the record's REF that
        # hold a SNP of the enriched genome (every visible SNP record)
        self.code = []
        for i, pl in enumerate(occ):
            if pl[1] == 0 and not (p0 <= pl[0] < iend):
                if pl[0] in win.snp_at and win.snp_at[pl[0]]:
                    self.code.append((i, pl[0]))
        self.crrna = crrna
        self.aligned = aligned
        self.cols = base_columns(aligned, strand) if aligned else None


def sample_calls(win, sample, local):
    """{key: (normalized tokens, kind)} of a sample's calls among ``local``."""
    out = {}
    for k in local:
        nk = win.norm(k, sample)
        if nk is not None:
            out[k] = nk
    return out


def on_slot(call, h, ikind, ploidy):
    """True / False: whether slot h holds the alt of ``call`` (normalized tokens, kind),
    read against the indel's own phase; None where it cannot be put on a slot."""
    toks, kind = call
    if kind == ikind and kind != "unphased" and len(toks) == ploidy:
        if toks[h] == "1":
            return True
        if toks[h] == ".":
            return None
        return False
    if all(t == "1" for t in toks):
        return True  # homozygous: on every slot
    return None if "1" in toks else False


def other_relation(main, key):
    """How a record stands to the row's record: 'apart', 'alternative' (another allele
    of its own anchor: a phasing artefact), 'deleted' (only on bases it deletes),
    'kept' (a SNP on a base an unnormalized record keeps past its POS) or 'complex'
    (overlapping it at another position)."""
    p0 = main[0] - 1
    iend = p0 + len(main[1])
    k0 = key[0] - 1
    kend = k0 + len(key[1])
    if kend <= p0 or k0 >= iend:
        return "apart"
    if key[0] == main[0] and (is_snp(key) or (kept_bases(key[1], key[2]) >= 1
                                              and kept_bases(main[1], main[2]) >= 1)):
        return "alternative"
    gone = gone_of(main)
    if is_snp(key):
        if gone and gone[0] <= k0 <= gone[1]:
            return "deleted"
        return "kept"
    if gone and not changes_span(key, k0, kend, gone):
        return "deleted"
    return "complex"


class Judge(object):
    """One target's verdict for any sample: CRISPRme's documented view (``model``) and
    the sample's true haplotypes (``truth``)."""

    def __init__(self, win, ctx, tgt):
        self.win = win
        self.ctx = ctx
        self.t = tgt
        self.main = tgt.main
        n = tgt.n
        w = max(len(self.main[1]), len(self.main[2]))
        self.lo = max(0, tgt.span[0] - 2 * n - w - 40)
        self.hi = min(len(win.refseq), tgt.span[1] + 2 * n + w + 40)
        self.local = [k for k in win.local_keys(self.lo, self.hi) if k != self.main]
        self.mgone = gone_of(self.main)
        self.rel = {k: other_relation(self.main, k) for k in self.local}
        self.cache = {}

    def signature(self, sample):
        icall = self.win.norm(self.main, sample)
        calls = tuple(sorted((k, tuple(c[0]), c[1])
                             for k, c in sample_calls(self.win, sample, self.local).items()))
        return (tuple(icall[0]), icall[1]) if icall else None, calls

    def judge(self, sample):
        sig = self.signature(sample)
        res = self.cache.get(sig)
        if res is None:
            res = self._judge(sig)
            self.cache[sig] = res
        return res

    # -- the verdict --------------------------------------------------------- #
    def _judge(self, sig):
        """{'model': ('exact', holds) | ('putative', allowed, holds or None) |
        ('none', why), 'truth': True / False / None}."""
        if sig[0] is None:
            return {"model": ("none", "no-record"), "truth": False, "tkind": "none"}
        itoks, ikind = list(sig[0][0]), sig[0][1]
        calls = {k: (list(tk), kd) for k, tk, kd in sig[1]}
        mnv = is_mnv(self.main)
        # an MNV's carriers come from the indel log, phase unknown (METHODS section 4)
        m_itoks, m_ikind = (["0", "1"], "unphased") if mnv else (itoks, ikind)
        ploidy = len(m_itoks)
        islots = [h for h, x in enumerate(m_itoks) if x == "1"]
        if mnv:
            blocked, why, unsure = set(), None, False
        else:
            blocked, why, unsure = self._breakers(calls, m_itoks, m_ikind, islots)
        mslots = [h for h in islots if h not in blocked]
        seen = self._seen(calls)
        truth, tkind = self._truth(calls, itoks, ikind)
        if not mslots:
            model = ("none", why or "no-slot")
        elif not seen and not unsure:
            # no SNP alt in the target: the reference bases, whatever the phase
            model = ("exact", self._window_holds([self.main]))
        elif (not mnv and not unsure and m_ikind != "unphased"
              and self._exact_ok(seen, mslots, m_ikind, ploidy)):
            model = ("exact", any(
                self._window_holds([self.main] + self._slot_snps(calls, h, m_ikind, ploidy))
                for h in mslots))
        else:
            model = self._putative(seen, ploidy)
        return {"model": model, "truth": truth, "tkind": tkind}

    def _breakers(self, calls, itoks, ikind, islots):
        """(blocked slots, why, unsure): the slots CRISPRme takes off this target and
        whether a call it cannot put on a slot makes the sample PUTATIVE."""
        t = self.t
        span = t.span
        ploidy = len(itoks)
        blocked = set()
        why = None
        unsure = False
        for k, call in calls.items():
            r = self.rel[k]
            toks, kind = call
            if r == "complex" and is_store_indel(k):
                if not changes_span(k, span[0], span[1], self.mgone):
                    continue
            elif r == "kept" and k not in self.win.shadowed:
                if (k[0] - 1, 0, 0) not in t.occ:
                    continue
            elif r == "apart" and is_store_indel(k):
                if not changes_span(k, span[0], span[1], self.mgone):
                    continue
            else:
                continue
            if r != "apart" and kind == "unphased":
                continue  # a multi-allelic heterozygote: one allele per haplotype
            if not (ikind != "unphased" and kind == ikind and len(toks) == ploidy):
                unsure = True
                continue
            for h in islots:
                if h in blocked or toks[h] != "1":
                    continue
                if r == "apart" and self._keeps(calls, h, ikind, ploidy):
                    continue
                blocked.add(h)
                why = why or {"complex": "complex", "kept": "kept"}.get(r, "second-indel")
        return blocked, why, unsure

    def _seen(self, calls):
        """(base index, key, tokens, kind) of every SNP alt the sample carries at the
        target's code columns, as CRISPRme's stores give them (the last record of a
        split site only; an alt of a multi-ALT line it cannot tie, as unphased)."""
        out = []
        for i, p in self.t.code:
            for k in self.win.snp_at.get(p, ()):
                if k in self.win.shadowed or k not in calls:
                    continue
                toks, kind = calls[k]
                if k in self.win.untieable:
                    toks, kind = ["1"] + ["0"] * (len(toks) - 1), "unphased"
                out.append((i, k, toks, kind))
        return out

    @staticmethod
    def _exact_ok(seen, slots, ikind, ploidy):
        for _i, _k, toks, kind in seen:
            if kind != ikind or len(toks) != ploidy:
                return False
            if any(toks[h] not in ("0", "1") for h in slots):
                return False
        for h in slots:
            cols = [i for i, _k, toks, _kind in seen if toks[h] == "1"]
            if len(cols) != len(set(cols)):
                return False  # two alts at one position on one slot
        return True

    def _slot_snps(self, calls, h, ikind, ploidy):
        out = []
        for k, call in calls.items():
            if is_snp(k) and k not in self.win.shadowed and self.rel[k] == "apart" \
                    and on_slot(call, h, ikind, ploidy) is True:
                out.append(k)
        return out

    def _window_holds(self, keys):
        """Whether CRISPRme's fake-indel window, written with ``keys`` (the record and a
        slot's SNP alleles), holds the target at its place."""
        seq, pl, _c = apply_keys(self.win.refseq, keys, self.lo, self.hi, self.main)
        return reads_as(seq, pl, self.t.probe, self.t.first, self.t.last, self.t.core)

    def _putative(self, seen, ploidy):
        """CRISPRme's PUTATIVE resolution: whether the target's bases at the code columns
        are among the sample's own alleles there, and -- when the choice needs no CFD --
        whether they are its strongest target."""
        t = self.t
        allowed = True
        determined = True
        holds = True
        for i, p in t.code:
            ref = self.win.refseq[p]
            carried = set()
            copies = 0
            pl = ploidy
            unknown = False
            for ii, k, toks, _kind in seen:
                if ii != i:
                    continue
                carried.add(k[2])
                copies += sum(1 for x in toks if x == "1")
                pl = max(pl, len(toks))
                unknown = unknown or any(x not in ("0", "1") for x in toks)
            if not carried or unknown or copies < pl:
                carried.add(ref)
            row_b = t.probe[i]
            if row_b not in carried:
                allowed = False
            pick = self._strongest(i, carried)
            if pick is None:
                determined = False
            elif pick == "" or pick != row_b:
                holds = False
        return ("putative", allowed, (holds and allowed) if determined else None)

    def _strongest(self, i, cands):
        """The base METHODS section 4 gives a PUTATIVE carrier at forward base i: the
        guide's where it carries it, the only candidate, or None when the choice needs
        the CFD; '' when no candidate satisfies the PAM there."""
        t = self.t
        if t.cols is None:
            return None
        col = t.cols[i]
        tb = (lambda b: b) if t.strand == "+" else comp
        cands_t = sorted({tb(b) for b in cands})
        n = len(t.crrna)
        p5 = self.ctx.pam_5prime
        first_pam = 0 if p5 else n - self.ctx.pam_len
        if first_pam <= col < first_pam + self.ctx.pam_len:
            m = self.ctx.pam_motif[col - first_pam]
            ok = [b for b in cands_t if iupac_match(m, b)]
            if not ok:
                return ""
            return tb(ok[0]) if len(ok) == 1 else None
        g = t.crrna[col].upper()
        if g != "-" and g in cands_t:
            return tb(g)
        if len(cands_t) == 1:
            return tb(cands_t[0])
        return None

    def _keeps(self, calls, h, ikind, ploidy):
        """CRISPRme's reading of a slot whose second store indel changes the target's
        span: the slot's bases (the record, the slot's store indels it can reach, its own
        SNP alleles) read as the window's target with the same alleles, from the first
        base on or up to the last, through the record (METHODS section 4). Unreadable
        bases (an unknown call where it reads, overlapping variants, a SNP alt on a base
        of a second indel's REF) do not keep it."""
        t = self.t
        win = self.win
        span = t.span
        mine = [k for k in calls if self.rel.get(k) == "apart" and is_store_indel(k)
                and on_slot(calls[k], h, ikind, ploidy) is True
                and changes_span(k, span[0], span[1], self.mgone)]
        margin = (t.n + max(len(self.main[1]), len(self.main[2]))
                  + sum(abs(len(k[1]) - len(k[2])) for k in mine))
        flo, fhi = max(0, span[0] - margin), min(len(win.refseq), span[1] + margin + 1)
        near = [k for k in win.local_keys(flo, fhi)
                if k != self.main and is_store_indel(k)
                and other_relation(self.main, k) == "apart"]
        snp_on, snp_unknown = [], set()
        for k in win.local_keys(flo, fhi):
            if not is_snp(k) or k in win.shadowed or k not in calls:
                continue
            if other_relation(self.main, k) != "apart":
                continue
            v = on_slot(calls[k], h, ikind, ploidy)
            if v is True:
                snp_on.append(k)
            elif v is None:
                snp_unknown.add(k[0] - 1)
        wkeys = [self.main] + snp_on
        wseq, wpl, _c = apply_keys(win.refseq, wkeys, flo, fhi, self.main)
        a = bisect.bisect_left(wpl, t.first)
        want = wseq[a:a + t.n]
        if len(want) != t.n:
            return False
        for side in (0, 1):
            held = set(mine)
            bad = False
            for k in near:
                if k in held or k not in calls:
                    continue
                reach = ((k[0] - 2 + len(k[1]) >= span[0]) if side == 0
                         else (k[0] - 1 <= span[1]))
                if not reach:
                    continue
                v = on_slot(calls[k], h, ikind, ploidy)
                if v is None:
                    bad = True
                    break
                if v:
                    held.add(k)
            if bad:
                continue
            on_ref = {sk[0] - 1 for sk in snp_on} | snp_unknown
            if any(k[0] - 1 <= q < k[0] - 1 + len(k[1]) for k in held for q in on_ref):
                continue  # a SNP on a base of a second indel's REF: two records of one base
            refs = sorted((k[0] - 1, k[0] - 1 + len(k[1])) for k in [self.main] + list(held))
            if any(b[0] < a_[1] for a_, b in zip(refs, refs[1:])):
                continue  # two REFs share a base
            mseq, mpl, conf = apply_keys(win.refseq, wkeys + list(held), flo, fhi,
                                         self.main)
            if conf:
                continue
            r = reading(mseq, mpl, t.n, side, t.first, t.last)
            if r is not None and r[1] == want and t.core <= set(r[2]):
                return True
        return False

    def _truth(self, calls, itoks, ikind):
        """Whether the sample's TRUE haplotypes hold the target at its place: True,
        False, or None when a slot holds records overlapping the row's record (no one
        sequence). A sample whose calls on the target cannot be put on a slot is judged
        on its alleles instead (the target's bases at every SNP are alleles it carries)."""
        t = self.t
        ploidy = len(itoks)
        islots = [h for h, x in enumerate(itoks) if x == "1"]
        touching = [k for k in calls
                    if k[0] - 1 <= t.span[1] + 1 and k[0] - 1 + len(k[1]) >= t.span[0] - 1]
        clean = ikind != "unphased" and all(
            calls[k][1] == ikind and len(calls[k][0]) == ploidy
            and all(x in ("0", "1") for x in calls[k][0]) for k in touching)
        if not clean or not islots:
            return self._alleles(calls, itoks), "alleles"
        anyc = False
        for h in islots:
            keys = [self.main]
            complex_ = False
            for k, (toks, kind) in calls.items():
                if kind != ikind or len(toks) != ploidy or toks[h] != "1":
                    continue
                r = self.rel[k]
                if r in ("alternative", "deleted"):
                    continue
                if r == "complex":
                    if changes_span(k, t.span[0], t.span[1], self.mgone):
                        complex_ = True
                    continue
                keys.append(k)
            if complex_:
                anyc = True
                continue
            seq, pl, _conf = apply_keys(self.win.refseq, keys, self.lo, self.hi, self.main)
            if reads_as(seq, pl, t.probe, t.first, t.last, t.core):
                return True, "haplotype"
        return (None if anyc else False), "haplotype"

    def _alleles(self, calls, itoks):
        """The target's base at every SNP position of it is an allele the sample carries
        (any record, the reference unless it is homozygous for alts there)."""
        t = self.t
        for i, p in t.code:
            row_b = t.probe[i]
            carried = set()
            copies = 0
            pl = len(itoks)
            for k in self.win.snp_at.get(p, ()):
                if k not in calls:
                    continue
                toks, _kind = calls[k]
                carried.add(k[2])
                copies += sum(1 for x in toks if x == "1")
                pl = max(pl, len(toks))
            if copies < pl:
                carried.add(self.win.refseq[p])
            if row_b not in carried:
                return False
        return True


def explain(judge, sample):
    """Why a sample's true haplotypes and CRISPRme's documented view part: the records on
    its slots that change the target and that CRISPRme does not read ('last-alt-wins':
    an earlier SNP record of a split site; 'two-indels': a second indel-path record, an
    MNV included; 'complex-allele': a record overlapping the row's)."""
    win, t = judge.win, judge.t
    if win.norm(judge.main, sample) is None:
        return set()
    causes = set()
    for k, (toks, _kind) in sample_calls(win, sample, judge.local).items():
        if "1" not in toks:
            continue
        rel = judge.rel.get(k, other_relation(judge.main, k))
        if rel in ("alternative", "deleted"):
            continue
        if not changes_span(k, t.span[0], t.span[1], judge.mgone):
            continue
        if k in win.shadowed:
            causes.add("last-alt-wins")
        elif rel == "complex" or rel == "kept":
            causes.add("complex-allele")
        elif not is_snp(k):
            causes.add("two-indels")
    return causes


def row_target(win, ctx, row, named, main):
    """The row's target on the reference with exactly its named variants: (Target,
    problem class, detail)."""
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    width = max(len(main[1]), len(main[2])) if main else 1
    anchor = (main[0] - 1) if main else row["start"]
    lo = max(0, min(row["start"], anchor) - n - width - 60)
    hi = min(len(win.refseq), max(row["start"], anchor + width) + 2 * n + width + 60)
    seq, places, conflicts = apply_keys(win.refseq, named, lo, hi, main)
    if conflicts:
        return None, "variants-overlap", (
            "the named variants cannot all be on one haplotype: %s overlaps another"
            % ",".join("%d_%s_%s" % k for k in conflicts))
    hits = []
    k = seq.find(probe)
    while k >= 0:
        occ = places[k:k + n]
        hits.append((occ[0][0], holds_main(occ, main) if main else True, occ))
        k = seq.find(probe, k + 1)
    for st, held, occ in hits:
        if st == row["start"] and held:
            if main is None:
                return occ, None, None
            return Target(main, probe, row["strand"], occ, row["crrna"], row["alt"],
                          win), None, None
    through = sorted({st for st, held, _o in hits if held})
    if through:
        return None, "start-not-on-variant-sequence", (
            "the target is the reference with its named variants%s starting at %s, not at "
            "%d" % (", through its record," if main else "",
                    ",".join(map(str, through[:6])), row["start"]))
    if any(st == row["start"] for st, _h, _o in hits):
        return None, "target-holds-no-indel-base", (
            "the target at %d holds no base its record %d_%s_%s changes (a flank target)"
            % (row["start"], main[0], main[1], main[2]))
    return None, "target-not-on-variant-sequence", (
        "the target is not the reference with its named variants near %d: a base no "
        "named variant explains" % row["start"])


def alignment_checks(row, ctx, rep, gating):
    """The checks of a row's own strings: its guide, its aligned spacer, the recount, the
    bulge type, the budget, the PAM, the REF text. False when the row cannot be read
    further."""
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
        return
    text = list(span if row["strand"] == "+" else rc(span))
    for i, c in enumerate(row["alt"]):
        if c == "-":
            text.insert(i, "-")
    if row["ref"].upper() != "".join(text).upper():
        rep.add(gating, "ref-text-mismatch", row,
                "REF %s, the reference at %d gapped as the target is %s"
                % (row["ref"], row["start"], "".join(text)))


def named_records(row, win, rep, gating):
    """The row's variants as VCF keys, or None after a finding."""
    named = []
    for t in [x for x in row["var"].split(",") if x and x != "NA"]:
        try:
            c, p, r, a = parse_token(t)
        except (ValueError, IndexError):
            rep.add(gating, "variant-unparseable", row, t)
            return None
        key = (p, r, a)
        if c != row["chrom"] or not win.is_record(key):
            rep.add(gating, "variant-not-in-vcf", row,
                    "%s is not a record of the VCF (forward strand)" % t)
            return None
        if key in win.shadowed:
            rep.add(gating, "variant-shadowed", row,
                    "%s is an earlier SNP record of a split site, which CRISPRme's stores "
                    "never hold" % t)
        named.append(key)
    return named


CARRIER_TEXT = {
    "carrier-without-indel": "listed without carrying %s",
    "carrier-without-target": "listed, and no haplotype of theirs holds the target at its "
                              "place with %s",
    "carrier-without-variants": "listed, and they do not carry every named variant on the "
                                "slot of %s",
    "variant-missing": "listed, and their slot of %s holds a variant inside the target "
                       "the row does not name",
    "carrier-without-allele": "listed (PUTATIVE) for an allele they do not carry, with %s",
    "carrier-off-by-rule": "listed although each of their slots of %s holds a record "
                           "overlapping it (CRISPRme takes such a slot off)",
    "carrier-missing": "hold the target at its place with %s, are listed by CRISPRme's "
                       "documented view, and are not on the row",
    "mnv-not-read-on-carrier": "listed with an MNV on their slot that changes the target "
                               "(MNVs are not in the phased indel store), with %s",
    "limit-last-alt-wins": "differ from the row only by an earlier SNP record of a split "
                           "site (CRISPRme reads the last record), with %s",
    "limit-two-indels": "hold the target only with a second indel-path record on the "
                        "slot of %s",
    "limit-complex-allele": "hold records overlapping %s on its slot (taken off by rule)",
    "limit-putative-carrier": "are PUTATIVE carriers of %s (phase unknown to CRISPRme), "
                              "listed for the strongest target their alleles allow",
}


def check_row(row, win, ctx, rep, gating=True):
    """Every soundness check of one indel/MNV row block."""
    if not alignment_checks(row, ctx, rep, gating):
        return
    ref_text_check(row, win, rep, gating)
    named = named_records(row, win, rep, gating)
    if named is None:
        return
    mains = [k for k in named if not is_snp(k)]
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
    tgt, cls, detail = row_target(win, ctx, row, named, main)
    if tgt is None:
        rep.add(gating, cls, row, detail)
        return
    for key in named:
        if key != main and (key[0] - 1, 0, 0) not in tgt.occ:
            rep.add(gating, "variant-outside-target", row,
                    "%s_%d_%s_%s changes no base of the target (%d..%d)"
                    % (row["chrom"], key[0], key[1], key[2], tgt.span[0], tgt.span[1]))
    if not win.samples:
        return  # a sites-only VCF: no genotype to judge a carrier by (noted once by main)
    judge = Judge(win, ctx, tgt)
    listed = set(samples)
    pool = set(win.carry.get(main, ())) | listed
    out = defaultdict(list)
    for s in sorted(pool):
        cls = carrier_verdict(judge, s, s in listed, named, main)
        if cls is not None:
            out[cls].append(s)
    tok = "%s_%d_%s_%s" % (row["chrom"], main[0], main[1], main[2])
    for (g, cls), who in sorted(out.items(), key=lambda kv: kv[0][1]):
        detail = "%d sample(s) %s: %s" % (len(who), CARRIER_TEXT.get(cls, "%s") % tok,
                                          ",".join(who[:12]))
        rep.add(g and gating, cls, row, detail, samples=",".join(who))


def carrier_verdict(judge, s, is_listed, named, main):
    """(gating, class) for one sample of a row, or None. The true haplotypes decide;
    CRISPRme's documented view only turns a disagreement into a limit."""
    v = judge.judge(s)
    model, truth = v["model"], v["truth"]
    if model == ("none", "no-record"):
        return (True, "carrier-without-indel") if is_listed else None
    if model[0] == "exact":
        m_holds = model[1]
    elif model[0] == "putative":
        m_holds = model[2]
    else:
        m_holds = False
    if is_listed:
        if truth:
            return None
        if truth is None:
            if model[0] == "none" and model[1] == "complex":
                return (True, "carrier-off-by-rule")
            return None  # a slot of overlapping records: no one sequence to hold
        causes = explain(judge, s)
        if model[0] == "putative" and model[1]:
            if "last-alt-wins" in causes:
                return (False, "limit-last-alt-wins")
            return (False, "limit-putative-carrier")
        if model[0] == "putative":
            return (True, "carrier-without-allele")
        if m_holds:
            if "last-alt-wins" in causes:
                return (False, "limit-last-alt-wins")
            if "two-indels" in causes and any(
                    is_mnv(k) and s in judge.win.carry.get(k, {}) for k in judge.local):
                return (True, "mnv-not-read-on-carrier")
            return (True, "carrier-without-target")
        return (True, diagnose(judge, s, named, main))
    if truth is None and model[0] == "none" and model[1] in ("complex", "kept"):
        return (False, "limit-complex-allele")  # no one sequence; CRISPRme takes it off
    if not truth:
        return None
    if m_holds:
        return (True, "carrier-missing")
    if m_holds is None:
        return None  # PUTATIVE, and its strongest target needs the CFD to choose
    if model[0] == "putative" and v["tkind"] == "alleles":
        return None  # phase unknown: listed for its strongest target, as documented
    causes = explain(judge, s)
    if model[0] == "none" and model[1] in ("complex", "kept"):
        return (False, "limit-complex-allele")
    if model[0] == "none" and model[1] == "second-indel":
        return (False, "limit-two-indels")
    if "last-alt-wins" in causes:
        return (False, "limit-last-alt-wins")
    if "two-indels" in causes:
        return (False, "limit-two-indels")
    if model[0] == "putative":
        return (False, "limit-putative-carrier")
    return (True, "carrier-missing")


def diagnose(judge, sample, named, main):
    """The class of a listed carrier that neither its haplotypes nor CRISPRme's view
    give the target."""
    win = judge.win
    calls = sample_calls(win, sample, judge.local)
    icall = win.norm(main, sample)
    if icall is None:
        return "carrier-without-indel"
    itoks, ikind = icall
    if ikind == "unphased":
        return "carrier-without-allele"
    snps = [k for k in named if k != main]
    for h, x in enumerate(itoks):
        if x != "1":
            continue
        if all(k in calls and calls[k][1] == ikind and h < len(calls[k][0])
               and calls[k][0][h] == "1" for k in snps):
            return "variant-missing" if snps or carries_more(judge, sample, h) else \
                "carrier-without-target"
    return "carrier-without-variants" if snps else "carrier-without-target"


def carries_more(judge, sample, h):
    """Whether the sample's slot h holds another variant inside the target."""
    t = judge.t
    for k, (toks, _kind) in sample_calls(judge.win, sample, judge.local).items():
        if h < len(toks) and toks[h] == "1" and judge.rel.get(k) == "apart" \
                and changes_span(k, t.span[0], t.span[1], judge.mgone):
            return True
    return False


def check_snp_row(row, win, ctx, rep):
    """The same checks on a SNP-only row of a per-sample run, report-only (CHANGELOG
    2.8.1 Known limits: a SNP row lists carriers from SNP genotypes only)."""
    if not alignment_checks(row, ctx, rep, False):
        return
    ref_text_check(row, win, rep, False)
    named = named_records(row, win, rep, False)
    if named is None:
        return
    occ, cls, detail = row_target(win, ctx, row, named, None)
    if occ is None:
        rep.add(False, cls, row, detail)
        return
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    lo = max(0, occ[0][0] - 2 * len(probe) - 60)
    hi = min(len(win.refseq), occ[-1][0] + 2 * len(probe) + 60)
    local = win.local_keys(lo, hi)
    span = (min(pl[0] for pl in occ), max(pl[0] for pl in occ))
    out = defaultdict(list)
    for s in samples_of(row):
        calls = sample_calls(win, s, local)
        kinds = {c[1] for c in calls.values()}
        if not calls or kinds - {"phased", "haploid"} or len(kinds) > 1:
            continue  # no slot to read: CRISPRme's PUTATIVE case
        ploidy = max(len(c[0]) for c in calls.values())
        ok = False
        why = None
        for h in range(ploidy):
            keys = [k for k, (toks, _kd) in calls.items() if h < len(toks) and toks[h] == "1"]
            if not all(k in keys for k in named):
                continue
            seq, pl, _c = apply_keys(win.refseq, keys, lo, hi)
            if reads_as(seq, pl, probe, occ[0], occ[-1], set()):
                ok = True
                break
            inside = [k for k in keys if k not in named
                      and changes_span(k, span[0], span[1])]
            if any(not is_snp(k) for k in inside):
                why = why or "limit-snp-row-indel"
            elif any(k in win.shadowed for k in inside):
                why = why or "limit-last-alt-wins"
        if not ok:
            out[why or "carrier-without-target"].append(s)
    text = {"limit-snp-row-indel": "hold an indel (or MNV) inside the target on the slot "
                                   "of the row's SNPs: a SNP row's carriers come from SNP "
                                   "genotypes only (CHANGELOG 2.8.1, Known limits)",
            "limit-last-alt-wins": "hold an earlier SNP record of a split site inside the "
                                   "target (CRISPRme reads the last record)",
            "carrier-without-target": "hold no haplotype with this target and every named "
                                      "variant"}
    for cls, who in sorted(out.items()):
        rep.add(False, cls, row, "%d listed sample(s) %s: %s"
                % (len(who), text[cls], ",".join(who[:12])), samples=",".join(who))


# --------------------------------------------------------------------------- #
# completeness
# --------------------------------------------------------------------------- #
def completeness(win, ctx, rep, reported, region, popmode):
    """The real-data brute force over the indel and MNV records of ``region``."""
    lo_r, hi_r = region
    for guide in ctx.guides:
        full = ctx.full(guide)
        reach = len(full) + max(ctx.bdna, ctx.brna) + 2
        for key in list(win.keys):
            if is_snp(key) or not (lo_r <= key[0] - 1 < hi_r):
                continue
            carriers = win.carry[key]
            if not carriers:
                continue
            lo = key[0] - 1 - reach - len(key[1])
            hi = key[0] - 1 + len(key[1]) + reach
            if lo < 0 or hi > len(win.refseq):
                rep.add(False, "limit-contig-end", {
                    "table": "completeness", "chrom": win.chrom, "start": key[0] - 1,
                    "var": "%s_%d_%s_%s" % (win.chrom, key[0], key[1], key[2])},
                    "the record is within one target of the contig's end")
                continue
            comp_record(win, ctx, rep, reported, key, guide, full, lo, hi, popmode)


def comp_record(win, ctx, rep, reported, key, guide, full, lo, hi, popmode):
    chrom = win.chrom
    carriers = win.carry[key]
    local = [k for k in win.local_keys(lo, hi) if k != key]
    sig = defaultdict(list)
    for k in local:
        for s, tp in win.carry[k].items():
            if s in carriers:
                sig[s].append((k,) + (tuple(tp[0]), tp[1]))
    groups = defaultdict(list)
    for s in carriers:
        groups[(tuple(carriers[s][0]), carriers[s][1], tuple(sig.get(s, ())))].append(s)
    complex_who = set()          # carriers whose slot of the record holds an overlapping one
    exact = defaultdict(set)     # (strand, start, probe, variants) -> samples
    exact_t = {}                 # the same key -> its Target (the window's occurrence)
    putative = defaultdict(set)  # (strand, start) -> samples
    limits = []                  # (class, sample, start, strand, detail)
    seen_t = {}
    tok = "%s_%d_%s_%s" % (chrom, key[0], key[1], key[2])
    for (itoks, ikind, recs), members in groups.items():
        rep_s = members[0]
        for seq, places, label, wseq, wpl in haplotypes(win, key, rep_s, lo, hi):
            if label == "complex":
                complex_who.update(members)
                continue
            hits = seen_t.get(seq)
            if hits is None:
                hits = list(scan_sequence(full, seq, ctx.mm, ctx.bdna, ctx.brna, ctx.cap,
                                          ctx.pam_len, ctx.pam_5prime))
                seen_t[seq] = hits
            for strand, off, sa, ta, mm, dn, rn in hits:
                width = len(ta.replace("-", ""))
                occ = places[off:off + width]
                if len(occ) < width or not holds_main(occ, key):
                    continue
                probe = seq[off:off + width]
                t = Target(key, probe, strand, occ, sa, ta, win)
                # CRISPRme finds the target on its window (the record and the slot's own
                # SNP alleles): read it there from the target's first base on or up to
                # its last, through the record, and judge that occurrence
                if wseq is not None and wseq != seq:
                    for side in (0, 1):
                        r = reading(wseq, wpl, width, side, t.first, t.last)
                        if r is not None and r[1] == probe and t.core <= set(r[2]):
                            t = Target(key, probe, strand, r[2], sa, ta, win)
                            break
                judge = Judge(win, ctx, t)
                for s in members:
                    v = judge.judge(s)
                    model = v["model"]
                    if label == "slot" and v["truth"] is False:
                        continue  # not this sample's sequence (another slot's)
                    if model[0] == "none":
                        why = model[1]
                        if why in ("complex", "kept"):
                            limits.append(("limit-complex-allele", s, t, sa, ta, ""))
                        elif why == "second-indel":
                            limits.append(("limit-two-indels", s, t, sa, ta,
                                           "a second indel on the slot changes the "
                                           "window's bases"))
                        continue
                    if model[0] == "exact":
                        if model[1]:
                            dk = (strand, t.start, probe, variants_of(win, t, s, key))
                            exact[dk].add(s)
                            exact_t.setdefault(dk, t)
                            continue
                        causes = explain(judge, s)
                        if "two-indels" in causes:
                            limits.append(("limit-two-indels", s, t, sa, ta, ""))
                        elif "last-alt-wins" in causes:
                            limits.append(("limit-last-alt-wins", s, t, sa, ta, ""))
                        elif "complex-allele" in causes:
                            limits.append(("limit-complex-allele", s, t, sa, ta, ""))
                        else:
                            limits.append(("model-gap", s, t, sa, ta, ""))
                        continue
                    putative[(strand, t.start, len(probe))].add(s)
    # resolve the demands against the tables
    full_len = len(full)
    fold_sets = {}
    for (strand, start, probe, var), who in sorted(exact.items()):
        cl = cluster_of(start, strand, len(probe), full_len, ctx.pam_5prime)
        put_here = set()
        for (pst, pstart, plen), pw in putative.items():
            if cluster_of(pstart, pst, plen, full_len, ctx.pam_5prime) == cl:
                put_here |= pw
        for s in sorted(who):
            if popmode:
                if not reported.locus(start):
                    rep.add(True, "missing-locus", {
                        "table": "completeness", "chrom": chrom, "start": start,
                        "strand": strand, "alt": probe, "var": ",".join(sorted(var))},
                        "no integrated row near %d (population-level: one row per "
                        "cluster)" % start, samples=s)
                break
            if (s, strand, start, probe, var) in reported.exact:
                continue
            # the merge folds rows of one Cluster Position, variant list and carrier set:
            # CRISPRme's carrier set, so it is computed by the documented view over every
            # carrier of the record (a carrier it lists for a sequence its haplotype lacks,
            # an earlier split-site allele unseen, belongs to the row as CRISPRme writes it)
            if fold_sets.get((strand, start, probe, var)) is None:
                fold_sets[(strand, start, probe, var)] = model_carriers(
                    win, ctx, exact_t[(strand, start, probe, var)], carriers)
            cm, pm = fold_sets[(strand, start, probe, var)]
            if reported.folded(s, cl, var, cm | set(who), pm | put_here):
                continue
            rep.add(True, "missing-row", {
                "table": "completeness", "chrom": chrom, "start": start, "strand": strand,
                "alt": probe, "var": ",".join(sorted(var))},
                "%s has this target at %d (+ reading %s) and no row reports it, nor one "
                "at that start naming the same variants with the target's carriers"
                % (s, start, probe), samples=s)
    for (strand, start, plen), who in sorted(putative.items()):
        cl = cluster_of(start, strand, plen, full_len, ctx.pam_5prime)
        for s in sorted(who):
            if popmode:
                if not reported.locus(start):
                    rep.add(True, "missing-locus", {
                        "table": "completeness", "chrom": chrom, "start": start,
                        "strand": strand, "var": tok},
                        "no integrated row near %d" % start, samples=s)
                break
            if reported.listed_at(s, cl, tok):
                continue
            rep.add(True, "missing-putative-row", {
                "table": "completeness", "chrom": chrom, "start": start, "strand": strand,
                "var": tok},
                "%s (PUTATIVE: its phase is unknown there) holds an in-budget target at %d "
                "and is on no row there naming %s" % (s, start, tok), samples=s)
    # the limits, one line per target (its carriers listed), and one per record for the
    # carriers whose slot holds a record overlapping it
    grouped = {}
    for cls, s, t, sa, ta, why in limits:
        k = (cls, t.strand, t.start, t.probe)
        if k not in grouped:
            grouped[k] = [set(), sa, ta, why]
        grouped[k][0].add(s)
    text = {"limit-two-indels": "the target needs a second indel-path record on the slot "
                                "(one record per fake-indel window)",
            "limit-last-alt-wins": "the target needs an earlier SNP record of a split site "
                                   "(only the last is searched)",
            "limit-complex-allele": "the slot holds a record overlapping %s" % tok,
            "model-gap": "the documented view differs from the haplotype for no documented "
                         "reason (a checker defect?)"}
    for (cls, strand, start, probe), (who, sa, ta, why) in sorted(grouped.items()):
        who = sorted(who)
        rep.add(cls == "model-gap", cls, {
            "table": "completeness", "chrom": chrom, "start": start, "strand": strand,
            "crrna": sa, "alt": ta, "var": tok},
            "%d sample(s): %s: %s" % (len(who), why or text[cls], ",".join(who[:12])),
            samples=",".join(who))
    if complex_who:
        who = sorted(complex_who)
        rep.add(False, "limit-complex-allele", {
            "table": "completeness", "chrom": chrom, "start": key[0] - 1, "var": tok},
            "%d sample(s) hold %s on a slot with a record overlapping it at another position "
            "(one complex allele split in two): no one sequence there, and CRISPRme takes "
            "the slot off the record's targets: %s" % (len(who), tok, ",".join(who[:12])),
            samples=",".join(who))


def model_carriers(win, ctx, t, carriers):
    """(exact, putative): the record's carriers CRISPRme's documented view lists for the
    target ``t`` -- exactly, or as PUTATIVE ones whose alleles allow it."""
    judge = Judge(win, ctx, t)
    exact, putative = set(), set()
    for s in carriers:
        model = judge.judge(s)["model"]
        if model == ("exact", True):
            exact.add(s)
        elif model[0] == "putative" and model[1] and model[2] is not False:
            putative.add(s)
    return exact, putative


def variants_of(win, t, sample, key):
    """The variant list CRISPRme writes for a sample's exact slot: its record and the
    SNP alts of the slot at the target's code columns (forward)."""
    out = {"%s_%d_%s_%s" % (win.chrom, key[0], key[1], key[2])}
    for i, p in t.code:
        b = t.probe[i]
        if b != win.refseq[p]:
            for k in win.snp_at.get(p, ()):
                if k[2] == b and k not in win.shadowed:
                    out.add("%s_%d_%s_%s" % (win.chrom, k[0], k[1], k[2]))
    return frozenset(out)


def haplotypes(win, key, sample, lo, hi):
    """The sequences to enumerate for one carrier of ``key``, as (seq, places, label,
    window seq, window places): its true slots when its calls there are phased ('slot',
    with CRISPRme's window for that slot: the record and the slot's own SNP alleles),
    else the record with its homozygous alleles and with every SNP alt it carries
    ('putative': the sites CRISPRme lists it at; its window is the sequence itself)."""
    icall = win.norm(key, sample)
    if icall is None:
        return
    itoks, ikind = icall
    local = [k for k in win.local_keys(lo, hi) if k != key]
    calls = sample_calls(win, sample, local)
    ploidy = len(itoks)
    clean = ikind != "unphased" and all(
        c[1] == ikind and len(c[0]) == ploidy and all(x in ("0", "1") for x in c[0])
        for c in calls.values())
    if clean and not is_mnv(key):
        for h, x in enumerate(itoks):
            if x != "1":
                continue
            keys = [key]
            wkeys = [key]
            skip = False
            for k, (toks, kind) in calls.items():
                if toks[h] != "1":
                    continue
                r = other_relation(key, k)
                if r in ("alternative", "deleted"):
                    continue
                if r == "complex":
                    skip = True
                    break
                keys.append(k)
                if r == "apart" and is_snp(k) and k not in win.shadowed:
                    wkeys.append(k)
            if skip:
                yield None, None, "complex", None, None  # no one sequence: listed
                continue
            seq, places, _conf = apply_keys(win.refseq, keys, lo, hi, key)
            wseq, wpl, _c = apply_keys(win.refseq, wkeys, lo, hi, key)
            yield seq, places, "slot", wseq, wpl
        return
    homo = [key] + [k for k, (toks, kind) in calls.items()
                    if all(x == "1" for x in toks) and other_relation(key, k) == "apart"
                    and is_snp(k)]
    allsnp = [key] + [k for k, (toks, kind) in calls.items()
                      if is_snp(k) and "1" in toks and other_relation(key, k) == "apart"]
    seen = set()
    for keys in (homo, allsnp):
        seq, places, _conf = apply_keys(win.refseq, keys, lo, hi, key)
        if seq in seen:
            continue
        seen.add(seq)
        yield seq, places, "putative", None, None


def cluster_of(start, strand, probe_len, full_len, pam_5prime):
    """CRISPRitz's Cluster Position of an alignment: its start, moved to the PAM's fixed
    end -- less the RNA bulges and plus the DNA bulges on the + strand with a 3' PAM (and
    the - strand with a 5' PAM), the start itself otherwise. Every alignment of one PAM
    site shares it, which is what the merge folds on."""
    if (strand == "+") != pam_5prime:
        return start + probe_len - full_len
    return start


class Reported(object):
    """What the tables report, for completeness."""

    def __init__(self, rows, ctx, starts=None, width=23):
        self.exact = set()
        self.site = defaultdict(list)    # (cluster, variants) -> [carrier sets]
        self.at = set()                  # (sample, cluster, record token)
        self.starts = sorted(starts or ())
        self.width = width
        self.ctx = ctx
        for r in rows:
            t = r["alt"].replace("-", "").upper()
            probe = t if r["strand"] == "+" else rc(t)
            carriers = frozenset(samples_of(r))
            var = frozenset(x for x in r["var"].split(",") if x)
            cl = self.cluster(r["start"], r["strand"], len(probe), r["spacer_col"])
            for s in carriers:
                self.exact.add((s, r["strand"], r["start"], probe, var))
                for v in var:
                    if is_indel_token(v):
                        self.at.add((s, cl, v))
            self.site[(cl, var)].append(carriers)

    def cluster(self, start, strand, probe_len, spacer_col=None):
        full = len(spacer_col) if spacer_col else self.width
        return cluster_of(start, strand, probe_len, full, self.ctx.pam_5prime)

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

    def locus(self, start):
        i = bisect.bisect_left(self.starts, start - self.width)
        return i < len(self.starts) and self.starts[i] <= start + self.width


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
    mode_run, alt_written = run_kind(args.results)
    if args.analysis != "auto":
        mode_run = args.analysis
        alt_written = False if mode_run == "population-level" else alt_written
    # one row per cluster (no alternative alignments): completeness is checked per locus
    popmode = alt_written is False
    rows = []
    starts = {}
    for path, label in tables:
        rows.extend(read_rows(path, label, chrom, lo, hi,
                              starts=starts if label == "integrated" else None))
    indel_rows = [r for r in rows if any(is_indel_token(t) for t in r["var"].split(","))]
    snp_rows = [r for r in rows if not any(is_indel_token(t) for t in r["var"].split(","))]
    sys.stderr.write("check_indel_rows: %d indel/MNV row block(s), %d SNP-only, %s run%s, "
                     "from %s\n" % (len(indel_rows), len(snp_rows), mode_run or "unknown",
                                    " without alternative alignments" if popmode else "",
                                    ", ".join(os.path.basename(p) for p, _l in tables)))
    if mode_run == "population-level" and args.snp_rows == "report":
        sys.stderr.write("check_indel_rows: a population-level run's SNP rows are a "
                         "worst-case union; --snp-rows report is ignored\n")
        args.snp_rows = "off"
    ctx = Ctx(args, pam_motif, pam_len, pam_5prime, guides)
    ref = Reference(args.reference)
    ref.keep_all = args.jobs > 1  # loaded once in this process, shared with the workers
    vcfs_all = [Vcf(p) for p in vcf_paths]
    if not any(v.samples for v in vcfs_all):
        sys.stderr.write("check_indel_rows: the VCF(s) hold no genotypes (a sites-only panel): "
                         "carriers are not checked, only the rows themselves, and completeness has no "
                         "haplotype to enumerate\n")
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

    check_rows = indel_rows + (snp_rows if args.snp_rows == "report" else [])
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
    rep = Report()
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
    width = len(ctx.full(guides[0]))
    for ch in chroms:
        here = [r for r in indel_rows if r["chrom"] == ch]
        st = []
        for b in BLOCKS:
            st.extend(starts.get((ch, b), ()))
        reported[ch] = Reported(here, ctx, starts=st, width=width)

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
