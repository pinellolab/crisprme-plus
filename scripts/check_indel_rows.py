#!/usr/bin/env python3
"""An independent check of CRISPRme's indel (and MNV) off-target rows on real data.

    check_indel_rows.py --results <results dir> --reference <genome dir or FASTA>
                        --vcf <bgzipped VCF> [--vcf ...] | --vcf-dir <dir>
                        --pam <PAM file> --guide <guide or guide file>
                        --mm 4 --bDNA 1 --bRNA 1 [--max-total-edits 6]
                        [--samples <samplesID file> ...]
                        [--chrom chr22 | --region chr22:1-500000]
                        [--mode soundness|completeness|both] [--snp-rows report|off]
                        [--jobs N] [--known known.tsv] [--out failures.tsv]

WHAT IT IS. A second opinion on the indel path, built only from what a user has: the
reference FASTA, the VCF(s) the search read, the samples list, the guide, the PAM file,
the search budget and the result folder. It imports nothing but Python's standard
library and no code of CRISPRme's variant or indel machinery
(``PostProcess/indel_snp_cis.py``, ``analisi_indels_NNN.py``, ``overlay_indel_snps*``,
``resultIntegrator``, ``remove_contiguous_samples``, the CRISPRitz enricher, the
dictionaries, the indel log or the fake contigs), and none of it is copied here. It
reads the VCF with its own BGZF + tabix reader, writes each sample's haplotypes itself
from the VCF 4.x definition of a record, keeps its own map from haplotype bases to
reference coordinates and aligns the guide with its own exhaustive aligner. The CFD
tables CRISPRme ships (``PostProcess/cfd_*_scores.json``) are read as data.

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
(limit-overlapping-records) wherever a target of the record reads (one in budget reads
with the overlapping alleles left out). A sample's alleles are on known slots when its
calls are phased ('|', or haploid). Where a call that touches the target is unphased or
has a missing allele, which of its alleles share the record's haplotype is unknown: when
the record's slot is known and no placement of those calls lets it hold the target (a
phased second indel inside it, say), or every placement does, the slot is judged so;
otherwise the truth is CRISPRme's documented PUTATIVE rule (METHODS section 4): the
sample is listed for the strongest target its own alleles allow -- the guide's base
wherever it carries it; elsewhere, and at the PAM jointly, the alleles that maximize the
alignment's CFD among those the motif allows (a product: where it is 0 for every choice,
a column no choice changes scoring 0, every choice is a strongest one) -- and never for
an allele it does not carry.

THE DOCUMENTED VIEW ONLY CLASSIFIES. Where CRISPRme's output and the truth part, the
rules CRISPRme documents for reading a carrier (METHODS section 4, CHANGELOG 2.8.1 Known
limits), re-implemented here from that text, say whether the difference is a documented
limit (counted and listed under its own class, non-gating) or not (gating):
  * the fake-indel window holds one indel-path record (an indel or an MNV) and the SNPs of
    the enriched genome, which keeps only the last of several SNP records at one position
    (limit-last-alt-wins);
  * a second phased indel inside the target is read on the slot's bases, from the
    target's first base on or up to its last base, through the record, before any call of
    unknown phase sends the carrier to the PUTATIVE rule; an unreadable or different
    reading takes the slot off, and a target that needs two indel-path records is never
    reported (limit-two-indels) -- unless the window's own reading holds its bases from
    the same first base (a repeat unit gained or lost puts them back), where CRISPRme
    keeps the carrier and completeness demands the row;
  * another allele at the record's own position is a phasing artefact, a variant only on
    bases the record deletes changes nothing, and an overlapping indel at another position
    takes the slot off (limit-overlapping-records; a listed carrier the rule takes off is
    carrier-off-by-rule, gating);
  * the carriers of an MNV, of an alt of a multi-ALT SNP line whose number the genotypes
    cannot tell, and of a second indel with an unphased call are resolved PUTATIVE
    (limit-putative-carrier where that differs from their known haplotype, or where the
    rule ties the target with another and the row leaves the carrier off);
  * an MNV beside the record on the slot is in neither store CRISPRme reads (the phased
    indel store keeps records whose REF and ALT lengths differ, the SNP tiers SNVs): a
    target that needs it needs two indel-path records (limit-two-indels), but a carrier
    listed for a sequence such an MNV changes is listed for a sequence it does not have,
    which nothing documents (variant-missing, gating).
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
    * every listed carrier is a sample of the panel (--samples) and holds it, and no
      carrier is missing (CARRIERS).
  An ALT-origin block that names no variant is reported (alt-without-variant). A row that
  names an earlier SNP record of a split site, which the documented limit says CRISPRme's
  stores do not hold, is a note (variant-shadowed, non-gating): the record is real and the
  target is checked with it.

  COMPLETENESS (no real indel off-target is missed -- the real-data brute force). For
  every indel and MNV record of the VCF in the checked region: each carrier's haplotypes
  around it are written (its phased slots; a call of unknown phase is tried both ways),
  identical sequences are aligned once, every in-budget alignment on both strands with the
  PAM that holds a base the record changes (the junction rule) is enumerated and mapped to
  the reference, and each (sample, target) must be reported for that sample at that site
  by its best alignment -- as an integrated row or an alternative alignment whose rows of
  the same site, variant list and carrier set are no worse than it -- or merged into a row
  of the same site by the documented merge into a row no worse than it. A target that
  touches a call of unknown phase is demanded as the strongest target the sample's alleles
  allow there (the PUTATIVE rule). A run that writes no alternative alignments keeps one
  row per merge cluster -- a population-level run, or one made with --ci-cd-test, as the
  documented chr22 self-test (.alt_mode = no-alt-alignments) -- so there each target's
  cluster must be reported with a fewest-mm+b block no worse than it and a highest-CFD
  block scoring at least its CFD.

  NOT CHECKED: the CFD_score_REF/ALT, CFD_risk_score and CRISPR-Bulge columns (and which
  row the CRISPR-Bulge block keeps), Variant_info_spacer+PAM, PAM_creation, the seed
  counts, Variant_MAF and Variant_rsID, Other_motifs, Not_found_in_REF, the annotations;
  a row's own CFD_score is not compared (the merge checks recompute each fold row's CFD
  from its aligned strings).

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

  CARRIERS (soundness). The panel is the samples of --samples (the samplesID files the
  search read), else of the run's .sampleID.txt (CRISPRme writes there the samples of
  those files), else every sample of the VCF(s): CRISPRme lists no carrier outside it. The
  pool is the listed samples and every panel carrier of the row's record. A listed sample
  outside the panel fails (carrier-not-in-panel). A listed sample whose haplotypes do not
  hold the target at its start with every named variant, or that is listed for a weaker
  target than its strongest (PUTATIVE), fails unless the documented view lists it (then a
  limit, with its cause); a sample whose haplotype holds the target at the start with the
  named variants, or whose strongest PUTATIVE target it is, fails when the row omits it,
  unless the documented view leaves it off (a limit). The merge never moves a carrier off
  its sequence's row (same position, variants and aligned target merge carriers; equal
  carrier sets fold to the best alignment), so the row itself must list it. A
  SITES-ONLY VCF (eight columns, no sample) genotypes nobody, and a sites-only dataset's
  rows list no carrier (CHANGELOG 2.8.1; its post-analysis is the population-level one, its
  rows the worst case of each window): so a row that lists no carrier for a record a
  sites-only VCF holds is that dataset's, and its strings, target and variants are checked,
  its carriers not (limit-no-genotypes, one line per row block). The decision is the
  record's, not the run's: beside a genotyped dataset, as the default index with
  hg38_mega, the genotyped dataset's rows list their carriers and are checked as above, and
  a record both hold is demanded for its genotyped carriers by completeness. A row that
  lists no carrier for a record no sites-only VCF holds fails (no-carrier-listed), and a
  sample listed for a record only sites-only VCFs hold carries nothing
  (carrier-without-indel).

  Merge (completeness). The pipeline clusters a guide's rows on a chromosome with Cluster
  Positions within --merge bp (3) of the previous, chained, whatever the strand.
  CRISPRitz's Cluster Position is the start moved to the PAM's fixed end (less the RNA
  bulges and plus the DNA bulges on the + strand with a 3' PAM, the - strand with a 5'
  one), which every alignment of one PAM site shares. Inside a cluster, variant rows of
  one (Cluster Position, variant list) -- and for an indel/MNV row one aligned target --
  merge their carriers, and the rows of one site, variant list and carrier set fold into
  the first of them in each file's order: the fewest file sorts on the sorting criteria in
  turn (--sorting-criteria, crisprme.py's default mm+bulges,mm), the CFD file on the CFD
  (as the analysis writes it, to three decimals) and then on the same criteria
  (merge_close_targets_cfd.sh hands both files --sorting-criteria; a later CRISPRme branch
  hands the CFD file --sorting-criteria-scoring, so where a run's command line names one, a
  CFD tie no worse under either is taken as no worse). Two alignments of one
  target share its start, variants, carriers and site, so they always fold. The site is the
  cluster position on the haplotype carrying the record: the Cluster Position plus the
  rank of the row's first base among the bases an insertion adds (0 on any other base), as
  the fix branch's merge compares rows that start inside one insertion (c3a75fe; at
  d7a6fff it compared the Cluster Position alone, which folds two targets of one
  insertion). CRISPRme's fold key holds no strand, while the alignments of one target
  share theirs: rows of the other strand are other targets, which a fold across strands
  would lose, so a fold group is taken on one strand. A fold group is NO WORSE than a
  target when its fewest-mm+b block's key under the sorting criteria is at most the least
  key of the target's alignments, and its highest-CFD block, its CFD recomputed from its
  own aligned strings, scores at least the target's best CFD at the three decimals the
  merge sorts on -- and, where the two surely print the same three decimals, has a key no
  greater than the least key of the target's alignments of that CFD. (Not compared: the
  merge's keys after the criteria -- the Cluster Position, the carrier count, the aligned
  text -- which only order alignments equal in all of these.) So a (sample, target) is
  satisfied by a row reporting
  exactly it with its variants whose fold group (the rows of its site naming those
  variants with its carrier set) is no worse than it, or by the rows of its site naming
  the same variants whose carrier set is the target's own as CRISPRme writes it (every
  carrier the documented view lists for it, and nobody else but carriers it resolves
  PUTATIVE) and that are no worse than it. A PUTATIVE target is the carrier's own: it is
  held to the alignments under which the rule picks it for that carrier (the choice is
  the alignment's highest CFD), and one the rule ties with others to its own best. The
  verdicts depend on nothing but the input (no set's order). In a run with one row per
  cluster, the cluster of
  a target is the integrated row with a block within --merge bp of its Cluster Position,
  or of one the in-budget alignments of the enriched reference (the SNP alleles its stores
  hold: the last record of a split site) and of the record's window chain to it, each
  within --merge bp of the next (every record the chain reaches is read, up to 3 kb
  either side); there is no other fallback. That row must be as good as the target in
  both blocks: no more edits by the first sorting criterion (the merge picks between a
  cluster's reference and variant rows on mm+bulges, or the CFD, alone) and at least its
  CFD. A row's CFD is always recomputed from its aligned strings:
  the merge sorts on the CFD of a row's own alignment (CFD_score_ALT), while CFD_score
  prints the higher of its REF and ALT scores.

LIMITS (non-gating, counted per class, one TSV line per row block or per (record,
target)): limit-last-alt-wins, limit-two-indels, limit-overlapping-records,
limit-putative-carrier (all above), limit-contig-end (a record within one target of a
contig's end, or a row whose REF text runs past it), limit-n-bases (a stretch of N in a
record's window: no target is aligned across it), limit-no-genotypes (a row of a
sites-only dataset, which lists no carrier: the row itself is checked, one line per
block; and in completeness one line per chromosome counting the records only sites-only
VCFs hold, which no sample carries), and in the SNP-row report limit-snp-row-indel and
limit-snp-row-unphased; variant-shadowed is a note. Two non-gating classes are this
checker's own COVERAGE GAPS, not CRISPRme limits: coverage-putative-cap (a carrier with
more than 10 positions of unknown phase within a target's reach of a record: its targets
that touch one are not enumerated, 2^k haplotypes each) and coverage-unindexed-long-ref
(a VCF read without a .tbi holds a record whose REF is longer than the 20 kb an unindexed
read looks back: a deletion that long, starting further before a region, is not read
there; give the VCF a .tbi).
Rows on a chromosome no VCF covers, or with no reference sequence, and table lines shorter
than their header are reported (gating): nothing is skipped silently.

OUTPUT. A TSV of every finding with its reason class, gating flag and evidence, a
per-class summary on stdout, exit status 1 on any gating finding (2 on a usage error).
``--known`` lists findings already reported upstream: they turn non-gating
(known-defect), and a listed one that no longer shows is reported (known-not-seen).
SNP-only rows of a per-sample run can be checked in report-only mode
(``--snp-rows report``: their strings, target and listed carriers, and the phased
carriers a SNP row of their site leaves out); population-level SNP rows are a documented
worst-case union and are out of scope, as are a sites-only dataset's SNP rows, which are
population-level whatever the search's mode (one limit-no-genotypes line each). Without
--max-total-edits the cap is the run's own
(its .Params.txt), else complete-search's 6; without --samples the panel is the run's
.sampleID.txt, else every sample of the VCF(s). Rows are keyed by guide: one guide's row
stands for no demand of another.

VCF ACCESS. A .tbi beside a VCF makes every read a seek. Without one, each VCF is read
once into an in-memory block index (nothing is written), which costs one decompression of
the file; a canonical install (one VCF per chromosome, 28 GB) is then feasible with
--jobs. An unindexed read starts 20 kb before its region, so give a .tbi for a panel
with deletions longer than that (HPRC's chr22 has 582; completeness lists each VCF whose
reads met one, coverage-unindexed-long-ref). Several VCFs are read as one panel: a
genotyped one beside a sites-only one is the default index with hg38_mega. A sites-only
VCF's split sites are not read as shadowed: CRISPRme's mega index holds earlier records
of some (its rows name them), so which it holds cannot be read off the VCF's order.
``--jobs N`` splits the work by region over N processes.
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
NO_VARIANT = ("", "NA", "n")


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

    def lines_at(self, voffset):
        """(virtual offset, line) of each line from a virtual offset on: a line's virtual
        offset is its first byte's (the block's compressed offset << 16 | the byte's
        offset in the block), as a tabix chunk names its ends."""
        self.fh.seek(voffset >> 16)
        pos = voffset & 0xFFFF
        tail, tail_v = b"", None
        while True:
            at = self.fh.tell()
            block = self._member()
            if block is None:
                break
            while True:
                nl = block.find(b"\n", pos)
                if nl < 0:
                    if pos < len(block):
                        if not tail:
                            tail_v = (at << 16) | pos
                        tail += block[pos:]
                    break
                if tail:
                    line, v = tail + block[pos:nl], tail_v
                    tail, tail_v = b"", None
                else:
                    line, v = block[pos:nl], (at << 16) | pos
                yield v, line.decode("utf-8", "replace")
                pos = nl + 1
            pos = 0
        if tail:
            yield tail_v, tail.decode("utf-8", "replace")

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
        self.long_ref = 0      # the longest REF an unindexed read met (0: none longer)
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
                # each chunk is read from its first line to its end (a virtual offset, the
                # end exclusive), and no line twice: the chunks come sorted, and a line an
                # earlier chunk's read reached is not yielded again
                done = -1
                for cb, ce in self.tbi.chunks(chrom, beg, end):
                    for v, line in bg.lines_at(cb):
                        if v >= ce:
                            break
                        if v <= done:
                            continue
                        done = v
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
                    if len(r["ref"]) > self.LONGEST:
                        # such a deletion is not read where a region starts more than
                        # LONGEST after it: noted (coverage-unindexed-long-ref)
                        self.long_ref = max(self.long_ref, len(r["ref"]))
                    yield r
        finally:
            bg.close()

    def _parse(self, line, chrom, beg, end):
        """A record as a dict; the genotype columns stay unsplit (only the lines a window
        keeps are split). A sites-only line (eight columns, no FORMAT) is a record with no
        genotypes."""
        if not line or line.startswith("#"):
            return None
        head, sep, tail = line.rstrip("\n").partition("\t")
        if head != chrom or not sep:
            return None
        f = tail.split("\t", 8)
        if len(f) < 7:
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
    (GT must be the first FORMAT key, as VCF 4.x requires when it is present). A sample
    whose name is None (outside the panel) is left out."""
    out = {}
    if not tail or not fmt.startswith("GT"):
        return out
    fields = tail.split("\t")
    if fmt != "GT":
        fields = [g.partition(":")[0] for g in fields]
    homref = _HOMREF
    for i in [i for i, g in enumerate(fields[:len(names)]) if g not in homref]:
        if names[i] is None:
            continue
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
    (METHODS section 4: their carriers are resolved as unphased ones).

    ``panel`` (a set, or None for every sample): the samples of the search's panel; the
    calls of any other sample are not read (CRISPRme lists no carrier outside it)."""

    def __init__(self, refseq, chrom, beg, end, vcfs, needed=None, panel=None):
        self.chrom = chrom
        self.refseq = refseq
        self.beg = max(0, beg)
        self.end = end
        self.needed = needed
        self.panel = panel
        self.vcfs = vcfs
        self.samples = []
        self.alleles = []
        self.by_key = defaultdict(list)
        self.line_lo = []       # line id -> 0-based start of its REF
        self.line_hi = []
        self.line_alleles = []  # line id -> [Allele or None per ALT]
        self.calls = []         # line id -> {sample: (allele numbers, phased)}
        self.line_sites_only = []  # line id -> whether its VCF has no sample columns
        self.shadowed = set()
        self.untied = set()
        self.snv_at = defaultdict(list)  # reference coordinate -> SNV alleles there
        self._order = None
        seen = set()
        for v in vcfs:
            names = [s if panel is None or s in panel else None for s in v.samples]
            for s in names:
                if s is not None and s not in seen:
                    seen.add(s)
                    self.samples.append(s)
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
                self.line_sites_only.append(not v.samples)
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
                if len(lids) < 2 or not v.samples:
                    # a sites-only VCF's split sites are not read as shadowed: the mega
                    # index CRISPRme ships holds earlier records of some (its rows name
                    # them, as chr22:25093750 C>G before C>T), so which ones it holds
                    # cannot be read off the order of the VCF given
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

    def sites_only(self, key):
        """Whether a sites-only VCF (no sample columns) holds the record ``key``: the
        rows a sites-only dataset gives of it list no carrier (CHANGELOG 2.8.1)."""
        return any(self.line_sites_only[a.line] for a in self.by_key.get(key, ()))

    def genotyped(self, key):
        """Whether a VCF with sample columns holds the record ``key``."""
        return any(not self.line_sites_only[a.line] for a in self.by_key.get(key, ()))

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


def guide_of(spacer_col, pam_len, pam_5prime):
    """The guide of a row, from its Spacer+PAM column (the guide and N's over the PAM)."""
    spc = spacer_col.upper()
    return spc[pam_len:] if pam_5prime else spc[:len(spc) - pam_len]


def row_cfd(cfd_model, crrna, aligned, pam_len):
    """The CFD of a row's alignment, from its own aligned strings (None without tables):
    what the analysis writes in the column the merge sorts on (CFD_score_ALT of the
    results; their CFD_score is the higher of the REF and ALT scores)."""
    if cfd_model is None or not crrna or aligned in NO_VARIANT or len(crrna) != len(
            aligned):
        return None
    return cfd_model.score(crrna, aligned.upper(), pam_len)


def read_rows(path, table, chrom=None, lo=None, hi=None, margin=0, pam_5prime=False,
              loci=None, bad=None, pam_len=3, cfd_model=None, novar=None):
    """Every block (dict) of a results table that names a variant. Rows whose start is
    within [lo - margin, hi + margin) are read; those within [lo, hi) are marked for the
    soundness check ('own'). ``loci`` (a dict) gets, per (guide, chromosome), every
    integrated line's (guide, Cluster Position of each block, the fewest-mm+b block's
    edits under each sorting criterion, the highest-CFD block's CFD -- recomputed from
    its aligned strings when ``cfd_model`` is given, the value the merge sorts on, else
    as printed --, line) for the per-cluster completeness check; ``bad`` (a list) gets
    the lines shorter than the header, ``novar`` (a list) the ALT-origin blocks of
    [lo, hi) that name no variant."""
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
            guide = (guide_of(row[col["Spacer+PAM"]], pam_len, pam_5prime)
                     if "Spacer+PAM" in col else "")
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
                cfd_col = col.get("CFD_score_(%s)" % b)
                block = {
                    "table": table, "line": line_no, "block": b, "chrom": ch,
                    "spacer_col": row[col["Spacer+PAM"]], "guide": guide,
                    "cfd": row[cfd_col] if cfd_col is not None else "",
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
                }
                if var in NO_VARIANT:
                    # an ALT row is made by variants and must name them
                    if origin == "alt" and novar is not None and block["own"]:
                        novar.append(block)
                    continue
                out.append(block)
            if loci is not None and line_blocks:
                try:
                    crit = edits(int(row[col["Mismatches_(fewest_mm+b)"]]),
                                 int(row[col["Bulges_(fewest_mm+b)"]]))
                except (KeyError, ValueError):
                    crit = {}
                try:
                    cfd = float(row[col["CFD_score_(highest_CFD)"]])
                except (KeyError, ValueError):
                    cfd = None
                try:
                    # the merge sorts on the CFD of the row's own alignment (its
                    # CFD_score_ALT); CFD_score prints the higher of REF and ALT
                    b = "highest_CFD"
                    alt = row[col["Aligned_protospacer+PAM_ALT_(%s)" % b]]
                    shown = (alt if alt not in NO_VARIANT
                             else row[col["Aligned_protospacer+PAM_REF_(%s)" % b]])
                    calc = row_cfd(cfd_model, row[col["Aligned_spacer+PAM_(%s)" % b]],
                                   shown, pam_len)
                    if cfd is not None and calc is not None:
                        cfd = calc
                except KeyError:
                    pass
                loci.setdefault((guide, ch), []).append(
                    (guide, [c for _b, c in line_blocks], crit, cfd, line_no))
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
    """(analysis mode, alternative alignments written, merge window, the sorting
    criteria of the fewest block (--sorting-criteria, crisprme.py's default
    'mm+bulges,mm' when the command line gives none), the scoring criteria the command
    line gives (--sorting-criteria-scoring, else None), total-edits cap or None) of a
    run, from its own markers. A population-level run writes no alternative alignments
    by design, nor does a run made with --ci-cd-test (the documented chr22 self-test):
    both keep one row per cluster (.alt_mode = no-alt-alignments). The cap is the one the
    search used (.Params.txt's Max_total_edits, else --max-total-edits in
    .command_line.txt)."""
    mode = alt = cap = scoring = None
    merge, criterion = 3, "mm+bulges,mm"
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
            criterion = m.group(1)
        m = re.search(r"--sorting-criteria-scoring\s+(\S+)", v)
        if m:
            scoring = m.group(1)
        m = re.search(r"--max-total-edits\s+(-?\d+)", v)
        if m:
            cap = int(m.group(1))
    for name in (".Params.txt", "Params.txt"):
        p = os.path.join(resdir, name)
        if os.path.exists(p):
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) == 2 and f[0] == "Max_total_edits":
                        try:
                            cap = int(f[1])
                        except ValueError:
                            pass
            break
    if mode == "population-level" and alt is None:
        alt = False
    return mode, alt, merge, criterion, scoring, cap


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
        self.panel = None     # the search's samples (a set), None for every VCF sample

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
            # the record is real, so the row is not wrong for naming it (its target is
            # checked against it); what it contradicts is a documented limit: a note
            rep.add(False, "variant-shadowed", row,
                    "%s is an earlier SNP record of a split site, which the documented "
                    "limit says CRISPRme's stores do not hold (CHANGELOG 2.8.1 Known "
                    "limits): this index holds it" % t)
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
        if unk and h is not None:
            # the record's slot is known: what its known alleles say comes first. When no
            # placement of the calls of unknown phase lets the slot hold the target (a
            # phased second indel inside it, say), or every placement does, the slot is
            # judged; only when they decide is the PUTATIVE rule the truth
            can, must = placements(win, tgt, main_al, present, unk, named)
            if can is not None and (not can or must):
                holds, unk = bool(must), []
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


def placements(win, tgt, main_al, present, unk, named, cap=PUTATIVE_CAP):
    """(can, must) for a known slot of the record: whether some placement of the calls of
    unknown phase ``unk`` (each of a line's alleles on the slot, or none) lets the slot
    read the target at its start through the record with every named variant, and
    whether every placement does. (None, None) when a placement gives overlapping alleles
    near the target (no sequence) or there are too many lines to try."""
    by_line = defaultdict(list)
    for a in unk:
        by_line[a.line].append(a)
    if len(by_line) > cap:
        return None, None
    want = {a.key for a in named}
    can, must = False, True
    for combo in itertools.product(*[[None] + v for v in by_line.values()]):
        chosen = [a for a in combo if a is not None]
        hap = make_hap(win.refseq, [main_al] + present + chosen, tgt.lo, tgt.hi)
        if hap.clashes(tgt.span[0] - 1, tgt.span[1] + 1):
            return None, None
        keys = {a.key for a in present} | {a.key for a in chosen}
        ok = want <= keys and bool(reads_at(hap, tgt.start, tgt.probe, main_al))
        can = can or ok
        must = must and ok
    return can, must


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
    the options that maximize the CFD of the whole alignment. CFD is a product of one
    factor per spacer column and the PAM's factor, so where the best product is above 0
    each column takes its own best options; where it is 0 (a column no choice changes
    scores 0, as a DNA bulge whose base the table gives 0, or every option of a column
    does, or every PAM the motif allows) every choice ties at 0 and each is a strongest
    one. Without CFD (it does not apply) every option ties."""
    if tgt.crrna is None or tgt.aligned is None:
        return None
    n_cols = len(tgt.crrna)
    cols = base_columns(tgt.aligned, tgt.strand)
    pam = set(ctx.pam_cols(n_cols))
    n_spacer = n_cols - ctx.pam_len
    positions = Cfd.positions(tgt.crrna, n_spacer)
    ori = (lambda b: b) if tgt.strand == "+" else comp
    cfd = ctx.cfd
    allowed = {}
    scores = {}       # base index -> {option: its column's CFD factor} (a choice to make)
    pam_codes = []
    code_cols = set()
    fixed = 1.0       # the factors no choice changes
    for i, _q in tgt.codes:
        opts = opts_by_code[i]
        col = cols[i]
        code_cols.add(col)
        if col in pam:
            pam_codes.append(i)
            continue
        g = tgt.crrna[col].upper()
        o_opts = {ori(b): b for b in opts}
        if g != "-" and g in o_opts:
            allowed[i] = {o_opts[g]}   # the guide's base wins (a factor of 1)
            continue
        if cfd is None:
            allowed[i] = set(opts)
            continue
        scores[i] = {b: cfd.factor(g, ori(b), positions[col]) for b in opts}
    if cfd is not None:
        for col in range(n_spacer):
            if col not in code_cols:
                fixed *= cfd.factor(tgt.crrna[col], tgt.aligned[col], positions[col])
    pam_ok = None
    pam_top = 1.0
    combos = []
    if pam_codes:
        target_bases = list(tgt.aligned.replace("-", "").upper())
        # target_bases is in the target's orientation: base index i (+ order) sits at
        # position i on the + strand, n-1-i on the - strand
        def at(i):
            return i if tgt.strand == "+" else len(target_bases) - 1 - i
        for combo in itertools.product(*[sorted(opts_by_code[i]) for i in pam_codes]):
            tb = list(target_bases)
            for i, b in zip(pam_codes, combo):
                tb[at(i)] = ori(b)
            pseq = "".join(tb[:ctx.pam_len] if ctx.pam_5prime else tb[-ctx.pam_len:])
            if not all(iupac_match(m, x) for m, x in zip(ctx.pam_motif, pseq)):
                continue
            score = cfd.pam.get(pseq[-2:], 0.0) if cfd is not None else 0.0
            combos.append((score, combo))
        pam_top = max((s for s, _c in combos), default=0.0)
    elif cfd is not None:
        tail = tgt.aligned.replace("-", "").upper()
        pam_top = cfd.pam.get((tail[-ctx.pam_len:])[-2:], 0.0)
    best = fixed * pam_top
    for sc in scores.values():
        best *= max(sc.values())
    tie = cfd is None or best <= 0.0
    for i, sc in scores.items():
        top = max(sc.values())
        allowed[i] = set(sc) if tie else {b for b, v in sc.items() if v >= top - 1e-12}
    if pam_codes:
        pam_ok = ({c for _s, c in combos} if tie
                  else {c for s, c in combos if s >= pam_top - 1e-12})
    determined = all(len(v) == 1 for v in allowed.values()) and (
        not pam_codes or len(pam_ok) == 1)
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
    verdicts = []     # per slot: (the rules list it, an MNV on it changes the target)
    put_mnv = []      # the same flag for each slot the rules resolve PUTATIVE
    putative = None   # a call of unknown phase (or an untied alt) on a slot the rules keep
    for h in slots:
        snvs, seconds = [], []
        off = None
        unsure = puts = None
        # an MNV is in no store CRISPRme reads for an indel's carrier (the phased indel
        # store keeps records whose REF and ALT lengths differ, the SNP tiers SNVs): the
        # slot is read without it, whatever it changes
        mnv_here = False
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
                    elif a.kind == "indel":
                        off = "complex"
                    # an overlapping MNV is not read: the slot stays as without it
                    continue
                if a.kind == "snv":
                    if not touches(a, s0, s1):
                        continue
                    if r == "?":
                        puts = puts or "putative"
                        continue
                    if a in win.untied:
                        puts = puts or "putative-untied"
                        continue
                    if a in win.shadowed:
                        causes.add("last-alt-wins")
                        continue
                    snvs.append(a)
                elif a.kind == "indel":
                    if r == "?":
                        if changes_target(a, tgt):
                            unsure = "putative"
                        continue
                    seconds.append(a)
                elif changes_target(a, tgt):  # an MNV, not read
                    if r == "yes":
                        mnv_here = True
                    else:
                        causes.add("second-indel")
        if unsure:
            # a second indel inside the target with an unphased call: listed PUTATIVE,
            # even when a phased second indel takes the slot off (METHODS section 4)
            causes.add(unsure)
            return _documented_putative(win, ctx, tgt, sample, causes)
        # the slot's phased variants come first (CRISPRme takes a slot off before it
        # resolves anything PUTATIVE): an overlapping record at another position, and a
        # second indel inside the target whose reading fails
        if off:
            causes.add(off)
            verdicts.append((False, mnv_here))
            continue
        inside = [a for a in seconds if changes_target(a, tgt)]
        if not inside:
            if puts:
                putative = putative or puts
                put_mnv.append(mnv_here)
                continue
            window = make_hap(win.refseq, [main_al] + snvs, tgt.lo, tgt.hi)
            verdicts.append((bool(reads_at(window, tgt.start, tgt.probe, main_al)),
                             mnv_here))
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
        if ok and puts:
            putative = putative or puts   # read through, then resolved PUTATIVE
            put_mnv.append(mnv_here)
            continue
        verdicts.append((ok, mnv_here))
    # 'mnv-unseen': the rules list the carrier only through slots read without an MNV
    # that changes the target there -- a sequence the carrier does not have, which no
    # documented rule lists (gating); a slot the rules leave off whose MNV changes the
    # target needs a second indel-path record for it (the two-indel limit)
    if putative is not None:
        # a slot the rules keep has a call of unknown phase in the target: PUTATIVE
        causes.add(putative)
        if all(put_mnv) and not any(v for v, _m in verdicts):
            causes.add("mnv-unseen")
        return _documented_putative(win, ctx, tgt, sample, causes)
    listing = [m for v, m in verdicts if v]
    if listing:
        if all(listing):
            causes.add("mnv-unseen")
        return "lists", causes
    if any(m for _v, m in verdicts):
        causes.add("second-indel")
    return "off", causes


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
            ("mnv-unseen", "limit-two-indels"),   # only where the carrier is left off
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
    "carrier-not-in-panel": "listed with %s, and not samples of the search's panel "
                            "(--samples, else the run's .sampleID.txt)",
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


MNV_UNSEEN = ("an MNV on their slot changes the target; CRISPRme reads no MNV beside the "
              "record (its phased indel store keeps records whose REF and ALT lengths "
              "differ, its SNP tiers SNVs), and no documented rule lists a carrier for a "
              "sequence it does not have")


def carrier_verdict(win, ctx, tgt, sample, listed, lids):
    """(gating, class[, detail]) for one sample of a row, or None when it agrees."""
    tv = truth(win, tgt, sample, lids)
    kind = tv[0]
    if kind == "none":
        if listed and not win.genotyped(tgt.main.key):
            return (True, "carrier-without-indel",
                    "only sites-only VCF(s) hold the record: nobody's genotype carries it, "
                    "and a sites-only dataset's rows list no carrier")
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
                # CRISPRme resolves it PUTATIVE and its rule ties this target with
                # another: it may be listed on the other (completeness looks there)
                return False, "limit-putative-carrier", (
                    "documented: %s; a PUTATIVE tie with another target"
                    % ",".join(sorted(causes)))
            return unexplained_or(causes)
        if verdict in ("lists", "maybe"):
            if "mnv-unseen" in causes:
                return True, "variant-missing", MNV_UNSEEN
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
            if "mnv-unseen" in causes:
                return True, "variant-missing", MNV_UNSEEN
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
    tgt, cls, detail = row_target(win, row, named, main)
    if tgt is None:
        rep.add(gating, cls, row, detail)
        return
    for a in named:
        if a is not main and not any(s is a for s in tgt.src):
            rep.add(gating, "variant-outside-target", row,
                    "%s changes no base of the target (%d..%d)"
                    % (a.token(win.chrom), tgt.span[0], tgt.span[1]))
    samples = samples_of(row)
    tok = main.token(win.chrom)
    if not samples:
        if win.sites_only(main.key):
            # a row of a sites-only dataset (no genotypes): CRISPRme documents no carrier
            # on it, so its strings, target and variants are all there is to check. The
            # decision is the record's, not the run's: beside a genotyped dataset (the
            # default index with hg38_mega) that dataset's rows still list carriers
            rep.add(False, "limit-no-genotypes", row,
                    "%s is a record of a sites-only VCF (no sample columns), whose "
                    "dataset's rows list no carrier (CHANGELOG 2.8.1): the row's strings, "
                    "target and variants are checked, no carrier is%s"
                    % (tok, "; the record's carriers in the genotyped VCF(s) are "
                            "demanded by completeness" if win.genotyped(main.key) else ""))
            return
        rep.add(gating, "no-carrier-listed", row, "the row lists no carrier")
        return
    lids = local_lines(win, tgt)
    listed = set(samples)
    out = defaultdict(list)
    why = defaultdict(set)
    if ctx.panel is not None:
        # CRISPRme lists no carrier outside the search's panel
        for s in sorted(listed - ctx.panel):
            out[(True, "carrier-not-in-panel")].append(s)
        listed &= ctx.panel
    pool = set(listed)
    for a in win.by_key.get(main.key, ()):
        pool.update(win.carriers(a))
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
    for (g, cls), who in sorted(out.items(), key=lambda kv: kv[0][1]):
        extra = sorted(why.get((g, cls), ()))
        detail = "%d sample(s) %s%s: %s" % (
            len(who), CARRIER_TEXT.get(cls, "%s") % tok,
            " (%s%s)" % ("; ".join(extra[:4]), "; ..." if len(extra) > 4 else "")
            if extra else "", ",".join(who[:12]))
        rep.add(g and gating, cls, row, detail, samples=",".join(who))


def check_snp_row(row, win, ctx, rep, snp_listed=None):
    """The same checks on a SNP-only row of a per-sample run, report-only (CHANGELOG
    2.8.1 Known limits: a SNP row lists carriers from SNP genotypes only): its own
    strings, its target, every listed carrier, and the carriers it leaves out -- a sample
    whose phased slot holds exactly the named SNPs inside the target is listed on a row
    of this site naming them (``snp_listed``: (guide, chromosome, Cluster Position,
    variants) -> the samples SNP rows there list; such rows merge their carriers). A row
    that lists no carrier and names only records a sites-only VCF holds is a sites-only
    dataset's, whose SNP post-analysis is the population-level one whatever the search's
    mode (CHANGELOG 2.8.1): a worst-case union, out of scope, so one line says so."""
    if not alignment_checks(row, ctx, rep, False):
        return
    ref_text_check(row, win, rep, False)
    named = named_alleles(row, win, rep, False)
    if named is None:
        return
    if not samples_of(row) and named and all(win.sites_only(a.key) for a in named):
        rep.add(False, "limit-no-genotypes", row,
                "a SNP row of a sites-only dataset (no carrier; %s held by a sites-only "
                "VCF): its post-analysis is the population-level one, a documented "
                "worst-case union, which the SNP-row report leaves out"
                % ",".join(a.token(win.chrom) for a in named))
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
        if ctx.panel is not None and s not in ctx.panel:
            out["carrier-not-in-panel"].append(s)
            continue
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
    # the carriers it leaves out (phased calls only: an unphased carrier is listed for
    # every subset it may hold, METHODS section 4, which is not checked here)
    var = frozenset(x for x in row["var"].split(",") if x and x != "NA")
    cp = cluster_of(row["start"], row["strand"], n, len(row["spacer_col"]),
                    ctx.pam_5prime)
    site = (snp_listed or {}).get((row.get("guide", ""), row["chrom"], cp, var), set())
    want = {a.key for a in named}
    pool = set()
    for a in named:
        for x in win.by_key.get(a.key, ()):
            pool.update(win.carriers(x))
    pool -= set(samples_of(row)) | set(site)
    for s in sorted(pool):
        calls = {lid: win.calls[lid][s] for lid in lids if s in win.calls[lid]}
        ploidies = {len(c[0]) for c in calls.values()}
        touching = [lid for lid in calls if win.line_lo[lid] <= span[1]
                    and win.line_hi[lid] - 1 >= span[0]]
        if len(ploidies) != 1 or any(
                (not calls[lid][1] and len(set(calls[lid][0])) > 1)
                or None in calls[lid][0] for lid in touching):
            continue
        for h in range(ploidies.pop()):
            snvs = [a for lid in touching for a in win.line_alleles[lid]
                    if a is not None and a.kind == "snv" and calls[lid][0][h] == a.number
                    and touches(a, span[0], span[1])]
            if {a.key for a in snvs} != want or any(a in win.shadowed for a in snvs):
                continue
            if reads_at(make_hap(win.refseq, snvs, lo, hi), row["start"], probe):
                out["carrier-missing"].append(s)
                break
    text = {"limit-snp-row-indel": "hold an indel (or MNV) inside the target on the slot "
                                   "of the row's SNPs: a SNP row's carriers come from SNP "
                                   "genotypes only (CHANGELOG 2.8.1, Known limits)",
            "carrier-missing": "are not listed, nor on another row of this site naming "
                               "these SNPs, and their phased slot holds exactly these SNPs "
                               "inside the target, which it reads (report-only)",
            "limit-last-alt-wins": "hold an earlier SNP record of a split site inside the "
                                   "target (CRISPRme reads the last record)",
            "limit-overlapping-records": "hold overlapping alleles inside the target",
            "limit-snp-row-unphased": "have a call of unknown phase inside the target: a "
                                      "SNP row lists every subset a carrier may hold "
                                      "(METHODS section 4), not checked here",
            "carrier-without-target": "hold no haplotype with this target and every named "
                                      "variant",
            "carrier-not-in-panel": "are not samples of the search's panel (--samples, "
                                    "else the run's .sampleID.txt)"}
    for cls, who in sorted(out.items()):
        rep.add(False, cls, row, "%d %ssample(s) %s: %s"
                % (len(who), "" if cls == "carrier-missing" else "listed ", text[cls],
                   ",".join(who[:12])), samples=",".join(who))


# --------------------------------------------------------------------------- #
# completeness
# --------------------------------------------------------------------------- #
class Best(object):
    """The best of a target's alignments: ``best`` the fewest mismatches+bulges,
    mismatches and bulges (each on its own), ``cfd`` the highest CFD (None without CFD),
    and ``keys`` every alignment's (mismatches, bulges, CFD), which the merge's sort keys
    are read from: what the rows standing for the target after the merge must match.
    ``real`` is False for a summary of several alignments (its minima are no one
    alignment's, so they give no key)."""

    __slots__ = ("best", "cfd", "keys")

    def __init__(self, best, cfd, real=True):
        self.best = dict(best)
        self.cfd = cfd
        self.keys = {(best["mm"], best["bulges"], cfd)} if real else set()

    def add(self, best, cfd, real=True):
        for c, v in best.items():
            if v < self.best[c]:
                self.best[c] = v
        if cfd is not None and (self.cfd is None or cfd > self.cfd):
            self.cfd = cfd
        if real:
            self.keys.add((best["mm"], best["bulges"], cfd))


class Demand(Best):
    """A target a carrier must be reported for: exact (its known haplotype) or the
    strongest of its PUTATIVE choices (any one of ``alts`` when the rule leaves a tie).
    ``best`` and ``cfd`` are the best over its alignments (and its alternatives), ``own``
    each alternative's own Best (by strand, start, + reading and variants), ``aligned``
    each alternative's alignments (a Target each), and ``k`` the rank of its first base
    among the bases an insertion adds (0 on any other base)."""

    __slots__ = ("kind", "alts", "k", "own", "aligned")

    def __init__(self, kind, alts, best, cfd, k):
        Best.__init__(self, best, cfd, real=False)   # add_alt notes each alignment
        self.kind = kind
        self.alts = alts  # [(Target, variant tokens frozenset)]
        self.k = k
        self.own = {}
        self.aligned = {}

    @property
    def total(self):
        return self.best["mm+bulges"]

    def add(self, best, cfd, k=None):
        Best.add(self, best, cfd, real=False)
        if k is not None:
            self.k = min(self.k, k)

    def add_alt(self, t, toks, best, cfd):
        """One alignment of alternative ``t`` (its own best, and the demand's)."""
        key = (t.strand, t.start, t.probe, toks)
        self.keys.add((best["mm"], best["bulges"], cfd))
        b = self.own.get(key)
        if b is None:
            self.own[key] = Best(best, cfd)
        else:
            b.add(best, cfd)
        ts = self.aligned.setdefault(key, [])
        if not any(x.crrna == t.crrna and x.aligned == t.aligned for x in ts):
            ts.append(t)

    def alignments(self, t, toks):
        """The alignments of alternative ``t`` (Targets: a PUTATIVE choice, and so
        CRISPRme's documented view of a carrier, can depend on the alignment)."""
        return self.aligned.get((t.strand, t.start, t.probe, toks)) or [t]

    def best_of(self, t, toks):
        """The Best a row standing for alternative ``t`` must match: its own alignments'
        (a PUTATIVE target the rule ties with others is compared with itself)."""
        return self.own.get((t.strand, t.start, t.probe, toks), self)


def edits(mm, bulges):
    """An alignment's value under each sorting criterion of the fewest block."""
    return {"mm+bulges": mm + bulges, "mm": mm, "bulges": bulges}


def completeness(win, ctx, rep, reported, region, popmode, stats=None):
    """The real-data brute force over the indel and MNV records of ``region``. ``stats``
    counts the records only sites-only VCFs hold (no sample carries them: nothing to
    demand, and they are counted, not skipped in silence)."""
    lo_r, hi_r = region
    for gi, guide in enumerate(ctx.guides):
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
                if gi == 0 and stats is not None and not win.genotyped(al.key):
                    stats["sites_only"] = stats.get("sites_only", 0) + 1
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


def record_reads(win, ctx, full, cache, main_al, rest, lo, hi):
    """Whether a target of the record reads on a slot holding it and ``rest`` (the slot's
    other alleles, those overlapping the record left out: the one-record reading)."""
    hap = make_hap(win.refseq, [main_al] + rest, lo, hi)
    if hap.clashes(main_al.lo, main_al.hi - 1):
        return False
    for _st, off, _sa, ta, _mm, _dn, _rn in scan_cached(cache, ctx, full, hap.seq, []):
        width = len(ta.replace("-", ""))
        if off + width <= len(hap.seq) and holds_edit(hap, off, width, main_al):
            return True
    return False


def one_record_reading(win, main_al, alleles, t, toks, lo, hi):
    """Whether the fake-indel window's reading of a slot -- the record and the slot's SNV
    alleles, no second indel-path record -- holds the target's bases from its first base
    on, with the same variants: CRISPRme then lists the carrier for that row (METHODS
    section 4: a unit a tandem repeat gains or loses leaves the bases as they are)."""
    hap = make_hap(win.refseq, [main_al] + [a for a in alleles if a.kind == "snv"], lo, hi)
    for i in reads_at(hap, t.start, t.probe, main_al):
        sp = span_of(hap.place[i:i + t.n])
        if hap.clashes(sp[0] - 1, sp[1] + 1):
            continue
        if demand_tokens(win, main_al, hap.place[i:i + t.n], hap.src[i:i + t.n],
                         [])[0] == toks:
            return True
    return False


def comp_record(win, ctx, rep, reported, al, carriers, guide, full, lo, hi, popmode,
                reach=None):
    """The demands of one record. ``reach``: the reference coordinates a target holding
    the record's edit can cover; an unknown call outside it cannot change such a target,
    so only the calls inside it are tried both ways."""
    reach = reach or (lo, hi)
    chrom = win.chrom
    lids = win.lines_in(lo, hi)
    groups = defaultdict(list)
    for s in sorted(carriers):   # a set: sorted, so that nothing depends on hash order
        _a, call = main_call(win, al.key, s)
        groups[(call, win.signature(s, lids))].append(s)
    cache = {}
    skipped = []
    demands = defaultdict(set)   # (kind, key of the demand) -> samples
    demand_of = {}
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
                    # the record overlaps another allele of the slot: the VCF gives no
                    # sequence there. It is counted where a target of the record reads
                    # on the slot with the overlapping alleles left out (the one-record
                    # reading CRISPRme's window makes); a record with no target either
                    # way has nothing left unjudged
                    over = {id(a) for z in own for a in z[2] if a is not main_al}
                    rest = [a for a in present + chosen if id(a) not in over]
                    if record_reads(win, ctx, full, cache, main_al, rest, lo, hi):
                        clash_who.update(members)
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
                    # demanded too: a row reporting it satisfies it, and a miss is the
                    # documented limit (one record per fake-indel window) -- unless the
                    # window's own reading holds the same bases from the same first base,
                    # where CRISPRme documents listing the carrier: demanded exactly
                    kind = "exact"
                    if others and not one_record_reading(win, main_al, present + chosen,
                                                         t, toks, lo, hi):
                        kind = "second"
                    key = (kind, strand, t.start, probe, toks)
                    demands[key].update(members)
                    keep_demand(demand_of, key, t, toks, edits(mm, dn + rn),
                                target_cfd(ctx, sa, ta))
        if need_putative:
            putative_targets(win, ctx, rep, al, main_al, s, members, lids, lo, hi, full,
                             known, cache, skipped, demands, demand_of, capped, reach)
    resolve(win, ctx, rep, reported, al, carriers, lids, demands, demand_of, popmode,
            guide)
    if clash_who:
        who = sorted(clash_who)
        how = sorted(clash_how)
        rep.add(False, "limit-overlapping-records", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d sample(s) hold %s on a slot whose alleles overlap where a target of it "
            "reads (one in budget reads there with the overlapping alleles left out), so "
            "the VCF gives no one sequence there and they are not demanded (overlapping "
            "on one slot: %s%s): %s"
            % (len(who), tok, "; ".join(how[:4]), "; ..." if len(how) > 4 else "",
               ",".join(who[:12])), samples=",".join(who))
    if capped:
        who = sorted(capped)
        rep.add(False, "coverage-putative-cap", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "this checker's own coverage gap, not a CRISPRme limit: %d sample(s) have more "
            "than %d positions of unknown phase within a target's reach of %s, and their "
            "targets there that touch one are not enumerated (2^k haplotypes each): %s"
            % (len(who), PUTATIVE_CAP, tok, ",".join(who[:12])), samples=",".join(who))
    if skipped:
        rep.add(False, "limit-n-bases", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d window(s) around %s hold an N and are not aligned (no base pairs with N)"
            % (len(skipped), tok))


def keep_demand(demand_of, key, t, toks, best, cfd):
    """Record one alignment of a demanded target: the target is one sequence at one
    start (or two copies of it one insertion holds), and its many alignments differ in
    edits and CFD, so the demand keeps the best of each (what the row standing for it
    must match)."""
    d = demand_of.get(key)
    if d is None:
        d = demand_of[key] = Demand(key[0], [(t, toks)], best, cfd, t.place[0][1])
    else:
        d.add(best, cfd, t.place[0][1])
    d.add_alt(t, toks, best, cfd)


def clash_text(win, zones):
    return {"+".join(a.token(win.chrom) for a in z[2]) for z in zones}


def target_cfd(ctx, sa, ta):
    if ctx.cfd is None:
        return None
    return ctx.cfd.score(sa, ta.upper(), ctx.pam_len)


def slot_allows(win, t, main_al, named, slot_sets):
    """Whether some known slot of the record can hold target ``t`` for some placement of
    its calls of unknown phase (placements); True when that cannot be told."""
    for present, unknown in slot_sets:
        unk = [a for a in unknown if touches(a, t.span[0], t.span[1])]
        can, _must = placements(win, t, main_al, present, unk, named)
        if can is None or can:
            return True
    return False


def putative_targets(win, ctx, rep, al, main_al, s, members, lids, lo, hi, full, known,
                     cache, skipped, demands, demand_of, capped, reach=None):
    """The PUTATIVE demands of one carrier group: every in-budget target its own alleles
    allow with the record (any choice at each SNV position: the documented rule ignores
    phase), grouped by alignment; the strongest choice of each is demanded. With the
    record's slot known, only alignments that touch a call of unknown phase are, and only
    choices some placement of those calls lets a slot of the record hold (its known
    alleles, a phased second indel say, can rule a choice out)."""
    main_call_ = main_call(win, al.key, s)[1]
    unknown_all = []
    slot_sets = []
    if known:
        for h in record_slots(main_call_, main_al.number)[0]:
            present_h, unknown_h = slot_alleles(win, s, main_al, main_call_, h, lids)
            unknown_all.extend(unknown_h)
            slot_sets.append((present_h, unknown_h))
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
    groups = defaultdict(dict)  # alignment -> {probe: (Target, mm, bulges)}
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
                g[probe] = (t, mm, dn + rn)
    # the strongest choice depends on the alignment (its CFD), so the alignments under
    # which a target is this group's strongest are the group's own: a demand is the
    # targets the rule picks (one, or a tie) together with every alignment that picks
    # them, and carriers share it only when their alignments are the same (pooling
    # another group's alignments would hold this one to a target it is not listed for)
    picks = defaultdict(list)   # (strand, start, + readings, variants) -> [(aid, alts)]
    for aid, members_t in groups.items():
        alts = []
        for probe, (t, mm, bul) in members_t.items():
            if not putative_check(win, ctx, t, s)[1]:
                continue
            snvs = [x for x in t.src if x is not None and x is not main_al
                    and x.kind == "snv"]
            if known and not slot_allows(win, t, main_al, snvs, slot_sets):
                continue
            toks = frozenset([al.token(win.chrom)] + [x.token(win.chrom) for x in snvs])
            alts.append((t, toks, mm, bul))
        if not alts:
            continue
        picks[(aid[0], alts[0][0].start, tuple(sorted(t.probe for t, *_x in alts)),
               tuple(sorted(tuple(sorted(x)) for _t, x, *_y in alts)))].append((aid, alts))
    for (strand, start, probes, toks_all), found in picks.items():
        key = ("putative", strand, start, probes, toks_all,
               frozenset(aid for aid, _alts in found))
        demands[key].update(members)
        for _aid, alts in found:
            best = {c: min(edits(mm, bul)[c] for _t, _x, mm, bul in alts)
                    for c in ("mm+bulges", "mm", "bulges")}
            best_cfd = (max((target_cfd(ctx, t.crrna, t.aligned) or 0.0)
                            for t, *_x in alts) if ctx.cfd else None)
            k = min(t.place[0][1] for t, *_x in alts)
            if key not in demand_of:
                demand_of[key] = Demand("putative", [(t, x) for t, x, *_y in alts], best,
                                        best_cfd, k)
            else:
                demand_of[key].add(best, best_cfd, k)
            d = demand_of[key]
            for t, x, mm, bul in alts:
                d.add_alt(t, x, edits(mm, bul), target_cfd(ctx, t.crrna, t.aligned))


def resolve(win, ctx, rep, reported, al, carriers, lids, demands, demand_of, popmode,
            guide):
    """Each demand against the tables; a miss is classified by the documented view."""
    chrom = win.chrom
    full_len = len(ctx.full(guide))
    fold_sets = {}

    def first_k(r):
        return reported.first_k(win, r)

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
            ok, why = reported.locus(guide, chrom, cl, d, lambda: cluster_chain(
                win, ctx, guide, cl, reported.merge))
            if d.kind == "second":
                if not ok:
                    rep.add(False, "limit-two-indels", row,
                            "%d sample(s) hold this target only with a second "
                            "indel-path record on the slot (one record per fake-indel "
                            "window), and its cluster %s: %s"
                            % (len(who), why, ",".join(who[:12])), samples=",".join(who))
                continue
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
        worse = {}
        listed_worse = set()   # listed for the target exactly, but not by its best
        # the rows that report an alternative exactly: their carriers, and whether their
        # fold group (the rows of their site naming those variants with those carriers,
        # which the merge folds into the best by each block's order: two alignments of
        # one target always fold) is worse than it
        exact = [reported.exact_groups(guide, t, toks, d.best_of(t, toks), first_k)
                 for t, toks in d.alts]
        for s in who:
            # a target that needs two indel-path records is reported, for the sample,
            # when a row lists it there under either record's window (the documented
            # limit is that it is not: its alignment is not compared)
            if d.kind == "second" and any(reported.has(guide, s, t.strand, t.start,
                                                       t.probe) for t, _x in d.alts):
                continue
            hit = False
            for groups in exact:
                for who_e, bad in groups:
                    if s not in who_e:
                        continue
                    if bad is None:
                        hit = True
                        break
                    listed_worse.add(s)
                    worse.setdefault(s, "the rows of its site naming its variants with the "
                                        "carriers of the row that lists them for it: %s"
                                     % bad)
                if hit:
                    break
            if hit:
                continue
            for t, toks in d.alts:
                hcl = cluster_of(t.start, t.strand, t.n, full_len, ctx.pam_5prime) + (
                    t.place[0][1] if d.kind == "putative" else d.k)
                # the carriers CRISPRme's documented rules list for the target, under each
                # of its alignments (a PUTATIVE choice can depend on it)
                sets = []
                for ta in d.alignments(t, toks):
                    fk = (ta.strand, ta.start, ta.probe, toks, ta.crrna, ta.aligned)
                    if fk not in fold_sets:
                        fold_sets[fk] = documented_carriers(win, ctx, ta, carriers, lids)
                    if fold_sets[fk] not in sets:
                        sets.append(fold_sets[fk])
                for cm, pm in sets:
                    ok, why = reported.folded(guide, s, hcl, t.strand, toks, cm | set(who),
                                              pm, d.best_of(t, toks), first_k)
                    if ok:
                        hit = True
                        break
                    if why and s not in worse:
                        worse[s] = why
                if hit:
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
            if verdict in ("lists", "maybe") and s in listed_worse:
                # listed for the target, by an alignment worse than its best: the merge
                # kept the worse of its alignments (or lost the best)
                out[(True, "missing-best-alignment", "")].append(s)
            elif verdict == "lists":
                out[(True, "missing-row" if d.kind == "exact" else "missing-putative-row",
                     "")].append(s)
            elif verdict == "maybe":
                if reported.listed_at(guide, s, cl, al.token(chrom)):
                    # CRISPRme's rule ties this target with another and the sample is
                    # listed at the site naming the record: likely on the other one
                    out[(False, "limit-putative-carrier",
                         "%s; a PUTATIVE tie, listed at the site on another row"
                         % ",".join(sorted(causes)))].append(s)
                    continue
                out[(True, "missing-putative-row", "")].append(s)
            else:
                g, cls = unexplained_or(causes)[:2]
                out[(g, cls, ",".join(sorted(causes)))].append(s)
        for (g, cls, why), ss in sorted(out.items(), key=lambda kv: kv[0][1]):
            text = ("hold%s this target at %d (+ reading %s, best %d edits%s)%s and no row "
                    "reports it for them no worse than that: neither one of exactly this "
                    "target and variants whose rows of its site with its carriers are no "
                    "worse, nor a fold into rows of its site naming the same variants with "
                    "the target's carriers%s"
                    % ("s" if g else "", t0.start, t0.probe, d.total,
                       "" if d.cfd is None else ", CFD %.4f" % d.cfd,
                       " as the strongest PUTATIVE target" if d.kind == "putative" else "",
                       "" if g else " (documented: %s)" % why))
            if g:  # a failure is one line per (sample, target)
                for x in ss:
                    rep.add(True, cls, row, "%s %s%s" % (
                        x, text, " (%s)" % worse[x] if x in worse else ""), samples=x)
            else:
                rep.add(False, cls, row, "%d sample(s) %s: %s"
                        % (len(ss), text, ",".join(ss[:12])), samples=",".join(ss))


def documented_for(win, ctx, d, s, lids):
    """The documented verdict on a demand: 'lists' if CRISPRme's rules report any of its
    alternatives for the sample under any of its alignments, 'maybe' if they leave the
    choice open."""
    seen = set()
    best = "off"
    for t0, toks in d.alts:
        for t in d.alignments(t0, toks):
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
    list for the target (under its alignment ``t``), and those they resolve PUTATIVE
    without fixing the choice."""
    listed, maybe = set(), set()
    lids = win.lines_in(t.lo, t.hi)
    cache = {}
    for s in sorted(carriers):
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


# --------------------------------------------------------------------------- #
# a merge cluster's extent (one row per cluster): the alignments CRISPRme chains
# --------------------------------------------------------------------------- #
def _feasible_sets(matches, target, n, m, max_mm, max_dna, max_rna, max_total, pam_first,
                   pam_last, pam_5prime):
    """_feasible on a target whose bases are sets (IUPAC: any allele of a position): an
    alignment exists if some choice of each position's alleles has one in budget. Used
    only to find where a merge cluster can extend, so it is kept apart from the aligner's
    hot path."""
    inf = max_mm + 1
    cost = {(0, 0): 0}
    i = 0
    for i in range(n + 1):
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
                ok = not target[j].isdisjoint(matches[i])
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


def enriched_windows(win, lo, hi):
    """[lo, hi) of the reference, and of each indel-path record's haplotype there, with
    every SNV allele of the window CRISPRme's stores hold as an option at its position
    (not an earlier record of a split site: only the last is kept): the enriched genome
    and the fake-indel windows, read as CRISPRitz reads an IUPAC code. (base sets, places)
    each. A merge cluster holds the rows of every record, so every record that touches
    the stretch is read (the alleles of one multi-allelic site put their targets a few
    bases apart)."""
    records = [None]
    seen = set()
    for lid in win.lines_in(lo, hi):
        for a in win.line_alleles[lid]:
            if a is not None and a.kind != "snv" and a.key not in seen:
                seen.add(a.key)
                records.append(a)
    out = []
    for al in records:
        hap = make_hap(win.refseq, [al] if al is not None else [], lo, hi)
        sets = []
        for b, pl, src in zip(hap.seq, hap.place, hap.src):
            opts = {b}
            if src is None and pl[1] == 0:
                opts.update(a.body for a in win.snv_at.get(pl[0], ())
                            if a not in win.shadowed)
            sets.append(frozenset(opts))
        out.append((sets, hap.place))
    return out


def feasible_cps(ctx, full, sets, places):
    """The Cluster Positions of the PAM sites of ``sets`` (both strands) where some
    alignment of the guide is in budget."""
    n = len(full)
    pam_len, p5 = ctx.pam_len, ctx.pam_5prime
    motif = full[:pam_len] if p5 else full[n - pam_len:]
    proto = n - pam_len
    matches = [_MATCHES.get(c.upper(), frozenset()) for c in full]
    pam_first, pam_last = (0, pam_len - 1) if p5 else (n - pam_len, n - 1)
    motif_sets = [_MATCHES.get(c.upper(), frozenset()) for c in motif]
    size = len(sets)
    rc_sets = [frozenset(comp(x) for x in s) for s in reversed(sets)]
    out = set()
    for strand, ss in (("+", sets), ("-", rc_sets)):
        for q in range(0, size - pam_len + 1):
            if any(ss[q + j].isdisjoint(motif_sets[j]) for j in range(pam_len)):
                continue
            for width in range(proto - ctx.brna, proto + ctx.bdna + 1):
                beg = (q - width) if not p5 else q
                end = (q + pam_len) if not p5 else (q + pam_len + width)
                if beg < 0 or end > size:
                    continue
                piece = ss[beg:end]
                if any("N" in x for x in piece):
                    continue
                if not _feasible_sets(matches, piece, n, end - beg, ctx.mm, ctx.bdna,
                                      ctx.brna, ctx.cap, pam_first, pam_last, p5):
                    continue
                off = beg if strand == "+" else size - end
                out.add(cluster_of(places[off][0], strand, end - beg, n, p5))
    return out


def chain_from(c0, cps, merge):
    """The run of Cluster Positions, each within ``merge`` of the next, that holds c0."""
    pts = sorted(set(cps) | {c0})
    i = pts.index(c0)
    lo = hi = i
    while lo > 0 and pts[lo] - pts[lo - 1] <= merge:
        lo -= 1
    while hi + 1 < len(pts) and pts[hi + 1] - pts[hi] <= merge:
        hi += 1
    return pts[lo:hi + 1]


def cluster_chain(win, ctx, guide, c0, merge, step=150, most=3000):
    """The Cluster Positions CRISPRme's merge can chain with c0: those of the in-budget
    alignments, on either strand, of the enriched reference and of every indel-path
    record's window (enriched_windows), each within ``merge`` bp of the next, from c0 on
    (CRISPRme clusters a guide's rows on a chromosome so, whatever their strand). The
    window grows while the chain reaches its edge (up to ``most`` bp either side)."""
    full = ctx.full(guide)
    lo, hi = c0 - step, c0 + step
    wide = False
    while True:
        cps = set()
        a = max(0, lo - len(full) - 10)
        b = min(len(win.refseq), hi + len(full) + 10)
        if not wide and (win.needed is not None or a < win.beg or b > win.end):
            # the chain reaches past what the window read (a tile's reads a pad around
            # the tile): every record it can still reach is read, once
            wide = True
            win = Window(win.refseq, win.chrom,
                         min(win.beg, max(0, c0 - most - step - len(full) - 10)),
                         max(win.end, c0 + most + step + len(full) + 10), win.vcfs,
                         panel=win.panel)
        for sets, places in enriched_windows(win, a, b):
            cps |= feasible_cps(ctx, full, sets, places)
        members = chain_from(c0, cps, merge)
        if (members[0] - merge > lo and members[-1] + merge < hi) or hi - lo >= 2 * most:
            return members
        lo -= step
        hi += step


class Entry(object):
    """One reported row block, for completeness: its carriers, Cluster Position, + reading,
    edits under each sorting criterion and CFD -- recomputed from its own aligned strings
    when the CFD tables apply (``printed`` keeps the table's value)."""

    __slots__ = ("row", "carriers", "cp", "block", "crit", "cfd", "printed", "probe")

    def __init__(self, row, carriers, cp, probe=None, cfd_model=None, pam_len=3):
        self.row = row
        self.carriers = carriers
        self.cp = cp
        self.probe = probe
        self.block = row["block"]
        self.crit = {}
        try:
            mm, bul = int(row["mm"]), int(row["bulges"])
            self.crit = edits(mm, bul)
        except (KeyError, ValueError):
            pass
        try:
            self.printed = float(row.get("cfd", ""))
        except ValueError:
            self.printed = None
        calc = row_cfd(cfd_model, row.get("crrna", ""), row.get("alt", ""), pam_len)
        self.cfd = calc if calc is not None else self.printed


def cfd3(x):
    """A CFD as the analysis writes it ("{:.3f}"), the value the merge sorts on."""
    return float("%.3f" % x)


def cfd_below(a, b):
    """Whether CFD ``a`` (a row's) sorts below CFD ``b`` (a target's) in the merge, which
    sorts on three decimals. A difference of float rounding alone is none: a value within
    1e-9 of a rounding boundary (x.xxx5) may print either way, as products taken in
    another order do, so the row's is rounded up there and the target's down."""
    return cfd3(a + 1e-9) < cfd3(b - 1e-9) and b - a > 1e-9


def cfd_same(a, b):
    """Whether CFDs ``a`` and ``b`` surely print the same three decimals (neither within
    1e-9 of a rounding boundary): only then does the merge's CFD file go on to its next
    sort key for them."""
    return len({cfd3(a - 1e-9), cfd3(a + 1e-9), cfd3(b - 1e-9), cfd3(b + 1e-9)}) == 1


def crit_key(crits, mm, bulges):
    """An alignment's key in the merge's sort on the sorting criteria ``crits``."""
    e = edits(mm, bulges)
    return tuple(e[c] for c in crits)


def crit_text(crits, key):
    return ", ".join("%d %s" % (v, c) for v, c in zip(key, crits))


def row_first_k(win, row):
    """The rank of a reported row's first base among the bases its insertion adds (0 when
    it is any other base), read from the row's own target on the reference with its named
    variants: its cluster position on the haplotype is its Cluster Position plus that.
    Two places hold the same target only in a repeat it lies wholly in; the first is
    taken, as CRISPRme takes it (indel_snp_cis.row_first_base, on the fix branch). None
    when the row's target does not read at its own start through its record (a row
    soundness fails, as 2.8.1's rows placed past an insertion's anchor): such a row stands
    for no other alignment."""
    named = []
    for t in row["var"].split(","):
        try:
            _c, p, r, a = parse_token(t)
        except (ValueError, IndexError):
            return None
        found = win.by_key.get((p, r, a))
        if not found:
            return None
        named.append(found[0])
    mains = [a for a in named if a.kind != "snv"]
    if len(mains) != 1:
        return None
    main = mains[0]
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    width = main.hi - main.lo
    lo = max(0, min(row["start"], main.lo) - n - width - 60)
    hi = min(len(win.refseq), max(row["start"], main.hi) + 2 * n + width + 60)
    hap = make_hap(win.refseq, named, lo, hi)
    if hap.clash:
        return None
    hits = reads_at(hap, row["start"], probe, main)
    return hap.place[hits[0]][1] if hits else None


class Reported(object):
    """What the tables report, for completeness: every row block that names a variant, by
    guide -- its carriers, its Cluster Position, its block's edits and CFD -- and every
    integrated line's clusters, edits and CFD (runs with one row per cluster)."""

    def __init__(self, rows, ctx, loci=None, merge=3, criterion="mm+bulges,mm",
                 scoring=None):
        self.exact = set()
        self.exact_any = set()           # (guide, sample, strand, start, + reading)
        self.by_var = defaultdict(list)  # (guide, variants) -> [Entry]
        self.by_target = defaultdict(list)  # (guide, strand, start, + reading, variants)
        self.at = set()                  # (guide, sample, cluster, record token)
        self.ctx = ctx
        self.merge = merge
        # the sorting criteria, in order: the merge sorts the fewest file on them and the
        # CFD file on the CFD, then on them (merge_close_targets_cfd.sh hands both files
        # --sorting-criteria); the first alone is compared where the reference row and
        # the variant row of one cluster are chosen between (on mm+bulges or CFD only)
        self.criteria = tuple(c for c in criterion.split(",")
                              if c in ("mm+bulges", "mm", "bulges")) or ("mm+bulges", "mm")
        self.criterion = self.criteria[0]
        # --sorting-criteria-scoring: a CRISPRme whose merge gives it to the CFD file (a
        # later branch's does; 2.8.x's never reads it) breaks a CFD tie on it instead, so
        # a tie no worse under either is taken as no worse
        self.scoring = tuple(c for c in (scoring or "").split(",")
                             if c in ("mm+bulges", "mm", "bulges")) or None
        self._k = {}
        self._groups = {}   # (guide, variants) -> {(strand, carriers, site): [Entry]}
        self.loci = {}
        for gk, lines in (loci or {}).items():
            pos = sorted((c, k) for k, ln in enumerate(lines) for c in ln[1])
            self.loci[gk] = (pos, lines)
        for r in rows:
            t = r["alt"].replace("-", "").upper()
            probe = t if r["strand"] == "+" else rc(t)
            carriers = frozenset(samples_of(r))
            var = frozenset(x for x in r["var"].split(",") if x)
            g = r.get("guide", "")
            cl = cluster_of(r["start"], r["strand"], len(probe), len(r["spacer_col"]),
                            ctx.pam_5prime)
            for s in carriers:
                self.exact.add((g, s, r["strand"], r["start"], probe, var))
                self.exact_any.add((g, s, r["strand"], r["start"], probe))
                for v in var:
                    if is_indel_token(v):
                        self.at.add((g, s, cl, v))
            e = Entry(r, carriers, cl, probe, ctx.cfd, ctx.pam_len)
            self.by_var[(g, var)].append(e)
            self.by_target[(g, r["strand"], r["start"], probe, var)].append(e)

    def has(self, guide, s, strand, start, probe, var=None):
        """Whether the sample is listed for the target with the variants ``var`` (any
        variants: None)."""
        if var is None:
            return (guide, s, strand, start, probe) in self.exact_any
        return (guide, s, strand, start, probe, var) in self.exact

    def first_k(self, win, row):
        key = (row["table"], row["line"], row["block"])
        if key not in self._k:
            self._k[key] = row_first_k(win, row)
        return self._k[key]

    def group(self, guide, var, carriers, site, strand, first_k):
        """The rows that stand for one target's alignments after the merge's fold: those
        naming ``var`` with the carrier set ``carriers`` on ``strand`` whose target starts
        at haplotype cluster position ``site`` (a row whose target is not at its own start
        belongs to none). CRISPRme's fold key holds no strand, but the alignments of one
        target share it: a row of the other strand is another target, which a fold would
        lose, so it stands for none of this one's."""
        groups = self._groups.get((guide, var))
        if groups is None:
            groups = self._groups[(guide, var)] = defaultdict(list)
            for e in self.by_var.get((guide, var), ()):
                k = first_k(e.row)
                if k is not None:
                    groups[(e.row["strand"], e.carriers, e.cp + k)].append(e)
        return groups.get((strand, carriers, site), [])

    def exact_groups(self, guide, t, toks, best, first_k):
        """[(carriers, why)] for every row block reporting target ``t`` exactly (its
        strand, start, + reading and the variants ``toks``): its carrier set, and why its
        fold group -- every row of its site naming those variants with that carrier set,
        which the merge folds into the best by each block's order, so the group holds the
        target's best alignment if the run reported it -- is worse than ``best`` (None:
        no worse)."""
        out = []
        seen = {}
        for e in self.by_target.get((guide, t.strand, t.start, t.probe, toks), ()):
            k = first_k(e.row)
            key = (e.carriers, None if k is None else e.cp + k)
            if key not in seen:
                if k is None:
                    seen[key] = "its row's target is not at its own start"
                else:
                    seen[key] = self.worse(self.group(guide, toks, e.carriers, e.cp + k,
                                                      t.strand, first_k), best)
                out.append((e.carriers, seen[key]))
        return out

    def folded(self, guide, sample, hcl, strand, var, exact_who, putative_who, d,
               first_k):
        """(ok, why): a row of the same site naming the same variants whose carriers are
        the target's own, and no worse than it in each block: the merge folds those rows
        into the first of them by each file's order (_fold_same_carriers), which worse()
        compares on the merge's sort keys. The site
        is the cluster position on the haplotype carrying the record, ``hcl``: the
        Cluster Position plus the rank of the first base among an insertion's added
        bases (first_k of a row), which is the Cluster Position itself on any row that
        does not start inside an insertion. (On d7a6fff the merge compares rows at their
        Cluster Position alone, so it folds two targets that start on different bases of
        one insertion; the fix branch's c3a75fe compares them on the haplotype.) The rows
        are on the target's ``strand``: CRISPRme's fold key holds none, and a fold across
        strands loses another target."""
        allowed = set(exact_who) | set(putative_who)
        groups = defaultdict(list)
        for e in self.by_var.get((guide, var), ()):
            if e.cp > hcl or sample not in e.carriers or e.row["strand"] != strand:
                continue
            if not (exact_who <= e.carriers <= allowed):
                continue
            k = first_k(e.row)
            if k is None or e.cp + k != hcl:
                continue   # another site, or a row whose target is not at its start
            groups[e.carriers].append(e)
        why = None
        for es in groups.values():
            bad = self.worse(es, d)
            if bad is None:
                return True, None
            why = why or bad
        return False, why

    def worse(self, es, d):
        """None when one fold group's rows ``es`` are no worse than ``d`` (a Best) by the
        merge's sort keys in each block, else why. The fewest file sorts on the sorting
        criteria in order (mm+bulges, then mm, by default); the CFD file on the CFD at
        the three decimals the analysis writes, then on the same criteria. So the
        fewest-mm+b block's key must be at most the least key of the target's alignments,
        and the highest-CFD block's CFD at least the target's best -- and, where the two
        surely print the same three decimals, its key at most the least key of the
        target's alignments of that CFD. A group the run lists has a row in each block
        (each file keeps one per group). A row's CFD is recomputed from its aligned
        strings. Not compared: the keys after the criteria (the Cluster Position, the
        carrier count, the aligned text), which only order alignments equal in all of
        these."""
        crits = self.criteria
        few_rows = [e for e in es if e.block == "fewest_mm+b"]
        cfd_rows = [e for e in es if e.block == "highest_CFD"]
        if not few_rows or not cfd_rows:
            return "no %s block of its site lists its carriers" % (
                "fewest-mm+b" if not few_rows else "highest-CFD")
        keys = sorted(d.keys, key=lambda k: (k[0], k[1], -1.0 if k[2] is None else k[2]))
        few = [tuple(e.crit[c] for c in crits) for e in few_rows
               if all(c in e.crit for c in crits)]
        if few and keys:
            rk = min(few)
            tk = min(crit_key(crits, mm, b) for mm, b, _c in keys)
            if rk > tk:
                if rk[0] > tk[0]:
                    return ("the row of its site with its carriers has a fewest-mm+b block "
                            "of %d (%s), the target %d" % (rk[0], crits[0], tk[0]))
                return ("the row of its site with its carriers has a fewest-mm+b block of "
                        "%s, the target %s (worse by the merge's next sorting criterion, "
                        "--sorting-criteria %s)"
                        % (crit_text(crits, rk), crit_text(crits, tk), ",".join(crits)))
        cfd = [e for e in cfd_rows if e.cfd is not None]   # none: no CFD column
        if cfd and d.cfd is not None:
            top = max(cfd, key=lambda e: e.cfd)
            if cfd_below(top.cfd, d.cfd):
                return ("the row of its site with its carriers has a highest-CFD block of "
                        "%.3f (%.6f, printed %s), the target %.3f (%.6f)"
                        % (cfd3(top.cfd), top.cfd, top.row.get("cfd", "") or "nothing",
                           cfd3(d.cfd), d.cfd))
            why = None
            for crs in [crits] + ([self.scoring] if self.scoring and
                                  self.scoring != crits else []):
                tied = [crit_key(crs, mm, b) for mm, b, c in keys
                        if c is not None and cfd_same(c, d.cfd)]
                if not tied:
                    return None
                tk = min(tied)
                ok = False
                for e in cfd:
                    if cfd_below(e.cfd, d.cfd):
                        continue
                    if (not cfd_same(e.cfd, d.cfd) or not all(c in e.crit for c in crs)
                            or tuple(e.crit[c] for c in crs) <= tk):
                        ok = True   # a higher (or not surely equal) CFD, or no worse key
                        break
                if ok:
                    return None
                rk = min(tuple(e.crit[c] for c in crs) for e in cfd
                         if cfd_same(e.cfd, d.cfd))
                why = why or ("the row of its site with its carriers has a highest-CFD block "
                              "of %.3f with %s, the target %.3f with %s (a tie at three "
                              "decimals, which the merge breaks on --sorting-criteria %s)"
                              % (cfd3(top.cfd), crit_text(crs, rk), cfd3(d.cfd),
                                 crit_text(crs, tk), ",".join(crs)))
            return why
        return None

    def listed_at(self, guide, sample, cluster, tok):
        return (guide, sample, cluster, tok) in self.at

    def near(self, pos, cps):
        out = set()
        for c in cps:
            i = bisect.bisect_left(pos, (c - self.merge, -1))
            while i < len(pos) and pos[i][0] <= c + self.merge:
                out.add(pos[i][1])
                i += 1
        return out

    def locus(self, guide, chrom, c, d, chain):
        """(ok, why) for a target of Cluster Position c in a run with one row per merge
        cluster. The cluster's row is an integrated line with a block within --merge bp
        of c, or of a Cluster Position chained to c (chain(): the in-budget alignments
        CRISPRme can chain, each within --merge bp of the next); none is a miss. It must
        have a fewest-mm+b block no worse than the target by the first sorting criterion
        and a highest-CFD block of at least its CFD at the three decimals the merge sorts
        on (the row's CFD recomputed from its aligned strings)."""
        pos, lines = self.loci.get((guide, chrom), ([], []))
        near = self.near(pos, [c])
        if not near:
            members = chain()
            near = self.near(pos, members)
            if not near:
                return False, (
                    "no integrated row lies in its cluster (Cluster Positions %d..%d: the "
                    "in-budget alignments of the reference, its SNP alleles and the indel "
                    "records there that chain with it within %d bp)"
                    % (members[0], members[-1], self.merge))
        why = []
        crit = self.criterion
        for k in sorted(near):
            _g, _clusters, r_crit, r_cfd, ln = lines[k]
            v = r_crit.get(crit) if r_crit else None
            if v is not None and v > d.best[crit]:
                why.append("the row at line %d has a fewest-mm+b block of %d (%s)"
                           % (ln, v, crit))
                continue
            if d.cfd is not None and r_cfd is not None and cfd_below(r_cfd, d.cfd):
                why.append("the row at line %d has a highest-CFD block of %.3f (%.6f), "
                           "the target %.3f" % (ln, cfd3(r_cfd), r_cfd, cfd3(d.cfd)))
                continue
            return True, ""
        return False, "; ".join(why)

# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def read_panel(paths):
    """The sample IDs of samplesID files (the first column of each line; '#' lines and a
    SAMPLE_ID header are skipped)."""
    out = set()
    for path in paths:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                sid = s.split()[0]
                if sid != "SAMPLE_ID":
                    out.add(sid)
    return out


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


def row_span(r, pad):
    """The reference stretch a soundness check of row ``r`` reads records over: its start
    and target and the REF span of every record it names, ``pad`` either side (a long
    deletion's far junction, and the carriers' alleles there, lie past the start by the
    deletion's length)."""
    lo = r["start"]
    hi = lo + len(r["alt"].replace("-", ""))
    for t in r["var"].split(","):
        try:
            _c, p, ref, _a = parse_token(t)
        except (ValueError, IndexError):
            continue
        lo = min(lo, p - 1)
        hi = max(hi, p - 1 + len(ref))
    return max(0, lo - pad), hi + pad


def run_unit(unit):
    """One work unit: ('sound' | 'comp', chrom, lo, hi, rows). Returns (the findings, the
    unit's counts: {'sites_only': records only sites-only VCFs hold})."""
    kind, ch, tlo, thi, rows = unit
    ctx = G["ctx"]
    rep = Report()
    stats = {"chrom": ch}
    refseq = G["ref"].seq(ch)
    pad = 2 * (len(ctx.guides[0]) + ctx.pam_len) + 1000
    vcfs = G["vcfs_for"][ch]
    if kind == "sound":
        needed = []
        for a, b in sorted(row_span(r, pad) for r in rows):
            if needed and a <= needed[-1][1]:
                needed[-1] = (needed[-1][0], max(needed[-1][1], b))
            else:
                needed.append((a, b))
        win = Window(refseq, ch, needed[0][0], needed[-1][1], vcfs, needed=needed,
                     panel=ctx.panel)
        for r in sorted(rows, key=lambda x: x["start"]):
            if any(is_indel_token(x) for x in r["var"].split(",")):
                check_row(r, win, ctx, rep)
            else:
                check_snp_row(r, win, ctx, rep, G.get("snp_listed", {}))
    else:
        win = Window(refseq, ch, max(0, tlo - pad), thi + pad, vcfs, panel=ctx.panel)
        # a record of the tile whose REF runs far (a long deletion) has its targets'
        # haplotypes written past what was read: read again, to its far end
        one = max(len(ctx.full(g)) for g in ctx.guides) + max(ctx.bdna, ctx.brna) + 2
        far = max([a.hi for a in win.alleles if a.kind != "snv" and tlo <= a.lo < thi]
                  or [0]) + 2 * one + 20   # the far end of what comp_record writes
        if far > win.end:
            win = Window(refseq, ch, max(0, tlo - pad), far + pad, vcfs, panel=ctx.panel)
        completeness(win, ctx, rep, G["reported"][ch], (tlo, thi), G["popmode"], stats)
    for v in vcfs:
        if v.long_ref:
            stats.setdefault("long_ref", {})[v.path] = v.long_ref
    return rep.rows, stats


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
    p.add_argument("--max-total-edits", type=int, default=None,
                   help="the cap on mismatches + bulges the search used (default: the "
                        "run's own, from its .Params.txt or .command_line.txt, else 6, "
                        "complete-search's default; a negative value: no cap)")
    p.add_argument("--samples", action="append", default=[],
                   help="a samplesID file the search read (repeatable: the panel is their "
                        "union; the sample ID is the first column, '#' lines are skipped). "
                        "Default: the run's .sampleID.txt (CRISPRme writes there the "
                        "samples of those files), else every sample of the VCF(s)")
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
    mode_run, alt_written, merge, criterion, scoring, run_cap = run_kind(args.results)
    # the total-edits cap: the one given, else the run's own, else complete-search's 6
    if args.max_total_edits is None:
        if run_cap is None:
            sys.stderr.write("check_indel_rows: the run does not say its --max-total-edits "
                             "(no .Params.txt): taking 6, complete-search's default\n")
            run_cap = 6
        args.max_total_edits = run_cap
    elif run_cap is not None and run_cap != args.max_total_edits:
        sys.stderr.write("check_indel_rows: WARNING: --max-total-edits %d, while the run "
                         "searched with %d\n" % (args.max_total_edits, run_cap))
    if args.max_total_edits < 0:
        args.max_total_edits = None   # no cap
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
    novar = []
    for path, label in tables:
        rows.extend(read_rows(path, label, chrom, lo, hi, margin, pam_5prime,
                              loci=loci if label == "integrated" else None, bad=bad,
                              pam_len=pam_len, cfd_model=ctx.cfd, novar=novar))
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
    no_genotypes = not any(v.samples for v in vcfs_all)
    # the panel: the samplesID files given, else the run's .sampleID.txt, else every
    # sample of the VCF(s); CRISPRme lists no carrier outside it
    panel_from = args.samples
    if not panel_from and os.path.isfile(os.path.join(args.results, ".sampleID.txt")):
        panel_from = [os.path.join(args.results, ".sampleID.txt")]
    if panel_from and not no_genotypes:
        try:
            ctx.panel = read_panel(panel_from)
        except OSError as e:
            sys.stderr.write("ERROR: cannot read the samples list: %s\n" % e)
            return 2
        in_vcf = {s for v in vcfs_all for s in v.samples}
        if not ctx.panel & in_vcf:
            sys.stderr.write("ERROR: no sample of %s is in the VCF(s)\n"
                             % ", ".join(panel_from))
            return 2
        sys.stderr.write("check_indel_rows: the panel is the %d sample(s) of %s (%d of the "
                         "VCF's %d)\n" % (len(ctx.panel), ", ".join(panel_from),
                                          len(ctx.panel & in_vcf), len(in_vcf)))
    elif not no_genotypes:
        sys.stderr.write("check_indel_rows: no samples list (--samples, or the run's "
                         ".sampleID.txt): every sample of the VCF(s) is a panel sample\n")
    if no_genotypes:
        sys.stderr.write("check_indel_rows: the VCF(s) hold no genotypes (a sites-only "
                         "panel): a row's carriers are not checked, only the row itself, "
                         "and completeness has no haplotype to enumerate\n")
    elif any(not v.samples for v in vcfs_all):
        sys.stderr.write("check_indel_rows: %d of the %d VCF(s) hold no genotypes (a "
                         "sites-only dataset beside a genotyped one): a row that lists no "
                         "carrier for a record such a VCF holds is that dataset's, and its "
                         "strings, target and variants are checked, not its carriers\n"
                         % (sum(1 for v in vcfs_all if not v.samples), len(vcfs_all)))
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
    for r in (novar if args.mode in ("soundness", "both") else ()):
        rep.add(True, "alt-without-variant", r,
                "REF/ALT_origin is alt and Variant_info_genome is %r: an ALT row is made "
                "by variants and must name them (the block is not checked further)"
                % r["var"])
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
        reported[ch] = Reported(here, ctx, loci={k: v for k, v in loci.items()
                                                  if k[1] == ch},
                                merge=merge, criterion=criterion, scoring=scoring)
    # the samples SNP rows list per site and variant list (they merge their carriers
    # there), for the SNP-row report's missing carriers
    snp_listed = defaultdict(set)
    if args.snp_rows == "report":
        for r in snp_rows:
            t = r["alt"].replace("-", "")
            snp_listed[(r["guide"], r["chrom"], cluster_of(
                r["start"], r["strand"], len(t), len(r["spacer_col"]), pam_5prime),
                frozenset(x for x in r["var"].split(",") if x and x != "NA"))].update(
                    samples_of(r))

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
    G.update(ctx=ctx, ref=ref, vcfs_for=vcfs_for, reported=reported, popmode=popmode,
             snp_listed=snp_listed)
    done = 0
    # each unit's findings are appended to <out>.partial as it ends, so a long run that is
    # stopped keeps what it found; the sorted TSV replaces it at the end
    partial = None
    if args.out:
        partial = open(args.out + ".partial", "w", newline="")
        pw = csv.DictWriter(partial, fieldnames=COLS, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
        pw.writeheader()

    sites_only = defaultdict(int)   # chromosome -> records only sites-only VCFs hold
    long_ref = {}                   # unindexed VCF -> the longest REF its reads met

    def keep(result):
        found, stats = result
        rep.rows.extend(found)
        sites_only[stats["chrom"]] += stats.get("sites_only", 0)
        for path, n in stats.get("long_ref", {}).items():
            long_ref[path] = max(n, long_ref.get(path, 0))
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
    for ch in sorted(c for c, n in sites_only.items() if n):
        rep.add(False, "limit-no-genotypes", {"table": "completeness", "chrom": ch},
                "%d indel/MNV record(s) of %s%s are held only by sites-only VCF(s) (no "
                "sample columns): no sample carries them, so there is no haplotype to "
                "enumerate and nothing to demand (a sites-only dataset's rows list no "
                "carrier; soundness checks those rows themselves)"
                % (sites_only[ch], ch, "" if lo is None else " in the region"))
    for path in sorted(long_ref):
        rep.add(False, "coverage-unindexed-long-ref", {"table": ""},
                "this checker's own coverage gap, not a CRISPRme limit: %s has no .tbi "
                "and holds a record whose REF spans %d bp; a read without an index starts "
                "%d bp before its region, so a deletion longer than that which starts "
                "further before a region is not read there (its carriers' haplotypes there "
                "lack it). Index the VCF (tabix -p vcf, on a copy or a link in a folder of "
                "your own) and run again" % (path, long_ref[path], Vcf.LONGEST))
    if args.known:
        modes = {"soundness", "completeness"} if args.mode == "both" else {args.mode}
        apply_known(rep, load_known(args.known), chroms, modes)
    # every column in the key: the same input gives the same TSV, whatever order the work
    # units ended in
    rep.rows.sort(key=lambda r: (r["chrom"], str(r["start"]).zfill(12), r["class"],
                                 r["table"], str(r["line"]), r["block"], r["strand"],
                                 r["target"], r["variants"], r["samples"], r["detail"],
                                 r["gating"], r["spacer"]))
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
