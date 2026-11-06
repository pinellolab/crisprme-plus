#!/usr/bin/env python3
"""An independent check of CRISPRme's indel (and MNV) off-target rows on real data.

    check_indel_rows.py --results <results dir> --reference <genome dir or FASTA>
                        --vcf <bgzipped VCF> [--vcf ...] --pam <PAM file>
                        --mm 4 --bDNA 1 --bRNA 1 [--max-total-edits 6]
                        [--chrom chr22 | --region chr22:1-500000]
                        [--mode soundness|completeness|both] [--snp-rows report|off]
                        [--out failures.tsv]

WHAT IT IS. A second opinion on the indel path, built only from what a user has: the
reference FASTA, the VCF(s) the search read, the guide, the PAM, the search budget and
the result folder. It shares no code with CRISPRme's variant or indel machinery
(``PostProcess/indel_snp_cis.py``, ``analisi_indels_NNN.py``, ``overlay_indel_snps*``,
``resultIntegrator``, ``remove_contiguous_samples``, the CRISPRitz enricher, the
dictionaries, the indel log or the fake contigs): it builds every sample's two
haplotypes itself from the VCF, keeps its own map from haplotype offsets to reference
coordinates, and aligns the guide with its own exhaustive dynamic-programming aligner.
The only CRISPRme file it may read is a scoring table, as data (``--cfd-tables``); it
never imports one.

``crisprme.py validate-test`` compares the engine's SITES (chromosome, start, strand) on
the reference and on the SNV-only 1000 Genomes chr22 genome, so an indel row's
coordinate, its aligned strings, its REF/ALT text, its variants, its carriers and the
alternative-alignment table are never checked against an independent truth on real data.
This is that check.

TWO MODES.

  SOUNDNESS (every reported indel/MNV row is real). For each row block of the integrated
  results (highest-CFD, fewest mm+b, highest CRISPR-Bulge) and every alternative
  alignment whose variant list names an indel or MNV:
    * the target (gaps removed, PAM included) is on a listed carrier's haplotype, at the
      reported coordinate, under the convention below;
    * the aligned spacer/target pair is a valid alignment of the guide, and its
      mismatches, bulge types and sizes and their total recount exactly, within budget;
    * the PAM matches the searched motif;
    * the REF text is the reference read from the row's start, gapped as the target;
    * every named variant is a VCF record, forward strand, on that carrier's haplotype
      and inside the target's span, and every variant of that haplotype that changes the
      target is named;
    * every listed carrier carries the target (phased: on one haplotype; unphased: it
      carries every named allele -- CRISPRme's documented PUTATIVE rule);
    * no carrier is missing: a panel sample whose haplotype holds this target at this
      site is listed on this row, or on another reported row of the same site that names
      the same variants (the merge folds those).

  COMPLETENESS (no real indel off-target is missed -- the real-data brute force). For
  every indel and MNV record of the VCF inside the checked region: build each carrier's
  haplotype window, deduplicate identical sequences, enumerate every in-budget alignment
  on both strands with the PAM that holds at least one base the record changes, map it to
  the reference, and require that each (sample, target, site) is reported for that sample.

DEFINITIONS THIS CHECKS (fix/282-indel-alignment; METHODS.md section 4, "SNP+indel
off-targets", and the merge of ``remove_contiguous_samples.py``).

  Coordinate. A row starts at the reference coordinate of its target's FIRST ALIGNED
  BASE -- the leftmost on either strand, 0-based as CRISPRme writes it. A base an
  insertion adds has no reference coordinate: such a row starts at the insertion's
  ANCHOR, the last reference base the added bases follow (VCF POS, 1-based; POS-1
  0-based). Up to 2.8.1 it started at that base's place on the haplotype, 1 to
  len(ALT)-len(REF) bases past the anchor.

  REF text. ``Aligned_protospacer+PAM_REF`` is the reference read from the row's start
  for as many bases as the target holds (an RNA bulge is a gap, so one base fewer per
  RNA bulge), in the target's orientation, with the target's gap columns put back at the
  same column indices. Letter case is not checked (it marks mismatches).

  Bulges. A gap in the aligned spacer is a DNA bulge (an extra genomic base); a gap in
  the aligned target is an RNA bulge. No gap may fall inside the PAM, nor at the first or
  last column.

  Junction rule (which targets are indel targets). A target holds its indel when it holds
  a base the ALT adds or replaces; for a pure deletion, when it holds the bases on BOTH
  sides of the deleted ones -- it spans the junction. A target that holds no base of the
  indel lies wholly in a flank: its sequence is reference (or SNP) sequence and it is not
  an indel row (2.8.1 drops those).

  Merge / superseded. The pipeline clusters rows by guide and chromosome with
  ``Cluster_Position`` within ``--merge`` bp (default 3) of the previous row, chained.
  Inside a cluster, variant rows with the same (cluster position, variant list) -- and,
  for an indel/MNV row, the same aligned target -- merge their carriers into one row; of
  the alignments left naming the same variants with the same carriers, the best by the
  file's own order is kept and the others are dropped. The best row of the cluster goes
  to ``integrated_results.tsv``, the rest to the alternative alignments (a
  population-level run writes none). So a (sample, target) may legitimately be absent
  only when another reported row of the same cluster lists that sample under the same
  variants.

  PUTATIVE rule (unphased). A carrier whose genotypes cannot be put on a haplotype -- an
  unphased call, a missing allele, no phased indel genotype, an MNV (the phased indel
  store keeps only records whose REF and ALT differ in length) -- is listed for the
  strongest target its own alleles allow, and is checked as carrying every allele the row
  names rather than on one haplotype.

OUTPUT. A TSV of every failing row with its reason class and evidence, a per-class
summary on stdout, and exit status 1 on any gating failure. Documented limits of 2.8.1
are counted and listed as their own classes and never gate. SNP-only rows of a per-sample
run can be checked in report-only mode (``--snp-rows report``): population-level SNP rows
are a documented worst-case union and are out of scope.
"""

import argparse
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


def rc(s):
    return s.translate(COMP)[::-1]


def iupac_match(motif_base, base):
    return base.upper() in IUPAC.get(motif_base.upper(), "")


# --------------------------------------------------------------------------- #
# BGZF + tabix: read a VCF region without pysam (an index is optional)
# --------------------------------------------------------------------------- #
class Bgzf(object):
    """Just enough BGZF to start reading at a tabix virtual offset. A BGZF file is a
    series of gzip members, so ``zlib`` with a gzip window decodes each in turn."""

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
        """Yield decoded lines from a virtual offset (coffset << 16 | uoffset). A
        plain-gzip VCF has no index, so it is simply streamed."""
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

    def _member(self):
        """One BGZF block: its gzip header carries a BC extra field holding the block's
        total size, so the next block starts right after it."""
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
    """A bgzipped VCF read by region, with its genotypes. Our own reader: no pysam, no
    CRISPRme code."""

    def __init__(self, path):
        self.path = path
        self.samples = []
        self.contigs = set()
        self.tbi = None
        for ext in (".tbi", ".csi"):
            if os.path.exists(path + ext) and ext == ".tbi":
                try:
                    self.tbi = Tabix(path + ext)
                except Exception:  # noqa: BLE001 - fall back to streaming
                    self.tbi = None
        bg = Bgzf(path)
        try:
            for line in bg.lines(0):
                if line.startswith("##"):
                    if line.startswith("##contig=<ID="):
                        self.contigs.add(line[13:].split(",")[0].split(">")[0])
                    continue
                if line.startswith("#CHROM"):
                    self.samples = line.rstrip("\n").split("\t")[9:]
                break
        finally:
            bg.close()

    def contig(self, chrom):
        """What this VCF calls ``chrom``. CRISPRme names a variant after the reference
        FASTA's contig (``chr22``), and a panel VCF may call it ``22`` -- the chr22
        1000 Genomes VCF of the self-test does."""
        names = [chrom]
        if chrom.startswith("chr"):
            names.append(chrom[3:])
        else:
            names.append("chr" + chrom)
        if self.tbi is not None:
            for n in names:
                if n in self.tbi.refs:
                    return n
            return chrom
        for n in names:
            if n in self.contigs:
                return n
        return chrom

    def records(self, chrom, beg, end):
        """Every record whose REF span overlaps [beg, end) (0-based), as a dict."""
        chrom = self.contig(chrom)
        bg = Bgzf(self.path)
        try:
            if self.tbi is not None:
                spans = self.tbi.chunks(chrom, beg, end)
                if not spans:
                    return
                for cb, ce in spans:
                    for line in bg.lines(cb):
                        if bg.fh.tell() > (ce >> 16) + (1 << 17):
                            break
                        r = self._parse(line, chrom, beg, end)
                        if r == "past":
                            break
                        if r:
                            yield r
            else:
                for line in bg.lines(0):
                    r = self._parse(line, chrom, beg, end)
                    if r == "past":
                        break
                    if r:
                        yield r
        finally:
            bg.close()

    def _parse(self, line, chrom, beg, end):
        """A record as (pos1, REF, [ALT], raw genotype tail). The genotype columns are
        left unsplit: a record of a 4,000-sample panel is 20 kB, and only the records a
        window really needs are split (Window._carriers)."""
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
    def __init__(self, path):
        self.path = path
        self.cache = {}

    def _file(self, chrom):
        if os.path.isdir(self.path):
            for name in ("%s.fa" % chrom, "%s.fasta" % chrom, "%s.fa.gz" % chrom):
                p = os.path.join(self.path, name)
                if os.path.exists(p):
                    return p
            raise FileNotFoundError("no FASTA for %s under %s" % (chrom, self.path))
        return self.path

    def seq(self, chrom):
        if chrom in self.cache:
            return self.cache[chrom]
        p = self._file(chrom)
        op = gzip.open if p.endswith(".gz") else open
        parts = []
        with op(p, "rt") as fh:
            want = False
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    want = (name == chrom) or os.path.isfile(self.path)
                    continue
                if want:
                    parts.append(line.strip())
        s = "".join(parts).upper()
        self.cache = {chrom: s}  # one chromosome at a time
        return s


# --------------------------------------------------------------------------- #
# haplotypes, built here from the VCF
# --------------------------------------------------------------------------- #
def gt_tokens(gt):
    """(tokens, phased). '.' and missing alleles come back as they are."""
    g = gt.split(":")[0]
    if "|" in g:
        return g.split("|"), True
    if "/" in g:
        return g.split("/"), False
    return [g], g not in (".",)


class Window(object):
    """The variants of one region and the haplotypes they make, per sample.

    ``alleles``: one (pos1, ref, alt) triple per VCF ALT, in VCF order.
    ``carry[(pos1, ref, alt)]`` -> {sample: (tokens, phased)} for samples carrying it.
    """

    HOMREF = ("0|0", "0/0", "0", "", ".", ".|.", "./.")

    def __init__(self, ref, chrom, beg, end, vcfs, pad=400, needed=None):
        """``needed``: sorted (lo, hi) reference intervals whose records need their
        genotypes read. A panel record is 4,000 genotype fields, so a soundness pass --
        which only looks at the few hundred bases around each reported row -- reads the
        genotypes of a small fraction of the window's records. None: all of them."""
        self.chrom = chrom
        self.beg = max(0, beg - pad)
        self.end = end + pad
        self.ref = ref
        self.refseq = ref.seq(chrom)
        self.alleles = []
        self.carry = {}
        self.alt_index = {}
        self.samples = []
        self.needed = needed
        self._sorted = None
        seen = set()
        for v in vcfs:
            for s in v.samples:
                if s not in seen:
                    seen.add(s)
                    self.samples.append(s)
            names = v.samples
            for r in v.records(chrom, self.beg, self.end):
                if not self._wanted(r["pos1"] - 1, r["pos1"] - 1 + len(r["ref"])):
                    continue
                carriers = None
                for k, alt in enumerate(r["alts"]):
                    if alt in (".", "*", "<NON_REF>") or not alt:
                        continue
                    key = (r["pos1"], r["ref"], alt)
                    if carriers is None:
                        carriers = self._carriers(r["tail"], names)
                    tok = str(k + 1)
                    who = {s: tp for s, tp in carriers.items() if tok in tp[0]}
                    if key not in self.carry:
                        self.alleles.append(key)
                        self.carry[key] = {}
                        self.alt_index[key] = tok
                    self.carry[key].update(who)
        self.by_sample = defaultdict(list)
        for key in self.alleles:
            tok = self.alt_index[key]
            for s, (toks, phased) in self.carry[key].items():
                self.by_sample[s].append((key, toks, phased, tok))
        for s in self.by_sample:
            self.by_sample[s].sort(key=lambda e: (e[0][0], len(e[0][1]), e[0][2]))

    def _wanted(self, lo, hi):
        if self.needed is None:
            return True
        import bisect
        i = bisect.bisect_right(self.needed, (lo, float("inf")))
        for j in (i - 1, i):
            if 0 <= j < len(self.needed):
                a, b = self.needed[j]
                if lo < b and hi > a:
                    return True
        return False

    @classmethod
    def _carriers(cls, tail, names):
        """{sample: (tokens, phased)} for the samples of this record that are not
        homozygous reference. Nearly every field of a panel VCF is '0|0', so the scan is
        one string compare per sample and only the rest is parsed."""
        out = {}
        if not tail:
            return out
        fields = tail.split("\t")
        homref = cls.HOMREF
        for i, g in enumerate(fields):
            if g in homref or i >= len(names):
                continue
            toks, phased = gt_tokens(g)
            if any(t not in ("0", ".") for t in toks):
                out[names[i]] = (toks, phased)
        return out

    def is_record(self, pos1, ref, alt):
        return (pos1, ref.upper(), alt.upper()) in self.carry

    def carries_allele(self, sample, pos1, ref, alt):
        return sample in self.carry.get((pos1, ref.upper(), alt.upper()), {})

    def slot_alleles(self, sample, h, lo=None, hi=None):
        """The alleles on haplotype ``h`` of ``sample``, in position order: phased calls
        only. Returns None when the sample has an unphased or missing call among them
        (no haplotype can be read). ``lo``/``hi`` restrict it to the records whose REF
        span overlaps [lo, hi), which is what a target's window needs."""
        out = []
        for key, toks, phased, tok in self.by_sample.get(sample, ()):
            if lo is not None and (key[0] - 1 + len(key[1]) <= lo or key[0] - 1 >= hi):
                continue
            if not phased or h >= len(toks):
                return None
            if any(not t.isdigit() for t in toks):
                return None
            if toks[h] == tok:
                out.append(key)
        return out

    def local_keys(self, lo, hi):
        """The records whose REF span overlaps [lo, hi), in position order. A window's
        carriers are then read record by record (a handful of dict lookups) instead of
        walking each sample's whole variant list."""
        import bisect
        if self._sorted is None:
            self._sorted = sorted(self.alleles, key=lambda k: (k[0], len(k[1]), k[2]))
            self._starts = [k[0] - 1 for k in self._sorted]
            self._widest = max([len(k[1]) for k in self._sorted] or [1])
        i = bisect.bisect_left(self._starts, lo - self._widest)
        j = bisect.bisect_left(self._starts, hi)
        return [k for k in self._sorted[i:j]
                if k[0] - 1 + len(k[1]) > lo and k[0] - 1 < hi]

    def slots_in(self, sample, local):
        """(slot 0 alleles, slot 1 alleles) of ``sample`` among ``local``, or None when
        one of its calls there is unphased or missing: no haplotype can be read."""
        out = ([], [])
        for k in local:
            tp = self.carry[k].get(sample)
            if tp is None:
                continue
            toks, phased = tp
            if not phased or len(toks) < 2 or any(not t.isdigit() for t in toks):
                return None
            tok = self.alt_index[k]
            for h in (0, 1):
                if toks[h] == tok:
                    out[h].append(k)
        return out

    def phased(self, sample):
        for _key, _toks, phased, _idx in self.by_sample.get(sample, ()):
            if not phased:
                return False
        return True

    def haplotype(self, sample, h, lo, hi):
        """(sequence, coords) of haplotype ``h`` over the reference window [lo, hi):
        ``coords[i]`` is the 0-based reference coordinate of base ``i``, or None for a
        base an insertion adds. Overlapping records: the first by position wins, the
        others are skipped (a phasing artefact; CRISPRme treats such a slot the same)."""
        keys = self.slot_alleles(sample, h)
        if keys is None:
            return None, None
        return self.apply(keys, lo, hi)

    def apply(self, keys, lo, hi):
        out, coords = [], []
        i = lo
        for pos1, ref, alt in sorted(keys, key=lambda k: (k[0], len(k[1]), k[2])):
            p0 = pos1 - 1
            if p0 + len(ref) <= lo or p0 >= hi:
                continue
            if p0 < i:
                continue  # overlaps one already applied
            out.append(self.refseq[i:p0])
            coords.extend(range(i, p0))
            c = 0
            while c < min(len(ref), len(alt)) and ref[c] == alt[c]:
                c += 1
            out.append(alt[:c])
            coords.extend(range(p0, p0 + c))
            rest = alt[c:]
            out.append(rest)
            if len(alt) == len(ref):          # an MNV: each ALT base stands in a REF place
                coords.extend(range(p0 + c, p0 + len(alt)))
            elif len(alt) < len(ref):         # a deletion: the rest replaces REF bases
                coords.extend(range(p0 + c, p0 + c + len(rest)))
            else:                             # an insertion: the added bases have none
                keep = max(0, len(ref) - c)
                coords.extend(range(p0 + c, p0 + c + keep))
                coords.extend([None] * (len(rest) - keep))
            i = p0 + len(ref)
        out.append(self.refseq[i:hi])
        coords.extend(range(i, hi))
        return "".join(out), coords


def changes_bases(pos1, ref, alt, lo, hi):
    """True iff the record changes the reference bases lo..hi (0-based, inclusive) of a
    haplotype carrying it. The bases REF and ALT share at the start (the VCF anchor) and
    at the end are kept, so they change nothing; bases an insertion adds sit between two
    reference bases, so the insertion changes the span when both sides are in it."""
    a = pos1 - 1
    c = 0
    while c < min(len(ref), len(alt)) and ref[c] == alt[c]:
        c += 1
    t = 0
    while t < min(len(ref), len(alt)) - c and ref[-1 - t] == alt[-1 - t]:
        t += 1
    if len(ref) - t > c:  # bases a+c .. a+len(ref)-t-1 replaced or deleted
        if a + c <= hi and a + len(ref) - t - 1 >= lo:
            return True
    if len(alt) - t > c:  # bases inserted between a+c-1 and a+c
        return lo <= a + c - 1 and a + c <= hi
    return False


def holds_record(seg, pos1, ref, alt):
    """The junction rule read on the target's own bases: ``seg`` is the reference
    coordinate of each base of the target (None for a base an insertion adds), in
    reference order. Exact where a span cannot be: a target holding only the anchor and
    the bases an insertion adds has the span (anchor, anchor)."""
    a = pos1 - 1
    c = 0
    while c < min(len(ref), len(alt)) and ref[c] == alt[c]:
        c += 1
    t = 0
    while t < min(len(ref), len(alt)) - c and ref[-1 - t] == alt[-1 - t]:
        t += 1
    if len(alt) > len(ref):  # an insertion: the target holds a base it adds
        if any(x is None for x in seg):
            reals = [x for x in seg if x is not None]
            return not reals or min(reals) <= a + c - 1 <= max(reals) or True
        return False
    if len(ref) > len(alt):  # a deletion: the target spans the junction
        lo_b, hi_b = a + c - 1, a + len(ref) - t
        keep = [x for x in seg if x is not None]
        if len(alt) - t > c:  # a complex record: it also adds bases
            if any(x is None for x in seg):
                return True
        return any(x == lo_b and y == hi_b for x, y in zip(keep, keep[1:]))
    # an MNV: the target covers a base the ALT replaces
    return any(x is not None and a + c <= x <= a + len(ref) - t - 1 for x in seg)


def holds_indel(first, last, pos1, ref, alt):
    """The junction rule: does a target covering reference bases first..last (0-based,
    inclusive, on the haplotype's own reading) hold a base of this record? A base the ALT
    adds or replaces, or -- for a pure deletion -- the bases on both sides of the deleted
    ones. Anything else lies wholly in a flank."""
    a = pos1 - 1
    c = 0
    while c < min(len(ref), len(alt)) and ref[c] == alt[c]:
        c += 1
    t = 0
    while t < min(len(ref), len(alt)) - c and ref[-1 - t] == alt[-1 - t]:
        t += 1
    added_from = a + c
    if len(ref) == len(alt):  # an MNV replaces bases in place
        return first <= a + len(ref) - t - 1 and last >= added_from
    if len(alt) > len(ref):   # an insertion: the added bases sit between two ref bases
        return first <= added_from - 1 and last >= added_from
    # a deletion: the target must span the junction
    return first < added_from and last >= added_from


# --------------------------------------------------------------------------- #
# our own exhaustive aligner
# --------------------------------------------------------------------------- #
def align_pairs(spacer, target, max_mm, max_dna, max_rna, max_total, pam_len,
                pam_5prime):
    """Every alignment of ``spacer`` (the guide + PAM motif, IUPAC) to ``target`` (plain
    bases) within the budget.

    A gap in the spacer is a DNA bulge (an extra genomic base); a gap in the target is an
    RNA bulge. The PAM columns hold no gap in either string, and the spacer's first and
    last columns are never a gap (a leading DNA bulge or a trailing one is the same
    alignment one base along). A leading RNA bulge is reported by CRISPRme and is kept.
    Yields (aligned spacer, aligned target, mismatches, dna bulges, rna bulges)."""
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
        if j < m:  # a match or mismatch column
            ok = iupac_match(spacer[i], target[j])
            if ok or not in_pam:
                sa.append(spacer[i])
                ta.append(target[j])
                walk(i + 1, j + 1, mm + (0 if ok else 1), dn, rn, sa, ta)
                sa.pop()
                ta.pop()
        # a DNA bulge: an extra genomic base, its column before the guide base i. It may
        # sit immediately before the PAM (i == pam_first on a 3' PAM) but never inside it
        dna_ok = (i <= pam_first) if not pam_5prime else (i > pam_last)
        if dna_ok and dn < max_dna and j < m and len(sa) > 0:
            sa.append("-")
            ta.append(target[j])
            walk(i, j + 1, mm, dn + 1, rn, sa, ta)
            sa.pop()
            ta.pop()
        if not in_pam and rn < max_rna:  # an RNA bulge: the guide base i over a gap
            sa.append(spacer[i])
            ta.append("-")
            walk(i + 1, j, mm, dn, rn + 1, sa, ta)
            sa.pop()
            ta.pop()

    walk(0, 0, 0, 0, 0, [], [])
    return out


def pam_sites(motif, seq, pam_len):
    """Every offset where ``seq`` matches the PAM motif."""
    out = []
    for q in range(0, len(seq) - pam_len + 1):
        seg = seq[q:q + pam_len]
        if all(iupac_match(m, b) for m, b in zip(motif, seg)):
            out.append(q)
    return out


def scan_sequence(spacer, seq, max_mm, max_dna, max_rna, max_total, pam_len, pam_5prime):
    """Every in-budget alignment of the guide on ``seq``, both strands, anchored on the
    PAM (which holds no gap, so it is a fixed window of the sequence). Yields
    (strand, offset of the leftmost aligned base on the + reading of ``seq``,
    aligned spacer, aligned target, mm, dna bulges, rna bulges)."""
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
            name = next(n for n in z.namelist() if n.endswith(".tsv"))
            return io.StringIO(z.read(name).decode())
    return open(path, newline="")


def is_indel_token(tok):
    p = tok.rsplit("_", 3)
    return len(p) == 4 and (len(p[2]) > 1 or len(p[3]) > 1)


def parse_token(tok):
    p = tok.rsplit("_", 3)
    return p[0], int(p[1]), p[2].upper(), p[3].upper()


def read_rows(path, table, chrom=None, lo=None, hi=None):
    """Every (block, fields) of a results table. ``table`` labels the source."""
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
            if chrom and row[col["Chromosome"]] != chrom:
                continue
            for b in BLOCKS:
                try:
                    start = int(row[col["Start_coordinate_(%s)" % b]])
                except (KeyError, ValueError):
                    continue
                if lo is not None and not (lo <= start < hi):
                    continue
                f = {
                    "table": table, "line": line_no, "block": b,
                    "chrom": row[col["Chromosome"]],
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
                    "origin": row[col["REF/ALT_origin_(%s)" % b]],
                    "var": row[col["Variant_info_genome_(%s)" % b]],
                    "samples": row[col["Variant_samples_(%s)" % b]],
                }
                f["cluster"] = start
                if f["origin"] != "alt" or f["var"] in ("", "NA"):
                    continue
                out.append(f)
    return out


def result_tables(resdir):
    """(path, label) of the integrated results and the alternative alignments."""
    out = []
    for p in sorted(glob.glob(os.path.join(resdir, "*integrated_results.tsv"))):
        out.append((p, "integrated"))
    if not out:
        for p in sorted(glob.glob(os.path.join(resdir, "*integrated_results.zip"))):
            out.append((p, "integrated"))
    alt = sorted(glob.glob(os.path.join(resdir, "*all_results_with_alternative_alignments.tsv")))
    if not alt:
        alt = sorted(glob.glob(os.path.join(resdir, "*all_results_with_alternative_alignments.zip")))
    for p in alt:
        out.append((p, "alternative"))
    return out


# --------------------------------------------------------------------------- #
# the checks
# --------------------------------------------------------------------------- #
class Report(object):
    def __init__(self):
        self.rows = []
        self.counts = defaultdict(int)

    def add(self, gating, cls, row, detail):
        self.counts[cls] += 1
        self.rows.append({
            "gating": "yes" if gating else "no", "class": cls,
            "table": row.get("table", ""), "line": row.get("line", ""),
            "block": row.get("block", ""), "chrom": row.get("chrom", ""),
            "start": row.get("start", ""), "strand": row.get("strand", ""),
            "spacer": row.get("crrna", ""), "target": row.get("alt", ""),
            "variants": row.get("var", ""), "samples": row.get("samples", "")[:300],
            "detail": detail,
        })

    def gating(self):
        return sum(1 for r in self.rows if r["gating"] == "yes")

    def write(self, path):
        cols = ["gating", "class", "table", "line", "block", "chrom", "start", "strand",
                "spacer", "target", "variants", "samples", "detail"]
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            for r in self.rows:
                w.writerow(r)


def recount(crrna, target, pam_len, pam_5prime):
    """(mismatches, dna bulges, rna bulges, problem) of an aligned pair as CRISPRme
    writes it: the spacer's gaps are DNA bulges, the target's are RNA bulges, case is
    ignored, and the PAM columns of the spacer are the motif (N/R/...)."""
    if len(crrna) != len(target):
        return 0, 0, 0, "aligned spacer and target differ in length"
    n = len(crrna)
    pam_cols = set(range(0, pam_len)) if pam_5prime else set(range(n - pam_len, n))
    mm = dn = rn = 0
    for i, (g, d) in enumerate(zip(crrna, target)):
        if g == "-" and d == "-":
            return 0, 0, 0, "a column with a gap in both strings"
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


def check_row(row, win, args, rep, pam_motif, pam_len, pam_5prime, cluster_index,
              gating=True, kind="indel"):
    """Every soundness check of one row block. Appends to ``rep``."""
    crrna, alt = row["crrna"], row["alt"]
    stripped = alt.replace("-", "").upper()
    bad = False
    if set(stripped) - set("ACGT"):
        rep.add(gating, "iupac-left", row, "the target holds %s"
                % "".join(sorted(set(stripped) - set("ACGT"))))
        return
    # 1. the alignment recounts
    mm, dn, rn, problem = recount(crrna, alt, pam_len, pam_5prime)
    if problem:
        rep.add(gating, "alignment-invalid", row, problem)
        return
    try:
        r_mm, r_bul, r_tot = int(row["mm"]), int(row["bulges"]), int(row["total"])
    except ValueError:
        rep.add(gating, "alignment-invalid", row, "unreadable counts")
        return
    if (mm, dn + rn) != (r_mm, r_bul) or r_tot != r_mm + r_bul:
        rep.add(gating, "counts-mismatch", row,
                "reported %d mm + %d bulges = %d; the alignment has %d mm, %d DNA and "
                "%d RNA bulges" % (r_mm, r_bul, r_tot, mm, dn, rn))
        bad = True
    want = ",".join([t for t, c in (("RNA", rn), ("DNA", dn)) if c] ) or "X"
    if row["btype"] not in (want, "X", "") and not (dn == rn == 0):
        rep.add(gating, "bulge-type-mismatch", row,
                "reported %s, the alignment has %d DNA and %d RNA bulges"
                % (row["btype"], dn, rn))
        bad = True
    if mm > args.mm or dn > args.bdna or rn > args.brna or (
            args.max_total_edits is not None and mm + dn + rn > args.max_total_edits):
        rep.add(gating, "budget-exceeded", row,
                "%d mm, %d DNA, %d RNA bulges against %d/%d/%d, cap %s"
                % (mm, dn, rn, args.mm, args.bdna, args.brna, args.max_total_edits))
        bad = True
    # 2. the PAM
    pam_seq = stripped[:pam_len] if pam_5prime else stripped[-pam_len:]
    if not all(iupac_match(m, b) for m, b in zip(pam_motif, pam_seq)):
        rep.add(gating, "pam-mismatch", row, "PAM %s against the motif %s"
                % (pam_seq, pam_motif))
        bad = True
    if row["pam"] and row["pam"] != "NA" and row["pam"].upper() != pam_seq:
        rep.add(gating, "pam-column-mismatch", row,
                "the PAM column says %s, the target has %s" % (row["pam"], pam_seq))
        bad = True
    # 3. the REF text is the reference read from the row's start, gapped as the target
    refseq = win.refseq
    span = refseq[row["start"]:row["start"] + len(stripped)]
    if len(span) == len(stripped):
        text = list(span if row["strand"] == "+" else rc(span))
        for i, c in enumerate(alt):
            if c == "-":
                text.insert(i, "-")
        if row["ref"].upper() != "".join(text).upper():
            rep.add(gating, "ref-text-mismatch", row,
                    "REF %s, the reference at %d gapped as the target is %s"
                    % (row["ref"], row["start"], "".join(text)))
            bad = True
    # 4. the variants
    toks = [t for t in row["var"].split(",") if t and t != "NA"]
    named = []
    for t in toks:
        try:
            c, p, r, a = parse_token(t)
        except (ValueError, IndexError):
            rep.add(gating, "variant-unparseable", row, t)
            bad = True
            continue
        named.append((p, r, a))
        if not win.is_record(p, r, a):
            rep.add(gating, "variant-not-in-vcf", row,
                    "%s is not a record of the VCF (forward strand)" % t)
            bad = True
    samples = [s for s in row["samples"].split(",") if s and s not in ("NA", "n")]
    if not samples:
        rep.add(gating, "no-carrier-listed", row, "the row lists no carrier")
        return
    # 5. the target is on a listed carrier's haplotype at the reported coordinate
    lo = max(0, row["start"] - 2 * len(stripped) - 400)
    hi = min(len(refseq), row["start"] + 2 * len(stripped) + 400)
    probe = stripped if row["strand"] == "+" else rc(stripped)
    placed = False
    carrier_ok = []
    for s in samples:
        hit = find_on_sample(win, s, probe, row["start"], lo, hi, named)
        if hit:
            carrier_ok.append((s, hit))
            placed = True
    if not placed:
        starts = set()
        for s in samples[:40]:
            starts |= all_starts(win, s, probe, lo, hi)
        rep.add(gating, "start-not-on-carrier-haplotype" if starts
                else "target-not-on-carrier-haplotype", row,
                "no listed carrier has the target at %d; its first base is at %s on "
                "their haplotypes" % (row["start"], sorted(starts)[:8] or "no place"))
        bad = True
    # 6. each listed carrier really has it. Only where the row is placed at all: when no
    #    carrier has the target at the reported start, every carrier fails for that one
    #    reason and the per-carrier lines would only restate it.
    for s in (samples if placed else ()):
        if win.phased(s):
            if not any(c == s for c, _h in carrier_ok):
                if find_on_sample(win, s, probe, row["start"], lo, hi, None):
                    rep.add(gating, "carrier-without-variants", row,
                            "%s has the sequence but not with every variant the row "
                            "names on one haplotype" % s)
                else:
                    rep.add(gating, "carrier-without-target", row,
                            "%s is listed and no haplotype of it has the target there" % s)
                bad = True
        else:  # the documented PUTATIVE rule: it must carry every named allele
            for p, r, a in named:
                if not win.carries_allele(s, p, r, a):
                    rep.add(gating, "carrier-without-allele", row,
                            "%s is listed (unphased) without %s_%d_%s_%s"
                            % (s, row["chrom"], p, r, a))
                    bad = True
    # 7. the named variants are inside the target's span, and none that changes it is
    #    left out
    if carrier_ok:
        s, (h, coords) = carrier_ok[0]
        reals = [c for c in coords if c is not None]
        if reals:
            flo, fhi = min(reals), max(reals)
            for p, r, a in named:
                if not holds_record(list(coords), p, r, a) and not changes_bases(
                        p, r, a, flo, fhi):
                    rep.add(gating, "variant-outside-target", row,
                            "%s_%d_%s_%s changes no base of the target (%d..%d)"
                            % (row["chrom"], p, r, a, flo, fhi))
                    bad = True
            keys = win.slot_alleles(s, h) or []
            for p, r, a in keys:
                if (p, r, a) in named:
                    continue
                if changes_bases(p, r, a, flo, fhi):
                    rep.add(gating, "variant-missing", row,
                            "%s_%d_%s_%s is on %s's haplotype and changes the target, "
                            "but the row does not name it" % (row["chrom"], p, r, a, s))
                    bad = True
    # 8. no carrier missing: a panel sample whose haplotype holds this target here and
    #    is on no reported row of the same cluster under the same variants
    if args.check_carriers and carrier_ok:
        listed = set(samples)
        for other in cluster_index.get((row["chrom"], row["strand"], stripped), ()):
            listed |= other
        # a sample whose haplotype holds this target must carry the row's indel (the
        # target holds a base of it), so only the record's own carriers can be missing
        pool = set()
        for p, r, a in named:
            if len(r) > 1 or len(a) > 1:
                pool |= set(win.carry.get((p, r, a), ()))
        missing = []
        for s in sorted(pool):
            if s in listed or not win.phased(s):
                continue
            if find_on_sample(win, s, probe, row["start"], lo, hi, named):
                missing.append(s)
        if missing:
            rep.add(gating, "carrier-missing", row,
                    "%d panel sample(s) hold this target here and are on no row of the "
                    "site: %s" % (len(missing), ",".join(sorted(missing)[:10])))
            bad = True
    return bad


def all_starts(win, sample, probe, lo, hi):
    """Where the target's first base sits on the sample's haplotypes, under the
    convention (the anchor for a base an insertion adds)."""
    out = set()
    for h in (0, 1):
        hap, coords = win.haplotype(sample, h, lo, hi)
        if hap is None:
            continue
        k = hap.find(probe)
        while k >= 0:
            out.add(first_coord(coords, k))
            k = hap.find(probe, k + 1)
    return {c for c in out if c is not None}


def first_coord(coords, k):
    """The reference coordinate a target starting at haplotype offset ``k`` reports: the
    base's own, or -- for a base an insertion adds -- the anchor, the last reference base
    before it."""
    if coords[k] is not None:
        return coords[k]
    for q in range(k - 1, -1, -1):
        if coords[q] is not None:
            return coords[q]
    return None


def find_on_sample(win, sample, probe, start, lo, hi, named):
    """(slot, coordinates of the target's bases) where ``sample`` holds ``probe`` with
    its first base reported at ``start`` -- and, when ``named`` is given, on a haplotype
    carrying every named variant."""
    for h in (0, 1):
        keys = win.slot_alleles(sample, h)
        if keys is None:
            continue
        if named and not all((p, r, a) in keys for p, r, a in named):
            continue
        hap, coords = win.apply(keys, lo, hi)
        k = hap.find(probe)
        while k >= 0:
            if first_coord(coords, k) == start:
                return h, coords[k:k + len(probe)]
            k = hap.find(probe, k + 1)
    return None


# --------------------------------------------------------------------------- #
# completeness
# --------------------------------------------------------------------------- #
def completeness(win, args, rep, pam_motif, pam_len, pam_5prime, reported, spacer,
                 region):
    """The real-data brute force: for every indel / MNV record in the region, every
    in-budget alignment on each carrier's haplotype that holds a base the record changes
    must be reported for that carrier."""
    lo_r, hi_r = region
    n = len(spacer)
    reach = n + max(args.bdna, args.brna) + 2
    limits = defaultdict(list)
    claimed = set()
    local_cache = {}
    for key in win.alleles:
        pos1, ref, alt = key
        if len(ref) == 1 and len(alt) == 1:
            continue  # a SNV: not an indel-path record
        if not (lo_r <= pos1 - 1 < hi_r):
            continue
        carriers = win.carry[key]
        if not carriers:
            continue
        lo = max(0, pos1 - 1 - reach - len(ref))
        hi = min(len(win.refseq), pos1 - 1 + len(ref) + reach)
        if lo == 0 or hi == len(win.refseq):
            limits["limit-contig-end"].append("%d_%s_%s" % (pos1, ref, alt))
            continue
        local = local_cache.get((lo, hi))
        if local is None:
            local = win.local_keys(lo, hi)
            local_cache.clear()
            local_cache[(lo, hi)] = local
        groups = {}
        for s in carriers:
            slots = win.slots_in(s, local)
            if slots is None:
                limits["limit-unphased-carrier"].append("%s@%d" % (s, pos1))
                continue
            for h in (0, 1):
                if key in slots[h]:
                    groups.setdefault(tuple(slots[h]), []).append((s, h))
        for keys, who in groups.items():
            hap, coords = win.apply(keys, lo, hi)
            for strand, off, sa, ta, mm, dn, rn in scan_sequence(
                    spacer, hap, args.mm, args.bdna, args.brna,
                    args.max_total_edits if args.max_total_edits is not None else 99,
                    pam_len, pam_5prime):
                width = len(ta.replace("-", ""))
                seg = list(coords[off:off + width])
                if len(seg) < width:
                    continue
                if not holds_record(seg, pos1, ref, alt):
                    continue  # a flank target: the reference / SNP search reports it
                start = first_coord(list(coords), off)
                if start is None:
                    continue
                target = hap[off:off + width]
                probe = (target if strand == "+" else rc(target)).upper()
                total = mm + dn + rn
                # a second indel (or MNV) of this haplotype inside the target: the
                # fake-indel window holds one indel, so 2.8.1 does not generate this
                # candidate at all (CHANGELOG 2.8.1, Known limits; METHODS section 8)
                reals = [x for x in seg if x is not None]
                flo = min(reals) if reals else start
                fhi = max(reals) if reals else start
                second = [k for k in keys
                          if k != (pos1, ref, alt) and (len(k[1]) > 1 or len(k[2]) > 1)
                          and changes_bases(k[0], k[1], k[2], flo, fhi)]
                if second:
                    limits["limit-two-indels"].append(
                        "%s:%d %s with %s" % (win.chrom, start,
                                              "%d_%s_%s" % (pos1, ref, alt),
                                              ",".join("%d_%s_%s" % k for k in second)))
                    continue
                for s, _h in who:
                    if (s, strand, start, probe) in reported["exact"]:
                        continue
                    if (s, strand, start, probe) in claimed:
                        continue  # the same alignment, reached from another record
                    claimed.add((s, strand, start, probe))
                    # the merge keeps one alignment per site, variants and carriers: a
                    # reported row of the same site, for this sample, that is at least
                    # as good supersedes this one
                    near = [reported["best"].get((s, strand, st))
                            for st in range(start - args.merge, start + args.merge + 1)]
                    near = [t for t in near if t is not None]
                    if near and min(near) <= total:
                        continue
                    cls = "missing-better-alignment" if near else "missing-row"
                    detail = ("%s has this target at %d (%d mm, %d DNA, %d RNA bulges, "
                              "%d edits); %s" % (
                                  s, start, mm, dn, rn, total,
                                  ("the best row of the site for it has %d edits"
                                   % min(near)) if near
                                  else "no row of either table reports it"))
                    rep.add(True, cls, {
                        "table": "completeness", "chrom": win.chrom, "start": start,
                        "strand": strand, "crrna": sa, "alt": ta,
                        "var": "%s_%d_%s_%s" % (win.chrom, pos1, ref, alt),
                        "samples": s}, detail)
    for cls, items in limits.items():
        rep.counts[cls] += len(items)
        rep.rows.append({
            "gating": "no", "class": cls, "table": "completeness", "line": "",
            "block": "", "chrom": win.chrom, "start": "", "strand": "", "spacer": "",
            "target": "", "variants": "", "samples": "",
            "detail": "%d case(s); first: %s" % (len(items), "; ".join(items[:5])),
        })


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def read_pam(path):
    """(motif, length, 5-prime?) from a CRISPRme PAM file: a padded motif and its PAM
    length, e.g. 'NNNNNNNNNNNNNNNNNNNNNRG 3'."""
    with open(path) as fh:
        text = fh.read().split()
    motif, k = text[0].upper(), int(text[1])
    if k < 0:
        return motif[:-k], -k, True
    return motif[-k:], k, False


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--results", required=True, help="a CRISPRme result folder")
    p.add_argument("--reference", required=True, help="the reference genome dir or FASTA")
    p.add_argument("--vcf", action="append", default=[], help="a bgzipped VCF (repeatable)")
    p.add_argument("--vcf-dir", action="append", default=[], help="a folder of VCFs")
    p.add_argument("--pam", required=True, help="the PAM file the search used")
    p.add_argument("--guide", default=None, help="the spacer (default: the table's)")
    p.add_argument("--mm", type=int, required=True)
    p.add_argument("--bDNA", dest="bdna", type=int, required=True)
    p.add_argument("--bRNA", dest="brna", type=int, required=True)
    p.add_argument("--max-total-edits", type=int, default=6)
    p.add_argument("--merge", type=int, default=3, help="the search's --merge threshold")
    p.add_argument("--chrom", default=None)
    p.add_argument("--region", default=None, help="chrom:start-end (1-based, inclusive)")
    p.add_argument("--mode", default="both",
                   choices=("soundness", "completeness", "both"))
    p.add_argument("--snp-rows", default="off", choices=("off", "report"),
                   help="also check SNP-only rows of a per-sample run, report-only")
    p.add_argument("--check-carriers", action="store_true", default=True)
    p.add_argument("--no-check-carriers", dest="check_carriers", action="store_false")
    p.add_argument("--out", default=None, help="the TSV of failures")
    p.add_argument("--tile", type=int, default=50000,
                   help="how many bases of a chromosome one VCF read covers")
    p.add_argument("--progress", type=int, default=0, help="a line every N tiles")
    p.add_argument("--max-rows", type=int, default=0, help="stop after N row blocks (debug)")
    args = p.parse_args(argv)

    if args.max_total_edits is not None and args.max_total_edits < 0:
        args.max_total_edits = None
    pam_motif, pam_len, pam_5prime = read_pam(args.pam)
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

    ref = Reference(args.reference)
    rep = Report()
    tables = result_tables(args.results)
    if not tables:
        sys.stderr.write("ERROR: no results table under %s\n" % args.results)
        return 2
    rows = []
    for path, label in tables:
        rows.extend(read_rows(path, label, chrom, lo, hi))
    indel_rows = [r for r in rows
                  if any(is_indel_token(t) for t in r["var"].split(","))]
    snp_rows = [r for r in rows
                if not any(is_indel_token(t) for t in r["var"].split(","))]
    sys.stderr.write("check_indel_rows: %d indel/MNV row block(s), %d SNP-only, from %s\n"
                     % (len(indel_rows), len(snp_rows),
                        ", ".join(os.path.basename(p) for p, _l in tables)))
    if args.max_rows:
        indel_rows = indel_rows[:args.max_rows]

    spacer = args.guide
    if spacer is None:
        # the table's Spacer+PAM column is the guide padded with N over the PAM
        for r in rows:
            s = r["spacer_col"]
            spacer = (s[:-pam_len] if not pam_5prime else s[pam_len:])
            break
    if spacer is None:
        if args.mode == "soundness":
            # a run with no variant row at all (a reference-only or tightly capped
            # search): there is nothing to check, and no guide to read from a row
            print("check_indel_rows: no variant row in the results; nothing to check")
            print("check_indel_rows: PASS")
            return 0
        sys.stderr.write("ERROR: completeness needs --guide when the results hold no "
                         "variant row to read one from\n")
        return 2
    full = (pam_motif + spacer) if pam_5prime else (spacer + pam_motif)

    # group the work by chromosome, then by tile: one VCF read and one haplotype build
    # serve every row of a tile, and no window is ever the whole chromosome
    by_chrom = defaultdict(list)
    for r in indel_rows + (snp_rows if args.snp_rows == "report" else []):
        by_chrom[r["chrom"]].append(r)
    vcfs_all = [Vcf(p) for p in vcf_paths]
    tile = args.tile
    if not any(v.tbi is not None for v in vcfs_all):
        # without a tabix index a region read is a full stream, so one tile per
        # chromosome: the VCF is read once instead of once per tile
        tile = 1 << 30
        sys.stderr.write("check_indel_rows: no .tbi beside the VCF(s); reading each "
                         "chromosome in one pass\n")
    pad = 2 * len(full) + 400

    chroms = sorted(set(by_chrom))
    if args.mode in ("completeness", "both"):
        if chrom:
            chroms = sorted(set(chroms) | {chrom})
        elif not chroms:
            sys.stderr.write("ERROR: completeness needs --chrom or --region\n")
            return 2
    for ch in chroms:
        vcfs = [v for v in vcfs_all
                if (v.tbi is None or v.contig(ch) in v.tbi.refs)]
        if not vcfs:
            continue
        rows_here = by_chrom.get(ch, [])
        # which (sample, strand, start, target) each site reports, for the
        # merge/superseded rule and for completeness
        reported = {"exact": set(), "best": {}}
        cluster_index = defaultdict(list)
        for r in rows_here:
            t = r["alt"].replace("-", "").upper()
            carriers = {s for s in r["samples"].split(",") if s and s not in ("NA", "n")}
            try:
                total = int(r["total"])
            except ValueError:
                total = 99
            for s in carriers:
                reported["exact"].add((s, r["strand"], r["start"], t))
                k = (s, r["strand"], r["start"])
                if total < reported["best"].get(k, 99):
                    reported["best"][k] = total
            cluster_index[(r["chrom"], r["strand"], t)].append(carriers)
        tiles = {}
        for r in rows_here:
            tiles.setdefault(r["start"] // tile, []).append(r)
        if args.mode in ("completeness", "both"):
            clo = lo if lo is not None else 0
            chi = hi if hi is not None else len(ref.seq(ch))
            for t in range(clo // tile, chi // tile + 1):
                tiles.setdefault(t, [])
        done = 0
        for t in sorted(tiles):
            tlo, thi = t * tile, (t + 1) * tile
            needed = None
            if args.mode == "soundness":
                spans = sorted((max(0, r["start"] - pad), r["start"] + pad)
                               for r in tiles[t])
                needed = []
                for a, b in spans:
                    if needed and a <= needed[-1][1]:
                        needed[-1] = (needed[-1][0], max(needed[-1][1], b))
                    else:
                        needed.append((a, b))
            # a tile that is the whole chromosome (no tabix index) still only needs the
            # records of --region, so clamp the VCF window to it
            wlo = max(0, tlo - pad, (lo - pad) if lo is not None else 0)
            whi = min(thi + pad, (hi + pad) if hi is not None else thi + pad)
            win = Window(ref, ch, wlo, whi, vcfs, pad=0, needed=needed)
            if args.mode in ("soundness", "both"):
                for r in sorted(tiles[t], key=lambda x: x["start"]):
                    kind = "indel" if any(is_indel_token(x)
                                          for x in r["var"].split(",")) else "snp"
                    check_row(r, win, args, rep, pam_motif, pam_len, pam_5prime,
                              cluster_index, gating=(kind == "indel"), kind=kind)
            if args.mode in ("completeness", "both"):
                rlo = max(tlo, lo) if lo is not None else tlo
                rhi = min(thi, hi) if hi is not None else thi
                if rlo < rhi:
                    completeness(win, args, rep, pam_motif, pam_len, pam_5prime,
                                 reported, full, (rlo, rhi))
            done += 1
            if args.progress and done % args.progress == 0:
                sys.stderr.write("  %s tile %d/%d, %d finding(s)\n"
                                 % (ch, done, len(tiles), len(rep.rows)))
                sys.stderr.flush()

    if args.out:
        rep.write(args.out)
    print("check_indel_rows: %d finding(s), %d gating" % (len(rep.rows), rep.gating()))
    for cls in sorted(rep.counts):
        print("  %-34s %d" % (cls, rep.counts[cls]))
    if not rep.rows:
        print("check_indel_rows: PASS")
    return 1 if rep.gating() else 0


if __name__ == "__main__":
    sys.exit(main())
