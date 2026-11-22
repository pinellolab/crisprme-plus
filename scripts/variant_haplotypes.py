#!/usr/bin/env python3
"""Reading a VCF and a reference FASTA, and writing the haplotypes they describe.

This is the reconstruction core that ``scripts/check_indel_rows.py`` was built on and
hardened with (the SBDSP1 real-data investigation): its own BGZF + tabix reader, its own
reference FASTA reader, VCF 4.x's definition of what an ALT allele does to a haplotype,
and the map from a haplotype's bases back to reference coordinates that makes an indel's
shift readable. It lives in a file of its own so the checker and
``crisprme.py extract-flanks`` (``PostProcess/extract_flanks.py``) write a sample's
haplotype by the SAME code: two readers of one VCF that disagree would be two different
answers to "what does this individual's genome hold here?", and one of them would send a
user's validation oligos to the wrong sequence.

It is the checker's own code, moved whole and unchanged. It keeps the two properties that
make the checker an independent second opinion, and ANY CHANGE HERE MUST KEEP THEM:

  * STANDARD LIBRARY ONLY. No pysam, no pandas, no third-party anything: a user with a
    Python and a VCF can run it.
  * NO CRISPRme VARIANT OR INDEL MACHINERY. Nothing here imports or copies
    ``PostProcess/indel_snp_cis.py``, ``site_alleles.py``, ``analisi_indels_NNN.py``,
    ``overlay_indel_snps*``, ``resultIntegrator``, ``remove_contiguous_samples``, the
    CRISPRitz enricher, the dictionaries, the indel log or the fake contigs. The truth
    here is the VCF as VCF 4.x defines it, not as CRISPRme's search read it.

What it does NOT hold: the checker's view of how a given CRISPRme/CRISPRitz release read a
VCF (``Window``, ``Semantics``, the dictionary and shadowing rules, the FILTER readings of
the several stores), its aligner, its CFD model and its reporting. Those model CRISPRme's
behaviour and belong to the checker alone. A caller that wants the plain VCF truth --
every ALT a genotype puts on a slot, in place of its REF bases -- builds it from
``Allele`` and ``make_hap`` here.

Coordinates are 0-based half-open internally; a VCF POS is 1-based, as the file writes it.
"""

import bisect
import gzip
import os
import struct
import zlib

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
# the FILTER values CRISPRitz's enricher reads (Enrichment/enricher.py and enricher.cpp:
# `if line[6] not in ('PASS', '.')`, the line skipped; PostProcess/validate_inputs.py's
# ENRICHER_PASS_VALUES): no line outside them is in the enriched genome, its dictionaries,
# its registry and genotype store, its indel log or its fake-indel windows
PASS_FILTERS = ("PASS", ".")

def rc(s):
    return s.translate(COMP)[::-1]


def comp(b):
    return b.translate(COMP)

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

    def holds(self, chrom):
        """Whether this VCF holds a record of ``chrom`` (a plain-gzip file that cannot
        tell from its first line is read until it finds one)."""
        known = self.covers(chrom)
        if known is not None:
            return known
        name = self.contig(chrom)
        with gzip.open(self.path, "rt") as fh:
            return any(line.startswith(name + "\t") for line in fh)

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
        genotypes. Its FILTER is kept as written (Window reads it: PASS_FILTERS). REF and
        ALT are read in upper case, as VCF 4.x reads them, as CRISPRitz 2.10.0's enricher
        reads them, and as CRISPRme 2.10.0 names an indel or MNV record; the ALT column as
        written is kept too (``alt_text``), for an SNV's ALT written in lower case, which
        CRISPRitz 2.9.1 and earlier do not code (Window)."""
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
                "alts": [a.upper() for a in f[3].split(",")], "alt_text": f[3],
                "filter": f[5], "info": f[6], "format": f[7] if len(f) > 7 else "",
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
# ALT alleles and genotypes, as VCF 4.x defines them
# --------------------------------------------------------------------------- #
class Allele(object):
    """One ALT allele of a VCF line. On a haplotype that carries it, ALT stands in place
    of the REF bases [lo, hi) (0-based). The bases REF and ALT share at either end stay
    where they are, so the allele's edit is what lies between them: it replaces the
    reference bases [edit_lo, edit_hi) with ``body``."""

    __slots__ = ("key", "lo", "hi", "line", "number", "edit_lo", "edit_hi", "body",
                 "kind", "_bases", "searched", "filtered")

    def __init__(self, pos1, ref, alt, line, number):
        self.searched = True    # False: an indel CRISPRme's search skips (FILTER) and only
                                # its phased indel store holds (Window)
        self.filtered = None    # the FILTER of a line outside PASS_FILTERS that is kept
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


def meets_ref(a, al):
    """Whether allele ``a``'s REF span shares a base with record ``al``'s: on one slot the
    two would describe those bases twice (no sequence there)."""
    return a.lo < al.hi and a.hi > al.lo


# --------------------------------------------------------------------------- #
# genotypes: which slot of a sample holds which allele
# --------------------------------------------------------------------------- #
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


def merge_call(a, b):
    """One sample's call of one SNV ALT (1, every other called allele 0) that two lines
    give, ``a`` the earlier line's: the one call CRISPRitz 2.10.0's enricher writes into the
    allele's one dictionary entry, which CRISPRme reads (CRISPRitz's CHANGELOG and its
    docs/MULTIALLELIC_ENRICHMENT.md section 2.3). Phased alike -- as many alleles, both
    phased ('|', or haploid) --, a slot holds the ALT where either call puts it there, else
    a called allele (0) where either has one; otherwise, as unphased slots cannot be
    matched up, the call with more copies of the ALT, the earlier one on a tie, as it is."""
    (xa, pa), (xb, pb) = a, b
    if pa and pb and len(xa) == len(xb):
        return tuple(1 if 1 in (x, y) else (0 if (x, y) != (None, None) else None)
                     for x, y in zip(xa, xb)), True
    return b if xb.count(1) > xa.count(1) else a


def merge_calls(into, more):
    """Merge ``more`` ({sample: call} of a later line of one SNV REF>ALT) into ``into``."""
    for s, c in more.items():
        into[s] = merge_call(into[s], c) if s in into else c


def split_call(call, number):
    """A call of a line of several ALTs as a line of one ALT writes it for its allele
    ``number`` (CRISPRitz's enricher's form of a dictionary entry it merges): that allele
    1, every other called allele 0, a missing one missing."""
    alle, phased = call
    return tuple(None if x is None else (1 if x == number else 0) for x in alle), phased


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
