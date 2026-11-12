#!/usr/bin/env python3
"""Unit tests for ``scripts/check_indel_rows.py``, the independent real-data check of
CRISPRme's indel (and MNV) off-target rows.

The checker must CATCH a wrong row and PASS a right one, so the suite is built both
ways: hand-built true rows that must raise nothing gating, and mutation tests that
corrupt one field of such a row -- the start shifted by 1 and by an insertion's offset,
one target base changed, a carrier swapped or added or dropped, a named variant dropped,
a cross-haplotype SNP added, the aligned spacer rewritten, a phantom of the shape Manuel
Tognon reported (an insertion's only carrier listed under SNP alleles only other
individuals carry), an MNV row and an unphased row each with a phantom base, an unphased
carrier listed on a weaker target than its strongest, an indel block labelled REF origin,
a population-level cluster whose row is worse than a carrier's target, a merge that keeps
the worse of two alignments of one site or of one target (fewest block, CFD block, a CFD
only the three printed decimals tell apart: three real rows of chr22; one worse only by the
merge's next sorting criterion), a fold across strands, a carrier listed for a sequence an
MNV on its slot changes, a carrier outside the samples list, an ALT
block that lost its variants, a cluster dropped beside another cluster's row, a row
placed past an insertion's anchor, two targets inside one insertion folded into one, a
long deletion's far junction, a second indel that keeps the target's bases, another
guide's row, a SNP named for a carrier whose slot overlaps the record and holds it nowhere,
a slot read without its phasing artefact that lacks the target, a listed carrier whose
slot has no one sequence and that METHODS' rules take off, a carrier left off whose
slot's overlapping alleles change no base of the target, a population-level cluster
joined only through a target that needs two indel records, a run read in another mode
than the one it ran -- and assert that each one is reported, with its class. The truth is
the VCF: a slot whose alleles overlap gives no sequence, so it is read without the
alleles METHODS section 4 reads it without (the one reading that holds the record), and
on every reading of two other alleles that overlap each other where they can change the
target, judged where all readings agree, and otherwise judged by CRISPRme's rules alone
(a listed carrier they take off fails) or counted (limit-overlapping-records), its named
variants always required on a slot that can be the record's; CRISPRme's
rules (a split multi-allelic site, two indels in one target, a second indel
read up to the target's last base, the PUTATIVE carriers of an MNV) only turn a
disagreement into a non-gating limit, and a tie their PUTATIVE rule leaves is a limit line,
not silence. The PUTATIVE truth is the rule on the whole alignment's CFD (a product of 0
ties every choice), after what a record's known slot says. The driver must never pass
silently (a chromosome no VCF covers, one FASTA of several contigs, an indexed VCF read
once per record, a sites-only VCF alone or beside a genotyped one, a deletion longer than
an unindexed read looks back, a short table line, a window holding N), must refuse a
configuration it does not model (two genotyped datasets that share samples, inputs other
than those the run records it searched with, a search that says it read no genotype
store: its log, or a --per-sample search run population-level), must read a VCF's FILTER as CRISPRme's
stores do (a record outside PASS and '.' is no record of the search, but a genotyped
dataset's indels stay on their carriers' haplotypes, as the phased indel store keeps
them), must read a VCF line of several ALTs as CRISPRme's readers do (an alt whose number
they cannot tie, beside a '*' or an indel of its line too, is PUTATIVE; a listing of an SNV
merged into its first, from a line of several ALTs too), must not read a run's allele coding
off a '*' line, and must read a target at its fake-indel window's edge by its site, and the
checker must hold no CRISPRme code.

Everything is built here: a small reference, a phased VCF (plain gzip, or BGZF with a
tabix index written by this file's own writer), and rows in the column layout of
``integrated_results.tsv``. No CRISPRme module is imported, as the checker imports none.
"""

import contextlib
import gzip
import io
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
import zlib
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                                "scripts"))
import check_indel_rows as C  # noqa: E402

GUIDE = "ACTGAAATCTGTAAGCAGGC"
SAMPLES = ["S%02d" % i for i in range(1, 9)]
PRE = "GGCATTACGGTTACCATGCAAT"
POST = "ACGTTGCAATTGGCCAATTGCA"
LONG = "TTGACCGGTATCAGCTAGCCTTAGGCATCGATGCTAGCATGACTGATCGTACGTAGCTAGT"


def other(b):
    return "A" if b != "A" else "C"


def al(pos1, ref, alt):
    return C.Allele(pos1, ref, alt, 0, 1)


# --------------------------------------------------------------------------- #
# writers: FASTA, VCF (plain gzip or BGZF + tabix), results tables
# --------------------------------------------------------------------------- #
def write_fasta(path, seq, name="chr21"):
    with open(path, "w") as fh:
        fh.write(">%s\n" % name)
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def vcf_text(chrom, records, samples=SAMPLES):
    body = ["##fileformat=VCFv4.2", "##contig=<ID=%s>" % chrom,
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for rec in records:
        pos1, ref, alt, gts = rec[:4]
        filt = rec[4] if len(rec) > 4 else "PASS"
        row = [chrom, str(pos1), ".", ref, alt, ".", filt, "AC=1", "GT"]
        row += [gts.get(s, "0|0") for s in samples]
        body.append("\t".join(row))
    return "\n".join(body) + "\n"


def write_vcf(path, chrom, records, samples=SAMPLES):
    """``records``: (pos1, ref, alt, {sample: gt}[, FILTER]), written as a plain gzip
    stream (FILTER PASS where a record gives none)."""
    with gzip.open(path, "wt") as fh:
        fh.write(vcf_text(chrom, records, samples))


def write_sites_only_vcf(path, chrom, records):
    """A sites-only VCF (eight columns: no FORMAT, no sample), as hg38_mega's: the
    records' (pos1, ref, alt[, FILTER]) (PASS where a record gives none; a fourth item that
    is not a string is ignored)."""
    body = ["##fileformat=VCFv4.2", "##contig=<ID=%s>" % chrom,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for r in records:
        filt = r[3] if len(r) > 3 and isinstance(r[3], str) else "PASS"
        body.append("\t".join([chrom, str(r[0]), ".", r[1], r[2], ".", filt, "AF=0.1"]))
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(body) + "\n")


def _bgzf_block(data):
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = comp.compress(data) + comp.flush()
    bsize = 18 + len(payload) + 8
    head = struct.pack("<4BI2BH2BHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize - 1)
    return head + payload + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))


def _reg2bin(beg, end):
    end -= 1
    for sh, base in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> sh == end >> sh:
            return base + (beg >> sh)
    return 0


def write_bgzf_tabix(path, text, block=200):
    """``text`` as BGZF in small blocks (so a region read has to seek) and a tabix index
    beside it, both written here with the standard library."""
    raw = text.encode()
    out = bytearray()
    starts = []  # (uncompressed offset of the block, compressed offset)
    for i in range(0, len(raw), block):
        starts.append((i, len(out)))
        out += _bgzf_block(raw[i:i + block])
    out += _bgzf_block(b"")
    with open(path, "wb") as fh:
        fh.write(bytes(out))

    def voff(u):
        k = max(j for j, (ui, _c) in enumerate(starts) if ui <= u)
        return (starts[k][1] << 16) | (u - starts[k][0])

    refs = {}
    order = []
    u = 0
    for line in text.split("\n"):
        nxt = u + len(line.encode()) + 1
        if line and not line.startswith("#"):
            f = line.split("\t")
            chrom, beg = f[0], int(f[1]) - 1
            end = beg + len(f[3])
            if chrom not in refs:
                refs[chrom] = ({}, {})
                order.append(chrom)
            bins, lin = refs[chrom]
            b = _reg2bin(beg, end)
            bins.setdefault(b, []).append((voff(u), voff(nxt) if nxt < len(raw) else
                                           (len(out) - 28) << 16))
            for w in range(beg >> 14, ((end - 1) >> 14) + 1):
                lin.setdefault(w, voff(u))
        u = nxt
    names = b"".join(n.encode() + b"\x00" for n in order)
    idx = bytearray(b"TBI\x01")
    idx += struct.pack("<8i", len(order), 2, 1, 2, 0, ord("#"), 0, len(names)) + names
    for n in order:
        bins, lin = refs[n]
        idx += struct.pack("<i", len(bins))
        for b, chunks in sorted(bins.items()):
            idx += struct.pack("<Ii", b, 1) + struct.pack("<QQ", chunks[0][0], chunks[-1][1])
        top = max(lin) + 1 if lin else 0
        idx += struct.pack("<i", top)
        last = 0
        for w in range(top):
            last = lin.get(w, last)
            idx += struct.pack("<Q", last)
    with gzip.open(path + ".tbi", "wb") as fh:
        fh.write(bytes(idx))


HEAD = (["Spacer+PAM", "Chromosome"]
        + sum(([
            "Start_coordinate_(%s)" % b, "Strand_(%s)" % b, "Aligned_spacer+PAM_(%s)" % b,
            "Aligned_protospacer+PAM_REF_(%s)" % b, "Aligned_protospacer+PAM_ALT_(%s)" % b,
            "PAM_(%s)" % b, "Mismatches_(%s)" % b, "Bulges_(%s)" % b,
            "Mismatches+bulges_(%s)" % b, "Bulge_type_(%s)" % b,
            "REF/ALT_origin_(%s)" % b, "Variant_info_genome_(%s)" % b,
            "Variant_samples_(%s)" % b] for b in C.BLOCKS), []))
OFF = {"start": 0, "strand": 1, "crrna": 2, "ref": 3, "alt": 4, "pam": 5, "mm": 6,
       "bul": 7, "total": 8, "btype": 9, "origin": 10, "var": 11, "samples": 12}


def row_fields(chrom, start, strand, crrna, ref, alt, pam, mm, bul, btype, var, samples,
               spacer=GUIDE + "NNN"):
    out = [spacer, chrom]
    for _b in C.BLOCKS:
        out += [str(start), strand, crrna, ref, alt, pam, str(mm), str(bul),
                str(mm + bul), btype, "alt", var, samples]
    return out


def set_field(row, **kw):
    row = list(row)
    per = (len(row) - 2) // len(C.BLOCKS)
    for key, value in kw.items():
        for b in range(len(C.BLOCKS)):
            row[2 + b * per + OFF[key]] = str(value)
    return row


def write_results(folder, rows, name="G+NRG_hg38+v_4+2_integrated_results.tsv"):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, name)
    with open(path, "w") as fh:
        fh.write("\t".join(HEAD) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return path


def tok(chrom, rec):
    return "%s_%d_%s_%s" % (chrom, rec[0], rec[1], rec[2])


def fill_markers(res):
    """The markers complete-search writes for every run since 2.8.0, where the test wrote
    none of its own: .search_mode (population-level where the test made the run write no
    alternative alignments, or gave it --population-level; else per-sample) and .alt_mode
    (the mode's default). Returns the files it wrote."""
    sm, am = os.path.join(res, ".search_mode"), os.path.join(res, ".alt_mode")
    cl = os.path.join(res, ".command_line.txt")
    made = []
    alt = None
    if os.path.exists(am):
        with open(am) as fh:
            alt = fh.read().strip() == "alt-alignments"
    if not os.path.exists(sm):
        pop = alt is False
        if os.path.exists(cl):
            with open(cl) as fh:
                pop = pop or "--population-level" in fh.read().split()
        with open(sm, "w") as fh:
            fh.write("population-level" if pop else "per-sample")
        made.append(sm)
    if alt is None:
        with open(sm) as fh:
            ps = fh.read().strip() == "per-sample"
        with open(am, "w") as fh:
            fh.write("alt-alignments" if ps else "no-alt-alignments")
        made.append(am)
    return made


def run_markers(res, mode, alt=None):
    """A run's .search_mode and .alt_mode, as complete-search writes them: the mode it ran
    and whether it wrote alternative alignments (by default per sample, not
    population-level)."""
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, ".search_mode"), "w") as fh:
        fh.write(mode)
    if alt is None:
        alt = mode == "per-sample"
    with open(os.path.join(res, ".alt_mode"), "w") as fh:
        fh.write("alt-alignments" if alt else "no-alt-alignments")


class Case(object):
    """One built world: a reference, a VCF, a results folder and the argv to check it."""

    def __init__(self, tmp, seq, records, rows, chrom="chr21", pam="N" * 20 + "NRG 3"):
        self.tmp = tmp
        self.chrom = chrom
        os.makedirs(os.path.join(tmp, "Genomes"), exist_ok=True)
        write_fasta(os.path.join(tmp, "Genomes", chrom + ".fa"), seq, chrom)
        self.vcf = os.path.join(tmp, "v.vcf.gz")
        write_vcf(self.vcf, chrom, records)
        self.pam = os.path.join(tmp, "pam.txt")
        with open(self.pam, "w") as fh:
            fh.write(pam + "\n")
        self.res = os.path.join(tmp, "Results")
        write_results(self.res, rows)
        # a new run: no record of an earlier one's mode in its folder (fill_markers writes
        # this one's while run() checks it, where the test wrote none)
        for name in (".search_mode", ".alt_mode", ".command_line.txt"):
            if os.path.exists(os.path.join(self.res, name)):
                os.remove(os.path.join(self.res, name))

    def run(self, mode="soundness", extra=(), budget=("4", "1", "1", "6"), guide=GUIDE):
        made = fill_markers(self.res)
        try:
            return self._run(mode, extra, budget, guide)
        finally:
            for path in made:
                os.remove(path)

    def _run(self, mode, extra, budget, guide):
        out = os.path.join(self.tmp, "fail.tsv")
        if os.path.exists(out):
            os.remove(out)
        argv = ["--results", self.res, "--reference", os.path.join(self.tmp, "Genomes"),
                "--vcf", self.vcf, "--pam", self.pam, "--mm", budget[0], "--bDNA",
                budget[1], "--bRNA", budget[2], "--max-total-edits", budget[3],
                "--mode", mode, "--chrom", self.chrom, "--guide", guide,
                "--out", out] + list(extra)
        found = []
        rc = C.main(argv)
        if os.path.exists(out):
            with open(out) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    f = dict(zip(head, line.rstrip("\n").split("\t")))
                    found.append((f["gating"], f["class"], f["detail"], f["samples"]))
        return rc, found


def results_table(case):
    return os.path.join(case.res, "G+NRG_hg38+v_4+2_integrated_results.tsv")


def gating(found):
    return [f for f in found if f[0] == "yes"]


def classes(found):
    return {f[1] for f in found}


# --------------------------------------------------------------------------- #
# the worlds
# --------------------------------------------------------------------------- #
def world_deletion(tmp):
    """A deletion restores a perfect target: one extra C after GUIDE[9] is deleted.
    S01 and S02 carry it, S03 does not."""
    seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
    pos1 = len(PRE) + 10
    rec = (pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
    start = len(PRE)
    var = tok("chr21", rec)
    ref_text = seq[start:start + 23]
    rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                       "TGG", 0, 0, "X", var, "S01,S02")]
    return Case(tmp, seq, [rec], rows), dict(
        seq=seq, pos1=pos1, start=start, var=var, ref_text=ref_text, rows=rows, rec=rec)


def world_insertion(tmp):
    """An insertion holds the target's first bases: the row starts at the anchor."""
    tail = GUIDE[6:] + "TGG"
    seq = PRE + tail + POST
    anchor0 = len(PRE) - 1
    added = "TTCAT" + GUIDE[:6]
    rec = (anchor0 + 1, seq[anchor0], seq[anchor0] + added, {"S04": "1|0"})
    var = tok("chr21", rec)
    ref_text = seq[anchor0:anchor0 + 23]
    rows = [row_fields("chr21", anchor0, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                       "TGG", 0, 0, "X", var, "S04")]
    return Case(tmp, seq, [rec], rows), dict(
        seq=seq, anchor0=anchor0, var=var, ref_text=ref_text, rows=rows, added=added)


def world_indel_snp(tmp):
    """A deletion and a spacer SNP the guide needs, in cis on S05 and in trans on S06."""
    seq = PRE + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST
    dpos = len(PRE) + 10
    spos = len(PRE) + 6
    recs = [(spos, "C", GUIDE[5], {"S05": "1|0", "S06": "0|1"}),
            (dpos, GUIDE[9] + "C", GUIDE[9], {"S05": "1|0", "S06": "1|0"})]
    start = len(PRE)
    var = "%s,%s" % (tok("chr21", recs[1]), tok("chr21", recs[0]))
    ref_text = seq[start:start + 23]
    rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                       "TGG", 0, 0, "X", var, "S05")]
    return Case(tmp, seq, recs, rows), dict(seq=seq, start=start, var=var, rows=rows,
                                            recs=recs)


def world_mnv(tmp, carriers=None):
    """An MNV AC>GU... restores the guide; S01 carries it phased (CRISPRme reads an MNV's
    carriers with the phase unknown: PUTATIVE)."""
    seq = PRE + GUIDE[:8] + "AC" + GUIDE[10:] + "TGG" + POST
    pos1 = len(PRE) + 9
    rec = (pos1, "AC", GUIDE[8:10], carriers or {"S01": "1|0"})
    start = len(PRE)
    var = tok("chr21", rec)
    rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                       GUIDE + "TGG", "TGG", 0, 0, "X", var, "S01")]
    return Case(tmp, seq, [rec], rows), dict(seq=seq, start=start, var=var, rows=rows)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="chkidl_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# right rows raise nothing
# --------------------------------------------------------------------------- #
class TrueRows(Base):
    def assertClean(self, case, mode="soundness", **kw):
        rc, found = case.run(mode=mode, **kw)
        self.assertEqual((rc, gating(found)), (0, []), found)
        return found

    def test_deletion_row_passes(self):
        case, _ = world_deletion(self.tmp)
        self.assertEqual(self.assertClean(case), [])

    def test_insertion_row_at_the_anchor_passes(self):
        case, _ = world_insertion(self.tmp)
        self.assertEqual(self.assertClean(case), [])

    def test_indel_plus_cis_snp_row_passes(self):
        case, _ = world_indel_snp(self.tmp)
        self.assertEqual(self.assertClean(case), [])

    def test_minus_strand_row_passes(self):
        fwd = C.rc(GUIDE + "TGG")
        seq = PRE + fwd[:9] + "A" + fwd[9:] + POST
        pos1 = len(PRE) + 9
        rec = (pos1, fwd[8] + "A", fwd[8], {"S07": "1|0"})
        start = len(PRE)
        rows = [row_fields("chr21", start, "-", GUIDE + "NNN", C.rc(seq[start:start + 23]),
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S07")]
        self.assertEqual(self.assertClean(Case(self.tmp, seq, [rec], rows)), [])

    def test_homozygous_carrier_passes(self):
        _case, info = world_deletion(self.tmp)
        rec = (info["pos1"], GUIDE[9] + "C", GUIDE[9], {"S01": "1|1", "S02": "0|1"})
        self.assertClean(Case(self.tmp, info["seq"], [rec], info["rows"]))

    def test_a_multi_alt_line_passes(self):
        """GC -> G,GA on one line, S01 1|2: the deletion on slot 0."""
        case, info = world_deletion(self.tmp)
        write_vcf(case.vcf, "chr21", [(info["pos1"], GUIDE[9] + "C",
                                       GUIDE[9] + "," + GUIDE[9] + "A",
                                       {"S01": "1|2", "S02": "0|1"})])
        self.assertClean(case)

    def test_minus_strand_insertion_at_the_anchor_passes(self):
        fwd = C.rc(GUIDE + "TGG")
        seq = PRE + fwd[5:] + POST
        anchor0 = len(PRE) - 1
        rec = (anchor0 + 1, seq[anchor0], seq[anchor0] + "TT" + fwd[:5], {"S03": "0|1"})
        rows = [row_fields("chr21", anchor0, "-", GUIDE + "NNN",
                           C.rc(seq[anchor0:anchor0 + 23]), GUIDE + "TGG", "TGG", 0, 0,
                           "X", tok("chr21", rec), "S03")]
        case = Case(self.tmp, seq, [rec], rows)
        self.assertClean(case)
        self.assertClean(case, mode="completeness", budget=("0", "0", "0", "0"))

    def test_a_pam_made_of_inserted_bases_passes(self):
        seq = PRE + GUIDE + "CAC" + POST
        anchor0 = len(PRE) + len(GUIDE) - 1
        rec = (anchor0 + 1, GUIDE[-1], GUIDE[-1] + "TGG", {"S01": "1|0"})
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01")]
        case = Case(self.tmp, seq, [rec], rows)
        self.assertClean(case)
        self.assertClean(case, mode="completeness", budget=("0", "0", "0", "0"))

    def test_a_deletion_joining_the_pam_passes(self):
        seq = PRE + GUIDE + "TAAAGG" + POST
        pos1 = len(PRE) + len(GUIDE) + 1
        rec = (pos1, "TAAA", "T", {"S01": "0|1"})
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01")]
        case = Case(self.tmp, seq, [rec], rows)
        self.assertClean(case)
        self.assertClean(case, mode="completeness", budget=("0", "0", "0", "0"))

    def test_a_cas12a_insertion_at_the_anchor_passes(self):
        """5' PAM TTTV: the insertion supplies the PAM and the first spacer bases."""
        seq = PRE + GUIDE[3:] + POST
        anchor0 = len(PRE) - 1
        rec = (anchor0 + 1, seq[anchor0], seq[anchor0] + "GG" + "TTTA" + GUIDE[:3],
               {"S04": "1|0"})
        rows = [row_fields("chr21", anchor0, "+", "NNNN" + GUIDE, seq[anchor0:anchor0 + 24],
                           "TTTA" + GUIDE, "TTTA", 0, 0, "X", tok("chr21", rec), "S04",
                           spacer="NNNN" + GUIDE)]
        case = Case(self.tmp, seq, [rec], rows, pam="TTTV" + "N" * 20 + " -4")
        self.assertClean(case)
        self.assertClean(case, mode="completeness", budget=("0", "0", "0", "0"))

    def test_a_dna_bulge_next_to_a_5prime_pam_passes(self):
        """CRISPRme writes NNNN-CCTTG... for Cas12a: a DNA bulge right after the PAM."""
        seq = PRE + "TTTA" + "G" + GUIDE[:10] + "C" + GUIDE[10:] + POST
        pos1 = len(PRE) + 5 + 10
        rec = (pos1, GUIDE[9] + "C", GUIDE[9], {"S02": "1|0"})
        start = len(PRE)
        crrna = "NNNN-" + GUIDE
        alt = "TTTAG" + GUIDE
        rows = [row_fields("chr21", start, "+", crrna, seq[start:start + 25], alt, "TTTA",
                           0, 1, "DNA", tok("chr21", rec), "S02", spacer="NNNN" + GUIDE)]
        case = Case(self.tmp, seq, [rec], rows, pam="TTTV" + "N" * 20 + " -4")
        self.assertClean(case)

    def test_an_mnv_row_passes(self):
        case, _ = world_mnv(self.tmp)
        self.assertClean(case)

    def test_an_unphased_row_passes(self):
        """Unphased carriers of the deletion and of the SNP the guide needs: PUTATIVE,
        listed for the strongest target their own alleles allow."""
        case, info = world_indel_snp(self.tmp)
        recs = [(r[0], r[1], r[2], {"S05": "0/1"}) for r in info["recs"]]
        self.assertClean(Case(self.tmp, info["seq"], recs, info["rows"]))

    def test_a_missing_allele_away_from_the_target_passes(self):
        """F7: '1|.' at a SNP 30 bp outside the target does not make the carrier wrong."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST + LONG
        pos1 = len(PRE) + 10
        snp = len(PRE) + 24 + 30
        recs = [(pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"}),
                (snp, seq[snp - 1], other(seq[snp - 1]), {"S02": "1|."})]
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]),
                           "S01,S02")]
        self.assertClean(Case(self.tmp, seq, recs, rows))

    def test_an_unnormalized_deletion_passes_and_is_demanded(self):
        """N1: GCA>GA (C deleted, A shared) is the deletion GC>G."""
        case, info = world_deletion(self.tmp)
        seq, p = info["seq"], info["pos1"]
        nxt = seq[p + 1]
        rec = (p, GUIDE[9] + "C" + nxt, GUIDE[9] + nxt, {"S01": "1|0", "S02": "0|1"})
        rows = [set_field(info["rows"][0], var=tok("chr21", rec))]
        case = Case(self.tmp, seq, [rec], rows)
        self.assertClean(case)
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertIn("missing-row", classes(found), found)

    def test_a_phasing_artefact_at_the_records_position_keeps_the_carrier(self):
        """F8 / oracle site P: on S45's slot another insertion and a SNP share the
        insertion's POS (alternative alleles): it stays listed with the window's
        sequence."""
        seq = PRE + GUIDE[:10] + GUIDE[13:] + "TGG" + POST
        pos1 = len(PRE) + 10
        recs = [(pos1, GUIDE[9], GUIDE[9] + GUIDE[10:13], {"S01": "1|0", "S02": "0|1"}),
                (pos1, GUIDE[9], GUIDE[9] + "AAA", {"S01": "1|0"}),
                (pos1, GUIDE[9], other(GUIDE[9]), {"S01": "1|0"})]
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]),
                           "S01,S02")]
        self.assertClean(Case(self.tmp, seq, recs, rows))

    def test_a_tandem_second_deletion_keeps_the_carrier_unnamed(self):
        """Oracle site W: a 1-bp deletion in the G run past the PAM leaves the target's
        bases (the run puts a G back): S01 keeps the row, which names only the first."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGGG" + "A" + POST
        d1 = len(PRE) + 10
        d2 = len(PRE) + 22
        recs = [(d1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"}),
                (d2, "TG", "T", {"S01": "1|0"})]
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]),
                           "S01,S02")]
        self.assertClean(Case(self.tmp, seq, recs, rows))


# --------------------------------------------------------------------------- #
# each corruption of a right row is caught
# --------------------------------------------------------------------------- #
class Mutations(Base):
    def _mutate(self, build, **kw):
        case, info = build(self.tmp)
        write_results(case.res, [set_field(info["rows"][0], **kw)])
        return case.run()

    def assertCaught(self, rc_found, *want):
        rc, found = rc_found
        self.assertEqual(rc, 1, found)
        self.assertTrue(set(want) & classes(gating(found)), found)

    def test_start_shifted_by_one_is_caught(self):
        _c, info = world_deletion(self.tmp)
        self.assertCaught(self._mutate(world_deletion, start=info["start"] + 1),
                          "start-not-on-variant-sequence")

    def test_start_shifted_by_the_insertion_offset_is_caught(self):
        """The 2.8.1 defect: the row at the first base's place on the haplotype."""
        _c, info = world_insertion(self.tmp)
        self.assertCaught(self._mutate(world_insertion, start=info["anchor0"] + 6),
                          "start-not-on-variant-sequence")

    def test_start_moved_to_the_base_after_the_insertion_is_caught(self):
        _c, info = world_insertion(self.tmp)
        self.assertCaught(self._mutate(world_insertion, start=info["anchor0"] + 1),
                          "start-not-on-variant-sequence")

    def test_minus_strand_insertion_shift_is_caught(self):
        fwd = C.rc(GUIDE + "TGG")
        seq = PRE + fwd[5:] + POST
        anchor0 = len(PRE) - 1
        rec = (anchor0 + 1, seq[anchor0], seq[anchor0] + "TT" + fwd[:5], {"S03": "0|1"})
        rows = [row_fields("chr21", anchor0 + 3, "-", GUIDE + "NNN",
                           C.rc(seq[anchor0 + 3:anchor0 + 26]), GUIDE + "TGG", "TGG", 0, 0,
                           "X", tok("chr21", rec), "S03")]
        self.assertCaught(Case(self.tmp, seq, [rec], rows).run(),
                          "start-not-on-variant-sequence")

    def test_one_target_base_changed_is_caught(self):
        alt = GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:] + "TGG"
        self.assertCaught(self._mutate(world_deletion, alt=alt, mm=1, total=1),
                          "target-not-on-variant-sequence")

    def test_counts_not_recounting_is_caught(self):
        self.assertCaught(self._mutate(world_deletion, mm=2, total=2), "counts-mismatch")

    def test_a_pam_the_search_never_matched_is_caught(self):
        self.assertCaught(self._mutate(world_deletion, alt=GUIDE + "TGA", pam="TGA"),
                          "pam-mismatch")

    def test_ref_text_read_at_the_wrong_place_is_caught(self):
        _c, info = world_deletion(self.tmp)
        wrong = info["seq"][info["start"] + 3:info["start"] + 26]
        self.assertCaught(self._mutate(world_deletion, ref=wrong), "ref-text-mismatch")

    def test_bulge_type_x_on_a_bulged_row_is_caught(self):
        """N3: a DNA bulge reported with Bulge_type X."""
        seq = PRE + GUIDE + "TGG" + POST
        pos1 = len(PRE) + 10
        rec = (pos1, GUIDE[9], GUIDE[9] + "A", {"S01": "1|0"})
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE[:10] + "-" + GUIDE[10:] + "NNN",
                           seq[start:start + 24], GUIDE[:10] + "A" + GUIDE[10:] + "TGG",
                           "TGG", 0, 1, "X", tok("chr21", rec), "S01")]
        self.assertCaught(Case(self.tmp, seq, [rec], rows).run(), "bulge-type-mismatch")

    def test_a_gap_in_the_pam_is_caught(self):
        self.assertCaught(self._mutate(world_deletion, crrna=GUIDE + "N-NN",
                                       alt=GUIDE + "TAGG"), "alignment-invalid")

    def test_an_aligned_spacer_that_is_not_the_guide_is_caught(self):
        """F5: the spacer rewritten to match a mismatching target, 0 mismatches."""
        bad = GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:]
        seq = PRE + bad[:10] + "C" + bad[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        rec = (pos1, bad[9] + "C", bad[9], {"S01": "1|0"})
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", bad + "NNN", seq[start:start + 23],
                           bad + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01")]
        self.assertCaught(Case(self.tmp, seq, [rec], rows).run(), "spacer-not-guide")

    def test_a_row_of_another_guide_is_caught(self):
        case, info = world_deletion(self.tmp)
        rc, found = case.run(guide="TTTTAAAATTTTAAAATTTT")
        self.assertEqual(rc, 1)
        self.assertIn("guide-unknown", classes(found), found)

    # --- carriers -------------------------------------------------------- #
    def test_a_carrier_swapped_for_a_non_carrier_is_caught(self):
        self.assertCaught(self._mutate(world_deletion, samples="S03"),
                          "carrier-without-indel")
        rc, found = self._mutate(world_deletion, samples="S03")
        self.assertIn("carrier-missing", classes(found), found)

    def test_a_carrier_added_beside_the_real_ones_is_caught(self):
        self.assertCaught(self._mutate(world_deletion, samples="S01,S02,S03"),
                          "carrier-without-indel")

    def test_a_missing_carrier_is_caught(self):
        rc, found = self._mutate(world_deletion, samples="S01")
        self.assertCaught((rc, found), "carrier-missing")
        self.assertTrue(any("S02" in f[3] for f in found if f[1] == "carrier-missing"))

    def test_a_missing_carrier_is_not_hidden_by_the_same_target_elsewhere(self):
        """F6: two copies of one locus 200 bp apart, each with the deletion. S02 carries
        both; row A lists S01 only, row B (the other copy) S02."""
        unit = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + "ACGTTGCAATTGG"
        spacer = LONG * 3
        seq = "AAAACCCCGGGG" + unit + spacer + unit + "CCAATTGCAGGA"
        p1 = 12 + len(PRE) + 10
        p2 = 12 + len(unit) + len(spacer) + len(PRE) + 10
        d = GUIDE[9] + "C"
        recs = [(p1, d, GUIDE[9], {"S01": "1|0", "S02": "1|0"}),
                (p2, d, GUIDE[9], {"S02": "1|0"})]
        s1 = 12 + len(PRE)
        s2 = 12 + len(unit) + len(spacer) + len(PRE)
        rows = [row_fields("chr21", s1, "+", GUIDE + "NNN", seq[s1:s1 + 23], GUIDE + "TGG",
                           "TGG", 0, 0, "X", tok("chr21", recs[0]), "S01"),
                row_fields("chr21", s2, "+", GUIDE + "NNN", seq[s2:s2 + 23], GUIDE + "TGG",
                           "TGG", 0, 0, "X", tok("chr21", recs[1]), "S02")]
        rc, found = Case(self.tmp, seq, recs, rows).run()
        self.assertTrue(any(f[1] == "carrier-missing" and "S02" in f[3] for f in found),
                        found)

    # --- variants -------------------------------------------------------- #
    def test_a_named_variant_dropped_is_caught(self):
        _c, info = world_indel_snp(self.tmp)
        only_indel = info["var"].split(",")[0]
        self.assertCaught(self._mutate(world_indel_snp, var=only_indel),
                          "target-not-on-variant-sequence")

    def test_a_variant_that_is_not_a_record_is_caught(self):
        _c, info = world_deletion(self.tmp)
        self.assertCaught(self._mutate(world_deletion, var=info["var"] + ",chr21_9_A_AAAA"),
                          "variant-not-in-vcf")

    def test_a_cross_haplotype_snp_added_to_the_row_is_caught(self):
        """S05 carries the SNP in cis, S06 in trans: S06 on the cis row is the 2.8.0
        cross-person pairing in miniature."""
        self.assertCaught(self._mutate(world_indel_snp, samples="S05,S06"),
                          "carrier-without-variants", "carrier-without-target")

    def test_an_unnamed_variant_of_a_listed_carrier_is_caught(self):
        """A listed carrier whose slot holds a SNP inside the target the row omits."""
        case, info = world_deletion(self.tmp)
        snp = info["start"] + 3 + 1
        recs = [info["rec"], (snp, info["seq"][snp - 1], other(info["seq"][snp - 1]),
                              {"S01": "1|0"})]
        rc, found = Case(self.tmp, info["seq"], recs, info["rows"]).run()
        self.assertEqual(rc, 1)
        self.assertIn("variant-missing", classes(found), found)

    def test_the_phantom_manuel_reported_is_caught(self):
        """An insertion's only carrier listed for a sequence pairing it with SNP alleles
        two other samples carry, nobody carrying any two of the three: the 2.8.0 row at
        chr22:22322768 in miniature."""
        seq = PRE + GUIDE[6:] + "TGG" + POST
        anchor0 = len(PRE) - 1
        wrong = list(GUIDE[:6])
        wrong[2] = other(wrong[2])
        added = "TTCAT" + "".join(wrong)
        recs = [(anchor0 + 1, seq[anchor0], seq[anchor0] + added, {"S04": "1|0"}),
                (len(PRE) + 2, seq[len(PRE) + 1], other(seq[len(PRE) + 1]), {"S05": "1|0"}),
                (len(PRE) + 5, seq[len(PRE) + 4], other(seq[len(PRE) + 4]), {"S06": "0|1"})]
        rows = [row_fields("chr21", anchor0, "+", GUIDE + "NNN", seq[anchor0:anchor0 + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]), "S04")]
        self.assertCaught(Case(self.tmp, seq, recs, rows).run(),
                          "target-not-on-variant-sequence")

    # --- MNV and unphased rows (F4) -------------------------------------- #
    def test_an_mnv_row_with_a_phantom_base_is_caught(self):
        bad = list(GUIDE + "TGG")
        bad[2] = other(bad[2])
        self.assertCaught(self._mutate(world_mnv, alt="".join(bad), mm=1, total=1),
                          "target-not-on-variant-sequence")

    def test_an_mnv_row_shifted_by_one_is_caught(self):
        _c, info = world_mnv(self.tmp)
        s = info["start"] + 1
        self.assertCaught(self._mutate(world_mnv, start=s, ref=info["seq"][s:s + 23]),
                          "start-not-on-variant-sequence")

    def test_an_unphased_row_with_a_phantom_base_is_caught(self):
        case, info = world_deletion(self.tmp)
        rec = (info["pos1"], GUIDE[9] + "C", GUIDE[9], {"S01": "0/1", "S02": "0/1"})
        bad = list(GUIDE + "TGG")
        bad[3] = other(bad[3])
        rows = [set_field(info["rows"][0], alt="".join(bad), mm=1, total=1)]
        self.assertCaught(Case(self.tmp, info["seq"], [rec], rows).run(),
                          "target-not-on-variant-sequence")

    def test_an_unphased_carrier_listed_for_an_allele_it_lacks_is_caught(self):
        """Oracle site H: homozygous for a SNP that spoils the target, so it has no
        reference allele there and is never listed for the perfect target."""
        case, info = world_deletion(self.tmp)
        snp = info["start"] + 4
        recs = [(info["pos1"], GUIDE[9] + "C", GUIDE[9], {"S01": "0/1", "S02": "0/1"}),
                (snp, info["seq"][snp - 1], other(info["seq"][snp - 1]), {"S01": "1/1"})]
        rc, found = Case(self.tmp, info["seq"], recs, info["rows"]).run()
        self.assertEqual(rc, 1, found)
        # S01 is homozygous for the SNP: its record's haplotype is known (the record and
        # the SNP), and it holds the SNP the row does not name
        self.assertTrue(classes(gating(found)) & {"carrier-without-allele",
                                                  "variant-missing"}, found)

    def test_a_complex_allele_carrier_listed_is_caught(self):
        """F8: TG>T and G>GT overlap on one slot: CRISPRme takes that slot off."""
        case, info = world_deletion(self.tmp)
        p = info["pos1"]
        recs = [info["rec"], (p + 1, "C", "CT", {"S01": "1|0"})]
        rc, found = Case(self.tmp, info["seq"], recs, info["rows"]).run()
        self.assertEqual(rc, 1, found)
        self.assertIn("carrier-off-by-rule", classes(found), found)


# --------------------------------------------------------------------------- #
# documented limits are counted under their own class, never as failures
# --------------------------------------------------------------------------- #
class Limits(Base):
    def test_last_alt_wins_listed_carrier_is_a_limit(self):
        """F3: S01 carries the earlier of two SNP records at one position inside the
        target; CRISPRme reads the later one and lists S01 with the reference base."""
        case, info = world_deletion(self.tmp)
        snp = info["start"] + 4
        b = info["seq"][snp - 1]
        recs = [info["rec"], (snp, b, other(b), {"S01": "1|0"}),
                (snp, b, "G" if b != "G" and other(b) != "G" else "T", {"S03": "1|0"})]
        rc, found = Case(self.tmp, info["seq"], recs, info["rows"]).run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-last-alt-wins", classes(found), found)

    def test_last_alt_wins_in_completeness_is_a_limit(self):
        """A target that needs the earlier record's alt is not demanded: a limit."""
        seq = PRE + GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        snp = len(PRE) + 5
        ref_b = seq[snp - 1]
        third = next(x for x in "ACGT" if x not in (ref_b, GUIDE[4]))
        recs = [(pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"}),
                (snp, ref_b, GUIDE[4], {"S01": "1|0"}),
                (snp, ref_b, third, {"S02": "1|0"})]
        case = Case(self.tmp, seq, recs, [])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-last-alt-wins", classes(found), found)

    def test_two_indels_in_one_target_are_a_limit(self):
        seq = PRE + GUIDE[:5] + "T" + GUIDE[5:10] + "C" + GUIDE[10:] + "TGG" + POST
        recs = [(len(PRE) + 5, GUIDE[4] + "T", GUIDE[4], {"S01": "1|0"}),
                (len(PRE) + 11, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"})]
        rc, found = Case(self.tmp, seq, recs, []).run(mode="completeness",
                                                      budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-two-indels", classes(found), found)

    def test_a_target_starting_inside_one_insertion_and_holding_a_deletion(self):
        """F13: insertion A adds the target's first bases, deletion B removes a base
        inside it, on one slot: the documented limit, not a missing row."""
        seq = PRE + GUIDE[4:12] + "C" + GUIDE[12:] + "TGG" + POST + LONG
        anchor_a = len(PRE) - 1
        recs = [(anchor_a + 1, seq[anchor_a], seq[anchor_a] + "CCA" + GUIDE[:4],
                 {"S01": "1|0"}),
                (len(PRE) + 8, GUIDE[11] + "C", GUIDE[11], {"S01": "1|0"})]
        case = Case(self.tmp, LONG + seq, [(r[0] + len(LONG), r[1], r[2], r[3])
                                           for r in recs], [])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-two-indels", classes(found), found)

    def test_a_complex_allele_carrier_left_off_is_not_missing(self):
        """F8 (HG01029's shape): the overlapping records take the slot off; a row
        without that carrier is right."""
        case, info = world_deletion(self.tmp)
        p = info["pos1"]
        recs = [info["rec"], (p + 1, "C", "CT", {"S02": "0|1"})]
        rows = [set_field(info["rows"][0], samples="S01")]
        case = Case(self.tmp, info["seq"], recs, rows)
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-overlapping-records", classes(found), found)
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "limit-overlapping-records" and "S02" in f[3]
                            for f in found), found)

    def test_an_mnv_carrier_known_phased_is_putative(self):
        """An MNV carrier whose phased VCF call puts a far-end SNP on the other slot is
        PUTATIVE to CRISPRme (MNVs have no phased store): listed on its strongest target
        it is a limit, not a failure."""
        seq = PRE + GUIDE[:8] + "AC" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 9
        far = len(PRE) + 1
        recs = [(pos1, "AC", GUIDE[8:10], {"S01": "1|0"}),
                (far, GUIDE[0], other(GUIDE[0]), {"S01": "0|1"})]
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]), "S01")]
        rc, found = Case(self.tmp, seq, recs, rows).run()
        self.assertEqual((rc, gating(found)), (0, []), found)


class SnpRowsReport(Base):
    def test_a_snp_row_carrier_with_an_indel_in_the_target_is_a_defect_report(self):
        """--snp-rows report: since 2.10.0 a per-sample SNP row does not list a carrier
        whose haplotype holding its SNPs also holds an indel of the phased store that
        changes the target (METHODS section 4, Track B): listing S02, whose slot of the
        SNP deletes a base of the target, is a CRISPRme defect, reported
        (snp-row-indel-carrier), never gating."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        snp = len(PRE) + 6
        dele = len(PRE) + 12
        recs = [(snp, other(GUIDE[5]), GUIDE[5], {"S01": "1|0", "S02": "1|0"}),
                (dele, seq[dele - 1:dele + 1], seq[dele - 1], {"S02": "1|0"})]
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]),
                           "S01,S02")]
        case = Case(self.tmp, seq, recs, rows)
        rc, found = case.run(extra=("--snp-rows", "report"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual(classes(found), {"snp-row-indel-carrier"}, found)
        self.assertTrue(all(f[3] == "S02" for f in found), found)
        # the row without S02, as CRISPRme 2.10.0 writes it: nothing to report, and S02,
        # left off, is not missing (its slot's indel changes the target)
        write_results(case.res, [set_field(rows[0], samples="S01")])
        rc, found = case.run(extra=("--snp-rows", "report"))
        self.assertEqual((rc, found), (0, []), found)

    def test_the_documented_exceptions_of_the_snp_row_drop(self):
        """The carriers 2.10.0 keeps on a SNP row though an indel or MNV of their slot
        changes the target: an MNV (in no store CRISPRme reads: limit-snp-row-mnv), an
        indel at the named SNP's own position (another allele of that base, a phasing
        artefact: limit-overlapping-records) and an unphased indel call
        (limit-snp-row-unphased). Reported, never gating; none is a defect."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        snp = len(PRE) + 6
        q = len(PRE) + 12          # 1-based: an MNV of two bases inside the target
        mnv = (q, seq[q - 1:q + 1], other(seq[q - 1]) + other(seq[q]), {"S02": "1|0"})
        art = (snp, seq[snp - 1], seq[snp - 1] + "T", {"S03": "1|0"})
        dele = (len(PRE) + 15, seq[len(PRE) + 14:len(PRE) + 16], seq[len(PRE) + 14],
                {"S04": "0/1"})
        recs = [(snp, other(GUIDE[5]), GUIDE[5],
                 {"S01": "1|0", "S02": "1|0", "S03": "1|0", "S04": "1|0"}), art, mnv, dele]
        recs.sort(key=lambda r: (r[0], len(r[1]), r[2]))
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X",
                           tok("chr21", (snp, other(GUIDE[5]), GUIDE[5])),
                           "S01,S02,S03,S04")]
        case = Case(self.tmp, seq, recs, rows)
        rc, found = case.run(extra=("--snp-rows", "report"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        who = {f[1]: f[3] for f in found}
        self.assertEqual(who, {"limit-snp-row-mnv": "S02",
                               "limit-overlapping-records": "S03",
                               "limit-snp-row-unphased": "S04"}, found)

    def test_a_repeat_unit_lost_at_the_targets_edge_keeps_the_snp_row_carrier(self):
        """Track B reads a carrier's slot through a repeat at the target's edge: S02's
        slot of the SNP also loses one CAG unit of a run the VCF writes at its left end,
        inside the target. Where the run stops at the target's end (the control), the
        deletion changes the target's bases: listing S02 is reported
        (snp-row-indel-carrier). Where the run goes on past the target's end, the bases
        read from the target's first base on are the target's: 2.10.0 keeps S02 (nothing
        to report), and S03, holding the SNP alone and left off, is missing."""
        unit = "CAG"
        start = len(PRE)
        snp = start + 6                  # 1-based POS of the SNP (the target's 6th base)
        recs_snp = (snp, other(GUIDE[5]), GUIDE[5])
        for runs_on, guide, body, after in (
                (False, GUIDE[:11] + unit * 3, GUIDE[:11] + unit * 3, "TGG" + "CAGTGG"),
                (True, GUIDE[:8] + unit * 4, GUIDE[:8] + unit * 4, unit * 3)):
            ref_t = body[:5] + other(body[5]) + body[6:]
            seq = PRE + ref_t + after + POST
            pam = (after[:3])
            d1 = start + len(body) - 3 * (4 if runs_on else 3)   # 0-based: the run's start
            dele = (d1, seq[d1 - 1] + unit, seq[d1 - 1], {"S02": "1|0"})
            self.assertEqual(seq[d1:d1 + 3], unit)
            recs = [recs_snp + ({"S01": "1|0", "S02": "1|0", "S03": "1|0"},), dele]
            row = row_fields("chr21", start, "+", guide + "NNN", seq[start:start + 23],
                             guide + pam, pam, 0, 0, "X", tok("chr21", recs_snp),
                             "S01,S02" if not runs_on else "S01,S02,S03",
                             spacer=guide + "NNN")
            tmp = os.path.join(self.tmp, "w%d" % runs_on)
            os.makedirs(tmp)
            case = Case(tmp, seq, recs, [row])
            rc, found = case.run(guide=guide, extra=("--snp-rows", "report"))
            want = {("snp-row-indel-carrier", "S02")} if not runs_on else set()
            got = {(f[1], f[3]) for f in found if f[1] != "carrier-missing"}
            self.assertEqual(got, want, (runs_on, found))
            write_results(case.res, [set_field(row, samples="S01,S02")])
            rc, found = case.run(guide=guide, extra=("--snp-rows", "report"))
            self.assertEqual({f[3] for f in found if f[1] == "carrier-missing"},
                             {"S03"}, (runs_on, found))

    def test_a_carrier_merged_from_another_alignment_of_the_site_is_listed_rightly(self):
        """CRISPRme's merge keys a SNP row by its Cluster Position and variants, with no
        alignment, so the carriers of a site's SNP rows of one variant list are joined:
        the row is a DNA-bulge alignment whose first base S02's slot changes (an SNV the
        shorter alignment of the same site does not reach); S02 holds that shorter
        alignment, which names the same SNP -- listed rightly, nothing to report. S03,
        listed and carrying nothing, is reported."""
        g = GUIDE
        x = other(g[5])
        seq = PRE + g[0] + "T" + g[1:5] + x + g[6:] + "TGG" + POST
        s0 = len(PRE)
        snp = (s0 + 7, x, g[5])                    # the guide's base at its column 5
        first = (s0 + 1, g[0], other(g[0]))         # only the bulged alignment reads it
        recs = [first + ({"S02": "1|0"},), snp + ({"S01": "1|0", "S02": "1|0"},)]
        t = g[0] + "T" + g[1:] + "TGG"
        bulged = row_fields("chr21", s0, "+", g[0] + "-" + g[1:] + "NNN", seq[s0:s0 + 24],
                            t, "TGG", 0, 1, "DNA", tok("chr21", snp), "S01,S02,S03")
        case = Case(self.tmp, seq, recs, [bulged])
        rc, found = case.run(extra=("--snp-rows", "report"))
        self.assertEqual({(f[1], f[3]) for f in found}, {("carrier-without-target", "S03")},
                         found)

    def test_an_unphased_carrier_listed_for_a_subset_of_its_variants_is_documented(self):
        """METHODS section 4: an unphased individual (a '/' call, homozygous too) is listed
        for every non-empty subset of the variants it carries in the target, also one that
        leaves out a variant it is homozygous for: limit-snp-row-unphased, not
        carrier-without-target."""
        g = GUIDE
        x, y = other(g[5]), other(g[9])
        seq = PRE + g[:5] + x + g[6:9] + y + g[10:] + "TGG" + POST
        a = (len(PRE) + 6, x, g[5])
        b = (len(PRE) + 10, y, g[9])
        recs = [a + ({"S01": "1/1"},), b + ({"S01": "1/1"},)]
        s0 = len(PRE)
        t = g[:9] + y + g[10:] + "TGG"     # a only: leaves out b, which S01 is 1/1 for
        row = row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], mark(g + "NNN", t),
                         "TGG", 1, 0, "X", tok("chr21", a), "S01")
        rc, found = Case(self.tmp, seq, recs, [row]).run(extra=("--snp-rows", "report"))
        self.assertEqual({(f[1], f[3]) for f in found}, {("limit-snp-row-unphased", "S01")},
                         found)

    def test_a_slot_the_documented_reading_cannot_read_leaves_the_snp_row(self):
        """METHODS section 4 (2.10.0): a slot whose store indel changes the target stays
        on the SNP row only where its own bases still read the target, from the target's
        first base on (from the first base the haplotype keeps, where an indel removes
        it) or up to its last. In a GA run that goes on past the target, S02's slot of
        the SNP loses the run's first unit, the target's first two bases: read from the
        next base on, its bases are the target's -- kept, nothing to report. An unphased
        insertion inside the target beside it leaves the slot unreadable (its slot cannot
        be told): CRISPRme drops S02, and listing it is reported; S03, holding the SNP
        alone and left off, is missing."""
        lo = len(PRE) + 1
        run = list("GA" * 14)
        run[10] = "T"                        # the reference breaks the run; the SNP mends it
        seq = PRE + "T" + "".join(run) + POST
        g = "GA" * 10                        # the target GAGA... and its PAM GAG (NRG)
        snp = (lo + 11, "T", "G")
        dele = (lo, "TGA", "T")              # POS = the base before the run: its 1st unit
        ins = (lo + 7, seq[lo + 6], seq[lo + 6] + "TT")
        t = "GA" * 11 + "G"
        s0 = lo
        self.assertEqual((seq[s0 + 10], seq[lo - 1:lo + 2]), ("T", "TGA"))
        for with_ins in (False, True):
            recs = [dele + ({"S02": "1|0"},),
                    snp + ({"S01": "1|0", "S02": "1|0", "S03": "0|1"},)]
            if with_ins:
                recs.append(ins + ({"S02": "0/1"},))
            recs.sort(key=lambda r: r[0])
            row = row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], t, t[-3:], 0,
                             0, "X", tok("chr21", snp), "S01,S02", spacer=g + "NNN")
            tmp = os.path.join(self.tmp, "w%d" % with_ins)
            os.makedirs(tmp)
            case = Case(tmp, seq, recs, [row])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
            keeps = C.snp_row_keeps(win, "S02", 0, s0, t, [win.by_key[snp][0]],
                                    win.lines_in(0, len(seq)))
            self.assertEqual(keeps, not with_ins)
            rc, found = case.run(guide=g, extra=("--snp-rows", "report"))
            got = {(f[1], f[3]) for f in found}
            want = {("carrier-missing", "S03")}
            if with_ins:
                want.add(("snp-row-indel-carrier", "S02"))
            self.assertEqual(got, want, (with_ins, found))

    def test_a_snp_row_naming_a_variant_outside_its_target_is_reported(self):
        """R9-F5 (the review's adv_r9a): the report runs the variant-outside-target check
        of an indel row: a SNP row naming an SNV 10 bases past its PAM, its carrier holding
        both, is reported (non-gating); the carrier, who holds the target with its own
        SNP, is not."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        snp = len(PRE) + 6
        far = len(PRE) + 23 + 10
        recs = [(snp, other(GUIDE[5]), GUIDE[5], {"S01": "1|0"}),
                (far + 1, seq[far], other(seq[far]), {"S01": "1|0"})]
        s = len(PRE)
        var = "%s,%s" % (tok("chr21", recs[0]), tok("chr21", recs[1]))
        rows = [row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], GUIDE + "TGG",
                           "TGG", 0, 0, "X", var, "S01")]
        rc, found = Case(self.tmp, seq, recs, rows).run(extra=("--snp-rows", "report"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual(classes(found), {"variant-outside-target"}, found)


class Edges(Base):
    def test_a_deletion_reaching_into_the_window_is_applied(self):
        """N3: a deletion that starts before the window's lower edge removes the bases it
        covers inside the window."""
        seq = "ACGT" * 10
        hap = C.make_hap(seq, [al(3, "GTACG", "G")], 5, 12)
        self.assertEqual(hap.clash, [])
        self.assertEqual(hap.seq, seq[7:12])
        self.assertEqual(hap.place[0], (7, 0))

    def test_a_record_within_a_target_of_the_contig_end_is_listed(self):
        seq = GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + "ACGT"
        rec = (10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"})
        rc, found = Case(self.tmp, seq, [rec], []).run(mode="completeness",
                                                       budget=("0", "0", "0", "0"))
        self.assertEqual(gating(found), [], found)
        self.assertIn("limit-contig-end", classes(found), found)

    def test_two_guides_are_each_checked(self):
        """Rows of two guides in one run: each row is read against its own guide."""
        case, info = world_deletion(self.tmp)
        g2 = "TTTTAAAATTTTAAAATTTT"
        rc, found = case.run(extra=("--guide", g2))
        self.assertEqual((rc, gating(found)), (0, []), found)


# --------------------------------------------------------------------------- #
# completeness
# --------------------------------------------------------------------------- #
class Completeness(Base):
    TIGHT = ("0", "0", "0", "0")

    def test_a_reported_target_satisfies_completeness(self):
        case, _ = world_deletion(self.tmp)
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_dropped_row_is_caught(self):
        case, _ = world_deletion(self.tmp)
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual(rc, 1)
        self.assertIn("missing-row", classes(found), found)

    def test_a_dropped_carrier_is_caught(self):
        case, info = world_deletion(self.tmp)
        write_results(case.res, [set_field(info["rows"][0], samples="S01")])
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual(rc, 1)
        self.assertTrue(any("S02" in f[3] for f in found if f[1] == "missing-row"), found)

    def test_a_flank_target_is_not_demanded(self):
        seq = PRE + "GT" + "CA" + GUIDE + "TGG" + POST
        rec = (len(PRE) + 1, "GT", "G", {"S01": "1|0"})
        rc, found = Case(self.tmp, seq, [rec], []).run(mode="completeness",
                                                       budget=self.TIGHT)
        self.assertEqual(gating(found), [], found)

    def test_a_row_with_other_variants_does_not_supersede(self):
        """N2: a reported row of the same sample and site naming another record does not
        stand for the missing one."""
        case, info = world_deletion(self.tmp)
        far = len(info["seq"]) - 5
        seq = info["seq"]
        rec2 = (far, seq[far - 1] + seq[far], seq[far - 1], {"S01": "0|1", "S02": "1|0"})
        recs = [info["rec"], rec2]
        rows = [set_field(info["rows"][0], var=tok("chr21", rec2))]
        case = Case(self.tmp, seq, recs, rows)
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertIn("missing-row", classes(found), found)

    def test_a_folded_alignment_is_satisfied_by_its_cluster_row(self):
        """An RNA-bulge alignment of a perfect site shares its Cluster Position and its
        carriers: the merge folds it into the perfect row."""
        case, info = world_deletion(self.tmp)
        rc, found = case.run(mode="completeness", budget=("0", "1", "1", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_the_fold_compares_crisprmes_own_carrier_sets(self):
        """The merge folds rows of one site, variant list and carrier set, and the
        carrier set is CRISPRme's: S02 holds an earlier record of a split site at the
        target's first base, which CRISPRme cannot see, so CRISPRme lists it for the
        perfect target as for the RNA-bulge one (which S02 holds), and the RNA-bulge
        alignment folds into the perfect row reported. (chr1:213641941 of the
        genome-wide README run is this case.)"""
        _case, info = world_deletion(self.tmp)
        seq, start = info["seq"], info["start"]
        first = start + 1  # 1-based POS of the target's first base
        recs = [info["rec"], (first, seq[start], "C", {"S02": "0|1"}),
                (first, seq[start], "G", {"S03": "1|0"})]
        case = Case(self.tmp, seq, recs, info["rows"])
        rc, found = case.run(mode="completeness", budget=("0", "1", "1", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_an_unphased_carrier_must_be_listed_at_the_site(self):
        case, info = world_indel_snp(self.tmp)
        recs = [(r[0], r[1], r[2], {"S05": "0/1"}) for r in info["recs"]]
        case = Case(self.tmp, info["seq"], recs, [])
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual(rc, 1, found)
        self.assertTrue(classes(gating(found)) & {"missing-putative-row", "missing-row"},
                        found)
        case = Case(self.tmp, info["seq"], recs, info["rows"])
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_population_level_run_is_checked_per_locus(self):
        """F9: a population-level run keeps one row per cluster; completeness asks only
        that the locus be reported."""
        case, info = world_deletion(self.tmp)
        run_markers(case.res, "population-level")
        write_results(case.res, [set_field(info["rows"][0], samples="S01")])
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", budget=self.TIGHT)
        self.assertEqual(rc, 1)
        self.assertIn("missing-locus", classes(found), found)


# --------------------------------------------------------------------------- #
# review round 2: the truth is the VCF, the documented view only classifies
# --------------------------------------------------------------------------- #
def world_strongest(tmp, gts):
    """A deletion and the spacer SNP that supplies the guide's base, with ``gts`` for
    S05 (the deletion's and the SNP's calls); the weak row lists S05 with the
    reference base at the SNP (1 mismatch), the strong row with the SNP (0)."""
    seq = PRE + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST + LONG
    dpos, spos = len(PRE) + 10, len(PRE) + 6
    snp = (spos, "C", GUIDE[5], {"S05": gts[1]})
    dele = (dpos, GUIDE[9] + "C", GUIDE[9], {"S05": gts[0]})
    start = len(PRE)
    ref_text = seq[start:start + 23]
    strong = row_fields("chr21", start, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG", "TGG",
                        0, 0, "X", "%s,%s" % (tok("chr21", dele), tok("chr21", snp)), "S05")
    weak = row_fields("chr21", start, "+", GUIDE + "NNN", ref_text,
                      GUIDE[:5] + "C" + GUIDE[6:] + "TGG", "TGG", 1, 0, "X",
                      tok("chr21", dele), "S05")
    return seq, [snp, dele], strong, weak


class ReviewRound2(Base):
    def test_a_putative_carrier_on_a_weaker_target_is_caught(self):
        """F2 (case A): unphased S05 carries the guide's base at the SNP: METHODS section
        4 lists it for that target only, never for the reference base there."""
        seq, recs, strong, weak = world_strongest(self.tmp, ("0/1", "0/1"))
        rc, found = Case(self.tmp, seq, recs, [strong, weak]).run()
        self.assertEqual(rc, 1, found)
        self.assertEqual(classes(gating(found)), {"carrier-not-strongest"}, found)
        rc, found = Case(self.tmp, seq, recs, [weak]).run()
        self.assertEqual(classes(gating(found)), {"carrier-not-strongest"}, found)
        rc, found = Case(self.tmp, seq, recs, [weak]).run(
            mode="completeness", budget=("1", "0", "0", "1"))
        self.assertIn("missing-putative-row", classes(gating(found)), found)
        rc, found = Case(self.tmp, seq, recs, [strong]).run(
            mode="both", budget=("1", "0", "0", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_phased_carrier_with_a_far_unknown_call_is_demanded_exactly(self):
        """F3 (case B): S05's slot 0 holds the deletion and the SNP in cis; one call 6 bp
        past the PAM is unphased or missing. The target touches no call of unknown phase,
        so its slot-0 target (0 mm) is demanded, and the right row passes."""
        seq = PRE + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST + LONG
        dpos, spos, bpos = len(PRE) + 10, len(PRE) + 6, len(PRE) + 15
        far = len(PRE) + 24 + 6
        for gt_far in ("0/1", "1|.", "0|1"):
            recs = [(spos, "C", GUIDE[5], {"S05": "1|0"}),
                    (dpos, GUIDE[9] + "C", GUIDE[9], {"S05": "1|0"}),
                    (bpos, seq[bpos - 1], other(seq[bpos - 1]), {"S05": "0|1"}),
                    (far, seq[far - 1], other(seq[far - 1]), {"S05": gt_far})]
            rc, found = Case(self.tmp, seq, recs, []).run(
                mode="completeness", budget=("0", "0", "0", "0"))
            self.assertIn("missing-row", classes(gating(found)), (gt_far, found))
            var = "%s,%s" % (tok("chr21", recs[1]), tok("chr21", recs[0]))
            start = len(PRE)
            rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                               GUIDE + "TGG", "TGG", 0, 0, "X", var, "S05")]
            for mm in ("0", "1", "2"):
                rc, found = Case(self.tmp, seq, recs, rows).run(
                    mode="both", budget=(mm, "0", "0", mm))
                self.assertEqual((rc, gating(found)), (0, []), (gt_far, mm, found))

    def test_region_reads_the_rows_of_its_records(self):
        """F4 (case C): a --region whose left edge falls between a right row's start and
        its record's POS."""
        case, info = world_deletion(self.tmp)
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"),
                             extra=["--region", "chr21:%d-%d" % (info["pos1"] - 2, 200)])
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"),
                             extra=["--region", "chr21:%d-%d" % (info["pos1"] - 2, 200)])
        self.assertIn("missing-row", classes(found), found)

    def test_a_same_position_snp_on_a_listed_slot_is_counted(self):
        """F1 (case D): a SNP at the insertion's own POS on its slot: the VCF describes
        the anchor twice, so there is no sequence to judge by; METHODS section 4 reads the
        slot without it (a phasing artefact), and the listed carrier, whose slot so read
        holds the target, is a note of its own (overlap-reading-agrees, R8-N1);
        completeness counts it (limit-overlapping-records), never passed over in
        silence."""
        tail = GUIDE[6:] + "TGG"
        seq = LONG + PRE + tail + POST
        anchor0 = len(LONG) + len(PRE) - 1
        a = seq[anchor0]
        g2 = a + GUIDE[1:]
        rec = (anchor0 + 1, a, a + GUIDE[1:6], {"S04": "1|0"})
        snp = (anchor0 + 1, a, other(a), {"S04": "1|0"})
        rows = [row_fields("chr21", anchor0, "+", g2 + "NNN", seq[anchor0:anchor0 + 23],
                           g2 + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S04",
                           spacer=g2 + "NNN")]
        case = Case(self.tmp, seq, [rec, snp], rows)
        rc, found = case.run(guide=g2)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual(classes(found), {"overlap-reading-agrees"}, found)
        self.assertTrue(any("S04" in f[3] for f in found), found)
        rc, found = case.run(guide=g2, mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-overlapping-records", classes(found), found)

    def test_an_indel_block_labelled_ref_is_checked(self):
        """N1 (case E): a shifted indel block whose origin column says ref."""
        case, info = world_deletion(self.tmp)
        write_results(case.res, [set_field(info["rows"][0], start=info["start"] + 1,
                                           origin="ref")])
        rc, found = case.run()
        self.assertEqual(rc, 1)
        self.assertTrue({"origin-not-alt", "start-not-on-variant-sequence"}
                        <= classes(gating(found)), found)

    def test_a_target_ending_inside_a_second_insertion_is_the_two_indel_limit(self):
        """Genome-wide (chr10:98432871, chr11:105091858 ...): the target's last bases are
        added by a second insertion anchored at its last reference coordinate. It needs two
        indel records, which no fake-indel window holds: the documented limit, not an
        unexplained miss."""
        seq = PRE + GUIDE[:10] + GUIDE[11:] + "T" + "AAAA" + POST
        k = len(PRE) + 10           # 1-based POS of GUIDE[9], the record's anchor
        t = len(PRE) + 10 + 9 + 1   # 1-based POS of the T the second insertion follows
        recs = [(k, GUIDE[9], GUIDE[9] + GUIDE[10], {"S01": "1|0"}),
                (t, "T", "TGG", {"S01": "1|0"})]
        case = Case(self.tmp, seq, recs, [])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-two-indels", classes(found), found)

    def test_a_second_indel_read_up_to_the_last_base_is_the_documented_limit(self):
        """A carrier that lost a CA repeat unit at the target's left edge holds the row's
        bases up to its last base, not from its start: its own haplotype holds the target
        two bases to the left. CRISPRme documents listing it (METHODS section 4); the
        truth disagrees, and the disagreement is a limit, not a pass and not a failure."""
        g = "CACATTGTCAGTACGATCGA"
        seq = LONG + "GTCA" + g[:10] + "C" + g[10:] + "TGG" + POST
        s0 = len(LONG) + 4
        d1 = (s0 + 10, g[9] + "C", g[9], {"S01": "1|0", "S02": "1|0"})
        d2 = (s0, seq[s0 - 1:s0 + 2], seq[s0 - 1], {"S01": "1|0"})
        rows = [row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], g + "TGG", "TGG",
                           0, 0, "X", tok("chr21", d1), "S01,S02", spacer=g + "NNN")]
        case = Case(self.tmp, seq, [d1, d2], rows)
        rc, found = case.run(guide=g)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "limit-two-indels" and f[3] == "S01" for f in found),
                        found)
        rc, found = case.run(guide=g, mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-two-indels", classes(found), found)

    def test_a_cluster_row_worse_than_the_target_is_caught(self):
        """F5: one row per cluster (no alternative alignments). The indel row's line is
        replaced by a worse alignment of the same cluster (the reference with the C as a
        DNA bulge, 1 edit): the cluster no longer holds a row as good as S01's target."""
        case, info = world_deletion(self.tmp)
        run_markers(case.res, "population-level")
        s0, seq = info["start"], info["seq"]
        worse = row_fields("chr21", s0, "+", GUIDE[:10] + "-" + GUIDE[10:] + "NNN",
                           seq[s0:s0 + 24], seq[s0:s0 + 24], "TGG", 0, 1, "DNA", "NA", "NA")
        worse = set_field(worse, origin="ref", var="NA", samples="NA")
        write_results(case.res, [worse])
        rc, found = case.run(mode="completeness", budget=("0", "1", "0", "1"))
        self.assertEqual(rc, 1, found)
        self.assertIn("missing-locus", classes(gating(found)), found)
        write_results(case.res, [worse, info["rows"][0]])
        rc, found = case.run(mode="completeness", budget=("0", "1", "0", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_the_cluster_check_takes_the_targets_best_alignment(self):
        """F5 on real data: chr22:25443200-25443340 of hg38 and the insertion
        chr22:25443277 T>TA, whose README-guide target at 25443269 (-) has many alignments
        (3 mm + 2 bulges at best, 4 + 2 the first the aligner meets). The review replaced
        its population-level row by a worse one of the cluster (6 edits): the cluster row
        must be as good as the target's best alignment."""
        g = "CTAACAGTTGCTTTTATCAC"
        seq = ("GCAACCTTTATGTTCTTCAGTAGGCGAATGTACCAACAAACTTTGTTACACTGATCAATGGAATATTATC"
               "CAGTGATAAAAAACCATTGAGATACATTTGGATTCGATTTATTAATCTCAGTTACCTACTATTCTGACAT")
        rec = (77, "T", "TA", {"S01": "1|0"})
        worse = row_fields("chr21", 69, "-", "CTAACAGTTG-CTTTTATCACNNN",
                           "-TATCAATGGTTTTTTATCACTGG", "NA", "TGG", 4, 2, "RNA,DNA", "NA",
                           "NA", spacer=g + "NNN")
        worse = set_field(worse, origin="ref", alt="NA")
        case = Case(self.tmp, seq, [rec], [worse])
        run_markers(case.res, "population-level")
        rc, found = case.run(mode="completeness", guide=g)
        hits = [f for f in gating(found) if f[1] == "missing-locus"
                and "Cluster Position 69 " in f[2]]
        self.assertTrue(hits and "(5 edits" in hits[0][2], found)

    def test_the_cluster_check_compares_cfd(self):
        rep = C.Reported([], C.Ctx(argparse_ns(), "NRG", 3, False, [GUIDE]),
                         loci={(GUIDE, "chr21"): [(GUIDE, [100, 101, 100], C.edits(2, 0),
                                                   0.300, 7)]})

        def d(total, cfd):
            return C.Demand("exact", [], C.edits(total, 0), cfd, 0)

        def alone(wide=False):
            return [102], False, set()
        self.assertTrue(rep.locus(GUIDE, "chr21", 102, d(3, 0.2), alone)[0])
        self.assertFalse(rep.locus(GUIDE, "chr21", 102, d(1, 0.2), alone)[0])  # fewer edits
        self.assertFalse(rep.locus(GUIDE, "chr21", 102, d(3, 0.5), alone)[0])  # higher CFD
        # F2 (review 3): a row past --merge bp is the cluster's only through alignments
        # that chain to it, never as the nearest row
        self.assertFalse(rep.locus(GUIDE, "chr21", 120, d(3, 0.2),
                                   lambda wide: ([120], False, set()))[0])
        self.assertEqual(rep.locus(GUIDE, "chr21", 120, d(3, 0.2), lambda wide: (
            [104, 107, 110, 113, 116, 118, 120], False, set()))[::3], (True, set()))
        # a chain cut at its limit is said to be so (the caller lists a coverage gap)
        self.assertEqual(rep.locus(GUIDE, "chr21", 120, d(3, 0.2),
                                   lambda wide: ([120], True, set()))[::2], (False, True))
        # R7-F2: a cluster that reaches its row only through rows CRISPRme may write or
        # not is said to, with why (the caller counts it, never passes it as a cluster row)
        self.assertEqual(rep.locus(GUIDE, "chr21", 120, d(3, 0.2), lambda wide: (
            ([104, 107, 110, 113, 116, 118, 120], False, {"capped"}) if wide
            else ([120], False, set())))[::3], (True, {"capped"}))
        self.assertFalse(rep.locus(GUIDE, "chr21", 120, d(3, 0.2), lambda wide: (
            ([120], False, set())))[0])
        # N3: another guide's row is not this guide's cluster
        self.assertFalse(rep.locus("TTTTAAAATTTTAAAATTTT", "chr21", 102, d(3, 0.2),
                                   alone)[0])

    def test_a_near_row_that_starts_inside_an_insertion_can_be_the_next_sites(self):
        """Since 2.10.0 the rows that start inside an insertion at one anchor are placed
        together, join one site and chain none (remove_contiguous_samples, rule of
        decision (d)): a row within --merge bp of a target that starts inside an insertion
        can be the next site's, and the target's own row further along the chain of every
        row. The gate's repeat-rich population-level run: the target of CFD 0.040 at
        Cluster Position 22,263,588 (+, VEGFA site 3) beside the new chr22:22,263,589 row of
        CFD 0.030 (an insertion's row), its site's row the chr22:22,263,577 one of CFD
        0.102, eleven bases before. It is counted under limit-insertion-site, not judged; a
        near row that names no insertion is the target's cluster's, as before, and a far row
        no better leaves it a miss. A near row that names an insertion and starts on a
        reference base is chained like any other: the target's site's, a miss (up to review
        r1 of lane M it was let through, M-R1-N1); and a near row whose block naming the
        insertion lies further than --merge bp from the target is no such row either."""
        ctx = C.Ctx(argparse_ns(), "NRG", 3, False, [GUIDE])

        def d(total, cfd):
            return C.Demand("exact", [], C.edits(total, 0), cfd, 0)

        def chain(wide=False):
            return [89, 92, 95, 98, 100], False, set()

        def inside(blk):
            return blk["inside"]
        for name, blocks, far_cfd, want in (
                ("inside", [(101, {"inside": True})], 0.102, (True, {"insertion-site"})),
                ("no insertion", [], 0.102, (False, set())),
                ("inside, the far row no better", [(101, {"inside": True})], 0.020,
                 (False, set())),
                ("on a reference base", [(101, {"inside": False})], 0.102, (False, set())),
                ("inside, its block far", [(110, {"inside": True})], 0.102,
                 (False, set()))):
            with self.subTest(name):
                far = (GUIDE, [89, 89, 89], C.edits(2, 0), far_cfd, 7, [])
                near = (GUIDE, [101, 101, 110], C.edits(3, 0), 0.030, 8, blocks)
                rep = C.Reported([], ctx, loci={(GUIDE, "chr21"): [far, near]})
                self.assertEqual(rep.locus(GUIDE, "chr21", 100, d(3, 0.040), chain,
                                           inside=inside)[::3], want)

    def test_an_integrated_line_lists_its_blocks_naming_an_insertion(self):
        # read_rows' per-line entry for the cluster check: its last field lists the
        # (Cluster Position, block) of each block of the line naming an insertion (an ALT
        # longer than its REF), which locus reads the first base of (row_first_k)
        r1 = row_fields("chr21", 69, "+", GUIDE + "NNN", GUIDE + "TGG", GUIDE + "TGG", "TGG",
                        0, 0, "X", "chr21_70_T_TA", "S01")
        r2 = row_fields("chr21", 80, "+", GUIDE + "NNN", GUIDE + "TGG", GUIDE + "TGG", "TGG",
                        0, 0, "X", "chr21_81_TA_T", "S01")
        p = write_results(os.path.join(self.tmp, "ins_res"), [r1, r2])
        loci = {}
        C.read_rows(p, "integrated", "chr21", None, None, 0, False, loci=loci)
        got = [[(cp, blk["start"], blk["var"]) for cp, blk in ln[5]]
               for lines in loci.values() for ln in lines]
        self.assertEqual([len(g) for g in got], [3, 0])
        self.assertEqual({(start, var) for _cp, start, var in got[0]}, {(69, "chr21_70_T_TA")})

    def test_a_short_table_line_is_reported(self):
        """N4: a results line with fewer fields than its header is not skipped silently."""
        case, info = world_deletion(self.tmp)
        write_results(case.res, [info["rows"][0], info["rows"][0][:5]])
        rc, found = case.run()
        self.assertEqual(rc, 1)
        self.assertIn("row-unreadable", classes(found), found)

    def test_a_window_with_ns_is_counted(self):
        """N4: a stretch of N next to a record: the windows across it are counted."""
        case, info = world_deletion(self.tmp)
        seq = info["seq"] + "NNNN" + "ACGTACGTACGTACGT" + "TGG" + LONG
        case = Case(self.tmp, seq, [info["rec"]], info["rows"])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-n-bases", classes(found), found)

    def test_the_checker_holds_no_crisprme_code(self):
        """F1: no import of CRISPRme's variant or indel modules, none of the helper names of
        PostProcess/indel_snp_cis.py or site_alleles.py, and no function of site_alleles
        named as the source of a rule (a regression guard against copying: R-B1's
        merge_call was site_alleles.merge_gt under another name)."""
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                               "scripts", "check_indel_rows.py")) as fh:
            src = fh.read()
        for mod in ("indel_snp_cis", "analisi_indels", "overlay_indel", "resultIntegrator",
                    "remove_contiguous_samples", "observed_haplotypes", "cfd_score",
                    "bulge_cfd_penalty", "site_alleles", "tier0_compile", "tier1_genotypes",
                    "simple_analysis_registry", "population_summary", "new_simple_analysis",
                    "haplotype_prioritizer"):
            self.assertNotRegex(src, r"(?m)^\s*(from|import)\s+%s" % mod)
        for name in ("kept_bases", "shared_tail", "changes_span", "target_ref_span",
                     "holds_indel", "local_haplotype", "target_pairs", "biallelic_site",
                     "resolve_carriers", "worst_case_bases", "_for_alt", "reference_start",
                     "window_place", "holds_target_bases", "snp_target_held", "merge_gt",
                     "merge_carriers", "read_records", "per_alt_entries", "per_alt",
                     "alt_tokens", "for_alt", "gt_tokens", "two_alts_on_a_slot"):
            self.assertNotRegex(src, r"\bdef %s\b" % name)
        self.assertNotRegex(src, r"site_alleles\.(?!py\b)\w")


# --------------------------------------------------------------------------- #
# review round 3: folds no worse, clusters by chaining, PUTATIVE on the whole CFD,
# long deletions, second indels that keep the bases, limits only where a target reads
# --------------------------------------------------------------------------- #
HEAD_CFD = HEAD + ["CFD_score_(%s)" % b for b in C.BLOCKS]


def write_results_cfd(folder, rows, cfds, name="G+NRG_hg38+v_4+2_integrated_results.tsv"):
    """Rows with a CFD column per block (``cfds``: one [CFD per block] per row)."""
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, name), "w") as fh:
        fh.write("\t".join(HEAD_CFD) + "\n")
        for r, c in zip(rows, cfds):
            fh.write("\t".join(list(r) + ["%.3f" % x for x in c]) + "\n")


def blocks_row(blocks):
    """One results line whose three blocks come from three row_fields() lines (BLOCKS
    order): a line whose blocks report different alignments."""
    per = (len(blocks[0]) - 2) // len(C.BLOCKS)
    out = list(blocks[0][:2])
    for k, b in enumerate(blocks):
        out += b[2 + k * per:2 + (k + 1) * per]
    return out


def mark(crrna, target):
    """A target as CRISPRme shows it: lower case where it mismatches the spacer."""
    return "".join(t if (g in "N-" or t == "-" or g == t) else t.lower()
                   for g, t in zip(crrna, target))


def world_inserted_targets(tmp, copy_mm=False):
    """An insertion whose added bases hold the guide's first bases: the target
    GUIDE + TGG starts at its 29th added base, and its RNA-bulge alignment before the
    same PAM at the 30th. With ``copy_mm`` the added bases also hold, from their first,
    the guide with one mismatch and its own PAM: another target starting inside the
    insertion, at the same Cluster Position (the anchor's), at another place on the
    haplotype."""
    first = (GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:] if copy_mm else
             "GATTACAGATTACAGATTACAGAT") + "TGG"
    added = first + "TTCAT" + GUIDE[:6]
    tail = GUIDE[6:] + "TGG"
    seq = LONG + PRE + tail + POST
    anchor0 = len(LONG) + len(PRE) - 1
    rec = (anchor0 + 1, seq[anchor0], seq[anchor0] + added, {"S01": "1|0", "S02": "0|1"})
    var = tok("chr21", rec)
    ref = seq[anchor0:anchor0 + 23]
    t1 = row_fields("chr21", anchor0, "+", GUIDE + "NNN", ref, GUIDE + "TGG", "TGG", 0, 0,
                    "X", var, "S01,S02")
    t2 = row_fields("chr21", anchor0, "+", GUIDE + "NNN", "-" + seq[anchor0:anchor0 + 22],
                    "-" + GUIDE[1:] + "TGG", "TGG", 0, 1, "RNA", var, "S01,S02")
    t3 = row_fields("chr21", anchor0, "+", GUIDE + "NNN", ref,
                    mark(GUIDE + "NNN", first), "TGG", 1, 0, "X", var, "S01,S02")
    return seq, rec, t1, t2, t3


class ReviewRound3(Base):
    # --- F1: a fold keeps the best alignment --------------------------------- #
    def test_a_fold_that_keeps_the_worse_alignment_is_caught(self):
        """F1: the perfect deletion target and the same PAM site's 1 mm + 1 DNA-bulge
        alignment have one variant and one carrier set, so the merge folds them into the
        best; a table holding only the worse one has lost the best (fewest block)."""
        case, info = world_deletion(self.tmp)
        seq, s = info["seq"], info["start"]
        crrna = GUIDE[0] + "-" + GUIDE[1:] + "NNN"
        worse = row_fields("chr21", s - 1, "+", crrna, seq[s - 1:s + 23],
                           seq[s - 1].lower() + GUIDE + "TGG", "TGG", 1, 1, "DNA",
                           info["var"], "S01,S02")
        write_results(case.res, [worse])
        rc, found = case.run(mode="completeness")
        self.assertEqual(rc, 1, found)
        self.assertTrue(any(f[1] == "missing-row" and "S01" == f[3] and
                            "fewest-mm+b block of 2" in f[2] for f in found), found)
        write_results(case.res, info["rows"])
        rc, found = case.run(mode="completeness")
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_fold_whose_cfd_block_keeps_a_lower_cfd_is_caught(self):
        """F1, highest-CFD block: a target with one mismatch at the guide's first base
        (CFD 1.0) and the same PAM's RNA-bulge alignment that skips that base (CFD 0.9),
        one edit each: the CFD block must keep the first."""
        seq = PRE + "G" + GUIDE[1:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        rec = (pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
        s = len(PRE)
        var = tok("chr21", rec)
        t1 = row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23],
                        "g" + GUIDE[1:] + "TGG", "TGG", 1, 0, "X", var, "S01,S02")
        t2 = row_fields("chr21", s + 1, "+", GUIDE + "NNN", "-" + seq[s + 1:s + 23],
                        "-" + GUIDE[1:] + "TGG", "TGG", 0, 1, "RNA", var, "S01,S02")
        case = Case(self.tmp, seq, [rec], [])
        budget = ("1", "0", "1", "1")
        # the CFD block holds the 0.9 alignment, the fewest block too (a tie on edits)
        write_results_cfd(case.res, [t2], [[0.9, 0.9, 0.9]])
        rc, found = case.run(mode="completeness", budget=budget)
        self.assertEqual(rc, 1, found)
        self.assertTrue(any(f[1] == "missing-row" and "highest-CFD block of 0.900" in f[2]
                            for f in found), found)
        write_results_cfd(case.res, [blocks_row([t1, t2, t2])], [[1.0, 0.9, 0.9]])
        rc, found = case.run(mode="completeness", budget=budget)
        self.assertEqual((rc, gating(found)), (0, []), found)

    # --- the fold compares rows on the haplotype (fix branch c3a75fe) ---------- #
    def test_bulged_alignments_inside_one_insertion_fold_on_the_haplotype(self):
        """The perfect target and its RNA-bulge alignment before the same PAM start on
        the 29th and 30th bases an insertion adds: both start at the anchor, their
        Cluster Positions differ by the bulge, and on the haplotype they are one site
        (c3a75fe folds them into the perfect row): a table with that row alone is
        right, as one with both (d7a6fff kept them apart)."""
        seq, rec, t1, t2, _t3 = world_inserted_targets(self.tmp)
        for rows in ([t1], [t1, t2]):
            case = Case(self.tmp, seq, [rec], rows)
            rc, found = case.run(mode="both", budget=("0", "0", "1", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (len(rows), found))

    def test_two_targets_inside_one_insertion_are_two_sites(self):
        """The added bases hold a second target (one mismatch) from their first base:
        both targets start at the anchor, at one Cluster Position, at two places on the
        haplotype. d7a6fff's merge folded the second into the first (one carrier set);
        c3a75fe lists both, and a table without the second misses it."""
        seq, rec, t1, _t2, t3 = world_inserted_targets(self.tmp, copy_mm=True)
        case = Case(self.tmp, seq, [rec], [t1])
        rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"))
        self.assertEqual(rc, 1, found)
        self.assertEqual({f[3] for f in found if f[1] == "missing-row"}, {"S01", "S02"},
                         found)
        case = Case(self.tmp, seq, [rec], [t1, t3])
        rc, found = case.run(mode="both", budget=("1", "0", "0", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_row_placed_past_the_anchor_stands_for_no_target(self):
        """2.8.1 placed a row whose first base is one an insertion adds at that base's
        place on the haplotype, the anchor plus its rank: its Cluster Position is then the
        target's cluster position on the haplotype. Its own target is not at its start, so
        it stands for no target in a fold, and completeness misses the target."""
        seq = LONG + PRE + GUIDE[6:] + "TGG" + POST + LONG
        anchor0 = len(LONG) + len(PRE) - 1
        rec = (anchor0 + 1, seq[anchor0], seq[anchor0] + "TTCAT" + GUIDE[:6],
               {"S04": "1|0"})
        right = row_fields("chr21", anchor0, "+", GUIDE + "NNN", seq[anchor0:anchor0 + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S04")
        s = anchor0 + 6   # GUIDE's first base is the 6th base added ("TTCAT" + ...)
        past = set_field(right, start=s, ref=seq[s:s + 23])
        for rows, want in (([right], []), ([past], ["S04"])):
            case = Case(self.tmp, seq, [rec], rows)
            rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
            self.assertEqual([f[3] for f in gating(found) if f[1] == "missing-row"], want,
                             found)

    # --- F2: a cluster is what chains ------------------------------------------ #
    def test_a_dropped_cluster_is_not_hidden_by_a_row_nearby(self):
        """F2: one row per cluster. Site A (the deletion's perfect target) is dropped and
        a perfect reference row sits 63 bp away, another cluster: A is missing."""
        seq_a = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG"
        seq = seq_a + LONG[:40] + GUIDE + "TGG" + POST
        rec = (len(PRE) + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
        s_a, s_b = len(PRE), len(seq_a) + 40
        row_a = row_fields("chr21", s_a, "+", GUIDE + "NNN", seq[s_a:s_a + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01,S02")
        row_b = set_field(row_fields("chr21", s_b, "+", GUIDE + "NNN", seq[s_b:s_b + 23],
                                     GUIDE + "TGG", "TGG", 0, 0, "X", "NA", "NA"),
                          origin="ref")
        for rows, want in (([row_a, row_b], 0), ([row_b], 1)):
            case = Case(self.tmp, seq, [rec], rows)
            run_markers(case.res, "population-level")
            rc, found = case.run(mode="completeness", budget=("0", "1", "1", "1"))
            self.assertEqual(rc, want, found)
            if want:
                self.assertIn("missing-locus", classes(gating(found)), found)

    def test_a_cluster_row_reached_by_chaining_alignments_satisfies(self):
        """F2: in a tandem repeat the guide's targets are 3 bp apart, so CRISPRme's merge
        chains a long run of them into one cluster, whose row can lie past --merge bp of
        a target: the repeat's own alignments link them. A stretch that breaks the run
        makes two clusters, and the row of one is not the other's."""
        g = "AGGAGGAGGAGGAGGAGGAG"
        left = "TTTCCTTTCATCTTTC"
        for broken in (False, True):
            run_a = "AGG" * 8
            gap = "TTCTTTCTTCTTTTCCTTTC" if broken else ""
            run_b = "AGG" * 14
            seq = LONG + left + run_a + gap + run_b + "TTTCTTTC" + POST
            s_row = len(LONG) + len(left)
            p = s_row + len(run_a) + len(gap) + 3 * 7 - 1   # the G ending run B's 7th unit
            assert seq[p:p + 4] == "GAGG", seq[p:p + 4]
            rec = (p + 1, "GAGG", "G", {"S01": "1|0"})
            row = set_field(row_fields("chr21", s_row, "+", g + "NNN",
                                       seq[s_row:s_row + 23], seq[s_row:s_row + 23],
                                       seq[s_row + 20:s_row + 23], 0, 0, "X", "NA", "NA",
                                       spacer=g + "NNN"), origin="ref")
            case = Case(self.tmp, seq, [rec], [row])
            run_markers(case.res, "population-level")
            rc, found = case.run(mode="completeness", guide=g, budget=("0", "0", "0", "0"))
            if broken:
                self.assertIn("missing-locus", classes(gating(found)), found)
            else:
                self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_cluster_chains_through_the_other_alleles_of_a_site(self):
        """F2: each allele of a multi-allelic site is a record, and each can put the
        guide's target a few bases from the others' (HPRC's alleles that fill an hg38 gap
        do): CRISPRme's cluster holds the rows of all of them, so the chain follows every
        record of the window. Three alleles put the target 2 bp apart; the cluster's one
        row is the third's, 4 bp from the first's target."""
        seq = LONG + PRE + "A" + "T" * 60 + POST + LONG
        a0 = len(LONG) + len(PRE)                   # 0-based: the shared first base
        recs, starts = [], []
        for k, shift in enumerate((0, 2, 4)):
            alt = "A" + "C" * (10 - shift) + GUIDE + "TGG" + "C" * 10
            recs.append((a0 + 1, "A" + "T" * 60, alt, {"S0%d" % (k + 1): "1|0"}))
            starts.append(a0 + 1 + 10 - shift)
        row = row_fields("chr21", starts[2], "+", GUIDE + "NNN",
                         seq[starts[2]:starts[2] + 23], GUIDE + "TGG", "TGG", 0, 0, "X",
                         tok("chr21", recs[2]), "S03")
        case = Case(self.tmp, seq, recs, [row])
        run_markers(case.res, "population-level")
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    # --- F3: the PUTATIVE rule as CRISPRme documents it ------------------------ #
    def _zero_cfd_world(self, bulge, listed):
        """A deletion target whose alignment holds a DNA bulge at guide position 7 (base
        ``bulge``) and an unphased SNV G>A at the column of the guide's C (position 16),
        S05 unphased at both; the row lists S05 with ``listed`` there."""
        g = GUIDE
        seq = PRE + g[:6] + bulge + g[6:10] + "T" + g[10:15] + "G" + g[16:] + "TGG" + POST
        dpos1, spos1 = len(PRE) + 11, len(PRE) + 18
        drec = (dpos1, g[9] + "T", g[9], {"S05": "0/1"})
        srec = (spos1, "G", "A", {"S05": "0/1"})
        s = len(PRE)
        crrna = g[:6] + "-" + g[6:] + "NNN"
        hap = g[:6] + bulge + g[6:15] + listed + g[16:] + "TGG"
        var = tok("chr21", drec) + ("," + tok("chr21", srec) if listed == "A" else "")
        rows = [row_fields("chr21", s, "+", crrna, seq[s:s + len(hap)], mark(crrna, hap),
                           "TGG", 1, 1, "DNA", var, "S05")]
        return Case(self.tmp, seq, [drec, srec], rows)

    def test_putative_bases_tie_when_the_alignment_scores_zero(self):
        """F3 (a): METHODS section 4 picks the alleles that maximize the alignment's CFD.
        A DNA bulge of G at guide position 7 scores 0, so every choice at the SNV scores 0
        and each is a strongest one (CRISPRme keeps the reference base); with a C there
        (0.059) the SNV's A (a match) is the strongest, and G is weaker."""
        m = C.Cfd(C.default_cfd_dir())
        self.assertEqual(m.factor("-", "G", 7), 0.0)
        self.assertGreater(m.factor("-", "C", 7), 0.0)
        for bulge, listed, ok in (("G", "G", True), ("G", "A", True), ("C", "A", True),
                                  ("C", "G", False)):
            for mode in ("soundness", "completeness"):
                rc, found = self._zero_cfd_world(bulge, listed).run(
                    mode=mode, budget=("1", "1", "0", "2"))
                if ok:
                    self.assertEqual((rc, gating(found)), (0, []), (bulge, listed, found))
                else:
                    self.assertEqual(rc, 1, (bulge, listed, found))
                    self.assertTrue(classes(gating(found)) & {"carrier-not-strongest",
                                                              "missing-putative-row"})

    def test_a_phased_second_indel_rules_a_slot_out_before_an_unphased_call(self):
        """F3 (b): S05's indel slot also holds a phased deletion of the guide's third
        base, inside the target, and an unphased SNV lies in the target: whatever the
        SNV's phase, S05's slot does not read the perfect target, and CRISPRme takes the
        slot off before resolving anything PUTATIVE (METHODS section 4)."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        main = (len(PRE) + 10, GUIDE[9] + "C", GUIDE[9], {"S05": "1|0", "S06": "1|0"})
        sec = (len(PRE) + 2, GUIDE[1] + GUIDE[2], GUIDE[1], {"S05": "1|0"})
        q = len(PRE) + 16
        snp = (q + 1, GUIDE[15], other(GUIDE[15]), {"S05": "0/1"})
        s = len(PRE)
        rows = [row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], GUIDE + "TGG",
                           "TGG", 0, 0, "X", tok("chr21", main), "S06")]
        case = Case(self.tmp, seq + LONG, [sec, main, snp], rows)
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
        # listing S05 there is wrong
        write_results(case.res, [set_field(rows[0], samples="S05,S06")])
        rc, found = case.run(mode="soundness", budget=("1", "0", "0", "1"))
        self.assertEqual(rc, 1, found)

    # --- F4: a long deletion's far junction ------------------------------------ #
    def test_a_long_deletion_reads_the_records_at_its_far_junction(self):
        """F4: the row names a 1,500-bp deletion and a cis SNV 4 bp past its far
        junction; the records there are read (soundness), and a work unit whose tile
        ends inside the deletion writes its carriers' haplotypes with them
        (completeness)."""
        import random
        rnd = random.Random(7)
        dele = "".join(rnd.choice("ACGT") for _ in range(1500))
        seq = PRE + GUIDE[:10] + dele + GUIDE[10:] + "TGG" + POST + LONG
        a = len(PRE) + 9
        drec = (a + 1, GUIDE[9] + dele, GUIDE[9], {"S01": "1|0"})
        q = a + 1 + 1500 + 4
        snp = (q + 1, GUIDE[14], other(GUIDE[14]), {"S01": "1|0"})
        s = len(PRE)
        tgt = GUIDE[:14] + other(GUIDE[14]).lower() + GUIDE[15:] + "TGG"
        rows = [row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], tgt, "TGG", 1, 0,
                           "X", "%s,%s" % (tok("chr21", drec), tok("chr21", snp)), "S01")]
        case = Case(self.tmp, seq, [drec, snp], rows)
        for extra in ((), ("--tile", "100")):
            rc, found = case.run(mode="both", budget=("1", "0", "0", "1"), extra=extra)
            self.assertEqual((rc, gating(found)), (0, []), (extra, found))
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"),
                             extra=("--tile", "100"))
        self.assertEqual([f[3] for f in gating(found)], ["S01"], found)

    # --- F5: a second indel that keeps the bases ------------------------------ #
    def test_a_second_indel_that_keeps_the_bases_is_demanded_exactly(self):
        """F5: S01's slot gains an AGG unit inside the target in an (AGG)n run that runs
        past the PAM: read from the target's first base, the bases are the window's, and
        CRISPRme documents keeping S01 on the row. S01 is its only carrier: a table
        without the row has lost it."""
        g = "ACTGATTCTGTACCAGGAGG"
        seq = LONG + PRE + g[:5] + "C" + g[5:] + "AGG" + "AGGAGGC" + POST
        s0 = len(LONG) + len(PRE)
        d1 = (s0 + 5, g[4] + "C", g[4], {"S01": "1|0"})
        d2 = (s0 + 15, g[13], g[13] + "AGG", {"S01": "1|0"})
        rows = [row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], g + "AGG", "AGG",
                           0, 0, "X", tok("chr21", d1), "S01", spacer=g + "NNN")]
        case = Case(self.tmp, seq, [d1, d2], rows)
        rc, found = case.run(mode="both", guide=g, budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", guide=g, budget=("0", "0", "0", "0"))
        self.assertEqual(rc, 1, found)
        self.assertIn("missing-row", classes(gating(found)), found)

    # --- F6: a limit only where a target reads -------------------------------- #
    def test_overlapping_alleles_with_no_target_are_no_limit(self):
        """F6: a deletion far from any target of the guide, with a SNV on a base it
        deletes on the same slot: the record has no target with or without the SNV, so
        nothing is left unjudged and no limit is listed."""
        seq = LONG + PRE + "ACGTTGCA" + POST + LONG
        p = len(LONG) + len(PRE) + 2
        recs = [(p, seq[p - 1:p + 2], seq[p - 1], {"S01": "1|0"}),
                (p + 1, seq[p], other(seq[p]), {"S01": "1|0"})]
        rc, found = Case(self.tmp, seq, recs, []).run(mode="completeness")
        self.assertEqual((rc, found), (0, []), found)

    # --- N1, N2, N3, N4, N5 ------------------------------------------------------ #
    def test_a_sites_only_vcf_checks_the_rows_themselves(self):
        """N1: an eight-column VCF (no FORMAT, no genotypes): its records are read, the
        row's strings, target and variants are checked; it lists no carrier, as a
        sites-only dataset's rows do (CHANGELOG 2.8.1), and a sample it lists is one
        no genotype carries (R5-F1)."""
        case, info = world_deletion(self.tmp)
        rec = info["rec"]
        write_sites_only_vcf(case.vcf, "chr21", [rec])
        write_results(case.res, [set_field(info["rows"][0], samples="NA")])
        rc, found = case.run(mode="both")
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual(sorted(f[1] for f in found), ["limit-no-genotypes"] * 4, found)
        write_results(case.res, [set_field(info["rows"][0], samples="S01,S02")])
        rc, found = case.run(mode="both")
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-without-indel", "S01,S02")}, found)
        write_results(case.res, [set_field(info["rows"][0], start=info["start"] + 1,
                                           samples="NA")])
        rc, found = case.run()
        self.assertIn("start-not-on-variant-sequence", classes(gating(found)), found)

    def test_a_tie_the_documented_rule_leaves_is_listed_not_silent(self):
        """N2: without CFD (none applies) every PUTATIVE choice the guide's base does not
        decide ties. S01 and S02 carry the MNV (CRISPRme resolves MNV carriers PUTATIVE);
        S01 also carries, phased on its other slot, a SNV whose two alleles both
        mismatch the guide. S01's haplotype holds the reference base there, the row of
        that target lists S02 only, and S01 is on the SNV's row: CRISPRme's rule ties the
        two, so S01's absence is a limit line, never silence."""
        x = other(GUIDE[2])
        y = next(b for b in "ACGT" if b not in (x, GUIDE[2]))
        seq = PRE + GUIDE[:2] + x + GUIDE[3:8] + "AC" + GUIDE[10:] + "TGG" + POST
        mnv = (len(PRE) + 9, "AC", GUIDE[8:10], {"S01": "1|0", "S02": "1|0"})
        snv = (len(PRE) + 3, x, y, {"S01": "0|1"})
        s = len(PRE)
        rows = [row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23],
                           mark(GUIDE + "NNN", GUIDE[:2] + b + GUIDE[3:] + "TGG"), "TGG",
                           1, 0, "X", var, who)
                for b, var, who in ((x, tok("chr21", mnv), "S02"),
                                    (y, "%s,%s" % (tok("chr21", mnv), tok("chr21", snv)),
                                     "S01"))]
        case = Case(self.tmp, seq, [snv, mnv], rows)
        nocfd = ("--cfd-tables", os.path.join(self.tmp, "none"))
        rc, found = case.run(budget=("1", "0", "0", "1"), extra=nocfd)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "limit-putative-carrier" and f[3] == "S01" and
                            "tie" in f[2] for f in found), found)
        rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"), extra=nocfd)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "limit-putative-carrier" and "S01" in f[3]
                            for f in found), found)

    def test_another_guides_row_satisfies_no_demand(self):
        """N3: two guides; the deletion's target is perfect for the first and one
        mismatch for the second, and the table reports it only for the second."""
        g2 = GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:]
        case, info = world_deletion(self.tmp)
        s, seq = info["start"], info["seq"]
        row2 = row_fields("chr21", s, "+", g2 + "NNN", seq[s:s + 23],
                          mark(g2 + "NNN", GUIDE + "TGG"), "TGG", 1, 0, "X", info["var"],
                          "S01,S02", spacer=g2 + "NNN")
        write_results(case.res, [row2])
        rc, found = case.run(mode="both", budget=("1", "0", "0", "1"),
                             extra=("--guide", g2))
        self.assertEqual(rc, 1, found)
        self.assertEqual({f[3] for f in gating(found) if f[1] == "missing-row"},
                         {"S01", "S02"}, found)
        write_results(case.res, [row2, info["rows"][0]])
        rc, found = case.run(mode="both", budget=("1", "0", "0", "1"),
                             extra=("--guide", g2))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_the_snp_row_report_names_a_carrier_it_leaves_out(self):
        """N4: a SNP row listing S01 while S02's phased slot holds the same SNV in the
        target: reported (report-only), and never gating."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        snp = len(PRE) + 6
        recs = [(snp, other(GUIDE[5]), GUIDE[5], {"S01": "1|0", "S02": "0|1"})]
        s = len(PRE)
        rows = [row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], GUIDE + "TGG",
                           "TGG", 0, 0, "X", tok("chr21", recs[0]), "S01")]
        rc, found = Case(self.tmp, seq, recs, rows).run(extra=("--snp-rows", "report"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "carrier-missing" and f[3] == "S02" for f in found),
                        found)

    def test_the_total_edits_cap_is_the_runs_own(self):
        """N5: without --max-total-edits the cap is the one the run searched with (its
        .Params.txt), else complete-search's 6."""
        case, info = world_deletion(self.tmp)
        s, seq = info["start"], info["seq"]
        two = row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23],
                         mark(GUIDE + "NNN", GUIDE[:3] + other(GUIDE[3]) + GUIDE[4:6] +
                              other(GUIDE[6]) + GUIDE[7:] + "TGG"), "TGG", 2, 0, "X",
                         info["var"], "S01,S02")
        write_results(case.res, [two])
        argv = ["--results", case.res, "--reference", os.path.join(self.tmp, "Genomes"),
                "--vcf", case.vcf, "--pam", case.pam, "--mm", "4", "--bDNA", "1",
                "--bRNA", "1", "--mode", "soundness", "--chrom", "chr21", "--guide", GUIDE,
                "--out", os.path.join(self.tmp, "o.tsv")]
        self.assertNotIn("budget-exceeded", open_findings(C.main(argv), self.tmp))
        with open(os.path.join(case.res, ".Params.txt"), "w") as fh:
            fh.write("Mismatches\t4\nDNA\t1\nRNA\t1\nMax_total_edits\t1\n")
        self.assertIn("budget-exceeded", open_findings(C.main(argv), self.tmp))
        # R7-F4's sweep: a cap other than the run's is not the search's: refused
        self.assertEqual(C.main(argv + ["--max-total-edits", "6"]), 2)
        self.assertEqual(C.main(argv + ["--max-total-edits", "1"]), 1)


# --------------------------------------------------------------------------- #
# review round 4: one target's alignments fold, MNVs CRISPRme cannot see, the panel
# --------------------------------------------------------------------------- #
README = "CTAACAGTTGCTTTTATCAC"
# hg38 chr22:17580151-17580290 (1-based) and the deletion chr22:17580229 ATCG>A: the
# README guide's target at chr22:17580213 (+) of fix282_NRG_ps line 404 (24 carriers),
# whose alignments include 3 mm + 2 bulges (CFD 0.000967) and 4 + 2 (CFD 0.001505)
REAL_404 = (17580150, "TATATGTTATTGTTCAACATAGTTTAATTGTTTGGGTTCATCTTCTCTAACTAC"
                      "AAAATTAAGTTTACAGTTCTCTTTATCGTCCCAAGTTTCTTCCTTTGTGAAAGGCC"
                      "CACTATCATCTACTCTTACCACTTCCAGGT",
            (17580229, "ATCG", "A"), 17580213, "+",
            ("CTAACAGTTGC-TTTTATCACNNN", "tTtACAGTT-CTcTTTATCcCAAG"),    # 4 + 2, 0.002
            ("CTAACAGTTGCT-TTTATCACNNN", "tTtACAGTT-CTCTTTATCcCAAG"))    # 3 + 2, 0.001
# hg38 chr22:40428111-40428250 and the deletion chr22:40428184 TA>T: the target at
# chr22:40428174 (-) of line 87 (23 carriers): two 4 + 2 alignments, CFD 0.032 and 0.012
REAL_87 = (40428110, "TTCTGTGAGCCATCCTAGCAAATTATCAAGCCCAAGGAGAGGGTATGGGAACCCC"
                     "AATCTATAGCTAGTTGGTTAAAAGCAAATGTAAAATAACCTGGGGCTTGGGACCAG"
                     "CATGGGAAGTGGGGCTAAGTCTTGTGGGA",
           (40428184, "TA", "T"), 40428174, "-",
           ("CTAACAGTTGCTTTTATCA-CNNN", "tT-ACAtTTGCTTTaAcCAACTAG"),     # 0.032
           ("CTAACAGTTGCTTTTATC-ACNNN", "-TtACAtTTGCTTTaAcCAACTAG"))     # 0.012


def real_row(seq, start, strand, crrna, alt, var, samples, spacer=README + "NNN"):
    """A block of a real table: its counts and REF text (the reference read from the start,
    in the target's orientation, gapped as the target) from the strings."""
    n = len(alt.replace("-", ""))
    span = seq[start:start + n]
    text = list(span if strand == "+" else C.rc(span))
    for i, ch in enumerate(alt):
        if ch == "-":
            text.insert(i, "-")
    mm, dn, rn, _p = C.recount(crrna, alt, 3, False)
    bt = ",".join([t for t, c in (("RNA", rn), ("DNA", dn)) if c]) or "X"
    return row_fields("chr22", start, strand, crrna, "".join(text), alt,
                      alt.replace("-", "")[-3:].upper(), mm, dn + rn, bt, var, samples,
                      spacer=spacer)


class ReviewRound4(Base):
    # --- R4-F1: two alignments of one target fold into the best ---------------- #
    def _real(self, world):
        w0, seq, (pos1, ref, alt), start, strand, a, b = world
        rec = (pos1 - w0, ref, alt, {"S01": "1|0", "S02": "0|1"})
        var = tok("chr22", rec)   # the window's own coordinates
        case = Case(self.tmp, seq, [rec], [], chrom="chr22")
        rows = [real_row(seq, start - w0, strand, x[0], x[1], var, "S01,S02")
                for x in (a, b)]
        return case, start - w0, rows

    def _target_lines(self, found, start):
        return [f for f in gating(found) if "this target at %d " % start in f[2]]

    def test_a_fewest_block_holding_the_worse_alignment_of_one_target_is_caught(self):
        """R4-F1 (mutation M1 of the review, on fix282_NRG_ps line 404): the fewest-mm+b
        block of the target's line takes the CFD block's 4 mm + 2 bulge alignment, so
        its 3 + 2 alignment is in no table. The sample is still listed for the target,
        but the rows of its site with its carriers are worse than it."""
        case, s, (a, b) = self._real(REAL_404)
        write_results_cfd(case.res, [blocks_row([a, b, b])], [[0.002, 0.001, 0.001]])
        rc, found = case.run(mode="completeness", guide=README)
        self.assertEqual(self._target_lines(found, s), [], found)
        write_results_cfd(case.res, [blocks_row([a, a, a])], [[0.002, 0.002, 0.002]])
        rc, found = case.run(mode="completeness", guide=README)
        hits = self._target_lines(found, s)
        self.assertEqual(rc, 1, found)
        self.assertEqual(sorted((f[1], f[3]) for f in hits),
                         [("missing-best-alignment", "S01"),
                          ("missing-best-alignment", "S02")], hits)
        self.assertIn("fewest-mm+b block of 6 (mm+bulges), the target 5", hits[0][2])

    def test_a_cfd_block_holding_the_lower_cfd_alignment_of_one_target_is_caught(self):
        """R4-F1 (mutation M2b, line 87): the highest-CFD block takes the fewest block's
        alignment of the same target, CFD 0.012 for 0.032."""
        case, s, (a, b) = self._real(REAL_87)
        write_results_cfd(case.res, [blocks_row([a, b, b])], [[0.032, 0.012, 0.012]])
        rc, found = case.run(mode="completeness", guide=README)
        self.assertEqual(self._target_lines(found, s), [], found)
        write_results_cfd(case.res, [blocks_row([b, b, b])], [[0.012, 0.012, 0.012]])
        rc, found = case.run(mode="completeness", guide=README)
        hits = self._target_lines(found, s)
        self.assertEqual({(f[1], f[3]) for f in hits},
                         {("missing-best-alignment", "S01"),
                          ("missing-best-alignment", "S02")}, hits)
        self.assertIn("highest-CFD block of 0.012", hits[0][2])

    def test_a_cfd_block_lower_only_below_the_old_tolerance_is_caught(self):
        """R4-F4 (mutation M2, line 404): the CFD block takes the 3 + 2 alignment, CFD
        0.000967 (printed 0.001) for 0.001505 (0.002). The merge sorts on the three
        decimals the analysis writes, so it would have kept 0.002: the checker compares
        each row's CFD, recomputed from its own strings, at those three decimals (a
        tolerance of 0.0006 on the printed values let it pass)."""
        case, s, (a, b) = self._real(REAL_404)
        write_results_cfd(case.res, [blocks_row([b, b, b])], [[0.001, 0.001, 0.001]])
        rc, found = case.run(mode="completeness", guide=README)
        hits = self._target_lines(found, s)
        self.assertEqual({(f[1], f[3]) for f in hits},
                         {("missing-best-alignment", "S01"),
                          ("missing-best-alignment", "S02")}, hits)
        self.assertIn("highest-CFD block of 0.001 (0.000967", hits[0][2])
        # a printed CFD the row's own strings do not give does not hide it
        write_results_cfd(case.res, [blocks_row([b, b, b])], [[0.002, 0.001, 0.001]])
        rc, found = case.run(mode="completeness", guide=README)
        self.assertTrue(self._target_lines(found, s), found)

    def test_a_cluster_row_is_scored_by_its_own_alignment(self):
        """R4-F4, one row per cluster (line 87's target, population-level): the merge
        sorts on the CFD of a row's own alignment (CFD_score_ALT); the results' CFD_score
        prints the higher of the REF and ALT scores, so a cluster row whose REF scores
        0.032 and whose own alignment 0.012 is worse than the target's 0.032."""
        case, s, (a, b) = self._real(REAL_87)
        run_markers(case.res, "population-level")
        for row, want in ((a, []), (b, ["missing-locus"])):
            write_results_cfd(case.res, [blocks_row([row, row, row])], [[0.032] * 3])
            rc, found = case.run(mode="completeness", guide=README)
            self.assertEqual([f[1] for f in gating(found)
                              if "Cluster Position %d " % s in f[2]], want, found)

    def test_only_the_worst_alignment_of_the_perfect_target_is_caught(self):
        """R4-F1 (the review's isolated case): every target of the deletion's site is
        reported by its best alignment except the perfect one, reported only by its
        worst (6 edits); the right table, with its 0-edit alignment, passes."""
        case, info = world_deletion(self.tmp)
        seq, s = info["seq"], info["start"]
        probe = GUIDE + "TGG"
        als = sorted(C.align_pairs(GUIDE + "NRG", probe, 4, 1, 1, 6, 3, False),
                     key=lambda x: (x[2] + x[3] + x[4], x))
        worst = als[-1]
        self.assertEqual(worst[2] + worst[3] + worst[4], 6)
        row = real_row(seq, s, "+", worst[0][:-3] + "NNN", mark(worst[0], worst[1]),
                       info["var"], "S01,S02", spacer=GUIDE + "NNN")
        row[1] = "chr21"
        write_results(case.res, [row])
        rc, found = case.run(mode="completeness")
        hits = [f for f in gating(found) if "this target at %d " % s in f[2]]
        self.assertEqual(sorted((f[1], f[3]) for f in hits),
                         [("missing-best-alignment", "S01"),
                          ("missing-best-alignment", "S02")], found)
        write_results(case.res, info["rows"])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_putative_target_is_held_to_its_own_carriers_alignments(self):
        """The oracle's site at chr21:8000 (synthetic genome, CRISPRme's own rows, S19 and
        S20 renamed S01 and S02): both are unphased for the deletion, S01 also for the
        SNV that gives the guide's A at the target's sixth base. Under the alignment
        without bulges S01's strongest target is the one with that A (0 edits, listed);
        the one with C is its strongest only under the bulged alignments (2 edits), which
        list it, while the 1-mismatch alignment is S02's alone. A demand that pooled
        S02's alignments held S01 to 1 edit, and which alignment it kept depended on the
        hash seed: every seed must give the same, right answer."""
        import subprocess
        seq = ("CGAGTCGTAAATGACGTGACTTGTCCCATTAATCACGTATTTGTGACCGCGAGGCGTCGAGTTGGCTGTT"
               "AGATCGCCGCCCCTCGAATTTAGTGAAATAACTGACATCTCGTAAGCAGGCTGGTGCAGTGGAACCGAACTCTC"
               "GCACCCAATGATGTATATGAGCTACACCATACCATCATTACTACATATCATCTTAT")
        snv = (106, "C", "A", {"S01": "0/1"})
        dele = (110, "TC", "T", {"S01": "0/1", "S02": "0/1"})
        d, v = tok("chr21", dele), tok("chr21", snv)
        r0 = "ACTGAcATCTcgtAagcaGgCTG"

        def blk(start, crrna, ref, alt, mm, bul, bt, var, who):
            return row_fields("chr21", start, "+", crrna, ref, alt, alt.replace("-", "")[-3:],
                              mm, bul, bt, var, who)
        perfect = blk(100, GUIDE + "NNN", r0, GUIDE + "TGG", 0, 0, "X", v + "," + d, "S01")
        one = blk(100, GUIDE + "NNN", r0, "ACTGAcATCTGTAAGCAGGCTGG", 1, 0, "X", d, "S02")
        snp1 = blk(100, "ACTGAAATCT-GTAAGCAGGCNNN", "ACTGAcATCTCGTAAGCAGGCTGG",
                   "ACTGAAATCTCGTAAGCAGGCTGG", 0, 1, "DNA", v, "S01")
        snp2 = blk(100, "ACTGAAATC-TGTAAGCAGGCNNN", "ACTGAcATCTcGTAAGCAGGCTGG",
                   "ACTGAAATCTcGTAAGCAGGCTGG", 1, 1, "DNA", v, "S01")
        bulged = blk(100, "A-CTGAAATCTGTAAGCAGGCNNN", "AC-TGAcATCTcgtAagcaGgCTG",
                     "AC-TGAcATCTGTAAGCAGGCTGG", 1, 2, "RNA,DNA", d, "S01,S02")
        two = blk(100, "ACTGAA-ATCTGTAAGCAGGCNNN", "ACTG-ACATCTcgtAagcaGgCTG",
                  "ACTG-ACATCTGTAAGCAGGCTGG", 0, 2, "RNA,DNA", d, "S01,S02")
        left = blk(99, "A-CTGAAATCTGTAAGCAGGCNNN", "AACTGAcATCTcgtAagcaGgCTG",
                   "AACTGAcATCTGTAAGCAGGCTGG", 1, 1, "DNA", d, "S01,S02")
        case = Case(self.tmp, seq, [snv, dele], [perfect])
        write_results(case.res, [blocks_row([one, snp1, one]), blocks_row([bulged, one, snp1]),
                                 blocks_row([snp2, two, left])],
                      name="G+NRG_hg38+v_4+2_all_results_with_alternative_alignments.tsv")
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                              "scripts", "check_indel_rows.py")
        argv = ["--results", case.res, "--reference", os.path.join(self.tmp, "Genomes"),
                "--vcf", case.vcf, "--pam", case.pam, "--guide", GUIDE, "--mm", "4",
                "--bDNA", "1", "--bRNA", "1", "--max-total-edits", "6", "--chrom", "chr21",
                "--mode", "completeness"]
        answers = set()
        run_markers(case.res, "per-sample")
        for seed in range(6):
            out = os.path.join(self.tmp, "seed%d.tsv" % seed)
            env = dict(os.environ, PYTHONHASHSEED=str(seed))
            p = subprocess.run([sys.executable, script] + argv + ["--out", out], env=env,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with open(out) as fh:
                answers.add((p.returncode, fh.read()))
        self.assertEqual(len(answers), 1, answers)
        rc, text = answers.pop()
        self.assertEqual(rc, 0, text)

    def test_a_putative_alternative_is_compared_with_its_own_best(self):
        """R4-F1: a PUTATIVE target the rule ties with another is compared with its own
        alignments' best, not the best of its alternatives."""
        t1 = C.Target.__new__(C.Target)
        t2 = C.Target.__new__(C.Target)
        for t, p in ((t1, "AAAA"), (t2, "AAAC")):
            t.strand, t.start, t.probe, t.crrna, t.aligned = "+", 10, p, "AAAA", p
        x = frozenset(["chr21_11_A_AT"])
        d = C.Demand("putative", [(t1, x), (t2, x)], C.edits(1, 0), 0.9, 0)
        d.add_alt(t1, x, C.edits(1, 0), 0.9)
        d.add_alt(t2, x, C.edits(3, 0), 0.2)
        d.add_alt(t2, x, C.edits(2, 1), 0.1)
        self.assertEqual((d.best_of(t1, x).best["mm+bulges"], d.best_of(t1, x).cfd),
                         (1, 0.9))
        self.assertEqual((d.best_of(t2, x).best["mm+bulges"], d.best_of(t2, x).cfd),
                         (3, 0.2))
        self.assertEqual(d.best_of(t2, x).best["bulges"], 0)
        self.assertTrue(C.cfd_below(0.000967, 0.001505))
        self.assertFalse(C.cfd_below(0.0024, 0.0016))     # both print 0.002
        self.assertFalse(C.cfd_below(0.0015, 0.0015 + 1e-15))
        # products taken in another order can print a x.xxx5 value either way (genome-wide,
        # 14 of 2.4 million REF and SNP blocks print 0.010 for 0.0105 and the like)
        self.assertFalse(C.cfd_below(0.0105 - 1e-12, 0.0108))
        self.assertFalse(C.cfd_below(0.0102, 0.0105 + 1e-12))
        self.assertTrue(C.cfd_below(0.0104, 0.0106))

    # --- R4-F2: an MNV CRISPRme cannot see ------------------------------------- #
    def _mnv_world(self, second, listed):
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dele = (len(PRE) + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
        p2 = len(PRE) + 1 + 14 + 1   # guide bases 15 and 16 (1-based POS)
        r2 = seq[p2 - 1:p2 + 1]
        if second == "mnv":
            rec2 = (p2, r2, other(r2[0]) + other(r2[1]), {"S01": "1|0"})
        else:
            rec2 = (p2, r2[0], other(r2[0]), {"S01": "1|0"})
        s = len(PRE)
        row = row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], GUIDE + "TGG",
                         "TGG", 0, 0, "X", tok("chr21", dele), listed)
        return Case(self.tmp, seq, [dele, rec2], [row])

    def test_a_carrier_listed_for_a_sequence_its_mnv_changes_is_caught(self):
        """R4-F2: S01 carries the deletion that makes the target and, phased on the same
        slot, an MNV on two of its bases. CRISPRme reads no MNV beside the record, but
        nothing documents listing S01 for the sequence without it: gating, as for an
        SNV there (the review found it passing as the two-indel limit)."""
        for second in ("mnv", "snv"):
            rc, found = self._mnv_world(second, "S01,S02").run(mode="soundness")
            self.assertEqual(rc, 1, (second, found))
            self.assertEqual({(f[1], f[3]) for f in gating(found)},
                             {("variant-missing", "S01")}, (second, found))
            self.assertNotIn("limit-two-indels", classes(found), found)
            if second == "mnv":
                self.assertIn("an MNV on their slot changes the target", found[0][2])

    def test_a_target_that_needs_the_mnv_is_the_two_indel_limit(self):
        """R4-F2: the target S01's slot holds, with the MNV, needs two indel-path records
        (CRISPRme's window holds one): not demanded as a row (limit-two-indels); the row
        listing S02 alone passes both checks."""
        case = self._mnv_world("mnv", "S02")
        rc, found = case.run(mode="both", budget=("2", "0", "0", "2"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-two-indels", classes(found), found)

    # --- R4-F3: the samples list ----------------------------------------------- #
    def test_the_panel_is_the_samples_list(self):
        """R4-F3: CRISPRme lists no carrier outside the samplesID panel. S02 carries the
        deletion but is not in the panel: a row listing S01 alone is right, one listing
        S02 is not; the panel is --samples, else the run's .sampleID.txt."""
        case, info = world_deletion(self.tmp)
        only1 = set_field(info["rows"][0], samples="S01")
        write_results(case.res, [only1])
        rc, found = case.run(mode="both", budget=("0", "0", "0", "0"))
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-missing", "S02"), ("missing-row", "S02")}, found)
        panel = os.path.join(self.tmp, "panel.txt")
        with open(panel, "w") as fh:
            fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
            fh.write("".join("%s\tP\tSP\tfemale\n" % s for s in SAMPLES if s != "S02"))
        rc, found = case.run(mode="both", budget=("0", "0", "0", "0"),
                             extra=["--samples", panel])
        self.assertEqual((rc, gating(found)), (0, []), found)
        shutil.copy(panel, os.path.join(case.res, ".sampleID.txt"))
        rc, found = case.run(mode="both", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, info["rows"])   # S01,S02
        rc, found = case.run(mode="soundness")
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-not-in-panel", "S02")}, found)
        with open(panel, "w") as fh:
            fh.write("NA12878\n")
        os.remove(os.path.join(case.res, ".sampleID.txt"))   # (R7: else refused for that)
        rc, _found = case.run(mode="soundness", extra=["--samples", panel])
        self.assertEqual(rc, 2)

    # --- R4-F5, F6, F7 ---------------------------------------------------------- #
    def test_an_alt_block_naming_no_variant_is_reported(self):
        """R4-F5: an ALT-origin block that lost its variant list is not skipped."""
        case, info = world_deletion(self.tmp)
        lost = set_field(info["rows"][0], var="NA")
        write_results(case.res, [info["rows"][0], lost])
        rc, found = case.run(mode="soundness")
        self.assertEqual(rc, 1, found)
        self.assertEqual([f[1] for f in gating(found)], ["alt-without-variant"] * 3, found)
        ref_row = set_field(lost, origin="ref", samples="NA")
        write_results(case.res, [info["rows"][0], ref_row])
        rc, found = case.run(mode="soundness")
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_the_putative_cap_is_a_coverage_gap_of_the_checker(self):
        """R4-F6: a carrier with more than 10 calls of unknown phase within a target's
        reach of a record: its targets that touch one are not enumerated, which is this
        checker's own coverage gap, listed under its own class, not a CRISPRme limit."""
        case, info = world_deletion(self.tmp)
        seq = info["seq"]
        recs = [info["rec"]] + [(q + 1, seq[q], other(seq[q]), {"S01": "0/1"})
                                for q in range(36, 47)]
        case = Case(self.tmp, seq, recs, info["rows"])
        rc, found = case.run(mode="completeness")
        caps = [f for f in found if f[1] == "coverage-putative-cap"]
        self.assertEqual([(f[0], f[3]) for f in caps], [("no", "S01")], found)
        self.assertIn("coverage gap, not a CRISPRme limit", caps[0][2])
        self.assertFalse(any(f[1].startswith("limit-putative-cap") for f in found))

    def test_the_cluster_chain_reads_only_the_alleles_crisprme_keeps(self):
        """R4-F7: the enriched genome keeps the last record of a split site, so the
        cluster chain reads no earlier one."""
        tmp = tempfile.mkdtemp(prefix="chkidl_")
        try:
            seq = "ACGT" * 20
            vcf = os.path.join(tmp, "v.vcf.gz")
            write_vcf(vcf, "chr21", [(5, "A", "C", {"S01": "1|0"}),
                                     (5, "A", "G", {"S02": "0|1"}),
                                     (9, "A", "T", {"S01": "1|0"})])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
            ctx = argparse_ns_ctx()     # population-level SNP rows: every allele kept
            union = [x for x in C.snp_row_sets(win, ctx, 0, 20) if len(x[4]) > 1]
            self.assertEqual(len(union), 1, union)
            self.assertEqual(union[0][4], frozenset("AG"))
            self.assertEqual(union[0][8], frozenset("AT"))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class ReviewRound5(Base):
    # --- R5-F1: a sites-only VCF beside a genotyped one ------------------------ #
    def _mixed(self, deletion_in_genotyped, sites_first=False):
        """world_deletion with two VCFs, as a search of the default index with hg38_mega
        reads them: a genotyped one (the eight samples; the deletion's carriers S01 and
        S02 when ``deletion_in_genotyped``, else only an unrelated SNV of S03) and a
        sites-only one holding the deletion. Returns (case, info, the other VCF's
        argv)."""
        case, info = world_deletion(self.tmp)
        seq = info["seq"]
        recs = [(5, seq[4], other(seq[4]), {"S03": "0|1"})]
        if deletion_in_genotyped:
            recs.append(info["rec"])
        write_vcf(case.vcf, "chr21", recs)
        sites = os.path.join(self.tmp, "sites.vcf.gz")
        write_sites_only_vcf(sites, "chr21", [info["rec"]])
        if sites_first:
            genotyped, case.vcf = case.vcf, sites
            return case, info, ["--vcf", genotyped]
        return case, info, ["--vcf", sites]

    def test_a_sites_only_record_beside_a_genotyped_vcf_lists_no_carrier(self):
        """R5-F1 (the review's adv_r5a case 1): the record is in the sites-only VCF alone,
        so its row is the sites-only dataset's and lists no carrier, which CRISPRme
        documents (CHANGELOG 2.8.1): the row is checked, a non-gating limit-no-genotypes
        line says its carriers are not, and completeness counts the record. Whichever
        VCF comes first. A sample it lists carries nothing; a wrong start still fails."""
        for sites_first in (False, True):
            case, info, extra = self._mixed(False, sites_first)
            write_results(case.res, [set_field(info["rows"][0], samples="")])
            rc, found = case.run(mode="both", extra=extra)
            self.assertEqual((rc, gating(found)), (0, []), (sites_first, found))
            lines = [f for f in found if f[1] == "limit-no-genotypes"]
            self.assertEqual(len(lines), 4, found)   # three blocks, one completeness count
            self.assertTrue(any("1 indel/MNV record(s) of chr21 are held only by "
                                "sites-only VCF(s)" in f[2] for f in lines), lines)
            write_results(case.res, [set_field(info["rows"][0], samples="S01")])
            rc, found = case.run(extra=extra)
            self.assertEqual({(f[1], f[3]) for f in gating(found)},
                             {("carrier-without-indel", "S01")}, (sites_first, found))
            self.assertIn("only sites-only VCF(s) hold the record", gating(found)[0][2])
            write_results(case.res, [set_field(info["rows"][0], samples="",
                                               start=info["start"] + 1)])
            rc, found = case.run(extra=extra)
            self.assertIn("start-not-on-variant-sequence", classes(gating(found)), found)
            self.assertNotIn("limit-no-genotypes", classes(found), found)

    def test_a_record_both_vcfs_hold_keeps_its_genotyped_carriers(self):
        """R5-F1: a record the genotyped VCF holds too. The genotyped dataset's row lists
        its carriers, and they are checked (one left off is missing). The sites-only
        dataset's carrier-less row of it is a limit, and completeness still demands the
        carriers' own row: a table with only the carrier-less row misses them."""
        for sites_first in (False, True):
            case, info, extra = self._mixed(True, sites_first)
            rc, found = case.run(mode="both", extra=extra)
            self.assertEqual((rc, found), (0, []), (sites_first, found))
            write_results(case.res, [set_field(info["rows"][0], samples="S01")])
            rc, found = case.run(extra=extra)
            self.assertEqual({(f[1], f[3]) for f in gating(found)},
                             {("carrier-missing", "S02")}, (sites_first, found))
            bare = set_field(info["rows"][0], samples="NA")
            write_results(case.res, [info["rows"][0], bare])
            rc, found = case.run(mode="both", extra=extra)
            self.assertEqual((rc, gating(found)), (0, []), found)
            self.assertTrue(all("demanded by completeness" in f[2] for f in found), found)
            write_results(case.res, [bare])
            rc, found = case.run(mode="completeness", extra=extra)
            self.assertEqual({(f[1], f[3]) for f in gating(found)
                              if "this target at %d " % info["start"] in f[2]},
                             {("missing-row", "S01"), ("missing-row", "S02")}, found)

    def test_a_carrier_less_row_of_a_genotyped_record_still_fails(self):
        """R5-F1: the decision is the record's: a row that lists no carrier for a record
        only the genotyped VCF holds is no sites-only dataset's row."""
        case, info = world_deletion(self.tmp)
        sites = os.path.join(self.tmp, "sites.vcf.gz")
        write_sites_only_vcf(sites, "chr21", [(5, info["seq"][4], other(info["seq"][4]))])
        write_results(case.res, [set_field(info["rows"][0], samples="")])
        rc, found = case.run(extra=["--vcf", sites])
        self.assertEqual(sorted(f[1] for f in gating(found)), ["no-carrier-listed"] * 3,
                         found)

    def test_an_earlier_record_of_a_split_site_is_a_note(self):
        """R5-F1 follow-up: a row naming an earlier record of a split site, which the
        documented limit says the stores never hold. The record is real and the target
        reads with it, so naming it is no wrong row: a non-gating note (variant-shadowed),
        a genotyped VCF's and, on an indel row, a sites-only VCF's alike -- the fake-indel
        windows hold the last record whose FILTER passes (R10-F1: the mega's indel rows
        name it at every split site they name on chr22, the last PASS one at
        chr22:25093750, whose last record C>T is LOW_VQSLOD; its SNP rows, which name
        earlier records, are no indel rows)."""
        case, info = world_indel_snp(self.tmp)
        recs = info["recs"]
        snv = recs[0]
        later = (snv[0], snv[1], other(snv[2]) if other(snv[2]) != snv[1] else "T", {})
        geno = [snv, later, recs[1]]
        write_vcf(case.vcf, "chr21", geno)
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual(sorted(f[1] for f in found), ["variant-shadowed"] * 3, found)
        write_sites_only_vcf(case.vcf, "chr21", geno)
        write_results(case.res, [set_field(info["rows"][0], samples="")])
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual(sorted(f[1] for f in found),
                         ["limit-no-genotypes"] * 3 + ["variant-shadowed"] * 3, found)

    def test_a_sites_only_datasets_snp_row_is_out_of_the_snp_report(self):
        """R5-F1 follow-up: beside a genotyped dataset, a sites-only one's SNP rows are
        population-level rows (a worst-case union, out of scope): the SNP-row report says
        so in one line per block, and still reads the genotyped dataset's rows."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        snp = len(PRE) + 6
        rec = (snp, other(GUIDE[5]), GUIDE[5], {"S01": "1|0"})
        start = len(PRE)
        bare = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                          GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "")
        case = Case(self.tmp, seq, [(5, seq[4], other(seq[4]), {"S03": "0|1"})], [bare])
        sites = os.path.join(self.tmp, "sites.vcf.gz")
        write_sites_only_vcf(sites, "chr21", [rec])
        rc, found = case.run(extra=("--snp-rows", "report", "--vcf", sites))
        self.assertEqual(sorted(f[1] for f in found), ["limit-no-genotypes"] * 3, found)
        self.assertIn("population-level one", found[0][2])
        write_vcf(case.vcf, "chr21", [rec])   # the genotyped dataset holds it, S01 carries
        listed = set_field(bare, samples="S01")
        write_results(case.res, [listed])
        rc, found = case.run(extra=("--snp-rows", "report", "--vcf", sites))
        self.assertEqual((rc, found), (0, []), found)

    # --- R5-F2: "no worse" is the merge's whole sort key ------------------------ #
    def _two_targets_one_site(self):
        """The review's adv_r5b world: one deletion's site holds a 24-mer (0 mm + 1 DNA
        bulge: 1 edit, 0 mm) and a 23-mer (1 mm: 1 edit, 1 mm) with the same variants and
        carriers, so the merge folds them, and its fewest file, sorted on mm+bulges and
        then mm (crisprme.py's default --sorting-criteria), keeps the 24-mer."""
        core = "AC" + "C" + GUIDE[2:12] + "G" + GUIDE[12:] + "TGG"
        seq = PRE + core + POST
        dpos1 = len(PRE) + 13
        anchor0 = dpos1 - 1
        rec = (dpos1, seq[anchor0] + "G", seq[anchor0], {"S01": "1|0", "S02": "0|1"})
        var = tok("chr21", rec)
        hap = seq[:anchor0 + 1] + seq[anchor0 + 2:]
        full = GUIDE + "NRG"
        als = ([(len(PRE), a) for a in C.align_pairs(full, hap[len(PRE):len(PRE) + 24],
                                                     4, 1, 1, 6, 3, False)]
               + [(len(PRE) + 1, a) for a in C.align_pairs(
                   full, hap[len(PRE) + 1:len(PRE) + 24], 4, 1, 1, 6, 3, False)])
        cfd = C.Cfd(C.default_cfd_dir())

        def few_key(x):
            _s, (_sa, _ta, mm, dn, rn) = x
            return (mm + dn + rn, mm)

        def cfd_key(x):
            _s, (sa, ta, mm, dn, rn) = x
            return (-C.cfd3(cfd.score(sa, ta, 3)), mm + dn + rn, mm)

        def block(x):
            start, (sa, ta, mm, dn, rn) = x
            ref = list(seq[start:start + len(ta.replace("-", ""))])
            for i, ch in enumerate(ta):
                if ch == "-":
                    ref.insert(i, "-")
            bt = ",".join(t for t, c in (("RNA", rn), ("DNA", dn)) if c) or "X"
            return row_fields("chr21", start, "+", sa.replace("NRG", "NNN"), "".join(ref),
                              mark(sa, ta), ta[-3:], mm, dn + rn, bt, var, "S01,S02")
        right_few = min(als, key=few_key)
        right_cfd = min(als, key=cfd_key)
        tied = [x for x in als if few_key(x)[0] == few_key(right_few)[0]
                and few_key(x) > few_key(right_few)]
        return seq, rec, block, right_few, right_cfd, tied

    def test_a_fewest_block_worse_by_the_next_sorting_criterion_is_caught(self):
        """R5-F2 (adv_r5b): the fewest-mm+b block holds the 23-mer, tied with the 24-mer
        on mm+bulges and worse on mm. The merge's next key decides, so the 24-mer is lost
        (and the 22-mer the site's RNA-bulge alignment holds, which folds there too). With
        --sorting-criteria mm+bulges alone the tie is the merge's to keep."""
        seq, rec, block, right_few, right_cfd, tied = self._two_targets_one_site()
        self.assertEqual((right_few[0], right_few[1][2:]), (len(PRE), (0, 1, 0)))
        self.assertEqual(tied[0][1][2:], (1, 0, 0))
        case = Case(self.tmp, seq, [rec], [blocks_row([block(right_cfd), block(right_few),
                                                       block(right_cfd)])])
        rc, found = case.run(mode="both")
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, [blocks_row([block(right_cfd), block(tied[0]),
                                             block(right_cfd)])])
        rc, found = case.run(mode="completeness")
        hits = [f for f in gating(found) if "this target at %d " % len(PRE) in f[2]]
        self.assertEqual(sorted((f[1], f[3]) for f in hits),
                         [("missing-best-alignment", "S01"),
                          ("missing-best-alignment", "S02")], found)
        self.assertIn("fewest-mm+b block of 1 mm+bulges, 1 mm, the target 1 mm+bulges, "
                      "0 mm (worse by the merge's next sorting criterion", hits[0][2])
        with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
            fh.write("crisprme.py complete-search --sorting-criteria mm+bulges\n")
        rc, found = case.run(mode="completeness")
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_cfd_tie_is_broken_by_the_sorting_criteria(self):
        """R5-F2: the CFD file sorts on the CFD at three decimals, then on the sorting
        criteria: a fold group whose highest-CFD block ties the target's best CFD at
        three decimals with more edits than the target's alignment of that CFD has kept
        the worse; one of a higher CFD, or within a rounding step of a boundary, has
        not."""
        reported = C.Reported([], argparse_ns_ctx(), criterion="mm+bulges,mm")

        def entry(block, mm, bul, cfd):
            return C.Entry({"block": block, "mm": str(mm), "bulges": str(bul),
                            "cfd": "%.6f" % cfd, "crrna": "", "alt": ""}, frozenset(), 0)
        d = C.Best(C.edits(2, 0), 0.5001)
        d.add(C.edits(4, 0), 0.4)
        few = entry("fewest_mm+b", 2, 0, 0.5001)
        self.assertIsNone(reported.worse([few, entry("highest_CFD", 2, 0, 0.5001)], d))
        why = reported.worse([few, entry("highest_CFD", 2, 1, 0.5003)], d)
        self.assertIn("a tie at three decimals, which the merge breaks on "
                      "--sorting-criteria mm+bulges,mm", why)
        # a CRISPRme whose CFD merge reads --sorting-criteria-scoring (here mm: the row's 2
        # mismatches tie the target's) breaks the tie on that: either reading is accepted
        lenient = C.Reported([], argparse_ns_ctx(), criterion="mm+bulges,mm", scoring="mm")
        self.assertIsNone(lenient.worse([few, entry("highest_CFD", 2, 1, 0.5003)], d))
        self.assertIsNotNone(lenient.worse([few, entry("highest_CFD", 3, 1, 0.5003)], d))
        # a second CFD row of the group, of a lower CFD, changes nothing either way
        self.assertIsNotNone(reported.worse([few, entry("highest_CFD", 2, 1, 0.5003),
                                             entry("highest_CFD", 0, 0, 0.3)], d))
        # the run's own command line names the criteria the merge read
        res = os.path.join(self.tmp, "run")
        os.makedirs(res)
        with open(os.path.join(res, ".command_line.txt"), "w") as fh:
            fh.write("crisprme.py complete-search --sorting-criteria mm,bulges "
                     "--sorting-criteria-scoring mm+bulges --max-total-edits 5\n")
        self.assertEqual(C.run_kind(res)[3:6], ("mm,bulges", "mm+bulges", 5))
        self.assertIsNone(reported.worse([few, entry("highest_CFD", 3, 1, 0.501)], d))
        d2 = C.Best(C.edits(2, 0), 0.5005)   # prints 0.500 or 0.501: not surely a tie
        self.assertIsNone(reported.worse([few, entry("highest_CFD", 3, 1, 0.5005)], d2))
        # a summary of several alignments gives no key: only real alignments do
        d3 = C.Demand("exact", [], C.edits(0, 0), 0.9, 0)
        self.assertEqual(d3.keys, set())

    # --- R5-F4: a fold across strands loses a target ---------------------------- #
    def test_a_fold_across_strands_is_a_lost_target(self):
        """R5-F4: a deletion makes one 23-mer that reads as a target on both strands at
        one Cluster Position (a perfect + target, a 1-mismatch - target) with the same
        variants and carriers. CRISPRme's fold key holds no strand, so a merge that
        folds them keeps one row for two targets: the other strand's is reported lost,
        not accepted as an alignment of the kept one."""
        x = "CCAGTACGATCAGATCGTACTGG"     # x and rc(x) differ at their 12th base only
        g = x[:20]
        seq = PRE + x[:15] + "G" + x[15:] + POST
        s = len(PRE)
        rec = (s + 15, x[14] + "G", x[14], {"S01": "1|0", "S02": "0|1"})
        var = tok("chr21", rec)
        ref = seq[s:s + 23]
        plus = row_fields("chr21", s, "+", g + "NNN", ref, x, "TGG", 0, 0, "X", var,
                          "S01,S02", spacer=g + "NNN")
        minus = row_fields("chr21", s, "-", g + "NNN", C.rc(ref), mark(g + "NNN", C.rc(x)),
                           "TGG", 1, 0, "X", var, "S01,S02", spacer=g + "NNN")
        case = Case(self.tmp, seq, [rec], [plus, minus])
        budget = ("1", "0", "0", "1")
        rc, found = case.run(mode="both", budget=budget, guide=g)
        self.assertEqual((rc, found), (0, []), found)
        write_results(case.res, [plus])
        rc, _found = case.run(mode="completeness", budget=budget, guide=g)
        lines = read_tsv(os.path.join(self.tmp, "fail.tsv"))
        self.assertEqual(rc, 1, lines)
        self.assertEqual(sorted((r["class"], r["strand"], r["samples"]) for r in lines),
                         [("missing-row", "-", "S01"), ("missing-row", "-", "S02")], lines)

    # --- R5-F5: the read windows ---------------------------------------------------- #
    def test_the_cluster_chain_reads_the_records_it_reaches(self):
        """R5-F5.1: a population-level cluster chains through a gap only the SNVs of the
        enriched genome bridge; a window that stops before them (a tile's) is read again
        as far as the chain can reach."""
        seq = "A" * 50 + "G" * 400 + "ATAT" + "G" * 400 + "A" * 50
        tmp = self.tmp
        vcf = os.path.join(tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", [(p + 1, seq[p], "G", {"S01": "0|1"})
                                 for p in range(450, 454)])
        import argparse
        ns = argparse.Namespace(mm=0, bdna=0, brna=0, max_total_edits=0)
        ctx = C.Ctx(ns, "NGG", 3, False, ["G" * 20], None)
        win = C.Window(seq, "chr21", 0, 440, [C.Vcf(vcf)])
        self.assertEqual(win.snv_at, {})
        members, cut, _why = C.cluster_chain(win, ctx, "G" * 20, 300, 3)
        self.assertGreater(max(members), 454 + 300, members[-5:])
        self.assertFalse(cut)

    def test_a_long_ref_an_unindexed_read_meets_is_listed(self):
        """R5-F5.2: without a .tbi a read starts Vcf.LONGEST before its region, so a
        longer deletion that starts further before is not read there: completeness lists
        each VCF whose reads met such a record (its own coverage gap, non-gating); with a
        .tbi there is nothing to list."""
        import random
        rnd = random.Random(7)
        dele = "".join(rnd.choice("ACGT") for _ in range(1500))
        seq = PRE + GUIDE[:10] + dele + GUIDE[10:] + "TGG" + POST + LONG
        a = len(PRE) + 9
        drec = (a + 1, GUIDE[9] + dele, GUIDE[9], {"S01": "1|0"})
        s = len(PRE)
        row = row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], GUIDE + "TGG", "TGG",
                         0, 0, "X", tok("chr21", drec), "S01")
        case = Case(self.tmp, seq, [drec], [row])
        bg = os.path.join(self.tmp, "idx.vcf.gz")
        write_bgzf_tabix(bg, vcf_text("chr21", [drec]))
        case.vcf = bg
        old = C.Vcf.LONGEST
        try:
            C.Vcf.LONGEST = 1000
            rc, found = case.run(mode="both", budget=("0", "0", "0", "0"))
            self.assertEqual((rc, gating(found)), (0, []), found)
            self.assertNotIn("coverage-unindexed-long-ref", classes(found))
            os.remove(bg + ".tbi")
            rc, found = case.run(mode="both", budget=("0", "0", "0", "0"))
            self.assertEqual((rc, gating(found)), (0, []), found)
            lines = [f for f in found if f[1] == "coverage-unindexed-long-ref"]
            self.assertEqual(len(lines), 1, found)
            self.assertIn("holds a record whose REF spans 1501 bp", lines[0][2])
        finally:
            C.Vcf.LONGEST = old


def world_push(tmp):
    """R6-F1 (A), the shape of the review's five tandem-repeat sites (chr22:17799944 and
    others): an (AGG)n run. R inserts two AGG units after the run's one reference unit
    (S02's slot 1; S01's slot 0); D inserts one unit before the run, on S01's slot 0 only.
    Read from the target's first base, S01's slot has the target's bases, but D lengthened
    the stretch and pushed R's second unit out of it."""
    g = "ACTGATTCTGTACC" + "AGGAGG"
    s0 = len(LONG) + len(PRE)
    seq = LONG + PRE + "ACTGATTCTGTACC" + "AGG" + "TTTCATTTCA" + POST
    d = (s0 + 14, "C", "CAGG", {"S01": "1|0"})
    r = (s0 + 17, "G", "GAGGAGG", {"S01": "1|0", "S02": "0|1"})

    def row(who):
        return row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], g + "AGG", "AGG",
                          0, 0, "X", tok("chr21", r), who, spacer=g + "NNN")
    return seq, [d, r], row, g, s0


def world_snv_bases(tmp):
    """R6-F1 (B), the shape of chr22:32030641-45: R deletes the T of P+GT|GGA; S (A>G) on
    S02's slot with R makes the target P+G|GG(G). S03's slot holds R and D, a G inserted
    after the run's last G, which puts the same bases there without S."""
    g = "TTGACCGGTATCAGCTAGCG"   # P (19 bases) + R's anchor G
    s0 = len(LONG) + len(PRE)
    seq = LONG + PRE + g[:19] + "G" + "T" + "GG" + "A" + "CATTTCAGTTAC" + POST
    a = s0 + 19                  # R's anchor, 0-based
    r = (a + 1, "GT", "G", {"S02": "1|0", "S03": "1|0"})
    snv = (a + 5, "A", "G", {"S02": "1|0"})
    d = (a + 4, "G", "GG", {"S03": "1|0"})
    var = "%s,%s" % (tok("chr21", r), tok("chr21", snv))
    row = row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], g + "GGG", "GGG", 0, 0,
                     "X", var, "S02", spacer=g + "NNN")
    return seq, [r, snv, d], row, g


def world_bridge(rows_kind):
    """R6-F2, the shape of the review's population-level line 646: an (AGG)n run whose
    reference lacks one base of the run (ideal offset 4), which R puts back (S01). Its
    target T1 (ideal 0, one mismatch) is at Cluster Position s; X (ideal 3) needs R and
    the SNV S (S02's, on no slot of R) and is at s+3; B, a SNP row of S (ideal 6), at s+5:
    T1 and B are two clusters, which only X, a target no carrier holds, would join."""
    g = "AGGAGGAGGAGGAGGAGGAG"
    ideal = list("AGG" * 16)
    ideal[5] = "T"      # a mismatch of T1 and X
    ideal[24] = "T"     # S (T>A) fixes it: a mismatch of X and B
    s0 = len(LONG) + len(PRE)
    seq = LONG + PRE + "".join(ideal[:4] + ideal[5:]) + "TTTCATTTCA" + POST
    r = (s0 + 4, seq[s0 + 3], seq[s0 + 3] + "G", {"S01": "1|0"})
    sp = s0 + 23
    snv = (sp + 1, "T", "A", {"S02": "0|1"})
    rows = []
    if "T1" in rows_kind:
        alt = "".join(ideal[0:23])
        rows.append(row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23],
                               mark(g + "NNN", alt), alt[-3:], 1, 0, "X", tok("chr21", r),
                               "S01", spacer=g + "NNN"))
    if "B" in rows_kind:
        b0 = s0 + 5
        alt = seq[b0:sp] + "A" + seq[sp + 1:b0 + 23]
        rows.append(row_fields("chr21", b0, "+", g + "NNN", seq[b0:b0 + 23], alt, alt[-3:],
                               0, 0, "X", tok("chr21", snv), "S02", spacer=g + "NNN"))
    return seq, [r, snv], rows, g


class ReviewRound6(Base):
    def _ctx(self, guide, mm=0):
        import argparse
        ns = argparse.Namespace(mm=mm, bdna=0, brna=0, max_total_edits=mm)
        return C.Ctx(ns, "NRG", 3, False, [guide], None)

    def _row_target(self, case, guide):
        win = C.Window(C.Reference(os.path.join(case.tmp, "Genomes")).seq("chr21"),
                       "chr21", 0, 10 ** 6, [C.Vcf(case.vcf)])
        row = C.read_rows(results_table(case), "integrated")[0]
        named = C.named_alleles(row, win, C.Report(), True)
        main = [x for x in named if x.kind != "snv"][0]
        tgt = C.row_target(win, row, named, main)[0]
        return win, tgt, C.local_lines(win, tgt)

    # --- R6-F1 (A): a second insertion that pushes the record's bases out --------- #
    def test_an_upstream_insertion_that_pushes_the_records_bases_out_lists_no_one(self):
        """R6-F1 (A): CRISPRme reads S01's slot through R (METHODS section 4), and D
        pushed R's second unit out of the stretch read from the target's first base (up to
        its last base the bases differ): S01 is not listed for the row, and its own target
        needs two indel-path records. 4c271fb demanded it (carrier-missing, missing-row);
        a row that lists S01 names R alone for bases D makes there (variant-missing)."""
        seq, recs, row, g, _s0 = world_push(self.tmp)
        b0 = ("0", "0", "0", "0")
        case = Case(self.tmp, seq, recs, [row("S02")])
        rc, found = case.run(mode="both", guide=g, budget=b0)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "limit-two-indels" and f[3] == "S01" for f in found),
                        found)
        win, tgt, lids = self._row_target(case, g)
        self.assertEqual(C.documented(win, self._ctx(g), tgt, "S01", lids)[0], "off")
        self.assertEqual(C.truth(win, tgt, "S01", lids)[:2], ("exact", False))
        self.assertEqual(C.truth(win, tgt, "S02", lids)[:2], ("exact", True))
        rc, found = Case(self.tmp, seq, recs, [row("S01,S02")]).run(guide=g, budget=b0)
        self.assertEqual([(f[1], f[3]) for f in gating(found)][:1],
                         [("variant-missing", "S01")], found)
        rc, found = Case(self.tmp, seq, recs, [row("S01")]).run(mode="both", guide=g,
                                                                 budget=b0)
        self.assertIn(("carrier-missing", "S02"), [(f[1], f[3]) for f in gating(found)])
        self.assertIn(("missing-row", "S02"), [(f[1], f[3]) for f in gating(found)])

    def test_the_core_of_a_reading_is_the_records_bases_and_the_base_past_them(self):
        """record_core: an insertion's anchor and added bases the target holds and, where
        the target goes on, the base after them; a deletion's anchor and junction."""
        ins = al(10, "A", "ACGT")   # anchor 9, adds C G T at (9, 1..3)
        place = [(7, 0), (8, 0), (9, 0), (9, 1), (9, 2), (9, 3), (10, 0), (11, 0)]
        src = [None, None, None, ins, ins, ins, None, None]
        self.assertEqual(C.record_core(place, src, ins),
                         {(9, 0), (9, 1), (9, 2), (9, 3), (10, 0)})
        self.assertEqual(C.record_core(place[:5], src[:5], ins), {(9, 0), (9, 1), (9, 2)})
        dele = al(10, "ACG", "A")   # anchor 9, deletes 10 and 11
        self.assertEqual(
            C.record_core([(8, 0), (9, 0), (12, 0), (13, 0)], [None] * 4, dele),
            {(9, 0), (12, 0)})

    # --- R6-F1 (B): a carrier's own SNP alleles make the window's target ---------- #
    def test_a_second_indel_giving_an_snvs_bases_lists_no_carrier_without_the_snv(self):
        """R6-F1 (B): METHODS section 4 resolves the window's target with each carrier's
        own SNP alleles and names those it holds: S03 holds no S, so its window target is
        another sequence, whatever D (a G inserted in the G run) puts at S's place. 4c271fb
        read S03's bases alone and listed it, so its fold carrier sets demanded S03 of
        every row of the target; S03 is not a carrier CRISPRme lists, nor required."""
        seq, recs, row, g = world_snv_bases(self.tmp)
        case = Case(self.tmp, seq, recs, [row])
        rc, found = case.run(mode="both", guide=g, budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        win, tgt, lids = self._row_target(case, g)
        ctx = self._ctx(g)
        self.assertEqual(C.documented(win, ctx, tgt, "S03", lids)[0], "off")
        self.assertEqual(C.documented(win, ctx, tgt, "S02", lids)[0], "lists")
        self.assertEqual(C.documented_carriers(win, ctx, tgt, {"S02", "S03"}, lids),
                         (frozenset({"S02"}), frozenset()))

    def test_a_carrier_listed_against_its_haplotype_vetoes_no_fold(self):
        """R6-F1, the fold: a carrier the documented rules list for a target its own
        haplotypes do not hold at the row's start (S01, read up to the target's last base
        after it lost a CA unit at the target's left edge) may be on the target's row, and
        is not required there; S02, whose haplotype and rules agree, is."""
        g = "CACATTGTCAGTACGATCGA"
        seq = LONG + "GTCA" + g[:10] + "C" + g[10:] + "TGG" + POST
        s0 = len(LONG) + 4
        d1 = (s0 + 10, g[9] + "C", g[9], {"S01": "1|0", "S02": "1|0"})
        d2 = (s0, seq[s0 - 1:s0 + 2], seq[s0 - 1], {"S01": "1|0"})
        rows = [row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], g + "TGG", "TGG",
                           0, 0, "X", tok("chr21", d1), "S01,S02", spacer=g + "NNN")]
        case = Case(self.tmp, seq, [d1, d2], rows)
        win, tgt, lids = self._row_target(case, g)
        ctx = self._ctx(g)
        self.assertEqual(C.documented(win, ctx, tgt, "S01", lids)[0], "lists")
        self.assertEqual(C.truth_lists(win, ctx, tgt, "S01", lids), "no")
        self.assertEqual(C.documented_carriers(win, ctx, tgt, {"S01", "S02"}, lids),
                         (frozenset({"S02"}), frozenset({"S01"})))

    # --- R6-F2: the cluster chain ------------------------------------------------- #
    def test_a_cluster_chains_through_no_target_that_no_carrier_holds(self):
        """R6-F2: one row per cluster. T1 and B are two clusters of CRISPRme's merge: X,
        between them, needs R and S together, which no haplotype holds, so the per-carrier
        resolution gives it no row. A table that lost T1's row has lost a cluster; 4c271fb
        chained T1 to B through X (every SNV allele an option on R's window) and passed."""
        for kind, want in ((("T1", "B"), 0), (("B",), 1)):
            seq, recs, rows, g = world_bridge(kind)
            case = Case(self.tmp, seq, recs, rows)
            run_markers(case.res, "population-level")
            with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
                fh.write("input_command\tcrisprme.py complete-search --population-level\n")
            rc, found = case.run(mode="completeness", guide=g, budget=("1", "0", "0", "1"))
            self.assertEqual(rc, want, found)
            if want:
                self.assertEqual([(f[1], f[3]) for f in gating(found)],
                                 [("missing-locus", "S01")], found)

    def test_the_snp_rows_a_cluster_holds_follow_the_analysis_mode(self):
        """R6-F2, the same hole on the reference: a population-level run's SNP rows are
        the worst case of each alignment over every allele (each dataset apart), a
        per-sample run's (one made with --no-alt-alignments keeps one row per cluster too)
        only the haplotypes individuals hold. Two SNVs on different samples: one union in population-level,
        each sample's own in per-sample, never both together."""
        seq = "ACGT" * 20
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", [(5, "A", "C", {"S01": "1|0"}),
                                 (9, "A", "T", {"S02": "0/1"}),
                                 (13, "A", "G", {"S01": "1|1"})])
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
        ctx = self._ctx("A" * 20)
        pop = C.snp_row_sets(win, ctx, 0, 20)
        self.assertIn(frozenset("AC"), {x[4] for x in pop})
        self.assertTrue(any(x[4] == frozenset("AC") and x[8] == frozenset("AT")
                            for x in pop))
        ctx.snp_union = False
        per = C.snp_row_sets(win, ctx, 0, 20)
        got = {(x[4], x[8], x[12]) for x in per}
        a, c, t, gg = (frozenset(b) for b in "ACTG")
        self.assertEqual(got, {(a, a, a), (c, a, gg), (a, a, gg), (a, frozenset("AT"), a)})

    # --- R6-F3: two genotyped datasets that share samples ------------------------- #
    def test_two_genotyped_vcfs_that_share_samples_are_refused(self):
        """R6-F3: the same records and calls in two genotyped VCFs (two datasets of one
        panel's individuals, as HPRC beside HGDP + 1000 Genomes): each dataset is resolved
        on its own genotypes and the merge joins their rows, which this checker does not
        model, so it refuses the run (exit 2), in either mode; datasets with no sample in
        common are checked (a dropped carrier is still caught)."""
        case, info = world_indel_snp(self.tmp)
        v2 = os.path.join(self.tmp, "v2.vcf.gz")
        write_vcf(v2, "chr21", info["recs"])
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, extra=["--vcf", v2])
            self.assertEqual((rc, found), (2, []), (mode, found))
        v3 = os.path.join(self.tmp, "v3.vcf.gz")
        write_vcf(v3, "chr21", info["recs"], samples=["X%02d" % i for i in range(1, 9)])
        rc, found = case.run(extra=["--vcf", v3])
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, [set_field(info["rows"][0], samples="S06")])
        rc, found = case.run(extra=["--vcf", v3])
        self.assertIn("carrier-missing", classes(gating(found)), found)

    # --- R6-F4: no known-slot cap; the record's REF span is no choice ------------- #
    def test_calls_of_unknown_phase_on_the_records_ref_span_cap_nothing(self):
        """R6-F4: S01 (the deletion phased) and S02 (unphased) have 11 unphased SNVs on
        the bases the deletion removes: an allele there makes no sequence on the
        deletion's slot and changes nothing CRISPRme reads (METHODS section 4). 4c271fb
        capped both (the known slot's 11 calls, and 11 PUTATIVE choices) and a table
        without the row passed; both are demanded now."""
        junk = "CATCGATTCAGC"
        seq = PRE + GUIDE[:10] + junk + GUIDE[10:] + "TGG" + POST
        rec = (len(PRE) + 10, GUIDE[9] + junk, GUIDE[9], {"S01": "1|0", "S02": "0/1"})
        recs = [rec] + [(len(PRE) + 11 + k, junk[k], other(junk[k]),
                         {"S01": "0/1", "S02": "0/1"}) for k in range(1, 12)]
        s = len(PRE)
        rows = [row_fields("chr21", s, "+", GUIDE + "NNN", seq[s:s + 23], GUIDE + "TGG",
                           "TGG", 0, 0, "X", tok("chr21", rec), "S01,S02")]
        rc, found = Case(self.tmp, seq, recs, rows).run(mode="both",
                                                         budget=("0", "0", "0", "0"))
        self.assertEqual((rc, found), (0, []), found)
        rc, found = Case(self.tmp, seq, recs, []).run(mode="completeness",
                                                       budget=("0", "0", "0", "0"))
        self.assertEqual(sorted((f[1], f[3]) for f in gating(found)),
                         [("missing-putative-row", "S02"), ("missing-row", "S01")], found)
        self.assertNotIn("coverage-putative-cap", classes(found))

    def test_a_cluster_past_the_chains_reach_is_a_coverage_gap_not_a_miss(self):
        """The sweep: a population-level cluster in a 7 kb (AGG)n run is one merge
        cluster, whose one row can lie past the 3 kb either side the chain follows: a
        deletion of one unit in the middle has its row 3.6 kb away. The checker cannot
        tell then (coverage-chain-limit, non-gating), where it reported missing-locus."""
        g = "AGGAGGAGGAGGAGGAGGAG"
        run = "AGG" * 2400
        seq = LONG + "TTTC" + run + "TTTC" + POST
        s0 = len(LONG) + 4
        mid = s0 + 3 * 1200                 # a unit boundary in the middle of the run
        rec = (mid, seq[mid - 1:mid + 3], seq[mid - 1], {"S01": "1|0"})
        row = set_field(row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23],
                                   seq[s0:s0 + 23], seq[s0 + 20:s0 + 23], 0, 0, "X", "NA",
                                   "NA", spacer=g + "NNN"), origin="ref")
        case = Case(self.tmp, seq, [rec], [row])
        run_markers(case.res, "population-level")
        rc, found = case.run(mode="completeness", guide=g, budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("coverage-chain-limit", classes(found), found)

    # --- R6-F5: a row that lists only limit carriers ------------------------------ #
    def test_a_limit_line_says_how_many_of_its_carriers_a_row_lists(self):
        """R6-F5: S01's slot holds two other alleles at the insertion's own position (a
        phasing artefact): the VCF gives no sequence, so its target is demanded of no one
        and its limit line stays. The line says how many of its carriers a row naming the
        record lists, so a table that lost the row they were on gives other findings
        (non-gating: CRISPRme's reading of such a slot is not modeled)."""
        seq = PRE + GUIDE[:10] + GUIDE[13:] + "TGG" + POST
        pos1 = len(PRE) + 10
        recs = [(pos1, GUIDE[9], GUIDE[9] + GUIDE[10:13], {"S01": "1|0", "S02": "0|1"}),
                (pos1, GUIDE[9], GUIDE[9] + "AAA", {"S01": "1|0"}),
                (pos1, GUIDE[9], other(GUIDE[9]), {"S01": "1|0"})]
        start = len(PRE)
        lines = {}
        for who in ("S01,S02", "S02"):
            rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                               GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[0]), who)]
            rc, found = Case(self.tmp, seq, recs, rows).run(mode="completeness")
            self.assertEqual((rc, gating(found)), (0, []), found)
            lines[who] = [f[2] for f in found if f[1] == "limit-overlapping-records"
                          and "hold %s on" % tok("chr21", recs[0]) in f[2]]
        self.assertEqual(len(lines["S02"]), 1, lines)
        self.assertIn("1 of them listed on a row naming the record", lines["S01,S02"][0])
        self.assertIn("0 of them listed on a row naming the record", lines["S02"][0])

    def test_a_target_only_an_unknown_call_makes_is_still_demanded(self):
        """R6-F4: with no placement of a known slot's calls of unknown phase enumerated,
        a target that exists only with one of them (an unphased SNV that makes the PAM)
        is still the PUTATIVE rule's: S01's deletion is phased, the SNV unphased."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TCG" + POST
        rec = (len(PRE) + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"})
        p = len(PRE) + 10 + 1 + 10 + 1     # the PAM's middle base, C>G
        snv = (p + 1, "C", "G", {"S01": "0/1"})
        rc, found = Case(self.tmp, seq, [rec, snv], []).run(mode="completeness",
                                                            budget=("0", "0", "0", "0"))
        self.assertEqual([(f[1], f[3]) for f in gating(found)],
                         [("missing-putative-row", "S01")], found)


def argparse_ns_ctx():
    """A Ctx for Reported (a 3' NGG PAM, no CFD tables: rows keep their printed CFD)."""
    import argparse
    ns = argparse.Namespace(mm=4, bdna=1, brna=1, max_total_edits=6)
    return C.Ctx(ns, "NRG", 3, False, [GUIDE], None)


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]


def open_findings(_rc, tmp):
    with open(os.path.join(tmp, "o.tsv")) as fh:
        return fh.read()


def argparse_ns():
    import argparse
    return argparse.Namespace(mm=4, bdna=1, brna=1, max_total_edits=6)


# --------------------------------------------------------------------------- #
# the driver never passes silently
# --------------------------------------------------------------------------- #
def world_phantom(art=None, s01="1|0", s_on=None):
    """R7-F1: the deletion world, with S at spacer position 3 (S02's, 1|1, unless
    ``s_on`` puts it on S01 too) and ``art``, an allele overlapping the deletion on S01:
    'pos' (a SNV at its own position: a phasing artefact), 'deleted' (a SNV on the base it
    deletes), 'complex' (an insertion at the next position). Returns (seq, records, the
    row naming the deletion and S with S's base written in, the row naming the deletion
    alone, the deletion)."""
    seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
    pos1 = len(PRE) + 10
    s_pos0 = len(PRE) + 3
    s_alt = other(seq[s_pos0])
    rec = (pos1, GUIDE[9] + "C", GUIDE[9], {"S01": s01, "S03": "0|1"})
    snp_gts = {"S02": "1|1"}
    if s_on:
        snp_gts["S01"] = s_on
    snp = (s_pos0 + 1, seq[s_pos0], s_alt, snp_gts)
    recs = [rec, snp]
    if art == "pos":
        recs.append((pos1, GUIDE[9], other(GUIDE[9]), {"S01": "1|0"}))
    elif art == "deleted":
        recs.append((pos1 + 1, "C", "A", {"S01": "1|0"}))
    elif art == "complex":
        recs.append((pos1 + 1, "C", "CT", {"S01": "1|0"}))
    start = len(PRE)
    with_s = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                        GUIDE[:3] + s_alt.lower() + GUIDE[4:] + "TGG", "TGG", 1, 0, "X",
                        "%s,%s" % (tok("chr21", rec), tok("chr21", snp)), "S01")
    alone = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                       GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01,S03")
    return seq, recs, with_s, alone, rec


def world_bridge_overlap(rows_kind, r2=None):
    """R7-F2 (the review's adv_r7b): world_bridge's shape, with X held by S01's slot of
    R and S. ``r2`` adds an allele on that slot overlapping R: 'complex' (a deletion of
    R's anchor base: METHODS takes the slot off, so X has no row and T1 (S03's) and B are
    two clusters), 'pos' (a SNV at R's own position: a phasing artefact, read without it,
    so X keeps its row)."""
    g = "AGGAGGAGGAGGAGGAGGAG"
    ideal = list("AGG" * 16)
    ideal[5] = "T"
    ideal[24] = "T"
    s0 = len(LONG) + len(PRE)
    seq = LONG + PRE + "".join(ideal[:4] + ideal[5:]) + "TTTCATTTCA" + POST
    r = (s0 + 4, seq[s0 + 3], seq[s0 + 3] + "G", {"S01": "1|0", "S03": "1|0"})
    sp = s0 + 23
    snv = (sp + 1, "T", "A", {"S02": "0|1", "S01": "1|0"})
    recs = [r, snv]
    if r2 == "complex":
        recs.append((s0 + 3, seq[s0 + 2:s0 + 4], seq[s0 + 2], {"S01": "1|0"}))
    elif r2 == "pos":
        recs.append((s0 + 4, seq[s0 + 3], other(seq[s0 + 3]), {"S01": "1|0"}))
    rows = []
    if "T1" in rows_kind:
        alt = "".join(ideal[0:23])
        rows.append(row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23],
                               mark(g + "NNN", alt), alt[-3:], 1, 0, "X", tok("chr21", r),
                               "S03", spacer=g + "NNN"))
    if "B" in rows_kind:
        b0 = s0 + 5
        alt = seq[b0:sp] + "A" + seq[sp + 1:b0 + 23]
        rows.append(row_fields("chr21", b0, "+", g + "NNN", seq[b0:b0 + 23], alt, alt[-3:],
                               0, 0, "X", tok("chr21", snv), "S01,S02", spacer=g + "NNN"))
    return seq, recs, rows, g, s0


class ReviewRound7(Base):
    def _ctx(self, guide, mm=0, cfd=None):
        import argparse
        ns = argparse.Namespace(mm=mm, bdna=0, brna=0, max_total_edits=mm)
        return C.Ctx(ns, "NRG", 3, False, [guide], cfd)

    def _popmode(self, case):
        run_markers(case.res, "population-level")
        with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
            fh.write("input_command\tcrisprme.py complete-search --population-level\n")

    # --- R7-F1: a carrier whose slot's alleles overlap the record ------------------ #
    def test_a_phantom_snv_on_a_carrier_whose_slot_overlaps_is_caught(self):
        """R7-F1 (the review's W1, W1b, W1c): the row names R and S and lists S01 alone,
        who holds S on no slot. S01's slot of R also holds an allele overlapping R (a SNV
        at R's position, a SNV on the base R deletes; R homozygous with that allele on one
        slot): 'no sequence' excused it (limit-overlapping-records); a slot that holds a
        named variant on no reading lacks the target whatever its overlaps."""
        for art, s01 in (("pos", "1|0"), ("deleted", "1|0"), ("pos", "1|1"),
                         ("complex", "1|0")):
            seq, recs, with_s, _alone, _rec = world_phantom(art, s01)
            rc, found = Case(self.tmp, seq, recs, [with_s]).run()
            self.assertEqual(rc, 1, (art, s01, found))
            self.assertTrue({(f[1], f[3]) for f in gating(found)} <= {
                ("carrier-without-variants", "S01"), ("carrier-off-by-rule", "S01")},
                (art, s01, found))
        # the control: S on S01's slot of R, the same row is right (a note: the slot, read
        # as METHODS section 4 reads it, holds the target; R8-N1)
        seq, recs, with_s, _alone, _rec = world_phantom("pos", "1|0", s_on="1|0")
        rc, found = Case(self.tmp, seq, recs, [with_s]).run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual({(f[1], f[3]) for f in found},
                         {("overlap-reading-agrees", "S01")}, found)

    def test_a_slot_read_without_its_artefact_must_hold_the_target(self):
        """R7-F1: METHODS section 4 reads a slot with a phasing artefact (another allele
        at the record's position) without it, and S01 stays listed for the target of
        that reading. S01's slot holds S, inside the target: listed for the target without
        S (the row naming R alone) it holds a variant the row does not name."""
        seq, recs, with_s, alone, _rec = world_phantom("pos", "1|0", s_on="1|0")
        rc, found = Case(self.tmp, seq, recs, [alone]).run()
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("variant-missing", "S01")}, found)
        # a variant only on a base the record deletes changes nothing: the same reading
        seq, recs, with_s, alone, _rec = world_phantom("deleted", "1|0", s_on="1|0")
        rc, found = Case(self.tmp, seq, recs, [alone]).run()
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("variant-missing", "S01")}, found)
        rc, found = Case(self.tmp, seq, recs, [with_s]).run()
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_carrier_off_by_rule_vetoes_its_fold(self):
        """R7-F1's sweep: S01's slot of the deletion holds an insertion overlapping it at
        another position (a complex allele: METHODS takes the slot off), so a fold group
        that lists S01 is not the target's own; S03, whose slot is clean, is required. A
        slot read without a phasing artefact that holds the target is the truth's, so it
        may be listed (it is never required)."""
        seq, recs, _with_s, alone, rec = world_phantom("complex", "1|0")
        case = Case(self.tmp, seq, recs, [alone])
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
        row = C.read_rows(results_table(case), "integrated")[0]
        named = C.named_alleles(row, win, C.Report(), True)
        tgt = C.row_target(win, row, named, named[0])[0]
        lids = C.local_lines(win, tgt)
        ctx = self._ctx(GUIDE)
        self.assertEqual(C.truth_lists(win, ctx, tgt, "S01", lids), "clash")
        self.assertEqual(C.documented(win, ctx, tgt, "S01", lids)[0], "off")
        self.assertEqual(C.documented_carriers(win, ctx, tgt, {"S01", "S03"}, lids),
                         (frozenset({"S03"}), frozenset()))
        # a slot read without a phasing artefact that holds the target: the truth lists
        # it (never required: not a gating demand), whatever the documented view says
        seq, recs, _with_s, alone, rec = world_phantom("pos", "1|0")
        case = Case(self.tmp, seq, recs, [alone])
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
        row = C.read_rows(results_table(case), "integrated")[0]
        named = C.named_alleles(row, win, C.Report(), True)
        tgt = C.row_target(win, row, named, named[0])[0]
        lids = C.local_lines(win, tgt)
        self.assertEqual(C.truth_lists(win, ctx, tgt, "S01", lids), "resolved")
        self.assertEqual(C.documented_carriers(win, ctx, tgt, {"S01", "S03"}, lids),
                         (frozenset({"S03"}), frozenset({"S01"})))

    def test_an_overlap_explains_no_disagreement_outside_its_own_slot(self):
        """R7-F1's sweep: a slot's overlap is judged where it is (the clash), so the causes
        it gives never turn another disagreement into a limit."""
        self.assertEqual(C.unexplained_or({"artefact"}), (True, "unexplained-disagreement"))
        self.assertEqual(C.unexplained_or({"complex", "deleted"}),
                         (True, "unexplained-disagreement"))
        self.assertEqual(C.unexplained_or({"artefact", "putative"})[:2],
                         (False, "limit-putative-carrier"))

    def test_a_slot_with_no_reading_and_an_unknown_call_demands_nothing(self):
        """R7-F1's sweep: S01's slot of the deletion holds a complex allele and S01 has an
        unphased SNV in the target. The slot has no reading (METHODS takes it off), so no
        PUTATIVE target is demanded of it -- 'cannot tell' made one, whose miss the
        documented view could only call unexplained."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        q = len(PRE) + 4
        recs = [(pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"}),
                (pos1 + 1, "C", "CT", {"S01": "1|0"}),
                (q + 1, seq[q], other(seq[q]), {"S01": "0/1"})]
        rc, found = Case(self.tmp, seq, recs, []).run(mode="completeness")
        self.assertEqual((rc, gating(found)), (0, []), found)
        # the record unphased, the complex allele homozygous: it is on the record's slot,
        # whichever that is (the perfect target, which only the record makes, is no
        # PUTATIVE demand; the insertion alone makes none in this budget)
        recs = [(pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "0/1"}),
                (pos1 + 1, "C", "CT", {"S01": "1/1"})]
        rc, found = Case(self.tmp, seq, recs, []).run(mode="completeness",
                                                      budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-overlapping-records", classes(found), found)

    # --- R7-F2: the cluster chain through a slot whose alleles overlap ------------- #
    def test_a_cluster_chains_through_no_slot_a_complex_allele_takes_off(self):
        """R7-F2 (the review's adv_r7b): X is held only by S01's slot of R, which also holds
        a deletion of R's anchor base (a complex allele: METHODS takes the slot off), so
        CRISPRme writes no row for X and T1 and B are two clusters: a table that lost T1's
        row has lost a cluster. Without that deletion, or with a phasing artefact in its
        place (the slot read without it, as METHODS reads it), X joins them."""
        for r2, kind, want in (("complex", ("T1", "B"), 0), ("complex", ("B",), 1),
                               (None, ("B",), 0), ("pos", ("B",), 0)):
            seq, recs, rows, g, _s0 = world_bridge_overlap(kind, r2)
            case = Case(self.tmp, seq, recs, rows)
            self._popmode(case)
            rc, found = case.run(mode="completeness", guide=g, budget=("1", "0", "0", "1"))
            self.assertEqual(rc, want, (r2, kind, found))
            if want:
                self.assertEqual([(f[1], f[3]) for f in gating(found)],
                                 [("missing-locus", "S03")], found)

    def test_the_rows_of_an_overlap_are_written_left_open_or_none(self):
        """R7-F2: the Cluster Positions a record gives the chain from a slot whose alleles
        overlap it: those of the slot read without a phasing artefact (CRISPRme writes
        them), none for a complex allele (it takes the slot off), and those of the
        one-record reading where nothing documents the overlap (a SNV on an MNV's second
        base) as rows CRISPRme may write or not."""
        g = "AGGAGGAGGAGGAGGAGGAG"
        for r2, cps, maybe in ((None, {0, 3}, {}), ("complex", {0}, {}),
                               ("pos", {0, 3}, {})):
            seq, recs, _rows, _g, s0 = world_bridge_overlap(("B",), r2)
            case = Case(self.tmp, seq, recs, [])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
            ctx = self._ctx(g, mm=1)
            r = win.by_key[recs[0][:3]][0]
            e = C.record_demands(win, ctx, r, {"S01", "S03"}, ctx.full(g),
                                 max(0, r.lo - 80), r.hi + 80,
                                 (r.edit_lo - 26, r.edit_hi + 26))
            self.assertEqual(({c - s0 for c in e["cps"]},
                              {c - s0: w for c, w in e["cps_open"].items()}),
                             (cps, maybe), r2)
        case, info = world_mnv(self.tmp)
        mnv_pos1 = len(PRE) + 9
        recs = [(mnv_pos1, "AC", GUIDE[8:10], {"S01": "1|0"}),
                (mnv_pos1 + 1, "C", "T", {"S01": "1|0"})]
        case = Case(self.tmp, info["seq"], recs, [])
        win = C.Window(info["seq"], "chr21", 0, len(info["seq"]), [C.Vcf(case.vcf)])
        ctx = self._ctx(GUIDE)
        r = win.by_key[recs[0][:3]][0]
        e = C.record_demands(win, ctx, r, {"S01"}, ctx.full(GUIDE), 0, len(info["seq"]),
                             (r.edit_lo - 26, r.edit_hi + 26))
        self.assertEqual((e["cps"], e["cps_open"]),
                         (frozenset(), {info["start"]: frozenset({"overlap"})}), e)

    def test_a_target_through_an_snv_the_stores_read_otherwise_is_an_open_row(self):
        """R7-F2's sweep: S01's target holds the earlier record of a split SNV site, which
        CRISPRme's stores do not hold (they read the last): it is a row CRISPRme may write
        or not, so it chains a cluster only as such."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        q = len(PRE) + 4
        b = seq[q]
        x1 = other(b)
        x2 = next(c for c in "ACGT" if c not in (b, x1))
        recs = [(pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"}),
                (q + 1, b, x1, {"S01": "1|0"}), (q + 1, b, x2, {"S02": "1|0"})]
        case = Case(self.tmp, seq, recs, [])
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
        ctx = self._ctx(GUIDE, mm=1)
        r = win.by_key[recs[0][:3]][0]
        e = C.record_demands(win, ctx, r, {"S01"}, ctx.full(GUIDE), 0, len(seq),
                             (r.edit_lo - 26, r.edit_hi + 26))
        start = len(PRE)   # the Cluster Position of a + target with a 3' PAM, no bulge
        self.assertEqual(e["cps_open"].get(start), frozenset({"last-alt-wins"}), e)
        self.assertNotIn(start, e["cps"])

    # --- R7-N1: a capped window's greedy row (per-sample, classic path) ------------ #
    def test_a_capped_windows_greedy_row_is_an_open_row_of_a_per_sample_chain(self):
        """R7-N1: under --per-sample the classic dictionary path (an index built without
        --samplesID, as the self-test's) still writes a capped window's greedy worst-case
        row (METHODS section 5, known limit): eleven SNVs, each its own sample's, make the
        guide's target only together. The per-sample chain holds it as a row CRISPRme may
        write or not (a cluster that reaches its row only through it is counted); ten are
        not a capped window."""
        for n, want in ((11, True), (10, False)):
            site = list(GUIDE + "TGG")
            cols = list(range(0, 2 * n, 2))[:n] if n <= 10 else list(range(0, 20, 2)) + [19]
            for i in cols:
                site[i] = other(site[i])
            seq = PRE + "".join(site) + POST
            samples = ["T%02d" % i for i in range(len(cols))]
            vcf = os.path.join(self.tmp, "v%d.vcf.gz" % n)
            write_vcf(vcf, "chr21", [(len(PRE) + i + 1, site[i], GUIDE[i], {smp: "0|1"})
                                     for i, smp in zip(cols, samples)], samples=samples)
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
            ctx = self._ctx(GUIDE)
            ctx.snp_union = False
            rows, capped, _drop, _filt, _deg, _grd = C.snp_row_cps(win, ctx, GUIDE, ctx.full(GUIDE), 0,
                                                len(seq))
            cp = len(PRE)
            self.assertNotIn(cp, rows, n)
            self.assertEqual(cp in capped, want, n)
            ctx.snp_union = True
            rows, capped, _drop, _filt, _deg, _grd = C.snp_row_cps(
                C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)]), ctx, GUIDE,
                ctx.full(GUIDE), 0, len(seq))
            self.assertIn(cp, rows, n)   # population-level: the union is a row
            self.assertEqual(capped, set(), n)

    # --- R7-F3: a haploid call against a diploid one -------------------------------- #
    def test_a_haploid_snv_beside_a_diploid_record_leaves_the_reference(self):
        """R7-F3 (the review's adv_r7c): S01 carries the deletion diploid (1|0) and a SNV
        haploid ('1') whose alt is a mismatch. METHODS section 4 resolves a haploid call
        against a diploid one PUTATIVE, its copies counted against the record's: the
        reference stays possible, so the perfect target is S01's strongest and the row
        listing S01 for it is right."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        rec = (pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"})
        snp = (len(PRE) + 6, GUIDE[5], other(GUIDE[5]), {"S01": "1"})
        start = len(PRE)
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01")]
        case = Case(self.tmp, seq, [rec, snp], rows)
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
        q = len(PRE) + 5
        self.assertEqual(C.options(win, "S01", q), {other(GUIDE[5])})
        self.assertEqual(C.options(win, "S01", q, ploidy=2), {other(GUIDE[5]), GUIDE[5]})
        # a phantom: the row lists S01 for the SNV's alt it carries, and for the guide's
        # base too -- but not for a base it does not carry
        bad = [set_field(rows[0], alt=GUIDE[:5] + "t" + GUIDE[6:] + "TGG")]
        if other(GUIDE[5]) != "T":
            write_results(case.res, bad)
            rc, found = case.run(budget=("1", "0", "0", "1"))
            self.assertEqual(rc, 1, found)

    # --- R7-F4: the run's own record of its inputs ------------------------------------ #
    def test_a_run_given_other_vcfs_than_its_datasets_is_refused(self):
        """R7-F4: the run's .Params.txt names the variant datasets it searched (Genome_idx):
        a check given fewer VCFs of a chromosome (one of two datasets that share samples
        passes or fails by which) or more (a dataset the search never read) is refused;
        one VCF per dataset is checked."""
        case, info = world_deletion(self.tmp)
        params = os.path.join(case.res, ".Params.txt")
        with open(params, "w") as fh:
            fh.write("Genome_idx\tNRG_2_hg38+one,NRG_2_hg38+two\n")
        rc, _found = case.run()
        self.assertEqual(rc, 2)
        with open(params, "w") as fh:
            fh.write("Genome_idx\tNRG_2_hg38+one\n")
        self.assertClean(case)
        other_vcf = os.path.join(self.tmp, "w.vcf.gz")
        write_vcf(other_vcf, "chr21", [info["rec"][:3] + ({"T01": "1|0"},)],
                  samples=["T01", "T02"])
        rc, _found = case.run(extra=["--vcf", other_vcf])
        self.assertEqual(rc, 2)
        with open(params, "w") as fh:
            fh.write("Genome_idx\tNRG_2_hg38+one,NRG_2_hg38+two\n")
        rc, found = case.run(extra=["--vcf", other_vcf])
        self.assertNotEqual(rc, 2, found)

    def test_inputs_other_than_the_searchs_are_refused(self):
        """R7-F4's sweep: the budget, the PAM, the guides and the samples list the run
        records (.Params.txt, .guides.txt, .sampleID.txt) must be the check's: any other is
        refused (exit 2), never half-checked."""
        case, _info = world_deletion(self.tmp)
        params = os.path.join(case.res, ".Params.txt")
        for text, want in (("Mismatches\t3\n", 2), ("DNA\t0\n", 2), ("RNA\t2\n", 2),
                           ("Max_total_edits\t5\n", 2), ("Pam\tNGG\n", 2),
                           ("Mismatches\t4\nDNA\t1\nRNA\t1\nMax_total_edits\t6\nPam\tNRG\n",
                            0)):
            with open(params, "w") as fh:
                fh.write(text)
            rc, _found = case.run()
            self.assertEqual(rc, want, text)
        with open(os.path.join(case.res, ".guides.txt"), "w") as fh:
            fh.write("TTTTAAAATTTTAAAATTTTNNN\n")
        self.assertEqual(case.run()[0], 2)
        with open(os.path.join(case.res, ".guides.txt"), "w") as fh:
            fh.write(GUIDE + "NNN\n")
        self.assertClean(case)
        # some of the run's guides: a check of those (rows are keyed by guide); the rows
        # of the others are counted, not checked (they are no row of a guide given)
        g2 = "TTTTAAAATTTTAAAATTTT"
        with open(os.path.join(case.res, ".guides.txt"), "w") as fh:
            fh.write(GUIDE + "NNN\n" + g2 + "NNN\n")
        rows = [_info["rows"][0], set_field(_info["rows"][0], start=_info["start"] + 1)]
        rows[1][0] = g2 + "NNN"
        write_results(case.res, rows)
        found = self.assertClean(case)
        self.assertEqual([(f[1], f[2][:22]) for f in found],
                         [("guide-not-given", "3 row block(s) of the ")], found)
        write_results(case.res, _info["rows"])
        sid = os.path.join(case.res, ".sampleID.txt")
        with open(sid, "w") as fh:
            fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
            fh.write("".join("%s\tP\tSP\tfemale\n" % x for x in SAMPLES))
        panel = os.path.join(self.tmp, "panel.txt")
        with open(panel, "w") as fh:
            fh.write("".join("%s\n" % x for x in SAMPLES[:4]))
        self.assertEqual(case.run(extra=["--samples", panel])[0], 2)
        shutil.copy(sid, panel)
        self.assertClean(case, extra=["--samples", panel])
        # an analysis mode other than the one the run's command line says
        with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
            fh.write("input_command\tcrisprme.py complete-search --per-sample\n")
        self.assertEqual(case.run(extra=["--analysis", "population-level"])[0], 2)
        self.assertClean(case, extra=["--analysis", "per-sample"])

    def assertClean(self, case, mode="soundness", **kw):
        rc, found = case.run(mode=mode, **kw)
        self.assertEqual((rc, gating(found)), (0, []), found)
        return found


def world_union_bridge():
    """R8-F1: a run of AGG whose target at ideal offset 0 needs the MNV R (S03's), at 3
    (X) both PAM bases that two SNVs of two samples restore (S01's, S02's), and at 6 (B)
    one of them: X is a row of a population-level run (the worst case of each alignment
    over the dataset's alleles) and of no per-sample one (no individual holds both), so
    only population-level do T1 (Cluster Position s) and B (s+6) share a merge cluster."""
    g = "AGGAGGAGGAGGAGGAGGAG"
    ideal = "AGG" * 16
    run = list(ideal)
    run[1] = run[2] = "T"      # R (TT>GG) restores them: T1
    run[24] = "C"              # X's PAM (GAG): S01's SNV restores the A
    run[25] = "T"              # and S02's the last G
    s0 = len(LONG) + len(PRE)
    seq = LONG + PRE + "".join(run) + "TTTCATTTCA" + POST
    r = (s0 + 2, "TT", "GG", {"S03": "1|0"})
    snv1 = (s0 + 25, "C", "A", {"S01": "1|0"})
    snv2 = (s0 + 26, "T", "G", {"S02": "0|1"})
    b0 = s0 + 6
    alt = ideal[6:29]
    row_b = row_fields("chr21", b0, "+", g + "NNN", seq[b0:b0 + 23], alt, alt[-3:], 0, 0,
                       "X", "%s,%s" % (tok("chr21", snv1), tok("chr21", snv2)), "S01,S02",
                       spacer=g + "NNN")
    return seq, [r, snv1, snv2], [row_b], g


def world_overlap_slot(kind):
    """R8-F2 (the review's adv_r8a, adv_r8b, adv_r8d): the deletion world (R, the C after
    GUIDE[9]; S01 1|0, S02 0|1), with alleles of S01's slot that overlap each other:
    'second' a second deletion D2 inside the target and a SNV on the base D2 deletes;
    'before' two deletions that both delete the base right before the target (they change
    no base of it); 'before+snv' those and an unnamed SNV inside the target; 'split' both
    alts of a split SNV site inside the target (A>C and A>T phased onto the slot)."""
    lp = LONG + PRE
    seq = lp + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
    start = len(lp)
    rec = (start + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
    recs = [rec]
    d0 = start + 15
    q = start + 14
    if kind == "second":
        recs += [(d0 + 1, seq[d0:d0 + 2], seq[d0], {"S01": "1|0"}),
                 (d0 + 2, seq[d0 + 1], other(seq[d0 + 1]), {"S01": "1|0"})]
    if kind.startswith("before"):
        recs += [(start - 5, seq[start - 6:start], seq[start - 6], {"S01": "1|0"}),
                 (start - 1, seq[start - 2:start], seq[start - 2], {"S01": "1|0"})]
    if kind == "before+snv":
        recs.append((q + 1, seq[q], other(seq[q]), {"S01": "1|0"}))
    if kind == "split":
        b = seq[q]
        x1 = other(b)
        x2 = next(c for c in "ACGT" if c not in (b, x1))
        recs += [(q + 1, b, x1, {"S01": "1|0"}), (q + 1, b, x2, {"S01": "1|0"})]
    recs.sort(key=lambda x: x[0])

    def row(samples):
        return row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                          GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), samples)
    return seq, recs, row, start


class ReviewRound8(Base):
    # --- R8-F1: the analysis mode is the one the run ran (.search_mode) ------------- #
    def test_the_mode_is_the_runs_own_search_mode(self):
        """R8-F1: complete-search writes the mode it ran to .search_mode, whatever chose
        it; --ci-cd-test chooses none (complete-test's searches take the default,
        population-level on the self-test's index). A --per-sample on what has no
        genotypes runs population-level, which .search_mode says; markers no run writes
        together (per-sample beside --population-level, .alt_mode against the command
        line's flag) are refused."""
        res = os.path.join(self.tmp, "run")

        def kind(search=None, alt=None, cmd=None):
            shutil.rmtree(res, ignore_errors=True)
            os.makedirs(res)
            for name, text in ((".search_mode", search), (".alt_mode", alt),
                               (".command_line.txt", cmd)):
                if text is not None:
                    with open(os.path.join(res, name), "w") as fh:
                        fh.write(text)
            out = C.run_kind(res)
            return out[0], out[1], out[6]
        cicd = "input_command\tcrisprme.py complete-search --mm 4 --ci-cd-test"
        self.assertEqual(kind("population-level", "no-alt-alignments", cicd),
                         ("population-level", False, ".search_mode"))
        self.assertEqual(kind(None, "no-alt-alignments", cicd)[:2], (None, False))
        self.assertEqual(kind("population-level", "no-alt-alignments",
                              "crisprme.py complete-search --per-sample"),
                         ("population-level", False, ".search_mode"))
        self.assertEqual(kind(None, None, "crisprme.py complete-search --per-sample"),
                         ("per-sample", True, "command line"))
        self.assertEqual(kind(None, None, "crisprme.py complete-search --population-level"),
                         ("population-level", False, "command line"))
        self.assertEqual(kind(None, None, "crisprme.py complete-search --per-sample "
                                          "--no-alt-alignments")[:2], ("per-sample", False))
        self.assertEqual(kind(None, None, None), (None, None, None))
        for bad in (("per-sample", None, "crisprme.py complete-search --population-level"),
                    ("per-sample", "alt-alignments",
                     "crisprme.py complete-search --no-alt-alignments")):
            with self.assertRaises(ValueError):
                kind(*bad)

    def test_a_ci_cd_test_run_is_checked_in_the_mode_it_ran(self):
        """R8-F1: a run without alternative alignments made with --ci-cd-test, as the
        self-test's: its .search_mode (population-level there) says how its SNP rows chain
        a merge cluster. T1's row is gone and only B's is left: population-level, X joins
        them (a worst-case row of both SNVs), so the cluster holds B, as good as T1; a
        per-sample run writes no X, so T1's cluster has lost its row. The run that records
        no mode is refused unless --analysis gives one; one given against the run's
        .search_mode is refused."""
        seq, recs, rows, g = world_union_bridge()
        budget = ("1", "0", "0", "1")
        for search, cmd, want in (
                ("population-level", "--ci-cd-test", 0),
                ("per-sample", "--per-sample --no-alt-alignments --ci-cd-test", 1)):
            case = Case(self.tmp, seq, recs, rows)
            run_markers(case.res, search, alt=False)
            with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
                fh.write("input_command\tcrisprme.py complete-search %s\n" % cmd)
            rc, found = case.run(mode="completeness", guide=g, budget=budget)
            self.assertEqual(rc, want, (search, found))
            if want:
                self.assertEqual([(f[1], f[3]) for f in gating(found)],
                                 [("missing-locus", "S03")], found)
        # no mode recorded (no .search_mode, no mode flag): refused, unless given
        case = Case(self.tmp, seq, recs, rows)
        with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
            fh.write("no-alt-alignments")
        with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
            fh.write("input_command\tcrisprme.py complete-search --ci-cd-test\n")
        argv = ["--results", case.res, "--reference", os.path.join(self.tmp, "Genomes"),
                "--vcf", case.vcf, "--pam", case.pam, "--mm", "1", "--bDNA", "0",
                "--bRNA", "0", "--max-total-edits", "1", "--chrom", "chr21", "--guide", g,
                "--out", os.path.join(self.tmp, "o.tsv")]
        self.assertEqual(C.main(argv + ["--mode", "completeness"]), 2)
        self.assertEqual(C.main(argv + ["--mode", "soundness", "--snp-rows", "report"]), 2)
        self.assertEqual(C.main(argv + ["--mode", "soundness"]), 0)
        self.assertEqual(C.main(argv + ["--mode", "completeness", "--analysis",
                                        "population-level"]), 0)
        self.assertEqual(C.main(argv + ["--mode", "completeness", "--analysis",
                                        "per-sample"]), 1)
        # a run that records nothing of its alternative alignments: completeness refused
        os.remove(os.path.join(case.res, ".alt_mode"))
        os.remove(os.path.join(case.res, ".command_line.txt"))
        self.assertEqual(C.main(argv + ["--mode", "completeness", "--analysis",
                                        "population-level"]), 2)
        # an --analysis against the run's .search_mode
        run_markers(case.res, "population-level")
        self.assertEqual(C.main(argv + ["--mode", "soundness", "--analysis",
                                        "per-sample"]), 2)

    # --- R8-F2: a slot whose other alleles overlap each other ----------------------- #
    def test_a_listed_carrier_the_documented_rules_take_off_fails(self):
        """R8-F2 (a, c): S01's slot has no one sequence where the target reads -- a second
        deletion inside it and a SNV on the base it deletes; both alts of a split SNV site
        -- and every reading of it but the one with neither lacks the target. METHODS
        section 4 lists no such carrier (a slot that cannot be read through its second
        indel; the last record of the split site): listed, it fails."""
        for kind, cause in (("second", "second-indel"), ("split", "last-alt-wins")):
            seq, recs, row, _start = world_overlap_slot(kind)
            rc, found = Case(self.tmp, seq, recs, [row("S01,S02")]).run()
            self.assertEqual(rc, 1, (kind, found))
            self.assertEqual({(f[1], f[3]) for f in gating(found)},
                             {("carrier-off-by-rule", "S01")}, (kind, found))
            self.assertTrue(all(cause in f[2] for f in gating(found)), found)
            rc, found = Case(self.tmp, seq, recs, [row("S02")]).run()
            self.assertEqual((rc, gating(found)), (0, []), (kind, found))
            self.assertEqual({f[1] for f in found}, {"limit-overlapping-records"}, found)

    def test_a_carrier_with_no_reading_is_counted_never_demanded(self):
        """R8-F2's sweep: S01's slot of R holds a second deletion inside the target and a
        SNV on the base it deletes (no one sequence there; a target through that deletion
        needs two indel-path records): completeness demands nothing of S01 and counts it
        (limit-overlapping-records), never leaving it out in silence; S02, on a clean
        slot, is demanded."""
        seq, recs, _row, _start = world_overlap_slot("second")
        rc, found = Case(self.tmp, seq, recs, []).run(mode="completeness",
                                                      budget=("1", "0", "0", "1"))
        self.assertEqual({f[3] for f in gating(found)}, {"S02"}, found)
        self.assertIn(("limit-overlapping-records", "S01"),
                      {(f[1], f[3]) for f in found if f[0] == "no"}, found)

    def test_a_split_sites_last_record_lists_its_carrier(self):
        """R8-F2: both alts of a split SNV site on S01's slot (no reading): a row naming the
        last record, as the documented rule reads the slot, may list S01 (counted); a row
        naming the earlier one shows the run's stores hold it, against the documented
        limit (variant-shadowed), so the rule is no model there and nothing gates."""
        seq, recs, row, start = world_overlap_slot("split")
        q = start + 14
        split = [r for r in recs if r[0] == q + 1]
        rec = recs[[r[0] for r in recs].index(start + 10)]
        for which, want in ((split[-1], {"limit-overlapping-records"}),
                            (split[0], {"limit-overlapping-records", "variant-shadowed"})):
            alt = GUIDE[:13] + which[2].lower() + GUIDE[14:] + "TGG"
            r = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23], alt,
                           "TGG", 1, 0, "X", "%s,%s" % (tok("chr21", rec), tok("chr21", which)),
                           "S01")
            rc, found = Case(self.tmp, seq, recs, [r]).run()
            self.assertEqual((rc, gating(found)), (0, []), (which, found))
            self.assertEqual({f[1] for f in found}, want, (which, found))

    def test_an_overlap_that_changes_no_target_base_is_judged(self):
        """R8-F2 (b, d): two deletions of S01's slot overlap right before the target and
        change none of its bases, so every reading of them holds the same target: S01 is
        judged as on a clean slot. Listed beside an unnamed SNV inside the target, it fails
        (variant-missing); holding the target and left off, it fails (carrier-missing), and
        completeness demands its target."""
        seq, recs, row, start = world_overlap_slot("before+snv")
        rc, found = Case(self.tmp, seq, recs, [row("S01,S02")]).run()
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("variant-missing", "S01")}, found)
        seq, recs, row, start = world_overlap_slot("before")
        case = Case(self.tmp, seq, recs, [row("S02")])
        rc, found = case.run()
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-missing", "S01")}, found)
        self.assertClean(Case(self.tmp, seq, recs, [row("S01,S02")]))
        write_results(case.res, [])
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual(sorted((f[1], f[3]) for f in gating(found)),
                         [("missing-row", "S01"), ("missing-row", "S02")], found)

    def test_zone_readings(self):
        """R8-F2: the readings of an overlap zone are its alleles no two of which meet, none
        included; a zone none of whose alleles is named or touches what the slot reads
        from the target's start leaves every verdict the same."""
        a = C.Allele(10, "ACG", "A", 0, 1)    # deletes 10, 11 (0-based)
        b = C.Allele(11, "C", "T", 1, 1)      # on a base a deletes
        c = C.Allele(12, "G", "GT", 2, 1)     # meets a, not b
        self.assertEqual(sorted(len(x) for x in C.zone_readings([a, b, c])),
                         [0, 1, 1, 1, 2])
        self.assertIn((b, c), C.zone_readings([a, b, c]))

    # --- R8-F3: a slot read without its artefact is enumerated as a clean slot ----- #
    def test_a_slot_read_without_its_artefact_chains_no_two_record_target(self):
        """R8-F3 (the review's adv_r8c): X needs R and a deletion D2 of S01's slot, two
        indel-path records (CRISPRme writes no such row), and a SNV at R's own position on
        that slot is a phasing artefact METHODS reads it without. Read so, the slot is a
        clean one: X needs a second indel-path record, so it chains no cluster, and a
        table that lost S03's T1 has lost T1's cluster."""
        g = "AGGAGGAGGAGGAGGAGGAG"
        ideal = list("AGG" * 16)
        ideal[5] = "T"
        ideal[24] = "T"
        s0 = len(LONG) + len(PRE)
        seq = LONG + PRE + "".join(ideal[:4] + ideal[5:]) + "TTTCATTTCA" + POST
        r = (s0 + 4, seq[s0 + 3], seq[s0 + 3] + "G", {"S01": "1|0", "S03": "1|0"})
        sp = s0 + 23
        snv = (sp + 1, "T", "A", {"S02": "0|1"})
        d2 = (sp, seq[sp - 1:sp + 1], seq[sp - 1], {"S01": "1|0"})
        b0 = s0 + 5
        alt = seq[b0:sp] + "A" + seq[sp + 1:b0 + 23]
        row_b = row_fields("chr21", b0, "+", g + "NNN", seq[b0:b0 + 23], alt, alt[-3:], 0,
                           0, "X", tok("chr21", snv), "S02", spacer=g + "NNN")
        t1 = "".join(ideal[0:23])
        row_t1 = row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23], mark(g + "NNN", t1),
                            t1[-3:], 1, 0, "X", tok("chr21", r), "S03", spacer=g + "NNN")
        art = (s0 + 4, seq[s0 + 3], other(seq[s0 + 3]), {"S01": "1|0"})
        for extra, rows, want in (([], [row_t1, row_b], 0), ([], [row_b], 1),
                                  ([art], [row_t1, row_b], 0), ([art], [row_b], 1)):
            case = Case(self.tmp, seq, sorted([r, snv, d2] + extra), rows)
            run_markers(case.res, "population-level")
            rc, found = case.run(mode="completeness", guide=g, budget=("1", "0", "0", "1"))
            self.assertEqual(rc, want, (extra, found))
            if want:
                self.assertEqual([(f[1], "S03" in f[3].split(",")) for f in gating(found)],
                                 [("missing-locus", True)], found)
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
        ctx = ReviewRound7._ctx(self, g, mm=1)
        rr = win.by_key[r[:3]][0]
        e = C.record_demands(win, ctx, rr, {"S01", "S03"}, ctx.full(g),
                             max(0, rr.lo - 80), rr.hi + 80, (rr.edit_lo - 26, rr.edit_hi + 26))
        self.assertEqual(({c - s0 for c in e["cps"]}, e["cps_open"]), ({0}, {}))

    # --- R8-N2: a check of some of a run's guides says which it left out ------------ #
    def test_a_completeness_check_of_some_guides_names_the_others(self):
        """R8-N2: completeness of a subset of the run's guides says, per guide left out,
        that none of its targets is demanded."""
        case, _info = world_deletion(self.tmp)
        g2 = "TTTTAAAATTTTAAAATTTT"
        with open(os.path.join(case.res, ".guides.txt"), "w") as fh:
            fh.write(GUIDE + "NNN\n" + g2 + "NNN\n")
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual([(f[1], f[2][:44]) for f in found if f[1] == "guide-not-given"],
                         [("guide-not-given",
                           "completeness demands none of the targets of ")], found)

    def assertClean(self, case, mode="soundness", **kw):
        rc, found = case.run(mode=mode, **kw)
        self.assertEqual((rc, gating(found)), (0, []), found)
        return found


def world_anchor_snv(kind):
    """R9-F1 (the review's adv_r9b2 / adv_r9c): an SNV A at an indel's own POS, its alt the
    guide's base at the anchor (whose reference base mismatches the guide), and a call of
    unknown phase inside the target, so that the truth is the PUTATIVE rule. CRISPRme's
    fake-indel window writes the ALT verbatim, the anchor included, and codes flank bases
    only: the anchor is never a PUTATIVE choice, whatever A is. 'deletion': S01 holds the
    deletion D, A and an SNV X unphased (0/1 each); 'insertion': S01 holds the insertion
    1|0, A phased on its other haplotype (0|1) and an SNV Y with a missing allele (.|1).
    Returns (seq, records, the right row, its start, the anchor's coordinate)."""
    G = GUIDE
    if kind == "deletion":
        anchor_ref, x_ref = other(G[9]), other(G[3])
        seq = PRE + G[:3] + x_ref + G[4:9] + anchor_ref + "C" + G[10:] + "TGG" + POST
        p_del, p_x = len(PRE) + 10, len(PRE) + 4
        D = (p_del, anchor_ref + "C", anchor_ref, {"S01": "0/1"})
        A = (p_del, anchor_ref, G[9], {"S01": "0/1", "S03": "0|1"})
        X = (p_x, x_ref, G[3], {"S01": "0/1"})
        recs = sorted([D, A, X], key=lambda r: (r[0], len(r[1])))
        s = len(PRE)
        tgt = G[:9] + anchor_ref + G[10:] + "TGG"
        row = row_fields("chr21", s, "+", G + "NNN", seq[s:s + 23],
                         mark(G + "NNN", tgt), "TGG", 1, 0, "X",
                         "%s,%s" % (tok("chr21", D), tok("chr21", X)), "S01")
        return seq, recs, row, s, p_del - 1
    pre, post = LONG + PRE, POST + LONG
    X0 = other(G[0])
    y_ref = other(G[12])
    seq = pre + X0 + G[5:12] + y_ref + G[13:] + "TGG" + post
    a0 = len(pre)
    I = (a0 + 1, X0, X0 + G[1:5], {"S01": "1|0"})
    A = (a0 + 1, X0, G[0], {"S01": "0|1", "S03": "1|0"})
    Y = (a0 + 1 + 7 + 1, y_ref, G[12], {"S01": ".|1"})
    recs = sorted([I, A, Y], key=lambda r: (r[0], -len(r[2])))
    tgt = X0 + G[1:] + "TGG"
    row = row_fields("chr21", a0, "+", G + "NNN", seq[a0:a0 + 23], mark(G + "NNN", tgt),
                     "TGG", 1, 0, "X", "%s,%s" % (tok("chr21", I), tok("chr21", Y)), "S01")
    return seq, recs, row, a0, a0


def world_sites_only_chain(genotyped_r=True, anchor_snv=False, t1_row=False):
    """R9-F2 / R9-F4 (the review's adv_r9f / adv_r9e): an (AGG)n run, guide AGGAGG...,
    NRG, a population-level run without alternative alignments beside a sites-only VCF.
    The reference lacks ideal[4] (R, an insertion anchored on ideal[3], puts it back);
    ideal[20] = T (S, T>G: S02's SNP row B, 0 edits, the table's one row). T1 (R's target,
    S01's, Cluster Position s) and B are one merge cluster only through W, the sites-only
    dataset's worst-case row of R's window with S (both in the sites-only VCF), at s+3.
    With ``anchor_snv``: ideal[3] = C (one mismatch on T1; budget 1 mismatch), and the
    sites-only VCF holds A, an SNV at R's own POS whose alt is the guide's base: a window
    coded at its anchor would hold a target CRISPRme's window (anchor uncoded) cannot."""
    g = "AGGAGGAGGAGGAGGAGGAG"
    ideal = list("AGG" * 16)
    if anchor_snv:
        ideal[3] = "C"
    ideal[20] = "T"
    s0 = len(LONG) + len(PRE)
    seq = LONG + PRE + "".join(ideal[:4] + ideal[5:]) + "TTTCATTTCA" + POST
    anchor = s0 + 3
    R = (anchor + 1, seq[anchor], seq[anchor] + "G", {"S01": "1|0"})
    S = (s0 + 20, "T", "G", {"S02": "0|1"})
    b0 = s0 + 5
    balt = seq[b0:s0 + 19] + "G" + seq[s0 + 20:b0 + 23]
    row_b = row_fields("chr21", b0, "+", g + "NNN", seq[b0:b0 + 23], balt, balt[-3:], 0, 0,
                       "X", tok("chr21", S), "S02", spacer=g + "NNN")
    rows = [row_b]
    if t1_row:
        t1 = "".join(ideal[0:23])
        rows.insert(0, row_fields("chr21", s0, "+", g + "NNN", seq[s0:s0 + 23],
                                  mark(g + "NNN", t1), t1[-3:], int(anchor_snv), 0, "X",
                                  tok("chr21", R), "S01", spacer=g + "NNN"))
    recs = [S] + ([R] if genotyped_r else [])
    mega = [R[:3], (anchor + 1, "C", "A")] if anchor_snv else [R[:3], S[:3]]
    budget = ("1", "0", "0", "1") if anchor_snv else ("0", "0", "0", "0")
    return g, seq, recs, rows, mega, budget, R


class ReviewRound9(Base):
    # --- R9-F1: no PUTATIVE choice on the record's own REF span ------------------------ #
    def test_the_putative_rule_has_no_choice_on_the_records_own_span(self):
        """R9-F1: the anchor of D carries no code column (CRISPRme's window writes the ALT
        verbatim, and an SNV at the record's POS is the other allele of a multi-allelic
        heterozygote, on the other haplotype): the row CRISPRme documents -- D and X, the
        reference anchor, 1 mismatch, listing S01 -- is right (it failed carrier-not-
        strongest three times), and one that writes A's alt at the anchor is no sequence
        of D and X. The insertion world (A phased on S01's other haplotype) likewise."""
        for kind in ("deletion", "insertion"):
            seq, recs, row, s, anchor = world_anchor_snv(kind)
            tmp = os.path.join(self.tmp, kind)
            os.makedirs(tmp)
            case = Case(tmp, seq, recs, [row])
            self.assertClean(case)
            rows = C.read_rows(results_table(case), "integrated")
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(case.vcf)])
            named = [win.by_key[(p, r, a)][0] for p, r, a in
                     (C.parse_token(t)[1:] for t in rows[0]["var"].split(","))]
            main = next(a for a in named if a.kind != "snv")
            tgt = C.row_target(win, rows[0], named, main)[0]
            self.assertTrue(all(q != anchor for _i, q in tgt.codes), (kind, tgt.codes))
            self.assertTrue(tgt.codes, kind)       # the SNV in the flank is still a code
            G = GUIDE
            alt = row[2 + C.BLOCKS.index("highest_CFD") * 13 + 4]
            i = anchor - s
            wrong = alt[:i] + (G[i] if kind == "deletion" else G[0]) + alt[i + 1:]
            write_results(case.res, [set_field(row, alt=wrong, mm=0)])
            rc, found = case.run()
            self.assertEqual(rc, 1, (kind, found))

    def test_the_best_alignments_putative_demand_is_kept(self):
        """R9-F1, completeness: the best alignment's PUTATIVE demand no longer disappears
        (no choice at the anchor), so the target reported only by a worse alignment is
        missing-best-alignment, and the best alignment's row satisfies it."""
        seq, recs, row, s, _anchor = world_anchor_snv("deletion")
        G = GUIDE
        tgt = G[:9] + other(G[9]) + G[10:] + "TGG"
        als = sorted(C.align_pairs(G + "NRG", tgt, 4, 1, 1, 6, 3, False),
                     key=lambda x: (x[2] + x[3] + x[4], x[2]))
        best = als[0]
        worse = next(a for a in als if a[2] + a[3] + a[4] > sum(best[2:5]))

        def mkrow(a):
            sa, ta, mm, dn, rn = a
            ref = list(seq[s:s + len(ta.replace("-", ""))])
            for k, c in enumerate(ta):
                if c == "-":
                    ref.insert(k, "-")
            bt = ",".join([t for t, c in (("RNA", rn), ("DNA", dn)) if c]) or "X"
            return row_fields("chr21", s, "+", sa[:-3] + "NNN", "".join(ref),
                              mark(sa[:-3] + "NNN", ta), ta.replace("-", "")[-3:], mm,
                              dn + rn, bt, row[2 + 11], "S01")
        case = Case(self.tmp, seq, recs, [mkrow(worse)])
        rc, found = case.run(mode="completeness")
        self.assertTrue(any(f[1] == "missing-best-alignment" and f[3] == "S01"
                            for f in gating(found)), found)
        write_results(case.res, [mkrow(best)])
        rc, found = case.run(mode="completeness")
        self.assertFalse(any(f[1] == "missing-best-alignment" for f in found), found)
        # (S01's targets of other alignments -- under which its strongest bases differ --
        # are demanded too, and this one-row table lacks them: only the best one's is
        # asked about here)
        self.assertFalse(any("best %d edits" % sum(best[2:5]) in f[2]
                             for f in gating(found)), found)

    # --- R9-F2 / R9-F4: a record's chain beside a sites-only dataset ----------------- #
    def test_a_records_chain_holds_its_sites_only_rows_whoever_reads_it_first(self):
        """R9-F2 (adv_r9f): a record both a genotyped and a sites-only VCF hold chains
        the sites-only dataset's worst-case rows of its window too, whichever record's
        work reads it first: the right table (B only, as good as T1) passes; R in the
        sites-only VCF alone passes too (control)."""
        for genotyped_r in (True, False):
            g, seq, recs, rows, mega, budget, _R = world_sites_only_chain(genotyped_r)
            tmp = os.path.join(self.tmp, str(genotyped_r))
            os.makedirs(tmp)
            case = Case(tmp, seq, recs, rows)
            run_markers(case.res, "population-level", alt=False)
            so = os.path.join(tmp, "mega.vcf.gz")
            write_sites_only_vcf(so, "chr21", mega)
            rc, found = case.run(mode="completeness", guide=g, budget=budget,
                                 extra=("--vcf", so))
            self.assertEqual((rc, gating(found)), (0, []), (genotyped_r, found))

    def test_a_sites_only_window_codes_no_base_of_its_record(self):
        """R9-F4 (adv_r9d, adv_r9e): a sites-only dataset's window writes its record's ALT
        verbatim and codes flank bases only, so an SNV at the record's own POS is in no
        row of it: such a window gives no Cluster Position there, and a dropped T1 that
        only such a phantom would chain to B is missing-locus; T1's row, reported,
        passes."""
        g, seq, recs, rows, mega, budget, R = world_sites_only_chain(anchor_snv=True)
        so = os.path.join(self.tmp, "mega.vcf.gz")
        write_sites_only_vcf(so, "chr21", mega)
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(so)])

        class A(object):
            mm, bdna, brna, max_total_edits = 1, 0, 0, 1
        ctx = C.Ctx(A(), "NRG", 3, False, [g])
        rd = win.by_key[R[:3]][0]
        cps = C.record_cps(win, ctx, g, rd)[0]
        anchor = R[0] - 1
        self.assertNotIn(anchor - 3 + 3, cps)    # X, the anchor-coded phantom (s + 3)
        case = Case(self.tmp, seq, recs, rows)
        run_markers(case.res, "population-level", alt=False)
        rc, found = case.run(mode="completeness", guide=g, budget=budget,
                             extra=("--vcf", so))
        self.assertEqual(rc, 1, found)
        self.assertIn("missing-locus", classes(gating(found)), found)
        _g, _seq, _recs, rows2, _m, _b, _R = world_sites_only_chain(anchor_snv=True,
                                                                     t1_row=True)
        write_results(case.res, rows2)
        rc, found = case.run(mode="completeness", guide=g, budget=budget,
                             extra=("--vcf", so))
        self.assertEqual((rc, gating(found)), (0, []), found)

    # --- R9-F3: a search without bulges records Genome_idx None ---------------------- #
    def test_a_search_without_bulges_reads_its_datasets_from_its_vcf_list(self):
        """R9-F3: crisprme.py writes 'Genome_idx None' for every search without bulges,
        whatever it searched; its datasets are its --vcf list's lines, read in its working
        directory (the folder holding Results/<output>). One dataset, one VCF: checked
        (it was refused, 'the reference alone'); a list of two: refused; no --vcf: the
        reference alone, refused; a list not found: the VCFs given are taken."""
        case, info = world_deletion(self.tmp)
        work = os.path.join(self.tmp, "work")
        os.makedirs(os.path.join(work, "Results"))
        dest = os.path.join(work, "Results", "b0")
        shutil.move(case.res, dest)
        case.res = dest
        with open(os.path.join(dest, ".Params.txt"), "w") as fh:
            fh.write("Genome_idx\tNone\nMismatches\t0\nDNA\t0\nRNA\t0\n"
                     "Max_total_edits\t0\n")
        b0 = ("0", "0", "0", "0")

        def cmd(text):
            with open(os.path.join(dest, ".command_line.txt"), "w") as fh:
                fh.write("input_command\t/opt/conda/bin/crisprme.py complete-search "
                         "--genome Genomes/hg38 %s--mm 0 --bDNA 0 --bRNA 0 --output b0 "
                         "--per-sample\n" % text)
        lv = os.path.join(work, "lv.txt")
        with open(lv, "w") as fh:
            fh.write("VCFs/hg38_one/\n")
        cmd("--vcf lv.txt ")
        self.assertEqual(C.run_datasets(dest)[0], ["hg38_one"])
        for mode in ("soundness", "completeness"):
            self.assertClean(case, mode=mode, budget=b0)
        with open(lv, "w") as fh:
            fh.write("VCFs/hg38_one\nVCFs/hg38_two\n")
        rc, _found = case.run(budget=b0)
        self.assertEqual(rc, 2)
        cmd("")
        self.assertEqual(C.run_datasets(dest)[0], [])
        rc, _found = case.run(budget=b0)
        self.assertEqual(rc, 2)
        cmd("--vcf elsewhere/lv.txt ")
        self.assertIsNone(C.run_datasets(dest)[0])
        self.assertClean(case, budget=b0)

    # --- the PUTATIVE rule of a search that scores no CFD ---------------------------- #
    def test_without_cfd_the_putative_choice_is_the_reference_then_acgt(self):
        """METHODS section 4: a search that scores no CFD (here a 19-nt guide) takes, where
        the guide's base does not decide, the reference base where the carrier's alleles
        allow it, else the first of A, C, G, T they allow. S01 (unphased deletion and SNV,
        neither allele the guide's) is listed for the reference base; a row listing it
        for the SNV's alt fails carrier-not-strongest, and with the row lost completeness
        demands the reference-base target. With the CFD tables missing for a search that
        scores CFD, the choice stays open (ReviewRound3's tie test)."""
        g = GUIDE[1:]
        r = other(g[4])
        b = next(x for x in "ACGT" if x not in (r, g[4]))
        seq = PRE + g[:4] + r + g[5:9] + "C" + g[9:] + "TGG" + POST
        D = (len(PRE) + 9, g[8] + "C", g[8], {"S01": "0/1", "S02": "1|0"})
        X = (len(PRE) + 5, r, b, {"S01": "0/1"})
        s = len(PRE)

        def row(base, var, who):
            t = g[:4] + base + g[5:] + "TGG"
            return row_fields("chr21", s, "+", g + "NNN", seq[s:s + 22], mark(g + "NNN", t),
                              "TGG", 1, 0, "X", var, who, spacer=g + "NNN")
        right = [row(r, tok("chr21", D), "S01,S02")]
        case = Case(self.tmp, seq, [X, D], right)
        budget = ("1", "0", "0", "1")
        self.assertClean(case, budget=budget, guide=g)
        self.assertClean(case, mode="completeness", budget=budget, guide=g)
        wrong = [row(r, tok("chr21", D), "S02"),
                 row(b, "%s,%s" % (tok("chr21", X), tok("chr21", D)), "S01")]
        write_results(case.res, wrong)
        rc, found = case.run(budget=budget, guide=g)
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-not-strongest", "S01"), ("carrier-missing", "S01")},
                         found)
        rc, found = case.run(mode="completeness", budget=budget, guide=g)
        self.assertIn("missing-putative-row", classes(gating(found)), found)
        # where the reference base is not an allele of S01 (1/2: two alts), the first of
        # A, C, G, T the alleles allow, in the target's orientation
        b2 = next(x for x in "ACGT" if x not in (r, g[4], b))
        X2 = (len(PRE) + 5, r, "%s,%s" % (b, b2), {"S01": "1/2"})
        first = min(b, b2)
        case = Case(self.tmp, seq, [X2, D],
                    [row(first, "%s,%s" % (tok("chr21", (X2[0], r, first)),
                                           tok("chr21", D)), "S01"),
                     row(r, tok("chr21", D), "S02")])
        self.assertClean(case, budget=budget, guide=g)

    def test_without_cfd_the_choice_reads_bases_in_the_targets_orientation(self):
        """The no-CFD order (the reference base, else A, C, G, T) is read in the target's
        orientation, as CRISPRme's code columns are: on the - strand the + alleles A and
        T read T and A, so A (+ T) comes first."""
        class Ctx0(object):
            cfd = None
            cfd_applies = False
            pam_len = 3
            pam_5prime = False
            pam_motif = "NRG"

            @staticmethod
            def pam_cols(n):
                return range(n - 3, n)

        class T0(object):
            pass
        t = T0()
        t.codes = [(5, 105)]
        t.crrna = GUIDE + "NNN"
        t.strand = "-"
        t.aligned = GUIDE + "TGG"
        refseq = "G" * 200
        # - strand: base 5 faces the guide's G (column 17); + A and T read T and A
        per = C.strongest(Ctx0(), t, {5: {"A", "T"}}, refseq)[0]
        self.assertEqual(per[5], {"T"})
        t.strand = "+"     # base 5 faces the guide's A (column 5): C and T read as they are
        per = C.strongest(Ctx0(), t, {5: {"C", "T"}}, refseq)[0]
        self.assertEqual(per[5], {"C"})
        per = C.strongest(Ctx0(), t, {5: {"C", "G"}}, refseq)[0]
        self.assertEqual(per[5], {"G"})     # the reference base where allowed

    # --- Track B in the per-sample chain --------------------------------------------- #
    def test_a_per_sample_snp_row_whose_slot_holds_a_changing_indel_is_open(self):
        """2.10.0 drops a per-sample SNP row's carrier whose slot of its SNPs holds an
        indel that changes the target (METHODS section 4): in a per-sample chain the SNP
        target of such a slot is a row CRISPRme may write or not (a repeat read through
        keeps the carrier), never a certain one -- unless another slot holds the SNP
        without the indel."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        snp = (len(PRE) + 6, other(GUIDE[5]), GUIDE[5])
        dele = len(PRE) + 12
        drec = (dele, seq[dele - 1:dele + 1], seq[dele - 1])
        cp = C.cluster_of(len(PRE), "+", 23, 23, False)
        for extra, sure in (({}, False), ({"S03": "1|0"}, True)):
            vcf = os.path.join(self.tmp, "v%d.vcf.gz" % sure)
            write_vcf(vcf, "chr21", [snp + (dict({"S02": "1|0"}, **extra),),
                                     drec + ({"S02": "1|0"},)])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])

            class A(object):
                mm, bdna, brna, max_total_edits = 0, 0, 0, 0
            ctx = C.Ctx(A(), "NRG", 3, False, [GUIDE])
            ctx.snp_union = False
            rows, _capped, dropped, _filt, _deg, _grd = C.snp_row_cps(win, ctx, GUIDE, ctx.full(GUIDE), 0,
                                                   len(seq))
            self.assertEqual((cp in rows, cp in dropped), (sure, not sure), extra)
            ctx.snp_union = True    # population-level: unchanged
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
            rows, _capped, dropped, _filt, _deg, _grd = C.snp_row_cps(win, ctx, GUIDE, ctx.full(GUIDE), 0,
                                                   len(seq))
            self.assertEqual((cp in rows, dropped), (True, set()), extra)

    def assertClean(self, case, mode="soundness", **kw):
        rc, found = case.run(mode=mode, **kw)
        self.assertEqual((rc, gating(found)), (0, []), found)
        return found


def world_filter_pam(filt):
    """R10-F1 (review-r10's LowQual world, site A): world_deletion's perfect target (D on
    S01's slot 0 and S02's slot 1) and P, a G>A SNV at the PAM's last base on S02's D
    slot, its FILTER ``filt``. CRISPRitz's enricher skips a line whose FILTER is not
    PASS or '.': with P LowQual no store CRISPRme reads holds it, and the row lists S02
    (CRISPRme does); with P PASS it does not."""
    seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
    p = len(PRE)
    D = (p + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
    P = (p + 24, "G", "A", {"S02": "0|1"}, filt)
    row = row_fields("chr21", p, "+", GUIDE + "NNN", seq[p:p + 23], GUIDE + "TGG", "TGG",
                     0, 0, "X", tok("chr21", D), "S01,S02")
    return seq, [D, P], row, D, P


def world_filter_second(filt):
    """R10-F1, the phased indel store: world_deletion's D (S01's slot 0, S02's slot 1) and
    J, a 1-bp deletion inside the perfect target on S01's D slot, its FILTER ``filt``.
    CRISPRme's phased indel store (build_indel_genotypes) keeps every indel of the VCF,
    whatever its FILTER, and its carrier resolution reads it as a second indel: S01 is
    taken off D's perfect target, which lists S02 alone (as the oracle's site O, with its
    second deletion LowQual, lists S22 alone)."""
    seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
    p = len(PRE)
    D = (p + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
    J = (p + 4, GUIDE[3:5], GUIDE[3], {"S01": "1|0"}, filt)
    row = row_fields("chr21", p, "+", GUIDE + "NNN", seq[p:p + 23], GUIDE + "TGG", "TGG",
                     0, 0, "X", tok("chr21", D), "S02")
    return seq, sorted([D, J], key=lambda r: r[0]), row, D, J


class ReviewRound10(Base):
    # --- R10-F1: the FILTER CRISPRitz's enricher reads ---------------------------------- #
    def test_a_window_reads_filter_as_crisprme_s_stores_do(self):
        """R10-F1: a genotyped line outside PASS and '.' (the exact strings CRISPRitz's
        enricher keeps: 'PASS;LowQual' is skipped too) keeps only its indels, unsearched
        (the phased indel store's: on their carriers' haplotypes, carrying no row);
        its SNVs and MNVs are in no store CRISPRme reads. A sites-only line keeps its
        SNVs (hg38_mega's index holds them) and drops its indels and MNVs."""
        seq = PRE + GUIDE + "TGG" + POST
        p = len(PRE)
        geno = os.path.join(self.tmp, "g.vcf.gz")
        recs = [(p + 2, seq[p + 1], "T" if seq[p + 1] != "T" else "A", {"S01": "1|0"},
                 "LowQual"),
                (p + 4, seq[p + 3:p + 5], seq[p + 3], {"S01": "1|0"}, "LowQual"),
                (p + 6, seq[p + 5], seq[p + 5] + "T," + other(seq[p + 5]),
                 {"S02": "1|2", "S03": "2|0"}, "PASS;LowQual"),
                (p + 8, seq[p + 7:p + 9], other(seq[p + 7]) + other(seq[p + 8]),
                 {"S04": "1|0"}, "VQSRTrancheSNP99.00to100.00"),
                (p + 11, seq[p + 10], other(seq[p + 10]), {"S05": "1|0"}, "."),
                (p + 13, seq[p + 12:p + 14], seq[p + 12], {"S06": "0|1"}, "PASS")]
        write_vcf(geno, "chr21", recs)
        sites = os.path.join(self.tmp, "s.vcf.gz")
        write_sites_only_vcf(sites, "chr21", [
            (p + 15, seq[p + 14], other(seq[p + 14]), "LOW_VQSLOD"),
            (p + 17, seq[p + 16:p + 18], seq[p + 16], "ExcHet"),
            (p + 19, seq[p + 18], seq[p + 18] + "G", "PASS")])
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(geno), C.Vcf(sites)])
        key = lambda r: (r[0], r[1], r[2])  # noqa: E731
        # the genotyped LowQual SNV: in no store, read as the reference
        self.assertNotIn(key(recs[0]), win.by_key)
        self.assertEqual(win.filtered_keys[key(recs[0])], "LowQual")
        # the genotyped LowQual deletion: on S01's haplotype, unsearched, no carrier
        d = win.by_key[key(recs[1])][0]
        self.assertEqual((d.searched, d.filtered), (False, "LowQual"))
        self.assertEqual(win.carriers(d), [])
        self.assertEqual(C.main_call(win, d.key, "S01"), (None, None))
        self.assertFalse(win.searched(d.key))
        self.assertTrue(C.on_slot(((1, 0), True), 0, win.calls[d.line]["S01"], d.number)
                        == "yes")
        # 'PASS;LowQual' is not PASS: the insertion stays (unsearched), the SNV alt goes
        ins = win.by_key[(p + 6, seq[p + 5], seq[p + 5] + "T")][0]
        self.assertFalse(ins.searched)
        self.assertNotIn((p + 6, seq[p + 5], other(seq[p + 5])), win.by_key)
        self.assertIsNone(win.line_alleles[ins.line][1])
        # the genotyped non-PASS MNV: in no store
        self.assertNotIn(key(recs[3]), win.by_key)
        # '.' and PASS: read as ever
        self.assertTrue(win.by_key[key(recs[4])][0].searched)
        self.assertIsNone(win.by_key[key(recs[4])][0].filtered)
        self.assertEqual(win.carriers(win.by_key[key(recs[5])][0]), ["S06"])
        # sites-only: the LOW_VQSLOD SNV kept (searched, flagged), the ExcHet deletion
        # dropped, the PASS insertion kept
        snv = win.by_key[(p + 15, seq[p + 14], other(seq[p + 14]))][0]
        self.assertEqual((snv.searched, snv.filtered), (True, "LOW_VQSLOD"))
        self.assertNotIn((p + 17, seq[p + 16:p + 18], seq[p + 16]), win.by_key)
        self.assertIn((p + 17, seq[p + 16:p + 18], seq[p + 16]), win.filtered_keys)
        self.assertTrue(win.sites_only((p + 19, seq[p + 18], seq[p + 18] + "G")))

    def test_a_snp_the_search_skips_leaves_its_carrier_listed(self):
        """R10-F1 (review-r10's site A): P LowQual, the row listing S02 -- as CRISPRme
        lists it -- passes in both modes (it failed variant-missing); P PASS, the same
        row fails for S02; a row naming P fails variant-filtered: no store CRISPRme's
        search reads holds it."""
        seq, recs, row, D, P = world_filter_pam("LowQual")
        case = Case(self.tmp, seq, recs, [row])
        self.assertClean(case)
        found = self.assertClean(case, mode="completeness")
        self.assertNotIn("limit-filtered-record", classes(found), found)  # an SNV
        tmp2 = os.path.join(self.tmp, "pass")
        os.makedirs(tmp2)
        seq, recs2, row2, _D, _P = world_filter_pam("PASS")
        case2 = Case(tmp2, seq, recs2, [row2])
        rc, found = case2.run()
        self.assertEqual(rc, 1, found)
        self.assertTrue(all(f[3] == "S02" for f in gating(found)), found)
        write_results(case.res, [set_field(row, var="%s,%s" % (tok("chr21", D),
                                                                 tok("chr21", P)))])
        rc, found = case.run()
        self.assertEqual(rc, 1, found)
        self.assertIn("variant-filtered", classes(gating(found)), found)

    def test_an_indel_the_search_skips_stays_on_its_carriers_haplotypes(self):
        """R10-F1, the phased indel store: J LowQual. CRISPRme lists S02 alone on D's
        perfect target (its store reads J on S01's slot: as the oracle's sites O and Q
        with their second deletions LowQual), and that row passes in both modes; a
        blanket drop of non-PASS records would demand S01 there. J's own targets are
        demanded of no one (no fake-indel window holds it): one limit-filtered-record
        line. Listing S01 fails; a row naming J fails variant-filtered."""
        seq, recs, row, _D, J = world_filter_second("LowQual")
        case = Case(self.tmp, seq, recs, [row])
        self.assertClean(case)
        found = self.assertClean(case, mode="completeness")
        lines = [f for f in found if f[1] == "limit-filtered-record"]
        self.assertEqual(len(lines), 1, found)
        self.assertTrue(lines[0][2].startswith("1 indel/MNV record(s)"), lines)
        self.assertFalse(any(f[3] == "S01" and f[1] in ("missing-row", "carrier-missing")
                             for f in found), found)
        write_results(case.res, [set_field(row, samples="S01,S02")])
        rc, found = case.run()
        self.assertEqual(rc, 1, found)
        self.assertTrue(all(f[3] == "S01" for f in gating(found)), found)
        write_results(case.res, [set_field(row, var=tok("chr21", J), samples="S01")])
        rc, found = case.run()
        self.assertEqual(rc, 1, found)
        self.assertIn("variant-filtered", classes(gating(found)), found)

    def test_track_b_reads_an_indel_the_search_skips(self):
        """R10-F1, Track B (review-r10's LowQual world, the oracle's Y1): a SNP row's
        carrier whose slot of the SNP also holds an insertion inside the target is dropped
        by CRISPRme whatever the insertion's FILTER (its phased indel store keeps it): the
        row listing S02 and S03 -- S01, in cis, dropped -- is nothing to report, LowQual
        or PASS."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        p = len(PRE)
        snp = (p + 6, other(GUIDE[5]), GUIDE[5], {"S01": "1|0", "S02": "1|0", "S03": "1|0"})
        start = p
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", snp), "S02,S03")
        for filt in ("LowQual", "PASS"):
            tmp = os.path.join(self.tmp, filt)
            os.makedirs(tmp)
            ins = (p + 13, seq[p + 12], seq[p + 12] + "TT", {"S01": "1|0", "S02": "0|1"},
                   filt)
            case = Case(tmp, seq, [snp, ins], [row])
            rc, found = case.run(extra=("--snp-rows", "report"))
            self.assertEqual((rc, found), (0, []), (filt, found))

    def test_a_sites_only_records_filter(self):
        """R10-F1, sites-only: the hg38_mega index's indel log holds exactly its PASS
        indels, so a sites-only indel outside PASS and '.' is in no row: a carrier-less row
        naming it fails (variant-filtered), and completeness counts it as
        limit-filtered-record, not limit-no-genotypes. The dataset's fake-indel windows
        code what CRISPRitz's overlay codes, its SNVs whose FILTER passes (the mega's indel
        rows name none other): a carrier-less row of a PASS sites-only deletion through a
        LOW_VQSLOD SNV fails variant-filtered, through the same SNV PASS it passes, and
        through an SNV of another dataset it fails variant-not-in-dataset."""
        case, info = world_deletion(self.tmp)
        seq, start, rec = info["seq"], info["start"], info["rec"]
        q = start + 2                         # a protospacer base, reference = GUIDE[2]
        y_alt = other(GUIDE[2])
        X = (start + 40, seq[start + 39:start + 41], seq[start + 39], "ExcHet")
        tgt = GUIDE[:2] + y_alt + GUIDE[3:] + "TGG"
        Y = (q + 1, seq[q], y_alt)
        mega_row = row_fields("chr21", start, "+", GUIDE + "NNN", info["ref_text"],
                              mark(GUIDE + "NNN", tgt), "TGG", 1, 0, "X",
                              "%s,%s" % (info["var"], tok("chr21", Y)), "")
        write_results(case.res, info["rows"] + [mega_row])
        so = os.path.join(self.tmp, "mega.vcf.gz")
        for filt, want in (("PASS", None), ("LOW_VQSLOD", "variant-filtered")):
            write_sites_only_vcf(so, "chr21", sorted([rec[:3], Y + (filt,), X],
                                                     key=lambda r: r[0]))
            rc, found = case.run(extra=("--vcf", so))
            if want is None:
                self.assertEqual((rc, gating(found)), (0, []), (filt, found))
            else:
                self.assertEqual(rc, 1, (filt, found))
                self.assertEqual(classes(gating(found)), {want}, (filt, found))
        # Y only in the genotyped VCF: no record of the sites-only dataset's window
        write_sites_only_vcf(so, "chr21", sorted([rec[:3], X], key=lambda r: r[0]))
        write_vcf(case.vcf, "chr21", sorted([rec, Y + ({"S05": "1|0"},)],
                                            key=lambda r: r[0]))
        rc, found = case.run(extra=("--vcf", so))
        self.assertEqual(classes(gating(found)), {"variant-not-in-dataset"}, found)
        write_vcf(case.vcf, "chr21", [rec])
        write_results(case.res, info["rows"])
        found = self.assertClean(case, mode="completeness", extra=("--vcf", so))
        lines = [f for f in found if f[1] == "limit-filtered-record"]
        self.assertEqual(len(lines), 1, found)
        self.assertFalse(any(f[1] == "limit-no-genotypes" for f in found), found)
        # a carrier-less row naming X (its strings are D's row's: X is read first)
        x_row = set_field(info["rows"][0], var=tok("chr21", X), samples="")
        write_results(case.res, info["rows"] + [x_row])
        rc, found = case.run(extra=("--vcf", so))
        self.assertEqual(rc, 1, found)
        self.assertIn("variant-filtered", classes(gating(found)), found)

    def test_a_sites_only_snv_outside_pass_and_the_cluster_chain(self):
        """R10-F1, the cluster chain. A sites-only dataset's SNP rows read its SNVs
        whatever their FILTER (hg38_mega's SNP rows name such SNVs; an index CRISPRitz's
        enricher builds does not hold them): a SNP row through one is a row CRISPRme may
        write or not, apart from the certain ones (snp_row_cps' fourth part). Its
        fake-indel windows code its PASS SNVs only: world_sites_only_chain's W, the
        dataset's worst-case row of R's window with S, is no row with S LOW_VQSLOD, so
        T1, which only W chains to B, is missed in the table holding B alone
        (missing-locus) and passes with its own row."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:] + "TGG" + POST
        p = len(PRE)
        so = os.path.join(self.tmp, "so.vcf.gz")

        class A(object):
            mm, bdna, brna, max_total_edits = 0, 0, 0, 0
        ctx = C.Ctx(A(), "NRG", 3, False, [GUIDE])
        cp = C.cluster_of(p, "+", 23, 23, False)
        for filt in ("PASS", "LOW_VQSLOD"):
            write_sites_only_vcf(so, "chr21", [(p + 6, other(GUIDE[5]), GUIDE[5], filt)])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(so)])
            rows, _capped, _drop, filt_cps, _deg, _grd = C.snp_row_cps(win, ctx, GUIDE, ctx.full(GUIDE),
                                                           0, len(seq))
            self.assertEqual((cp in rows, cp in filt_cps), (filt == "PASS",
                                                            filt != "PASS"), filt)
        g, seq, recs, rows, mega, budget, _R = world_sites_only_chain(True)
        case = Case(self.tmp, seq, recs, rows)
        run_markers(case.res, "population-level", alt=False)
        write_sites_only_vcf(so, "chr21", [mega[0], mega[1] + ("LOW_VQSLOD",)])
        rc, found = case.run(mode="completeness", guide=g, budget=budget,
                             extra=("--vcf", so))
        self.assertEqual(rc, 1, found)
        self.assertIn("missing-locus", classes(gating(found)), found)
        _g, _seq, _recs, rows2, _m, _b, _R2 = world_sites_only_chain(True, t1_row=True)
        write_results(case.res, rows2)
        rc, found = case.run(mode="completeness", guide=g, budget=budget,
                             extra=("--vcf", so))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_sites_only_window_codes_the_last_passing_record_of_a_split_site(self):
        """R10-F1, a sites-only dataset's fake-indel window (sites_only_cps): CRISPRitz's
        overlay codes the enriched genome's SNVs, the last record of a split site whose
        FILTER passes. D (a PASS deletion) has a perfect target only with E, the earlier
        of two SNV records at one flank base: with the later record L PASS, the window
        holds L, and no row of D's window is there; with L LOW_VQSLOD, E is the last record
        the window holds, and the row is."""
        G = GUIDE
        x = other(G[2])
        seq = PRE + G[:2] + x + G[3:10] + "C" + G[10:] + "TGG" + POST
        p = len(PRE)
        D = (p + 10, G[9] + "C", G[9])
        E = (p + 3, x, G[2])
        y = next(b for b in "ACGT" if b not in (x, G[2]))
        so = os.path.join(self.tmp, "so.vcf.gz")

        class A(object):
            mm, bdna, brna, max_total_edits = 0, 0, 0, 0
        ctx = C.Ctx(A(), "NRG", 3, False, [G])
        cp = C.cluster_of(p, "+", 23, 23, False)
        for filt, held in (("PASS", False), ("LOW_VQSLOD", True)):
            write_sites_only_vcf(so, "chr21", [E, (p + 3, x, y, filt), D])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(so)])
            d = win.by_key[D][0]
            got, _guarded = C.sites_only_cps(win, ctx, ctx.full(G), d, 0, len(seq))
            self.assertEqual(cp in got, held, (filt, sorted(got)))

    def test_the_density_guard_withholds_a_sites_only_datasets_dense_rows(self):
        """The sites-only density guard (2.10.0; the run's *.sites_only_density_guard.bed):
        a sites-only dataset's candidate whose window spells 2^bits sequences or more is
        withheld from its post-analysis, so its rows -- SNP and indel -- are rows CRISPRme
        does not write: open (limit-density-guard), never ones a chain counts on. A
        genotyped dataset's rows are not guarded; without the guard they are rows."""
        G = GUIDE
        x = other(G[2])
        seq = PRE + G[:2] + x + G[3:10] + "C" + G[10:] + "TGG" + POST
        p = len(PRE)
        D = (p + 10, G[9] + "C", G[9])      # a deletion: its window's target is perfect
        E = (p + 3, x, G[2])                # with E's allele at a flank base
        # two more SNVs of the dataset in the target (past the C D deletes, seq[p + 11 + k]
        # is G[10 + k]): the window spells 2^3 sequences
        F = (p + 15, G[13], other(G[13]))
        H = (p + 17, G[15], other(G[15]))
        so = os.path.join(self.tmp, "so.vcf.gz")
        write_sites_only_vcf(so, "chr21", [E, F, H, D])

        class A(object):
            mm, bdna, brna, max_total_edits = 0, 0, 0, 0
        ctx = C.Ctx(A(), "NRG", 3, False, [G])
        cp = C.cluster_of(p, "+", 23, 23, False)
        for bits, want in ((None, (True, False)), (4, (True, False)), (3, (False, True)),
                           (-1, (True, False))):
            ctx.guard_bits = bits
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(so)])
            d = win.by_key[D][0]
            got, guarded = C.sites_only_cps(win, ctx, ctx.full(G), d, 0, len(seq))
            self.assertEqual((cp in got, cp in guarded), want, bits)
            rows, opened = C.chain_cps(win, ctx, G, ctx.full(G), 0, len(seq))
            self.assertEqual(cp in rows, want[0], bits)
            self.assertEqual("guard" in opened.get(cp, set()), want[1], bits)
        # a SNP row of the dataset: a target three of whose bases its SNVs restore (2^3)
        site = list(G + "TGG")
        for i in (3, 8, 12):
            site[i] = other(site[i])
        seq2 = PRE + "".join(site) + POST
        snvs = [(p + i + 1, site[i], G[i]) for i in (3, 8, 12)]
        write_sites_only_vcf(so, "chr21", snvs)
        for bits, want in ((3, (False, True)), (4, (True, False)), (None, (True, False))):
            ctx.guard_bits = bits
            win = C.Window(seq2, "chr21", 0, len(seq2), [C.Vcf(so)])
            rows, _c, _d, _f, _g, guarded = C.snp_row_cps(win, ctx, G, ctx.full(G), 0,
                                                          len(seq2))
            self.assertEqual((cp in rows, cp in guarded), want, bits)
        # a genotyped dataset's own rows are not guarded
        gt = os.path.join(self.tmp, "gt.vcf.gz")
        write_vcf(gt, "chr21", [(q, r, a, {"S01": "1|1"}) for q, r, a in snvs])
        ctx.guard_bits = 3
        win = C.Window(seq2, "chr21", 0, len(seq2), [C.Vcf(gt)])
        rows, _c, _d, _f, _g, guarded = C.snp_row_cps(win, ctx, G, ctx.full(G), 0, len(seq2))
        self.assertIn(cp, rows)
        self.assertEqual(guarded, set())

    def test_a_run_says_its_guard_in_its_beds(self):
        """run_guard_bits: the bits a run's sites-only density guard wrote on the second line
        of its BEDs; None for a run that wrote none."""
        self.assertIsNone(C.run_guard_bits(self.tmp))
        with open(os.path.join(self.tmp, "x_chr21.snp" + C.GUARD_SUFFIX), "w") as fh:
            fh.write("## CRISPRme sites-only density guard: ...\n"
                     "## bits=20\tdataset=hg38_mega\tchrom=chr21\tkind=snp\tcandidates=4"
                     "\twithheld=1\n#chrom\tstart\tend\tguide\tcandidates_withheld\t"
                     "max_bits\n")
        self.assertEqual(C.run_guard_bits(self.tmp), 20)

    # --- R10-F2: --per-sample without the genotype store ------------------------------ #
    def test_a_per_sample_run_without_its_genotype_store_is_refused(self):
        """R10-F2: complete-search runs --per-sample population-level, and records it
        (.search_mode population-level beside --per-sample), where a genotyped dataset has
        no genotype store and no classic dictionaries (download --what index
        --no-genotypes): with a genotyped VCF given the run is refused (exit 2). Not so
        with --population-level or no mode flag (which record nothing that tells such an
        install apart), nor where every VCF is sites-only (the other case of that pair)."""
        case, info = world_deletion(self.tmp)

        def run(flag, vcf=None):
            run_markers(case.res, "population-level", alt=False)
            with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
                fh.write("input_command\t/opt/conda/bin/crisprme.py complete-search "
                         "--genome Genomes/hg38 --vcf list_vcf.txt --mm 4 --bDNA 1 "
                         "--bRNA 1 --output r %s\n" % flag)
            if vcf is not None:
                saved, case.vcf = case.vcf, vcf
            try:
                return case.run()
            finally:
                if vcf is not None:
                    case.vcf = saved
                for name in (".search_mode", ".alt_mode", ".command_line.txt"):
                    os.remove(os.path.join(case.res, name))
        rc, _found = run("--per-sample")
        self.assertEqual(rc, 2)
        self.assertTrue(C.per_sample_downgraded(case.res) is False)  # (files removed)
        for flag in ("--population-level", ""):
            rc, found = run(flag)
            self.assertEqual((rc, gating(found)), (0, []), (flag, found))
        so = os.path.join(self.tmp, "so.vcf.gz")
        write_sites_only_vcf(so, "chr21", [info["rec"][:3]])
        write_results(case.res, [set_field(info["rows"][0], samples="")])
        rc, found = run("--per-sample", vcf=so)
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_run_whose_log_says_it_had_no_genotype_store_is_refused(self):
        """R10-F2, any mode: the run's log_verbose.txt (the pipeline's stdout, in the
        result folder) records a genotyped dataset searched without its genotype store
        ('Registry-only install for chr21: ... no genotype tier to resolve carriers.', as
        new_simple_analysis.py prints it), so with --population-level or no mode flag too
        the run is refused; the same line of an aggregate registry (hg38_mega's, marked
        '(aggregate info_af registry; ...') is not one, and a folder without the log is
        checked (the PRECONDITION)."""
        case, _info = world_deletion(self.tmp)
        log = os.path.join(case.res, "log_verbose.txt")
        plain = ("Registry-only install for chr21: emitting variant off-targets with "
                 "degraded (NA) Samples -- no genotype tier to resolve carriers.\n")
        agg = ("Registry-only install for chr21: emitting variant off-targets with "
               "degraded (NA) Samples -- no genotype tier to resolve carriers (aggregate "
               "info_af registry; sample-less dict ignored).\n")
        for flag in ("--population-level", "", "--per-sample"):
            for text, refused in ((plain, True), (agg, False), (None, False)):
                if text is None:
                    if os.path.exists(log):
                        os.remove(log)
                else:
                    with open(log, "w") as fh:
                        fh.write("Analysis of chr21\n" + text)
                run_markers(case.res, "population-level" if flag != "--per-sample"
                            else "per-sample", alt=flag == "--per-sample")
                with open(os.path.join(case.res, ".command_line.txt"), "w") as fh:
                    fh.write("input_command\t/opt/conda/bin/crisprme.py complete-search "
                             "--genome Genomes/hg38 --output r %s\n" % flag)
                rc, found = case.run()
                self.assertEqual(rc, 2 if refused else 0, (flag, text, found))
        self.assertEqual(C.registry_only_logged(case.res), None)
        with open(log, "w") as fh:
            fh.write(plain + agg)
        self.assertEqual(C.registry_only_logged(case.res), ["chr21"])

    def assertClean(self, case, mode="soundness", **kw):
        rc, found = case.run(mode=mode, **kw)
        self.assertEqual((rc, gating(found)), (0, []), found)
        return found


class Driver(Base):
    def test_a_chromosome_no_vcf_covers_is_reported(self):
        """F11: an indexed VCF of chr21 and a row on chr20."""
        case, info = world_deletion(self.tmp)
        row = set_field(info["rows"][0], var=info["var"].replace("chr21", "chr20"))
        row[1] = "chr20"
        write_results(case.res, [row])
        write_fasta(os.path.join(self.tmp, "Genomes", "chr20.fa"), info["seq"], "chr20")
        bg = os.path.join(self.tmp, "idx.vcf.gz")
        write_bgzf_tabix(bg, vcf_text("chr21", [info["rec"]]))
        out = os.path.join(self.tmp, "o.tsv")
        rc = C.main(["--results", case.res, "--reference", os.path.join(self.tmp, "Genomes"),
                     "--vcf", bg, "--pam", case.pam, "--mm", "4", "--bDNA", "1", "--bRNA",
                     "1", "--mode", "soundness", "--guide", GUIDE, "--out", out])
        self.assertEqual(rc, 1)
        with open(out) as fh:
            self.assertIn("row-not-checked", fh.read())

    def test_one_fasta_of_several_contigs(self):
        """F12: chr20 before chr21 in one FASTA must not shift chr21."""
        case, info = world_deletion(self.tmp)
        fa = os.path.join(self.tmp, "both.fa")
        with open(fa, "w") as fh:
            fh.write(">chr20\n" + "ACGT" * 50 + "\n>chr21\n" + info["seq"] + "\n>chrX\nAAAA\n")
        out = os.path.join(self.tmp, "o.tsv")
        rc = C.main(["--results", case.res, "--reference", fa, "--vcf", case.vcf, "--pam",
                     case.pam, "--mm", "4", "--bDNA", "1", "--bRNA", "1", "--mode",
                     "soundness", "--chrom", "chr21", "--guide", GUIDE, "--out", out])
        with open(out) as fh:
            self.assertEqual(rc, 0, fh.read())

    def test_the_bgzf_tabix_reader_matches_the_stream(self):
        """N4: a BGZF file in small blocks with a tabix index: region reads seek and find
        what a full stream finds, and the checker catches a mutation through it."""
        case, info = world_deletion(self.tmp)
        recs = [info["rec"]] + [(p, info["seq"][p - 1], other(info["seq"][p - 1]),
                                 {"S0%d" % (1 + p % 7): "1|0"})
                                for p in range(3, len(info["seq"]) - 2, 4)]
        recs.sort(key=lambda r: r[0])
        bg = os.path.join(self.tmp, "idx.vcf.gz")
        write_bgzf_tabix(bg, vcf_text("chr21", recs))
        plain = os.path.join(self.tmp, "plain.vcf.gz")
        write_vcf(plain, "chr21", recs)
        vi, vp = C.Vcf(bg), C.Vcf(plain)
        self.assertIsNotNone(vi.tbi)
        for lo, hi in ((0, 10), (20, 40), (35, 36), (0, 1000)):
            a = [(r["pos1"], r["ref"], r["alts"]) for r in vi.records("chr21", lo, hi)]
            b = [(r["pos1"], r["ref"], r["alts"]) for r in vp.records("chr21", lo, hi)]
            self.assertEqual(a, b, (lo, hi))
        # the in-memory block index of an unindexed BGZF file gives the same
        os.remove(bg + ".tbi")
        vb = C.Vcf(bg)
        self.assertIsNone(vb.tbi)
        for lo, hi in ((0, 10), (20, 40), (0, 1000)):
            a = [(r["pos1"], r["ref"]) for r in vb.records("chr21", lo, hi)]
            b = [(r["pos1"], r["ref"]) for r in vp.records("chr21", lo, hi)]
            self.assertEqual(a, b, (lo, hi))
        write_bgzf_tabix(bg, vcf_text("chr21", recs))
        write_results(case.res, [set_field(info["rows"][0], samples="S01")])
        out = os.path.join(self.tmp, "o.tsv")
        rc = C.main(["--results", case.res, "--reference", os.path.join(self.tmp, "Genomes"),
                     "--vcf", bg, "--pam", case.pam, "--mm", "4", "--bDNA", "1", "--bRNA",
                     "1", "--mode", "soundness", "--chrom", "chr21", "--guide", GUIDE,
                     "--out", out])
        self.assertEqual(rc, 1)

    def test_a_region_read_yields_each_record_once(self):
        """A long deletion crossing a 16-kb tile edge lives in a coarser tabix bin, so a
        query of the next tile gets its one-line chunk, then (past the records of the tile
        before, which the query does not ask for) the tile's own chunk. Reading a chunk
        past its end gave every record of the window twice, and a record twice on one slot
        reads as overlapping alleles: each chunk is read to its end, each line once."""
        recs = [(p, "A", "C", {"S01": "0|1"}) for p in range(16001, 16381, 10)]
        recs += [(16300, "A" + "C" * 200, "A", {"S02": "1|0"})]
        recs += [(p, "A", "G", {"S03": "1|0"}) for p in range(16401, 17001, 10)]
        recs.sort(key=lambda r: r[0])
        bg = os.path.join(self.tmp, "long.vcf.gz")
        write_bgzf_tabix(bg, vcf_text("chr21", recs), block=120)
        plain = os.path.join(self.tmp, "long_plain.vcf.gz")
        write_vcf(plain, "chr21", recs)
        vi, vp = C.Vcf(bg), C.Vcf(plain)
        self.assertGreater(len(vi.tbi.chunks("chr21", 16450, 16600)), 1)
        for lo, hi in ((16450, 16600), (16390, 16700), (16000, 17000), (16100, 16200)):
            a = [(r["pos1"], r["ref"]) for r in vi.records("chr21", lo, hi)]
            b = [(r["pos1"], r["ref"]) for r in vp.records("chr21", lo, hi)]
            self.assertEqual(a, b, (lo, hi))

    def test_a_known_finding_turns_non_gating(self):
        case, info = world_deletion(self.tmp)
        write_results(case.res, [])
        known = os.path.join(self.tmp, "known.tsv")
        with open(known, "w") as fh:
            fh.write("class\tchrom\tstart\tstrand\tvariants\tsamples\treason\n")
            for s in ("S01", "S02"):
                fh.write("missing-row\tchr21\t%d\t+\t%s\t%s\tan engine loss under study\n"
                         % (info["start"], info["var"], s))
            fh.write("missing-row\tchr21\t5\t+\tchr21_6_A_AT\tS09\tfixed since\n")
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"),
                             extra=("--known", known))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("known-defect", classes(found))
        self.assertIn("known-not-seen", classes(found))

    def test_jobs_give_the_same_findings(self):
        case, info = world_deletion(self.tmp)
        write_results(case.res, [set_field(info["rows"][0], samples="S01")])
        r1 = case.run(mode="both", budget=("0", "0", "0", "0"))
        r4 = case.run(mode="both", budget=("0", "0", "0", "0"),
                      extra=("--jobs", "3", "--tile", "20"))
        self.assertEqual(r1[0], r4[0])
        self.assertEqual(sorted(r1[1]), sorted(r4[1]))
        # the findings written as units end are replaced by the sorted TSV
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "fail.tsv")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "fail.tsv.partial")))


# --------------------------------------------------------------------------- #
# the checker's own primitives
# --------------------------------------------------------------------------- #
class Pieces(unittest.TestCase):
    def test_places_of_an_insertion(self):
        seq = "ACGT" * 10
        hap = C.make_hap(seq, [al(5, "A", "ATTT")], 0, 12)
        self.assertEqual(hap.seq, seq[:5] + "TTT" + seq[5:12])
        self.assertEqual(hap.place[4], (4, 0))
        self.assertEqual(hap.place[5:8], [(4, 1), (4, 2), (4, 3)])
        self.assertEqual(hap.place[8], (5, 0))
        self.assertEqual(hap.clash, [])

    def test_places_of_a_deletion_an_mnv_a_shared_tail_and_a_complex_allele(self):
        seq = "ACGT" * 10
        hap = C.make_hap(seq, [al(5, "ACG", "A")], 0, 12)
        self.assertEqual(hap.seq, seq[:5] + seq[7:12])
        self.assertEqual(hap.place[4:6], [(4, 0), (7, 0)])
        hap = C.make_hap(seq, [al(5, "AC", "TG")], 0, 12)
        self.assertEqual((hap.seq[4:6], hap.place[4:6]), ("TG", [(4, 0), (5, 0)]))
        # GTA>GA at POS 3: the shared A keeps its own coordinate
        hap = C.make_hap(seq, [al(3, "GTA", "GA")], 0, 8)
        self.assertEqual(hap.seq, "ACGACGT")
        self.assertEqual(hap.place[2:4], [(2, 0), (4, 0)])
        # A>TCG: T in the A's place, C and G added after it (the last replaced base)
        hap = C.make_hap(seq, [al(5, "A", "TCG")], 0, 8)
        self.assertEqual(hap.seq[4:7], "TCG")
        self.assertEqual(hap.place[4:7], [(4, 0), (4, 1), (4, 2)])

    def test_an_alleles_window_is_its_bases_placed_inside(self):
        """make_hap writes an allele's edit by slices (HPRC holds 86-kb insertions): the
        slice for a window is exactly the bases whose place lies in it."""
        import random
        rnd = random.Random(7)
        ref = "".join(rnd.choice("ACGT") for _ in range(400))
        for _ in range(3000):
            p = rnd.randint(50, 300)
            r = ref[p - 1:p - 1 + rnd.randint(1, 30)]
            alt = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(1, 30)))
            if rnd.random() < 0.5:
                alt = r[0] + alt
            if alt == r:
                continue
            a = al(p, r, alt)
            lo = rnd.randint(0, 380)
            hi = rnd.randint(lo, 400)
            i, j = a.window(lo, hi)
            want = [k for k, pl in enumerate(a.places()) if lo <= pl[0] < hi]
            self.assertEqual(list(range(i, j)), want, (a.key, lo, hi))
            hap = C.make_hap(ref, [a], lo, hi)
            self.assertEqual(hap.place, sorted(hap.place))  # bisection relies on it

    def test_alleles_whose_ref_spans_meet_are_a_clash(self):
        """VCF 4.x: two alleles of one haplotype describing one base give no sequence."""
        seq = "ACGT" * 10
        for other in (al(6, "CG", "C"), al(5, "A", "G"), al(6, "C", "T"), al(5, "A", "AT")):
            hap = C.make_hap(seq, [al(5, "ACG", "A"), other], 0, 12)
            self.assertEqual(len(hap.clash), 1, other.key)
            self.assertEqual(hap.seq, seq[:12])  # neither is written
        hap = C.make_hap(seq, [al(5, "A", "AT"), al(6, "C", "T")], 0, 12)
        self.assertEqual(hap.clash, [])  # side by side, not overlapping

    def test_the_junction_rule(self):
        seq = "ACGTACGTACGTACGTACGT"
        d = al(6, "CG", "C")     # deletes the G at 0-based 6
        hap = C.make_hap(seq, [d], 0, 20)
        i = hap.at(5)[0]
        self.assertTrue(C.holds_edit(hap, i - 2, 4, d))       # holds 5 and 7 side by side
        self.assertFalse(C.holds_edit(hap, i + 1, 4, d))      # 7.. only
        self.assertFalse(C.holds_edit(hap, i - 3, 4, d))      # ..5 only
        ins = al(6, "C", "CAAA")
        hap = C.make_hap(seq, [ins], 0, 20)
        j = hap.at(5)[1]                                      # the first added base
        self.assertTrue(C.holds_edit(hap, j + 2, 3, ins))
        self.assertFalse(C.holds_edit(hap, j - 3, 3, ins))
        m = al(6, "CG", "TT")
        hap = C.make_hap(seq, [m], 0, 20)
        self.assertTrue(C.holds_edit(hap, 4, 2, m))
        self.assertFalse(C.holds_edit(hap, 7, 3, m))

    def test_touches(self):
        self.assertTrue(C.touches(al(100, "GT", "G"), 99, 101))
        self.assertFalse(C.touches(al(100, "GT", "G"), 101, 120))
        self.assertTrue(C.touches(al(100, "G", "GAAA"), 98, 100))    # adds between 99, 100
        self.assertFalse(C.touches(al(100, "G", "GAAA"), 100, 120))  # adds before 100
        self.assertFalse(C.touches(al(100, "G", "GAAA"), 90, 99))    # adds after 99
        self.assertTrue(C.touches(al(100, "A", "C"), 99, 99))

    def test_genotypes(self):
        self.assertEqual(C.parse_gt("0|1"), ((0, 1), True))
        self.assertEqual(C.parse_gt("1/0"), ((1, 0), False))
        self.assertEqual(C.parse_gt("1|."), ((1, None), True))
        self.assertEqual(C.parse_gt("1"), ((1,), True))
        self.assertEqual(C.record_slots(((0, 1), True), 1), ([1], True))
        self.assertEqual(C.record_slots(((0, 1), False), 1), ([None], False))
        self.assertEqual(C.record_slots(((1, 1), False), 1), ([0, 1], True))
        main = ((1, 0), True)
        self.assertEqual(C.on_slot(main, 0, ((1, 0), True), 1), "yes")
        self.assertEqual(C.on_slot(main, 0, ((0, 1), True), 1), "no")
        self.assertEqual(C.on_slot(main, 0, ((0, 1), False), 1), "?")
        self.assertEqual(C.on_slot(main, 0, ((None, 1), True), 1), "?")
        self.assertEqual(C.on_slot(main, 0, ((1, 1), False), 1), "yes")
        self.assertEqual(C.on_slot(main, 0, ((1,), True), 1), "?")    # haploid vs diploid

    def test_cfd(self):
        """The CFD tables read as data reproduce scores CRISPRme printed: a perfect NGG
        site, one mismatch, a DNA and an RNA bulge."""
        cfd = C.Cfd(C.default_cfd_dir())
        self.assertAlmostEqual(cfd.score(GUIDE + "NNN", GUIDE + "TGG"), 1.0)
        self.assertAlmostEqual(cfd.score(GUIDE + "NNN", GUIDE + "TAG"), 0.259259259, 6)
        # the oracle's alternative alignment chr21:8000 (CFD 0.455, printed)
        self.assertAlmostEqual(cfd.score("ACTGAAATCTGTAAGCAGGCNNN",
                                         "ACTGACATCTGTAAGCAGGCTGG"), 0.455, 3)
        self.assertAlmostEqual(cfd.score("A-CTGAAATCTGTAAGCAGGCNNN",
                                         "AC-TGACATCTGTAAGCAGGCTGG"), 0.423, 3)

    def test_cluster_position(self):
        # + strand, 3' PAM: an RNA bulge (one base fewer) or a DNA bulge moves the start
        self.assertEqual(C.cluster_of(1001, "+", 22, 23, False), 1000)
        self.assertEqual(C.cluster_of(999, "+", 24, 23, False), 1000)
        self.assertEqual(C.cluster_of(1001, "-", 22, 23, False), 1001)
        self.assertEqual(C.cluster_of(1001, "-", 22, 24, True), 999)

    def test_the_aligner_finds_a_perfect_target_and_its_bulges(self):
        full = GUIDE + "NRG"
        seq = "TTTT" + GUIDE + "TGG" + "TTTT"
        hits = list(C.scan_sequence(full, seq, 4, 1, 1, 6, 3, False))
        perfect = [h for h in hits if h[4:] == (0, 0, 0)]
        self.assertEqual([(h[0], h[1], h[3]) for h in perfect], [("+", 4, GUIDE + "TGG")])
        for _st, _off, sa, ta, _m, _d, _r in hits:
            self.assertNotIn("-", sa[-3:])
            self.assertNotIn("-", ta[-3:])

    def test_the_aligner_reports_a_dna_bulge_next_to_the_pam(self):
        full = GUIDE + "NRG"
        seq = "TTTT" + GUIDE + "A" + "TGG" + "TTTT"
        hits = [h for h in C.scan_sequence(full, seq, 4, 1, 1, 6, 3, False)
                if h[2].endswith("-NRG")]
        self.assertTrue(hits)
        self.assertEqual(hits[0][5], 1)

    def test_recount(self):
        self.assertEqual(C.recount("CTAACAGTTGCTTTTATCAC-NNN",
                                   "-TAACtGTTtCTgTTATCAgTAAG", 3, False), (4, 1, 1, None))
        self.assertIsNotNone(C.recount(GUIDE + "N-N", GUIDE + "TGG", 3, False)[3])
        self.assertIsNotNone(C.recount("-" + GUIDE[1:] + "NNN", "A" + GUIDE[1:] + "TGG",
                                       3, False)[3])

    def test_the_pam_file_and_the_guides(self):
        tmp = tempfile.mkdtemp(prefix="chkidl_")
        try:
            p3 = os.path.join(tmp, "nrg.txt")
            with open(p3, "w") as fh:
                fh.write("N" * 20 + "NRG 3\n")
            self.assertEqual(C.read_pam(p3), ("NRG", 3, False))
            p5 = os.path.join(tmp, "ttt.txt")
            with open(p5, "w") as fh:
                fh.write("TTTV" + "N" * 20 + " -4\n")
            self.assertEqual(C.read_pam(p5), ("TTTV", 4, True))
            g = os.path.join(tmp, "guides.txt")
            with open(g, "w") as fh:
                fh.write(GUIDE + "NNN\n\nTTTTAAAATTTTAAAATTTTNNN\n")
            self.assertEqual(C.read_guides([g], 3, False), [GUIDE, "TTTTAAAATTTTAAAATTTT"])
            self.assertEqual(C.read_guides(["NNNN" + GUIDE], 4, True), [GUIDE])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_a_split_site_shadows_its_earlier_records(self):
        tmp = tempfile.mkdtemp(prefix="chkidl_")
        try:
            seq = "ACGT" * 20
            vcf = os.path.join(tmp, "v.vcf.gz")
            write_vcf(vcf, "chr21", [(5, "A", "C", {"S01": "1|0"}),
                                     (5, "A", "G", {"S02": "0|1"}),
                                     (9, "A", "C,G", {"S01": "1|2"})])
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
            self.assertEqual({a.key for a in win.shadowed}, {(5, "A", "C")})
            # every carrier of the two alts of POS 9 holds both: the genotypes cannot
            # tell CRISPRme which allele number is which
            self.assertEqual({a.key for a in win.untied}, {(9, "A", "C"), (9, "A", "G")})
            g = win.by_key[(9, "A", "G")][0]
            self.assertEqual(win.call(g, "S01"), ((1, 2), True))
            self.assertEqual(g.number, 2)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# a run of CRISPRme 2.10.0: what its own records say of how it read its data
# --------------------------------------------------------------------------- #
G21 = "GATTACAGTCCATGCATGGCT"   # a 21-nt guide (SaCas9: NNGRRT, a 27-nt target)
SA_PAM = "N" * 21 + "NNGRRT 6"


def write_bin(path, fmt, flags=0, guide_len=20, pam="NRG", budget=2):
    """A CRISPRitz .bin's opening bytes, as its index format writes them (CRISPRitz's
    docs/INDEX_FORMAT.md): format 2 its magic and header (flags: 4 = a target holds Q),
    format 1 its leaf count and guide length; a few body bytes after."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        if fmt == 2:
            fh.write(b"CZI\xfe" + struct.pack("=6i", 2, 28 + len(pam), flags, guide_len,
                                              len(pam), budget) + pam.encode())
        else:
            fh.write(struct.pack("=ii", 100, guide_len))
        fh.write(b"\0" * 64)


def write_indel_log(path, vcf, window):
    """The indel log CRISPRitz's add-variants writes of ``vcf`` (log<chrom>.txt of
    Dictionaries/log_indels_<dataset>), as far as the checker reads it: one row per
    indel-path ALT of a PASS or '.' record (an indel or an MNV), its CHR field holding its
    fake-indel window's start and end, ``window`` = (before, after) reference bases around
    the record's REF (POS - 1 - before, POS - 1 + len(REF) + after)."""
    before, after = window
    rows, k = [], 0
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[6] not in ("PASS", "."):
                continue
            pos, ref = int(f[1]), f[3]
            for alt in f[4].split(","):
                if alt.startswith("<") or alt == "*" or not (
                        len(alt) != len(ref) or (len(alt) > 1 and len(ref) > 1)):
                    continue
                k += 1
                rows.append("%s_%d-%d_%d\tS01\t.\t0.1\t%s_%d_%s_%s\t0,1\tN" % (
                    f[0], pos - 1 - before, pos - 1 + len(ref) + after, k, f[0], pos, ref,
                    alt))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n" +
                 "".join(r + "\n" for r in rows))


def as_2100(case, log="", var_flags=4, var_format=2, indels_format=2, enriched=None,
            budget=2, guide_len=20, pam="NRG", window=(25, 27)):
    """Lay ``case`` out as a run of CRISPRme 2.10.0 of dataset v: its results folder moved
    to work/Results/run, with .Params.txt (Genome_idx), .indexes_searched.txt and (``log``
    not None) log_verbose.txt; the variant index and the INDELS index it names in
    work/genome_library, a .bin of the format given each (``var_format`` None: no variant
    index, as for a search without bulges); ``enriched``: the variant genome, written as
    work/Genomes/hg38+v/chr21.enriched.fa; ``window`` (None: none): the (before, after)
    of the fake-indel windows the dataset's indel log, work/Dictionaries/log_indels_v,
    records (write_indel_log): 25/27, an enrichment of CRISPRitz 2.9.1 or earlier, unless
    a test says otherwise (30/30: 2.10.0's)."""
    work = os.path.join(case.tmp, "work")
    lib = os.path.join(work, "genome_library")
    res = os.path.join(work, "Results", "run")
    os.makedirs(os.path.dirname(res), exist_ok=True)
    shutil.move(case.res, res)
    case.res = res
    var = "%d-%s_%d_hg38+v" % (guide_len, pam, budget)
    if var_format is not None:
        write_bin(os.path.join(lib, var, "%s_chr21_1.bin" % pam), var_format, var_flags,
                  guide_len, pam, budget)
    write_bin(os.path.join(lib, var + "_INDELS", "%s_%d_fakechr21" % (pam, budget),
                           "%s_fakechr21_1.bin" % pam), indels_format, 0, guide_len, pam,
              budget)
    with open(os.path.join(res, ".Params.txt"), "w") as fh:
        fh.write("Genome_ref\thg38\nGenome_idx\t%s_%d_hg38+v\n" % (pam, budget))
    with open(os.path.join(res, ".indexes_searched.txt"), "w") as fh:
        fh.write("Reference_index\t%d-%s_%d_hg38 (in %s): built for the PAM %s\n"
                 % (guide_len, pam, budget, lib, pam))
        fh.write("Variant_index:v\t%s (in %s): built for the PAM %s on the 3' side of a "
                 "%d-nt guide, with up to %d bulge(s) of each type; its indels: %s_INDELS\n"
                 % (var, lib, pam, guide_len, budget - 1, var))
    if log is not None:
        with open(os.path.join(res, "log_verbose.txt"), "w") as fh:
            fh.write(log)
    if enriched is not None:
        os.makedirs(os.path.join(work, "Genomes", "hg38+v"), exist_ok=True)
        write_fasta(os.path.join(work, "Genomes", "hg38+v", "chr21.enriched.fa"), enriched)
    if window is not None:
        write_indel_log(os.path.join(work, "Dictionaries", "log_indels_v",
                                     "log%s.txt" % case.chrom), case.vcf, window)
    return case


def world_window_edge(tmp, bulged=False):
    """R10-N1 in miniature: an insertion of a T after the anchor A makes a SaCas9 target
    (G21 + ACGAA + the T, PAM NNGRRT) that starts at its fake-indel window's first base, 25
    bases before the anchor (the reference reads ACGAAC there: no PAM). S01 carries it.
    ``bulged``: the target holds one base more than the guide (a DNA bulge), and starts a
    base before the window."""
    head = G21[:10] + ("C" if bulged else "") + G21[10:]
    left = head + "ACGAA"
    seq = PRE + left + "CTTGGACCATGG" + POST
    anchor = len(PRE) + len(left)          # VCF POS of the anchor A
    rec = (anchor, "A", "AT", {"S01": "1|0"})
    return seq, [rec], len(PRE)


def world_site_edge(kind):
    """R2B-F2's worlds: a SaCas9 target (G21, NNGRRT) of S01's record, whose fake-indel
    window starts 25 bases before the record's first REF base, at PRE's end. CRISPRitz
    keeps a site by its PAM: its guide region, the guide's 21 bases from the PAM, must lie
    in the window, and on an index of format 1 its pad too, the budget's bases past it.
      'rna'  a deletion AC>A; the target, one RNA bulge (26 bases), starts 2 bases in, so
             its site's guide region starts 1 base in: format 1's 2-base pad reaches the
             separator N, format 2 holds it;
      'dna'  a deletion GC>G; the target, one DNA bulge (28 bases), starts 1 base in, its
             site's guide region 2 bases in: both formats hold it;
      'mnv'  an MNV CA>TG whose first ALT base ends the PAM; the target, one RNA bulge,
             starts at the window's first base, so its site's guide region starts on the
             separator N: no index holds it.
    Returns (reference, records, budget of the search)."""
    if kind == "rna":
        seq = PRE + "TT" + G21[:10] + G21[11:] + "CCGACAT" + POST + LONG
        rec = (len(PRE) + 26, "AC", "A", {"S01": "1|0"})
        budget = ("0", "0", "1", "1")
    elif kind == "dna":
        seq = PRE + "T" + G21[:10] + "T" + G21[10:] + "CCGCAAT" + POST + LONG
        rec = (len(PRE) + 26, "GC", "G", {"S01": "1|0"})
        budget = ("0", "1", "0", "1")
    else:
        seq = PRE + G21[:10] + G21[11:] + "CCGAACA" + POST + LONG
        rec = (len(PRE) + 26, "CA", "TG", {"S01": "1|0"})
        budget = ("0", "0", "1", "1")
    assert seq[rec[0] - 1:rec[0] - 1 + len(rec[1])] == rec[1]
    return seq, [rec], budget


class ReleaseTwoTen(Base):
    """2.10.0: the run's records say how its CRISPRme and CRISPRitz read the data."""

    def test_an_opt_out_run_is_refused(self):
        """R11-F1: CRISPRME_INDEL_SNP=0 (its search notes it in log_verbose.txt): every SNP
        code reads as the reference, an indel row lists every carrier of its indel, and the
        per-carrier resolution this checker models does not apply: refused (exit 2), in
        every mode; the same run without the note is checked."""
        case, _ = world_indel_snp(self.tmp)
        as_2100(case, log="STARTING ANALYSIS FOR v\nNOTE: the SNP+indel co-occurrence is off "
                          "(CRISPRME_INDEL_SNP=0): every SNP code of the fake-indel genome of "
                          "v reads as the reference base, so ...\n")
        for mode in ("soundness", "completeness", "both"):
            rc, _found = case.run(mode=mode)
            self.assertEqual(rc, 2, mode)
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("STARTING ANALYSIS FOR v\n")
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_store_of_2_10_0_holds_no_indel_the_search_skips(self):
        """The store's FILTER rule (2.10.0): J, a LowQual deletion inside D's target on
        S01's D slot, is in no store, so S01 is listed for D's perfect target beside S02. A
        row listing S02 alone fails (S01 missing); the row listing both passes. With the
        log's note of a store written before 2.10.0 (which keeps J), S02 alone passes, as
        before; a run of 2.10.0 with no log_verbose.txt cannot say: refused."""
        seq, recs, row, _D, J = world_filter_second("LowQual")
        case = as_2100(Case(self.tmp, seq, recs, [row]))
        rc, found = case.run()
        self.assertEqual(rc, 1, found)
        self.assertTrue(gating(found) and all(f[3] == "S01" for f in gating(found)), found)
        write_results(case.res, [set_field(row, samples="S01,S02")])
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        rc, found = case.run(mode="completeness")
        self.assertEqual((rc, gating(found)), (0, []), found)
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("NOTE [indel-snp]: the phased indel genotype store of v on chr21 was "
                     "written before CRISPRme 2.10.0: it holds ...\n")
        write_results(case.res, [row])
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        os.remove(os.path.join(case.res, "log_verbose.txt"))
        rc, _found = case.run()
        self.assertEqual(rc, 2)

    def test_the_filtered_record_line_says_what_the_store_did(self):
        """The limit-filtered-record line says what a genotyped dataset's indel outside PASS
        and '.' was to the run: on its carriers' haplotypes where the run read a store of
        every FILTER (an earlier CRISPRme's), the reference where it read one of 2.10.0."""
        seq, recs, row, _D, _J = world_filter_second("LowQual")
        os.makedirs(os.path.join(self.tmp, "old"))
        case = Case(os.path.join(self.tmp, "old"), seq, recs, [row])
        rc, found = case.run(mode="completeness")
        lines = [f[2] for f in found if f[1] == "limit-filtered-record"]
        self.assertEqual(rc, 0, found)
        self.assertEqual(len(lines), 1, found)
        self.assertIn("stay on their carriers' haplotypes", lines[0])
        os.makedirs(os.path.join(self.tmp, "new"))
        case = as_2100(Case(os.path.join(self.tmp, "new"), seq, recs,
                            [set_field(row, samples="S01,S02")]))
        rc, found = case.run(mode="completeness")
        lines = [f[2] for f in found if f[1] == "limit-filtered-record"]
        self.assertEqual(rc, 0, found)
        self.assertEqual(len(lines), 1, found)
        self.assertIn("are read as the reference too", lines[0])
        self.assertNotIn("stay on their carriers' haplotypes", lines[0])

    def test_every_allele_of_a_split_site_is_read(self):
        """CRISPRitz 2.10.0 codes every allele of a split site: S01's target through the
        earlier record's ALT (the guide's base) is demanded -- a missing row fails --,
        where on a genome of the last record (2.9.1) it is limit-last-alt-wins. The variant
        index's Q flag says so; failing that, its enriched genome's code at the site."""
        seq = PRE + GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1, snp = len(PRE) + 10, len(PRE) + 5
        ref_b = seq[snp - 1]
        third = next(x for x in "ACGT" if x not in (ref_b, GUIDE[4]))
        recs = [(pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"}),
                (snp, ref_b, GUIDE[4], {"S01": "1|0"}),
                (snp, ref_b, third, {"S02": "1|0"})]
        budget = ("0", "0", "0", "0")
        old = Case(os.path.join(self.tmp, "old") if os.makedirs(
            os.path.join(self.tmp, "old")) is None else None, seq, recs, [])
        rc, found = old.run(mode="completeness", budget=budget)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-last-alt-wins", classes(found), found)
        code = {frozenset("AC"): "M", frozenset("AG"): "R", frozenset("AT"): "W",
                frozenset("CG"): "S", frozenset("CT"): "Y", frozenset("GT"): "K",
                frozenset("ACG"): "V", frozenset("ACT"): "H", frozenset("AGT"): "D",
                frozenset("CGT"): "B"}
        for tag, kw in (("q", dict(var_flags=4)),
                        ("fa", dict(var_flags=0, enriched=seq[:snp - 1] + code[frozenset(
                            (ref_b, GUIDE[4], third))] + seq[snp:]))):
            tmp = os.path.join(self.tmp, tag)
            os.makedirs(tmp)
            case = as_2100(Case(tmp, seq, recs, []), **kw)
            rc, found = case.run(mode="completeness", budget=budget)
            self.assertEqual(rc, 1, (tag, found))
            self.assertIn("missing-row", classes(gating(found)), (tag, found))
            self.assertNotIn("limit-last-alt-wins", classes(found), (tag, found))
        # the enriched genome codes the last record alone (and the index flags no Q): a
        # genome of 2.9.1's reading, the limit again
        tmp = os.path.join(self.tmp, "last")
        os.makedirs(tmp)
        case = as_2100(Case(tmp, seq, recs, []), var_flags=0, enriched=seq[:snp - 1]
                       + code[frozenset((ref_b, third))] + seq[snp:])
        rc, found = case.run(mode="completeness", budget=budget)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-last-alt-wins", classes(found), found)
        # neither tells (no Q, no enriched genome, no four-base position): refused
        tmp = os.path.join(self.tmp, "none")
        os.makedirs(tmp)
        case = as_2100(Case(tmp, seq, recs, []), var_flags=0)
        rc, _found = case.run(mode="completeness", budget=budget)
        self.assertEqual(rc, 2)

    def test_the_reading_taken_is_printed(self):
        """The check prints how it read the run: of a run of 2.10.0 whose VCF has no split
        SNP site, that nothing told the allele coding apart (it changes nothing there), not
        the last record's reading; of one with a split site and a variant index that flags
        Q, every allele; of a run without .indexes_searched.txt, an earlier CRISPRme's; and
        where soundness has no row to check, that there was nothing to read it for."""
        def said(case, **kw):
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                case.run(**kw)
            return [x for x in err.getvalue().splitlines() if "read as a run of" in x]

        plain = os.path.join(self.tmp, "plain")
        os.makedirs(plain)
        case, _ = world_indel_snp(plain)
        self.assertEqual(said(case), ["check_indel_rows: read as a run of CRISPRme 2.8.1 or "
                                      "earlier (every dataset: the last record of a split "
                                      "site, store of every FILTER, INDELS format 1, "
                                      "fake-indel windows 25/27)"])
        as_2100(case, var_flags=0, window=(30, 30))
        self.assertEqual(said(case), ["check_indel_rows: read as a run of CRISPRme 2.10.0 or "
                                      "later (v on chr21: no split SNP site to tell, store "
                                      "of PASS and '.', INDELS format 2, fake-indel windows "
                                      "30/30)"])
        seq = PRE + GUIDE[:4] + other(GUIDE[4]) + GUIDE[5:10] + "C" + GUIDE[10:] + "TGG" + POST
        snp, ref_b = len(PRE) + 5, seq[len(PRE) + 4]
        third = next(x for x in "ACGT" if x not in (ref_b, GUIDE[4]))
        recs = [(len(PRE) + 10, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0"}),
                (snp, ref_b, GUIDE[4], {"S01": "1|0"}), (snp, ref_b, third, {"S02": "1|0"})]
        os.makedirs(os.path.join(self.tmp, "split"))
        case = as_2100(Case(os.path.join(self.tmp, "split"), seq, recs, []))
        self.assertEqual(said(case), ["check_indel_rows: read as a run of CRISPRme 2.10.0 or "
                                      "later (no chromosome to read it for)"])
        self.assertEqual(said(case, mode="completeness", budget=("0", "0", "0", "0")),
                         ["check_indel_rows: read as a run of CRISPRme 2.10.0 or later (v on "
                          "chr21: every allele, store of PASS and '.', INDELS format 2, "
                          "fake-indel windows 25/27)"])

    def printed_reading(self, window):
        """The stderr of a check of a run of 2.10.0 (dataset v on chr21; its variant index
        codes every allele, its store 2.10.0's), whose indel log records fake-indel windows
        of ``window``."""
        tmp = os.path.join(self.tmp, "w%d" % window[0])
        os.makedirs(tmp)
        case, _ = world_indel_snp(tmp)
        as_2100(case, window=window)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            case.run()
        return err.getvalue()

    def test_the_documented_reading_is_the_printed_one(self):
        """W-R1-F3: the CHANGELOG and RELEASE_CHECKLIST (gate 4c) quote the reading the check
        prints (on stderr) of a run on an index this release built as the check prints it,
        its fake-indel windows 30/30 included; gate 4c's check of that reading passes such a
        run and fails one of other windows, or of windows its records do not tell."""
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        err = self.printed_reading((30, 30))
        said = [x for x in err.splitlines() if "read as a run of" in x]
        self.assertEqual(len(said), 1, err)
        printed = said[0].split("(v on chr21: ", 1)[1][:-1]
        self.assertTrue(printed.endswith("INDELS format 2, fake-indel windows 30/30"), printed)
        with open(os.path.join(repo, "docs", "RELEASE_CHECKLIST.md")) as fh:
            checklist = fh.read()
        quote = " ".join(re.sub(r"\n#\s+", " ", checklist).split())
        m = re.search(r'"read as a run of CRISPRme 2\.10\.0 or later \(\.\.\. ([^"]*)\)"', quote)
        self.assertIsNotNone(m, "RELEASE_CHECKLIST 4c quotes no reading")
        self.assertEqual(m.group(1), printed)
        with open(os.path.join(repo, "CHANGELOG.md")) as fh:
            changelog = " ".join(fh.read().split())
        m = re.search(r"`read as a run of CRISPRme 2\.10\.0 or later \(<dataset> on <chrom>: "
                      r"([^`]*)\)`", changelog)
        self.assertIsNotNone(m, "the CHANGELOG quotes no reading")
        self.assertEqual(m.group(1), printed)
        checks = [x for x in checklist.splitlines() if 'grep -o "fake-indel windows' in x]
        self.assertEqual(len(checks), 2, checks)   # (per chromosome, and genome-wide)
        gate = os.path.join(self.tmp, "gate")
        os.makedirs(gate)
        told = self.printed_reading((25, 27))
        for text, passes in ((err, True), (told, False),
                             (err.replace("windows 30/30", "windows not told"), False)):
            for name in ("indelrows_per-sample_soundness.err", "indelrows_gw_readme.err"):
                with open(os.path.join(gate, name), "w") as fh:
                    fh.write(text)
            for line in checks:
                out = subprocess.run(["bash", "-c", line], capture_output=True, text=True,
                                     env={"GATE": gate, "mode": "per-sample",
                                          "check": "soundness", "g": "readme:guide.txt",
                                          "PATH": os.environ.get("PATH", "/usr/bin:/bin")})
                self.assertEqual(out.stdout.strip(), "" if passes else "GATE_FAILED",
                                 (line, text[-300:]))

    def test_a_chromosome_without_an_indel_reads_no_window(self):
        """W-R2-N2: a chromosome whose VCF holds no indel record has no fake-indel window (its
        indel log no row, its INDELS index none recorded): its reading says so, "fake-indel
        windows none: no indel record", where it said "not told" as of a run whose records do
        not tell them, and the run is checked; gate 4c passes a genome-wide check that reads it
        beside chromosomes of 30/30 (it failed every such check of a panel with a chromosome
        free of indels), and still fails one of windows not told."""
        tmp = os.path.join(self.tmp, "snp_only")
        os.makedirs(tmp)
        seq = PRE + GUIDE[:5] + "C" + GUIDE[6:] + "TGG" + POST
        case = as_2100(Case(tmp, seq, [(len(PRE) + 6, "C", GUIDE[5], {"S05": "1|0"})], []),
                       window=(30, 30))
        err = io.StringIO()
        with contextlib.redirect_stderr(err):   # (checked: every chromosome the VCF holds)
            rc, found = case.run("completeness")
        none = err.getvalue()
        self.assertEqual((rc, gating(found)), (0, []), none)
        said = [x for x in none.splitlines() if "read as a run of" in x]
        self.assertEqual(len(said), 1, none)
        self.assertTrue(said[0].endswith("INDELS format 2, fake-indel windows none: no indel "
                                         "record)"), said[0])
        # a sites-only VCF (the mega's) with no indel record there reads none too (W-R3-N2)
        write_sites_only_vcf(case.vcf, "chr21", [(len(PRE) + 6, "C", GUIDE[5])])
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc, found = case.run("completeness")
        sites = [x for x in err.getvalue().splitlines() if "read as a run of" in x]
        self.assertEqual(len(sites), 1, err.getvalue())
        self.assertTrue(sites[0].endswith("fake-indel windows none: no indel record)"), sites[0])
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(repo, "docs", "RELEASE_CHECKLIST.md")) as fh:
            checks = [x for x in fh.read().splitlines() if 'grep -o "fake-indel windows' in x]
        self.assertEqual(len(checks), 2, checks)
        gate = os.path.join(self.tmp, "gate")
        os.makedirs(gate)
        told = self.printed_reading((30, 30))
        for text, passes in ((told + none, True), (none, False),
                             (told + none.replace("none: no indel record", "not told"), False)):
            for name in ("indelrows_per-sample_soundness.err", "indelrows_gw_readme.err"):
                with open(os.path.join(gate, name), "w") as fh:
                    fh.write(text)
            for line in checks:
                out = subprocess.run(["bash", "-c", line], capture_output=True, text=True,
                                     env={"GATE": gate, "mode": "per-sample",
                                          "check": "soundness", "g": "readme:guide.txt",
                                          "PATH": os.environ.get("PATH", "/usr/bin:/bin")})
                self.assertEqual(out.stdout.strip(), "" if passes else "GATE_FAILED",
                                 (line, text[-300:]))

    def test_two_alts_of_one_position_on_a_slot_are_putative(self):
        """Every allele coded, a slot can hold two ALTs of one position (split records
        phased one by one): no haplotype, and CRISPRme resolves the carrier PUTATIVE
        (METHODS section 4, decision B4), the rule it resolves a call of unknown phase by,
        which so is the truth there (R-B2). S01's D slot holds the deletion and both ALTs
        of a spacer position, one of them the guide's base: listed for the guide's base,
        the row passes; a row that leaves S01 out fails (carrier-missing: it was counted
        as an overlap, never required), and one listing S01 for the other ALT's weaker
        target fails (carrier-not-strongest)."""
        seq = PRE + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        third = next(x for x in "ACGT" if x not in ("C", GUIDE[5]))
        recs = [(spos, "C", GUIDE[5], {"S01": "1|0", "S02": "1|0"}),
                (spos, "C", third, {"S01": "1|0"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"})]
        start = len(PRE)
        var = "%s,%s" % (tok("chr21", recs[2]), tok("chr21", recs[0]))
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         GUIDE + "TGG", "TGG", 0, 0, "X", var, "S01,S02")
        case = as_2100(Case(self.tmp, seq, recs, [row]))
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        write_results(case.res, [set_field(row, samples="S02")])
        rc, found = case.run()
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)}, {("carrier-missing", "S01")},
                         found)
        weak = list(GUIDE + "TGG")
        weak[5] = third.lower()
        var3 = "%s,%s" % (tok("chr21", recs[2]), tok("chr21", recs[1]))
        row3 = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                          "".join(weak), "TGG", 1, 0, "X", var3, "S01")
        write_results(case.res, [row, row3])
        rc, found = case.run()
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-not-strongest", "S01")}, found)

    def test_two_alts_of_one_slot_leave_the_reference_out(self):
        """The PUTATIVE rule counts a carrier's copies against the record's call (METHODS
        section 4): S01's D slot holds the deletion and two ALTs of a spacer position whose
        reference base is the guide's, so its two ALT copies fill its two slots and the
        reference is no option there, as CRISPRme 2.10.0 resolves it (candidates C and G,
        not A). Listed for the stronger ALT's target, the row passes; listed for the weaker
        one's (completeness then misses the stronger one's row), or for the perfect target,
        it fails."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        alts = [x for x in "ACGT" if x != GUIDE[5]][:2]
        recs = [(spos, GUIDE[5], alts[0], {"S01": "1|0"}),
                (spos, GUIDE[5], alts[1], {"S01": "1|0"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"})]
        start = len(PRE)
        ref23 = seq[start:start + 23]
        dtok = tok("chr21", recs[2])

        def run(samples_perfect, alt_k=None):
            tmp = os.path.join(self.tmp, "%s_%s" % (samples_perfect, alt_k))
            os.makedirs(tmp)
            rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref23, GUIDE + "TGG", "TGG",
                               0, 0, "X", dtok, samples_perfect)]
            if alt_k is not None:
                weak = list(GUIDE + "TGG")
                weak[5] = alts[alt_k].lower()
                rows.append(row_fields("chr21", start, "+", GUIDE + "NNN", ref23,
                                       "".join(weak), "TGG", 1, 0, "X", "%s,%s" % (
                                           dtok, tok("chr21", recs[alt_k])), "S01"))
            case = as_2100(Case(tmp, seq, recs, rows))
            # (one mismatch, no bulge: the rows written are every target in that budget)
            return [case.run(mode=m, budget=("1", "0", "0", "1"))
                    for m in ("soundness", "completeness")]

        verdicts = []
        for k in (0, 1):
            (rc_s, found_s), (rc_c, found_c) = run("S02", k)
            verdicts.append((rc_s, sorted({(f[1], f[3]) for f in gating(found_s)}),
                             rc_c, sorted({(f[1], f[3]) for f in gating(found_c)})))
        self.assertEqual(sorted(verdicts),
                         [(0, [], 0, []),
                          (1, [("carrier-not-strongest", "S01")], 1,
                           [("missing-putative-row", "S01")])], verdicts)
        (rc, found), _ = run("S01,S02")
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-without-allele", "S01")}, found)

    def test_two_records_of_one_snv_are_one_allele(self):
        """Every allele coded, two lines of one SNV REF>ALT are one allele, their calls
        merged as CRISPRitz 2.10.0's enricher merges them into the allele's one dictionary
        entry (its CHANGELOG, docs/MULTIALLELIC_ENRICHMENT.md section 2.3; R-B1): phased
        alike, slot by slot, a slot holding the ALT if either line puts it there; otherwise
        the call with more copies of the ALT, the earlier one on a tie, as it is. On the
        genome of the last record (no semantics given) they stay two lines."""
        seq = "ACGT" * 20
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", [(5, "A", "C", {"S01": "1|0", "S02": "1|0", "S03": "0/1",
                                                "S04": "1|0", "S05": "1", "S06": "0/1"}),
                                 (5, "A", "C", {"S01": "0|1", "S02": "1|0", "S03": "1|0",
                                                "S04": "0/1", "S05": "0|1", "S06": "1|1"})])
        sem = C.Semantics(True, {(None, None): {"all_alleles": True, "store_filter": True,
                                                "indels": (2, 2)}}, {})
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)], sem=sem)
        self.assertEqual(len(win.by_key[(5, "A", "C")]), 1)
        a = win.by_key[(5, "A", "C")][0]
        self.assertEqual({s: win.call(a, s) for s in ("S01", "S02", "S03", "S04", "S05",
                                                      "S06")},
                         {"S01": ((1, 1), True), "S02": ((1, 0), True), "S03": ((0, 1), False),
                          "S04": ((1, 0), True), "S05": ((1,), True), "S06": ((1, 1), True)})
        self.assertEqual(C.merge_call(((1, None), True), ((None, 0), True)), ((1, 0), True))
        self.assertEqual(C.merge_call(((None, None), True), ((None, 0), True)),
                         ((None, 0), True))
        old = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
        self.assertEqual(len(old.by_key[(5, "A", "C")]), 2)

    def test_a_phased_call_beside_an_unphased_one_of_its_allele_stays_phased(self):
        """R-B1, through the check: S01 is 1|0 on the first line of a SNV and 0/1 on a
        second line of the same REF>ALT, both beside the deletion D on its first slot.
        CRISPRitz 2.10.0's dictionary entry keeps 1|0 (the calls are not phased alike, and
        their copies tie), so CRISPRme lists S01 for D's target through the SNV, its slot's
        own sequence. Read with the calls made unphased (CRISPRme's own reader merge,
        site_alleles.merge_gt, which never sees such an entry), the check demanded the
        PUTATIVE rule's strongest target, the reference base there, and failed the row."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        alt_b = other(GUIDE[5])
        recs = [(spos, GUIDE[5], alt_b, {"S01": "1|0"}),
                (spos, GUIDE[5], alt_b, {"S01": "0/1"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"})]
        start = len(PRE)
        ref23 = seq[start:start + 23]
        dtok = tok("chr21", recs[2])
        snp = list(GUIDE + "TGG")
        snp[5] = alt_b.lower()
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref23, GUIDE + "TGG", "TGG",
                           0, 0, "X", dtok, "S02"),
                row_fields("chr21", start, "+", GUIDE + "NNN", ref23, "".join(snp), "TGG", 1,
                           0, "X", "%s,%s" % (dtok, tok("chr21", recs[0])), "S01")]
        case = as_2100(Case(self.tmp, seq, recs, rows))
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))

    def test_the_last_record_of_one_snv_drops_the_earlier_lines_carriers(self):
        """R-B3: on a genome of the last record (CRISPRitz 2.9.1 and earlier), whose
        dictionary entry of a position was the last line's, two lines of one SNV REF>ALT
        are read so too: the earlier line's carriers are lost (S05, here), its ALT coded
        by the later line (S03). The SNV restores the guide's base, beside the deletion D.
        CRISPRme lists S03 for the perfect target and S05 for D's target on the reference:
        the documented limit (limit-last-alt-wins), no failure; a row that leaves S03 out
        fails. S04, on both lines and the same slot, holds the ALT once: listed with S03,
        it agrees; left out, it fails (it was a slot of two records of one base)."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        ref_b = other(GUIDE[5])
        recs = [(spos, ref_b, GUIDE[5], {"S05": "1|0", "S04": "1|0"}),
                (spos, ref_b, GUIDE[5], {"S03": "1|0", "S04": "1|0"}),
                (dpos, GUIDE[9] + "C", GUIDE[9],
                 {"S02": "1|0", "S03": "1|0", "S04": "1|0", "S05": "1|0"})]
        start = len(PRE)
        ref23 = seq[start:start + 23]
        dtok = tok("chr21", recs[2])
        on_ref = list(GUIDE + "TGG")
        on_ref[5] = ref_b.lower()

        def rows(perfect):
            return [row_fields("chr21", start, "+", GUIDE + "NNN", ref23, GUIDE + "TGG", "TGG",
                               0, 0, "X", "%s,%s" % (dtok, tok("chr21", recs[1])), perfect),
                    row_fields("chr21", start, "+", GUIDE + "NNN", ref23, "".join(on_ref),
                               "TGG", 1, 0, "X", dtok, "S02,S05")]

        case = Case(self.tmp, seq, recs, rows("S03,S04"))
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
            self.assertIn("limit-last-alt-wins", classes(found), (mode, found))
        for perfect, gone in (("S04", "S03"), ("S03", "S04")):
            write_results(case.res, rows(perfect))
            rc, found = case.run(budget=("1", "0", "0", "1"))
            self.assertEqual(rc, 1, found)
            self.assertEqual({(f[1], f[3]) for f in gating(found)}, {("carrier-missing", gone)},
                             found)

    def test_two_lines_of_one_snv_count_their_copies_once(self):
        """R-B3, the PUTATIVE rule: two lines of one SNV REF>ALT that put the ALT on one slot
        of S06 (1|0 on both) give one copy of it, not two. S06's call of the deletion D is
        unphased, so it is listed for the strongest target its alleles allow: the guide's
        base at the SNV's position (its other slot holds the reference there), D's perfect
        target, as CRISPRme resolves it from the last line: the row agrees with S06, no line
        of either mode names it. Counted twice, the copies left the reference out, and the
        row listing S06 failed (or, beside the shadowed line, it was a PUTATIVE limit)."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        alt_b = other(GUIDE[5])
        recs = [(spos, GUIDE[5], alt_b, {"S06": "1|0"}),
                (spos, GUIDE[5], alt_b, {"S06": "1|0"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S06": "0/1", "S02": "1|0"})]
        start = len(PRE)
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[2]), "S02,S06")
        case = Case(self.tmp, seq, recs, [row])
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
            self.assertFalse([f for f in found if "S06" in f[3].split(",")], (mode, found))

    def test_the_alts_crisprme_cannot_tie_are_untied(self):
        """R2B-F1: CRISPRme's readers tie an alt of a line of several ALTs to the number every
        carrier of it holds, taking off it the numbers the position's other alts settle
        (METHODS section 4), over the alts of one character the dictionary lists there,
        a '*' among them, an indel's ALT not, and every sample of the VCF, whatever the
        panel. POS 5: C, G and T all tied (C's carriers share 1, which leaves G 2 and T 3);
        POS 9: C beside a '*' whose carriers hold both numbers, untied; POS 13: a carrier of
        the '*' alone settles 2, C tied; POS 17: C beside an insertion, untied; POS 21: C and
        G of one 1|2 carrier, untied; POS 25: S09, outside the panel, settles G's 2, C tied.
        Read line by line, the alts no carrier holds alone were untied (POS 5), and an
        SNV beside a '*' or an indel never (POS 9, 17)."""
        seq = "ACGT" * 20
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        samples = SAMPLES + ["S09"]
        with gzip.open(vcf, "wt") as fh:
            fh.write(vcf_text("chr21", [
                (5, "A", "C,G,T", {"S01": "1|2", "S02": "1|3"}),
                (9, "A", "C,*", {"S03": "1|2"}),
                (13, "A", "C,*", {"S03": "1|2", "S04": "0|2"}),
                (17, "A", "C,AT", {"S05": "1|2"}),
                (21, "A", "C,G", {"S06": "1|2"}),
                (25, "A", "C,G", {"S07": "1|2", "S09": "0|2"})], samples))
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)], panel=set(SAMPLES))
        self.assertEqual(sorted(a.key for a in win.untied),
                         [(9, "A", "C"), (17, "A", "C"), (21, "A", "C"), (21, "A", "G")])
        sem = C.Semantics(True, {(None, None): {"all_alleles": True, "store_filter": True,
                                                "indels": (2, 2)}}, {})
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)], panel=set(SAMPLES), sem=sem)
        self.assertEqual(sorted(a.key for a in win.untied),
                         [(9, "A", "C"), (17, "A", "C"), (21, "A", "C"), (21, "A", "G")])

    def test_an_snv_beside_a_star_or_an_indel_on_its_line_is_putative(self):
        """R2B-F1 through the check: the SNV T of a line G>T,* or G>T,GA (a '*' or an
        insertion as its second ALT), S01 1|2, beside the deletion D on S01's first slot.
        CRISPRme's readers list the SNV alone among the alts that tie numbers (no entry of
        an indel; the '*''s carriers are the SNV's), so it cannot tell which slot holds T
        and resolves S01 PUTATIVE: the strongest target its alleles allow, the guide's base
        at T's position, D's perfect target beside S02 (CRISPRitz 2.10.0 and the real CLI
        list so). That row is the documented limit (limit-putative-carrier), no failure in
        either mode; read with T's number off the VCF, it failed (variant-missing: S01's
        slot of D holds T), and so did S01's target through T, never reported."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST + LONG
        dpos, spos = len(PRE) + 10, len(PRE) + 4
        ref_b = GUIDE[3]
        alt_b = other(ref_b)
        start = len(PRE)
        for second in ("*", ref_b + "A"):
            tmp = os.path.join(self.tmp, "star" if second == "*" else "ins")
            os.makedirs(tmp)
            recs = [(spos, ref_b, alt_b + "," + second, {"S01": "1|2"}),
                    (dpos, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "1|0"})]
            row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                             GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[1]),
                             "S01,S02")
            case = as_2100(Case(tmp, seq, recs, [row]))
            for mode in ("soundness", "completeness"):
                rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
                self.assertEqual((rc, gating(found)), (0, []), (second, mode, found))
                if mode == "soundness":
                    self.assertEqual({(f[1], f[3]) for f in found},
                                     {("limit-putative-carrier", "S01")}, (second, found))

    def test_a_putative_carriers_target_is_listed_where_the_rule_picks_it(self):
        """R2B-F1, CRISPRme's own rows (the reviewer's world W, through the real CLI of
        2.10.0, 4 + 1 + 1, per sample; chr21:900-1249 of its reference, its coordinates less
        900): S01 holds the deletion D (110 TC>T) and, by the VCF, T of 104 G>T,* on D's slot,
        a 1-mismatch target; CRISPRme cannot tie T's number and resolves S01 PUTATIVE, so it
        lists S01 on D's perfect row (the guide's G at 104 under the alignment without
        bulges) and with T only under the bulged alignment that pairs T with the guide's T.
        The target through T, reported for S01 under a worse alignment than its best, is
        the documented limit (limit-putative-carrier), as the perfect row is; the check
        failed it (missing-best-alignment)."""
        seq = (
            "CACATTAGTCCCTTGTATGCAGGCGGTATCGGACGGCGCCCACACCTTGGAGGTATCCAGCGCAAGGCGC"
            "CATATCCGTACCTTACTATCGCGCGAACTTACTGAAATCTCGTAAGCAGGCTGGCATCTATACGTCAGTC"
            "CTAAACATAGCGAGCATTTCGCAGATGGGTCTCCGACGGTACCCCAAGGGTCGTTACCGACGCCGGGACG"
            "CCGCATATAAAGGTACGCCCGACCATTATACAGGTAGCCATCTGCGTCTGACATCGCATTTGAAACCCAG"
            "TAGGTACTGCCTTAGTTGCACTCCTAACTCATGTTAACGGACTTACGGGCACTAGCTTCTTACTGCCCTC"
        )
        recs = [(102, "CTG", "C", {"S01": "0|1"}), (104, "G", "T,*", {"S01": "1|2"}),
                (110, "TC", "T", {"S01": "1|0", "S02": "1|0"})]

        def blk(start, crrna, ref, alt, mm, bul, bt, var, who):
            return row_fields("chr21", start, "+", crrna, ref, alt,
                              alt.replace("-", "")[-3:], mm, bul, bt, var, who)
        perfect = blk(100, "ACTGAAATCTGTAAGCAGGCNNN", "ACTGAAATCTcgtAagcaGgCTG",
                      "ACTGAAATCTGTAAGCAGGCTGG", 0, 0, "X", "chr21_110_TC_T", "S01,S02")
        alt = [
            [blk(100, "A-CTGAAATCTGTAAGCAGGCNNN", "ACtg-AAATCTcgtAagcaGgCTG",
                 "ACtT-AAATCTGTAAGCAGGCTGG", 1, 2, "RNA,DNA", "chr21_104_G_T,chr21_110_TC_T",
                 "S01"),
             blk(100, "ACTG-AAATCTGTAAGCAGGCNNN", "ACT-GAAATCTcgtAagcaGgCTG",
                 "ACT-GAAATCTGTAAGCAGGCTGG", 0, 2, "RNA,DNA", "chr21_110_TC_T", "S02"),
             blk(100, "ACTGAAATCT-GTAAGCAGGCNNN", "ACTGAAATCTCGTAAGCAGGCTGG",
                 "ACTtAAATCTCGTAAGCAGGCTGG", 1, 1, "DNA", "chr21_104_G_T", "S01")],
            [blk(100, "A-CTGAAATCTGTAAGCAGGCNNN", "ACtg-AAATCTcgtAagcaGgCTG",
                 "ACtg-AAATCTGTAAGCAGGCTGG", 2, 2, "RNA,DNA", "chr21_110_TC_T", "S02"),
             blk(100, "ACTG-AAATCTGTAAGCAGGCNNN", "ACT-GAAATCTcgtAagcaGgCTG",
                 "ACT-TAAATCTGTAAGCAGGCTGG", 0, 2, "RNA,DNA", "chr21_104_G_T,chr21_110_TC_T",
                 "S01"),
             blk(101, "ACTGAAATCTGTAAGCAGGCNNN", "ctg-AAATCTcgtAagcaGgCTG",
                 "ctT-AAATCTGTAAGCAGGCTGG", 2, 1, "RNA", "chr21_104_G_T,chr21_110_TC_T",
                 "S01")],
            [blk(100, "ACTGAAATCTGTAAGCAGGCNNN", "ACTGAAATC-tcgtaagcagGCT",
                 "ACaaAtcTC-GTAAGCAGGCTGG", 4, 1, "RNA", "chr21_102_CTG_C", "S01"),
             blk(100, "ACTGAAATCT-GTAAGCAGGCNNN", "ACTGAAATCTCGTAAGCAGGCTGG",
                 "ACTtAAATCTCGTAAGCAGGCTGG", 1, 1, "DNA", "chr21_104_G_T", "S01"),
             blk(99, "ACTG-AAATCTGTAAGCAGGCNNN", "tactGAAATCTcgtAagcaGgCTG",
                 "tactGAAATCTGTAAGCAGGCTGG", 4, 1, "DNA", "chr21_110_TC_T", "S02")],
            [blk(101, "ACTGAAATC-TGTAAGCAGGCNNN", "ctg-AAATCTcGTAAGCAGGCTGG",
                 "ctT-AAATCTcGTAAGCAGGCTGG", 3, 2, "RNA,DNA", "chr21_104_G_T", "S01"),
             blk(99, "ACTGAAATCT-GTAAGCAGGCNNN", "-tactgAaaTCtcgtaagcagGCT",
                 "-tacAAATCTCGTAAGCAGGCTGG", 3, 2, "RNA,DNA", "chr21_102_CTG_C", "S01"),
             blk(98, "ACTGAAATCT-GTAAGCAGGCNNN", "ttactgAaaTCtcgtaagcagGCT",
                 "ttacAAATCTCGTAAGCAGGCTGG", 4, 1, "DNA", "chr21_102_CTG_C", "S01")]]
        case = as_2100(Case(self.tmp, seq, recs, [perfect]))
        write_results(case.res, [blocks_row(b) for b in alt],
                      name="G+NRG_hg38+v_4+2_all_results_with_alternative_alignments.tsv")
        run_markers(case.res, "per-sample")
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode)
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
            self.assertIn(("limit-putative-carrier", "S01"), {(f[1], f[3]) for f in found},
                          (mode, found))

    def test_a_line_of_several_alts_merges_its_listing_of_an_snv(self):
        """R2B-F3: every allele coded, CRISPRitz 2.10.0's enricher merges a listing of an
        SNV REF>ALT into its first listing at the position wherever it is listed, on a line
        of several ALTs too, and CRISPRme reads that one entry: its calls written as a line
        of one ALT writes them, a slot holding the ALT where either puts it there when they
        are phased alike, else the call with more copies, the earlier on a tie. POS 5: A>C,T
        after A>C: S07's 1/0 leaves its 0|1 as it is, S06's 1|1 joins; T keeps its line and
        its number 2. POS 9: A>C,T before A>C: C gets a line of its own, S01's 1|0 and 0|1
        give 1|1, and T stays on the first line, 2 on S01's second slot. On the genome of the
        last record they stay two listings."""
        seq = "ACGT" * 20
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", [(5, "A", "C", {"S07": "0|1"}),
                                 (5, "A", "C,T", {"S07": "1/0", "S08": "2|0", "S06": "1|1"}),
                                 (9, "A", "C,T", {"S01": "1|2", "S02": "2|0"}),
                                 (9, "A", "C", {"S01": "0|1", "S03": "1|0"})])
        sem = C.Semantics(True, {(None, None): {"all_alleles": True, "store_filter": True,
                                                "indels": (2, 2)}}, {})
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)], sem=sem)
        self.assertEqual({k: len(v) for k, v in win.by_key.items()},
                         {(5, "A", "C"): 1, (5, "A", "T"): 1, (9, "A", "C"): 1, (9, "A", "T"): 1})
        a5, t5 = win.by_key[(5, "A", "C")][0], win.by_key[(5, "A", "T")][0]
        self.assertEqual({s: win.call(a5, s) for s in ("S06", "S07")},
                         {"S06": ((1, 1), True), "S07": ((0, 1), True)})
        self.assertEqual((t5.number, win.call(t5, "S08")), (2, ((2, 0), True)))
        a9, t9 = win.by_key[(9, "A", "C")][0], win.by_key[(9, "A", "T")][0]
        self.assertEqual((a9.number, win.line_alleles[a9.line], win.call(a9, "S01"),
                          win.call(a9, "S03")), (1, [a9], ((1, 1), True), ((1, 0), True)))
        self.assertEqual((t9.number, win.call(t9, "S01"), win.line_alleles[t9.line]),
                         (2, ((1, 2), True), [None, t9]))
        self.assertEqual(sorted(a.key for a in win.untied), [])
        old = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
        self.assertEqual((len(old.by_key[(5, "A", "C")]), len(old.by_key[(9, "A", "C")])),
                         (2, 2))

    def test_a_listing_merged_from_a_line_of_several_alts_keeps_its_phase(self):
        """R2B-F3 through the check (the reviewer's world W2): C>A, S07 0|1, then C>A,G,
        S07 1/0 and S08 2|0, at a spacer position whose reference base C is a mismatch and A
        the guide's base, beside the deletion D on the first slot of S07 and S06. CRISPRitz
        2.10.0's entry of C>A keeps S07's 0|1 (the calls are not phased alike, and their
        copies tie), so CRISPRme lists S07 with S06 on D's one-mismatch target, C at that
        position, which the check now passes in both modes. Read on its own, the second
        line's 1/0 sent S07 to the PUTATIVE rule, whose strongest target, A's perfect one,
        was demanded and S07's listing failed."""
        seq = PRE + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST + LONG
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        third = next(x for x in "ACGT" if x not in ("C", GUIDE[5]))
        recs = [(spos, "C", GUIDE[5], {"S07": "0|1"}),
                (spos, "C", GUIDE[5] + "," + third, {"S07": "1/0", "S08": "2|0"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S07": "1|0", "S06": "1|0"})]
        start = len(PRE)
        on_ref = list(GUIDE + "TGG")
        on_ref[5] = "c"
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         "".join(on_ref), "TGG", 1, 0, "X", tok("chr21", recs[2]), "S06,S07")
        case = as_2100(Case(self.tmp, seq, recs, [row]))
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
            self.assertFalse([f for f in found if "S07" in f[3].split(",")], (mode, found))

    def test_a_run_that_cannot_tell_how_it_read_a_repeated_allele_is_refused(self):
        """A position whose lines are C>A and C>A,G (the reviewer's world W, its four-allele
        position left out) is read otherwise by a genome of every allele (one allele, the
        listings merged) and by one of the last record (the earlier line's carriers lost),
        and nothing in it tells them apart: both genomes code C, A and G there. A run of
        2.10.0 whose variant index flags no Q, whose VCF has no other split site, is refused
        (exit 2), as where a biallelic split site cannot be read; the index's Q flag tells,
        and the run is checked. The check sought only split sites of biallelic SNV lines,
        and read such a run as coding the last record."""
        seq = PRE + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST + LONG
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        recs = [(spos, "C", GUIDE[5], {"S07": "0|1"}),
                (spos, "C", GUIDE[5] + ",G", {"S07": "1/0", "S08": "2|0"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S07": "1|0", "S06": "1|0"})]
        enriched = seq[:spos - 1] + "V" + seq[spos:]   # C, A and G on either reading
        start = len(PRE)
        on_ref = list(GUIDE + "TGG")
        on_ref[5] = "c"
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         "".join(on_ref), "TGG", 1, 0, "X", tok("chr21", recs[2]), "S06,S07")
        for flags, want in ((0, 2), (4, 0)):
            tmp = os.path.join(self.tmp, "q%d" % flags)
            os.makedirs(tmp)
            case = as_2100(Case(tmp, seq, recs, [row]), var_flags=flags, enriched=enriched)
            rc, found = case.run(budget=("1", "0", "0", "1"))
            self.assertEqual(rc, want, (flags, found))

    def test_a_star_line_tells_nothing_of_the_coding(self):
        """R2B-F4: a line whose ALT is a '*' (a spanning deletion's) beside an SNV line at
        one position is no split SNP site: no IUPAC code holds a '*' and CRISPRitz codes
        nothing for it, so the genome codes the SNV there on either reading. Twelve such
        positions before the VCF's one split site took the check's whole sample of sites,
        and a run of 2.10.0 whose variant index flags no Q (and no four-base position)
        was read as coding the last record of a split site; now the split site tells, and
        the run is read as coding every allele."""
        seq = "ACGT" * 40
        code = {frozenset("AC"): "M", frozenset("AG"): "R", frozenset("AT"): "W",
                frozenset("CG"): "S", frozenset("CT"): "Y", frozenset("GT"): "K",
                frozenset("ACT"): "H"}
        recs, enriched = [], list(seq)
        for i in range(12):
            q = 4 * i + 5          # an A
            recs += [(q, "A", "*", {"S01": "0|1"}), (q, "A", "C", {"S01": "1|0"})]
            enriched[q - 1] = code[frozenset("AC")]
        recs += [(101, "A", "C", {"S02": "1|0"}), (101, "A", "T", {"S03": "1|0"})]
        enriched[100] = code[frozenset("ACT")]
        case = as_2100(Case(self.tmp, seq, recs, []), var_flags=0, enriched="".join(enriched))
        vcfs = [C.Vcf(case.vcf)]
        self.assertEqual([x[0] for x in C.split_sites(vcfs[0], "chr21")], [100])
        sem, why = C.run_semantics(case.res, vcfs, ["chr21"], ["v"], "hg38")
        self.assertIsNone(why)
        self.assertTrue(sem.all_alleles(vcfs[0], "chr21"))

    def test_a_target_at_the_window_edge_follows_the_indels_index_format(self):
        """R10-N1: a 27-nt SaCas9 target that starts at its fake-indel window's first base.
        An INDELS index of format 1 drops it (the separator N in its pad):
        limit-window-edge, not demanded; one of format 2 holds it: demanded, a missing row
        fails. A run that says neither (no index found, no format-1 NOTE) is refused."""
        seq, recs, _start = world_window_edge(self.tmp)
        budget = ("2", "0", "0", "2")
        for fmt in (1, 2):
            tmp = os.path.join(self.tmp, "f%d" % fmt)
            os.makedirs(tmp)
            case = as_2100(Case(tmp, seq, recs, [], pam=SA_PAM), indels_format=fmt, budget=1,
                           guide_len=21, pam="NNGRRT")
            rc, found = case.run(mode="completeness", budget=budget, guide=G21)
            if fmt == 1:
                self.assertEqual((rc, gating(found)), (0, []), found)
                self.assertIn("limit-window-edge", classes(found), found)
            else:
                self.assertEqual(rc, 1, found)
                self.assertIn("missing-row", classes(gating(found)), found)
        tmp = os.path.join(self.tmp, "gone")
        os.makedirs(tmp)
        case = as_2100(Case(tmp, seq, recs, [], pam=SA_PAM), budget=1, guide_len=21,
                       pam="NNGRRT")
        shutil.rmtree(os.path.join(tmp, "work", "genome_library"))
        rc, _found = case.run(mode="completeness", budget=budget, guide=G21)
        self.assertEqual(rc, 2)

    def test_a_bulged_target_at_the_window_edge_is_read_by_its_site(self):
        """R2B-F2: whether an index holds a target at its fake-indel window's edge is its
        site's: the guide region (the guide's length from the PAM) and, on format 1, its
        pad. An RNA-bulged target whose own bases leave format 1's pad clear while its
        site's guide region does not is limit-window-edge on format 1, demanded on format 2;
        a DNA-bulged one whose own first base lies in the pad while its site's pad is clear
        is demanded on both (a missing row fails); one whose site's guide region starts on
        the separator N is held by no index (limit-window-reach). Read from the target's
        end base, the first was demanded of format 1, the second left to format 1 as an
        edge, the third demanded of format 2."""
        want = {("rna", 1): (0, "limit-window-edge"), ("rna", 2): (1, "missing-row"),
                ("dna", 1): (1, "missing-row"), ("dna", 2): (1, "missing-row"),
                ("mnv", 1): (0, "limit-window-reach"), ("mnv", 2): (0, "limit-window-reach")}
        for (kind, fmt), (rc_want, cls) in sorted(want.items()):
            seq, recs, budget = world_site_edge(kind)
            tmp = os.path.join(self.tmp, "%s%d" % (kind, fmt))
            os.makedirs(tmp)
            case = as_2100(Case(tmp, seq, recs, [], pam=SA_PAM), indels_format=fmt, budget=2,
                           guide_len=21, pam="NNGRRT")
            rc, found = case.run(mode="completeness", budget=budget, guide=G21)
            self.assertEqual(rc, rc_want, (kind, fmt, found))
            self.assertIn(cls, classes(gating(found) if rc_want else found), (kind, fmt, found))
            if not rc_want:
                self.assertEqual(gating(found), [], (kind, fmt, found))

    def test_a_target_past_the_window_is_not_demanded(self):
        """A target that runs past its record's fake-indel window (here a DNA-bulged
        SaCas9 target, 28 bases from a base before the window to the inserted T) is held
        by no index, whatever its format: limit-window-reach, never demanded."""
        seq, recs, _start = world_window_edge(self.tmp, bulged=True)
        case = as_2100(Case(self.tmp, seq, recs, [], pam=SA_PAM), budget=2, guide_len=21,
                       pam="NNGRRT")
        rc, found = case.run(mode="completeness", budget=("2", "1", "0", "3"), guide=G21)
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertIn("limit-window-reach", classes(found), found)

    def test_the_fake_indel_window_is_the_runs_own(self):
        """The fake-indel window is read from the run's indel log (each row's window start
        and end): on windows of 25/27 (CRISPRitz 2.9.1) a DNA-bulged SaCas9 target that
        starts a base before the window is held by no index (limit-window-reach, not
        demanded); on windows of 30/30 (2.10.0) it is in the window, demanded, and a missing
        row fails. Read as 25/27 whatever the log, it was never demanded."""
        seq, recs, _start = world_window_edge(self.tmp, bulged=True)
        budget = ("2", "1", "0", "3")
        for window, want in (((25, 27), (0, "limit-window-reach")),
                             ((30, 30), (1, "missing-row"))):
            tmp = os.path.join(self.tmp, "w%d" % window[0])
            os.makedirs(tmp)
            case = as_2100(Case(tmp, seq, recs, [], pam=SA_PAM), budget=2, guide_len=21,
                           pam="NNGRRT", window=window)
            rc, found = case.run(mode="completeness", budget=budget, guide=G21)
            self.assertEqual(rc, want[0], (window, found))
            self.assertIn(want[1], classes(gating(found) if want[0] else found),
                          (window, found))
            if not want[0]:
                self.assertEqual(gating(found), [], found)
                line = [f[2] for f in found if f[1] == "limit-window-reach"][0]
                self.assertIn("(25 bases before the record's first REF base, 27 after its "
                              "last), or whose site's guide region does (a target of a 27-nt "
                              "guide and PAM with a DNA bulge, or one with an RNA bulge at the "
                              "window's first base, say)", line)

    def test_the_reach_example_is_the_windows(self):
        """The example a limit-window-reach line gives is the shortest target past the run's
        own windows: L + b past min(LEFT + 2, RIGHT + 1), a 27-nt SaCas9 target with a DNA
        bulge on 25/27, a 29-nt one with 3 on 30/30 (with one it fits), any 33-nt one."""
        self.assertEqual(C.reach_example((25, 27), 27),
                         "a target of a 27-nt guide and PAM with a DNA bulge")
        self.assertEqual(C.reach_example((30, 30), 29),
                         "a target of a 29-nt guide and PAM with 3 DNA bulges")
        self.assertEqual(C.reach_example((30, 30), 23),
                         "a target of a 23-nt guide and PAM with 9 DNA bulges")
        self.assertEqual(C.reach_example((30, 30), 33), "a target of a 33-nt guide and PAM")

    def test_a_run_whose_records_say_no_window_is_refused(self):
        """A run of 2.10.0 whose working directory holds no indel log of the dataset, and
        whose INDELS index records no window (index_meta's indel_window), does not say the
        fake-indel windows its targets were searched in: refused (exit 2) where its VCF has
        indel records. With the INDELS sub-index's record of 30/30 it is read so (the
        DNA-bulged target is demanded)."""
        seq, recs, _start = world_window_edge(self.tmp, bulged=True)
        budget = ("2", "1", "0", "3")
        case = as_2100(Case(self.tmp, seq, recs, [], pam=SA_PAM), budget=2, guide_len=21,
                       pam="NNGRRT", window=None)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc, _found = case.run(mode="completeness", budget=budget, guide=G21)
        self.assertEqual(rc, 2)
        self.assertIn("the fake-indel windows of v on chr21", err.getvalue())
        sub = os.path.join(self.tmp, "work", "genome_library", "21-NNGRRT_2_hg38+v_INDELS",
                           "NNGRRT_2_fakechr21")
        with open(os.path.join(sub, ".index_meta.json"), "w") as fh:
            json.dump({"index_format": 2, "indel_window": [30, 30]}, fh)
        rc, found = case.run(mode="completeness", budget=budget, guide=G21)
        self.assertEqual(rc, 1, found)
        self.assertIn("missing-row", classes(gating(found)), found)

    def test_the_refusal_says_why_the_windows_are_not_told(self):
        """W-R1-N4.5: the refusal of a run whose records do not tell its fake-indel windows
        says why: no indel log of the dataset in the working directory; one with no row;
        one whose rows hold two widths, which no one enrichment writes; one with a row of
        no indel log's layout."""
        seq, recs, _start = world_window_edge(self.tmp, bulged=True)
        budget = ("2", "1", "0", "3")
        case = as_2100(Case(self.tmp, seq, recs, [], pam=SA_PAM), budget=2, guide_len=21,
                       pam="NNGRRT", window=None)
        log = os.path.join(self.tmp, "work", "Dictionaries", "log_indels_v", "logchr21.txt")
        head = "CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
        row = "chr21_%d-%d_%d\tS01\t.\t0.1\tchr21_%d_A_AT\t0,1\tN\n"
        pos = recs[0][0]
        for body, said in (
                (None, "the run's working directory holds no indel log of it, "
                       "Dictionaries/log_indels_v/logchr21.txt[.gz]"),
                (head, "its indel log %s holds no row" % log),
                (head + row % (pos - 26, pos + 27, 1, pos)
                 + row % (pos - 31, pos + 30, 2, pos),
                 "%s holds fake-indel windows of 25/27 and 30/30 bases" % log),
                (head + "chr21\tS01\t.\t0.1\tchr21_%d_A_AT\t0,1\tN\n" % pos,
                 "%s holds a row of no indel log's layout" % log)):
            if body is not None:
                os.makedirs(os.path.dirname(log), exist_ok=True)
                with open(log, "w") as fh:
                    fh.write(body)
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                rc, _found = case.run(mode="completeness", budget=budget, guide=G21)
            self.assertEqual(rc, 2, err.getvalue())
            self.assertIn("the fake-indel windows of v on chr21 (%s" % said, err.getvalue())

    def test_indel_log_window(self):
        """indel_log_window reads a log's windows from its rows (CHR <chrom>_<start>-<end>_<id>,
        indel <chrom>_<POS>_<REF>_<ALT>), a chromosome name holding '_' or '-' too; a log
        of two widths is no one enrichment's (ValueError); one with no row says none."""
        path = os.path.join(self.tmp, "logchrUn_K-1.txt")
        rows = ["chrUn_K-1_%d-%d_%d\tS01\t.\t0.1\tchrUn_K-1_%d_%s_%s\t0,1\tN"
                % (p - 31, p - 1 + len(r) + 30, k, p, r, a)
                for k, (p, r, a) in enumerate(((100, "A", "AT"), (200, "ACG", "A"),
                                               (300, "AC", "TG")))]
        write_text = lambda body: open(path, "w").write(  # noqa: E731
            "CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n" + body)
        write_text("".join(r + "\n" for r in rows))
        self.assertEqual(C.indel_log_window(path), (30, 30))
        write_text(rows[0] + "\n" + "chr1_74-127_9\tS01\t.\t0.1\tchr1_100_A_AT\t0,1\tN\n")
        with self.assertRaises(ValueError):
            C.indel_log_window(path)
        write_text("")
        self.assertIsNone(C.indel_log_window(path))

    def test_window_cut(self):
        """window_cut on its own: 25 reference bases before the record's first REF base and
        27 after its last; the site's guide region, the guide's length from the PAM
        (``full_n`` 27: SaCas9's 21-nt guide and NNGRRT), inside them on either format; the
        edge rule only on an INDELS index of format 1, on the site's guide region's far end
        (+ strand, 3' PAM: the left end), within the index's budget of the window's edge
        (R2B-F2: the site, not the target's end base, which an RNA bulge moves in and a DNA
        bulge out)."""
        class Win(object):
            chrom = "chr21"
            vcfs = ["v"]
            line_vcf = [0]
        win = Win()
        al_ = C.Allele(100, "A", "AT", 0, 1)   # first REF base 99 (0-based)
        ctx = type("Ctx", (), {"bdna": 0, "brna": 0, "pam_5prime": False})()
        win.sem = None
        place = [(q, 0) for q in range(74, 100)] + [(99, 1)]   # 74 = 99 - 25
        self.assertIsNone(C.window_cut(win, ctx, al_, place, "+", 27))
        self.assertEqual(C.window_cut(win, ctx, al_, [(73, 0)] + place, "+", 27), "reach")
        self.assertIsNone(C.window_cut(win, ctx, al_, [(99, 1)] + [(q, 0) for q in
                                                                range(100, 127)], "+", 27))
        self.assertEqual(C.window_cut(win, ctx, al_, [(99, 1)] + [(q, 0) for q in
                                                               range(100, 128)], "+", 27),
                         "reach")   # 127: past the 27 bases after the REF's last, 99
        # one RNA bulge, the target at the window's first base: its site's guide region
        # starts a base before it, on the separator N -- no index holds it
        rna0 = [(q, 0) for q in range(74, 99)] + [(99, 1)]
        self.assertEqual(C.window_cut(win, ctx, al_, rna0, "+", 27), "reach")
        # one RNA bulge, the target 2 bases in; one DNA bulge, the target a base in
        rna2 = [(q, 0) for q in range(76, 100)] + [(99, 1), (100, 0)]
        dna1 = [(q, 0) for q in range(75, 100)] + [(99, 1), (100, 0), (101, 0)]
        self.assertEqual((len(rna0), len(rna2), len(dna1)), (26, 26, 28))
        for fmt, pad, want in ((1, 1, ("edge", None, None)), (1, 2, ("edge", "edge", None)),
                               (2, 2, (None, None, None))):
            win.sem = C.Semantics(True, {(None, None): {"indels": (fmt, pad)}}, {})
            self.assertEqual(tuple(C.window_cut(win, ctx, al_, t, "+", 27)
                                   for t in (place, rna2, dna1)), want, (fmt, pad))
            # the - strand's guide is on the right: the left end is its PAM's
            self.assertIsNone(C.window_cut(win, ctx, al_, place, "-", 27))
            self.assertEqual(C.window_cut(win, ctx, al_, rna0, "+", 27), "reach")
        # the - strand, 3' PAM: the far end on the right, 27 bases after the REF's last (99);
        # a target of 26 bases whose last is 2 in: an RNA-bulged one of a 27-nt site, or a
        # 26-nt site's
        win.sem = C.Semantics(True, {(None, None): {"indels": (1, 2)}}, {})
        right = [(99, 1)] + [(q, 0) for q in range(100, 125)]
        self.assertEqual(C.window_cut(win, ctx, al_, right, "-", 27), "edge")
        self.assertIsNone(C.window_cut(win, ctx, al_, right, "-", 26))
        # the run's windows of 30/30 (CRISPRitz 2.10.0, its indel log): 30 bases before the
        # REF's first base (99) and 30 after its last
        win.sem = C.Semantics(True, {(None, None): {"indels": (2, 2), "window": (30, 30)}}, {})
        wide = [(q, 0) for q in range(69, 100)] + [(99, 1)]          # 69 = 99 - 30
        self.assertIsNone(C.window_cut(win, ctx, al_, wide, "+", 27))
        self.assertEqual(C.window_cut(win, ctx, al_, [(68, 0)] + wide, "+", 27), "reach")
        self.assertIsNone(C.window_cut(win, ctx, al_, [(99, 1)] + [(q, 0) for q in
                                                                range(100, 130)], "+", 27))
        self.assertEqual(C.window_cut(win, ctx, al_, [(99, 1)] + [(q, 0) for q in
                                                               range(100, 131)], "+", 27),
                         "reach")
        self.assertIsNone(C.window_cut(win, ctx, al_, [(q, 0) for q in range(73, 99)]
                                       + [(99, 1)], "+", 27))         # rna0 of 25/27

    def test_a_run_of_2_10_0_writes_no_capped_greedy_row(self):
        """A per-sample run of 2.10.0 reports a window over CRISPRME_IUPAC_CAP by the
        haplotypes its carriers have, on the classic dictionaries too: no greedy row a
        cluster can chain through (capped_sets), which a run of 2.8.1 may write."""
        seq = "ACGTTGCAAG" * 8
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", [(3 + 2 * i, seq[2 + 2 * i], other(seq[2 + 2 * i]),
                                  {"S01": "1|0"}) for i in range(C.IUPAC_CAP + 2)])
        ctx = type("Ctx", (), {"snp_union": False, "sem": None})()
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
        self.assertTrue(C.capped_sets(win, ctx, 0, len(seq)))
        ctx.sem = C.Semantics(True, {}, {})
        self.assertEqual(C.capped_sets(win, ctx, 0, len(seq)), [])

    def test_max_degeneracy_opens_a_snp_row_over_its_cap(self):
        """--max-degeneracy K (2.10.0; the run's .Params.txt Max_degeneracy): a SNP row
        resting on more than K ambiguous columns is cut back or left out (METHODS section
        5), so a population-level alignment over more than K positions of several options
        is a row CRISPRme may write or not (snp_row_cps' fifth set), never one the chain
        counts on; without the cap it is a row."""
        site = list(GUIDE + "TGG")
        for i in (3, 8):
            site[i] = other(site[i])
        seq = PRE + "".join(site) + POST
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", [(len(PRE) + i + 1, site[i], GUIDE[i], {"S01": "0|1"})
                                 for i in (3, 8)])

        class A(object):
            mm, bdna, brna, max_total_edits = 0, 0, 0, 0
        ctx = C.Ctx(A(), "NRG", 3, False, [GUIDE])
        ctx.snp_union = True
        cp = C.cluster_of(len(PRE), "+", 23, 23, False)
        for k, want in ((None, (True, False)), (1, (False, True)), (2, (True, False))):
            ctx.max_degeneracy = k
            win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
            rows, _c, _d, _f, deg, _g = C.snp_row_cps(win, ctx, GUIDE, ctx.full(GUIDE), 0, len(seq))
            self.assertEqual((cp in rows, cp in deg), want, k)

    def test_the_semantics_of_a_run(self):
        """run_semantics: no .indexes_searched.txt, a run before 2.10.0 (the last record,
        stores of every FILTER, INDELS format 1); with it, the INDELS index's format and
        budget from its .bin, the store's rule from the log (an older store's note), none
        with no log."""
        seq, recs, row, _D, _J = world_filter_second("PASS")
        case = Case(self.tmp, seq, recs, [row])
        vcfs = [C.Vcf(case.vcf)]
        sem, why = C.run_semantics(case.res, vcfs, ["chr21"], None, "hg38")
        self.assertIsNone(why)
        self.assertFalse(sem.release2100)
        self.assertEqual((sem.all_alleles(vcfs[0], "chr21"), sem.store_filter(vcfs[0], "chr21"),
                          sem.indels(vcfs[0], "chr21")), (False, False, (1, None)))
        self.assertTrue(sem.capped)
        as_2100(case, log="", indels_format=2, budget=3)
        sem, why = C.run_semantics(case.res, vcfs, ["chr21"], ["v"], "hg38")
        self.assertTrue(sem.release2100)
        self.assertFalse(sem.capped)
        self.assertEqual((sem.all_alleles(vcfs[0], "chr21"), sem.store_filter(vcfs[0], "chr21"),
                          sem.indels(vcfs[0], "chr21")), (True, True, (2, 3)))
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("NOTE [indel-snp]: the phased indel genotype store of v on chr21 was "
                     "written before CRISPRme 2.10.0: ...\n")
        sem, _why = C.run_semantics(case.res, vcfs, ["chr21"], ["v"], "hg38")
        self.assertFalse(sem.store_filter(vcfs[0], "chr21"))
        self.assertTrue(sem.values[("v", "chr21")]["store"])
        os.remove(os.path.join(case.res, "log_verbose.txt"))
        sem, _why = C.run_semantics(case.res, vcfs, ["chr21"], ["v"], "hg38")
        self.assertIsNone(sem.values[("v", "chr21")]["store_filter"])
        self.assertIsNone(sem.values[("v", "chr21")]["store"])
        lib = os.path.join(case.tmp, "work", "genome_library")
        shutil.rmtree(os.path.join(lib, "20-NRG_3_hg38+v_INDELS"))
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("NOTE: 20-NRG_3_hg38+v_INDELS is in index format 1 (built by CRISPRitz "
                     "2.9.1 or earlier): it lacks ...\n")
        sem, _why = C.run_semantics(case.res, vcfs, ["chr21"], ["v"], "hg38")
        self.assertEqual(sem.indels(vcfs[0], "chr21"), (1, 3))

    def test_every_chromosome_reads_the_store_from_the_runs_own_log(self):
        """W-R1-F1: run_semantics reads each chromosome's store FILTER rule from the run's
        log_verbose.txt, whatever chromosomes come before it. The indel log it reads a
        chromosome's fake-indel windows from (Dictionaries/log_indels_<dataset>) never takes
        that log's place. With a log_verbose.txt that notes no older store, every chromosome
        reads 2.10.0's store, chr21 after chr20 (which has no indel log) too. With one that
        notes chr21's store as older, chr21 alone reads it so. Without log_verbose.txt none
        is told, chr22 after chr21 (which has an indel log) too."""
        seq, recs, row, _D, _J = world_filter_second("PASS")
        case = Case(self.tmp, seq, recs, [row])
        vcfs = [C.Vcf(case.vcf)]
        as_2100(case, log="", indels_format=2, budget=3)
        chroms = ["chr20", "chr21", "chr22"]

        def stores():
            sem, _why = C.run_semantics(case.res, vcfs, chroms, ["v"], "hg38")
            return [sem.values[("v", ch)]["store_filter"] for ch in chroms]

        self.assertEqual(stores(), [True, True, True])
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("NOTE [indel-snp]: the phased indel genotype store of v on chr21 was "
                     "written before CRISPRme 2.10.0: ...\n")
        self.assertEqual(stores(), [True, False, True])
        os.remove(os.path.join(case.res, "log_verbose.txt"))
        self.assertEqual(stores(), [None, None, None])

    def add_chromosome(self, case, chrom, records):
        """A second chromosome of the run's dataset v: its reference beside chr21's, and a
        VCF of its own in the same folder as v.vcf.gz (read as v, the run's one dataset).
        Returns the VCF's path."""
        seq = PRE + GUIDE + "TGG" + POST + PRE
        write_fasta(os.path.join(case.tmp, "Genomes", chrom + ".fa"), seq, chrom)
        path = os.path.join(os.path.dirname(case.vcf), "v_%s.vcf.gz" % chrom)
        write_vcf(path, chrom, [(p, seq[p - 1:p - 1 + n], alt, gts, filt)
                                for p, n, alt, gts, filt in records])
        return path

    def check_all(self, case, vcfs):
        """main in completeness mode without --chrom (every chromosome the VCFs hold), as a
        genome-wide check runs: (exit status, stderr)."""
        out = os.path.join(case.tmp, "all.tsv")
        argv = ["--results", case.res, "--reference", os.path.join(case.tmp, "Genomes"),
                "--pam", case.pam, "--mm", "4", "--bDNA", "1", "--bRNA", "1",
                "--max-total-edits", "6", "--mode", "completeness", "--guide", GUIDE,
                "--out", out]
        for v in vcfs:
            argv += ["--vcf", v]
        made = fill_markers(case.res)
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
                rc = C.main(argv)
        finally:
            for p in made:
                os.remove(p)
        return rc, err.getvalue()

    def test_a_check_of_every_chromosome_reads_the_runs_own_log(self):
        """W-R1-F1 through main, in completeness mode without --chrom. The run's log notes
        chr21's store as one written before 2.10.0 (which keeps J, so the row listing S02
        alone is right). With a second VCF of the dataset for chr20, which sorts first and
        has no indel log, chr21 is still read as that older store and the run passes."""
        seq, recs, row, _D, _J = world_filter_second("LowQual")
        case = as_2100(Case(self.tmp, seq, recs, [row]))
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("NOTE [indel-snp]: the phased indel genotype store of v on chr21 was "
                     "written before CRISPRme 2.10.0: it holds ...\n")
        v20 = self.add_chromosome(case, "chr20", [(len(PRE) + 3, 1, "A" if GUIDE[2] != "A"
                                                   else "C", {"S01": "1|0"}, "PASS")])
        rc, err = self.check_all(case, [case.vcf, v20])
        self.assertEqual(rc, 0, err)
        self.assertRegex(err, r"v on chr21: [^;]*store of every FILTER")
        self.assertRegex(err, r"v on chr20: [^;]*store of PASS and '\.'")

    def test_a_run_without_its_log_is_refused_on_every_chromosome(self):
        """W-R1-F1 through main, the other way: a run of 2.10.0 with no log_verbose.txt does
        not say which indels its store held, so a chr22 VCF holding a LowQual indel is
        refused (exit 2), though chr21's indel log sorts before it; chr21, whose indels all
        pass, is not."""
        seq, recs, row, _D, _J = world_filter_second("PASS")
        case = as_2100(Case(self.tmp, seq, recs, [row]), log=None)
        v22 = self.add_chromosome(case, "chr22", [(len(PRE) + 5, 2, GUIDE[4],
                                                   {"S03": "1|0"}, "LowQual")])
        rc, err = self.check_all(case, [case.vcf, v22])
        self.assertEqual(rc, 2, err)
        self.assertIn("which indels the phased indel genotype store of v held on chr22", err)
        self.assertNotIn("held on chr21", err)


# --------------------------------------------------------------------------- #
# stage 2a's review round 3 (s2a-review-r3-A and -B), fixed in stage 2b's lane C
# --------------------------------------------------------------------------- #
class StageTwoARoundThree(Base):
    """A record written in lower case, a target whose alignments a PUTATIVE co-carrier
    splits between two carrier sets, an indel written on two VCF lines, the narrowing of
    an alt's number by a line of one ALT, the checker's independence of CRISPRme's code,
    and a split site whose last line is a '*' on a genome of the last record."""

    def test_a_record_written_in_lower_case_is_found_on_its_rows(self):
        """R3B-F1: a row names a record as its VCF writes it, and HPRC writes 1,546 indel
        records of chr22 in lower case. The checker reads every record in upper case, and
        looked a row's variants up by its raw tokens: the right row of a lower-case
        deletion was never found, and every carrier of every target of it failed
        (missing-row). Read in upper case, the row passes, and one that leaves a carrier
        out still fails."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        pos1 = len(PRE) + 10
        rec = (pos1, (GUIDE[9] + "C").lower(), GUIDE[9].lower(), {"S01": "1|0", "S02": "0|1"})
        start = len(PRE)
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", rec), "S01,S02")
        self.assertIn("_tc_t", row[2 + OFF["var"]])
        case = Case(self.tmp, seq, [rec], [row])
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
        write_results(case.res, [set_field(row, samples="S01")])
        rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"))
        self.assertEqual({(f[1], f[3]) for f in gating(found)}, {("missing-row", "S02")},
                         found)

    def test_an_indel_written_on_two_lines_is_one_allele(self):
        """R3B-F3 (the reviewer's world WD in miniature): the deletion D on two VCF lines of
        one REF>ALT, D1 then D2, beside the SNV q that gives the guide's base. Two records
        of one REF>ALT are one allele (decision B2; CRISPRme 2.10.0's readers of the phased
        indel store merge them as CRISPRitz merges an SNV's): S01 (D1 1|0, D2 1|1, q 0|1)
        holds D on both slots, the perfect target on slot 1; S05 (D1 0|1, D2 1|0, q 1|0)
        on slot 0; S02 (D1 1|0 alone) the one-mismatch target. The check read D's call off
        the first line that lists the sample, so S01's and S05's slot of q was no slot of
        D, and the perfect row failed for them (carrier-without-variants)."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        ref_b = other(GUIDE[5])
        dref, dalt = GUIDE[9] + "C", GUIDE[9]
        recs = [(spos, ref_b, GUIDE[5], {"S01": "0|1", "S05": "1|0"}),
                (dpos, dref, dalt, {"S01": "1|0", "S02": "1|0", "S05": "0|1"}),
                (dpos, dref, dalt, {"S01": "1|1", "S05": "1|0"})]
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(self._vcf(seq, recs))])
        d = win.by_key[(dpos, dref, dalt)]
        self.assertEqual(len(d), 1)
        self.assertEqual({s: win.call(d[0], s) for s in ("S01", "S02", "S05")},
                         {"S01": ((1, 1), True), "S02": ((1, 0), True), "S05": ((1, 1), True)})
        start = len(PRE)
        ref23 = seq[start:start + 23]
        dtok = tok("chr21", recs[1])
        on_ref = list(GUIDE + "TGG")
        on_ref[5] = ref_b.lower()
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref23, GUIDE + "TGG", "TGG", 0,
                           0, "X", "%s,%s" % (dtok, tok("chr21", recs[0])), "S01,S05"),
                row_fields("chr21", start, "+", GUIDE + "NNN", ref23, "".join(on_ref), "TGG",
                           1, 0, "X", dtok, "S01,S02,S05")]
        case = Case(self.tmp, seq, recs, rows)
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
        # the carriers only the first line lists are carriers: S02 left off fails
        write_results(case.res, [rows[0], set_field(rows[1], samples="S01,S05")])
        rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"))
        self.assertEqual({(f[1], f[3]) for f in gating(found)}, {("missing-row", "S02")},
                         found)

    def _vcf(self, seq, recs, name="w.vcf.gz", samples=SAMPLES):
        path = os.path.join(self.tmp, name)
        with gzip.open(path, "wt") as fh:
            fh.write(vcf_text("chr21", recs, samples))
        return path

    def test_a_line_of_one_alt_settles_no_number_of_a_line_of_several(self):
        """R3A-1, the documented view: CRISPRitz 2.10.0 lists a position's lines of one ALT
        (and the listing it merges from two lines, written as one of a line of one ALT)
        beside a line of several ALTs, whose raw calls number its own ALTs. METHODS section
        4: an alt whose carriers hold only 0 and 1 settles no number. POS 5: A>C (S01 1|0)
        beside A>G,AT (S02 1|2: G on slot 0, the insertion on slot 1): G cannot be tied
        (untied: PUTATIVE). POS 9: A>C (S03 1|0) beside A>T,C (S04 1|2): C's listing is
        merged into A>C's, T untied. POS 13: one line A>C,G,T (S05 1|2, S06 1|3), all tied.
        The narrowing took the 1 of the line of one ALT off G and off T, which tied G to 2 (the
        insertion's number) and T to 2 (C's), reading each on the other slot."""
        seq = "ACGT" * 20
        vcf = self._vcf(seq, [(5, "A", "C", {"S01": "1|0"}), (5, "A", "G,AT", {"S02": "1|2"}),
                              (9, "A", "C", {"S03": "1|0"}), (9, "A", "T,C", {"S04": "1|2"}),
                              (13, "A", "C,G,T", {"S05": "1|2", "S06": "1|3"})])
        sem = C.Semantics(True, {(None, None): {"all_alleles": True, "store_filter": True,
                                                "indels": (2, 2)}}, {})
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)], sem=sem)
        self.assertEqual(sorted(a.key for a in win.untied), [(5, "A", "G"), (9, "A", "T")])
        # a genome of the last record holds the last line's entries alone: G untied there
        # too, and both alts of A>T,C, S04's 1|2 telling neither
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(vcf)])
        self.assertEqual(sorted(a.key for a in win.untied),
                         [(5, "A", "G"), (9, "A", "C"), (9, "A", "T")])

    def test_the_checker_copies_no_crisprme_code(self):
        """R3B-F4: no run of four code lines or more (comments and docstrings left out,
        whitespace collapsed; a run of imports aside) is the same in the checker and in a
        Python file of CRISPRme's: every file of the repository's tree, in every folder
        (PostProcess and its supportFunctions and fixtures, scripts, seq_script, pages, the
        oracle and benchmark scripts of test/, plot_generation_paper, paper, crisprme.py),
        the unit tests (test_*.py) aside. Window._untie was site_alleles.alt_tokens with
        its narrowing loop copied line for line, which the test of function names let
        through, as R-B1's merge_call (site_alleles.merge_gt) before it. C-R1-2: the test
        read only the top level of PostProcess, scripts, seq_script and the root, where
        the docs said every file."""
        root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            os.pardir))
        checker = os.path.join(root, "scripts", "check_indel_rows.py")
        mine = code_lines(checker)
        others = []
        for d, subdirs, names in os.walk(root):
            subdirs[:] = sorted(x for x in subdirs
                                if not x.startswith(".") and x != "__pycache__")
            for name in sorted(names):
                path = os.path.join(d, name)
                if (name.endswith(".py") and not name.startswith("test_")
                        and os.path.isfile(path) and os.path.abspath(path) !=
                        os.path.abspath(checker)):
                    others.append(path)
        rel = {os.path.relpath(p, root) for p in others}
        for want in ("PostProcess/site_alleles.py", "crisprme.py",
                     "PostProcess/supportFunctions/loadSample/associateSample.py",
                     "test/indel_snp_oracle/oracle.py"):
            self.assertIn(want, rel)
        found = []
        for path in others:
            for n, i, j, text in shared_runs(mine, code_lines(path), 4):
                found.append("%d lines: checker line %d = %s line %d: %s"
                             % (n, i, os.path.relpath(path, root), j, " / ".join(text)))
        self.assertEqual(found, [])

    def test_a_star_as_the_last_line_of_a_split_site_on_a_genome_of_the_last_record(self):
        """Stage 2a's open item: on a genome of CRISPRitz 2.9.1 the dictionary entry of a
        position is its last line's, a line of '*' alone too (a spanning deletion's), while
        the genome keeps the code the earlier lines wrote (a '*' adds no base): C>G (S03
        1|0, the guide's base, on the slot of the deletion D) then C>* (S04 0|1). The SNV
        is in no store, so CRISPRme 2.8.x lists S03 for D's target on the reference base:
        the documented last-record limit (limit-last-alt-wins), no failure. The check read
        the line before the '*' as the last, and failed the row. On a genome of every
        allele the '*' shadows nothing."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 6
        ref_b = other(GUIDE[5])
        recs = [(spos, ref_b, GUIDE[5], {"S03": "1|0"}),
                (spos, ref_b, "*", {"S04": "0|1"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S02": "1|0", "S03": "1|0"})]
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(self._vcf(seq, recs))])
        self.assertEqual({a.key for a in win.shadowed}, {(spos, ref_b, GUIDE[5])})
        sem = C.Semantics(True, {(None, None): {"all_alleles": True, "store_filter": True,
                                                "indels": (2, 2)}}, {})
        win = C.Window(seq, "chr21", 0, len(seq), [C.Vcf(self._vcf(seq, recs))], sem=sem)
        self.assertEqual(win.shadowed, set())
        start = len(PRE)
        on_ref = list(GUIDE + "TGG")
        on_ref[5] = ref_b.lower()
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         "".join(on_ref), "TGG", 1, 0, "X", tok("chr21", recs[2]), "S02,S03")
        os.makedirs(os.path.join(self.tmp, "old"))
        case = Case(os.path.join(self.tmp, "old"), seq, recs, [row])
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode, budget=("1", "0", "0", "1"))
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
            self.assertIn("limit-last-alt-wins", classes(found), (mode, found))
        # a run of 2.10.0 is read otherwise there, so where nothing in its records tells
        # how its genome was coded (no Q flag, no split site of two bases), it is refused;
        # its index's Q flag tells, and S03's perfect target is demanded
        for flags, want in ((0, 2), (4, 1)):
            tmp = os.path.join(self.tmp, "q%d" % flags)
            os.makedirs(tmp)
            case = as_2100(Case(tmp, seq, recs, [row]), var_flags=flags)
            rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"))
            self.assertEqual(rc, want, (flags, found))
            if want == 1:
                self.assertEqual({(f[1], f[3]) for f in gating(found)},
                                 {("missing-row", "S03")}, found)

    def test_a_putative_co_carrier_splits_a_targets_alignments_between_two_carrier_sets(self):
        """R3B-F2, CRISPRme's own rows (world WF2, the reviewer's WQ site A with three of its
        carriers, through the real CLI of 2.10.0, 4 + 1 + 1, per sample; chr21:900-1249 of its
        reference, coordinates less 900, its S12 written S08 here): S03 holds the deletion R
        (110 TC>T) and T at 104 on one slot, a target CRISPRme reports by two alignments, the
        1-mismatch one (CFD 0.364) and a bulged one (1 mismatch, a DNA and an RNA bulge, CFD
        0.457). S08 (R 1|0, A 1/0, T 0/1) is PUTATIVE: T where the bulged alignment pairs the
        column with the guide's T, A on the 1-mismatch one, so the merge, which keeps one
        entry per carrier set (METHODS section 4), lists the 1-mismatch alignment for {S03}
        and the bulged one for {S03, S08}. S03's best by edits stands in one group, its best
        CFD in the other: judged one group at a time, both failed (missing-best-alignment,
        and S07's A target, two ALTs on its slot, likewise); judged together, as the merge
        keeps them, the rows pass."""
        seq = ("CGCGAGGAAATCGCGCTAGGGCTTCTCGTTATAGGTTATAGGTAAGGTATATTCCGCAATTACGCTCCAAGATG"
               "CAGCAAACTCACCTAACGTATGTTGCACTCAAATCTCGTAAGCAGGCTGGTAGCAGAGTTTCAATTCCGCCCC"
               "GGGAGTTCTAGCCGAACTTGACCGTCCTTGCCAGGGGGCGAACTATGGACGCGCACGTAGACGACCCGCAACAC"
               "ATCCCGACTCATTAGGCCTCACATCGGCCATAAGCGGCAGGGTAGTAGGCCACGAAAGCCGGAACGGTTTCACA"
               "AGTACCACCTACGAACACGCAGATGACTCTACGACGCATTTTATGCAGTGGGCGA")
        recs = [(104, "C", "A", {"S07": "1|0", "S08": "1/0"}), (104, "C", "G", {"S07": "1|0"}),
                (104, "C", "T", {"S03": "1|0", "S08": "0/1"}),
                (110, "TC", "T", {"S03": "1|0", "S07": "1|0", "S08": "1|0"})]
        rows = [
            ('integrated', [
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTGAAATCTGTAAGCAGGCTGG', 'TGG', 0, 0, 'X', 'chr21_104_C_G,chr21_110_TC_T', 'S07'),
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTGAAATCTGTAAGCAGGCTGG', 'TGG', 0, 0, 'X', 'chr21_104_C_G,chr21_110_TC_T', 'S07'),
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTGAAATCTGTAAGCAGGCTGG', 'TGG', 0, 0, 'X', 'chr21_104_C_G,chr21_110_TC_T', 'S07'),
            ]),
            ('alternative', [
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTaAAATCTGTAAGCAGGCTGG', 'TGG', 1, 0, 'X', 'chr21_104_C_A,chr21_110_TC_T', 'S08'),
                (100, '+', 'ACTGAAATCT-GTAAGCAGGCNNN', 'ACTcAAATCTCGTAAGCAGGCTGG', 'ACTGAAATCTCGTAAGCAGGCTGG', 'TGG', 0, 1, 'DNA', 'chr21_104_C_G', 'S07'),
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTaAAATCTGTAAGCAGGCTGG', 'TGG', 1, 0, 'X', 'chr21_104_C_A,chr21_110_TC_T', 'S08'),
            ]),
            ('alternative', [
                (100, '+', 'A-CTGAAATCTGTAAGCAGGCNNN', 'ACtc-AAATCTcgtAagcaGgCTG', 'ACtT-AAATCTGTAAGCAGGCTGG', 'TGG', 1, 2, 'RNA,DNA', 'chr21_104_C_T,chr21_110_TC_T', 'S03,S08'),
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTaAAATCTGTAAGCAGGCTGG', 'TGG', 1, 0, 'X', 'chr21_104_C_A,chr21_110_TC_T', 'S08'),
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTtAAATCTGTAAGCAGGCTGG', 'TGG', 1, 0, 'X', 'chr21_104_C_T,chr21_110_TC_T', 'S03'),
            ]),
            ('alternative', [
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTtAAATCTGTAAGCAGGCTGG', 'TGG', 1, 0, 'X', 'chr21_104_C_T,chr21_110_TC_T', 'S03'),
                (100, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ACTcAAATCTcgtAagcaGgCTG', 'ACTtAAATCTGTAAGCAGGCTGG', 'TGG', 1, 0, 'X', 'chr21_104_C_T,chr21_110_TC_T', 'S03'),
                (100, '+', 'ACTGAAATCT-GTAAGCAGGCNNN', 'ACTcAAATCTCGTAAGCAGGCTGG', 'ACTGAAATCTCGTAAGCAGGCTGG', 'TGG', 0, 1, 'DNA', 'chr21_104_C_G', 'S07'),
            ]),
            ('alternative', [
                (100, '+', 'A-CTGAAATCTGTAAGCAGGCNNN', 'ACtc-AAATCTcgtAagcaGgCTG', 'ACta-AAATCTGTAAGCAGGCTGG', 'TGG', 2, 2, 'RNA,DNA', 'chr21_104_C_A,chr21_110_TC_T', 'S07'),
                (100, '+', 'ACTGA-AATCTGTAAGCAGGCNNN', 'ACT-cAAATCTcgtAagcaGgCTG', 'ACT-AAAATCTGTAAGCAGGCTGG', 'TGG', 0, 2, 'RNA,DNA', 'chr21_104_C_A,chr21_110_TC_T', 'S07,S08'),
                (100, '+', 'ACTGAAATCT-GTAAGCAGGCNNN', 'ACTcAAATCTCGTAAGCAGGCTGG', 'ACTaAAATCTCGTAAGCAGGCTGG', 'TGG', 1, 1, 'DNA', 'chr21_104_C_A', 'S07,S08'),
            ]),
            ('alternative', [
                (100, '+', 'ACTGAAATC-TGTAAGCAGGCNNN', 'ACTcAAATCTcGTAAGCAGGCTGG', 'ACTGAAATCTcGTAAGCAGGCTGG', 'TGG', 1, 1, 'DNA', 'chr21_104_C_G', 'S07'),
                (100, '+', 'ACTG-AAATCTGTAAGCAGGCNNN', 'ACT-CAAATCTcgtAagcaGgCTG', 'ACT-TAAATCTGTAAGCAGGCTGG', 'TGG', 0, 2, 'RNA,DNA', 'chr21_104_C_T,chr21_110_TC_T', 'S03,S08'),
                (100, '+', 'ACTGAAATCT-GTAAGCAGGCNNN', 'ACTcAAATCTCGTAAGCAGGCTGG', 'ACTtAAATCTCGTAAGCAGGCTGG', 'TGG', 1, 1, 'DNA', 'chr21_104_C_T', 'S08'),
            ]),
            ('alternative', [
                (100, '+', 'ACTGAAATC-TGTAAGCAGGCNNN', 'ACTcAAATCTcGTAAGCAGGCTGG', 'ACTaAAATCTcGTAAGCAGGCTGG', 'TGG', 2, 1, 'DNA', 'chr21_104_C_A', 'S07,S08'),
                (100, '+', 'ACTGAAATCT-GTAAGCAGGCNNN', 'ACTcAAATCTCGTAAGCAGGCTGG', 'ACTaAAATCTCGTAAGCAGGCTGG', 'TGG', 1, 1, 'DNA', 'chr21_104_C_A', 'S07,S08'),
                (101, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ctc-AAATCTcgtAagcaGgCTG', 'ctT-AAATCTGTAAGCAGGCTGG', 'TGG', 2, 1, 'RNA', 'chr21_104_C_T,chr21_110_TC_T', 'S03,S08'),
            ]),
            ('alternative', [
                (100, '+', 'ACTGAAATC-TGTAAGCAGGCNNN', 'ACT-cAAatCTcgtAagcaGgCTG', 'ACT-AAAatCTGTAAGCAGGCTGG', 'TGG', 2, 2, 'RNA,DNA', 'chr21_104_C_A,chr21_110_TC_T', 'S07,S08'),
                (100, '+', 'ACTGAAATCT-GTAAGCAGGCNNN', 'ACTcAAATCTCGTAAGCAGGCTGG', 'ACTtAAATCTCGTAAGCAGGCTGG', 'TGG', 1, 1, 'DNA', 'chr21_104_C_T', 'S08'),
                (99, '+', 'ACTGA-AATCTGTAAGCAGGCNNN', 'cactcAAATCTcgtAagcaGgCTG', 'cactAAAATCTGTAAGCAGGCTGG', 'TGG', 4, 1, 'DNA', 'chr21_104_C_A,chr21_110_TC_T', 'S07,S08'),
            ]),
            ('alternative', [
                (101, '+', 'ACTGAAATC-TGTAAGCAGGCNNN', 'ctc-AAATCTcGTAAGCAGGCTGG', 'ctT-AAATCTcGTAAGCAGGCTGG', 'TGG', 3, 2, 'RNA,DNA', 'chr21_104_C_T', 'S08'),
                (100, '+', 'AC-TGAAATCTGTAAGCAGGCNNN', 'ACTc-AAATCTcgtAagcaGgCTG', 'ACTa-AAATCTGTAAGCAGGCTGG', 'TGG', 1, 2, 'RNA,DNA', 'chr21_104_C_A,chr21_110_TC_T', 'S07'),
                (101, '+', 'ACTGAAATCTGTAAGCAGGCNNN', 'ctc-AAATCTcgtAagcaGgCTG', 'cta-AAATCTGTAAGCAGGCTGG', 'TGG', 3, 1, 'RNA', 'chr21_104_C_A,chr21_110_TC_T', 'S07'),
            ]),
        ]
        lines = {"integrated": [], "alternative": []}
        for table, blocks in rows:
            lines[table].append(blocks_row([row_fields("chr21", *b) for b in blocks]))
        case = as_2100(Case(self.tmp, seq, recs, lines["integrated"]))
        write_results(case.res, lines["alternative"],
                      name="G+NRG_hg38+v_4+2_all_results_with_alternative_alignments.tsv")
        run_markers(case.res, "per-sample")
        for mode in ("soundness", "completeness"):
            rc, found = case.run(mode=mode)
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
        # the reviewer's mutation m9: S08 added to the 1-mismatch T rows, which its strongest
        # target on that alignment (A) is not, fails (and judged one group at a time, S03's
        # failure went away there, the wrong row passing)
        bad = {"integrated": [], "alternative": []}
        for table, blocks in rows:
            bad[table].append(blocks_row([row_fields("chr21", *(
                b[:10] + ("S03,S08",) if b[4] == "ACTtAAATCTGTAAGCAGGCTGG" and b[7] == 0
                else b)) for b in blocks]))
        write_results(case.res, bad["alternative"],
                      name="G+NRG_hg38+v_4+2_all_results_with_alternative_alignments.tsv")
        rc, found = case.run(mode="both")
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("carrier-not-strongest", "S08")}, found)

    def test_two_lines_of_several_alts_at_one_position_are_refused(self):
        """R3A-1's remaining case: two lines that each list several ALTs at one position (here
        a spacer base of the deletion D's target, C>A,T then C>G,CT). Where every allele is
        coded, CRISPRitz 2.10.0 lists both lines' one-base ALTs in one dictionary value, each
        with its line's own allele numbers, which CRISPRme's readers cannot tell apart
        (METHODS section 4), so what CRISPRme lists for their carriers is not modeled: the
        run is refused (exit 2), with the position, and writes no findings. On a genome of
        the last record (CRISPRitz 2.9.1) the dictionary holds the last line alone: it is
        checked. The check read the two lines as two numberings it knew, and judged
        CRISPRme against them."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 4
        ref_b = GUIDE[3]
        alts = [x for x in "ACGT" if x != ref_b]
        recs = [(spos, ref_b, alts[0] + "," + alts[1], {"S03": "1|2"}),
                (spos, ref_b, alts[2] + "," + ref_b + "T", {"S04": "1|2"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})]
        start = len(PRE)
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[2]), "S01,S02")
        os.makedirs(os.path.join(self.tmp, "new"))
        case = as_2100(Case(os.path.join(self.tmp, "new"), seq, recs, [row]))
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc, found = case.run(mode="both", budget=("1", "0", "0", "1"))
        self.assertEqual((rc, found), (2, []), err.getvalue())
        self.assertIn("chr21:%d holds two or more VCF lines that each list several ALTs"
                      % spos, err.getvalue())
        os.makedirs(os.path.join(self.tmp, "old"))
        case = Case(os.path.join(self.tmp, "old"), seq, recs, [row])
        rc, found = case.run(mode="both", budget=("1", "0", "0", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_refused_run_stops_at_the_first_unit_that_meets_it(self):
        """C-R1-6: the run of the test above is refused as soon as a work unit meets the
        position, the other units not checked (here one unit of 30 bases, the first in the
        order units are run); the refusal came only after every unit had been checked, the
        whole check time of a genome-wide run."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 4
        ref_b = GUIDE[3]
        alts = [x for x in "ACGT" if x != ref_b]
        recs = [(spos, ref_b, alts[0] + "," + alts[1], {"S03": "1|2"}),
                (spos, ref_b, alts[2] + "," + ref_b + "T", {"S04": "1|2"}),
                (dpos, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})]
        start = len(PRE)
        row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                         GUIDE + "TGG", "TGG", 0, 0, "X", tok("chr21", recs[2]), "S01,S02")
        case = as_2100(Case(self.tmp, seq, recs, [row]))
        ran = []
        real = C.run_unit

        def counted(unit):
            ran.append(unit[:4])
            return real(unit)
        err = io.StringIO()
        with mock.patch.object(C, "run_unit", counted), contextlib.redirect_stderr(err):
            rc, found = case.run(mode="completeness", budget=("1", "0", "0", "1"),
                                 extra=("--tile", "30"))
        self.assertEqual((rc, found), (2, []), err.getvalue())
        self.assertIn("chr21:%d holds two or more VCF lines" % spos, err.getvalue())
        self.assertEqual(ran, [("comp", "chr21", 0, 30)])

    def test_an_snv_whose_alt_is_in_lower_case_is_refused(self):
        """The sweep of C-R1-1: CRISPRitz codes no lower-case ALT (a lower-case ALT is not a
        SNP ALT, every release), while its dictionary lists the carriers; read in upper case
        as VCF 4.x reads it, the check put the SNV on its carriers' haplotypes and demanded
        targets the search cannot find. Not modeled: the run is refused (exit 2) with the
        position, on every release. A lower-case REF with an upper-case ALT is coded (the
        enricher reads the reference's base) and is checked, as is an indel written in lower
        case (CRISPRme names it in upper case)."""
        seq = PRE + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + POST
        dpos, spos = len(PRE) + 10, len(PRE) + 4
        ref_b = GUIDE[3]
        alt_b = next(x for x in "ACGT" if x != ref_b)
        start = len(PRE)
        for name, snv, refused in (("low", (spos, ref_b, alt_b.lower(), {"S03": "1|0"}), True),
                                   ("ref", (spos, ref_b.lower(), alt_b, {"S03": "1|0"}), False)):
            recs = [snv, (dpos, (GUIDE[9] + "C").lower(), GUIDE[9].lower(),
                          {"S01": "1|0", "S02": "0|1"})]
            row = row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                             GUIDE + "TGG", "TGG", 0, 0, "X",
                             "chr21_%d_%s_%s" % (dpos, GUIDE[9] + "C", GUIDE[9]), "S01,S02")
            for release in ("2100", "281"):
                tmp = os.path.join(self.tmp, name + release)
                os.makedirs(tmp)
                case = Case(tmp, seq, recs, [row])
                if release == "2100":
                    case = as_2100(case)
                err = io.StringIO()
                with contextlib.redirect_stderr(err):
                    rc, found = case.run(mode="both", budget=("1", "0", "0", "1"))
                if refused:
                    self.assertEqual((rc, found), (2, []), (name, release, err.getvalue()))
                    self.assertIn("chr21:%d holds an SNV whose ALT the VCF writes in lower "
                                  "case" % spos, err.getvalue())
                else:
                    self.assertEqual((rc, gating(found)), (0, []), (name, release, found))

    def test_a_snp_row_lists_an_untied_alts_carrier_as_an_unphased_one(self):
        """R3A-1, the SNP rows (the reviewer's world Y, site F): S06 is 1|2 on a line C>A,CT
        (A on slot 0, the insertion on slot 1) and 1|0 on r, an SNV of the target; A's
        number cannot be tied (every carrier of A holds 1 and 2), so CRISPRme reads S06's A
        as one copy on an unknown slot and lists it, as an unphased individual, on every
        SNP row of a subset of its alleles (METHODS section 4): the A-only row too, which
        neither of its slots holds. The SNP-row report counts that listing as the
        documented limit (limit-snp-row-unphased), not as a carrier without the target."""
        seq = PRE + GUIDE[:5] + other(GUIDE[5]) + GUIDE[6:12] + other(GUIDE[12]) + \
            GUIDE[13:] + "TGG" + POST
        qpos, rpos = len(PRE) + 6, len(PRE) + 13
        recs = [(qpos, other(GUIDE[5]), GUIDE[5] + "," + other(GUIDE[5]) + "T",
                 {"S06": "1|2"}),
                (rpos, other(GUIDE[12]), GUIDE[12], {"S06": "1|0", "S07": "0|1"})]
        start = len(PRE)
        a_only = list(GUIDE + "TGG")
        a_only[12] = other(GUIDE[12]).lower()
        rows = [row_fields("chr21", start, "+", GUIDE + "NNN", seq[start:start + 23],
                           "".join(a_only), "TGG", 1, 0, "X",
                           tok("chr21", (qpos, other(GUIDE[5]), GUIDE[5])), "S06")]
        case = as_2100(Case(self.tmp, seq, recs, rows))
        rc, found = case.run(extra=("--snp-rows", "report"), budget=("1", "0", "0", "1"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertEqual({(f[1], f[3]) for f in found},
                         {("limit-snp-row-unphased", "S06")}, found)

    def test_a_run_whose_indels_were_read_without_their_phased_store_is_refused(self):
        """The sweep of C-R2-1: a search that found no phased indel store for a genotyped
        dataset's chromosome (an index built with CRISPRME_INDEL_SNP=0, or a store whose
        compile failed) reads each indel's carriers from the indel log, phase unknown, and
        says so in its log_verbose.txt. The checker models the store's phased reading, so
        it judged such a run's right rows wrong (a world of a deletion on two lines, the
        store removed: 43 gating findings, every one the unphased reading): refused (exit
        2), in every mode; the note of another chromosome or dataset refuses nothing."""
        case, _ = world_indel_snp(self.tmp)
        note = ("NOTE [indel-snp]: no phased indel genotype store of %s on %s: each indel's "
                "carriers are its indel log's, their phase unknown, ...\n")
        as_2100(case, log="STARTING ANALYSIS FOR v\n" + note % ("v", "chr21"))
        for mode in ("soundness", "completeness", "both"):
            rc, _found = case.run(mode=mode)
            self.assertEqual(rc, 2, mode)
        for other_run in (note % ("v", "chr22"), note % ("w", "chr21"), ""):
            with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
                fh.write("STARTING ANALYSIS FOR v\n" + other_run)
            rc, found = case.run()
            self.assertEqual((rc, gating(found)), (0, []), (other_run, found))

    def test_a_run_of_2_10_0_without_its_log_is_refused_where_its_vcf_holds_indels(self):
        """C-R3-4: only the log_verbose.txt of a run of 2.10.0 says that the run read a
        genotyped dataset's indels with no phased store, or that its SNP+indel co-occurrence
        was off (the NOTEs above). A run whose log is gone was read as a default run that
        read its store: the review's world WDNS (wd.py, its store removed) without its log
        gave 43 gating findings, each of them the reading the NOTE refuses. Such a run
        cannot say how it read its indels: refused (exit 2) in every mode where the VCF
        holds indel or MNV records whose FILTER passes; the same run with its log is
        checked, a run of 2.10.0 without its log whose VCF holds no such record is checked,
        and so is a run before 2.10.0, whose log records neither (as before)."""
        case, _ = world_indel_snp(self.tmp)
        as_2100(case, log=None)
        for mode in ("soundness", "completeness", "both"):
            rc, _found = case.run(mode=mode)
            self.assertEqual(rc, 2, mode)
        with open(os.path.join(case.res, "log_verbose.txt"), "w") as fh:
            fh.write("STARTING ANALYSIS FOR v\n")
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)
        # an MNV is a record of the indel path too; an SNV alone is not (completeness:
        # a chromosome is read with no row on it)
        seq = PRE + GUIDE[:5] + "C" + GUIDE[6:] + "TGG" + POST
        snv = (len(PRE) + 6, "C", GUIDE[5], {"S05": "1|0"})
        mnv = (len(PRE) + 6, "CA", GUIDE[5] + "T", {"S06": "1|0"})
        assert seq[mnv[0] - 1:mnv[0] + 1] == mnv[1]
        for recs, refused in (([snv], False), ([snv, mnv], True)):
            tmp = tempfile.mkdtemp(dir=self.tmp)
            other = as_2100(Case(tmp, seq, recs, []), log=None)
            rc, found = other.run(mode="completeness")
            self.assertEqual(rc, 2 if refused else 0, (recs, found))
        # a run before 2.10.0 (no .indexes_searched.txt) without its log: checked
        case, _ = world_indel_snp(os.path.join(self.tmp, "old"))
        rc, found = case.run()
        self.assertEqual((rc, gating(found)), (0, []), found)


# (GGT)4 + G11: a guide of a GGT repeat and a run of G, and its PAM GGG
REPEAT_TARGET = "GGTGGTGGTGGTGGGGGGGGGGG"
REPEAT_GUIDE = REPEAT_TARGET[:20]


def world_repeat_places():
    """M-R3-N2 (chr22:26,153,539 (+), VEGFA site 3, 2 DNA + 2 RNA bulges) in miniature, with
    no mismatch and no bulge: an insertion of two units of a GGT repeat, A>AGGTGGT, before
    the reference's GGTGGTGGT and a run of G, and an SNV that turns the third reference
    unit's T into a G. The target (GGT)4 G11 is on every haplotype of the insertion, at one
    of two places: from the 4th added base on the haplotype of the insertion alone (S01,
    S02; and S05, whose SNV is on its other slot), where the target holds the third unit's
    T as the reference has it; from the 1st added base on the haplotype of the insertion
    and the SNV (S03, S04), where the target holds the SNV, its 15th base. S06 holds the SNV
    alone. Both rows start at the insertion's anchor, with the same bases: the insertion's
    row lists S01, S02 and S05, the row of the insertion and the SNV S03 and S04 (METHODS
    section 4: the haplotypes left with one sequence form one row, which names the SNP alts
    they hold). Returns (reference, records, the insertion's row, the other row)."""
    seq = LONG + PRE + "A" + "GGTGGTGGT" + "G" * 11 + POST
    anchor0 = len(LONG) + len(PRE)
    assert seq[anchor0] == "A" and seq[anchor0 + 9] == "T"
    ins = (anchor0 + 1, "A", "AGGTGGT",
           {"S01": "1|0", "S02": "0|1", "S03": "1|0", "S04": "0|1", "S05": "1|0"})
    snv = (anchor0 + 10, "T", "G", {"S03": "1|0", "S04": "0|1", "S05": "0|1", "S06": "1|1"})
    ref = seq[anchor0:anchor0 + 23]
    row_ins = row_fields("chr21", anchor0, "+", REPEAT_GUIDE + "NNN", ref, REPEAT_TARGET,
                         "GGG", 0, 0, "X", tok("chr21", ins), "S01,S02,S05",
                         spacer=REPEAT_GUIDE + "NNN")
    row_both = row_fields("chr21", anchor0, "+", REPEAT_GUIDE + "NNN", ref, REPEAT_TARGET,
                          "GGG", 0, 0, "X", "%s,%s" % (tok("chr21", ins), tok("chr21", snv)),
                          "S03,S04", spacer=REPEAT_GUIDE + "NNN")
    return seq, [ins, snv], row_ins, row_both


class RepeatPlaces(Base):
    """A haplotype holds a row's target where its SNP alts in the target are the ones the
    row names (M-R3-N2). In a repeat, a target that starts inside an insertion starts at
    the anchor whichever added base it starts on, and another haplotype can read the same
    bases there from another added base through an SNV the row does not name: that reading
    is another row's target, the row naming the SNV too."""

    def test_two_places_of_one_target_in_a_repeat_are_two_rows(self):
        """CRISPRme's two rows are right: S03 and S04 hold the target from the 1st added
        base through the SNV, and are listed on the row that names it; the row naming the
        insertion alone, whose target their haplotype holds from no added base without the
        SNV, does not list them. The check read their reading of the same bases at the
        anchor as the insertion's row's target and demanded them there (carrier-missing,
        gating, on each block), at any budget."""
        seq, recs, row_ins, row_both = world_repeat_places()
        case = Case(self.tmp, seq, recs, [row_ins, row_both])
        for budget in (("0", "0", "0", "0"), ("4", "1", "1", "6"), ("4", "2", "2", "6")):
            rc, found = case.run(mode="soundness", budget=budget, guide=REPEAT_GUIDE)
            self.assertEqual((rc, gating(found)), (0, []), (budget, found))
        rc, found = case.run(mode="both", budget=("0", "0", "0", "0"), guide=REPEAT_GUIDE)
        self.assertEqual((rc, gating(found)), (0, []), found)

    def test_a_carrier_listed_for_its_target_at_the_other_place_is_caught(self):
        """CRISPRme listing S03 and S04 on the insertion's row (whose bases their haplotype
        holds at the anchor only through the SNV the row does not name), with no row naming
        the SNV, lists them for a target they do not have: soundness fails them on that row
        (variant-missing: their slot holds a variant inside the target the row does not
        name), which it let pass as long as it read their other place as the row's, and
        completeness misses their own target."""
        seq, recs, row_ins, _row_both = world_repeat_places()
        wrong = set_field(row_ins, samples="S01,S02,S03,S04,S05")
        case = Case(self.tmp, seq, recs, [wrong])
        rc, found = case.run(mode="soundness", budget=("0", "0", "0", "0"),
                             guide=REPEAT_GUIDE)
        self.assertEqual(rc, 1, found)
        self.assertEqual({(f[1], f[3]) for f in gating(found)},
                         {("variant-missing", "S03,S04")}, found)
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"),
                             guide=REPEAT_GUIDE)
        self.assertEqual({f[3] for f in gating(found) if f[1] == "missing-row"},
                         {"S03", "S04"}, found)

    def test_a_named_snv_outside_the_reading_is_no_listing_either(self):
        """The row's named SNV on the slot, but not in the slot's reading of the row's bases:
        the run of G is one base short of the target (an A after ten G), and an SNV A>G at
        that base (s1) completes it, the target then read from the insertion's 4th added base
        (the row of the insertion and s1: S01). The SNV of the third unit's T (s2) makes the
        target from the 1st added base, which ends before s1's base (the row of the
        insertion and s2: S03). S07 holds the insertion, s1 and s2 on one slot: it reads the
        bases from the 1st added base through s2, not through s1 (whose base its target does
        not reach), and from the 4th added base not at all (s2 is in the way). S07 is on the
        row naming s2, and not on the row naming s1, which the check demanded it on for the
        s1 its slot holds elsewhere (carrier-missing)."""
        seq = LONG + PRE + "A" + "GGTGGTGGT" + "G" * 10 + "A" + POST
        anchor0 = len(LONG) + len(PRE)
        assert seq[anchor0 + 9] == "T" and seq[anchor0 + 20] == "A"
        ins = (anchor0 + 1, "A", "AGGTGGT",
               {"S01": "1|0", "S02": "0|1", "S03": "1|0", "S07": "1|0"})
        s1 = (anchor0 + 21, "A", "G", {"S01": "1|0", "S07": "1|0"})
        s2 = (anchor0 + 10, "T", "G", {"S03": "1|0", "S07": "1|0"})
        ref = seq[anchor0:anchor0 + 23]

        def row(var, samples):
            return row_fields("chr21", anchor0, "+", REPEAT_GUIDE + "NNN", ref, REPEAT_TARGET,
                              "GGG", 0, 0, "X", ",".join(tok("chr21", v) for v in var),
                              samples, spacer=REPEAT_GUIDE + "NNN")
        rows = [row([ins, s1], "S01"), row([ins, s2], "S03,S07")]
        case = Case(self.tmp, seq, [ins, s1, s2], rows)
        for mode in ("soundness", "both"):
            rc, found = case.run(mode=mode, budget=("0", "0", "0", "0"), guide=REPEAT_GUIDE)
            self.assertEqual((rc, gating(found)), (0, []), (mode, found))
        # S07 listed on the row naming s1 is listed for a target it does not have
        case = Case(self.tmp, seq, [ins, s1, s2], [row([ins, s1], "S01,S07"), rows[1]])
        rc, found = case.run(mode="soundness", budget=("0", "0", "0", "0"),
                             guide=REPEAT_GUIDE)
        self.assertEqual(rc, 1, found)
        self.assertEqual({f[3] for f in gating(found)}, {"S07"}, found)


def code_lines(path):
    """[(line number, text)] of a Python file's code lines: comments, docstrings (a string
    that is a statement) and blank lines left out, whitespace collapsed."""
    import tokenize
    with open(path, encoding="utf-8", errors="replace") as fh:
        src = fh.read()
    skip = set()
    prev = None
    try:
        for t in tokenize.generate_tokens(io.StringIO(src).readline):
            if t.type == tokenize.STRING and (prev is None or prev.type in (
                    tokenize.INDENT, tokenize.DEDENT, tokenize.NEWLINE)):
                skip.update(range(t.start[0], t.end[0] + 1))
            if t.type == tokenize.COMMENT:
                skip.add(("c", t.start[0], t.start[1]))
            if t.type not in (tokenize.NL, tokenize.COMMENT):
                prev = t
    except (tokenize.TokenError, IndentationError, SyntaxError):
        pass   # a file Python cannot read: every line counts as code from where it stopped
    cut = {k[1]: k[2] for k in skip if isinstance(k, tuple)}
    out = []
    for i, line in enumerate(src.splitlines(), 1):
        if i in skip:
            continue
        text = " ".join(line[:cut.get(i, len(line))].split())
        if text:
            out.append((i, text))
    return out


def shared_runs(a, b, k):
    """The runs of ``k`` or more consecutive lines of ``a`` that ``b`` holds consecutively
    too (code_lines), a run of imports aside: [(length, line in a, line in b, texts)]."""
    at = {}
    for j, (_n, text) in enumerate(b):
        at.setdefault(text, []).append(j)
    out = []
    done = set()
    for i, (_n, text) in enumerate(a):
        for j in at.get(text, ()):
            if (i, j) in done:
                continue
            n = 0
            while i + n < len(a) and j + n < len(b) and a[i + n][1] == b[j + n][1]:
                done.add((i + n, j + n))
                n += 1
            texts = [x for _l, x in a[i:i + n]]
            if n >= k and not all(x.startswith(("import ", "from ")) for x in texts):
                out.append((n, a[i][0], b[j][0], texts))
    return out


if __name__ == "__main__":
    unittest.main(verbosity=2)
