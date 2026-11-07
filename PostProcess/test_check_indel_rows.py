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
the worse of two alignments of one site (fewest block or CFD block), a cluster dropped
beside another cluster's row, a row placed past an insertion's anchor, two targets inside
one insertion folded into one, a long deletion's far junction, a second indel that keeps
the target's bases, another guide's row -- and assert that each one is reported, with its
class. The truth is the VCF: a slot whose alleles overlap
gives no sequence and is counted (limit-overlapping-records), never judged by CRISPRme's
rules; those rules (a split multi-allelic site, two indels in one target, a second indel
read up to the target's last base, the PUTATIVE carriers of an MNV) only turn a
disagreement into a non-gating limit, and a tie their PUTATIVE rule leaves is a limit line,
not silence. The PUTATIVE truth is the rule on the whole alignment's CFD (a product of 0
ties every choice), after what a record's known slot says. The driver must never pass
silently (a chromosome no VCF covers, one FASTA of several contigs, an indexed VCF read
once per record, a sites-only VCF, a short table line, a window holding N), and the
checker must hold no CRISPRme code.

Everything is built here: a small reference, a phased VCF (plain gzip, or BGZF with a
tabix index written by this file's own writer), and rows in the column layout of
``integrated_results.tsv``. No CRISPRme module is imported, as the checker imports none.
"""

import gzip
import os
import shutil
import struct
import sys
import tempfile
import unittest
import zlib

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
    for pos1, ref, alt, gts in records:
        row = [chrom, str(pos1), ".", ref, alt, ".", "PASS", "AC=1", "GT"]
        row += [gts.get(s, "0|0") for s in samples]
        body.append("\t".join(row))
    return "\n".join(body) + "\n"


def write_vcf(path, chrom, records, samples=SAMPLES):
    """``records``: (pos1, ref, alt, {sample: gt}), written as a plain gzip stream."""
    with gzip.open(path, "wt") as fh:
        fh.write(vcf_text(chrom, records, samples))


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

    def run(self, mode="soundness", extra=(), budget=("4", "1", "1", "6"), guide=GUIDE):
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
    def test_a_snp_row_carrier_with_an_indel_in_the_target_is_the_documented_limit(self):
        """--snp-rows report: a SNP row lists carriers from SNP genotypes only, so a
        carrier whose haplotype holding the SNP also holds an indel inside the target is
        listed (CHANGELOG 2.8.1, Known limits). Reported, never gating."""
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
        self.assertIn("limit-snp-row-indel", classes(found), found)
        self.assertTrue(any("S02" in f[3] for f in found), found)


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
        with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
            fh.write("no-alt-alignments\n")
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
        the anchor twice, so there is no sequence to judge by; the carrier is listed under
        limit-overlapping-records, never passed over in silence."""
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
        self.assertIn("limit-overlapping-records", classes(found), found)
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
        with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
            fh.write("no-alt-alignments\n")
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
        with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
            fh.write("no-alt-alignments\n")
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

        def alone():
            return [102]
        self.assertTrue(rep.locus(GUIDE, "chr21", 102, d(3, 0.2), alone)[0])
        self.assertFalse(rep.locus(GUIDE, "chr21", 102, d(1, 0.2), alone)[0])  # fewer edits
        self.assertFalse(rep.locus(GUIDE, "chr21", 102, d(3, 0.5), alone)[0])  # higher CFD
        # F2 (review 3): a row past --merge bp is the cluster's only through alignments
        # that chain to it, never as the nearest row
        self.assertFalse(rep.locus(GUIDE, "chr21", 120, d(3, 0.2), lambda: [120])[0])
        self.assertTrue(rep.locus(GUIDE, "chr21", 120, d(3, 0.2),
                                  lambda: [104, 107, 110, 113, 116, 118, 120])[0])
        # N3: another guide's row is not this guide's cluster
        self.assertFalse(rep.locus("TTTTAAAATTTTAAAATTTT", "chr21", 102, d(3, 0.2),
                                   alone)[0])

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
        """F1: no import of CRISPRme's variant or indel modules, and none of the helper
        names of PostProcess/indel_snp_cis.py (a regression guard against copying)."""
        src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                                "scripts", "check_indel_rows.py")).read()
        for mod in ("indel_snp_cis", "analisi_indels", "overlay_indel", "resultIntegrator",
                    "remove_contiguous_samples", "observed_haplotypes", "cfd_score",
                    "bulge_cfd_penalty"):
            self.assertNotRegex(src, r"(?m)^\s*(from|import)\s+%s" % mod)
        for name in ("kept_bases", "shared_tail", "changes_span", "target_ref_span",
                     "holds_indel", "local_haplotype", "target_pairs", "biallelic_site",
                     "resolve_carriers", "worst_case_bases", "_for_alt", "reference_start",
                     "window_place", "holds_target_bases"):
            self.assertNotRegex(src, r"\bdef %s\b" % name)


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
            with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
                fh.write("no-alt-alignments\n")
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
            with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
                fh.write("no-alt-alignments\n")
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
        with open(os.path.join(case.res, ".alt_mode"), "w") as fh:
            fh.write("no-alt-alignments\n")
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
        row's strings, target and variants are checked, the carriers are not."""
        case, info = world_deletion(self.tmp)
        rec = info["rec"]
        with gzip.open(case.vcf, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr21>\n#CHROM\tPOS\tID\tREF\tALT"
                     "\tQUAL\tFILTER\tINFO\nchr21\t%d\t.\t%s\t%s\t.\tPASS\tAF=0.1\n"
                     % (rec[0], rec[1], rec[2]))
        for samples in ("NA", "S01,S02"):
            write_results(case.res, [set_field(info["rows"][0], samples=samples)])
            rc, found = case.run(mode="both")
            self.assertEqual((rc, gating(found)), (0, []), found)
            self.assertIn("limit-no-genotypes", classes(found))
        write_results(case.res, [set_field(info["rows"][0], start=info["start"] + 1)])
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
        self.assertNotIn("budget-exceeded", open_findings(
            C.main(argv + ["--max-total-edits", "6"]), self.tmp))


def open_findings(_rc, tmp):
    with open(os.path.join(tmp, "o.tsv")) as fh:
        return fh.read()


def argparse_ns():
    import argparse
    return argparse.Namespace(mm=4, bdna=1, brna=1, max_total_edits=6)


# --------------------------------------------------------------------------- #
# the driver never passes silently
# --------------------------------------------------------------------------- #
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
        self.assertEqual(rc, 0, open(out).read())

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


if __name__ == "__main__":
    unittest.main(verbosity=2)
