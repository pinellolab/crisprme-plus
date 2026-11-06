#!/usr/bin/env python3
"""Unit tests for ``scripts/check_indel_rows.py``, the independent real-data check of
CRISPRme's indel (and MNV) off-target rows.

The checker must CATCH a wrong row and PASS a right one, so the suite is built both
ways: hand-built true rows that must raise nothing gating, and mutation tests that
corrupt one field of such a row -- the start shifted by 1 and by an insertion's offset,
one target base changed, a carrier swapped or added or dropped, a named variant dropped,
a cross-haplotype SNP added, the aligned spacer rewritten, a phantom of the shape Manuel
Tognon reported (an insertion's only carrier listed under SNP alleles only other
individuals carry), an MNV row and an unphased row each with a phantom base -- and
assert that each one is reported, with its class. CRISPRme's documented limits (a split
multi-allelic site, two indels in one target, a complex allele split into overlapping
records, a phasing artefact, a tandem-repeat second indel) must come out under their own
non-gating class, never as a failure, and the driver must never pass silently (a
chromosome no VCF covers, one FASTA of several contigs, an indexed VCF).

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
        self.assertIn("carrier-without-allele", classes(found), found)

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
        self.assertIn("limit-complex-allele", classes(found), found)
        rc, found = case.run(mode="completeness", budget=("0", "0", "0", "0"))
        self.assertEqual((rc, gating(found)), (0, []), found)
        self.assertTrue(any(f[1] == "limit-complex-allele" and "S02" in f[3]
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
        hap, pl, conf = C.apply_keys(seq, [(3, "GTACG", "G")], 5, 12)
        self.assertEqual(conf, [])
        self.assertEqual(hap, seq[7:12])
        self.assertEqual(pl[0], (7, 0, 0))

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
        """The merge folds rows of one Cluster Position, variant list and carrier set, and
        the carrier set is CRISPRme's: S02 holds an earlier record of a split site at the
        target's first base, which CRISPRme cannot see, so CRISPRme lists it for the
        perfect target as for the RNA-bulge one, and the two fold into the row reported.
        (chr1:213641941 of the genome-wide README run is this case.)"""
        _case, info = world_deletion(self.tmp)
        seq, start = info["seq"], info["start"]
        first = start + 1  # 1-based POS of the target's first base
        recs = [info["rec"], (first, seq[start], "C", {"S02": "0|1"}),
                (first, seq[start], "G", {"S03": "1|0"})]
        rows = [row_fields("chr21", start + 1, "+", GUIDE + "NNN",
                           "-" + seq[start + 1:start + 23], "-" + GUIDE[1:] + "TGG", "TGG",
                           0, 1, "RNA", info["var"], "S01,S02")]
        case = Case(self.tmp, seq, recs, rows)
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
        hap, pl, conf = C.apply_keys(seq, [(5, "A", "ATTT")], 0, 12, (5, "A", "ATTT"))
        self.assertEqual(hap, seq[:5] + "TTT" + seq[5:12])
        self.assertEqual(pl[4], (4, 0, 0))
        self.assertEqual(pl[5:8], [(4, 2, 0), (4, 2, 1), (4, 2, 2)])
        self.assertEqual(pl[8], (5, 0, 0))
        self.assertEqual(conf, [])

    def test_places_of_a_deletion_an_mnv_and_a_shared_tail(self):
        seq = "ACGT" * 10
        hap, pl, _c = C.apply_keys(seq, [(5, "ACG", "A")], 0, 12)
        self.assertEqual(hap, seq[:5] + seq[7:12])
        self.assertEqual(pl[4:6], [(4, 0, 0), (7, 0, 0)])
        hap, pl, _c = C.apply_keys(seq, [(5, "AC", "TG")], 0, 12)
        self.assertEqual((hap[4:6], pl[4:6]), ("TG", [(4, 0, 0), (5, 0, 0)]))
        # GCA>GA at POS 3: the shared A keeps its own coordinate
        hap, pl, _c = C.apply_keys(seq, [(3, "GTA", "GA")], 0, 8)
        self.assertEqual(hap, "ACGACGT")
        self.assertEqual(pl[2:4], [(2, 0, 0), (4, 0, 0)])

    def test_overlapping_records_conflict(self):
        seq = "ACGT" * 10
        _h, _p, conf = C.apply_keys(seq, [(5, "ACG", "A"), (6, "CG", "C")], 0, 12)
        self.assertEqual(conf, [(6, "CG", "C")])

    def test_the_junction_rule(self):
        # GT>G at POS 100: the junction is between 99 and 101
        d = (100, "GT", "G")
        self.assertTrue(C.holds_main([(98, 0, 0), (99, 0, 0), (101, 0, 0)], d))
        self.assertFalse(C.holds_main([(101, 0, 0), (102, 0, 0)], d))
        self.assertFalse(C.holds_main([(97, 0, 0), (98, 0, 0), (99, 0, 0)], d))
        i = (100, "G", "GAAA")
        self.assertTrue(C.holds_main([(99, 2, 2), (100, 0, 0)], i))
        self.assertFalse(C.holds_main([(98, 0, 0), (99, 0, 0)], i))
        m = (100, "AC", "TG")
        self.assertTrue(C.holds_main([(99, 0, 0), (100, 0, 0)], m))
        self.assertFalse(C.holds_main([(101, 0, 0), (102, 0, 0)], m))

    def test_changes_span(self):
        self.assertTrue(C.changes_span((100, "GT", "G"), 99, 101))
        self.assertFalse(C.changes_span((100, "GT", "G"), 102, 120))
        self.assertTrue(C.changes_span((100, "G", "GAAA"), 99, 100))
        self.assertFalse(C.changes_span((100, "G", "GAAA"), 100, 120))
        self.assertFalse(C.changes_span((100, "GTT", "G"), 99, 101, gone=(100, 101)))

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

    def test_readings_from_either_end(self):
        """A haplotype that lost a repeat unit at the target's left edge holds the target
        up to its last base, not from its first."""
        seq = "ACAC" + GUIDE + "TGG" + "TTTT"
        n = len(GUIDE) + 3
        hap, pl, _c = C.apply_keys(seq, [(2, "CAC", "C")], 0, len(seq))
        first, last = (4, 0, 0), (4 + n - 1, 0, 0)
        self.assertTrue(C.reads_as(hap, pl, seq[4:4 + n], first, last, set()))
        self.assertIsNotNone(C.reading(hap, pl, n, 1, first, last))

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
            self.assertEqual(win.shadowed, {(5, "A", "C")})
            self.assertEqual(win.norm((9, "A", "G"), "S01"), (["0", "1"], "phased"))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
