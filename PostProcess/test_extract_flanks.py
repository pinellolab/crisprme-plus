#!/usr/bin/env python3
"""Unit tests for ``PostProcess/extract_flanks.py``, the flanking-sequence writer behind
``crisprme.py extract-flanks``.

The command hands a user the sequence they will order validation oligos from, so the one
thing it must never do is give a sequence that haplotype does not hold. Everything here
is built so the right answer is known BY HAND: a small reference written in this file, a
phased VCF written in this file, and sites placed where the arithmetic is checkable --

  * a site with no variant in its window (one entry, the reference, nothing else),
  * a SNV beside the target and a SNV inside it (the second one makes a haplotype that
    does NOT read the row's protospacer, and must say so),
  * an INSERTION before the target (that haplotype's flank is longer and the target's
    offset moves right by exactly the inserted bases),
  * a DELETION before the target (shorter, and the offset moves left),
  * an insertion INSIDE the target (the carrier reads the row's protospacer exactly; the
    reference does not),
  * two variants in one window (a slot carrying both is a third sequence),
  * a site at a contig start and at a contig end (truncated, said, never a crash),
  * an unphased call and a call with a missing allele (REPORTED, never written as a
    haplotype), and a slot whose records overlap (the VCF gives no sequence there),
  * identical haplotypes deduplicated by sequence, with the carriers of each.

Then a MUTANT for each of those rules: the module is patched to a plausible wrong
implementation (the offset read off the reference rather than the haplotype, an unphased
call read as phased, no deduplication, the window cut by arithmetic rather than by
coordinate, a variant named for a flank it does not change), and the test asserts that
the expectation above then FAILS -- a test that cannot fail is not a test.

The VCFs are written both as plain gzip and as BGZF (this file writes BGZF itself), so
both of the reader's paths are exercised.

Standard library only: no CRISPRme module beside the one under test is imported, and the
shared reader it uses (``scripts/variant_haplotypes.py``) imports none.
"""

import ast
import contextlib
import gzip
import io
import os
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
import zlib
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import extract_flanks as F  # noqa: E402

SAMPLES = ["S1", "S2", "S3", "S4", "S5"]
# a made-up 23-mer (20 spacer + NGG), not any real guide: this suite needs a target the
# arithmetic can be checked against, not a biological one
TARGET = "GTCAGTACCTGATCGGATCAAGG"
TARGET_AT = 400                          # 0-based start on the reference below
CONTIG_LEN = 900


def reference_sequence():
    """A 900 bp contig with the target at TARGET_AT, and no second copy of it (so a
    reading anchored at the row's start is the only reading)."""
    rng = random.Random(20261122)
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(CONTIG_LEN))
        seq = seq[:TARGET_AT] + TARGET + seq[TARGET_AT + len(TARGET):]
        if seq.count(TARGET) == 1 and seq.count(_rc(TARGET)) == 0:
            return seq


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


REF = reference_sequence()


# --------------------------------------------------------------------------- #
# writing the inputs
# --------------------------------------------------------------------------- #
def write_fasta(folder, chrom=None, seq=None, width=60):
    """A per-chromosome FASTA folder, as a CRISPRme genome folder is."""
    os.makedirs(folder, exist_ok=True)
    chrom, seq = chrom or "chr1", seq if seq is not None else REF
    with open(os.path.join(folder, "%s.fa" % chrom), "w") as fh:
        fh.write(">%s\n" % chrom)
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")
    return folder


def vcf_text(records, samples=SAMPLES, chrom="chr1"):
    """A VCF 4.2 text. ``records``: (pos1, ref, alt, {sample: GT}, filter)."""
    out = ["##fileformat=VCFv4.2", "##contig=<ID=%s,length=%d>" % (chrom, CONTIG_LEN),
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for pos1, ref, alt, gts, filt in sorted(records, key=lambda r: r[0]):
        row = [chrom, str(pos1), ".", ref, alt, ".", filt, ".", "GT"]
        row.extend(gts.get(s, "0|0") for s in samples)
        out.append("\t".join(row))
    return "\n".join(out) + "\n"


def write_vcf_gzip(path, records, samples=SAMPLES, chrom="chr1"):
    with gzip.open(path, "wt") as fh:
        fh.write(vcf_text(records, samples, chrom))
    return path


def _bgzf_block(payload):
    """One BGZF member: a deflate stream in a gzip member carrying its size in BC."""
    body = zlib.compressobj(6, zlib.DEFLATED, -15).compress(payload)
    body += zlib.compressobj(6, zlib.DEFLATED, -15).flush()
    comp = zlib.compress(payload, 6)[2:-4]
    bsize = 12 + 6 + len(comp) + 8
    head = (b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6)
            + b"BC" + struct.pack("<H", 2) + struct.pack("<H", bsize - 1))
    return head + comp + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


def write_vcf_bgzf(path, records, samples=SAMPLES, chrom="chr1", per_block=2):
    """The same VCF as BGZF, split over several blocks, with no .tbi: the reader's
    block-index path rather than its plain-gzip one."""
    text = vcf_text(records, samples, chrom).encode()
    lines = text.split(b"\n")
    out = io.BytesIO()
    chunk = []
    for line in lines:
        chunk.append(line)
        if len(chunk) >= per_block:
            out.write(_bgzf_block(b"\n".join(chunk) + b"\n"))
            chunk = []
    if chunk:
        out.write(_bgzf_block(b"\n".join(chunk)))
    out.write(_bgzf_block(b""))   # the EOF block
    with open(path, "wb") as fh:
        fh.write(out.getvalue())
    return path


def write_samples(path, samples=SAMPLES):
    with open(path, "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\n")
        for s in samples:
            fh.write("%s\tPOP\n" % s)
    return path


def write_targets(path, rows, curated=True, block="highest_CFD", spacer=TARGET,
                  ref_aligned="NA"):
    """A results table. ``rows``: (chrom, start0, strand, aligned protospacer+PAM).

    The curated shape is the report's top-N download down to the name it gives the
    guide: ``Aligned_spacer+PAM``, the guide written against the row's alignment. That
    table has no plain ``Spacer+PAM`` column at all, and it is the input of both
    documented examples, so the fixture carries what the real file carries.
    """
    with open(path, "w") as fh:
        if curated:
            fh.write("rank\tChromosome\tPosition\tStrand\tAligned_spacer+PAM\t"
                     "Aligned_protospacer+PAM\tCFD\n")
            for i, (chrom, start, strand, aligned) in enumerate(rows, 1):
                fh.write("%d\t%s\t%d\t%s\t%s\t%s\t0.5000\n"
                         % (i, chrom, start, strand, spacer, aligned))
        else:
            head = ["Spacer+PAM", "Chromosome"]
            for b in F.BLOCKS:
                head += ["Start_coordinate_(%s)" % b, "Strand_(%s)" % b,
                         "Aligned_protospacer+PAM_REF_(%s)" % b,
                         "Aligned_protospacer+PAM_ALT_(%s)" % b]
            fh.write("\t".join(head) + "\n")
            for chrom, start, strand, aligned in rows:
                row = [TARGET, chrom]
                for b in F.BLOCKS:
                    row += ([str(start), strand, ref_aligned, aligned] if b == block
                            else ["0", "+", "NA", "NA"])
                fh.write("\t".join(row) + "\n")
    return path


# --------------------------------------------------------------------------- #
# running
# --------------------------------------------------------------------------- #
class Case(unittest.TestCase):
    """One run of the command over inputs this file writes, read back as a table."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="flanks_")
        self.addCleanup(shutil.rmtree, self.dir, True)

    def path(self, name):
        return os.path.join(self.dir, name)

    def go(self, rows, records=(), window=40, bgzf=False, samples=SAMPLES, chrom="chr1",
           seq=None, curated=True, spacer=TARGET, ref_aligned="NA", **kw):
        genome = write_fasta(self.path("ref"), chrom, REF if seq is None else seq)
        vcf = None
        if records:
            vcf = self.path("v.vcf.gz")
            (write_vcf_bgzf if bgzf else write_vcf_gzip)(vcf, records, samples, chrom)
        table = write_targets(self.path("t.tsv"), rows, curated=curated, spacer=spacer,
                              ref_aligned=ref_aligned)
        out = self.path("out")
        summary = F.run(targets=table, reference=genome, vcfs=[vcf] if vcf else [],
                        samples=[write_samples(self.path("s.txt"), samples)],
                        out=out, window=window, **kw)
        return summary, read_tsv(out + ".tsv"), read_fasta(out + ".fa")


def text_of(path, mode="r"):
    with open(path, mode) as fh:
        return fh.read()


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_fasta(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = ""
            elif name:
                out[name] += line.strip()
    return out


def entries(table, row_id=1):
    return [r for r in table if int(r["row_id"]) == row_id
            and r["haplotype"] != "unresolved"]


def unresolved(table, row_id=1):
    return [r for r in table if int(r["row_id"]) == row_id
            and r["haplotype"] == "unresolved"]


def marked(row, fasta):
    """The bases a row's own offset+length mark inside its own FASTA record: what the
    user would cut out of the file to read the protospacer."""
    seq = fasta[row["fasta_id"]]
    off, n = int(row["target_offset"]), int(row["target_len"])
    return seq[off:off + n]


# the variants used across the cases, all well inside the 40 bp flank of TARGET_AT
SNV_BESIDE = (TARGET_AT - 10, "C")          # 10 bases before the target
SNV_IN_TARGET = (TARGET_AT + 10, "T")       # inside the protospacer
INS_BEFORE = (TARGET_AT - 20, "TTTTT")      # 5 bases added before the target
DEL_BEFORE = (TARGET_AT - 25, 4)            # 4 bases removed before the target


def snv(pos0, base, gts, filt="PASS"):
    ref = REF[pos0]
    assert base != ref
    return (pos0 + 1, ref, base, gts, filt)


def insertion(pos0, added, gts, filt="PASS"):
    """An insertion anchored on the reference base at pos0."""
    return (pos0 + 1, REF[pos0], REF[pos0] + added, gts, filt)


def deletion(pos0, n, gts, filt="PASS"):
    """n reference bases after pos0 removed."""
    return (pos0 + 1, REF[pos0:pos0 + n + 1], REF[pos0], gts, filt)


ROW = [("chr1", TARGET_AT, "+", TARGET)]


# --------------------------------------------------------------------------- #
# the rules
# --------------------------------------------------------------------------- #
class ReferenceOnlySite(Case):
    def test_no_variant_in_the_window_gives_one_entry(self):
        _s, table, fasta = self.go(ROW)
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref"])
        self.assertEqual(got[0]["n_haplotypes"], "1")
        self.assertEqual(fasta[got[0]["fasta_id"]], REF[TARGET_AT - 40:TARGET_AT + 63])
        self.assertEqual(marked(got[0], fasta), TARGET)
        self.assertEqual(got[0]["target_found"], "yes")

    def test_a_variant_outside_the_window_changes_nothing(self):
        far = snv(TARGET_AT - 300, "A" if REF[TARGET_AT - 300] != "A" else "G",
                  {"S1": "1|0"})
        _s, table, _f = self.go(ROW, [far])
        self.assertEqual([r["haplotype"] for r in entries(table)], ["ref"])

    def test_a_record_reaching_into_the_flank_that_changes_no_base_is_not_named(self):
        """An MNV whose REF span crosses the flank's left edge while the base it
        rewrites lies outside: the flank IS the reference there, and naming the record
        for it would read as a difference the user would design around."""
        flo = TARGET_AT - 40
        ref4 = REF[flo - 2:flo + 2]
        alt4 = ("A" if ref4[0] != "A" else "C") + ref4[1:]
        rec = (flo - 1, ref4, alt4, {"S1": "1|0"}, "PASS")
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref"])
        self.assertEqual(fasta[got[0]["fasta_id"]], REF[flo:TARGET_AT + 63])
        self.assertEqual(got[0]["variants"], "")


class SnvSites(Case):
    def test_a_snv_beside_the_target_makes_a_second_sequence(self):
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "1|1"})
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref", "hap1"])
        # S1 slot 1, S3/S4/S5 both slots read the reference; S1 slot 0 and both of S2's
        # carry the SNV
        self.assertEqual(got[1]["carriers"], "S1:0,S2:0,S2:1")
        self.assertEqual(got[1]["n_carriers"], "3")
        self.assertEqual(got[0]["carriers"], "S1:1")
        alt = fasta[got[1]["fasta_id"]]
        want = list(REF[TARGET_AT - 40:TARGET_AT + 63])
        want[SNV_BESIDE[0] - (TARGET_AT - 40)] = SNV_BESIDE[1]
        self.assertEqual(alt, "".join(want))
        # the target itself is untouched, so both haplotypes read the row's protospacer
        self.assertEqual([r["target_found"] for r in got], ["yes", "yes"])
        self.assertEqual(marked(got[1], fasta), TARGET)
        self.assertEqual(got[1]["variants"],
                         "chr1_%d_%s_%s" % (SNV_BESIDE[0] + 1, REF[SNV_BESIDE[0]],
                                            SNV_BESIDE[1]))

    def test_a_snv_inside_the_target_is_not_the_rows_protospacer(self):
        rec = snv(*SNV_IN_TARGET, gts={"S1": "1|0"})
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        self.assertEqual([r["target_found"] for r in got], ["yes", "no"])
        self.assertIn("protospacer-not-on-this-haplotype", got[1]["notes"])
        # the marked stretch is still the same piece of genome, with the ALT in it
        self.assertEqual(len(marked(got[1], fasta)), len(TARGET))
        self.assertNotEqual(marked(got[1], fasta), TARGET)
        self.assertEqual(marked(got[1], fasta)[10], SNV_IN_TARGET[1])


class IndelsMoveTheTarget(Case):
    def test_an_insertion_before_the_target_moves_its_offset_right(self):
        rec = insertion(*INS_BEFORE, gts={"S1": "1|0"})
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        self.assertEqual(got[0]["target_offset"], "40")
        self.assertEqual(got[1]["target_offset"], "45")       # + the 5 bases added
        self.assertEqual(int(got[1]["flank_len"]) - int(got[0]["flank_len"]), 5)
        self.assertEqual(marked(got[1], fasta), TARGET)       # still the row's target
        self.assertEqual(got[1]["target_found"], "yes")

    def test_a_deletion_before_the_target_moves_its_offset_left(self):
        rec = deletion(*DEL_BEFORE, gts={"S1": "1|0"})
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        self.assertEqual(got[0]["target_offset"], "40")
        self.assertEqual(got[1]["target_offset"], "36")       # - the 4 bases removed
        self.assertEqual(int(got[0]["flank_len"]) - int(got[1]["flank_len"]), 4)
        self.assertEqual(marked(got[1], fasta), TARGET)

    def test_an_insertion_inside_the_target_is_read_on_its_carrier(self):
        """An ALT-origin row: the protospacer exists only where the insertion is."""
        added = "GGG"
        rec = insertion(TARGET_AT + 5, added, gts={"S1": "1|0", "S2": "0|1"})
        alt_target = TARGET[:6] + added + TARGET[6:]
        _s, table, fasta = self.go([("chr1", TARGET_AT, "+", alt_target)], [rec])
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"][0]
        self.assertEqual(hap["target_found"], "yes")
        self.assertEqual(marked(hap, fasta), alt_target)
        self.assertEqual(hap["carriers"], "S1:0,S2:1")
        # the target covers len(TARGET) reference bases, not len(alt_target)
        self.assertEqual(int(hap["target_end"]) - int(hap["target_start"]) + 1,
                         len(TARGET))
        # the reference does not hold it, and says so rather than marking it anyway
        ref_row = got[0]
        self.assertEqual(ref_row["target_found"], "no")
        self.assertEqual(marked(ref_row, fasta), TARGET)

    def test_a_deletion_inside_the_target_shortens_what_is_marked(self):
        """A haplotype that does NOT read the row's protospacer still has to mark where
        the site IS on it. A 4 bp deletion inside the target leaves 19 of its bases
        there, so target_len is 19 and the marked bases are those 19 -- not the 23 the
        reference span is wide, which would mark four bases PAST the site and read as
        the protospacer to anyone cutting the FASTA by this offset and length."""
        rec = deletion(TARGET_AT + 5, 4, {"S1": "1|0"})
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"][0]
        self.assertEqual(hap["target_found"], "no")
        self.assertEqual(hap["target_len"], "19")
        self.assertEqual(marked(hap, fasta), TARGET[:6] + TARGET[10:])
        self.assertEqual(hap["target_sequence"], TARGET[:6] + TARGET[10:])
        # the row's own span on the reference is unchanged: 23 coordinates
        self.assertEqual(int(hap["target_end"]) - int(hap["target_start"]) + 1, 23)
        self.assertEqual(got[0]["target_len"], "23")


class ClashesAwayFromTheTarget(Case):
    """A slot whose records overlap has no sequence WHERE THEY DO. That is a fact about
    those bases, not about the row: the target's reference span, and so every coordinate
    the row prints, must be the same whether or not such a pair sits elsewhere in the
    window read around the target."""

    # the row's protospacer: an A inserted after the 10th base of the target, so it is
    # 23 bases long over 22 reference coordinates
    PROBE = TARGET[:10] + "A" + TARGET[10:22]
    CLASH_AT = TARGET_AT + 40        # 18 bases past the protospacer's end

    def records(self, clash):
        out = [insertion(TARGET_AT + 9, "A", {"S1": "1|0"})]
        if clash:
            # two records whose REF spans meet, both on S1's slot 0
            at = self.CLASH_AT
            out += [(at + 1, REF[at:at + 3], REF[at], {"S1": "1|0"}, "PASS"),
                    (at + 2, REF[at + 1], "A" if REF[at + 1] != "A" else "C",
                     {"S1": "1|0"}, "PASS")]
        return out

    def go_with(self, clash):
        # a 10 bp window: the clash is outside the flank, and inside the wider stretch
        # read to size the target
        _s, table, fasta = self.go([("chr1", TARGET_AT, "+", self.PROBE)],
                                   self.records(clash), window=10, max_haplotypes=0)
        return table, fasta

    def test_the_target_span_is_the_same_with_and_without_it(self):
        plain, _f = self.go_with(False)
        with_clash, _f2 = self.go_with(True)
        keys = ("target_start", "target_end", "flank_start", "flank_end")
        ref_plain = [r for r in entries(plain) if r["haplotype"] == "ref"][0]
        ref_clash = [r for r in entries(with_clash) if r["haplotype"] == "ref"][0]
        self.assertEqual([ref_plain[k] for k in keys], [ref_clash[k] for k in keys])
        # and they are the ones the insertion gives: 22 reference coordinates
        self.assertEqual(int(ref_plain["target_end"]) - int(ref_plain["target_start"]) + 1,
                         22)
        self.assertNotIn("protospacer-not-located", ref_clash["notes"])

    def test_the_haplotype_that_holds_the_site_is_still_written(self):
        table, fasta = self.go_with(True)
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"]
        self.assertEqual([r["target_found"] for r in hap], ["yes"])
        self.assertEqual(hap[0]["carriers"], "S1:0")
        self.assertEqual(marked(hap[0], fasta), self.PROBE)

    def test_a_clash_the_reading_itself_crosses_is_not_read_through(self):
        """The other side of the rule: two records that overlap INSIDE the target's own
        span leave those bases unknown, so nothing is read there and the row says so."""
        at = TARGET_AT + 3
        recs = [(at + 1, REF[at:at + 3], REF[at], {"S1": "1|0"}, "PASS"),
                (at + 2, REF[at + 1], "A" if REF[at + 1] != "A" else "C",
                 {"S1": "1|0"}, "PASS")]
        _s, table, _f = self.go([("chr1", TARGET_AT, "+", self.PROBE)], recs,
                                window=10, max_haplotypes=0)
        ref_row = [r for r in entries(table) if r["haplotype"] == "ref"][0]
        self.assertIn("protospacer-not-located", ref_row["notes"])


class SeveralVariantsInOneWindow(Case):
    def test_a_slot_with_both_is_its_own_sequence(self):
        a = snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "1|0"})
        b = snv(TARGET_AT + 30, "A" if REF[TARGET_AT + 30] != "A" else "C",
                {"S2": "1|0", "S3": "1|0"})
        _s, table, fasta = self.go(ROW, [a, b])
        got = entries(table)
        self.assertEqual(len(got), 4)   # reference, a alone, b alone, a+b
        by_carrier = {r["carriers"]: r for r in got}
        self.assertIn("S2:0", by_carrier)                       # both variants
        self.assertEqual(by_carrier["S2:0"]["variants"].count(","), 1)
        self.assertEqual(by_carrier["S1:0"]["variants"].count(","), 0)
        # every emitted sequence differs from every other
        self.assertEqual(len({fasta[r["fasta_id"]] for r in got}), 4)


class ContigEnds(Case):
    def test_a_target_near_the_contig_start_is_truncated_and_said(self):
        seq = TARGET + REF[len(TARGET):]
        _s, table, fasta = self.go([("chr1", 0, "+", TARGET)], seq=seq)
        got = entries(table)
        self.assertIn("truncated-at-contig-start", got[0]["notes"])
        self.assertEqual(got[0]["flank_start"], "1")
        self.assertEqual(got[0]["target_offset"], "0")
        self.assertEqual(marked(got[0], fasta), TARGET)

    def test_a_target_near_the_contig_end_is_truncated_and_said(self):
        start = CONTIG_LEN - len(TARGET)
        seq = REF[:start] + TARGET
        _s, table, fasta = self.go([("chr1", start, "+", TARGET)], seq=seq)
        got = entries(table)
        self.assertIn("truncated-at-contig-end", got[0]["notes"])
        self.assertEqual(got[0]["flank_end"], str(CONTIG_LEN))
        self.assertEqual(marked(got[0], fasta), TARGET)
        self.assertEqual(int(got[0]["target_offset"]) + len(TARGET),
                         int(got[0]["flank_len"]))


class WindowEdges(Case):
    """A record whose REF span STRADDLES the edge of the flank. The flank is a genomic
    interval, so what goes into it is the part of that record's edit lying INSIDE the
    interval -- not the whole edit (the sequence would then be longer than the interval
    it claims, and every offset past it would be wrong) and not nothing (the bases the
    carrier holds there would read as the reference)."""

    WINDOW = 40

    def straddling(self, at, width=4):
        """A ``width``-base MNV at ``at``, every base changed, so there is no shared
        prefix or suffix to trim: its whole REF span is its edit."""
        ref = REF[at:at + width]
        alt = ref.translate(str.maketrans("ACGT", "CATG"))
        self.assertEqual(len([1 for a, b in zip(ref, alt) if a == b]), 0)
        return (at + 1, ref, alt, {"S1": "1|0"}, "PASS"), alt

    def test_a_record_over_the_right_edge_puts_in_only_the_bases_inside(self):
        fhi = TARGET_AT + len(TARGET) + self.WINDOW          # the flank's end, 0-based
        rec, alt = self.straddling(fhi - 2)                  # 2 bases in, 2 bases out
        _s, table, fasta = self.go(ROW, [rec], window=self.WINDOW)
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"][0]
        want = REF[TARGET_AT - self.WINDOW:fhi - 2] + alt[:2]
        self.assertEqual(fasta[hap["fasta_id"]], want)
        # the interval is 103 reference coordinates, so the sequence is 103 bases: a
        # carrier of an MNV holds no more of them than the reference does
        self.assertEqual(hap["flank_len"], got[0]["flank_len"])
        self.assertEqual(int(hap["flank_len"]), len(want))
        self.assertEqual(marked(hap, fasta), TARGET)

    def test_a_record_over_the_left_edge_puts_in_only_the_bases_inside(self):
        flo = TARGET_AT - self.WINDOW                        # the flank's start
        rec, alt = self.straddling(flo - 2)                  # 2 bases out, 2 bases in
        _s, table, fasta = self.go(ROW, [rec], window=self.WINDOW)
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"][0]
        want = alt[2:] + REF[flo + 2:TARGET_AT + len(TARGET) + self.WINDOW]
        self.assertEqual(fasta[hap["fasta_id"]], want)
        self.assertEqual(hap["flank_len"], got[0]["flank_len"])
        self.assertEqual(marked(hap, fasta), TARGET)


class UnknownGenotypes(Case):
    def test_an_unphased_heterozygous_call_writes_no_haplotype(self):
        rec = snv(*SNV_BESIDE, gts={"S1": "0/1", "S2": "1|0"})
        _s, table, _f = self.go(ROW, [rec])
        self.assertNotIn("S1", ",".join(r["carriers"] for r in entries(table)))
        bad = unresolved(table)
        self.assertEqual(len(bad), 1)
        self.assertEqual(bad[0]["carriers"], "S1")
        self.assertTrue(bad[0]["notes"].startswith("unphased-or-missing-genotype:"))
        self.assertIn("chr1_%d_" % (SNV_BESIDE[0] + 1), bad[0]["notes"])
        self.assertEqual(bad[0]["flank_len"], "")   # no sequence is invented for it

    def test_a_missing_allele_writes_no_haplotype(self):
        rec = snv(*SNV_BESIDE, gts={"S1": ".|1"})
        _s, table, _f = self.go(ROW, [rec])
        self.assertEqual([r["carriers"] for r in unresolved(table)], ["S1"])

    def test_an_unphased_homozygote_is_known_and_is_written(self):
        """1/1 says which slots hold it however it is written: both."""
        rec = snv(*SNV_BESIDE, gts={"S1": "1/1"})
        _s, table, _f = self.go(ROW, [rec])
        self.assertEqual(unresolved(table), [])
        self.assertEqual(entries(table)[1]["carriers"], "S1:0,S1:1")

    def test_overlapping_records_on_one_slot_have_no_sequence(self):
        """A deletion and a SNV on a base it removes, phased onto one slot: the VCF
        describes those bases twice, so that slot is reported, not guessed."""
        d = deletion(TARGET_AT - 25, 4, {"S1": "1|0"})
        s = snv(TARGET_AT - 23, "A" if REF[TARGET_AT - 23] != "A" else "C", {"S1": "1|0"})
        _s, table, _f = self.go(ROW, [d, s])
        bad = unresolved(table)
        self.assertEqual([r["carriers"] for r in bad], ["S1:0"])
        self.assertTrue(bad[0]["notes"].startswith("overlapping-records:"))
        # only the copy the clash is on is lost: S1's other copy is plain reference
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref"])
        self.assertEqual(got[0]["carriers"], "S1:1")
        # ONE sample has a record in this window, on both lists: its written copy and
        # its unresolvable one. n_samples counts samples, so it counts it once.
        self.assertEqual(got[0]["n_samples"], "1")

    def test_a_haploid_call_is_one_slot_and_is_written(self):
        """A male chrX is coded haploid (a bare ``0`` / ``1``). Such a sample has ONE
        copy here, and it is known: it is written, with one slot."""
        rec = snv(*SNV_BESIDE, gts={"S1": "1", "S2": "0"})
        _s, table, _f = self.go(ROW, [rec])
        self.assertEqual(unresolved(table), [])
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref", "hap1"])
        self.assertEqual(got[1]["carriers"], "S1:0")
        self.assertEqual(got[1]["n_carriers"], "1")

    def test_a_sample_haploid_at_one_record_and_diploid_at_another_is_reported(self):
        """Mixed ploidy in one window: which of its calls share a copy is not in the
        VCF, so neither copy is written -- reading the haploid call onto BOTH slots
        would invent a carrier, and a haplotype nobody holds."""
        one = snv(*SNV_BESIDE, gts={"S1": "1", "S2": "1|0"})
        two = snv(TARGET_AT - 15, "A" if REF[TARGET_AT - 15] != "A" else "C",
                  {"S1": "1|0"})
        _s, table, _f = self.go(ROW, [one, two])
        bad = unresolved(table)
        self.assertEqual([r["carriers"] for r in bad], ["S1"])
        got = entries(table)
        self.assertNotIn("S1", ",".join(r["carriers"] for r in got))
        # S2's slot 0 is the only haplotype written, and the reference is the other
        self.assertEqual([r["haplotype"] for r in got], ["ref", "hap1"])
        self.assertEqual(got[1]["carriers"], "S2:0")
        self.assertEqual(got[0]["n_samples"], "2")


class Deduplication(Case):
    def test_identical_haplotypes_are_one_entry_with_every_carrier(self):
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "0|1", "S3": "1|1", "S4": "1|0"})
        _s, table, fasta = self.go(ROW, [rec])
        got = entries(table)
        self.assertEqual(len(got), 2)
        self.assertEqual(got[1]["carriers"], "S1:0,S2:1,S3:0,S3:1,S4:0")
        self.assertEqual(got[1]["n_carriers"], "5")
        self.assertEqual(len(fasta), 2)

    def test_two_variant_sets_that_spell_one_sequence_are_one_entry(self):
        """The same inserted base written at two anchors of a homopolymer: one sequence,
        so one oligo, so one entry -- and the entry says its variant list is a union."""
        run_at = _homopolymer_before_target()
        a = insertion(run_at, REF[run_at], {"S1": "1|0"})
        b = insertion(run_at + 1, REF[run_at + 1], {"S2": "1|0"})
        self.assertEqual(REF[run_at], REF[run_at + 1])
        _s, table, fasta = self.go(ROW, [a, b])
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"]
        self.assertEqual(len(hap), 1)
        self.assertEqual(hap[0]["carriers"], "S1:0,S2:0")
        self.assertIn("one-sequence-from-more-than-one-variant-set", hap[0]["notes"])
        self.assertEqual(hap[0]["variants"].count(","), 1)


def _homopolymer_before_target():
    """A coordinate before the target where two reference bases are the same: inserting
    that base at either of them spells one sequence."""
    for i in range(TARGET_AT - 35, TARGET_AT - 5):
        if REF[i] == REF[i + 1]:
            return i
    raise AssertionError("the test reference has no repeated base before the target")


class TooManyHaplotypes(Case):
    """A 400 bp window over a panel of thousands holds dozens of haplotypes. The cap
    keeps the output readable WITHOUT hiding anything: the true total stays in
    n_haplotypes and the notes say how many were left out."""

    def records(self):
        """Five SNVs before the target, each on a sample of its own: five haplotypes."""
        out = []
        for i, sample in enumerate(SAMPLES):
            at = TARGET_AT - 30 + i * 3
            out.append(snv(at, "A" if REF[at] != "A" else "C", {sample: "1|0"}))
        return out

    def test_zero_writes_every_haplotype(self):
        _s, table, _f = self.go(ROW, self.records(), max_haplotypes=0)
        got = entries(table)
        self.assertEqual(len(got), 6)                 # the reference and five
        self.assertEqual(got[0]["n_haplotypes"], "6")
        self.assertNotIn("haplotypes-not-written", got[0]["notes"])

    def test_the_cap_keeps_the_count_and_says_what_it_left_out(self):
        _s, table, fasta = self.go(ROW, self.records(), max_haplotypes=2)
        got = entries(table)
        self.assertEqual(len(got), 3)                 # the reference and two
        self.assertEqual(got[0]["n_haplotypes"], "6")  # the TRUE total, still
        self.assertIn("haplotypes-not-written:3", got[0]["notes"])
        self.assertIn("3 carrier haplotypes", got[0]["notes"])
        # no note holds the column's own separator, so the column stays splittable
        for row in table:
            for note in row["notes"].split(";"):
                self.assertNotIn(";", note)
        self.assertEqual(len(fasta), 3)

    def test_the_haplotypes_that_hold_the_off_target_come_first(self):
        """A row that exists only on an insertion: its carrier must be hap1, however
        many commoner haplotypes the window has, or the cap could drop the one sequence
        the user needs."""
        added = "GGG"
        alt_target = TARGET[:6] + added + TARGET[6:]
        recs = self.records() + [insertion(TARGET_AT + 5, added, {"S5": "0|1"})]
        _s, table, fasta = self.go([("chr1", TARGET_AT, "+", alt_target)], recs,
                                   max_haplotypes=1)
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref", "hap1"])
        self.assertEqual(got[1]["target_found"], "yes")
        self.assertEqual(got[1]["carriers"], "S5:1")
        self.assertEqual(marked(got[1], fasta), alt_target)

    def test_n_samples_counts_only_the_panel_samples_with_a_record_here(self):
        far = snv(TARGET_AT - 300, "A" if REF[TARGET_AT - 300] != "A" else "G",
                  {"S4": "1|0"})
        near = snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "1|1"})
        _s, table, _f = self.go(ROW, [far, near])
        got = entries(table)
        self.assertEqual(got[0]["n_samples"], "2")    # S1 and S2; not S4, not S3/S5
        self.assertEqual(got[0]["carriers"], "S1:1")  # S1's other copy reads as the ref


class Formats(Case):
    def test_an_integrated_table_is_read_by_its_block(self):
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0"})
        _s, table, _f = self.go(ROW, [rec], curated=False)
        self.assertEqual(len(entries(table)), 2)
        self.assertEqual(entries(table)[0]["guide"], TARGET)

    def test_the_curated_download_names_the_guide_aligned_and_it_is_still_read(self):
        """The report's top-N download -- the documented input -- has no plain
        ``Spacer+PAM`` column: the guide is ``Aligned_spacer+PAM``, written against the
        row's alignment, so a DNA-bulge row carries the alignment's gap in it. The gap
        belongs to the row, not to the guide, so the column holds the same guide on
        every line of a file and several guides' files merge on it."""
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0"})
        _s, table, _f = self.go(ROW, [rec])
        self.assertTrue(table)
        for row in table:
            self.assertEqual(row["guide"], TARGET)
        # the same row with the alignment's gap in the guide cell
        _s2, dashed, _f2 = self.go(ROW, [rec], spacer=TARGET[:4] + "-" + TARGET[4:])
        self.assertEqual([r["guide"] for r in dashed], [TARGET] * len(dashed))

    def test_an_integrated_row_is_read_on_its_alt_alignment(self):
        """An integrated row carries BOTH alignments. The ALT one is the sequence the
        site is, where a variant makes it -- 398 of the 1,000 rows of a real shortlist
        have both -- so it is the one read; the REF one stands in only where the row
        has no ALT alignment (a reference off-target, written NA)."""
        added = "GGG"
        alt_target = TARGET[:6] + added + TARGET[6:]
        rec = insertion(TARGET_AT + 5, added, gts={"S1": "1|0"})
        _s, table, fasta = self.go([("chr1", TARGET_AT, "+", alt_target)], [rec],
                                   curated=False, ref_aligned=TARGET)
        got = entries(table)
        hap = [r for r in got if r["haplotype"] != "ref"][0]
        self.assertEqual(hap["target_found"], "yes")        # the ALT alignment is read
        self.assertEqual(marked(hap, fasta), alt_target)
        self.assertEqual(got[0]["target_found"], "no")      # ...and the reference is not
        # the other way round: no ALT alignment, so the REF one is what the row has
        _s2, table2, fasta2 = self.go([("chr1", TARGET_AT, "+", "NA")], [rec],
                                      curated=False, ref_aligned=TARGET)
        ref2 = entries(table2)[0]
        self.assertEqual(ref2["target_found"], "yes")
        self.assertEqual(marked(ref2, fasta2), TARGET)

    def test_another_block_of_an_integrated_table_is_taken_by_name(self):
        genome = write_fasta(self.path("ref"))
        t = write_targets(self.path("t.tsv"), ROW, curated=False, block="fewest_mm+b")
        out = self.path("out")
        F.run(targets=t, reference=genome, out=out, block="fewest_mm+b")
        got = entries(read_tsv(out + ".tsv"))
        self.assertEqual(got[0]["target_start"], str(TARGET_AT + 1))
        # the default block of that file says start 0: taking it would be another site
        F.run(targets=t, reference=genome, out=out + "2")
        self.assertEqual(entries(read_tsv(out + "2.tsv"))[0]["target_start"], "1")

    def test_a_zipped_table_is_read_and_an_ambiguous_zip_is_refused(self):
        import zipfile
        plain = write_targets(self.path("t.tsv"), ROW)
        one = self.path("one.zip")
        with zipfile.ZipFile(one, "w") as zf:
            zf.write(plain, "integrated_results.tsv")
        genome = write_fasta(self.path("ref"))
        F.run(targets=one, reference=genome, out=self.path("z"), window=40)
        self.assertEqual(entries(read_tsv(self.path("z.tsv")))[0]["target_start"],
                         str(TARGET_AT + 1))
        many = self.path("many.zip")
        with zipfile.ZipFile(many, "w") as zf:
            zf.write(plain, "data/top1000.tsv")
            zf.write(plain, "data/mmb_le_3.tsv")
        with self.assertRaises(F.FlankError) as e:
            F.run(targets=many, reference=genome, out=self.path("z2"))
        self.assertIn("data/top1000.tsv", str(e.exception))
        self.assertIn("unzip it", str(e.exception))

    def test_the_bgzf_reader_path_gives_the_same_answer(self):
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "1|1"})
        _s, plain, _f = self.go(ROW, [rec], bgzf=False)
        _s2, bgz, _f2 = self.go(ROW, [rec], bgzf=True)
        self.assertEqual(plain, bgz)

    def test_a_minus_strand_row_is_written_on_the_plus_strand(self):
        row = [("chr1", TARGET_AT, "-", _rc(TARGET))]
        _s, table, fasta = self.go(row)
        got = entries(table)
        self.assertEqual(marked(got[0], fasta), TARGET)             # plus strand
        self.assertEqual(got[0]["protospacer_on_strand"], _rc(TARGET))
        self.assertEqual(got[0]["target_found"], "yes")


class Determinism(Case):
    def test_two_runs_write_the_same_bytes(self):
        rec = [snv(*SNV_BESIDE, gts={"S3": "1|0", "S1": "0|1"}),
               insertion(*INS_BEFORE, gts={"S2": "1|0", "S5": "1|1"})]
        genome = write_fasta(self.path("ref"))
        vcf = write_vcf_gzip(self.path("v.vcf.gz"), rec)
        t = write_targets(self.path("t.tsv"), ROW)
        s = write_samples(self.path("s.txt"))
        first = {}
        for name in ("a", "b"):
            F.run(targets=t, reference=genome, vcfs=[vcf], samples=[s],
                  out=self.path(name))
            for ext in (".fa", ".tsv"):
                first.setdefault(ext, []).append(text_of(self.path(name) + ext, "rb"))
        self.assertEqual(first[".fa"][0], first[".fa"][1])
        self.assertEqual(first[".tsv"][0], first[".tsv"][1])

    def test_haplotypes_come_by_carrier_count_then_sequence(self):
        a = snv(*SNV_BESIDE, gts={"S1": "1|0"})
        b = snv(TARGET_AT + 30, "A" if REF[TARGET_AT + 30] != "A" else "C",
                {"S2": "1|1", "S3": "1|1"})
        _s, table, _f = self.go(ROW, [a, b])
        got = [r for r in entries(table) if r["haplotype"] != "ref"]
        self.assertEqual([int(r["n_carriers"]) for r in got], [4, 1])
        self.assertEqual([r["haplotype"] for r in got], ["hap1", "hap2"])


class Anonymised(Case):
    def test_no_sample_identity_is_written(self):
        rec = [snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "0/1"})]
        _s, table, _f = self.go(ROW, rec, anonymise=True)
        for row in table:
            self.assertEqual(row["carriers"], "")
        text = text_of(self.path("out.tsv")) + text_of(self.path("out.fa"))
        for s in SAMPLES:
            self.assertNotIn(s, text)
        # the counts stay: they are what a report needs, and name nobody
        # (S1 is 1|0: its slot 1 is the reference, its slot 0 the SNV; S2 is unphased)
        self.assertEqual([r["n_carriers"] for r in entries(table)], ["1", "1"])

    def test_a_fasta_name_never_carries_a_sample_even_without_the_flag(self):
        rec = [snv(*SNV_BESIDE, gts={"S1": "1|0"})]
        _s, _t, fasta = self.go(ROW, rec)
        self.assertEqual(sorted(fasta), ["crisprme_row1_hap1", "crisprme_row1_ref"])


class FilteredRecords(Case):
    def test_a_record_outside_pass_is_left_out_and_counted(self):
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0"}, filt="LowQual")
        _s, table, _f = self.go(ROW, [rec])
        self.assertEqual([r["haplotype"] for r in entries(table)], ["ref"])
        self.assertIn("filtered-records-in-window:1", entries(table)[0]["notes"])

    def test_include_filtered_puts_it_back_and_marks_it(self):
        """With the flag ON those bases are IN the sequence the user orders oligos
        against, which is the case that most needs saying: the note says they were
        written (not merely nearby), and the record carries its FILTER in the variants
        column, the one place the entry names its records."""
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0"}, filt="LowQual")
        _s, table, _f = self.go(ROW, [rec], include_filtered=True)
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref", "hap1"])
        self.assertIn("filtered-records-written:1", got[0]["notes"])
        self.assertNotIn("filtered-records-in-window", got[0]["notes"])
        token = "chr1_%d_%s_%s|LowQual" % (SNV_BESIDE[0] + 1, REF[SNV_BESIDE[0]],
                                           SNV_BESIDE[1])
        self.assertEqual(got[1]["variants"], token)

    def test_a_pass_record_carries_no_filter_in_its_token(self):
        """The marking is for the non-PASS ones only: a PASS record reads as before,
        so the column stays joinable with the search's own variant tokens."""
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0"})
        _s, table, _f = self.go(ROW, [rec], include_filtered=True)
        got = entries(table)
        self.assertEqual(got[1]["variants"],
                         "chr1_%d_%s_%s" % (SNV_BESIDE[0] + 1, REF[SNV_BESIDE[0]],
                                            SNV_BESIDE[1]))
        self.assertNotIn("filtered-records", got[0]["notes"])


class WhatTheReferenceEntryClaims(Case):
    """``reference-of-every-sample-without-a-variant-here`` is a claim about the panel:
    that every one of its samples reads as the reference across this flank. It may be
    made only where the samples were READ. Where no VCF covers the contig -- chrY, an
    _alt, a _random, a chrUn: 46 contigs and 54 rows of one real 1,000-row shortlist --
    or where no VCF was read at all, the sequence is the same and the claim is not."""

    def two_contigs(self, **kw):
        """A run whose table names chr1 (in the VCF) and chr2 (in no VCF)."""
        genome = write_fasta(self.path("ref"), "chr1")
        write_fasta(genome, "chr2")
        vcf = write_vcf_gzip(self.path("v.vcf.gz"),
                             [snv(*SNV_BESIDE, gts={"S1": "1|0"})])
        table = write_targets(self.path("t.tsv"),
                              [("chr1", TARGET_AT, "+", TARGET),
                               ("chr2", TARGET_AT, "+", TARGET)])
        said = []
        F.run(targets=table, reference=genome,
              samples=[write_samples(self.path("s.txt"))], out=self.path("out"),
              window=40, log=said.append, vcfs=[vcf], **kw)
        rows = {r["chrom"]: r for r in read_tsv(self.path("out.tsv"))
                if r["haplotype"] == "ref"}
        return rows, said

    def test_a_contig_no_vcf_covers_is_not_called_invariant(self):
        rows, said = self.two_contigs()
        self.assertIn(F.REF_NOTE_CHECKED, rows["chr1"]["notes"])
        self.assertIn(F.REF_NOTE_NO_CONTIG, rows["chr2"]["notes"])
        self.assertNotIn(F.REF_NOTE_CHECKED, rows["chr2"]["notes"])
        # the two sequences are the same; only what was checked differs
        self.assertEqual(rows["chr1"]["flank_len"], rows["chr2"]["flank_len"])
        self.assertEqual(rows["chr2"]["n_samples"], "0")
        # and the contigs nobody was read for are named once, in the log
        named = [m for m in said if F.REF_NOTE_NO_CONTIG in m]
        self.assertEqual(len(named), 1, said)
        self.assertIn("chr2", named[0])
        self.assertNotIn("chr1,", named[0])

    def test_reference_only_claims_nothing_about_any_sample(self):
        rows, said = self.two_contigs(reference_only=True)
        for chrom, row in rows.items():
            with self.subTest(chrom=chrom):
                self.assertIn(F.REF_NOTE_NO_VCF, row["notes"])
                self.assertNotIn(F.REF_NOTE_CHECKED, row["notes"])
        # no VCF was read, so no contig is "the one contig that was missing"
        self.assertEqual([m for m in said if F.REF_NOTE_NO_CONTIG in m], [])

    def test_a_contig_a_vcf_covers_with_no_variant_here_still_says_so(self):
        """The claim itself is unchanged where it is earned: a covered contig whose
        window holds no record is every sample's reference, and says it."""
        far = snv(TARGET_AT - 300, "A" if REF[TARGET_AT - 300] != "A" else "G",
                  {"S1": "1|0"})
        _s, table, _f = self.go(ROW, [far])
        got = entries(table)
        self.assertEqual([r["haplotype"] for r in got], ["ref"])
        self.assertIn(F.REF_NOTE_CHECKED, got[0]["notes"])


# --------------------------------------------------------------------------- #
# the inputs: resolved from the run, or named
# --------------------------------------------------------------------------- #
class RunRecords(Case):
    def make_run(self, command, params=True, sample_ids=True):
        """A results folder as complete-search leaves one: <workdir>/Results/<job>."""
        work = self.path("work")
        job = os.path.join(work, "Results", "job1")
        os.makedirs(job)
        write_fasta(os.path.join(work, "Genomes", "hg38"))
        os.makedirs(os.path.join(work, "VCFs", "ds1"))
        os.makedirs(os.path.join(work, "samplesIDs"))
        write_vcf_gzip(os.path.join(work, "VCFs", "ds1", "ds1.chr1.vcf.gz"),
                       [snv(*SNV_BESIDE, gts={"S1": "1|0"})])
        write_samples(os.path.join(work, "samplesIDs", "ds1.samplesID.txt"))
        with open(os.path.join(work, "vcf.config.txt"), "w") as fh:
            fh.write("ds1\n")
        with open(os.path.join(work, "samples.config.txt"), "w") as fh:
            fh.write("ds1.samplesID.txt\n")
        with open(os.path.join(job, ".command_line.txt"), "w") as fh:
            fh.write("input_command\t" + command + "\n")
        if params:
            with open(os.path.join(job, ".Params.txt"), "w") as fh:
                fh.write("Genome_selected\thg38\nGenome_idx\tNRR_2_hg38+ds1\n")
        if sample_ids:
            write_samples(os.path.join(job, ".sampleID.txt"))
        shutil.copy(write_targets(self.path("t.tsv"), ROW),
                    os.path.join(job, "g_integrated_results.tsv"))
        return job

    CMD = ("crisprme.py complete-search --genome Genomes/hg38 --vcf vcf.config.txt "
           "--samplesID samples.config.txt --output job1")

    def test_a_shortlist_is_read_with_the_runs_own_inputs(self):
        """--results for the inputs, --targets for the table: the report's top-N read
        against the genome and VCFs the run actually searched."""
        job = self.make_run(self.CMD)
        short = write_targets(self.path("short.tsv"), ROW)
        summary = F.run(results=job, targets=short, out=self.path("out"), window=40)
        self.assertEqual(summary["table"], short)
        self.assertTrue(summary["reference"].endswith(os.path.join("Genomes", "hg38")))
        self.assertEqual(len(summary["vcfs"]), 1)
        self.assertEqual([r["haplotype"] for r in entries(read_tsv(self.path("out.tsv")))],
                         ["ref", "hap1"])

    def test_the_run_gives_its_genome_vcfs_and_samples(self):
        job = self.make_run(self.CMD)
        lines = []
        summary = F.run(results=job, out=self.path("out"), window=40,
                        log=lines.append)
        self.assertTrue(summary["reference"].endswith(os.path.join("Genomes", "hg38")))
        self.assertEqual(len(summary["vcfs"]), 1)
        self.assertTrue(summary["samples"][0].endswith(".sampleID.txt"))
        self.assertEqual(len(entries(read_tsv(self.path("out.tsv")))), 2)
        self.assertTrue(any("the run's --genome" in m for m in lines))

    def test_a_results_folder_without_an_integrated_table_is_named(self):
        job = self.make_run(self.CMD)
        os.remove(os.path.join(job, "g_integrated_results.tsv"))
        with self.assertRaises(F.FlankError) as e:
            F.run(results=job, out=self.path("out"))
        self.assertIn("no integrated results table", str(e.exception))
        self.assertIn("--targets", str(e.exception))

    def test_vcfs_the_run_names_but_that_are_gone_are_named_not_guessed(self):
        job = self.make_run(self.CMD)
        work = self.path("work")
        os.remove(os.path.join(work, "VCFs", "ds1", "ds1.chr1.vcf.gz"))
        with self.assertRaises(F.FlankError) as e:
            F.run(results=job, out=self.path("out"))
        self.assertIn("ds1", str(e.exception))
        self.assertIn("--vcf-dir", str(e.exception))
        self.assertIn("--index-path", str(e.exception))
        # and --reference-only still works, saying only what it can
        summary = F.run(results=job, out=self.path("out"), reference_only=True,
                        window=40)
        self.assertEqual(summary["vcfs"], [])
        self.assertEqual([r["haplotype"] for r in entries(read_tsv(self.path("out.tsv")))],
                         ["ref"])

    def test_a_genome_the_run_names_but_that_is_gone_is_named(self):
        job = self.make_run(self.CMD)
        shutil.rmtree(os.path.join(self.path("work"), "Genomes", "hg38"))
        with self.assertRaises(F.FlankError) as e:
            F.run(results=job, out=self.path("out"))
        self.assertIn("Genomes/hg38", str(e.exception).replace(os.sep, "/"))
        self.assertIn("--reference", str(e.exception))

    def test_a_run_that_records_no_genome_is_named(self):
        job = self.make_run("crisprme.py complete-search --output job1")
        with self.assertRaises(F.FlankError) as e:
            F.run(results=job, out=self.path("out"))
        self.assertIn("does not record its reference genome", str(e.exception))

    def test_a_chromosome_the_reference_lacks_is_named(self):
        with self.assertRaises(F.FlankError) as e:
            self.go([("chrZZ", TARGET_AT, "+", TARGET)])
        self.assertIn("chrZZ", str(e.exception))
        self.assertIn("no contig", str(e.exception))

    def test_two_vcfs_that_share_a_sample_are_refused(self):
        genome = write_fasta(self.path("ref"))
        a = write_vcf_gzip(self.path("a.vcf.gz"), [snv(*SNV_BESIDE, gts={"S1": "1|0"})])
        b = write_vcf_gzip(self.path("b.vcf.gz"),
                           [snv(TARGET_AT + 30, "A" if REF[TARGET_AT + 30] != "A" else "C",
                                {"S1": "1|0"})])
        with self.assertRaises(F.FlankError) as e:
            F.run(targets=write_targets(self.path("t.tsv"), ROW), reference=genome,
                  vcfs=[a, b], samples=[write_samples(self.path("s.txt"))],
                  out=self.path("out"), window=40)
        self.assertIn("two VCFs of chr1", str(e.exception))

    def test_a_samples_list_of_another_dataset_is_refused(self):
        genome = write_fasta(self.path("ref"))
        vcf = write_vcf_gzip(self.path("v.vcf.gz"), [snv(*SNV_BESIDE, gts={"S1": "1|0"})])
        other = write_samples(self.path("s.txt"), ["X1", "X2"])
        with self.assertRaises(F.FlankError) as e:
            F.run(targets=write_targets(self.path("t.tsv"), ROW), reference=genome,
                  vcfs=[vcf], samples=[other], out=self.path("out"))
        self.assertIn("different datasets", str(e.exception))

    def test_a_sites_only_vcf_is_said_to_list_no_carrier(self):
        genome = write_fasta(self.path("ref"))
        vcf = self.path("v.vcf.gz")
        with gzip.open(vcf, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
                     "\tINFO\nchr1\t%d\t.\t%s\tC\t.\tPASS\t.\n"
                     % (SNV_BESIDE[0] + 1, REF[SNV_BESIDE[0]]))
        lines = []
        F.run(targets=write_targets(self.path("t.tsv"), ROW), reference=genome,
              vcfs=[vcf], out=self.path("out"), window=40, log=lines.append)
        self.assertTrue(any("list no carrier" in m for m in lines))
        got = entries(read_tsv(self.path("out.tsv")))
        self.assertEqual([r["haplotype"] for r in got], ["ref"])
        # ...and the reference entry does not claim the panel either: there were no
        # genotypes to read, so no sample was checked
        self.assertIn(F.REF_NOTE_NO_GENOTYPES, got[0]["notes"])
        self.assertNotIn(F.REF_NOTE_CHECKED, got[0]["notes"])

    def test_bad_arguments_are_refused_by_name(self):
        for kw, text in (({"out": None}, "--out"),
                         ({"window": -1}, "--window"),
                         ({"block": "nope"}, "--block")):
            args = {"targets": write_targets(self.path("t.tsv"), ROW),
                    "reference": write_fasta(self.path("ref")), "out": self.path("o")}
            args.update(kw)
            with self.assertRaises(F.FlankError) as e:
                F.run(**args)
            self.assertIn(text, str(e.exception))

    def test_a_start_past_the_contig_is_a_line_not_a_crash(self):
        _s, table, _f = self.go([("chr1", CONTIG_LEN + 50, "+", TARGET)])
        self.assertEqual(table[0]["haplotype"], "unresolved")
        self.assertIn("start-past-the-contig-end", table[0]["notes"])
        self.assertIn(str(CONTIG_LEN), table[0]["notes"])

    def test_an_unreadable_row_comes_back_as_a_line_not_a_crash(self):
        t = self.path("t.tsv")
        with open(t, "w") as fh:
            fh.write("rank\tChromosome\tPosition\tStrand\tAligned_protospacer+PAM\n")
            fh.write("1\tchr1\tnot-a-number\t+\t%s\n" % TARGET)
            fh.write("2\tchr1\t%d\t+\t%s\n" % (TARGET_AT, TARGET))
        F.run(targets=t, reference=write_fasta(self.path("ref")), out=self.path("out"),
              window=40)
        table = read_tsv(self.path("out.tsv"))
        self.assertEqual(table[0]["notes"], "unreadable-start")
        self.assertEqual(table[0]["haplotype"], "unresolved")
        self.assertEqual(len(entries(table, row_id=2)), 1)


class Cli(Case):
    def test_the_module_runs_standalone_and_prints_its_help(self):
        script = os.path.join(HERE, "extract_flanks.py")
        p = subprocess.run([sys.executable, script], capture_output=True, text=True)
        self.assertEqual(p.returncode, 1)
        self.assertIn("extract-flanks", p.stderr)
        self.assertIn("--window", p.stderr)
        t = write_targets(self.path("t.tsv"), ROW)
        ok = subprocess.run(
            [sys.executable, script, "--targets", t, "--reference",
             write_fasta(self.path("ref")), "--out", self.path("out"), "--window", "40"],
            capture_output=True, text=True)
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertTrue(os.path.isfile(self.path("out.fa")))
        bad = subprocess.run([sys.executable, script, "--targets", self.path("nope.tsv"),
                              "--out", self.path("o")], capture_output=True, text=True)
        self.assertEqual(bad.returncode, 2)
        self.assertTrue(bad.stderr.startswith("ERROR: "))
        self.assertNotIn("Traceback", bad.stderr)

    def test_the_crisprme_subcommand_is_wired_and_reads_the_same_inputs(self):
        cli = os.path.join(REPO, "crisprme.py")
        p = subprocess.run([sys.executable, cli, "extract-flanks", "--help"],
                           capture_output=True, text=True, cwd=self.dir)
        self.assertEqual(p.returncode, 1)
        self.assertIn("extract-flanks", p.stderr)
        self.assertIn("--window", p.stderr)
        general = subprocess.run([sys.executable, cli], capture_output=True, text=True,
                                 cwd=self.dir)
        self.assertIn("crisprme.py extract-flanks", general.stderr)
        rec = write_vcf_gzip(self.path("v.vcf.gz"),
                             [snv(*SNV_BESIDE, gts={"S1": "1|0"})])
        ok = subprocess.run(
            [sys.executable, cli, "extract-flanks",
             "--targets", write_targets(self.path("t.tsv"), ROW),
             "--reference", write_fasta(self.path("ref")), "--vcf", rec,
             "--samples", write_samples(self.path("s.txt")),
             "--out", self.path("out"), "--window", "40"],
            capture_output=True, text=True, cwd=self.dir)
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertEqual([r["haplotype"] for r in entries(read_tsv(self.path("out.tsv")))],
                         ["ref", "hap1"])
        for flag, text in (("--window", "nope"), ("--block", "nope"),
                           ("--max-rows", "nope")):
            bad = subprocess.run(
                [sys.executable, cli, "extract-flanks", "--targets",
                 self.path("t.tsv"), "--out", self.path("o"), flag, text],
                capture_output=True, text=True, cwd=self.dir)
            self.assertEqual(bad.returncode, 1, flag)
            self.assertTrue(bad.stderr.startswith("Error: "), bad.stderr)
            self.assertIn(flag, bad.stderr)
        none = subprocess.run(
            [sys.executable, cli, "extract-flanks", "--out", self.path("o")],
            capture_output=True, text=True, cwd=self.dir)
        self.assertEqual(none.returncode, 1)
        self.assertIn("--targets", none.stderr)


# --------------------------------------------------------------------------- #
# the shared reader, and the checker it came from
# --------------------------------------------------------------------------- #
class SharedWithTheChecker(unittest.TestCase):
    """The reconstruction is ONE piece of code: the release gate's notion of a haplotype
    and the user's validation sequence must not drift apart. And it must stay what made
    the checker an independent second opinion."""

    MODULE = os.path.join(REPO, "scripts", "variant_haplotypes.py")
    RECORDS = os.path.join(REPO, "scripts", "run_records.py")

    def test_the_tool_and_the_checker_read_the_same_module(self):
        self.assertTrue(os.path.isfile(self.MODULE))
        self.assertTrue(os.path.isfile(self.RECORDS))
        src = text_of(os.path.join(REPO, "scripts", "check_indel_rows.py"))
        self.assertIn("from variant_haplotypes import", src)
        tool = text_of(os.path.join(HERE, "extract_flanks.py"))
        self.assertIn("from variant_haplotypes import", tool)

    def test_the_shared_modules_import_only_the_standard_library(self):
        for path, want in ((self.MODULE, {"bisect", "gzip", "os", "struct", "zlib"}),
                           (self.RECORDS, {"gzip", "os"})):
            tree = ast.parse(text_of(path))
            names = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names.update(a.name.split(".")[0] for a in node.names)
                elif isinstance(node, ast.ImportFrom):
                    names.add((node.module or "").split(".")[0])
            self.assertEqual(names, want, path)

    def test_the_checker_still_imports_no_crisprme_module(self):
        tree = ast.parse(text_of(os.path.join(REPO, "scripts",
                                              "check_indel_rows.py")))
        names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom):
                names.add((node.module or "").split(".")[0])
        self.assertEqual(names - set(sys.stdlib_module_names),
                         {"run_records", "variant_haplotypes"})

    def test_the_checker_still_runs_standalone(self):
        p = subprocess.run([sys.executable,
                            os.path.join(REPO, "scripts", "check_indel_rows.py"),
                            "--help"], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0)
        self.assertIn("--results", p.stdout)

    def test_the_tool_names_the_shared_module_when_it_is_missing(self):
        source = text_of(os.path.join(HERE, "extract_flanks.py"))
        self.assertIn("variant_haplotypes.py", source)
        self.assertIn("run_records.py", source)
        self.assertIn("ships scripts/", source)


# --------------------------------------------------------------------------- #
# the mutants: each rule above, broken on purpose, must break its test
# --------------------------------------------------------------------------- #
@contextlib.contextmanager
def expect_failure(test, why):
    """The block must raise AssertionError: the rule's test is load-bearing."""
    try:
        yield
    except AssertionError:
        return
    raise AssertionError("the mutant (%s) was NOT caught: the rule is untested" % why)


class Mutants(Case):
    def test_an_offset_read_off_the_reference_breaks_the_insertion_rule(self):
        def flat(hap, t_lo, t_hi, probe):
            return 40, len(probe), True        # the reference's own offset, always
        with mock.patch.object(F, "locate_target", flat):
            with expect_failure(self, "the offset ignores the haplotype's own places"):
                IndelsMoveTheTarget(
                    "test_an_insertion_before_the_target_moves_its_offset_right"
                ).debug()

    def test_an_unphased_call_read_as_phased_breaks_the_unphased_rule(self):
        def guess(call, number, h, ploidy):
            answer = F.on_slot(((-1,) * max(ploidy, 1), True), h, call, number)
            return "yes" if answer == "?" else answer
        with mock.patch.object(F, "slot_of", guess):
            with expect_failure(self, "an unphased call is written as a haplotype"):
                UnknownGenotypes(
                    "test_an_unphased_heterozygous_call_writes_no_haplotype").debug()

    def test_no_deduplication_breaks_the_dedup_rule(self):
        original = F.Entry

        class PerSlot(original):
            __slots__ = ("tag",)

        def entry(label, seq):
            e = PerSlot(label, seq)
            e.tag = object()
            return e
        real_site = F.site_entries

        def no_dedup(refseq, chrom, row, variants, window, include_filtered,
                     max_haplotypes=0, **kw):
            out, bad, facts = real_site(refseq, chrom, row, variants, window,
                                        include_filtered, max_haplotypes, **kw)
            split = [out[0]]
            for e in out[1:]:
                for carrier in e.carriers:
                    one = original(e.label, e.seq)
                    one.offset, one.length, one.found = e.offset, e.length, e.found
                    one.carriers = [carrier]
                    one.variants = set(e.variants)
                    split.append(one)
            for i, e in enumerate(split[1:], 1):
                e.label = "hap%d" % i
            return split, bad, facts
        with mock.patch.object(F, "site_entries", no_dedup):
            with expect_failure(self, "every slot gets its own entry"):
                Deduplication(
                    "test_identical_haplotypes_are_one_entry_with_every_carrier").debug()

    def test_sizing_a_flank_by_arithmetic_breaks_the_deletion_rule(self):
        real = F.make_hap

        def padded(refseq, alleles, lo, hi):
            """As if every haplotype were as long as the reference stretch."""
            hap = real(refseq, alleles, lo, hi)
            want = min(hi, len(refseq)) - lo
            hap.seq = (hap.seq + refseq[lo:lo + want])[:want]
            return hap
        with mock.patch.object(F, "make_hap", padded):
            with expect_failure(self, "the flank is sized by arithmetic, not by what "
                                      "the haplotype holds"):
                IndelsMoveTheTarget(
                    "test_a_deletion_before_the_target_moves_its_offset_left").debug()

    def test_naming_a_variant_that_changes_nothing_breaks_the_variant_rule(self):
        with mock.patch.object(F, "touches", lambda a, lo, hi: True):
            with expect_failure(self, "a variant outside the flank is named for it"):
                ReferenceOnlySite(
                    "test_a_record_reaching_into_the_flank_that_changes_no_base_"
                    "is_not_named").debug()

    def test_ignoring_a_clash_breaks_the_overlapping_records_rule(self):
        real = F.make_hap

        def no_clash(refseq, alleles, lo, hi):
            hap = real(refseq, alleles, lo, hi)
            hap.clash = []
            return hap
        with mock.patch.object(F, "make_hap", no_clash):
            with expect_failure(self, "overlapping records are written anyway"):
                UnknownGenotypes(
                    "test_overlapping_records_on_one_slot_have_no_sequence").debug()

    def test_marking_the_protospacer_as_found_everywhere_breaks_the_found_rule(self):
        real = F.locate_target

        def always(hap, t_lo, t_hi, probe):
            off, n, _found = real(hap, t_lo, t_hi, probe)
            return off, n, True
        with mock.patch.object(F, "locate_target", always):
            with expect_failure(self, "every haplotype claims the row's protospacer"):
                SnvSites(
                    "test_a_snv_inside_the_target_is_not_the_rows_protospacer").debug()

    def test_a_cap_that_hides_what_it_dropped_breaks_the_cap_rule(self):
        real_site = F.site_entries

        def quiet(refseq, chrom, row, variants, window, include_filtered,
                  max_haplotypes=0, **kw):
            out, bad, facts = real_site(refseq, chrom, row, variants, window,
                                        include_filtered, max_haplotypes, **kw)
            facts["notes"] = {n for n in facts["notes"]
                              if not n.startswith("haplotypes-not-written")}
            facts["n_haplotypes"] = len(out)
            return out, bad, facts
        with mock.patch.object(F, "site_entries", quiet):
            with expect_failure(self, "the cap drops haplotypes without saying so"):
                TooManyHaplotypes(
                    "test_the_cap_keeps_the_count_and_says_what_it_left_out").debug()

    def test_leaking_a_sample_into_a_fasta_name_breaks_the_anonymise_rule(self):
        with mock.patch.object(F, "fasta_id",
                               lambda row_id, label: "row%d_%s_S1" % (row_id, label)):
            with expect_failure(self, "a sample name reaches a FASTA record name"):
                Anonymised("test_no_sample_identity_is_written").debug()

    def test_dropping_a_slot_for_a_clash_anywhere_breaks_the_target_span_rule(self):
        """The rule before: a slot with ANY overlapping pair in the window read around
        the target was skipped, wherever the pair was."""
        real = F.reads_at

        def blind(hap, coord, probe, al=None):
            return [] if hap.clash else real(hap, coord, probe, al)
        with mock.patch.object(F, "reads_at", blind):
            with expect_failure(self, "a clash far from the target loses the reading"):
                ClashesAwayFromTheTarget(
                    "test_the_target_span_is_the_same_with_and_without_it").debug()

    def test_a_window_that_takes_one_base_too_many_breaks_the_edge_rule(self):
        """``min(stand, hi - edit_lo)`` + 1: a record straddling the right edge then
        puts one base more into the flank than the interval has coordinates."""
        def wide(self, lo, hi):
            width = self.edit_hi - self.edit_lo
            stand = min(width, len(self.body))
            a = max(0, lo - self.edit_lo)
            b = max(a, min(stand, hi - self.edit_lo + 1))
            if len(self.body) > stand and lo <= self.edit_lo + width - 1 < hi:
                b = len(self.body)
            return a, b
        with mock.patch.object(F.Allele, "window", wide):
            with expect_failure(self, "the window takes a base from past its edge"):
                WindowEdges(
                    "test_a_record_over_the_right_edge_puts_in_only_the_bases_inside"
                ).debug()

    def test_a_target_len_from_the_reference_span_breaks_the_deletion_rule(self):
        """``t_hi - t_lo`` instead of the span mapped through the haplotype: on a
        carrier of a deletion inside the target it marks bases past the site."""
        real = F.locate_target

        def flat(hap, t_lo, t_hi, probe):
            off, length, found = real(hap, t_lo, t_hi, probe)
            return (off, length, found) if found else (off, t_hi - t_lo, found)
        with mock.patch.object(F, "locate_target", flat):
            with expect_failure(self, "the target's length is read off the reference"):
                IndelsMoveTheTarget(
                    "test_a_deletion_inside_the_target_shortens_what_is_marked").debug()

    def test_a_haploid_call_read_onto_every_slot_breaks_the_ploidy_rule(self):
        """``all(x == number for x in alle)`` without the ploidy guard: a haploid call
        beside a diploid one then covers both slots, inventing a carrier."""
        real = F.slot_of

        def loose(call, number, h, ploidy):
            alle, _phased = call
            return "yes" if all(x == number for x in alle) else real(call, number, h,
                                                                     ploidy)
        with mock.patch.object(F, "slot_of", loose):
            with expect_failure(self, "a haploid call is read onto both slots"):
                UnknownGenotypes(
                    "test_a_sample_haploid_at_one_record_and_diploid_at_another_is_"
                    "reported").debug()

    def test_preferring_the_ref_alignment_breaks_the_integrated_alignment_rule(self):
        """``cell('aligned_ref') or cell('aligned')``: the row is then read on the
        sequence the reference has, not on the one the variant makes."""
        real = F.table_columns

        def ref_first(header, block):
            cols, shape = real(header, block)
            if cols.get("aligned_ref") is not None:
                cols["aligned"] = cols["aligned_ref"]
            return cols, shape
        with mock.patch.object(F, "table_columns", ref_first):
            with expect_failure(self, "the REF alignment is read instead of the ALT"):
                Formats("test_an_integrated_row_is_read_on_its_alt_alignment").debug()

    def test_a_reference_entry_that_always_claims_the_panel_breaks_the_claim_rule(self):
        """The note before: every reference entry said every sample reads as the
        reference here, including on a contig no VCF covers."""
        with mock.patch.object(F, "REF_NOTE_NO_CONTIG", F.REF_NOTE_CHECKED):
            with expect_failure(self, "an unchecked contig claims the whole panel"):
                WhatTheReferenceEntryClaims(
                    "test_a_contig_no_vcf_covers_is_not_called_invariant").debug()

    def test_counting_a_clashing_sample_twice_breaks_the_n_samples_rule(self):
        """``len(written) + len(unresolved)``: a sample with one copy written and one
        unresolvable is in both, and the column then runs past the panel."""
        real_site = F.site_entries

        def twice(refseq, chrom, row, variants, window, include_filtered,
                  max_haplotypes=0, **kw):
            out, bad, facts = real_site(refseq, chrom, row, variants, window,
                                        include_filtered, max_haplotypes, **kw)
            facts["n_samples"] = facts["n_samples"] + len(
                {w.split(":")[0] for w, _r in bad})
            return out, bad, facts
        with mock.patch.object(F, "site_entries", twice):
            with expect_failure(self, "a sample on both lists is counted twice"):
                UnknownGenotypes(
                    "test_overlapping_records_on_one_slot_have_no_sequence").debug()


if __name__ == "__main__":
    unittest.main(verbosity=2)
