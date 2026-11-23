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
import csv
import gzip
import io
import os
import random
import re
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
# the GUIDE as a results table writes it: the spacer, then one N per PAM base. That run
# of N's is how the PAM's side is read off a table (a trailing run is a 3' PAM, a
# leading one a 5' PAM), and its length with the gaps taken out is the nominal
# spacer+PAM length a fixed-total window is measured against -- so the fixture carries
# what the real file carries.
SPACER = TARGET[:20] + "NNN"             # 20 spacer + a 3-nt 3' PAM: SpCas9
SPACER_5P = "NNNN" + TARGET[:19]         # 4-nt 5' PAM then 19 spacer: Cas12a-shaped
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


def write_targets(path, rows, curated=True, block="highest_CFD", spacer=SPACER,
                  ref_aligned="NA", extra=()):
    """A results table. ``rows``: (chrom, start0, strand, aligned protospacer+PAM).

    The curated shape is the report's top-N download down to the name it gives the
    guide: ``Aligned_spacer+PAM``, the guide written against the row's alignment. That
    table has no plain ``Spacer+PAM`` column at all, and it is the input of both
    documented examples, so the fixture carries what the real file carries.

    ``extra``: (column name, [value per row]) pairs appended to the curated shape, for
    the tests that filter on a score column.
    """
    names = [c for c, _v in extra]
    with open(path, "w") as fh:
        if curated:
            fh.write("\t".join(["rank", "Chromosome", "Position", "Strand",
                                "Aligned_spacer+PAM", "Aligned_protospacer+PAM", "CFD"]
                               + names) + "\n")
            for i, (chrom, start, strand, aligned) in enumerate(rows, 1):
                fh.write("\t".join([str(i), chrom, str(start), strand, spacer, aligned,
                                    "0.5000"]
                                   + [str(v[i - 1]) for _c, v in extra]) + "\n")
        else:
            head = ["Spacer+PAM", "Chromosome"]
            for b in F.BLOCKS:
                head += ["Start_coordinate_(%s)" % b, "Strand_(%s)" % b,
                         "Aligned_spacer+PAM_(%s)" % b,
                         "Aligned_protospacer+PAM_REF_(%s)" % b,
                         "Aligned_protospacer+PAM_ALT_(%s)" % b]
            fh.write("\t".join(head + names) + "\n")
            for i, (chrom, start, strand, aligned) in enumerate(rows, 1):
                row = [spacer, chrom]
                for b in F.BLOCKS:
                    row += ([str(start), strand, spacer, ref_aligned, aligned]
                            if b == block else ["0", "+", spacer, "NA", "NA"])
                fh.write("\t".join(row + [str(v[i - 1]) for _c, v in extra]) + "\n")
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
           seq=None, curated=True, spacer=SPACER, ref_aligned="NA", extra=(), **kw):
        genome = write_fasta(self.path("ref"), chrom, REF if seq is None else seq)
        vcf = None
        if records:
            vcf = self.path("v.vcf.gz")
            (write_vcf_bgzf if bgzf else write_vcf_gzip)(vcf, records, samples, chrom)
        table = write_targets(self.path("t.tsv"), rows, curated=curated, spacer=spacer,
                              ref_aligned=ref_aligned, extra=extra)
        out = self.path("out")
        if window is not None:
            kw["window"] = window
        summary = F.run(targets=table, reference=genome, vcfs=[vcf] if vcf else [],
                        samples=[write_samples(self.path("s.txt"), samples)],
                        out=out, **kw)
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
        self.assertEqual(entries(table)[0]["guide"], SPACER)

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
            self.assertEqual(row["guide"], SPACER)
        # the same row with the alignment's gap in the guide cell
        _s2, dashed, _f2 = self.go(ROW, [rec], spacer=SPACER[:4] + "-" + SPACER[4:])
        self.assertEqual([r["guide"] for r in dashed], [SPACER] * len(dashed))

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

    def test_a_minus_strand_row_is_written_on_the_plus_strand_when_asked(self):
        row = [("chr1", TARGET_AT, "-", _rc(TARGET))]
        _s, table, fasta = self.go(row, orientation="genome")
        got = entries(table)
        self.assertEqual(marked(got[0], fasta), TARGET)             # plus strand
        self.assertEqual(got[0]["protospacer_on_strand"], _rc(TARGET))
        self.assertEqual(got[0]["target_found"], "yes")
        self.assertEqual(got[0]["orientation"], "genome")
        self.assertEqual(got[0]["emitted_strand"], "+")
        self.assertEqual(fasta[got[0]["fasta_id"]],
                         REF[TARGET_AT - 40:TARGET_AT + 63])


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


def _biopython_is_installed():
    try:
        import Bio.Seq  # noqa: F401
    except ImportError:
        return False
    return True


class Cli(Case):
    """The command line, run as a user runs it.

    ``crisprme.py`` imports Biopython at module load (``from Bio.Seq import Seq``) and
    the unit-tests job installs the light deps only -- not Biopython, which is why the
    other suites exec crisprme.py's definitions with Bio stubbed rather than import it
    (``test_bulge_cap._load_crisprme``, ``test_annotation_tmpdirs._load_crisprme``).
    These tests must run the REAL process, so where Biopython is absent they put that
    same stub (``Seq = object``) on the subprocess's path instead of skipping: the CLI
    wiring of a new subcommand is exactly what CI has to hold. ``extract-flanks``
    reaches ``Seq`` nowhere -- only ``complete-search``'s guide extraction does -- and
    where Biopython IS installed the real module is imported and nothing is stubbed.
    """

    def setUp(self):
        super().setUp()
        if _biopython_is_installed():
            return
        stub = os.path.join(self.dir, "_bio_stub")
        os.makedirs(os.path.join(stub, "Bio"))
        open(os.path.join(stub, "Bio", "__init__.py"), "w").close()
        with open(os.path.join(stub, "Bio", "Seq.py"), "w") as fh:
            fh.write("class Seq(object):\n    pass\n")
        path = os.pathsep.join([stub] + ([os.environ["PYTHONPATH"]]
                                         if os.environ.get("PYTHONPATH") else []))
        patch = mock.patch.dict(os.environ, {"PYTHONPATH": path})
        patch.start()
        self.addCleanup(patch.stop)

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

    def test_the_subcommand_passes_the_new_flags_through(self):
        """The anchored window, the column mode and the filter, reached the way a user
        reaches them: through crisprme.py, not through run()."""
        cli = os.path.join(REPO, "crisprme.py")
        import cb_operating_point as cbop
        thr = cbop.sensitive_threshold()
        rows = [("chr1", TARGET_AT, "+", TARGET)] * 2
        table = write_targets(self.path("t.tsv"), rows,
                              extra=[("CRISPR-Bulge", [thr * 10, thr / 10]),
                                     ("CRISPR-Bulge_ref", [thr * 10, thr * 10])])
        ok = subprocess.run(
            [sys.executable, cli, "extract-flanks", "--targets", table,
             "--reference", write_fasta(self.path("ref")),
             "--upstream", "3", "--downstream", "6",
             "--annotate-table", self.path("ann.tsv"), "--alt-list",
             "--filter", "CRISPR_Bulge >= CB_SENSITIVE and "
                         "CRISPR_Bulge_ref >= CB_SENSITIVE",
             "--out", self.path("out")],
            capture_output=True, text=True, cwd=self.dir)
        self.assertEqual(ok.returncode, 0, ok.stderr)
        got = entries(read_tsv(self.path("out.tsv")))
        self.assertEqual([r["row_id"] for r in got], ["1"])
        self.assertEqual(got[0]["flank_len"], "32")
        ann = read_tsv(self.path("ann.tsv"))
        self.assertEqual(len(ann), 2)
        self.assertEqual(ann[0]["flank_ref_len"], "32")
        self.assertEqual(ann[1]["flank_notes"], "filtered-out")
        self.assertIn("flank_alt_all", ann[0])
        # and the presets reach the subcommand's help
        helped = subprocess.run([sys.executable, cli, "extract-flanks", "--help"],
                                capture_output=True, text=True, cwd=self.dir)
        self.assertIn("PRESETS (--preset <name>)", helped.stderr)
        self.assertIn("--upstream", helped.stderr)
        for flag, text in (("--upstream", "nope"), ("--downstream", "-4"),
                           ("--length-policy", "nope"), ("--orientation", "nope"),
                           ("--pam-side", "nope")):
            bad = subprocess.run(
                [sys.executable, cli, "extract-flanks", "--targets", table,
                 "--out", self.path("o"), flag, text],
                capture_output=True, text=True, cwd=self.dir)
            self.assertEqual(bad.returncode, 1, flag)
            self.assertTrue(bad.stderr.startswith("Error: "), bad.stderr)
            self.assertIn(flag, bad.stderr)


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
# the anchored window: the padding is named for the end of the target it sits
# beyond, not for a genomic side
# --------------------------------------------------------------------------- #
# a 5'-PAM row: the protospacer+PAM is still 23 bases of genome at the same place, but
# the guide's N's are at its FRONT, so the two paddings land on the other genomic sides
FIVE_PRIME = dict(spacer=SPACER_5P)


class AnchoredWindows(Case):
    def test_the_default_is_a_32mer_on_a_mismatch_only_site(self):
        """3 bases beyond the protospacer + the 23-base target + 6 beyond the PAM, and
        the target at offset 3. Every number of it is in the output: nobody has to
        count, and nobody has to remember the command line."""
        _s, table, fasta = self.go(ROW, window=None)
        got = entries(table)[0]
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 3:TARGET_AT + 29])
        self.assertEqual(got["flank_len"], "32")
        self.assertEqual(got["target_offset"], "3")
        self.assertEqual(got["target_len"], "23")
        self.assertEqual(marked(got, fasta), TARGET)
        self.assertEqual(got["window_spec"], "up3+down6")
        self.assertEqual((got["pad_upstream"], got["pad_downstream"]), ("3", "6"))
        self.assertEqual(got["nominal_target_len"], "23")
        self.assertEqual(got["target_span_ref"], "23")
        self.assertEqual(got["bulge_type"], "none")
        self.assertEqual(got["length_policy"], "anchored")

    def test_the_two_forms_of_one_window_agree(self):
        """--window W IS --upstream W --downstream W, so the two cannot drift apart."""
        rec = snv(*SNV_BESIDE, gts={"S1": "1|0", "S2": "1|1"})
        _s, sym, _f = self.go(ROW, [rec], window=40)
        _s2, anchored, _f2 = self.go(ROW, [rec], window=None, upstream=40, downstream=40)
        self.assertEqual([r["flank_len"] for r in sym],
                         [r["flank_len"] for r in anchored])
        self.assertEqual([r["target_sequence"] for r in sym],
                         [r["target_sequence"] for r in anchored])
        self.assertEqual(sym[0]["window_spec"], "window40")
        self.assertEqual(anchored[0]["window_spec"], "up40+down40")

    def test_giving_both_forms_is_refused_rather_than_resolved(self):
        with self.assertRaises(F.FlankError) as e:
            self.go(ROW, window=40, upstream=3)
        self.assertIn("--window", str(e.exception))
        self.assertIn("one form or the other", str(e.exception))

    def test_the_fasta_header_says_which_way_round_and_in_what_window(self):
        """The FASTA is used on its own, so it must not need the TSV to be read right:
        a reverse-complemented record that does not say so is a wrong oligo."""
        row = [("chr1", TARGET_AT, "-", _rc(TARGET))]
        _s, table, _f = self.go(row, window=None)
        head = {}
        with open(self.path("out.fa")) as fh:
            for line in fh:
                if line.startswith(">"):
                    head = dict(f.split("=", 1) for f in line.split()[1:])
                    break
        self.assertEqual(head["orientation"], "guide")
        self.assertEqual(head["emitted_strand"], "-")
        self.assertEqual(head["window_spec"], "up3+down6")
        self.assertEqual(head["target_offset"], "3")

    def test_a_padding_of_zero_is_the_target_itself(self):
        _s, table, fasta = self.go(ROW, window=None, upstream=0, downstream=0)
        got = entries(table)[0]
        self.assertEqual(fasta[got["fasta_id"]], TARGET)
        self.assertEqual(got["flank_len"], "23")
        self.assertEqual(got["target_offset"], "0")


class StrandAndPamSide(Case):
    def test_a_minus_strand_row_pads_in_the_guides_orientation(self):
        """The same 3 and 6, read the same way round: on a minus row the 6 bases beyond
        the PAM are at the genomic LEFT, and the emitted sequence is the reverse
        complement of that interval, so the target is still at offset 3."""
        row = [("chr1", TARGET_AT, "-", _rc(TARGET))]
        _s, table, fasta = self.go(row, window=None)
        got = entries(table)[0]
        self.assertEqual(got["emitted_strand"], "-")
        self.assertEqual(got["orientation"], "guide")
        # the genomic interval: 6 to the left of the target, 3 to its right
        self.assertEqual((got["flank_start"], got["flank_end"]),
                         (str(TARGET_AT - 6 + 1), str(TARGET_AT + 23 + 3)))
        self.assertEqual(fasta[got["fasta_id"]],
                         _rc(REF[TARGET_AT - 6:TARGET_AT + 26]))
        self.assertEqual(got["target_offset"], "3")
        self.assertEqual(marked(got, fasta), _rc(TARGET))
        self.assertEqual(got["protospacer_on_strand"], _rc(TARGET))

    def test_a_five_prime_pam_puts_the_paddings_on_the_other_genomic_sides(self):
        """Cas12a reads PAM then protospacer, so the target's two ends swap: the 6
        bases beyond the PAM are now at the genomic LEFT of a plus-strand row, which is
        where the 3 were under a 3' PAM. Read from the guide's own run of N's."""
        _s, table, fasta = self.go(ROW, window=None, **FIVE_PRIME)
        got = entries(table)[0]
        self.assertEqual(got["pam_side"], "5prime")
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 6:TARGET_AT + 26])
        self.assertEqual(got["target_offset"], "6")
        self.assertEqual(marked(got, fasta), TARGET)
        # ...and the 3' case is the mirror of it, at the same site
        _s2, three, f3 = self.go(ROW, window=None)
        self.assertEqual(f3[entries(three)[0]["fasta_id"]],
                         REF[TARGET_AT - 3:TARGET_AT + 29])

    def test_the_pam_side_comes_from_the_pam_file_given(self):
        """A PAM file carries the side in the SIGN of its length: the one unambiguous
        source, and the one the run's own command line names."""
        five = self.path("5p.txt")
        with open(five, "w") as fh:
            fh.write("TTTN" + "N" * 19 + " -4\n")
        self.assertEqual(F.read_pam_file(five), ("TTTN", "5prime", 4))
        three = self.path("3p.txt")
        with open(three, "w") as fh:
            fh.write("N" * 20 + "NRG 3\n")
        self.assertEqual(F.read_pam_file(three), ("NRG", "3prime", 3))
        # a guide with no N's says nothing; the file still does
        _s, table, fasta = self.go(ROW, window=None, spacer=TARGET, pam_file=five)
        got = entries(table)[0]
        self.assertEqual(got["pam_side"], "5prime")
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 6:TARGET_AT + 26])

    def test_an_asymmetric_window_with_no_readable_pam_side_is_refused(self):
        """Never assume SpCas9: with 3 and 6 on a table that does not say which end is
        the PAM, writing the window either way round is a guess, so it is refused and
        the three ways to answer are named."""
        with self.assertRaises(F.FlankError) as e:
            self.go(ROW, window=None, spacer=TARGET)
        self.assertIn("--pam-file", str(e.exception))
        self.assertIn("--pam-side", str(e.exception))
        self.assertIn("--window", str(e.exception))

    def test_pam_side_alone_answers_a_table_that_cannot(self):
        """A guide with no N's says nothing about the side; --pam-side does, with no
        motif needed, and it is what the row then reports."""
        _s, table, fasta = self.go(ROW, window=None, spacer=TARGET,
                                   pam_side="5prime")
        got = entries(table)[0]
        self.assertEqual(got["pam_side"], "5prime")
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 6:TARGET_AT + 26])
        # the column is still read for the length and the bulges: it is only the SIDE
        # that a guide written without its N's cannot give
        self.assertEqual(got["nominal_target_len"], "23")
        self.assertEqual(got["bulge_type"], "none")

    def test_a_symmetric_window_does_not_need_the_pam_side(self):
        """With the same padding at both ends the side changes nothing, so a table that
        cannot say is still served."""
        _s, table, fasta = self.go(ROW, window=40, spacer=TARGET)
        got = entries(table)[0]
        self.assertEqual(got["pam_side"], "unknown")
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 40:TARGET_AT + 63])


class ThePamSideOfACopiedResultsFolder(Case):
    """The branch of ``pam_side_from_spacer`` a results folder without its PAMs/ folder
    takes -- the common case for a run copied or downloaded on its own, because
    .Params.txt records the PAM's MOTIF and not which end of the target it is on. The
    side is then read off the table's own ``Aligned_spacer+PAM`` with the motif's length
    known, and getting it backwards puts both paddings on the opposite genomic side of
    every row, which is the one failure the whole design exists to prevent.
    """

    # a guide whose first position is degenerate (an N) as well as its 3-nt PAM: with
    # the motif's LENGTH known the two ends are still told apart, and without it they
    # are not. The real tables carry such guides, and this is the branch they take.
    LEADING_N = "N" + TARGET[1:20] + "NNN"

    def test_the_side_is_read_off_the_spacer_with_the_motifs_length(self):
        self.assertEqual(F.pam_side_from_spacer(SPACER, 3), (F.PAM_3PRIME, 3))
        self.assertEqual(F.pam_side_from_spacer(SPACER_5P, 4), (F.PAM_5PRIME, 4))
        # the length is what tells a degenerate guide position from the PAM
        self.assertEqual(F.pam_side_from_spacer(self.LEADING_N, 3), (F.PAM_3PRIME, 3))
        self.assertEqual(F.pam_side_from_spacer("NNNN" + TARGET[4:] + "N", 4),
                         (F.PAM_5PRIME, 4))
        # ...and where it cannot tell them apart it says so rather than guess
        self.assertIsNone(F.pam_side_from_spacer("NNN" + TARGET[3:20] + "NNN", 3))
        self.assertIsNone(F.pam_side_from_spacer(TARGET, 3))
        self.assertIsNone(F.pam_side_from_spacer("N" * 23, 3))
        # a DNA bulge writes a gap into the guide; the run of N's is read without it
        self.assertEqual(F.pam_side_from_spacer(SPACER[:10] + "-" + SPACER[10:], 3),
                         (F.PAM_3PRIME, 3))

    def run_folder(self, spacer, motif="NRG"):
        """A results folder whose PAM file is NOT there -- the run records the motif in
        .Params.txt and names a --pam the folder was copied without."""
        job = self.path("job1")
        os.makedirs(job)
        with open(os.path.join(job, ".command_line.txt"), "w") as fh:
            fh.write("input_command\tcrisprme.py complete-search --pam "
                     "PAMs/20bp-NRG-SpCas9.txt --output job1\n")
        with open(os.path.join(job, ".Params.txt"), "w") as fh:
            fh.write("Genome_selected\thg38\nPam\t%s\n" % motif)
        return job, write_targets(self.path("t.tsv"), ROW, spacer=spacer)

    def test_the_run_gives_the_motif_and_the_table_gives_the_side(self):
        job, table = self.run_folder(self.LEADING_N)
        said = []
        F.run(results=job, targets=table, reference=write_fasta(self.path("ref")),
              reference_only=True, out=self.path("o"), log=said.append)
        got = entries(read_tsv(self.path("o.tsv")))[0]
        fasta = read_fasta(self.path("o.fa"))
        self.assertEqual(got["pam_side"], "3prime")
        # the 3 beyond the protospacer at the genomic left, the 6 beyond the PAM at the
        # right: the whole point of reading the side rather than assuming it
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 3:TARGET_AT + 29])
        self.assertEqual(got["target_offset"], "3")
        self.assertTrue(any("its length from the run's Params.txt" in m for m in said),
                        said)

    def test_a_five_prime_run_read_the_same_way_lands_on_the_other_sides(self):
        job, table = self.run_folder("NNNN" + TARGET[4:] + "N", motif="TTTV")
        F.run(results=job, targets=table, reference=write_fasta(self.path("ref")),
              reference_only=True, out=self.path("o"))
        got = entries(read_tsv(self.path("o.tsv")))[0]
        self.assertEqual(got["pam_side"], "5prime")
        self.assertEqual(read_fasta(self.path("o.fa"))[got["fasta_id"]],
                         REF[TARGET_AT - 6:TARGET_AT + 26])

    def test_a_side_that_cannot_be_read_is_refused_not_assumed(self):
        """The same folder with a guide whose two ends cannot be told apart: the run
        says the motif, nothing says the side, and an asymmetric window is refused."""
        job, table = self.run_folder("NNN" + TARGET[3:20] + "NNN")
        with self.assertRaises(F.FlankError) as e:
            F.run(results=job, targets=table, reference=write_fasta(self.path("ref")),
                  reference_only=True, out=self.path("o"))
        self.assertIn("--pam-side", str(e.exception))


class BulgeFootprint(Case):
    """A bulge makes the target's genomic footprint differ from the nominal spacer+PAM
    length -- longer by one per DNA bulge, shorter by one per RNA bulge -- which is
    exactly what a fixed 32-mer assumes away. Both policies, and both said per row."""

    DNA_SPACER = SPACER[:5] + "-" + SPACER[5:]      # the GUIDE carries the gap

    def dna_row(self):
        return [("chr1", TARGET_AT, "+", REF[TARGET_AT:TARGET_AT + 24])]

    def rna_row(self):
        # 23 alignment columns, one of them a gap in the PROTOSPACER: 22 genomic bases
        bases = REF[TARGET_AT:TARGET_AT + 22]
        return [("chr1", TARGET_AT, "+", bases[:11] + "-" + bases[11:])]

    def test_a_dna_bulge_lengthens_the_window_and_says_so(self):
        _s, table, fasta = self.go(self.dna_row(), window=None,
                                   spacer=self.DNA_SPACER)
        got = entries(table)[0]
        self.assertEqual(got["bulge_type"], "DNA:1")
        self.assertEqual(got["target_span_ref"], "24")
        self.assertEqual(got["nominal_target_len"], "23")
        self.assertEqual(got["flank_len"], "33")          # 3 + 24 + 6
        self.assertEqual((got["pad_upstream"], got["pad_downstream"]), ("3", "6"))
        self.assertIn("target-span-not-nominal:24/23", got["notes"])
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 3:TARGET_AT + 30])

    def test_an_rna_bulge_shortens_it(self):
        _s, table, fasta = self.go(self.rna_row(), window=None)
        got = entries(table)[0]
        self.assertEqual(got["bulge_type"], "RNA:1")
        self.assertEqual(got["target_span_ref"], "22")
        self.assertEqual(got["flank_len"], "31")          # 3 + 22 + 6
        self.assertEqual(got["target_found"], "yes")
        self.assertIn("target-span-not-nominal:22/23", got["notes"])

    def test_fixed_total_holds_the_length_and_says_which_padding_moved(self):
        for rows, spacer, pads in ((self.dna_row(), self.DNA_SPACER, ("3", "5")),
                                   (self.rna_row(), SPACER, ("3", "7"))):
            _s, table, fasta = self.go(rows, window=None, spacer=spacer,
                                       length_policy="fixed-total")
            got = entries(table)[0]
            self.assertEqual(got["flank_len"], "32", got["notes"])
            self.assertEqual((got["pad_upstream"], got["pad_downstream"]), pads)
            self.assertEqual(got["length_policy"], "fixed-total")
            self.assertIn("padding-moved-to-hold-the-total", got["notes"])
            self.assertEqual(len(fasta[got["fasta_id"]]), 32)

    def test_a_padding_is_never_pushed_below_zero(self):
        """A bulge longer than the padding cannot make a negative one: the rest comes
        off the other end, and the row still says what it emitted."""
        long_spacer = SPACER[:5] + "---" + SPACER[5:]
        rows = [("chr1", TARGET_AT, "+", REF[TARGET_AT:TARGET_AT + 26])]
        _s, table, fasta = self.go(rows, window=None, spacer=long_spacer,
                                   upstream=1, downstream=1, length_policy="fixed-total")
        got = entries(table)[0]
        self.assertEqual(got["bulge_type"], "DNA:3")
        self.assertEqual((got["pad_upstream"], got["pad_downstream"]), ("0", "0"))
        self.assertEqual(got["flank_len"], "26")
        self.assertEqual(len(fasta[got["fasta_id"]]), 26)


# --------------------------------------------------------------------------- #
# the column mode: the input table back out, with the sequences on the row
# --------------------------------------------------------------------------- #
class ColumnMode(Case):
    def annotate(self, rows, records=(), curated=True, **kw):
        genome = write_fasta(self.path("ref"))
        vcf = write_vcf_gzip(self.path("v.vcf.gz"), records) if records else None
        table = write_targets(self.path("t.tsv"), rows, curated=curated,
                              spacer=kw.pop("spacer", SPACER))
        out, ann = self.path("out"), self.path("ann.tsv")
        summary = F.run(targets=table, reference=genome,
                        vcfs=[vcf] if vcf else [],
                        samples=[write_samples(self.path("s.txt"))] if vcf else [],
                        out=out, annotate_table=ann, **kw)
        return (summary, read_tsv(ann), read_tsv(out + ".tsv"),
                read_fasta(out + ".fa"))

    def test_the_curated_download_comes_back_with_the_sequences_on_the_row(self):
        rec = snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                  {"S1": "1|0"})
        _s, ann, _tsv, _fa = self.annotate(ROW, [rec])
        self.assertEqual(len(ann), 1)
        row = ann[0]
        # every original column is still there, in its own order
        self.assertEqual(row["rank"], "1")
        self.assertEqual(row["Aligned_protospacer+PAM"], TARGET)
        self.assertEqual(row["flank_ref"], REF[TARGET_AT - 3:TARGET_AT + 29])
        self.assertEqual(row["flank_ref_len"], "32")
        self.assertEqual(row["flank_ref_offset"], "3")
        self.assertEqual(row["flank_alt_count"], "1")
        self.assertEqual(len(row["flank_alt"]), 32)
        self.assertNotEqual(row["flank_alt"], row["flank_ref"])
        self.assertEqual(row["flank_window"], "up3+down6")
        self.assertEqual(row["flank_row_id"], "1")

    def test_an_integrated_table_too(self):
        _s, ann, _tsv, _fa = self.annotate(ROW, curated=False)
        row = ann[0]
        self.assertEqual(row["Spacer+PAM"], SPACER)
        self.assertEqual(row["flank_ref"], REF[TARGET_AT - 3:TARGET_AT + 29])
        self.assertEqual(row["flank_alt"], "")
        self.assertEqual(row["flank_alt_count"], "0")

    def test_several_alts_are_counted_and_kept_not_dropped(self):
        """A row with more than one distinct ALT sequence: the columns carry the
        primary one and the COUNT, and the sequences are in the companion files,
        reachable by flank_row_id -- all of them here, two being under the default cap
        of ten; nothing is crammed into a cell."""
        recs = [snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                    {"S1": "1|0"}),
                snv(TARGET_AT + 25, "A" if REF[TARGET_AT + 25] != "A" else "C",
                    {"S2": "1|0"})]
        summary, ann, tsv, fasta = self.annotate(ROW, recs)
        row = ann[0]
        self.assertEqual(row["flank_alt_count"], "2")
        self.assertEqual(summary["multi_alt_rows"], 1)
        companion = [r for r in tsv if r["haplotype"].startswith("hap")]
        self.assertEqual(len(companion), 2)
        self.assertEqual({r["row_id"] for r in companion}, {row["flank_row_id"]})
        # the column carries hap1 itself -- the same bytes the FASTA holds -- and the
        # other one is in the companion files, not lost
        self.assertEqual(row["flank_alt"], fasta[companion[0]["fasta_id"]])
        self.assertEqual({fasta[r["fasta_id"]] for r in companion} - {row["flank_alt"]},
                         {fasta[companion[1]["fasta_id"]]})
        # the extras are NOT in the cell by default
        self.assertNotIn("|", row["flank_alt"])
        self.assertNotIn("flank_alt_all", row)

    def test_alt_list_puts_them_all_in_one_cell_when_asked(self):
        recs = [snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                    {"S1": "1|0"}),
                snv(TARGET_AT + 25, "A" if REF[TARGET_AT + 25] != "A" else "C",
                    {"S2": "1|0"})]
        _s, ann, _tsv, _fa = self.annotate(ROW, recs, alt_list=True)
        every = ann[0]["flank_alt_all"].split("|")
        self.assertEqual(len(every), 2)
        self.assertEqual(len(set(every)), 2)
        self.assertIn(ann[0]["flank_alt"], every)

    THREE = [snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                 {"S1": "1|0", "S2": "1|0", "S3": "1|0"}),      # 3 carriers: hap1
             snv(TARGET_AT + 25, "A" if REF[TARGET_AT + 25] != "A" else "C",
                 {"S4": "1|1"}),                                # 2 carriers: hap2
             snv(TARGET_AT + 27, "A" if REF[TARGET_AT + 27] != "A" else "C",
                 {"S5": "1|0"})]                                # 1 carrier:  hap3

    def test_the_alt_list_is_in_haplotype_order_hap1_first(self):
        """The cell is positional: element k IS hapK of the companion TSV, and the
        first element IS flank_alt. A reader pairing the cell with the hap names has
        nothing else to go on."""
        _s, ann, tsv, fasta = self.annotate(ROW, self.THREE, alt_list=True)
        haps = sorted((r for r in tsv if r["haplotype"].startswith("hap")),
                      key=lambda r: int(r["haplotype"][3:]))
        self.assertEqual([r["haplotype"] for r in haps], ["hap1", "hap2", "hap3"])
        every = ann[0]["flank_alt_all"].split("|")
        self.assertEqual(every, [fasta[r["fasta_id"]] for r in haps])
        self.assertEqual(every[0], ann[0]["flank_alt"])
        self.assertEqual(len(set(every)), 3)       # three sequences, three elements

    def test_the_alt_carriers_column_counts_the_primary_alts_carriers(self):
        """flank_alt_carriers is how many haplotypes carry THE SEQUENCE IN flank_alt.
        The reference entry has carriers of its own -- the samples whose other copy
        reads as the reference here -- and they are a different number of different
        people."""
        _s, ann, tsv, _fa = self.annotate(ROW, self.THREE)
        ref = [r for r in tsv if r["haplotype"] == "ref"][0]
        hap1 = [r for r in tsv if r["haplotype"] == "hap1"][0]
        self.assertEqual(hap1["n_carriers"], "3")
        self.assertEqual(ref["n_carriers"], "4")
        self.assertEqual(ann[0]["flank_alt_carriers"], hap1["n_carriers"])
        self.assertNotEqual(ann[0]["flank_alt_carriers"], ref["n_carriers"])

    def test_the_cap_never_makes_the_count_or_the_notes_lie(self):
        """--max-haplotypes keeps the file readable; it must not make the table say a
        row has fewer ALT sequences than it has. The count is the true total, the
        delimited list holds what was written, and the notes say what was left out."""
        recs = [snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                    {"S1": "1|0"}),
                snv(TARGET_AT + 25, "A" if REF[TARGET_AT + 25] != "A" else "C",
                    {"S2": "1|0"})]
        _s, ann, _tsv, _fa = self.annotate(ROW, recs, alt_list=True, max_haplotypes=1)
        row = ann[0]
        self.assertEqual(row["flank_alt_count"], "2")                  # the true total
        self.assertEqual(len(row["flank_alt_all"].split("|")), 1)      # what was written
        self.assertIn("haplotypes-not-written:1", row["flank_notes"])

    def test_the_original_columns_come_back_byte_for_byte(self):
        """The user's table must come back as THEIR table: the original text is copied,
        not round-tripped through a CSV reader that would eat a quote character, and
        the line count is the same whatever the filter did."""
        genome = write_fasta(self.path("ref"))
        table = self.path("t.tsv")
        head = ('rank\tChromosome\tPosition\tStrand\tAligned_spacer+PAM\t'
                'Aligned_protospacer+PAM\tNote')
        body = ['1\tchr1\t%d\t+\t%s\t%s\t"quoted, with a comma"' % (TARGET_AT, SPACER,
                                                                    TARGET),
                '2\tchr1\t%d\t+\t%s\t%s\t' % (TARGET_AT, SPACER, TARGET)]
        with open(table, "w") as fh:
            fh.write(head + "\n" + "\n".join(body) + "\n")
        F.run(targets=table, reference=genome, out=self.path("o"),
              annotate_table=self.path("a.tsv"))
        with open(self.path("a.tsv")) as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith(head + "\t"))
        for original, written in zip(body, lines[1:]):
            self.assertTrue(written.startswith(original + "\t"), written)
        self.assertIn('"quoted, with a comma"', lines[1])

    def test_a_row_past_max_rows_says_so_rather_than_reading_as_empty(self):
        """--max-rows stops the READING (a genome-wide table has millions of lines), so
        the lines past it have no flank -- and must say which of the two reasons it is,
        or an empty cell reads as "this site has no ALT sequence"."""
        genome = write_fasta(self.path("ref"))
        table = write_targets(self.path("t.tsv"),
                              [("chr1", TARGET_AT, "+", TARGET)] * 3)
        F.run(targets=table, reference=genome, out=self.path("o"), max_rows=1,
              annotate_table=self.path("a.tsv"))
        got = read_tsv(self.path("a.tsv"))
        self.assertEqual(len(got), 3)
        self.assertEqual(got[0]["flank_ref_len"], "32")
        for row in got[1:]:
            self.assertEqual(row["flank_ref"], "")
            self.assertEqual(row["flank_notes"], "not-read(--max-rows)")

    def test_a_column_it_would_overwrite_is_refused(self):
        genome = write_fasta(self.path("ref"))
        table = self.path("t.tsv")
        with open(table, "w") as fh:
            fh.write("rank\tChromosome\tPosition\tStrand\tAligned_spacer+PAM\t"
                     "Aligned_protospacer+PAM\tflank_ref\n")
            fh.write("1\tchr1\t%d\t+\t%s\t%s\tmine\n" % (TARGET_AT, SPACER, TARGET))
        with self.assertRaises(F.FlankError) as e:
            F.run(targets=table, reference=genome, out=self.path("o"),
                  annotate_table=self.path("a.tsv"))
        self.assertIn("flank_ref", str(e.exception))
        self.assertIn("overwrite", str(e.exception))

    def test_writing_over_the_input_table_is_refused(self):
        genome = write_fasta(self.path("ref"))
        table = write_targets(self.path("t.tsv"), ROW)
        with self.assertRaises(F.FlankError) as e:
            F.run(targets=table, reference=genome, out=self.path("o"),
                  annotate_table=table)
        self.assertIn("never modified in place", str(e.exception))


class RowFilters(Case):
    """--filter: which rows get flanks, as an expression over the table's own columns."""

    EXPR = "CRISPR_Bulge >= CB_SENSITIVE and CRISPR_Bulge_ref >= CB_SENSITIVE"

    def scored(self, alt, ref):
        rows = [("chr1", TARGET_AT, "+", TARGET)] * len(alt)
        return rows, [("CRISPR-Bulge", alt), ("CRISPR-Bulge_ref", ref)]

    def test_both_scores_above_the_sensitive_threshold(self):
        """The worked example: the pair only where the ALT and the REF alignment are
        both at or above the CRISPR-Bulge sensitive operating point."""
        import cb_operating_point as cbop
        thr = cbop.sensitive_threshold()
        rows, extra = self.scored([thr * 10, thr / 10, thr * 10, thr * 10],
                                  [thr * 10, thr * 10, thr / 10, thr])
        summary, table, _f = self.go(rows, window=None, extra=extra,
                                     row_filter=self.EXPR)
        self.assertEqual({int(r["row_id"]) for r in table}, {1, 4})
        self.assertEqual(summary["filter_counts"]["kept"], 2)
        self.assertEqual(summary["filter_counts"]["read"], 4)

    def test_the_threshold_is_the_operating_points_not_a_literal(self):
        """CB_SENSITIVE follows cb_operating_point, so a filter written today still
        means 'the sensitive threshold' after a re-fit or a different recall target."""
        import cb_operating_point as cbop
        rows, extra = self.scored([2.0e-3], [2.0e-3])
        at_95 = cbop.CB_RECALL_POINTS[95][0]
        at_90 = cbop.CB_RECALL_POINTS[90][0]
        self.assertLess(at_95, 2.0e-3)
        self.assertGreater(at_90, 2.0e-3)
        _s, kept, _f = self.go(rows, window=None, extra=extra, row_filter=self.EXPR)
        self.assertEqual(len(entries(kept)), 1)
        with mock.patch.dict(os.environ, {"CRISPRME_CB_RECALL_TARGET": "90"}):
            with self.assertRaises(F.FlankError) as e:
                self.go(rows, window=None, extra=extra, row_filter=self.EXPR)
        self.assertIn("kept none", str(e.exception))

    def test_a_missing_value_does_not_pass_and_is_counted(self):
        """A cell the table does not fill is not zero: the row is not kept, and how
        many such rows there were is said rather than left to be noticed."""
        rows, extra = self.scored([0.5, 0.5], ["NA", 0.5])
        summary, table, _f = self.go(rows, window=None, extra=extra,
                                     row_filter=self.EXPR)
        self.assertEqual({int(r["row_id"]) for r in table}, {2})
        self.assertEqual(summary["filter_counts"]["filter_missing"], 1)

    def test_a_negated_test_does_not_let_the_missing_rows_through(self):
        """The rule is one rule: a row whose value the expression READ is missing is
        not kept, whichever way round the test is written. Otherwise 'not above the
        threshold' would quietly select every row nobody measured -- the opposite of
        what a validation shortlist wants. has() asks about presence without reading
        the value, and short-circuits, so the escape hatch still works."""
        rows, extra = self.scored([0.5, "NA", 1e-9], [0.5, 0.5, 0.5])
        _s, table, _f = self.go(rows, window=None, extra=extra,
                                row_filter="not (CRISPR_Bulge >= CB_SENSITIVE)")
        self.assertEqual({int(r["row_id"]) for r in table}, {3})
        _s2, with_hatch, _f2 = self.go(
            rows, window=None, extra=extra,
            row_filter='not has("CRISPR-Bulge") or CRISPR_Bulge < CB_SENSITIVE')
        self.assertEqual({int(r["row_id"]) for r in with_hatch}, {2, 3})

    def test_a_name_the_table_does_not_have_is_refused_by_name(self):
        rows, extra = self.scored([0.5], [0.5])
        with self.assertRaises(F.FlankError) as e:
            self.go(rows, window=None, extra=extra, row_filter="CRISPR_Bulgé > 0")
        self.assertIn("has no column", str(e.exception))
        self.assertIn("CRISPR-Bulge", str(e.exception))     # the columns are listed

    def test_the_exact_name_is_always_reachable(self):
        rows, extra = self.scored([0.5, 1e-9], [0.5, 0.5])
        _s, table, _f = self.go(rows, window=None, extra=extra,
                                row_filter='col("CRISPR-Bulge") >= CB_SENSITIVE')
        self.assertEqual({int(r["row_id"]) for r in table}, {1})

    def test_an_expression_that_is_not_a_row_test_is_refused(self):
        rows, extra = self.scored([0.5], [0.5])
        for bad, why in ((" __import__('os').system('true')", "may only call"),
                         ("CFD.real > 0", "may not use"),
                         ("[c for c in (1,)]", "may not use"),
                         ("abs > 1", "is a function"),
                         ('col("nope") > 1', 'has no column "nope"'),
                         ("CFD >", "is not an expression")):
            with self.assertRaises(F.FlankError) as e:
                self.go(rows, window=None, extra=extra, row_filter=bad)
            self.assertIn(why, str(e.exception), bad)

    def test_a_filtered_row_keeps_its_line_in_the_table_unless_dropped(self):
        rows, extra = self.scored([0.5, 1e-9], [0.5, 0.5])
        genome = write_fasta(self.path("ref"))
        table = write_targets(self.path("t.tsv"), rows, extra=extra)
        F.run(targets=table, reference=genome, out=self.path("o"),
              row_filter=self.EXPR, annotate_table=self.path("a.tsv"))
        kept = read_tsv(self.path("a.tsv"))
        self.assertEqual(len(kept), 2)
        self.assertEqual(kept[1]["flank_notes"], "filtered-out")
        self.assertEqual(kept[1]["flank_ref"], "")
        F.run(targets=table, reference=genome, out=self.path("o2"),
              row_filter=self.EXPR, annotate_table=self.path("a2.tsv"),
              filter_drop=True)
        self.assertEqual(len(read_tsv(self.path("a2.tsv"))), 1)

    def test_a_dash_is_the_value_of_a_column_that_has_one(self):
        """'-' is CRISPRme's "not measured" in a score column and the MINUS STRAND in
        Strand (and an intergenic Gene_region, and a PAM_creation that did not happen).
        A filter over the table's own columns has to serve both: compared with a string
        the cell is its own text, so these select rows rather than report absence."""
        rows = [("chr1", TARGET_AT, "+" if i % 2 else "-", TARGET if i % 2
                 else _rc(TARGET)) for i in range(4)]
        extra = [("Gene_region", ["-", "intron", "-", "exon"])]
        summary, table, _f = self.go(rows, window=None, extra=extra,
                                     row_filter='Strand == "-"')
        self.assertEqual({int(r["row_id"]) for r in table}, {1, 3})
        self.assertEqual(summary["filter_counts"]["filter_missing"], 0)
        _s2, both, _f2 = self.go(rows, window=None, extra=extra,
                                 row_filter='Strand in ("+", "-")')
        self.assertEqual({int(r["row_id"]) for r in both}, {1, 2, 3, 4})
        _s3, gene, _f3 = self.go(rows, window=None, extra=extra,
                                 row_filter='Gene_region == "-"')
        self.assertEqual({int(r["row_id"]) for r in gene}, {1, 3})

    def test_a_dash_where_a_number_belongs_is_still_missing(self):
        """The other half of the same rule, and the one the worked example depends on:
        '-' in a SCORE column is not a value to compare with a threshold. The row is not
        kept and it is counted -- on one real 1,000-row shortlist, the 14 rows whose two
        CRISPR-Bulge columns are both '-'."""
        rows, extra = self.scored([0.5, "-", 0.5], [0.5, 0.5, "."])
        summary, table, _f = self.go(rows, window=None, extra=extra,
                                     row_filter=self.EXPR)
        self.assertEqual({int(r["row_id"]) for r in table}, {1})
        self.assertEqual(summary["filter_counts"]["filter_missing"], 2)
        # and a disjunction does not drop a row that passes one half because the OTHER
        # half's cell holds a '-' that is a value
        rows2 = [("chr1", TARGET_AT, "+", TARGET)] * 3
        extra2 = [("Score", [0.9, 0.1, 0.1]), ("Gene_region", ["intron", "-", "exon"])]
        _s2, kept, _f2 = self.go(rows2, window=None, extra=extra2,
                                 row_filter='Score > 0.5 or Gene_region == "-"')
        self.assertEqual({int(r["row_id"]) for r in kept}, {1, 2})


class JoiningBackToTheInputTable(Case):
    """The join key of every output is the input's own data line number. The reader and
    the --annotate-table writer have to count lines THE SAME WAY, or a user reads one
    site's sequence off another site's row with nothing to see."""

    def table_with(self, cells):
        """A 4-row curated table, one site per row, 30 bases apart, with ``cells`` as
        the last column -- written as given, quote characters and all."""
        path = self.path("t.tsv")
        with open(path, "w") as fh:
            fh.write("\t".join(["rank", "Chromosome", "Position", "Strand",
                                "Aligned_spacer+PAM", "Aligned_protospacer+PAM",
                                "Note"]) + "\n")
            for i, cell in enumerate(cells):
                fh.write("\t".join([str(i + 1), "chr1", str(100 + 30 * i), "+", SPACER,
                                    TARGET, cell]) + "\n")
        return path

    def at(self, position, window=20):
        return REF[position - window:position + len(TARGET) + window]

    def test_a_quote_in_a_cell_does_not_move_a_sequence_onto_another_row(self):
        """A quote opened in one cell and closed in another: csv's default dialect
        joins those lines into ONE row, while the writer counts them as two, and every
        row after it then carries the next site's sequence. A TSV has no quoting: the
        table is read one row per line."""
        table = self.table_with(['plain', '"gene X', 'end here"', 'plain'])
        F.run(targets=table, reference=write_fasta(self.path("ref")), window=20,
              out=self.path("o"), reference_only=True,
              annotate_table=self.path("a.tsv"))
        got = read_tsv(self.path("a.tsv"))
        self.assertEqual(len(got), 4)
        for i, row in enumerate(got):
            self.assertEqual(row["flank_row_id"], str(i + 1))
            self.assertEqual(row["flank_ref"], self.at(100 + 30 * i), row["Position"])
        # the quote characters themselves come back untouched
        self.assertEqual(got[1]["Note"], '"gene X')

    def test_a_cell_holding_a_newline_is_a_line_not_a_second_row_s_sequence(self):
        """What Excel's tab-delimited export writes for a cell with a line break in it.
        Those bytes are not a row of the table: the line they leave behind comes back
        saying so, and every real row keeps its own sequence."""
        table = self.table_with(['"first line', 'plain', 'plain', 'plain'])
        with open(table) as fh:
            raw = fh.read()
        with open(table, "w") as fh:          # close the quote on the NEXT line
            fh.write(raw.replace('"first line', '"first\nsecond"', 1))
        summary = F.run(targets=table, reference=write_fasta(self.path("ref")),
                        window=20, out=self.path("o"), reference_only=True,
                        annotate_table=self.path("a.tsv"))
        tsv = read_tsv(self.path("o.tsv"))
        self.assertEqual(summary["rows"], 5)             # the stray line is a line
        for i in range(4):
            seq = [r for r in tsv if r["row_id"] == str(i + 1 + (i > 0))]
            self.assertEqual(seq[0]["flank_start"], str(100 + 30 * i - 20 + 1))
        stray = [r for r in tsv if r["notes"].endswith("no-chromosome")]
        self.assertEqual(len(stray), 1)
        self.assertEqual(stray[0]["haplotype"], "unresolved")

    def test_a_reader_and_a_writer_that_disagree_write_no_table(self):
        """The tripwire under the rule above: if the two ever count differently again,
        the column mode refuses rather than hand over a table whose every number looks
        right. Nothing is left behind to be read by mistake."""
        table = self.table_with(["plain"] * 4)
        header = ["rank", "Chromosome", "Position", "Strand", "Aligned_spacer+PAM",
                  "Aligned_protospacer+PAM", "Note"]
        out = self.path("a.tsv")
        with self.assertRaises(F.FlankError) as e:
            F.write_annotated_table(table, out, {}, header,
                                    counts={"read": 3, "truncated": False,
                                            "last_line": 3})
        self.assertIn("4 data line(s)", str(e.exception))
        self.assertIn("read as 3 row(s)", str(e.exception))
        self.assertFalse(os.path.exists(out))
        # the same call with the true count writes the table
        F.write_annotated_table(table, out, {}, header,
                                counts={"read": 4, "truncated": False, "last_line": 4})
        self.assertEqual(len(read_tsv(out)), 4)


class AHaplotypeThatDeletesTheWindow(Case):
    """A deletion that covers the whole window: that genome holds NO bases between the
    flank's two coordinates. Real at a short window (one entry of 4,017 on a 1,000-row
    shortlist at 3/6; none of 10,354 of the same shortlist at --window 200)."""

    def go_deleted(self, **kw):
        # 40 reference bases removed, over the whole 32-base window of the row; plus an
        # ordinary SNV carrier, so the row has a second ALT and an order to get wrong
        recs = [deletion(TARGET_AT - 5, 40, {"S1": "1|0"}),
                snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                    {"S2": "1|0", "S3": "1|1"})]
        return self.go(ROW, recs, window=None, max_haplotypes=0, **kw)

    def test_the_row_says_the_window_is_deleted_rather_than_only_a_zero_length(self):
        _s, table, _f = self.go_deleted()
        gone = [r for r in table if r["flank_len"] == "0"
                and r["haplotype"] != "unresolved"]
        self.assertEqual(len(gone), 1)
        self.assertEqual(gone[0]["n_carriers"], "1")
        self.assertIn(F.NOTE_WINDOW_DELETED, gone[0]["notes"])

    def test_no_fasta_record_is_written_without_bases(self):
        """A header followed by a blank line is not a sequence: it is how every record
        after it gets read as the wrong haplotype."""
        said = []
        summary, table, fasta = self.go_deleted(log=said.append)
        self.assertEqual(summary["window_deleted"], 1)
        # and the log says why the FASTA holds one record fewer than the TSV has entries
        self.assertTrue(any(F.NOTE_WINDOW_DELETED in m for m in said), said)
        self.assertEqual(len(fasta), summary["sequences"])
        self.assertTrue(all(fasta.values()), fasta)
        gone = [r for r in table if r["flank_len"] == "0"
                and r["haplotype"] != "unresolved"][0]
        self.assertEqual(gone["fasta_id"], "")      # nothing points at a missing record
        self.assertNotIn("crisprme_row1_%s" % gone["haplotype"], fasta)
        # the entries that DO have bases are all there
        self.assertEqual(sorted(fasta), ["crisprme_row1_hap1", "crisprme_row1_ref"])

    def test_the_alt_list_keeps_one_element_per_haplotype(self):
        """Element k of the cell is hapK. An empty element is dropped by the obvious
        split, and then every sequence after it is read as the wrong haplotype."""
        genome = write_fasta(self.path("ref"))
        vcf = write_vcf_gzip(self.path("v.vcf.gz"),
                             [deletion(TARGET_AT - 5, 40, {"S1": "1|0"}),
                              snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A"
                                  else "C", {"S2": "1|0", "S3": "1|1"})])
        F.run(targets=write_targets(self.path("t.tsv"), ROW), reference=genome,
              vcfs=[vcf], samples=[write_samples(self.path("s.txt"))],
              out=self.path("o"), annotate_table=self.path("a.tsv"), alt_list=True,
              max_haplotypes=0)
        row = read_tsv(self.path("a.tsv"))[0]
        haps = sorted((r for r in read_tsv(self.path("o.tsv"))
                       if r["haplotype"].startswith("hap")),
                      key=lambda r: int(r["haplotype"][3:]))
        every = row["flank_alt_all"].split("|")
        self.assertEqual(len(every), len(haps))
        self.assertEqual(len([p for p in every if p]), len(haps))   # none is empty
        self.assertEqual(every[0], row["flank_alt"])
        self.assertEqual(every[1], F.SEQ_DELETED)
        self.assertEqual(row["flank_notes"].split(";")[-1],
                         "%s:%s" % (F.NOTE_WINDOW_DELETED, haps[1]["haplotype"]))

    def test_the_primary_alt_cell_is_never_blank_for_a_deleted_window(self):
        """An empty flank_alt means "this row has no ALT sequence". If the primary ALT
        were the deleted one, the cell would say that and be wrong."""
        cells = F.annotation_rows([
            {"row_id": 1, "haplotype": "ref", "_seq": "ACGT", "flank_len": 4,
             "target_offset": 0, "target_len": 4, "n_haplotypes": 2, "n_carriers": 0,
             "notes": "", "window_spec": "up3+down6", "orientation": "guide",
             "emitted_strand": "+"},
            {"row_id": 1, "haplotype": "hap1", "_seq": "", "flank_len": 0,
             "target_offset": 0, "target_len": 0, "n_haplotypes": 2, "n_carriers": 7,
             "notes": F.NOTE_WINDOW_DELETED, "window_spec": "up3+down6",
             "orientation": "guide", "emitted_strand": "+"}])[1]
        self.assertEqual(cells["flank_alt"], F.SEQ_DELETED)
        self.assertEqual(cells["flank_alt_len"], 0)
        self.assertIn("%s:hap1" % F.NOTE_WINDOW_DELETED, cells["flank_notes"])


class Presets(Case):
    def test_a_preset_names_a_window_and_a_flag_still_wins(self):
        _s, table, fasta = self.go(ROW, window=None, preset="amplicon")
        got = entries(table)[0]
        self.assertEqual(got["window_spec"], "window200")
        self.assertEqual(fasta[got["fasta_id"]],
                         REF[TARGET_AT - 200:TARGET_AT + 223])
        # the flag beats the preset
        _s2, table2, _f2 = self.go(ROW, window=40, preset="amplicon")
        self.assertEqual(entries(table2)[0]["window_spec"], "window40")

    def test_a_preset_file_adds_one_without_a_code_change(self):
        """The requirement that makes this general: another assay is a section in a
        file, not an edit to the module."""
        path = self.path("my.conf")
        with open(path, "w") as fh:
            fh.write("[mine]\nupstream = 5\ndownstream = 1\n"
                     "description = my own assay\n")
        _s, table, fasta = self.go(ROW, window=None, preset="mine",
                                   preset_files=[path])
        got = entries(table)[0]
        self.assertEqual(got["window_spec"], "up5+down1")
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 5:TARGET_AT + 24])
        self.assertIn("mine", F.preset_help(F.load_presets([path])))

    def test_a_preset_may_carry_the_filter_and_the_cap_too(self):
        """A preset is a name for a SET of flags, not for a window alone -- and
        --max-haplotypes 10 (the default's own value) still beats a preset that says
        something else, because "not given" is told from "given the default"."""
        import cb_operating_point as cbop
        thr = cbop.sensitive_threshold()
        path = self.path("my.conf")
        with open(path, "w") as fh:
            fh.write("[assay]\nupstream = 3\ndownstream = 6\nmax_haplotypes = 1\n"
                     "filter = CRISPR_Bulge >= CB_SENSITIVE\n")
        rows = [("chr1", TARGET_AT, "+", TARGET)] * 2
        extra = [("CRISPR-Bulge", [thr * 10, thr / 10])]
        recs = [snv(TARGET_AT - 2, "A" if REF[TARGET_AT - 2] != "A" else "C",
                    {"S1": "1|0"}),
                snv(TARGET_AT + 25, "A" if REF[TARGET_AT + 25] != "A" else "C",
                    {"S2": "1|0"})]
        _s, table, _f = self.go(rows, recs, window=None, extra=extra,
                                preset="assay", preset_files=[path])
        self.assertEqual({int(r["row_id"]) for r in table}, {1})      # the filter
        self.assertEqual(len([r for r in table
                              if r["haplotype"].startswith("hap")]), 1)   # the cap
        _s2, more, _f2 = self.go(rows, recs, window=None, extra=extra,
                                 preset="assay", preset_files=[path],
                                 max_haplotypes=10)
        self.assertEqual(len([r for r in more
                              if r["haplotype"].startswith("hap")]), 2)

    def test_the_presets_are_listed_by_the_help(self):
        text = F.help_text()
        for name in F.load_presets():
            self.assertIn(name, text)
        self.assertIn("PRESETS (--preset <name>)", text)

    def test_a_flag_wins_over_the_presets_other_form_of_the_window(self):
        """"A flag on the command line always wins over the preset" is one rule, and it
        holds across the two FORMS of a window: --window W IS --upstream W --downstream
        W, so a preset's window is its two ends and a flag replaces the end it names.
        Refusing the pair would make the preset the documentation tells a user to write
        (upstream/downstream, as in flank_presets.conf) unusable with --window."""
        path = self.path("my.conf")
        with open(path, "w") as fh:
            fh.write("[mine]\nupstream = 3\ndownstream = 6\n")
        # the anchored preset + the symmetric flag
        _s, table, fasta = self.go(ROW, window=200, preset="mine", preset_files=[path])
        got = entries(table)[0]
        self.assertEqual(got["window_spec"], "window200")
        self.assertEqual(fasta[got["fasta_id"]], REF[TARGET_AT - 200:TARGET_AT + 223])
        # the symmetric preset + one anchored flag: the flag replaces that end only
        _s2, one, f2 = self.go(ROW, window=None, upstream=5, preset="amplicon")
        self.assertEqual(entries(one)[0]["window_spec"], "up5+down200")
        self.assertEqual(f2[entries(one)[0]["fasta_id"]],
                         REF[TARGET_AT - 5:TARGET_AT + 223])
        # ...and the shipped anchored preset with a symmetric flag, the other way round
        _s3, three, _f3 = self.go(ROW, window=50, preset="oligo")
        self.assertEqual(entries(three)[0]["window_spec"], "window50")

    def test_two_forms_from_one_source_are_still_refused(self):
        """What the refusal is actually for: one source saying two different things.
        On the command line, or inside one preset -- not a flag overriding a preset."""
        with self.assertRaises(F.FlankError) as e:
            self.go(ROW, window=40, upstream=3)
        self.assertIn("one form or the other", str(e.exception))
        path = self.path("bad.conf")
        with open(path, "w") as fh:
            fh.write("[muddled]\nwindow = 10\nupstream = 3\n")
        with self.assertRaises(F.FlankError) as e:
            self.go(ROW, window=None, preset="muddled", preset_files=[path])
        self.assertIn("preset [muddled]", str(e.exception))
        self.assertIn("one form or the other", str(e.exception))

    def test_an_unknown_preset_and_an_unknown_key_are_refused(self):
        with self.assertRaises(F.FlankError) as e:
            self.go(ROW, window=None, preset="nope")
        self.assertIn("no preset", str(e.exception))
        self.assertIn("amplicon", str(e.exception))     # the ones there are, named
        path = self.path("bad.conf")
        with open(path, "w") as fh:
            fh.write("[mine]\nupstraem = 5\n")
        with self.assertRaises(F.FlankError) as e:
            F.load_presets([path])
        self.assertIn("upstraem", str(e.exception))
        self.assertIn("is not a preset key", str(e.exception))


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


@contextlib.contextmanager
def expect_break(test, why):
    """The block must not PASS: the rule's test either fails an assertion or the run it
    makes refuses outright. Which of the two depends only on whether the wrong answer is
    still a runnable one -- a filter that now matches nothing, a window with no readable
    PAM side -- and either way the suite is not silent about the mutant, which is what
    these tests are for."""
    try:
        yield
    except (AssertionError, F.FlankError):
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

    # ---- the anchored window, the strand, the PAM's side, the bulge, the columns,
    # the filter: one mutant each
    def test_an_always_symmetric_window_breaks_the_default_rule(self):
        """The window before this change: one number, the same at both ends. A 32-mer
        is 3 and 6, so a symmetric reading of them is a different sequence."""
        real = F.Geometry.edges

        def symmetric(self_, strand, span, nominal):
            left, right, up, down, adj = real(self_, strand, span, nominal)
            return self_.upstream, self_.upstream, up, down, adj
        with mock.patch.object(F.Geometry, "edges", symmetric):
            with expect_failure(self, "both paddings take the upstream number"):
                AnchoredWindows("test_the_default_is_a_32mer_on_a_mismatch_only_site"
                                ).debug()

    def test_ignoring_the_strand_breaks_the_minus_strand_rule(self):
        """The paddings mapped to the genomic sides as if every row were on the plus
        strand: on a minus row the 3 and the 6 then swap ends."""
        real = F.Geometry.edges

        def plus_only(self_, strand, span, nominal):
            return real(self_, "+", span, nominal)
        with mock.patch.object(F.Geometry, "edges", plus_only):
            with expect_failure(self, "a minus row is padded as if it were plus"):
                StrandAndPamSide(
                    "test_a_minus_strand_row_pads_in_the_guides_orientation").debug()

    def test_assuming_spcas9_breaks_the_five_prime_pam_rule(self):
        """The thing this must never do: read every run as a 3' PAM. On a Cas12a run
        that puts both paddings on the wrong genomic side."""
        with mock.patch.object(F, "pam_side_from_spacer",
                               lambda spacer, motif_len=None: (F.PAM_3PRIME, 3)):
            with expect_failure(self, "a 5' PAM is read as a 3' one"):
                StrandAndPamSide(
                    "test_a_five_prime_pam_puts_the_paddings_on_the_other_genomic_sides"
                ).debug()

    def test_a_footprint_taken_as_nominal_breaks_the_bulge_rule(self):
        """The assumption the 32-mer makes: protospacer+PAM is always the nominal
        number of genomic bases. It is not, wherever the alignment carries a bulge."""
        with mock.patch.object(F, "bulge_counts", lambda spacer, proto: (0, 0)):
            real = F.Geometry.edges

            def nominal_span(self_, strand, span, nominal):
                return real(self_, strand, nominal if nominal else span, nominal)
            with mock.patch.object(F.Geometry, "edges", nominal_span):
                with expect_failure(self, "the footprint is assumed to be nominal"):
                    BulgeFootprint(
                        "test_fixed_total_holds_the_length_and_says_which_padding_moved"
                    ).debug()

    def test_an_alt_column_filled_from_the_reference_breaks_the_column_rule(self):
        """The quiet disaster of the column mode: the ALT cell carrying the reference
        sequence, which is the very thing a user must not order oligos from."""
        real = F.annotation_rows

        def ref_in_alt(out_lines):
            rows = real(out_lines)
            for cells in rows.values():
                cells["flank_alt"] = cells["flank_ref"]
            return rows
        with mock.patch.object(F, "annotation_rows", ref_in_alt):
            with expect_failure(self, "the ALT column holds the reference"):
                ColumnMode(
                    "test_the_curated_download_comes_back_with_the_sequences_on_the_row"
                ).debug()

    def test_counting_only_what_was_written_breaks_the_alt_count_rule(self):
        """The count must be the row's TRUE number of distinct ALT sequences, not the
        number that fitted: a user sizes extra assay targets from it."""
        real = F.annotation_rows

        def written_only(out_lines):
            rows = real(out_lines)
            for cells in rows.values():
                cells["flank_alt_count"] = 1 if cells["flank_alt"] else 0
            return rows
        with mock.patch.object(F, "annotation_rows", written_only):
            with expect_failure(self, "the ALT count is what was written"):
                ColumnMode(
                    "test_several_alts_are_counted_and_kept_not_dropped").debug()

    def test_a_missing_value_read_as_zero_breaks_the_filter_rule(self):
        """An empty cell read as 0 passes a 'greater than' test it should fail and
        fails a 'less than' one it should pass: either way the user gets a row they
        did not ask for, silently."""
        real = F.RowFilter._cell

        def as_zero(self_, i, row, state):
            text = row[i].strip() if i < len(row) else ""
            if text in F._MISSING_CELLS:
                return 1.0            # a number, and nothing said about it
            return real(self_, i, row, state)
        with mock.patch.object(F.RowFilter, "_cell", as_zero):
            with expect_failure(self, "a missing value is read as a number"):
                RowFilters("test_a_missing_value_does_not_pass_and_is_counted").debug()

    def test_a_hardcoded_threshold_breaks_the_operating_point_rule(self):
        """CB_SENSITIVE written as the number it happens to be today: the filter then
        keeps saying 1.9e-4 after the operating point moves."""
        real = F.cb_sensitive_point

        def frozen():
            _thr, recall = real()
            return 1.9e-4, recall
        with mock.patch.object(F, "cb_sensitive_point", frozen):
            with expect_failure(self, "the threshold is a literal, not the point"):
                RowFilters(
                    "test_the_threshold_is_the_operating_points_not_a_literal").debug()

    def test_a_preset_that_beats_the_flag_breaks_the_precedence_rule(self):
        """A preset must supply a default, never override what the user typed."""
        real = F.resolve_window

        def preset_wins(upstream, downstream, window, preset, presets):
            if preset:
                return real(None, None, None, preset, presets)
            return real(upstream, downstream, window, preset, presets)
        with mock.patch.object(F, "resolve_window", preset_wins):
            with expect_failure(self, "the preset beats the flag"):
                Presets("test_a_preset_names_a_window_and_a_flag_still_wins").debug()

    # ---- the join key, the deleted window, the '-' that is a value, the two forms of
    # a window, the PAM side a copied results folder reads: one mutant each
    def test_reading_the_table_as_csv_breaks_the_join_rule(self):
        """csv's default dialect: one quote character anywhere in the table and a row
        spans two lines for the reader while the writer counts two -- every row after
        it then carries the next site's sequence, with nothing wrong to see."""
        with mock.patch.object(F, "table_reader",
                               lambda lines: csv.reader(lines, delimiter="\t")):
            with expect_break(self, "a quoted cell moves a row's sequence"):
                JoiningBackToTheInputTable(
                    "test_a_quote_in_a_cell_does_not_move_a_sequence_onto_another_row"
                ).debug()

    def test_an_empty_fasta_record_breaks_the_deleted_window_rule(self):
        """A header followed by a blank line, for a haplotype that holds no bases
        there: the records then outnumber the sequences and a reader pairing them by
        position is one out from there on."""
        with mock.patch.object(F, "has_sequence",
                               lambda line: line["haplotype"] != "unresolved"):
            with expect_failure(self, "a record is written with no bases"):
                AHaplotypeThatDeletesTheWindow(
                    "test_no_fasta_record_is_written_without_bases").debug()

    def test_an_empty_element_in_the_alt_list_breaks_the_one_per_haplotype_rule(self):
        """The same thing in the cell: an empty element between two sequences is
        dropped by the obvious split, and every later element is then the wrong hap."""
        real = F.annotation_rows

        def blank(out_lines):
            rows = real(out_lines)
            for cells in rows.values():
                if F.ANNOTATE_ALT_LIST in cells:
                    cells[F.ANNOTATE_ALT_LIST] = cells[F.ANNOTATE_ALT_LIST].replace(
                        F.SEQ_DELETED, "")
            return rows
        with mock.patch.object(F, "annotation_rows", blank):
            with expect_failure(self, "an ALT element is written empty"):
                AHaplotypeThatDeletesTheWindow(
                    "test_the_alt_list_keeps_one_element_per_haplotype").debug()

    def test_a_dash_read_as_missing_breaks_the_column_value_rule(self):
        """Both sides of a comparison read as numbers, whatever the expression asked
        for: Strand == "-" then selects nothing and blames missing data."""
        with mock.patch.object(
                F, "_read_pair",
                lambda a, b, state: (F._read_value(a, state), F._read_value(b, state))):
            with expect_break(self, "a '-' that is a value reads as missing"):
                RowFilters("test_a_dash_is_the_value_of_a_column_that_has_one").debug()

    def test_a_string_cell_read_as_a_value_breaks_the_missing_rule(self):
        """And the other way: the markers taken as ordinary text everywhere would let
        a row through on a score nobody measured."""
        with mock.patch.object(F._Cell, "as_value",
                               lambda self_, state: self_.number
                               if self_.number is not None else self_.text):
            with expect_break(self, "a missing score reads as its own text"):
                RowFilters("test_a_missing_value_does_not_pass_and_is_counted").debug()

    def test_refusing_a_flag_beside_a_preset_breaks_the_precedence_rule(self):
        """The refusal before this change: --preset <anything> --window W was read as
        the user giving both forms, which made the preset the docs tell a user to write
        unusable with the flag the docs say always wins."""
        real = F.resolve_window

        def same_source(upstream, downstream, window, preset, presets):
            if preset and (window is not None or upstream is not None
                           or downstream is not None):
                raise F.FlankError("give one form or the other")
            return real(upstream, downstream, window, preset, presets)
        with mock.patch.object(F, "resolve_window", same_source):
            with expect_break(self, "a flag beside a preset is refused"):
                Presets("test_a_flag_wins_over_the_presets_other_form_of_the_window"
                        ).debug()

    def test_the_pam_side_read_backwards_breaks_the_copied_folder_rule(self):
        """The branch a results folder without its PAMs/ takes, read backwards: the 3
        and the 6 then land on the opposite genomic side of every row."""
        real = F.pam_side_from_spacer

        def swapped(spacer, motif_len=None):
            answer = real(spacer, motif_len)
            if answer and motif_len:
                side, n = answer
                return (F.PAM_5PRIME if side == F.PAM_3PRIME else F.PAM_3PRIME), n
            return answer
        with mock.patch.object(F, "pam_side_from_spacer", swapped):
            with expect_failure(self, "the PAM's side is read backwards"):
                ThePamSideOfACopiedResultsFolder(
                    "test_the_run_gives_the_motif_and_the_table_gives_the_side").debug()

    def test_ignoring_the_motifs_length_breaks_the_pam_side_rule(self):
        """And the same branch never taken: with the motif's length thrown away, a
        guide with a degenerate first position can no longer be read at all."""
        real = F.pam_side_from_spacer
        with mock.patch.object(F, "pam_side_from_spacer",
                               lambda spacer, motif_len=None: real(spacer, None)):
            with expect_failure(self, "the motif's length is not used"):
                ThePamSideOfACopiedResultsFolder(
                    "test_the_side_is_read_off_the_spacer_with_the_motifs_length"
                ).debug()

    def test_an_alt_list_in_the_wrong_order_breaks_the_positional_rule(self):
        """Element k is hapK, and nothing else in the file says so."""
        real = F.annotation_rows

        def backwards(out_lines):
            rows = real(out_lines)
            for cells in rows.values():
                if F.ANNOTATE_ALT_LIST in cells:
                    cells[F.ANNOTATE_ALT_LIST] = F.ALT_LIST_SEP.join(
                        reversed(cells[F.ANNOTATE_ALT_LIST].split(F.ALT_LIST_SEP)))
            return rows
        with mock.patch.object(F, "annotation_rows", backwards):
            with expect_failure(self, "the ALT list is written backwards"):
                ColumnMode("test_the_alt_list_is_in_haplotype_order_hap1_first").debug()

    def test_alt_carriers_taken_from_the_reference_breaks_the_carrier_rule(self):
        """flank_alt_carriers must count the carriers of the sequence in flank_alt.
        The reference entry has carriers of its own: different people, different
        number, and ordering assay material off it is ordering the wrong amount."""
        real = F.annotation_rows

        def from_ref(out_lines):
            rows = real(out_lines)
            for line in out_lines:
                if line["haplotype"] == "ref" and line["row_id"] in rows:
                    rows[line["row_id"]]["flank_alt_carriers"] = line["n_carriers"]
            return rows
        with mock.patch.object(F, "annotation_rows", from_ref):
            with expect_failure(self, "the ALT carriers come from the reference"):
                ColumnMode(
                    "test_the_alt_carriers_column_counts_the_primary_alts_carriers"
                ).debug()


class WhatTheDocsPromiseAboutTheHaplotypeCap(unittest.TestCase):
    """Every surface that describes the companion files says the cap applies to them.

    README, CHANGELOG and this module's own docstring each said a row "keeps every one"
    of its distinct ALT sequences in the companion FASTA/TSV and that "none is dropped",
    with ``--max-haplotypes`` mentioned only as something that trims the ``--alt-list``
    CELL. It does not: ``site_entries`` truncates ``alts`` to ``--max-haplotypes``
    (default 10) BEFORE returning, so the TSV and the FASTA are capped too -- which is
    what ``test_the_cap_keeps_the_count_and_says_what_it_left_out`` asserts above, three
    sequences written of six at ``--max-haplotypes 2``. Someone ordering oligos off the
    companion files would have been short on exactly the busiest sites (11 rows of the
    1,000-row shortlist at the 32-mer, and most of them at ``--window 200``, where the
    median site has 22 distinct ALT sequences). METHODS § 7.1, which calls its last
    column "sequences written", was the one surface that already had it right.
    """

    # every shipped text that describes this command: the README, METHODS, the
    # changelog entry, the risk-summary document, the module's own docstring (which is
    # also what `extract-flanks --help` prints), and crisprme.py's subcommand help
    SURFACES = ("README.md", "METHODS.md", "CHANGELOG.md",
                os.path.join("docs", "offtarget_risk_summary_METHODS.md"),
                os.path.join("PostProcess", "extract_flanks.py"),
                "crisprme.py")

    def text(self, name):
        with open(os.path.join(REPO, name), encoding="utf-8") as fh:
            return fh.read()

    def test_the_promise_is_not_hiding_from_the_guard_behind_markdown(self):
        """The guard above is only worth as much as its regex. Three occurrences were
        invisible to the plain-substring version because Markdown writes the phrase as
        "one flank per **distinct observed haplotype**" and the line wraps inside it, so
        this pins that the emphasis-tolerant pattern still finds the plain form and that
        every surface the tool is described on is in the list."""
        self.assertRegex("one flank per distinct observed haplotype", self.PROMISE)
        self.assertRegex("one flank per **distinct observed\nhaplotype**", self.PROMISE)
        self.assertRegex("one flank per DISTINCT observed haplotype", self.PROMISE)
        found = {name: len(self.PROMISE.findall(self.text(name)))
                 for name in self.SURFACES}
        for name, n in found.items():
            with self.subTest(surface=name):
                self.assertTrue(n, f"{name} no longer describes extract-flanks at all")

    def test_no_surface_promises_that_every_alt_sequence_is_always_written(self):
        # assertNotIn would print the whole file -- a megabyte of README around the one
        # phrase -- and bury its own message, so this asserts on the boolean
        for name in self.SURFACES:
            with self.subTest(surface=name):
                self.assertFalse(
                    "none is dropped" in self.text(name).lower(),
                    f"{name} still promises that no ALT sequence is dropped; the "
                    f"default --max-haplotypes {F.DEFAULT_MAX_HAPLOTYPES} drops them")

    def test_every_surface_names_the_cap_and_how_to_turn_it_off(self):
        self.assertEqual(F.DEFAULT_MAX_HAPLOTYPES, 10)
        n = F.DEFAULT_MAX_HAPLOTYPES
        # the default is written both ways across these files and both are fine; what is
        # not fine is a surface that gives no default at all
        spellings = (f"default {n}", f"{n} by default", f"default is {n}")
        for name in self.SURFACES:
            body = self.text(name).lower()
            with self.subTest(surface=name):
                self.assertIn("--max-haplotypes 0", body,
                              f"{name} does not say how to write every haplotype")
                self.assertTrue(any(sp in body for sp in spellings),
                                f"{name} gives the cap no default: none of {spellings}")

    def test_the_help_a_user_reads_first_says_the_cap_applies(self):
        """`extract-flanks --help` is where most people meet this. Its opening summary
        said "one flank per distinct observed haplotype" with the cap mentioned only
        further down, among the flags."""
        help_text = F.help_text()
        self.assertIn("--max-haplotypes 0", help_text)
        self.assertIn(f"default {F.DEFAULT_MAX_HAPLOTYPES}", help_text.lower())
        # the PARAGRAPH that makes the promise has to carry the qualification -- the
        # usage block names the flag a few lines above, which is not the same thing
        paras = [p for p in re.split(r"\n\s*\n", help_text)
                 if "distinct observed haplotype" in p]
        self.assertEqual(len(paras), 1, "the promise is made in more than one place")
        self.assertIn(
            "--max-haplotypes 0", paras[0],
            "the paragraph promising one flank per distinct haplotype does not say "
            "that only --max-haplotypes of them are written")

    # the phrase as it is WRITTEN: Markdown puts ** inside it and the line wraps in the
    # middle of it, and a plain-substring version of this regex missed three occurrences
    # for exactly that reason -- README 2.2.12, METHODS 7.1 and the risk-summary document
    PROMISE = re.compile(r"flank per[\s*_`]+distinct[\s*_`]+observed[\s*_`]+haplotype",
                         re.I)

    def test_every_place_that_makes_the_promise_qualifies_it(self):
        """Per OCCURRENCE, not per file -- which is how five of seven survived.

        A file-level check ("README names --max-haplotypes somewhere") passes while the
        sentence a reader actually meets says none of it. The phrase "one flank per
        distinct observed haplotype" is the promise; wherever it is made, the cap has to
        be within reach of it. 400 characters is about a paragraph.
        """
        reach = 400
        for name in self.SURFACES:
            body = self.text(name)
            hits = list(self.PROMISE.finditer(body))
            for m in hits:
                line = body[:m.start()].count("\n") + 1
                with self.subTest(surface=name, line=line):
                    self.assertIn(
                        "--max-haplotypes", body[m.end():m.end() + reach],
                        f"{name}:{line} promises one flank per distinct haplotype and "
                        f"does not name --max-haplotypes within {reach} characters")

    def test_the_cap_really_does_reach_the_two_companion_files(self):
        """The claim above is only worth making because the code behaves this way: if
        the cap ever stops reaching the FASTA and the TSV, these docs become wrong in
        the other direction and this test says so."""
        src = self.text(os.path.join("PostProcess", "extract_flanks.py"))
        self.assertIn("alts = alts[:max_haplotypes]", src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
