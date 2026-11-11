#!/usr/bin/env python
"""The real indel post-analysis (analisi_indels_NNN.py) on SNP-overlaid targets.

A hermetic, engine-free fixture: a small reference chromosome with planted SpCas9
(NRG) sites, a phased panel (indels and SNPs, cis / trans / unphased), the fake-indel
windows cut as CRISPRitz add-variants cuts them and SNP-overlaid by
overlay_indel_snps, the engine's targets on them (IUPAC codes in the DNA column, a
code matching when any allele fits), and the dictionaries a search reads: the
Tier-0 registry, the Tier-1 genotype store, the phased indel store and the indel log,
compiled with the repo's own writers. The SHIPPED analisi_indels_NNN.py is then run
on them and its rows checked:

  * every reported base comes from the listed carrier's own haplotype: the reported
    sequence (gaps removed) is on a phased slot of every CONFIRMED carrier, at the
    reported site, and every PUTATIVE carrier holds each allele it uses;
  * every PAM satisfies NRG, and no IUPAC code is left;
  * the planted sites give exactly the rows their carriers have (a PAM SNP in cis
    and in trans, a chimera only another sample's SNP would complete, a cis spacer
    SNP, an unphased carrier, a minus-strand PAM SNP and a minus-strand spacer SNP, a
    second variant on the same slot: an indel inside the target, inserted next to the
    indel's anchor or overlapping the indel (one complex allele), a SNP on a base a
    record keeps past its POS, two indels at one POS inside the target that together
    would leave the reference (no one haplotype holds both); and one that changes
    nothing: a SNP or a nested indel on bases the indel deletes, another allele (an
    indel or a SNP) at the indel's own POS, a SNP on a kept base the target does not
    hold; several SNP alts on one VCF line; a target that lies wholly in a flank of the
    window and holds no base of the indel);
  * the SNP+indel co-occurrence companion names the SNPs a row requires, CONFIRMED
    or PUTATIVE, PAM-window SNPs included, its joint frequency counting haplotype
    copies;
  * a classic dictionary install (no Tier-1 store: the per-sample SNP dictionary), a
    sites-only panel (no genotypes), no phased indel store (which the run notes in its
    log), CRISPRME_INDEL_SNP=0 and --max-total-edits; an SNV whose ALT is written in lower
    case, which CRISPRitz does not code, is no allele on a Tier-1 or a classic install;
  * a record written on two VCF lines, in one case and in two, is one allele: an indel
    read from the store, an MNV (which no store holds) and an indel without the store
    read from the indel log, each window's rows listing the carriers of both lines.

On 2.8.0 / aa9c54b the chimera and the PAM-less rows (an R at the last PAM base read
as A) are reported, and a classic dictionary install leaves the codes in the rows.

Needs pandas and numpy (the analyzer imports them); intervaltree is stubbed when it
is not installed. Run: python -m unittest test_indel_overlaid_rows -v
"""

import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import build_indel_genotypes as big
import overlay_indel_snps as ov
import tier0_registry as t0
import tier1_genotypes as t1
from bulge_cfd_penalty import cfd_score_bulge_aware

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYZER = os.path.join(HERE, "analisi_indels_NNN.py")
CHROM = "chrT"
VCF = "synth"
GUIDE = "ACTGAAATCTGTAAGCAGGC"
CODE = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
        frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
        frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
        frozenset("ACG"): "V"}
ALLELES = {v: set(k) for k, v in CODE.items()}
COMP = str.maketrans("ACGTRYSWKMBDHVacgtryswkmbdhv", "TGCAYRSWMKVHDBtgcayrswmkvhdb")
SAMPLES = ["S%02d" % i for i in range(1, 32)]
AN = 2 * len(SAMPLES)
SAMPLE_META = {s: ("synth", "POP", "female") for s in SAMPLES}


def _gnu_wc():
    try:
        return subprocess.run(["wc", "--version"], capture_output=True).returncode == 0
    except OSError:
        return False


_SKIP_WC = "analisi_indels_NNN.py reads the line count from GNU `wc -l`, as the pipeline image has"

_STUB_INTERVALTREE = '''
class Interval(tuple):
    """Hashed and compared by (begin, end), as intervaltree's Interval is."""
    def __new__(cls, begin, end, data=None):
        return tuple.__new__(cls, (begin, end, data))
    begin = property(lambda self: self[0])
    end = property(lambda self: self[1])
    data = property(lambda self: self[2])
    def __hash__(self):
        return hash((self[0], self[1]))
    def __eq__(self, other):
        return (self[0], self[1]) == (other[0], other[1])
    def __lt__(self, other):
        return (self[0], self[1]) < (other[0], other[1])

class IntervalTree(object):
    def __init__(self):
        self._iv = []
    def __setitem__(self, key, data):
        self._iv.append(Interval(key.start, key.stop, data))
    def __getitem__(self, point):
        return set(iv for iv in self._iv if iv.begin <= point < iv.end)
'''


def rc(s):
    return s.translate(COMP)[::-1]


class World(object):
    """The reference, the panel and the planted sites."""

    def __init__(self):
        rng = random.Random(28101)
        self.seq = [rng.choice("ACGT") for _ in range(6000)]
        self.snps = []    # (pos1, ref, alt, rsid, {sample: gt})
        self.indels = []  # (pos1, ref, alt, rsid, {sample: gt})
        self.alt_token = {}   # (pos1, alt) -> its allele index, on a multi-allelic line
        self.extra_codes = {}  # pos1 -> a code the enriched genome holds with no SNP
        g = GUIDE
        # A (+): a deletion of CAC restores a perfect target, PAM TGG; a G>A SNP at the
        # PAM's last base in cis (S02), in trans (S03) and with no indel (S04).
        a0 = 100
        self._put(a0, g[:10] + "CAC" + g[10:] + "TGG")
        self._indel("A", a0 + 10, g[9] + "CAC", g[9], {"S01": "1|0", "S02": "1|0", "S03": "1|0"})
        self._snp("A", a0 + 10 + 3 + 10 + 2 + 1, "G", "A", {"S02": "1|0", "S03": "0|1", "S04": "1|0"})
        # C (+): a chimera. With its deletion S05's haplotype is 4 mismatches + one
        # more at column 6 (C); a C>A SNP there would complete a 4-mismatch target,
        # but only S06 (no deletion) carries it.
        c0 = 300
        xc = "tgac" + g[4:6] + "C" + g[7:10]
        self._put(c0, xc.upper() + "T" + g[10:] + "AGG")
        self._indel("C", c0 + 10, xc[9].upper() + "T", xc[9].upper(), {"S05": "1|0"})
        self._snp("C", c0 + 6 + 1, "C", "A", {"S06": "0|1"})
        # D (+): a C>A spacer SNP the guide needs, in cis with the deletion in S07,
        # absent in S08.
        d0 = 500
        self._put(d0, g[:3] + "C" + g[4:10] + "G" + g[10:] + "CGG")
        self._indel("D", d0 + 10, g[9] + "G", g[9], {"S07": "1|0", "S08": "0|1"})
        self._snp("D", d0 + 3 + 1, "C", "G", {"S07": "1|0", "S09": "1|0"})
        # E (+): the same, unphased: S10 carries both, phase unknown -> PUTATIVE.
        e0 = 700
        self._put(e0, g[:3] + "C" + g[4:10] + "G" + g[10:] + "TGG")
        self._indel("E", e0 + 10, g[9] + "G", g[9], {"S10": "0/1", "S11": "0/1"})
        self._snp("E", e0 + 3 + 1, "C", "G", {"S10": "0/1"})
        # F (-): an insertion of TT restores the target on the minus strand; a C>T
        # SNP at the PAM's last base (forward C, first base of 'CCT') in cis (S11) and
        # in trans (S12).
        f0 = 900
        pre = rc(g[14:] + "AGG")
        post = rc(g[:12])
        self._put(f0, pre + post)
        self._indel("F", f0 + len(pre), pre[-1], pre[-1] + "TT", {"S11": "1|0", "S12": "1|0"})
        self._snp("F", f0 + 1, "C", "T", {"S11": "1|0", "S12": "0|1"})
        # G (+): a deletion plus an A>G SNP that CREATES the PAM (TGA -> TGG), in cis in
        # S03 and in trans in S04; the co-occurrence companion must name it.
        g0 = 1100
        self._put(g0, g[:10] + "C" + g[10:] + "TGA")
        self._indel("G", g0 + 10, g[9] + "C", g[9], {"S03": "1|0", "S04": "1|0"})
        self._snp("G", g0 + 10 + 1 + 10 + 2 + 1, "A", "G", {"S03": "1|0", "S04": "0|1"})
        # O (+), no SNP: a deletion restores a perfect target (S01, S02), but S01's
        # haplotype also carries a second deletion inside it (on the same slot), so it
        # does not have that target: the one-indel fake window cannot show it.
        o0 = 1250
        self._put(o0, g[:10] + "C" + g[10:] + "TGG")
        self._indel("O", o0 + 10, g[9] + "C", g[9], {"S01": "1|0", "S02": "0|1"})
        self._indel("O2", o0 + 3, g[2] + g[3], g[2], {"S01": "1|0"})
        # P (+), no SNP: inserting the PAM after the spacer's last base (the anchor) makes
        # a perfect target that ends inside the insertion; S09's slot also holds an A
        # inserted right before the anchor (P2), inside the target, S08's the other slot.
        p0 = 1400
        self._put(p0, g + "T" * 30)
        self._indel("P", p0 + 20, g[19], g[19] + "TGG", {"S08": "1|0", "S09": "1|0"})
        self._indel("P2", p0 + 19, g[18], g[18] + "A", {"S08": "0|1", "S09": "1|0"})
        # Q (+), no SNP: the same insertion; on S12's slot a deletion of the anchor and
        # the two bases before it (Q2) overlaps it, one complex allele split into two
        # records that the window does not hold. S10 carries both unphased: a
        # multi-allelic heterozygote, one allele per haplotype, so it stays listed.
        q0 = 1550
        self._put(q0, g + "T" * 30)
        self._indel("Q", q0 + 20, g[19], g[19] + "TGG", {"S07": "0|1", "S10": "0/1", "S12": "0|1"})
        self._indel("Q2", q0 + 17, g[16:20], g[16], {"S10": "0/1", "S12": "0|1"})
        # R (+): the same insertion, with a C>T SNP at its own POS, on one slot in S05 and
        # on the other in S06: two records at one POS are alternative alleles of that base
        # (a haplotype holding both would be one record, C>TTGG), so one slot phased with
        # both is a phasing artefact, and both samples stay listed.
        r0 = 1630
        self._put(r0, g + "T" * 30)
        self._indel("R", r0 + 20, g[19], g[19] + "TGG", {"S05": "1|0", "S06": "1|0"})
        self._snp("R", r0 + 20, g[19], "T", {"S05": "1|0", "S06": "0|1"})
        # RK (+): a record left unnormalized, GCA>GC, keeps two bases; a C>T SNP on the second
        # is in cis in S05, in trans in S06, unphased in S07. The overlay leaves the ALT
        # alone, so the window shows the reference base there: S05's slot does not have
        # the window's sequence, and an unphased carrier stays.
        rk = 4200
        self._put(rk, g + "A" + "TGG" + "T" * 20)
        self._indel("RK", rk + 19, g[18:20] + "A", g[18:20],
                    {"S05": "1|0", "S06": "1|0", "S07": "1|0"})
        self._snp("RK", rk + 20, g[19], "T", {"S05": "1|0", "S06": "0|1", "S07": "0/1"})
        # S (+): a deletion restores a perfect target; S03's slot also holds a SNP on a
        # base it deletes, which is not in the window's sequence: S03 stays listed.
        s0 = 1720
        self._put(s0, g[:10] + "CAG" + g[10:] + "TGG" + "T" * 20)
        self._indel("S", s0 + 10, g[9] + "CAG", g[9], {"S03": "1|0"})
        self._snp("S", s0 + 12, "A", "T", {"S03": "1|0"})
        # T (+): a longer deletion, and on S04's slot a shorter one nested in it (T2): the
        # outer window is S04's sequence, the inner one is not.
        t0 = 1830
        self._put(t0, g[:10] + "CAGTA" + g[10:] + "TGG" + "T" * 20)
        self._indel("T", t0 + 10, g[9] + "CAGTA", g[9], {"S04": "1|0"})
        self._indel("T2", t0 + 12, "AGT", "A", {"S04": "1|0"})
        # U (+): GT>G, 2 bp, then a target a C>A SNP makes perfect, wholly in the deletion's
        # downstream flank: S13 holds both on one slot, S14 the deletion, S15 the SNP. It
        # is S13's and S15's SNP target; on the window it holds no base of the deletion.
        u0 = 2000
        self._put(u0, "GGTCA" + g[:5] + "C" + g[6:] + "TGG")
        self._indel("U", u0 + 2, "GT", "G", {"S13": "1|0", "S14": "1|0"})
        self._snp("U", u0 + 11, "C", "A", {"S13": "1|0", "S15": "1|0"})
        # V, V2 (+): two insertions at one POS and anchor, both on S16's slot 0, inside the
        # target: no haplotype holds both, so that is a phasing artefact and S16 stays on
        # each window; S17 and S18 carry one of them.
        v0 = 2200
        self._put(v0, g[:10] + g[13:] + "TGG" + "A" * 30)
        self._indel("V", v0 + 10, g[9], g[9] + g[10:13], {"S16": "1|0", "S17": "1|0"})
        self._indel("V2", v0 + 10, g[9], g[9] + "AAA", {"S16": "1|0", "S18": "0|1"})
        # V3, V4 (+): two deletions at one POS and anchor (G>., two As and one A of five), on
        # S16's slot 0: V4's window would read V3 as changing a base it covers.
        v3 = 2300
        self._put(v3, g[:4] + "AAAAA" + g[7:] + "TGG")
        self._indel("V3", v3 + 4, "GAA", "G", {"S16": "1|0", "S17": "0|1"})
        self._indel("V4", v3 + 4, "GA", "G", {"S16": "1|0", "S18": "1|0"})
        # W (+): a C>A,T SNP on one VCF line (raw genotypes under each alt), A the guide's
        # base; a deletion. S14's deletion slot holds A, S15's T (2|1, deletion on slot 0),
        # S13 is 2|2, S16 has neither.
        w0 = 2400
        self._put(w0, g[:5] + "C" + g[6:10] + "G" + g[10:] + "CGG")
        self._indel("W", w0 + 10, g[9] + "G", g[9],
                    {"S13": "1|0", "S14": "0|1", "S15": "1|0", "S16": "1|0"})
        self._snp("W", w0 + 6, "C", "A", {"S14": "0|1", "S15": "2|1"})
        self._snp("W", w0 + 6, "C", "T", {"S13": "2|2", "S15": "2|1"})
        self.alt_token[(w0 + 6, "A")] = "1"
        self.alt_token[(w0 + 6, "T")] = "2"
        # X (-): a spacer SNP the guide needs on the minus strand: forward G>T (the target
        # reads C>A), in cis in S17, in trans in S18.
        x0 = 2600
        site = g[:5] + "C" + g[6:10] + "C" + g[10:] + "TGG"
        fwd = rc(site)
        self._put(x0, fwd)
        self._indel("X", x0 + 13, fwd[12:14], fwd[12], {"S17": "1|0", "S18": "1|0"})
        self._snp("X", x0 + 19, "G", "T", {"S17": "1|0", "S18": "0|1"})
        # Y (+): a target with a SNP code; S19's slot also holds a second deletion inside
        # it (Y2), so S19 does not have the window's sequence; S20 does.
        y0 = 2800
        self._put(y0, g[:5] + "C" + g[6:10] + "C" + g[10:] + "TGG")
        self._indel("Y", y0 + 10, g[9] + "C", g[9], {"S19": "1|0", "S20": "1|0"})
        self._indel("Y2", y0 + 3, g[2] + g[3], g[2], {"S19": "1|0"})
        self._snp("Y", y0 + 6, "C", "A", {"S19": "1|0", "S20": "1|0"})
        # Z (+): a carrier homozygous for the deletion and the SNP the target needs: two
        # cis copies, one cis sample.
        z0 = 3000
        self._put(z0, g[:3] + "C" + g[4:10] + "G" + g[10:] + "CGG")
        self._indel("Z", z0 + 10, g[9] + "G", g[9], {"S21": "1|1"})
        self._snp("Z", z0 + 4, "C", "G", {"S21": "1|1"})
        # HH (+): an unphased carrier homozygous for a SNP that spoils the target (S22):
        # it has no reference allele there, so never the perfect target; S23 has no SNP.
        h0 = 3200
        self._put(h0, g[:10] + "C" + g[10:] + "TGG")
        self._indel("HH", h0 + 10, g[9] + "C", g[9], {"S22": "0/1", "S23": "0/1"})
        self._snp("HH", h0 + 4, g[3], "T", {"S22": "1/1"})
        # AA (+): a spacer SNP in cis, and on the second base a record keeps (CTG>CT) an
        # unphased SNP: one allele per haplotype, so it does not take S24 off a slot, which
        # stays exact (CONFIRMED).
        a2 = 3400
        self._put(a2, g[:3] + "C" + g[4:10] + "G" + g[10:] + "TGG")
        self._indel("AA", a2 + 9, g[8:10] + "G", g[8:10], {"S24": "1|0"})
        self._snp("AA", a2 + 4, "C", "G", {"S24": "1|0"})
        self._snp("AA", a2 + 10, g[9], "C", {"S24": "0/1"})
        # A3 (+): an insertion after two kept bases (TC>TC...), a target starting at the first
        # inserted base; S25's slot also holds a SNP on the second kept base, which the
        # target does not hold: S25 stays listed.
        a3 = 3600
        self._put(a3, "TC" + g[12:] + "TGG" + "A" * 20)
        self._indel("A3", a3 + 1, "TC", "TC" + g[:12], {"S25": "1|0", "S26": "1|0"})
        self._snp("A3", a3 + 2, "C", "A", {"S25": "1|0"})
        # B2 (+): a code column neither of whose alleles is the guide's: an unphased
        # carrier gets the allele with the higher CFD, not simply the reference base.
        b2 = 3800
        weak, strong = sorted("ACT", key=lambda b: cfd_score_bulge_aware(
            g + "NNN", g[:3] + b.lower() + g[4:] + "CGG"))[:2]
        self.b2_strong = strong
        self._put(b2, g[:3] + weak + g[4:10] + "G" + g[10:] + "CGG")
        self._indel("B2", b2 + 10, g[9] + "G", g[9], {"S27": "0/1", "S28": "1|0"})
        self._snp("B2", b2 + 4, weak, strong, {"S27": "0/1"})
        # C2 (+): a code that does not hold the reference base (never expected: the overlay
        # copies the enriched genome): the column is resolved without genotypes, and the
        # run says so.
        c2 = 4000
        self._put(c2, g[:3] + "T" + g[4:10] + "C" + g[10:] + "TGG")
        self._indel("C2", c2 + 10, g[9] + "C", g[9], {"S29": "1|0"})
        self.extra_codes[c2 + 4] = "R"
        # A second indel inside the target is read on the slot's bases: it takes a carrier
        # off only when it changes them. PR (+): a deletion restores a perfect target, PAM
        # TGG and one more G; on S12's slot a 1-bp deletion in the G run (PR2) leaves the
        # target's bases as they are (the run puts a G back); on S13's a spacer deletion
        # (PR3) changes them; S14 holds PR2 on its other slot. On S30's slot PR2 comes with
        # a T>C SNP on its anchor, the PAM's first base: two records of one base, a phasing
        # artefact, so the slot has no one sequence there and stays off; S19 holds the SNP
        # and PR alone (PAM CGG).
        pr = 4400
        self._put(pr, g[:10] + "C" + g[10:] + "TGGG" + "A")
        self._indel("PR", pr + 10, g[9] + "C", g[9], {"S11": "1|0", "S12": "1|0", "S13": "1|0",
                                                       "S14": "1|0", "S19": "1|0", "S30": "1|0"})
        self._indel("PR2", pr + 22, "TG", "T", {"S12": "1|0", "S14": "0|1", "S30": "1|0"})
        self._indel("PR3", pr + 3, g[2:4], g[2], {"S13": "1|0"})
        self._snp("PR", pr + 22, "T", "C", {"S19": "1|0", "S30": "1|0"})
        # RP (+): the same deletion; on S15's slot a CA insertion written inside the
        # target's 5' AC repeat (a record left unnormalized: bcftools norm moves it out of
        # the target) holds the target up to its last base, not from its first.
        rp = 4600
        self._put(rp, "ACAC" + g[:10] + "C" + g[10:] + "TGG")
        self._indel("RP", rp + 14, g[9] + "C", g[9], {"S15": "0|1", "S17": "1|0"})
        self._indel("RP2", rp + 5, "A", "ACA", {"S15": "0|1"})
        # TL (+): the PAM is the first unit of a TGG repeat; a left-aligned deletion of one
        # unit takes the PAM, and the next unit puts it back (S18 on one slot, S20 on both).
        # S31 holds both too, and upstream of the target an unphased deletion (TLU), which no
        # stretch read from the target's first base reaches: it does not stop that reading.
        tl = 4800
        self._put(tl, g[:10] + "C" + g[10:] + "TGG" + "TGG" + "A")
        self._indel("TL", tl + 10, g[9] + "C", g[9],
                    {"S18": "1|0", "S19": "1|0", "S20": "1|0", "S31": "1|0"})
        self._indel("TL2", tl + 21, g[19] + "TGG", g[19], {"S18": "1|0", "S20": "1|1", "S31": "1|0"})
        self._indel("TLU", tl - 8, "".join(self.seq[tl - 9:tl - 7]), self.seq[tl - 9], {"S31": "0/1"})
        # SH (+): the same, with a G>A SNP on the base the deletion brings into the PAM: on
        # S21's slot the target has no PAM; on S28's other slot it does not matter; S24's
        # SNP is unphased, so its slot's bases cannot be read: it stays off.
        sh = 5000
        self._put(sh, g[:10] + "C" + g[10:] + "TGG" + "TGG" + "A")
        self._indel("SH", sh + 10, g[9] + "C", g[9],
                    {"S21": "1|0", "S23": "1|0", "S24": "1|0", "S28": "1|0"})
        self._indel("SH2", sh + 21, g[19] + "TGG", g[19],
                    {"S21": "1|0", "S24": "1|0", "S28": "1|0"})
        self._snp("SH", sh + 27, "G", "A", {"S21": "1|0", "S24": "0/1", "S28": "0|1"})
        # CO (+): a target with a SNP code, the guide's base an alt in cis in S25, and the
        # deletion of the PAM's unit on S25's and S26's slots: both keep their rows.
        co = 5200
        self._put(co, g[:5] + "C" + g[6:10] + "C" + g[10:] + "TGG" + "TGG" + "A")
        self._indel("CO", co + 10, g[9] + "C", g[9], {"S25": "1|0", "S26": "1|0", "S27": "1|0"})
        self._indel("CO2", co + 21, g[19] + "TGG", g[19], {"S25": "1|0", "S26": "1|0"})
        self._snp("CO", co + 6, "C", "A", {"S25": "1|0"})
        # XT (+): a deletion restores a perfect target; inside it, at one POS, a tandem
        # duplication of two bases (XT2, C>CTG) and the deletion of those two bases (XT3,
        # CTG>C), each of which changes the target's bases (the review's X1). S02's slot holds
        # both with the first: two alleles of one site, which no one haplotype holds, so its
        # bases cannot be read and it stays off, as S03 (XT2) and S04 (XT3) do; S01 holds the
        # deletion alone. (At f1b75ba XT2 and XT3 added up to the reference, and S02 was listed.)
        xt = 5400
        self._put(xt, g[:13] + "T" + g[13:] + "TGG" + "A")
        self._indel("XT", xt + 13, g[12] + "T", g[12],
                    {"S01": "1|0", "S02": "1|0", "S03": "1|0", "S04": "1|0"})
        self._indel("XT2", xt + 2, g[1], g[1:4], {"S02": "1|0", "S03": "1|0"})
        self._indel("XT3", xt + 2, g[1:4], g[1], {"S02": "1|0", "S04": "1|0"})
        self.ref = "".join(self.seq)

    def _put(self, start0, s):
        for i, c in enumerate(s):
            self.seq[start0 + i] = c

    def _indel(self, site, pos1, ref, alt, gts):
        assert "".join(self.seq[pos1 - 1:pos1 - 1 + len(ref)]) == ref, (site, ref)
        self.indels.append((pos1, ref, alt, "rs_%s_indel" % site, gts))

    def _snp(self, site, pos1, ref, alt, gts):
        assert self.seq[pos1 - 1] == ref, (site, pos1, ref, self.seq[pos1 - 1])
        self.snps.append((pos1, ref, alt, "rs_%s_snp" % site, gts))

    # -- the fake-indel genome, cut and overlaid as the build does --------------------
    def fake_windows(self):
        """[(indel, start_position, window)] -- 25 bp of plain reference, the ALT, then
        27 bp downstream (CRISPRitz add-variants), SNP-overlaid (overlay_indel_snps)."""
        enriched = list(self.ref)
        alleles = {}
        for pos1, ref, alt, _rs, _g in self.snps:
            if alt != alt.upper():
                continue  # CRISPRitz codes no ALT written in lower case
            alleles.setdefault(pos1, {ref}).add(alt)
        for pos1, bases in alleles.items():
            enriched[pos1 - 1] = CODE[frozenset(bases)]
        for pos1, code in self.extra_codes.items():
            enriched[pos1 - 1] = code
        enriched = "".join(enriched)
        out = []
        for ind in self.indels:
            pos1, ref, alt = ind[:3]
            s = pos1 - 26
            sub = self.ref[s:s + 25] + alt + self.ref[s + 25 + len(ref):s + 25 + len(ref) + 27]
            out.append((ind, s, ov.overlay_sub_fasta(sub, s, ref, alt, enriched)))
        return out

    def haplotypes(self, sample):
        """The sample's two haplotypes as {slot: (sequence, real 0-based pos per base)};
        an unphased genotype's alt goes on slot 0 (the check below only reads phased
        carriers' slots)."""
        out = {}
        for h in (0, 1):
            variants = []
            for pos1, ref, alt, _rs, gts in self.snps + self.indels:
                gt = gts.get(sample)
                if gt is None:
                    continue
                toks = gt.replace("/", "|").split("|")
                if h < len(toks) and toks[h] == self.alt_token.get((pos1, alt), "1"):
                    variants.append((pos1, ref, alt))
            seq, pos, i = [], [], 0
            for pos1, ref, alt in sorted(variants):
                p0 = pos1 - 1
                if p0 < i:
                    continue  # inside a deletion already applied
                seq.extend(self.ref[i:p0])
                pos.extend(range(i, p0))
                seq.extend(alt)
                pos.extend([p0] + [None] * (len(alt) - 1) if len(alt) > len(ref)
                           else range(p0, p0 + len(alt)))
                i = p0 + len(ref)
            seq.extend(self.ref[i:])
            pos.extend(range(i, len(self.ref)))
            out[h] = ("".join(seq), pos)
        return out


def engine_targets(windows, max_mm=4, max_total=None):
    """The engine's X (no-bulge) NRG targets on the overlaid fake genome: a code
    matches the guide or the PAM when any of its alleles fits; mismatches lower
    case; past --mm or --max-total-edits no target, as in the engine. Returns (fake
    chromosome, rows)."""
    fake, rows = [], []
    starts = []
    for _ind, _s, win in windows:
        if fake:
            fake.append("N")
        starts.append(sum(len(x) for x in fake))
        fake.append(win)
    fake = "".join(fake)

    def fits(base, want):
        return want in ALLELES.get(base.upper(), {base.upper()})

    for strand in "+-":
        for f0 in range(len(fake) - 22):
            seg = fake[f0:f0 + 23]
            if "N" in seg:
                continue
            shown = seg if strand == "+" else rc(seg)
            pam = shown[20:]
            if not (fits(pam[1], "A") or fits(pam[1], "G")) or not fits(pam[2], "G"):
                continue
            dna, mm = [], 0
            for g, b in zip(GUIDE, shown[:20]):
                if fits(b, g):
                    dna.append(b.upper())
                else:
                    dna.append(b.lower())
                    mm += 1
            if mm > max_mm or (max_total is not None and mm > max_total):
                continue
            rows.append(["X", GUIDE + "NNN", "".join(dna) + pam.upper(), "fake" + CHROM,
                         str(f0), str(f0), strand, str(mm), "0", str(mm)])
    return fake, starts, rows


class _Fixture(object):
    """Write the world as a data folder and run the shipped analyzer on it."""

    def __init__(self, tmp, world, tiers=("registry", "tier1", "indelgt"), env=None,
                 max_total=None, indel_af=None, store_marker=True):
        self.tmp = tmp
        self.world = world
        d = os.path.join(tmp, "Dictionaries")
        os.makedirs(d)
        windows = world.fake_windows()
        fake, starts, rows = engine_targets(windows, max_total=max_total)
        self.rows_in = rows
        with open(os.path.join(tmp, CHROM + ".fa"), "w") as fh:
            fh.write(">%s\n%s\n" % (CHROM, world.ref))
        logd = os.path.join(d, "log_indels_" + VCF)
        os.makedirs(logd)
        with open(os.path.join(logd, "log%s.txt" % CHROM), "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            for k, ((pos1, ref, alt, rs, gts), s, win) in enumerate(windows):
                fh.write("%s_%d-%d_%d\t%s\t%s\t%s\t%s_%d_%s_%s\t%d,%d\t%s\n" % (
                    CHROM, s, s + len(win), k + 1, ",".join(sorted(gts)), rs,
                    "%.6g" % (sum(g.count("1") for g in gts.values()) / float(AN)),
                    CHROM, pos1, ref, alt, starts[k], starts[k] + len(win),
                    world.ref[s:s + len(win)]))
        recs = [(p, r, a, rs, dict(g)) for p, r, a, rs, g in world.snps]
        if "registry" in tiers:
            rd = os.path.join(d, "registry_" + VCF)
            os.makedirs(rd)
            t0.compile_registry_panel(recs, SAMPLE_META, None, t0.autosomal_ploidy,
                                      os.path.join(rd, "reg_%s.bin" % CHROM),
                                      os.path.join(rd, "reg_%s.idx" % CHROM))
        if "tier1" in tiers:
            gd = os.path.join(d, "genotypes_" + VCF)
            os.makedirs(gd)
            t1.compile_genotypes([(p, a, dict(g), r) for p, r, a, _rs, g in world.snps],
                                 t1.build_sample_axis(SAMPLE_META),
                                 os.path.join(gd, "gt_%s.bin" % CHROM),
                                 os.path.join(gd, "gt_%s.idx" % CHROM))
        if "dict" in tiers:
            dd = os.path.join(d, "dictionaries_" + VCF)
            os.makedirs(dd)
            entries = {}
            for p, r, a, rs, g in world.snps:
                entry = "%s;%s,%s;%s;%s" % (
                    ",".join("%s:%s" % kv for kv in sorted(g.items())), r, a, rs,
                    "%.6g" % (sum(x.count("1") for x in g.values()) / float(AN)))
                key = "%s,%d" % (CHROM, p)
                entries[key] = entries[key] + "$" + entry if key in entries else entry
            with gzip.open(os.path.join(dd, "my_dict_%s.json.gz" % CHROM), "wt") as fh:
                json.dump(entries, fh)
        if "indelgt" in tiers:
            igd = os.path.join(d, "indel_genotypes_" + VCF)
            os.makedirs(igd)
            with gzip.open(os.path.join(igd, "gt_indel_%s.tsv.gz" % CHROM), "wt") as fh:
                for p, r, a, _rs, g in world.indels:
                    if not big._is_indel(r, a):
                        continue  # an MNV is in no store, as compile_indel_genotypes keeps them
                    fh.write("%d_%s_%s\t%s\n" % (p, r, a, ",".join(
                        "%s:%s" % kv for kv in sorted(g.items()))))
            if store_marker:  # (as build_indel_genotypes marks a store of 2.10.0)
                with open(os.path.join(igd, "gt_indel_%s.tsv.gz.done" % CHROM), "w") as fh:
                    fh.write("%d\n%s\n" % (len(world.indels), big.FILTER_MARK))
        if indel_af is not None:
            iad = os.path.join(d, "indel_af_" + VCF)
            os.makedirs(iad)
            with gzip.open(os.path.join(iad, "indel_af_%s.tsv.gz" % CHROM), "wt") as fh:
                fh.write("#pos\tref\talt\tAF_panel\tAF_max\n")
                for p, r, a, _rs, _g in world.indels:
                    fh.write("%d\t%s\t%s\t%s\t%s\n" % (p, r, a, indel_af, indel_af))
        with open(os.path.join(tmp, "pam.txt"), "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        with open(os.path.join(tmp, "guides.txt"), "w") as fh:
            fh.write(GUIDE + "NNN\n")
        cluster = os.path.join(tmp, "targets.cluster.txt")
        with open(cluster, "w") as fh:
            fh.write("#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
                     "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\t"
                     "Annotation Type\tReal Guide\tID\tAF\tSNP\n")
            for r in rows:
                fh.write("\t".join(r + ["n", "n", "n", "n", GUIDE, "n", "n", "n"]) + "\n")
        self.out = os.path.join(tmp, "out", "fake" + CHROM)
        os.makedirs(os.path.dirname(self.out))
        argv = [sys.executable, ANALYZER, "", cluster, self.out, logd,
                os.path.join(tmp, "pam.txt"), "4", os.path.join(tmp, CHROM + ".fa"),
                os.path.join(tmp, "guides.txt"), "1", "1"]
        if max_total is not None:
            argv.append(str(max_total))
        run_env = {k: v for k, v in os.environ.items()
                   if k not in ("MAMBA_EXE", "CONDA_EXE", "CRISPRME_CONDA_EXE")}
        run_env["PATH"] = os.path.dirname(sys.executable) + ":/usr/bin:/bin"
        run_env.update(env or {})
        try:
            import intervaltree  # noqa: F401
        except ImportError:
            stub = os.path.join(tmp, "stub")
            os.makedirs(stub)
            with open(os.path.join(stub, "intervaltree.py"), "w") as fh:
                fh.write(_STUB_INTERVALTREE)
            run_env["PYTHONPATH"] = stub + os.pathsep + run_env.get("PYTHONPATH", "")
        self.proc = subprocess.run(argv, cwd=HERE, env=run_env, capture_output=True,
                                   text=True, timeout=600)

    def rows(self):
        """The CFD-ranked output rows: (crRNA, DNA, start, strand, mm, samples, rsID, AF,
        SNP, CFD)."""
        out = []
        with open(self.out + ".bestCFD_INDEL.txt") as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                out.append({"crRNA": c[1], "dna": c[2], "start": int(c[4]), "strand": c[6],
                            "mm": int(c[7]), "total": int(c[9]), "pam_gen": c[10],
                            "samples": c[12], "rsid": c[15], "af": c[16], "var": c[17],
                            "cfd": float(c[20])})
        return out

    def cooc(self):
        path = self.out + ".indel_snp_cooc.tsv"
        if not os.path.exists(path):
            return None
        with open(path) as fh:
            next(fh)
            return sorted(set(tuple(line.rstrip("\n").split("\t")) for line in fh))


class _Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_overlaid_rows_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.world = World()

    def run_analyzer(self, **kw):
        fx = _Fixture(self.tmp, self.world, **kw)
        self.assertEqual(fx.proc.returncode, 0, fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
        return fx

    def by_site(self, rows):
        """{site letter: [rows]} by the indel token each row starts with."""
        sites = {}
        for r in rows:
            pos1 = int(r["var"].split(",")[0].split("_")[1])
            site = next(rs for p, _r, _a, rs, _g in self.world.indels if p == pos1)
            sites.setdefault(site.split("_")[1], []).append(r)
        return sites

    def fake_start(self, fx, site):
        """The fake-genome offset of a site's window (from the indel log)."""
        log = os.path.join(self.tmp, "Dictionaries", "log_indels_" + VCF,
                           "log%s.txt" % CHROM)
        pos1 = next(p for p, _r, _a, rs, _g in self.world.indels
                    if rs == "rs_%s_indel" % site)
        with open(log) as fh:
            for line in fh:
                c = line.split("\t")
                if c[0] != "CHR" and int(c[4].split("_")[1]) == pos1:
                    return int(c[5].split(",")[0])
        raise KeyError(site)

    def assert_on_carrier_haplotypes(self, rows, putative_ok=()):
        """Every listed carrier holds the reported sequence: phased carriers on one
        of their slots, at the reported strand; for the samples in putative_ok each
        allele the row names must be one they carry."""
        for r in rows:
            seq = r["dna"].replace("-", "").upper()
            self.assertTrue(set(seq) <= set("ACGT"), r)
            self.assertIn(seq[-2:], ("AG", "GG"), r)
            for s in r["samples"].split(","):
                if s in putative_ok:
                    continue
                haps = self.world.haplotypes(s)
                found = any(seq in (h if r["strand"] == "+" else rc(h))
                            for h, _pos in haps.values())
                self.assertTrue(found, "%s does not carry %s (%s)" % (s, r["dna"], r["var"]))


@unittest.skipUnless(_gnu_wc(), _SKIP_WC)
class FullTiers(_Base):
    """A genotyped, phased panel: registry + Tier-1 + phased indel store."""

    def test_rows_are_the_carriers_own_haplotypes(self):
        fx = self.run_analyzer()
        rows = fx.rows()
        self.assertTrue(rows)
        # S16 holds two alleles at one POS on a slot: the fixture's haplotypes apply the first
        # variant of a position only. (S07, S22, S27 and S29 are unphased or unresolved at
        # some site: their rows are pinned below.)
        self.assert_on_carrier_haplotypes(
            rows, putative_ok=("S07", "S10", "S11", "S16", "S22", "S27", "S29", "S31"))
        site = self.by_site(rows)
        g = GUIDE
        best = lambda rs: max(rs, key=lambda r: r["cfd"])  # noqa: E731
        # A: the PAM SNP in cis destroys S02's PAM; S01 and S03 keep TGG
        a = best(site["A"])
        self.assertEqual((a["dna"], a["samples"], a["cfd"]), (g + "TGG", "S01,S03", 1.0))
        self.assertFalse(any(r["dna"].endswith("A") for r in site["A"]), site["A"])
        # C: the chimera (S05's deletion + S06's SNP) is not reported
        self.assertNotIn("C", site, site.get("C"))
        # D: the cis SNP gives S07 a perfect target, named; S08 has one mismatch
        d = {r["samples"]: r for r in site["D"] if r["mm"] == 0 or r["samples"] == "S08"}
        self.assertEqual(d["S07"]["dna"], g + "CGG")
        self.assertEqual(d["S07"]["var"].split(","),
                         ["%s_510_TG_T" % CHROM, "%s_504_C_G" % CHROM])
        self.assertEqual(d["S07"]["rsid"], "rs_D_indel,rs_D_snp")
        self.assertEqual(len(d["S07"]["af"].split(",")), 2)
        self.assertEqual((d["S08"]["mm"], d["S08"]["var"]), (1, "%s_510_TG_T" % CHROM))
        # E: unphased -- the worst case of S10's own alleles, S11 never gets the SNP
        e = {r["samples"]: r for r in site["E"]}
        self.assertEqual(e["S10"]["dna"], g + "TGG")
        self.assertIn("S11", e)
        self.assertEqual(e["S11"]["dna"][3], "c")
        # F: the minus-strand PAM SNP in cis (S11) destroys the PAM
        f = best(site["F"])
        self.assertEqual((f["dna"], f["samples"], f["cfd"]), (g + "AGG", "S12", 1.0))
        # G: the PAM-creating SNP in cis gives S03 a perfect target; S04 has no PAM
        gg = site["G"]
        self.assertEqual([(r["dna"], r["samples"], r["cfd"]) for r in gg],
                         [(g + "TGG", "S03", 1.0)])
        self.assertEqual(gg[0]["var"].split(",")[1], "%s_1124_A_G" % CHROM)
        self.assertEqual(gg[0]["pam_gen"], "TGG")  # the reference has no PAM there
        # O: no SNP code, but S01's slot also holds the second deletion: S02 only
        self.assertEqual([(r["dna"], r["samples"]) for r in site["O"]], [(g + "TGG", "S02")])
        self.assertNotIn("O2", site, site.get("O2"))
        # P: the target ends inside the insertion, and S09's slot has an A inserted
        # next to the anchor: S08 only; S09 is never listed for P2's window either
        self.assertEqual([(r["dna"], r["samples"]) for r in site["P"]], [(g + "TGG", "S08")])
        self.assertFalse([r for r in site.get("P2", ()) if "S09" in r["samples"].split(",")])
        # Q: S12's slot holds the overlapping deletion: listed on neither window
        self.assertEqual([(r["dna"], r["samples"]) for r in site["Q"]],
                         [(g + "TGG", "S07,S10")])
        self.assertFalse([r for r in site.get("Q2", ()) if "S12" in r["samples"].split(",")])
        # R: a SNP at the indel's own POS is another allele: S05 stays, on the same slot
        self.assertEqual([(r["dna"], r["samples"]) for r in site["R"]],
                         [(g + "TGG", "S05,S06")])
        # RK: a SNP on a base the record keeps past its POS, on the same slot in S05: off
        self.assertEqual([(r["dna"], r["samples"]) for r in site["RK"]],
                         [(g + "TGG", "S06,S07")])
        # S: a SNP on a base the deletion removes changes nothing: S03 stays
        self.assertEqual([(r["dna"], r["samples"]) for r in site["S"]], [(g + "TGG", "S03")])
        # T: the nested deletion is inside what the outer one removes: S04 stays there,
        # and is never listed for the inner window
        self.assertEqual([(r["dna"], r["samples"]) for r in site["T"]], [(g + "TGG", "S04")])
        self.assertFalse([r for r in site.get("T2", ()) if "S04" in r["samples"].split(",")])
        # U: a target wholly in the deletion's flank holds no base of it: no indel row,
        # although the engine reported it on the window (it is S13's and S15's SNP target)
        self.assertTrue([x for x in fx.rows_in if x[2].upper().replace("M", "A") == g + "TGG"
                         and 0 <= int(x[4]) - self.fake_start(fx, "U") < 80],
                        "the engine must report the flank target")
        self.assertNotIn("U", site, site.get("U"))
        # V..V4: another allele at one POS and anchor takes no carrier off a window
        anchor = g[9]
        v = sorted((r["var"], r["dna"], r["samples"]) for r in rows
                   if r["var"].startswith("%s_2210_" % CHROM) and r["mm"] <= 2)
        self.assertEqual(v, sorted([
            ("%s_2210_%s_%s" % (CHROM, anchor, anchor + g[10:13]), g + "TGG", "S16,S17"),
            ("%s_2210_%s_%sAAA" % (CHROM, anchor, anchor),
             g[:10] + "aa" + g[12:] + "TGG", "S16,S18")]))
        v = sorted((r["var"], r["samples"]) for r in rows
                   if r["var"].startswith("%s_2304_" % CHROM))
        self.assertEqual(v, [("%s_2304_GAA_G" % CHROM, "S16,S17"),
                             ("%s_2304_GA_G" % CHROM, "S16,S18")])
        # W: several SNP alts on one VCF line: each slot gets the allele it holds
        w = sorted((r["dna"][5], r["samples"], r["var"]) for r in site["W"])
        self.assertEqual(w, [
            ("A", "S14", "%s_2410_TG_T,%s_2406_C_A" % (CHROM, CHROM)),
            ("c", "S16", "%s_2410_TG_T" % CHROM),
            ("t", "S13,S15", "%s_2410_TG_T,%s_2406_C_T" % (CHROM, CHROM))])
        # X: a minus-strand SNP is named forward, as the VCF has it; the trans carrier's
        # mismatch counts in Mismatches and in mismatches + bulges (the engine counted the
        # code as a match)
        x = {r["samples"]: r for r in site["X"]}
        self.assertEqual(x["S17"]["dna"], g + "TGG")
        self.assertEqual(x["S17"]["var"].split(","),
                         ["%s_2613_%s_%s" % (CHROM, self.world.ref[2612:2614],
                                             self.world.ref[2612]),
                          "%s_2619_G_T" % CHROM])
        self.assertEqual((x["S18"]["mm"], x["S18"]["total"]), (1, 1))
        self.assertEqual(x["S18"]["dna"][5], "c")
        raw_x = [q for q in fx.rows_in if q[2].upper() == g[:5] + "M" + g[6:] + "TGG"]
        self.assertTrue(raw_x and all(q[7] == "0" and q[9] == "0" for q in raw_x), raw_x)
        # Y: the second-indel rule on a target with a SNP code: S19 is on no row of it
        self.assertEqual([(r["dna"], r["samples"]) for r in site["Y"]], [(g + "TGG", "S20")])
        # Z: 1|1 at both: one row, one sample
        self.assertEqual([(r["dna"], r["samples"]) for r in site["Z"]], [(g + "CGG", "S21")])
        # HH: 1/1 at a SNP that spoils the target: S22 never gets the reference base
        hh = {r["samples"]: r for r in site["HH"]}
        self.assertEqual((hh["S22"]["dna"][3], hh["S22"]["mm"]), ("t", 1))
        self.assertEqual(hh["S22"]["var"].split(",")[1], "%s_3204_G_T" % CHROM)
        self.assertEqual((hh["S23"]["dna"], hh["S23"]["cfd"]), (g + "TGG", 1.0))
        # A3: a SNP on an anchor the target does not hold: S25 stays
        self.assertEqual([(r["dna"], r["samples"]) for r in site["A3"]],
                         [(g + "TGG", "S25,S26")])
        # B2: S27's worst case takes the allele with the higher CFD, S28 has the reference
        b2 = {r["samples"]: r for r in site["B2"]}
        self.assertEqual(b2["S27"]["dna"][3], self.world.b2_strong.lower())
        self.assertNotEqual(b2["S28"]["dna"][3], self.world.b2_strong.lower())
        # C2: a code without the reference base: resolved over the code, and reported
        self.assertEqual([(r["dna"], r["samples"]) for r in site["C2"]], [(g + "TGG", "S29")])
        self.assertIn("did not map onto the reference base", fx.proc.stdout)
        # PR..CO: a second indel inside the target takes a carrier off only when it changes
        # the target's bases (decided on the position alone, as at 9c6b4046, S12, S15,
        # S18, S20, S25, S26 and S28 were taken off too)
        listed = lambda rs, s: [r for r in rs if s in r["samples"].split(",")]  # noqa: E731
        self.assertEqual([r["samples"] for r in site["PR"] if r["dna"] == g + "TGG"],
                         ["S11,S12,S14"])
        self.assertFalse(listed(site["PR"], "S13"))  # its spacer deletion changes the bases
        self.assertEqual([r["samples"] for r in site["PR"] if r["dna"] == g + "CGG"], ["S19"])
        self.assertFalse(listed(site["PR"], "S30"))  # a SNP on the second indel's anchor
        rp = best(site["RP"])  # the target is on S15's slot up to its last base
        self.assertEqual((rp["dna"], rp["samples"]), (g + "TGG", "S15,S17"))
        tl = best(site["TL"])
        self.assertEqual((tl["dna"], tl["samples"]), (g + "TGG", "S18,S19,S20,S31"))
        sh = best(site["SH"])
        self.assertEqual((sh["dna"], sh["samples"]), (g + "TGG", "S23,S28"))
        self.assertFalse(listed(site["SH"], "S21") + listed(site["SH"], "S24"))
        co = {r["samples"]: r for r in site["CO"] if r["mm"] <= 1}
        self.assertEqual(sorted(co), ["S25", "S26,S27"])
        self.assertEqual(co["S25"]["dna"], g + "TGG")
        self.assertEqual(co["S25"]["var"].split(","),
                         ["%s_5210_TC_T" % CHROM, "%s_5206_C_A" % CHROM])
        self.assertEqual(co["S26,S27"]["dna"], g[:5] + "c" + g[6:] + "TGG")
        # XT: an insertion and a deletion at one POS on S02's slot, added together, are the
        # reference, but no one haplotype holds both: S02 is not listed
        self.assertEqual([(r["dna"], r["samples"]) for r in site["XT"]], [(g + "TGG", "S01")])
        # every recomputed CFD, mismatch count and total agrees with its alignment
        for r in rows:
            self.assertAlmostEqual(cfd_score_bulge_aware(r["crRNA"], r["dna"]), r["cfd"],
                                   places=3)
            self.assertEqual(r["mm"], sum(1 for c in r["dna"][:20] if c.islower()))
            self.assertEqual(r["total"], r["mm"])  # the fixture's targets have no bulge

    def test_the_cooccurrence_companion(self):
        fx = self.run_analyzer()
        cooc = fx.cooc()
        got = sorted(set((c[1], c[6], c[7], c[8], c[10], c[11]) for c in cooc))
        self.assertEqual(got, sorted([
            ("1110", "1124", "rs_G_snp", "CONFIRMED", "1", "S03"),  # a PAM-window SNP
            ("510", "504", "rs_D_snp", "CONFIRMED", "1", "S07"),
            ("710", "704", "rs_E_snp", "PUTATIVE", "1", "S10"),
            ("2410", "2406", "rs_W_snp", "CONFIRMED", "1", "S14"),  # one alt of C>A,T
            ("2613", "2619", "rs_X_snp", "CONFIRMED", "1", "S17"),
            ("2810", "2806", "rs_Y_snp", "CONFIRMED", "1", "S20"),
            ("3010", "3004", "rs_Z_snp", "CONFIRMED", "1", "S21"),
            ("3409", "3404", "rs_AA_snp", "CONFIRMED", "1", "S24"),  # its kept-base SNP is /
            ("4010", "4004", ".", "PUTATIVE", "1", "S29"),  # the code without a record
            ("5210", "5206", "rs_CO_snp", "CONFIRMED", "1", "S25"),  # and a second indel
        ]))
        # the joint frequency counts cis copies: 1|1 is two, on the panel's alleles
        z = [c for c in cooc if c[1] == "3010"]
        self.assertEqual(sorted(set(c[9] for c in z)), ["%.6g" % (2.0 / AN)])
        d = [c for c in cooc if c[1] == "510"]
        self.assertEqual(sorted(set(c[9] for c in d)), ["%.6g" % (1.0 / AN)])

    def test_the_total_edits_cap(self):
        """With --max-total-edits 0 only the perfect haplotypes stay."""
        fx = self.run_analyzer(max_total=0)
        rows = fx.rows()
        self.assertTrue(rows)
        self.assertTrue(all(r["total"] == 0 for r in rows), rows)

    def test_an_indel_written_in_two_cases_is_one_record(self):
        """C-R1-1: site A's deletion written on two VCF lines in two cases, S01 on the
        upper-case line and S02, S03 on the lower-case one (HPRC's vcfwave writes 140
        indels of chr22 so: TG>T beside tg>t). CRISPRitz cuts a fake-indel window for each
        line and its log names each as written; the record is one allele, so both windows'
        rows name it in upper case with the carriers of both lines (S02's PAM SNP is in cis,
        S03's in trans), and the merge folds the two copies into one row. Read as written,
        the upper-case window's rows listed S01 alone and the lower-case one's, named
        chrT_110_gcac_g, S03 alone: one target's carriers split over two rows."""
        w = self.world
        i = next(k for k, ind in enumerate(w.indels) if ind[3] == "rs_A_indel")
        pos1, ref, alt, _rs, _gts = w.indels[i]
        w.indels[i] = (pos1, ref, alt, "rs_A_indel", {"S01": "1|0"})
        w.indels.insert(i + 1, (pos1, ref.lower(), alt.lower(), "rs_A_indel_lc",
                                {"S02": "1|0", "S03": "1|0"}))
        fx = self.run_analyzer()
        rows = [r for r in fx.rows() if int(r["var"].split(",")[0].split("_")[1]) == pos1]
        self.assertTrue(rows)
        self.assertEqual({r["var"].split(",")[0] for r in rows},
                         {"%s_%d_%s_%s" % (CHROM, pos1, ref, alt)})
        best = [r for r in rows if r["cfd"] == 1.0]
        self.assertEqual(sorted((r["dna"], r["samples"], r["rsid"]) for r in best),
                         [(GUIDE + "TGG", "S01,S03", "rs_A_indel"),
                          (GUIDE + "TGG", "S01,S03", "rs_A_indel_lc")])

    def test_an_mnv_on_two_lines_is_one_record(self):
        """C-R2-1: an MNV M (TT>GA, the guide's bases 4-5, a perfect target with PAM TGG)
        written on two VCF lines, in one case and in two: S01 on both (1|0, 0|1), S02 on the
        first (0|1), S03 on the second (1|0). An SNV at the target's base 13 (A>C), S04's
        alone, puts a code in M's window, so its target is resolved per carrier
        (_overlaid_rows). The phased indel store holds no MNV, so the carriers are the indel
        log's, phase unknown, each line's window listing its own; the two lines are one
        allele, so each window's row lists the carriers of both. The carriers were cached by
        the record's name, which both windows share (in one case, and since the name is
        read in upper case in two): the second window read the first's, and S03 was on no
        row of M."""
        g = GUIDE
        m0 = 5600
        for lower in (False, True):
            with self.subTest(second_line="tt>ga" if lower else "TT>GA"):
                w = World()
                w._put(m0, g[:3] + "TT" + g[5:] + "TGG")
                w.ref = "".join(w.seq)
                w.indels.append((m0 + 4, "TT", g[3:5], "rs_M_mnv", {"S01": "1|0", "S02": "0|1"}))
                ref2, alt2 = ("tt", g[3:5].lower()) if lower else ("TT", g[3:5])
                w.indels.append((m0 + 4, ref2, alt2, "rs_M_mnv_2", {"S01": "0|1", "S03": "1|0"}))
                w._snp("M", m0 + 13, g[12], "C", {"S04": "1|0"})
                tmp = tempfile.mkdtemp(prefix="indel_overlaid_rows_mnv_")
                self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
                fx = _Fixture(tmp, w)
                self.assertEqual(fx.proc.returncode, 0,
                                 fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
                self.assertTrue([x for x in fx.rows_in if "M" in x[2]],
                                "the SNV must put a code in the MNV's targets")
                rows = [r for r in fx.rows()
                        if int(r["var"].split(",")[0].split("_")[1]) == m0 + 4]
                self.assertEqual({r["var"].split(",")[0] for r in rows},
                                 {"%s_%d_TT_GA" % (CHROM, m0 + 4)})
                best = [r for r in rows if r["cfd"] == 1.0]
                self.assertEqual(sorted((r["dna"], r["samples"], r["rsid"]) for r in best),
                                 [(g + "TGG", "S01,S02,S03", "rs_M_mnv"),
                                  (g + "TGG", "S01,S02,S03", "rs_M_mnv_2")])


@unittest.skipUnless(_gnu_wc(), _SKIP_WC)
class OtherInstalls(_Base):
    def test_classic_dictionary_install(self):
        """No Tier-1 store and no registry: the per-sample SNP dictionary gives the
        genotypes; the rows match the full tiers' rows."""
        full = self.run_analyzer().rows()
        tmp2 = tempfile.mkdtemp(prefix="indel_overlaid_rows_dict_")
        self.addCleanup(shutil.rmtree, tmp2, ignore_errors=True)
        fx = _Fixture(tmp2, self.world, tiers=("dict", "indelgt"))
        self.assertEqual(fx.proc.returncode, 0, fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
        rows = fx.rows()
        key = lambda r: (r["dna"], r["start"], r["strand"], r["samples"], r["var"], r["cfd"])  # noqa: E731
        self.assertEqual(sorted(map(key, rows)), sorted(map(key, full)))
        self.assertIsNone(fx.cooc())  # the companion needs the registry, as before

    def test_no_phased_indel_store(self):
        """The indel log's carriers, no phase: a sample holding the SNP gets its worst
        case, PUTATIVE, on whichever slot the SNP is (S07 is 1|0 at the SNP). The run says
        so in its log, once per chromosome (the real-data check refuses such a run: it
        models the store's phased reading); with the store it says nothing
        (test_a_store_of_an_earlier_release_is_noted)."""
        fx = self.run_analyzer(tiers=("registry", "tier1"))
        self.assertEqual(fx.proc.stdout.count(
            "NOTE [indel-snp]: no phased indel genotype store of %s on %s:" % (VCF, CHROM)), 1,
            fx.proc.stdout[-2000:])
        rows = fx.rows()
        d = {r["samples"]: r for r in self.by_site(rows)["D"]}
        self.assertEqual(d["S07"]["dna"], GUIDE + "CGG")
        self.assertEqual(d["S08"]["dna"][3], "c")
        lines = [c for c in fx.cooc() if c[1] == "510"]
        self.assertEqual(sorted(set(c[8] for c in lines)), ["PUTATIVE"])

    def test_an_indel_on_two_lines_without_the_store_is_one_record(self):
        """C-R2-1's sweep: with no phased indel store the carriers of every indel are the
        log's, as an MNV's are. Site A's deletion on two VCF lines, in one case and in two:
        S01 on the first, S02 and S03 on the second (their PAM SNP phased, the deletion's
        phase unknown: each gets the strongest target its own alleles allow, the PAM TGG).
        Each window's row lists the carriers of both lines; read from a cache keyed by the
        record's name, the second window listed the first's (S01 alone), and S02 and S03
        were on no row of the deletion."""
        w = self.world
        i = next(k for k, ind in enumerate(w.indels) if ind[3] == "rs_A_indel")
        pos1, ref, alt, _rs, _gts = w.indels[i]
        for lower in (False, True):
            with self.subTest(second_line=(ref.lower() if lower else ref) + ">"
                              + (alt.lower() if lower else alt)):
                world = World()
                world.indels[i] = (pos1, ref, alt, "rs_A_indel", {"S01": "1|0"})
                world.indels.insert(i + 1, (pos1, ref.lower() if lower else ref,
                                            alt.lower() if lower else alt, "rs_A_indel_2",
                                            {"S02": "1|0", "S03": "1|0"}))
                tmp = tempfile.mkdtemp(prefix="indel_overlaid_rows_nostore_")
                self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
                fx = _Fixture(tmp, world, tiers=("registry", "tier1"))
                self.assertEqual(fx.proc.returncode, 0,
                                 fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
                rows = [r for r in fx.rows()
                        if int(r["var"].split(",")[0].split("_")[1]) == pos1]
                self.assertEqual({r["var"].split(",")[0] for r in rows},
                                 {"%s_%d_%s_%s" % (CHROM, pos1, ref, alt)})
                best = [r for r in rows if r["cfd"] == 1.0]
                self.assertEqual(sorted((r["dna"], r["samples"], r["rsid"]) for r in best),
                                 [(GUIDE + "TGG", "S01,S02,S03", "rs_A_indel"),
                                  (GUIDE + "TGG", "S01,S02,S03", "rs_A_indel_2")])

    def test_an_snv_alt_written_in_lower_case_is_not_read(self):
        """C-R3-1 (the review's world WLOWALT): a deletion R (TC>T) makes a target that an
        SNV q makes perfect (C>G, the guide's base 4; the reference C). q is written on two
        VCF lines, C>G (S01, S03) and C>g (S04, S05). CRISPRitz codes no ALT written in
        lower case: its genome codes q from C>G alone, and its dictionary lists both lines.
        S05 holds R on both slots and g on the first, which reads the reference C there, so
        R's perfect target is S01's alone (R and G on its second slot). A Tier-1 install's
        indel path reads the ACGT alts only; a classic install's reader of the dictionary
        read the ALT in upper case, g as G, and listed S05 on the perfect target too."""
        g = GUIDE
        m0 = 5600
        w = World()
        w._put(m0, g[:3] + "C" + g[4:10] + "C" + g[10:] + "TGG")
        w.ref = "".join(w.seq)
        w.indels.append((m0 + 10, g[9] + "C", g[9], "rs_R_indel",
                         {"S01": "1|1", "S02": "1|0", "S03": "1|0", "S04": "0|1", "S05": "1|1"}))
        w._snp("R", m0 + 4, "C", "G", {"S01": "0|1", "S03": "0|1"})
        w._snp("R", m0 + 4, "C", "g", {"S04": "1|0", "S05": "1|0"})
        got = {}
        for install, tiers in (("tier1", ("registry", "tier1", "indelgt")),
                               ("classic", ("dict", "indelgt"))):
            tmp = tempfile.mkdtemp(prefix="indel_overlaid_rows_lc_alt_%s_" % install)
            self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
            fx = _Fixture(tmp, w, tiers=tiers)
            self.assertEqual(fx.proc.returncode, 0,
                             fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
            got[install] = sorted(
                (r["dna"], r["samples"], r["cfd"]) for r in fx.rows()
                if int(r["var"].split(",")[0].split("_")[1]) == m0 + 10)
            best = [(dna, samples) for dna, samples, cfd in got[install] if cfd == 1.0]
            self.assertEqual(best, [(g + "TGG", "S01")], install)
        self.assertEqual(got["classic"], got["tier1"])

    def test_a_store_of_an_earlier_release_is_noted(self):
        """A phased indel store written before 2.10.0 (its .done marker without the FILTER
        rule) may hold the indels of records whose FILTER the search skips: the search says
        so in its log, once per chromosome, and reads it as before; one of 2.10.0 says
        nothing."""
        note = ("NOTE [indel-snp]: the phased indel genotype store of %s on %s was written "
                "before CRISPRme 2.10.0" % (VCF, CHROM))
        fx = self.run_analyzer()
        self.assertNotIn("NOTE [indel-snp]", fx.proc.stdout)
        tmp2 = tempfile.mkdtemp(prefix="indel_overlaid_rows_old_store_")
        self.addCleanup(shutil.rmtree, tmp2, ignore_errors=True)
        old = _Fixture(tmp2, self.world, store_marker=False)
        self.assertEqual((old.proc.returncode, old.proc.stderr), (0, ""))
        self.assertEqual(old.proc.stdout.count(note), 1, old.proc.stdout[-2000:])
        key = lambda r: (r["dna"], r["start"], r["strand"], r["samples"], r["var"], r["cfd"])  # noqa: E731
        self.assertEqual(sorted(map(key, old.rows())), sorted(map(key, fx.rows())))

    def test_sites_only_joint_frequency(self):
        """No genotypes: the joint frequency of a PAM-creating SNP and the indel is bounded
        by the rarer of the two (the indel's AF from its per-dataset sidecar)."""
        fx = self.run_analyzer(tiers=("registry",), indel_af="0.002")
        lines = [c for c in fx.cooc() if c[1] == "1110"]
        self.assertTrue(lines)
        self.assertEqual(sorted(set((c[8], c[9], c[10], c[11]) for c in lines)),
                         [("PUTATIVE", "0.002", "NA", "NA")])

    def test_sites_only_panel(self):
        """No genotypes at all: one PUTATIVE worst-case row per target, the log's
        carriers, a valid PAM, no code left."""
        fx = self.run_analyzer(tiers=("registry",))
        rows = fx.rows()
        self.assertTrue(rows)
        for r in rows:
            seq = r["dna"].replace("-", "").upper()
            self.assertTrue(set(seq) <= set("ACGT"), r)
            self.assertIn(seq[-2:], ("AG", "GG"), r)
        a = max(self.by_site(rows)["A"], key=lambda r: r["cfd"])
        self.assertEqual((a["dna"], a["samples"]), (GUIDE + "TGG", "S01,S02,S03"))
        self.assertTrue(all(c[8] == "PUTATIVE" and c[10] == "NA" for c in fx.cooc() or []))

    def test_opting_out_reads_the_reference(self):
        """CRISPRME_INDEL_SNP=0 on an overlaid genome: each code is the reference
        base, as on a fake-indel genome without the overlay. A target that holds no base
        of its indel gives no row there either: that rule does not depend on the gate."""
        fx = self.run_analyzer(env={"CRISPRME_INDEL_SNP": "0"})
        rows = fx.rows()
        site = self.by_site(rows)
        a = max(site["A"], key=lambda r: r["cfd"])
        self.assertEqual((a["dna"], a["samples"], a["var"]),
                         (GUIDE + "TGG", "S01,S02,S03", "%s_110_TCAC_T" % CHROM))
        for r in rows:
            self.assertEqual(len(r["var"].split(",")), 1, r)  # no SNP is named
        self.assertIsNone(fx.cooc())
        # U: the engine reports the target in the deletion's flank; it is no indel row
        self.assertTrue([x for x in fx.rows_in
                         if x[2].upper().replace("M", "A") == GUIDE + "TGG"
                         and 0 <= int(x[4]) - self.fake_start(fx, "U") < 80],
                        "the engine must report the flank target")
        self.assertNotIn("U", site, site.get("U"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
