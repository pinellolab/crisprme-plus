#!/usr/bin/env python3
"""Synthetic oracle for SNP+indel off-targets, checked through the real CLI.

  oracle.py generate <workdir>
      A small reference chromosome, a phased VCF (with a few unphased calls), its
      samplesID file, the NRG PAM, the guide and spec.json, with SpCas9 sites that
      need an indel, a SNP, or both, planted in cis, in trans, unphased, as two
      indels on one haplotype, as a multi-allelic SNP on one VCF line, with a SNP at
      the far end of the target, with a target in an indel's flank that does not
      hold the indel, and with an equal-length multi-base record (an MNV).
  oracle.py check <workdir> <PostProcess dir> <results dir> [<results dir> ...]
      Builds every sample's two haplotypes from the VCF (a brute-force truth that
      shares no code with CRISPRme's resolution) and checks every row whose variant
      list holds an indel or an MNV, in each results dir's integrated results AND its
      alternative alignments, all three blocks:
        * the PAM is a concrete NRG and no IUPAC code is left;
        * CFD, recomputed from the row's own alignment, is the reported CFD;
        * the mismatch count is the alignment's, and mismatches + bulges is their sum;
        * the REF text is the reference at the reported start, gapped as the target;
        * every variant it names is a record of the VCF, forward strand;
        * every listed carrier really has the reported sequence: a phased carrier on
          one of its haplotypes, at the reported site, spanning the indel the row
          names, with every named variant on that haplotype; an unphased one carries
          every allele the row names;
      then the planted sites' expectations (who has which target, at what CFD, under
      which variants), and the SNP+indel co-occurrence companion's CONFIRMED lines
      (their carriers hold the indel and the named SNP alts on one haplotype, their
      joint frequency counts haplotype copies). Exit status 1 on any failure.

Up to 2.8.0 the indel post-analysis read each IUPAC code of the SNP-overlaid
fake-indel genome column by column, without genotypes: 5 of the first 7 sites below
came out wrong (sites A, B, D a PAM read as TGA/AGA, sites E and E2 phantom 0-edit
targets no sample has, and site F missing a carrier's real 1-mismatch target).

Stdlib only, plus the release's own CFD scorer (PostProcess/bulge_cfd_penalty.py).
"""

import csv
import glob
import io
import json
import os
import random
import sys
import zipfile

G = "ACTGAAATCTGTAAGCAGGC"  # the guide
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
CODES = set("RYSWKMBDHVryswkmbdhv")
NSAMPLES = 53


def rc(s):
    return s.translate(COMP)[::-1]


def other(b, avoid=()):
    for c in "ACGT":
        if c != b and c not in avoid:
            return c
    raise ValueError(b)


# --------------------------------------------------------------------------- #
# generate
# --------------------------------------------------------------------------- #
def generate(wd):
    random.seed(7)
    length = 16000
    seq = [random.choice("ACGT") for _ in range(length)]
    variants = []  # dicts: site, kind, pos1, ref, alt, gts{sample: 'a|b'}
    samples = ["S%02d" % i for i in range(1, NSAMPLES + 1)]
    expect = {}

    def put(start0, s):
        for i, c in enumerate(s):
            seq[start0 + i] = c

    def var(site, kind, pos1, ref, alt, gts):
        assert "".join(seq[pos1 - 1:pos1 - 1 + len(ref)]) == ref, (site, kind, pos1, ref)
        variants.append(dict(site=site, kind=kind, pos1=pos1, ref=ref, alt=alt, gts=gts))

    def multi(site, kind, pos1, ref, alts, raw):
        """One VCF line with several alts. The truth holds one biallelic variant per alt
        (a slot holds this alt or not); the VCF line keeps the raw genotypes."""
        for k, alt in enumerate(alts):
            one = {}
            for s, g in raw.items():
                sep = "/" if "/" in g else "|"
                toks = g.split(sep)
                if str(k + 1) in toks:
                    one[s] = sep.join("1" if t == str(k + 1) else "0" for t in toks)
            var(site, kind, pos1, ref, alt, one)
            variants[-1]["line"] = dict(alts=alts, raw=raw, index=k)

    # A (+): a deletion of CAC restores a perfect target (PAM TGG); a G>A SNP at the
    # PAM's last base in cis (S02: no PAM), in trans (S03) and with no indel (S04).
    a0 = 1000
    put(a0, G[:10] + "CAC" + G[10:] + "TGG")
    var("A", "del", a0 + 10, G[9] + "CAC", G[9], {"S01": "1|0", "S02": "1|0", "S03": "1|0"})
    var("A", "pamsnp", a0 + 10 + 3 + 10 + 2 + 1, "G", "A",
        {"S02": "1|0", "S03": "0|1", "S04": "1|0"})
    expect["A"] = {"perfect": ["S01", "S03"], "never": ["S02"]}

    # B (-): a 2-bp insertion restores the target on the minus strand; a C>T SNP at
    # the PAM's last base (forward C) in cis (S06) and in trans (S07).
    b0 = 2000
    pre = rc(G[14:] + "AGG")
    post = rc(G[:12])
    put(b0, pre + post)
    var("B", "ins", b0 + len(pre), pre[-1], pre[-1] + "TT",
        {"S05": "1|0", "S06": "1|0", "S07": "1|0"})
    var("B", "pamsnp", b0 + 1, "C", "T", {"S06": "1|0", "S07": "0|1", "S08": "1|0"})
    expect["B"] = {"perfect": ["S05", "S07"], "never": ["S06"]}

    # C (+): a control deletion with no SNP.
    c0 = 3000
    put(c0, G[:6] + "GG" + G[6:] + "CGG")
    var("C", "del", c0 + 6, G[5] + "GG", G[5], {"S09": "1|0"})
    expect["C"] = {"perfect": ["S09"]}

    # D (+): a PAM-creating A>G SNP (TGA -> TGG) in cis with a deletion (S10).
    d0 = 4000
    put(d0, G[:10] + "C" + G[10:] + "TGA")
    var("D", "del", d0 + 10, G[9] + "C", G[9], {"S10": "1|0", "S11": "1|0"})
    var("D", "pamsnp", d0 + 10 + 1 + 10 + 2 + 1, "A", "G",
        {"S10": "1|0", "S11": "0|1", "S12": "1|1"})
    expect["D"] = {"perfect": ["S10"], "never": ["S11"], "cooc": ["S10"]}

    def dense(site, start0, pam, del_sample, snp_a, snp_b, pam_snp=None):
        # X differs from G at every position but 9; an extra C after X[9] is deleted
        # in del_sample; every other position has a SNP whose ALT is the guide's base,
        # split between two samples that do not carry the deletion.
        x = []
        for i, g in enumerate(G):
            if i == 9:
                x.append(g)
            elif i == 10:
                x.append(other(g, avoid=("C",)))
            else:
                x.append(other(g))
        x = "".join(x)
        put(start0, x[:10] + "C" + x[10:] + pam)
        var(site, "del", start0 + 10, x[9] + "C", x[9], {del_sample: "1|0"})
        for i in range(20):
            if i == 9:
                continue
            off = i if i < 10 else i + 1
            var(site, "snp%02d" % i, start0 + off + 1, x[i], G[i],
                {(snp_a if i < 10 else snp_b): "1|0"})
        if pam_snp:
            var(site, "pamsnp", start0 + 21 + 2 + 1, pam[2], pam_snp[0], {pam_snp[1]: "0|1"})

    # E, E2 (+): every protospacer base carries a SNP whose ALT is the guide base, but
    # no haplotype carries them together, and none with the deletion: no target.
    dense("E", 5000, "TGG", "S13", "S14", "S15")
    dense("E2", 6000, "AGG", "S16", "S14", "S15", pam_snp=("A", "S14"))
    expect["E"] = {"none_above": 0.5}
    expect["E2"] = {"none_above": 0.5}

    # F (+): a protospacer C>A SNP the guide needs, in cis with a deletion (S17:
    # perfect) and in trans (S18: one mismatch, CFD 0.4545).
    f0 = 7000
    put(f0, G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG")
    var("F", "del", f0 + 10, G[9] + "C", G[9], {"S17": "1|0", "S18": "1|0"})
    var("F", "snp05", f0 + 6, "C", "A", {"S17": "1|0", "S18": "0|1"})
    expect["F"] = {"perfect": ["S17"], "one_mm": ["S18"], "cooc": ["S17"]}

    # U (+): F's geometry with unphased calls: S19 carries the deletion and the SNP,
    # phase unknown (PUTATIVE, the worst case of its own alleles); S20 only the
    # deletion (one mismatch).
    u0 = 8000
    put(u0, G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG")
    var("U", "del", u0 + 10, G[9] + "C", G[9], {"S19": "0/1", "S20": "0/1"})
    var("U", "snp05", u0 + 6, "C", "A", {"S19": "0/1"})
    expect["U"] = {"perfect": ["S19"], "one_mm": ["S20"], "cooc_putative": ["S19"]}

    # O (+): two deletions inside one target. O1 alone leaves a 1-DNA-bulge target
    # (S22); O1 + O2 on one haplotype (S21) make a perfect one that a one-indel
    # fake-indel window cannot hold, so S21 is listed on neither window's rows.
    o0 = 9000
    put(o0, G[:5] + "T" + G[5:10] + "C" + G[10:] + "TGG")
    var("O", "del2", o0 + 5, G[4] + "T", G[4], {"S21": "1|0"})
    var("O", "del", o0 + 11, G[9] + "C", G[9], {"S21": "1|0", "S22": "1|0"})
    expect["O"] = {"listed": ["S22"], "never": ["S21"], "never_any": ["S21"]}

    # M (+), M2 (-): the indel restores a perfect target. S24's and S26's slot also holds
    # a SNP at the target's far-end base, so their target has one mismatch there. An
    # RNA-bulge alignment at the same site does not reach that base and lists them too,
    # under the same variant list: the merge step must not put them on the perfect row.
    m0 = 10000
    put(m0, G[:10] + "C" + G[10:] + "TGG")
    var("M", "del", m0 + 10, G[9] + "C", G[9], {"S23": "1|0", "S24": "1|0"})
    var("M", "farsnp", m0 + 1, G[0], other(G[0]), {"S24": "1|0"})
    expect["M"] = {"perfect": ["S23"], "never_perfect": ["S24"], "one_mm": ["S24"]}
    m2 = 10300
    pre = rc(G[14:] + "AGG")
    post = rc(G[:12])
    put(m2, pre + post)
    var("M2", "ins", m2 + len(pre), pre[-1], pre[-1] + "TT", {"S25": "1|0", "S26": "1|0"})
    far = m2 + len(pre) + len(post) - 1
    var("M2", "farsnp", far + 1, seq[far], "G", {"S26": "1|0"})  # a CFD below 1 there
    expect["M2"] = {"perfect": ["S25"], "never_perfect": ["S26"], "one_mm": ["S26"]}

    # N (+): F's geometry with the spacer SNP as one VCF line C>A,T; A is the guide's base.
    # S28's deletion slot holds A (perfect), S29's holds T (2|1, its deletion on slot 1),
    # S27 is 2|2 (no A, no C), S30 has no SNP (C).
    n0 = 10600
    put(n0, G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG")
    var("N", "del", n0 + 10, G[9] + "C", G[9],
        {"S27": "1|0", "S28": "0|1", "S29": "1|0", "S30": "1|0"})
    multi("N", "snp05", n0 + 6, "C", ["A", "T"], {"S27": "2|2", "S28": "0|1", "S29": "2|1"})
    expect["N"] = {"perfect": ["S28"], "never_perfect": ["S27", "S29", "S30"],
                   "one_mm": ["S27", "S29", "S30"],
                   "named": {"S28": ["chr21_%d_C_A" % (n0 + 6)],
                             "S27": ["chr21_%d_C_T" % (n0 + 6)],
                             "S29": ["chr21_%d_C_T" % (n0 + 6)]},
                   "bare": ["S30"]}

    # N2 (+): the same with S31 the only carrier of A, 2|1: its carrier genotypes alone
    # do not say which token is A; S32's T (2|2) does. Its deletion slot holds A.
    n2 = 10900
    put(n2, G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG")
    var("N2", "del", n2 + 10, G[9] + "C", G[9], {"S31": "0|1", "S32": "1|0"})
    multi("N2", "snp05", n2 + 6, "C", ["A", "T"], {"S31": "2|1", "S32": "2|2"})
    expect["N2"] = {"perfect": ["S31"], "never_perfect": ["S32"], "one_mm": ["S32"]}

    # L (+): a target in the downstream flank of a deletion, holding no base of it: GT>G,
    # then 2 bp, then a target a C>A SNP makes perfect. S33 carries both on one slot, S36
    # only the deletion, S34 and S35 only the SNP. It is a SNP off-target of S33, S34,
    # S35; the fake-indel window's row (deletion + SNP, S33) must not take its place.
    l0 = 11200
    put(l0, "GGTCA" + G[:5] + "C" + G[6:] + "TGG")
    var("L", "del", l0 + 2, "GT", "G", {"S33": "1|0", "S36": "1|0"})
    var("L", "snp05", l0 + 5 + 5 + 1, "C", "A", {"S33": "1|0", "S34": "1|0", "S35": "0|1"})
    expect["L"] = {"no_indel_rows": True, "snp_perfect": ["S33", "S34", "S35"]}

    # L2 (+): the same on an insertion: A>ACC, then 3 bp, then a perfect reference target.
    # L3 (+): a perfect reference target whose last base is the anchor of an insertion:
    # the target holds no inserted base, so its sequence is the reference's.
    l2 = 11500
    put(l2, "GATCA" + G + "TGG")
    var("L2", "ins", l2 + 2, "A", "ACC", {"S33": "0|1", "S36": "0|1"})
    expect["L2"] = {"no_indel_rows": True}
    l3 = 11800
    put(l3, G + "TGG" + "TTCAGTTAC")
    var("L3", "ins", l3 + 23, "G", "GCAT", {"S34": "0|1", "S35": "1|0"})
    expect["L3"] = {"no_indel_rows": True}

    # R (-): F's geometry on the minus strand. The forward SNP is G>T (the target reads
    # C>A), in cis in S37, in trans in S38: named forward, as the VCF has it.
    r0 = 12100
    site = G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG"
    put(r0, rc(site))
    fwd = rc(site)
    # the deleted C is forward index 13 (G), its anchor the base left of it
    var("R", "del", r0 + 13, fwd[12:14], fwd[12], {"S37": "1|0", "S38": "1|0"})
    var("R", "snp05", r0 + 19, "G", "T", {"S37": "1|0", "S38": "0|1"})
    expect["R"] = {"perfect": ["S37"], "one_mm": ["S38"], "cooc": ["S37"],
                   "named": {"S37": ["chr21_%d_G_T" % (r0 + 19)]}, "bare": ["S38"]}

    # Q (+): a target holding a SNP code and, on S39's slot, a second deletion inside it.
    # The window keeps the T its second deletion removes, so S39 does not have the
    # window's sequence and is on no row; S40, with the SNP and no second deletion, is.
    q0 = 12400
    put(q0, G[:2] + "T" + G[2:5] + "C" + G[6:10] + "C" + G[10:] + "TGG")
    var("Q", "del2", q0 + 2, G[1] + "T", G[1], {"S39": "1|0"})
    var("Q", "snp05", q0 + 7, "C", "A", {"S39": "1|0", "S40": "1|0"})
    var("Q", "del", q0 + 11, G[9] + "C", G[9], {"S39": "1|0", "S40": "1|0"})
    expect["Q"] = {"listed": ["S40"], "never_any": ["S39"]}

    # H (+): an unphased carrier homozygous for a G>T SNP that spoils the target: it has
    # no reference allele there, so it is never listed for the perfect target. S42 has
    # the deletion and no SNP.
    h0 = 12700
    put(h0, G[:10] + "C" + G[10:] + "TGG")
    var("H", "del", h0 + 10, G[9] + "C", G[9], {"S41": "0/1", "S42": "0/1"})
    var("H", "snp03", h0 + 4, "G", "T", {"S41": "1/1"})
    expect["H"] = {"perfect": ["S42"], "never_perfect": ["S41"], "one_mm": ["S41"],
                   "named": {"S41": ["chr21_%d_G_T" % (h0 + 4)]}}

    # J (+): a carrier homozygous for the deletion and the SNP the target needs: two cis
    # copies, so the joint frequency is 2 / AN, with one cis sample.
    j0 = 13000
    put(j0, G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG")
    var("J", "del", j0 + 10, G[9] + "C", G[9], {"S43": "1|1"})
    var("J", "snp05", j0 + 6, "C", "A", {"S43": "1|1"})
    expect["J"] = {"perfect": ["S43"], "cooc": ["S43"], "cooc_copies": 2}

    # I (+): a long insertion holding a whole perfect target: it spans no reference base,
    # yet it needs the indel (the window's ALT), so it is an indel row of the carrier.
    i0 = 13300
    var("I", "ins", i0 + 1, seq[i0], seq[i0] + "CAT" + G + "TGG" + "CAT", {"S44": "1|0"})
    expect["I"] = {"perfect": ["S44"]}

    # P (+): an insertion restores a perfect target; on S45's slot 0 another insertion and a
    # SNP share its POS and anchor: alternative alleles, so a phasing artefact, and S45
    # stays on each window. S46 carries the first insertion alone.
    p0 = 13600
    put(p0, G[:10] + G[13:] + "TGG")
    var("P", "ins", p0 + 10, G[9], G[9] + G[10:13], {"S45": "1|0", "S46": "0|1"})
    var("P", "ins2", p0 + 10, G[9], G[9] + "AAA", {"S45": "1|0"})
    var("P", "snp", p0 + 10, G[9], other(G[9]), {"S45": "1|0"})
    expect["P"] = {"perfect": ["S45", "S46"]}

    # K (+), K2 (-): an equal-length multi-base record (an MNV, AC>TG on the guide's
    # strand) restores a perfect target; CRISPRitz puts it in a fake-indel window, as it
    # does an indel. S48's and S50's slots both hold a SNP at the target's far-end base,
    # so their MNV haplotype has one mismatch there, whatever the MNV's phase. An
    # RNA-bulge alignment at the same site does not reach that base and lists them under
    # the MNV alone: the merge step must not put them on the perfect row.
    k0 = 13900
    put(k0, G[:9] + "AC" + G[11:] + "TGG")
    var("K", "mnv", k0 + 10, "AC", G[9:11], {"S47": "1|0", "S48": "1|0"})
    var("K", "farsnp", k0 + 1, G[0], other(G[0]), {"S48": "1|1"})
    expect["K"] = {"perfect": ["S47"], "never_perfect": ["S48"], "one_mm": ["S48"]}
    k2 = 14200
    fwd = rc(G[:9] + "AC" + G[11:] + "TGG")
    put(k2, fwd)
    var("K2", "mnv", k2 + 13, fwd[12:14], rc(G[9:11]), {"S49": "1|0", "S50": "1|0"})
    var("K2", "farsnp", k2 + 23, fwd[22], "G", {"S50": "1|1"})  # a CFD below 1 there
    # (no one-mismatch row is expected for S50: that alignment of its haplotype ties on
    # mismatches + bulges with a DNA-bulge one with none, which each file keeps)
    expect["K2"] = {"perfect": ["S49"], "never_perfect": ["S50"]}

    # W (+): a deletion restores a perfect target whose PAM TGG is followed by one more G;
    # on S51's slot a 1-bp deletion in that G run leaves the target's bases as they are
    # (the run puts a G back), so S51 has the target; on S52's slot a spacer deletion
    # changes them, and S52 has none. A second indel inside a target is read on the
    # carrier's bases, not on its position.
    w0 = 14500
    put(w0, G[:10] + "C" + G[10:] + "TGGG" + "A")
    var("W", "del", w0 + 10, G[9] + "C", G[9], {"S51": "1|0", "S52": "1|0"})
    var("W", "pamrun", w0 + 22, "TG", "T", {"S51": "1|0"})
    var("W", "spacer", w0 + 3, G[2:4], G[2], {"S52": "1|0"})
    expect["W"] = {"perfect": ["S51"], "never_any": ["S52"]}

    # V (+): the same deletion, the target's 5' end inside an AC repeat; on S53's slot a CA
    # insertion written inside the target (a record left unnormalized: bcftools norm moves
    # it upstream, out of the target) leaves the target on the haplotype up to its last
    # base, not from its first.
    v0 = 14800
    put(v0, "ACAC" + G[:10] + "C" + G[10:] + "TGG")
    var("V", "del", v0 + 14, G[9] + "C", G[9], {"S53": "0|1"})
    var("V", "ins", v0 + 5, "A", "ACA", {"S53": "0|1"})
    expect["V"] = {"perfect": ["S53"]}

    # (No split multi-allelic SNP records here: CRISPRitz add-variants keeps only the
    # last of two records at one position in the enriched genome and the dictionary,
    # so the other alt is invisible to every stage. PostProcess/
    # test_indel_haplotype_resolution.py covers a position whose source knows both.
    # No MNV carrier heterozygous for a SNP its target covers either: the phased indel
    # store holds only records whose REF and ALT differ in length, so an MNV's phase is
    # unknown to the post-analysis and such a carrier gets its PUTATIVE worst case.)

    variants.sort(key=lambda v: (v["pos1"], len(v["ref"]), v["alt"]))
    for d in ("Genomes/synth", "VCFs/synth_vcf", "samplesIDs", "PAMs"):
        os.makedirs(os.path.join(wd, d), exist_ok=True)
    ref = "".join(seq)
    with open(os.path.join(wd, "Genomes/synth/chr21.fa"), "w") as fh:
        fh.write(">chr21\n")
        for i in range(0, length, 60):
            fh.write(ref[i:i + 60] + "\n")
    with open(os.path.join(wd, "VCFs/synth_vcf/synth.chr21.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr21,length=%d>\n" % length)
        fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n')
        fh.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n')
        fh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(samples) + "\n")
        an = 2 * len(samples)
        for v in variants:
            line = v.get("line")
            if line is not None:
                if line["index"]:
                    continue  # the first alt's line holds them all
                gts = [line["raw"].get(s, "0|0") for s in samples]
                acs = [sum(g.count(str(k + 1)) for g in gts) for k in range(len(line["alts"]))]
                fh.write("chr21\t%d\tsyn_%s_%s\t%s\t%s\t.\tPASS\tAC=%s;AN=%d;AF=%s\tGT\t%s\n" % (
                    v["pos1"], v["site"], v["kind"], v["ref"], ",".join(line["alts"]),
                    ",".join(str(a) for a in acs), an,
                    ",".join("%.6g" % (a / an) for a in acs), "\t".join(gts)))
                continue
            gts = [v["gts"].get(s, "0|0") for s in samples]
            ac = sum(g.count("1") for g in gts)
            fh.write("chr21\t%d\tsyn_%s_%s\t%s\t%s\t.\tPASS\tAC=%d;AN=%d;AF=%.6g\tGT\t%s\n" % (
                v["pos1"], v["site"], v["kind"], v["ref"], v["alt"], ac, an, ac / an,
                "\t".join(gts)))
    with open(os.path.join(wd, "samplesIDs/synth.samplesID.txt"), "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for i, s in enumerate(samples):
            fh.write("%s\tPOP%d\tSUP%d\tfemale\n" % (s, i % 3, i % 2))
    with open(os.path.join(wd, "PAMs/20bp-NRG-SpCas9.txt"), "w") as fh:
        fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
    with open(os.path.join(wd, "guides.txt"), "w") as fh:
        fh.write(G + "NNN\n")
    with open(os.path.join(wd, "list_vcf.txt"), "w") as fh:
        fh.write("synth_vcf\n")
    with open(os.path.join(wd, "list_samplesID.txt"), "w") as fh:
        fh.write("synth.samplesID.txt\n")
    sites = {}
    for v in variants:
        sites.setdefault(v["site"], [v["pos1"], v["pos1"]])
        sites[v["site"]][0] = min(sites[v["site"]][0], v["pos1"])
        sites[v["site"]][1] = max(sites[v["site"]][1], v["pos1"] + len(v["ref"]))
    with open(os.path.join(wd, "spec.json"), "w") as fh:
        json.dump(dict(guide=G, samples=samples, variants=variants, expect=expect,
                       sites=sites), fh, indent=1)
    print("oracle: wrote %d variants over %d samples, reference %d bp"
          % (len(variants), len(samples), length))


# --------------------------------------------------------------------------- #
# check
# --------------------------------------------------------------------------- #
class Truth(object):
    def __init__(self, wd):
        with open(os.path.join(wd, "spec.json")) as fh:
            self.spec = json.load(fh)
        with open(os.path.join(wd, "Genomes/synth/chr21.fa")) as fh:
            self.ref = "".join(l.strip() for l in fh if not l.startswith(">")).upper()
        self.gts = {}  # sample -> [(variant, tokens, phased)]
        for v in self.spec["variants"]:
            for s, g in v["gts"].items():
                phased = "|" in g
                self.gts.setdefault(s, []).append((v, g.replace("/", "|").split("|"), phased))

    def phased(self, sample):
        return all(p for _v, _t, p in self.gts.get(sample, ()))

    def haplotype(self, sample, h, first=()):
        """(sequence, reference 0-based position of each base, None for inserted). Of two
        records at one position on a slot (a phasing artefact) the first applies, those
        in ``first`` before any other."""
        vs = sorted((v for v, toks, _p in self.gts.get(sample, ()) if toks[h] == "1"),
                    key=lambda v: (v["pos1"], (v["pos1"], v["ref"], v["alt"]) not in first))
        out, pos, i = [], [], 0
        for v in vs:
            p0 = v["pos1"] - 1
            if p0 < i:
                continue
            out.extend(self.ref[i:p0])
            pos.extend(range(i, p0))
            r, a = v["ref"], v["alt"]
            out.append(a[0])
            pos.append(p0)
            if len(a) > len(r):
                out.extend(a[1:])
                pos.extend([None] * (len(a) - 1))
            elif len(a) == len(r):
                out.extend(a[1:])
                pos.extend(range(p0 + 1, p0 + len(a)))
            i = p0 + len(r)
        out.extend(self.ref[i:])
        pos.extend(range(i, len(self.ref)))
        return "".join(out), pos

    def occurrences(self, sample, seq, strand, start, slack=40, first=()):
        """(slot, reference position of each base) of every occurrence of ``seq``
        (strand-oriented) on the sample's haplotypes, its leftmost base within ``slack``
        bp of the reference coordinate ``start``."""
        out = []
        for h in (0, 1):
            hap, pos = self.haplotype(sample, h, first)
            probe = seq if strand == "+" else rc(seq)
            k = hap.find(probe)
            while k >= 0:
                seg = pos[k:k + len(probe)]
                anchor = next((p for p in seg if p is not None), None)
                if anchor is None:  # wholly inside an insertion: the base before it
                    anchor = next((p for p in reversed(pos[:k]) if p is not None), None)
                if anchor is not None and abs(anchor - start) <= slack:
                    out.append((h, seg))
                k = hap.find(probe, k + 1)
        return out

    def carries_at(self, sample, seq, strand, start, slack=40, first=()):
        """True iff one haplotype of the sample holds ``seq`` near ``start``."""
        return bool(self.occurrences(sample, seq, strand, start, slack, first))

    def on_slot(self, sample, h, pos1, ref, alt):
        """True iff the sample's slot ``h`` holds the variant."""
        for v, toks, _p in self.gts.get(sample, ()):
            if (v["pos1"], v["ref"], v["alt"]) == (pos1, ref, alt) and toks[h] == "1":
                return True
        return False

    def holds_all(self, sample, seq, strand, start, tokens):
        """True iff the sample has ``seq`` near ``start`` on a slot holding every named
        variant, the occurrence spanning each named indel: an inserted base (no reference
        position) or a jump across the deleted bases."""
        for h, seg in self.occurrences(sample, seq, strand, start, first=tokens):
            ok = True
            for pos1, ref, alt in tokens:
                if not self.on_slot(sample, h, pos1, ref, alt):
                    ok = False
                    break
                if len(ref) > len(alt):
                    lo, hi = pos1 - 1, pos1 - 1 + len(ref)  # anchor; first base kept after it
                    if not any(a == lo and b == hi for a, b in zip(seg, seg[1:])):
                        ok = False
                        break
                elif len(alt) > len(ref) and None not in seg:
                    ok = False
                    break
            if ok:
                return True
        return False

    def is_record(self, pos1, ref, alt):
        return any((v["pos1"], v["ref"], v["alt"]) == (pos1, ref, alt)
                   for v in self.spec["variants"])

    def carries_allele(self, sample, pos1, ref, alt):
        for v, toks, _p in self.gts.get(sample, ()):
            if (v["pos1"], v["ref"], v["alt"]) == (pos1, ref, alt) and "1" in toks:
                return True
        return False

    def site_of(self, start):
        for site, (lo, hi) in self.spec["sites"].items():
            if lo - 60 <= start <= hi + 60:
                return site
        return None


BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _is_indel_token(tok):
    """A variant of the indel path: a REF or ALT longer than one base (an indel, or an
    equal-length multi-base record, an MNV, which CRISPRitz also puts in a fake-indel
    window)."""
    p = tok.rsplit("_", 3)
    return len(p) == 4 and (len(p[2]) > 1 or len(p[3]) > 1)


def _open_table(path):
    """A results table, or the one TSV inside its .zip (a run without --ci-cd-test
    keeps the alternative alignments zipped only)."""
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith(".tsv"))
            return io.StringIO(z.read(name).decode())
    return open(path, newline="")


def alt_rows(path):
    """Every (block, field dict) of a results table with a variant."""
    csv.field_size_limit(1 << 30)
    out = []
    with _open_table(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        head = next(r)
        col = {n: i for i, n in enumerate(head)}
        for row in r:
            for b in BLOCKS:
                try:
                    f = {
                        "origin": row[col["REF/ALT_origin_(%s)" % b]],
                        "start": int(row[col["Start_coordinate_(%s)" % b]]),
                        "strand": row[col["Strand_(%s)" % b]],
                        "crrna": row[col["Aligned_spacer+PAM_(%s)" % b]],
                        "ref": row[col["Aligned_protospacer+PAM_REF_(%s)" % b]],
                        "alt": row[col["Aligned_protospacer+PAM_ALT_(%s)" % b]],
                        "mm": int(row[col["Mismatches_(%s)" % b]]),
                        "bulges": int(row[col["Bulges_(%s)" % b]]),
                        "total": int(row[col["Mismatches+bulges_(%s)" % b]]),
                        "var": row[col["Variant_info_genome_(%s)" % b]],
                        "samples": row[col["Variant_samples_(%s)" % b]],
                    }
                except (KeyError, ValueError):
                    continue
                # the CRISPR-Bulge block carries no CFD of its own alignment
                key = "CFD_score_ALT_(%s)" % b
                f["cfd"] = float(row[col[key]]) if key in col else None
                if f["origin"] != "alt" or f["var"] in ("", "NA"):
                    continue
                f["block"] = b
                out.append(f)
    return out


def indel_rows(path):
    """Every (block, field dict) of a results table whose variant list holds an indel."""
    return [f for f in alt_rows(path)
            if any(_is_indel_token(t) for t in f["var"].split(","))]


def check(wd, pp, result_dirs):
    sys.path.insert(0, pp)
    from bulge_cfd_penalty import cfd_score_bulge_aware

    truth = Truth(wd)
    expect = truth.spec["expect"]
    fails = []

    def fail(msg):
        fails.append(msg)
        print("FAIL  " + msg)

    for rd in result_dirs:
        name = os.path.basename(os.path.normpath(rd))
        per_sample = "per" in name
        alt = glob.glob(os.path.join(rd, "*all_results_with_alternative_alignments.tsv"))
        tables = (glob.glob(os.path.join(rd, "*integrated_results.tsv")) + (alt or glob.glob(
            os.path.join(rd, "*all_results_with_alternative_alignments.zip"))))
        if not tables:
            fail("%s: no results table" % rd)
            continue
        rows = []
        for t in tables:
            rows.extend(indel_rows(t))
        print("oracle: %s: %d indel row blocks" % (name, len(rows)))
        snp_only = [f for t in glob.glob(os.path.join(rd, "*integrated_results.tsv"))
                    for f in alt_rows(t)
                    if not any(_is_indel_token(x) for x in f["var"].split(","))]
        seen = {}  # site -> [(samples set, cfd, mm, names)]
        for f in rows:
            where = "%s %s:%d(%s) %s %s" % (name, f["block"], f["start"], f["strand"],
                                            f["alt"], f["var"])
            stripped = f["alt"].replace("-", "")
            if set(stripped) - set("ACGTacgt"):
                fail("%s: IUPAC code left" % where)
                continue
            pam = stripped[-3:].upper()
            if not (pam[1] in "AG" and pam[2] == "G"):
                fail("%s: PAM %s is not NRG" % (where, pam))
            mm = sum(1 for c in f["alt"][:-3] if c.islower())
            if mm != f["mm"]:
                fail("%s: %d mismatches reported, the alignment has %d" % (where, f["mm"], mm))
            if f["total"] != f["mm"] + f["bulges"]:
                fail("%s: %d mismatches + bulges reported, the row has %d + %d"
                     % (where, f["total"], f["mm"], f["bulges"]))
            # the REF text is the reference at the row's start: as many bases as its target
            # holds, in its orientation, gapped as its target (an RNA bulge is a gap)
            ref_at = truth.ref[f["start"]:f["start"] + len(stripped)]
            ref_at = list(ref_at if f["strand"] == "+" else rc(ref_at))
            for i, c in enumerate(f["alt"]):
                if c == "-":
                    ref_at.insert(i, "-")
            if f["ref"].upper() != "".join(ref_at):
                fail("%s: REF %s is not the reference at the row's start (%s)"
                     % (where, f["ref"], "".join(ref_at)))
            if f["cfd"] is not None:
                cfd = cfd_score_bulge_aware(f["crrna"], f["alt"])
                if abs(round(cfd, 3) - f["cfd"]) > 1e-9:
                    fail("%s: CFD %.3f reported, the alignment scores %.3f"
                         % (where, f["cfd"], cfd))
            samples = [s for s in f["samples"].split(",") if s and s != "NA"]
            toks = [t.split("_") for t in f["var"].split(",")]
            named = [(int(t[-3]), t[-2], t[-1]) for t in toks]
            for n in named:
                if not truth.is_record(*n):
                    fail("%s: %s is not a record of the VCF" % (where, "_".join(map(str, n))))
            for s in samples:
                if truth.phased(s):
                    if not truth.carries_at(s, stripped.upper(), f["strand"], f["start"],
                                            first=named):
                        fail("%s: %s is listed but no haplotype of it has the sequence"
                             % (where, s))
                    elif not truth.holds_all(s, stripped.upper(), f["strand"], f["start"], named):
                        fail("%s: no haplotype of %s has the sequence spanning the indel with "
                             "every variant the row names" % (where, s))
                else:
                    for t in toks:
                        if not truth.carries_allele(s, int(t[-3]), t[-2], t[-1]):
                            fail("%s: %s is listed (unphased) without %s"
                                 % (where, s, "_".join(t)))
            site = truth.site_of(f["start"])
            if site:
                seen.setdefault(site, []).append((set(samples), f["cfd"], f["mm"], f["var"]))
        for site, ex in sorted(expect.items()):
            got = seen.get(site, [])
            listed = set().union(*[s for s, _c, _m, _v in got]) if got else set()
            if ex.get("no_indel_rows") and got:
                fail("%s: site %s: indel row(s) on a target that holds no base of the indel: %s"
                     % (name, site, got))
            if "snp_perfect" in ex:
                lo, hi = truth.spec["sites"][site]
                best = [set(f["samples"].split(",")) for f in snp_only
                        if lo - 60 <= f["start"] <= hi + 60 and f["block"] == "highest_CFD"
                        and f["cfd"] is not None and f["cfd"] >= 0.999]
                if not any(set(ex["snp_perfect"]) <= b for b in best):
                    fail("%s: site %s: no CFD-1 SNP row in the integrated results lists %s"
                         % (name, site, ex["snp_perfect"]))
            for s in ex.get("never_perfect", ()):
                if any(s in ss and c is not None and c >= 0.999 for ss, c, _m, _v in got):
                    fail("%s: site %s: %s is listed on a CFD-1 row it does not have"
                         % (name, site, s))
            for s, want in sorted(ex.get("named", {}).items()) if per_sample else ():
                if not any(s in ss and all(w in v.split(",") for w in want)
                           for ss, _c, _m, v in got):
                    fail("%s: site %s: no row lists %s under %s" % (name, site, s, want))
            for s in ex.get("bare", ()) if per_sample else ():
                if not any(s in ss and len(v.split(",")) == 1 for ss, _c, _m, v in got):
                    fail("%s: site %s: no row lists %s under the indel alone" % (name, site, s))
            for s in ex.get("perfect", ()):
                if not any(s in ss and c is not None and c >= 0.999 for ss, c, _m, _v in got):
                    fail("%s: site %s: no CFD-1 row lists %s" % (name, site, s))
            for s in ex.get("one_mm", ()) if per_sample else ():
                if not any(s in ss and m == 1 for ss, _c, m, _v in got):
                    fail("%s: site %s: no 1-mismatch row lists %s" % (name, site, s))
            for s in ex.get("listed", ()):
                if s not in listed:
                    fail("%s: site %s: %s is listed on no row" % (name, site, s))
            for s in ex.get("never", ()):
                if any(s in ss and c is not None and c >= 0.5 for ss, c, _m, _v in got):
                    fail("%s: site %s: %s is listed for a strong target it does not have"
                         % (name, site, s))
            for s in ex.get("never_any", ()):
                if s in listed:
                    fail("%s: site %s: %s is listed although no indel window holds its "
                         "haplotype there" % (name, site, s))
            if "none_above" in ex:
                hi = [c for _ss, c, _m, _v in got if c is not None and c > ex["none_above"]]
                if hi:
                    fail("%s: site %s: indel row(s) scored %s, no sample has such a target"
                         % (name, site, hi))
        # the co-occurrence companion
        cooc = []
        for path in glob.glob(os.path.join(rd, "*indel_snp_cooc.tsv")):
            with open(path) as fh:
                next(fh, None)
                cooc.extend(l.rstrip("\n").split("\t") for l in fh if l.strip())
        for c in cooc:
            if c[8] != "CONFIRMED":
                continue
            ipos = int(c[1])
            snps = [int(p) for p in c[6].split(";")]
            for s in c[11].split(","):
                ok = False
                for h in (0, 1):
                    hv = [v for v, toks, _p in truth.gts.get(s, ()) if toks[h] == "1"]
                    pos = {v["pos1"] for v in hv}
                    if ipos in pos and all(p in pos for p in snps):
                        ok = True
                if not ok:
                    fail("%s: co-occurrence %s lists %s, who does not carry the indel and "
                         "SNP %s on one haplotype" % (name, c[1], s, c[6]))
        for site, ex in sorted(expect.items()):
            lo, hi = truth.spec["sites"][site]
            for s in ex.get("cooc", ()):
                if not any(c[8] == "CONFIRMED" and lo <= int(c[1]) <= hi and s in c[11].split(",")
                           for c in cooc):
                    fail("%s: site %s: no CONFIRMED co-occurrence lists %s" % (name, site, s))
            if "cooc_copies" in ex:
                want = ex["cooc_copies"] / (2.0 * NSAMPLES)
                for c in cooc:
                    if lo <= int(c[1]) <= hi and abs(float(c[9]) - want) > 1e-5:
                        fail("%s: site %s: joint frequency %s, not %d copies of %d"
                             % (name, site, c[9], ex["cooc_copies"], 2 * NSAMPLES))
            for s in ex.get("cooc_putative", ()):
                if not any(c[8] == "PUTATIVE" and lo <= int(c[1]) <= hi and s in c[11].split(",")
                           for c in cooc):
                    fail("%s: site %s: no PUTATIVE co-occurrence lists %s" % (name, site, s))
    print("oracle: %d failure(s)" % len(fails))
    return 1 if fails else 0


def main(argv):
    if len(argv) >= 3 and argv[1] == "generate":
        generate(argv[2])
        return 0
    if len(argv) >= 5 and argv[1] == "check":
        return check(argv[2], argv[3], argv[4:])
    sys.stdout.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
