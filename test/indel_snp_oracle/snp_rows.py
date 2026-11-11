"""The SNP+indel oracle's SNP-row part (2.10.0): sites Y1-Y4, and a per-sample check of
every SNP-only row's carriers against their haplotypes. oracle.py calls ``plant`` from
``generate`` and ``check`` from its ``check``.

A SNP row's target is the reference with the SNP alleles of a haplotype: the enriched
genome holds SNPs only, so up to 2.8.1 a SNP row's carriers came from SNP genotypes alone
and it listed every carrier of its SNPs, also one whose haplotype holding them has an indel
inside the target, though that haplotype does not have the row's sequence (an HGDP + 1000
Genomes row: chr22:21052638 -, its SNP 21052641 G>A on the same 1,448 haplotypes as
21052640 T>TCACC). Since 2.10.0 a per-sample SNP row does not list such a carrier.

Stdlib only.
"""

import glob
import os

# the sites' place on the oracle's reference and the samples they use (S70-S78): clear
# of the other sites' and samples, so the oracle can grow on either side
Y_SITES = (17000, 17300, 17600, 17900)
LENGTH = 21000


def plant(put, var, expect, seq, guide, rc):
    """Plant sites Y1-Y4 with oracle.generate's ``put`` and ``var``: a SNP target of the
    enriched genome, a C>A SNP at its 6th base making it perfect, and an indel inside it
    on some haplotypes."""
    g = guide
    y1, y2, y3, y4 = Y_SITES
    # Y1 (+): a 2-bp insertion after the 13th base, in cis (S70), in trans (S71),
    # homozygous (S72); S73 has the SNP alone; S74 both unphased (its slot unknown: it
    # stays, PUTATIVE).
    put(y1, g[:5] + "C" + g[6:] + "TGG")
    var("Y1", "snp05", y1 + 6, "C", "A",
        {"S70": "1|0", "S71": "1|0", "S72": "1|0", "S73": "1|0", "S74": "0/1"})
    var("Y1", "ins", y1 + 13, seq[y1 + 12], seq[y1 + 12] + "TT",
        {"S70": "1|0", "S71": "0|1", "S72": "1|1", "S74": "0/1"})
    expect["Y1"] = {"snp_listed": ["S71", "S73", "S74"], "snp_never": ["S70", "S72"]}
    # Y2 (+): a phased SNP carrier whose insertion call is unphased (S75): its slot cannot
    # be tied to the SNP's, so it stays listed, PUTATIVE, no longer CONFIRMED.
    put(y2, g[:5] + "C" + g[6:] + "TGG")
    var("Y2", "snp05", y2 + 6, "C", "A", {"S75": "1|0"})
    var("Y2", "ins", y2 + 13, seq[y2 + 12], seq[y2 + 12] + "TT", {"S75": "0/1"})
    expect["Y2"] = {"snp_listed": ["S75"], "snp_phase": "PUTATIVE"}
    # Y3 (+): the PAM TGG followed by one more G, and on S76's SNP slot a 1-bp deletion in
    # that G run: the run puts a G back, the target's bases are unchanged, S76 stays.
    put(y3, g[:5] + "C" + g[6:] + "TGGG" + "A")
    var("Y3", "snp05", y3 + 6, "C", "A", {"S76": "1|0"})
    var("Y3", "pamrun", y3 + 21, "TG", "T", {"S76": "1|0"})
    expect["Y3"] = {"snp_listed": ["S76"]}
    # Y4 (-): the SNP (forward G>T) always in cis with a 2-bp insertion inside the target
    # (S77, S78): no haplotype has the SNP target, so it has no row at all.
    fwd = rc(g[:5] + "C" + g[6:] + "TGG")
    put(y4, fwd)
    var("Y4", "snp05", y4 + 18, "G", "T", {"S77": "1|0", "S78": "0|1"})
    var("Y4", "ins", y4 + 9, fwd[8], fwd[8] + "AA", {"S77": "1|0", "S78": "0|1"})
    expect["Y4"] = {"snp_never": ["S77", "S78"]}


def snp_row_at(truth, sample, seq, strand, start, named, rc):
    """True iff a haplotype of the sample holding every SNP a SNP-only row names has
    ``seq`` (strand-oriented) with its first base (the leftmost, forward strand) at the
    reference coordinate ``start``: that haplotype has the row's target there."""
    probe = seq if strand == "+" else rc(seq)
    for h in (0, 1):
        if not all(truth.on_slot(sample, h, *n) for n in named):
            continue
        hap, pos = truth.haplotype(sample, h)
        k = hap.find(probe)
        while k >= 0:
            if pos[k] == start:
                return True
            k = hap.find(probe, k + 1)
    return False


def check(name, rd, tables, truth, expect, fail, alt_rows, is_indel_token, rc):
    """In a per-sample results dir ``rd`` (its tables, the integrated results and the
    alternative alignments), every SNP-only row (no indel or MNV named) must list only
    carriers that have its target: a phased one on a haplotype holding every SNP it names,
    the target's first base at the row's start; an unphased one must carry every SNP it
    names. A row naming two alts of one position, which no haplotype holds, fails: before
    2.10.0 the SNP path read site N2's lone 2|1 carrier of an alt on a VCF line of several
    ALTs as carrying the other alt too and listed it on such rows; since 2.10.0 every reader
    reads a genotype against its own alt's allele index where the genotypes tell it (else the
    carrier is read without phase), except where two lines of several ALTs share a position:
    an alt of one can then be tied to an index of the other's and read on the wrong slot
    (METHODS section 4; the real-data check refuses a run that holds such a position). Then
    the Y sites' expectations (snp_listed: on a CFD-1 SNP row; snp_never: on no SNP row;
    snp_phase: the phase confirmation of the site's SNP rows)."""
    seen = {}  # site -> [(samples set, cfd)]
    for t in tables:
        for f in alt_rows(t):
            if any(is_indel_token(x) for x in f["var"].split(",")):
                continue
            where = "%s %s:%d(%s) %s %s" % (name, f["block"], f["start"], f["strand"],
                                            f["alt"], f["var"])
            stripped = f["alt"].replace("-", "").upper()
            named = [(int(t2.split("_")[-3]), t2.split("_")[-2], t2.split("_")[-1])
                     for t2 in f["var"].split(",")]
            samples = [x for x in f["samples"].split(",") if x and x != "NA"]
            if len({n[0] for n in named}) < len(named):
                fail("%s: a SNP row names two alts of one position, which no haplotype "
                     "holds (listing %s)" % (where, ",".join(samples)))
                continue
            for smp in samples:
                if truth.phased(smp):
                    if not snp_row_at(truth, smp, stripped, f["strand"], f["start"], named,
                                      rc):
                        fail("%s: %s is listed on a SNP row but no haplotype of it holding "
                             "its SNPs has the sequence at its start" % (where, smp))
                elif not all(truth.carries_allele(smp, *n) for n in named):
                    fail("%s: %s is listed (unphased) without every SNP the row names"
                         % (where, smp))
            site = truth.site_of(f["start"])
            if site:
                seen.setdefault(site, []).append((set(samples), f["cfd"]))
    for site, ex in sorted(expect.items()):
        for s in ex.get("snp_listed", ()):
            if not any(s in ss and c is not None and c >= 0.999 for ss, c in seen.get(site, [])):
                fail("%s: site %s: no CFD-1 SNP row lists %s" % (name, site, s))
        for s in ex.get("snp_never", ()):
            if any(s in ss for ss, _c in seen.get(site, [])):
                fail("%s: site %s: %s is listed on a SNP row although its haplotype holding "
                     "the SNP has an indel inside the target" % (name, site, s))
        if "snp_phase" in ex:
            lo, hi = truth.spec["sites"][site]
            got = set()
            for path in glob.glob(os.path.join(rd, "*phase_confirmation.tsv")):
                with open(path) as fh:
                    for line in fh:
                        c = line.rstrip("\n").split("\t")
                        if line.startswith("#") or len(c) < 7:
                            continue
                        if lo - 60 <= int(c[1]) <= hi + 60:
                            got.add(c[6])
            if got != {ex["snp_phase"]}:
                fail("%s: site %s: the SNP rows' phase confirmation is %s, not %s"
                     % (name, site, sorted(got), ex["snp_phase"]))
