#!/usr/bin/env python
"""The carrier gate: which variant combinations of a search's rows real haplotypes carry.

A population-level search reports, per variant window, the worst possible off-target: its
alignment takes, column by column, the allele that matches the guide best, and its row lists
every sample that carries any of those alleles. Two or more alleles taken that way can sit on
no one haplotype: no individual has the sequence the row reports. A search of a sites-only
dataset (the mega panel) has no genotypes at all, and its rows are built the same way.

The gate reads each row's variants (Variant_info_genome, for each of the three alignments of
the row) in a phased genotype store -- a genotyped dataset's own (Dictionaries/genotypes_<vcf>/
for the SNPs, indel_genotypes_<vcf>/ for the indels), the phased HGDP + 1000 Genomes store for
a sites-only dataset (the genotype oracle, opt-in) -- and gives the alignment a tier:

  A   two or more variants that one phased haplotype of the panel carries together, or that
      a haplotype reaches the alignment's CFD without: a variant the alignment takes which
      leaves the CFD the results report, the alignment's mismatches and its PAM as they are
      (an allele at a PAM base the search leaves free, a mismatch for a mismatch the table
      scores alike) is not needed for the row's score, so a haplotype carrying the rest has
      the row's off-target at the row's score. The CFD says this of the alignments it
      scores (the highest-CFD one, which the report ranks by, and the fewest-mismatches
      one), not of the one ranked by CRISPR-Bulge score, which the gate does not read.
      States: "one haplotype" when a single haplotype carries them all, "phase not
      known" when only samples whose genotypes do not say which haplotype holds which
      variant carry them all, "same score" (with "same score, phase not known" and "same
      score as the reference") when a haplotype reaches the alignment's CFD without every
      variant;
  B   one variant: the panel or the dataset's allele frequency says that some individual
      carries it. It cannot be a combination nobody carries, which is all B says;
  C1  two or more variants, each one in the panel, that no haplotype of it carries together.
      E, the number of panel haplotypes expected to carry them together if they were
      independent, says how much the zero means: below the threshold (3 by default) the
      state is "too rare to tell". Or variants (one or more) every haplotype carrying which
      also holds an indel inside the target, so that none of them has the row's sequence
      (state "indel in the target");
  C2  two or more variants, one or more not in the panel: it cannot say. The safeguards'
      frequency bound (the product of the variants' frequencies >= 1e-6 and the rarest
      >= 0.001) is given with them. In a search of a sites-only dataset (the mega: no
      genotypes, so its rows are the worst case of each window, the immunoglobulin loci's
      combinations of 10 to 17 variants among them) a C2 alignment below the bound is
      ranked after the others, as a C1 one is (state "below the frequency bound": plan
      decision D6, demote and label, never drop); elsewhere the bound is information only.
      With no phased panel to read (a sites-only dataset searched without the genotype
      oracle) every alignment of two or more variants is C2, checked against the bound
      alone.

A haplotype holds a row only when it holds the row's target: on a phased slot that carries
the row's SNPs, an indel of the store that changes the bases the target covers means the slot
does not have the row's sequence, as a per-sample SNP row reads its carriers (METHODS section
4) -- one whose REF holds a listed SNP's base at another position too (the slot holds two
records of that base, so its bases cannot be read). One at a listed SNP's own position is
another allele of that base (a phasing artefact: it takes nothing off, and the slot holds two
alleles of one position, "phase not known"), and one that leaves the target's bases where
they were, read on the reference from the target's first base on or up to its last (a
tandem repeat's unit gained or lost), changes nothing. An unphased call keeps the carrier
("phase not known"), as a slot that holds two ALTs at one position does (the PUTATIVE rule
for two ALTs on one slot).

The gate changes no row and drops none: it writes a companion TSV (<results>_carrier_gate.tsv,
one line per alignment that uses a variant) and a sidecar the report reads
(.carrier_gate.json). The report ranks the C1 rows, and the C2 rows below the frequency bound
of a sites-only dataset, after the others, in a table of their own, with the factual label of
each tier. In a per-sample search the gate reads only the rows that name no carrier (a
sites-only dataset's, which are the population-level analysis's in every search): the
others name their cis carriers.

STDLIB ONLY (with the repository's tier readers).
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import math
import os
import sys
from collections import namedtuple

import build_mega_indel_af as bia
import bulge_cfd_penalty as bcp
import indel_snp_cis as isc
import pam_filter as pf
import population_summary as ps
import tier0_compile as t0c
import tier0_registry as t0
import tier1_genotypes as t1

TIER_A, TIER_B, TIER_C1, TIER_C2 = "A", "B", "C1", "C2"
ONE_HAPLOTYPE = "one haplotype"
PHASE_NOT_KNOWN = "phase not known"
SAME_SCORE = "same score"
SAME_SCORE_PHASE = "same score, phase not known"
SAME_SCORE_REF = "same score as the reference"
SAME_SCORE_STATES = (SAME_SCORE, SAME_SCORE_PHASE, SAME_SCORE_REF)
# the counts of the sidecar hold the tiers, and this key the A rows of a "same score" state
SAME_SCORE_COUNT = "A_same_score"
TOO_RARE = "too rare to tell"
INDEL_IN_TARGET = "indel in the target"
NOT_IN_PANEL = "not in the panel"
# D6: a C2 alignment of a sites-only search whose frequencies fail the bound
BELOW_BOUND = "below the frequency bound"
# the counts of the sidecar hold the tiers, and this key the C2 rows below the bound
BELOW_BOUND_COUNT = "C2_below_bound"

# D2: below this many expected haplotypes a zero says little ("too rare to tell"). A
# parameter (complete-search --gate-min-expected); 0 gives every C1 row the same state.
MIN_EXPECTED = 3.0
# The safeguards' frequency bound (plan section 2.3, decision D6): the product of the
# variants' frequencies and the rarest of them. A C2 alignment of a sites-only search below
# it is ranked after the others (BELOW_BOUND); elsewhere it is shown only.
FREQ_BOUND_PRODUCT = 1e-6
FREQ_BOUND_RAREST = 0.001


def below_bound(af_product, af_min):
    """True when frequencies fail the bound: the product below FREQ_BOUND_PRODUCT or the
    rarest below FREQ_BOUND_RAREST; False when they are not known (None)."""
    if af_product is None or af_min is None:
        return False
    return af_product < FREQ_BOUND_PRODUCT or af_min < FREQ_BOUND_RAREST

# the oracle a sites-only search is checked against (download --what genotype-oracle)
ORACLE_PANEL = "hg38_HGDP1kGP"

# the integrated results' three alignments of a row, by the suffix of their columns
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
# the alignment the report's second ranking takes, which has no CFD column: the gate reads
# no CRISPR-Bulge score (the model is not run again), so the "same score" states, which the
# CFD decides, are not its to give
CB_BLOCK = BLOCKS[2]
_BLOCK_COLUMNS = ("Start_coordinate", "Strand", "Aligned_spacer+PAM",
                  "Aligned_protospacer+PAM_ALT", "PAM", "REF/ALT_origin",
                  "Variant_info_genome", "Variant_MAF")

# the alignment CFD is defined for: a 20-nt spacer and a 3-base 3' PAM (cfd_score)
SPACER_PAM_LEN = 23
CFD_PAM_LEN = 3
# the integrated results report a CFD of three decimals: two alignments have the CFD the row
# reports when they round to the same three decimals, which is what the label says of them
CFD_DECIMALS = 3
# a minus-strand alignment is written on the other strand, base by base
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

COMPANION_SUFFIX = "_carrier_gate.tsv"
META_NAME = ".carrier_gate.json"
COMPANION_COLUMNS = (
    "Row", "Block", "Spacer+PAM", "Chromosome", "Start", "Strand", "Variants", "Tier",
    "State", "Panel", "Panel_haplotypes", "Haplotypes", "Samples", "Expected",
    "AF_product", "AF_min", "Label", "Scoring_variants",
)

# the search's PAM (pam_filter: the motif the PAM file asks for, and whether it is 3'), read
# by the gate over an alignment's PAM region
Pam = namedtuple("Pam", "motif three_prime")

# one alignment of a row of the integrated results, as the gate reads it ("pam": the
# search's PAM, None where the gate was not given it)
Alignment = namedtuple(
    "Alignment", "row block guide start strand crrna aligned pam_len variants maf cfd pam",
    defaults=(None,))


# --------------------------------------------------------------------------- #
# Variants of a row
# --------------------------------------------------------------------------- #
def parse_variants(field):
    """The variants of a Variant_info_genome cell, [(chrom, pos, REF, ALT)] in the order it
    lists them, each once. A token is chrom_pos_REF_ALT (a chromosome name may hold '_'); a
    reference alignment's cell ("NA", "n", empty) has none."""
    out = []
    for tok in str(field).split(","):
        tok = tok.strip()
        if not tok or tok in ("NA", "n", ".", "-"):
            continue
        parts = tok.rsplit("_", 3)
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts
        try:
            v = (chrom, int(pos), ref.upper(), alt.upper())
        except ValueError:
            continue
        if v not in out:
            out.append(v)
    return out


def is_snp(v):
    return len(v[2]) == 1 and len(v[3]) == 1


def target_span(start, aligned):
    """(lo, hi): the reference bases (0-based, inclusive) an alignment of a SNP-only row
    covers. Start_coordinate is the 0-based leftmost base on either strand, and every
    character of the aligned protospacer + PAM but an RNA-bulge gap ('-') is a genome base."""
    n = sum(1 for ch in str(aligned) if ch != "-")
    if n <= 0:
        return None
    return int(start), int(start) + n - 1


# --------------------------------------------------------------------------- #
# Which variants of an alignment its score needs
# --------------------------------------------------------------------------- #
def alignment_cfd(crrna, aligned, pam_len=CFD_PAM_LEN):
    """The CFD of one alignment, from its aligned spacer+PAM and protospacer+PAM
    (bulge_cfd_penalty, as the search scores a row), or None where CFD is not defined: a
    spacer or a PAM of another length, or two strings that do not align."""
    if (pam_len != CFD_PAM_LEN or len(crrna) != len(aligned)
            or len(crrna) - crrna.count("-") != SPACER_PAM_LEN):
        return None
    return bcp.cfd_score_bulge_aware(crrna.upper(), aligned.upper(), pam_len=pam_len)


def variant_columns(start, aligned, strand, variants):
    """[(column, reference base)] of each variant of a SNP-only alignment in its aligned
    protospacer+PAM, read as target_span reads it (every character but an RNA-bulge gap is
    a genome base, the leftmost first; a minus-strand alignment is the reverse complement).
    None when a variant falls outside the alignment or the alignment does not take its
    allele there."""
    cols = [i for i, ch in enumerate(aligned) if ch != "-"]
    out = []
    for _chrom, pos, ref, alt in variants:
        off = (pos - 1) - int(start)
        if not (0 <= off < len(cols)):
            return None
        i = cols[len(cols) - 1 - off] if strand == "-" else cols[off]
        if strand == "-":
            ref, alt = COMPLEMENT.get(ref, ref), COMPLEMENT.get(alt, alt)
        if aligned[i] != alt:
            return None
        out.append((i, ref))
    return out


def same_cfd(one, other):
    """Whether two alignments have the CFD the integrated results report of them: the same
    three decimals (CFD_DECIMALS), which is what a "same score" label says of them."""
    return (one is not None and other is not None
            and round(one, CFD_DECIMALS) == round(other, CFD_DECIMALS))


def matches(guide_base, base):
    """Whether a column of an alignment matches the guide: its base, or the N of a PAM
    the guide's side writes as N, which any base matches."""
    return guide_base == "N" or guide_base == base


def read_pam(path):
    """The search's PAM (Pam) from its PAM file, as the post-analysis reads it
    (pam_filter._requested_pam: the motif and the side it sits on), or None."""
    try:
        motif, three_prime = pf._requested_pam(path)
    except (OSError, ValueError, IndexError):
        return None
    return Pam(motif.upper(), three_prime) if motif else None


def pam_asks(aligned, pam_len, pam):
    """{column: the IUPAC character the search's PAM asks there} over an alignment's PAM
    region -- its last ``pam_len`` characters, or its first for a 5' PAM, as the
    post-analysis reads it. The guide's side of an alignment writes its PAM as N, which
    every base matches, so a column of the region is read against the search's PAM
    instead: an allele the PAM needs is one the row needs, however the CFD scores it.
    Without the search's PAM (``pam`` None, or a motif of another length) the gate cannot
    say which bases make a PAM, and asks for "" -- which no base gives -- at either end's
    ``pam_len`` columns."""
    n = len(aligned)
    if pam_len <= 0 or pam_len > n:
        return {}
    if pam is None or len(pam.motif) != pam_len:
        return {c: "" for c in list(range(pam_len)) + list(range(n - pam_len, n))}
    cols = range(n - pam_len, n) if pam.three_prime else range(pam_len)
    return {c: pam.motif[i] for i, c in enumerate(cols)}


def pam_gives(base, asks):
    """Whether ``base`` makes a PAM at a column that asks for the IUPAC character ``asks``
    (pam_filter, as the post-analysis gates a target's PAM); "" is given by no base."""
    return bool(asks) and pf._pam_can_match(base, asks)


def scoring_variants(crrna, aligned, placed, pam_len=CFD_PAM_LEN, shown=None, pam=None):
    """The indexes of the variants of an alignment its score needs: those whose allele,
    taken back to the reference base (``placed``: variant_columns), change the CFD the row
    reports, or turn a column that matches the guide into one that does not, or the other
    way round (the row's mismatches are the alignment's, and one more would be another
    alignment, outside the search's budget), or leave the PAM region a sequence the
    search's PAM (``pam``) does not take, where there is no off-target at all. The others
    leave the row's score and its mismatches as they are -- an allele at a PAM base the
    search leaves free, a mismatch for a mismatch the table scores alike -- so a haplotype
    carrying the rest alone has the row's off-target at the row's score. Every index when
    they change the CFD together (the PAM's factor is read over its bases at once), and
    None when the gate cannot read the alignment's CFD, or ``shown`` (the row's CFD) is not
    the CFD it reads."""
    whole = alignment_cfd(crrna, aligned, pam_len)
    if not whole:  # another geometry, or a score of 0 that no allele can lower
        return None
    if shown is not None and not same_cfd(whole, shown):
        return None
    asks = pam_asks(aligned, pam_len, pam)
    keep = [k for k, (col, ref) in enumerate(placed)
            if not same_cfd(alignment_cfd(crrna, aligned[:col] + ref + aligned[col + 1:],
                                          pam_len), whole)
            or matches(crrna[col], ref) != matches(crrna[col], aligned[col])
            or (col in asks and not pam_gives(ref, asks[col]))]
    if len(keep) == len(placed):
        return keep
    rest = list(aligned)
    for k, (col, ref) in enumerate(placed):
        if k not in keep:
            rest[col] = ref
    if not same_cfd(alignment_cfd(crrna, "".join(rest), pam_len), whole):
        return list(range(len(placed)))
    return keep


def score_breakers(panel, alignment, variants):
    """{sample: genotype} of each allele of the panel inside the alignment, other than the
    ones the row takes, that changes the CFD the row reports or the alignment's
    mismatches, or leaves the PAM region a sequence the search's PAM does not take: a slot
    holding one does not have the row's off-target at the row's score, however the row's
    own variants sit on it."""
    aligned = alignment.aligned.upper()
    whole = alignment_cfd(alignment.crrna, aligned, alignment.pam_len)
    if whole is None:
        return []
    minus = alignment.strand == "-"
    cols = [i for i, ch in enumerate(aligned) if ch != "-"]
    own = {(v[1], v[3]) for v in variants}
    asks = pam_asks(aligned, alignment.pam_len, alignment.pam)
    out = []
    for off, col in enumerate(cols):
        pos = int(alignment.start) + (len(cols) - 1 - off if minus else off) + 1
        for alt, gts in panel.alleles_at(pos).items():
            if (pos, alt) in own or len(alt) != 1:
                continue
            base = COMPLEMENT.get(alt, alt) if minus else alt
            if base == aligned[col]:
                continue
            if (not same_cfd(alignment_cfd(alignment.crrna,
                                           aligned[:col] + base + aligned[col + 1:],
                                           alignment.pam_len), whole)
                    or matches(alignment.crrna[col], base)
                    != matches(alignment.crrna[col], aligned[col])
                    or (col in asks and not pam_gives(base, asks[col]))):
                out.append(gts)
    return out


# --------------------------------------------------------------------------- #
# Slots: the haplotypes of one sample that hold a set of alleles
# --------------------------------------------------------------------------- #
def _tokens(gt):
    """(allele tokens, kind) of a genotype read against one allele ("1" = this allele):
    kind "phased", "unphased" or "haploid" (indel_snp_cis._gt_tokens)."""
    return isc._gt_tokens(gt)


def _exact(parsed):
    """The number of slots when every genotype of ``parsed`` ([(tokens, kind)]) can be
    put on a slot: all phased with one ploidy, or all haploid. None otherwise."""
    kinds = {kind for _toks, kind in parsed}
    lens = {len(toks) for toks, _kind in parsed}
    if len(lens) != 1 or not (kinds == {"phased"} or kinds == {"haploid"}):
        return None
    return lens.pop()


def support(variant_gts, breakers=(), rivals=()):
    """Which haplotypes of a panel hold a set of alleles together.

    Args:
      variant_gts: one {sample: genotype} per allele of the set, each genotype read against
        its own allele ("1|0": the first slot holds it).
      breakers: {sample: genotype} of each variant whose allele, on a slot, takes the row's
        target off it (an indel that changes the target's bases).
      rivals: {sample: genotype} of each other allele at a position of the set: a slot that
        holds two alleles of one position is a phasing artefact.

    Returns (haplotypes, cis, unknown, carrying): the number of phased slots that hold every
    allele and no breaker, the samples with such a slot, the samples that carry every
    allele with no such slot but one the genotypes cannot place (an unphased call, of an
    allele or of a breaker, or two alleles of one position on one slot), and the number of
    phased slots that hold every allele, breakers or not. A sample that carries them all
    only on different slots (trans) is in none."""
    haps = carrying = 0
    cis, unknown = set(), set()
    if not variant_gts:
        return haps, cis, unknown, carrying
    common = set(variant_gts[0])
    for gts in variant_gts[1:]:
        common &= set(gts)
    for s in sorted(common):
        parsed = [_tokens(gts[s]) for gts in variant_gts]
        n = _exact(parsed)
        if n is None:
            unknown.add(s)
            continue
        slots = [h for h in range(n) if all(toks[h] == "1" for toks, _kind in parsed)]
        if not slots:
            continue
        carrying += len(slots)
        held, unsure = 0, False
        for h in slots:
            state = "held"
            for rv in rivals:
                g = rv.get(s)
                if g is None:
                    continue
                rt = _tokens(g)
                if _exact(parsed + [rt]) == n and rt[0][h] == "1":
                    state = "artefact"
                    break
            if state == "artefact":
                unsure = True
                continue
            for br in breakers:
                g = br.get(s)
                if g is None:
                    continue
                bt = _tokens(g)
                if _exact(parsed + [bt]) == n:
                    if bt[0][h] == "1":
                        state = "broken"
                        break
                elif "1" in bt[0]:
                    state = "unsure"
            if state == "held":
                held += 1
            elif state == "unsure":
                unsure = True
        if held:
            haps += held
            cis.add(s)
        elif unsure:
            unknown.add(s)
    return haps, cis, unknown, carrying


# --------------------------------------------------------------------------- #
# A panel's phased store on one chromosome
# --------------------------------------------------------------------------- #
class IndelIndex(object):
    """The indels of a phased indel store (build_indel_genotypes: one gzipped line per
    indel, "<pos>_<REF>_<ALT>\t<sample>:<genotype>,..."), by position, with their carriers
    read on demand. ``wanted``: the 0-based reference positions a search's rows touch; only
    an indel whose REF (or the base past it, where an insertion goes) meets one of them is
    kept, so a chromosome's whole store is not held for a few targets. None keeps all."""

    def __init__(self, path, wanted=None):
        self._toks = {}
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                key, _, toks = line.rstrip("\n").partition("\t")
                pos, _, rest = key.partition("_")
                ref, _, alt = rest.partition("_")
                try:
                    pos = int(pos)
                except ValueError:
                    continue
                if wanted is not None and not any(
                        p in wanted for p in range(pos - 1, pos + len(ref))):
                    continue
                self._toks[(pos, ref, alt)] = toks
        self.indels = sorted(self._toks)
        self.starts = [p for p, _r, _a in self.indels]
        self.maxref = max((len(r) for _p, r, _a in self.indels), default=1)
        self._gts = {}

    def carriers(self, pos, ref, alt):
        key = (pos, ref, alt)
        gts = self._gts.get(key)
        if gts is None:
            gts = {}
            for tok in self._toks.get(key, "").split(","):
                s, _, g = tok.partition(":")
                if s:
                    gts[s] = g
            self._gts[key] = gts
        return gts

    def has(self, pos, ref, alt):
        return (pos, ref, alt) in self._toks

    def near(self, lo, hi):
        """Every indel whose REF meets the reference bases lo..hi (0-based, inclusive) or
        whose POS is in it."""
        k = bisect.bisect_left(self.starts, lo - self.maxref + 1)
        while k < len(self.indels) and self.indels[k][0] - 1 <= hi:
            p, r, a = self.indels[k]
            k += 1
            if p - 1 + len(r) <= lo - 1:
                continue
            yield p, r, a

    def at(self, pos):
        k = bisect.bisect_left(self.starts, pos)
        while k < len(self.indels) and self.indels[k][0] == pos:
            yield self.indels[k]
            k += 1


def changing_indels(index, lo, hi):
    """[(POS, REF, ALT, {sample: genotype})] of each indel of ``index`` that changes the
    reference bases lo..hi (0-based, inclusive) of a haplotype carrying it
    (indel_snp_cis.changes_span): a target's, once for all its haplotypes."""
    out = []
    if index is None:
        return out
    for p, r, a in index.near(lo, hi):
        if isc.changes_span(p, r, a, lo, hi):
            out.append((p, r, a, index.carriers(p, r, a)))
    return out


def breakers_of(changing, alleles, lo, hi, seq_at=None, own=()):
    """{sample: genotype} of each of ``changing`` (changing_indels) that takes a target's
    bases lo..hi off a haplotype holding the SNP alleles ``alleles`` ({POS: base}), as a
    per-sample SNP row's carrier is read (METHODS section 4): not one of ``own``, not one at
    a SNP's position (another allele of that base: a slot phased with both is a phasing
    artefact, which support() reads among the rivals), and, given the reference
    (``seq_at``), not one that leaves the target's bases where they were (holds_target: a
    tandem repeat's unit gained or lost) -- but one whose REF holds a SNP's base at another
    position whatever the reference says: a slot holding both holds two records of that
    base, whose bases cannot be read. Without the reference, every other indel that
    changes the bases."""
    out = []
    for p, r, a, g in changing:
        if (p, r, a) in own or p in alleles:
            continue
        over = any(p < q <= p + len(r) - 1 for q in alleles)
        if not over and seq_at is not None and holds_target(seq_at, lo, hi, alleles,
                                                            (p, r, a)):
            continue
        out.append(g)
    return out


def holds_target(seq_at, lo, hi, alleles, indel):
    """Whether a haplotype that holds the SNP alleles ``alleles`` ({POS: base}) and the
    indel ``indel`` (POS, REF, ALT) still has the target's bases lo..hi (0-based, inclusive):
    the reference there with those alleles, read on the haplotype from the target's first
    base on, or up to its last base, as the indel rows read a second indel
    (indel_snp_cis.target_pairs). In a tandem repeat a unit gained or lost at the target's
    edge leaves its bases where they were. ``seq_at(a, b)``: the reference bases a..b-1
    (0-based, upper case; shorter at a chromosome's end). The bases the indel's REF and ALT
    share at the start or at the end stay reference bases."""
    p, r, a = indel
    p0 = p - 1
    n = hi - lo + 1
    target = "".join(alleles.get(q + 1, seq_at(q, q + 1)) for q in range(lo, hi + 1))
    c = isc.kept_bases(r, a)
    t = isc.shared_tail(r, a)
    gone_from, gone_to = p0 + c, p0 + len(r) - t  # reference bases the indel takes out
    put = a[c:len(a) - t].upper()                  # and the bases it puts in their place
    pad = n + len(r) + len(a)
    first = max(0, min(lo, p0) - pad)
    ref = seq_at(first, max(hi, p0 + len(r)) + 1 + pad)
    hap = []  # (reference position or None, base)
    for q in range(first, first + len(ref)):
        if q == gone_from:
            hap.extend((None, b) for b in put)
        if gone_from <= q < gone_to:
            continue
        hap.append((q, alleles.get(q + 1, ref[q - first])))
    if gone_from >= first + len(ref):
        hap.extend((None, b) for b in put)
    where = {q: i for i, (q, _b) in enumerate(hap) if q is not None}
    bases = "".join(b for _q, b in hap)
    i, j = where.get(lo), where.get(hi)
    return ((i is not None and bases[i:i + n] == target)
            or (j is not None and j - n + 1 >= 0 and bases[j - n + 1:j + 1] == target))


def breaking_indels(index, lo, hi, snps, own=(), seq_at=None):
    """{sample: genotype} of each indel of ``index`` that takes the target's bases lo..hi
    (0-based, inclusive) off a haplotype holding the SNPs ``snps`` ([(pos, REF, ALT)]):
    breakers_of."""
    return breakers_of(changing_indels(index, lo, hi), {p: a for p, _r, a in snps}, lo, hi,
                       seq_at, own)


class PanelChrom(object):
    """One panel's phased genotypes on one chromosome: the SNP store (Tier-1) and the indel
    store, either absent. ``carriers`` gives a variant's carriers with each genotype read
    against that variant ("1|0"), or None for a variant the panel does not hold.
    ``wanted``: the positions the rows touch, for the indel store (IndelIndex); ``seq_at``:
    the chromosome's reference bases (holds_target), or None."""

    def __init__(self, name, dictionaries, chrom, role="own", wanted=None, seq_at=None):
        self.name = name
        self.seq_at = seq_at
        self.role = role
        self.chrom = chrom
        self.gt = None
        self.indels = None
        self.ref_mismatches = 0
        gdir = os.path.join(dictionaries, f"genotypes_{name}")
        gbin = os.path.join(gdir, f"gt_{chrom}.bin")
        gidx = os.path.join(gdir, f"gt_{chrom}.idx")
        if os.path.isfile(gbin) and os.path.isfile(gidx):
            self.gt = t1.GenotypeReader(gbin, gidx)
        ipath = os.path.join(dictionaries, f"indel_genotypes_{name}", f"gt_indel_{chrom}.tsv.gz")
        if os.path.isfile(ipath):
            self.indels = IndelIndex(ipath, wanted)
        self._sites = {}
        self._ploidy = {}
        self._an = None
        if self.gt is not None:
            axis = self.gt._axis
            ploidy_of = t0c.ploidy_of_for_chrom(chrom)
            for i in range(len(axis)):
                _db, _sp, sex = axis.sample_meta(i)
                self._ploidy[axis.sample_id(i)] = ploidy_of(axis.sample_id(i), sex)
            if t0c._norm_chrom(chrom) == "X":
                self._an = t0c.chrx_par_panels(axis, ps.Panel)
            else:
                panel = ps.Panel(axis, ploidy_of)
                self._an = lambda _pos, _p=panel: _p
        elif self.indels is not None:
            self._ploidy = None

    def close(self):
        if self.gt is not None:
            self.gt.close()

    def haplotypes(self, pos):
        """The panel's haplotypes at ``pos`` (its AN: chrX males are haploid outside PAR)."""
        if self._an is None:
            return 0
        return self._an(pos).global_[0]

    def _snp_site(self, pos):
        site = self._sites.get(pos)
        if site is None:
            rec = {}
            refs = {}
            for alt in self.gt.alts_at(pos):
                rec[alt] = dict(self.gt.carriers(pos, alt) or [])
                refs[alt] = (self.gt.ref(pos, alt) or "").upper()
            site = (isc.biallelic_site(rec) if rec else {}, refs)
            self._sites[pos] = site
        return site

    def carriers(self, v):
        _chrom, pos, ref, alt = v
        if is_snp(v):
            if self.gt is None:
                return None
            site, refs = self._snp_site(pos)
            if alt not in site:
                return None
            if refs.get(alt) and refs[alt] != ref:
                self.ref_mismatches += 1
                return None
            return site[alt]
        if self.indels is None or not self.indels.has(pos, ref, alt):
            return None
        return self.indels.carriers(pos, ref, alt)

    def alleles_at(self, pos):
        """{ALT: {sample: genotype}} of the panel's SNP alleles at ``pos``."""
        if self.gt is None:
            return {}
        site, _refs = self._snp_site(pos)
        return site

    def rivals(self, v):
        """{sample: genotype} of each other allele at the variant's position: the other ALTs
        of a SNP site, and the indels that start there."""
        _chrom, pos, ref, alt = v
        out = []
        if self.gt is not None:
            site, _refs = self._snp_site(pos)
            out.extend(gts for a, gts in site.items() if not (is_snp(v) and a == alt))
        if self.indels is not None:
            for p, r, a in self.indels.at(pos):
                if (p, r.upper(), a.upper()) != (pos, ref, alt):
                    out.append(self.indels.carriers(p, r, a))
        return out

    def allele_freq(self, v, gts):
        """The variant's frequency in the panel: its copies over the panel's haplotypes."""
        an = self.haplotypes(v[1])
        if not an:
            return 0.0
        copies = 0
        for s, g in gts.items():
            toks = _tokens(g)[0]
            cap = self._ploidy.get(s, len(toks)) if self._ploidy else len(toks)
            copies += sum(1 for t in toks[:max(cap, 1)] if t == "1")
        return copies / an

    def samples(self):
        return len(self.gt._axis) if self.gt is not None else 0

    def superpopulations(self):
        """{superpopulation: samples} of the panel (the Tier-1 axis keeps it per sample)."""
        out = {}
        if self.gt is None:
            return out
        axis = self.gt._axis
        for i in range(len(axis)):
            _db, sp, _sex = axis.sample_meta(i)
            out[sp] = out.get(sp, 0) + 1
        return out


# --------------------------------------------------------------------------- #
# Tiers
# --------------------------------------------------------------------------- #
class Verdict(object):
    """The tier of one alignment: ``tier`` (A, B, C1, C2), ``state`` (its sub-state or
    ""), the panel that decided it and its haplotypes (``panel_haplotypes``), the haplotypes
    that hold the variants together (``haplotypes``), the samples to name (``samples``),
    ``expected`` (E), the frequency bound's ``af_product`` / ``af_min``, ``missing``,
    the number of variants no panel holds, and ``kept``, the variants the alignment's CFD
    needs where a haplotype reaches it without every one of them (scoring_variants)."""

    __slots__ = ("tier", "state", "panel", "panel_haplotypes", "haplotypes", "samples",
                 "expected", "af_product", "af_min", "missing", "n_variants", "kept")

    def __init__(self, tier, state="", panel=None, panel_haplotypes=0, haplotypes=None,
                 samples=(), expected=None, af_product=None, af_min=None, missing=0,
                 n_variants=0, kept=()):
        self.tier = tier
        self.state = state
        self.panel = panel
        self.panel_haplotypes = panel_haplotypes
        self.haplotypes = haplotypes
        self.samples = tuple(samples)
        self.expected = expected
        self.af_product = af_product
        self.af_min = af_min
        self.missing = missing
        self.n_variants = n_variants
        self.kept = tuple(kept)

    def label(self):
        """The factual text the report shows (ASCII: it goes into the TSV downloads too)."""
        hap_of = f"{(self.panel_haplotypes or 0):,} {self.panel} haplotypes"
        if self.tier == TIER_A:
            if self.state in SAME_SCORE_STATES:
                rest = self.n_variants - len(self.kept)
                others = (f"the other variant does not change it" if rest == 1
                          else f"the other {rest} variants do not change it")
                if self.state == SAME_SCORE_REF:
                    return (f"the same CFD as the reference here: none of the "
                            f"{self.n_variants} variants changes it")
                carried = " and ".join(self.kept)
                if self.state == SAME_SCORE_PHASE:
                    n = len(self.samples)
                    return (f"the same CFD in {n} {self.panel} sample{'' if n == 1 else 's'} "
                            f"carrying {carried}, phase not known: {others}")
                return (f"the same CFD on {self.haplotypes:,} of {hap_of}, which carry "
                        f"{carried}: {others}")
            if self.state == PHASE_NOT_KNOWN:
                n = len(self.samples)
                return (f"carried together by {n} {self.panel} sample{'' if n == 1 else 's'}"
                        f", phase not known")
            return f"together on {self.haplotypes:,} of {hap_of}"
        if self.tier == TIER_B:
            if self.state == NOT_IN_PANEL:
                return "single variant (not in the panel)"
            return "single variant"
        if self.tier == TIER_C1:
            if self.state == INDEL_IN_TARGET:
                n = self.haplotypes or 0
                return (f"{'the' if n == 1 else 'all'} {n:,} {self.panel} haplotype"
                        f"{'' if n == 1 else 's'} carrying "
                        f"{'it' if self.n_variants == 1 else 'them together'} "
                        f"{'has' if n == 1 else 'have'} an indel in the target")
            text = f"not together on any of {hap_of}"
            if self.state == TOO_RARE:
                return text + f"; too rare to tell (expected {self.expected:.2g})"
            return text + f" (expected {self.expected:.3g})"
        if self.panel_haplotypes is None:
            # no phased panel read at all (a sites-only search without the oracle)
            text = "not checked against phased haplotypes (no panel)"
        else:
            text = (f"not checked: {self.missing} of {self.n_variants} variants not in "
                    f"{self.panel or 'the panel'}")
        if self.af_product is not None:
            inside = not below_bound(self.af_product, self.af_min)
            text += (f"; frequency product {self.af_product:.1e}, rarest {self.af_min:.1e}"
                     f" ({'within' if inside else 'below'} the 1e-06 / 0.001 bound"
                     + (", ranked after the others" if self.state == BELOW_BOUND else "")
                     + ")")
        return text


def panel_support(variants, panels, span=None, quick_single=False, also=None):
    """([(panel, genotypes, haplotypes, cis, unknown, carrying)] of each panel that holds
    every variant (support), the fewest variants a panel misses, and -- with
    ``quick_single`` -- the first panel holding a single variant with nothing in the target
    to take a slot off it, whose slots need not be read). ``also(panel)``: more alleles
    that take a slot off the set (score_breakers)."""
    full, missing, quick = [], len(variants), None
    snps_only = all(is_snp(v) for v in variants)
    for panel in panels:
        gts = [panel.carriers(v) for v in variants]
        absent = sum(1 for g in gts if not g)
        missing = min(missing, absent)
        if absent:
            continue
        breakers = []
        if snps_only and span is not None:
            breakers = breaking_indels(panel.indels, span[0], span[1],
                                       [(v[1], v[2], v[3]) for v in variants],
                                       seq_at=panel.seq_at)
        if also is not None:
            breakers = breakers + also(panel)
        if quick_single and len(variants) == 1 and not breakers:
            quick = panel
            break
        rivals = []
        for v in variants:
            rivals.extend(panel.rivals(v))
        haps, cis, unknown, carrying = support(gts, breakers, rivals)
        full.append((panel, gts, haps, cis, unknown, carrying))
    return full, missing, quick


def _token(v):
    return "{0}_{1}_{2}_{3}".format(*v)


def same_score_verdict(variants, panels, span, alignment):
    """The verdict of an alignment whose CFD a real haplotype reaches without every
    variant it takes (tier A, states "same score", "same score, phase not known", "same
    score as the reference"), or None: no such haplotype, every variant changes the CFD,
    the alignment is the one ranked by CRISPR-Bulge score (CB_BLOCK: the CFD does not say
    what that score is, and the gate does not run the model again, so an alignment no
    haplotype carries is C1 there), or the gate cannot read it (``alignment``: an
    Alignment, or None)."""
    if (alignment is None or alignment.block == CB_BLOCK or len(variants) < 2
            or not all(is_snp(v) for v in variants)):
        return None
    placed = variant_columns(alignment.start, alignment.aligned.upper(), alignment.strand,
                             variants)
    if placed is None:
        return None
    keep = scoring_variants(alignment.crrna, alignment.aligned.upper(), placed,
                            alignment.pam_len, alignment.cfd, alignment.pam)
    if keep is None or len(keep) == len(variants):
        return None
    pos = variants[0][1]
    sub = [variants[k] for k in keep]
    if not sub:
        # the reference sequence scores the same here, so every haplotype has the row's
        # off-target at the row's score
        panel = max(panels, key=lambda p: p.haplotypes(pos))
        return Verdict(TIER_A, SAME_SCORE_REF, panel.name, panel.haplotypes(pos),
                       n_variants=len(variants))
    full, _missing, _quick = panel_support(
        sub, panels, span, also=lambda p: score_breakers(p, alignment, variants))
    if not full:
        return None
    panel, _gts, haps, cis, unknown, _carrying = max(full, key=lambda f: (f[2], len(f[4])))
    kept = [_token(v) for v in sub]
    if haps:
        return Verdict(TIER_A, SAME_SCORE, panel.name, panel.haplotypes(pos), haps,
                       sorted(cis), n_variants=len(variants), kept=kept)
    if unknown:
        return Verdict(TIER_A, SAME_SCORE_PHASE, panel.name, panel.haplotypes(pos), 0,
                       sorted(unknown), n_variants=len(variants), kept=kept)
    return None


def classify(variants, panels, span=None, min_expected=MIN_EXPECTED, af_of=None,
             alignment=None):
    """The tier of one alignment's variants against the panels (PanelChrom).

    ``span`` (lo, hi): the target's reference bases, for an alignment of SNPs only: a slot
    that holds them with an indel changing those bases does not hold the target. ``af_of``:
    the dataset's frequency of a variant (its registry; an indel's, its indel frequencies),
    or None, for the frequency bound of a multi-variant alignment. ``alignment``: the row's
    Alignment, for a combination no haplotype carries: where a haplotype reaches the
    alignment's CFD without every variant the tier is A, not C1 (same_score_verdict), as
    the row's score is one a real haplotype has. A variant is in a panel when the panel
    holds its POS, REF and ALT; the panels are the stores a search reads (a row does not
    say which dataset it comes from): the first with evidence decides, the one with the
    most haplotypes holding the variants together for A, the most expected for C1."""
    k = len(variants)
    if k == 0:
        return None
    afs = None
    if k >= 2 and af_of is not None:
        afs = [af_of(v) for v in variants]
        if any(a is None for a in afs):
            afs = None
    af_product = math.prod(afs) if afs else None
    af_min = min(afs) if afs else None
    pos = variants[0][1]
    full, missing, quick = panel_support(variants, panels, span, quick_single=True)
    if quick is not None:
        # one variant, and nothing in the target to take a slot off it
        return Verdict(TIER_B, panel=quick.name, panel_haplotypes=quick.haplotypes(pos),
                       haplotypes=None, n_variants=1)
    if not full:
        if k == 1:
            return Verdict(TIER_B, NOT_IN_PANEL, n_variants=1)
        named = panels[0].name if len(panels) == 1 else None
        return Verdict(TIER_C2, panel=named, missing=missing, af_product=af_product,
                       af_min=af_min, n_variants=k)
    panel, gts, haps, cis, unknown, _carrying = max(full, key=lambda f: (f[2], len(f[4])))
    hap_total = panel.haplotypes(pos)
    if haps or unknown:
        if k == 1:
            return Verdict(TIER_B, panel=panel.name, panel_haplotypes=hap_total,
                           haplotypes=haps, n_variants=1)
        if haps:
            expected = hap_total * math.prod(
                panel.allele_freq(v, g) for v, g in zip(variants, gts))
            return Verdict(TIER_A, ONE_HAPLOTYPE if haps == 1 else "", panel.name,
                           hap_total, haps, sorted(cis), expected, af_product, af_min,
                           n_variants=k)
        return Verdict(TIER_A, PHASE_NOT_KNOWN, panel.name, hap_total, 0, sorted(unknown),
                       None, af_product, af_min, n_variants=k)
    same = same_score_verdict(variants, panels, span, alignment)
    if same is not None:
        return same
    panel, gts, _h, _c, _u, carrying = max(full, key=lambda f: f[5])
    if carrying:
        # haplotypes carry them (together), each with an indel that changes the target
        return Verdict(TIER_C1, INDEL_IN_TARGET, panel.name, panel.haplotypes(pos),
                       carrying, n_variants=k)
    expected, panel = max(
        ((p.haplotypes(pos) * math.prod(p.allele_freq(v, gv) for v, gv in zip(variants, g)), p)
         for p, g, _h, _c, _u, _t in full), key=lambda e: e[0])
    state = TOO_RARE if (min_expected and expected < min_expected) else ""
    return Verdict(TIER_C1, state, panel.name, panel.haplotypes(pos), 0, (), expected,
                   af_product, af_min, n_variants=k)


# --------------------------------------------------------------------------- #
# The gate over a search's integrated results
# --------------------------------------------------------------------------- #
def _read_rows(results, pam=None, unlisted_only=False):
    """{chromosome: [Alignment]} of every alignment of the integrated results that uses a
    variant (row: its 1-based data line; cfd: the CFD the row shows, where the alignment
    has one -- the ranking by CRISPR-Bulge score has no CFD column; pam: the search's
    PAM). ``unlisted_only``: only the alignments that name no carrier (their
    Variant_samples empty), a sites-only dataset's, in a per-sample search whose other
    rows name their cis carriers."""
    by_chrom = {}
    with open(results) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        gi, ci = ix["Spacer+PAM"], ix["Chromosome"]
        cols = []
        for block in BLOCKS:
            try:
                cols.append((block, [ix[f"{c}_({block})"] for c in _BLOCK_COLUMNS],
                             ix.get(f"CFD_score_({block})"),
                             ix.get(f"Variant_samples_({block})")))
            except KeyError:
                continue
        for row, line in enumerate(fh, start=1):
            f = line.rstrip("\n").split("\t")
            for block, (si, sti, ci_, ai, pi, oi, vi, mi), fi, smi in cols:
                if vi >= len(f) or str(f[oi]).strip().lower() != "alt":
                    continue
                if unlisted_only and smi is not None and smi < len(f) \
                        and f[smi].strip() not in ("", "NA", "n", ".", "-"):
                    continue
                variants = parse_variants(f[vi])
                if not variants:
                    continue
                cfd = None
                if fi is not None and fi < len(f):
                    try:
                        cfd = float(f[fi])
                    except ValueError:
                        cfd = None
                by_chrom.setdefault(f[ci], []).append(Alignment(
                    row, block, f[gi], f[si], f[sti], f[ci_].upper(), f[ai],
                    len(f[pi].strip()), variants, f[mi], cfd, pam))
    return by_chrom


def _maf_values(field):
    out = []
    for tok in str(field).split(","):
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out


def _gate_chromosome(job):
    """Tier every alignment of one chromosome's rows: (chrom, lines, counts, panel facts)."""
    chrom, rows, dictionaries, panel_specs, datasets, min_expected, genome, bound = job
    seq_at = None
    fasta = os.path.join(genome, f"{chrom}.fa") if genome else None
    if fasta and os.path.isfile(fasta):
        with open(fasta) as fh:
            ref = "".join(ln.strip() for ln in fh if not ln.startswith(">")).upper()
        seq_at = lambda a, b, _s=ref: _s[max(a, 0):max(b, 0)]  # noqa: E731
    wanted = set()
    for aln in rows:
        wanted.update(v[1] - 1 for v in aln.variants)
        if all(is_snp(v) for v in aln.variants):
            try:
                span = target_span(aln.start, aln.aligned)
            except ValueError:
                span = None
            if span is not None:
                wanted.update(range(span[0], span[1] + 1))
    panels = [PanelChrom(name, root, chrom, role, wanted, seq_at)
              for name, root, role in panel_specs]
    panels = [p for p in panels if p.gt is not None or p.indels is not None]
    registries = []
    for name in datasets:
        rb = os.path.join(dictionaries, f"registry_{name}", f"reg_{chrom}.bin")
        ri = os.path.join(dictionaries, f"registry_{name}", f"reg_{chrom}.idx")
        if os.path.isfile(rb) and os.path.isfile(ri):
            try:
                registries.append(t0.RegistryReader(rb, ri))
            except Exception:  # noqa: BLE001 - a frequency for the bound only
                pass

    # an indel's, of a sites-only dataset's indel frequencies (build_mega_indel_af: AF_max,
    # as its registry gives for a SNP); a genotyped dataset's rows carry theirs (MAF)
    at_indels = {v[1] for aln in rows for v in aln.variants if not is_snp(v)}
    indel_afs = []
    for name in datasets:
        path = os.path.join(dictionaries, f"indel_af_{name}", f"indel_af_{chrom}.tsv.gz")
        if at_indels and os.path.isfile(path):
            try:
                indel_afs.append(bia.IndelAfReader(path, at_indels))
            except Exception:  # noqa: BLE001 - a frequency for the bound only
                pass

    def af_of(v):
        if not is_snp(v):
            for reader in indel_afs:
                rec = reader.lookup(v[1], v[2], v[3])
                if rec and "AF_max" in rec:
                    return rec["AF_max"]
            return None
        for reg in registries:
            groups = reg.lookup(v[1], v[3])
            if groups and t0.GLOBAL_GROUP_ID in groups:
                return groups[t0.GLOBAL_GROUP_ID].allele_freq()
        return None

    lines, counts, cache = [], {}, {}
    for aln in rows:
        variants = aln.variants
        span = None
        if all(is_snp(v) for v in variants):
            try:
                span = target_span(aln.start, aln.aligned)
            except ValueError:
                span = None
        # the alignment decides the tier too (same_score_verdict), so two rows share a
        # verdict only when they take the same alleles over the same aligned bases, and
        # are both ranked by a score the gate reads (the CRISPR-Bulge one is not)
        key = (tuple(variants), span, aln.crrna, aln.aligned, aln.cfd,
               aln.block == CB_BLOCK)
        verdict = cache.get(key)
        if verdict is None:
            fallback = _maf_values(aln.maf)
            if panels:
                freq_of = af_of if registries or indel_afs else None
                verdict = classify(variants, panels, span, min_expected, freq_of, aln)
            else:
                # no phased panel at all: one variant is B, more are C2 (no panel)
                verdict = Verdict(TIER_B if len(variants) == 1 else TIER_C2,
                                  NOT_IN_PANEL if len(variants) == 1 else "",
                                  panel_haplotypes=None if len(variants) > 1 else 0,
                                  missing=len(variants), n_variants=len(variants))
                if len(variants) >= 2 and (registries or indel_afs):
                    afs = [af_of(v) for v in variants]
                    if all(a is not None for a in afs):
                        verdict.af_product, verdict.af_min = math.prod(afs), min(afs)
            if verdict.af_product is None and len(variants) >= 2 \
                    and len(fallback) == len(variants):
                verdict.af_product = math.prod(fallback)
                verdict.af_min = min(fallback)
            # D6: in a search of a sites-only dataset a C2 combination below the frequency
            # bound is ranked after the others; never dropped
            if bound and verdict.tier == TIER_C2 and below_bound(verdict.af_product,
                                                                 verdict.af_min):
                verdict.state = BELOW_BOUND
            cache[key] = verdict
        if aln.block == BLOCKS[0]:
            counts[verdict.tier] = counts.get(verdict.tier, 0) + 1
            if verdict.state in SAME_SCORE_STATES:
                counts[SAME_SCORE_COUNT] = counts.get(SAME_SCORE_COUNT, 0) + 1
            if verdict.state == BELOW_BOUND:
                counts[BELOW_BOUND_COUNT] = counts.get(BELOW_BOUND_COUNT, 0) + 1
        lines.append((aln.row, BLOCKS.index(aln.block), "\t".join((
            str(aln.row), aln.block, aln.guide, chrom, str(aln.start), aln.strand,
            ",".join(_token(v) for v in variants),
            verdict.tier, verdict.state or "-", verdict.panel or "-",
            str(verdict.panel_haplotypes or "NA"),
            "NA" if verdict.haplotypes is None else str(verdict.haplotypes),
            ",".join(verdict.samples) or "-",
            "NA" if verdict.expected is None else "%.4g" % verdict.expected,
            "NA" if verdict.af_product is None else "%.4g" % verdict.af_product,
            "NA" if verdict.af_min is None else "%.4g" % verdict.af_min,
            verdict.label(), ",".join(verdict.kept) or "-",
        ))))
    facts = {p.name: {"samples": p.samples(), "haplotypes": p.haplotypes(1) or 2 * p.samples(),
                      "role": p.role, "indels": p.indels is not None,
                      "superpopulations": p.superpopulations(),
                      "ref_mismatches": p.ref_mismatches}
             for p in panels if p.gt is not None}
    for p in panels:
        p.close()
    for reg in registries:
        try:
            reg.close()
        except Exception:  # noqa: BLE001
            pass
    return chrom, lines, counts, facts


def run_gate(results, out_path, dictionaries, panel_specs, datasets,
             min_expected=MIN_EXPECTED, threads=1, genome=None, pam=None, bound=False,
             unlisted_only=False):
    """Write the companion of ``results`` to ``out_path``; return (counts, panel facts).
    ``panel_specs``: [(panel name, Dictionaries folder holding its stores, role)];
    ``genome``: the reference genome folder (<chrom>.fa), for holds_target; ``pam``: the
    search's PAM (Pam), read over an alignment's PAM region; ``bound``: the search has a
    sites-only dataset, whose C2 rows below the frequency bound are ranked after the others
    (D6); ``unlisted_only``: only the alignments that name no carrier (_read_rows)."""
    by_chrom = _read_rows(results, pam, unlisted_only)
    jobs = [(chrom, rows, dictionaries, panel_specs, datasets, min_expected, genome, bound)
            for chrom, rows in sorted(by_chrom.items(), key=lambda kv: -len(kv[1]))]
    if threads > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(min(threads, len(jobs))) as pool:
            done = pool.map(_gate_chromosome, jobs, chunksize=1)
    else:
        done = [_gate_chromosome(job) for job in jobs]
    lines, counts, facts = [], {}, {}
    for _chrom, chrom_lines, chrom_counts, chrom_facts in done:
        lines.extend(chrom_lines)
        for tier, n in chrom_counts.items():
            counts[tier] = counts.get(tier, 0) + n
        for name, fact in chrom_facts.items():
            seen = facts.get(name)
            if seen is None:
                facts[name] = dict(fact)
                continue
            seen["ref_mismatches"] += fact["ref_mismatches"]
            seen["haplotypes"] = max(seen["haplotypes"], fact["haplotypes"])
            seen["indels"] = seen["indels"] or fact["indels"]
    lines.sort(key=lambda t: (t[0], t[1]))
    tmp = out_path + ".tmp"
    with open(tmp, "w") as out:
        out.write(
            "# CRISPRme carrier gate. One line per alignment of a row of the integrated results\n"
            "# that uses a variant (Row: the row's line, header excluded; Block: which alignment).\n"
            "# Tier A: two or more variants one phased haplotype of the panel carries together,\n"
            "# or, of an alignment ranked by CFD, that a haplotype reaches its CFD without\n"
            "# (Scoring_variants: those its CFD needs, which that haplotype carries; the others\n"
            "# leave the CFD the results report, the alignment's mismatches and its PAM as they\n"
            "# are).\n"
            "# Tier B: one variant. Tier C1: two or more variants, each in the panel, that no\n"
            "# haplotype of it carries together (or a variant whose every carrying haplotype has\n"
            "# an indel in the target): the report ranks these after every other row. Tier C2:\n"
            "# a variant not in the panel, so the panel cannot say. Expected: the panel\n"
            "# haplotypes expected to carry the variants together if they were independent.\n"
            "# AF_product / AF_min: the frequency bound (>= 1e-06 / >= 0.001). In a search of\n"
            "# a sites-only dataset a C2 alignment below it (State: below the frequency bound)\n"
            "# is ranked after every other row too; elsewhere it is shown only.\n"
        )
        out.write("#" + "\t".join(COMPANION_COLUMNS) + "\n")
        for _row, _b, text in lines:
            out.write(text + "\n")
    os.replace(tmp, out_path)
    return counts, facts


# --------------------------------------------------------------------------- #
# Which panels a search's gate reads
# --------------------------------------------------------------------------- #
def _registry_manifest(dictionaries, name):
    """The manifest of a dataset's registry (any chromosome's), or {}."""
    rdir = os.path.join(dictionaries, f"registry_{name}")
    try:
        for f in sorted(os.listdir(rdir)):
            if f.endswith(".idx"):
                with open(os.path.join(rdir, f)) as fh:
                    return json.load(fh)
    except (OSError, ValueError):
        pass
    return {}


def dataset_kind(dictionaries, name):
    """"genotyped" when the dataset has a phased SNP store, "sites-only" when its registry
    says it has no genotypes, "untiered" when it has the classic per-sample dictionaries
    alone (an index built without --samplesID), else "other"."""
    if os.path.isdir(os.path.join(dictionaries, f"genotypes_{name}")):
        return "genotyped"
    man = _registry_manifest(dictionaries, name)
    if man.get("aggregation") == "info_af" or man.get("data_type") == "sites-only":
        return "sites-only"
    if os.path.isdir(os.path.join(dictionaries, f"dictionaries_{name}")):
        return "untiered"
    return "other"


def find_oracle(dictionaries, given=None, panel=ORACLE_PANEL):
    """(name, Dictionaries folder) of the phased store a sites-only search is checked
    against, or None: ``given`` (CRISPRME_GENOTYPE_TIER: a folder holding genotypes_<panel>/,
    or that folder itself), else the install's own Dictionaries/genotypes_<panel>/."""
    for cand in ([given] if given else []) + [dictionaries]:
        cand = os.path.abspath(cand)
        base = os.path.basename(os.path.normpath(cand))
        if base.startswith("genotypes_") and os.path.isdir(cand):
            return base[len("genotypes_"):], os.path.dirname(cand)
        if os.path.isdir(os.path.join(cand, f"genotypes_{panel}")):
            return panel, cand
    return None


def plan_panels(dictionaries, datasets, mode, prioritize, oracle=None):
    """(status, [(panel name, Dictionaries folder, role)]) of a search's gate.

    ``mode``: "population-level" or "per-sample"; ``prioritize``: "auto" (the default: a
    population-level search of a genotyped dataset is checked against its own store),
    "on" (--prioritize-observed: a sites-only dataset is checked against the oracle too,
    in either mode), "off". Status: "run", "off", "per-sample" (nothing to check: the
    per-sample rows name their cis carriers), "bound" (a sites-only dataset with no phased
    panel to read -- not asked for the oracle, or asked and the oracle not installed: its
    rows are checked against the frequency bound alone, D6), "untiered" (a genotyped
    dataset with the classic dictionaries alone: the gate reads the phased store
    build-index-only --samplesID writes), "no-panel" (no phased store). A per-sample search
    of a sites-only dataset beside a genotyped one reads the genotyped datasets' stores for
    the sites-only rows (gate_scope: the rows that name no carrier), as a population-level
    one does. (Up to the plan's D6 a sites-only dataset not asked for the oracle had
    status "opt-in", and one asked without it "no-oracle": nothing was run.)"""
    if prioritize == "off":
        return "off", []
    kinds = {d: dataset_kind(dictionaries, d) for d in datasets}
    genotyped = [d for d in datasets if kinds[d] == "genotyped"]
    sites_only = [d for d in datasets if kinds[d] == "sites-only"]
    panels = []
    if mode == "population-level" or prioritize == "on" or sites_only:
        panels.extend((d, dictionaries, "own") for d in genotyped)
    if prioritize == "on" and sites_only:
        found = find_oracle(dictionaries, oracle)
        if found is not None and not any(
                os.path.realpath(os.path.join(root, f"genotypes_{n}"))
                == os.path.realpath(os.path.join(found[1], f"genotypes_{found[0]}"))
                for n, root, _r in panels):
            panels.append((found[0], found[1], "oracle"))
    if panels:
        return "run", panels
    if sites_only:
        return "bound", []
    if mode == "per-sample" and genotyped:
        return "per-sample", []
    if mode == "population-level" and any(kinds[d] == "untiered" for d in datasets):
        return "untiered", []
    return "no-panel", []


def gate_scope(dictionaries, datasets, mode, prioritize):
    """(bound, unlisted_only) of a search's gate: ``bound``, the search has a sites-only
    dataset, whose C2 rows below the frequency bound are ranked after the others (D6);
    ``unlisted_only``, a per-sample search of one, whose gate reads only the rows that
    name no carrier (the sites-only dataset's: the others name their cis carriers)."""
    if prioritize == "off":
        return False, False
    sites_only = any(dataset_kind(dictionaries, d) == "sites-only" for d in datasets)
    return sites_only, sites_only and mode == "per-sample"


def oracle_missing(dictionaries, datasets, prioritize, oracle=None):
    """True when --prioritize-observed asks for the genotype oracle a sites-only dataset
    is checked against, and it is not installed."""
    return (prioritize == "on"
            and any(dataset_kind(dictionaries, d) == "sites-only" for d in datasets)
            and find_oracle(dictionaries, oracle) is None)


def companion_path(results):
    base = results[:-len("_integrated_results.tsv")] if results.endswith(
        "_integrated_results.tsv") else os.path.splitext(results)[0]
    return base + COMPANION_SUFFIX


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="The carrier gate: tier a search's variant combinations against "
                    "phased haplotypes (writes <results>_carrier_gate.tsv)")
    parser.add_argument("--results", required=True, help="the integrated results TSV")
    parser.add_argument("--dictionaries", required=True, help="the install's Dictionaries/")
    parser.add_argument("--vcf-list", required=True,
                        help="the search's list of VCF folders ('_': none)")
    parser.add_argument("--mode", default="population-level",
                        choices=("population-level", "per-sample"))
    parser.add_argument("--prioritize-observed", default="auto", choices=("auto", "on", "off"))
    parser.add_argument("--oracle", default="",
                        help="a folder holding genotypes_<panel>/ (CRISPRME_GENOTYPE_TIER)")
    parser.add_argument("--min-expected", type=float, default=MIN_EXPECTED)
    parser.add_argument("--genome", default="",
                        help="the search's reference genome folder (<chrom>.fa): an indel "
                             "that leaves a target's bases (a repeat unit) takes no slot off")
    parser.add_argument("--pam", default="",
                        help="the search's PAM file: an allele a row's PAM needs is one its "
                             "score needs. Without it every allele in a PAM region is")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--meta", default=None, help="the sidecar JSON (default: "
                        "<results folder>/" + META_NAME + ")")
    args = parser.parse_args(argv)
    datasets = []
    try:
        with open(args.vcf_list) as fh:
            for line in fh:
                name = os.path.basename(os.path.normpath(line.strip()))
                if line.strip() and line.strip() != "_" and name not in datasets:
                    datasets.append(name)
    except OSError:
        pass
    meta_path = args.meta or os.path.join(os.path.dirname(os.path.abspath(args.results)),
                                          META_NAME)
    status, panels = plan_panels(args.dictionaries, datasets, args.mode,
                                 args.prioritize_observed, args.oracle or None)
    bound, unlisted = gate_scope(args.dictionaries, datasets, args.mode,
                                 args.prioritize_observed)
    no_oracle = oracle_missing(args.dictionaries, datasets, args.prioritize_observed,
                               args.oracle or None)
    meta = {"status": status, "min_expected": args.min_expected, "panels": [],
            "counts": {}, "companion": None, "mode": args.mode,
            "prioritize_observed": args.prioritize_observed,
            "frequency_bound": bound, "rows": "no carrier listed" if unlisted else "all",
            "oracle_missing": no_oracle}
    if no_oracle:
        print("Carrier gate: --prioritize-observed asked to check the sites-only rows against "
              f"the phased {ORACLE_PANEL} genotypes, which are not installed "
              f"(Dictionaries/genotypes_{ORACLE_PANEL}/, or CRISPRME_GENOTYPE_TIER). Install "
              "them with: crisprme.py download --what genotype-oracle", flush=True)
    if status in ("run", "bound"):
        out = companion_path(args.results)
        counts, facts = run_gate(args.results, out, args.dictionaries, panels, datasets,
                                 args.min_expected, max(1, args.threads),
                                 args.genome or None,
                                 read_pam(args.pam) if args.pam else None, bound, unlisted)
        meta["counts"] = counts
        meta["companion"] = os.path.basename(out)
        for name, root, role in panels:
            fact = facts.get(name, {"role": role})
            fact["name"] = name
            fact["store"] = os.path.join(root, f"genotypes_{name}")
            meta["panels"].append(fact)
        named = ", ".join(f"{p['name']} ({p.get('samples', 0):,} samples)"
                          for p in meta["panels"])
        same = counts.get(SAME_SCORE_COUNT, 0)
        below = counts.get(BELOW_BOUND_COUNT, 0)
        scope = (" (the rows that name no carrier: the sites-only dataset's)"
                 if unlisted else "")
        bound_text = (f", {below:,} of them below the frequency bound (product of the "
                      f"variants' frequencies >= 1e-06, rarest >= 0.001) and ranked after the "
                      f"others" if bound else "")
        if status == "run":
            print(f"Carrier gate: the rows' variant combinations{scope} checked against the "
                  f"phased genotypes of {named}: A {counts.get(TIER_A, 0):,} ({same:,} of "
                  f"them on a haplotype that reaches the row's CFD without every variant), "
                  f"B {counts.get(TIER_B, 0):,}, C1 {counts.get(TIER_C1, 0):,} (ranked after "
                  f"the others), C2 {counts.get(TIER_C2, 0):,}{bound_text} (by the "
                  f"highest-CFD alignment of each row); {os.path.basename(out)}", flush=True)
        else:
            print(f"Carrier gate: no phased panel to check the sites-only rows{scope} "
                  f"against (complete-search --prioritize-observed reads the phased "
                  f"{ORACLE_PANEL} genotypes): their combinations of two or more variants "
                  f"were checked against the frequency bound alone: C2 "
                  f"{counts.get(TIER_C2, 0):,}{bound_text}, B {counts.get(TIER_B, 0):,} "
                  f"(one variant) (by the highest-CFD alignment of each row); "
                  f"{os.path.basename(out)}", flush=True)
        mism = sum(p.get("ref_mismatches", 0) for p in meta["panels"])
        if mism:
            print(f"Carrier gate: {mism:,} variant lookup(s) found the panel's REF different "
                  f"from the row's: those variants were taken as not in the panel (the "
                  f"dataset and the panel were not normalized alike there)", flush=True)
    elif status == "untiered":
        print("Carrier gate: not run: the search's genotyped datasets have the classic "
              "per-sample dictionaries alone, not the phased genotype store it reads "
              "(Dictionaries/genotypes_<vcf>/, which build-index-only --samplesID writes, and "
              "a downloaded index has); the rows are reported unchecked.", flush=True)
    tmp = meta_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    os.replace(tmp, meta_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
