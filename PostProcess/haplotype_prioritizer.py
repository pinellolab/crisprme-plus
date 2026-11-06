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

  A   two or more variants that one phased haplotype of the panel carries together. States:
      "one haplotype" when a single haplotype does, "phase not known" when only samples whose
      genotypes do not say which haplotype holds which variant carry them all;
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
      >= 0.001) is given with them, as information, and applied to nothing.

A haplotype holds a row only when it holds the row's target: on a phased slot that carries
the row's SNPs, an indel of the store that changes the bases the target covers (an indel
another of the slot's alleles makes impossible -- one at a listed SNP's own position, or one
whose REF covers it -- is a phasing artefact and changes nothing) means the slot does not have
the row's sequence. An unphased call keeps the carrier ("phase not known"), as a slot that
holds two ALTs at one position does (the PUTATIVE rule for two ALTs on one slot).

The gate changes no row and drops none: it writes a companion TSV (<results>_carrier_gate.tsv,
one line per alignment that uses a variant) and a sidecar the report reads
(.carrier_gate.json). The report ranks the C1 rows after the others, in a tier of their own,
with the factual label of each tier. The same helpers drop, in the per-sample SNP
post-analysis (new_simple_analysis.py), a carrier slot that holds such an indel.

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

import indel_snp_cis as isc
import population_summary as ps
import tier0_compile as t0c
import tier0_registry as t0
import tier1_genotypes as t1

TIER_A, TIER_B, TIER_C1, TIER_C2 = "A", "B", "C1", "C2"
ONE_HAPLOTYPE = "one haplotype"
PHASE_NOT_KNOWN = "phase not known"
TOO_RARE = "too rare to tell"
INDEL_IN_TARGET = "indel in the target"
NOT_IN_PANEL = "not in the panel"

# D2: below this many expected haplotypes a zero says little ("too rare to tell"). A
# parameter (complete-search --gate-min-expected); 0 gives every C1 row the same state.
MIN_EXPECTED = 3.0
# The safeguards' frequency bound, shown for a multi-variant row and applied to nothing
# (D6): the product of the variants' frequencies and the rarest of them.
FREQ_BOUND_PRODUCT = 1e-6
FREQ_BOUND_RAREST = 0.001

# the oracle a sites-only search is checked against (download --what genotype-oracle)
ORACLE_PANEL = "hg38_HGDP1kGP"

# the integrated results' three alignments of a row, by the suffix of their columns
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
_BLOCK_COLUMNS = ("Start_coordinate", "Strand", "Aligned_protospacer+PAM_ALT",
                  "REF/ALT_origin", "Variant_info_genome", "Variant_MAF")

COMPANION_SUFFIX = "_carrier_gate.tsv"
META_NAME = ".carrier_gate.json"
COMPANION_COLUMNS = (
    "Row", "Block", "Spacer+PAM", "Chromosome", "Start", "Strand", "Variants", "Tier",
    "State", "Panel", "Panel_haplotypes", "Haplotypes", "Samples", "Expected",
    "AF_product", "AF_min", "Label",
)


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


def slots_holding(variant_gts, others=(), breakers=()):
    """(held, broken) for ONE sample: of its slots that hold every allele of
    ``variant_gts`` ([genotype], each read against its allele) and none of ``others`` (its
    other alleles in the window: such a slot is another haplotype), the number with no
    breaker (an indel that changes the target, ``breakers``: [genotype]) and the number with
    one. None when the genotypes cannot put the alleles on slots (an unphased or a missing
    call, haploid against diploid): such a carrier stays listed. A breaker whose call cannot
    be put on a slot takes no slot off."""
    parsed = [_tokens(g) for g in variant_gts]
    n = _exact(parsed)
    if n is None:
        return None
    other = [_tokens(g) for g in others]
    if other and _exact(parsed + other) != n:
        return None
    held = broken = 0
    for h in range(n):
        if not all(toks[h] == "1" for toks, _kind in parsed):
            continue
        if any(toks[h] == "1" for toks, _kind in other):
            continue
        hit = False
        for g in breakers:
            bt = _tokens(g)
            if _exact(parsed + [bt]) == n and bt[0][h] == "1":
                hit = True
                break
        if hit:
            broken += 1
        else:
            held += 1
    return held, broken


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


def breaking_indels(index, lo, hi, snps, own=()):
    """{sample: genotype} of each indel of ``index`` that changes the reference bases
    lo..hi (0-based, inclusive) of a haplotype carrying it (indel_snp_cis.changes_span),
    but those of ``own`` and those a slot cannot hold with the SNPs ``snps`` ([(pos, REF,
    ALT)]): an indel at a SNP's own position, or one whose REF covers it, is another allele
    of that base, and a slot phased with both is a phasing artefact."""
    out = []
    if index is None:
        return out
    snp_pos = {p for p, _r, _a in snps}
    for p, r, a in index.near(lo, hi):
        if (p, r, a) in own:
            continue
        if any(p <= q <= p + len(r) - 1 for q in snp_pos):
            continue
        if not isc.changes_span(p, r, a, lo, hi):
            continue
        out.append(index.carriers(p, r, a))
    return out


class PanelChrom(object):
    """One panel's phased genotypes on one chromosome: the SNP store (Tier-1) and the indel
    store, either absent. ``carriers`` gives a variant's carriers with each genotype read
    against that variant ("1|0"), or None for a variant the panel does not hold.
    ``wanted``: the positions the rows touch, for the indel store (IndelIndex)."""

    def __init__(self, name, dictionaries, chrom, role="own", wanted=None):
        self.name = name
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
    ``expected`` (E), the frequency bound's ``af_product`` / ``af_min`` and ``missing``,
    the number of variants no panel holds."""

    __slots__ = ("tier", "state", "panel", "panel_haplotypes", "haplotypes", "samples",
                 "expected", "af_product", "af_min", "missing", "n_variants")

    def __init__(self, tier, state="", panel=None, panel_haplotypes=0, haplotypes=None,
                 samples=(), expected=None, af_product=None, af_min=None, missing=0,
                 n_variants=0):
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

    def label(self):
        """The factual text the report shows (ASCII: it goes into the TSV downloads too)."""
        hap_of = f"{self.panel_haplotypes:,} {self.panel} haplotypes"
        if self.tier == TIER_A:
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
        text = (f"not checked: {self.missing} of {self.n_variants} variants not in "
                f"{self.panel or 'the panel'}")
        if self.af_product is not None:
            inside = (self.af_product >= FREQ_BOUND_PRODUCT
                      and self.af_min >= FREQ_BOUND_RAREST)
            text += (f"; frequency product {self.af_product:.1e}, rarest {self.af_min:.1e}"
                     f" ({'within' if inside else 'below'} the 1e-06 / 0.001 bound)")
        return text


def classify(variants, panels, span=None, min_expected=MIN_EXPECTED, af_of=None):
    """The tier of one alignment's variants against the panels (PanelChrom).

    ``span`` (lo, hi): the target's reference bases, for an alignment of SNPs only: a slot
    that holds them with an indel changing those bases does not hold the target. ``af_of``:
    the dataset's frequency of a variant (its registry), or None, for the frequency bound of
    a multi-variant alignment. A variant is in a panel when the panel holds its POS, REF and
    ALT; the panels are the stores a search reads (a row does not say which dataset it comes
    from): the first with evidence decides, the one with the most haplotypes holding the
    variants together for A, the most expected for C1."""
    k = len(variants)
    if k == 0:
        return None
    snps_only = all(is_snp(v) for v in variants)
    afs = None
    if k >= 2 and af_of is not None:
        afs = [af_of(v) for v in variants]
        if any(a is None for a in afs):
            afs = None
    af_product = math.prod(afs) if afs else None
    af_min = min(afs) if afs else None
    pos = variants[0][1]
    full, missing = [], k
    for panel in panels:
        gts = [panel.carriers(v) for v in variants]
        absent = sum(1 for g in gts if not g)
        missing = min(missing, absent)
        if absent:
            continue
        breakers = []
        if snps_only and span is not None:
            breakers = breaking_indels(panel.indels, span[0], span[1],
                                       [(v[1], v[2], v[3]) for v in variants])
        if k == 1 and not breakers:
            # one variant, and nothing in the target to take a slot off it
            return Verdict(TIER_B, panel=panel.name, panel_haplotypes=panel.haplotypes(pos),
                           haplotypes=None, n_variants=1)
        rivals = []
        for v in variants:
            rivals.extend(panel.rivals(v))
        haps, cis, unknown, carrying = support(gts, breakers, rivals)
        full.append((panel, gts, haps, cis, unknown, carrying))
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
def _read_rows(results):
    """{chromosome: [(row, block, guide, start, strand, aligned, variants, maf)]} of every
    alignment of the integrated results that uses a variant (row: its 1-based data line)."""
    by_chrom = {}
    with open(results) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        gi, ci = ix["Spacer+PAM"], ix["Chromosome"]
        cols = []
        for block in BLOCKS:
            try:
                cols.append((block, [ix[f"{c}_({block})"] for c in _BLOCK_COLUMNS]))
            except KeyError:
                continue
        for row, line in enumerate(fh, start=1):
            f = line.rstrip("\n").split("\t")
            for block, (si, sti, ai, oi, vi, mi) in cols:
                if vi >= len(f) or str(f[oi]).strip().lower() != "alt":
                    continue
                variants = parse_variants(f[vi])
                if not variants:
                    continue
                by_chrom.setdefault(f[ci], []).append(
                    (row, block, f[gi], f[si], f[sti], f[ai], variants, f[mi]))
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
    chrom, rows, dictionaries, panel_specs, datasets, min_expected = job
    wanted = set()
    for _row, _block, _guide, start, _strand, aligned, variants, _maf in rows:
        wanted.update(v[1] - 1 for v in variants)
        if all(is_snp(v) for v in variants):
            try:
                span = target_span(start, aligned)
            except ValueError:
                span = None
            if span is not None:
                wanted.update(range(span[0], span[1] + 1))
    panels = [PanelChrom(name, root, chrom, role, wanted) for name, root, role in panel_specs]
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

    def af_of(v):
        if not is_snp(v):
            return None
        for reg in registries:
            groups = reg.lookup(v[1], v[3])
            if groups and t0.GLOBAL_GROUP_ID in groups:
                return groups[t0.GLOBAL_GROUP_ID].allele_freq()
        return None

    lines, counts, cache = [], {}, {}
    for row, block, guide, start, strand, aligned, variants, maf in rows:
        span = None
        if all(is_snp(v) for v in variants):
            try:
                span = target_span(start, aligned)
            except ValueError:
                span = None
        key = (tuple(variants), span)
        verdict = cache.get(key)
        if verdict is None:
            if panels:
                fallback = _maf_values(maf)
                bound = af_of if registries else None
                verdict = classify(variants, panels, span, min_expected, bound)
                if verdict.af_product is None and len(variants) >= 2 \
                        and len(fallback) == len(variants):
                    verdict.af_product = math.prod(fallback)
                    verdict.af_min = min(fallback)
            else:
                verdict = Verdict(TIER_B if len(variants) == 1 else TIER_C2,
                                  NOT_IN_PANEL if len(variants) == 1 else "",
                                  missing=len(variants), n_variants=len(variants))
            cache[key] = verdict
        if block == BLOCKS[0]:
            counts[verdict.tier] = counts.get(verdict.tier, 0) + 1
        lines.append((row, BLOCKS.index(block), "\t".join((
            str(row), block, guide, chrom, str(start), strand,
            ",".join(f"{c}_{p}_{r}_{a}" for c, p, r, a in variants),
            verdict.tier, verdict.state or "-", verdict.panel or "-",
            str(verdict.panel_haplotypes or "NA"),
            "NA" if verdict.haplotypes is None else str(verdict.haplotypes),
            ",".join(verdict.samples) or "-",
            "NA" if verdict.expected is None else "%.4g" % verdict.expected,
            "NA" if verdict.af_product is None else "%.4g" % verdict.af_product,
            "NA" if verdict.af_min is None else "%.4g" % verdict.af_min,
            verdict.label(),
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
             min_expected=MIN_EXPECTED, threads=1):
    """Write the companion of ``results`` to ``out_path``; return (counts, panel facts).
    ``panel_specs``: [(panel name, Dictionaries folder holding its stores, role)]."""
    by_chrom = _read_rows(results)
    jobs = [(chrom, rows, dictionaries, panel_specs, datasets, min_expected)
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
            "# Tier A: two or more variants one phased haplotype of the panel carries together.\n"
            "# Tier B: one variant. Tier C1: two or more variants, each in the panel, that no\n"
            "# haplotype of it carries together (or a variant whose every carrying haplotype has\n"
            "# an indel in the target): the report ranks these after every other row. Tier C2:\n"
            "# a variant not in the panel, so the panel cannot say. Expected: the panel\n"
            "# haplotypes expected to carry the variants together if they were independent.\n"
            "# AF_product / AF_min: the frequency bound (>= 1e-06 / >= 0.001), shown only.\n"
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
    per-sample rows name their cis carriers), "opt-in" (sites-only only, not asked),
    "untiered" (a genotyped dataset with the classic dictionaries alone: the gate reads the
    phased store build-index-only --samplesID writes), "no-panel" (no phased store),
    "no-oracle" (asked, the oracle not installed)."""
    if prioritize == "off":
        return "off", []
    kinds = {d: dataset_kind(dictionaries, d) for d in datasets}
    genotyped = [d for d in datasets if kinds[d] == "genotyped"]
    sites_only = [d for d in datasets if kinds[d] == "sites-only"]
    panels = []
    if mode == "population-level" or prioritize == "on":
        panels.extend((d, dictionaries, "own") for d in genotyped)
    if prioritize == "on" and sites_only:
        found = find_oracle(dictionaries, oracle)
        if found is None:
            if not panels:
                return "no-oracle", []
        elif not any(os.path.realpath(os.path.join(root, f"genotypes_{n}"))
                     == os.path.realpath(os.path.join(found[1], f"genotypes_{found[0]}"))
                     for n, root, _r in panels):
            panels.append((found[0], found[1], "oracle"))
    if panels:
        return "run", panels
    if mode == "per-sample" and genotyped:
        return "per-sample", []
    if sites_only and prioritize != "on":
        return "opt-in", []
    if mode == "population-level" and any(kinds[d] == "untiered" for d in datasets):
        return "untiered", []
    return "no-panel", []


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
    meta = {"status": status, "min_expected": args.min_expected, "panels": [],
            "counts": {}, "companion": None, "mode": args.mode,
            "prioritize_observed": args.prioritize_observed}
    if status == "run":
        out = companion_path(args.results)
        counts, facts = run_gate(args.results, out, args.dictionaries, panels, datasets,
                                 args.min_expected, max(1, args.threads))
        meta["counts"] = counts
        meta["companion"] = os.path.basename(out)
        for name, root, role in panels:
            fact = facts.get(name, {"role": role})
            fact["name"] = name
            fact["store"] = os.path.join(root, f"genotypes_{name}")
            meta["panels"].append(fact)
        named = ", ".join(f"{p['name']} ({p.get('samples', 0):,} samples)"
                          for p in meta["panels"])
        print(f"Carrier gate: the rows' variant combinations checked against the phased "
              f"genotypes of {named}: A {counts.get(TIER_A, 0):,}, B {counts.get(TIER_B, 0):,}, "
              f"C1 {counts.get(TIER_C1, 0):,} (ranked after the others), C2 "
              f"{counts.get(TIER_C2, 0):,} (by the highest-CFD alignment of each row); "
              f"{os.path.basename(out)}", flush=True)
        mism = sum(p.get("ref_mismatches", 0) for p in meta["panels"])
        if mism:
            print(f"Carrier gate: {mism:,} variant lookup(s) found the panel's REF different "
                  f"from the row's: those variants were taken as not in the panel (the "
                  f"dataset and the panel were not normalized alike there)", flush=True)
    elif status == "no-oracle":
        print("Carrier gate: --prioritize-observed asked to check the sites-only rows against "
              f"the phased {ORACLE_PANEL} genotypes, which are not installed "
              f"(Dictionaries/genotypes_{ORACLE_PANEL}/, or CRISPRME_GENOTYPE_TIER): the "
              "rows are reported unchecked. Install them with: crisprme.py download --what "
              "genotype-oracle", flush=True)
    elif status == "opt-in":
        print("Carrier gate: not run on a sites-only dataset unless asked "
              "(--prioritize-observed).", flush=True)
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
