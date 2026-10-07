# Why re-phasing mattered, and phased vs. mega

Paper narrative / internal notes.

**Guide** `ACTGAAATCTGTAAGCAGGC` · SpCas9 (NRG PAM) · hg38 · CFD + CRISPR-Bulge scores · all
functional/cancer annotations.

> Cohort note: the diversity cohort is **HGDP** (Human Genome Diversity Project), jointly called
> with the 1000 Genomes Project (1kGP) in the gnomAD HGDP+1kGP callset.

This document makes **two points**, each from a dedicated search:

- **Point 1 (search 6,1,1) — re-phasing was necessary, two ways.** (a) The old default mixed a phased
  1kGP with an *unphased* HGDP, forcing HGDP co-occurrences to PUTATIVE; re-phasing both together
  fixes that. (b) We could **not** simply use the *published* SHAPEIT5-phased release, because it drops
  singletons (MAC ≥ 2) — so we re-phased ourselves and **kept the singletons**, which matter for
  off-targets.
- **Point 2 (search 4,1,1) — mega over-calls putative off-targets.** A broad sites-only panel (mega)
  cannot confirm cis, so its top-ranked sites are inflated with putative (cis-unconfirmed) calls — a
  false-positive–prone ranking that the phased index cleans up.

---

## Numbers at a glance (reported + top-1000, by CFD and by CRISPR-Bulge)

For each comparison, two views: **how many off-targets are reported in total** (and how many reach each
score tier under each scorer), and **the top-1000** (overlap between the two indices + the Observed
breakdown), ranked **by CFD** and **by CRISPR-Bulge**.

### Comparison 1 (6,1,1) — phased HGDP+1kGP vs unphased hybrid

**1A · Reported off-targets (all), by score tier**

| index | total reported | CFD ≥0.5 | CFD ≥0.2 | CFD ≥0.1 | CRISPR-Bulge ≥0.5 | ≥0.2 | ≥0.1 |
|---|--:|--:|--:|--:|--:|--:|--:|
| **phased (new)** | 5,996,092 | 417 | 13,857 | 63,178 | 4 | 8 | 15 |
| hybrid (old) | 5,540,022 | 347 | 12,275 | 56,169 | 4 | 7 | 14 |

**1B · Top-1000 (overlap between the two indices + Observed), by scorer**

| ranked by | overlap | phased: carrier / ref / putative | hybrid: carrier / ref / putative |
|---|--:|---|---|
| **CFD** | **793 / 1000** | 478 / 522 / 0 | 388 / 612 / 0 |
| **CRISPR-Bulge** \* | 100 / 1000 | 453 / 547 / 0 | 527 / 473 / 0 |

### Comparison 2 (4,1,1) — phased HGDP+1kGP vs mega (sites-only)

**2A · Reported off-targets (all), by score tier**

| index | total reported | CFD ≥0.5 | CFD ≥0.2 | CFD ≥0.1 | CRISPR-Bulge ≥0.5 | ≥0.2 | ≥0.1 |
|---|--:|--:|--:|--:|--:|--:|--:|
| **phased** | 213,545 | 159 | 2,116 | 6,476 | 5 | 9 | 15 |
| mega (sites-only) | 184,449 | 220 | 2,061 | 5,785 | **33** | **53** | **73** |

**2B · Top-1000 (overlap + Observed), by scorer**

| ranked by | overlap | phased: carrier / ref / putative | mega: carrier / ref / putative / obs-AF |
|---|--:|---|---|
| **CFD** | 558 / 1000 | 578 / 422 / **0** | 0 / 399 / **273** / 328 |
| **CRISPR-Bulge** \* | 81 / 1000 | 521 / 479 / **0** | 0 / 369 / **186** / 445 |

> Mega's **high-CRISPR-Bulge counts are inflated** (33/53/73 at ≥0.5/≥0.2/≥0.1 vs the phased index's
> 5/9/15) for the same reason as its CFD=1.0 top putatives: the combinatorial variant stacks score high
> under *both* scorers. The mega-over-call story therefore holds under CFD **and** CRISPR-Bulge — 273
> (CFD) / 186 (CB) putative in its top-1000, vs **0** for the phased index either way.
>
> \* **CRISPR-Bulge top-1000 caveat.** CRISPR-Bulge produces very few high-scoring sites (only ~15 at
> CB ≥ 0.1 genome-wide), so its "top-1000" is dominated by a large near-zero tie region where rank
> order is near-arbitrary — the low CB overlap (100, 81) is tie-break noise, **not** index
> disagreement. CFD (which spreads scores) is the meaningful stability metric; the CB view is included
> for completeness and for the putative-load contrast, which *is* meaningful.

---

## POINT 1 — Re-phasing was necessary (search 6,1,1)

Phased HGDP+1kGP vs the old hybrid, identical search: **mm6 + 1 DNA + 1 RNA bulge**, `--per-sample`.

### 1a. Phasing the unphased HGDP: PUTATIVE → CONFIRMED

The previously shipped default (`NRG_3_hg38+hg38_1000G2021_HGDP`) is a **hybrid**: 1kGP phased, **HGDP
unphased**. A co-occurring off-target (two nearby variants that must sit on the *same* chromosome copy,
**in cis**) is reported **CONFIRMED** only when cis can be proven in a real individual — which
**requires phased genotypes**. On the hybrid, every HGDP-involving co-occurrence collapses to
**PUTATIVE** with no named carrier.

| phase_confirmation | CONFIRMED | PUTATIVE | % CONFIRMED |
|---|---|---|---|
| Hybrid (old, HGDP unphased) | 8,393,667 | 22,515,383 | 27.2% |
| **Re-phased (new)** | **36,319,150** | **971,014** | **97.4%** |

PUTATIVE collapses **23×**, CONFIRMED rises **4.3×**. Ground-truthed: off-target `chr13:100001388(−)`
(SNPs `chr13_100001399_A_T` + `chr13_100001406_G_A`, joint AF 9/8182) — HGDP sample
`LP6005443-DNA_E02` is phased **`0|1`/`0|1`** = both ALTs on one haplotype (cis), confirmed in the VCF;
PUTATIVE-only on the hybrid.

### 1b. Keeping singletons: why we re-phased instead of using the published phased release

A natural shortcut would be to use gnomAD's **already-phased** HGDP+1kGP release
(`phased_haplotypes_v2`). We could **not**: that release applies a post-phasing **MAC ≥ 2 filter**
(`remove_singletons.py`) that **removes every singleton** (AC = 1 — variants private to a single
haplotype). Its allele-count spectrum has **zero AC = 1**. SHAPEIT5 itself does *not* drop singletons
(`phase_rare` phases them); the drop is a separate release step the gnomAD paper's methods do not
mention. We reproduced the SHAPEIT5 pipeline **skipping that step**, recovering:

- **59.9 M** autosomal singletons + **1,597,916** chrX singletons (public release: **0**),
  at ≥99.978% genotype concordance and ≤0.58% switch error vs the public backbone.

**Why that matters for off-targets — measured in this very search:** of the 5,996,092 off-targets,
**1,090,787 (18.2%) are driven by a recovered singleton**, including **13,684 at CFD ≥ 0.1** and
**2,972 at CFD ≥ 0.2**. Every one of these would be **invisible** with the published singleton-dropped
release. Singletons create off-targets *private to one individual* — exactly the private risk a
per-sample / clinical off-target screen exists to catch. Using the "phased release from the original
paper" would silently discard ~1.09 M off-targets (thousands of them scoring in an actionable range).

**Net of Point 1:** re-phasing is not cosmetic. (a) It unlocks CONFIRMED-cis + named carriers across
HGDP (27% → 97% confirmed); (b) doing it *ourselves, singleton-inclusive* avoids losing ~18% of all
off-targets that the ready-made phased release would have dropped.

---

## POINT 2 — Mega over-calls putative off-targets (search 4,1,1)

Mega (5 sources: 1kGP + HGDP + gnomAD v4.1 + TOPMed + All-of-Us; **sites-only** — aggregate allele
frequencies, **no genotypes**) vs the phased HGDP+1kGP, identical search: **mm4 + 1 DNA + 1 RNA bulge**.

> Mega can appear **only** in this comparison, not in Point 1: with no sample genotypes there is
> nothing to phase or unphase.

Because mega aggregates a very large, multi-source variant set but **cannot resolve cis**, it "matches"
a guide against many *possible* variant combinations — far more than co-occur in any real haplotype.
Those surface as **PUTATIVE** calls, a fraction of which are genuine false positives (combinations no
individual carries). The phased, genotyped index confirms which combinations are actually real, giving
a higher-precision ranking.

Top-1000 ranked sites (by CFD):

| top-1000 | CONFIRMED/named-carrier | reference | PUTATIVE (cis-unconfirmed) |
|---|---|---|---|
| **Phased HGDP+1kGP** | **578** (named carriers) | 422 | **0** |
| **Mega (sites-only)** | 0 (cannot resolve) | 399 | **273** + 328 AF-only "observed" |

- **27.3% of mega's top-1000 are PUTATIVE** (cis-unconfirmed, false-positive–prone), vs **0%** for the
  phased index — every phased top site is either reference or a confirmed named-carrier off-target.
- Detection breadth: mega 184,449 vs phased 213,545 off-targets; top-1000 overlap 558/1000.

### Cross-match: how many of mega's top putatives are real? (an estimate)

The mega putative top sites are **CFD = 1.0 perfect-match off-targets fabricated by stacking many
rare ALT alleles into one 23 bp window** — a median of **14 variants** per site, up to **33** (min 2).
We looked each one up in the phased, genotyped panel (6,1,1 per-sample oracle, 3.62 M CONFIRMED-cis
loci over 4,091 individuals) to ask: does any real person actually carry this combination in cis?

Of mega's **228 distinct** putative top-1000 loci:

| outcome in the 4,091 phased individuals | loci | median variants/site |
|---|---|---|
| **CONFIRMED real** (a real carrier exists) — *mega recovers a true site* | **38 (17%)** | **2** |
| **Not confirmable** (189 absent from panel + 1 present-but-never-cis) | **190 (83%)** | **15** |

**Estimated false-positive rate ≈ 83%** of mega's top-ranked putatives. The split is the whole story:
the sites that are real are **simple** co-occurrences (median **2** variants); the ones that cannot be
confirmed are **implausible stacks** (median **15**, up to 33 variants in cis) that no haplotype
realistically carries.

> **This is an estimate, stated explicitly.** "Not confirmable" means *not observed in cis in our
> 4,091-sample phased panel* — 189 of the 190 involve variants absent from HGDP+1kGP (they come from
> mega's much larger aggregate sources: gnomAD v4.1, TOPMed, All-of-Us), so "absent here" is not
> *proof* of falsehood. But a 15-variant (up to 33) cis stack in 23 bp is combinatorially implausible
> at **any** cohort size, so the overwhelming majority are genuine false positives, not merely
> unobserved. The point: **mega recovers the real sites** (the 38 simple ones) **but pays for its
> breadth with an ~83% putative-false-positive load at the top of the ranking**, which the phased
> index does not.

### The flip side — real rare single-SNP sites the phased index confirms

The phased index's value is not only pruning mega's stacks; it also **surfaces real off-targets created
by a single rare SNP and names the carrier** — hits a reference-only search misses entirely and a
sites-only panel can only guess at. In the phased 6,1,1 search: **2,374,162** single-SNP
carrier-backed off-targets, of which **1,870,135 are rare** (MAF ≤ 0.001), **5,085 actionable**
(CFD ≥ 0.2), and **1,062,519 driven by a singleton** (invisible to the singleton-dropped release,
cf. Point 1b). These are unambiguously real (one variant — no cis question) and tied to a named
individual.

**Net of Point 2:** mega's breadth is a useful worst-case detection screen and it *does* recover the
real sites — but ~83% of its top-ranked putatives are combinatorial stacks no individual carries, so
its ranking is dominated by false positives at the top. The phased index removes that load (0% putative
in its top-1000) **and** adds ~1.87 M real, carrier-backed single-rare-SNP off-targets — a cleaner,
higher-precision, individually-actionable ranking.

---

## Take-home for the paper

1. **Re-phasing was necessary and had to be done in-house, singleton-inclusive.** It converts 22.5 M
   PUTATIVE → CONFIRMED cis with named HGDP carriers (27% → 97%), *and* recovers ~61.5 M singletons that
   drive 18.2% of all off-targets (thousands actionable) — all lost by the ready-made phased release.
2. **A sites-only mega panel over-calls putative off-targets.** 27% of mega's top-1000 are
   cis-unconfirmed (vs 0% for the phased index); cross-matched against real genotypes, **~83% of its
   top putatives cannot be confirmed in 4,091 individuals** — they are implausible stacks (median 15
   variants). Mega *does* recover the real sites (the simple, median-2-variant ones), but its ranking
   is top-loaded with false positives; the phased index gives a higher-precision, carrier-backed
   ranking and additionally confirms ~1.87 M real single-rare-SNP off-targets with named carriers.
3. The top single-site hits themselves are **stable** phased-vs-hybrid (793/1000 shared, 767 identical
   CFD), so making the phased index the default is safe.

## Methods / reproducibility

- **Re-phasing:** gnomAD HGDP+1kGP dense callset, SHAPEIT5 (`phase_common` scaffold + `phase_rare`,
  skipping the release's `remove_singletons.py` MAC≥2 step), chrX PAR/non-PAR correct ploidy.
- **Index:** `NRG_3_hg38+hg38_HGDP1kGP`, 4,091 samples (929 HGDP + 3,162 1kGP), 125.2M SNPs + 17.2M
  indels. HuggingFace `lucapinello/crisprme-data`; default in CRISPRme+ **v2.7.0**
  (`pinellolab/crisprme:v2.7.0`).
- **Searches:** Point 1 = 6,1,1 `--per-sample` (phased vs hybrid); Point 2 = 4,1,1 (phased vs mega).
  Different budgets by design — each conclusion is budget-independent.
- **Artifacts:** cluster `/srv/local/lp698_PAPER_phasing_importance/` (both full reports, top-1000
  tables, phase tally); build recipe `seq_script/merge_panels/hgdp1kgp_build.sh`; methods
  `docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`.
