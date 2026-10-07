# Why re-phasing mattered, and phased vs. mega

Paper narrative / internal notes.

**Guide** `ACTGAAATCTGTAAGCAGGC` · SpCas9 (NRG PAM) · hg38 · CFD + CRISPR-Bulge scores · all
functional/cancer annotations.

## 1. Background & motivation

The previously shipped default (`NRG_3_hg38+hg38_1000G2021_HGDP`) is a **hybrid**: the 1000
Genomes portion is phased, but the **HGDP portion is unphased**. CRISPRme reports a co-occurring
off-target — two nearby variants that must sit together on one chromosome copy (**in cis**) for the
off-target to form — as **CONFIRMED** only when it can prove the variants are in cis in a real
individual, which **requires phased genotypes**. On the hybrid index, every HGDP-involving
co-occurrence therefore collapses to **PUTATIVE** ("possible, but cis unproven"), with no named HGDP
carrier. We re-phased the gnomAD HGDP+1kGP callset end-to-end with SHAPEIT5 so **both** cohorts are
phased, rebuilt the index (`NRG_3_hg38+hg38_HGDP1kGP`, 4,091 samples), and re-ran the identical search.

**There are two distinct comparisons:**

1. **The re-phasing comparison** — phased vs unphased (same cohorts; phasing is the variable).
2. **Mega vs the latest phased index** — broad sites-only vs focused phased.

> **Why mega is only in comparison 2:** mega is **sites-only** (aggregate allele frequencies, *no
> sample genotypes at all*), so there is nothing to phase or unphase — it cannot participate in the
> phased-vs-unphased comparison.

---

## 2. Comparison 1 — the re-phasing comparison (phased vs unphased) — **6,1,1**

Phased HGDP+1kGP vs the unphased hybrid, identical search: **mm6 + 1 DNA + 1 RNA bulge**, `--per-sample`.

### 2a. Headline — phasing's real value is CONFIRMED cis

| phase_confirmation | CONFIRMED | PUTATIVE | % CONFIRMED |
|---|---|---|---|
| Hybrid (old) | 8,393,667 | 22,515,383 | 27.2% |
| **Phased (new)** | **36,319,150** | **971,014** | **97.4%** |

PUTATIVE collapses **23×**, CONFIRMED rises **4.3×**. The residual ~0.97M PUTATIVE is *correct* —
genuinely unobserved cis combinations. Autosomes ~99.6–99.9% CONFIRMED; chrX improves 5× (106k → 526k).

### 2b. Do the top off-targets change? Mostly no — the switch is safe

Top-1000 single-site off-targets (by CFD): **793/1000 identical**; of the 793 shared, **767 (96.7%)
have identical CFD** (max Δ 0.22). The highest-risk, actionable predictions do **not** move.

### 2c. What *are* the 207 differences?

| | # | What they are |
|---|---|---|
| **Phased-only** top-1000 | 207 | **All variant-driven, all with named carriers** (CFD 0.42–0.75) — the richer QC'd gnomAD panel surfaces more real variant-driven off-targets, with their carriers, into the top ranks |
| **Hybrid-only** top-1000 | 207 | **75 reference** + **132 old-panel variant** (CFD 0.41–0.68). The 75 reference sites are **not lost** — detected identically in both, just re-ranked below #1000 by the phased index's extra high-scoring variant sites. The 132 variant sites are alleles specific to the raw 1000G+HGDP panel, not in gnomAD's QC'd callset |

So the differences are **panel-composition** effects (which variants each callset contains), **not
phasing** effects. Phasing does not add/remove sites or change CFD — it changes the *cis interpretation*.

### 2d. Ground truth — cis confirmed in the VCF

Off-target `chr13:100001388(−)` uses co-occurring SNPs `chr13_100001399_A_T` + `chr13_100001406_G_A`
(joint AF 9/8182 = 0.0011). HGDP sample `LP6005443-DNA_E02` is phased **`0|1` / `0|1`** — both ALT
alleles on the **same haplotype (cis)**, confirmed directly in the re-phased VCF. On the unphased
hybrid, PUTATIVE only.

---

## 3. Comparison 2 — mega vs the latest phased index — **4,1,1**

Mega (5 sources: 1000G + HGDP + gnomAD v4.1 + TOPMed + All-of-Us; **sites-only**) vs phased HGDP+1kGP,
identical search: **mm4 + 1 DNA + 1 RNA bulge**.

- **Detection breadth:** mega 184,449 vs phased 213,545 off-targets.
- **Top-1000 overlap:** 558/1000 sites shared.
- **Carrier resolution — the decisive difference:**

| top-1000 | named carriers | breakdown |
|---|---|---|
| **Phased** | **578** | 578 named-carrier + 422 reference (0 unresolved) |
| **Mega** | **0** | 399 reference + 328 observed (AF only) + 273 putative |

Mega is sites-only: it gives a worst-case allele-frequency bound but **cannot tell you who carries a
site or whether co-occurring variants are in cis**. Its top-1000 is a useful broad detection screen
but is **not actionable at the individual/population level**. The phased index gives named carriers +
confirmed cis for the majority of its top sites.

> **Edit-budget note.** Comparison 1 (re-phasing) was run at **6,1,1**; comparison 2 (mega) at
> **4,1,1**. If a single budget is preferred for the paper, the mega comparison can be re-run at
> **6,1,1** to match (the phased 6,1,1 search already exists). Mega is tractable at 6,1,1 (1-bulge);
> only the dense 2+2-bulge worst case is heavy.

---

## 4. Take-home for the paper

- Re-phasing's value is the **completeness/correctness of co-occurrence interpretation** (22.5M →
  0.97M PUTATIVE; named HGDP carriers; confirmed cis) — **not** a reshuffling of the top off-targets.
- The top off-targets are **stable** across the phased/unphased switch → making the phased index the
  default is **safe**.
- Phased ≫ sites-only mega for carrier-level prioritization; mega is a broad worst-case detection
  screen only.

## 5. Methods / reproducibility

- **Re-phasing:** gnomAD HGDP+1kGP callset, SHAPEIT5 (phase_common scaffold + phase_rare for
  singletons), chrX PAR/non-PAR correct ploidy. Validation: ≥99.978% GT concordance, switch error
  ≤0.58% vs the public backbone, 59.9M singletons recovered.
- **Index:** `NRG_3_hg38+hg38_HGDP1kGP`, 4,091 samples (929 HGDP + 3,162 1000 Genomes), 125.2M SNPs +
  17.2M indels. Published on HuggingFace `lucapinello/crisprme-data`; shipped as the default in
  CRISPRme+ **v2.7.0**.
- **Artifacts:** cluster `/srv/local/lp698_PAPER_phasing_importance/` (both full reports, top-1000
  tables, phase tally); build recipe `seq_script/merge_panels/hgdp1kgp_build.sh`; methods
  `docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`.
