# Genome-wide SNP + indel co-occurrence summary (2.5.1)

Summary of the SNP+indel *cis* co-occurrence output of the CRISPRme+ 2.5.x indel+SNP
feature, from the complete genome-wide feature-on run (`v2_dev_cooc`). This is the
canonical GW 2.5.1 cooc reference; it supersedes the interim `gw_cooc_25` recovery sweep
(which was reuse-based and never completed chr1).

## Run
| Field | Value |
|---|---|
| Guide | `TGCTTGGTCGGCACTGATAG` (NRG, SpCas9) |
| Index | `NRG_3_hg38+hg38_1000G2021_HGDP` (1000G 2021 @30× + HGDP, feature-on) |
| Search | mm 5, bDNA 2, bRNA 2, `--max-total-edits 6` (from-HF-index, dev 2.5.x) |
| Total off-targets | 407,222 |
| Report | `v2_dev_cooc_report.zip` — bundles `data/indel_snp_cooc.tsv` + a "SNP + indel cis co-occurrences" section |

> Note: at `--max-total-edits 6` a 5mm+2-bulge indel alignment (7 edits) is pruned; raise
> the cap to keep the deepest indel co-occurrences. This run captures every co-occurrence
> within the 6-edit budget.

## Co-occurrence results
- **2,729 co-occurrence rows** across **2,261 distinct off-target loci** (the extra rows are
  bulge-alignment variants of the same locus, ±1–2 bp).
- **Phase split:** **843 CONFIRMED** (cis proven by phasing) / **1,886 PUTATIVE**.
  - **All 843 CONFIRMED contain only 1000G-2021 (phased) samples** — 0 HGDP, as expected
    (HGDP is genotyped-**unphased**, so it can never be CONFIRMED-cis).
  - Of the 1,886 PUTATIVE: **1,449 are pure-HGDP** (unphased) and **136 carry a phased 1000G
    sample** but span a phase-uncertain boundary (so they stay PUTATIVE despite a phased panel).
- **Joint allele-frequency distribution** (the frequency of the indel+SNP haplotype):

  | joint AF | rows |
  |---|---|
  | ≥ 0.5 (common / near-fixed) | 20 |
  | 0.05 – 0.5 | 64 |
  | 0.005 – 0.05 | 179 |
  | < 0.005 (rare) | 2,466 (90%) |

  So the vast majority of co-occurrences are **rare, individual-specific** haplotypes; a
  small tail are common population haplotypes.
- **Highest-AF co-occurrences** (all PUTATIVE — mixed 1000G+HGDP panels, ~4,000 carriers):
  - `chr8:14511210` `G>GT` + SNP `rs55905580` — joint AF **0.884**
  - `chr8:9458936` `A>ATG` + SNP `rs56350029` — joint AF 0.865
  - `chr10:50227336` `T>TGG` + SNP `rs1192968381` — joint AF 0.818
  - `chr10:53124384` `TG>T` + SNP `rs35603270` — joint AF 0.736
- Highest-AF **CONFIRMED** (phased 1000G): `chrX:621958` `A>AC` + rs77899810 — joint AF 0.428.

## Interpretation
- The co-occurrence feature identifies off-targets that require an indel **and** a SNP on the
  same haplotype — a class 2.4.0 cannot represent. On this guide it accounts for ~1.2% of the
  105,471 new off-targets on the 24 shared primary contigs (≈1,262 cooc-proximal loci, of
  which only 14 are net-new beyond the panel expansion); the overwhelming majority (~96.3%)
  of the 2.5.x-vs-2.4.0 gain is the **denser 2021 panel** (variant-created off-targets), not
  the cooc feature. See `docs/RELEASE_CHECKLIST_2.5.1.md` for the confound-free delta.
- **Frequency ≠ cut strength.** The high-AF co-occurrences above are *common variants that
  create a weak off-target*; most cooc loci sit at the low-CFD margin. The clinically
  actionable subset (high CFD **and** co-occurring) is small — the report's validation panel
  and CFD tiers flag these; the cooc TSV is the exhaustive per-sample record.
- **`--fast` does not change this output** — the `indel_snp_cooc.tsv` is byte-identical
  between `--fast` and non-`--fast` runs (the cis-phasing pass is independent of the SNP
  worst-possible collapse). See METHODS §8.
