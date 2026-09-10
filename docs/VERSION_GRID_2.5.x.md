# CRISPRme+ 2.5.x version grid

How the 2.5.x line evolved. **2.5.1 was a development-only line** (never tagged); its features
ship in the **2.5.2** release alongside the new sites-only co-occurrence work. So a user goes
from **2.5.0** (tagged) directly to **2.5.2** (tagged); the 2.5.1 column shows where each piece
first landed on `dev`.

## Capability grid

| Capability | 2.5.0 | 2.5.1 *(dev only)* | 2.5.2 |
|---|:---:|:---:|:---:|
| Dictionary-less variant search (Tier-0 registry + Tier-1 genotype store) | ✅ | ✅ | ✅ |
| Observed-haplotype enumeration (no phantom off-targets) | ✅ | ✅ | ✅ |
| **SNP + indel** co-occurrence, genotyped (CONFIRMED phased / PUTATIVE unphased, per-sample carriers + joint AF) | ✅ default-on | ✅ | ✅ |
| Genotyped index `NRG_3_hg38+hg38_1000G2021_HGDP` on HF | ✅ | ✅ | ✅ |
| COSMIC + ENCODE SCREEN v4 + GENCODE + DHS annotation | ✅ | ✅ | ✅ |
| Shareable HTML off-target report | ✅ | ✅ | ✅ |
| High-variant-density cap + `--max-total-edits` | ✅ | ✅ | ✅ |
| **Two-pass `--fast` mode** (worst-possible reps for dense panels; exact CFD, screen-grade CRISTA) | ❌ | ✅ opt-in | ✅ |
| **All-source "mega" index** (1000G-2021 + HGDP + gnomAD + TOPMed + AoU; sites-only, per-dataset AF + `AF_max`) | ❌ | ⚠️ SNP-only (indels not searchable) | ✅ **searchable indels genome-wide** |
| Per-dataset **indel AF** companion (`indel_af.tsv`) | ❌ | ✅ emitted | ✅ emitted + **bundled** into report.zip |
| **SNP + SNP** co-occurrence companion (`snp_snp_cooc.tsv`) | ❌ | ❌ | ✅ **new** |
| **SNP + indel PUTATIVE** on sites-only panels (no genotypes → min-AF joint bound) | ❌ | ❌ | ✅ **new** |
| **Lossless dense-region** worst-case haplotype on sites-only (`CRISPRME_LOSSLESS_DENSE`) | ❌ | ❌ | ✅ default-on (sites-only) |
| Three-rung confidence model (CONFIRMED / PUTATIVE co-carrier / PUTATIVE estimated) | partial | partial | ✅ full (both dimensions, both indices) |
| `radar_chart` tolerant of IUPAC in sites-only DNA (report never crashes) | ❌ | ❌ | ✅ fixed |

Legend: ✅ present · ⚠️ present-but-limited · ❌ absent.

## The two production indices (2.5.2)

| | `NRG_3_hg38+hg38_1000G2021_HGDP` | `NRG_3_hg38+hg38_mega` |
|---|---|---|
| Sources | 1000 Genomes 2021 + HGDP | 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us |
| Data model | **sample-level genotypes** | **sites-only** (aggregate allele frequency) |
| Haplotypes | **observed** (real individuals) | **putative** (co-located variants) |
| Co-occurrence confidence | CONFIRMED (phased) / PUTATIVE (unphased) with **named carriers** + exact joint AF | **PUTATIVE** with conservative **min-AF** bound (no carriers) |
| Provenance | per-panel AF | **per-dataset AF** (`AF_1000G2021 … AF_AoU`) + `AF_max` |
| SNP+SNP + SNP+indel | ✅ | ✅ |
| Panel AN | 8262 (4131 samples) | sites-only (no AN; AF from sources) |

## Empirical progression (chr22 development validation, guide CTAAC, NRG, mm4/b1/b1)

Measured while validating 2.5.2 (see the chr22 progression + the genome-wide mega verify):

| Output | Genotyped (observed) | Mega (sites-only) |
|---|---|---|
| SNP+SNP co-occurrence rows | 3,329 (35 CONFIRMED + 3,294 PUTATIVE) | 139 PUTATIVE |
| SNP+indel co-occurrence rows | 14 (2 CONFIRMED + 11 PUTATIVE) | 6,502 PUTATIVE |
| indel_af rows | — | 6,824 (per-dataset AF) |

Genome-wide mega verify (all 24 chr, guide CTAAC): **839,897 SNP+indel + 10,753 SNP+SNP
co-occurrences (all PUTATIVE), 868,724 indel_af rows, indel off-targets on all 24 chromosomes,
0 pseudo-sample leaks.** Provenance spread: 785K single-source → 8.7K reported by all 5 datasets.

## Genome-wide clean-room (2.5.2, new-user from scratch: apptainer-pull v2.5.2 + HF-download)

TRAC guide `CTCTCAGCTGGTACACGGCA`, NRG. Two configs on separate hosts (see
`RELEASE_REPORT_2.5.2.md` for the full tables):

- **Mega (sites-only) validated at both configs.** ml008 mm6/b2/b2: 318,455 SNP+indel + 3,453
  SNP+SNP (**all PUTATIVE**) + 329,980 indel_af, 550K off-targets. ml007 mm4/1/1: 38,222 SNP+indel
  + 327 SNP+SNP PUTATIVE + 39,358 indel_af. Genome-wide sites-only feature set holds, heavyweight.
- **Genotyped observed-haplotype path validated** (ml007 mm4/1/1, all 4 cells passed): **124
  CONFIRMED SNP+SNP + 18 CONFIRMED SNP+indel with named carriers**, 119K off-targets, report.zip
  — the phased-cis class the mega cannot assert.
- **slow vs `--fast` (genotyped):** slow = observed CONFIRMED haplotypes (with carriers); `--fast`
  = more worst-possible representative rows, PUTATIVE (no per-sample cis) — its documented tradeoff.
- **Resource note:** 4× genotyped GW at mm6/b2/b2 in parallel over-subscribed RAM on one host
  (one cell's post-analysis worker OOM-hung); run fewer genotyped GW searches concurrently, or
  lower `CRISPRME_POSTPROC_MAX_WORKERS`. Single-search / lighter-config genotyped completes fine.
- **New-user friction:** anonymous HF downloads hit 429 rate limits; an HF token (or a
  retry-on-429 in the downloader — candidate fix) avoids it.
