---
license: other
license_name: mgb-open-access-license-1.0
license_link: https://github.com/pinellolab/crisprme-plus/blob/main/LICENSE
pretty_name: CRISPRme+ reference data & precomputed indexes
tags:
  - crispr
  - off-target
  - genomics
  - crisprme
  - reference-data
viewer: false
---

# CRISPRme+ reference data & precomputed indexes

A fast, reliable mirror of the reference resources CRISPRme+ needs, plus
**precomputed indexes** — so users can skip the slow raw-data download, the
multi-hour variant **enrichment**, and the index build entirely.

Pulled automatically by `crisprme.py download` (and the web **Settings → Data
Manager**). See the CRISPRme+ docs for the full data layout and quickstart:
https://github.com/pinellolab/crisprme-plus

## Why this exists

A from-scratch, genome-wide CRISPRme+ run is dominated by **data movement +
enrichment**, not the search itself:

- Downloading the 1000 Genomes VCFs from EBI ran at **~1.4 MB/s (~3 h)**; the UCSC
  genome download is similarly throttled. **Hugging Face measured ~6.5× faster**
  from the same machine.
- The variant **enrichment** stage was **~93 % of the ~13.5 h wall-clock**
  (single-threaded per chromosome).

Shipping **precomputed indexes** removes both the download and the
enrichment/index-build for the common case, taking a genome-wide setup from many
hours to minutes.

## Layout

```
crisprme-data/
├── indexes/                        # PRECOMPUTED search indexes (genome_library tarballs)
│   ├── NRG_3_hg38.tar.gz                       # reference-only (no variants)
│   ├── NRG_3_hg38+hg38_HGDP1kGP.tar.gz         # ← recommended variant DEFAULT
│   ├── NRG_3_hg38+hg38_HPRC.tar.gz             # HPRC pangenome panel
│   ├── NRG_3_hg38+hg38_mega.tar.gz             # sites-only aggregate (worst-case fallback)
│   └── genotypes_hg38_{HGDP1kGP,HPRC}.tar.gz   # per-sample genotype companions
├── build-inputs/                   # PRE-ENRICHED kits (fast rebuild / new PAM: download --what build-inputs)
│   ├── hg38_HGDP1kGP.tar.gz (+ .SHA256SUMS)
│   └── hg38_mega.tar.gz     (+ .SHA256SUMS)
├── vcfs/                           # SOURCE VCFs — build any index from scratch / with a new PAM
│   └── hg38_{HGDP1kGP,HPRC,mega}/         # download --what vcf --dataset <name>
├── genomes/hg38/                   # per-chromosome FASTA mirror (chr1.fa … chrX.fa)
├── annotations/                    # functional + cancer-driver tracks (see below)
├── pams/                           # 20bp-NRG-SpCas9, 20bp-NGG-SpCas9, 20bp-NNN-NO-PAM
└── samplesIDs/                     # per-dataset sample-id lists
```

## Precomputed indexes

All indexes are **NRG** PAM (the SpCas9 `NGG`/`NGA` superset, the conservative
choice for off-target screening), scored with **CFD + CRISPR-Bulge**.

| Index | Variants | Phasing | Notes |
|---|---|---|---|
| `NRG_3_hg38` | none | — | reference-only |
| **`NRG_3_hg38+hg38_HGDP1kGP`** | HGDP + 1000 Genomes (4,091 samples = 925 HGDP + 3,166 1000G) | **genotyped-phased** (SHAPEIT5 re-phased) | **recommended variant default** — CONFIRMED *cis* co-occurrence + named carriers |
| `NRG_3_hg38+hg38_HPRC` | HPRC pangenome (232 samples incl. CHM13) | genotyped-phased | pangenome panel |
| `NRG_3_hg38+hg38_mega` | aggregate (1000G + HGDP + gnomAD + TOPMed + All of Us) | sites-only | worst-case fallback; a listed haplotype may not exist in any individual |

The two genotyped panels ship a `genotypes_hg38_*` companion that enables
`--per-sample` resolution (CONFIRMED *cis* + exact joint AF + named carriers).

**Annotations:** `dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz` (functional +
cancer, with COSMIC gated behind a licence attestation),
`intogen_drivers.hg38.bed.gz` (IntOGen cancer drivers, CC0, default-on), and
`gencode.protein_coding.bed.gz` (closest-gene).

## Usage

The indexes are fetched automatically:

```bash
# one-shot: reference + recommended variant default + annotations + PAMs
crisprme.py download --what all --path .

# or a single index
crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP --path .
```

Direct download also works:

```bash
hf download lucapinello/crisprme-data --repo-type dataset --local-dir crisprme-data
```

## License & citation

Reference data are redistributed under the terms of their original sources
(1000 Genomes, HGDP, HPRC, gnomAD, GENCODE/ENCODE, IntOGen, UCSC hg38, etc.).
**CRISPRme+ is licensed under the MGB Open Access License 1.0 (MGBOA 1.0)** — see
https://github.com/pinellolab/crisprme-plus/blob/main/LICENSE.

If you use CRISPRme, please cite:
Cancellieri S, *et al.* Human genetic diversity alters off-target outcomes of
therapeutic gene editing. *Nat Genet* 55, 34–43 (2023).
doi:10.1038/s41588-022-01257-y

Project: https://github.com/pinellolab/crisprme-plus
