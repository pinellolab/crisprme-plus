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
│   ├── NRR_3_hg38{,+hg38_HGDP1kGP,+hg38_HPRC,+hg38_mega}.tar.gz   # NRR PAM
│   ├── NGG_3_hg38{,+hg38_HGDP1kGP}.tar.gz      # NGG PAM (cannot serve NAG off-targets)
│   ├── NNN_3_hg38{,+hg38_HGDP1kGP}.tar.gz      # PAM-LESS 3-nt flank — serves ANY 3-nt PAM
│   ├── genotypes_hg38_{HGDP1kGP,HPRC}.tar.gz   # per-sample genotype companions
│   └── SHA256SUMS                              # checksum per archive (sha256sum -c, inside indexes/)
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

**NRG** is the recommended PAM (the SpCas9 `NGG` + `NAG` superset — the conservative
choice for off-target screening); NRG and NRR are published on all three variant
panels, NGG and a PAM-less index on HGDP + 1000 Genomes only. All are scored with
**CFD + CRISPR-Bulge**. Every index here stores a **3-nt PAM region**, which is
what any CRISPRme release reads for a 3-nt PAM, so all of them work with older
clients.

### Variant panels (the axis that matters most)

| Index | Variants | Phasing | Notes |
|---|---|---|---|
| `NRG_3_hg38` | none | — | reference-only |
| **`NRG_3_hg38+hg38_HGDP1kGP`** | HGDP + 1000 Genomes (4,091 samples = 925 HGDP + 3,166 1000G) | **genotyped-phased** (SHAPEIT5 re-phased) | **recommended variant default** — CONFIRMED *cis* co-occurrence + named carriers |
| `NRG_3_hg38+hg38_HPRC` | HPRC pangenome (232 samples incl. CHM13) | genotyped-phased | pangenome panel |
| `NRG_3_hg38+hg38_mega` | aggregate (1000G + HGDP + gnomAD + TOPMed + All of Us) | sites-only | worst-case fallback; a listed haplotype may not exist in any individual |

The two genotyped panels ship a `genotypes_hg38_*` companion that enables
`--per-sample` resolution (CONFIRMED *cis* + exact joint AF + named carriers).

### Other PAMs

| PAM | Indexes | When to use |
|---|---|---|
| `NRR` | `NRR_3_hg38`, `+hg38_HGDP1kGP`, `+hg38_HPRC`, `+hg38_mega` | broader than NRG (`NAA`/`NAG`/`NGA`/`NGG`); a wider, more conservative screen |
| `NGG` | `NGG_3_hg38`, `+hg38_HGDP1kGP` | canonical SpCas9 PAM only. **NGG cannot serve NAG off-targets**, so NRG is preferred for a safety screen |
| `NNN` (PAM-less) | `NNN_3_hg38`, `+hg38_HGDP1kGP` | one index for ANY 3-nt PAM, selected at search time with `--prune-pam`. Larger (71 GB / 110 GB downloads, GB = 10⁹ bytes) but removes the per-PAM rebuild |

**`NRR_3_hg38+hg38_mega` was replaced for CRISPRme+ 2.8.1.** The archive uploaded on
2026-10-27 shipped an empty indel index and empty indel logs, so no indel off-target
was ever reported from it, and installing it next to `NRG_3_hg38+hg38_mega` emptied
that index's indel logs as well (an install replaces each `Dictionaries/` folder).
If you downloaded it before the replacement, remove
`genome_library/NRR_3_hg38+hg38_mega` and `genome_library/NRR_3_hg38+hg38_mega_INDELS`
and download it again; that also restores NRG mega's logs. Both mega archives now
carry identical `Dictionaries/`.

**`NRR_3_hg38+hg38_HPRC` was replaced too**, for a smaller reason: its registry lacked
the stamp that makes a search on a phased panel default to per-sample analysis, so it
ran population-level unless `--per-sample` was passed (and, installed next to
`NRG_3_hg38+hg38_HPRC`, did the same to that one). Its index is unchanged; re-download
it the same way to get the per-sample default back.

**`NNN_3_hg38` and `NNN_3_hg38+hg38_HGDP1kGP` were rebuilt for CRISPRme+ 2.8.1**, with
the genomic-`N` flank fix in CRISPRitz 2.9.1. For a **constrained** search PAM (`NGG`,
`NRG`, `NRR`, i.e. every intended use of a PAM-less index) the archives uploaded on
2026-10-25 already gave the same results as the rebuilt ones and as the dedicated PAM
index (verified genome-wide, per-sample, byte-for-byte against
`NRG_3_hg38+hg38_HGDP1kGP`), so those searches do not change. The rebuild matters only
when the *search* PAM is the fully unconstrained `20bp-NNN-NO-PAM`: the old archives
under-reported sites whose PAM bases are assembly `N`. If you search that way, remove
`genome_library/NNN_3_hg38`, `genome_library/NNN_3_hg38+hg38_HGDP1kGP` and
`genome_library/NNN_3_hg38+hg38_HGDP1kGP_INDELS`, then download again.

**Wider PAM-less flanks** (flank > 3, e.g. `NNNNNNNN_3`) are **not** readable by a
CRISPRme release that predates `--prune-pam`, so they are published separately in
[`lucapinello/crisprme-data-general-index`](https://huggingface.co/datasets/lucapinello/crisprme-data-general-index)
and never land in this repo.

**Annotations:** `dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz` (functional +
cancer, with COSMIC gated behind a licence attestation),
`intogen_drivers.hg38.bed.gz` (IntOGen cancer drivers, CC0, default-on), and
`gencode.protein_coding.bed.gz` (closest-gene).

## Usage

The indexes are fetched automatically:

```bash
# one-shot: genome + annotations + PAMs + sample lists + the NRG_3_hg38 reference index
crisprme.py download --what all --path .

# then the recommended variant default (or any other index by name)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP --path .
```

Direct download of single files also works (`crisprme.py download` additionally unpacks and
installs them). Name the files: the whole repo is about 432 GB.

```bash
hf download lucapinello/crisprme-data indexes/NRG_3_hg38+hg38_HGDP1kGP.tar.gz \
  indexes/SHA256SUMS --repo-type dataset --local-dir crisprme-data
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
