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
│   ├── NRR_3_hg38.tar.gz                       # reference-only (no variants), NRR: the default PAM since CRISPRme+ 2.10.0
│   ├── NRR_3_hg38+hg38_HGDP1kGP.tar.gz         # ← recommended variant DEFAULT
│   ├── NRR_3_hg38+hg38_HPRC.tar.gz             # HPRC pangenome panel
│   ├── NRR_3_hg38+hg38_mega.tar.gz             # sites-only aggregate (worst-case fallback)
│   ├── NRG_3_hg38{,+hg38_HGDP1kGP,+hg38_HPRC,+hg38_mega}.tar.gz   # NRG PAM (the default before 2.10.0)
│   ├── NGG_3_hg38{,+hg38_HGDP1kGP}.tar.gz      # NGG PAM (cannot serve NAG off-targets)
│   ├── NNN_3_hg38{,+hg38_HGDP1kGP}.tar.gz      # PAM-LESS 3-nt flank — serves ANY 3-nt 3′ PAM
│   ├── genotypes_hg38_{HGDP1kGP,HPRC}.tar.gz   # per-sample genotype companions
│   └── SHA256SUMS                              # checksum per archive (sha256sum -c --ignore-missing, inside indexes/)
├── build-inputs/                   # PRE-ENRICHED kits (fast rebuild / new PAM: download --what build-inputs)
│   ├── hg38_HGDP1kGP.tar.gz (+ .SHA256SUMS)
│   └── hg38_mega.tar.gz     (+ .SHA256SUMS)
├── vcfs/                           # SOURCE VCFs — build any index from scratch / with a new PAM
│   └── hg38_{HGDP1kGP,HPRC,mega}/         # download --what vcf --dataset <name>
├── genomes/hg38/                   # hg38 FASTA mirror, one file per sequence (download --what genome unpacks it)
│   ├── chr1.fa.gz … chrY.fa.gz, chrM.fa.gz     # gzipped; the alt, random and unplaced contigs too (455 files)
│   └── chr22.fa                                # chr22 also unpacked: complete-test and setup fetch it as is
├── annotations/                    # functional + cancer-driver tracks (see below)
├── pams/                           # 20bp-NRG-SpCas9, 20bp-NGG-SpCas9, 20bp-NNN-NO-PAM (CRISPRme ships 20bp-NRR-SpCas9)
├── samplesIDs/                     # per-dataset sample-id lists
├── guides/                         # sg1617.txt: the example guide CTAACAGTTGCTTTTATCACNNN (download does not fetch it)
└── README.md
```

## Precomputed indexes

**NRR** is the recommended PAM since CRISPRme+ 2.10.0 (the SpCas9 `NGG`, `NAG`, `NGA` and
`NAA` superset — the conservative choice for off-target screening; it was **NRG**, `NGG` +
`NAG`, before); NRR and NRG are published on all three variant panels, NGG and a PAM-less
index on HGDP + 1000 Genomes only. CRISPRme+ 2.10.0 and later search an NRR index for NGG,
NAG, NRG or NGA, and an NRG one for NGG or NAG, pruning it to that PAM by itself. All are scored with
**CFD + CRISPR-Bulge**. Every index here stores a **3-nt PAM region**, which is
what any CRISPRme release reads for a 3-nt PAM, so any CRISPRme release can read
them; use the rebuilt `NNN_3` archives with CRISPRme+ 2.8.1 or later (see below).

### Variant panels (the axis that matters most)

| Index | Variants | Phasing | Notes |
|---|---|---|---|
| `NRR_3_hg38` | none | — | reference-only |
| **`NRR_3_hg38+hg38_HGDP1kGP`** | HGDP + 1000 Genomes (4,091 samples = 925 HGDP + 3,166 1000G) | **genotyped-phased** (SHAPEIT5 re-phased) | **recommended variant default** — CONFIRMED *cis* co-occurrence + named carriers |
| `NRR_3_hg38+hg38_HPRC` | HPRC pangenome (232 samples incl. CHM13) | genotyped-phased | pangenome panel |
| `NRR_3_hg38+hg38_mega` | aggregate (1000G + HGDP + gnomAD + TOPMed + All of Us) | sites-only | worst-case fallback; a listed haplotype may not exist in any individual |

The two genotyped panels ship a `genotypes_hg38_*` companion that enables
`--per-sample` resolution (CONFIRMED *cis* + exact joint AF + named carriers).

**The mega for CRISPRme+ 2.10.0** is rebuilt from the same merged VCFs (`vcfs/hg38_mega`) by
`seq_script/merge_panels/mega_rebuild.sh`, with CRISPRitz 2.10.0 (index format 2, every allele
of a split site coded): its registry and its indel frequencies hold exactly the records its
genome codes, those whose FILTER is `PASS` or `.` (the archives above, of 2.8.1, list every
record of the VCF, 2,442,989 SNVs outside `PASS` among them), one allele once (of an allele on
two records, the more frequent record's: 8 long deletions genome-wide), and its indel log carries
each indel's `AF_max`, so an indel row shows its frequency. Searched with CRISPRme+ 2.10.0, a
sites-only dataset's candidates in windows of 2^20 sequences or more (the immunoglobulin loci
above all) are withheld and listed in the dense-region files, and noted on the rows there
(`--no-sites-only-density-guard` analyses them all), and its
rows of two or more variants whose frequencies fail the bound (product under 1e-06, or the
rarest under 0.001) are ranked after the others; with the guard a search of it took 1.37 (chr22;
1.46 in an earlier measurement) and 1.57 (chr14) times the shipped mega's, at 4 mismatches and
1 + 1 bulges. The rebuilt
archives are published with the release, in the location of the 2.10.0 indexes; the ones above
stay for earlier releases.

### Other PAMs

| PAM | Indexes | When to use |
|---|---|---|
| `NRG` | `NRG_3_hg38`, `+hg38_HGDP1kGP`, `+hg38_HPRC`, `+hg38_mega` | the default before CRISPRme+ 2.10.0 (`NAG`/`NGG`); a narrower, faster screen than NRR |
| `NGG` | `NGG_3_hg38`, `+hg38_HGDP1kGP` | canonical SpCas9 PAM only. **NGG cannot serve NAG or NGA off-targets**, so NRR is preferred for a safety screen |
| `NNN` (PAM-less) | `NNN_3_hg38`, `+hg38_HGDP1kGP` | one index for ANY 3-nt 3′ PAM, selected at search time with `--prune-pam`. Larger (71 GB / 110 GB downloads, GB = 10⁹ bytes) but removes the per-PAM rebuild |

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
the genomic-`N` flank fix in CRISPRitz 2.9.1, and replaced here on 2026-10-30. The
archives uploaded on 2026-10-25 lack every target whose 3-nt flank holds an assembly
`N`. CRISPRme+ 2.8.1 never reports a target with an `N` inside its PAM, so a search
with a **3-nt PAM** (`NGG`, `NRG`, `NRR`, or the unconstrained `20bp-NNN-NO-PAM`) gives
the same report on the old and the rebuilt archives (for NRG verified genome-wide,
per-sample, byte-for-byte against `NRG_3_hg38+hg38_HGDP1kGP`). A **shorter 3′ PAM**
searched with `--prune-pam`, such as `20bp-NG-SpCas9NG`, does change: the old archives
miss each site whose base just past the PAM is an `N` (at an assembly-gap edge, and in
the indel index, whose windows are separated by `N`s), and the rebuilt ones report it,
as a dedicated index for that PAM does. If you downloaded either archive before the
replacement and search it that way, remove `genome_library/NNN_3_hg38`,
`genome_library/NNN_3_hg38+hg38_HGDP1kGP` and
`genome_library/NNN_3_hg38+hg38_HGDP1kGP_INDELS`, then download again.
Use the rebuilt archives with CRISPRme+ 2.8.1 or later. 2.8.0 keeps the targets with
an `N` in their PAM that the rebuild adds, and a search that reaches one fails in
post-analysis or, on the CRISPRitz 2.9.0 that the 2.8.0 image builds and with RNA
bulges only, reports it as a 0-mismatch, 0-bulge off-target whose sequence stops at the
`N`. If you searched them with 2.8.0, run that search again with 2.8.1.

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
# one-shot: genome + annotations + PAMs + sample lists + the default reference index
# (NRR_3_hg38 since CRISPRme+ 2.10.0, NRG_3_hg38 before)
crisprme.py download --what all --path .

# then the recommended variant default (or any other index by name)
crisprme.py download --what index --index-name NRR_3_hg38+hg38_HGDP1kGP --path .
```

Direct download of single files also works (`crisprme.py download` additionally unpacks and
installs them). Name the files: the whole repo is about 432 GB.

```bash
hf download lucapinello/crisprme-data indexes/NRR_3_hg38+hg38_HGDP1kGP.tar.gz \
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
