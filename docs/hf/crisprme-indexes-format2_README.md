---
license: other
license_name: mgb-open-access-license-1.0
license_link: https://github.com/pinellolab/crisprme-plus/blob/main/LICENSE
pretty_name: CRISPRme+ 2.10 indexes (index format 2)
tags:
  - crispr
  - off-target
  - genomics
  - crisprme
  - reference-data
viewer: false
---

# CRISPRme+ 2.10 indexes — index format 2

The precomputed indexes of **CRISPRme+ 2.10.0** and later, built by **CRISPRitz 2.10.0** in
**index format 2**, with each variant panel's Dictionaries and genotype store once, and the
enriched-genome kits they were built from. `crisprme.py download --what index --index-name
<name>` fetches an index from here, with what it needs, for any name of this repo.

## Why this is a separate repo

**These indexes are not readable by earlier releases.** An index of format 2 opens and ends
with a magic and records its build parameters in a header; CRISPRme+ 2.8.1 refuses one ("this
CRISPRme reads index format 1 only") and 2.8.0 aborts on one. The repos earlier releases read,
[`lucapinello/crisprme-data`](https://huggingface.co/datasets/lucapinello/crisprme-data) and
[`lucapinello/crisprme-data-general-index`](https://huggingface.co/datasets/lucapinello/crisprme-data-general-index),
keep their archives (index format 1) unchanged, so an older client never finds an index here
(decision C2). CRISPRme+ 2.10.0 reads both formats; a download of an archive of 2.8.1 by its
own name names the index here that replaces it.

**Versions:** CRISPRme+ and CRISPRitz share their major.minor: CRISPRme+ 2.10.z with CRISPRitz
2.10.z (Bioconda pins `crispritz >=2.10.0,<2.11`; the Docker image the exact CRISPRitz it was
tested with). Every install route of 2.10.0 installs the matching CRISPRitz.

## Names

An index is named for its guide length and the side of its PAM, in the order of the target on
the strand: `20-NRR_3_hg38` is a 20-nt guide with a 3′ PAM, `TTTN-23_3_hg38` a 23-nt guide
with a 5′ PAM; `+<dataset>` marks a variant index (`20-NRR_3_hg38+hg38_HGDP1kGP`). `3` is N, a
budget of up to 2 DNA and 2 RNA bulges, the web's ceiling. The 2.8.1 names (`NRR_3_hg38`,
...) say neither length nor side.

## Indexes

| Index | Guide | PAM | Reference | HGDP + 1000 Genomes | HPRC | mega |
|---|---|---|---|---|---|---|
| SpCas9 **NRR** (NAA, NAG, NGA, NGG; the default) | 20 nt | 3′ | `20-NRR_3_hg38` | `20-NRR_3_hg38+hg38_HGDP1kGP` | `20-NRR_3_hg38+hg38_HPRC` | `20-NRR_3_hg38+hg38_mega` |
| SpCas9 **NGG** (the fast canonical search) | 20 nt | 3′ | `20-NGG_3_hg38` | `20-NGG_3_hg38+hg38_HGDP1kGP` | | |
| SpCas9 **NRG** (NAG + NGG) | 20 nt | 3′ | `20-NRG_3_hg38` | `20-NRG_3_hg38+hg38_HGDP1kGP` | | |
| **General** (PAM-less, 8-nt flank: any 3′ PAM of up to 8 nt, pruned) | 20 nt (21 at 1 + 1) | 3′ | `20-NNNNNNNN_3_hg38` | `20-NNNNNNNN_3_hg38+hg38_HGDP1kGP` | | |
| SaCas9 **NNGRRN** (serves NNGRRT, pruned) | 22 nt | 3′ | `22-NNGRRN_3_hg38` | `22-NNGRRN_3_hg38+hg38_HGDP1kGP` | | |
| Cas12a **TTTN** (serves TTTV, pruned) | 20 nt | 5′ | `TTTN-20_3_hg38` | `TTTN-20_3_hg38+hg38_HGDP1kGP` | | |
| Cas12a **TTTN** | 23 nt | 5′ | `TTTN-23_3_hg38` | `TTTN-23_3_hg38+hg38_HGDP1kGP` | | |

An index serves a search of a PAM its own holds, pruned (an NRR index serves NGG, NAG, NGA and
NRG; NNGRRN serves NNGRRT; TTTN serves TTTV), and the general index any 3′ PAM of up to 8 nt.
NRG on HPRC and mega is served by the NRR index there. HGDP + 1000 Genomes is genotyped and
phased (4,091 samples), HPRC genotyped and phased (232 assembly-derived genomes), mega
sites-only (5 sources, allele frequencies only).

**`NNN_3` is not built again for 2.10.0.** The general index serves every search `NNN_3`
served, at the same speed (a 3-nt and an 8-nt PAM-less index cost the same: every position is
a site), with the same rows; 2.8.1's `NNN_3` archives stay in `lucapinello/crisprme-data`.

The list, with how each archive is built and published, is the index set,
[`seq_script/rebuild/index_set_2.10.0.tsv`](https://github.com/pinellolab/crisprme-plus/blob/main/seq_script/rebuild/index_set_2.10.0.tsv).

## Layout

```
crisprme-indexes-format2/
├── indexes/
│   ├── 20-NRR_3_hg38.tar.gz                           # a reference index: the index alone
│   ├── 20-NRR_3_hg38+hg38_HGDP1kGP.tar.gz             # a variant index: it and its _INDELS index
│   ├── ...                                            # (the indexes above)
│   ├── dictionaries-bundle_hg38_HGDP1kGP.tar.gz       # the panel's Dictionaries, once
│   ├── dictionaries-bundle_hg38_HGDP1kGP.json         # what the bundle holds (its fingerprint)
│   ├── genotypes_hg38_HGDP1kGP.tar.gz (+ .json)       # the panel's Tier-1 genotype store, once
│   ├── indel_genotypes_hg38_HGDP1kGP.tar.gz (+ .json) # its phased indel genotypes, for the genotype oracle
│   ├── ... (the same for hg38_HPRC; hg38_mega has a bundle, no genotype store)
│   └── SHA256SUMS                                     # one SHA-256 per .tar.gz
├── build-inputs/
│   ├── hg38_HGDP1kGP.tar.gz (+ .SHA256SUMS)           # the enriched-genome kits these indexes were built from
│   ├── hg38_HPRC.tar.gz (+ .SHA256SUMS)
│   └── hg38_mega.tar.gz (+ .SHA256SUMS)
└── README.md
```

**Each panel's Dictionaries once.** A variant archive holds the index, its `_INDELS` index and a
`manifest.json`, which names the panel's **Dictionaries bundle**
(`dictionaries-bundle_<dataset>`: the Tier-0 registry, the indel logs, the phased indel
genotypes, mega's indel allele frequencies, and the samplesID lists) and its **genotype store**
(`genotypes_<dataset>`) by name and fingerprint. `download` fetches each once per panel: a
second index of the panel finds them installed, the same, and fetches its own archive alone.
`--no-genotypes` leaves the genotype store out. (The archives of 2.8.1 each carried their own
copy of the Dictionaries, 16.5 GB for every HGDP + 1000 Genomes archive.)

## Usage

```bash
# the genome, annotations, PAMs and sample lists, and the default reference index
# (20-NRR_3_hg38, from here)
crisprme.py download --what all --path .

# the recommended variant-aware default: phased HGDP + 1000 Genomes, SpCas9 NRR
crisprme.py download --what index --index-name 20-NRR_3_hg38+hg38_HGDP1kGP --path .

# another index of the same panel: its archive alone (the bundle and store are installed)
crisprme.py download --what index --index-name 20-NGG_3_hg38+hg38_HGDP1kGP --path .
```

Two installs of one panel must read one build of its Dictionaries: install an index of 2.8.1
(the older names) in another working directory. A download refuses to put an archive of
2.8.1's Dictionaries over a bundle of this repo, and says which installed indexes read
Dictionaries of another build where a bundle replaces them.

## Verifying

From inside a downloaded `indexes/`: `sha256sum -c --ignore-missing SHA256SUMS`.
