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
keep their archives (index format 1) unchanged, so an older client never finds an index here.
CRISPRme+ 2.10.0 reads both formats; a download of an archive of 2.8.1 by its own name names the
index here that replaces it.

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
NRG on HPRC and mega is served by the NRR index there, pruned: the same rows, more slowly than
a dedicated NRG index (about twice the wall time in the engine's search of chr22's reference
sequence). HGDP + 1000 Genomes is genotyped and
phased (4,091 samples), HPRC genotyped and phased (232 assembly-derived genomes; HPRC Release 2's
Minigraph-Cactus graph v2.1, `hprc-v2.1-mc-grch38`, where 2.8.1's HPRC archives are of v2.0), mega
sites-only (5 sources, allele frequencies only).

**Contigs.** The reference indexes and those of HGDP + 1000 Genomes hold all 455 sequences of
hg38 (its alt, random and unplaced contigs and chrM included), so a variant search finds every
site a reference search finds; those of HPRC and mega hold its 24 chromosomes, chr1-22, X and Y.
An alt contig repeats a region of a primary chromosome (the HLA region, for one) and no VCF of
these panels covers it, so a site there can be reported twice, the second time as a reference
row on the alt contig.

**What these indexes hold.** Every allele of a split multi-allelic site is coded (the
four-allele symbol Q where the alleles cover A, C, G and T), each indel is searched in a
fake-indel window of 30 reference bases on either side (25 and 27 in 2.8.1's archives), and
an index and its Dictionaries hold the records whose FILTER is `PASS` or `.`, those CRISPRitz's
enricher reads. The mega is rebuilt from
the same merged VCFs: its registry and indel frequencies hold exactly the records its genome
codes, and each indel carries its `AF_max`.

**Searching them.** The default budget is 4 mismatches + 1 DNA + 1 RNA bulge. Above 1 + 1 a
search of a variant index here warns: the variant search of the densest windows (the
immunoglobulin loci) grows by orders of magnitude, to hours or out of memory; raise the
mismatches first (CRISPRme 2.10.1 is planned to search the samples' own haplotypes there, with
these indexes as they are). A search of mega withholds by default its windows of 2^20
sequences or more, listing them in the dense-region files (`--no-sites-only-density-guard`
analyses them), and ranks its combinations whose frequencies fail the safeguards' bound after
the others.

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
│   ├── hg38_HGDP1kGP.members.md5                      # the md5 of each file the kit holds, by its path in it
│   ├── hg38_HPRC.tar.gz (+ .SHA256SUMS, .members.md5)
│   ├── hg38_mega.tar.gz (+ .SHA256SUMS, .members.md5)
│   └── hg38_HPRC/                                     # HPRC's normalised VCFs, which its index was built from
│       ├── hprc.chr1.norm.vcf.gz (+ .tbi)             # one a chromosome: chr1-22, X and Y
│       ├── ...
│       └── SHA256SUMS                                 # one SHA-256 per VCF and .tbi
└── README.md
```

**The kits and their inputs.** A kit (`build-inputs/<dataset>.tar.gz`) holds the enriched
genome and the fake-indel genome of its panel, `Genomes/hg38+<dataset>/` and
`Genomes/hg38+<dataset>_INDELS/`, from which `build-index-only` builds an index of any PAM of the
panel without its VCFs (`crisprme.py download --what build-inputs --dataset <dataset>`). Its
member list, `<dataset>.members.md5`, gives the md5 of each file the kit holds: the tarball is
not byte-reproducible (tar records each file's time, owner and listing order), so its
`SHA256SUMS` says that a download is the published file, and the list what the kit holds, the
same for any packing of the same files. `build-inputs/hg38_HPRC/` holds the VCFs the HPRC index
reads: HPRC Release 2's Minigraph-Cactus graph v2.1 on GRCh38 (`hprc-v2.1-mc-grch38`, the VCF
decomposed with `vcfwave`), split per chromosome and normalised as the panel's recipe does it,
`hprc.<chr>.norm.vcf.gz` and its `.tbi` (`lucapinello/crisprme-data`'s `vcfs/hg38_HPRC` holds
the v2.0 files 2.8.1's HPRC archives were built from). HGDP + 1000 Genomes' VCFs are
`lucapinello/crisprme-data`'s `vcfs/hg38_HGDP1kGP`, as they are, and mega's index is rebuilt
from that repo's `vcfs/hg38_mega`.

**Each panel's Dictionaries once.** A variant archive holds the index, its `_INDELS` index and a
`manifest.json`, which names the panel's **Dictionaries bundle**
(`dictionaries-bundle_<dataset>`: the Tier-0 registry, the indel logs, the phased indel
genotypes, mega's indel allele frequencies, and the samplesID lists) and its **genotype store**
(`genotypes_<dataset>`) by name and fingerprint. `download` fetches each once per panel: a
second index of the panel finds them installed, the same, and fetches its own archive (and its
PAM's reference index where that is not installed: `20-NGG_3_hg38` for
`20-NGG_3_hg38+hg38_HGDP1kGP`), not the bundle and the store again.
`--no-genotypes` leaves the genotype store out. (The archives of 2.8.1 each carried their own
copy of the Dictionaries, 16.5 GB for every HGDP + 1000 Genomes archive.)

## Usage

```bash
# the genome, annotations, PAMs and sample lists, and the default reference index
# (20-NRR_3_hg38, from here)
crisprme.py download --what all --path .

# the recommended variant-aware default: phased HGDP + 1000 Genomes, SpCas9 NRR
crisprme.py download --what index --index-name 20-NRR_3_hg38+hg38_HGDP1kGP --path .

# another index of the same panel: its archive, and 20-NGG_3_hg38 if not installed (the
# bundle and store are installed already)
crisprme.py download --what index --index-name 20-NGG_3_hg38+hg38_HGDP1kGP --path .
```

Two installs of one panel must read one build of its Dictionaries: install an index of 2.8.1
(the older names) in another working directory. A download refuses to put an archive of
2.8.1's Dictionaries over a bundle of this repo, and says which installed indexes read
Dictionaries of another build where a bundle replaces them, and where the panel's enrichment in
place (`Genomes/hg38+<panel>`, which 2.8.1's build or kit left) is of another build: a search
with no bulges scans it in place of the index, and refuses it, so remove it, or put this repo's
kit of the panel in its place.

## Verifying

- An index or a companion: from inside a downloaded `indexes/`, `sha256sum -c --ignore-missing
  SHA256SUMS`.
- A kit: beside the downloaded tarball, `sha256sum -c <dataset>.SHA256SUMS`; once unpacked,
  `md5sum -c <dataset>.members.md5` in the folder it was unpacked in (the one that holds its
  `Genomes/`).
- HPRC's VCFs: from inside a downloaded `build-inputs/hg38_HPRC/`, `sha256sum -c --ignore-missing
  SHA256SUMS`.
