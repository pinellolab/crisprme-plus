---
license: other
license_name: mgb-open-access-license-1.0
license_link: https://github.com/pinellolab/crisprme-plus/blob/main/LICENSE
pretty_name: CRISPRme+ general (wide-flank PAM-less) indexes
tags:
  - crispr
  - off-target
  - genomics
  - crisprme
  - reference-data
viewer: false
---

# CRISPRme+ general indexes — wide-flank, PAM-less

One index that serves **many PAMs**. A *general* index stores the raw genomic
flank instead of a fixed PAM, so a single build answers searches for any 3′ PAM no
longer than that flank — `NGG`, `NRG`, `NRR`, `NGA`, `NAG`, `NNGRRT` and so on —
with the PAM applied at search time by `--prune-pam` (3′ PAMs only: see the
orientation limit below).

Pulled by `crisprme.py download --what index --index-name <name>`, which routes
to this repo automatically for any flank wider than 3.

## Why this is a separate repo

**These indexes are not backward compatible.** A CRISPRme release that predates
`--prune-pam` reads exactly as many stored PAM bytes per leaf as the *search* PAM
is long. Handed a wider stored flank it reads the leaf stream at the wrong width,
desynchronises at the first leaf, and crashes — and because a crashed search used
to be reported as a clean empty result, that failure was silent.

CRISPRme+ **2.8.1** (CRISPRitz **2.9.1**) refuses such a mismatch up front, naming
both widths, and resolves the width correctly under `--prune-pam`. Older clients
cannot, so these archives are kept out of
[`lucapinello/crisprme-data`](https://huggingface.co/datasets/lucapinello/crisprme-data),
which older defaults point at. The split is derived from the flank width at
publish time rather than kept as a list of names, so it cannot be forgotten.

**Minimum versions: CRISPRme+ 2.8.1 / CRISPRitz 2.9.1.**

## Layout

```
crisprme-data-general-index/
├── indexes/
│   ├── NNNNNNNN_3_hg38.tar.gz                   # reference-only, 8-nt flank
│   ├── NNNNNNNN_3_hg38+hg38_HGDP1kGP.tar.gz     # HGDP + 1000 Genomes (4,091 samples, phased)
│   ├── genotypes_hg38_HGDP1kGP.tar.gz           # per-sample genotype companion of the variant index
│   └── SHA256SUMS                               # checksum per archive
└── README.md
```

## Indexes

| Index | Flank | Serves | Variants |
|---|---|---|---|
| `NNNNNNNN_3_hg38` | 8 nt | any 3′ PAM of length ≤ 8 | none (reference-only) |
| `NNNNNNNN_3_hg38+hg38_HGDP1kGP` | 8 nt | any 3′ PAM of length ≤ 8 | HGDP + 1000 Genomes, genotyped-phased (SHAPEIT5 re-phased) |

The variant archive bundles its `_INDELS` companion. Its per-sample genotype
store, `genotypes_hg38_HGDP1kGP.tar.gz`, is published here as well, byte-identical
to the copy in the main repo: `download` fetches the genotype companion and the
reference index from the repo the index itself came from, so without a copy here
a per-sample search on this index would run with degraded Samples.

**Orientation limit.** These indexes serve a **3′** PAM only. They store the flank
on the 3′ side of the guide, so a 5′ PAM (Cas12a-style, such as `TTTV`) is not served
at any width: use a dedicated index, which `crisprme.py build-index-only` builds.

## Usage

```bash
# the genome, annotations, PAMs and sample lists (plus the NRG_3_hg38 reference index)
crisprme.py download --what all --path .

# fetch the general index (routes here automatically on flank width)
crisprme.py download --what index --index-name NNNNNNNN_3_hg38+hg38_HGDP1kGP --path .

# then search with whatever PAM you need -- no rebuild
crisprme.py complete-search --genome Genomes/hg38 \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --guide guide.txt --pam PAMs/20bp-NRG-SpCas9.txt \
  --mm 4 --bDNA 1 --bRNA 1 --prune-pam --index-path genome_library \
  --annotation Annotations/dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz \
  --output my_search
```

Run the search from the folder you downloaded into: `download` writes
`list_vcf.txt` and `list_samplesID.txt` there, and `--prune-pam` needs
`--index-path`. `guide.txt` holds one spacer per line, padded with an `N` for each
PAM position (`NNN` for a 3-nt PAM).

`--prune-pam` is **CLI-only**: the web form does not pass the flag, so a general
index cannot be used from the web interface.

## Equivalence claim

A PAM-less index of flank *K* searched with `--prune-pam` and a 3′ PAM of length
*L ≤ K* produces the **same report as a dedicated index built for that PAM**, but for
what a flank wider than the PAM cannot index at a contig's very end and start: a + strand
PAM with 1 to *K − L* nt after it to the end, and a − strand PAM in the first *K − L* nt,
are not indexed, and one exactly *K − L* nt from the start is reported on the + strand at
another position. On hg38 at *K* = 8 that is 145 + 101 NGG and 308 + 241 NRG sites, none
on a primary chromosome; CRISPRme+ says so when a search reads such a flank.
Measured genome-wide at `NRG` 4 mismatches + 1 DNA bulge + 1 RNA bulge over HGDP +
1000 Genomes with per-sample resolution: the flank-8 arm, the flank-3 arm and the
dedicated `NRG_3` arm give **byte-identical** `integrated_results.tsv`
files, 182,369 rows each (the dedicated arm's licence-gated COSMIC column, empty on its
fresh install, set aside); that guide has no site at a contig's end. The archives here
were built by CRISPRitz 2.9.1, which keeps a target whose flank holds an `N`: an index of
this flank built by an earlier release misses those sites, and CRISPRme+ 2.8.1 warns when
a search reads one.

## Verifying a download

```bash
sha256sum -c --ignore-missing SHA256SUMS   # from inside indexes/
```

## License & citation

Released under the **MGB Open Access License 1.0** (see `license_link`).
Commercial use requires a separate licence. Cite CRISPRme+ and the underlying
variant panels (HGDP, 1000 Genomes) as described in the
[CRISPRme+ repository](https://github.com/pinellolab/crisprme-plus).
