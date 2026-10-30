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
flank instead of a fixed PAM, so a single build answers searches for any PAM no
longer than that flank — `NGG`, `NRG`, `NRR`, `NGA`, `NAG`, `TTTV`, `NNGRRT` and
so on — with the PAM applied at search time by `--prune-pam`.

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
│   └── SHA256SUMS                               # checksum per archive
└── README.md
```

## Indexes

| Index | Flank | Serves | Variants |
|---|---|---|---|
| `NNNNNNNN_3_hg38` | 8 nt | any PAM of length ≤ 8 | none (reference-only) |
| `NNNNNNNN_3_hg38+hg38_HGDP1kGP` | 8 nt | any PAM of length ≤ 8 | HGDP + 1000 Genomes, genotyped-phased (SHAPEIT5 re-phased) |

Both carry the matching `_INDELS` companion, and the variant index is accompanied
by the per-sample genotype store already published in the main repo.

**Orientation limit.** A wider flank currently serves a **3′** PAM only. For a 5′
PAM (Cas12a-style) the engine refuses a *shorter* PAM off a wider flank, because
the emitted sequence is trimmed at the end while a 5′ pad sits at the front — so
the reported PAM region would carry pad bases. An equal-width 5′ PAM works. Use a
dedicated index for a shorter 5′ PAM.

## Usage

```bash
# fetch the general index (routes here automatically on flank width)
crisprme.py download --what index --index-name NNNNNNNN_3_hg38+hg38_HGDP1kGP --path .

# then search with whatever PAM you need -- no rebuild
crisprme.py complete-search --genome Genomes/hg38 \
  --vcf vcf_config.txt --guide guide.txt --pam PAMs/20bp-NRG-SpCas9.txt \
  --bMax 2 --mm 4 --bDNA 1 --bRNA 1 --prune-pam \
  --annotation ... --output my_search
```

`--prune-pam` is **CLI-only**: the web form does not pass the flag, so a general
index cannot be used from the web interface.

## Equivalence claim

A PAM-less index of flank *K* searched with `--prune-pam` and a PAM of length
*L ≤ K* produces the **same report as a dedicated index built for that PAM**.
Verified genome-wide at `NRG` 4 mismatches + 1 DNA bulge + 1 RNA bulge over HGDP +
1000 Genomes with per-sample resolution: the flank-8 arm, the flank-3 arm and the
dedicated `NRG_3` arm give **byte-identical** `integrated_results.tsv`
(189,950 rows), with every indel row justified by its own window.

## Verifying a download

```bash
sha256sum -c SHA256SUMS            # from inside indexes/
```

## License & citation

Released under the **MGB Open Access License 1.0** (see `license_link`).
Commercial use requires a separate licence. Cite CRISPRme+ and the underlying
variant panels (HGDP, 1000 Genomes) as described in the
[CRISPRme+ repository](https://github.com/pinellolab/crisprme-plus).
