[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crisprme/README.html)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/pinellolab/crisprme-plus?include_prereleases)
![Conda](https://img.shields.io/conda/dn/bioconda/crisprme)
![license](https://img.shields.io/badge/license-MGBOA%201.0-lightgrey)

<p align="center">
  <img src="assets/readme/crisprme-logo.png" alt="CRISPRme" width="700"/>
</p>

# CRISPRme+ 
CRISPRme is a comprehensive tool designed for thorough off-target assessment in 
CRISPR-Cas systems. It is available as a command-line interface and an offline tool
with a locally deployable web interface. CRISPRme integrates human genetic variant 
datasets with orthogonal genomic annotations to predict and prioritize potential 
off-target sites at scale. CRISPRme accounts for single-nucleotide variants (SNVs) 
and indels, considers *bona fide* haplotypes, and allows for spacer:protospacer 
mismatches and bulges, making it well-suited for both population-wide and personal 
genome analyses. CRISPRme automates the entire workflow, from data download to 
executing the search, and delivers detailed reports complete with tables and figures 
through an interactive web-based interface.

### ✨ Highlights

- **Dictionary-less variant-aware search** — a compact allele-frequency registry + genotype store ship with the index, so variant off-target search (with allele frequencies, rsIDs, and per-dataset provenance) runs out of the box. ([methods](METHODS.md#1-variant-aware-dictionary-less-data-model))
- **Co-occurring off-targets** — off-targets that need two nearby variants together on the same chromosome copy (*in cis*, i.e. the same haplotype), so one individual carries both — **SNP+indel** (on by default; opt out with `CRISPRME_INDEL_SNP=0`) or **SNP+SNP** — are detected and reported **CONFIRMED** (phasing proves cis) or **PUTATIVE** (co-occurrence possible but cis unproven — unphased genotyped panels still name the carriers; sites-only panels like `mega` give a conservative min-AF bound with no carriers), with joint AF, and surfaced in the report. ([methods](METHODS.md#4-haplotype-scanning-observed-haplotype-enumeration))
- **Complementary production indexes** — the **recommended default** is the phased **HGDP + 1000 Genomes** index (`NRG_3_hg38+hg38_HGDP1kGP`: 4,091 samples, jointly SHAPEIT5-re-phased → CONFIRMED cis + named per-sample carriers throughout), with the **HPRC pangenome** and the sites-only five-source **mega** alternatives. All carry searchable indels genome-wide with SNP+SNP / SNP+indel co-occurrence — see the [index comparison table below](#-quickstart--web-interface-in-docker-no-conda-no-giant-build) and [`docs/PRECOMPUTED_INDEXES.md`](docs/PRECOMPUTED_INDEXES.md). (Per-source retrieval recipes for 1000G 2019/2021, HGDP, gnomAD v4.1, TOPMed and All-of-Us are in [`paper/reproducibility/variant-datasets/`](paper/reproducibility/variant-datasets/).)
- **Context-aware analysis mode: per-sample on phased indexes by default; `--population-level` to opt out** — as of 2.8.0 `complete-search` auto-selects **per-sample** (CONFIRMED cis phasing + named carriers + exact joint AF) on a genotyped + **PHASED** index (e.g. the HGDP + 1000 Genomes default) and **population-level** (worst-possible representatives per variant window, lossless for detection, tractable on any panel) otherwise (sites-only / unphased / hybrid / reference-only). Add **`--population-level`** (CLI) or pick **Population-level** in the web form to force the fast worst-case screen on a phased index; add **`--per-sample`** to force resolution on any index (it cannot name carriers on a sites-only index). This is a **genotype-resolution** dial, not a speed dial — for a single guide the two modes take about the same time; the resolution cost shows up only on dense/aggregate panels. The report states which mode was used. ([methods](METHODS.md#population-level-analysis-default-and---per-sample-genotype-resolution))
- **Lean default output** — the report, tables, mismatch/bulge matrix and validation panel are built from the **best alignment per locus**. The exhaustive *alternative-alignments* file (`..._all_results_with_alternative_alignments.tsv`, the non-best alignments per locus — expensive on large searches) is emitted by mode: **off** in population-level runs (`--population-level` on a genotyped index, or the default on sites-only/unphased indexes), **on** in per-sample runs (the default on genotyped+phased indexes, or under `--per-sample`) — on **both the CLI and the web**; override on the CLI with **`--alt-alignments`** / **`--no-alt-alignments`**. Each index also self-describes its `data_type` (sites-only / genotyped-unphased / genotyped-phased / hybrid), shown in the web variant-dataset selector.
- **Cancer-gene annotations: IntOGen by default, COSMIC by licence** — off-targets in cancer-driver genes are flagged via **IntOGen** (CC0, on by default); **COSMIC** (Cancer Gene Census) is licence-gated and **excluded by default** — enable it in **Settings** or with `crisprme.py cosmic-license enable`. Alongside ENCODE SCREEN v4, GENCODE and DHS. ([methods](METHODS.md#6-functional-annotation-of-off-targets))
- **Shareable off-target report** — every run auto-generates a self-contained HTML report (summary, calibrated per-scorer risk card, precision-recall curves, **recommended validation panel**, annotated top hits, per-tier downloads, legend). The panel keys bulges on **CRISPR-Bulge at a principled sensitive threshold** — the **95%-recall operating point** on a 40-guide held-out benchmark (`CB ≥ 1.9×10⁻⁴`; selectable to 90 % / 95 % / 99 % recall via `CRISPRME_CB_RECALL_TARGET`) — plus a high-CFD rescue for no-/single-bulge sites and all perfect/near-perfect and functional / cancer-gene hits, so the shipped panel *is* the recommended default set. ([methods](METHODS.md#7-shareable-off-target-assessment-report))
- **One-command Docker web interface + indexes on demand** — no conda, no giant build; pull reference data + prebuilt indexes from the HuggingFace CDN (`crisprme.py download`), or **build & publish your own** dict-less indexes (`build-index-only` / `publish-index --dictless`), and manage genomes / indexes / VCFs / annotations / PAMs from the browser **Data Manager**. ([protocol](docs/PRECOMPUTED_INDEXES.md))
- **Personal-assembly search (CLI + web)** — `assembly-search` scans a fully assembled diploid genome (paternal + maternal haplotype FASTAs) directly, lifts each haplotype to hg38, and reconciles them (homozygous- / heterozygous-equivalent + haplotype-non-mappable sites) — no VCF inference. Each individual is stored as one self-contained `Assemblies/<individual>/` folder + a `metadata.json` (fetch an HPRC individual, or bring your own as a `.zip`/`.tar.gz` bundle); `--assembly-individual <name>` runs a registered one from the CLI. Available on the website (a "Personal assembly" tab + a Settings card) with a combined report. An optional `--gene_annotation` adds the nearest gene + signed distance (`Gene` / `Gene_distance_kb`) to every hg38-mapped site. ([web guide](docs/crisprme_web_interface_user_guide.md#searching-a-personal-assembly-genome-assembly-search))
- **Also** — merged multi-dataset VCF panels (per-dataset provenance + per-haplotype phasing), PAM-geometry-aware (pamless) indexes, bounded-complexity controls (`--max-total-edits` + a high-variant-density skip), and optional email notifications.

### What changed vs stock CRISPRme (migration guide)

Coming from the original CRISPRme? The baseline below is **CRISPRme 2.1.14** — the last
stable release of the original [`pinellolab/CRISPRme`](https://github.com/pinellolab/CRISPRme)
before CRISPRme+ (installable with `mamba install crisprme`). Here are the behavior deltas
that matter most:

| Topic | CRISPRme 2.1.14 (last stable before CRISPRme+) | CRISPRme+ (now) | What to do |
|---|---|---|---|
| Off-target ML score | CFD + **CRISTA** | CFD + **CRISPR-Bulge** (OrensteinLab, *NAR* 2024; MIT — more accurate, esp. on bulges), run in a dedicated `cbulge` env | Nothing: `install_from_source.sh` / the Docker image provision it; or create it later with `crisprme.py scorer-env create` |
| Analysis mode | single search, **always** per-sample genotype-resolved (named carriers; no mode flag) | **context-aware default** — per-sample on genotyped+phased panels (CONFIRMED cis + named carriers + exact joint AF), population-level on sites-only/unphased; **`--population-level`** opts out of per-sample, **`--per-sample`** forces it | Add `--per-sample` to reproduce the old per-sample output ([methods](METHODS.md#population-level-analysis-default-and---per-sample-genotype-resolution)) |
| SNP+SNP / SNP+indel co-occurrence | not detected | **detected + reported** — **CONFIRMED** (phasing proves cis) / **PUTATIVE**, with joint AF; on by default | Opt out with `CRISPRME_INDEL_SNP=0` ([methods](METHODS.md#4-haplotype-scanning-observed-haplotype-enumeration)) |
| Alternative-alignments TSV | always written | **off in population-level runs, on in per-sample runs** (per the context-aware default; both CLI and web) | Force either way on the CLI with `--alt-alignments` / `--no-alt-alignments` ([methods](METHODS.md#7-shareable-off-target-assessment-report)) |
| Variant datasets / default index | build-your-own from `--vcf` configs (e.g. 1000G + HGDP) | **precomputed, downloadable** indexes on HuggingFace; recommended default **`NRG_3_hg38+hg38_HGDP1kGP`** (phased HGDP + 1000 Genomes, 4,091 samples → CONFIRMED cis + named carriers), with `HPRC` / sites-only `mega` alternatives | `crisprme.py download --what index --index-name …` (Quickstart step 2) ([details](docs/PRECOMPUTED_INDEXES.md)) |
| License | AGPL-3.0 for academic/non-commercial use **plus a separate commercial license required for commercial use** (dual-licensed — AGPL never covered commercial use) | **MGB Open Access License 1.0** — free for non-commercial, non-revenue-generating, academic use; commercial use **likewise requires a separate commercial license** | See [§6 License](#6-license) |
| Report | CFD / CRISTA summaries | adds an **`Observed`** column (`reference` / `N carrier(s)` / `observed` / `putative`), CRISPR-Bulge panels, and a shareable HTML report | Read it in the report ([methods](METHODS.md#7-shareable-off-target-assessment-report)) |

> **Note on the retired flags (CRISPRme+ 2.5.x users only).** `--fast` / `--full` existed
> only in the intermediate CRISPRme+ 2.5.x line — never in stock CRISPRme 2.1.14 — and were
> removed in 2.5.4 (they are **not** aliases). `--per-sample` is the replacement for `--full`.
> `complete-search` now prints a **WARNING** if a saved script passes the retired
> `--fast`/`--full`, so a pre-2.5.4 command no longer silently changes the analysis mode;
> mode is now the context-aware default unless `--population-level`/`--per-sample` is given.

### Which analysis mode + which index should I use?

- **Analysis mode.** The **default (no flag) is context-aware**: per-sample on a
  genotyped+PHASED index (`HGDP1kGP` — the default; also `HPRC`) resolving
  CONFIRMED cis + named carriers + exact joint AF, or population-level otherwise
  (sites-only/unphased/hybrid/reference-only; lossless for detection, works on any index).
  Add **`--population-level`** to opt out of per-sample on a phased index (e.g. a very
  dense / high-edit / genome-wide search); add **`--per-sample`** to force resolution. The
  two flags are mutually exclusive. It is a genotype-resolution dial, not a speed dial.
- **Index.** You only need one; pick by what you're after (full detail in the index
  table in the [Quickstart](#-quickstart--web-interface-in-docker-no-conda-no-giant-build)
  below):

  | Index | Best for |
  |---|---|
  | `HGDP1kGP` *(recommended default)* | most cases — phased HGDP + 1000 Genomes, 4,091 samples → CONFIRMED cis + named carriers throughout |
  | `HPRC` | assembly / pangenome variation |
  | `mega` | **escalate here only** for a widest-provenance worst-case screen (sites-only, all co-occurrences PUTATIVE, no carriers) |

> ⚠️ **Note**  
> The original public CRISPRme web service is no longer available.  
> All functionalities are now accessible via the command-line interface or the 
locally hosted web interface. See [`docs/DOCKER_QUICKSTART.md`](docs/DOCKER_QUICKSTART.md) for a quick guide on
deploying the web interface locally.

## ⚡ Quickstart — web interface in Docker (no conda, no giant build)

### Easiest: one-line install → a clickable app (no terminal after this)

Not comfortable with the command line? Run **one** command once, and it installs a
clickable **CRISPRme** app (into Applications / Start Menu) that manages everything
— you never touch the terminal again. Requires [Docker Desktop](https://www.docker.com/products/docker-desktop/).

**macOS / Linux** (Terminal):
```bash
curl -fsSL https://pinellolab.github.io/crisprme-plus/install.sh | bash
```
**Windows** (PowerShell):
```powershell
irm https://pinellolab.github.io/crisprme-plus/install.ps1 | iex
```
Then open **CRISPRme** (Applications on macOS; Desktop / Start Menu on Windows) — a
small window with three buttons: **Start**, **Update**, **Stop**. Click **Start**: the
first time it downloads the reference + variant data automatically (~45 GB download, once;
~85 GB on disk once installed — keep ~100 GB free; a reference-only setup is far smaller, ~25 GB), then
opens the web interface at http://localhost:8080; every Start after that is instant.
(Memory: the **whole-genome *variant-aware* search** needs **~64 GB of RAM** — set that
under Docker Desktop → Settings → Resources → Memory. **On Windows with the WSL 2 backend
there is no memory slider**; Windows caps the WSL 2 VM at roughly half your RAM, so set
`memory=` in `C:\Users\<you>\.wslconfig`, run `wsl --shutdown`, and verify with
`(docker info --format '{{.MemTotal}}') / 1073741824` — see the
[Docker Quickstart](docs/DOCKER_QUICKSTART.md#1-install-docker-one-time). Many current
laptops ship with 36–128 GB and can run it; otherwise use a workstation/HPC.
**Reference-only or single-chromosome / target-region variant** searches fit **~16 GB**,
so any laptop handles those.)
Details: [`install/`](install/).

> **GPU acceleration (optional).** The default image `pinellolab/crisprme:latest` is
> **GPU-capable**: on an NVIDIA host, add `--gpus all` to `docker run` and pass
> `--compute-backend cuda` to run the CRISPR-Bulge scorer on the GPU; it **falls back to
> CPU automatically** when no GPU is present. The **amd64** image bundles CUDA; **arm64**
> (Apple Silicon) and Docker-on-Mac are always CPU. The clickable installers auto-add
> `--gpus all` when Docker exposes a GPU. GPU pays off on large/dense searches; a tiny
> single-guide run is often faster on CPU due to init overhead. (See the
> [`--compute-backend {cpu,cuda,metal}`](#input-arguments) flag.)

---

### Or run the Docker commands yourself

New to CRISPRme? Get the point-and-click web interface running in a few commands
— just install Docker, then:

> **Allocate memory first.** In **Docker Desktop → Settings → Resources → Memory**:
> **16 GB** is enough for reference-only / single-chromosome runs, but the
> genome-wide 1000G+HGDP variant search shown here peaks around **68–75 GB**, so
> give Docker **64 GB+** for it (32 GB will OOM). On **Windows with the WSL 2 backend there is no memory slider** — set `memory=` in `C:\Users\<you>\.wslconfig` instead, then `wsl --shutdown`; see the [Docker Quickstart](docs/DOCKER_QUICKSTART.md#1-install-docker-one-time). (On Linux/Docker Engine the
> container uses host memory directly.) See [System Requirements](#0-system-requirements).

```bash
mkdir -p ~/crisprme && cd ~/crisprme    # your data folder — ALL shared inputs + per-search Results live here (see docs/DATA_LAYOUT.md)
# 1) fast-download the reference data + reference index (minutes, HuggingFace CDN)
#    (this does NOT include the variant index — that is step 2)
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py download --what all --path /DATA
# 2) grab the prebuilt SpCas9 (NRG = NAG+NGG) indexes so no long index build is needed:
#    the reference index, and the compact variant-aware hg38 + HGDP + 1000
#    Genomes index (the recommended default; phased HGDP + 1000 Genomes, 4,091 samples
#    jointly SHAPEIT5-re-phased → CONFIRMED cis + named per-sample carriers throughout,
#    including HGDP's seven genetic regions; the broadest fully-phased panel of real
#    observed haplotypes)
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py download --what index --index-name NRG_3_hg38 --path /DATA
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP --path /DATA
#    ...or the HPRC pangenome index (232 phased assembly-derived genomes incl. CHM13,
#    CONFIRMED cis + named carriers, pangenome-specific variation, searchable indels):
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py download --what index --index-name NRG_3_hg38+hg38_HPRC --path /DATA
#    ...or the all-source sites-only "mega" index (5 datasets, per-dataset AF + AF_max,
#    searchable indels, SNP+SNP & SNP+indel PUTATIVE co-occurrence):
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py download --what index --index-name NRG_3_hg38+hg38_mega --path /DATA
# 3) launch the web interface, then open http://127.0.0.1:8080
#    (on an NVIDIA host add `--gpus all` here to run the CRISPR-Bulge scorer on the GPU)
docker run --rm -v "${PWD}:/DATA" -w /DATA -p 8080:8080 -it pinellolab/crisprme:latest crisprme.py web-interface
```

> **`:latest` vs a pinned tag.** These commands use `pinellolab/crisprme:latest` — the
> auto-updating, GPU-capable default image. For a fully reproducible run, pin the exact
> release instead, e.g. `pinellolab/crisprme:v2.8.0`.

> **On HPC with Apptainer/Singularity.** Use **`apptainer run`**, which invokes the image
> entrypoint so `crisprme.py` is on `PATH`, binding your data dir to `/DATA` — e.g.
> `apptainer run --bind "${PWD}:/DATA" --pwd /DATA docker://pinellolab/crisprme:v2.8.0 crisprme.py --version`
> (add `crisprme.py web-interface` / `complete-search …` the same way). Note that
> `apptainer exec … crisprme.py …` **bypasses** the entrypoint, so `crisprme.py` is not found on
> `PATH` — prefer `apptainer run` (or prefix `bash -lc`).

**The prebuilt SpCas9 (NRG = NAG+NGG) variant indexes.** All carry searchable
indels genome-wide with SNP+SNP / SNP+indel co-occurrence; download whichever you need with
`--what index --index-name <name>` (or grab more later from the web **Settings**):

| Index (`--index-name`) | Variant sources | Resolution | Best for |
|---|---|---|---|
| **HGDP + 1000 Genomes** — `NRG_3_hg38+hg38_HGDP1kGP` *(recommended default)* | HGDP (925) + 1000 Genomes (3,166), genotyped, **jointly SHAPEIT5-re-phased** — 4,091 samples | **CONFIRMED cis** throughout; named per-sample carriers + exact joint AF (including HGDP's seven genetic regions, with correct chrX PAR allele frequencies) | most use cases; the broadest fully-phased panel of real observed haplotypes |
| **HPRC pangenome** — `NRG_3_hg38+hg38_HPRC` | HPRC Release 2 Minigraph-Cactus (`hprc-v2.0-mc-grch38`) — 232 assembly-derived genomes incl. CHM13, phased | CONFIRMED cis + named carriers; captures **pangenome-specific** variation | assembly-derived / pangenome variation |
| **mega (sites-only)** — `NRG_3_hg38+hg38_mega` | 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us (aggregate) | **sites-only** — union of allele frequencies, no genotypes; every multi-variant/co-occurring off-target is **PUTATIVE** with a min-AF bound and **no named carriers** | a conservative worst-case screen when you must not miss a rare allele from ANY of five databases (escalate here after a genotyped-index run) |

> **Which variant index should I pick?** You only need **one** — if unsure, use the
> **recommended default HGDP + 1000 Genomes** (`NRG_3_hg38+hg38_HGDP1kGP`, downloaded in step 2):
> phased HGDP + 1000 Genomes, 4,091 samples jointly SHAPEIT5-re-phased → CONFIRMED cis +
> named per-sample carriers throughout, including HGDP's seven genetic regions and with
> correct chrX PAR allele frequencies — the broadest fully-phased panel of real observed
> haplotypes, covering most use cases. Pick **HPRC** (`+hg38_HPRC`) for pangenome / assembly-derived
> variation. **Escalate to mega** (`+hg38_mega`) only for a widest-provenance worst-case
> screen across five databases — it is **sites-only**, so every co-occurring off-target is
> **PUTATIVE** (a worst-case reconstruction that may or may not exist in any real
> individual) with a min-AF bound and **no named carriers**; run it after a genotyped-index
> run when you need maximal allele-frequency breadth. Grab the others any time from the web
> **Settings** or by re-running the download with a different `--index-name`.

**Full step-by-step (with variants, more indexes, troubleshooting):
[`docs/DOCKER_QUICKSTART.md`](docs/DOCKER_QUICKSTART.md).**

### …or run the same search from the command line

Prefer the CLI? After the download steps above (1 & 2), run the variant-aware
example search and generate the shareable report — no web UI needed:

```bash
# a genome-wide SpCas9 (NRG) search over the hg38 + HGDP + 1000 Genomes jointly re-phased index, up to 4
# mismatches + 1 DNA / 1 RNA bulge (a 6-edit budget), with --per-sample genotype resolution
# (named carriers + CONFIRMED cis + exact joint allele frequency), plus rsIDs
# and functional annotations. (Downloaded the +hg38_HGDP1kGP index in step 2 above.)
echo "CTAACAGTTGCTTTTATCACNNN" > guide.txt
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt --guide guide.txt \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --annotation Annotations/dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 4 --bDNA 1 --bRNA 1 --per-sample --output my_search --thread 8
# ^ the recommended default search: up to 4 mismatches + 1 DNA and 1 RNA bulge (a 6-edit
#   budget, matching the default --max-total-edits 6, so nothing is pruned). This sits in the
#   CRISPR-Bulge scorer's validated single-bulge domain. Raise the caps for a deeper search,
#   but also raise --max-total-edits to their sum or alignments over the budget are PRUNED
#   (e.g. --mm 6 --bDNA 2 --bRNA 2 is a 10-edit budget → add --max-total-edits 10 to keep all).
#   Uses --per-sample here to resolve per-sample genotypes — CONFIRMED cis phasing throughout
#   (the whole HGDP + 1000 Genomes panel is jointly re-phased) + named carrier samples + exact joint
#   allele frequency — recommended for genotyped panels / clinical validation. On HGDP1kGP (a
#   genotyped+phased index) per-sample is already the default; pass --population-level to switch to
#   the worst-possible representatives (removes the per-haplotype enumeration wall that makes
#   dense/aggregate panels intractable). --per-sample is a genotype-resolution mode,
#   not a speed mode: for a single guide the runtimes are comparable, but it can be intractable on
#   dense/aggregate panels; on a sites-only index (mega) it cannot resolve carriers and is inert.
#   COSMIC cancer annotations ship in the bundle but are OFF by default (licence-gated) — enable
#   with `crisprme.py cosmic-license enable`; IntOGen cancer-driver flags are on by default.

# build the self-contained, shareable HTML report (report.html + a data/ folder)
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest crisprme.py generate-report \
  --result-dir Results/my_search
# -> open Results/my_search/<jobid>_report.zip, then report.html
```

`download --what all/index` writes `list_vcf.txt` and `list_samplesID.txt` for you,
so the search above works out of the box on a fresh install.

## Table Of Contents

0 [System Requirements](#0-system-requirements)
<br>1 [Installation](#1-installation)
<br>&nbsp;&nbsp;1.1 [Install CRISPRme via Conda/Mamba](#11-install-crisprme-via-condamamba)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.1 [Installing Conda or Mamba](#111-installing-conda-or-mamba)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.2 [Installing CRISPRme](#112-installing-crisprme)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.3 [Updating CRISPRme](#113-updating-crisprme)
<br>&nbsp;&nbsp;1.2 [Install CRISPRme via Docker](#12-install-crisprme-via-docker)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.2.1 [Installing Docker](#121-installing-docker)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.2.2 [Pulling the CRISPRme Docker Image](#122-pulling-the-crisprme-docker-image)
<br>&nbsp;&nbsp;1.3 [Install CRISPRme from source (without Bioconda)](#13-install-crisprme-from-source-without-bioconda)
<br>2 [Usage](#2-usage)
<br>&nbsp;&nbsp;2.1 [Directory Structure](#21-directory-structure)
<br>&nbsp;&nbsp;2.2 [CRISPRme Functions](#22-crisprme-functions)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.1 [Complete Search](#221-complete-search)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.2 [Complete Test](#222-complete-test)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.3 [Off-target sites validation Test](#223-off-target-sites-validation-test)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.4 [Targets Integration](#224-targets-integration)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.5 [GNOMAD Converter](#225-gnomad-converter)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.6 [Generate Personal Card](#226-generate-personal-card)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.7 [Setup Legacy Database](#227-setup-legacy-database)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.8 [Local Web Interface](#228-local-web-interface)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.9 [Generate Report](#229-generate-report)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.10 [Reference-Index and Data-Distribution Commands](#2210-reference-index-and-data-distribution-commands)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2.2.11 [Targets Scoring](#2211-targets-scoring)
<br>3 [Test](#3-test)
<br>&nbsp;&nbsp;3.1 [Quick Test](#31-quick-test)
<br>&nbsp;&nbsp;3.2 [Detailed Test](#32-detailed-test)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 3.2.1 [Single Chromosome Test](#321-single-chromosome-test)
<br>&nbsp;&nbsp;&nbsp;&nbsp; 3.2.2 [Full Genome Test](#322-full-genome-test)
<br>4 [Citation](#4-citation)
<br>5 [Contacts](#5-contacts)
<br>6 [License](#6-license)
<br>• [Methods &amp; validation](METHODS.md) — the self-contained technical reference (variant-aware data model, haplotype scanning, worst-case scoring, annotation, report, scoring assumptions &amp; limitations)

## 0 System Requirements

Memory scales with the search — pick your tier (measured figures):

| Search | RAM | Disk |
|---|---|---|
| **Reference-only / single-chromosome** (toy or reference-genome runs) | **8–16 GB** | ~22 GB (reference-only install) |
| **Whole-genome, variant-aware** (e.g. HGDP + 1000 Genomes (`NRG_3_hg38+hg38_HGDP1kGP`); the web/Quickstart default) | **64 GB minimum, ~80 GB for high `--max-total-edits`** | ~85 GB (full variant-aware install) |

The variant-aware genome-wide search is the memory-heavy one: the CRISPRitz search
process alone peaks around **68–75 GB** loading and traversing the variant index, so
**16 GB or 32 GB will OOM a real genome-wide variant search** — allocate 64 GB+.
Reference-only runs are light (a few GB). To estimate the cost of a large search up
front, see `docs/SCALABILITY_ANALYSIS.md`.

## 1 Installation

> **Which version do I get?** For **CRISPRme+ (2.8.0, this release)** use **Docker**
> (the [Quickstart](#-quickstart--web-interface-in-docker-no-conda-no-giant-build) above, or
> §1.2). **Conda/Bioconda currently installs the stable 2.1.x line (Python 3.8), not the
> 2.8.0 line** — use it only if you specifically want the stable release. If in doubt,
> use Docker.

This section outlines the steps to install CRISPRme, tailored to suit different 
operating systems. Select the method that best matches your setup:

- [Install CRISPRme via Docker (compatible with all operating systems — recommended for 2.8.0)](#12-install-crisprme-via-docker)

- [Install CRISPRme from source (§1.3)](#13-install-crisprme-from-source-without-bioconda) — the CRISPRme+ (2.8.0) line on Python 3.11

- [Install CRISPRme via Conda/Mamba (Linux; installs the stable 2.1.x line)](#11-install-crisprme-via-condamamba)

### 1.1 Install CRISPRme via Conda/Mamba (stable 2.1.x line only)
---

Conda/Bioconda installs **only the stable 2.1.x line** (Python 3.8) — *not* the
CRISPRme+ (2.8.x) line. Use this path only if you specifically want the stable
release; for CRISPRme+ 2.8.0 use [source (§1.3)](#13-install-crisprme-from-source-without-bioconda)
or [Docker](#12-install-crisprme-via-docker).

This section is organized into three subsections to guide you through the installation 
and maintenance of CRISPRme:

- [Installing Conda or Mamba](#111-installing-conda-or-mamba):
  <br>This subsection provides step-by-step instructions to install either 
  Conda or Mamba. Begin here if you do not have these package managers installed 
  on your machine.

- [Installing CRISPRme](#112-installing-crisprme):
  <br>Once you have Conda or Mamba set up, proceed to this subsection for detailed 
  instructions on creating the CRISPRme environment and installing the necessary 
  dependencies.

- [Updating CRISPRme](#113-updating-crisprme):
  <br>Learn how to update an existing CRISPRme installation to the latest version, 
  ensuring access to new features and bug fixes.

#### 1.1.1 Installing Conda or Mamba
---

Before installing CRISPRme, ensure you have either Conda or Mamba installed on 
your machine. Based on recommendations from the Bioconda community, we highly 
recommend using Mamba over Conda. Mamba is a faster, more efficient drop-in replacement 
for Conda, leveraging a high-performance dependency solver and components 
optimized in C++.

**Step1: Install `Conda` or `Mamba`**

- To install `Conda`, refer to the official installation guide:
  <br>[Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

- To install `Mamba`, refer to the official installation guide: 
  <br>[Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

**Step 2: Configure Bioconda Channels**

Once `Mamba` is installed, configure it to use Bioconda and related channels by 
running the following one-time setup commands:
```bash
mamba config --add channels bioconda
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```
> **Note:** If you prefer to use `Conda`, replace `mamba` with `conda` in the 
commands above. (Per current Bioconda guidance we do **not** add the `defaults`
channel, and `conda-forge` is given the highest priority.)

By completing these steps, your system will be fully prepared for installing CRISPRme.

#### 1.1.2 Installing CRISPRme
---

> **CRISPRme+ runs on Python 3.11 and installs from source** — the build
> compiles CRISPRitz 2.9.0 and installs both tools into a conda environment. A native
> Bioconda package for the CRISPRme+ line is **in preparation**; until it lands, the Bioconda
> `crisprme` package installs the last **stable 2.1.x** line (Python 3.8), **not** the
> CRISPRme+ (2.8.x) line.

To create the CRISPRme+ conda environment, follow **[1.3 Install CRISPRme from source](#13-install-crisprme-from-source-without-bioconda)**
(`git clone` → `mamba env create -f environment.yml` (Python 3.11) → `bash install_from_source.sh`),
or use the **[Docker quickstart](docs/DOCKER_QUICKSTART.md)** for the fastest path.

Verify the installation:

```bash
crisprme.py --version
crisprme.py
```

#### 1.1.3 Updating CRISPRme
---

> **Updating CRISPRme+ (2.8.0)?** Conda/Bioconda cannot update you to the CRISPRme+
> line — update via the [source build (§1.3)](#13-install-crisprme-from-source-without-bioconda)
> or [Docker](#12-install-crisprme-via-docker). The steps below update only the stable
> 2.1.x Bioconda line.

To update an existing **stable 2.1.x** installation using `Mamba` or `Conda`, follow
the steps below:

**Step 1: Check the Latest Version**

Visit the CRISPRme README to identify the latest version of the tool.

**Step 2: Update CRISPRme**

Run the following command in your terminal, replacing <latest_version> with the 
desired version number:
```bash
mamba install crisprme=<latest_version>
```

This updates within the **stable 2.1.x** Bioconda line (latest is `crisprme=2.1.14`):
```bash
mamba install crisprme=2.1.14
```
If you're using `Conda`, replace `mamba` with `conda` in the commands above.

**Step 3: Verify the Update**

After the update completes, ensure the installation was successful by checking the 
version:
```bash
crisprme.py --version
```
If the displayed version matches the one you installed, the update was successful.

### 1.2 Install CRISPRme via Docker
---

This section is organized into two subsections to guide you through the setup of 
**CRISPRme** using Docker:  

- [Installing Docker](#121-installing-docker): 
<br>Provides step-by-step instructions for installing Docker on your system, 
ensuring compatibility with all operating systems, including Linux, macOS, and Windows.  

- [Pulling the CRISPRme Docker Image](#122-pulling-the-crisprme-docker-image): 
  <br>Explains how to download the prebuilt CRISPRme Docker image to set up a 
  containerized environment for seamless execution.  

Follow the subsections in order if Docker is not yet installed on your machine. 
If Docker is already installed, skip to the second subsection.

#### 1.2.1 Installing Docker 
---

MacOS and Windows users are encouraged to install [Docker](https://www.docker.com/get-started) 
to use CRISPRme. Linux users may also choose Docker for convenience and compatibility.

Docker provides tailored distributions for different operating systems. Follow the 
official Docker installation guide specific to your OS:

- [MacOS Installation Guide](https://docs.docker.com/docker-for-mac/install/)
- [Windows Installation Guide](https://docs.docker.com/docker-for-windows/install/)
- [Linux Installation Guide](https://docs.docker.com/engine/install/ubuntu/)

**Apple Silicon (M1/M2/M3) support**

The CRISPRme Docker image is now published as a multi-architecture image with a
native `linux/arm64` build, so it runs **natively on Apple Silicon Macs with no
emulation** (no Rosetta / QEMU, no "platform mismatch" warning). Both the
`crispritz` dependency (native `linux-aarch64` bioconda build) and `crisprme`
itself (a `noarch` package) install per-architecture without any source
compilation. Docker automatically selects the correct architecture when you
pull or run the image, so no special flags are required.

**Docker Desktop memory allocation (important)**

CRISPRme's searches are memory intensive. By default Docker Desktop caps the
memory available to containers, which can cause searches to fail with
out-of-memory errors. Before running large analyses, raise the limit in
**Docker Desktop → Settings → Resources → Memory** — **16 GB** is fine for
reference-only runs, but a genome-wide **variant** search needs **64 GB+** (it peaks
at ~68–75 GB). On **Windows with the WSL 2 backend there is no memory slider** — set `memory=` in `C:\Users\<you>\.wslconfig` instead, then `wsl --shutdown`; see the [Docker Quickstart](docs/DOCKER_QUICKSTART.md#1-install-docker-one-time). CRISPRme prints a non-fatal
`WARNING:` at the start of a `complete-search` / `complete-test` run when it
detects less than 32 GB of total memory, reminding you to increase this
allocation. That threshold sits below the 64 GB a genome-wide variant search needs, so a
run that starts without a warning may still not have enough memory.
(You can override the threshold with `CRISPRME_MIN_GB` or skip the
check entirely with `CRISPRME_SKIP_MEM_CHECK=1`.)

**Linux-Specific Post-Installation Steps**

If you're using Linux, additional configuration steps are required:
1. Create the Docker Group:
  ```bash
  sudo groupadd docker
  ```

2. Add Your User to the Docker Group:
  ```bash
  sudo usermod -aG docker $USER
  ```
  Repeat this command for any additional users you want to include in the Docker 
  Group.

3. Restart Your Machine
  <br>Log out and log back in, or restart your machine to apply the changes.

**Testing Docker Installation**

Once Docker is installed, verify the setup by opening a terminal window and typing:
```bash
docker run hello-world
```

If Docker is installed correctly, you should see output like this:
```
Hello from Docker!
This message shows that your installation appears to be working correctly.

To generate this message, Docker took the following steps:
 1. The Docker client contacted the Docker daemon.
 2. The Docker daemon pulled the "hello-world" image from the Docker Hub.
 3. The Docker daemon created a new container from that image, which runs the executable that produces this output.
 4. The Docker daemon streamed this output to the Docker client, which displayed it on your terminal.

For more examples and ideas, visit:
 https://docs.docker.com/get-started/
```

#### 1.2.2 Pulling the CRISPRme Docker Image
---

After installing Docker, you can download the prebuilt CRISPRme Docker image by 
running the following command in a terminal:
```bash
docker pull pinellolab/crisprme:latest
```

This command retrieves the latest pre-built CRISPRme image from Docker Hub and sets 
it up on your system, ensuring all required dependencies and configurations are 
included. `:latest` is the auto-updating default image; it runs on **CPU** out of the
box (fully supported — this is the default scoring path), and a **GPU can optionally
accelerate the CRISPR-Bulge scorer** on an NVIDIA host (add `--gpus all` + `--compute-backend
cuda`). For a fully reproducible setup pin the exact release instead, e.g.
`docker pull pinellolab/crisprme:v2.8.0`.

Once the download is complete, the CRISPRme Docker image will be ready for use. 
To confirm the image is successfully installed, you can list all available Docker 
images by typing:
```bash
docker images
```

Look for an entry similar to the following:
```
REPOSITORY          TAG       IMAGE ID       CREATED        SIZE
pinellolab/crisprme   latest    <image_id>     <timestamp>    ~5 GB (amd64; ~1.8 GB on Apple Silicon)
```

You are now ready to run CRISPRme using Docker.

### 1.3 Install CRISPRme from source (without Bioconda)

Use this to run the CRISPRme+ line (Python 3.11 + Dash 2.x) from source — the recommended path until it is published to Bioconda, or for development. It installs the runtime dependencies into a conda environment, **builds CRISPRitz 2.9.0 from source**, and installs CRISPRme from the checkout — using the same layout the Bioconda/Docker builds use, so `crisprme.py` and `crispritz.py` end up on your `PATH` and resolve their support files correctly.

**Prerequisites:** `conda`/`mamba`, `git`, and internet access. A C++ compiler with OpenMP and every Python dependency are provided by the environment file below (no `apt`/system packages required).

```bash
# 1. clone the repository (the CRISPRme+ development line lives on the main branch)
git clone https://github.com/pinellolab/crisprme-plus.git
cd crisprme-plus

# 2. create + activate the runtime environment (pinned deps from environment.yml)
mamba env create -f environment.yml
mamba activate crisprme

# 3. build CRISPRitz 2.9.0 from source and install both tools into the env
bash install_from_source.sh

# 4. verify (both tools are now on your PATH)
crisprme.py --version
command -v crispritz.py
```

The `install_from_source.sh` script compiles the CRISPRitz C++ binaries, then copies `crisprme.py`/`crispritz.py` into `$CONDA_PREFIX/bin` and their support trees into `$CONDA_PREFIX/opt/…`, and creates the dedicated `cbulge` conda env for the CRISPR-Bulge off-target scorer (set `CRISPRME_SKIP_SCORER_ENV=1` to defer it — create it later with `crisprme.py scorer-env create`). `CRISPRITZ_REF` defaults to the **v2.9.0** tag; override it with `CRISPRITZ_REF=<tag> bash install_from_source.sh` to build a different CRISPRitz revision.

**Smoke-test the install** (downloads a small chr22 dataset, runs the full pipeline, and compares against the committed ground truth):

```bash
mkdir crisprme_test && cd crisprme_test
crisprme.py complete-test --chrom chr22 --thread 4
crisprme.py validate-test --chrom chr22
```

Notes:
- This is exactly the recipe the project `Dockerfile` automates — if you prefer full isolation, build the image instead (Section 1.2).
- The `assembly-search` subcommand lifts coordinates through chain files with the `pyliftover` library (included in `environment.yml`); no external liftOver binary is required.
- To launch the web interface after installing: `crisprme.py web-interface` (serves on port 8080 — open http://localhost:8080).

## 2 Usage

CRISPRme is a tool designed for variant- and haplotype-aware CRISPR off-target 
analysis. It integrates robust functionalities for off-target detection, 
variant-aware search, and result analysis. The tool also includes a 
user-friendly graphical interface, which can be deployed locally to streamline 
its usage.

> 📄 **Preparing inputs?** See [`docs/INPUT_FORMATS.md`](docs/INPUT_FORMATS.md)
> for PAM/guide file formats, chromosome naming, VCF and annotation
> requirements, bulge reporting, and how to estimate expensive searches.

### 2.1 Directory Structure
---

CRISPRme operates within a specific directory structure to manage input data and 
outputs efficiently. Everything lives inside **one user-owned working folder** (the CLI
`--path`, Docker `-v $PWD:/DATA`, or Apptainer `--bind $PWD:/DATA --pwd /DATA`). The
subdirectories below are **shared inputs** reused by every search; each run writes its
own outputs under `Results/<job>/`. Most of these are created for you by `download` /
first launch. See [docs/DATA_LAYOUT.md](docs/DATA_LAYOUT.md) for the full reads-vs-writes
map. The main subdirectories:

- **Genomes** 
    - **Purpose**: Stores reference genomes.
    - **Structure**: Each reference genome resides in its own subdirectory.
    - **Requirements**: The genome must be split into separate files, each 
    representing a single chromosome.

- **VCFs**
    - **Purpose**: Contains variant data in VCF format.
    - **Structure**: Similar to the Genomes directory, each dataset has a 
    dedicated subdirectory with VCF files split by chromosome.
    - **Requirements**: Files must be compressed using `bgzip` (with a `.gz`
    extension).

- **samplesIDs** 
    - **Purpose**: Lists the sample identifiers corresponding to the VCF 
    datasets.
    - **Structure**: Tab-separated files, one for each VCF dataset, specifying 
    the sample IDs.
    
- **Annotations** 
    - **Purpose**: Provides genome annotation data.
    - **Format**: Annotation files must be in BED format.

- **PAMs** 
    - **Purpose**: Specifies the Protospacer Adjacent Motif (PAM) sequences for 
    off-target search.
    - **Format**: Text files containing PAM sequences.

- **genome_library** 
    - **Purpose**: Holds the built/downloaded **search indexes** (`NRG_3_hg38`,
    `NRG_3_hg38+hg38_HGDP1kGP`, …) — the largest, most load-bearing directory, where
    `download --what index` / `build-index-only` place their output.

- **Dictionaries** 
    - **Purpose**: The variant **registry / genotype / `log_indels` tiers** that
    variant- and indel-aware post-analysis depends on (fetched alongside a variant index).

- **Assemblies** 
    - **Purpose**: Per-individual personal assemblies used by `assembly-search`
    (one folder per individual; `download --what assembly` populates it).

The directory organization required by CRISPRme is illustrated below:
<p align="center">
  <img src="https://github.com/pinellolab/crisprme-plus/blob/main/assets/readme/crisprme_dirtree.png" alt="crisprme_dirtree.png", width=400/>
</p>

### 2.2 CRISPRme Functions
---

> **Running the examples.** Each example below is the bare `crisprme.py <command> …`.
> To run it in **Docker**, prefix it with
> `docker run --rm -v "${PWD}:/DATA" -w /DATA -i pinellolab/crisprme:latest`
> (add `-p 8080:8080` for `web-interface`; on an NVIDIA host add `--gpus all` to run the
> CRISPR-Bulge scorer on the GPU). Pin `:v2.8.0` instead of `:latest` for reproducibility.
> From a **source / Conda** install, run it as-is inside the activated `crisprme`
> environment.

This section provides a comprehensive overview of CRISPRme's core functions, 
detailing each feature, the required input data and formats, and the resulting 
outputs. The following is a summary of CRISPRme's key features:

- [**Complete Search**](#221-complete-search) (`complete-search`)
  <br>Executes a genome-wide off-targets
  search across both reference and variant datasets (if specified), conducts 
  Cutting Frequency Determination (CFD) and CRISPR-Bulge analyses (if applicable), and 
  identifies candidate targets.

- **Assembly Search** (`assembly-search`)
  <br>Searches a personal / HPRC diploid assembly (paternal + maternal haplotype
  FASTAs) directly for off-targets, lifting each haplotype to hg38 and reconciling
  them — no VCF inference. See the
  [assembly-search web guide](docs/crisprme_web_interface_user_guide.md#searching-a-personal-assembly-genome-assembly-search).

- [**Generate Report**](#229-generate-report) (`generate-report`)
  <br>Builds a self-contained, easily shareable off-target assessment report for a
  completed run — a single **slim** ZIP bundling an offline `report.html` (embedded
  plots, the top-1000 off-target table, a recommended validation panel) and the
  per-tier curated TSVs, kept small and focused for sharing. Add `--report-full` to
  also bundle the heavy data (complete `integrated_results.tsv.gz`, raw co-occurrence
  dump, high-variant-density BED). It runs automatically at the end of every
  `complete-search`, and can be re-run standalone on any `Results/<jobid>` folder.

- [**Complete Test**](#222-complete-test) (`complete-test`)
  <br>Tests CRISPRme pipeline on a small input dataset or the full genome, 
  enabling users to validate the tool's functionality before performing 
  large-scale analyses.

- [**Off-target sites validation Test**](#223-off-target-sites-validation-test) (`validate-test`)
  <br>Validates off-target sites generated by the Complete Test workflow by
  comparing CRISPRme predictions against brute-force ground-truth alignments
  derived from 1000 Genomes variant data.

- [**Targets Integration**](#224-targets-integration) (`targets-integration`) 
  <br>Combines *in silico* predicted targets with experimental data to create a 
  finalized target panel.

- [**Targets Scoring**](#2211-targets-scoring) (`targets-scoring`) 
  <br>Scores a supplied list of off-target sites (e.g. Cas-OFFinder output) with the 
  CRISPRme+ CFD + CRISPR-Bulge calculations and an optional per-site variant carrier 
  overlay — no genome search.

- [**GNOMAD Converter**](#225-gnomad-converter) (`gnomAD-converter`) 
  <br>Transforms GNOMAD VCFs (`vcf.bgz` format) into a format compatible with 
  CRISPRme. The function supports VCFs from GNOMAD v3.1, v4.0, and v4.1, 
  including *joint* VCFs.

- [**Generate Personal Card**](#226-generate-personal-card) (`generate-personal-card`) 
  <br>Generates a personalized summary for a specific sample, identifying all 
  private off-targets unique to that individual.

- [**Setup Legacy Database**](#227-setup-legacy-database) (`setup`)
  <br>Initializes the CRISPRme legacy database by downloading reference genomes,
  variant datasets, PAM definition files, and all associated resources
  originally distributed through the CRISPRme web server. The downloaded
  resources can be reused across analyses, eliminating the need for repeated
  downloads and simplifying environment setup.

- [**Web Interface**](#228-local-web-interface) (`web-interface`)
  <br>Launches CRISPRme's interactive web interface, allowing users to manage 
  and execute tasks directly via a local browser.

#### 2.2.1 Complete Search
---

The **Complete Search** function performs an exhaustive variant- and 
haplotype-aware off-target analysis, leveraging the provided reference genome 
and variant datasets to deliver comprehensive results. This feature integrates 
all critical stages of the CRISPRme pipeline, encompassing off-target 
identification, functional annotation, and detailed reporting.

**Quickstart (batteries install).** After the Docker
[quickstart](docs/DOCKER_QUICKSTART.md) download, a variant-aware search is a
single command — `download --what index` already wrote the `list_vcf.txt` /
`list_samplesID.txt` the search reads, and the guide is auto-padded to the PAM:

```bash
printf '%s\n' CTAACAGTTGCTTTTATCACNNN > my_guide.txt   # 20-nt spacer + NNN for the PAM positions
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:latest \
  crisprme.py complete-search \
    --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
    --guide my_guide.txt --vcf list_vcf.txt --samplesID list_samplesID.txt \
    --annotation Annotations/dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz \
    --gene_annotation Annotations/gencode.protein_coding.bed.gz \
    --mm 4 --bDNA 1 --bRNA 1 --output my_search --thread 4
# ^ the recommended default: 4 mismatches + 1 DNA + 1 RNA bulge = a 6-edit budget, matching
#   the default --max-total-edits 6, so nothing (incl. SNP+indel co-occurrences) is pruned.
#   Raise the per-type caps for a deeper search, but raise --max-total-edits to their sum too.
# shareable report: crisprme.py generate-report --result-dir Results/my_search
```

The example below shows **every** argument (a classic setup with a custom VCF
config file and base-editing options); see the
[CLI guide](docs/crisprme_data_setup_051826.md) for each one.

Usage Example for the Complete Search function:
  ```bash
  crisprme.py complete-search \
    --genome Genomes/hg38 \
    --vcf vcf_config.1000G.HGDP.txt \
    --guide sg1617.txt \
    --pam PAMs/20bp-NRG-SpCas9.txt \
    --annotation Annotations/dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz \
    --gene_annotation Annotations/gencode.protein_coding.bed.gz \
    --samplesID samplesIDs.1000G.HGDP.txt \
    --be-window 4,8 \
    --be-base A,G \
    --mm 4 \
    --bDNA 1 \
    --bRNA 1 \
    --merge 3 \
    --sorting-criteria-scoring mm+bulges \
    --sorting-criteria mm,bulges \
    --output sg1617-NRG-1000G-HGDP \
    --thread 8
  ```
  (Each argument is documented under **Input Arguments** below.)

##### Input Arguments
---

Below is a detailed list of the input arguments required or optionally used by 
the Complete Search function. Each parameter is explained to ensure clarity in 
its purpose and usage:

**General Parameters**
- `--help`
  <br>Displays the help message with usage details and exits. Useful for quickly 
  referencing all available options.

- `--output` (*Required*)
  <br>Specifies the name of the output directory where all results from the 
  analysis will be saved. This directory will be created within the `Results` 
  directory.

- `--thread` (*Optional - Default: 8*)
  <br>Defines the number of CPU threads to use for parallel computation. 
  Increasing the number of threads can speed up analysis on systems with 
  multiple cores.

- `--debug` (*Optional*)
  <br>Runs the tool in debug mode. 

**Input Data Parameters** 

- `--genome` (*Required*)
  <br>Path to the directory containing the reference genome in FASTA format. 
  Each chromosome must be in a separate file (e.g., `chr1.fa`, `chr2.fa`, etc.). 

- `--vcf` (*Optional*)
  <br>Path to text config file listing the directories containing VCF files to 
  be integrated into the analysis. When provided, CRISPRme conducts variant- 
  and haplotype-aware searches. If not specified, the tool searches only on the 
  reference genome. 

- `--guide` 
  <br>Path to a text file containing one or more guide RNA sequences (one per 
  line) to search for in the input genome and variants. This argument cannot be
  used together with `--sequence`.

- `--sequence`
  <br>Path to a FASTA file listing guide RNA sequences. This argument is an 
  alternative to `--guide` and cannot be used simultaneously.

- `--pam` (*Required*)
  <br>Path to a text file specifying the PAM sequence(s) required for the 
  search. The file should define the PAM format (e.g., `NGG` for SpCas9).

**Annotation Parameters (*Optional*)**

- `--annotation`
  <br>Path to a BED file containing genomic annotations, such as regulatory 
  regions (e.g., DNase hypersensitive sites, enhancers, promoters). These 
  annotations provide functional context for identified off-targets.

- `--gene_annotation`
  <br>Path to a BED file of gene information (e.g. GENCODE protein-coding
  annotations). Used for a **nearest-feature** lookup (not an overlap): it adds the
  closest gene and its signed distance to every scored off-target — reported as
  `Annotation_closest_gene_name`, `Annotation_closest_gene_ID` and
  `Annotation_closest_gene_distance_(kb)` in the full results (curated report tables
  show `Gene` + `Gene_distance_kb`; distance `0.0` = overlapping, negative = the
  other side of the feature). Available for both `complete-search` **and**
  `assembly-search` (where sites with no hg38 coordinate get no closest gene).
  Independent of `--annotation`; needs BEDOPS (`closest-features`/`sort-bed`).

**Base Editing Parameters (*Optional*)**

- `--be-window`
  <br>Specifies the editing window for base editors, defined as start and stop 
  positions relative to the guide RNA (1-based, comma-separated). This defines 
  the region of interest for base-editing analysis.

- `--be-base`
  <br>Defines the target nucleotide(s) for base editing. This is only used when 
  base editing functionality is needed.

**Sample-Specific Parameters (*Optional*)**

- `--samplesID` 
  <br>Path to a text config file listing sample identifiers (one per line)
  corresponding to VCF datasets. This enables sample-specific off-target 
  analyses. Mandatory if `--vcf` is specified.

**Search and Merging Parameters**

- `--mm` (*Required*)
  <br>Maximum number of mismatches allowed during off-target identification.

- `--bDNA` (*Optional*, default `1`)
  <br>Maximum allowable DNA bulge size. If omitted (together with `--bRNA`), it defaults to
  **1** (one DNA bulge), bounded by the bulge depth the installed index supports — the
  recommended single-bulge default (the CRISPR-Bulge ML scorer's validated domain). Pass it
  explicitly for a deeper search. If no bulge-capable index is installed the search stays
  bulge-free.

- `--bRNA` (*Optional*, default `1`)
  <br>Maximum allowable RNA bulge size. See `--bDNA` — defaults to one RNA bulge when neither
  is given.

- `--max-total-edits` (*Optional*, default `6`)
  <br>Cap on the combined number of edits (mismatches + DNA/RNA bulges) considered
  per candidate off-target, and — when `--bDNA/--bRNA` are omitted — the budget the
  per-type bulge caps are derived from. Lower values speed up dense-variant searches;
  raise it to recover targets with many simultaneous edits. A **positive** value is
  required: set it to `mm + bDNA + bRNA` (the sum of the per-type caps) to keep every
  requested-depth alignment; **negative values are rejected**. (This is a different
  knob from the `CRISPRME_IUPAC_CAP` environment variable, where `-1` disables the
  per-window greedy IUPAC cap — see the high-variant-density section below.)
  <br>NOTE: for a **variant-created** off-target you must also pass `--vcf` (or use a
  precomputed variant index); a reference-only search will not surface it regardless of
  this value.

- `--per-sample` (*Optional*)
  <br>Force per-sample genotype resolution into observed haplotypes for the SNP
  post-analysis. Reports each observed haplotype with CONFIRMED cis phasing, named carrier
  samples, and exact joint allele frequency. On a genotyped+PHASED index (HGDP + 1000
  Genomes (the default), or HPRC) this is already the context-aware default; on
  a sites-only index (`mega`) it cannot resolve carriers (no genotypes) and
  `complete-search` warns. Can be slow / intractable on dense or aggregate panels. To opt
  out of per-sample on a phased index use **`--population-level`**, which runs one
  worst-possible representative per variant window (exact worst-case CFD, dataset-level
  allele frequencies), lossless for detection and tractable on any panel. It replaces the
  removed `--full` flag; `--per-sample` is a genotype-resolution dial, not a speed dial
  (see METHODS §5).

- `--compute-backend {cpu,cuda,metal}` (*Optional*)
  <br>ML-scorer compute device: `cpu` (default) | `cuda` (NVIDIA/Linux GPU) | `metal`
  (Apple Silicon, self-tests and falls back to CPU); degrades to CFD-only if the cbulge
  scorer env is absent. In Docker, the default `pinellolab/crisprme:latest` image is
  GPU-capable — add `--gpus all` to `docker run` on an NVIDIA host and pass
  `--compute-backend cuda`; it falls back to CPU when no GPU is present. See
  [GPU acceleration](#2210-reference-index-and-data-distribution-commands) for details.

> **Scratch space (`$CRISPRME_TMPDIR`).** A genome-wide run needs tens of GB of
> temporary scratch. On HPC nodes and in containers `/tmp` is often a small
> RAM-backed `tmpfs`, which can make a run fail with a `No space left` / temp-file
> error. Set `CRISPRME_TMPDIR` to a roomy, writable directory on your data disk and
> CRISPRme routes **all** temp there (its own + the `sort-bed`/`bcftools` children):
> `export CRISPRME_TMPDIR=/big/disk/crisprme_tmp` (or, in Docker,
> `-e CRISPRME_TMPDIR=/DATA/crisprme_tmp`). It is created if missing. See
> [Troubleshooting](docs/TROUBLESHOOTING.md) for this and other environment issues
> (disk, memory, `/tmp`, GPU, degenerate-PAM+bulges).

- `--alt-alignments` / `--no-alt-alignments` (*Optional*)
  <br>Emit (or skip) the exhaustive alternative-alignments file
  (`..._all_results_with_alternative_alignments.tsv` — the non-best alignments at each
  locus). The default follows the analysis mode: **off** for the default
  population-level run (the report is built from the best alignment per locus, so this
  is a large, slow-to-build companion that nothing in the report reads), **on** under
  `--per-sample`. (`assembly-search` reconciles haplotypes from `integrated_results.tsv`
  only, so it does not need this file either.) Suppressing it does not change the report
  or any reported number.

- `--merge` (*Optional - Default: 3*)
  <br>Defines the window size (in base pairs) used to merge closely spaced 
  off-targets. Pivot targets are selected based on the highest score (e.g., 
  CFD, CRISPR-Bulge) or criteria defined by the `--sorting-criteria`. 

- `--sorting-criteria-scoring` (*Optional - Default: `mm+bulges`*)
  <br>Specifies sorting criteria for merging when using CFD/CRISPR-Bulge scores. 
  Options include:
    - mm: Number of mismatches.
    - bulges: Total bulge size.
    - mm+bulges: Combined mismatches and bulges. 

- `--sorting-criteria` (*Optional - Default: `mm+bulges,mm`*)
  <br>Sorting criteria used when CFD/CRISPR-Bulge scores are unavailable. Options are 
  similar to `--sorting-criteria-scoring` but tailored for simpler analyses. 

- `--keep-intermediates` (*Optional*)
  <br>Keep the bulky pre-integration working files. By default, once the finalized
  `*.integrated_results.tsv` and the report are produced, `complete-search` removes the
  intermediate tables it no longer needs — the merged non-best `*.altMerge.txt`,
  the raw `*_all_results_with_alternative_alignments.tsv` (the compressed `.zip` is kept),
  and the per-chromosome `*.population_summary.tsv` — reclaiming on the order of ~15 GB on a
  genome-wide run. **Always kept** (regardless of this flag): `*.integrated_results.tsv` and
  its `.zip`, the `*.bestMerge.txt` best-alignments table (the web results page greps it live
  for the per-position cluster drill-down, so it must persist for the job to stay browsable),
  the `*_report.zip`, the alternative-alignments `.zip`, the per-chromosome co-occurrence /
  high-variant-density sidecars (so a later `generate-report --report-full` can rebuild those
  bundles), the per-guide/per-sample summaries, `imgs/`, and the SQLite `.<job>.db` web
  results-page index. Pass `--keep-intermediates` to retain everything (e.g. for debugging or
  custom downstream parsing).

**Note 1**: Ensure compatibility between input files and genome builds (e.g., 
hg38 or hg19) to avoid alignment issues.

**Note 2**: Optional arguments can be omitted when not applicable, but required 
arguments must always be specified.

##### Output Data Overview
---

The Complete Search function generates a comprehensive suite of reports, 
detailing the identified and prioritized targets, along with statistical and 
graphical summaries. These outputs are essential for interpreting the search 
results and understanding the impact of genetic diversity on CRISPR off-target.

**Output Off-targets Files Description**

1. `*.integrated_results.tsv`
    
    - **Contents**: A detailed file containing the top targets (`*.bestMerge.txt`) 
    enriched with annotations from the input files. Includes:
        - Gene proximity of off-targets — `Gene_distance_kb` plus a `Gene_region`
        column locating each off-target relative to the nearest protein-coding gene
        (`CDS`, `5'UTR`, `3'UTR`, `exon`, `intron`, or `intergenic`; `intergenic`
        iff `Gene_distance_kb > 0`)
        - Overlaps with regulatory elements
    
    - **Purpose**: Integrates functional genomic context into the prioritization 
    of off-targets.

2. `*.all_results_with_alternative_alignments.tsv` 
    
    - **Contents**: Comprehensive listing of all identified targets, including 
    alternative alignments. Annotated with:
        - Gene proximity
        - Overlaps with regulatory elements
    
    - **Purpose**: Facilitates a full exploration of CRISPRme's off-target 
    predictions and their functional relevance.
    - **Note**: emitted **by analysis mode** — off in population-level runs
    (`--population-level` on genotyped indexes, or the default on sites-only/unphased
    indexes), on in per-sample runs (the default on genotyped+phased indexes, or under
    `--per-sample`); force either way with `--alt-alignments` / `--no-alt-alignments`.
    Nothing in the report reads it.

**Co-occurrence companion outputs** (the default-on SNP+SNP / SNP+indel co-occurrence
feature; both are bundled into the report ZIP):

2a. `*.indel_snp_cooc.tsv`

    - **Contents**: One row per off-target requiring **both** a nearby SNP **and** an
    indel in the same protospacer window — chrom, indel pos/ref/alt, off-target start,
    strand, SNP positions + rsIDs, phase (`CONFIRMED`/`PUTATIVE`), joint AF, cis-carrier
    count, and carrier samples.

    - **Purpose**: Surfaces SNP+indel co-occurring off-targets, which the classic
    two-independent-passes search could not see.

2b. `*.snp_snp_cooc.tsv`

    - **Contents**: One row per off-target requiring **≥2** nearby SNP alt alleles
    together — Chromosome, Position, Direction, crRNA, DNA, SNP positions, rsIDs, Phase,
    `MinAF_bound`, N carriers, and carriers.

    - **Purpose**: Surfaces SNP+SNP co-occurring off-targets.

    Both companions are joinable to the `bestMerge` tables by
    `(Chromosome, Position, Direction, crRNA, DNA)`.

**Guide-Specific Summary Files**

These files summarize off-target statistics per guide sequence based on 
different sorting criteria.
    
3. `*.summary_by_guide.<guide-sequence>_CFD.txt` 
    
    - **Contents**: Summarizes off-target counts per guide using the CFD score
    as the primary sorting criterion (data derived from 
    `*.integrated_results.tsv`). Includes counts of:
        - Targets by bulge type (DNA, RNA).
        - Mismatch number and bulge size.
        - Targets in the reference genome, variant genome, and those caused by 
        PAM creation due to variants.

    - **Purpose**: Provides insight into the distribution and characteristics of 
    off-targets prioritized by CFD score.
    
4. `*.summary_by_guide.<guide-sequence>_CRISPR-Bulge.txt`

    - **Contents**: Summarizes off-target counts per guide using the CRISPR-Bulge score
    as the primary sorting criterion (data derived from 
    `*.integrated_results.tsv`). Includes counts of:
        - Targets by bulge type (DNA, RNA).
        - Mismatch number and bulge size.
        - Targets in the reference genome, variant genome, and those caused by 
        PAM creation due to variants.

    - **Purpose**: Provides insight into the distribution and characteristics of 
    off-targets prioritized by CRISPR-Bulge score.

5. `*.summary_by_guide.<guide-sequence>_fewest.txt`

    - **Contents**: Summarizes off-target counts per guide using the fewest 
    mismatches and bulges as the sorting criterion.

    - **Purpose**: Highlights off-targets that are closest to perfect matches, 
    providing an alternative prioritization method.

**Sample-Specific Summary Files**

These files focus on off-targets unique to individual samples and their 
populations.

6. `*.summary_by_samples.<guide-sequence>_CFD.txt`

    - **Contents**: Counts of private off-targets per sample, sorted by CFD 
    score. Reports targets:
        - Private to the sample.
        - Found in the population or superpopulation.
        - Resulting from PAM creation due to a variant.

    - **Purpose**: Quantifies sample-specific off-targets and their broader 
    population impact.

7. `*.summary_by_samples.<guide-sequence>_CRISPR-Bulge.txt`
    
    - **Contents**: Similar to the CFD-based sample summary but uses the CRISPR-Bulge 
    score for sorting.

8. `*.summary_by_samples.<guide-sequence>_fewest.txt`
    - **Contents**: Summarizes private off-targets using the fewest mismatches 
    and bulges as the sorting criterion.

**Graphical Output**

9. `imgs` directory 

    - **Contents**: Contains visual representations of the top 1000 targets 
      based on CFD score, CRISPR-Bulge score, and Fewest mismatches and bulges. Images 
      include:
        - Bar plots showing the distribution of targets across populations and 
        bulge types.
        - Graphical summaries illustrating the impact of genetic variants on 
        mismatches, bulge size, and scores.

    - **Purpose**: Facilitates easy interpretation and presentation of CRISPRme 
    results.

#### 2.2.2 Complete Test
---

The **Complete Test** module verifies the correct installation and functionality 
of CRISPRme by running the built-in example guide against a **prebuilt index that 
you have already downloaded** (or built). It does **not** build an index or 
auto-download all the data for you: CRISPRme never auto-builds an index, so if no 
prebuilt index is present, `complete-test` prints guidance to download one (e.g. 
`crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`) and 
exits cleanly.

Once an index is available, `complete-test` runs the example search end-to-end 
and produces the same set of outputs as Complete Search, confirming that the 
installed tool works as expected. You can limit the test to a specific chromosome 
(e.g., chromosome 22) with `--chrom`, significantly reducing runtime and resource 
usage, or omit `--chrom` to run across the entire genome.

This module is suited for:

- **Initial Installation Verification**: Test whether CRISPRme has been 
installed correctly and is functioning as expected before running actual 
analyses. 

- **System Configuration Testing**: Validate the compatibility of CRISPRme with 
different system configurations (e.g., varying thread counts, datasets).

- **Dataset Evaluation**: Confirm that the variant-aware index you downloaded
(e.g. the default `NRG_3_hg38+hg38_HGDP1kGP`) searches correctly end-to-end.

Usage Example for the Complete Test function (run after downloading a prebuilt
index, e.g. `crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`):
  ```bash
  crisprme.py complete-test \
    --chrom chr22
  ```

##### Input Arguments
---

Below is a detailed list of the input arguments required by the Complete Test
function, including detailed explanations and default behaviors:

**General Parameters**

- `--help`
  <br>Displays the help message with details about the available options and 
  exits the program.

- `--thread` (*Optional - Default: 4*)
  <br>Specifies the number of threads to use for the computation, allowing for 
  parallel processing.

- `--debug` (*Optional*)
  <br>Runs the tool in debug mode. 

**Input Data Parameters (*Optional*)**

- `--chrom`
  <br>Specifies the chromosome to be used in the CRISPRme test. The chromosome 
  identifier must be prefixed with chr (e.g., chr22). When this argument is 
  provided, the test will be limited to the specified chromosome, enabling 
  faster execution for targeted validation. *Default behavior*: Run the test 
  across the entire genome.

> **Note:** `complete-test` requires a prebuilt index to already be present.
> CRISPRme never auto-builds an index; download one first with
> `crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`.
> If no index is found, `complete-test` prints this guidance and exits cleanly.

> **Two modes.** With no `--vcf_dataset`, `complete-test` runs the **modern smoke** above —
> the example guide against whichever prebuilt index you downloaded (recommended: the default
> `NRG_3_hg38+hg38_HGDP1kGP`). There is also a **legacy build-from-raw** mode,
> `complete-test --vcf_dataset {1000G|HGDP|1000G+HGDP}`, which downloads raw VCFs and builds a
> small test index on the fly; it exists mainly to feed the brute-force `validate-test` gate
> (see [§2.2.3](#223-validate-test)), which is supported only for a **1000 Genomes** run.

##### Output Data Overview
---

The results generated by the Complete Test function are saved in a subdirectory 
within the Results directory, named `crisprme-test-out`.

The output folder contains the same set of files produced by the 
Complete Search functionality, including:

- Detailed target reports that prioritize off-target candidates based on 
various criteria.

- Summaries categorized by guides and samples.

- Graphical representations of results, such as bar plots and impact 
assessments of genetic variants on target metrics.

For a detailed description of the contents of these files, please refer to the 
corresponding subsection in [Complete Search](#221-complete-search).

Finally, the off-target sites generated by the Complete Test workflow can be
validated using the Off-Target Sites Validation Test functionality, which
compares CRISPRme predictions against brute-force ground-truth results.
For additional details on this validation step, refer to the
[Off-Target Sites Validation Test section](#223-off-target-sites-validation-test)
in the documentation.

#### 2.2.3 Off-Target Sites Validation Test
---

The Validate Test functionality provides a systematic validation of
off-target sites generated by the Complete Test workflow in CRISPRme.
This module evaluates the correctness and completeness of CRISPRme’s
variant-aware off-target predictions by comparing them against results obtained
through a brute-force search and alignment strategy, serving as a ground-truth
reference.

The validation step is designed to ensure robustness and reliability of the
pipeline, particularly when working with large-scale population variant data.

> ⚠️ **Important Requirements**
>
> * The `complete-test` functionality **must be executed beforehand**
> * Validation is **supported only** when `complete-test` was run using
>   **1000 Genomes Project variant data**

Usage Example for the Off-Target Sites Validation Test function:
  ```bash
  crisprme.py validate-test \
    --chrom chr22
  ```

If no chromosome is specified, validation is performed across all chromosomes
available from the prior `complete-test` execution.

##### Input Arguments
---

Below is a detailed list of the input arguments supported by the Off-Target Sites 
Validation Test function, including their purpose and default behavior:

**General Parameters**

- `--help`
  <br>Displays the help message with details about the available options and
  exits the program.

- `--debug` (*Optional*)
  <br>Runs the tool in debug mode.

**Validation Parameters**

- `--chrom` (*Optional*)
  <br>Restricts validation to a specific chromosome (e.g., chr22).
  If not provided, validation is performed across all chromosomes processed
  during the `complete-test` run.

##### Output Data Overview
---

The **Off-Target Sites Validation Test** functionality does not generate new 
result files. Instead, it performs an **in-place validation** of the off-target 
sites produced by the **Complete Test** workflow and reports the outcome 
directly through standard output and error streams.

At completion, the Off-Target Sites Validation Test produces a validation 
summary that includes:

- Confirmation of a successful validation when the CRISPRme and brute-force
  off-target sites match exactly.
- Detailed diagnostics in case of discrepancies, including:
  - Total number of off-target sites in each set
  - Number of sites missing from CRISPRme predictions
  - Number of extra sites detected only by CRISPRme
  - Example mismatched sites to aid debugging

#### 2.2.4 Targets Integration
---

The **Targets Integration** function enables the seamless combination of 
computationally predicted off-targets identified by CRISPRme with 
experimentally validated off-targets, such as those obtained from GUIDE-seq, 
CIRCLE-seq, or other high-throughput methods. This integration enhances the 
interpretability and reliability of CRISPRme’s results by merging empirical 
data with *in silico* predictions.

This module supports integration with user-provided datasets in `BED` format, 
enabling flexibility in validation sources. Targets Integration 

Usage Example for the Targets Integration function:

> **Note:** the `--output` folder must already exist — `targets-integration`
> errors out if it does not. Create it first with `mkdir -p`.

  ```bash
  mkdir -p integrated_targets_dir
  crisprme.py targets-integration \
    --targets results.integrated_results.tsv \
    --empirical_data empirical_data.bed \
    --output integrated_targets_dir
  ```

##### Input Arguments
---

Below is a detailed list of the input arguments required by the Target 
Integration function, including detailed explanations and default behaviors:

**General Parameters**

- `--help`
  <br>Displays the help message with details about the available options and 
  exits the program.

**Input Data Parameters**

- `--targets` (*Required*)
  <br>Specifies the file containing targets identified and processed from a 
  CRISPRme search. This file should include predicted off-target data such as 
  mismatch counts, bulge sizes, and scores (e.g., CFD, CRISPR-Bulge).

- `--empirical_data` (*Required*)
  <br>Path to a `BED` file containing empirically validated off-targets. This 
  file typically includes genomic coordinates and additional metadata from 
  experiments such as GUIDE-seq or CIRCLE-seq. Ensure that the `BED` file 
  format adheres to standard conventions for compatibility.

- `--output` (*Required*)
  <br>Name of the directory where the resulting integrated targets file will be 
  saved. 

##### Output Data Overview
---

The Targets Integration module generates output files that merge 
computationally predicted targets from CRISPRme with experimentally validated 
off-targets provided in the BED file. The integrated data is stored in the 
specified output directory and includes the `integrated_results.tsv`. This 
tab-separated file contains the combined list of off-target sites. Each entry 
represents a merged record from CRISPRme predictions and empirical data, where 
overlaps are identified based on genomic coordinates.

#### 2.2.5 GNOMAD Converter
--- 

The **GNOMAD Converter** function is a utility designed to preprocess and 
convert GNOMAD VCF files into a format compatible with CRISPRme for off-target 
analysis. This tool facilitates the inclusion of population-level genetic 
variation data from GNOMAD into CRISPRme.

The converter currently supports VCF files from the following GNOMAD versions:
- **v3.1**
- **v4.0**
- **v4.1**, including joint VCF files (exomes + genomes).

Since individual sample data are not available in GNOMAD VCFs, the tool treats each
GNOMAD **population as a pseudo-individual sample**. This supports population-based
statistical analyses in CRISPRme (such as identifying population-specific off-targets),
but is not suitable for applications requiring individual-level granularity.

**Note 1**: For studies requiring sample-specific statistics, GNOMAD is not 
recommended due to its population-level nature.

**Note 2**: The GNOMAD Converter function is particularly useful for creating 
GNOMAD-compatible datasets formatted for CRISPRme's population-aware off-target 
analysis.

To ensure a smooth conversion process, sample ID files compatible with GNOMAD 
VCFs must be provided. These files are available for download from the CRISPRme 
GitHub repository:

- [Sample IDs file for GNOMAD v3.1 and v4.0](https://github.com/pinellolab/crisprme-plus/blob/main/test/data/samplesIDs/samplesIDs.gnomad.v40.txt)

- [Sample IDs file for GNOMAD v4.1](https://github.com/pinellolab/crisprme-plus/blob/main/test/data/samplesIDs/samplesIDs.gnomad.v41.txt)

The conversion process preserves all variant information necessary for CRISPRme
analyses, including allele frequencies and genotypes (if applicable).

Usage Example for the GNOMAD Converter function:
  ```bash
  crisprme.py gnomAD-converter \
    --gnomAD_VCFdir gnomad_vcf_dir \
    --samplesID samplesIDs.gnomad.v41.txt \
    --keep \
    --threads 4
  ```

##### Input Arguments
---

Below is a detailed list of the input arguments required by the GNOMAD 
Converter function, including detailed explanations and default behaviors:

**General Parameters**

- `--help`
  <br>Displays the help message with details about the available options and 
  exits the program.

- `--threads` (*Optional - Default: 8*)
  <br>Specifies the number of threads to use for the conversion, allowing for 
  parallel processing.

- `--debug` (*Optional*)
  <br>Runs the tool in debug mode. 

**Input Data Parameters**

- `--gnomAD_VCFdir` (*Required*)
  <br>Specifies the directory containing the gnomAD VCF files to be processed 
  and converted into a format compatible with CRISPRme.

- `--samplesID` (*Required*)
  <br>Path to a text file containing the sample IDs used during the conversion 
  process (the per-population sample lists linked above).

- `--joint` (*Optional*)
  <br>Use this flag if GNOMAD VCFs being processed are joint VCFs, such as 
  GNOMAD v4.1 joint variant files. *Default Behavior*: Assumes the input VCFs 
  are not joint.

- `--keep` (*Optional*)
  <br>Use this flag to retain all variants during the conversion, regardless of
  the FILTER field value. *Default behavior*: Excludes variants that do not have 
  a `PASS` value in the FILTER field.

- `--multiallelic` (*Optional*)
  <br>Indicates whether to merge variants at the same genomic position into a 
  single multiallelic record. *Default behavior*: Keeps variants as biallelic 
  and does not merge them.

##### Output Data Overview
---

The output of the GNOMAD Converter function consists of converted VCF files 
formatted to be compatible with CRISPRme. These files are stored in the same 
directory as the input VCF files. The conversion process ensures that the output 
adheres to CRISPRme's input specifications for population-level analysis. Below 
are details about the generated data:

**File Naming Conventions**

- **Multiallelic Variant Merging** (`--multiallelic`)
  <br>If multiallelic entries were generated by merging variants at the same 
  position, the filename will include a tag such as `*.multiallelic.*`, 
  `*.biallelic.*` otherwise.

- **Joint Variant Files** (`--joint`)
  <br>If joint VCFs were processed, the filenames will be labeled accordingly, 
  for example, `*.joint.*`.

**Content of the Converted VCFs**

- **Population Representation**
  <br>The pseudo-individual sample columns (one per GNOMAD population) are written
  into each output VCF, so CRISPRme can run its population-based statistical analysis.

- **Variant Quality**
  <br>If the `--keep` flag is used, all variants from the input VCF are included,
  regardless of their quality as indicated in the FILTER field. Without the 
  `--keep` flag, only variants with a `PASS` in the FILTER field are retained.

- **Allele Representation**
  <br>By default, the converter preserves biallelic representation, creating one 
  row per variant. If the `--multiallelic` flag is used, variants at the same 
  position are merged into multiallelic entries.

- **Compatible Structure**
  <br>The output files are structured to align with CRISPRme's population-aware 
  off-target analysis, ensuring seamless integration into the tool's pipeline.

#### 2.2.6 Generate Personal Card
---

The **Generate Personal Card** functionality creates a **sample-specific** 
report, referred to as **personal card**. This report details all off-targets 
identified by CRISPRme for a given sample's unique genomic sequence. The report 
accounts for private variants (genetic differences specific to the sample not 
shared with other populations or individuals) and their potential impact on 
off-target editing outcomes.

This feature is particularly useful in scenarios where personalized 
gene-editing strategies are required, such as analyzing how private genetic 
variations influence the efficacy and safety of CRISPR-based interventions. The 
personal card provides the following key insights: 

- Identification of off-target sequences present exclusively in the sample due 
to private variants. This allows researchers to evaluate potential risks or 
opportunities unique to the individual’s genome. 

- Detailed information for each input guide, showing how private variants 
affect off-target sequences.

This functionality is a critical tool for advancing personalized medicine and 
precision genome editing, enabling researchers and clinicians to tailor 
CRISPR-based solutions to an individual’s unique genetic makeup.

Usage Example for the Generate Personal Card function:
  ```bash
  crisprme.py generate-personal-card \
    --result_dir Results/sg1617.4.1.1 \
    --guide_seq CTAACAGTTGCTTTTATCACNNN \
    --sample_id NA21129
  ```

##### Input Arguments
--- 

Below is a detailed list of the input arguments required by the Generate 
Personal Card function, including detailed explanations and default behaviors:

**General Parameters**

- `--help`
  <br>Displays the help message with details about the available options and 
  exits the program.

- `--debug` (*Optional*)
  <br>Runs the tool in debug mode.

**Input Data Parameters**

- `--result_dir` (*Required*)
  <br>Specifies the directory containing the CRISPRme search results. The tool 
  extracts the relevant targets from the reports available in this directory. 
  Ensure this path includes all necessary output files generated by CRISPRme 
  during the analysis (`*.integrated_results.*`).

- `--guide_seq` (*Required*)
  <br>The sgRNA sequence for which the sample-specific targets will be 
  extracted. This argument ensures that the generated personal card is tailored 
  to a specific guide of interest, enabling targeted analysis.

- `--sample_id` (*Required*)
  <br>The identifier of the sample for which the personal card will be created. 
  This ID corresponds to the unique genetic profile being analyzed. Ensure the 
  sample ID matches the format used in the input data to avoid discrepancies.

##### Output Data Overview
---

The Generate Personal Card functionality produces sample-specific outputs, 
allowing researchers to assess how private genetic variants influence CRISPR 
off-target activity. 

**Note 1**: All output files include the sample ID in their file names for easy 
identification and traceability (e.g., `*.<sample_id>.*`).

**Note 2**: Outputs are tailored to the specified guide sequence and sample ID, 
ensuring precise and personalized reporting.

**Output Off-targets Files description**

- `*.<sample_id>.<guide_seq>.private_targets.tsv`
  <br>A detailed table reporting all off-target sequences specific to the 
  sample. Extracted from the *.integrated_results.* file within the input 
  directory. The report is generated in the results directory specified via 
  `--result_dir`.

**Graphical Output**

- `imgs` directory
  <br>The function generates plots illustrating the effect of private genetic 
  variants on the sample-specific targets. Displays changes in the CFD and  
  CRISPR-Bulge scores, and number of Mismatches and Bulges highlighting off-target 
  risk influenced by genetic variants.

#### 2.2.7 Setup Legacy Database
---

The Setup Legacy Database functionality initializes the CRISPRme legacy
database by downloading all reference genomes, variant datasets, PAM
definition files, and associated resources originally distributed through the
CRISPRme web server.

This functionality provides a convenient way to reproduce the original CRISPRme
environment and ensures that all required resources are available locally
before running analyses. Once downloaded, the resources can be reused across
multiple analyses without requiring additional downloads, reducing setup time
and improving reproducibility.

This feature is particularly useful when deploying CRISPRme on new systems,
setting up shared computational environments, or preparing offline analysis
workflows. The setup process provides the following key benefits:

- Automated download and organization of all required reference resources,
eliminating the need for manual retrieval and configuration.
- Support for downloading either the complete database or chromosome-specific
resources, allowing users to reduce storage requirements when working on
targeted analyses.
- Reuse of downloaded resources across multiple CRISPRme runs, avoiding
redundant downloads and ensuring consistent analysis environments.

This functionality serves as the foundation for CRISPRme analyses by providing
all reference data and auxiliary resources required by downstream workflows.

Usage Example for the Setup Legacy Database function:

```bash
crisprme.py setup \
  --path CRISPRme_database
```

##### Input Arguments
---

Below is a detailed list of the input arguments supported by the Setup Legacy
Database function, including detailed explanations and default behaviors:

**General Parameters**

- `--help`
  <br>Displays the help message with details about the available options and
  exits the program.
- `--debug` (Optional)
  <br>Runs the tool in debug mode.

**Database Configuration Parameters**

- `--path` (Optional)
  <br>Specifies the directory where the CRISPRme legacy database and all
  associated resources will be downloaded and installed. If not provided, the
  current working directory is used.
- `--chrom` (Optional)
  <br>Downloads resources associated with a specific chromosome only (for
  example, chr22). This option can significantly reduce download time and
  storage requirements when testing or performing chromosome-restricted
  analyses. By default, resources for all chromosomes are downloaded.
- `--force` (Optional)
  <br>Forces the download of all resources even if an existing CRISPRme
  database is already detected in the specified installation directory. This
  option is useful when repairing incomplete installations or refreshing
  previously downloaded data.

##### Output Data Overview
---

The Setup Legacy Database functionality creates a local CRISPRme database
containing all resources required for downstream analyses.

**Note 1**: The exact directory structure depends on the resources available
from the CRISPRme legacy repository and the selected download options.

**Note 2**: Downloaded resources are preserved and can be reused across
multiple CRISPRme analyses without requiring additional setup steps.

**Downloaded Resources**

- Reference genome files
  <br>Genome assemblies required by CRISPRme for off-target identification and
  genomic coordinate annotation.
- Variant datasets
  <br>Population-scale variant collections used to perform variant-aware
  off-target searches and population analyses.
- PAM definition files
  <br>Configuration files describing PAM sequences and nuclease-specific
  targeting rules supported by CRISPRme.
- Auxiliary resources
  <br>Additional metadata, indexes, and support files required by CRISPRme
  workflows and analysis modules.

**Output Directory**

- Installation directory specified by `--path`
  <br>Contains the complete CRISPRme legacy database organized in a structure
  compatible with all CRISPRme analysis workflows.

#### 2.2.8 Local Web Interface
---

The **Web Interface** module offers a user-friendly, locally hosted graphical 
user interface (GUI) for CRISPRme. This interface enables users to execute CRISPRme
workflows and explore results interactively through a web browser, without requiring
any external web service. This interface provides a complete interactive environment
for running CRISPRme analyses locally.

The GUI allows users to submit CRISPRme jobs directly through the web interface. 
Users can upload input files, configure parameters, and monitor progress with 
ease.

The interface provides an intuitive way to explore the output files generated by 
CRISPRme. Users can filter targets by criteria such as mismatch count, bulge 
size, or scores (e.g., CFD or CRISPR-Bulge). The interface includes dynamic plots and 
charts. Results are presented in a structured, easy-to-navigate format, linking 
data to relevant genomic annotations. The web interface runs as a local server, 
ensuring data privacy and fast response times. Users can access it via their 
preferred web browser, with compatibility confirmed for: Google Chrome, Mozilla 
Firefox, and Safari. All functionalities are self-contained, eliminating the 
need for an internet connection. This is particularly useful for secure
environments or systems without reliable internet access.

Usage example for the Web Interface function:

  ```bash
  crisprme.py web-interface
  ```

> 🌐 **Accessing the interface from another machine (SSH tunnel).** The server
> binds IPv4 `0.0.0.0:8080`. When forwarding it over SSH from a remote host, use
> the explicit IPv4 address as the forward target — `localhost` may resolve to
> IPv6 (`::1`) and be refused:
> ```bash
> ssh -N -L 8080:127.0.0.1:8080 user@remote-host   # then open http://localhost:8080
> ```

##### Input Arguments
---

Below is a detailed list of the input arguments required by the Web Interface 
function:

- `--help`
  <br>Displays the help message and exits.

- `--debug` (*Optional*)
  <br>Launches the local server in debug mode. This mode enables verbose logging, 
  which is useful for troubleshooting and diagnosing issues during setup or 
  operation. It provides detailed error messages and runtime information in the 
  console.

##### Output Data Overview
---

This function does not generate output data in the form of files or reports. 
Instead, it serves to launch a local graphical user interface, which allows 
users to interactively explore and run CRISPRme analyses. All results are 
displayed dynamically within the web interface itself, offering an interactive 
experience for viewing CRISPRme data and outputs (see [Section 2.2.1](#221-complete-search)
for details).

#### 2.2.9 Generate Report
---

The **Generate Report** function (`generate-report`) builds a self-contained,
easily shareable off-target assessment report for a completed CRISPRme run. It
produces a single ZIP (`<jobid>_report.zip`) in which **`report.html` is the only
top-level file — open it** — while every bundled data file sits under a **`data/`**
subfolder. `report.html` is fully offline (plots embedded as base64 PNGs, the
top-1000 table + CSS inline, no external dependencies), and every filename it
mentions links to its bundled file under `data/`. It is a portable digest of the
full interactive website result, aimed at sharing off-target predictions (for
example, to design a targeted-NGS / rhAMP-Seq confirmation panel).

**The report is SLIM by default.** It bundles the human-facing layer — the
top-1000 tables (both scorers), the ~100-site validation panel, the non-empty
per-tier curated CFD/edit-distance TSVs, and the co-occurrence companion tables —
which keeps the shareable ZIP small (single-digit MB on a typical run) and focused
on what a collaborator needs to seed a confirmation panel. The bulky, redundant
layers — the complete `integrated_results.tsv` (gzipped), the raw SNP+SNP
co-occurrence dump, the variant-created per-tier table, and the high-variant-density
BED — are **omitted by default**; pass **`--report-full`** to additionally bundle
them (gzipped) for a fully self-contained archive. The omitted data is never lost:
it stays in the result folder (and `complete-search` keeps the per-chromosome
sidecars so `--report-full` can rebuild those bundles even after an intermediate
cleanup). On the **web interface**, tick **"Include the full report"** under *Report
options* on the submission page for the same effect.

This report is generated **automatically at the end of every `complete-search`**
run, so most users never call it directly. Run it standalone only to
(re)generate the ZIP for an existing result folder.

Usage Example for the Generate Report function:
  ```bash
  crisprme.py generate-report \
    --result-dir Results/my-job
  ```

##### How to read the report
---

Unzip and open **`report.html`** (the only file at the top level; everything else
is under `data/`). Reading top to bottom:

- **1. Summary** — run parameters and an off-target-by-mismatch/bulge matrix. The
  **Analysis inputs & criteria** box's *Variants included* line states the
  genotyped panel size (e.g. **4,091 individuals** for HGDP + 1000 Genomes,
  `NRG_3_hg38+hg38_HGDP1kGP`) and the number of **SNPs** and **indels** searched.
- **⚠ Red banner ("no unambiguous on-target")** — appears when a guide has **more
  than one perfect (0-mismatch / 0-bulge) match**: there is then no *a-priori*
  on-target, so every perfect-match site is listed, flagged `Perfect_match = Yes`,
  and forced to the top of the validation panel. (One perfect match = a normal
  amber "presumed on-target" note instead.)
- **2. Recommended validation panel** — CFD / CRISPR-Bulge / edit-distance threshold
  tables plus a **hybrid worst-case ~100-site panel** (`data/panel_top100.tsv`) to
  seed a targeted-NGS / rhAMP-Seq confirmation assay. Among equally-severe
  candidates the panel **prioritizes sites observed in ≥1 real individual** (see the
  `Observed` column) — real haplotypes are preferred over worst-case reconstructions
  — without ever displacing a strictly-worse site, so worst-case coverage is kept.
- **`Observed` column** — whether an off-target is supported by **at least one real
  individual**:

  | Value | Meaning | When shown |
  |---|---|---|
  | `reference` | Present in the reference genome — carried by essentially every individual | Reference off-targets |
  | `N carrier(s)` | A variant site carried by N named individuals | Genotyped panel (never on a sites-only panel like `mega`, which has no carrier roster) |
  | `observed` | A single variant with non-zero allele frequency — carried by ≥1 individual by definition, even where the carrier cannot be named | Any variant-aware run, incl. sites-only panels |
  | `putative` | A multi-variant combination whose cis co-occurrence cannot be confirmed | Any variant-aware run |

  Under **`--per-sample`** the carrier list is the exact cis carriers, so `putative`
  cleanly marks reconstructions no one actually carries. On a sites-only panel (e.g.
  `mega`) the panel tie-break is inert but the column is still shown (values `observed`
  / `putative`); it is dropped only for a reference-only run.
- **`Gene_region` column** — where the off-target sits relative to the nearest
  protein-coding gene: `CDS`, `5'UTR`, `3'UTR`, `exon`, `intron`, or `intergenic`.
  It is coherent with `Gene_distance_kb` — `intergenic` iff the distance is `> 0`.
  Present in the curated report table, all download TSVs, and the raw
  `integrated_results.tsv`.
- **MAF column** — a value of **`1e-05` is a display floor** meaning "present but
  frequency effectively 0" (used so a zero-frequency allele still renders on the
  log-scale plots), **not** a measured frequency.
- **Key off-target table columns** — in the per-site table (and every download TSV):
  - **Chromosome / Position / Direction (`+`/`-`)** — where the off-target is and on which strand.
  - **crRNA vs DNA (aligned)** — the guide and the genomic sequence aligned; `-` marks a bulge.
  - **Mismatches**, **Bulge_type** (`DNA` / `RNA` / `-`) and **Bulge_Size** — the alignment's edit make-up.
  - **Variant / rsID** — the variant(s) that create a *variant* off-target (`chrom;pos;ref;alt`,
    plus rsID when known); empty for reference off-targets.
  - **Samples** — the named carrier individuals (populated under `--per-sample` on a genotyped index).
  - **CFD** — the primary specificity score. Designed for ≤3–4 mismatches; it **declines
    non-linearly** as mismatches grow, so do **not** compare CFD across searches that used
    different thresholds or variant datasets.
  - **CRISPR-Bulge score** — a complementary deep-learning activity score (Yaish & Orenstein,
    *NAR* 2024) on its **own** scale (bands **0.5 / 0.2 / 0.1**). It was trained on off-targets
    with at most one bulge, so any site that needs **≥ 2 bulges is out-of-domain and shows `-`**
    (N/A) — read CFD for those sites (a `-` means "not scored", not "no risk").
- **Off-target risk summary (top of the report).** The top of the report **replaces the legacy
  `100/(100+Σ)` aggregated specificity score** with a per-scorer **risk card**: for CFD and for
  CRISPR-Bulge separately (never blended), the count of ≤1-gap off-target sites in each score band and a
  **calibrated expected number of true off-targets** with a 90% range, plus a qualitative note that a
  deployed guide's own off-target propensity (unknown, not given a specific multiplier) can shift the
  count further. The estimates come from an isotonic calibration fit leave-one-guide-out on 359
  experimentally confirmed off-targets across 46 guides (held-out AUPR 0.168 [0.13–0.24] CFD /
  0.613 [0.56–0.70] CRISPR-Bulge; CRISPR-Bulge held out on 287/40). A "How to read this" pop-up explains
  the two scorers and shows the precision-recall curves. Full method, the four benchmark studies, and
  provenance (Zenodo 10.5281/zenodo.10902445) in `docs/offtarget_risk_summary_METHODS.md`.
- **`data/observed_haplotypes.tsv`** — for every site reached by **two or more distinct observed
  variant haplotypes**, this companion lists **each** haplotype with its variants, rsIDs, mm/bulges,
  **both** scores, origin, allele frequency and **named carriers**. The main table shows one
  representative per site (highest score, ties → most carriers); this companion is where you see the
  *full* set of real haplotypes a dense/multi-haplotype site contributes (it serves the CFD- and the
  CRISPR-Bulge-ranked views alike).
- **`--per-sample` reports only off-targets a real individual carries.** In a high-variant-density
  window CRISPRme does **not** fabricate a synthetic cross-sample "worst-case" haplotype under
  `--per-sample`; it reports the real observed haplotypes (exact carriers) and flags the dense region in
  `high_variant_density_regions.bed`. A population-level run (`--population-level`, or the default on a
  sites-only/unphased index) still reports the conservative worst-case representative, as it should for a screen.

##### Input Arguments
---

- `--result-dir` (*Required unless `--integrated-results` is given*)
  <br>A CRISPRme result folder (`Results/<jobid>`). The integrated-results TSV,
  `.Params.txt`, `.version.txt` and the default output location are resolved from
  here.

- `--integrated-results` (*Optional*)
  <br>An explicit `integrated_results` TSV (`.tsv` or `.tsv.gz`) to report on,
  instead of resolving it from `--result-dir`.

- `--samplesID-dir` (*Optional*)
  <br>Directory of samplesID files used for the superpopulation breakdown plot;
  auto-detected near the install when omitted.

- `--output` (*Optional*)
  <br>Output ZIP path. Defaults to `<result-dir>/<jobid>_report.zip`.

- `--report-full` (*Optional*)
  <br>Build the **full** (non-slim) report. In addition to the default human-facing
  layer, bundles the bulky redundant data — the complete `integrated_results.tsv.gz`,
  the raw `snp_snp_cooc.tsv.gz`, the `variant_created` per-tier table, and the
  `high_variant_density_regions.bed.gz` — so the ZIP is a fully self-contained
  archive. Omit it for the compact shareable default. `complete-search` accepts the
  same `--report-full` to make its auto-generated report full; the web submission
  page exposes it as the *Report options → "Include the full report"* checkbox.

#### 2.2.10 Reference-Index and Data-Distribution Commands

Bulge-enabled searches build a CRISPRitz index of the reference genome. `complete-search` does this automatically on the first run and reuses it afterward, but four commands let you control it explicitly — `build-index-only`, `download`, `publish-index`, and `scorer-env` — useful for batch jobs, shared clusters, and skipping the build via prebuilt indexes.

**`build-index-only`** — build the reusable reference index (into `genome_library/<PAM>_<bulges+1>_<genome>`) without running a search. Pass the same `--genome`/`--pam`/`--bDNA`/`--bRNA` you will search with:

```bash
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt --bDNA 2 --bRNA 2 --thread 16 --path /data/crisprme
```

This writes `genome_library/NRG_3_hg38` (the folder number is `max(bDNA,bRNA)+1`, so 2 bulges → `_3_`), which is the same index shipped precomputed on HuggingFace.

Add `--vcf` (a VCF dataset directory, e.g. `VCFs/hg38_HGDP1kGP`, **or** a config file listing one VCF folder name under `VCFs/` — the same `--vcf` form `complete-search` accepts) to also pre-build the **variant-aware** index (genome enrichment + SNP/indel indexes). Add `--samplesID` (a listing file, one samplesID filename per line under `samplesIDs/`; a combined panel lists each per-db samplesID file) to make it **self-complete**: with `--samplesID` the build ALSO emits the dict-less Tier-0 registry + Tier-1 genotype tiers and, for a merged panel, writes the combined `samplesIDs/<vcf>.samplesID.txt` list. Without `--samplesID` you get a dicts-only index (no fast post-analysis).

```bash
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt --bDNA 2 --bRNA 2 --vcf VCFs/hg38_HGDP1kGP --samplesID samplesIDs.config.txt --path /data/crisprme
```

To build a **different PAM** (NGG, NNGRRT, TTTN, …) on the shipped re-phased **HGDP+1kGP** or **mega** variant data — reusing the published enriched genome (fast) or the source VCFs — see [*Building a new-PAM index for a published dataset*](docs/PRECOMPUTED_INDEXES.md#building-a-new-pam-index-for-a-published-dataset-hgdp1kgp-mega).

**`complete-search --index-path <dir>`** — reuse a prebuilt/staged index library instead of building one under the working directory. A missing matching index is a hard error (rather than a silent rebuild), which is what you want on a read-only shared mount.

**`download`** — fetch reference data from a HuggingFace dataset repository over HF's CDN (typically much faster than the FTP/UCSC sources). `--what` selects the component: `genome`, `annotations`, `pams`, `samples`, `index` (with `--index-name`), `assembly` (with `--sample <HPRC id>`), `vcf` (with `--dataset` — raw per-dataset VCFs, only if the repo hosts a `vcfs/<dataset>/` folder; the default repo ships the standard panels as prebuilt `index`es instead, so this is for custom/self-hosted repos), or `all` (genome + annotations + PAMs + samples; on hg38 `all` also fetches the default `NRG_3_hg38` reference index so a reference scan works with no on-demand build). For `--what index`, `--no-genotypes` skips the big separate Tier-1 genotype store; the Tier-0 registry it keeps still corrects allele frequencies and rsIDs, but **per-individual `Samples` and variant off-target rows require the genotype tier** — omit `--no-genotypes` for full dict-less variant-aware output. Default repo `lucapinello/crisprme-data`, overridable with `--hf-repo` or `CRISPRME_HF_REPO`:

```bash
crisprme.py download --what all --path /data/crisprme
crisprme.py download --what index --index-name NRG_3_hg38 --path /data/crisprme
crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP --path /data/crisprme           # recommended default (phased HGDP + 1000 Genomes, 4,091 samples)
crisprme.py download --what assembly --sample HG01255 --path /data/crisprme                            # an HPRC individual's diploid assembly (for assembly-search)
```

**`publish-index`** — upload a locally built index to a HuggingFace dataset repository so other machines can skip the build (needs an HF write token via `--token` or `HF_TOKEN`). Add `--dictless` for a variant index to drop the ~152 GB per-sample SNP dicts (the registry + genotype tiers replace them; indel logs kept), upload the genotype store as a separate companion, and bundle the samplesID lists so the index is self-complete:

```bash
crisprme.py publish-index --index genome_library/NRG_3_hg38
crisprme.py publish-index --index genome_library/NRG_3_hg38+hg38_HGDP1kGP --dictless
```

**`scorer-env {create,check,update,list,doctor}`** — manages the dedicated `cbulge` conda env the CRISPR-Bulge scorer runs in (created automatically by `install_from_source.sh` unless deferred with `CRISPRME_SKIP_SCORER_ENV=1`):

```bash
crisprme.py scorer-env create
crisprme.py scorer-env doctor
```

The ML off-target scorer is **CRISPR-Bulge** (Yaish & Orenstein, *NAR* 2024) — a TensorFlow 2.13.x GRU ensemble — running isolated in the `cbulge` env so its heavy pins never touch the main environment. **CFD remains the primary score**; if the `cbulge` env is absent the run completes **CFD-only** (the ML column reads `-1`/N/A). **GPU acceleration:** the default image `pinellolab/crisprme:latest` is GPU-capable — on an NVIDIA host add `--gpus all` to `docker run` and pass [`--compute-backend cuda`](#input-arguments) to run the CRISPR-Bulge scorer on the GPU; it falls back to CPU automatically when no GPU is present. The amd64 image bundles CUDA; arm64 (Apple Silicon) and Docker-on-Mac are always CPU. The clickable installers auto-add `--gpus all` when Docker exposes a GPU. GPU pays off on large/dense searches; a tiny single-guide run is often faster on CPU due to init overhead. The env is **architecture-aware and scores identically across CPUs**: it installs the same TensorFlow 2.13.1 per-arch — conda `tensorflow-cpu`/CUDA on x86-64, `tensorflow-macos` on Apple Silicon, and the PyPI `tensorflow==2.13.1` aarch64 wheel on Linux ARM (validated x86-64 vs aarch64 agreement to max |Δ| = 3.6e-7, zero threshold-band changes), so the multi-arch Docker image scores identically on both. On networks that block `conda.anaconda.org`, point the env solve at a mirror with `CONDA_CHANNEL_BASE=https://prefix.dev` (or `CRISPRME_CONDA_CHANNEL_BASE`; the Docker build accepts `--build-arg CONDA_CHANNEL_BASE=…`).

See the companion data-setup guide (`docs/crisprme_data_setup_051826.md`, Sections 2d and 3.5) for the end-to-end workflow.

#### 2.2.11 Targets Scoring
---

The **Targets Scoring** function scores a user-supplied **list** of candidate
off-target sites — for example Cas-OFFinder output, or any table of aligned
guide/off-target pairs — **without running a genome search**. Each site is
annotated with the same two scores a normal CRISPRme+ report uses: the **CFD**
score and the **CRISPR-Bulge** ML score (the exact same calculations, so results
are directly comparable to a full CRISPRme+ analysis). This is useful for scoring
sites nominated by another tool — e.g. non-canonical-PAM hits — with CRISPRme+'s
calculations when building a dedicated index is not warranted. Bulge-containing
sites are scored with CFD's position-aware bulge penalty (as of v2.8.0; a
PAM-proximal bulge is penalized, not ignored) — supply the aligned pair in the
in-place convention (crRNA and off-target the same length, each bulge a `-` at its
position). CFD still cannot *rank* bulges; use the CRISPR-Bulge column for those
(see `docs/offtarget_risk_summary_METHODS.md` §3).

With `--index-path` pointing at an installed variant index, each site also gets a
**variant carrier overlay**: the variants overlapping its protospacer+PAM window
and, per population dataset, their allele frequency and carrier count — i.e. in how
many individuals the site is polymorphic (read from the index's
`Dictionaries/registry_<vcf>` tier).

Usage example (Cas-OFFinder output, with the carrier overlay):

  ```bash
  crisprme.py targets-scoring \
    --targets casoffinder_output.txt \
    --casoffinder \
    --index-path /path/to/index_install \
    --output scored_sites.tsv
  ```

For a generic TSV instead of Cas-OFFinder output, map the columns explicitly:

  ```bash
  crisprme.py targets-scoring \
    --targets sites.tsv \
    --columns guide=crRNA,dna=DNA,chrom=Chr,pos=Pos,strand=Strand,mm=Mismatches \
    --output scored_sites.tsv
  ```

##### Input Arguments
---

- `--targets` (*Required*) — the input site list (Cas-OFFinder output or a TSV).
- `--output` (*Required*) — output path (`.tsv`, or `.xlsx` with `--xlsx`).
- `--casoffinder` — parse the input as Cas-OFFinder output (schema auto-detected);
  omit and use `--columns` for a generic TSV.
- `--columns` — generic-TSV column map, e.g.
  `guide=crRNA,dna=DNA,chrom=Chr,pos=Pos,strand=Strand,mm=Mismatches` (at least
  `guide` and `dna` are required).
- `--index-path` — an installed variant index; adds the per-site carrier overlay.
- `--device` — CRISPR-Bulge scorer compute device: `cpu` (default) | `gpu`/`auto`
  (auto-pick `cuda`→`metal`→`cpu`) | `cuda` (NVIDIA/Linux GPU) | `metal` (Apple Silicon);
  same vocabulary as `complete-search --compute-backend`, and it falls back to CPU when
  the requested accelerator is unavailable.
- `--xlsx` — write an `.xlsx` workbook instead of TSV.
- `--pam-length` — PAM length in bp (default `3`).
- `--position-base` — `0`- or `1`-based input positions (default `0` for
  Cas-OFFinder, `1` for a generic TSV).

> **Note:** pure scoring (CFD + CRISPR-Bulge) needs no index; the carrier overlay
> requires an installed variant index via `--index-path`. Mismatch-only sites are
> scored exactly; bulge rows must follow CRISPRme's alignment convention (bulge
> columns as `-`).

## 3 Test

This section covers CRISPRme testing. This ensures that the software is 
correctly installed, configured, and ready for use. This section covers two 
testing options tailored to user needs:

- [Quick Test](#31-quick-test)
  <br>A lightweight test to verify the installation and basic functionality.

- [Detailed Test](#32-detailed-test) 
  <br>A comprehensive pipeline test, replicating a full-scale CRISPRme analysis.

For persistent issues, refer to the 
[CRISPRme GitHub Issues Page](https://github.com/pinellolab/crisprme-plus/issues) 
or contact the maintainers.

### 3.1 Quick Test
---

The Quick Test provides a simple and efficient way to confirm that CRISPRme is 
correctly installed and operational on your system. It ensures that the 
software is accessible, dependencies are properly configured, and the version 
matches the expected release.

**Step 1: Verify CRISPRme Installation**

Open a terminal and execute the following command to check the software version:

  ```bash
  crisprme.py --version
  ```

If the output displays the correct software version (e.g., `v2.8.0`), CRISPRme 
is successfully installed and ready for use.

**Step 2: Access CRISPRme Help Menu**

To explore the functionalities and input parameters of CRISPRme, use the 
`--help` flag:

  ```bash
  crisprme.py
  ```

The help menu provides detailed descriptions of CRISPRme's features and usage 
instructions, enabling users to familiarize themselves with available options 
and workflows.

**Why Perform the Quick Test?**

The Quick Test is highly recommended:

1. Immediately after installation to ensure the setup is correct.

2. Before running full-scale analyses to verify dependencies and environment 
configurations.

If the Quick Test completes without errors, you can confidently proceed to detailed analyses using CRISPRme’s powerful tools and features.

### 3.2 Detailed Test
---

For a more comprehensive validation of CRISPRme’s functionality, this detailed 
test offers a full-scale test of the main CRISPRme pipeline, specifically its 
Complete Search module. This test provides two distinct options for testing: a 
quicker, focused test on a single chromosome and a more detailed, exhaustive 
test across the entire genome:

- **Single Chromosome Test**
  <br>This test runs a search for potential off-targets on a single chromosome, 
  such as chromosome 22, using a prebuilt variant-aware index you have already 
  downloaded. This provides a fast way to check CRISPRme's ability to process 
  variant data and identify off-targets for a specific chromosome.

- **Full Genome Test**
  <br>This test checks CRISPRme's ability to search for potential off-targets 
  across the entire human genome using a prebuilt index. It runs a search using 
  the example guide RNA, incorporating both Gencode and ENCODE annotations to 
  ensure comprehensive results.

Upon successful completion of the Complete Test, users may optionally
perform an additional off-target validation step using the
Off-Target Sites Validation Test functionality. This step compares 
CRISPRme-predicted off-targets against brute-force ground-truth alignments 
derived from 1000 Genomes variant data, providing a rigorous correctness check 
of the pipeline.

**Why Perform the Detailed Test?**

The Detailed Test is ideal for users who wish to:

- Fully assess CRISPRme’s performance across the genome.
- Perform a more comprehensive check on CRISPRme’s ability to handle 
large-scale data and complex analyses.
- Validate CRISPRme’s off-target predictions against brute-force results using
the Validate Test functionality.

Successful completion of the detailed test, optionally followed by the
Validate Test, confirms the full functionality of CRISPRme, ensuring it is ready 
to handle large datasets and complex genetic analysis tasks.

#### 3.2.1 Single Chromosome Test
---

To run the quicker test on chromosome 22 against a prebuilt index you have
already downloaded (e.g.
`crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`),
execute the following commands:

  ```bash
  crisprme.py complete-test \
    --chrom chr22
  ```

After completion, off-target sites generated by this test can be validated using:
```bash
crisprme.py validate-test --chrom chr22
```

#### 3.2.2 Full Genome Test
---

To run the detailed test across the entire genome against a prebuilt index you
have already downloaded, execute the following commands:

  ```bash
  crisprme.py complete-test
  ```

Once the full genome test has completed, validation can be performed across all
chromosomes using:
```bash
crisprme.py validate-test
```

## 4 Citation
If you use CRISPRme in your research, please cite our [paper](https://rdcu.be/c1GYQ):

Cancellieri S, Zeng J, Lin LY, Tognon M, Nguyen MA, Lin J, Bombieri N, Maitland 
SA, Ciuculescu MF, Katta V, Tsai SQ, Armant M, Wolfe SA, Giugno R, Bauer DE, 
Pinello L. Human genetic diversity alters off-target outcomes of therapeutic 
gene editing. Nat Genet. 2023 Jan;55(1):34-43. 
doi: [10.1038/s41588-022-01257-y](https://doi.org/10.1038/s41588-022-01257-y).
Epub 2022 Dec 15. PMID: 36522432; PMCID: PMC10272994.

## 5 Contacts

- Luca Pinello
  <br>lpinello@mgh.harvard.edu

- Rosalba Giugno
  <br>rosalba.giugno@univr.it

- Daniel Bauer
  <br>bauer@bloodgroup.tch.harvard.edu

## 6 License

CRISPRme+ is licensed under the **MGB Open Access License 1.0 (MGBOA 1.0)**
([full terms](LICENSE)), which permits its use **only for non-commercial,
non-revenue-generating, academic purposes**.

### Commercial licensing

> **Draft summary — for guidance only.** The description below states the *intended*
> commercial-licensing model; the binding terms are those of the separate commercial
> license agreement. Licensing is administered with **Mass General Brigham
> Innovation** (which manages the underlying intellectual property). This text is not
> legal advice and is subject to their review.

**Any commercial or for-profit use requires a separate commercial license.** The
requirement is triggered by the commercial use of CRISPRme itself or of any result,
report, or other output it produces — for example in a clinical trial, diagnostic,
product, or drug-development program — **regardless of who ran the software.** In
particular, having an academic or non-profit party run CRISPRme on your behalf does
not confer commercial rights: the free academic license grants no commercial rights
to the outputs.

The model is meant to be fair and to keep CRISPRme widely usable:

| Scenario | Needs a commercial license? |
|---|---|
| Academic / non-profit user, for their own non-commercial, non-revenue-generating research and teaching | No — free under MGBOA 1.0 (no fee, no time limit) |
| Commercial / for-profit use of CRISPRme or of a result it produces (clinical trial, diagnostic, product, drug-development program) | Yes — covers **running** CRISPRme and the commercial use of what it produces. A result generated under a valid commercial license may then be used and transferred with the associated program (e.g. after an acquisition or out-licensing) with **no** repeat fee merely to *use* it |
| Service provider / CRO running CRISPRme to produce results for third parties | Yes — a commercial (service-provider) license |

**No warranty.** CRISPRme and its outputs are provided **"AS IS", without warranty
of any kind**; predictions are computational, may contain false positives and false
negatives, and are **not a substitute for experimental validation**. Outputs may
**change** as the software, its methods, or the underlying data are updated or
improved. To the maximum extent permitted by law, the authors, contributors, and
their institutions accept **no liability** for any use of, or reliance on, CRISPRme
or its outputs. (The same disclaimer appears in every generated report.)

For a commercial license or any licensing question, please contact **Luca Pinello**
at lpinello@mgh.harvard.edu.
