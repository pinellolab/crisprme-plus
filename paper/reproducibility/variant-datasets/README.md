# Variant datasets — CRISPRme+ paper

This directory documents how the population variant datasets used in the **CRISPRme+** paper were retrieved (or built) and prepared for CRISPRme's population-aware off-target nomination. Each dataset lives in its own `hg38_<dataset>/` subfolder containing a **dataset-specific README** and, where the data are publicly retrievable, a **retrieval / build script**. All datasets are on **GRCh38 / hg38**.

---

## Datasets at a glance

| Folder | Dataset & release | Data type | Access | Retrieval / build |
|---|---|---|---|---|
| [`hg38_1000G_2019`](hg38_1000G_2019/) | 1000 Genomes Phase 3 — GRCh38 2019 biallelic (`20190312_biallelic_SNV_and_INDEL`) | Genotyped, phased — 2,548 samples | Open (IGSR/EBI) | `retrieve_1000G_2019.sh` — download + MD5 |
| [`hg38_1000G_2021`](hg38_1000G_2021/) | 1000 Genomes high-coverage NYGC (`20201028_3202_phased`) | Genotyped, phased — 3,202 samples, ~30× | Open (IGSR/EBI) | `retrieve_1000G_2021.sh` — download + MD5 |
| [`hg38_HGDP`](hg38_HGDP/) | Human Genome Diversity Project (`hgdp_wgs.20190516`) | Genotyped, high-coverage — 929 samples | Open (Wellcome Sanger) | `retrieve_HGDP.sh` — download + MD5 |
| [`hg38_gnomAD`](hg38_gnomAD/) | gnomAD v4.1 joint, sites-only | Aggregate allele frequencies — no genotypes | Open (Google Cloud) | `retrieve_gnomAD.sh` — download + convert + filter |
| [`hg38_TOPMed`](hg38_TOPMed/) | TOPMed aggregate | Aggregate — single pseudo-sample | **Controlled** — provided directly | `rewrite_vcf.py` — normalise (no download) |
| [`hg38_AoU`](hg38_AoU/) | All of Us aggregate (July 2025 snapshot) | Aggregate — single pseudo-sample | Public API — no bulk download | `retrieve_AoU.sh` — scrape + build |

---

## The datasets

### 1000 Genomes Project — 2019 biallelic reanalysis (`hg38_1000G_2019`)
1000 Genomes Project Phase 3 samples called natively on GRCh38 (March 2019 biallelic SNV+INDEL release). 2,548 phased samples across 26 populations / 5 super-populations; low-coverage source data. Openly downloaded from IGSR/EBI and MD5-verified. → [details](hg38_1000G_2019/README.md)

### 1000 Genomes Project — high-coverage NYGC (`hg38_1000G_2021`)
The ~30× re-sequencing of the expanded 1000 Genomes cohort (3,202 samples = 2,504 unrelated + 698 related forming 602 trios), phased with SHAPEIT2-duohm (chrX with Eagle2). Openly downloaded from IGSR/EBI and MD5-verified against the release's own manifest. → [details](hg38_1000G_2021/README.md)

### Human Genome Diversity Project (`hg38_HGDP`)
High-coverage WGS across 54 diverse populations (929 samples), GATK-called and VQSR-filtered. Openly downloaded from the Wellcome Sanger Institute and MD5-verified. → [details](hg38_HGDP/README.md)

### gnomAD v4.1 joint (`hg38_gnomAD`)
An aggregate, **sites-only** allele-frequency resource (730,947 exomes + 76,215 genomes; 10 ancestry groups; no individual genotypes). Openly downloaded from the gnomAD public Google Cloud bucket (MD5 from the object metadata), converted to a CRISPRme-compatible form with the built-in `gnomAD-converter`, and filtered to MAF > 0.001. → [details](hg38_gnomAD/README.md)

### TOPMed (`hg38_TOPMed`)
NHLBI TOPMed aggregate allele frequencies as a single `TopMed` pseudo-sample. TOPMed is **controlled-access**: these files were provided directly to the authors (not downloadable) and normalised into a clean VCF, then filtered to MAF > 0.001. Not reproducible by public download. → [details](hg38_TOPMed/README.md)

### All of Us (`hg38_AoU`)
NIH All of Us aggregate allele frequencies as a single `AllOfUs` pseudo-sample. All of Us has no bulk VCF release; the data are **built** by scraping the public Data Browser API, assembling per-chromosome VCFs, and filtering to MAF > 0.001. → [details](hg38_AoU/README.md)

---

## Retrieval modes

The datasets fall into three retrieval categories, reflected in each folder's script (or absence of one):

* **Open direct download** — 1000G 2019, 1000G 2021, HGDP, gnomAD. A `retrieve_<dataset>.sh` script downloads the files and **verifies every one by MD5** against a known-good manifest, resuming interrupted transfers and retrying transient failures.
* **Controlled access** — TOPMed. The data cannot be downloaded openly; they were provided directly and are normalised in place. Reproducing them requires independent TOPMed access.
* **Public API, no bulk download** — All of Us. The aggregate frequencies are scraped from the Data Browser API and assembled locally.

## Allele-frequency filtering (aggregate resources)

The three **aggregate** resources, namely gnomAD, TOPMed, All of Us, are filtered to **MAF > 0.001** (`bcftools view -i 'INFO/AF > 0.001'`). This common-variant subset is the variant of each aggregate dataset used in the paper. The genotyped panels (1000G, HGDP) are used unfiltered.

## Combined multi-source panel

  For the paper, the genotyped panels (1000G, HGDP), and the aggregate collections can be merged into a single multi-source panel so CRISPRme enriches the reference and scans **once** rather than per dataset. 

<!-- TODO: Add merged dataset description here  -->

## Common conventions

* **Format.** All variant data are bgzip-compressed VCFs, split by chromosome.
* **Reference / contigs.** GRCh38 / hg38 throughout. Contig naming is `chr`-prefixed in every dataset **except** 1000G 2019 (unprefixed, e.g. `22`)CRISPRme handles the matching internally.
* **Aggregate resources** carry a single pseudo-sample (gnomAD's ancestry groups; `TopMed`; `AllOfUs`) with `INFO` allele frequencies (`AF`/`AC`/`AN`[/`HOM`]) and no real genotypes (not available from source).
* **Sample metadata.** Each dataset needs a `samplesID` file mapping samples (or pseudo-samples) to populations for CRISPRme's population-aware statistics; see the per-dataset READMEs for where each comes from.

## Required tools

* `wget` (or `curl`) and one of `md5sum` / `md5` / `openssl` — download + verify.
* `bcftools` and HTSlib (`bgzip`, `tabix`) — sort / filter / index.
* `python3` with `pysam` (gnomAD-converter, TOPMed `rewrite_vcf.py`, AoU build) and, for the AoU scrape, `requests` + `pandas`.
* `crisprme.py` — for the gnomAD-converter (and downstream index building).
* A GRCh38 / hg38 reference FASTA — required to build the AoU VCFs.

## Reproducibility summary

* **Fixed-file, checksum-verified.** 1000G 2019, 1000G 2021, HGDP and gnomAD are pinned to specific releases and verified per file by MD5; a successful run yields byte-for-byte the same downloads on any machine.
* **Controlled, deterministic transform.** TOPMed is not re-downloadable, but the normalisation of the provided files is deterministic (no public checksums).
* **Live source, dated snapshot.** All of Us is scraped from a live API whose contents grow with each data release; the build is deterministic downstream of the scrape, but the scraped set is a point-in-time snapshot (dated in the dataset README), not a fixed upstream release.
