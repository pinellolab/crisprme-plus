# All of Us (AoU) — aggregate variant frequencies, scrape-built (`hg38_AoU`)

Aggregate allele-frequency data from the NIH **All of Us Research Program**, on **GRCh38 / hg38**. All of Us offers **no bulk VCF download**: the aggregate frequencies are only reachable through the Data Browser's public API. This dataset is therefore **built** by scraping that API per chromosome, assembling the results into an aggregate VCF with a single `AllOfUs` pseudo-sample, and filtering to MAF > 0.001. Only the filtered set is distributed (August 2026 snapshot).

---

## Summary

| Property | Value |
|---|---|
| Dataset | All of Us (AoU) Research Program, aggregate allele frequencies |
| Program | NIH All of Us (Choi et al. 2024) |
| Reference build | GRCh38 / hg38 |
| Data type | Aggregate, sites-only — single pseudo-sample, **no individual genotypes** |
| Access | Public Data Browser API (no authentication); no bulk file download |
| Distributed set | `more_than_001/` — variants with MAF > 0.001 (same threshold as gnomAD / TOPMed) |
| Variant types | SNVs and short INDELs; aggregate, no individual genotypes |
| File format | bgzip-compressed VCF (`VCFv4.2`) + tabix index |
| Files | `<chrom>.allofus.unphased.sorted.vcf.gz` (+ `.tbi`) |
| Chromosomes | chr1–chr22, chrX, chrY (24) |
| INFO fields | `AF`, `AC`, `AN`, `HOM` |
| Sample column | one pseudo-sample, `AllOfUs` |
| Contig naming | `chr`-prefixed |
| Build | scrape → CSV → VCF → MAF-filter (see [Source](#source)) |

---

## Background

*All of Us*, part of the NIH, is building a large, diverse US-based genomic health database. Traditional genomic databases are often biased toward European ancestry; *All of Us* enrolls more than half of its participants from racial and ethnic groups historically underrepresented in biomedical research. Its aggregate variant frequencies are exposed through the public Data Browser, but only via an interactive API (there is no downloadable per-chromosome VCF release).

---

## Data description

### Samples and populations

The dataset contains **no individual-level genotypes**. It carries a **single aggregate pseudo-sample** named `AllOfUs`, with population-scale allele counts and frequencies in `INFO`. The Data Browser API exposes only global aggregates, so, unlike gnomAD, these files are **not stratified by ancestry group**, despite the underlying cohort being large and diverse.

### Variant types

* **SNVs** and **short INDELs** (from the API's `variantType`).
* Each API variant is a single `REF`/`ALT` pair → one biallelic record; multiallelic sites appear as **separate biallelic records** at the same position.
* There are **no individual genotypes** — a single aggregate pseudo-sample carries `AC`/`AN`/`AF`/`HOM` only; its `GT` is synthesised from `AC`/`AN`/`HOM` and carries no phase information (phased vs unphased does not apply).
* The distributed set is filtered to **MAF > 0.001** (see [Source](#source)),
  matching the threshold applied to gnomAD.

### Reference genome and coordinates

* **GRCh38 / hg38**, 1-based coordinates.
* **Contig names are `chr`-prefixed.** `csv_to_vcf.py` prefixes the API's unprefixed chromosome labels and takes the `##contig` lines (with lengths) from the reference FASTA; records are then position-sorted with `bcftools sort`.

### Allele-frequency (INFO) fields

The built VCFs carry a minimal, clean header — exactly four `INFO` fields:

| Field | Number | Type | Description |
|---|---|---|---|
| `AF` | A | Float | Allele frequency |
| `AC` | A | Integer | Allele count |
| `AN` | 1 | Integer | Allele number |
| `HOM` | 1 | Integer | Homozygote count |

`FILTER` is always set to `PASS`; `FORMAT` is `GT` only; the single sample column is `AllOfUs`. `AN` reflects the AoU cohort and is roughly uniform across well-covered regions, dropping where coverage is incomplete (see [Known limitations](#known-limitations)). The pseudo-sample `GT` is synthesised from `AC`/`AN`/`HOM` (`1/1` only when every allele is ALT and homozygotes exist, `0/1` otherwise).

### Record counts

Per-chromosome counts of the distributed MAF > 0.001 set:

| Chrom | Records | Chrom | Records | Chrom | Records |
|---|---|---|---|---|---|
| chr1 | 2,765,787 | chr9 | 1,525,649 | chr17 | 1,074,713 |
| chr2 | 2,887,236 | chr10 | 1,759,200 | chr18 | 976,070 |
| chr3 | 2,305,415 | chr11 | 1,644,534 | chr19 | 866,083 |
| chr4 | 2,440,306 | chr12 | 1,643,868 | chr20 | 816,130 |
| chr5 | 2,198,202 | chr13 | 1,236,699 | chr21 | 540,514 |
| chr6 | 2,177,273 | chr14 | 1,101,401 | chr22 | 539,931 |
| chr7 | 2,045,924 | chr15 | 1,013,562 | chrX | 1,285,829 |
| chr8 | 1,853,860 | chr16 | 1,145,364 | chrY | 63,357 |

### Sample and population metadata

These files use a single pseudo-sample, `AllOfUs` (no real individuals, no per-population breakdown). For a CRISPRme run, a matching `samplesID` file with a single `AllOfUs` row (with `SEX = n`) is required. Unlike gnomAD, no such file is bundled in the CRISPRme test data, so it must be supplied alongside the VCFs.

### Provenance & population-aware search notes

- **Phasing:** aggregate (sites-only, no genotypes) — the single pseudo-sample `GT` is *synthesised* from `AC`/`AN`/`HOM`, not a real observed genotype, so "phased" vs "unphased" does not apply.
- **indel+SNP co-occurrence:** **not applicable** (aggregate, no per-sample genotypes) — CRISPRme cannot tag an indel+SNP haplotype as CONFIRMED or PUTATIVE here because there are no genotypes to establish co-occurrence; each record is treated as an independent site-level allele frequency.
- **Conventions:** MAF filter **`>0.001`** (`INFO/AF > 0.001`; the distributed `more_than_001/` set); contigs `chr`-prefixed.
- **Caveats:** the single pseudo-sample is named **`AllOfUs`** (matching the build script `csv_to_vcf.py`, `header.add_sample("AllOfUs")`, and verified in the distributed VCFs via `bcftools query -l`), so the accompanying `samplesID` must contain a single **`AllOfUs`** row (`SEX = n`) to match. The dataset is a **moving snapshot** (live API, August 2026) with no fixed checksum.

---

## Source

All of Us exposes its aggregate variant frequencies only through the **public Data Browser API** (no authentication, no bulk download):

```
https://public.api.researchallofus.org/v1/genomics/search-variants
```

Because there is no file release, the dataset is **built** by the pipeline in this directory, orchestrated by **`retrieve_AoU.sh`**:

```bash
REFERENCE=/path/to/hg38.fa ./retrieve_AoU.sh          # all 24 chromosomes
REFERENCE=/path/to/hg38.fa ./retrieve_AoU.sh chr22    # a subset
```

Per chromosome it runs four phases:

1. **Scrape** (`scrape_aou.py`) — walk the chromosome in 1 Mb windows, paginate the API, and save each response page as raw JSON. Resumable (checkpointed progress + a `.complete` marker) and rate-limited.
2. **Combine** (`json_to_csv.py`) — merge the JSON pages into one de-duplicated CSV (by `variantId`), keeping the API field names.
3. **Convert** (`csv_to_vcf.py`) — build an aggregate `VCFv4.2` with the single `AllOfUs` pseudo-sample and `AF`/`AC`/`AN`/`HOM` INFO, then `bcftools sort` + `tabix` → `full/` (an intermediate, **not distributed**).
4. **Filter** — `bcftools view -i 'INFO/AF > 0.001'` on the full VCF, writing the distributed **`more_than_001/`** set (bgzipped + indexed).

## Reproducibility notes

* **Built from a live API, not a fixed file.** Reproducibility depends on the Data Browser's current contents: the AoU aggregate grows with each data release, so a re-scrape at a later date may yield a different variant set. There is no fixed upstream release or checksum to pin against.
* **Deterministic downstream of the scrape.** Given the same scraped CSV and reference, the CSV → VCF → filter steps produce the same output (fixed four-field INFO, single `AllOfUs` sample, reference-ordered sort, `INFO/AF > 0.001`).
* **Only the filtered set is distributed.** The `full/` VCFs are an unshared build intermediate; the distributed dataset is `more_than_001/` only.

---

## Known limitations

* **Moving source.** The API is a live service, not a versioned file release; the returned set can change between AoU data releases (see Reproducibility notes).
* **Aggregate only.** A single `AllOfUs` pseudo-sample, with no individual genotypes and no ancestry stratification (contrast gnomAD's per-group AFs).
* **Coverage gaps in hard regions.** Large low-/zero-density stretches occur throughout the callset, concentrated in centromeric, pericentromeric and other highly repetitive, low-mappability regions where short-read WGS loses sensitivity (e.g. the chr2 centromere, 90–96 Mb). This is not AoU-specific.
* **Slow to build.** The scrape is paginated and politely rate-limited; a whole-genome build takes many hours — run inside tmux/screen.

---

## References

* Choi SH, Wang X, Rosenthal EA, et al. *Genomic data in the All of Us Research Program.* Nature. 2024;627(8003):340–346.
* All of Us Research Program. *Genomic Quality Report.* Data Browser Documentation (2025).
