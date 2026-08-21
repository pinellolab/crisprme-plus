# TOPMed — Trans-Omics for Precision Medicine, aggregate (controlled-access) (`hg38_TOPMed`)

Aggregate allele-frequency data from the **NHLBI TOPMed** program, on **GRCh38 / hg38**. TOPMed is **controlled-access**: unlike 1000G / HGDP / gnomAD, these data cannot be downloaded from a public source. The files here were **provided directly to the authors** and normalised into a clean, CRISPRme-compatible form (a single aggregate pseudo-sample carrying `AC`/`AN`/`AF`/`HOM`), then filtered to MAF > 0.001. Because there is no open source, this directory has **no retrieval script** — only the rewrite script that documents how the distributed VCFs are organised.

---

## Summary

| Property | Value |
|---|---|
| Dataset | TOPMed (Trans-Omics for Precision Medicine), aggregate allele frequencies |
| Program | NHLBI TOPMed whole-genome sequencing (Taliun et al. 2021) |
| Reference build | GRCh38 / hg38 |
| Data type | Aggregate — single pseudo-sample; allele counts/frequencies, **no individual genotypes** |
| Access | **Controlled** — provided to the authors; not a public download |
| Distributed set | `more_than_001/` — variants with MAF > 0.001 (same threshold as gnomAD) |
| Variant types | SNVs and short INDELs (multiallelic ALTs preserved), unphased |
| File format | bgzip-compressed VCF (`VCFv4.2`) + tabix index |
| Files | `chr<chrom>.topmed.unphased.sorted.vcf.gz` (+ `.tbi`) |
| Chromosomes | chr1–chr22 and chrX (23 files) |
| INFO fields | `AF`, `AC`, `AN`, `HOM` |
| Sample column | one pseudo-sample, `TopMed` |
| Contig naming | `chr`-prefixed |
| Integrity | No public checksums (provided files) |

---

## Background

TOPMed (NHLBI) is a large, ancestrally diverse whole-genome sequencing program (Taliun et al. 2021 described 53,831 genomes). Individual-level TOPMed data are **dbGaP controlled-access**; NHLBI also exposes an aggregate allele-frequency resource through the **BRAVO** browser (sites and frequencies only, no individual genotypes). The data tracked here are aggregate allele frequencies obtained under controlled access — not a public release download.

---

## Data description

### Samples and populations

The distributed VCFs contain **no individual-level genotypes**. They carry a **single aggregate pseudo-sample** named `TopMed`, with population-scale allele counts and frequencies in `INFO`. Unlike gnomAD, these files are **not stratified by ancestry group** — there is one global aggregate only. (The underlying TOPMed cohort is large and diverse, but that structure is not exposed in these files.)

### Variant types

* **SNVs** and **short INDELs**; multiallelic ALT alleles are preserved (`REF` + comma-separated `ALT`).
* Genotypes are **unphased** (the single pseudo-sample's `GT` uses `/`).
* The distributed set is filtered to **MAF > 0.001** (see [Source](#source)), matching the threshold applied to gnomAD.

### Reference genome and coordinates

* **GRCh38 / hg38**, 1-based coordinates.
* **Contig names are `chr`-prefixed.** `rewrite_vcf.py` normalises contigs to the reference and **sorts** records by the reference's contig order, then position; the `##contig` lines (with lengths) are taken from the reference FASTA.

### Allele-frequency (INFO) fields

`rewrite_vcf.py` emits a deliberately minimal, clean header — exactly four `INFO` fields:

| Field | Number | Type | Description |
|---|---|---|---|
| `AF` | A | Float | Allele frequency (per ALT) |
| `AC` | A | Integer | Allele count (per ALT) |
| `AN` | 1 | Integer | Allele number |
| `HOM` | 1 | Integer | Homozygote count |

`FILTER` is always set to `PASS`; `FORMAT` is `GT` only; the single sample column is `TopMed`. `AC` / `AN` are clamped to the signed 32-bit range, and any missing `AF` / `AC` / `AN` / `HOM` defaults to 0.

### Record counts

Per-chromosome counts of the distributed MAF > 0.001:

| Chrom | Records | Chrom | Records | Chrom | Records |
|---|---|---|---|---|---|
| chr1 | 1,706,024 | chr9 | 951,207 | chr17 | 649,907 |
| chr2 | 1,655,778 | chr10 | 1,054,209 | chr18 | 597,573 |
| chr3 | 1,387,926 | chr11 | 1,021,062 | chr19 | 506,347 |
| chr4 | 1,371,425 | chr12 | 995,908 | chr20 | 547,544 |
| chr5 | 1,242,757 | chr13 | 780,957 | chr21 | 338,550 |
| chr6 | 1,236,542 | chr14 | 633,343 | chr22 | 365,069 |
| chr7 | 1,201,536 | chr15 | 617,555 | chrX | 575,796 |
| chr8 | 1,047,915 | chr16 | 657,660 | | |

### Sample and population metadata

These files use a single pseudo-sample, `TopMed` (no real individuals, no per-population breakdown). For a CRISPRme run, a matching `samplesID` file with a single `TopMed` row (with `SEX = n`) is required. Unlike gnomAD, no such file is bundled in the CRISPRme test data, so it must be supplied alongside the VCFs.

---

## Source

TOPMed is **controlled-access**, so there is no open download and **no retrieval script** in this directory:

* **Individual-level data** are hosted on **dbGaP** under controlled access.
* **Aggregate allele frequencies** are browsable via **BRAVO** ([bravo.sph.umich.edu](https://bravo.sph.umich.edu)) but are not offered as a bulk per-chromosome VCF download.

The VCFs here were **provided directly to the authors** and then normalised into the clean, CRISPRme-compatible form described above by **`rewrite_vcf.py`**:

```bash
python rewrite_vcf.py input.vcf.gz reference.fa output.vcf.gz
```

It reads the provided TOPMed VCF, keeps only `AF` / `AC` / `AN` / `HOM`, rebuilds
a clean `VCFv4.2` header with a single `TopMed` pseudo-sample, `chr`-prefixes and
sorts records against the reference, sets `FILTER=PASS`, and writes a
bgzip-compressed, tabix-indexed VCF. The distributed **`more_than_001/`** set is
this output restricted to **MAF > 0.001** (`INFO/AF > 0.001`), the same
allele-frequency threshold used for gnomAD.

## Reproducibility notes

* **Not reproducible by public download.** TOPMed access is controlled; the
  source data cannot be re-fetched from an open URL. Reproducing this dataset
  requires obtaining TOPMed access (dbGaP / the program) independently.
* **Deterministic transformation.** Given the same provided input VCF and
  reference, `rewrite_vcf.py` produces the same normalised output (fixed
  four-field INFO, single `TopMed` sample, reference-ordered sort); the
  distributed set is that output filtered to MAF > 0.001.
* **No public checksums.** As provided (non-public) files, there is no upstream
  MD5 manifest to verify against.

---

## Known limitations

* **Controlled access / not redistributable.** These files cannot be shared or
  re-downloaded openly; other researchers need their own TOPMed access.
* **Aggregate only.** A single `TopMed` pseudo-sample, with no individual
  genotypes and no ancestry stratification (contrast gnomAD's per-group AFs).
* **Freeze / version not recorded here.** The exact TOPMed freeze, and the exact
  chromosome list of the distributed files, are not independently verifiable from
  this repository — confirm against the provided data.
* **`HOM` may be uninformative.** If the source lacks homozygote counts, `HOM`
  defaults to 0.
* **Coverage gaps.** Like any short-read aggregate call set, reduced sensitivity
  in centromeric / segmental-duplication / low-complexity regions is expected.

---

## References

* Taliun D, Harris DN, Kessler MD, et al. *Sequencing of 53,831 diverse genomes
  from the NHLBI TOPMed Program.* Nature. 2021;590(7845):290–299.
* BRAVO variant browser (NHLBI / University of Michigan):
  [bravo.sph.umich.edu](https://bravo.sph.umich.edu)
