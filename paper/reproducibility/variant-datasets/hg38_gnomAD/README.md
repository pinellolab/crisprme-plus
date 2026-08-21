# Genome Aggregation Database (gnomAD) — v4.1 joint, sites-only (`hg38_gnomAD`)

Population **allele-frequency** data from **gnomAD v4.1 joint** (exome + genome
combined), aligned to **GRCh38 / hg38**. gnomAD is an **aggregate, sites-only**
resource: the VCFs carry population allele counts and frequencies in the `INFO`
column and contain **no individual-level genotypes**. CRISPRme+ uses it (via the
gnomAD-converter) for population-aware off-target frequencies. The retrieval
script downloads the raw v4.1 joint sites files and prepares them for CRISPRme by
converting them and allele-frequency-filtering to MAF > 0.001.

---

## Summary

| Property | Value |
|---|---|
| Dataset | Genome Aggregation Database (gnomAD), v4.1 joint (exome + genome) |
| Release | v4.1 (April 2024; files dated 2024-05-01) |
| Reference build | GRCh38 / hg38 (`assembly=gnomAD_GRCh38`) |
| Data type | Sites-only aggregate — allele frequencies, **no individual genotypes** |
| Aggregate cohort | 730,947 exomes + 76,215 genomes (807,162 individuals) |
| Ancestry groups | 10 (afr, ami, amr, asj, eas, fin, nfe, mid, sas, remaining) |
| Variant types | SNVs and short INDELs (multiallelic, not split) |
| File format | bgzip-compressed VCF (`.vcf.bgz`, `VCFv4.2`) |
| Files | `gnomad.joint.v4.1.sites.<chrom>.vcf.bgz` |
| Chromosomes | chr1–chr22, chrX, chrY (24 files) |
| Contig naming | `chr`-prefixed |
| Size | ~817 GB total (chr1 ~72 GB; smallest chrY ~0.8 GB) |
| Source | gnomAD public bucket on Google Cloud Storage |
| Access | Open / unrestricted |
| Integrity | Per-file MD5 checksums (GCS object metadata; see [Reproducibility notes](#reproducibility-notes)) |

---

## Background

gnomAD (Broad Institute) aggregates large-scale exome and genome sequencing data
to provide population allele frequencies for variant interpretation. **v4.1**
(released April 2024) combines **730,947 exomes and 76,215 genomes** (including UK
Biobank) into a single "joint" call set — larger and more diverse than earlier
gnomAD releases — and is called natively on GRCh38.

---

## Data description

### Samples and populations

gnomAD has **no individual samples** in the usual sense: it is an aggregate of
**730,947 exomes + 76,215 genomes (807,162 individuals)**, released as sites-only
allele-frequency summaries. Frequencies are stratified across **10 genetic
ancestry groups**:

* **afr** — African / African-American
* **ami** — Amish
* **amr** — Admixed American
* **asj** — Ashkenazi Jewish
* **eas** — East Asian
* **fin** — Finnish
* **nfe** — Non-Finnish European
* **mid** — Middle Eastern
* **sas** — South Asian
* **remaining** — individuals not assigned to a group above

CRISPRme treats each ancestry group as a pseudo-sample (see
[Sample and population metadata](#sample-and-population-metadata)).

### Variant types

* **SNVs** and **short INDELs**.
* **Sites-only**: each record has the 8 fixed VCF columns only
  (`#CHROM POS ID REF ALT QUAL FILTER INFO`) — no `FORMAT`/genotype columns
  (confirmed: `bcftools query -l` → 0 samples).
* **Multiallelic sites are not split** in the raw files; normalise with
  `bcftools norm -m -any` if strictly biallelic records are required.
* The `FILTER` column is meaningful and frequently non-`PASS` (e.g.
  `GENOMES_FILTERED`, `AC0`); choose a filtering policy accordingly.

### Reference genome and coordinates

* Called natively on **GRCh38** (`##contig` lines carry
  `assembly=gnomAD_GRCh38`). Positions follow the 1-based VCF convention.
* **Contig names are `chr`-prefixed** (`chr1`, …, `chrX`, `chrY`), matching the
  UCSC hg38 reference directly.

### Allele-frequency (INFO) fields

Allele frequencies are the core content. The primary field is **`AF_joint`** (the
combined exome + genome frequency), with matching `AC_joint` / `AN_joint`, and:

* per-dataset sub-fields `AF_exomes` / `AF_genomes` (with their `AC_*` / `AN_*`);
* per-ancestry-group fields `AF_joint_<grp>` (e.g. `AF_joint_afr`,
  `AF_joint_nfe`, `AF_joint_eas`);
* sex-stratified (`_XX` / `_XY`) and homozygote (`nhomalt_*`) breakdowns.

CRISPRme's gnomAD-converter reads `AF_joint` (its `--joint` mode) as the allele
frequency. For example, the record `chrY:2781489 C>T` carries `AC_joint`,
`AN_joint`, and `AF_joint` together with the per-dataset and per-ancestry
breakdowns.

### Record counts

Variant counts per chromosome in the final prepared set — after conversion and
the MAF filter (phase C, `INFO/AF > 0.001`):

| Chrom | Records | Chrom | Records | Chrom | Records |
|---|---|---|---|---|---|
| chr1 | 2,700,493 | chr9 | 1,513,459 | chr17 | 1,005,968 |
| chr2 | 2,889,343 | chr10 | 1,735,466 | chr18 | 954,406 |
| chr3 | 2,399,643 | chr11 | 1,652,134 | chr19 | 813,882 |
| chr4 | 2,441,038 | chr12 | 1,632,968 | chr20 | 779,641 |
| chr5 | 2,204,436 | chr13 | 1,221,538 | chr21 | 506,578 |
| chr6 | 2,178,072 | chr14 | 1,108,948 | chr22 | 511,366 |
| chr7 | 2,030,890 | chr15 | 1,017,179 | chrX | 1,354,225 |
| chr8 | 1,883,580 | chr16 | 1,097,964 | chrY | 36,598 |

### Sample and population metadata

Because gnomAD is aggregate-only, the "samples" CRISPRme uses are the ancestry
groups themselves, treated as pseudo-individuals. CRISPRme ships this mapping as
`samplesIDs.gnomad.v41.txt` (in the repository's `test/data/samplesIDs/`); a copy
is **bundled next to the retrieval script** in this folder, which the pipeline
uses automatically. It is a tab-separated file
(`#SAMPLE_ID  POPULATION_ID  SUPERPOPULATION_ID  SEX`) listing the 10 ancestry
groups above with `SEX = n` (not applicable). A separate
`samplesIDs.gnomad.v41.txt` covers gnomAD v3.1 / v4.1. These files are consumed
by the gnomAD-converter, which emits one pseudo-sample column per ancestry group.

---

## Source

The VCFs are served from the **gnomAD public bucket on Google Cloud Storage**:

```
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/
```

The data are openly accessible and require no credentials. The file dates
(2024-05-01) match gnomAD v4.1's April 2024 release, consistent with a direct
download rather than a local re-processing. The retrieval script downloads the
24 `gnomad.joint.v4.1.sites.<chrom>.vcf.bgz` files from this bucket.

## Reproducibility notes

* **Deterministic and version-pinned.** The dataset is fixed by the v4.1 release
  and by the per-file MD5 checksums embedded in the retrieval script
  `retrieve_gnomAD.sh`.
* **Checksums from the object store.** Each MD5 is the object's own MD5 as
  reported by Google Cloud Storage (`x-goog-hash: md5=…`, base64-decoded to hex),
  so integrity is verifiable against the upstream source of truth. GCS reports
  identical stored and served content lengths (no transcoding), so a downloaded
  file's MD5 matches the stored digest.
* **Open access, automated end-to-end.** The source needs no credentials, and
  beyond the download the retrieval script now performs the two preparation steps
  that turn the sites-only files into CRISPRme inputs (no manual conversion
  required):
  * **Convert** (phase B) — `crisprme.py gnomAD-converter --joint` turns the
    sites-only VCFs into CRISPRme-compatible, biallelic, PASS-filtered VCFs with
    one pseudo-sample column per ancestry group (reading `INFO/AF_joint`).
  * **MAF-filter** (phase C) — `bcftools view -i 'INFO/AF > 0.001'` on the
    converted files writes the `more_than_001/` subset (bgzipped + indexed), the
    variant of this dataset used in the paper's combined multi-source panel.
  The converter reads the `samplesIDs.gnomad.v41.txt` bundled next to the script,
  so the pipeline runs with no extra inputs.

---

## Known limitations

* **Aggregate, sites-only.** There is no individual-level data anywhere in this
  dataset; it cannot support per-sample analyses. CRISPRme uses ancestry groups
  as pseudo-samples, which enables population-level statistics only.
* **Multiallelic sites not split** in the raw files (see
  [Variant types](#variant-types)).
* **Very large.** ~817 GB for the 24 files; fetch subsets where possible and
  ensure ample disk.
* **Reduced sensitivity in hard regions.** Like any short-read aggregate call
  set, gnomAD is expected to have reduced sensitivity / coverage gaps in
  centromeric, segmental-duplication, and other low-complexity regions (not
  independently quantified here).

---

## References

* Chen S, Francioli LC, Goodrich JK, et al. *A genomic mutational constraint map
  using variation in 76,156 human genomes.* Nature. 2024;625(7993):92–100.
* gnomAD v4.1 release notes:
  [gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1)
* gnomAD downloads: [gnomad.broadinstitute.org/downloads](https://gnomad.broadinstitute.org/downloads)
