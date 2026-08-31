# 1000 Genomes Project — GRCh38 high-coverage re-sequencing, 2021 (`hg38_1000G_2021`)

High-coverage (~30×) whole-genome re-sequencing of the **1000 Genomes Project (1KGP)** cohort, called natively on **GRCh38 / hg38**. This 2021 release provides deep Illumina WGS for **3,202 individuals** (2,504 unrelated plus 698 related forming 602 trios), with fully phased genotypes, a substantially more sensitive, GRCh38-native successor to the low-coverage 2019 reanalysis (tracked separately as `hg38_1000G_2019`). It is suitable for fine-mapping, imputation, and haplotype-based analyses.

---

## Summary

| Property | Value |
|---|---|
| Dataset | 1000 Genomes Project — high-coverage NYGC re-sequencing |
| Release identifier | `20201028_3202_phased` |
| Reference build | GRCh38 / hg38 |
| Coverage | ~30× WGS |
| Samples | 3,202 (2,504 unrelated + 698 related; 602 trios) |
| Populations | 26 populations across 5 super-populations (AFR, AMR, EAS, EUR, SAS) |
| Variant types | SNVs and short INDELs (GATK HaplotypeCaller, joint-genotyped, site-filtered) |
| Genotypes | Phased (autosomes SHAPEIT2-duohmm; chrX Eagle2 v2) |
| File format | bgzip-compressed VCF (`VCFv4.1`) |
| Files | autosomes `CCDG_14151_B01_GRM_WGS_2020-08-05_<chrom>.filtered.shapeit2-duohmm-phased.vcf.gz`; chrX `…chrX.filtered.eagle2-phased.v2.vcf.gz` |
| Chromosomes | chr1–chr22 (default); chrX available but excluded by default |
| Contig naming | `chr`-prefixed (`chr1`, …, `chrX`) |
| Source | International Genome Sample Resource (IGSR) @ EBI |
| Access | Open / unrestricted |
| Integrity | Per-file MD5 checksums (official manifest; see [Reproducibility notes](#reproducibility-notes)) |

---

## Background

The original 1000 Genomes Project Phase 3 was sequenced at low coverage (~7.4×) on the GRCh37 assembly. While foundational for human population genetics, the low depth limited rare-variant sensitivity and the GRCh37 coordinates restricted compatibility with modern pipelines.

The **2021 high-coverage re-sequencing** was carried out by the **New York Genome Center (NYGC)** and distributed through IGSR, motivated by the need to:

* dramatically improve rare-variant detection through high-depth (~30×) sequencing;
* provide a GRCh38-native call set without lift-over artifacts;
* leverage family structure (trios) for improved phasing via SHAPEIT2-duohmm;
* expand the cohort with 698 related individuals, enabling trio-based quality control.

---

## Data description

### Samples and populations

The dataset includes **3,202 individuals** across **26 populations** grouped into
five continental super-populations:

* **AFR** — African
* **AMR** — Admixed American
* **EAS** — East Asian
* **EUR** — European
* **SAS** — South Asian

The 3,202 total comprises 2,504 unrelated individuals (the classic 1KGP Phase 3 set) plus 698 additional related individuals forming 602 parent-offspring trios used to aid phasing. The related individuals should be accounted for in population-structure analyses (see [Known limitations](#known-limitations)). The sample → population mapping is described under [Sample and population metadata](#sample-and-population-metadata) below.

### Variant types

* **SNVs** and **short INDELs**.
* Called *de novo* at high coverage with **GATK HaplotypeCaller** (GVCF mode), followed by joint genotyping across all samples and site-level filtering (the `filtered` call set).
* Genotypes are **phased**: autosomes with **SHAPEIT2-duohmm** (trio-aware, for longer, more accurate haplotype blocks); chrX with **Eagle2 (v2)**.

### Reference genome and coordinates

* Aligned to **GRCh38 / hg38** with **BWA-MEM** (NYGC functional-equivalence protocol) and called natively on GRCh38. Positions follow the 1-based VCF convention.
* **Contig names are `chr`-prefixed** (`##contig=<ID=chr21>`; records use `CHROM=chr21`), matching the UCSC hg38 reference directly.
* Chromosomes covered: **chr1–chr22** and **chrX**; no chrY or chrM files are included in this release.

### Allele-frequency (INFO) fields

Each record carries allele counts and frequencies in the `INFO` column: the global `AC`, `AN`, `AF`, per-super-population frequencies (`AF_AFR`, `AF_AMR`, `AF_EAS`, `AF_EUR`, `AF_SAS`), and per-population allele-count breakdowns (`AC_Het` / `AC_Hom`, globally and per population). For example, the first record of chr21 is `chr21:5030578 C>T` with `FILTER=PASS`.

### Record counts

Per-chromosome variant counts (from the release's tabix indexes, via `bcftools index`):

| Chrom | Records | Chrom | Records | Chrom | Records |
|---|---|---|---|---|---|
| chr1 | 5,769,087 | chr9 | 3,169,328 | chr17 | 2,075,523 |
| chr2 | 6,095,976 | chr10 | 3,499,286 | chr18 | 1,965,907 |
| chr3 | 4,986,824 | chr11 | 3,425,446 | chr19 | 1,672,929 |
| chr4 | 4,878,537 | chr12 | 3,335,036 | chr20 | 1,647,102 |
| chr5 | 4,539,890 | chr13 | 2,512,948 | chr21 | 1,004,437 |
| chr6 | 4,317,093 | chr14 | 2,294,933 | chr22 | 1,070,401 |
| chr7 | 4,140,924 | chr15 | 2,111,611 | chrX | 2,858,184 |
| chr8 | 3,888,893 | chr16 | 2,366,114 | | |

### Sample and population metadata

The sample → population → super-population mapping is **not** stored in the VCFs. For this release the authoritative source is the official 3,202-sample pedigree / population file `20130606_g1k_3202_samples_ped_population.txt`, in the parent `1000G_2504_high_coverage/` directory at IGSR/EBI ([link](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt)). It is space-delimited, one row per sample, with columns `FamilyID`, `SampleID`, `FatherID`, `MotherID`, `Sex` (1 = male, 2 = female), `Population`, `Superpopulation`. A CRISPRme `samplesID` file (`#SAMPLE_ID  POPULATION_ID  SUPERPOPULATION_ID  SEX`) can be derived from it as a separate preparation step.

### Provenance & population-aware search notes

- **Phasing:** phased (autosomes SHAPEIT2-duohmm, trio-aware; chrX Eagle2 v2).
- **indel+SNP co-occurrence:** **CONFIRMED-capable.** Because SNVs and indels are jointly phased onto the same haplotypes, an off-target that requires an indel together with one or more nearby SNVs can be validated as co-occurring on a single observed haplotype and reported as **CONFIRMED** (not merely PUTATIVE).
- **Conventions:** MAF filter **none** (site-level QC only — rare and singleton alleles are retained, no allele-frequency cutoff); contigs `chr`-prefixed.
- **Caveats:**
  - The 698 related samples **inflate the `AN` denominator** relative to the 2,504-unrelated panel; account for relatedness when reporting or interpreting allele frequencies (use the unrelated subset if unbiased population frequencies are needed).
  - **chrX is excluded by default:** it is Eagle2-phased with haploid males encoded as `0/1` (not diploid-safe), so the default run is autosomes (chr1–chr22) only. Re-enable chrX only if your pipeline handles the haploid-male encoding.
  - No chrY or chrM in this release (see [Known limitations](#known-limitations)).

---

## Source

The VCFs are distributed by the **International Genome Sample Resource (IGSR)**, hosted at EBI, in the high-coverage re-sequencing `working/` tree:

```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
```

The data are openly accessible and require no credentials. The per-chromosome files keep their upstream NYGC names (prefixed `CCDG_14151_B01_GRM_WGS_2020-08-05_`), and the retrieval script downloads them under those original names so that checksums match the release's own manifest (`phased-manifest_July2021.tsv`).

## Reproducibility notes

* **Deterministic and version-pinned.** The dataset is fixed by its release identifier (`20201028_3202_phased`) and by the per-file MD5 checksums embedded in the retrieval script `retrieve_1000G_2021.sh`. A successful run yields byte-for-byte the same files on any machine.
* **Checksums from the official manifest.** The MD5s are taken verbatim from the release's own `phased-manifest_July2021.tsv`, so integrity is verifiable against the upstream source of truth.
* **No manual or credentialed steps for download.** The source is openly accessible; the download is fully automated. Two follow-up preparation steps are required before CRISPRme use (see [Known limitations](#known-limitations)): renaming the files to a CRISPRme-compatible pattern, and building the `samplesID` file.

---

## Known limitations

* **No chrY or chrM.** This release does not include Y-chromosome or mitochondrial variant files.
* **Working-directory snapshot.** Files were retrieved from the IGSR `working/` directory (pre-publication snapshot dated 2020-10-28 extended to include all SNVs and INDELs on 2022-11-18); the canonical published release also includes structural variants.

---

## References

1. Byrska-Bishop M, Evani US, Zhao X, et al. *High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios.* Cell. 2022;185(18):3426–3440.e19. doi:10.1016/j.cell.2022.08.004. PMID: 36055201; PMCID: PMC9560234.
2. Lowy-Gallego E, Fairley S, Zheng-Bradley X, Ruffier M, Clarke L, Flicek P; 1000 Genomes Project Consortium. *Variant calling on the GRCh38 assembly with the data from phase three of the 1000 Genomes Project.* Wellcome Open Res. 2019;4:50. doi:10.12688/wellcomeopenres.15126.2. PMID: 32175479; PMCID: PMC7059836.
