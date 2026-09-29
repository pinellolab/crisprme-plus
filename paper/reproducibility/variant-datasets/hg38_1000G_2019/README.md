# 1000 Genomes Project — GRCh38 2019 biallelic reanalysis (`hg38_1000G_2019`)

Phased, population-scale germline variants from the **1000 Genomes Project (1KGP) Phase 3** sample panel, called directly on the **GRCh38 / hg38** assembly (the March 2019 *biallelic SNV + INDEL* release).

> **Legacy dataset.** For CRISPRme+ this low-coverage (~7×) 2019 reanalysis is **superseded by the high-coverage (~30×) 2021 release** (`hg38_1000G_2021`), which is the 1000 Genomes source used for CRISPRme+ indexes. The 2019 set is retained here for paper reproducibility and is **not** part of the shipped/merged CRISPRme+ index.

---

## Summary

| Property | Value |
|---|---|
| Dataset | 1000 Genomes Project Phase 3 — GRCh38 reanalysis |
| Release identifier | `20190312_biallelic_SNV_and_INDEL` |
| Reference build | GRCh38 / hg38 (`GRCh38_full_analysis_set_plus_decoy_hla.fa`) |
| Samples | 2,548 |
| Populations | 26 populations across 5 super-populations (AFR, AMR, EAS, EUR, SAS) |
| Variant types | Biallelic SNVs and short INDELs |
| Genotypes | Phased (SHAPEIT2-integrated) |
| File format | bgzip-compressed VCF (`VCFv4.3`) |
| Files | `ALL.<chrom>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz` |
| Chromosomes | chr1–chr22, chrX (23 files; no chrY / chrM) |
| Contig naming | **bare, no `chr` prefix** (`22`, `X`, …) |
| Source | International Genome Sample Resource (IGSR) @ EBI |
| Access | Open / unrestricted |
| Integrity | Per-file MD5 checksums (see [Reproducibility notes](#reproducibility-notes)) |

---

## Background

The original 1000 Genomes Project Phase 3 call set was produced on the GRCh37 (hg19) assembly. To make the panel usable with modern GRCh38-based pipelines, the project reanalysed the Phase 3 samples **on GRCh38**, realigning the reads and calling variants *de novo* on the new assembly rather than lifting the old GRCh37 coordinates over. The genotypes were then integrated and phased with SHAPEIT2. The result is the `shapeit2_integrated_snvindels_v2a_27022019` call set released on 2019-03-12.

This directory tracks the **biallelic SNV + INDEL** subset of that release: one phased VCF per chromosome, restricted to biallelic sites.

---

## Data description

### Samples and populations

The phased VCFs contain **2,548 samples** (standard 1KGP `HG*` / `NA*` identifiers) drawn from **26 populations** grouped into five continental super-populations:

* **AFR** — African
* **AMR** — Admixed American
* **EAS** — East Asian
* **EUR** — European
* **SAS** — South Asian

The sample → population → super-population mapping is described under [Sample and population metadata](#sample-and-population-metadata) below.

### Variant types

* **Biallelic SNVs** and **short INDELs** only (one ALT allele per record).
* Produced with **multiple variant callers** whose call sets were integrated before final genotyping.
* Genotypes are **phased** (`GT` uses the `|` separator, integrated with SHAPEIT2), suitable for haplotype-aware analyses.

### Reference genome and coordinates

* Aligned and called against **GRCh38** (`GRCh38_full_analysis_set_plus_decoy_hla.fa`, as recorded in the VCF `##reference` header). Positions follow the 1-based VCF convention.
* **Contig names carry no `chr` prefix** (`##contig=<ID=22>`; records use `CHROM=22`, `X`, …). This differs from the UCSC hg38 reference, whose contig are `chr`-prefixed. Any step that combines these VCFs with the UCSC hg38 FASTA (for example, building a combined multi-source panel) must first rename the contigs to the `chr`-prefixed form; CRISPRme's own enrichment handles this matching internally.

### Allele-frequency (INFO) fields

Each record carries allele counts and frequencies in the `INFO` column, including the global `AC`, `AN`, `AF` and per-super-population frequencies (`AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF`). For example, the first record of chr22 is `22:10516173 A>G` with `AF=0.02`.

### Record counts

Per-chromosome variant counts (from the release's tabix indexes, via `bcftools index`):

| Chrom | Records | Chrom | Records | Chrom | Records |
|---|---|---|---|---|---|
| chr1 | 6,191,833 | chr9 | 3,384,360 | chr17 | 2,209,149 |
| chr2 | 6,790,551 | chr10 | 3,874,259 | chr18 | 2,189,529 |
| chr3 | 5,641,493 | chr11 | 3,881,791 | chr19 | 1,738,824 |
| chr4 | 5,477,810 | chr12 | 3,745,465 | chr20 | 1,817,492 |
| chr5 | 5,115,036 | chr13 | 2,760,845 | chr21 | 1,045,269 |
| chr6 | 4,863,337 | chr14 | 2,548,903 | chr22 | 1,059,079 |
| chr7 | 4,511,408 | chr15 | 2,301,453 | chrX | 106,963 |
| chr8 | 4,425,449 | chr16 | 2,548,920 | | |

In this biallelic release chrX carries far fewer records than the autosomes (106,963); interpret chrX-based statistics with that in mind.

### Sample and population metadata

The sample → population → super-population mapping is **not** stored in the VCFs; two upstream sources are relevant, depending on the need:

* **To run CRISPRme**, use CRISPRme's own sample-ID list, `test/data/samplesIDs/samplesIDs.1000G.txt` in the CRISPRme repository (installed as `hg38_1000G.samplesID.txt` and pulled automatically by the tool's setup). 
* **For the underlying population definitions** (the 26 populations and five super-populations, with their sample assignments) refer to the standard 1000 Genomes / IGSR sample metadata, available through the IGSR data portal ([internationalgenome.org](https://www.internationalgenome.org/data-portal/)).

### Provenance & population-aware search notes

- **Status:** **legacy** — superseded by the high-coverage 2021 release (`hg38_1000G_2021`) for CRISPRme+ indexes; retained for reproducibility.
- **Phasing:** phased (SHAPEIT2-integrated; `GT` uses the `|` separator).
- **indel+SNP co-occurrence:** **CONFIRMED-capable** (phased) — a co-occurring SNV and indel can be resolved onto the same observed sample haplotype and reported as CONFIRMED rather than dropped. (In practice CRISPRme+ uses the 2021 release for this.)
- **Conventions:** MAF filter **none** (rare and singleton alleles retained; sensitivity bounded only by the ~7× source coverage); contigs **UNPREFIXED** (`22`, `X`, …) — rename to `chr`-prefixed before combining with the UCSC hg38 FASTA.

---

## Source

The VCFs are distributed by the **International Genome Sample Resource (IGSR)**, hosted at EBI, under the Phase 3 GRCh38 reanalysis release:

```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
```

The same tree is reachable over `ftp://` on the same host.

## Reproducibility notes

* **Deterministic and version-pinned.** The dataset is fixed by its release identifier (`20190312_biallelic_SNV_and_INDEL`) and by the per-file MD5 checksums embedded in the retrieval script. A successful run yields byte-for-byte the same files on any machine.
* **No manual or credentialed steps.** The source is openly accessible; the retrieval is fully automated and needs no login or data-use agreement.

---

## Known limitations

* **Low-coverage source data.** The Phase 3 samples were sequenced at low coverage (~7×); sensitivity for rare variants is lower than in the more recent high-coverage (30×) 1000 Genomes releases (e.g. the 2020/2021 NYGC 3,202-sample call set, tracked separately as `hg38_1000G_2021`).
* **Biallelic subset only.** Multiallelic sites are represented as biallelic records in this release; analyses needing the full multiallelic representation should consult the complete project release.

---

## References

1. Zheng-Bradley X, Streeter I, Fairley S, Richardson D, Clarke L, Flicek P; 1000 Genomes Project Consortium. *Alignment of 1000 Genomes Project reads to reference assembly GRCh38.* GigaScience. 2017;6(7):1–8. doi:10.1093/gigascience/gix038. PMID: 28531267; PMCID: PMC5522380.
2. Lowy-Gallego E, Fairley S, Zheng-Bradley X, Ruffier M, Clarke L, Flicek P; 1000 Genomes Project Consortium. *Variant calling on the GRCh38 assembly with the data from phase three of the 1000 Genomes Project.* Wellcome Open Res. 2019;4:50. doi:10.12688/wellcomeopenres.15126.2. PMID: 32175479; PMCID: PMC7059836.
