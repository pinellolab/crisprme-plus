# 1000 Genomes Project — GRCh38 2019 biallelic reanalysis (`hg38_1000G_2019`)

Phased, population-scale germline variants from the **1000 Genomes Project
(1KGP) Phase 3** sample panel, called directly on the **GRCh38 / hg38**
assembly (the March 2019 *biallelic SNV + INDEL* release). This is the default
1000 Genomes resource used by CRISPRme+ for population-aware off-target
nomination.

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
| Integrity | Per-file MD5 checksums (see [retrieval](#how-the-data-are-retrieved)) |

---

## Background

The original 1000 Genomes Project Phase 3 call set was produced on the GRCh37
(hg19) assembly. To make the panel usable with modern GRCh38-based pipelines,
the project reanalysed the Phase 3 samples **on GRCh38** — realigning the reads
and calling variants *de novo* on the new assembly rather than lifting the old
GRCh37 coordinates over. The genotypes were then integrated and phased with
SHAPEIT2. The result is the `shapeit2_integrated_snvindels_v2a_27022019`
call set released on 2019-03-12.

This directory tracks the **biallelic SNV + INDEL** subset of that release: one
phased VCF per chromosome, restricted to biallelic sites.

---

## Data description

### Samples and populations

The phased VCFs contain **2,548 samples** (standard 1KGP `HG*` / `NA*`
identifiers) drawn from **26 populations** grouped into five continental
super-populations:

* **AFR** — African
* **AMR** — Admixed American
* **EAS** — East Asian
* **EUR** — European
* **SAS** — South Asian

The sample → population → super-population mapping is **not** part of these VCFs;
CRISPRme reads it from a separate `samplesID` metadata file (obtained
independently of this retrieval script — see [Related files](#related-files)).

### Variant types and phasing

* **Biallelic SNVs** and **short INDELs** only (one ALT allele per record).
* Genotypes are **phased** (`GT` uses the `|` separator), suitable for
  haplotype-aware analyses.

### Reference genome and coordinates

* Aligned and called against **GRCh38**
  (`GRCh38_full_analysis_set_plus_decoy_hla.fa`, as recorded in the VCF
  `##reference` header). Positions follow the 1-based VCF convention.
* **Contig names carry no `chr` prefix** (`##contig=<ID=22>`; records use
  `CHROM=22`, `X`, …). This differs from the UCSC hg38 reference, whose contigs
  are `chr`-prefixed. Any step that combines these VCFs with the UCSC hg38 FASTA
  (for example, building a combined multi-source panel) must first rename the
  contigs to the `chr`-prefixed form; CRISPRme's own enrichment handles this
  matching internally.

### Allele-frequency (INFO) fields

Each record carries allele counts and frequencies in the `INFO` column,
including the global `AC`, `AN`, `AF` and per-super-population frequencies
(`AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF`). For example, the first
record of chr22 is `22:10516173 A>G` with `AF=0.02`.

---

## Source

The VCFs are distributed by the **International Genome Sample Resource (IGSR)**,
hosted at EBI, under the Phase 3 GRCh38 reanalysis release:

```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
```

The same tree is reachable over `ftp://` on the same host. The data are openly
accessible and require no credentials. This is also the exact release that
CRISPRme's built-in `crisprme.py setup` downloads for its default 1000 Genomes
data, so the resource reproduced here matches the tool's shipped configuration.

## Reproducibility notes

* **Deterministic and version-pinned.** The dataset is fixed by its release
  identifier (`20190312_biallelic_SNV_and_INDEL`) and by the per-file MD5
  checksums embedded in the retrieval script. A successful run yields byte-for-
  byte the same files on any machine.
* **No manual or credentialed steps.** The source is openly accessible; the
  retrieval is fully automated and needs no login or data-use agreement.
* **Provenance.** The URLs and checksums mirror CRISPRme's own
  `setup_legacy_database.py` / `utils.py`, so this reproducibility bundle stays
  consistent with the tool's default data.

---

## Known limitations

* **Low-coverage source data.** The Phase 3 samples were sequenced at low
  coverage (~7×); sensitivity for rare variants is lower than in the more recent
  high-coverage (30×) 1000 Genomes releases (e.g. the 2020/2021 NYGC 3,202-sample
  call set, tracked separately as `hg38_1000G_2021`).
* **Biallelic subset only.** Multiallelic sites are represented as biallelic
  records in this release; analyses needing the full multiallelic representation
  should consult the complete project release.

---

## References

1. Zheng-Bradley X, Streeter I, Fairley S, Richardson D, Clarke L, Flicek P;
   1000 Genomes Project Consortium. *Alignment of 1000 Genomes Project reads to
   reference assembly GRCh38.* GigaScience. 2017;6(7):1–8.
   doi:10.1093/gigascience/gix038. PMID: 28531267; PMCID: PMC5522380.
2. Lowy-Gallego E, Fairley S, Zheng-Bradley X, Ruffier M, Clarke L, Flicek P;
   1000 Genomes Project Consortium. *Variant calling on the GRCh38 assembly with
   the data from phase three of the 1000 Genomes Project.* Wellcome Open Res.
   2019;4:50. doi:10.12688/wellcomeopenres.15126.2. PMID: 32175479;
   PMCID: PMC7059836.
