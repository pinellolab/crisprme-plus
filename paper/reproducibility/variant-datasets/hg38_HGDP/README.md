# Human Genome Diversity Project — GRCh38 high-coverage WGS (`hg38_HGDP`)

High-coverage whole-genome sequencing (WGS) variant calls from the **Human
Genome Diversity Project (HGDP)**, aligned and called against **GRCh38 / hg38**
(the 2019-05-16 release). HGDP documents human genetic diversity across
geographically and ethnolinguistically diverse populations and — together with
the 1000 Genomes Project — is one of the two population resources CRISPRme+ uses
for population-aware off-target nomination. This is CRISPRme's own canonical HGDP
source: the same file set retrieved by the tool's built-in setup, and the HGDP
data used in the original CRISPRme publication (Cancellieri et al., 2023).

---

## Summary

| Property | Value |
|---|---|
| Dataset | Human Genome Diversity Project (HGDP) — high-coverage WGS |
| Release identifier | `hgdp_wgs.20190516` ("full" call set) |
| Reference build | GRCh38 / hg38 (`GRCh38_full_analysis_set_plus_decoy_hla.fa`) |
| Samples | 929 (autosomes); 927 (chrX); 610 (chrY) |
| Populations | 54 populations across 7 geographic regions |
| Variant types | SNVs and short INDELs (GATK HaplotypeCaller + VQSR); multiallelic, not pre-split |
| File format | bgzip-compressed VCF (`VCFv4.2`) |
| Files | `hgdp_wgs.20190516.full.<chrom>.vcf.gz` (upstream release names) |
| Chromosomes | chr1–chr22, chrX (default); chrY also available |
| Contig naming | `chr`-prefixed (`chr1`, …, `chrX`, `chrY`) |
| Source | Wellcome Sanger Institute |
| Access | Open / unrestricted |
| Integrity | Per-file MD5 checksums (see [Reproducibility notes](#reproducibility-notes)) |

---

## Background

The HGDP was established to capture the full spectrum of human genetic diversity
across diverse ethnolinguistic groups. The high-coverage WGS reanalysis released
on 2019-05-16 (Bergström et al., 2020) provides deep-sequencing data spanning 54
populations, complementing cohorts such as the 1000 Genomes Project by including
deeply divergent lineages from Africa, Oceania, and the Americas. Its high
coverage — in contrast to the low-coverage 1000 Genomes Phase 3 data — yields
high genotype accuracy and strong sensitivity for rare variants.

---

## Data description

### Samples and populations

* **929 samples** on the autosomes (chr1–22), consistent across autosomes
  (verifiable with `bcftools query -l`). Public HGDP metadata describes these as
  spanning **54 populations** across **seven broad geographic regions** (Africa,
  Europe, the Middle East, Central/South Asia, East Asia, Oceania, the Americas).
* **927 samples** on chrX — two fewer than the autosomes, reflecting two samples
  with documented sex-chromosome aneuploidy (XO, i.e. Turner syndrome).
* **610 samples** on chrY — present in the male subset only.

HGDP includes several "bottlenecked" or isolated populations; expect higher
homozygosity and distinct linkage-disequilibrium patterns in these groups (see
[Known limitations](#known-limitations)).

### Variant types

* **SNVs and short INDELs**, called with **GATK (HaplotypeCaller)** and filtered
  with **VQSR** (the `VQSLOD` score and the usual GATK annotations —
  `BaseQRankSum`, `FS`, `InbreedingCoeff`, `ExcessHet`, … — are present in the
  `INFO` column).
* **Multiallelic sites are not pre-split**: a record may carry a comma-separated
  `ALT`. Downstream tools that assume strictly biallelic records must first run
  `bcftools norm -m -any`.

### Reference genome and coordinates

* Aligned and called against **GRCh38**
  (`GRCh38_full_analysis_set_plus_decoy_hla.fa`, as recorded in the VCF
  `##reference` header). Positions follow the 1-based VCF convention.
* **Contig names are `chr`-prefixed** (`##contig=<ID=chr1>`, …; records use
  `CHROM=chr22`). They therefore match the UCSC hg38 reference directly — in
  contrast to the 1000 Genomes 2019 VCFs, whose contigs are unprefixed and must
  be renamed before they can be combined with the same reference.

### Allele-frequency (INFO) fields

Each record carries allele counts and frequency in the `INFO` column (`AC`, `AN`,
`AF`, `DP`), alongside the GATK/VQSR quality annotations noted above. For example,
the first record of chr22 is `chr22:10510212 A>T` with `AF=0.00295858`,
`FILTER=PASS`.

### Record counts

Per-chromosome variant counts (verified with `bcftools index -n`):

| Chrom | Records | Chrom | Records | Chrom | Records |
|---|---|---|---|---|---|
| chr1 | 6,330,165 | chr9 | 3,424,990 | chr17 | 2,221,758 |
| chr2 | 6,387,900 | chr10 | 3,704,986 | chr18 | 2,127,418 |
| chr3 | 5,264,762 | chr11 | 3,641,107 | chr19 | 1,759,281 |
| chr4 | 5,152,732 | chr12 | 3,554,651 | chr20 | 1,820,421 |
| chr5 | 4,782,993 | chr13 | 2,813,331 | chr21 | 1,086,522 |
| chr6 | 4,495,706 | chr14 | 2,400,086 | chr22 | 1,185,008 |
| chr7 | 4,381,365 | chr15 | 2,223,509 | chrX | 2,654,850 |
| chr8 | 4,068,606 | chr16 | 2,483,073 | chrY | 132,457 |

### Sample and population metadata

The sample → population mapping is **not** stored in this directory; two distinct
upstream sources are relevant, depending on the need:

* **To run CRISPRme**, use CRISPRme's own sample-ID list,
  `test/data/samplesIDs/samplesIDs.HGDP.txt` in the CRISPRme repository (pulled
  automatically by the tool's setup). This is the file the CRISPRme pipeline
  actually consumes to match VCF samples. It is served from the `main` branch of
  the pinellolab repository — a moving target, not a pinned commit/tag, matching
  what CRISPRme's own setup does.
* **For the real population / sex / QC metadata** (population name, latitude and
  longitude, continental region, sex, sequencing coverage and QC), the official
  Sanger release ships its own metadata file alongside the VCFs, at
  `https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/hgdp_wgs.20190516.metadata.txt`.
  It is tab-separated, one row per sample, with columns: `sample`, `library`,
  `sample_accession`, `source`, `library_type`, `population`, `latitude`,
  `longitude`, `region`, `sex`, `coverage`, `freemix`, `capmq`,
  `insert_size_average`, `array_non_reference_discordance`,
  `library_alias_ENA`. Note this content differs from CRISPRme's ID list.

---

## Source

The VCFs are distributed by the **Wellcome Sanger Institute**, from the official
HGDP WGS release dated 2019-05-16:

```
https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/
```

The same tree is reachable over `ftp://` on the same host. The data are openly
accessible and require no credentials, and the per-chromosome files keep their
upstream release names (`hgdp_wgs.20190516.full.<chrom>.vcf.gz`, not renamed
locally). This is also the exact release that CRISPRme's built-in setup
downloads for its default HGDP data (`setup_legacy_database.py`,
`VCF_HGDP_SERVER = "ngs.sanger.ac.uk"`), so the resource reproduced here matches
the tool's shipped configuration.

## Reproducibility notes

* **Deterministic and version-pinned.** The dataset is fixed by its release
  identifier (`hgdp_wgs.20190516`) and by the per-file MD5 checksums embedded in
  the retrieval script `retrieve_HGDP.sh`. A successful run yields byte-for-byte
  the same files on any machine.
* **No manual or credentialed steps.** The source is openly accessible; the
  retrieval is fully automated and needs no login or data-use agreement.
* **Provenance.** The URLs and checksums mirror CRISPRme's own
  `setup_legacy_database.py` / `utils.py` (`MD5HGDP`), so this reproducibility
  bundle stays consistent with the tool's default data.

---

## Known limitations

* **Multiallelic sites are not pre-split.** Records may carry a comma-separated
  `ALT`; normalise with `bcftools norm -m -any` before any strictly biallelic
  analysis (see [Variant types](#variant-types)).
* **Chromosome sample-count differences.** chrX (927) and chrY (610) contain
  fewer samples than the autosomes (929), reflecting sex-chromosome aneuploidy
  and the male-only chrY subset (see [Samples and populations](#samples-and-populations)).
* **Isolated / bottlenecked populations.** Several HGDP groups are strongly
  bottlenecked; higher homozygosity and distinct LD patterns are expected and
  should be considered when interpreting population-level statistics.

---

## References

* Bergström A, et al. *Insights into human genetic variation and population
  history from 929 high-coverage genome sequences.* Science. 2020;367(6484):eaay5012.
* Cancellieri S, Zeng J, Lin LY, et al. *Human genetic diversity alters
  off-target outcomes of therapeutic gene editing.* Nat Genet. 2023;55(1):34–43.
  — the original CRISPRme publication; uses 1000 Genomes + HGDP (this dataset)
  for population-aware off-target nomination.
* International Genome Sample Resource (IGSR):
  [HGDP Data Collection Portal](https://www.internationalgenome.org/data-portal/data-collection/HGDP).
