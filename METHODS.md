# CRISPRme+ — Methods

This document describes the methods introduced in **CRISPRme+** (the 2.3/2.4
line), intended as a self-contained technical reference and as source material
for the Methods section of the manuscript. It focuses on what is **new or
changed** relative to the original CRISPRme (Cancellieri, Zeng, Lin et al.,
*Nature Genetics* 2023); the core enumeration of candidate off-targets by
CRISPRitz and the CFD / CRISTA scoring are unchanged unless stated otherwise.

Sections:
1. [Variant-aware, dictionary-less data model](#1-variant-aware-dictionary-less-data-model)
2. [Homogenization and merging of population VCFs](#2-homogenization-and-merging-of-population-vcfs)
3. [Allele-frequency estimation](#3-allele-frequency-estimation)
4. [Haplotype scanning: observed-haplotype enumeration for phased, unphased, and mixed data](#4-haplotype-scanning-observed-haplotype-enumeration)
5. [Search-space control for high-variant-density regions](#5-search-space-control-for-high-variant-density-regions)
6. [Functional annotation of off-targets](#6-functional-annotation-of-off-targets)
7. [Shareable off-target assessment report](#7-shareable-off-target-assessment-report)

Throughout, "protospacer window" means the genomic interval spanned by a
candidate off-target's protospacer plus PAM (and any bulges), i.e. the interval
in which an overlapping genetic variant can create, destroy, or modify an
off-target site.

---

## 1. Variant-aware, dictionary-less data model

### Motivation
The original CRISPRme makes off-target search variant-aware by *enriching* the
reference genome with population variants and, for reporting, by storing — for
every variant position — the genotype of **every sample** in a per-chromosome
JSON "dictionary" (`my_dict_<chrom>.json`). For the combined 1000 Genomes +
HGDP panel this per-sample store is ≈**152 GB**. That per-sample detail is only
required for *sample-attribution* features (which individuals / populations
carry a given off-target). Finding the variant off-targets themselves and
computing their allele frequency and rsID does **not** require per-sample
genotypes — only the variant alleles and their frequency at each candidate
position.

### Two-tier compact representation
CRISPRme+ replaces the monolithic per-sample dictionary with two tiers:

- **Tier-0 registry (`tier0_registry.py`).** A compact, memory-mapped binary
  index. Each record is one `(position, alternative allele)` pair and stores the
  aggregate counts needed for off-target detection and frequency reporting:
  allele count (AC), allele number (AN), carrier / homozygote / called-individual
  counts, and the rsID — resolved **per group**, where a group is a database
  (e.g. `1000G`), a database × super-population cell (e.g. `1000G::EUR`), and a
  deduplicated **global** aggregate. Records are a fixed-width, sorted array, so
  a lookup is an `O(log n)` binary search directly on the `mmap` — the file is
  never parsed into Python objects. Field widths are chosen per file to minimize
  size (positions as `u32`, counts as the smallest width that cannot overflow).

- **Tier-1 genotype store (`tier1_genotypes.py`).** A compact per-sample carrier
  representation consulted **only** when sample-level attribution is requested
  (the "Samples" column, personal risk cards, per-sample summaries). It
  reconstructs, for a `(position, alt)`, the list of carrier samples and their
  phased or unphased genotypes.

Population summaries (per-database, per-super-population, and global frequency
distributions) are **first-class** and computed at build time, so they are
available without touching Tier-1; per-individual queries remain lazy.

### Compact on-disk encoding of the Tier-0 registry
Physically, the registry is three contiguous sections: a fixed-width **record
array** (16 bytes per record — position, reference base, alternate base, the
record's group count, and offsets into the two following sections), a **group
blob** (the per-group AC / AN / carrier / homozygote / called-allele counts), and
a deduplicated **string pool** (rsIDs). The group blob dominates the file
(~75–80 %), and within it the allele-number and called-allele columns are
near-constant across the millions of records on a chromosome — the genotyped
panel is essentially the same size everywhere — so the raw layout is highly
redundant.

CRISPRme+ therefore stores the registry in a **block-compressed** format: each of
the three sections is partitioned into fixed-size blocks (4,096 records) that are
individually `zlib`-compressed, with a small sparse index of block offsets carried
in the JSON manifest. A lookup bisects that sparse index to the one covering
block, decompresses it (kept in a tiny LRU so clustered/adjacent lookups stay
warm), and reads it exactly as before — so the public reader API and the
`O(log n)` random access by genomic position are unchanged, while the
near-constant count columns compress away. On the shipped 1000 Genomes + HGDP
registry this is ≈3.5× on disk (≈7.5 GB → ≈2.1 GB after extraction). The reader is
fully backward-compatible (it still reads an uncompressed registry byte-for-byte),
and an existing registry is re-encoded **losslessly** by a verbatim block
re-chunking (`transcode_registry`) with no VCF re-parse — the records, counts,
rsIDs and lookups are identical, only the container changes.

### Out-of-the-box variant search
Because Tier-0/Tier-1 are small, they are shipped **with the pre-built index**
(compressed, read on the fly). A user who downloads a variant index can
therefore run variant-aware searches — with correct allele frequencies and
rsIDs — without ever materializing the 152 GB of per-sample dictionaries. The
enriched reference genome and the raw VCFs are needed only at *build* time.

---

## 2. Homogenization and merging of population VCFs

CRISPRme+ supports combining several population resources (e.g. 1000 Genomes,
HGDP, gnomAD, HPRC, All-of-Us) into a single variant panel. Correct frequencies
across heterogeneous sources require careful homogenization.

### Homogenization
Each source VCF is normalized so that variant records are directly comparable:

- **Multiallelic decomposition.** Multiallelic sites are split into biallelic
  records (`bcftools norm -m-`) and left-aligned against the reference. The
  carried alternative allele is recoded to the canonical `1` token, so a record
  is unambiguously per-alt. This is the convention the Tier-0 aggregator assumes
  (`alt_index="1"`): a genotype token equal to `1` is *this record's* alt, `0`
  or any foreign alt index (`2`, …) is treated as a called non-carrier allele,
  and `.` is missing. This prevents a `1|2` genotype from being miscounted as
  two copies of alt `1`.
- **Consistent coordinates and contigs** (reference build, chromosome naming)
  across datasets.

### Merging and per-dataset provenance
Datasets are merged (`bcftools merge`) into a combined panel, and allele
frequencies are recomputed on the merged multiallelic records (this depends on
the multiallelic-AF fix in CRISPRitz PR #36). Crucially, CRISPRme+ **preserves
dataset provenance**: allele frequencies are reported per native dataset label
(1000G vs HGDP vs gnomAD are never conflated), and a combined global frequency
is reported over the union panel. This matters because the same variant can have
very different frequencies across ancestries, and a merge must not silently
average them away.

### Panel definition (the genotyped sample set)
The **allele number (AN)**, i.e. the denominator of every allele frequency, is
the number of *genotyped* alleles in the panel. A subtle but consequential point
is that a dataset's sample roster (`samplesID`) may list more samples than are
actually genotyped in its VCF (for the phased 1000 Genomes VCF used here, the
roster lists 3,500 samples but only 2,548 are genotyped). Counting the un-
genotyped "phantom" samples as called reference inflates AN and biases every
allele frequency low. CRISPRme+ therefore defines the panel from the samples
**actually present in the VCF** (VCF-filtered `samplesID`), giving the correct
AN (here 2×(2,548 + 929) = 6,954 for combined 1000G+HGDP autosomes).

---

## 3. Allele-frequency estimation

For a candidate variant off-target, the reported minor/alternate allele
frequency is

  AF = AC / AN

computed **per group** (global, per-database, per-super-population) from the
Tier-0 registry. The aggregator is ploidy-aware: on autosomes every sample
contributes ploidy 2; on chrX-nonPAR a male contributes ploidy 1 (and is never
double-counted as homozygous), and on chrY females contribute ploidy 0 (they add
no phantom alleles to the denominator). Missing genotypes are treated as
reference for the panel denominator ("missing-as-ref"), so AN is the full
genotyped-panel called-allele count rather than the alleles among carriers only.
Homozygous carriers contribute 2 to AC; heterozygous carriers contribute 1 —
so AF is a true allele frequency, not a carrier-individual frequency.

For a multi-variant haplotype off-target, the reported frequency is bounded by
the **rarest** contributing allele (the haplotype can be no more frequent than
its least-frequent variant), and the exact carriers of the haplotype (Section 4)
give the tightest available estimate.

---

## 4. Haplotype scanning: observed-haplotype enumeration

### The problem
When more than one genetic variant falls inside a single protospacer window, the
off-target that a real genome presents depends on **which combination of those
variants co-occurs on the same DNA molecule** (haplotype). With *k* variant
positions in the window there are up to 2^*k* possible allele combinations, but
the vast majority never occur in any real individual. Two failure modes must be
avoided:

- **False positives (phantom off-targets).** Enumerating a worst-case
  combination attributed to the *union* of all carriers describes a haplotype
  that **no single individual carries**. In dense, low-complexity regions this
  can stack many variants into a fictional high-scoring off-target. (Empirically
  we observed windows where the union stacked 10–28 variants while the maximum
  *real* cis haplotype carried by any individual was 2.)
- **False negatives (dropped haplotypes).** Conversely, collapsing or capping
  combinations can drop a genuine multi-variant haplotype that a real individual
  carries, under-reporting a real off-target.

### Observed-haplotype enumeration
CRISPRme+ enumerates, for each protospacer window, exactly the variant-sets
(haplotypes) that occur in **at least one real individual**, using the genotypes
in the Tier-1 store (`observed_haplotypes.py`). The number of enumerated
haplotypes is bounded by ≈2×(number of carriers) rather than 2^*k*, so it is
efficient, and each enumerated haplotype carries its **exact set of carriers**
(which in turn yields its exact frequency, Section 3). Cross-individual chimeras
— combinations assembled from variants carried by *different* people — are
**excluded** by construction.

### Behavior by phasing status
The confidence attached to an enumerated haplotype depends on the phasing of the
input genotypes:

- **Phased data** (alleles separated by `|`). The cis/trans configuration is
  known, so a per-individual, per-chromosome variant-set is a **confirmed** cis
  haplotype: that individual carries exactly those variants together on one
  chromosome.
- **Unphased data** (alleles separated by `/`). The cis/trans configuration is
  unknown, so **every** cis/trans arrangement of the individual's carried variants
  is possible. CRISPRme+ therefore reports **every non-empty subset** of that
  individual's variant-set as a **putative** haplotype (each flagged as
  unconfirmed, none dropped). Enumerating the subsets — not only the maximal union
  — is essential because a variant can *break* an off-target as well as create one:
  e.g. a variant that disrupts the PAM must be droppable, so that a sub-combination
  which keeps the reference allele at that position (and is a genuine off-target
  under one possible phasing) is not hidden by the all-variants union. The subsets
  are the individual's **own** variants only (never cross-individual chimeras); the
  scoring/PAM/mismatch-budget gates prune subsets that are not in-budget PAM-valid
  targets, and identical subsets carried by different individuals are deduplicated.
  This is the sensitivity-first choice — we never omit a haplotype a real
  individual could plausibly carry — bounded by `CRISPRME_IUPAC_CAP`: an individual
  carrying more variants than the cap in a single window falls back to the union
  (the combinatorial blow-up is confined to that individual, and the dense window
  is surfaced in the high-variant-density BED).
- **Mixed data** (e.g. a merge of phased and unphased sources, or block-phased
  VCFs from WhatsHap/GATK/HapCUT2). Genotypes are handled conservatively:
  same-individual variant-sets are enumerated; combinations that would require
  assuming cis across an unknown or different **phase set (PS)** are reported as
  putative rather than confirmed. Phase-set awareness prevents wrongly fusing
  `1|0` at position A and `1|0` at position B when A and B lie in different
  phasing blocks and are therefore not known to be in cis. When phase-set
  information is unavailable, the whole-chromosome statistical phasing of
  resources like 1000 Genomes and HGDP is assumed (documented), and the
  conservative putative labeling absorbs the residual uncertainty.

Confirmed and putative haplotypes are reported distinctly, so a reviewer can
weight them appropriately.

### Locus completeness
For every candidate window, CRISPRme+ additionally emits the **reference**
off-target (the site as it appears in the reference genome, independent of any
variant), so that a locus is never dropped merely because it lacks a productive
variant haplotype. Variant off-targets are reported alongside the reference
off-target rather than replacing it.

### Relationship to the original method
The original CRISPRme performs the multi-variant combination step with a greedy
cap and, on the dictionary-less path, could additionally mis-select the phased
vs. unphased branch. CRISPRme+ replaces the greedy union with per-individual
observed-haplotype enumeration driven directly by the genotype tier, which
removes the phantom off-targets, restores dropped real haplotypes, and yields
exact carriers (hence exact frequencies) for every reported combination.

---

## 5. Search-space control for high-variant-density regions

Variant-aware search can explode combinatorially where an IUPAC-dense,
low-complexity region coincides with a permissive search (many mismatches/bulges,
minimal PAM constraint, unphased genotypes). CRISPRme+ bounds this with three
complementary controls:

- **`max_total_edits`.** The total number of edits (mismatches + bulges) of the
  **reconstructed** reference/alternate alignment is enforced against the
  user-requested budget, so a reported off-target never silently exceeds the
  stated edit distance.
- **High-variant-density cap.** Windows exceeding a configurable variant-count
  threshold (`CRISPRME_IUPAC_CAP`) fall back to a bounded procedure instead of full
  2ᵏ enumeration, so a single pathological window cannot dominate runtime or memory.
  Crucially, the bounded procedure still emits a **greedy minimum-mismatch
  representative** for the window — at each variant column it takes the allele that
  most lowers the mismatch count, which (mismatches being additive per column) is the
  exact argmin over all 2ᵏ combinations, i.e. the window's worst-case off-target.
  So a capped window **always surfaces at least one off-target row**; the cap trades
  exhaustive per-haplotype enumeration for a single conservative representative, never
  the whole region.
- **Density reporting + transparency.** Every window that triggers the cap is written
  to a `high_variant_density_regions.bed` sidecar (region span, variant count,
  carriers, and the full IUPAC protospacer), and each affected off-target carries a
  `High_variant_density_region` column in `integrated_results.tsv` noting that a
  greedy worst-case alignment is reported and additional haplotype alignments may
  exist, with the full IUPAC sequence so a user can dig into them. Nothing is silently
  truncated — the bound, the representative, and the alternatives are all auditable.

Together these keep genome-wide variant search tractable while making any bound
that was applied explicit and reviewable, and guaranteeing no region is dropped.

---

## 6. Functional annotation of off-targets

Each reported off-target is annotated with its genomic context by intersecting
its coordinates with a 4-column BED (`chrom  start  end  label`), where the label
is suffixed by its source. `resultIntegrator.py` buckets labels by suffix into
dedicated columns; CRISPRme+ ships an updated annotation set (ENCODE **SCREEN
v4**) and adds a **COSMIC** cancer-gene column:

- **GENCODE** — gene-model context (`exon`, `CDS`, `UTR`, `transcript`,
  `start_codon`/`stop_codon`; `intergenic` otherwise), plus nearest-gene name
  and distance.
- **DHS** — DNase I hypersensitive (open-chromatin) sites, labeled by tissue /
  organ system.
- **ENCODE SCREEN v4 cCREs** — candidate cis-regulatory elements: promoter-like
  (`PLS`), proximal / distal enhancer-like (`pELS` / `dELS`), and
  chromatin-accessible / TF classes (`CA-CTCF`, `CA-H3K4me3`, `CA-TF`, `CA`,
  `TF`).
- **COSMIC (Cancer Gene Census)** — whether the off-target falls in a curated
  cancer gene, tagged by confidence **tier** (Tier 1: extensive curated causal
  evidence; Tier 2: strong but less-curated) and documented **role**
  (`oncogene`, `TSG`/tumor-suppressor, `fusion`). This flag is particularly
  relevant for therapeutic and pre-IND assessment, where an off-target in a
  known cancer gene warrants scrutiny.

Annotations are managed as an enable/disable set and applied automatically to
every search; the shareable report includes a plain-language legend for every
annotation value.

---

## 7. Shareable off-target assessment report

CRISPRme+ auto-generates a **self-contained HTML report** at the end of every
search (`generate_report.py`), bundled as a ZIP with the underlying tables. The
report is designed to be opened by a non-specialist (e.g. a reviewer or
collaborator) with no software beyond a web browser: plots are inlined as base64
PNGs, the table and styles are inline, and there are no external references.

It contains: a run summary and a mismatch × bulge count matrix; a graphical
report of reference vs. variant off-target scores (CFD and CRISTA, under
multiple rankings); a reference-vs-population origin breakdown; a **recommended
validation panel** (a hybrid worst-case top-N shortlist selected by combining
CFD, CRISTA, and edit-distance floors, with the selection logic stated
explicitly); per-threshold **downloads** sharing one curated, spreadsheet-ready
column schema; a scrollable **top-1000** table with the functional annotations;
and the annotation legend of Section 6. Allele frequencies can be omitted
(`--no-maf`) for runs where they are not yet finalized, so the site set and
scores can be shared without misleading frequency values.

A guide with **more than one perfect genomic match** (0 mismatch, 0 bulge) has no
a-priori on-target — each is an equally-efficient candidate cut site, and a
perfect-match *off*-target is the highest-risk class. Every 0-mismatch site is
therefore forced to the top of the validation panel (never truncated), flagged in
a `Perfect_match` column, and called out in a prominent warning banner (red when
there are several, listing the sites; amber "presumed on-target" when there is
exactly one) — in both the report and the interactive web results page.

The report ZIP places `report.html` at the top level with all data files under a
`data/` subfolder, and its *Variants included* line states the genotyped panel
size and the number of SNPs and indels searched. These database counts are read
from a build-time `Dictionaries/registry_<vcf>/variant_count.json` sidecar
(`n_records` SNPs + `n_indels`) written when the dict-less index is built, so the
report needs no VCF access at report time (for the shipped combined 1000G+HGDP
panel: 106,664,924 SNPs). Reported allele frequencies are AC/AN over the genotyped
panel; a variant present in the panel but whose source allele frequency is exactly
0 (e.g. a secondary allele of a multiallelic site) is shown at a **display floor of
1×10⁻⁵** so it still renders on the log-scale plots — this is a plotting floor,
read as "present, frequency effectively 0", not a measured frequency.

## 8. Off-target scoring, assumptions and limitations

**Scoring models.** Candidate off-targets are scored with two independent,
previously published models, used **unchanged** from their original definitions:
**CFD** (Cutting Frequency Determination; Doench *et al.*, *Nat. Biotechnol.*
2016) and **CRISTA** (Abadi *et al.*, *PLoS Comput. Biol.* 2017). Both return a
value in [0, 1]; higher means more likely to be cut. They are reported
side-by-side because they can disagree, and the recommended validation panel is
deliberately model-agnostic (a site worst by *any* metric is included), so no
single model gates the shortlist. CRISPRme+ does not re-train or modify either
model; it applies them to the aligned protospacer+PAM of each candidate.

**Domain-of-validity caveat (important).** CFD was trained on **single-nucleotide
mismatches** and was not designed to score DNA/RNA bulges (insertions/deletions).
CRISPRme+ nonetheless reports a CFD value for bulge-containing alignments, which
is an **extrapolation beyond the model's training domain**; for gapped/bulge
sites, CRISTA (which models indels) is the more appropriate score, and both
should be read as relative risk indicators rather than calibrated probabilities.
The CFD/CRISTA threshold tiers in the report are **model-relative** (CRISTA's
cut points differ from CFD's because the two scores are on different scales).

**Assumptions.** (i) Results are relative to the chosen **reference assembly** and
its coordinates. (ii) The variant panel is only as representative as the input
databases — **1000G + HGDP is broad but not exhaustive**, and a variant absent
from the panel cannot generate a variant-created off-target. (iii) Allele
frequencies are **AC/AN over the genotyped panel** (see §3); a frequency-only
database contributes its reported population AF directly, so a single combined
number can mix a genotyped denominator with imported summary statistics. (iv)
Phasing is resolved per haplotype from the genotypes (confirmed vs putative,
§4); unphased or cross-phase-set haplotypes are reported as **putative**.

**Limitations.** CRISPRme+ does **not** model somatic/mosaic variants, copy-number
or large structural variants, epigenetic state beyond the supplied annotations, or
chromatin accessibility as a cutting determinant. **SNP+indel co-occurrence** —
by default the SNP and indel searches run as two independent passes (SNPs on the
IUPAC-enriched genome, indels on a plain-reference fake-indel genome), so an
off-target requiring **both** a nearby SNP **and** an indel in the same protospacer
is not reported. This limitation is **lifted by an opt-in, experimental gate**
(`CRISPRME_INDEL_SNP=1`, off by default; classic builds byte-identical): when
enabled, the build overlays SNP IUPAC codes onto the fake-indel genome and compiles
a phased indel genotype tier, and post-analysis reports CONFIRMED-cis (phased) /
PUTATIVE (unphased) indel+SNP haplotypes with per-sample carriers + joint AF.
The `Max_total_edits` value is a
**search cap on the variant-collapsed (IUPAC) genome**; individual variant-expanded
alignments may exceed it (the report surfaces the observed maximum). A reported MAF
of `1e-05` is a **display floor** for a source-AF of 0, not a measured frequency
(§3, §7). These are computational predictions and are **not a substitute for
experimental off-target validation** (e.g. GUIDE-seq / CIRCLE-seq / targeted
amplicon or rhAMP-Seq sequencing of the recommended panel).

**Validation (search-engine correctness).** The fast index-based search has been
verified against an **independent brute-force ground truth** — an exhaustive
dynamic-programming alignment implemented by a *different* method (and a
dependency-free Rust re-implementation as a second independent check), driven from
a benchmark registry (`test/benchmark/benchmarks.json`) via `crisprme.py
validate-test`. This is a one-time correctness check, not a per-search step:

- **SpCas9 (sg1617, NGG, mm 4 / DNA-bulge 1 / RNA-bulge 1)** on the hg38 + 1000G
  variant-enriched genome — CRISPRme reproduces the brute-force reference
  **exactly (3,495 off-targets, 0 differences)**.
- **enAsCas12a (HBG1/HBG2 clinical guide, TTTV 5′ PAM)** on chr22 — **953
  off-targets**, matched against the brute-force reference.
- The registered benchmarks pass end-to-end (`complete-test` → `validate-test`),
  including the default single-`--max-total-edits` web-mode cases (4/4 locally);
  `validate-test` exits non-zero on any mismatch, and the round-trip runs in CI.
- Separately, the high-variant-density greedy decomposition (§5) was validated as
  provably exact against brute-force argmin on **4,000/4,000 random cases** (both
  strands, with/without bulges), with PAM-creating-variant cases reproducing full
  enumeration (variant attribution identical).

This establishes that the engine **does not miss** off-targets relative to
exhaustive search. It does **not** validate the *scoring* models' predictive
accuracy against experimental cleavage — a regulatory-grade package should add a
retrospective comparison of CFD/CRISTA ranking to experimental off-target assays
(e.g. GUIDE-seq / CIRCLE-seq / targeted amplicon or rhAMP-Seq).

---

*Software: CRISPRme+ (`pinellolab/crisprme-plus`). This document tracks the
methods as of the 2.4.0 line; see the CHANGELOG and the referenced source files
for implementation detail.*
