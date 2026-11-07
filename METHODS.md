# CRISPRme+ — Methods

This document describes the methods introduced in **CRISPRme+** (the 2.3 through
2.8.x line; current shipped release 2.8.1), intended as a self-contained technical reference and as source material
for the Methods section of the manuscript. It focuses on what is **new or
changed** relative to the original CRISPRme (Cancellieri, Zeng, Lin et al.,
*Nature Genetics* 2023); the core enumeration of candidate off-targets by
CRISPRitz and the primary CFD scoring are unchanged unless stated otherwise. The
machine-learning off-target score is **CRISPR-Bulge** (Yaish & Orenstein,
*Nucleic Acids Res.* 2024); see §8.

**Licence.** CRISPRme+ is released under the **MGB Open Access License 1.0** — free
for non-commercial, non-revenue-generating academic use; commercial use requires a
separate commercial licence (for **both** CRISPRme+ and stock CRISPRme). See
[README §6 License](README.md#6-license).

Sections:
1. [Variant-aware, dictionary-less data model](#1-variant-aware-dictionary-less-data-model)
2. [Homogenization and merging of population VCFs](#2-homogenization-and-merging-of-population-vcfs)
3. [Allele-frequency estimation](#3-allele-frequency-estimation)
4. [Haplotype scanning: observed-haplotype enumeration for phased, unphased, and mixed data](#4-haplotype-scanning-observed-haplotype-enumeration)
5. [Search-space control for high-variant-density regions](#5-search-space-control-for-high-variant-density-regions)
6. [Functional annotation of off-targets](#6-functional-annotation-of-off-targets)
7. [Shareable off-target assessment report](#7-shareable-off-target-assessment-report)
8. [Off-target scoring, assumptions and limitations](#8-off-target-scoring-assumptions-and-limitations)
9. [PAM-less indexes and dynamic PAM filtering](#9-pam-less-indexes-and-dynamic-pam-filtering)

Throughout, "protospacer window" means the genomic interval spanned by a
candidate off-target's protospacer plus PAM (and any bulges), i.e. the interval
in which an overlapping genetic variant can create, destroy, or modify an
off-target site.

**Terminology (used throughout).**
- **cis** — two variants are *in cis* when they lie on the **same DNA molecule** (the
  same haplotype / chromosome copy) in one individual, as opposed to *in trans* (on the
  two homologous copies). A co-occurring off-target requires its variants **in cis**.
- **observed haplotype** — a variant-set that occurs in **at least one real individual**
  (§4), as opposed to a synthetic worst-case reconstruction assembled from the union of
  all carriers.
- **CONFIRMED** — phased genotypes prove the variants are carried together in cis (exact
  carriers + joint AF). **PUTATIVE** — co-occurrence is possible but cis is unproven
  (genotyped-unphased panels still name the co-carriers; sites-only panels give a
  conservative min-AF bound with no carriers). See §4 for the three-rung model.
- **worst-possible representative** — a single per-window off-target row that
  upper-bounds the score/severity of every realizable haplotype in that window, emitted
  by the population-level analysis instead of enumerating each haplotype (§5).

**Analysis-mode decision rule.** The **default** (no flag) is **context-aware**:
per-sample (CONFIRMED cis phasing + named carriers + exact joint AF) on a
genotyped+PHASED index, or a population-level worst-possible screen (works on any
index, lossless for detection) otherwise (sites-only/unphased/hybrid/reference-only).
Force the fast population-level mode with **`--population-level`**; force resolution
with **`--per-sample`** (mutually exclusive) — the latter is inert on a sites-only
index (`mega`). The genotyped+phased panels are `HGDP1kGP` — the recommended default,
catalog name `NRG_3_hg38+hg38_HGDP1kGP`: 925 HGDP + 3,166 1000 Genomes = 4,091 samples,
jointly SHAPEIT5 re-phased — or `HPRC`. See
[§5, Search-space control for high-variant-density regions](#context-aware-analysis-mode-per-sample-genotype-resolution-or-population-level-analysis).

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
Physically, the registry is three contiguous sections:

- a fixed-width **record array** — 16 bytes per record (position, reference base,
  alternate base, the record's group count, and offsets into the two following sections);
- a **group blob** — the per-group AC / AN / carrier / homozygote / called-allele counts;
- a deduplicated **string pool** — rsIDs.

The group blob dominates the file (~75–80 %), and within it the allele-number and
called-allele columns are near-constant across the millions of records on a chromosome
(the genotyped panel is essentially the same size everywhere), so the raw layout is
highly redundant.

**Raw layout (default).** CRISPRme+ ships the registry **uncompressed (raw)**: the three
sections are memory-mapped and read by direct `O(log n)` bisection on genomic position,
with **no per-lookup decompression**. This is the fastest option — a dense variant search
does millions of registry lookups (one per candidate off-target window, scattered across a
chromosome), and measurement showed raw lookups are **~2× faster** than the compressed path
even with a warm cache. All shipped indexes (the genotyped 1000 Genomes + HGDP and the
sites-only mega) use this raw format, for consistency and lookup speed ("speed over space").
The registry has shipped raw since 2.5.3 (still raw as of 2.8.1).

**Optional `zlib` compression.** A block-compressed format is available for
space-constrained storage. Each section is partitioned into 4,096-record blocks,
individually compressed (~3.5× smaller on disk), with a sparse block index in the manifest.
The reader is fully backward-compatible and reads either format byte-for-byte identically.
When a compressed registry is used, a decompressed-block **LRU cache** amortizes
decompression to one pass per block: profiling had shown block `zlib` decompression reaching
**~71 %** of a dense post-analysis with a tiny cache, and the right-sized cache plus an O(1)
LRU cut a dense chr22 SNP post-analysis **326 s → ~89 s (~3.66×)** (measured for the
8→512-block cache step; the speedup is cache-dependent, and the default 4,096-block cache
extends it to the largest contigs) — but the raw default avoids that cost entirely.

| Knob | Default | Effect |
|---|---|---|
| `CRISPRME_REGISTRY_COMPRESS` | `0` (raw) | set `=1` at build time to write the `zlib` block-compressed format |
| `transcode_registry` | — | re-encode an existing raw registry losslessly to compressed (or back), no VCF re-parse |
| `CRISPRME_REGISTRY_CACHE_BLOCKS` | `4096` | decompressed-block LRU cache size (each block ~64 KB; fills only to the blocks a search touches) |

### Out-of-the-box variant search
Because Tier-0/Tier-1 are small, they are shipped **with the pre-built index**
(the Tier-0 registry raw/uncompressed by default, memory-mapped and read on the
fly; zlib block-compression is opt-in). A user who downloads a variant index can
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

  *Known limit (engine).* The enricher writes one IUPAC code a position and one
  SNP entry a position, rebuilt from the reference base by each record, so at a
  position with several biallelic SNP records only the **last** record's ALT is
  searched: an earlier record's ALT is in neither the enriched genome nor the
  dictionaries nor the Tier-0 registry nor the Tier-1 genotype store, and an
  off-target that needs it is not found, or is reported at a lower score, without
  its carriers (a VCF line that lists several ALTs itself is coded with all of them).
  It is not new in 2.8.1: every release behaves the same (the SNP dictionaries have
  always kept one record a position, the enriched genome since CRISPRitz 2.6.0).
  On the HGDP + 1000 Genomes panel, 2,486,165 of its 127,666,429 SNP records
  (1.95 %), at 2,441,139 positions, are not coded for this reason. Per guide (NRG, 4
  mismatches + 1 DNA + 1 RNA bulge, per-sample, genome-wide) this leaves 3 to 7 sites
  whose true CFD puts them in the top-1000 out of it or understated there, and at most
  one site reaching CFD ≥ 0.5 in the five guides tested (the README's, ACTGA, TRAC,
  EMX1 and VEGFA site 3; the one is EMX1's chr17:12179615, 11 haplotypes, CFD 0.508);
  and 52 to 324 rows per guide (99 % of them at CFD < 0.1) name an individual who
  carries such an ALT, read as the reference base, for a sequence they do not have.
  On chr22, in the normalized source VCF: 45,837 of 1,757,775 SNP records (2.6 %), at
  43,389 positions, are not the last at their position; 5,104 of them have AC ≥ 41
  (AF ≥ 0.5 % of 8,182 haplotypes) and 1,566 AC ≥ 818 (about 10 %). Checked directly
  at six such positions (the three with the highest AC on chr22, two at AC 82 and 257,
  and chr20:12296886): each earlier ALT, up to AC = 8,179 (`rs807755` T>A at
  chr22:19142014), is absent from the enriched genome's code and from all three
  stores, and only the last record's ALT (AC 1 to 3) is held. It will be fixed, with
  rebuilt indexes, in the next release.
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
AN (here 2×(2,548 + 929) = 6,954 for combined 1000G+HGDP autosomes). (This
illustrative example uses the earlier 1000 Genomes 2019 + HGDP merge; the shipped
default is the re-phased HGDP1kGP, 4,091 samples.)

### Two panel modes (genotyped vs sites-only)
The merge above is the **genotyped, cis-capable** mode (`merge_vcf_panels.sh` /
`build_combined_panel.sh`): it keeps per-sample genotypes, recomputes a **pooled**
`INFO/AF` (AC/AN over the union of genotyped samples, `bcftools +fill-tags`), and
feeds the genotype-counting Tier-0 registry. It is correct only when *every* merged
source is genotyped, and it is what enables indel+SNP cis reconstruction.

The **default** shipped index, **`NRG_3_hg38+hg38_HGDP1kGP`**, is built in this genotyped mode
from the gnomAD **v3.1.2 HGDP+1kGP** jointly-called callset (GRCh38; 4,091 samples = 925 HGDP +
3,166 1000 Genomes), which CRISPRme+ **re-phases with SHAPEIT5** (`phase_common` scaffold +
`phase_rare`). Re-phasing deliberately **keeps singletons** — the public `phased_haplotypes_v2`
release drops variants at MAC ≥ 2, discarding exactly the rare alleles most relevant to a
conservative off-target screen. Because *both* cohorts are co-phased, the whole panel is
**genotyped-phased**, so HGDP off-targets report CONFIRMED cis with named carriers (the earlier
`1000G2021_HGDP` hybrid, which mixed phased 1000G with *unphased* HGDP, could only call the HGDP
portion PUTATIVE). chrX uses a 3-way PAR split with non-PAR **males treated as haploid**
(`bcftools +fixploidy`), giving the correct allele-number denominator (non-PAR AN = 5,990, PAR
AN = 8,182). The full re-phasing protocol, chrX ploidy handling, and validation are documented in
`docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`.

CRISPRme+ also ships an **HPRC pangenome panel** (`hprc_build.sh`) in this same
genotyped, cis-capable mode. It is a *single* source rather than a merge: the HPRC
Minigraph-Cactus pangenome — the **HPRC Release 2 graph `hprc-v2.0-mc-grch38`** (GRCh38) —
decomposed to a per-chromosome VCF (`vg deconstruct` → `vcfbub` → `vcfwave`), 232
assembly-derived genomes (231 individuals + the CHM13 graph backbone). Because the genotypes are **phased** (each sample is two assembled
haplotypes), it supports CONFIRMED cis co-occurrence and named per-sample carriers
like a short-read 1000 Genomes + HGDP panel, while capturing pangenome-specific variation
absent from short-read call sets. The graph-coverage half-missingness (`1|.`) is
handled by recomputing `AF/AC/AN` from the genotypes (`bcftools +fill-tags`), so the
allele-number denominator excludes unassembled alleles.

> **Assembly-search (a complementary, VCF-free approach).** For an individual whose
> genome is fully assembled, CRISPRme+ also offers `assembly-search`, which sidesteps
> variant inference entirely: it searches the two haplotype assemblies (paternal +
> maternal FASTAs) **directly**, lifts each haplotype's off-targets to hg38 via a
> liftOver chain, and reconciles them — a site found on both haplotypes is
> homozygous-equivalent, on one is heterozygous-equivalent, and one with no hg38
> equivalent is *haplotype-non-mappable* (invisible to any reference-based search). A
> haplotype with **no** off-targets at all is a valid outcome, not an error — it
> contributes nothing, so every locus on the other haplotype is heterozygous-equivalent.
> This is the assembly-level counterpart to the vcfwave-derived HPRC index above: the
> index scans one enriched reference for population-level HPRC variation, while
> assembly-search scans a specific individual's assembled genome. It is available on the
> CLI (`crisprme.py assembly-search`) and the web interface.

CRISPRme+ also ships a **sites-only "mega" panel** (`merge_mega_sites.sh`) that
merges heterogeneous **aggregate** resources — 1000 Genomes 2021, HGDP, gnomAD v4.1,
TOPMed, and All-of-Us — where genotypes are unavailable or meaningless (gnomAD is
frequency-only, TOPMed distributes `AN=0`, All-of-Us is a single aggregate
pseudo-sample). Because there is no honest pooled AC/AN across such sources, and no
shared samples to reconstruct cross-source haplotypes, the mega:

1. **normalizes** each source (`bcftools norm -m -any -f REF`: split multiallelics to
   biallelic and left-align, so the same variant from two sources is represented
   identically and merges rather than duplicating; `AF` is `Number=A`, so the correct
   per-alt frequency is carried without reading genotypes);
2. applies a uniform **MAF > 0.001** filter and **strips genotypes** (`view -G`);
3. keeps each source's frequency verbatim as `AF_<source>` and, after
   `bcftools merge -m none`, annotates a per-site **`AF_max`** — the maximum
   `AF_<source>` at that site — as the global summary frequency (no pooled `AF`).

The mega's Tier-0 registry is built **directly from these frequencies**
(`compile_registry_from_info_af`, not from genotypes): each source becomes one
database group with allele count `AC = round(AF · AN_nom)` and `AN = AN_nom` (twice
the source's nominal sample size — so the reported allele frequency reproduces the
source AF exactly, to within `0.5/AN`, and the reported AN is the source's true
cohort size), and the GLOBAL group carries `AF_max`. Per-individual carrier and
homozygote counts do not exist for aggregate data, and none is estimated: the registry
stores 0 for both (the frequency itself is exact), and its manifest marks the panel as
aggregate (`aggregation` `info_af`, and `aggregate_af_only` for each dataset), so the
population-summary companion gives them as NA, not 0 (§7), while `targets-scoring`'s
carrier overlay prints the stored 0 as the carrier count. The Tier-0 registry stores SNPs
only (single-base ref/alt); indels (~24 % of merged sites) are surfaced as off-targets
through the fake-indel genome and **are** per-dataset frequency-annotated via a separate
`indel_af` sidecar (per-dataset `AF_<dataset>` + `AF_max`), which the report reads for the
mega panel.

The two modes are therefore complementary: the genotyped panel gives phased,
cis-capable frequencies over a curated sample set; the mega gives one-scan,
frequency-annotated coverage across the widest set of population resources.

Each registry manifest (`reg_<chrom>.idx`) records the index's **`data_type`** —
`sites-only`, `genotyped-unphased`, `genotyped-phased`, or `hybrid` — plus a
`phased` flag per database and overall, observed from the genotypes at build time.
This lets tooling (e.g. the web variant-dataset selector) identify an index's type
and phasing without scanning the multi-GB per-sample dictionary; a missing flag on
an older index falls back to genotype-store presence. The shipped indexes are
`HGDP1kGP` = **genotyped-phased** (the re-phased default, above), `HPRC` =
**genotyped-phased**, and `mega` = **sites-only**;
`hybrid` is the fourth possible value, for a merged panel mixing phased and unphased
sources.

### Re-indexing for a different PAM

An enriched reference genome is **PAM-independent**: the per-position IUPAC overlay and the
Dictionaries/registry encode the *variants*, not the nuclease. Only the CRISPRitz search index
(`genome_library/<PAM>_<N>_<ref>+<vcf>`) is PAM-specific. CRISPRme+ therefore publishes the
PAM-independent **enriched-genome kit** for the flagship datasets (under `build-inputs/<dataset>.tar.gz`
on the HuggingFace data repo), so a user targeting a different nuclease (NGG, NNGRRT, TTTN, …)
rebuilds *only* the index off the published enriched genome rather than re-deriving the variants.
The complete workflow — downloading the kit, indexing a new PAM, and reproducing the kit from
source (including the re-phasing above and the mega merge/filters) — is in
`docs/PRECOMPUTED_INDEXES.md`.

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

The Tier-0 registry holds SNPs only. An indel has, in either mode, the `AF` of the
index's indel log, which CRISPRitz `add-variants` copied from the VCF's `INFO/AF` when
it built the index, as rounded there; one whose `AF` the log gives as 0, or not at all,
has no frequency. The log holds none for a record whose `INFO` has `AF` at another place
than the VCF's first record (a known limit of `add-variants`, §8), and none for any
indel of the sites-only mega, whose VCF has no `INFO/AF` (only `AF_<dataset>` and
`AF_max`, §2): a mega indel's per-dataset frequencies are in the `indel_af.tsv`
companion.

For a multi-variant haplotype off-target, the reported frequency is bounded by
the **rarest** contributing allele (the haplotype can be no more frequent than
its least-frequent variant), and the exact carriers of the haplotype (Section 4)
give the tightest available estimate.

An index without the dict-less tiers (one built without `--samplesID`, such as the one the
benchmark `complete-test` builds) has no Tier-0 registry: a SNP has the `AF` of the classic
per-sample dictionaries, which `add-variants` also copied from the VCF's `INFO/AF`, as
rounded there, with the same limit as an indel's (§8). The benchmark's 1000 Genomes VCF
gives two decimals, so a rare variant reads 0 and has no frequency, and a multi-variant
row's minimum leaves it out. On the NGG index the README's self-test builds
(`NGG_2_hg38_chr22+hg38_1000G`, Section 5's known limit), a `--per-sample` search of the
README's guide with the `20bp-NGG-SpCas9` PAM, 4 mismatches, 1 DNA and 1 RNA bulge gives 154
of its 189 variant rows a blank MAF (3 of 184 with the tiers), and four two-variant rows the
frequency of their commoner variant (0.09 to 0.86) for a haplotype that one to five of the
VCF's 2,548 individuals carry. The report of a search on such an index says so, in either
mode; rebuild the index with `build-index-only --samplesID`, or use a prebuilt one.

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
**excluded** by construction. A SNP row's target is the reference with the SNP alleles
of a haplotype (the enriched genome holds SNPs only), so it is not the sequence of a
haplotype that also holds an indel changing the target's bases: that haplotype's own
target there is the indel row's, when it is in budget. Since 2.10.0 a per-sample SNP row
does not list such a carrier. The carriers are read against the phased indel genotype
store the indel post-analysis reads (§4, SNP+indel), as an indel row's second indel is:
- an indel changes the target when it deletes or replaces a base of the target's span
  (protospacer and PAM, a DNA-bulge base included) or inserts bases between two of them.
  One whose deleted bases all lie outside the span, its junction at or past the target's
  edge, or that inserts before the target's first base or after its last, leaves the
  target's bases as they are;
- whether the bases change is read on the haplotype's own bases, written with every indel
  of the store it holds there and its own SNP alleles: in a tandem repeat a unit gained or
  lost at the target's edge leaves them as they are (from the first base on, or up to the
  last), and the carrier stays. The bases must be read through the SNPs the row names: a
  haplotype whose indel pushes a named SNP out of the stretch, or deletes its base, has its
  target thanks to the indel, not the SNP;
- a phased carrier is decided by its haplotype: one whose indel is on the other haplotype
  stays. One homozygous for an indel that changes the bases holds it on every haplotype and
  leaves the row (one homozygous for an indel that leaves them as they are stays, as above).
  Otherwise one whose indel call or SNP call is unphased cannot be tied to a haplotype and
  stays, PUTATIVE, as every carrier on an unphased panel is. An indel at the position of a
  SNP the row names is another allele of that base (a VCF writes a haplotype holding both as
  one record), so a haplotype phased with both is a phasing artefact, and it stays; a
  haplotype whose bases cannot be read (an indel there whose slot is unknown, two records
  of one base) leaves the row;
- a row left with no carrier is not reported. The population-level analysis (§5) is
  unchanged: its rows are worst-case composites. With `CRISPRME_INDEL_SNP=0`, the opt-out of
  the SNP+indel co-occurrence, the store is not read, and a SNP row lists its carriers as up
  to 2.8.1.

On chr22 of the HGDP + 1000 Genomes panel, the README guide, NRG, 4 mismatches + 1 + 1, per
sample, up to 2.8.1 3 of 628 highest-CFD SNP rows of the integrated results and 7 of 399 of
the alternative alignments listed such carriers (4,234 listings), and 7 of them had no other
carrier: chr22:21052638 − listed 1,292 carriers of 21052641 G>A, which sits on the same 1,448
haplotypes as 21052640 T>TCACC. Genome-wide, 280 of 65,184 highest-CFD SNP rows (28,700
listings, 13,525 of them carriers homozygous for the indel), and 208 had no other carrier (none
at CFD ≥ 0.2, highest 0.112).

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
  individual could plausibly carry — bounded by `CRISPRME_IUPAC_CAP` (default **10**): an individual
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

### Sites-only panels and co-occurrence without genotypes
On a **sites-only / aggregate panel** (the mega index, or any download without a
genotype tier) there are no per-sample genotypes to reconstruct cis, but the
variants co-located in a window still form **putative** haplotypes. CRISPRme+
reports these too, forming a three-rung confidence model:

1. **CONFIRMED** — phased genotypes prove the variants are carried together in
   cis (e.g. 1000 Genomes); exact carriers and joint frequency.
2. **PUTATIVE (co-carrier)** — genotyped but unphased (e.g. HGDP); the individuals
   who carry all the variants are known (a both-carrier count), but cis is unproven.
3. **PUTATIVE (estimated)** — sites-only; no genotypes at all. The variants are
   known to segregate in the population at their marginal allele frequencies, so the
   joint cis frequency is reported as a **conservative upper bound = the minimum
   participating marginal AF** (a cis haplotype can never be more frequent than its
   rarest allele; no LD assumed).

This applies to **both** co-occurrence dimensions: **SNP+SNP** (an off-target that
requires ≥2 nearby SNP alt alleles together; `snp_snp_cooc.tsv`) and **SNP+indel**
(`indel_snp_cooc.tsv`). On a genotyped panel both are emitted CONFIRMED/PUTATIVE with
carriers; on a sites-only panel both fall to rung 3 (PUTATIVE, min-AF, no carriers) —
so the co-occurrence signal is never silently dropped for lack of genotypes.

### SNP+indel off-targets: each carrier's own haplotype
Off-targets that need an indel are searched on the **fake-indel genome**: one short window
per indel (25 bp of reference, the indel's ALT, then reference), onto which the build
writes the IUPAC code of every SNP of the enriched genome. The search engine matches a code
when **any** of its alleles fits the guide (or the PAM), so one raw indel target stands for
every combination of the window's SNP alleles at once; only the genotypes say which of them
exist. The indel post-analysis therefore resolves each such target **per carrier
haplotype** (`indel_snp_cis.resolve_target`):

- every haplotype that carries the indel gets **its own allele** at each code column — the
  SNP alt it carries, else the reference base — and the haplotypes left with the same
  sequence form **one reported row**, listing only their individuals and naming the SNP alts
  they hold next to the indel (Variant info, rsID, MAF). A SNP with several alts on one VCF
  line keeps its raw genotype under each alt (2|1 under A and under T), so each genotype is
  first read against its own alt; where the genotypes cannot tell which index is which (every
  carrier of two alts is 1|2) no slot can be tied to the alt and the carrier is resolved as an
  unphased one;
- a sequence whose PAM no longer satisfies the searched motif, or whose mismatches,
  recounted on the concrete sequence, exceed `--mm` or `--max-total-edits`, is **dropped**:
  a raw target that no individual carries gives no row;
- an individual whose genotypes cannot put an allele on the indel's haplotype (an unphased
  call, a missing allele, no phased indel genotype, a haploid call against a diploid one)
  gets the **strongest target its own alleles allow** on each alignment (the guide's base
  wherever it carries it and the alignment pairs the position with a guide base; elsewhere, a
  DNA bulge included, and at the PAM jointly, the alleles that maximize CFD among those the
  motif allows; a search that scores no CFD, with a 5′ PAM or a guide or PAM not 20 and 3 nt
  long, takes there the reference base where its alleles and the motif allow it, else the first
  of A, C, G, T they allow), PUTATIVE in the co-occurrence companion: it is never listed for an
  allele it does not carry. The carriers of an equal-length multi-base record (an MNV, AC>TG),
  which CRISPRitz puts in a fake-indel window as it does an indel, are always in this case: the
  phased indel genotype store keeps only records whose REF and ALT differ in length. With
  **no genotypes at all** (a sites-only panel such as `mega`) the row
  holds that worst case over the code's alleles, keeps the carriers the indel log lists,
  and its co-occurrence is rung 3 above (PUTATIVE, min-AF bound); with no SNP genotypes
  (no Tier-1 store and no classic dictionaries, as after
  `download --what index --no-genotypes`) the row is that worst case too, under every
  carrier of the indel the log lists, PUTATIVE;
- the fake-indel window holds **one indel**, so a haplotype that carries a second indel
  inside the target may not have the window's sequence there, and whether it has is read on
  its bases, not on the indel's position: the reference around the target is written out
  with every indel of the phased store the haplotype carries there and its own SNP alleles,
  and the haplotype is listed when the target's bases are on it at the target's place, from
  the target's first base on or up to its last base, through the indel (not in a stretch
  that another of its variants lengthened, pushing the indel out of it). In a tandem repeat
  a unit gained or lost at the target's edge leaves the bases as they are (the repeat puts
  the same bases back), so such a carrier keeps its row; one whose second indel changes the
  bases is not listed for it, on any indel row (§8, residuals), nor one whose bases cannot
  be read either way: each way meets an unphased or missing SNP call, an unphased indel
  that does not itself change the target, or two of its variants that overlap, two records
  of one base: two other indels whose REFs share a base (any two
  at one position do, even a duplication and the deletion of the same bases, which together
  would leave the reference), or a SNP on a base of the second indel's REF. A reading from
  the target's first base on meets no variant upstream of that base, one up to its last base
  none downstream of it. A second indel inserted right next to the indel's anchor base is
  inside the target. A carrier whose second indel inside the target has an unphased call
  cannot be tied to a haplotype, so it stays listed, PUTATIVE, as every carrier on an
  unphased panel is; this holds even when another, phased second indel on its indel
  haplotype changes the target, so such a carrier is listed for a row its indel haplotype
  does not have (§8, residuals).
  On the same phased slot, an indel overlapping the first at another position (one complex
  allele the VCF splits into overlapping records; the window keeps the reference there)
  takes the haplotype off the target's rows whatever its bases. A variant only on bases the
  indel deletes is not in the window's sequence and changes nothing. Another allele at the
  indel's own position, an indel or a SNP, is an alternative to it (a haplotype holding both
  would be one record), so a slot phased with both is a phasing artefact. And an unphased
  carrier of two overlapping records is a multi-allelic heterozygote, one allele per
  haplotype. These three stay listed;
- a target of the window that holds **no base of the indel** — none the ALT adds or replaces,
  and not the two bases on either side of the deleted ones — lies wholly in a flank and gives
  **no indel row**, with or without a code: its sequence is reference (or SNP) sequence, which
  the reference and SNP searches report at its true coordinate for its true carriers, while
  the indel shifts its coordinate on the window;
- a classic dictionary install (no Tier-1 store, e.g. an index built without `--samplesID`)
  reads the same genotypes from the per-sample SNP dictionary; with `CRISPRME_INDEL_SNP=0`
  at search time every code reads as the reference base, as on a fake-indel genome built
  without the overlay.

The co-occurrence companion `indel_snp_cooc.tsv` lists the SNP alts the indel
post-analysis's rows **require** — the guide's own base at a spacer column, or a base that
creates the PAM where the reference does not — with those rows' carriers: CONFIRMED when
every carrier was resolved exactly, PUTATIVE otherwise. It is written before the merge
picks each cluster's alignments, one line per engine alignment, required SNP set and phase
(an alignment that needs no SNP writes no line), so a line can name an alignment the merge
did not keep: a population-level run, which writes no alternative alignments, can have
lines with no row in its tables. The resolution does not depend on the analysis mode (§5).
A target with no code (every target of a classic, non-overlaid
fake-indel genome) that holds the indel needs no resolution: its sequence is the one the
engine found. An indel row's carriers belong to its sequence and its alignment: two alignments
of one site can name the same variants with different carriers (an RNA-bulge alignment does not
reach the base a SNP sits on, so it lists that SNP's carriers under the indel alone; and a
PUTATIVE carrier can have different bases on two alignments of one target string, the guide's
base where one pairs the position with a guide base and the CFD's choice, or the reference base
in a search that scores no CFD, where the other leaves it as a DNA bulge), and the merge of a
cluster's rows, which otherwise keeps one row per position and variant list, keeps one per
carrier set there, so that no carrier is added to a sequence it does not have, nor to an
alignment it was not resolved on (the merge keys an indel row by its aligned spacer as well as
by its target, since 2.10.0: up to 2.8.1 the alignments of one target string were one entry
whatever their carriers, and in the SNP+indel oracle's site U the unphased S19 was listed with
the reference base on the alignment where its alleles give the guide's base); of the alignments
with one carrier set it keeps the best by the order that picks each block's representative (the
score first in the CFD and CRISPR-Bulge blocks, then `--sorting-criteria`). Carrier sets are
compared whole, so that every block keeps one row count whatever its order. So an alignment
whose carriers are all listed on a better alignment of the same target and variants still keeps
an entry when the better one also lists a PUTATIVE carrier resolved to another base on the
weaker one: with a SNP whose base the weaker alignment leaves as a DNA bulge, an unphased
carrier takes the guide's base on the better alignment and the CFD's choice on the weaker, so
only the carriers homozygous for the SNP, or phased with it in cis, are listed on the weaker
one with the SNP. That entry is an extra alternative alignment, its carriers and score right;
none occurs in the chr22 searches checked on HGDP + 1000 Genomes, HPRC and an unphased copy of
HPRC. An MNV row is merged as an indel row is.

An indel row is placed as the reference and SNP rows of the same alignment are: it starts at
its target's first aligned base, the leftmost on either strand, carried from the window to
the chromosome, and its REF alignment holds the reference bases the target covers from
there, in the target's orientation and gapped as the target (an RNA bulge is a gap, so the
REF covers one base fewer per RNA bulge). Up to 2.8.0 an indel row with only RNA bulges
started one base per bulge early on the + strand (on the - strand with a 5' PAM), and the
REF alignment of an indel row with both bulge types had no gap and one base too many.

**A row that starts inside an insertion.** A first base that an insertion adds has no reference
coordinate: such a row starts at the reference base the inserted bases follow, the insertion's
anchor (its VCF POS, 1-based; POS − 1 in the 0-based `Start_coordinate`), as the VCF writes the
insertion and as the target's reference span already counts it. Up to 2.8.1 it started at that
base's place on the haplotype, 1 to len(ALT) − len(REF) bases past the anchor, where its REF
alignment was read, from reference bases the target does not hold. Alignments that start on
bases of one insertion therefore share the anchor. The merge takes a guide's rows on one
chromosome for one site while the cluster position of each (its start, less its RNA bulges and
plus its DNA bulges on the + strand with a 3′ PAM and on the − strand with a 5′ one) is within
`--merge` bases (3 by default) of the one before, on either strand and whatever their bases,
and keeps one row to represent the site, the others being its alternative alignments:
alignments that start on bases of one insertion are now one site, and such a row can join the
site of a reference or SNP row that starts near the anchor (on the HGDP + 1000 Genomes index,
the README guide's chr22:21052642 (−) of 2.8.1, 4 bases from the reference row at 21052638 (−),
is now 21052639, one of that row's alternative alignments). Within a site, the alignments that
name the same variants with the same carriers are one entry, the best, as on every row: they
are one target's alignments with different bulges. The merge compares indel rows there at their
cluster position on the haplotype carrying the indel, which on a row that starts inside an
insertion is past the anchor by the rank of its first base among the inserted ones, read from
the row's bases (`indel_snp_cis.row_first_base`). Two alignments whose cluster positions on
that haplotype differ are two targets of it, as many bases apart as 2.8.1 placed them, and both
are listed, as the site's row or one of its alternative alignments (a population-level run
lists only each site's row); one target's alignments with different bulges before one PAM fold
into the best, as up to 2.8.1, also when they start on different bases an insertion adds. Two
copies of one target inside an insertion (a repeat it holds) make the same row at the anchor,
its bases, REF alignment and variants the same, and that row is listed once, with the carriers
of both, where 2.8.1 listed each copy at its own start. A bulged alignment of the second copy,
which 2.8.1 folded into that copy's row, then has no row of its own copy to fold into and is
listed as one of the site's alternative alignments (the SNP+indel oracle's X2). With the TRAC
guide (CTCTCAGCTGGTACACGGCA) on that index, two alleles of chr22:50193929 that insert 96 bases,
a 48-base unit twice, with one carrier each, hold the CFD-0 alignment CTACTCAGgaGGgACA-GGCtGAG
(−) in both copies, which 2.8.1 listed at 50193950 and 50193998 and 2.10.0 lists once per
allele, at 50193928. Such a row's ALT alignment begins mid-insertion: on the + strand its first
bases are the last bases the insertion adds that the target holds (TTAT, the 7th to 10th of
G>GCAATCATTAT, in a 2.8.0 alternative alignment at chr22:22322768). Its REF alignment holds the
reference bases from the anchor on, as many as the target holds, so it begins with the anchor
base, which the target does not hold, and past the insertion it is out of register with the ALT
by the number of inserted bases the target holds, less one. On the − strand both are shown
reverse-complemented, in the target's orientation: the ALT ends mid-insertion (with a 3′ PAM,
its PAM holds those inserted bases) and the REF ends with the anchor's complement.
`CFD_score_REF` and the REF seed and non-seed counts are read from that REF alignment.

**`Variant_info_spacer+PAM`.** Each variant is numbered as on a SNP row, by its distance from
the target's first aligned base to the variant's POS base, the first base being 1: on
the + strand counted from the left, POS − start; on the − strand from the right,
len(`Aligned_spacer+PAM`) + 1 − (POS − start), less the deleted bases for a deletion, its REF
and ALT reverse-complemented. The number shown is its absolute value: a deletion that
reaches past the target's far end makes it negative (TRAC guide, chr22:48378297 (−), 21
bases deleted at 48378315: 24 + 1 − 18 − 21 = −14, shown 14, as in 2.8.1). A variant
before the target's first base is 1 on the + strand and
the last on the − strand. On an indel row the distance is counted along the haplotype carrying
the indel, the target's own bases: the bases an insertion adds count, those a deletion removes
do not. So a SNP past an indel is len(ALT) − len(REF) bases further along the target than along
the reference, and a row that starts inside an insertion is counted from its first base, the
inserted one, not from its start: in the SNP+indel oracle's site Z3 the SNP is the target's
13th base. The indel itself is numbered at its POS base, the anchor, which comes before the
bases it changes: on a row that starts inside an insertion it is 1 on the + strand and the last
on the − strand. Up to 2.8.1 the distance was counted along the reference: a SNP past an
insertion was len(ALT) − len(REF) bases short of its place on the + strand (long on the −
strand), one past a deletion as many long (short). So was a SNP on a row starting inside an
insertion, which 2.8.1 placed at its first base's place on the haplotype (Z3: 4 for 13), and
one whose POS was before that place was taken for a variant before the target's first base, 1
on the + strand and the last on the − strand (chr22:21052639 (−) above: its SNP, numbered 22,
was numbered 24). On the − strand each RNA bulge adds one to the numbers, on SNP rows as well:
chr22:21052639 has one, and its SNP is the target's 21st base.

**`PAM_creation`.** It holds the ALT's PAM when the variants created it and NA when the
reference holds a PAM there. On a SNP row the PAM sits on reference bases, and the predicate
reads the reference bases under it. On an indel row the PAM's place on the reference depends
on which of the target's bases are the indel's, and in a repeat that is not fixed: the VCF
writes an insertion or a deletion at the left end of its repeat (bcftools norm left-aligns it),
and the same haplotype holds the same bases with the inserted sequence rotated through the
repeat and written further right, or with the deletion's junction moved through it. Since
2.10.0 the PAM is judged on every placement of the indel that keeps the row's haplotype, the
reference with the row's own variants (its indel and the SNP alts it names), so the answer
does not depend on where in its repeat the VCF writes the indel:

1. a placement that leaves the whole PAM on reference bases (no inserted base, no deletion's
   junction inside it): the reference bases it sits on there. Without a SNP alt of the row in
   the PAM they are the PAM's own bases: NA. A SNP alt is read as the reference base it
   replaces, as on a SNP row; a PAM that two placements put on the reference (before the
   indel and past it, when the repeat spans it) is not created if either reading holds the
   motif;
2. else a placement that leaves the protospacer's base that touches the PAM on a reference
   base: the reference PAM beside that base, on the PAM's side, the PAM's place beside that
   protospacer. A deletion's PAM always has one, its bases all being reference bases; an
   insertion's has one when some placement puts the inserted bases at or past the PAM's
   first base (the PAM at the target's right end, forward: the + strand with a 3′ PAM, the −
   strand with a 5′ one) or ending at or before the PAM's last (the PAM at its left end);
3. else, an insertion that holds that base, with the protospacer's end and the PAM's start,
   wherever it is placed: no place on the reference, created.

A row whose REF alignment is its ALT has its target, PAM included, on the reference at its
start: NA (since 2.10.0). A PAM with a bulge in it keeps the REF alignment's bases, and so does
a PAM holding an ALT base that stands in a REF base's place, an MNV's or a complex record's, as
up to 2.8.1 (for an MNV, which moves no base, those are the reference bases under the PAM). The
rule reads
the row's haplotype: a carrier's other indels are not part of it (a carrier whose second indel
changes the target's bases is not listed on the row, above), and a repeat that runs past the
target is read on the reference.

Examples, SNP+indel oracle (`test/indel_snp_oracle`): a tandem duplication of a 27-base unit
holding a whole target, the copy written left of the reference's (PR, PM): the copy's target
starts inside the insertion, its PAM wholly inserted bases, and placed one unit right the
insertion leaves that PAM on the reference's own: NA. Where an insertion after the protospacer
copies the reference's own PAM (T, C6), NA likewise. C>CCAG after its anchor on the − strand
(QC), whose first base the reference repeats there, holds the protospacer's first forward base
in both its places: created. Not a repeat, C>CCA after an anchor C that starts a − strand row
(QN): the protospacer's touching base is a reference base, and the reference holds CTC beside it
(GAG on −): NA; and with an inserted TGG before the reference's AGG on the + strand (T2), or
TTCCA inserted between the reference's CCG and its protospacer on the − strand (QW, a wholly
inserted PAM), the reference holds a valid PAM beside the very protospacer: NA. CA>C inside the
PAM CCG on the − strand (QD) leaves the reference's CAG (CTG on −) beside the protospacer:
created. Inserted bases holding the protospacer's last four bases and the PAM on the + strand
(QA): created.

In 2.8.1 a PAM on reference bases was read on them, through the indel, and any other on the REF
alignment's bases, read from the row's start. That answered the question only
where the protospacer runs from the row's start to the PAM: on the + strand with a 3′ PAM (the −
strand with a 5′ one), and on a row that starts on a reference base. On the − strand with a 3′
PAM (the + strand with a 5′ one) the start is the PAM's far end, and the bases read were
len(ALT) − len(REF) bases from the reference PAM beside the protospacer: across a deletion's
junction they included the deleted bases, the PAM of another protospacer one base on or more
(README guide, HGDP + 1000 Genomes: chr22:21056178 (−), CT>C at 21056179, read CTC, GAG on −,
where the reference holds TCT, AGA on −, beside its protospacer: now AGG); on a row starting on
the anchor of an insertion inside the PAM, the anchor and the bases after it (chr3:154367344 (−),
C>CCT, read AGG as created where the reference holds GGG beside the protospacer: now NA). A row
that starts inside an insertion started at its first base's place on the haplotype, past the
anchor, and the bases read were reference bases the target does not hold (chr22:21052639 (−):
TGG, where the reference has CAG beside its protospacer: now NA; with NRR, chr22:44115947 (+),
wholly inside a 37- or a 33-base insertion: now TAA, at 44115944 and 44115946). And every PAM
holding inserted bases was read where the VCF writes the insertion, so a repeat unit inserted
next to a reference PAM read as created, in 2.10.0's first fixes (README 48 of 75 such
alignments genome-wide; NRR chr22:25415044 (−): CAA, the same alignment as its site's reference
row; TRAC guide chr22:49196564 (−) and 50193928 (−): now NA, as 2.8.1 read them). Genome-wide,
2.8.1's tables (2.8.1's code gives back every value in them) with each row's PAM_creation
computed by 2.10.0's code, highest-CFD alignments: README guide, per-sample, 97 go NA → created
(70 rows and 27 alternative alignments: a PAM across a deletion's junction on the − strand 56,
CFD at most 0.083; a PAM holding inserted bases on a row that starts inside the insertion 27, on
a row that starts on a reference base 14) and 4 created → NA (3 rows and 1 alternative
alignment); population-level 70 and 3 rows; TRAC guide (CTCTCAGCTGGTACACGGCA), per-sample, 109
(67 and 42: deletions 29, CFD at most 0.039; rows inside an insertion 58; on a reference base 22)
and 5 (2 and 3).

From 2.5.0 to 2.8.0 each code was instead read **column by column without genotypes** —
the guide's base when it was an allele, else the code's first allele. That reported one
individual's indel together with SNP alleles that only other individuals carry, scored as
if real (on the default HGDP + 1000 Genomes index, a guide's report listed at its fifth
highest CFD, 0.708, a site whose only listed carrier is homozygous reference at the SNP the
displayed sequence uses), understated mismatch counts, and read an R at NRG's last PAM base
as A (a PAM the search never matched; CFD ≈ 0). 2.8.1 fixes this; results of earlier
releases are corrected only by running the search again.

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

- **`--max-total-edits` (default 6).** An independent cap on the combined number of
  edits (mismatches + DNA/RNA bulges) per alignment, pruned inside the TST search. It
  defaults to **6** — matching the default search depth of **4 mismatches + 1 DNA
  bulge + 1 RNA bulge** — regardless of the requested `--mm`/`--bDNA`/`--bRNA`, so a
  permissive request (e.g. `--mm 6 --bDNA 2 --bRNA 2`, budget 10) still drops any
  alignment needing more than 6 combined edits unless raised. `complete-search` prints
  a WARNING when the cap is below the requested budget; set `--max-total-edits` to
  `mm+bDNA+bRNA` (the sum of the requested per-type caps) to keep all requested-depth
  off-targets. It bounds the
  variant-density combinatorial explosion. (`assembly-search` instead defaults it to
  the full `mm+bDNA+bRNA` budget — a resolved diploid genome has no IUPAC-variant
  lattice, so it needs no variant-density cap.)
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

On a **genotyped** panel the observed-haplotype enumerator (§4) already emits every
carried multi-variant haplotype exactly, so the greedy representative is only a
tractability fallback for pathological windows. On a **sites-only** panel there are no
carriers to enumerate, so a dense window emits only the greedy min-mismatch
representative — which can be a strict subset of a genuinely co-located haplotype.
`CRISPRME_LOSSLESS_DENSE` closes this gap for sites-only panels: it additionally emits
the **full co-located variant union** for the window (the maximal PUTATIVE haplotype),
bounded by the carrier-free union rather than the 2ᵏ lattice and gated by the same
mismatch/PAM budget so no over-budget or PAM-invalid row is produced. It is **on by
default** (so a sites-only panel fulfils "don't miss a region") but **scoped to the
sites-only path** — the effect requires `registry_only_mode`, so it is byte-identical
for every genotyped / legacy install (the genotyped path is already lossless via the
observed enumerator). Set `CRISPRME_LOSSLESS_DENSE=0` to opt out. A per-sample union on
the genotyped path is intentionally *not* done — it would risk trans-as-cis phantoms.

Under `--per-sample` a capped window emits no greedy representative at all (since 2.7.10):
the representative unions the carriers of alleles that different individuals carry, so it
can be a haplotype no single one has, and the window is reported by its observed haplotypes,
its reference row and its BED line. That is the observed-haplotype path, which reads the
index's dict-less tiers. **Known limit:** an index without them (one built without
`--samplesID`, such as the one the benchmark `complete-test` builds) goes the classic
dictionary path, which still reports a capped window by its greedy worst-case alignment,
flagged in `High_variant_density_region`, whose listed carriers (the carriers of any allele
it takes) need not carry every variant in it. On the NGG index the README's self-test builds
(`complete-test`'s benchmark mode: `NGG_2_hg38_chr22+hg38_1000G`, chr22 of the benchmark's
phased 1000 Genomes VCF, built without `--samplesID`), the README's guide
(`CTAACAGTTGCTTTTATCAC`) with the `20bp-NGG-SpCas9` PAM, 4 mismatches, 1 DNA and 1 RNA bulge,
searched per sample at the default `CRISPRME_IUPAC_CAP`, gives 5 such rows (8 to 13 variants
each, CFD at most 0.001). They name 102 individuals as carriers 114 times (some on more than
one row), and each time the carrier lacks one of its row's variants: 342 of the run's 56,730
SNP-row carrier listings (a row lists its carriers once for each of its three alignments).
After `build-index-only --samplesID` the same search returns exactly those 5 rows fewer, and
every carrier it lists has all of its row's variants on one phased haplotype. `complete-search` and
the report say so (its Search mode row; and, when no dataset of the search has the tiers, such
a row of two or more SNPs reads `N carriers (>=1 of K vars, poss. chimera)` in its `Observed`
column, as in a population-level analysis, not `N cis carriers`: the report cannot tell which
dataset a row comes from, so beside a dataset with the tiers such a row still reads
`N cis carriers`, and its legend says why); rebuild the index with
`build-index-only --samplesID`, or use a prebuilt one. It is fixed in the next release.

### Context-aware analysis mode: per-sample genotype resolution or population-level analysis

The controls above bound any *single* window, but a **dense panel** (many merged
sources) or a **sites-only aggregate panel** can present so many variant-dense windows
that even the observed-haplotype enumeration of Section 4 becomes intractable — measured
at **49 h+ without completing** on a 4×-density 1000G+HGDP panel. The two analysis modes
are therefore not really a *speed* dial but a **genotype-resolution** dial: per-sample
resolves each individual's observed haplotype (only possible, and only meaningful, when the
panel carries per-sample genotypes), while the population-level analysis reports
worst-possible off-targets (works on any panel). As of **2.8.0** the mode is **context-aware
by default**: on a genotyped+PHASED index CRISPRme+ runs the exact observed-haplotype
enumeration of Section 4 (per-sample) by default, and on a sites-only / unphased / hybrid /
reference-only index it runs the two-pass population-level analysis by default (propagated to
the whole post-analysis via the internal `CRISPRME_FAST_MODE` env var). **`--population-level`**
forces the worst-possible screen (the opt-out from per-sample on a phased index); **`--per-sample`**
forces the observed-haplotype enumeration (the two are mutually exclusive). The 2.5.4 two-pass
population-level machinery is now the opt-out path on a phased index and the default elsewhere.
Instead of enumerating the 2ᵏ IUPAC haplotype lattice per window, the population-level analysis
emits a small fixed set of **worst-possible representatives** per window:

- **Pass 1 — score-free find.** The window's per-position IUPAC allele sets yield a
  **minimum-edit** representative whose edit distance `D` (the additive-per-column argmin)
  **lower-bounds every realizable haplotype**. A window is therefore dropped only when `D`
  already exceeds the requested budget — detection stays **lossless** (no locus is lost),
  while the whole 2ᵏ expansion is skipped.
- **Pass 2 — worst-case score.** Each surviving window emits (i) its **reference**
  alignment, (ii) the **minimum-edit** representative, and (iii) the **maximum-CFD**
  representative (Section 8). This collapses the per-sample lattice to O(1) rows per window
  while preserving the window's worst-case scores; per-sample phasing becomes an
  *annotation* rather than a dependency, and rows are tagged **PUTATIVE** (a synthetic
  worst case, not an observed haplotype).

The min-edit and max-CFD representatives are distinct because CFD is position-weighted (Section
8): the fewest-mismatch haplotype is often **not** the highest-scoring one, so both are emitted
so the reported worst case is never understated. The population-level analysis is validated to
be **lossless for locus detection and non-understating for the worst-case score** against the
exact per-sample enumeration path on a real chr22 1000G-2021+HGDP slice (0 CFD under-reports; it
in fact surfaces *stronger* worst cases at 182 loci that per-sample enumeration misses), and it
emits ~1.9× fewer rows on that 1× slice, a factor that grows with density — turning the
otherwise-intractable 4× panel into a tractable run. **The trade-off is that per-sample
carriers, CONFIRMED cis phasing and exact joint allele frequency are not computed in the
population-level analysis** (rows are PUTATIVE worst cases); the launch-time message, the web
**Analysis mode** control and the report's **Search mode** row all say which mode ran, so it is
never silent. The default therefore follows the panel. On a genotyped+PHASED index (such as the
default HGDP + 1000 Genomes one) it is per-sample, the phased haplotype resolution that
*confirmatory / pre-IND* runs need, and `--population-level` is the fast *screen* for a very
dense, high-edit or genome-wide search. On a sites-only, unphased, hybrid or reference-only
index it is the population-level analysis (the only meaningful mode on sites-only panels), and
`--per-sample` forces resolution on a genotyped one. On sites-only (aggregate) panels there are
no per-sample genotypes to recover — verified on the released **mega** index: `--per-sample`
yields **zero CONFIRMED rows and no cis carriers** in either mode (there is no genotype store to
enumerate), and the `indel_snp_cooc.tsv` / `indel_af.tsv` companions are **byte-identical**
across modes. A sites-only dataset has no genotypes for `--per-sample` to resolve, so as of 2.8.1
its SNP post-analysis is the two-pass population-level one whatever the search's mode, dataset by
dataset, from its registry's `data_type` (for a registry with none, its aggregation: an INFO-AF
registry is sites-only): `--per-sample` on sites-only datasets alone **is** the
population-level analysis, and `complete-search` **warns**; a phased dataset searched beside a
sites-only one (per-sample by default) gives per-sample rows, the sites-only one population-level
rows, and `complete-search` says so. (Before 2.8.1 such a dataset went through the per-sample path,
whose emission joins variants only through the carriers they share, and it has none: it lost
multi-variant worst cases, whole sites among them, with *fewer* rows than the population-level path:
integrated 3,065 vs 3,224 rows, SNP+SNP co-occurrence 4 vs 22 PUTATIVE rows.) The web form
**disables** the per-sample option when a sites-only index is selected.

**Measured population-level vs `--per-sample` behavior (adversarially verified).** Running both
modes on the then-released 1000 Genomes 2021 + HGDP genotyped index for one guide (identical
guide/PAM/thresholds; only `--per-sample` differs; population-level was then the default on every
index) makes the trade concrete and confirms the guarantees hold end-to-end:

- *Detection is lossless at the **window** level.* Every off-target window `--per-sample`
  reports is present in the population-level analysis, which in fact reports **more** loci
  (3,287 vs 3,103) because it emits worst-possible representatives. A small number of windows
  are anchored a few bases apart between the two modes — `--per-sample` anchors each *observed*
  haplotype, population-level its worst-possible representative — so a strict
  coordinate-equality comparison flags them as "per-sample-only" (24 loci here, all 1–10 bp
  shifts of a window population-level did detect). Compared by cluster they match exactly and
  population-level is a strict superset. Losslessness is therefore a **window/cluster** property,
  not byte-exact coordinate identity.
- *The worst-case CFD bound holds exactly.* Across all 3,079 shared loci, the population-level
  CFD is **≥** the per-sample CFD with **zero** violations; it is strictly higher only where it
  assumes worst-possible co-occurrence (e.g. CFD 0.354 vs 0.113 at a 3-SNP window), never lower.
- *`--per-sample` corrects — and **tightens** — the carrier sets.* Population-level lists a
  PUTATIVE **union** of every sample carrying *any* contributing variant in the worst-possible
  window; `--per-sample` prunes to the actually-observed cis haplotype's carriers.
  Counterintuitively `--per-sample` therefore names **fewer** samples (1,447 vs 1,642 distinct
  here), not more — the population-level list deliberately over-includes non-cis carriers (a
  conservative screening bias). The SNP+SNP co-occurrence companion is all-PUTATIVE at
  population level (0 CONFIRMED) and gains **106 CONFIRMED phased-cis rows** with exact joint AF
  and single named carriers under `--per-sample` (e.g. a locus that is a 175-carrier PUTATIVE
  union at population level resolves to a single confirmed carrier under `--per-sample`).
- *On a **genotyped** index population-level does not blank carriers.* The `NA`-carrier behavior
  is specific to sites-only panels; a genotyped panel still surfaces dataset-level carriers in the
  population-level analysis, so the real trade is **PUTATIVE-union vs CONFIRMED-exact-cis**, not
  named-vs-`NA`.
- *For a single guide the two modes take essentially the same wall time* (~13–14 min) — the
  modes differ in genotype resolution, not intrinsic speed; population-level's advantage
  materializes only on the dense/aggregate-panel enumeration wall, not sparse single-guide runs.

The measured numbers above, side by side (same guide/PAM/thresholds; only `--per-sample` differs):

| Metric | Population-level | `--per-sample` |
|---|---|---|
| **Genotyped index (1000G 2021 + HGDP), one guide** | | |
| Off-target loci reported | 3,287 (superset; 3,079 shared) | 3,103 |
| Worst-case CFD on the 3,079 shared loci | ≥ per-sample, **0** violations | baseline |
| Distinct carrier samples named | 1,642 (PUTATIVE union, over-includes non-cis) | 1,447 (exact cis carriers) |
| SNP+SNP co-occurrence CONFIRMED rows | 0 (all PUTATIVE) | **106** (phased cis, exact joint AF + named carrier) |
| Single-guide wall time | ~13–14 min | ~13–14 min |
| **Sites-only mega index, before 2.8.1** | | |
| Integrated rows | **3,224** (more complete screen) | 3,065 (registry-only fallback) |
| SNP+SNP co-occurrence PUTATIVE rows | **22** | 4 |
| CONFIRMED rows / cis carriers | 0 / none (no genotype store) | 0 / none (no genotype store) |

(On the sites-only mega, `--per-sample` cannot resolve carriers and emitted *fewer* rows than the
population-level analysis before 2.8.1. It now runs that analysis, whose rows are the left column's,
and `complete-search` warns; the web form disables per-sample for a sites-only index.)

**Scope of the analysis mode.** The mode toggles only the **SNP** post-analysis (`--per-sample`
re-enables the 2^k IUPAC haplotype lattice / observed-haplotype enumeration that the
population-level analysis collapses). The **indel** post-analysis is single-threaded and
ML-scoring-bound (CRISPR-Bulge), and is **unaffected by the mode** — a dense indel search pays
the full indel cost regardless (parallelizing that path is a follow-up). Correspondingly, the
`indel_snp_cooc.tsv` companion is **byte-identical** in both modes (§8).

---

## 6. Functional annotation of off-targets

Each reported off-target is annotated with its genomic context by intersecting
its coordinates with a 4-column BED (`chrom  start  end  label`), where the label
is suffixed by its source. `resultIntegrator.py` buckets labels by suffix into
dedicated columns; CRISPRme+ ships an updated annotation set (ENCODE **SCREEN
v4**) and adds cancer-gene columns — **IntOGen** (CC0) cancer-driver genes on by
default and a licence-gated **COSMIC** Cancer Gene Census column:

- **GENCODE** — gene-model context (`exon`, `CDS`, `UTR`, `transcript`,
  `start_codon`/`stop_codon`; `intergenic` otherwise), plus nearest-gene name
  and distance. The `Gene_region` column records where the off-target sits
  relative to the nearest protein-coding gene — `CDS`, `5'UTR`, `3'UTR`,
  `exon`, `intron`, or `intergenic` — coherent with `Gene_distance_kb`
  (`intergenic` iff the distance is > 0). `Gene_region` is present in the
  curated report table, all download TSVs, and the raw `integrated_results.tsv`.
- **DHS** — DNase I hypersensitive (open-chromatin) sites, labeled by tissue /
  organ system.
- **ENCODE SCREEN v4 cCREs** — candidate cis-regulatory elements: promoter-like
  (`PLS`), proximal / distal enhancer-like (`pELS` / `dELS`), and
  chromatin-accessible / TF classes (`CA-CTCF`, `CA-H3K4me3`, `CA-TF`, `CA`,
  `TF`).
- **IntOGen (cancer driver genes)** — whether the off-target falls in a gene
  reported as a cancer **driver** by IntOGen (Integrative OncoGenomics), a
  compendium of computationally-identified drivers across tumour cohorts. This is
  the **default** cancer-gene flag (`Annotation_INTOGEN`): IntOGen's current release
  is **CC0** (public domain), so it is licence-free for any use — the open complement
  to the licence-gated COSMIC.
- **COSMIC (Cancer Gene Census)** — whether the off-target falls in a curated
  cancer gene, tagged by confidence **tier** (Tier 1: extensive curated causal
  evidence; Tier 2: strong but less-curated) and documented **role**
  (`oncogene`, `TSG`/tumor-suppressor, `fusion`). This flag is particularly
  relevant for therapeutic and pre-IND assessment, where an off-target in a
  known cancer gene warrants scrutiny. Unlike IntOGen it is **licence-gated** (below).

**Cancer-gene annotations: IntOGen by default, COSMIC by licence.** CRISPRme+ ships two
cancer-gene flags with different licences and defaults:

| Flag | Source / licence | Default | Report column |
|---|---|---|---|
| **IntOGen** | Integrative OncoGenomics driver compendium; **CC0** (public domain), legally clean for academic and commercial use | **default-ON** (`intogen_drivers.hg38.bed.gz`, built by `seq_script/build_intogen_annotation.py` mapping the CC0 driver symbols to GENCODE gene spans) | `Annotation_INTOGEN` |
| **COSMIC** | Cancer Gene Census (Genome Research Ltd / Wellcome Sanger); free for academic / non-commercial research, **commercial use requires a licence** ([terms](https://www.cosmickb.org/terms/)) | **default-stripped** — the built-in bundle bakes COSMIC into a single BED, so CRISPRme+ strips its rows from the active annotation and no values appear in the output | `Annotation_COSMIC` |

**Opt-in mechanic.** The user enables COSMIC **once** — a Settings checkbox on the web
(with a link to the terms) or `crisprme.py cosmic-license enable` on the CLI (both persist to
`Annotations/.cosmic_license.json`; the CLI `enable` requires an explicit confirmation).
Toggling it changes the active-annotation cache signature, so the next search rebuilds the
annotation with or without COSMIC accordingly.

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
report of reference vs. variant off-target scores (CFD and CRISPR-Bulge, under
multiple rankings); a reference-vs-population origin breakdown; a **recommended
validation panel** (a hybrid worst-case top-N shortlist selected by combining
CFD, CRISPR-Bulge, edit-distance, and pure-bulge floors, with the selection logic
stated explicitly); per-threshold **downloads** sharing one curated, spreadsheet-ready
column schema; a scrollable **top-1000** table with the functional annotations;
and the annotation legend of Section 6. Allele frequencies can be omitted
(`--no-maf`) for runs where they are not yet finalized, so the site set and
scores can be shared without misleading frequency values.

The report and all of its tables are built from the **best alignment per locus**
(`bestMerge`). The exhaustive **alternative-alignments** file
(`..._all_results_with_alternative_alignments.tsv`, the non-best alignments at each
locus) is a separate, optional artifact whose size grows combinatorially with the
edit budget; it is emitted by mode — **off** for a population-level analysis, **on**
for a per-sample one (the default on a phased genotyped index, or forced with
`--per-sample`) — with `--alt-alignments` / `--no-alt-alignments` to force either way.
(`assembly-search` reconciles its two haplotypes from each per-haplotype
`integrated_results.tsv` alone, so it does not require this file.)
Suppressing it does not change the report or any reported number.

A guide with **more than one perfect genomic match** (0 mismatch, 0 bulge) has no
a-priori on-target — each is an equally-efficient candidate cut site, and a
perfect-match *off*-target is the highest-risk class. Every 0-mismatch site is
therefore forced to the top of the validation panel (never truncated), flagged in
a `Perfect_match` column, and called out in a prominent warning banner (red when
there are several, listing the sites; amber "presumed on-target" when there is
exactly one) — in both the report and the interactive web results page.

An **`Observed`** column records whether each off-target is supported by at least
one real individual:

| Value | Meaning | When shown |
|---|---|---|
| `reference` | present in the reference genome | always → carried by essentially every individual |
| `N carrier(s)` | named carriers | genotyped panel only |
| `observed` | a single AF>0 variant, carried by ≥1 individual even where the carriers can't be named | any panel |
| `putative` | a multi-variant combination whose cis co-occurrence can't be confirmed | any panel |

On a genotyped index the carrier roster is
populated in both analysis modes; under `--per-sample` it is pruned to the **exact
cis carriers** (but for a dense window's greedy alignment on an index without the
dict-less tiers, §5), so `putative` cleanly marks reconstructions that no one actually
carries — separating genuine haplotypes from worst-case combinations. In a
population-level analysis a multi-variant SNP roster is the union of each variant's
carriers and reads `N carriers (>=1 of K vars, poss. chimera)`; an indel (or MNV) row's
roster is resolved per carrier in both modes (§4), so it reads `N carrier(s)`: CONFIRMED
where each carrier is placed on its phased haplotype, PUTATIVE where the phase is unknown
(an unphased call, an MNV). On an install without SNP genotypes (`download --what index
--no-genotypes`, §4), a row an indel and a SNP make together lists every carrier of the
indel, with or without the SNP; its phase is in `indel_snp_cooc.tsv`.

The
recommended validation panel uses this to **prioritize observed sites** among
equally-severe candidates: `Observed` enters the ordering only as the first
tie-break *after* worst-case severity, so a strictly-worse (higher-severity)
reconstructed site is never displaced and worst-case coverage is preserved, while
genuine real-carrier haplotypes of equal severity float to the top of the
shortlist. On a sites-only panel (e.g. the mega index) there is no per-sample
roster, so the `N carrier(s)` value never appears and the panel tie-break has no
effect (every site is `observed` or `putative`); the column itself is still shown.
It is dropped only for a reference-only run with no variant off-targets.

The report ZIP places `report.html` at the top level with all data files under a
`data/` subfolder, and its *Variants included* line states the genotyped panel
size and the number of SNPs and indels searched. These database counts are read
from a build-time `Dictionaries/registry_<vcf>/variant_count.json` sidecar
(`n_records` SNPs + `n_indels`) written when the dict-less index is built, so the
report needs no VCF access at report time (e.g. the earlier 1000 Genomes + HGDP
merge recorded 106,664,924 SNPs). Reported SNP allele frequencies are AC/AN over the
genotyped panel (on an index without the dict-less tiers, the `AF` of the classic per-sample
dictionaries, §3, as the report of a search there says), and an indel's is the `AF` of the
index's indel log; `add-variants` copied both from the VCF's `INFO/AF` (§3), and §8 gives
when they hold none. A variant whose allele frequency is 0 or missing (e.g. an allele that
no individual of the panel carries) has no frequency: `integrated_results.tsv` lists it as
`NA`, and the report leaves it out of a row's MAF, the lowest frequency of the row's
variants, which is blank when no variant of the row has one. No floor is applied: a
frequency below 1×10⁻⁵ is reported as it is, and the plots, whose log scale is the site's
rank, use a frequency only to size the points. For an **aggregate (sites-only) panel** (e.g.
the all-source mega index, §2), only allele frequencies exist — there are no per-individual
genotypes — so no carrier or homozygote count is estimated: the registry holds 0 for both
(§2), which the population-summary companion (`*.population_summary.tsv`, kept with
`--keep-intermediates`) gives as **NA**, and a per-dataset **`indel_af.tsv`** companion
carries indel allele frequencies by source alongside the SNP+indel co-occurrence table.

## 8. Off-target scoring, assumptions and limitations

**Scoring models.** Candidate off-targets are scored with two independent,
previously published models, used **unchanged** from their original definitions:
**CFD** (Cutting Frequency Determination; Doench *et al.*, *Nat. Biotechnol.*
2016), the **primary** score, and **CRISPR-Bulge** (Yaish & Orenstein, *Nucleic
Acids Res.* 2024, doi:10.1093/nar/gkae428), a deep-learning score built for gapped
off-targets: a TensorFlow GRU-embedding five-model ensemble (MIT-licensed) that
was **trained explicitly on bulge-containing off-targets**. On the authors'
Refined-TrueOT benchmark it reaches full-set AUPR **0.498** and, on the case that
matters most here, **bulge-only AUPR 0.245**, while scoring ≈ 14.5k
off-targets/s at 16 CPU threads (a ~7.8 min genome-wide pass). Both scores return
a value in [0, 1]; higher means more likely to be cut. They are reported
side-by-side because they can disagree, and the recommended validation panel is
deliberately model-agnostic (a site worst by *any* metric is included), so no
single model gates the shortlist. CRISPRme+ does not re-train or modify either
model; CFD is applied as its published lookup tables (CRISPRme+'s own
implementation, `PostProcess/cfd_score.py`, Doench et al. 2016) and CRISPR-Bulge
as its published, pinned weights, to the aligned protospacer+PAM of each candidate.

**Domain-of-validity caveat (important).** Each score has a domain where it is
trustworthy and a boundary past which CRISPRme+ refuses to fabricate a value:

- **CFD on gaps is a relative risk indicator, not a calibrated indel probability.**
  CFD was trained on **single-nucleotide mismatches**, but it does *not* ignore gaps:
  the CFD gap table (`cfd_mismatch_scores.json`) carries **152 position-specific gap entries**
  (DNA-bulge `r-:dX` and RNA-bulge `rX:d-` at protospacer positions 2–20), so a bulge is penalized by a
  per-position factor and **two bulges multiply** — a PAM-distal bulge is near-neutral while a
  seed-proximal one is strongly penalizing. **As of v2.8.0 these factors are applied at each
  bulge's true spacer position** (the extra genomic base of a DNA bulge is removed so the 20
  guide bases keep their Cas9 positions, then the gap factor is applied; the two undefined
  terminal columns are extrapolated — lenient at the PAM-distal 5' tip, abolished at the
  PAM-proximal junction). Before v2.8.0 a DNA-bulge realignment discarded the gap, so DNA
  bulges — and any bulge at a terminal column — were scored as if absent (a PAM-proximal bulge
  scored like a perfect match); mismatch-only and in-range RNA-bulge scores are unchanged by the
  fix. Even corrected, a CFD value on a gapped alignment is only a **relative** risk indicator
  (held-out bulge-only AP ≈ 0.01) — **CRISPR-Bulge is the calibrated bulge ranker**; see
  `docs/offtarget_risk_summary_METHODS.md` §3.
- **CRISPR-Bulge is the model built for gaps — for exactly one bulge.** It natively scores
  a mismatched protospacer carrying **exactly one 1-bp DNA or RNA bulge**; its trained
  weights bake in a fixed input width.
- **≥ 2 bulges are out of domain and nulled to N/A.** An alignment carrying **≥ 2 bulges is
  outside CRISPR-Bulge's domain** and is deliberately reported as N/A (sentinel `-1`, never
  sorted or thresholded) rather than fabricated from an out-of-domain extrapolation. The
  model's training corpus of ~4.6M measured off-targets contains **zero ≥ 2-bulge examples**,
  so any 2-bulge output would be pure extrapolation — the honest choice is N/A plus the
  model-agnostic safety nets.
- **Nulled ≥ 2-bulge sites are *not dropped*.** They remain scored by CFD's position-specific
  bulge penalty and are covered by two validation-panel safety nets. First, an unconditional
  **pure-bulge** rule: a 0-mismatch site with ≥ 1 bulge is always hard-included. Second, a
  ≥ 2-bulge site that is high on **CFD alone** (CFD ≥ 0.5, where CFD position-penalizes each
  bulge but is an uncalibrated extrapolation outside the single-bulge domain, and CRISPR-Bulge is
  N/A) is **kept in the recommended panel but *demoted*** — ranked
  after every reliable hard-include (sites close by sequence, or high CFD at ≤ 1 bulge in its
  calibrated domain, or a clean bulge). We demote rather than drop because a ≥ 2-bulge
  off-target can still be real — the benchmark simply contains no ≥ 2-bulge examples, so its
  risk is *unknown*, not zero — and only the user can decide what to assay. Demoting keeps the
  trustworthy candidates first without ever silently removing an unknown-risk site, and the
  count of such demoted sites is surfaced in the panel note even on a CFD-only run.

The N/A sentinel is uniform across the report: a `-` means *not scored*, not *no risk* — it is
never sorted or thresholded as if it were a low score. The threshold tiers in the report are
**model-relative** (CRISPR-Bulge's cut points differ from CFD's because the two scores are on
different scales). The CRISPR-Bulge **sensitive** tier — the operating point the recommended
panel keys on — is the **95%-recall threshold** on the 40-guide held-out benchmark
(**CB ≥ 1.9×10⁻⁴**; point estimate, guide-bootstrap 95% CI [5.8×10⁻⁶, 3.1×10⁻³]; selectable to
90% / 99% via `CRISPRME_CB_RECALL_TARGET`), not an ad-hoc `> 0`: CRISPR-Bulge scores almost every
site non-zero, so `> 0` would list the whole genome. CRISPR-Bulge scores are written at full
precision (`{:.6g}`) so this sub-5×10⁻⁴ operating point is representable. See
`docs/offtarget_risk_summary_METHODS.md` §3.5 for the recall↔validation-size trade (162 off-target
validation loci at 95% recall for the ACTGA/NRG worked example).

**Top-of-report off-target risk summary (replaces the aggregated specificity score).**
The legacy single `100/(100+Σ)` specificity score is replaced by a per-scorer risk card. For
each scorer, shown side by side and never blended, the card reports the count of ≤1-gap
off-target sites in each score band and a **calibrated expected number of true off-targets**
with a coverage-valid 90% prediction interval. The calibration is an isotonic
score-to-probability map fit leave-one-guide-out on an independent benchmark of 359
experimentally confirmed off-targets across 46 guides (four published external-evaluation sets
distributed with the CRISPR-Bulge repository: Tsai 2015, Refined TrueOT, Chen 2017, Listgarten
2018; Zenodo 10.5281/zenodo.10902445). Held-out AUPR on that benchmark is **0.168 (CFD)** and
**0.613 (CRISPR-Bulge)**, with guide-cluster-bootstrap 90% CIs of 0.13 to 0.24 and 0.56 to
0.70. The reported range is a two-level guide-cluster bootstrap plus a negative-binomial
between-guide overdispersion term, validated per-held-out-guide to about 90% coverage. Because
both scorers are sequence-based, the calibration applies to any ≤1-gap alignment, reference or
variant allele; ≥2-gap sites and the per-sample carrier-incidence dimension are out of scope.
The card also flags, qualitatively, that a deployed guide's own off-target propensity is
unknown and can shift the absolute count substantially in either direction — usually a larger
uncertainty than the statistical range, and deliberately not given a specific multiplier. The card restricts CFD counts to the ≤1-gap calibration domain, whereas
the recommended validation panel counts candidates over all gaps; a bridge note in the panel
states this so the two are not read as contradictory. Full method, benchmark provenance,
precision-recall curves, and the independent statistician, clinician, and genome-editor
sign-off are in `docs/offtarget_risk_summary_METHODS.md`.

**Compute backend and modular environment.** CFD is pure-Python and always runs
in-process on the CPU. CRISPR-Bulge (TensorFlow) runs in a **dedicated per-scorer
conda environment** (`cbulge`) driven by a persistent worker subprocess, so its
heavy, version-pinned dependency stack is fully isolated from the main CRISPRme+
environment (whose legacy scikit-learn / numpy pins have been dropped and
modernized accordingly). The scorer honors `--compute-backend {cpu,cuda,metal}`
(`gpu` and `auto` are accepted aliases: `auto` picks cuda→metal→cpu): **cpu** is
the default and the numerical oracle; **cuda** (conda-forge `tensorflow=2.13`
CUDA build) is validated bit-close to CPU on NVIDIA A100 (~1e-5 agreement) for
Linux/Docker/GPU runs; **metal** (Apple Silicon) is guarded by a **load-time
numerical self-test** — `tensorflow-metal` silently miscomputes this GRU kernel
(collapsing every prediction to 1.0), so on detecting the miscompute the scorer
transparently falls back to a CPU device context, producing **bit-identical**
results. With the env or its weights missing or broken, `complete-search` stops
before the search, naming the fix, unless `CRISPRME_SKIP_SCORER_ENV=1` is set, which
runs it **CFD-only**; `targets-scoring` and searches started from the web interface do
not check first, and complete **CFD-only**. A CFD-only run has the CRISPR-Bulge score at
`-1` and a banner saying so in its log. The released **amd64** Docker image bundles the CUDA
TensorFlow build, so `docker run --gpus all ... --compute-backend cuda` runs
CRISPR-Bulge on an NVIDIA GPU with transparent CPU fallback; the **arm64** image
stays CPU.

The scorer environment is **architecture-aware** and yields **numerically
identical** CRISPR-Bulge scores across CPU architectures: it installs the same
TensorFlow 2.13.x per-arch (2.13.1 on x86-64/aarch64; 2.13.0 via tensorflow-macos
on Apple Silicon) — the conda-forge `tensorflow-cpu`/CUDA builds on
x86-64, `tensorflow-macos` on Apple Silicon, and the PyPI `tensorflow==2.13.1`
aarch64 wheel on Linux ARM (conda-forge ships `tensorflow-cpu=2.13` only for
x86-64). A cross-architecture check on a 254-off-target panel (150 with bulges)
found x86-64 vs aarch64 scores agree to a maximum absolute difference of
3.6 × 10⁻⁷ with zero threshold-band changes, so the multi-architecture
(amd64 + arm64) Docker image scores identically on both. On networks that block
`conda.anaconda.org`, the environment is built from a mirror (e.g. prefix.dev) via
`CONDA_CHANNEL_BASE` / `CRISPRME_CONDA_CHANNEL_BASE` (the Docker build exposes a
`CONDA_CHANNEL_BASE` build argument).

**Worst-case scoring in the population-level analysis.** In the population-level
analysis (§5), each window is represented by worst-possible rows rather than every haplotype,
so the *scores* attached to those rows are defined as worst cases over the window's allele
combinations. The two scores earn different trust:

- **CFD is the exact worst case — the safe actionable gate.** CFD factorizes as a product of
  per-position maxima times a **joint two-base PAM factor**, so the maximum over all
  combinations is found by a per-position argmax plus a bounded brute-force over the PAM region
  (the joint factor is why a naïve per-column greedy is insufficient). This exact maximizer is
  validated bit-for-bit against an independent factorized oracle and against the exact
  per-sample enumeration path: **zero CFD under-reports** on a real chr22 1000G+HGDP slice
  (and, on the legacy dict / aggregate-panel path, catching cases where the fewest-mismatch
  allele scores materially *lower* CFD than another carried allele, up to a threshold-crossing
  0.14). A whole-index comparison reproduces this end-to-end — across all shared loci the
  population-level CFD ≥ the per-sample CFD with **zero** violations, strictly higher only where
  it assumes worst-possible co-occurrence (the per-index numbers are in the measured comparison
  in §5).
- **The CRISPR-Bulge ML score is best-effort — a screen, not an action gate.** It is a
  non-factorizable neural ensemble, so its worst case is taken as the maximum over the emitted
  representatives rather than an exhaustive per-haplotype search. This is tight exactly where
  decisions are made: high-ML-score off-targets are low-edit, and the min-edit + max-CFD
  representatives already span the low-edit shell, so strong sites are surfaced at full or
  greater strength and any under-report is confined to the weak tail.

So by default **CFD is the safe actionable gate** (exact worst case) and the **ML score is a
screen**. A **guaranteed per-haplotype ML worst case** is available by running `--per-sample`;
this is the screening-vs-confirmatory two-tier split of Section 5.

**SNP+indel co-occurrence is unaffected by the analysis mode.** The default analysis collapses
only the *SNP* worst-possible representative emission in `integrated_results.tsv` (Section 5);
the SNP+indel co-occurrence companion (`indel_snp_cooc.tsv`) is produced by the indel
post-analysis' cis phasing pass over the genotype tiers, which the mode does not touch. Measured
on the complete genome-wide matrix (2021 panel, same guide, default vs `--per-sample`): the two
`indel_snp_cooc.tsv` files are **byte-identical** (same MD5, 2,729 rows, 843 CONFIRMED / 1,886 PUTATIVE, full
per-sample `cis_samples` and joint-AF in both). So per-sample cis attribution — which individual
carries the indel and SNP together — is preserved identically in both analysis modes. Measured
again with 2.8.1, whose companion holds only the SNPs the indel post-analysis's rows require, on
the HGDP + 1000 Genomes index: the two modes give the same lines for each of two guides (339 and
363 distinct lines, all CONFIRMED).

**Two co-occurrence companions.** Alongside `indel_snp_cooc.tsv` (SNP+indel),
`snp_snp_cooc.tsv` reports **SNP+SNP** co-occurrences — off-targets that require ≥2
nearby SNP alt alleles together. Both use the same three-rung confidence model (§4):
CONFIRMED (phased cis, exact carriers + joint AF), PUTATIVE co-carrier (genotyped
unphased), and PUTATIVE estimated (sites-only, min marginal AF as a conservative
upper bound, no carriers). On a sites-only panel the SNP+indel companion likewise
falls back to a PUTATIVE min-AF row rather than emitting nothing, so a co-occurrence
is never dropped merely for lack of genotypes. Both companions are bundled into the
report ZIP with a confirmed-count summary.

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
chromatin accessibility as a cutting determinant. **SNP+indel co-occurrence** is
searched **by default** as of 2.5.0: the build overlays SNP IUPAC codes onto the
fake-indel genome and compiles a phased indel genotype tier, and post-analysis reports
off-targets that require **both** a nearby SNP **and** an indel in the same protospacer,
each resolved on its carriers' own haplotypes (§4; since 2.8.1) and tagged CONFIRMED-cis
(phased) / PUTATIVE (unphased) with per-sample carriers and joint allele frequency. (The
pre-2.5.0 behavior — two independent passes, SNPs on the IUPAC-enriched genome and indels
on a plain-reference fake-indel genome — remains available by disabling the integration;
classic dict builds are byte-identical, and searches run the two passes with the fixes of 2.8.1
that apply to every search, so their results can differ from an earlier release's; three of
those fixes change indel rows (§4): a target that holds no base of its indel gives no row, an
indel row with an RNA bulge starts at its first aligned base with its REF gapped as its target,
and of the alignments of one site that name the same variants with the same carriers the merge
keeps the best, not the first it reads.) One **residual**: the indel search materializes
**one indel per fake contig**, so an off-target requiring **≥ 2 co-occurring cis indels
within a single protospacer** is not generated as a candidate — a pre-existing
single-indel-search property, independent of the analysis mode; since 2.8.1 an individual
whose haplotype carries such a second indel inside a target is not listed for that target
when the second indel changes the target's bases, read on the haplotype's own bases (a
repeat unit gained or lost at the target's edge leaves them, and the individual stays
listed: §4), or when those bases cannot be read (an unphased or missing SNP call, an
unphased indel that does not itself change them, or two of its variants that overlap: §4);
a carrier whose second indel inside the target is unphased cannot be tied to a haplotype
and stays listed, PUTATIVE, even when another, phased second indel on its indel haplotype
changes the target, nor when, as one complex allele split into overlapping records,
an indel overlaps it at another position. This residual is **low-frequency**: raw
multi-indel cis co-occurrence is dominated by STR/VNTR repeats (which off-target analysis should
soft-mask), falling to **~1–2%** of indel loci after
repeat-masking and deduplication, and the **genuinely-missed** off-targets are **~0.1–0.2%**
of indel off-targets — all at the edit-budget ceiling (the weakest, ≈0-CFD tier). Windows
carrying ≥ 2 cis indels **can** be flag-all'd for conservative (lossless over-reporting)
treatment — a design option the min-edit primitive supports (proven lossless in
`test_twopass_lynchpin_counterexamples`), **not yet wired into the production indel search**.
Two SNP records at one position (a split multi-allelic site) keep only the last record's
alt in the enriched genome and its dictionaries (a CRISPRitz `add-variants` property), so
a carrier of an earlier record's alt is read as reference there, in the SNP and the indel
post-analysis alike, however common that alt: on the HGDP + 1000 Genomes panel
chr8:4174067 T>C (8,170 copies) is read as reference because T>G (1 copy) is the later
record. That is 1.95 % of the panel's SNP records, whose effect §2 gives as measured; it
will be fixed, with rebuilt indexes, in the next release. `add-variants` also looks for the
`AF` key in the VCF's first record only, and reads every record's `AF` at that place in
`INFO`: a record that has it at another place gets none, in the index's indel log and in its
classic per-sample dictionaries. Its indel then has no frequency on any index (§3), nor its
SNP on an index without the dict-less tiers (the registry, which gives a SNP's frequency with
them, counts AC/AN from the genotypes, or on the mega reads each source's `AF_<dataset>` by
its key), and the report's MAF leaves it out. A VCF merged with
`seq_script/merge_panels/merge_vcf_panels.sh` has such records, since a site of one source has
fewer keys before `AF` than a site of two: the retired `NRG_3_hg38+hg38_1000G2021_HGDP` index
(1000 Genomes 2021 + HGDP, built from such a merge, retired in 2.7.8) has no `AF` for
14,567,429 of its 15,583,931 indels (224,981 of the 232,384 of chr22). Of the shipped indexes only HPRC's chrY is affected,
whose VCF has `AF` at five places in `INFO`, after one to five other keys: its log has none for
91,950 of its 102,665 indels. The HGDP1kGP index and HPRC's other chromosomes have an `AF` for
each of their 17,192,381 and 12,000,599 indels, as the self-test's 1000 Genomes index (chr22)
has for its 64,528 indels and 993,880 SNPs; the mega's VCF has no `INFO/AF` at all (§2), so
none of its 13,853,761 indels has one, and their per-dataset frequencies are in the
`indel_af.tsv` companion. 2.8.0 behaves the same; it will be fixed in the next CRISPRitz
release (HPRC's chrY log needs rebuilding with it). The phased indel genotype store
leaves out equal-length multi-base records (MNVs), so their carriers are resolved with the
phase unknown, PUTATIVE (§4): the HGDP + 1000 Genomes index has none, HPRC's chr22 has 9,079
among its 195,813 indel-search records. Since 2.10.0 a per-sample SNP row does not list a
carrier whose haplotype holding its SNPs also holds an indel that changes the target's bases
(§4); up to 2.8.1 it did, and so it does with `CRISPRME_INDEL_SNP=0`. What remains: an MNV
inside the target, which the phased indel store does not hold (the HGDP + 1000 Genomes index
has none); a carrier whose indel call is unphased, which stays listed, PUTATIVE; a dense
window reported by its greedy worst case on an index without the dict-less tiers (§5), whose
carriers need not carry every variant in it; and the population-level analysis, whose rows
are worst-case composites.
On a row of several variants, the `Variant_rsID` and `Variant_MAF` lists of the results
tables are not in the order of its `Variant_info_genome`: the merge of a cluster's rows
(`remove_contiguous_samples.py`) sorts the three lists each on its own, as text, and keeps
each value once, so the n-th rsID or MAF need not be the n-th variant's, and two variants
of one frequency show it once. On ACTGA (NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per
sample, genome-wide, HGDP + 1000 Genomes), 953 of 1,967 rows of two or more SNPs and 118 of
293 rows of an indel and one or more SNPs list their rsIDs in another order than their
variants (in the highest-CFD columns). The report's MAF, the lowest frequency of the list,
does not depend on it; the rsIDs its `Variant` column gives after the variants, and the rsIDs
and AFs of `observed_haplotypes.tsv` and of a site's haplotype table in the report, are in
the list's order. 2.8.0 behaves the same, and it will be fixed in the next release.
The `Max_total_edits` value is a
**search cap on the variant-collapsed (IUPAC) genome**; individual variant-expanded SNP
alignments may exceed it (the report surfaces the observed maximum), while an indel row,
recounted on its carriers' own sequence, is dropped past it as past `--mm` (§4). No floor is
applied to a MAF: a variant whose allele frequency is 0 or missing has none, which a
multi-variant row's minimum leaves out, and a reported MAF below `1e-05` is the frequency
itself (§3, §7). These are computational predictions and are **not a substitute for
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
  including the single-`--max-total-edits` web-mode cases at 3 and at 5 total edits
  (all six on chr22); `validate-test` exits non-zero on any mismatch, and the whole
  round-trip runs in CI.
- Separately, the high-variant-density greedy decomposition (§5) was validated as
  provably exact against brute-force argmin on **4,000/4,000 random cases** (both
  strands, with/without bulges), with PAM-creating-variant cases reproducing full
  enumeration (variant attribution identical).
- The **two-pass population-level analysis** (§5) was validated against the exact per-sample
  enumeration path on a real chr22 1000G-2021+HGDP slice: **lossless locus detection** and a
  **non-understating worst-case bound** (0 CFD under-reports; 182 loci where it surfaces a
  *stronger* worst case). Its exact worst-case-CFD maximizer is additionally cross-checked on
  4,000 random windows against an independent factorized CFD oracle (agreement to the raw
  double, including the joint-PAM case), and on the legacy dict / aggregate-panel path against a
  real CFD-scored multiallelic fixture; the ML score's (CRISPR-Bulge) best-effort bound is
  described under *Worst-case scoring in the population-level analysis* above.

This establishes that the engine **does not miss** off-targets relative to
exhaustive search. It does **not** validate the *scoring* models' predictive
accuracy against experimental cleavage — a regulatory-grade package should add a
retrospective comparison of CFD/CRISPR-Bulge ranking to experimental off-target
assays (e.g. GUIDE-seq / CIRCLE-seq / targeted amplicon or rhAMP-Seq).

---

## 9. PAM-less indexes and dynamic PAM filtering

A CRISPRitz index is built for one PAM. The PAM is matched against the genome at build
time, only the surviving sites enter the ternary search tree, and the PAM region itself is
stored as a flank beside each guide. One index therefore answers one PAM: searching NGG
and NRG over the same genome means building and storing two indexes.

A **PAM-less** index is built with an all-`N` PAM (`PAMs/20bp-NNN-NO-PAM.txt`), which
constrains nothing, so every guide-length window enters the tree and the stored flank is
the raw genomic sequence. `--prune-pam` then supplies the real PAM at search time and
discards a leaf whose stored flank cannot match it — at the moment the leaf is reached,
before it is written out.

**Validity.** A search PAM of length *L* is served by a stored flank of width *K* when
*L ≤ K* and the flank is on the PAM's side of the guide. The flank lies where the PAM file
the index was built from puts the PAM (the 3′ side for `PAMs/20bp-NNN-NO-PAM.txt` and every
published PAM-less index), and the engine reads it on the side the search's PAM file names,
so a 5′ search of a 3′-built index finds nothing at all. A dedicated index is named after its
PAM by PAM files of both sides alike (`N…NTTTV 4` and the 5′ Cas12a file both build
`TTTV_<N>_<genome>`), and read on the wrong side it reports sites that are not there. CRISPRme
therefore takes any index only for a PAM on the side recorded in its `.pam_build` (written by
`build-index-only`, and by a search that builds an index). A dedicated index without one is used
as before. A PAM-less index without one is taken for a 3′ build at the published widths, 3 and 8
(the K = 8 archive has no `.pam_build`), or under the search PAM's own name, and at no other
width: a search before 2.8.1 built one of any width from a 5′ PAM file and recorded nothing. A
search builds its index under a name an index for the other side may already have, and stops
rather than search that one; `build-index-only` does not record a PAM file over an index
recorded for the other side, nor a 3′ one over an unrecorded PAM-less index of another width.
The search PAM is
widened to *K* with `N`s on the distal side (the end for a 3′ PAM, the front for a 5′ PAM),
and a leaf is rejected when, at any position the real PAM covers, the stored and requested
IUPAC codes share no nucleotide. Pad positions are not tested at all: they state no
requirement, and a stored genomic `N` intersects nothing, so testing one would reject a
target that a dedicated index of the real width keeps. At *L = K*, on the same side, the
pad is empty and the geometry is exactly a dedicated index's.

**Equivalence.** Pruning changes only which leaves are written out, never how a guide is
aligned to the genome, so it cannot move an alignment or a score. Measured genome-wide
over HGDP + 1000 Genomes (4,091 samples, guide `ACTGAAATCTGTAAGCAGGC`, NRG, 4 mismatches +
1 DNA bulge + 1 RNA bulge, per-sample): a 3-nt-flank PAM-less index with `--prune-pam`, an
8-nt-flank PAM-less index with `--prune-pam`, and a dedicated NRG index produce
**byte-identical** reports: the integrated results are identical, 182,369 rows each (MD5
`c76f3490` once the licence-gated COSMIC column, empty on the dedicated arm's fresh install, is
set aside; the two PAM-less arms are identical with it, MD5 `c4586f7a`). The same holds for
NRR: 386,012 rows each (MD5 `f810dd76`) over a dedicated index, both PAM-less indexes with
`--prune-pam`, and the 3-nt one without it. That holds at *L = K*, and at *K > L* for a guide
with no site at a contig's very end or start, as this one has none. The build scans a contig
for a PAM only where a whole window of *K* bases fits, so with *K > L* a + strand PAM with 1 to
*K − L* bases after it to the contig's end, and a − strand PAM in the contig's first *K − L*
bases, are not indexed, and a − strand PAM exactly *K − L* bases from the start is reported on
the + strand at another coordinate. On hg38 at *K* = 8 that is 145 + 101 NGG and 308 + 241 NRG
sites (about 7.5 × 10⁻⁷ of all), none on a primary chromosome for those PAMs; in two
whole-genome searches the 8-nt flank lost one row each against a dedicated index (27,215 of
27,216 rows for a guide picked for its perfect match on `chr11_JH159136v1_alt`, and 216,977 of
216,978 for VEGFA site 3 at 5 mismatches). A search that reads a flank wider than its PAM says
so on stdout and in the run's `Params.txt`. N runs inside a contig are not affected by the
width.

**Genomic `N`.** This equivalence needs engine support, not just a filter. A genomic `N`
is stored as a code that intersects nothing, so it fails a PAM match at every position; a
PAM-less index of width *K* stores a wider flank than a dedicated index of a shorter PAM
and so is more exposed to one. Three separate gates acted on that — the build-time PAM
scan, the discard of any window containing an `N`, and the prune test itself — and the
effect is concentrated on the synthetic indel genome, which is `N`-separated roughly every
56 bases: a width-31 window admitted 34,690,315 start positions on one contig against a
width-26 window's 41,454,105, so **16.3%** were invisible to a wider index. Real assembly
`N` sits in megabase gaps where nothing matches, which is why reference and SNP searches
were never affected. CRISPRitz now keeps a target whose stored PAM/flank region holds an
`N` and lets the search PAM decide, exactly as a dedicated index does.

An index built before CRISPRitz 2.9.1 (CRISPRme 2.8.0 and earlier) lacks every target whose
flank holds an `N`, and searched with `--prune-pam` for a PAM narrower than its flank it misses
the sites whose flank holds an `N` past the PAM (3 of 36 rows on a synthetic genome, 2.6 % of
the indel rows on chr22). Its `.pam_build` has the same format, but its trees tell: one built
since stores an `N` (a 0 nibble) in the flank of a target beside one, and no earlier build does.
A search that reads such an index wider than its PAM warns on stdout (it reads the smallest
tree of a contig whose sequence holds an `N`, and of the INDELS index one of several indel
windows), and goes on.

**Contig ends and N runs.** These limits hold for every index, dedicated or PAM-less, at any
width. An index of budget *N* holds a window of the guide, the PAM and *N* more bases, and the
build drops a window that holds an `N` or runs past the contig's start, so it cannot report a
site whose guide lies within *N − 1* bases of a contig start, a contig end or an `N` run, even a
perfect match; a search without bulges, which scans the FASTA, reports it. A variant search
without bulges on an install without its whole enriched genome (a download) searches the variant
index instead, so its SNP variant rows have this limit too: of such a site it reports the
reference version, from the scan of the reference FASTA, when that is within the search's
mismatches, and not the variant version or its carriers. So an index of a larger budget is not a
strict superset of a smaller one: a search that asks for more bulges can return fewer such
sites. For `NRG_3_hg38` that is 546 + strand and 504 − strand PAM sites over hg38's 1,580
stretches between `N` runs; of the 1,050 such NRG sites (535 on primary contigs), none lies
within 4 mismatches of the nine guides checked (TRAC, `ACTGAAATCTGTAAGCAGGC`, BCL11A sg1617,
EMX1, HBG1, CCR5, VEGFA site 3, FANCF and HEK4), and 9 of 20,000 random 20-mers have one at up
to 4 mismatches. The exposure that matters is a guide whose own on-target sits at such a
boundary: its report then holds no perfect match. Two more engine limits at contig ends, as
measured on hg38: a − strand PAM at a contig's first base is reported on the + strand at base
*N* (11 NGG and 37 NRG contigs open so, all alt, random or chrUn contigs; the primary
chromosomes open with `N`), and a + strand PAM on a contig's last bases is found by no search,
indexed or not (8 NGG and 31 NRG contigs end so, none a primary chromosome, and no fake-indel
contig, which all end in `N`). They need an engine change, and the indexes rebuilt with it.

**Cost.** A PAM-less index is much larger, because no site is excluded at build time: the
reference index is 71 GB against 9 GB for a dedicated NRG index, and the HGDP + 1000
Genomes variant index 110 GB against 28 GB (download sizes, GB = 10⁹ bytes). Pruning
happens at leaf time, so a PAM-less search is not faster than a dedicated one — the gain
is that one index answers many PAMs.

**Guide length.** An index name carries no guide length; each tree records the one it was
built for, as the second int32 of its `.bin` (written by `buildTST` since 2019, whoever built
the index). The engine reads that length only to offset coordinates and searches any guide it
is given, so past what the tree stores it gives wrong alignments or loses targets, with exit 0.
A search therefore takes an index only where

    built + N >= guide + bDNA + (1 if bRNA > 0, else 0)

with N the number in the index's name. We measured the rule with the real CLI against an index
dedicated to each guide, on a synthetic genome with no site near a contig end or an N run: of
246 cells (built 20, 22 or 23 nt; N 2 to 4; guide 20 to 24 nt; dedicated, all-`N`, and all-`N`
under `--prune-pam`; up to 3 bulges), the 136 inside the rule are byte-identical to the
dedicated index and the 110 outside it all differ, up to an empty result. Another 70 cells, of
variant searches (SNP and indel rows, with no bulges too, through the index) and of a 5′ PAM,
fall on the same side of the rule. So a published 20-nt N = 3 index serves a 21-nt guide at up
to 1 bulge of each type, a 22-nt guide with one bulge of one type, a 23-nt guide with no bulge
through the index, and no longer guide. An index built for this very guide is preferred to one
that merely fits. With none that fits, a search builds its own under the name
`<PAM>_<N>_<genome>`, or stops with an ERROR that names both lengths when an index of that
name, built for another guide, is in the way; `build-index-only` does the same and leaves that
index and its `.pam_build` as they are.
A longer build serves a shorter guide, but, like any deeper index, it does not hold a site
next to an N run or a contig end that a shallower one does.

**Limits.** `--prune-pam` requires `--index-path` and is command-line only: the web
interface never passes the flag. (Its PAM picker prefers a PAM whose geometry matches the
installed index, but falls back to offering every PAM when none matches, so it does not
itself prevent an unusable pairing: the search does not take an all-`N` index of another
width without the flag, or one built for the other side, and when no index serves it and
neither the dataset's VCFs nor its enriched genome are there to build one, it stops with an
ERROR that names what was asked, the variant indexes installed and what to do.) A shorter PAM
off a wider flank is supported for a 3′ PAM only; for a 5′ PAM the engine requires the flank width to equal the
PAM length, and the index must have been built from a 5′ all-`N` PAM file, as its
`.pam_build` records (above). The published
PAM-less indexes are 3′ builds, so they serve no 5′ PAM: a Cas12a search needs a dedicated
index. A flank wider than
3 nt is not readable by CRISPRme versions before this one, so those indexes are published
in a separate repository rather than alongside the standard ones.

---

*Software: CRISPRme+ (`pinellolab/crisprme-plus`). This document tracks the
methods as of the 2.8.x line (CRISPRme+ 2.8.1) (default SNP+indel co-occurrence, context-aware
analysis mode — per-sample on genotyped+PHASED indexes / population-level otherwise, with opt-out
`--population-level` / forced `--per-sample`, and the CRISPR-Bulge ML off-target
scorer in a dedicated conda environment with a `cpu|cuda|metal` compute backend); see the
CHANGELOG and the referenced source files for implementation detail.*
