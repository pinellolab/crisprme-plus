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
catalog name `NRR_3_hg38+hg38_HGDP1kGP` (SpCas9 NRR, the default PAM since 2.10.0; NRG
before): 925 HGDP + 3,166 1000 Genomes = 4,091 samples,
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

  *Known limit (engine), up to CRISPRitz 2.9.1.* The enricher of CRISPRitz 2.9.1 and
  earlier writes one IUPAC code a position and one SNP entry a position, rebuilt from
  the reference base by each record, so at a position with several biallelic SNP
  records an index it built, as every index published until the rebuild for CRISPRme
  2.10.0, searches only the **last** record's ALT: an earlier record's ALT is in neither
  the enriched genome nor the dictionaries nor the Tier-0 registry nor the Tier-1
  genotype store, and an off-target that needs it is not found, or is reported at a
  lower score, without its carriers (a VCF line that lists several ALTs itself is coded
  with all of them). Every release up to 2.8.1 behaves the same (the SNP dictionaries
  kept one record a position, the enriched genome since CRISPRitz 2.6.0), and so does
  2.10.0 on such an index.
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
  stores, and only the last record's ALT (AC 1 to 3) is held. CRISPRitz 2.10.0, the
  one CRISPRme 2.10.0 installs, codes every ALT of such a site (`Q`, any base, where
  the reference base and its ALTs cover A, C, G and T) and lists each in the
  dictionaries with its own record's carriers, rsID and AF, so an index it builds, as
  every one rebuilt for CRISPRme 2.10.0, holds them all. Two records of one REF>ALT at
  a position (362 alleles in the published panels, 351 of them pure duplicates) are one
  allele there: CRISPRitz 2.10.0 lists it once, with the carriers of both, a sample both
  list carrying the ALT on every haplotype slot where either puts it when its two calls
  are phased alike (as many alleles, every separator `|`: `0|1` and `1|0` give `1|1`),
  and otherwise keeping the call with more copies of the ALT, the earlier record's on a
  tie, as it is (CRISPRitz's `docs/MULTIALLELIC_ENRICHMENT.md` §2.3). So too where a VCF
  line of several ALTs lists it, whose calls of it are first written as a split record's
  (that ALT `1`, every other called allele `0`: a `1|2` reads `0|1` for its second ALT).
  Up to 2.9.1 the later record's entry replaced the earlier one, whose carriers were lost.
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

The **default** shipped index, **`NRR_3_hg38+hg38_HGDP1kGP`** (the SpCas9 NRR PAM's, the
default since 2.10.0; the NRG and NGG indexes of the panel are built from the same data), is
built in this genotyped mode
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
sources. A registry compiled from a genotyped panel's dictionaries also lists, in its
manifest (`four_base_positions`), the positions whose reference base and SNP ALTs hold all
four bases, an ALT that no sample carries counted: CRISPRitz 2.10.0 codes such a position Q,
any base, and codes an ALT that no sample carries as any other, but the registry holds no
record of such an ALT (it has no group). The SNP post-analysis of a genotyped panel keeps, at a
target's column, only the ALTs the column's code holds: its dictionary and the tiers compiled
from it list what the enrichment coded, so an ALT the code does not hold is another
position's. It reads a Q column's ALTs only at such a position (a Q at a bulge boundary can be
a neighbouring position's), so on an index with the tiers alone it takes the position from
this list. A sites-only registry is read from its VCF, every record of it, and lists ALTs the
enrichment never coded: on an index built by CRISPRitz 2.9.1 or earlier, which codes the last
PASS record of a split site alone (Homogenization, above), the site's earlier records, and on
any index non-PASS records. Each is the panel's own allele at its own position, so the
post-analysis of a sites-only dataset substitutes every ALT its registry lists at a variant
column, whatever the column's code, a Q included, as up to 2.8.1. A target that only such an
ALT brings within the search's budget is still not found on such an index; a mega index
rebuilt by CRISPRitz 2.10.0 codes every PASS ALT.

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
than the VCF's first record on an index built by CRISPRitz 2.9.1 or earlier (a known limit
of its `add-variants`, §8), and none for any
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
(`complete-test`'s benchmark mode: `20-NGG_2_hg38_chr22+hg38_1000G`, Section 5), a
`--per-sample` search of the README's guide with the `20bp-NGG-SpCas9` PAM, 4 mismatches, 1
DNA and 1 RNA bulge gives 149 of its 184 variant rows a blank MAF (3 of 184 with the tiers),
and four two-variant rows the frequency of their commoner variant (0.09 to 0.86) for a
haplotype that one to five of the VCF's 2,548 individuals carry. The report of a search on
such an index says so, in either mode; rebuild the index with `build-index-only --samplesID`,
or use a prebuilt one.

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
- a haploid call — a male's on chrX outside PAR — is his one haplotype's, so his calls in a
  window are cis by construction (they make one enumerated haplotype, not every subset of
  them, below) and an indel of his that changes the target takes him off the row, whatever
  else that haplotype holds;
- a row left with no carrier is not reported. The same holds on the classic per-sample
  dictionaries. The population-level analysis (§5) is unchanged: its rows are worst-case
  composites, whose variant combinations the carrier gate checks (§5). With
  `CRISPRME_INDEL_SNP=0`, the opt-out of the SNP+indel co-occurrence, the store is not read,
  and a SNP row lists its carriers as up to 2.8.1.

On chr22 of the HGDP + 1000 Genomes panel, the README guide, NRG, 4 mismatches + 1 + 1, per
sample, up to 2.8.1 3 of 628 highest-CFD SNP rows of the integrated results and 7 of 399 of
the alternative alignments listed such carriers (4,234 listings), and 7 of them had no other
carrier: chr22:21052638 − listed 1,292 carriers of 21052641 G>A, which sits on the same 1,448
haplotypes as 21052640 T>TCACC. Genome-wide, 280 of 65,184 highest-CFD SNP rows (28,700
listings, 13,525 of them carriers homozygous for the indel), and 208 had no other carrier (none
at CFD ≥ 0.2, highest 0.112). Over five guides on chr22, per sample, 2.2% of SNP
row-blocks (990 of 45,921) named such a carrier up to 2.8.1, and 545 of them had no listed
carrier that holds the sequence.

### Behavior by phasing status
The confidence attached to an enumerated haplotype depends on the phasing of the
input genotypes:

- **Phased data** (alleles separated by `|`). The cis/trans configuration is
  known, so a per-individual, per-chromosome variant-set is a **confirmed** cis
  haplotype: that individual carries exactly those variants together on one
  chromosome. A **haploid** call (one allele, no separator: a male's on chrX
  outside PAR) is read the same way since 2.10.0 — it is his one haplotype's, so
  there is nothing to phase — and a missing call (`.`) carries nothing.
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
  is surfaced in the high-variant-density BED). A haplotype holds one allele of a position,
  so no subset and no union takes two ALTs of one: where the individual carries two or three
  ALTs of a position (a `1/2` heterozygote, or a phased slot holding two), there is one union
  per choice of ALT while that makes at most 16, and past that one union, holding at each such
  position the ALT nearest the guide (one the PAM allows, then the fewest mismatches).
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
per indel (25 bp of reference, the indel's ALT, then 27 bp of reference, an `N` between two
windows), onto which the build writes the IUPAC code of every SNP of the enriched genome.
A target is searched in its window only where its own bases and its site lie there, the
site being its PAM and the guide's length of bases from it, as CRISPRitz indexes a site, not
an alignment (an RNA bulge leaves the site's far end past the target's first base, a DNA
bulge the target's first base past the site's). The window holds every SpCas9 site (23 nt)
that holds a base of its record, with its targets of up to 3 DNA bulges, and a 27-nt site
(SaCas9's 21-nt guide and `NNGRRT`, SaCas9-KKH, Cas12a `TTTV`/`TTTN` at 23 nt) only from its
first base on: a known limit for such motifs, whose target with a DNA bulge, from an
insertion's first added base, needs 26 bases before it (from 3% to 20% of the indel-holding
alignment slots of those enzymes at 1 to 3 DNA bulges; SaCas9 4.2% at one), and whose site
ending on an MNV's first replaced base starts a base before the window. Such targets cannot
be searched, whatever the index. An index of format 1 (CRISPRitz 2.9.1) also lost a target
whose site's guide region ends within its budget of the window's edge, the `N` in its pad:
every budget of a 27-nt motif, so for SaCas9 the sites that start within the budget of a
window's first base (at it, for a search without bulges: 2,103 carrier targets at 6 records
on chr22 of HGDP + 1000 Genomes for one guide, 2 mismatches); an index of format 2 (§9)
holds them. The search engine matches a code
when **any** of its alleles fits the guide (or the PAM), so one raw indel target stands for
every combination of the window's SNP alleles at once; only the genotypes say which of them
exist. The indel post-analysis therefore resolves each such target **per carrier
haplotype** (`indel_snp_cis.resolve_target`):

- every haplotype that carries the indel gets **its own allele** at each code column — the
  SNP alt it carries, else the reference base — and the haplotypes left with the same
  sequence form **one reported row**, listing only their individuals and naming the SNP alts
  they hold next to the indel (Variant info, rsID, MAF). A SNP with several alts on one VCF
  line keeps its raw genotype under each alt (2|1 under A and under T), so each genotype is
  first read against its own alt's allele index: the index every carrier of the alt holds,
  less, where that leaves several, the indexes the position's other alts settle (each alt of
  one character the dictionary lists there, the `*` of a spanning deletion among them; an
  indel's ALT is in no entry); where the genotypes still cannot tell (every carrier of two
  alts is 1|2, or of an SNV and an indel of its line) no slot can be tied to the alt and the
  carrier is resolved as an unphased one;
- a sequence whose PAM no longer satisfies the searched motif, or whose mismatches,
  recounted on the concrete sequence, exceed `--mm` or `--max-total-edits`, is **dropped**:
  a raw target that no individual carries gives no row;
- an individual whose genotypes cannot put an allele on the indel's haplotype (an unphased
  call, a missing allele, no phased indel genotype, a haploid call against a diploid one, two
  ALTs of one position on the indel's slot, which is no haplotype) gets the **strongest target
  its own alleles allow** on each alignment (the guide's base wherever it carries it and the
  alignment pairs the position with a guide base; elsewhere, a DNA bulge included, and at the
  PAM jointly, the alleles that maximize CFD among those the motif allows; a search that
  scores no CFD, with a 5′ PAM or a guide or PAM not 20 and 3 nt long, takes there the
  reference base where its alleles and the motif allow it, else the first of A, C, G, T they
  allow), PUTATIVE in the co-occurrence companion: it is never listed for an allele it does
  not carry, its ALT copies at a position counted against its call (the reference stays
  possible unless they fill its slots: two ALT copies of a diploid leave it out). The carriers
  of an equal-length multi-base record (an MNV, AC>TG),
  which CRISPRitz puts in a fake-indel window as it does an indel, are always in this case: the
  phased indel genotype store keeps only records whose REF and ALT differ in length (and,
  since 2.10.0, whose FILTER is `PASS` or `.`, as the enricher reads: §8). With
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
creates the PAM where the reference does not — with those rows' carriers (a line counts them
all and lists the first 50 by name): CONFIRMED when every carrier was resolved exactly,
PUTATIVE otherwise. It is written before the merge picks each cluster's alignments, one line
per engine alignment, required SNP set and phase (an alignment that needs no SNP writes no
line), so a line can name an alignment the merge did not keep: a population-level run,
which writes no alternative alignments, can have lines with no row in its tables. The
resolution does not depend on the analysis mode (§5).
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
of both, where 2.8.1 listed each copy at its own start; its `PAM_creation` is created when
either copy's PAM is (below). A bulged alignment of the second copy,
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
   motif. A tandem duplication of the PAM's bases is this case: the "inserted" PAM is the
   reference's own, shifted, and it stays NA;
2. else, a PAM that every placement makes wholly of the bases the insertion adds: an
   additional PAM, which the reference does not have there: created, even where the reference
   holds a valid PAM beside the protospacer (decision D4, 2.10.0: a new PAM within the
   search's mismatches and bulges is reported as one; before D4 such a PAM went to step 3, and
   a valid reference PAM beside the protospacer read NA, since the site existed before the
   insertion);
3. else a placement that leaves the protospacer's base that touches the PAM on a reference
   base: the reference PAM beside that base, on the PAM's side, the PAM's place beside that
   protospacer. A deletion's PAM always has one, its bases all being reference bases; an
   insertion's has one when some placement puts the inserted bases at or past the PAM's
   first base (the PAM at the target's right end, forward: the + strand with a 3′ PAM, the −
   strand with a 5′ one) or ending at or before the PAM's last (the PAM at its left end): a
   PAM partly made of inserted bases, in that placement;
4. else, an insertion that holds that base, with the protospacer's end and the PAM's start,
   wherever it is placed: no place on the reference, created.

A complex record replaces REF bases and adds bases or loses some (A>T and inserted bases after
it: HPRC writes many insertions so). Its first min(len(REF), len(ALT)) changed bases stand in
the places of the REF bases they replace, as an MNV's, and are read as those REF bases, as a
SNP alt of the row is; its other ALT bases are added after them, or its other REF bases
deleted, and are placed as an insertion's or a deletion's are. So a PAM made of a complex record's added bases is judged as one made of an
insertion's: a target wholly inside them, its PAM and the protospacer's end added bases, is
created unless a placement of those bases through their repeat puts its PAM on the reference
(up to 2.8.1 such a PAM was read on the REF alignment's window, reference bases the target
does not hold, NA where they happen to hold a valid PAM). On HPRC's chr22 (8,022 complex
records with a net insertion, 1,055 with a net deletion), with the TRAC guide, three CFD-0
alignments wholly inside the added bases of three complex records whose added bases repeat
the reference after them are not created (chr22:37509610 (−), 42134819 (−) and 42135364 (−):
placed further right, their PAM is the reference's own GGG).

A row whose REF alignment is its ALT has its target, PAM included, on the reference at its
start: NA (since 2.10.0). A PAM with a bulge in it keeps the REF alignment's bases, and so does
a PAM holding an MNV's ALT bases, as up to 2.8.1 (an MNV moves no base: those are the reference
bases under the PAM). The rule reads the row's haplotype: a carrier's other indels are not part
of it (a carrier whose second indel changes the target's bases is not listed on the row,
above), and a repeat that runs past the target is read on the reference. A row that stands for
two copies of one target inside one insertion (above) is created when either copy's PAM is:
each copy is judged on its own, and a placement of the insertion can put one copy's PAM on the
reference's own (not created) while the copy the insertion adds has none (TRAC guide,
chr22:50193928 (−): a 96-base allele holds the CFD-0 alignment CTACTCAGgaGGgACA-GGCtGAG twice,
at 50193950 and 50193998 on its haplotype, where 2.8.1 listed each; placed further right, the
insertion leaves the left copy's PAM on the reference's own CTC, GAG on −, and no placement
takes the right copy's off it: the row reads GAG). The merge keeps that value whichever copy it
meets first.

Examples, SNP+indel oracle (`test/indel_snp_oracle`): a tandem duplication of a 27-base unit
holding a whole target, the copy written left of the reference's (PR, PM): the copy's target
starts inside the insertion, its PAM wholly inserted bases, and placed one unit right the
insertion leaves that PAM on the reference's own: NA. Where an insertion after the protospacer
copies the reference's own PAM (T, C6), NA likewise. C>CCAG after its anchor on the − strand
(QC), whose first base the reference repeats there, holds the protospacer's first forward base
in both its places: created. Not a repeat, C>CCA after an anchor C that starts a − strand row
(QN): the PAM is the anchor and two added bases, the protospacer's touching base a reference
base, and the reference holds CTC beside it (GAG on −): NA. An inserted TGG before the
reference's AGG on the + strand (T2), or TTCCA inserted between the reference's CCG and its
protospacer on the − strand (QW): the PAM is inserted bases wherever the insertion sits, an
additional PAM, created (TGG), though the reference holds a valid PAM beside the very
protospacer (decision D4; NA before it). CA>C inside the
PAM CCG on the − strand (QD) leaves the reference's CAG (CTG on −) beside the protospacer:
created. Inserted bases holding the protospacer's last four bases and the PAM on the + strand
(QA): created. Complex records: A>T on the + strand and C>T on the − strand, each with 33
bases added after it that hold a whole target (QX, QY): created, though the REF alignment,
read from the replaced base where such a row starts, holds AGG (+) and CCG (CGG on −) in the
PAM's columns; G>CCA on the − strand, the row starting on the replaced base, its PAM CCA (TGG
on −) that base and the two added ones, beside a reference CTG (CAG on −) (QZ): NA; TG>C on
the − strand, the T replaced and the G lost, the junction inside the PAM CCA (QV): the
reference holds TGA (TCA on −) beside the protospacer, created.

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
wholly inside a 37- or a 33-base insertion: now TAA, at 44115944 and 44115946). With such
rows at the anchor, reading the PAM where the VCF writes the insertion would not do either: a
repeat unit inserted next to a reference PAM would read as created, where a placement further
right puts that PAM on the reference's own bases (with NRR, chr22:25415044 (−), the same
alignment as its site's reference row; TRAC guide, chr22:49196564 (−) and the row of
50193928 (−): NA, as 2.8.1 read them). Genome-wide,
2.8.1's tables (2.8.1's code gives back every value in them) with each row's PAM_creation
computed by 2.10.0's code, highest-CFD alignments: README guide, per-sample, 99 go NA → created
(70 rows and 29 alternative alignments: a PAM across a deletion's junction on the − strand 56,
CFD at most 0.083; a PAM holding inserted bases on a row that starts inside the insertion 27, on
a row that starts on a reference base 16) and 4 created → NA (3 rows and 1 alternative
alignment); population-level 70 and 3 rows; TRAC guide (CTCTCAGCTGGTACACGGCA), per-sample, 113
(67 and 46: deletions 29, CFD at most 0.039; rows inside an insertion 58; on a reference base 26)
and none created → NA. Of these, decision D4 (step 2) accounts for 11 highest-CFD blocks, every
one a PAM wholly inserted beside a reference protospacer that held a valid PAM, at CFD at most
0.008: 2 README alternative alignments (chr11:35248792 (+), CFD 0), none population-level, and
on TRAC 2 rows (chr17:77419671 (−), 0.008; chr2:159346 (−), 0.003) and 7 alternative
alignments (chr1:246083493 (+), chr20:63226848 (+), chr4:183404212 (−) three,
chr5:157864188 (+) two), 5 of which 2.8.1 also read created; nothing else changes.

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
  lattice, so it needs no variant-density cap.) The default bulge budget is **1 DNA + 1
  RNA bulge**: a search given neither `--bDNA` nor `--bRNA` takes it (complete-search, up to
  the bulges an index it can use or build holds; assembly-search, since 2.10.0; the web
  form's simple and Advanced modes, whose Advanced panel offers up to 2 of each).
- **High-variant-density cap.** Dense windows, those spelling more sequences than a
  configurable budget, 2 to the power `CRISPRME_IUPAC_CAP` (default 2¹⁰), fall back to a
  bounded procedure instead of full enumeration, so a single pathological window cannot
  dominate runtime or memory. The number of sequences a window spells is the product, over
  its columns, of the bases each letter of the variant genome stands for: 2 for a
  two-allele column (R, Y, S, W, K, M), 3 for a three-allele one (B, D, H, V), 4 for a
  four-allele one (`Q`), and 1 for a base or an `N` (an assembly gap, which matches
  nothing). Its log₂ is the window's bits, F. Up to 2.8.1 the cap counted a window's
  ambiguous columns: the same number for a window of two-allele columns, but blind to the
  1.58 bits of a three-allele column and the 2 of a four-allele one, which coding every
  allele of a split multi-allelic site makes common, above all in the immunoglobulin loci
  (on chr22 of HGDP + 1000 Genomes, the variant search of five guides, NRG, 4 mismatches +
  1 DNA + 1 RNA bulge, returns 4,372,248 raw targets with CRISPRitz 2.10.0, 96.6% of them
  from the IGL locus, chr22:22.00-22.95 Mb, where 2.9.1 returned 206,883). Crucially, the
  bounded procedure still emits a **greedy minimum-mismatch representative** for the
  window — at each variant column it takes the allele that most lowers the mismatch count,
  which (mismatches being additive per column) is the exact argmin over all 2ᵏ
  combinations, i.e. the window's worst-case off-target. So a capped window **surfaces at
  least one off-target row**, unless `--max-degeneracy` (below) cuts that representative
  back past the search's limits; the cap trades exhaustive per-haplotype enumeration for a
  single conservative representative, never the whole region.
- **Hyperdense windows, flagged whatever the cap.** A window spelling 2²⁰ sequences or
  more is hyperdense: it matches nearly any guide somewhere, as uninformative as a gap, yet
  the engine reports it, as it does any window it aligns to the guide within the search's
  mismatches and bulges (each of the 381 hyperdense windows the README guide's search of
  chr22 lists on the published index is aligned with 2 to 4 columns that do not admit the
  guide's base). Such a window is listed in the BED below whatever `CRISPRME_IUPAC_CAP` is,
  read from its own letters (so a window whose only variant letters are four-allele
  positions is listed too), and the report's `High_complexity_region` column gives its count
  (`24 in window (2^24.0 sequences)`).
  Among the indexable NRG windows of chr22 (holding no `N`, their letters able to spell NRG
  on either strand) the reference genome has none, the variant genome of HGDP + 1000
  Genomes published with CRISPRitz 2.9.1 has 1,763 (the worst 2^23.00), and the one
  CRISPRitz 2.10.0 writes from the same VCF, every allele coded (`Q` at its 2,448
  four-allele positions), 5,103 (2,266 of them at 2³⁰ or more; the worst 2^44.75, at
  chr22:22,881,189 in IGL). The threshold, 2²⁰, is a judgment call
  (`window_density.HYPERDENSE_BITS`, the one constant to tune). At the default cap every
  hyperdense window is dense as well. Within a cap raised to its bits, or with
  `CRISPRME_IUPAC_CAP=-1`, a `--per-sample` search enumerates its combinations as any other
  window's, so its rows are its carriers' own haplotypes: its BED line says no greedy
  worst case stands for it, and its rows are not noted as one.
- **Four-allele positions.** Where the reference base and three ALTs cover A, C, G and T,
  the variant genome writes `Q` (CRISPRitz 2.10.0), which matches any base; `N` stays an
  assembly gap. The post-analysis lists the four-allele positions of the windows the
  search reported in `<out>.four_allele_positions.bed`, and `High_variant_density_region`
  names them on the row whose window holds one (`four_allele_positions=chr:pos`, 1-based;
  the report: `1 four-allele position`).
- **Density reporting + transparency.** Every dense window, and every hyperdense one, is
  written to a `high_variant_density_regions.bed` sidecar (region span; variant count,
  every variant column of the window, a four-allele one too; carriers; the full IUPAC
  protospacer; and, since 2.10.0, its bits, its number of four-allele positions, and
  whether a greedy worst-case alignment stands for it, as it does in a population-level
  search), and each affected off-target carries a `High_variant_density_region` column in
  `integrated_results.tsv` giving the window's variant count and sequences and the full
  IUPAC sequence so a user can dig into them, and noting how the window was reported:
  where a greedy worst-case alignment stands for it, that one is reported and additional
  haplotype alignments may exist; under `--per-sample` (since 2.10.0) that the haplotypes
  its carriers have are reported (below), but for a row that names no carrier, a sites-only
  dataset's, whose rows are the population-level analysis's in either mode and keep the
  first note. Nothing is silently truncated — the bound, the representative, and the
  alternatives are all auditable. `complete-search` says per chromosome how many windows
  were dense and how many hyperdense (stdout, `log_verbose.txt`).
- **The ambiguous columns an alignment rests on, D, and `--max-degeneracy`.** D counts the
  columns of a reported alignment where the reference base would not fit (a mismatch
  against the guide, or a base the PAM does not admit) and the row's allele does. Few
  alignments resting on two or more are carried by anyone: on chr22 of HGDP + 1000 Genomes
  (the published 2.9.1 variant genome; five guides, NRG, 4 mismatches + 1 DNA + 1 RNA
  bulge), a real phased haplotype carries 96.8% of the alignments with D = 1 (44,315 of
  45,789), 5.6% with D = 2 (1,091 of 19,561), 1.1% with D = 3 (78 of 7,369) and 0.02% with
  D of 4 or more (10 of 59,208). `complete-search --max-degeneracy K` (off by default)
  caps D for the SNP rows: a worst-case alignment (population-level, or a dense window's
  greedy one) is cut back to the strongest resting on at most K columns (a valid PAM
  first, then the fewest mismatches or the highest CFD), and any other row resting on more
  is left out, unless its carriers hold its variants together on one phased haplotype
  (a CONFIRMED observed haplotype, or a phased dictionary's combination). The reference
  row keeps the locus where the reference alignment is within the search's limits; a site
  that only such alignments reach goes, and a capped window can be left with no row. It
  is a precision knob, not a gate: an alignment resting on two columns is rare but real
  for some. `complete-search` says per chromosome how many were cut back and how many left
  out, and the report states the cap. An indel row is resolved per carrier haplotype (§4)
  and is not cut. On chr22 of the published index, the README
  guide's population-level search (4 mismatches + 1 DNA + 1 RNA bulge) reports 2,197 rows;
  681 of its 1,257 SNP variant rows rest on two or more columns, the second-highest by CFD
  (0.299, in IGL) on eleven. With `--max-degeneracy 1` it reports 1,572: 625 sites go, those
  only such alignments reach (9 more move by 1 or 2 nt, to where their cut-back alignment
  starts), the rows lying in a dense window fall from 281 to 2, and 47 to 49 of the 100
  sites with the highest CFD leave the top 100, as the ties at the 100th CFD fall. Its
  per-sample search changes no row, its multi-variant rows being CONFIRMED haplotypes.

Together these keep genome-wide variant search tractable while making any bound
that was applied explicit and reviewable, and, unless `--max-degeneracy` is given,
guaranteeing no region is dropped.

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
its reference row and its BED line (which says that no greedy worst case stands for it).
That is the observed-haplotype path, which reads the index's dict-less tiers; since 2.10.0 a
capped window on an index without them (one built without `--samplesID`, such as the one the
benchmark `complete-test` builds) goes the same path, reading the same genotypes in the
classic per-sample dictionaries. Up to 2.8.1 such an index's classic dictionary path still
reported a capped window by its greedy worst-case alignment, flagged in
`High_variant_density_region`, whose listed carriers (the carriers of any allele it takes)
need not carry every variant in it. On the NGG index the README's self-test builds
(`complete-test`'s benchmark mode: `20-NGG_2_hg38_chr22+hg38_1000G`, named
`NGG_2_hg38_chr22+hg38_1000G` up to 2.8.1; chr22 of the benchmark's phased 1000 Genomes VCF,
built without `--samplesID`), the README's guide (`CTAACAGTTGCTTTTATCAC`) with the
`20bp-NGG-SpCas9` PAM, 4 mismatches, 1 DNA and 1 RNA bulge, searched per sample at the
default `CRISPRME_IUPAC_CAP`, gave 5 such rows with 2.8.1 (8 to 13 variants each, CFD at most
0.001). They named 102 individuals as carriers 114 times (some on more than one row), and
each time the carrier lacked one of its row's variants: 342 of the run's 56,730 SNP-row
carrier listings (a row lists its carriers once for each of its three alignments).
After `build-index-only --samplesID` the same search returned exactly those 5 rows fewer, and
every carrier it lists has all of its row's variants on one phased haplotype; since 2.10.0 the
index without the tiers gives the same rows (560 against 565, the five rows above gone). The
2.8.1 report said so of a 2.8.1 search's results, and still does when it reads them again (a
search since 2.10.0 writes the carrier gate's sidecar, `.carrier_gate.json`, in either mode).

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

### The carrier gate: a population-level row's variants against phased haplotypes

A population-level row is a window's worst case: its alignment takes, column by column, the
allele that matches the guide best, and on a genotyped panel it lists every sample that carries
any of those alleles. Two or more of them can sit on no one haplotype, so that no individual has
the sequence the row reports, and such a row, being the best-matching, ranks first. A sites-only
panel (`mega`) has no genotypes at all, and its rows are built the same way. Since 2.10.0 the
search checks each row's variant combination against phased haplotypes after it has written the
integrated results (`haplotype_prioritizer.py`), and the report ranks by the result. Detection is
unchanged: the gate adds and drops no row, and the integrated results are byte for byte those of
the search without it (`--no-prioritize-observed`).

- **Which haplotypes.** A genotyped dataset's own phased genotypes (`Dictionaries/genotypes_<vcf>/`
  for the SNPs, `indel_genotypes_<vcf>/` for the indels), by default in a population-level
  search. A sites-only dataset's rows are checked, on request (`--prioritize-observed`; the web
  form's check box), against the phased HGDP + 1000 Genomes genotypes (4,091 samples, 8,182
  haplotypes on an autosome): the genotype oracle, which an installed
  `NRG_3_hg38+hg38_HGDP1kGP` index brings, or `download --what genotype-oracle` installs
  alone (`CRISPRME_GENOTYPE_TIER` names another folder). Without it the search says so and the
  rows are reported unchecked. A per-sample search needs no gate: its rows name their cis
  carriers. A row does not say which dataset it comes from, so the gate reads every phased
  store the search has: a variant is in a panel when the panel holds its POS, REF and ALT (a
  REF that differs is counted, and said, as a normalization difference).
- **Tiers**, for each of a row's three alignments (highest CFD, fewest mismatches + bulges,
  highest CRISPR-Bulge):
  - **A**: two or more variants that one phased haplotype of the panel carries together — on
    one slot of a phased sample, with no indel of the store on that slot that takes the
    target's bases off it, as for the per-sample rows (§4), read on the search's reference
    and the row's variants: one that changes the target's bases does, unless they are where
    they were (a repeat unit gained or lost at the target's edge), and so does one whose REF
    holds a variant's base at another position (two records of one base, whose bases cannot
    be read); one at a variant's own position, another allele of that base, takes nothing
    off (a phasing artefact). States: `one haplotype` when a single one does; `phase not
    known` when only samples whose genotypes cannot say which haplotype holds which variant
    carry them all (an unphased call, or two ALTs of one position on one slot, the PUTATIVE
    rule; an indel at a variant's own position on the slot is such an ALT).

    **A also where a haplotype reaches the alignment's CFD without every variant.** The
    worst case takes the best-matching allele in every column, and some of those alleles
    leave the CFD the results report (three decimals), the alignment's mismatches and its
    PAM as they are: one at a PAM base the search's PAM leaves free, or a mismatch for a
    mismatch the CFD table scores alike. The gate scores the alignment from its aligned
    spacer+PAM and protospacer+PAM (`bulge_cfd_penalty`, as the search scores a row), takes
    each allele in turn back to its reference base, and keeps the variants that change the
    CFD, the alignment's mismatches, or its PAM — the guide's side of an alignment writes
    its PAM as `N`, which every base matches, so the PAM region is read against the
    search's own PAM instead (its PAM file, `pam_filter`, as the post-analysis gates a
    target's PAM): without an allele the PAM needs there is no off-target at the site at
    all. If a phased haplotype carries the variants kept, it has the row's off-target at
    the row's score, and the row is A (states `same score`, `same score, phase not known`,
    and `same score as the reference` where no variant changes the CFD), with the variants
    the score needs in the companion's `Scoring_variants` column and in the label. Were it
    C1 the whole site would leave the top views, though a named individual carries an
    off-target there with the row's score. Checked only where the gate reads the
    alignment's CFD: a 20-nt spacer with a 3-base 3' PAM, the row's own CFD where the
    alignment has one, every variant a SNP inside the alignment, and the alleles the
    alignment takes. A CRISPR-Bulge score is not recomputed (the model is not run again),
    and an allele the CFD never reads does change it, so this state is the CFD's to give,
    on the alignments it scores: the alignment ranked by CRISPR-Bulge score, which has no CFD
    column, is C1 where no haplotype carries its combination, and its row leaves the views
    ranked by that score.
    Run without the search's PAM the gate cannot read a PAM region, and takes every allele
    in one for an allele the row needs.
  - **B**: one variant. The panel, or the frequency the dataset gives it, says that some
    individual carries it, so it cannot be a combination nobody carries; that is all B says
    (call quality is another question).
  - **C1**: two or more variants, each in the panel, that none of its haplotypes carries
    together and none of which reaches the alignment's CFD without (above); or variants (one
    or more, a single one too) every phased haplotype carrying which also holds an
    indel that takes the target's bases off it, as A reads one (state `indel in the
    target`). A zero means something only where
    some were expected: E, the panel haplotypes expected to carry the variants together if
    they were independent (the panel's haplotypes times the product of the variants' panel
    frequencies; linkage is not modelled), is given, and below a threshold (3 by default,
    `--gate-min-expected`; 0 turns it off) the state is `too rare to tell`. Ranked after the
    others all the same.
  - **C2**: two or more variants, one or more of them not in the panel: absence of evidence,
    not evidence of absence, as the combination may be a real haplotype of an ancestry the
    panel has few samples of. Kept in its place. The safeguards' frequency bound (the
    product of the variants' frequencies at least 1e-06 and the rarest at least 0.001) is
    given as information and applied to nothing: a SNP's frequency is the dataset's registry's
    (on `mega`, `AF_max`), an indel's that of its indel frequencies (`indel_af_<vcf>/`), or
    else the row's own.
- **Output.** `..._carrier_gate.tsv` beside the integrated results, one line per alignment
  that uses a variant (the row's line, the alignment, the variants, tier, state, panel, the
  haplotypes that carry them together and their samples, E, the frequency product and
  rarest frequency, and the label), and `.carrier_gate.json` (how the gate ran, its panels
  with their samples by superpopulation, and its tiers by each row's highest-CFD alignment).
- **Report.** A C1 row is ranked after every other row in the views ranked by the score of
  the alignment that is C1. A row whose highest-CFD alignment is C1 leaves the top-1000 by
  CFD, panels (a) and (b) of Section 2, the validation panel's CFD thresholds and the risk
  card's CFD bands, and is listed in a table of its own after the top-1000, by CFD, with all
  of them in `unobserved_worst_case.tsv`; a row whose highest CRISPR-Bulge alignment is C1
  leaves the top-1000 by CRISPR-Bulge score, panels (c) and (d), the panel's CRISPR-Bulge
  thresholds and the card's CRISPR-Bulge bands. The two alignments of a row need not be the
  same one, so a row of the demoted table can keep its place in the views ranked by
  CRISPR-Bulge score, and the other way round. The summary matrix, the validation panel's
  totals (its off-targets, perfect matches and variant-created off-targets, beside a row
  that counts the C1 rows), `run_manifest.json` and the population plot count every row. A
  demoted site can still have a real off-target there with fewer of the variants and a lower
  score, which a `--per-sample` search of a genotyped panel reports as its own row. The
  `Observed` column of a row the gate checked gives its tier's factual label (`together on N
  of H <panel> haplotypes`, `the same CFD on N of H <panel> haplotypes, which carry
  <variants>`, `not together on any of H <panel> haplotypes (expected E)`, `not checked: K
  of N variants not in <panel>`), with no "chimera" word; a row of one variant reads its
  carriers as before, unless every haplotype carrying it has an indel in the target, when it
  reads that label (`all N <panel> haplotypes carrying it have an indel in the target`) and
  is demoted with the others. The inputs card has a **Carrier gate** row:
  the panels, the counts by tier, and the panel's samples by superpopulation (on HGDP +
  1000 Genomes: AFR 986, AMR 549, CSA 782, EAS 816, EUR 771, MID 157, OCE 30), as a
  haplotype found only in an ancestry with few samples there can be missing from it. The
  carriers it names are the panel's individuals, not the patient's.

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
cis carriers**, so `putative` cleanly marks reconstructions that no one actually
carries — separating genuine haplotypes from worst-case combinations. In a
population-level analysis a multi-variant SNP roster is the union of each variant's
carriers; where the carrier gate ran (§5) the cell gives the row's tier instead (`together
on N of H <panel> haplotypes`, `the same CFD on N of H <panel> haplotypes, which carry
<variants>`, `not together on any of H <panel> haplotypes (expected E)`,
`not checked: K of N variants not in <panel>`), and without it the cell reads `N carriers
(>=1 of K vars, poss. chimera)`. A row of one variant reads its carriers either way, unless
the gate found that every haplotype carrying it has an indel that takes the target's bases
off it (`all N <panel> haplotypes carrying it have an indel in the target`, C1); an indel
(or MNV) row's
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
`-1` and says so, a search in its `log_verbose.txt` and `targets-scoring` on stdout. The
released **amd64** Docker image bundles the CUDA TensorFlow build, so
`docker run --gpus all ... --compute-backend cuda` runs CRISPR-Bulge on an NVIDIA GPU with
transparent CPU fallback; the **arm64** image stays CPU.

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
`indel_snp_cooc.tsv` files are **byte-identical** (same MD5, 2,729 rows, 843 CONFIRMED / 1,886 PUTATIVE, the
same per-sample `cis_samples` and joint-AF in both). So per-sample cis attribution — which individual
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
On an index CRISPRitz 2.9.1 or earlier built (every one published until the rebuild for
CRISPRme 2.10.0), two SNP records at one position (a split multi-allelic site) keep only the
last record's alt in the enriched genome and its dictionaries (a property of its
`add-variants`), so a carrier of an earlier record's alt is read as reference there, in the
SNP and the indel post-analysis alike, however common that alt: on the HGDP + 1000 Genomes
panel chr8:4174067 T>C (8,170 copies) is read as reference because T>G (1 copy) is the later
record. That is 1.95 % of the panel's SNP records, whose effect §2 gives as measured.
CRISPRitz 2.10.0 codes every alt of such a site and lists each in the dictionaries (§2), so an
index it builds, as every one rebuilt for CRISPRme 2.10.0, holds them all. `add-variants` up
to CRISPRitz 2.9.1 also looks for the `AF` key in the VCF's first record only, and reads every
record's `AF` at that place in `INFO`: a record that has it at another place gets none, in the
index's indel log and in its classic per-sample dictionaries. Its indel then has no frequency
on any index such a CRISPRitz built (§3), nor its SNP on one without the dict-less tiers (the
registry, which gives a SNP's frequency with them, counts AC/AN from the genotypes, or on the
mega reads each source's `AF_<dataset>` by its key), and the report's MAF leaves it out. A VCF
merged with `seq_script/merge_panels/merge_vcf_panels.sh` has such records, since a site of one
source has fewer keys before `AF` than a site of two: the retired
`NRG_3_hg38+hg38_1000G2021_HGDP` index (1000 Genomes 2021 + HGDP, built from such a merge,
retired in 2.7.8) has no `AF` for 14,567,429 of its 15,583,931 indels (224,981 of the 232,384
of chr22). Of the indexes published until the rebuild only HPRC's chrY is affected,
whose VCF has `AF` at five places in `INFO`, after one to five other keys: its log has none for
91,950 of its 102,665 indels. The HGDP1kGP index and HPRC's other chromosomes have an `AF` for
each of their 17,192,381 and 12,000,599 indels, as the self-test's 1000 Genomes index (chr22)
has for its 64,528 indels and 993,880 SNPs; the mega's VCF has no `INFO/AF` at all (§2), so
none of its 13,853,761 indels has one, and their per-dataset frequencies are in the
`indel_af.tsv` companion. CRISPRitz 2.10.0 finds each record's `AF` by its key, so an index it
builds, as every one rebuilt for CRISPRme 2.10.0, has the `AF` of each record that has one.
The phased indel genotype store
leaves out equal-length multi-base records (MNVs), so their carriers are resolved with the
phase unknown, PUTATIVE (§4): the HGDP + 1000 Genomes index has none, HPRC's chr22 has 9,079
among its 195,813 indel-search records. Since 2.10.0 it holds only the indels of records whose
FILTER is `PASS` or `.`, the records CRISPRitz's enricher reads (no other is in a fake-indel
window, a dictionary or a tier); up to 2.8.1 it held every indel of the VCF, so on a VCF with
records of another FILTER a carrier was taken off an indel row, for a second indel in its
target, or off a SNP row, for an indel no window holds. A store written so is written again by
`build-index-only` where its VCFs are (its `.done` marker does not record the FILTER rule), and
a search that reads one says so in its `log_verbose.txt`; the panels published hold no such
record. Since 2.10.0 a per-sample SNP row does not list a
carrier whose haplotype holding its SNPs also holds an indel that changes the target's bases
(§4); up to 2.8.1 it did, and so it does with `CRISPRME_INDEL_SNP=0`. What remains: an MNV
inside the target, which the phased indel store does not hold (the HGDP + 1000 Genomes index
has none); a carrier whose indel call or SNP call is unphased, which stays listed, PUTATIVE;
and the population-level analysis, whose rows are worst-case composites (the carrier gate
ranks a combination no haplotype carries after the others, §5).
On a row of several variants, the `Variant_rsID` and `Variant_MAF` lists of the results
tables are in the order of its `Variant_info_genome`, one rsID and one AF a variant, and two
variants of one frequency give it each: the merge of a cluster's rows
(`remove_contiguous_samples.py`) keeps each variant with its own rsID and AF. A variant that
two datasets of a search give with two rsIDs or frequencies lists both pairs, one after the
other; an rsID field holding two IDs (a VCF ID `rs1;rs2`) reads `rs1;rs2`, one item. Up to
2.8.1 the merge sorted the three lists each on its own, as text, and kept each value once: on
ACTGA (NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per sample, genome-wide, HGDP + 1000 Genomes)
953 of 1,967 rows of two or more SNPs and 118 of 293 rows of an indel and one or more SNPs
listed their rsIDs in another order than their variants. The report's MAF, the lowest
frequency of the list, did not depend on it.
An indel row gives its record's first ID, and so does each ALT of a VCF line that lists
several ALTs: CRISPRitz `add-variants` writes the first ID alone into the indel log for every
indel record, whatever the number of its ALTs, and into the per-sample dictionary for each
ALT of a line of several, so no reader of either can recover the others. A SNP of a record of
one ALT keeps the record's whole ID column, which CRISPRme gives as one item (`rs1;rs2`).
So where the ID column lists one ID per ALT (`rs1;rs2`, as `bcftools norm -m+` and
`gnomAD-converter --multiallelic` join them) the rows of the second and later ALTs give the
first's (`rs1`); and a record of a VCF split with `bcftools norm -m -any`, as the published
panels are, keeps the line's whole ID column, which its SNP rows give whole and its indel rows
by the first ID. The mega's VCF lists two IDs for 68,017 of the 222,283 PASS indel records of
chr22 (`22:22652455:TTATG:T;rs1372713682`, whose rows give `22:22652455:TTATG:T`), the
HGDP + 1000 Genomes panel's for 8 of its 240,557; the sites-only (mega) registry builder,
which names the mega's SNPs, gives each ALT of a line its own ID where the ID column lists one
per ALT. The input validator says so where it finds a line of several ALTs (in each VCF's
first 1,000 records on every run, anywhere with `--full_input_validate`).
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

`validate-test` compares the engine's **sites** — (chromosome, start, strand) — on the
reference and on the **SNV-only** 1000 Genomes chr22 enriched genome. It therefore says
nothing about the **indel path**: the brute-force references are built from a FASTA that
holds no indel, so the fake-indel genome, an indel row's coordinate, its aligned
spacer/target strings, its REF/ALT text, its named variants, its carriers and the whole
alternative-alignments table are outside the comparison, and so is the per-carrier truth
(the enriched genome collapses every individual into IUPAC codes). The synthetic SNP+indel
oracle (`test/indel_snp_oracle`) checks those on a planted synthetic chromosome, through the
real CLI, but only there.

**Validation (the indel path on real data).** `scripts/check_indel_rows.py` is an
independent checker of the indel path on a real panel, added in 2.10.0. It imports only
Python's standard library and no code of CRISPRme's variant or indel machinery, and copies
none: its inputs are what a user has — the reference FASTA, the VCF(s) the search read, the
samples list (the samplesID files; by default the run's `.sampleID.txt`, which CRISPRme
writes from them), the guide, the PAM file, the search budget and the result folder (the CFD
tables are read as data) — and it writes each sample's haplotypes itself from the VCF 4.x
definition of a record, keeps its own map from haplotype bases to reference coordinates, and
aligns with its own exhaustive aligner. **The truth is the VCF, as the search read it.** A
slot is the reference with every ALT allele its genotype puts there of the records CRISPRme
reads (FILTER, below), and it holds a row's target when it reads
the target's bases from the row's start through the row's record: each base the record
puts in the target at its place, and the base after them (a second insertion between the
start and the record lengthens the stretch and pushes some of them out of it, which makes
another target, one that needs two indel-path records). Two alleles of one slot whose REF
spans meet — a SNP at an insertion's own position, a SNP on a base a deletion removes, one
complex allele split into overlapping records, two ALTs at one position phased onto one
haplotype — describe one stretch twice, so the VCF gives no sequence there. Where every
allele that meets the record's REF is one §4 reads the slot without (another allele at the
record's own position, an allele only on bases it deletes, an MNV) and no allele kept
meets one left out, the slot is read without them: the one reading that holds the record,
and so the one a row naming the record and not them can mean. A listed carrier must hold
the target there, and one that does is a note (`overlap-reading-agrees`); one the row
leaves off while its slot holds it there (§4: it stays listed) is counted
(`limit-overlapping-records`). Any other overlap of the record (an indel at another
position, which §4 takes off, or a SNP there, which §4 does not cover) leaves the slot no
reading. Two other alleles of the slot whose REFs meet (an overlap zone) leave it no
sequence only where they can change what it says of the target — one of them is named, or
changes a base read from the target's start —: every reading of such zones (any of their
alleles no two of which meet, none included) is tried, and where all agree the slot is
judged so; a zone wholly before or after the stretch read changes none of its bases. Where
the readings part, the slot has no reading. There the VCF holds no truth, and §4's rules
are the only model: a listed carrier they take off the target for a documented cause (a
complex allele, a second indel the slot cannot be read through, the last record of a
split site, the PUTATIVE rule) fails (`carrier-off-by-rule`), any other is counted.
Whatever its overlaps, a slot that holds a named variant on no reading (neither on it nor
of unknown phase there) lacks the target. Where a call that touches the target has unknown
phase (unphased, a missing allele, or fewer copies than the record's call — a haploid call
against a diploid one), or the record's slot holds two SNV ALTs of one position that the
search's stores both read (split records phased one by one, every allele coded: the slot has
no one sequence in the VCF there, and §4 resolves such a carrier PUTATIVE, decision B4), the
record's known slot is judged first — when no placement of
those calls lets it hold the target (a phased second indel inside it, say), or every
placement does — and otherwise the truth is the PUTATIVE rule of §4: the strongest target
the sample's own alleles allow on the row's alignment, the guide's base wherever it
carries it and the alignment pairs the position with a guide base, elsewhere (a DNA bulge
included, and at the PAM jointly) the alleles that maximize the alignment's CFD, its copies
counted against the record's call (a haploid ALT call beside a diploid record leaves the
reference possible); CFD is a product, so where it is 0 whatever the choice (a DNA bulge
whose base the table scores 0, say) every choice is a strongest one, and §4 orders no CFD
tie, so a tie is left open; a search that scores no CFD (a 5′ PAM, or a guide or PAM not
20 and 3 nt long) takes the reference base where the alleles and the motif allow it, else
the first of A, C, G, T they allow, read in the target's orientation. The rule chooses
only at the fake-indel window's code columns, the target's reference bases off the
record's own REF span: the window writes the ALT verbatim, the anchor included, and codes
flank bases only, so a SNV the VCF holds at the record's position (the other allele of a
multi-allelic heterozygote, one per haplotype) is never a choice there.
CRISPRitz's enricher reads only the VCF lines whose FILTER is `PASS` or `.` (those exact
strings; `PostProcess/validate_inputs.py`'s full check says the others "will be excluded
from enrichment"): no other line
is in the enriched genome, its dictionaries, registry and genotype store, the indel log or
the fake-indel windows. The phased indel store (`build_indel_genotypes.py`) holds the same
records since 2.10.0; one written before streamed the VCF without reading FILTER and kept
every indel, and §4's carrier resolution reads the store: a second indel or an overlapping
one on a carrier's slot takes the carrier off a target, and a SNP row drops a carrier for an
indel of its slot, whatever that indel's FILTER (the SNP+indel oracle, the second deletions
of its sites O and Q and the insertion of Y1 marked `LowQual`, lists its carriers so on such
a store). So of a genotyped dataset a line outside `PASS` and `.` is read as the reference,
but for its indels where the run read a store written before 2.10.0 (its `log_verbose.txt`
says so: the releases, below), which stay on their carriers' haplotypes as that store has
them; either way no fake-indel window holds one, so a row naming one fails
(`variant-filtered`) and completeness demands none of its targets (one
`limit-filtered-record` line per chromosome counts them), and its SNVs and MNVs are read as
the reference: a carrier of such a SNP stays listed for a target the SNP would change (the
oracle's site A, its PAM SNP marked `LowQual`, lists its cis carrier). Of a sites-only
dataset such a line's indels and MNVs are dropped (the `hg38_mega` index's indel log holds
exactly the 222,283 `PASS` indels of its chr22 VCF, none of its 40,840 others), and its SNVs
are read by its SNP rows whatever their FILTER: `hg38_mega`'s registry holds all 838,025
SNVs of its chr22 VCF, the 78,116 outside `PASS` among them, and its SNP rows name some
(chr22:11574550 C>G, `LOW_VQSLOD`), while an index CRISPRitz's enricher builds from a
sites-only VCF holds none of them, and a run does not say which it read, so a SNP row
through one is a row CRISPRme may write or not (below). Its fake-indel windows code what
CRISPRitz's overlay codes from the enriched genome: its SNVs whose FILTER passes, every
allele of a split site on a genome of CRISPRitz 2.10.0 and the last record of one on a genome
of 2.9.1 or earlier (on chr22, the published `hg38_mega`'s indel rows name the last record at
279 split sites and, at chr22:25093750, the `PASS` C>G before the `LOW_VQSLOD` C>T, and no
other), so a row of a sites-only dataset (no carrier) must name SNVs of that dataset's VCF
whose FILTER passes (`variant-not-in-dataset`, `variant-filtered`; on a genome of 2.9.1 an
earlier record of a split site is the note `variant-shadowed`).
**Releases.** Where CRISPRme's and CRISPRitz's releases read the data otherwise, the check
reads which one did from the run's own records, and refuses (exit 2) a run whose records do
not tell where its data makes it matter. A run of CRISPRme 2.10.0 or later writes
`.indexes_searched.txt` (the indexes that served it); one without it is an earlier
CRISPRme's, which read indexes CRISPRitz 2.9.1 or earlier built (format 1, the last record
of a split site) and a phased indel store of every FILTER. Of a run of 2.10.0: the variant
genome codes every allele of a split site where the variant index that served it flags Q in
its `.bin` headers, or where its enriched genome (`Genomes/<genome>+<dataset>/` of the run's
working directory) codes an earlier record's ALT at a split site of the VCF, and the last
record alone where that index is of format 1, or the enriched genome holds the last ALT
alone, or an index of format 2 flags no Q while the VCF has a position of four bases (a
line whose ALT is a `*` is no record of such a site: no code holds it); with every allele
coded, two listings of one SNV REF>ALT at a position, on a line of several ALTs too, are one
allele, their calls merged as CRISPRitz 2.10.0 merges them into the allele's one dictionary
entry, which the search reads (§2: taken from CRISPRitz's documentation, not from CRISPRme's
code), and a slot holding two ALTs of one position is resolved PUTATIVE (§4), the truth there
(above); with the last record alone coded, an earlier line of the same REF>ALT is an earlier
record of the site, whose carriers that dictionary lost (its entry was the last line's:
`limit-last-alt-wins`). The INDELS index's format is that of its `.bin` files, else the
run's format-1 NOTE of it (§9). The phased indel store holds only the records whose FILTER passes but where the run's
`log_verbose.txt` notes a store written before 2.10.0 for the dataset and chromosome (without
the log, a run whose genotyped VCF holds an indel record of another FILTER is refused). A
per-sample run of 2.10.0 writes no greedy row of a window over `CRISPRME_IUPAC_CAP`, so none
joins a cluster; its `--max-degeneracy K` (§5) cuts back or leaves out a SNP row resting on
more than K ambiguous columns, so in a cluster chain an alignment over more than K positions
of several options is a row CRISPRme may write or not (`limit-max-degeneracy`); an indel row
is not cut. The carrier gate (§5) ranks and labels population-level rows in its companion
and the report and leaves the integrated results and their alternative alignments as they
are, so it changes nothing the check reads. A run made with the SNP+indel co-occurrence off
(`CRISPRME_INDEL_SNP=0`, which a search of 2.10.0 notes in its `log_verbose.txt`) is
refused: every SNP code of its fake-indel genomes read as the reference, its indel rows name
no SNP and list every carrier of their indel, and its SNP rows keep a carrier whose
haplotype holds an indel in the target, none of the per-carrier resolution modeled here (a
run of an earlier CRISPRme made so records nothing that tells it).
CRISPRme's documented reading of a carrier (§4: the one-indel window, resolved with each
carrier's own SNP alleles, so that a row lists a carrier only where those give its target;
the second indel read from either end of the target through the indel, each base the indel
puts in the target at its place, before any unphased call; the phasing artefact and the
complex allele, the PUTATIVE carriers of MNVs and of multi-ALT alts the genotypes cannot
tie by §4's rule, the last record of a split site) is re-implemented from this text only to
classify a disagreement as a documented limit (counted, non-gating, a PUTATIVE tie included)
or a failure. An MNV beside the record is in no store CRISPRme reads (§4): a target that needs it
needs two indel-path records (the two-indel limit), while a carrier listed for a sequence
such an MNV on its slot changes is listed for a sequence it does not have, which nothing
documents, and fails. It runs two ways:

- **Soundness.** Every row block of the integrated results (the highest-CFD, fewest mm+b
  and highest-CRISPR-Bulge blocks) and of the alternative alignments whose variant list
  names an indel or an MNV must be labelled ALT, name one indel-path record (and SNPs),
  all records of the VCF the search read (one whose FILTER it skips is
  `variant-filtered`); its aligned spacer must be the guide with its PAM, its
  mismatches, bulge types and sizes and totals must recount within the budget, no gap may
  fall in the PAM, the PAM must match the searched motif and the REF text must be the
  reference read from the row's start, gapped as the target; its target (gaps removed, PAM
  included) must be the reference with exactly its named variants, read from the reported
  start (the first aligned base's reference coordinate, the insertion's anchor for a base
  an insertion adds, the last replaced base for one a complex record adds), and hold a
  base of its record. Each listed carrier must be a sample of the panel and hold the target
  on a haplotype at that start with every named variant — or, its phase unknown there, be
  listed for its strongest target; a carrier whose slot's alleles overlap the record is
  held to the slot's one reading where §4 gives it, and to the named variants on a slot
  that can be the record's always, and is otherwise counted — and each carrier of the
  record whose haplotype holds it (or whose strongest it is) must be listed: the merge
  never moves a carrier off its sequence's row (its key for an indel or MNV row holds the
  aligned target and, since 2.10.0, the aligned spacer, so a PUTATIVE carrier stays on the
  alignment it was resolved on). An ALT-origin block that names no variant is reported. A
  sites-only VCF genotypes nobody, and a sites-only dataset's rows list no carrier (§4): a
  row that lists none for a record a sites-only VCF holds is checked but for its carriers,
  each record deciding for itself, so that beside a genotyped dataset (the default index
  with `hg38_mega`) the genotyped dataset's rows are checked in full. A row that names an
  earlier SNP record of a split site, which the documented limit says the stores do not
  hold, is a note: the record is real and the target is checked with it.
- **Completeness (the real-data brute force).** For every indel and MNV record of the VCF
  in the checked region, each carrier's known slots around it are written with their
  known alleles (a target that touches a call of unknown phase is the PUTATIVE rule's, and
  no other target depends on where such a call is), identical sequences are aligned once,
  and every in-budget alignment on both strands with the PAM that holds a base the record
  changes (for a pure deletion, the bases on both sides of it — a target wholly in a flank
  is reference or SNP sequence and gives no indel row) is mapped to the reference and
  required to be reported for that sample by rows no worse than its best alignment by the
  merge's own sort keys — a fewest-mm+b block whose key under `--sorting-criteria`
  (mm+bulges, then mm, by default) is at most the least of the target's alignments', and a
  highest-CFD block of at least its CFD that, where the two print the same three decimals,
  has a key no greater than the least of the target's alignments of that CFD (the CFD file
  sorts on the CFD, then on the same criteria; where a run names
  `--sorting-criteria-scoring`, which a later branch hands the CFD file, a tie no worse
  under either is accepted) — since the merge keeps, of the rows of one site, variant list
  and carrier set, the first by each file's order (§4): a row with exactly that target and
  variants whose fold group (the rows of its site naming those variants with its carrier
  set, which two alignments of one target always share) is no worse, or rows of its site
  naming the same variants with the target's own carrier set that are no worse — every
  carrier both the documented reading and its own haplotypes list for it, and besides them
  only carriers either lists or the reading leaves open (the documented reading only
  classifies: a carrier it would list against the carrier's haplotypes vetoes no fold; a
  carrier whose slot overlaps the record is never required, and allowed where its slot,
  read as §4 reads it, holds the target, or, a slot with no reading, where the reading
  lists it or leaves it open). The fold key holds no strand, while one target's alignments
  share theirs, so a fold group is taken on one strand: a fold across strands loses another
  target. The merge's later keys (the Cluster Position, the carrier count, the aligned text)
  only order alignments equal in all of these and are not compared. The merge sorts on the
  CFD of a row's own alignment, written to three decimals (the results' `CFD_score_ALT`;
  their `CFD_score` prints the higher of the REF and ALT scores), so each row's CFD is
  recomputed from its aligned strings and compared at three decimals. A PUTATIVE target is
  the carrier's own: it is held to the alignments under which the rule picks it for that
  carrier, and a target the rule ties with another to its own best. The site is the cluster
  position on the haplotype carrying the record: CRISPRitz's Cluster Position (the start
  moved to the PAM's fixed end, shared by every alignment of one PAM site) plus the rank of
  the row's first base among the bases an insertion adds, as the merge compares rows that
  start inside one insertion; a row whose own target is not at its start stands for no
  other alignment. A target that needs a second indel-path record is demanded, as the
  window's target of the same bases, start and variants, exactly when CRISPRme's reading of
  the slot against that target holds (from its first base on or up to its last, through the
  record: it then keeps the carrier), and is the two-indel limit otherwise. A target
  touching a call of unknown phase is demanded as the sample's strongest (a base on the
  record's own REF span is no choice, as above). A target that an overlap zone of the
  carrier's slot changes on some reading and not on another is counted, not demanded. A
  run that writes no alternative alignments (`.alt_mode` = `no-alt-alignments`: by default
  a population-level run, as the chr22 self-test's are, whatever `--ci-cd-test` says)
  keeps one row per merge cluster: there the target's cluster — the integrated row with a
  block within `--merge` bp of its Cluster Position, or of one that the rows CRISPRme can
  write chain to it, each within `--merge` bp of the next, as the merge clusters: its
  reference rows, its SNP rows (in the mode the run ran, its `.search_mode`:
  population-level the worst case of each alignment over each dataset's SNP alleles that
  CRISPRme's stores hold, per sample only the haplotypes individuals hold) and the targets
  some carrier of an indel-path record holds (the per-carrier resolution gives no row for a
  raw target no carrier holds, nor for one that needs two indel-path records; a slot whose
  overlap §4 reads it without is read so, as a clean slot, one §4 takes off gives none; a
  sites-only dataset's rows are the worst case of each alignment of its window, the
  record's ALT written verbatim and a SNV coded only off its REF span, and a record that a
  genotyped and a sites-only dataset both hold gives the rows of both, whichever record's
  work reads it first), every record the chain reaches read, up to 3 kb either side (a
  chain still going there with no row as good as the target is the checker's own coverage
  gap, `coverage-chain-limit`, not a miss) — must be reported with a fewest-mm+b block of
  no more edits than the target's best alignment by the first sorting criterion (the merge
  picks between a cluster's reference and variant rows on mm+bulges, or the CFD, alone)
  and a highest-CFD block of at least its CFD (each row's recomputed from its own aligned
  strings, the value the merge sorts on). The rows CRISPRme may write or not, which no
  documented rule settles (the one-record reading of an overlap §4 does not cover, a target
  an overlap zone of two other alleles of a slot changes, on any reading of the zone, a
  target of a slot read without its overlap that touches a call of unknown phase, a target
  through an SNV its stores read otherwise, a per-sample run's capped windows' greedy rows,
  which §5's known limit writes on the classic dictionary path only, a per-sample SNP
  target of a slot that also holds an indel touching it, which 2.10.0 writes only where its
  reading of the slot through a repeat keeps the carrier, and a sites-only dataset's row
  through one of its SNVs outside `PASS` and `.`), are followed only when the rows it
  writes reach no row as good as the target, and never one that needs a second indel-path
  record: a cluster reached only through them is counted under their cause's class
  (`limit-capped-window`, the checker's `coverage-snp-row-indel` and
  `limit-filtered-record` for the last three), neither passed nor failed. A record only
  sites-only VCFs hold is carried by nobody: it is counted, not demanded; so is one whose
  FILTER the search skips (`limit-filtered-record`).

The check's inputs must be the search's: where the run's own record of them says otherwise
(its `.Params.txt` budget, PAM and variant datasets — one VCF per dataset must hold each
chromosome checked; a search without bulges records `Genome_idx None` whatever it searched,
and its datasets are then the lines of the `--vcf` list its command line names, read in its
working directory —, a guide its `.guides.txt` does not hold, a `--samples` other than its
`.sampleID.txt`, an `--analysis` other than its `.search_mode`, which records the mode the
search ran whatever chose it), the run is refused (exit 2), as is a search of two genotyped
datasets that share samples, a run whose records contradict each other, and one that does
not record what the check needs (completeness: whether it wrote alternative alignments and,
where it wrote none, its mode, which `--analysis` can give). It models a search that read
each genotyped dataset's genotypes (its genotype store, or classic dictionaries), as every
index CRISPRme installs or builds with `--samplesID` does; an install without them
(`download --what index --no-genotypes`) searches population-level and lists every carrier
of an indel on a SNP+indel row, with or without its SNPs (§4). Such a search records it:
in any mode its `log_verbose.txt`, the pipeline's stdout kept in the result folder, holds
the SNP post-analysis's `Registry-only install for <chrom>: ...` line for the dataset
(`hg38_mega`'s aggregate registry prints that line too, marked `(aggregate info_af
registry`, and is not counted), and one made with `--per-sample` writes `.search_mode`
population-level beside `--per-sample`, which `complete-search` does only there or where
every dataset it searches is sites-only; given a genotyped VCF, the check refuses either
(exit 2). Without that line (a result folder without its `log_verbose.txt`, or a dataset
its SNP post-analysis never read), a run made with `--population-level` or no mode flag
records nothing else that tells it apart (its `.search_mode`, `.command_line.txt` and
`.Params.txt` read as a search that read the store; the warning goes to the terminal only):
such a run is not to be checked.

The SNP-only rows of a per-sample run can be checked report-only (`--snp-rows report`):
their strings, their target and that each named variant lies in it, their listed carriers,
and the phased carriers a row of their site leaves out, each carrier read as §4 reads a
SNP row's carriers — an individual with an SNV call of unknown phase in the target is
listed for every non-empty subset of its variants there (`limit-snp-row-unphased`); the
merge joins the carriers of a site's SNP rows of one variant list (its key for a SNP row
is the Cluster Position and the variants, no alignment), so a carrier may be listed for
another alignment of the site its slot holds; and since 2.10.0 a carrier whose slot of the
row's SNPs holds an indel of the phased store that changes the target is dropped unless its
own bases still read the target, from the first base on or up to the last, through the
row's SNPs (kept too where that indel's call cannot be tied to the slot, PUTATIVE, or where
it sits at a named SNP's position, a phasing artefact), and a slot whose bases cannot be
read (an indel near it whose slot is unknown, two records of one base, a base whose allele
its calls leave unknown) leaves the row. Listing a carrier the rule drops is reported
(`snp-row-indel-carrier`), an MNV the store does not hold is `limit-snp-row-mnv`.
Population-level SNP rows are a documented worst-case union and are out of scope.

It does not check the `CFD_score_REF`/`_ALT`, `CFD_risk_score` and CRISPR-Bulge columns or
which row the CRISPR-Bulge block keeps, `Variant_info_spacer+PAM`, `PAM_creation`, the seed
counts, `Variant_MAF`, `Variant_rsID`, `Other_motifs`, `Not_found_in_REF` or the
annotations, nor a row's own `CFD_score` against its alignment. The panel is the samples
list; CRISPRme's input validator refuses one that leaves out a sample of the VCF, so on a
run it accepts it is every sample of the VCF.

The documented limits — the last record of a split multi-allelic site, two indel-path
records in one target, a PUTATIVE carrier listed for its strongest target, a site within
one target of a contig's end, a record whose FILTER the search skips — and the slots the
VCF leaves without a sequence are listed
with their carriers under their own class and do not fail the check, as do the checker's
own coverage gaps (a carrier the PUTATIVE rule resolves with more than 10 positions within
a target's reach of the record, off its REF span, where its alleles allow two bases, whose
targets touching a call of unknown phase it does not enumerate, `coverage-putative-cap`; a
population-level cluster that runs past the 3 kb either side its chain follows,
`coverage-chain-limit`; a slot whose overlap zones near a record have more readings than
it tries for the targets only such readings hold, which it counts and never demands,
`coverage-overlap-cap`; a per-sample cluster reached only through a SNP target of a slot
with an indel touching it, `coverage-snp-row-indel`; a VCF read without a `.tbi` that holds
a deletion longer than the 20 kb such a read looks back, `coverage-unindexed-long-ref`); a
row on a chromosome no VCF covers, a table line shorter than its header and a window
holding N are reported. A target all of whose carriers are documented limits is demanded
of no one, so a row that lists only them is not missed when it is lost: each completeness
limit line says how many of its carriers a row naming the record lists, and that count
changes then, without failing the check (the documented reading of those slots, "stays
listed" and "the last record", is not modeled into a gating demand).

Measured on 2.10.0 itself (its code over the 2.8.1 image with CRISPRitz 2.10.0; every
index built by `build-index-only` with it, of format 2; each search from a clean
environment, as the HPC route runs it), chr22 of the HGDP + 1000 Genomes panel (4,091
samples), the README guide, NRG, 4 mismatches, 1 DNA and 1 RNA bulge: soundness and
completeness pass in both analysis modes (the two, with the per-sample run's SNP-row
report, about 2.5 min with `--jobs 16`), its limits 8 two-indel targets, 4 records whose
carriers' slots hold overlapping alleles where a target reads and the checker's coverage
gap for one carrier's slot near two records (`coverage-putative-cap`: in the IGL locus, it
holds two ALTs of one position at more than 10 positions there), and no split-site
target, every allele of a split site being coded; so do NRR and the TRAC guide per sample
(19 and 20 two-indel targets, 8 overlapping records each), SaCas9's 21-nt
`CACTTTGGGAGGCCGAGGCGG` with `NNGRRT` per sample at 2 mismatches without bulges (the 2,103
carrier targets at a fake-indel window's first base reported; 22 two-indel targets, 7
overlapping records) and at 4 + 1 + 1 (the targets of 16 records that run past their
fake-indel window, a 27-nt target with a DNA bulge, held by no index and counted for their
2,121 carrier listings, `limit-window-reach`; two-indel targets at 66 records, overlapping
alleles at 22), three repeat-rich guides (VEGFA site 3 `GGTGAGTGAGTGTGTGCGTG`, EMX1
`GAGTCCGAGCAGAAGAAGAA`, VEGFA site 1 `GGGTGGGGGGAGTTTGCTCC`) per sample and
population-level (per sample 8,289 notes where a listed carrier's slot, read as §4 reads
it, holds the target, 1,082 row blocks and 1,457 records with overlapping alleles, 1,284
two-indel targets; population-level 160 notes, 67 row blocks and the same 1,457 records,
271 two-indel targets), the HPRC chr22 panel (9,125 MNVs, 4.5% of its calls with a missing
allele) and a copy of it with every fourth sample unphased, both modes (3 row blocks of an
MNV whose carrier CRISPRme resolves PUTATIVE, as it documents; 6 two-indel targets; 17 and
26 windows across N; the checker's coverage gap for 21 records of the unphased copy,
carriers resolved PUTATIVE with more than 10 positions of two possible bases within a
target's reach), and the default index beside the sites-only `hg38_mega`, both modes
(2,254,806 and 666 row blocks of `hg38_mega` list no carrier and are checked as rows, each
naming SNVs of `hg38_mega`'s VCF whose FILTER passes, per sample 615,803 lines of the
alternative alignments at one (AT)n repeat, chr22:24,238,700-24,238,900, where `hg38_mega`
lists several SNVs and indels at almost every position; completeness counts the 127,735
indel records of chr22 only `hg38_mega` holds and the 32,178 the search skips for their
FILTER; 19 and 9 two-indel targets, 8 overlapping records); each read as a run of 2.10.0
from its own records (every allele, a store of `PASS` and `.`, an INDELS index of format
2), as is the README per-sample search with 2.10.0's code on the indexes of CRISPRitz
2.9.1 (the last record of a split site, a store of every FILTER, which its log notes,
format 1), which passes too, its limits those above and 4 `limit-last-alt-wins` lines. The
genotyped panels of these configurations have no FILTER other than `.` (HGDP + 1000
Genomes, HPRC) or `PASS` (the self-test's 1000 Genomes); the SNP+indel oracle with the
second deletions of its sites O and Q, the insertion of Y1 and an MNV marked `LowQual`,
searched through the real CLI, passes in both modes, and fails where the run read a store
of every FILTER, as one written before 2.10.0 is (14 and 6 gating findings: S21 and S39
left off); the same world searched with `CRISPRME_INDEL_SNP=0` is refused in both modes
(the check of 2.8.2's code gave 166 and 77 gating findings there). On the documented chr22
self-test (`complete-test`, as CI runs it) soundness passes every benchmark, the 5′ Cas12a
PAM included, and completeness of `cas9_sg1617` (population-level, its `.search_mode`)
passes: CI's step, replayed verbatim, passes there (about 4 min on 4 processes). The
SNP+indel oracle passes in both profiles and both modes, with no known defect left to
list. The negative control fails as before: on the released 2.8.1's searches of the README
guide soundness reports exactly the four rows per mode whose first aligned base is one an
insertion adds, and nothing else, and completeness reports their 1,306 carrier targets
missing at the right coordinate per sample (and, population-level, the one such row placed
more than `--merge` bp past its cluster), the run read as one of 2.8.1. Measured with the
check of 2.8.2's code (rel/2.8.2's post-analysis over the released 2.8.1 image, indexes of
format 1), which gives those 2.8.1 searches the same findings: completeness of
`cas12a_hbg` on the self-test passes, removing the 5 indel rows of `cas9_sg1617` is caught
at each of their 5 clusters, and CI's step fails on 2.8.1's self-test; the README guide
and review-r9's ACTGAAATCTGTAAGCAGGC without bulges (4 mismatches, both modes; the
second's 3 indel row blocks) pass. Corrupting rows of the README per-sample run — the
start moved by one base, a REF base, the bulge count, a non-carrier added, a carrier
dropped, the strand flipped, a target base — is caught on 8 of 8 rows of each class in the
integrated table and in the alternative alignments. On the 2.8.0 run soundness reports
10,093 row blocks per sample, among them the alternative alignment at chr22:22322768 that
paired an insertion's only carrier with SNP alleles other individuals carry and two rows
that list a carrier whose slot holds a second deletion inside the target and a SNV on the
base it deletes, which §4 cannot read through and lists on no row, and completeness 1,425
carrier targets, 109 of them met only through a fold into a phantom row, 5 listed only by
rows whose site has no fewest-mm+b block with their carriers, and 2 whose target (5 edits,
3 mismatches) the rows of its site with their carriers hold only as another target of 5
edits and 4 mismatches, which the merge's fewest file, sorting on mm+bulges and then mm,
ranks after it. Genome-wide, soundness of the 2.8.1 runs of the README guide (per sample
and population-level) and of TRAC reports only that inside-insertion class (1,111, 682 and
1,332 row blocks, 70 to 100 s each with `--jobs 24`), and completeness of the README run
reports the same class again, as carrier targets missing at the right coordinate (63,651
carrier targets of 348 records, each a carrier of a 2.8.1 row of the same record placed
past its anchor, and 4 whose target at the anchor lost its highest-CFD block to such a row
at its Cluster Position), and nothing else that gates (about an hour on 40 processes of a
loaded node: a background job, not the release gate).
The SNP-row report finds no listing the documented rules do not explain on 2.10.0's
per-sample runs above, but beside `hg38_mega`, where 144 SNP row blocks (65 lines of the
alternative alignments, at 51 sites of chr22; 60 blocks at 23 sites on 2.8.2's code) name
only an HGDP + 1000 Genomes SNP and list its carriers for a target that holds more variant
bases than that: the merge, whose key for a SNP row is the Cluster Position and the
variant list, keeps `hg38_mega`'s worst-case alignment of the site under them. On 2.8.1
and d7a6fff its 30 row blocks list carriers whose slot of the SNP holds an indel that
changes the target, which 2.10.0 drops. On the runs of 2.8.2's code the checker's CFD of a
row's own alignment was the `CFD_score_ALT` the runs printed on every one of their 46,626
indel row blocks, and on all but 14 of about 2.4 million SNP and REF blocks (each a value
at a x.xxx5 rounding boundary, which the comparison gives the benefit of the doubt); the
same input gives the same findings byte for byte (and did on a fresh run of 2.10.0). It
found three CRISPRme defects that 2.10.0 fixes. The merge keyed an indel row by Cluster
Position, variant list and aligned target, not the aligned spacer, so an unphased carrier
resolved on an alignment where its SNP falls in a DNA bulge (CFD picks the reference base
there) was listed on another alignment of the same aligned target whose column there faces
the guide's base, which the carrier carries — a weaker target than its strongest (the
SNP+indel oracle's site U, S19; and in a synthetic world of partly phased carriers
searched through the real CLI, a carrier phased for the deletion and unphased for the SNP,
which 2.10.0's code passes). A search whose targets reach 27 nt (SaCas9's 21-nt guide and
6-nt PAM) lost, on an INDELS index of format 1, those that start at a fake-indel window's
first base (chr22, `CACTTTGGGAGGCCGAGGCGG` with `NNGRRT`, 2 mismatches, per sample: 2,103
carrier targets at 6 records, 2,088 of them the insertion chr22:17325232 A>AT's; an SpCas9
target, 24 nt at most here, never starts there holding its record), which an index of
format 2 holds: the same search of indexes 2.10.0 built reports them all, its row at
chr22:17325206 (+) listing the insertion's 2,088 carriers. And the phased indel store kept
every indel whatever its FILTER, so a carrier was taken off a target, and off a SNP row,
for an indel no fake-indel window of the search holds (the oracle's sites O, Q and Y1
above). It also shows a reading no text documents: `hg38_mega`'s SNP rows name SNVs whose
FILTER is not `PASS` (on the published index also earlier records of split sites: 58 split
sites of chr22, 57 of them a `PASS` record before a later `PASS` one), where its indel
rows and indel log follow CRISPRitz's enricher. The variant genome CRISPRitz 2.10.0 builds
from that VCF holds no such SNV, but the registry lists every ALT of a position and a
sites-only dataset's SNP row substitutes each: the fresh per-sample search beside it names
one of the 74,025 SNVs of chr22 outside `PASS` that no passing record of either VCF holds
in 1 integrated row and 3,434 alternative alignments (chr22:12810728 (−),
chr22_12810729_T_C, `ExcHet`, beside the `PASS` T>A the genome codes), the
population-level one in none.
`test/indel_snp_oracle/run_oracle.sh` runs it beside `oracle.py check` in both profiles,
`validate-benchmarks.yml` runs it on the chr22 self-test output after `validate-test`, and
`docs/RELEASE_CHECKLIST.md` makes it a release gate.

This establishes that the engine **does not miss** off-targets relative to
exhaustive search. It does **not** validate the *scoring* models' predictive
accuracy against experimental cleavage — a regulatory-grade package should add a
retrospective comparison of CFD/CRISPR-Bulge ranking to experimental off-target
assays (e.g. GUIDE-seq / CIRCLE-seq / targeted amplicon or rhAMP-Seq).

---

## 9. PAM-less indexes and dynamic PAM filtering

A CRISPRitz index is built for one PAM. The PAM is matched against the genome at build
time, only the surviving sites enter the ternary search tree, and the PAM region itself is
stored as a flank beside each guide. An index therefore answers its own PAM and, pruned
(below), a PAM its own holds, and no other: an NRG index serves an NGG search, but an NGG
index cannot serve an NRG one, so searching NRG and NGA over the same genome means building
and storing two dedicated indexes, or one NRR index, or a PAM-less one.

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
therefore takes any index only for a PAM on the side it was built for: the side the header of
its `.bin` files records, for an index of format 2 (below), else the one recorded in its
`.pam_build` (written by `build-index-only`, and by a search that builds an index). A dedicated
index of format 1 without one is used
as before. A PAM-less index without one is taken for a 3′ build at the published widths, 3 and 8
(the K = 8 archive has no `.pam_build`), or under the search PAM's own name, and at no other
width: a search before 2.8.1 built one of any width from a 5′ PAM file and recorded nothing. A
search, and `build-index-only`, build an index under a name that says its side (since 2.10.0;
**Index names** below), so an index of the older name for the other side is left beside it, and
never searched for this side; `build-index-only` does not record a PAM file over an index
recorded for the other side, nor a 3′ one over an unrecorded PAM-less index of another width.
The search PAM is
widened to *K* with `N`s on the distal side (the end for a 3′ PAM, the front for a 5′ PAM),
and a leaf is rejected when, at any position the real PAM covers, the stored and requested
IUPAC codes share no nucleotide. Pad positions are not tested at all: they state no
requirement, and a stored genomic `N` intersects nothing, so testing one would reject a
target that a dedicated index of the real width keeps. At *L = K*, on the same side, the
pad is empty and the geometry is exactly a dedicated index's.

**A dedicated index of a broader PAM.** A dedicated index stores, as its flank, the concrete
bases of every site of its own PAM, so it holds every site of a PAM its own holds at each
position: an NRR index those of NGG, NAG, NGA, NAA and NRG; an NRG index those of NGG and
NAG. `--prune-pam` narrows it to the search PAM the same way (the two are as long, so there is
no pad), and since 2.10.0 CRISPRme takes such an index for a search, under `--prune-pam`, as
it takes a PAM-less one: the search PAM's own index first, then a dedicated one of a PAM that
holds it (of two, the narrower PAM, NRG before NRR for NGG), then a PAM-less one. Of the
PAM-less ones, the narrowest flank, as before. The side rule above applies to it as to a
PAM-less index: a record, or else a 3′ build at the published widths.

**Index names.** An index was named `<PAM>_<N>_<genome>[+<dataset>]`, which says neither the
guide length it was built for nor the side of its PAM, so a 3′ and a 5′ build of one motif, or
builds for two guide lengths, took one name, and the second stopped on the first. Since 2.10.0
an index is built under a name that carries both, in the order of the target on the strand:
`<guide>-<PAM>_<N>_…` for a 3′ PAM (`20-NRR_3_hg38`), `<PAM>-<guide>_<N>_…` for a 5′ one
(`TTTV-23_3_hg38`), its INDELS index `<name>_INDELS`. `build-index-only` and a search name what
they build so (the engine writes the older name, in a build folder, and the index is renamed
into place). The name records the side as a `.pam_build` does, ahead of one, so such an index
needs none: CRISPRme takes it only for a search of the side its name says, at any flank width
the rules above allow on that side. It reads the guide length from the trees as before, and of
one index under both names takes the newer. An index under the older name is read as before:
`build-index-only` keeps one in place that records this side and this guide length (or no side,
by the rules above), and builds beside one that records another, whose refusal of 2.8.1 said
to build in another working directory; one whose side its record leaves open (an all-`N` one
of a width but 3 and 8, for a 3′ PAM file) is refused, with the name of each side that
records it.

**Pruning by itself.** Since 2.10.0 a search prunes, as if `--prune-pam` was given, whenever
no installed index serves it unpruned and one of these does: a wider all-`N` one, or a
dedicated one of a PAM that holds the search's. It says so first, naming the index that needs
it ("no installed reference index of hg38 serves PAM NGG unpruned, and NRR_3_hg38 serves it
pruned to it at leaf time: searching with --prune-pam, as if it was given"): a search from the
command line on its console, one from the web, which has none, in its `log_verbose.txt`. The
run's record of the indexes it searched (`.indexes_searched.txt` in the result folder, beside
the parameters: `Reference_index`, `Variant_index:<dataset>`, `Prune_pam auto`;
`indexes_searched.txt` in the report's zip) and the report's summary say which index served
it. Before, the search stopped for want of a variant index (and said to pass the flag), or
built a reference index of its PAM from the genome: hours on hg38. Where an index serves the
search unpruned (the PAM's own, or an all-`N` one of its width), the search runs as before,
unpruned; pass `--prune-pam` to prune it anyway. A search with no bulges opens no reference
index, so it neither prunes nor builds one (below).

**Equivalence.** Pruning changes only which leaves are written out, never how a guide is
aligned to the genome, so it cannot move an alignment or a score. Measured genome-wide
over HGDP + 1000 Genomes (4,091 samples, guide `ACTGAAATCTGTAAGCAGGC`, NRG, 4 mismatches +
1 DNA bulge + 1 RNA bulge, per-sample): a 3-nt-flank PAM-less index with `--prune-pam`, an
8-nt-flank PAM-less index with `--prune-pam`, and a dedicated NRG index produce
**byte-identical** reports: the integrated results are identical, 182,369 rows each (MD5
`c76f3490` once the licence-gated COSMIC column, empty on the dedicated arm's fresh install, is
set aside; the two PAM-less arms are identical with it, MD5 `c4586f7a`). The same holds for
NRR: 386,012 rows each (MD5 `f810dd76`) over a dedicated index, both PAM-less indexes with
`--prune-pam`, and the 3-nt one without it. Those indexes are of format 1, built by CRISPRitz
2.9.1 or earlier (below). For them it holds at *L = K*, and at *K > L* for a guide
with no site at a contig's very end or start, as this one has none. The build scans a contig
for a PAM only where a whole window of *K* bases fits, so with *K > L* a + strand PAM with 1 to
*K − L* bases after it to the contig's end, and a − strand PAM in the contig's first *K − L*
bases, are not indexed. On hg38 at *K* = 8 that is 145 + 101 NGG and 308 + 241 NRG
sites (about 7.5 × 10⁻⁷ of all), none on a primary chromosome for those PAMs; in two
whole-genome searches the 8-nt flank lost one row each against a dedicated index (27,215 of
27,216 rows for a guide picked for its perfect match on `chr11_JH159136v1_alt`, and 216,977 of
216,978 for VEGFA site 3 at 5 mismatches). A − strand PAM exactly *K − L* bases from the start
is indexed in the window at the contig's first base, which CRISPRitz 2.9.1's search read as a
window of the + strand, as it read that of a − strand PAM at a contig's first base (below): it
reported that PAM on the + strand at another coordinate. CRISPRitz 2.10.0's search reads it
right, in an index of either format. A search that reads a flank of format 1 wider than its PAM
says so: from the command line on stdout and in the run's `Params.txt`, from the web in its
`log_verbose.txt`; under CRISPRitz 2.9.1 it also says that PAM is reported on the + strand. N
runs inside a contig are not affected by the width. An index of
format 2 stores a PAM-less window whose flank runs off the contig, its
bases off the contig as `N`, which matches no PAM base: with `--prune-pam` it returns the rows of
a dedicated index of format 2 **built for the same guide length and budget**, for any 3′ PAM of
*L ≤ K*, contig ends included (CRISPRitz's `docs/INDEX_FORMAT.md`, section 8; through the real CLI
on a synthetic genome of sites at contig ends and N runs, `test/index_format2/check_boundary_rows.py`,
an NGG search off an 8-nt flank gives the dedicated index's 13 rows, where of format 1 it gives
5 of 6). Of two indexes of format 2 built for one guide length, the one of the larger budget
returns a superset of the other's rows. Neither holds across guide lengths (below).

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
A search that reads such an index wider than its PAM warns of it beside its note of the flank
(it reads the smallest tree of a contig whose sequence holds an `N`, and of the INDELS index
one of several indel windows), and goes on.

**Contig ends and N runs.** CRISPRitz writes an index in one of two formats (its
`docs/INDEX_FORMAT.md`). Up to 2.9.1 it writes **format 1**, which every index published until
the rebuild for CRISPRme 2.10.0 is; from 2.10.0, the CRISPRitz CRISPRme 2.10.0 installs, it writes
**format 2**, whose `.bin` files open with a header that records the format, the build PAM, its
side, the guide length and the budget. The limits of this paragraph are those of format 1, for
every index, dedicated or PAM-less, at any width. An index of format 1 of budget *N* holds a window of the guide, the PAM and *N* more bases, and the
build drops a window that holds an `N` or runs past the contig's start, so it cannot report a
site whose guide lies within *N − 1* bases of a contig start, a contig end or an `N` run, even a
perfect match; a search without bulges, which scans the FASTA, reports it (and, since 2.10.0,
neither looks for nor builds a reference index, nor a variant index where the whole enriched
genome is there to scan: it built `<PAM>_1_<genome>`, which nothing read). A variant search
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
boundary: its report then holds no perfect match. Two more limits of CRISPRitz 2.9.1 at contig
ends, as measured on hg38: a − strand PAM at a contig's first base is reported on the + strand
at base *N* (11 NGG and 37 NRG contigs open so, all alt, random or chrUn contigs; the primary
chromosomes open with `N`), and a + strand PAM on a contig's last bases is found by no search,
indexed or not (8 NGG and 31 NRG contigs end so, none a primary chromosome, and no fake-indel
contig, which all end in `N`). CRISPRitz 2.10.0 lifts all three. Its search reads a − strand PAM
at a contig's first base on the − strand, in an index of either format; its scan reaches a PAM on
a contig's last base, also without an index; and an index of format 2 keeps every site whose guide
lies on the contig with no `N` in it, the pad of *N* bases being allowed to run off the contig or
into an `N`. An alignment is then reported only where every base it uses exists, so every row is
the contig's own bases. On hg38 chr22 and chr1 a format-2 NGG index returns every row of the
format-1 one and 1,827 and 5,512 more, all next to an `N` run (CRISPRitz's `docs/INDEX_FORMAT.md`,
section 11). An index of format 1 keeps the first and the last limit until it is built again
(a search without bulges, which scans the FASTA, finds the last PAM of a contig at once), and a
search of one by CRISPRitz 2.10.0 says so in its `log_verbose.txt`. CRISPRme 2.8.1 refuses an index
of format 2, and 2.8.0 aborts on one; CRISPRme 2.10.0 searches it where the installed CRISPRitz
reads it (2.10.0 and later: the `INDEX_FORMAT` of its `crispritz.py`), and refuses it, naming
both, under CRISPRitz 2.9.1. Which sites a contig holds, which CRISPRme checks to tell a lost
`.bin` from a contig with no site, follows the format the index was built in, read from its own
`.bin` files; an index CRISPRme builds records it, with the build PAM, side, flank and guide
length, in its `.index_meta.json` (a sub-index with no `.bin` has no other record of it).

**Cost.** A PAM-less index is much larger, because no site is excluded at build time: the
reference index is 71 GB against 9 GB for a dedicated NRG index, and the HGDP + 1000
Genomes variant index 110 GB against 28 GB (download sizes, GB = 10⁹ bytes). Pruning
happens at leaf time, so a PAM-less search is not faster than a dedicated one — the gain
is that one index answers many PAMs.

**Guide length.** An index name of the older form carries no guide length (one of the 2.10.0
form says it, **Index names** above); each tree records the one it was built for, as the second
int32 of its `.bin` in format 1 (written by `buildTST` since 2019, whoever built the index),
and in its header in format 2, which is what CRISPRme reads. The engine reads that length only
to offset coordinates and searches any guide it
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
that merely fits. With none that fits, a search builds its own under a name that says its guide
length (`<guide>-<PAM>_<N>_<genome>`, **Index names** above), beside an index of the older name
`<PAM>_<N>_<genome>` built for another guide, which it stopped on before 2.10.0, as
`build-index-only` did; both leave that index and its `.pam_build` as they are. One of the new
name whose trees were built for another guide than its name says stops either, with an ERROR
that names both lengths and the name that records what it is.
A longer build serves a shorter guide, but, like any deeper index of format 1, it does not hold a
site next to an N run or a contig end that a shallower one does. The rule was measured on
indexes of format 1; on one of format 2 the rows inside it are a dedicated index's away from
contig ends and N runs, but not there: an index built for 20-nt guides holds sites whose 21st
base would lie off the contig or on an `N`, which an index built for 21-nt guides does not, and,
with bulges, a few alignments differ at interior sites too (CRISPRitz's `docs/INDEX_FORMAT.md`,
section 5.1, measures both). So every equivalence above holds between indexes built for the same
guide length; for the rows of a guide length, search an index built for it, which the search
prefers.

**Limits.** The web interface never passes `--prune-pam` itself; its searches prune by
themselves as above. (Before 2.10.0 the flag also required `--index-path`.) (Its PAM list
prefers a PAM whose geometry matches the installed index, but falls back to offering every PAM
when none matches, so it does not itself prevent an unusable pairing; the PAM picker of a
PAM-less index offers only the PAMs that index serves. The search takes a wider all-`N` index
for a 3′ PAM only pruned, which it then does by itself, and never one built for the other side,
and when no index serves it and neither the dataset's VCFs nor its enriched genome are there to
build one, it stops with an ERROR that names what was asked, the variant indexes installed and
what to do.) A shorter PAM
off a wider flank is supported for a 3′ PAM only; for a 5′ PAM the engine requires the flank width to equal the
PAM length, and the index must have been built from a 5′ all-`N` PAM file, as its name,
the header of its `.bin` files or its `.pam_build` records (above). The published
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
