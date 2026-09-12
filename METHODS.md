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

CRISPRme+ therefore ships the registry **uncompressed (raw)** by default: the three
sections are memory-mapped and read by direct `O(log n)` bisection on genomic
position, with **no per-lookup decompression**. This is the fastest option — a dense
variant search does millions of registry lookups (one per candidate off-target window,
scattered across a chromosome), and measurement showed raw lookups are **~2× faster**
than the compressed path even with a warm cache. All shipped indexes (the genotyped
1000 Genomes + HGDP and the sites-only mega) use this raw format, for consistency and
lookup speed ("speed over space").

An optional **`zlib` block-compressed** format is available for space-constrained
storage (opt in with `CRISPRME_REGISTRY_COMPRESS=1` at build time, or re-encode an
existing raw registry losslessly with `transcode_registry` — no VCF re-parse): each
section is partitioned into 4,096-record blocks, individually compressed (~3.5× smaller
on disk), with a sparse block index in the manifest. The reader is fully
backward-compatible and reads either format byte-for-byte identically. When a
compressed registry is used, a decompressed-block **LRU cache** (default 4,096 blocks,
`CRISPRME_REGISTRY_CACHE_BLOCKS`; each block ~64 KB, filling only to the blocks a search
touches) amortizes decompression to one pass per block: profiling had shown block
`zlib` decompression reaching **~71 %** of a dense post-analysis with a tiny cache, and
the right-sized cache plus an O(1) LRU cut a dense chr22 SNP post-analysis **326 s → 43 s
(~7.6×)** — but the raw default avoids that cost entirely, which is why 2.5.3 ships raw.

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

### Two panel modes (genotyped vs sites-only)
The merge above is the **genotyped, cis-capable** mode (`merge_vcf_panels.sh` /
`build_combined_panel.sh`): it keeps per-sample genotypes, recomputes a **pooled**
`INFO/AF` (AC/AN over the union of genotyped samples, `bcftools +fill-tags`), and
feeds the genotype-counting Tier-0 registry. It is correct only when *every* merged
source is genotyped, and it is what enables indel+SNP cis reconstruction.

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
homozygote counts do not exist for aggregate data, so they are reported as
Hardy–Weinberg expectations from the allele frequency (the frequency itself is exact;
carrier/hom are flagged as estimates). The registry stores SNPs only (single-base
ref/alt); indels (~24 % of merged sites) are surfaced as off-targets through the
fake-indel genome but are not yet frequency-annotated in this mode.

The two modes are therefore complementary: the genotyped panel gives phased,
cis-capable frequencies over a curated sample set; the mega gives one-scan,
frequency-annotated coverage across the widest set of population resources.

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

### Population-level analysis (default) and `--per-sample` genotype resolution

The controls above bound any *single* window, but a **dense panel** (many merged
sources) or a **sites-only aggregate panel** can present so many variant-dense windows
that even the observed-haplotype enumeration of Section 4 becomes intractable — measured
at **49 h+ without completing** on a 4×-density 1000G+HGDP panel. The two analysis modes
are therefore not really a *speed* dial but a **genotype-resolution** dial: the default
reports population-level worst-possible off-targets (works on any panel), while
`--per-sample` resolves each individual's observed haplotype (only possible, and only
meaningful, when the panel carries per-sample genotypes). So as of 2.5.4 CRISPRme+ runs a
**two-pass population-level analysis by default** (propagated to the whole post-analysis via
the internal `CRISPRME_FAST_MODE` env var); **`--per-sample`** opts into the exact
observed-haplotype enumeration of Section 4. Instead of enumerating the 2ᵏ IUPAC haplotype
lattice per window, the population-level analysis emits a small fixed set of
**worst-possible representatives** per window:

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

The min-edit and max-CFD representatives are distinct because CFD is position-weighted
(Section 8): the fewest-mismatch haplotype is often **not** the highest-scoring one, so
both are emitted so the reported worst case is never understated. The population-level
analysis is validated to be **lossless for locus detection and non-understating for the
worst-case score** against the exact per-sample enumeration path on a real chr22
1000G-2021+HGDP slice (0 CFD under-reports; it in fact surfaces *stronger* worst cases at
182 loci that per-sample enumeration misses), and it collapses ~1.9× fewer rows on that 1×
slice, growing with density — turning the otherwise-intractable 4× panel into a tractable
run. **The trade-off is that per-sample carriers, CONFIRMED cis phasing and exact joint
allele frequency are not computed in the default mode** (rows are PUTATIVE worst cases); the
launch-time message, the web **Analysis mode** control and the report's **Search mode** row
all state this explicitly, so it is never silent. This yields a **two-tier workflow**: the
default population-level analysis for routine, high-density, or aggregate-panel *screening*
(and the only meaningful mode on sites-only panels), and `--per-sample` (byte-identical to
the pre-2.5.3 enumeration path) for *confirmatory / pre-IND* runs on a genotyped panel where
per-sample phased haplotype resolution is required. On sites-only (aggregate) panels there
are no per-sample genotypes to recover — verified on the released **mega** index:
`--per-sample` yields **zero CONFIRMED rows and no cis carriers** in either mode (there is no
genotype store to enumerate), and the `indel_snp_cooc.tsv` / `indel_af.tsv` companions are
**byte-identical** across modes. `--per-sample`'s SNP path in fact emits *fewer* worst-possible
windows (integrated 3,065 vs 3,224 rows; SNP+SNP co-occurrence 4 vs 22 PUTATIVE rows), so the
population-level analysis — with its conservative worst-possible over-listing — is the more
complete screen there. Accordingly `--per-sample` is a **no-op (with a warning) on a sites-only
index**, and the web form **disables** the per-sample option when a sites-only index is selected.

**Measured population-level vs `--per-sample` behavior (adversarially verified).** Running both
modes on the released 1000G-2021 + HGDP genotyped index for one guide (identical
guide/PAM/thresholds; only `--per-sample` differs) makes the trade concrete and confirms the
guarantees hold end-to-end:

- *Detection is lossless at the **window** level.* Every off-target window `--per-sample`
  reports is present in the default analysis; the default in fact reports **more** loci
  (3,287 vs 3,103) because it emits worst-possible representatives. A small number of windows
  are anchored a few bases apart between the two modes — `--per-sample` anchors each *observed*
  haplotype, the default anchors its worst-possible representative — so a strict
  coordinate-equality comparison flags them as "per-sample-only" (24 loci here, all 1–10 bp
  shifts of a window the default did detect). Compared by cluster they match exactly and the
  default is a strict superset. Losslessness is therefore a **window/cluster** property, not
  byte-exact coordinate identity.
- *The worst-case CFD bound holds exactly.* Across all 3,079 shared loci, the population-level
  CFD is **≥** the per-sample CFD with **zero** violations; it is strictly higher only where it
  assumes worst-possible co-occurrence (e.g. CFD 0.354 vs 0.113 at a 3-SNP window), never lower.
- *`--per-sample` corrects — and **tightens** — the carrier sets.* The default lists a PUTATIVE
  **union** of every sample carrying *any* contributing variant in the worst-possible window;
  `--per-sample` prunes to the actually-observed cis haplotype's carriers. Counterintuitively
  `--per-sample` therefore names **fewer** samples (1,447 vs 1,642 distinct here), not more — the
  default list deliberately over-includes non-cis carriers (a conservative screening bias). The
  SNP+SNP co-occurrence companion is all-PUTATIVE by default (0 CONFIRMED) and gains **106
  CONFIRMED phased-cis rows** with exact joint AF and single named carriers under `--per-sample`
  (e.g. a locus that is a 175-carrier PUTATIVE union by default resolves to a single confirmed
  carrier under `--per-sample`).
- *On a **genotyped** index the default does not blank carriers.* The `NA`-carrier behavior is
  specific to sites-only panels; a genotyped panel still surfaces dataset-level carriers in the
  default analysis, so the real trade is **PUTATIVE-union vs CONFIRMED-exact-cis**, not
  named-vs-`NA`.
- *For a single guide the two modes take essentially the same wall time* (~13–14 min) — the
  modes differ in genotype resolution, not intrinsic speed; the default's advantage materializes
  only on the dense/aggregate-panel enumeration wall, not sparse single-guide runs.

**Scope of the analysis mode.** The mode toggles only the **SNP** post-analysis (`--per-sample`
re-enables the 2^k IUPAC haplotype lattice / observed-haplotype enumeration that the default
collapses). The **indel** post-analysis is single-threaded and CRISTA-scoring-bound, and is
**unaffected by the mode** — a dense indel search pays the full indel cost regardless
(parallelizing that path is a follow-up). Correspondingly, the `indel_snp_cooc.tsv` companion is
**byte-identical** in both modes (§8).

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

**Cancer-gene annotations: IntOGen by default, COSMIC by licence.** For cancer-gene
context CRISPRme+ ships **IntOGen (CC0)** enabled by default (`intogen_drivers.hg38.bed.gz`,
built by `seq_script/build_intogen_annotation.py` mapping the CC0 driver symbols to
GENCODE gene spans) — legally clean for academic and commercial use. **COSMIC** (Genome Research Ltd /
Wellcome Sanger) is free for academic / non-commercial research but its
**commercial use requires a licence** ([terms](https://www.cosmickb.org/terms/)).
Because the built-in bundle bakes COSMIC into a single BED, CRISPRme+ **strips the
COSMIC rows from the active annotation by default**, so no `Annotation_COSMIC`
values appear in the output. The user opts in **once** — a Settings checkbox on the
web (with a link to the terms) or `crisprme.py cosmic-license enable` on the CLI
(both persist to `Annotations/.cosmic_license.json`); the CLI `enable` requires an
explicit confirmation. Toggling it changes the active-annotation cache signature, so
the next search rebuilds the annotation with or without COSMIC accordingly.

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
read as "present, frequency effectively 0", not a measured frequency. For an
**aggregate (sites-only) panel** (e.g. the all-source mega index, §2), only allele
frequencies exist — there are no per-individual genotypes — so carrier and homozygote
frequencies are rendered **NA** rather than fabricated, and a per-dataset **`indel_af.tsv`**
companion carries indel allele frequencies by source alongside the SNP+indel co-occurrence
table.

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

**Worst-case scoring in the population-level analysis.** In the default population-level
analysis (§5), each window is represented by worst-possible rows rather than every haplotype,
so the *scores* attached to those rows are defined as worst cases over the window's allele
combinations. **CFD is the exact worst case.** CFD factorizes as a product of per-position
maxima times a **joint two-base PAM factor**, so the maximum over all combinations is found by
a per-position argmax plus a bounded brute-force over the PAM region (the joint factor is why a
naïve per-column greedy is insufficient). This exact maximizer is validated bit-for-bit against
an independent factorized oracle and against the exact per-sample enumeration path — **zero CFD
under-reports** on a real chr22 1000G+HGDP slice (and, on the legacy dict / aggregate-panel
path, catching cases where the fewest-mismatch allele scores materially *lower* CFD than
another carried allele, up to a threshold-crossing 0.14). A whole-index comparison of the two
modes on the released 1000G-2021 + HGDP panel reproduces this end-to-end: across all **3,079
shared off-target loci, the population-level CFD ≥ the per-sample CFD with zero violations** (it
is strictly higher only where it assumes worst-possible co-occurrence). **CRISTA is
best-effort.** CRISTA is a non-factorizable RandomForest, so its worst case is taken as the
maximum over the emitted representatives rather than an exhaustive per-haplotype search.
Measured against the exact path (chr22 1000G-2021+HGDP), this approximation is tight exactly
where decisions are made: **every off-target with CRISTA ≥ 0.2 is reported at full or greater
strength** (the population-level analysis even surfaces *more* actionable sites than per-sample
enumeration), and under-reporting is **bounded to ≤ 0.04 and confined to the sub-0.19 weak
tail** (median gap 0.006, no threshold crossings) — structurally, because high-CRISTA
off-targets are low-edit and the min-edit + max-CFD representatives already span the low-edit
shell. **At genome-wide scale the CRISTA tail is heavier than the chr22 slice:** across the full
genome ~5 % of CRISTA ≥ 0.2 loci can drop below 0.2 in the default analysis (largest observed
gap ~0.12), whereas **CFD had zero ≥ 0.2 losses**. So by default CFD is a safe actionable gate
but **CRISTA is a screen**, not an action gate. A **guaranteed per-haplotype CRISTA worst case**
is available by running `--per-sample`; this is the screening-vs-confirmatory two-tier split of
Section 5.

**SNP+indel co-occurrence is unaffected by the analysis mode.** The default analysis collapses
only the *SNP* worst-possible representative emission in `integrated_results.tsv` (Section 5);
the SNP+indel co-occurrence companion (`indel_snp_cooc.tsv`) is produced by the indel
post-analysis' cis phasing pass over the genotype tiers, which the mode does not touch. Measured
on the complete genome-wide matrix (2021 panel, same guide, default vs `--per-sample`): the two
`indel_snp_cooc.tsv` files are **byte-identical** (same MD5, 2,729 rows, 843 CONFIRMED / 1,886 PUTATIVE, full
per-sample `cis_samples` and joint-AF in both). So per-sample cis attribution — which individual
carries the indel and SNP together — is preserved identically in both analysis modes.

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
tagged CONFIRMED-cis (phased) / PUTATIVE (unphased) with per-sample carriers and joint
allele frequency. (The pre-2.5.0 behavior — two independent passes, SNPs on the
IUPAC-enriched genome and indels on a plain-reference fake-indel genome — remains available
by disabling the integration; classic dict builds are byte-identical.) One **residual**:
the indel search materializes **one indel per fake contig**, so an off-target requiring
**≥ 2 co-occurring cis indels within a single protospacer** is not generated as a candidate
— a pre-existing single-indel-search property, independent of the analysis mode. This is a
**low-frequency** case: raw multi-indel cis co-occurrence is dominated by STR/VNTR repeats
(which off-target analysis should soft-mask), falling to **~1–2%** of indel loci after
repeat-masking and deduplication, and the **genuinely-missed** off-targets are **~0.1–0.2%**
of indel off-targets — all at the edit-budget ceiling (the weakest, ≈0-CFD tier). Windows
carrying ≥ 2 cis indels **can** be flag-all'd for conservative (lossless over-reporting)
treatment — a design option the min-edit primitive supports (proven lossless in
`test_twopass_lynchpin_counterexamples`), **not yet wired into the production indel search**. The `Max_total_edits` value is a
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
- The **two-pass population-level analysis** (§5) was validated against the exact per-sample
  enumeration path on a real chr22 1000G-2021+HGDP slice: **lossless locus detection** and a
  **non-understating worst-case bound** (0 CFD under-reports; 182 loci where it surfaces a
  *stronger* worst case). Its exact worst-case-CFD maximizer is additionally cross-checked on
  4,000 random windows against an independent factorized CFD oracle (agreement to the raw
  double, including the joint-PAM case), and on the legacy dict / aggregate-panel path against a
  real CFD-scored multiallelic fixture; the CRISTA best-effort bound is the measurement
  reported under *Worst-case scoring in the population-level analysis* above.

This establishes that the engine **does not miss** off-targets relative to
exhaustive search. It does **not** validate the *scoring* models' predictive
accuracy against experimental cleavage — a regulatory-grade package should add a
retrospective comparison of CFD/CRISTA ranking to experimental off-target assays
(e.g. GUIDE-seq / CIRCLE-seq / targeted amplicon or rhAMP-Seq).

---

*Software: CRISPRme+ (`pinellolab/crisprme-plus`). This document tracks the
methods as of the 2.5.x line (default SNP+indel co-occurrence, two-pass population-level
analysis with opt-in `--per-sample` genotype resolution); see the CHANGELOG and the referenced
source files for implementation detail.*
