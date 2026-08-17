# Design: dictless variant post-analysis (v2)

> Status: **design for review (v2)** — not yet implemented. Supersedes the v1 tiered
> sketch. v2 incorporates the requirements worked out 2026-08-17: off-target allele
> frequency is a **SNP-combination (haplotype) property**, reporting must be
> **global (deduplicated) + per-database + per-subpopulation**, both as **allele**
> and **carrier (individual)** frequency, and multi-database merges must not blow up
> disk. Implementation is gated on the golden-output diff in §12.

---

## 1. Problem

Variant post-analysis (`new_simple_analysis.py`) reads per-chromosome SNP dictionaries
`Dictionaries/dictionaries_<vcf>/my_dict_<chrom>.json` mapping
`"<chrom>,<pos>" -> $-joined "<samples>;<ref,alt>;<rsID>;<AF>"` where `<samples>` is a
per-sample `sampleID:genotype` list. For 1000G+HGDP these total **~152 GB** (gzipped
~31 GB since alpha.27) and need ~3–26 GB RAM/worker. They are the reason the batteries
variant index is large and post-analysis is slow. The alpha.27 gzip change bought space,
**not speed** (still a full per-sample parse, now with decompress overhead).

Measured on a real shipped chr22 dict (1.55 GB raw, 219 MB gzip, 993,880 positions):

| Quantity | Value |
|---|---|
| Per-position metadata (pos, ref, alt, rsID, AF), gzip | **4.0 MB = 0.26%** |
| Per-sample genotype lists | **99.74%** |

## 2. What the user actually needs

For each reported off-target, its **allele frequency** *and* **carrier (individual)
frequency**, in three coherent views:

- **Global** — across **all individuals in all loaded databases, deduplicated** (each
  unique individual counted once).
- **Per database** — computed within each database's own individuals (1000G, HGDP, …).
- **Per subpopulation** — within a database, using that database's **native** subpop
  labels (1000G: AFR/AMR/EAS/EUR/SAS; HGDP: AFRICA/AMERICA/CSA/EAST_ASIA/EUROPE/
  MIDDLE_EAST/OCEANIA; future gnomAD/UKB: their own).

Plus, out of the box: off-target detection, PAM-creation, and the rsID / per-position MAF
columns. Per-**individual** carrier lists (personal cards) remain available but are the
"nuanced, on-demand" layer, not the headline.

## 3. The correctness crux — off-target frequency is a COMBINATION property

An off-target is created/modified when a candidate genomic window matches the guide well
enough *given the alleles an individual carries there*. Its frequency is therefore the
frequency of the **specific combination of alleles** (the haplotype) that produces the
off-target — **not** any single SNP's AF, and **not** a product of per-SNP AFs (the SNPs
can be in linkage disequilibrium).

Two consequences drive the whole design:

1. **Off-targets are guide-dependent** → unknown at index-build time → we **cannot**
   precompute per-off-target frequencies into the shipped index. They must be computed at
   **post-analysis**, for the specific targets a search finds.
2. The efficiency win is therefore **not** "precompute everything into a tiny table." It
   is: compute the combination frequency **at query time but SPARSELY** — only at the
   handful of positions the found off-targets touch (dozens–thousands) — instead of
   streaming 152 GB across the whole genome.

**Single- vs multi-variant off-targets** (the load-bearing split):

- **Single-variant** off-target (one SNP flips one base — the common case): its frequency
  **is** that position's allele/carrier frequency. Fully derivable from a per-position
  table (Tier 0), per database and per subpopulation, exactly.
- **Multi-variant** off-target (2+ SNPs must co-occur): frequency = co-occurrence of that
  exact combination among individuals → needs **per-sample genotypes at those positions**.
  Not derivable from marginal per-position AFs.

**Phasing / cis.** A multi-SNP off-target on a single DNA molecule requires the alts in
**cis** (same haplotype). Phased panels (1000G) → exact haplotype counting; unphased →
either assume-cis (upper bound) or report genotype co-occurrence. We carry a per-database
`phased` flag and state the mode used in the output. (This matches what
`iupac_decomposition` / the haplotype loop already grapple with today.)

## 4. Architecture — Tier 0 registry + sparse per-database genotype tier

### Tier 0 — position registry (ALWAYS shipped, tiny)

`Dictionaries/registry_<vcf>/reg_<chrom>.bin` (+ `.idx`): position-sorted binary of, per
variant position:

- `pos, ref, alt, rsID`
- `AF_global` over the **deduplicated union** of all individuals in the databases bundled
  in this index (precomputed at build time; see §5/§6 for dedup)
- `AF[database][subpopulation]` — a **sparse** map (skip zero-carrier groups), database and
  subpop as small integer codes resolved via the manifest taxonomy
- optional `AC/AN` alongside each AF so downstream can recompute or combine safely

Loaded fully into RAM per chromosome (MBs). Powers: (a) off-target detection ("is there a
variant here"), (b) PAM-creation (ref vs alt), (c) the rsID/MAF columns, and (d) **exact
allele + carrier frequency for single-variant off-targets**, in all three views. No
per-sample data. A `manifest.json` carries `{databases:[{name, sample_count, phased,
subpops:[...] }], chroms, total_deduped_individuals}`.

### Genotype tier — per-database, position-indexed, sparse, compact

`Dictionaries/genotypes_<db>/gt_<chrom>.<ext>` (+ position index), **one block per
database** (NOT a merged matrix). Per variant position, store only what's needed to count
combinations by group:

- the **alt-carrier set** as sample indices (sparse — most samples are ref), or a 2-bit
  genotype packing when a locus is dense; whichever is smaller per record
- for phased data, per-haplotype carrier bits (to count cis combinations)
- genotypes keyed by a **canonical sample ID** so the same individual across databases can
  be deduplicated for the global view

Random-access by position (block-compressed + offset index), read **only at the found
off-target positions**. So the store can be large on disk yet cost <100 MB RAM and a few
seeks per off-target at query time.

### Runtime selection (in `new_simple_analysis.py`)

```
detection, rsID/MAF, single-variant AF (global/per-db/per-subpop)  -> Tier 0 only (always)
multi-variant off-target frequency + carrier aggregation           -> sparse genotype-tier
  reads at the off-target's positions, per database, then:
    per-database   = count carriers within that db's samples
    per-subpop     = bucket carriers by that db's subpop labels
    global         = union carriers across dbs, DEDUP by canonical sample ID
personal per-individual carrier list                               -> same reads, not aggregated
genotype tier absent but source VCF+tbi present                    -> tabix fetch fallback
genotype tier absent and no VCF                                    -> single-variant exact;
    multi-variant flagged "combination frequency needs sample tier"
```

## 5. Frequency semantics (exact)

For an off-target using variant positions `P = {p1..pk}` with required alleles `a1..ak`:

- **carrier (individual) frequency** in group `G` = |{ individuals in `G` carrying the
  combination }| / |G|. "Carrying" = phased: some haplotype has all `ai` in cis; unphased:
  individual is (het/hom) at every `pi` for `ai` (assume-cis upper bound, flagged).
- **allele (haplotype) frequency** in group `G` = |{ haplotypes in `G` with all `ai` in
  cis }| / (2·|G|) — phased only; for unphased we report carrier frequency and mark
  allele frequency as not-well-defined.
- `k == 1` (single-variant): both reduce to the position's AC/AN in `G` → **Tier 0**.
- `k >= 2`: computed from the genotype tier reads at `P`.

**Group `G` instances:** each database (within-db samples), each subpopulation (within-db,
native labels), and **global** (deduplicated union of all individuals across databases).

## 6. Multi-database & merging

- **Never merge into one cross-database AF by pooling raw.** Report per database + per
  subpop with native labels; subpop taxonomies are not reconciled across databases.
- **Global view = deduplicated union.** Correct **iff** individuals are deduped (each
  unique person once). Dedup needs a canonical cross-database sample ID. 1000G and HGDP are
  **disjoint cohorts** (distinct IDs) → union is clean and `AF_global` ≈ today's single
  `AF` column, made explicit. When adding an **overlapping** resource (gnomAD ⊇ 1000G/HGDP,
  UKB): dedup by ID; if IDs can't be matched, **exclude** that database from the global
  number or mark it an **upper bound**. Label the number `AF_global (N=<deduped>)`.
- **Interpretation caveat (surfaced in output/docs):** pooled cohorts are a *convenience
  sample*, not a random population draw; `AF_global` is "fraction of sampled individuals,"
  biased by cohort ancestry composition (same caveat gnomAD's global AF carries).
- **Provenance.** The `samplesID` files are `#SAMPLE_ID POPULATION_ID SUPERPOPULATION_ID
  SEX` with **no database column**, and the web's combined `hg38_1000G_HGDP.samplesID.txt`
  concatenates them (provenance lost). Build-time aggregation reads the **per-database**
  samplesID files (before concatenation); a `DATABASE` column is added to the combined
  samplesID for runtime.
- **Disk composition is linear.** Each database is its own Tier-0 AF contribution +
  genotype block. Adding a database = add blocks + a taxonomy entry; no cross-product, no
  re-encode of existing databases.

## 7. Disk & performance

| | Today (per-sample dicts) | v2 default (Tier 0) | v2 with genotype tier |
|---|---|---|---|
| Always-shipped artifact (genome, 1000G+HGDP) | ~40–50 GB (gz) | **~0.3–0.6 GB** | same + optional genotype tier |
| RAM / worker | 3–26 GB | **<100 MB** | <100 MB (sparse seeks) |
| Single-variant off-target AF (all 3 views) | full parse | **Tier 0, exact** | Tier 0, exact |
| Multi-variant off-target AF | full parse | flag / VCF-tabix | **sparse seeks at k positions** |

Tier 0 grows vs v1 (~0.15 GB) because it now carries per-(database×subpop) + global AF, but
stays sub-GB genome-wide. The genotype tier is far smaller than the 152 GB JSON (sparse
carrier indices / 2-bit packing) and is read only at found positions.

## 8. Build / publish / download

- Build (`crisprme.py build_variant_index`, after `add-variants`): compile
  `SNPs_genome/*.json` → (a) Tier 0 `registry_<vcf>/reg_<chrom>.bin` with per-(db×subpop) +
  deduped-global AF (aggregated from the per-database `samplesID` labels), and (b) per-db
  `genotypes_<db>/gt_<chrom>` compact stores. The raw JSON becomes a deleted intermediate.
- Publish/download (`crisprme_hf.py`): **Tier 0 rides in the main index tarball** (tiny →
  batteries stays self-sufficient for detection + single-variant AF); **each database's
  genotype tier publishes as a separate optional `genotypes_<db>.tar.gz`**. Default index
  download stays small; a user opts into genotype tiers per database as needed.
- Self-sufficiency check requires only Tier 0.

## 9. Read-site changes (contract-preserving)

- Keep `retrieveFromDict`'s 5-tuple `(snp_list, sample_list, rsID_list, AF_list,
  snp_info_list)`; fill metadata + single-variant AF from Tier 0; `sample_list` and
  multi-variant combination counts from the genotype provider (empty when not needed).
- `_load_dict_targeted` → `_open_registry` + a genotype-provider factory; `haplotype_check`
  from the manifest `phased` flag, not a scan.
- Summaries: `process_summaries.py` / `PopulationDistribution` consume the new per-(db×
  subpop) + global aggregates directly (Tier 0 for single-variant; genotype-tier
  aggregation for multi-variant) instead of re-deriving from a per-sample Samples column.
- Legacy `my_dict_<chrom>.json[.gz]` path retained as a detection-based fallback.

## 10. Backward compatibility

5-tuple contract preserved → downstream columns untouched. `_open_registry` falls back to
the legacy dict when no registry is present, so existing installs keep working (no flag
day). Output gains explicit `AF_global`/per-db/per-subpop + carrier columns; existing
columns keep their meaning.

## 11. Rollout

1. **Tier 0 registry** (compiler + reader + single-variant AF path, all three views) with
   legacy-dict fallback for multi-variant/NEED_SAMPLES. Ships the shrink + single-variant
   exactness first, smallest correctness surface.
2. **Sparse genotype tier** (compact per-db store + position index) for multi-variant
   combination AF + carrier aggregation + dedup-global; make it the default sample provider.
3. **tabix-VCF fallback** for the genotype tier when source VCFs are present but the store
   isn't shipped.
4. Extend to indel logs (`analisi_indels_NNN.py`); retire legacy dict shipping (keep
   `--legacy-dicts` one release).

## 12. Golden-output diff gate (mandatory before landing each phase)

Assert byte-identical results (old per-sample path vs new) for one guide on:

- a **phased** dataset (1000G) and an **unphased** dataset;
- targets that include **multi-variant haplotype** off-targets (not just single-SNP);
- the `PopulationDistribution` summary **and** all three AF views (global/per-db/per-subpop),
  both allele and carrier frequency;
- the **CPS1** example (chr2:210,530,658, rs114518452) as a regression anchor.

## 13. Open questions / risks

- **Allele vs carrier frequency for unphased multi-variant** off-targets: confirm the
  assume-cis upper bound is the desired behavior vs reporting only carrier frequency.
- **Cross-database ID canonicalization** for dedup when gnomAD/UKB arrive (ID schemes,
  cryptic relatedness) — may need an explicit sample-map artifact.
- **Genotype-tier encoding**: sparse carrier indices vs 2-bit matrix vs roaring bitmaps vs
  reusing the VCF+tbi directly — pick per measured size/seek cost (see §7; open to the
  optimization review).
- **Multiallelic / MNV / overlapping-variant** positions: confirm per-alt AC/AF and the
  combination logic compose correctly.
- **Biggest risk:** phasing/haplotype correctness → mitigated by the §12 golden diff on a
  phased and an unphased dataset with real multi-variant targets.
