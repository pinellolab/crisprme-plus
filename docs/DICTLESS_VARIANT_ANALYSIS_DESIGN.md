# Design: dictless variant post-analysis (v3)

> Status: **design for review (v3)** — not yet implemented. v3 incorporates the
> geneticist + systems/optimization review of v2 (both "sound-with-changes") and two
> product decisions (2026-08-17): **(1)** the frequency feature is scoped to
> variant-**created/strengthened** off-targets (weaken/abolish/ref-already-off-target are
> labeled, not given a misleading number); **(2)** ship the extra headline outputs
> (max-subpop AF, homozygous-carrier count, observed-in-≥1-genome, absolute carrier N).
> Implementation is gated on the two-part diff in §12.

---

## 1. Problem

Variant post-analysis (`new_simple_analysis.py`) reads per-chromosome SNP dictionaries
`Dictionaries/dictionaries_<vcf>/my_dict_<chrom>.json` mapping
`"<chrom>,<pos>" -> $-joined "<samples>;<ref,alt>;<rsID>;<AF>"` (`<samples>` = per-sample
`sampleID:genotype`). For 1000G+HGDP these total **~152 GB** (gz ~31 GB since alpha.27),
needing ~3–26 GB RAM/worker. alpha.27's gzip bought space, **not speed** (still a full
per-sample parse). This design makes them tiny **and** fast, and fixes several latent
frequency-correctness bugs on the way.

**Scale (measured, corrected from v2).** The shipped chr22 1000G dict has **993,880**
positions; chr22 ≈ 1.64% of hg38 → **~60M variant positions genome-wide** (1000G alone;
the combined 1000G+HGDP panel, ~4,430 samples, is larger). All v2 size numbers were ~20×
low; §7 is recomputed against ~60M.

## 2. What the user needs (outputs)

For each reported off-target, both **allele frequency** and **carrier (individual)
frequency**, plus absolute **carrier count (N)**, in three views:

- **Global** — across **all individuals in all loaded databases, deduplicated**.
- **Per database** — within each database's own individuals (1000G, HGDP, …).
- **Per subpopulation** — within a database, using that database's **native** labels.

Plus these headline columns (product decision, all cheap given AC/AN + the genotype tier):
- **max-subpopulation AF** and **which subpopulation** (the most decision-relevant number
  for ancestry-specific risk);
- **homozygous-carrier count/frequency** (biallelic risk);
- **observed-in-≥1-genome** boolean;
- **absolute carrier N** alongside each frequency.

Detection, PAM-creation, and rsID/MAF remain out-of-the-box. Per-**individual** carrier
lists (personal cards) stay available as the on-demand layer.

## 3. Correctness crux — off-target frequency is a COMBINATION property, with a direction

An off-target's frequency is the frequency of the **specific allele combination
(haplotype)** that produces it — not any single SNP's marginal AF, nor a product (LD), nor
min/max of per-SNP AFs. (This already fixes a real bug: `CRISPRme_plots_personal.py` takes
`min()` over per-SNP AFs for a multi-variant target — neither the haplotype AF nor a valid
bound.)

Two structural consequences:
1. **Off-targets are guide-dependent** → unknown at index-build time → per-off-target
   frequencies **cannot** be precomputed into the index; they are computed at
   **post-analysis**, for the found targets.
2. So the win is: compute the combination frequency **at query time, sparsely** — only at
   the handful of positions the found off-targets touch — not by streaming 152 GB.

### 3.1 Direction of effect (product decision: created/strengthened only)

A variant can (create | strengthen | weaken | abolish) an off-target, at the (seed | PAM).
The risk-relevant frequency differs by direction (alt-frequency for create/strengthen;
**ref**-frequency = 1−AF for weaken/abolish; ~1.0 for a reference-only off-target). To avoid
ever reporting a misleading number, **v3 scopes the frequency feature to
create/strengthen** (the classic CRISPRme question: "which variants introduce new/worse
off-targets, and how common are they"):
- Every off-target gets a `variant_effect ∈ {created, strengthened, weakened, abolished,
  reference}` label, derived from CFD_ref vs CFD_alt / `Var_uniq` (already computed today).
- **Frequency + carrier columns are populated only for `created`/`strengthened`** (the alt
  combination is the risk allele). For `weakened`/`abolished`/`reference`, frequency fields
  are `N/A` with the effect label, so no one reads an alt-AF as the risk when the variant is
  protective. (Reporting direction-correct ref-frequencies for those is a documented future
  extension; it also requires surfacing ref-only/abolished candidates the search drops today.)

### 3.2 Single- vs multi-variant (the load-bearing split)

- **Single-variant** created/strengthened off-target (one alt at one position — the common
  case): its frequency **is** that (position, alt)'s allele/carrier frequency → exact from
  **Tier 0**, per db and per subpop.
- **Multi-variant** (≥2 alts must co-occur in cis): frequency = co-occurrence of the exact
  combination → needs per-sample genotypes at those positions → **genotype tier**.

### 3.3 Phasing / cis

A multi-SNP off-target on one molecule needs the alts in **cis**. Phased (1000G) → exact
haplotype counting; **unphased** → report **carrier frequency as an upper bound**
(assume-cis) *with* a lower bound (`max(0, Σ single-SNP carrier freq − (k−1))`), and mark
the **haplotype allele frequency undefined** (not identifiable without phase). Carry a
per-database `phased` flag; state the mode in output.

## 4. Architecture — Tier 0 registry (mmap) + sparse per-database genotype tier

### Tier 0 — position registry (ALWAYS shipped, mmap'd)

`Dictionaries/registry_<vcf>/reg_<chrom>.bin` + `reg_<chrom>.idx`. **mmap'd, not parsed
into a dict** (largest chromosome ~5M positions ≈ 100–180 MB binary; a parsed dict would be
2–4× that and blow the RAM budget). Layout:
- **primary record array**, fixed-width, sorted by pos, keyed per **(position, specific
  alt)** (multiallelic → one record per alt): `pos:u32, ref:u8, alt:u8, rsid_off:u32,
  group_blob_off:u32, flags:u8` (flags: strand-irrelevant here; alt encodes the base).
- a **parallel sorted-pos array** for O(log n) binary search / mmap lookup.
- a per-record **group blob** (variable-length, offset-referenced) holding **AC/AN integers
  per (database × subpopulation)** and the **deduplicated-global AC/AN** — sparse (only
  groups with carriers). Store **AC/AN, not AF** (exact for rare alleles; lets downstream
  recombine groups; AF = AC/AN on read). Size AN for the largest cohort (u32).
- an **rsID string pool** (offset-referenced), shared per chromosome.

Powers, from Tier 0 alone: detection, PAM-creation (ref vs alt), rsID/MAF, and **exact
allele + carrier + hom + max-subpop frequencies for single-variant created/strengthened
off-targets**, all three views. **No per-sample data.** `manifest.json`:
`{databases:[{name, sample_count, phased, has_genotypes, subpops:[…]}], chroms,
global_sample_axis_size}`.

### Genotype tier — per-database, position-indexed, sparse, compact (OPTIONAL)

`Dictionaries/genotypes_<db>/gt_<chrom>.zst` + `.idx`, **one block per database** (not a
merged matrix). Per **(position, alt)** record, a **per-record HYBRID** encoding chosen by
min-size at build (the carrier distribution is bimodal — median 2 carriers, mean 139, p99
~2,300, max ~2,550 of ~3,500):
- `tag:u8` then EITHER **delta-sorted varint carrier indices** (rare loci, ~2 B/carrier) OR
  a **2-bit-per-sample bitplane** (dense loci, N/4 B);
- carrier indices are into the **build-time global deduplicated sample axis** (§6) so a db's
  carriers map into their slice and cross-db dedup is free;
- when **phased**, store 2-bit **per-haplotype** codes so cis combinations are countable
  (reproduces the existing haplotype set-arithmetic); when unphased, store genotype
  (het/hom) only.

**Small 4–8 KB zstd compression blocks** (NOT htslib's 64 KB BGZF): a single off-target
touches ~a few hundred bytes; 64 KB blocks give ~173× decompress amplification, 4–8 KB give
~11–22×. `.idx` maps `first_pos → block voffset`. Read **only at the found off-target
positions**. Store large on disk, cost <100 MB RAM + a few small-block seeks per off-target.

### Runtime selection (`new_simple_analysis.py`)

```
detection, rsID/MAF, single-variant created/strengthened AF (all 3 views + hom + max-subpop)
                                              -> Tier 0 only (always; mmap)
multi-variant created/strengthened AF/carriers-> genotype tier: read (pos,alt) records at the
   off-target's k positions, intersect carrier sets (cis for phased), then bucket:
     per-database  = popcount within that db's slice
     per-subpop    = popcount within that subpop's mask
     global        = popcount over the deduped global axis (dedup free by construction)
personal per-individual carrier list          -> same reads, not aggregated
db has_genotypes=false (gnomAD)                -> single-variant AF exact; multi-variant marked
                                                  "not computable (frequency-only database)"
genotype tier absent but source VCF+tbi present-> tabix fallback (secondary; NOT primary store)
```
Query-time combination counting is cheap — measured ~5 µs/off-target for the bitset
intersect + all-group popcounts (bitsets over ~4.4k samples = ~70 u64 words).

## 5. Frequency semantics (exact, AN- and ploidy-aware)

For a created/strengthened off-target needing (position, alt) pairs `P = {(p1,a1)…(pk,ak)}`
in group `G`:
- **allele (haplotype) frequency** = `AC_G / AN_G` over **called alleles** (AN, not 2·|G| —
  excludes missing `./.`). k=1 → the record's AC/AN in `G`. k≥2 (phased) → count haplotypes
  carrying all `ai` in cis. **Unphased k≥2 → undefined** (report carrier bound instead).
- **carrier (individual) frequency** = `carriers_G / called_individuals_G`, denominator
  excludes individuals missing at **any** `pi`. Policy: **drop-if-any-missing** (state it).
- **ploidy/sex-aware:** on chrX non-PAR / chrY, **males contribute 1 allele** to AC and AN
  (ploidy 1) — do **not** `hap*2` (the current indel code's `hap = hap*2` double-counts and
  is a bug to fix, not reproduce). Use the `SEX` column already in `samplesID`. PAR regions
  are diploid.
- **strand:** registry alleles are on the reference **plus** strand; the reader
  reverse-complements the required off-target base to plus-strand before matching (mirrors
  `reverse_complement_table` in `iupac_decomposition`).
- **multiallelic:** `P` carries the **specific alt**; AC/AN and carrier sets are **per-alt**,
  never site-total non-ref.
- **derived columns:** homozygous-carrier count (both/all copies carry the combination),
  max-subpop AF (+label) = argmax over that db's subpops, observed = carriers_global ≥ 1.

**AF is re-derived from genotypes + samplesID at build time**, per (db×subpop) and global —
**not** the source VCF INFO `AF` (which is pan-cohort/ release-specific and not
subpop-decomposed; used only as a build-time cross-check).

## 6. Multi-database & merging

- Report **per database + per subpop (native labels)**; never a pooled cross-db AF.
- **Global = deduplicated union.** Build a single **global deduplicated sample axis** (u32)
  at build time — the union of all databases' individuals, each canonical ID once; each db's
  carrier indices map into their slice. Then the global view's dedup is **free** (popcount
  over the axis). 1000G/HGDP are disjoint (clean union; `AF_global` ≈ today's single `AF`,
  made explicit). Overlapping resources (gnomAD ⊇ 1000G/HGDP, UKB): dedup by canonical ID via
  a build-time sample-map; if unmatchable, exclude that db from global or mark upper bound.
- **Interpretation caveat** (in output + docs): pooled cohorts are a convenience sample, not
  a random population draw; `AF_global (N=<deduped>)` = fraction of sampled individuals,
  biased by cohort composition (gnomAD carries the same caveat).
- **Relatedness:** AF/carrier estimates assume approximate independence of sampled genomes;
  1000G has trios/related samples, HGDP has related individuals. Where a db ships a
  related/unrelated list, prefer the unrelated subset for frequencies (or note the
  independence assumption). Record which subset was used.
- **Provenance:** `samplesID` is `#SAMPLE_ID POPULATION_ID SUPERPOPULATION_ID SEX` with **no
  database column**; the combined file just concatenates. Aggregate from the **per-database**
  samplesID files at build; add a `DATABASE` column to the combined samplesID for runtime.
- **Capability flags** per db in the manifest: `{has_genotypes, phased,
  has_relatedness_subset}`. Frequency-only dbs (gnomAD, genotype-less): single-variant AF
  exact from Tier 0; multi-variant explicitly "not computable."
- **Disk composition linear:** add a db = add its Tier-0 AC/AN contribution + genotype block
  + a taxonomy entry; no cross-product, no re-encode of existing dbs.

## 7. Disk & performance (recomputed at ~60M positions)

| | Today | v3 default (Tier 0, mmap) | v3 + genotype tier |
|---|---|---|---|
| Always-shipped (genome, 1000G+HGDP) | ~40–50 GB gz | **~1–2 GB** | + optional ~5–10 GB/db |
| RAM / worker | 3–26 GB | **<100 MB** (page cache, evictable) | <100 MB (small-block seeks) |
| Resident/chrom | full parse | tens of MB (mmap, largest chrom ~100–180 MB on disk) | same |
| Single-variant AF (3 views) | full parse | **Tier 0, exact** | Tier 0, exact |
| Multi-variant AF | full parse | flag (or tabix) | **sparse seeks at k positions, ~5 µs** |

Still a 20–40× shrink of the always-shipped artifact and RAM bottleneck gone. Re-run sizing
on the **actual combined 1000G+HGDP** dict (more positions + ~4,430 samples) before locking §7.

## 8. Build / publish / download

- Build (`crisprme.py build_variant_index`, after `add-variants`): **single streaming pass
  per chromosome** over `my_dict_<chrom>.json[.gz]` (ijson, bounded RAM), **parallel across
  24 chromosomes** (~1–2 core-hours genome-wide), emitting Tier 0 (`registry_<vcf>/`, per-alt
  AC/AN per (db×subpop)+global, computed from genotypes+samplesID with ploidy/sex/AN
  correctness) and the per-db `genotypes_<db>/` hybrid stores. Raw JSON becomes a deleted
  intermediate.
- Publish/download (`crisprme_hf.py`): **Tier 0 rides the main index tarball** (small →
  batteries self-sufficient for detection + single-variant AF); add a **`genotypes`
  component** (`_COMPONENT_PREFIXES`/`_COMPONENT_LOCALDIR`) publishing **`genotypes_<db>.tar.gz`
  separately/optionally**. Self-sufficiency check requires only Tier 0.

## 9. Read-site changes (contract-preserving)

- Keep `retrieveFromDict`'s 5-tuple; fill metadata + single-variant AF from Tier 0;
  `sample_list` + multi-variant combination counts from the genotype provider (empty when not
  needed). Preserve the phased **level-subtraction** logic (`new_simple_analysis.py:387-395`,
  the top correctness risk) exactly.
- `_load_dict_targeted` → `_open_registry` (mmap) + genotype-provider factory;
  `haplotype_check` from the manifest `phased` flag.
- `process_summaries.py` / `PopulationDistribution` consume the new per-(db×subpop) + global
  aggregates directly. Legacy `my_dict_<chrom>.json[.gz]` retained as detection-based fallback.

## 10. Backward compatibility

5-tuple contract preserved → downstream columns unchanged in meaning; new columns
(`variant_effect`, AF_global/per-db/per-subpop, carrier N, hom count, max-subpop, observed)
are additive. `_open_registry` falls back to the legacy dict when no registry is present (no
flag day).

## 11. Rollout

1. **Tier 0** (compiler + mmap reader + single-variant created/strengthened AF path, all
   views + derived columns), legacy fallback for multi-variant/NEED_SAMPLES. Ships the shrink
   + single-variant exactness + the AN/ploidy/per-alt fixes first.
2. **Genotype tier** (hybrid encoding, small-block zstd, global sample axis) for multi-variant
   combination AF + carrier aggregation + dedup-global; default sample provider.
3. **tabix-VCF fallback** when source VCFs present but the store isn't shipped.
4. Extend to indel logs (`analisi_indels_NNN.py`); retire legacy dict shipping (`--legacy-dicts`
   one release).

## 12. Two-part diff gate (mandatory per phase — supersedes v2's byte-identical gate)

The legacy path doesn't compute allele frequency and has AN/ploidy/multiallelic bugs, so a
byte-identical gate would force reproducing bugs. Split it:
- **REPRODUCE gate** — byte-identical on fields whose semantics are unchanged: detection,
  `Samples` set membership, rsID, `Var_uniq`, single-variant **carrier** frequency where the
  legacy value is already correct (autosomal, no missingness).
- **CORRECTNESS gate** — new/fixed fields (AN-corrected AF, ploidy/sex on chrX/Y, per-alt
  multiallelic, per-(db×subpop)+global, hom/max-subpop) validated against an **independent
  oracle** (`bcftools +fill-tags` / `plink --freq` on the same VCF subset), on a **phased
  (1000G)** and an **unphased** dataset, including **multi-variant haplotype** off-targets and
  the **CPS1** anchor (chr2:210,530,658, rs114518452).

## 13. Remaining open items (technical, no product decision pending)

- Cross-database ID canonicalization + relatedness subset artifact (needed when gnomAD/UKB
  arrive; non-issue for disjoint 1000G/HGDP now).
- Confirm multiallelic/MNV/overlapping-variant composition end-to-end on the **combined**
  enriched index (the chr22 1000G sample was pre-split by `bcftools norm -m-`; the combined
  index may carry `$`-joined multiallelic records).
- Lock §7 numbers against the actual combined 1000G+HGDP dict (positions + ~4,430 samples).
- **Top risk:** phasing/haplotype correctness + the level-subtraction reimplementation →
  mitigated by the §12 CORRECTNESS gate on phased + unphased data with real multi-variant targets.
