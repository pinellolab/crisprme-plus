# Design: dictless variant post-analysis (replace the 152 GB per-sample dicts)

> Status: **design for review** (parallel exploration, not yet implemented). The
> shipping fix for the batteries gap is the separate "compress dicts + read on the
> fly" change (branch `feat/compressed-dicts-on-the-fly`, alpha.27). THIS document
> is the follow-on that makes the dicts tiny + fast rather than merely compressed.

## Context
The variant post-analysis reads per-chromosome SNP dictionaries
(`Dictionaries/dictionaries_<vcf>/my_dict_<chrom>.json`) that map
`"<chrom>,<pos>" -> $-joined records of "<samples>;<ref,alt>;<rsID>;<AF>"`, where
`<samples>` is a per-sample `sampleID:genotype` list. For 1000G+HGDP these total
**~152 GB** — the reason the batteries variant index isn't self-sufficient and why
post-analysis needs ~3–26 GB RAM/worker.

## The measurement that drives the design
Measured on a real shipped chr22 dict (1.55 GB uncompressed, 219 MB gzip, 993,880
variant positions):

| Quantity | Value |
|---|---|
| Per-position metadata only (pos, ref, alt, rsID, AF), gzip | **4.0 MB = 0.26% of the dict** |
| Per-sample genotype lists | **99.74% of the dict** |

The columns that must work out of the box (off-target detection, PAM-creation,
MAF/rsID) need only the **0.26%** metadata. Per-sample attribution needs the 99.74%
— exactly what we make **optional/lazy**. Genome-wide the metadata tier is
**~100–150 MB** vs **152 GB** today.

## Options considered
| | (A) tabix-VCF on the fly | (B) tiered: compact registry + optional per-sample | (C) position-indexed BGZF binary |
|---|---|---|---|
| Tiny always-shipped MAF/rsID/PAM artifact | no | **yes (~4 MB/chr)** | no (~dict size) |
| Batteries self-sufficient w/o source VCFs | **breaks** | **kept** | kept |
| Disk vs today | ≈ source VCF (50–100+ GB) | **~0.15 GB + optional sample tier** | ≈ today |
| Per-sample attribution | on-the-fly tabix | lazy tier (store or tabix) | always present |
| Default-run RAM / speed | <100 MB / medium | **<100 MB / fast** | <100 MB / fast |

**Recommendation: (B) tiered, with (A) as the lazy sample-tier fallback and (C)'s
seek trick for the sample tier.** (A) alone is a non-starter (breaks batteries); (C)
doesn't deliver the tiny always-shipped artifact.

## Recommended architecture — three artifacts per dataset
- **Tier 0 — position registry (ALWAYS shipped, ~4 MB/chr):**
  `Dictionaries/registry_<vcf>/reg_<chrom>.bin` (+ `.idx`) — sorted binary of
  `(pos, ref, alt, af, rsid)`, loaded fully into RAM per chromosome (~MBs). Powers
  off-target "is this position variant", PAM-creation (ref+alt), and the rsID/AF/SNP
  columns. **No per-sample data.** Plus a `manifest.json` with
  `{sample_count, phased, chroms}` (`phased` replaces the on-the-fly `|`/`/` scan).
- **Tier 1 — per-sample genotype store (OPTIONAL, lazy, published separately):**
  `Dictionaries/genotypes_<vcf>/gt_<chrom>.bgz` (+ offset `.gti`), BGZF-blocked,
  position-sorted; payload is the **existing** `samples;ref,alt;rsID;AF` value string
  (byte-identical) so `retrieveFromDict`'s parser is unchanged — O(1) seek instead of
  streaming a 13 GB JSON.
- **Tier 1-alt — tabix VCF mode (no store shipped):** when the sample tier is absent
  but source `.vcf.gz`+`.tbi` exist, resolve carriers via `pysam.VariantFile.fetch`
  (reuse `convert_gnomAD_vcfs.py:150-194`, `annotation.py:47-73`).

Runtime selection (new, in `new_simple_analysis.py`):
```
NEED_SAMPLES = personal_cards OR summaries OR Variant_samples requested
not NEED_SAMPLES        -> Tier 0 only            (default: fast, tiny)
elif genotypes_<vcf>/   -> Tier 0 + Tier 1 (seek)
elif VCF+tbi present    -> Tier 0 + Tier 1-alt (tabix)
else                    -> Tier 0 only + warn "sample attribution unavailable"
```

## Read-site changes (file:line, contract-preserving)
- `new_simple_analysis.py:157-199` `retrieveFromDict` — keep the exact 5-tuple
  `(snp_list, sample_list, rsID_list, AF_list, snp_info_list)`; fill all but
  `sample_list` from Tier 0; `sample_list` from the active provider (`[]` when not
  `NEED_SAMPLES` — already a supported branch). `iupac_decomposition` haplotype loop
  (`:238-288`) is **untouched**.
- `new_simple_analysis.py:823-867` `_load_dict_targeted` — **replaced** by
  `_open_registry` (+ a gt-provider factory); `haplotype_check` comes from the
  manifest, not a genotype scan. Keep the legacy `.json[.gz]` path as a fallback.
- `:792-820` `_collect_needed_dict_keys` (needed positions) + `:916-931` bootstrap —
  updated to open the registry + choose the gt provider.
- Consumers unchanged as long as the 5-tuple + `Samples`/`Variant_samples_*` columns
  are preserved: `change_headers_bestMerge.py:46-129`, `process_summaries.py:118,153`,
  `generate_sample_card.py:18,61-71`, `CRISPRme_plots_personal.py`. Tier-0-only runs
  leave the Samples column empty → summaries/cards degrade gracefully (the intended
  "optional tier").

## Build / publish changes
- `crisprme.py:1811-1863` (`build_variant_index` STEP 1): after `add-variants`, compile
  `SNPs_genome/*.json` → `registry_<vcf>/reg_<chrom>.bin` (Tier 0) + optional
  `genotypes_<vcf>/gt_<chrom>.bgz` (Tier 1, reuse `pysam.BGZFile`); the raw JSON
  becomes a deleted intermediate.
- `crisprme_hf.py` `_make_index_tarball`/`publish_index` (`:416-579`): **registry tier
  goes in the main index tarball** (tiny → batteries stays self-sufficient); the
  **genotype tier publishes as a separate optional `genotypes_<vcf>.tar.gz`**. Default
  index download drops from ~40–50 GB of dicts to **~0.15 GB registry**. Extract/route
  (`:370-385`) always installs `registry_<vcf>/`; opt-in fetch for `genotypes_<vcf>/`.
  Self-sufficiency check (`crisprme.py:554-567`) requires only the registry.

## Backward compatibility
5-tuple contract preserved → downstream untouched. `_open_registry` falls back to the
legacy `my_dict_<chrom>.json[.gz]` when no registry is present, so old installs keep
working — non-breaking, detection-based (not a flag day).

## Size / speed vs today
| | Today | Tiered default (Tier 0) | With sample tier |
|---|---|---|---|
| Shipped dicts (genome, 1000G+HGDP) | ~40–50 GB | **~0.1–0.15 GB** | +optional ~40–50 GB (separate) |
| chr22 artifact | 219 MB | **4 MB** | 4 MB + ~219 MB on demand |
| RAM / worker | 3–26 GB | **<100 MB** | <100 MB |
| Default-run speed | stream ≤26 GB | **~50–100× faster** | ≈ today when samples needed |

**Net: ~300× smaller always-shipped artifact, ~50–100× faster default post-analysis,
RAM bottleneck gone.** Per-sample attribution costs the same but only when requested.

## Phased rollout (recommend build Phase 1 now)
1. **Tier 0 registry** — compiler + reader + `retrieveFromDict` metadata path, **legacy
   dict fallback retained** for `NEED_SAMPLES` runs. Ships the ~300× shrink with the
   smallest correctness surface. **Gate on a golden-output diff** (one guide, old path
   vs new, assert identical `.bestCFD`/`.bestmmblg`), on a phased (1000G) + unphased set.
2. **Tabix VCF sample-tier fallback** (Tier 1-alt) — per-sample attribution without
   shipping the store, when VCFs are present.
3. **BGZF genotype store** (Tier 1) + separate optional publish/download; make it the
   default sample provider; retire legacy dict shipping (keep `--legacy-dicts` one release).
4. Extend the tiering to indel logs (`analisi_indels_NNN.py`).

**Biggest risk:** phasing/haplotype correctness → mitigate with the golden diff on a
phased and an unphased dataset.
