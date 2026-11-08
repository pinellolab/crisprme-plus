# CRISPRme+ "mega" all-source merged index — design + build status

**The sites-only, widest-provenance fallback among the published CRISPRme+ indexes
(alongside the HGDP+1kGP default, HPRC, and the reference-only index).**
This document describes the mega index. An earlier 1000G-2021 + HGDP genotyped/cis
hybrid was superseded by the re-phased HGDP1kGP default and retired in v2.7.8. The mega merges **all five** variant sources into
one **sites-only** index so a user can scan once against the union of common
variation across cohorts, with **per-dataset allele frequencies** + a **global
`AF_max`**. Individual-level (genotype) data is intentionally discarded (user
decision, 2026-09-01). NRG (NRG_3) only for now.

Sources (`/data/pinello/SHARED_DATA/CRISPRme_data/variants_datasets_20260827/`):
1000G-2021, HGDP, gnomAD v4.1, TOPMed, All-of-Us. MAF > 0.001 for **all** (uniform
threshold, user decision).

## Why the mega needed new code (the crux)
The 2.5.0 Tier-0 registry derives per-database AF by **counting per-sample
genotypes** (`compile_registry` / `compile_registry_panel`). The mega is
**sites-only** (`bcftools view -G` strips genotypes) — there are no genotypes to
count. The aggregate datasets' existing workaround (`convert_gnomAD_vcfs.py`
synthesises a *presence* genotype per subpop) yields a subpop-presence proxy, not
the true AF. For a paper-grade mega we want the **true** per-dataset AF. So the
mega reads the pre-computed `AF_<dataset>` INFO fields directly — a new registry
path. (User chose this over presence-genotypes or a no-AF sites-only index.)

## Build pipeline

### 1. Merge (per chrom) — VALIDATED on chr22, running GW
`/srv/local/lp698/mega_gw.sh` (resume-safe `.done` markers, `xargs -P 6`). Per
dataset, per chrom:
```
bcftools norm -m -any -f <chrom>.fa <src.vcf.gz> -Ou     # split multiallelics; AF is Number=A so no GT read
  | bcftools view -e "AF<=0.001 || AF>=0.999" -G -Ou      # MAF>0.001 filter; strip genotypes
  | bcftools annotate -x "^INFO/AF" -Ou                   # keep ONLY INFO/AF (NB: "^INFO/AF", NOT "INFO,^INFO/AF")
  | bcftools annotate --rename-annots <ren> -Oz -o <tag>.<chrom>.vcf.gz   # INFO/AF -> AF_<dataset>
```
Then `bcftools merge -m none` the present datasets → union sites, each carrying its
own `AF_<dataset>` (co-existing, not blended). Then `AF_max` (query→awk max over
present per-dataset AF→`annotate -c ...,INFO/AF_max`). Output:
`mega_gw/mega.<chrom>.afmax.vcf.gz`.

chr22 validation: 1,101,148 union sites; all 5 `AF_<ds>` present; 649,646 sites
(59%) in ≥2 datasets; `AF_max` = max per row (0 violations). Datasets missing a
contig (1000G/TOPMed have no chrY) are simply absent from that chrom's merge.

**GW merge DONE** (2026-09-02, ~1h40m, `mega_gw/mega.<chrom>.afmax.vcf.gz`):
**66,913,681 union sites** across 24 chroms (chr1 5.52M, chrX 2.44M, chr22 1.10M,
chrY 187K). Verified: all covered datasets' `AF_<ds>` + `AF_max` present per chrom;
chrY correctly carries only HGDP/gnomAD/AoU (no 1000G/TOPMed). Indels are a stable
**~24%** of sites (chr10 23.6%, chr22 23.9%).

Nominal cohort sizes (for the registry `AN_nom = 2N`): 1000G-2021 **3202**, HGDP
**929**, gnomAD **807162** (joint), TOPMed **53831**, AoU **~535662** (chr22 max
AN 1,071,324; Aug 2026 snapshot).

### 2. SNP registry (INFO-AF) — CODE DONE + tested + validated on real chr22
`PostProcess/tier0_registry.py::compile_registry_from_info_af` (+ helper
`_af_only_counts`) reuses `_write_registry` → **binary format byte-identical**
to `compile_registry` (read by the shipped `RegistryReader`, no reader change).
Per site: each dataset → one db-level group with `AC=round(AF*AN_nom)`,
`AN=AN_nom` (AF exact to ±0.5/AN — <8e-5 for 1000G, <3e-7 for gnomAD); the GLOBAL
group carries **`AF_max`** over the pooled nominal panel. Individual carrier/hom
counts do not exist for aggregate data, and none is estimated: they are stored as
**0** (`n_called` = the nominal cohort size); **AF is exact** (see *Carrier/hom
presentation* below). `aggregation` tag = `"info_af"`, taxonomy carries
`aggregate_af_only: True`.

Driver: `PostProcess/build_mega_registry.py` — dependency-free gzip VCF reader →
`reg_<chrom>.bin/.idx`. **SNP-only**: the Tier-0 registry packs single-char
ref/alt; indels are split out + counted (logged, not silently dropped).

Real chr22 run: **838,025 SNPs** compiled, **263,123 indels (24%) skipped**;
read-back of chr22:10516173 A>G is exact — HGDP AF 0.0307 (src 0.0305), TOPMed
0.0956 (src 0.0956), **global 0.0956 = AF_max**.

**GW registry DONE** (2026-09-02, `registry_mega/reg_<chrom>.{bin,idx}`, 24 chroms):
**51,174,004 SNPs** (+15.7M indels skipped, 23.5%). Read-back exact across chroms —
chr1:10147 C>G → AoU 0.056 / gnomAD 0.004 / global=AF_max 0.056; chrX:10072 G>T →
HGDP 0.065 / TOPMed 0.075 / global=AF_max 0.075.

Tests: `test_registry_from_info_af.py` (13), `test_build_mega_registry.py` (7) —
AF round-trip, GLOBAL=AF_max, sparsity, multiallelic, tiny-AF floor, SNP/indel
split, carrier/hom 0, RegistryReader round-trip. Existing `test_tier0_registry` (20) +
`test_tier0_compile` (19) still green.

### 3. Enriched genome + fake-indel contigs — DE-RISKED (sites-only OK)
`crispritz.py add-variants` on a sites-only (8-column, no FORMAT/samples) VCF
**parses cleanly** — "Variants Extraction and Processing" completes, no
sample-column dependency. (Standalone-run path/naming quirks are handled by
crisprme.py's own orchestration.) Naming contract: the enricher derives the genome
FASTA from the VCF **filename**, so per-chrom mega VCFs feed as `chr<N>.vcf`.

### 4. Assembly (HYBRID) — REMAINING (morning, needs SIF deploy of this code)
Naming: `--vcf <dir>` sets `vcf_name = basename(dir)`, which drives every path
(`Genomes/{ref}+{vcf_name}`, `genome_library/{pam}_{bMax+1}_{ref}+{vcf_name}`,
`Dictionaries/registry_{vcf_name}`). So place the 24 mega VCFs (named `chr<N>.vcf.gz`
to match the genome contig fastas — the enricher derives the genome file from the
VCF filename) in `VCFs/hg38_mega/`.

Steps (chr22 smoke first, then GW):
1. `crisprme.py build-index-only --genome <hg38dir> --pam <NRG_3 pam> --bDNA 2
   --bRNA 2 --vcf VCFs/hg38_mega/ --path <workdir>` **WITHOUT `--samplesID`** (mega
   has no samples). Produces the enriched SNP genome + `_INDELS` fake contigs +
   CRISPRitz index + legacy dicts. Omitting `--samplesID` means NO genotype-based
   registry is emitted (nothing empty to override).
2. Drop in the INFO-AF registry as `Dictionaries/registry_hg38_mega/reg_<chrom>.{bin,idx}`
   (rename `registry_mega/` → `registry_hg38_mega/` to match `vcf_name`), + a
   `variant_count.json` (per-dataset SNP counts).
3. Publish to HF (NRG_3, enriched genome + index + registry; **no** genotype store,
   **no** indel-GT store — aggregate data has no phase).

**Open questions the smoke must answer** (empirical — cannot resolve by reading):
- Does the dict-less **search prefer `registry_<vcf>/`** over the (empty, no-sample)
  dicts when built without `--samplesID`? If the search defaults to dict-based, the
  registry drop-in won't be consulted → need to force the dict-less path.
- Are the `_INDELS` fake contigs built + searchable without an indel-GT store
  (indel off-targets found, just AF-unannotated)?
- Exact NRG_3 PAM file + bulge geometry (match the shipped NRG index).

## Indel per-dataset AF (decided: build now)
The Tier-0 registry is SNP-only, so indels (~24% of sites) get their per-dataset AF
from a **separate sites-only sidecar**:
- `PostProcess/build_mega_indel_af.py` builds `indel_af_<chrom>.tsv.gz` (gzipped TSV,
  `pos ref alt AF_<ds>.. AF_max`, atomic write) from the merged VCF's INDEL records,
  and exposes `IndelAfReader.lookup(pos, ref, alt) -> {ds: af, "AF_max": af}`.
- `analisi_indels_NNN.py` loads the sidecar (sibling `indel_af_<vcf>/`, guarded
  no-op if absent) and emits a companion `<output>.indel_af.tsv` (one row per indel
  off-target: chrom, off-target start, indel pos/ref/alt, per-dataset AF, AF_max) —
  mirrors the cooc companion, never touches the fixed indel columns/scores.
- `merge_mega_sites.sh` also writes a plain `INFO/AF = AF_max`, after the site's
  `AF_<src>` keys and `AF_max`. On an index built from its output the **global** AF_max
  would reach an indel off-target through the existing fake-indel MAF path just where
  the indel's site has as many sources as the VCF's first record: CRISPRitz
  `add-variants` reads each record's `AF` at the place in `INFO` where the first record
  has it (a known limit, METHODS section 8), and this `AF`'s place moves with the number
  of sources at the site, so every other indel would get none. That is CRISPRitz 2.9.1
  and earlier: 2.10.0 finds each record's `AF` by its key, so on an index it builds
  every indel would get it. The as-run genome-wide
  merge that made the shipped mega's VCFs, `mega_gw_merge.sh` (read by
  `mega_build_indels.sh`), writes only `AF_<src>` and `AF_max`, with no plain `INFO/AF`,
  so the shipped mega's indel logs hold no AF (for none of its 13,853,761 indels,
  222,283 of them on chr22) and a row of an indel alone has no MAF there; the sidecar
  gives each indel's **per-dataset** AF and its `AF_max`, in the `indel_af.tsv`
  companion.
- Tests: `test_build_mega_indel_af.py`. Remaining: surface `indel_af.tsv` in the
  report (same merge pattern as the cooc companion) + validate end-to-end in the
  assembly smoke (the `analisi_indels_NNN.py` hook is guarded but untested in a live
  indel search).

## Carrier/hom presentation (decided: AF-only)
Aggregate data has no genotypes, so per-individual carrier/hom counts are **0** in
the registry (never a fabricated Hardy-Weinberg estimate). The taxonomy sets
`aggregate_af_only: True` and the manifest `aggregation: "info_af"`, from which the
population-summary companion (`*.population_summary.tsv`) renders every carrier/hom field
as NA and keeps AF + AF_max (`test_aggregate_carrier_na.py`); the report shows no
carrier/hom count from the registry. `targets-scoring`'s carrier overlay prints the stored 0.

## Other notes
- **No pooled global AF** (ΣAC/ΣAN not computable: gnomAD AF-only, TOPMed AN=0,
  cohort overlap) — the global group uses **`AF_max`**, per user.
- **Review fixes applied** (adversarial pass): reject non-finite / out-of-range AF
  before it `OverflowError`s a whole chromosome build; count multiallelic + no-AF
  drops (nothing silently discarded); uppercase the SNP gate; guard `n_called==0`;
  document that a `--datasets` subset recomputes GLOBAL.
- **Hardening round** (second workflow, 17 confirmed latent findings — 0 BLOCKER/MAJOR):
  the indel store now DEDUPS duplicate `(pos,ref,alt)` keys by **max** per-dataset AF
  (`bcftools norm` emits ~19 such keys genome-wide in VNTR regions; previously the
  reader silently kept last-wins); both builders honor `--chrom` (multi-contig
  safety — the registry is keyed by (pos,alt) with no chrom); the merge awk emits
  `.` not `0` for an all-AF-missing row; `_find_vcf` fails on ambiguity instead of
  silent `head -1`; chrom-token grep anchored to the basename; registry skips
  phantom empty-group records; `IndelAfReader` validates its header + parses AF_max
  by position; the indel hook warns-once instead of pure silent swallow.
- **Full 24-chrom verification** (`verify_mega_chrom.py`): registry + indel sidecar
  PARTITION all 66.9M sites (SNP registry XOR indel store, no gap/overlap) and
  sampled AFs round-trip on **every** chromosome — ALL PASS.

## Files
- `PostProcess/tier0_registry.py` — `compile_registry_from_info_af`, `_af_only_counts`
- `PostProcess/build_mega_registry.py` — merged-VCF → registry driver (CLI + `build()`)
- `PostProcess/test_registry_from_info_af.py`, `test_build_mega_registry.py`
- ml007: `/srv/local/lp698/mega_gw.sh` (GW merge), `mega_gw/` (outputs), `mega_dryrun/` (chr22 proof)
