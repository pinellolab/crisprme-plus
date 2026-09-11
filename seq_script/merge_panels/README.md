# Combined multi-source VCF panels for CRISPRme (single-scan enrichment)

There are **two** merge modes, for two different purposes. Pick by whether the
sources are all genotyped and whether you need cross-variant cis:

| | **Mode 1 — genotyped panel** | **Mode 2 — sites-only "mega"** |
|---|---|---|
| Script | `merge_vcf_panels.sh` + `build_combined_panel.sh` | `merge_mega_sites.sh` + `build_mega_registry.py` |
| Sources | fully genotyped (1000G-2021, HGDP) | any, incl. aggregate (gnomAD AF-only, TOPMed AN=0, AoU pseudo-sample) |
| Genotypes | **kept** (per-sample GT) | **stripped** (`view -G`) |
| AF reported | **pooled** `AF` (`+fill-tags` AC/AN over union) + `AF_<src>` | per-source `AF_<src>` + **`AF_max`** (pooled AC/AN not computable) |
| Registry | genotype-counted (`compile_registry`, via `--samplesID`) | INFO-AF (`compile_registry_from_info_af`, no samples) |
| indel+SNP cis | CONFIRMED (1000G phased) / PUTATIVE (HGDP) | **N/A** (aggregate, no phase) |
| MAF filter | none | **> 0.001** (uniform) |

**Status:** Mode 1 validated **genome-wide** (1000G + HGDP), both NGG and pamless
NNN indexes built from one enrichment, pooled AF confirmed. Mode 2 merge + registry
validated **genome-wide** (all 5 sources, 66.9M sites, 51.2M SNPs); hybrid index
assembly (enrich + drop-in registry) pending. See §"Mode 2" below and
`docs/DESIGN_mega_index.md`.

---
## Mode 1 — genotyped cis-capable panel (1000G-2021 + HGDP)

## Motivation
CRISPRme enriches a reference with a VCF's variants and scans once. Running N
resources (1000G, HGDP, AllOfUs, TOPMed, gnomAD, …) as N separate datasets means N
enrichments + N scans + a post-hoc result merge. **Merging the resources into one
panel first → enrich once, scan once** — a large, growing time saving. The catch is
doing the merge without losing information we care about (per-source allele
frequency and which dataset a variant came from).

## What the merge produces (provenance, per variant)
| INFO field | Meaning |
|---|---|
| `AF` | **Pooled** allele frequency, recomputed as AC/AN across ALL merged samples (`bcftools +fill-tags`). This is what CRISPRme's enricher reports. |
| `AF_<SRC>` | Each source's **original** AF, verbatim; `.` where that source lacked the variant. |
| `SRC` | Comma-list of sources containing the variant: `1000G,HGDP` (shared), `HGDP` (HGDP-unique), … Derived from which `AF_<SRC>` are present. |

So provenance is recoverable two ways: implicitly (`AF_1000G` present ⇔ in 1000G) and
explicitly (`SRC`). A variant **unique** to one source is immediately identifiable.

## How it works (see `merge_vcf_panels.sh`)
Per chromosome, for each source: normalize contigs to `chr`-prefixed (to match the
hg38 reference), rename `INFO/AF → INFO/AF_<SRC>`, index. Then
`bcftools merge -m none` the sources, `bcftools +fill-tags -- -t AF` to compute the
pooled `AF`, and derive `INFO/SRC` from a fast sites-only pass over the `AF_<SRC>`
fields. `build_combined_panel.sh` runs this over all chromosomes (parallel,
resumable), assembles the panel `samplesID` (union), and calls
`crisprme.py build-index-only` to enrich + index.

## Validation

### Genome-wide (1000G + HGDP) — the full build
- **All 23 chromosomes merged** (chr1–22 + chrX), 3,477 samples each, AF-provenance
  intact (`AF` pooled + `AF_1000G`/`AF_HGDP` + `SRC`).
- **Both combined indexes built from ONE shared enrichment:**
  `NGG_3_hg38+hg38_1000G_HGDP` (8.7 GB + 730 MB INDELS) and the pamless
  `NNN_3_hg38+hg38_1000G_HGDP` (144 GB + 15 GB INDELS).
- **Pooled AF correct at scale:** `chr22:10516173 → AF=0.024992` (`AF_1000G=0.02`,
  `AF_HGDP=0.030541`, `SRC=1000G,HGDP`).
- **Why this matters — the earlier plain merge was wrong:** a naive
  `bcftools merge` (no `+fill-tags`) leaves shared 1000G∩HGDP variants carrying only
  ONE source's AF (e.g. the HGDP 0.0305 instead of the pooled 0.0250 above) — 15/1000
  sampled sites mismatched. This tooling recomputes the pooled AF, so the shipped
  combined index reports correct frequencies.
- **max-edits tolerance:** the pamless index is searched with `--max-total-edits 5`
  (CRISPRitz ≥2.8.1); this bounds the pamless search-space explosion in
  variant-dense regions (validated separately).

### chr22 detail (provenance cross-check)
- **Samples:** 3,477 = 2,548 (1000G) + 929 (HGDP); cohorts disjoint (0 shared).
- **Record provenance cross-checks exactly:** SRC = 572,013 `1000G` + 697,942 `HGDP`
  + 487,066 `1000G,HGDP` = 1,757,021 total; and 572,013+487,066 = **1,059,079**
  (= 1000G source records) and 697,942+487,066 = **1,185,008** (= HGDP source records).
- **Pooled AF correct:** e.g. chr22:10516173 `AC=156 AN=6242 → AF=0.024992`, with
  `AF_1000G=0.02`, `AF_HGDP=0.030541` preserved.
- **Enricher reads the pooled AF:** CRISPRitz #36's exact-key `AF` match selects
  `AF=0.024992`, correctly skipping the earlier `AF_1000G`/`AFR_AF` fields (the old
  2-char-prefix match would have grabbed `AF_1000G`).
- **End-to-end:** the panel enriches hg38_chr22 and builds the pamless `NNN_3` index
  (+ INDELS) without error. [`e2e_afmerge.sh`]

## Trade-offs (intentional, documented)
- **Phasing is not preserved across sources.** Cohorts are disjoint, so cross-source
  multi-variant haplotypes have no real carrier; single-variant off-targets (the vast
  majority) are fully correct. The single-scan speed win is the point.
- **Sites-only sources** (gnomAD, TOPMed: no per-sample genotypes) contribute variant
  positions + `AF_<SRC>` but no samples; the pooled `AF` then reflects only the
  genotyped sources. Add them as extra `SOURCES` entries.
- `bcftools merge -m none` can emit same-POS records for multiallelic sites; #36's
  enricher guards the AF-count mismatch this can create.

## Dependencies
- `bcftools >= 1.18` with the `+fill-tags` plugin; `htslib` (`bgzip`/`tabix`).
- **CRISPRitz #36** (enricher AF/FILTER robustness) for enrichment (exact-key pooled
  AF read + multiallelic AF-count guard) and **#37/#38** (`--max-edits`) for the
  bounded pamless search — both shipped in **CRISPRitz ≥ 2.8.1**.
- **Scratch space:** the per-chromosome intermediates are multi-GB and a genome-wide
  N-way-parallel run writes tens–hundreds of GB of temp. `merge_vcf_panels.sh` now
  puts its work dir on `$TMPDIR` (fall back: the output volume) rather than `/tmp` —
  point `TMPDIR` at a volume with ample free space (a genome-wide 1000G+HGDP run
  needs well over 100 GB of temp + ~160 GB for the final pamless index).

## Reproduce
```bash
# one chromosome (configure SOURCES in the script):
BCFTOOLS=bcftools ./merge_vcf_panels.sh /path/to/cwd/VCFs/hg38_1000G_HGDP chr22
./validate_merge.sh /path/to/cwd/VCFs/hg38_1000G_HGDP/merged.chr22.vcf.gz 1000G HGDP
# genome-wide + enrich + index:
./build_combined_panel.sh /path/to/crisprme_working_dir hg38_1000G_HGDP 20bp-NNN-NO-PAM 2 2 4
```

## To compare with Ann
- Merge mechanics: `-m none` vs `-m both`/normalization; how you handle multiallelic
  records and same-POS duplicates.
- Provenance: is `SRC` (+ `AF_<SRC>`) the same shape you use, or a different tag?
- Pooled vs per-source AF: do you recompute pooled AF, keep per-source, or both?
- Sites-only sources (gnomAD/TOPMed) and phasing policy.

---
## Mode 2 — sites-only "mega" all-source panel (`merge_mega_sites.sh`)

The mega merges **all five** sources — 1000G-2021, HGDP, gnomAD v4.1, TOPMed,
All-of-Us — into one **sites-only** index so a user scans once against the union of
common variation, with each source's frequency kept separate plus a global max. It
exists because Mode 1's pooled AF is only meaningful when every source is genotyped;
gnomAD is AF-only, TOPMed ships `AN=0`, and AoU is a single pseudo-sample, so there
is no honest pooled AC/AN and no cross-source cis. Individual-level data is
intentionally discarded.

### What it produces (per variant)
| INFO field | Meaning |
|---|---|
| `AF_<src>` | each source's **original** per-alt AF; absent where the source lacks the variant |
| `AF_max` | max over the present `AF_<src>` — the global frequency the mega registry's GLOBAL group carries |

There is **no pooled `AF`** (unlike Mode 1): ΣAC/ΣAN is not computable across these
sources, and averaging would be dishonest.

### How it works (`merge_mega_sites.sh`, per chromosome)
Per source: `bcftools norm -m -any -f <chr>.fa` (split multiallelics **and**
left-align, so the same indel from two sources merges instead of duplicating — AF is
`Number=A` so the per-alt AF is carried with no genotype read) → MAF filter
`-e "AF<=0.001 || AF>=0.999"` (keeps 0.001<AF<0.999) → `view -G` (strip genotypes) →
`annotate -x "^INFO/AF"` (keep **only** AF; note the `^` — `"INFO,^INFO/AF"` would
strip everything) → rename `INFO/AF → AF_<src>`. Then `bcftools merge -m none` the
sources (distinct `AF_<src>` co-exist at union sites) and add `AF_max` (query → awk
max → `annotate`). Genome-wide: `xargs -P N` over `chr1..chrY`.

Then `PostProcess/build_mega_registry.py` reads the merged `mega.<chr>.afmax.vcf.gz`
and builds the Tier-0 registry via `tier0_registry.compile_registry_from_info_af`
(each source → a db group with `AC=round(AF·AN_nom)`, `AN=AN_nom`=2×cohort-N so AF is
exact and AN is meaningful; GLOBAL group = `AF_max`). The registry is **SNP-only**
(the Tier-0 format packs single-char ref/alt); indels (~24% of sites) are counted and
routed to the fake-indel path, not annotated with AF in v1.

### Validation (genome-wide)
- Merge: 24 chroms, **66,913,681 union sites**; chr22 1.10M, chr1 5.52M, chrX 2.44M,
  chrY 187K (chrY correctly carries only HGDP/gnomAD/AoU — 1000G/TOPMed have none);
  `AF_max = max(AF_<src>)` verified (0 violations).
- Registry: **51,174,004 SNPs** (+15.7M indels routed out, 23.5%); read-back exact —
  chr1:10147 C>G → gnomAD 0.004 / AoU 0.056 / global=`AF_max` 0.056; chrX:10072 G>T →
  HGDP 0.065 / TOPMed 0.075 / global=`AF_max` 0.075.
- 13 unit tests: `PostProcess/test_registry_from_info_af.py`, `test_build_mega_registry.py`.

### Reproduce
```bash
# one chromosome (configure SOURCES + REF dir in the script / env):
CRISPRME_REF_DIR=/path/to/Genomes/hg38 ./merge_mega_sites.sh /path/to/VCFs/hg38_mega chr22
# genome-wide (6-way parallel):
printf '%s\n' $(seq 1 22) X Y | sed 's/^/chr/' \
  | CRISPRME_REF_DIR=/path/to/Genomes/hg38 xargs -P 6 -I{} ./merge_mega_sites.sh /path/to/VCFs/hg38_mega {}
# registry per chrom:
python3 ../../PostProcess/build_mega_registry.py --vcf /path/to/VCFs/hg38_mega/mega.chr22.afmax.vcf.gz \
  --chrom chr22 --out-dir /path/to/Dictionaries/registry_hg38_mega
```
Then the hybrid index assembly (enrich via `build-index-only` **without**
`--samplesID`, then drop `registry_hg38_mega/` in) — see `docs/DESIGN_mega_index.md`.

### Gaps / open items
- **Indel per-dataset AF absent in v1** (SNP-only registry). Indel off-targets are
  still discovered via fake-indel contigs, just un-annotated. Fast-follow: a
  sites-only indel-AF store.
- **Carrier/hom counts are HWE expectations** from AF (AF itself is exact) — the
  registry taxonomy sets `aggregate_af_only: True` so the report can mark them.
