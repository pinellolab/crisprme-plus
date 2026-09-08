# Release checklist — CRISPRme+ 2.5.1 (STAGED — do not execute without sign-off)

Prepared for review. The code + docs are on `dev` and validated; this lists the scope,
the evidence, the exact steps, and the gated decisions. **No tag / image build / HF publish
has been done.** v2.5.0 is already tagged; this is the next release (`crisprme.py` = `2.5.1-dev`).

## Scope (what 2.5.1 adds over v2.5.0)
- **Two-pass `--fast` mode** (opt-in) for dense/aggregate panels — worst-possible reps instead
  of the 2^k IUPAC lattice. CFD exact/conservative, CRISTA screen-grade; non-`--fast` byte-identical.
- **All-source "mega" sites-only index** (1000G-2021 + HGDP + gnomAD + TOPMed + AoU) with
  per-dataset AF + `AF_max`, and the `indel_af.tsv` report companion + carrier/hom "NA" for
  aggregate groups.
- **`download` ships the reference index** alongside a variant index.
- **Performance:** CRISTA model cached (load once) + get_features N-guard (~1.6× on get_features),
  both byte-identical.
- **Fixes:** post-analysis Pool deadlock at high `--thread`; pre-flight guard for a variant
  search with no prebuilt index; exact worst-case CFD on the legacy/mega fast path.

## Validation evidence (done)
- **Clean-room real-user e2e** (`cleanroom_250`, #178): TRAC guide, 2021 panel, log_error empty,
  cooc in report.
- **V1/V2/V3 GW matrix:** V1 (2.4.0/2019) 325,530 → V2 (2.5.x/2021) 407,222 → V3 (`--fast`) 416,500.
  Delta computed on the **24 shared primary contigs** (V1's index also carries 10,620
  alt/random/decoy-contig rows absent from V2's primary-only index — excluded to avoid a
  contig-set confound). V1→V2 gain = **105,471 new off-targets: ~96.3% panel-driven** (2021
  variant-created) **/ ~1.2% feature** (cooc; only 14 net-new beyond the panel), ~3.7%
  high-edit boundary churn. **5 material reference off-targets (CFD 0.22–0.38, all mm+b 6–7
  with bulges) are genuinely absent in V2** — the denser 2021 variant-enriched index masks
  the reference-allele k-mer at those loci (compounded by the `--max-total-edits` default);
  none are strong/high-risk.
- **`--fast` GW-lossless:** V3 locus-level superset of V2; **0 of 1,458 CFD≥0.2 loci lost or
  demoted**; CFD exact-or-conservative; CRISTA screen-grade near 0.2 only.
- **cooc:** byte-identical fast vs non-fast; 843 CONFIRMED / 1,886 PUTATIVE GW (see
  `docs/GW_COOC_SUMMARY_2.5.1.md`).
- **CRISTA perf:** byte-identical to the shipped model (score `0.9060622093602125`), model 1× vs 3× load.

## Release steps (execute via `/release-crisprme` when approved)
1. Land `dev` → `main` (release is cut from `main`; `main` must be green).
2. `scripts/prepare_release.py 2.5.1` — drops `-dev`, bumps version everywhere (`crisprme.py`,
   Dockerfile, recipe).
3. Move `CHANGELOG.md` `[Unreleased]` → `[2.5.1]` dated section + update the link footer.
4. Create GitHub release/tag `v2.5.1` → triggers `.github/workflows/docker-multiarch.yml`
   (needs `DOCKERHUB_USERNAME`/`DOCKERHUB_TOKEN`).
5. SIF: `apptainer pull` the published image for ml007/ml008.
6. Bioconda: in progress — **blocked** on CRISPRitz 2.8.1 for py3.11; does not block the
   Docker/GitHub release.

## Gated decisions / open items (need your call)
- **Go/no-go on the tag** — everything else is a button-press after that.
- **Index republish:** NOT required for the 2.5.1 *code* release — the shipped HF indices
  (2021 feature-on, mega, old 2019) are compatible. Republish only if you want the perf/doc
  changes reflected in a bundled index (they don't change index format).
- **`--fast` indel parallelization (~cores×):** DEFERRED — needs a risky `analisi_indels`
  `__main__`-guard refactor + real-run validation. Not in 2.5.1; tracked for a follow-up.
- **#185 mega `indel_af` live demo:** reclassified — the mega has no *searchable* indels (empty
  fake-indel genome, SNP-only build), so the companion is code-verified but not live-exercisable
  on the mega; not a blocker (see `mega-indel-stores-empty-not-searchable`).

## NOT to do without explicit approval
Tagging, image build, HF publish, `main` merge, `prepare_release.py`. All gated.
