# Changelog

All notable changes to CRISPRme are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When cutting a release, move items out of `[Unreleased]` into a new dated
version section and update the link-reference footer. See `docs/RELEASING.md`
and the `release-crisprme` skill.

## [Unreleased]

## [2.3.0] - 2026-08-18

The dict-less variant-analysis engine: replaces the ~152 GB per-sample SNP
dictionaries with compact, memory-mapped tiers, adds population-level off-target
summaries, and corrects the allele-frequency column. Backward-compatible — a
dict-based install still works unchanged.

### Added
- **Compact dict-less variant post-analysis.** A tiny, always-shipped **Tier-0
  registry** (per-(pos,alt) AC/AN + per-(db×subpop)/global allele counts + rsID,
  mmap+bisect) powers off-target detection and corrected allele frequencies out of
  the box; an optional **Tier-1 genotype store** reconstructs the exact per-sample
  Samples column. Together they replace the 152 GB per-sample dicts with ~7 GB
  (registry) + ~22 GB (genotype tier), random-access — SNP post-analysis no longer
  streams 152 GB of per-sample JSON. `build-index-only` emits the tiers; the search
  auto-detects and uses them (falling back to dicts when present).
- **Combination-aware population-summary output** — a `<output>.population_summary.tsv`
  companion per variant off-target with per-database, per-superpopulation and global
  allele/carrier frequencies, max-subpopulation AF (+label), homozygote counts and
  absolute Ns, with **dataset provenance preserved** (never conflates 1000G vs HGDP).
  Phased datasets use exact cis co-occurrence; unphased report an assume-cis upper
  bound + a labeled lower bound.
- **Dict-less HuggingFace distribution.** `publish-index --dictless` ships the
  registry in the main index tarball + the genotype tier as a separate optional
  `genotypes_<vcf>.tar.gz`; `download` fetches the genotype tier by default
  (`--no-genotypes` to skip). Published `NRG_3_hg38-dictless+hg38_1000G_HGDP` on
  `lucapinello/crisprme-data` (the dict-based `NRG_3`/`NGG_3` indexes remain the
  default).

### Fixed
- **Corrected allele frequencies.** The AF column — empty/mis-polarized for ~95 % of
  variants in the dict format (a documented 2.2.0 limitation) — is re-derived from
  AC/AN over the full panel.
- **Out-of-the-box dict-less install.** The skip-enrichment gates
  (`submit_job_automated_new_multiple_vcfs.sh`, `validate_inputs.py`) now accept a
  `registry_<vcf>/` tier in place of `dictionaries_<vcf>/`; `download` installs a
  `-dictless`-named index under its search-resolvable canonical name and synthesizes
  the combined `<vcf>.samplesID.txt` from the per-db files. Verified end-to-end: a
  fresh dict-less download finds the CPS1 off-target (chr2:210530658, CFD 0.947,
  rs114518452) with Samples reconstructed from the genotype tier and a populated
  population summary.

## [2.2.0] - 2026-08-18

First stable release of the CRISPRme+ 2.2.0 line (Python 3.11 / Dash 2.x). It
consolidates the `2.2.0-alpha.27`–`alpha.30` pre-releases (see those sections
below for the full history); the entries here are the changes since `alpha.30`.

### Fixed
- **Phased multi-SNP haplotype off-targets are no longer under-reported.** In the
  phased post-analysis path — the default for phased datasets such as 1000G/HGDP —
  `iupac_decomposition` failed to assemble the full-haplotype off-target when one
  sample carried **≥4 co-occurring variants in a single protospacer window**: the
  true worst target (fully substituted, lowest-mismatch/highest-CFD) was silently
  under-reported or dropped. The in-loop level-0 subtraction (a dedup device that
  starved the deeper combination layers) is replaced with a **deferred peel** that
  runs after the full lattice is built, so the maximal *cis* combination forms while
  per-haplotype attribution stays deduplicated. Validated byte-identical for every
  previously-correct case (unphased, phased ≤3-variant, single-SNP, >cap greedy)
  plus a real-data genome-wide no-regression diff. (#41; a long-standing defect
  inherited from classic CRISPRme, tracked there at pinellolab/CRISPRme#175.)

### Changed
- `environment.yml` now includes `ijson`, so from-source installs get the streaming
  low-RAM SNP-dict reader that the Docker image already had (without it, SNP
  post-analysis falls back to a whole-file `json.load`; results are identical, but
  RAM use and the OOM-guard estimate are higher).

## [2.2.0-alpha.30] - 2026-08-17

### Fixed
- **Indel off-targets were silently dropped since v2.1.9.** `post_analisi_indel.sh`
  subset the per-chromosome indel targets with `grep -F -w $chrom`, but the indel
  search runs on a per-chromosome *fake* genome (`pool_search_indels.py` searches and
  names its targets `fake<chrom>`), so every indel target row's Chromosome column is
  `fake<chrom>` (e.g. `fakechr22`) — which `grep -F -w chr22` can never match (`-w`'s
  left word boundary fails because `chr22` is preceded by the word character `e`). The
  per-chromosome subsets came back empty, the indel post-analysis processed nothing,
  and **all indel off-targets were dropped with a clean exit**. Both grep lines now
  match on `"$fakechrom"` (`-w` still prevents `fakechr2` matching `fakechr22`); the
  `NF >= 10` malformed-line guard that the same refactor had dropped is restored; and
  an empty subset now prints a WARNING to stdout (stderr is fatal in this pipeline) so
  the failure can never be silent again. A hermetic regression test
  (`test_indel_chrom_subset.py`) is added and wired into CI. SNP post-analysis is
  unaffected (SNP targets use real chromosome names). Thanks to @munchr-gene1 for the
  report and diagnosis (#172).

## [2.2.0-alpha.29] - 2026-08-16

### Fixed
- **Variant search without `--annotation` no longer crashes on a read-only install
  (apptainer SIF).** `post_process.sh` wrote the annotation intermediate BEDs
  (`<annot>.tmp.bed` / `.tmp.sorted.bed`) next to the annotation file; with no
  `--annotation`/`--gene_annotation` that file is the shipped empty placeholder
  `PostProcess/vuoto.txt` inside the (read-only) install dir, so the writes failed
  and the run died at "Integrating results". Temps now go to the run's writable
  output dir. As defense-in-depth, `resultIntegrator.py` guards the closest-gene
  distance parse so a missing/empty annotation field is treated as "no annotation"
  (columns stay `NA`) instead of raising `ValueError: could not convert string to
  float: ''`. The Docker/web path (real annotation, writable FS) is byte-identical.

### Added
- **Pre-flight check for unpadded guides.** A CLI `--guide` file whose guides lack
  the PAM's N-padding (e.g. the bare 20 nt `CTAACAGTTGCTTTTATCAC` instead of
  `CTAACAGTTGCTTTTATCACNNN`) previously ran the full multi-hour search then crashed
  deep in `radar_chart_dict_generator.py` (`IndexError`). The lightweight validator
  now catches this up front with an exact fix message (append/prepend N's), handling
  both 3' (SpCas9) and 5' (Cas12a) PAM orientations. The web is unaffected (it
  auto-pads guides). Covered by 10 new `check_guide_pam_consistency` tests; the whole
  `test_validate_inputs` suite is now run in the unit-tests CI (previously it wasn't).

## [2.2.0-alpha.28] - 2026-08-16

### Fixed
- **Batteries variant search no longer aborts at pre-flight validation.** The
  lightweight input validator (`validate_inputs.run_lightweight`, always-on before
  `complete-search`) hard-required `VCFs/<dataset>/` to exist and aborted otherwise —
  but a batteries-included install (`download --what index`) ships the precomputed
  variant index + per-sample dictionaries *without* the multi-GB source VCFs. The
  validator now recognizes that case: when the VCF dataset directory is absent but the
  dictionaries (`Dictionaries/dictionaries_<dataset>/` + `log_indels_<dataset>/`) are
  present, it passes the check (source VCFs genuinely not required), mirroring the
  skip-enrichment logic in `submit_job_automated_new_multiple_vcfs.sh`. A genuinely
  missing dataset (no dictionaries either) still fails loudly. Completes the alpha.27
  batteries fix: `download --index-name NRG_3_hg38+hg38_1000G_HGDP` + a variant search
  now runs end-to-end on a fresh install. Covered by two new `run_lightweight` tests.

## [2.2.0-alpha.27] - 2026-08-16

### Fixed
- **Batteries variant search is now self-sufficient (the index-only download can run
  a variant search).** The variant post-analysis needs per-sample SNP dictionaries
  (~152 GB for 1000G+HGDP) that weren't shipped, so a fresh `download --index-name
  <variant>` install failed at genome enrichment (no source VCFs). The dicts are now
  shipped **gzipped** (152 GB → ~28 GB; indel logs 18 GB → 3.3 GB) and read **on the
  fly**, mirroring the existing bgzipped-VCF pattern — no 150 GB on-disk decompress.
  - `new_simple_analysis.py` / `analisi_indels_NNN.py`: prefer `*.json.gz` / `*.txt.gz`
    (ijson streams through gzip; json.load + line reads via `gzip.open`), fall back to
    plain — backward-compatible. Bash wrappers unchanged (Python resolves the extension).
  - `pool_search_indels.py` / `pool_post_analisi_indel.py` + the `submit_job` batteries
    chrom glob: strip an optional `.gz` from `log*.txt(.gz)`.
  - `build-index-only` gzips the per-chromosome dicts it produces (pigz-preferred).
  - Verified end-to-end by `validate-benchmarks` (build gzips → search reads `.gz`) and
    a new `test_dict_gzip` unit test (publish bundles + download routes `.gz` dicts).

## [2.2.0-alpha.26] - 2026-08-16

### Changed
- **Default SpCas9 PAM is now NRG (NAG + NGG), not NGG.** The web form and the
  "Load Example" default, plus the CLI/download docs, now default to `NRG`
  (`20bp-NRG-SpCas9`) and the `NRG_3_*` indexes. NRG matches SpCas9's broad
  recognition, so **variant-created NAG off-targets are found out of the box** —
  e.g. the CPS1 off-target from the CRISPRme paper (rs114518452: chr2:210530658,
  2 mismatches + 1 RNA bulge, variant-created `TGG` PAM, CFD 0.947). Verified
  end-to-end on a from-scratch `NRG_3_hg38+hg38_1000G_HGDP` index. `_default_pam`
  prefers NRG, falls back to NGG, then the first available PAM. The `NGG_3_*`
  indexes remain available for canonical-NGG-only searches. New `NRG_3_*` indexes
  are published on `lucapinello/crisprme-data`.

### Fixed
- **Off-target search now uses all requested cores.** `submit_job` split
  `$ncpus` in half and ran reference + variant searches in parallel; the small
  reference genome finishes first, so the long variant-search tail ran on only
  half the cores. Now runs the two searches sequentially at full `$ncpus`
  (reference 44.6 → 25.6 min in testing; `searchTST` 24 → ~45 cores).
- **Freshly-built / downloaded indexes now search correctly.** `build-index-only`
  writes `.pam_build`/`.display_label` sidecars inside the index dir, but
  CRISPRitz's `-index` search requires the dir to contain only `.bin` files
  ("only .bin files" error). Added a symlink-only `.bin` view so the sidecars no
  longer break the search, regardless of how the index was built or downloaded.

## [2.2.0-alpha.25] - 2026-08-15

### Changed
- **Contacts page tidied.** Developers are now listed with the Pinello lab first
  (Luca Pinello), then the Verona/InfOmics group (Manuel Tognon et al.), then the
  Bauer lab (Linda Lin, Daniel Bauer). Comments are directed to Luca Pinello only
  (single email), and the page now says bugs should be filed as a GitHub issue.
  Also fixed a missing space after the comma in the per-lab author lists.

## [2.2.0-alpha.24] - 2026-08-15

### Changed
- **Submission page: consolidated to exactly three font sizes (Large / Medium /
  Small).** The form previously mixed ~7 inline sizes on a 10px root, which left
  helper notes rendering *larger* than the inputs they described. It now uses one
  clear hierarchy: Large `2.2rem` for the numbered section headers ("Select
  gRNA/genome/PAM/thresholds"), Medium `1.5rem` for every interactive control
  (spacer textarea, dropdowns, Load Example / Submit buttons, the "Maximum edits"
  label and form body), and Small `1.25rem` for helper text, notes and slider
  tick marks. Verified via full-page + high-DPI screenshots — no clipping or
  overflow. Cosmetic only; no change to search behavior.
- **"Max edits" slider maximum is now derived, not a magic number.** The slider
  topped out at a hardcoded 5; it now tops out at the app's mismatch ceiling
  (`MAX_MMS - 1` = 6). Rationale: mismatches are the non-index-limited edit
  dimension and the dominant search cost, so the total-edits knob should reach as
  far as the mismatch dropdown does. The installed index caps *bulges* only
  (`bMax = N - 1`, ≤2 in practice) — that limit is applied separately where
  simple-mode bulges are derived and must not shrink the mismatch-driven total.
  Default stays 4.

## [2.2.0-alpha.23] - 2026-08-15

### Changed
- **"Max edits" is now a self-sufficient single knob — bulges are derived from it.**
  Previously a variant off-target that needs a bulge (e.g. the CRISPRme paper's CPS1
  off-target rs114518452: 2 mismatches + 1 RNA bulge) was only found if bulges were
  explicitly enabled — the CLI defaulted to bulge-free, so `--max-total-edits` alone
  found no bulge targets. Now:
  - **CLI:** when `--bDNA/--bRNA` are omitted, they are derived from `--max-total-edits`,
    capped by the bulge depth the installed index supports (so `complete-search --vcf …
    --max-total-edits 4` finds up to 2mm+2bulge patterns with no bulge flags, reusing the
    installed index — never forcing a heavier build; explicit `--bDNA/--bRNA` still win).
  - **Web:** in simple mode the per-type bulge caps now track `min(max_edits, index_cap)`
    instead of a fixed 2/2, so the slider is genuinely the only control needed.
- **Consistent default cap of 4** across web slider, Load Example, CLI `--max-total-edits`,
  and the shell fallback (was 5). Four surfaces the CPS1 off-target (needs 3 total edits)
  with headroom and is the validated performance sweet spot (cap 6 was pathologically slow
  on genome-wide combined searches).
- Docs: `--bDNA/--bRNA` are documented as optional (derived from `--max-total-edits`), not
  "Required"; the variant-aware requirement for variant-created off-targets is noted.

## [2.2.0-alpha.22] - 2026-08-15

### Added
- **History page now reports the Max edits used.** The results-history table only showed
  Mismatches / DNA bulge / RNA bulge; it now has a "Max edits" column showing the total
  mismatches+bulges cap (the simple-mode slider value, or the per-type sum with an
  "(advanced)" note). Jobs predating this field show "-".

## [2.2.0-alpha.21] - 2026-08-15

### Fixed
- **Graphical Reports: every variant off-target now shows its blue "Alternative" point.**
  The top-1000 scatter drew a blue point only when a variant had a known allele frequency;
  a variant whose MAF was missing (e.g. the frequency is not annotated from the external
  AF source) was treated like a reference-only site and its blue point was suppressed. On
  the combined 1000G+HGDP panel that hid ~93% of the variant points, so the plot looked
  "all red." Now a blue point is drawn for any off-target carried by >= 1 sample, sized by
  allele frequency when known and a default size when not; only true reference-only sites
  get no blue point. This reproduces the CRISPRme Nature Genetics figure — e.g. the CPS1
  off-target rs114518452 (chr2:210530658, CFD 0.021 -> 0.947, a variant-created off-target)
  appears as the top blue point. Also fixed in the per-sample Personal Risk Card plots.

## [2.2.0-alpha.20] - 2026-08-15

### Added
- **Graphical Reports: a "top 1000 by variant effect" plot.** Alongside the existing
  score-ranked top-1000 scatter, the CFD and CRISTA reports now also show a companion
  plot ranked by the size of the variant-induced score change (`|ALT - REF|`), so the
  variants that actually change the off-target score are foregrounded instead of being
  buried among the many that leave it unchanged. The score-ranked plot is unchanged; the
  fewest-mm+b report stays score-only. (`CRISPRme_plots.py`, shown in `results_page.py`.)

## [2.2.0-alpha.19] - 2026-08-15

### Fixed
- **Graphical Reports: the top-1000 scatter now shows the variant (ALT) points.** The blue
  ALT markers were drawn ~10x smaller than the red REF markers at the same x, so they were
  swallowed and the plot looked "all red." They now have a floored size + a black edge (and
  stay hidden for reference-only sites), so a variant that changes the CFD/CRISTA score is
  visible at the tip of its arrow. (`CRISPRme_plots.py`, `CRISPRme_plots_personal.py`.)
- **Graphical Reports: the CRISTA and fewest-mm+b top-1000 plots are generated again.** The
  plot generator sorted the CRISTA plot on a non-existent column
  (`CFD_score_(highest_CRISTA)`) and crashed with a silently-swallowed `KeyError` right after
  the CFD image, so only the CFD plot appeared. Fixed to `CRISTA_score_(highest_CRISTA)`.
- **Reference-only searches no longer break the Graphical Reports tab.** The population-barplot
  callback `Output` targeted a container the layout omitted for reference-only searches, which
  could fail the whole callback (radar chart included). The container is now always present
  (hidden for reference-only).
- **Personal Risk Card now shows the candidate sites.** The detail table was built only from
  *private* targets (unique to the sample), so it was empty for the common case of a sample
  with shared-but-not-private candidates. It now shows both a **Personal candidates** table
  (all off-targets the sample carries) and the **Private** table.
- `update_images_tabs` now guards an empty table selection (no crash on initial render) and
  applies the non-SpCas9 -> "fewest" criterion remap the rest of the page uses; fixed a radar
  chart `except` handler that assigned an unused variable.

## [2.2.0-alpha.18] - 2026-08-15

### Fixed
- **INDEL post-analysis: forward-ported the boundary/dense-variant crash guards.** The
  alpha.17 CRISTA window guard was only in the SNP path (`new_simple_analysis.py`); the
  parallel INDEL path (`analisi_indels_NNN.py`) still assumed a full 29 nt window, so an
  out-of-range indel coordinate could crash post-analysis. The `!= 29` null-guard is now
  mirrored in both INDEL CRISTA blocks.
- **INDEL CFD scoring no longer KeyErrors on an unexpected PAM.** `analisi_indels_NNN.py`
  used a raw `pam_scores[pam]` where the SNP path uses the guarded `pam_scores.get(pam,
  0.0)`; now matched (issue-#94 class).
- **CLI `--sorting-criteria` default now matches its own `--help` and the web.** It
  defaulted to `mm+bulges` but documented (and the web uses) `mm+bulges,mm`; the
  fewest/total tie-break now ranks tied targets identically on both paths.

### Changed
- **CLI and web now share the same default total-edits cap (5).** The web "Maximum edits"
  slider defaulted to 3 while the CLI `--max-total-edits` defaults to 5; the web slider
  (and Load Example) now default to 5 so a default web search and a default CLI search
  cover the same edit space.

## [2.2.0-alpha.17] - 2026-08-14

### Fixed
- **Post-analysis no longer crashes on variant-dense / chromosome-boundary targets.**
  On the combined 1000G+HGDP index the example search (and any dense-variant search)
  could abort in post-analysis with a `ZeroDivisionError` (CRISTA scoring divided by the
  length of a 29 nt window that had been stripped empty of IUPAC codes) or an
  `IndexError` in `iupac_decomposition` (a reference window truncated at a chromosome
  boundary was shorter than the target). Both came from assuming a genome slice is always
  full length. Now an un-scoreable CRISTA window is nulled (score −1, as for windows with
  `N`) and IUPAC positions past a truncated reference are skipped, with belt-and-suspenders
  divide-by-zero guards in `CRISTA_score.py`. Verified end-to-end on the batteries
  1000G+HGDP data: the example completes and its Personal Risk Card renders.

### Changed
- **Base editing is off by default after "Load Example."** The example now loads a
  standard (non-base-editing) search; the base-editor toggle stays on "No".
- **Larger helper text on the submission form.** The small gray guidance texts in the
  submission cards were ~10.5 px (the page root font is 10 px); enlarged to ~14 px so
  they are readable. The page's base font is unchanged.

## [2.2.0-alpha.16] - 2026-08-14

### Fixed
- **Result title is now mode-aware and no longer self-contradictory.** In simple mode
  the title listed the per-type mismatch/DNA/RNA caps *and* the total "Max edits" cap,
  which contradict when the total cap is tighter (e.g. `Mismatches 4 - DNA bulges 1 -
  RNA bulges 1 - Max edits 1`). Now a simple-mode search shows only `Max edits
  (mismatches + bulges) N` (the single control the user set), while an advanced search
  shows the explicit per-type caps. `change_url` records `Threshold_mode` (simple/
  advanced) in `.Params.txt`; older jobs are inferred (a total cap tighter than the
  per-type sum ⇒ simple) and pre-alpha.15 jobs fall back to the per-type title.
- **Dedup reuse wrote a malformed completion marker.** When a resubmission reused a
  prior clean result, the marker was written `Job\nDone` (newline) instead of
  `Job\tDone` (tab), so the dedup check and the status page's completion detection
  would not recognize the reused job as finished. Fixed to write `Job\tDone`.

## [2.2.0-alpha.15] - 2026-08-14

### Fixed
- **The result title now shows the "Max edits" cap that actually bound the search.**
  In simple mode the "Maximum edits (mismatches + bulges)" slider is the binding
  constraint, but the Result Summary title only listed the per-type mismatch/DNA/RNA
  caps — which can be looser than the total cap (and, after "Load Example", are values
  the user never explicitly chose), so the title overstated how broad the search was.
  The total-edits cap is now recorded in `.Params.txt` (`Max_total_edits`) and appended
  to the title as "Max edits (mismatches + bulges) N". Backward-compatible: results for
  jobs created before this field simply omit the clause.

### Changed
- **"Load Example" now sets the Max-edits slider (to 4) explicitly** instead of
  leaving it at the default, so the example is driven by the same total-edits cap the
  report displays and the per-type caps it fills in (4/1/1) are consistent with it.
- Docs: added mentions of the optional **email notification** (README, Docker
  quickstart, landing pages) and clarified that **Personal Risk Cards** appear only
  for variant-aware searches (landing pages, in-app help, quickstart, user guide).

## [2.2.0-alpha.14] - 2026-08-14

### Fixed
- **A resubmitted search no longer resurfaces a previously-FAILED job.** Results are
  deduplicated by search parameters; the dedup reused any matching prior job, so
  re-running a search whose earlier attempt had failed (e.g. left behind by a bug
  since fixed) returned that stale failure — surfaced in the UI as "The selected
  result encountered some errors, please remove it and try to submit again." Dedup
  now skips prior results that did not finish cleanly (non-empty `log_error.txt` or
  no `Job Done`) and runs the search fresh instead.

### Changed
- Docs: dropped the stale "peak memory is in post-processing and can exceed
  32 GB / ~64–100 GB" guidance (the alpha.12 streamed SNP-dict load removed that
  post-analysis spike); the 16 / 32 / 64 GB tiers stay.

## [2.2.0-alpha.13] - 2026-08-14

### Changed
- **Post-analysis performance** (follow-up to the streaming dict-load fix):
  - Added the `yajl` C library + `cffi` so `ijson` uses its fast `yajl2_cffi` C
    backend instead of the pure-python fallback (much faster streaming of the SNP
    dictionary).
  - The post-analysis worker cap no longer sizes per-worker RAM from the (large) SNP
    dictionary file: with streaming, per-worker peak is small, so it uses a modest
    fixed estimate and runs chromosomes in parallel again (the SNP path checks that
    `ijson` is importable; the INDEL path never loaded the big dict). Previously it
    over-estimated and serialized to a single worker.

## [2.2.0-alpha.12] - 2026-08-14

### Fixed
- **On-demand INDELS index built with a doubled `<PAM>_<N>` prefix** (regression
  from the alpha.10 indels-index-resolution change). crispritz `index-genome`
  already auto-prepends `<PAM>_<bMax>_` to the output dir, but the builder also
  added it manually, producing `NGG_2_NGG_2_..._INDELS`; detection/search look for
  the single-prefix name, so the indel search produced no output and the run failed
  with "off-targets search failed". Removed the manual prefix (let crispritz add
  it). This unblocks `validate-benchmarks` (complete-test + validate-test now match
  the brute-force ground truth: 2/2 light benchmarks). Batteries searches were
  unaffected (they use a correctly-named pre-built `_INDELS`).
- **Genome-wide variant post-analysis OOM ("Killed … EmptyDataError").** On a
  genome-wide 1000G+HGDP search, each per-chromosome worker `json.load()`ed the
  whole chromosome SNP dictionary — a chr2 dict is ~13 GB on disk and ~26 GB in
  RAM — so even one worker could exhaust a normal machine (and the concurrency
  guard defaulted to a fixed 64 GB budget / 4 GB-per-worker, over-subscribing on
  smaller hosts). `new_simple_analysis.py` now loads **only the dictionary entries
  the search's off-targets actually query** (streamed with `ijson`), dropping peak
  RAM from ~26 GB to <0.1 GB on the real chr2 dict (byte-identical results). The
  worker-cap guards now **detect the machine's real RAM** (host + cgroup limit),
  **estimate per-worker RAM from the actual dictionary sizes**, and print a clear
  "needs ~N GB / may be OOM-killed" warning instead of the cryptic cascade. Added
  `ijson` to the image.

### Changed
- Docs: corrected the in-app manual and web guide (removed the stale "up to 100
  spacers / first-100-sequences" copy; annotation is Settings-managed not a Step-3
  dropdown; email is configured in Settings, not "the server notifies").

## [2.2.0-alpha.11] - 2026-08-14

### Fixed
- **Graphical Reports population barplot never rendered.** The pipeline writes
  `populations_distribution_<guide>_<N>total_<score>_new.png`, but the results
  page read the name without the `_new` suffix, so every mismatch/bulge combo
  showed "No result found for this combination". Added the missing suffix.
- **Graphical Reports defaulted to an empty combo.** The "up to N edits" selector
  defaulted to `0` (on-target only), which is empty for typical sparse guides. It
  now defaults to the full budget (mm + bulges) — the cumulative barplot the
  pipeline only emits when it has data — so a plot shows immediately.
- **"Load Example" left the base-editing window blank.** The `be-window-start/stop`
  dropdowns start with no options, so the example's window values were dropped
  before `update_base_editing_dropdown` populated the options. Load Example now
  emits the window options atomically with the values (and passes ints, matching
  the option type), so it fills 4/8 correctly.

## [2.2.0-alpha.10] - 2026-08-14

### Added
- **Email notifications are now configurable in the app** (Settings → Email
  notifications, local mode): SMTP host/port/SSL, sender, and an app password are
  saved to a local `.email.json` (0600, password never echoed back to the
  browser). `send_mail.py` reads that config (falling back to the
  `CRISPRME_SMTP_*` env vars) and is now best-effort — it never writes to stderr
  or raises, so a mail-server problem can't fail an already-completed search.

### Changed
- **Search form polish**: the CRISPRme+ logo now sits above the first step; the
  *Load Example* button moved to the bottom of the *Select gRNA* box; larger,
  more legible fonts throughout the form; the spacer help text explains that
  multiple spacers can be entered (and points to the CLI for bulk runs); the
  *Notify me by email* box explains the one-time Settings setup.

### Fixed
- **Variant indel indexing/search now resolves the INDELS companion from the
  actually-resolved variant index** instead of reconstructing an exact
  `<true_pam>_<bMax>_..._INDELS` name. A batteries/precomputed index whose N or
  PAM differs from the request (e.g. the shipped `NGG_3` served against a
  `bMax=2` search, or a pamless `NNN` index) was not matched, so the pipeline
  attempted an impossible rebuild and crashed in `pool_index_indels.py`
  (`FileNotFoundError: .../Genomes/<ref>+<vcf>_INDELS/`). The INDELS index, its
  search path, and the on-demand builder now all key off the index's own
  `<PAM>_<N>` prefix (consistent with `build-index-only`/`publish`/`download`),
  so newly built-and-uploaded indexes work generally. When the companion is
  genuinely absent and cannot be rebuilt (no source genome), the run now fails
  with a clear "re-download the index" message instead of a raw traceback.

### Added
- **Annotation manager** in Settings → Data Manager: a per-genome
  *Manage annotations (enable / disable)* checklist backed by a persisted
  `.enabled.json` manifest, plus format validation on every uploaded annotation
  BED (≥4 tab columns, integer `start ≤ end`, whitespace-free label, size cap;
  malformed files are rejected with a clear message and nothing is written).
  Enabled tracks are merged on demand into a single cached active annotation and
  applied automatically to every search.

### Changed
- Annotations are now applied automatically from the enabled set rather than
  chosen per search. The annotation dropdown was removed from the search form
  (Step 3 is now just email + job name); the built-in
  ENCODE cCREs (SCREEN) + DHS + GENCODE bundle is enabled by default so
  annotations work out of the box.
- The **Maximum edits** slider is floored at 1 (the on-target / 0-edit hit is
  always reported at any setting), avoiding an empty-result edge case on large
  precomputed indexes.
- Search-submission form: consistent left indentation for the PAM and
  variant-dataset dropdowns and the submit/example buttons.

## [2.2.0] - 2026-08-06

### Added
- Graphical **Settings / Data Manager** page in the web interface: add reference
  genomes (UCSC by assembly name — e.g. the pig `susScr11` — HuggingFace, or a
  direct URL), precomputed indexes (download from HuggingFace or build locally
  from an installed genome + PAM), VCF datasets (HuggingFace or register an
  existing server folder), annotations (BED upload), and nucleases/PAMs (a small
  form). New data lands in the local data folder and is auto-discovered by the
  search form. Long operations run as detached jobs on a dedicated executor with
  live progress, so they never starve the search slots. Mutations are local-mode
  only; publishing an index to the shared HuggingFace repo is maintainer-only.
  `download --what genome` gained `--source {hf,ucsc,url}` (+ `--url`) so the CLI
  and web share one non-human-genome download path.
- Python 3.11 modernization: pipeline fixes for pandas 2.x / matplotlib 3.x and
  a Dash 1.x → 2.x web-app migration, plus a Python-3.11 Docker image built from
  source (CRISPRitz 2.8.1) (#131).
- `assembly-search` subcommand: off-target search on a personal diploid genome
  assembly (two haplotypes, no VCF), reconciled to hg38 via liftOver (#113).
- Reference-index UX: `build-index-only` pre-builds the reusable CRISPRitz
  reference index without running a search, and `complete-search --index-path`
  reuses a prebuilt/staged index library (a missing index is a hard error rather
  than a silent rebuild).
- HuggingFace data distribution: `download` fetches reference data (genome,
  annotations, PAMs, sample IDs, VCFs, precomputed indexes) from a HuggingFace
  dataset repository over its CDN, and `publish-index` uploads a locally built
  index for reuse. Default repo `lucapinello/crisprme-data`, overridable via
  `--hf-repo` / `CRISPRME_HF_REPO`. `setup`/`complete-test` also try HuggingFace
  first and fall back transparently to the original UCSC/EBI/Sanger sources
  (#140, #141).
- `complete-search --max-total-edits N`: cap the total edits (mismatches +
  bulges) per reported alignment; over-cap targets are dropped right after the
  search, shrinking intermediate files and post-analysis time (#107).

### CI
- New `unit tests` workflow: fast, hermetic byte-compile + network-free HF/index
  unit tests on every code PR.
- New `web e2e (playwright)` workflow: builds the py3.11 image, serves the web
  app, and drives Chromium to assert every Dash 2.x page renders (no blank pages
  / JS errors).
- `validate-benchmarks` gained a new-subcommand dispatch + unit-test smoke step.

### Changed
- Clearer failure reporting: when a search fails, CRISPRme now prints *which
  stage* failed (from the per-stage log) and the last lines of the error log,
  instead of only "run failed — see log_error.txt". Makes failures actionable
  for non-expert users.

### Fixed
- Web interface (Dash 2.x) hardening, from a full Playwright stress test of the
  running app:
  - The web server no longer crashes on a from-source install. Dash 2.x's
    `app.run()` lets the `HOST` environment variable override the host argument,
    and the from-source conda env sets `HOST` to a non-bindable compiler build
    triple (`x86_64-conda-linux-gnu`); the server now forces the intended host/
    port so it binds correctly.
  - The **Query Genomic Region** and **Personal Risk Cards** result tabs no
    longer return HTTP 500. Both callbacks type-checked their inputs before the
    "no click yet" guard, so Dash's initial (empty) render raised a `TypeError`;
    the guard now runs first, and Filter/Generate with nothing selected is a
    graceful no-op.
  - Removed the dead cross-origin "skeleton" stylesheet (blocked by browsers on
    every page); the layout already uses the Bootstrap grid.
  - The nuclease dropdown collapses case-variant duplicate PAM files so each
    nuclease is listed once.
- Zero-hit searches now complete cleanly with an empty result instead of
  aborting. A search that finds no off-targets (e.g. a very stringent
  guide/parameter combination) previously failed part-way through post-analysis
  ("off-targets post-analysis (reference) failed", then a cascade through the
  rsID / summary / integration steps, all of which assumed at least one target).
  Added a zero-target guard to the reference SNP post-analysis (mirroring the
  INDELs one), made `remove_n_and_dots.py` tolerate a header-only report, and
  added a high-level short-circuit that emits an empty-but-valid result and exits
  0 when no off-targets are found. Verified end-to-end on ml007 (empty result
  exits 0; a normal with-hits search is unaffected).

### Documentation
- New `docs/DOCKER_QUICKSTART.md` + a README quickstart callout: a few-command,
  no-conda / no-410 GB path to the web interface for non-experts — fast HF data
  download, a prebuilt index, then `web-interface` in the browser, with
  troubleshooting and "install more indexes as needed".
- Data-setup guide: documented what `setup` produces (including the combined
  1000G+HGDP config files), the HuggingFace fast-download path, and a new
  "Prebuild, reuse, and share the reference index" section.
- README: added a from-source install path and a reference-index /
  data-distribution commands section.
- Documented the variant PAM behaviour in `docs/INPUT_FORMATS.md`: PAM *creation*
  (variants that add an off-target) is reported; PAM *disruption* (variants that
  remove a reference PAM for carriers) is a known, intentional gap — with the
  rationale (disruption lowers rather than raises predicted risk and would need
  sample-specific, non-scorable semantics).

## [2.1.13] - 2026-08-04

### Added
- Robust `complete-test` downloads: HTTP retry with linear backoff plus
  checksum-verified resume, so a dropped connection or truncated transfer from a
  slow/flaky host retries (and skips already-verified files) instead of aborting
  a multi-hour run (#136).
- Memory-bounded post-analysis: per-chromosome SNP/indel workers are capped to a
  RAM budget (default 64 GB, overridable via `CRISPRME_MAX_MEM_GB` and
  `CRISPRME_POSTPROC_WORKER_GB`), preventing the peak-RAM spike observed on
  genome-wide 1000G runs (#136).

### Fixed
- Post-analysis worker-count diagnostic is now written to stdout (`log_verbose`)
  instead of stderr. The caller treats a non-empty stderr log as a fatal
  post-analysis failure, so the informational message previously aborted every
  run with a false "post-analysis (snps) failed" even when the analysis
  succeeded (#136).

## [2.1.12] - 2026-08-01

### Added
- Cas9 and Cas12a validation benchmark examples with precomputed brute-force
  ground-truth references and an extensible registry
  (`test/benchmark/benchmarks.json`); `validate-test` compares CRISPRme
  off-targets against them (#116).
- A fast, dependency-free Rust port of the brute-force ground-truth generator
  (`test/benchmark/rust/`) that produces output identical to the Python
  generator (#116).
- `validate-benchmarks` GitHub Actions CI that runs the full
  `complete-test` → `validate-test` round-trip on Linux for every registered
  benchmark (#116, #120).
- Native Apple Silicon (`linux/arm64`) Docker image via a multi-arch `buildx`
  workflow; the multi-arch image is published to Docker Hub on release tags
  (#121).
- Low-memory startup warning (`PostProcess/memory_check.py`) advising Docker
  Desktop users to raise the memory limit when total RAM is below ~32 GB (#121).
- Search-space budget estimator and guard (`PostProcess/search_budget.py`) that
  warns before resource-explosive runs, with a companion analysis in
  `docs/SCALABILITY_ANALYSIS.md` (#118).
- Release tooling: a `release-crisprme` Claude Code skill, `docs/RELEASING.md`,
  this `CHANGELOG.md`, and `scripts/prepare_release.py` (#117); plus docs for the
  fast Rust generator and using the release skill (#122).
- Input-format hardening: `--gene_annotation` is now auto-sorted like
  `--annotation`; clear errors for multi-line PAM files and malformed PAM
  filenames; and a non-fatal warning when a degenerate PAM motif is combined with
  bulges (the CRISPRitz #105 crash, fixed in CRISPRitz 2.7.1) (#125).
- Documentation: `docs/INPUT_FORMATS.md` (PAM/Cas12a/VCF/chromosome/annotation/
  bulge guidance, #124), `docs/ISSUE_AUDIT.md` (a review of all open + closed
  issues, #126), and `docs/ROADMAP.md` (the 2.1.12 → 2.1.13 plan, #128).

### Fixed
- Multiple bugs affecting new VCF dataset processing, concurrent runs, and
  gnomAD handling (#96).
- Guard the CFD PAM-score lookup against non-ATCG PAM bases present in real hg38,
  which previously crashed the run with a `KeyError` (#94, #125).
- Guard the post-analysis intermediate reads against non-UTF-8 bytes
  (`UnicodeDecodeError`) so a stray byte no longer aborts a run (#52, #125).
- `validate-test` now exits non-zero when a benchmark mismatches or the
  `complete-test` output is missing, instead of passing silently (#116).
- `complete-test` downloads 1000 Genomes VCFs over HTTPS instead of FTP,
  unblocking FTP-restricted and CI networks (#116).

### Changed
- Synced the `Dockerfile` `crisprme_version` pin to the released version, which
  had lagged behind `crisprme.py` and the Bioconda recipe (#117).

## [2.1.11] - 2026-07-06

### Fixed
- Consistent `.gz` suffix for compressed annotation files, removing a mismatch
  between file content and file name in the standard and personal-annotation
  pipelines (`crisprme.py`).
- `_check_personal_annotation` now strips a `.gz` suffix before processing so
  downstream steps see the expected filename.
- Removed double-compression of annotation files, preventing corrupted or
  unreadable annotation outputs.

### Changed
- Refactored per-chromosome test genome-directory handling: replaced
  `ensure_hg38_directory` with `_assign_genome_directory_name` (whole-genome
  builds in `<dest>/hg38`, single-chromosome builds in `<dest>/hg38_<chrom>`).
- `download_genome_data` now returns the resolved genome directory path directly.
- Renamed downloaded sample-ID files to include genome/dataset context
  (`hg38_{dataset}.samplesID.txt`) via a centralized mapping.
- Reduced default web-interface example-run parameters (max mismatches 6 → 4,
  DNA bulges 2 → 1, RNA bulges 2 → 1) for a faster example run.

### Added
- New `_write_sg1617_file` helper to generate an `sg1617` guide fixture during
  legacy database setup.
- New guides `docs/crisprme_data_setup_051826.md` (CLI setup/usage) and
  `docs/crisprme_web_interface_user_guide.md` (web interface walkthrough).
- "Setup Legacy Database" section in `README.md`.

### Removed
- Unused `_bgzip_ann_data` helper and stale typing imports/docstrings.

## [2.1.10] - 2026-05-29

### Added
- New `setup` command for automated download/installation of CRISPRme reference
  resources directly from the command line, removing the need for manual
  downloads from the CRISPRme website.

### Changed
- Updated documentation and container configuration to support the new
  installation workflow and improve reproducibility across deployments.

## [2.1.9] - 2026-01-16

### Added
- `validate-test` functionality: validates `complete-test` off-target
  predictions against brute-force ground-truth alignments derived from
  1000 Genomes variant data (PR #92).
- Chromosome-level validation support (single chromosome or genome-wide).

### Changed
- Strengthened the `complete-test` pipeline and refined benchmarking utilities
  for correctness and reproducibility with population-scale variant data.

## [2.1.8] - 2025-12-10

### Fixed
- Targeted fixes and stability improvements across annotation handling,
  off-target search execution, test coverage, and report generation.
- Fixed and updated the `complete-test` workflow to align with the latest
  internal logic and outputs (PR #86).

### Changed
- Updated Docker distribution and testing environment for consistent runtime
  behavior (PR #85).

## [2.1.7] - 2025-03-10

### Added
- Support for exome-based gnomAD VCF files.

### Fixed
- Corrected handling of indels in off-target reporting (previously missing or
  misreported).
- Fixed the visualization of the Personal Card tab in the GUI and website.
- Closed file streams after creating graphical reports to prevent memory issues.

### Changed
- Improved contiguous-target merge logic (always reports the leftmost off-target
  in a cluster).
- Improved error-tracing logic.

## [2.1.6] - 2024-11-27

### Added
- Support for gnomAD v4.1 VCFs, including joint variant files.
- New `complete-test` function for quick single-chromosome and full-genome
  off-target testing.

### Fixed
- Corrected indel handling in off-target reporting to eliminate false positives.
- Fixed the alternative-alignments report to list all possible alignments.
- Resolved Personal Card tab visualization issues in the GUI and website.
- Addressed a memory-overflow issue when merging contiguous targets.

### Changed
- Upgraded the DockerHub image with the latest fixes.

[Unreleased]: https://github.com/pinellolab/crisprme-plus/compare/v2.3.0...HEAD
[2.3.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.0
[2.2.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.2.0
[2.2.0-alpha.30]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.2.0-alpha.30
[2.1.13]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.13
[2.1.12]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.12
[2.1.11]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.11
[2.1.10]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.10
[2.1.9]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.9
[2.1.8]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.8
[2.1.7]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.7
[2.1.6]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.6
