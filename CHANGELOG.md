# Changelog

All notable changes to CRISPRme are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When cutting a release, move items out of `[Unreleased]` into a new dated
version section and update the link-reference footer. See `docs/RELEASING.md`
and the `release-crisprme` skill.

## [Unreleased]

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

[Unreleased]: https://github.com/pinellolab/crisprme-plus/compare/v2.1.13...HEAD
[2.1.13]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.13
[2.1.12]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.12
[2.1.11]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.11
[2.1.10]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.10
[2.1.9]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.9
[2.1.8]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.8
[2.1.7]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.7
[2.1.6]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.6
