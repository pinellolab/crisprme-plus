# CRISPRme+ release checklist

A living, run-every-release gate. It has two halves:

1. **Gates** — what must be true before a new `vX.Y.Z` tag.
2. **Regression traps** — concrete things we have gotten wrong before, each paired with the
   check that catches it. **Append to this list whenever a release pass finds a new problem**,
   so the same class of bug can never ship twice.

> Versioning is **stingy** — stay in the `2.x` line (patch bumps) unless the change is a true
> CRISPRme-3 break. Never cut a tag on anything less than all-green CI.

---

## 1. Gates

### A. Release-surface consistency (the invariant)
Every user-facing surface must agree on **version · default index · scorer · licence · analysis-mode · data-layout**. The surfaces:
`README.md` · `METHODS.md` · everything under `docs/` **including the GitHub-Pages website `docs/index.html`** and the `.md` pages it links · installation (`README` Conda/Docker/source, `docs/DOCKER_QUICKSTART.md`, `install_from_source.sh`, `Dockerfile`, `environment.yml`, the Bioconda recipe) · **the Hugging Face dataset card** (`lucapinello/crisprme-data` README).

- [ ] **Version sweep** — `crisprme.py` `version`, `Dockerfile` (CRISPRitz ARG), the tag/image `vX.Y.Z`, and every version string in README + `docs/` read the SAME `X.Y.Z`. Grep gate: no stale `2.5.x/2.6.x/2.7.x` strings or `v2.6.x` image tags linger outside historical/CHANGELOG/archived-doc context.
- [ ] **Default index** — every surface names `NRG_3_hg38+hg38_HGDP1kGP` (SHAPEIT5 re-phased, 4,091 = 925 HGDP + 3,166 1000 Genomes) as the recommended **variant** default; `NRG_3_hg38` is the reference-only default. No surface proposes a retired/alternative index as *the* default.
- [ ] **No retired index advertised** — grep the whole repo (code help strings, scripts, tests, docs) + the HF card for retired names; they may appear ONLY in a past-tense/retirement note. Retired set: `NRG_3_hg38+hg38_1000G2021_HGDP` (hybrid), any `*-dictless*`, `NGG_3_hg38`, `genotypes_hg38_1000G_HGDP`.
- [ ] **Scorer** — CFD + CRISPR-Bulge everywhere; zero mentions of CRISTA or azimuth (removed).
- [ ] **Licence** — MGBOA 1.0 everywhere (README, METHODS, `LICENSE`, HF card frontmatter `license: other` + `license_name`/`license_link`, report footer). The migration table states commercial use needs a separate licence for **both** CRISPRme+ and stock CRISPRme.
- [ ] **Analysis mode** — population-level default vs opt-in `--per-sample`; zero mentions of the removed `--fast`/`--full`.
- [ ] **HF catalog matches reality** — `crisprme_hf.list_available_downloads` / the Settings dropdowns / the docs list EXACTLY the indexes present on HF (`list_repo_files` on `lucapinello/crisprme-data`). Nothing listed that 404s; nothing on HF that is undocumented.
- [ ] **Every documented `download --what X` example resolves on HF** — for each `download --what {genome,index,vcf,annotations,pams,samples,assembly,all}` command shown in ANY doc (README, DOCKER_QUICKSTART, index.html, DATA_LAYOUT, SETTINGS, web guide, legacy setup), confirm the HF path it fetches actually exists (`list_repo_files`). A documented command that targets unhosted data is a broken step (see trap T-21: `--what vcf --dataset 1000G|HGDP` pointed at a `vcfs/` folder that was never uploaded).

### B. Functional gates (must pass)
- [ ] **Clean-room new-user e2e** on current `main` (overlay into the latest SIF, or a fresh Docker pull): `download --what all` → a real search → report renders with closest-gene columns + CONFIRMED cis + named carriers from HGDP1kGP. The gold-standard gate.
- [ ] **Full CI green on the release PR** — `build`, `unit-tests`, `validate-benchmarks`, `web-e2e`, `Build multi-arch image`. ⚠️ A PR touching only `pages/**` or docs does **not** trigger `unit-tests`/`validate-benchmarks` (path filters — see trap T-09); the version bump touches `crisprme.py`, which re-triggers the full suite, so confirm they actually ran.
- [ ] **Web walkthrough (Playwright)** — home search submit → `/load`; Settings download index / add annotation / download genome / build index / delete genome; assembly-search submit. Screenshots archived.
- [ ] **CLI test suite** — `PostProcess/test_*.py` + `test/web/*.py` pass inside the shipped env.

### C. Release mechanics (at/after tag)
- [ ] Bump `crisprme.py` version; move `CHANGELOG [Unreleased]` → `[X.Y.Z]` dated + fresh `[Unreleased]` + fix link footer.
- [ ] Tag → `docker-multiarch` CI publishes `pinellolab/crisprme:vX.Y.Z` + `:latest`; verify `crisprme.py --version` == `X.Y.Z` and CRISPRitz `v2.9.0` inside the image.
- [ ] Bioconda: refresh the recipe/Gmail draft + `sha256` of the new release tarball; README's bioconda-status line matches reality.
- [ ] HF dataset card re-checked against the shipped catalog (part of surface A).
- [ ] Ping Ann (@anncir1) to review any web/assembly-search changes on the PR.

---

## 2. Regression traps (things we got wrong before → the check that catches it)

| # | What bit us | The check that catches it |
|---|---|---|
| T-01 | Version strings drifted across surfaces (v2.7.8: README/docs still said 2.6.2 / 2.7.1 / 2.5.5 while shipped was 2.7.7). | Gate A "version sweep" grep, every release. |
| T-02 | Docs/scripts/tests kept advertising indexes **after they were deleted from HF** (hybrid `1000G2021_HGDP`, `*-dictless`, `NGG_3_hg38`). | Gate A "no retired index advertised" repo-wide + HF-card grep. |
| T-03 | The **HF dataset card** was never part of the consistency sweep → it still declared AGPL, a fictional index layout, and the wrong `pinellolab/` repo-id. | HF card is now an explicit surface in Gate A. |
| T-04 | Web **first-paint** seeded the variant default via `_preferred_variant()` **without the genome**, so the phased-preference was skipped and it could pick a non-HGDP1kGP panel. | Both call sites pass the genome; Gate A "default index" + a web walkthrough confirming the proposed default. |
| T-05 | `crisprme.py` path-correction hard-stripped a 3-char `"bin"` suffix → the `opt/` copy spawned by the web got a garbage PostProcess path (`ModuleNotFoundError: validate_inputs`), breaking **every** web assembly-search. | Layout-detection (PostProcess-dir probe) instead of suffix math; web assembly-search e2e in the walkthrough. |
| T-06 | Hidden web dropdowns used `dcc.Dropdown` whose `value` didn't round-trip → assembly-search submitted `None` and no job was submittable. | Use `dcc.Store` for hidden state; assembly-search submit e2e. |
| T-07 | Web search didn't pass `--gene_annotation` → the gencode companion was discarded, so **closest-gene columns were missing** from web jobs. | Closest-gene columns asserted in the web-job report during the walkthrough. |
| T-08 | A single multi-chromosome custom VCF produced a cryptic `Genomes/<g>/.fa` failure (the enricher needs per-chrom files); a SNP-only VCF aborted the whole build at the indel step. | `_ensure_per_chromosome_vcfs` auto-split + graceful 0-indel skip; custom-VCF build test. |
| T-09 | CI **path filters** mean a `pages/**`-only or docs-only PR silently skips `unit-tests` + `validate-benchmarks` — a PR can look "green" without running them. | Gate B note: confirm the full suite actually *ran* (the version-bump PR touches `crisprme.py`, which re-triggers it). |
| T-10 | `apptainer exec` bypasses the image entrypoint, so `crisprme.py` isn't on PATH. | Use `apptainer run` (or call `/opt/conda/bin/python` directly) in all HPC/test invocations + docs. |
| T-11 | HF uploads/deletes 401 when the token comes from the env; `HfFolder` was removed in new `huggingface_hub`. | Pass `token=get_token()` explicitly; delete via `create_commit([CommitOperationDelete…])`. |
| T-12 | Editing a bash script **while instances of it are running** corrupts the running jobs (broke chr7/chr9 in a GW re-phase). | Never edit a running script — copy to a new filename and run that. |
| T-13 | A fresh `complete-search` with an empty index-path builds the search index but **not** the `Dictionaries/` tiers → variant/indel post-analysis crashes. | New-user path = `download` (or `build-index-only`) FIRST, then `complete-search --index-path`; documented in the quickstart. |
| T-14 | `samplesID` files that over-list samples (vs VCF-genotyped) silently **deflate every reported AF**. | Pre-flight `check_samplesid_overlisting` WARN; AN denominator sanity check. |
| T-15 | Indel off-targets were silently dropped by a `grep -w $chrom` vs `fakechr` mismatch in `post_analisi_indel.sh`. | `grep -F -w "$fakechrom"`; indel off-targets asserted present in the e2e. |
| T-16 | Read-only `vuoto.txt` in the SIF made an **unconditional** `truncate` fail and abort the pipeline. | Guard writes/`truncate` on existence + writability before touching a file. |
| T-17 | A test run via `apptainer exec --cleanenv` + explicit `/opt/conda/bin/python` gives a minimal PATH that EXCLUDES `/opt/conda/bin`, so external-binary-dependent tests (BEDOPS `sort-bed`/`closest-features`, `bcftools`) spuriously "fail" even though the binary ships in the image. | When running the suite out-of-entrypoint, pass `--env PATH=/opt/conda/bin:/usr/bin:/bin` (or use `apptainer run`); confirm a reported missing-binary failure isn't just PATH before treating it as a release bug. |
| T-18 | The same fact stated in multiple README/website sections drifted apart — the default-index sample counts disagreed across sections (`929/3162` vs the correct `925/3166`; a stale `3,477` retired-hybrid panel size in the report-summary example). | The section-by-section docs/website review (see Gate B) + a grep that every default-index count/panel-size agrees; single-source the number and point other sections at it. |
| T-19 | A public install page handed users the WRONG package's subcommands — the website Conda path told pre-2.2.0 stock-Bioconda users to run CRISPRme+-only `download`/`web-interface`, which that package lacks → the documented step errors. | Every install path's example commands must match the subcommands the package it just installed actually has (stock Bioconda = legacy `setup`; CRISPRme+ Docker/source = `download`/`web-interface`). |
| T-20 | Overlay testing bound `crisprme.py` to `/opt/conda/opt/crisprme/crisprme.py`, but the EXECUTED wrapper is `/opt/conda/bin/crisprme.py` (a COPY, not a symlink) — so current-main `crisprme.py` wrapper changes (help text, arg defaults, new subcommands) were silently NOT exercised; tests ran the image's old wrapper. | To test current-main `crisprme.py` in a SIF, bind BOTH `--bind <src>/crisprme.py:/opt/conda/bin/crisprme.py` AND the `opt/crisprme/{PostProcess,pages}` trees. (Modules/pages load from `opt/`; the entrypoint runs the `bin/` copy.) |
| T-21 | A documented download command pointed at data NOT hosted on HF — `download --what vcf --dataset 1000G\|HGDP` appeared as a working example on 7 doc surfaces (+ a `~16 GB` size claim) + the web Settings action, but the raw per-dataset `vcfs/` folders were never uploaded, so the command errors. Doc examples outran what was actually published. | Gate A "every documented `download --what X` resolves on HF" (`list_repo_files`). Only hosted paths may appear as copy-pasteable examples; for anything fetched-from-HF, verify the target exists, and make the CLI error point to the real path (`--what index` / `build-inputs/`). |
| T-22 | A documented command instructed a flag value the shipped validator REJECTS — METHODS told users `--max-total-edits -1` to "disable" the cap, but the validator errors on any negative (`crisprme.py:1834`), so the documented command aborts; a stale in-code comment reinforced it. | For every flag VALUE shown in docs (especially sentinels like `-1`/`0`/`all`), confirm the shipped validator accepts it; keep the doc, the in-code comment, and the validator in agreement. |
| T-23 | METHODS claimed a capability was "not yet" shipped (mega indels "not yet frequency-annotated") though the code + CHANGELOG showed it ships — stale "not yet / forthcoming" claims understate shipped capability in an IND reference. | Grep METHODS + docs for "not yet / forthcoming / will land / in preparation" and verify each against the shipped code/CHANGELOG before tagging (also a fresh-user-audit heuristic). |
| T-24 | A per-release doc carried HARDCODED pinned image tags `pinellolab/crisprme:vX.Y.Z` in copy-paste run commands (the web-interface user guide, 5×) instead of `:latest`, so every release silently left the guide pointing at the PREVIOUS image. | In the version sweep, grep for hardcoded `:vX.Y.Z` tags in runnable doc commands — each must either be on the bump list OR use `:latest`. A bump that edits only `crisprme.py`+`Dockerfile` (CI derives the image tag from the git tag) MISSES ~27 hardcoded `X.Y.Z` doc strings → triage every `grep -rn "<OLD-version>"` hit into (a) current-version claim → bump, (b) pinned run-command tag → bump/`:latest`, (c) historical/"retired in"/"RESOLVED in"/CHANGELOG → keep. |
| T-25 | The only container-wrapped CLI-search instruction for HPC used `apptainer exec` (bypasses the entrypoint → `crisprme.py: not found`) while the same docs warn to use `apptainer run`, and NO correct `apptainer run complete-search`/`generate-report` worked example existed anywhere. | `grep -rn "apptainer exec.*crisprme.py"` must return 0 outside an explicit "wrong-way" example; at least one full worked `apptainer run ... crisprme.py complete-search` + `generate-report` script must exist for the HPC/CLI user. |
| T-26 | The clickable-app / web-form DEFAULT first search selected the 64 GB whole-genome variant index, OOM-killing a typical laptop (exit 137, empty log) with no clickable recovery; and a RAM escape hatch the docs promised ("~16 GB single-chromosome variant search") had no actual WEB-FORM control. | New-user gate: the first-run default must fit ~16 GB (reference-only) with an explicit variant-aware upgrade prompt, OR the app reads container-allocated memory and warns-before-Submit below the measured peak, AND an OOM/exit-137 dialog points to raise-memory / reference-only / single-chromosome. Any RAM escape hatch the docs promise must have a real web-form control, not only a CLI recipe. |

**Gate A sub-checks added from this pass** (wire into section A above):
- Sentinel-value multi-surface sweep (T-22): `grep -rn "disables the cap\|-1.*cap"` must find `-1` only on the ENV VAR `CRISPRME_IUPAC_CAP` (accepts it), never on the CLI flag `--max-total-edits` (rejects negatives). The T-22 fix once landed on METHODS/CHANGELOG/comment but missed README:884.
- Default-index count single-source (T-18): `grep -rnE "929|3,?162|3,?477"` across docs/ must return ONLY raw-source / mega-source-list / retired / 2019-AN-example contexts; any "929 HGDP + 3,162 1kGP" DEFAULT split is the stale pair (correct = 925 + 3,166; 925 = post-QC, 929 = raw source count).
- No retired-index copy-paste even under a Historical banner: a whole-file banner is not cover for a verbatim dead `download --index-name <retired>` command — replace the inline name with the live equivalent.
- Help-text value resolvability: every VALUE in a flag's `--help` must resolve against the DEFAULT repo/config or carry an inline "custom only" caveat (e.g. `--what vcf` help listing `1000G, HGDP`).
- Scripted HF-card + catalog check (make mandatory, not skipped offline): the HF `indexes/` tree must equal the documented 5-index set and the card must carry zero retired names + MGBOA frontmatter.

> **How to use T-rows:** every one is a real incident. When a release pass uncovers a new one,
> add a row here with the symptom and the check, and wire the check into Gate A or B above.
