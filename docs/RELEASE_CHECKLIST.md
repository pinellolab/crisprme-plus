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

> **How to use T-rows:** every one is a real incident. When a release pass uncovers a new one,
> add a row here with the symptom and the check, and wire the check into Gate A or B above.
