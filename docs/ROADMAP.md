# CRISPRme+ / CRISPRitz roadmap

What is planned after CRISPRme+ 2.10.0 and CRISPRitz 2.10.0, which share their major.minor from
2.10 on ([`RELEASING.md`](RELEASING.md#version-policy)). Nothing here is in 2.10.0. Deferred
ideas with no release named are in [`FUTURE_WORK.md`](FUTURE_WORK.md); the classic-CRISPRme plan
this file held before is kept at its end, for the record.

## CRISPRme+ 2.10.1 (planned)

- **The haplotype scan.** In the dense regions of a variant genome that codes every allele (the
  immunoglobulin loci above all), a search above 1 DNA + 1 RNA bulge expands so many IUPAC
  combinations that it runs for hours or out of memory, and 2.10.0 warns of it
  ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §14). 2.10.1 is to search there the samples' own
  haplotype sequences, indels included, in place of the IUPAC search: its design was checked on
  IGL and IGH, up to 4 mismatches + 2 + 2 bulges, where the lines of the IUPAC search that a real
  haplotype can give are exactly the scan's. The engine's hooks for it (a region skip, keep-only,
  the skipped-leaf dump and a capability report) ship in CRISPRitz 2.10.0, which 2.10.1 runs on,
  and the indexes of 2.10.0 serve it as they are: no index is rebuilt.
- **One density pass before the search**, defined by the genome and shared by the scan, the speed
  steps of 2.10.0 and the mega's density guard (in 2.10.0 the dense-region file is written from
  each search's own raw targets).
- **Known limits of 2.10.0 the scan or the release fixes**
  ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §15): a carrier listed on an indel row for a SNP
  that a second indel of its haplotype pushes out of its own target (the scan reads each carrier's
  own bases); one target written as two rows, or twice, in a microsatellite (the merge).
- **The real-data check** (`scripts/check_indel_rows.py`) checks `Variant_info_spacer+PAM` (each
  variant's number in its target) and redundant alternative alignments, which it does not today
  ([`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md), step 4c).

## Open, no release named

- An unphased call of two ALTs at one position (`1/2`) is listed on a row that needs the
  reference base there, and on a classic install an unphased carrier only on the row of all its
  variants ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §15).
- The advice a search prints for an INDELS index it cannot use, on a sites-only dataset whose
  fake-indel genome is not in place, removes the enriched genome before a `build-index-only`
  that refuses the dataset ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §12 says what to do
  instead); a sites-only branch would say to put the kit back.
- Messages: `pool_index_indels.py`'s refusal of a fake-indel genome and a log of two enrichments
  could name the paths to remove, as `build-index-only` does; the sites-only advice of
  `indel_window` could be shorter, pointing to TROUBLESHOOTING §12; a few printed `rm -rf`
  commands quote paths by hand, so a path holding an apostrophe cannot be pasted.
- The real-data check reads a whole enriched chromosome for each SNV written in lower case: slow
  on a user's VCF written so (none of the published panels holds one).

---

## Historical: the classic-CRISPRme plan (2.1.12 → 2.1.13)

Kept for the record. Its items shipped (Python 3.11, Dash 2, #107, #112, #113) or were superseded
(#108).

This document records the release strategy: a **stable 2.1.12** on the current
Python (so the Dash web app keeps working) with every known bug fixed, followed
by a **2.1.13** modernization + new-features release (Python 3.11, Dash upgrade,
dependency bumps, deferred features).

---

### 2.1.12 — stability release (Python 3.8)

**Goal:** a stable, fully-tested release with all known bugs fixed. Deliberately
**stays on Python 3.8** so the Dash web interface and its dependency stack keep
working unchanged — **no modernization risk**.

**Included (merged to `main`):**
- VCF / gnomAD / concurrent-run fixes (#96); annotation `.gz` handling (#109);
  malformed-cluster-row guard (#114).
- Cas9 + Cas12a validation benchmarks + brute-force ground truth + Rust generator,
  and the `validate-benchmarks` CI (#116, #120).
- Scalability guard `search_budget.py` + analysis (#118).
- Native Apple-Silicon (arm64) multi-arch Docker image + low-memory startup
  warning (#121).
- Input hardening (#125): non-ATCG PAM guard (**#94**), `--gene_annotation`
  sorting, multi-line-PAM + PAM-filename validation, a `UnicodeDecode` guard
  (**#52**), and a **#105** guard (clear warning on degenerate-PAM + bulge).
- Docs: `INPUT_FORMATS.md`, `ISSUE_AUDIT.md`, `RELEASING.md`, release skill.
- Full audit of all 8 open + 58 closed issues — nothing closed by mistake.

**CRISPRitz side (2.7.1):** the #105 heap-overflow fix (PR #19), a CI safety net
(PR #20), and the stale-`.pyc` cleanup (#17). CRISPRme 2.1.12 runs on the existing
**crispritz 2.7.0** with the CRISPRme-side #105 guard, so it is not blocked on the
crispritz/Bioconda chain. Repinning to crispritz 2.7.1 can fold in at tag time if
Bioconda is ready, else it moves to 2.1.13.

**Explicitly deferred out of 2.1.12 (→ 2.1.13):** the input validator (#112, has a
`samplesID` false-positive that aborts valid 1000G runs), Python 3.11, the Dash
upgrade, and the assembly-search feature (#113).

---

### 2.1.13 — modernization + new features

#### A. Python 3.11 migration
Both tools are currently pinned to Python 3.8 (EOL). Assessed difficulty: **MEDIUM**.
- **CRISPRitz → 2.8.0:** drop stale `.pyc` (#17, done in PR #21); bump the pins;
  **adopt CRISPR-HAWK's `scores/azimuth`** (a modern, sklearn-1.1-loadable model) to
  restore on-target Doench scoring — the only hard part, since CRISPRitz's own
  azimuth model is a sklearn-0.21 pickle; **fix the broken standalone Dockerfile**
  (fails today at `conda install python=3.8`). Groundwork: PRs #21, `azimuth-crisprhawk`.
  The off-target (CFD) pipeline already runs on 3.11 (tested).
- **CRISPRme → 2.2.0:** pipeline fixes — `DataFrame.append`→`pd.concat`
  (`process_summaries.py`), `seaborn-poster`→`seaborn-v0_8-poster` (4 sites),
  bare `python`→`sys.executable` (`generate_sample_card.py`, `crisprme.py`); pin
  `numpy<2`; repin to the crispritz-py3.11 build. The ML off-target scorer is now
  **CRISPR-Bulge** (TensorFlow, no sklearn pickle), so no sklearn forward-compat shim
  is needed on the CRISPRme side. Groundwork: PR `py311-core`. **Prerequisite:** a
  crispritz **Python-3.11 Bioconda build** (compiled; linux-64 + linux-aarch64) must
  exist first.

#### B. Dash web-app modernization (the long pole)
The web app is written for **Dash 1.x** (`dash 1.10`, `flask 1.1`, `werkzeug 1.0`,
`itsdangerous<2.0`, and the retired `dash-core-components`/`dash-html-components`).
Migrate to **Dash 2.x**: rewrite 11 files' imports to `from dash import dcc, html,
Input, Output, State`, `app.run_server`→`app.run`, and bump Flask/Werkzeug/
itsdangerous/plotly. Not covered by the benchmark CI — **validate with Playwright**
click-through (as in the original audit).

#### C. Deferred features / bugs
- **#112 input validator** — re-enable after fixing the `samplesID` sample-check
  false-positive (it currently `ERROR`s on the valid 1000G workflow). It's the
  FDA-requested pre-flight check.
- **#113 assembly-search** subcommand (re-test on current `main` + review).
- **#107** true `--max-edits` (total mm+bulges) filter.
- **#108** deep scalability fix (stream `best*` via `sort -m`; flush at cluster
  boundaries) — see `docs/SCALABILITY_ANALYSIS.md`.
- Separate **DNA/RNA bulge count** columns in the output.
- CRISPRitz **#13** (non-human/scaffold genome crash — check if the #105 fix
  resolves it) and the multi-motif / additional-PAM improvements.

#### D. Pickled-model modernization (cross-cutting)
This now applies only to **azimuth (CRISPRitz)**, which ships a model pickled under
sklearn ~0.21 (paths removed in 0.22). Fix via a `sys.modules` alias shim and/or
re-export, and **validate scores numerically** against the benchmark ground truth
(a clean unpickle is not sufficient). On the CRISPRme side the ML off-target scorer
is now **CRISPR-Bulge** (TensorFlow, no sklearn pickle), so no such shim is needed.

---

### Order of operations (priority: stability first)

1. **Ship CRISPRme 2.1.12** — finalize docs `#124` + changelog `#123`, tag `v2.1.12`
   (triggers the arm64 Docker publish + Bioconda autobump). Old Python, working Dash,
   all bugs fixed.
2. **CRISPRitz 2.7.1** (parallel) — Manuel reviews/merges #19 (#105) + #20 (CI) +
   the #17 cleanup → tag `v2.7.1` → Bioconda → optionally repin CRISPRme.
3. **2.1.13 modernization** — Python 3.11 (CRISPRitz 2.8.0 → CRISPRme 2.2.0) + Dash
   2.x + deferred features, each gated by the benchmark CI (and Playwright for the web).

### Validation gates
- `validate-benchmarks` (byte-identical Cas9 + Cas12a off-targets, incl. CFD **and**
  CRISPR-Bulge) gates every CRISPRme change.
- CRISPRitz CI (build + AddressSanitizer #105-regression + native arm64).
- Playwright click-through for the Dash web app.
- Multi-arch Docker build (linux/amd64 + linux/arm64) for both images.
