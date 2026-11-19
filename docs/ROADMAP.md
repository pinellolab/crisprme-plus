# CRISPRme+ / CRISPRitz roadmap

What is planned after CRISPRme+ 2.10.0 and CRISPRitz 2.10.0, which share their major.minor from
2.10 on ([`RELEASING.md`](RELEASING.md#version-policy)). Nothing here is in 2.10.0. Deferred
ideas with no release named are in [`FUTURE_WORK.md`](FUTURE_WORK.md); the classic-CRISPRme plan
this file held before is kept at its end, for the record.

## CRISPRme+ 2.10.1 (planned)

2.10.1 is to run on CRISPRitz 2.10 and on the indexes of 2.10.0: no index is to be built again
(HPRC's Dictionaries are, for the MNVs' phase). Each item below is a limit of 2.10.0 that METHODS
§10 and [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §14 and §15 state with its measured size, or a
gap the release's gates found.

- **The MNVs' phase (HPRC).** No store the search reads holds an MNV, so on HPRC v2.1 a carrier
  whose haplotype holding an indel row's indel also holds an MNV in the target is listed for a
  sequence it does not have ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §15): 1,163 listings at 94
  sites in the README guide's genome-wide search per sample, 3,334 on SNP rows. 2.10.1 is to store
  each MNV with its haplotype, HPRC's Dictionaries built again and its index as it is, and to read
  it as a second variant of the slot, which takes such a carrier off that row. A target held only
  through an MNV and another indel or MNV stays unreported, as one held through two indels does
  (below); in 171 of the 1,163 listings the haplotype holds another target within the search's
  budget there (METHODS §10), and each such target, held through the indel and the MNV, would then
  go unreported too.
- **The haplotype scan.** In the dense regions of a variant genome that codes every allele (the
  immunoglobulin loci above all), a search above 1 DNA + 1 RNA bulge expands so many IUPAC
  combinations that it runs for hours or out of memory, and 2.10.0 warns of it
  ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §14). 2.10.1 is to search there, per sample, the
  samples' own haplotype sequences, indels included, in place of the IUPAC search (one indel per
  fake-indel window, as the search writes them; a population-level search keeps the IUPAC search
  there): its design was checked on IGL and IGH, up to 4 mismatches + 2 + 2 bulges, where the lines
  of the IUPAC search that a real haplotype can give are exactly the scan's. The engine's hooks for
  it (a region skip, keep-only, the skipped-leaf dump and a capability report) ship in CRISPRitz
  2.10.0, which 2.10.1 runs on, and the indexes of 2.10.0 serve it as they are: no index is rebuilt.
- **One density pass before the search**, defined by the genome and shared by the scan, the speed
  steps of 2.10.0 and the mega's density guard (in 2.10.0 the dense-region file is written from
  each search's own raw targets).
- **A faster population-level post-analysis of the dense loci.** In the release gate, 2.10.0's
  genome-wide population-level search of the README guide with NRG on HGDP + 1000 Genomes took 2.20
  times 2.8.1's time, its search phases about 2.8.1's: a population-level post-analysis reads no
  haplotype, so the haplotype pre-check skips none of the dense loci's windows, the immunoglobulin
  loci among them (METHODS §10). 2.10.1 is to make that post-analysis faster, and to time the
  search again (C9, with C8 and C10, the per-sample ones beside it), with the post-analysis changes
  made after the gate.
- **Known limits of 2.10.0 the release fixes** ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §15):
  one target written as two rows, or twice, in a microsatellite (the merge); and a SNP row's carrier
  left off beside an indel that adds or deletes more than the 500 bases the store is read around a
  target, less the target's length (one carrier on two rows of the release gate's HPRC search), the
  store to be read as far as a haplotype's own indels need.
- **One rule for a complex record inside an indel's deleted bases.** CRISPRme takes a carrier off
  an indel row where its slot holds a second record, inside the first's deleted or replaced bases,
  that puts bases of its own there (METHODS §4), and the real-data check reads that slot without it
  (31 pairs of records on one slot in HPRC v2.1, 249 slots): 2.10.1 is to decide which reading is
  right and make the other say the same.
- **The real-data check** (`scripts/check_indel_rows.py`): to check `Variant_info_spacer+PAM` (each
  variant's number in its target) and redundant alternative alignments, which it does not today
  ([`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md), step 4c); to hold less memory (one process holds
  the whole alternative-alignment table: 106.5 to 116.1 GiB on the release gate's largest searches
  genome-wide, 124.6 GiB with a SNP-row report); to check a search beside the mega on chrY, which
  only the mega's VCFs hold (it refuses the whole genome there today, so the gate checked each
  other chromosome); to check SNP rows for completeness, so that it sees a SNP row the rule takes
  off whole; to name every record a missed target needs in its completeness line, a record past a
  long deletion included; and to tell a run's genome where it refuses one today.
- **CRISPR-Bulge scores of near-ties.** A CRISPR-Bulge score can move by about 1e-5 with how a
  chromosome's rows are batched, or between the amd64 and arm64 images
  ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §5), which flips the order of near-ties in the
  highest-CRISPR-Bulge column: the scoring is to be made batch-independent, or the gates to compare
  with a tolerance.
- **The web app**: temporary folders by name in a shared `CRISPRME_TMPDIR` that is not sticky; the
  History page's reads of a results folder's files through the safe-file reader, and a `Results/`
  folder that is a link (it answers HTTP 500 today); the reads of headers and PAM files; and its
  stop under Apptainer, where a SIGINT sent to the `crisprme.py` process alone (as a script or a
  service manager sends it) left the worker process of the Settings page's jobs running and
  holding the port (a Ctrl-C in its terminal, which reaches every process, stops them all; Docker
  is unaffected, as a container's end stops every process in it).
- **The rebuild's tooling**: an INDELS index of an exact budget refused early, not late; the
  completeness of `$REF` checked; the checklist's lessons of the 2.10.0 rebuild
  ([`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md), D and E).

## Open, no release named

- Two indels on one haplotype: a fake-indel window holds one indel, so a target that a haplotype
  holds only through two indels is not reported for its carriers, and a carrier whose second indel
  pushes a SNP the row names out of its own target is listed on that row
  ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §15; METHODS §10 gives the release gate's counts).
  The haplotype scan does not lift it: it writes, as the search does, one indel per window.
- An unphased call of two ALTs at one position (`1/2`) is listed on a row that needs the
  reference base there, and on a classic install a carrier with an unphased call in a VCF that is
  otherwise phased only on the row of all its variants ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md)
  §15). No index of 2.10.0 meets either: none of its panels holds a `1/2` call (measured on their
  VCFs: METHODS §10), and every published variant index of a genotyped panel is built with
  `--samplesID`; an index one builds oneself can.
- The advice a search prints for an INDELS index it cannot use, on a sites-only dataset whose
  fake-indel genome is not in place, removes the enriched genome before a `build-index-only`
  that refuses the dataset ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §12 says what to do
  instead); a sites-only branch would say to put the kit back.
- Messages: `pool_index_indels.py`'s refusal of a fake-indel genome and a log of two enrichments
  could name the paths to remove, as `build-index-only` does; the sites-only advice of
  `indel_window` could be shorter, pointing to TROUBLESHOOTING §12; a few printed `rm -rf`
  commands quote paths by hand, so a path holding an apostrophe cannot be pasted.

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
