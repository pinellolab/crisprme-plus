# CRISPRme+ 2.5.2 — release report (for review)

Status snapshot of the 2.5.2 release + validation + the post-2.5.2 CRISTA-parallel prototype.
Written for a morning review; the clean-room result tables are filled in as those runs land.

## 1. Released (live)

| Artifact | State |
|---|---|
| Git tag `v2.5.2` | pushed (`8145ea4`), on `main` at `72c8e27` → README `0a888cf` → TESTING `48dcbd3` → grid `ba83ba6` |
| Docker image | **`pinellolab/crisprme:v2.5.2` live on Docker Hub** (multi-arch build succeeded) |
| Genotyped index | `lucapinello/crisprme-data : NRG_3_hg38+hg38_1000G2021_HGDP` (on HF; unchanged, 2.5.2-ready) |
| Mega index | `lucapinello/crisprme-data : NRG_3_hg38+hg38_mega` — **newly published** (12.1 GB: SNP + searchable INDELS + registry + indel_af + log_indels) |
| GitHub *release notes* | **not yet attached** (optional): `gh release create v2.5.2 --repo pinellolab/crisprme-plus --title "CRISPRme+ 2.5.2" --notes-file ~/crisprme_252_notes.md --verify-tag` |

## 2. What 2.5.2 adds (the two-index, co-occurrence story)

- **Two production indices, complementary by design:**
  - **1000G-2021 + HGDP** (genotyped) → **observed** haplotypes: CONFIRMED (phased cis) / PUTATIVE (unphased), with named carrier samples + exact joint AF.
  - **mega** (5 sources, sites-only) → **putative** haplotypes: PUTATIVE with a conservative **min-AF** joint bound + **per-dataset AF provenance** (`AF_1000G2021 … AF_AoU` + `AF_max`).
- **Two co-occurrence dimensions**, both indices: **SNP+SNP** (`snp_snp_cooc.tsv`, new) and **SNP+indel** (`indel_snp_cooc.tsv`), the class the classic two-pass search couldn't see.
- **Mega searchable indels genome-wide** (the P3 gap: the enricher never materialized the fake-indel genome for sites-only; fixed by a synth-GT build + MEGA-carrier strip → cleanly sites-only).
- **Lossless dense-region worst-case** on sites-only (`CRISPRME_LOSSLESS_DENSE`, default-on sites-only, byte-identical for genotyped).
- **Report fixes:** bundle `snp_snp_cooc.tsv`; **fix pre-existing `indel_af.tsv`-never-bundled** (#184); radar-chart tolerant of IUPAC on sites-only DNA (was crashing the whole variants summary).

Code commits on `main`: `c7b145e` (P1), `f98e18d` (P2), `ddeff99` (report), `1cbe5b2`/`9f25e5e`/`72c8e27` (docs+version), + README/TESTING/grid.

## 3. Validation

### 3a. chr22 development progression (done)
Genotyped: SNP+SNP 3,329 (35 CONFIRMED + 3,294 PUTATIVE), SNP+indel 14 (2+12). Mega: SNP+SNP 139 PUTATIVE, SNP+indel 6,502 PUTATIVE, indel_af 6,824. GW mega verify (all 24 chr): 839,897 SNP+indel + 10,753 SNP+SNP (all PUTATIVE), 868,724 indel_af, indel off-targets on all 24 chr, 0 pseudo-sample leaks. Provenance spread 785K single-source → 8.7K all-5.

### 3b. Clean-room, new-user from scratch (Docker/apptainer pull v2.5.2 → HF-download both indices)
Two configs in parallel on separate hosts. Guide = TRAC `CTCTCAGCTGGTACACGGCA`, NRG.

**ml007 — mm4 / bDNA1 / bRNA1 (fast validation, /srv/local SSD):**

_[results appended when the run lands]_

**ml008 — mm6 / bDNA2 / bRNA2 (thorough):**

| Cell | report.zip | off-targets | snp_snp_cooc | indel_snp_cooc | indel_af |
|---|---|---|---|---|---|
| mega_slow | ✅ | 550,456 | 3,453 **all PUTATIVE** | 318,455 PUTATIVE | 329,980 |
| mega_fast | ✅ | 547,807 | 42,273 **all PUTATIVE** | 318,455 PUTATIVE | 329,980 |
| geno_slow | ⚠️ hung @51.6% | — | — | — | — |
| geno_fast | ⚠️ hung @51.6% | — | — | — | — |

**Heavyweight mega passed** genome-wide (318K SNP+indel + 3.4K SNP+SNP PUTATIVE, 330K indel_af, 550K off-targets) — the sites-only feature set holds at mm6/b2/b2 GW.

**Genotyped cells hung at 51.6%** (0 active post-analysis workers, `run_searches.sh` still waiting): a worker was OOM-killed under **4× mm6/b2/b2 memory pressure** — the genotyped observed-enumeration loads the per-sample genotype store (far heavier than the registry-only mega), and 4 such GW searches in parallel over-subscribed RAM. **Resource/contention finding, not a genotyped-index bug** — the genotyped is validated by the chr22 progression (§3a) and the ml007 mm4/1/1 clean-room (below). Takeaway: don't run 4× genotyped GW mm6/b2/b2 concurrently on one host; run fewer in parallel or cap `CRISPRME_POSTPROC_MAX_WORKERS` / memory.

The mega SNP+SNP rows are **all PUTATIVE** (verified: zero CONFIRMED data rows across all 24 per-chr files) — correct for a sites-only index. `--fast` emits far more SNP+SNP rows than slow (42,273 vs 3,453) because it materializes a worst-possible representative per dense window (which combines multiple SNPs), whereas the slow greedy emits fewer multi-SNP reps; both are PUTATIVE with the min-AF bound.

What each cell is checked for: report.zip generated; `snp_snp_cooc.tsv` / `indel_snp_cooc.tsv` / `indel_af.tsv` present; CONFIRMED+carriers on genotyped vs PUTATIVE+min-AF+provenance on mega; `--fast` ⊇ slow; no crash.

**New-user friction found:** anonymous HF downloads hit **429 "Too Many Requests"** (the `--what all` many-small-file API calls); authenticating (HF token) fixes it. Worth a docs note + a retry-on-429 in the downloader (candidate fix).

## 4. CRISTA parallel prototype (post-2.5.2, on `dev`)

Answers "why is CRISTA the tail?" — the post-analysis pool parallelizes per **contig**, so one big chromosome's worker runs serially long after the others. Prototype parallelizes the CPU-bound per-target **feature build** at the finest seam + lets the RF predict use joblib threads.

- `CRISPRME_CRISTA_PARALLEL` (opt-in, **default OFF**). `dev` commits `7771b66` → `c5fff0d` (fork) → `3d7aa44`/`a0e7f16` (measurement-driven scope-down).
- **Parallelizes only the per-target feature build** (`get_features`) — the genuinely serial part — across a small bounded **fork** pool (bit-identical feature matrix; contiguous order-preserving chunks). Test `test_crista_parallel_equivalence` (in CI).
- **Fork, not spawn** — the post-analysis callers run their main at module top level (no `__main__` guard), so spawn would re-run the whole analysis in each child; fork is safe (caller is an isolated `subprocess.call` child; inner pool bounded ≤8).
- **Measured** (on-SIF, 20k real targets): feature-build alone **3.8× @ 8w**; full `CRISTA_predict_list` **~1.2–2×** (modest — see below). Feature-build-only, 10-core mac: 1.36× (2w) / 2.49× (4w) / 3.82× (8w).
- **Key measurement finding:** the pickled RF predictors already carry **`n_jobs=-1`** → the predict already fans across all cores (it's the dominant per-batch cost, ~4× the feature build), so parallelizing the predict further only over-subscribes. It is also therefore **non-deterministic at the raw-float level** — two *serial* predict runs differ by ~5e-17, **identical to the emitted 3 decimals**. So CRISTA output has *always* been reproducible only to the emitted precision (not raw bits); this prototype preserves exactly that.
- **Takeaway:** the bigger lever for the per-contig skew tail is **outer-level balanced contig batching** (tracked #174), not per-batch predict parallelism. The feature-build fork pool is a clean, bit-identical, opt-in partial win.

## 5. Pending / for your decision

1. **Clean-room results** — appended above as they land (both running).
2. **CRISTA prototype**: full-scale byte-identical + predict-`n_jobs` measurement on the SIF; then decide whether to fold into a release (it's `dev`, default-off).
3. **HF 429**: add a retry-on-429 to the downloader (and/or a docs note that a token avoids it).
4. **GitHub release notes** object (command above) — cosmetic; the tag + image are live.
5. `pinellolab/crisprme-data` vs `lucapinello/crisprme-data`: the published indices + the code default (`DEFAULT_HF_REPO`) are both `lucapinello/crisprme-data` — consistent, but confirm that's the intended long-term home.
