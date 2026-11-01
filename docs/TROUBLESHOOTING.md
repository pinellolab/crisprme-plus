# CRISPRme Troubleshooting

Practical fixes for problems you may hit while following the
[Docker Quickstart](DOCKER_QUICKSTART.md), the README, or the
[Web Interface guide](crisprme_web_interface_user_guide.md). Most are environment
issues (disk, `/tmp`, memory, GPU) rather than CRISPRme bugs.

If something here is wrong or missing, please open an issue — the goal is that a
first-time user following the docs never gets stuck.

---

## 1. "No space left" or temp-file errors during a run (small or full `/tmp`)

**Symptom.** A search fails partway with `No space left on device`, a
`sort-bed`/`bgzip`/`bcftools` error, or a `FileNotFoundError` on a file under `/tmp`.
Common on HPC nodes and inside containers, where `/tmp` is often a small
RAM-backed `tmpfs` (a few GB) while a genome-wide run needs tens of GB of scratch.

**Fix — point CRISPRme's temp at a roomy, writable directory with `CRISPRME_TMPDIR`.**
This single knob redirects *all* temporary files (CRISPRme's own, plus the
`sort-bed` / `bcftools` / `bgzip` children it spawns) — it is exported as `TMPDIR`
for the whole run. Point it at a path on the same big disk as your data/results:

```bash
# conda / source install
export CRISPRME_TMPDIR=/path/to/big/disk/crisprme_tmp
crisprme.py complete-search ...
```

```bash
# Docker: put the temp dir on the mounted data volume and tell the container to use it
docker run --rm -v "${PWD}:/DATA" -w /DATA \
  -e CRISPRME_TMPDIR=/DATA/crisprme_tmp \
  pinellolab/crisprme:latest crisprme.py complete-search ...
```

The directory is created if it does not exist. If it is not writable, CRISPRme
warns and falls back to the default (`$TMPDIR`, else `/tmp`). You can also just set
the standard `TMPDIR` — CRISPRME_TMPDIR takes precedence when both are set.

---

## 2. A genome-wide search needs a lot of disk — plan for it

A full install + one variant-aware genome-wide run is **tens of GB**:

| Item | Approx size |
|---|---|
| Reference genome (hg38) | ~3 GB |
| Variant-aware index (e.g. `NRG_3_hg38+hg38_HGDP1kGP`) | ~40 GB |
| Its genotype store + dictionaries (`--per-sample`) | ~50 GB |
| Per-run results (grows with mismatches/bulges; `--per-sample` is larger) | GB–tens of GB |

Keep everything on a disk with **~150 GB free** before a first variant-aware run.
Use `CRISPRME_TMPDIR` (Section 1) on the same roomy disk.

**Reclaiming space:** by default `complete-search` removes bulky intermediates when
it finishes (keeping the finalized results + report). Pass `--keep-intermediates`
to retain them (for debugging / custom parsing).

---

## 3. (Historical) PAM `NRG` (or other partially-degenerate PAM) + bulges warning

Some releases before 2.7.10 printed this at the start of a search:

> `WARNING: PAM motif 'NRG' is a partially-degenerate, odd-length motif used with
> bulges. This combination can crash the underlying search engine (issue #105 ...)`

CRISPRme+ no longer prints it. The crash it warned about was fixed in CRISPRitz 2.7.1,
and every supported install pins a later CRISPRitz (2.9.1 for CRISPRme+ 2.8.1). Only if
you run CRISPRme against an older, unpinned CRISPRitz and a search **aborts with a memory
error** while using a degenerate PAM (`NRG`, `NNN`, …) *with* bulges, either:

- retry with bulges set to 0 (`--bDNA 0 --bRNA 0`), or
- use a fully-specified PAM (e.g. `NGG` instead of `NRG`).

---

## 4. Per-sample runs need memory (and a genotyped index)

Per-sample analysis resolves genotypes across the whole panel (e.g. 4,091 individuals
for the HGDP+1000G index) and is both slower and more memory-hungry than the
population-level analysis. On a large panel budget **tens of GB of RAM** per run. It is
the default on a genotyped + phased index (HGDP + 1000 Genomes, HPRC), and `--per-sample`
forces it on any other genotyped one.

- `--per-sample` requires a **genotyped** index. It no-ops (with a warning) on a
  sites-only index such as `NRG_3_hg38+hg38_mega`.
- Named carriers / personal risk cards come from the **genotype store** bundled with
  the index (`download --what index`, unless you pass `--no-genotypes`). If carriers
  show as degraded/absent, the store was skipped — re-download without `--no-genotypes`.
- If you only need site-level results, pass `--population-level` — it is much lighter.
  On a phased index leaving out `--per-sample` is not enough, since per-sample is that
  index's default; population-level is the default only on a sites-only, unphased,
  hybrid or reference-only index.

---

## 5. GPU vs CPU for the CRISPR-Bulge scorer

The ML scorer (CRISPR-Bulge) runs on CPU by default and **always falls back to CPU**
when no GPU is available — a CPU-only machine is fully supported. To use a GPU:

- Docker: add `--gpus all` on an NVIDIA host and pass `--compute-backend cuda`.
- The accepted devices are `cpu | cuda | metal | auto` (`complete-search
  --compute-backend`; `targets-scoring --device` accepts the same set).
- If the `cbulge` scorer environment is absent, the run completes **CFD-only** (the
  CRISPR-Bulge column reads `-1`/N/A) rather than failing.

GPU helps on large/dense searches; a tiny single-guide run is often faster on CPU
due to initialization overhead.

---

## 6. Download / data setup

- **Get everything a reference scan needs:** `crisprme.py download --what all` fetches
  the genome, annotations, PAMs, sample lists, and the default reference index.
- **Add a variant-aware index:** `crisprme.py download --what index --index-name
  NRG_3_hg38+hg38_HGDP1kGP` (the recommended phased HGDP + 1000 Genomes default).
- The public repository ships **prebuilt indexes** that already encode the variants, so
  you do not need the VCFs to search the shipped datasets. It also hosts their source
  VCFs (`hg38_HGDP1kGP`, `hg38_HPRC`, `hg38_mega`), for building an index yourself:
  `crisprme.py download --what vcf --dataset hg38_HGDP1kGP`.
- A partial/interrupted download is re-fetched and structurally verified on the next
  run; set `CRISPRME_VERIFY_DEEP=1` to force a deep re-check, or `SKIP_VERIFY=1` to
  skip verification.

---

## 7. Running several searches at once on one install

Independent `complete-search` jobs can run in parallel, but give each job its own
**output directory** (`--output`), and prefer a separate `CRISPRME_TMPDIR` per job if
`/tmp` is tight. Jobs correctly share the read-only reference genome + index. (A
long-standing collision where two jobs stripping the same COSMIC-containing annotation
raced on a temp file is fixed in 2.7.10.)

---

## 8. Very high `--thread` can stall post-analysis

The post-analysis pool is bounded by `--thread`. On many-core nodes, extremely high
values (e.g. `--thread 200`) can stall the per-chromosome pool. A value in the
**16–32** range is a good balance on large machines; raising it further rarely helps
and can hurt.

---

## 9. Assembly-search (personal genome) coordinate lifting

`assembly-search` lifts off-targets between a personal assembly and hg38 using the
MIT `pyliftover` library (no external `liftOver` binary required). Pass explicit chain
files with `--chain-paternal` / `--chain-maternal`; sites that do not map, or map
ambiguously to more than one chain, are reported as unmapped.
