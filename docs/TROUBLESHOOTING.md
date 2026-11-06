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

A home CRISPRme cannot write in (Apptainer `--no-home`, or `--home` on a read-only folder)
needs nothing more: matplotlib then keeps its configuration and font
cache in a temporary folder there (`MPLCONFIGDIR`), removed when the command exits (a command
killed by a signal, such as `scancel` or a timeout, leaves it in the temporary folder). When
`~/.cache` cannot be written, that folder is `XDG_CACHE_HOME` too, so the caches kept under it
go there as well (fontconfig's, written when matplotlib builds its font cache, and those of any
other program a command runs), but not Hugging Face's folder: `HF_HOME` stays
`~/.cache/huggingface`, so a login saved there with `hf auth login` is still found
(`publish-index` takes it when given no token, and a download on a conda or source install
sends it). A `MPLCONFIGDIR`, `XDG_CACHE_HOME` or `HF_HOME` you set is used as it is.

`docker run --read-only` makes `/tmp` read-only as well, so it needs a writable temporary folder:
add `--tmpfs /tmp`, or set `CRISPRME_TMPDIR` to a folder of the mounted data folder (or to that
folder itself: `-w /DATA -e CRISPRME_TMPDIR=/DATA`). Without one,
the CRISPR-Bulge scorer and the report have nowhere to write, and complete-search stops before
it starts, saying so (until 2.9.0 the search ended with CRISPR-Bulge scores of -1, no report zip,
and only a WARNING in its error log).

---

## 2. A genome-wide search needs a lot of disk — plan for it

The default variant-aware install (`download --what all`, then
`download --what index --index-name NRR_3_hg38+hg38_HGDP1kGP`) takes **~150 GB** on disk
(measured: 150.0 GB) before any results:

| Item | Approx size |
|---|---|
| Reference genome (hg38) | ~3 GB |
| Reference + variant-aware indexes (`NRR_3_hg38`, `NRR_3_hg38+hg38_HGDP1kGP` and its indel index) | ~94 GB |
| Its `Dictionaries/`: genotype store ~30 GB, allele-frequency registry ~8.5 GB, indel genotype store ~7.7 GB, indel logs ~6.5 GB (fetched by default; per-sample analysis, this index's default, reads the store) | ~52 GB |
| Annotations | ~0.1 GB |
| Per-run results (grows with mismatches/bulges; larger for per-sample runs, the default on a phased index) | GB–tens of GB |

The install peaks around 161 GB while the archives unpack. Keep everything on a
disk with **~200 GB free** before a first variant-aware run.
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
  sites-only index such as `NRR_3_hg38+hg38_mega`.
- Named carriers / personal risk cards come from the **genotype store** bundled with
  the index (`download --what index`, unless you pass `--no-genotypes`). Without it
  (and without the classic per-sample dictionaries) a search runs population-level,
  naming the carriers of an indel alone (the index holds them), and `--per-sample`
  falls back to it with a `WARNING` that names the dataset: re-download without
  `--no-genotypes` to get carriers.
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
- If the `cbulge` scorer environment or its weights are missing or broken,
  `complete-search` stops before it searches and names the fix
  (`crisprme.py scorer-env create`; `crisprme.py scorer-env doctor` says what is wrong).
  Set `CRISPRME_SKIP_SCORER_ENV=1` to run it **CFD-only** anyway. `targets-scoring` and
  searches started from the web interface do not check first: with the env absent they
  complete **CFD-only**. A CFD-only run has the CRISPR-Bulge column at `-1` and says so
  (a search in its `log_verbose.txt`, `targets-scoring` on stdout).

GPU helps on large/dense searches; a tiny single-guide run is often faster on CPU
due to initialization overhead.

---

## 6. Download / data setup

- **Get everything a reference scan needs:** `crisprme.py download --what all` fetches
  the genome, annotations, PAMs, sample lists, and the default reference index.
- **Add a variant-aware index:** `crisprme.py download --what index --index-name
  NRR_3_hg38+hg38_HGDP1kGP` (the recommended phased HGDP + 1000 Genomes default, of the
  default SpCas9 NRR PAM).
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

---

## 10. Hidden build folders left after a container stop

**Symptom.** Storage counts more under *Indexes* than the installed indexes hold, and
the data folder has hidden `genome_library/.build-index-only-*` (or
`Genomes/.build-index-only-enrich-*`, `Dictionaries/.build-index-only-dicts-*`)
folders; the Settings page names them when a build ends.

**Cause.** `build-index-only` builds in such a folder and removes it on every way
out but a kill. A `docker stop` while a Settings build runs kills it. A later build
removes the folder of a run that is gone, but one made in another container (each
`docker run --rm` has a new host name) cannot be checked, so it is left in place,
with a line on stdout that says so.

**Fix.** Once no CRISPRme container is building, remove them:
`rm -rf genome_library/.build-index-only-* Genomes/.build-index-only-enrich-* Dictionaries/.build-index-only-dicts-*`
from your data folder.

## 11. A search stops and says no index serves it

Each of these names the index and what to do; none of them removes anything.

- **`no variant index here serves this search (...)`.** The dataset's indexes are for another PAM or
  fewer bulges, and neither its VCFs nor its enriched genome are here to build one. A download holds
  its own index only: `NRR_3_hg38+hg38_HGDP1kGP` searches NRR, and the PAMs NRR holds (NGG, NAG,
  NGA, NRG), pruned to them (and, with a 20-nt guide, up to 2 bulges of each type). Search with what
  the message lists, or install another index ([Precomputed indexes](PRECOMPUTED_INDEXES.md)).
- **`index '...' is not installed, and this run cannot build it`.** The library is on a read-only
  mount, has no write permission or is full, and lacks the index this PAM and these bulges need. Make
  it writable, or run `build-index-only` where it is, or search with no bulges, which opens no index.
- **`... was built for a 20-nt guide, and this search has a 23-nt guide ...: it stores too few bases
  per site`.** This search builds its own index under a name that says its guide length (since
  2.9.0: `23-NGG_2_hg38`), and an index of that name is there, whose trees were built for another
  guide (a folder renamed by hand). Name it for the guide it was built for, as the ERROR says (`mv
  .../23-NGG_2_hg38 .../20-NGG_2_hg38`), or move it aside. An index of the older name
  (`NGG_2_hg38`) built for another guide is left as it is, and the search builds its own beside it.
  Any published index serves a 21-nt guide at up to 1 bulge of each type.
