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
it starts, saying so (until 2.10.0 the search ended with CRISPR-Bulge scores of -1, no report zip,
and only a WARNING in its error log).

---

## 2. A genome-wide search needs a lot of disk — plan for it

The default variant-aware install (`download --what all`, then
`download --what index --index-name 20-NRR_3_hg38+hg38_HGDP1kGP`) takes **~150 GB** on disk
(measured: 150.0 GB, with 2.8.1's `NRR_3_hg38+hg38_HGDP1kGP`) before any results:

| Item | Approx size |
|---|---|
| Reference genome (hg38) | ~3 GB |
| Reference + variant-aware indexes (`20-NRR_3_hg38`, `20-NRR_3_hg38+hg38_HGDP1kGP` and its indel index) | ~94 GB |
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
and every supported install pins a later CRISPRitz, of CRISPRme's own major.minor (2.10.0 for
CRISPRme+ 2.10.0). Only if you run CRISPRme against an older, unpinned CRISPRitz and a search
**aborts with a memory error** while using a degenerate PAM (`NRG`, `NNN`, …) *with* bulges, either:

- retry with bulges set to 0 (`--bDNA 0 --bRNA 0`), or
- use a fully-specified PAM (e.g. `NGG` instead of `NRG`).

---

## 4. Per-sample runs need memory (and a genotyped index)

Per-sample analysis resolves genotypes across the whole panel (e.g. 4,091 individuals
for the HGDP+1000G index) and reads the panel's genotype store: on a large panel budget
**tens of GB of RAM** per run. It is the default on a genotyped + phased index (HGDP +
1000 Genomes, HPRC), and `--per-sample` forces it on any other genotyped one.

- `--per-sample` requires a **genotyped** index. It no-ops (with a warning) on a
  sites-only index such as `20-NRR_3_hg38+hg38_mega`.
- Named carriers / personal risk cards come from the **genotype store** bundled with
  the index (`download --what index`, unless you pass `--no-genotypes`). Without it
  (and without the classic per-sample dictionaries) a search runs population-level,
  naming the carriers of an indel alone (the index holds them), and `--per-sample`
  falls back to it with a `WARNING` that names the dataset: re-download without
  `--no-genotypes` to get carriers.
- If you only need site-level results, pass `--population-level`. Which mode is faster
  depends on the panel. On the indexes CRISPRitz 2.10.0 builds, which code every allele, a
  population-level search reads no haplotype, so nothing skips the raw targets of the dense
  windows of the immunoglobulin loci, which a per-sample search skips when no sample can have
  them: on HGDP + 1000 Genomes (4,091 samples) the README guide's population-level search of
  chr22 took 23 to 32 min against 8 to 11 min per sample. On a small panel it can be the faster:
  on HPRC (232 samples) 9 to 12 min against 14 to 17 min per sample (16 threads; METHODS.md
  §5). These were measured during 2.10.0's development, on fake-indel windows of 25/27; on the
  release's 30/30 the integration's gate (8 threads a search, with many searches at once) took
  31 against 10 min on HGDP + 1000 Genomes and 9 against 14 min on HPRC, and beside mega 39
  against 17 min (README, "What the default search costs").
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
  20-NRR_3_hg38+hg38_HGDP1kGP` (the recommended phased HGDP + 1000 Genomes default, of the
  default SpCas9 NRR PAM, for a 20-nt guide).
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

On a machine other people use, keep the data folder, the output folders and the temporary
folder (`$TMPDIR`, or `CRISPRME_TMPDIR` when set) writable by you alone (a umask of 022 or
077, or `chmod go-w` on them; a shared temporary folder such as `/tmp` will do while it is
sticky). Anyone who can write a folder a search or the web app works in can do more there
than change the inputs and the results:

- stop the search, by making first a name it makes there (the post-analysis's
  `<output>.pieces` folder, for one);
- rename a folder the search made there and put their own in its place (a rename needs only
  write access to the folder that holds it), and put links in it, which the search then
  writes through, over files of yours elsewhere (the web app follows none: below);
- put lines and names in the files it reads there (its inputs, and the files it writes and
  reads back), which older parts of the pipeline still pass to shell commands; and, in a
  temporary folder that is not sticky, replace the modules TensorFlow writes there and then
  imports.

What holds wherever you run it, since 2.10.0:

- Python looks for no module in the folder CRISPRme is run from: every `python -c`, `-m` or
  `-` it runs is isolated (`-I`), and `--debug` keeps the `PostProcess/` beside the
  `crisprme.py` that runs. Up to 2.10.0 the scorer's health check, which every search runs
  first from the folder you run `crisprme.py` in, loaded a module put there, and `--debug` took
  that folder's `PostProcess/` for its own;
- no file CRISPRme writes is loaded as code (they are JSON or text; TensorFlow's own
  modules, above, and the web app's table cache, in a folder of the server's own, below,
  aside), and the web app's Settings Retry re-runs a job with the command the server
  launched it with, kept in its memory and run without a shell. Up to 2.10.0 it ran through
  a shell the command a job's `Settings/<job>/cmd.json` held;
- the values a review found going from a job's files into shell commands as the shell's own
  text are quoted or checked now: the guides of the guide file, the indel targets' header,
  the names of the CFD graph files, and the guide and the annotation of the results' region
  page;
- the temporary files and folders it makes with `tempfile.mkdtemp`, `mkstemp` or `mktemp -d`
  take names no one else could make first. Where they are made, in the temporary folder or
  (the annotation temp folders, without `CRISPRME_TMPDIR`) in the working directory, only a
  sticky or private folder keeps others from replacing them; and the intermediate files a
  search writes in its output folder keep fixed names (such as `<file>.tmp`);
- the web app follows no link another user put in its data folder. Every file it reads,
  serves or writes there (a search's logs, parameters and results, the downloads of
  `/Results/<path>`, the Settings jobs' logs, the uploads, the email, concurrent searches
  and COSMIC licence settings, the annotation manifest, the sample lists, the PAM files, an
  index's label) is reached one folder at a time with `O_NOFOLLOW` (`openat`), and must be a
  regular file with no other name: not a link, a FIFO, a device or, in a folder others can
  write, a hard link. A folder that one of its files names (the personal-assembly results
  page reads the two haplotypes' folders the combined job's `.Params.txt` names) must be a
  folder of `Results/` named alone: a link, a `..` or a path there is refused. A file it
  writes whole is made under a name of its own (`O_EXCL`) and renamed over its name, so a
  link there is replaced, not written through, and a log is opened with `O_NOFOLLOW`; a
  program it runs on a job's files (the plots, the radar chart, the region page's lookup)
  reads copies of them or its standard input, and the zips are made in Python. A folder it
  makes for files of its own in a folder others can write (a program's outputs in a job's
  folder, an uploaded genome or personal-assembly bundle unpacked in `Genomes/` or
  `Assemblies/`, an annotation sorted) is used through its descriptor alone (Linux's
  `/proc/self/fd`), never by its name again: renamed away, with a link put in its place, its
  name leads nowhere the web app writes. A link is followed only where no one else can have
  put it: a link of yours in a folder only you can write (or a sticky one of yours), such as
  a genome linked to a shared copy, a VCF dataset registered from a server folder, or
  `Results/` or `Settings/` linked onto a volume; in a folder your group can write, the web
  app refuses those too, and below `Results/` and `Settings/` it follows none. A Settings
  job's folder is writable by you alone, and the cache of the results tables, pickles that
  Flask-Caching loads, is in a folder of the server's own in the temporary folder (in memory
  when that folder is not private or sticky), no longer in `Cache/` of the data folder, which
  can be deleted. Up to 2.10.0 the web app followed such links: it showed and served files of
  yours to anyone who could reach it (it listens on every network interface, `0.0.0.0`), the
  results of any search of yours too, through the haplotype folder a combined job's
  `.Params.txt` named; it wrote through them, a chunked upload's data, the email password and
  a job's logs among what it wrote, and through a link put in place of a folder it had made:
  an uploaded archive unpacked into a folder of yours (a bundle's `.ssh/` into your home),
  and every file of a folder of yours moved into a job's `imgs/`, which it then served; it
  copied a JSON file of yours that a link named into the data folder, readable by all, when
  the concurrent searches setting was saved; and another user's pickle in `Cache/` ran as you
  when it was loaded.

What a shared server still needs, with all of the above:

- the search, and the other commands the web app starts (a download, an index build), still
  follow links in the folders they write and write through them, as the list above says: keep
  the data folder, the output folders and the temporary folder writable by you alone, or the
  web app refuses the links of yours there as well;
- the web app has no login: anyone who can reach its port can start searches and download,
  upload and delete data in Settings. On a machine others can reach, let only those you trust
  reach it (a firewall, or Docker's `-p 127.0.0.1:8080:8080` and an SSH tunnel: README,
  Local Web Interface);
- a folder that others can write through an ACL its mode bits do not show (macOS, NFSv4)
  counts as yours: the web app follows your links there.

`test_safe_files` (in the unit-test job) tests how the web app opens a file of its data folder
and uses a folder it makes there, `test_web_no_links` (in the web image's step) the pages that
read, serve and write one, through the real Dash server, each with a link another user put
there or a folder of the web app's renamed and a link put in its place, and
`test_settings_progress` the Settings page's, without Dash. A test (`test_no_code_from_files`,
in the unit-test job and the web image's) keeps out of CRISPRme what a scan of the source can
show: loads that run code (pickle and the like, `eval`, `exec`, YAML loads but `safe_load`), a
temporary name made by name only, and a `python -c`, `-m` or `-` that is not isolated. It also
runs `process_summaries.py` and `assemble_cfd_graphs.py` on a line or a name another user could
put in a search's folder and, in the web image's step, `crisprme.py --debug` from a folder
holding someone else's `PostProcess/`. The other fixes above have tests of their own:
`test_scorer_env` (the scorer's health check, its self-test and the look-up of its env's
python), `test_settings_progress` (the Retry), `test_indel_empty_chromosome` (the indel targets'
header) and `test_web_shell_safety` (the region page). The scan cannot follow every value read
back into a shell command, which is one more reason to keep these folders yours.

---

## 8. Very high `--thread` can stall post-analysis

The post-analysis pool is bounded by `--thread` (and by memory; at most 32 workers unless
`CRISPRME_POSTPROC_MAX_WORKERS` says otherwise). Since 2.10.0 it runs each chromosome in
pieces over those workers, so a dense locus no longer runs on one core
(`CRISPRME_POSTPROC_PIECES=0` runs one process a chromosome, as before). On many-core
nodes, extremely high values (e.g. `--thread 200`) stalled the per-chromosome pool of
earlier releases. A value in the **16–32** range is a good balance on large machines;
raising it further rarely helps and can hurt.

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

**Fix.** Once no CRISPRme container is building, remove them, from your data folder:

```bash
find genome_library -mindepth 1 -maxdepth 1 -name '.build-index-only-*' -exec rm -rf -- {} +
find Genomes -mindepth 1 -maxdepth 1 -name '.build-index-only-enrich-*' -exec rm -rf -- {} +
find Dictionaries -mindepth 1 -maxdepth 1 -name '.build-index-only-dicts-*' -exec rm -rf -- {} +
```

A `find` of each name pattern does the same pasted into bash or into zsh, whichever of the
three folders holds such a folder. Up to 2.8.1 this section gave one `rm -rf` of the three
patterns, which zsh, the default shell of macOS, expands itself: where one of them matched
nothing, it stopped the whole command (`no matches found`) and removed none of them.

## 11. A search stops and says no index serves it

Each of these names the index and what to do; none of them removes anything.

- **`no variant index here serves this search (...)`.** The dataset's indexes are for another PAM or
  fewer bulges, and neither its VCFs nor its enriched genome are here to build one. A download holds
  its own index only: `20-NRR_3_hg38+hg38_HGDP1kGP` searches NRR, and the PAMs NRR holds (NGG, NAG,
  NGA, NRG), pruned to them (and, with a 20-nt guide, up to 2 bulges of each type). Search with what
  the message lists, or install another index ([Precomputed indexes](PRECOMPUTED_INDEXES.md)).
- **`index '...' is not installed, and this run cannot build it`.** The library is on a read-only
  mount, has no write permission or is full, and lacks the index this PAM and these bulges need. Make
  it writable, or run `build-index-only` where it is, or search with no bulges, which opens no index.
- **`... was built for a 20-nt guide, and this search has a 23-nt guide ...: it stores too few bases
  per site`.** This search builds its own index under a name that says its guide length (since
  2.10.0: `23-NGG_2_hg38`), and an index of that name is there, whose trees were built for another
  guide (a folder renamed by hand). Name it for the guide it was built for, as the ERROR says (`mv
  .../23-NGG_2_hg38 .../20-NGG_2_hg38`), or move it aside. An index of the older name
  (`NGG_2_hg38`) built for another guide is left as it is, and the search builds its own beside it.
  Any published index serves a 21-nt guide at up to 1 bulge of each type.

## 12. A build or a search refuses to enrich a dataset anew

**Symptom.** `build-index-only`, or a variant search, stops before it builds anything with
`No finished enrichment of <dataset> is in place (...), so ... would enrich it anew with the
installed CRISPRitz, whose fake-indel windows are 30/30; but an earlier enrichment of it, whose
records tell windows of 25/27 (...), left what was built from it: ...`.

**Cause.** The dataset's enriched genome (`Genomes/<genome>+<dataset>`) is not there (removed to
save space, or never installed: an index download holds none), while its dict-less tiers or
indexes are, built from an enrichment of CRISPRitz 2.9.1 or earlier (CRISPRme 2.8.1, fake-indel
windows of 25/27). A new enrichment of CRISPRitz 2.10.0 differs from it in its SNPs too (2.9.1
enriched only the last ALT of a site written on several records), and a search would read the
stores of the two together, with carriers or rows wrong or lost.

**Fix.** Remove what the message lists (it gives the `rm -rf` command, each path absolute), then
run `build-index-only` with the dataset's VCFs and its `--samplesID`, which enriches it anew and
builds its tiers and index; or put the earlier enrichment back and search it as it is (its windows
hold targets of up to 27 bases, guide, PAM and DNA bulges counted, 26 at an MNV). Beside the stores
of an enrichment of the installed CRISPRitz's windows, a dataset is enriched anew as before.
`build-index-only` run on an enrichment of 2.9.1's windows in place says how to enrich it with
2.10.0's in the same way (METHODS §4). The refusal goes by what the earlier enrichment's records
tell (its indel logs, an INDELS index of it, or a variant index of it in index format 1): with
all of those removed, nothing tells its windows, and its tiers or a variant index of format 2
left in place are not caught.

**A sites-only dataset (`hg38_mega`).** The message reads `... would enrich it anew from its
VCFs; but hg38_mega is sites-only: its VCFs name no sample, and CRISPRitz writes a fake-indel
window only for an indel a sample carries, so an enrichment of them holds none, while its indel
logs or its INDELS indexes record windows of ...`, whatever their width. The mega's windows come
from its enriched-genome kit, whose indels a synthetic sample carries
(`seq_script/merge_panels/mega_rebuild.sh`, the mega's recipe of 2.10.0, builds the whole kit
so): enriched anew from its VCFs, it would lose every indel off-target. Put its kit back
(`Genomes/hg38+hg38_mega` and `_INDELS`: download the kit,
[Precomputed indexes](PRECOMPUTED_INDEXES.md)) and build or search as before; for the
windows of 2.10.0, install its archive and kit of CRISPRme 2.10.0 in place of what the message
lists, or build the kit again that way. Its registry, which comes from its VCFs alone (their
INFO AF fields), is never in that list.

**Building the mega's kit again.** The recipe builds the dataset whole, its VCFs, registry
and indel frequencies too, in a working directory that holds none of it: it refuses one that
holds the dataset's VCFs, or a part it builds, and makes nothing there. So run it in a new
working directory, then put what it built in place of every part of the dataset in the data
folder: the commands below remove each such part by name first, whatever a message's own
`rm -rf` removed or left, namely the kit in place and all that was built from it (its genomes,
dictionaries, indel logs, tiers and every index of the dataset, which the messages name for
removal too) and the dataset's VCFs, registry and indel frequencies, which the recipe wrote
anew. Remove any copy of the dataset's indexes in a library you search with `--index-path`
too. The commands run in bash whatever shell you paste them into, bash or zsh (the
default shell of macOS, which reads a pasted `#` as a word, not a comment, and stops a command
at a pattern that matches nothing): `bash -euo pipefail -s ... <<'EOF'` hands them to bash as
they are written, which stops at the first command that fails (an unset or empty variable
among them): nothing in the data folder is removed or moved unless the recipe has built the
whole dataset and passed its checks, and every part they move in is there. A mistyped `SRC`,
a missing chromosome FASTA or a disk that fills during the recipe stops them with the data
folder as it was; remove `NEW`, which the recipe refuses once it holds a part, and run them
again. Set, in the shell you paste them into (no `export` needed: they are passed to bash as
its arguments), `DATA` the data folder, `NEW` a new folder outside it on the same file system
(so that each move is a rename, which needs no space), `SRC` a folder of the merged VCFs the
mega is built from, `mega.<chrom>.afmax.vcf.gz`, each an absolute path (`crisprme.py download
--what vcf --dataset hg38_mega --path <dir>` fetches them into `<dir>/VCFs/hg38_mega`, `<dir>`
being another folder than `NEW`, which the recipe refuses once it holds the dataset's VCFs; a
2.8.1 install's `VCFs/hg38_mega` holds them, so you can give it as `SRC`: it is removed only
after the recipe has built the dataset from it), and `REPO` the CRISPRme repository of 2.10.0,
whose `PAMs/` holds the recipe's PAM file (a data folder that CRISPRme 2.8.1 set up has none of
NRR). Run the commands one of two ways. In a shell inside the image (`docker run -it ...
bash`, or `apptainer shell`), `REPO` is `/opt/conda/opt/crisprme` and the recipe's `RUN` prefix
stays empty. From the host, with the recipe's programs run in the container, export `RUN` first,
the container prefix the recipe's header describes, and set `REPO` to a checkout of the
repository of 2.10.0 on the host: the recipe runs the scripts of `$REPO/PostProcess` through
`RUN`, so the container must see `REPO`, `DATA`, `NEW` and `SRC` at the paths the host gives
them (bind them so). Set to the image's path on the host, `REPO` stops the commands at their
first `cp`, the data folder unchanged:

```bash
bash -euo pipefail -s "$DATA" "$NEW" "$SRC" "$REPO" <<'EOF'
# bash, whichever shell this is pasted into: it stops at the first command that fails, so the
# data folder changes only after the recipe has built the whole dataset
DATA=${1:?DATA is not set} NEW=${2:?NEW is not set} SRC=${3:?SRC is not set} \
  REPO=${4:?REPO is not set}
# a new working directory: the FASTAs of the chromosomes the mega is built on and no other
# (chr1-22, X and Y), and the PAM file, the repository's
mkdir -p "$NEW/Genomes/hg38" "$NEW/PAMs"
for c in chr{1..22} chrX chrY; do ln -s "$DATA/Genomes/hg38/$c.fa" "$NEW/Genomes/hg38/"; done
cp "$REPO/PAMs/20bp-NRR-SpCas9.txt" "$NEW/PAMs/"
# (its input /dev/null: bash reads these commands from its own input, the recipe must not)
PAM=PAMs/20bp-NRR-SpCas9.txt bash "$REPO/seq_script/merge_panels/mega_rebuild.sh" "$NEW" "$SRC" \
  < /dev/null
# each part the recipe built that the moves below take, there before anything is removed
for p in Genomes/hg38+hg38_mega Genomes/hg38+hg38_mega_INDELS Dictionaries/log_indels_hg38_mega \
  Dictionaries/registry_hg38_mega Dictionaries/indel_af_hg38_mega VCFs/hg38_mega; do
  [ -d "$NEW/$p" ] || { echo "no $NEW/$p: the data folder is left as it was" >&2; exit 1; }
done
ls -d "$NEW"/genome_library/*_hg38+hg38_mega "$NEW"/genome_library/*_hg38+hg38_mega_INDELS
mkdir -p "$DATA/Genomes" "$DATA/genome_library" "$DATA/Dictionaries" "$DATA/VCFs" "$DATA/samplesIDs"
# in the data folder, every part of the dataset the recipe built (each it refuses where it
# runs): the kit in place and all that was built from it, then the dataset's VCFs, registry
# and indel frequencies, which the recipe wrote anew
rm -rf "$DATA/Genomes/hg38+hg38_mega" "$DATA/Genomes/hg38+hg38_mega_INDELS" \
  "$DATA/Dictionaries/dictionaries_hg38_mega" "$DATA/Dictionaries/log_indels_hg38_mega" \
  "$DATA/Dictionaries/genotypes_hg38_mega" "$DATA/Dictionaries/indel_genotypes_hg38_mega" \
  "$DATA/VCFs/hg38_mega" "$DATA/Dictionaries/registry_hg38_mega" \
  "$DATA/Dictionaries/indel_af_hg38_mega"
# and the dataset's indexes, of every PAM (<PAM>_<N>_hg38+hg38_mega and its _INDELS), a find of
# each name pattern: none of them, or one kind alone, is no error (a hidden folder is a build's)
find "$DATA/genome_library" -mindepth 1 -maxdepth 1 -name '*_[0-9]_hg38+hg38_mega' \
  -not -name '.*' -exec rm -rf -- {} +
find "$DATA/genome_library" -mindepth 1 -maxdepth 1 -name '*_[0-9]_hg38+hg38_mega_INDELS' \
  -not -name '.*' -exec rm -rf -- {} +
# what the recipe built, moved in (not its reference index, genome_library/<guide>-<PAM>_<N>_hg38)
mv "$NEW/Genomes/hg38+hg38_mega" "$NEW/Genomes/hg38+hg38_mega_INDELS" "$DATA/Genomes/"
mv "$NEW"/genome_library/*_hg38+hg38_mega "$NEW"/genome_library/*_hg38+hg38_mega_INDELS \
  "$DATA/genome_library/"
mv "$NEW/Dictionaries/log_indels_hg38_mega" "$NEW/Dictionaries/registry_hg38_mega" \
  "$NEW/Dictionaries/indel_af_hg38_mega" "$DATA/Dictionaries/"
mv "$NEW/VCFs/hg38_mega" "$DATA/VCFs/"
[ -e "$DATA/samplesIDs/hg38_mega.samplesID.txt" ] \
  || mv "$NEW/samplesIDs/hg38_mega.samplesID.txt" "$DATA/samplesIDs/"
EOF
```

The recipe's last step, `verify_mega_chrom.py`, prints `ALL PASS` for every chromosome. The
index it built is of its PAM; for another, run `build-index-only --pam <it>` on the data folder,
which reuses the kit. The name patterns of the `find` take an index of up to 9 bulges of a type
(`_[0-9]_`, as the messages' own removal does): one built with `--bDNA` or `--bRNA` of 9 or more
(`_10_` and up), which no build of the web or of the README makes, is removed by hand. Each
message's removal and these commands do the same pasted into bash or into zsh.

**A message that says to remove the mega's enriched genome.** Where a search finds the mega's
INDELS index incomplete, or of other windows than its indel logs, and the fake-indel genome
(`Genomes/hg38+hg38_mega_INDELS`) is not in place while the dataset's VCFs are, its ERROR ends
`Rebuild it: rm -rf ... <data>/Genomes/hg38+hg38_mega && crisprme.py build-index-only ...` (or
`Rebuild them:`). For `hg38_mega` do not run it: `build-index-only` refuses the sites-only
dataset (above), after that `rm` has removed its enriched genome. Put the kit back instead (its
download, or the recipe above).

## 13. A build or a search says a fake-indel genome and its indel logs come from two enrichments

**Symptom.** `build-index-only`, before it builds anything, or a variant search that would
build an index from the enrichment in place (the INDELS index it lacks beside a variant index in
place, say), before it puts the SNP codes on the fake-indel genome, stops with
`The fake-indel genome Genomes/<genome>+<dataset>_INDELS and the indel logs ... come from two
enrichments: on <chrom>, the window of ... ends elsewhere there`; or `download --what index`
warns that `the enrichment of <dataset> in place ... is another one than the archive's`. Nothing
is changed: the SNP overlay, which places each code by the logs, refuses such a genome too.

**Cause.** The dataset's indel logs, dictionaries and tiers came one way (a download of an index
installs them) and its enriched and fake-indel genomes another (an enrichment made here, by
another CRISPRitz, whose windows were another width). An index built from the enrichment in place
would be read with the other enrichment's logs: every indel row misplaced.

**Fix.** Put the enrichment that goes with the logs in place of this one: remove it, as the
message says (`rm -rf` of the enriched and the fake-indel genome), and install the published
index's enriched-genome kit ([Precomputed indexes](PRECOMPUTED_INDEXES.md)); or remove all of the
dataset that comes from an enrichment (the message lists it) and run `build-index-only` with its
VCFs and `--samplesID`, which enriches it anew. For a sites-only dataset (`hg38_mega`) the second
way is to build its kit again with a synthetic sample (`seq_script/merge_panels/mega_rebuild.sh`,
the enriched genome, the fake-indel genome and their indexes in one build), as the message says:
an enrichment of its VCFs would hold no fake-indel window. The recipe runs in a new working
directory, and what it built goes in place of all of the dataset that comes from an enrichment
(the message lists it) and of the dataset's VCFs, registry and indel frequencies, which it
writes anew: the commands, which remove each of those parts by name once the recipe has built
the whole dataset, and change nothing in the data folder where it fails, are in section 12,
"Building the mega's kit again" (they run in bash, pasted into bash or into zsh). The `rm -rf`
of the messages, each path quoted, does the same in either shell.

## 14. A search above 1 DNA + 1 RNA bulge warns, then runs for hours or runs out of memory

**Symptom.** `complete-search` prints `WARNING [complete-search]: 4 mismatch(es) + 2 DNA + 1 RNA
bulge(s) on a variant index that codes every allele (index format 2, CRISPRitz 2.10.0). Above
the 1 DNA + 1 RNA default, ...` (its `.Params.txt` keeps it, and its `log_verbose.txt` says it
again, a web search's too); the search then runs for many hours, or is killed for memory (exit
137, an empty `log_error.txt`).

**Cause.** The indexes of 2.10.0 code every allele of a split site, so the densest windows of
the variant genome, in the immunoglobulin loci (IGL on chr22, IGH on chr14, IGK on chr2), spell
millions of sequences, and each bulge beyond one of a type multiplies the alignments the search
writes there: on IGL alone, at the default `--max-total-edits` of 6, 0.76 million SNP raw
targets at 4 + 1 + 1 and 33.7 million at 4 + 2 + 2 (28 GB of memory), where the indel search
wrote 204 million (157 GB). No step after the search removes that cost.

**Fix.** Keep 1 DNA + 1 RNA bulge, the default, and raise `--mm` first for a wider search; or
keep `--max-total-edits` at 3 or less, which bounds the alignments; or search the reference
alone, or a chromosome that holds none of those loci. CRISPRme 2.10.1 is planned to search, in
those regions, the samples' own haplotypes, indels included, in place of the variant genome's
combinations (the haplotype scan); the indexes of 2.10.0 serve it as they are.

## 15. A row that looks wrong: the known limits of 2.10.0

Each is rare and changes no other row; [METHODS.md](../METHODS.md) §4 gives the rules behind
them.

- **A carrier listed for a SNP that its own target does not hold.** On an indel row, a carrier
  whose haplotype holds a second indel that pushes one of the row's named SNPs out of the
  target's bases is listed on the row all the same. On 1000 Genomes chr22, a per-sample search of
  seven guides at 4 mismatches + 2 DNA + 2 RNA bulges lists one such carrier (NA20515, at VEGFA
  site 1), and none of the 17 gate searches of 2.10.0 does. The haplotype scan planned for
  2.10.1 reads each carrier's own bases.
- **One target written as two rows, or twice, in a microsatellite.** In a tandem repeat that
  holds a whole target, the merge places each row by its own bases, and the target's
  alignments can fit different places of the repeat: they can then go to two rows, or the
  target be written twice (in a model of microsatellites at 1 + 1 bulges and `--merge` 3, 134
  targets of 10,000 split and 15 written twice). No row is lost or wrong, and no search of real
  data checked has such a target. To be fixed in 2.10.1.
- **An unphased carrier of two ALTs on a row that needs the reference there.** A carrier whose
  call is unphased (`/`) is listed, PUTATIVE, on every subset of its variants a target can hold.
  For a call of two ALTs at one position (`1/2`), that includes a row that needs the reference
  base there, which its genotype rules out: both of its slots hold an ALT. The rows are the same
  for a line of several ALTs and for its split form.
- **On a classic install, an unphased carrier only on the row of all its variants.** A
  per-sample search of an index built without `--samplesID` (a classic install, as the
  benchmark `complete-test` builds its own) lists a carrier with an unphased call in a phased VCF
  only on the row of all its variants of a window together, where a dict-less install lists it
  on the row of each subset too. Build the index with `--samplesID`, or download a published
  one.
