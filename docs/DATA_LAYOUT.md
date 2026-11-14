# Where CRISPRme+ stores your data

CRISPRme+ reads and writes **everything under a single data folder that you choose and own** on
your own machine. Nothing important is kept inside the container — so you can always browse, back
up, move, or share any part of it with your normal file manager.

## One folder — you pick it, you own it

| How you run it | The data folder is… |
|---|---|
| **CLI** | `--path <DIR>` on every command (defaults to the current directory). All data + results land under `<DIR>`. |
| **Docker** | `docker run -v "$PWD:/DATA" -w /DATA …` — your host `$PWD` is mounted as `/DATA` in the container. |
| **Apptainer / Singularity** | `apptainer run --bind "$PWD:/DATA" --pwd /DATA …` — same idea. |
| **One-click app / web** | asks *where to store the data* on first launch (a ~70 GB download → ~150 GB on disk; keep ~200 GB free — a reference-only setup is far smaller, ~44 GB) and **remembers** it. |

Because it is a normal host directory, you have full read/write access to it at all times.

## What's inside: shared inputs vs. per-search results

```
<your data folder>/                 e.g.  ~/crisprme     (in the container: /DATA)
│
├─ SHARED INPUTS ── downloaded or built ONCE, reused by EVERY search ──────────────
│  ├─ Genomes/          reference-genome FASTAs (hg38/, mm39/, …)        download --what genome
│  ├─ genome_library/   the search INDEXES — reference + variant-aware    download --what index  /  build-index-only
│  ├─ VCFs/             raw variant datasets (one folder per dataset)      your own, or download --what vcf
│  ├─ Dictionaries/     variant registry + per-sample genotype tiers      (companions to the variant indexes)
│  ├─ samplesIDs/       sample-ID lists for the variant panels
│  ├─ Annotations/      functional + cancer annotations (DHS, ENCODE SCREEN, GENCODE, IntOGen, COSMIC)
│  ├─ PAMs/             nuclease PAM definitions (NRR, the default; NGG, NRG, NNGRRT, …)
│  ├─ Assemblies/       personal / HPRC diploid assemblies, one folder per individual   download --what assembly
│  ├─ LiftoverFiles/    chain + chromAlias files used to reconcile a personal assembly to GRCh38
│  ├─ haplotype_genomes/     one link per personal-assembly haplotype, named after it      made by assembly-search
│  └─ haplotype_alignments/  the alignments of each individual's two haplotypes          made by assembly-search
│
└─ YOUR RESULTS ── one folder per search run; this is where OUTPUTS go ─────────────
   └─ Results/
      ├─ <job>/                               a complete-search run
      └─ <job>_paternal / _maternal / _combined   an assembly-search run (per-haplotype + merged)
```

**Shared inputs** are your reusable "library": download or build them once and *every* search points
at them. They are the large part (the indexes under `genome_library/` are the biggest — tens of GB).

Searches add to them, once: a search builds the index it needs in `genome_library/` when it is
missing (not in a library given with `--index-path`, which a search only reads: a reference index
from there, and a variant index from there first, then from `genome_library/`). A variant dataset's
dictionaries (`Dictionaries/`) come with its index download (`download --what index`), or are built
with its enriched genome (`Genomes/<genome>+<dataset>/`) and its index by `build-index-only --vcf`;
`complete-search --vcf` refuses a dataset that has none yet and says how to get them. From the web,
a search without bulges of a dataset that has only its VCFs (the form then offers no bulges)
enriches the genome and builds the dictionaries and the index of its indels itself (a search
without bulges scans the variant genome, and builds no index of it). A personal-assembly
search (`assembly-search`) builds an index for each haplotype, named after the individual and the
haplotype plus a hash of the haplotype's folder (e.g. `genome_library/20-NRR_2_HG01255_paternal_5aeb0fc6`
and `20-NRR_2_HG01255_maternal_067bdb25`, each the size of a reference index; before 2.10.0 named
without the guide length, `NRG_2_HG01255_paternal_5aeb0fc6`), searches each haplotype
through a link of that name in `haplotype_genomes/`, and keeps the alignments of the two haplotypes
in `haplotype_alignments/`. Deleting an individual (Settings → Data Manager → Delete → *Personal
assembly*) also removes its haplotype indexes and links; the alignments stay, and
`haplotype_alignments/` can be deleted by hand: a later search makes them again.

**Results** is where *your* outputs are written — **one subfolder per run**, named from the job id and
your `--output` name. Each holds the shareable **`report.html`** + its `data/` folder, the off-target
tables, and the run logs (`log.txt`, `log_verbose.txt`, `log_error.txt`). Deleting a `Results/` subfolder
never touches your shared inputs.

## What a search reads vs. what it writes

| A search **READS** (shared inputs) | A search **WRITES** |
|---|---|
| a genome (`Genomes/`) + its matching index (`genome_library/`), a PAM (`PAMs/`), and — for a variant-aware search — a dataset (`VCFs/` + `Dictionaries/` + `samplesIDs/`) and annotations (`Annotations/`). For a personal-genome search: an individual (`Assemblies/` + `LiftoverFiles/`). | its outputs **only** under `Results/<job>/`; and, once, the shared inputs it needs that are missing: the index (`genome_library/`), and for a personal-genome search each haplotype's index, link and alignments (`genome_library/`, `haplotype_genomes/`, `haplotype_alignments/`). A variant dataset's dictionaries come from its download or from `build-index-only` before its first `complete-search --vcf` (a web search can build them: see above). |

So the rule is simple: **inputs are shared and reused (a search adds to them only what it is
missing, once: see above); every run's outputs go under `Results/`.**

## Managing the folder
- **Move everything:** move the whole data folder; nothing is stored elsewhere.
- **Free space / add data:** use the web **Settings → Data Manager** (download, build, delete, with a
  size-aware confirm) or the command-line equivalents — see the CLI↔web table in
  [`SETTINGS_DATA_MANAGER.md`](SETTINGS_DATA_MANAGER.md).
- **Back up:** the shared inputs can always be re-downloaded/rebuilt, so the irreplaceable part is
  your `Results/` — back that up.

> *Transient only:* while fetching from HuggingFace the downloader may use a short-lived cache; the
> final extracted data always lands in your data folder.
