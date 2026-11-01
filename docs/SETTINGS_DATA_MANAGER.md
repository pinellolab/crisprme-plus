# Settings / Data Manager — add data from the browser

CRISPRme's web interface has a **Settings / Data Manager** page where you add
everything a search needs — reference genomes, indexes, variant (VCF) datasets,
annotations and PAMs — without touching the command line. Whatever you add lands
in your local data folder and shows up in the search form automatically.
(For exactly where everything is stored — the shared inputs reused by every search
vs. where each run's results go — see [`DATA_LAYOUT.md`](DATA_LAYOUT.md).)

Open it from the **"Settings / Data Manager"** button on the home page, or go to
`/settings` directly (e.g. `http://localhost:8080/settings`).

> Data management runs only in **local mode** (how you run CRISPRme on your own
> machine / Docker). On the public web server the page is read-only, because it
> would change shared data.

## What you can do

### Add a reference genome
Type an assembly name and pick a source:
- **UCSC** — downloads the assembly by name from UCSC (e.g. the pig genome
  `susScr11`). CRISPRme handles both UCSC layouts (per-chromosome archive or a
  single multi-FASTA) automatically: a single multi-FASTA download is split into one
  file per sequence, which is the layout every search requires.
- **HuggingFace** — pulls a genome we host on the CRISPRme data repo.
- **Direct URL** — any `.fa.gz` / `.tar.gz` link.
- **Upload** — pick a local genome file (`.fa`, `.fa.gz` or `.tar.gz`); it is sent
  in chunks, so there is no browser size limit. It is laid out as a Direct URL
  download is: every FASTA of an archive goes to the top of the genome folder (a UCSC
  `chromFa.tar.gz` holds them in `chroms/`), and a FASTA holding several sequences is
  split into one file per sequence. An upload with the name of an installed genome
  replaces it.

The same thing from the CLI:
```bash
crisprme.py download --what genome --ref susScr11 --source ucsc
```

### Add an index (or build one)
Bulge-enabled searches need a genome **index**.
- **Download** a ready-made index from HuggingFace by name (e.g. `NRG_3_hg38`).
- **Build locally** from an installed genome + PAM + bulge counts.
- **Pre-build a variant-aware index**: pick a VCF dataset in the build panel and
  CRISPRme enriches the genome with the variants and indexes it ahead of time, so
  your first variant-aware search is fast instead of slow. The build passes the
  dataset's sample list (`samplesIDs/<dataset>.samplesID.txt`, see *Add a VCF
  dataset*) as `--samplesID` when it names samples the VCFs genotype (a sites-only
  VCF, such as the `hg38_mega` ones, has none and is built without it), so it also
  emits the dict-less Tier-0 registry + Tier-1 genotype tiers (and, for a merged panel, the
  combined samplesID list): per-sample analysis and the dataset's data type on the
  search form need them, and the index is self-complete for fast post-analysis. (CLI
  equivalent:
  `crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt --bDNA 2 --bRNA 2 --vcf VCFs/hg38_HGDP1kGP --samplesID samplesIDs.config.txt`.)
- A build stops with an error, and puts no index in place, when the genome folder
  holds no `.fa` file at its top, holds a FASTA of several sequences, or gives an
  index with no `.bin` file: every search would refuse such an index.

### Add a VCF dataset
Variant datasets are large, so they are **fetched server-side** rather than
uploaded from the browser: from HuggingFace, an existing server folder already on
the machine (registered in place), or a chunked upload for a local file. Pick the
reference genome the VCF is called against on the card first: the search form
offers the dataset for that genome, whatever the dataset's name. Your VCF does
**not** need special preparation — a single multi-chromosome `.vcf.gz` is
**auto-split per chromosome** when you build an index from it (its `.tbi` is
optional: CRISPRme indexes a bgzip or plain-gzip file itself), and a SNP-only VCF
(no indels) builds a complete SNP variant index.

A search reads the dataset's samples from `samplesIDs/<dataset>.samplesID.txt`
(tab-separated `#SAMPLE_ID POPULATION_ID SUPERPOPULATION_ID SEX`). When the dataset
has none, adding it here writes one from the VCF header: exactly the samples the VCF
genotypes, with the dataset as their population and `unknown` superpopulation and
sex. To give their populations and sex, put your own list there before you build
the index; a list you put there is never changed. The search form refuses a variant
search whose sample list is missing, rather than letting it fail at its end.

### Add an annotation
Upload a `.bed` file (annotations are small), then enable it under *Manage
annotations*. That card offers a BED for the genome its file name names (e.g.
`mytrack.hg38.bed` for hg38), and a BED whose file name names no installed genome
for every genome: enable it for the one its coordinates are on.

### Define a nuclease / PAM
Fill in a short form — name, guide length, motif (e.g. `NGG`, `TTTV`), and whether
the PAM is at the 3′ (Cas9-like) or 5′ (Cas12a-like) end — and CRISPRme writes the
PAM file for you.

## Storage and cleanup
The right-hand **Storage** panel shows how much disk each category
(genomes / VCFs / indexes / annotations) is using, the total, and how much free
disk you have. Each download card also states the space it typically needs, so
you can check before starting.

The **Remove installed data** card lets you delete any downloaded genome, index,
VCF dataset or annotation to free space — with a confirmation prompt. Anything you
delete can be downloaded again later.

## Command-line equivalents (everything the web does)

The web app is a GUI over the same `crisprme.py` commands — **anything you can do in
the browser you can do on the command line.** The data all lands in the same layout
(`Genomes/`, `genome_library/`, `VCFs/`, `Annotations/`, `PAMs/`, `Assemblies/`), so a
CLI-built install and a web-built install are interchangeable.

| In the web | On the command line |
|---|---|
| Run an off-target **search** | `crisprme.py complete-search --genome Genomes/<g> --pam PAMs/<p>.txt --guide guides.txt --mm <n> --bDNA <n> --bRNA <n> [--vcf list_vcf.txt --samplesID list_samplesID.txt] --annotation <bed> --output <name>` |
| **Personal-assembly search** | `crisprme.py assembly-search --genome-paternal … --genome-maternal … --chain-paternal … --chain-maternal … --chrom-alias-paternal … --chrom-alias-maternal … --guide guides.txt --pam PAMs/<p>.txt --mm <n> --bDNA <n> --bRNA <n> --output <name>` |
| **Download a reference genome** (UCSC) | `crisprme.py download --what genome --ref <assembly> --source ucsc` |
| **Download a precomputed index** | `crisprme.py download --what index --index-name <name>` |
| **Download a VCF dataset** (HuggingFace) | `crisprme.py download --what vcf --dataset <name>` — the raw source VCFs; the default repo hosts `hg38_HGDP1kGP`, `hg38_HPRC` and `hg38_mega`. A search needs only the prebuilt **index** (use *Download a precomputed index*); the VCFs are for building one yourself |
| **Fetch an individual from HPRC** | `crisprme.py download --what assembly --sample <id>`  *(e.g. `HG01255`)* |
| First-run **reference data + annotations** | `crisprme.py download --what all` |
| **Create (build) an index** — reference or variant | `crisprme.py build-index-only --genome Genomes/<g> --pam PAMs/<p>.txt --bDNA <n> --bRNA <n> [--vcf VCFs/<dataset>]` |
| **Add a VCF dataset** (register a folder) | place the `.vcf.gz` folder under `VCFs/<name>/` and pass `--vcf VCFs/<name>` (a single multi-chromosome VCF is auto-split; SNP-only is fine) |
| **Add an annotation** (BED) | place the `.bed` under `Annotations/` and pass `--annotation Annotations/<file>.bed` to a search |
| **Define a nuclease / PAM** | write the one-line PAM file under `PAMs/<name>.txt` (`<IUPAC_motif> <signed_position>`; see `docs/INPUT_FORMATS.md`) and pass `--pam PAMs/<name>.txt` |
| **COSMIC licence** enable/disable | `crisprme.py cosmic-license {enable|disable|status}` |
| CRISPR-Bulge **scorer environment** | `crisprme.py scorer-env {create|check|update|list|doctor}` |
| **Publish an index** to HuggingFace | `crisprme.py publish-index --index genome_library/<name> [--dictless]` |
| **Remove installed data** | delete the folder/file (e.g. `rm -r genome_library/<name>` · `Genomes/<name>` · `Assemblies/<name>`) |
| Shareable **report** | produced by each search; regenerate with `crisprme.py generate-report …` |

Run any subcommand with `--help` for its full option list.

## How long operations behave
Downloads and index builds run in the background on a dedicated slot (so they
never slow down searches) and report **live progress** in the Activity panel.
When one finishes, the installed-data and storage panels update automatically and
the new data is immediately available in the search form.
