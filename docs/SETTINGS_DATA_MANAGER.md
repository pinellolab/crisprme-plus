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

A name you give on this page (a genome's, a VCF dataset's, a nuclease's, an HPRC sample's),
the name of a file you upload here (but a personal-assembly bundle's, which is named after its
top folder) and the name a Direct URL ends in may hold no blank and none of `|`, `;`, `'`, `"`,
`$`, `(`, `)`, `&`, `<`, `>`, a backquote, `#`, `*`, `?` or `[`, and may not start with `-` or
`.`: CRISPRme takes these names into shell commands as they are. The page refuses such a name,
and says why. A personal-assembly bundle is refused, saying why, when the name of a FASTA file
of its genome folders holds a blank or any of `|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`, `>`,
a backquote or `#`: a search takes that name, as its chromosome's, into shell commands too.

### Add a reference genome
Type an assembly name and pick a source:
- **UCSC** — downloads the assembly by name from UCSC (e.g. the pig genome
  `susScr11`). CRISPRme handles both UCSC layouts (per-chromosome archive or a
  single multi-FASTA) automatically: a single multi-FASTA download is split into one
  file per sequence, which is the layout every search requires.
- **HuggingFace** — pulls a genome we host on the CRISPRme data repo.
- **Direct URL** — any `.fa.gz` / `.tar.gz` link. Its FASTA files go to the top of the
  genome folder, a FASTA of several sequences is split into one file per sequence, one of
  one sequence is named after it, and each header is written as `>id` (as for an upload,
  below). An archive with two FASTA files of one name (`a/chr1.fa` and `b/chr1.fa`) is
  refused, as by an upload, and nothing is installed: if it holds two genomes (the two
  haplotypes of an assembly), install each under a name of its own. A download replaces a
  genome installed under its name only once it is laid out: one that is refused leaves
  that genome as it was. As for an upload, a `.fasta` file goes in as a `.fa` file (the
  files an index and a search read), so `chr1.fa` and `chr1.fasta` are two of one name, and
  the hidden files of an archive are left out (the `._chr1.fa` that macOS `tar` adds beside
  `chr1.fa`). Before 2.10.0 the two went in as they were.
  Like an upload, it is refused while indexes or variant genomes built on a genome of that
  name are installed, and the page says so before anything is fetched: the link can hold
  another genome than the one they were built from. Download it under another name, or
  remove them first.
- **Upload** — pick a local genome file (`.fa`, `.fa.gz` or `.tar.gz`); it is sent
  in chunks, so there is no browser size limit. It is laid out as every search needs
  it: every FASTA of an archive goes to the top of the genome folder (a UCSC
  `chromFa.tar.gz` holds them in `chroms/`) as a `.fa` file (`chr1.fasta` as
  `chr1.fa`: an index is built from the `.fa` files only), a FASTA holding several
  sequences is split into one file per sequence, and a FASTA of one sequence is named
  after it (`hg38_chr22.fa` holding `>chr22` goes in as `chr22.fa`: a search takes each
  file for the sequence it is named after). Each header is written as `>id`, with no
  description, and every line with no carriage return, the last one too: an NCBI, GENCODE
  or Ensembl FASTA
  (`>NC_000913.3 Escherichia coli ...`, `>chr1 1`) goes in as `>NC_000913.3`, `>chr1`. An
  id the pipeline cannot take (one that holds `|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`,
  `>`, a backquote or `#`, as `>gi|49175990|ref|NC_000913.2|` does, or that starts with `-`
  or `.`) is refused, by the upload when it writes that file and otherwise by the build:
  rename the sequence and upload it again.
  Hidden files are left out (the
  `._chr1.fa` that macOS `tar` adds beside `chr1.fa` holds its file attributes, not a
  sequence). An archive with two FASTA files that would go in under one name
  (`a/chr1.fa` and `b/chr1.fa`, `chr1.fa` and `chr1.fasta`, or two files of one
  sequence) is refused. An upload with the name of an installed genome replaces it.
  It is refused while indexes or variant genomes built on a genome of that name are
  installed: every search would go on using them, with the sequences and positions of
  the genome they were built from.
  Remove them first (an index under *Remove installed data*, a
  `Genomes/<genome>+<dataset>` variant genome from your data folder), or upload the
  genome under another name.

The same thing from the CLI:
```bash
crisprme.py download --what genome --ref susScr11 --source ucsc
```

### Add an index (or build one)
Bulge-enabled searches need a genome **index**.
- **Download** a ready-made index from HuggingFace by name (e.g. `NRR_3_hg38`, of the default PAM).
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
  `crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt --bDNA 2 --bRNA 2 --vcf VCFs/hg38_HGDP1kGP --samplesID samplesIDs.config.txt`.)
  The tiers take memory: at panel scale ~150-180 GB a chromosome, by default up to 4
  chromosomes at once (`CRISPRME_TIER_WORKERS` in the container's environment sets how
  many). A worker that runs out sends the chromosomes left to one at a time. If one dies
  then too, the build fails with a message: give the container more memory and build
  again, which writes the rest, or follow the message to build from the command line
  without the tiers (this page has no way to). Until one of those is done, searches of
  the dataset are refused.
- A build stops with an error, and puts no index in place, when the genome folder
  holds no `.fa` file at its top, holds a FASTA of several sequences, or gives an
  index with no `.bin` file: every search would refuse such an index. It stops as
  well for a FASTA of one sequence not named after it (`hg38_chr22.fa` holding
  `>chr22`), whose off-targets every search would find and then drop: rename it after
  its sequence (`python PostProcess/genome_layout.py --split <genome folder>` renames
  them all) and build again. It stops as well for a header that is more than `>id`, or
  lines that end in a carriage return anywhere in a file, in a folder put there by hand or
  installed before 2.8.1 (a search stops too, for a header, or carriage returns in a file's
  first two lines): every search would lose results, often all of them. One that only
  ends a file (on its last line, or after it) is no exception: with a PAM site in the last
  bases of a sequence, the build warned, and a search of the index, or of a variant index
  built on the genome, could crash or report other rows. The
  same command writes them as `>id`; then remove each index built from the folder as it
  was. It stops as well, and so does a search, for a sequence id that holds `*`, `?` or
  `[` and, read as a shell pattern, matches the name of another sequence of the folder
  (`HLA-A*02:10` matches `HLA-A*24:02:10`, as 10 of the HLA sequences of the GATK GRCh38
  FASTA match others): a search, which takes each name into shell commands as it is, would
  fail or lose that sequence's off-targets. Rename such a sequence, its header and its
  file, with none of those characters, or remove it; then remove each index built from the
  folder and build again. It stops too
  when the reference index already there holds no `.bin` (one built before from such
  a folder), or a `.bin` named after more than an id: remove that index and build again.

### Add a VCF dataset
Variant datasets are large, so they are **fetched server-side** rather than
uploaded from the browser: from HuggingFace, an existing server folder already on
the machine (registered in place), or a chunked upload for a local file. An upload
takes a dataset's `.vcf.gz` files and their `.tbi` or `.csi` index, into a dataset
added by upload or download: a server folder registered as a dataset is read where
it is, and nothing is uploaded into it (copy the file into that folder). Pick the
reference genome the VCF is called against on the card first: the search form
offers the dataset for that genome, whatever the dataset's name. Your VCF does
**not** need special preparation — a single multi-chromosome `.vcf.gz` is
**auto-split per chromosome** when you build an index from it (its `.tbi` is
optional: CRISPRme indexes a bgzip or plain-gzip file itself), and a SNP-only VCF
(no indels) builds a complete SNP variant index.

A search reads the dataset's samples from `samplesIDs/<dataset>.samplesID.txt`
(tab-separated `#SAMPLE_ID POPULATION_ID SUPERPOPULATION_ID SEX`). When the dataset
has none, adding it here (an upload or a server folder), or building its index here,
writes one from the VCF header: exactly the samples the VCF genotypes, each with the
population, superpopulation and sex an installed sample list gives it (a downloaded
panel's, or one you wrote), and with the dataset as its population and `unknown`
superpopulation and sex when no list gives it a superpopulation. Each time you add a
file to the dataset or build its index here, the list is written again if that gives
another one (a file with other samples, a sample list installed or removed since),
and it goes when you remove the dataset here. To give its samples' populations and
sex yourself, put your own list there before you build the index (or edit the one
written here): a list you put there or edit is never changed. A list put there, or
written there again, once the variant index is built reaches its dict-less tiers only
after you remove `Dictionaries/registry_<dataset>` and
`Dictionaries/genotypes_<dataset>` from your data folder and build the index again
here: until then a build leaves the tiers as they are, from the list they were built
with. The search form refuses a variant search whose sample list is missing, rather
than letting it fail at its end, and the note under its bulge options says so before
you submit, with how to get one. A variant index built here stays when you remove its
dataset here, and serves a search again once you add the dataset again under the same
name, which writes the list. Until then the search form offers the dataset, and refuses
its search, if the dataset's name begins with its genome's and an underscore (as
`hg38_<name>`), and does not offer it otherwise (as `my_pig_vcf`). A sites-only
VCF (no sample columns) gets no list: a search reports its variants only through a
per-dataset allele-frequency registry, as the shipped `hg38_mega` index has.

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
| **Publish an index** to HuggingFace: command line only (the web app does not publish; in maintainer mode Settings says so) | `crisprme.py publish-index --index genome_library/<name> [--dictless]` |
| **Remove installed data** | delete the folder/file (e.g. `rm -r genome_library/<name>` · `Genomes/<name>` · `Assemblies/<name>`) |
| Shareable **report** | produced by each search; regenerate with `crisprme.py generate-report …` |

Run any subcommand with `--help` for its full option list.

## How long operations behave
Downloads and index builds run in the background on a dedicated slot (so they
never slow down searches) and report **live progress** in the Activity panel.
When one finishes, the installed-data and storage panels update automatically and
the new data is immediately available in the search form.
