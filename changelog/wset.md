<!-- wset track (web Settings / data manager): [2.8.1] entries, merged by the docs track. -->

**For the docs track: two amendments to an entry already in `## [2.8.1]`, then new entries.** A
variant build in the web Settings page now passes the dataset's sample list as `--samplesID` (the
first entry under Fixed below), and the ERROR that `build-index-only` and the search script give for
tiers still marked unfinished now says so: run `build-index-only` again to its end "with <dataset>
and its --samplesID (or build it again in the web Settings page, which passes the dataset's sample
list)". The entry **"A `build-index-only --vcf` run stopped while it writes the dict-less tiers or
the indel genotype store no longer leaves a part of them that searches use as the whole"** (under
`### Fixed`, line 614 of `CHANGELOG.md` at aa9c54b) still says twice that a Settings build passes no
`--samplesID`. Replace (line numbers are those of `CHANGELOG.md` at aa9c54b):

1. lines 636-637:

   > `build-index-only` again to its end (for the tiers, from the command line, with their
   > `--samplesID`, which a build from the web Settings page passes none of) or, for the tiers and the

   with

   > `build-index-only` again to its end (for the tiers, with their `--samplesID`, or built again in
   > the web Settings page) or, for the tiers and the

2. lines 647-648:

   > with `--samplesID` was stopped writing, which a run without it (a build from the web Settings
   > page) does not write: it exited 0, which the Settings page showed as done, while every search of

   with

   > with `--samplesID` was stopped writing, which a run without it (as a build from the web
   > Settings page was before this release) does not write: it exited 0, which the Settings page
   > showed as done, while every search of

The rest of that entry stays as it is. The new entries below go where their headings say: the first
two under `### Changed`, the others under `### Fixed`. The first entry under Changed gives the
tables as they are with the fix of the tier path's REF columns (the tierref track,
`PostProcess/new_simple_analysis.py`), with which those columns are as the dictionaries give
them. If that fix does not land in 2.8.1, add to that entry, after "and an indel's is still the
VCF's;": "on a SNP row the reference alignment has no lower-case mismatches, so its REF seed and
non-seed counts leave the mismatches out (they count only the row's bulges);".

### Changed
- **A search of a VCF dataset built in the web Settings page reads its dict-less tiers, so its
  tables change.** A variant build in Settings now passes the dataset's sample list as `--samplesID`
  when its VCFs have samples (see Fixed), so it also writes the dict-less tiers
  (`Dictionaries/registry_<dataset>/` and `Dictionaries/genotypes_<dataset>/`). Searches of an index
  built (or built again) in Settings now read them, as searches of a command-line `--samplesID`
  build and of a downloaded index do, where they read the dataset's dictionaries before. On a phased
  dataset a search now runs per-sample by default, as on the default index (*Population-level* on
  the search form is the faster worst-possible screen), and a per-sample search names on each row
  only the samples one of whose haplotypes carries all of the row's variants (of two samples with
  the same two SNPs, one that carries both on one haplotype is named on the row of both SNPs, and
  one that carries one on each haplotype is not, where a population-level search names both).
  Measured on a phased 6-sample panel, with a population-level search of one index on its
  dictionaries and then on its tiers, the tables change as follows: an indel and a SNP that one
  haplotype carries are reported as one off-target, for the samples one of whose haplotypes carries
  both (the dictionaries gave the SNP alone, aligned with a bulge, or, for the fewest mismatches and
  bulges, the indel alone, for every sample that carries either; a sample that carries them on
  different haplotypes is no longer named at that site, and the per-sample summary and the report's
  threshold counts change with that row); a SNP's allele frequency is computed from the genotypes of
  the samples the list names, to six significant digits (`0.333333`), where the dictionaries gave
  the VCF's `AF` (`0.3333`), and an indel's is still the VCF's; `Other_motifs` counts the alignments
  the tiers give at a site: one more at the site of that indel and SNP (their row), and one fewer at
  a site where the dictionaries also listed the reference alignment as a variant one, for the
  carriers of a SNP that destroys its PAM; and rows come in another order. The result folder also
  holds the SNP + indel and SNP + SNP co-occurrence and phase-confirmation tables, and the shareable
  report gives the dataset's SNP and indel counts and its SNP + indel cis co-occurrences. A dataset
  built in Settings before this release has no tiers until its variant index is built again there,
  and its searches read its dictionaries until then, as before.
- **A genome uploaded in the web Settings page under the name of an installed genome replaces it,
  and is refused while indexes built on a genome of that name are installed.** The upload went into
  that genome's folder, beside the files already there (and so into the shared copy, when that
  folder is a link to one). It is now built in a hidden folder and replaces the installed genome
  only once whole; a link is replaced and the copy it points to left alone. It is refused, at its
  first chunk, while indexes or variant genomes built on a genome of that name are installed
  (`genome_library/<PAM>_<N>_<genome>[+<dataset>]`, `Genomes/<genome>+<dataset>`): every search
  and build of the new genome would go on using them, with the sequences and positions of the one
  they were built from. Remove them first (an index under *Remove installed data*, a
  `Genomes/<genome>+<dataset>` folder from your data folder), or upload the genome under another
  name; to add the very genome an installed index was built on, use *Download genome*, which this
  does not refuse.

### Fixed
- **A VCF dataset added in the web Settings page can be searched from the browser, per sample too.**
  A web variant search reads the dataset's sample list, `samplesIDs/<dataset>.samplesID.txt`, only
  after the search and post-analysis. Nothing in Settings wrote one for a dataset added there, and
  the Settings build passed no `--samplesID`, so a dataset added and built only in Settings ran its
  whole search and then failed (`ERROR: samples IDs processing failed on dataset ...`), was
  labelled `[sites-only (allele frequencies)]` on the search form, and could not run per-sample.
  Adding a dataset in Settings (an upload, or a registered server folder), or building its variant
  index there, now writes that list from the VCF headers when there is none, in the 4-column layout
  of the shipped lists: exactly the samples the VCFs genotype (a command-line search refuses a list
  that names a sample no VCF genotypes, #46, or that leaves out one they genotype), each with the
  population, superpopulation and sex an installed sample list gives it (a downloaded panel's, or
  yours; of the lists that give it a superpopulation, the last by file name), and with the dataset
  as its population and `unknown` superpopulation and sex when no list gives it a superpopulation
  (the shareable report reads every list for a sample's superpopulation, and a list that gave
  `unknown` to a panel's samples turned them unknown in the reports of the panel's searches).
  The list is written again when a file added to the dataset, or a build of its index, gives another
  one (a file with other samples, a sample list installed or removed since), and it goes when the
  dataset is removed in Settings. A list from anywhere else (one you wrote or edited, or a
  download's) is never changed: the page tells its own by the SHA-256 it keeps of it in the hidden
  `samplesIDs/.<dataset>.samplesID.sha256`. The Settings build passes the list, as a hidden one-line
  listing, `samplesIDs/.<dataset>.samplesID.list`, so it also builds the dict-less tiers (which
  searches of its index then read: see Changed): the dataset reads `[genotyped, phased]`
  (`[genotyped, unphased]` for an unphased VCF) and per-sample analysis is available. It passes a
  list only when it names samples the VCFs genotype: a sites-only dataset (no sample columns, as
  the `hg38_mega` VCFs) is built without one, as before, so the per-dataset allele-frequency
  registry put in place for it is left as it is. Removing the dataset in Settings removes the two
  hidden files too. Checked in the browser on a phased 6-sample dataset uploaded and built only in
  Settings, once per chromosome and once as a single multi-chromosome file named without the
  genome: the per-sample search runs to its end and names the same carriers, with the same allele
  frequencies, as the same panel built from the command line with `--samplesID`. Submit now refuses
  a variant search whose sample list is missing, with a message that says how to get one, rather
  than letting it fail at its end (for a sites-only dataset, which has no samples to list, it says
  so, and that a search reports such a VCF's variants only through a per-dataset allele-frequency
  registry, as the shipped `hg38_mega` index has); and a genotyped dataset with no tiers yet reads
  `[genotyped, samples not indexed]`, not sites-only.
  **If you added a dataset in Settings,** build its variant index again there; to give its
  samples' populations and sex (the tiers count a sample of unknown sex as diploid on chrX and
  chrY), put your own list at `samplesIDs/<dataset>.samplesID.txt` first. A list put there, or
  written there again, once the variant index is built reaches its tiers only after you remove
  `Dictionaries/registry_<dataset>` and `Dictionaries/genotypes_<dataset>` from your data folder
  and build the index again: until then a build leaves the tiers as they are, from the list they
  were built with, while a search's tables and report take the list that is there.
- **The Settings VCF upload records the reference genome you pick, and the search form offers the
  dataset for it.** The chunked uploader looked the genome up with selectors Dash 2's dropdown does
  not have and sent none, so `VCFs/.<dataset>.refgenome` was never written; and the search form
  offered only datasets whose name holds the genome's, so one named like the card's own example
  (`my_pig_vcf`) never appeared, even once its variant index was built. The uploader now sends the
  dropdown's value; an upload into a dataset with no recorded genome is refused until one is picked
  (the card already said it is required); and the search form offers a dataset for the genome its
  marker names, with the bulges its index supports. A dataset with no marker is offered by its name,
  as before, and now also for the genome of an enriched genome built from it
  (`Genomes/<genome>+<dataset>`).
- **A genome uploaded in Settings is laid out as every search needs it.** A FASTA of several
  sequences was stored as it was, so every search on it was refused with a command-line-only fix,
  and a UCSC `chromFa.tar.gz` kept its `chroms/` folder and an archive of `.fasta` files their
  names, so their index held no `.bin`; all were listed as installed. An upload is now laid out as
  every search needs it: every FASTA of an archive goes to the top of the genome folder as a `.fa`
  file (`chr1.fasta` as `chr1.fa`: an index is built from the `.fa` files only), a FASTA of several
  sequences is split into one file per sequence, and a FASTA of one sequence is named after it: a
  one-chromosome subset (`hg38_chr22.fa` holding `>chr22`) or a bacterial genome went in as
  `<genome>.fa`, and as a search takes each file for the sequence it is named after, every search of
  it found its off-targets and then dropped them all (`none survived post-analysis`), though its
  index build said done. (A genome downloaded by Direct URL or from UCSC as one FASTA of one
  sequence is named after it the same way, and `PostProcess/genome_layout.py --split` renames such a
  file too.) Hidden files of an archive are left out: the `._chr1.fa` that macOS `tar` adds beside
  `chr1.fa` (its file attributes) was installed too, and every search read it as a chromosome and
  failed, though the build said done. An archive with two FASTA files that would go in under one
  name (`chr1.fa` and `chr1.fasta`, or two files of one sequence) is refused. (An upload under the
  name of an installed genome now replaces it: see Changed.)
  **A genome uploaded before** this fix that holds a multi-sequence FASTA, a `chroms/` folder,
  `.fasta` files, `._` files or a FASTA not named after its one sequence: remove it in Settings with
  the indexes built on it (*Remove installed data*), and any `Genomes/<genome>+<dataset>` variant
  genome from your data folder, then upload it again. (An index built from a multi-sequence FASTA
  holds the sites of every sequence under the name of the first: kept, it gave a chr2 off-target as
  `chr1:37386` instead of `chr2:17381`.)
- **`build-index-only` fails on a genome folder it cannot index, and never puts in place an index
  with no `.bin`.** `index-genome` indexes the `.fa` files at the top of the genome folder only and
  exits 0 with an empty index when there are none (a `chroms/` folder, `.fasta` names, gzipped
  FASTAs); `build-index-only` put that index in place, and the Settings page showed the build as
  done, though every search refuses an index with no `.bin`. It now ends with an ERROR (exit status
  1) that says what to do: before it builds anything, for a folder with no `.fa` file at its top,
  with a FASTA of several sequences (which every search refuses too; the message says to remove the
  indexes built from the folder while it held that file) or with a FASTA of one sequence not named
  after it (`walkone.fa` holding `>chr1`, whose off-targets every search found and then dropped; the
  message gives the `mv`, or `genome_layout.py --split`, which renames them all); and, putting
  nothing in place, when the engine writes no `.bin` for the reference or the variant genome. A
  reference index with no `.bin` that an earlier run left is no longer taken for present (the run
  said `Reference index already present` and the Settings page showed the build done, though no
  search uses it): the run says to remove it. A genome of one FASTA per sequence, each named after
  its sequence, builds exactly as before.
- **A single multi-chromosome `.vcf.gz` without its `.tbi` builds a variant index.** The
  documented auto-split read the file by region through its index, so a file uploaded alone (the
  Settings uploader takes one file at a time) failed with `Could not load .tbi index ...
  Auto-splitting the VCF dataset into per-chromosome files failed`. A file with no `.tbi` or `.csi`
  is now indexed first, in a scratch folder (a registered server folder is never written to); a
  plain-gzip file, which tabix cannot index, is first recompressed there in bgzip format. What
  tabix and bcftools print goes to stdout (a warning of theirs on stderr fails a pipeline stage),
  and an unsorted file fails with their reason and the `bcftools sort` to run. The split folder now
  gets the same permissions as any other folder made there, not the `0700` of the temporary folder
  it came from. Per-chromosome datasets are not touched, as before.
- **No false "Samples degraded ... Re-download this index" alert after a Settings variant build.**
  The alert (new in this release) is meant for an index download whose genotype store did not
  install, but it matched any `genotype store` text, so the line every variant build prints
  (`[indel-snp] phased indel genotype store: ...`), and a publish's, raised it: wrong advice, as a
  local build has nothing to re-download. It now matches the download's own notes only, and is
  checked for `Download index` jobs only.
- **The Settings page says when a stopped build left its build directory behind.** A `docker stop`
  during a Settings build kills it before it can clean up, and each later container, with a new host
  name, cannot tell whether that build still runs, so it leaves the hidden directory (which can be
  as large as the finished index or enrichment), with a line on stdout that the page showed at most
  in the full log of a failed build. When a Settings job ends, done or failed, the page now names
  each such directory still there, with its size, and says to delete it once no other CRISPRme
  container is building. The Docker quickstart and `docs/TROUBLESHOOTING.md` give the command.
- **An annotation BED whose file name names no genome can be enabled; a second upload of the same
  bytes is added.** The upload said "Enable it under 'Manage annotations'", but that card offered
  only BEDs whose file name holds the genome's, so `walk_regions.bed` was offered for no genome. A
  BED whose file name names no installed genome is now offered for every genome, to enable for the
  one its coordinates are on; one that names an installed genome is still offered for that genome
  only, and the card says so. The upload widget is emptied after each file, so a next file with the
  same bytes (a copy under another name) is added, where nothing happened.
- **Settings: the "Started — progress is shown in the Activity panel" note is no longer shown as an
  error, and goes when its job ends;** the VCF card no longer lists "a URL" as a source it does not
  offer, and says that a single multi-chromosome file's `.tbi` is optional and where the dataset's
  sample list is; the web interface guide no longer says that custom VCFs and PAMs need the command
  line and a server restart, nor that an annotation upload takes a `.bed.gz`.
- **The shareable report's "By superpopulation" chart keeps a sample's superpopulation when a
  sample list read after it gives none.** The report reads every list of `samplesIDs/` and took a
  sample's superpopulation from the last list (by file name) that names it, so a list that gives
  it as `unknown` or `NA` (one you wrote for your own dataset of some of a panel's samples, or one
  the Settings page wrote before a panel naming them was installed) turned those samples unknown in
  the reports of the panel's searches. A list that gives no superpopulation no longer replaces one
  another list gives; where lists give two, the last one's is kept, as before.
- **`Other_motifs` and the alternative alignments count and list each reference alignment once.**
  The merge step of a search, which keeps each site's best alignment and lists the others as
  alternative alignments, took a site's reference rows as they came, and two steps before it give it
  the same row more than once: each variant dataset's post-analysis carries every row of the one
  reference search, and on the dict-less tiers (searches of a `--samplesID` build, of a downloaded
  index such as the default one, and of a variant index built in the web Settings page) the
  post-analysis also gives the reference alignment of each variant window again. So a search that
  read the tiers counted each reference alignment over a SNP twice in `Other_motifs` (`444` where
  the dictionaries gave `222`), a search of two variant datasets counted every reference alignment
  once per dataset (`491` where one dataset gave `245`), and a per-sample search's file of
  alternative alignments listed each of them as many times, a copy of the site's best alignment
  among them (on a 6-sample panel, 2446 rows for 1888 alternative alignments on the tiers, and 3779
  rows for the panel's 1887 when a second dataset was searched with it). Each reference alignment
  is now kept once: `Other_motifs`, the `#Seq_in_cluster` columns of the `bestMerge` table, the
  alternative alignments and the results page's table change with it; the best alignment of each
  site, its scores, carriers and allele frequencies are as before, and a search of one dataset on
  its dictionaries, or of the reference alone, is unchanged. In the table of a site's haplotypes
  (`observed_haplotypes.tsv`, and the report's), a haplotype's CRISPR-Bulge score is taken from the
  rows of the alternative alignments, which pair the three criteria's lists by position, so it can
  change with them: in the two-dataset search, 2 of its 5 rows did. Run a search again for these
  counts and lists.
