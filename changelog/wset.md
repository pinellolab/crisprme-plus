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
  and its searches read its dictionaries until then, as before. Writing the tiers takes memory that
  a build without them did not need: at panel scale ~150-180 GB a chromosome, by default up to 4
  chromosomes at once (~700 GB, `CRISPRME_TIER_WORKERS` sets how many). A build that runs out of it
  fails with a message that says what to do (see Fixed).
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
  were built with. A search takes its samples' populations and superpopulations from the list that
  is there, and a SNP's allele frequency from the tiers, so on chrX and chrY that frequency counts
  each sample by the sex the tiers were built with.
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
  file too.) Each header is written as `>id` (see the entry on those). Hidden files of an archive
  are left out: the `._chr1.fa` that macOS `tar` adds beside `chr1.fa` (its file attributes) was
  installed too, and every search read it as a chromosome and failed, though the build said done. An
  archive with two FASTA files that would go in under one name (`chr1.fa` and `chr1.fasta`, or two
  files of one sequence) is refused. (An upload under the name of an installed genome now replaces
  it: see Changed.)
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
  indexes built from the folder while it held that file), with a FASTA of one sequence not named
  after it (`walkone.fa` holding `>chr1`, whose off-targets every search found and then dropped; the
  message gives the `mv`, or `genome_layout.py --split`, which renames them all) or with a header of
  more than `>id` (see the entry on those); and, putting nothing in place, when the engine writes no
  `.bin` for the reference or the variant genome. A reference index with no `.bin` that an earlier
  run left is no longer taken for present (the run said `Reference index already present` and the
  Settings page showed the build done, though no search uses it): the run says to remove it. A
  genome of one FASTA per sequence, each named after its sequence and with a header of `>id` (see
  the entry on those), builds exactly as before.
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
  the same alignment more than once: each variant dataset's post-analysis carries every row of the
  one reference search, and on the dict-less tiers (searches of a `--samplesID` build, of a
  downloaded index such as the default one, and of a variant index built in the web Settings page)
  the post-analysis also gives the reference alignment of a variant window on the plus strand again.
  So a search that read the tiers counted each reference alignment over a SNP twice in
  `Other_motifs` (`444` where the dictionaries gave `222`), a search of two variant datasets counted
  every reference alignment once per dataset (`491` where one dataset gave `245`), and a per-sample
  search's file of alternative alignments listed each of them as many times, a copy of the site's
  best alignment among them. On chr22 of the default index, a per-sample search gave 6152 rows for
  its 4941 alternative alignments, and 280 of its 1940 sites a higher `Other_motifs` (280 of 2356 at
  population level); on a 6-sample panel searched with a second dataset, 3779 rows for 1887. Each
  reference alignment is now kept once. Rows equal in every column but `#Seq_in_cluster` and the two
  score columns are one alignment: each dataset scores the rows itself, and CRISPR-Bulge can give
  one alignment two values (`4.3390600000000004e-23` and `4.33904e-23`), of which the higher stays.
  What changes is `Other_motifs` in the integrated table and on the results page, the
  `#Seq_in_cluster` columns of the `bestMerge` table, the file of alternative alignments and the
  `observed_haplotypes_raw` file made from it. The best alignment of each site, its scores, carriers
  and allele frequencies are as before, `observed_haplotypes.tsv` and the report do not change, and
  a search of one dataset on its dictionaries, or of the reference alone, is unchanged (a search
  with a sites-only dataset changes further: see the entry on those). Run a search again for these
  counts and lists.
- **A search with a sites-only dataset lists each site once and loses none; a VCF with no `AF` keeps
  a variant's best alignment.** The merge step of a search sorts each criterion's rows by guide,
  chromosome, position and score before it makes one row of each site, and its sort counted columns
  by blanks, so an empty column moved the later ones up. A row of a sites-only dataset (`hg38_mega`,
  alone or beside another dataset, since 2.5.2) has an empty `Samples`, so it went by its rsID and
  not its guide: such rows came before or after all the others, and the merge, which starts a new
  site only when the position rises, took unrelated sites for one. A site was listed more than once;
  a variant row was marked `Not_found_in_REF` (a `VARIANT` off-target in the report) though its site
  has a reference row; a site could lose its row altogether; and in a per-sample search a row of
  alternative alignments could pair the three criteria's alignments across sites, so
  `observed_haplotypes.tsv` gave a site another site's variants. The sorts now read the
  tab-separated columns, and a variant both datasets give is one row, with the genotyped dataset's
  carriers and both allele frequencies. On chr22 (guide `CTAACAGTTGCTTTTATCACNNN`, 4 mismatches, 1+1
  bulges), `hg38_mega` alone goes from 1621 rows to 1384, rows within 3 bp of another from 556 to 0,
  `Not_found_in_REF` rows from 659 to 422, and it lists 50 sites it did not (the best at
  chr22:50091459, CFD 0.022); the report's `REFERENCE` and `VARIANT` off-targets go from 766 and 659
  to 962 and 422, and its CFD >= 0.05 net from 53 sites to 46. The default index with `hg38_mega`
  goes at population level from 3015 rows to 2520 (within 3 bp 997 to 4, `Not_found_in_REF` 2053 to
  1558, 21 sites more); per sample, with the next entry's fix, from 2406 rows to 2137 (671 to 2,
  1444 to 1175, 73 sites more), and its rows of alternative alignments with the fewest or the
  CRISPR-Bulge block more than 100 bp from the CFD block go from 28 and 27 to 0, and its haplotypes
  with a variant more than 1 kb from their site from 26 to 0. A row with only its `AF` empty (a VCF
  with no `AF`, which the input check warns about: every variant on the dictionaries, every indel on
  the tiers) stayed at its site but went by its REF score, so the merge could keep a lower-scoring
  alignment of a variant: on a 6-sample panel with no `AF`, the perfect match a deletion makes at
  chr2:10000 read CFD 0.9 with a bulge and now reads 1.0 with none, and the report's perfect matches
  go from 3 to 4. A search whose rows have no empty column (the default index alone, any VCF with
  `AF`, the reference alone) sorts as before. **If you searched a sites-only dataset or a VCF with
  no `AF`,** run the search again.
- **A sites-only dataset goes the population-level way in a per-sample search too, as the command
  line said it did, and loses no multi-variant worst case there.** A sites-only dataset
  (`hg38_mega`) has no genotypes for a per-sample search to resolve, and `complete-search` said that
  `the population-level analysis runs instead`. It did not: the dataset's rows went through the
  per-sample path, which joins two variants only through carriers they share, and it has none, so a
  site that takes two variants together lost its row, and others scored lower, with no warning.
  Installing the default index and `hg38_mega` lists both for a search
  (`docs/PRECOMPUTED_INDEXES.md`), and the two then ran that way with no mode flag. Now each dataset
  goes by its registry (its `data_type`, or for a registry that has none its aggregation): a
  sites-only dataset's post-analysis is the population-level one whatever the search's mode.
  `--per-sample` on sites-only datasets alone is the population-level analysis (the same files as
  with no flag), with a WARNING on stdout. Beside a genotyped dataset, whose rows stay per-sample,
  `complete-search` says which goes which way. With no mode flag a phased dataset beside a
  sites-only one now runs per-sample whichever name sorts first; with the sites-only name first it
  ran population-level. On chr22 (guide `CTAACAGTTGCTTTTATCACNNN`, 4 mismatches, 1+1 bulges),
  `hg38_mega --per-sample` goes from 1318 rows to 1384, the table of `hg38_mega` with no flag, byte
  for byte, and lists 63 sites it did not; the default index with `hg38_mega`, per-sample as it runs
  with no flag, goes from 2077 rows to 2137 (`Not_found_in_REF` 1115 to 1175, alternative alignments
  5863 to 6226, `observed_haplotypes.tsv` 1556 rows to 1796). It lists 57 sites it did not and loses
  none, 14 sites score higher and none lower (chr22:34270481 from CFD 0.173 to 0.233, and
  chr22:22881199 reads 0.274 where chr22:22881195 read 0.262), and the report's CFD >= 0.05 and >=
  0.2 tables each gain one site. Such rows have no carriers and read `putative` in the report's
  Observed column (`observed` for one variant). A search with no sites-only dataset, `hg38_mega`
  with no flag and any search at `--population-level` are unchanged, byte for byte. **If you ran a
  per-sample search with a sites-only dataset,** run it again. Covered by
  `PostProcess/test_sites_only_per_sample.py`.
- **A command-line search of two datasets is no longer refused because only one has its VCFs on
  disk.** The input check of `complete-search` compared the samples of all the `--samplesID` lists
  with the VCFs on disk (#46). With `hg38_walk`'s VCFs there and `hg38_walkB` installed as an index,
  it stopped with `--samplesID lists 3 sample(s) genotyped in no VCF (B4, B5, B6)`, and told you to
  filter a list that was right; so too beside a sites-only dataset, on disk or not, whose VCFs have
  no sample columns and whose list names `1000G2021` and `HGDP`. Each step that leads there is
  documented (`download --what vcf` for one dataset, several indexes in one list, your own cohort
  beside a dataset). Each list is now compared on its own: up to as many lists as the search names
  such datasets, those that name no sample of the VCFs on disk, are left out. A list that names some
  of them and others they lack still stops the search, with those others named; so does a list of no
  dataset when every dataset's VCFs are there. Measured on a 6-sample panel, the search that was
  refused gives the same table as with both datasets' VCFs on disk, and as with neither. A web
  search does not run this check and is not affected. Covered by
  `PostProcess/test_samplesid_overlisting.py`.
- **A variant build that runs out of memory as it writes its tiers ends with a message, and no
  longer hangs.** `build-index-only --samplesID`, and in this release every variant build in the web
  Settings page, writes the dict-less tiers and the phased indel genotype store in forked workers,
  one chromosome each. The pool they ran in did not see a worker die: when the OOM killer took one,
  the run waited for good, at 0% CPU, with nothing on stdout or stderr. In Settings the build read
  `running` for ever, with its buttons locked and every later job queued behind it, until the
  container was restarted; the same build finished in 2.8.0, which wrote no tiers. A worker that
  dies now sends the chromosomes left to one worker at a time, with a WARNING on stdout, and the
  build goes on if they fit. If a worker dies then too, the run ends with exit status 1 and an ERROR
  that says memory is the likely cause, names `CRISPRME_TIER_WORKERS` (`CRISPRME_INDELGT_WORKERS`
  for the store), and says how to go without what it could not write; Settings shows the build as
  failed with that message. A run with one worker, or of a dataset of one chromosome, which did not
  hang but died at exit status 137 with no message, now ends with that ERROR too. What it could not
  write stays marked unfinished, so a search refuses it, and the next run writes the rest; no
  variant index is built by that run. To build with no tiers, remove them as the message says and
  run `build-index-only` from the command line with no `--samplesID`: the Settings page has no way
  to. No worker outlives a run that is killed outright either, one that waits for its next job looks
  for the run every few seconds. Measured on a synthetic panel of 3000 samples and a SNP every 3 bp
  on 3 x 50 kb: with a limit of 1.5 GB the Settings build reads `failed` with the message after 5 s,
  and with 6 GB, where three workers at once are too many, the build finishes one chromosome at a
  time. A build with memory enough writes the same files as before, byte for byte, and so does that
  one.
- **A genome whose FASTA headers carry a description, or whose lines end in CRLF, no longer indexes
  as done and then loses every search's results.** CRISPRitz names a chromosome after the whole
  header line, and the pipeline splits on blanks, so the FASTAs NCBI, GENCODE, Ensembl and GATK ship
  went wrong, each in its own way, with the same sequence that gives 7 sites under a plain `>chr1`.
  An NCBI `>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome` or a GATK
  header: exit status 0, `No off-targets found` (in the browser, `Finished: no off-targets were
  found`). GENCODE `>chr1 1`: 7 sites, but with `chr1 1` for a chromosome and `NA` for every
  annotation. Ensembl `>1 dna:chromosome ...`: the search fails at integration (`Start coordinate is
  too large ... BEDOPS`). CRLF: 451 KB of warnings on stderr at exit status 0, and `No off-targets
  found`. A variant build on `>chr1 1` failed with `Indels indexing failed`. An upload in Settings,
  a Direct URL or UCSC download and `genome_layout.py --split` now write each header as `>id`, and
  every line with no CR, also for a FASTA already named after its sequence; a file that needs none
  of this is not touched. Each of those genomes then gives the sites and table of the same sequence
  under a plain header, through the browser and the command line, and a variant build on one gives
  the same table too. A folder left with such a file is refused before anything is built or
  searched: `build-index-only`, `complete-search` and a web search end with an ERROR that gives the
  `--split` command and says to remove the indexes built from the folder as it was (`--split` lists
  them). Such an index is refused as well, where its `.bin` names show it: it holds the old names
  once the FASTA is written as it should be, and every search of it went on losing results, and a
  build took it as present. An index built from a FASTA with a clean header line and CRLF after it
  has `.bin` names that do not show it, so it is not caught: remove it as the ERROR says. A sequence
  id that holds `|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`, `>`, a backquote or `#` is refused too,
  at an upload, a download, `--split`, a build and a search (`gi|49175990|ref|NC_000913.2|`, a PanSN
  `HG01255#1#CM086702.1`; measured character by character through the real command line, `.`, `:`,
  `_`, `-`, `=`, `,`, `+`, `@`, `%`, `^`, `~`, `!`, `*`, `?` and brackets work): a search of one
  failed in a shell command (`/bin/sh: Syntax error`), or, for `#`, with `none survived
  post-analysis`. A file whose lines end in CR alone is refused by `--split`, which says to convert
  it, as its first line holds the whole file. **If you have such a genome,** run `python
  PostProcess/genome_layout.py --split <genome folder>`, remove the indexes it names, and build
  again, or remove it in Settings and upload it again. A search of a genome of plain headers is
  unchanged. Covered by `PostProcess/test_genome_layout.py` and
  `PostProcess/test_build_index_only_layout.py`.
