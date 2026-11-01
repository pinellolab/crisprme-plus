<!-- [2.8.1] entries of the "wset" track (web Settings / data manager), to merge into
CHANGELOG.md. All go under "### Fixed".
For the docs track: a Settings variant build now passes the dataset's sample list as
--samplesID (first entry below), so the build-index-only / search ERROR for unfinished
tiers no longer says "a build from the web Settings page passes none": it says to run
build-index-only again to its end "with <dataset> and its --samplesID (or build it again in
the web Settings page, which passes the dataset's sample list)". The two places in
CHANGELOG.md's existing [2.8.1] entry on unfinished parts that say a Settings build passes
no --samplesID ("from the command line, with their `--samplesID`, which a build from the
web Settings page passes none of" and "a run without it (a build from the web Settings
page) does not write") should say the same: "(for the tiers, with their `--samplesID`, or
built again in the web Settings page)" and "a run without it (as a build from the web
Settings page was before this release)".
The first entry's last sentence describes the tables of a search that reads the tiers as
they are at this commit: if the fix of their REF columns (round 3 review F3:
PostProcess/new_simple_analysis.py writes the REF alignment of an observed or fast-mode row
without lowercase mismatches, so its REF seed and non-seed counts read 0) lands in 2.8.1,
drop its clause "and on a SNP row the reference alignment has no lowercase mismatches and
its REF seed and non-seed counts read 0, as on every search that reads the tiers". -->

### Fixed
- **A VCF dataset added in the web Settings page can be searched from the browser, per sample too.**
  A web variant search reads the dataset's sample list, `samplesIDs/<dataset>.samplesID.txt`, only
  after the search and post-analysis. Nothing in Settings wrote one and the Settings build passed no
  `--samplesID`, so a dataset added and built only in Settings ran its whole search and then failed
  (`ERROR: samples IDs processing failed`), was labelled `[sites-only (allele frequencies)]` on the
  search form, and could not run per-sample. Adding a dataset in Settings (an upload, or a
  registered server folder) now writes that list from the VCF headers when there is none, in the
  4-column layout of the shipped lists: exactly the samples the VCFs genotype (a listed sample no
  VCF genotypes would count as hom-ref and deflate every allele frequency), each with the
  population, superpopulation and sex an installed sample list gives it (a downloaded panel's, or
  yours), and with the dataset as its population and `unknown` superpopulation and sex when no list
  names it (the shareable report reads every list for a sample's superpopulation, and a list that
  gave `unknown` to a panel's samples turned them unknown in the reports of the panel's searches).
  The list is written again when a file added to the dataset, or a build of its index, gives another
  one (a file with other samples, a sample list installed or removed since), and it goes when the
  dataset is removed in Settings. A list from anywhere else (one you wrote or edited, or a
  download's) is never changed. The Settings build passes the list, so it also builds the dict-less
  tiers: the dataset reads `[genotyped, phased]` and per-sample analysis is available. It passes a
  list only when it names samples the VCFs genotype: a sites-only dataset (no sample columns, as the
  `hg38_mega` VCFs) is built without one, as before, so the per-dataset allele-frequency registry
  put in place for it is left as it is. Checked in the browser on a phased 6-sample dataset uploaded
  and built only in Settings, once per chromosome and once as a single multi-chromosome file named
  without the genome: the per-sample search runs to its end and names the same carriers, with the
  same allele frequencies, as the same panel built from the command line with `--samplesID`. Submit
  now refuses a variant search whose sample list is missing, with a message that says how to get
  one, rather than letting it fail at its end (for a sites-only dataset, which has no samples to
  list, it says so, and that a search reports such a VCF's variants only through a per-dataset
  allele-frequency registry, as the shipped `hg38_mega` index has); and a genotyped dataset with no
  tiers yet reads `[genotyped, samples not indexed]`, not sites-only. Searches of an index built (or
  built again) in Settings now read its dict-less tiers, as searches of a command-line `--samplesID`
  build and of a downloaded index do, where they read the dataset's dictionaries before, so their
  tables change: an indel and a SNP that one haplotype carries are reported as one off-target, for
  the samples that carry both (the dictionary path reported the SNP alone, aligned with a bulge, for
  every carrier of the SNP), allele frequencies have six decimals, rows come in another order, and
  on a SNP row the reference alignment has no lowercase mismatches and its REF seed and non-seed
  counts read 0, as on every search that reads the tiers.
  **If you added a dataset in Settings,** build its variant index again there; to give its
  samples' populations and sex, put your own list at `samplesIDs/<dataset>.samplesID.txt` first.
- **The Settings VCF upload records the reference genome you pick, and the search form offers the
  dataset for it.** The chunked uploader looked the genome up with selectors Dash 2's dropdown does
  not have and sent none, so `VCFs/.<dataset>.refgenome` was never written; and the search form
  offered only datasets whose name holds the genome's, so one named like the card's own example
  (`my_pig_vcf`) never appeared, even once its variant index was built. The uploader now sends the
  dropdown's value; an upload into a dataset with no recorded genome is refused until one is picked
  (the card already said it is required); and the search form offers a dataset for the genome its
  marker names, with the bulges its index supports. A dataset with no marker is still offered by its
  name, or by an enriched genome built from it.
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
  name (`chr1.fa` and `chr1.fasta`, or two files of one sequence) is refused. It is built in a
  hidden folder and replaces an installed genome of the same name only once whole (a genome folder
  that is a link to a shared copy is replaced and the copy left alone, where the upload used to
  write into it). An upload is refused, at its first chunk, while indexes or variant genomes built
  on a genome of that name are installed: every search and build of the new genome would go on using
  them, with the sequences and positions of the one they were built from.
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
  1) that says what to do, before it builds anything, for a folder with no `.fa` file at its top,
  with a FASTA of several sequences (which every search refuses too; the message says to remove the
  indexes built from the folder while it held that file) or with a FASTA of one sequence not named
  after it (`walkone.fa` holding `>chr1`, whose off-targets every search found and then dropped; the
  message gives the `mv`, or `genome_layout.py --split`, which renames them all), and when the
  engine writes no `.bin` for the reference or the variant genome. A reference index with no `.bin`
  that an earlier run left is no longer taken for present (the run said
  `Reference index already present` and the Settings page showed the build done, though no search
  uses it): the run says to remove it. A genome of one FASTA per sequence builds exactly as before.
- **A single multi-chromosome `.vcf.gz` without its `.tbi` builds a variant index.** The
  documented auto-split read the file by region through its index, so a file uploaded alone (the
  Settings uploader takes one file at a time) failed with `Could not load .tbi index ...
  Auto-splitting the VCF dataset into per-chromosome files failed`. A file with no `.tbi` or `.csi`
  is now indexed first, in a scratch folder (a registered server folder is never written to); a
  plain-gzip file, which tabix cannot index, is recompressed with bgzip there first. What tabix and
  bcftools print goes to stdout (a warning of theirs on stderr fails a pipeline stage), and an
  unsorted file fails with their reason and the `bcftools sort` to run. The split folder now gets the
  same permissions as any other folder made there, not the `0700` of the temporary folder it came
  from. Per-chromosome datasets are not touched, as before.
- **No false "Samples degraded ... Re-download this index" alert after a Settings variant build.**
  The alert (new in this release) is meant for an index download whose genotype store did not
  install, but it matched any `genotype store` text, so the line every variant build prints
  (`[indel-snp] phased indel genotype store: ...`), and a publish's, raised it: wrong advice, as a
  local build has nothing to re-download. It now matches the download's own notes only, and is
  checked for `Download index` jobs only.
- **The Settings page says when a stopped build left its build directory behind.** A `docker stop`
  during a Settings build kills it before it can clean up, and each later container, with a new host
  name, cannot tell whether that build still runs, so it leaves the hidden directory (as large as
  the finished index or enrichment) with a line on stdout the page did not show. When a Settings job
  ends, done or failed, the page now names each such directory still there, with its size, and says
  to delete it once no other CRISPRme container is building. The Docker quickstart and
  `docs/TROUBLESHOOTING.md` give the command.
- **An annotation BED whose file name names no genome can be enabled; a second upload of the same
  bytes is added.** The upload said "Enable it under 'Manage annotations'", but that card offered
  only BEDs whose file name holds the genome's, so `walk_regions.bed` was offered for no genome. A
  BED whose file name names no installed genome is now offered for every genome, to enable for the
  one its coordinates are on; one that names a genome is still offered for that genome only, and the
  card says so. The upload widget is emptied after each file, so a next file with the same bytes (a
  copy under another name) is added, where nothing happened.
- **Settings: the "Started — progress is shown in the Activity panel" note is no longer shown as an
  error, and goes when its job ends;** the VCF card no longer lists "a URL" as a source it does not
  offer, and says that a single multi-chromosome file's `.tbi` is optional and where the dataset's
  sample list is; the web interface guide no longer says that custom VCFs and PAMs need the command
  line and a server restart.
- **The shareable report's "By superpopulation" chart keeps a sample's superpopulation when a
  sample list read after it gives none.** The report reads every list of `samplesIDs/` and took a
  sample's superpopulation from the last list (by file name) that names it, so a list that gives
  it as `unknown` or `NA` (one you wrote for your own dataset of some of a panel's samples, or one
  the Settings page wrote before a panel naming them was installed) turned those samples unknown in
  the reports of the panel's searches. A list that gives no superpopulation no longer replaces one
  another list gives; where lists give two, the last one's is kept, as before.
