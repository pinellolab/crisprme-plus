<!-- [2.8.1] entries of the "wset" track (web Settings / data manager), to merge into
CHANGELOG.md. All go under "### Fixed". -->

### Fixed
- **A VCF dataset added in the web Settings page can be searched from the browser, per sample
  too.** A web variant search reads the dataset's sample list, `samplesIDs/<dataset>.samplesID.txt`,
  only after the search and post-analysis. Nothing in Settings wrote one and the Settings build
  passed no `--samplesID`, so a dataset added and built only in Settings ran its whole search and
  then failed (`ERROR: samples IDs processing failed`), was labelled `[sites-only (allele
  frequencies)]` on the search form, and could not run per-sample. Adding a dataset in Settings (an
  upload, or a registered server folder) now writes that list from the VCF headers when there is
  none, in the 4-column layout of the shipped lists: exactly the samples the VCFs genotype (a listed
  sample no VCF genotypes would count as hom-ref and deflate every allele frequency), with the
  dataset as their population and `unknown` superpopulation and sex. A list from anywhere else (one
  you wrote, or a download's) is never changed. The Settings build passes the list, so it also
  builds the dict-less tiers: the dataset reads `[genotyped, phased]` and per-sample analysis is
  available. Checked in the browser on a phased 6-sample dataset uploaded and built only in
  Settings, once per chromosome and once as a single multi-chromosome file named without the
  genome: the per-sample search runs to its end and names the same carriers, with the same allele
  frequencies, as the same panel built from the command line with `--samplesID`. Submit now
  refuses a variant search whose sample list is missing, with a message that says how to get one,
  rather than letting it fail at its end; and a genotyped dataset with no tiers yet reads
  `[genotyped, samples not indexed]`, not sites-only. **If you added a dataset in Settings,** build
  its variant index again there; to give its samples' populations and sex, put your own list at
  `samplesIDs/<dataset>.samplesID.txt` first.
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
  and a UCSC `chromFa.tar.gz` kept its `chroms/` folder, so its index held no `.bin`; both were
  listed as installed. An upload is now laid out as a Direct URL download is: every FASTA of an
  archive goes to the top of the genome folder, and a FASTA of several sequences is split into one
  file per sequence. It is built in a hidden folder and replaces an installed genome of the same name
  only once whole (a genome folder that is a link to a shared copy is replaced and the copy left
  alone, where the upload used to write into it). **A genome uploaded before** this fix that holds a
  multi-sequence FASTA or a `chroms/` folder: remove it in Settings and upload it again.
- **`build-index-only` fails on a genome folder it cannot index, and never puts in place an index
  with no `.bin`.** `index-genome` indexes the `.fa` files at the top of the genome folder only and
  exits 0 with an empty index when there are none (a `chroms/` folder, `.fasta` names, gzipped
  FASTAs); `build-index-only` put that index in place, and the Settings page showed the build as
  done, though every search refuses an index with no `.bin`. It now ends with an ERROR (exit status
  1) that says what to do, before it builds anything, for a folder with no `.fa` file at its top or
  with a FASTA of several sequences (which every search refuses too), and when the engine writes no
  `.bin` for the reference or the variant genome. A genome of one FASTA per sequence builds exactly
  as before.
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
