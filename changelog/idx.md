<!-- idx track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. The [2.8.1] entry
that says "An INDELS index that `complete-search` builds from a genome it enriched itself (the
on-demand build) still has no SNP codes, as before" stops being true with these fixes: drop that
sentence. The entry below on the SNP codes needs the indel track's fix for a dataset with no
registry. -->

### Fixed

- **A 0-bulge variant search on a downloaded index no longer loses every SNP-path variant row.**
  With `--bDNA 0 --bRNA 0` the variant search scanned `Genomes/<genome>+<dataset>/`, which a
  download does not ship. The engine printed `ERROR! ... is not a directory.` on stdout and exited
  0, so the run went on with the reference and indel rows alone, wrote a full report, and exited 0
  with an empty error log. **Affected:** every such search of a downloaded index without the
  enriched genome, on the CLI and in the web form's Advanced panel, in 2.8.0 and back to at least
  2.3.3 (searches with bulges, and installs with the enriched genome, were right). **What to do:**
  run those searches again. Measured on hg38 with the default index (TRAC guide, 4 mismatches, NRG,
  per-sample): 483 rows before, 722 after, which is byte-identical to a scan of the enriched genome.
  The 240 new sites are all variant sites, and one site listed before is now reported 1 nt over, at
  a variant alignment that scores higher; sites found only in the variant genome go from 80 to 319;
  sites at CFD >= 0.5 from 7 to 11, >= 0.2 from 33 to 50, >= 0.05 from 156 to 233; and 207 of the
  sites already listed change too (a variant alignment that scores higher, or another motif counted).
  The search now takes the installed variant index when any chromosome's enriched FASTA is missing,
  with the same `-var`, `--max-edits` and `--prune-pam` a search with bulges passes; with the whole
  enriched genome in place it scans it, as before.
- **A search the engine never ran now fails the run.** `crispritz.py` says what stops it before it
  searches on stdout and exits 0, so a search that never ran passed for one that found nothing.
  After each reference and variant search, the run now stops with an ERROR when the targets file is
  missing. A search that finds nothing writes a file with its header alone, and still ends in an
  empty result with exit 0.
- **A search with bulges no longer ends empty when its library cannot be written to and lacks the
  index it needs.** On a read-only mount, without write permission or with no inode left, the failed
  lock read as "another run built it", so the search went on with a folder that is not there: the
  engine said so on stdout and exited 0, and the run ended `No off-targets found ... completed
  successfully` (web: "Finished: no off-targets were found"), or, with the reference index in place
  and the variant index short, reported the reference and indel rows alone. **Affected:** every
  such search, in 2.8.0 too. It now stops with one ERROR that names the index and the library and
  says what to do, shown on the web job's banner. A library that holds the right index still serves
  the search, byte-identically, and a 0-bulge search, which opens no index there, still runs.
- **A guide that does not fit an index's stored window is no longer searched on it.** An index name
  carries no guide length, and the engine searches any guide: past what the index stores, it
  reports wrong alignments or loses sites at exit 0, or, on a site planted at a contig start, stops
  with `Non-numeric start coordinate`. A 20-nt search at default flags builds `<PAM>_2`, and the next
  23-nt search took it: on real chr22, 181 rows against 1003, the on-target no longer a perfect
  match. On the web, the same, with no message. `build-index-only` with a 23-nt PAM file over a 20-nt
  index said "already present" and rewrote `.pam_build` with the wrong geometry. **Affected:** a
  search of a guide longer than its index was built for, in 2.8.0 too; every published index is
  20-nt, so a 20-nt search is not (nothing changes for it), a 21-nt one at up to 1 bulge of each
  type on an N=3 index is not, and any longer one is. **What to do:** run those searches again;
  they build their own index, or say how to clear the way. The length comes from the index's own
  `.bin` header, so no index needs rebuilding. A search takes an index only where built + N >= guide +
  bDNA + (1 if bRNA > 0), measured through the real CLI on 246 cells: every cell inside is
  byte-identical to a dedicated index, every cell outside differs. When an index of the same name
  is in the way it stops with one ERROR naming the index, both lengths and what to do; the index is
  never removed. An index built for this very guide is now preferred to one that merely fits, so a
  20-nt search beside a 23-nt `NGG_2` and a 20-nt `NGG_3` takes `NGG_3`. `build-index-only` stops
  the same way, and leaves the record of an index built for another guide. The web form's bulge
  caps and PAM list follow the same rule, and two PAM files that read the same show their guide
  length ("Cas12a · TTTV (23 nt)").
- **A guide longer than its PAM file's spacer is refused.** The search kept the spacer's length of
  it and dropped the rest with no message: a 23-nt guide with `20bp-NGG` was searched without its
  last 3 nt (2 rows against 35, and none at 0 bulges, at exit 0), under the full guide's name. The
  pre-flight now says which bases would go and what to do. Padding with extra `N`s drops nothing
  and stays accepted. CLI `--guide` only: the web form and `--sequence` never pass such a guide.
- **A search that enriches a dataset itself no longer loses targets that need an indel and a SNP
  beside it, and no longer keeps a blind INDELS index.** Only `build-index-only` put the SNP codes
  on the fake-indel genome. A search that found the dataset's VCFs and no index for its PAM or
  bulges (a download plus `--what vcf`, or a renamed genome folder), and the web's first 0-bulge
  search of a dataset that has only its VCFs, enriched the genome without them, so a cis SNP+indel
  site read 2 mismatches, CFD 0.096, listing a carrier that does not hold both, against 1 mismatch,
  CFD 0.353 on a built index; the INDELS index stayed in the library, so every later run gave the
  same, and `build-index-only` said "already present". Exit 0, no message. **Affected:** every such
  search since 2.5.0. Now the search script overlays the codes before the fake-indel genome goes
  into place, and again before it builds an INDELS index from one that is there, and writes the
  phased indel genotype store where it enriches. Measured on the walk genome: the on-demand search
  is byte-identical to a search of a built index, in all 41 result files. The first web search gets
  the same row, carrier and scores (it needs the indel track's fix, in this release, to resolve the
  codes without tiers); it has no registry, so its report says no co-occurrence table is written
  rather than "unaffected", and the CLI warning says the same. **What to do:** an INDELS index built
  by such a search stays blind until it is rebuilt. Run `build-index-only` with the dataset's VCFs and
  your PAM: it now codes a fake-indel genome it finds without codes, removes the INDELS indexes of
  every PAM and budget of that dataset built before the codes, says which, and builds its own (a
  later search rebuilds the others). A genotyped dataset with no `Dictionaries/indel_genotypes_<vcf>/`
  (an index built before 2.5.0) now gets a `WARNING` that says so and what to run.
- **A downloaded variant index asked for a PAM or bulge budget it does not hold now says so.**
  The search tried to enrich a dataset with no VCFs, and failed in `mv: cannot stat ...` four times
  and `Genome enrichment failed`, after a pre-flight that printed `[OK] ... precomputed variant index
  ... present`, leaving empty folders; the web form promised that a first search would prepare the
  genome, and its job ended "please remove it and try to submit again". This is likely the first
  slip on the default index (searching `NRG_3_hg38+hg38_HGDP1kGP` with NGG), and on a K=8 index
  without `--prune-pam`. **Affected:** a search with no index for its PAM or bulges and neither the
  dataset's VCFs nor its enriched genome installed, in 2.8.0 too. It now stops before it does
  anything, on the CLI before the pre-flight's `[OK]` and on the web job's banner, with one ERROR that
  names the genome, dataset, PAM, bulges and guide length asked for, the variant indexes installed,
  and what to do: search with the bulges an index of this PAM holds, `--prune-pam --index-path
  genome_library` where a wider all-N index serves it, or `download` or `build-index-only` another
  index. Web serving of a wider all-N index is not new here: the web still needs an all-N index of
  exactly the PAM's width. An install with the enriched genome (the build-inputs kit) or the VCFs
  still builds the index and is byte-identical to a dedicated one. The README no longer says a
  PAM-less index must be searched from the CLI: only a wider one must.
- **With no `--bDNA`/`--bRNA`, a search no longer becomes a 0-bulge search with no message.** The
  default of 1 bulge of each type is capped by the indexes a search can reach, but that cap looked
  only at the working directory's `genome_library/` and `Genomes/<name>`. A library at
  `--index-path` (variant, K=8 with `--prune-pam`), or a `--genome` folder elsewhere (the
  shared-mount recipe), gave a cap of 0, and the run was byte-identical to a 0-bulge search: 33 rows
  against 38 with bulges, with a site that needs a bulge reported at a worse alignment or not at all
  (2.8.0 too). **What to do:** run those searches again, or pass the bulges. The cap now reads
  `--index-path` and the `--genome` folder given, so these searches take their 1 bulge, and where no
  usable index holds bulges (K=8 without `--prune-pam`, an index of another PAM) a `WARNING` says
  the search runs with none, and what to pass.
- **`complete-test` works on any index you installed, and no longer breaks when the K=8 index sits
  beside the default.** Its smoke wrote the example guide and a 3' PAM from the index's name, so on
  a PAM wider than 3 nt it stopped with `guide ... is 23 nt but the PAM expects 28 ... Fix: append 8
  'N'`, about a file it rewrites each run; it took the first name in sorted order, so
  `NNNNNNNN_3_hg38` shadowed `NRG_3_hg38`, on the install `DOCKER_QUICKSTART` describes; and on a 5'
  index it searched the wrong side, which under `--chrom` read `PASSED ... holds none of its
  off-targets` with sites there. Each failed run also left a PAM file such as
  `PAMs/20bp-TTTV-SpCas9.txt`, which the web then offered as a nuclease. **Affected:** the
  documented self-test on a K=8, SaCas9, Cas12a or other non-SpCas9-class index, in 2.8.0 too; a
  leftover `PAMs/20bp-<motif>-SpCas9.txt` of such a motif can be removed. It now takes a 3' index of
  a 3-nt PAM first, and an all-N index wider than 3 nt last; searches a wider 3' all-N index for NRG
  with `--prune-pam`, byte-identical to a dedicated NRG index (measured on the synthetic genome,
  and on the published K=8 indexes on chr22); and searches any other for its own PAM, on the side
  its `.pam_build` records (3' without one; for a dedicated index, what the PAM files of its motif
  in `PAMs/` say), with the example guide padded for it. It writes nothing into `PAMs/`: it uses the
  file there that holds the line it needs, or a hidden one for the run. With another nuclease's PAM
  an empty result on the whole genome passes too, as it does on one chromosome. An index whose side
  cannot be told stops in one line. Results on NGG, NRG, NNN and a 2-nt PAM are unchanged.
- **A variant search on a genotyped index without `--samplesID` no longer ends in a `KeyError` after
  the whole search.** The help marks the flag optional and the run warned "using an empty sample
  set", a case meant for sites-only panels; on a genotyped index the rows name carriers, and
  `process_summaries.py` looked each up in an empty table at "Creating images", then deleted the
  results: 65 to 78 minutes on hg38 for a `KeyError`. An empty or a partial list on a folder with no
  VCFs did the same, and so did `complete-test` with a sample list under another name (2.8.0 too).
  **What to do:** run those searches again. Without the flag the search now uses
  `samplesIDs/<dataset>.samplesID.txt` for each dataset that has one, which an index download
  installs, and says so: every result file but its own list `.list_samplesID.txt` is byte-identical
  to passing that list, in both analysis modes. With no such file the search still ends with the
  same `integrated_results`, with empty by-sample and population tables, and the warning says so.
  A carrier a partial list does not name has no row in those tables.
- **`build-index-only` records the folder an index of a genome outside `Genomes/` holds, and a
  Direct URL download no longer replaces a genome something is built on.** A search takes a
  genome's index by the folder's name alone. `complete-search` records, in an index it builds of a
  folder outside `Genomes/`, which folder it holds (`.genome_source`); `build-index-only` did not.
  So after `build-index-only --genome projA/ref`, a search of another folder named `ref` reported
  the results of `projA/ref`, byte for byte, at exit 0; `build-index-only` of a folder named like
  one in `Genomes/` built that one's index from it, which every search of the `Genomes/` folder
  then took; a recorded folder whose FASTA files had changed was "already present"; and the search
  of a personal-assembly haplotype by its own path (`Assemblies/<individual>/<haplotype>/genome`)
  was refused on the index `build-index-only` had just built for it, which it had said would be
  reused. `download --what genome --source url`, the Settings page's Direct URL, put another
  genome's FASTA files under an installed name and kept its indexes, which every search then used,
  with the old genome's sequences and positions. **Affected:** `build-index-only` of a folder
  outside `Genomes/` and a Direct URL download over a genome with an index, in 2.8.0 too; the
  refused haplotype search is new in 2.8.1. **What to do:** a search run on such an index, or on a
  genome downloaded over another, searched the other folder's sequence: remove that index and run
  it again. `build-index-only` now records the folder, refuses one named like a folder in
  `Genomes/` or like another folder an index records, and builds an index again when its folder's
  FASTA files have changed (`Removing index ...`), as `complete-search` does, with the same
  messages. Of a folder in `Genomes/` it records nothing, as before, so a downloaded index, and a
  genome renamed together with its indexes, are served as before. A Direct URL download under the
  name of a genome with an index or a variant genome installed is refused before anything is
  fetched, by the rule a Settings upload has (remove them, or use another name), and the Settings
  page says so at once. `--source ucsc` and the HuggingFace download give the same assembly by
  name, and are not refused.
- **A variant search reads the variant index at `--index-path`.** `--index-path` served the
  reference index alone: the variant index and its INDELS index were looked for in the working
  directory's `genome_library/` only, so a variant search of a complete library kept elsewhere
  (a shared, read-only install) stopped with `mv: cannot stat` and `Genome enrichment failed`
  when the dataset's VCFs were absent, and with them it enriched the genome and built a second
  variant index in the working directory, which every later run there took, also when the one at
  `--index-path` was newer. **Affected:** `complete-search --vcf --index-path <library>` of a
  library other than the working directory's, in 2.8.0 too; the web passes no `--index-path`.
  **What to do:** a variant index such a run built in the working directory can be removed. The
  variant index, with its INDELS index beside it, is now read from the `--index-path` library first
  and then from the working directory's, and the run says which library it read
  (`Variant index library: ...` in `log_verbose.txt`). Measured on the synthetic genome: a library
  elsewhere with no VCFs here, the same with the VCFs here, a K=8 library with `--prune-pam`, and
  a stale index of the same name in the working directory are each byte-identical to the search
  with that library in place, and build nothing. A search builds no index in an `--index-path`
  library: a variant index there without its INDELS index stops the run with a message.
  Dictionaries, sample lists and the genome are still read from the working directory, and a
  search with the library in the working directory, with or without `--index-path`, is unchanged.
- **A search refuses an index that lacks a chromosome's trees, instead of losing that
  chromosome's sites.** A search took an index holding one `.bin` for a whole one: a partial copy,
  a build that an earlier release left stopped between chromosomes, or an index of a genome folder
  that gained a FASTA since, gave no row on the missing chromosome, at exit 0 and with no message,
  on the CLI and the web, and `build-index-only` said "already present". Measured on the synthetic
  genome: a reference index without chr2's tree gave 6 rows instead of 7, a variant index without
  chr1's another result, and a FASTA added after the build was not searched. **Affected:** searches
  of such an index, in 2.8.0 too; a stock install is not (every index on the lab's hg38 install, 32
  reference and variant indexes of NRG, NGG, NRR, NNN and K=8, holds every tree). **What to do:**
  nothing, unless a search now stops on an index: remove it and build it again, or download it
  again, as the message says. Before it opens an index (with bulges, or for a variant search
  without bulges and without the whole enriched genome), the search now checks that it holds the
  trees of every FASTA of the genome folder, and a variant index those of every chromosome of its
  dataset, and stops naming the FASTAs that lack them; `build-index-only` does the same before it
  says "already present". A FASTA with no site of the index's PAM (all N, shorter than a window, no
  PAM) has no tree in a healthy index and is not counted, by the engine's own rule for a site, and
  neither is a variant chromosome the genome folder has no FASTA of. A gap in a chromosome's
  numbered trees, and a full last tree of a chromosome that holds more sites, count as missing
  trees too; a chromosome's last trees lost after a tree that is not full cannot be told from the
  files. The check reads the first line of each FASTA, and a FASTA in full only where it has no
  tree or a full last tree: under a second on hg38.
- **An install without the genotype store runs the population-level analysis, rather than a
  per-sample one that loses off-targets.** On a phased dataset the CLI chose the per-sample
  analysis from the registry's data type alone, so after `download --what index --no-genotypes`,
  or a download whose genotype store did not come, it ran per-sample with no genotypes to read:
  it joined no two SNPs of one site, dropped the off-targets that only two SNPs of a sample
  together reach, and still printed "CONFIRMED cis phasing + named carrier samples". `--per-sample`
  did the same. Measured on the synthetic genome with a site whose PAM only two cis SNPs make (1
  mismatch, CFD 0.6, carrier S4 on a full install): the no-store default reported it at a 2-mismatch
  bulged alignment, CFD 0.095. **Affected:** CLI searches of such an install, in 2.8.0 too; the web
  already offers population-level alone there. **What to do:** run them again, or install the store
  (`download --what index` without `--no-genotypes`) for carriers. A dataset with a registry and
  neither a genotype store with anything in it nor the classic per-sample dictionaries now runs
  population-level, by default and under `--per-sample`, with one `WARNING` on stdout that names
  the dataset and the download. On the synthetic genome the no-store default is then byte-identical
  to its `--population-level` run and reports that site at 1 mismatch, CFD 0.6, with no carriers; it
  lists the same sites with the same scores as a full install's `--population-level` run, without
  the carrier columns. A full install, and one with the classic dictionaries, are unchanged. The
  download help, the README, docs/PRECOMPUTED_INDEXES.md and docs/TROUBLESHOOTING.md say what such
  an install runs.
- **An index is taken only for a search on the side of the guide it was built for, a dedicated
  one too, and an all-N index that records no side is taken for a 3' PAM only at the published
  widths.** 2.8.1 takes an all-N index of a 3' PAM's width, and one with no `.pam_build` for a 3'
  build. A search before 2.8.1 records no side, so after a 5' all-N search (`23bp-NNNN-NO-PAM`,
  which builds `NNNN_2_<genome>`) a 3' search of that width (`NGCG`) read the 5' index: on the
  synthetic genome it found no site against 2 on a dedicated index, which 2.8.0 built. This is new
  in 2.8.1. The side check also covered all-N indexes alone, while a dedicated index is named
  after its PAM by PAM files of both sides (`N...NTTTV 4` and the 5' `TTTV-Cas12a` both build
  `TTTV_<N>_<genome>`): a 5' TTTV search of a 3' build reported nothing, or a site that is in no
  chromosome, at exit 0, as in 2.8.0, and `build-index-only` with the other side's file recorded
  that side over it. **Affected:** a search beside an all-N index of a width other than 3 and 8
  that a search before 2.8.1 built, and a search beside a dedicated index built from the other
  side's PAM file; no published index and no shipped PAM file. **What to do:** run those searches
  again. Now every index's recorded side must be the search's side, and a search records the side
  of every index it builds, dedicated ones too (build-index-only always did). An all-N index with no
  record serves a 3' PAM at widths 3 and 8 (the published `NNN_3` and the K=8 archive, which has
  no `.pam_build`) or under the search PAM's own name, as before; at any other width the search
  builds its dedicated index, as 2.8.0 did, and under `--index-path` it stops and says why. A
  dedicated index with no record is used as before. Where the index of the build name is in the
  way the search stops with a message, and `build-index-only` refuses another side's PAM file over
  a recorded index, and a 3' file over an unrecorded all-N index of another width (a 5' file is
  recorded, as the messages say to do). Measured on the synthetic genome: the NGCG search after the
  5' NNNN search, and the 5' TTTV search beside a 3' `TTTV_2` and a 5' `TTTV_3`, now equal a
  dedicated index; NRG on an unrecorded `NNN_3`, and on an unrecorded K=8 with `--prune-pam`, is
  unchanged. A dedicated index built for the other side's PAM by a search before 2.8.1 records
  nothing and cannot be told apart: remove it if you built one from such a PAM file.
