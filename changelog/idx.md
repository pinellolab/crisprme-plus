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
  the search runs with none, and what to pass. A variant index at `--index-path` is still read from
  `genome_library/` alone: with the dataset's VCFs the search rebuilds it there, as it does with
  explicit flags, and without them it stops and says so.
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
