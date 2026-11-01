<!-- emptyidx track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

### Fixed

- **A variant search no longer fails when a chromosome's indels hold no site of the PAM** (a
  regression of this release's CRISPRitz pin, like the SNP-only chromosome above, with another
  trigger). `build-index-only --vcf` builds one indel sub-index per chromosome from the short
  windows around that chromosome's indels. When none of them holds a site of the PAM, `buildTST`
  writes no `.bin` there, and CRISPRitz 2.9.1 refuses to search an index directory without one, as
  it should for a broken one. The search therefore stopped with exit status 1, reported at
  "Post-analysis SNPs" with `ERROR: no .bin index files found in '.../TTTV_2_fakechr3/'`, where
  2.8.0 returned the correct report. **Affected:** variant indexes built for a PAM whose sites are
  rare enough, on a dataset with few indels on some chromosome: the shipped Cas12a (TTTV, TTTN) and
  SpCas9-VRER (NGCG) PAMs on a test genome whose chr3 carries a single deletion; NRG, NGG and SaCas9
  (NNGRRT) had a site there. The same failure followed a search that rebuilt a missing INDELS index
  itself. `pool_index_indels.py` now writes `.no_targets` into such a sub-index, after checking by
  `buildTST`'s own rule (the check build-index-only already makes) that no window holds a target,
  and the search skips that chromosome as it skips one without indels. The record travels with the
  index, so a copied or downloaded index needs nothing else. An index built before this fix has no
  record: the search then checks the fake-indel genome the index was built from
  (`Genomes/<genome>+<dataset>_INDELS/`) by the same rule, for the PAM the index was built with
  (for a PAM-less index searched for another PAM, such as `NNN_2` serving NRG, the one its
  `.pam_build` records), and skips the chromosome when no window holds a target. Otherwise (a
  window holds one, so the index is incomplete, or that genome is not there) the search stops with
  an ERROR that says which and how to rebuild the index, where it used to show the engine's
  message: the command `rm -rf ... && crisprme.py build-index-only ...` with the PAM file the index
  was built with (for a PAM-less index, a file of its PAM in `PAMs/` or beside the search's PAM
  file), and, for a downloaded index, also the download command. Where there is no such PAM file,
  the ERROR gives the steps in words instead of a command: remove the INDELS index and run the
  search again, which rebuilds it from the fake-indel genome, or, without that genome, the PAM line
  to write into a PAM file for build-index-only. Measured through the real CLI: the TTTV, TTTN and
  NGCG searches now exit 0, both on an index built with the fix and on one built before it, and
  every result file is the same as that of the same code run on CRISPRitz 2.9.0's search engine,
  which has no such refusal (but for timestamps, and the engine version `run_manifest.json`
  records); for TTTV the result tables are also the same as 2.8.0's. NRG, NGG and NNGRRT searches,
  NNN searches and `--prune-pam` searches over 3- and 8-nt-flank PAM-less indexes give
  byte-identical results to before, and the indexes build-index-only writes are the same file for
  file, apart from the new record where a sub-index holds no target. Each rebuild the ERROR gives,
  followed as given, leads to the same result files as an intact index (measured for an NRG index,
  and for NRG served by the 3- and 8-nt-flank PAM-less indexes, with and without the fake-indel
  genome). **No published index is affected:** in the copies of the published indexes on the build
  machine, every indel sub-index holds a `.bin`, except in the superseded `NRR_3_hg38+hg38_mega`
  archive (see its entry above), whose own indel logs are empty, so a search skips those
  chromosomes as before. **What to do:** nothing for an index built with 2.8.1, or for one built
  before whose fake-indel genome is still in place; otherwise follow the ERROR. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.
- **A variant search with bulges no longer loses the indel off-targets of a chromosome whose indel
  sub-index is missing** (in 2.8.0 too). Such a chromosome was searched in its fake-indel genome
  without bulges, a fallback meant for searches without bulges: the run exited 0 with no message,
  and the bulged indel off-targets of that chromosome were missing from the report. It happened in
  two ways. (1) A search that rebuilt a missing INDELS index for a PAM-less index serving another
  PAM (`NNN_<N>` serving NRG, or `NNNNNNNN_<N>` under `--prune-pam`) built it with the search's PAM
  file, so it went into `<that PAM>_<N>_<genome>+<dataset>_INDELS` (for example
  `NRG_2_..._INDELS/NNN_2_fake<chr>`), where the search does not look; every later such search did
  the same. (2) An INDELS index copied in part. Measured through the real CLI on a test genome
  (NRG, 4 mismatches, 1 DNA and 1 RNA bulge): way 1 lost 5 of the 2446 rows with alternative
  alignments and 190 of the 191 indel+SNP co-occurrence rows of chr3 (2447 -> 2442 rows on 2.8.0),
  way 2 lost 3 rows with alternative alignments, and the integrated results changed with both. Now
  (1) that rebuild uses the PAM the index was built with (from its `.pam_build`, else of the
  search's guide length), and is byte-identical to the INDELS index build-index-only builds, so the
  search's result files are those of an intact index; and (2) a search with bulges stops at a
  chromosome whose indel log, or fake-indel genome, shows indels but which has no sub-index, with
  an ERROR that says how to rebuild the INDELS index (as above). A search without bulges still
  searches the fake-indel genome, and a chromosome without indels is skipped as before. Normal runs
  are unchanged: a build gives every chromosome with indels its sub-index, as it does in each of
  the 50 INDELS indexes on the build machine that have their indel logs beside them. **Affected:**
  searches with bulges through a PAM-less index serving another PAM whose INDELS index was missing,
  and searches with an INDELS index copied in part. **What to do:** a job was affected by (1) when
  it allowed bulges and its `log_verbose.txt` shows `Variant index present (NNN_...)` (or
  `NNNNNNNN_...`) for another PAM, followed by `Indexing INDELs in ...` lines: run it again with
  2.8.1, and remove what it left under the other PAM's name, the `NNN_<N>_fake<chr>` (or
  `NNNNNNNN_<N>_fake<chr>`) directories in `genome_library/<that PAM>_<N>_<genome>+<dataset>_INDELS`,
  with that directory when nothing else is in it. For (2), a search now says which index is
  incomplete; once it is rebuilt, run again the searches made with it. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.
