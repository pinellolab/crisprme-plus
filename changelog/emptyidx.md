<!-- emptyidx track: [2.8.1] entry for CHANGELOG.md, to be merged by the docs track. -->

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
  (`Genomes/<genome>+<dataset>_INDELS/`) by the same rule and skips the chromosome when no window
  holds a target. Otherwise (a window holds one, so the index is incomplete; that genome is not
  there; or a PAM-less index is searched for another PAM) the search stops with an ERROR that says
  which, and gives the command that rebuilds the index (`rm -rf ... && crisprme.py build-index-only
  ...`) or, for a downloaded index, downloads it again, where it used to show the engine's message.
  Measured through the real CLI: the TTTV, TTTN and NGCG searches now exit 0, both on an index built
  with the fix and on one built before it, and every result file is the same as that of the same
  code run on CRISPRitz 2.9.0's search engine, which has no such refusal (but for timestamps, and
  the engine version `run_manifest.json` records); for TTTV the result tables are also the same as
  2.8.0's. NRG, NGG and NNGRRT searches, NNN searches and `--prune-pam` searches over 3- and
  8-nt-flank PAM-less indexes give byte-identical results to before, and the indexes
  build-index-only writes are the same file for file, apart from the new record where a sub-index
  holds no target. The printed rebuild command works as given. **No published index is affected:**
  in the copies of the published indexes on the build machine, every indel sub-index holds a `.bin`,
  except in the superseded `NRR_3_hg38+hg38_mega` archive (see its entry above), whose own indel
  logs are empty, so a search skips those chromosomes as before. **What to do:** nothing for an
  index built with 2.8.1, or for one built before whose fake-indel genome is still in place;
  otherwise run the command the ERROR gives. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.
