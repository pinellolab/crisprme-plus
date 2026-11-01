### Fixed
- **Building a second PAM or bulge budget on an enriched dataset no longer refuses every search of
  that dataset while it runs.** When `build-index-only` found the fake-indel genome already in place
  (`Genomes/<genome>+<vcf>_INDELS`), it marked it unfinished for the whole SNP overlay and rewrote
  every FASTA, even when no base changed, as on a genome this release had coded already. Meanwhile
  every variant search of the dataset was refused, from the command line and from the web ("The
  search was not run"), including searches whose INDELS index was in place and that never read the
  fake-indel genome. Genome-wide this lasts minutes, by extrapolation from a 50 Mb test. The overlay
  in place now leaves a FASTA none of whose bases changes as it is, and marks the genome only just
  before the first FASTA that changes. A mark that a stopped run left is still cleared once the
  overlay ends. A search now checks the mark only where it would read the fake-indel genome, that
  is, when it has to build its own INDELS index. Such a search is still refused: by
  `complete-search` before it starts, and by a search launched from the web at that step, before it
  builds anything. A search whose INDELS index is in place (for its PAM, with enough bulges) runs.
  As before, if a search built an INDELS index from the genome while the codes went on (one that got
  past its check before the mark was set), `build-index-only` removes that index once the codes are
  on and builds it again; a search that finds that index in place before then uses it. For the same
  reason, `build-index-only` with `CRISPRME_INDEL_SNP=0` no longer builds an INDELS index from a
  fake-indel genome that a stopped run left half coded. Such an index would be blind to the SNPs on
  the FASTAs not coded yet, and a search that finds an INDELS index in place no longer looks at the
  mark. The run now stops before building, with an ERROR saying to run `build-index-only` again with
  the overlay on, and its end-of-build ERROR for that mark says which searches refuse it. Search
  results are unchanged, and so are the files a build leaves, except that a fake-indel FASTA
  already in place that needs no SNP code is no longer rewritten with the same sequence on one
  line. Covered by `PostProcess/test_build_index_only_atomic.py` and
  `PostProcess/test_index_resolution.py`.
- **`complete-test --chrom` now works on a variant index that has no INDELS index.**
  `build-index-only` builds no INDELS index for a dataset without indels, and a build stopped before
  its last step leaves the variant index without one. In both cases the one-chromosome smoke failed
  with "INDELS index ... is missing and cannot be built ... re-download it", for an index that was
  never published, and a traceback, although the whole-genome test passed. The smoke's per-run view
  now holds `Genomes/<genome>+<vcf>_INDELS` with that chromosome's fake-indel genome. The search
  builds a throwaway INDELS index for that one chromosome inside the view, or finds no indel there
  to index. The rows match the whole-genome test's rows on that chromosome, and nothing is left in
  `genome_library`. Covered by `PostProcess/test_complete_test.py`.
- **`complete-test --chrom` no longer passes on a dataset that every search refuses.** The
  one-chromosome smoke runs `complete-search` from a per-run view of the working directory, and the
  view held no `Genomes/<genome>+<vcf>`. So when `build-index-only` had left the enrichment or the
  fake-indel genome marked unfinished, the smoke ended `complete-test PASSED`, while the
  whole-genome test and every real search of the dataset were refused. For a marked dict-less tier
  or indel genotype store, the refusal named a path inside the view, which was gone once the run
  ended, and its `rm -rf` advice pointed there. The smoke now checks the same unfinished parts that
  `complete-search` refuses, in the working directory itself, before it makes the view. Each marked
  part gets one ERROR line with its real path, a symlink followed. Nothing is searched, and no
  `Results/crisprme-test-out` is left behind to block the next run. `complete-search` also names a
  dict-less tier, indel genotype store, enrichment or fake-indel genome that is a symlink by the path
  it points to. Separately, any failed `complete-test` now ends with one line instead of
  `complete-test`'s own `OSError` traceback, with exit status 1 as before. That line names the
  search's error log only when one was written, under its real folder name
  (`Results/crisprme-test-out*/log_error.txt`), and otherwise says to see the messages above, since
  a refusal before the search starts writes no log. Covered by `PostProcess/test_complete_test.py`.
- **`targets-scoring` refuses a registry that `build-index-only` has not finished writing.** The
  carrier overlay (`--index-path`) used `Dictionaries/registry_<vcf>` even while it was marked
  unfinished. For a chromosome whose tiers were not written yet it reported 0 overlapping variants,
  with no warning, where `complete-search`, the web search and `publish-index` all refuse that mark.
  It now stops with a one-line `Error:` saying to let that run finish, or to run `build-index-only`
  again to its end with its `--samplesID`. A missing registry is now also reported in one line
  instead of a traceback. Output for a finished registry is unchanged. Covered by
  `PostProcess/test_targets_scoring.py`.

### Changed
- **`--per-sample` on an index without the dict-less tiers now says what it reports.** An index
  built without `--samplesID` has only the classic per-sample dictionaries, and no
  `Dictionaries/registry_<vcf>/` or `genotypes_<vcf>/`. This covers the index that benchmark-mode
  `complete-test` builds and a build from the web Settings page. On such an index, `complete-search
  --per-sample` still promised "CONFIRMED cis phasing, named carrier samples and exact joint AF",
  but wrote no phase-confirmation or SNP+SNP co-occurrence table. The banner now says that this
  dataset's results name the carrier samples, with the allele frequencies the VCF gives, but have no
  CONFIRMED/PUTATIVE cis phasing, no SNP+SNP co-occurrence and no exact joint AF. It also says to
  rebuild with `build-index-only --samplesID` or download a prebuilt index. The banner of an index
  with the tiers, and every result, are unchanged. Covered by
  `PostProcess/test_build_index_only_atomic.py`.
- **`complete-test` now says when it ignores `--vcf_dataset`.** The benchmark mode, which
  `--vcf_dataset` configures and `validate-test` checks, runs only with
  `CRISPRME_ALLOW_ONDEMAND_BUILD=1` set. A Docker user who put the variables in front of `docker
  run` set them on the host only, and `complete-test` silently ran the prebuilt-index smoke
  instead. It now prints one line on stdout saying so, and that in Docker the variable goes with
  `-e`. The run and its exit status are unchanged. The one-line hint that `validate-test` gives when
  there is nothing to check now says the same: `(Docker: pass them with -e)`. Covered by
  `PostProcess/test_complete_test.py` and `PostProcess/test_validate_test.py`.
