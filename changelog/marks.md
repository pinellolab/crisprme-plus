<!-- marks track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the docs track: five amendments to entries already in `## [2.8.1]`, then new entries.**
Most of this track's work narrows or completes behavior that is itself new in 2.8.1, which 2.8.0
never had: the unfinished mark of the fake-indel genome, the one-chromosome view of
`complete-test --chrom`, and the one-line hint of `validate-test`. An entry of its own would
describe a failure no release had, and the entries already there would contradict what 2.8.1
does. So each is an amendment: replace the quoted text (line numbers are those of `CHANGELOG.md`
at aa9c54b) with the text after it. All five entries are under `### Fixed`.

1. **"The SNP overlay can no longer leave a fake-indel FASTA cut short"**, lines 573-574. Replace

   > Either way each FASTA is written in a hidden directory beside the original and renamed over it,
   > with the same contents, mode and owner, and an overlay that fails

   with

   > Either way each FASTA is written in a hidden directory beside the original and renamed over it,
   > with the same contents, mode and owner (one in place none of whose bases changes, as on a
   > dataset this release enriched, is left as it is), and an overlay that fails

2. **"A `build-index-only --vcf` run stopped after the enrichment no longer leaves the INDELS index
   blind to the SNPs beside an indel"**, lines 597-605. Replace

   > it now gets them as the run begins, whenever its INDELS index is still to be built, and while
   > they go on the fake-indel genome is marked unfinished, as the tiers are (next entry), so a
   > search refuses it until a `build-index-only` run has finished them; an INDELS index that a
   > search built from it meanwhile (one that got past its check of the mark before the run set it)
   > is removed once they are on, and built again; that search's own results come from it, with exit
   > status 0, so run such a search again once the `build-index-only` run has finished. On success
   > the files, their contents, stdout lines (the overlay's still print just before the INDELS index
   > is built) and exit status are as before, but for a re-run on a fake-indel genome that carries
   > its codes already: its overlay changes no base, and where its lines said

   with

   > it now gets them as the run begins, whenever its INDELS index is still to be built. From just
   > before the first FASTA that changes until the overlay ends, the fake-indel genome is marked
   > unfinished, as the tiers are (next entry), and a search that would build an INDELS index from
   > it (the variant index it searches has none yet) refuses it until a `build-index-only` run has
   > finished the codes. `complete-search` refuses it before it starts, unless the dataset has
   > another variant index that serves the search with its INDELS index in place: the search then
   > goes on, and if the variant index it picks has no INDELS index, it fails as it comes to build
   > one, before it builds anything (`CRISPRme run failed`, with that ERROR in its `log_error.txt`),
   > as a search launched from the web does (next entry). A search whose INDELS index is in place
   > does not read that genome, and the mark does not hold it up. A fake-indel genome none of whose
   > FASTAs changes, as when a second PAM or bulge budget is built on a dataset this release
   > enriched, is neither marked nor rewritten. If a search built, from that genome meanwhile, the
   > INDELS index of the PAM and bulge budget that the `build-index-only` run builds (one that got
   > past its check of the mark before the run set it), that index is removed once the overlay ends,
   > and built again; that search's own results come from it, with exit status 0, so run such a
   > search again once the `build-index-only` run has finished. An INDELS index of another PAM or
   > bulge budget that a search built that way is not removed: remove it
   > (`genome_library/<PAM>_<N>_<genome>+<vcf>_INDELS`) once that run has finished, and the next
   > search or `build-index-only` run with that PAM and budget builds it again. With
   > `CRISPRME_INDEL_SNP=0`, which puts no SNP code on, `build-index-only` does not build an INDELS
   > index from a fake-indel genome that a stopped run left marked: it stops with an ERROR saying to
   > run it again with the overlay on (the default). On success the files, their contents, stdout
   > lines (the overlay's still print just before the INDELS index is built) and exit status are as
   > before, but for a re-run on a fake-indel genome that carries its codes already: its overlay
   > changes no base and rewrites none of its FASTAs, and where its lines said

   The rest of the entry stays as it is.

3. **"A `build-index-only --vcf` run stopped while it writes the dict-less tiers or the indel
   genotype store no longer leaves a part of them that searches use as the whole"**, lines 638-639.
   Replace

   > The same goes for the enriched genome and the fake-indel genome marked unfinished (above).

   with

   > The same goes for the enriched genome marked unfinished, and for the fake-indel genome marked
   > unfinished (above) when the search would build its INDELS index from it. `complete-search`
   > names each marked part by its real path, the one a symlink there points to.

   Line 631 ("every search ends with an ERROR while the mark is there") is about the tiers and the
   indel genotype store, and stays as it is.

4. **"`complete-test --chrom chr22` searched the whole genome, and the source installer started on
   macOS"**, lines 362-364. Replace

   > to its `fake<chrom>` INDELS sub-index and to its indel log, so that the search does not load
   > another chromosome's search or post-analysis data (the dataset-wide tier folders are linked
   > whole). Nothing is copied,

   with

   > to its `fake<chrom>` INDELS sub-index, to its fake-indel genome and to its indel log, so that
   > the search does not load another chromosome's search or post-analysis data (the dataset-wide
   > tier folders are linked whole). For a variant index without an INDELS index (that of a dataset
   > without indels, or one whose build stopped before it), the search builds a throwaway one for
   > that chromosome inside the folder, or finds no indel there to index, and leaves none in
   > `genome_library`. Before it makes the folder, the smoke refuses each part of its dataset that
   > `build-index-only` has marked unfinished and `complete-search` would refuse (see below), in one
   > ERROR line that names its real path, with nothing searched and no `Results/crisprme-test-out`
   > left behind. Nothing is copied,

5. **"The documented install self-test tested nothing, and `validate-test` then failed with a
   traceback"**, lines 282-283. Replace

   > naming both commands for the chromosome it was asked to check (the chr22 pair without
   > `--chrom`), and exits 1;

   with

   > naming both commands for the chromosome it was asked to check (the chr22 pair without
   > `--chrom`) and saying that in Docker their variables go with `-e`, and exits 1;

### Fixed
- **A failed `complete-test` now ends with one line, which names an error log only where one was
  written.** It ended with `crisprme.py`'s own `OSError` traceback, which pointed to
  `Results/crisprme-test-out/log_error.txt`: a search refused before it starts writes no log there,
  and the benchmark mode writes one for each benchmark, in `Results/crisprme-test-out_<name>`. It
  now ends with one line naming each error log, not empty, that this run's searches wrote
  (`complete-test` refuses to start while a `Results/crisprme-test-out*` folder of an earlier run is
  there), or else saying to see the messages above, and exits 1, as before. Covered by
  `PostProcess/test_complete_test.py`.
- **`targets-scoring --index-path` refuses a dict-less registry that `build-index-only` has not
  finished writing.** The carrier overlay read `Dictionaries/registry_<vcf>` as it found it, also
  while a `build-index-only` run was writing it chromosome by chromosome, or after one had stopped
  partway, and for each chromosome not written yet it reported 0 overlapping variants, with no
  warning. It now refuses a registry marked unfinished (see the dict-less tiers entry), as
  `complete-search` and `publish-index` do: it stops with a one-line `Error:` saying to let that run
  finish, or to run `build-index-only` again to its end with its `--samplesID`. A
  registry that a run of an earlier release left partway carries no mark and is not detected. A
  missing registry is now also reported in one line instead of a traceback. Output for a finished
  registry is unchanged. Covered by `PostProcess/test_targets_scoring.py`.

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
  `-e`. The run and its exit status are unchanged. Covered by `PostProcess/test_complete_test.py`.
