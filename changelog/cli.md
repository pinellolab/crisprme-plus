<!-- cli track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the docs track: one amendment to an entry already in `## [2.8.1]`, then seven new
entries.** The entry "A failed `complete-test` now ends with one line, which names an error log
only where one was written" (under `### Fixed`, lines 1164-1171 of `CHANGELOG.md` at 38a7c4e)
says, in parentheses, why the logs it names are this run's; the `complete-test` re-run entry below
changes that reason, as each mode now runs beside the other's folders. Replace

> (`complete-test` refuses to start while a `Results/crisprme-test-out*` folder of an earlier run
> is there)

with

> (`complete-test` refuses to start while a folder of an earlier run of its mode is there, and
> names only that mode's: `Results/crisprme-test-out_<name>` in the benchmark mode,
> `Results/crisprme-test-out` in the smoke)

All seven entries are pre-existing (2.8.0 did the same; the card's line is the live card's) and
change no search result.

### Fixed

- **A failed `complete-search` ends on its summary of the failure, with exit status 1 and no
  Python traceback.** It printed where and why the run failed (`CRISPRme run failed.`, the stage,
  the last lines of `log_error.txt`, the log's path), then raised an `OSError`, so every failed
  search, the refusal of an index a newer CRISPRitz built among them, went on with `Traceback (most
  recent call last): ...` and an `OSError` line that only repeated the log's path.
  `assembly-search` added a traceback of its own when a haplotype's search failed; it now ends in
  one line, `Error: Haplotype search failed for --output <name>! See <log> for details`, which the
  web's assembly job reads as the job's failure, as it read the traceback. Exit statuses are
  unchanged, and a failed search still removes its annotation temp folders. Pre-existing, in 2.8.0
  too. Covered by `PostProcess/test_build_index_only_atomic.py`,
  `PostProcess/test_annotation_tmpdirs.py` and `PostProcess/test_assembly_haplotype_genomes.py`.
- **The documented self-test no longer reports a HuggingFace failure for files the data repo does
  not host.** `CRISPRME_ALLOW_ONDEMAND_BUILD=1 crisprme.py complete-test --chrom chr22
  --vcf_dataset 1000G` always printed `HuggingFace fetch failed for
  ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz (404 Client Error ...
  Entry Not Found for url ...); falling back to original source`, then downloaded the VCF from EBI
  and passed. `complete-test` tries the Hugging Face data repo first for the chromosome's FASTA,
  the VCFs and the sample lists it downloads (the annotations come from GitHub), and
  `lucapinello/crisprme-data` holds no raw 1000 Genomes or HGDP VCF (its `vcfs/` holds the
  source VCFs of the three shipped panels), nor `samplesIDs.HGDP.txt`, nor the unpacked FASTA of
  any chromosome but chr22, so `--vcf_dataset HGDP` and `--chrom chr21` printed it too. A file the
  Hub answers it does not hold is now said in one line on stdout, `<file> is not on HuggingFace
  (lucapinello/crisprme-data); downloading it from its original source` (`from UCSC` for a
  FASTA). Any other failure to fetch a file, such as no connection or a checksum that does not
  match, is still a warning on stderr, word for word as before. What is downloaded, and from
  where, is unchanged. Pre-existing, in 2.8.0 too. Covered by `PostProcess/test_complete_test.py`.
- **`setup` no longer reports a HuggingFace failure for files the data repo does not host
  either.** `crisprme.py setup`, which downloads the legacy database of the CRISPRme web server,
  tries the Hugging Face data repo first for the chromosome's FASTA, the 1000 Genomes and HGDP
  VCFs and the two sample lists, and printed `HuggingFace fetch failed for <file> (404 Client
  Error ... Entry Not Found for url ...); falling back to original source` (`to UCSC` for a FASTA)
  on stderr for each one the repo does not hold: every VCF, `samplesIDs.HGDP.txt` and the
  unpacked FASTA of any chromosome but chr22, so `setup --chrom chr22` printed it three times,
  and a whole-genome `setup` 47 times. It now says each in one line on stdout, as
  `complete-test` does: `<file> is not on HuggingFace (lucapinello/crisprme-data); downloading it
  from its original source` (`from UCSC` for a FASTA). Any other failure to fetch a file is still
  a warning on stderr, word for word as before. What is downloaded, and from where, is unchanged.
  The setup guide, `docs/crisprme_data_setup_051826.md`, said `setup` pulls the genome and the
  variants from their original sources, a single chromosome's FASTA from UCSC; it now says `setup`
  takes two files from the data repo, `chr22.fa` for `--chrom chr22` and `samplesIDs.1000G.txt`,
  and everything else from UCSC, EBI, Sanger or GitHub. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_setup_legacy_database.py`, a new step of `unit-tests.yml`.
- **The `lucapinello/crisprme-data` card says what its `genomes/hg38/` holds.** Its layout, kept
  from the live card, gave the folder as a per-chromosome FASTA mirror, `chr1.fa … chrX.fa`. The
  repo holds one gzipped FASTA per hg38 sequence (`chr1.fa.gz` … `chrY.fa.gz`, `chrM.fa.gz` and
  the alt, random and unplaced contigs: 455 files, which `download --what genome` unpacks) and one
  unpacked FASTA, `chr22.fa`, which `complete-test` and `setup` fetch as is. The layout now lists
  both, as the repo's file list gave them on 2026-11-04. The live card is corrected when the
  release uploads `docs/hf/crisprme-data_README.md`.
- **`complete-test` run again over an earlier run's results names every folder in its way, and
  exits 2.** It named only the first such folder in sorted order, so after a benchmark-mode run it
  named `Results/crisprme-test-out_cas12a_hbg`, one of six, and a user who deleted it was stopped
  five more times; it stopped the prebuilt-index smoke on those six folders too, though the smoke
  writes only `crisprme-test-out`; and it exited 0, which a script or a CI step reads as a pass
  though nothing ran. It now names, in one line, every folder there of the run's own mode (the
  benchmark mode's `crisprme-test-out_<benchmark>`, the smoke's `crisprme-test-out`) and exits
  with status 2 before it downloads or searches anything; `crisprme.py` passes the 2 on, as for a
  refused `--chrom`. A folder of the other mode no longer stops a run, and the `CRISPRme complete
  test encountered an error` line of a failed run names only its own mode's error logs. README
  §2.2.2 says how to run a mode again. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_complete_test.py`.
- **A search's `log_verbose.txt` keeps every `Searching for INDELs in` and `Search ended for
  INDELs in` line.** `pool_search_indels.py` printed them from its pool's workers without flushing
  them, and leaving the pool terminates each worker it finds idle, with the lines it had not yet
  written: a search's log missed 3 to 23 of them, a different number each run.
  `pool_index_indels.py`, which builds a missing INDELS index for a search and for
  `build-index-only`, lost its `Indexing INDELs in` lines the same way. Both flush them now; no
  result changes. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.
- **A source install gets the web stack the Docker image and the Bioconda package run.**
  `environment.yml`, which `install_from_source.sh` installs into and which says it mirrors the
  Dockerfile's validated dependency set, left `flask`, `flask-caching`, `flask-compress` and
  `werkzeug` unpinned, and `cachelib` out. The Dockerfile and the Bioconda recipe pin them (`flask`
  3.1.3, `flask-caching` 2.4.1, `flask-compress` 1.24, `werkzeug` 3.1.8, `cachelib` 0.17.0), as
  an unpinned `flask-caching` once drifted to a release that crashed the web app at import; solved
  for linux-64 on 2026-11-03, the file took `flask-caching` 2.5.1, `flask-compress` 1.26 and
  `werkzeug` 3.1.9. It now pins the five as the Dockerfile does, and a dry-run solve of it for
  linux-64 gives the very builds of the release image. An environment created from the earlier
  file keeps the versions it got: create it again to get these. Pre-existing, in 2.8.0 too.
  Covered by `PostProcess/test_install_from_source.py`.
