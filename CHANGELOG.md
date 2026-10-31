# Changelog

All notable changes to CRISPRme are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When cutting a release, move items out of `[Unreleased]` into a new dated
version section and update the link-reference footer. See `docs/RELEASING.md`
and the `release-crisprme` skill.

## [Unreleased]

## [2.8.1] - 2026-10-31

### Added
- **`--prune-pam`: one index can serve many PAMs.** A *PAM-less* index stores the raw genomic flank
  instead of a matched PAM, and `--prune-pam` applies the PAM you pass during the search, discarding
  a leaf whose stored flank cannot match it. One index therefore answers any PAM no longer than its
  flank, instead of one index per PAM. Verified genome-wide over HGDP + 1000 Genomes (4,091 samples,
  NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per-sample): a 3-nt-flank PAM-less index, an 8-nt-flank
  PAM-less index, and a dedicated NRG index produce **byte-identical** reports — the same 189,950
  rows, every indel row justified by its own variant window. Pruning changes only which leaves are
  written out, never how a guide is aligned, so it cannot move an alignment or a score. The cost is
  size, since nothing is excluded at build time: 110 GB against 28 GB for the HGDP + 1000 Genomes
  variant index. Requires `--index-path`, and is command-line only — the web interface never passes
  the flag. A shorter PAM off a wider flank is served for a 3' PAM only; a 5' PAM requires the
  flank width to equal the PAM length, and the engine refuses any other width pairing rather
  than mis-parsing the index. The `NNN_3` indexes are published alongside the dedicated
  ones; a flank wider than 3 nt is not readable by earlier releases, so those live in a separate
  repository. Needs a CRISPRitz build that supports `--prune-pam`.
  See [`METHODS.md` §9](METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering).
- **A second dataset repo for wide-flank PAM-less indexes.** A PAM-less index with a 3-nt flank
  stores what a dedicated 3-nt PAM index stores, so any release can read it and `NNN_3` stays in
  `lucapinello/crisprme-data`. A wider flank cannot be read by a release predating `--prune-pam`, so
  those are published to `lucapinello/crisprme-data-general-index`, which older defaults never point
  at. `download` and `publish-index` route on the flank width parsed from the index name — not on a
  list of names, so it cannot be forgotten at publish time — `download` falls back to the default
  repo when the routed one lacks the archive (and, if both lack it, says it looked in both), and the
  catalog listing merges both repos, taking each name's row from the copy `download` would fetch, so
  what is offered is what can be fetched. If the default repo cannot be read the catalog is empty,
  as before, rather than offering only the slower wide-flank indexes. Dedicated-PAM and `NNN_3`
  indexes are listed and downloaded from the default repo exactly as in 2.8.0, never looked up in
  the new one. `--hf-repo` / `CRISPRME_HF_REPO` still override. `download` takes
  an index's genotype companion and reference index from that index's own repo, so the new repo
  carries its own copy of `genotypes_hg38_HGDP1kGP.tar.gz`, byte-identical to the main repo's.
- **`indexes/SHA256SUMS` in both dataset repos**, written by `scripts/hf_sha256sums.py` from the
  checksums the Hub already records (a Git LFS object is addressed by its SHA-256; the oids match
  the independently computed `build-inputs/*.SHA256SUMS`), so `sha256sum -c SHA256SUMS` verifies a
  download without anyone re-hashing 400 GB. `--staged` folds in the checksums of archives about to
  be uploaded; regenerating after the upload and diffing proves the Hub holds exactly what was
  staged. It refuses to write an empty file, which would verify nothing and read as a pass.

- **Unrecognized `complete-search` flags are now rejected.** `complete-search` parses `sys.argv` by
  hand, so a flag nothing reads was simply ignored: `--prune-pams` instead of `--prune-pam` ran a
  genome-wide search with pruning quietly off, finished without complaint, and was detectable only by
  diffing the results against a dedicated index. Unknown flags now fail up front, with a suggestion
  for a near miss.
- **Settings: configurable number of concurrent search jobs.** The web app ran at most 2 searches at
  once (a hard-coded legacy limit); this is now adjustable from the Settings page (or the
  `CRISPRME_MAX_CONCURRENT_JOBS` environment variable), clamped to `[1, 16]`. The **default stays a
  conservative 2** — each genome-wide search is memory-heavy (~tens of GB of RAM) — so a laptop / the
  public server is never accidentally oversubscribed; raise it only on a machine with enough RAM and
  cores. Extra jobs beyond the limit simply queue. Applies on the next app restart. When
  `CRISPRME_MAX_CONCURRENT_JOBS` is set it takes precedence: the Settings card says so and its input
  is disabled, and a value that is not a whole number from 1 to 16 is ignored or clamped, with a
  one-line warning at startup.
- **Settings: Hugging Face download progress bar.** Index/annotation downloads now show percent,
  size, ETA, and phase in the Settings page, backed by an env-gated progress file
  (`CRISPRME_PROGRESS_FILE`); the progress path is a no-op when the variable is unset, so the CLI
  happy path is byte-for-byte unchanged.

### Changed
- **CRISPRitz pin moved to v2.9.1**, the tag carrying `--prune-pam` plus eight engine fixes: a
  genomic `N` in a stored PAM/flank region no longer drops the target; a non-`.bin` sidecar in the
  index directory no longer crashes the search; a PAM/flank width mismatch is refused up front
  instead of mis-parsing the leaf stream, and `buildTST`/`searchTST`/`searchBruteForce` now have
  their exit status checked, so a crashed search can no longer be reported as an empty result; a
  multiallelic site no longer compounds its indel alleles into one fake contig; the compiled C++
  `enricher` no longer drops a SNP wherever an indel shares its position; `enricher.py` keeps the
  real SNP at a multiallelic site whose ALTs include `*` (it left the base unchanged; the C++ port
  already coded it); and `add-variants` now fails when an enricher fails, instead of exiting 0 with
  that chromosome's enriched genome missing, which the index build would then have indexed as it
  stood (the Python enricher, moreover, no longer needs a `python` on PATH, and no longer races
  its sibling workers to create their shared output folder, which could lose a chromosome); and
  `buildTST` skips a contig too short to hold one target, such as the header-only FASTA the
  enricher writes for a chromosome with no indel, instead of crashing on it.
  Applies to the Docker build (`crispritz_ref`) and `install_from_source.sh`
  (`CRISPRITZ_REF`). The `*` case needs a multiallelic VCF, and every CRISPRme-supplied VCF is
  split, so it moves no published index.
  **No published index is affected** by the C++ enricher fix: the image built only the three search
  binaries, so enrichment ran through `enricher.py`, confirmed by direct measurement on the
  artifacts (chr22:10684299 carries a SNP and an indel as separate records; every enriched
  `chr22.enriched.fa` on the build host holds the IUPAC `M` there, not a plain `C`). A **source**
  install is the exception -- `install_from_source.sh` has installed the compiled binary for some
  time, so an enriched genome produced that way before CRISPRitz 2.9.1 is missing the SNP at every
  position that also carries an indel, and should be re-enriched.
- **The compiled `enricher` is now installed in the image, making it the default enrichment path.**
  It is built beside the three search binaries and copied to both places `crispritz.py`'s
  `_enricher_command()` checks, matching what `install_from_source.sh` already did and removing a
  silent divergence between the two install routes. Measured **1.82x** faster than `enricher.py`
  through the real `add-variants` CLI on the FULL chr22 of HGDP + 1000 Genomes -- 1,998,332 records,
  4,091 samples, 2.7 GB of output: 1,776 s vs 974 s -- on the stage that dominates a from-scratch
  genome-wide build. Gated on CRISPRitz 2.9.1's byte-identity fixes, verified at that same scale:
  all four enricher outputs match `enricher.py` hash for hash, as do the 17 synthetic cases of
  `verify_enricher.sh`, the SNP+indel combinations among them. The same holds off Linux/amd64:
  natively on macOS arm64 (Apple clang, libc++) the compiled port reproduces those full-chr22
  hashes exactly, and so does `enricher.py` under pandas 3 (the recipe leaves pandas unpinned);
  inside the arm64 image, the binary it ships reproduces them too.
- **When a dedicated index for the PAM is installed beside a PAM-less one, the dedicated index is
  now searched.** The index scans settle every choice by rule: the requested PAM first, then the
  smallest sufficient bulge budget, then the narrowest PAM-less flank. 2.8.0 kept the smallest
  budget first, so with `NRG_4` and `NNN_3` installed it searched `NNN_3`; 2.8.1 searches `NRG_4`,
  which gives the same report faster (genome-wide, the search stage took 44-51 min on a dedicated
  NRG index against 2 h 59 min on `NNN_3` under `--prune-pam`). No rule depends on glob order, so
  the choice cannot depend on the shell's locale either, which matters now that PAM-less indexes
  of two widths can sit side by side: `NNNNNNNN_3_hg38` globs before `NNN_3_hg38` under
  `LC_ALL=C` and after it under `en_US.UTF-8`. `PostProcess/test_index_resolution.py` checks each
  rule with the wrong candidate globbing first, and runs the whole shipped script against a stub
  engine to check that each scan hands the search the index it picked.

- **`build-index-only` makes a missing `--samplesID` impossible to miss.** Building a genotyped
  variant index without `--samplesID` (or with no usable entries) skips the dict-less Tier-0
  registry, so the index lacks fast post-analysis, the report's variant-count panel and the web's
  data-type labels, and cannot be published. The one-line notice is now a red banner on stderr
  that names the consequence and the fix (re-run with `--samplesID`). Behaviour is unchanged.

- **Install docs: Bioconda and macOS.** The README and website said Bioconda installs the old 2.1.x
  line; it has served CRISPRme+ (Python 3.11, Linux) since 2.7.9. The channel setup now uses
  `conda config` (Mamba 2 rejects `mamba config --add`), and the install names its version
  (`crisprme=2.8.1 python=3.11`), since an unpinned solve can silently fall back to a years-old
  release; `scorer-env create` and `doctor` are part of it. macOS is Docker-only: Bioconda has no
  macOS CRISPRitz, and the source build needs GCC (`buildTST` uses GCC's parallel-mode STL).

### Fixed
- **One indel off-target could be reported on every chromosome, each time with that chromosome's real
  rsID, allele frequency and carrier samples.** Every CRISPRitz targets file carries a `#` header, and
  `post_analisi_indel.sh` copies the merged indel file's first line as the header of each
  per-chromosome subset — it is the only stage that does so, which is why only the indel leg was
  affected. Two filters read that header as a target and dropped it:
  `pam_filter.py` guarded headers with a column count, but the CRISPRitz header has all 10 columns, so
  its `DNA` column was taken for a PAM region and rejected; and the `--max-total-edits` cap compared
  `$10` (`"Total"`) against the numeric cap as strings. With the header gone, `head -1` returned the
  first *data* row, and that single row was prepended to all ~23 subsets. Nothing rejected it, because
  every chromosome's fake genome starts at coordinate 0 with near-uniform windows, so the indel lookup
  (keyed on the fake position alone) resolved it against whichever chromosome was running and overwrote
  the `Chromosome` column with that one. One raw search row therefore became one correct plus N phantom
  report rows, each indistinguishable from a genuine annotated variant off-target. Measured genome-wide
  on `ACTGAAATCTGTAAGCAGGC` at NRG 4+1+1 over HGDP+1kGP: 8 phantom rows on a dedicated NRG index and 11
  on a pamless index, against an identical 13,253-row justified core; because the injected row is
  whatever lands at line 1, the count varied with the index. `pam_filter.py` and the edits cap now keep
  `#` lines, `post_analisi_indel.sh` prepends only an actual header, and the indel post-analysis rejects
  and counts rows carrying another contig's label. The injected row also moved one legitimate row's
  CRISPR-Bulge score by one unit in the last place: scoring runs in fixed-size batches, so an extra row
  changes a batch's contents and with it the order of the reductions inside it. That score now agrees
  across indices. The SNP leg was never affected. Covered by
  `PostProcess/test_indel_header_injection.py`.
- **A chromosome with SNPs but no indels no longer aborts the whole search** (a regression in this
  release's CRISPRitz pin, caught before shipping). The enricher builds each chromosome's fake-indel
  genome from the indels some sample carries, so a chromosome with none gets an empty genome;
  `buildTST` cannot index it and leaves that `_INDELS` sub-index without a `.bin`, and CRISPRitz 2.9.1
  -- rightly, for a broken index -- refuses to search such a directory. One SNP-only chromosome
  therefore failed the entire run, reported as an SNP post-analysis failure, where 2.8.0 had silently
  returned the correct answer: no indel off-target there. Reproduced through the real CLI on chr21
  (23,235 SNPs, no indel) + chr22 (2,891 indels): 2.8.0 exits 0 with 284 rows, 2.8.1 without the fix
  exits 1 with none. `pool_index_indels.py` no longer hands an empty fake genome to `buildTST`, and
  `pool_search_indels.py` skips a chromosome whose indel log records no indel, tolerating the missing
  part file for that chromosome only. With the fix the same run exits 0 with a report byte-identical
  to 2.8.0's. An index that is empty although its log lists indels is still a failure, and still
  reported as one. Covered by `PostProcess/test_indel_empty_chromosome.py`.
- **`NRR_3_hg38+hg38_mega` is republished: the 2026-10-27 archive found no indel off-target at
  all.** Measured in the archive itself (downloaded; SHA-256 equal to the Hub's): its `_INDELS` holds
  24 sub-indexes and **0** `.bin` files, and its `log_indels_hg38_mega` logs are header-only. It was
  indexed from a mega fake-indel genome that predates the sites-only fix (#191) -- every chromosome's
  FASTA a bare header -- by a script that bypassed `seq_script/merge_panels/mega_build_indels.sh`;
  `buildTST` crashed on each empty genome and nothing checked. Every NRR mega search therefore
  reported reference and SNP off-targets only, and on 2.8.0 said nothing. It also damaged
  `NRG_3_hg38+hg38_mega` in a shared working directory: an install replaces each `Dictionaries/`
  directory, so it swapped NRG mega's 419 MB of indel logs for the empty ones, and NRG mega then lost
  its indel off-targets too. The replacement's `_INDELS` is rebuilt with CRISPRitz 2.9.1 from the
  2026-10-10 mega kit genome -- which reproduces the published NRG mega `_INDELS` byte-for-byte, 32 of
  32 files -- giving 67 `.bin` across 24 sub-indexes, and its `Dictionaries/` + `samplesIDs/` are
  taken from the published NRG mega archive, so both mega archives now install identical
  dictionaries. Its SNP index is unchanged. Measured genome-wide (sites-only, 4+1+1): the published
  archive gives **0** indel off-target rows on 2.8.0, the replacement 2,081,489 on 2.8.1; and on
  each of the 24 indel sub-indexes, the replacement's NRR targets filtered to NRG equal the
  published NRG mega's exactly (1,131,494). **If you downloaded `NRR_3_hg38+hg38_mega` before
  2.8.1,** remove `genome_library/NRR_3_hg38+hg38_mega{,_INDELS}` and download it again; that also
  restores NRG mega's indel logs where the two share a working directory. Every other published
  variant archive was checked the same way and is complete. A new release gate (T-40) checks this in
  the archive for every variant index.
- **`NRR_3_hg38+hg38_HPRC` is republished so that it defaults to per-sample analysis, like
  `NRG_3_hg38+hg38_HPRC`.** Its index is complete, but its registry was built before the
  `data_type`/`phased` stamps were added, and `complete-search` picks the default mode from that stamp
  (per-sample only on `genotyped-phased`). An NRR HPRC search therefore ran population-level --
  no named carriers, no CONFIRMED cis phasing -- on a phased panel, without saying so, and installing
  it next to NRG HPRC replaced NRG HPRC's stamped registry and changed its default too. The
  replacement changes only `Dictionaries/`, taken from the published NRG HPRC archive; the registry
  data itself is unchanged. Passing `--per-sample` explicitly was never affected.
- **Both `NNN_3` archives are republished, rebuilt with the CRISPRitz 2.9.1 builder.** The archives
  uploaded on 2026-10-25 predate the genomic-`N` flank fix, so a search that used the unconstrained
  `20bp-NNN-NO-PAM` as its *search* PAM missed sites whose PAM bases are assembly `N`. Any
  constrained PAM was already exact: the downloaded `NNN_3_hg38+hg38_HGDP1kGP`, searched with NRG
  `--prune-pam`, reproduced the dedicated NRG report byte-for-byte genome-wide. The rebuild uses the
  old build's own inputs and changes the index only where the fix predicts: 1,222 of the reference
  index's 1,628 bins (+145,540 B) and 1,221 of the variant index's (+106,849 B), on exactly the 66
  contigs that contain `N` and on none of the other 389; and the whole indel index, whose 23 contigs
  are cut by `N` every ~56 bases by construction (222 bins become 245, +2.1 GB). Built independently
  on two nodes, the reference index is byte-identical. Searched genome-wide by 2.8.1 (NRG
  `--prune-pam`, 4 mm + 1 DNA + 1 RNA bulge, per-sample), the rebuilt pair matches the gate's report
  in every column but `Annotation_COSMIC` (189,950 rows, 0 unjustified indel rows), which differs
  only because that run used freshly downloaded annotations, where COSMIC stays off until its licence
  is attested. The released 2.8.0, searching chr22, gives byte-identical reports on the old and the
  new reference index, with NRG and with NNN. The variant archive grows from 108.4 GB to 109.5 GB.
  Nothing to do unless you search `NNN_3` with the unconstrained PAM; then remove
  `genome_library/NNN_3_hg38`, `genome_library/NNN_3_hg38+hg38_HGDP1kGP` and its `_INDELS`, and
  download again (`download` skips an index that is already installed).
- **A PAM of any length other than 3 could be handed the 3-nt PAM-less index.** The three scans
  that pick an index (the reference scan, the variant scan and the precomputed-variant detection)
  globbed the requested PAM and then a literal `NNN_*`, whatever the PAM's length, and kept the
  smallest sufficient bulge budget. So with `NNN_3` installed (published since 2026-10-25) and no
  dedicated index winning, a search with SaCas9's 6-nt `NNGRRT`, SpCas9-NG's `NG` or Cas12a's
  `TTTV` ran on an index whose stored flank cannot serve it. Reproduced on the released 2.8.0 with a
  synthetic genome holding three planted SaCas9 sites: with `NNN_2` installed the run reported
  `Reference index present (NNN_2_synth)` and failed with `Error on reading the PAM, input char is
  not an IUPAC nucleotide`, where an empty library builds `NNGRRT_2` and finds all three. 2.8.1,
  in the same layout, builds `NNGRRT_2` and returns the same three rows. The scans now share one
  resolver that takes a PAM-less stand-in only when its name is all `N` at the PAM's width (or,
  under `--prune-pam` and for a 3′ PAM, wider), so a dedicated PAM that merely begins with `N`s,
  such as the shipped `NNNRRT`, is never taken for one either. Covered by
  `PostProcess/test_index_resolution.py`.
- **`complete-test` passed a flag that did not exist.** The brute-force benchmark command included
  `--bmax`, which `complete-search` has never had — it recomputes bMax from the per-type budgets — so
  the argument was discarded unread. Harmless until unknown flags started being rejected, which then
  broke every benchmark. The flag allowlist was checked against the help text and against what the
  parser reads, but not against the commands this repo builds, so a dead flag in a caller could only
  surface by running the benchmarks; that gap is now a static test.
- **An index still being built could be selected as if it were finished.** A build takes
  `<dir>.lock`, and `index-genome` then creates the index directory and writes one `.bin` per
  chromosome, so for the whole build a half-written index is indistinguishable from a complete one
  by name. A run starting in that
  window printed `Reference index present` and searched a **partial** index — fewer off-targets, no
  error, exit 0 — and did so even when a complete, sufficient index sat right beside it. The scan now
  skips any candidate whose build lock is held, and any directory holding no `.bin` at all, and the
  run waits for that build instead, under `--index-path` as well. A symlinked index directory, as a
  library shared between working directories uses, still counts: the check follows the link. Index
  choice on a finished library is unchanged.
- **A build lock left by a run that died stopped every later run that needed the same build.** A run
  that builds a shared output (an enriched genome, a genome index, an INDELS index) holds
  `<output>.lock` while it writes it, and every other run that needs that output waits for the lock
  to go. In 2.8.0 that wait had no end: a run killed mid-build (SIGKILL, an OOM kill, a container
  stop, a reboot) left its lock behind, and a later run needing the same output waited on it
  forever, or, for an index, took the partial one for finished (above). Each lock now holds an owner
  record (host, PID namespace, PID, start time), refreshed every minute while the build runs. When
  the owner is gone on this host (it died, is a zombie, or its PID now belongs to another process),
  the next run takes the lock over and builds the output again from scratch; from another host or
  container, where the PID cannot be checked, it does so once the record has gone 15 minutes without
  a refresh. With `--index-path`, which never builds a reference index, the run stops instead with
  an ERROR that names the partial index and the command that clears it. A build stopped on this host
  (Ctrl-Z, `kill -STOP`, a SLURM suspend) is waited for however long. A run ended by Ctrl-C, a
  SIGTERM (`scancel`, a SLURM time limit) or a hangup leaves its lock for that takeover rather than
  releasing it, since a released lock would hand on a half-written output as finished. A build
  stopped for 15 minutes or more where no other run can check it (in another container or on
  another host) can lose its lock; resumed, it ends with an ERROR once its build step returns and
  uses none of what it built, so a library shared across hosts or containers must not have a build
  suspended that long. While it waits, a run notes in `log_verbose.txt` who holds the lock and what
  to do. A lock with no record, as 2.8.0 takes, is never broken by age, since a 2.8.0 run may still
  be building under it: the run waits and says every five minutes how to clear it. Covered by
  `PostProcess/test_index_resolution.py`.
- **The RAM thread guard sized a search served by a PAM-less index as if it were a dedicated one.**
  It chose its memory model from the requested PAM, so an NRG bulge search served by an all-N index
  (`NNN_3`, or a wider flank under `--prune-pam`) was sized as constrained (4 GB base, 0.8 GB per
  thread) instead of PAM-less (12 GB base, 2.2 GB per thread), and could be given more threads than
  the host's memory holds, and be OOM-killed. A bulge search is now sized by the index it actually
  searches; a search with no bulges, which scans the genome FASTA and opens no index, is sized as
  before. Only the thread count changes, never the results.

- **A downloaded variant index warned that its genotype store was missing right after installing
  it.** The "no `genotypes_<vcf>.tar.gz` companion ... per-sample Samples degraded" note was written
  as the `else:` of the shared-reference-index `try`/`except`, so it fired whenever *that* download
  succeeded — an unrelated event. It was wrong in both directions: false on a complete download, and
  silent when the companion really was absent but the reference index was already installed, so its
  block never ran. The note now follows the companion fetch's own result, and distinguishes
  "genuinely absent" from "the fetch errored" — the latter already reports itself, and calling it
  absent would be wrong.

- **A search that legitimately found nothing was reported as a failure.** Three places sized a
  CRISPRitz targets file with `-s` or `wc -l` to decide whether the search had found anything. Those
  are equivalent to counting data rows only while the file has no header — which was true only
  because the filter bug above was stripping it. With the header correctly preserved, a zero-hit
  search left a header-only file, the zero-hit short-circuit saw a non-empty file, and the run
  aborted with "the search found 1 raw off-target(s) but none survived post-analysis". The two bugs
  had been cancelling. All three sites now count non-header rows through one helper, which also
  handles a file that has no header rather than blind-skipping its first line. Caught by the Cas12a
  brute-force benchmark, whose chr22 result is genuinely empty.

- **`--prune-pam` could run with pruning silently off.** `crispritz.py` builds its search argv
  explicitly from the options it recognizes, so a CRISPRitz build without `--prune-pam` does not
  reject the flag — it drops it, and the search runs unpruned over the whole PAM-less index. The
  post-search PAM filter still enforces the requested PAM, so the report came out right while the run
  read an index many times larger than it needed. `complete-search` now checks that the installed
  CRISPRitz can prune and fails with an explanation if not, rather than assuming it from the pin,
  which is overridable through `CRISPRITZ_REF` or an engine already on `PATH`.

- **The report manifest could stamp a hardcoded version.** The run manifest in the shareable report
  fell back to a literal `"2.8.0"` when a run recorded no version. It now resolves the real version,
  and reports `unknown` rather than a plausible-but-wrong number — this manifest is the IND
  traceability record.

- **The report manifest named the wrong search-engine version.** `crispritz_note` carried a hardcoded
  `(CRISPRitz v2.8.3)`, a literal that outlived three CRISPRitz pins, so every shipped manifest named
  an engine the install had never run. The manifest now carries a `crispritz_version` field resolved
  from the `crispritz.py` on `PATH` — the engine that actually ran — and `null` when it cannot be
  determined, with no version in the prose. The test that should have caught this was asserting the
  stale literal; it now asserts the property instead.

- **History page no longer crashes when a job is queued.** With every search slot busy, a queued
  job's `log.txt` exists (a queue note is written to it) but has no `Job\tStart` line yet;
  `read_job_info()` raised `StopIteration`, which aborted the entire History-table build — so running
  and completed jobs could not be opened until the queued job started. `read_job_info()` now returns
  `None` for a not-yet-started or missing log (the caller falls back to the `.Params.txt` mtime, the
  same stand-in the assembly path uses), and the per-job summary loop is guarded so a single not-ready
  or malformed job directory can never take down the whole table. (Reported by M. Tognon.)
- **Index/annotation download install is atomic and always structurally verified.** A partial or
  failed download never touches the final install location, and a failure surfaces a copy-paste log
  plus a Retry action instead of a silent or corrupt state. Two genome/VCF download behaviours change
  from 2.8.0, which wrote straight into `Genomes/<ref>` and `VCFs/<ds>`: a folder that is a symbolic
  link (for example to a shared copy) is replaced by the downloaded copy, with a NOTE, and the shared
  copy is left untouched; a folder that is a mount point (including a bind mount) is refused before
  the download starts, saying what to do instead, so `download --what all` stops at that component.
- **Pipeline failure summaries are binary-safe.** `summarize_pipeline_failure()` now reads `log.txt`
  and `log_error.txt` with `errors="replace"`, so a single non-UTF-8 byte in a tool's captured
  stderr no longer raises `UnicodeDecodeError` and masks the real error.

- **`scorer-env create` failed whenever conda was the env manager, so a Bioconda or source install
  could not run a search.** The CRISPR-Bulge env was created with
  `conda create -y -n cbulge python=3.10 -c conda-forge tensorflow-cpu=2.13 …`. conda parses the
  package specs as one positional argument, so the specs after `-c` were rejected as "unrecognized
  arguments" (conda 25.3 and 26.9 alike; mamba and micromamba accept either order). conda is the
  manager whenever the env was activated with `conda activate` (it sets `CONDA_EXE`), which covers a
  Bioconda install and `install_from_source.sh`: the scorer env was never created, and every search
  then refused to start (or ran CFD-only under `CRISPRME_SKIP_SCORER_ENV=1`). The command now puts
  every option before the package specs. The Docker image builds the env with micromamba and was not
  affected. Found by installing `crisprme` 2.8.0 from Bioconda the way a new user would (fresh
  Miniforge, conda 26.9.4); with the fix, `scorer-env create`, its compute self-test and real
  searches all pass.

### Security
- **The web Settings page could be pointed outside `Settings/` and made to re-run a recorded
  command.** Retry, Show-full-log and the job-progress refresh built a path from the job id the
  browser sends back, so an id such as `../../x` reached outside `Settings/`, and Retry re-ran the
  `cmd.json` it found there through the shell. The page now accepts only ids of the exact form it
  issues and only directories inside `Settings/`, and Retry re-runs only a job its log shows as
  failed, never a maintainer-only publish while maintainer mode is off.
- **Genome archives are checked before they are extracted.** A genome `.tar.gz` uploaded in
  Settings, or fetched by `download --what genome --source url` (the Settings "Direct URL" source),
  was extracted with no check on its member names, so a crafted archive could write anywhere the app
  can (as root in the Docker image), including a job folder that Retry would then run. Every member
  is now checked first: an absolute or `..` path, a link, a device or a FIFO refuses the whole
  archive and nothing is written. Real genome archives (UCSC hg38 and mm39 chromFa) are unaffected.

## [2.8.0] - 2026-10-23

### Added
- **Reference-sequence score columns + aligned spacer in the curated tables.** The top-1000 tables and
  every curated download TSV gain `CFD_ref` and `CRISPR-Bulge_ref` — the score of the REFERENCE sequence
  at each site, shown beside the ALT/worst-case score so a reviewer sees how much a variant moves each
  score — and an `Aligned_spacer+PAM` column placed immediately before `Aligned_protospacer+PAM`, so the
  guide and protospacer read as a pairwise alignment (bulge columns as `-`). `CRISPR-Bulge_ref` is gated
  like the CRISPR-Bulge column (shown only when CRISPR-Bulge was computed).
- **CRISPR-Bulge sensitive tier (benchmark 95%-recall operating point) in the Section-4 threshold
  table + a `crispr_bulge_gt_0.tsv` download.** The CRISPR-Bulge candidate-count table shows
  `≥ 0.5 / ≥ 0.2 / ≥ 1.9e-4 (95% recall)` — the orphan `≥ 0.1` tier is replaced by the **sensitive
  operating point** the recommended panel keys on: the threshold that recovers ~95% of validated
  off-targets on the 40-guide held-out benchmark (point 1.9e-4; guide-bootstrap 95% CI [5.8e-6, 3.1e-3];
  selectable via `CRISPRME_CB_RECALL_TARGET={90|95|99}`, single source of truth in `cb_operating_point`).
  The matching subset is bundled as a download, consistent with the risk-card sensitive row. The file
  key stays `crispr_bulge_gt_0.tsv` (stable download/web identifier).
- **CRISPR-Bulge scores are written at full precision (`{:.6g}`) instead of 3 decimals.** The sensitive
  default (1.9e-4) is below the old 3-decimal display floor (5e-4), which would have silently rounded
  ~7% of validated off-targets to `0.000` and dropped them from the sensitive tier (recall 0.93→1.00 on
  the benchmark). Full precision across the population/per-sample and indel scorers makes the operating
  point representable. CFD display is unchanged (3dp; its tiers are ≥0.2/≥0.5, far above the floor).

### Changed
- **CFD now penalizes bulges by position (Doench et al., 2016 gap matrix), instead of ignoring them.**
  CFD's score tables include position-weighted DNA-bulge (`r-:dX`) and RNA-bulge (`rX:d-`) factors, but
  the DNA-bulge factors were never applied: the scorer realigned a DNA-bulge target by trimming the 5'
  end, which collapsed the gap to a table position that has no entry, so every DNA bulge — and every
  bulge at the two undefined terminal columns — contributed a factor of 1.0 (no penalty). A bulge
  immediately adjacent to the PAM (maximally disruptive) therefore scored the same as no bulge
  (e.g. a chr19 site scored CFD 0.87 with a PAM-proximal DNA bulge, while CRISPR-Bulge scored it 0).
  CFD now scores the alignment as displayed — the extra genomic base of a DNA bulge is removed so the
  20 guide bases keep their true Cas9 positions (this also fixes a mispairing in the old 5'-trim), and
  the position-weighted gap factor is applied at the bulge's real position, with the two undefined
  terminal columns extrapolated (lenient at the PAM-distal 5' tip, abolished at the PAM-proximal
  junction). Mismatch-only sites and in-range RNA bulges are unchanged (bitwise); only previously
  unpenalized bulges move. The same position-aware scoring is applied across every CFD path (the
  population and per-sample analyzers, the indel analyzer, the two-pass worst-case screen, the
  standalone target scorer, and the factorized fast-mode scorer). A ≥2-bp bulge (opt-in `--bDNA`/
  `--bRNA ≥ 2`, outside the single-bulge calibration domain) applies the position factor per looped-out
  base; see METHODS §3. The deployment calibration (isotonic P(real), AUPR, precision/recall) was refit
  on the re-scored 2.94M-site benchmark and the report's PR curves regenerated, including a new
  bulge-stratified panel embedded in the risk-card modal. Measured effect: a per-site correctness fix,
  not a bulge-ranking gain — overall AUPR unchanged (0.168; the benchmark is mismatch-dominated) and
  CFD's bulge-only AP stays at the floor (≈0.01 vs CRISPR-Bulge ≈0.25), but removing the PAM-proximal
  over-call drops ~9,000 spurious high-CFD bulge calls at CFD≥0.2 (precision 0.006→0.007; at CFD≥0.5
  0.05→0.085). It reinforces the recommended-panel logic (CRISPR-Bulge ranks bulges; CFD rescues a
  single-bulge site only at high confidence) rather than changing it.
- **Context-aware default analysis mode: per-sample on a genotyped+phased index.** With no mode flag,
  `complete-search` now auto-selects **per-sample** (CONFIRMED cis phasing + named carrier samples +
  exact joint AF) when the variant index is genotyped **and phased** (e.g. the default HGDP+1000G
  panel), and **population-level** otherwise (sites-only / unphased / hybrid / reference-only). The
  index `data_type` is read from the Phase-2.1 registry manifest. New **`--population-level`** flag
  forces the fast worst-case screen (opt out of per-sample on a phased index — recommended on a very
  dense / high-edit / genome-wide search); `--per-sample` still forces resolution; the two are mutually
  exclusive. The web Analysis-mode selector mirrors this (defaults to per-sample on a phased panel). A
  one-line notice prints the auto-selected mode and how to override. Rationale: phased panels are built
  for exactly this resolution, so the richest, most honest result is the default where it is tractable.
- **The recommended validation panel (report Section 4) selects on a principled sensitive operating
  point (benchmark 95%-recall), not an ad-hoc `> 0`.** The panel hard-includes
  `(CRISPR-Bulge ≥ 1.9e-4) ∪ (mm+b ≤ 2) ∪ (CFD ≥ 0.2 at 0 bulges) ∪ (CFD ≥ 0.5 at 1 bulge) ∪
  (clean 0-mismatch bulge) ∪ must-includes` — making the shipped `panel_top100.tsv` *be* the recommended
  default panel rather than a CFD-centric subset. The CRISPR-Bulge sensitive threshold is the **95%-recall
  operating point** on the 40-guide held-out benchmark (CB cannot rank bulges: bulge=1 AP ≈ 0.02 vs ≈ 0.25,
  so it is the primary axis); the CFD rescue is sensitive off-bulge (≥ 0.2) but requires high confidence
  (≥ 0.5) on single-bulge sites where CFD is an extrapolation; ≥2-bulge CFD-only sites are kept-and-demoted
  (uncalibrated for both scores), never dropped. One consistent threshold message now runs across the risk
  card, the "how to read this" modal, the Section-4 methods note, the score glossary, and METHODS §3.5.
  Why a principled cutoff and not literal `> 0`: CRISPR-Bulge assigns a tiny non-zero score to ~every site
  (99.9% of benchmark negatives fall below 1.9e-4), so `> 0` would list essentially the whole genome.
  95% recall is the sensible trade — for the ACTGA/NRG worked example (189,957 candidate loci) it yields
  162 off-target validation loci (90% → 58; 99% → 1,158; the shipped crispr_bulge_gt_0.tsv tier); pushing from 95% to 99% recall is ~7× the sites because the
  last few percent of true off-targets sit in the score's noise floor. Framed as a risk/validation-size optimization, not a
  recall guarantee; recall above ~95% is reachable (selectable) but the list grows steeply.

### Fixed
- **Reference-only search on a read-only install no longer dies on the `vuoto.txt` placeholder.**
  The serial annotation-index pre-build (the flock-race fix) called `load_annotation_bed()` directly
  on the annotation file, bypassing the empty-annotation guard; for a reference-only search that file
  is the empty `vuoto.txt` placeholder, which `pysam.tabix_index()` bgzips in place and removes —
  consuming it, so the later gene-annotation step's `cp vuoto.txt` failed and aborted the search (and
  the write fails outright on a truly read-only install / Apptainer SIF). `load_annotation_bed` now
  returns `None` for `vuoto.txt` (never indexes the placeholder) and the pre-build skips it
  (`test_annotation.test_vuoto_placeholder_never_indexed_or_consumed`; release trap T-38).
- **The "Observed" cell no longer reports a misleading raw carrier count for population-level,
  multi-variant off-targets.** A dense-window site that stacks N variants into one worst-case
  protospacer previously showed e.g. `128 carriers`, which reads as "128 people carry this
  off-target" — but in the population-level (unphased) default that number is the **union** of
  samples carrying ≥1 of the N variants, and the stacked combination may exist on no single
  haplotype (a possible chimera; the row's own low joint AF gives it away). It now reads
  `128 carriers (>=1 of 12 vars, poss. chimera)`. A single-variant site still reads `N carriers`, and
  under `--per-sample` (phased, exact cis set) it reads `N cis carriers`. Applies to the in-report
  table and the download TSVs (one shared cell).
- **The per-site "Observed" modal/caption no longer calls population-level reconstructions
  "observed haplotypes."** In the population-level (unphased) default a site reached by several
  variant combinations now reads `N variant combinations (worst-case)` (modal), `N other variant
  combinations` (link), and the caption spells out that these are reconstructed from the observed
  variants and are unphased — "run `--per-sample` to resolve the true cis haplotypes and carriers."
  Under `--per-sample` the wording stays `observed haplotypes` (genuinely resolved/phased). Mirrors
  the carrier-count fix one level up.
- **`--per-sample` searches no longer fail at the annotation stage on Docker Desktop for Mac/Windows.**
  The post-search step annotates three target files in parallel, each building the same tabix index;
  `flock` is non-functional on VM-backed bind mounts (virtiofs/9p → ENOTTY), so the concurrent builds
  raced, corrupted the index read, and tripped the stderr-fatal gate — killing the whole search at the
  final stage with a misleading "file may be corrupt, re-download". Fix: pre-build the tabix index once,
  serially, before the parallel annotators (they then only read a complete index), and tolerate a
  `flock` that errors on such filesystems. Linux/cluster (working `flock`) was unaffected.
- **Transient truncated reads of the annotation no longer abort a search.** `load_annotation_bed` now
  retries the BGZF-EOF guard + index build + open with a short backoff (filesystem read-after-write lag
  on virtiofs/NFS), declaring the file corrupt only after the view fails to settle (+ 4 regression tests,
  now wired into CI).

### Documentation
- Troubleshooting entries for two CLI first-run messages a newcomer can hit: a guide needing `NNN`
  PAM-padding (23 nt), and the "output folder not empty" a validation-aborted run leaves behind.

### Removed
- Retired the standalone 1000G-2021 index (`NRG_3_hg38+hg38_1000G2021` + genotype companion) —
  superseded by the default HGDP1kGP (which includes the jointly re-phased 1000G-2021 samples);
  deleted from HuggingFace (26 GB reclaimed), rebuildable from the retained `hg38_1000G_2021`
  source. The `hg38_1000G_2021` source dataset remains (a mega + HGDP1kGP input).

## [2.7.10] - 2026-10-22

### Added
- **Top-of-report off-target risk summary (replaces the aggregated specificity score).** The report
  now opens with a per-scorer risk card: for CFD and for CRISPR-Bulge separately (never blended), the
  count of ≤1-gap off-target sites in each score band and a **calibrated expected number of true
  off-targets** with a coverage-valid 90% range, plus a qualitative note that a deployed guide's own
  off-target propensity (unknown, and not given a specific multiplier) can shift the count further.
  Estimates come from
  an isotonic score-to-probability calibration fit leave-one-guide-out on 359 experimentally confirmed
  off-targets across 46 guides (Tsai 2015, Refined TrueOT, Chen 2017, Listgarten 2018; Zenodo
  10.5281/zenodo.10902445); held-out AUPR 0.168 CFD / 0.613 CRISPR-Bulge; the range is a two-level
  guide-cluster bootstrap with a negative-binomial between-guide overdispersion term, validated
  per-held-out-guide to about 90% coverage. A CSS-only "How to read this" pop-up (no JavaScript, so the
  report stays self-contained) explains the two scorers and shows the precision-recall curves. New
  `PostProcess/offtarget_calibration.py`, `PostProcess/risk_card.py`, and `PostProcess/calibration/`
  (bundled isotonic maps + bootstrap ensemble + per-band overdispersion). The recommended validation
  panel reconciles with the card: its CFD tiers use the SAME cut points (≥ 0.5, ≥ 0.2), each showing
  the all-gaps assay count with the ≤1-gap calibrated subset (the card's domain) in parentheses, plus
  a bridge note; the wider CFD ≥ 0.05 net stays a bundled download referenced in the note. The MM/B
  matrix legend moved to a full-width band below the summary grid (flowing in columns) for a more
  compact page. Full method, provenance, and the statistician/clinician/genome-editor sign-off in
  `docs/offtarget_risk_summary_METHODS.md`.

- **Own implementation of the CFD score** (Doench et al., *Nat. Biotechnol.* 2016) in
  `PostProcess/cfd_score.py`, with the mismatch and PAM penalty tables as reviewable JSON
  (`cfd_mismatch_scores.json`, 392 keys; `cfd_pam_scores.json`, 16 keys) and a committed generator
  (`scripts/generate_cfd_tables.py`, with a `--check` table-validation mode). It is now **the**
  production CFD scorer: the SNP, INDEL, fast-mode and `targets-scoring` paths all score through
  `cfd_score.cfd` via the position-aware `bulge_cfd_penalty.cfd_score_bulge_aware`, and the previous
  per-module scorer (three `calc_cfd` copies + their helpers) and the binary `.pkl` penalty tables
  have been removed. Output is **bit-for-bit identical** to prior releases' CFD: a committed golden
  fixture (`cfd_golden.tsv.gz`, 6,132 cases spanning mismatches, DNA/RNA/terminal/PAM-proximal bulges,
  every gap key, IUPAC and non-canonical PAM) is reproduced exactly (raw IEEE-754) by
  `test_cfd_cleanroom_equiv`, which also checks clean-room internal consistency and the published
  Doench anchor values.

- **`crisprme.py targets-scoring` — score a supplied LIST of off-target sites** (no genome search).
  Takes Cas-OFFinder output (`--casoffinder`, schema auto-detected) or a generic TSV with an explicit
  column map (`--columns guide=..,dna=..,chrom=..,pos=..`) and annotates each site with the **same**
  CFD and CRISPR-Bulge calculations a normal CRISPRme+ report uses (the shared clean-room
  `cfd_score.cfd` via `cfd_score_bulge_aware`; CRISPR-Bulge via the `cbulge` env, `-1.0`
  when absent). With `--index-path` it adds a per-site **variant carrier overlay**: the variants
  overlapping each site's protospacer+PAM window and, per population dataset, their allele frequency +
  carrier count (read from the index's `Dictionaries/registry_<vcf>` tier). Enables scoring
  externally-nominated sites (e.g. non-canonical-PAM hits) without building a new index. New module
  `PostProcess/targets_scoring.py` + `PostProcess/test_targets_scoring.py` (wired into unit-tests CI).

- **Full-report opt-in surfaced on every entry point.** The slim report (above) can now be upgraded to
  the complete self-contained bundle from anywhere: `complete-search --report-full` makes the
  automatically-generated end-of-search report full (threaded via `CRISPRME_REPORT_FULL` → the shell's
  `generate_report.py --report-full`), the standalone `generate-report` CLI gained the `--report-full`
  flag (its `main()` now parses and threads it — previously only the Python `build_report(report_full=…)`
  API honored it), and the **web submission page** gained a *Report options → "Include the full report"*
  checkbox (present on every search tab, including Personal assembly) that sets the same env var. New
  `test_cli_main_report_full_flag_threads_to_bundle` covers the CLI→bundle path.

- **`observed_haplotypes.tsv` report companion — every in-budget haplotype per multi-haplotype site.**
  The per-cluster representative in `integrated_results` necessarily collapses a dense/multi-haplotype
  site to ONE headline row; this bundled companion lists **each** distinct observed variant haplotype
  at any site with ≥2 of them — its variant set, rsIDs, mm/bulges, **both** scores (CFD **and**
  CRISPR-Bulge), origin, AF, and named carriers — so a reviewer sees the full set a site contributes
  (serving the CFD-ranked **and** the CRISPR-Bulge-ranked views). Report-layer only (built in
  `generate_report` from the representatives + the alt-alignments file, merged per `(site, variant-set)`
  across axes); bundled in the SLIM report. New `PostProcess/observed_haplotypes_companion.py`.

### Changed
- **Performance: observed-haplotypes sidecar cuts the per-sample report render ~35%.** The report's
  clickable observed-haplotype modals were built by re-scanning the whole (30+ GB) alt-alignments
  table every render, even though only the alt-origin rows matter. The search now pre-reduces that
  table once (`make_oh_sidecar.py`, new STEP 8a) into a small `*_observed_haplotypes_raw.tsv` sidecar
  (alt-origin rows only — a strict superset of what the report keeps), and the report reads the sidecar
  instead. Measured on a 6,1,1 `--per-sample` run: 51.3M → 2.93M rows (31.5 GB → 8.9 GB), render
  ~25 → ~17 min, and `observed_haplotypes.tsv` is **byte-identical** (the report re-applies its exact
  per-row filter + clustering) with the prioritized top-1000 unchanged. A small zipped sidecar is kept
  after cleanup so standalone `generate-report` re-renders stay fast; older result dirs with no sidecar
  fall back to the full alt scan. New `test_make_oh_sidecar.py` + a byte-identity guard.

- **Recommended validation panel: domain-aware CFD hard-include with ≥2-bulge demotion.** The panel
  now hard-includes a site on CFD only within CFD's calibrated domain (CFD ≥ 0.5 **at ≤1 bulge**). A
  site that needs **≥2 bulges** and is high on CFD alone — where CFD is an uncalibrated over-call and
  CRISPR-Bulge is N/A — is **kept in the panel but *demoted*** (ranked after every reliable
  hard-include) rather than forced to the top on an untrustworthy score, and it is never dropped: a
  ≥2-bulge off-target can still be real (the benchmark contains no ≥2-bulge examples, so its risk is
  *unknown*, not zero) and only the user can decide what to assay. The count of such demoted sites is
  surfaced in the panel note even on a CFD-only run; the unconditional pure-bulge hard-include
  (0 mismatches, ≥1 bulge) is unchanged. The Section-7 CFD glossary tooltip and the in-report Methods
  note were updated to match, with regression tests pinning both the saturated-regime demote behaviour
  and the glossary wording.

- **`complete-search` now cleans up bulky intermediate files by default** (new `--keep-intermediates`
  to retain them). Once the finalized `*_integrated_results.tsv` and the report are produced, the
  pre-integration working files — the merged `bestMerge`/`altMerge` tables, the raw
  `all_results_with_alternative_alignments.tsv` (the `.zip` is kept), and the per-chromosome
  population summaries — are removed, reclaiming ~15 GB on a genome-wide run (measured on a 5+1+1
  `--per-sample` run). **Kept either way:** `integrated_results.tsv`/`.zip`, the `report.zip`, the
  alt-alignments `.zip`, the per-chromosome co-occurrence / high-variant-density sidecars (so a later
  `generate-report --report-full` can still rebuild those bundles), the per-guide/per-sample
  summaries, and the sqlite `.db` web results-page index (removing it would break browsing the job in
  the web interface). Gated by `CRISPRME_KEEP_INTERMEDIATES`.

- **`assembly-search` coordinate lifting now uses the `pyliftover` library** instead of an external
  liftOver binary. `pyliftover` replaces `ucsc-liftover` in `environment.yml`/`Dockerfile`, so the
  image no longer carries that binary and `assembly-search` works out of the box. Point-lifting
  behaviour is unchanged (per-endpoint target coordinate + strand; a point covered by more than one
  chain, or none, is left unmapped); chain files are still passed explicitly via `--chain-paternal` /
  `--chain-maternal` (nothing is fetched). Added a `THIRD_PARTY.md` attribution file for bundled
  components (CRISPR-Bulge, impg, pyliftover — all MIT).

- **The shareable report (`report.zip`) is now SLIM by default.** It bundles the human-facing layer
  — `report.html`, the curated top-1000 (both scorers), the validation panel, the per-tier CFD/
  CRISPR-Bulge/mmb subsets, the indel+SNP co-occurrence, and the run manifest/Params — which stays a
  few MB and easy to share. The bulky full-data files (the complete raw `integrated_results.tsv.gz`,
  the full `snp_snp_cooc`, `variant_created`, and `high_variant_density_regions`) are now bundled ONLY
  with the new **`--report-full`** flag (gzipped when bundled — e.g. a ~300 MB `snp_snp_cooc` → ~23 MB);
  they always remain in the result directory regardless, so nothing is lost. On a real genome-wide
  5+1+1 `--per-sample` run this took the default report from **274 MB → ~5 MB**; `--report-full`
  reproduces the complete bundle. The slim report degrades its HTML gracefully (no download links to
  the excluded files). `crisprme.py generate-report --report-full` opts into the full bundle.

- **The slim report now size-gates the observed-haplotypes companion.** The
  `observed_haplotypes.tsv` sidecar (one row per distinct observed variant haplotype at every
  multi-haplotype off-target site) is tiny on a normal run and keeps shipping in the slim bundle; on an
  extreme dense search it can grow to ~GB and would otherwise bloat the slim/shareable report. It is
  now bundled in the slim report only while it stays under a size threshold (default ~25 MB, overridable
  via `CRISPRME_OH_SLIM_MAX_MB`); over the threshold it moves to `--report-full`, exactly like the other
  bulky full-data files. The full file is always written to the result directory regardless, so nothing
  is lost. On a 6+1+1 `--per-sample` genome-wide run this takes the slim report from ~442 MB to ~20 MB.
  The clickable per-site top-1000 modals (both scorers, every haplotype for each shown site) are
  unchanged either way; only the top-level download link/callout/note and the modal's download button
  degrade (the per-site haplotype table and a short "available via `--report-full` / in the result
  directory" note remain) when the companion is not in the slim bundle.

- **Web Contacts page:** added Ann Cirincione under the Pinello Lab, and reworded "Please send any
  comments to Luca Pinello only" → "…to Luca Pinello" (dropped "only").

- **Web form hides the "Analysis mode" selector for a personal-assembly search.** Population-level vs
  per-sample resolves genotypes across a multi-sample variant *panel*; a personal-assembly search
  works on one individual's own two haplotypes (no panel to resolve), so the control is irrelevant
  there. The tab-toggle callback now hides the Analysis-mode block when the Personal-assembly tab is
  active and restores it for the reference-genome tab.

- **Web form shows an annotation note.** The per-search annotation selector was removed earlier, which
  obscured that searches are still annotated by default. The form now carries a short note that results
  are annotated with the functional regions enabled in Settings → Annotations (DHS, ENCODE SCREEN cCREs,
  GENCODE, COSMIC by default) plus the closest gene, and that tracks can be toggled — or a custom
  annotation added — there. (No new control; discoverability only.)

- **Install-walkthrough corrections (#92, Ann Cirincione).** Cluster SSH port-forwarding guidance
  (the server binds `0.0.0.0` with no auth; both tunnel forms given inline); per-OS Docker memory
  (macOS slider / Windows `.wslconfig` / Linux no cap); and a **backwards `--max-total-edits` rule** —
  the docs said it "must be ≥ mm+bDNA+bRNA," which is actually how you *disable* the cap (default is 6
  = 4mm+1DNA+1RNA; setting it below the sum prunes bulged sites, by design). Plus reconciled the
  reference-only disk figure (~25 GB), de-hardcoded a `:v2.7.9` pin, and refreshed stale
  `1000G2021` → `HGDP1kGP` references.

### Removed
- **Removed the stale "odd-length degenerate PAM + bulges can crash the search engine (issue #105)"
  runtime warning.** It fired on the *default* PAM (`NRG` is odd-length and carries the
  partially-degenerate `R`), so it cried wolf on essentially every default search — about a CRISPRitz
  C++ crash that was fixed upstream. CRISPRme+ hard-pins CRISPRitz 2.9.0 on every install surface
  (Docker `crispritz_ref=v2.9.0`, `install_from_source.sh`, Bioconda `crispritz >=2.9.0`), so no
  supported install can run a vulnerable engine. The `docs/` note is reworded to historical.

### Fixed
- **Observed-haplotypes companion no longer blows up memory on dense per-sample searches.**
  The build loaded the entire `all_results_with_alternative_alignments.tsv` (30+ GB on a
  dense run) and materialized every row, exploding to ~140 GB. It now streams the file in
  200k-row chunks and, for the slim report, keeps only rows near the top-1000 displayed
  sites (whose per-site modals are the only consumers of the haplotype map), with a
  vectorized pre-filter so the giant file is scanned, not materialized. Peak memory drops to
  tens of GB and the pass is ~2.5x faster; the full all-sites companion is built only under
  `--report-full`. Output is byte-identical to the unscoped build for the displayed sites.

- **A samplesID file with a trailing blank line (or any blank/comment/short row) no longer aborts the
  run.** (Reported independently by @ManuelTgn, PR #105, from a real 4091-sample HGDP+1kGP run.)
  The per-sample summary and specificity-radar steps (`process_summaries.py`,
  `radar_chart_dict_generator.py`) indexed samplesID columns without guarding against short rows, so an
  ordinary trailing newline in the `--samplesID` file crashed the job with `IndexError: list index out
  of range` at the "Creating images" stage — **after** the full multi-hour search had completed. The
  samplesID consumers now skip blank, `#`-comment, and short rows; `submit_job` strips blank lines when
  building the working `.sampleID.txt` (central fix for every consumer); and the pre-flight validator
  now checks the resolved samplesID association files up front and WARNs (never errors) if a *data* row
  has fewer than the expected 4 columns (`SAMPLE_ID`, `POPULATION_ID`, `SUPERPOPULATION_ID`, `SEX`),
  naming the file and an example line so a genuinely malformed samplesID is caught before the search
  rather than hours into it. +4 `test_validate_inputs` cases.

- **Indel off-targets now carry their cis-phased co-occurring SNP(s) in the `Variant` column.** When an
  indel off-target's true haplotype is an indel **plus** one or more co-occurring, cis-phased SNPs, the
  `Variant` column (and the `observed_haplotypes` companion, which reads the same field) previously
  showed the indel **only** — the co-occurring SNP(s) were recorded solely in the per-chromosome
  `*.indel_snp_cooc.tsv` sidecar, never merged into the off-target row. The INDEL post-analysis now
  appends the co-SNP forward-strand variant tokens (`chrom_pos_ref_alt`, matching the SNP-path format so
  `resultIntegrator` parses and strand-corrects them identically), with the `rsID` and `MAF` lists kept
  1:1 token-aligned, and narrows the `Samples` column to the **exact cis carriers** (individuals carrying
  the indel AND every required SNP on the same phased slot). **Gated to CONFIRMED genotype-phased
  co-occurrence only** — the gate is the co-occurrence *phase* (derived from the index's on-disk phased
  genotype stores), **not** the search mode, so the enrichment applies equally to a population-level and
  a `--per-sample` search on a genotyped-phased index, and cannot fire on a sites-only index (no
  genotype store → that row is byte-identical to before). A PUTATIVE (unphased / sites-only / mega) call
  is left sidecar-only and does not touch the row, preserving the anti-chimera contract (never assert a
  haplotype no individual is confirmed to carry). Verified end-to-end on a genome-wide 4+1+1 `--per-sample` run: 369 indel
  off-targets gained their co-SNP(s) with CFD byte-identical and the carrier roster narrowed to the
  exact cis set (e.g. chr10:105008040, carrier HG03711: `chr10_105008061_ATG_A` →
  `chr10_105008045_G_C,chr10_105008061_ATG_A`); SNP-only rows, PUTATIVE rows, and non-cis co-locations
  (indel and SNP carried by different individuals — e.g. the top-1000 rank-13 at chr18:12386282, whose
  worst-case alignment uses a SNP its indel carrier lacks) are left untouched. Because the corrected cis
  carrier count feeds the most-carriers representative tie-break, at a few dense multi-allele indel/STR
  loci a different — still real — observed haplotype becomes the cluster representative (no score
  change; this removes a pre-existing chimeric over-count where a worst-case indel alignment had used a
  cis SNP that most of its reported carriers do not carry).

- **`--per-sample` no longer over-states dense-window off-targets (the "chimera" over-call).** In a
  high-variant-density window (>`CRISPRME_IUPAC_CAP` co-located variants) the per-sample path emitted a
  synthetic **cross-sample** worst-case representative whose carrier roster was the UNION of different
  individuals' best-per-column alleles — a haplotype **no single individual carries** — which, being
  highest-CFD, became the per-cluster headline and over-stated both risk and carrier count (observed in
  the wild as a 15-SNP, "47-carrier", CFD-0.57 row at chr14:106770868 that nobody actually carries;
  the real observed haplotypes there are ≤2 variants, CFD ≈ 0). Per the per-sample contract — report an
  off-target only when ≥1 real individual carries it — that synthetic representative is no longer
  emitted; the dense region stays visible via the high-variant-density BED + the reference row. **Only**
  `--per-sample` is affected; the population-level default (which *should* report the conservative
  worst case) and sites-only panels are byte-identical. *Expected when diffing a `--per-sample` run
  against 2.7.9:* on dense-window variant sites the reported `highest_CFD` representative — and hence its
  `Start_coordinate` / `CFD_score` / `CFD_ref` and carrier roster — shifts to a **real** carrier
  haplotype (the chimera is gone), so a handful of such rows differ. This is the intended effect, not a
  scoring change: `CFD_alt` and all CFD arithmetic are byte-identical to 2.7.9 (validated on a
  genome-wide run — CFD identical on 189,119/189,132 shared alignments; the differing rows are exactly
  these double-bulge dense-window sites whose chimera representative was replaced by a real-carrier one).

- **Per-cluster representative tie-break → most carriers.** Among off-targets that tie EXACTLY on
  score and severity, the representative shown in `integrated_results` now prefers the haplotype carried
  by the most individuals (applies to the CFD, fewest-mm+b, and CRISPR-Bulge passes); a no-op for the
  unique-best ~99.7% of clusters.

- **CRISPR-Bulge aggregate specificity sanitized.** The CRISPR-Bulge `-1.0` out-of-domain sentinel
  (sites with ≥2 bulges, outside the model's domain) was summed into the aggregate `100/(100+Σ)`,
  driving the CRISPR-Bulge specificity negative; it is now excluded from the sum (no-op for CFD).

- **Indel-derived off-targets now get a `PAM_creation` value (#94, Ann Cirincione).** The field was
  only ever filled on the SNP path, so a PAM created by an insertion/deletion surfaced as `NA` even
  when it was the top-scoring site. The indel analyzer now mirrors the SNP predicate (if any reference
  base in the PAM window is incompatible with the motif, the variant created the PAM; store the ALT
  PAM). Writes the existing `PAM_gen` slot in place (no column shift); skipped when the window holds an
  IUPAC-undefined character (e.g. a bulge `-`). `calc_cfd` untouched.

- **`download --what index` no-ops on a *complete* install but still heals a partial one (#93, Ann
  Cirincione).** Re-running the documented download re-fetched the full tarball (~8.86 GB). It now
  skips the fetch when the install is complete, so a repeat costs no network call — but only when
  complete: a variant index needs its `_INDELS` companion, the shared reference index, and (unless
  `--no-genotypes`) its Tier-1 genotype store. Anything missing falls through to re-fetch just the
  absent piece, preserving the documented "no variant off-targets / Samples degraded → re-run
  download" recovery (and the `--no-genotypes`-then-add-genotypes workflow).

- **No more spurious `BiopythonWarning` at CLI startup.** `dash-bootstrap-components`' sdist drops a
  `pyproject.toml` into `site-packages/`, which Biopython misreads as "importing from a source tree"
  and warns about on every `crisprme.py` invocation. The Docker image already strips that stray file;
  conda/source installs did not, so `crisprme.py --version` printed the warning. Now silenced in two
  places: `crisprme.py` filters that one message before importing Bio (covers any install), and
  `install_from_source.sh` removes the stray file (matching the Dockerfile).

- **Web dropdown options no longer overlap for long labels.** The variant-index dropdown's long entries
  overflowed their fixed row height and collided with the next entry after repeated opens. The
  `variant-dataset` dropdown now reserves room for wrapped labels (`optionHeight`), and a CSS rule lets
  option text wrap and centers it within the row.

- **`targets-scoring` DNA-bulge CFD** now matches CRISPRme's reported CFD convention for in-place
  (Cas-OFFinder / report-style) alignments. The first cut front-trimmed the first `bulge_size`
  columns (a crispritz raw-format assumption), which misaligns an in-place gap that isn't at the
  5' end. It now DELETES the gap column(s) (the looped-out DNA base) from both strands and scores
  20-vs-20, so CFD's position weights stay guide-indexed. Verified against a real CRISPRme report:
  mismatch + RNA-bulge rows reproduce the reported CFD exactly; DNA-bulge CFD is the standard
  bulged-base-excluded value (CRISPRme's internal DNA-bulge gap placement isn't recoverable from an
  in-place alignment, so it stays a documented caveat). Updated `test_dna_bulge_delete_gap`.

- **Web "Load Example" no longer prunes bulged off-targets.** The Load Example button set the
  "Maximum edits" total-edit cap to 4 while loading 4 mm + 1 DNA + 1 RNA bulge, so `--max-total-edits=4`
  silently dropped every 4-mismatch-plus-bulge site (bulges had no headroom at the full mismatch
  budget). The cap is now 6 (= mm + bDNA + bRNA, the recommended default; matches the slider's own
  default and the CLI), and two inconsistent fallback defaults (5, 4) were aligned to 6. +1 regression
  test (`test/web/test_max_edits_default.py`, wired into web-e2e).

- **Web low-memory banner no longer nags a correctly-provisioned host.** A host set to the
  recommended ~64 GB reports a few GB less inside a container (GiB-vs-GB + VM reservation), so the
  banner fired even when Docker was set to 64 GB. Added a detection margin (`_DETECT_MARGIN_GB`), so
  the warning only fires for a clearly under-provisioned host; +1 unit test.

- **Removed the stale "download the offline version" note** from the web form footer. It was a
  hosted-web-service-era leftover — now that the service is retired, a user running the local web
  interface already *is* the full/offline version, so the note only caused confusion.


## [2.7.9] - 2026-10-11

### Added
- **Web low-memory pre-flight guard.** A new-user review found the clickable-Mac path dead-ends: the
  default first Submit is a genome-wide variant-aware search (peaks ~64 GB) that silently OOM-kills
  (exit 137) a typical 16-36 GB laptop. The web form now shows a **pre-Submit warning banner** when the
  host's detectable RAM is below the variant-aware whole-genome threshold (`memory_check.memory_status`,
  threshold `VARIANT_WGS_MIN_GB=64`, honoring `CRISPRME_MIN_GB`/`CRISPRME_SKIP_MEM_CHECK`), pointing to
  the escape hatches (reference-only index, single-chromosome/region, or raise Docker memory). Silent
  when RAM is undetectable; +7 unit tests. (A dedicated exit-137 message on the load page is a planned
  fast-follow; the search log already carries the `warn_low_memory` warning.)
- **`docs/RELEASE_CHECKLIST.md`** — a living release gate: the pre-tag consistency + functional +
  mechanics checks, plus a "regression traps" table that pairs every bug we've shipped (or nearly
  shipped) with the check that catches it. Appended to whenever a release pass finds a new problem.
- **`download --what assembly --sample <id>`** — fetch an HPRC individual (both haplotypes' genome +
  chain + chromAlias into `Assemblies/<id>/`) from the command line — the CLI equivalent of the web
  Data Manager's *Fetch from HPRC*. Closes the last CLI↔web gap; `docs/SETTINGS_DATA_MANAGER.md` now
  carries a full **CLI ↔ web parity** table (every web action maps to a `crisprme.py` command, so a
  CLI-built and a web-built install are interchangeable).
- **Build a variant-aware index from a single user-provided VCF.** `build-index-only --vcf` (and the
  web **Settings → Add a VCF dataset → Build index**) now **auto-splits** a single multi-chromosome
  VCF into the per-chromosome files the enricher (`crispritz add-variants`) requires — so a normal
  user VCF "just works". Previously a single multi-chromosome VCF failed with a cryptic
  empty-chromosome error (`Genomes/<g>/.fa`). (No-op for datasets already named per-chromosome, e.g.
  the shipped indexes.)
- **SNP-only custom VCFs now build a complete index.** A variant dataset with no indels previously
  aborted the whole build at the indel-indexing step (`Indels indexing failed`, empty fake-indel
  contigs); the indels index is now **skipped gracefully** and the SNP variant index is produced
  (SNP+indel co-occurrence simply does not apply to an indel-free dataset).
- **`docs/DATA_LAYOUT.md`** — a clear reference for *where your data lives*: the single user-owned
  data folder (CLI `--path`, Docker `-v …:/DATA`, Apptainer `--bind …:/DATA`), the **shared inputs**
  reused by every search (`Genomes/`, `genome_library/`, `VCFs/`, `Dictionaries/`, `samplesIDs/`,
  `Annotations/`, `PAMs/`, `Assemblies/`, `LiftoverFiles/`) vs. **where each run's outputs go**
  (`Results/<job>/`). Linked from the README quickstart, the Docker quickstart, and the Data-Manager
  guide.

### Changed
- **Final docs polish.** Reworded the web low-memory banner to actionable, web-achievable advice
  (switch to the reference-only index, or raise Docker memory — dropped the CLI-only "single
  chromosome" suggestion that has no web control). Restructured the densest METHODS sections into
  bullets/tables (§1 raw-vs-compressed encoding, §5 population-level comparison, §6 cancer
  annotations, §7 `Observed` column, §8 domain-of-validity + worst-case scoring). Reconciled the disk
  figures to one framing everywhere (~45 GB download → ~85 GB on disk → ~100 GB free; reference-only
  ~25 GB). Standardized the example guide on the full `CTAACAGTTGCTTTTATCAC` (what the web *Load
  Example* button loads). Added an SSH port-forward note to the website Singularity/Apptainer tab for
  reaching the web UI on a headless compute node.
- **CI + test hygiene.** Removed the redundant `docker-image.yml` workflow (single-arch, no-push,
  `checkout@v3` EOL — superseded by `docker-multiarch.yml`); wired the two Dash results-page unit
  tests (`test_results_page_assembly`, `test_results_summary_and_status`) into `web-e2e`; updated the
  `test_complete_test` mock index from `1000G2021` to the shipped default `HGDP1kGP`. The brute-force
  `validate-test` gate (exact set-equality vs an exhaustive aligner, fails the build on any mismatch)
  was audited and confirmed enforced on every PR + push + weekly cron.
- **complete-test docs (§2.2.2).** Download examples now use the default `NRG_3_hg38+hg38_HGDP1kGP`;
  added a note distinguishing the modern no-arg smoke (runs against whichever index you downloaded)
  from the legacy `--vcf_dataset` build-from-raw mode that feeds the `validate-test` gate (1000G only).
- **README/website readability pass (new-user experience).** Applied the docs-review nice-to-haves
  as surgical edits: Quickstart density trimmed (deduplicated the 3× RAM/default-index repeats, index
  pitch collapsed to the comparison table); Highlights five-index run-on → concise lead + table
  pointer; the mode/index decision bullet → a tight 2-column table; the report's `Observed`-column
  cascade and the §6 commercial-licence scenarios → tables; gnomAD-converter "pseudo-individuals"
  de-duplicated; `assembly-search` added to the function list; §1.1 reframed as the stable-2.1.x
  path with the CRISPRme+ update guidance first; Docker GPU wording qualified (CPU is the default,
  GPU optionally accelerates CRISPR-Bulge); `CRISPRITZ_REF` default (v2.9.0) + `http://localhost:8080`
  documented; a few internal-consistency nits (NGG→NRG example label, "three"→"four" commands,
  report list renumbered). No facts/commands/versions changed.

### Fixed
- **Final METHODS + linked-docs review pass (multi-agent, 2 blockers + consistency fixes):**
  - **`--max-total-edits -1`** — METHODS told users to pass `-1` to keep all off-targets, but the
    shipped validator rejects negatives (the search aborts). Dropped the `-1` option and fixed the
    stale in-code comment so doc + comment + validator agree (set it to `mm+bDNA+bRNA`).
  - **Mega-indel allele frequency** — METHODS said mega indels were "not yet frequency-annotated";
    they ARE (per-dataset `AF_<dataset>` + `AF_max` via the `indel_af` sidecar, read by the report).
    Scoped the "SNPs only" statement to the Tier-0 registry.
  - **`download --what vcf` deprecated in docs** — the raw per-dataset VCFs (`vcfs/1000G/`, `vcfs/HGDP/`)
    were never hosted on HF, so the command errored. Removed the misleading examples across README,
    DOCKER_QUICKSTART, the website, the web guide, DATA_LAYOUT and the legacy setup doc; clarified that
    the standard panels ship as prebuilt **indexes** (whose bundled genotype store covers sample-level
    analysis — raw VCFs are not needed) and that `--what vcf` is for custom/self-hosted repos. The CLI
    error now points at `--what index` / `build-inputs/` / `VCFs/<dataset>/`.
  - Consistency/accuracy: METHODS intro names the full default index + composition and the shipped
    version; stale present-tense "2.5.3 ships raw" / chrX "forthcoming" anchors reframed/resolved; a
    registry-cache speedup number reconciled to the source code (326→~89 s, 3.66×); the legacy-setup
    CRISPRitz-version self-contradiction collapsed to `v2.9.0`; three stale `INPUT_FORMATS` code line
    citations corrected; `PRECOMPUTED_INDEXES` `<vcf>` token defined once; `CRISPRME_IUPAC_CAP` default
    (10) stated. Also fixed a README download example that repeated `HGDP1kGP` under a "1000G-only" label.
- **Docs/website review pass (multi-agent section-by-section review, 3 blockers closed):**
  - README default-index per-source split said `HGDP (929) + 1000 Genomes (3162)` — contradicted the
    correct `925 HGDP + 3,166 1000 Genomes = 4,091` elsewhere in the same file; corrected.
  - README report-Summary example cited `3,477 individuals for 1000G+HGDP` (the retired-hybrid
    panel size) — corrected to `4,091` for `NRG_3_hg38+hg38_HGDP1kGP`.
  - Website Conda/Mamba path told pre-2.2.0 stock-Bioconda users to run CRISPRme+-only subcommands
    (`download`, `web-interface`) the stock package lacks — now routes them to the legacy `setup`
    and points the fast download/web path at the CRISPRme+ Docker/source install.
  - Also: §2.1 directory list was missing `genome_library` / `Dictionaries` / `Assemblies` and the
    one-working-folder framing (added, links DATA_LAYOUT.md); `sampleIDs` → `samplesIDs`; website
    hero now names the shipped default; the first-launch `~45 GB` figure is reconciled with the
    `~100 GB` disk budget; version-specific "unreleased 2.7.8 / in preparation" prose made
    version-agnostic; "dict-less" download jargon relabelled "variant-aware".
- **Web default variant panel now honours the re-phased default.** The first-paint variant-dataset
  seed (`_preferred_variant`) was called without the genome, so its genotyped-**phased** preference +
  sample-count tiebreak were skipped — with two or more phased panels installed it could seed a
  non-default panel. Both call sites now pass the genome, so the default page consistently proposes
  the SHAPEIT5 re-phased `NRG_3_hg38+hg38_HGDP1kGP` as the recommended variant default (still fully
  name-free: the preference is read from each index's `data_type` manifest).
- **No surface advertises a retired index.** Swept the repo for indexes removed from HuggingFace in
  v2.7.8 and repointed the few remaining references at the live catalog: the `seq_script/merge_panels/`
  README (retirement banner + repoint to `HGDP1kGP`), the `onekg2021_build.sh` header (reworded from
  "the default index" to an *alternative* single-source panel), a `pages_utils` docstring example, and
  the `test_variant_dropdown` docstring (clarifies the two-source `1000G_HGDP` fixture is a generic
  discovery-mechanism stand-in, not a shipped index). Archived 2.5.x docs already carry retirement
  banners and are unlinked from the live site.
- **Web Settings / Data-Manager (3 fresh-user-QA findings):**
  - Clicking **Download index** on an already-installed index re-downloaded all its GB and blocked
    the single-slot job queue; it now short-circuits with *"Index '…' is already installed — nothing
    to download."*
  - An in-progress download's `.stage_<genome>` staging dir was listed as an **installed genome** (and
    offered in the build / search / delete dropdowns); dot-dirs (staging + assembly markers) are now
    excluded from the installed-genome listing.
  - Launching a **Download / Build** left the inline feedback blank (job status appeared only in the
    Activity panel); every job launch now confirms *"Started — progress is shown in the Activity
    panel."*
- **Web-submitted `assembly-search` was broken** (`ModuleNotFoundError: No module named
  'validate_inputs'` at startup, so every personal-assembly search failed immediately). Two causes,
  both fixed: (1) crisprme.py's install-path correction hard-stripped a 3-char `"bin"` suffix, so when
  the web spawned its `opt/` copy (`<prefix>/opt/crisprme/crisprme.py`) the PostProcess path was
  garbled — it now detects whether PostProcess is a direct sibling (opt/ copy or source checkout) or
  reached via `opt/` (bin copy); (2) the job was launched with `--debug`, whose cwd-relative
  PostProcess path also missed (the job's cwd is the Results output dir) — `--debug` now only switches
  to that path when it exists, and the web assembly-search command no longer passes `--debug`. Verified
  live on a real HPRC individual (both haplotypes search + `--gene_annotation` closest-gene wired).
- **`assembly-search` (and `complete-search`) aborted on a read-only install** (Apptainer/Singularity
  SIF) with `truncate: … vuoto.txt: Read-only file system` (issue #143): the results-integration step
  stripped a trailing newline from the "no empirical data" placeholder, which is the shipped read-only
  `PostProcess/vuoto.txt`. It now only truncates a real, non-empty, writable file. (Docker's writable
  layer masked this; it affected HPC SIF runs.)
- **Web-submitted `assembly-search` crashed in the combine step** with `shutil.SameFileError`: the web
  writes the guide straight into `Results/<job>_combined/.guides.txt`, then `write_combined_guides_file`
  copied the guide onto that same path. It is now a no-op when source and destination are the same
  file. With this, a web assembly-search runs **end-to-end** — both haplotype searches → hg38
  reconciliation → the combined report with per-haplotype CFD + closest-gene columns (verified live on
  HPRC HG01255: 60,745 off-targets, `Annotation_closest_gene_*` present).
- **Docs/consistency precision pass (full release-gate QA, 161 PASS / 14 WARN / 2 FAIL).** Corrected
  the HGDP+1kGP cohort split to **925 HGDP + 3,166 1000 Genomes** (the shipped BCF header-derived
  split; still 4,091) across README / METHODS / DOCKER_QUICKSTART / PRECOMPUTED_INDEXES;
  `gnomAD-converter --help` now lists **v3.1, v4.0, and v4.1 (incl. joint VCFs)** (it already
  supported v4.1); documented the `scorer-env list` action; fixed the Dockerfile CRISPRitz comment
  (→ v2.9.0), the README `download --what` narrative (added `assembly`), the DESIGN_mega "two
  predefined indexes" framing (now the published set), the TensorFlow "2.13.1" → "2.13.x per-arch"
  wording (Apple-Silicon is 2.13.0), and drift-prone code-line citations in INPUT_FORMATS.
- **`complete-test` / `setup` help clarified.** `complete-test`'s "no prebuilt index" pointer now
  recommends the default `NRG_3_hg38+hg38_HGDP1kGP` (was `1000G2021`); the `--vcf_dataset`
  (1000G | HGDP | 1000G+HGDP) help now states it is the **legacy build-from-raw-VCF test panel**
  (1000 Genomes 2019 + HGDP 2019), explicitly **not** the shipped re-phased HGDP1kGP default and
  **not** the retired `NRG_3_hg38+hg38_1000G2021_HGDP` index. `setup --help` now points to the
  faster modern `download --what all` path (setup is the legacy from-scratch FTP download).

## [2.7.8] - 2026-10-11

### Added
- **Build a new-PAM index for a published dataset.** The PAM-independent enriched-genome "kits"
  for the re-phased `NRG_3_hg38+hg38_HGDP1kGP` and sites-only `NRG_3_hg38+hg38_mega` datasets are
  now published on HuggingFace (`build-inputs/hg38_HGDP1kGP.tar.gz`, `build-inputs/hg38_mega.tar.gz`,
  each with `SHA256SUMS`). A new `docs/PRECOMPUTED_INDEXES.md` section documents downloading a kit
  and indexing **any** nuclease/PAM (NGG, NNGRRT, TTTN, …) off it, plus the full source-to-kit
  reproduction (SHAPEIT5 re-phasing + data provenance for HGDP+1kGP; the five-source merge + MAF
  filters for mega). `METHODS.md` gains a re-phasing-provenance note and a "Re-indexing for a
  different PAM" subsection; `README.md` and the data-setup guide cross-link the new section.

### Fixed
- **Web `assembly-search` could not be submitted.** The hidden genome/chain/chromAlias holders
  auto-populated from the Individual picker were `dcc.Dropdown`s, which silently drop a value that
  is not one of their options; the picker stores cwd-relative *paths*, so the submit handler read
  `None` for all six and always reported "missing fields". They are now `dcc.Store` components
  (hold arbitrary strings), so a web-submitted assembly-search launches correctly.
- **Web `assembly-search` did not add closest-gene columns.** The submit handler discarded the
  gencode companion returned by `build_active_annotation`; it now passes `--gene_annotation`
  (the bundled `gencode.protein_coding.bed`) when the built-in annotation bundle is enabled, so a
  web-submitted job gets the `Annotation_closest_gene_*` columns (parity with complete-search).
- **`synth_sites_gt.py` wrote plain gzip instead of BGZF**, so `tabix` rejected the output
  ("not BGZF") and the documented sites-only (mega) fake-indel reproduction recipe failed. It now
  writes BGZF via pysam (with a BGZF-magic-byte test assertion).
- **`seq_script/merge_panels/mega_build_indels.sh` never placed the fake-indel genome** into the
  index directory — its `strip_and_place` glob (`*fakechr<N>*`) did not match the enricher's
  `fake_hg38_mega_<chr>` directory name. Search was unaffected (it uses the index bins), but the
  PAM-independent rebuild genome was missing; the glob is fixed.
- **Migration-guide license row was misleading.** The "What changed vs stock CRISPRme" table framed
  CRISPRme 2.1.14 as merely "AGPL-3.0 (copyleft)", implying AGPL alone permitted commercial use.
  Stock CRISPRme was dual-licensed — AGPL for academic/non-commercial use **plus** a separate
  commercial license for commercial use. The row now states that both CRISPRme and CRISPRme+ require
  a separate commercial license for commercial use; only the academic-use license changed
  (AGPL-3.0 → MGBOA 1.0).

### Changed
- **Documentation + website consistency sweep.** Version strings were unified to `2.7.8` across
  `README.md`, the website (`docs/index.html`), `docs/DOCKER_QUICKSTART.md`, the web-interface
  user guide, `docs/PRECOMPUTED_INDEXES.md`, and `install_from_source.sh` (several surfaces had
  drifted to `2.6.x` / `2.7.1` / `2.5.5`). The legacy data-setup guide now points to the modern
  `crisprme.py download` flow for first-time setup.
- **Fresh-user consistency pass.** A full install→run→build→interpret walkthrough surfaced and
  fixed: the **default index** is now consistently **HGDP + 1000 Genomes** across README, the
  website, the Docker quickstart, and the web-interface guide (several spots still named the old
  1000G-2021 default); the website RAM figure is corrected to **64 GB** for a genome-wide
  variant-aware search and now surfaces the **one-line installer**; a CLI-guide example that
  dropped the PAM `N`-padding (and would prune all off-targets) is fixed; the web guide gains an
  **"Understanding your results"** reference (CONFIRMED vs PUTATIVE, carriers, key columns); and
  a `docs/FUTURE_WORK.md` captures deferred enhancements.
- **One-click installer moved into this repo.** `install.sh` / `install.ps1` (previously served
  from the legacy `pinellolab/CRISPRme` repo) now live in `docs/` here and are served at
  `https://pinellolab.github.io/crisprme-plus/install.{sh,ps1}`, so the release owns its own
  installer. Their default first-run download is switched from the old `1000G-2021` index to the
  recommended **HGDP + 1000 Genomes** default, and the stated size corrected to **~45 GB** (the old
  "~75 GB" was wrong). The one-liner URLs in the README, website, and `install/README.md` are
  updated to match. (The old URL should be left as a redirect for existing bookmarks.)

### Removed
- **Retired ~71 GB of stale/superseded precomputed indexes from HuggingFace.** Deleted the
  superseded hybrid `NRG_3_hg38+hg38_1000G2021_HGDP` (+ its genotype companion), which is replaced
  by the fully re-phased `HGDP1kGP` default; the orphaned v2.3.0 dict-less lineage
  (`NRG_3_hg38-dictless+hg38_1000G_HGDP`, `NGG_3_hg38-dictless+hg38_1000G_HGDP`,
  `genotypes_hg38_1000G_HGDP`); and the unused reference-only `NGG_3_hg38`. The hosted catalog is
  now the four production indexes (`HGDP1kGP` default, `1000G2021`, `HPRC`, `mega`) plus
  reference-only `NRG_3_hg38`; all user-facing docs were scrubbed to match (catalog entries,
  download commands, and every worked example repointed to a hosted index).

## [2.7.7] - 2026-10-10

### Added
- **Closest-gene annotation for `assembly-search`** (`--gene_annotation`, #66). A
  nearest-feature (BEDOPS `closest-features`) lookup adds the closest gene + signed
  distance to every hg38-mapped off-target of a personal diploid assembly — `Gene` +
  `Gene_distance_kb` in the curated report tables / web result page, and
  `Annotation_closest_gene_name` / `_ID` / `_distance_(kb)` in the full TSV
  (distance `0.0` = overlapping, negative = the other side; haplotype-non-mappable
  sites have no hg38 coordinate and so get no closest gene). Validated end-to-end on a
  real diploid assembly (HG01255 chr19 → 952/960 sites annotated) and in the web UI.
- **Variant-dataset retrieval recipes** under `paper/reproducibility/variant-datasets/`
  for 1000G 2019/2021, HGDP, gnomAD v4.1, TOPMed and All-of-Us (paper reproducibility,
  docs-only; #45).

### Changed
- Docs: README (`--gene_annotation` nearest-feature columns + assembly-search support +
  variant-datasets pointer) and the web-interface user guide (assembly result columns).

## [2.7.6] - 2026-10-09

Stability release candidate. One correctness hardening on top of v2.7.5; no change
to the search engine, default index, or Docker/Bioconda toolchain.

### Changed
- **`--samplesID` over-listing is now a hard ERROR (was a warning).** A `--samplesID`
  that lists samples genotyped in *none* of the provided VCFs silently inflates the
  panel AN (phantom samples are counted as hom-ref) and deflates every reported
  allele frequency (Samples / AF / MAF) — no crash, just wrong numbers (#46). Since
  there is no legitimate reason to list a sample absent from every VCF, pre-flight
  validation now refuses to run and points to the fix (a VCF-filtered samplesID, or
  the `registry_fix_an` / genotyped-panel transform). The shipped default indices
  (built with correct panels) and index-path searches are unaffected — this only
  guards fresh `--genome --vcf --samplesID` builds with a mis-listed panel.

### Fixed
- **Default 1000G samplesID over-listed its VCF panel.** The shipped
  `samplesIDs.1000G.txt` listed 3500 samples while the complete-test 1000G VCF
  (`ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019`) genotypes only 2548 —
  the 952 phantom samples were counted as hom-ref, deflating every reported AF
  ~1.37× (#46), and (with the hard ERROR above) blocked the default 1000G
  `complete-test`. The panel is now filtered to exactly the 2548 VCF-genotyped
  samples (HuggingFace data repo + in-repo fallback + `MD5SAMPLES` pin updated).
  HGDP was verified already correct (929 ↔ 929).

## [2.7.5] - 2026-10-09

Stability release: off-target scoring correctness (CFD audit), assembly-search
coordinate + report fixes, and test-suite hardening. No change to the search
engine, the default index, or the Docker/Bioconda toolchain.

### Fixed
- **CFD scorer — INDEL-path N-in-PAM divergence.** `analisi_indels_NNN.calc_cfd`
  had drifted from the canonical SNP-path `new_simple_analysis.calc_cfd`: an `N`
  in the scored 2 bp PAM kept a non-zero CFD on the INDEL path, whereas the SNP
  path zeroes a non-canonical PAM (`pam_scores.get(pam, 0.0)`, issue-#94 guard);
  it also carried a dead `if "N" == sl` branch. The INDEL twin is now byte-for-byte
  equivalent to the SNP path (an `N` genomic base is unknown → zeroing is
  consistent and conservative), pinned by a new regression test
  (`test_calc_cfd_twins_agree.py`). A full CFD correctness audit confirming the
  core is bit-identical to the Doench-2016/CRISPOR reference is documented in
  `docs/analyses/cfd_correctness_audit.md`. (#63)
- **assembly-search — hg38 off-target coordinates were off by one and lost
  strand.** Off-targets are now lifted via both end bases as stranded records,
  with orientation resolved and non-mappable sites refused rather than
  mis-placed; two edge-case crashes (contig-edge `impg`, empty lift) fixed. (#62)
- **assembly-search report — stale curated-output flags + inconsistent zip
  schema.** The module-level curation flags (`_DROP_MAF`, `_PRESENT_ANN_KINDS`,
  `_HAS_VARIANTS`, `_DROP_KINDS`) are now set on entry to *both* report builders,
  so an assembly report rendered after a complete-search can no longer inherit
  stale flags; and every table in an export zip now shares one schema (decided
  once from the pooled frame). Complete-search report output is unchanged
  (verified byte-identical). (#64)
- **Test suite — collection-order-dependent failures.** `test_bulge_cap.py` no
  longer leaks its import-time `app`/`dash`/library stubs into `sys.modules`, so a
  later-collected test's `from app import …` resolves the real module regardless
  of collection order (fixes two order-dependent `TestAnnotationLogNote` failures).

### Changed
- **assembly-search report columns.** The combined annotation column is split
  into per-annotation columns matching the complete-search report, and the new
  haplotype-context columns (`Haplotype_origin`, `hg38_orientation`,
  `hg38_lift_confidence`, `n_copies`/`copy_loci`, `hg38_end`) are surfaced;
  assembly-irrelevant (e.g. VCF-specific) columns are dropped. (#64)

## [2.7.4] - 2026-10-08

### Changed
- **Pinned CRISPRitz to v2.9.0** (`Dockerfile` + `install_from_source.sh`). CRISPRitz
  v2.9.0 is now a **pure enumeration engine** — all scoring, annotation and reporting
  were removed from CRISPRitz (they live in CRISPRme+), which drops its azimuth
  `scikit-learn==1.1.3 / numpy==1.24.4 / pandas==2.0.3 / scipy==1.10.1` pins. This has
  **no functional effect on CRISPRme+** (which only uses CRISPRitz's `index-genome` /
  `search` / `add-variants` engine — verified byte-identical, golden-search regression
  on x86_64 + arm64) and unblocks publishing CRISPRme+ on **Bioconda** (the `crisprme`
  package can now depend on `crispritz` without the old, unsolvable ML version pins).

## [2.7.3] - 2026-10-08

### Changed
- **Pinned CRISPRitz to v2.8.4** (`Dockerfile` + `install_from_source.sh`), which carries the
  CRISPRitz #44 enricher fix: on **multiallelic** records the `add-variants` enricher previously
  matched a sample to an ALT allele by substring-searching the raw genotype text and stopping at
  the first hit, so (1) a record with ≥10 ALTs could bind a sample to an allele it does not carry
  (`0|12` matched allele `1`) and (2) a sample carrying two different ALTs (`2|1`) had its second
  allele dropped. The genotype is now parsed into integer alleles and the sample is recorded under
  every ALT it carries. CRISPRme+'s shipped indices are `bcftools norm -m -any` pre-split to
  biallelic (unaffected); this corrects per-sample carrier assignment for **user-supplied
  multiallelic VCFs**. The biallelic path and the enriched FASTA are unchanged.
- **The CRISPR-Bulge scorer is now a guaranteed part of a source install, not best-effort.**
  `install_from_source.sh` builds the `cbulge` env **and compute-self-tests it** (`scorer-env
  create && scorer-env check`), and **aborts the install** if either fails — mirroring the Docker
  image (`build_scorer_envs=1`, failure fatal). Opt out deliberately with
  `CRISPRME_SKIP_SCORER_ENV=1`.

### Fixed
- **A missing/broken CRISPR-Bulge scorer env no longer silently degrades a search to CFD-only.**
  `complete-search` now **fails fast up front** (before the long search) with the one command that
  fixes it when the scorer env is not ready, instead of discovering it per-chromosome in
  post-analysis. A deliberate CFD-only run is still possible with `CRISPRME_SKIP_SCORER_ENV=1`,
  which proceeds with a single loud notice. The scorer-unavailable banner is written to **stdout,
  not stderr** (the post-analysis pipeline treats any stderr byte as fatal, so the old banner
  turned a deliberate CFD-only run into a failed run — supersedes #59). +1 regression test.

## [2.7.2] - 2026-10-08

### Fixed
- **#46 (co-occurrence path): Tier-1 genotype axis was not filtered to VCF-genotyped samples.**
  The phantom-sample filter (`genotyped_samples`) was wired on the Tier-0 registry
  (`compile_from_dict`) but not on the Tier-1 genotype store (`compile_genotypes_from_dict`), so an
  over-listing samplesID inflated the Tier-1 sample axis / AN and **deflated the per-sample
  CONFIRMED-cis joint-AF denominator**. The same `genotyped_samples` set is now threaded through to
  `tier0_compile.build_sample_meta` (already filter-capable); no-op when `None` (batteries install /
  source VCF absent). Build-path only — correctly-built shipped indices are unaffected (e.g.
  HGDP1kGP's header-derived 4,091-sample panel has no phantoms). Added
  `test_genotyped_samples_filter_drops_phantom_from_axis`.

## [2.7.1] - 2026-10-07

### Fixed
- **Report: `IntOGen_cancer_driver` column rendered `-` in every row.** `_curated_cell()` had no
  `intogen` branch, so the value fell through to the trailing `else: v = None` and displayed as `-`
  for every off-target. Added the missing `elif kind == "intogen"` branch (regenerating the curated
  top-1000 goes from 0/1000 to 45/1000 populated IntOGen cells; spot-checked on BRCA2 / EPHA3 /
  PTPRK coding-region off-targets). Added a guard that every `kind` declared in `CURATED_COLUMNS` is
  matched by a branch in `_curated_cell`. (#58, @anncir1)
- **CLI: `generate-report --no-maf` was silently ignored.** The documented top-level
  `crisprme.py generate-report` wrapper hand-parsed only its four value-flags, so the boolean
  `--no-maf` was dropped with no error and the MAF column was still emitted. It is now parsed
  (`"--no-maf" in args`) and threaded into `build_report(drop_maf=...)`, listed in the help, and
  covered by a test asserting the MAF column is actually removed.
- **Web: per-cluster view/download silently dropped all 24 CRISPR-Bulge columns.**
  `change_headers_bestMerge.py` reordered the 3-scorer `.bestMerge.txt` through a fixed `new_order`
  list that predates CRISPR-Bulge and contains no `CRISPR_BULGE_*` names, so the pandas column
  selection silently discarded the entire CRISPR-Bulge block (no crash — all listed names were
  present) from the genomic-region cluster download. The block is now preserved (no-op on a
  CFD-only input); new `test_change_headers_bestmerge.py` locks it in.
- **Web: the Reset button never cleared the threshold dropdown.** The reset callback's Output
  targeted the prop `"value "` (trailing space) instead of `"value"`, so Dash silently never wrote
  it. Fixed, plus `test_dash_output_props.py` guards every `pages/` Output/Input/State prop against
  stray-whitespace / empty names (the same silent-inert class as the IntOGen column).

### Added
- Coverage tests for the "declared-but-dead" bug class (an output field/flag that is declared but
  has no populating/wiring code, silently rendering a placeholder/no-op in every case): the
  CURATED_COLUMNS↔`_curated_cell` guard (#58), the `--no-maf` drop assertion, the CRISPR-Bulge
  per-cluster preservation test, and the Dash callback-prop whitespace guard.

## [2.7.0] - 2026-10-07

### Added
- **Scorer-env build is registry-driven + compute-verified (robust release).** The Docker image
  now provisions **and COMPUTE-self-tests every env in the `scorer_env.SCORER_ENVS` registry**
  (via `scorer_env.build_all`) instead of hardcoding `cbulge` — the modular scorer registry is the
  single source of truth, so adding a scorer env automatically builds it into the image. The build
  runs a real compute self-test (`scorer_env.selftest`: scores a perfect-match vs a mismatched
  target and asserts real, in-`[0,1]`, non-degenerate scores) rather than an import-only check, so
  a broken / miscomputing / weight-less scorer **fails the image build** instead of silently
  shipping (catches the tensorflow-metal GRU-miscompute class). `crisprme.py scorer-env
  check`/`doctor` now run the same compute self-test.
- **New recommended default variant index: `NRG_3_hg38+hg38_HGDP1kGP`** — the gnomAD HGDP + 1000
  Genomes callset (4,091 samples = 929 HGDP + 3,162 1000 Genomes), jointly **SHAPEIT5-re-phased** so
  **both** cohorts are phased: HGDP off-targets now get **CONFIRMED cis** + named per-sample carriers
  (not PUTATIVE-from-phasing), covering HGDP's seven genetic regions. Published to HuggingFace
  (`indexes/NRG_3_hg38+hg38_HGDP1kGP.tar.gz` + `genotypes_hg38_HGDP1kGP.tar.gz`); SNP+indel
  searchable; supersedes the hybrid `NRG_3_hg38+hg38_1000G2021_HGDP`.
- **PAR-aware chrX ploidy.** chrX allele frequencies now use the correct ploidy per region — diploid
  in the pseudoautosomal regions (PAR1/PAR2, AN=8182 for a 4,091-sample panel) and haploid-male in
  the non-PAR region (AN=5990) — in both the registry build and the search-time k≥2 co-occurrence
  denominators. Fixes a ~1.37× PAR allele-frequency inflation; autosomes/chrY unchanged.

### Changed
- **A missing/broken ML scorer env is now LOUD.** A search that degrades to CFD-only because the
  CRISPR-Bulge env is absent/broken prints an un-missable banner to **stdout + stderr** (was a
  one-line stderr note buried in `log_error`) — so a CFD-only run is never silent.
- **Default variant index flipped to `NRG_3_hg38+hg38_HGDP1kGP`** (was `NRG_3_hg38+hg38_1000G2021`).
  The web UI picks it automatically (richest genotyped-phased panel, by sample count — no hardcoded
  name); CLI help + docs recommend it first. `1000G2021` remains the 1000-Genomes-only alternative;
  `1000G2021_HGDP`, `HPRC`, and `mega` remain available.

## [2.6.2] - 2026-09-28

### Added
- **`assembly-search`: direct haplotype-vs-haplotype reconciliation (impg / minimap2).** The two
  personal haplotypes are now reconciled directly against each other (cached, threaded pairwise
  alignments), independent of `liftOver` — more accurate two-haplotype off-target reconciliation.
  `minimap2` / `impg` are checked at run start (`check_impg_available()`). (#53)
- **`assembly-search`: hg38 annotation of reconciled sites (`--annotation`).** Reconciled
  off-target sites can be annotated against the hg38 functional / gene tracks; a failed or missing
  annotation degrades gracefully and never discards an otherwise-finished two-haplotype search. (#55)

### Fixed
- **Non-human genome robustness.** Genome-layout detection and empty-result handling fixes
  surfaced while running CRISPRme on a non-human (pig) genome. (#56)
- **Restored the brute-force benchmark CI gate (`validate-benchmarks`).** The v2.6.0
  `complete-test` rewrite made `run_crisprme_test` always require a prebuilt index and run a
  single example-guide smoke, which silently dropped the CI's download + per-benchmark
  on-demand-build flow — so `validate-test` had no predictions to compare against the committed
  brute-force ground truth and the gate failed. `run_crisprme_test` now branches: the default
  (user) path keeps the no-auto-build prebuilt-index smoke, while the CI path (when
  `CRISPRME_ALLOW_ONDEMAND_BUILD` is set, as the workflow does) downloads the test data and runs
  every registered benchmark with on-demand index build, producing the per-benchmark output dirs
  `validate-test` needs. The brute-force comparison against CRISPRme's predicted off-targets is
  exercised again.
- **Guardrails so this class of regression can't recur silently.** Added `test_complete_test.py`
  — a fast, hermetic dispatch test (in the push-to-main unit suite) that pins the two-mode
  contract: CI mode (`CRISPRME_ALLOW_ONDEMAND_BUILD` set) downloads + builds each registered
  benchmark on demand with no `--index-path`, and user mode requires a prebuilt index and never
  downloads/builds. It fails against the broken rewrite. Also made `validate-benchmarks` (the
  brute-force gate) run on **direct pushes to `main`** touching the pipeline — not just PRs +
  the weekly cron — so a fast-forward release can no longer bypass it (the gap that let the
  breakage sit unnoticed since mid-September).

## [2.6.1] - 2026-09-28

### Added
- **The default image is now GPU-capable — one image, both modes.** `pinellolab/crisprme:latest`
  (and `:<tag>`) builds its **amd64** layer with the conda-forge CUDA TensorFlow
  (`tensorflow=2.13=cuda*`, validated correct + fast on an NVIDIA A100), so the CRISPR-Bulge scorer
  runs on the GPU under `docker run --gpus all …` and **falls back to CPU transparently** when no
  GPU is visible (the compute-backend device guard; byte-identical results) — no flag needed for
  CPU. The **arm64** layer stays CPU (CUDA has no ARM build), so Apple Silicon is unaffected. A
  single `docker pull pinellolab/crisprme:latest` therefore works everywhere and uses the GPU when
  present; the only cost is a larger amd64 image (bundles cudatoolkit/cudnn). Build the Dockerfile
  with `--build-arg scorer_backend=cpu` for a lean CPU-only image if size matters. Selection is via
  the Dockerfile `scorer_backend=auto` default keyed on `TARGETARCH`.
- **`Gene_region` column — exon/intron/UTR granularity for the gene annotation, on every
  report surface.** States in one word where an off-target sits relative to the nearest
  protein-coding gene: `CDS` (coding exon), `5'UTR`, `3'UTR`, `exon` (non-coding-transcript
  exon), `intron` (inside a gene body but not an exon), or `intergenic`. Computed from the
  existing `Annotation_GENCODE` feature-set (no new pipeline pass) and kept coherent with the
  `Gene` / `Gene_distance_kb` columns — `intergenic` is pinned to `Gene_distance_kb > 0`, so
  the region never contradicts the distance. Surfaced **consistently everywhere** from ONE
  shared helper (`utils.gene_region_class`): the raw `integrated_results.tsv`
  (`Annotation_gene_region`, written by the integrator), the in-report curated table, and
  every download TSV (all tiers, panel, top1000, variant_created); documented in the report's
  annotation legend. The raw `GENCODE` column still keeps the full overlapping feature-set.
  Report curation also derives it on the fly, so it renders on pre-existing
  `integrated_results.tsv` files too.

### Fixed
- **GPU image now actually builds its CUDA scorer env (and a broken scorer env can no longer
  ship silently).** The `cbulge` CUDA solve failed on GPU-less builders (e.g. GitHub CI) with
  `__cuda missing`, because the micromamba base activation exports `CONDA_OVERRIDE_CUDA=""` and
  `scorer_env` used `setdefault` (which will not replace an existing empty value). Now the
  override is force-set to `11.8` when empty/absent, so the CUDA TensorFlow build resolves on any
  builder (the image still falls back to CPU at runtime when no GPU is present). Additionally, the
  Dockerfile scorer-env build step is now **FATAL** — a failed or incomplete `cbulge` env fails
  the image build (verified by importing TensorFlow) instead of being swallowed by a warning and
  shipping an image with no scorer.

## [2.6.0] - 2026-09-28

### Changed
- **New default search: 4 mismatches + 1 DNA + 1 RNA bulge (6-edit budget).** The web form
  and CLI now default to `--mm 4 --bDNA 1 --bRNA 1` with `--max-total-edits 6` (was a 6/2/2
  ceiling with a 4-edit budget). This keeps the default in the CRISPR-Bulge scorer's validated
  **single-bulge** domain and, crucially, sizes the total-edit budget to the params so a bulged
  off-target at the full mismatch budget is no longer silently pruned. Simple mode (Advanced
  panel closed) fixes the per-type bulge cap at 1; the CLI likewise defaults omitted
  `--bDNA/--bRNA` to a single bulge of each type. Open the Advanced panel / pass explicit
  `--bDNA/--bRNA` (and raise `--max-total-edits` to their sum) for deeper searches.
- **Per-sample availability + label are index-aware on the web.** The "Per-sample" analysis
  option is auto-enabled only for genotyped panels (grayed out + reverted for sites-only like
  `mega`), and its label now states what it buys on the chosen index: **CONFIRMED cis + named
  carriers on a phased panel**, **carriers with PUTATIVE cis on an unphased panel**, or the
  hybrid mix.
- **CRISPR-Bulge is the machine-learning off-target scorer.**
  The legacy vendored RandomForest (scikit-learn) ML scorer is retired in favor of CRISPR-Bulge
  (Yaish & Orenstein, *NAR* 2024; TensorFlow GRU ensemble; MIT), which is more accurate
  — ~2.1× overall and ~5.4× on bulge-containing off-targets by AUPR on Refined-TrueOT —
  and faster per off-target on CPU. **CFD remains the primary score** (fast, in-process,
  pure-Python); CRISPR-Bulge is the secondary ML score. To keep the main environment free
  of a heavy TensorFlow pin, the scorer runs in its **own dedicated conda env** (`cbulge`),
  invoked over a persistent-worker batch protocol; if that env is absent the run completes
  on CFD with a single warning (the ML column is `-1`). Reports and the web UI now label
  the column/plots/legend "CRISPR-Bulge" and calibrate on its own [0,1] scale
  (0.5 / 0.2 / 0.1 bands); the physical `CRISPR_BULGE_score_*` column name and `.bestCRISPR_BULGE`
  file are retained as stable internal identifiers so older reports still render.
  The ensemble is a single-bulge model, so off-targets whose alignment needs ≥2 bulges (or a
  ≥2-bp bulge) are NOT ML-scored — they are nulled to `-1` and carried by CFD (the primary
  score) + mismatch/bulge counts. (An earlier build reduced such alignments to a single bulge;
  an adversarial review showed that reduction could corrupt the alignment and, on a non-monotonic
  model, was not provably conservative, so it was removed in favor of the honest null. Each
  off-target locus is still ML-scored on its best ≤1-bulge alignment when one exists.)
- **Modernized the main environment now that the legacy ML scorer's hard pins are gone.**
  `scikit-learn` is removed entirely (nothing in the main env imported it after the legacy
  RandomForest scorer was retired);
  `numpy` 1.24→≥1.26, `scipy` 1.10→≥1.11, `pandas` 2.0→≥2.1, and the `matplotlib-base` ceiling
  is dropped. `numpy` is capped `<2` and `pandas` `<3` as deliberate guards (the pysam /
  CRISPRitz C-ABI, and pandas 3.0 copy-on-write) — both are validated fast-follows. The
  full unit suite passes on the rebuilt, un-pinned environment (783 tests on the
  fresh-from-`environment.yml` source install).
- **The CRISPR-Bulge scorer env is architecture-aware, with identical scores on x86-64 and
  arm64.** conda-forge ships `tensorflow-cpu=2.13` only for x86-64 (aarch64 has just 2.18/2.19),
  so the `cbulge` env is now built per-arch behind the same name: conda `tensorflow-cpu=2.13` on
  x86-64, the conda CUDA build on GPU hosts, and — new — the **same TensorFlow 2.13.1** from the
  official PyPI aarch64 wheel on Linux ARM (base deps from conda, TF from pip, mirroring the
  Apple-Silicon Metal path). Same version + weights → validated **identical**: on a 254-pair panel
  (150 with bulges) x86-64 vs aarch64 scores agree to max |Δ| = 3.6e-7 with **zero** differences at
  display precision and **zero** threshold-band flips. The multi-arch Docker image therefore ships
  the full scorer on both arches (previously arm64 would have degraded to CFD-only). The scorer-env
  solve also honors the same `CONDA_CHANNEL_BASE` mirror as the main env (e.g. `https://prefix.dev`),
  so `scorer-env create` works on networks that block `conda.anaconda.org`.
- **Alternative-alignments output is now mode-driven (off by default in the
  population-level analysis).** The `..._all_results_with_alternative_alignments.tsv`
  dump — the *non-best* alignments per locus — grows combinatorially with the edit
  budget and its annotation/sort/integration was the single most expensive
  post-analysis step on a large search, yet the report, summary tables, mismatch/bulge
  matrix and top-100 panel are all built from the *best* alignment per locus and never
  read it. It is now emitted only under `--per-sample`; the default population-level
  run drops the non-best rows right after the merge, skipping the expensive step. New
  `--alt-alignments` / `--no-alt-alignments` flags override either way.
  `assembly-search` no longer forces it on either (its haplotype reconciliation reads
  only `*_integrated_results.tsv`), so its per-haplotype searches run lean too. (A worst-case
  6mm/2+2/max-8 search on the 5-source mega had ballooned this file past 160 GB and
  made the run take days; the default population path no longer pays that cost while
  the report is byte-identical.)
- **License: CRISPRme+ is now the MGB Open Access License 1.0 (MGBOA 1.0).** The
  software moves from AGPL-3.0 to MGBOA 1.0 (non-commercial, non-revenue-generating
  academic use; commercial use requires a license). `LICENSE` carries the full MGBOA
  text; the README, GitHub Pages footer, and the generated-report footer are updated
  to match (the report footer uses the MGB-approved "License and permitted use" +
  disclaimer language). The main `pinellolab/CRISPRme` repo is relicensed to match.
- **Report: the off-target Mismatch/Bulge matrix is now resolved by bulge TYPE
  (DNA vs RNA).** The single matrix splits each row into DNA-bulge and RNA-bulge
  counts (derived from the guide-aligned gap positions), so off-targets are countable
  by bulge type; the coarse total-bulge-only matrix is removed (the DNA/RNA table
  includes the 0-bulge row, so it is a strict superset and keeps the perfect-match
  highlight).
- **Report: variant identifiers are shown as `chrom;pos;ref;alt`** (with the rsID as a
  companion), instead of preferring the rsID — normalizing the previously mixed
  rsID / `pos:ref:alt` / `chr_pos_ref_alt` forms; multi-SNP haplotypes are joined.
- **Report: the high-variant-density column now reads "N in window"** (was
  "Yes (N var)"), and the legend spells out that this counts total variants in the
  surrounding window — distinct from the `Variant` column, which lists the variant(s)
  that create the individual off-target.
- **Personal-assembly storage is now one self-contained folder per individual with a
  single `metadata.json`.** `Assemblies/<individual>/{paternal,maternal}/genome/*.fa`
  plus each haplotype's `*_vs_GRCh38.chain.gz` and `*.chromAlias.txt`, described by a
  relocatable `metadata.json` (paths relative to the folder). This replaces the flat
  layout that scattered files across `Genomes/` + `LiftoverFiles/` with six hidden
  `.assembly_individual` marker files. Benefits: deleting an individual removes one
  folder and reclaims all its space (the old delete only removed markers), clear
  separation between individuals, and portable bring-your-own bundles. The old flat
  layout is still discovered for backward compatibility. New shared module
  `PostProcess/personal_assembly.py` (+ `test_personal_assembly.py`).

### Added
- **`assembly-search` annotates reconciled sites against hg38.** New optional
  `--annotation` flag (same file format, COSMIC licence gate and IntOGen handling as
  `complete-search`'s) adds an `Annotation` column to `combined_hg38.tsv` for every site
  with an hg38 coordinate (found on both haplotypes, or on one and mappable); sites with
  no hg38 equivalent get none. The web form applies the annotations enabled in
  Settings -> Annotations, like `complete-search`. In the results table the column is
  filterable and sortable like every other one. The annotation describes the *lifted*
  hg38 locus, not the haplotype's own sequence.
- **`scorer-env` command group + modular scorer environment.** `crisprme.py scorer-env
  {create,check,update,list,doctor}` manages the dedicated `cbulge` conda env (detects
  micromamba/mamba/conda, provisions the pinned CRISPR-Bulge source+weights, health-checks
  imports/weights, records state under `Annotations/.scorer_env.json`). The Docker image
  and `install_from_source.sh` create it at install time; set `build_scorer_envs=0` /
  `CRISPRME_SKIP_SCORER_ENV=1` to defer. `CRISPRME_SCORER_ENV` points CRISPRme at an
  existing/shared scorer env instead of creating its own.
- **`--compute-backend {cpu,gpu,auto,cuda,metal}` for the ML scorer** (default `cpu`;
  GPU is opt-in and always degrades gracefully to CPU). `cuda` uses a CUDA TensorFlow build
  on NVIDIA/Linux; `metal` uses tensorflow-metal on Apple Silicon **guarded by a load-time
  numerical self-test** — tensorflow-metal silently miscomputes this model's GRU kernel, so
  the scorer detects the miscompute and transparently falls back to a CPU device context
  (bit-identical scores) rather than emit wrong ones. A CPU thread cap
  (`CRISPRME_SCORER_THREADS`, default `min(16, cores)`) prevents TensorFlow oversubscription
  on many-core nodes. **CUDA is validated on an NVIDIA A100** (Linux): the GRU is numerically
  correct (~1e-5 of CPU) and scoring runs ~14.7k off-targets/s — ~1.6× the CPU rate, so a
  genome-wide worst-case (6.77M off-targets) scores in ~8 min. The Docker image accepts a
  `--build-arg scorer_backend=gpu` to bake in the CUDA scorer env (run with
  `docker run --gpus all … --compute-backend cuda`); CFD stays on the CPU alongside it.
- **New recommended default variant index: single-source, phased
  `NRG_3_hg38+hg38_1000G2021`** (published to HuggingFace). 1000 Genomes 2021,
  3,202 samples, fully phased (`data_type=genotyped-phased`) — so every reported
  co-occurrence is CONFIRMED cis with named per-sample carriers throughout (no
  PUTATIVE-from-phasing): a clean, simple default of real observed haplotypes.
  `NRG_3_hg38+hg38_1000G2021_HGDP` remains available as the broader-coverage option
  (adds HGDP's 929 individuals; hybrid confidence). The docs and CLI download hints
  now point at the single-source default first.
- **`--alt-alignments` / `--no-alt-alignments`** (`complete-search`) — explicitly
  emit or suppress the exhaustive alternative-alignments file, overriding the
  mode-driven default (off in population-level, on under `--per-sample`).
- **Index `data_type` / `phased` metadata in the registry manifest.** Each
  `reg_<chrom>.idx` now records `data_type` (`sites-only` / `genotyped-unphased` /
  `genotyped-phased` / `hybrid`) and a `phased` flag (per-database + overall),
  observed from the genotypes at build time — so tooling can identify an index's
  type without scanning the multi-GB per-sample dictionary. The variant-dataset
  selector in the web UI reads it to label each dataset (falling back to
  genotype-store presence for indexes whose manifest predates the flag). The three
  shipped indexes (`1000G2021_HGDP` = hybrid, `HPRC` = genotyped-phased, `mega` =
  sites-only) were backfilled in place on HuggingFace.
- **Report: an `Observed` column flagging off-targets supported by ≥1 real
  individual.** Values are `reference` (present in the reference genome, so carried
  by essentially every individual), `N carrier(s)` (a variant site carried by N
  named individuals in a genotyped panel), `observed` (a single AF>0 variant —
  carried by ≥1 individual by definition, even where carriers can't be named), or
  `putative` (a multi-variant combination whose cis co-occurrence cannot be
  confirmed). On a genotyped index the carrier roster is populated in both analysis
  modes; under `--per-sample` it is pruned to the EXACT cis carriers, so the column
  cleanly separates genuine haplotypes from worst-case reconstructions. The
  recommended validation panel now **prioritizes observed (real-carrier) sites**
  among equally-severe candidates — as the first tie-break *after* worst-case
  severity, so a strictly-worse reconstructed site is never displaced (worst-case
  coverage preserved), and the panel summary reports how many selected sites are
  observed. On a sites-only panel (e.g. `mega`) the `N carrier(s)` value and the
  tie-break are inactive (no per-sample roster) but the column is still shown; it is
  dropped only for a reference-only run.
- **`assembly-search --assembly-individual <name>`** — resolve both haplotypes'
  genome/chain/chromAlias from a registered `Assemblies/<name>/` bundle instead of
  passing the six paths explicitly (the six flags still work, incl. the legacy layout).
- **`assembly-search --max-total-edits <N>`** — override the total mismatch+bulge
  budget for the per-haplotype searches. The default is unchanged (the full
  `mm + bDNA + bRNA` budget — a personal assembly is a plain resolved genome with no
  IUPAC-variant lattice, so it does not need `complete-search`'s variant-density cap of
  4), but the override lets you tighten or widen it explicitly. Also fixes reuse
  detection so a re-run with the same parameters reuses the prior per-haplotype search.
- **`assembly-search` run history + self-describing jobs.** Each assembly-search job
  now records its parameters (individual, guide, PAM, edit budget, mode) alongside the
  results, and the web **History** page lists past assembly-search runs next to the
  standard searches so they can be reopened.
- **Bring-your-own personal assembly** — the "Add a personal assembly" Settings card has a
  file-picker to upload a `.zip`/`.tar.gz` bundle (one top-level `<individual>/` folder + a
  valid `metadata.json`); it is extracted, validated, and registered
  (`personal_assembly.import_archive`). Extraction is path-traversal-safe (rejects
  absolute / `..` members and tar symlink/device entries, and verifies every member stays
  within the staging directory).
- **Single canonical installer** — the clickable-app `install.sh`/`install.ps1` now live
  only in the served GitHub Pages repo (`pinellolab/CRISPRme` `docs/`); the drift-prone
  duplicate under `install/` was removed and replaced with a pointer README. Docs clarify
  that only the whole-genome variant-aware search needs ~64 GB RAM (many current laptops
  have it); reference-only / single-chromosome / region searches fit ~16 GB.
- **Direct haplotype-vs-haplotype reconciliation for `assembly-search`.** A
  subset of each haplotype's non-mappable (no hg38 equivalent) off-target
  predictions are now resolved by aligning the two haplotype genomes
  directly to each other (`minimap2`/`impg`), independent of `liftOver`.
  A prediction is only merged into the new `both_haplotype_private`
  category when both a single unambiguous alignment hit and mutual
  (bidirectional) reciprocal confirmation hold — real, measured directional
  asymmetry between haplotypes required this over a simpler one-directional
  check. `minimap2`/`impg` are new required dependencies for
  `assembly-search` (`check_impg_available()` mirrors the existing
  `check_liftover_available()` gate). 92/92 tests, real-data verified. The
  live results page and the static combined report both surface the new
  category: the "Both haplotypes" haplotype-coverage bar now has a real
  unmapped segment (previously hardcoded to 0, since the category didn't
  exist yet), and the live page's summary stats gained a "Found in both
  haplotypes, non-mappable to hg38" card. Reports from before this
  category existed render unaffected (the new count defaults to 0).

### Removed
- **The legacy vendored RandomForest ML scorer and azimuth are removed.** The 276 MB vendored
  RandomForest model (predictors archive, `dnaShape.pkl`, and its scorer module), the whole
  `PostProcess/azimuth/` tree (already inert / never invoked), and their now-unneeded
  scikit-learn 1.1.3 / numpy 1.24 pin are deleted. CRISPR-Bulge is now the ML scorer (see
  *Changed*); CFD's small lookup pickles (`mismatch_score.pkl`, `PAM_scores.pkl`) are
  version-agnostic and unaffected. `--scorer` is gone (the scorer is always CRISPR-Bulge); the
  internal `CRISPR_BULGE_*` column names are kept for back-compat only.
- **Legacy per-file personal-assembly upload in Settings.** The old "Add a personal
  assembly" card that took the six individual files (two genomes, two chains, two
  chromAlias) one at a time is gone; personal assemblies are now added as a single
  self-contained `.zip`/`.tar.gz` bundle (see *Bring-your-own personal assembly*
  above), matching the one-folder-per-individual storage layout. The CLI still accepts
  the six explicit paths for scripted use.

### Fixed
- **A genome delivered as one multi-sequence FASTA file no longer silently produces an
  EMPTY result.** CRISPRme identifies chromosomes by FASTA *file name*, one sequence per
  file (as the README requires), but CRISPRitz happily indexes and searches a single
  multi-sequence file -- so a UCSC assembly with no per-chromosome archive (e.g. the pig
  `susScr11`, `bigZips/susScr11.fa.gz`) ran to completion and then reported "No
  off-targets found" even with tens of thousands of raw hits, because the whole genome
  was treated as one "chromosome" named after the file. Now: `download` /
  the Data Manager split such a download into one `<sequence>.fa` per sequence; the
  input validator and the pipeline itself refuse a multi-sequence genome folder up front
  with the exact fix (`python PostProcess/genome_layout.py --split <genome_dir>`, which
  also repairs an existing folder -- rebuild any index made from the old layout); and the
  pipeline no longer reports success when the search found raw hits but none survived
  post-analysis (it fails with an explanatory error instead). New
  `PostProcess/genome_layout.py` (+ `test_genome_layout.py`); the layout check reads one
  header line per file, so it adds no measurable time to a search.
- **A search that genuinely finds nothing now says so on the job status page** (green
  "Finished: no off-targets were found ...", later steps marked "Not needed") instead of
  leaving "Merging Targets" and the rest as "To do" indefinitely.
- **The Result Summary matrix no longer shows bulge rows `0..5` for a reference-only
  search that only reaches 2 bulges.** The per-guide count file always carries a
  reference block and an all-zero variant block; the page only split the two blocks for
  searches with variants. A reference-only job now shows just its reference block.
- Removed two Dash dev-mode prop-type warnings on the results page (invalid extra keys in
  a `DataTable` `css` entry; a hidden radio's options passed as `[labels, values]`
  instead of one list of `{label, value}`). No visible change.
- **`complete-test` no longer relies on on-demand index building.** The smoke test ran a
  variant `complete-search` that expected an index to be built on demand — a path that was
  gated off (it produced incomplete dict-less tiers), so `complete-test` failed. CRISPRme does
  **not** build an index automatically: `complete-test` now checks for a prebuilt index and, if
  none is present, prints clear guidance to **download** one (`crisprme.py download --what index
  --index-name …`) and exits cleanly without downloading test data or building anything — the
  user decides which index(es) to fetch (as the web instructions already recommend). When an
  index is present it is used via `--index-path`.
- **Web interface now boots on the modernized environment.** The Dash cache was configured
  with the Flask-Caching **1.x** short-name `CACHE_TYPE: "filesystem"`; Flask-Caching **2.x**
  (pulled in by the un-pinned, modernized web stack) removed the short-name aliases and expects
  the backend class name, so `filesystem` raised `ImportStringError` at import time and the whole
  web app failed to start (`crisprme.py web-interface` and the served site). `CACHE_TYPE` is now
  `"FileSystemCache"`, which resolves on Flask-Caching 2.x (and remains valid on 1.10+).
- **Web results page no longer crashes (or silently drops a plot) when the CRISPR-Bulge
  filter is selected.** The filter dropdown value is the hyphen form `CRISPR-Bulge`, but the
  cached top-1000 / personal / private lolliplot files are written with the underscore token
  `CRISPR_BULGE`; the results page interpolated the hyphen value straight into the filename, so
  the personal-card page raised `ValueError('Personal and Private Lolliplots not found')` and
  the per-guide top-1000 plot vanished. The filename token is now normalized
  (`CRISPR-Bulge`→`CRISPR_BULGE`); `CFD`/`fewest` are unaffected.
- **The report now states, as a live count, how many off-targets are ≥2-bulge (N/A for
  CRISPR-Bulge).** The recommended-panel note reads "N of M off-targets are ≥2-bulge and
  therefore show `-` (N/A) in the CRISPR-Bulge column — ranked by CFD and edit distance
  instead, not dropped," so a mostly-empty ML column on a bulge-heavy run reads as expected
  behavior rather than missing data.
- **Variant search no longer aborts during report generation** (population-distribution
  step). `process_summaries.py` writes the per-criterion sidecars using the filter-criterion
  token (`..._CRISPR-Bulge.txt`, hyphen), but `submit_job` read the population-distribution
  input — and error-path cleanup globs — with the underscore identifier form
  (`..._CRISPR_BULGE.txt`); the mismatch raised `FileNotFoundError`, wrote to stderr, and the
  `[ -s $logerror ]` guard aborted the whole run with "population distribution plots creation
  failed" on **every** search. The consumers now use the hyphen form that `process_summaries`
  actually emits.
- **`build-index-only --vcf` now accepts the same config-file form as `complete-search`.**
  `complete-search --vcf` takes a file listing VCF folder names (under `VCFs/`), but
  `build-index-only --vcf` accepted only a bare directory path — so the same flag meant two
  different things. `build-index-only` now accepts **either** a VCF directory (unchanged) **or**
  a one-line config file naming a folder under `VCFs/`, with a clear error if the listing names
  more than one dataset (build each separately / merge first). Help text + README updated.
- **`complete-search`'s "build the index first" hint now shows the correct flags.** The error
  told users to run `build-index-only ... --bMax <N>`, but that command reads `--bDNA`/`--bRNA`
  (and ignores `--bMax`), so following the hint verbatim silently built nothing. The message now
  says `--bDNA <N> --bRNA <N>`.
- **Report ≥2-bulge N/A wording is precise.** The scores legend said such off-targets are
  "shown as -1"; the report tables actually show `–` (only the raw `integrated_results.tsv`
  records `-1`). The legend now states both.
- **`assembly-search` no longer crashes when a haplotype has zero off-targets.** A
  haplotype whose search finds no hits on the assembly is a legitimate result, not an
  error; reconciliation now treats it as an empty set — every locus on the other haplotype
  becomes heterozygous-equivalent (`<other>_only`), and the empty side contributes 0
  non-mappable — instead of aborting on the missing `*_integrated_results.tsv`. Unsupported
  multi-guide runs (more than one results file) still error.
- **`assembly-search` haplotype reconciliation reads only `*_integrated_results.tsv`.**
  It previously also outer-merged in `*_all_results_with_alternative_alignments.tsv`
  and re-clustered by CFD, on the theory that `integrated_results.tsv` alone
  was missing real distinct sites. Re-derived: that theory was wrong — the
  union's "duplicates" (up to 233 rows at one locus) were alt-file-only
  representations of an already-counted site, not additional real
  off-targets, and `integrated_results.tsv` already has exactly one row per
  real physical cluster by construction. The representative-selection
  criterion is also now consistently fewest-mismatches+bulges (matching
  `merge_contiguous_targets.py`'s own choice) rather than highest-CFD at
  every call site. Real genome-wide effect on a real dataset: net one fewer
  `both` row (a duplicate collapsing), zero effect elsewhere. 58/58 tests.

## [2.5.5] - 2026-09-13

### Added
- **HPRC pangenome production index** (`NRG_3_hg38+hg38_HPRC`) — a third prebuilt SpCas9
  NRG variant index, published on HuggingFace alongside the genotyped 1000G-2021 + HGDP and
  the sites-only mega. It enriches hg38 with the HPRC Minigraph-Cactus pangenome decomposed
  to a per-chromosome VCF (`vg deconstruct` → `vcfbub` → `vcfwave`; 232 assembly-derived
  genomes incl. CHM13). Genotypes are **phased**, so it supports CONFIRMED cis co-occurrence
  and named per-sample carriers, and it captures pangenome-specific variation absent from
  short-read panels. Reproducible recipe: `seq_script/merge_panels/hprc_build.sh` (RAW
  registry + SNP+indel, consistent with the other two indexes). METHODS §merge documents it.
- **`assembly-search` web interface** (PR #48) — the diploid personal-genome search
  (paternal + maternal assembly FASTAs, liftOver-reconciled to hg38; no VCF) is now available
  on the website: a "Personal assembly" tab under Select genome, a parallel status-polling
  `/load` page, and a dedicated results view (`result_page_assembly`) that mirrors
  complete-search's layout for two-haplotype data (haplotype-coverage stats, reconciled
  off-target table, per-haplotype non-mappable sites, Graphical Reports). A combined
  report.zip is produced at parity with complete-search (`build_combined_report`), and a
  Settings "Add a personal assembly" card fetches an HPRC Release-2 individual
  (`download_hprc_assembly.py`) or accepts uploaded assemblies.

### Changed
- Pinned **`dash-bootstrap-components<2`** in both `environment.yml` (PR #48) **and the
  `Dockerfile`** — the Docker build installs deps directly (not from `environment.yml`), so the
  image pin is what actually reaches users. dbc ≥ 2.0.0 requires React 18 / Dash ≥ 3; this app
  targets React 16 / Dash 2.x, and the assembly-search `dbc.Tabs` do not render under dbc 2.x.
  The published v2.5.5 image ships dbc 1.7.x so the web assembly-search feature renders (the CLI
  is unaffected).

### Fixed
- `assembly-search` now passes `--max-total-edits {mm+bDNA+bRNA}` into each haplotype's
  `complete-search` (PR #47) — previously it silently used the default cap of 4, pruning
  mismatch+bulge off-targets.
- `download_hprc_assembly.py`: the `--register` path used a bare `from pages_utils import …`
  (→ `ModuleNotFoundError`, the only mode the Settings card invokes) — fixed to
  `from pages.pages_utils import …`; and the multi-GB scratch `.fa.gz` is now removed after a
  successful split and on MD5-mismatch abort (was leaking ~2 GB per haplotype).

### Verified
- HPRC index validated end-to-end (`--per-sample`): 68,057 CONFIRMED phased-cis off-targets,
  855 SNP+indel and 140 SNP+SNP co-occurring rows, with named cis carriers ground-truth-verified
  against the source genotypes (e.g. HG03579 carries chr10:12896898 `TG>T` + chr10:12896900
  `T>A` in cis, joint AF 1/464). Assembly-search validated e2e (HG01255 chr19, merged code):
  reconciliation + combined report.zip. Full web click-through: tabs render (dbc<2), Load
  Example, Per-sample gray-out survives the assembly-tab round-trip, Settings cards, assembly
  results page. HPRC download card exercised against real HPRC S3 (92 contigs, register, leak-free).

## [2.5.4] - 2026-09-11

### Added
- **IntOGen (CC0) cancer-driver genes as the default cancer-gene annotation.** With COSMIC
  now licence-gated (below), CRISPRme+ ships **IntOGen** — a compendium of cancer *driver*
  genes whose current release is **CC0** (public domain; free for academic AND commercial
  use + redistribution) — **enabled by default** as the licence-free cancer-gene flag. New
  report column `Annotation_INTOGEN` (mirrors `_COSMIC`), a default-on `intogen_drivers.hg38.bed.gz`
  track (`seq_script/build_intogen_annotation.py` maps the CC0 driver symbols to GENCODE gene
  spans), report legend, Settings label. So a default search still flags cancer-driver
  off-targets — legally — and users with a COSMIC licence add the curated Census on top.
- **COSMIC cancer annotations are now licence-gated and EXCLUDED by default.** COSMIC
  (Genome Research Ltd / Wellcome Sanger) is free for academic / non-commercial research
  but **commercial use requires a licence** ([terms](https://www.cosmickb.org/terms/)).
  Because the built-in hg38 bundle bakes COSMIC in, CRISPRme+ now **strips COSMIC features
  from the annotation (and therefore the output) by default**, on both the web and the CLI.
  The user opts in **once** (the choice persists in `Annotations/.cosmic_license.json`):
  a **Settings → "COSMIC cancer annotations (licence)"** checkbox on the web (with a link to
  the terms), or `crisprme.py cosmic-license {enable|disable|status}` on the CLI (`enable`
  asks for confirmation; `--accept` for non-interactive/Docker use). New module
  `PostProcess/cosmic_license.py` (attestation flag + row-level `_COSMIC` filter);
  `build_active_annotation` and `complete-search` consult it. (`PostProcess/test_cosmic_license.py`.)

### Changed
- **Renamed the analysis-mode control: default is now "population-level"; `--full` → `--per-sample`.**
  The measurements showed the two modes differ in **genotype resolution**, not intrinsic speed
  (on a single guide the two take the same wall time; on a sites-only panel the opt-in mode does
  nothing). So "fast/full" was the wrong framing. The CLI now runs a **population-level analysis**
  by default (no flag) and exposes an opt-in **`--per-sample`** that resolves per-sample genotypes
  into observed haplotypes (CONFIRMED cis + named carriers + exact joint AF). The web form's mode
  control is relabeled **"Analysis mode: Population-level (default) / Per-sample"**, and the
  report's Search-mode note + the `.search_mode` marker use the new names. `--per-sample` requires
  a **genotyped** panel: on a sites-only index it is a **no-op with a warning**, and the web form
  **disables** the per-sample option when a sites-only index is selected. The internal
  `CRISPRME_FAST_MODE` env var (post-analysis plumbing) is unchanged and not user-facing.

### Removed
- **`--fast` / `--full` flags removed (clean break, no aliases).** They shipped only in 2.5.3 and
  are replaced by the population-level default + `--per-sample` above. Legacy `.search_mode`
  markers (`fast`/`full`) written by 2.5.3 are still recognized when rendering old result dirs.

### Verified
- **New-user clean-room on the released v2.5.3 (image + both published RAW indexes).**
  `apptainer pull docker://pinellolab/crisprme:v2.5.3` reports `v2.5.3`; a fresh
  `download` of **both** production indexes confirms the Tier-0 registry ships
  **uncompressed (`codec=0`, RAW)** as intended (`registry_hg38_1000G2021_HGDP` and
  `registry_hg38_mega`); a fast-default `complete-search` on **each** completes `rc=0`
  and produces a `report.zip` carrying the **Search-mode note + fast caveat**. Confirms
  the two-index design end-to-end: genotyped (CONFIRMED-capable) vs mega (sites-only,
  PUTATIVE + min-AF), and that the **mega now searches indels** (32-bin fake-indel
  genome; 387 SNP+indel co-occurrence + 393 `indel_af` rows with 5-source per-dataset AF).
- **Measured fast-vs-`--full` comparison on both indexes (adversarially verified).**
  Genotyped 1000G-2021+HGDP: detection is **lossless at the window level** (fast a superset;
  the 24 "full-only" loci are 1–10 bp anchor shifts of windows fast detected); **`fast_CFD ≥
  full_CFD` with 0 violations across 3,079 shared loci**; `--full` **tightens** carrier sets
  (fast over-lists a PUTATIVE union → 1,642 distinct samples; full prunes to observed cis →
  1,447) and turns 0 CONFIRMED co-occurrence rows into **106 CONFIRMED phased-cis** rows with
  exact joint AF + named carriers. Sites-only mega: `--full` yields **0 CONFIRMED / no
  carriers** in either mode and the `indel_snp_cooc` / `indel_af` companions are
  **byte-identical** — so `--full` is meaningful only on a genotyped panel. Single-guide
  runtimes are near-equal (fast ~13.8 min vs full ~13.4 min); fast's speedup is a
  dense/aggregate-panel effect, not a per-guide one.

### Fixed
- **New-user docs pinned stale artifacts.** `README.md` copy-paste commands referenced the
  previous image `v2.5.2` → bumped to `v2.5.3`; `docs/DOCKER_QUICKSTART.md` was a v2.4.0-era
  doc with a **dead index name** (`NRG_3_hg38-dictless+hg38_1000G_HGDP`, whose download would
  fail) → fixed to `NRG_3_hg38+hg38_1000G2021_HGDP`, image tags bumped to `v2.5.3`, the mega
  index added, the "Alpha release" banner replaced, and a fast-vs-full explanation added.

### Documentation
- **METHODS §5/§8 + README** now document the *measured* fast-vs-`--full` behavior above
  (window-level losslessness, the exact worst-case CFD bound with 0 violations, carrier
  tightening, and the verified sites-only inertness of `--full`), so the fast/full contract
  is documented from data rather than design intent alone.

## [2.5.3] - 2026-09-11

### Performance
- **All indexes ship the Tier-0 registry uncompressed (raw) by default — "speed over space."**
  Registry lookups (millions per dense search, scattered across a chromosome) are now direct
  memory-mapped reads with **no per-lookup decompression** — measured **~2× faster** than the
  compressed path even with a warm cache. The genotyped build previously hard-coded zlib
  (`compress=True`); it now defaults to raw, matching the sites-only mega index, so both
  production indexes use **one consistent format**. Optional zlib block-compression (~3.5×
  smaller on disk) remains available via `CRISPRME_REGISTRY_COMPRESS=1` at build time (or
  `transcode_registry`); the reader reads either format byte-for-byte identically. The block
  cache + O(1) LRU below now serve that opt-in compressed path.
- **Tier-0 registry block-cache: ~3.7× faster dense variant post-analysis, genome-wide, byte-identical.**
  The genotyped v3 registry is zlib block-compressed (`~3.6×` smaller on disk); profiling a
  dense per-contig SNP post-analysis showed the cost had moved to lookup-time decompression —
  **`zlib.decompress` was 71% of the run** — because the decompressed-block LRU was only 8
  blocks while a dense IUPAC decomposition touches positions across the whole chromosome, so
  hot blocks were re-decompressed ~10⁶ times. The cache default is raised **8 → 4096 blocks**
  (tunable via `CRISPRME_REGISTRY_CACHE_BLOCKS`); the LRU only grows to the blocks actually
  touched, so a block being ~64 KB, this "caches the whole contig registry" for any single
  contig. **Measured** on a dense chr22 SNP post-analysis: **326 s → 89 s (3.66×), output
  byte-identical** (same md5). Sized from measured block counts — genotyped chr22 = 339
  blocks (~22 MB), chr1 (largest) = 1,868 blocks (~122 MB) — so **4096 holds every genotyped
  contig fully**, extending the speedup genome-wide (not just small chromosomes), and the cap
  bounds a pathological panel at ~256 MB/reader. The **mega sites-only index is uncompressed
  (`codec=RAW`) → unaffected** (harmless no-op). Also corrects the earlier "the ML scorer is the
  tail" assumption — the ML scorer is ~4% of the run; the registry decompression was the real tail.
- **O(1) registry block-cache LRU (`OrderedDict`).** Profiling the *indel* post-analysis with
  the enlarged cache exposed a second tail the small cache had hidden: the LRU was a plain
  list, so every cache hit did an O(cache-size) `list.remove()` — **43 s / 26 %** of a dense
  indel post-analysis at ~15 M hits. Replaced with an `OrderedDict` (`move_to_end` /
  `popitem`, all O(1)); byte-identical LRU semantics. **Indel post-analysis 167 s → 121 s**,
  and it also roughly halved the SNP path (cache=512: 89 s list-LRU → 43 s O(1)). A clean
  cache sweep shows the hot working set is tiny (plateau at ~128 blocks), so the combined
  cache + O(1)-LRU fix is **326 s → 43 s ≈ 7.6×** on a dense chr22 SNP post-analysis, and the
  4096 default is generous headroom that fills only to touched blocks (~8-40 MB in practice).

### Changed
- **Fast mode is now the DEFAULT search behavior; `--full` opts into exact enumeration.**
  The SNP variant post-analysis now reports one worst-possible representative per variant
  window by default (exact worst-case CFD; the ML scorer is a best-effort screen), instead of
  enumerating every observed haplotype — so a search stays tractable on dense / aggregate
  panels out of the box (the full enumeration was measured 49 h+ without completing on a
  4×-density panel). **Trade-off: per-sample carriers, CONFIRMED cis phasing and exact joint
  allele frequency are NOT computed in fast mode** — pass `--full` (CLI) or pick **Full** in
  the web search form for that per-sample resolution (recommended for genotyped panels /
  clinical validation). `--fast` is retained as a deprecated no-op alias. The web search form
  exposes a **Search mode** control with an inline explanation, and the shareable report now
  carries a **Search mode** row stating which mode was used and its implications. SNP+indel
  co-occurrence is unchanged by the search mode. A launch-time message and the report note
  make the carrier trade-off explicit (never silent). Mechanics unchanged: the resolved mode
  propagates through the post-analysis subprocess tree via `CRISPRME_FAST_MODE`
  (`0`=full, `1`=fast); a direct `new_simple_analysis.py` call still defaults to full, so unit
  tests and the byte-identical guarantee for the enumeration path are unaffected.

### Added
- **Parallel ML scoring (opt-in, default OFF)** — prototype
  that removes the per-contig ML-scoring tail (one large chromosome's worker running
  serially long after the others finish). The per-contig post-analysis pool is unchanged; this
  parallelizes the CPU-bound per-target **feature build** (`get_features`) — the genuinely
  serial part — at the finest seam across a small bounded **fork** pool.
  The model stays in the parent (children only build features). **Bit-identical feature
  matrix by construction** (a feature row is a pure per-target function; contiguous,
  order-preserving chunks) and default OFF ⇒ output + process count unchanged. Fork (not spawn)
  because the post-analysis entry points run their main at module top level. Measured
  feature-build speedup ≈ 3.8× at 8 workers. The model's
  `predict` is left untouched — the pickled predictors already carry `n_jobs=-1` (all cores), so
  it is already parallel (and, as a result, non-deterministic at the raw-float level: two serial
  runs differ ~5e-17, identical to the emitted 3 decimals — the ML scorer's output has always been
  reproducible only to the emitted precision, which this preserves).

### Fixed
- **HuggingFace download resilience against rate limits (HTTP 429).** A first-time,
  unauthenticated `crisprme download` of the many-small-file indices could hit HF's
  "Too Many Requests" and abort the whole download. The downloader now retries transient
  429 / 5xx responses with exponential backoff (honoring a sane `Retry-After` header) around
  both the snapshot fetch and the size-verification pass; `snapshot_download` resumes, so a
  retry only re-pulls what's missing. 4xx / auth errors still fail fast. Tunable via
  `CRISPRME_HF_MAX_RETRIES` (set to `1` to disable). Authenticating with an HF token still
  avoids the limit entirely.

## [2.5.2] - 2026-09-09

### Added
- **Two-pass fast mode (`complete-search --fast`, opt-in)** for dense / aggregate variant
  panels where per-haplotype enumeration is intractable (measured **49 h+ without
  completing** on a 4×-density 1000G+HGDP panel). Instead of enumerating the 2^k IUPAC
  haplotype lattice per window, it emits a small set of **worst-possible representatives**
  (reference + minimum-edit + maximum-CFD), propagated to the whole post-analysis via
  `CRISPRME_FAST_MODE`. Validated **lossless for locus detection** and **non-understating
  for the worst-case score** against the slow full-enumeration path on a real chr22
  1000G-2021+HGDP slice (0 CFD under-reports; surfaces *stronger* worst cases at 182 loci),
  and **confirmed genome-wide** (V2 non-`--fast` vs V3 `--fast` on the 2021 panel: V3 is a
  locus-level superset, **0 of V2's 1,458 CFD≥0.2 loci lost or demoted**, CFD exact-or-
  conservative, the ML scorer screen-grade only near the 0.2 line). The default (non-`--fast`) path
  is **byte-identical**. See `docs/DESIGN_2.5.1_two_pass_fast_mode.md` and METHODS §5/§8.
- **All-source "mega" sites-only index (5 datasets).** A new merged panel — 1000 Genomes
  2021, HGDP, gnomAD v4.1, TOPMed, All-of-Us — built directly from each source's aggregate
  `INFO/AF` (these carry no shared samples to reconstruct cross-source haplotypes), with
  per-dataset AF plus a per-site **`AF_max`** (the maximum across sources) as the reported
  global frequency. `merge_mega_sites.sh` + `build_mega_registry.py` compile the Tier-0
  registry directly from the frequencies (byte-identical binary; SNPs, with an indel-AF
  companion store), and the enricher tolerates sites-only VCFs. See METHODS §2/§3.
- **Report companions for aggregate panels.** The report now bundles the per-dataset
  **`indel_af.tsv`** companion (indel allele frequencies by source) alongside the SNP+indel
  co-occurrence TSV, and renders carrier / homozygote frequencies as **"NA"** for aggregate
  (sites-only) groups — where only allele frequencies exist and per-individual carrier counts
  are undefined — instead of fabricating numbers.
- **`download` always ships the shared reference index** alongside a variant index, so a
  reference-genome scan (or an on-demand variant-index rebuild) works immediately after a
  download with no separate build step.
- **SNP+SNP co-occurrence companion (`<output>.snp_snp_cooc.tsv`).** The SNP-side analogue
  of the SNP+indel co-occurrence companion: one joinable row per emitted variant off-target
  that USES ≥2 co-occurring SNP alt alleles, recorded from BOTH the observed-haplotype
  enumerator (genotyped panels → `CONFIRMED`/`PUTATIVE`, exact phase) and the registry-only /
  capped finalizer (sites-only panels → `PUTATIVE`). Reports the phase, a conservative
  `MinAF_bound` (the minimum participating marginal AF — a valid upper bound on the joint cis
  AF for both phases), and the observed carriers (`NA` on the sites-only path). Gated on a
  variant registry being present, so a legacy dict install is byte- and allocation-identical.
- **SNP+indel PUTATIVE co-occurrence on sites-only panels (mega).** The SNP+indel
  co-occurrence emission previously wrote a row only when cis carriers were countable, so on
  an aggregate sites-only panel (no per-sample genotypes) an indel + used SNP alt that
  co-occur at a locus produced no row. It now emits a `PUTATIVE` row with the conservative
  `min`-AF joint bound (over the indel `AF_max` + each used SNP's registry AF) and `NA`
  carrier fields. Gated on genotype ABSENCE, so a genotyped panel's real no-cis-carrier call
  is never turned into a phantom row.
- **Lossless dense-region worst-case haplotype (`CRISPRME_LOSSLESS_DENSE`, default ON for
  sites-only).** In a capped dense window the min-mismatch representative can be a strict subset
  of a genuine carried haplotype, so a sites-only off-target requiring ≥4 co-occurring variants
  on one haplotype could be dropped. The registry-only (sites-only) path therefore also emits
  the full co-located variant union — the maximal PUTATIVE haplotype — as an extra
  representative, bounded (not the 2^k lattice) and gated by the finalizer's own mm/PAM budget
  so nothing over-budget or PAM-invalid is emitted. The effect is **scoped to sites-only**
  (`registry_only_mode`), so it is **byte-identical for every genotyped / legacy / dict install**
  (those never enter the branch, and the genotyped path is already lossless via the observed
  enumerator). It is **on by default** so the sites-only panel fulfils "don't miss a region";
  set `CRISPRME_LOSSLESS_DENSE=0` to opt out (sites-only reverts to a greedy representative
  only). A per-sample union on the genotyped path is intentionally *not* done — it would risk
  phantom trans-as-cis haplotypes.

### Performance
- **ML scoring: load the model once + skip eager per-pentamer work.** The 276 MB
  RandomForest ensemble was re-`pickle.load()`ed on every scoring batch (both the SNP and
  INDEL post-analysis); it is now cached at module scope and loaded once per process.
  Separately, `get_features` no longer runs an eager per-pentamer `re.sub` + `random.choice`
  on N-free windows (guarded behind `if "N" in ...`), which also removes a latent
  nondeterminism from the hot path. Both are **byte-identical** at the reported score
  precision (validated against the shipped model); `get_features` micro-benchmarks ~1.6×
  faster. The dense-panel INDEL post-analysis remains dominated by single-threaded per-target
  ML-scorer compute (profiled scoring-bound, not enumeration-bound); parallelizing it is tracked
  separately.

### Fixed
- **Post-analysis no longer deadlocks at high thread counts.** A genome-wide variant
  post-analysis launched with a large `--thread` (e.g. 200 on a 200-core host) could hang
  indefinitely: a `multiprocessing.Pool` with hundreds of workers makes every worker inherit
  every other worker's internal sentinel pipe at fork, so the pool's join never finalizes
  (all per-chromosome tasks complete, the parent hangs in `pool.map` forever). Both the SNP
  and indel post-analysis pools now cap the worker count (default 32, override
  `CRISPRME_POSTPROC_MAX_WORKERS`) — no speed loss (post-analysis is I/O-bound, one
  subprocess per worker) and no deadlock.
- **Mega (all-source) index hardening.** Resolved 17 adversarial-review findings in the
  sites-only aggregate index path and de-duplicated colliding indel registry keys.
- **Clear pre-flight error when a variant search has no prebuilt index.** `complete-search
  --vcf` without a prebuilt index used to fall back to a legacy on-demand enrichment that
  builds a classic per-sample dict WITHOUT the dict-less tiers (registry/genotype) and
  could leave the indel coordinate map (`log_indels`) empty — crashing the variant/indel
  post-analysis ~30 min in with a cryptic `FileNotFoundError: log<chrom>.txt`. It now
  detects the missing tiers per dataset up front and stops immediately with an actionable
  message pointing to `download --what index` (recommended) or `build-index-only` first.
  Override for advanced/legacy on-demand builds with `CRISPRME_ALLOW_ONDEMAND_BUILD=1`.
- **Warn when SNP+indel co-occurrence is requested but the registry tier is absent.** The
  pre-flight guard accepts either the classic per-sample dict OR the dict-less registry, but
  co-occurrence needs the **registry** specifically. A classic-dict-only index (e.g. a local
  `build-index-only` run **without `--samplesID`**) previously ran cooc as a **silent no-op**;
  `complete-search` now prints a WARNING at launch naming the datasets missing the registry
  tier, so the missing co-occurrence output is no longer silent.
- **Fast-mode worst-case CFD is now exact on every path.** CFD is position-weighted, so the
  fewest-mismatch representative does not maximize CFD; fast mode now also emits the exact
  **maximum-CFD** representative (per-position argmax + bounded brute-force for the joint
  two-base PAM factor). On the **legacy dict / aggregate-panel (mega) path** this closes a
  **threshold-crossing** under-report (up to **0.14** CFD). Cross-checked against an
  independent factorized CFD oracle (4,000 random windows, joint-PAM case included) and a
  real CFD-scored multiallelic fixture.

### Notes
- In fast mode, **CFD is the exact worst case; the ML scorer is best-effort** (a non-factorizable
  RandomForest). On a chr22 1000G-2021+HGDP slice every ML-scorer ≥ 0.2 off-target was reported
  at full strength (under-reporting ≤ 0.04, sub-0.19 tail). **Genome-wide the ML-scorer tail is
  heavier:** ~5 % of ML-scorer ≥ 0.2 loci can drop below 0.2 under `--fast` (largest gap ~0.12),
  while **CFD had zero ≥ 0.2 losses** (V2-vs-V3 GW matrix). So `--fast` CFD is a safe
  actionable gate but **the ML scorer is a screen** — run **without** `--fast` for a guaranteed
  per-haplotype ML-scorer worst case / an ML-scorer action gate.
- **`--fast` scope + index/version compatibility.** `--fast` accelerates only the SNP
  post-analysis; the indel post-analysis is single-threaded and **unaffected** by it (and the
  `indel_snp_cooc.tsv` companion is byte-identical with/without `--fast`). The feature-on 2021
  index requires CRISPRme **≥ 2.5.0** — stock 2.4.0's indel post-analysis fails on it (2.4.0
  users use the 2019 / `-dictless` indices). The **mega** index currently searches **SNPs
  only** (its fake-indel genome is not built); searchable indels + PUTATIVE
  co-occurrence-without-genotypes (SNP+SNP and SNP+indel possible haplotypes from AF) are
  planned for **2.5.2**.
- **Multi-indel residual (documented):** the indel search materializes one indel per fake
  contig, so an off-target needing **≥ 2 co-occurring cis indels in one protospacer** is not
  a candidate (pre-existing single-indel-search property). Low-frequency: **~1–2%** of indel
  loci after repeat-masking/dedup; genuinely-missed off-targets **~0.1–0.2%**, all at the
  weakest (≈0-CFD) edit-budget ceiling.
- **`--fast` does NOT change SNP+indel co-occurrence.** `--fast` collapses only the *SNP*
  representative emission in `integrated_results.tsv`; the `indel_snp_cooc.tsv` companion is
  produced by the indel post-analysis' cis-phasing pass over the genotype tiers, which `--fast`
  does not touch. Measured on the completed genome-wide matrix (2021 panel, `--fast` vs
  non-`--fast`, guide TGCTTGGTCGGCACTGATAG): the two `indel_snp_cooc.tsv` files are
  **byte-identical** (same MD5, 2,729 rows, 843 CONFIRMED / 1,886 PUTATIVE, full per-sample
  `cis_samples` + joint-AF in both). Per-sample cis attribution is preserved in fast runs.
  (Corrects an earlier note that claimed `--fast` collapsed cooc to representatives only.)

## [2.5.0] - 2026-09-01

### Added
- **SNP+indel co-occurring off-targets — now enabled by default** (opt out with
  `CRISPRME_INDEL_SNP=0`; previously experimental/opt-in). Detects off-targets that
  require BOTH a nearby SNP **and** an indel on the same haplotype within a
  protospacer window — a class the classic two-pass search could not see. The build
  compiles a per-chromosome **phased indel genotype store** and overlays SNP IUPAC
  codes onto the fake-indel genome, so the indel search surfaces SNP+indel
  haplotypes; post-analysis reports **CONFIRMED-cis** (phased) and **PUTATIVE**
  (unphased) co-occurrences with per-sample carriers + joint allele frequency, and
  the shareable report bundles `indel_snp_cooc.tsv` and surfaces a
  "SNP + indel cis co-occurrences (N confirmed-cis)" section.
- **New 1000G-2021 + HGDP NRG variant index** published to HuggingFace
  (`NRG_3_hg38+hg38_1000G2021_HGDP`): 3,202 high-coverage phased 1000G samples +
  929 HGDP. `download --what index` yields a searchable, variant-ready dict-less
  install (registry + genotype + phased indel-genotype tiers bundled).
- Shareable report: publication-grade inline **vector SVG** figures.

### Changed
- **Dict-less observed-haplotype enumerator** replaces the legacy 2^k IUPAC lattice
  + greedy cap on the variant post-analysis path: every haplotype real individuals
  carry is emitted (sensitivity-first, **no truncation**), and dense/hypervariable
  windows are logged to `high_variant_density_regions.bed` for visibility rather
  than dropped — resolving the phased multi-SNP haplotype under-report on the
  default path. (Exhaustive worst-case coverage of *unobserved* haplotypes in dense
  windows remains a documented, flagged limitation; a lossless two-pass fix is planned.)
- `--max-total-edits` pruning now emits a WARNING and persists a `Pruning_note`
  into the run params so the "deep off-targets were pruned" signal reaches the
  report bundle.

### Fixed
- **Atomic per-chromosome build writes + integrity-checked resume** for the Tier-0
  registry, Tier-1 genotype, and phased indel-genotype stores: a SIGKILL/OOM
  mid-write no longer leaves a truncated file that a restart accepts as "done" (a
  silent-truncation class that had shipped a short chr1 indel-genotype store).
  Writers use tmp + `os.replace`; resume verifies integrity, not mere presence.
- Report: high-variant-density region count is **deduplicated** (was inflated by
  duplicate per-chromosome rows), and the co-occurrence TSV is deduplicated (the
  writer emits a row per alignment pass). Population-summary companion **degrades
  gracefully** (emits an NA-frequency row) instead of dropping multi-variant rows
  on a registry-only install.
- Pre-flight guard + build-time auto-fix so an over-listing `--samplesID` can no
  longer inflate the AN denominator (which deflated every reported allele frequency).

## [2.4.0] - 2026-08-20

### Added
- **Dict-less variant index is the default distribution.** The shipped SpCas9
  (NRG) variant index (`NRG_3_hg38-dictless+hg38_1000G_HGDP`, 1000G + HGDP)
  replaces the 152 GB per-sample SNP dictionaries with the compact Tier-0
  allele registry + Tier-1 genotype store — a ~16 GB download that searches out
  of the box with correct allele frequencies, rsIDs, per-individual Samples, and
  population summaries. See `METHODS.md`.
- **Shareable off-target report** (`generate-report`, auto-run at the end of
  `complete-search`): a single self-contained, branded `report.html` with an
  annotation / COSMIC legend plus a top-N table for validation. The web results
  page now shows a one-click **Download Report (.zip)** button in place of the
  three raw result-file links, and the CLI prints an "open the report" hint.
- **COSMIC cancer-gene annotation on by default.** The built-in functional-region
  bundle is now `dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz` (ENCODE SCREEN v4
  cCREs + COSMIC Cancer Gene Census tiers); results carry an `Annotation_COSMIC`
  column. Legacy bundles remain supported via automatic fallback.
- **Observed-haplotype enumerator** (dict-less multi-variant path): multi-variant
  off-targets are enumerated only for haplotypes actually carried by sampled
  individuals (CONFIRMED = phased, PUTATIVE = unphased / mixed), removing phantom
  combinations. Gated on the Tier-1 genotype store.
- `download --what all` on hg38 also fetches the default reference index
  (`NRG_3_hg38`) so a reference-genome scan works with no on-demand index build.
- **Turn-key CLI search after a download.** `download --what index` now writes the
  `list_vcf.txt` / `list_samplesID.txt` list files (the ones the web form builds
  per-search) at the install root, so `complete-search --vcf list_vcf.txt
  --samplesID list_samplesID.txt` works out of the box — no hand-written config.
  A "Run a search from the command line" section is added to the Docker
  quickstart, README and the landing site, with the `generate-report` step.
- Manuscript-level `METHODS.md` (dict-less model, VCF homogenization / merging,
  AF estimation, phased / unphased / mixed haplotype scanning, variant-density
  control, annotations, report) and README "what's new" pointers.
- **Perfect-match (0-mismatch) multiplicity flagging.** A guide with more than one
  perfect genomic match (0 mismatch, 0 bulge) has no a-priori on-target -- each is
  an equally-efficient candidate cut site. The report **and** the web results page
  now force EVERY 0-mm site to the top of the validation panel (never truncated),
  add a `Perfect_match` column to the curated tables, and show a prominent red
  "no unambiguous on-target" banner listing the sites (an amber "presumed
  on-target" note when there is exactly one). Previously all 0-mm sites were
  bucketed as "the on-target" and excluded from the off-target set and the panel
  with no warning -- so a perfect-match off-target, the highest-risk class, was the
  one thing the report never surfaced.

### Fixed
- **Web job launch crashed on special characters in the job name.** A job name
  with shell / path metacharacters (e.g. `SBDS(T>C)`) was interpolated raw into a
  `mkdir` run through `/bin/sh`, aborting the launch (`Syntax error: "("
  unexpected`, HTTP 500 in `change_url`). The name is now sanitized to a
  filesystem/shell-safe charset, and the results directory + queue file are
  created without a shell (`os.makedirs` / `open`).
- **Combined allele-frequency denominator in the shipped dict-less index.** The
  registry's global panel size was built from an over-listed sample list
  (AN = 8858) rather than the VCF-resident panel (1000G 2548 + HGDP 929 = 3477 →
  AN = 6954), deflating combined MAF by ~27 %. The published registry now reports
  the correct combined AF (AC / 6954); allele counts were already correct. A
  singleton reads 1 / 6954 = 1.44e-4.
- **Annotation sorting is read-only-install + concurrent-job safe.**
  `_sort_annotation` / `_process_personal_annotation` build the sorted / merged
  annotation in a per-invocation temp dir instead of decompressing and
  overwriting files in the shared install `Annotations/` directory.
- **Dict-less unphased multi-variant off-targets: enumerate sub-combinations.**
  The observed-haplotype enumerator reported only the maximal-*union* variant set
  for an unphased individual, which hides a valid off-target when a variant
  *breaks* the target (e.g. one that disrupts the PAM) — the sub-combination that
  keeps the reference allele at that column (a real off-target under one possible
  phasing) was dropped. Unphased (and cross-phase-set) individuals now enumerate
  every non-empty subset of their own carried variants (never cross-individual
  chimeras); the mismatch/PAM gates prune out-of-budget subsets and the enumerator
  dedups shared subsets, bounded by `CRISPRME_IUPAC_CAP`. Recovers the missed
  off-targets with correct MAF + samples, with zero change to previously-reported
  rows or to the phased (confirmed-cis) path.
- **Registry-only (`--no-genotypes`) installs no longer drop variant off-target
  sites.** With a Tier-0 registry but no per-sample dict and no genotype tier, the
  variant finalizer's "≥1 carrier" guard silently discarded every variant-created
  off-target. Such sites are now emitted as degraded rows (Samples = `NA`) that
  still carry the creating-variant's registry rsID / AF, so the risk is surfaced
  without per-sample resolution. Gated on the registry-only mode, so legacy,
  registry+dict, and dict-less-with-genotypes installs are byte-identical.
- **`complete-test` sample-ID download.** The github fallback URL was built with a
  double slash (`…/data//samplesIDs/…`, a redirect the downloader could stumble
  on); it now normalizes to a single slash.
- **`download` prefers the corrected dict-less variant index.** Requesting a
  variant index by its bare (legacy dict-based) name now transparently resolves to
  the `-dictless` companion when one exists on the repo — both install to the same
  folder, so old links / muscle-memory get the smaller, combined-AF index.
- **Annotation sorting survives a small or full `/tmp`.** `_sort_annotation` /
  `_process_personal_annotation` now stage intermediates on a roomy writable disk
  (`$CRISPRME_TMPDIR`, else the working/output dir) instead of the default `/tmp`,
  which is tiny on many HPC/container setups.
- **Dict-less never misses a high-variant-density region.** In a capped dense
  window the observed path emitted only the real observed haplotypes (per-individual
  unions) and logged the region to `*.high_variant_density_regions.bed`; when those
  unions all exceeded the mismatch budget the region surfaced ONLY in the BED with no
  off-target ROW in the results (the dict path always reports a greedy representative
  there). The observed path now also emits the same greedy min-mismatch representative
  (argmin over the window's alleles) per capped region, so every dense region always
  surfaces ≥1 off-target row at parity with the dict path (verified genome-wide).
- **Dense-region transparency.** A new `High_variant_density_region` column in
  `integrated_results.tsv` flags each off-target that falls in a dense window,
  noting that a greedy worst-case alignment is reported and additional haplotype
  alignments may exist, and gives the FULL IUPAC protospacer (every variant column
  as its ref+alt ambiguity code, also added to the `.bed`) so a user can see all
  alleles in the window and dig into the other possible alignments.
- **Deterministic off-target representative.** When two alignments of a cluster
  tied exactly on CFD and mismatches+bulges, the reported "highest-CFD
  representative" fell to enumeration order, so the dict and dict-less paths (which
  enumerate a locus's alleles in different orders) could show a different
  representative *coordinate* for a small fraction of clusters (identical cluster,
  carriers and CFD). The greedy dense-region representative and the cluster
  representative selector now break exact ties by an intrinsic key (smaller alt
  allele; aligned-DNA + genomic position), so the representative is
  order-independent — a no-op on every already-unique cluster.
- **Dict-less correctly omits dict phantom off-targets (verified).** A genome-wide
  dict-vs-dict-less comparison leaves a handful of variant off-targets that the
  legacy dict path reports and the dict-less path does not. These are NOT missed
  regions: they are the dict's cross-variant *chimeras* — an off-target attributed
  to carriers who, in their actual (phased) haplotype, also carry an additional
  variant in the window that breaks the alignment, so no real individual carries the
  off-target. The dict-less observed enumerator (which reconstructs each individual's
  real haplotype and every sub-combination) correctly emits every off-target a real
  haplotype creates and drops these phantoms — it is *more* accurate than the dict
  here, not less (see §3/§4 of `METHODS.md`, and #118/#120).
- Report polish: combined-panel MAF is spelled out in the legend, the Section-3
  population bar no longer counts the on-target, duplicate reference rows are
  collapsed, and the section-numbering docstring matches the emitted sections.
- Doc/CLI accuracy: the Docker quick-start + web-interface guide use the v2.4.0
  image and the `-dictless` index name; `--max-total-edits` help states its real
  default (4); the web "Load Example" picks the installed default PAM; and a note
  documents pointing `HF_HOME` at a large disk on small-home systems.

## [2.3.3] - 2026-08-18

### Fixed
- **Bulge dropdowns could not be raised above 0 on a downloaded dict-less index.**
  The search form caps the selectable bulges at `min(index_max_bulges(reference),
  index_max_bulges(variant))`. A dict-less variant tarball ships the variant index
  but not the reference index, so the reference term was read as `0` and forced the
  cap — and therefore the bulge dropdowns — to `0`. The reference term is now
  **buildable-aware** (`reference_bulge_capacity`): the reference index is built
  on demand from the downloaded `Genomes/<ref>`, so its capacity is taken from the
  resolved index's `<PAM>_<N>` prefix instead of an absent folder. The default
  SpCas9 install now correctly offers up to 2 DNA / 2 RNA bulges.
- **Default search parameters restored to MM 6 / bDNA 2 / bRNA 2 with max-total-edits 4.**
  These are the intended out-of-the-box SpCas9 defaults (mismatches 6, both bulges 2,
  capped at 4 total edits); a prior change had drifted them.
- **A 0-bulge variant search failed on a dict-less / batteries install.** With
  `bDNA=0, bRNA=0` the indel branch ran a non-indexed CRISPRitz search of the
  **raw** fake-indel genome (`Genomes/<ref>+<vcf>_INDELS/`), which a dict-less
  install does not ship (build-time only) — producing missing `fakechr*.targets.txt`
  files and a `tail/head/rm/sed` cascade in `post_analisi_indel.sh`. The 0-bulge
  indel search now routes through the **indexed** indel genome
  (`genome_library/<idx>_INDELS/`, which the tarball does ship), matching the
  ≥1-bulge path. No behavior change for dict-based installs.
- **Population-distribution post-analysis crashed on out-of-range totals**
  (`max() arg is an empty sequence`). `populations_distribution.py` /
  `process_summaries.py` now guard the empty case (`max(..., default=0)`, an empty
  "no targets" figure) and the driver loop bound is computed from `mm + bDNA + bRNA`
  (the true edit budget) instead of the index `N`, so no row index falls outside the
  populated range.

## [2.3.2] - 2026-08-18

### Fixed
- **Self-built dict-less indexes are now self-complete for samplesIDs.**
  `build-index-only` (with `--samplesID`) writes the per-dataset and the combined
  `<vcf>.samplesID.txt` into the install (reusing the one union implementation), and
  `publish-index` bundles them into the index tarball — so a fresh `download
  --index-name <X>` + search works with no separate `--what samples`/`--what all`
  dependency (the download-side synthesis added in 2.3.1 becomes a pure fallback).
  No-op for single-dataset, reference, and dict-based indexes.

### Documentation
- **Rewrote the index-build docs for the dict-less flow.** `docs/PRECOMPUTED_INDEXES.md`
  now covers `build-index-only --vcf --samplesID` (**required** to emit the Tier-0
  registry + Tier-1 genotype tiers — without it you get a dicts-only index),
  `publish-index --dictless`, the separate `genotypes_<vcf>.tar.gz`, `download
  --no-genotypes`, the `-dictless` canonical install-name behavior, and the updated
  `manifest.json` fields (with a caveat that name-parsed fields are reference-index
  only). Added the variant-index build path to README §2.2.9, `crisprme_data_setup`
  §3.5, and `SETTINGS_DATA_MANAGER.md`; fixed an NGG→NRG example nit.

## [2.3.1] - 2026-08-18

### Fixed
- **Dict-less `download --what index` is now self-sufficient.** A standalone
  `download --what index` for a merged variant index (e.g.
  `NRG_3_hg38-dictless+hg38_1000G_HGDP`) now fetches the per-dataset `samplesID`
  lists it needs before synthesizing the combined `<vcf>.samplesID.txt`, so it no
  longer requires a prior `download --what all`/`--what samples`. No-op when those
  files are already present (the `--what all` → `--what index` path is unchanged).
- **Restored the `unit-tests` CI gate.** A colon in a workflow step name
  (`(#172: fakechrom match)`) made GitHub Actions reject `unit-tests.yml` at startup
  (a 0 s "workflow file issue"), so the unit-tests job had silently not run since the
  v2.3.0 merge (the tests themselves pass locally). Quoting the name fixes the parse;
  the full 22-step unit suite runs in CI again.

## [2.3.0] - 2026-08-18

The dict-less variant-analysis engine: replaces the ~152 GB per-sample SNP
dictionaries with compact, memory-mapped tiers, adds population-level off-target
summaries, and corrects the allele-frequency column. Backward-compatible — a
dict-based install still works unchanged.

### Added
- **Compact dict-less variant post-analysis.** A tiny, always-shipped **Tier-0
  registry** (per-(pos,alt) AC/AN + per-(db×subpop)/global allele counts + rsID,
  mmap+bisect) powers off-target detection and corrected allele frequencies out of
  the box; an optional **Tier-1 genotype store** reconstructs the exact per-sample
  Samples column. Together they replace the 152 GB per-sample dicts with ~7 GB
  (registry) + ~22 GB (genotype tier), random-access — SNP post-analysis no longer
  streams 152 GB of per-sample JSON. `build-index-only` emits the tiers; the search
  auto-detects and uses them (falling back to dicts when present).
- **Combination-aware population-summary output** — a `<output>.population_summary.tsv`
  companion per variant off-target with per-database, per-superpopulation and global
  allele/carrier frequencies, max-subpopulation AF (+label), homozygote counts and
  absolute Ns, with **dataset provenance preserved** (never conflates 1000G vs HGDP).
  Phased datasets use exact cis co-occurrence; unphased report an assume-cis upper
  bound + a labeled lower bound.
- **Dict-less HuggingFace distribution.** `publish-index --dictless` ships the
  registry in the main index tarball + the genotype tier as a separate optional
  `genotypes_<vcf>.tar.gz`; `download` fetches the genotype tier by default
  (`--no-genotypes` to skip). Published `NRG_3_hg38-dictless+hg38_1000G_HGDP` on
  `lucapinello/crisprme-data` (the dict-based `NRG_3`/`NGG_3` indexes remain the
  default). Fetch it with `download --what all` first, then `download --what index`,
  so the per-dataset `samplesID` lists are on disk when the combined
  `<vcf>.samplesID.txt` is auto-synthesized (a standalone `--what index` without a
  prior `samples` fetch leaves the combined file absent — hardening tracked as a
  follow-up).

### Fixed
- **Corrected allele frequencies.** The AF column — empty/mis-polarized for ~95 % of
  variants in the dict format (a documented 2.2.0 limitation) — is re-derived from
  AC/AN over the full panel.
- **Out-of-the-box dict-less install.** The skip-enrichment gates
  (`submit_job_automated_new_multiple_vcfs.sh`, `validate_inputs.py`) now accept a
  `registry_<vcf>/` tier in place of `dictionaries_<vcf>/`; `download` installs a
  `-dictless`-named index under its search-resolvable canonical name and synthesizes
  the combined `<vcf>.samplesID.txt` from the per-db files. Verified end-to-end: a
  fresh dict-less download finds the CPS1 off-target (chr2:210530658, CFD 0.947,
  rs114518452) with Samples reconstructed from the genotype tier and a populated
  population summary.

## [2.2.0] - 2026-08-18

First stable release of the CRISPRme+ 2.2.0 line (Python 3.11 / Dash 2.x). It
consolidates the `2.2.0-alpha.27`–`alpha.30` pre-releases (see those sections
below for the full history); the entries here are the changes since `alpha.30`.

### Fixed
- **Phased multi-SNP haplotype off-targets are no longer under-reported.** In the
  phased post-analysis path — the default for phased datasets such as 1000G/HGDP —
  `iupac_decomposition` failed to assemble the full-haplotype off-target when one
  sample carried **≥4 co-occurring variants in a single protospacer window**: the
  true worst target (fully substituted, lowest-mismatch/highest-CFD) was silently
  under-reported or dropped. The in-loop level-0 subtraction (a dedup device that
  starved the deeper combination layers) is replaced with a **deferred peel** that
  runs after the full lattice is built, so the maximal *cis* combination forms while
  per-haplotype attribution stays deduplicated. Validated byte-identical for every
  previously-correct case (unphased, phased ≤3-variant, single-SNP, >cap greedy)
  plus a real-data genome-wide no-regression diff. (#41; a long-standing defect
  inherited from classic CRISPRme, tracked there at pinellolab/CRISPRme#175.)

### Changed
- `environment.yml` now includes `ijson`, so from-source installs get the streaming
  low-RAM SNP-dict reader that the Docker image already had (without it, SNP
  post-analysis falls back to a whole-file `json.load`; results are identical, but
  RAM use and the OOM-guard estimate are higher).

## [2.2.0-alpha.30] - 2026-08-17

### Fixed
- **Indel off-targets were silently dropped since v2.1.9.** `post_analisi_indel.sh`
  subset the per-chromosome indel targets with `grep -F -w $chrom`, but the indel
  search runs on a per-chromosome *fake* genome (`pool_search_indels.py` searches and
  names its targets `fake<chrom>`), so every indel target row's Chromosome column is
  `fake<chrom>` (e.g. `fakechr22`) — which `grep -F -w chr22` can never match (`-w`'s
  left word boundary fails because `chr22` is preceded by the word character `e`). The
  per-chromosome subsets came back empty, the indel post-analysis processed nothing,
  and **all indel off-targets were dropped with a clean exit**. Both grep lines now
  match on `"$fakechrom"` (`-w` still prevents `fakechr2` matching `fakechr22`); the
  `NF >= 10` malformed-line guard that the same refactor had dropped is restored; and
  an empty subset now prints a WARNING to stdout (stderr is fatal in this pipeline) so
  the failure can never be silent again. A hermetic regression test
  (`test_indel_chrom_subset.py`) is added and wired into CI. SNP post-analysis is
  unaffected (SNP targets use real chromosome names). Thanks to @munchr-gene1 for the
  report and diagnosis (#172).

## [2.2.0-alpha.29] - 2026-08-16

### Fixed
- **Variant search without `--annotation` no longer crashes on a read-only install
  (apptainer SIF).** `post_process.sh` wrote the annotation intermediate BEDs
  (`<annot>.tmp.bed` / `.tmp.sorted.bed`) next to the annotation file; with no
  `--annotation`/`--gene_annotation` that file is the shipped empty placeholder
  `PostProcess/vuoto.txt` inside the (read-only) install dir, so the writes failed
  and the run died at "Integrating results". Temps now go to the run's writable
  output dir. As defense-in-depth, `resultIntegrator.py` guards the closest-gene
  distance parse so a missing/empty annotation field is treated as "no annotation"
  (columns stay `NA`) instead of raising `ValueError: could not convert string to
  float: ''`. The Docker/web path (real annotation, writable FS) is byte-identical.

### Added
- **Pre-flight check for unpadded guides.** A CLI `--guide` file whose guides lack
  the PAM's N-padding (e.g. the bare 20 nt `CTAACAGTTGCTTTTATCAC` instead of
  `CTAACAGTTGCTTTTATCACNNN`) previously ran the full multi-hour search then crashed
  deep in `radar_chart_dict_generator.py` (`IndexError`). The lightweight validator
  now catches this up front with an exact fix message (append/prepend N's), handling
  both 3' (SpCas9) and 5' (Cas12a) PAM orientations. The web is unaffected (it
  auto-pads guides). Covered by 10 new `check_guide_pam_consistency` tests; the whole
  `test_validate_inputs` suite is now run in the unit-tests CI (previously it wasn't).

## [2.2.0-alpha.28] - 2026-08-16

### Fixed
- **Batteries variant search no longer aborts at pre-flight validation.** The
  lightweight input validator (`validate_inputs.run_lightweight`, always-on before
  `complete-search`) hard-required `VCFs/<dataset>/` to exist and aborted otherwise —
  but a batteries-included install (`download --what index`) ships the precomputed
  variant index + per-sample dictionaries *without* the multi-GB source VCFs. The
  validator now recognizes that case: when the VCF dataset directory is absent but the
  dictionaries (`Dictionaries/dictionaries_<dataset>/` + `log_indels_<dataset>/`) are
  present, it passes the check (source VCFs genuinely not required), mirroring the
  skip-enrichment logic in `submit_job_automated_new_multiple_vcfs.sh`. A genuinely
  missing dataset (no dictionaries either) still fails loudly. Completes the alpha.27
  batteries fix: `download --index-name NRG_3_hg38+hg38_1000G_HGDP` + a variant search
  now runs end-to-end on a fresh install. Covered by two new `run_lightweight` tests.

## [2.2.0-alpha.27] - 2026-08-16

### Fixed
- **Batteries variant search is now self-sufficient (the index-only download can run
  a variant search).** The variant post-analysis needs per-sample SNP dictionaries
  (~152 GB for 1000G+HGDP) that weren't shipped, so a fresh `download --index-name
  <variant>` install failed at genome enrichment (no source VCFs). The dicts are now
  shipped **gzipped** (152 GB → ~28 GB; indel logs 18 GB → 3.3 GB) and read **on the
  fly**, mirroring the existing bgzipped-VCF pattern — no 150 GB on-disk decompress.
  - `new_simple_analysis.py` / `analisi_indels_NNN.py`: prefer `*.json.gz` / `*.txt.gz`
    (ijson streams through gzip; json.load + line reads via `gzip.open`), fall back to
    plain — backward-compatible. Bash wrappers unchanged (Python resolves the extension).
  - `pool_search_indels.py` / `pool_post_analisi_indel.py` + the `submit_job` batteries
    chrom glob: strip an optional `.gz` from `log*.txt(.gz)`.
  - `build-index-only` gzips the per-chromosome dicts it produces (pigz-preferred).
  - Verified end-to-end by `validate-benchmarks` (build gzips → search reads `.gz`) and
    a new `test_dict_gzip` unit test (publish bundles + download routes `.gz` dicts).

## [2.2.0-alpha.26] - 2026-08-16

### Changed
- **Default SpCas9 PAM is now NRG (NAG + NGG), not NGG.** The web form and the
  "Load Example" default, plus the CLI/download docs, now default to `NRG`
  (`20bp-NRG-SpCas9`) and the `NRG_3_*` indexes. NRG matches SpCas9's broad
  recognition, so **variant-created NAG off-targets are found out of the box** —
  e.g. the CPS1 off-target from the CRISPRme paper (rs114518452: chr2:210530658,
  2 mismatches + 1 RNA bulge, variant-created `TGG` PAM, CFD 0.947). Verified
  end-to-end on a from-scratch `NRG_3_hg38+hg38_1000G_HGDP` index. `_default_pam`
  prefers NRG, falls back to NGG, then the first available PAM. The `NGG_3_*`
  indexes remain available for canonical-NGG-only searches. New `NRG_3_*` indexes
  are published on `lucapinello/crisprme-data`.

### Fixed
- **Off-target search now uses all requested cores.** `submit_job` split
  `$ncpus` in half and ran reference + variant searches in parallel; the small
  reference genome finishes first, so the long variant-search tail ran on only
  half the cores. Now runs the two searches sequentially at full `$ncpus`
  (reference 44.6 → 25.6 min in testing; `searchTST` 24 → ~45 cores).
- **Freshly-built / downloaded indexes now search correctly.** `build-index-only`
  writes `.pam_build`/`.display_label` sidecars inside the index dir, but
  CRISPRitz's `-index` search requires the dir to contain only `.bin` files
  ("only .bin files" error). Added a symlink-only `.bin` view so the sidecars no
  longer break the search, regardless of how the index was built or downloaded.

## [2.2.0-alpha.25] - 2026-08-15

### Changed
- **Contacts page tidied.** Developers are now listed with the Pinello lab first
  (Luca Pinello), then the Verona/InfOmics group (Manuel Tognon et al.), then the
  Bauer lab (Linda Lin, Daniel Bauer). Comments are directed to Luca Pinello only
  (single email), and the page now says bugs should be filed as a GitHub issue.
  Also fixed a missing space after the comma in the per-lab author lists.

## [2.2.0-alpha.24] - 2026-08-15

### Changed
- **Submission page: consolidated to exactly three font sizes (Large / Medium /
  Small).** The form previously mixed ~7 inline sizes on a 10px root, which left
  helper notes rendering *larger* than the inputs they described. It now uses one
  clear hierarchy: Large `2.2rem` for the numbered section headers ("Select
  gRNA/genome/PAM/thresholds"), Medium `1.5rem` for every interactive control
  (spacer textarea, dropdowns, Load Example / Submit buttons, the "Maximum edits"
  label and form body), and Small `1.25rem` for helper text, notes and slider
  tick marks. Verified via full-page + high-DPI screenshots — no clipping or
  overflow. Cosmetic only; no change to search behavior.
- **"Max edits" slider maximum is now derived, not a magic number.** The slider
  topped out at a hardcoded 5; it now tops out at the app's mismatch ceiling
  (`MAX_MMS - 1` = 6). Rationale: mismatches are the non-index-limited edit
  dimension and the dominant search cost, so the total-edits knob should reach as
  far as the mismatch dropdown does. The installed index caps *bulges* only
  (`bMax = N - 1`, ≤2 in practice) — that limit is applied separately where
  simple-mode bulges are derived and must not shrink the mismatch-driven total.
  Default stays 4.

## [2.2.0-alpha.23] - 2026-08-15

### Changed
- **"Max edits" is now a self-sufficient single knob — bulges are derived from it.**
  Previously a variant off-target that needs a bulge (e.g. the CRISPRme paper's CPS1
  off-target rs114518452: 2 mismatches + 1 RNA bulge) was only found if bulges were
  explicitly enabled — the CLI defaulted to bulge-free, so `--max-total-edits` alone
  found no bulge targets. Now:
  - **CLI:** when `--bDNA/--bRNA` are omitted, they are derived from `--max-total-edits`,
    capped by the bulge depth the installed index supports (so `complete-search --vcf …
    --max-total-edits 4` finds up to 2mm+2bulge patterns with no bulge flags, reusing the
    installed index — never forcing a heavier build; explicit `--bDNA/--bRNA` still win).
  - **Web:** in simple mode the per-type bulge caps now track `min(max_edits, index_cap)`
    instead of a fixed 2/2, so the slider is genuinely the only control needed.
- **Consistent default cap of 4** across web slider, Load Example, CLI `--max-total-edits`,
  and the shell fallback (was 5). Four surfaces the CPS1 off-target (needs 3 total edits)
  with headroom and is the validated performance sweet spot (cap 6 was pathologically slow
  on genome-wide combined searches).
- Docs: `--bDNA/--bRNA` are documented as optional (derived from `--max-total-edits`), not
  "Required"; the variant-aware requirement for variant-created off-targets is noted.

## [2.2.0-alpha.22] - 2026-08-15

### Added
- **History page now reports the Max edits used.** The results-history table only showed
  Mismatches / DNA bulge / RNA bulge; it now has a "Max edits" column showing the total
  mismatches+bulges cap (the simple-mode slider value, or the per-type sum with an
  "(advanced)" note). Jobs predating this field show "-".

## [2.2.0-alpha.21] - 2026-08-15

### Fixed
- **Graphical Reports: every variant off-target now shows its blue "Alternative" point.**
  The top-1000 scatter drew a blue point only when a variant had a known allele frequency;
  a variant whose MAF was missing (e.g. the frequency is not annotated from the external
  AF source) was treated like a reference-only site and its blue point was suppressed. On
  the combined 1000G+HGDP panel that hid ~93% of the variant points, so the plot looked
  "all red." Now a blue point is drawn for any off-target carried by >= 1 sample, sized by
  allele frequency when known and a default size when not; only true reference-only sites
  get no blue point. This reproduces the CRISPRme Nature Genetics figure — e.g. the CPS1
  off-target rs114518452 (chr2:210530658, CFD 0.021 -> 0.947, a variant-created off-target)
  appears as the top blue point. Also fixed in the per-sample Personal Risk Card plots.

## [2.2.0-alpha.20] - 2026-08-15

### Added
- **Graphical Reports: a "top 1000 by variant effect" plot.** Alongside the existing
  score-ranked top-1000 scatter, the CFD and ML-scorer reports now also show a companion
  plot ranked by the size of the variant-induced score change (`|ALT - REF|`), so the
  variants that actually change the off-target score are foregrounded instead of being
  buried among the many that leave it unchanged. The score-ranked plot is unchanged; the
  fewest-mm+b report stays score-only. (`CRISPRme_plots.py`, shown in `results_page.py`.)

## [2.2.0-alpha.19] - 2026-08-15

### Fixed
- **Graphical Reports: the top-1000 scatter now shows the variant (ALT) points.** The blue
  ALT markers were drawn ~10x smaller than the red REF markers at the same x, so they were
  swallowed and the plot looked "all red." They now have a floored size + a black edge (and
  stay hidden for reference-only sites), so a variant that changes the CFD / ML-scorer score is
  visible at the tip of its arrow. (`CRISPRme_plots.py`, `CRISPRme_plots_personal.py`.)
- **Graphical Reports: the ML-scorer and fewest-mm+b top-1000 plots are generated again.** The
  plot generator sorted the ML-scorer plot on a non-existent column
  (`CFD_score_(highest_CRISPR_BULGE)`) and crashed with a silently-swallowed `KeyError` right after
  the CFD image, so only the CFD plot appeared. Fixed to `CRISPR_BULGE_score_(highest_CRISPR_BULGE)`.
- **Reference-only searches no longer break the Graphical Reports tab.** The population-barplot
  callback `Output` targeted a container the layout omitted for reference-only searches, which
  could fail the whole callback (radar chart included). The container is now always present
  (hidden for reference-only).
- **Personal Risk Card now shows the candidate sites.** The detail table was built only from
  *private* targets (unique to the sample), so it was empty for the common case of a sample
  with shared-but-not-private candidates. It now shows both a **Personal candidates** table
  (all off-targets the sample carries) and the **Private** table.
- `update_images_tabs` now guards an empty table selection (no crash on initial render) and
  applies the non-SpCas9 -> "fewest" criterion remap the rest of the page uses; fixed a radar
  chart `except` handler that assigned an unused variable.

## [2.2.0-alpha.18] - 2026-08-15

### Fixed
- **INDEL post-analysis: forward-ported the boundary/dense-variant crash guards.** The
  alpha.17 ML-scorer window guard was only in the SNP path (`new_simple_analysis.py`); the
  parallel INDEL path (`analisi_indels_NNN.py`) still assumed a full 29 nt window, so an
  out-of-range indel coordinate could crash post-analysis. The `!= 29` null-guard is now
  mirrored in both INDEL ML-scorer blocks.
- **INDEL CFD scoring no longer KeyErrors on an unexpected PAM.** `analisi_indels_NNN.py`
  used a raw `pam_scores[pam]` where the SNP path uses the guarded `pam_scores.get(pam,
  0.0)`; now matched (issue-#94 class).
- **CLI `--sorting-criteria` default now matches its own `--help` and the web.** It
  defaulted to `mm+bulges` but documented (and the web uses) `mm+bulges,mm`; the
  fewest/total tie-break now ranks tied targets identically on both paths.

### Changed
- **CLI and web now share the same default total-edits cap (5).** The web "Maximum edits"
  slider defaulted to 3 while the CLI `--max-total-edits` defaults to 5; the web slider
  (and Load Example) now default to 5 so a default web search and a default CLI search
  cover the same edit space.

## [2.2.0-alpha.17] - 2026-08-14

### Fixed
- **Post-analysis no longer crashes on variant-dense / chromosome-boundary targets.**
  On the combined 1000G+HGDP index the example search (and any dense-variant search)
  could abort in post-analysis with a `ZeroDivisionError` (ML scoring divided by the
  length of a 29 nt window that had been stripped empty of IUPAC codes) or an
  `IndexError` in `iupac_decomposition` (a reference window truncated at a chromosome
  boundary was shorter than the target). Both came from assuming a genome slice is always
  full length. Now an un-scoreable ML-scorer window is nulled (score −1, as for windows with
  `N`) and IUPAC positions past a truncated reference are skipped, with belt-and-suspenders
  divide-by-zero guards in the ML scorer module. Verified end-to-end on the batteries
  1000G+HGDP data: the example completes and its Personal Risk Card renders.

### Changed
- **Base editing is off by default after "Load Example."** The example now loads a
  standard (non-base-editing) search; the base-editor toggle stays on "No".
- **Larger helper text on the submission form.** The small gray guidance texts in the
  submission cards were ~10.5 px (the page root font is 10 px); enlarged to ~14 px so
  they are readable. The page's base font is unchanged.

## [2.2.0-alpha.16] - 2026-08-14

### Fixed
- **Result title is now mode-aware and no longer self-contradictory.** In simple mode
  the title listed the per-type mismatch/DNA/RNA caps *and* the total "Max edits" cap,
  which contradict when the total cap is tighter (e.g. `Mismatches 4 - DNA bulges 1 -
  RNA bulges 1 - Max edits 1`). Now a simple-mode search shows only `Max edits
  (mismatches + bulges) N` (the single control the user set), while an advanced search
  shows the explicit per-type caps. `change_url` records `Threshold_mode` (simple/
  advanced) in `.Params.txt`; older jobs are inferred (a total cap tighter than the
  per-type sum ⇒ simple) and pre-alpha.15 jobs fall back to the per-type title.
- **Dedup reuse wrote a malformed completion marker.** When a resubmission reused a
  prior clean result, the marker was written `Job\nDone` (newline) instead of
  `Job\tDone` (tab), so the dedup check and the status page's completion detection
  would not recognize the reused job as finished. Fixed to write `Job\tDone`.

## [2.2.0-alpha.15] - 2026-08-14

### Fixed
- **The result title now shows the "Max edits" cap that actually bound the search.**
  In simple mode the "Maximum edits (mismatches + bulges)" slider is the binding
  constraint, but the Result Summary title only listed the per-type mismatch/DNA/RNA
  caps — which can be looser than the total cap (and, after "Load Example", are values
  the user never explicitly chose), so the title overstated how broad the search was.
  The total-edits cap is now recorded in `.Params.txt` (`Max_total_edits`) and appended
  to the title as "Max edits (mismatches + bulges) N". Backward-compatible: results for
  jobs created before this field simply omit the clause.

### Changed
- **"Load Example" now sets the Max-edits slider (to 4) explicitly** instead of
  leaving it at the default, so the example is driven by the same total-edits cap the
  report displays and the per-type caps it fills in (4/1/1) are consistent with it.
- Docs: added mentions of the optional **email notification** (README, Docker
  quickstart, landing pages) and clarified that **Personal Risk Cards** appear only
  for variant-aware searches (landing pages, in-app help, quickstart, user guide).

## [2.2.0-alpha.14] - 2026-08-14

### Fixed
- **A resubmitted search no longer resurfaces a previously-FAILED job.** Results are
  deduplicated by search parameters; the dedup reused any matching prior job, so
  re-running a search whose earlier attempt had failed (e.g. left behind by a bug
  since fixed) returned that stale failure — surfaced in the UI as "The selected
  result encountered some errors, please remove it and try to submit again." Dedup
  now skips prior results that did not finish cleanly (non-empty `log_error.txt` or
  no `Job Done`) and runs the search fresh instead.

### Changed
- Docs: dropped the stale "peak memory is in post-processing and can exceed
  32 GB / ~64–100 GB" guidance (the alpha.12 streamed SNP-dict load removed that
  post-analysis spike); the 16 / 32 / 64 GB tiers stay.

## [2.2.0-alpha.13] - 2026-08-14

### Changed
- **Post-analysis performance** (follow-up to the streaming dict-load fix):
  - Added the `yajl` C library + `cffi` so `ijson` uses its fast `yajl2_cffi` C
    backend instead of the pure-python fallback (much faster streaming of the SNP
    dictionary).
  - The post-analysis worker cap no longer sizes per-worker RAM from the (large) SNP
    dictionary file: with streaming, per-worker peak is small, so it uses a modest
    fixed estimate and runs chromosomes in parallel again (the SNP path checks that
    `ijson` is importable; the INDEL path never loaded the big dict). Previously it
    over-estimated and serialized to a single worker.

## [2.2.0-alpha.12] - 2026-08-14

### Fixed
- **On-demand INDELS index built with a doubled `<PAM>_<N>` prefix** (regression
  from the alpha.10 indels-index-resolution change). crispritz `index-genome`
  already auto-prepends `<PAM>_<bMax>_` to the output dir, but the builder also
  added it manually, producing `NGG_2_NGG_2_..._INDELS`; detection/search look for
  the single-prefix name, so the indel search produced no output and the run failed
  with "off-targets search failed". Removed the manual prefix (let crispritz add
  it). This unblocks `validate-benchmarks` (complete-test + validate-test now match
  the brute-force ground truth: 2/2 light benchmarks). Batteries searches were
  unaffected (they use a correctly-named pre-built `_INDELS`).
- **Genome-wide variant post-analysis OOM ("Killed … EmptyDataError").** On a
  genome-wide 1000G+HGDP search, each per-chromosome worker `json.load()`ed the
  whole chromosome SNP dictionary — a chr2 dict is ~13 GB on disk and ~26 GB in
  RAM — so even one worker could exhaust a normal machine (and the concurrency
  guard defaulted to a fixed 64 GB budget / 4 GB-per-worker, over-subscribing on
  smaller hosts). `new_simple_analysis.py` now loads **only the dictionary entries
  the search's off-targets actually query** (streamed with `ijson`), dropping peak
  RAM from ~26 GB to <0.1 GB on the real chr2 dict (byte-identical results). The
  worker-cap guards now **detect the machine's real RAM** (host + cgroup limit),
  **estimate per-worker RAM from the actual dictionary sizes**, and print a clear
  "needs ~N GB / may be OOM-killed" warning instead of the cryptic cascade. Added
  `ijson` to the image.

### Changed
- Docs: corrected the in-app manual and web guide (removed the stale "up to 100
  spacers / first-100-sequences" copy; annotation is Settings-managed not a Step-3
  dropdown; email is configured in Settings, not "the server notifies").

## [2.2.0-alpha.11] - 2026-08-14

### Fixed
- **Graphical Reports population barplot never rendered.** The pipeline writes
  `populations_distribution_<guide>_<N>total_<score>_new.png`, but the results
  page read the name without the `_new` suffix, so every mismatch/bulge combo
  showed "No result found for this combination". Added the missing suffix.
- **Graphical Reports defaulted to an empty combo.** The "up to N edits" selector
  defaulted to `0` (on-target only), which is empty for typical sparse guides. It
  now defaults to the full budget (mm + bulges) — the cumulative barplot the
  pipeline only emits when it has data — so a plot shows immediately.
- **"Load Example" left the base-editing window blank.** The `be-window-start/stop`
  dropdowns start with no options, so the example's window values were dropped
  before `update_base_editing_dropdown` populated the options. Load Example now
  emits the window options atomically with the values (and passes ints, matching
  the option type), so it fills 4/8 correctly.

## [2.2.0-alpha.10] - 2026-08-14

### Added
- **Email notifications are now configurable in the app** (Settings → Email
  notifications, local mode): SMTP host/port/SSL, sender, and an app password are
  saved to a local `.email.json` (0600, password never echoed back to the
  browser). `send_mail.py` reads that config (falling back to the
  `CRISPRME_SMTP_*` env vars) and is now best-effort — it never writes to stderr
  or raises, so a mail-server problem can't fail an already-completed search.

### Changed
- **Search form polish**: the CRISPRme+ logo now sits above the first step; the
  *Load Example* button moved to the bottom of the *Select gRNA* box; larger,
  more legible fonts throughout the form; the spacer help text explains that
  multiple spacers can be entered (and points to the CLI for bulk runs); the
  *Notify me by email* box explains the one-time Settings setup.

### Fixed
- **Variant indel indexing/search now resolves the INDELS companion from the
  actually-resolved variant index** instead of reconstructing an exact
  `<true_pam>_<bMax>_..._INDELS` name. A batteries/precomputed index whose N or
  PAM differs from the request (e.g. the shipped `NGG_3` served against a
  `bMax=2` search, or a pamless `NNN` index) was not matched, so the pipeline
  attempted an impossible rebuild and crashed in `pool_index_indels.py`
  (`FileNotFoundError: .../Genomes/<ref>+<vcf>_INDELS/`). The INDELS index, its
  search path, and the on-demand builder now all key off the index's own
  `<PAM>_<N>` prefix (consistent with `build-index-only`/`publish`/`download`),
  so newly built-and-uploaded indexes work generally. When the companion is
  genuinely absent and cannot be rebuilt (no source genome), the run now fails
  with a clear "re-download the index" message instead of a raw traceback.

### Added
- **Annotation manager** in Settings → Data Manager: a per-genome
  *Manage annotations (enable / disable)* checklist backed by a persisted
  `.enabled.json` manifest, plus format validation on every uploaded annotation
  BED (≥4 tab columns, integer `start ≤ end`, whitespace-free label, size cap;
  malformed files are rejected with a clear message and nothing is written).
  Enabled tracks are merged on demand into a single cached active annotation and
  applied automatically to every search.

### Changed
- Annotations are now applied automatically from the enabled set rather than
  chosen per search. The annotation dropdown was removed from the search form
  (Step 3 is now just email + job name); the built-in
  ENCODE cCREs (SCREEN) + DHS + GENCODE bundle is enabled by default so
  annotations work out of the box.
- The **Maximum edits** slider is floored at 1 (the on-target / 0-edit hit is
  always reported at any setting), avoiding an empty-result edge case on large
  precomputed indexes.
- Search-submission form: consistent left indentation for the PAM and
  variant-dataset dropdowns and the submit/example buttons.

## [2.2.0] - 2026-08-06

### Added
- Graphical **Settings / Data Manager** page in the web interface: add reference
  genomes (UCSC by assembly name — e.g. the pig `susScr11` — HuggingFace, or a
  direct URL), precomputed indexes (download from HuggingFace or build locally
  from an installed genome + PAM), VCF datasets (HuggingFace or register an
  existing server folder), annotations (BED upload), and nucleases/PAMs (a small
  form). New data lands in the local data folder and is auto-discovered by the
  search form. Long operations run as detached jobs on a dedicated executor with
  live progress, so they never starve the search slots. Mutations are local-mode
  only; publishing an index to the shared HuggingFace repo is maintainer-only.
  `download --what genome` gained `--source {hf,ucsc,url}` (+ `--url`) so the CLI
  and web share one non-human-genome download path.
- Python 3.11 modernization: pipeline fixes for pandas 2.x / matplotlib 3.x and
  a Dash 1.x → 2.x web-app migration, plus a Python-3.11 Docker image built from
  source (CRISPRitz 2.8.1) (#131).
- `assembly-search` subcommand: off-target search on a personal diploid genome
  assembly (two haplotypes, no VCF), reconciled to hg38 via liftOver (#113).
- Reference-index UX: `build-index-only` pre-builds the reusable CRISPRitz
  reference index without running a search, and `complete-search --index-path`
  reuses a prebuilt/staged index library (a missing index is a hard error rather
  than a silent rebuild).
- HuggingFace data distribution: `download` fetches reference data (genome,
  annotations, PAMs, sample IDs, VCFs, precomputed indexes) from a HuggingFace
  dataset repository over its CDN, and `publish-index` uploads a locally built
  index for reuse. Default repo `lucapinello/crisprme-data`, overridable via
  `--hf-repo` / `CRISPRME_HF_REPO`. `setup`/`complete-test` also try HuggingFace
  first and fall back transparently to the original UCSC/EBI/Sanger sources
  (#140, #141).
- `complete-search --max-total-edits N`: cap the total edits (mismatches +
  bulges) per reported alignment; over-cap targets are dropped right after the
  search, shrinking intermediate files and post-analysis time (#107).

### CI
- New `unit tests` workflow: fast, hermetic byte-compile + network-free HF/index
  unit tests on every code PR.
- New `web e2e (playwright)` workflow: builds the py3.11 image, serves the web
  app, and drives Chromium to assert every Dash 2.x page renders (no blank pages
  / JS errors).
- `validate-benchmarks` gained a new-subcommand dispatch + unit-test smoke step.

### Changed
- Clearer failure reporting: when a search fails, CRISPRme now prints *which
  stage* failed (from the per-stage log) and the last lines of the error log,
  instead of only "run failed — see log_error.txt". Makes failures actionable
  for non-expert users.

### Fixed
- Web interface (Dash 2.x) hardening, from a full Playwright stress test of the
  running app:
  - The web server no longer crashes on a from-source install. Dash 2.x's
    `app.run()` lets the `HOST` environment variable override the host argument,
    and the from-source conda env sets `HOST` to a non-bindable compiler build
    triple (`x86_64-conda-linux-gnu`); the server now forces the intended host/
    port so it binds correctly.
  - The **Query Genomic Region** and **Personal Risk Cards** result tabs no
    longer return HTTP 500. Both callbacks type-checked their inputs before the
    "no click yet" guard, so Dash's initial (empty) render raised a `TypeError`;
    the guard now runs first, and Filter/Generate with nothing selected is a
    graceful no-op.
  - Removed the dead cross-origin "skeleton" stylesheet (blocked by browsers on
    every page); the layout already uses the Bootstrap grid.
  - The nuclease dropdown collapses case-variant duplicate PAM files so each
    nuclease is listed once.
- Zero-hit searches now complete cleanly with an empty result instead of
  aborting. A search that finds no off-targets (e.g. a very stringent
  guide/parameter combination) previously failed part-way through post-analysis
  ("off-targets post-analysis (reference) failed", then a cascade through the
  rsID / summary / integration steps, all of which assumed at least one target).
  Added a zero-target guard to the reference SNP post-analysis (mirroring the
  INDELs one), made `remove_n_and_dots.py` tolerate a header-only report, and
  added a high-level short-circuit that emits an empty-but-valid result and exits
  0 when no off-targets are found. Verified end-to-end on ml007 (empty result
  exits 0; a normal with-hits search is unaffected).

### Documentation
- New `docs/DOCKER_QUICKSTART.md` + a README quickstart callout: a few-command,
  no-conda / no-410 GB path to the web interface for non-experts — fast HF data
  download, a prebuilt index, then `web-interface` in the browser, with
  troubleshooting and "install more indexes as needed".
- Data-setup guide: documented what `setup` produces (including the combined
  1000G+HGDP config files), the HuggingFace fast-download path, and a new
  "Prebuild, reuse, and share the reference index" section.
- README: added a from-source install path and a reference-index /
  data-distribution commands section.
- Documented the variant PAM behaviour in `docs/INPUT_FORMATS.md`: PAM *creation*
  (variants that add an off-target) is reported; PAM *disruption* (variants that
  remove a reference PAM for carriers) is a known, intentional gap — with the
  rationale (disruption lowers rather than raises predicted risk and would need
  sample-specific, non-scorable semantics).

## [2.1.13] - 2026-08-04

### Added
- Robust `complete-test` downloads: HTTP retry with linear backoff plus
  checksum-verified resume, so a dropped connection or truncated transfer from a
  slow/flaky host retries (and skips already-verified files) instead of aborting
  a multi-hour run (#136).
- Memory-bounded post-analysis: per-chromosome SNP/indel workers are capped to a
  RAM budget (default 64 GB, overridable via `CRISPRME_MAX_MEM_GB` and
  `CRISPRME_POSTPROC_WORKER_GB`), preventing the peak-RAM spike observed on
  genome-wide 1000G runs (#136).

### Fixed
- Post-analysis worker-count diagnostic is now written to stdout (`log_verbose`)
  instead of stderr. The caller treats a non-empty stderr log as a fatal
  post-analysis failure, so the informational message previously aborted every
  run with a false "post-analysis (snps) failed" even when the analysis
  succeeded (#136).

## [2.1.12] - 2026-08-01

### Added
- Cas9 and Cas12a validation benchmark examples with precomputed brute-force
  ground-truth references and an extensible registry
  (`test/benchmark/benchmarks.json`); `validate-test` compares CRISPRme
  off-targets against them (#116).
- A fast, dependency-free Rust port of the brute-force ground-truth generator
  (`test/benchmark/rust/`) that produces output identical to the Python
  generator (#116).
- `validate-benchmarks` GitHub Actions CI that runs the full
  `complete-test` → `validate-test` round-trip on Linux for every registered
  benchmark (#116, #120).
- Native Apple Silicon (`linux/arm64`) Docker image via a multi-arch `buildx`
  workflow; the multi-arch image is published to Docker Hub on release tags
  (#121).
- Low-memory startup warning (`PostProcess/memory_check.py`) advising Docker
  Desktop users to raise the memory limit when total RAM is below ~32 GB (#121).
- Search-space budget estimator and guard (`PostProcess/search_budget.py`) that
  warns before resource-explosive runs, with a companion analysis in
  `docs/SCALABILITY_ANALYSIS.md` (#118).
- Release tooling: a `release-crisprme` Claude Code skill, `docs/RELEASING.md`,
  this `CHANGELOG.md`, and `scripts/prepare_release.py` (#117); plus docs for the
  fast Rust generator and using the release skill (#122).
- Input-format hardening: `--gene_annotation` is now auto-sorted like
  `--annotation`; clear errors for multi-line PAM files and malformed PAM
  filenames; and a non-fatal warning when a degenerate PAM motif is combined with
  bulges (the CRISPRitz #105 crash, fixed in CRISPRitz 2.7.1) (#125).
- Documentation: `docs/INPUT_FORMATS.md` (PAM/Cas12a/VCF/chromosome/annotation/
  bulge guidance, #124), `docs/ISSUE_AUDIT.md` (a review of all open + closed
  issues, #126), and `docs/ROADMAP.md` (the 2.1.12 → 2.1.13 plan, #128).

### Fixed
- Multiple bugs affecting new VCF dataset processing, concurrent runs, and
  gnomAD handling (#96).
- Guard the CFD PAM-score lookup against non-ATCG PAM bases present in real hg38,
  which previously crashed the run with a `KeyError` (#94, #125).
- Guard the post-analysis intermediate reads against non-UTF-8 bytes
  (`UnicodeDecodeError`) so a stray byte no longer aborts a run (#52, #125).
- `validate-test` now exits non-zero when a benchmark mismatches or the
  `complete-test` output is missing, instead of passing silently (#116).
- `complete-test` downloads 1000 Genomes VCFs over HTTPS instead of FTP,
  unblocking FTP-restricted and CI networks (#116).

### Changed
- Synced the `Dockerfile` `crisprme_version` pin to the released version, which
  had lagged behind `crisprme.py` and the Bioconda recipe (#117).

## [2.1.11] - 2026-07-06

### Fixed
- Consistent `.gz` suffix for compressed annotation files, removing a mismatch
  between file content and file name in the standard and personal-annotation
  pipelines (`crisprme.py`).
- `_check_personal_annotation` now strips a `.gz` suffix before processing so
  downstream steps see the expected filename.
- Removed double-compression of annotation files, preventing corrupted or
  unreadable annotation outputs.

### Changed
- Refactored per-chromosome test genome-directory handling: replaced
  `ensure_hg38_directory` with `_assign_genome_directory_name` (whole-genome
  builds in `<dest>/hg38`, single-chromosome builds in `<dest>/hg38_<chrom>`).
- `download_genome_data` now returns the resolved genome directory path directly.
- Renamed downloaded sample-ID files to include genome/dataset context
  (`hg38_{dataset}.samplesID.txt`) via a centralized mapping.
- Reduced default web-interface example-run parameters (max mismatches 6 → 4,
  DNA bulges 2 → 1, RNA bulges 2 → 1) for a faster example run.

### Added
- New `_write_sg1617_file` helper to generate an `sg1617` guide fixture during
  legacy database setup.
- New guides `docs/crisprme_data_setup_051826.md` (CLI setup/usage) and
  `docs/crisprme_web_interface_user_guide.md` (web interface walkthrough).
- "Setup Legacy Database" section in `README.md`.

### Removed
- Unused `_bgzip_ann_data` helper and stale typing imports/docstrings.

## [2.1.10] - 2026-05-29

### Added
- New `setup` command for automated download/installation of CRISPRme reference
  resources directly from the command line, removing the need for manual
  downloads from the CRISPRme website.

### Changed
- Updated documentation and container configuration to support the new
  installation workflow and improve reproducibility across deployments.

## [2.1.9] - 2026-01-16

### Added
- `validate-test` functionality: validates `complete-test` off-target
  predictions against brute-force ground-truth alignments derived from
  1000 Genomes variant data (PR #92).
- Chromosome-level validation support (single chromosome or genome-wide).

### Changed
- Strengthened the `complete-test` pipeline and refined benchmarking utilities
  for correctness and reproducibility with population-scale variant data.

## [2.1.8] - 2025-12-10

### Fixed
- Targeted fixes and stability improvements across annotation handling,
  off-target search execution, test coverage, and report generation.
- Fixed and updated the `complete-test` workflow to align with the latest
  internal logic and outputs (PR #86).

### Changed
- Updated Docker distribution and testing environment for consistent runtime
  behavior (PR #85).

## [2.1.7] - 2025-03-10

### Added
- Support for exome-based gnomAD VCF files.

### Fixed
- Corrected handling of indels in off-target reporting (previously missing or
  misreported).
- Fixed the visualization of the Personal Card tab in the GUI and website.
- Closed file streams after creating graphical reports to prevent memory issues.

### Changed
- Improved contiguous-target merge logic (always reports the leftmost off-target
  in a cluster).
- Improved error-tracing logic.

## [2.1.6] - 2024-11-27

### Added
- Support for gnomAD v4.1 VCFs, including joint variant files.
- New `complete-test` function for quick single-chromosome and full-genome
  off-target testing.

### Fixed
- Corrected indel handling in off-target reporting to eliminate false positives.
- Fixed the alternative-alignments report to list all possible alignments.
- Resolved Personal Card tab visualization issues in the GUI and website.
- Addressed a memory-overflow issue when merging contiguous targets.

### Changed
- Upgraded the DockerHub image with the latest fixes.

[Unreleased]: https://github.com/pinellolab/crisprme-plus/compare/v2.8.1...HEAD
[2.8.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.8.1
[2.8.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.8.0
[2.7.10]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.10
[2.7.9]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.9
[2.7.8]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.8
[2.7.7]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.7
[2.7.6]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.6
[2.7.5]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.5
[2.7.4]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.4
[2.7.3]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.3
[2.7.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.2
[2.7.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.1
[2.7.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.0
[2.6.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.6.2
[2.6.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.6.1
[2.6.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.6.0
[2.5.5]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.5
[2.5.4]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.4
[2.5.3]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.3
[2.5.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.2
[2.5.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.0
[2.4.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.4.0
[2.3.3]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.3
[2.3.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.2
[2.3.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.1
[2.3.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.0
[2.2.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.2.0
[2.2.0-alpha.30]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.2.0-alpha.30
[2.1.13]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.13
[2.1.12]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.12
[2.1.11]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.11
[2.1.10]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.10
[2.1.9]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.9
[2.1.8]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.8
[2.1.7]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.7
[2.1.6]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.6
