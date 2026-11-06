<!-- defaults track (n290): CRISPRme 2.9.0 entries for CHANGELOG.md. NRR as the default PAM and
index, the index-flexibility items of plan 2.7 (E1), index names that carry the guide length and
the side (E2), and the clear fixes of the backlog. -->

**For the merge, not for the changelog.** The install sizes below (150 GB on disk, a 70 GB
download) were measured on ml003 on 2026-11-06 with the two documented downloads; the peak while
the archives unpack is given where it was measured. Every published NRR archive is the 2.8.1 one:
the rebuild lane's archives change these sizes again (by 1-2 %, plan 2.5), so re-measure with them.

### Changed

- **SpCas9 NRR (`NAA`, `NAG`, `NGA`, `NGG`) is the default PAM and index; it was NRG (`NAG`,
  `NGG`).** The broadest SpCas9 screen of the three published ones, and the one recommended for an
  IND package: an NRR search also reports the `NGA` and `NAA` sites, and the variant-created
  ones, that an NRG search leaves out, as NRG reports the `NAG` ones NGG leaves out. The PAM file
  `PAMs/20bp-NRR-SpCas9.txt` ships with CRISPRme, and a download of the PAM files (`download --what
  pams`, and so `--what all`) installs it where the data repo has none (the repo's `pams/` holds the
  NRG, NGG and PAM-less files). The web form selects it when it is installed (then NRG, then NGG, as
  before), `download --what all` fetches the reference index `NRR_3_hg38`, and every surface that
  names the recommended index names `NRR_3_hg38+hg38_HGDP1kGP` (and its HPRC and mega panels):
  the README's and the website's commands, the Docker Quickstart, the web user guide, the CLI guide,
  PRECOMPUTED_INDEXES, METHODS, the Hugging Face cards, RELEASE_CHECKLIST's default-index gate, the
  help of `download`, `setup` and `complete-test`, complete-test's advice when no index is
  installed, and the advice of a search and of the report to download a tiered index.
  complete-test's smoke searches a PAM-less 8-nt-flank index for NRR (it searched it for NRG). The
  NRG and NGG indexes stay published and searchable. What it costs: the default install is a ~70 GB
  download (it was ~45 GB) and takes ~150 GB on disk (it took ~98 GB; measured 150,048,109,201
  bytes: indexes 94.2 GB, `Dictionaries/` 52.4 GB, genome 3.3 GB), and a genome-wide search takes
  about twice as long (2 h 26 min against 1 h 14 min for `ACTGAAATCTGTAAGCAGGC` on HGDP + 1000
  Genomes, per-sample, 4 mismatches + 1 + 1 bulge, 32 threads, CRISPRme 2.8.1), at no higher peak
  memory (66 GB, against 82 GB summed over the NRG run's processes). Covered by
  `PostProcess/test_defaults_and_index_flex.py`.
- **`--prune-pam` no longer needs `--index-path`**, and a search passes it by itself (see Added).
  The refusals that told a user to run again with `--prune-pam --index-path` now say, where the
  engine cannot prune, that an index serves the search pruned and the installed CRISPRitz cannot do
  it; the advice to run without `--index-path` no longer says to drop `--prune-pam` too.

### Added

- **A dedicated index serves the PAMs its own holds.** An `NRR` index stores every site of `NGG`,
  `NAG`, `NGA`, `NAA` and `NRG`, an `NRG` index every site of `NGG` and `NAG`, as concrete bases in
  their stored flank; under `--prune-pam` the engine writes out only those of the search's PAM, at
  leaf time, as it does an all-`N` index's. A search takes such an index after its PAM's own and
  before a PAM-less one (of two, the narrower PAM: `NRG` before `NRR` for `NGG`), and only for a
  PAM on the side of the guide it was built for. It was refused with `--index-path`, and a search
  without it built an index of its PAM (20 to 35 minutes for hg38; a variant search stopped). One
  rule, `index_fit.pam_serves` and the search script's `_pam_holds`, in the search script,
  `crisprme.py`, `unfinished_marks` and the web form. Through the real engine (the release gate's
  `test/index_flex/check_index_flex.py`): `NGG` from an `NRG` or an `NRR` index, and `NRG` from an
  `NRR` index, give the same results as a dedicated index of the PAM.
- **A search prunes by itself.** When no installed index serves a search unpruned and one serves it
  pruned (a dedicated index of a PAM that holds the search's, or an all-`N` one wider than the PAM),
  the search prunes, as if `--prune-pam` was given, and says so first: on the console (`NOTE
  [complete-search]: no installed index serves PAM NGG unpruned, and NRR_3_hg38+hg38_HGDP1kGP (...)
  serves it pruned to it at leaf time: searching with --prune-pam, as if it was given. A dedicated
  NGG index is searched faster: ...`) and in its log. The web form, which never passes the flag,
  gets the same: the search script decides it for the search it runs, and the form's bulge cap and
  PAM list count those indexes. Before, a search that only a K = 8 index served stopped for want
  of a variant index (or, for the reference index, built one), and said to pass the flag. Where an
  index serves the search unpruned, as an all-`N` one of the PAM's width does, nothing changes.
- **The run's parameters and the report name the index that served the search.** `.Params.txt`'s
  `Genome_idx` says what the search asked for, which another index may serve (`NGG_2_hg38` served
  by `NNNNNNNN_3_hg38`), or none at all (IDX-21). The search script now records the indexes it
  took: `Reference_index` and `Variant_index:<dataset>`, each with its library, the PAM it was
  built for and its side of the guide, the guide length its trees were built for, its bulges and
  how the search narrows it to the search's PAM, the variant one with its INDELS index, or `none`
  for a search with no bulges, which scans the genome; and `Prune_pam`, `auto` or `given`. The
  report's summary shows them, a run of an earlier release has none.
- **Index names may carry the guide length and the side of the PAM.** `<PAM>_<N>_<genome>` says
  neither, so a 3′ and a 5′ build of one motif, or builds for two guide lengths (20- and 23-nt
  Cas12a), take one name. A name may now carry both, in the order of the target on the strand:
  `20-NRR_3_hg38` for a 20-nt guide and a 3′ PAM, `TTTV-23_3_hg38` for a 5′ PAM and a 23-nt guide.
  Every reader takes both forms (the search script, `crisprme.py`, the web form, `download` and its
  router, which reads a PAM-less name's flank in either form, `publish-index`'s manifest,
  complete-test, the index coverage and N-flank checks): an index of the new form only for a search
  of the side its name says, and of one index under both names, the newer one. The INDELS index of
  a variant index so named is `<name>_INDELS`; one a search builds for it is moved there. Builds
  still write the older name. Covered by `PostProcess/test_defaults_and_index_flex.py`, and through
  the real engine by `check_index_flex.py`: a download renamed to the new form gives the same
  results, and one named for a 5′ PAM is not taken for a 3′ search.
- **The web form's PAM picker for a PAM-less index** (plan 2.7, DYNAMIC_PAM_INDEX_COMPAT §8). When an
  all-`N` index is installed for the genome (or the dataset) selected, the form says what it serves
  (any PAM of its guide length on its side, up to its flank) and offers the PAMs of common nucleases
  it serves (SpCas9 NRR, recommended for IND, then NRG and NGG, SpG, SpRY, VQR, VRER, xCas9, ...),
  those it does not serve disabled with the reason, and a field for a PAM of one's own, checked when
  the field is left. "Use this PAM" writes that PAM's file into `PAMs/` and selects it; the search
  prunes the index to it by itself. Only a PAM the index serves, and a nuclease of the table, is
  written. Covered by `PostProcess/test_defaults_and_index_flex.py` and
  `PostProcess/test_web_bulge_notes.py` (in the image).

### Fixed

- **A search with no bulges no longer builds an index.** It scans the genome's FASTA files, and a
  variant search the enriched genome's when it is here whole, so it opens no index; it built the
  reference index `<PAM>_1_<genome>` all the same (hours on hg38), and a variant one, which nothing
  read, and under `--index-path` it stopped on a library that held no reference index. Its INDELS
  index, which the indel search reads at any bulges, is still found or built as before. The web's
  progress page still reads the indexing step as done. A personal-assembly search with no bulges
  no longer builds a reference index per haplotype (R4-4).
- **Load Example selects the PAM and the genome it computed.** It picked the installed default PAM
  and then returned "20bp-NRG-SpCas9" and "hg38" whatever was installed (FR :125).
- **The web form's bulge cap read an index built for the other side of the guide as serving.** Its
  mirror of the search script's choice had no side rule (F2): it now takes an index only for a PAM on
  the side its `.pam_build` (or its name) records, as the search does.
- **`crisprme.py setup` exits non-zero when its download fails.** It caught the legacy script's
  failure, wrote `Legacy database setup exited with error code 1` with no end of line, and exited 0
  (with no network, a `requests` traceback and status 0).
- **setup's test-data URLs no longer hold `//`.** The base URL ended in `/` and each URL added
  another (`.../test/data//samplesIDs/...`); GitHub answered 307 and the download worked.
- **`generate-report` says why a folder has nothing to report.** On a search that found no
  off-target, and on an assembly-search's combined folder, it ended on a `FileNotFoundError`
  traceback (OF-12). It now says, in one line, that the folder holds no integrated results (and
  why), or that assembly-search writes the combined report there itself.
- **A PAM-less search records its nuclease as `NO-PAM` from the command line too.** complete-search
  took the third field of the PAM file's name, `NO`, where the web reads `NO-PAM` (OF-13).
- **environment.yml installs pigz**, which the Docker image has and `publish-index` packs an index
  with on every core (install P9).
