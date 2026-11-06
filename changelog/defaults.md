<!-- defaults track (n290): CRISPRme 2.9.0 entries for CHANGELOG.md. NRR as the default PAM and
index, the index-flexibility items of plan 2.7 (E1), index names that carry the guide length and
the side (E2), and the clear fixes of the backlog. -->

**For the merge, not for the changelog.** The install sizes below (150 GB on disk, a 161 GB
peak, a 70 GB download) were measured on ml003 on 2026-11-06 with the two documented downloads,
twice, the peak by sampling `du -sb` every 10 s. Every published NRR archive is the 2.8.1 one: the
rebuild lane's archives change these sizes again (by 1-2 %, plan 2.5), so re-measure with them,
and change `PostProcess/test_documented_sizes.py` and every surface it reads together.

**The backlog: what this track fixed, and what it leaves, and why** (list_282.md and plan
section 4, for the merge). Fixed here, each with its test (the entries below): IDX-17a, IDX-19,
IDX-21, IDX-22, cli5 D7, OF-12, OF-13, FR :125, F2, RV1-F5, RV4-1, RV4-3, M2, M3, CI-1, WEB-9 (and
the non-archive bundle answered 500), C11/P7, emptyidx F6, F7 (and the genome-layout refusal's
banner, WEB-3's class), F8, suite P4, skips D7, P15.3, install P9, R4-4's unused NRG_1 index of
a 0-bulge assembly search, setup's exit status, progress lines and URLs (FR2 cli a, c, e), the
checklist's HPC route (`--cleanenv`, a fresh sandbox), backlog 4.4 (the installers). Found done
already: P12 (complete-test with no index exits 2 since 2.8.1). Left:
- Design questions, or decisions taken for another release: all of section 4.1 (two cis indels,
  CRISPR-Bulge batch noise, the one-row drops, phasing truth, phase sets, #77, the calibration
  edge, the HPRC and HG01255 runs); the report's A6 rsID/MAF order (with the multi-allelic work,
  plan 2.2), P7 (what the slim report shows), the 0-1-edit sites of the shareable report, the
  CRISPR-Bulge of reconciled assembly rows; WEB-7 (plan 2.6), WEB-4, WEB-11 (a default), the
  annotation copies of a web search (the #97 trade), Settings builds without the tiers' RAM; the
  stderr gate of the pipeline stages (exit codes instead, 44 gates); build locks across hosts.
- Other tracks' code: the served-mode security items (the security track's round 3, A1); the
  report's RV1-F1/F2, D6, the 1e-5 wording, P8, mega IUPAC rows, Risk Card sort, text N1 (on its
  branch); the engine's AF-by-first-record and 5' last-window items, the Bioconda enricher, PKG1
  (format 2: the engine and fmt tracks); download's completeness check of the dictionaries (its
  heal path refetches the whole archive: the HF track's); contigs with "_" (bio11.3: the indel
  pipeline's naming, with its own tests); R2-F1, M5, R4-6 (the build tracks' locks and marks).
- Not reproduced, or needing a browser walk this track had no time for: P14 (the results page's
  Group by now opens on CFD_score), N1-N8, M6, OF-11, the base-editing dropdown's 500s, history
  and Settings hashed haplotype names, the Settings build form's genome list, W1.
- Test hygiene left: suite P5 (tests that import the web app write `PostProcess/Cache/`, the app's
  cache folder, relative to the working directory by design: moving it is the web track's call;
  the folder is ignored by git), bio11.2 (a flake recorded with no detail), the T-31 flock trigger
  (could not be reproduced on the Mac), R2-F4/R2-F5 and wset N3 coverage, the 5'-PAM synthetic case
  and a published test archive (audit additions, sized small-medium).
- Docs and wording nits of 4.6 and 4.8 (RV1-F3, N4, F4, F5, RV4-2, the 2.7.10 entry, the FFR
  D3R3 items): the docs track's, as they touch other tracks' changelog fragments.

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
  the README's and the website's commands, the one-click installers (`docs/install.sh`,
  `docs/install.ps1`), the Docker Quickstart, the web user guide, the CLI guide,
  PRECOMPUTED_INDEXES, METHODS, the Hugging Face cards, RELEASE_CHECKLIST's default-index gate, the
  help of `download`, `setup` and `complete-test`, complete-test's advice when no index is
  installed, the advice of a search and of the report to download a tiered index, and the notes
  of the panel build recipes (`seq_script/merge_panels/`) and of the paper's dataset table.
  complete-test's smoke searches a PAM-less 8-nt-flank index for NRR (it searched it for NRG). The
  NRG and NGG indexes stay published and searchable. What it costs: the default install is a ~70 GB
  download (it was ~45 GB) and takes ~150 GB on disk (it took ~98 GB; measured 150,048,109,201
  bytes: indexes 94.2 GB, `Dictionaries/` 52.4 GB, genome 3.3 GB), with a peak of 161 GB while the
  archives unpack, so every surface that sizes it (the installers, the README, the website, the
  Quickstart, TROUBLESHOOTING, DATA_LAYOUT, install/README) asks for ~200 GB free (it asked for
  ~150 GB); a reference-only setup takes ~44 GB (it took ~25 GB). A genome-wide search takes
  about twice as long (2 h 26 min against 1 h 14 min for `ACTGAAATCTGTAAGCAGGC` on HGDP + 1000
  Genomes, per-sample, 4 mismatches + 1 + 1 bulge, 32 threads, CRISPRme 2.8.1), at no higher peak
  memory (66 GB, against 82 GB summed over the NRG run's processes). Covered by
  `PostProcess/test_defaults_and_index_flex.py`.
- **`--prune-pam` no longer needs `--index-path`**, and a search passes it by itself (see Added).
  The refusals that told a user to run again with `--prune-pam --index-path` now say, where the
  engine cannot prune, that an index serves the search pruned and the installed CRISPRitz cannot do
  it; the advice to run without `--index-path` no longer says to drop `--prune-pam` too.
- `PostProcess/pool_search_indels.py` drops `_index_pam`, a helper nothing called (emptyidx F8);
  the results page's sample tab no longer computes paths it never reads (wsearch N8); and
  `remove_contiguous_samples.py` reads a variant token as an indel's in one place, for the
  representative's tie-break and the carrier resolution alike (RV3-D).

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
  `NRR` index, give the same results as a dedicated index of the PAM. On chr22 with HGDP + 1000
  Genomes (two guides, 4 mismatches + 1 + 1 bulge, per-sample, 16 threads), an `NGG` search that
  the `NRR` index alone served, pruned by itself, wrote the same 2,023 rows, byte for byte, as one
  of the dedicated `NGG` index, in 397 s against 295 s, at 9.6 GB peak memory against 2.6 GB
  (CRISPRitz 2.9.1); with CRISPRitz 2.10.0 both wrote those rows again.
- **A search prunes by itself.** When no installed index serves a search unpruned and one serves it
  pruned (a dedicated index of a PAM that holds the search's, or an all-`N` one wider than the PAM),
  the search prunes, as if `--prune-pam` was given, and says so first, naming the index that needs
  it: on the console (`NOTE [complete-search]: no installed variant index of hg38+hg38_HGDP1kGP
  serves PAM NGG unpruned, and NRR_3_hg38+hg38_HGDP1kGP (...) serves it pruned to it at leaf time:
  searching with --prune-pam, as if it was given. A dedicated NGG variant index is searched faster:
  ...`; in a library of the PAM's own variant index and another PAM's reference index, `no
  installed reference index of hg38`) and in its log. The web form, which never passes the flag,
  gets the same: the search script decides it for the search it runs, and the form's bulge cap and
  PAM list count those indexes. Before, a search that only a K = 8 index served stopped for want
  of a variant index (or, for the reference index, built one), and said to pass the flag. Where an
  index serves the search unpruned, as an all-`N` one of the PAM's width does, nothing changes. On
  chr22 with HGDP + 1000 Genomes (CRISPRitz 2.10.0), an `NGG` search whose variant index was
  `NGG`'s own and whose reference index only the `NRR` one, under its older name or as
  `20-NRR_3_hg38` with no `.pam_build`, said `no installed reference index of hg38`, and wrote
  the same 2,023 rows, byte for byte, as the dedicated `NGG` indexes.
- **The run and its report name the index that served the search.** `.Params.txt`'s `Genome_idx`
  says what the search asked for, which another index may serve (`NGG_2_hg38` served by
  `NNNNNNNN_3_hg38`), or none at all (IDX-21). The search script now records the indexes it took,
  in `.indexes_searched.txt` in the result folder, beside the parameters (which the web compares
  byte for byte to reuse an earlier job, so they stay as the request wrote them): `Reference_index`
  and `Variant_index:<dataset>`, each with its library, the PAM it was built for and its side of
  the guide, the guide length its trees were built for, its bulges and how the search narrows it
  to the search's PAM, the variant one with its INDELS index, or `none` for a search with no
  bulges, which scans the genome; and `Prune_pam auto` for a search that pruned by itself
  (complete-search records `Prune_pam`, `auto` or `given`, in its parameters). The report's
  summary shows them, and its zip bundles the record as `indexes_searched.txt`; a run of an
  earlier release has none.
- **Index names may carry the guide length and the side of the PAM.** `<PAM>_<N>_<genome>` says
  neither, so a 3′ and a 5′ build of one motif, or builds for two guide lengths (20- and 23-nt
  Cas12a), take one name. A name may now carry both, in the order of the target on the strand:
  `20-NRR_3_hg38` for a 20-nt guide and a 3′ PAM, `TTTV-23_3_hg38` for a 5′ PAM and a 23-nt guide.
  Every reader takes both forms (the search script, `crisprme.py`, the web form, `download` and its
  router, which reads a PAM-less name's flank in either form, `publish-index`'s manifest,
  complete-test, the index coverage and N-flank checks): an index of the new form only for a search
  of the side its name says, which records the side as a `.pam_build` does, ahead of one (such an
  index needs none, at any flank width), and of one index under both names, the newer one. The
  INDELS index of
  a variant index so named is `<name>_INDELS`; one a search builds for it is moved there. Builds
  still write the older name. Covered by `PostProcess/test_defaults_and_index_flex.py` (the
  search script, `crisprme.py`, the web and complete-test take the same of every name form, side
  and record), and through the real engine by `check_index_flex.py`: a download renamed to the new
  form gives the same results, one named for a 5′ PAM is not taken for a 3′ search, and a
  PAM-less index named for its side with no `.pam_build` serves a search of that side (a 3′ one of
  width 4 an NGG search, pruned; a 5′ one a 5′ PAM of its width), which the search script refused
  as an unrecorded index that crisprme.py and the web offered.
- **The web form's PAM picker for a PAM-less index** (plan 2.7, DYNAMIC_PAM_INDEX_COMPAT §8). When an
  all-`N` index is installed for the genome (or the dataset) selected, the form says what it serves
  (any PAM of its guide length on its side, up to its flank) and offers the PAMs of common nucleases
  it serves (SpCas9 NRR, recommended for IND, then NRG and NGG, SpG, SpRY, VQR, VRER, xCas9, ...),
  those it does not serve disabled with the reason, and a field for a PAM of one's own, checked when
  the field is left. "Use this PAM" writes that PAM's file into `PAMs/` and selects it; the search
  prunes the index to it by itself. Only a PAM the index serves, and a nuclease of the table, is
  written. A server online writes nothing into the data folder, as its Settings page does not:
  there the picker selects the PAM file installed for the PAM, whatever its name, offers only the
  nucleases whose PAM file is installed, and says that it adds none for any other. Covered by
  `PostProcess/test_defaults_and_index_flex.py` and `PostProcess/test_web_bulge_notes.py` (in the
  image).

### Fixed

- **A search with no bulges no longer builds an index.** It scans the genome's FASTA files, and a
  variant search the enriched genome's when it is here whole, so it opens no index; it built the
  reference index `<PAM>_1_<genome>` all the same (hours on hg38), and a variant one, which nothing
  read, and under `--index-path` it stopped on a library that held no reference index. Its INDELS
  index, which the indel search reads at any bulges, is still found or built as before. The web's
  progress page still reads the indexing step as done. A personal-assembly search with no bulges
  no longer builds a reference index per haplotype (R4-4). On chr22 with HGDP + 1000 Genomes, a
  search with no bulges of `NGCG`, which no installed index holds, ran in 67 s and built only
  `NGCG_1_hg38+hg38_HGDP1kGP_INDELS`; its record of the indexes searched says `Reference_index
  none` and `Variant_index:hg38_HGDP1kGP none`, each with why.
- **A search its input check refuses leaves its output folder empty.** A guide without its
  PAM's `N` padding, for one, is refused before the search starts, but the run had written
  `.command_line.txt`, `.version.txt`, `Params.txt` and `guides.txt` into `Results/<output>/`,
  so the corrected command, with the same `--output`, stopped on `Output folder ... not
  empty!`. Found running this release's chr22 searches.
- **Load Example selects the PAM and the genome it computed.** It picked the installed default PAM
  and then returned "20bp-NRG-SpCas9" and "hg38" whatever was installed (FR :125).
- **The web form's bulge cap read an index built for the other side of the guide as serving.** Its
  mirror of the search script's choice had no side rule (F2): it now takes an index only for a PAM on
  the side its `.pam_build` (or its name) records, as the search does.
- **`crisprme.py setup` exits non-zero when its download fails.** It caught the legacy script's
  failure, wrote `Legacy database setup exited with error code 1` with no end of line, and exited 0
  (with no network, a `requests` traceback and status 0).
- **setup's progress lines come as they are printed.** Into a pipe or a log file, the download
  script's progress lines came at its end, after its errors, so a merged log read out of order.
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
- **The one-click apps start nothing on half the data after a failed first download.** The
  Windows app's first Start ran `download --what all`, the variant index's download and the app
  with `;` between them, so a failed download went on to the next and started the app on half the
  data; each step now runs only when the one before it succeeded (`$LASTEXITCODE`: Windows
  PowerShell 5.1 has no `&&`) (backlog 4.4). And both apps, the macOS one too, chose between
  setting up and opening the app by `Genomes/` alone, which the download makes in its first
  minutes, so the next Start after a failed download opened the app on part of the data. The
  setup now marks the data folder (`.setup_unfinished`) until both downloads have ended well: a
  Start after a failure sets up again, picking up what is there, and the failure says to click
  Start again to finish setting up. Checked with the macOS app's AppleScript under `osascript`
  and `install.ps1` under PowerShell 7, each with a stub `docker` that fails the index download,
  then not; an install of an earlier installer, with `Genomes/` and no mark, opens as before.
- **A search no longer edits the user's `--vcf` and `--samplesID` lists.** The search script
  appended a newline to each, for `while read`, and cut it off again at its end: a run that
  ended early (no off-target, an error, Ctrl-C) left the newline in them, and a read-only
  `--vcf` list failed the write, on stderr, which fails the stage (RV4-1). It reads copies in
  the output folder (`.vcf_list_searched.txt`, `.samplesID_list_searched.txt`).
- **A `--genome` folder or a working directory with a space in its path is refused up front**,
  by complete-search and by build-index-only, each in one line that names the command: the
  search script takes both paths into its commands unquoted, and failed later with a path that
  is no folder's (IDX-19). The working directory is checked first; for the genome folder, the
  line says to move it or to link it from the working directory under a name without one (`ln
  -s '<folder>' <working directory>/Genomes/<a name without spaces>`), which the search takes.
- **The Settings genome upload says what to do with an archive of gzipped FASTA files.** It was
  refused with `no .fa or .fasta file in X`; it now names the `.fa.gz` files it holds and says
  to gunzip them, or to upload one `.fa.gz` a sequence (M3).
- **RELEASE_CHECKLIST's HPC route runs from a clean host environment and a fresh sandbox.**
  Apptainer passes the host's environment into the container (a login shell's `CONDA_*`
  reached it on ml008), and steps 2-4b reused the writable sandbox step 1 had left its caches
  in (list_282.md).
- **The web results page no longer warns of "Multiple perfect matches" for a search of several
  guides with one perfect match each.** It counted the perfect-match sites of every guide
  together, so two guides with an on-target each got the red banner that no on-target can be
  told (RV1-F5). It counts them per guide, as the report's banner does, and names the guide.
- **A 5′ search refused for want of an index names the PAM-less index whose flank is too wide
  for it.** A search skips, without a word, a PAM-less index whose flank is wider than a 5′ PAM
  (the engine reads a 5′ PAM off a flank of its own width only), so the refusal named none of
  the indexes there, or blamed the index's recorded side (cli5 D7). The refusal under
  `--index-path`, the search script's refusal of a variant search, and complete-search's
  pre-flight now say, for example, that `NNNNNNNN_2_hg38 is a PAM-less index whose flank is 8
  nt, and TTTV is a 5' PAM of 4 nt`, and that the search needs `NNNN_2_hg38` or `TTTV_2_hg38`.
- **docs/DYNAMIC_PAM_INDEX_COMPAT.md says what the release does.** Its CLI, web and PAM-picker
  sections called planned what 2.9.0 now does, and described a CLI that took any PAM inside an
  index's own and passed `--prune-pam` by itself, which 2.8.1 did not (IDX-22); they now say
  what ships and how it differs from the design, and that a search no installed index serves
  builds an index of its PAM where it can, as before (it said every one was refused). METHODS
  §9 no longer says that one index answers one PAM, nor that a search takes no all-`N` index of
  another width without the flag.
- **validate-benchmarks' read-only smoke can fail.** It made the code folder read-only with
  `chmod -R a-w` in a container that runs as root, which writes through permission bits with
  `CAP_DAC_OVERRIDE`, so a write into the code folder passed and the step could not fail (CI-1).
  Its container drops the capability (`--cap-drop DAC_OVERRIDE`) and works in a folder of its
  own; RELEASE_CHECKLIST T-38 says so. Run on `pinellolab/crisprme:v2.8.1` with this tree: a
  search script planted to write into its own folder fails the step (`Permission denied`), and
  passes it without the flag.
- **A Settings genome upload of another kind of file is refused on its first chunk, and every
  refusal of an upload is answered as the client's.** A file that is not FASTA, an NCBI
  Datasets `.zip` for one, was taken for a FASTA file and refused once whole, for its first line,
  the archive's bytes, with HTTP 500 `finalize failed: ...`, as was every refusal of what an
  upload holds (WEB-9). An upload whose name is not `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`,
  `.tar.gz` or `.tgz`, or whose first bytes are not what its name says (gzip, or FASTA text), is
  refused before the rest of it streams in, and a refusal at its end (an archive of gzipped FASTA
  files, a truncated gzip, a VCF upload into a registered folder) is answered 400 with its reason;
  so is a personal-assembly bundle that is no archive (`not a gzip file`). A failure of the
  server is still 500.
- **validate-test downloads the brute-force references of its own release.** They came from the
  main branch, and are md5-checked against the registry the release ships, so the day main
  regenerated one, validate-test of every installed release (Bioconda leaves the references out
  of its package) failed (C11/P7). They come from the release's tag (`v<version>`, which
  `crisprme.py validate-test` passes on), or from main for a version with no tag yet.
- **`PostProcess/bench_cbulge.py` does nothing on import.** It changed the working folder and the
  module path of whatever imported it, to a developer's path, or failed to import where that path
  is not (RV4-3); it takes the CRISPR-Bulge repository from `--repo` or `CBULGE_REPO` when run.
- **The report's "Panel: N individuals" counts the samples the search ran with.** It counted every
  list in `samplesIDs/`: the datasets not searched, and the lists Settings writes (M2).
- **build-index-only with no `--samplesID` no longer calls the tiers in place missing.** Run as a
  search's ERROR says to run it, to rebuild an INDELS index, it left the dict-less tiers of the
  build before as they were, and said, in a loud WARNING on stderr, that the index was MISSING
  them, INCOMPLETE (emptyidx F6). It says, on stdout, that they are kept; with none in place, or
  ones a stopped run left unfinished, it warns as before.
- **The web status page says why the search was not run when an INDELS sub-index is empty or
  missing, or the genome folder is not as a search needs it.** Each got the general "remove it and
  try to submit again", though a new submission fails the same way (emptyidx F7, and the genome
  layout refusal of the final 2.8.1 review); the banner now gives the ERROR, and for the genome
  folder what `genome_layout.py` says of it, with the command that fixes it.
- **Settings shows an index's display name.** The Installed indexes list and the Remove card named
  an index by its folder alone, never by the display name given it (`build-index-only --name`, or
  the Display name of a build there), which the search form shows (P15.3): both show it beside the
  folder.
- **Six function-based test modules run under `python -m unittest`.** Each reported `Ran 0 tests`
  and `OK` there, a false pass (suite P4); a `load_tests` hook collects every test function. CI
  and the release checklist still run them as scripts too.
- **Four web test modules skip their classes where Dash is not installed.** They raised SkipTest
  at import, which `python -m unittest <module>` reports as a traceback and an error (skips D7).
- **With no `--bDNA`/`--bRNA`, a search under `--index-path` takes the bulges that library
  serves.** complete-search chose 1 bulge of each type when the genome's FASTA files were there
  to build a reference index from, or the working directory's library held one, though under
  `--index-path` a search reads the reference index from that library alone and builds none
  there: given a library with no index of the PAM, it then stopped on `no matching reference index
  under --index-path` (plan 4.8 RV2-5). It searches with no bulges, which needs no index, and its
  WARNING says to install an index of the PAM in that library, or to run without
  `--index-path`.
- **A search leaves no hidden `.binview/` folder in its results.** The engine reads an index
  folder that holds `.bin` files alone, so a search hands it a folder of links to the trees of an
  index with a sidecar (a download's, and since 2.8.1 every index a search or `build-index-only`
  builds), which stayed in the result folder (RV3-4, RV3-6). It goes as each search ends, whether
  the search found targets or failed.
- **Submit makes no job folder for base-editing values the form does not offer.** With base
  editing on, a window bound that is no whole number, or a nucleotide that is no base, was refused
  (by `int()` and an `assert`) only after the job folder and its files were made, so each such
  request left a folder that never ran, as a Max edits value did before 2.8.1 (secfu RV1-N1). The
  values are checked with the others, first.
- **The Personal assembly tab offers every PAM installed and the app's bulges, whatever the
  Reference tab's genome and variant dataset.** Its PAM list, its bulge options and the note
  under them followed the hidden dropdowns of the Reference tab (OF-11): with a VCF dataset
  selected there that had no index, it offered 0 bulges only, with a note about that dataset's
  index, and with an NRG variant index, NRG alone, though an assembly search indexes each
  haplotype itself; its Submit took the same counts. It offers every installed PAM and 0 to 2
  bulges of each type, its note says what its search builds, and the PAM picker of a PAM-less
  index is not shown there.
- **Settings' Build an index locally pairs a VCF dataset with its genome.** It offered every
  genome for a dataset whose genome Settings recorded when the dataset was added, and a build of
  the dataset on another genome ended in the enricher's crash (`std::out_of_range`, SIGABRT;
  FV :328). Choosing such a dataset selects its genome, the one the form then offers, and a
  build of the dataset on another is refused, naming its genome.
- **A job an identical search is sent to keeps its log.** Submit gave the finished job's "Job
  Done" line today's date (its results are kept 3 days from it) by writing the log back
  without its last line, the report's, and with a second "Job Done" and no line end
  (FV :335); it changes the date of that line alone.
- **The search form offers a variant dataset installed as an index alone under its own genome,
  and by its name.** It read the dataset off the index's name without the genome before the
  `+`, so an index of another genome, of a dataset named after this one, was offered too
  (FR :278), and it showed each underscore of the dataset's name as " + ", as for a merged
  panel of two datasets (`my + vcf` for `hg38_my_vcf`; FR :286).
- **A Ctrl-C ends a command in one line.** crisprme.py ended in a `KeyboardInterrupt` traceback
  (install P6); it says `crisprme.py <command> interrupted (Ctrl-C)`, for complete-search with
  the folder the run wrote into, and still ends by the signal (status 130 in a shell). The
  pipeline's own stages, which the same Ctrl-C stops, leave theirs in the run's
  `log_error.txt`, as before.
- **A source or Bioconda install's self-test prints no warning the image does not.** Every
  anonymous download of the public data printed "You are sending unauthenticated requests to
  the HF Hub", a server warning that huggingface_hub logs and the Docker image hides
  (`HF_HUB_VERBOSITY=error`): crisprme.py sets the same level, unless the variable is set
  already. And `install_from_source.sh` printed git's advice about a detached HEAD as it checked
  out the CRISPRitz tag (install P10).
- **`--sorting-criteria-scoring` is used.** The merge of close targets chose the criteria of the
  CFD and CRISPR-Bulge merges by a variable nothing set, so every merge took
  `--sorting-criteria` and `--sorting-criteria-scoring` did nothing (wset N1). The CFD and
  CRISPR-Bulge merges take it now; its default is what they took, `mm+bulges,mm` (it said
  `mm+bulges`), so a search with neither flag, and every web search, merges as before.
- **build-index-only refuses, before it enriches anything, a VCF the enrichment cannot read.**
  `crispritz add-variants` takes each VCF's chromosome from the last part of its file name that
  holds `chr` and reads the genome's FASTA file of that name, so a multi-chromosome VCF whose
  chromosomes are named without `chr` (`1` on a `chr1` genome, `Contig1`, `NC_000913.3`), which
  the build split under those names, a VCF named so by chance (`nochr.vcf.gz`), and one of a
  chromosome the genome has no FASTA of (another genome's dataset) each ended in `ERROR: cannot
  open <genome>/<name>.fa` and "Genome enrichment (add-variants) failed", the split one after
  splitting the dataset in place (wset N2). Each is one ERROR that says why and what to do, and
  the dataset is left as it was.
- **The Settings uploader shows a refusal in red.** An upload refused (`Error: ...`), one that
  failed, and "Choose a file first." were grey, as the progress and the success are, so a
  refusal read as if all were well (final confirmation of 2.8.1); they take the red of the
  Settings cards' messages.
- **The advice for a merged panel's missing sample list names the list that is missing.** Submit
  makes a merged panel's list (`hg38_1000G_HGDP`'s) from the lists of the datasets it merges;
  with one of them missing, and no VCF folder of the panel, the note and Submit said to add the
  VCF dataset `hg38_1000G_HGDP` again in Settings, a dataset never added (FR :98). They name the
  missing list (`samplesIDs/hg38_HGDP.samplesID.txt`), and how to get it.
- **A genome folder holding a FASTA file of one sequence named otherwise is refused before the
  search.** A search takes each file for the sequence it is named after, so it found the
  off-targets of `walkone.fa` holding `>chr1`, dropped them all, and ended saying that none
  survived post-analysis and to check the folder with `genome_layout.py`, whose check, which
  the search runs first, passed the file; complete-search's input check only warned that
  results "may be dropped" (M1). The check, the search that runs it first (a web search too)
  and the input check now refuse such a `.fa` file, saying to name it after its sequence with
  `genome_layout.py --split`, as build-index-only already did; a `.fasta` file, which no
  search reads, is still only named.
- **`crisprme.py web-interface --port <N>` serves on another port.** The server listened on
  8080 alone, so where another program held it, as on a shared node, an Apptainer or source
  install could not start it (Docker maps another with `-p`; C19). Its start message, and the
  link of a notification e-mail, name the port.
- **`docker stop` stops a CRISPRme container at once.** Its first process, crisprme.py, had no
  handler for the SIGTERM a stop sends, so `docker stop` waited its 10 s and killed the
  container (exit 137; C16): the image's stop signal is SIGINT, which crisprme.py takes as a
  Ctrl-C (the web interface stops its server and exits 0). Measured with the 2.8.1 image, run
  with each signal: 10.2 s against 0.5 s.
- **The Docker image is built on a pinned base.** `FROM mambaorg/micromamba` took whatever image
  was latest the day of the build (C20); the Dockerfile names the tag and the digest of the one
  2.8.1's image was built on (micromamba 2.10.0, Debian 13).
- **A search refused for want of a reference index under `--index-path` names an all-`N` index
  there that records no side, and how to use it.** Such an index of a width but 3 and 8 is not
  taken for a 3′ PAM (a search before 2.8.1 built one of any width from a 5′ PAM file), and
  `build-index-only` records no side over it, so one really built for a 3′ PAM could only be
  rebuilt (R6B RV1-6); without `--prune-pam` the refusal did not name it at all. It names it
  either way, and says that a name of the 2.9.0 form records the side: `mv
  .../NNNN_3_hg38 .../20-NNNN_3_hg38`, for the guide length its trees were built for, after
  which the search takes it.
- **complete-search stops at once when no temporary folder can be written.** With none
  (`docker run --read-only` without `--tmpfs /tmp`), a search ended with exit 0, CRISPR-Bulge
  scores of -1, no report zip and a WARNING in its error log alone (final confirmation of
  2.8.1; TROUBLESHOOTING §1 said so). It says, before it makes anything, which folders it
  tried, and to set `CRISPRME_TMPDIR` to a writable one or give the container a `/tmp`.
- **complete-test is not stopped by the empty results folder a refused search leaves.**
  complete-search makes its output folder, and a refusal of its pre-flight leaves it empty: the
  next complete-test stopped on it, "Complete-test already run once", until the user removed it
  (R2-F2). An empty folder holds no earlier run's results, and complete-search takes one: only
  a folder with something in it stops complete-test now.
- **The base-editing window of BED lines answers no server error.** The form recomputes the
  window's options at each change of the guides, and for genomic coordinates it raised with no
  genome selected (a fresh data folder after Load Example), on a chromosome the genome has no
  FASTA of, on coordinates past its end and on a line short of its three fields: each an HTTP 500
  and a callback error in the browser's console (FR :89, :189). It offers no window for them, as
  for any line Submit then refuses, saying why.
- **`complete-test --chrom` fetches a chromosome but chr22 from Hugging Face too.** The data repo
  holds every hg38 chromosome gzipped (`genomes/hg38/<chrom>.fa.gz`) and chr22 unpacked as well,
  and complete-test asked for the unpacked file alone, so every other chromosome came from UCSC
  (FR :725). It takes the gzipped one from the Hub when the unpacked one is not there, and UCSC
  stays its source when neither is: chr21 came from the Hub, 47,644,190 bytes, as from UCSC.
