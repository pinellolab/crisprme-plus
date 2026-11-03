<!-- docs3 track: entries for CHANGELOG.md [2.8.1], to be merged by the docs track. The last
section holds notes for the merge: an amendment to an entry already in [2.8.1], the
corrections this track made in CHANGELOG.md itself, and where the EQUIVALENCE_NUMBERS
placeholders the release manager fills are. -->

### Fixed

- **`--personal_annotation` reads a bgzip- or gzip-compressed BED, as its `--help` asked.** The
  help said the file "must be compressed using bgzip", but since 2.2.0 the personal BED has been
  read as plain text, so a bgzipped one reached the sort garbled and every search given one
  stopped before it began, with or without `--annotation` (`No tabs/spaces found at line 1 ...
  Error: Sorting BED file failed`, exit 1; the message never named compression). A compressed
  file, told by its content and not its name, is now decompressed first, and a damaged one stops
  the search with one line that says so. `--help` and `docs/INPUT_FORMATS.md` now say that a
  plain, bgzipped or gzipped file works. Measured through the real CLI on a toy genome: a
  bgzipped or gzipped personal BED, alone or with a plain or bgzipped `--annotation`, now gives
  the plain BED's results byte for byte, where it stopped at exit 1; with a plain BED the results
  are unchanged, byte for byte. Covered by `PostProcess/test_annotation_tmpdirs.py`.
- **`complete-search --help` and the docs said a search without the CRISPR-Bulge scorer env
  completes on CFD only.** Since 2.7.3 `complete-search` refuses to start while the `cbulge` env
  or its weights are missing or broken, naming the fix (`crisprme.py scorer-env create`;
  `crisprme.py scorer-env doctor` says what is wrong), and runs CFD-only, with the CRISPR-Bulge
  columns at `-1`, only under `CRISPRME_SKIP_SCORER_ENV=1`, as README §1.1.2 says. The
  `--compute-backend` help, its README entry, the README's scorer paragraph,
  `docs/TROUBLESHOOTING.md` (which added "rather than failing") and METHODS now say so, and that
  `targets-scoring` and a search started from the web interface do not check first: with the env
  absent they complete CFD-only, with the CRISPR-Bulge column at `-1` and a banner in the log.
  Nothing else changes. Covered by `PostProcess/test_complete_search_flags.py` (the help) and
  `PostProcess/test_annotation_tmpdirs.py` (the refusal).
- **The installers and the docs give the measured image and install sizes.** The one-line
  installers said the engine image is ~800 MB: `docker pull pinellolab/crisprme:latest`
  downloads 4.84 GB on amd64 (Windows, Intel Macs, most Linux hosts) and 1.81 GB on Apple
  Silicon, summed from the registry manifests, and unpacks to ~7 GB or more. They now say a
  ~5 GB download, ~1.8 GB on Apple Silicon. They also asked for a data folder with ~45 GB free,
  which is the download; the README, the website, `docs/DATA_LAYOUT.md` and the Quickstart said
  ~85 GB on disk and ~100 GB free, and `docs/TROUBLESHOOTING.md` ~3 + 40 + 50 GB. A default
  variant-aware install (`download --what all`, then
  `download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`) measures 97.9 GB, and peaks
  around 103 GB while the genotype store unpacks beside its archive, so ~100 GB free was full
  before the first result. Every one of them now says ~98 GB on disk and ~150 GB free, the
  Quickstart's breakdown lists the indel genotype store and the indel logs (14 GB it left out),
  and `install/README.md` names the figures to keep in sync. Covered by
  `PostProcess/test_documented_sizes.py`, a new step of `unit-tests.yml`.
- **The `lucapinello/crisprme-data` card said every index there works with older clients.** Its
  overview said so, while the same card, `docs/INDEX_PROVENANCE.md` and this changelog say a
  2.8.0 search of the rebuilt `NNN_3` archives fails or reports malformed rows. It now says any
  release can read them and to use the rebuilt `NNN_3` archives with 2.8.1 or later. The live
  card is corrected when the release uploads `docs/hf/crisprme-data_README.md`.
- **`docs/DYNAMIC_PAM_INDEX_COMPAT.md` read as a description of 2.8.1.** The page, new in 2.8.1,
  is the design written before `--prune-pam` was built, and described what does not ship: a
  `complete-search` that adds `--prune-pam` itself when the PAM is narrower than the index, a web
  PAM picker, an `index_meta.json` sidecar on every PAM-less index, and a `pinellolab` repository
  for PAM-less indexes. A status note under its title now says what 2.8.1 does: the flag is
  explicit, needs `--index-path` and is command-line only, and only `complete-test` passes it
  itself, to search a 3′ all-N index wider than 3 nt such as the K = 8 one; an index serves only
  its own PAM or, all-N, a narrower one; no `index_meta.json` is written or read (the PAM file
  is recorded in a `.pam_build` by `build-index-only` and by a search that builds an index);
  `NNN_3` is in `lucapinello/crisprme-data` and the 8-nt-flank indexes in
  `lucapinello/crisprme-data-general-index`. The sections that describe the rest are marked
  planned, and its example of the engine's subset rule no longer has `NRG` serve `NGA`.

### Changed

- **Release docs (maintainers): the CRISPRitz pin of the Bioconda recipe is bounded, and the
  version sweep compares like with like.** `docs/RELEASING.md` and the `release-crisprme` skill
  ask for `crispritz >=X,<Y`, below the next CRISPRitz that changes the index format or the
  enricher: `>=2.9.1,<2.9.2` for 2.8.1, which reads index format 1, where the CRISPRitz planned as
  2.9.2 writes format 2 and enriches multi-allelic sites another way. With a bare `>=` the solver
  would install that engine beside 2.8.1, which cannot read the indexes it builds:
  `build-index-only` fails, and a bulge search refuses such an index when it is reused. The
  checklist's version sweep asked the Dockerfile's `crispritz_ref`
  (v2.9.1) to read CRISPRme's version (2.8.1), which no correct release could pass; it now asks it
  to equal `install_from_source.sh`'s `CRISPRITZ_REF` default, which
  `PostProcess/test_install_from_source.py` checks.

### Notes for the merge

- Amend the entry "`docs/RELEASE_CHECKLIST.md` Gate B gives the local gate" (Changed) with: "Its
  step 4b runs, on the release image, amd64 and arm64, the engine checks CI runs only with the
  branch overlaid on the published image and its older CRISPRitz: `test/index_flex/
  check_index_flex.py`, whose K = 8 cells must give the same results as a dedicated index (on an
  engine without `--prune-pam` they pass by refusing), `test/indel_snp_oracle/run_oracle.sh` and
  `scripts/check_alt_pam.py` on step 4's tables. `PostProcess/test_ci_workflows.py` pins the step
  and its greps to what the checks print." (Trap T-42 in the checklist.)
- Already corrected in CHANGELOG.md on this branch, against what the code does at 9c6b4046: the
  index-choice rules of "When a dedicated index for the PAM is installed beside a PAM-less one"
  (the guide-length rule, second of four); the no-bulge sizing of "The RAM thread guard sized a
  search served by a PAM-less index" (a variant search of a download is sized by its variant
  index); "A search of a VCF dataset built in the web Settings page reads its dict-less tiers"
  (the SNP + indel row is the same on the dictionaries since the indel fix; no extra
  `Other_motifs` alignment or new row order on the tiers; a dataset not rebuilt changes too);
  "`scorer-env create` failed whenever conda was the env manager" (a `complete-search` refused, a
  web search ran CFD-only); and the `--prune-pam` equivalence in Added, the header-injection entry
  and the `NNN_3` republish entry (no release-candidate counts).
- The placeholder `EQUIVALENCE_NUMBERS`, in double braces inside an HTML comment beside the claim
  (grep the repository for it), marks where the final genome-wide
  three-arm run's numbers go, once the last results fix is in: the `--prune-pam` entry under
  Added, METHODS §9's "Equivalence" paragraph, and the "Equivalence claim" of
  `docs/hf/crisprme-data-general-index_README.md`. The claim reads correctly without them.
