# [2.8.1] entries from the tests track

All five go under **### Fixed** of `## [2.8.1]` in CHANGELOG.md.

### Fixed
- **A search with no samples no longer writes into the install's empty `vuoto.txt`.**
  `complete-search` without `--samplesID` (a reference-only search, or a variant-aware one
  without it) hands the pipeline the shared empty placeholder `PostProcess/vuoto.txt` as its
  samples file. The pipeline added a newline to it as it started and removed it only at the
  very end, so every such run that ended earlier left one more newline in it for good: a
  search that found no off-targets, one that failed, or one stopped with Ctrl-C (measured
  through the CLI: 0 bytes, then 4 after one failed run and three with no hits). This
  happened on every install the run could write to: Docker (as root), source installs, and
  Bioconda, where the file is a hard link into the conda package cache, so the stray bytes
  also reached any env later created from that cache. A read-only install was never written,
  and a search started from the web interface uses a samples file of its own, so it never
  touched the placeholder. No result changed: the pipeline skips blank lines, and the other
  places that use the placeholder check its name, not its contents. Reference-only and
  variant-aware searches without `--samplesID` report exactly what they reported before
  (every result file compared, except the plot PDFs, which differ between any two runs). The
  pipeline now leaves the placeholder alone, telling it apart by its name as it already does
  for the annotation placeholder. An install that already has the extra newlines needs
  nothing; to empty the file anyway, run
  `: > "$CONDA_PREFIX/opt/crisprme/PostProcess/vuoto.txt"`. Present in 2.8.0 too. Covered by
  `PostProcess/test_index_resolution.py`.
- **`scorer-env check` no longer lists `pkg_resources=?`.** `pkg_resources` has no version of
  its own, so a healthy scorer env showed `pkg_resources=?`, which read like a missing
  package. It now shows the version of setuptools, which ships it (for example
  `pkg_resources=80.10.2 (setuptools)`; CRISPR-Bulge needs setuptools below 81 for it), or
  `ok` where no version can be read. The same value goes to `Annotations/.scorer_env.json`.
  The status, the other versions and the compute self-test are unchanged. Covered by
  `PostProcess/test_scorer_env.py`.
- **`validate-test` no longer prints pandas `FutureWarning`s, and runs under pandas 3.** Each
  benchmark it compared printed `FutureWarning: Series.__getitem__ treating keys as positions
  is deprecated` (two in the documented chr22 self-test): the brute-force reference's columns
  were read by position from rows labelled by column name. pandas 3 removes that fallback, and
  there `validate-test` stopped with `KeyError: 0`; the shipped environments keep pandas below
  3, so a released install saw only the warnings. The columns are now read by position
  explicitly. The comparison, its report and its exit code are unchanged. Covered by
  `PostProcess/test_validate_test.py`.
- **A source install no longer leaves an untracked `Annotations/` folder in the checkout.**
  `install_from_source.sh` runs `scorer-env create` and `scorer-env check` in the checkout,
  and they save the scorer env's state to `Annotations/.scorer_env.json` there. The
  checkout's `.gitignore` covered every other data folder that CRISPRme creates in the folder
  it runs in, but not `Annotations/`, so `git status` showed it after every install. It is now
  ignored (the top-level one only: `test/data/Annotations/` holds tracked test files). Nothing
  to do; the file is a diagnostic record. Covered by `PostProcess/test_install_from_source.py`.
- **The unit tests pass inside the shipped image.** The build-lock tests in
  `test_index_resolution.py` listed processes with `ps`, which the image does not have (no
  procps): run there, the module ended with 44 errors and 6 failures. They now read `/proc`,
  as the pipeline's own lock code does (`ps` only where there is no `/proc`), and pass in the
  image as on CI. `test_crispr_bulge_fidelity.py` changed to a developer's folder as it was
  imported, so `python -m unittest test_crispr_bulge_fidelity` failed everywhere. It now skips,
  saying why, where TensorFlow or its dataset is absent, and in the cbulge env with the
  dataset it runs the fidelity check as before (exit 0 on a match); its docstring says how to
  fetch the one dataset file it needs. Searches are not affected.
