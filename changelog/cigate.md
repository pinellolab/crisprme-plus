<!-- cigate track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the docs track: amendments to one entry already in `## [2.8.1]`, then new entries.** The
entry "The documented install self-test tested nothing, and `validate-test` then failed with a
traceback" (`### Fixed`, CHANGELOG.md lines 360-384 at a5379e0) gives the commands and counts this
track changes. Replace the quoted text with the text after it, after marks' amendment 5 to the same
entry, whose "their variables go" then becomes "its variable goes".

1. Replace

   > installer now give the commands CI runs: `CRISPRME_ALLOW_ONDEMAND_BUILD=1 CRISPRME_SKIP_HEAVY=1
   > crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G`, then `crisprme.py validate-test
   > --chrom chr22`.

   with

   > installer now give the commands CI runs: `CRISPRME_ALLOW_ONDEMAND_BUILD=1 crisprme.py
   > complete-test --chrom chr22 --vcf_dataset 1000G`, then `crisprme.py validate-test --chrom chr22`.

2. Replace

   > covers chr22 only (the Cas12a one and both single-"n edits" ones) is now skipped on any other

   with

   > covers chr22 only (the Cas12a one and the single-"n edits" ones) is now skipped on any other

3. Replace

   > (223 extra sites) where it now ends `1 passed, 0 failed, 3 skipped`. Validation on chr22 is
   > unchanged. In Docker, the only route on macOS, the two variables must be passed with `-e`: set
   > in front of `docker run` they stay on the host, and the run tests nothing. README §1.2 and

   with

   > (223 extra sites) where it now ends `1 passed, 0 failed, 5 skipped`. Validation on chr22 is
   > unchanged. In Docker, the only route on macOS, the variable must be passed with `-e`: set in
   > front of `docker run` it stays on the host, and the run tests nothing. README §1.2 and

### Changed
- **The documented self-test now runs every benchmark, and it ends `6 passed, 0 failed, 0 skipped`.**
  **Affected:** anyone who runs `complete-test`'s benchmark mode and `validate-test`, as README
  §1.2, §1.3, §2.2.2, §2.2.3 and §3.2.1, `docs/DOCKER_QUICKSTART.md`, both commands' `--help`, the
  installer's last line and `validate-benchmarks.yml` give them. All of these set
  `CRISPRME_SKIP_HEAVY=1`, on the claim that the bulge-2 index builds need more than 16 GB of RAM,
  so each left out 2 of its 4 benchmarks. That claim was wrong: the run peaks at 3.66 GB of RAM,
  and it passes in a container held to 16 GiB, 4 CPUs and no swap. **What changes:** the commands
  no longer set it, and CI no longer does, so every benchmark runs, and CI fails if one did not.
  All six took about 25 minutes at 8 threads on an Apple Silicon Mac in testing, where the two
  light ones took about 12. **What to do:** nothing. The variable stays as an opt-out for a very
  small machine: `complete-test` says what it leaves out, and `validate-test` ends `2 passed, 0
  failed, 4 skipped`. A folder that ran the old command has no output for the others:
  `validate-test` there ends the same way, saying `complete-test targets not found`. Covered by
  `PostProcess/test_validate_test.py`, `PostProcess/test_complete_test.py` and
  `PostProcess/test_ci_workflows.py`.
- **Two new benchmarks that one lost site fails, `cas9_sg1617_n5` and `cas12a_hbg_n5`.** The two
  bulge-2 benchmarks could hardly catch a lost site: `cas12a_hbg_simple`'s brute-force reference is
  a header with no rows, so that case compares empty with empty, and `cas9_sg1617_simple`'s holds 3
  rows at 2 sites. The same guides at `--max-total-edits 5`, on the same indexes, have 381 and 92
  sites on chr22, and CRISPRme matches brute force on them, on CRISPRitz 2.9.1 and on the 2.9.0 CI
  overlays. No result was lost: this was a gap in the test's power. Their references are in
  `test/benchmark/max-edits-5/`, which Docker, source and Bioconda installs all hold, so they need no
  network. Covered by `PostProcess/test_validate_test.py`.
- **CI runs what it says it runs** (contributors). `unit-tests.yml`'s job limit goes from 10 to 20
  minutes: the hosted job took 523 s at the first release head of 2.8.1, and the steps have grown by
  a quarter since, which predicts about 640 s, past the limit (a prediction: no job was run to
  measure it). `web-e2e.yml` gets a step that runs, in the image, the tests that skip in the light
  job for want of BEDOPS, minimap2, impg or Dash: `test_assembly_annotate`, which no workflow named,
  the impg truth tests, the flag-gate behaviour tests, the `post_process.sh` test, whose own step in
  the light job tested nothing, and `test_index_resolution` (368 tests, about 100 s). It fails on any
  skip but the two the image cannot help. `validate-benchmarks.yml`'s smoke step runs the real
  Hugging Face download test, which no job ran. The path filters of `unit-tests` and `web-e2e` now
  start, on a change to any file one of their tests reads, a job that runs that test, as found by
  replaying each step with every file read logged: `unit-tests` lists every file its tests read,
  and `web-e2e` every file read by the tests only it runs (the docs and scripts its in-image step
  reads through `test_complete_search_flags`'s scans start `unit-tests`, which runs the same
  scans). Before, a pull request that changed only `pages/`, `scripts/`, `docs/`, the README,
  `test/benchmark/`, `test/data/`, `assets/`, or a script outside `PostProcess/` such as one in
  `seq_script/`, did not start `unit-tests`, and one that changed only `PostProcess/`,
  `seq_script/`, `PAMs/` or `download_hprc_assembly.py` did not start `web-e2e`. Covered by
  `PostProcess/test_ci_workflows.py`, which `unit-tests.yml` runs.
- **`docs/RELEASE_CHECKLIST.md` Gate B gives the local gate** (maintainers): the exact commands that
  leave nothing skipped, on the release image, amd64 and arm64, with the three skips it allows
  named. The six function-based test modules run as scripts there, since `python -m unittest` alone
  reports `Ran 0 tests` and `OK` for them.

### Fixed
- **Tests that failed where nothing was wrong** (contributors; no search or report changes).
  `test_build_index_only_atomic`'s setgid test errored in any unprivileged Apptainer run, where
  `os.getgroups()` gives 65534 for every group the user namespace does not map: it now skips there,
  saying so. `test_generate_report.TestAssemblyExtraColumns` failed 3 of its tests when a combined
  report had been built before it in the same process, as under `python -m unittest discover`: it
  no longer depends on that. `test_complete_search_flags` and `test_post_process_tempdir` failed on
  a read-only code directory such as a SIF: they now pass there, and the first no longer leaves the
  data folders beside the code on a writable install.
