<!-- defaults track (n290): CRISPRme 2.9.0 entries for CHANGELOG.md. NRR as the default PAM and
index, the index-flexibility items of plan 2.7 (E1), index names that carry the guide length and
the side (E2), and the clear fixes of the backlog. -->

**For the merge, not for the changelog.** The install sizes below (150 GB on disk, a 161 GB
peak, a 70 GB download) were measured on ml003 on 2026-11-06 with the two documented downloads,
twice, the peak by sampling `du -sb` every 10 s. Every published NRR archive is the 2.8.1 one: the
rebuild lane's archives change these sizes again (by 1-2 %, plan 2.5), so re-measure with them,
and change `PostProcess/test_documented_sizes.py` and every surface it reads together.

**For the merge, not for the changelog (E2 and the rebuild).** `build-index-only` names an index
for its guide length and side (`20-NRR_3_hg38+hg38_HGDP1kGP`), and `publish-index` publishes it
under its folder's name. The rebuild's archives, built with this release, take that name unless
renamed at staging; every command and document that downloads by the published name
(`download --index-name NRR_3_hg38+hg38_HGDP1kGP`, the installers, the self-test) then needs
`download` to find an archive of the new name for the older one, or the rebuild to publish under
the older name, which every reader takes: the rebuild's choice (plan 2.5), with plan 3 C2's new
location.

**The backlog: what this track fixed, what it found fixed, and what it leaves, and why** (every item
of `_r281/r109/list_282.md` and of the plan's section 4, for the merge; the ids are the sources', a
JSON line where the plan cites one, as FV :335).

*Fixed here, each with its test (the entries below):* IDX-17 (a), (c) and (e), and (b) by E1 (a
search with no bulges reads no index); IDX-19; IDX-21; IDX-22; cli5 D7; OF-11; OF-12; OF-13; FR :98;
FR :125; FR :89 and :189 (the base-editing window's HTTP 500s); FR :278; FR :286; FV :328 (the
Settings build form's genome list); FV :335; F2 (the web's bulge-cap side rule); RV1-F5; RV4-1;
RV4-3; M1; M2; M3; CI-1; WEB-9 and the bundle that is no archive; C11/P7; C14; C16; C19; C20;
emptyidx F6, F7 (with the genome-layout refusal's banner) and F8; suite P4; skips D7; P15.3; install
P6; install P10 (the Hub's warning, git's advice); install P9 (pigz); R2-F2; R4-4 (the unused index
of an assembly search with no bulges); R6B RV1-6, RV2-5, RV3-4 and RV3-6; LP
unterminated-last-line-dropped; M6; wset N1 and N2; wsearch N8; RV3-D; secfu RV1-N1; the chunked
uploader's grey refusals and the search with no temporary folder (the final confirmation of 2.8.1);
setup's exit status, progress lines and URLs (FR2 cli a, c, e); FR :725 (complete-test's chromosomes
from the Hub); the checklist's HPC route (`--cleanenv`, a fresh sandbox); the installers' first-run
chain (plan 4.4); the Personal assembly tab's bulge limits (the final confirmation of 2.8.1, with
OF-11); R6B direct-url-fasta-and-hidden-files (plan 4.8: a Direct URL archive is laid out as an
upload); of the focused review of f8dbe73 (the last bullet of list_282.md's final confirmation on
bfddb59), INPUT_FORMATS §4's first lines (the length alone, a third field, a blank line: refused as
it says) and CRISPR-Bulge's autograph `.pyc` files, left in the temporary folder by every scoring
run; the stale test comment of plan 4.5's coverage gaps (test_web_bulge_notes, FR's
W3-none-selected-changed: a comment, with no test of its own); of plan 4.8's docs and wording,
RV1-WEB2-3 (the web user guide's ALT rows), RV1-WEB2-5 (the web's refusal of a VCF folder with a
file named for no chromosome beside others), D3R3-2 (a CFD-only run's CRISPR-Bulge risk score),
D3R3-3 (targets-scoring's banner, on stdout), LP RV4-cooc-carrier-cap (the SNP+indel companion names
50 carriers at most) and D3R3-4 (a test docstring's count of the release candidate: a docstring,
with no test of its own); WEB-4 (the Graphical Reports selector's totals); of plan 4.2's report
items, RV1-F1 (a search of two or more datasets), RV1-F2 (a haplotype's report) and cli5 D6 (a run
with no score); FR :163 (download's check of an index's dictionaries); R2-F1 (an INDELS index built
while the fake-indel genome changed); R4-6 (a haplotype result reused after its FASTA files
changed); M5 (by E2's builds in the build lock: a search puts an index under its name only whole,
its `.bin` files in it, so the ERROR of a folder with no `.bin` can no longer name one a search is
building, as `test_index_resolution`'s
`test_a_build_leaves_nothing_under_its_name_and_touches_no_index_of_the_older_one` shows; a search
of 2.8.1 or earlier, which builds in place, run in the same library at once, still can); wsearch N3
(the plots' perfect-match legend with no star); FV :342's first two, whole: (1) the per-sample
option of a dataset not indexed yet, and (2) both its halves, the sites-only note and Submit's
refusal that cited `docs/DESIGN_mega_index.md`, with every piece of the report's advice that
named a command a web search cannot pass: both sentences that say to re-run the search
per-sample (the Search mode row and the caption of the variant combinations), the commands that
get the tiers, and every sentence that names `--report-full` (the variant-created list and
subset, the raw integrated results, the dense windows, the observed-haplotypes companion the
slim zip leaves out); its other three are R4-4's names, P14 and N7, each recorded here. Of the audit's test
additions (plan 4.5), the leftover INDELS folder (build-index-only's "Indels index already present"
of an INDELS index without its sub-indexes: a behaviour, not a test). Of R6B: the cigate items F1
(the scan guard that saw no skip of a method), F2 and F3 (the wiring guards a YAML comment
satisfied, and the logs of a job its time limit cancels) and F4 (the SKIP_HEAVY wording), with FFR
RV1-WEB2-6 (a test helper's open file); emptyidx RV1-2, RV2-2 and R3-1 (T-40's "always a broken
build"). Plan 4.1's FUTURE_WORK figure (~17–19%, a raw count, beside METHODS §8's).

*Found fixed already (checked on this branch):* P8 (fixed in 2.8.1: the tier path writes the REF
protospacer's mismatches in lower case, as the dictionary path does, CHANGELOG [2.8.1], covered by
`test_tier_ref_columns`, which unit-tests.yml runs); P12 (exit 2 since 2.8.1); t31 P6 (C3 of 2.8.1:
every folder named, exit 2, and the banner says reference when the smoke runs reference-only);
R2-F3; R4-5; WEB-8; of plan 4.5's test additions from the audit, the smoke test that knows the
geometry (2.8.1's I7, with IDX-9, OF-9 and cli5 D8: the prebuilt-index smoke searches each index's
own PAM on its side, the guide padded for it); cli5 D5 (= IDX-7); IDX-20 (2.8.1 runs a store-less
install population-level and says why); F3 of the docs review (a wrong-side all-N variant index is
named by the search script's refusal since 2.8.1); F6 of the docs review (DYNAMIC_PAM §1 says a 3′
PAM); install P4 (Biopython's warning); install P10's pandas FutureWarning (validate.py reads by
position, run in the image with `-W error`); the environment.yml web pins; the self-test's 404 lines
(t31 P9's noise); RV4-4 and R4-7 (web-e2e's in-image step runs test_index_resolution; no `ps` there
is its one allowed skip); text N1 to N3; P14 (the results page's Group by opens on CFD_score); the
older list's first-paint label (Gate B walk 2b and the 38a7c4e walk's step 9 passed); the stderr
fixtures of a runner with less memory than memory_check asks for (list_282.md, "Found 2026-11-04":
CRISPRME_SKIP_MEM_CHECK): with memory_check reading the 15.6 GB of GitHub's runner, so that every
complete-search and complete-test that does not skip it warns, every step of unit-tests.yml passes
in the CI image, and web-e2e's in-image tests in the image (b648227's fix of the one fixture that
compared it was the one needed). And fixed in 2.8.1: `/load`'s Show full log (`show_full_log` and
`show_assembly_full_log` take the job folder through `results_job_dir`); F4 of the docs review (the
release skill says that the only pin is `crisprme.py`); the [2.8.0] date (2026-10-24); the
CHANGELOG's unit-tests.yml 20-minute sentence (it gives GitHub's 14 min 50 s and 13 min 55 s); the
1e-5 wording (the README, METHODS and the report say that no floor is applied);
DOC-PAMCREATE-LEFT-COUNTS (the PAM_creation entry counts the other 55, 186 in all); M4 (the entry
names the hidden samplesIDs files; whether its disclosure goes under Fixed or Changed was the docs
track's call); C18, the FutureWarning (install P10's, above); LP RV4-3 (its case, a sites-only
dataset searched with a phased one, runs per-sample since 2.8.1, with no population-level banner; of
sites-only datasets alone, the banner's advice is for a genotyped panel, which the review found not
wrong). D3R3-7 is moot: no `{{EQUIVALENCE_NUMBERS}}` placeholder is left in the tree to fill.

*Left: design questions, or the PI's or an owner's call.*
- Plan 4.1, every item but FUTURE_WORK's figure (fixed, above): two cis indels in one protospacer;
  CRISPR-Bulge's batch noise; the one-row drops; phasing against trio truth; phase sets; #77; the
  calibration edge; the HPRC and HG01255 runs. Each is science or a measurement, not a code fix.
  And its documented limits, which the plan lists with no fix planned: the code-fix option it
  gives the fourth, `resolve_carriers` following phase (FFR OLD-UNPHASED-OVERRIDE), would change
  who an indel row lists on a partly phased panel, an option the plan did not take.
- IDX-17 (d) (a dataset downloaded to another folder: the search reads its dictionaries, sample
  list and VCF folder from the working directory, and the input check says its VCF folder is not
  there): `--index-path` covering a dataset's data is plan 2.7's "medium" item, not taken by E1.
- R6B RV3-7 (a variant search under `--index-path` with no variant index in either library builds
  one in the working directory): R6B leaves it to the owner, the docs describe what it does, and
  refusing would stop a search the dataset's VCFs can serve (plan 2.7's `--index-path` row,
  "medium").
- FR :171 and FR :207 (a dataset removed in Settings is offered only if its name begins with its
  genome's, and a stale page then names the wrong one): 2.8.1 documents the rule; offering more
  removed datasets would offer searches Submit refuses (their list went with them).
- WEB-7 (the analysis mode on the results page, History and the Risk Cards: plan 2.6 designs what
  each shows, with the carrier gate's display), WEB-11 (Hide Reference Targets ticked by default
  since 2021: a default); P7 (what the
  slim report shows); the 0-1-edit sites of the shareable report and N7; the CRISPR-Bulge of
  reconciled assembly rows; the annotation copies of a web search (the #97 trade), and Submit's
  stall while it makes them; Settings builds without the tiers' RAM, and two things that follow,
  each documented in 2.8.1 (CHANGELOG, SETTINGS_DATA_MANAGER): a dataset built there before 2.8.1
  has no tiers until it is built there again, and a sample list put there after its build reaches
  them only once its registry and genotypes are removed and it is built again (a build that
  replaced tiers in place by itself would take the RAM above unasked); W1 (deliberately not
  added); FR :107 (deliberate); the stderr gate of the pipeline stages (exit codes instead: 44
  gates); the build locks of plan 4.4, each 2.8.1's design, whose CHANGELOG says why and what to
  do: a lock held from another host or container is taken over only after 15 minutes without a
  refresh, a 2.8.0 lock, which has no owner record, is never broken by age (a 2.8.0 run may still
  be building under it), and a build directory left in another PID namespace or on another host
  is not reclaimed (its owner cannot be checked from here); the undetected damage of
  earlier releases (each needs a stamp those releases never wrote, and each CHANGELOG entry gives
  its manual fix); the genome-folder name collisions of plan 4.3; source installs enriched before
  CRISPRitz 2.9.1 (documented, re-enrich); FV :180 (the Hub's cache keeps a copy in $HOME: it spares
  a re-download, and TROUBLESHOOTING says HF_HOME moves it); the dash interrupt of plan 4.4 (a
  signal to crisprme.py alone, which leaves the pipeline running and keeps its annotation
  temporary folders)
  and c3's signal nits (a process-scan race, SIGHUP to the caller's
  group, two signals): c3's signal design (the image now stops by SIGINT); `targets-scoring` and web
  searches without the scorer env (documented in 2.8.1: refusing them, as complete-search does,
  or warning on the web form, decides what an install without that env can still do, CFD-only
  searches); a variant search of a `Genomes/` link of
  another name (which name a variant genome takes, the link's or its target's, with installs holding
  either); R6B RV1-5 (a partial `--samplesID` list drops carriers: "owner's call"), dash-alone-id
  ("none needed"), f4-wrong-list-named ("owner call"), renamed-ids-stale-index and
  refused-download-drops-installed-genome ("for the owner to decide"); R4-4's hashed haplotype index
  names in Settings and History (each haplotype's own index, named so since 2.8.1, one name a
  haplotype: whether those lists name the individual instead is a choice of what they show); the
  results pages that open only
  job-id-form folders and the region page's refused chromosome ids (by design, the security round 2
  decisions b and c); setup re-downloading the annotations, with the setup guide's A2e (rv2-3: which
  form setup keeps, bgzipped as a search reads or plain, is one choice the guide then follows).
- The served-mode names of plan 4.3 (a genome, PAM or dataset name placed by hand in the data
  folder, or through Settings before 2.8.1, reaches the pipeline's shell commands): plan 3 A1
  leaves the PI to choose between restricting the names a search takes and recording the limit.
  Such a name is the operator's: Settings, which refuses one since 2.8.1, is off on a served
  instance.
- The classic path's MAF of a row one of whose variants has no frequency (FR :438, :531): 2.8.1
  documents the rule, in the report's MAF note (a variant with no frequency is left out of a
  row's MAF, which is blank when none has one); a blank MAF, or a joint one, for such a row is a
  choice of what the column means.
- The IUPAC letters of a sites-only (mega) row's DNA: no haplotype resolves its codes into one
  allele, and printing one would print a sequence no individual is known to carry; how such a
  row shows its alleles is plan 2.6's degeneracy flags, a design.
- Contigs with "_" in their name (bio11.3, and around it R6B emptyidx's RV1-4, RV2-3 and R3-2 to
  R3-4, plan 4.8's list: no word at search time that the bulged indel off-targets of such a
  contig were not searched, the raw shell errors of a search without its fake-indel genome, the
  crash of a variant there in the integrator, and the "Indexing INDELs in KI270302v1" line of a
  build; and the empty `_INDELS` of a dataset whose every indel is on such contigs, which
  download refuses, as T-40 now says, above): every reader of the fake-indel genome
  (the INDELS build, the SNP overlay, the indel post-analysis, the search script) takes the
  chromosome off the engine's folder name `fake_<vcf>_<chrom>` after its last `_`; whether the
  engine's naming changes, with index format 2, or every reader does, is a design question.
  INPUT_FORMATS documents the rename that avoids it, and no published index has such a contig.
  (R3-5 is emptyidx F6, fixed here, and R3-6 a wording of 2.8.1's CHANGELOG, below.)
- INT N2 (an indel-free dataset installed without its enriched genome: every search says to
  download its INDELS index again): telling a dataset with no indel from an INDELS index lost or
  never built needs a record of the dataset's indels that an install may not hold, a design of
  the INDELS check.
- Plan 2.7's roadmap beyond E1 and E2 (list_282.md's index-flexibility list): `index_meta.json`,
  serving other guide lengths (E3 decided: pair them), the 5′ wider-flank prune (a CRISPRitz
  change), migrating an old fake-indel genome, a stale-index stamp, a web setting for a shared
  library, per-sample on a tierless index from the web, other PAMs on a variant panel, a 2.8.0-built
  5′ all-N index, the builder version on PAM-less indexes: features the E1 decision did not take;
  the 5′ PAM-less flank-4 indexes are the rebuild's (E4, E5).
- Of the focused review of f8dbe73, the two items not fixed. A read-only home's Hugging Face cache
  (`HF_TOKEN_PATH` in place of `HF_HOME`): where the Hub's cache goes when the home cannot hold it
  is a design question, as `hf_hub_download` keeps each file whole in that cache before `hf_fetch`
  copies it, so a temporary folder removed at exit would hold each download twice, in /tmp or
  CRISPRME_TMPDIR, while 2.8.1 keeps a saved login found and falls back to the file's original
  source, which works. SIGTERM leaving matplotlib's temporary folder: a SIGTERM
  handler changes how every command, and the pipeline under it, stops (c3's signal design, above),
  and TROUBLESHOOTING §1 says that a command killed by a signal leaves it in the temporary folder.
- list_282.md's parked design work: the two-pass dense search (its min-edit lower bound breaks
  for two or more cis indels and equal-length MNVs: the density track's) and last-ALT-wins at
  split multi-allelic sites (B1 decided that every allele is coded: the alleles and engine
  tracks').
- wsearch N1 (the `/result` region route answers HTTP 500): not a path fix. Fixing the relative
  path, tried here on a real job, stops the 500 but shows no row, as the route reads the columns of
  a pre-2.x bestMerge and an annotation file no search writes (an awk error too); removing the
  route, which nothing links to, or rewriting it on the integrated results is a choice. N5 (a
  finished job whose last steps warn shows the failure text until its log's rename, its last
  step, as its status page does): whether a log of warnings alone reads as done on both pages is
  a choice.

*Left: the work of this release's other lanes, the rebuild and the integration.* A6 (the rsID
and MAF lists of a multi-variant row, sorted apart from its variants) and results F2 (a SNP row
that lists a carrier whose haplotype holds an indel in the target): fixed by the alleles and gate
lanes (`changelog/alleles.md`, `changelog/gate.md`); the engine's AF-by-first-record, its PAMs at a
contig's last bases (the 5′ last-window item), its contig-end and N-run windows, the Bioconda
enricher (its recipe installs it), PKG1 and the on-demand index check: the engine and fmt lanes'
(CRISPRitz 2.10.0, index format 2, whose builds check what they write), with the CRISPRitz items
of plan 4.4's added block; the mega indels' AF: the rebuild's (plan 2.5); the chr22 ALT-PAM step
of validate-benchmarks, which cannot fail on the split multi-allelic sites (R6B indel
RV2-chr22-gate-inert): the release gates', whose G4 (plan 2.10) checks every allele by brute force
on the rebuild's machine (plan 2.5's phase 5); the 2.8.2 correctness items of list_282.md
(fix/282-indel-alignment, test/282-indel-realdata-check, the SNP rows whose carriers all have an
indel): their branches, forward-ported at the integration (the plan's status log).

*Left: not reproduced, or nothing to change.* RV1-F6 (it names no item of its own: it confirms items
of the list, each named here); the radar chart's guide-length guard (#39, closed as user error);
the "observed haplotypes" wording of FR :558 (kept on purpose); wsearch N2 (the Graphical Reports
caption states the perfect-match stars and the plots drawn before 2.8.1 in general terms: accurate
as it is, as its review found); the tie order of a Personal Risk Card's plots (text r1 F3: no value
changes, and the web's top-1000 plots and the report's agree, as its review measured); C17 (a
banner) and the older list's validate-test hint and scope (no
detail kept in the sources: the list names them alone); t31 P4 (VirtioFS: not reproduced on
38a7c4e); the bulge note counting a merged panel's list as made (FR2: only a failed write at Submit,
which then says why, parts them); the older list's code comments (no behaviour); bio11.2 (a flake
recorded with no detail); the T-31 flock trigger (flock works on Docker Desktop and on Linux; no
ENOTTY host to test on); the conda-build warnings of install P10 (the Bioconda recipe's, not this
repository's); the README's apptainer note, which list_282.md's HPC item would have say so too
(its "consider"): the host environment that reached the container changed nothing, as the image
sets its own `MAMBA_*`, `PATH` and `LANG`, and the checklist runs clean so that its gate is a new
user's run; R6B cr-past-first-lines (build-index-only reads each file whole and refuses it since
2.8.1; a search that builds its own index checks two lines a file, and reading every file whole on
each search, seconds on hg38, is the owner's trade to weigh).

*Left: tests and CI.* Suite P5 (tests that import the web app write `PostProcess/Cache/`, the app's
cache folder relative to the working directory, ignored by git: where the app keeps its cache is a
choice of its layout); R2-F4, R2-F5, wset N3, RV1-F4, N6, text F2, the value tests of plan 4.5
(resultIntegrator's columns, add_risk_score's fields, the CB sentinels, resolve_phased's
manifest), the recovery route of a missing INDELS sub-index, the audit's test additions (a 5′-PAM
synthetic case with brute-force truth, a published test archive) and LP RV4-1 (a registry that
says sites-only by its aggregation alone, as the real mega's does): coverage of code that works,
where a new test must fail before its fix and none has a failing behaviour to fix. The audit's
`complete-test --index-name` and a per-kind example guide (plan 4.5, from the audit's roadmap):
features of the smoke, not a fix of it, which since 2.8.1 (I7) tests the index it ranks first in
that index's own geometry, and takes a run with no site for a pass only on another nuclease's
index. R6B cigate F3's retry loop (each of up to three tries runs every benchmark, against the
job's 120 minutes): the job ends red, never with a false pass, and its logs now go up when the
limit cancels it (above); whether a try that wrote results is tried again is the workflow's
choice. FFR D3R3-6 (the fragment test skips fenced code, which merge_changelog.py reads): no
fragment holds any, so the two read every fragment alike.

*Left: docs and wording.* RV1-F3, N4, F5, RV4-2, the 2.7.10 entry, R6B's three changelog wordings
and R6B emptyidx RV1-3, RV2-4 and R3-6 (the "(or `NNNNNNNN_...`)" of the F3 entry's What to do,
which no job of an earlier release can show): each is the wording of an entry of the released
CHANGELOG (2.8.1's, the CFD entry 2.7.10's), whose correction this track leaves to the release
manager.

*Left: public actions.* t31 P9 (host the 1000G chr22 VCF on the Hub), the HF Xet report, and the
deprecated pinellolab/CRISPRme site's installers (another repository): each needs Luca's go.

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
  it; the advice to run without `--index-path` no longer says to drop `--prune-pam` too. The
  warning of an all-`N` index built before CRISPRitz 2.9.1, which holds no target whose flank has
  an `N`, said to "search without `--prune-pam` and `--index-path`, which builds a dedicated
  index": such a search now prunes by itself from that very index and builds none, so the advice
  led back to the same warning. It names the two routes that serve the search whole: the index
  built again, or a dedicated index of the PAM beside it, which a search takes first. And the NOTE
  of a search that prunes by itself offers a published index "if it is published", as the refusals
  do — none is published for a genome of one's own, or for `TTTV`, `NGA` and most other PAMs.
- **An index is named for its guide length and the side of its PAM (E2).** `<PAM>_<N>_<genome>`
  says neither, so a 3′ and a 5′ build of one motif, or builds for two guide lengths (Cas12a
  `TTTV` at 20 and 23 nt, as the rebuild adds), took one name, and the second build stopped on
  the first ("... is already there ... The index for that PAM file has the same name: build it in
  another working directory (--path)"), as a search did. `build-index-only` and a search now
  build an index under a name that carries both, in the order of the target on the strand:
  `20-NRR_3_hg38` for a 20-nt guide and a 3′ PAM, `TTTV-23_3_hg38` for a 5′ PAM and a 23-nt
  guide, the variant index (`20-NRR_3_hg38+hg38_HGDP1kGP`) and its INDELS index
  (`<name>_INDELS`) too. The engine writes the older name, in a build folder (a search's in the
  index's build lock), and the index is renamed into place: an index of the older name in place
  is never touched (in place, the engine removed a folder of the name it writes before writing
  it). Every reader takes both forms (the search script, `crisprme.py`, the web form, `download`
  and its router, which reads a PAM-less name's flank in either form, `publish-index`'s
  manifest, complete-test, the index coverage, FASTA-header and N-flank checks): an index of the
  new form only for a search of the side its name says, which records the side as a `.pam_build`
  does, ahead of one (such an index needs none, at any flank width), and of one index under both
  names, the newer one. An index of the older name (every one published so far, or built by an
  earlier release) is used as before. `build-index-only` keeps one in place that is the index it
  builds (its `.bin` files built for the PAM file's guide length, its `.pam_build` recording that
  side, or none: a dedicated index, an all-`N` one for a 5′ file, or one at the published widths,
  3 and 8), and builds beside one that records another, saying so (`.../NNNN_2_hg38 is an index
  for a PAM on the 3' side of the guide (as its .pam_build records): it is left as it is, and
  this index is built beside it, as NNNN-20_2_hg38.`). An unrecorded all-`N` one of another
  width, for a 3′ file, is refused as before, the ERROR now naming the `mv` to each side's name,
  which records it; so is a folder of the new name that its trees or its record contradict, with
  the name of what it is. A search's refusal of an unrecorded index for a 5′ PAM names that `mv`
  too (it said to run `build-index-only` with that file, which still records it). An index is
  published under its folder's name, so one built by this release under the new one. An earlier
  release does not read the new form, and builds its own. The documented builds and publishes name
  the folder a build writes: PRECOMPUTED_INDEXES' dict-less flow, which the rebuild follows, said
  that its build wrote `NRR_3_hg38+hg38_HGDP1kGP`, and its `publish-index` of that folder stopped
  on `Index directory ... does not exist`, and the release checklist's build said
  `NNN_3_hg38+hg38_HGDP1kGP`; the panel recipes (`seq_script/merge_panels/`) print the index their
  build wrote, of either form, as the CRISPRme they run names it, and the mega INDELS recipe's
  stop on an unexpected name says which to set (`20-NRG_3` from 2.9.0 on). On chr22 with HGDP +
  1000 Genomes, the documented build wrote `20-NRR_3_hg38`, `20-NRR_3_hg38+hg38_HGDP1kGP` and its
  `_INDELS`, with nothing on stderr; the documented `publish-index` of the variant one went past
  its folder to the token it asks for (none was given), where the older name stops on the missing
  folder; the HGDP + 1000 Genomes recipe's lines after its build printed that folder; and the
  page's search under `--index-path` took both indexes. Covered by
  `PostProcess/test_build_index_only_atomic.py` (both sides and two guide lengths built in one
  library, the variant and INDELS indexes too, an index of the older name kept or built beside,
  each refusal and its `mv`), `test_index_resolution.py` (a search builds its own under the new
  name beside an index of the older name for the other side or another guide length, which it
  stopped on), `test_defaults_and_index_flex.py` (every reader takes both forms alike; every
  documented build and publish names the folder of its build, and each recipe prints the one its
  build wrote) and, through
  the real engine, `check_index_flex.py`: the 3′ and 5′ `NNNN` indexes, and the 20- and 23-nt
  `NGG` ones with their variant and INDELS indexes, built in one library, each searched with the
  results of that index alone, byte for byte; an index under the older name gives the same
  results; a 5′ build beside a 3′ index of the older name, each side's search taking its own; one
  named for a 5′ PAM is not taken for a 3′ search; and a PAM-less index named for its side with no
  `.pam_build` serves a search of that side (a 3′ one of width 4 an `NGG` search, pruned; a 5′ one
  a 5′ PAM of its width).
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
  bulges, which scans the genome (and the variant genome, where it is here whole: the record then
  names neither index, whatever is installed, only the INDELS index the indel search reads); and
  `Prune_pam auto` for a search that pruned by itself, whichever decided it — complete-search,
  which tells the search script that it took the flag (it records `Prune_pam`, `auto` or `given`,
  in the parameters too), or the script, for a web search. The report's
  summary shows them, and its zip bundles the record as `indexes_searched.txt`; a run of an
  earlier release has none.
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
  search with no bulges of `NGCG`, which no installed index holds, built only the index of its
  indels, `20-NGCG_1_hg38+hg38_HGDP1kGP_INDELS`, in 50 s with nothing on stderr, where 2.8.1 built
  `NGCG_1_hg38` and `NGCG_1_hg38+hg38_HGDP1kGP` as well (54 s); its record of the indexes searched
  says `Reference_index none` and `Variant_index:hg38_HGDP1kGP none`, each with why, the second
  naming that index of its indels. The record says `none` for the variant index of such a search
  wherever its enriched genome is whole, an installed index of the search's PAM beside it or not:
  it named that index, with its bulges, as the index searched, and the report's summary with it,
  though only the genome was scanned (every local or Settings build keeps its enriched genome, so
  a web search with no bulges, "Max edits" at 0 bulges, is the common case).
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
- **A PAM file whose first line is not the motif and the length, one space apart, is refused as
  INPUT_FORMATS §4 says.** complete-search checked only that the last field of the first line is a
  number: a line of the length alone made `Results/<output>` and failed the input check later, one
  with a third field (`N...NGG 3 4`) searched with that field for the PAM's length, and a blank
  first line, with the motif on the next, was refused as a file of two motif lines (the focused
  review of f8dbe73). Each gets the one-line refusal of a first line that is not `<motif>
  <length>`, before anything is made.
- **CRISPR-Bulge's scoring leaves nothing in the temporary folder.** TensorFlow's autograph writes
  the code it generates to a module in `$TMPDIR` and removes it at exit, but importing it left
  `__pycache__/__autograph_generated_file*.pyc` there, one a run: in the data folder itself under
  `CRISPRME_TMPDIR=/DATA` (the focused review of f8dbe73). The scorer's worker writes no bytecode
  now (`PYTHONDONTWRITEBYTECODE`, as a user's own setting is kept): a read-only container's search
  with `-e CRISPRME_TMPDIR=/DATA` leaves `/DATA` as it found it, its results those of 2.8.1.
- **A Settings Direct URL download is laid out as an upload is.** A `.fasta` file goes in as a
  `.fa` file, and the hidden files of an archive are left out: a `.fasta` file that needed
  nothing else stayed so, which no index or search reads, and the `._chr1.fa` that macOS `tar`
  adds beside `chr1.fa` went in, with which a search stops (plan 4.8, R6B
  direct-url-fasta-and-hidden-files). `chr1.fa` and `chr1.fasta` in one archive are two of one
  name, refused as an upload refuses them. docs/SETTINGS_DATA_MANAGER.md said so as a limit.
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
  build of the dataset on another is refused, naming its genome. The form reads that record for
  a dataset it offers alone: a request can send any value. Covered by
  `PostProcess/test_web_bulge_notes.py` and `test_web_shell_safety.py` (in the image).
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
  `mm+bulges`), so a search with neither flag, and every web search, merges as before. On chr22
  with HGDP + 1000 Genomes (three guides, NGG, 4 mismatches + 1 + 1 bulge, per-sample), a search
  with neither flag wrote the same 2,023 rows, byte for byte, as this branch before the change,
  and one with `--sorting-criteria-scoring bulges` 82 of them otherwise; both with nothing on
  stderr.
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
  tried, and to set `CRISPRME_TMPDIR` to a writable one or give the container a `/tmp`. A
  folder set with `CRISPRME_TMPDIR` or `TMPDIR` is taken, the working folder included (`docker
  run --read-only -w /DATA -e CRISPRME_TMPDIR=/DATA` runs to its report, its rows those of
  2.8.1); the working folder is not taken when nothing names it, as the scorer and the report
  do not run in it.
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
- **A FASTA whose last line has no line end is refused, and `genome_layout.py --split` ends it
  with one.** The index build drops a last line with no LF, with no error: a sequence written on
  one line lost all its off-targets, and a 60-wrapped chr1 whose last line held a guide's
  on-target gave no row with 2.8.1, and a search with bulges failed (no `.bin`); `--split`,
  which every genome download and Settings upload runs, deleted lone CRs past the first lines,
  which left one long line with no LF (LP unterminated-last-line-dropped). The check a search
  runs first, complete-search's input check and build-index-only refuse such a file (one byte
  a file), saying to run `--split`, which writes a CR between bases as a line end and ends
  every file it writes with LF; after it, the same search found the on-target, with and
  without bulges.
- **The advice for a dataset with no sample list names a VCF file that cannot be read.** A
  Settings build writes the list from the header of every VCF file of the dataset, and writes
  none when one has no readable header (a file cut short, or no VCF); the note and Submit said
  all the same to build the index in Settings, "which writes the list" (M6). They name the file,
  and say to add the dataset again with it whole, or to add the list.
- **A variant search under `--index-path` whose library holds no reference index it can use
  stops before it enriches the genome.** The search script stops on such a library, but only
  after its first step, the enrichment of the genome with the dataset's variants, hours on
  hg38 (IDX-17 e): complete-search stops first, in one line, when that step would run (the
  dataset's VCFs here, no enriched genome, no variant index that serves the search). A
  reference search, which reaches the script's refusal at once and says what the library holds,
  is left to it, and so is a library where a run is building an index, which the script waits
  for. It names the index to download as the script's refusal now does (below).
- **A search refused for want of a reference index under `--index-path` names one to download
  that is published.** The refusal said to download `<PAM>_<N>_<genome>` with N the search's
  bulges + 1: `NRR_2_hg38` for a search of 1 DNA and 1 RNA bulge, the default, where every index
  published is of N = 3, so the download it advised stopped on `No index 'NRR_2_hg38' found`. It
  names the one of N = 3 for up to 2 bulges (`NRR_3_hg38`), "if it is published", as the refusal
  for want of a variant index does, and its other route, to run without `--index-path`, names
  what such a search does since 2.9.0: it takes an index of the working directory's library,
  pruned to the PAM where only such a one serves it, and builds one from the genome where none
  does (it said that it builds one). On chr22 with HGDP + 1000 Genomes, both refusals of a search of
  1 + 1 bulge under an empty `--index-path` name `NRR_3_hg38`, which the Hub holds, and the Hub
  holds no `NRR_2_hg38` (`download` stops on it), nor an `NGCG_3_hg38`. Covered by
  `PostProcess/test_index_resolution.py` and `test_bulge_cap.py`.
- **A variant index of `--index-path` whose dictionaries are not in the working directory is
  named.** complete-search said "No prebuilt variant index found" for a dataset whose variant
  index the library held (IDX-17 c): it says that the index is there, that a search reads the
  dataset's dictionaries, sample list and VCF folder from the working directory, and to run the
  search from the folder that holds them.
- **`crisprme.py` with no argument makes no folder.** It prints the help, and made the seven
  folders of a data folder (`Annotations/`, `Dictionaries/`, `Genomes/`, `PAMs/`, `Results/`,
  `VCFs/`, `samplesIDs/`) wherever one ran it to read it (C14); README §1.1.2, which said so,
  says that it makes none. Covered by `PostProcess/test_complete_search_flags.py`.
- **The web form names the VCF file a search cannot read in a folder of several.** A search of a
  dataset with no variant index for the PAM prepares the variant genome from its VCF folder, which
  takes each file's chromosome from a `.chrN.` part of its name, and the note under the bulges and
  Submit refused a folder with any file named otherwise, saying that "the dataset's VCF is a single
  multi-chromosome file", also of `walk.chr1.vcf.gz` and `walk.chr2.vcf.gz` beside `merged.vcf.gz`
  (RV1-WEB2-5). They say that the folder holds `merged.vcf.gz`, whose name says no chromosome; a
  dataset of one such file is refused as before. Covered by `PostProcess/test_web_bulge_notes.py`
  (in the image).
- **The docs say what a run writes.** The web user guide defined an ALT row (`REF/ALT_origin`) as
  one "created by a variant allele", which the results page's caption contradicts: it is a variant
  alignment, of a site a variant creates (`Not_found_in_REF` `y`) or of a reference site it
  changes (RV1-WEB2-3). METHODS put targets-scoring's CFD-only banner "in its log", which it
  prints on stdout, where a search's `log_verbose.txt` holds it (D3R3-3). The README and
  complete-search's help said that "the CRISPR-Bulge columns read -1" in a CFD-only run, whose
  CRISPR-Bulge risk score, ALT less REF, reads 0 (D3R3-2). And no text said that a line of the
  SNP+indel companion (`indel_snp_cooc.tsv`) lists the first 50 of its carriers by name, which
  its count takes all of: the README and METHODS say so (LP RV4-cooc-carrier-cap). Covered by
  `PostProcess/test_defaults_and_index_flex.py`.
- **The Graphical Reports selector offers the totals of mismatches and bulges the search has.**
  It offered 0 to the mismatches + `Max_bulges`, an index dimension (max(DNA, RNA) + 1 as the
  web writes it, max(DNA, RNA) as complete-search does), and opened on the largest, where the
  search draws the population plots of 0 to its mismatches + DNA + RNA bulges (WEB-4): a web
  search of 4 mismatches and no bulge opened on 5, which has no plot ("No result found for this
  combination of mismatches and bulges"), and a web search of 4 + 2 + 2 could not select its
  sites of 8, nor a complete-search of 4 + 1 + 1 its sites of 6. It offers 0 to mismatches +
  DNA + RNA, up to the cap on them together (Max edits, `--max-total-edits`), and opens on that; a
  job written before the DNA and RNA lines of its parameters is shown as before. On chr22 with
  HGDP + 1000 Genomes, a complete-search of 4 + 1 + 1 drew the plots of totals 4 to 6, and its
  selector offered 0 to 5 and opened on 5: it offers 0 to 6 and opens on 6. Covered by
  `PostProcess/test_web_results_links.py` (in the image).
- **The report of a search of two or more variant datasets names them all.** A search of a
  `--vcf` list of two names its results file after the last dataset alone, and its report
  named that one in its header, inputs card and run manifest, with that dataset's registry
  panel and counts (RV1-F1). The report reads the datasets from the copy of the list the
  search keeps in its result folder: it names each (`1000G, hgdpx`), credits no one registry
  with the run's counts (its panel is the samples the search ran with, each under its
  dataset), and says that each variant's frequency is that of its own dataset's panel, as
  each dataset is searched with its own. A search of one dataset, or a result folder of an
  earlier release, is reported as before: `generate-report` of a chr22 search of HGDP + 1000
  Genomes gives the same labels and registry counts with its list of one dataset, and names both
  datasets, crediting no registry, with a list of two. Covered by
  `PostProcess/test_generate_report.py`.
- **The report of a search with no variant dataset fits a personal haplotype too.** It said
  `none (reference genome only)` and "Sites in the reference genome", also in the report of
  each haplotype of assembly-search, a search of a personal genome (RV1-F2): it says `none
  (no variant dataset searched)`, and "Sites in the genome searched (hg38)", with the name the
  search gives the genome: `WALK01_paternal_575e5055` in the paternal report of an assembly-search
  of the walkthrough's bundle. Covered by `PostProcess/test_generate_report.py`.
- **The report of a search that computed no score says so in place of its score tables.** A
  Cas12a search, or one of any PAM on the 5′ side, computes neither CFD nor CRISPR-Bulge, and
  its report printed the risk card's CFD bands and Section 4's CFD tiers, each 0, and the
  CRISPR-Bulge and CFD sentences of the recommended panel: read as "no site at risk", where
  the panel was chosen by mismatches and bulges alone (cli5 D6). The card and Section 4 say
  that neither score is computed for the nuclease and that the panel is chosen by mismatches
  and bulges, and the panel's note says how; assembly-search's combined report too. A scored
  search's report is unchanged. A Cas12a (`TTTV`) search of chr22 with no bulges reports so.
  Covered by `PostProcess/test_generate_report.py`.
- **`download --what index` brings again the dictionaries an installed variant index has
  lost.** A search of a variant index whose dictionaries were gone
  (`Dictionaries/registry_<dataset>` or `dictionaries_<dataset>`, and `log_indels_<dataset>`) said
  to download the index again,
  and complete-search's check said to download it, but the download took the index for
  installed, as its folder, its INDELS index, the reference index and the genotype store were
  there, and did nothing (FR :163). It counts the dictionaries a search reads as part of the
  install, and brings them with the rest of the index's archive. The search's advice names the
  command, and says that in the web app, whose Settings page downloads no index that is
  installed, the index is removed there first. With the Hub offline, a download of
  `NRR_3_hg38+hg38_HGDP1kGP` over an install of it without the dictionaries said "already
  installed" before, and goes to fetch it now. Covered by `PostProcess/test_crisprme_hf.py` and
  `test_index_resolution.py`.
- **assembly-search searches a haplotype again when its FASTA files changed since its
  result.** A run of the same `--output` reused a finished haplotype search whose parameters
  and genome folder were the run's, whatever the folder held: after a bundle was imported
  again under the same path with other FASTA files, the result of the old sequence was
  reconciled (R4-6). Each haplotype's result records the genome it was searched from (the
  folder and the fingerprint of its FASTA files, as an index records it), and a run reuses it
  only while the genome is that one; a result of CRISPRme 2.8.1 or earlier, which records
  none, is searched again. An assembly-search of the walkthrough's bundle, run again after a base
  of the maternal on-target was changed, reused the paternal result and searched the maternal
  haplotype again (one perfect match, not two), where the code before reused both. Covered by
  `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search keeps no INDELS index it built while build-index-only put the SNP codes on the
  fake-indel genome.** build-index-only marks that genome unfinished only just before the
  first FASTA it changes, so a search past its check of the mark before then built its INDELS
  index from the FASTAs not coded yet, and every later run kept it, blind to a SNP beside an
  indel on them (R2-F1). An INDELS index during whose build a FASTA of the genome changed, or
  after which the mark is there, is not put in place: the search stops, saying to let that
  run finish and to run the search again, which builds it from the coded genome. Covered by
  `PostProcess/test_index_resolution.py`.
- **The top-1000 plots of a run with no score name no perfect match they do not star.** Where
  the score was not computed (a nuclease of no SpCas9 family, or a run without the scorer),
  its column holds -1, and the stars of the perfect matches were drawn there, below the axes,
  while the legend listed them (wsearch N3). A perfect match is starred, and in the legend,
  only at a score of its own; the plot of the count of edits stars them as before. The
  results page's plots and the Personal Risk Cards' alike: the CFD plot of the Cas12a search above,
  whose one site is its on-target, names no perfect match. Covered by
  `PostProcess/test_top_plots.py`.
- **The web's per-sample option, its word on a sites-only dataset and the report's advice say
  what a web user can do.** The search form disabled Per-sample for a genotyped dataset whose
  index is not built yet (listed as "genotyped, samples not indexed") saying that it "needs a
  genotyped panel"; the note under its bulge options, and Submit's refusal, said of a sites-only
  dataset that only an allele-frequency registry reports its variants, "as the shipped hg38_mega
  index has (docs/DESIGN_mega_index.md)", a file of the repository a web user cannot open; and
  the report of a population-level search said to re-run it with `--per-sample` alone, a flag a
  web search cannot pass (FV :342), as it named commands alone to get the tiers of a dataset
  without them, and to bundle what its slim zip leaves out (the variant-created list and subset,
  the complete raw integrated results, the windows near the top-ranked sites, an
  observed-haplotypes companion over the size gate). The option
  says that the dataset's samples are not indexed yet, and to build its index in Settings; the
  note and the refusal name the hg38_mega index published on HuggingFace, which Settings
  downloads (under Download a prebuilt index); the report names the web form's Per-sample
  analysis beside the flag in both sentences that tell the reader to re-run the search per-sample
  (the Search mode row, and the caption of a population-level search's variant combinations),
  Settings' Build an index locally and Download a prebuilt index
  beside the commands that get the tiers, and the form's Include the full report beside
  `--report-full` in every sentence that names it (the website's example report with them).
  Covered by
  `PostProcess/test_web_shell_safety.py`, `test_web_bulge_notes.py` and
  `test/web/test_settings_vcf_datasets.py` (in the image), and `test_generate_report.py`.
- **build-index-only refuses an INDELS index in place that lacks the sub-indexes it builds.** It
  said "Indels index already present", then "Variant index ready" at exit 0, of any INDELS index
  under its name: one copied in part, or one where a search before 2.8.1 had put the sub-indexes
  of a PAM-less index it served (`NNN_2_fakechr1` in `NRG_2_<genome>+<dataset>_INDELS`), each of
  which every search with bulges refuses (the index audit's leftover INDELS folder, plan 4.5).
  It now looks for the sub-index of each chromosome whose fake-indel genome holds indels (none is
  due for a contig with an underscore in its name, or a genome of headers alone), and stops with
  an ERROR that names the ones missing, what the folder holds, and how to build it again. With
  the real engine on a synthetic dataset, a run over that folder of `NNN_2_fakechr1` and
  `NNN_2_fakechr2`, or with `NRG_2_fakechr2` gone, ended 0 before and ends 1 now, and a search
  with bulges refuses the same folder; once it is removed, the build makes it whole. On chr22 of
  HGDP1kGP, the published `NRR_3_hg38+hg38_HGDP1kGP` with its `NRR_3_fakechr22` renamed
  `NNN_3_fakechr22`, or gone, is refused the same way, as a search with bulges refuses it, and is
  present again once named back. Covered by `PostProcess/test_build_index_only_atomic.py`.
- **CI keeps the logs of a job its time limit cancels, and the checks of what the workflows run
  read no comment and see every skip.** validate-benchmarks uploaded its logs on `failure()`
  alone, so a job that up to three tries of every benchmark took past its 120 minutes, which
  cancels it, kept no log of what failed (R6B cigate F3): it uploads them on `cancelled()` too.
  `test_ci_workflows`' checks that every test module is named in a workflow, and that each
  function-based one runs as a script, read the workflows' comments too, where unit-tests.yml's
  paths comment names a dozen modules: with the line that runs `test_dash_output_props` gone, the
  first still passed (R6B cigate F2). Its check that the light job runs the flag-gate scans looked
  at a skip of their classes alone, and passed with one of their tests skipped (R6B cigate F1).
  Both fail on those now. And `test_history_page_assembly` reads the page with its file closed
  (it printed six ResourceWarnings: FFR RV1-WEB2-6). Covered by `PostProcess/test_ci_workflows.py`.
- **The docs say what the self-test's opt-out saves, where an empty `_INDELS` comes from, and how
  common two indels in cis are.** The README, the benchmark README and the registry called
  `CRISPRME_SKIP_HEAVY=1` an opt-out "for very small machines", which reads as one for a machine
  short of memory: the run's peak is the variant enrichment, which every benchmark needs, so it
  saves time on a slow machine, not memory, as they now say (R6B cigate F4). RELEASE_CHECKLIST's
  archive gate (T-40) said that a present-but-empty `_INDELS` is always a broken build; a dataset
  whose every indel is on a contig with an underscore in its name gets one from a build that ends
  well, and the gate names it beside the broken build, neither to be published (R6B emptyidx
  RV1-2, RV2-2, R3-1). FUTURE_WORK §3 called two or more indels in cis common, "~17–19% of real
  indel off-targets": it says that this is the raw count, mostly repeats, and gives METHODS §8's
  ~1–2% of indel loci after repeat-masking (plan 4.1). Covered by
  `PostProcess/test_ci_workflows.py` and `test_defaults_and_index_flex.py`.
