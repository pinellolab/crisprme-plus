# CRISPRme+ release checklist

A living, run-every-release gate. It has two halves:

1. **Gates** — what must be true before a new `vX.Y.Z` tag.
2. **Regression traps** — concrete things we have gotten wrong before, each paired with the
   check that catches it. **Append to this list whenever a release pass finds a new problem**,
   so the same class of bug can never ship twice.

> Versioning is **stingy** — stay in the `2.x` line (patch bumps) unless the change is a true
> CRISPRme-3 break. Never cut a tag on anything less than all-green CI.

---

## 1. Gates

### A. Release-surface consistency (the invariant)
Every user-facing surface must agree on **version · default index · scorer · licence · analysis-mode · data-layout**. The surfaces:
`README.md` · `METHODS.md` · everything under `docs/` **including the GitHub-Pages website `docs/index.html`** and the `.md` pages it links · installation (`README` Conda/Docker/source, `docs/DOCKER_QUICKSTART.md`, `install_from_source.sh`, `Dockerfile`, `environment.yml`, the Bioconda recipe) · **the Hugging Face dataset card** (`lucapinello/crisprme-data` README).

- [ ] **Version sweep** — `crisprme.py` `version`, the tag/image `vX.Y.Z`, and every CRISPRme version string in README + `docs/` read the SAME `X.Y.Z`. The `Dockerfile` pins no CRISPRme version (its header comment names the release: bump it with the rest, as `scripts/prepare_release.py` leaves it); its `ARG crispritz_ref` is the CRISPRitz pin (v2.9.1 for 2.8.1) and must equal `install_from_source.sh`'s `CRISPRITZ_REF` default, the lower bound of the Bioconda recipe's bounded `crispritz` pin ([`RELEASING.md`](RELEASING.md)) and the CRISPRitz version the docs name. Grep gate: no stale `2.5.x/2.6.x/2.7.x` strings or `v2.6.x` image tags linger outside historical/CHANGELOG/archived-doc context.
- [ ] **Default index** — every surface names `NRG_3_hg38+hg38_HGDP1kGP` (SHAPEIT5 re-phased, 4,091 = 925 HGDP + 3,166 1000 Genomes) as the recommended **variant** default; `NRG_3_hg38` is the reference-only default. No surface proposes a retired/alternative index as *the* default.
- [ ] **No retired index advertised** — grep the whole repo (code help strings, scripts, tests, docs) + the HF card for retired names; they may appear ONLY in a past-tense/retirement note. Retired set: `NRG_3_hg38+hg38_1000G2021` (+ `genotypes_hg38_1000G2021`, superseded by the HGDP1kGP default), `NRG_3_hg38+hg38_1000G2021_HGDP` (hybrid), any `*-dictless*`, `genotypes_hg38_1000G_HGDP`. (`NGG_3_hg38` was retired once and **re-published 2026-10-25**: it is the shared reference index that `download --what index` co-fetches for `NGG_3_hg38+hg38_HGDP1kGP`, so it is live and must stay documented.)
- [ ] **Scorer** — CFD + CRISPR-Bulge everywhere; zero mentions of CRISTA or azimuth (removed).
- [ ] **Licence** — MGBOA 1.0 everywhere (README, METHODS, `LICENSE`, HF card frontmatter `license: other` + `license_name`/`license_link`, report footer). The migration table states commercial use needs a separate licence for **both** CRISPRme+ and stock CRISPRme.
- [ ] **Analysis mode** — the context-aware default (per-sample on a genotyped + phased index such as the default HGDP1kGP one, population-level otherwise) with the `--per-sample` / `--population-level` overrides; zero mentions of the removed `--fast`/`--full`.
- [ ] **Published-index provenance recorded** — every archive on HF is accounted for in
  [`INDEX_PROVENANCE.md`](INDEX_PROVENANCE.md) with its stored flank width (which decides
  whether an older client can read it) and whether the current CRISPRitz builder would
  reproduce its content. A new archive, or a CRISPRitz change touching `mainParallel.cpp`
  or `pamGenerator.cpp`, means re-running that comparison. `SHA256SUMS` present under
  `indexes/` in both dataset repos.
- [ ] **Every published variant archive can find indels — checked in the ARCHIVE** (T-40). For
  each `indexes/*+*.tar.gz` (download it, or stage it before upload), list its members and
  confirm: every `<index>_INDELS/<sub-index>/` holds at least one `.bin`, or the `.no_targets` record of a build that found no target, unless that chromosome's
  `Dictionaries/log_indels_<vcf>/log<chr>.txt.gz` is header-only (a SNP-only chromosome), and the
  `log_indels_<vcf>` logs are not all header-only. A present-but-empty `_INDELS` is always a broken
  build: a dataset with no indels gets no `_INDELS` at all. Also confirm every archive built on a
  panel bundles the `Dictionaries/` and `samplesIDs/` the published NRG archive of that panel
  installs, **by content, member by member** (SHA-256 of every member; an install REPLACES each
  directory, so a mismatch means whichever archive was installed last wins). Names and sizes are
  not enough: equal sizes passed while 46 files differed. A `.gz` member may differ only in the
  gzip header's 4-byte MTIME (bytes 4-7); compare those with that field masked, and treat any
  other difference as a failure. Applies to both dataset repos — the K = 8 archive replaces the
  same `Dictionaries/` folders as the NRG index it sits beside.
- [ ] **HF catalog matches reality** — `crisprme_hf.list_available_downloads` / the Settings dropdowns / the docs list EXACTLY the indexes present on HF, in **both** dataset repos (`list_repo_files` on `lucapinello/crisprme-data` and on `lucapinello/crisprme-data-general-index`; the Settings catalog merges the two). Nothing listed that 404s; nothing on HF that is undocumented; `genotypes_hg38_HGDP1kGP.tar.gz` has the same SHA-256 in both repos.
- [ ] **Every documented `download --what X` example resolves on HF** — for each `download --what {genome,index,vcf,annotations,pams,samples,assembly,all}` command shown in ANY doc (README, DOCKER_QUICKSTART, index.html, DATA_LAYOUT, SETTINGS, web guide, legacy setup), confirm the HF path it fetches actually exists (`list_repo_files`). A documented command that targets unhosted data is a broken step (see trap T-21: `--what vcf --dataset 1000G|HGDP` pointed at a `vcfs/` folder that was never uploaded).

### B. Functional gates (must pass)
- [ ] **Clean-room new-user e2e** on current `main` (overlay into the latest SIF, or a fresh Docker pull): `download --what all` → a real search → report renders with closest-gene columns + CONFIRMED cis + named carriers from HGDP1kGP. The gold-standard gate.
- [ ] **From-scratch RC image test** (stronger than overlay; closes T-20). Build a BAKED RC image from the release branch and run the full new-user flow against it with **no code binds**. On a Docker host: `docker build` the branch → run `download --what all` + `download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP` + pop-level + `--per-sample` 4/1/1 + `generate-report`. **On a restricted/corporate network that blocks `conda.anaconda.org` (the local/Mac-build SSL failure), pass `--build-arg CONDA_CHANNEL_BASE=https://prefix.dev`** (the Dockerfile threads it into both the main and scorer-env conda solves; see T-37). GitHub-Actions CI builds on GitHub runners (open network) so it uses the default channel fine. On an **apptainer-only node (no Docker, e.g. ml003)**: `apptainer build --fakeroot rc.sif` from a def that bootstraps the latest released SIF and, in `%setup`, copies the branch source into `/opt/conda/opt/crisprme` + `/opt/conda/bin/crisprme.py` (mirrors the Dockerfile `COPY .`); invoke via `apptainer exec --env PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin` (this SIF has no runscript — see T-33). This exercises the BAKED `bin/` wrapper + baked PostProcess/calibration, not a runtime overlay.
- [ ] **`--per-sample` on Docker Desktop for Mac** (virtiofs bind mount) — not only Linux. The annotation flock race (T-31) is filesystem-specific; Linux/cluster never reproduces it.
- [ ] **Report threshold-message consistency** (T-32) — the shipped `panel_top100.tsv` IS the recommended default panel `(CB≥1.9e-4, 95%-recall) ∪ (mm+b≤2) ∪ (CFD≥0.2@0-bulge) ∪ (CFD≥0.5@1-bulge) ∪ clean-bulge ∪ must-includes` (a superset of the previous panel); the Section-4 CRISPR-Bulge tier table carries the sensitive `≥1.9e-4 (95% recall)` row + `crispr_bulge_gt_0.tsv` (stable file key); and the risk card, the "how to read this" modal, the Section-4 note, the score glossary, and METHODS §3.5 tell ONE threshold story. The new REF/aligned columns (`CFD_ref`, `CRISPR-Bulge_ref`, `Aligned_spacer+PAM`) appear in the in-report top-1000 table AND every download TSV.
- [ ] **Full CI green on the release PR** — `build`, `unit-tests`, `validate-benchmarks`, `web-e2e`, `Build multi-arch image`. ⚠️ A PR touching only `pages/**` or docs does **not** trigger `validate-benchmarks` (path filters — see trap T-09); the version bump touches `crisprme.py`, which re-triggers the full suite, so confirm they actually ran.
- [ ] **Web walkthrough (Playwright)** — home search submit → `/load`; Settings download index / add annotation / download genome / build index / delete genome; assembly-search submit. Screenshots archived.
- [ ] **CLI test suite, with nothing left skipped** (T-41) — the local gate below, on the release image, amd64 and arm64.

#### The local gate: nothing stays skipped

Run it on the release image with no overlay, on amd64 Linux and again on the Mac (arm64): CI builds arm64 and tests nothing on it, and `validate-benchmarks` overlays the published image, so it never sees the release code with the release's CRISPRitz. So step 4b runs here, on the release's own engine, the three engine checks CI runs only on that published image: `test/index_flex/check_index_flex.py`, `test/indel_snp_oracle/run_oracle.sh` and `scripts/check_alt_pam.py`. It needs Docker, git-lfs, and network to Docker Hub, Hugging Face, EBI, GitHub, MCR, PyPI (step 5) and the Debian package mirror (step 2). 4 CPUs and 16 GB are enough; it took 36 min on the Mac and 55 min on a busy amd64 node, step 4 being 24 and 34 of them. Step 4b adds about 16 min on the Mac (index_flex 14, the oracle under 2) and 30 to 40 min on a busy amd64 node (index_flex 28 to 37, the oracle 3); its index_flex and oracle need nothing from step 4 and can run beside it. Step 4c, the independent real-data indel check, needs a chr22 HGDP + 1000 Genomes install and the genome-wide per-sample runs of the README guide and TRAC; the two chr22 searches are about 5 min each on 8 threads, the four chr22 checks a few minutes in all with `--jobs 24` (soundness 10 to 15 s; completeness of the whole chromosome 100 to 180 s depending on the node's load, about 40 CPU-minutes), and genome-wide soundness under 2 min per guide with `--jobs 24` once the VCFs have a `.tbi` (each unindexed VCF is first read once into an in-memory index).

**Three skips are allowed, and only these**, each run elsewhere in the gate: `no TensorFlow here` (`test_crispr_bulge_fidelity`, the CRISPR-Bulge fidelity gate, runs in its own env in step 3), `en_US.UTF-8 is not installed` (`test_index_resolution`'s locale tests, step 2) and `no ps here to compare` (`test_index_resolution`, runs in the `unit-tests` job, whose runner has `ps`).

**Six function-based modules must run as scripts**: `test_determinism`, `test_indel_genotypes`, `test_indel_snp_bruteforce`, `test_indel_snp_cis`, `test_indel_snp_overlay`, `test_samplesid_overlisting`. Under `python -m unittest` each reports `Ran 0 tests` and `OK`, a false pass; as scripts their runners call every test and raise on a failure. Step 1's loop does it. Run each module in a process of its own, as the CI steps do.

```bash
IMG=pinellolab/crisprme:vX.Y.Z      # or an RC: docker build -t crisprme:rc . (T-37 if conda.anaconda.org is blocked)
SRC=$PWD                            # the release checkout (for test/web)
GATE=$PWD/gate; rm -rf "$GATE"; mkdir -p "$GATE/heavy"
docker run --rm "$IMG" crisprme.py scorer-env check   # must print: compute self-test: OK

# 1. Every module, the function-based ones as scripts, the hermetic test/web scripts, the HF network test.
docker run --rm -e CRISPRME_RUN_HF_TESTS=1 -v "$GATE:/gate" -w /opt/conda/opt/crisprme/PostProcess "$IMG" bash -c '
  fail=0
  for f in test_*.py; do m=${f%.py}
    python -m unittest "$m" -v > "/gate/$m.log" 2>&1 || { fail=1; echo "FAIL $m"; }
    if grep -q "^Ran 0 tests" "/gate/$m.log"; then python "$f" > "/gate/$m.script.log" 2>&1 || { fail=1; echo "FAIL $m (script)"; }; fi
  done
  for t in ../test/web/test_*.py; do case "$t" in *_e2e.py) continue;; esac
    python -u "$t" > "/gate/web_$(basename "$t" .py).log" 2>&1 || { fail=1; echo "FAIL web $t"; }
  done
  exit $fail' || echo GATE_FAILED
# Must print nothing: any other skip is a test that did not run.
grep -hE "skipped ['\"]" "$GATE"/test_*.log | grep -vE "no TensorFlow here|en_US.UTF-8 is not installed|no ps here to compare"

# 2. The en_US arm of the locale tests.
docker run --rm -w /opt/conda/opt/crisprme/PostProcess "$IMG" bash -c '(apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq locales) >/dev/null 2>&1; localedef -i en_US -f UTF-8 en_US.UTF-8; python -m unittest test_index_resolution.ChoiceIsLocaleIndependent -v' 2>&1 | tee "$GATE/en_us.log"
grep -q "^OK$" "$GATE/en_us.log" || echo GATE_FAILED   # Ran 4 tests, OK, no skip

# 3. CRISPR-Bulge fidelity, at the commit the image provisions (a 524 MB LFS download).
PIN=$(docker run --rm -w /opt/conda/opt/crisprme/PostProcess "$IMG" python -c 'import scorer_env; print(scorer_env.CBULGE_PIN)')
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/OrensteinLab/CRISPR-Bulge.git cb
git -C cb checkout "$PIN"
git -C cb lfs install --local
git -C cb lfs pull --include files/datasets.zip
unzip -q cb/files/datasets.zip 'datasets/Refined_TrueOT.csv' -d cb/files
docker run --rm -v "$PWD/cb:/cb" -e CBULGE_REPO=/cb -w /opt/conda/opt/crisprme/PostProcess "$IMG" micromamba run -n cbulge python -m unittest test_crispr_bulge_fidelity -v 2>&1 | tee "$GATE/fidelity.log"
grep -q "^OK$" "$GATE/fidelity.log" || echo GATE_FAILED   # Ran 2 tests, OK, no skip; [fidelity] ... MATCH

# 4. Every benchmark, with no CRISPRME_SKIP_HEAVY.
docker run --rm -e CRISPRME_ALLOW_ONDEMAND_BUILD=1 -v "$GATE/heavy:/DATA" -w /DATA "$IMG" crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G --thread 8
docker run --rm -v "$GATE/heavy:/DATA" -w /DATA "$IMG" crisprme.py validate-test --chrom chr22 2>&1 | tee "$GATE/validate.out"
grep -q "Validation summary: 6 passed, 0 failed, 0 skipped" "$GATE/validate.out" || echo GATE_FAILED

# 4b. The engine checks, on the release's CRISPRitz: which index serves which search, the SNP+indel
#     oracle, and the ALT-PAM and no-N gate on step 4's tables. index_flex runs in the container's own
#     /tmp, as CI runs it (its read-only cells chown files to uid 65534). On an engine without
#     --prune-pam it prints INDEX FLEX PASS too, taking its refusal branch: so grep for the K=8 cells.
docker run --rm "$IMG" python3 /opt/conda/opt/crisprme/test/index_flex/check_index_flex.py /tmp/index_flex 2>&1 | tee "$GATE/flex.log"
grep -q "^INDEX FLEX PASS" "$GATE/flex.log" && grep -q "^PASS  K=8 with --prune-pam: same results as a dedicated index" "$GATE/flex.log" && grep -q "^PASS  complete-test on K=8: same results as on NRG" "$GATE/flex.log" || echo GATE_FAILED
docker run --rm -v "$GATE:/gate" "$IMG" bash /opt/conda/opt/crisprme/test/indel_snp_oracle/run_oracle.sh /gate/oracle 2>&1 | tee "$GATE/oracle.log"
grep -q "^ORACLE PASS" "$GATE/oracle.log" || echo GATE_FAILED
python3 "$SRC/scripts/check_alt_pam.py" "$GATE"/heavy/Results/crisprme-test-out_*/*integrated_results.tsv || echo GATE_FAILED

# 4c. The independent real-data indel check (scripts/check_indel_rows.py), on a REAL panel:
#     validate-test compares only the engine's sites on the reference and the SNV-only chr22
#     1000 Genomes genome, so no indel row's coordinate, REF text, variants or carriers are
#     ever checked against an independent truth (test/benchmark/README.md). The checker
#     writes each sample's haplotypes from the VCF itself (every allele of a slot) and
#     aligns with its own exhaustive aligner; it imports and copies no CRISPRme code and
#     needs only Python's standard library (the repository's CFD tables are read as data).
#     $D22 is a chr22-only HGDP + 1000 Genomes install (Genomes/hg38, VCFs/hg38_HGDP1kGP,
#     genome_library, Dictionaries), $VCF its chr22 VCF (a .tbi beside it makes reads seeks;
#     without one the VCF is read once into an in-memory index, nothing is written).
#     Run the README guide, NRG, 4 mismatches + 1 DNA + 1 RNA bulge, in BOTH modes first:
for mode in per-sample population-level; do
  docker run --rm -v "$D22:/DATA" -w /DATA "$IMG" crisprme.py complete-search --genome Genomes/hg38 \
    --pam PAMs/20bp-NRG-SpCas9.txt --guide guide.txt --vcf list_vcf.txt --samplesID list_samplesID.txt \
    --annotation Annotations/dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz \
    --gene_annotation Annotations/gencode.protein_coding.bed.gz --index-path genome_library \
    --mm 4 --bDNA 1 --bRNA 1 --thread 8 --output "gate_${mode}" --"$mode"
  #   soundness (every indel/MNV row real) and completeness (no real indel off-target
  #   missed: per sample on the per-sample run, where every carrier target must be
  #   reported for its carrier by rows of its site and strand no worse than its best
  #   alignment by the merge's sort keys (--sorting-criteria, then the CFD at three
  #   decimals and the criteria again) -- the fold group of the row that lists it, or the
  #   rows a fold puts it under; per merge cluster on the population-level run, which keeps
  #   one row per cluster -- its row must have no more edits than each target and at least
  #   its CFD, at the three decimals the merge sorts on; the checker reads the mode from
  #   Results/*/.alt_mode, the cap from .Params.txt, the criteria from .command_line.txt
  #   and the panel from .sampleID.txt, or give --samples)
  for check in soundness completeness; do
    python3 "$SRC/scripts/check_indel_rows.py" --results "$D22/Results/gate_${mode}" \
      --reference "$D22/Genomes/hg38" --vcf "$VCF" --pam "$D22/PAMs/20bp-NRG-SpCas9.txt" \
      --guide "$D22/guide.txt" --mm 4 --bDNA 1 --bRNA 1 --max-total-edits 6 --chrom chr22 \
      --mode "$check" --jobs 24 --out "$GATE/indelrows_${mode}_${check}.tsv" || echo GATE_FAILED
  done
done
#     Genome-wide SOUNDNESS of the README guide and of TRAC: $GW is a genome-wide HGDP + 1000
#     Genomes install whose Results/readme and Results/trac are the release's per-sample
#     searches of guide.txt (README) and guide_trac.txt (TRAC), NRG, 4 + 1 + 1. --vcf-dir
#     is the panel's VCF folder; give it a .tbi per VCF (in a folder you can write: link
#     the VCFs there and run `tabix -p vcf`), or let the checker index each one in memory
#     (one decompression of the file, about 28 GB for this panel: the README run took
#     388 s and 8.3 GB so on 24 processes of a node at load 250 to 300, against 89 s and
#     4.5 GB with a .tbi).
for g in readme:guide.txt trac:guide_trac.txt; do
  python3 "$SRC/scripts/check_indel_rows.py" --results "$GW/Results/${g%%:*}" \
    --reference "$GW/Genomes/hg38" --vcf-dir "$GW/VCFs/hg38_HGDP1kGP" \
    --pam "$GW/PAMs/20bp-NRG-SpCas9.txt" --guide "$GW/${g#*:}" --mm 4 --bDNA 1 --bRNA 1 \
    --max-total-edits 6 --mode soundness --jobs 24 --out "$GATE/indelrows_gw_${g%%:*}.tsv" \
    || echo GATE_FAILED
done
#     The gate is the exit status: 0 gating findings. Non-gating lines are documented
#     CRISPRme limits (limit-last-alt-wins, limit-two-indels, limit-putative-carrier,
#     limit-contig-end), slots whose VCF alleles overlap where a target of their record
#     reads, so the VCF gives it no sequence (limit-overlapping-records: 3 records on
#     chr22), windows across N (limit-n-bases), a sites-only dataset's rows, which list no
#     carrier (limit-no-genotypes: none here; to check a search of the default index beside
#     hg38_mega, give both VCFs, --vcf "$VCF" --vcf <mega chr22 VCF>), the note
#     variant-shadowed, and the checker's own coverage gaps (coverage-putative-cap: a
#     carrier resolved PUTATIVE with more than 10 positions near a record where its
#     alleles allow two bases, whose targets touching a call of unknown phase are not
#     enumerated, none on this phased panel; coverage-chain-limit: a population-level
#     cluster longer than the 3 kb either side its chain follows;
#     coverage-unindexed-long-ref: a VCF without a .tbi holding a deletion longer than
#     20 kb, which an unindexed read may miss: give it a .tbi), each listed with its
#     samples and target: read them, and say in the release notes if a class grew. A
#     known upstream defect can be listed with --known (class, chrom, start, strand,
#     variants, samples, reason): it turns non-gating as known-defect, and a listed one
#     that stops showing is reported as known-not-seen. Exit 2 is no verdict: a usage
#     error, or a configuration the checker refuses (two genotyped datasets that share
#     samples; check each dataset's search alone). Genome-wide COMPLETENESS (at most
#     about 52 CPU-hours per guide: 53 to 78 min on 40 processes of a loaded node, 6 GB, a
#     TSV of about 31 MB) is not a release gate; chr22 above stands for it.

# 5. The browser tests. On a Mac: drop --network host and use CRISPRME_WEB_URL=http://host.docker.internal:8090.
docker run -d --rm --name gate-web -p 8090:8080 "$IMG" crisprme.py web-interface
until curl -sf -o /dev/null http://localhost:8090/index; do sleep 3; done
for t in test_web_e2e test_settings_upload_e2e; do
  docker run --rm --network host -v "$SRC/test/web:/t:ro" -e CRISPRME_WEB_URL=http://localhost:8090 mcr.microsoft.com/playwright/python:v1.47.0-jammy bash -c "pip install -q playwright==1.47.0 && python -u /t/$t.py" || echo GATE_FAILED
done
docker rm -f gate-web

# 6. All four workflows green on the release PR, and confirm in each run log that unit-tests and validate-benchmarks actually ran (T-09).
```

Nothing above may print `GATE_FAILED` or `FAIL`. When `benchmarks.json` or a reference changes, also regenerate the references with `test/benchmark/generate_references.py` (see its README): no CI runs it.

**On an HPC node without Docker**, a plain SIF is read-only, so use a writable sandbox. One more skip is allowed there, in step 1: `no group to give genome_library` (an unprivileged user namespace maps no second group; the Docker run above tests it).

```bash
apptainer build --sandbox sb crisprme.sif && mkdir -p "sb$PWD" gate/heavy tmp locale hf
# step 1: save its loop as gate_step1.sh with /gate replaced by $PWD/gate, then
apptainer exec --writable --no-home --bind "$PWD" --pwd /opt/conda/opt/crisprme/PostProcess --env TMPDIR="$PWD/tmp" --env CRISPRME_RUN_HF_TESTS=1 --env HF_HOME="$PWD/hf" sb /usr/local/bin/_entrypoint.sh bash "$PWD/gate_step1.sh"
grep -hE "skipped ['\"]" gate/test_*.log | grep -vE "no TensorFlow here|en_US.UTF-8 is not installed|no ps here to compare|no group to give genome_library"
# step 2 (needs the host's /usr/share/i18n):
apptainer exec --bind /usr/share/i18n:/usr/share/i18n --bind "$PWD" sb /usr/local/bin/_entrypoint.sh localedef -i en_US -f UTF-8 --no-archive "$PWD/locale/en_US.UTF-8"
apptainer exec --env LOCPATH="$PWD/locale" --bind "$PWD" --pwd /opt/conda/opt/crisprme/PostProcess sb /usr/local/bin/_entrypoint.sh python -m unittest test_index_resolution.ChoiceIsLocaleIndependent -v
# step 3 (home bound: micromamba writes ~/.cache/mamba), after the clone of step 3 above, with
PIN=$(apptainer exec --pwd /opt/conda/opt/crisprme/PostProcess sb /usr/local/bin/_entrypoint.sh python -c 'import scorer_env; print(scorer_env.CBULGE_PIN)')
apptainer exec --env CBULGE_REPO="$PWD/cb" --bind "$PWD" --pwd /opt/conda/opt/crisprme/PostProcess sb /usr/local/bin/_entrypoint.sh micromamba run -n cbulge python -m unittest test_crispr_bulge_fidelity -v
# step 4:
apptainer run --bind "$PWD/gate/heavy:/DATA" --pwd /DATA --env CRISPRME_ALLOW_ONDEMAND_BUILD=1 sb crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G --thread 8
apptainer run --bind "$PWD/gate/heavy:/DATA" --pwd /DATA sb crisprme.py validate-test --chrom chr22
# step 4b (the same greps as above; unprivileged, index_flex's read-only cells chmod instead of chown):
apptainer exec --bind "$PWD" --pwd "$PWD" sb /usr/local/bin/_entrypoint.sh python3 /opt/conda/opt/crisprme/test/index_flex/check_index_flex.py "$PWD/gate/index_flex" 2>&1 | tee gate/flex.log
grep -q "^INDEX FLEX PASS" gate/flex.log && grep -q "^PASS  K=8 with --prune-pam: same results as a dedicated index" gate/flex.log && grep -q "^PASS  complete-test on K=8: same results as on NRG" gate/flex.log || echo GATE_FAILED
apptainer exec --bind "$PWD" --pwd "$PWD" sb /usr/local/bin/_entrypoint.sh bash /opt/conda/opt/crisprme/test/indel_snp_oracle/run_oracle.sh "$PWD/gate/oracle" 2>&1 | tee gate/oracle.log
grep -q "^ORACLE PASS" gate/oracle.log || echo GATE_FAILED
apptainer exec --bind "$PWD" sb /usr/local/bin/_entrypoint.sh python3 /opt/conda/opt/crisprme/scripts/check_alt_pam.py "$PWD"/gate/heavy/Results/crisprme-test-out_*/*integrated_results.tsv || echo GATE_FAILED
# step 4c: the same commands as above, with `apptainer run --bind ... sb crisprme.py complete-search`
# for the two searches. check_indel_rows.py itself needs only the host's python3 (stdlib only);
# --jobs N splits it over N processes.
```

### C. Release mechanics (at/after tag)
- [ ] Bump `crisprme.py` version; move `CHANGELOG [Unreleased]` → `[X.Y.Z]` dated + fresh `[Unreleased]` + fix link footer.
- [ ] Tag → `docker-multiarch` CI publishes `pinellolab/crisprme:vX.Y.Z` + `:latest`; verify `crisprme.py --version` == `X.Y.Z` and that `crispritz.py`'s `VERSION` inside the image equals the Dockerfile's `crispritz_ref` (v2.9.1 for 2.8.1).
- [ ] Bioconda: refresh the recipe/Gmail draft + `sha256` of the new release tarball; README's bioconda-status line matches reality.
- [ ] HF dataset card re-checked against the shipped catalog (part of surface A).
- [ ] Ping Ann (@anncir1) to review any web/assembly-search changes on the PR.

---

## 2. Regression traps (things we got wrong before → the check that catches it)

| # | What bit us | The check that catches it |
|---|---|---|
| T-01 | Version strings drifted across surfaces (v2.7.8: README/docs still said 2.6.2 / 2.7.1 / 2.5.5 while shipped was 2.7.7). | Gate A "version sweep" grep, every release. |
| T-02 | Docs/scripts/tests kept advertising indexes **after they were deleted from HF** (hybrid `1000G2021_HGDP`, `*-dictless`, `NGG_3_hg38` — the last since re-published, 2026-10-25, as the reference companion of `NGG_3_hg38+hg38_HGDP1kGP`). | Gate A "no retired index advertised" repo-wide + HF-card grep. |
| T-03 | The **HF dataset card** was never part of the consistency sweep → it still declared AGPL, a fictional index layout, and the wrong `pinellolab/` repo-id. | HF card is now an explicit surface in Gate A. |
| T-04 | Web **first-paint** seeded the variant default via `_preferred_variant()` **without the genome**, so the phased-preference was skipped and it could pick a non-HGDP1kGP panel. | Both call sites pass the genome; Gate A "default index" + a web walkthrough confirming the proposed default. |
| T-05 | `crisprme.py` path-correction hard-stripped a 3-char `"bin"` suffix → the `opt/` copy spawned by the web got a garbage PostProcess path (`ModuleNotFoundError: validate_inputs`), breaking **every** web assembly-search. | Layout-detection (PostProcess-dir probe) instead of suffix math; web assembly-search e2e in the walkthrough. |
| T-06 | Hidden web dropdowns used `dcc.Dropdown` whose `value` didn't round-trip → assembly-search submitted `None` and no job was submittable. | Use `dcc.Store` for hidden state; assembly-search submit e2e. |
| T-07 | Web search didn't pass `--gene_annotation` → the gencode companion was discarded, so **closest-gene columns were missing** from web jobs. | Closest-gene columns asserted in the web-job report during the walkthrough. |
| T-08 | A single multi-chromosome custom VCF produced a cryptic `Genomes/<g>/.fa` failure (the enricher needs per-chrom files); a SNP-only VCF aborted the whole build at the indel step. | `_ensure_per_chromosome_vcfs` auto-split + graceful 0-indel skip; custom-VCF build test. |
| T-09 | CI **path filters** mean a `pages/**`-only or docs-only PR silently skips `validate-benchmarks` (and, until 2.8.1, `unit-tests`) — a PR can look "green" without running them. | Gate B note: confirm the full suite actually *ran* (the version-bump PR touches `crisprme.py`, which re-triggers it). `test_ci_workflows` pins the paths `unit-tests` and `web-e2e` must list. |
| T-10 | `apptainer exec` bypasses the image entrypoint, so `crisprme.py` isn't on PATH. | Use `apptainer run` (or call `/opt/conda/bin/python` directly) in all HPC/test invocations + docs. |
| T-11 | HF uploads/deletes 401 when the token comes from the env; `HfFolder` was removed in new `huggingface_hub`. | Pass `token=get_token()` explicitly; delete via `create_commit([CommitOperationDelete…])`. |
| T-12 | Editing a bash script **while instances of it are running** corrupts the running jobs (broke chr7/chr9 in a GW re-phase). | Never edit a running script — copy to a new filename and run that. |
| T-13 | A fresh `complete-search` with an empty index-path builds the search index but **not** the `Dictionaries/` tiers → variant/indel post-analysis crashes. | New-user path = `download` (or `build-index-only`) FIRST, then `complete-search --index-path`; documented in the quickstart. |
| T-14 | `samplesID` files that over-list samples (vs VCF-genotyped) silently **deflate every reported AF**. | Pre-flight `check_samplesid_overlisting` WARN; AN denominator sanity check. |
| T-15 | Indel off-targets were silently dropped by a `grep -w $chrom` vs `fakechr` mismatch in `post_analisi_indel.sh`. | `grep -F -w "$fakechrom"`; indel off-targets asserted present in the e2e. |
| T-16 | Read-only `vuoto.txt` in the SIF made an **unconditional** `truncate` fail and abort the pipeline. | Guard writes/`truncate` on existence + writability before touching a file. |
| T-17 | A test run via `apptainer exec --cleanenv` + explicit `/opt/conda/bin/python` gives a minimal PATH that EXCLUDES `/opt/conda/bin`, so external-binary-dependent tests (BEDOPS `sort-bed`/`closest-features`, `bcftools`) spuriously "fail" even though the binary ships in the image. | When running the suite out-of-entrypoint, pass `--env PATH=/opt/conda/bin:/usr/bin:/bin` (or use `apptainer run`); confirm a reported missing-binary failure isn't just PATH before treating it as a release bug. |
| T-18 | The same fact stated in multiple README/website sections drifted apart — the default-index sample counts disagreed across sections (`929/3162` vs the correct `925/3166`; a stale `3,477` retired-hybrid panel size in the report-summary example). | The section-by-section docs/website review (see Gate B) + a grep that every default-index count/panel-size agrees; single-source the number and point other sections at it. |
| T-19 | A public install page handed users the WRONG package's subcommands — the website Conda path told pre-2.2.0 stock-Bioconda users to run CRISPRme+-only `download`/`web-interface`, which that package lacks → the documented step errors. | Every install path's example commands must match the subcommands the package it just installed actually has (stock Bioconda = legacy `setup`; CRISPRme+ Docker/source = `download`/`web-interface`). |
| T-20 | Overlay testing bound `crisprme.py` to `/opt/conda/opt/crisprme/crisprme.py`, but the EXECUTED wrapper is `/opt/conda/bin/crisprme.py` (a COPY, not a symlink) — so current-main `crisprme.py` wrapper changes (help text, arg defaults, new subcommands) were silently NOT exercised; tests ran the image's old wrapper. | To test current-main `crisprme.py` in a SIF, bind BOTH `--bind <src>/crisprme.py:/opt/conda/bin/crisprme.py` AND the `opt/crisprme/{PostProcess,pages}` trees. (Modules/pages load from `opt/`; the entrypoint runs the `bin/` copy.) |
| T-21 | A documented download command pointed at data NOT hosted on HF — `download --what vcf --dataset 1000G\|HGDP` appeared as a working example on 7 doc surfaces (+ a `~16 GB` size claim) + the web Settings action, but the raw per-dataset `vcfs/` folders were never uploaded, so the command errors. Doc examples outran what was actually published. | Gate A "every documented `download --what X` resolves on HF" (`list_repo_files`). Only hosted paths may appear as copy-pasteable examples; for anything fetched-from-HF, verify the target exists, and make the CLI error point to the real path (`--what index` / `build-inputs/`). |
| T-22 | A documented command instructed a flag value the shipped validator REJECTS — METHODS told users `--max-total-edits -1` to "disable" the cap, but the validator errors on any negative (`crisprme.py:1834`), so the documented command aborts; a stale in-code comment reinforced it. | For every flag VALUE shown in docs (especially sentinels like `-1`/`0`/`all`), confirm the shipped validator accepts it; keep the doc, the in-code comment, and the validator in agreement. |
| T-23 | METHODS claimed a capability was "not yet" shipped (mega indels "not yet frequency-annotated") though the code + CHANGELOG showed it ships — stale "not yet / forthcoming" claims understate shipped capability in an IND reference. | Grep METHODS + docs for "not yet / forthcoming / will land / in preparation" and verify each against the shipped code/CHANGELOG before tagging (also a fresh-user-audit heuristic). |
| T-24 | A per-release doc carried HARDCODED pinned image tags `pinellolab/crisprme:vX.Y.Z` in copy-paste run commands (the web-interface user guide, 5×) instead of `:latest`, so every release silently left the guide pointing at the PREVIOUS image. | In the version sweep, grep for hardcoded `:vX.Y.Z` tags in runnable doc commands — each must either be on the bump list OR use `:latest`. A bump that edits only `crisprme.py`+`Dockerfile` (CI derives the image tag from the git tag) MISSES ~27 hardcoded `X.Y.Z` doc strings → triage every `grep -rn "<OLD-version>"` hit into (a) current-version claim → bump, (b) pinned run-command tag → bump/`:latest`, (c) historical/"retired in"/"RESOLVED in"/CHANGELOG → keep. |
| T-25 | The only container-wrapped CLI-search instruction for HPC used `apptainer exec` (bypasses the entrypoint → `crisprme.py: not found`) while the same docs warn to use `apptainer run`, and NO correct `apptainer run complete-search`/`generate-report` worked example existed anywhere. | `grep -rn "apptainer exec.*crisprme.py"` must return 0 outside an explicit "wrong-way" example; at least one full worked `apptainer run ... crisprme.py complete-search` + `generate-report` script must exist for the HPC/CLI user. |
| T-26 | The clickable-app / web-form DEFAULT first search selected the 64 GB whole-genome variant index, OOM-killing a typical laptop (exit 137, empty log) with no clickable recovery; and a RAM escape hatch the docs promised ("~16 GB single-chromosome variant search") had no actual WEB-FORM control. | New-user gate: the first-run default must fit ~16 GB (reference-only) with an explicit variant-aware upgrade prompt, OR the app reads container-allocated memory and warns-before-Submit below the measured peak, AND an OOM/exit-137 dialog points to raise-memory / reference-only / single-chromosome. Any RAM escape hatch the docs promise must have a real web-form control, not only a CLI recipe. |
| T-27 | A web-stack smoke that imported `flask_caching.backends.filesystem` FAILED even on the correct, working image (flask-caching 2.4.1) — the module is `filesystemcache`, and the app deliberately uses the CLASS NAME `CACHE_TYPE: "FileSystemCache"` (app.py), NOT the legacy lowercase `filesystem` path (which raises `ImportStringError` by design). Testing the wrong path gives a false "broken web app" verdict (and I nearly shipped it as a Bioconda recipe test that would have red-ed CI). | To smoke-test the web/Flask stack (image or recipe), import what the APP uses: `python -c "from flask_caching.backends import FileSystemCache; import cachelib"` — never `import flask_caching.backends.filesystem`. Verify against what `app.py` actually configures (`CACHE_TYPE`), not an assumed module path. |
| T-28 | The web **Load Example** button set the "Maximum edits" total-edit cap to **4** while loading mm=4/bDNA=1/bRNA=1 — so `--max-total-edits=4` silently PRUNED every 4-mismatch-plus-bulge off-target (bulges had zero headroom at the full mismatch budget). The slider's own fresh-load default was correctly 6, but Load Example and two fallback defaults (5, 4) were all inconsistent with it — a user clicking Load Example got a tighter search than the documented default without any indication. | Every default edit-budget value (slider default, Load-Example cap, and all `else` fallbacks) must equal the recommended default `mm + bDNA + bRNA` (= 6 for 4/1/1); a total cap **below** that sum prunes bulged sites. Regression test `test/web/test_max_edits_default.py` asserts the Load-Example cap == mm+bDNA+bRNA == 6. Also: memory-banner / edit-budget / any "default" the UI seeds must be checked against the CLI + docs recommended default, not an independent magic number. |
| T-29 | The web **"Analysis mode"** selector (population-level vs per-sample) stayed visible on the **Personal-assembly** tab, where it is meaningless — per-sample resolves genotypes across a multi-sample variant *panel*, but a personal-assembly search works on one individual's own two haplotypes (no panel to resolve), and assembly-search never reads `search-mode`. A user on that tab saw a toggle that does nothing. | Web controls that only apply to the variant-PANEL path (Analysis mode = population-level/per-sample) must be hidden when the Personal-assembly tab is active. The `genome-mode-tabs` active-tab callback (`toggle_assembly_search_panel`) toggles the `analysis-mode-container` style along with the submit-field swap. When adding any new panel-only control, wire it into this same tab toggle (and vice-versa: anything genuinely shared — guide/PAM/mm/bulges/annotations — stays visible on both tabs). |
| T-30 | A runtime WARNING ("odd-length partially-degenerate PAM + bulges can crash the engine, issue #105, fixed in CRISPRitz 2.7.1") kept firing on the **default** search long after the fix shipped: `NRG` is odd-length and carries the partially-degenerate `R`, so the gate (`bMax>0 and pam_len%2==1 and any(c in "WRYSKMBDHV")`) tripped on every default run — crying wolf about a crash that cannot occur once CRISPRitz is pinned past the fix. A warning that references "fixed in CRISPRitz X" but never checks the installed version becomes a permanent false alarm once X ships. | Any runtime warning/guard that references "fixed in CRISPRitz/dependency X" must be re-evaluated whenever the pinned dependency crosses X: if the pin is hard (`>=X` on all surfaces — Docker build arg, `install_from_source.sh`, Bioconda), remove the warning (no supported install can hit it) rather than leave it firing. Grep for "crash"/"issue #"/"fixed in CRISPRitz" in `crisprme.py` + `pages/` each release and confirm each still describes a *reachable* condition on the pinned stack. |
| T-31 | A `--per-sample` search aborted at the FINAL annotation stage on **Docker Desktop for Mac/Windows**: the post-search step annotates three target files in PARALLEL, each building the same tabix index; `flock` is non-functional on VM-backed bind mounts (virtiofs/9p → ENOTTY), so the concurrent index builds raced, corrupted the read, and tripped the stderr-fatal gate with a misleading "file may be corrupt, re-download". Linux/cluster (working `flock`) was unaffected, so CI + ml003/ml008 never reproduced it — only a real Mac Docker Desktop run did. | Pre-build the tabix index ONCE, serially, before the parallel annotators (they then only read a complete index); tolerate a `flock` that errors on such filesystems; `load_annotation_bed` retries the BGZF-EOF guard + index build + open with backoff. Regression: `PostProcess/test_annotation.py` (8 mock tests, wired into CI). Each release, run at least one `--per-sample` search on **Docker Desktop for Mac** (virtiofs bind mount), not only on Linux — the race is filesystem-specific. |
| T-32 | The report told **two different threshold stories**: Section 4's recommended-panel PROSE said "validate every CRISPR-Bulge ≥1.9e-4 (95%-recall) site + the high-CFD rescue," but the shipped `panel_top100.tsv` was still built CFD-centric (hard-include `CFD≥0.5 at ≤1 bulge`, CRISPR-Bulge only via the fill ranking), and the CRISPR-Bulge tier TABLE bottomed out at an orphan `≥0.1` with **no sensitive row** — so the panel the user downloaded didn't implement the panel the report recommended, and the sensitive (CB≥1.9e-4, 95%-recall) operating point appeared in the risk card but not the Section-4 table/download. | The panel now hard-includes `(CRISPR-Bulge≥1.9e-4, 95%-recall) ∪ (mm+b≤2) ∪ (CFD≥0.2 at 0 bulges) ∪ (CFD≥0.5 at 1 bulge) ∪ clean-0mm-bulge ∪ must-includes` (the shipped panel = the recommended default panel; a strict **superset** of the old one — verify `in OLD not NEW == 0` on real data), and the CB tier table shows `≥0.5 / ≥0.2 / ≥1.9e-4 (95% recall)` + a `crispr_bulge_gt_0.tsv` download. ONE threshold message across the risk card, the "how to read this" modal, the Section-4 methods note, the score glossary, and METHODS §3.5. Each release: a threshold-consistency sweep (grep every CFD/CRISPR-Bulge cut point across `cb_operating_point.py` (single source of truth) + `generate_report.py` + `risk_card.py` + `docs/` + the collaborator rationale doc — they must agree) + `test_generate_report` assertions (the sensitive (95%-recall) tier ROW renders, panel superset, Option-C floors `PANEL_RESCUE_CFD=0.2`/`PANEL_HIGH_CFD=0.5`). |
| T-33 | The docs tell Apptainer users `apptainer run crisprme.sif crisprme.py complete-search …` ("use run, not exec" — T-10/T-25), but that only works if the SIF carries a **runscript** (from the Docker `CMD ["crisprme.py"]`) AND `/opt/conda/bin` is on the container `PATH`. A SIF obtained by any route OTHER than `apptainer pull docker://…` (e.g. built from a docker-archive/sandbox, or a locally-baked RC) can DROP the runscript and leak the host `PATH`, so `apptainer run` no-ops ("no runscript defined") and bare `apptainer exec crisprme.py` fails ("command not found") — the documented HPC command then doesn't run. | At Gate C, verify against a **fresh `apptainer pull docker://pinellolab/crisprme:vX.Y.Z`**: `apptainer inspect --runscript` references `crisprme.py`, `apptainer run crisprme.sif crisprme.py --version` prints `vX.Y.Z`, and `$PATH` inside includes `/opt/conda/bin`. When testing a SIF built another way (e.g. a baked RC on an apptainer-only node), invoke with `apptainer exec --env PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin crisprme.py …` and do NOT treat the missing runscript as a code bug — it's a SIF-provenance artifact of how that SIF was created, not the shipped image. |
| T-34 | The mega sites-only INDEL build driver (`seq_script/merge_panels/mega_build_indels.sh`) advertised an overridable `PAM=${PAM:-…}` but then HARDCODED `--pam …/20bp-NRG-SpCas9.txt` (line 57) and the `NRG_3_hg38+hg38_mega_INDELS` destination dir names — so any non-NRG mega indel rebuild aborted (missing NRG file) or silently built the WRONG PAM under the NRG path. The three phased per-dataset scripts threaded `$PAM` correctly; only the mega driver drifted. | A new-PAM build of ANY shipped index must honor `--pam $PAM` end-to-end: `grep -rn "20bp-NRG\|NRG_3_" seq_script/` must find only overridable defaults (`${PAM:-…}` / `${IDXNAME:-NRG_3}`), never a hardcoded literal in a `--pam`/dest path. Driver now parameterizes `IDXNAME` + guards loudly if crispritz produces a dir the name didn't expect. |
| T-35 | Multiple candidate source-VCF dirs exist per dataset (current + `archive/` + a dated `_20260827` snapshot + a decoy `hg38_combined_five_datasets`); publishing the WRONG one as "exact provenance" would mislabel what built the shipped index. The mega decoy had **1.83× too many records** (122.4M vs the shipped 66.9M) and a different INFO scheme. | Before publishing ANY build input, VERIFY it matches the shipped index: sample list (`bcftools query -l` == samplesID, md5 where possible) + record/indel counts (`bcftools index -n` sum) == the installed index's `Dictionaries/registry_<vcf>/variant_count.json`. A count/sample MISMATCH blocks publication. Canonical matched paths are recorded in `paper/reproducibility/variant-datasets/` + §3 below. |
| T-36 | CFD's DNA-bulge half of the Doench gap matrix was DEAD: the DNA-bulge caller front-trimmed the alignment, collapsing the gap to a table position with no key → `except: pass` free-pass (factor 1.0), so every DNA bulge (and terminal bulge) scored like no bulge — a PAM-proximal bulge scored as high as a perfect match (chr19 CFD 0.87 vs CRISPR-Bulge 0). | CFD now scores the displayed alignment with the position-weighted gap factor at the bulge's true position (terminal columns extrapolated); mismatch-only + in-range RNA bulges are bitwise-unchanged (`test_bulge_cfd_penalty`, `test_cfd_cleanroom_equiv`). The deployment calibration was refit on the re-scored benchmark + PR curves regenerated (incl. a bulge-stratified panel marking the report's own thresholds). Each release that touches scoring: re-verify the chr19-class PAM-proximal bulge is penalized AND mismatch-only CFD is unchanged, and that `PostProcess/calibration/*.npz` + `pr_curves.png` were regenerated from the current scorer. |
| T-37 | A local/Mac `docker build` of the image failed at the conda solve with an SSL/connection error because the **local/corporate network blocks `conda.anaconda.org`** — so the from-scratch RC image + the T-31 Mac `--per-sample` test couldn't be built on-site, and it read as a code/image defect rather than a network block. | The Dockerfile parameterizes the channel (`ARG CONDA_CHANNEL_BASE`, default `https://conda.anaconda.org`); on such a network build with `docker build --build-arg CONDA_CHANNEL_BASE=https://prefix.dev …` (prefix.dev mirrors conda-forge + bioconda; `scorer_env._channel_args()` + the Dockerfile both honor it). CI on GitHub runners is unaffected (default channel reachable). Documented on the Gate-B RC-build step. |
| T-38 | RE-OCCURRENCE of the T-16 read-only-`vuoto.txt` class via a NEW path: the 2.8.0 annotation **serial index pre-build** (`submit_job…:984`, the flock-race fix) called `annotation.load_annotation_bed($annotation_file)` DIRECTLY, bypassing `annotation.py main()`'s `empty` guard. For a **reference-only search** `$annotation_file` is the empty placeholder `vuoto.txt`, and `load_annotation_bed`→`_reindex`→`pysam.tabix_index()` **bgzips it in place and REMOVES the plain file** → the later gene-annotation `cp vuoto.txt` fails ("cannot stat"), aborting the search. It went green on `main` and failed only on the branch (the pre-build was new). SUBTLE CI trap: the `validate-benchmarks` RO-smoke `chmod -R a-w`s the code dir to mimic Apptainer, but **Docker runs as root and root ignores the perm bits**, so the in-place bgzip SUCCEEDS in Docker (consuming vuoto.txt) where a REAL read-only SIF would fail the write — same end result, opaque cause. | `load_annotation_bed` now returns `None` for `vuoto.txt` (never lock/reindex/bgzip the placeholder; callers already take the no-annotation path, `annotate_target` accepts `None`), AND the shell pre-build skips it (`if basename != vuoto.txt`). `test_annotation.test_vuoto_placeholder_never_indexed_or_consumed` pins it (returns None, `tabix_index` never called, no `.gz`/`.tbi`/`.tbi.lock`). On any change to annotation/index-prebuild code: reproduce the **reference-only search on a read-only install** with the REAL CI base image (`docker run pinellolab/crisprme:latest` + overlay the branch, chmod a-w, reference-only search) — do NOT trust an old local SIF (an outdated base image masks the regression), and remember Docker-root bypasses `chmod a-w`. |
| T-39 | CFD was scored by code inherited from the original CRISPRme: three `calc_cfd` copies with their `revcom` / `get_mm_pam_scores` helpers and the binary `mismatch_score.pkl` / `PAM_scores.pkl` tables. The copies had already drifted apart once (the INDEL path's `N`-in-PAM divergence), and the tables were binary pickles rather than reviewable data. 2.8.0 replaced them with CRISPRme+'s own implementation (`PostProcess/cfd_score.py`, Doench 2016, reviewable JSON tables), bit-for-bit identical in output. | Gate A sub-check "CFD clean-room provenance" below: the grep finds no inherited body or `.pkl`, `test_cfd_cleanroom_equiv` reproduces the frozen `cfd_golden.tsv.gz` (6,132 cases) bit for bit, and `scripts/generate_cfd_tables.py --check` validates the JSON tables. |
| T-40 | `NRR_3_hg38+hg38_mega` was published (2026-10-27) with an **empty** `_INDELS` (24 sub-indexes, 0 `.bin`) and header-only `log_indels`, so every indel search on it found nothing — and installing it next to `NRG_3_hg38+hg38_mega` REPLACED that index's correct 419 MB indel logs with the empty ones. It was built by an ad-hoc script from a fake-indel genome that predates the sites-only fix (#191) instead of through `seq_script/merge_panels/mega_build_indels.sh`; `buildTST` crashed on each empty genome and nothing checked. The release's own index audit verified builder content and readability, not indel completeness. The same archive check found `NRR_3_hg38+hg38_HPRC`'s registry predating the `data_type` stamp, so it silently defaulted to population-level on a phased panel; both came from ml007 build dirs that missed the later registry backfill. | Gate A "every published variant archive can find indels", run on the archive itself; build any new-PAM mega index only through `mega_build_indels.sh` or from its kit genome. A search cannot catch this archive: 2.8.1 fails loudly on a `.bin`-less sub-index only when that chromosome's log lists indels, and here the logs were empty too, so it reads exactly like a dataset with no indels. The archive check is the safeguard. |
| T-41 | A skipped test looks like a pass. Every documented self-test and `validate-benchmarks.yml` set `CRISPRME_SKIP_HEAVY=1` on the claim that the bulge-2 index builds need more than 16 GB of RAM, never measured (peak 3.66 GB; passes in 16 GiB, 4 CPUs, no swap); the two benchmarks so skipped could hardly catch a lost site anyway (one reference empty, the other 2 sites). The light CI job's step for `post_process.sh` ran one test and skipped it, 8 other tests skipped there too, and one module was in no workflow: all passed in the image, where no job ran them. Six function-based modules report `Ran 0 tests`, `OK` under `python -m unittest`. | Gate B's local gate: only the three named skips are allowed, the six modules run as scripts, every benchmark runs. A new claim that something does not fit is measured before it becomes a skip. `test_ci_workflows` pins the workflows and the docs' expected output to the registry's count. |
| T-42 | The engine checks ran on the old engine only. `validate-benchmarks.yml` runs `test/index_flex/check_index_flex.py`, the SNP+indel oracle (`test/indel_snp_oracle/run_oracle.sh`) and `scripts/check_alt_pam.py` with the branch overlaid on the published `:latest`, so for 2.8.1 on CRISPRitz 2.9.0 while the release ships 2.9.1: index_flex took its refusal branch for both K = 8 cells and still printed INDEX FLEX PASS, and the oracle never saw the compiled enricher the release ships. 2.8.1 had two regressions that only the engine bump caused (a chromosome with SNPs but no indels, a chromosome whose indels hold no site of the PAM). | Gate B's local gate, step 4b: the three run on the release image (amd64 and arm64), and the K = 8 cells must read "same results as". `test_ci_workflows` pins the step and its greps to what the scripts print. |

**Gate A sub-checks added from this pass** (wire into section A above):
- Read-only reference-only search (T-38, re-occurrence of T-16): a reference-only `complete-search` (no `--annotation`, so the `vuoto.txt` placeholder is used) on a READ-ONLY install must emit `integrated_results` and NOT consume/need to write `vuoto.txt`. Covered by the `validate-benchmarks` RO-smoke; any NEW code that indexes/pre-builds an annotation must skip the `vuoto.txt` placeholder. Validate against the REAL current base image, not a stale local SIF (Docker-root bypasses `chmod a-w`; a true read-only FS or an Apptainer SIF is the faithful test).
- Build-input provenance + publication (T-35, §3): the HF `vcfs/<dataset>/` tree hosts EXACTLY the matched build inputs of the 3 shipped variant panels (`hg38_HGDP1kGP`, `hg38_HPRC`, `hg38_mega`, the panels every published variant index is built on), each verified (sample list + record counts) against its shipped `variant_count.json`; every documented `download --what vcf --dataset <X>` resolves (closes T-21). Provenance table lives in `paper/reproducibility/variant-datasets/`.
- Scoring-change re-verify (T-36): on any scorer change, `test_bulge_cfd_penalty` + `test_cfd_cleanroom_equiv` green, calibration npz + PR curves regenerated, and the report's threshold markers match METHODS §3.
- CFD clean-room provenance (T-39, 2.8.0): CFD is CRISPRme+'s own implementation (`PostProcess/cfd_score.py`, Doench 2016) scored via `cfd_score_bulge_aware`; there is **no inherited `calc_cfd`/`revcom`/`get_mm_pam_scores` body and no `mismatch_score.pkl`/`PAM_scores.pkl`** anywhere in the shipped tree (production OR tests). Gate on each release: `grep -rnE "def (calc_cfd|revcom|get_mm_pam_scores)|(mismatch_score|PAM_scores)\.pkl" PostProcess scripts` returns nothing; `test_cfd_cleanroom_equiv` reproduces the committed golden fixture (`cfd_golden.tsv.gz`) bit-for-bit; `python scripts/generate_cfd_tables.py --check` validates the JSON tables (392 mismatch + 16 PAM, anchors). The golden fixture is the frozen historical oracle — the clean-room must reproduce it exactly, bit-for-bit.
- Sentinel-value multi-surface sweep (T-22): `grep -rn "disables the cap\|-1.*cap"` must find `-1` only on the ENV VAR `CRISPRME_IUPAC_CAP` (accepts it), never on the CLI flag `--max-total-edits` (rejects negatives). The T-22 fix once landed on METHODS/CHANGELOG/comment but missed README:884.
- Default-index count single-source (T-18): `grep -rnE "929|3,?162|3,?477"` across docs/ must return ONLY raw-source / mega-source-list / retired / 2019-AN-example contexts; any "929 HGDP + 3,162 1kGP" DEFAULT split is the stale pair (correct = 925 + 3,166; 925 = post-QC, 929 = raw source count).
- No retired-index copy-paste even under a Historical banner: a whole-file banner is not cover for a verbatim dead `download --index-name <retired>` command — replace the inline name with the live equivalent.
- Help-text value resolvability: every VALUE in a flag's `--help` must resolve against the DEFAULT repo/config or carry an inline "custom only" caveat (e.g. `--what vcf` help listing `1000G, HGDP`).
- Scripted HF-card + catalog check (make mandatory, not skipped offline): the HF `indexes/` tree of each dataset repo must equal the set [`INDEX_PROVENANCE.md`](INDEX_PROVENANCE.md) documents for it, and the card must carry zero retired names + MGBOA frontmatter.

> **How to use T-rows:** every one is a real incident. When a release pass uncovers a new one,
> add a row here with the symptom and the check, and wire the check into Gate A or B above.

---

## 3. Building & publishing a new index (new dataset, or a new PAM for an existing dataset)

Anyone can reproduce a shipped index or build a new one from the **public build inputs** — the source
VCFs + samplesID + PAMs are hosted on HuggingFace (`crisprme.py download --what vcf --dataset <name>`)
and every dataset's exact upstream provenance + retrieval/QC recipe is in
[`paper/reproducibility/variant-datasets/`](../paper/reproducibility/variant-datasets/). The shipped
indexes, their matched source VCFs, and the exact build command are the table there; `docs/PRECOMPUTED_INDEXES.md`
lists the downloadable indexes.

### What a PAM change does and doesn't cost
A CRISPRitz index is built for ONE PAM (encoded in the name `<MOTIF>_<N>_<genome>+<enriched>`, e.g.
`NRG_3_…`). To search a different PAM you rebuild ONLY the index, not the data:

- **Reused verbatim (PAM-independent):** the enriched genome (`Genomes/hg38+hg38_<ds>/` + `_INDELS/`),
  ALL `Dictionaries/` tiers (registry / genotypes / log_indels), `samplesIDs/`, `variant_count.json`.
- **Rebuilt (PAM-specific):** only `genome_library/<MOTIF>_<N>_hg38+hg38_<ds>{,_INDELS}` and the reference
  `genome_library/<MOTIF>_<N>_hg38`. The reference index is dataset-independent → build it ONCE per PAM
  and reuse across every variant index with that PAM.

### Add a new PAM file
One line: `<20 N's><PAM-motif> <signed-pos>` — `pos > 0` = 3′ PAM (last |pos| chars, e.g. SpCas9);
`pos < 0` = 5′ PAM (first |pos| chars, e.g. Cas12a). Filename MUST be `<len>-<motif>-<name>.txt` (≥3
hyphen fields) or the parser errors. Examples in `PAMs/`: `20bp-NRG-SpCas9.txt` (`…NRG 3`),
`20bp-NNN-NO-PAM.txt` (`…NNN 3`, the PAM-less index that captures non-canonical PAMs / SpRY-like enzymes).
> **CFD/CRISPR-Bulge caveat for a new PAM:** CFD is computed only for a 20-nt guide + 3-nt 3′ PAM; it
> down-weights the off-target PAM by Doench's **SpCas9** ladder (NGG 1.0, NAG 0.26, NCG 0.11, … most
> others 0). A 4-nt/2-nt/5′ PAM builds + searches but returns NO CFD; and for a non-SpCas9 enzyme
> (SpRY/SpCas9-NG/Cas12a) both scores' PAM terms are SpCas9-specific — the index **finds** the sites,
> but rank by edit distance + the enzyme's own PAM. See METHODS §3.

### Build (the real user command — no shortcuts)
```bash
# 1. get the public inputs (reference bundle + the dataset's source VCFs + samplesID)
crisprme.py download --what all
crisprme.py download --what vcf --dataset hg38_HGDP1kGP        # the matched build input from HF
# 2. build the index for your PAM (enriches the genome, builds the index, builds the Dictionaries);
#    --samplesID takes a LISTING of samplesID file names, not the samplesID file itself
printf 'hg38_HGDP1kGP.samplesID.txt\n' > hgdp1kgp.samplesID.listing.txt
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NNN-NO-PAM.txt \
    --bDNA 2 --bRNA 2 --thread 16 --vcf VCFs/hg38_HGDP1kGP --samplesID hgdp1kgp.samplesID.listing.txt
# -> genome_library/NNN_3_hg38+hg38_HGDP1kGP{,_INDELS}  (reuses the enriched genome + Dictionaries on a rebuild)
# 3. search it
crisprme.py complete-search --pam PAMs/20bp-NNN-NO-PAM.txt --index-path genome_library … --guide guides.txt …
```
The mega sites-only indel half uses `seq_script/merge_panels/mega_build_indels.sh` (honors `PAM=`/`IDXNAME=`
since T-34); the three phased per-dataset scripts (`{hgdp1kgp,onekg2021,hprc}_build.sh`) thread `--pam`
directly.

### Release checklist for a NEW published index
- [ ] **Provenance verified (T-35)** — the source VCF matches the shipped/intended index (sample list +
      record counts vs `variant_count.json`); the canonical matched path recorded in the reproducibility folder.
- [ ] **Inputs on HF** — `vcfs/<dataset>/` + samplesID + the PAM file uploaded; `download --what vcf --dataset <ds>`
      resolves (closes T-21); size + SHA256 recorded.
- [ ] **PAM honored end-to-end (T-34)** — no hardcoded PAM/dest literals in the build path; index named `<MOTIF>_<N>_…`.
- [ ] **New index documented** — added to `docs/PRECOMPUTED_INDEXES.md` + the reproducibility table (dataset,
      PAM, upstream, exact command, checksums) + the HF dataset card; `BUILDING_INDEXES`/§3 build command shown.
- [ ] **Smoke** — a real `complete-search` on the new index renders a report; for a non-NGG PAM confirm the
      CFD/CRISPR-Bulge caveat is surfaced (NO-SCORES where expected), not silently wrong.
- [ ] **Surface consistency (Gate A)** — version/default-index/scorer/licence unchanged; the new index is
      ADDED, not advertised as *the* default unless that's the decision.
