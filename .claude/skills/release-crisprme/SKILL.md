---
name: release-crisprme
description: Cut a new CRISPRme release and keep GitHub, Docker, and Bioconda in sync. Use when the user wants to release, publish, tag, or bump the version of CRISPRme; update the Bioconda recipe; or update the changelog for a new version. Covers the version-consistency pre-flight, CHANGELOG.md, the GitHub release/tag, the Bioconda recipe update (autobump or manual PR), and final version verification.
---

# Release CRISPRme

Cut a CRISPRme release and keep three artifacts in sync:

```
git tag vX.Y.Z ──▶ GitHub release tarball ──▶ Bioconda meta.yaml (sha256 of that tarball)
        │                                              │
        └── version pinned in crisprme.py ─────────────┘
```

Version-sync points (ALL must equal the same `X.Y.Z`):
- `crisprme.py` → `version = "X.Y.Z"` (also sets `__version__`)
- `Dockerfile` → its header comment only, `# CRISPRme+ X.Y.Z — Python 3.11 image, built FROM SOURCE.`,
  which pins nothing and which the helper leaves: bump it by hand. There is NO
  `ARG crisprme_version`. The image is built FROM SOURCE (`COPY . ${PREFIX}/opt/crisprme/`), so
  the only pin is `crisprme.py`. The Dockerfile does carry `ARG crispritz_ref=vA.B.C`, the
  CRISPRitz pin and never CRISPRme's version: the exact CRISPRitz tag the release was tested
  with, which a release bumps when CRISPRitz moves; it equals `install_from_source.sh`'s
  `CRISPRITZ_REF` default. The recipe's bounded pin is `crispritz >=X.Y.0,<X.(Y+1)`, whose lower
  bound can be below it (`docs/RELEASING.md`, Version policy).
- git tag / GitHub release → `vX.Y.Z`
- Bioconda `recipes/crisprme/meta.yaml` → `{% set version = "X.Y.Z" %}` + fresh `sha256` + `build.number: 0`

Version policy (Luca, 2026-11-07; `docs/RELEASING.md`): CRISPRme and CRISPRitz share their
major.minor (2.10, then 2.11, ...) and their patch numbers are independent. CRISPRme X.Y.z
requires CRISPRitz X.Y.*: the recipe pins `crispritz >=X.Y.0,<X.(Y+1)`, and `crispritz_ref` is a
CRISPRitz X.Y tag. A breaking change in either package (an index format, the enricher's output)
bumps the minor of BOTH. So a release whose minor differs from its `crispritz_ref`'s is wrong:
stop and fix it, or release CRISPRitz first. When CRISPRitz gets a patch (decision D5, Luca,
2026-11-08): `crispritz_ref` and `CRISPRITZ_REF` pin the exact CRISPRitz tag the release was
tested with; the recipe's lower bound stays X.Y.0 and moves up only when a CRISPRitz patch fixes
something CRISPRme depends on (then to that patch, e.g. `>=2.10.1,<2.11`). So the lower bound may
be below `crispritz_ref`, never above it.

Bioconda source URL is `https://github.com/pinellolab/crisprme-plus/archive/refs/tags/v{{ version }}.tar.gz`,
so the recipe pulls the **GitHub release tarball** — its sha256 is the hash of that exact tarball.

Run all commands from the CRISPRme repo root. Replace `X.Y.Z` with the real version (no leading `v` except in tags/URLs).

---

## Step 0 — Pre-flight: version-consistency check

Confirm where the repo currently stands and that nothing is half-bumped.

```bash
echo "crisprme.py: $(grep -E '^version *= *' crisprme.py)"
echo "CRISPRitz pin: $(grep -E '^ARG crispritz_ref=' Dockerfile)"
echo "Dockerfile header: $(grep -E '^# CRISPRme\+ ' Dockerfile)"
echo "latest tag:  $(git tag --sort=-v:refname | head -1)"
echo "latest release: $(gh release list --repo pinellolab/crisprme-plus -L1)"
curl -sL https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/crisprme/meta.yaml | grep -E 'set version|sha256|number:'
# the version policy: the release X.Y.Z and its CRISPRitz pin share X.Y
echo "CRISPRitz pin's major.minor, which must be the X.Y of the release: $(grep -E '^ARG crispritz_ref=' Dockerfile | sed -E 's/.*=v?([0-9]+\.[0-9]+)\..*/\1/')"
```

Expected before a release: `crisprme.py`, the Dockerfile's header comment, the latest tag, and the recipe all show the **previous** version, consistent with each other, and the Dockerfile's `crispritz_ref` (a CRISPRitz version, never CRISPRme's) equals `install_from_source.sh`'s `CRISPRITZ_REF` default. If they disagree, fix that first. **Do not proceed with inconsistent state.**

## Step 1 — Bump in-repo version pins (helper script)

The helper edits ONLY `crisprme.py` (it finds no `ARG crisprme_version` in the source-built `Dockerfile` and says so), is idempotent, and prints the recipe + changelog scaffolds. It leaves the release the Dockerfile's header comment names, so bump that by hand (the `sed` below). The sha256 fetch only works AFTER the tag exists, so run it now with `--no-download`, then again after Step 3 to get the hash.

```bash
python scripts/prepare_release.py X.Y.Z --no-download
sed -i.bak -E 's/^# CRISPRme\+ [0-9]+\.[0-9]+\.[0-9]+ /# CRISPRme+ X.Y.Z /' Dockerfile && rm Dockerfile.bak
git --no-pager diff --stat        # must show ONLY crisprme.py + Dockerfile (its header comment)
git --no-pager diff crisprme.py Dockerfile
```

Rollback if wrong: `git checkout crisprme.py Dockerfile`.

## Step 2 — Update CHANGELOG.md

Move items from `## [Unreleased]` into a new dated section, using the script's scaffold as a template. Keep-a-Changelog categories: Added / Changed / Deprecated / Removed / Fixed / Security.

```bash
python scripts/prepare_release.py X.Y.Z --no-download | sed -n '/CHANGELOG.md scaffold/,/releases\/tag/p'
```

Edit `CHANGELOG.md`:
- Add `## [X.Y.Z] - YYYY-MM-DD` with the curated entries.
- Leave a fresh empty `## [Unreleased]` at the top.
- Update the link-reference footer:
  `[X.Y.Z]: https://github.com/pinellolab/crisprme-plus/releases/tag/vX.Y.Z`
  and repoint `[Unreleased]` to `.../compare/vX.Y.Z...HEAD`.

## Step 3 — Commit, then create the GitHub release + tag

Commit the version bump + changelog on the release branch/`main` (per repo policy), push, then create the release. `gh release create` creates the annotated tag `vX.Y.Z` for you.

```bash
git add crisprme.py Dockerfile CHANGELOG.md
git commit -m "Release vX.Y.Z"
git push origin HEAD

# Notes come straight from the changelog section you just wrote:
awk '/^## \[X.Y.Z\]/{f=1;next} /^## \[/{f=0} f' CHANGELOG.md > /tmp/crisprme_notes.md

gh release create vX.Y.Z \
  --repo pinellolab/crisprme-plus \
  --title "CRISPRme vX.Y.Z" \
  --notes-file /tmp/crisprme_notes.md \
  --target main            # or the branch you released from
```

Verify the tarball is now published:

```bash
curl -sIL https://github.com/pinellolab/crisprme-plus/archive/refs/tags/vX.Y.Z.tar.gz | head -1  # expect 200
```

Rollback: `gh release delete vX.Y.Z --repo pinellolab/crisprme-plus --cleanup-tag` (only if the release must be withdrawn before Bioconda picks it up).

## Step 4 — Update the Bioconda recipe

First compute the sha256 of the now-published tarball (this is the recipe's `source.sha256`):

```bash
python scripts/prepare_release.py X.Y.Z    # no --no-download; prints sha256 + meta.yaml diff
# or directly:
curl -sL https://github.com/pinellolab/crisprme-plus/archive/refs/tags/vX.Y.Z.tar.gz | shasum -a 256
```

You have two routes. **Prefer autobump; fall back to a manual PR.**

### Route A — Autobump bot (preferred, low-effort)

Bioconda runs an autobump bot that periodically scans GitHub releases and opens
`Update crisprme to X.Y.Z by BiocondaBot` PRs on its own. Wait ~a day, then check:

```bash
gh pr list --repo bioconda/bioconda-recipes --search "crisprme in:title" --state all -L 5
```

If a PR exists, review that it has the right version + sha256 + `build.number: 0` and the bounded `crispritz` pin this release needs (the bot moves only the version and the hash), wait for CI/lint (green), then a Bioconda member or you (once approved) comment on the PR:

```
@BiocondaBot please merge
```

To nudge the bot to scan immediately instead of waiting, comment on any recipe PR/issue:

```
@BiocondaBot autobump crisprme
```

### Route B — Manual PR (reliable, when you don't want to wait)

```bash
# one-time: fork bioconda/bioconda-recipes to your account
gh repo fork bioconda/bioconda-recipes --clone --remote
cd bioconda-recipes
git checkout master && git pull upstream master
git checkout -b crisprme-X.Y.Z
```

Edit `recipes/crisprme/meta.yaml`:
- `{% set version = "X.Y.Z" %}`  ← bump
- `source.sha256:` ← the hash from above
- `build.number: 0`  ← RESET to 0 (only INCREMENT, keeping version, for a same-version rebuild)
- `crispritz >=X.Y.0,<X.(Y+1)`  ← CRISPRme X.Y.z requires CRISPRitz X.Y.* (the version policy: `>=2.10.0,<2.11` for 2.10.z; `docs/RELEASING.md`); keep it BOUNDED, and move it with the minor (within a minor, raise the lower bound only to a CRISPRitz patch that fixes something CRISPRme depends on: D5)
- Do NOT edit `source.url` — it uses `{{ version }}` and resolves automatically.

```bash
git commit -am "Update crisprme to X.Y.Z"
git push -u origin crisprme-X.Y.Z
gh pr create --repo bioconda/bioconda-recipes \
  --title "Update crisprme to X.Y.Z" \
  --body "Bump crisprme to vX.Y.Z. sha256 recomputed from the GitHub release tarball. build.number reset to 0."
```

Useful @BiocondaBot PR comments:
- `@BiocondaBot please update` — sync your PR branch with upstream master.
- `@BiocondaBot please add label` — request the `please review & merge` label.
- `@BiocondaBot please fetch artifacts` — links to the built test packages/containers.
- `@BiocondaBot please merge` — after lint + CI are green and the PR is approved, uploads packages and merges.

Wait for CI/lint to pass. Fix any lint findings the bot reports.

## Step 5 — Verify GitHub ↔ Bioconda ↔ Docker match

Do this AFTER the Bioconda PR merges and the package is on the channel (can take an hour+).

```bash
# GitHub / in-repo
grep -E '^version *= *' crisprme.py
grep -E '^ARG crispritz_ref=' Dockerfile   # CRISPRme's own version is pinned only in crisprme.py
grep -E '^# CRISPRme\+ ' Dockerfile        # ...and named in the Dockerfile's header comment
git tag --sort=-v:refname | head -1

# Bioconda published package
conda search -c bioconda crisprme | tail -3          # or: mamba search -c bioconda crisprme

# Optional: build the Docker image and check it installs the new version
# (the image is built from source, so this does NOT wait on Bioconda)
docker build -t crisprme:X.Y.Z . && \
  docker run --rm crisprme:X.Y.Z crisprme.py --version   # expect vX.Y.Z
```

crisprme.py, the Dockerfile's header comment, the git tag and the bioconda version must read `X.Y.Z`; the Dockerfile's `crispritz_ref` reads the CRISPRitz pin (`vA.B.C`). Done.

---

## Failure / rollback quick reference
- Wrong in-repo bump: `git checkout crisprme.py Dockerfile CHANGELOG.md`.
- Wrong/premature GitHub release: `gh release delete vX.Y.Z --repo pinellolab/crisprme-plus --cleanup-tag` and redo.
- sha256 mismatch on the Bioconda PR: GitHub can regenerate archive tarballs; recompute the sha256 from the live tarball and update `meta.yaml`, then push.
- Deps changed but version unchanged: keep `version`, INCREMENT `build.number` (do not reset to 0).
- Docker build fails building CRISPRitz: check `ARG crispritz_ref` points at a tag that exists.
