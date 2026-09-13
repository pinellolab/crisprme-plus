# CRISPRme+ one-line install (clickable app)

For a scientist who does **not** want to use the terminal: one command installs a
clickable **CRISPRme** app that manages Docker, downloads the reference data, and opens
the web interface. After the install command, no terminal is ever needed again.

## For users

**macOS / Linux** (Terminal, once):
```bash
curl -fsSL https://pinellolab.github.io/CRISPRme/install.sh | bash
```

**Windows** (PowerShell, once):
```powershell
irm https://pinellolab.github.io/CRISPRme/install.ps1 | iex
```

Then open **CRISPRme** from Applications (macOS) or the Desktop / Start Menu
(Windows) — a small window with **three buttons: Start / Update / Stop**. Click
**Start**: the first time it **asks where to store the data** (~45 GB) — accept the
default in your home folder, or pick another disk/folder if your home disk is low
on space; the choice is **remembered** for every future launch. It then downloads
the reference + variant data automatically (once) and opens the web interface at
http://localhost:8080; every Start after that is instant. (Manage additional data
from the web app's own Settings.)

**Memory:** most searches are light, but the **whole-genome variant-aware search**
(across 1000G + HGDP) needs **~64 GB of RAM** — give Docker Desktop that much under
Settings → Resources → Memory. Many current laptops ship with 36–128 GB and can run
it; otherwise use a workstation or HPC. **Reference-only or single-chromosome /
target-region variant searches fit ~16 GB**, so any laptop can do those. Requires
[Docker Desktop](https://www.docker.com/products/docker-desktop/) (the installer
guides you if it is missing).

## Maintainers — canonical source

**There is a single canonical copy of the installer, and it is NOT in this repo.**
The two scripts are **served from GitHub Pages** and their source of truth is:

- `docs/install.sh` and `docs/install.ps1` in the **`pinellolab/CRISPRme`** repository
  (served at `https://pinellolab.github.io/CRISPRme/install.{sh,ps1}` — the exact URLs
  the one-liners above fetch).

This directory previously kept a **duplicate** of those two scripts, which silently
drifted out of sync (it lagged to an old image tag). To prevent that, the duplicate has
been removed — edit the installer only in the `pinellolab/CRISPRme` repo.

- The apps are thin managers that shell out to `docker compose`/`docker run` and open the
  browser; the engine is the published `pinellolab/crisprme` image. **Bump the image tag in
  those two scripts at every release** — or, better, point them at the floating
  `pinellolab/crisprme:latest` tag (the release CI already pushes `:latest`), so the
  installer never goes stale and the **Update** button always fetches the current release.
- MVP is unsigned (the installer strips the macOS quarantine flag since a user-run script
  placed the app); a later phase can ship a signed/notarized Tauri build.
