# CRISPRme+ one-line install (clickable app)

For a scientist who does **not** want to use the terminal: one command installs a
clickable **CRISPRme** app that manages Docker, downloads reference data, and opens
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
**Start**: the first time it downloads the reference + variant data automatically
(once), then opens the web interface at http://localhost:8080; every Start after
that is instant. (Manage additional data from the web app's own Settings.)

**Requirements:** Docker Desktop (the installer guides you if it is missing) and,
for the genome-wide *variant* search, 64 GB of RAM (reference-only fits 16 GB).

## For maintainers

- `install.sh` (macOS/Linux) generates `~/Applications/CRISPRme.app` with
  `osacompile` (built into macOS — no build toolchain), writes `~/CRISPRme/`
  (compose file + `crisprme-data/`), and strips the quarantine flag so the app
  opens with **no Gatekeeper warning** (it was placed by a script the user ran).
- `install.ps1` (Windows) writes a WinForms app script to `%LOCALAPPDATA%\CRISPRme`
  and Desktop + Start-Menu shortcuts.
- Both apps are thin managers that shell out to `docker compose` and open the
  browser; the CRISPRme engine is the published `pinellolab/crisprme:v2.4.0` image.
- This is the MVP (unsigned; the installer strips quarantine). A later phase can
  replace the app with a small signed/notarized Tauri build for standalone
  distribution.
- **Hosting:** these two files are published to GitHub Pages at
  `https://pinellolab.github.io/CRISPRme/` so the one-liners above resolve.
