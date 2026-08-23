===============================================================================
  CRISPRme+  —  run the point-and-click web interface (no terminal needed)
===============================================================================

CRISPRme+ predicts CRISPR off-targets, including ones created by human genetic
variants, and produces a shareable report. This folder lets you run it by
double-clicking — no command line.

-------------------------------------------------------------------------------
BEFORE YOU START (one time)
-------------------------------------------------------------------------------
1. Install Docker Desktop (free):
      Mac / Windows:  https://www.docker.com/products/docker-desktop/
   Open Docker Desktop once and let it finish starting.

2. Give Docker enough memory:
      Docker Desktop  ->  Settings  ->  Resources  ->  Memory
        - 16 GB  is enough for reference-only searches
        - 64 GB  is required for the genome-wide VARIANT search (1000G + HGDP)
   (A typical laptop has 8-32 GB; the variant search needs a workstation/server.)

-------------------------------------------------------------------------------
HOW TO USE (double-click, in order)
-------------------------------------------------------------------------------
Mac users:  the files ending in  .command
Windows users:  the files ending in  .bat

  STEP 1 (once):   "1 - Download data"
                   Downloads the reference genome + index (~25 GB). One time.

  STEP 1b (once, optional): "1b - Download variant index"
                   Only if you want VARIANT-aware search and have 64 GB RAM.
                   Downloads the 1000G + HGDP index (~16 GB download, ~60 GB on disk).

  STEP 2 (each time): "2 - Start CRISPRme"
                   Starts CRISPRme and opens your browser at
                        http://localhost:8080
                   Keep the little window that appears OPEN while you work.
                   To stop CRISPRme, close that window.

-------------------------------------------------------------------------------
IN THE BROWSER
-------------------------------------------------------------------------------
- Click "Load Example" to try it, or paste your own 20-nt guide.
- The default PAM is 20bp-NRG-SpCas9 (covers NGG + NAG).
- Submit, wait for the job, then open the Results page and download the report
  (a ZIP whose report.html is the file you open).

-------------------------------------------------------------------------------
NOTES
-------------------------------------------------------------------------------
- All data and results stay in the "crisprme-data" folder next to these files,
  on your own computer. Nothing is uploaded.
- First run of "2 - Start CRISPRme" may take a minute while Docker warms up.
- Mac: if double-clicking a .command file is blocked, right-click it -> Open
  (the first time only), or run:  chmod +x "*.command"
- No 64 GB machine? Use reference-only mode (skip step 1b); it still finds the
  classic (non-variant) off-targets. For genome-wide variant runs, use an HPC /
  core-facility / cloud instance.

Full documentation: https://github.com/pinellolab/crisprme-plus
===============================================================================
