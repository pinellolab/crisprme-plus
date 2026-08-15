# CRISPRme Docker Quickstart — run the web interface in a few commands

> ⚠️ **Alpha release.** CRISPRme+ is a preview of the next major version; interfaces, data
> layouts, and results may change between builds. For production or clinical/critical work,
> use the frozen stable line **[CRISPRme 2.1.14](https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.14)**.

This is the fastest way to get CRISPRme running with its point-and-click **web
interface**, with **no conda, no compiling, and no 410 GB download**. You copy a
few commands, open a browser, fill in a short form, and read the results.

Everything runs inside Docker, so the only thing you install is Docker itself.

---

## 1. Install Docker (one time)

- **Mac / Windows:** install **Docker Desktop** — https://docs.docker.com/get-started/
- **Linux:** install **Docker Engine** — https://docs.docker.com/engine/install/

Then, in Docker Desktop → **Settings → Resources**, give Docker enough memory:
**16 GB** is fine for a first run / reference-only searches, but the default
genome-wide 1000G+HGDP variant search is memory-intensive — give it **at least
32 GB (64 GB recommended)**.

> **Disk:** the batteries-included setup below needs **≈ 250 GB free**. The download
> is ~44 GB, but the combined 1000G+HGDP variant index **expands to ~190 GB on disk**
> — its per-sample variant dictionaries (for population/sample annotation) are
> ~170 GB uncompressed. A reference-only setup is far smaller (~15 GB).

Check Docker works, then pull the CRISPRme+ image:

```bash
docker run --rm hello-world

# Pull the current CRISPRme+ alpha (multi-arch: Apple Silicon + Intel/Linux).
docker pull pinellolab/crisprme:v2.2.0-alpha.24
```

> **Already have an older image?** Docker does **not** re-download a tag you already
> have — run `docker pull pinellolab/crisprme:v2.2.0-alpha.24` again to update. Skipping this makes an
> old image error with `download is not an allowed command`.

## 2. Make a folder to hold your data and results

Everything CRISPRme downloads and produces will live here (so it is kept between
runs). Pick any folder:

```bash
mkdir -p ~/crisprme && cd ~/crisprme
```

> In every command below, `-v "${PWD}:/DATA"` shares *this* folder with the
> container as `/DATA`. That is how your downloads and results are saved to your
> computer instead of disappearing when the container stops.

## 3. Download the reference data (one time, minutes — not hours)

This pulls the human genome, annotations, PAM files and sample lists from the
CRISPRme HuggingFace mirror (a fast CDN). It replaces the old multi-hour `setup`:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:v2.2.0-alpha.24 \
  crisprme.py download --what all --path /DATA
```

## 4. Download a prebuilt index (one time, minutes — skips a long build)

Bulge-enabled searches need a genome **index**. Building it yourself takes ~10
minutes of CPU; instead, download the ready-made SpCas9 (NGG) index:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:v2.2.0-alpha.24 \
  crisprme.py download --what index --index-name NGG_3_hg38 --path /DATA
```

For a variant-aware search (what the default web search uses), also download the
1000 Genomes + HGDP enriched index:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:v2.2.0-alpha.24 \
  crisprme.py download --what index --index-name NGG_3_hg38+hg38_1000G_HGDP --path /DATA
```

The web interface picks these up automatically — a search that uses the NGG PAM
with up to 2 bulges will reuse them instead of rebuilding. (Need a different
nuclease? See **"Installing more indexes"** at the bottom.)

## 5. (Optional, advanced) Add the raw 1000 Genomes VCFs

CRISPRme's superpower is finding off-targets created by genetic variants — and the
`NGG_3_hg38+hg38_1000G_HGDP` index you downloaded in step 2 **already** makes the
default web search variant-aware. You do **not** need the raw VCFs for that.

Download the raw 1000 Genomes variant set (~16 GB) only for CLI sample-level
analyses / personal risk cards:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:v2.2.0-alpha.24 \
  crisprme.py download --what vcf --dataset 1000G --path /DATA
```

## 6. Start the web interface

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA -p 8080:8080 -it \
  pinellolab/crisprme:v2.2.0-alpha.24 crisprme.py web-interface
```

`-p 8080:8080` connects the app inside the container to your browser. Leave this
running, and open:

**http://127.0.0.1:8080**

(If you started it on a remote server, use that server's address instead of
`127.0.0.1`.)

## 7. Run a search in the browser

1. Enter a **guide/spacer** sequence (e.g. `CTAACAGTTGCTTTTATCAC`).
2. Choose the **PAM** (e.g. `20bp-NGG-SpCas9`) and the **genome** (`hg38`).
3. Leave the default **Maximum edits** slider (3) for a quick search — it caps the
   total mismatches + DNA/RNA bulges — or open **Advanced options** to set
   mismatches / DNA bulges / RNA bulges individually (e.g. 4 / 1 / 1). The default
   **1000G+HGDP** variant index is pre-selected, so the search is variant-aware out
   of the box.
4. Give the job a name and click **Submit**. (Optionally tick **Notify me by
   email** to be sent a results link when the job finishes — set your SMTP details
   once under **Settings → Email notifications**.)
5. Watch the live status page; when it finishes, open the **Results** to explore
   the off-targets, scores and plots. Per-sample tabs (Summary by Sample and
   **Personal Risk Cards**) appear only when the search is variant-aware.

Functional annotations (ENCODE cCREs / SCREEN + DHS + GENCODE) are applied
automatically — the built-in bundle is enabled by default, so there is nothing to
pick on the search form. To add your own annotation BEDs or turn tracks on/off, use
**Settings → Data Manager → Manage annotations** (local mode only).

Your results are saved on your computer under `~/crisprme/Results/<job name>/`.

To stop the web server, press **Ctrl+C** in the terminal.

---

## Running on HPC with Singularity / Apptainer

On a cluster you usually cannot run Docker, but Singularity/Apptainer runs the
same image with no root and no daemon. Most clusters already have it (the command
is `apptainer`, or `singularity` on older systems — they are interchangeable).

```bash
# 1. build the image once (a ~2 GB .sif file; no root needed)
apptainer pull crisprme.sif docker://pinellolab/crisprme:v2.2.0-alpha.24

# 2. download data + a prebuilt index into a working folder
mkdir -p ~/crisprme && cd ~/crisprme
apptainer run --bind "${PWD}:/DATA" --pwd /DATA crisprme.sif \
  crisprme.py download --what all --path /DATA
apptainer run --bind "${PWD}:/DATA" --pwd /DATA crisprme.sif \
  crisprme.py download --what index --index-name NGG_3_hg38 --path /DATA
apptainer run --bind "${PWD}:/DATA" --pwd /DATA crisprme.sif \
  crisprme.py download --what index --index-name NGG_3_hg38+hg38_1000G_HGDP --path /DATA
# optional (advanced): the raw 1000G VCFs (~16 GB) — only for CLI sample-level
# analyses / personal risk cards; the index above already makes the web search
# variant-aware
apptainer run --bind "${PWD}:/DATA" --pwd /DATA crisprme.sif \
  crisprme.py download --what vcf --dataset 1000G --path /DATA

# 3. launch the web interface
apptainer run --bind "${PWD}:/DATA" --pwd /DATA crisprme.sif \
  crisprme.py web-interface
# then open http://127.0.0.1:8080 (or the node's address on a cluster)
```

Two Singularity-specific notes:
- **Use `apptainer run`, not `apptainer exec`.** `run` executes the image's
  entrypoint, which activates the conda environment so `crisprme.py` is on the
  PATH; `exec` skips it and you get `crisprme.py: not found`.
- **Networking is the host's.** Apptainer shares the host network, so there is no
  `-p` port mapping — the app is reachable directly at port **8080**, which must be
  free on that node. On a shared login node, run on a compute node / interactive
  session instead.

Everything else (the web form, results, downloading more indexes) is identical to
the Docker instructions above.

---

## Installing more indexes (as you need them)

An index is specific to a **PAM + bulge count + genome**. Download whichever you
need by its exact published name — for example the pamless variant index:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:v2.2.0-alpha.24 crisprme.py download --what index --index-name NNN_3_hg38+hg38_1000G_HGDP --path /DATA
```

To see which indexes are published, browse the dataset repository
(`lucapinello/crisprme-data`, folder `indexes/`) — the names there are exactly
what you pass to `--index-name`.

If an index you need is **not** published, the web interface will **not** build it
on the fly: it tells you no matching index is installed and asks you to install
one first. Pre-build it from the command line with `crisprme.py build-index-only`
(same `--genome`/`--pam`/`--bDNA`/`--bRNA` you plan to search with) — see
[`PRECOMPUTED_INDEXES.md`](PRECOMPUTED_INDEXES.md) — then it is reused on every
later search.

Or do all of this **from the browser**: the web interface has a **Settings /
Data Manager** page (the button on the home page, or `/settings`) where you can
add genomes, indexes, VCF datasets, annotations and PAMs — including downloading
a UCSC assembly by name, pre-building a variant-aware index, checking disk usage,
and deleting data you no longer need. See
[`SETTINGS_DATA_MANAGER.md`](SETTINGS_DATA_MANAGER.md).

## Troubleshooting

- **The genome/PAM dropdowns are empty** → you skipped step 3. Run
  `crisprme.py download --what all --path /DATA` (inside Docker, as above).
- **Browser can't connect to `127.0.0.1:8080`** → make sure the `docker run`
  in step 6 includes `-p 8080:8080` and is still running.
- **Results disappeared after closing the container** → make sure every command
  includes `-v "${PWD}:/DATA"` so files are written to your folder, not the
  throwaway container.
- **A search ran out of memory** → raise Docker's memory limit (step 1), or
  start with a smaller search (fewer mismatches/bulges, or one chromosome).
- **`Bind for 0.0.0.0:8080 failed: port is already allocated`** → another
  container is already using port 8080. Stop it (`docker ps`, then
  `docker stop <id>`), or map a different host port with `-p 8081:8080` and open
  http://127.0.0.1:8081.
- **Search finished but the results table is empty (or shows no variant
  off-targets)** →
  1. Thresholds too strict — raise **Maximum edits**, and if you opened Advanced
     options check the DNA and RNA bulges are not both `0`.
  2. Reference-only selected — keep the **1000G+HGDP** option (pre-selected by
     default) to get variant off-targets.
  3. Variant index not installed — re-run
     `crisprme.py download --what index --index-name NGG_3_hg38+hg38_1000G_HGDP --path /DATA`.
  4. Confirm success: `Results/<name>/log_error.txt` is empty and
     `*.integrated_results.tsv` is non-empty.

For the full web-interface reference (every form field, all the result tabs),
see [`crisprme_web_interface_user_guide.md`](crisprme_web_interface_user_guide.md).
