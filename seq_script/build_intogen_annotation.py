#!/usr/bin/env python3
"""Build the IntOGen (CC0) cancer-driver-gene annotation BED for CRISPRme+.

IntOGen (Integrative OncoGenomics) publishes a compendium of cancer *driver* genes.
The CURRENT release is CC0 (public domain) -> free for academic AND commercial use +
redistribution, so CRISPRme+ ships it as a **default-on** cancer-gene annotation (the
licence-free complement to COSMIC, which stays opt-in behind the licence gate).

IMPORTANT: only the CURRENT release is CC0. Releases 2016.5-2023.05.31 were CC-BY-NC
(non-commercial). Always build from the current CC0 drivers file
(https://www.intogen.org/download -> IntOGen-Drivers-<date>.zip, CC0).

SHIPPED-INDEX PROVENANCE: the intogen_drivers.hg38.bed.gz published for CRISPRme+ v2.5.4
was built from the IntOGen **2024.09.20** release (zip label; internal release date
2024-06-18), file Compendium_Cancer_Genes.tsv (633 unique driver gene symbols), LICENSE.txt
= CC0 1.0. GENCODE input: the shipped Annotations/gencode.protein_coding.bed.gz. Rebuilding
from a newer IntOGen release will drift the gene set — re-pin this line when you do.

This maps each driver gene SYMBOL to its genomic interval via the GENCODE annotation
CRISPRme already ships (gene_name= in the attribute column), and writes a BED whose
4th column is ``<SYMBOL>_INTOGEN`` -- the label CRISPRme's resultIntegrator parses into
the ``Annotation_INTOGEN`` report column (mirrors the ``_COSMIC`` convention).

Usage:
    python build_intogen_annotation.py \
        --drivers   drivers.tsv \                 # from IntOGen-Drivers-<date>.zip (CC0)
        --gencode   gencode.protein_coding.bed.gz \  # shipped GENCODE (gene_name= in col9)
        --out       intogen_drivers.hg38.bed.gz      # bgzips; also writes .tbi if tabix present

Notes:
- --drivers: any TSV with a gene-symbol column (auto-detected: SYMBOL / GENE / Gene /
  gene / symbol). One row per (gene, cohort/cancer-type) is fine -- symbols are deduped.
- Symbols absent from GENCODE (e.g. non-coding / alias-only) are reported as unmapped.
- Output is coordinate-sorted; requires ``bgzip`` (and optionally ``tabix``) on PATH.
"""

import argparse
import gzip
import os
import re
import subprocess
import sys

# Primary assembly contigs only (chr1..chr22, chrX/Y/M, with or without the "chr"
# prefix). GENCODE places some genes on ALT/scaffold contigs too (e.g. KN.../GL.../KQ...
# or chr*_alt); anchoring a driver gene there would put its interval on a contig the
# search never reports off-targets on, so it would silently never match. Restrict the
# gene-span computation to primary contigs.
_PRIMARY_CHROM = re.compile(r"^(chr)?([0-9]{1,2}|X|Y|M|MT)$")


def _open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def load_driver_symbols(drivers_tsv):
    """Unique driver gene symbols from an IntOGen drivers TSV (symbol column auto-detected)."""
    with _open(drivers_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lower = [h.strip().lower() for h in header]
        col = None
        for cand in ("symbol", "gene", "gene_symbol", "hugo_symbol", "genename"):
            if cand in lower:
                col = lower.index(cand)
                break
        if col is None:
            raise SystemExit(
                f"Could not find a gene-symbol column in {drivers_tsv} header: {header}"
            )
        symbols = set()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if col < len(parts) and parts[col]:
                symbols.add(parts[col].strip())
    return symbols


def gene_intervals_from_gencode(gencode_bed):
    """Map gene_name -> (chrom, min_start, max_end) by spanning all rows of each gene.

    Reads the ``gene_name=<SYMBOL>`` attribute from the GFF-style 9th column of the
    shipped GENCODE BED, spanning every feature row so each gene gets one interval.
    """
    spans = {}  # symbol -> [chrom, start, end]
    with _open(gencode_bed) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            if not _PRIMARY_CHROM.match(f[0]):
                continue  # skip ALT/scaffold placements (would never match search targets)
            attrs = f[-1]
            i = attrs.find("gene_name=")
            if i < 0:
                continue
            sym = attrs[i + len("gene_name="):].split(";", 1)[0].strip()
            if not sym:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            chrom = f[0]
            cur = spans.get(sym)
            if cur is None:
                spans[sym] = [chrom, start, end]
            elif cur[0] == chrom:  # same chromosome: widen the span
                if start < cur[1]:
                    cur[1] = start
                if end > cur[2]:
                    cur[2] = end
    return spans


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drivers", required=True, help="IntOGen CC0 drivers TSV")
    ap.add_argument("--gencode", required=True, help="shipped GENCODE bed(.gz) with gene_name=")
    ap.add_argument("--out", required=True, help="output intogen_drivers.hg38.bed.gz")
    args = ap.parse_args()

    symbols = load_driver_symbols(args.drivers)
    sys.stderr.write(f"[intogen] {len(symbols)} unique driver gene symbols\n")
    spans = gene_intervals_from_gencode(args.gencode)
    sys.stderr.write(f"[intogen] {len(spans)} gene intervals in GENCODE\n")

    rows, unmapped = [], []
    for sym in symbols:
        span = spans.get(sym)
        if span is None:
            unmapped.append(sym)
            continue
        chrom, start, end = span
        rows.append((chrom, start, end, f"{sym}_INTOGEN"))
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    sys.stderr.write(
        f"[intogen] mapped {len(rows)} / {len(symbols)} symbols; "
        f"{len(unmapped)} unmapped (first 10: {sorted(unmapped)[:10]})\n"
    )
    if not rows:
        raise SystemExit("[intogen] no rows mapped -- check the drivers/gencode inputs")

    plain = args.out[:-3] if args.out.endswith(".gz") else args.out
    with open(plain, "w") as out:
        for chrom, start, end, name in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    # bgzip (required) + tabix (optional)
    subprocess.run(["bgzip", "-f", plain], check=True)
    if args.out != plain + ".gz":
        os.replace(plain + ".gz", args.out)
    try:
        subprocess.run(["tabix", "-p", "bed", "-f", args.out], check=True)
        sys.stderr.write(f"[intogen] wrote {args.out} + .tbi ({len(rows)} rows)\n")
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.stderr.write(f"[intogen] wrote {args.out} ({len(rows)} rows; tabix skipped)\n")


if __name__ == "__main__":
    main()
