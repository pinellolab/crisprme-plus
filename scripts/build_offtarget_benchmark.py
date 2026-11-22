#!/usr/bin/env python3
"""Build the external off-target benchmark TSV the calibration and PR figures use.

Sources: the CRISPR-Bulge repository's external *evaluation* data (Yaish &
Orenstein, NAR 2024, doi:10.1093/nar/gkae428) -- **not** its CHANGE-seq training
data:

  * ``Refined_TrueOT.csv``            already binary (editing rate >0.1% => active)
  * ``Tsai_2015_GUIDE_seq.csv``       read counts, labelled by the repo rule
  * ``Chen_2017_GUIDE_seq.csv``       read counts, labelled by the repo rule
  * ``Listgarten_2018_GUIDE_seq.csv`` read counts, labelled by the repo rule

Repo label rule for the read-count GUIDE-seq files (``OT_deep_score_src/
dataset_utilities.py`` ``load_dataset``, ``read_threshold=100``): keep rows with
reads >= 100 (ACTIVE) or reads == 0 (INACTIVE); drop ambiguous 0 < reads < 100 and
off-targets containing a literal N.

``bulge_type``: a gap in the off-target is an RNA bulge, a gap in the sgRNA is a
DNA bulge.

``cb_heldout`` is the leakage flag and the **only** flag that means "held out":
1 when the guide is not in the CRISPR-Bulge training guide set (CHANGE-seq union
FullGUIDE-seq), else 0. CFD has no training set, so every row is CFD-fair
regardless of ``cb_heldout``; anything that evaluates CRISPR-Bulge must restrict
itself to ``cb_heldout == 1``.

This builder used to live outside the repository and used to emit a second flag,
``in_domain = 1 if n_bulge <= 1``, described as CRISPR-Bulge's single-bulge input
domain. Every row of this benchmark has ``n_bulge <= 1``, so the column was
constant 1 and a filter on it dropped nothing -- which is exactly what happened:
a precision-recall panel titled "(in-domain)" filtered on it, dropped none of the
2,940,588 rows, and so reported CRISPR-Bulge on the 6 guides it was trained on.
A flag that cannot be false is worse than no flag, so it is gone; use
``cb_heldout``.

Usage::

    python3 scripts/build_offtarget_benchmark.py \\
        --datasets <crispr-bulge-checkout>/files/datasets --out benchmark.tsv
"""
import argparse
import csv
import os
import sys

READ_THRESHOLD = 100
COLUMNS = ("guide", "offtarget_seq", "chrom", "pos", "strand", "n_mismatch",
           "n_bulge", "bulge_type", "label", "source", "cb_heldout", "crrna_aligned")
READCOUNT_SOURCES = (("Tsai_2015", "Tsai_2015_GUIDE_seq.csv"),
                     ("Chen_2017", "Chen_2017_GUIDE_seq.csv"),
                     ("Listgarten_2018", "Listgarten_2018_GUIDE_seq.csv"))


def _load_list(datasets, name):
    with open(os.path.join(datasets, name)) as fh:
        return {line.strip() for line in fh if line.strip()}


def bulge_type(offtarget, sgrna):
    """RNA bulge = gap in the off-target; DNA bulge = gap in the crRNA."""
    off_gap, sg_gap = offtarget.count("-"), sgrna.count("-")
    if off_gap and not sg_gap:
        return "RNA"
    if sg_gap and not off_gap:
        return "DNA"
    if off_gap and sg_gap:
        return "MIXED"
    return "X0"


class Benchmark:
    """Accumulates deduplicated sites; on a label conflict the active row wins."""

    def __init__(self, training_guides):
        self.training = training_guides
        self.rows = []
        self._seen = {}

    def emit(self, source, guide, off_aln, sg_aln, chrom, start, strand, n_mm, n_b, label):
        rec = dict(guide=guide, offtarget_seq=off_aln, chrom=chrom, pos=start,
                   strand=strand, n_mismatch=n_mm, n_bulge=n_b,
                   bulge_type=bulge_type(off_aln, sg_aln), label=label, source=source,
                   cb_heldout=0 if guide in self.training else 1, crrna_aligned=sg_aln)
        key = (guide, chrom, str(start), strand, off_aln)
        if key in self._seen:
            i = self._seen[key]
            if label == 1 and self.rows[i]["label"] == 0:
                self.rows[i] = rec
            return
        self._seen[key] = len(self.rows)
        self.rows.append(rec)


def _common(row):
    return (row["sgRNA"].strip(), row["Align.off-target"], row["Align.sgRNA"].strip(),
            row["chrom"].strip(), row["Align.chromStart"].strip(),
            row["Align.strand"].strip(), int(row["Align.#Mismatches"]),
            int(row["Align.#Bulges"]))


def build(datasets, read_threshold=READ_THRESHOLD):
    csv.field_size_limit(10 ** 7)
    bench = Benchmark(_load_list(datasets, "CHANGEseq_sgRNAs_list.csv")
                      | _load_list(datasets, "FullGUIDEseq_sgRNAs_list.csv"))
    with open(os.path.join(datasets, "Refined_TrueOT.csv")) as fh:
        for row in csv.DictReader(fh):
            if "N" in row["Align.off-target"].replace("-", ""):
                continue
            g, off, sg, chrom, start, strand, n_mm, n_b = _common(row)
            bench.emit("Refined_TrueOT", g, off, sg, chrom, start, strand,
                       n_mm, n_b, int(row["label"]))
    for source, fname in READCOUNT_SOURCES:
        with open(os.path.join(datasets, fname)) as fh:
            for row in csv.DictReader(fh):
                if "N" in row["Align.off-target"].replace("-", ""):
                    continue
                reads = float(row["reads"])
                if reads == 0:
                    label = 0
                elif reads >= read_threshold:
                    label = 1
                else:
                    continue  # ambiguous 0 < reads < threshold
                g, off, sg, chrom, start, strand, n_mm, n_b = _common(row)
                bench.emit(source, g, off, sg, chrom, start, strand, n_mm, n_b, label)
    return bench.rows


def write(rows, out):
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in COLUMNS])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--datasets", required=True,
                    help="CRISPR-Bulge checkout's files/datasets directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--read-threshold", type=int, default=READ_THRESHOLD)
    args = ap.parse_args(argv)
    rows = build(args.datasets, args.read_threshold)
    write(rows, args.out)
    held = sum(r["cb_heldout"] for r in rows)
    pos = sum(r["label"] for r in rows)
    sys.stderr.write(f"wrote {len(rows):,} rows ({pos:,} positive; {held:,} CRISPR-Bulge "
                     f"held-out) to {args.out}\n")


if __name__ == "__main__":
    sys.exit(main())
