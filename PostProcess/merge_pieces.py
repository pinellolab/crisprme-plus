#!/usr/bin/env python
"""Merge the pieces of one chromosome's post-analysis into the files one process writes
(2.10.0, option P).

pool_post_analisi_snp.py and pool_post_analisi_indel.py cut a chromosome's clustered raw
targets into contiguous pieces of lines and run new_simple_analysis.py (or
analisi_indels_NNN.py) once a piece, over the threads of the search (CRISPRME_POSTPROC_PIECE=1).
Each line is analysed on its own, so a piece writes exactly the rows one process writes for
its lines, in their order. This script joins them into the chromosome's files, byte for byte
those of one process:

* the CFD and fewest-edits files: the pieces' rows in piece order, one header;
* the CRISPR-Bulge file: the pieces leave its rows unscored (<piece>.cbulge_rows.jsonl, each
  row with its input line). One process scores them in batches: whenever its rows reach
  cbulge_rows.BATCH_ROWS, after the input line that reached it, and the rest at the end; the
  model's floats depend on which rows share a batch, so they are scored here in exactly those
  batches, with the same function (cbulge_rows.score_rows);
* the SNP companions (population summary, phase confirmation, SNP+SNP co-occurrence): one
  process keeps one row a variant off-target, the first it meets; the pieces' rows in piece
  order, the first of each identity kept;
* the indel companions (SNP+indel co-occurrence, the sites-only indel AF): rows in piece order;
* <out>.high_variant_density_regions.bed: one line a dense window (start, end, guide), the
  line resultIntegrator reads first where pieces both list a window (new_simple_analysis.
  _record_dense_window), sorted;
* <out>.four_allele_positions.bed: the sorted union; <out>.CFDGraph.txt: the sum.

Then it writes on stdout (log_verbose.txt) a line on the merge and the notes one process
ends with (analysis_notes): the pieces write none, but leave their counts in
<piece>.piece.json; the dense windows and the four-allele positions are counted in the
merged files.

Usage: merge_pieces.py snp|indel <output prefix> <piece prefix>...
"""

import json
import os
import sys

import analysis_notes
import cbulge_rows

_SNP_COMPANIONS = {
    "population_summary.tsv": ("Chromosome", "Position", "Direction", "crRNA", "DNA", "SNP"),
    "phase_confirmation.tsv": ("Chromosome", "Position", "Direction", "crRNA", "DNA", "SNP"),
    "snp_snp_cooc.tsv": ("Chromosome", "Position", "Direction", "crRNA", "DNA"),
}
_INDEL_COMPANIONS = ("indel_snp_cooc.tsv", "indel_af.tsv")


def _existing(paths):
    return [p for p in paths if os.path.exists(p)]


def concat_rows(dst, sources):
    """``sources``' rows after their first line, under the first one's header line."""
    with open(dst, "w") as out:
        for i, src in enumerate(sources):
            with open(src) as fh:
                head = fh.readline()
                if i == 0:
                    out.write(head)
                for line in fh:
                    out.write(line)


def merge_companion(dst, sources, identity):
    """A companion TSV of the pieces: the first source's block of '#' lines (the last one
    names the columns), then every source's rows in order, the first of each ``identity``
    (column names) kept; with ``identity`` None, every row."""
    seen = set()
    head = None
    idx = None
    with open(dst, "w") as out:
        for src in sources:
            with open(src) as fh:
                lines = fh.readlines()
            n_h = 0
            while n_h < len(lines) and lines[n_h].startswith("#"):
                n_h += 1
            if head is None:
                head = lines[:n_h]
                out.writelines(head)
                if identity is not None and head:
                    cols = head[-1].lstrip("#").rstrip("\n").split("\t")
                    idx = [cols.index(c) for c in identity]
            for line in lines[n_h:]:
                if idx is not None:
                    f = line.rstrip("\n").split("\t")
                    key = tuple(f[i] if i < len(f) else None for i in idx)
                    if key in seen:
                        continue
                    seen.add(key)
                out.write(line)


def merge_plain_companion(dst, sources):
    """A companion whose one header line is its first: the rows in order."""
    concat_rows(dst, sources)


def _bed_record(line):
    f = line.rstrip("\n").split("\t")
    key = (int(f[1]), int(f[2]), f[3])
    # resultIntegrator reads (n_variants, iupac, bits, greedy != "0") and keeps the first
    # sorted: the same order picks the line kept here
    rep = (f[4], f[6], f[7], f[9] != "0")
    return key, rep, line if line.endswith("\n") else line + "\n"


def merge_hvdr(dst, sources):
    """The high-variant-density BED of the pieces: one line a window, the one with the
    smallest (n_variants, iupac, bits, greedy), sorted by start, end and guide."""
    best = {}
    header = None
    for src in sources:
        with open(src) as fh:
            for line in fh:
                if line.startswith("#"):
                    header = header or line
                    continue
                if not line.strip():
                    continue
                key, rep, text = _bed_record(line)
                old = best.get(key)
                if old is None or rep < old[0]:
                    best[key] = (rep, text)
    with open(dst, "w") as out:
        if header is not None:
            out.write(header)
        for key in sorted(best):
            out.write(best[key][1])


def merge_four_allele(dst, sources):
    """The four-allele positions of the pieces: their sorted union."""
    header = None
    pos = {}
    for src in sources:
        with open(src) as fh:
            for line in fh:
                if line.startswith("#"):
                    header = header or line
                    continue
                f = line.rstrip("\n").split("\t")
                pos[(f[0], int(f[1]))] = line if line.endswith("\n") else line + "\n"
    with open(dst, "w") as out:
        if header is not None:
            out.write(header)
        for key in sorted(pos):
            out.write(pos[key])


def merge_cfd_graph(dst, sources):
    """The pieces' CFD graphs, summed (each a header and rows of integers)."""
    header = None
    total = None
    for src in sources:
        with open(src) as fh:
            h = fh.readline()
            header = header or h
            rows = [[int(x) for x in line.rstrip("\n").split("\t")] for line in fh
                    if line.strip()]
        total = rows if total is None else [[a + b for a, b in zip(r, s)]
                                            for r, s in zip(total, rows)]
    with open(dst, "w") as out:
        out.write(header or "")
        for r in total or []:
            out.write("\t".join(str(x) for x in r) + "\n")


def score_cbulge(dst, header_line, row_files, do_scores, batch_rows=cbulge_rows.BATCH_ROWS,
                 score_fn=None):
    """Write the CRISPR-Bulge file: ``header_line``, then the pieces' unscored rows
    (``row_files``, JSON lines of [input line, row], in order) scored in the batches one
    process forms: rows gather line by line, and once they reach ``batch_rows`` after a
    line, the batch is scored; the rest at the end. Returns the number of rows."""
    score = score_fn or (lambda rows: cbulge_rows.score_rows(rows, do_scores))
    n = 0

    def flush(rows, out):
        for t in score(rows):
            # remove count of tmp_mms
            t.pop(-2)
            out.write("\t".join(t) + "\t" + str(0) + "\n")

    with open(dst, "w") as out:
        out.write(header_line)
        buf = []
        for src in row_files:
            last = None
            with open(src) as fh:
                for line in fh:
                    ln, row = json.loads(line)
                    if last is not None and ln != last and len(buf) >= batch_rows:
                        flush(buf, out)
                        buf = []
                    last = ln
                    buf.append(row)
                    n += 1
            # a piece ends after its last line
            if len(buf) >= batch_rows:
                flush(buf, out)
                buf = []
        if buf:
            flush(buf, out)
    return n


def _do_scores(pieces):
    for p in pieces:
        meta = p + ".piece.json"
        if os.path.exists(meta):
            with open(meta) as fh:
                return bool(json.load(fh).get("do_scores", True))
    return True


def merge_snp(out, pieces):
    best_cfd = [p + ".bestCFD.txt" for p in pieces]
    with open(best_cfd[0]) as fh:
        header_line = fh.readline()
    concat_rows(out + ".bestCFD.txt", best_cfd)
    concat_rows(out + ".bestmmblg.txt", [p + ".bestmmblg.txt" for p in pieces])
    score_cbulge(out + ".bestCRISPR_BULGE.txt", header_line,
                 _existing(p + ".cbulge_rows.jsonl" for p in pieces), _do_scores(pieces))
    merge_cfd_graph(out + ".CFDGraph.txt", _existing(p + ".CFDGraph.txt" for p in pieces))
    merge_hvdr(out + ".high_variant_density_regions.bed",
               _existing(p + ".high_variant_density_regions.bed" for p in pieces))
    four = _existing(p + ".four_allele_positions.bed" for p in pieces)
    if four:
        merge_four_allele(out + ".four_allele_positions.bed", four)
    for name, identity in _SNP_COMPANIONS.items():
        srcs = _existing(p + "." + name for p in pieces)
        if srcs:
            merge_companion(out + "." + name, srcs, identity)


def merge_indel(out, pieces):
    best_cfd = [p + ".bestCFD_INDEL.txt" for p in pieces]
    with open(best_cfd[0]) as fh:
        header_line = fh.readline()
    concat_rows(out + ".bestCFD_INDEL.txt", best_cfd)
    concat_rows(out + ".bestmmblg_INDEL.txt", [p + ".bestmmblg_INDEL.txt" for p in pieces])
    score_cbulge(out + ".bestCRISPR_BULGE_INDEL.txt", header_line,
                 _existing(p + ".cbulge_rows.jsonl" for p in pieces), _do_scores(pieces))
    for name in _INDEL_COMPANIONS:
        srcs = _existing(p + "." + name for p in pieces)
        if srcs:
            merge_plain_companion(out + "." + name, srcs)
    acfd = _existing(p + ".acfd.txt" for p in pieces)
    if acfd:
        with open(acfd[0]) as src, open(out + ".acfd.txt", "w") as dst:
            dst.write(src.read())


def _metas(pieces):
    """The pieces' <piece>.piece.json, in piece order (those that hold notes' counts)."""
    out = []
    for p in pieces:
        path = p + ".piece.json"
        if os.path.exists(path):
            with open(path) as fh:
                meta = json.load(fh)
            if "notes" in meta:
                out.append(meta)
    return out


def _data_lines(path):
    """The lines of ``path`` that are neither blank nor '#' lines (0 for no file)."""
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith("#"))


def _all_lines(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def snp_notes(out, pieces):
    """The notes one process writes at the end of a chromosome's SNP post-analysis
    (analysis_notes), written from the pieces' counts and the merged files at ``out``."""
    metas = _metas(pieces)
    if not metas:
        return []
    chrom, first = metas[0]["chrom"], metas[0]["notes"]
    counts = analysis_notes.snp_counts()
    for meta in metas:
        analysis_notes.add_counts(counts, meta["notes"]["counts"])
    letters = []
    bed = out + ".high_variant_density_regions.bed"
    if os.path.exists(bed):
        with open(bed) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    letters.append(line.split("\t")[6])
    cap = first["iupac_cap"]
    n_over, n_hyper = analysis_notes.dense_window_counts(letters, cap)
    return (analysis_notes.density_lines(chrom, os.path.basename(out),
                                         _data_lines(out + ".four_allele_positions.bed"),
                                         n_over, n_hyper, cap)
            + analysis_notes.snp_lines(chrom, counts, first["precheck_on"],
                                       first["fast_mode"], first["max_degeneracy"]))


def indel_notes(pieces):
    """The notes one process writes at the end of a chromosome's indel post-analysis
    (analysis_notes), written from the pieces' counts."""
    metas = _metas(pieces)
    if not metas:
        return []
    counts = analysis_notes.indel_counts()
    for meta in metas:
        analysis_notes.add_counts(counts, meta["notes"]["counts"])
    return analysis_notes.indel_lines(metas[0]["chrom"], counts)


def merge_line(kind, out, pieces):
    """The log's line on a merge: the pieces, the rows, the companions."""
    metas = _metas(pieces)
    chrom = metas[0]["chrom"] if metas else "?"
    best = out + (".bestCFD.txt" if kind == "snp" else ".bestCFD_INDEL.txt")
    companions = []
    for name in (_SNP_COMPANIONS if kind == "snp" else _INDEL_COMPANIONS):
        path = out + "." + name
        if os.path.exists(path):
            # the SNP companions' header lines are '#' lines; an indel one's is its first
            rows = _data_lines(path) if kind == "snp" else max(0, _all_lines(path) - 1)
            companions.append("%s (%d row(s))" % (name, rows))
    return ("Post-analysis %s of %s: merged its %d piece(s), %d row(s)%s"
            % (kind.upper(), chrom, len(pieces), _data_lines(best),
               "; " + ", ".join(companions) if companions else ""))


def main(argv):
    if len(argv) < 4 or argv[1] not in ("snp", "indel"):
        sys.stderr.write("usage: merge_pieces.py snp|indel <output prefix> <piece prefix>...\n")
        return 2
    kind, out, pieces = argv[1], argv[2], argv[3:]
    (merge_snp if kind == "snp" else merge_indel)(out, pieces)
    analysis_notes.write([merge_line(kind, out, pieces)]
                         + (snp_notes(out, pieces) if kind == "snp" else indel_notes(pieces)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
