#!/usr/bin/env python3
"""Emit the observed-haplotypes sidecar from the alt-alignments table.

The report's ``build_observed_haplotypes_companion`` builds the clickable per-site
observed-haplotype modals by scanning the (often 30+ GB) ``*_all_results_with_
alternative_alignments.tsv`` once per render -- yet it keeps only rows whose
highest-CFD OR highest-CRISPR_BULGE axis is ``alt``-origin. This script does that
row filter ONCE, at search time, writing a small ``*_observed_haplotypes_raw.tsv``
sidecar so the report reads a few MB instead of the whole alt file.

Safety: the sidecar is a strict SUPERSET of the rows the companion keeps -- it
filters only on origin (``alt`` on either axis), NOT on the variant-set realness or
the top-1000/merge-window scoping the report still applies. The report re-applies
its exact per-row filter, site clustering, and >=2-haplotype test on the sidecar, so
``observed_haplotypes.tsv`` is byte-identical whether built from the sidecar or the
full alt file. Rows are copied VERBATIM (raw lines, no re-quoting), so a sidecar row
is byte-identical to its alt-file row.

If the two origin columns are absent (an older/odd header), the whole file is copied
(safe degrade -- the report then filters everything itself, as before).

Usage: make_oh_sidecar.py <alt_alignments.tsv> <out_observed_haplotypes_raw.tsv>
"""
import sys

CFD_ORIGIN = "REF/ALT_origin_(highest_CFD)"
CB_ORIGIN = "REF/ALT_origin_(highest_CRISPR_BULGE)"


def _origin_indices(header_line):
    """Return (i_cfd, i_cb) 0-based column indices, or (None, None) if absent."""
    cols = header_line.rstrip("\n").rstrip("\r").split("\t")
    try:
        return cols.index(CFD_ORIGIN), cols.index(CB_ORIGIN)
    except ValueError:
        return None, None


def _keep(line, i_cfd, i_cb):
    """True if the highest-CFD OR highest-CB origin column equals 'alt'
    (case/space-insensitive) -- the superset the report refines further."""
    parts = line.split("\t")
    if i_cfd < len(parts) and parts[i_cfd].strip().lower() == "alt":
        return True
    if i_cb < len(parts) and parts[i_cb].strip().lower() == "alt":
        return True
    return False


def write_sidecar(alt_path, out_path):
    """Stream alt_path -> out_path keeping alt-origin rows. Returns (total, kept)."""
    total = kept = 0
    with open(alt_path, "r", newline="") as fin, open(out_path, "w", newline="") as fout:
        header = fin.readline()
        if not header:
            return 0, 0
        fout.write(header)
        i_cfd, i_cb = _origin_indices(header)
        if i_cfd is None:
            # columns absent -> copy verbatim (report filters as before)
            for line in fin:
                total += 1
                fout.write(line)
            return total, total
        for line in fin:
            total += 1
            if _keep(line, i_cfd, i_cb):
                fout.write(line)
                kept += 1
    return total, kept


def main(argv):
    if len(argv) != 3:
        sys.stderr.write("usage: make_oh_sidecar.py <alt.tsv> <out.tsv>\n")
        return 2
    total, kept = write_sidecar(argv[1], argv[2])
    frac = (kept / total) if total else 0.0
    sys.stderr.write(
        "make_oh_sidecar: total=%d kept=%d (%.2f%%)\n" % (total, kept, 100.0 * frac))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
