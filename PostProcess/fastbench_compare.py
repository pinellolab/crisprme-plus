#!/usr/bin/env python3
"""Compare a FAST-mode vs SLOW-mode new_simple_analysis.py output (same cluster input).

Proves, on real data, the two 2.5.1 fast-mode correctness properties. The unit is an
off-target LOCUS (chrom, cluster_position, direction), and fast mode is WORST-POSSIBLE:
it emits ONE strongest off-target per window, collapsing weaker per-haplotype variant
off-targets into the (stronger) reference off-target at that locus -- so the right
criterion is over ALL rows (variant OR reference), not variant-specific:
  1. LOSSLESS LOCUS COVERAGE: every off-target locus the slow enumeration reports is
     ALSO covered by fast (as a variant OR a reference row). Fast is typically a SUPERSET
     -- worst-possible surfaces formable-but-unobserved combinations too.
  2. WORST-POSSIBLE BOUND: at each shared locus fast's STRONGEST off-target (min Total =
     mismatches + bulges, over all its rows) is <= slow's strongest (fast never
     understates the worst case).
The variant-specific breakdown is printed for diagnostics (a large variant-only gap that
is fully covered by fast reference rows is EXPECTED, not a miss). Usage:
    fastbench_compare.py <slow.bestmmblg.txt> <fast.bestmmblg.txt>
"""
import sys


def _load(path):
    hdr = None
    idx = {}
    var_min = {}    # locus -> min Total over VARIANT rows (diagnostic)
    loc_min = {}    # locus -> min Total over ALL rows (worst-case at the locus)
    ref_loci = set()
    var_loci = set()
    all_loci = set()
    n_rows = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or (hdr is None and "Chromosome" in line):
                hdr = line.rstrip("\n").lstrip("#").split("\t")
                for name in ("Chromosome", "Cluster_Position", "Direction", "Total", "SNP"):
                    if name in hdr:
                        idx[name] = hdr.index(name)
                continue
            if hdr is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx.values()):
                continue
            n_rows += 1
            locus = (parts[idx["Chromosome"]], parts[idx["Cluster_Position"]],
                     parts[idx["Direction"]])
            all_loci.add(locus)
            snp = parts[idx["SNP"]].strip()
            try:
                total = int(parts[idx["Total"]])
            except ValueError:
                total = None
            if total is not None and (locus not in loc_min or total < loc_min[locus]):
                loc_min[locus] = total
            if snp in ("n", "NA", "", "."):
                ref_loci.add(locus)
            else:
                var_loci.add(locus)
                if total is not None and (locus not in var_min or total < var_min[locus]):
                    var_min[locus] = total
    return {"n_rows": n_rows, "ref_loci": ref_loci, "var_loci": var_loci,
            "all_loci": all_loci, "var_min": var_min, "loc_min": loc_min}


def main():
    slow = _load(sys.argv[1])
    fast = _load(sys.argv[2])
    print("=== FAST vs SLOW new_simple_analysis comparison ===")
    print("rows:              slow=%d  fast=%d  (collapse %.1fx)"
          % (slow["n_rows"], fast["n_rows"],
             (slow["n_rows"] / fast["n_rows"]) if fast["n_rows"] else float("nan")))
    print("off-target loci:   slow=%d  fast=%d" % (len(slow["all_loci"]), len(fast["all_loci"])))
    print("  variant loci:    slow=%d  fast=%d" % (len(slow["var_loci"]), len(fast["var_loci"])))
    print("  reference loci:  slow=%d  fast=%d" % (len(slow["ref_loci"]), len(fast["ref_loci"])))

    # diagnostic: variant-specific gap (EXPECTED to be covered by fast reference rows)
    var_missed = slow["var_loci"] - fast["var_loci"]
    covered_by_ref = var_missed & fast["ref_loci"]
    print("diagnostic: slow variant loci not a fast VARIANT locus: %d "
          "(of which covered by a fast REFERENCE row: %d)"
          % (len(var_missed), len(covered_by_ref)))

    # (1) lossless LOCUS coverage: every slow locus present in fast (variant OR reference)
    absent = slow["all_loci"] - fast["all_loci"]
    print("NO-MISS (slow off-target loci absent from fast entirely): %d" % len(absent))
    for m in list(absent)[:10]:
        print("   ABSENT:", m)

    # (2) worst-possible: fast's strongest (min Total over all rows) <= slow's, per locus
    understated = []
    for locus, s_total in slow["loc_min"].items():
        f_total = fast["loc_min"].get(locus)
        if f_total is not None and f_total > s_total:
            understated.append((locus, s_total, f_total))
    print("WORST-POSSIBLE violations (fast strongest weaker than slow): %d" % len(understated))
    for u in understated[:10]:
        print("   UNDERSTATED (locus, slow_min, fast_min):", u)

    ok = (len(absent) == 0 and len(understated) == 0)
    print("\nRESULT:", "PASS (lossless locus coverage + worst-possible bound)" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
