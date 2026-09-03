#!/usr/bin/env python3
"""Compare a FAST-mode vs SLOW-mode new_simple_analysis.py output (same cluster input).

Proves, on real data, the two 2.5.1 fast-mode properties:
  1. LOSSLESS DETECTION: every variant off-target LOCUS (chrom, cluster_position,
     direction) the slow enumeration reports is ALSO present in the fast output
     (fast is a superset -- worst-possible may surface a formable-but-unobserved combo).
  2. WORST-POSSIBLE: at each shared locus the fast representative's minimum Total
     (mismatches + bulges) is <= the slow path's minimum Total (fast never understates
     the worst case).
Also reports row counts + the reference-locus parity. Usage:
    fastbench_compare.py <slow.bestmmblg.txt> <fast.bestmmblg.txt>
"""
import sys


def _load(path):
    hdr = None
    idx = {}
    var_min = {}    # locus -> min Total over variant rows
    ref_loci = set()
    var_loci = set()
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
            snp = parts[idx["SNP"]].strip()
            try:
                total = int(parts[idx["Total"]])
            except ValueError:
                total = None
            if snp in ("n", "NA", "", "."):
                ref_loci.add(locus)
            else:
                var_loci.add(locus)
                if total is not None:
                    if locus not in var_min or total < var_min[locus]:
                        var_min[locus] = total
    return {"n_rows": n_rows, "ref_loci": ref_loci, "var_loci": var_loci,
            "var_min": var_min}


def main():
    slow = _load(sys.argv[1])
    fast = _load(sys.argv[2])
    print("=== FAST vs SLOW new_simple_analysis comparison ===")
    print("rows:              slow=%d  fast=%d  (collapse %.1fx)"
          % (slow["n_rows"], fast["n_rows"],
             (slow["n_rows"] / fast["n_rows"]) if fast["n_rows"] else float("nan")))
    print("variant loci:      slow=%d  fast=%d"
          % (len(slow["var_loci"]), len(fast["var_loci"])))
    print("reference loci:    slow=%d  fast=%d"
          % (len(slow["ref_loci"]), len(fast["ref_loci"])))

    # (1) lossless detection: every slow variant locus present in fast
    missed = slow["var_loci"] - fast["var_loci"]
    print("NO-MISS (slow variant loci absent from fast): %d" % len(missed))
    if missed:
        for m in list(missed)[:10]:
            print("   MISSED:", m)

    # reference locus parity
    ref_missed = slow["ref_loci"] - fast["ref_loci"]
    print("reference loci in slow absent from fast: %d" % len(ref_missed))

    # (2) worst-possible: fast min Total <= slow min Total at each shared locus
    understated = []
    for locus, s_total in slow["var_min"].items():
        f_total = fast["var_min"].get(locus)
        if f_total is not None and f_total > s_total:
            understated.append((locus, s_total, f_total))
    print("WORST-POSSIBLE violations (fast Total > slow Total): %d" % len(understated))
    for u in understated[:10]:
        print("   UNDERSTATED:", u)

    ok = (len(missed) == 0 and len(understated) == 0)
    print("\nRESULT:", "PASS (lossless + worst-possible)" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
