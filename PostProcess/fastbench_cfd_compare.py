#!/usr/bin/env python3
"""Worst-case CFD check: fast's strongest CFD per locus must be >= slow's.

After the max-CFD-representative fix, fast mode should report a CFD at each off-target
LOCUS that is >= the slow enumeration's max CFD there (fast never UNDER-states the
worst-possible CFD -- the safety-relevant direction). Reads two *.bestCFD.txt outputs.
Usage: fastbench_cfd_compare.py <slow.bestCFD.txt> <fast.bestCFD.txt>
"""
import sys


def _load(path):
    hdr, idx, loc_max = None, {}, {}
    for line in open(path):
        if line.startswith("#") or (hdr is None and "Chromosome" in line):
            hdr = line.rstrip("\n").lstrip("#").split("\t")
            for n in ("Chromosome", "Cluster_Position", "Direction", "CFD"):
                if n in hdr:
                    idx[n] = hdr.index(n)
            continue
        if hdr is None:
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(idx.values()):
            continue
        loc = (p[idx["Chromosome"]], p[idx["Cluster_Position"]], p[idx["Direction"]])
        try:
            c = float(p[idx["CFD"]])
        except ValueError:
            continue
        if loc not in loc_max or c > loc_max[loc]:
            loc_max[loc] = c
    return loc_max


def main():
    slow, fast = _load(sys.argv[1]), _load(sys.argv[2])
    # tolerance: floats formatted to 3 decimals in the output; allow 1 ulp of that.
    TOL = 5e-4
    under = []
    for loc, s in slow.items():
        f = fast.get(loc)
        if f is not None and f + TOL < s:
            under.append((loc, round(s, 3), round(f, 3)))
    print("loci compared:            %d (slow) vs %d (fast)" % (len(slow), len(fast)))
    print("loci where fast CFD < slow CFD (worst-case UNDER-report): %d" % len(under))
    for u in sorted(under, key=lambda x: x[1] - x[2], reverse=True)[:15]:
        print("   UNDER (locus, slow_cfd, fast_cfd):", u)
    # also report how many loci fast strictly RAISED (the fix working)
    raised = sum(1 for loc, s in slow.items()
                 if fast.get(loc) is not None and fast[loc] > s + TOL)
    print("loci where fast CFD > slow CFD (fix surfaced a stronger worst case): %d" % raised)
    print("\nRESULT:", "PASS (no CFD under-report)" if not under else "FAIL")
    sys.exit(0 if not under else 1)


if __name__ == "__main__":
    main()
