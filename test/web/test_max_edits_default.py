#!/usr/bin/env python3
"""Regression: the web 'Maximum edits' total-edit cap must match the recommended
default (4 mm + 1 DNA + 1 RNA = 6), so bulged sites at the full mismatch budget are
never pruned. The 'Load Example' button once set the slider to 4 while loading
mm=4/bDNA=1/bRNA=1, which capped --max-total-edits at 4 and dropped every 4-mismatch
+ bulge off-target. Guards that (and the fresh-load slider default)."""
import sys
sys.path.insert(0, "/opt/conda/opt/crisprme")
from pages import main_page  # noqa: E402

_fail = []
def check(c, m):
    print(("PASS " if c else "FAIL "), m); sys.stdout.flush()
    if not c: _fail.append(m)

# Load Example: the returned per-type caps and the max-edits cap must be consistent
# (cap == mm + dna + rna) and equal the recommended 6.
r = main_page.load_example_data(1)
mm, dna, rna, maxedits = r[4], r[5], r[6], r[13]
check((mm, dna, rna) == (4, 1, 1), f"Load Example loads 4 mm + 1 DNA + 1 RNA (got {mm},{dna},{rna})")
check(maxedits == mm + dna + rna,
      f"Load Example max-edits cap == mm+dna+rna, so no bulged site is pruned (got {maxedits} vs {mm+dna+rna})")
check(maxedits == 6, f"Load Example max-edits == 6 (the recommended default) (got {maxedits})")

print(f"\n{'FAILED' if _fail else 'OK'}: {len(_fail)} failure(s)")
sys.stdout.flush(); sys.stderr.flush()
import os
os._exit(1 if _fail else 0)
