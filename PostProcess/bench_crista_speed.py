#!/usr/bin/env python
"""
Phase-0 CRISTA speed baseline (the incumbent scorer) — run in the crisprme env
(SIF) so the 276 MB RandomForest is present. Mirrors bench_cbulge.py's batch
sizes so the OT/sec numbers are directly comparable and the REPLACE gate
("CRISPR-Bulge CPU throughput >= CRISTA") can be decided.

    apptainer exec <sif> bash -lc '
      export PATH=/opt/conda/bin:$PATH
      cd /opt/conda/opt/crisprme/PostProcess
      python bench_crista_speed.py --sizes 1000,10000,100000'
"""
import argparse
import json
import os
import resource
import sys
import time


def peak_rss_mb():
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v / 1024.0 if sys.platform.startswith("linux") else v / (1024.0 * 1024.0)


def make_batch(n, seed=0):
    """Synthetic aligned SpCas9 pairs + 29-nt genomic context (23-mer +/- 3bp)."""
    import random

    rng = random.Random(seed)
    bases = "ACGT"
    sg, off, ctx = [], [], []
    guide = "".join(rng.choice(bases) for _ in range(20))
    for _ in range(n):
        proto = list(guide)
        # 0-5 mismatches
        for _ in range(rng.randint(0, 5)):
            j = rng.randrange(20)
            proto[j] = rng.choice(bases)
        proto = "".join(proto)
        pam = rng.choice(bases) + "GG"
        sg_aln = guide + "AGG"          # sgRNA aligned (guide + canonical NGG placeholder)
        off_aln = proto + pam           # 23-mer off-target aligned
        left = "".join(rng.choice(bases) for _ in range(3))
        right = "".join(rng.choice(bases) for _ in range(3))
        genomic29 = left + off_aln + right   # 3 + 23 + 3 = 29
        sg.append(sg_aln)
        off.append(off_aln)
        ctx.append(genomic29)
    return sg, off, ctx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="1000,10000,100000")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--out", default="bench_crista_out.json")
    args = ap.parse_args()

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import CRISTA_score

    # load-once (unzip + unpickle the 276 MB RF)
    t0 = time.time()
    _sg, _off, _ctx = make_batch(64, seed=1)
    _ = CRISTA_score.CRISTA_predict_list(_sg, _off, _ctx)  # warm-up + model load
    load_s = time.time() - t0
    print(f"[crista] warm-up + model-load: {load_s:.2f}s")

    sizes = [int(x) for x in args.sizes.split(",") if x]
    results = []
    for n in sizes:
        sg, off, ctx = make_batch(n, seed=n)
        lat = []
        for _ in range(args.repeats):
            t0 = time.time()
            y = CRISTA_score.CRISTA_predict_list(sg, off, ctx)
            lat.append(time.time() - t0)
            assert len(y) == n
        best = min(lat)
        rec = {
            "n": n,
            "best_s": best,
            "ots_per_sec": n / best,
            "peak_rss_mb": round(peak_rss_mb(), 1),
        }
        results.append(rec)
        print(f"[crista] n={n:>8}  best={best:8.3f}s  {rec['ots_per_sec']:>12,.0f} OT/s  RSS={rec['peak_rss_mb']:.0f}MB")

    if results:
        biggest = max(results, key=lambda r: r["n"])
        gw = 6_770_000 / biggest["ots_per_sec"]
        print(f"[crista] projected GW (6.77M OT): {gw/60:.1f} min ({biggest['ots_per_sec']:,.0f} OT/s)")

    with open(args.out, "w") as fh:
        json.dump({"kind": "crista_speed", "model_load_s": load_s, "results": results}, fh, indent=2)


if __name__ == "__main__":
    main()
