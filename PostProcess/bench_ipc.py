#!/usr/bin/env python
"""
Phase-0 IPC / subprocess-boundary overhead probe for the modular-env design.

Spawns ONE persistent scorer_worker.py in its dedicated env via `mamba run`,
feeds it batches over stdin/stdout JSON, and reports:
  * one-time worker spawn + model-load latency (amortized per run/chromosome)
  * per-batch round-trip time at {1k,10k,100k}
  * the round-trip MINUS pure compute (the serialization + pipe cost) so we can
    confirm the boundary is negligible vs GRU inference at 100k-row batches.

Runs in the LIGHT parent env (only stdlib) — that's the whole point: the main
CRISPRme process never imports TensorFlow.

    python bench_ipc.py --env cbulge_cpu --sizes 1000,10000,100000
"""
import argparse
import json
import os
import subprocess
import sys
import time

REPO = os.environ.get("CBULGE_REPO", "/srv/local/lp698/cbulge_bench/CRISPR-Bulge")
MAMBA = os.environ.get("MAMBA_BIN", "/data/pinello/SHARED_SOFTWARE/miniforge3/bin/mamba")


def make_pairs(n, seed=0):
    import random

    rng = random.Random(seed)
    bases = "ACGT"
    guide = "".join(rng.choice(bases) for _ in range(20)) + "AGG"
    pairs = []
    for _ in range(n):
        off = list(guide)
        for _ in range(rng.randint(0, 5)):
            off[rng.randrange(20)] = rng.choice(bases)
        pairs.append([guide, "".join(off)])
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--env", default="cbulge_cpu")
    ap.add_argument("--sizes", default="1000,10000,100000")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--out", default="bench_ipc_out.json")
    args = ap.parse_args()

    worker = os.path.join(REPO, "scorer_worker.py")
    # Resolve the env's python ONCE and invoke it directly. `mamba run` captures/
    # buffers stdio by default, which breaks a streaming stdin/stdout protocol
    # (handshake works, first batch deadlocks/broken-pipe). Invoking the env python
    # by absolute path is both the fix and what the production runner should do.
    env_python = subprocess.check_output(
        [MAMBA, "run", "-n", args.env, "python", "-c", "import sys;print(sys.executable)"],
        text=True,
    ).strip().splitlines()[-1]
    print(f"[ipc] env python: {env_python}")
    cmd = [env_python, "-u", worker]
    env = dict(os.environ, CBULGE_REPO=REPO, CUDA_VISIBLE_DEVICES="")

    t_spawn = time.time()
    proc = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr,
        text=True, bufsize=1, env=env,
    )
    # first line = readiness handshake (includes model-load time inside the worker)
    ready_line = proc.stdout.readline()
    spawn_plus_load = time.time() - t_spawn
    ready = json.loads(ready_line)
    print(f"[ipc] spawn+load handshake: {spawn_plus_load:.2f}s  (worker load_s={ready.get('load_s'):.2f}s)")

    results = []
    for n in [int(x) for x in args.sizes.split(",") if x]:
        pairs = make_pairs(n, seed=n)
        payload = json.dumps({"pairs": pairs})
        rtts = []
        for _ in range(args.repeats):
            t0 = time.time()
            proc.stdin.write(payload + "\n")
            proc.stdin.flush()
            resp = proc.stdout.readline()
            rtt = time.time() - t0
            obj = json.loads(resp)
            assert len(obj["scores"]) == n, (len(obj["scores"]), n)
            rtts.append(rtt)
        best = min(rtts)
        # serialize-only cost (encode+decode the same payload, no IPC/compute)
        t0 = time.time()
        _ = json.loads(json.dumps({"scores": [0.0] * n}))
        _ = json.dumps({"pairs": pairs})
        ser_s = time.time() - t0
        rec = {"n": n, "rtt_best_s": best, "ots_per_sec": n / best, "serialize_s": ser_s}
        results.append(rec)
        print(f"[ipc] n={n:>8}  rtt={best:8.3f}s  {rec['ots_per_sec']:>12,.0f} OT/s  serialize~{ser_s*1000:.1f}ms")

    proc.stdin.write(json.dumps({"cmd": "quit"}) + "\n")
    proc.stdin.flush()
    proc.wait(timeout=30)

    with open(args.out, "w") as fh:
        json.dump({"kind": "ipc", "spawn_plus_load_s": spawn_plus_load,
                   "worker_load_s": ready.get("load_s"), "results": results}, fh, indent=2)
    print(f"[ipc] wrote {args.out}")


if __name__ == "__main__":
    main()
