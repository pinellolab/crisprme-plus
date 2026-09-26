#!/usr/bin/env python
"""
Persistent CRISPR-Bulge scoring worker — the modular-env boundary prototype.

Launched ONCE per run/chromosome inside its dedicated conda env, e.g.
    mamba run -n cbulge_cpu python scorer_worker.py
then fed newline-delimited JSON batches on stdin, one response line per request:

    request : {"pairs": [[sg_aligned, off_aligned], ...]}\n
    response: {"scores": [float, ...]}\n          (order preserved; invalid -> -1.0)
    "quit"  : {"cmd": "quit"}\n  -> worker exits 0

Models load once at startup (amortized across every batch). Scoring itself is
delegated to crispr_bulge_score.CRISPR_BULGE_predict_list so this worker and the
in-process path share ONE implementation.
"""
import json
import os
import sys
import time

REPO = os.environ.get("CBULGE_REPO", "/srv/local/lp698/cbulge_bench/CRISPR-Bulge")
if REPO not in sys.path:
    sys.path.insert(0, REPO)
# make crispr_bulge_score importable regardless of cwd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _log(msg):
    sys.stderr.write(f"[worker] {msg}\n")
    sys.stderr.flush()


def main():
    device = "gpu" if ("--gpu" in sys.argv or os.environ.get("CRISPRME_COMPUTE_BACKEND") == "gpu") else "cpu"

    # The scoring libs pollute stdout: build_sequence_features prints "The features
    # sizes are ..." and Keras prints "1/1 [====]" progress bars to fd 1, which would
    # corrupt the JSON line protocol. So dup fd 1 to a private protocol channel, then
    # redirect fd 1 -> fd 2 (stderr) so ALL library noise is harmless. Only clean JSON
    # goes to the real stdout the parent reads.
    proto = os.fdopen(os.dup(1), "w")
    os.dup2(2, 1)

    def send(obj):
        proto.write(json.dumps(obj) + "\n")
        proto.flush()

    import crispr_bulge_score as cb

    t0 = time.time()
    cb.load_models(device)
    _log(f"ready: 5 models loaded in {time.time()-t0:.2f}s (device={device})")
    send({"ready": True, "load_s": time.time() - t0})

    # readline() loop (NOT `for line in sys.stdin`, which read-ahead-buffers and
    # deadlocks a request/response protocol).
    while True:
        line = sys.stdin.readline()
        if not line:
            break  # EOF
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except Exception as e:
            send({"error": f"bad json: {e}"})
            continue
        if req.get("cmd") == "quit":
            _log("quit")
            return
        pairs = req.get("pairs", [])
        if pairs:
            sg = [p[0] for p in pairs]
            off = [p[1] for p in pairs]
            scores = cb.CRISPR_BULGE_predict_list(sg, off, device=device)
        else:
            scores = []
        send({"scores": scores})


if __name__ == "__main__":
    main()
