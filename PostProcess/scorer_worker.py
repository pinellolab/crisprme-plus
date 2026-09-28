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

# make crispr_bulge_score importable regardless of cwd. NOTE: we intentionally do
# NOT put CBULGE_REPO on sys.path here — crispr_bulge_score.load_models() does that
# via a validated path (single, guarded injection site).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _log(msg):
    sys.stderr.write(f"[worker] {msg}\n")
    sys.stderr.flush()


def main():
    # Quiet TensorFlow's own C++ logs BEFORE it is imported (crispr_bulge_score imports it).
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")   # errors only, no INFO/WARNING spam

    # Pass the requested backend through verbatim; compute_backend.resolve_backend()
    # (called inside load_models) does detection + fallback for cpu|gpu|cuda|metal|auto.
    # --gpu is a legacy alias kept for back-compat with older spawn sites.
    device = os.environ.get("CRISPRME_COMPUTE_BACKEND") or ("gpu" if "--gpu" in sys.argv else "cpu")

    # The scoring libs are NOISY on BOTH stdout and stderr: build_sequence_features prints
    # "The features sizes are ...", Keras prints "1/1 [====]" progress bars to fd 1, and
    # TF/Keras emit retracing WARNINGs + C++ logs to fd 2. That output is harmless in
    # isolation, but the CRISPRme post-analysis stage treats ANY bytes on the subprocess's
    # stderr as fatal (`[ -s $logerror ]`), so routing this noise to stderr would fail every
    # real search's post-analysis (CRISTA was quiet; TF is not). So: dup fd 1 to a PRIVATE
    # protocol channel, then send BOTH fd 1 and fd 2 to /dev/null. The worker communicates
    # exclusively over the proto channel (ready/error + score responses), so silencing the
    # real fds loses nothing the parent needs. Set CRISPRME_SCORER_DEBUG=1 to keep the
    # library noise on the real stderr for troubleshooting.
    proto = os.fdopen(os.dup(1), "w")
    if os.environ.get("CRISPRME_SCORER_DEBUG"):
        os.dup2(2, 1)  # legacy: library stdout -> stderr (visible, but fatal to post-analysis)
    else:
        _devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(_devnull, 1)   # library stdout noise -> /dev/null
        os.dup2(_devnull, 2)   # library stderr noise (TF logs, retracing warnings) -> /dev/null

    def send(obj):
        proto.write(json.dumps(obj) + "\n")
        proto.flush()

    import crispr_bulge_score as cb

    t0 = time.time()
    try:
        cb.load_models(device)
    except Exception as e:
        # fail-fast: tell the parent we're NOT ready so it disables cleanly instead
        # of waiting out the handshake timeout
        _log(f"model load failed: {e}")
        send({"ready": False, "error": str(e)[:300]})
        return
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
