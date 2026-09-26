"""Modular scorer-runner — dispatch a batch to a scorer, in-process or in its own env.

The main CRISPRme process carries no TensorFlow. To score with CRISPR-Bulge it
launches ONE persistent worker (``scorer_worker.py``) inside the dedicated ``cbulge``
conda env and streams batches to it over a newline-JSON protocol:

    request : {"pairs": [[sg_aligned, off_aligned], ...]}\\n
    response: {"scores": [float, ...]}\\n      (order preserved; invalid row -> -1.0)

Contract (uniform, mirrors the CRISTA/CFD list API):

    with ScorerRunner("cbulge", device="cpu") as r:
        scores = r.predict(sg_list, off_list)      # list[float], len == len(sg_list)

Design points (Phase-0 findings baked in):
  * The worker is spawned via the env's python ABSOLUTE PATH (never ``mamba run`` —
    that captures/buffers stdio and breaks the streaming protocol).
  * ONE worker per run/chromosome amortizes the ~4.5 s model spawn+load.
  * GRACEFUL DEGRADATION: if the env/worker is missing or dies, the runner disables
    itself, warns ONCE, and returns -1.0 sentinels — a scoring run is never killed by
    a missing optional scorer (the caller can fall back to CRISTA/CFD).

Depends only on the stdlib + ``scorer_env`` (also stdlib-only).
"""

import json
import math
import os
import select
import subprocess
import sys
import time

import scorer_env


def _pos_float(env_name, default):
    """Parse a positive float env var; fall back to default on absent/<=0/garbage."""
    try:
        v = float(os.environ.get(env_name, default))
        return v if v > 0 else float(default)
    except (TypeError, ValueError):
        return float(default)


# how long to wait for one batch response before declaring the worker hung
_BATCH_TIMEOUT_S = _pos_float("CRISPRME_SCORER_TIMEOUT", "1800")
# how long to wait for the startup handshake (model spawn+load)
_HANDSHAKE_TIMEOUT_S = _pos_float("CRISPRME_SCORER_LOAD_TIMEOUT", "600")
# split very large requests so neither side buffers an unbounded JSON blob
_MAX_BATCH = max(1, int(_pos_float("CRISPRME_SCORER_MAX_BATCH", "100000")))
# reject a single protocol line larger than this (a healthy worker sends ~tens of MB
# for a 100k batch; anything far bigger means a broken/runaway worker)
_MAX_LINE_BYTES = int(_pos_float("CRISPRME_SCORER_MAX_LINE_MB", "512")) * 1024 * 1024


class ScorerRunner:
    def __init__(self, scorer=scorer_env.DEFAULT_ENV, device=None, worker=None):
        self.scorer = scorer
        self.device = (device or os.environ.get("CRISPRME_COMPUTE_BACKEND", "cpu")).lower()
        self.worker = worker or os.path.join(os.path.dirname(os.path.abspath(__file__)), "scorer_worker.py")
        self.proc = None
        self.disabled = False
        self._reason = ""
        self._warned = False
        self.load_s = None

    # -- lifecycle ---------------------------------------------------------
    def _disable(self, reason):
        self.disabled = True
        self._reason = reason
        if not self._warned:
            sys.stderr.write(f"[scorer-runner] {self.scorer} disabled: {reason}. "
                             f"Falling back (scores = -1).\n")
            self._warned = True

    def _ensure_worker(self):
        if self.proc is not None or self.disabled:
            return
        py = scorer_env.env_python(self.scorer)
        if not py:
            self._disable(f"env '{self.scorer}' not found "
                          f"(create with: crisprme.py scorer-env create)")
            return
        env = dict(os.environ)
        env["CRISPRME_COMPUTE_BACKEND"] = self.device
        env.setdefault("PYTHONUNBUFFERED", "1")
        # point the worker at the provisioned CRISPR-Bulge source+weights
        env.setdefault("CBULGE_REPO", scorer_env.default_cbulge_repo())
        cmd = [py, "-u", self.worker] + (["--gpu"] if self.device == "gpu" else [])
        try:
            self.proc = subprocess.Popen(
                cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr,
                text=True, bufsize=1, env=env,
            )
        except Exception as e:
            self._disable(f"failed to spawn worker: {e}")
            return
        # readiness handshake (worker sends one JSON line when models are loaded)
        line = self._readline(_HANDSHAKE_TIMEOUT_S)
        if not line:
            self._disable("worker did not report ready (spawn/load timeout or crash)")
            self._kill()
            return
        try:
            ready = json.loads(line)
        except Exception:
            self._disable(f"bad handshake from worker: {line[:120]!r}")
            self._kill()
            return
        if not ready.get("ready"):
            self._disable(f"worker reported not-ready: {ready.get('error', 'unknown')}")
            self._kill()
            return
        self.load_s = ready.get("load_s")

    def _readline(self, timeout_s):
        """Blocking readline with a timeout; returns '' on timeout/EOF."""
        if self.proc is None or self.proc.stdout is None:
            return ""
        deadline = time.monotonic() + timeout_s
        fd = self.proc.stdout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return ""
            r, _, _ = select.select([fd], [], [], remaining)
            if not r:
                return ""
            line = fd.readline()
            if not line:
                return ""  # EOF: worker closed stdout / exited (don't busy-loop)
            if len(line) > _MAX_LINE_BYTES:
                return ""  # runaway/broken worker line -> treat as failure
            return line

    def _kill(self):
        if self.proc is not None:
            proc = self.proc
            try:
                proc.kill()
            except Exception:
                pass
            # close the pipes so their fds aren't leaked (stderr is our sys.stderr)
            for pipe in (proc.stdin, proc.stdout):
                try:
                    if pipe is not None:
                        pipe.close()
                except Exception:
                    pass
            self.proc = None

    def close(self):
        if self.proc is not None:
            proc = self.proc
            try:
                proc.stdin.write(json.dumps({"cmd": "quit"}) + "\n")
                proc.stdin.flush()
                proc.wait(timeout=15)
            except Exception:
                try:
                    proc.kill()
                except Exception:
                    pass
            finally:
                # explicitly close the pipes so their fds are not leaked if the
                # worker died mid-protocol (stderr is our own sys.stderr; leave it)
                for pipe in (proc.stdin, proc.stdout):
                    try:
                        if pipe is not None:
                            pipe.close()
                    except Exception:
                        pass
                self.proc = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # -- scoring -----------------------------------------------------------
    def predict(self, sg_list, off_list):
        """Score aligned (sgRNA, off-target) pairs -> list[float] (len == inputs).

        Never raises: on any failure the runner disables and returns -1.0 sentinels
        so the caller can fall back to another scorer.
        """
        n = len(sg_list)
        if n != len(off_list):
            raise ValueError("sg_list and off_list differ in length")
        if n == 0:
            return []
        # split oversized requests so neither side buffers an unbounded JSON blob;
        # a disabled runner makes later chunks return sentinels (sticky), no respawn
        if n > _MAX_BATCH:
            out = []
            for i in range(0, n, _MAX_BATCH):
                out.extend(self.predict(sg_list[i:i + _MAX_BATCH], off_list[i:i + _MAX_BATCH]))
            return out

        self._ensure_worker()
        if self.disabled or self.proc is None:
            return [-1.0] * n
        req = json.dumps({"pairs": [[sg_list[i], off_list[i]] for i in range(n)]})
        try:
            self.proc.stdin.write(req + "\n")
            self.proc.stdin.flush()
        except Exception as e:
            self._disable(f"worker write failed: {e}")
            self._kill()
            return [-1.0] * n
        line = self._readline(_BATCH_TIMEOUT_S)
        if not line:
            self._disable("worker timed out or crashed during scoring")
            self._kill()
            return [-1.0] * n
        # validate the whole response defensively: a malformed/garbage reply must
        # degrade to sentinels, never raise (contract: predict() never throws)
        try:
            scores = json.loads(line).get("scores", [])
            if not isinstance(scores, list):
                raise ValueError(f"'scores' is {type(scores).__name__}, not a list")
            if len(scores) != n:
                self._disable(f"worker returned {len(scores)} scores for {n} inputs")
                self._kill()
                return [-1.0] * n
            out = []
            for s in scores:
                f = float(s)  # raises on str/None/bool-ish garbage
                # -1.0 is the worker's own per-row failure sentinel; keep it. Any
                # other out-of-[0,1] / NaN / Inf value is corrupt -> fail the batch.
                if f != -1.0 and (math.isnan(f) or math.isinf(f) or f < 0.0 or f > 1.0):
                    raise ValueError(f"out-of-range score {f!r}")
                out.append(f)
            return out
        except Exception as e:
            self._disable(f"bad response from worker: {e}")
            self._kill()
            return [-1.0] * n


# ---------------------------------------------------------------------------
# Module-level singleton convenience (what the pipeline seam will call in P4)
# ---------------------------------------------------------------------------
_SINGLETON = None


def get_runner(scorer=scorer_env.DEFAULT_ENV, device=None):
    """Return a process-wide singleton runner (one persistent worker per run)."""
    global _SINGLETON
    if _SINGLETON is None:
        _SINGLETON = ScorerRunner(scorer=scorer, device=device)
    return _SINGLETON


def CRISPR_BULGE_predict_list(sg_list, off_list, device=None):
    """Drop-in batch scorer mirroring CRISTA_predict_list's shape (no 29-nt ctx),
    routed through the persistent env-worker singleton. Graceful -1.0 on failure."""
    return get_runner("cbulge", device=device).predict(sg_list, off_list)


def close_runner():
    global _SINGLETON
    if _SINGLETON is not None:
        _SINGLETON.close()
        _SINGLETON = None
