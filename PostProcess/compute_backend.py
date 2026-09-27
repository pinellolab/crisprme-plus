"""Compute-backend selection for the optional GPU acceleration.

CRISPRme+ runs CPU by default; GPU is opt-in and, crucially, PORTABLE across the two
GPU stacks scientists actually have:

  * ``cuda``  — NVIDIA GPUs on Linux (clusters, workstations): TensorFlow CUDA build
                for the CRISPR-Bulge scorer; cupy for a vectorized CFD.
  * ``metal`` — Apple-Silicon Macs (M1/M2/M3): TensorFlow via the tensorflow-metal
                PluggableDevice for the scorer. NOTE: cupy is CUDA-only, so there is no
                Metal cupy — CFD stays on CPU (it is a tiny lookup + per-position product,
                already fast; the GPU win on Mac is the neural scorer, not CFD).
  * ``cpu``   — the mandatory, always-available fallback.

One env var drives it: ``CRISPRME_COMPUTE_BACKEND`` = ``cpu`` (default) | ``gpu``/``auto``
(auto-pick cuda→metal→cpu) | ``cuda`` | ``metal``. A requested accelerator that isn't
present degrades to CPU with a single warning (never hard-fails) — CPU-mandatory.

Dependency-free at import (stdlib only); heavy libs (tensorflow/cupy) are imported lazily
by the callers, guided by resolve_backend().
"""

import os
import platform
import shutil
import sys

CPU, CUDA, METAL = "cpu", "cuda", "metal"


def _is_apple_silicon():
    return platform.system() == "Darwin" and platform.machine() == "arm64"


def _has_cuda():
    # lightweight proxy: the NVIDIA runtime tool is on PATH, or CUDA devices are pinned.
    if shutil.which("nvidia-smi"):
        return True
    v = os.environ.get("CUDA_VISIBLE_DEVICES")
    return bool(v) and v not in ("", "-1")


def _warn(msg):
    sys.stderr.write(f"[compute-backend] {msg}\n")


def resolve_backend(requested=None):
    """Return the EFFECTIVE backend ('cpu'|'cuda'|'metal') after detection + fallback.

    ``requested`` (or $CRISPRME_COMPUTE_BACKEND) is one of cpu | gpu/auto | cuda | metal.
    A requested accelerator that isn't available degrades to CPU with a warning.
    """
    req = (requested or os.environ.get("CRISPRME_COMPUTE_BACKEND", "cpu") or "cpu").lower()
    if req in ("cpu", ""):
        return CPU
    if req == CUDA:
        if _has_cuda():
            return CUDA
        _warn("cuda requested but no NVIDIA GPU detected; using CPU")
        return CPU
    if req == METAL:
        if _is_apple_silicon():
            return METAL
        _warn("metal requested but this is not an Apple-Silicon Mac; using CPU")
        return CPU
    if req in ("gpu", "auto"):
        if _has_cuda():
            return CUDA
        if _is_apple_silicon():
            return METAL
        _warn("gpu requested but no CUDA/Metal device detected; using CPU")
        return CPU
    _warn(f"unknown CRISPRME_COMPUTE_BACKEND={req!r}; using CPU")
    return CPU


def pre_import_env(backend):
    """Set env that must be in place BEFORE `import tensorflow`. For CPU we hide any
    NVIDIA GPU via CUDA_VISIBLE_DEVICES='' (harmless on Mac). cuda/metal leave it alone
    (Metal is selected post-import via tf.config; there is no CUDA env for it)."""
    if backend == CPU:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""


def _cpu_thread_cap():
    """How many intra-op threads the scorer should use on CPU. TF's default grabs ALL
    logical cores, which on a big shared node (e.g. a 256-core cluster box) OVERSUBSCRIBES
    and runs ~40% slower than a modest cap (measured: 6.0k vs 8.3k OT/s at cap=16). Honor
    an explicit override, else cap at 16."""
    for var in ("CRISPRME_SCORER_THREADS", "OMP_NUM_THREADS"):
        v = os.environ.get(var, "")
        if v.isdigit() and int(v) > 0:
            return int(v)
    try:
        return min(16, os.cpu_count() or 8)
    except Exception:
        return 8


def configure_tf(tf, backend):
    """Apply thread + device configuration AFTER `import tensorflow`. Returns the backend
    actually in effect ('cpu' if the requested accelerator turned out to have no visible
    GPU device).

    On Metal, tensorflow-metal exposes the Apple GPU as a 'GPU' physical device, so the
    same code path (enable memory growth, keep visible) drives both cuda and metal; for
    cpu we hide all GPU devices so TF stays on the CPU even if a device is present.
    """
    # Cap CPU thread parallelism first (must precede any op execution). Harmless under
    # cuda/metal -- it only bounds CPU-side ops -- and crucial for the CPU path + the
    # Metal->CPU numerical fallback, which otherwise oversubscribe a many-core host.
    try:
        tf.config.threading.set_intra_op_parallelism_threads(_cpu_thread_cap())
        tf.config.threading.set_inter_op_parallelism_threads(2)
    except Exception:
        pass
    try:
        gpus = tf.config.list_physical_devices("GPU")
    except Exception:
        gpus = []
    if backend == CPU:
        if gpus:
            try:
                tf.config.set_visible_devices([], "GPU")
            except Exception:
                pass
        return CPU
    if not gpus:
        _warn(f"{backend} requested but TensorFlow sees no GPU device; using CPU")
        return CPU
    for g in gpus:  # memory growth is a no-op / best-effort on unified-memory Metal
        try:
            tf.config.experimental.set_memory_growth(g, True)
        except Exception:
            pass
    return backend


def array_module(backend):
    """Return the array module for a vectorized CFD path: cupy under CUDA, else numpy.
    Metal has no cupy, so CFD runs on numpy (CPU) there -- it is not the bottleneck."""
    if backend == CUDA:
        try:
            import cupy as cp  # noqa
            return cp
        except Exception:
            _warn("cuda backend but cupy unavailable; CFD uses numpy (CPU)")
    import numpy as np
    return np
