"""Dedicated conda environments for the ML off-target scorers (modular scorer envs).

Heavy/vendored-model scorers (CRISPR-Bulge = TensorFlow) run in their OWN conda
env so their pinned deps never co-solve with CRISPRme's main stack. This module is
the lifecycle manager for those envs, used by:

  * the installer (Dockerfile / install_from_source.sh)  -> create the env
  * ``crisprme.py scorer-env {create,check,update,list,doctor}``  -> manage/diagnose
  * ``scorer_runner``  -> locate the env's python to launch a persistent worker

Design notes (all learned in the Phase-0 benchmark):
  * The env manager DIFFERS by install path: **micromamba** inside the Docker image
    (``FROM mambaorg/micromamba``), **mamba/conda** for source installs. We detect
    whichever is available rather than hard-coding one.
  * The CRISPR-Bulge encoder import chain pulls xgboost + catboost + pkg_resources
    even for the NN-only path, and setuptools>=81 dropped pkg_resources -> pin
    ``setuptools<81``. numpy is pinned to 1.23.5 (CRISPR-Bulge requirement).
  * State is persisted to ``<data>/Annotations/.scorer_env.json`` (atomic tmp +
    os.replace), mirroring ``cosmic_license.py``. Reads never raise.

Dependency-free (stdlib only) so both the CLI and the Dash web layer can import it.
Actual scoring happens in ``crispr_bulge_score`` / ``scorer_worker`` inside the env.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Environment specifications
# ---------------------------------------------------------------------------
# Versions here are the ones validated end-to-end in Phase 0 (fidelity bit-identical
# to upstream on tf 2.13 CPU; A100 speed on tf 2.12 CUDA). numpy=1.23.5 + setuptools<81
# are hard requirements; xgboost/catboost are pulled by the encoder import chain.
_COMMON = ["numpy=1.23.5", "pandas", "scikit-learn", "xgboost", "catboost", "setuptools<81"]

SCORER_ENVS: Dict[str, dict] = {
    "cbulge": {
        "description": "CRISPR-Bulge off-target scorer (TensorFlow GRU ensemble)",
        "channels": ["conda-forge"],
        "python": "3.10",
        # CPU is mandatory; the GPU variant is opt-in and swaps in the CUDA TF build.
        "cpu_packages": ["tensorflow-cpu=2.13"] + _COMMON,
        "gpu_packages": ['tensorflow=2.12=cuda*'] + _COMMON,
        # import probe run INSIDE the env to confirm the scorer can load
        "probe_imports": [
            "tensorflow", "numpy", "pandas", "sklearn", "xgboost", "catboost", "pkg_resources",
        ],
    },
}

DEFAULT_ENV = "cbulge"

# ---------------------------------------------------------------------------
# Env-manager detection
# ---------------------------------------------------------------------------


def detect_env_manager() -> Optional[Tuple[str, str]]:
    """Return (absolute_exe, kind) for the available conda env manager, or None.

    Preference order: an explicit override, then micromamba (Docker image), then
    mamba, then conda. ``kind`` is one of 'micromamba' | 'mamba' | 'conda'.
    """
    override = os.environ.get("CRISPRME_CONDA_EXE")
    candidates: List[Tuple[Optional[str], str]] = []
    if override:
        candidates.append((override, _kind_of(override)))
    # env vars set by the respective installers (kind derived from the actual binary:
    # miniforge sets MAMBA_EXE to a 'mamba', the Docker image to 'micromamba')
    mamba_exe = os.environ.get("MAMBA_EXE")
    if mamba_exe:
        candidates.append((mamba_exe, _kind_of(mamba_exe)))
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        candidates.append((conda_exe, _kind_of(conda_exe)))
    for name in ("micromamba", "mamba", "conda"):
        candidates.append((shutil.which(name), _kind_of(name)))
    for exe, kind in candidates:
        if exe and os.path.exists(exe) and os.access(exe, os.X_OK):
            return os.path.abspath(exe), kind
        if exe and shutil.which(exe):
            return shutil.which(exe), kind
    return None


def _kind_of(path_or_name: str) -> str:
    base = os.path.basename(path_or_name)
    if "micromamba" in base:
        return "micromamba"
    if "mamba" in base:
        return "mamba"
    return "conda"


# ---------------------------------------------------------------------------
# Env location / existence
# ---------------------------------------------------------------------------


def _run(cmd: List[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def env_python(name: str) -> Optional[str]:
    """Absolute path to the env's python, or None if the env is absent.

    Resolved via ``<manager> run -n <name> python -c ...`` once (authoritative,
    location-independent). Callers should invoke the worker via THIS absolute path,
    NOT ``<manager> run`` — the latter captures/buffers stdio and breaks the
    streaming worker protocol (Phase-0 finding).
    """
    mgr = detect_env_manager()
    if not mgr:
        return None
    exe, _ = mgr
    try:
        cp = _run([exe, "run", "-n", name, "python", "-c", "import sys;print(sys.executable)"])
        if cp.returncode == 0:
            path = cp.stdout.strip().splitlines()[-1].strip() if cp.stdout.strip() else ""
            if path and os.path.exists(path):
                return path
    except Exception:
        pass
    # fallback: derive from the manager root prefix
    for prefix in _candidate_roots(exe):
        cand = os.path.join(prefix, "envs", name, "bin", "python")
        if os.path.exists(cand):
            return cand
    return None


def _candidate_roots(exe: str) -> List[str]:
    roots = []
    for var in ("MAMBA_ROOT_PREFIX", "CONDA_ROOT", "CONDA_PREFIX"):
        v = os.environ.get(var)
        if v:
            roots.append(v)
    # exe is typically <root>/bin/<manager>
    roots.append(os.path.dirname(os.path.dirname(os.path.abspath(exe))))
    return roots


def env_exists(name: str) -> bool:
    return env_python(name) is not None


# ---------------------------------------------------------------------------
# Create / update
# ---------------------------------------------------------------------------


def _packages(spec: dict, gpu: bool) -> List[str]:
    return spec["gpu_packages"] if gpu else spec["cpu_packages"]


def build_create_command(name: str, gpu: bool = False) -> Optional[List[str]]:
    """Return the argv to create the env, or None if no manager is available."""
    spec = SCORER_ENVS[name]
    mgr = detect_env_manager()
    if not mgr:
        return None
    exe, _ = mgr
    channels = []
    for ch in spec["channels"]:
        channels += ["-c", ch]
    return (
        [exe, "create", "-y", "-n", name, "python=" + spec["python"]]
        + channels
        + _packages(spec, gpu)
    )


def create_env(name: str = DEFAULT_ENV, gpu: bool = False, force: bool = False,
               stream: bool = True) -> Tuple[bool, str]:
    """Create the scorer env (idempotent). Returns (ok, message).

    If the env already exists and ``force`` is False, this is a no-op success.
    GPU builds set CONDA_OVERRIDE_CUDA so the CUDA TF build resolves off-GPU hosts.
    """
    if name not in SCORER_ENVS:
        return False, f"unknown scorer env '{name}'"
    if env_exists(name) and not force:
        return True, f"env '{name}' already exists (use --force to recreate)"
    cmd = build_create_command(name, gpu=gpu)
    if cmd is None:
        return False, ("no conda env manager found (need micromamba/mamba/conda on "
                       "PATH, or set CRISPRME_CONDA_EXE)")
    if force and env_exists(name):
        exe = cmd[0]
        _run([exe, "env", "remove", "-y", "-n", name])
    env = dict(os.environ)
    if gpu:
        env.setdefault("CONDA_OVERRIDE_CUDA", "12.0")
    sys.stderr.write(f"[scorer-env] creating '{name}' ({'gpu' if gpu else 'cpu'}): {' '.join(cmd)}\n")
    if stream:
        proc = subprocess.run(cmd, env=env)
        ok = proc.returncode == 0
        return ok, ("created" if ok else f"create failed (exit {proc.returncode})")
    cp = _run(cmd, env=env)
    return cp.returncode == 0, (cp.stdout + cp.stderr)[-4000:]


def update_env(name: str = DEFAULT_ENV, gpu: bool = False) -> Tuple[bool, str]:
    """Repair/update: create if missing, else re-assert the package set (idempotent)."""
    if not env_exists(name):
        return create_env(name, gpu=gpu)
    spec = SCORER_ENVS[name]
    mgr = detect_env_manager()
    if not mgr:
        return False, "no conda env manager found"
    exe, _ = mgr
    channels = []
    for ch in spec["channels"]:
        channels += ["-c", ch]
    cmd = [exe, "install", "-y", "-n", name] + channels + _packages(spec, gpu)
    env = dict(os.environ)
    if gpu:
        env.setdefault("CONDA_OVERRIDE_CUDA", "12.0")
    sys.stderr.write(f"[scorer-env] updating '{name}': {' '.join(cmd)}\n")
    proc = subprocess.run(cmd, env=env)
    ok = proc.returncode == 0
    return ok, ("updated" if ok else f"update failed (exit {proc.returncode})")


# ---------------------------------------------------------------------------
# Health check
# ---------------------------------------------------------------------------

OK, WARN, ERROR = "ok", "warn", "error"


def health_check(name: str = DEFAULT_ENV) -> dict:
    """Probe env health. Returns a structured record (never raises):

        {status: 'ok'|'warn'|'error', issues: [(level, msg), ...],
         python_version, versions: {pkg: ver}, weights_present: bool, ...}

    'error' = env or a required import is missing (scorer unusable).
    'warn'  = env is fine but weights/source aren't staged yet (P4 fetches them).
    """
    spec = SCORER_ENVS.get(name)
    rec: dict = {"env": name, "status": OK, "issues": [], "versions": {},
                 "python_version": None, "weights_present": False}
    if spec is None:
        rec["status"] = ERROR
        rec["issues"].append((ERROR, f"unknown scorer env '{name}' "
                                     f"(known: {', '.join(sorted(SCORER_ENVS))})"))
        return rec

    mgr = detect_env_manager()
    if not mgr:
        rec["status"] = ERROR
        rec["issues"].append((ERROR, "no conda env manager (micromamba/mamba/conda) found"))
        return rec
    rec["manager"] = {"exe": mgr[0], "kind": mgr[1]}

    py = env_python(name)
    if not py:
        rec["status"] = ERROR
        rec["issues"].append((ERROR, f"env '{name}' does not exist "
                                     f"(create with: crisprme.py scorer-env create)"))
        return rec
    rec["env_python"] = py

    # import probe INSIDE the env (module list injected as a literal; no %-format
    # collisions with the inner code)
    probe = (
        "import json,sys\n"
        "mods=" + repr(list(spec["probe_imports"])) + "\n"
        "vers={}\n"
        "bad=[]\n"
        "for m in mods:\n"
        "    try:\n"
        "        mod=__import__(m)\n"
        "        vers[m]=getattr(mod,'__version__','?')\n"
        "    except Exception as e:\n"
        "        bad.append([m,str(e)])\n"
        "pv='.'.join(str(x) for x in sys.version_info[:3])\n"
        "print(json.dumps({'py':pv,'vers':vers,'bad':bad}))\n"
    )
    try:
        cp = _run([py, "-c", probe])
        payload = json.loads(cp.stdout.strip().splitlines()[-1]) if cp.stdout.strip() else {}
    except Exception as e:
        rec["status"] = ERROR
        rec["issues"].append((ERROR, f"import probe failed to run: {e}"))
        return rec

    rec["python_version"] = payload.get("py")
    rec["versions"] = payload.get("vers", {})
    for mod, err in payload.get("bad", []):
        rec["status"] = ERROR
        rec["issues"].append((ERROR, f"import '{mod}' failed in env: {err.splitlines()[0][:160]}"))

    # weights / source presence (P4 will fetch from HF; for now via CBULGE_REPO)
    repo = os.environ.get("CBULGE_REPO", "")
    if repo and os.path.isdir(os.path.join(repo, "OT_deep_score_src")):
        rec["weights_present"] = True
        rec["cbulge_repo"] = repo
    else:
        if rec["status"] != ERROR:
            rec["status"] = WARN
        rec["issues"].append((WARN, "CRISPR-Bulge source/weights not found "
                                    "(set CBULGE_REPO; P4 will fetch from HuggingFace on first use)"))
    return rec


# ---------------------------------------------------------------------------
# Persisted state  (mirrors cosmic_license.py)
# ---------------------------------------------------------------------------

SCORER_ENV_FILE = ".scorer_env.json"


def state_path(annotations_dir: str) -> str:
    return os.path.join(annotations_dir, SCORER_ENV_FILE)


def get_scorer_env_state(annotations_dir: str) -> dict:
    try:
        with open(state_path(annotations_dir)) as fh:
            data = json.load(fh)
            return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def set_scorer_env_state(annotations_dir: str, record: dict) -> bool:
    """Persist a health/state record atomically. Best-effort: this is a diagnostic
    cache, so a write failure (read-only dir, concurrent race) must never crash the
    caller. Uses a PROCESS-UNIQUE temp file so concurrent writers on a shared
    (e.g. NFS) Annotations dir don't clobber each other's temp. Returns True on success.
    """
    payload = {"version": 1}
    payload.update(record)
    try:
        os.makedirs(annotations_dir, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=annotations_dir, prefix=".scorer_env.", suffix=".tmp")
        try:
            with os.fdopen(fd, "w") as fh:
                json.dump(payload, fh)
            os.replace(tmp, state_path(annotations_dir))
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
        return True
    except OSError:
        return False


def render_health(rec: dict) -> str:
    """Human-friendly one-block summary of a health record."""
    lines = [f"Scorer env '{rec.get('env', '?')}': {rec.get('status', '?').upper()}"]
    if rec.get("manager"):
        lines.append(f"  manager: {rec['manager']['kind']} ({rec['manager']['exe']})")
    if rec.get("env_python"):
        lines.append(f"  python:  {rec.get('python_version')}  ({rec['env_python']})")
    if rec.get("versions"):
        key = ", ".join(f"{k}={v}" for k, v in rec["versions"].items())
        lines.append(f"  packages: {key}")
    for level, msg in rec.get("issues", []):
        lines.append(f"  [{level.upper()}] {msg}")
    return "\n".join(lines)
