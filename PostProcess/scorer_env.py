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
import platform
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
        # Keep CPU and GPU on the SAME TF version (2.13) so the model behaves identically
        # on both — the conda-forge cuda build (tensorflow=2.13=cuda*, e.g. cuda118) is the
        # exact stack validated correct + fast on an A100 (14.7k OT/s, GRU numerically correct;
        # unlike tensorflow-metal it does NOT miscompute). It also runs on CPU-only hosts.
        "cpu_packages": ["tensorflow-cpu=2.13"] + _COMMON,
        "gpu_packages": ['tensorflow=2.13=cuda*'] + _COMMON,
        # Apple-Silicon Metal variant: conda-forge has no Metal TensorFlow, so the base
        # deps come from conda (NO tensorflow) and the Metal TF comes from pip
        # (tensorflow-macos + the tensorflow-metal PluggableDevice). hdf5 for h5py.
        "metal_packages": ["hdf5"] + _COMMON,
        "metal_pip": ["tensorflow-macos==2.13.0", "tensorflow-metal==1.0.1"],
        # Linux aarch64 variant: conda-forge ships no tensorflow-cpu=2.13 for aarch64
        # (only 2.18/2.19), but PyPI publishes the official manylinux2014_aarch64 wheel of
        # the full tensorflow==2.13.1. So the base deps come from conda (NO tensorflow) and
        # the SAME TF version used on x86-64 (2.13.1) is pip-installed -- keeping the model
        # numerically identical across arches. The full 'tensorflow' package runs CPU-only on
        # a GPU-less ARM host (benign "no CUDA" warnings); there is no CUDA/Metal on aarch64.
        "linux_aarch64_packages": ["hdf5"] + _COMMON,
        "linux_aarch64_pip": ["tensorflow==2.13.1"],
        # import probe run INSIDE the env to confirm the scorer can load
        "probe_imports": [
            "tensorflow", "numpy", "pandas", "sklearn", "xgboost", "catboost", "pkg_resources",
        ],
    },
}

# The scorer conda env name. Overridable via CRISPRME_SCORER_ENV so a site can point
# CRISPRme at an existing / shared scorer env (e.g. a lab-wide 'cbulge_cpu') instead of
# creating its own. Resolved once at import; env_python() resolves ANY existing env by
# name (it need not appear in SCORER_ENVS).
DEFAULT_ENV = os.environ.get("CRISPRME_SCORER_ENV", "cbulge")

# ---------------------------------------------------------------------------
# CRISPR-Bulge source + weights (the model itself)
# ---------------------------------------------------------------------------
# The 5-model ensemble weights (~53 MB) are committed inside the CRISPR-Bulge repo
# (regular .h5 files; only the 524 MB datasets.zip is Git-LFS). So we provision the
# model by cloning the repo at a PINNED commit with LFS smudge disabled — the same
# pinned-source pattern crispritz uses. Overridable via env for air-gapped installs.
CBULGE_URL = os.environ.get("CRISPRME_CBULGE_URL", "https://github.com/OrensteinLab/CRISPR-Bulge.git")
CBULGE_PIN = os.environ.get("CRISPRME_CBULGE_PIN", "3eddcd5bfcaff00b2bdf29425116ec9756ace870")


def default_cbulge_repo() -> str:
    """Resolve where the CRISPR-Bulge source+weights live (no I/O).

    Precedence: $CBULGE_REPO > $CRISPRME_CBULGE_HOME > `<prefix>/opt/CRISPR-Bulge`
    derived from this file at `<prefix>/opt/crisprme/PostProcess` (sibling of the
    crisprme + crispritz opt trees).
    """
    for var in ("CBULGE_REPO", "CRISPRME_CBULGE_HOME"):
        v = os.environ.get(var)
        if v:
            return v
    here = os.path.dirname(os.path.abspath(__file__))       # <prefix>/opt/crisprme/PostProcess
    opt = os.path.dirname(os.path.dirname(here))            # <prefix>/opt
    return os.path.join(opt, "CRISPR-Bulge")


def source_present(repo: Optional[str] = None) -> bool:
    repo = repo or default_cbulge_repo()
    return os.path.isdir(os.path.join(repo, "OT_deep_score_src"))


def provision_source(dest: Optional[str] = None, pin: Optional[str] = None,
                     force: bool = False) -> Tuple[bool, str]:
    """Clone the CRISPR-Bulge repo at the pinned commit (LFS smudge OFF) if absent.

    Idempotent: a no-op success if the source is already present (unless force).
    Returns (ok, message). Requires `git` on PATH.
    """
    dest = dest or default_cbulge_repo()
    pin = pin or CBULGE_PIN
    if source_present(dest) and not force:
        return True, f"source present at {dest}"
    if not shutil.which("git"):
        return False, "git not found on PATH (cannot fetch CRISPR-Bulge source)"
    if force and os.path.isdir(dest):
        shutil.rmtree(dest, ignore_errors=True)
    env = dict(os.environ, GIT_LFS_SKIP_SMUDGE="1")
    sys.stderr.write(f"[scorer-env] provisioning CRISPR-Bulge @ {pin[:10]} -> {dest}\n")
    cp = _run(["git", "clone", CBULGE_URL, dest], env=env)
    if cp.returncode != 0:
        return False, f"git clone failed: {cp.stderr[-400:]}"
    cp = _run(["git", "-C", dest, "checkout", "--quiet", pin], env=env)
    if cp.returncode != 0:
        return False, f"git checkout {pin[:10]} failed: {cp.stderr[-400:]}"
    if not source_present(dest):
        return False, f"clone completed but OT_deep_score_src missing under {dest}"
    # drop .git history to slim the footprint (~179 MB -> ~53 MB source+weights);
    # re-provision with force re-clones from scratch
    shutil.rmtree(os.path.join(dest, ".git"), ignore_errors=True)
    return True, f"provisioned at {dest}"


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


def _variant(gpu: bool) -> str:
    """The TF variant to build: 'metal' on Apple-Silicon Macs (the only GPU stack there),
    'linux_aarch64' on Linux ARM (conda-forge has no tensorflow-cpu=2.13 there, so the same
    TF 2.13.1 comes from the PyPI aarch64 wheel), else 'cuda' when gpu is requested, else 'cpu'.
    On a Mac the ONE metal env serves both CPU and GPU runs -- the compute backend picks the
    device at runtime. aarch64 is always CPU (no CUDA/Metal there)."""
    sysname, machine = platform.system(), platform.machine()
    if sysname == "Darwin" and machine == "arm64":
        return "metal"
    if sysname == "Linux" and machine in ("aarch64", "arm64"):
        return "linux_aarch64"
    return "cuda" if gpu else "cpu"


def _conda_packages(spec: dict, variant: str) -> List[str]:
    return {
        "metal": spec.get("metal_packages"),
        "cuda": spec.get("gpu_packages"),
        "linux_aarch64": spec.get("linux_aarch64_packages"),
    }.get(variant, spec["cpu_packages"])


def _pip_packages(spec: dict, variant: str) -> Optional[List[str]]:
    """Packages a variant installs via pip AFTER the conda env (TensorFlow builds that
    conda-forge doesn't provide): the Metal TF on Apple Silicon, the PyPI aarch64 TF wheel
    on Linux ARM. None for the conda-only variants (cpu/cuda)."""
    return {
        "metal": spec.get("metal_pip"),
        "linux_aarch64": spec.get("linux_aarch64_pip"),
    }.get(variant)


def _channel_args(spec: dict) -> List[str]:
    """``-c`` args for the env's channels, honoring a mirror base so sites whose network
    blocks conda.anaconda.org can build the scorer env. If ``CRISPRME_CONDA_CHANNEL_BASE``
    or ``CONDA_CHANNEL_BASE`` is set (e.g. ``https://prefix.dev``), a bare channel name like
    ``conda-forge`` is rewritten to ``<base>/conda-forge`` -- the SAME convention the
    Dockerfile uses. Full-URL channels are passed through untouched; with no base set the
    plain name is used (which still honors any global condarc ``channel_alias``)."""
    base = (os.environ.get("CRISPRME_CONDA_CHANNEL_BASE")
            or os.environ.get("CONDA_CHANNEL_BASE") or "").rstrip("/")
    args: List[str] = []
    for ch in spec["channels"]:
        args += ["-c", f"{base}/{ch}" if (base and "://" not in ch) else ch]
    return args


def build_create_command(name: str, gpu: bool = False, variant: Optional[str] = None) -> Optional[List[str]]:
    """Return the argv to create the conda env, or None if no manager is available."""
    spec = SCORER_ENVS[name]
    mgr = detect_env_manager()
    if not mgr:
        return None
    exe, _ = mgr
    variant = variant or _variant(gpu)
    return (
        [exe, "create", "-y", "-n", name, "python=" + spec["python"]]
        + _channel_args(spec)
        + _conda_packages(spec, variant)
    )


def _create_conda_env(name, gpu, force, stream, variant) -> Tuple[bool, str]:
    """Create the conda env (no source provisioning). For the Metal variant, the Metal
    TensorFlow (tensorflow-macos + tensorflow-metal) is pip-installed after the conda env,
    since conda-forge has no Metal TF build."""
    if env_exists(name) and not force:
        return True, f"env '{name}' already exists"
    cmd = build_create_command(name, gpu=gpu, variant=variant)
    if cmd is None:
        return False, ("no conda env manager found (need micromamba/mamba/conda on "
                       "PATH, or set CRISPRME_CONDA_EXE)")
    if force and env_exists(name):
        _run([cmd[0], "env", "remove", "-y", "-n", name])
    env = dict(os.environ)
    if variant == "cuda" and not env.get("CONDA_OVERRIDE_CUDA"):
        # Force-set (NOT setdefault): the micromamba base-env activation exports
        # CONDA_OVERRIDE_CUDA="" (empty), which setdefault() would keep -> the conda CUDA
        # TensorFlow build then fails to solve on a GPU-less builder ("__cuda missing"),
        # e.g. GitHub CI. Treat empty as unset. 11.8 = TF 2.13's CUDA; the built image
        # still runs on CPU when no GPU is visible (compute-backend device guard).
        env["CONDA_OVERRIDE_CUDA"] = "11.8"
    sys.stderr.write(f"[scorer-env] creating '{name}' ({variant}): {' '.join(cmd)}\n")
    proc = subprocess.run(cmd, env=env) if stream else _run(cmd, env=env)
    if proc.returncode != 0:
        detail = "" if stream else (proc.stdout + proc.stderr)[-2000:]
        return False, f"env create failed (exit {proc.returncode}) {detail}"
    # Variants whose TensorFlow isn't on conda-forge (Metal on Apple Silicon; the aarch64
    # wheel on Linux ARM) pip-install it into the fresh env after the conda step.
    spec = SCORER_ENVS[name]
    pip_pkgs = _pip_packages(spec, variant)
    if pip_pkgs:
        exe = cmd[0]
        pip_cmd = [exe, "run", "-n", name, "pip", "install"] + pip_pkgs
        sys.stderr.write(f"[scorer-env] pip-installing TensorFlow ({variant}): {' '.join(pip_pkgs)}\n")
        pp = subprocess.run(pip_cmd) if stream else _run(pip_cmd)
        if pp.returncode != 0:
            return False, f"env created but TensorFlow pip install failed (exit {pp.returncode})"
    return True, f"env created ({variant})"


def create_env(name: str = DEFAULT_ENV, gpu: bool = False, force: bool = False,
               stream: bool = True) -> Tuple[bool, str]:
    """Create the scorer env AND provision its model source (idempotent).

    Two independent steps: (1) the conda env (CPU / CUDA / Metal variant, auto-picked by
    platform -- Metal on Apple Silicon), (2) the CRISPR-Bulge source+weights (pinned clone).
    Both must succeed.
    """
    if name not in SCORER_ENVS:
        return False, f"unknown scorer env '{name}'"
    variant = _variant(gpu)
    ok_env, msg_env = _create_conda_env(name, gpu, force, stream, variant)
    if not ok_env:
        return False, msg_env
    ok_src, msg_src = provision_source(force=force)
    return ok_src, f"{msg_env}; {msg_src}"


def update_env(name: str = DEFAULT_ENV, gpu: bool = False) -> Tuple[bool, str]:
    """Repair/update: create if missing, else re-assert the package set (idempotent)."""
    if not env_exists(name):
        return create_env(name, gpu=gpu)
    spec = SCORER_ENVS[name]
    mgr = detect_env_manager()
    if not mgr:
        return False, "no conda env manager found"
    exe, _ = mgr
    variant = _variant(gpu)
    cmd = [exe, "install", "-y", "-n", name] + _channel_args(spec) + _conda_packages(spec, variant)
    env = dict(os.environ)
    if variant == "cuda" and not env.get("CONDA_OVERRIDE_CUDA"):
        # Force-set (NOT setdefault): the micromamba base-env activation exports
        # CONDA_OVERRIDE_CUDA="" (empty), which setdefault() would keep -> the conda CUDA
        # TensorFlow build then fails to solve on a GPU-less builder ("__cuda missing"),
        # e.g. GitHub CI. Treat empty as unset. 11.8 = TF 2.13's CUDA; the built image
        # still runs on CPU when no GPU is visible (compute-backend device guard).
        env["CONDA_OVERRIDE_CUDA"] = "11.8"
    sys.stderr.write(f"[scorer-env] updating '{name}' ({variant}): {' '.join(cmd)}\n")
    proc = subprocess.run(cmd, env=env)
    if proc.returncode != 0:
        return False, f"update failed (exit {proc.returncode})"
    pip_pkgs = _pip_packages(spec, variant)
    if pip_pkgs:
        pp = subprocess.run([exe, "run", "-n", name, "pip", "install", "-U"] + pip_pkgs)
        if pp.returncode != 0:
            return False, f"conda update ok but TF pip update failed (exit {pp.returncode})"
    return True, f"updated ({variant})"


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
        # A CRISPRME_SCORER_ENV override makes DEFAULT_ENV a user/site-provided
        # CRISPR-Bulge env that isn't one of our named specs. Probe it against the
        # 'cbulge' spec (same import contract) instead of declaring it unknown. An
        # arbitrary unknown name (not the configured default) still errors as unknown.
        if name == DEFAULT_ENV:
            spec = SCORER_ENVS.get("cbulge")
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

    # model source + weights presence (pinned clone; see provision_source)
    repo = default_cbulge_repo()
    rec["cbulge_repo"] = repo
    if source_present(repo):
        rec["weights_present"] = True
    else:
        if rec["status"] != ERROR:
            rec["status"] = WARN
        rec["issues"].append((WARN, f"CRISPR-Bulge source/weights not found at {repo} "
                                    "(provision with: crisprme.py scorer-env create)"))
    return rec


def selftest(name: str = DEFAULT_ENV) -> tuple:
    """COMPUTE-based self-test: actually score known (sgRNA, off-target) pairs inside
    the env and assert the result is real and non-degenerate. Returns (ok, msg).

    This catches failure modes that an import-only probe (``health_check``) MISSES:
      * missing/corrupt model weights or a bad CBULGE_REPO -> all scores -1 (disabled);
      * a silently-miscomputing backend (e.g. the tensorflow-metal GRU bug that
        returned all 1.0) -> degenerate, all-identical scores.
    A perfect-match and a mismatched target MUST produce distinct, in-[0,1] scores.
    Build/CI/`scorer-env check` run this so a broken scorer env can never ship or
    silently degrade a search to CFD-only.
    """
    py = env_python(name)
    if not py:
        return False, ("env '%s' does not exist (create with: "
                       "crisprme.py scorer-env create)" % name)
    here = os.path.dirname(os.path.abspath(__file__))
    # perfect match vs a multi-mismatch off-target (23 nt incl. PAM). A correct scorer
    # gives a high score for the match and a clearly lower one for the mismatch.
    code = (
        "import json,sys\n"
        "sys.path.insert(0, " + repr(here) + ")\n"
        "import crispr_bulge_score as s\n"
        "sg =['GTAACGGCAGACTTCTCCACAGG','GTAACGGCAGACTTCTCCACAGG']\n"
        "off=['GTAACGGCAGACTTCTCCACAGG','GTCACGGCTGACTACTCCACAGG']\n"
        "print(json.dumps(s.CRISPR_BULGE_predict_list(sg, off)))\n"
    )
    try:
        cp = _run([py, "-c", code])
        r = json.loads(cp.stdout.strip().splitlines()[-1]) if cp.stdout.strip() else None
    except Exception as e:  # noqa: BLE001
        return False, "self-test failed to run: %s" % e
    if not r or any(x is None for x in r):
        return False, "self-test returned no scores: %r" % (r,)
    if all(x == -1 for x in r):
        return False, ("scorer is DISABLED/broken (all scores -1) -- missing env, "
                       "weights, or a bad CBULGE_REPO: %r" % (r,))
    if any((x < 0.0 or x > 1.0) for x in r):
        return False, "scores out of [0,1] (bad model/output): %r" % (r,)
    if len({round(x, 4) for x in r}) == 1:
        return False, ("degenerate all-identical scores -> suspect a miscomputing "
                       "backend (cf. the tensorflow-metal GRU bug): %r" % (r,))
    return True, "scorer computes real, distinct scores: %r" % (r,)


def build_all(gpu: bool = False, selftest_each: bool = True, stream: bool = True) -> tuple:
    """Provision (and COMPUTE-self-test) EVERY scorer env in the ``SCORER_ENVS``
    registry. The registry is the SINGLE SOURCE OF TRUTH for the modular scorer
    setup, so adding an env there automatically provisions it in the Docker image /
    any install -- the build never drifts from the registry. Returns (ok, report).

    ``selftest_each`` runs the compute self-test (not just an import probe) on each
    env so a broken/miscomputing scorer FAILS the build instead of silently shipping.
    """
    report = []
    names = list(SCORER_ENVS)
    report.append("building %d registered scorer env(s): %s"
                  % (len(names), ", ".join(names) or "(none)"))
    failed = []
    for n in names:
        ok, msg = create_env(n, gpu=gpu, stream=stream)
        report.append("%s create: %s" % (n, msg))
        if not ok:
            failed.append(n)
            continue
        if selftest_each:
            ok, msg = selftest(n)
            report.append("%s self-test: %s" % (n, msg))
            if not ok:
                failed.append(n)
    ok_all = not failed
    report.append("ALL scorer envs OK" if ok_all
                  else ("FAILED scorer env(s): " + ", ".join(failed)))
    return ok_all, report


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
