"""HuggingFace-backed fast download/upload of CRISPRme reference data and
precomputed indexes.

CRISPRme's reference bundle (genome, annotations, PAMs, sample-ID files, variant
VCFs) and its precomputed CRISPRitz indexes can be hosted on a HuggingFace
dataset repository and fetched over HF's CDN, which is typically far faster and
more reliable than the original UCSC/FTP sources. This module is the thin client
for that: it resolves the target repo, downloads requested components into the
canonical CRISPRme directory layout (decompressing where needed), and — on the
publish side — uploads a locally built index back to the repo.

Design notes
------------
* Everything network-facing goes through ``huggingface_hub``. The dependency is
  imported lazily so that importing this module (e.g. for the pure-logic unit
  tests) never requires it; only the functions that actually touch HF do.
* Public downloads need no token. Uploads (``publish_index``) need a write token,
  read from ``--token`` or the ``HF_TOKEN`` / ``HUGGING_FACE_HUB_TOKEN`` env vars.
* The remote layout mirrors the local one so mapping is trivial:

      genomes/<ref>/<chrom>.fa[.gz]      -> Genomes/<ref>/<chrom>.fa
      annotations/<file>                 -> Annotations/<file>
      PAMs/<file>                        -> PAMs/<file>
      samplesIDs/<file>                  -> samplesIDs/<file>
      VCFs/<dataset>/<file>.vcf.gz       -> VCFs/<dataset>/<file>.vcf.gz
      indexes/<index_name>.tar.gz        -> genome_library/<index_name>/

  Large FASTA is stored gzipped and decompressed after download; VCFs are kept
  bgzipped (CRISPRme consumes them that way); indexes travel as a single
  ``.tar.gz`` and are unpacked into ``genome_library/``.
"""

from typing import Dict, List, Optional
import os
import sys
import io
import json
import gzip
import shutil
import subprocess
import tarfile
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime, timezone

# Default HuggingFace dataset repo. Overridable per-invocation with --hf-repo or
# the CRISPRME_HF_REPO environment variable (e.g. switch to a pinellolab org
# repo once it exists).
DEFAULT_HF_REPO = "lucapinello/crisprme-data"

# PAM-less indexes with a flank WIDER than 3 nt live in their own dataset repo.
# A 3-nt flank is stored exactly like a dedicated 3-nt PAM index, so any CRISPRme
# version can read it and NNN_3 ships in the default repo. A wider flank cannot be
# read correctly by a release that predates --prune-pam, so keeping those out of the
# repo older defaults point at is what stops an old client from downloading an index
# it would misread. Routing is by flank width, not by name, so it cannot be forgotten
# at publish time.
GENERAL_INDEX_HF_REPO = "lucapinello/crisprme-data-general-index"

# A flank of this width or narrower is readable by every release, so it stays in the
# default repo.
_COMPAT_FLANK_WIDTH = 3

# Components understood by the download command and their remote sub-paths.
# This layout is the single source of truth for the HF dataset tree; the
# transparent setup/complete_test fast-path (utils.hf_fetch) targets the same
# paths (genomes/<ref>/, vcfs/<dataset>/, samplesIDs/, annotations/, indexes/).
_COMPONENT_PREFIXES = {
    "genome": "genomes",
    "annotations": "annotations",
    "pams": "pams",
    "samples": "samplesIDs",
    "vcf": "vcfs",
    "index": "indexes",
    # build-inputs/<dataset>.tar.gz = the PRE-ENRICHED genome kit (the intermediate build
    # path): staging it lets build-index-only skip the (slow) enrichment and only index,
    # so a new-PAM build is just a re-index. The raw-source path is `vcf` above.
    "build-inputs": "build-inputs",
}

# Local destination directory for each component (relative to the working dir).
_COMPONENT_LOCALDIR = {
    "genome": "Genomes",
    "annotations": "Annotations",
    "pams": "PAMs",
    "samples": "samplesIDs",
    "vcf": "VCFs",
    "index": "genome_library",
    "build-inputs": "Genomes",  # the enriched genome kit unpacks under Genomes/
}


# Publish-only markers a build may embed in the REF segment of an index dir name
# (<pam>_<N>_<ref>[+<vcf>]). They are NOT part of the search convention
# (ref segment == genome-folder basename), so download must strip them from the
# INSTALL dir name or the search cannot resolve the index (GAP 3).
_DICTLESS_MARKERS = ("-dictless",)


def canonical_index_name(index_name: str) -> str:
    """Strip a publish-only marker (e.g. ``-dictless``) from the REF segment of
    ``<pam>_<N>_<ref>[+<vcf>]`` so the installed dir matches the search convention
    (ref segment == genome-folder basename). NO-OP when already canonical.

    Only the ``base`` (``<pam>_<N>_<ref>``) segment before the FIRST ``+`` is
    mutated; the ``<vcf>`` segment after ``+`` is preserved verbatim (it equals
    ``vcf_name`` and names the shared ``genotypes_<vcf>`` companion, so it must
    never change). A reference-only name (no ``+``) or an already-canonical name
    passes through unchanged.
    """
    base, plus, vcf = index_name.partition("+")  # vcf untouched
    for mark in _DICTLESS_MARKERS:
        base = base.replace(mark, "")
    return base + plus + vcf


def _ensure_search_lists(workdir, vcf_name):
    """Make a downloaded variant index searchable from the CLI out of the box.

    The web form builds the per-search VCF / samplesID lists automatically; a CLI
    ``crisprme.py complete-search`` instead reads them from the ``--vcf`` /
    ``--samplesID`` list files. A fresh ``download --what index`` installs the
    index, registry and combined samplesID but NOT those small list files, so a
    CLI search previously required the user to hand-write them (and to know the
    internal dataset name). Create/extend them here: append this dataset to
    ``list_vcf.txt`` and its combined samplesID to ``list_samplesID.txt`` at the
    install root. Idempotent (each dataset added once), so installing several
    indexes accumulates them. Guarded, STDOUT-only, never fatal -- an optional
    convenience file, exactly like the combined samplesID.
    """
    try:
        entries = [("list_vcf.txt", vcf_name)]
        sid_name = f"{vcf_name}.samplesID.txt"
        if os.path.isfile(os.path.join(workdir, "samplesIDs", sid_name)):
            entries.append(("list_samplesID.txt", sid_name))
        for fname, value in entries:
            path = os.path.join(workdir, fname)
            lines = []
            if os.path.isfile(path):
                with open(path) as handle:
                    lines = [ln.strip() for ln in handle if ln.strip()]
            if value not in lines:
                lines.append(value)
                with open(path, "w") as handle:
                    handle.write("\n".join(lines) + "\n")
    except OSError as exc:
        sys.stdout.write(
            f"NOTE: could not write CLI search lists for '{vcf_name}' ({exc}); "
            f"create list_vcf.txt / list_samplesID.txt manually for CLI searches.\n"
        )


def synthesize_combined_samplesid(workdir, vcf_name, ref="hg38"):
    """Generate ``samplesIDs/<vcf_name>.samplesID.txt`` for a merged dataset when
    it is missing but the per-db component files are present (GAP 2).

    A merged variant index (e.g. ``hg38_1000G_HGDP``) searches ONE vcf folder and
    the search expects a single combined ``samplesIDs/<vcf_name>.samplesID.txt``,
    but the published artifact ships only the per-db lists
    (``hg38_1000G.samplesID.txt`` + ``hg38_HGDP.samplesID.txt``). This unions the
    per-db DATA rows (header-less, dedup by SAMPLE_ID, stable dataset order),
    mirroring ``pages/main_page._ensure_samplesid`` byte-for-byte so CLI == web.

    Returns the written path, or ``None`` (NO-OP) when the combined file already
    exists, the dataset is single (non-merged), or a component is missing.
    """
    sdir = os.path.join(workdir, "samplesIDs")
    target = os.path.join(sdir, f"{vcf_name}.samplesID.txt")
    if os.path.isfile(target):
        return None  # already present -> NO-OP (never touch a shipped file)
    # derive dataset tokens via the SHARED helper so the ensure-step
    # (_ensure_perdb_samplesids) and this synthesizer never disagree on which
    # per-db files matter. Empty -> single (non-merged) dataset -> nothing to do.
    dbs = _derive_perdb_datasets(vcf_name, ref)
    if not dbs:
        return None  # single (non-merged) dataset -> nothing to synthesize
    comps = [os.path.join(sdir, f"{ref}_{c}.samplesID.txt") for c in dbs]
    if not comps or not all(os.path.isfile(c) for c in comps):
        return None  # cannot synthesize (don't write a half-union)
    seen, rows = set(), []
    for c in comps:  # stable dataset order (e.g. 1000G then HGDP)
        with open(c) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                key = line.split("\t", 1)[0]
                if key not in seen:
                    seen.add(key)
                    rows.append(line.rstrip("\n"))
    os.makedirs(sdir, exist_ok=True)
    with open(target, "w") as out:
        out.write("\n".join(rows) + "\n")
    return target


def resolve_repo(cli_repo: Optional[str] = None) -> str:
    """Resolve the HuggingFace repo id: CLI value > env var > built-in default."""
    if cli_repo:
        return cli_repo
    return os.environ.get("CRISPRME_HF_REPO", DEFAULT_HF_REPO)


def index_flank_width(index_name: str) -> Optional[int]:
    """Stored flank width of a PAM-less index, or None if it is not PAM-less.

    An index directory/archive is ``<PAM>_<N>_<ref>[+<ref>_<vcf>]``. A PAM-less
    index declares an all-N PAM, and its width is that PAM's length.
    """
    if not index_name:
        return None
    pam = canonical_index_name(index_name).split("_", 1)[0]
    if pam and set(pam) == {"N"}:
        return len(pam)
    return None


def repo_overridden(cli_repo: Optional[str] = None) -> bool:
    """True when the caller named a repo explicitly (flag or env var)."""
    return bool(cli_repo) or bool(os.environ.get("CRISPRME_HF_REPO"))


def _index_archive_exists(repo: str, index_name: str) -> bool:
    """Best-effort check that ``indexes/<index_name>.tar.gz`` is in ``repo``.

    Returns False on any error (offline, repo absent, dependency missing) so a
    caller treats "cannot tell" the same as "not there" and falls back.
    """
    try:
        hf = _require_hf()
        want = f"{_COMPONENT_PREFIXES['index']}/{canonical_index_name(index_name)}.tar.gz"
        files = hf.HfApi().list_repo_files(repo, repo_type="dataset")
        return want in set(files)
    except Exception:  # noqa: BLE001
        return False


def resolve_index_repo(index_name: str, cli_repo: Optional[str] = None) -> str:
    """Which repo holds (or should hold) this index.

    An explicit ``--hf-repo`` / ``CRISPRME_HF_REPO`` always wins. Otherwise a
    PAM-less index whose flank is wider than a dedicated 3-nt PAM index routes to
    :data:`GENERAL_INDEX_HF_REPO`; everything else to :data:`DEFAULT_HF_REPO`.
    """
    if repo_overridden(cli_repo):
        return resolve_repo(cli_repo)
    width = index_flank_width(index_name)
    if width is not None and width > _COMPAT_FLANK_WIDTH:
        return GENERAL_INDEX_HF_REPO
    return DEFAULT_HF_REPO


def resolve_token(cli_token: Optional[str] = None) -> Optional[str]:
    """Resolve an HF token for uploads.

    Order: explicit ``--token`` > ``HF_TOKEN`` > ``HUGGING_FACE_HUB_TOKEN`` >
    a cached ``huggingface-cli login`` token. The cached-login fallback means
    that on a machine already logged into HuggingFace, ``publish-index`` works
    with no extra flags or env vars (matching how other HF tools behave).
    """
    tok = (
        cli_token
        or os.environ.get("HF_TOKEN")
        or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    )
    if tok:
        return tok
    try:  # cached `huggingface-cli login` token, if any (no network)
        from huggingface_hub import get_token

        return get_token()
    except Exception:  # huggingface_hub absent or no cached token
        return None


def _require_hf():
    """Import huggingface_hub lazily, with an actionable error if it is missing."""
    try:
        import huggingface_hub  # noqa: F401

        return huggingface_hub
    except ImportError as e:  # pragma: no cover - exercised only without the dep
        raise ImportError(
            "The 'huggingface_hub' package is required for HuggingFace "
            "downloads/uploads but is not installed. Install it with "
            "'conda install -c conda-forge huggingface_hub' or "
            "'pip install huggingface_hub'."
        ) from e


def list_available_downloads(
    component: str, repo: Optional[str] = None
) -> List[Dict]:
    """List items available to download from the HuggingFace dataset repo.

    Returns ``[{"name": str, "size": int}, ...]`` for the given component:
    ``genome``/``vcf`` -> the assembly/dataset folders under ``genomes/`` or
    ``vcfs/`` (size = sum of that folder's files); ``index`` -> the
    ``indexes/<name>.tar.gz`` archives. Returns ``[]`` on any error (offline,
    missing repo, dependency absent) so callers can fall back gracefully.

    For ``index`` with no repo named, both the default repo and the wide-flank
    PAM-less repo are listed and merged, each row taken from the copy ``download``
    would fetch: every default-repo row (a legacy name is fetched only there, and a
    wide-flank name falls back there), and a wide-flank-repo row only for a name
    routed to that repo, where it replaces the default repo's copy. A legacy-named
    archive that exists only in the wide-flank repo is therefore not offered, since
    ``download`` would never look there for it. The wide-flank repo being
    unreadable just drops its rows; the default repo being unreadable returns
    ``[]``, the offline signal, rather than a catalog of only wide-flank indexes.
    """
    if component == "index" and not repo_overridden(repo):
        try:
            default_rows = _list_available_downloads_raw("index", DEFAULT_HF_REPO)
        except Exception:
            return []
        merged: Dict[str, Dict] = {row["name"]: row for row in default_rows}
        for row in _list_available_downloads_one("index", GENERAL_INDEX_HF_REPO):
            if resolve_index_repo(row["name"]) == GENERAL_INDEX_HF_REPO:
                merged[row["name"]] = row
        return sorted(merged.values(), key=lambda row: row["name"])
    return _list_available_downloads_one(component, repo)


def _list_available_downloads_one(
    component: str, repo: Optional[str] = None
) -> List[Dict]:
    """One repo's listing (``[]`` if that repo cannot be read)."""
    try:
        return _list_available_downloads_raw(component, repo)
    except Exception:
        return []


def _list_available_downloads_raw(
    component: str, repo: Optional[str] = None
) -> List[Dict]:
    """One repo's listing; raises on any error, so a caller can tell a repo that
    cannot be read from one that holds nothing."""
    if component not in ("genome", "vcf", "index"):
        return []
    hf = _require_hf()
    repo = resolve_repo(repo)
    api = hf.HfApi()
    prefix = _COMPONENT_PREFIXES[component]
    if component == "index":
        items = []
        gt_present = set()  # vcf names that have a genotypes_<vcf>.tar.gz companion
        for entry in api.list_repo_tree(
            repo, repo_type="dataset", path_in_repo=prefix
        ):
            size = getattr(entry, "size", None)
            path = getattr(entry, "path", "")
            if size is None or not path.endswith(".tar.gz"):
                continue
            name = os.path.basename(path)[: -len(".tar.gz")]
            # The separate Tier-1 genotype companions live under the SAME
            # indexes/ prefix but are NOT selectable indexes themselves — they
            # ride along with their index on download. Record which vcfs have
            # one (so the index rows can advertise it) and skip them as rows.
            if name.startswith("genotypes_"):
                gt_present.add(name[len("genotypes_"):])
                continue
            items.append({"name": name, "size": size})
        # ADDITIVE metadata: mark whether each variant index has a genotype
        # companion available (has_genotypes). Consumers that only read
        # name/size are unaffected (extra keys are ignored).
        for it in items:
            vcf = it["name"].partition("+")[2]
            it["has_genotypes"] = bool(vcf) and vcf in gt_present
        return sorted(items, key=lambda d: d["name"])
    # genome / vcf: folders under prefix/, sum their files' sizes
    sizes: Dict[str, int] = {}
    for entry in api.list_repo_tree(
        repo, repo_type="dataset", path_in_repo=prefix, recursive=True
    ):
        size = getattr(entry, "size", None)
        path = getattr(entry, "path", "")
        if size is None:  # a folder entry, not a file
            continue
        rel = path[len(prefix) + 1 :]
        top = rel.split("/")[0]
        if top:
            sizes[top] = sizes.get(top, 0) + size
    return [{"name": k, "size": v} for k, v in sorted(sizes.items())]


def _catalog_component_size(
    component: str, name: str, repo: Optional[str] = None
) -> int:
    """Best-effort download size (bytes) of a catalog item, 0 if unknown.

    Reuses :func:`list_available_downloads` (which already returns per-item
    sizes) so the progress bar has a denominator up-front. Any failure (offline,
    missing item) -> 0, so the bar simply shows transferred bytes with no total."""
    try:
        for item in list_available_downloads(component, repo):
            if item.get("name") == name:
                return int(item.get("size", 0) or 0)
    except Exception:  # noqa: BLE001 - sizing is best-effort
        return 0
    return 0


def decompress_gz(path: str, remove_source: bool = True) -> str:
    """Gunzip ``path`` (``*.gz`` -> ``*``) and return the decompressed path.

    Used for FASTA that is stored gzipped on HF to save space/bandwidth but must
    be plain text for CRISPRme/CRISPRitz. VCFs are intentionally NOT run through
    this (they must stay bgzipped).
    """
    if not path.endswith(".gz"):
        return path
    out = path[:-3]
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    if remove_source:
        os.remove(path)
    return out


# --- download progress -------------------------------------------------------
# The web Settings download used to show only a spinner + a static "running..."
# while a multi-GB fetch ran for minutes; the underlying huggingface_hub tqdm
# progress goes to stderr (log_error.txt) and was never parsed. This emits a
# small JSON progress file the web poller renders as a real bar. It is
# ENV-GATED: without CRISPRME_PROGRESS_FILE set, every function here is a no-op
# and the happy path is byte-for-byte unchanged. All failures are swallowed —
# progress reporting must never break a download.


class _ProgressReporter:
    """Write a throttled JSON progress file for the web Settings poller.

    The file path comes from the ``CRISPRME_PROGRESS_FILE`` environment variable
    (set by the web launcher, see settings_page._run_settings_job). When it is
    unset this is a pure no-op — there is zero behaviour change on the CLI / the
    current happy path. The JSON schema the poller reads::

        {"phase": str, "bytes_done": int, "bytes_total": int, "pct": float,
         "eta_seconds": int|null, "detail": str, "updated": <epoch float>}

    Writes are throttled (``min_interval`` seconds, default 1.0) and atomic
    (write a temp sibling + ``os.replace``) so a concurrent reader never sees a
    half-written file.
    """

    def __init__(self, path: Optional[str] = None, min_interval: float = 1.0):
        if path is None:
            path = os.environ.get("CRISPRME_PROGRESS_FILE") or None
        self.path = path
        self.min_interval = min_interval
        self.phase = ""
        self.detail = ""
        self.bytes_total = 0
        self.bytes_done = 0
        self._base_done = 0  # bytes finished in PRIOR stages (multi-part index)
        self._start = time.time()
        self._last_write = 0.0

    @property
    def enabled(self) -> bool:
        return bool(self.path)

    def set_total(self, total: int) -> None:
        """Set the known total byte size for the current stage (0 == unknown)."""
        try:
            self.bytes_total = int(total or 0)
        except (TypeError, ValueError):
            self.bytes_total = 0

    def set_phase(self, phase: str, detail: str = "", carry_done: bool = False) -> None:
        """Start a new phase. ``carry_done`` keeps the running byte total across
        stages (a multi-part index: index 1/3 -> genotypes 2/3 -> reference 3/3);
        otherwise the per-stage byte counter resets."""
        self.phase = phase
        self.detail = detail
        if carry_done:
            self._base_done = self.bytes_done
        else:
            self._base_done = 0
            self.bytes_done = 0
        self.write(force=True)

    def add_bytes(self, n: int) -> None:
        """Accumulate ``n`` newly-transferred bytes and write (throttled)."""
        try:
            self.bytes_done += int(n or 0)
        except (TypeError, ValueError):
            return
        self.write()

    def set_done(self, done: int) -> None:
        """Set the absolute bytes-done for the current stage (decompress poll)."""
        try:
            self.bytes_done = self._base_done + int(done or 0)
        except (TypeError, ValueError):
            return
        self.write()

    def write(self, force: bool = False) -> None:
        if not self.path:
            return
        now = time.time()
        if not force and (now - self._last_write) < self.min_interval:
            return
        self._last_write = now
        total = self.bytes_total
        done = min(self.bytes_done, total) if total else self.bytes_done
        pct = (100.0 * done / total) if total else 0.0
        eta = None
        elapsed = now - self._start
        if total and done > 0 and elapsed > 0:
            rate = done / elapsed
            if rate > 0:
                eta = int(max(0, (total - done) / rate))
        payload = {
            "phase": self.phase,
            "bytes_done": int(done),
            "bytes_total": int(total),
            "pct": round(pct, 1),
            "eta_seconds": eta,
            "detail": self.detail,
            "updated": now,
        }
        try:
            tmp = self.path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(payload, fh)
            os.replace(tmp, self.path)
        except OSError:
            pass  # progress reporting is best-effort; never raise


def _make_progress_tqdm(reporter: "_ProgressReporter"):
    """Build a ``tqdm_class`` subclass that feeds a :class:`_ProgressReporter`.

    ``huggingface_hub.snapshot_download`` accepts a ``tqdm_class`` and drives it
    with per-file byte updates; we subclass tqdm so each ``update(n)`` also
    accumulates into the reporter. Across the many-small-file index tree tqdm
    creates several bars, so we sum ``n`` across ALL instances (hence the
    reporter, not per-bar totals). Returns ``None`` when tqdm is unavailable so
    the caller falls back to a plain (progress-free) download.
    """
    try:
        from tqdm.auto import tqdm as _tqdm
    except Exception:  # noqa: BLE001 - tqdm absent -> no custom bar, plain download
        return None

    class _ProgressTqdm(_tqdm):  # pragma: no cover - exercised via a fake in tests
        def update(self, n=1):
            try:
                reporter.add_bytes(n or 0)
            except Exception:  # noqa: BLE001 - never let reporting break the bar
                pass
            return super().update(n)

    return _ProgressTqdm


def _dir_bytes(path: str) -> int:
    """Total bytes of regular files under ``path`` (no symlink following)."""
    total = 0
    for root, _dirs, files in os.walk(path):
        for fn in files:
            fp = os.path.join(root, fn)
            try:
                if not os.path.islink(fp):
                    total += os.path.getsize(fp)
            except OSError:
                pass
    return total


@contextmanager
def _decompress_progress(
    reporter: Optional["_ProgressReporter"], out_dir: str, phase: str
):
    """Context manager: while the body runs a local decompress/extract into
    ``out_dir``, a background thread polls its growing on-disk size and pushes it
    to ``reporter`` as phase=``phase``. A strict no-op when progress is disabled —
    no thread is started, so the CLI path is unchanged. The byte total already
    reflects the known download size (set by the caller), so the growing extracted
    size is a reasonable denominator for the decompress stage too."""
    if reporter is None or not reporter.enabled:
        yield
        return
    import threading

    reporter.set_phase(phase, carry_done=False)
    stop_evt = threading.Event()

    def _poll():
        while not stop_evt.is_set():
            reporter.set_done(_dir_bytes(out_dir))
            stop_evt.wait(1.0)

    th = threading.Thread(target=_poll, daemon=True)
    th.start()
    try:
        yield
    finally:
        stop_evt.set()
        th.join(timeout=2.0)
        reporter.set_done(_dir_bytes(out_dir))


# --- atomic directory install -----------------------------------------------
# A download + decompress that writes straight into the final data-tree location
# (Genomes/<ref>/, VCFs/<ds>/, ...) leaves a half-populated dir on any mid-way
# failure (interrupt, disk-full, flaky mount) — the search then sees a
# discoverable-but-incomplete genome/dataset and fails cryptically later. Build
# into a hidden staging dir on the SAME filesystem, run the verifier, then
# os.replace the whole dir into place as the LAST step. On ANY exception the
# staging dir is discarded so the final location is never touched by a partial.


def _atomic_install_dir(dest: str, build, verify=None) -> str:
    """Populate ``dest`` atomically via a same-filesystem staging dir.

    ``build(stage_dir)`` writes the component into a fresh hidden staging dir
    (created as a sibling of ``dest`` so ``os.replace`` is an atomic rename, not a
    cross-device copy). ``verify(stage_dir)`` — if given — may raise to reject a
    bad staging dir BEFORE anything is swapped in. Only after a clean build+verify
    is the staged dir swapped into ``dest`` (any pre-existing ``dest`` is moved
    aside first and deleted after the swap; rolled-forward if a prior run died
    mid-swap). On ANY exception the staging dir is removed and ``dest`` is left
    exactly as it was — no partial is ever installed. Returns ``dest``.

    A ``dest`` that is a symlink (e.g. Genomes/hg38 linked to a shared copy) is
    replaced by the new real directory: the old link is removed with ``os.unlink``,
    never rmtree'd through, and the shared target is never written to or deleted."""
    parent = os.path.dirname(dest.rstrip("/")) or "."
    os.makedirs(parent, exist_ok=True)
    base = os.path.basename(dest.rstrip("/"))
    stage = os.path.join(parent, f".{base}.staging")
    backup = os.path.join(parent, f".{base}.replaced")
    # roll forward an interrupted prior swap: a backup with no live dest means a
    # previous run died between the two renames below -> restore it first.
    if os.path.isdir(backup) and not os.path.lexists(dest):
        os.rename(backup, dest)
    _remove_moved_aside(backup)
    shutil.rmtree(stage, ignore_errors=True)
    os.makedirs(stage)
    try:
        build(stage)
        if verify is not None:
            verify(stage)
        # clean build -> swap in as the last step (same-FS renames == atomic). Both
        # renames are inside this try, so a failed swap still discards the stage.
        link_target = os.readlink(dest) if os.path.islink(dest) else None
        if os.path.lexists(dest):
            os.rename(dest, backup)
        try:
            os.rename(stage, dest)
        except BaseException:
            if os.path.lexists(backup) and not os.path.lexists(dest):
                os.rename(backup, dest)  # restore on a failed final swap
            raise
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)  # discard the partial; dest intact
        raise
    _remove_moved_aside(backup)
    if link_target is not None:
        sys.stdout.write(
            f"NOTE: {dest} was a symbolic link to {link_target}; replaced the link "
            f"with the downloaded copy and left {link_target} untouched.\n"
        )
    return dest


def _remove_moved_aside(path: str) -> None:
    """Delete a moved-aside ``dest``: rmtree a real directory, but only unlink a
    symlink (or file), so a linked shared copy is never followed. Best-effort."""
    if os.path.islink(path) or os.path.isfile(path):
        try:
            os.unlink(path)
        except OSError:
            pass
    else:
        shutil.rmtree(path, ignore_errors=True)


# --- download integrity -----------------------------------------------------
# A truncated/interrupted download (partial file, flaky mount, disk full, a proxy
# that returns a short body) used to install silently and only blow up later at
# search time with a cryptic htslib error ("EOF marker is absent … could not open
# index"). We verify every downloaded file and re-fetch the bad ones before they
# are moved into place.

# A valid bgzipped file (*.bed.gz, *.vcf.gz, *.tbi, …) ends with this 28-byte
# empty-BGZF-block "EOF marker"; its absence is htslib's definition of truncated.
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
# A regular gzip stream (FASTA .fa.gz, index .tar.gz) ends with an 8-byte trailer
# (CRC32 + ISIZE); it has magic 1f 8b at the start. These are cheap structural
# sanity checks; a full-stream scan (interior corruption) is opt-in (see below).
_GZIP_MAGIC = b"\x1f\x8b"


class DownloadIntegrityError(RuntimeError):
    """A downloaded file failed integrity checks and could not be repaired."""


def _has_bgzf_eof(path: str) -> bool:
    try:
        if os.path.getsize(path) < len(_BGZF_EOF):
            return False
        with open(path, "rb") as fh:
            fh.seek(-len(_BGZF_EOF), os.SEEK_END)
            return fh.read(len(_BGZF_EOF)) == _BGZF_EOF
    except OSError:
        return False


def _gzip_stream_ok(path: str) -> bool:
    """Full decompress-scan — catches interior corruption, but reads the whole
    (possibly multi-GB) file, so it is opt-in via CRISPRME_VERIFY_DEEP=1."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:  # noqa: BLE001 — any read/decompress error == corrupt
        return False


def _tar_stream_ok(path: str) -> bool:
    try:
        with tarfile.open(path) as tf:
            tf.getmembers()
        return True
    except Exception:  # noqa: BLE001 — any error opening the archive == corrupt
        return False


def verify_download_file(path: str, expected_size: Optional[int] = None,
                         deep: Optional[bool] = None) -> Optional[str]:
    """Return None if ``path`` looks intact, else a short human reason.

    Cheap by default: the file exists, is non-empty, matches HF's advertised
    size (the reliable truncation catch), and passes a per-type structural check
    (BGZF EOF marker for *.bed.gz/*.vcf.gz/*.tbi; gzip magic for other *.gz). Set
    ``CRISPRME_VERIFY_DEEP=1`` (or deep=True) to also scan the full gzip/tar
    stream for interior corruption.
    """
    if deep is None:
        deep = os.environ.get("CRISPRME_VERIFY_DEEP") == "1"
    if not os.path.isfile(path):
        return "missing"
    size = os.path.getsize(path)
    if size == 0:
        return "empty (0 bytes)"
    if expected_size is not None and size != expected_size:
        return f"size {size:,} B != expected {expected_size:,} B (truncated/partial)"
    low = path.lower()
    if low.endswith((".bed.gz", ".vcf.gz", ".tbi", ".csi", ".bgz")):
        if not _has_bgzf_eof(path):
            return "BGZF EOF marker absent (truncated)"
        if deep and not _gzip_stream_ok(path):
            return "gzip stream corrupt"
    elif low.endswith((".tar.gz", ".tgz")):
        with open(path, "rb") as fh:
            if fh.read(2) != _GZIP_MAGIC:
                return "not a gzip stream (corrupt/partial)"
        if deep and not _tar_stream_ok(path):
            return "tar archive corrupt"
    elif low.endswith(".gz"):
        with open(path, "rb") as fh:
            if fh.read(2) != _GZIP_MAGIC:
                return "not a gzip stream (corrupt/partial)"
        if deep and not _gzip_stream_ok(path):
            return "gzip stream corrupt"
    return None


def _remote_sizes(hf, repo: str, rel_paths: List[str],
                  token: Optional[str]) -> Dict[str, int]:
    """Best-effort {remote_path: size} from HF metadata (the ground truth for a
    truncation check). Any failure (old hub, offline) -> {} = structural-only."""
    sizes: Dict[str, int] = {}
    try:
        for i in range(0, len(rel_paths), 256):  # chunk large trees
            chunk = rel_paths[i:i + 256]
            for info in hf.get_paths_info(
                repo, chunk, repo_type="dataset", token=token
            ):
                size = getattr(info, "size", None)
                if size is not None:
                    sizes[info.path] = size
    except Exception:  # noqa: BLE001 — verification degrades gracefully
        return {}
    return sizes


def _staged_files(local_dir: str) -> List[str]:
    """Relative paths of downloaded files under ``local_dir`` (skip hub cache)."""
    out: List[str] = []
    for root, dirs, files in os.walk(local_dir):
        dirs[:] = [d for d in dirs if d != ".cache"]  # hf_hub metadata dir
        for fn in files:
            rel = os.path.relpath(os.path.join(root, fn), local_dir)
            if not rel.startswith(".cache" + os.sep):
                out.append(rel)
    return out


# HTTP statuses worth retrying: 429 (rate limit) + transient 5xx from the CDN.
_HF_RETRY_STATUS = (429, 500, 502, 503, 504)


def _hf_retryable(exc: BaseException) -> bool:
    """True when ``exc`` is a transient HF HTTP error worth retrying — a 429
    rate-limit or a 5xx. Prefers the structured ``response.status_code`` and
    falls back to matching the status in the message (some hub errors don't
    expose ``.response``)."""
    status = getattr(getattr(exc, "response", None), "status_code", None)
    if status in _HF_RETRY_STATUS:
        return True
    msg = str(exc)
    return (
        "Too Many Requests" in msg
        or "rate limit" in msg.lower()
        or "Service Unavailable" in msg
        or any((" %d " % c) in (" %s " % msg) for c in _HF_RETRY_STATUS)
    )


def _hf_retry_after(exc: BaseException) -> Optional[float]:
    """Return the server-advertised ``Retry-After`` (seconds) if present + sane,
    so we wait exactly as long as the rate-limiter asks instead of guessing."""
    headers = getattr(getattr(exc, "response", None), "headers", None) or {}
    raw = headers.get("Retry-After") or headers.get("retry-after")
    if not raw:
        return None
    try:
        secs = float(raw)
    except (TypeError, ValueError):
        return None
    # cap a hostile/absurd value so a bad header can't wedge the download
    return secs if 0 < secs <= 300 else None


def _hf_call_with_retry(fn, *args, _what: str = "download",
                        _max_attempts: int = 6, **kwargs):
    """Call ``fn(*args, **kwargs)``, retrying transient HF HTTP errors (429
    rate-limit + 5xx) with exponential backoff. Anonymous downloads of the
    many-small-file indices hit 429 easily; without this the whole download
    aborts on the first rate-limit. ``snapshot_download`` resumes (already-present
    files are skipped) so a retry is cheap. Honors CRISPRME_HF_MAX_RETRIES; set
    it to 1 to disable retrying."""
    try:
        _max_attempts = max(
            1, int(os.environ.get("CRISPRME_HF_MAX_RETRIES", _max_attempts))
        )
    except (TypeError, ValueError):
        pass
    delay = 4.0
    for attempt in range(1, _max_attempts + 1):
        try:
            return fn(*args, **kwargs)
        except Exception as e:  # noqa: BLE001 - re-raised below if not retryable
            if attempt >= _max_attempts or not _hf_retryable(e):
                raise
            wait = _hf_retry_after(e)
            if wait is None:
                wait = delay
            status = getattr(getattr(e, "response", None), "status_code", "?")
            sys.stderr.write(
                f"HF {_what}: transient error (HTTP {status}); retry "
                f"{attempt}/{_max_attempts - 1} in {wait:.0f}s "
                f"(set CRISPRME_HF_MAX_RETRIES=1 to disable, or authenticate "
                f"with an HF token to avoid rate limits)...\n"
            )
            sys.stderr.flush()
            time.sleep(wait)
            delay = min(delay * 2, 60.0)
    raise RuntimeError("unreachable")  # pragma: no cover


def _hf_snapshot(repo: str, allow_patterns: List[str], local_dir: str,
                 token: Optional[str] = None, retries: int = 2,
                 reporter: Optional["_ProgressReporter"] = None) -> str:
    """Download the files matching ``allow_patterns`` from ``repo`` into
    ``local_dir`` (flat, following the repo tree), VERIFY each one against HF's
    advertised size + a structural check, and re-fetch any truncated/corrupt file
    (up to ``retries`` extra passes) before returning. Raises
    ``DownloadIntegrityError`` if a file can't be made whole. Returns
    ``local_dir``. Set CRISPRME_SKIP_VERIFY=1 to bypass (not recommended).

    When ``reporter`` is an enabled :class:`_ProgressReporter`, a custom
    ``tqdm_class`` is passed to ``snapshot_download`` so transferred bytes feed
    the web progress file. The extra kwarg is omitted entirely when progress is
    off (or tqdm is unavailable), keeping the happy path byte-for-byte unchanged.
    """
    hf = _require_hf()
    os.makedirs(local_dir, exist_ok=True)
    extra = {}
    if reporter is not None and reporter.enabled:
        tqdm_cls = _make_progress_tqdm(reporter)
        if tqdm_cls is not None:
            extra["tqdm_class"] = tqdm_cls
    _hf_call_with_retry(
        hf.snapshot_download,
        repo_id=repo,
        repo_type="dataset",
        allow_patterns=allow_patterns,
        local_dir=local_dir,
        token=token,
        _what="snapshot",
        **extra,
    )
    if os.environ.get("CRISPRME_SKIP_VERIFY") == "1":
        return local_dir

    rels = _staged_files(local_dir)
    # _remote_sizes returns {} when HF's size API fails (old hub / offline). That
    # is NOT a reason to skip verification: when sizes.get(rel) is None below,
    # verify_download_file skips ONLY the size comparison and STILL runs the
    # BGZF-EOF / gzip-magic / tar structural checks. So a truncated bgzf is caught
    # even with no advertised sizes.
    sizes = _hf_call_with_retry(
        _remote_sizes, hf, repo, rels, token, _what="size-check"
    )
    for attempt in range(retries + 1):
        bad = {}
        for rel in _staged_files(local_dir):
            reason = verify_download_file(
                os.path.join(local_dir, rel), sizes.get(rel)
            )
            if reason:
                bad[rel] = reason
        if not bad:
            return local_dir
        if attempt == retries:
            detail = "\n".join(f"  - {r}: {why}" for r, why in sorted(bad.items()))
            raise DownloadIntegrityError(
                f"{len(bad)} file(s) from '{repo}' are still corrupt after "
                f"{retries + 1} download attempt(s):\n{detail}\n"
                "This is usually a flaky network/mount or low disk space. Free up "
                "space, check the connection, and re-run the same download command "
                "(it resumes and re-fetches only the bad files)."
            )
        sys.stderr.write(
            f"Integrity check: re-downloading {len(bad)} truncated/corrupt "
            f"file(s) (attempt {attempt + 2}/{retries + 1})...\n"
        )
        sys.stderr.flush()
        for rel in bad:  # drop the bad copy so the re-fetch cannot be skipped
            try:
                os.remove(os.path.join(local_dir, rel))
            except OSError:
                pass
        _hf_call_with_retry(
            hf.snapshot_download,
            repo_id=repo,
            repo_type="dataset",
            allow_patterns=list(bad),  # exact remote paths of the bad files
            local_dir=local_dir,
            token=token,
            force_download=True,
            _what="snapshot-refetch",
        )
    return local_dir


def _fetch_genotypes_companion(
    repo: str,
    workdir: str,
    vcf_name: str,
    remote_prefix: str,
    token: Optional[str] = None,
    reporter: Optional["_ProgressReporter"] = None,
) -> bool:
    """Fetch + extract the SEPARATE Tier-1 ``genotypes_<vcf>.tar.gz`` companion
    into ``<workdir>/Dictionaries/`` so ``genotypes_<vcf>/`` lands as a sibling of
    ``dictionaries_``/``registry_`` (the exact path the search resolver reads).

    Returns ``True`` when the store was installed, ``False`` when the repo has no
    such companion (a registry-only / detection-only publish). Absence is NOT an
    error: the search falls back gracefully (degraded Samples). Uses the same
    atomic-swap staging as the main index install; all diagnostics go to STDOUT
    (stderr is fatal downstream)."""
    remote_name = f"genotypes_{vcf_name}.tar.gz"
    patterns = [f"{remote_prefix}/{remote_name}"]
    staging = os.path.join(workdir, ".hf_stage_genotypes")
    shutil.rmtree(staging, ignore_errors=True)
    _hf_snapshot(repo, patterns, staging, token, reporter=reporter)
    gt_tarball = os.path.join(staging, remote_prefix, remote_name)
    if not os.path.isfile(gt_tarball):
        shutil.rmtree(staging, ignore_errors=True)
        return False
    dst_dicts_root = os.path.join(workdir, "Dictionaries")
    os.makedirs(dst_dicts_root, exist_ok=True)
    extract_tmp = os.path.join(dst_dicts_root, f".extract_genotypes_{vcf_name}")
    shutil.rmtree(extract_tmp, ignore_errors=True)
    os.makedirs(extract_tmp)
    try:
        with _decompress_progress(
            reporter, extract_tmp, phase="Extracting genotypes"
        ), tarfile.open(gt_tarball) as tf:
            tf.extractall(extract_tmp)
        staged = os.path.join(extract_tmp, "Dictionaries")
        # tolerate an archive rooted at genotypes_<vcf>/ directly (no Dictionaries/
        # prefix) as well as the canonical Dictionaries/genotypes_<vcf>/ layout.
        src_root = staged if os.path.isdir(staged) else extract_tmp
        installed_any = False
        for sub in os.listdir(src_root):
            src = os.path.join(src_root, sub)
            if not os.path.isdir(src):
                continue
            dst = os.path.join(dst_dicts_root, sub)
            backup = os.path.join(dst_dicts_root, f".{sub}.replaced")
            if os.path.isdir(backup) and not os.path.exists(dst):
                os.rename(backup, dst)
            shutil.rmtree(backup, ignore_errors=True)
            if os.path.exists(dst):
                os.rename(dst, backup)
            os.rename(src, dst)
            shutil.rmtree(backup, ignore_errors=True)
            installed_any = True
        return installed_any
    finally:
        shutil.rmtree(extract_tmp, ignore_errors=True)
        shutil.rmtree(staging, ignore_errors=True)


def _derive_perdb_datasets(vcf_name: str, ref: str = "hg38") -> List[str]:
    """Derive the per-db dataset tokens of a MERGED variant index, EXACTLY the way
    :func:`synthesize_combined_samplesid` does (lines 116-120), so the ensure-step
    and the synthesizer never diverge on which per-db files matter.

    ``vcf_name`` is the segment after the first ``+`` of an index name (e.g.
    ``hg38_1000G_HGDP``). Strip a leading ``<ref>_`` then split on ``_``. Returns
    ``[]`` for a SINGLE (non-merged) dataset (no ``_`` after the ref-strip), which
    is the exact guard the synthesizer no-ops on — a single dataset needs no
    combined file, and its per-db list IS the samplesID the search uses.
    """
    dataset = vcf_name[len(ref) + 1:] if vcf_name.startswith(ref + "_") else vcf_name
    if "_" not in dataset:
        return []  # single (non-merged) dataset -> nothing to ensure/synthesize
    return dataset.split("_")


def _ensure_perdb_samplesids(
    repo: str,
    workdir: str,
    vcf_name: str,
    ref: str = "hg38",
    token: Optional[str] = None,
) -> None:
    """Ensure the per-db samplesID files for a MERGED variant index are present in
    ``<workdir>/samplesIDs/`` before :func:`synthesize_combined_samplesid` runs,
    fetching any that are missing from HF's ``samplesIDs/`` component (GAP 2b).

    Rationale: ``synthesize_combined_samplesid`` unions the per-db lists
    (``<ref>_<db>.samplesID.txt``) into the combined
    ``<ref>_<db1>_<db2>...samplesID.txt`` the search's ``--samplesID`` expects, but
    ONLY when those per-db files are already on disk. They arrive via
    ``download --what all`` / ``--what samples`` (the 'samples' component), NOT via
    ``download --what index``. A standalone ``--what index`` therefore had no per-db
    files and the synthesis silently no-op'd. This fetches the genuinely-missing
    ones — reusing the exact HF snapshot mechanism (_hf_snapshot + exact-path
    allow_patterns) the 'samples' component / genotypes companion already use — so
    the standalone path works too.

    Strict NO-OP when every per-db file is already present (the ``--what all`` path:
    zero extra network calls, no staging dir, no STDOUT noise) and for a single
    (non-merged) dataset. Fully guarded: any network/HF error is reported to STDOUT
    and swallowed (stderr is fatal to CRISPRme post-analysis) — degrading to today's
    silent no-op but with a diagnostic; the caller never crashes over an optional
    convenience file.
    """
    try:
        dbs = _derive_perdb_datasets(vcf_name, ref)
        if not dbs:
            return  # single-dataset index -> nothing to ensure (matches synthesizer)
        sdir = os.path.join(workdir, "samplesIDs")
        # Cheap short-circuit: if the combined target already exists, the synthesizer
        # will no-op at its first guard anyway -> nothing to ensure.
        if os.path.isfile(os.path.join(sdir, f"{vcf_name}.samplesID.txt")):
            return
        remote_prefix = _COMPONENT_PREFIXES["samples"]
        missing = [
            db for db in dbs
            if not os.path.isfile(os.path.join(sdir, f"{ref}_{db}.samplesID.txt"))
        ]
        if not missing:
            return  # all per-db files already present (the --what all path) -> NO-OP
        patterns = [f"{remote_prefix}/{ref}_{db}.samplesID.txt" for db in missing]
        staging = os.path.join(workdir, ".hf_stage_samplesids")
        shutil.rmtree(staging, ignore_errors=True)
        try:
            _hf_snapshot(repo, patterns, staging, token)
            os.makedirs(sdir, exist_ok=True)
            for db in missing:
                fname = f"{ref}_{db}.samplesID.txt"
                staged = os.path.join(staging, remote_prefix, fname)
                if os.path.isfile(staged):
                    # plain move, no decompress — same as the flat 'samples' path.
                    shutil.move(staged, os.path.join(sdir, fname))
                else:
                    # not published in the repo: warn + continue. The synthesizer
                    # then no-ops gracefully at its own missing-component guard.
                    sys.stdout.write(
                        f"NOTE: could not fetch per-dataset samplesID "
                        f"'{remote_prefix}/{fname}' from {repo}; combined samplesID "
                        f"may not be synthesized.\n"
                    )
        finally:
            shutil.rmtree(staging, ignore_errors=True)
    except Exception as exc:  # network/HF/dep error -> non-fatal, STDOUT only
        sys.stdout.write(
            f"NOTE: could not ensure per-dataset samplesID files for '{vcf_name}' "
            f"({exc}); combined samplesID may not be synthesized.\n"
        )


def download_component(
    component: str,
    workdir: str,
    repo: Optional[str] = None,
    ref: str = "hg38",
    dataset: Optional[str] = None,
    index_name: Optional[str] = None,
    token: Optional[str] = None,
    genotypes: bool = True,
) -> str:
    """Download one CRISPRme reference component from HF into the canonical
    local layout under ``workdir``.

    Args:
        component: one of genome|annotations|pams|samples|vcf|index.
        workdir: working directory the CRISPRme dir-tree lives under.
        repo: HF repo id (defaults via :func:`resolve_repo`).
        ref: reference-genome name (for ``genome``/``index``), e.g. "hg38".
        dataset: variant dataset name for ``vcf`` (e.g. "1000G", "HGDP").
        index_name: precomputed index directory name for ``index``
            (e.g. "NGG_2_hg38").

    Returns:
        The local destination directory populated by the download.
    """
    if component not in _COMPONENT_PREFIXES:
        raise ValueError(
            f"Unknown component '{component}'. Expected one of: "
            f"{', '.join(sorted(_COMPONENT_PREFIXES))}."
        )
    routed = component == "index" and bool(index_name) and not repo_overridden(repo)
    routed_from = None  # the wide-flank repo an index fell back from, if it did
    if routed:
        # Same routing as publish; the fallback probe waits until the index is known
        # not to be installed (below), so a complete install costs no network call.
        repo = resolve_index_repo(index_name, repo)
    else:
        repo = resolve_repo(repo)
    remote_prefix = _COMPONENT_PREFIXES[component]
    local_dir = os.path.join(workdir, _COMPONENT_LOCALDIR[component])

    # Web progress (env-gated via CRISPRME_PROGRESS_FILE; a pure no-op on the CLI /
    # the current happy path). Seed the known total up-front from the HF catalog so
    # the bar has a denominator from byte one; refined per-stage for the multi-part
    # index path below.
    reporter = _ProgressReporter()
    if reporter.enabled:
        _name = {"genome": ref, "vcf": dataset, "index": index_name}.get(component)
        if component in ("genome", "vcf", "index") and _name:
            reporter.set_total(_catalog_component_size(component, _name, repo))
        reporter.set_phase(f"Downloading {component}", _name or "", carry_done=False)

    if component == "genome":
        patterns = [f"{remote_prefix}/{ref}/*"]
        staging = os.path.join(workdir, ".hf_stage_genome")
        _hf_snapshot(repo, patterns, staging, token, reporter=reporter)
        src = os.path.join(staging, remote_prefix, ref)
        if not os.path.isdir(src):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No genome '{ref}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{ref}/'). Check --ref, or the repo "
                f"may not have this genome uploaded yet."
            )
        dest = os.path.join(local_dir, ref)
        files = sorted(os.listdir(src))
        # The move + gunzip below has no progress bar and, for a full genome
        # (~hundreds of per-contig FASTAs, several GB), runs silently for minutes
        # after "Fetching N files: 100%" — which reads as "stuck". Announce it.
        n_gz = sum(1 for f in files if f.endswith(".gz"))
        sys.stderr.write(
            f"Decompressing + staging {len(files)} file(s)"
            + (f" ({n_gz} gzip)" if n_gz else "")
            + f" into {dest} — no progress bar; a full genome can take a few minutes...\n"
        )
        sys.stderr.flush()

        def _build_genome(stage_dir):
            # decompress into the hidden staging dir; it is swapped into Genomes/<ref>/
            # only after this returns cleanly, so an interrupted/failed decompress
            # never leaves a half-populated genome folder in place.
            with _decompress_progress(
                reporter, stage_dir, phase="Decompressing genome"
            ):
                for fn in files:
                    moved = shutil.move(
                        os.path.join(src, fn), os.path.join(stage_dir, fn)
                    )
                    if moved.endswith(".gz"):
                        decompress_gz(moved)

        try:
            _atomic_install_dir(dest, _build_genome)
        finally:
            shutil.rmtree(staging, ignore_errors=True)
        return dest

    if component == "vcf":
        if not dataset:
            raise ValueError("--dataset is required to download a 'vcf' component")
        patterns = [f"{remote_prefix}/{dataset}/*"]
        staging = os.path.join(workdir, ".hf_stage_vcf")
        _hf_snapshot(repo, patterns, staging, token, reporter=reporter)
        src = os.path.join(staging, remote_prefix, dataset)
        if not os.path.isdir(src):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No raw VCF dataset '{dataset}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{dataset}/'). The default repo hosts "
                f"prebuilt variant-aware INDEXES (with their genotype stores for "
                f"sample-level analysis), not raw per-dataset VCFs. For variant-aware "
                f"search use:  crisprme.py download --what index --index-name "
                f"NRG_3_hg38+hg38_HGDP1kGP  (the genotype store it ships covers carriers "
                f"and per-sample output). To build your own index from raw VCFs, place "
                f"them under VCFs/{dataset}/ yourself, or use --hf-repo to point at a "
                f"repo that hosts '{remote_prefix}/{dataset}/'."
            )
        dest = os.path.join(local_dir, dataset)

        def _build_vcf(stage_dir):
            for fn in sorted(os.listdir(src)):
                shutil.move(os.path.join(src, fn), os.path.join(stage_dir, fn))  # keep bgzip

        try:
            _atomic_install_dir(dest, _build_vcf)
        finally:
            shutil.rmtree(staging, ignore_errors=True)
        return dest

    if component == "build-inputs":
        if not dataset:
            raise ValueError(
                "--dataset is required for 'build-inputs' (e.g. hg38_HGDP1kGP). This fetches "
                "the PRE-ENRICHED genome kit so build-index-only skips enrichment and only "
                "indexes (fast / PAM-swap path); use --what vcf to build from the raw source VCFs."
            )
        patterns = [f"{remote_prefix}/{dataset}.tar.gz", f"{remote_prefix}/{dataset}.SHA256SUMS"]
        staging = os.path.join(workdir, ".hf_stage_binputs")
        _hf_snapshot(repo, patterns, staging, token)
        tgz = os.path.join(staging, remote_prefix, f"{dataset}.tar.gz")
        if not os.path.isfile(tgz):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No build-inputs kit '{dataset}.tar.gz' in '{repo}' (looked for "
                f"'{remote_prefix}/{dataset}.tar.gz'). Pre-enriched kits exist for a SUBSET of "
                f"datasets; build from the raw source VCFs instead:  crisprme.py download "
                f"--what vcf --dataset {dataset}  then  build-index-only --vcf VCFs/{dataset} ..."
            )
        # Optional checksum verification (the kit ships a <dataset>.SHA256SUMS sidecar).
        sums = os.path.join(staging, remote_prefix, f"{dataset}.SHA256SUMS")
        if os.path.isfile(sums):
            import hashlib
            want = {}
            for line in open(sums):
                parts = line.split()
                if len(parts) >= 2:
                    want[os.path.basename(parts[-1])] = parts[0].lower()
            exp = want.get(f"{dataset}.tar.gz")
            if exp:
                h = hashlib.sha256()
                with open(tgz, "rb") as fh:
                    for blk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(blk)
                if h.hexdigest().lower() != exp:
                    shutil.rmtree(staging, ignore_errors=True)
                    raise ValueError(
                        f"build-inputs kit '{dataset}.tar.gz' SHA256 mismatch vs its "
                        f"SHA256SUMS — re-download (the transfer may be truncated)."
                    )
                sys.stderr.write(f"build-inputs: {dataset}.tar.gz SHA256 verified.\n")
        # Extract the enriched genome into Genomes/. Robust to the archive's internal
        # layout: unpack to a temp dir, then move every 'hg38+hg38_<ds>...' top dir into
        # Genomes/ (the kit may be rooted at that dir OR nested under Genomes/).
        # (tarfile is imported at module scope; do NOT re-import locally — a function-local
        # import would shadow it for the whole function and break the index/genome branches.)
        extract = os.path.join(staging, "_x")
        os.makedirs(extract, exist_ok=True)
        extract_real = os.path.realpath(extract)
        with tarfile.open(tgz, "r:gz") as tf:
            # Validate members before extracting: a build-inputs kit can come from a
            # user-chosen --hf-repo, so reject non-regular/dir entries (symlink/hardlink/
            # device) and any member that would resolve OUTSIDE the staging dir (path
            # traversal). Same guard as personal_assembly.py's bundle unpack.
            for m in tf.getmembers():
                if not (m.isreg() or m.isdir()):
                    raise ValueError(
                        f"build-inputs: refusing unsafe tar entry (not a regular file "
                        f"or directory): {m.name!r}")
                dest = os.path.realpath(os.path.join(extract, m.name))
                if dest != extract_real and not dest.startswith(extract_real + os.sep):
                    raise ValueError(
                        f"build-inputs: refusing tar entry outside the staging dir: {m.name!r}")
            tf.extractall(extract)
        os.makedirs(local_dir, exist_ok=True)
        moved = []
        for root, dirs, _files in os.walk(extract):
            for d in list(dirs):
                if d.startswith(ref + "+") or d.startswith("hg38+hg38_"):
                    src = os.path.join(root, d)
                    dst = os.path.join(local_dir, d)
                    if os.path.exists(dst):
                        shutil.rmtree(dst, ignore_errors=True)
                    shutil.move(src, dst)
                    moved.append(d)
                    dirs.remove(d)  # don't descend into a moved dir
        shutil.rmtree(staging, ignore_errors=True)
        if not moved:
            raise ValueError(
                f"build-inputs kit '{dataset}.tar.gz' unpacked but contained no enriched-genome "
                f"dir (expected 'hg38+hg38_{dataset}/'). Report this kit as malformed."
            )
        sys.stderr.write(
            f"build-inputs: staged {', '.join(moved)} into {local_dir} — build-index-only will "
            f"now SKIP enrichment and only index. Pair with the PAM-independent Dictionaries "
            f"(download --what index --index-name <...> for the shipped index, or build-index-only "
            f"--samplesID to regenerate them).\n"
        )
        return local_dir

    if component == "index":
        if not index_name:
            raise ValueError("--index-name is required to download an 'index' component")
        # Already installed -> NO-OP, mirroring the guard the variant path already
        # applies to its companion reference-index fetch below. Checked BEFORE the
        # dict-less probe so a COMPLETE install costs no network call at all.
        # canonical_index_name() is a pure string function and strips a "-dictless"
        # publish marker, so this resolves to the same folder whichever name the
        # caller passed.
        #
        # Skip ONLY when the install is COMPLETE, so re-running `download --what
        # index` still HEALS a partial/interrupted install (the documented recovery
        # for "no variant off-targets / Samples degraded") rather than no-opping over
        # it. A reference-only index has no companions. A VARIANT index is complete
        # only with ALL of: its _INDELS companion (renamed into place by a SEPARATE
        # os.rename, so a main dir without it is an interrupted install); the shared
        # reference index (<pam>_<N>_<ref>, which travels once -- not inside each
        # variant tarball); and, unless --no-genotypes, its Tier-1 genotype store
        # (Dictionaries/genotypes_<vcf>/). Any missing piece falls through to the
        # fetch/heal path below, which re-fetches only what is absent.
        installed_name = canonical_index_name(index_name)
        installed_dir = os.path.join(local_dir, installed_name)
        _vcf = installed_name.partition("+")[2]
        _ref_index_dir = os.path.join(local_dir, installed_name.partition("+")[0])
        _genotype_dir = os.path.join(workdir, "Dictionaries", f"genotypes_{_vcf}")
        _complete = os.path.isdir(installed_dir) and (
            not _vcf
            or (
                os.path.isdir(f"{installed_dir}_INDELS")
                and os.path.isdir(_ref_index_dir)
                and (not genotypes or os.path.isdir(_genotype_dir))
            )
        )
        if _complete:
            print(
                f"Index {installed_name} is already installed -> {installed_dir} "
                f"(remove that folder to re-download)",
                flush=True,
            )
            return installed_dir
        # Fall back to the default repo if the routed one does not have it, so an
        # index published before the split still resolves.
        if routed and repo != DEFAULT_HF_REPO and not _index_archive_exists(repo, index_name):
            routed_from, repo = repo, DEFAULT_HF_REPO
        # Dict-less is the default in 2.4.0. If the caller asked for a VARIANT
        # index by its bare (legacy dict-based) name but a corrected `-dictless`
        # companion exists on the repo, prefer it: the `-dictless` marker is
        # stripped at install time (canonical_index_name), so both names install
        # to the SAME search-resolvable folder, and the dict-less tarball is the
        # smaller, combined-AF one that 2.4.0 ships. Best-effort — any lookup
        # failure falls through to the requested name (legacy behaviour), so this
        # never makes a previously-working download fail.
        if "+" in index_name and "-dictless" not in index_name:
            dictless_name = index_name.replace("+", "-dictless+", 1)
            try:
                if _require_hf().HfApi().file_exists(
                    repo,
                    f"{remote_prefix}/{dictless_name}.tar.gz",
                    repo_type="dataset",
                    token=token,
                ):
                    sys.stderr.write(
                        f"Preferring dict-less variant index '{dictless_name}' "
                        f"(smaller, corrected combined-AF) over legacy "
                        f"'{index_name}'; both install to the same folder.\n"
                    )
                    index_name = dictless_name
            except Exception:
                pass  # existence check failed — use the requested name as-is
        patterns = [f"{remote_prefix}/{index_name}.tar.gz"]
        staging = os.path.join(workdir, ".hf_stage_index")
        # multi-part phase labels: a VARIANT index installs in 3 stages (main index
        # tarball -> genotypes companion -> shared reference index); a reference-only
        # index is a single download. Re-seed the total to the ACTUAL (possibly
        # -dictless) tarball's size so the bar denominator matches what is fetched.
        _is_variant = "+" in index_name
        if reporter.enabled:
            reporter.set_total(
                _catalog_component_size("index", index_name, repo)
                or reporter.bytes_total
            )
            reporter.set_phase(
                "Downloading index 1/3" if _is_variant else "Downloading index",
                index_name,
                carry_done=False,
            )
        _hf_snapshot(repo, patterns, staging, token, reporter=reporter)
        tarball = os.path.join(staging, remote_prefix, f"{index_name}.tar.gz")
        if not os.path.isfile(tarball):
            shutil.rmtree(staging, ignore_errors=True)
            where = (
                f"neither in HuggingFace repo '{routed_from}' (where indexes with a "
                f"flank wider than 3 nt are published) nor in the fallback '{repo}'"
                if routed_from else f"in HuggingFace repo '{repo}'"
            )
            raise ValueError(
                f"No index '{index_name}' found {where} "
                f"(looked for '{remote_prefix}/{index_name}.tar.gz'). Check "
                f"--index-name, or build it locally with 'crisprme.py "
                f"build-index-only' — the repo may not have it uploaded yet."
                + (f" If '{routed_from}' was only unreachable, retry, or name it "
                   f"with --hf-repo {routed_from}." if routed_from else "")
            )
        os.makedirs(local_dir, exist_ok=True)
        # GAP 3: the remote tarball filename AND the archive's internal top-level
        # dir are the ORIGINAL index_name (that is how it was published), so the
        # fetch (above) + the extract-validation below keep using index_name. But
        # the INSTALL dir name must be search-resolvable: strip a publish-only
        # marker (e.g. "-dictless") from the ref segment so the search convention
        # (<pam>_<N>_<ref>+<vcf>, ref == genome-folder basename) resolves it. This
        # is a NO-OP for a canonical name (install_name == index_name).
        install_name = canonical_index_name(index_name)
        # Extract into a HIDDEN staging dir on the SAME filesystem as
        # genome_library, validate, then atomically rename into place. This keeps
        # the install atomic: a partial/failed extract never leaves a discoverable
        # <index_name>/ dir (get_available_indexes skips dot-prefixed dirs), so the
        # UI only ever sees a complete index.
        extract_tmp = os.path.join(local_dir, f".extract_{index_name}")
        shutil.rmtree(extract_tmp, ignore_errors=True)
        os.makedirs(extract_tmp)
        try:
            label = None
            with tarfile.open(tarball) as tf:
                # surface the provenance manifest (if present) without unpacking it
                # alongside the index; extract only the index folder members
                members = [m for m in tf.getmembers() if m.name != "manifest.json"]
                with _decompress_progress(
                    reporter, extract_tmp, phase="Extracting index"
                ):
                    tf.extractall(extract_tmp, members=members)
                try:
                    mf = tf.extractfile("manifest.json")
                    if mf is not None:
                        meta = json.load(mf)
                        sys.stderr.write(
                            f"Downloaded index {index_name} "
                            f"(built {meta.get('created_at', 'unknown')})\n"
                        )
                        label = meta.get("display_label")
                except (KeyError, json.JSONDecodeError):
                    pass  # no/invalid manifest — the index itself is what matters
            # validate: the archive must contain the expected index folder
            produced = os.path.join(extract_tmp, index_name)
            if not os.path.isdir(produced) or not os.listdir(produced):
                raise ValueError(
                    f"Downloaded archive for '{index_name}' did not contain the "
                    f"expected index folder — nothing installed."
                )
            # a variant index (name has '+') must ship a non-empty _INDELS
            # companion; a truncated/corrupt archive missing it would otherwise
            # install a half-index whose indel search silently fails.
            if "+" in index_name:
                indels_produced = os.path.join(extract_tmp, index_name + "_INDELS")
                if not os.path.isdir(indels_produced) or not os.listdir(indels_produced):
                    raise ValueError(
                        f"Downloaded archive for variant index '{index_name}' is "
                        f"missing its _INDELS companion (incomplete/corrupt download) "
                        f"— nothing installed."
                    )
            # restore the friendly display name into the staged index before swap
            if label:
                try:
                    with open(os.path.join(produced, ".display_label"), "w") as lf:
                        lf.write(label)
                except OSError:
                    pass
            # GAP 3 diagnostic: announce a rename (STDOUT only; stderr is fatal
            # downstream). Never fires for a canonical name (NO-OP).
            if install_name != index_name:
                sys.stdout.write(
                    f"Installing index '{index_name}' under search-resolvable "
                    f"name '{install_name}' (stripped publish marker from ref "
                    f"segment)\n"
                )
            # atomically swap the index and its _INDELS companion into place. The
            # SOURCE dir in extract_tmp is the ORIGINAL name (that is what the
            # tarball contains); the DESTINATION is the canonical name so the
            # search can resolve it. For a canonical name these are identical, so
            # the paths (and the backup name) are byte-for-byte today's behavior.
            for sub in (index_name, index_name + "_INDELS"):
                src = os.path.join(extract_tmp, sub)
                if not os.path.isdir(src):
                    continue  # reference-only indexes have no _INDELS companion
                dst_name = canonical_index_name(sub)
                dst = os.path.join(local_dir, dst_name)
                backup = os.path.join(local_dir, f".{dst_name}.replaced")  # hidden -> unlisted
                # roll forward an interrupted prior swap: if a backup exists but the
                # live dir is gone, a previous run died between the two renames below
                # -> restore it before we would otherwise delete it (avoids losing
                # the only surviving copy of the index).
                if os.path.isdir(backup) and not os.path.exists(dst):
                    os.rename(backup, dst)
                shutil.rmtree(backup, ignore_errors=True)
                if os.path.exists(dst):
                    os.rename(dst, backup)  # move any existing index aside (same FS)
                os.rename(src, dst)  # move the new one into place (atomic, same FS)
                shutil.rmtree(backup, ignore_errors=True)
            # a variant index bundles its per-sample dictionaries (see publish_index):
            # install them into <workdir>/Dictionaries/ so the variant search finds
            # them without the source VCFs (same atomic swap as the index dirs above).
            staged_dicts = os.path.join(extract_tmp, "Dictionaries")
            if os.path.isdir(staged_dicts):
                dst_dicts_root = os.path.join(workdir, "Dictionaries")
                os.makedirs(dst_dicts_root, exist_ok=True)
                for sub in os.listdir(staged_dicts):
                    src = os.path.join(staged_dicts, sub)
                    if not os.path.isdir(src):
                        continue
                    dst = os.path.join(dst_dicts_root, sub)
                    backup = os.path.join(dst_dicts_root, f".{sub}.replaced")
                    if os.path.isdir(backup) and not os.path.exists(dst):
                        os.rename(backup, dst)
                    shutil.rmtree(backup, ignore_errors=True)
                    if os.path.exists(dst):
                        os.rename(dst, backup)
                    os.rename(src, dst)
                    shutil.rmtree(backup, ignore_errors=True)
            # ADDITIVE: a self-complete variant index (published after the samplesID
            # emit fix) bundles its samplesID lists under samplesIDs/ (see
            # publish_index / _make_index_tarball). Install them into
            # <workdir>/samplesIDs/ so `download --what index` searches WITHOUT a
            # separate `--what samples`. INSTALL-ONLY-IF-ABSENT: never clobber a
            # user's / `--what all`-fetched samplesID (a stale bundled file must not
            # shadow a freshly fetched per-db list). This lands BEFORE the
            # _ensure_perdb_samplesids + synthesize_combined_samplesid calls below,
            # which then strictly no-op. Tiny text files -> a plain guarded copy is
            # enough (no atomic swap needed). Absent for a legacy/reference index ->
            # no-op, download-side 2.3.1 fallback still covers a merged index.
            staged_sids = os.path.join(extract_tmp, "samplesIDs")
            if os.path.isdir(staged_sids):
                dst_sids_root = os.path.join(workdir, "samplesIDs")
                os.makedirs(dst_sids_root, exist_ok=True)
                for fn in sorted(os.listdir(staged_sids)):
                    src = os.path.join(staged_sids, fn)
                    if not os.path.isfile(src):
                        continue
                    dst = os.path.join(dst_sids_root, fn)
                    if os.path.isfile(dst):
                        continue  # never overwrite a shipped/fresher samplesID
                    shutil.move(src, dst)
        finally:
            shutil.rmtree(extract_tmp, ignore_errors=True)
            shutil.rmtree(staging, ignore_errors=True)
        # ADDITIVE Tier-1: after the main tarball (which may now carry
        # registry_<vcf>/), ALSO fetch + extract the separate genotypes_<vcf>.tar.gz
        # into <workdir>/Dictionaries/ so genotypes_<vcf>/ lands as a sibling of
        # dictionaries_/registry_ — the exact path _resolve_genotype_paths reads.
        # Default ON; --no-genotypes (genotypes=False) skips it (detection-only,
        # degraded Samples). A missing companion (registry-only publish) is NOT an
        # error: warn to STDOUT and continue — search falls back gracefully. A
        # dictless index (no dictionaries_<vcf>/ in the main tarball) is also FINE:
        # the registry + genotype tiers cover it.
        vcf_name = index_name.partition("+")[2]
        if vcf_name:
            if not genotypes:
                sys.stdout.write(
                    f"Skipping Tier-1 genotype store for '{index_name}' "
                    f"(--no-genotypes): off-target detection works, but per-sample "
                    f"Samples will be degraded until genotypes_{vcf_name}/ is present.\n"
                )
            else:
                if reporter.enabled:
                    reporter.set_total(
                        _catalog_component_size(
                            "index", f"genotypes_{vcf_name}", repo
                        )
                    )
                    reporter.set_phase(
                        "Downloading genotypes 2/3", f"genotypes_{vcf_name}",
                        carry_done=False,
                    )
                try:
                    got = _fetch_genotypes_companion(
                        repo, workdir, vcf_name, remote_prefix, token,
                        reporter=reporter,
                    )
                except Exception as exc:  # network/extract error -> non-fatal
                    got = None  # None = errored and already reported below
                    sys.stdout.write(
                        f"NOTE: could not fetch Tier-1 genotype store "
                        f"'genotypes_{vcf_name}.tar.gz' ({exc}); continuing without "
                        f"it — off-target detection works, Samples degraded.\n"
                    )
                if got:
                    sys.stdout.write(
                        f"Installed Tier-1 genotype store genotypes_{vcf_name}/ "
                        f"into {os.path.join(workdir, 'Dictionaries')}\n"
                    )
                elif got is False:
                    # Genuinely absent (a registry-only publish), as distinct from a
                    # fetch that errored -- that case reported itself above, and
                    # saying "has no companion" for it would be wrong.
                    sys.stdout.write(
                        f"NOTE: index '{index_name}' has no genotypes_{vcf_name}.tar.gz "
                        f"companion in {repo} (registry-only/detection publish); "
                        f"off-target detection works, per-sample Samples degraded.\n"
                    )
            # ALWAYS-SHIP the shared reference index. A variant index tarball
            # (<pam>_<N>_<ref>+<vcf>) does NOT bundle the reference index
            # (<pam>_<N>_<ref>) -- it is common to every variant dataset, so it
            # travels once, not inside each variant archive. Fetch it here so a
            # fresh download is IMMEDIATELY searchable: the search needs BOTH the
            # variant AND the reference index (the latter for reference off-targets).
            # No-op if it is already installed. Non-fatal on failure -- the search
            # falls back to building it on demand from the raw genome.
            ref_index_name = install_name.partition("+")[0]  # <pam>_<N>_<ref>
            if ref_index_name and not os.path.isdir(
                os.path.join(local_dir, ref_index_name)
            ):
                if reporter.enabled:
                    reporter.set_total(
                        _catalog_component_size("index", ref_index_name, repo)
                    )
                    reporter.set_phase(
                        "Downloading reference 3/3", ref_index_name,
                        carry_done=False,
                    )
                try:
                    # the recursive call builds its own reporter from the same env
                    # progress file, so it keeps updating the bar for this stage.
                    download_component(
                        "index", workdir, repo=repo, ref=ref,
                        index_name=ref_index_name, token=token, genotypes=False,
                    )
                    sys.stdout.write(
                        f"Also installed the shared reference index "
                        f"'{ref_index_name}' (needed for reference off-targets).\n"
                    )
                except Exception as exc:  # non-fatal: build-on-demand fallback exists
                    sys.stdout.write(
                        f"NOTE: could not fetch the reference index "
                        f"'{ref_index_name}' ({exc}); the search will build it on "
                        f"demand from the raw genome (a shared base component).\n"
                    )
        # GAP 2: a merged variant index searches ONE vcf folder and the search
        # expects a single combined samplesIDs/<vcf_name>.samplesID.txt, but the
        # published artifact ships only the per-db lists. Synthesize the combined
        # file from those when it is missing (NO-OP when already present / single
        # dataset / a component is absent). Non-fatal, STDOUT-only (stderr fatal).
        # vcf_name is unaffected by GAP 3's canonicalization (marker is only in the
        # base segment), so the file is named for the canonical installed index.
        if vcf_name:
            # GAP 2b: a standalone `download --what index` never runs the 'samples'
            # component, so the per-db samplesID lists the synthesizer unions may be
            # absent (whereas `--what all` fetches them first). Ensure the per-db
            # files for THIS index's datasets are present — fetching any missing from
            # HF's samplesIDs/ component — BEFORE synthesizing. Internally guarded;
            # STDOUT-only; strict NO-OP when they are already present (`--what all`).
            try:
                _ensure_perdb_samplesids(repo, workdir, vcf_name, ref, token)
            except Exception as exc:  # belt-and-suspenders (helper is guarded too)
                sys.stdout.write(
                    f"NOTE: could not ensure per-dataset samplesID files for "
                    f"'{vcf_name}' ({exc}); combined samplesID may not be "
                    f"synthesized.\n"
                )
            try:
                _sid = synthesize_combined_samplesid(workdir, vcf_name, ref=ref)
            except OSError as exc:
                _sid = None
                sys.stdout.write(
                    f"NOTE: could not synthesize combined samplesID for "
                    f"'{vcf_name}' ({exc}); provide "
                    f"samplesIDs/{vcf_name}.samplesID.txt manually if needed.\n"
                )
            if _sid:
                sys.stdout.write(
                    f"Generated combined samplesID {os.path.basename(_sid)} "
                    f"from per-dataset lists in "
                    f"{os.path.join(workdir, 'samplesIDs')}\n"
                )
            # Make the index searchable from the CLI out of the box: name this
            # dataset + its samples in the --vcf / --samplesID list files the web
            # builds per-search but a CLI complete-search reads from disk.
            _ensure_search_lists(workdir, vcf_name)
        return os.path.join(local_dir, install_name)

    # flat components: annotations, pams, samples
    patterns = [f"{remote_prefix}/*"]
    staging = os.path.join(workdir, f".hf_stage_{component}")
    _hf_snapshot(repo, patterns, staging, token)
    src = os.path.join(staging, remote_prefix)
    if not os.path.isdir(src):
        shutil.rmtree(staging, ignore_errors=True)
        raise ValueError(
            f"No '{component}' files found in HuggingFace repo '{repo}' "
            f"(looked for '{remote_prefix}/'). The repo may not have this "
            f"component uploaded yet."
        )
    os.makedirs(local_dir, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        moved = shutil.move(os.path.join(src, fn), os.path.join(local_dir, fn))
        if component == "annotations" and moved.endswith(".gz") and not moved.endswith(".bed.gz"):
            decompress_gz(moved)
    shutil.rmtree(staging, ignore_errors=True)
    return local_dir


def _make_index_tarball(
    tarball: str,
    index_dir: str,
    index_name: str,
    indels_dir: str,
    manifest: Dict,
    dict_dirs: Optional[List[str]] = None,
    samplesid_files: Optional[List[str]] = None,
) -> None:
    """Builds the one-file-per-index ``.tar.gz`` (index + _INDELS + dicts +
    samplesID lists + manifest).

    Prefers ``pigz`` (parallel gzip) via GNU ``tar`` so a large index compresses
    across all cores instead of one — a 170GB pamless index tars in minutes with
    pigz vs ~a day with single-threaded gzip. Falls back to Python's (single-
    threaded) ``tarfile`` when ``pigz``/``tar`` are unavailable. Both paths write
    an identical archive layout: ``<index_name>/``, optional ``<index_name>_INDELS/``,
    each ``dict_dirs`` entry stored under ``Dictionaries/<basename>/`` (the per-sample
    variant dictionaries a variant search needs at post-analysis time — bundling them
    lets a downloaded index search WITHOUT the source VCFs), each ``samplesid_files``
    entry stored under ``samplesIDs/<basename>`` (the combined + per-db samplesID
    lists a MERGED variant search's ``--samplesID`` needs — bundling them makes
    ``download --what index`` self-complete, no separate ``--what samples``), and
    ``manifest.json`` at the archive root (see download_component's extractor).

    Both the samplesID files (small text) and manifest.json are staged into a
    single temp dir so the two code paths share ONE staging mechanism.
    """
    manifest_bytes = json.dumps(manifest, indent=2).encode()
    pigz = shutil.which("pigz")
    tar = shutil.which("tar")
    # stage the small text artifacts (manifest.json + samplesID lists) in a temp
    # dir so BOTH paths add them by the same arcnames the extractor expects:
    #   <tmpd>/manifest.json          -> manifest.json (archive root)
    #   <tmpd>/samplesIDs/<basename>  -> samplesIDs/<basename>
    tmpd = tempfile.mkdtemp(prefix="hfpub_")
    try:
        with open(os.path.join(tmpd, "manifest.json"), "wb") as mf:
            mf.write(manifest_bytes)
        staged_sids = []  # arcnames (relative to tmpd) actually staged
        for sf in (samplesid_files or []):
            if not (os.path.isfile(sf) and os.path.getsize(sf) > 0):
                continue  # skip absent/empty (legacy build) -> fewer bundled files
            sdir = os.path.join(tmpd, "samplesIDs")
            os.makedirs(sdir, exist_ok=True)
            arc = os.path.join("samplesIDs", os.path.basename(sf))
            shutil.copy2(sf, os.path.join(tmpd, arc))
            staged_sids.append(arc)
        if pigz and tar:
            try:
                parent = os.path.dirname(index_dir)
                # -C <dir> <name> stores <name> under its basename == the arcname
                inputs = ["-C", parent, index_name]
                if os.path.isdir(indels_dir):
                    inputs += ["-C", parent, os.path.basename(indels_dir)]
                for dp in (dict_dirs or []):
                    # store as Dictionaries/<basename>/... so the extractor routes it
                    # to <workdir>/Dictionaries/ (not genome_library/)
                    inputs += ["-C", os.path.dirname(os.path.dirname(dp)),
                               os.path.join("Dictionaries", os.path.basename(dp))]
                # samplesID lists + manifest.json both come from the staging tmpd
                for arc in staged_sids:
                    inputs += ["-C", tmpd, arc]
                inputs += ["-C", tmpd, "manifest.json"]
                subprocess.check_call(
                    [tar, "--use-compress-program", pigz, "-cf", tarball, *inputs]
                )
                return
            except (subprocess.CalledProcessError, OSError) as exc:
                # this tar may not support --use-compress-program (busybox / very
                # old), or pigz died mid-stream -> drop any partial tarball and fall
                # through to the single-threaded Python path so publish succeeds.
                sys.stderr.write(
                    f"Note: parallel (pigz) compression failed ({exc}); falling back "
                    f"to single-threaded gzip.\n"
                )
                if os.path.exists(tarball):
                    os.remove(tarball)
        else:
            sys.stderr.write(
                "Note: pigz/tar not found — compressing with single-threaded gzip "
                "(slow for large indexes). Install pigz for parallel compression.\n"
            )
        # fallback: single-threaded Python tarfile (identical layout)
        with tarfile.open(tarball, "w:gz") as tf:
            tf.add(index_dir, arcname=index_name)  # -> genome_library/<name>/ on unpack
            if os.path.isdir(indels_dir):
                tf.add(indels_dir, arcname=os.path.basename(indels_dir))
            for dp in (dict_dirs or []):  # -> Dictionaries/<basename>/ on unpack
                tf.add(dp, arcname=os.path.join("Dictionaries", os.path.basename(dp)))
            for arc in staged_sids:  # -> samplesIDs/<basename> on unpack
                tf.add(os.path.join(tmpd, arc), arcname=arc)
            info = tarfile.TarInfo(name="manifest.json")
            info.size = len(manifest_bytes)
            tf.addfile(info, io.BytesIO(manifest_bytes))
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def _make_genotypes_tarball(tarball: str, genotypes_dir: str) -> None:
    """Build the SEPARATE Tier-1 genotype-store archive ``genotypes_<vcf>.tar.gz``.

    Layout mirrors ``_make_index_tarball``'s dict routing: the single member is
    ``Dictionaries/<basename>/...`` so the download extractor drops it straight
    into ``<workdir>/Dictionaries/`` alongside ``dictionaries_``/``registry_`` —
    the exact sibling path ``_resolve_genotype_paths`` reads. Prefers pigz/tar,
    falls back to single-threaded Python tarfile (identical archive layout). This
    is the big (~22GB genome) optional artifact, uploaded separately from the main
    index tarball so a detection-only install can skip it (``--no-genotypes``)."""
    arcname = os.path.join("Dictionaries", os.path.basename(genotypes_dir))
    parent_of_dictionaries = os.path.dirname(os.path.dirname(genotypes_dir))
    pigz = shutil.which("pigz")
    tar = shutil.which("tar")
    if pigz and tar:
        try:
            subprocess.check_call(
                [tar, "--use-compress-program", pigz, "-cf", tarball,
                 "-C", parent_of_dictionaries, arcname]
            )
            return
        except (subprocess.CalledProcessError, OSError) as exc:
            # diagnostics to STDOUT (stderr is fatal to the post-analysis pipeline)
            sys.stdout.write(
                f"Note: parallel (pigz) compression of genotype store failed "
                f"({exc}); falling back to single-threaded gzip.\n"
            )
            if os.path.exists(tarball):
                os.remove(tarball)
    else:
        sys.stdout.write(
            "Note: pigz/tar not found — compressing the genotype store with "
            "single-threaded gzip (slow; ~22GB genome). Install pigz for parallel "
            "compression.\n"
        )
    with tarfile.open(tarball, "w:gz") as tf:
        tf.add(genotypes_dir, arcname=arcname)


def publish_index(
    index_dir: str,
    repo: Optional[str] = None,
    token: Optional[str] = None,
    display_name: Optional[str] = None,
    dictless: bool = False,
) -> str:
    """Tar a locally built genome_library index and upload it to HF.

    Args:
        index_dir: path to the ``genome_library/<index_name>`` directory.
        repo: HF repo id (defaults via :func:`resolve_repo`).
        token: HF write token (defaults via :func:`resolve_token`).
        display_name: optional human-friendly label for the index.
        dictless: when ``True``, EXCLUDE the per-sample SNP dictionaries
            (``dictionaries_<vcf>/``) from the main tarball — the additive
            registry (Tier-0) + genotype (Tier-1) tiers replace them — while
            KEEPING ``log_indels_<vcf>/`` (indel post-analysis still needs the
            indel logs; the tiers are SNP-only). When ``False`` (default),
            publishing is BYTE-FOR-BYTE unchanged from the classic path: the
            SNP dicts are bundled exactly as before; the registry tier, if
            present, is simply ADDED.

    In BOTH modes, the small Tier-0 ``registry_<vcf>/`` dir is added to the main
    tarball when it exists (it powers out-of-the-box off-target detection +
    corrected AF/rsID), and — when a ``genotypes_<vcf>/`` dir exists — a SEPARATE
    ``genotypes_<vcf>.tar.gz`` companion is produced and uploaded to the same repo
    under the same ``indexes/`` prefix (the big optional Tier-1 artifact).

    SECONDARY (future builds): the download side now synthesizes the combined
    ``samplesIDs/<vcf_name>.samplesID.txt`` for a MERGED index when it is missing
    (see :func:`synthesize_combined_samplesid`), so the already-published artifact
    works with zero manual steps. The cleaner long-term fix is to EMIT that
    combined file at build time (build-index-only already has the ordered
    ``{db: samplesID_path}`` map for the panel) and upload it under ``samplesIDs/``
    (component ``samples``); the download-side generation then becomes a pure
    fallback. This publish path is unchanged; the note is intentional.

    Returns:
        The remote path (``indexes/<index_name>.tar.gz``) the main index was
        uploaded to.
    """
    index_dir = os.path.abspath(index_dir.rstrip("/"))
    if not os.path.isdir(index_dir):
        raise ValueError(f"Index directory {index_dir} does not exist")
    index_name = os.path.basename(index_dir)
    # Route a wide-flank PAM-less index to its own repo, so an older client that
    # only knows the default repo cannot download an index it would misread.
    repo = resolve_index_repo(index_name, repo)
    token = resolve_token(token)
    if not token:
        raise ValueError(
            "An HF write token is required to publish an index. Provide --token "
            "or set HF_TOKEN (never commit it)."
        )
    hf = _require_hf()
    # self-describing provenance travels inside the tarball (adopted from the
    # #141 design): PAM / bulge / genome parsed from the folder name
    # <PAM>_<bMax+1>_<ref>, plus a UTC timestamp.
    manifest = {"name": index_name, "created_at": datetime.now(timezone.utc).isoformat()}
    parts = index_name.rsplit("_", 2)
    if len(parts) == 3:
        manifest["pam"], manifest["index_bmax"], manifest["genome"] = (
            parts[0],
            parts[1],
            parts[2],
        )
    # optional human-friendly name: explicit arg wins, else a .display_label sidecar
    # written at build time; travels in the manifest so download can restore it.
    if not display_name:
        sidecar = os.path.join(index_dir, ".display_label")
        if os.path.isfile(sidecar):
            try:
                display_name = open(sidecar).read().strip()
            except OSError:
                display_name = None
    if display_name:
        manifest["display_label"] = display_name
    # ONE tarball per index: bundle the indel-genome companion (<name>_INDELS) into
    # the same archive so an index is a single file on HF (indels are part of the
    # index, not a separate artifact). download extracts both dirs into genome_library.
    indels_dir = index_dir + "_INDELS"
    # A variant index (<ref>+<vcf>) needs its per-sample dictionaries at SEARCH time
    # (SNP + indel post-analysis). Bundle them so `download --what index` yields a
    # self-sufficient variant index that searches WITHOUT the multi-GB source VCFs —
    # submit_job skips genome enrichment when the precomputed index + these dicts are
    # present. Reference indexes (no '+') have no dictionaries.
    vcf_name = index_name.partition("+")[2]
    dict_dirs = []
    genotypes_dir = None
    if vcf_name:
        cwd = os.path.dirname(os.path.dirname(index_dir))  # .../genome_library/<name> -> workdir
        dicts_root = os.path.join(cwd, "Dictionaries")
        snp_name = f"dictionaries_{vcf_name}"
        indel_name = f"log_indels_{vcf_name}"
        registry_name = f"registry_{vcf_name}"
        genotypes_name = f"genotypes_{vcf_name}"
        indel_gt_name = f"indel_genotypes_{vcf_name}"
        # DEFAULT (classic) main-tarball members: SNP dicts + indel logs, exactly
        # as before. DICTLESS: drop the (152GB) per-sample SNP dicts but KEEP the
        # indel logs (the tiers are SNP-only). Missing-dir handling is unchanged:
        # a dir is bundled only when present + non-empty.
        wanted = [indel_name] if dictless else [snp_name, indel_name]
        for d in wanted:
            p = os.path.join(dicts_root, d)
            if os.path.isdir(p) and os.listdir(p):
                dict_dirs.append(p)
        missing = set(wanted) - {os.path.basename(p) for p in dict_dirs}
        if missing:
            sys.stderr.write(
                f"WARNING: variant index '{index_name}' is missing dictionaries "
                f"{sorted(missing)} under {dicts_root}; the "
                f"published index will NOT be searchable without the source VCFs. "
                f"Build the index in the same working dir so the dicts are present.\n"
            )
        # ADDITIVE Tier-0: always include the small registry_<vcf>/ in the MAIN
        # tarball when it exists (out-of-the-box off-target detection + AF/rsID).
        # Absent => behave exactly as today (no error). It is routed under
        # Dictionaries/ by _make_index_tarball just like the other dict dirs.
        reg_p = os.path.join(dicts_root, registry_name)
        if os.path.isdir(reg_p) and os.listdir(reg_p):
            dict_dirs.append(reg_p)
            manifest["has_registry"] = True
        # ADDITIVE (mega / aggregate index): bundle the sites-only INDEL AF sidecar
        # (indel_af_<vcf>/) in the MAIN tarball so indel off-targets from an aggregate
        # all-source panel carry per-dataset AF on download -- exactly like the SNP
        # registry does for SNPs (analisi_indels_NNN reads it from Dictionaries/ on
        # install). Small gzipped TSVs; absent (genotyped index) => no-op.
        indel_af_p = os.path.join(dicts_root, f"indel_af_{vcf_name}")
        if os.path.isdir(indel_af_p) and os.listdir(indel_af_p):
            dict_dirs.append(indel_af_p)
            manifest["has_indel_af"] = True
        # ADDITIVE (SNP+indel feature): bundle the phased indel-genotype store in the
        # MAIN tarball so CONFIRMED-cis SNP+indel co-occurrence survives download
        # (indel_genotypes_<vcf>/ -> Dictionaries/ on install, where analisi_indels
        # reads it). It is << the per-sample SNP genotype store; absent (feature-off
        # index) => no-op. Already-gzipped .tsv.gz, so it adds ~its on-disk size.
        indel_gt_p = os.path.join(dicts_root, indel_gt_name)
        if os.path.isdir(indel_gt_p) and os.listdir(indel_gt_p):
            dict_dirs.append(indel_gt_p)
            manifest["has_indel_genotypes"] = True
        # ADDITIVE Tier-1: the big genotype store travels as a SEPARATE tarball
        # (see below), never inside the main archive.
        gt_p = os.path.join(dicts_root, genotypes_name)
        if os.path.isdir(gt_p) and os.listdir(gt_p):
            genotypes_dir = gt_p
    # ADDITIVE samplesID self-completeness: bundle the samplesID lists this variant
    # index's datasets need — the combined <vcf_name>.samplesID.txt (emitted at
    # build time by crisprme.py._emit_combined_samplesid) plus the per-db
    # <ref>_<db>.samplesID.txt lists — INTO the main tarball, routed to
    # <workdir>/samplesIDs/ on download. This makes `download --what index` yield a
    # searchable install with NO separate `--what samples`/`--what all`; download's
    # 2.3.1 _ensure_perdb_samplesids + synthesize_combined_samplesid then no-op
    # (belt-and-suspenders). Each file is guarded by exists+non-empty; a legacy
    # index built before the emit fix simply bundles fewer files (or none) and the
    # download-side fallback still covers it. Reference-only indexes (no vcf_name)
    # skip this block entirely -> byte-for-byte-unchanged tarball.
    samplesid_files: List[str] = []
    if vcf_name:
        try:
            # ref token for the strict <ref>_<db>.samplesID.txt convention: parsed
            # from the index-dir name's REF segment, CANONICALIZED to drop a
            # publish-only marker (e.g. '-dictless'). NOTE: manifest["genome"] is a
            # name-split artifact (index_name.rsplit('_',2)) that is only meaningful
            # for a REFERENCE index; for a merged variant index it is the last VCF
            # token (e.g. 'HGDP'), NOT the ref -- so we do NOT use it here. Canonicalize
            # the whole name then re-parse the ref, so the bundled per-db files are
            # named 'hg38_<db>.samplesID.txt' and download (--ref hg38) resolves them.
            _canon = canonical_index_name(index_name)
            _canon_base = _canon.partition("+")[0]  # <pam>_<N>_<ref>
            _cparts = _canon_base.rsplit("_", 2)
            ref_token = _cparts[2] if len(_cparts) == 3 else manifest.get("genome", "hg38")
            sdir = os.path.join(cwd, "samplesIDs")
            candidates = [os.path.join(sdir, f"{vcf_name}.samplesID.txt")]  # combined
            for db in _derive_perdb_datasets(vcf_name, ref_token):  # [] => single ds
                candidates.append(os.path.join(sdir, f"{ref_token}_{db}.samplesID.txt"))
            for c in candidates:
                if os.path.isfile(c) and os.path.getsize(c) > 0:
                    samplesid_files.append(c)
        except Exception as exc:  # additive convenience -> never fail a publish
            sys.stdout.write(
                f"NOTE: could not collect samplesID lists for '{index_name}' "
                f"({exc}); publishing without bundled samplesIDs (download-side "
                f"synthesis still covers a merged index if the per-db lists are on "
                f"HF).\n"
            )
            samplesid_files = []
    # record the publish shape so a consumer (list_available_downloads / download)
    # can tell a dictless index + a genotype companion apart without unpacking.
    manifest["dictless"] = bool(dictless and vcf_name)
    manifest["has_genotypes"] = genotypes_dir is not None
    manifest["has_samplesids"] = bool(samplesid_files)
    tarball = f"{index_dir}.tar.gz"
    _make_index_tarball(
        tarball, index_dir, index_name, indels_dir, manifest, dict_dirs,
        samplesid_files=samplesid_files,
    )
    remote_path = f"indexes/{index_name}.tar.gz"
    api = hf.HfApi()
    api.upload_file(
        path_or_fileobj=tarball,
        path_in_repo=remote_path,
        repo_id=repo,
        repo_type="dataset",
        token=token,
    )
    os.remove(tarball)
    # SEPARATE upload of the big Tier-1 genotype store (only when it exists), to
    # the SAME repo + indexes/ prefix as the main tarball, named after the vcf so
    # download can find it deterministically from the index name.
    if genotypes_dir is not None:
        gt_tarball = os.path.join(os.path.dirname(index_dir),
                                  f"genotypes_{vcf_name}.tar.gz")
        _make_genotypes_tarball(gt_tarball, genotypes_dir)
        gt_remote = f"indexes/genotypes_{vcf_name}.tar.gz"
        api.upload_file(
            path_or_fileobj=gt_tarball,
            path_in_repo=gt_remote,
            repo_id=repo,
            repo_type="dataset",
            token=token,
        )
        os.remove(gt_tarball)
        sys.stdout.write(
            f"Published Tier-1 genotype store to {repo}:{gt_remote}\n"
        )
    return remote_path
