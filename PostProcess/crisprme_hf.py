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
from datetime import datetime, timezone

# Default HuggingFace dataset repo. Overridable per-invocation with --hf-repo or
# the CRISPRME_HF_REPO environment variable (e.g. switch to a pinellolab org
# repo once it exists).
DEFAULT_HF_REPO = "lucapinello/crisprme-data"

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
}

# Local destination directory for each component (relative to the working dir).
_COMPONENT_LOCALDIR = {
    "genome": "Genomes",
    "annotations": "Annotations",
    "pams": "PAMs",
    "samples": "samplesIDs",
    "vcf": "VCFs",
    "index": "genome_library",
}


def resolve_repo(cli_repo: Optional[str] = None) -> str:
    """Resolve the HuggingFace repo id: CLI value > env var > built-in default."""
    if cli_repo:
        return cli_repo
    return os.environ.get("CRISPRME_HF_REPO", DEFAULT_HF_REPO)


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
    """
    if component not in ("genome", "vcf", "index"):
        return []
    try:
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
    except Exception:
        return []


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


def _hf_snapshot(repo: str, allow_patterns: List[str], local_dir: str,
                 token: Optional[str] = None) -> str:
    """Download the files matching ``allow_patterns`` from ``repo`` into
    ``local_dir`` (flat, following the repo tree). Returns ``local_dir``."""
    hf = _require_hf()
    os.makedirs(local_dir, exist_ok=True)
    hf.snapshot_download(
        repo_id=repo,
        repo_type="dataset",
        allow_patterns=allow_patterns,
        local_dir=local_dir,
        token=token,
    )
    return local_dir


def _fetch_genotypes_companion(
    repo: str,
    workdir: str,
    vcf_name: str,
    remote_prefix: str,
    token: Optional[str] = None,
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
    _hf_snapshot(repo, patterns, staging, token)
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
        with tarfile.open(gt_tarball) as tf:
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
    repo = resolve_repo(repo)
    remote_prefix = _COMPONENT_PREFIXES[component]
    local_dir = os.path.join(workdir, _COMPONENT_LOCALDIR[component])

    if component == "genome":
        patterns = [f"{remote_prefix}/{ref}/*"]
        staging = os.path.join(workdir, ".hf_stage_genome")
        _hf_snapshot(repo, patterns, staging, token)
        src = os.path.join(staging, remote_prefix, ref)
        if not os.path.isdir(src):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No genome '{ref}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{ref}/'). Check --ref, or the repo "
                f"may not have this genome uploaded yet."
            )
        dest = os.path.join(local_dir, ref)
        os.makedirs(dest, exist_ok=True)
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
        for fn in files:
            moved = shutil.move(os.path.join(src, fn), os.path.join(dest, fn))
            if moved.endswith(".gz"):
                decompress_gz(moved)
        shutil.rmtree(staging, ignore_errors=True)
        return dest

    if component == "vcf":
        if not dataset:
            raise ValueError("--dataset is required to download a 'vcf' component")
        patterns = [f"{remote_prefix}/{dataset}/*"]
        staging = os.path.join(workdir, ".hf_stage_vcf")
        _hf_snapshot(repo, patterns, staging, token)
        src = os.path.join(staging, remote_prefix, dataset)
        if not os.path.isdir(src):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No VCF dataset '{dataset}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{dataset}/'). Check --dataset, or "
                f"the repo may not have this dataset uploaded yet."
            )
        dest = os.path.join(local_dir, dataset)
        os.makedirs(dest, exist_ok=True)
        for fn in sorted(os.listdir(src)):
            shutil.move(os.path.join(src, fn), os.path.join(dest, fn))  # keep bgzip
        shutil.rmtree(staging, ignore_errors=True)
        return dest

    if component == "index":
        if not index_name:
            raise ValueError("--index-name is required to download an 'index' component")
        patterns = [f"{remote_prefix}/{index_name}.tar.gz"]
        staging = os.path.join(workdir, ".hf_stage_index")
        _hf_snapshot(repo, patterns, staging, token)
        tarball = os.path.join(staging, remote_prefix, f"{index_name}.tar.gz")
        if not os.path.isfile(tarball):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No index '{index_name}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{index_name}.tar.gz'). Check "
                f"--index-name, or build it locally with 'crisprme.py "
                f"build-index-only' — the repo may not have it uploaded yet."
            )
        os.makedirs(local_dir, exist_ok=True)
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
            # atomically swap the index and its _INDELS companion into place
            for sub in (index_name, index_name + "_INDELS"):
                src = os.path.join(extract_tmp, sub)
                if not os.path.isdir(src):
                    continue  # reference-only indexes have no _INDELS companion
                dst = os.path.join(local_dir, sub)
                backup = os.path.join(local_dir, f".{sub}.replaced")  # hidden -> unlisted
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
                try:
                    got = _fetch_genotypes_companion(
                        repo, workdir, vcf_name, remote_prefix, token
                    )
                except Exception as exc:  # network/extract error -> non-fatal
                    got = False
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
                else:
                    sys.stdout.write(
                        f"NOTE: index '{index_name}' has no genotypes_{vcf_name}.tar.gz "
                        f"companion in {repo} (registry-only/detection publish); "
                        f"off-target detection works, per-sample Samples degraded.\n"
                    )
        return os.path.join(local_dir, index_name)

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
) -> None:
    """Builds the one-file-per-index ``.tar.gz`` (index + _INDELS + dicts + manifest).

    Prefers ``pigz`` (parallel gzip) via GNU ``tar`` so a large index compresses
    across all cores instead of one — a 170GB pamless index tars in minutes with
    pigz vs ~a day with single-threaded gzip. Falls back to Python's (single-
    threaded) ``tarfile`` when ``pigz``/``tar`` are unavailable. Both paths write
    an identical archive layout: ``<index_name>/``, optional ``<index_name>_INDELS/``,
    each ``dict_dirs`` entry stored under ``Dictionaries/<basename>/`` (the per-sample
    variant dictionaries a variant search needs at post-analysis time — bundling them
    lets a downloaded index search WITHOUT the source VCFs), and ``manifest.json`` at
    the archive root (see download_component's extractor).
    """
    manifest_bytes = json.dumps(manifest, indent=2).encode()
    pigz = shutil.which("pigz")
    tar = shutil.which("tar")
    if pigz and tar:
        # stage manifest.json in a temp dir so tar can add it at the archive root
        # with the same "manifest.json" name the extractor expects. All three -C
        # paths are absolute, so the multi -C chaining is order-independent.
        tmpd = tempfile.mkdtemp(prefix="hfpub_")
        try:
            with open(os.path.join(tmpd, "manifest.json"), "wb") as mf:
                mf.write(manifest_bytes)
            parent = os.path.dirname(index_dir)
            # -C <dir> <name> stores <name> under its basename == the desired arcname
            inputs = ["-C", parent, index_name]
            if os.path.isdir(indels_dir):
                inputs += ["-C", parent, os.path.basename(indels_dir)]
            for dp in (dict_dirs or []):
                # store as Dictionaries/<basename>/... so the extractor routes it to
                # <workdir>/Dictionaries/ (not genome_library/)
                inputs += ["-C", os.path.dirname(os.path.dirname(dp)),
                           os.path.join("Dictionaries", os.path.basename(dp))]
            inputs += ["-C", tmpd, "manifest.json"]
            subprocess.check_call(
                [tar, "--use-compress-program", pigz, "-cf", tarball, *inputs]
            )
            return
        except (subprocess.CalledProcessError, OSError) as exc:
            # this tar may not support --use-compress-program (busybox / very old),
            # or pigz died mid-stream -> drop any partial tarball and fall through
            # to the single-threaded Python path so publish still succeeds.
            sys.stderr.write(
                f"Note: parallel (pigz) compression failed ({exc}); falling back "
                f"to single-threaded gzip.\n"
            )
            if os.path.exists(tarball):
                os.remove(tarball)
        finally:
            shutil.rmtree(tmpd, ignore_errors=True)
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
        info = tarfile.TarInfo(name="manifest.json")
        info.size = len(manifest_bytes)
        tf.addfile(info, io.BytesIO(manifest_bytes))


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

    Returns:
        The remote path (``indexes/<index_name>.tar.gz``) the main index was
        uploaded to.
    """
    index_dir = os.path.abspath(index_dir.rstrip("/"))
    if not os.path.isdir(index_dir):
        raise ValueError(f"Index directory {index_dir} does not exist")
    index_name = os.path.basename(index_dir)
    repo = resolve_repo(repo)
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
        # ADDITIVE Tier-1: the big genotype store travels as a SEPARATE tarball
        # (see below), never inside the main archive.
        gt_p = os.path.join(dicts_root, genotypes_name)
        if os.path.isdir(gt_p) and os.listdir(gt_p):
            genotypes_dir = gt_p
    # record the publish shape so a consumer (list_available_downloads / download)
    # can tell a dictless index + a genotype companion apart without unpacking.
    manifest["dictless"] = bool(dictless and vcf_name)
    manifest["has_genotypes"] = genotypes_dir is not None
    tarball = f"{index_dir}.tar.gz"
    _make_index_tarball(tarball, index_dir, index_name, indels_dir, manifest, dict_dirs)
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
