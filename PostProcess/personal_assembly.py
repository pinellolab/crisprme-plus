#!/usr/bin/env python3
"""Shared, dependency-light helpers for the *personal-assembly* on-disk layout.

A personal (diploid) assembly is stored as **one self-contained folder per
individual**, with a single ``metadata.json`` describing every file:

    Assemblies/
      <individual>/
        metadata.json
        paternal/
          genome/            <- one plain .fa per contig (CRISPRme Genomes layout)
            chr1.fa … chrN.fa
          <assembly>_vs_GRCh38.chain.gz
          <assembly>.chromAlias.txt
        maternal/
          genome/  …
          <assembly>_vs_GRCh38.chain.gz
          <assembly>.chromAlias.txt

Why a folder + metadata rather than the older flat ``Genomes/<name>_<hap>/`` +
``LiftoverFiles/<file>`` + hidden ``.assembly_individual`` markers:

  * **Delete one individual = remove one folder** (``rm -rf Assemblies/<id>/``),
    which actually reclaims the multi-GB genome files -- the old marker-only
    delete left them on disk.
  * **Clear separation** -- one individual never interleaves with another's
    files or with reference genomes/indexes.
  * **Portable / bring-your-own** -- a user can share or import a personal
    assembly as a ``.zip``/``.tar.gz`` of the same ``<individual>/`` structure
    plus a valid ``metadata.json``; nothing else is needed to register it.

All paths inside ``metadata.json`` are **relative to the individual's folder**,
so a bundle is relocatable. This module has **no third-party dependencies**
(stdlib only) so it can be imported from the CLI (``crisprme.py``), the HPRC
downloader (``download_hprc_assembly.py``), and the web app (``pages/*``) alike.

The web layer additionally understands the **legacy** flat/marker layout for
backward compatibility (see ``pages.pages_utils.installed_assemblies``); this
module deals only with the folder+metadata bundle format.
"""
from __future__ import annotations

import datetime
import json
import os
import shutil
import tarfile
import zipfile
from typing import Dict, List, Optional

import safe_tar  # (PostProcess/, the standard library only)

ASSEMBLIES_DIR = "Assemblies"       # top-level folder holding one dir per individual
METADATA_NAME = "metadata.json"     # the single source of truth inside each individual dir
BUNDLE_FORMAT = "crisprme-personal-assembly"
BUNDLE_FORMAT_VERSION = 1
HAPLOTYPES = ("paternal", "maternal")


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
def individual_dir(cwd: str, individual: str) -> str:
    """Absolute path of ``Assemblies/<individual>/`` (basename-guarded)."""
    return os.path.join(cwd, ASSEMBLIES_DIR, os.path.basename(individual))


def metadata_path(cwd: str, individual: str) -> str:
    return os.path.join(individual_dir(cwd, individual), METADATA_NAME)


# ---------------------------------------------------------------------------
# metadata.json read / write
# ---------------------------------------------------------------------------
def read_metadata(cwd: str, individual: str) -> Optional[Dict]:
    """Parse ``Assemblies/<individual>/metadata.json``; ``None`` if absent/invalid."""
    try:
        with open(metadata_path(cwd, individual)) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return None
    if not isinstance(meta, dict) or meta.get("format") != BUNDLE_FORMAT:
        return None
    return meta


def _now_iso() -> str:
    return datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")


def write_haplotype(
    cwd: str,
    individual: str,
    haplotype: str,
    entry: Dict,
    source: Optional[str] = None,
    reference: str = "GRCh38",
) -> None:
    """Record one haplotype in the individual's ``metadata.json``, merging into
    any existing file (so fetching paternal then maternal accumulates rather
    than clobbers). ``entry`` holds that haplotype's paths **relative to the
    individual folder** plus provenance, e.g.::

        {"assembly_name": "HG01255_paternal_CM086702.1",
         "genome_dir": "paternal/genome",
         "chain": "paternal/HG01255_paternal_CM086702.1_vs_GRCh38.chain.gz",
         "chromalias": "paternal/HG01255_paternal_CM086702.1.chromAlias.txt",
         "n_contigs": 92, "fasta_md5": "…"}
    """
    if haplotype not in HAPLOTYPES:
        raise ValueError(f"haplotype must be one of {HAPLOTYPES}, got {haplotype!r}")
    d = individual_dir(cwd, individual)
    os.makedirs(d, exist_ok=True)
    meta = read_metadata(cwd, individual) or {
        "format": BUNDLE_FORMAT,
        "format_version": BUNDLE_FORMAT_VERSION,
        "individual": os.path.basename(individual),
        "reference": reference,
        "created_utc": _now_iso(),
        "haplotypes": {},
    }
    if source and not meta.get("source"):
        meta["source"] = source
    meta.setdefault("haplotypes", {})[haplotype] = entry
    meta["updated_utc"] = _now_iso()
    tmp = metadata_path(cwd, individual) + ".part"
    with open(tmp, "w") as fh:
        json.dump(meta, fh, indent=2, sort_keys=False)
        fh.write("\n")
    os.replace(tmp, metadata_path(cwd, individual))


# ---------------------------------------------------------------------------
# Resolution (metadata -> concrete, existence-checked paths)
# ---------------------------------------------------------------------------
def _inside(bundle_dir: str, rel) -> bool:
    """Whether ``rel``, a path ``metadata.json`` gives, names a file or folder inside
    the bundle folder ``bundle_dir`` once '..' and links are resolved. metadata.json is
    the bundle's own content, a client's for an archive imported in Settings: an
    absolute path or '../' in it named any file or folder of the server, which the
    individual was then searched with."""
    if not isinstance(rel, str):
        return False
    root = os.path.realpath(bundle_dir)
    return os.path.commonpath([root, os.path.realpath(os.path.join(bundle_dir, rel))]) == root


def _resolve_haplotype(cwd: str, individual: str, hap_meta: Dict) -> Optional[Dict]:
    """Turn one haplotype's relative metadata into cwd-relative, existence-checked
    paths. Returns ``None`` if the metadata is missing required keys, or names a file
    or folder outside the individual's folder (a bundle imported before 2.8.1, whose
    import did not check that)."""
    root = individual_dir(cwd, individual)
    genome_rel = hap_meta.get("genome_dir")
    chain_rel = hap_meta.get("chain")
    chromalias_rel = hap_meta.get("chromalias")
    if not (genome_rel and chain_rel and chromalias_rel):
        return None
    if not all(_inside(root, rel) for rel in (genome_rel, chain_rel, chromalias_rel)):
        return None
    genome_abs = os.path.join(root, genome_rel)
    chain_abs = os.path.join(root, chain_rel)
    chromalias_abs = os.path.join(root, chromalias_rel)

    def _rel(p: str) -> str:
        return os.path.relpath(p, cwd)

    def _genome_ok(p: str) -> bool:
        return os.path.isdir(p) and any(f.endswith(".fa") for f in os.listdir(p)) if os.path.isdir(p) else False

    present = _genome_ok(genome_abs) and os.path.isfile(chain_abs) and os.path.isfile(chromalias_abs)
    assembly_name = hap_meta.get("assembly_name") or os.path.basename(genome_rel)
    return {
        # display strings (kept truthy so completeness/gap checks work)
        "genome": assembly_name,
        "chain": os.path.basename(chain_rel),
        "chromalias": os.path.basename(chromalias_rel),
        # cwd-relative resolved paths for the search engine
        "genome_path": _rel(genome_abs),
        "chain_path": _rel(chain_abs),
        "chromalias_path": _rel(chromalias_abs),
        "assembly_name": assembly_name,
        "present": present,
    }


def resolve(cwd: str, individual: str) -> Optional[Dict]:
    """One individual's registration dict (bundle layout), or ``None`` if there
    is no ``metadata.json`` for it. Shape matches
    ``pages.pages_utils.installed_assemblies`` entries, with extra ``layout`` /
    ``root`` / ``source`` / ``*_path`` fields. The individual is named after its
    folder, as ``--assembly-individual`` and Settings' Remove take it, not after
    metadata.json's ``individual``: that is the bundle's own content, a client's
    for an archive imported in Settings, and naming another individual listed the
    bundle in that one's place, or one with no folder, which Remove then took for
    a pairing of the legacy layout, removing paths named after the bundle's
    assembly names (outside the data folder too)."""
    meta = read_metadata(cwd, individual)
    if meta is None:
        return None
    haps = meta.get("haplotypes", {}) or {}
    out = {
        "individual": os.path.basename(individual),
        "layout": "bundle",
        "root": os.path.relpath(individual_dir(cwd, individual), cwd),
        "source": meta.get("source"),
        "paternal": None,
        "maternal": None,
    }
    for hap in HAPLOTYPES:
        if hap in haps:
            out[hap] = _resolve_haplotype(cwd, individual, haps[hap])

    def _complete(h: Optional[Dict]) -> bool:
        return bool(h) and bool(h.get("present"))

    out["complete"] = _complete(out["paternal"]) and _complete(out["maternal"])
    return out


def bundle_haplotype(genome_dir: str, haplotype: Optional[str] = None) -> Optional[Dict]:
    """The bundle haplotype whose genome folder is ``genome_dir``, found from the
    folder alone, so it works for the explicit ``--genome-paternal`` /
    ``--genome-maternal`` paths (the web passes those) as well as for
    ``--assembly-individual``.

    Looks for a bundle ``metadata.json`` in each folder above ``genome_dir`` (a
    bundle may live outside this working directory's ``Assemblies/``) and returns
    the haplotype whose ``genome_dir`` resolves to the same folder::

        {"root": <bundle folder>, "folder": <its name>, "individual": ...,
         "haplotype": "paternal"|"maternal", "assembly_name": ... or None}

    When both haplotypes name that folder, ``haplotype`` wins. ``None`` when the
    folder is no bundle haplotype's genome (a legacy ``Genomes/<name>`` folder, or
    any other folder)."""
    real = os.path.realpath(genome_dir)
    level = os.path.dirname(real)
    while True:
        try:
            with open(os.path.join(level, METADATA_NAME)) as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            meta = None
        if isinstance(meta, dict) and meta.get("format") == BUNDLE_FORMAT:
            haps = meta.get("haplotypes") or {}
            for hap in [haplotype] + [h for h in HAPLOTYPES if h != haplotype]:
                hap_meta = haps.get(hap) if hap else None
                rel = hap_meta.get("genome_dir") if isinstance(hap_meta, dict) else None
                if rel and os.path.realpath(os.path.join(level, rel)) == real:
                    return {
                        "root": level,
                        "folder": os.path.basename(level),
                        "individual": meta.get("individual") or os.path.basename(level),
                        "haplotype": hap,
                        "assembly_name": hap_meta.get("assembly_name") or None,
                    }
        parent = os.path.dirname(level)
        if parent == level:
            return None
        level = parent


def haplotype_label(genome_dir: str, haplotype: Optional[str] = None) -> str:
    """Display name of one haplotype's genome, for the ``Genome_paternal`` /
    ``Genome_maternal`` lines of an assembly-search job's ``.Params.txt`` (the
    results header and History show them): the bundle's ``assembly_name``
    (``<individual>_<haplotype>`` when it records none), else the folder's own
    name, as in the legacy ``Genomes/<name>_<hap>`` layout. A bundle keeps every
    haplotype in a folder named ``genome``, so there the folder name says nothing."""
    hap = bundle_haplotype(genome_dir, haplotype)
    if hap is not None:
        label = hap["assembly_name"] or f"{hap['individual']}_{hap['haplotype']}"
    else:
        label = os.path.basename(os.path.normpath(genome_dir))
    # one field of a tab-separated file
    return " ".join(str(label).split())


def discover_bundles(cwd: str) -> List[Dict]:
    """All individuals under ``Assemblies/`` that carry a ``metadata.json``."""
    root = os.path.join(cwd, ASSEMBLIES_DIR)
    if not os.path.isdir(root):
        return []
    out = []
    for name in sorted(os.listdir(root)):
        if not os.path.isdir(os.path.join(root, name)):
            continue
        entry = resolve(cwd, name)
        if entry is not None:
            out.append(entry)
    return out


def delete_target(cwd: str, individual: str) -> Optional[str]:
    """Absolute path of the whole individual folder to remove, or ``None`` if
    this individual isn't a bundle (caller should fall back to legacy un-pair)."""
    d = individual_dir(cwd, individual)
    if os.path.isfile(os.path.join(d, METADATA_NAME)):
        return d
    return None


# ---------------------------------------------------------------------------
# Import a bring-your-own bundle archive (.zip / .tar.gz)
# ---------------------------------------------------------------------------
def _unusable_fasta_name(genome_dir: str) -> Optional[str]:
    """Why a FASTA file of the haplotype genome folder ``genome_dir`` cannot be
    searched under its name, else ``None``. A search takes each ``.fa`` file's name,
    less the extension, for its chromosome's, and its post-analysis takes that name
    into shell commands as it is (pool_post_analisi_snp.py refuses one it cannot
    take). A bundle brings its files named as its archive names them, a client's for
    one imported in Settings: a name holding a character a shell reads ran it when the
    individual was searched. Refused as in a sequence's name (genome_layout.py), and
    so is any blank; a real chromosome's file, named after its sequence, holds none."""
    from genome_layout import UNUSABLE_IN_ID

    for fname in sorted(os.listdir(genome_dir)):
        if not fname.endswith(".fa"):
            continue
        chrom = os.path.splitext(fname)[0]
        if any(c in UNUSABLE_IN_ID or c.isspace() for c in chrom):
            return (
                f"holds {fname!r}: a FASTA file's name must not hold a blank or any of "
                f"{' '.join(UNUSABLE_IN_ID)}, which a search takes into shell commands as "
                f"the chromosome's; name each file after the sequence it holds"
            )
    return None


def validate_bundle_dir(bundle_dir: str) -> Optional[str]:
    """Return an error string if ``bundle_dir`` is not a valid personal-assembly
    bundle (has a parseable ``metadata.json`` of a supported version, every file
    it references exists, inside the bundle, and its FASTA files are named as a
    search can take them), else ``None``."""
    mpath = os.path.join(bundle_dir, METADATA_NAME)
    try:
        with open(mpath) as fh:
            meta = json.load(fh)
    except OSError:
        return f"missing {METADATA_NAME} in the bundle"
    except ValueError:
        return f"{METADATA_NAME} is not valid JSON"
    if meta.get("format") != BUNDLE_FORMAT:
        return f"{METADATA_NAME} 'format' must be {BUNDLE_FORMAT!r}"
    if meta.get("format_version", 1) > BUNDLE_FORMAT_VERSION:
        return (
            f"bundle format_version {meta.get('format_version')} is newer than this "
            f"CRISPRme supports ({BUNDLE_FORMAT_VERSION}); please update CRISPRme"
        )
    haps = meta.get("haplotypes") or {}
    if not haps:
        return f"{METADATA_NAME} lists no haplotypes"
    for hap, hm in haps.items():
        if hap not in HAPLOTYPES:
            return f"unexpected haplotype {hap!r} (expected paternal/maternal)"
        genome_rel, chain_rel, ca_rel = hm.get("genome_dir"), hm.get("chain"), hm.get("chromalias")
        if not (genome_rel and chain_rel and ca_rel):
            return f"{hap} entry is missing genome_dir/chain/chromalias"
        for key, rel in (("genome_dir", genome_rel), ("chain", chain_rel), ("chromalias", ca_rel)):
            if not _inside(bundle_dir, rel):
                return f"{hap} {key} {rel!r} is not a path inside the bundle"
        gdir = os.path.join(bundle_dir, genome_rel)
        if not (os.path.isdir(gdir) and any(f.endswith(".fa") for f in os.listdir(gdir))):
            return f"{hap} genome dir {genome_rel!r} missing or has no .fa contigs"
        unusable = _unusable_fasta_name(gdir)
        if unusable:
            return f"{hap} genome dir {genome_rel!r} {unusable}"
        if not os.path.isfile(os.path.join(bundle_dir, chain_rel)):
            return f"{hap} chain file {chain_rel!r} not found in the bundle"
        if not os.path.isfile(os.path.join(bundle_dir, ca_rel)):
            return f"{hap} chromAlias file {ca_rel!r} not found in the bundle"
    return None


def _safe_members_top_dir(names: List[str]) -> Optional[str]:
    """The single top-level directory all archive entries live under, or None.

    Rejects, BEFORE any normalization, absolute paths (POSIX ``/x`` or Windows
    ``C:\\x``) and every form of ``..`` traversal -- so a member like
    ``/etc/passwd`` (which would otherwise survive an ``lstrip('/')``) or
    ``a/../..`` is refused rather than silently rewritten. Every member is
    checked, not just the first, so a pass here guarantees all entries are
    relative and confined under one top folder."""
    tops = set()
    for raw in names:
        n = raw.replace("\\", "/")
        # absolute paths: POSIX (leading /) or Windows drive (X:...)
        if n.startswith("/") or (len(n) >= 2 and n[1] == ":"):
            return None
        if (
            not n
            or n == ".."
            or n.startswith("../")
            or "/../" in n
            or n.endswith("/..")
        ):
            return None
        top = n.split("/", 1)[0]
        if not top or top == ".":
            return None
        tops.add(top)
    if len(tops) != 1:
        return None
    return tops.pop()


def import_archive(archive_path: str, cwd: str) -> str:
    """Extract a personal-assembly bundle archive into ``Assemblies/`` and
    validate it. The archive must contain exactly one top-level ``<individual>/``
    directory holding the bundle structure + ``metadata.json``. Returns the
    individual name on success; raises ``ValueError`` with a clear reason
    otherwise (and leaves nothing half-extracted)."""
    assemblies_root = os.path.join(cwd, ASSEMBLIES_DIR)
    os.makedirs(assemblies_root, exist_ok=True)

    if archive_path.endswith((".tar.gz", ".tgz")):
        opener = lambda: tarfile.open(archive_path, "r:gz")  # noqa: E731
        namer = lambda a: a.getnames()  # noqa: E731
        is_tar = True
    elif archive_path.endswith(".zip"):
        opener = lambda: zipfile.ZipFile(archive_path)  # noqa: E731
        namer = lambda a: a.namelist()  # noqa: E731
        is_tar = False
    else:
        raise ValueError("bundle must be a .zip or .tar.gz/.tgz archive")

    with opener() as arch:
        names = namer(arch)
        top = _safe_members_top_dir(names)
        if top is None:
            raise ValueError(
                "archive must contain exactly one top-level <individual>/ folder "
                "and no absolute or '..' paths"
            )
        dest = os.path.join(assemblies_root, os.path.basename(top))
        if os.path.exists(dest):
            raise ValueError(
                f"an assembly named {os.path.basename(top)!r} already exists -- "
                "delete it first or rename the bundle's top folder"
            )
        staging = dest + ".importing"
        if os.path.exists(staging):
            shutil.rmtree(staging, ignore_errors=True)
        os.makedirs(staging)
        try:
            # Defense-in-depth: _safe_members_top_dir already rejected absolute
            # and '..' members, but also (1) refuse tar symlink/hardlink/device
            # entries -- a symlink member could otherwise redirect a later
            # write outside staging -- and (2) verify every member still
            # resolves inside staging before extracting anything.
            if is_tar:
                for m in arch.getmembers():
                    if not (m.isreg() or m.isdir()):
                        raise ValueError(
                            f"unsupported archive entry {m.name!r} "
                            "(only regular files and directories are allowed)"
                        )
            staging_real = os.path.realpath(staging)
            for n in names:
                target = os.path.realpath(os.path.join(staging, n))
                if target != staging_real and not target.startswith(staging_real + os.sep):
                    raise ValueError(f"archive member escapes the extraction directory: {n!r}")
            # extract into staging/<top>/… then validate before publishing: a tar
            # archive's files with safe modes, owned by this user, not with the modes
            # and owners it gives them (safe_tar); zipfile sets neither
            if is_tar:
                safe_tar.extract(arch, staging)
            else:
                arch.extractall(staging)
            extracted = os.path.join(staging, top)
            err = validate_bundle_dir(extracted)
            if err:
                raise ValueError(f"invalid personal-assembly bundle: {err}")
            os.replace(extracted, dest)
        finally:
            shutil.rmtree(staging, ignore_errors=True)
    return os.path.basename(top)
