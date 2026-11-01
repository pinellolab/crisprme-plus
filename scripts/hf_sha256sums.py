#!/usr/bin/env python3
"""hf_sha256sums.py - Write the ``SHA256SUMS`` for a folder of a Hugging Face dataset repo.

A file stored through Git LFS is addressed by its SHA-256 (the LFS ``oid``), so the
checksum of every archive under ``indexes/`` is already recorded by the Hub and needs
no download. This prints them in ``sha256sum`` format, one line per ``.tar.gz``
directly under the folder, sorted by name, so that from inside a downloaded
``indexes/``::

    sha256sum -c SHA256SUMS

verifies each archive. (Checked on ``build-inputs/``, whose ``.SHA256SUMS`` were
computed from the tarballs themselves: the oids match them.)

``--staged`` takes a ``sha256sum``-format file for archives about to be uploaded
(computed where they were built) and lets each of its lines replace or add an
entry. That gives the SHA256SUMS to upload WITH them; re-running without it once
the upload lands must print the same bytes, which is what proves the Hub holds what
was staged.

Exits non-zero when no archive is listed (an empty SHA256SUMS verifies nothing and
would read as a pass) or when an archive has no LFS checksum. Standard library only.

Usage:
    python scripts/hf_sha256sums.py lucapinello/crisprme-data > SHA256SUMS
    python scripts/hf_sha256sums.py lucapinello/crisprme-data --staged staged.sha256
    python scripts/hf_sha256sums.py lucapinello/crisprme-data-general-index --token "$HF_TOKEN"
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "https://huggingface.co/api/datasets/{repo}/tree/{revision}/{path}?expand=true"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _get_json(url: str, token: str | None):
    """GETs ``url``; returns (decoded JSON, URL of the next page or None)."""
    req = urllib.request.Request(url)
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(req, timeout=120) as resp:
        nxt = None
        # the tree endpoint paginates through an RFC 8288 Link header
        for part in (resp.headers.get("Link") or "").split(","):
            m = re.match(r'\s*<([^>]+)>\s*;\s*rel="next"', part)
            if m:
                nxt = m.group(1)
        return json.load(resp), nxt


def hub_checksums(repo: str, path: str, revision: str = "main", token: str | None = None) -> dict:
    """Maps each ``.tar.gz`` directly under ``path`` to its LFS SHA-256."""
    url = API.format(
        repo=repo,
        revision=urllib.parse.quote(revision, safe=""),
        path=urllib.parse.quote(path.strip("/")),
    )
    sums, missing = {}, []
    while url:
        entries, url = _get_json(url, token)
        for e in entries:
            if e.get("type") != "file" or not e["path"].endswith(".tar.gz"):
                continue
            name = e["path"].rsplit("/", 1)[-1]
            oid = (e.get("lfs") or {}).get("oid", "")
            if SHA256_RE.match(oid):
                sums[name] = oid
            else:
                missing.append(name)
    if missing:
        raise SystemExit(f"no LFS SHA-256 recorded for: {', '.join(sorted(missing))}")
    return sums


def read_staged(path: str) -> dict:
    """Parses a ``sha256sum``-format file (``<hex>  [*]<name>``) into {basename: hex}."""
    sums = {}
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip():
                continue
            parts = line.split(None, 1)
            if len(parts) != 2 or not SHA256_RE.match(parts[0]):
                raise SystemExit(f"{path}:{n}: not a sha256sum line: {line.rstrip()}")
            sums[os.path.basename(parts[1].strip().lstrip("*"))] = parts[0]
    return sums


def render(sums: dict) -> str:
    """``sha256sum`` format: hex, two spaces, name; sorted by name."""
    return "".join(f"{sums[name]}  {name}\n" for name in sorted(sums))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("repo", help="dataset repo id, e.g. lucapinello/crisprme-data")
    ap.add_argument("--path", default="indexes", help="folder in the repo (default: indexes)")
    ap.add_argument("--revision", default="main")
    ap.add_argument("--token", default=os.environ.get("HF_TOKEN"),
                    help="only needed for a private repo (default: $HF_TOKEN)")
    ap.add_argument("--staged", help="sha256sum file of archives about to be uploaded")
    args = ap.parse_args(argv)

    try:
        sums = hub_checksums(args.repo, args.path, args.revision, args.token)
    except urllib.error.HTTPError as exc:
        # the Hub answers 401 both for a private repo without a token and for one
        # that does not exist, so say both
        raise SystemExit(
            f"HTTP {exc.code} listing {args.repo}/{args.path}: the repo is private "
            f"(pass --token) or does not exist"
        )
    if args.staged:
        for name, hexsum in read_staged(args.staged).items():
            was = sums.get(name)
            verb = "adds" if was is None else ("keeps" if was == hexsum else "replaces")
            sys.stderr.write(f"staged {verb} {name}\n")
            sums[name] = hexsum
    if not sums:
        raise SystemExit(f"no .tar.gz under {args.repo}/{args.path}: refusing to write an empty SHA256SUMS")
    sys.stdout.write(render(sums))
    return 0


if __name__ == "__main__":
    sys.exit(main())
