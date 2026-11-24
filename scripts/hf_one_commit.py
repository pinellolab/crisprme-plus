#!/usr/bin/env python3
"""hf_one_commit.py - put every file of a flat folder into a Hugging Face dataset repo in
ONE commit.

Why not ``hf upload <folder>``: since huggingface_hub 1.19 ``upload_folder`` (which the CLI
calls) goes through a streamed pipeline when hf_xet is installed, and that pipeline closes
the running batch, and commits it, as soon as the batch is older than MAX_COMMIT_INTERVAL
(300 s). For six files whose first one is 29.55 GB, the archive's registration alone can
pass 300 s, and the upload then lands as "<title>" (the archive) plus "<title> (part 2)"
(the rest) - two commits, with a window in between in which the repo holds the new
companions and the old index archive, or the other way round. ``HfApi.create_commit`` has
no batching, no clock and no file count: it pre-uploads every file, then sends ONE commit
request, which is not retried (so a commit whose response was lost cannot be duplicated).

The write token is the one the machine is logged in with (``hf auth login``); HF_TOKEN in
the environment wins where it is set. No token is ever taken from the command line, and
none is printed.

``--parent`` is the branch head the SHA256SUMS was written from: the Hub refuses the commit
if anything else landed in between, so the folder's SHA256SUMS can never describe a repo
that has moved.

``--dry-run`` makes no Hub call at all beyond the login check: it prints the operations and
exits.

Usage:
  python3 hf_one_commit.py REPO FOLDER PATH_IN_REPO --message TITLE --parent OID [--dry-run]
"""
from __future__ import annotations

import argparse
import os
import re
import sys

OID = re.compile(r"^[0-9a-f]{40}$")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("repo", help="dataset repo id, e.g. lucapinello/crisprme_plus_indices")
    ap.add_argument("folder", help="a flat folder: every file in it goes up, nothing else")
    ap.add_argument("path_in_repo", help="the folder in the repo, e.g. indexes")
    ap.add_argument("--message", required=True, help="the commit title")
    ap.add_argument("--parent", required=True,
                    help="the 40-hex oid of the branch head this was prepared against")
    ap.add_argument("--revision", default="main")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the operations and stop; no commit, no upload")
    a = ap.parse_args(argv)

    if not OID.match(a.parent):
        sys.exit(f"--parent {a.parent!r} is not a 40-character commit oid")
    names = sorted(os.listdir(a.folder))
    odd = [n for n in names
           if n.startswith(".") or not os.path.isfile(os.path.join(a.folder, n))]
    if not names or odd:
        sys.exit(f"{a.folder}: not a flat folder of files ({odd or 'it is empty'}); "
                 f"nothing was uploaded")

    from huggingface_hub import CommitOperationAdd, HfApi
    from huggingface_hub.utils import get_token

    token = os.environ.get("HF_TOKEN") or get_token()
    if not token:
        sys.exit("This machine is not logged in to the Hub and HF_TOKEN is not set "
                 "(hf auth login). Nothing was uploaded.")

    ops, total = [], 0
    for n in names:
        local = os.path.join(a.folder, n)
        size = os.path.getsize(local)
        total += size
        ops.append(CommitOperationAdd(path_in_repo=f"{a.path_in_repo.strip('/')}/{n}",
                                      path_or_fileobj=local))
    sys.stdout.write(f"{len(ops)} file(s), {total:,} bytes, into {a.repo}:"
                     f"{a.path_in_repo.strip('/')}/ in one commit on {a.revision}\n")
    for op in ops:
        sys.stdout.write(f"  {op.path_in_repo}  "
                         f"{os.path.getsize(str(op.path_or_fileobj)):,} bytes\n")
    sys.stdout.write(f"parent commit required: {a.parent}\n")
    sys.stdout.write(f"title: {a.message}\n")
    sys.stdout.flush()
    if a.dry_run:
        sys.stdout.write("--dry-run: nothing was uploaded and no commit was made.\n")
        return 0

    info = HfApi(token=token).create_commit(
        repo_id=a.repo, repo_type="dataset", operations=ops,
        commit_message=a.message, revision=a.revision, parent_commit=a.parent)
    sys.stdout.write(f"ONE COMMIT: {info.oid}\n{info.commit_url}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
