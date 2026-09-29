#!/usr/bin/env python3
"""Scrape the All of Us (AoU) public Data Browser variant API for one chromosome.

The All of Us Research Program exposes an aggregate variant browser backed by a
public JSON API (``public.api.researchallofus.org``) that needs no
authentication. This script walks a single chromosome in fixed-size windows,
paginates the API within each window, and saves every response page as a raw
JSON file. It is **resumable**: progress is checkpointed to ``progress.json`` so
an interrupted run continues where it stopped, and a chromosome that has already
been fully scraped (marked by a ``.complete`` file) is skipped.

The API is aggregate and sites-only — each variant carries allele
counts/frequencies, not individual genotypes.

Usage
-----
    scrape_aou.py --chrom chr22 --length 50818468 --outdir json/chr22
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from typing import Dict, List

import requests

API_URL = "https://public.api.researchallofus.org/v1/genomics/search-variants"

# Browser-like headers; the endpoint is public but expects a JSON content type
# and the Data Browser origin.
HEADERS = {
    "User-Agent": "Mozilla/5.0",
    "Accept": "*/*",
    "Content-Type": "application/json",
    "Origin": "https://databrowser.researchallofus.org",
    "Referer": "https://databrowser.researchallofus.org/",
}

# The Data Browser sends a sortMetadata block; sorting by variantId ascending
# gives a stable pagination order across pages.
SORT_METADATA = {
    "variantId": {"sortActive": True, "sortDirection": "asc", "sortOrder": 1},
}


def normalize_chrom(chrom: str) -> str:
    """Return *chrom* in ``chr``-prefixed form.

    Parameters
    ----------
    chrom : str
        Chromosome label such as ``"1"``, ``"chr1"``, ``"X"`` or ``"chrX"``.

    Returns
    -------
    str
        The ``chr``-prefixed label (e.g. ``"chr1"``).
    """
    chrom = chrom.strip()
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def post_with_retries(
    payload: dict, max_retries: int = 5, timeout: int = 60
) -> requests.Response:
    """POST *payload* to the AoU API, retrying transient failures with backoff.

    Parameters
    ----------
    payload : dict
        JSON request body.
    max_retries : int, optional
        Maximum number of attempts before giving up (default 5).
    timeout : int, optional
        Per-request timeout in seconds (default 60).

    Returns
    -------
    requests.Response
        The successful (HTTP 200) response.

    Raises
    ------
    SystemExit
        If every attempt fails.
    """
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.post(API_URL, headers=HEADERS, json=payload, timeout=timeout)
            if resp.status_code == 200:
                return resp
            sys.stderr.write(
                f"HTTP {resp.status_code} (attempt {attempt}/{max_retries})\n"
            )
        except requests.exceptions.RequestException as exc:
            sys.stderr.write(f"request error: {exc} (attempt {attempt}/{max_retries})\n")
        time.sleep(random.uniform(1.0, 5.0))
    sys.stderr.write("Max retries exceeded; aborting.\n")
    sys.exit(1)


def load_progress(path: str) -> Dict[str, int]:
    """Load a checkpoint, or a fresh starting point if none exists.

    Parameters
    ----------
    path : str
        Path to ``progress.json``.

    Returns
    -------
    dict
        Mapping with integer keys ``window_index``, ``page_number`` and
        ``file_counter`` (all 0/1 defaults when the file is absent).
    """
    if os.path.isfile(path):
        with open(path) as handle:
            data = json.load(handle)
        return {
            "window_index": int(data.get("window_index", 0)),
            "page_number": int(data.get("page_number", 1)),
            "file_counter": int(data.get("file_counter", 0)),
        }
    return {"window_index": 0, "page_number": 1, "file_counter": 0}


def save_progress(path: str, window_index: int, page_number: int, file_counter: int) -> None:
    """Write the current checkpoint to *path* (``progress.json``).

    Parameters
    ----------
    path : str
        Destination path.
    window_index : int
        Index of the window to resume at.
    page_number : int
        Page number to resume at within that window.
    file_counter : int
        Number of JSON pages saved so far.
    """
    with open(path, "w") as handle:
        json.dump(
            {
                "window_index": window_index,
                "page_number": page_number,
                "file_counter": file_counter,
            },
            handle,
        )


def scrape_chromosome(
    chrom: str,
    length: int,
    outdir: str,
    window: int = 1_000_000,
    row_count: int = 50_000,
    max_retries: int = 5,
) -> int:
    """Scrape one chromosome window-by-window into raw JSON pages.

    Parameters
    ----------
    chrom : str
        Chromosome label (any accepted by :func:`normalize_chrom`).
    length : int
        Chromosome length in base pairs (windows are generated up to this).
    outdir : str
        Directory for the ``file_<N>.json`` pages, ``progress.json`` and the
        ``.complete`` marker. Created if absent.
    window : int, optional
        Window size in base pairs (default 1,000,000).
    row_count : int, optional
        Rows per API page (default 50,000).
    max_retries : int, optional
        Per-request retry budget (default 5).

    Returns
    -------
    int
        The total number of JSON pages saved for this chromosome.
    """
    chrom = normalize_chrom(chrom)
    os.makedirs(outdir, exist_ok=True)
    complete_marker = os.path.join(outdir, ".complete")
    if os.path.isfile(complete_marker):
        sys.stderr.write(f"{chrom}: already complete ({outdir}), skipping.\n")
        return load_progress(os.path.join(outdir, "progress.json"))["file_counter"]

    progress_path = os.path.join(outdir, "progress.json")
    progress = load_progress(progress_path)
    file_counter = progress["file_counter"]
    total_windows = math.ceil(length / window)

    for window_index in range(progress["window_index"], total_windows):
        base = window_index * window
        start = base + 1 if base != 0 else 0
        end = min(base + window, length)
        query = f"{chrom}:{start}-{end}"
        # Resume mid-window only for the first window of this run.
        page_number = (
            progress["page_number"] if window_index == progress["window_index"] else 1
        )

        while True:
            payload = {
                "query": query,
                "pageNumber": page_number,
                "rowCount": row_count,
                "sortMetadata": SORT_METADATA,
                "filterMetadata": None,
            }
            resp = post_with_retries(payload, max_retries=max_retries)
            try:
                items = resp.json().get("items") or []
            except ValueError:
                items = []
            if not items:  # no more pages in this window
                break

            file_counter += 1
            page_path = os.path.join(outdir, f"file_{file_counter}.json")
            with open(page_path, "w") as handle:
                handle.write(resp.text)
            sys.stderr.write(f"{query} p{page_number}: {len(items)} -> {page_path}\n")

            page_number += 1
            save_progress(progress_path, window_index, page_number, file_counter)
            time.sleep(random.uniform(1.0, 5.0))

        # Window finished; checkpoint the start of the next window.
        save_progress(progress_path, window_index + 1, 1, file_counter)

    # Mark the whole chromosome as complete.
    with open(complete_marker, "w") as handle:
        handle.write("")
    sys.stderr.write(f"{chrom}: complete, {file_counter} pages in {outdir}\n")
    return file_counter


def main() -> None:
    """Parse arguments and scrape a single chromosome."""
    parser = argparse.ArgumentParser(
        description="Scrape the All of Us public variant API for one chromosome."
    )
    parser.add_argument("--chrom", required=True, help="Chromosome, e.g. chr22 or 22")
    parser.add_argument("--length", required=True, type=int, help="Chromosome length (bp)")
    parser.add_argument("--outdir", required=True, help="Output directory for JSON pages")
    parser.add_argument("--window", type=int, default=1_000_000, help="Window size (bp)")
    parser.add_argument("--rowcount", type=int, default=50_000, help="Rows per API page")
    parser.add_argument("--max-retries", type=int, default=5, help="Per-request retries")
    args = parser.parse_args()

    scrape_chromosome(
        chrom=args.chrom,
        length=args.length,
        outdir=args.outdir,
        window=args.window,
        row_count=args.rowcount,
        max_retries=args.max_retries,
    )


if __name__ == "__main__":
    main()
