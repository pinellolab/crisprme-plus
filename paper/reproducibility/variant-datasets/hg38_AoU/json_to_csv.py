#!/usr/bin/env python3
"""Combine scraped All of Us JSON pages for one chromosome into a single CSV.

Reads every ``file_*.json`` produced by ``scrape_aou.py`` in a directory,
extracts the ``items`` array from each page, concatenates them, drops duplicate
variants (by ``variantId``, which can recur at window boundaries), and writes one
CSV. The CSV keeps the API's field names so downstream tools can select columns
by name rather than position.

Usage
-----
    json_to_csv.py --json-dir json/chr22 --out-csv csv/chr22.csv
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from typing import Dict, List

import pandas as pd

# Field order the AoU API returns for each variant item.
ITEM_COLUMNS: List[str] = [
    "variantId",
    "genes",
    "consequence",
    "variantType",
    "proteinChange",
    "clinicalSignificance",
    "alleleCount",
    "alleleNumber",
    "alleleFrequency",
    "homozygoteCount",
]


def collect_items(json_dir: str) -> List[Dict]:
    """Read all ``file_*.json`` pages in *json_dir* and return their items.

    Empty or unreadable files are skipped with a warning.

    Parameters
    ----------
    json_dir : str
        Directory containing the scraped ``file_<N>.json`` pages.

    Returns
    -------
    list of dict
        The concatenated ``items`` records across every page (order preserved).
    """
    rows: List[Dict] = []
    files = sorted(glob.glob(os.path.join(json_dir, "file_*.json")))
    for path in files:
        if os.path.getsize(path) == 0:
            sys.stderr.write(f"skip empty {path}\n")
            continue
        try:
            with open(path) as handle:
                data = json.load(handle)
        except json.JSONDecodeError as exc:
            sys.stderr.write(f"skip {path}: {exc}\n")
            continue
        rows.extend(data.get("items", []) or [])
    return rows


def main() -> None:
    """Parse arguments and write the combined, de-duplicated CSV."""
    parser = argparse.ArgumentParser(
        description="Combine scraped AoU JSON pages into one per-chromosome CSV."
    )
    parser.add_argument("--json-dir", required=True, help="Directory of file_*.json pages")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out_csv))
    os.makedirs(out_dir, exist_ok=True)

    rows = collect_items(args.json_dir)
    if not rows:
        pd.DataFrame(columns=ITEM_COLUMNS).to_csv(args.out_csv, index=False)
        sys.stderr.write(
            f"no items found in {args.json_dir}; wrote empty {args.out_csv}\n"
        )
        return

    df = pd.DataFrame(rows).drop_duplicates(subset=["variantId"])
    df.to_csv(args.out_csv, index=False)
    print(f"wrote {len(df)} variants -> {args.out_csv}")


if __name__ == "__main__":
    main()
