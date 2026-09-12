"""COSMIC licence gate for the functional-annotation bundle.

COSMIC (the Catalogue Of Somatic Mutations In Cancer, Genome Research Ltd /
Wellcome Sanger Institute) is free for academic / non-commercial research (with
registration), but **commercial use requires a paid licence** from GRL (via
QIAGEN); COSMIC data older than twelve months is released under CC BY-NC-SA 3.0
(non-commercial). Full terms: https://www.cosmickb.org/terms/ .

CRISPRme+'s built-in hg38 annotation bundle
(``dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz``) bakes COSMIC Cancer Gene
Census regions into a *single* BED; those rows carry a ``_COSMIC`` label suffix
and surface downstream as the ``Annotation_COSMIC`` report column. To stay
licence-safe out of the box, CRISPRme+ **excludes COSMIC from every search by
default** and only includes it after the user explicitly attests, in Settings,
that they hold a COSMIC licence appropriate for their use.

This module is intentionally dependency-free (stdlib only) so it can be imported
by both the CLI (``crisprme.py`` / ``PostProcess``) and the Dash web layer. It
holds the persisted attestation flag and the row-level COSMIC filter used to
strip COSMIC from the active annotation when the licence has not been attested.
"""

import gzip
import json
import os
from typing import Tuple

# stored under <data_dir>/Annotations/
COSMIC_LICENSE_FILE = ".cosmic_license.json"
# the label suffix that marks a COSMIC row in the combined annotation BED
COSMIC_TAG = "_COSMIC"


def license_path(annotations_dir: str) -> str:
    """Path to the per-installation COSMIC attestation file."""
    return os.path.join(annotations_dir, COSMIC_LICENSE_FILE)


def cosmic_enabled(annotations_dir: str) -> bool:
    """True only if the user has explicitly attested a COSMIC licence in Settings.

    Defaults to **False** (COSMIC excluded) when the flag file is missing or
    unreadable — the licence-safe default.
    """
    try:
        with open(license_path(annotations_dir)) as fh:
            return bool(json.load(fh).get("enabled", False))
    except (OSError, ValueError):
        return False


def get_cosmic_license(annotations_dir: str) -> dict:
    """Return the full attestation record ({'enabled', 'attestation', ...})."""
    try:
        with open(license_path(annotations_dir)) as fh:
            data = json.load(fh)
            return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def set_cosmic_license(annotations_dir: str, enabled: bool, attestation: str = "") -> None:
    """Persist the COSMIC attestation atomically (tmp + os.replace)."""
    os.makedirs(annotations_dir, exist_ok=True)
    tmp = license_path(annotations_dir) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(
            {"version": 1, "enabled": bool(enabled), "attestation": attestation or ""},
            fh,
        )
    os.replace(tmp, license_path(annotations_dir))


def _open_maybe_gz(path: str, mode: str):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def is_cosmic_row(line: str) -> bool:
    """A COSMIC row = its BED name column (4th, tab-separated) carries the
    ``_COSMIC`` suffix (e.g. ``TP53_COSMIC``). Falls back to a substring test if
    the row has fewer than 4 columns, so no COSMIC feature can slip through."""
    if COSMIC_TAG not in line:
        return False
    fields = line.rstrip("\n").split("\t")
    if len(fields) >= 4:
        return COSMIC_TAG in fields[3]
    return True


def bed_has_cosmic(bed_path: str) -> bool:
    """Cheap scan: does the annotation BED contain any COSMIC rows?"""
    try:
        with _open_maybe_gz(bed_path, "rt") as fh:
            for line in fh:
                if is_cosmic_row(line):
                    return True
    except OSError:
        pass
    return False


def strip_cosmic(in_bed: str, out_bed: str) -> Tuple[int, int]:
    """Write ``out_bed`` = ``in_bed`` with every COSMIC row removed.

    Comment/track/browser lines and non-COSMIC features are copied through
    verbatim. In/out compression is chosen by the ``.gz`` extension. Returns
    ``(kept, dropped)``.
    """
    kept = dropped = 0
    tmp = out_bed + ".tmp"
    # compression of the temp file must match the FINAL name, not the ".tmp" suffix
    dst_open = gzip.open(tmp, "wt") if out_bed.endswith(".gz") else open(tmp, "wt")
    with _open_maybe_gz(in_bed, "rt") as src, dst_open as dst:
        for line in src:
            if line and line[0] not in "#\n" and is_cosmic_row(line):
                dropped += 1
                continue
            dst.write(line)
            if line.strip() and line[0] != "#":
                kept += 1
    os.replace(tmp, out_bed)
    return kept, dropped
