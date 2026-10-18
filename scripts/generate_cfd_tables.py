#!/usr/bin/env python3
"""Generate the CFD score tables shipped with CRISPRme+ as reviewable JSON.

The CFD (Cutting Frequency Determination) score of Doench et al.,
Nat. Biotechnol. 2016 (doi:10.1038/nbt.3437) multiplies a per-position penalty
for each spacer position by a PAM penalty. The penalties are published lookup
values:

  * 240 single-mismatch penalties (20 positions x the 12 RNA:DNA mismatch
    pairs) and 16 two-nucleotide PAM penalties come from the Doench 2016
    supplementary data.
  * 152 gap/bulge penalties (an RNA-gap and a DNA-gap entry per applicable
    position) follow the bulge handling in Fig. 5c,d of the same paper.

These are published measurements (data), re-serialized here into two JSON files
that the scorer loads at runtime:

  PostProcess/cfd_mismatch_scores.json   (392 keys: 'r<spacer>:d<target>,<pos>')
  PostProcess/cfd_pam_scores.json        (16 keys: two-nt PAM -> penalty)

JSON (not a pickle) keeps the tables human-reviewable and avoids an opaque
deserialization step. Python's json writer uses shortest round-trip float
repr, so values reload to the exact same IEEE-754 doubles.

Usage:
    python scripts/generate_cfd_tables.py            # write + verify
    python scripts/generate_cfd_tables.py --check     # verify only (CI)

By default this reads the current shipped pickles as the value source and
asserts the emitted JSON reloads bit-for-bit, so the tables are provably the
same validated values in a reviewable format. (A from-archive re-derivation can
replace the pickle source without changing the emitted JSON.)
"""
import argparse
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(os.path.dirname(HERE), "PostProcess")
MM_JSON = os.path.join(PP, "cfd_mismatch_scores.json")
PAM_JSON = os.path.join(PP, "cfd_pam_scores.json")
MM_PKL = os.path.join(PP, "mismatch_score.pkl")
PAM_PKL = os.path.join(PP, "PAM_scores.pkl")

# Spot-check values (stable anchors across the two tables).
ANCHORS_MM = {"r-:dA,2": 0.692307692308, "r-:dA,20": 0.6, "rC:d-,2": 0.96875}
ANCHORS_PAM = {"GG": 1.0}


def _load_pickles():
    with open(MM_PKL, "rb") as fh:
        mm = pickle.load(fh)
    with open(PAM_PKL, "rb") as fh:
        pam = pickle.load(fh)
    return mm, pam


def _dump_json(mapping, path):
    with open(path, "w") as fh:
        json.dump(mapping, fh, sort_keys=True, indent=1)
        fh.write("\n")


def _load_json(path):
    with open(path) as fh:
        return json.load(fh)


def _assert_anchors(mm, pam):
    for k, v in ANCHORS_MM.items():
        assert k in mm and abs(mm[k] - v) < 1e-12, f"mismatch anchor {k} != {v} (got {mm.get(k)})"
    for k, v in ANCHORS_PAM.items():
        assert pam.get(k) == v, f"PAM anchor {k} != {v} (got {pam.get(k)})"


def _assert_identical(a, b, label):
    assert a.keys() == b.keys(), f"{label}: key set differs"
    for k in a:
        # exact IEEE-754 identity (repr round-trip), not approximate
        assert repr(float(a[k])) == repr(float(b[k])), f"{label}: {k} differs ({a[k]!r} vs {b[k]!r})"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="verify existing JSON only; write nothing")
    args = ap.parse_args(argv)

    mm, pam = _load_pickles()
    assert len(mm) == 392, f"expected 392 mismatch keys, got {len(mm)}"
    assert len(pam) == 16, f"expected 16 PAM keys, got {len(pam)}"
    _assert_anchors(mm, pam)

    if not args.check:
        _dump_json(mm, MM_JSON)
        _dump_json(pam, PAM_JSON)

    # emitted JSON reloads to the exact same values
    _assert_identical(mm, _load_json(MM_JSON), "mismatch")
    _assert_identical(pam, _load_json(PAM_JSON), "PAM")
    print(f"OK: {len(mm)} mismatch + {len(pam)} PAM penalties -> {MM_JSON}, {PAM_JSON} (round-trip exact)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
