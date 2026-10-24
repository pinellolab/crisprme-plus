#!/usr/bin/env python3
"""Validate (and canonically re-emit) the CFD score tables shipped with CRISPRme+.

The CFD (Cutting Frequency Determination) score of Doench et al., Nat. Biotechnol.
2016 (doi:10.1038/nbt.3437) multiplies a per-position penalty for each spacer
position by a PAM penalty. The penalties are published lookup values:

  * 240 single-mismatch penalties (20 positions x the 12 RNA:DNA mismatch pairs)
    and 16 two-nucleotide PAM penalties come from the Doench 2016 supplementary data.
  * 152 gap/bulge penalties (76 RNA-gap ``rX:d-,p`` + 76 DNA-gap ``r-:dX,p``) follow
    the bulge handling in Fig. 5c,d of the same paper.

These published measurements (data) are committed as two human-reviewable JSON files,
which the scorer (``PostProcess/cfd_score.py``) loads at runtime:

  PostProcess/cfd_mismatch_scores.json   (392 keys: 'r<spacer>:d<target>,<pos>')
  PostProcess/cfd_pam_scores.json        (16 keys: two-nt PAM -> penalty)

JSON (not a pickle) keeps the tables reviewable and avoids an opaque deserialization
step; Python's json writer uses shortest round-trip float repr, so values reload to
the exact same IEEE-754 doubles. The JSON IS the source of truth -- this script does
not derive it from any other artifact; it validates the committed tables' structure,
counts, value range and published anchor values, and (without --check) re-emits them
in canonical sorted form so formatting stays stable.

Usage:
    python scripts/generate_cfd_tables.py            # validate + canonically re-emit
    python scripts/generate_cfd_tables.py --check     # validate only (CI gate)
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(os.path.dirname(HERE), "PostProcess")
MM_JSON = os.path.join(PP, "cfd_mismatch_scores.json")
PAM_JSON = os.path.join(PP, "cfd_pam_scores.json")

# Spot-check values (stable anchors across the two tables).
ANCHORS_MM = {"r-:dA,2": 0.692307692308, "r-:dA,20": 0.6, "rC:d-,2": 0.96875}
ANCHORS_PAM = {"GG": 1.0}
_RNA = set("ACGU")
_DNA = set("ACGT")


def _load_json(path):
    with open(path) as fh:
        return json.load(fh)


def _dump_json(mapping, path):
    with open(path, "w") as fh:
        json.dump(mapping, fh, sort_keys=True, indent=1)
        fh.write("\n")


def _validate(mm, pam):
    # counts
    assert len(mm) == 392, f"expected 392 mismatch keys, got {len(mm)}"
    assert len(pam) == 16, f"expected 16 PAM keys, got {len(pam)}"
    # anchors
    for k, v in ANCHORS_MM.items():
        assert k in mm and abs(mm[k] - v) < 1e-12, f"mismatch anchor {k} != {v} (got {mm.get(k)})"
    for k, v in ANCHORS_PAM.items():
        assert pam.get(k) == v, f"PAM anchor {k} != {v} (got {pam.get(k)})"
    # structure + value range of mismatch keys: 'r<spacer>:d<target>,<pos>'
    single = rna_gap = dna_gap = 0
    for key, val in mm.items():
        assert 0.0 <= float(val) <= 1.0, f"mismatch penalty out of [0,1]: {key}={val}"
        head, pos = key.split(",")
        assert pos.isdigit() and 1 <= int(pos) <= 21, f"bad position in {key}"
        assert head[0] == "r" and head[2:4] == ":d", f"bad key shape {key}"
        spacer, target = head[1], head[4]
        if spacer == "-":
            dna_gap += 1
            assert target in _DNA, f"bad DNA-gap target in {key}"
        elif target == "-":
            rna_gap += 1
            assert spacer in _RNA, f"bad RNA-gap spacer in {key}"
        else:
            single += 1
            assert spacer in _RNA and target in _DNA, f"bad single-mm bases in {key}"
    assert single == 240, f"expected 240 single-mismatch keys, got {single}"
    assert rna_gap == 76 and dna_gap == 76, f"expected 76+76 gap keys, got {rna_gap}+{dna_gap}"
    # PAM keys: two-nt, value in [0,1]
    for key, val in pam.items():
        assert len(key) == 2 and set(key) <= _DNA, f"bad PAM key {key}"
        assert 0.0 <= float(val) <= 1.0, f"PAM penalty out of [0,1]: {key}={val}"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="validate committed JSON only; write nothing")
    args = ap.parse_args(argv)

    mm, pam = _load_json(MM_JSON), _load_json(PAM_JSON)
    _validate(mm, pam)
    # round-trip: the committed JSON reloads to the exact same values after a canonical dump
    assert mm == _load_json(MM_JSON) and pam == _load_json(PAM_JSON)

    if not args.check:
        _dump_json(mm, MM_JSON)   # canonical re-emit (idempotent)
        _dump_json(pam, PAM_JSON)
        _validate(_load_json(MM_JSON), _load_json(PAM_JSON))

    print(f"OK: {len(mm)} mismatch (240 single + 76+76 gap) + {len(pam)} PAM penalties "
          f"validated -> {MM_JSON}, {PAM_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
