"""Fix the Tier-0 registry AF denominator by removing phantom panel samples.

Background (see the dictless-vs-dict correctness audit): the shipped Tier-0
registry was compiled with the OVER-LISTED 1000G samplesID (3500 rows) even
though the phased 1000G VCF only genotypes 2548 of them. The 952 "phantom"
samples are listed in the panel but never appear in any dict, so the panel-AN
compiler (``aggregate_record_panel``) counts each of them as a called hom-ref
for EVERY record. Result: every group's AN is inflated by a CONSTANT (the
phantom ploidy sum for that group), so AF = AC / AN comes out uniformly too low
(GLOBAL autosomal AN 8858 instead of the correct 2*3477 = 6954). AC, n_carrier
and n_hom are UNAFFECTED (phantoms are never carriers).

Why a transform instead of a full dict reparse: reparsing the ~300 GB of
per-chromosome JSON dicts is memory-bandwidth bound and thrashes under
parallelism. But the fix is purely a denominator correction, and because the
phantom contribution is a per-group CONSTANT, we can patch the existing
registry directly:

  correct_AN[gid]      = existing_AN[gid]      - delta_an[gid]
  correct_n_called[gid] = existing_n_called[gid] - delta_nc[gid]

where delta_an[gid] / delta_nc[gid] are the phantom baselines. We DERIVE those
deltas empirically from the one chromosome we already reparsed correctly
(``reg_chr1``): for any (pos, alt) present in both the existing and the correct
registry, ``existing_AN[gid] - correct_AN[gid]`` is exactly delta_an[gid] --
identical for every record because the carrier-missingness reduction is the same
in both (same carriers, same genotypes). We assert that invariant, then apply
the deltas to the remaining AUTOSOMES (chr2..chr22 share the diploid phantom
baseline). chrX has a different ploidy model (haploid males) so its phantom
baseline differs; chrX is handled separately (reparsed solo), NOT by this
autosomal transform.

The patch operates on the group-entry blob linearly: every entry is a fixed
width record ``<code><AC><AN><n_carrier><n_hom><n_called>`` and the delta depends
ONLY on the group code, so we walk the blob and rewrite the AN and n_called
fields in place. Everything else (record array, string pool, manifest, group
codes, field widths) is byte-identical, so the output is a drop-in registry.

STDLIB ONLY (mmap/struct/json/argparse) -- same contract as tier0_registry.
"""

from __future__ import annotations

import argparse
import json
import mmap
import struct
import sys

import tier0_registry as t0r
from tier0_registry import (
    GLOBAL_GROUP_ID,
    SEP,
    _HEADER_STRUCT,
    _RECORD_STRUCT,
    HEADER_SIZE,
    RECORD_SIZE,
    MAGIC,
    _group_struct,
)


# --------------------------------------------------------------------------- #
# Delta derivation (existing vs a known-correct registry for the same chrom)
# --------------------------------------------------------------------------- #
def derive_deltas(existing_bin, existing_idx, correct_bin, correct_idx,
                  sample_limit=None):
    """Derive per-group (delta_an, delta_nc) from existing vs correct registry.

    For every (pos, alt) present in BOTH registries and every group present in
    both, computes:
        d_an = existing.AN - correct.AN
        d_nc = existing.n_called_indiv - correct.n_called_indiv
    and asserts that:
      * AC, n_carrier_indiv, n_hom_indiv are IDENTICAL (carriers must match --
        the whole premise is that phantoms are non-carriers),
      * d_an and d_nc are CONSTANT per group across all sampled records.

    Returns (delta_an, delta_nc): dict[gid] -> int. Groups seen only in one
    registry (e.g. a group with a carrier here but not there) are skipped for
    derivation but MUST still receive a delta -- see ``complete_deltas``.

    ``sample_limit`` caps how many correct-registry records are scanned (None =
    all). A few thousand already cover every group; we default to all for a
    genome-scale proof but allow capping for speed.
    """
    er = t0r.RegistryReader(existing_bin, existing_idx)
    cr = t0r.RegistryReader(correct_bin, correct_idx)
    delta_an = {}
    delta_nc = {}
    mismatches = []
    n = len(cr)
    scanned = 0
    try:
        for i in range(n):
            pos, alt = cr.record_key_at(i)
            cgroups = cr.lookup(pos, alt)
            egroups = er.lookup(pos, alt)
            if egroups is None:
                # record missing in existing registry -- cannot derive here.
                continue
            for gid, cc in cgroups.items():
                ec = egroups.get(gid)
                if ec is None:
                    continue
                # carriers MUST be identical -- this is the invariant that
                # proves the only difference is the phantom denominator.
                if (ec.AC != cc.AC or ec.n_carrier_indiv != cc.n_carrier_indiv
                        or ec.n_hom_indiv != cc.n_hom_indiv):
                    mismatches.append(
                        (pos, alt, gid,
                         (ec.AC, ec.n_carrier_indiv, ec.n_hom_indiv),
                         (cc.AC, cc.n_carrier_indiv, cc.n_hom_indiv)))
                    if len(mismatches) > 20:
                        raise ValueError(
                            "carrier counts differ between existing and correct "
                            "registry -- the difference is NOT purely the AN "
                            "denominator; refusing to transform. First 20:\n" +
                            "\n".join(repr(m) for m in mismatches[:20]))
                    continue
                d_an = ec.AN - cc.AN
                d_nc = ec.n_called_indiv - cc.n_called_indiv
                if gid in delta_an:
                    if delta_an[gid] != d_an or delta_nc[gid] != d_nc:
                        raise ValueError(
                            "delta not constant for group %r: had (an=%d,nc=%d) "
                            "now (an=%d,nc=%d) at pos=%d alt=%s -- the phantom "
                            "contribution must be constant; aborting"
                            % (gid, delta_an[gid], delta_nc[gid], d_an, d_nc,
                               pos, alt))
                else:
                    delta_an[gid] = d_an
                    delta_nc[gid] = d_nc
            scanned += 1
            if sample_limit is not None and scanned >= sample_limit:
                break
    finally:
        er.close()
        cr.close()
    if mismatches:
        raise ValueError(
            "carrier mismatches found (%d) -- aborting:\n%s"
            % (len(mismatches), "\n".join(repr(m) for m in mismatches[:20])))
    return delta_an, delta_nc


def summarize_deltas(delta_an, delta_nc):
    """Human-readable summary + sanity checks. Returns a list of warning strs."""
    warnings = []
    for gid in sorted(delta_an):
        # HGDP groups must have ZERO phantom delta (phantoms are all 1000G).
        if gid == "HGDP" or gid.startswith("HGDP" + SEP):
            if delta_an[gid] != 0 or delta_nc[gid] != 0:
                warnings.append(
                    "UNEXPECTED non-zero delta on HGDP group %r: an=%d nc=%d "
                    "(phantoms should be 1000G-only)"
                    % (gid, delta_an[gid], delta_nc[gid]))
    return warnings


# --------------------------------------------------------------------------- #
# Transform: patch AN + n_called in the group blob by group-code delta
# --------------------------------------------------------------------------- #
def transform_registry(in_bin, in_idx, delta_an, delta_nc, out_bin, out_idx):
    """Write a corrected registry: AN -= delta_an[gid], n_called -= delta_nc[gid].

    Walks the group-entry blob linearly (every entry is fixed-width and the delta
    depends only on the entry's group code), rewriting the AN (2nd count) and
    n_called (5th count) fields in place. AC/n_carrier/n_hom and everything else
    are preserved byte-for-byte, and the manifest (group codes, widths, offsets,
    string pool) is copied unchanged except an ``an_denominator_fixed`` note.

    Guards: refuses to write a group whose corrected AN < AC or AN <= 0 (would be
    a corrupt frequency), or an unknown group code.
    """
    with open(in_idx) as fh:
        manifest = json.load(fh)
    count_width = manifest["count_width"]
    code_width = manifest["code_width"]
    group_size = manifest["group_size"]
    group_codes = manifest["group_codes"]           # code index -> gid
    group_blob_off = manifest["group_blob_off"]
    string_pool_off = manifest["string_pool_off"]
    grp_struct = _group_struct(count_width, code_width)
    if grp_struct.size != group_size:
        raise ValueError("group_size mismatch: manifest %d vs struct %d"
                         % (group_size, grp_struct.size))

    # Map code -> (d_an, d_nc); a code (group) with no delta defaults to 0.
    code_delta = []
    for gid in group_codes:
        code_delta.append((delta_an.get(gid, 0), delta_nc.get(gid, 0)))

    with open(in_bin, "rb") as fh:
        buf = bytearray(fh.read())

    blob_len = string_pool_off - group_blob_off
    if blob_len % group_size != 0:
        raise ValueError("group blob length %d not a multiple of group_size %d"
                         % (blob_len, group_size))
    n_entries = blob_len // group_size

    patched = 0
    for k in range(n_entries):
        off = group_blob_off + k * group_size
        code, AC, AN, ncar, nhom, ncall = grp_struct.unpack_from(buf, off)
        if code >= len(code_delta):
            raise ValueError("group code %d out of range (%d groups) at entry %d"
                             % (code, len(code_delta), k))
        d_an, d_nc = code_delta[code]
        if d_an == 0 and d_nc == 0:
            continue
        new_an = AN - d_an
        new_ncall = ncall - d_nc
        if new_an <= 0 or new_an < AC or new_ncall < 0 or new_ncall < ncar:
            raise ValueError(
                "corrupt correction at entry %d code %d (gid=%s): AC=%d AN=%d->%d "
                "n_carrier=%d n_called=%d->%d"
                % (k, code, group_codes[code], AC, AN, new_an, ncar, ncall,
                   new_ncall))
        grp_struct.pack_into(buf, off, code, AC, new_an, ncar, nhom, new_ncall)
        patched += 1

    with open(out_bin, "wb") as fh:
        fh.write(buf)
    manifest = dict(manifest)
    manifest["an_denominator_fixed"] = True
    manifest["an_denominator_note"] = (
        "AN/n_called corrected by removing phantom panel samples "
        "(registry_fix_an transform)")
    with open(out_idx, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return {"n_entries": n_entries, "patched": patched}


# --------------------------------------------------------------------------- #
# Validation: transformed must equal a reference (correct) registry exactly
# --------------------------------------------------------------------------- #
def registries_equal(a_bin, a_idx, b_bin, b_idx, max_report=20):
    """Return (ok, diffs). Compares record keys + per-group Counts + rsid + ref.

    Used to prove transform(existing_chrN) == reparsed_correct_chrN.
    """
    ar = t0r.RegistryReader(a_bin, a_idx)
    br = t0r.RegistryReader(b_bin, b_idx)
    diffs = []
    try:
        if len(ar) != len(br):
            diffs.append("n_records: %d vs %d" % (len(ar), len(br)))
        n = min(len(ar), len(br))
        for i in range(n):
            ak = ar.record_key_at(i)
            bk = br.record_key_at(i)
            if ak != bk:
                diffs.append("key[%d]: %r vs %r" % (i, ak, bk))
                if len(diffs) >= max_report:
                    break
                continue
            pos, alt = ak
            ag = ar.lookup(pos, alt)
            bg = br.lookup(pos, alt)
            if set(ag) != set(bg):
                diffs.append("groups@%d,%s: %r vs %r"
                             % (pos, alt, sorted(ag), sorted(bg)))
            else:
                for gid in ag:
                    if ag[gid].as_tuple() != bg[gid].as_tuple():
                        diffs.append("counts@%d,%s,%s: %r vs %r"
                                     % (pos, alt, gid, ag[gid].as_tuple(),
                                        bg[gid].as_tuple()))
            if ar.rsid(pos, alt) != br.rsid(pos, alt):
                diffs.append("rsid@%d,%s: %r vs %r"
                             % (pos, alt, ar.rsid(pos, alt), br.rsid(pos, alt)))
            if ar.ref(pos, alt) != br.ref(pos, alt):
                diffs.append("ref@%d,%s" % (pos, alt))
            if len(diffs) >= max_report:
                break
    finally:
        ar.close()
        br.close()
    return (len(diffs) == 0, diffs)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _cli(argv=None):
    p = argparse.ArgumentParser(prog="registry_fix_an")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("derive", help="derive + print deltas from chr1")
    d.add_argument("--existing-bin", required=True)
    d.add_argument("--existing-idx", required=True)
    d.add_argument("--correct-bin", required=True)
    d.add_argument("--correct-idx", required=True)
    d.add_argument("--out-json", required=True, help="write deltas here")
    d.add_argument("--sample-limit", type=int, default=None)

    t = sub.add_parser("transform", help="apply deltas to a registry")
    t.add_argument("--in-bin", required=True)
    t.add_argument("--in-idx", required=True)
    t.add_argument("--deltas", required=True, help="json from derive")
    t.add_argument("--out-bin", required=True)
    t.add_argument("--out-idx", required=True)

    v = sub.add_parser("validate", help="assert two registries are identical")
    v.add_argument("--a-bin", required=True)
    v.add_argument("--a-idx", required=True)
    v.add_argument("--b-bin", required=True)
    v.add_argument("--b-idx", required=True)

    args = p.parse_args(argv)
    if args.cmd == "derive":
        da, dn = derive_deltas(args.existing_bin, args.existing_idx,
                               args.correct_bin, args.correct_idx,
                               sample_limit=args.sample_limit)
        warns = summarize_deltas(da, dn)
        with open(args.out_json, "w") as fh:
            json.dump({"delta_an": da, "delta_nc": dn}, fh, indent=2,
                      sort_keys=True)
        for gid in sorted(da):
            print("  %-16s  d_an=%-6d d_nc=%-6d" % (gid, da[gid], dn[gid]))
        for w in warns:
            print("WARNING:", w)
        print("wrote", args.out_json)
        return 0
    if args.cmd == "transform":
        with open(args.deltas) as fh:
            d = json.load(fh)
        stats = transform_registry(args.in_bin, args.in_idx,
                                   d["delta_an"], d["delta_nc"],
                                   args.out_bin, args.out_idx)
        print("transformed:", stats)
        return 0
    if args.cmd == "validate":
        ok, diffs = registries_equal(args.a_bin, args.a_idx,
                                     args.b_bin, args.b_idx)
        if ok:
            print("IDENTICAL")
            return 0
        print("DIFFERENCES (first %d):" % len(diffs))
        for x in diffs:
            print("  ", x)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli())
