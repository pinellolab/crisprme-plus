#!/usr/bin/env python3
"""
Turnkey driver: (re)generate the brute-force ground-truth reference(s) for every
benchmark in benchmarks.json, given a per-chromosome variant-enriched FASTA.

This is how new benchmark examples are produced/cached: add an entry to
benchmarks.json, run this driver on the enriched genome, host the resulting TSV,
and paste its printed MD5 back into the registry.

Usage:
    python generate_references.py --enriched-fasta chr22.enriched.fa --chrom chr22 \
        [--only cas12a_hbg] [--out-dir .]

Notes:
- The brute-force generator is a pure-Python exhaustive scanner (~minutes per Mb),
  so a full chromosome is a batch job (run on a server / in CI, not interactively).
- Run once per chromosome and concatenate the per-chromosome TSVs to build a
  genome-wide reference, then record the combined MD5 in benchmarks.json.
- The enriched FASTA MUST be the same genome CRISPRme searches (its add-variants
  output); a differently-enriched FASTA yields spurious hits in variant-dense
  regions.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enriched-fasta", required=True,
                    help="Per-chromosome variant-enriched FASTA (CRISPRme add-variants output)")
    ap.add_argument("--chrom", required=True, help="Chromosome name, e.g. chr22")
    ap.add_argument("--only", default=None, help="Generate only this benchmark 'name'")
    ap.add_argument("--out-dir", default=HERE, help="Directory for output TSVs")
    args = ap.parse_args()

    registry = json.load(open(os.path.join(HERE, "benchmarks.json")))
    global_th = registry.get("thresholds", {"mm": 4, "bDNA": 1, "bRNA": 1})
    generator = os.path.join(HERE, "generate_brute_force.py")

    for bench in registry["benchmarks"]:
        if args.only and bench["name"] != args.only:
            continue
        # Per-case thresholds mirror complete_test.py so the ground truth matches the
        # exact search this case runs (per-type budgets + optional single-n total cap).
        th = dict(global_th)
        th.update(bench.get("thresholds", {}))
        edits_cap = th.get("max_total_edits")  # None -> independent per-type budgets
        # When a total cap binds, no in-budget alignment can exceed it in any single
        # type, so capping each per-type budget at the total keeps the DP result
        # identical yet much faster (e.g. mm 6 with cap 3 explores only mm<=3).
        mm = min(th["mm"], edits_cap) if edits_cap is not None else th["mm"]
        bdna = min(th["bDNA"], edits_cap) if edits_cap is not None else th["bDNA"]
        brna = min(th["bRNA"], edits_cap) if edits_cap is not None else th["bRNA"]
        out = os.path.join(args.out_dir, f"{bench['name']}.{args.chrom}.tsv")
        cmd = [
            sys.executable, generator,
            "--fasta", args.enriched_fasta,
            "--rna", bench["guide_search"],
            "--max-mismatches", str(mm),
            "--max-dna-gaps", str(bdna),
            "--max-rna-gaps", str(brna),
            "--chrom", args.chrom,
            "--pam", bench.get("pam_concrete", ""),
            "--output", out,
        ]
        if edits_cap is not None:
            cmd += ["--max-total-edits", str(edits_cap)]
        if bench.get("pam_5prime"):
            cmd.append("--pam-5prime")
        sys.stderr.write(f"[{bench['name']}] {bench['nuclease']} {bench['guide_search']} "
                         f"mm={mm} bDNA={bdna} bRNA={brna}"
                         + (f" max-total-edits={edits_cap}" if edits_cap is not None else "")
                         + f" -> {out}\n")
        rc = subprocess.call(cmd)
        if rc != 0:
            sys.stderr.write(f"[{bench['name']}] generation FAILED (exit {rc})\n")
            sys.exit(rc)
        nrows = sum(1 for _ in open(out)) - 1
        sys.stderr.write(f"[{bench['name']}] done: {nrows} rows, md5({os.path.basename(out)})={md5(out)}\n")


if __name__ == "__main__":
    main()
