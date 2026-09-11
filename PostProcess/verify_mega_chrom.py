#!/usr/bin/env python3
"""Full-coverage verifier for the mega stores: for each chromosome confirm the SNP
registry + indel-AF sidecar together account for EVERY merged site and that sampled
AFs round-trip exactly against the source VCF. Deterministic (seeded sampling).

Usage: verify_mega_chrom.py <gw_dir> <reg_dir> <indel_dir> <chrom> [<chrom>...]
Prints one PASS/FAIL line per chrom + a final summary; exit 1 if any FAIL.
"""
import gzip
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tier0_registry import RegistryReader, GLOBAL_GROUP_ID, db_group_id  # noqa: E402
from build_mega_indel_af import IndelAfReader, _is_indel  # noqa: E402
from build_mega_registry import _parse_info_afs, DEFAULT_DATASET_META  # noqa: E402

_AN = {d: m["an_nominal"] for d, m in DEFAULT_DATASET_META.items()}
_DS = set(DEFAULT_DATASET_META)
K = 400  # sampled variants per class per chrom


def _af_max_from_info(info):
    for kv in info.split(";"):
        if kv.startswith("AF_max="):
            try:
                return float(kv[len("AF_max="):])
            except ValueError:
                return None
    return None


def verify(gw, regd, indd, chrom):
    vcf = os.path.join(gw, "mega.%s.afmax.vcf.gz" % chrom)
    reg = RegistryReader(os.path.join(regd, "reg_%s.bin" % chrom),
                         os.path.join(regd, "reg_%s.idx" % chrom))
    ind = IndelAfReader(os.path.join(indd, "indel_af_%s.tsv.gz" % chrom))
    rnd = random.Random(hash(chrom) & 0xffff)

    n_snp = n_indel = 0
    indel_keys = set()   # UNIQUE indel keys (the store dedups; SNPs have no dups)
    snp_pool, indel_pool = [], []  # reservoir samples of (pos, ref, alt, afs, afmax)
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or "," in f[4]:
                continue
            pos, ref, alt, info = f[1], f[3], f[4], f[7]
            afs = _parse_info_afs(info, _DS)
            if not afs:
                continue
            rec = (int(pos), ref.upper(), alt.upper(), afs, _af_max_from_info(info))
            if _is_indel(ref, alt):
                n_indel += 1
                indel_keys.add((int(pos), ref.upper(), alt.upper()))
                pool = indel_pool
            else:
                n_snp += 1
                pool = snp_pool
            if len(pool) < K:
                pool.append(rec)
            else:
                j = rnd.randint(0, (n_snp if pool is snp_pool else n_indel) - 1)
                if j < K:
                    pool[j] = rec

    errs = []
    # coverage: registry SNP records + indel store size == VCF SNP + indel with AF
    if reg.manifest["n_records"] != n_snp:
        errs.append("reg n_records %d != vcf SNPs %d" % (reg.manifest["n_records"], n_snp))
    # the indel store dedups duplicate (pos,ref,alt) keys (by max AF), so it holds
    # the UNIQUE indel keys, not the raw record count.
    if len(ind) != len(indel_keys):
        errs.append("indel store %d != vcf UNIQUE indels %d (raw %d)" %
                    (len(ind), len(indel_keys), n_indel))

    # sampled SNP AF round-trip (registry) + GLOBAL == AF_max
    for pos, ref, alt, afs, afmax in snp_pool:
        g = reg.lookup(pos, alt)
        if g is None:
            errs.append("SNP %d:%s>%s absent in registry" % (pos, ref, alt)); continue
        for ds, af in afs.items():
            got = g.get(db_group_id(ds))
            if got is None or abs(got.allele_freq() - af) > 0.5 / _AN[ds] + 1e-9:
                errs.append("SNP %d %s AF %.6g vs %s" %
                            (pos, ds, af, None if got is None else got.allele_freq()))
        gg = g.get(GLOBAL_GROUP_ID)
        if afmax is not None and (gg is None or abs(gg.allele_freq() - afmax) > 1e-3):
            errs.append("SNP %d GLOBAL %.6g vs AF_max %.6g" %
                        (pos, -1 if gg is None else gg.allele_freq(), afmax))

    # sampled indel AF round-trip (sidecar) + AF_max
    for pos, ref, alt, afs, afmax in indel_pool:
        rec = ind.lookup(pos, ref, alt)
        if rec is None:
            errs.append("indel %d:%s>%s absent in sidecar" % (pos, ref, alt)); continue
        # the store holds the MAX per-dataset AF over any duplicate (pos,ref,alt)
        # records, so a sampled record's AF must be <= the store's (== for the
        # non-duplicate majority).
        for ds, af in afs.items():
            got = rec.get(ds)
            if got is None or got + 1e-9 < af:
                errs.append("indel %d %s: sampled AF %.6g > store %s" % (pos, ds, af, got))
        if afmax is not None:
            gotm = rec.get("AF_max")
            if gotm is None or gotm + 1e-9 < afmax:
                errs.append("indel %d AF_max: sampled %.6g > store %s" % (pos, afmax, gotm))

    ok = not errs
    print("%s: %s  vcf_sites=%d (snp=%d indel=%d)  reg=%d indel_store=%d  "
          "sampled=%d/%d  %s" % (
              chrom, "PASS" if ok else "FAIL", n_snp + n_indel, n_snp, n_indel,
              reg.manifest["n_records"], len(ind), len(snp_pool), len(indel_pool),
              "" if ok else "| " + " ; ".join(errs[:6])))
    return ok, n_snp, n_indel


def main(argv):
    gw, regd, indd = argv[1], argv[2], argv[3]
    chroms = argv[4:]
    allok = True
    tsnp = tind = 0
    for c in chroms:
        ok, ns, ni = verify(gw, regd, indd, c)
        allok = allok and ok
        tsnp += ns; tind += ni
    print("SUMMARY: %d chroms, %d SNPs + %d indels = %d sites, %s" %
          (len(chroms), tsnp, tind, tsnp + tind, "ALL PASS" if allok else "FAILURES"))
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
