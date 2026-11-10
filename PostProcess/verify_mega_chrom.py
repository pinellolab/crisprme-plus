#!/usr/bin/env python3
"""Full-coverage verifier of a mega index's stores, against its source VCF read
independently: for each chromosome, the SNP registry and the indel-AF store hold exactly
the alleles CRISPRitz 2.10.0's enricher reads from the merged VCF, with its frequencies;
and, given them, the indel log, the enriched genome and the build input agree.

What the enricher reads (its docs/MULTIALLELIC_ENRICHMENT.md, re-derived here, not taken
from build_mega_registry.py, which this checks): the records whose FILTER is PASS or '.';
a SNP allele is an ALT of one upper-case A, C, G or T over a REF of one base; an indel-path
allele is one of a record whose REF or ALT is longer than one base, of another length or
both longer than one base, not symbolic, at POS LEFT or more (a window from the chromosome's
start: LEFT the window's left flank, 26 with the 25/27 window of CRISPRitz 2.10.0's enricher,
read from the indel log with --log); one REF>ALT on two records is one allele, with the more
frequent record's frequencies (the higher AF_max, the build input's INFO/AF; on a tie, the
earlier): the build input keeps that record alone (prepare_mega_build_vcf.py), as the
enricher, which merges two SNP records of one allele in its dictionaries, gives each indel
record a window and a log row of its own.

Checks, each a FAIL with its first examples:
  * registry: its records are exactly the SNP alleles (every one read back by position and
    ALT: none more, none less), each dataset's AF within 0.5/AN of the VCF's, the GLOBAL
    group's AF the record's AF_max;
  * indel store: its keys are exactly the indel-path alleles, each dataset's AF and AF_max
    the record's;
  * --log DIR (log_indels_<dataset>/): the indel log's alleles are the store's (one row per
    allele: a second record of one allele would be a second row), no SAMPLES listed (the
    presence genotype taken out), each row's AF the allele's AF_max;
  * --enriched DIR --reference DIR (Genomes/<genome>+<dataset>/, Genomes/<genome>/): the
    enriched genome codes, at each position, the IUPAC letter of the reference base and the
    SNP alleles there (Q for all four bases; a reference N or a byte that is no base or
    IUPAC code stays as it is; a reference IUPAC base adds nothing of its own), and differs
    from the reference nowhere else;
  * --build-vcf-dir DIR (VCFs/<dataset>/, mega.<chrom>.build.vcf.gz): the build input holds
    the VCF's PASS and '.' records, one record of each REF>ALT (the more frequent), each
    with INFO/AF its AF_max.

Usage: verify_mega_chrom.py <gw_dir> <reg_dir> <indel_dir> <chrom> [<chrom>...]
          [--log DIR] [--enriched DIR --reference DIR] [--build-vcf-dir DIR]
(<gw_dir>: the folder of the source mega.<chrom>.afmax.vcf.gz.) Prints one PASS/FAIL line
per chromosome and a summary; exit 1 if any FAIL.
"""
import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tier0_registry import RegistryReader, GLOBAL_GROUP_ID, db_group_id  # noqa: E402
from build_mega_indel_af import IndelAfReader  # noqa: E402
from build_mega_registry import DEFAULT_DATASET_META  # noqa: E402

_AN = {d: m["an_nominal"] for d, m in DEFAULT_DATASET_META.items()}
_DS = tuple(DEFAULT_DATASET_META)
_PASS = ("PASS", ".")
_BASES = ("A", "C", "G", "T")   # one base each: a tuple, so "CG" is no member
# the enricher's codes by the set of bases (Q: all four, CRISPRitz 2.10.0)
_IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
          frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
          frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
          frozenset("ACG"): "V", frozenset("ACGT"): "Q"}
_REFERENCE_CODES = set(_BASES) | set("RYSWKMBDHV")
# the enricher's left flank of a fake-indel window (25/27 in CRISPRitz 2.10.0's enricher): no
# window of an indel at POS < LEFT. --log reads it from the indel log (window_left_of_log).
LEFT = 26


def window_left_of_log(log_dir, chrom):
    """The left flank of the enricher's fake-indel windows, as the chromosome's indel log
    shows it: an indel's POS less its window's start (its CHR, <chrom>_<start>-<end>_<n>), or
    None when the log has no row to read."""
    for name in ("log%s.txt.gz" % chrom, "log%s.txt" % chrom):
        path = os.path.join(log_dir, name)
        if not os.path.isfile(path):
            continue
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    start = int(f[ix["CHR"]].rsplit("_", 2)[-2].split("-")[0])
                    return int(f[ix["indel"]].rsplit("_", 3)[1]) - start
                except (KeyError, IndexError, ValueError):
                    continue
    return None


def _info(info):
    out = {}
    for kv in info.split(";"):
        k, _s, v = kv.partition("=")
        out[k] = v
    return out


def _afs(info, k, n):
    """{dataset: AF} of ALT k of n (positive, at most 1), and AF_max (or None)."""
    kv = _info(info)
    afs = {}
    for ds in _DS:
        vals = kv.get("AF_" + ds, "").split(",")
        if len(vals) != n:
            continue
        try:
            v = float(vals[k])
        except ValueError:
            continue
        if 0.0 < v <= 1.0:
            afs[ds] = v
    mx = None
    vals = kv.get("AF_max", "").split(",")
    if len(vals) == n:
        try:
            mx = float(vals[k])
        except ValueError:
            mx = None
        if mx is not None and not 0.0 < mx <= 1.0:
            mx = None
    return afs, mx


def _rank(info, k, n):
    """How two records of one allele are told apart: the higher AF_max (the build input's
    INFO/AF), 0 when it is missing or unreadable, as the enricher ranks INFO/AF; on a tie,
    the earlier."""
    vals = _info(info).get("AF_max", "").split(",")
    try:
        return float(vals[k]) if len(vals) == n else 0.0
    except ValueError:
        return 0.0


def read_source(vcf, left=LEFT):
    """({(pos, REF, ALT): (afs, af_max)} of the SNP alleles, the same of the indel-path
    alleles, the number of PASS/'.' records, the number of them the build input keeps), as
    the enricher reads the VCF: one record of each REF>ALT, the more frequent (_rank),
    gives the allele's frequencies, and a record is left out when another one has each of
    its ALTs; ``left``, the enricher's left flank (no window of an indel at POS < left)."""
    best = {}     # (pos, REF, ALT) -> (rank, record, k, ALTs of the record, INFO)
    alleles = []  # the (pos, REF, ALT) of each PASS/'.' record's ALTs, in record order
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[6] not in _PASS:
                continue
            pos, ref, alts, info = int(f[1]), f[3], f[4].split(","), f[7]
            keys = []
            for k, a in enumerate(alts):
                key = (pos, ref, a)
                keys.append(key)
                rank = _rank(info, k, len(alts))
                if key not in best or rank > best[key][0]:
                    best[key] = (rank, len(alleles), k, alts, info)
            alleles.append(keys)
    snps, indels = {}, {}
    for (pos, ref, a), (_rank_of, _i, k, alts, info) in best.items():
        if len(ref) == 1 and a in _BASES:
            store = snps
        elif ((len(ref) > 1 or any(len(x) > 1 for x in alts)) and "<" not in ref
              and "<" not in a and (len(a) != len(ref) or (len(a) > 1 and len(ref) > 1))
              and pos >= left):
            store = indels
        else:
            continue
        afs, mx = _afs(info, k, len(alts))
        if afs:
            store[(pos, ref, a)] = (afs, mx)
    n_build = sum(1 for i, keys in enumerate(alleles)
                  if any(best[key][1:3] == (i, k) for k, key in enumerate(keys)))
    return snps, indels, len(alleles), n_build


def _registry_keys(reg):
    out = set()
    for i in range(len(reg)):
        out.add(reg.record_key_at(i))
    return out


def _fasta(path):
    with open(path) as fh:
        return "".join(ln.strip() for ln in fh if not ln.startswith(">")).upper()


def check_enriched(enr_path, ref_path, snps, errs):
    """The enriched genome against the reference and the SNP alleles."""
    enr, ref = _fasta(enr_path), _fasta(ref_path)
    if len(enr) != len(ref):
        errs.append("enriched genome %d bases, reference %d" % (len(enr), len(ref)))
        return 0
    alts = {}
    for pos, _r, a in snps:
        alts.setdefault(pos - 1, set()).add(a)
    want = {}
    for p, al in alts.items():
        b = ref[p]
        if b not in _REFERENCE_CODES:      # an N (a gap), or no base: stays
            continue
        bases = set(al) | ({b} if b in _BASES else set())
        if len(bases) >= 2:
            want[p] = _IUPAC[frozenset(bases)]
    bad = 0
    for p, (e, r) in enumerate(zip(enr, ref)):
        w = want.get(p, r)
        if e != w:
            bad += 1
            if bad <= 3:
                errs.append("enriched %d: %s (reference %s, expected %s)" % (p + 1, e, r, w))
    if bad:
        errs.append("enriched genome: %d position(s) differ from the expected code" % bad)
    return len(want)


def check_log(log_dir, chrom, indels, errs):
    """The indel log against the indel-path alleles: one row each, no sample, AF_max."""
    path = None
    for name in ("log%s.txt.gz" % chrom, "log%s.txt" % chrom):
        if os.path.isfile(os.path.join(log_dir, name)):
            path = os.path.join(log_dir, name)
    if path is None:
        errs.append("no indel log for %s in %s" % (chrom, log_dir))
        return 0
    opener = gzip.open if path.endswith(".gz") else open
    seen, rows, sampled, af_bad = {}, 0, 0, 0
    with opener(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows += 1
            _c, pos, r, a = f[ix["indel"]].rsplit("_", 3)
            key = (int(pos), r, a)
            seen[key] = seen.get(key, 0) + 1
            if f[ix["SAMPLES"]].strip():
                sampled += 1
            want = indels.get(key)
            try:
                af = float(f[ix["AF"]])
            except ValueError:
                af = None
            if want is not None and want[1] is not None and (
                    af is None or abs(af - want[1]) > 1e-9):
                af_bad += 1
                if af_bad <= 3:
                    errs.append("log %s AF %s, AF_max %s" % (f[ix["indel"]], f[ix["AF"]],
                                                            want[1]))
    if set(seen) != set(indels):
        more, less = set(seen) - set(indels), set(indels) - set(seen)
        errs.append("indel log alleles %d, expected %d: %d more (%s), %d missing (%s)" % (
            len(seen), len(indels), len(more), sorted(more)[:3], len(less), sorted(less)[:3]))
    twice = sum(1 for n in seen.values() if n > 1)
    if twice:
        errs.append("indel log: %d allele(s) on more than one row" % twice)
    if sampled:
        errs.append("indel log: %d row(s) list samples (the presence genotype left in)"
                    % sampled)
    return rows


def check_build_vcf(path, n_build, errs):
    """The build input: ``n_build`` records (the source's PASS/'.' ones, one of each
    REF>ALT), every one PASS or '.', INFO/AF its AF_max, and no REF>ALT on two of them."""
    if not os.path.isfile(path):
        errs.append("no build input %s" % path)
        return
    n, bad, seen, twice = 0, 0, set(), 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            n += 1
            kv = _info(f[7])
            if f[6] not in _PASS or ("AF_max" in kv and kv.get("AF") != kv["AF_max"]):
                bad += 1
                if bad <= 3:
                    errs.append("build input %s:%s FILTER %s AF %s AF_max %s" % (
                        f[0], f[1], f[6], kv.get("AF"), kv.get("AF_max")))
            for a in f[4].split(","):
                key = (f[1], f[3], a)
                twice += key in seen
                seen.add(key)
    if n != n_build:
        errs.append("build input %d records, the source's PASS/. records less those of a "
                    "REF>ALT another one has %d" % (n, n_build))
    if twice:
        errs.append("build input: %d REF>ALT(s) on more than one record" % twice)


def verify(gw, regd, indd, chrom, log=None, enriched=None, reference=None, build_dir=None):
    vcf = os.path.join(gw, "mega.%s.afmax.vcf.gz" % chrom)
    left = (window_left_of_log(log, chrom) if log else None) or LEFT
    snps, indels, n_pass, n_build = read_source(vcf, left)
    errs = []
    reg = RegistryReader(os.path.join(regd, "reg_%s.bin" % chrom),
                         os.path.join(regd, "reg_%s.idx" % chrom))
    have = _registry_keys(reg)
    want = {(p, a) for p, _r, a in snps}
    if have != want:
        more, less = have - want, want - have
        errs.append("registry %d records, expected %d: %d more (%s), %d missing (%s)" % (
            len(have), len(want), len(more), sorted(more)[:3], len(less), sorted(less)[:3]))
    bad = 0
    for (pos, _r, a), (afs, mx) in snps.items():
        g = reg.lookup(pos, a)
        if g is None:
            continue
        for ds, af in afs.items():
            got = g.get(db_group_id(ds))
            if got is None or abs(got.allele_freq() - af) > 0.5 / _AN[ds] + 1e-9:
                bad += 1
        gg = g.get(GLOBAL_GROUP_ID)
        if mx is not None and (gg is None or abs(gg.allele_freq() - mx) > 1e-3):
            bad += 1
    if bad:
        errs.append("registry: %d frequency(ies) differ from the VCF's" % bad)
    reg.close()

    ind = IndelAfReader(os.path.join(indd, "indel_af_%s.tsv.gz" % chrom))
    have_i = set(ind._d)
    want_i = {(p, r.upper(), a.upper()) for p, r, a in indels}
    if have_i != want_i:
        more, less = have_i - want_i, want_i - have_i
        errs.append("indel store %d, expected %d: %d more (%s), %d missing (%s)" % (
            len(have_i), len(want_i), len(more), sorted(more)[:3], len(less),
            sorted(less)[:3]))
    bad = 0
    for (pos, r, a), (afs, mx) in indels.items():
        rec = ind.lookup(pos, r, a)
        if rec is None:
            continue
        for ds, af in afs.items():
            if rec.get(ds) is None or abs(rec[ds] - af) > 1e-12:
                bad += 1
        if mx is not None and (rec.get("AF_max") is None or abs(rec["AF_max"] - mx) > 1e-12):
            bad += 1
    if bad:
        errs.append("indel store: %d frequency(ies) differ from the VCF's" % bad)

    extra = []
    if log:
        extra.append("log_rows=%d" % check_log(log, chrom, indels, errs))
    if enriched:
        coded = check_enriched(os.path.join(enriched, "%s.enriched.fa" % chrom),
                               os.path.join(reference, "%s.fa" % chrom), snps, errs)
        extra.append("coded_positions=%d" % coded)
    if build_dir:
        check_build_vcf(os.path.join(build_dir, "mega.%s.build.vcf.gz" % chrom), n_build, errs)
    if n_build != n_pass:
        extra.append("second_records_of_one_allele=%d" % (n_pass - n_build))
    if left != LEFT:
        extra.append("window_left=%d" % left)
    ok = not errs
    print("%s: %s  pass_records=%d snp_alleles=%d indel_alleles=%d registry=%d indel_store=%d"
          " %s %s" % (chrom, "PASS" if ok else "FAIL", n_pass, len(snps), len(indels),
                      len(have), len(have_i), " ".join(extra),
                      "" if ok else "| " + " ; ".join(errs[:8])))
    sys.stdout.flush()
    return ok, len(snps), len(indels)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gw_dir")
    ap.add_argument("reg_dir")
    ap.add_argument("indel_dir")
    ap.add_argument("chroms", nargs="+")
    ap.add_argument("--log", default=None)
    ap.add_argument("--enriched", default=None)
    ap.add_argument("--reference", default=None)
    ap.add_argument("--build-vcf-dir", default=None)
    args = ap.parse_args(argv)
    if args.enriched and not args.reference:
        ap.error("--enriched needs --reference")
    allok, tsnp, tind = True, 0, 0
    for c in args.chroms:
        ok, ns, ni = verify(args.gw_dir, args.reg_dir, args.indel_dir, c, args.log,
                            args.enriched, args.reference, args.build_vcf_dir)
        allok = allok and ok
        tsnp += ns
        tind += ni
    print("SUMMARY: %d chroms, %d SNP alleles + %d indel-path alleles, %s"
          % (len(args.chroms), tsnp, tind, "ALL PASS" if allok else "FAILURES"))
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
