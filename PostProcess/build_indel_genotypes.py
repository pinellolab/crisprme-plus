#!/usr/bin/env python
"""Phased indel genotype store (feature/indel-snp, gated).

The ``log_indels`` file lists each indel's carrier SAMPLES as a FLAT, UNPHASED
list, so it cannot prove that an indel and a nearby SNP sit on the SAME haplotype
(cis). CONFIRMED-cis SNP+indel reporting needs the indel's PHASED genotype per
carrier. The SNP side already has this (Tier-1 genotypes); this module produces
the indel-side equivalent.

It streams the build-time VCF (the same VCFs ``add-variants`` consumed), and for
every INDEL alt (len(REF) != len(alt)) emits, keyed by ``<pos>_<REF>_<ALT>``, the
per-carrier phased genotype NORMALIZED to a biallelic "this-indel vs not" call
(so it is directly comparable to a SNP's biallelic phased GT when checking cis).

STDLIB ONLY (gzip). Runs at build time; the reader loads one chromosome's indels
into a dict (indels are ~10x fewer than SNPs, so this is cheap). Gated behind
CRISPRME_INDEL_SNP; classic build untouched.
"""

import os
import sys
import gzip


def _is_indel(ref, alt):
    """True iff (ref, alt) is an indel (length differs) and not a symbolic allele."""
    return len(ref) != len(alt) and "<" not in alt and "*" not in alt


def normalize_gt_for_alt(gt, alt_idx):
    """Normalize a phased/unphased VCF genotype string to a biallelic
    "this-indel-alt vs everything-else" call, preserving the phase separator.

    ``gt`` is the GT subfield (e.g. "1|0", "0|1", "1|2", "2|2", "1/2", "0|0").
    ``alt_idx`` is the 1-based ALT index of the indel of interest. Each allele
    slot becomes "1" iff it equals ``alt_idx``, else "0". Returns the normalized
    2-slot string, or ``None`` if the sample carries the indel on NO slot (so it
    is not a carrier of THIS indel).

    Examples (alt_idx=1): "1|0"->"1|0", "0|1"->"0|1", "1|1"->"1|1", "0|0"->None.
    Examples (alt_idx=2): "1|2"->"0|1", "2|2"->"1|1", "1|0"->None, "1/2"->"0/1".
    A '/' (unphased) separator is preserved so the downstream cis check forces
    PUTATIVE (never CONFIRMED) for unphased samples.
    """
    gt = gt.split(":", 1)[0]  # GT is always the first FORMAT subfield
    if "|" in gt:
        sep, slots = "|", gt.split("|")
    elif "/" in gt:
        sep, slots = "/", gt.split("/")
    else:
        slots = [gt]  # haploid (e.g. chrY / male chrX)
        sep = ""
    tgt = str(alt_idx)
    norm = ["1" if s == tgt else "0" for s in slots]
    if "1" not in norm:
        return None  # not a carrier of this alt
    return sep.join(norm) if sep else norm[0]


def iter_indel_genotypes(vcf_path, keep_samples=None):
    """Stream a (b)gzipped or plain VCF and yield, per indel alt:
    ``(pos:int, ref:str, alt:str, {sample_id: normalized_gt})``.

    Multiallelic sites emit one tuple per indel ALT (SNP alts are skipped).
    ``keep_samples`` (a set) restricts the emitted carriers to that panel; None
    keeps all. Only carriers (normalize_gt_for_alt != None) are stored.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        sample_cols = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_cols = line.rstrip("\n").split("\t")[9:]
                continue
            if sample_cols is None:
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            pos, ref, alt_field = c[1], c[3], c[4]
            gts = c[9:]
            for a_i, alt in enumerate(alt_field.split(","), start=1):
                if not _is_indel(ref, alt):
                    continue
                carriers = {}
                for col, gt in enumerate(gts):
                    sid = sample_cols[col]
                    if keep_samples is not None and sid not in keep_samples:
                        continue
                    ng = normalize_gt_for_alt(gt, a_i)
                    if ng is not None:
                        carriers[sid] = ng
                if carriers:
                    yield int(pos), ref, alt, carriers


def compile_indel_genotypes(vcf_path, out_path, keep_samples=None):
    """Write the phased indel genotype store for one chromosome's VCF.

    Format (gzipped TSV, one line per indel alt):
        <pos>_<REF>_<ALT> \t sample:gt,sample:gt,...
    Returns the number of indel records written.
    """
    n = 0
    with gzip.open(out_path, "wt") as out:
        for pos, ref, alt, carriers in iter_indel_genotypes(vcf_path, keep_samples):
            toks = ",".join(f"{s}:{g}" for s, g in carriers.items())
            out.write(f"{pos}_{ref}_{alt}\t{toks}\n")
            n += 1
    return n


class IndelGenotypeReader(object):
    """Read a phased indel genotype store (one chromosome) into a dict keyed by
    ``<pos>_<REF>_<ALT>``. Small enough to hold in memory (indels << SNPs)."""

    __slots__ = ("_by_key",)

    def __init__(self, store_path):
        self._by_key = {}
        opener = gzip.open if store_path.endswith(".gz") else open
        with opener(store_path, "rt") as fh:
            for line in fh:
                key, _, toks = line.rstrip("\n").partition("\t")
                if not key:
                    continue
                self._by_key[key] = toks

    def carriers(self, pos, ref, alt):
        """Return "sample:gt" tokens (comma-joined) for the indel, or "" if absent."""
        return self._by_key.get(f"{pos}_{ref}_{alt}", "")

    def carriers_dict(self, pos, ref, alt):
        """Return {sample_id: normalized_gt} for the indel (empty if absent)."""
        toks = self.carriers(pos, ref, alt)
        out = {}
        for t in toks.split(","):
            if not t:
                continue
            s, _, g = t.partition(":")
            out[s] = g
        return out


def _load_panel(samplesid_path):
    """Load the set of panel sample IDs from a samplesID file (skip header)."""
    keep = set()
    opener = gzip.open if samplesid_path.endswith(".gz") else open
    with opener(samplesid_path, "rt") as fh:
        for i, line in enumerate(fh):
            s = line.strip()
            if not s or (i == 0 and (s.lower().startswith("sample") or "\t" in s)):
                # tolerate a header row; otherwise treat the first token as an ID
                pass
            keep.add(s.split("\t")[0])
    keep.discard("SAMPLE_ID")
    return keep


def main(argv):
    """CLI: build_indel_genotypes.py <vcf> <out_store.tsv.gz> [<samplesID>]

    Gated: no-op unless CRISPRME_INDEL_SNP is set.
    """
    if os.environ.get("CRISPRME_INDEL_SNP", "0") not in ("1", "true", "True", "yes"):
        return 0
    if len(argv) < 3:
        sys.stderr.write("usage: build_indel_genotypes.py <vcf> <out.tsv.gz> [samplesID]\n")
        return 2
    vcf_path, out_path = argv[1], argv[2]
    keep = _load_panel(argv[3]) if len(argv) > 3 and os.path.isfile(argv[3]) else None
    n = compile_indel_genotypes(vcf_path, out_path, keep_samples=keep)
    print(f"build_indel_genotypes: {n} phased indel records -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
