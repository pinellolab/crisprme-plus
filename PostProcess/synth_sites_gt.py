#!/usr/bin/env python3
"""Synthesize a presence-genotype column for a sites-only (aggregate) VCF.

Aggregate / sites-only panels (e.g. the "mega" merge: gnomAD, TOPMed, All-of-Us)
carry per-site allele frequencies but NO per-sample genotypes. The CRISPRitz
enricher's ``indel_to_fasta()`` only materializes a fake-indel contig (and its
``log_indels`` row) when at least one VCF *sample column* carries the indel --
it iterates ``line[9:]``. A sites-only VCF (``bcftools view -G``) has just 8
columns, so that gate never passes and the fake-indel genome comes out empty ->
indels are not searchable on the index.

Every record in the mega merge is *observed* (it survived a MAF threshold), so a
blanket presence genotype ``0/1`` on a single synthetic sample is correct and
lets the enricher build the fake-indel genome. REF/ALT/INFO (all the per-dataset
AF / AF_max fields) are preserved verbatim, so the downstream AF annotation is
unchanged.

A ``.vcf.gz`` output is written as **BGZF** (blocked gzip) via pysam, not plain
``gzip`` -- the enricher recipe `tabix`-indexes this file, and tabix REJECTS a
plain-gzip stream ("not BGZF"). (A non-``.gz`` output path is written uncompressed.)

Usage:  synth_sites_gt.py IN.vcf[.gz] OUT.vcf.gz [SAMPLE_NAME]
"""
import gzip
import sys


def _open_read(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _open_write(path):
    """Return (writer, is_bgzf). A ``.gz`` path -> pysam BGZFile (tabix-indexable);
    otherwise a plain text file."""
    if path.endswith(".gz"):
        try:
            from pysam import BGZFile
        except ImportError:  # pragma: no cover
            sys.exit(
                "synth_sites_gt: writing a BGZF .vcf.gz needs pysam "
                "(install pysam, or give a plain .vcf output path)"
            )
        return BGZFile(path, "wb"), True
    return open(path, "w"), False


def synth(inp, out, sample="MEGA"):
    n_records = 0
    w, is_bgzf = _open_write(out)

    def emit(s):
        w.write(s.encode() if is_bgzf else s)

    try:
        with _open_read(inp) as fh:
            for line in fh:
                if line.startswith("##"):
                    emit(line)
                    continue
                if line.startswith("#CHROM"):
                    emit(
                        '##FORMAT=<ID=GT,Number=1,Type=String,'
                        'Description="Synthesized presence genotype for a sites-only '
                        'aggregate panel (every record is observed by MAF threshold)">\n'
                    )
                    cols = line.rstrip("\n").split("\t")
                    # sites-only has 8 fixed columns (no FORMAT/sample); add both.
                    emit("\t".join(cols[:8] + ["FORMAT", sample]) + "\n")
                    continue
                if not line.strip():
                    continue
                row = line.rstrip("\n").split("\t")
                # keep the 8 fixed columns verbatim; append FORMAT=GT + presence 0/1
                emit("\t".join(row[:8] + ["GT", "0/1"]) + "\n")
                n_records += 1
    finally:
        w.close()
    return n_records


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: synth_sites_gt.py IN.vcf[.gz] OUT.vcf.gz [SAMPLE_NAME]")
    _sample = sys.argv[3] if len(sys.argv) > 3 else "MEGA"
    _n = synth(sys.argv[1], sys.argv[2], _sample)
    print("synth_sites_gt: wrote %d records with a '%s' presence genotype -> %s"
          % (_n, _sample, sys.argv[2]))
