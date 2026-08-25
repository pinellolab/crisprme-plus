#!/usr/bin/env python
"""Overlay SNP IUPAC codes onto the fake-indel genome (feature/indel-snp, gated).

CRISPRme's indel search runs against a "fake-indel genome" that ``crispritz
add-variants`` cuts from the PLAIN reference (one fake chromosome per real
chromosome; ~25 bp flanks around each indel with the indel ALT applied at
offset 25). Because those flanks are plain reference, an off-target that needs
BOTH an indel AND a nearby SNP is invisible to the indel search.

This module OVERLAYS the SNP ambiguity (IUPAC) codes from the already-produced
enriched genome (``Genomes/<ref>+<vcf>/<chrom>.enriched.fa``) onto the flank
(non-indel) bases of each fake sub-sequence, using the ``log_indels`` fake<->real
mapping. IUPAC is compact (1 char per SNP), so record LENGTHS are preserved and
the existing ``FAKEPOS`` / ``INDELS_tree`` remap stays valid; the CRISPRitz index
then matches SNP+indel targets for free (searchTST matches IUPAC codes).

Run AFTER ``add-variants`` and BEFORE indexing the ``_INDELS`` genome. No-op unless
``CRISPRME_INDEL_SNP`` is truthy (the feature gate).

Geometry reference: CRISPRitz enricher.py ``indel_to_fasta`` (25 bp left flank,
REF->ALT re.sub at offset 25, ``refseq``/``FAKEPOS`` = matched-length window,
``start_position = POS-26`` in the enricher's 0-based ``genomeStr`` coordinate).
"""

import os
import sys
import gzip

# The 15 non-ACGT / non-N IUPAC ambiguity codes the enricher writes for SNPs. A
# base is "SNP-carrying" in the enriched genome iff it is one of these.
_IUPAC_AMBIG = set("RYSWKMBDHVryswkmbdhv")

# The enricher's fixed left-flank width (enricher.py:302, start_position = POS-26,
# indel applied at sub_fasta[25:]).
LEFT_FLANK = 25


def map_fake_offset_to_real(k, start_position, ref_len, alt_len):
    """Map a 0-based offset ``k`` within a fake sub-sequence to its 0-based real
    reference position (enricher ``genomeStr`` coordinate), or ``None`` if ``k`` is
    an indel/ALT base with no reference position.

    Mirrors enricher.indel_to_fasta: ``sub_fasta = leftflank(25 plain-ref) +
    resub(REF->ALT applied)``.
      * k in [0, 25)                 -> LEFT FLANK:  start_position + k
      * k in [25, 25 + alt_len)      -> ALT/indel bases: None (skip; it's the variant)
      * k in [25 + alt_len, len)     -> DOWNSTREAM ref: start_position + k + (ref_len - alt_len)

    The downstream shift is delta = ref_len - alt_len (positive for deletions:
    downstream real positions run ahead of the fake offset; negative for
    insertions). Verified against real chr22 log rows (AT>A deletion, G>GATAA
    insertion) in test_indel_snp_overlay.py.

    COORDINATE CONVENTION (validated on real hg38_1000G chr22 -- see
    validate_indel_snp_coords.py): the returned ``real`` is a 0-based index into
    the enriched/reference sequence as read by read_enriched_chromosome (i.e.
    enriched[real] is the base). The per-sample SNP dict (my_dict_<chrom>.json) is
    keyed 1-BASED, so the Phase-3 indel post-analysis must look a SNP up at dict
    key ``real + 1`` (equivalently, tier0_registry.retrieve_5tuple uses chr_pos+1
    internally, so pass ``real`` as its 0-based chr_pos).
    """
    if k < LEFT_FLANK:
        return start_position + k
    if k < LEFT_FLANK + alt_len:
        return None  # the ALT allele (indel) bases — not a reference position
    return start_position + k + (ref_len - alt_len)


def overlay_sub_fasta(sub_fasta, start_position, ref, alt, enriched, enriched_offset=0):
    """Return ``sub_fasta`` with SNP IUPAC codes overlaid on its flank bases.

    For every offset whose real position carries an IUPAC ambiguity code in the
    ``enriched`` chromosome (and only where the fake base currently equals a plain
    nucleotide, so we never clobber the indel ALT or an ``N`` pad), substitute the
    IUPAC code. Length is unchanged (1 char -> 1 char).

    ``enriched_offset`` lets ``enriched`` be a SLICE rather than the whole
    chromosome: ``enriched[i]`` holds the base at absolute position
    ``i + enriched_offset`` (default 0 = full chromosome). Used by the CI fixture
    oracle, which ships only each indel's small enriched window.
    """
    ref_len, alt_len = len(ref), len(alt)
    n_enr = len(enriched) + enriched_offset
    out = list(sub_fasta)
    for k in range(len(sub_fasta)):
        # only overlay onto plain nucleotides of the fake flanks; leave ALT bases,
        # any existing ambiguity, and N pads untouched
        base = out[k]
        if base not in "ACGTacgt":
            continue
        real = map_fake_offset_to_real(k, start_position, ref_len, alt_len)
        if real is None or real < enriched_offset or real >= n_enr:
            continue
        e = enriched[real - enriched_offset]
        if e in _IUPAC_AMBIG:
            out[k] = e
    return "".join(out)


def _open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_enriched_chromosome(enriched_fa):
    """Read an enriched ``<chrom>.enriched.fa`` into a single 0-based sequence
    string (header stripped, newlines removed, uppercased) — same convention as
    the enricher's ``genomeStr`` and new_simple_analysis's ``genomeStr``."""
    with _open_text(enriched_fa) as fh:
        fh.readline()  # header (>chrom)
        return "".join(line.strip() for line in fh).upper()


def parse_log_indels(log_path):
    """Yield (start_position, ref, alt, fake_start, fake_end) per indel log row.

    log_indels columns: CHR(<chrom>_<start>-<end>_<id>), SAMPLES, rsID, AF,
    indel(<chrom>_<pos>_<REF>_<ALT>), FAKEPOS(<start>,<end>), refseq.
    ``start`` in CHR is the enricher's ``start_position`` (= POS-26).
    """
    with _open_text(log_path) as fh:
        header = fh.readline()  # CHR SAMPLES rsID AF indel FAKEPOS refseq
        if not header.startswith("CHR"):
            # no header (older layout): rewind by reprocessing this line too
            fh.seek(0)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6 or not c[0]:
                continue
            try:
                # CHR = <chrom>_<start>-<end>_<id>  -> start
                start_position = int(c[0].rsplit("_", 1)[0].split("-")[0].split("_")[-1])
                # indel = <chrom>_<pos>_<REF>_<ALT>
                ind = c[4].split("_")
                ref, alt = ind[-2], ind[-1]
                fs, fe = c[5].split(",")
                yield start_position, ref, alt, int(fs), int(fe)
            except (ValueError, IndexError):
                continue


def overlay_fake_chromosome(fake_fa, enriched_fa, log_path, out_fa=None):
    """Overlay SNP IUPAC codes onto one fake chromosome FASTA in place (or to
    ``out_fa``). Preserves the exact record/line structure and all lengths; only
    substitutes flank base identities. Returns the number of bases changed."""
    enriched = read_enriched_chromosome(enriched_fa)
    records = list(parse_log_indels(log_path))

    # Read the fake FASTA as header + one flat sequence (newlines removed). The
    # fake genome is a single record: sub_fasta_1 N sub_fasta_2 N ... with FAKEPOS
    # indexing into this flat (N-inclusive) sequence.
    with _open_text(fake_fa) as fh:
        header = fh.readline().rstrip("\n")
        flat = "".join(line.strip() for line in fh)
    seq = list(flat)

    changed = 0
    for start_position, ref, alt, fs, fe in records:
        sub = "".join(seq[fs:fe])
        overlaid = overlay_sub_fasta(sub, start_position, ref, alt, enriched)
        for i, ch in enumerate(overlaid):
            if seq[fs + i] != ch:
                seq[fs + i] = ch
                changed += 1

    out_fa = out_fa or fake_fa
    with open(out_fa, "w") as out:
        out.write(header + "\n")
        # one base per position; write the whole sequence on a single line (CRISPRitz
        # index strips newlines, so line wrapping is irrelevant to correctness)
        out.write("".join(seq) + "\n")
    return changed


def main(argv):
    """CLI: overlay_indel_snps.py <indels_genome_dir> <enriched_genome_dir> <log_indels_dir>

    Gated: no-op unless CRISPRME_INDEL_SNP is truthy. For each fake<chrom>.fa under
    <indels_genome_dir>/fake_<vcf>_<chrom>/, overlay using the matching enriched
    <chrom>.enriched.fa and log_indels/log<chrom>.txt[.gz].
    """
    if os.environ.get("CRISPRME_INDEL_SNP", "0") not in ("1", "true", "True", "yes"):
        return 0
    if len(argv) < 4:
        sys.stderr.write(
            "usage: overlay_indel_snps.py <indels_genome_dir> <enriched_dir> <log_indels_dir>\n"
        )
        return 2
    indels_dir, enriched_dir, log_dir = argv[1], argv[2], argv[3]
    total = 0
    # Progress goes to STDOUT: this runs inside the crisprme build where a subprocess
    # writing to stderr is treated as fatal (see crisprme.py build-index-only).
    for entry in sorted(os.listdir(indels_dir)):
        sub = os.path.join(indels_dir, entry)
        if not os.path.isdir(sub) or not entry.startswith("fake_"):
            continue
        chrom = entry.rsplit("_", 1)[-1]  # fake_<vcf>_<chrom> -> <chrom>
        fake_fa = os.path.join(sub, f"fake{chrom}.fa")
        enriched_fa = os.path.join(enriched_dir, f"{chrom}.enriched.fa")
        log_path = None
        for cand in (f"log{chrom}.txt", f"log{chrom}.txt.gz"):
            p = os.path.join(log_dir, cand)
            if os.path.isfile(p):
                log_path = p
                break
        if not (os.path.isfile(fake_fa) and os.path.isfile(enriched_fa) and log_path):
            print(f"overlay_indel_snps: skip {chrom} (missing inputs)", flush=True)
            continue
        n = overlay_fake_chromosome(fake_fa, enriched_fa, log_path)
        total += n
        print(f"overlay_indel_snps: {chrom} -> {n} SNP overlays", flush=True)
    print(f"overlay_indel_snps: {total} total SNP overlays", flush=True)
    return total


if __name__ == "__main__":
    sys.exit(main(sys.argv))
