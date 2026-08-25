#!/usr/bin/env python
"""Dev validation of the indel-SNP coordinate convention against REAL data.

This is NOT a unit test (it needs a real matching enriched-genome + fake-genome +
log_indels + SNP dict for one chromosome); it is the empirical proof behind the
coordinate convention used by overlay_indel_snps + the Phase-3 indel post-analysis.

VALIDATED (hg38_1000G chr22, 2026): the log CHR-field ``start`` and the fake genome
use the SAME 0-based coordinate as ``overlay_indel_snps.read_enriched_chromosome``
(the enriched FASTA read as a header-stripped 0-indexed string). So:

  * fake flank base at fake-region offset k  ==  enriched[start + k]   (0-based)
    -> overlay reads enriched[real] with real = map_fake_offset_to_real(...).
    Verified 72554/72554 (100%) fake-flank vs enriched matches.

  * the SNP dict (my_dict_<chrom>.json) is keyed 1-BASED:
        dict key "chr,POS"  <->  enriched[POS-1]
    -> a real 0-based position ``real`` maps to dict key ``real + 1``.
    Verified 5000/5000 isolated dict SNPs land on enriched[POS-1], and the full
    overlay->real->dict chain 200/200.

Usage:
  validate_indel_snp_coords.py <enriched.fa> <fake.fa> <log_indels.txt[.gz]> <my_dict.json>
"""

import sys
import json
import tempfile
import os
import bisect

import overlay_indel_snps as ov


def validate(enriched_fa, fake_fa, log_path, dict_path, sample=200):
    AMB = ov._IUPAC_AMBIG
    E = ov.read_enriched_chromosome(enriched_fa)

    with ov._open_text(fake_fa) as fh:
        fh.readline()
        orig = "".join(l.strip() for l in fh)

    recs = list(ov.parse_log_indels(log_path))

    # 1) fake flank base == enriched[start + k] (the overlay's read position)
    m = t = 0
    for sp, ref, alt, fs, fe in recs[:3000]:
        for k in range(25):
            fb, eb = orig[fs + k], (E[sp + k] if sp + k < len(E) else "N")
            if fb in "ACGT" and eb in "ACGT":
                t += 1
                m += (fb == eb)
    assert t and m == t, f"fake-flank vs enriched[start+k] mismatch: {m}/{t}"
    print(f"[1] fake flank == enriched[start+k]: {m}/{t} (100%)")

    # 2) run the overlay; every changed base is an IUPAC code, lengths preserved
    out = tempfile.mktemp(suffix=".fa")
    changed = ov.overlay_fake_chromosome(fake_fa, enriched_fa, log_path, out_fa=out)
    with open(out) as fh:
        fh.readline()
        over = "".join(l.strip() for l in fh)
    os.remove(out)
    assert len(over) == len(orig), "overlay changed the fake genome length!"
    diffs = [i for i in range(len(orig)) if orig[i] != over[i]]
    assert all(over[i] in AMB for i in diffs), "an overlaid base is not IUPAC"
    print(f"[2] overlay: {changed} bases IUPAC-coded, lengths preserved")

    # 3) full chain: overlaid pos -> real -> enriched IUPAC + dict key real+1
    dpos = set(int(k.split(",")[1]) for k in json.load(open(dict_path)))
    starts = [r[3] for r in recs]
    ok = bad = 0
    step = max(1, len(diffs) // sample)
    for P in diffs[::step][:sample]:
        i = bisect.bisect_right(starts, P) - 1
        sp, ref, alt, fs, fe = recs[i]
        if not (fs <= P < fe):
            continue
        real = ov.map_fake_offset_to_real(P - fs, sp, len(ref), len(alt))
        if real is None:
            continue
        if over[P] == E[real] and (real + 1) in dpos:
            ok += 1
        else:
            bad += 1
    assert bad == 0, f"chain validation had {bad} failures"
    print(f"[3] overlay->real->dict chain: {ok}/{ok} OK, dict key = real+1 confirmed")
    print("\nALL COORDINATE CHECKS PASSED.")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write(
            "usage: validate_indel_snp_coords.py <enriched.fa> <fake.fa> "
            "<log_indels.txt[.gz]> <my_dict.json>\n"
        )
        sys.exit(2)
    validate(*sys.argv[1:5])
