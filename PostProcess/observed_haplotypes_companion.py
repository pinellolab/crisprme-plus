"""ADDITIVE observed-haplotypes companion writer (CRISPRme+).

Surfaces, per OFF-TARGET SITE that carries >=2 distinct observed variant
haplotypes, EVERY such haplotype with its variant set, BOTH scores (CFD +
CRISPR-Bulge), its carriers, allele frequency, and phase. This is the per-sample
counterpart to the single per-cluster REPRESENTATIVE the ``integrated_results``
best-merge collapses: a site can be reached by several distinct genotype
combinations (haplotypes), each a separate experimental liability, but the
best-merge keeps only the highest-scoring one on each axis. The companion lists
them all so a reviewer sees the full set of observed-haplotype off-targets a site
contributes.

WHY A COMPANION FILE (not a new best-merge column): identical reasoning to the
other co-occurrence companions -- the finalized ``final_line`` addresses its
Reference / score-sentinel tail by NEGATIVE index in the CFD / CRISPR-Bulge
scorers, so appending a trailing column silently corrupts every variant row's
score, and the positional/name-list downstream consumers would desync or drop it.

DUAL-SCORE DESIGN: the report ships both a CFD-ranked and a CRISPR-Bulge-ranked
view, so every row carries BOTH scores. A haplotype that is the per-cluster CFD
representative supplies the CFD score; the same or a different haplotype that is
the per-cluster CRISPR-Bulge representative supplies the CRISPR-Bulge score; when
one axis never named a given haplotype its score reads "NA".

PHASE: ``integrated_results`` carries no phase column, so Phase is reported as
"NA" here (a follow-up may join a phase source; this writer never fabricates it).

GATE: the caller only records rows for a site with >=2 DISTINCT variant
haplotypes, and only writes the file when the row list is non-empty -- so a
single-haplotype or reference-only run writes NOTHING (byte-identical). Any error
is caught by the caller.

STDLIB ONLY.
"""

from __future__ import annotations

COMPANION_HEADER = (
    "#Chromosome\tSite_start\tStrand\tSpacer+PAM\tVariant_info_genome\t"
    "Aligned_protospacer+PAM\trsIDs\t"
    "Mismatches\tBulges\tBulge_type\tCFD\tCRISPR_Bulge\tOrigin\tAF\tPhase\t"
    "N_carriers\tCarriers"
)


def write_companion(out_path, rows):
    """Write the observed-haplotypes companion TSV.

    Args:
      out_path: destination path (``observed_haplotypes.tsv``).
      rows: list of dicts with keys Chromosome, Site_start, Strand, Spacer+PAM,
        Variant_info_genome, Aligned_protospacer+PAM, rsIDs, Mismatches, Bulges,
        Bulge_type, CFD, CRISPR_Bulge, Origin, AF, Phase, N_carriers, Carriers.

    Returns the number of data rows written.
    """
    with open(out_path, "w") as fh:
        fh.write(
            "# CRISPRme+ observed-haplotypes companion. One row per DISTINCT observed\n"
            "# variant haplotype at an off-target SITE that carries >=2 such haplotypes.\n"
            "# These are the per-sample haplotypes the per-cluster representative in\n"
            "# integrated_results collapses; each row carries BOTH scores so it serves\n"
            "# the CFD-ranked AND the CRISPR-Bulge-ranked report views.\n"
            "#   CFD / CRISPR_Bulge = score from whichever axis named this haplotype (NA if none).\n"
            "#   Origin             = always 'alt' (reference sites are not variant haplotypes).\n"
            "#   AF                 = variant allele frequency (NA if unknown).\n"
            "#   Phase              = NA (integrated_results carries no phase column).\n"
            "#   N_carriers/Carriers = observed carrier count + sample list (0/NA if none).\n"
        )
        fh.write(COMPANION_HEADER + "\n")
        n = 0
        for rec in rows:
            fh.write(
                "\t".join((
                    str(rec.get("Chromosome", ".")),
                    str(rec.get("Site_start", ".")),
                    str(rec.get("Strand", ".")),
                    str(rec.get("Spacer+PAM", ".")),
                    str(rec.get("Variant_info_genome", ".")),
                    str(rec.get("Aligned_protospacer+PAM", "NA")),
                    str(rec.get("rsIDs", "NA")),
                    str(rec.get("Mismatches", "NA")),
                    str(rec.get("Bulges", "NA")),
                    str(rec.get("Bulge_type", "NA")),
                    str(rec.get("CFD", "NA")),
                    str(rec.get("CRISPR_Bulge", "NA")),
                    str(rec.get("Origin", "alt")),
                    str(rec.get("AF", "NA")),
                    str(rec.get("Phase", "NA")),
                    str(rec.get("N_carriers", "0")),
                    str(rec.get("Carriers", "NA")),
                )) + "\n"
            )
            n += 1
    return n
