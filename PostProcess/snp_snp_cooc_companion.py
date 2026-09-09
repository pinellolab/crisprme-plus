"""ADDITIVE SNP+SNP co-occurrence companion writer (CRISPRme+ 2.5.2).

Surfaces, per emitted VARIANT off-target that USES >=2 co-occurring SNP alt alleles,
the SNP+SNP combination that makes the off-target real -- the SNP-side analogue of the
SNP+indel co-occurrence companion (``*.snp_indel_cooc.tsv``). One row per off-target,
keyed by the SAME identity columns the bestMerge uses (Chromosome, Position, Direction,
crRNA, DNA) so a downstream/web join is a straight key match.

WHY A COMPANION FILE (not a new bestMerge column): identical reasoning to
``phase_confirmation_companion`` -- the finalized ``final_line`` addresses its
Reference / ref-score-sentinel / tmp_pos_mms tail by NEGATIVE index in the CFD/CRISTA
scorers, so appending a trailing column silently corrupts every variant row's score,
and the positional/name-list downstream consumers would desync or drop it.

PHASE + FREQUENCY SEMANTICS:
  * Phase CONFIRMED -- every cis carrier reached this multi-SNP haplotype via a phased,
    same-phase-set path (genotyped path, observed enumerator). PUTATIVE -- >=1 carrier
    unphased / cross-phase-set, OR the registry-only (sites-only, e.g. mega) path where
    no per-sample genotypes exist to confirm cis.
  * MinAF_bound -- min over the participating SNPs' marginal allele frequencies. A cis
    haplotype can NEVER be more frequent than its rarest allele, so this is a valid
    conservative UPPER bound on the joint cis AF for BOTH phases (exact observed counts
    ride along in N_carriers / Carriers for the genotyped case). "." if unknown.
  * N_carriers / Carriers -- observed cis-carrier count + sample list (genotyped path);
    "NA" on the registry-only path (uncountable without genotypes).

GATE: the caller only records a row when an off-target uses >=2 distinct SNP positions,
and only writes the file when the row list is non-empty -- so on a single-SNP or legacy
run NOTHING is written (byte-identical). Any error is caught by the caller.

STDLIB ONLY.
"""

from __future__ import annotations

COMPANION_HEADER = (
    "#Chromosome\tPosition\tDirection\tcrRNA\tDNA\tSNP_positions\trsIDs\t"
    "Phase\tMinAF_bound\tN_carriers\tCarriers"
)


def write_companion(out_path, rows):
    """Write the SNP+SNP co-occurrence companion TSV.

    Args:
      out_path: destination path (``<outputFile>.snp_snp_cooc.tsv``).
      rows: list of dicts with keys Chromosome, Position, Direction, crRNA, DNA,
        SNP_positions, rsIDs, Phase, MinAF_bound, N_carriers, Carriers.

    Returns the number of data rows written.
    """
    with open(out_path, "w") as fh:
        fh.write(
            "# CRISPRme+ SNP+SNP co-occurrence companion. One row per variant off-target\n"
            "# that USES >=2 co-occurring SNP alt alleles (deduped by identity). Join FROM\n"
            "# bestMerge/bestCFD/bestCRISTA by (Chromosome, Position, Direction, crRNA, DNA).\n"
            "#   Phase CONFIRMED = phased same-phase-set cis haplotype (genotyped path).\n"
            "#   Phase PUTATIVE  = unphased / cross-phase-set / registry-only (sites-only).\n"
            "#   MinAF_bound     = min marginal AF = conservative upper bound on joint cis AF.\n"
            "#   N_carriers/Carriers = observed cis carriers (NA on the registry-only path).\n"
        )
        fh.write(COMPANION_HEADER + "\n")
        n = 0
        for rec in rows:
            fh.write(
                "\t".join((
                    str(rec.get("Chromosome", ".")),
                    str(rec.get("Position", ".")),
                    str(rec.get("Direction", ".")),
                    str(rec.get("crRNA", ".")),
                    str(rec.get("DNA", ".")),
                    str(rec.get("SNP_positions", ".")),
                    str(rec.get("rsIDs", ".")),
                    str(rec.get("Phase", ".")),
                    str(rec.get("MinAF_bound", ".")),
                    str(rec.get("N_carriers", "NA")),
                    str(rec.get("Carriers", "NA")),
                )) + "\n"
            )
            n += 1
    return n
