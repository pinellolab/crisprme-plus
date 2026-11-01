# SNP+indel oracle (synthetic, real CLI)

`run_oracle.sh <empty work dir>` runs inside a CRISPRme image. It generates a 10.5 kb
reference chromosome and a phased VCF of 22 samples (`oracle.py generate`), builds the
indexes and dictionaries with `build-index-only`, runs `complete-search` with
`--per-sample` and with `--population-level`, then checks every row whose variant list
holds an indel (`oracle.py check`) against each sample's two haplotypes built from the
VCF, and runs the ALT-PAM release gate (`scripts/check_alt_pam.py`) on both results.
It exits non-zero on any failure; CI runs it in `validate-benchmarks.yml`.

Every indel row of the integrated results and of the alternative alignments, in all
three blocks, must show a concrete NRG PAM; its CFD must be the CFD of its own
alignment; its mismatch count must be the alignment's; and every carrier it lists must
have the sequence: a phased carrier on one of its haplotypes at that site, an unphased
one every allele the row names. The planted sites then pin who has which target:

| Site | Planted | Expected |
|---|---|---|
| A | deletion + G>A at the PAM's last base, in cis (S02) and in trans (S03) | S01, S03 CFD 1 (TGG); S02 no target |
| B | minus-strand insertion + PAM SNP in cis (S06) and in trans (S07) | S05, S07 CFD 1 (AGG) |
| C | deletion, no SNP | S09 CFD 1 |
| D | deletion + a PAM-creating A>G SNP in cis (S10) | S10 CFD 1, CONFIRMED co-occurrence; S11 no target |
| E, E2 | a SNP the guide needs at every base, no haplotype carrying them together | no indel row above CFD 0.5 |
| F | deletion + a spacer SNP the guide needs, cis (S17) and trans (S18) | S17 CFD 1, CONFIRMED; S18 one mismatch |
| U | F's geometry, unphased | S19 CFD 1, PUTATIVE; S20 one mismatch |
| O | two deletions inside one target on S21's haplotype | S22 listed; S21 on no row |

Up to 2.8.0 the indel post-analysis read each IUPAC code of the SNP-overlaid fake-indel
genome column by column, without genotypes: on aa9c54b (the 2.8.1 release candidate)
the check fails at every site but the control C (121 failures).

Not covered: split multi-allelic SNP records at one position. CRISPRitz add-variants
keeps only the last of them in the enriched genome and in the dictionary, so the other
alt is invisible to every stage; `PostProcess/test_indel_haplotype_resolution.py` covers
a position whose genotype source knows both alts.
