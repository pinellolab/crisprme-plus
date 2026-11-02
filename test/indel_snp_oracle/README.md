# SNP+indel oracle (synthetic, real CLI)

`run_oracle.sh <empty work dir>` runs inside a CRISPRme image. It generates a 16 kb
reference chromosome and a phased VCF of 46 samples (`oracle.py generate`), builds the
indexes and dictionaries with `build-index-only`, runs `complete-search` with
`--per-sample` and with `--population-level`, then checks every row whose variant list
holds an indel (`oracle.py check`) against each sample's two haplotypes built from the
VCF, and runs the ALT-PAM release gate (`scripts/check_alt_pam.py`) on both results.
It exits non-zero on any failure; CI runs it in `validate-benchmarks.yml`.

Every indel row of the integrated results and of the alternative alignments, in all
three blocks, must show a concrete NRG PAM; its CFD must be the CFD of its own
alignment; its mismatch count must be the alignment's and its mismatches + bulges their
sum; every variant it names must be a record of the VCF, forward strand; and every
carrier it lists must have the sequence: a phased carrier on one of its haplotypes at
that site, spanning the indel, with every variant the row names on that haplotype; an
unphased one every allele the row names. The planted sites then pin who has which target:

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
| M, M2 | deletion (+), insertion (-), and on S24's, S26's slot a SNP at the target's far-end base | S23, S25 CFD 1 alone; S24, S26 on a one-mismatch row, off the perfect one (an RNA-bulge alignment that does not reach that base lists them under the same variants) |
| N, N2 | deletion + C>A,T on one VCF line, A the guide's base (2\|1, 2\|2, 0\|1) | the deletion slot's own allele: S28, S31 CFD 1; S27, S29, S32 T, one mismatch; S30 C |
| L, L2, L3 | a target in a deletion's downstream flank, one in an insertion's, one ending at an insertion's anchor | no indel row; L is the SNP off-target of its carriers, S33, S34, S35 |
| R | F's geometry on the minus strand | S37 CFD 1, the SNP named forward (G>T); S38 one mismatch |
| Q | a target with a SNP code and, on S39's slot, a second deletion inside it | S40 listed; S39 on no row |
| H | an unphased carrier homozygous for a SNP that spoils the target | S41 on the mismatch row, never the perfect one; S42 CFD 1 |
| J | a carrier homozygous for the deletion and for the SNP the target needs | S43 CFD 1; joint frequency 2 / AN |

Up to 2.8.0 the indel post-analysis read each IUPAC code of the SNP-overlaid fake-indel
genome column by column, without genotypes: on aa9c54b (the 2.8.1 release candidate)
the check fails at every site but the control C and the sites with no SNP (L2, L3 hold
no SNP but L3's row is there all the same).

Not covered: split multi-allelic SNP records at one position. CRISPRitz add-variants
keeps only the last of them in the enriched genome and in the dictionary, so the other
alt is invisible to every stage; `PostProcess/test_indel_haplotype_resolution.py` covers
a position whose genotype source knows both alts. (Several alts on one VCF line reach
every stage: sites N and N2.)
