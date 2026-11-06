# SNP+indel oracle (synthetic, real CLI)

`run_oracle.sh <empty work dir> [<oracle dir>] [cas12a]` runs inside a CRISPRme image. It
generates a 16.6 kb reference chromosome and a phased VCF of 59 samples (`oracle.py
generate`), builds the indexes and dictionaries with `build-index-only`, runs
`complete-search` with `--per-sample` and with `--population-level`, then checks every row
whose variant list holds an indel or an MNV (`oracle.py check`) against each sample's two
haplotypes built from the VCF, and runs the ALT-PAM and no-N release gate
(`scripts/check_alt_pam.py`) on both results. It exits non-zero on any failure; CI runs it
in `validate-benchmarks.yml`, with both profiles: SpCas9 (3' PAM NRG, the sites below) and
`cas12a` (5' PAM TTTV, the C sites further below).

Every such row of the integrated results and of the alternative alignments, in all
three blocks, must show a concrete PAM of the searched motif; its CFD must be the CFD of its
own alignment (3' PAM); its mismatch count must be the alignment's and its mismatches +
bulges their sum; its start must be the reference coordinate of its target's first base (the
leftmost, forward strand) on a listed carrier's haplotype, at an occurrence that holds every
variant the row names, or, when that base is one an insertion adds, the insertion's anchor
(its VCF POS); its REF text must be the reference at its start, the bases its target holds,
gapped as its target; every variant it names must be a record of the VCF, forward strand;
every carrier it lists must have the sequence: a phased carrier on one of its haplotypes at
that site, spanning the indel, with every variant the row names on that haplotype; an
unphased one every allele the row names. Its `Variant_info_spacer+PAM` must number each
variant as a SNP row numbers its variants on the reference, along the target on that
haplotype (the distance from the target's first base to the variant's POS base, the first
being 1; from the right on the - strand, len(aligned spacer) + 1 - distance, less a
deletion's deleted bases; a variant before the first base is 1 on the + strand, the last on
the - strand). Its `PAM_creation` must be the ALT's PAM when the variants created the PAM, NA
otherwise: a PAM on reference bases is created when the reference bases it sits on lack the
motif, one made wholly of inserted bases always, one holding some inserted bases when the
reference PAM beside the protospacer (on the PAM's side of the protospacer's base that touches
the PAM) lacks it, or when that base is inserted too. The planted sites then pin who has which
target:

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
| I | a long insertion holding a whole perfect target | S44 CFD 1 |
| P | an insertion, and on S45's slot another insertion and a SNP at its POS and anchor | S45 (a phasing artefact) and S46 CFD 1 |
| K, K2 | an equal-length multi-base record, an MNV (+, -), and on both of S48's, S50's slots a SNP at the target's far-end base | S47, S49 CFD 1 alone; S48, S50 off the perfect row (an RNA-bulge alignment that does not reach that base lists them under the MNV alone); S48 on a one-mismatch row |
| W | a deletion, the PAM TGG followed by a G; on S51's slot a 1-bp deletion in that G run, on S52's a spacer deletion | S51 CFD 1 (the run puts the G back: the target's bases are unchanged); S52 on no row |
| V | a deletion, the target's 5' end in an AC repeat; on S53's slot a CA insertion written inside the target (left unnormalized) | S53 CFD 1 (the target is on its haplotype up to its last base) |
| Z, Z2 | an insertion holding the target's first bases: the 5' end of the spacer (+), the PAM (-) | S54, S55 CFD 1, the row at the insertion's anchor; Z2's PAM, wholly inserted bases, created (TGG) |
| Z3 | Z's geometry and a SNP alt the guide needs past the insertion, cis (S56) and trans (S57) | S56 CFD 1 under both variants, CONFIRMED, at the anchor, the SNP the target's 13th base; S57 off the perfect row |
| Z4, Z5 | (-) the PAM's first two bases (forward CC) the last two an insertion adds, its third and the protospacer reference bases; beside the protospacer the reference holds CTA (TAG on -) at Z4, ATA at Z5 | S58, S59 CFD 1 at the anchor; PAM_creation NA at Z4 (the reference has the target with TAG), TGG at Z5 |

The `cas12a` profile (a 3 kb reference, 8 samples, 5' PAM TTTV, a 20-nt spacer; CFD is not
computed for a 5' PAM, so a perfect row is one with no mismatch and no bulge):

| Site | Planted | Expected |
|---|---|---|
| C1 | (+) an insertion holding the PAM and the spacer's first three bases; a SNP the spacer needs past it, cis (S01) and trans (S02) | S01 perfect at the anchor, the SNP the target's 15th base, the PAM created (TTTA); S02 one mismatch |
| C2, C2B | (+) the PAM's first two bases the last two an insertion adds; the reference before the spacer TTTA at C2, GTTA at C2B | S03 (per-sample: the reference holds the same target one base left, which represents the site in a population-level run) PAM_creation NA; S04 TTTA |
| C3 | (-) an insertion holding the spacer's 3' end (the target's leftmost bases); a SNP past it, cis (S05) and trans (S06) | S05 perfect at the anchor, the SNP the 16th base, the PAM (reference bases) not created; S06 one mismatch |
| C4 | (-) the PAM, at the right, wholly bases an insertion adds | S07 perfect, the PAM created (TTTA) |
| C5 | (+) a deletion in the spacer and a SNP past it, cis (S08) | S08 perfect, the SNP the 20th base |

Up to 2.8.0 the indel post-analysis read each IUPAC code of the SNP-overlaid fake-indel
genome column by column, without genotypes: on aa9c54b (the 2.8.1 release candidate)
the check fails at every site but L2 (at C, I, J and P only the REF check fails). On
3520eeb, before the fix of an RNA-bulge row's start, the REF check fails on 71 row
blocks: every RNA-bulge row on the + strand, which started one base early, and every
RNA,DNA row, whose REF text had no gap and one base too many; and S50 is listed on K2's
perfect row, in both modes. On d4730a0, which has that fix but not the merge keeping an
MNV row's carriers on its sequence, it fails at K and K2 in both modes: S48 and S50 are
listed on the perfect row (K's RNA-bulge row, no longer one base early, now meets the
site's other alignments). On 9c6b4046, which took a carrier off an indel row whenever its
slot held a second indel inside the target, decided on the indel's position, it fails at W
and V in both modes: S51 and S53, whose haplotypes hold the target, are on no CFD-1 row.
S51 is on no row at all, and S53 only on the CFD-0.9 RNA-bulge alignment that starts past
the insertion. On the released v2.8.1 image the check fails 74 times: the start at I, Z, Z2,
Z3, Z4 and Z5 (a row whose first base is one an insertion adds started at that base's place
on the haplotype, 3 to 6 bases past the anchor here, and its REF text was read there),
`Variant_info_spacer+PAM` at D, M2, R and Z3 (variants numbered along the reference: D's SNP
at the PAM's last base, past a deletion, read 24 for 23, R's 5 for 6, M2's 3 for 1, Z3's 4 for
13), and `PAM_creation` at Z5 (NA, read from the REF text at the haplotype place); and the
`cas12a` profile 47 times, at C1, C2, C2B, C3 and C5. With only the first 2.8.2 commit
(d7a6fff, rows at the anchor) it fails 51 and 24 times: `Variant_info_spacer+PAM` as above
(Z3: 8), and `PAM_creation` at I, Z2 (NA) and Z4 (TGG), read from the REF text at the anchor,
and C2. In the genome-wide HGDP + 1000 Genomes search of the README guide (per-sample, 4
mismatches + 1 + 1) 229 highest-CFD rows of 2.8.1 start inside an insertion, up to 229 bases
past the anchor; on 2.8.0, an alternative alignment of that guide started at chr22:22322768,
the 7th base of G>GCAATCATTAT at 22322762.

Not covered: split multi-allelic SNP records at one position. CRISPRitz add-variants
keeps only the last of them in the enriched genome and in the dictionary, so the other
alt is invisible to every stage; `PostProcess/test_indel_haplotype_resolution.py` covers
a position whose genotype source knows both alts. (Several alts on one VCF line reach
every stage: sites N and N2.) Nor an MNV carrier heterozygous for a SNP its target
covers: the phased indel genotype store holds only records whose REF and ALT differ in
length, so an MNV's phase is unknown to the post-analysis, and such a carrier gets the
strongest target its own alleles allow, PUTATIVE, as an unphased carrier does (site U).
