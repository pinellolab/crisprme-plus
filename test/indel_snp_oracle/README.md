# SNP+indel oracle (synthetic, real CLI)

`run_oracle.sh <empty work dir> [<oracle dir>] [cas12a]` runs inside a CRISPRme image. It
generates a 22 kb reference chromosome and a phased VCF of 89 samples (`oracle.py
generate`), builds the indexes and dictionaries with `build-index-only`, runs
`complete-search` with `--per-sample` and with `--population-level`, then checks every row
whose variant list holds an indel or an MNV (`oracle.py check`) against each sample's two
haplotypes built from the VCF, and the carriers of every per-sample SNP-only row ("SNP rows"
below), and runs the ALT-PAM and no-N release gate
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
unphased one every allele the row names, and its target must be the strongest its own alleles
allow on the row's alignment (PUTATIVE: the guide's base where the alignment pairs the position
with a guide base, else the alleles with the higher CFD; with a 5' PAM, fewer mismatches). Its
`Variant_info_spacer+PAM` must number each
variant as a SNP row numbers its variants on the reference, along the target on that
haplotype (the distance from the target's first base to the variant's POS base, the first
being 1; from the right on the - strand, len(aligned spacer) + 1 - distance, less a
deletion's deleted bases; a variant before the first base is 1 on the + strand, the last on
the - strand). Its `PAM_creation` must be the ALT's PAM when the variants created the PAM, NA
otherwise: never created when its REF alignment is its ALT (the reference holds the target at
its start); otherwise read on the haplotype the row's variants make (the reference with its
SNP alts and its indel), over every placement of the indel that keeps that haplotype, each
built whole (an insertion's bases rotated through its repeat, a deletion's junction moved
through it): created when no placement leaves the PAM on reference bases holding the motif,
and, where none leaves it on reference bases at all, when the reference PAM beside the
protospacer's base that touches the PAM lacks the motif in a placement that leaves that base
on the reference, or when no placement does (an MNV: the reference bases under the PAM; a
complex record, bases replaced and bases added or lost: its first min(len(REF), len(ALT))
changed bases in the places of the REF bases they replace, read as those, and its other bases
an insertion's or a deletion's). The planted sites then pin who has which target:

| Site | Planted | Expected |
|---|---|---|
| A | deletion + G>A at the PAM's last base, in cis (S02) and in trans (S03) | S01, S03 CFD 1 (TGG); S02 no target |
| T, T2 | (+) an insertion after a target's reference protospacer adds TGG: a copy of the reference's own TGG there (T, a tandem duplication), or before the reference's AGG (T2) | per-sample (the reference's own target represents each site in a population-level run): S62, S63 CFD 1, the row at the target's first base, a reference one; PAM_creation NA (at T a placement 3 bases right leaves the PAM on the reference's own TGG; at T2 the PAM is wholly inserted wherever the insertion sits, and the reference holds AGG beside the protospacer) |
| B | minus-strand insertion + PAM SNP in cis (S06) and in trans (S07) | S05, S07 CFD 1 (AGG) |
| C | deletion, no SNP | S09 CFD 1 |
| D | deletion + a PAM-creating A>G SNP in cis (S10) | S10 CFD 1, CONFIRMED co-occurrence; S11 no target |
| E, E2 | a SNP the guide needs at every base, no haplotype carrying them together | no indel row above CFD 0.5 |
| F | deletion + a spacer SNP the guide needs, cis (S17) and trans (S18) | S17 CFD 1, CONFIRMED; S18 one mismatch |
| U | F's geometry, unphased | S19 CFD 1, PUTATIVE, and on each alignment it is listed on, the strongest target its alleles allow there (C where the SNP's base is a DNA bulge, the guide's A elsewhere); S20 one mismatch |
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
| X, X2 | (+) an insertion holding two whole targets: the guide with TGG and, 25 bases on, with AGG (X), or with TGG again (X2, a repeat) | per-sample: S60 listed with both targets, both rows at the anchor, one the site's row and the other one of its alternative alignments; S61 on one row with the target (its two copies are the same row at the anchor) |
| PR, PM | (+, -) a tandem duplication of a 27-base unit holding a whole target, the copy written left of the reference's: the copy's target starts inside the insertion, its PAM wholly inserted | per-sample (the reference's target 3 bases on represents each site in a population-level run): S79, S80 CFD 1 at the anchor; PAM_creation NA (placed one unit right, the insertion leaves the PAM on the reference's own) |
| QC | (-) C>CCAG after the anchor C, the PAM's far end; the reference repeats its first base: in both places it holds the PAM's inner bases and the protospacer's first forward base | S81 CFD 1; PAM_creation TGG (no placement gives the PAM a place on the reference; the REF alignment's window, CCC from the anchor, is GGG on -) |
| QN | (-) C>CCA after the anchor C, no repeat; the reference holds CTC (GAG on -) beside the protospacer | per-sample (the reference's own target 2 bases left represents the site in a population-level run): S82 CFD 1; PAM_creation NA (the window, CGC, is GCG on -) |
| QD | (-) CA>C deletes the A of CCAG before the protospacer, its junction inside the PAM CCG (CGG on -) | S83 CFD 1; PAM_creation CGG (the reference holds CAG, CTG on -, beside the protospacer; the window from the row's start, CCA, is TGG on -, another protospacer's PAM) |
| QW | (-) TTCCA inserted after the anchor G: a target from its 3rd added base, its PAM the last three, and the reference's CCG (CGG on -) before the protospacer | per-sample (the reference's own target represents the site in a population-level run): S84 CFD 1 at the anchor; PAM_creation NA |
| QA | (+) G[16:] + TGG inserted after G[15], no repeat | S85 CFD 1; PAM_creation TGG (the protospacer's end and the PAM are inserted wherever the insertion sits; the window holds AGG) |
| QX | (+) a complex record, A>T and 33 bases added after it that hold a whole target from the 6th, its PAM and the protospacer's end added bases | S86 CFD 1, the row at the replaced A; PAM_creation TGG (as for the insertion of those bases after an anchor; the window, read from the replaced A, holds AGG) |
| QY | (-) its mirror, C>T and the bases added, the PAM (forward CCA) and the protospacer's first forward base added | S87 CFD 1; PAM_creation TGG (the window, CCG from the replaced C, is CGG on -) |
| QZ | (-) a complex record, G>CCA: the row starts on the replaced G, its PAM CCA (TGG on -) that base and the two added ones, the protospacer reference bases, and the reference holds CTG (CAG on -) beside it | S88 CFD 1, the site's row in both modes; PAM_creation NA (the window from the replaced G, GGC, is GCC on -) |
| QV | (-) a complex record, TG>C: the T replaced by C and the G lost, the junction inside the PAM CCA (TGG on -) | S89 CFD 1; PAM_creation TGG (the reference holds TGA, TCA on -, beside the protospacer; the window from the row's start, CTG, is CAG on -, with the lost G) |

The `cas12a` profile (a 3.8 kb reference, 11 samples, 5' PAM TTTV, a 20-nt spacer; CFD is not
computed for a 5' PAM, so a perfect row is one with no mismatch and no bulge):

| Site | Planted | Expected |
|---|---|---|
| C1 | (+) an insertion holding the PAM and the spacer's first three bases; a SNP the spacer needs past it, cis (S01) and trans (S02) | S01 perfect at the anchor, the SNP the target's 15th base, the PAM created (TTTA); S02 one mismatch |
| C2, C2B | (+) the PAM's first two bases the last two an insertion adds; the reference before the spacer TTTA at C2, GTTA at C2B | S03 (per-sample: the reference holds the same target one base left, which represents the site in a population-level run) PAM_creation NA; S04 TTTA |
| C3 | (-) an insertion holding the spacer's 3' end (the target's leftmost bases); a SNP past it, cis (S05) and trans (S06) | S05 perfect at the anchor, the SNP the 16th base, the PAM (reference bases) not created; S06 one mismatch |
| C4 | (-) the PAM, at the right, wholly bases an insertion adds | S07 perfect, the PAM created (TTTA) |
| C5 | (+) a deletion in the spacer and a SNP past it, cis (S08) | S08 perfect, the SNP the 20th base |
| C6 | (-) C4's geometry where the insertion copies the reference's own TAAA (TTTA on -) after the spacer, a tandem duplication | S09 perfect (per-sample: the reference's own target represents the site in a population-level run), the PAM not created |
| CB1 | (+) T>TTT written at the end of a T run (left unnormalized): the PAM the anchor, the two added T and a reference A | S10 perfect (per-sample: the reference's own target represents the site), the PAM not created (placed three bases left, as bcftools norm writes it, the insertion leaves TTTA on the reference; the REF alignment's window, TAGA from the anchor, held none) |
| CB3 | (-) the spacer's last four forward bases and the PAM (TAAA, TTTA on -) inserted after its 16th, no repeat | S11 perfect, the PAM created (TTTA; the window from the row's start holds TAAA, TTTA on -) |

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
the insertion. On the released v2.8.1 image the check fails 104 times: the start at I, Z,
Z2, Z3, Z4, Z5, X and X2 (a row whose first base is one an insertion adds started at that
base's place on the haplotype, 3 to 28 bases past the anchor here, and its REF text was read
there), `Variant_info_spacer+PAM` at D, M2, R and Z3 (variants numbered along the reference:
D's SNP at the PAM's last base, past a deletion, read 24 for 23, R's 5 for 6, M2's 3 for 1,
Z3's 4 for 13), `PAM_creation` at Z5 (NA, read from the REF text at the haplotype place),
and X2's target on two rows; and the `cas12a` profile 47 times, at C1, C2, C2B, C3 and C5.
With only the first 2.10.0 commit (d7a6fff, rows at the anchor) it fails 60 and 24 times:
`Variant_info_spacer+PAM` as above (Z3: 8), `PAM_creation` at I, Z2 and X (NA) and Z4 (TGG),
read from the REF text at the anchor, and C2; and at X the TGG target is on no row: it
shares the anchor, the variant and the carrier with the AGG one, and the merge folded the
two into one entry. On 43809d3, with the other 2.10.0 fixes but the merge still comparing a
site's alignments at their cluster position on the reference, that fails 3 times. There and
on c3a75fe, whose inserted-PAM rules judged every row, also one that starts on a reference
base, the PAMs of T and T2, and the `cas12a` profile's C6, read as created: 11 and 4
failures on 43809d3, 8 and 4 on c3a75fe. The released image and d7a6fff read those PAMs on
the REF alignment: NA, as now. In the genome-wide HGDP + 1000 Genomes search of the README
guide (per-sample, 4 mismatches + 1 + 1) 229 highest-CFD rows of 2.8.1 start inside an
insertion, up to 229 bases past the anchor; on 2.8.0, an alternative alignment of that guide
started at chr22:22322768, the 7th base of G>GCAATCATTAT at 22322762. From 75ff25b to 8d3f278
the check did not test the carriers a row lists, only that the row starts where a listed
carrier's occurrence does (75ff25b's rewrite of the start check dropped the test of every
carrier, and of every variant being a record of the VCF); it does again, and tests an unphased
carrier's target against the strongest its alleles allow on the row's alignment too. Up to
8d3f278 the merge keyed an indel row by its target alone, so two alignments of one target were
one entry, and that test fails twice, at site U, on every version before the fix (not counted
above): S19, listed with C on the alignment that makes the SNP's base a DNA bulge (CFD 0.014;
0 with A), went onto A-CTGAAATCTG…, which aligns that base with the guide's A, with C (CFD
0.423, where its alleles give 0.93), in the highest-CFD block, and onto the alignment at 7999
in the CRISPR-Bulge block.

Not covered: split multi-allelic SNP records at one position. CRISPRitz add-variants
keeps only the last of them in the enriched genome and in the dictionary, so the other
alt is invisible to every stage; `PostProcess/test_indel_haplotype_resolution.py` covers
a position whose genotype source knows both alts. (Several alts on one VCF line reach
every stage: sites N and N2.) Nor an MNV carrier heterozygous for a SNP its target
covers: the phased indel genotype store holds only records whose REF and ALT differ in
length, so an MNV's phase is unknown to the post-analysis, and such a carrier gets the
strongest target its own alleles allow, PUTATIVE, as an unphased carrier does (site U).

## SNP rows (2.10.0, `snp_rows.py`)

A SNP row's carriers came from SNP genotypes alone, so up to 2.8.1 it listed a carrier whose
haplotype holding its SNPs also holds an indel inside the target, though that haplotype does
not have the row's sequence. In the per-sample results, every SNP-only row (no indel or MNV
named) is now checked as the indel rows are: a phased carrier it lists must have its sequence
on a haplotype holding every SNP it names, the target's first base at the row's start, and an
unphased one must carry every SNP it names. (A row naming two alts of one position, which no
haplotype holds, is counted and not checked: site N2's lone `2|1` carrier of an alt on a
multi-allelic VCF line is read by the SNP path as carrying the other alt too, a defect of its
own.) Sites Y1-Y4 sit at 17000-17900 and use samples S70-S78 (the reference is 22 kb, the
panel 89 samples):

| Site | Planted | Expected (per sample) |
|---|---|---|
| Y1 | a SNP target (a C>A SNP makes it perfect) and a 2-bp insertion inside it: in cis (S70), in trans (S71), the insertion homozygous (S72), the SNP alone (S73), both unphased (S74) | S71, S73, S74 on the CFD-1 SNP row; S70, S72 on no SNP row |
| Y2 | Y1's geometry, the SNP phased and the insertion unphased (S75) | S75 on the CFD-1 SNP row, PUTATIVE |
| Y3 | Y1's SNP, the PAM TGG followed by a G, and on S76's SNP slot a 1-bp deletion in that G run | S76 on the CFD-1 SNP row (the run puts the G back) |
| Y4 | a minus-strand SNP target whose SNP is in cis with a 2-bp insertion inside it on every haplotype (S77, S78) | no SNP row lists them |

Without this fix (fix/282-indel-alignment at e4dbb9c) the check fails 62 times: the SNP-only
alignments of sites D, F, M, N, Q, R, J, Z3, Y1, Y3 and Y4 list the indel's cis carriers (57
row blocks), S70, S72, S77 and S78 are on a SNP row, and Y2's row is CONFIRMED; and the
`cas12a` profile 6 times: the SNP-only rows of C3 and C5 list S05 and S08, whose haplotypes
holding the SNP hold the insertion or the deletion. These failures come on top of the
indel-row failures counted above: the released v2.8.1 image fails those checks 166 times
(104 + 62) and the `cas12a` profile 53 times (47 + 6), and the whole oracle, with site U's two
failures of the check of an unphased carrier's target, 168 and 53 times (e4dbb9c 64 and 6,
8d3f278 2 and 0).

2.10.0's placement rule (a PAM holding inserted bases, or a deletion's junction, judged
wherever in its repeat the VCF writes the indel) added sites PR, PM, QC, QN, QD, QW and QA,
and CB1 and CB3 in the `cas12a` profile (samples S79-S85 and S10-S11; the reference grew to
21 kb and 3.8 kb). On 37d81de, before that rule, the oracle fails 42 and 10 times, all at those
sites (rows read where the VCF writes the indel: on their REF alignment's window from the row's
start, or, for a row starting inside an insertion, by the inserted-PAM rules there); on the
released v2.8.1 image 204 and 60 times.

The complex records' reading (a complex record's replaced bases read as the REF bases they
replace, its added or lost bases placed as an insertion's or a deletion's) added sites QX, QY,
QZ and QV (samples S86-S89; the reference grew to 22 kb). On 50803b6, which left a complex
record to its REF alignment's window, the SpCas9 oracle fails 31 times, all at those sites; on
37d81de, whose inserted-PAM rules already read QX and QY right (their rows start inside the
added bases), 57 times (the 42 above and 15 at QZ and QV); on the released v2.8.1 image 231
times. The `cas12a` profile is unchanged: 0 failures on 50803b6, 10 on 37d81de, 60 on v2.8.1.
