<!-- car track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the docs track: three amendments to entries already in `## [2.8.1]`, then new entries.**
Two of this track's fixes narrow behavior that is itself new in 2.8.1, which 2.8.0 never had:
the second-indel rule of the indel post-analysis, and the guide length read from an index's
`.bin` header. An entry of its own would describe a failure no release had, so each is an
amendment. The third takes the `PAM_creation` count out of the RNA-bulge entry: the new
`PAM_creation` entry below changes how that column is set. Replace the quoted text (line numbers
are those of `CHANGELOG.md` at 9c6b4046) with the text after it. All three entries are under
`### Fixed`.

1. **"An individual whose haplotype carries a second indel inside an indel off-target was listed
   for the sequence of the first one"**, lines 2460-2476: replace the whole entry with

   > - **An individual whose haplotype carries a second indel that changes an indel off-target was
   >   listed for the sequence of the first indel.** The fake-indel genome holds one indel per
   >   window, so such a haplotype does not have the window's sequence; on the synthetic oracle a
   >   sample carrying two deletions on one haplotype was listed on both windows' rows, for
   >   sequences neither of its haplotypes has. Whether a carrier whose haplotype holds a second
   >   indel inside the target has the row's sequence is read on its bases: the reference around the
   >   target is written with every indel of the phased indel store the haplotype carries there and
   >   its own SNP alleles, and the carrier is listed when the target's bases are on it at the
   >   target's place, from the target's first base on or up to its last, through the indel (not in
   >   a stretch that another of its variants lengthened, pushing the indel out of it). In a tandem
   >   repeat a unit gained or lost at the target's edge leaves the bases as they are (the repeat
   >   puts the same bases back), so such a carrier stays listed; one whose second indel changes
   >   them does not, nor one whose bases cannot be read either way: each way meets an unphased or
   >   missing SNP call, an unphased indel that does not itself change the target, or two of its
   >   variants that overlap, two records of one base: two other indels
   >   whose REFs share a base (any two at one position do, even a duplication and the deletion of
   >   the same bases, which together would leave the reference), or a SNP on a base of the second
   >   indel's REF. A reading from the first base on meets no variant upstream of it, one up to the
   >   last base none downstream. A carrier whose second indel inside the target has an unphased
   >   call cannot be tied to a haplotype and stays listed, PUTATIVE, as every carrier on an
   >   unphased panel is, even when another, phased second indel on its indel haplotype changes
   >   the target. Nor is one listed whose second indel, on the same phased slot,
   >   overlaps the first at another position (one complex allele that the VCF splits into
   >   overlapping records). A variant only on bases the indel deletes changes nothing. Another
   >   allele at the indel's own position, an indel or a SNP, is an alternative to it (a haplotype
   >   holding both would be one record), so a slot phased with both is a phasing artefact: on chr22
   >   of the HGDP + 1000 Genomes panel, two indels at one position share a slot for 40,679 carrier
   >   pairs against 2,463,324 on opposite slots, a SNP and an indel for 2,404 against 167,951. And
   >   an unphased carrier of two overlapping records is a multi-allelic heterozygote, one allele
   >   per haplotype. All of these stay listed. This needs the phased indel genotype store, which
   >   the build writes with the SNP overlay; a row left with no carrier is dropped. Its own target,
   >   which needs both indels, is still not searched (METHODS §8). On the search above, 693 targets
   >   with no SNP code lost such carriers. On chr22 of that panel (per sample, five guides),
   >   checked against haplotypes built from the phased VCF: of the 7,046 carrier listings the rule
   >   takes off the integrated results' highest-CFD block (96 rows are left with none, and
   >   dropped), 7,043 are of a haplotype that lacks the row's sequence within 40 bases of it, 2 of
   >   one that has it a few bases away without the indel (a repeat's copy, a target of its own) and
   >   1 of a slot with two overlapping records; no carrier left on an indel row lacks the row's
   >   sequence because of a second indel.

2. **"An indel off-target with an RNA bulge was reported one base early, and the reference sequence
   shown for one with both bulge types was a base off"**, lines 2514-2515. Replace

   > `PAM_creation` is now set on 815 and no longer set on 240 (all on the + strand, where the old
   > REF's PAM was read one base off), `CFD_score_REF` changed on 785, and `CFD_score` went down on 218

   with

   > `CFD_score_REF` changed on 785 (`PAM_creation` no longer reads this REF: the entry "An indel
   > row was flagged as creating a PAM the reference already has"), and `CFD_score` went down on 218

3. **"A guide that does not fit an index's stored window is no longer searched on it"**, lines
   2142-2143. Replace

   > The length comes from the index's own
   > `.bin` header, so no index needs rebuilding.

   with

   > The length comes from the index's own `.bin` header, so no index needs rebuilding. An index a
   > newer CRISPRitz built, in the format CRISPRitz 2.9.2 and later write (index format 2, whose
   > files open with a magic and the format version), holds no length this release reads: a search
   > that takes one, found or built with such a CRISPRitz installed, stops with one ERROR that says
   > what it is and what to do (use a CRISPRme that reads index format 2, or install CRISPRitz 2.9.1
   > and build the index again), and so does `build-index-only` with one in place that it would
   > keep (the reference, the variant or the INDELS index), or one the installed CRISPRitz writes.
   > Measured through the real CLI with that CRISPRitz built into the image, they had read its
   > format version as the guide length ("built for a 2-nt guide") and its magic as a leaf count
   > ("not a tree ..., as when the disk fills up"); with CRISPRitz 2.9.1 the same runs give the
   > same tables. Only a Bioconda install can pull such a CRISPRitz: the image and
   > `install_from_source.sh` pin 2.9.1.

### Fixed

- **An indel row was flagged as creating a PAM the reference already has.** The reference bases the
  `PAM_creation` test read for an indel row were its REF alignment's, read from the row's start; on
  the + strand with a 3' PAM (on the - strand with a 5' one, such as Cas12a's) a deletion or an
  insertion in the protospacer moves the PAM `len(REF) - len(ALT)` bases from there, so the test
  read other bases, and a PAM the reference already has was flagged as created: on the WALK panel
  the rows chr2:10000 (`rsW_J_del`) and chr3:20000 (`rsW_N_del` + `rsW_N14`) showed `PAM_creation`
  TGG over a reference TGG, counted in the web's Summary by Mismatches/Bulges and Summary by Sample
  (S2, S4, S6) and on the Personal Risk Card (S4 "PAM creation 2"). The test now reads the reference
  bases the PAM sits on, through the indel; a PAM that holds an ALT base or spans a deletion's
  junction keeps the REF alignment's bases, as before, and so does a PAM on the other side, which no
  indel inside the target moves. **Affected:** indel rows of every variant search since
  `PAM_creation` was first set on them (2.7.10, #94), in 2.8.0 too; only that column and the counts
  built from it change (the report's `PAM_creation`, the summaries, the Personal Risk Card, the
  web), not the detection, scores, carriers or ranking. Measured genome-wide on the default HGDP +
  1000 Genomes index (per sample, 4 mismatches + 1 DNA + 1 RNA bulge), on the highest-CFD block: for
  guide `ACTGAAATCTGTAAGCAGGC` (NRG) the flag is no longer set on 1,954 indel rows and is set on 1
  more; rebuilt from the reference through each row's variants, every one of the 1,834 of them that
  could be checked lies over a reference PAM and the one more over none; of the rows still
  flagged that could be checked, 11 hold a PAM the reference lacks and 120 a PAM that holds a
  base of the indel, flagged as before, and the other 55 (186 in all) have their PAM across the
  indel's junction or inside a long insertion, also flagged as before. For the README's guide, `CTAACAGTTGCTTTTATCAC`, 2,655 no longer set (2,561
  checked, all over a reference PAM; of the rows left that could be checked, 19 + 150), which on its report changes that column, and no
  other, on 15 rows of the top-1000 table and 2 of the recommended panel; for TRAC,
  `CTCTCAGCTGGTACACGGCA`, 1,124 no longer set and 1 more (1,047 checked; of the rows left that could be checked, 10 + 89); for
  `ACTGAAATCTGTAAGCAGGC` on the NRR index, 2,761 and 2 more (2,585 checked; of the rows left that could be checked, 15 + 196). **What
  to do:** run such searches again. Covered by `PostProcess/test_indel_pam_creation.py` (the WALK J
  and N rows among them) and `PostProcess/test_indel_row_coordinates.py`.

### Changed

- **A known limit, stated with its scale: a SNP row lists a carrier whose haplotype holding its SNPs
  also carries an indel inside the target.** A SNP row's carriers are resolved from SNP genotypes
  only, so an individual whose haplotype holding the row's SNPs also carries an indel inside the
  target is still listed, though its own sequence there is the indel's (reported on the indel row
  when it is in budget). On chr22 of the HGDP + 1000 Genomes panel (per sample, NRG, 4 mismatches +
  1 DNA + 1 RNA bulge, five guides) this affects 2.2% of SNP row-blocks (990 of 45,921, 194 in the
  integrated results), and 545 of them have no listed carrier that holds the sequence; with no
  bulges, 12 row-blocks. It only over-lists, no site is lost, and 2.8.0 behaves the same. METHODS
  sections 4 and 8 now say so; it is fixed in the next release.
