<!-- alleles track: [2.9.0] entries for CHANGELOG.md, to be merged with the other tracks'. -->

**For whoever merges.** These entries describe CRISPRme's side of CRISPRitz 2.10.0's every-allele
enrichment and its fifth symbol Q; the engine's own entries (every ALT of a split site coded, Q
written at a position of four alleles, the records of one REF>ALT merged in the dictionary) are the
engine track's. They were checked through the real CLI with CRISPRitz 2.9.1, with the every-allele
enricher before Q (`feat/multiallelic-enrich` d75c603) and with 2.10.0's (`feat/2.10.0` 25ee50c)
beside 2.9.1's `buildTST` and `searchTST`; a search of an index holding Q needs CRISPRme to read
index format 2, which 2.10.0 writes (the format track's): with that track's 087fe9b and fa6dbb7
merged in, on 2.10.0's engine (25ee50c's `add-variants`, `buildTST` and `searchTST`), MULTIALLELIC
section 6.6's kit passed, the review's synthetic set gave the same tables on the tiers alone as on
tiers and dictionary, and so did an HPRC kit at a position of four bases, one of which nobody
carries (under Added); with CRISPRitz 2.9.1 (the 2.8.1 image), a search of the published mega index
on chr22 gives 2.8.1's highest CFD at every site (under Fixed). CRISPRitz 2.10.0 still writes a
record's first ID alone into the indel log, for every indel record of one ALT or several, and into
the dictionary for each ALT of a VCF line listing several (`enricher.py` and `enricher.cpp`), which
no CRISPRme reader can undo: the known limit under Changed, METHODS section 8 and the validator's
warning say so. If the engine track writes an indel record's whole ID column, or gives each ALT its
own ID where the ID column lists one per ALT, change those three.

### Added

- **Q, the code of a variant position whose reference base and SNP ALTs hold all four bases, is read
  as every base wherever CRISPRme reads a variant code; N stays a gap that matches nothing.**
  CRISPRitz 2.10.0 writes Q at such a position (decision B1), where 2.9.1 coded the reference base
  and the ALT of the position's last split record alone or, for one VCF line listing three SNP ALTs,
  wrote N and lost every target through the position. The SNP post-analysis decomposes a Q column
  into each ALT's carriers and loads the dictionary keys of a target whose only code is Q (it took
  such a target for a reference row, with Q in its DNA column and no carrier); the PAM filter of a
  PAM-less or wider index keeps a target with Q in its PAM; the SNP overlay copies Q onto the
  fake-indel genome's flanks and the indel path resolves it per carrier; the site rule CRISPRme
  checks an index's trees against (`crispritz_bins`) reads a FASTA holding Q and matches it to every
  PAM base; and the population-level lower bound (`twopass_find`), the cluster sort and the guide
  extraction from a sequence take it. A Q stands for every base, so it cannot say which position it
  codes: at a bulge boundary, where the target's column can be read against a neighbouring position
  (#139), a Q column admits that position's ALTs only where they and its reference base hold all
  four bases, an ALT that no sample carries counted, as CRISPRitz codes it (a sites-only dataset's
  admits every ALT its registry lists, as each of its columns does: see Fixed), and of them only
  the bases: an ALT that is no base, such as the `*` of a spanning deletion (in a GATK joint call),
  which the dictionary of CRISPRitz 2.10.0 lists, is no allele of a Q, as of any other code (with
  the dictionary, a row substituted it, `*` in its DNA column). The dictionary lists such an ALT;
  the Tier-0 registry holds no record of it, and lists the positions of four bases in its manifest
  (`four_base_positions`), so that the tiers alone, the layout of every published index, read such
  a position as the dictionary does. Checked through the real CLI on an HPRC kit
  (hg38 chr22:1-17,950,000 and the normalised HPRC records up to 17,949,000; TAGGGAAGTTAGGGAGGTTA,
  NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per sample): of its 17 positions of four bases,
  chr22:17,877,081 holds T>A (8 copies), T>C (none) and T>G (46). On the tiers alone no row named
  it, and the table held 432,738 rows against the dictionary's 439,108; it is now the dictionary's,
  byte for byte, T>A on 648 rows (8 carriers) and T>G on 6,220 (44). Each of the search's 8,129 Q
  columns, 8,108 of them in an alignment with a bulge, is read against a position of four bases.
  Covered by `PostProcess/test_site_alleles.py`, with the real `new_simple_analysis.py` on tiers and
  dictionaries the repository's own compilers build.

### Changed

- **Two VCF records of one REF>ALT at a position are one allele: their carriers are merged, a
  haplotype slot carrying the ALT if either record puts it there, with the rsID and AF of the more
  frequent record (decision B2).** CRISPRitz 2.10.0 merges them in the dictionary; CRISPRme's
  dictionary readers (the Tier-0 and Tier-1 compiles, the dictionary path of the SNP
  post-analysis, the indel path) merge, by the same rule, two entries of one ALT, or a sample
  listed twice in one entry, that reach them from another enrichment. Up to 2.8.1 two entries of
  one ALT stopped the Tier-1 compile of the whole chromosome (`duplicate (pos, alt) record`), and
  up to CRISPRitz 2.9.1 the later record's entry replaced the earlier's. A sample that the merge
  leaves with two ALTs of the position on one slot is resolved without phase (see Fixed). Covered
  by `PostProcess/test_site_alleles.py`.
- **The input validator says what CRISPRme does with a VCF line that lists several ALTs, and with
  positions that hold several records, and looks for such lines on every run.** It warned only at 10
  ALTs or more, in the opt-in `--full_input_validate`, of an allele-index substring match the
  enricher stopped making in CRISPRitz 2.8.4 (CRISPRitz #44); it now warns on any line of two or
  more ALTs, in the full scan and in the first 1,000 records of each VCF on every run. A position's
  split records (`bcftools norm -m -any`, one ALT each, the published panels' form) no longer warn
  as dropped duplicates; two records of one REF>ALT warn that their carriers are merged; a line of
  several ALTs warns that each of them is given the line's first ID, and that split, each ALT's
  record holds the whole ID column, which a SNP's rows give whole and an indel's by its first ID
  (see the known limit below). `gnomAD-converter --multiallelic` says, in its help, on stdout and in
  the README, what it writes: a position's SNP records joined into one line and its indel records
  into another (`bcftools norm -m+`), with the converter's genotypes, which are unphased (`1/2`),
  and how CRISPRme reads them. Covered by `PostProcess/test_validate_inputs.py`.
- **A known limit, stated with its scale: an indel row gives its record's first ID, and so does each
  ALT of a VCF line that lists several ALTs.** CRISPRitz `add-variants` (2.10.0, as 2.9.1) writes
  the first ID alone into the indel log for every indel record, whatever the number of its ALTs, and
  into the per-sample dictionary for each ALT of a line of several; no reader of either can recover
  the others. A SNP of a record of one ALT keeps the record's whole ID column, which CRISPRme gives
  as one item (`rs1;rs2`). So where the ID column lists one ID per ALT (`rs1;rs2`, as
  `bcftools norm -m+` and `gnomAD-converter --multiallelic` join them), the rows of the second and
  later ALTs give the first's (`rs1`) as their `Variant_rsID`, with the tiers, with the dictionary
  or both; and a record of a split VCF (`bcftools norm -m -any`, the published panels' form) keeps
  the line's whole ID column, which its SNP rows give whole and its indel rows by the first ID. The
  mega's VCF lists two IDs for 68,017 of the 222,283 PASS indel records of chr22
  (`22:22652455:TTATG:T;rs1372713682`, whose rows give `22:22652455:TTATG:T`), the HGDP + 1000
  Genomes panel's for 8 of its 240,557; the sites-only (mega) registry builder, which names the
  mega's SNPs, gives each ALT of a line its own ID where the ID column lists one per ALT. Lines of
  several ALTs come from a user's VCF, HF's raw HPRC `vcfwave` VCFs (33,993 on chr22, whose IDs are
  graph IDs, not rsIDs), the 2019 HGDP data `complete-test` uses, and
  `gnomAD-converter --multiallelic`. METHODS section 8 and the input validator's warning on such a
  line say so. Checked through the real CLI with CRISPRitz 2.10.0: a line `rsB1;rsB2` A>G,C gives
  the rows of C `rsB1`; and the review's synthetic set, in all four layouts (tiers, tiers and
  dictionary, dictionary, population-level), gives a record `rsI1;rsI2` AA>A of one ALT `rsI1` on
  its indel rows and a record `rsS1;rsS2` A>C `rsS1;rsS2` on its SNP rows.

### Fixed

- **A VCF line that lists several ALTs is read one ALT at a time, each carrier's genotype against
  that ALT's own allele index (decision B3).** Such a line keeps its raw genotype (`1|2`) under
  every ALT its sample carries, and the readers took "1" for "this ALT": the Tier-0 registry
  counted the second ALT's copies as reference (on chr22 of HGDP + 1000 Genomes with its split
  records joined into such lines, 42,952 of its 43,389 positions of several SNP ALTs had an ALT
  with a wrong global AC: chr22:10,570,455 T, AC 2 in the registry against 2,977 of 8,182 in the
  VCF, and Variant_MAF was wrong on 3 to 10 rows per guide); the observed path guessed an ALT's
  index from its own carriers' most frequent token (`_alt_index_for_record`, wrong for 449 of
  45,837 ALTs, 2 or 3 carriers lost per search); and the dictionary path put a `1|2` sample's
  first ALT on both its slots. The allele index of an ALT is
  now the token all its carriers share, narrowed by what the position's other ALTs need; where the
  genotypes cannot tell (every carrier of two ALTs is `1|2`), the carrier holds one copy on an
  unknown slot and is resolved without phase (PUTATIVE); the population summary of a phased panel
  decides it without phase too (it read the unknown slot as slot 0, a cis copy of the row of each
  ALT of the line). The Tier-1 store holds each ALT's genotypes read so (`1|0` and `0|1` for a
  `1|2`), and a store of an earlier release is read the same way. The sites-only (mega) registry
  builder and its indel AF companion read such a line one ALT at a time, each with its own AF and,
  where the ID column lists one per ALT, its own ID; they skipped it whole (the dictionary gives
  each ALT the line's first ID: see the known limit under Changed). A position of split or biallelic
  records reads as before. These lines come from HF's raw HPRC `vcfwave` VCFs (33,993 on chr22), the
  2019 HGDP data `complete-test` uses, and `gnomAD-converter --multiallelic`. Covered by
  `PostProcess/test_site_alleles.py` and `PostProcess/test_build_mega_registry.py`.
- **A haplotype slot holding two ALTs of one position is no haplotype: the sample is PUTATIVE
  there, a carrier of each ALT's row, and no row names two ALTs of one position (decision B4).**
  With every ALT of a split site in the dictionary, a panel's split records phased one by one put
  two ALTs on one slot (672 sample-positions on chr22 of HGDP + 1000 Genomes): the per-sample
  analysis made it one CONFIRMED haplotype holding both, a row naming two alleles at one position
  (chr22:10,661,629 G and T for NA18869, CONFIRMED in 190 alignments), and listed the sample under
  neither ALT's row. The PUTATIVE branch now takes at most one ALT of a position in each
  combination, the union above `CRISPRME_IUPAC_CAP` variants included: one union per choice of ALT,
  at most 16 a sample, and past that one, holding at each such position the ALT nearest the guide
  (uncapped, NA19042, PUTATIVE in the IGL window chr22:22,881,180-22,881,205 and carrying two or
  three ALTs at 16 of its 26 positions, gave 8.5 million unions there, and 90 GB; now 331 distinct
  haplotypes in all for the window's 790 carriers, in 0.1 s). The population summary decides such
  a carrier without phase, as it now decides a carrier whose genotype at a variant of the row is
  unphased; the SNP+SNP co-occurrence companion counts positions, not variants, so two ALTs of one
  position are no co-occurrence; and a sites-only dataset's co-located union takes one ALT of a
  position (one union per choice, at most 16 a window). The phase label is per row, so a row such a
  sample shares reads PUTATIVE for all its carriers. Covered by
  `PostProcess/test_observed_haplotypes.py` and `PostProcess/test_site_alleles.py`.
- **The dictionary path keeps only the ALTs a target's own column holds, as the observed path does
  (#139); a sites-only dataset substitutes every ALT its registry lists, as up to 2.8.1.** The
  dictionary path, and a genotyped panel's registry read alone (`download --no-genotypes`),
  substituted every ALT listed at a position, an ALT the code the engine matched does not hold
  included, so the two paths could report different rows for one window: those sources list what the
  enrichment coded, so such an ALT is another position's. A sites-only (mega) registry is read from
  its VCF, every record of it, and lists ALTs the enrichment never coded: on an index built by
  CRISPRitz 2.9.1 or earlier, as every mega index published until the rebuild, the earlier records
  of a split site (the enrichment codes the last PASS record's), and on any index non-PASS records.
  Each is the panel's own allele at its own position, so a sites-only dataset substitutes every one
  at a variant column, whatever its code, a Q included, as 2.8.1 did. Applied there, the dictionary
  path's filter lowered the highest CFD of 114 sites in a search of the published
  `NRG_3_hg38+hg38_mega` (chr22; the README's guide, ACTGA, TRAC, EMX1 and VEGFA site 3; NRG, 4
  mismatches + 1 DNA + 1 RNA bulge, population-level; against 2.8.1): 110 for a PASS ALT the 2.9.1
  code left out (chr22:10,752,820 C>A, whose later record C>T is the coded one) and 4 for a non-PASS
  one, 22 of them from 0.1 or more and 2 from 0.5 or more (VEGFA site 3's chr22:22,893,912 from
  0.585 to 0.327, TRAC's chr22:22,900,049 from 0.524 to 0.127); this release's search gives 2.8.1's
  20,338 sites with 2.8.1's highest CFD at every one, and no row naming two ALTs of one position
  (2.8.1: 4 of the integrated table's), with nothing on stderr (CRISPRitz 2.9.1, the 2.8.1 image). A
  target that only such an ALT brings within the search's budget is still not found on that index:
  the mega index rebuilt by CRISPRitz 2.10.0 codes every PASS ALT (METHODS section 2). Covered by
  `PostProcess/test_site_alleles.py`.
- **On a row of several variants, the rsID and AF lists follow its variants: one rsID and one AF a
  variant, in the order of `Variant_info_genome`, and two variants of one frequency give it each.**
  The merge of a cluster's rows (`remove_contiguous_samples.py`) sorted the variant, rsID and AF
  lists each on its own, as text, and kept each value once (the known limit 2.8.1 stated: on ACTGA,
  953 of 1,967 rows of two or more SNPs out of order). It now keeps (variant, rsID, AF) together; a
  variant two datasets of a search give with two rsIDs or AFs lists both pairs, one after the other,
  and an rsID field holding two IDs (a VCF ID `rs1;rs2`, which the dictionary holds as `rs1,rs2`)
  reads `rs1;rs2`, one item. The report's MAF, the lowest AF of the row, is unchanged. Covered by
  `PostProcess/test_site_alleles.py`, `PostProcess/test_generate_report.py` and
  `PostProcess/test_representative_tiebreak.py`.
