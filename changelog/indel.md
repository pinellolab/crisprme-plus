<!-- indel track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

# [2.8.1] entries from the indel track (SNP+indel off-targets)

### Fixed

- **An indel off-target near a SNP was reported with SNP alleles its carriers do not have, and with
  PAMs the search never matched.** Since 2.5.0 the build writes the IUPAC code of every SNP onto the
  fake-indel genome the indel search runs on, and the search engine matches a code when any of its
  alleles fits, so one raw indel target stands for every combination of the SNP alleles around the
  indel. The indel post-analysis then turned each code into one base, column by column, without
  looking at any genotype: the guide's own base when it was one of the code's alleles, otherwise the
  code's first allele. A row therefore paired one individual's indel with SNP alleles that only
  other individuals carry, scored as if real, under the indel's carriers: on the default HGDP + 1000
  Genomes index, guide `ACTGAAATCTGTAAGCAGGC` (NRG, 4 mismatches + 1 DNA + 1 RNA bulge) reported
  chr18:12386282 at CFD 0.708, fifth in the report's top-1000 table and 18th in the recommended
  validation panel, as "1 cis carrier": NA19072, who is homozygous reference at the SNP the
  displayed sequence uses (its only carrier, NA20541, does not carry the deletion); on NA19072's
  real haplotype the site is 5 mismatches + 1 bulge, past the search's limit. Its mismatch counts
  were understated the same way. And a code at a PAM position always read as its first allele: an R
  (A/G) at NRG's last base was shown as A, a PAM the search never matched, with a CFD of about 0
  (1,088 rows of that search's highest-CFD block). On 9,255 of its 13,622 highest-CFD indel rows no
  listed carrier's real haplotype, built from the phased VCF, holds the reported sequence, mostly
  for this reason (the next two entries and METHODS §8 cover the others). Each carrier haplotype of
  the indel now gets its own allele at every code: the SNP alt it carries, else the reference base,
  several alts on one VCF line each read against its own allele index. Haplotypes left with the same
  sequence form one row that lists only their individuals and names the SNP alts they hold next to
  the indel (`Variant_info`, rsID, MAF), with its mismatches recounted, its CFD, CRISPR-Bulge and
  risk scores computed on that sequence, and its PAM creation flagged from it; a sequence whose PAM
  does not satisfy the searched motif, or that is past `--mm` or `--max-total-edits` (now passed to
  the indel post-analysis), is dropped, so a raw target no individual has gives no row. A row's
  carriers stay on its sequence through the merge of a cluster's rows, which would otherwise add
  those of another alignment of the site naming the same variants; so do those of a row of an
  equal-length multi-base record (an MNV, such as AC>TG), which the search treats as an indel. An
  individual whose genotypes cannot put an allele on the indel's haplotype (an unphased call, a
  missing allele, no phased indel genotype, a haploid call against a diploid one; every carrier of
  an MNV, which the phased indel genotype store leaves out) gets the strongest target its own
  alleles allow, and is never listed for an allele it does not carry; with no genotypes at all (a
  sites-only panel, `hg38_mega`), the row holds the strongest target the codes allow (its PAM
  always one the motif allows) under the carriers the
  indel log lists, as before. Both are PUTATIVE in the SNP+indel co-occurrence companion
  (`*.indel_snp_cooc.tsv`), which now lists the SNP alts the indel post-analysis's rows require, a
  SNP that creates the PAM included, with those rows' carriers (CONFIRMED when each was resolved on
  its phased haplotype), and no longer a SNP set no carrier's haplotype holds. The companion is
  written before the merge picks each cluster's alignments, one line per engine alignment, required
  SNP set and phase (an alignment that needs no SNP writes no line), so a population-level run,
  which writes no alternative alignments, can have lines with no row in its tables. A classic
  dictionary install (an
  index built without `--samplesID`, as the web Settings page builds them, so with no Tier-1
  genotype store) reads the genotypes from the per-sample SNP dictionary: there, the codes used to
  stay in the rows unresolved, scoring CFD 0, so the real SNP+indel targets were lost. With
  `CRISPRME_INDEL_SNP=0` at search time each code reads as the reference base, as on a fake-indel
  genome built without the overlay. A target with no code (every target of a fake-indel genome built
  without the overlay) that holds its indel, as the entry "An indel row could hold no base of its
  indel" defines it, is not changed by this fix. **Affected:**
  every search since 2.5.0 on an index that carries the overlay, which every published variant index
  and every index `build-index-only` builds do: from the web or the command line, per-sample and
  population-level (the indel post-analysis does not depend on the mode), with any PAM. A search that
  enriches its dataset on demand (`CRISPRME_ALLOW_ONDEMAND_BUILD=1`, as the benchmark mode of
  `complete-test`) builds its fake-indel genome without the overlay, so no code reached its
  post-analysis. Under a PAM whose motif allows the base the misreading produced (NRR allows A at its
  last position) the PAM stayed legal while the CFD collapsed, so nothing showed. Reference and
  SNP-only off-targets are not affected, though a cluster whose indel candidate goes away shows its
  next best in the results and its `Other_motifs` count changes. Measured on the search above, rerun
  genome-wide with the fixes of this entry and the next two: its highest-CFD indel rows went from
  13,622 to 5,047, of which 5,037 hold their sequence on a listed carrier's real haplotype, read
  literally from the phased VCF, and every one of those spans its indel; the other 10 hold a SNP at
  a split multi-allelic site, which the enriched genome and its dictionaries read as reference
  (METHODS §8), as does one of the 103 carriers of one more row; no ALT row shows a PAM NRG does
  not allow (3,581 over the table's three blocks before); chr18:12386282 is now the SNP row of its
  real carrier, NA20541, at CFD 0.14, out of the top-1000 table and the panel; that table lost 66
  indel rows and gained 2 (chr5:170503088 and chr10:103001704), from 88 to 24, the recommended
  panel went from 259 to 234 sites, the CFD >= 0.5 tier from 67 to 57, the CRISPR-Bulge sensitive
  tier from 162 to 156 and the <= 2 mismatches + bulges tier from 9 to 7; the SNP+indel companion
  went from 590 to 336 lines, all CONFIRMED and the same in both analysis modes; all 755,576
  reference and SNP-only candidates are unchanged, and so are the two perfect matches and the
  multiple-perfect-match banner.
  HPRC is the only published index with MNVs (chr22: 9,079 of its 195,813 indel-search records);
  the search above on HPRC's chr22 lists 8 MNV row blocks per block, and the MNV part of the merge
  fix changes none of their carriers.
  The search takes as long as before and peaks at the same memory, its search stage's; an indel
  post-analysis worker holds about 2% more. **What to do:** run such searches again; results of
  earlier releases are corrected only by running the search again. Covered by
  `PostProcess/test_indel_haplotype_resolution.py`, `PostProcess/test_indel_overlaid_rows.py`,
  `PostProcess/test_representative_tiebreak.py`, the synthetic SNP+indel oracle through the real CLI
  (`test/indel_snp_oracle`, run by `validate-benchmarks.yml`) and a release gate,
  `scripts/check_alt_pam.py` (no ALT row may show a PAM the searched motif does not allow), run on
  the chr22 benchmark outputs.
- **An individual whose haplotype carries a second indel inside an indel off-target was listed for
  the sequence of the first one.** The fake-indel genome holds one indel per window, so such a
  haplotype does not have the window's sequence; on the synthetic oracle a sample carrying two
  deletions on one haplotype was listed on both windows' rows, for sequences neither of its
  haplotypes has. It is no longer listed for them, nor when its second indel is inserted right next
  to the first, nor when, on the same phased slot, an indel overlaps the first at another position
  (one complex allele that the VCF splits into overlapping records). A variant only on bases the
  indel deletes changes nothing. Another allele at the indel's own position, an indel or a SNP, is
  an alternative to it (a haplotype holding both would be one record), so a slot phased with both is
  a phasing artefact: on chr22 of the HGDP + 1000 Genomes panel, two indels at one position share a
  slot for 40,679 carrier pairs against 2,463,324 on opposite slots, a SNP and an indel for 2,404
  against 167,951. And an unphased carrier of two overlapping records is a multi-allelic
  heterozygote, one allele per haplotype. All of these stay listed. This needs the phased indel
  genotype store, which the build writes with the SNP overlay; a row left with no carrier is
  dropped. Its own target, which needs both indels, is still not searched (METHODS §8). On the search
  above, 693 targets with no SNP code lost such carriers, and no row left lies off a listed
  haplotype because of a second indel.
- **An indel row could hold no base of its indel.** The engine reports every target on a fake-indel
  window, also one wholly in a flank, whose sequence is plain reference (or SNP) sequence that needs
  no indel. The post-analysis computed a row's Start with no shift for the indel and kept the row when
  that range met the indel's, so a target in the downstream flank of a deletion passed: a row with
  the indel's name and carriers, at a coordinate the indel shifts, that could take the real row's
  place in its cluster. On the search above, the top-1000 table showed at rank 628 chr1:69462019 as a
  2-carrier off-target of a deletion and a SNP, hiding the SNP row chr1:69462021, chr1_69462039_A_G,
  with its 3,459 carriers, so a common off-target read as a rare one (it is now at rank 588), and
  at rank 983 chr8:28788278, the sequence of the reference row chr8:28788272 reported again 6 bp off.
  A raw target now has to hold a base the ALT adds or replaces, or the two bases on either side of the
  deleted ones; one wholly inside a long insertion does, and stays. **Affected:** every search with a
  variant dataset that has indels, in every release, on any index, with or without the overlay and
  with `CRISPRME_INDEL_SNP=0` too. No off-target goes with these rows: on the search above, each of the
  1,566 distinct rows this removes has a reference or SNP row at its site that holds its sequence
  (1,560) or names the same SNP alts (6). **What to do:** run such searches again. Covered by the
  same tests.
- **An indel off-target with an RNA bulge was reported one base early, and the reference sequence
  shown for one with both bulge types was a base off.** For an indel row whose bulges are all RNA
  bulges, the indel post-analysis took the search engine's cluster position for the row's start on
  the + strand (on the - strand with a 5' PAM such as Cas12a's), so its `Start_coordinate` was one
  base per RNA bulge before its first aligned base, where the reference and SNP rows of the same
  alignment start, and its variant's position in `Variant_info_spacer+PAM` was off by as much. The
  reference sequence of an indel row (`Aligned_protospacer+PAM_REF`) was read with one base too many
  per RNA bulge and was given its gaps only on a row whose bulges are all RNA bulges: a row with an
  RNA and a DNA bulge (`RNA,DNA`) showed it without the gap and one base too long, out of register
  with the guide from the RNA bulge on, and with a 5' PAM an RNA-bulge row's lost its first PAM
  base. The columns computed from it followed: the REF seed and non-seed mismatches + bulges,
  `CFD_score_REF`, the CFD risk score, `PAM_creation`, and `CFD_score`, which is the higher of the
  REF and ALT scores. Each indel row now starts at its first aligned base and shows as REF the
  reference bases its target covers from there, gapped as its target. Its cluster position, which
  the merge of a cluster's rows sorts and groups on, moves with its start: an RNA-bulge indel row on
  the + strand no longer comes first among rows that tie exactly, and, as for SNP rows, the merge
  lists once the alignments of a site that name the same variants with the same carriers (the best
  of them: the next entry). Measured
  on the search above, rerun genome-wide with this fix alone: before it, 238 of the 5,047 indel rows
  of the highest-CFD block started one base early and 4,467 showed a REF sequence that is not the
  reference at their start; now each of the 5,045 shows the reference at its start. On those rows
  `PAM_creation` is now set on 815 and no longer set on 240 (all on the + strand, where the old
  REF's PAM was read one base off), `CFD_score_REF` changed on 785, and `CFD_score` went down on 218
  and up on 16 (chrX:42118270, 0.149 -> 0.017: the old value was the misaligned REF's, so the row
  left the report's CFD >= 0.05 table, as did 5 others). Four clusters now report the row of SNPs
  that has the same sequence, CFD and carriers as their indel row, and two indel off-targets
  that sat in a reference off-target's cluster because of their early start are rows of their own
  (chr3:55057499, CFD 0.087, 20 cis carriers; chr16:12067029, CFD 0.069, 1). The alternative
  alignments went from 722,564 to 722,347 rows, and 39 rows of the integrated results changed only
  their `Other_motifs` count. The report's tiers keep their sites, some one base later, except for
  those CFD >= 0.05 rows, the CFD >= 0.05 table gaining the two new rows and chr2:33157279 (whose
  REF now scores 0.1), and one tie at the top-1000 table's cut-off (CFD 0.188), now broken the other
  way. The SNP+indel companion keeps its 1,618 lines, one per alignment, the off-target start one
  base later on those rows; counted once each, as in the first entry and in the report, it went
  from 336 to 339 lines, all CONFIRMED, as some of those rows no longer share their line with
  another alignment of their site. The same search of the README's example guide,
  `CTAACAGTTGCTTTTATCAC` (the README search), keeps its 1,672 lines, 351 counted once each
  before and 363 now.
  Apart from `Other_motifs` counts and those four rows of SNPs, no reference or SNP-only row
  changes, nor do the perfect matches and the banner. The search takes as long as before, with the
  same peak memory. **Affected:** every search that allows RNA bulges (`--bRNA` 1 or more; the
  recommended default allows one) with a variant dataset that has indels, in every release, with or
  without the overlay and with `CRISPRME_INDEL_SNP=0` too; with a 5' PAM, every indel row with an
  RNA bulge. **What to do:** run such searches again. Covered by
  `PostProcess/test_indel_row_coordinates.py` (one and two RNA bulges, both bulge types, a DNA bulge
  and none, on both strands, with an NRG and a TTTV PAM, with `CRISPRME_INDEL_SNP=0` too) and the
  SNP+indel oracle, which now checks every indel row's REF sequence against the reference at its
  start.
- **Of the alignments of one site that name the same variants with the same carriers, the merge
  reported the first it read, not the best.** The merge of a cluster's rows lists such alignments of
  an indel row once and drops the others from both tables. It kept the first of them in the order
  its input is sorted, where alignments tied on score and on mismatches + bulges come in text order
  (an `RNA` alignment before an `RNA,DNA` one), so with the default `--sorting-criteria`
  (mismatches + bulges, then mismatches) the alignment reported could have one mismatch more and one
  bulge fewer than a dropped one, and in the fewest mismatches + bulges block another CFD. It now
  keeps the best of them by the order that picks each block's representative; among alignments
  that tie exactly there, the one that order's last tie-break puts first (the aligned sequence first
  in text order), whatever order they came in. Measured on the search above, this fix against the
  code before it: of its 182,367 integrated rows, 243 show another alignment of the same site,
  variants and carriers. In the highest-CFD block that is 56 rows, 12 with one mismatch fewer and
  one bulge more at the same CFD and 44 an alignment of equal rank; in the fewest mismatches +
  bulges block 232, 94 with one mismatch fewer and one bulge more (47 at another CFD) and 138 of
  equal rank (88 at another CFD, as that block's order does not look at CFD), and in 4 clusters the
  better indel alignment now comes before the cluster's SNP row, which moves unchanged to the
  alternative alignments; in the CRISPR-Bulge block 135 of equal rank. The alternative alignments
  keep their 722,347 rows, 191 of them changed in the same ways. Against the run before the
  previous entry's fix, no reported alignment has more mismatches: the 4 rows of the highest-CFD
  block and 18 of the fewest block that fix had moved to the worse alignment show the one reported
  before it again. The README search changes 323 of its 170,657 integrated rows the same way (of
  the 3 + 15 rows that fix had made worse, 16 show the alignment of before and 2 one of equal rank
  to it), and the search above at population level 241 of 211,830 (with no alternative alignments
  there, the 5 SNP rows that give the fewest block to a better indel alignment leave it). The
  report's tables keep their sites (one row of the <= 4 mismatches + bulges table now starts
  elsewhere); no carrier, rsID, AF, row count, reference or SNP-only row changes, the SNP+indel
  companion is the same, and the summaries and the raw observed-haplotypes sidecar follow the
  tables. **Affected:** every search with a variant dataset that has indels, in every release, with
  or without the overlay and with `CRISPRME_INDEL_SNP=0` too (until the previous entry's fix a +
  strand indel row with RNA bulges had a cluster position of its own and did not meet the site's
  other alignments). **What to do:** run such searches again.
  Covered by `PostProcess/test_representative_tiebreak.py`.
