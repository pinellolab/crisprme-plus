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
  (1,088 rows of that search's highest-CFD block). Of its 13,622 highest-CFD indel rows, 8,260 lie
  on no listed carrier's real haplotype (built from the phased VCF): for 6,572 nothing else explains
  it; in the others a listed carrier also holds a second indel there (1,317, the next entry) or a
  SNP at a split multi-allelic site (371, METHODS §8). Each carrier haplotype of the indel now gets
  its own allele at every code: the SNP alt it carries, else the reference base. Haplotypes left
  with the same sequence form one row that lists only their individuals and names the SNP alts they
  hold next to the indel (`Variant_info`, rsID, MAF), with its mismatches recounted, its CFD,
  CRISPR-Bulge and risk scores computed on that sequence, and its PAM creation flagged from it; a
  sequence whose PAM does not satisfy the searched motif, or that is past `--mm` or
  `--max-total-edits` (now passed to the indel post-analysis), is dropped, so a raw target no
  individual has gives no row. An individual whose genotypes cannot put an allele on the indel's
  haplotype (an unphased call, a missing allele, no phased indel genotype, a haploid call against a
  diploid one) gets the strongest target its own alleles allow, and is never listed for an allele it
  does not carry; with no genotypes at all (a sites-only panel, `hg38_mega`), the row holds the
  strongest target the codes allow (its PAM always one the motif allows) under the carriers the
  indel log lists, as before. Both are PUTATIVE in the SNP+indel co-occurrence companion
  (`*.indel_snp_cooc.tsv`), which now lists, per reported row, the SNP alts it requires, a SNP that
  creates the PAM included, with that row's carriers (CONFIRMED when each was resolved on its phased
  haplotype), and no longer a SNP set no carrier's haplotype holds. A classic dictionary install (an
  index built without `--samplesID`, as the web Settings page builds them, so with no Tier-1
  genotype store) reads the genotypes from the per-sample SNP dictionary: there, the codes used to
  stay in the rows unresolved, scoring CFD 0, so the real SNP+indel targets were lost. With
  `CRISPRME_INDEL_SNP=0` at search time each code reads as the reference base, as on a fake-indel
  genome built without the overlay. A target with no code, every target of a fake-indel genome built
  without the overlay, is reported as before. **Affected:** every search with a variant dataset that
  has indels since 2.5.0, from the web or the command line, per-sample and population-level (the
  indel post-analysis does not depend on the mode), with any PAM: every published variant index and
  every index `build-index-only` builds carry the overlay. Under a PAM whose motif allows the base
  the misreading produced (NRR allows A at its last position) the PAM stayed legal while the CFD
  collapsed, so nothing showed. Reference and SNP-only off-targets are not affected, though a
  cluster whose indel candidate goes away shows its next best in the results and its `Other_motifs`
  count changes. Measured on the search above, rerun genome-wide with the fix: its highest-CFD indel
  rows went from 13,622 to 5,311, of which 5,303 check out on a listed carrier's real haplotype and
  the other 8 hold a SNP at a split multi-allelic site, which the enriched genome and its
  dictionaries read as reference (METHODS §8); no ALT row shows a PAM NRG does not allow (3,581 over
  the table's three blocks before); chr18:12386282 is now the SNP row of its real carrier, NA20541,
  at CFD 0.14, out of the top-1000 table and the panel; that table lost 62 indel rows and gained 2
  (chr5:170503088 among them), the recommended panel went from 259 to 234 sites, the CFD >= 0.5 tier
  from 67 to 57, the CRISPR-Bulge sensitive tier from 162 to 156 and the <= 2 mismatches + bulges
  tier from 9 to 7; the SNP+indel companion went from 590 to 423 lines, all CONFIRMED and the same
  in both analysis modes; all 755,576 reference and SNP-only candidates are unchanged, and so are
  the two perfect matches and the multiple-perfect-match banner. The search takes as long as before
  and peaks at the same memory, its search stage's; an indel post-analysis worker holds about 3%
  more. **What to do:** run such searches again; results of earlier releases are corrected only by
  running the search again. Covered by `PostProcess/test_indel_haplotype_resolution.py`,
  `PostProcess/test_indel_overlaid_rows.py`, the synthetic SNP+indel oracle through the real CLI
  (`test/indel_snp_oracle`, run by `validate-benchmarks.yml`) and a release gate,
  `scripts/check_alt_pam.py` (no ALT row may show a PAM the searched motif does not allow), run on
  the chr22 benchmark outputs.
- **An individual whose haplotype carries a second indel inside an indel off-target was listed for
  the sequence of the first one.** The fake-indel genome holds one indel per window, so such a
  haplotype does not have the window's sequence; on the synthetic oracle a sample carrying two
  deletions on one haplotype was listed on both windows' rows, for sequences neither of its
  haplotypes has. It is no longer listed for them, nor when its second indel is inserted right next
  to the first, nor when, on the same phased slot, an indel overlaps the first or a SNP alt sits on
  its anchor base (one complex allele that the VCF splits into overlapping records). A variant only
  on bases the indel deletes changes nothing, and an unphased carrier of two overlapping records, a
  multi-allelic heterozygote, stays listed. This needs the phased indel genotype store, which the
  build writes with the SNP overlay; a row left with no carrier is dropped. Its own target, which
  needs both indels, is still not searched (METHODS §8). On the search above, 1,221 targets with no
  SNP code lost such carriers, and no row left lies off a listed haplotype because of a second
  indel.
