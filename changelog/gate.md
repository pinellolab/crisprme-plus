<!-- gate track: [2.9.0] entries for CHANGELOG.md (the carrier gate, plan 2.3; decisions D1, D2,
D3, D5, D6; the F2 and greedy-row fixes the 2.8.1 text promised for this release). -->

**For the merge, not for the changelog.** 2.8.1's `### Known issues` entries on SNP rows that
list a carrier whose haplotype holds an indel inside the target, and on the greedy rows of a
`--per-sample` search on an index without the dict-less tiers, both say "fixed in the next
release": the two `### Fixed` entries below are those fixes. The runs measured are of this
branch's code in the 2.8.1 release image (CRISPRitz 2.9.1, the 2.8.1 indexes), beside the same
searches with the 2.8.1 code; the harness check re-gated gate38a's arm C9 (the 2.8.1
population-level run of the README guide).

### Added

- **The carrier gate: a population-level row's variant combination is checked against phased
  haplotypes, and the report ranks a combination no haplotype carries after every other row, never
  dropping it.** A population-level row takes, per variant window, the alleles that match the guide
  best and lists every carrier of any of them, so two or more of them can sit on no one haplotype:
  no individual has the sequence the row reports, and, being the best-matching, it ranked first.
  After writing the integrated results the search now reads each row's variants (each of its three
  alignments) in the dataset's own phased genotypes, SNPs and indels
  (`Dictionaries/genotypes_<vcf>/`, `indel_genotypes_<vcf>/`;
  `PostProcess/haplotype_prioritizer.py`), and gives each a tier: **A**, two or more variants one
  phased haplotype carries together (states `one haplotype`, and `phase not known` where only
  unphased calls, or two ALTs of one position on one slot, show it), or that a haplotype reaches
  the alignment's CFD without (states `same score`, `same score, phase not known` and `same score
  as the reference`): an allele the worst case takes can leave the CFD the row reports (its three
  decimals), the alignment's mismatches and its PAM as they are -- one at a PAM base the search's
  PAM leaves free, a mismatch for a mismatch the
  CFD table scores alike -- and then a haplotype carrying the rest, and nothing else inside the
  alignment that changes them, has the row's off-target at the row's score (the companion's
  `Scoring_variants` column names the variants the score needs). The guide's side of an alignment
  writes its PAM as `N`, which every base matches, so the PAM region is read against the search's
  own PAM (the search gives the gate its PAM file): without an allele the PAM needs there is no
  off-target at the site at all. The CFD gives this state to the alignments it scores (the
  highest-CFD one, which the report ranks by, and the fewest-mismatches one), not to the one ranked
  by CRISPR-Bulge score, which the gate does not read (the model is not run again) and which an
  allele the CFD never reads does change; **B**, one variant; **C1**,
  two or more variants, each in the panel, that none of its haplotypes carries together nor
  reaches the row's CFD without, with E,
  the haplotypes expected to carry them together if they were independent (below 3, the state is
  `too rare to tell`), or variants, a single one too, every haplotype carrying which has an indel
  that takes the target's bases off it (`indel in the target`); **C2**, a variant not in the panel,
  which cannot say (the frequency bound, the product of the variants' frequencies >= 1e-06 and the
  rarest >= 0.001, is shown and applied to nothing). The gate changes no row: it writes
  `..._carrier_gate.tsv` beside the integrated results (one line per alignment that uses a variant)
  and `.carrier_gate.json`. The report ranks a C1 row after every other row in the views ranked by
  the score of the alignment that is C1: a row whose highest-CFD alignment is C1 leaves the
  top-1000 by CFD, panels (a) and (b) of Section 2, the validation panel's CFD thresholds and the
  risk card's CFD bands, and is listed in a table of its own after the top-1000
  (`unobserved_worst_case.tsv` holds them all; its caption
  counts the perfect matches among them); a row whose highest-CRISPR-Bulge alignment is C1 leaves
  the top-1000 by that score, panels (c) and (d), and the validation panel's CRISPR-Bulge
  thresholds, as the risk card's
  CRISPR-Bulge bands. The two alignments of a row need not be the same one, so a row of the
  demoted table can keep its place in the views ranked by CRISPR-Bulge score, and the other way
  round. The panel's totals (its off-targets, perfect matches and variant-created
  off-targets) and `run_manifest.json`'s counts are of every row, as the Section-1 matrix's, and a
  row of the panel's table says how many rows its thresholds, the panel and the risk card leave
  out. A, B and C2 rows keep their place, and the `Observed` column says what the gate found (see
  Changed). The inputs card has a **Carrier gate** row: its panels, its counts by tier, and the
  panel's samples by superpopulation (AFR 986, AMR 549, CSA 782, EAS 816, EUR 771, MID 157, OCE 30
  on HGDP + 1000 Genomes), as a haplotype found only in an ancestry with few samples there can be
  missing from it. On by default for a population-level search of a genotyped dataset;
  `--no-prioritize-observed` turns it off, `--gate-min-expected` sets the `too rare to tell`
  threshold (default 3; 0 turns it off). A per-sample search, whose rows name their cis carriers,
  needs none; a dataset with the classic dictionaries alone has no phased store to read, and the
  search and the report say so. On the genome-wide population-level search of the README guide
  (HGDP + 1000 Genomes, NRG, 4 mismatches + 1 DNA + 1 RNA bulge, 198,826 rows) the gate took 31 s
  in the search (35 s and 1.6 GB on 12 threads alone) and gave, by each row's highest-CFD
  alignment, A 4,155 (690 `one haplotype`, 1 `phase not known`, 2,574 `same score`, none `same
  score as the reference`), B 55,191, C1 40,231 (246 where 3 or more haplotypes were expected,
  39,781
  `too rare to tell`, 204 `indel in the target`, 172 of them of one variant) and C2 0. The report's
  top-1000 by CFD and validation panel then hold no row whose highest-CFD alignment is C1, and its
  own table all 40,231, while the
  panel's totals are the search's (198,826 sites, 1 perfect match, 91,289 variant-created
  off-targets) and 43,264 rows are C1 at their highest-CRISPR-Bulge alignment. The integrated
  results are byte for byte those of the same search with `--no-prioritize-observed`. Checked
  against all 8,182 phased haplotypes by a harness that rebuilds every alignment's window on each
  haplotype from the VCF (on the 2.8.1 run of the same search, whose rows are the same file): the
  only alignments left in tier A or B that list a combination no haplotype carries are the 2,574
  of a `same score` state (2,328 at the fewest-mismatches alignment, none at the CRISPR-Bulge one),
  which is what that state
  says; every other A or B alignment lists a combination a haplotype carries, and every C1
  alignment that a haplotype does carry is one with an indel in the target. The harness's stricter
  test (a haplotype that reproduces the alignment within its mismatches with a valid PAM) reaches
  2,572 of the 2,574, and the per-sample run of the same search reports the row's CFD at the locus
  of PS_SAME of them (PS_LOWER below it, PS_NONE with no row within three bases). Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **`--prioritize-observed` checks a sites-only dataset's rows against the phased HGDP + 1000
  Genomes genotypes, and `download --what genotype-oracle` installs them.** The `mega` panel has no
  genotypes, so its multi-variant rows could not be told apart. With `--prioritize-observed` (the
  web form's check box, offered for a sites-only panel once they are installed) the gate reads its
  rows in the phased HGDP + 1000 Genomes genotypes and indel genotypes (the genotype oracle), which
  an installed `NRG_3_hg38+hg38_HGDP1kGP` index already has, which
  `download --what genotype-oracle [--panel HGDP1kGP]` installs without the index (4 min and 28 GB
  on disk from the Hub), and which `CRISPRME_GENOTYPE_TIER` can name elsewhere. A variant the panel
  does not hold makes its row C2, kept in place, with its frequency bound (a SNP's frequency from
  the dataset's registry, an indel's from its indel frequencies, `AF_max` both). On chr22 of the
  mega panel (five guides, NRG, 4 + 1 + 1, 20,338 rows): A 137 (32 of them `same score`), B 4,419
  (1,568 of a variant not in
  the panel), C1 594 (58 where 3 or more haplotypes were expected, 520 `too rare to tell`, 16
  `indel in the target`, 13 of them of one variant), C2 3,519 (1,919 within the bound, 1,600 below
  it), the integrated results byte for byte those of the search without the flag. Without the
  oracle the search says so (a `WARNING` line on stdout), and the rows are reported unchecked.
  Until the Hub has the panel's phased indel genotypes (see Changed), the download says so in one
  line and installs the SNP genotypes, and an indel inside a target then takes no haplotype off.
  Covered by `PostProcess/test_haplotype_prioritizer.py` and
  `PostProcess/test_web_search_reuse.py`.

### Changed

- **The `Observed` column of a population-level row of two or more variants says what the
  carrier gate found, in place of `N carriers (>=1 of K vars, poss. chimera)`.** The row's
  carriers are those of any of its variants, so the cell could not say whether anyone carries
  them together, and "chimera" named a guess. Where the gate checked the row, the cell gives the
  tier of its highest-CFD alignment: `together on N of H <panel> haplotypes` or `carried
  together by N <panel> samples, phase not known` (A); `not together on any of H <panel>
  haplotypes (expected E)`, the same with `; too rare to tell`, or `all N <panel> haplotypes
  carrying them together have an indel in the target` (C1); `not checked: K of N variants not
  in <panel>`, with the variants' frequency product and rarest frequency against the bound where
  the dataset gives them (C2); `the same CFD on N of H <panel> haplotypes, which carry <variants>`
  where a haplotype reaches the row's CFD without every variant (A). A sites-only panel's row,
  `putative` up to 2.8.1, reads so with
  `--prioritize-observed`. A row of one variant reads its carriers as before, unless every
  haplotype carrying it has an indel that takes the target's bases off it, when it reads `all N
  <panel> haplotypes carrying it have an indel in the target` (C1, ranked after the others as
  well: 172 rows of the genome-wide run, 13 of mega chr22). A per-sample row, and the rows of a
  search the gate did not check (or of a 2.8.1 search, read again), read as before. The report's
  legend, the inputs card, README and METHODS list the labels. Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **`publish-index` uploads a variant index's phased indel genotypes on their own too**
  (`indexes/indel_genotypes_<vcf>.tar.gz`, beside `genotypes_<vcf>.tar.gz`; they stay in the
  index's archive as well), for `download --what genotype-oracle`, which installs them without
  the index. Covered by `PostProcess/test_haplotype_prioritizer.py`.
- **Under `--per-sample` a dense window's `High_variant_density_region` note says that the
  haplotypes its carriers have are reported, and so do the report's legend and its note on the
  dense windows near the top-ranked sites.** In either mode the note said, after the window's
  variant count,
  `a greedy worst-case alignment is reported here; additional haplotype alignments may exist`, the
  legend's `High_complexity_region` entry "Such dense windows are reported greedily", and the
  `--report-full` note "a single greedy worst-case alignment is reported for each". A per-sample
  search has reported a genotyped dataset's dense window by the haplotypes its carriers have on the
  dict-less tiers since 2.7.10, and does on the classic dictionaries since this release (see
  Fixed): its rows now say `the haplotypes its carriers have are reported here`, and so does the
  report. A row that names no carrier, a sites-only dataset's, whose rows are the population-level
  analysis's in either mode, keeps the greedy note, and the report says so of it where such a row
  is flagged; a population-level search reads as before. On chr22 of the HGDP + 1000 Genomes panel
  (five guides, NRG, 4 + 1 + 1, per sample) 50 rows carry the note, as many as with 2.8.1. Covered
  by `PostProcess/test_generate_report.py`.

### Fixed

- **A per-sample SNP row no longer names a carrier whose haplotype holding its SNPs also holds an
  indel that changes the target.** Its carriers were read from the SNP genotypes alone, so an
  individual whose only haplotypes holding the row's SNPs carried an indel changing the bases the
  target covers was listed, though its own sequence there is the indel's. The per-sample
  post-analysis now reads the index's phased indel genotypes, on the dict-less tiers and on the
  classic dictionaries alike, and takes such a carrier off; a row left with no carrier is not
  reported. The haplotype is read on the reference from the target's first base on and up to its
  last, as the indel rows read a second indel, so an indel that leaves the target's bases where
  they were (a tandem repeat's unit gained or lost) takes no carrier off; an indel at one of the
  row's SNPs, or over one (another allele of that base: a phasing artefact), does not either, nor
  does an unphased call. A haploid call -- a male's on chrX outside PAR -- is his one haplotype's,
  so his calls in a window are now read as that one haplotype (not as every subset of them, as an
  unphased sample's are) and an indel of his that changes the target takes him off the row,
  whatever else that haplotype holds. On the genome-wide per-sample search of the README guide
  (HGDP + 1000 Genomes, NRG, 4 + 1 + 1) that leaves PSX_ROWS. On chr22 of the HGDP + 1000 Genomes
  panel (five guides, NRG, 4 + 1 + 1,
  per sample), of the 22,628 SNP alignments of 2.8.1's 20,751 rows, 132 left their row (32,227
  carrier listings): 127 have no carrier left, and 5 keep 8 of their carriers but, tied with
  another alignment of their row, leave its place to the one with more carriers (the same SNP
  beside an indel, or another SNP) and are among the alternative alignments; at chr22:35633518(+)
  the row's highest-CFD alignment is now chr22_35633524_C_T with chr22_35633527_ATG_A (448
  carriers), and the SNP alone keeps 4 of 2.8.1's 451. 64 alignments list 4,214 carriers fewer; 15
  SNP alignments take the place of ones gone (another alignment of their row, or a row of their own
  where the one gone had merged two targets into one row), and the search returns 20,726 rows. Its
  SNP post-analysis took 18 min 45 s, the same search's with 2.8.1, run beside it, 17 min 2 s.
  Pre-existing, in 2.8.0 too. Covered by `PostProcess/test_haplotype_prioritizer.py`.
- **A `--per-sample` search of an index without the dict-less tiers no longer reports a dense
  window by its greedy worst-case alignment.** On the classic dictionaries a window denser than
  `CRISPRME_IUPAC_CAP` was reported by an alignment that took every allele matching the guide and
  listed the carriers of any of them, so that nobody need carry it: on the self-test's
  `NGG_2_hg38_chr22+hg38_1000G` the README guide (`20bp-NGG-SpCas9`, 4 mismatches + 1 DNA + 1 RNA
  bulge) gave 5 such rows (chr22:22410074, 22687072, 22687116, 22698190 and 22888934), naming 102
  individuals as carriers 114 times. Such a window now gives the haplotypes its carriers have, read
  from the dictionaries as from the tiers: the same search returns 560 rows against 565, the five
  gone. Three other rows differ: chr22:44864861 lists one carrier fewer (the fix above),
  chr22:42848613 one other motif fewer (a greedy alignment of its cluster), and chr22:50763084's
  CRISPR-Bulge score its last digit (3.61989e-21, 3.6199e-21: a score of another batch). Neither
  `complete-search` nor the report says any longer that a dense window can be such an alignment
  (but the report of a 2.8.1 search's results, read again), and neither does a per-sample row's
  `High_variant_density_region` note (see Changed). Covered by
  `PostProcess/test_haplotype_prioritizer.py`, `PostProcess/test_phased_haplotype.py` and
  `PostProcess/test_generate_report.py`.
