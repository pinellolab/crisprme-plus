<!-- density track: [2.9.0] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the merge, not for the changelog: three notes.**

1. The entries name `Q`, which the variant genome writes at a four-allele position from CRISPRitz
   2.10.0 (the engine track), and the dense threshold, 2^20 sequences, which the fifth-symbol
   measurement may tune (`window_density.DENSE_WINDOW_BITS`, the one constant): keep both as they
   land, here and in METHODS.md §5.
2. The chr22 counts come from the README guide's searches of chr22 of HGDP + 1000 Genomes in the
   2.8.1 image (rc281g) with v2.8.1's tree and with this track's (`NRG_3_hg38+hg38_HGDP1kGP`, 4
   mismatches + 1 DNA + 1 RNA bulge, on ml003); the index rebuilt with every allele changes them.
   The known limit the Changed entry names (a dense window's greedy alignment on an index without
   the dict-less tiers, under `--per-sample`) is to be fixed in this release by the carrier-gate
   track; if it is, drop that clause.
3. With CRISPRitz 2.10.0 (the engine track's 25ee50c, every allele coded, `Q` at all 2,448
   four-allele positions of chr22, uncapped), the search of the five guides on chr22 of HGDP + 1000
   Genomes (NRG, 4 mismatches + 1 DNA + 1 RNA bulge) reports 4,372,248 rows where the 2.9.1 genome
   gave 206,883, and 4,201,268 of them (96.1%; 375,834 of its 462,663 windows) lie in a window of
   2^20 sequences or more (the worst 2^46.34, a bulged window). The dense flag thus marks nearly
   all of what the fifth symbol adds, and a post-analysis would have that many rows to resolve:
   the "stricter density cap" of the fifth-symbol decision (an enricher-side cap on `Q`, or the
   engine's `--max-iupac`) is what bounds them; this track's cap and `--max-degeneracy` bound what
   one window gives, not how many windows there are. Measured with this track's
   `window_density.py` on the engine track's `idx_mm4b11.targets.txt`, on ml008. This CRISPRme
   (2.8.1's reader) refuses the 2.10.0 indexes (format 2), so no CRISPRme run of them is here.

### Added

- **A window that spells 2^20 sequences or more is flagged dense, whatever `CRISPRME_IUPAC_CAP`
  is.** A window's sequences are the product of the bases each letter of the variant genome stands
  for: 2 for a two-allele column (R, Y, S, W, K, M), 3 for a three-allele one (B, D, H, V), 4 for a
  four-allele one (`Q`), 1 for a base or an `N` (an assembly gap, which matches nothing). Such a
  window matches nearly any guide somewhere, as uninformative as a gap, yet the search reports it,
  each of its columns admitting the guide's base. It is listed in
  `<out>.high_variant_density_regions.bed`, which gains two columns, `bits` (log2 of the window's
  sequences) and `four_allele` (how many four-allele positions it holds), and the report's
  `High_complexity_region` column gives its count (`24 in window (2^24.0 sequences)`), as the
  integrated results' `High_variant_density_region` note does (`2^24.00 sequences`). The report's
  inputs card says how many rows lie in such a window, and `complete-search` says per chromosome
  how many windows were over `CRISPRME_IUPAC_CAP` and how many dense (stdout, `log_verbose.txt`).
  Among the indexable NRG windows of chr22, the reference genome has none, the variant genome of
  HGDP + 1000 Genomes published with CRISPRitz 2.9.1 has 1,223 (the worst 2^23.00), and one coding
  every allele of a split site 3,780 (the worst 2^36.45, at chr22:22,888,423 in the IGL locus). On
  the published index, the README guide's search of chr22 lists 565 windows over the cap, 381 of
  them dense (the worst 2^24.00), and 164 of the population-level report's 2,198 rows lie in one
  (1 of the per-sample report's 1,574). The threshold is one constant,
  `window_density.DENSE_WINDOW_BITS`; see
  [`METHODS.md` §5](METHODS.md#5-search-space-control-for-high-variant-density-regions).
- **A row whose window holds a four-allele position says so.** Where the reference base and three
  ALTs cover A, C, G and T, the variant genome writes `Q`, which matches any base (`N` stays an
  assembly gap). The post-analysis lists the four-allele positions of the windows the search
  reported in `<out>.four_allele_positions.bed` (0-based; written when there is one) and says how
  many on stdout; `High_variant_density_region` names those in a row's window
  (`four_allele_positions=<chrom>:<position>`, 1-based), the report's `High_complexity_region`
  column reads `1 four-allele position`, and its inputs card says how many rows have one. A window
  that holds a four-allele position and is not dense is no dense window to the report.
- **`complete-search --max-degeneracy K`: the most ambiguous columns a reported SNP row may rest
  on, off by default.** An ambiguous column a row rests on is one of the variant genome where the
  reference base would be a mismatch, or break the PAM, and the row's allele fits. Few alignments
  resting on two or more are carried by anyone: on chr22 of HGDP + 1000 Genomes (the published
  variant genome; five guides, NRG, 4 mismatches + 1 DNA + 1 RNA bulge) a real phased haplotype
  carries 96.8% of the alignments resting on one, 5.6% of those resting on two, 1.1% on three and
  0.02% on four or more. With the cap, a worst-case alignment (population-level, or a dense
  window's greedy one) is cut back to the strongest within it (a valid PAM first, then the fewest
  mismatches or the highest CFD), and any other row resting on more is left out, unless its
  carriers hold its variants together on one phased haplotype (a CONFIRMED observed haplotype, or a
  phased dictionary's combination); the reference row keeps the locus. `complete-search` says per
  chromosome how many alignments were cut back and how many rows left out, writes the cap to
  `Params.txt` (`Max_degeneracy`), and the report's search criteria state it. An indel row is
  resolved per carrier haplotype, as before, and is not cut. On chr22 of the published index, the
  README guide's population-level search with `--max-degeneracy 1` cut 20,918 worst-case alignments
  back and reports 1,573 rows instead of 2,198, every SNP variant row then resting on at most one
  column (of the 1,257 such rows without it, 681 rested on two or more): the 634 sites only such
  alignments reach go, and 49 of the 100 sites with the highest CFD leave the top 100 (the second
  and third, CFD 0.299 at chr22:22358032 in IGL and 0.275 at chr22:49464611, rest on 11 and 2
  columns). The per-sample search changes no row: its multi-variant rows are CONFIRMED haplotypes.

### Changed

- **`CRISPRME_IUPAC_CAP` counts bits, log2 of the sequences a window spells, where it counted the
  window's ambiguous columns.** Above it the post-analysis of an index without the dict-less tiers
  reports a window by its greedy worst-case alignment instead of every combination of its alleles,
  and every path lists the window in `high_variant_density_regions.bed`. The default, 10, caps a
  window of two-allele columns where the count of columns did: on the published HGDP + 1000 Genomes
  index, whose chr22 holds no three-allele column, the README guide's searches of chr22 report the
  same rows as with 2.8.1 (per-sample 1,574, population-level 2,198; their `bestMerge` tables are
  byte-identical, and their BEDs list the same 565 windows, 8,155 lines), the notes only giving the
  windows' sequences. A three-allele column weighs 1.58 bits and a four-allele one 2, where each
  counted 1: coding every allele of a split site makes them common, and the immunoglobulin loci
  explode (on a chr20 + chr22 install of HGDP + 1000 Genomes coded so, the variant search's raw
  targets went from 436,393 to 2,001,729, 81% of them in IGL), which the count of columns could not
  see. A VCF line that lists several ALTs writes a three-allele column with any release: on a small
  genome whose phased VCF has seven such lines in one window of the README guide (3^7 = 2^11.09
  sequences), 2.8.1 listed the window nowhere; it is listed now, and an index built without
  `--samplesID` reports it by its greedy worst-case alignment under `--per-sample` too, with the
  known limit of 2.8.1 for such an index: its listed carriers need not carry that alignment
  (METHODS.md §5). `-1` still disables the cap.
