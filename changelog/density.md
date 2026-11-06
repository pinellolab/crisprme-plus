<!-- density track: [2.9.0] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the merge, not for the changelog: three notes.**

1. The entries name `Q`, which the variant genome writes at a four-allele position from CRISPRitz
   2.10.0 (the engine track), and the hyperdense threshold, 2^20 sequences
   (`window_density.HYPERDENSE_BITS`, the one constant), which the fifth-symbol measurement
   (`_r281/MEASURE_fifth_symbol.md`) keeps at 20 bits: keep both as they land, here and in
   METHODS.md §5. "Dense" keeps the meaning 2.8.1's report, README, METHODS and help give it, a
   window over `CRISPRME_IUPAC_CAP`; "hyperdense" is one of 2^20 sequences or more.
2. The chr22 counts come from the README guide's searches of chr22 of HGDP + 1000 Genomes in the
   2.8.1 image (rc281g) with v2.8.1's tree and with this track's (`NRG_3_hg38+hg38_HGDP1kGP`, 4
   mismatches + 1 DNA + 1 RNA bulge, on ml003); the index rebuilt with every allele changes them.
   The known limit the Changed entry names (a dense window's greedy alignment on an index without
   the dict-less tiers, under `--per-sample`) is to be fixed in this release by the carrier-gate
   track; if it is, drop that clause. The chr22 window census (NRG windows of 2^20 sequences or
   more: reference 0, 2.9.1 1,763, 2.10.0 5,103) is the fifth-symbol measurement's (its §3),
   re-counted with this track's `window_density.py` on the same genomes on ml008: the published
   2.9.1 chr22 (`crisprme_canonical_2.8.0`) is byte-identical to the one measured, and the 2.10.0
   one is the engine track's 25ee50c enrichment. The plan's 1,223 and 3,780 (its §2.3) are another
   count, not reproduced here (on either strand 1,763 and 4,745, on the plus strand only 1,123 and
   3,234), and 3,780 is the four-base rule's genome, which no release writes.
3. With CRISPRitz 2.10.0 (the engine track's 25ee50c, every allele coded, `Q` at all 2,448
   four-allele positions of chr22, uncapped), the search of the five guides on chr22 of HGDP + 1000
   Genomes (NRG, 4 mismatches + 1 DNA + 1 RNA bulge) reports 4,372,248 rows where the 2.9.1 genome
   gave 206,883, and 4,201,268 of them (96.1%; 375,834 of its 462,663 windows) lie in a window of
   2^20 sequences or more (the worst 2^46.34, a bulged window): the hyperdense flag marks nearly
   all of what the fifth symbol adds. Measured with this track's `window_density.py` on the engine
   track's `idx_mm4b11.targets.txt`, on ml008. No cap upstream of CRISPRme will bound those rows:
   the fifth-symbol measurement (its §8) rejects an engine-side `--max-iupac`, which cannot tell a
   reference base from an ALT (at 20 bits it loses 21 IGL and 47 IGH loci a real haplotype
   carries, 14 and 28 of them loci 2.9.1 reports, so the superset gate G5 would fail), and keeps
   the bound in CRISPRme at 20 bits (its §7); the four-allele judge (`r109/fourallele_design.json`,
   point 6) chose to flag and gate in CRISPRme, not to cap in the enricher. So the "stricter density
   cap" of the fifth-symbol decision is this track's: `CRISPRME_IUPAC_CAP` counts bits, a `Q`
   weighing 2, and at its default of 10 every hyperdense window is over it, so the dictionary path,
   and any population-level search, reports each such window by its greedy worst-case alignment;
   per sample the tier path enumerates its observed haplotypes, as for any window; and
   `--max-degeneracy` cuts a worst case back. Two points stay open for the merge with the format
   and readers tracks: what the post-analysis costs per raw row on a 2.10.0 index, not measured here
   (this CRISPRme, 2.8.1's reader, refuses the 2.10.0 indexes, format 2), and the measurement's
   advice that the dense path read the observed haplotypes once per window, not once per row,
   which is not built here.

### Added

- **A window that spells 2^20 sequences or more is flagged hyperdense, whatever `CRISPRME_IUPAC_CAP`
  is.** A window's sequences are the product of the bases each letter of the variant genome stands
  for: 2 for a two-allele column (R, Y, S, W, K, M), 3 for a three-allele one (B, D, H, V), 4 for a
  four-allele one (`Q`), 1 for a base or an `N` (an assembly gap, which matches nothing). Such a
  window matches nearly any guide somewhere, as uninformative as a gap, yet the search reports it,
  each of its columns admitting the guide's base. It is listed in
  `<out>.high_variant_density_regions.bed`, read from its own letters (a window whose only variant
  letters are four-allele positions too), and the report's `High_complexity_region` column gives its
  count (`24 in window (2^24.0 sequences)`), as the integrated results'
  `High_variant_density_region` note does (`2^24.00 sequences`). The BED gains three columns, `bits`
  (log2 of the window's sequences), `four_allele` (how many four-allele positions it holds) and
  `capped` (1 when the window was over `CRISPRME_IUPAC_CAP` or the search population-level), and its
  `n_variants` counts every variant column of the window, a four-allele one too. Only a capped
  window's note says that a greedy worst-case alignment is reported there: within a cap raised to
  its bits, or with `CRISPRME_IUPAC_CAP=-1`, a `--per-sample` search enumerates a hyperdense
  window's combinations as any other's, and its rows read as in 2.8.1 (`N cis carriers`). The
  report's inputs card says how many rows lie in such a window, and `complete-search` says per
  chromosome how many windows were dense (over `CRISPRME_IUPAC_CAP`) and how many hyperdense
  (stdout, `log_verbose.txt`). Among the indexable NRG windows of chr22, the reference genome has
  none, the variant genome of HGDP + 1000 Genomes published with CRISPRitz 2.9.1 has 1,763 (the
  worst 2^23.00), and the one CRISPRitz 2.10.0 writes from the same VCF, every allele coded, 5,103
  (the worst 2^44.75, at chr22:22,881,189 in the IGL locus). On the published index, the README
  guide's search of chr22 lists 565 dense windows, 381 of them hyperdense (the worst 2^24.00), and
  164 of the population-level report's 2,198 rows lie in one (1 of the per-sample report's 1,574).
  The threshold is one constant, `window_density.HYPERDENSE_BITS`; see
  [`METHODS.md` §5](METHODS.md#5-search-space-control-for-high-variant-density-regions).
- **A row whose window holds a four-allele position says so.** Where the reference base and three
  ALTs cover A, C, G and T, the variant genome writes `Q`, which matches any base (`N` stays an
  assembly gap). The post-analysis lists the four-allele positions of the windows the search
  reported in `<out>.four_allele_positions.bed` (0-based; written when there is one) and says how
  many on stdout; `High_variant_density_region` names those in a row's window
  (`four_allele_positions=<chrom>:<position>`, 1-based), the report's `High_complexity_region`
  column reads `1 four-allele position`, and its inputs card says how many rows have one. A
  four-allele position alone makes no window dense to the report.
- **`complete-search --max-degeneracy K`, off by default: a reported SNP row rests on at most K
  ambiguous columns unless its carriers hold its variants together on one phased haplotype.** An
  ambiguous column a row rests on is one of the variant genome where the reference base would be a
  mismatch, or break the PAM, and the row's allele fits. Few alignments resting on two or more are
  carried by anyone: on chr22 of HGDP + 1000 Genomes (the published variant genome; five guides,
  NRG, 4 mismatches + 1 DNA + 1 RNA bulge) a real phased haplotype carries 96.8% of the alignments
  resting on one, 5.6% of those resting on two, 1.1% on three and 0.02% on four or more. With the
  cap, a worst-case alignment (population-level, or a dense window's greedy one) is cut back to the
  strongest within it (a valid PAM first, then the fewest mismatches or the highest CFD), and any
  other row resting on more is left out, unless its carriers hold its variants together that way (a
  CONFIRMED observed haplotype, or a phased dictionary's combination). The reference row keeps the
  locus where the reference alignment is within the search's limits; a site that only such
  alignments reach goes. `complete-search` says per chromosome how many alignments were cut back and
  how many rows left out, writes the cap to `Params.txt` (`Max_degeneracy`), and the report's search
  criteria state it. An indel row is resolved per carrier haplotype, as before, and is not cut. On
  chr22 of the published index, the README guide's population-level search with `--max-degeneracy 1`
  cut 20,918 worst-case alignments back and reports 1,573 rows instead of 2,198, every SNP variant
  row then resting on at most one column (of the 1,257 such rows without it, 681 rested on two or
  more): 625 sites go, those only such alignments reach (9 more move by 1 or 2 nt, to where their
  cut-back alignment starts), the rows lying in a dense window fall from 281 to 2, and 47 to 49 of
  the 100 sites with the highest CFD leave the top 100, as the ties at the 100th CFD fall (the
  second and third, CFD 0.299 at chr22:22358032 in IGL and 0.275 at chr22:49464611, rest on 11 and 2
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
  explode (on chr22 of HGDP + 1000 Genomes, five guides' variant search returns 4,372,248 raw
  targets with CRISPRitz 2.10.0, 96.6% of them in IGL, where 2.9.1 returned 206,883), which the
  count of columns could not see. A VCF line that lists several ALTs writes a three-allele column
  with any release: on a small genome whose phased VCF has seven such lines in one window of the
  README guide (3^7 = 2^11.09 sequences), 2.8.1 listed the window nowhere; it is listed now, and an
  index built without `--samplesID` reports it by its greedy worst-case alignment under
  `--per-sample` too, with the known limit of 2.8.1 for such an index: its listed carriers need not
  carry that alignment (METHODS.md §5). `-1` still disables the cap; a hyperdense window is listed
  all the same, not capped.
