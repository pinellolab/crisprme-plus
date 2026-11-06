<!-- gate track: [2.9.0] entries for CHANGELOG.md (the carrier gate, plan 2.3; decisions D1, D2,
D3, D5, D6; the F2 and greedy-row fixes the 2.8.1 text promised for this release). -->

**For the merge, not for the changelog.** 2.8.1's `### Known issues` entries on SNP rows that
list a carrier whose haplotype holds an indel inside the target, and on the greedy rows of a
`--per-sample` search on an index without the dict-less tiers, both say "fixed in the next
release": the two `### Fixed` entries below are those fixes. The numbers measured on the 2.8.1
population-level run are of gate38a's arm C9 (the 2.8.1 code), re-gated with this release's
`haplotype_prioritizer.py`; the others are of this branch's own runs.

### Added

- **The carrier gate: a population-level row's variant combination is checked against phased
  haplotypes, and the report ranks a combination no haplotype carries after every other row,
  never dropping it.** A population-level row takes, per variant window, the alleles that match
  the guide best and lists every carrier of any of them, so two or more of them can sit on no one
  haplotype: no individual has the sequence the row reports, and, being the best-matching, it
  ranked first. After writing the integrated results the search now reads each row's variants
  (each of its three alignments) in the dataset's own phased genotypes, SNPs and indels
  (`Dictionaries/genotypes_<vcf>/`, `indel_genotypes_<vcf>/`;
  `PostProcess/haplotype_prioritizer.py`), and gives each a tier: **A**, two or more variants one
  phased haplotype carries together (states `one haplotype`, and `phase not known` where only
  unphased calls, or two ALTs of one position on one slot, show it); **B**, one variant; **C1**,
  two or more variants, each in the panel, that none of its haplotypes carries together, with E,
  the haplotypes expected to carry them together if they were independent (below 3, the state
  is `too rare to tell`), or variants every haplotype carrying which has an indel that changes
  the target (`indel in the target`); **C2**, a variant not in the panel, which cannot say (the
  frequency bound, the product of the variants' frequencies >= 1e-06 and the rarest >= 0.001, is
  shown and applied to nothing). The gate changes no row: it writes `..._carrier_gate.tsv` beside
  the integrated results (one line per alignment that uses a variant) and `.carrier_gate.json`.
  The report leaves a row whose highest-CFD alignment is C1 out of the top-1000 tables, the plots
  of Section 2, the validation panel and the risk card, and lists it in a table of its own after
  the top-1000 (`unobserved_worst_case.tsv` holds them all); a row whose highest-CRISPR-Bulge
  alignment is C1 leaves the ranking by that score. A, B and C2 rows keep their place. The
  `Observed` column of a row of two or more variants gives its tier as a factual label
  (`together on N of H <panel> haplotypes`, `not together on any of H <panel> haplotypes
  (expected E)`, `not checked: K of N variants not in <panel>`), with no "chimera" word, and the
  inputs card a **Carrier gate** row: its panels, its counts by tier, and the panel's samples by
  superpopulation (AFR 986, AMR 549, CSA 782, EAS 816, EUR 771, MID 157, OCE 30 on HGDP + 1000
  Genomes), as a haplotype found only in an ancestry with few samples there can be missing from
  it. On by default for a population-level search of a genotyped dataset; `--no-prioritize-observed`
  turns it off, `--gate-min-expected` sets the `too rare to tell` threshold (default 3; 0 turns it
  off). A per-sample search, whose rows name their cis carriers, needs none. On the 2.8.1
  genome-wide population-level search of the README guide (HGDP + 1000 Genomes, NRG, 4
  mismatches + 1 DNA + 1 RNA bulge, 198,826 rows) it took 26 s on 24 threads (1.3 GB), left the
  integrated results' md5 as it was, and gave, by each row's highest-CFD alignment, A 1,581, B
  55,183, C1 42,813 (278 where 3 or more haplotypes were expected, 42,323 `too rare to tell`, 212
  `indel in the target`) and C2 0. Checked against all 8,182 phased haplotypes by a harness that
  rebuilds every alignment's window on each haplotype from the VCF: no alignment left in tier A
  or B lists a combination no haplotype carries (each of the three alignments), and every C1
  alignment that a haplotype does carry is one with an indel in the target. Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **`--prioritize-observed` checks a sites-only dataset's rows against the phased HGDP + 1000
  Genomes genotypes, and `download --what genotype-oracle` installs them.** The `mega` panel has
  no genotypes, so its multi-variant rows could not be told apart. With `--prioritize-observed`
  (the web form's check box, offered for a sites-only panel once they are installed) the gate
  reads its rows in the phased HGDP + 1000 Genomes genotypes and indel genotypes (the genotype
  oracle), which an installed `NRG_3_hg38+hg38_HGDP1kGP` index already has, which
  `download --what genotype-oracle [--panel HGDP1kGP]` installs without the index, and which
  `CRISPRME_GENOTYPE_TIER` can name elsewhere. A variant the panel does not hold makes its row
  C2, kept in place. Without the oracle the search says so, and the rows are reported unchecked.
  `publish-index` now also uploads a variant index's phased indel genotypes on their own
  (`indel_genotypes_<vcf>.tar.gz`), for the oracle. Covered by
  `PostProcess/test_haplotype_prioritizer.py` and `PostProcess/test_web_search_reuse.py`.

### Fixed

- **A per-sample SNP row no longer names a carrier whose haplotype holding its SNPs also holds an
  indel inside the target.** Its carriers were read from the SNP genotypes alone, so an
  individual whose only haplotypes holding the row's SNPs carried an indel changing the bases the
  target covers was listed, though its own sequence there is the indel's (2.2% of SNP row-blocks
  on chr22 of the HGDP + 1000 Genomes panel, five guides; 545 of them had no carrier holding the
  sequence). The per-sample post-analysis now reads the index's phased indel genotypes, on the
  dict-less tiers and on the classic dictionaries alike, and takes such a carrier off; a row left
  with no carrier is not reported. An indel at one of the row's SNPs, or over one, is another
  allele of that base, so a slot phased with both (a phasing artefact) keeps its carrier, as an
  unphased call of the SNP or of the indel does. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **A `--per-sample` search of an index without the dict-less tiers no longer reports a dense
  window by its greedy worst-case alignment.** On the classic dictionaries a window denser than
  `CRISPRME_IUPAC_CAP` was reported by an alignment that took every allele matching the guide and
  listed the carriers of any of them, so that nobody need carry it: on the self-test's
  `NGG_2_hg38_chr22+hg38_1000G` the README guide (`20bp-NGG-SpCas9`, 4 mismatches + 1 DNA + 1 RNA
  bulge) gave 5 such rows (chr22:22410074, 22687072, 22687116, 22698190 and 22888934), naming 102
  individuals as carriers 114 times. Such a window now gives the haplotypes its carriers have,
  read from the dictionaries as from the tiers: the same search returns 560 rows against 565, the
  five gone. The report and `complete-search` no longer say that a dense window can be such an
  alignment, but of a 2.8.1 search's results, read again. Covered by
  `PostProcess/test_haplotype_prioritizer.py`, `PostProcess/test_phased_haplotype.py` and
  `PostProcess/test_generate_report.py`.
