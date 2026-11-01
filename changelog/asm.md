# [2.8.1] entries from the asm track (personal-assembly search)

### Fixed

- **A personal-assembly search with bulges searched the maternal haplotype with the paternal
  haplotype's index, and a second individual with the first one's.** `assembly-search` runs one
  `complete-search` per haplotype, and `complete-search` names a genome, and the index it builds or
  reuses (`genome_library/<PAM>_<N>_<name>`), after the genome folder's name. A bundle assembly
  (fetched from HPRC, or imported in Settings) keeps each haplotype in
  `Assemblies/<individual>/<haplotype>/genome`, so every haplotype of every individual was searched
  as `genome`: the paternal search built `NRG_2_genome`, the maternal search logged `Reference index
  present (NRG_2_genome)` and searched the paternal sequence, and the reconciliation then reported
  each site as found on both haplotypes. Haplotype-private off-targets were missed without a warning.
  Two `--genome-paternal` / `--genome-maternal` folders with the same name collided the same way.
  Reproduced with the real CLI on a synthetic bundle (4 mismatches + 1 DNA + 1 RNA bulge, the web
  default): the two haplotypes' result files were byte-identical, and at chr1:35000 the maternal
  results showed the paternal 2-mismatch protospacer instead of the maternal 4-mismatch one; a second
  individual whose paternal haplotype carries that 4-mismatch sequence also got the 2-mismatch one.
  Each haplotype is now searched through a link named after it,
  `haplotype_genomes/<key>/<individual>_<haplotype>` in the working directory, so it builds and reuses
  its own index (`NRG_2_HG01255_paternal`, `NRG_2_HG01255_maternal`); a genome folder that is not in
  a bundle is named after the folder plus a short hash of its path. The per-haplotype result files
  carry that name too (`<guide>+NRG_HG01255_paternal+__4+2_integrated_results.tsv` instead of
  `<guide>+NRG_genome+__4+2_integrated_results.tsv`). Each index
  a haplotype search builds records which genome it holds (`.genome_source`, which CRISPRitz does not
  read), so an individual deleted and imported again with other content gets its index rebuilt
  instead of reused; an index without that record, built earlier or downloaded, is used as before. A
  haplotype result written by the colliding search is not reused by a later run with the same
  `--output`: it is searched again. With the fix the paternal results are unchanged byte for byte,
  the maternal haplotype reports its own sites, and a legacy `Genomes/<name>` pair (the layout before
  2.6.0) runs exactly as before. **Affected:** every `assembly-search` with bulges on a bundle
  assembly since 2.6.0, from the web or the CLI. The haplotype whose search built the shared index
  (usually the first paternal haplotype searched) has correct results; every other haplotype
  searched with that index does not. **What to do:** run those searches again. The index
  `genome_library/<PAM>_<N>_genome` they left (in Settings > Installed indexes, e.g. `NRG_2_genome`)
  is no longer used and can be deleted there. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A personal-assembly search without bulges returned no off-target for a bundle assembly, and then
  stopped with a traceback.** The 0-bulge search scanned `Genomes/<genome name>/` instead of the
  genome folder it was given. A bundle haplotype has no folder there, the search engine reported
  `ERROR! /DATA/Genomes/genome/ is not a directory.` but exited 0, and both haplotypes came back
  empty ("The search completed successfully; the result set is empty"); the combined report then
  failed on the empty set (`TypeError`), so `assembly-search` exited 1. `--bDNA` and `--bRNA` are
  optional and default to 0, so this was the documented minimal command
  (`assembly-search --assembly-individual <name> --guide ... --pam ... --mm ...`). The search now
  scans the folder it is given: on the synthetic bundle each haplotype reports its 9 sites and the
  run completes. Searches of a genome in `Genomes/` are unchanged (the same folder); reference and
  variant search outputs were compared byte for byte, with and without bulges. Covered by
  `PostProcess/test_assembly_haplotype_genomes.py`.
- **The personal-assembly results header and History named both haplotypes "genome".** They showed
  `Paternal genome: genome | Maternal genome: genome`: the name of the folder every bundle haplotype
  is kept in, which `assembly-search` wrote into the job's `.Params.txt` when it finished (the web
  form's own label was overwritten). Both now show the bundle's assembly name, e.g.
  `HG01255_paternal_<accession>` for an HPRC fetch, or `<individual>_<haplotype>` when the bundle
  records none, from the CLI and the web alike; a legacy `Genomes/<name>` folder keeps its name.
  Jobs that already finished keep the label they were written with.
- **Submitting the Personal assembly tab with no individual selected listed six fields the form does
  not show.** The warning listed "paternal genome", "maternal genome", the two chain files and the
  two chromAlias files, which the Individual picker fills. It now says "select an individual on the
  Personal assembly tab" once; any other missing input, such as the guide, is still listed.
- **The Personal assembly tab ignored the "Maximum edits" slider.** The web form never passed it on,
  so an assembly search always allowed mismatches + DNA + RNA bulges in total, whatever the slider
  showed (with the slider at 3 and Advanced options closed, the results header said "Max total
  edits: 6"). The slider now applies as on the Reference genome tab: with Advanced options closed it
  caps the total edits (`assembly-search --max-total-edits`), with them open the per-type caps
  govern. A search left at the defaults (slider 6 = 4 mismatches + 1 + 1 bulges) runs the same
  command as before.
