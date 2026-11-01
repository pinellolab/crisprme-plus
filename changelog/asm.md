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
  each site as found on both haplotypes. Haplotype-private off-targets were missed without a
  warning. Two `--genome-paternal` / `--genome-maternal` folders with the same name collided the
  same way. Reproduced with the real CLI on a synthetic bundle (4 mismatches + 1 DNA + 1 RNA bulge,
  the web default): the two haplotypes' result files were byte-identical, and at chr1:35000 the
  maternal results showed the paternal 2-mismatch protospacer instead of the maternal 4-mismatch
  one; a second individual whose paternal haplotype carries that 4-mismatch sequence also got the
  2-mismatch one. Each haplotype is now searched through a link named after it and its folder,
  `haplotype_genomes/<key>/<individual>_<haplotype>_<hash>` in the working directory, where `<hash>`
  is 8 hex digits of the folder's path, so it builds and reuses its own index
  (`NRG_2_HG01255_paternal_5aeb0fc6`, `NRG_2_HG01255_maternal_067bdb25`); a genome folder that is
  not in a bundle is named after the folder plus the same hash. The hash keeps these names apart
  from every genome in `Genomes/`: a bundle search never takes the index of a
  `Genomes/<individual>_<haplotype>` folder of the same individual (the layout before 2.6.0, or a
  custom upload), a search of that folder never takes the bundle's, and two bundles whose folder
  names differ only in characters a file name cannot keep (`HG.002`, `HG_002`) get an index each.
  The per-haplotype result files carry that name too
  (`<guide>+NRG_HG01255_paternal_5aeb0fc6+__4+2_integrated_results.tsv` instead of
  `<guide>+NRG_genome+__4+2_integrated_results.tsv`). Each index a haplotype search builds records
  the genome folder it was built from and the names, sizes and modification times of that folder's
  FASTA files (`.genome_source`, which CRISPRitz does not read), and a later search of that
  haplotype rebuilds the index when they have changed: a FASTA file added, removed, resized or
  modified at another time, as when an individual is deleted and imported again with other files. A
  FASTA file replaced by other content of the same size and modification time is not noticed, nor is
  a change made after a search that was killed before it ended (its container stopped, say): the
  record is written when the search ends, and an index without one is used as it is, as an index of
  a legacy `Genomes/<name>` folder built by an earlier CRISPRme, or downloaded, is used as before. A
  haplotype result written by the colliding search is not reused by a later run with the same
  `--output`: it is searched again. With the fix the paternal results are unchanged byte for byte
  apart from the genome's name (in file names, `.Params.txt` and `.command_line.txt`), the maternal
  haplotype reports its own sites, and a legacy `Genomes/<name>` pair (the layout before 2.6.0) runs
  exactly as before. Each haplotype's index is the size of a reference index; Settings > Delete >
  Personal assembly removes the individual's haplotype indexes and links with its folder, and its
  confirmation names those indexes and counts their size. **Affected:** every `assembly-search` with
  bulges on a bundle assembly since 2.6.0, from the web or the CLI. The haplotype whose search built
  the shared index (usually the first paternal haplotype searched) has correct results; every other
  haplotype searched with that index does not. **What to do:** run those searches again, and delete
  the index `genome_library/<PAM>_<N>_genome` they left (e.g. `NRG_2_genome`: Settings > Delete >
  `Index: NRG_2_genome`). `assembly-search` no longer uses it, but a `complete-search` with bulges
  of a folder named `genome` still finds it by that name: such a search of one haplotype of a bundle
  is refused while it is there (see below), and a search of any other folder named `genome` takes
  it. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search without bulges of a genome folder outside `Genomes/` came back empty, or scanned
  another genome.** The 0-bulge search reads the genome's FASTA files directly, and read them from
  `Genomes/<name of the --genome folder>/` instead of the folder it was given. Without such a folder
  the search engine reported `ERROR! /DATA/Genomes/<name>/ is not a directory.` but exited 0, and
  the search ended with "The search completed successfully; the result set is empty". With one, that
  folder was scanned in place of the given one (a `--genome other/hg38` whose chr1 carries a
  2-mismatch site at chr1:35000 reported the 4-mismatch protospacer of `Genomes/hg38` there). On the
  command line this hit `complete-search --genome <folder outside Genomes/>` with `--bDNA 0 --bRNA
  0`, and the same search without bulge flags when the library has no index of that genome, since
  the bulges then default to 0 (`complete-search --genome /data/mygenome --guide ... --pam ... --mm
  4`). It also hit every personal-assembly search without bulges whose haplotypes are kept outside
  `Genomes/`, as every bundle's are (`Assemblies/<individual>/<haplotype>/genome`): both haplotypes
  came back empty, and the combined report then failed on the empty set (`TypeError`, see below), so
  `assembly-search` exited 1; `--bDNA` and `--bRNA` are optional and default to 0, so this was the
  documented minimal command (`assembly-search --assembly-individual <name> --guide ... --pam ...
  --mm ...`). The search now scans the folder it is given: a copy of `Genomes/hg38` outside
  `Genomes/` reports the same 9 sites as `Genomes/hg38`, and on the synthetic bundle each haplotype
  reports its 9 sites. A search of a genome in `Genomes/`, which includes every search from the
  web's Reference genome tab, scans the same folder as before; reference and variant search outputs
  were compared byte for byte, with and without bulges. A search of such a folder with bulges takes
  the index of its name, and is now refused when that index, or a folder of that name in `Genomes/`,
  is another folder's (next entry). **Affected:** 0-bulge `complete-search` runs of a genome folder
  outside `Genomes/` from the command line, and 0-bulge personal-assembly searches of haplotypes
  outside `Genomes/` (every bundle assembly) from the CLI or the web. **What to do:** run them
  again. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search with bulges of a genome folder outside `Genomes/` could search another folder's
  sequence.** `complete-search` takes a genome's index (`genome_library/<PAM>_<N>_<name>`) and its
  variant genomes (`Genomes/<name>+<dataset>`) by the genome folder's name alone, whatever folder
  they were built from. With bulges, `--genome /data/projB/ref` took the index a search of
  `/data/projA/ref` had built and reported projA's sites (at chr1:35000 of the synthetic bundle,
  projA's 2-mismatch protospacer instead of projB's own 4-mismatch one), and `--genome other/hg38`
  took the index of `Genomes/hg38` and reported its sites, as the same command without bulges did by
  scanning `Genomes/hg38` itself (previous entry). Given no `--bDNA`/`--bRNA`, that search takes 1 +
  1 bulges when `Genomes/hg38` is there, so the plain command was hit too. The other way round, a
  search of `Genomes/<name>` took the index that a search of a folder of that name outside
  `Genomes/` had built. `complete-search` now refuses a search with bulges, or with `--vcf`, of a
  folder outside `Genomes/` whose name a folder in `Genomes/` has, and of any genome folder whose
  name an index recorded from another folder carries, with a message naming the folder or index in
  the way. The message says how to run the search again: with its `--bDNA`/`--bRNA` (given none, the
  search of a renamed folder, or of one whose index was removed, may take no bulges and miss every
  site that needs one), and, under `--index-path`, without it, since the pipeline builds no index
  there (and so without `--prune-pam`, which needs `--index-path`). An index that `complete-search`
  builds of a folder outside `Genomes/` records that folder and its FASTA files (`.genome_source`,
  as a personal-assembly haplotype's index does), and is built again when those files change. The
  record names a folder of the machine that built the index, so `publish-index` leaves it out of the
  archive and `download` drops it: a downloaded index records no folder. A search of a genome in
  `Genomes/` runs exactly as before and its index records nothing; an index that records nothing
  (built by an earlier CRISPRme, or downloaded) is taken as before, except by a search of one
  haplotype of a personal assembly (next entry), and so are the indexes of a library given with
  `--index-path` other than the working directory's `genome_library/`. The reference and variant
  search outputs of `Genomes/hg38` were compared byte for byte, with and without bulges. Two cases
  are not caught: two first searches of two such folders that start together (one of them searches
  the index the other builds), and a search from the web's Reference genome tab, which runs the
  pipeline without this check, of a `Genomes/<name>` added after a folder of that name outside
  `Genomes/` was searched with bulges. **Affected:** `complete-search` runs from the command line,
  with bulges or with `--vcf`, of a folder outside `Genomes/` named like a folder in `Genomes/` or
  like another folder searched with bulges in the same working directory, and searches, from the CLI
  or the web, of that `Genomes/` folder made after them. **What to do:** run them again, with a name
  of its own for each folder and with the bulges they ran with (`--bDNA 1 --bRNA 1` for one with
  bulges that gave none). An index built before 2.8.1 records no folder: where folders of one name
  outside `Genomes/` were searched with bulges, remove the indexes of that name (Settings > Delete >
  `Index: <PAM>_<N>_<name>`, or the folder `genome_library/<PAM>_<N>_<name>`), so that each search
  builds its own. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search with bulges of one haplotype of a personal assembly could search another haplotype's
  sequence.** Searched on its own, `complete-search --genome
  Assemblies/<individual>/<haplotype>/genome` names the genome, and its index, after that folder,
  `genome` for every haplotype of every bundle (first entry), so with bulges it took any
  `genome_library/<PAM>_<N>_genome` there: the one an `assembly-search` with bulges had left, or one
  that a search like this one of another haplotype had built. Given no `--bDNA`/`--bRNA`, the search
  takes 1 + 1 bulges when such an index is there, so the plain command was hit too. Reproduced on
  the synthetic bundle after an `assembly-search` of it: the maternal haplotype logged `Reference
  index present (NRG_2_genome)` and reported the paternal 2-mismatch protospacer at chr1:35000
  instead of its own 4-mismatch one. Such a search is now refused before it starts while an index of
  that name that it would take records no genome folder, as one built before 2.8.1 does, or while
  one records another folder (as any search is: previous entry). The message names the index and
  says to search the individual with `assembly-search`, which searches each haplotype on an index of
  its own, or to remove the index and run the search again, which then builds its own from the
  haplotype's folder. For either it gives the search's `--bDNA`/`--bRNA`: `assembly-search` takes no
  bulges by default, and, given none, the search of the haplotype may take none once the index is
  removed (on the synthetic bundle it then reported 9 sites instead of 17). A search without bulges
  and without `--vcf` is not refused: it scans the haplotype's folder (second entry).
  `assembly-search` and searches of any other folder are not refused this way. **Affected:**
  `complete-search` runs with bulges, from the command line, of one haplotype folder of a bundle, in
  a data folder where another haplotype had been searched with bulges. **What to do:** run them
  again with `assembly-search`, or delete the index (e.g. `NRG_2_genome`) and run them again, giving
  the bulges they ran with (`--bDNA 1 --bRNA 1` for one that gave none; `assembly-search`'s default
  is 0). Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A genome in `Genomes/` that is a link to a folder of another name failed in the post-analysis.**
  The pipeline resolved the link and named the genome after the link's target in one step and after
  the link in another, so a search of `Genomes/hg38lnk` -> `../shared/GRCh38_copy` stopped with
  `grep: .../crispritz_targets/GRCh38_copy_..._targets.txt: No such file or directory` and `ERROR:
  off-targets post-analysis (reference) failed` (exit 1), with or without bulges. The web runs the
  same pipeline on `Genomes/<name>`, so its searches of such a genome failed the same way. A link
  named like its target (`Genomes/hg38` -> `/shared/genomes/hg38`) was not affected. The genome now
  keeps the name of the folder the search is given (each personal-assembly haplotype is searched
  through such a link, named after the haplotype). The linked genome's results equal those of the
  same sequence in `Genomes/hg38`, apart from the genome's name in the file names. **Affected:**
  every reference search (no `--vcf`) of such a link, from the CLI or the web. **What to do:** run
  them again. A variant-aware search (`--vcf`) of such a link still fails, as before, at the
  enrichment (`mv: cannot stat '.../SNPs_genome/hg38lnk_enriched/'`, then `ERROR: Genome enrichment
  failed on <dataset>`): for it, name the link like its target folder, or put the folder itself in
  `Genomes/`.
- **A personal-assembly search stopped before its combined report when a haplotype had no
  off-target.** When the guide had no off-target on one haplotype (or none there that lifted to
  hg38), or none on either, `assembly-search` wrote the combined table and then failed with
  `TypeError: cannot concatenate object of type '<class 'numpy.float64'>'` (exit 1); on the web the
  job finished without a combined report. The combined table (`<job>_combined_hg38.tsv`) left out
  the columns of a haplotype without off-targets, and the report reads both haplotypes' columns. The
  table now always has both haplotypes' columns, blank for a haplotype without off-targets as on any
  row found on one haplotype only, and the combined report is built (with no off-target on either
  haplotype it says "No off-targets in this run."). A search with off-targets on both haplotypes
  writes the same table and report as before (compared byte for byte, apart from generation times
  and the random ids in the plots). **Affected:** such searches, from the CLI or the web, e.g.
  `--mm 0` with a guide that matches one haplotype only. **What to do:** run them again for the
  combined report. Covered by `PostProcess/test_assembly_reconcile.py`.
- **A search could run on the index of another genome whose name ends in the searched genome's
  name.** The pipeline picks an installed index by the genome's name,
  `genome_library/<PAM>_<N>_<name>`, but the pattern it matched also took
  `<PAM>_<N>_<prefix>_<name>`, the index of another genome: with `NRG_2_mini_hg38` installed, a
  search of `hg38` with bulges used it instead of building `hg38`'s own, and used it ahead of an
  installed `NRG_3_hg38` too, as the smaller bulge budget won; nothing in the output said so. The
  same happened between two individuals of the `Genomes/<individual>_<haplotype>` layout (before
  2.6.0) whose names end alike: reproduced with `assembly-search` of `X_WALK01` first, after which
  both of `WALK01`'s haplotypes logged `Reference index present (NRG_2_X_WALK01_...)` and reported
  `X_WALK01`'s sites. Only an index named for the genome itself is now taken. A library without two
  such names picks the same index as before, and reference and variant search outputs were compared
  byte for byte. **Affected:** searches with bulges, from the CLI or the web, of a genome `<name>`
  while `genome_library/` held an index of a genome named `<prefix>_<name>` (`NRG_2_mini_hg38`
  beside `hg38`; the legacy layout's `X_WALK01_*` beside `WALK01_*`). **What to do:** run them
  again; for `assembly-search`, with a new `--output`, since a completed haplotype result of the
  same `--output` is reused. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **Two searches that started together on a data folder without `genome_library/` could fail.** The
  pipeline makes `Results/`, `Dictionaries/`, `Genomes/` and `genome_library/` when they are
  missing, after checking for each; two runs that started at the same moment (two web jobs, or two
  `docker run`s) both found `genome_library/` missing, and the second to make it stopped with
  `mkdir: cannot create directory '/DATA/genome_library': File exists` and `ERROR: Genome enrichment
  failed on _` (exit 1). `crisprme.py` makes `Genomes/`, `Results/` and its other folders the same
  way, and could stop with `FileExistsError` there. A folder made meanwhile is now no error.
  **Affected:** searches started together on a data folder without those folders, such as the first
  ones before any index is downloaded or built there, from the CLI or the web. **What to do:** run
  the one that failed again. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **The personal-assembly results header, History and combined report named both haplotypes
  "genome".** The header and History showed `Paternal genome: genome | Maternal genome: genome`, and
  the combined report (`report.html` in `<job>_combined_report.zip`) the same in its Paternal genome
  and Maternal genome rows: the name of the folder every bundle haplotype is kept in, which
  `assembly-search` wrote into the job's `.Params.txt` when it finished (the web form's own label
  was overwritten). All three now show the bundle's assembly name (`assembly_name` in its
  `metadata.json`; for an HPRC fetch, the HPRC assembly name), or `<individual>_<haplotype>` when
  the bundle records none, from the CLI and the web alike; a legacy `Genomes/<name>` folder keeps
  its name. Jobs that already finished keep the label they were written with.
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
