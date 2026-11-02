<!-- idx track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

### Fixed

- **A 0-bulge variant search on a downloaded index no longer loses every SNP-path variant row.**
  With `--bDNA 0 --bRNA 0` the variant search scanned `Genomes/<genome>+<dataset>/`, which a
  download does not ship. The engine printed `ERROR! ... is not a directory.` on stdout and exited
  0, so the run went on with the reference and indel rows alone, wrote a full report, and exited 0
  with an empty error log. **Affected:** every such search of a downloaded index without the
  enriched genome, on the CLI and in the web form's Advanced panel, in 2.8.0 and back to at least
  2.3.3 (searches with bulges, and installs with the enriched genome, were right). **What to do:**
  run those searches again. Measured on hg38 with the default index (TRAC guide, 4 mismatches, NRG,
  per-sample): 483 rows before, 722 after, which is byte-identical to a scan of the enriched genome.
  The 240 new sites are all variant sites; sites found only in the variant genome go from 80 to 319;
  sites at CFD >= 0.5 from 7 to 11, >= 0.2 from 33 to 50, >= 0.05 from 156 to 233; and 207 of the
  sites already listed change too (a variant alignment that scores higher, or another motif counted).
  The search now takes the installed variant index when any chromosome's enriched FASTA is missing,
  with the same `-var`, `--max-edits` and `--prune-pam` a search with bulges passes; with the whole
  enriched genome in place it scans it, as before.
- **A search the engine never ran now fails the run.** `crispritz.py` says what stops it before it
  searches on stdout and exits 0, so a search that never ran passed for one that found nothing.
  After each reference and variant search, the run now stops with an ERROR when the targets file is
  missing. A search that finds nothing writes a file with its header alone, and still ends in an
  empty result with exit 0.
