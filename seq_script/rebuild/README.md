# The rebuild of the published indexes for CRISPRme 2.10.0

`index_set_2.10.0.tsv` is the **rebuild manifest**: every archive the rebuild for 2.10.0 builds
and publishes in [`lucapinello/crisprme-indexes-format2`](../../docs/hf/crisprme-indexes-format2_README.md),
the location of index format 2 that no release before 2.10.0 reads (decision C2), one line each,
with its panel, PAM file, PAM, side, spacer (guide) length, bMax (N, the name's), the build's
`--bDNA`/`--bRNA`, the archives of 2.8.1 it replaces, its build and publish commands, and a note.
A panel's companions (its Dictionaries bundle and genotype stores) have no command of their own
(`-`): the panel's first build writes them and its first `publish-index` uploads them.

| Who reads it | For |
|---|---|
| the rebuild | what to build, in which working directory, and how to publish it |
| the D7 readiness audit | `check_index_set.py`, before any build: the set against the code |
| `crisprme.py download` (`crisprme_hf.successor_index`) | which index of 2.10.0 replaces an archive of 2.8.1 a user asks for by name |
| the docs | `docs/PRECOMPUTED_INDEXES.md`'s table and the Hugging Face card list the same set |

## The set (approved 2026-11-08, "the lean set")

- SpCas9 **NRR** 20-nt (the default) on the reference, HGDP + 1000 Genomes, HPRC and mega;
- SpCas9 **NGG** 20-nt, the **general** PAM-less index (an 8-nt 3′ flank, 20-nt guide), SaCas9
  **NNGRRN** 22-nt (Manuel's CEP290 guides), Cas12a **TTTN** 5′ at 20 and 23 nt (every shipped
  Cas12a PAM and Manuel's Cas12a guides), and SpCas9 **NRG** 20-nt (Kris's TRAC reports), each on
  the reference and HGDP + 1000 Genomes;
- **NNN_3 is not built again**: the general index serves its searches at the same speed, and a
  download of `NNN_3_*` says so; 2.8.1's archives stay where they are;
- each variant panel's **Dictionaries bundle** and genotype store once
  (`dictionaries-bundle_<dataset>`, `genotypes_<dataset>`, `indel_genotypes_<dataset>`), and its
  enriched-genome kit (`build-inputs/<dataset>`).

Every index is of N = 3 (`--bDNA 2 --bRNA 2`), the web's ceiling of 2 bulges of each type, which
it serves for its own guide length (`index_fit.fits`).

## Rules the rebuild keeps

1. **One working directory per panel**, every index of the panel built in it (`$HGDP`, `$HPRC`,
   `$MEGA`; `$REF` for the reference indexes), each with hg38 in `Genomes/hg38`, one FASTA per
   contig, as 2.8.1's default index of each panel was built (the contigs of its archives):
   `$REF`'s and `$HGDP`'s hold all of hg38, as `crisprme.py download --what genome` installs it
   (455 contigs, chrM and the alt, random and chrUn ones among them, as 2.8.1's reference indexes
   and its NRR, NGG and NNN indexes of HGDP + 1000 Genomes hold; its NRG index of the panel held
   chr1-22, X and Y); `$HPRC`'s and `$MEGA`'s hold its chromosomes chr1-22, X and Y alone, the 24
   contigs of 2.8.1's indexes of those panels. mega's recipe builds every chromosome of the merged
   VCFs in `$MEGA_SRC` (those 24 in `vcfs/hg38_mega`) and refuses a genome folder that holds any
   other FASTA. Copy them from the reference's:

   ```bash
   mkdir -p $HGDP/Genomes && cp -r $REF/Genomes/hg38 $HGDP/Genomes/
   for d in $HPRC $MEGA; do
     mkdir -p $d/Genomes/hg38
     for c in chr{1..22} chrX chrY; do cp $REF/Genomes/hg38/$c.fa $d/Genomes/hg38/; done
   done
   ```

   The rebuild's script checks the genome folders before its first command (rule 5); on a test
   of a few chromosomes, `$HGDP`'s holds the FASTAs of the reference's, `$HPRC`'s and `$MEGA`'s
   those of them that are of chr1-22, X and Y, and mega's merged VCFs are of those (the HGDP +
   1000 Genomes and HPRC recipes build the chromosomes of their genome folder alone, whatever
   other chromosomes' sources `$RP` and `$HPRC_SRC` hold). The first build enriches the genome and
   writes the Dictionaries; each later PAM re-indexes the same enriched genome ("Enriched genome
   already present"). So the panel's indexes share one Dictionaries, whose fingerprint every one
   of their manifests records: `publish-index` uploads the bundle and the genotype store with the
   panel's first index, skips them for the others (the repo holds them, the same), and refuses to
   publish an index built against another build of them (`--replace-companions` replaces them,
   for a whole panel built again).
2. **Build with the release image only** (the CRISPRme and CRISPRitz the release ships), after the
   D7 readiness audit: no fix that can change index creation is pending. `$SIF` is that image,
   which the HGDP + 1000 Genomes and HPRC recipes run themselves, and `$RUN` the prefix that runs
   the commands' `crisprme.py`, `hf` and `pigz` in it, e.g. `RUN="apptainer exec --bind /srv/local
   $SIF /usr/local/bin/_entrypoint.sh"` (empty where the release is installed). The recipes bind
   `/srv/local` (`hgdp1kgp_build.sh`), and `/srv/local/lp698` and `/data/pinello`
   (`hprc_build.sh`), into the image: their working directories and inputs lie there, and `$RUN`
   binds the working directories and the repository's checkout (mega's recipe runs the
   repository's `PostProcess/` in the image).
3. **Publish each panel as a unit**, then write `indexes/SHA256SUMS` from the Hub's own
   checksums once its uploads are done, `python scripts/hf_sha256sums.py
   lucapinello/crisprme-indexes-format2 > SHA256SUMS`, and upload it as `indexes/SHA256SUMS`
   (`docs/hf/README.md`). `publish-index` packs each archive, uploads it and deletes it, so no
   archive is left to hash before its upload, and `--staged staged.sha256`, which checks that
   the Hub holds what was hashed where it was built, has nothing to read: this list lets a user
   check a download against the Hub, and what the Hub holds is checked on archives downloaded
   from it (`docs/RELEASE_CHECKLIST.md`: every published variant archive can find indels,
   checked in the archive).
4. **Upload the set's PAM files** to `pams/` of `lucapinello/crisprme-data`, which
   `download --what pams` installs: every `pam_file` of the set and the search PAMs its indexes
   serve (`22bp-NNGRRT-SaCas9.txt`, `21bp-NNGRRT-SaCas9.txt`, `20bp-TTTV-Cas12a.txt`,
   `23bp-TTTV-Cas12a.txt`), as `PAMs/` ships them (`docs/RELEASE_CHECKLIST.md`). Until then a
   working directory set up by `download` holds the default `20bp-NRR-SpCas9.txt` and the data
   repo's NGG, NNN and NRG files only, and the web offers no other PAM.
5. **The commands run from the repository's root**, each on its own: a `cd` is inside `( )`, so no
   command moves the next; a recipe is called by its repository path (`seq_script/merge_panels/`),
   with every input it would otherwise read from a path it was written with, the cluster's or its
   working directory's (`$RP`, `$META` of HGDP + 1000 Genomes; `$HPRC_SRC`, `$P1KG` of HPRC;
   `$MEGA_SRC` of mega), and stops, with a non-zero status, at a step that fails; each build
   first copies its PAM file from the repository's `PAMs/` into its working directory's `PAMs/`,
   which a working directory set up by `download` does not hold (rule 4), and reads it there; and
   every file a publish uploads is named by its working directory. `--commands` prints them as one
   script, which stops before its first command where one of these variables is not set or a
   working directory's genome folder is not as rule 1 says, and at the first command that fails.

## Commands

```bash
python seq_script/rebuild/check_index_set.py              # the check and its dry run (the D7 audit runs it): exit 1 on a problem
python seq_script/rebuild/check_index_set.py --commands > rebuild.sh   # the rebuild: every build command, then every publish command, in order
REF=... HGDP=... HPRC=... MEGA=... RUN=... SIF=... RP=... META=... HPRC_SRC=... P1KG=... MEGA_SRC=... bash rebuild.sh   # from the repository's root
```

The check holds when: every PAM file exists and gives the row's PAM, side and spacer; N = max(bDNA,
bRNA) + 1; every index's name is the one `build-index-only` gives its PAM file
(`index_fit.index_name`, guide length and side); every archive routes to the format-2 repo
(`crisprme_hf.resolve_index_repo`); every index serves its own guide at 2 + 2 bulges; every variant
index's reference index is in the set; every variant panel has one bundle and one kit, and a
genotyped one (HGDP + 1000 Genomes, HPRC) one genotype store and one store of phased indel
genotypes; and every index archive of 2.8.1 is replaced by exactly one archive of the set, the one
a download of it names; and every build and publish command is one bash reads (`bash -n`), a
companion having none, and keeps rule 5: each command reads only the variables above, and each of
them is read by one (a variable no command reads is an input a command lost), in its own panel's
working directory, with no `cd` outside `( )` and no repository path after one; every script it
calls is the repository's, and given every input it would otherwise read from a path it was written
with (one of the cluster's, or its working directory's: HGDP + 1000 Genomes' `RP=${RP:-$B}`); each
build copies its PAM file and reads that one; `crisprme.py`, `hf` and `pigz` run through `$RUN`;
every upload goes to the format-2 repo and no command passes a token (`HF_TOKEN` holds it); and
every variant index is published `--dictless`, so that each of a panel's indexes names one
Dictionaries bundle. Then the check runs `--commands`' script dry: every program it calls stubbed
(`crisprme.py`, `hf`, `pigz`, `tar`, `sha256sum`, a recipe; `apptainer` called directly fails),
from a copy of the repository's root made of links, in temporary working directories that hold
`Genomes/hg38` as rule 1 says (all of hg38 for the reference and HGDP + 1000 Genomes, chr1-22, X
and Y for HPRC and mega, and mega's merged VCFs of those) and what the builds leave. It must run to
its end, which it does only where every script it calls is there, every PAM file a build reads is
in place, every path a command hands a program or a recipe is there, and every index of the set is
published from its working directory; and it must write nothing into the repository's root.
`PostProcess/test_index_set.py` runs the check and the dry run in the unit tests, mutates the set
and its commands, and checks that each mutation is found: by the check and by the dry run where
both can see it (a recipe called after a `cd`, a PAM file not copied, a kit packed after a `cd`
that moves the next commands, a program of the image run on the host, a script the repository does
not have, a text the script cannot hold), by the dry run alone (a kit uploaded from the
repository's root, or packed into it), and by the check alone where the dry run cannot see it (an
input left to a recipe's default, a variable the rebuild does not set, or one no command reads,
another panel's working directory, a pending script the tree holds, an upload to another repo, a
token, a variant index published with its per-sample dictionaries, and the set's own errors); it
runs the script on genome folders that are not as rule 1 says, which stops before its first
command. `PostProcess/test_panel_recipes.py` runs the HGDP + 1000 Genomes and HPRC recipes on
stubs: each stops, with a non-zero status, at a step that fails (the samples, a chromosome's VCF,
the build), and a chromosome of the genome folder without its source stops it; on a test genome
of chr22, or of chr21 and chr22, beside a source folder of every chromosome's source, each builds
those chromosomes alone and reads no other's source.

The mega rows follow lane M's full rebuild of `hg38_mega` (decision D2: PASS/'.' records, AF =
AF_max, its registry and indel allele frequencies consistent), from the merged VCFs
`mega.<chr>.afmax.vcf.gz` (`$MEGA_SRC`: `vcfs/hg38_mega` of `lucapinello/crisprme-data`, or what
`seq_script/merge_panels/mega_gw_merge.sh` writes), whose recipe,
`seq_script/merge_panels/mega_rebuild.sh <install_dir> <merged_vcf_dir>`, that lane owns. Until
the integration of stage 2b brings it into this tree, the check names it as pending
(`check_index_set.PENDING_SCRIPTS`) and the dry run takes its call as done; once the tree holds
it, the check asks for its entry to go.
