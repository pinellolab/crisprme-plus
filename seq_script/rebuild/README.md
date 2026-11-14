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
  Cas12a PAM and Manuel's Cas12a guides), and SpCas9 **NRG** 20-nt (the PAM of the genome-wide
  TRAC-guide runs made with 2.8.1), each on the reference and HGDP + 1000 Genomes;
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
   contig, as most of 2.8.1's indexes were built (the contigs of their archives):
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
2. **Build with the release's CRISPRme and CRISPRitz** (the code the release ships), after the
   D7 readiness audit: no fix that can change index creation is pending. The indexes of 2.10.0
   were not built with the release image, which is built from the release's tags, after the
   indexes: they were built with an image of the image route, `rc281f.sif`'s environment with
   CRISPRme's commit `db50986` and CRISPRitz's `4ded5e6` compiled in it, whose sha256, commits and
   binaries `docs/INDEX_PROVENANCE.md` records (2.10.0). `$SIF` is that image,
   which the HGDP + 1000 Genomes and HPRC recipes run themselves, and `$RUN` the prefix that runs
   the commands' `crisprme.py`, `hf` and `pigz` in it, e.g. `RUN="apptainer exec --bind /srv/local
   $SIF /usr/local/bin/_entrypoint.sh"` (empty where the release is installed). The recipes bind
   `/srv/local` (`hgdp1kgp_build.sh`), and `/srv/local/lp698` and `/data/pinello`
   (`hprc_build.sh`), into the image: their working directories and inputs lie there, and `$RUN`
   binds the working directories and the repository's checkout (mega's recipe runs the
   repository's `PostProcess/` in the image).
3. **Publish each panel as a unit** into `lucapinello/crisprme-indexes-format2` (the publish
   half, once: The publish, below), which must exist first: `publish-index` creates no repo, so
   its owner creates the dataset repo, and pushes its card,
   `docs/hf/crisprme-indexes-format2_README.md` (`docs/hf/README.md`), before the first
   publish. Then write `indexes/SHA256SUMS` from the Hub's own
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
   every file a publish uploads is named by its working directory. `--commands build` and
   `--commands publish` print them as the two halves of the rebuild, a script each (Commands): each
   stops before its first command where one of the variables its commands read is not set (the
   build half also where a working directory's genome folder is not as rule 1 says), and at the
   first command that fails.

## Commands

The rebuild has two halves, a script each. The build half builds every archive of the set, and
publishes nothing. The publish half publishes them, and builds nothing: it is run once, after the
indexes have passed the release's gates (`docs/RELEASE_CHECKLIST.md`) and the PI has approved the
publish. The script of 2.10.0's audit held both, every publish command straight after the builds,
so a run published before any gate, wherever a token or a node's cached login could upload
(PP-2 of the release).

```bash
python seq_script/rebuild/check_index_set.py                          # the check and its dry run (the D7 audit runs it): exit 1 on a problem
python seq_script/rebuild/check_index_set.py --commands build > rebuild.sh     # the build half: every build command, in order; it publishes nothing
REF=... HGDP=... HPRC=... MEGA=... RUN=... SIF=... RP=... META=... HPRC_SRC=... P1KG=... MEGA_SRC=... bash rebuild.sh   # from the repository's root
# once the indexes have passed the gates and the publish is approved:
python seq_script/rebuild/check_index_set.py --commands publish > publish.sh   # the publish half: every publish command, in order
REF=... HGDP=... HPRC=... MEGA=... RUN=... bash publish.sh            # from the repository's root
```

Each line of a half is its row's command, as the manifest gives it, then `|| exit` with what
failed and a comment of the row's archive and note; around it, an `if` and an `else` that skip a
row done before, and a mark once it has run, each on lines of their own. The lines of the build
half, the rows' lines alone, are those the D7 readiness audit checked and the rebuild of 2.10.0's
indexes ran (md5 `7389f7e234d4e617860f8f9b34abdc49`, 19 lines; `check_index_set.half_lines` takes
them out, as `PostProcess/test_index_set.py` does).

**A re-run.** Each half skips a row it ran to its end before: once a row's command has run, the
half leaves a mark in the row's working directory, `.rebuild-2.10.0/<half>.<archive>` (a `/` as
`_`), and a re-run skips a row with its mark (an index's with the index still in `genome_library/`).
The build half also passes a kit whose `SHA256SUMS` holds (`sha256sum -c` in its working directory:
the list is written once the tarball is whole), and mega's build where every part its recipe makes
is in `$MEGA`: its index and `_INDELS` index, its enriched and fake-indel genomes, its
`log_indels_`, `registry_` and `indel_af_` Dictionaries, its samplesID, and the build VCF of each
chromosome of `$MEGA_SRC`, with no staging folder left (`check_index_set.WHOLE_RECIPES`). The
recipe refuses a working directory that holds any part of the dataset, so its line, run again once
done, exited 1 and stopped the script there, before any later row. A build run row by row (as
the D7 plan's `rows.sh` ran 2.10.0's) leaves no mark: there, the build half builds each index row
again, which the tools keep where it is in place ("already present", a few seconds each), and
passes the kits and mega. To build a row again, remove its mark and what it built. A mega build
that stopped partway is refused by its recipe until its parts are gone: move `$MEGA` aside, put
back its genome folder, and run the build half again.

**The publish.** The publish half uploads with the release image's own `hf` and `crisprme.py`
(`$RUN`), with a token given in the environment, `HF_TOKEN`, a write token of
`lucapinello/crisprme-indexes-format2` (a fine-grained token lists the repo, which its owner adds
once the repo exists), never on a command line, which `ps` shows. The image sets
`HF_HUB_DISABLE_IMPLICIT_TOKEN=1`, so its `hf upload` sends no token unless that is overridden in
it: under apptainer by `APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0`, which the script exports
itself (a host variable `HF_HUB_DISABLE_IMPLICIT_TOKEN=0` does not reach the image), and with
Docker by the `RUN` given (`docker run -e HF_TOKEN -e HF_HUB_DISABLE_IMPLICIT_TOKEN=0 ...`). It
exports `PYTHONNOUSERSITE=1` too, so that the image's Python reads no package of `~/.local`,
which apptainer's bind of `$HOME` puts on its path. Before its first upload it checks, and stops
(exit 2) saying what to do, where: `HF_TOKEN` is not set; the token does not reach the image, or
the image would not send it; the `hf` that `$RUN` runs is not the one installed beside its
`crisprme.py` (the image's own: a login shell's `PATH` can put `~/.local/bin/hf` first, so run it
with a clean `PATH`); the token does not authenticate (`hf auth whoami`); or the repo is not
there (`hf datasets info`: its owner creates it, with its card, first; `hf upload` would create
a repo that is not there, under the token's account). The kits go up as the manifest gives it,
`$RUN hf upload lucapinello/crisprme-indexes-format2 <kit> build-inputs/<kit> --repo-type dataset`,
the command of `hf` in `huggingface_hub` 2.4.0, the image's. On a node:

```bash
read -rs HF_TOKEN && export HF_TOKEN     # typed, not echoed; or read from a file only you can read
python seq_script/rebuild/check_index_set.py --commands publish > publish.sh
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin REF=... HGDP=... HPRC=... MEGA=... \
  RUN="apptainer exec --bind /srv/local $SIF /usr/local/bin/_entrypoint.sh" bash publish.sh
unset HF_TOKEN
```

Run against a stand-in for the Hub on 127.0.0.1 (`HF_ENDPOINT`; it answers the token check and
refuses every upload) in the image the rebuild image was made from (`rc281f.sif`,
`huggingface_hub` 2.4.0), on ml003 and ml007, a kit's upload as the manifest gives it sent no
token without the override and sent it with the script's exports, as `docs/hf/README.md`'s
card upload did; the preflight passed with a token and stopped without one, and stopped where
the token did not reach the image (`apptainer exec --cleanenv`).

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
every upload goes to the format-2 repo and no command passes a token (`HF_TOKEN` holds it); no
build command publishes, and no publish command builds; and every variant index is published
`--dictless`, so that each of a panel's indexes names one Dictionaries bundle. Then the check runs
the two halves dry, the build half, then the publish half: every program they call stubbed
(`crisprme.py`, `hf`, `pigz`, `tar`, `sha256sum`, a recipe; `apptainer` called directly fails),
from a copy of the repository's root made of links, in temporary working directories that hold
`Genomes/hg38` as rule 1 says (all of hg38 for the reference and HGDP + 1000 Genomes, chr1-22, X
and Y for HPRC and mega, and mega's merged VCFs of those) and what the builds leave. Each must run
to its end, which it does only where every script it calls is there, every PAM file a build reads
is in place, every path a command hands a program or a recipe is there, and every index of the set
is published from its working directory, by the publish half alone; and neither may write into
the repository's root. Run again, each must do nothing again; and the build half, with the marks
of the kits and of mega's build gone and mega's parts all there, must build nothing (A re-run).
`PostProcess/test_index_set.py` runs the check and the dry run in the unit tests, mutates the set
and its commands, and checks that each mutation is found: by the check and by the dry run where
both can see it (a recipe called after a `cd`, a PAM file not copied, a kit packed after a `cd`
that moves the next commands, a program of the image run on the host, a script the repository does
not have, a text the script cannot hold), by the dry run alone (a kit uploaded from the
repository's root, or packed into it; a half that does again what it did on a re-run, or runs
mega's recipe on a working directory that holds every part it makes), and by the check alone where
the dry run cannot see it (an input left to a recipe's default, a variable the rebuild does not
set, or one no command reads, another panel's working directory, a pending script the tree holds,
an upload to another repo, a token, a variant index published with its per-sample dictionaries,
and the set's own errors), and by both (a build command that publishes); it runs the build half on
genome folders that are not as rule 1 says, which stops before its first command, and on a `$MEGA`
that holds every part mega's recipe makes, whose line it passes (the recipe runs, and refuses it,
where a part is missing); and it checks that the build half's lines are the 19 the rebuild ran.
`PostProcess/test_panel_recipes.py` runs the HGDP + 1000 Genomes and HPRC recipes on stubs: each
stops, with a non-zero status, at a step that fails (the samples, a chromosome's VCF, the build),
and a chromosome of the genome folder without its source stops it; on a test genome of chr22, or
of chr21 and chr22, beside a source folder of every chromosome's source, each builds those
chromosomes alone and reads no other's source.

The HPRC rows read `$HPRC_SRC`: HPRC Release 2's Minigraph-Cactus graph **v2.1** on GRCh38, the
VCF HPRC decomposed with `vcfwave` (`hprc-v2.1-mc-grch38.wave.vcf.gz`, 2,291,014,302 bytes, md5
`63cc82deb0e010ab98719a9e31dfa352`), split per chromosome into `ALL.<chr>.HPRC.vcfwave.vcf.gz` and
its `.tbi` for chr1-22, X and Y, as `docs/PRECOMPUTED_INDEXES.md` gives the commands
(`bcftools view --regions <chr>`, then `tabix -p vcf`). The panel moved to v2.1 for 2.10.0
(decision of 2026-11-13): 2.8.1's HPRC archives were built from v2.0, split the same way (the
CHANGELOG gives what differs). The recipe reads the samples of `ALL.chr22.HPRC.vcfwave.vcf.gz`,
the same 232 in both, in the same order, so the panel's samplesID is the same, and stops at a
chromosome of its genome folder with no VCF.

The mega rows follow the full rebuild of `hg38_mega` (decision D2: PASS/'.' records, AF =
AF_max, its registry and indel allele frequencies consistent), from the merged VCFs
`mega.<chr>.afmax.vcf.gz` (`$MEGA_SRC`: `vcfs/hg38_mega` of `lucapinello/crisprme-data`, or what
`seq_script/merge_panels/mega_gw_merge.sh` writes), with its recipe,
`seq_script/merge_panels/mega_rebuild.sh <install_dir> <merged_vcf_dir>`; the dry run checks its
call as it checks every other recipe's. (A recipe of a branch not merged yet would be named in
`check_index_set.PENDING_SCRIPTS`: the check names it as pending and the dry run takes its call
as done, until the tree holds it and the check asks for its entry to go. None is, now.)
