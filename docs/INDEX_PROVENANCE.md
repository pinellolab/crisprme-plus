# Published-index provenance and validity

Which published index is valid, for which client, and why — with the measurement
behind each claim. Written for the 2.8.1 release, when CRISPRitz 2.9.1 changed
the index builder and a second, non-backward-compatible index family appeared.

> **The indexes of 2.10.0** are rebuilt by CRISPRitz 2.10.0 in index format 2 and published in
> their own repo, `lucapinello/crisprme-indexes-format2`, from the index set
> [`seq_script/rebuild/index_set_2.10.0.tsv`](../seq_script/rebuild/index_set_2.10.0.tsv), which
> gives each archive's build and publish commands (`seq_script/rebuild/README.md`); their
> provenance is the next section. The rest of this page covers the archives of 2.8.1, which stay
> in `lucapinello/crisprme-data` and `lucapinello/crisprme-data-general-index` for 2.8.1.

## The indexes of CRISPRme 2.10.0 (index format 2)

The rebuild for 2.10.0 began on 2026-11-14 (ml007 and ml008): the 26 archives of the index set,
16 indexes, each panel's companions and its enriched-genome kit. Every row ran from one image and
one commit, so an index's bytes depend on that image and on the thread count alone. What was
fixed before the first row is recorded here; what each row records once built is marked *after
the build*, to be filled in from the build's logs before the publish, and the size and SHA-256 of
each index and companion archive, which exist only once it is published, *after the publish*,
from the Hub (Each archive, below; the release checklist, Gate D, "Provenance").

### The build image: the image route, not the release image

The indexes were **not** built with the release image, which `seq_script/rebuild/README.md`'s
rule 2 names: the release image is built from the Dockerfile at the release's tags, and the
indexes are built before the tags exist (CRISPRitz's `v2.10.0` among them). They were built with
an image made by the image route of the D7 readiness audit: the environment of `rc281f.sif`, a
release candidate of 2.8.1 that was never published, with CRISPRme's tree replaced by the commit
below and CRISPRitz compiled inside it, from the commit below, with the Dockerfile's `g++` lines.
`rc281f.sif`'s Dockerfile equals v2.8.1's, and its own `crisprme.py` (md5 `4c527335…`, not v2.8.1's
`1fbc9c9f…`) is replaced; its conda environment is v2.8.1's (2.10.0's `environment.yml` adds
`pigz` alone, which the image holds). That image, and its
`.provenance.txt`, are the only artifacts that reproduce these indexes' bytes: keep them with the
release's records.

| | |
|---|---|
| Image | `crisprme-2.10.0-db50986.sif`, 4,513,169,408 bytes, sha256 `e95dd278ac76b17277528a89996815790186e86d9a67b21e3b03b3e180e06e16` |
| Made | 2026-11-14 at 02:18 UTC on ml007, unprivileged (apptainer 1.4.4), by the plan's `mkimage.sh`; copied to ml008 and ml003, the sha256 checked on each |
| Base | `rc281f.sif`, sha256 `c10c34f903b0fb28d8c40fabedf047b29c70028fe7073bd66a766006599fc116`, checked before the build |
| CRISPRme | commit `db50986cc54e4e39533be8b7725d94aa1664b2bc` of `feat/2.10.0` (version 2.10.0), as a `git archive` (sha256 `6ec147db0a1a0267aed5285543a59f3d16c0edb898c2182673b134979cb512e1`); the image's `/opt/conda/opt/crisprme` equals the archive, 548 files by md5, and its `crisprme.py` (md5 `10148d8b8b2c61fb399bbda70c53265d`) the commit's |
| CRISPRitz | commit `4ded5e62230fcea4a4031ae1bf431e64344b5309` of its `feat/2.10.0` (`v2.10.0` not yet tagged), as a `git archive` (sha256 `f92e5801f9864d5e7141ffc448960b561797dfe917a08318f7b88aba6fecf360`): VERSION 2.10.0, INDEX_FORMAT 2, INDEL_LEFT and INDEL_RIGHT 30; the enricher compiled without `-fopenmp`, as the Dockerfile's line compiles it |
| Binaries (md5) | `buildTST` `f6a1c29a246caff13d04efe5bc8eaf1a`, `searchTST` `b16b9e340556f1046bd359420bc0f535`, `searchBruteForce` `eefeb865acdacf5cdd32f84f646feee8`, `enricher` `b2a5d3a78c7713853d6ad259babe383e`, `crispritz.py` `070a462e50c24befde4791b791055912` |
| Tools | bcftools 1.24, ijson 3.6.0 (yajl2_cffi), huggingface_hub 2.4.0 |
| Environment | every row under `env -i`, with `PYTHONNOUSERSITE=1` (the image's Python packages, not those of `~/.local`, which the bind of `$HOME` puts first: ijson 3.5.1 there, which parses the dictionaries the same) |
| Checked | `crisprme.py --version` v2.10.0 in the image, on ml007 and on ml003 (apptainer 1.3.4); the scorer's self-test computes; the smoke (the reference index and HPRC on chr22, mega on chr3, chr4 and chr12) records `built_by` 2.10.0/2.10.0 |

At the release, check that the tagged sources that build indexes are those used:
`git diff --stat 4ded5e6 v2.10.0 -- sourceCode crispritz.py Makefile Makefile_conda conda` in
CRISPRitz prints nothing, and in CRISPRme `git diff db50986 v2.10.0 -- PAMs
seq_script/merge_panels seq_script/rebuild/index_set_2.10.0.tsv` prints only comments of the
recipes, which the docs wave after the build changes (the HPRC recipe's header says v2.1); the
19 row lines of the build half keep their md5 (`7389f7e234d4e617860f8f9b34abdc49`: those of
`check_index_set.py --commands build` after its `# build` line, the `if`, `else` and mark lines
around each left out, as `check_index_set.half_lines` takes them; the printed script, which holds
its checks and those lines too, has an md5 of its own); and `git diff db50986 v2.10.0 --
crisprme.py PostProcess` changes nothing `build-index-only` runs (the changes after `db50986` are
search, publish and checker code, each named in the CHANGELOG).

### Threads

The `.bin` files of an index depend on the thread count it was built with (`buildTST` sorts in
parallel, and the order it writes depends on the threads), its targets do not: in the D7 audit,
an index of HGDP + 1000 Genomes built again from its kit at 64 threads, the count it was built
with, was byte-identical (7 of 7 `.bin` files), and at 16 threads 6 of 7 differed; searches of
indexes built at 16 and at 64 threads gave identical rows. Every row ran with its recipe's count,
under `env -i`, so that no `THREADS` of a shell reached it, and none was changed on a re-run:

| Archives | Threads | Where it is set |
|---|---|---|
| The 7 reference indexes (`20-NRR_3_hg38`, ...) | 16 | the manifest's `--thread 16` |
| The 7 indexes of HGDP + 1000 Genomes | 64 | `seq_script/merge_panels/hgdp1kgp_build.sh`, `THREADS` |
| `20-NRR_3_hg38+hg38_HPRC` | 32 | `seq_script/merge_panels/hprc_build.sh`, `THREADS` |
| `20-NRR_3_hg38+hg38_mega` | 16 | `seq_script/merge_panels/mega_rebuild.sh`, `THREADS` |
| Every `_INDELS` sub-index | 1 | `PostProcess/pool_index_indels.py` (`-th 1` a chromosome) |

The tier emitter's workers (`CRISPRME_TIER_WORKERS`, 3 on ml007) and the enrichment's change no
file: each chromosome's tiers are written alone. A user who rebuilds an index from its kit
(`docs/PRECOMPUTED_INDEXES.md`, step 3) gets the published `.bin` files with these counts.

### Inputs

Each input was checked against a list before the build (the lists are kept with the build's
records; each list's own sha256 is given here):

| Input | What | List (sha256 of the list) |
|---|---|---|
| hg38 | the 455 FASTAs of `download --what genome`, md5 equal to HF's | `hg38_455.md5` (`58f53f0e9b9dd85cfa597c4925202909094b51871b0074778f2e5d3ca75017ee`) |
| HGDP + 1000 Genomes | 2.8.1's as-run normalised VCFs, `hgdp1kgp.<chr>.norm.vcf.gz` and `.tbi` (46 files), sha256 equal to HF's `vcfs/hg38_HGDP1kGP`, put in place of the recipe's VCF-PREP (which reproduces them record for record) | `hgdp_vcfs_281.sha256` (`fb91de7fe0e1299290196d0373d9eab9429db2db120a28544e53307103a15b74`) |
| | the re-phased BCFs the recipe reads first | `hgdp_rp.sha256` (`1891be055ad8a4d6d8c55cee1a31251b7f11ec98612b9a8c9fa166be63b5abf1`) |
| | gnomAD's metadata, `gnomad_meta_updated.tsv` | `meta.sha256` (`377b549762758217719a0715748b91ca2ec64a4bd4e816399a1d853f61a18946`) |
| HPRC | HPRC Release 2's Minigraph-Cactus graph **v2.1**, split per chromosome (below) | `hprc_src.sha256` (`83589242e6b7843bf51a657b43d223abb08221114b0e942e0276ec6e16e75350`) |
| | the samplesID of 1000 Genomes 2021, which names the HPRC samples' populations | `p1kg.sha256` (`fa27a591f71965d76372bf24f93a3a73589579db9993c0445a9e22a2a1a66352`) |
| mega | the merged VCFs `mega.<chr>.afmax.vcf.gz` and `.tbi` (48 files), sha256 equal to HF's `vcfs/hg38_mega` | `mega_src.sha256` (`a7b9671631fbd0e17f39b2ed7e8afbf74cb14b75d02229007dc75c70e3b6e6da`) |

**HPRC v2.1.** The source is `hprc-v2.1-mc-grch38.wave.vcf.gz`, the `vcfwave`-decomposed VCF of
the GRCh38 graph, from S3
`human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/`:
2,291,014,302 bytes, md5 `63cc82deb0e010ab98719a9e31dfa352`, its S3 multipart ETag
`3e33ef1c419591a0774722d0ca9e10c5-274` recomputed from the file (8 MiB parts) and equal; its `.tbi`
2,212,134 bytes. It is split per chromosome as `docs/PRECOMPUTED_INDEXES.md` gives the commands
(`bcftools view --regions <chr>`, then `tabix -p vcf`), into `ALL.<chr>.HPRC.vcfwave.vcf.gz` and
its `.tbi` for chr1-22, X and Y, whose sha256 are these (`hprc_src.sha256`):

```text
1e26e98e6c4c0cbef0b05b1acc115eb0fb681e05bd3730b96eab9f049fb1e689  ALL.chr1.HPRC.vcfwave.vcf.gz
34f2fc7636757989224a85229d357e831124d102d1936f3f6e66dfce86cec8a5  ALL.chr1.HPRC.vcfwave.vcf.gz.tbi
6cd9311fbb25769c326c4b022e1a4a7be89cf838946d529509502f1ca16a2902  ALL.chr10.HPRC.vcfwave.vcf.gz
49a50b8b07a37649ed2436d282c1805a51c4f1c03bcbe20637c064c6d7686b48  ALL.chr10.HPRC.vcfwave.vcf.gz.tbi
a0f1592a8cc8bcde05d5022f73cc032b19c3c17936a571a84fba595dee803076  ALL.chr11.HPRC.vcfwave.vcf.gz
c08dace5a5e7f42ee03b3bf8c2cb0d807da58f81496630116ee6f702d1acfccf  ALL.chr11.HPRC.vcfwave.vcf.gz.tbi
373b9bdefe5f7f7dffa22cdd5f9522a8bd9a94410bbcb0fbae6525a62024e493  ALL.chr12.HPRC.vcfwave.vcf.gz
5442b7953e590fdea5fda434fc2dadbcb50834ed4b47dd5b38849e3ab13a2e46  ALL.chr12.HPRC.vcfwave.vcf.gz.tbi
37d18433341612a79ee3a3c2c5170124d283fb58d6bd255776886b605ae1310b  ALL.chr13.HPRC.vcfwave.vcf.gz
86d6e1b41e90f34362e2f95808a249c9be41423ffdb03c2904165309aab82555  ALL.chr13.HPRC.vcfwave.vcf.gz.tbi
88d5b0152f02787ac4ac5694a1ef2624765f72bcb5b1fa2a1b888511a09723cc  ALL.chr14.HPRC.vcfwave.vcf.gz
55e683e3cadd7fd187632abb89fe12830820a94aec96a38e17173b93b935796f  ALL.chr14.HPRC.vcfwave.vcf.gz.tbi
b927c2efa7f641e2626b0bd3661f20a1921da911e401a985a8a2b5d5d10cb731  ALL.chr15.HPRC.vcfwave.vcf.gz
a95aacef2d200332804f2b254e20cdcda85ae699bf98224394be6d4be73c57d3  ALL.chr15.HPRC.vcfwave.vcf.gz.tbi
27ef9b089d51904a59a3b89f556adf37803aede6c03670f918f4bdfe48b2d282  ALL.chr16.HPRC.vcfwave.vcf.gz
6b602506ba43ef3cc3f39c4ac2b5dfbc6dc78549093155bcde532925a6402bfe  ALL.chr16.HPRC.vcfwave.vcf.gz.tbi
1208b7b3515e77329d1f023ddfdfd9f5d1587f86a739efaf2591572e9f0a40cf  ALL.chr17.HPRC.vcfwave.vcf.gz
a22eb700c22e7dd38ed01911ce1610f8a502e30074d04fc78c6df365ebbc92d4  ALL.chr17.HPRC.vcfwave.vcf.gz.tbi
8cf98d1bb7375433fb2c35c89ee1556034282a9016312c7be5632d7733781bc1  ALL.chr18.HPRC.vcfwave.vcf.gz
cf0b2ee832b2c5ff116cd8cf2b2a131a02c7062e5b9fbddd5840e9c6ec7322d4  ALL.chr18.HPRC.vcfwave.vcf.gz.tbi
4a2c3b86a58e0041122286fa3b8cf0215af12cbfaabc7379b965fb6dfc52fe4e  ALL.chr19.HPRC.vcfwave.vcf.gz
896640baccdb5555fbaa41760417232ceec5dc7653b3323c74d05e1a4fd2fe43  ALL.chr19.HPRC.vcfwave.vcf.gz.tbi
c0ca5f26e29ed5fb20a0c57a748fa7783df82bdb4cf885a8ec0846ef3b17ecbc  ALL.chr2.HPRC.vcfwave.vcf.gz
1e6f92f727ffc31a30c23751eade450c8c49394deef92c16b0dc200f151bd7f9  ALL.chr2.HPRC.vcfwave.vcf.gz.tbi
de170e59f32a19ec79bf8b6617e531ec6121576e8bbe39566ff82385656a69ad  ALL.chr20.HPRC.vcfwave.vcf.gz
811830c436f87380bf7d33cdb29c9f7afe44285641e06ff1b21b1948e1287150  ALL.chr20.HPRC.vcfwave.vcf.gz.tbi
40ca4681103bd64bb28c7b8d17c9461b31928e3ec78fd0e8751c9b70f1d9830c  ALL.chr21.HPRC.vcfwave.vcf.gz
7aaea0681329af709c40c7aa055d0e18b005b80e079f10ed6c224ad204b2685e  ALL.chr21.HPRC.vcfwave.vcf.gz.tbi
69f8e290918bb58b16c867a03ff171b6e465049bcdd3890607b3a7700887a812  ALL.chr22.HPRC.vcfwave.vcf.gz
0317f6a779ae87af9978737218c7a85a44570ed65480291fa1b7ba9c369d5424  ALL.chr22.HPRC.vcfwave.vcf.gz.tbi
77487d780463c5e1f06f441077bd4ac268ef4cfbfc2c1a6392b43eb9c4d9cbd8  ALL.chr3.HPRC.vcfwave.vcf.gz
19555a711b2d19bad25028b27588c85b9381592d88239e7e91f2d4ff2cd6c8b2  ALL.chr3.HPRC.vcfwave.vcf.gz.tbi
ced2e7b891c02fc17eb4436e8c1862e84042a08e121ef0fc4aee7dda57038606  ALL.chr4.HPRC.vcfwave.vcf.gz
e248942b0519bae03334a5e7fe835c0eedaab1a19198c74786ccb1a3312898ce  ALL.chr4.HPRC.vcfwave.vcf.gz.tbi
2cb125575563871a3ad9d75b0ba5a5c50d82d77c60348fcd0017bc2acc3ecddb  ALL.chr5.HPRC.vcfwave.vcf.gz
3029e9f5377f3d2b59690e1bd82ca9d2831a591479e8afb0207250d522f5f6e2  ALL.chr5.HPRC.vcfwave.vcf.gz.tbi
ee41eea0afb9ebca552621a7c901f7078f2571fcd320ad0201d035f352a4b0e9  ALL.chr6.HPRC.vcfwave.vcf.gz
bcccad92e449b8751b8fc79334c940868e03caae6be0fb3a896d93d76e3eb13f  ALL.chr6.HPRC.vcfwave.vcf.gz.tbi
02fe8b9932f1ac024910bf7f82467566bad7f95b58ec7c30ecdfa1b29c8d4e86  ALL.chr7.HPRC.vcfwave.vcf.gz
d78e083f59c78b53ebd8bf3c63be3c4c562a8da444815f773f7b9e09363f3fbd  ALL.chr7.HPRC.vcfwave.vcf.gz.tbi
97c4848f8e6cce7aea4229473771598f55911d27e30e23c9f72d99e3c13c0a15  ALL.chr8.HPRC.vcfwave.vcf.gz
c00b7bffab7431c05c12518dda7ea47642fb01680cadb9b20633751752e5325a  ALL.chr8.HPRC.vcfwave.vcf.gz.tbi
110446d9fa5807ad61cee06ec109a5459596a2b796610760308056cc70544eb4  ALL.chr9.HPRC.vcfwave.vcf.gz
5bbb2d84f2389705a11a46db4041f55d1c1a1759209358081bb33b66641e59b7  ALL.chr9.HPRC.vcfwave.vcf.gz.tbi
87be3d502253a6a140e9d849033fa3b3f869845cfc88147c12317915d1643771  ALL.chrX.HPRC.vcfwave.vcf.gz
5e9b1d32237dbefffb0ba90603228932f81aaaa36886eab356dbcd2115e3b4de  ALL.chrX.HPRC.vcfwave.vcf.gz.tbi
1504e59263583dfe020265ebf2025ae45515543f3f04d3f976a70f53b2fe9895  ALL.chrY.HPRC.vcfwave.vcf.gz
a5a7e44fb9c48129fd88b948ab2e099f8c42141e38fc95457cd2b7a91dfee5fd  ALL.chrY.HPRC.vcfwave.vcf.gz.tbi
```

Its 232 samples are v2.0's, in the same order, so the panel's samplesID equals 2.8.1's, byte for
byte (sha256 `9fc1ce527d5390beb2e00fa9541a012c74279ab8925411a9449b00e7056180ff`). 2.8.1's HPRC
archives were built from v2.0 (the CHANGELOG, 2.10.0, Changed, gives what differs). The index reads
the recipe's normalised output of these files (`hprc.<chr>.norm.vcf.gz`), which the publish hosts
with the panel's kit (`docs/RELEASE_CHECKLIST.md`, T-35).

### Each archive (after the build, and after the publish)

The node and the threads are fixed. From the build's logs (`rows.status`, each row's `.time`
and `.out`), once the build is done and before the publish: the wall time and the peak memory of
each row; and each kit's size and SHA-256, from the `<dataset>.SHA256SUMS` its row writes beside
the tarball, with its member list's own SHA-256 (`<dataset>.members.md5`, the md5 of each file
the kit holds, which `scripts/kit_members.py` writes, as the tarball itself is not
byte-reproducible; the list goes beside the kit, `build-inputs/<dataset>.members.md5`). The size
and SHA-256 of the 16 index archives and the 7 companions are known only once they are published:
`publish-index` packs each archive at the publish, with a `manifest.json` that records when it was
packed, uploads it and deletes it, so no build log holds them, and a tarball packed by hand would
not be the Hub's (the manifest's time, and each file's time that tar records). They are filled in
after the publish, from the Hub: each SHA-256 from `indexes/SHA256SUMS`, which
`scripts/hf_sha256sums.py` writes from the Hub's own checksums (the release checklist, E, step
5), and each size from the Hub's listing of `indexes/` (`list_repo_tree`).

| Archive | Node | Threads | Wall time | Peak memory | Size | SHA-256 |
|---|---|---|---|---|---|---|
| `20-NRR_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NGG_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NRG_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NNNNNNNN_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `22-NNGRRN_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `TTTN-20_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `TTTN-23_3_hg38` | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NRR_3_hg38+hg38_HGDP1kGP` (the panel's first: its companions) | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NGG_3_hg38+hg38_HGDP1kGP` | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NRG_3_hg38+hg38_HGDP1kGP` | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `20-NNNNNNNN_3_hg38+hg38_HGDP1kGP` | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `22-NNGRRN_3_hg38+hg38_HGDP1kGP` | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `TTTN-20_3_hg38+hg38_HGDP1kGP` | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `TTTN-23_3_hg38+hg38_HGDP1kGP` | ml007 | 64 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `dictionaries-bundle_hg38_HGDP1kGP`, `genotypes_hg38_HGDP1kGP`, `indel_genotypes_hg38_HGDP1kGP` | ml007 | (its first build's) | — | — | *after the publish* | *after the publish* |
| `build-inputs/hg38_HGDP1kGP` | ml007 | (pigz) | *after the build* | *after the build* | *after the build* | *after the build* |
| `20-NRR_3_hg38+hg38_HPRC` (its companions) | ml008 | 32 | *after the build* (the pilot, with its kit: 03:49-06:39 UTC) | *after the build* | *after the publish* | *after the publish* |
| `dictionaries-bundle_hg38_HPRC`, `genotypes_hg38_HPRC`, `indel_genotypes_hg38_HPRC` | ml008 | (its first build's) | — | — | *after the publish* | *after the publish* |
| `build-inputs/hg38_HPRC` | ml008 | (pigz) | *after the build* | *after the build* | *after the build* | *after the build* |
| `20-NRR_3_hg38+hg38_mega` (its bundle) | ml008 | 16 | *after the build* | *after the build* | *after the publish* | *after the publish* |
| `dictionaries-bundle_hg38_mega` | ml008 | (its build's) | — | — | *after the publish* | *after the publish* |
| `build-inputs/hg38_mega` | ml008 | (pigz) | *after the build* | *after the build* | *after the build* | *after the build* |

| Kit | Its member list | SHA-256 of the list |
|---|---|---|
| `build-inputs/hg38_HGDP1kGP` | `build-inputs/hg38_HGDP1kGP.members.md5` | *after the build* |
| `build-inputs/hg38_HPRC` | `build-inputs/hg38_HPRC.members.md5` | *after the build* |
| `build-inputs/hg38_mega` | `build-inputs/hg38_mega.members.md5` | *after the build* |

## The archives of CRISPRme 2.8.1 (index format 1)

Two separate questions per archive:

1. **Readability** — can an older CRISPRme release read it? Decided by the
   **stored flank width**, which equals the PAM length in the archive's name.
   `loadTST` consumes exactly `ceil(len(search_pam)/2)` stored PAM bytes per leaf,
   so a stored width that differs from the search PAM desynchronises the leaf
   stream at the first leaf.
2. **Content** — does the 2.9.1 builder produce the same index? Decided by whether
   the genomic-`N` flank fix (CRISPRitz `5c612c0`) changes that PAM's build.
3. **Completeness** — does the archive hold what a variant search needs? Checked in
   the archive itself: every variant archive was downloaded (SHA-256 equal to the
   Hub's LFS oid) and its members listed, for `_INDELS` bins per sub-index and the
   bundled `Dictionaries/`.

## Verdict

**Every index in `lucapinello/crisprme-data` stores a 3-nt PAM region, so every
one of them is readable by any CRISPRme release, 2.8.0 included.** Wider PAM-less
flanks are the only exception and live in a separate repo. Readable is not enough
for the rebuilt `NNN_3` archives: a 2.8.0 search of them can fail or report malformed
rows (*With 2.8.0*, below), so use them with 2.8.1 or later.

Sizes are download sizes as the Hugging Face API reports them, in GB (10⁹ bytes);
"Uploaded" is the date of the commit that last wrote the archive, so for the four
2.8.1 replacements it is the replacement's upload.

| Archive | Size | Uploaded | PAM | Readable by 2.8.0 | 2.9.1 builder content | `_INDELS` in the archive |
|---|---|---|---|---|---|---|
| `NRG_3_hg38` | 8.9 GB | 2026-08-25 | NRG | yes | bit-identical | — |
| `NRG_3_hg38+hg38_HGDP1kGP` | 27.7 GB | 2026-10-07 | NRG | yes | bit-identical | 23 sub-indexes, 34 bins |
| `NRG_3_hg38+hg38_HPRC` | 15.2 GB | 2026-09-21 | NRG | yes | bit-identical | 24, 129 |
| `NRG_3_hg38+hg38_mega` | 12.1 GB | 2026-09-21 | NRG | yes | bit-identical | 24, 32 |
| `NRR_3_hg38` | 20.1 GB | 2026-10-25 | NRR | yes | bit-identical | — |
| `NRR_3_hg38+hg38_HGDP1kGP` | 42.9 GB | 2026-10-25 | NRR | yes | bit-identical | 23, 70 |
| `NRR_3_hg38+hg38_HPRC` | 30.7 GB | 2026-10-30 | NRR | yes | bit-identical | 24, 248; **replaced** (the 2026-10-27 archive's registry metadata was stale) |
| `NRR_3_hg38+hg38_mega` | 25.2 GB | 2026-10-30 | NRR | yes | bit-identical | 24, 67; **replaced** (the 2026-10-27 archive had 24, 0: empty) |
| `NGG_3_hg38` | 3.7 GB | 2026-10-25 | NGG | yes | bit-identical | — |
| `NGG_3_hg38+hg38_HGDP1kGP` | 21.5 GB | 2026-10-25 | NGG | yes | bit-identical | 23, 23 |
| `NNN_3_hg38` | 71.5 GB | 2026-10-30 | none (all-N, width 3) | yes, but a search can fail or report malformed rows (*With 2.8.0*) | **rebuilt with it, replaced** (the 2026-10-25 archive differs) | — |
| `NNN_3_hg38+hg38_HGDP1kGP` | 109.5 GB | 2026-10-30 | none (all-N, width 3) | yes, but a search can fail or report malformed rows (*With 2.8.0*) | **rebuilt with it, replaced** (the 2026-10-25 archive differs) | 23, 245 (the 2026-10-25 archive: 23, 222) |
| `genotypes_hg38_HGDP1kGP` | 5.7 GB | 2026-10-07 | — | yes | not an index | — |
| `genotypes_hg38_HPRC` | 0.6 GB | 2026-10-27 | — | yes | not an index | — |

**2.8.1 replaced `NRR_3_hg38+hg38_mega`**, which shipped an empty indel index,
**`NRR_3_hg38+hg38_HPRC`**, whose registry lacked the stamp that selects the
per-sample default, and **both `NNN_3` archives**, rebuilt with the 2.9.1 builder (all
below). The replacements were uploaded on 2026-10-30, each at its archive's existing
path.

An install replaces each `Dictionaries/` folder of its panel, so every replacement and
the K = 8 archive were checked member by member (SHA-256) against what the published
NRG archive of the same panel installs. The mega and HPRC replacements are
byte-identical to it. The HGDP1kGP archives (the rebuilt `NNN_3` variant archive and
the K = 8 one) carry the published `NNN_3`'s dictionaries byte for byte, and those
match NRG's in content: `registry_` and `samplesIDs/` are byte-identical, and the 46
compressed `log_indels_`/`indel_genotypes_` files differ only in the gzip header's
4-byte timestamp (identical once it is masked), which no reader uses. The published
`NNN_3` has carried exactly these dictionaries since 2026-10-25, and the dedicated NRG
index reproduces its gate report reading them.

## The 10 dedicated-PAM archives: bit-identical, measured

The N-gate commit claims a dedicated index built by the patched builder is
bit-identical to one from the pristine builder. Measured rather than quoted —
same chr22, same PAM file, both builders, driven through
`crispritz.py index-genome`:

| PAM | pre-fix builder | post-fix builder | verdict |
|---|---|---|---|
| NRG | 276,820,413 B, 3 bins | 276,820,413 B, 3 bins | bit-identical |
| NGG | 136,238,465 B, 2 bins | 136,238,465 B, 2 bins | bit-identical |
| NRR | 533,448,996 B, 5 bins | 533,448,996 B, 5 bins | bit-identical |

Why: the relaxation applies only when the PAM states **no** constraint at any
position. `NGG`/`NRG`/`NRR` each constrain two positions, so a genomic `N`
(encoded `0000`, intersecting nothing) still fails the PAM match exactly as
before.

## The 2 `NNN_3` archives: rebuilt for 2.8.1

The archives published on 2026-10-25 were built before the fix. An all-N PAM is fully
unconstrained, so the relaxation applies to them, and 2.8.1 replaced both with a 2.9.1
build, uploaded on 2026-10-30.

**What the fix changes** (chr22, same inputs, both builders):

| | bytes | bins |
|---|---|---|
| pre-fix builder | 1,801,785,226 | 16 |
| post-fix builder | 1,801,789,651 | 16 |

A difference of **4,425 bytes in 1.8 GB (0.00025 %)** — the extra entries are
targets whose 3-nt PAM region holds a genomic `N` while the 20-nt guide region
does not, i.e. sites at the edge of an assembly gap.

**Reachability in the engine**, counted in CRISPRitz's own search output, before
CRISPRme's post-search PAM filter, on a genome with such a site planted on both strands
(the minus-strand case matters: a minus-strand 3′ PAM sits at the *front* of the
extracted window):

| search PAM | pre-fix index | post-fix index | |
|---|---|---|---|
| `20bp-NNN-NO-PAM` (unconstrained) | 1 row | **3 rows** | the two extra rows are the planted `NNN`-PAM sites |
| `20bp-NRG-SpCas9` (`--prune-pam`) | 1 row | 1 row | identical — genomic `N` satisfies neither `R` nor `G` |

**What a CRISPRme 2.8.1 report shows.** 2.8.1 never reports a target with an `N`
inside its PAM: the post-search PAM filter drops it, and under `--prune-pam` the
engine already does. Measured through `complete-search` (reference-only, 2
mismatches, 1 DNA and 1 RNA bulge, rows of `integrated_results.tsv`) on a synthetic
genome with one control site and four gap-edge sites (PAM region `GNN` and `AGN` on
the plus strand, `TGN` and `GNN` on the minus strand), with an `NNN_2` index from the
pre-fix builder (CRISPRitz `5c612c0^`) and one from the 2.9.1 builder:

| search PAM | pre-fix index | post-fix index | |
|---|---|---|---|
| `20bp-NNN-NO-PAM`, with or without `--prune-pam` | 3 rows | 3 rows | byte-identical; without `--prune-pam` the filter drops the post-fix index's 639 extra targets |
| `20bp-NRG-SpCas9`, without `--prune-pam` | 2 rows | 2 rows | byte-identical |
| `20bp-NG-SpCas9NG` (`--prune-pam`) | 1 row | **5 rows** | the pre-fix index misses the two sites whose PAM (`AG`, `TG`) is followed by an `N`, and a 1-bulge row at each of the other two; a dedicated `NG_2` index gives the same 5 rows, byte-identical |

So in a 2.8.1 search with **any 3-nt 3′ PAM** the 2026-10-25 archives already behaved
exactly like a fresh build. Confirmed at full scale by a 2.8.1 release candidate: the
published `NNN_3_hg38+hg38_HGDP1kGP`, downloaded from Hugging Face and searched with
`NRG --prune-pam` genome-wide at 4 mismatches + 1 DNA + 1 RNA bulge with per-sample
resolution, reproduces the dedicated-`NRG_3` report **byte-for-byte**. What they miss is a
**shorter 3′ PAM** searched with `--prune-pam`: a site whose base just past that PAM is an
assembly `N` is a real
target there, since the engine tests the real PAM positions and not the pad
([METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)), and a
dedicated index of that PAM keeps it.

**The rebuild.** Genome-wide, with the 2.9.1 `buildTST` (md5 `76febff1…`), from the
2026-10-25 build's own inputs: the same hg38 FASTAs, the same enriched genome, the
same fake-indel genome (the md5 of all 23 files checked before and at the build), the
same 27-byte PAM file, and `-th 16`, because the builder's parallel sort is not
stable across thread counts. `Dictionaries/` and `samplesIDs/` are the published
archive's own, byte for byte, so only the index directories change. Where they
changed, against the published archive's member list, bin by bin, by exact size:

| index | bins (published → rebuilt) | bins changed | bytes | contigs changed |
|---|---|---|---|---|
| reference `NNN_3_hg38` | 1,628 → 1,628 | 1,222 | +145,540 | 66: exactly the contigs that contain `N` |
| variant `NNN_3_hg38+hg38_HGDP1kGP` | 1,628 → 1,628 | 1,221 | +106,849 | the same 66 |
| `_INDELS` | 222 → 245 | all | +2,105,770,267 (+9 %) | all 23 fake-indel contigs |

The fix only adds a target whose PAM region overlaps assembly `N`, so it predicts
that no contig without an `N` changes. None did: 0 of the 389 `N`-free contigs
changed, in the reference or the variant index. Every fake-indel contig is cut by an
`N` roughly every 56 bases, so all of them change, and that is where nearly all the
growth is; a 3-nt PAM such as NRG still filters those targets out, while a shorter
3′ PAM under `--prune-pam` reports those whose `N` lies just past it.

**Reproducible.** ml008 and ml007 each built the reference index independently, from
the same binary and the same inputs: 1,628 of 1,628 bins byte-identical.

**Verified.** The rebuilt pair was searched by a 2.8.1 release candidate genome-wide with
`NRG --prune-pam` (4 mm + 1 DNA + 1 RNA bulge, per-sample, HGDP + 1000 Genomes), with
`pam_filter` dropping nothing on any of the three searches. Every column is
byte-identical to the gate's report except `Annotation_COSMIC`: that run used a fresh
`download --what annotations`, where COSMIC stays off until its licence is attested
(its log: "COSMIC EXCLUDED ... dropped 742 COSMIC feature(s)"), so the column reads
`NA` throughout. With that column removed, the two reports are byte-identical. For a
user still on 2.8.0, the released 2.8.0 image, unmodified, searching chr22 on the
rebuilt reference index and on the published one gives **byte-identical** reports,
both with NRG (the NNN index as a stand-in, how 2.8.0 serves a PAM it has no dedicated
index for; 1,488 rows) and with the unconstrained NNN PAM (8,108 rows). Both hold only
because no hit there has an `N` in its PAM: one that does makes 2.8.0 fail or
misreport.

**With 2.8.0.** Use the rebuilt archives with 2.8.1 or later. 2.8.0 reads them, but it
has no `--prune-pam` (it ignores the flag), and its post-search PAM filter treats a
target's `N` as matching any base, so it keeps the targets with an `N` in their PAM
that the rebuild adds. What a search that reaches one does then depends on the
CRISPRitz under 2.8.0. On 2.9.1 the post-analysis fails on the `N` (`cluster.dict.py`).
On 2.9.0, which 2.8.0's image and its `install_from_source.sh` build, `loadTST` leaves
the character of a stored `N` unset (2.9.1 writes `N`), and the engine writes the
target cut short at its `N`: 21 or 22 nt where the others have 23 or 24. A search with
DNA bulges then fails in post-analysis (`IndexError: string index out of range`, then
`adjust_cols.py`: `ValueError: list.remove(x): x not in list`); one with RNA bulges
only finishes, exit 0, and reports the short target as a 0-mismatch, 0-bulge
off-target. Measured with 2.8.0's own code (`crisprme.py`, `PostProcess/`, `pages/`)
on each engine (2.9.1 as the 2.8.1 image ships it; 2.9.0 compiled from its tag with
the 2.8.0 Dockerfile's `g++` line, `searchTST` md5 `c13cfcf0…`), on the synthetic
gap-edge genome and the two `NNN_2` indexes of *What a CRISPRme 2.8.1 report shows*
above, reference-only:

| search PAM and budget | pre-fix index (either engine) | post-fix index, CRISPRitz 2.9.1 | post-fix index, CRISPRitz 2.9.0 |
|---|---|---|---|
| `20bp-NRG-SpCas9`, 2 mm + 1 DNA + 1 RNA (with or without `--prune-pam`) | 2 rows | exit 1: `ValueError: substring not found`, `ERROR: off-targets post-analysis (reference) failed` | exit 1: `IndexError: string index out of range`, the same `ERROR` line |
| `20bp-NNN-NO-PAM`, 2 mm + 1 DNA + 1 RNA | 3 rows | exit 1, the same error | exit 1, the same error |
| `20bp-NRG-SpCas9`, `--mm 2 --bDNA 0 --bRNA 1` (with or without `--prune-pam`) | 2 rows | exit 1, the same error | **exit 0, 4 rows**: 3 of them short 0-mismatch, 0-bulge rows with CFD 0.259, one at `chrA:984` in place of the real 1-bulge alignment there (`GAGTCCGAGCAGAAGAAGA-AAG`, CFD 0.016) |
| `20bp-NNN-NO-PAM`, `--mm 2 --bDNA 0 --bRNA 1` | 3 rows | exit 1, the same error | **exit 0, 5 rows**, 4 of them short 0-mismatch, 0-bulge rows |
| either PAM, `--mm 2 --bDNA 0 --bRNA 0` | 1 row | 1 row | 1 row: a reference search without bulges scans the genome, not the index |

The `--mm 0` runs at `--bDNA 0 --bRNA 1` give the same counts. Genome-wide, the NRR
search behind 2.8.1's PAM-filter fix (CHANGELOG) meets 1 such reference target, 1
variant target and 93,820 indel targets, so a 2.8.0 variant search through the
rebuilt `_INDELS` is likely to reach one. If you searched the rebuilt `NNN_3` with
2.8.0, run that search again with 2.8.1.

**Staged as searched.** The staged archives hold the build that was searched: rsync
checked every transferred file against a whole-file checksum, every file is in the
archive at its path with its exact size, and the gzip stream's CRC held on unpacking.
After the upload every entry of `indexes/SHA256SUMS` equals the Hub's LFS SHA-256 for
that archive, the post-upload check `docs/hf/README.md` prescribes for proving that
the Hub holds what was staged.

## `NRR_3_hg38+hg38_mega`: published with an empty indel index

The two questions above say nothing about whether an archive actually *contains*
what an indel search needs, so every published variant archive was also checked in
the archive itself: downloaded, SHA-256 equal to the Hub's LFS oid, members listed.
One failed.

| | published 2026-10-27 | replacement (on the Hub since 2026-10-30) |
|---|---|---|
| `_INDELS` | 24 sub-indexes, **0** `.bin` | 24 sub-indexes, 67 `.bin` (5.29 GB) |
| `Dictionaries/log_indels_hg38_mega` | 24 header-only logs (72–73 B each) | 419.1 MB, as in NRG mega |
| `Dictionaries/indel_genotypes_hg38_mega` | 48 empty files | absent, as in NRG mega |
| `Dictionaries/registry_hg38_mega` | data identical to NRG mega's; index metadata predates the `data_type`/`phased` stamps | NRG mega's |
| SNP index `NRR_3_hg38+hg38_mega` | 346 `.bin` | unchanged |

**Cause.** `build_nrr_ml007.sh` indexed ml007's `mega_gw_build` fake-indel genome,
built 2026-09-02, before the sites-only fix (#191): the enricher writes an indel
contig only when a sample carries the indel, a sites-only VCF has no samples, so
every chromosome's FASTA there is a bare header. `buildTST` crashed on each empty
genome, `os.system` ignored it, and `publish_nrr_mh.sh` uploaded the result with
that build's equally empty `log_indels`. The repo's own driver for this,
`seq_script/merge_panels/mega_build_indels.sh`, adds a synthetic pseudo-sample at
every indel site precisely to avoid it, and was not used.

**Effect.** Every indel search on NRR mega found nothing: silently on 2.8.0; on 2.8.1
the empty logs make it read as a panel with no indels, so also silently. And an
install replaces each `Dictionaries/` directory, so installing NRR mega next to
`NRG_3_hg38+hg38_mega` replaced NRG mega's 419 MB of indel logs with the empty ones,
and NRG mega lost its indel off-targets in that working directory too. (The other two
differences have no effect found: the registry's `.bin` data is byte-identical and,
without the stamps, the CLI and the web both fall back to population-level/sites-only
for this panel; an empty indel genotype store changes nothing for a panel with no
per-sample store, whose co-occurrence rows are PUTATIVE either way.)

**Repair.** The PAM-independent mega fake-indel genome that `mega_build_indels.sh`
produced on 2026-10-10 (ml007 `mega_kit_stage`) was re-indexed with the 2.9.1
builder. As a provenance check it was indexed for NRG first: the result is
**byte-identical to the published NRG mega `_INDELS`** (32 of 32 files,
2,137,365,874 B), so it is that index's exact input. Indexed for NRR, it gives the
67 `.bin` above. The replacement archive is the real `publish_index` output
(upload intercepted) from the published SNP index, that `_INDELS`, and the
`Dictionaries/` + `samplesIDs/` extracted from the published NRG mega archive.

**Verified.** Genome-wide, sites-only, population-level, 4 mismatches + 1 DNA + 1 RNA
bulge, guide `ACTGAAATCTGTAAGCAGGC`, all three runs exit 0:

| run | indel off-target rows (per-dataset AF) | integrated rows |
|---|---|---|
| published NRR mega, CRISPRme 2.8.0 | **0** | 379,288 |
| replacement NRR mega, 2.8.1 release candidate | 2,081,489 | 400,459 |
| published NRG mega, 2.8.1 release candidate | 1,131,494 | 188,744 |

And at the engine: on each of the 24 `_INDELS` sub-indexes, the replacement's NRR
targets passed through `pam_filter.py` with NRG **equal** the published NRG mega's
targets, 24 of 24 chromosomes, 1,131,494 targets in all (NRR's `N[AG][AG]` contains
NRG's `N[AG]G`, so nothing else would be consistent). Each engine-level total equals
the pipeline's indel row count above, so the pipeline carried every one through.

## `NRR_3_hg38+hg38_HPRC`: complete, with a registry that set the wrong default

Its indel index is complete (24 sub-indexes, 248 bins) and its index directories
match the publish source member for member. But it was published from ml007's
`hprc_build`, whose registry predates the `data_type`/`phased` stamps the other
archives carry (`reg_*.idx` 1,584 B against NRG HPRC's 1,625 B). `crisprme.py`
chooses the default analysis mode from that stamp — per-sample when it reads
`genotyped-phased`, population-level otherwise — so an NRR HPRC search silently ran
population-level, without named carriers or CONFIRMED cis phasing, on a phased
panel. And an install replaces `Dictionaries/registry_hg38_HPRC`, so installing NRR
HPRC next to NRG HPRC changed NRG HPRC's default too. The replacement, uploaded on
2026-10-30, takes its `Dictionaries/` + `samplesIDs/` from the published NRG HPRC
archive and changes nothing else.

## The wide-flank (K = 8) family: a separate repo

`NNNNNNNN_3` stores 8 PAM bytes. A release predating `--prune-pam` reads 2 bytes
per leaf for a 3-nt PAM, desynchronises, and (before 2.9.1's refusal) crashed —
reported as a clean empty result. These therefore go to
`lucapinello/crisprme-data-general-index`, which no older default points at, with
the split derived from the flank width at publish time rather than from a list of
names.

In `lucapinello/crisprme-data-general-index` (Size and Uploaded as defined above):

| Archive | Size | Uploaded | Readable by 2.8.0 | 2.9.1 builder content |
|---|---|---|---|---|
| `NNNNNNNN_3_hg38` | 82.2 GB | 2026-10-30 | no (8-nt flank) | built with it (the 2026-10-28 build below) |
| `NNNNNNNN_3_hg38+hg38_HGDP1kGP` | 122.7 GB | 2026-10-30 | no (8-nt flank) | built with it (the 2026-10-28 build below) |
| `genotypes_hg38_HGDP1kGP` | 5.7 GB | 2026-10-30 | yes | not an index; the main repo's file, unchanged |

**The 2026-10-28 build was published**, not the 2026-10-26 one:

| source | built | status |
|---|---|---|
| `n8_fixed_install/genome_library/NNNNNNNN_3_hg38` (1,628 bins) | 2026-10-28 17:01 | **published** |
| `n8_fixed_install/genome_library/NNNNNNNN_3_hg38+hg38_HGDP1kGP` (1,628 bins) + `_INDELS` | 2026-10-28 19:56 / 20:00 | **published** |
| `nnn_fromscratch/genome_library/NNNNNNNN_3_*` | 2026-10-26 11:25 | **not published**: pre-N-gate-fix |

The repo carries its own `genotypes_hg38_HGDP1kGP.tar.gz`. `download` fetches
the genotype companion (and the shared reference index) from the repo the index
came from, with no fallback, so without a copy there a per-sample search on the
K = 8 variant index would run with degraded Samples. It is the main repo's file
uploaded unchanged: both repos' `SHA256SUMS` list the same SHA-256 for it.

A wide PAM-less index is the one case the N-gate fix materially changes, and it
bites hardest on the synthetic indel genome, which is `N`-separated roughly every
56 bases: a width-31 window saw 34,690,315 usable starts against a width-26
window's 41,454,105, so **16.3 % of start positions were invisible** to a
pre-fix wide PAM-less build.

**Provenance of the 2026-10-28 build.** `buildTST` compiles from
`mainParallel.cpp` alone (`Makefile_conda`), which includes `pamGenerator.cpp`.
Both are byte-identical between the ml008 tree that produced the index and the
committed 2.9.1 release:

```
1bb392601d1ecd9fee138451e991540b  pamGenerator.cpp
6854efd77db8cd128a3f124422986097  mainParallel.cpp
```

`searchOnTST.cpp` did change afterwards (the width-mismatch refusal) but is not
compiled into `buildTST`, so it cannot affect index content.

## Gaps

- **Closed at the 2026-10-30 upload: `SHA256SUMS` under `indexes/`.** Until then
  neither repo had one, although `build-inputs/` has shipped them since the kits
  were uploaded. Both repos now carry one, written by `scripts/hf_sha256sums.py`
  from the Hub's own LFS checksums (an LFS object is addressed by its SHA-256; the
  oids match both shipped `build-inputs/*.SHA256SUMS`): 14 archives in
  `lucapinello/crisprme-data` and 3 in `lucapinello/crisprme-data-general-index`,
  each entry equal to that archive's LFS SHA-256. The "is what is on Hugging Face
  what we validated?" question is now answered by checksum instead of by re-running
  a gate arm (`docs/hf/README.md`).
- **Nothing stamps the builder version into an index.** The provenance above was
  reconstructed from upload dates, build mtimes and source hashes. A
  `.pam_build`-style sidecar recording the CRISPRitz version and git SHA would
  make this a lookup. Not done in 2.8.1.
- **An indel genotype store built by an enricher at or below 2.9.0 from a
  multiallelic VCF** can hold a fabricated compound allele under a real variant's
  rsID (CRISPRitz `cfcb6a6`). Latent on every CRISPRme-supplied VCF, which are
  split — 0 of 1,998,332 chr22 records carry a comma in `ALT` — and reachable only
  through `gnomAD-converter --multiallelic`. The fix does not retroactively touch
  an already-built store, and nothing invalidates one.
- **At a position carrying more than one SNP record, only the last record's ALT
  reaches the enriched genome** (CRISPRitz up to 2.9.1, both enricher ports, since 2.6.0
  in 2021; the SNP dict kept one record's entry a position up to then; in CRISPRitz 2.9.1's
  Known issues). Normalising a VCF splits a multiallelic site into
  biallelic records, and each record rebuilds the position's IUPAC code from the
  reference base and overwrites the SNP dict entry, so an earlier record's ALT is in
  neither the enriched genome nor the dict, and an off-target that needs it is not
  found, or is reported at a lower score, without its carriers (a VCF line that lists
  several ALTs itself is coded with all of them). It applies to any index CRISPRitz 2.9.1
  or earlier enriched from a VCF with such sites, every archive above among them.
  Measured for HGDP + 1000 Genomes: 2,486,165 of the
  panel's 127,666,429 SNP records (1.95 %), at 2,441,139 positions, are not coded for
  this reason; per guide (NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per-sample,
  genome-wide) this leaves 3 to 7 sites whose true CFD puts them in the top-1000 out
  of it or understated there, and at most one site reaching CFD ≥ 0.5 in the five
  guides tested (EMX1's chr17:12179615, 11 haplotypes, CFD 0.508), and 52 to 324 rows
  per guide (99 % of them at CFD < 0.1) name an individual who carries such an ALT,
  read as the reference base, for a sequence they do not have. On chr22, 45,837 of
  1,757,775 SNP records (2.6 %), at 43,389 positions, lose their ALT this way,
  counted in the normalized source VCF
  (`hgdp1kgp.chr22.norm.vcf.gz`) as the records not last at their position; 5,104
  of them have AC ≥ 41 (AF ≥ 0.5 % of 8,182 haplotypes) and 1,566 AC ≥ 818 (about 10 %).
  Confirmed directly in the shipped data at six such positions, the three of chr22
  with the highest AC, two at AC 257 and 82, and chr20:12296886: the earlier ALT is
  absent from the enriched genome's IUPAC code, the Tier-0 registry, the Tier-1
  genotype store and the per-sample dictionary, and only the last record's ALT is
  held. Two examples: chr22:19142014 has `rs807755` T>A at AC = 8,179 and T>G at
  AC = 1, and the enriched genome reads `K` (G/T), the stores hold G alone;
  chr20:12296886 has `rs6109350` G>A at AC = 1,451 and G>C at AC = 2, the code is
  `S` (C/G), and the stores hold C alone. The fix changes the enriched genomes, and
  with them every affected index: CRISPRitz 2.10.0 codes every ALT of such a site (`Q`
  where the reference base and its ALTs cover A, C, G and T) and lists each in the dict,
  and every index is rebuilt with it for CRISPRme 2.10.0.
- **A record whose `INFO` has `AF` at another place than the VCF's first record gets
  no `AF`** (CRISPRitz `add-variants` up to 2.9.1, both enricher ports: it looks for the
  `AF` key in the first record only and reads every record's `AF` at that place). Of the
  published indexes it affects HPRC's chrY alone (the NRG and NRR archives install the
  same logs), whose VCF has `AF` at five places: its indel log has no `AF` for 91,950
  of its 102,665 indels, which the report then shows with no MAF; its SNPs, whose
  frequencies the registry counts from the genotypes, are not affected. HPRC's other
  chromosomes and HGDP1kGP have an `AF` for each of their 12,000,599 and 17,192,381
  indels. The mega's VCF has no `INFO/AF` at all, so none of its 13,853,761 indels has
  one there, by design: their per-dataset frequencies are in `indel_af.tsv`. CRISPRitz
  2.10.0 finds each record's `AF` by its key, and HPRC's chrY log is rebuilt with it for
  CRISPRme 2.10.0 (METHODS section 8).

## How to re-run these measurements

The builder comparison needs two engines from the same tree, one with the three
N-gate files reverted to `5c612c0^`, each driven through its own conda-style
layout (`<root>/bin/crispritz.py` + `<root>/opt/crispritz/`). Two gotchas on a
cluster node: `apptainer exec` does **not** activate the image's conda
environment, so call `/opt/conda/bin/python3` by absolute path (`bash -lc` also
re-reads the profile and drops conda's `PATH`), and `/srv/local` is visible inside
the container only if it is the working directory or is bound with
`-B /srv/local`. Guard every comparison against an empty result: two empty hash
listings compare equal, so a failed build reads as "bit-identical" unless the
script counts the files it hashed.
