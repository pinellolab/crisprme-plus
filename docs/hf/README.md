# Hugging Face dataset cards (tracked here, pushed on release)

`docs/RELEASE_CHECKLIST.md` lists the Hugging Face dataset card as a
release surface, but until 2.8.1 nothing in the repo held a copy — so it drifted
unnoticed: the live card documented 4 of the 12 index archives and opened with
"All indexes are NRG PAM", which stopped being true when NGG, NNN and NRR were
uploaded on 2026-10-25/27.

These files are the source of truth. Edit them here, review them in the PR, and
push them to Hugging Face as part of the release:

| File | Hugging Face dataset |
|---|---|
| `crisprme-indexes-format2_README.md` | [`lucapinello/crisprme-indexes-format2`](https://huggingface.co/datasets/lucapinello/crisprme-indexes-format2): the indexes of CRISPRme 2.10.0 and later (index format 2), each panel's Dictionaries bundle and genotype store once, and the enriched-genome kits (decision C2) |
| `crisprme-data_README.md` | [`lucapinello/crisprme-data`](https://huggingface.co/datasets/lucapinello/crisprme-data): the reference data of every release, and the index archives of 2.8.1 (index format 1) |
| `crisprme-data-general-index_README.md` | [`lucapinello/crisprme-data-general-index`](https://huggingface.co/datasets/lucapinello/crisprme-data-general-index): the general (8-nt flank) index archives of 2.8.1 |

Push with (requires a write token):

```bash
huggingface-cli upload lucapinello/crisprme-indexes-format2 \
  docs/hf/crisprme-indexes-format2_README.md README.md --repo-type dataset
huggingface-cli upload lucapinello/crisprme-data \
  docs/hf/crisprme-data_README.md README.md --repo-type dataset
huggingface-cli upload lucapinello/crisprme-data-general-index \
  docs/hf/crisprme-data-general-index_README.md README.md --repo-type dataset
```

**Gate:** after pushing, re-read the three live cards and diff them against these
files. A card that disagrees with `list_repo_files` is a release-blocking
inconsistency (`RELEASE_CHECKLIST.md` → "nothing on HF that is undocumented").

## The location of 2.10.0's indexes, and what it holds

`lucapinello/crisprme-indexes-format2` holds the index set of 2.10.0,
[`seq_script/rebuild/index_set_2.10.0.tsv`](../../seq_script/rebuild/index_set_2.10.0.tsv), every
archive of it and nothing else (`seq_script/rebuild/check_index_set.py` checks the set against
the code; `seq_script/rebuild/README.md` says how the rebuild follows it). `publish-index` routes an
index of format 2, or one named for its guide length and side, there by itself, and refuses
to publish an index of format 2 in the two repos of 2.8.1, which 2.8.1 reads. 2.10.0's `download`
and its Settings catalog read it by default; an archive of 2.8.1 is still fetched by its own
name from the repos of 2.8.1.

## `indexes/SHA256SUMS`

Each repo carries one, so a user can `sha256sum -c --ignore-missing SHA256SUMS` a download. The
Hub already records every LFS file's SHA-256, so the file is generated from that
metadata rather than by re-downloading the archives:

```bash
# archives about to be uploaded: hash them where they were built
( cd staged/indexes && sha256sum *.tar.gz ) > staged.sha256
python scripts/hf_sha256sums.py lucapinello/crisprme-data --staged staged.sha256 > SHA256SUMS
# upload the archives and SHA256SUMS together, then confirm the Hub holds what was staged:
python scripts/hf_sha256sums.py lucapinello/crisprme-data | diff - SHA256SUMS
```

That final diff must be empty: it is the check that what is on the Hub is
byte-for-byte what was validated. For a repo that does not exist yet, the staged
file alone is the whole list (`LC_ALL=C sort -k2 staged.sha256`, the order the
script writes).

## Each panel's Dictionaries and genotype store, once (2.10.0)

In `lucapinello/crisprme-indexes-format2` a variant archive holds the index, its `_INDELS`
index and `manifest.json` only. The panel's Dictionaries (`registry_`, `log_indels_`,
`indel_genotypes_`, mega's `indel_af_`, a classic publish's `dictionaries_`) and its samplesID
lists go once into `dictionaries-bundle_<vcf>.tar.gz`, its Tier-1 genotype store into
`genotypes_<vcf>.tar.gz` and its phased indel genotypes alone into `indel_genotypes_<vcf>.tar.gz`
(the genotype oracle), each with a record, `<name>.json`, of what it holds: a fingerprint of its
files (a gzip header's timestamp read as zeros, which no reader uses). The index's manifest
names the bundle and the store by name and fingerprint.

- **Publish.** `publish-index` of a panel's first index uploads its companions, then the
  index; of each later index, the index alone (the repo holds the companions, the same: "not
  uploaded again"). An index built against another build of them is refused before anything
  is uploaded, as the indexes already published read the ones in place; a panel built again as a
  whole replaces them with `--replace-companions`. So publish every index of a panel from the
  working directory it was built in.
- **Download.** `download --what index` installs the companions an index's manifest names once
  per panel: kept, and not fetched again, where the install holds them with that fingerprint
  (a marker, `Dictionaries/.<name>.json`, records them and the indexes installed with them);
  checked against the manifest (a repo whose companion is another build's is refused); and an
  archive of 2.8.1, which carries its own Dictionaries, is refused over a bundle of this repo,
  while a bundle that replaces another build's Dictionaries names the installed indexes that
  read them.

## The genotype companion travels with its index (the repos of 2.8.1)

`download --what index` fetches `genotypes_<vcf>.tar.gz` and the shared reference
index from the **same repo** the index came from. A variant index published to
`crisprme-data-general-index` therefore needs its own copy of
`genotypes_hg38_HGDP1kGP.tar.gz` there; upload the main repo's file unchanged
(same SHA-256) rather than re-packing it.

Conversely, `publish-index` of an index of format 1 re-packs and re-uploads that companion
whenever `Dictionaries/genotypes_<vcf>/` sits in the publishing working directory. One
companion serves every index built on the same panel, so republishing a single
index from such a directory replaces it for all of them. Publish from a working
directory that holds only `registry_`, `log_indels_` and `indel_genotypes_`, as
the 2026-10-25 uploads did. (An index of format 2 publishes it once: above.)

## Archives of 2.8.1 on one panel must bundle identical `Dictionaries/`

Every variant archive of 2.8.1 carries the panel's `Dictionaries/` (`registry_`,
`log_indels_`, …) and an install **replaces** each of those directories, so with
two indexes of one panel in a working directory, the one installed last decides
what both read. On 2026-10-27 `NRR_3_hg38+hg38_mega` shipped empty indel logs and
thereby emptied `NRG_3_hg38+hg38_mega`'s too (`INDEX_PROVENANCE.md`). Before
uploading a new PAM for an existing panel, compare its `Dictionaries/` member by
member with an archive already published for that panel, and take them from that
archive when in doubt. Compare content, not names and sizes: for 2.8.1 every member
was hashed (SHA-256), and equal sizes had hidden 46 compressed files that differed
only in the gzip header's timestamp (bytes 4-7, which no reader uses; identical once
masked). That difference is harmless; any other one is not. This holds across the two
repos too: a K = 8 install replaces the same folders as the NRG index beside it.
