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
| `crisprme-data_README.md` | [`lucapinello/crisprme-data`](https://huggingface.co/datasets/lucapinello/crisprme-data) |
| `crisprme-data-general-index_README.md` | [`lucapinello/crisprme-data-general-index`](https://huggingface.co/datasets/lucapinello/crisprme-data-general-index) |

Push with (requires a write token):

```bash
huggingface-cli upload lucapinello/crisprme-data \
  docs/hf/crisprme-data_README.md README.md --repo-type dataset
huggingface-cli upload lucapinello/crisprme-data-general-index \
  docs/hf/crisprme-data-general-index_README.md README.md --repo-type dataset
```

**Gate:** after pushing, re-read both live cards and diff them against these
files. A card that disagrees with `list_repo_files` is a release-blocking
inconsistency (`RELEASE_CHECKLIST.md` → "nothing on HF that is undocumented").

## `indexes/SHA256SUMS`

Each repo carries one, so a user can `sha256sum -c SHA256SUMS` a download. The
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

## The genotype companion travels with its index

`download --what index` fetches `genotypes_<vcf>.tar.gz` and the shared reference
index from the **same repo** the index came from. A variant index published to
`crisprme-data-general-index` therefore needs its own copy of
`genotypes_hg38_HGDP1kGP.tar.gz` there; upload the main repo's file unchanged
(same SHA-256) rather than re-packing it.

Conversely, `publish-index` re-packs and re-uploads that companion whenever
`Dictionaries/genotypes_<vcf>/` sits in the publishing working directory. One
companion serves every index built on the same panel, so republishing a single
index from such a directory replaces it for all of them. Publish from a working
directory that holds only `registry_`, `log_indels_` and `indel_genotypes_`, as
the 2026-10-25 uploads did.

## Archives on one panel must bundle identical `Dictionaries/`

Every variant archive carries the panel's `Dictionaries/` (`registry_`,
`log_indels_`, …) and an install **replaces** each of those directories, so with
two indexes of one panel in a working directory, the one installed last decides
what both read. On 2026-10-27 `NRR_3_hg38+hg38_mega` shipped empty indel logs and
thereby emptied `NRG_3_hg38+hg38_mega`'s too (`INDEX_PROVENANCE.md`). Before
uploading a new PAM for an existing panel, compare its `Dictionaries/` member by
member with an archive already published for that panel, and take them from that
archive when in doubt.
