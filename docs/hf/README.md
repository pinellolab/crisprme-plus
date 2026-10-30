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
