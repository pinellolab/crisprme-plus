# Dynamic-PAM pam-less indexes — design & backward-compatibility contract

> **Status: pre-implementation design (2026-10-24), not a description of CRISPRme+ 2.9.0.**
> What 2.9.0 ships: a search takes an index whose PAM is the search PAM, or all-N, or -- under
> `--prune-pam` -- a dedicated index of a PAM that holds the search PAM at every position (an
> `NRR` index serves NGG, NAG, NGA, NAA and NRG; an `NRG` one NGG and NAG; `NGA` is not inside
> `NRG`). An all-N flank as wide as the PAM (`NNN_3` for a 3-nt PAM) serves without the flag,
> narrowed by the post-search PAM filter; a wider flank (`NNNNNNNN_3`), for a 3′ PAM only, and
> a dedicated index of a broader PAM serve only pruned. `--prune-pam` no longer needs
> `--index-path`, and a search passes it by itself, and says so first, whenever no installed
> index serves it unpruned and one serves it pruned: on the command line and from the web
> interface, which launches the search script with pruning off. The run's record of the
> indexes it searched (`.indexes_searched.txt`, beside its parameters) and its report name the
> index that served the search. No `index_meta.json` sidecar is written or
> read. The nearest record is `.pam_build`, the PAM file an index was built from, which
> `build-index-only` writes into a reference or variant index, and a search into one it
> builds (the published 8-nt-flank archives have none). `NNN_3` is published in
> `lucapinello/crisprme-data` and the 8-nt-flank indexes in
> `lucapinello/crisprme-data-general-index`. The parts below marked *planned* are not in
> 2.9.0. What ships is described in
> [METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering) and the README.

Design contract for the CRISPRitz "on-the-fly PAM filter" + CRISPRme's consumption of it.
Drives the CRISPRitz PR (`feat/dynamic-pam-filter`) and the CRISPRme-side wiring. Nothing
here changes existing indexes, commands, or results — it is strictly additive and opt-in.

## 1. The feature
A pam-less index built with an all-N PAM of width **K** stores, per leaf, the concrete K
flanking bases of the whole-genome protospacer space (PAM-independent k-mer tree + variants +
bulges). At search time a restrictive PAM is applied as a **leaf-time skip-before-write**
filter (`--prune-pam`), so one index serves ANY 3' PAM of length `L ≤ K` via N-padding — no
per-PAM rebuilds. A 5' PAM needs an index built for the 5' side; the published ones are 3'
builds.

- **Validity rule:** a search PAM is valid for an index iff `L ≤ K` **and** at every position
  `search_bit[i] & build_bit[i] == search_bit[i]` (IUPAC subset). So an `NRG` index serves
  {NGG, NAG, NRG}; an all-N flank-K index serves any 3′ PAM ≤ K. You can only *narrow*.
  (This is the engine's rule. CRISPRme 2.8.1 narrowed only an all-N index; 2.9.0 narrows a
  dedicated one too, of a PAM as long as the search's: see the status note above.)
- **Size cost (measured, 2×):** the flank is 4-bit nibble-packed → `+2 bytes/leaf` for 3→8,
  i.e. **~+9%** (ref index ~134 GB → ~146 GB). Linear, not multiplicative.
- **Measured:** NNN+`--prune-pam NGG` == a dedicated NGG index, identical loci (chr22 cross-check,
  `L = K`). With `K > L` the flank does not index a + strand PAM with 1..K−L nt after it to the
  end of a contig, nor a − strand PAM in its first K−L nt, and reports one exactly K−L nt from the
  start on the + strand (on hg38 at K = 8: 145 + 101 NGG and 308 + 241 NRG sites, none on a
  primary chromosome; METHODS §9).

## 2. Backward-compatibility contract (two independent layers)

### Layer 1 — PHYSICAL isolation: old code cannot resolve pam-less indexes
Confirmed resolver facts (`PostProcess/crisprme_hf.py`):
- Default repo `DEFAULT_HF_REPO = "lucapinello/crisprme-data"` (overridable `--hf-repo` /
  `CRISPRME_HF_REPO`); indexes under `indexes/<name>/`, `name = {PAM}_{bMax+1}_{genome}`.
- `download --what index` **enumerates `indexes/` recursively** (`list_repo_files`) — so a
  same-repo subtree or new name WOULD be listed/fetchable by old code. A subtree is NOT safe.

**Therefore — the isolation boundary is FLANK WIDTH > 3, not "pam-less" per se.** The length-guard
(Layer 2) only bites on a width MISMATCH (search-PAM length != stored flank width). A **width-3**
pam-less index (e.g. `NNN_3`) searched with any 3-nt PAM has NO mismatch, so old code handles it as
an ordinary NNN index — it simply emits ALL 3-nt PAMs (per-PAM narrowing is what `--prune-pam`
adds). A width-3 NNN index is therefore format- AND width-safe for old code and **ships in the
standard repo**: CRISPRme 2.8.0 could download and search the first such archive (2026-10-25) with no
update (as an all-PAM / broadest survey — the top of the NGG->NRR->NNN ladder). The rebuilt archives
(2026-10-30) keep the targets with an `N` in their flank, and a 2.8.0 search of them can fail or
report malformed rows: use them with 2.8.1 or later. Only **variable/longer-flank** indexes (width > 3,
e.g. `NNNNNNNN_3` with an 8-nt flank) are unsafe for old code — a 3-nt search is a width mismatch
and old tools lack the guard — so those go in a **SEPARATE HF repo** (proposed
`pinellolab/crisprme-indexes-pamless`; 2.8.1 publishes them in
`lucapinello/crisprme-data-general-index`) the old default never points at. New code's resolver
knows both; `publish-index` targets the pam-less repo only for width>3 indexes.

### Layer 2 — RUNTIME version-gate: fail-safe even if an index is hand-moved
*Planned, not in 2.9.0 (which writes and reads no `index_meta.json`):* each pam-less index ships
a sidecar `index_meta.json`:
```json
{ "build_pam": "NNNNNNNN", "flank_len": 8, "pam_at_start": false,
  "min_crispritz_version": "<this release>", "min_crisprme_version": "<consuming release>" }
```
New tools read it, validate the requested PAM (subset + `L ≤ K`), and version-gate on load. The
new `searchTST` **PAM-length guard** errors cleanly on any width mismatch. Old tools lack the
guard, which is exactly why Layer 1 (they never obtain the file) is the primary defense; a
deliberate manual copy still fails loudly (flank-tagged name + length guard on any new tool).

**The sidecar is OPTIONAL — a CRISPRitz index is self-describing.** The build PAM + flank width are
encoded in the index name (`{PAM}_{bMax+1}_{genome}`) AND stored inside the index itself, so
`--prune-pam` validates a search PAM (subset + `L ≤ K`) directly against the index's own build PAM
with **no sidecar**. The current **width-3 NNN index needs NO sidecar** to be used with
`--prune-pam`. The sidecar only adds explicit min-version metadata + web-picker hints; the
`searchTST` length guard provides the actual hand-move safety. Ship a sidecar when convenient (it is
free and ignored by old tools), but nothing depends on it for correctness.

## 3. Format compatibility
No `.bin` format change — a pam-less index is a standard CRISPRitz index built with a longer
all-N PAM. New tools read legacy indexes with **byte-identical** default behavior. Naming is
flank-tagged + distinct (e.g. `NNNNNNNN_3_<genome>`, K=8) so it is never confused with a legacy
`NNN_3`. Adding the sidecar is invisible to old tools (they ignore unknown files).

## 4. CLI
- **CRISPRitz:** default call unchanged (byte-identical). `--prune-pam` (opt-in) applies the PAM
  file as a subset filter, validating subset + `L ≤ K` against the index name/sidecar.
- **CRISPRme** (as designed, 2.9.0 does this, with or without `--index-path`; in 2.8.1 the user
  passed `--prune-pam` with `--index-path`, and only an all-N index was narrowed):
  `complete-search --pam <searchPAM>` — `--pam` may be any valid subset of an installed index's
  `build_pam`. CRISPRme passes `--prune-pam` to crispritz by itself, and says so first, when no
  installed index serves the search unpruned and one serves it pruned (a dedicated index of a
  broader PAM, or an all-N flank wider than a 3′ PAM); an all-N flank as wide as the PAM serves
  unpruned, narrowed by the post-search PAM filter. A search for a PAM no installed index serves
  builds an index of that PAM where it can, as before: the reference index from the genome, the
  variant index from the dataset's VCFs or its enriched genome. Where it cannot, it is refused,
  naming the PAM and the indexes there: a reference index under `--index-path`, which a search
  never builds into, and a variant index of a dataset installed as a download (no VCFs, no
  enriched genome). `--pam == build_pam` is the legacy path (unchanged). Inline `--genome
  --vcf` build unchanged.

## 5. Web (2.9.0; in 2.8.1 the web interface never pruned)
The search script prunes by itself for a web search as for one on the command line. The PAM
picker is section 8: it reads the index's name and `.pam_build`, not a sidecar.

## 6. NNN index publish plan (*superseded*)
2.8.1 publishes `NNN_3` in `lucapinello/crisprme-data` and the 8-nt-flank indexes in
`lucapinello/crisprme-data-general-index`, with no `index_meta.json`. The plan was:
The genome-wide NNN HGDP+1kGP index (in progress) is published to the **new pam-less HF repo**,
flank-tagged name, with the `index_meta.json` sidecar (`build_pam`, `flank_len`, min versions) —
**not** the legacy `lucapinello/crisprme-data` repo. So no old client ever discovers it.

## 7. CHANGELOG / README wording (both repos) (*planned*)
State plainly: *existing indexes and commands are unchanged; pam-less (dynamic-PAM) indexes live
at a separate HF repo and require CRISPRitz ≥ `<X>` and CRISPRme ≥ `<Y>`.*

## 8. Web PAM picker (2.9.0, with the differences below)
2.9.0 ships this picker (`pages_utils.ENZYME_TABLE`, `pamless_presets`, `pamless_pam_problem`)
for the widest all-N index installed for the genome or dataset selected, with these differences:
the guide length and side come from the index name, `.pam_build` and the trees (`.bin`), not
from a sidecar; a preset is offered only for the index's own guide length (a 21-nt SaCas9 guide,
which a 20-nt index serves at up to 1 bulge of each type, is listed as not on this index); a 5′
PAM only as long as the flank; the custom field shows no site count; "Use this PAM" writes the
PAM file into `PAMs/` and selects it in the form's PAM list, and the search prunes by itself, so
no flag is passed; the recommended preset is SpCas9 NRR (the default PAM), not NGG.

The design, as written before it shipped:
Legacy per-PAM indexes keep the **existing flat dropdown** (the `guidelen-PAM-enzyme` list)
unchanged. Pam-less indexes get a new **index-scoped, data-driven** picker:

- A committed reference table `enzyme -> {guide_len, pam, side, notes}` (the content of the
  existing dropdown: SpCas9 20bp NGG, SpG NGN, SpRY NRN/NYN/NNN, VQR NGA, VRER NGCG, xCas9,
  S. canis NNG/NNGT, S. thermophilus NNAGAA, Nme-20bp NNNNGATT, Cj NNNNACA, Thermo NNNNCNAA,
  SaCas9 21bp NNGRRT, Nme2 24bp, Cas12a/CasX/Cas12b 5' TTTV/TTCN/ATTN, …).
- On index selection, FILTER that table by the index metadata: show a preset iff
  `guide_len == index.guide_len AND side == index.side AND len(pam) <= K AND pam ⊆ build_pam`.
  So the index's geometry prunes the menu — an un-servable enzyme never appears for an index.
  (A 20-nt / 3' / K=8 pam-less index serves the whole 20-nt 3'-PAM family, and 21-nt SaCas9 at up
  to 1 bulge of each type (METHODS §9 has the rule); it does NOT serve a longer guide such as Nme2's
  24 nt, or Cas12a/CasX/Cas12b (5' PAM) — those get their own index.)
- Layout: *SpCas9 & variants (recommended)* — **NGG (default)**, NRG, NRR, SpG, SpRY (NRN/NYN),
  VQR, xCas9, VRER; *Other 20-nt 3' enzymes (expandable)*; *Custom PAM…*; *Not available on this
  index (greyed, with reason + link to the matching index)*.
- Custom field: validate **on blur** — `length <= K` + subset + side — show `✓ valid (~N sites)`
  or `✗ <reason>` (e.g. "NYY isn't within this index's NRG-class sites"); no submit until valid.
- Chosen PAM → `--prune-pam <PAM>` to the search. Default preset = NGG.

## 9. Performance characteristics (measured)
A pam-less index is SLOWER + heavier per search than a dedicated per-PAM index, because the
`--prune-pam` filter is a leaf-time skip-before-write: it saves the OUTPUT + downstream
scoring/post-analysis, NOT the TST descent (the pam-less tree is traversed in full).

Measured (hg38 chr22:20-25Mb, EMX1, NGG 4mm/1/1, same host):
| search | off-targets | wall | peak RAM |
|---|---|---|---|
| dedicated NGG index | 772 | 1.3 s | 0.29 GB |
| pam-less NNN + `--prune-pam NGG` | 772 (identical) | 13.3 s | 3.36 GB |
| pam-less NNN, no filter | 9,736 | 15.6 s | 3.37 GB |

~10x slower / ~11x RAM for NGG. The penalty scales with PAM SELECTIVITY: largest for selective
PAMs (NGG ~10-16x, where a dedicated index is tiny), ~4x for NRR, and ~par for near-pamless
(no small dedicated alternative). Higher mm/bulge pam-less searches are proportionally heavier
(whole tree traversed regardless of the search PAM). **Guidance:** keep dedicated indexes as the
fast default for a fixed known PAM; use the pam-less index for PAM flexibility / variable /
arbitrary PAMs, accepting the per-search cost.

## Net
Not a breaking change: a new, opt-in index **type** + capability, quarantined to a new HF repo
(in 2.8.1 the wider flanks only), version-gated (*planned*), distinctly named; defaults and
legacy paths untouched. The pam-less index trades
per-search speed/RAM for one-index-any-PAM flexibility.
