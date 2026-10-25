# Dynamic-PAM pam-less indexes — design & backward-compatibility contract

Design contract for the CRISPRitz "on-the-fly PAM filter" + CRISPRme's consumption of it.
Drives the CRISPRitz PR (`feat/dynamic-pam-filter`) and the CRISPRme-side wiring. Nothing
here changes existing indexes, commands, or results — it is strictly additive and opt-in.

## 1. The feature
A pam-less index built with an all-N PAM of width **K** stores, per leaf, the concrete K
flanking bases of the whole-genome protospacer space (PAM-independent k-mer tree + variants
+ bulges). At search time a restrictive PAM is applied as a **leaf-time skip-before-write**
filter (`--prune-pam`), so one index serves ANY PAM of length `L ≤ K` via N-padding — no
per-PAM rebuilds.

- **Validity rule:** a search PAM is valid for an index iff `L ≤ K` **and** at every position
  `search_bit[i] & build_bit[i] == search_bit[i]` (IUPAC subset). So an `NRG` index serves
  {NGG, NAG, NGA, NRG}; an all-N flank-K index serves anything ≤ K. You can only *narrow*.
- **Size cost (measured, 2×):** the flank is 4-bit nibble-packed → `+2 bytes/leaf` for 3→8,
  i.e. **~+9%** (ref index ~134 GB → ~146 GB). Linear, not multiplicative.
- **Proven:** NNN+`--prune-pam NGG` == a dedicated NGG index, identical loci (chr22 cross-check).

## 2. Backward-compatibility contract (two independent layers)

### Layer 1 — PHYSICAL isolation: old code cannot resolve pam-less indexes
Confirmed resolver facts (`PostProcess/crisprme_hf.py`):
- Default repo `DEFAULT_HF_REPO = "lucapinello/crisprme-data"` (overridable `--hf-repo` /
  `CRISPRME_HF_REPO`); indexes under `indexes/<name>/`, `name = {PAM}_{bMax+1}_{genome}`.
- `download --what index` **enumerates `indexes/` recursively** (`list_repo_files`) — so a
  same-repo subtree or new name WOULD be listed/fetchable by old code. A subtree is NOT safe.

**Therefore:** pam-less / variable-flank indexes are published to a **SEPARATE HF repo**
(proposed `pinellolab/crisprme-indexes-pamless`) that the old default never points at. Old code
cannot list or download them. New code's resolver knows both: legacy repo for legacy names, the
pam-less repo for pam-less indexes; `publish-index` targets the pam-less repo for variable-flank
indexes.

### Layer 2 — RUNTIME version-gate: fail-safe even if an index is hand-moved
Each pam-less index ships a sidecar `index_meta.json`:
```json
{ "build_pam": "NNNNNNNN", "flank_len": 8, "pam_at_start": false,
  "min_crispritz_version": "<this release>", "min_crisprme_version": "<consuming release>" }
```
New tools read it, validate the requested PAM (subset + `L ≤ K`), and version-gate on load. The
new `searchTST` **PAM-length guard** errors cleanly on any width mismatch. Old tools lack the
guard, which is exactly why Layer 1 (they never obtain the file) is the primary defense; a
deliberate manual copy still fails loudly (flank-tagged name + length guard on any new tool).

## 3. Format compatibility
No `.bin` format change — a pam-less index is a standard CRISPRitz index built with a longer
all-N PAM. New tools read legacy indexes with **byte-identical** default behavior. Naming is
flank-tagged + distinct (e.g. `NNNNNNNN_3_<genome>`, K=8) so it is never confused with a legacy
`NNN_3`. Adding the sidecar is invisible to old tools (they ignore unknown files).

## 4. CLI
- **CRISPRitz:** default call unchanged (byte-identical). `--prune-pam` (opt-in) applies the PAM
  file as a subset filter, validating subset + `L ≤ K` against the index name/sidecar.
- **CRISPRme:** `complete-search --index-path <idx> --pam <searchPAM>` — `--pam` may be any valid
  subset of the index's `build_pam`; CRISPRme **auto-passes `--prune-pam` to crispritz when
  `--pam ⊊ build_pam`** (narrower than build) and validates with a clear error on an invalid PAM.
  `--pam == build_pam` is the legacy path (unchanged). Inline `--genome --vcf` build unchanged.

## 5. Web
On index selection, read the sidecar → a **PAM picker** populated with that index's valid set
(dropdown of common subsets + a custom field with live subset/length validation); default =
`build_pam`. A one-line hint states the index's supported PAM family and max length.

## 6. NNN index publish plan
The genome-wide NNN HGDP+1kGP index (in progress) is published to the **new pam-less HF repo**,
flank-tagged name, with the `index_meta.json` sidecar (`build_pam`, `flank_len`, min versions) —
**not** the legacy `lucapinello/crisprme-data` repo. So no old client ever discovers it.

## 7. CHANGELOG / README wording (both repos)
State plainly: *existing indexes and commands are unchanged; pam-less (dynamic-PAM) indexes live
at a separate HF repo and require CRISPRitz ≥ `<X>` and CRISPRme ≥ `<Y>`.*

## 8. Web PAM picker
Legacy per-PAM indexes keep the **existing flat dropdown** (the `guidelen-PAM-enzyme` list)
unchanged. Pam-less indexes get a new **index-scoped, data-driven** picker:

- A committed reference table `enzyme -> {guide_len, pam, side, notes}` (the content of the
  existing dropdown: SpCas9 20bp NGG, SpG NGN, SpRY NRN/NYN/NNN, VQR NGA, VRER NGCG, xCas9,
  S. canis NNG/NNGT, S. thermophilus NNAGAA, Nme-20bp NNNNGATT, Cj NNNNACA, Thermo NNNNCNAA,
  SaCas9 21bp NNGRRT, Nme2 24bp, Cas12a/CasX/Cas12b 5' TTTV/TTCN/ATTN, …).
- On index selection, FILTER that table by the index metadata: show a preset iff
  `guide_len == index.guide_len AND side == index.side AND len(pam) <= K AND pam ⊆ build_pam`.
  So the index's geometry prunes the menu — an un-servable enzyme never appears for an index.
  (A 20-nt / 3' / K=8 pam-less index serves the whole 20-nt 3'-PAM family; it does NOT serve
  SaCas9/Nme2 (different guide length) or Cas12a/CasX/Cas12b (5' PAM) — those get their own index.)
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
Not a breaking change: a new, opt-in index **type** + capability, quarantined to a new HF repo,
version-gated, distinctly named; defaults and legacy paths untouched. The pam-less index trades
per-search speed/RAM for one-index-any-PAM flexibility.
