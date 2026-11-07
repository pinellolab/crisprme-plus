# CRISPRme validation benchmark (brute-force ground truth)

`crisprme.py validate-test` checks a CRISPRme installation by comparing the
off-targets that `complete-test` predicts against an **independent brute-force
ground truth** generated with a different method (exhaustive dynamic-programming
alignment instead of CRISPRme's TST index search).

## Files
- `brute-force-1000G/brute_force_1000G.tsv` — Cas9 reference (SpCas9, NGG PAM,
  guide `sg1617`, thresholds mm=4 / DNA-bulge=1 / RNA-bulge=1) on the hg38 +
  1000G variant-enriched genome.
- `brute-force-cas12-1000G/brute_force_cas12_1000G.tsv` — Cas12a reference
  (HBG1/HBG2 clinical guide, TTTV PAM). **Currently chr22-scoped** (953 off-targets;
  `"scope": "chr22"` in `benchmarks.json`, so `validate-test` skips it on any other
  chromosome and on a whole-genome run); run `validate-test --chrom chr22`. Regenerate
  genome-wide with the driver below, and drop its `scope`.
- `brute-force-simple/` — references for the two **default single-"n edits" web
  mode** cases (`cas9_sg1617_simple`, `cas12a_hbg_simple`): chr22, per-type caps
  wide open (mm=6/bDNA=2/bRNA=2) with a binding `--max-total-edits 3`, exactly the
  default search slider. These are marked `"heavy": true`: they build a bulge-2
  (NGG_3/TTTV_3) index over the variant-enriched genome. They run by default, locally and
  in CI: the whole set fits a 16 GiB, 4-CPU runner with no swap (measured peak 3.66 GB of
  RAM). `CRISPRME_SKIP_HEAVY=1` leaves them out, an opt-out for very small machines.
- `max-edits-5/` — references for the same two guides at **5 total edits**
  (`cas9_sg1617_n5`, `cas12a_hbg_n5`, also `"heavy": true`): chr22, the same bulge-2
  indexes, a binding `--max-total-edits 5`, with **381 and 92 sites**, 271 and 55 of them
  with an alignment holding two bulges of one type, which the per-type cases (one of each)
  do not search for. They are there
  because the two cases above can hardly catch a lost site: `cas12a_hbg_simple`'s
  reference is a header with no rows and `cas9_sg1617_simple`'s holds 3 rows at 2 sites,
  so a search that lost its sites still passed them. One lost site fails these. The
  folder is deliberately not named `brute-force-*`: the Bioconda recipe leaves those
  folders out of its package (`validate.py` then downloads them from `reference_base_url`,
  the main branch), where a reference that lands with a release is not there until that
  release is merged. These ship in every install, Docker, source and Bioconda, and need
  no network. The committed files came from the Rust generator run over 1 Mb pieces of the
  chromosome overlapping by 200 nt, sorted by position. `generate_references.py --only
  <name> --generator rust/target/release/brute_force_gen` on the whole enriched chr22 (about
  1.6 h on one core each) gives the same rows, 1988 and 454 with 0 differences, in another
  order, so with another md5.
- `generate_brute_force.py` — the ground-truth generator (see credit below).
- `generate_references.py` — registry-driven driver (regenerates any benchmark).

## Search modes: per-type vs single-"n edits"
A benchmark validates one of two modes, driven by its (optional per-case)
`thresholds` in `benchmarks.json`:
- **Per-type** (the original cases): independent `mm` / `bDNA` / `bRNA` budgets,
  so a target may carry a DNA *and* an RNA bulge. No total cap.
- **Single-"n edits"** (the DEFAULT web slider): leave the per-type caps wide open
  and add `"max_total_edits": n`. The binding constraint is total
  mismatches+bulges ≤ `n` (CRISPRme's `--max-total-edits`, each bulge one edit).
  The generator takes the matching `--max-total-edits n`; `generate_references.py`
  caps each per-type budget at `n` when it binds (proven set-identical, just
  faster). Both `generate_brute_force.py` and the Rust port accept `--max-total-edits`
  and produce the same set.

**Validation status:** on chr22, `generate_brute_force.py` reproduces the committed
Cas9 reference **exactly** (3495 rows, 0 diff) on the identical enriched genome,
confirming the algorithm + PAM/edge-gap filters. The Cas12a reference was produced
by the same validated generator.

## What `validate-test` covers, and what it does not

`validate-test` compares **sites** — the (chromosome, start, strand) triples of the
engine's raw CRISPRitz target reports — on two genomes: the plain reference and the
**SNV-only** 1000 Genomes chr22 enriched genome. That is the right check for the search
engine: it proves the TST index misses no site an exhaustive aligner finds.

It does **not** check:

- **indel rows.** The references are built from the SNV-enriched FASTA, which holds no
  indel, so the fake-indel genome and everything the indel post-analysis derives from it
  are outside the comparison. A wrong indel coordinate, a wrong REF text, a wrong carrier
  list or a dropped indel off-target all pass it. (The chr22 1000 Genomes VCF the
  self-test uses has no indel in `chr22:22322700-22322850` at all, which is why the
  2.8.0 row reported at `chr22:22322768` could not be caught here.)
- **anything but the site.** The aligned spacer/target strings, the REF/ALT text, the
  PAM, the mismatch and bulge counts, the named variants, the carriers and the whole
  alternative-alignments table are never compared with an independent truth.
- **per-carrier truth.** The enriched genome collapses every individual into IUPAC codes,
  so the reference cannot say which sample carries which sequence.

**`scripts/check_indel_rows.py` is the check for that part**, on real data: it writes
each sample's haplotypes itself from the VCF the search read (every allele of a slot, as
VCF 4.x defines a record), keeps its own map from haplotype bases to reference
coordinates and aligns with its own exhaustive dynamic-programming aligner (it imports
and copies no CRISPRme code, only Python's standard library; the CFD tables are read as
data). The truth is the VCF: a slot holds a row's target when it reads its bases from
the row's start through the record (each base the record puts in the target at its
place); a slot whose alleles overlap gives no sequence, so it is read without the alleles
METHODS section 4 reads it without (a phasing artefact at the record's position, an
allele only on bases it deletes, an MNV), the one reading that holds the record, and is
otherwise counted (`limit-overlapping-records`) where a target of its record reads, a
named variant it holds on no reading always failing a listed carrier; where a call
touching a target has unknown phase (a haploid call against a diploid one included) the
record's known slot is judged first, then METHODS section 4's PUTATIVE rule (the
strongest target the sample's alleles allow, by the alignment's CFD, its copies counted
against the record's). CRISPRme's documented reading of a carrier (the one-indel window
resolved with the carrier's own SNP alleles, a second indel read from either end through
the indel, the last record of a split multi-allelic site, the PUTATIVE carriers of MNVs)
only classifies a disagreement, and a carrier it would list against the carrier's own
haplotypes vetoes no fold: every disagreement it checks
for is a line, gating or a listed limit, and none passes in silence (what it does not
check is listed below). **Soundness**: every indel/MNV row block of the integrated
results and of the alternative alignments must be the given guide aligned to the
reference with exactly its named variants, from its reported start, and its carriers must
be the panel samples whose haplotypes hold it (or whose strongest target it is), within
the samples list the search read; a row of a sites-only dataset (a VCF with no sample
columns, as `hg38_mega`, alone or beside a genotyped one) lists no carrier, and only its
own strings, target and variants are checked. **Completeness**: every in-budget alignment
on a carrier's haplotypes that holds a base of an indel or MNV must be reported for that
carrier by rows no worse than its best alignment by the merge's own sort keys (the
fewest-mm+b block's key under `--sorting-criteria`, mm+bulges then mm by default; the
highest CFD at the three decimals the merge sorts on, recomputed from each row's strings,
then those criteria): a row of exactly that target whose fold group -- the rows of its
site (the cluster position on the record's haplotype) and strand naming the same variants
with the same carriers, which the documented merge folds into the first by each block's
order -- is no worse, or such rows with the target's own carrier set; a run that writes
no alternative alignments -- by default a population-level one, as `complete-test`'s runs
are (their `.search_mode` says so; `--ci-cd-test` chooses no mode) -- keeps one row per
merge cluster, and there the row of the target's cluster (rows chained within `--merge`
bp by the rows CRISPRme can write, as the merge clusters: the reference rows, the SNP rows
of the mode the run ran -- population-level each alignment's worst case, per sample the
haplotypes individuals hold -- and the targets some carrier of an indel record holds, none
that needs two indel records) must be no worse than each target; a cluster that reaches
such a row only through rows CRISPRme may write or not, which no documented rule settles
(a per-sample run's capped window's greedy row, which only the classic dictionary path
writes; the one-record reading of an undocumented overlap; a target two overlapping
alleles of a slot change), is counted under its cause's class, not passed. A slot whose
VCF alleles overlap is read as METHODS section 4 reads it where it says how, and judged
on every reading of the rest (a carrier whose slot holds the target on all of them, or on
none, is judged so); where its readings part, a listed carrier METHODS' rules take off
the target fails, and any other is counted. It does not check the CFD_REF/ALT, CFD risk
and CRISPR-Bulge scores (nor which row the CRISPR-Bulge block keeps),
`Variant_info_spacer+PAM`, `PAM_creation`, the seed counts, MAF, rsID, `Other_motifs`,
`Not_found_in_REF` or the annotations, nor the merge's keys after the sorting criteria and
the CFD (the Cluster Position, the carrier count, the aligned text), which only order
alignments equal in all of these. Give each VCF a `.tbi` when the panel has deletions longer
than 20 kb (HPRC's does): a read without one looks back 20 kb, and the checker lists such a
VCF under its own coverage class. Two genotyped datasets that share samples are refused
(exit 2): CRISPRme resolves each dataset on its own genotypes and joins their rows, which
the checker does not model; a sites-only dataset beside a genotyped one is checked. So are
inputs other than the search's where the run's own record of them says so (its
`.Params.txt` budget, PAM and variant datasets -- one VCF per dataset for each chromosome
checked --, its `.guides.txt`, `.sampleID.txt` and `.search_mode`, the mode it ran), a run
whose records contradict each other, and one that does not record what the check needs
(whether it wrote alternative alignments and, where it wrote none, its mode, which
`--analysis` can give).
`test/indel_snp_oracle/run_oracle.sh` runs it beside `oracle.py check`,
`validate-benchmarks.yml` runs it on the chr22 self-test output after `validate-test`, and
`docs/RELEASE_CHECKLIST.md` (step 4c) makes it a release gate on a real HGDP + 1000 Genomes
chr22 search in both modes and genome-wide. Its unit and mutation tests are
`PostProcess/test_check_indel_rows.py`.

## Reproducible generation pipeline
```bash
# 1. Build the variant-enriched genome the same way CRISPRme does internally
#    (per-chromosome FASTA with 1000G SNVs encoded as IUPAC bases).
#    This is the `add-variants` step of a normal complete-search/complete-test.

# 2. Run the brute-force generator per chromosome, e.g. for Cas9 (sg1617):
python test/benchmark/generate_brute_force.py \
    --fasta <chrN.enriched.fa> \
    --rna CTAACAGTTGCTTTTATCACNGG \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --pam GG \
    --chrom chrN --output chrN.cas9.tsv

#    ...and for Cas12a (TTTV 5' PAM, 20 nt spacer): PAM first, --pam-5prime
python test/benchmark/generate_brute_force.py \
    --fasta <chrN.enriched.fa> \
    --rna TTTV<CAS12_SPACER> \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --pam TTTV --pam-5prime \
    --chrom chrN --output chrN.cas12.tsv

# Or drive both from the registry:
#   python test/benchmark/generate_references.py --enriched-fasta <chrN.enriched.fa> --chrom chrN

# 3. Concatenate per-chromosome TSVs, sort, and record the MD5 (BFTARGETSMD5 in
#    PostProcess/utils.py) so validate.py can integrity-check the download.
```

Thresholds are kept at **mm=4 / bulge 1+1** to keep both the search and the
brute-force generation cheap (matches `complete_test.py`). Runtime note: the
generator is a pure-Python exhaustive scanner (~minutes per Mb), so a full
chromosome is a batch job.

## Faster generation: the Rust port

For large or genome-wide runs, use the drop-in Rust generator in [`rust/`](rust/)
instead of the Python script. It is std-only (no external crates), roughly an
order of magnitude faster, and produces **byte-for-byte identical output**
(validated against the committed chr22 reference: 3495 rows, 0 missing, 0 extra).
It takes the **same flags** as `generate_brute_force.py`.

```sh
# build once (needs a Rust toolchain: https://rustup.rs)
cd test/benchmark/rust && cargo build --release

# Cas9 example — same flags as the Python generator
./target/release/brute_force_gen \
    --fasta chrN.enriched.fa --rna CTAACAGTTGCTTTTATCACNGG \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --chrom chrN --pam GG --output chrN.cas9.tsv
```

See [`rust/README.md`](rust/README.md) for the full flag reference and the
Cas12a (5' PAM) example. Because the output is identical, you can generate
references with either implementation interchangeably.

## Important: genome provenance
The brute force must be generated from **the exact same variant-enriched genome
CRISPRme searches** (its `add-variants` output). Using a differently-enriched
FASTA yields spurious/missing hits in variant-dense regions (an all-IUPAC window
matches almost anything). Pin the enriched genome build alongside the reference.

## Credit / license
The dynamic-programming brute-force checker is by **Benjamin Vyshedskiy**
(https://github.com/benjaminvyshedskiy/Dynamic_checker). `generate_brute_force.py`
vendors it with attribution; only output formatting (tab-delimited, underscore
column names) was adapted. Developed in the Pinello Lab; included with the
author's permission.
