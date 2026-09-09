# Testing CRISPRme+ 2.5.2 — guide for Ann & Manuel

Thanks for helping test the new release! This walks you through installing CRISPRme+ 2.5.2,
downloading the two new variant indices, running a couple of searches, and checking that the
new co-occurrence outputs look right. It should take ~30–60 min of hands-on time (plus
download/search wall-clock).

---

## What's new in 2.5.2 (what you're testing)

CRISPRme+ 2.5.2 ships **two complementary variant-aware indices** and new **co-occurrence**
reporting — off-targets that need *more than one* nearby variant on the same DNA molecule.

| Index | What it has | Co-occurrence reported as |
|---|---|---|
| **`NRG_3_hg38+hg38_1000G2021_HGDP`** (genotyped) | 1000 Genomes 2021 + HGDP, **sample-level genotypes** | **CONFIRMED** (phased cis) or **PUTATIVE** (unphased), with **named carrier samples** + exact joint allele frequency |
| **`NRG_3_hg38+hg38_mega`** (sites-only) | **5 sources** — 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us — **allele frequencies only** (no shared samples) | **PUTATIVE**, with a conservative **min-AF** joint bound + **per-dataset provenance** (which sources report each variant) |

Both indices report **two kinds** of co-occurring off-targets that the classic search could not see:
- **SNP + indel** on the same haplotype → `indel_snp_cooc.tsv`
- **SNP + SNP** (≥2 SNP alt alleles needed together) → `snp_snp_cooc.tsv`

The idea: with samples we report the **observed** haplotypes (CONFIRMED); without samples we
still report the **putative** haplotypes from allele frequency, so a co-occurring off-target is
never silently dropped just because a panel has no genotypes.

---

## 0. Install

### Option A — Docker (easiest; once the v2.5.2 image is published)
```bash
docker pull pinellolab/crisprme:v2.5.2
```
> If `v2.5.2` isn't on Docker Hub yet, the release tag is still building the image — use Option B
> in the meantime, or ping Luca.

### Option B — from source (works right now, `main` branch)
```bash
git clone https://github.com/pinellolab/crisprme-plus && cd crisprme-plus
mamba env create -f environment.yml
mamba activate crisprme-2.5.2
bash install_from_source.sh     # builds CRISPRitz 2.8.1 + puts crisprme.py on PATH
crisprme.py --version           # expect 2.5.2
```
(See README §1.3 for details / troubleshooting.)

---

## 1. Download the reference data + both indices

```bash
mkdir crisprme_test && cd crisprme_test
crisprme.py download --what all --path .                                              # reference genome + annotations
crisprme.py download --what index --index-name NRG_3_hg38+hg38_1000G2021_HGDP --path .  # genotyped index
crisprme.py download --what index --index-name NRG_3_hg38+hg38_mega            --path .  # mega sites-only index
```
> **Docker users:** prefix each command with
> `docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme:v2.5.2`
> and use `--path /DATA`.

The download writes `list_vcf.txt` / `list_samplesID.txt` for you (the CLI search reads these).

---

## 2. Search #1 — genotyped index (observed haplotypes)

Create `guides.txt` with your 20-mer + `NNN` PAM placeholder (one per line). Example (the TRAC guide):
```
CTCTCAGCTGGTACACGGCANNN
```
Run a genome-wide SpCas9 (NRG) search:
```bash
printf "hg38_1000G2021_HGDP\n"               > vcf_1000G.txt
printf "hg38_1000G2021_HGDP.samplesID.txt\n" > sid_1000G.txt

crisprme.py complete-search \
  --genome Genomes/hg38 \
  --vcf vcf_1000G.txt --samplesID sid_1000G.txt \
  --guide guides.txt --pam PAMs/20bp-NRG-SpCas9.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 6 --bDNA 2 --bRNA 2 --max-total-edits 6 \
  --output trac_1000G_HGDP --thread 16
```
(Substitute your actual annotation file names — `ls Annotations/` to see what downloaded.)

---

## 3. Search #2 — mega index (putative haplotypes)

Same guide, pointed at the mega index:
```bash
printf "hg38_mega\n"               > vcf_mega.txt
printf "hg38_mega.samplesID.txt\n" > sid_mega.txt

crisprme.py complete-search \
  --genome Genomes/hg38 \
  --vcf vcf_mega.txt --samplesID sid_mega.txt \
  --guide guides.txt --pam PAMs/20bp-NRG-SpCas9.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 6 --bDNA 2 --bRNA 2 --max-total-edits 6 \
  --output trac_mega --thread 16
```

---

## 4. What to check

For each run, open `Results/<output>/<jobid>_report.zip` → unzip → open **`report.html`** in a browser.
Also look at the bundled TSVs under the zip's `data/` folder (and in the `Results/<output>/` dir):

1. **The report renders** — summary, recommended validation panel, top-1000 table with annotations,
   and **co-occurrence sections** near the bottom.
2. **`snp_snp_cooc.tsv`** — off-targets needing ≥2 SNP alt alleles together.
   - Genotyped run: rows tagged **CONFIRMED** (phased cis) or **PUTATIVE**, with a **carrier sample list** + joint AF.
   - Mega run: rows tagged **PUTATIVE**, with a **`MinAF_bound`** and `NA` carriers.
3. **`indel_snp_cooc.tsv`** — SNP+indel co-occurrences, same CONFIRMED/PUTATIVE model.
4. **`indel_af.tsv`** (mega run) — each indel off-target's **per-dataset allele frequency**:
   columns `AF_1000G2021 AF_HGDP AF_gnomAD AF_TOPMed AF_AoU AF_max`. A `.` means "not reported by
   that source"; a number is that source's frequency. This is the **provenance** — which of the
   5 databases actually see each variant.
5. **`integrated_results.tsv`** — the main off-target table (aligned protospacer, mismatches/bulges,
   CFD + CRISTA scores, annotations, and the creating variant's rsID/AF/samples).

---

## 5. What "good" looks like

- ✅ Both searches complete and produce a `report.zip` (no crash / no error log).
- ✅ **Genotyped** report: CONFIRMED SNP+indel and SNP+SNP rows carry **named carrier samples** and an
  exact joint allele frequency; unphased ones are PUTATIVE.
- ✅ **Mega** report: SNP+indel and SNP+SNP rows are **PUTATIVE** with a **min-AF** bound and
  per-dataset provenance in `indel_af.tsv`; carrier columns read `NA` (correct — no genotypes).
- ✅ The two indices are **complementary**: the same locus may be CONFIRMED with carriers in the
  genotyped run and PUTATIVE (broader, 5-source frequency) in the mega run.

---

## 6. Feedback

Please note anything that looks off and send it back (or open an issue at
**github.com/pinellolab/crisprme-plus/issues**):
- crashes, error logs, or a search that never finishes;
- a missing or empty companion file (`snp_snp_cooc.tsv`, `indel_snp_cooc.tsv`, `indel_af.tsv`);
- report sections that are confusing or render wrong;
- allele frequencies / provenance / carrier lists that look implausible;
- anything in the docs (this file, the README, METHODS) that was unclear or didn't match reality.

Rough sanity numbers from our chr22 + genome-wide validation (so you know the ballpark): a
genome-wide NRG search on this guide produces thousands of off-targets and hundreds–thousands of
co-occurrence rows; the mega run reports **all** co-occurrences as PUTATIVE with min-AF bounds.
