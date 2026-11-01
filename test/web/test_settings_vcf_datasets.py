#!/usr/bin/env python3
"""A VCF dataset added in the web Settings page is built and searchable from the
browser (Gate B walk P3, P4).

A web variant search reads samplesIDs/<dataset>.samplesID.txt. Nothing in Settings
wrote one and its build never passed --samplesID, so a dataset added and built only
there ran a whole search and failed at its end ("samples IDs processing failed"),
read "sites-only" on the search form, and could not run per-sample. Runs inside the
image: the Settings page's own helpers, and the real build-index-only (CRISPRitz,
tabix, bcftools) on a small synthetic genome, run with the argv the page launches.

1. The upload writes the dataset's reference-genome marker and its sample list,
   from the VCF header; the Settings build passes that list (--samplesID), so the
   dict-less tiers are built and the dataset reads "genotyped, phased", with
   per-sample analysis available. A dataset uploaded per chromosome, and one
   uploaded as a single multi-chromosome file, alike.
2. The search form offers a dataset for the genome its marker names, whatever its
   name, with the bulges of its index; a genotyped dataset with no tiers reads
   "genotyped, samples not indexed", not "sites-only".
3. Submit refuses a variant search whose sample list is missing, instead of letting
   it fail at its end, and still runs one whose list is there.

Standalone (no pytest); exits non-zero on failure:

    docker run --rm -v "$PWD/test/web:/t" crisprme:web python -u /t/test_settings_vcf_datasets.py
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, "/opt/conda/opt/crisprme")

_fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL "), msg)
    sys.stdout.flush()
    if not cond:
        _fails.append(msg)


from pages import pages_utils  # noqa: E402
from pages import settings_page  # noqa: E402
from pages import main_page  # noqa: E402

D = tempfile.mkdtemp(prefix="crisprme_setvcf_")
for sub in ("Genomes/toyg", "PAMs", "VCFs", "samplesIDs", "Settings/uploads", "Annotations"):
    os.makedirs(os.path.join(D, sub))
for mod in (pages_utils, settings_page, main_page):
    mod.current_working_directory = D + "/"
PAM = "20bp-NRG-SpCas9"
with open(os.path.join(D, "PAMs", PAM + ".txt"), "w") as fh:
    fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")

# a 3 x 3 kb genome and a phased VCF of 6 samples, 4 SNPs a chromosome
random.seed(281)
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]
GENOME = {c: "".join(random.choice("ACGT") for _ in range(3000)) for c in ("chr1", "chr2", "chr3")}
for c, seq in GENOME.items():
    with open(os.path.join(D, "Genomes", "toyg", c + ".fa"), "w") as fh:
        fh.write(f">{c}\n{seq}\n")
RECORDS = {c: [] for c in GENOME}
for c, seq in GENOME.items():
    for pos in (500, 1200, 1900, 2600):
        ref = seq[pos - 1]
        alt = "ACGT"[("ACGT".index(ref) + 1) % 4]
        gts = "\t".join(random.choice(["0|0", "0|1", "1|0", "1|1"]) for _ in SAMPLES)
        RECORDS[c].append(f"{c}\t{pos}\trs_{c}_{pos}\t{ref}\t{alt}\t.\tPASS\tAF=0.5\tGT\t{gts}")
HEADER = ("##fileformat=VCFv4.2\n"
          + "".join(f"##contig=<ID={c},length=3000>\n" for c in GENOME)
          + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n")


def bgzip(text, path):
    with open(path, "wb") as out:
        subprocess.run(["bgzip", "-c"], input=text.encode(), stdout=out, check=True)


def upload_vcf(dataset, fname, records):
    """The chunked uploader's finished file, handed to _finalize_upload as the
    upload route does, with the reference genome picked on the card."""
    part = os.path.join(D, "Settings", "uploads", fname + ".part")
    bgzip(HEADER + "\n".join(records) + "\n", part)
    return settings_page._finalize_upload("vcf", part, fname, dataset, "toyg")


def settings_build(dataset):
    """The command the Settings 'Build' button launches for dataset, run as the
    Settings job runner runs it. Returns (argv, rc, stderr)."""
    launched = []
    settings_page.launch_settings_job = lambda argv, stage: launched.append(argv) or "set_TEST000000"
    settings_page.build_index(1, "toyg", PAM, 1, 1, dataset, None)
    argv = launched[-1]
    done = subprocess.run(["crisprme.py"] + argv + ["--path", D + "/"], cwd=D,
                          capture_output=True, text=True, timeout=900)
    return argv, done.returncode, done.stderr


ALL = [r for c in GENOME for r in RECORDS[c]]
UPLOADS = {
    "toyg_walk": [(f"walk.{c}.vcf.gz", RECORDS[c]) for c in GENOME],  # per chromosome
    "walkvcf": [("walk.vcf.gz", ALL)],  # one multi-chromosome file, named like my_pig_vcf
}
for dataset, files in UPLOADS.items():
    for fname, records in files:
        upload_vcf(dataset, fname, records)
    marker = os.path.join(D, "VCFs", f".{dataset}.refgenome")
    check(os.path.isfile(marker) and open(marker).read() == "toyg",
          f"{dataset}: the upload records its reference genome")
    sid = os.path.join(D, "samplesIDs", f"{dataset}.samplesID.txt")
    rows = open(sid).read().splitlines() if os.path.isfile(sid) else None
    check(rows == ["#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX"]
          + [f"{s}\t{dataset}\tunknown\tunknown" for s in SAMPLES],
          f"{dataset}: the upload writes its sample list from the VCF header "
          f"(got {rows and rows[:3]}...)")
    argv, rc, err = settings_build(dataset)
    check("--samplesID" in argv, f"{dataset}: the Settings build passes --samplesID ({argv})")
    check(rc == 0, f"{dataset}: the Settings build succeeds (rc {rc}; stderr {err[-1500:]!r})")
    value = dataset[len("toyg_"):] if dataset.startswith("toyg_") else dataset
    dt = pages_utils.variant_dataset_data_type("toyg", value)
    check(dt == "genotyped-phased", f"{dataset}: reads genotyped-phased (got {dt})")

# ---- the search form ---------------------------------------------------------------------
opts = {o["value"]: o["label"] for o in pages_utils.get_variant_dataset_options("toyg")}
check("walk" in opts and "walkvcf" in opts,
      f"the search form offers both, whatever their name (got {opts})")
for value in ("walk", "walkvcf"):
    check("[genotyped, phased]" in opts.get(value, ""), f"{value} label ({opts.get(value)})")
    check(pages_utils.variant_dataset_has_genotypes("toyg", value),
          f"{value}: per-sample analysis is available (a genotype store)")
    check(pages_utils.index_max_bulges("toyg", PAM, value) == 1,
          f"{value}: its index's bulges are offered "
          f"(got {pages_utils.index_max_bulges('toyg', PAM, value)})")
# a genotyped dataset never built, and a sites-only one
for name, header in (("toyg_fresh", HEADER),
                     ("toyg_agg", HEADER.replace("\tFORMAT\t" + "\t".join(SAMPLES), ""))):
    os.makedirs(os.path.join(D, "VCFs", name))
    bgzip(header, os.path.join(D, "VCFs", name, f"{name}.chr1.vcf.gz"))
opts = {o["value"]: o["label"] for o in pages_utils.get_variant_dataset_options("toyg")}
check(opts.get("fresh", "").endswith("[genotyped, samples not indexed]"),
      f"a genotyped dataset with no tiers is not labelled sites-only (got {opts.get('fresh')})")
check(opts.get("agg", "").endswith("[sites-only (allele frequencies)]"),
      f"a sites-only dataset still is (got {opts.get('agg')})")


def submit(variant):
    """check_input, as the Submit button runs it: (runs, messages)."""
    out = main_page.check_input(1, None, "toyg", PAM, "IP", "CTAACAGTTGCTTTTATCAC",
                                4, 0, 0, variant, False)
    return out[0] == 1, [li.children for li in out[-1].children[1].children]


for value in ("walk", "walkvcf", "ref"):
    ok, msgs = submit(value)
    check(ok and not msgs, f"Submit runs a search on {value} ({msgs})")
ok, msgs = submit("fresh")
check(not ok and any("No sample list for the variant dataset 'fresh'" in m
                     and "samplesIDs/toyg_fresh.samplesID.txt" in m for m in msgs),
      f"Submit refuses one whose sample list is missing (got {msgs})")

shutil.rmtree(D, ignore_errors=True)
print(f"\n{'FAILED' if _fails else 'OK'}: {len(_fails)} failure(s)")
sys.stdout.flush()
sys.stderr.flush()
# force-exit: importing the Dash app can leave non-daemon helper threads that would
# otherwise stall interpreter shutdown (as test_variant_dropdown.py does)
os._exit(1 if _fails else 0)
