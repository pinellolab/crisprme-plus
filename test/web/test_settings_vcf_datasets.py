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
   uploaded as a single multi-chromosome file, alike. Before its build, the search
   form offers a 0-bulge search of the first, which prepares the variant genome; of
   the second, which a search does not split per chromosome (the build does), it says
   to build it first, and Submit refuses that search with the same message.
2. The search form offers a dataset for the genome its marker names, whatever its
   name, with the bulges of its index; a genotyped dataset with no tiers reads
   "genotyped, samples not indexed", not "sites-only".
3. Submit refuses a variant search whose sample list is missing, instead of letting
   it fail at its end, and still runs one whose list is there; for a sites-only
   dataset it says why there is no list, not to build one (final review round 1, F4).
4. A sites-only dataset is built without a sample list, as the mega index is
   (docs/DESIGN_mega_index.md): a list of its pseudo-samples there, and the registry
   made from the VCFs' INFO allele frequencies (build_mega_registry.py), are left as
   they are (final review round 1, F1); and a list of only its header, which someone
   put there, survives an upload and a build (F3).

Standalone (no pytest); exits non-zero on failure:

    docker run --rm -v "$PWD/test/web:/t" crisprme:web python -u /t/test_settings_vcf_datasets.py
"""
import hashlib
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


def upload_vcf(dataset, fname, records, header=HEADER):
    """The chunked uploader's finished file, handed to _finalize_upload as the
    upload route does, with the reference genome picked on the card."""
    part = os.path.join(D, "Settings", "uploads", fname + ".part")
    bgzip(header + "\n".join(records) + "\n", part)
    return settings_page._finalize_upload("vcf", part, fname, dataset, "toyg")


def settings_build(dataset, bulges=1):
    """The command the Settings 'Build' button launches for dataset, run as the
    Settings job runner runs it. Returns (argv, rc, stderr)."""
    launched = []
    settings_page.launch_settings_job = lambda argv, stage: launched.append(argv) or "set_TEST000000"
    settings_page.build_index(1, "toyg", PAM, bulges, bulges, dataset, None)
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
    # before its build: a search prepares the variant genome from the VCFs, and does not
    # split one multi-chromosome file per chromosome (the build does), so that search
    # failed as it began; the form says to build it first, and Submit refuses it
    value = dataset[len("toyg_"):] if dataset.startswith("toyg_") else dataset
    note = main_page.update_bulge_note("toyg", PAM, value, 0, 0)
    out = main_page.check_input(1, None, "toyg", PAM, "IP", "CTAACAGTTGCTTTTATCAC",
                                4, 0, 0, value, False)
    msgs = [li.children for li in out[-1].children[1].children]
    if len(files) == 1:
        check("the dataset's VCF is a single multi-chromosome file" in note
              and "Build its index in Settings (Build an index locally) first" in note,
              f"{dataset}: before its build the note says to build it first (got {note!r})")
        check(out[0] is None and msgs == [note],
              f"{dataset}: and Submit refuses its search with the same message (got {msgs})")
    else:
        check("the first one prepares the variant genome" in note,
              f"{dataset}: before its build the note offers a 0-bulge search (got {note!r})")
        check(out[0] == 1 and not msgs, f"{dataset}: and Submit runs it ({msgs})")
    argv, rc, err = settings_build(dataset)
    check("--samplesID" in argv, f"{dataset}: the Settings build passes --samplesID ({argv})")
    check(rc == 0, f"{dataset}: the Settings build succeeds (rc {rc}; stderr {err[-1500:]!r})")
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
# a sites-only dataset: building it writes no list, so that is not the remedy
ok, msgs = submit("agg")
check(not ok and len(msgs) == 1 and msgs[0].startswith(
          "The variant dataset 'agg' is sites-only: its VCFs have no sample columns, so it "
          "has no sample list (samplesIDs/toyg_agg.samplesID.txt)")
      and "per-dataset allele-frequency registry" in msgs[0] and "Build its" not in msgs[0],
      f"Submit says why a sites-only dataset has no list (got {msgs})")
# and names the index with a registry as a web user gets it, not a file of the repository
# (it cited docs/DESIGN_mega_index.md: FV :342)
check(not ok and len(msgs) == 1
      and "(Settings downloads it, under Download a prebuilt index)" in msgs[0]
      and ".md" not in msgs[0],
      f"Submit's refusal of a sites-only dataset says what a web user can do (got {msgs})")

# ---- a list of only its header, put there for a sites-only dataset ----------------------
SITES = HEADER.replace("\tFORMAT\t" + "\t".join(SAMPLES), "")
sid = os.path.join(D, "samplesIDs", "toyg_aggh.samplesID.txt")
with open(sid, "w") as fh:
    fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
upload_vcf("toyg_aggh", "aggh.vcf.gz", ["\t".join(r.split("\t")[:8]) for r in ALL], header=SITES)
check(os.path.isfile(sid), "toyg_aggh: a list of only its header survives an upload of the dataset")
argv, rc, err = settings_build("toyg_aggh")
check("--samplesID" not in argv and rc == 0,
      f"toyg_aggh: built without a sample list (rc {rc}; {argv})")
check(os.path.isfile(sid) and open(sid).read() == "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n",
      "toyg_aggh: and survives its build")

# ---- a sites-only dataset built as the mega index is ----------------------------------
# build-index-only without --samplesID, the registry made from the VCFs' INFO allele
# frequencies dropped in, and the list the shipped index installs: one pseudo-sample a
# source dataset. A Settings build passed that list, its tier emission rewrote the
# registry from genotypes the VCFs do not have, and every variant off-target and
# allele frequency of a search was gone.
MEGA = "toyg_mega"
mdir = os.path.join(D, "VCFs", MEGA)
os.makedirs(mdir)
mheader = HEADER.replace("\tFORMAT\t" + "\t".join(SAMPLES), "").replace(
    "##FORMAT", '##INFO=<ID=AF_1000G2021,Number=A,Type=Float,Description="AF">\n'
    '##INFO=<ID=AF_HGDP,Number=A,Type=Float,Description="AF">\n'
    '##INFO=<ID=AF_max,Number=A,Type=Float,Description="AF">\n##FORMAT')
for c in GENOME:
    rows = []
    for rec in RECORDS[c]:
        a, b = round(random.uniform(0.05, 0.45), 2), round(random.uniform(0.05, 0.45), 2)
        rows.append("\t".join(rec.split("\t")[:7]
                              + [f"AF_1000G2021={a};AF_HGDP={b};AF_max={max(a, b)}"]))
    path = os.path.join(mdir, f"mega.{c}.vcf.gz")
    bgzip(mheader + "\n".join(rows) + "\n", path)
    subprocess.run(["tabix", "-p", "vcf", path], check=True)
recipe = subprocess.run(["crisprme.py", "build-index-only", "--genome", os.path.join(D, "Genomes", "toyg"),
                         "--pam", os.path.join(D, "PAMs", PAM + ".txt"), "--bDNA", "1", "--bRNA", "1",
                         "--vcf", mdir, "--path", D + "/"], cwd=D, capture_output=True, text=True,
                        timeout=900)
check(recipe.returncode == 0, f"{MEGA}: built as the mega recipe does (rc {recipe.returncode}; "
      f"stderr {recipe.stderr[-800:]!r})")
reg = os.path.join(D, "Dictionaries", f"registry_{MEGA}")
for c in GENOME:
    subprocess.run([sys.executable, "/opt/conda/opt/crisprme/PostProcess/build_mega_registry.py",
                    "--vcf", os.path.join(mdir, f"mega.{c}.vcf.gz"), "--chrom", c,
                    "--out-dir", reg, "--datasets", "1000G2021,HGDP"],
                   check=True, capture_output=True)
MEGA_LIST = "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n" + "".join(
    f"{d}\t{d}\t{d}\tNA\n" for d in ("1000G2021", "HGDP", "gnomAD", "TOPMed", "AoU"))
msid = os.path.join(D, "samplesIDs", f"{MEGA}.samplesID.txt")
with open(msid, "w") as fh:
    fh.write(MEGA_LIST)


def md5s(folder):
    return {f: hashlib.md5(open(os.path.join(folder, f), "rb").read()).hexdigest()
            for f in sorted(os.listdir(folder))}


before, dt_before = md5s(reg), pages_utils.variant_dataset_data_type("toyg", "mega")
check(sorted(before) == sorted(f"reg_{c}.{e}" for c in GENOME for e in ("bin", "idx")),
      f"{MEGA}: its INFO-AF registry is in place ({sorted(before)})")
argv, rc, err = settings_build(MEGA, bulges=2)
check("--samplesID" not in argv, f"{MEGA}: the Settings build passes no --samplesID ({argv})")
check(rc == 0, f"{MEGA}: the Settings build succeeds (rc {rc}; stderr {err[-1500:]!r})")
check(md5s(reg) == before, f"{MEGA}: its registry is left as it was ({md5s(reg)} != {before})")
check(not os.path.exists(os.path.join(D, "Dictionaries", f"genotypes_{MEGA}")),
      f"{MEGA}: no genotype store is made for it")
check(open(msid).read() == MEGA_LIST, f"{MEGA}: its list is left as it was")
dt = pages_utils.variant_dataset_data_type("toyg", "mega")
check(dt == dt_before == "sites-only", f"{MEGA}: still reads sites-only ({dt_before} -> {dt})")
check(not pages_utils.variant_dataset_has_genotypes("toyg", "mega"),
      f"{MEGA}: no per-sample analysis is offered for it")

shutil.rmtree(D, ignore_errors=True)
print(f"\n{'FAILED' if _fails else 'OK'}: {len(_fails)} failure(s)")
sys.stdout.flush()
sys.stderr.flush()
# force-exit: importing the Dash app can leave non-daemon helper threads that would
# otherwise stall interpreter shutdown (as test_variant_dropdown.py does)
os._exit(1 if _fails else 0)
