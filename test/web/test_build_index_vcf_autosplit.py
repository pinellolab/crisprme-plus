#!/usr/bin/env python3
"""build-index-only splits a single multi-chromosome .vcf.gz that has no .tbi
(Gate B walk P7).

The documented auto-split read the file by region (tabix -l, bcftools view -r),
which needs an index, so a .vcf.gz uploaded alone (the web Settings uploader takes
one file at a time) failed: "Could not load .tbi index ... Auto-splitting the VCF
dataset into per-chromosome files failed". It also left the dataset folder with the
0700 mode of the mkdtemp folder it came from, and let the tools write to stderr,
which fails a pipeline stage. Runs the real build-index-only in the image
(CRISPRitz, tabix, bcftools) on a small synthetic genome:

* bgzip with no .tbi, and plain gzip (which tabix cannot index): split, built, and
  nothing on stderr; the dataset folder gets the umask's mode; no scratch folder
  is left; the user's own folder (a registered server folder) is never written;
* with its .tbi: split as before;
* unsorted: the run fails, saying why and what to do, and the folder is left as it
  was.

Standalone (no pytest); exits non-zero on failure:

    docker run --rm -v "$PWD/test/web:/t" crisprme:web python -u /t/test_build_index_vcf_autosplit.py
"""
import gzip
import os
import random
import shutil
import stat
import subprocess
import sys
import tempfile

_fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL "), msg)
    sys.stdout.flush()
    if not cond:
        _fails.append(msg)


D = tempfile.mkdtemp(prefix="crisprme_autosplit_")
for sub in ("Genomes/toyg", "PAMs", "VCFs", "server"):
    os.makedirs(os.path.join(D, sub))
PAM = os.path.join(D, "PAMs", "20bp-NRG-SpCas9.txt")
with open(PAM, "w") as fh:
    fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
random.seed(7)
SAMPLES = ["S1", "S2", "S3", "S4"]
GENOME = {c: "".join(random.choice("ACGT") for _ in range(3000)) for c in ("chr1", "chr2", "chr3")}
for c, seq in GENOME.items():
    with open(os.path.join(D, "Genomes", "toyg", c + ".fa"), "w") as fh:
        fh.write(f">{c}\n{seq}\n")
RECORDS = []
for c, seq in GENOME.items():
    for pos in (600, 1500, 2400):
        ref = seq[pos - 1]
        alt = "ACGT"[("ACGT".index(ref) + 1) % 4]
        gts = "\t".join(random.choice(["0|0", "0|1", "1|0", "1|1"]) for _ in SAMPLES)
        RECORDS.append(f"{c}\t{pos}\trs_{c}_{pos}\t{ref}\t{alt}\t.\tPASS\tAF=0.5\tGT\t{gts}")
VCF = ("##fileformat=VCFv4.2\n"
       + "".join(f"##contig=<ID={c},length=3000>\n" for c in GENOME)
       + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
       + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n")


def bgzip(text, path):
    with open(path, "wb") as out:
        subprocess.run(["bgzip", "-c"], input=text.encode(), stdout=out, check=True)


def plain_gzip(text, path):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def build(dataset):
    """build-index-only with the dataset's sample list, as the Settings page runs it
    (without one, its INCOMPLETE VARIANT INDEX warning is on stderr)."""
    os.makedirs(os.path.join(D, "samplesIDs"), exist_ok=True)
    with open(os.path.join(D, "samplesIDs", f"{dataset}.samplesID.txt"), "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                 + "".join(f"{s}\tPOP\tSUP\tunknown\n" for s in SAMPLES))
    listing = os.path.join(D, f"samples_{dataset}.txt")
    with open(listing, "w") as fh:
        fh.write(f"{dataset}.samplesID.txt\n")
    done = subprocess.run(
        ["crisprme.py", "build-index-only", "--genome", os.path.join(D, "Genomes", "toyg"),
         "--pam", PAM, "--bDNA", "1", "--bRNA", "1", "--vcf", os.path.join(D, "VCFs", dataset),
         "--samplesID", listing, "--path", D], cwd=D, capture_output=True, text=True,
        timeout=900)
    return done.returncode, done.stdout, done.stderr


def tree(path):
    return sorted(os.listdir(path))


UMASK = os.umask(0)
os.umask(UMASK)
records = VCF + "\n".join(RECORDS) + "\n"
# a dataset uploaded as one file, without its .tbi (bgzip, then plain gzip); one with
# its .tbi; and one registered as a server folder (VCFs/<name> a link to it)
os.makedirs(os.path.join(D, "VCFs", "bgz"))
bgzip(records, os.path.join(D, "VCFs", "bgz", "bgz.vcf.gz"))
os.makedirs(os.path.join(D, "VCFs", "pgz"))
plain_gzip(records, os.path.join(D, "VCFs", "pgz", "pgz.vcf.gz"))
os.makedirs(os.path.join(D, "VCFs", "tbi"))
bgzip(records, os.path.join(D, "VCFs", "tbi", "tbi.vcf.gz"))
subprocess.run(["tabix", "-p", "vcf", os.path.join(D, "VCFs", "tbi", "tbi.vcf.gz")], check=True)
bgzip(records, os.path.join(D, "server", "srv.vcf.gz"))
os.symlink(os.path.join(D, "server"), os.path.join(D, "VCFs", "srv"))
server_before = tree(os.path.join(D, "server"))
for dataset in ("bgz", "pgz", "tbi", "srv"):
    rc, out, err = build(dataset)
    vdir = os.path.join(D, "VCFs", dataset)
    check(rc == 0, f"{dataset}: build-index-only succeeds (rc {rc}; stderr {err[-1500:]!r})")
    check(err == "", f"{dataset}: nothing on stderr (got {err[-600:]!r})")
    check([f for f in tree(vdir) if f.endswith(".vcf.gz")]
          == [f"{dataset}.{c}.vcf.gz" for c in ("chr1", "chr2", "chr3")],
          f"{dataset}: split per chromosome ({tree(vdir)})")
    check(stat.S_IMODE(os.stat(vdir).st_mode) == 0o777 & ~UMASK,
          f"{dataset}: the split folder has the umask's mode, not mkdtemp's 0700 "
          f"(got {oct(stat.S_IMODE(os.stat(vdir).st_mode))})")
    check(os.path.isdir(os.path.join(D, "genome_library", f"20-NRG_2_toyg+{dataset}")),
          f"{dataset}: the variant index is built (named for its 20-nt guide and 3' PAM)")
    split = sum(1 for line in subprocess.run(
        f"zcat {vdir}/*.chr*.vcf.gz | grep -v '^#'", shell=True, capture_output=True,
        text=True).stdout.splitlines())
    check(split == len(RECORDS), f"{dataset}: every record kept ({split} of {len(RECORDS)})")
check(not [f for f in tree(os.path.join(D, "VCFs")) if f.startswith(".vcfsplit")],
      "no scratch folder left in VCFs/")
check(tree(os.path.join(D, "server")) == server_before,
      "the registered server folder is never written to")

# an unsorted file: fails, saying why and what to do, and is left as it was
os.makedirs(os.path.join(D, "VCFs", "uns"))
bgzip(VCF + "\n".join(reversed(RECORDS)) + "\n", os.path.join(D, "VCFs", "uns", "uns.vcf.gz"))
rc, out, err = build("uns")
check(rc == 1 and "Auto-splitting the VCF dataset into per-chromosome files failed" in err
      and "Unsorted positions" in err and "Sort it by position (bcftools sort" in err,
      f"unsorted: fails, saying why and what to do (rc {rc}; {err[-700:]!r})")
check(tree(os.path.join(D, "VCFs", "uns")) == ["uns.vcf.gz"], "unsorted: left as it was")
check(not [f for f in tree(os.path.join(D, "VCFs")) if f.startswith(".vcfsplit")],
      "unsorted: no scratch folder left in VCFs/")

# chromosomes named without 'chr' (CHROM 1 on a chr1 genome; a contig such as Contig1 or
# NC_000913.3): the enrichment takes a file's chromosome from the 'chr' token of its name,
# so the split gave <base>.1.vcf.gz and the build failed on "cannot open .../.fa" (wset N2).
# It fails before the split, saying why and what to do, and the folder is left as it was.
RULE = ("the enrichment (crispritz add-variants) takes a VCF's chromosome from the part of "
        "its file name that holds 'chr' (<name>.chr<N>.vcf.gz) and reads the genome's FASTA "
        "file of that name. Rename them in the VCF to the names of the genome's sequences "
        "(bcftools annotate --rename-chrs <old new file>, e.g. 1 to chr1), sequences named "
        "with 'chr' too, then add it again")
for dataset, rename, bad in (("numbered", {"chr1": "1", "chr2": "2", "chr3": "3"}, "1, 2, 3"),
                             ("dotted", {"chr1": "NC_000913.3"}, "NC_000913.3")):
    renamed = records
    for old, new in rename.items():
        renamed = renamed.replace(f"ID={old},", f"ID={new},").replace(f"\n{old}\t", f"\n{new}\t")
    os.makedirs(os.path.join(D, "VCFs", dataset))
    bgzip(renamed, os.path.join(D, "VCFs", dataset, f"{dataset}.vcf.gz"))
    rc, out, err = build(dataset)
    check(rc == 1 and err.strip().splitlines()[-1:] == [
        f"Error: The VCF dataset cannot be enriched: {dataset}.vcf.gz holds the chromosome(s) "
        f"{bad}, whose names do not hold 'chr': {RULE}"],
          f"{dataset}: fails at once, saying why and what to do (rc {rc}; {err[-900:]!r})")
# a file whose name holds 'chr' is taken for the file of that chromosome, which the genome
# must have a FASTA file of: one named so by chance (nochr.vcf.gz, of three chromosomes), and
# one of a chromosome another genome has (chr4 on a genome of chr1 to chr3)
os.makedirs(os.path.join(D, "VCFs", "nochr"))
bgzip(records, os.path.join(D, "VCFs", "nochr", "nochr.vcf.gz"))
os.makedirs(os.path.join(D, "VCFs", "other"))
bgzip(records.replace("\nchr1\t", "\nchr4\t").replace("ID=chr1,", "ID=chr4,"),
      os.path.join(D, "VCFs", "other", "other.chr4.vcf.gz"))
for dataset, name, chrom in (("nochr", "nochr.vcf.gz", "nochr"),
                             ("other", "other.chr4.vcf.gz", "chr4")):
    rc, out, err = build(dataset)
    check(rc == 1 and err.strip().splitlines()[-1:] == [
        f"Error: The VCF dataset cannot be enriched: the name of {name} makes it a VCF of "
        f"chromosome {chrom} (the enrichment, crispritz add-variants, takes the last part of a "
        f"VCF's file name that holds 'chr' for its chromosome), and the genome "
        f"{os.path.join(D, 'Genomes', 'toyg')} has no {chrom}.fa. Name each VCF after the "
        f"chromosome it holds (<name>.chr<N>.vcf.gz), as the genome's FASTA file of it is "
        f"named, or build the index with the genome the dataset is of"],
          f"{dataset}: fails at once, saying why and what to do (rc {rc}; {err[-900:]!r})")
for dataset in ("numbered", "dotted", "nochr", "other"):
    names = [f for f in tree(os.path.join(D, "VCFs", dataset)) if not f.endswith(".tbi")]
    check(len(names) == 1 and names[0].startswith(dataset), f"{dataset}: left as it was")
    check(not [n for n in os.listdir(os.path.join(D, "genome_library"))
               if n.endswith(f"_toyg+{dataset}")], f"{dataset}: no variant index")
    check(not os.path.isdir(os.path.join(D, "Genomes", f"toyg+{dataset}")),
          f"{dataset}: no enriched genome")
check(not [f for f in tree(os.path.join(D, "VCFs")) if f.startswith(".vcfsplit")],
      "no chr: no scratch folder left in VCFs/")

shutil.rmtree(D, ignore_errors=True)
print(f"\n{'FAILED' if _fails else 'OK'}: {len(_fails)} failure(s)")
sys.exit(1 if _fails else 0)
