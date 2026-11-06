#!/usr/bin/env python3
"""Which index serves which search, through the real CLI and the real engine.

    check_index_flex.py <work dir, which it creates>

A synthetic phased genome (3 x 30 kb) with sites planted for the example guide, built with
build-index-only, and the layouts people actually have: a download, a library the run cannot
write to, an index for another PAM, budget or guide length, a K=8 one. Each cell compares a
search with the one on an index dedicated to it, byte for byte, or checks that it refuses
clearly. test_index_resolution.py and friends pin the same decisions with a stub engine: this
runs the engine, which answers what it cannot do on stdout with exit 0. Run inside an image
with the tree under test in place. A few minutes on 2 threads. Exits non-zero on any failure.

CI runs it on the published image, which until a release is out can hold an engine older than
--prune-pam. Those K=8 cells then check the search refuses, saying so, as it must on that engine.
"""

import glob
import gzip
import json
import os
import random
import shutil
import subprocess
import sys

G = "CTAACAGTTGCTTTTATCAC"
PAMS = {
    "20bp-NRG-SpCas9.txt": "N" * 20 + "NRG 3\n",
    "20bp-NRR-SpCas9.txt": "N" * 20 + "NRR 3\n",
    "20bp-NGG-SpCas9.txt": "N" * 20 + "NGG 3\n",
    "20bp-NGA-SpCas9VQR.txt": "N" * 20 + "NGA 3\n",
    "21bp-NGG-X.txt": "N" * 21 + "NGG 3\n",
    "23bp-NGG-X.txt": "N" * 23 + "NGG 3\n",
    "20bp-NNNNNNNN-NO-PAM.txt": "N" * 28 + " 8\n",
    "20bp-NGCG-SpCas9VRER.txt": "N" * 20 + "NGCG 4\n",
    "24bp-NNNN-5p.txt": "N" * 24 + " -4\n",
    "24bp-NNNN-3p.txt": "N" * 24 + " 4\n",
    "24bp-VNNN-5p.txt": "VNNN" + "N" * 20 + " -4\n",
}
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]
THREADS = "2"
failed = []


def say(cell, ok, why=""):
    print(f"{'PASS' if ok else 'FAIL'}  {cell}{': ' + why if why and not ok else ''}", flush=True)
    if not ok:
        failed.append(cell)


def make_seed(seed):
    """The genome, VCFs, sample list, PAM files and guides of a working directory."""
    r = random.Random(20261102)
    chroms = {c: [r.choice("ACGT") for _ in range(30000)] for c in ("chr1", "chr2", "chr3")}

    def plant(chrom, at, site):
        chroms[chrom][at:at + len(site)] = list(site)

    plant("chr1", 2000, G + "TGG")                                    # the on-target
    plant("chr1", 6000, G[:9] + "A" + G[9:] + "AGG")                  # a DNA bulge
    plant("chr1", 10000, G[:6] + "T" + G[7:] + "CGG")                 # a SNP makes it perfect
    plant("chr2", 12000, G + "TCA")                                   # a SNP makes its PAM
    plant("chr3", 20000, "CTACCAGTTGCATTTCATCACTGG")                  # an indel and a SNP in cis
    plant("chr1", 25000, G + "TGCG")                                  # the on-target of NGCG
    fold = os.path.join(seed, "Genomes", "gx")
    os.makedirs(fold)
    for name, seq in chroms.items():
        with open(os.path.join(fold, f"{name}.fa"), "w") as fh:
            fh.write(f">{name}\n{''.join(seq)}\n")
    rows = {  # (POS, REF, ALT, the six genotypes)
        "chr1": [(10007, chroms["chr1"][10006], "G", "1|0 0|0 1|1 0|0 0|1 0|0")],
        "chr2": [(12023, chroms["chr2"][12022], "G", "0|0 1|0 0|0 0|1 0|0 0|0")],
        "chr3": [(20011, "CA", "C", "0|0 0|0 0|0 1|0 1|0 0|0"),
                 (20016, "C", "T", "0|0 0|0 0|0 1|0 0|1 0|0")],
    }
    vcfs = os.path.join(seed, "VCFs", "gx_ds")
    os.makedirs(vcfs)
    for name, recs in rows.items():
        with gzip.open(os.path.join(vcfs, f"ds.{name}.vcf.gz"), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n"
                     f"##contig=<ID={name},length=30000>\n"
                     '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(SAMPLES) + "\n")
            for pos, ref, alt, gts in recs:
                fh.write(f"{name}\t{pos}\trs{pos}\t{ref}\t{alt}\t.\tPASS\tAF=0.25\tGT\t"
                         + "\t".join(gts.split()) + "\n")
    os.makedirs(os.path.join(seed, "samplesIDs"))
    with open(os.path.join(seed, "samplesIDs", "gx_ds.samplesID.txt"), "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for s in SAMPLES:
            fh.write(f"{s}\tP1\tS1\tmale\n")
    os.makedirs(os.path.join(seed, "PAMs"))
    for name, line in PAMS.items():
        with open(os.path.join(seed, "PAMs", name), "w") as fh:
            fh.write(line)
    for name, text in (("list_vcf.txt", "gx_ds\n"), ("list_samplesID.txt", "gx_ds.samplesID.txt\n"),
                       ("g20.txt", G + "NNN\n"), ("g21.txt", "A" + G + "NNN\n"),
                       ("g23.txt", "ATA" + G + "NNN\n"), ("g20_4.txt", G + "NNNN\n"),
                       ("g20_5p.txt", "NNNN" + G + "\n")):
        with open(os.path.join(seed, name), "w") as fh:
            fh.write(text)
    for d in ("Results", "Dictionaries", "Annotations", "genome_library"):
        os.makedirs(os.path.join(seed, d), exist_ok=True)


def cli(work, *args, user=None, env=None):
    """Run crisprme.py in ``work``; (exit status, stdout + stderr)."""
    cmd = ["crisprme.py", *args]
    if user is not None:
        cmd = ["setpriv", f"--reuid={user}", f"--regid={user}", "--clear-groups", *cmd]
    done = subprocess.run(cmd, cwd=work, capture_output=True, text=True,
                          env=dict(os.environ, HOME=work, **(env or {})))
    return done.returncode, done.stdout + done.stderr


def engine_prunes():
    """Whether the installed CRISPRitz has --prune-pam, read as the search script reads it."""
    with open(shutil.which("crispritz.py")) as fh:
        return "--prune-pam" in fh.read()


def search(work, output, guide="g20.txt", pam="20bp-NRG-SpCas9.txt", bulges=(1, 1), variant=True,
           extra=(), samples=True, user=None, env=None, vcf="list_vcf.txt"):
    args = ["complete-search", "--genome", "Genomes/gx", "--guide", guide, "--pam", f"PAMs/{pam}",
            "--mm", "4", "--output", output, "--thread", THREADS]
    if bulges is not None:
        args += ["--bDNA", str(bulges[0]), "--bRNA", str(bulges[1])]
    if variant:
        args += ["--vcf", vcf] + (["--samplesID", "list_samplesID.txt"] if samples else [])
    return cli(work, *args, *extra, user=user, env=env)


def result_log(work, output, name="log_error.txt"):
    """What a run wrote to its error log (or to another log of its folder)."""
    try:
        with open(os.path.join(work, "Results", output, name)) as fh:
            return fh.read()
    except OSError:
        return ""


def sites_only(work, dataset):
    """Stamp a dataset's registry sites-only, as the published mega panel's is."""
    registry = os.path.join(work, "Dictionaries", f"registry_{dataset}")
    for name in os.listdir(registry):
        if name.startswith("reg_") and name.endswith(".idx"):
            with open(os.path.join(registry, name)) as fh:
                manifest = json.load(fh)
            manifest["data_type"] = "sites-only"
            with open(os.path.join(registry, name), "w") as fh:
                json.dump(manifest, fh)


def table(work, output):
    """The integrated results of a run, or None."""
    found = glob.glob(os.path.join(work, "Results", output, "*integrated_results.tsv"))
    if not found:
        return None
    with open(found[0]) as fh:
        return fh.read()


def rows_at(text, chrom, start, block="highest_CFD"):
    """The rows of an integrated results table at chrom:start (its ``block`` columns), as
    dicts of the block's column names without the block suffix."""
    if not text:
        return []
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    suffix = f"_({block})"
    found = []
    for line in lines[1:]:
        row = dict(zip(head, line.split("\t")))
        if row.get("Chromosome") == chrom and row.get(f"Start_coordinate{suffix}") == str(start):
            found.append({k[:-len(suffix)] if k.endswith(suffix) else k: v for k, v in row.items()})
    return found


def the_cis_site(text):
    """Whether the cis SNP+indel site at chr3:20000 is the row the truth gives: one mismatch
    with the deletion and the SNP on one haplotype, CFD 0.353, carrier S4 alone (S5 has the
    deletion and the SNP on its two haplotypes, in trans)."""
    rows = rows_at(text, "chr3", 20000)
    return len(rows) == 1 and (
        rows[0]["Aligned_protospacer+PAM_ALT"] == "CTAcCAGTTGCTTTTATCACTGG"
        and rows[0]["CFD_score"] == "0.353" and rows[0]["Variant_samples"] == "S4"
        and rows[0]["Variant_rsID"] == "rs20011,rs20016")


def without_carriers(text):
    """The rows of an integrated results table without their Variant_samples columns, sorted:
    the same for two runs that differ only in the carriers they name and the order of their rows."""
    if not text:
        return None
    lines = text.rstrip("\n").split("\n")
    keep = [i for i, name in enumerate(lines[0].split("\t"))
            if not name.startswith("Variant_samples")]
    return sorted([line.split("\t")[i] for i in keep] for line in lines[1:])


def carriers(text):
    """{(chromosome, start): the set of the Variant_samples of a row's three blocks} for the rows
    of an integrated results table that hold a variant ('NA' is a reference row)."""
    found = {}
    lines = (text or "").rstrip("\n").split("\n")
    head = lines[0].split("\t")
    for line in lines[1:]:
        row = dict(zip(head, line.split("\t")))
        named = {v for k, v in row.items() if k.startswith("Variant_samples")}
        if named != {"NA"}:
            found[(row["Chromosome"], row["Start_coordinate_(highest_CFD)"])] = named
    return found


def cooc(work, output):
    """The (phase, cis carriers) of the lines of a run's SNP+indel co-occurrence companions."""
    found = set()
    for path in glob.glob(os.path.join(work, "Results", output, "*.indel_snp_cooc.tsv")):
        with open(path) as fh:
            lines = fh.read().splitlines()
        for line in lines[1:]:
            row = dict(zip(lines[0].split("\t"), line.split("\t")))
            found.add((row["phase"], row["cis_samples"]))
    return found


def build(root, seed, name, pam, bulges, variant=True):
    work = os.path.join(root, "lib", name)
    shutil.copytree(seed, work, symlinks=True)
    args = ["build-index-only", "--genome", "Genomes/gx", "--pam", f"PAMs/{pam}", "--bDNA",
            str(bulges), "--bRNA", str(bulges), "--thread", THREADS]
    if variant:
        args += ["--vcf", "list_vcf.txt", "--samplesID", "list_samplesID.txt"]
    rc, out = cli(work, *args)
    if rc:
        sys.exit(f"build-index-only {name} failed:\n{out}")
    return work


def copy_of(root, work, cell):
    to = os.path.join(root, "cell", cell)
    shutil.copytree(work, to, symlinks=True)
    return to


def as_download(work):
    """What `download --what index` leaves: no VCFs, no enriched genome."""
    for d in ("VCFs/gx_ds", "Genomes/gx+gx_ds", "Genomes/gx+gx_ds_INDELS",
              "Dictionaries/dictionaries_gx_ds"):
        shutil.rmtree(os.path.join(work, d), ignore_errors=True)
    return work


def identical(a, b):
    return a is not None and a == b


def main(root):
    seed = os.path.join(root, "seed")
    make_seed(seed)
    nrg2 = build(root, seed, "NRG_2", "20bp-NRG-SpCas9.txt", 1)
    nrr2 = build(root, seed, "NRR_2", "20bp-NRR-SpCas9.txt", 1)
    ngg2 = build(root, seed, "NGG_2", "20bp-NGG-SpCas9.txt", 1)
    k8 = build(root, seed, "K8_3", "20bp-NNNNNNNN-NO-PAM.txt", 2)
    g21 = build(root, seed, "g21_N2", "21bp-NGG-X.txt", 1, variant=False)
    ref3 = build(root, seed, "ref_NGG_3", "20bp-NGG-SpCas9.txt", 2, variant=False)

    # the truths: dedicated indexes, with everything in place
    full = copy_of(root, nrg2, "truth")
    search(full, "b11")
    search(full, "b00", bulges=(0, 0))
    b11, b00 = table(full, "b11"), table(full, "b00")
    say("truth: the example finds its off-targets at 1/1 and 0/0",
        bool(b11 and b00) and b11.count("\n") > b00.count("\n"))
    say("truth: the cis SNP+indel site is one mismatch, CFD 0.353, carrier S4 alone",
        the_cis_site(b11), str(rows_at(b11, "chr3", 20000)))

    # a 0-bulge variant search of a download serves from the installed index
    w = as_download(copy_of(root, nrg2, "download_0_bulges"))
    rc, out = search(w, "o", bulges=(0, 0))
    say("0 bulges, no enriched genome: same results", rc == 0 and identical(table(w, "o"), b00), out[-600:])

    # a search with bulges of the same download is what it was
    w = as_download(copy_of(root, nrg2, "download_1_bulge"))
    rc, out = search(w, "o")
    say("1/1 on a download: same results", rc == 0 and identical(table(w, "o"), b11), out[-600:])

    # no --samplesID: the list the download installed is taken
    w = as_download(copy_of(root, nrg2, "no_samples_flag"))
    rc, out = search(w, "o", samples=False)
    say("no --samplesID: same results", rc == 0 and identical(table(w, "o"), b11), out[-600:])
    os.remove(os.path.join(w, "samplesIDs", "gx_ds.samplesID.txt"))
    rc, out = search(w, "o2", samples=False)
    say("no --samplesID and no list: same results, no KeyError",
        rc == 0 and identical(table(w, "o2"), b11), out[-600:])
    # refused, it leaves its output folder empty: the corrected command, same --output, runs
    # (NGA, which NRG does not hold: an NGG search of NRG runs, pruned)
    w = as_download(copy_of(root, nrg2, "no_samples_flag_refused"))
    rc, out = search(w, "o", pam="20bp-NGA-SpCas9VQR.txt", samples=False)
    left = sorted(os.listdir(os.path.join(w, "Results", "o")))
    rc2, out2 = search(w, "o", samples=False)
    say("no --samplesID, refused: its output folder left empty, and the corrected command runs",
        rc != 0 and "no variant index here serves this search" in out and left == []
        and rc2 == 0 and identical(table(w, "o"), b11), f"left {left}: " + (out + out2)[-800:])

    # another PAM or budget than the download holds: one clear ERROR, and nothing made (NGA,
    # which NRG does not hold)
    w = as_download(copy_of(root, nrg2, "download_other_pam"))
    before = sorted(os.listdir(os.path.join(w, "Genomes")))
    rc, out = search(w, "o", pam="20bp-NGA-SpCas9VQR.txt")
    say("another PAM: refused, naming the PAM and the indexes there",
        rc != 0 and "no variant index here serves this search" in out and "PAM NGA" in out
        and "20-NRG_2_gx+gx_ds" in out and "cannot stat" not in out and "[OK] gx_ds:" not in out, out[-800:])
    say("another PAM: nothing made", sorted(os.listdir(os.path.join(w, "Genomes"))) == before)
    # a PAM the download's own holds (NGG in NRG, NGG and NRG in NRR): pruned to it by
    # itself, as a dedicated index of that PAM answers it (2.9.0)
    ngg_truth = copy_of(root, ngg2, "truth_ngg")
    search(ngg_truth, "b11", pam="20bp-NGG-SpCas9.txt")
    for index, pam, truth in ((nrg2, "20bp-NGG-SpCas9.txt", table(ngg_truth, "b11")),
                              (nrr2, "20bp-NGG-SpCas9.txt", table(ngg_truth, "b11")),
                              (nrr2, "20bp-NRG-SpCas9.txt", b11)):
        w = as_download(copy_of(root, index, f"download_{os.path.basename(index)}_{pam[5:8]}"))
        rc, out = search(w, "o", pam=pam)
        cell = f"{pam[5:8]} on a {os.path.basename(index)} download"
        if engine_prunes():
            say(f"{cell}: pruned by itself, same results as a dedicated index",
                rc == 0 and "searching with --prune-pam, as if it was given" in out
                and identical(table(w, "o"), truth), out[-600:])
        else:
            say(f"{cell} on an engine without --prune-pam: refused, saying what it needs",
                rc != 0 and "which the installed CRISPRitz cannot do" in out, out[-800:])
    w = as_download(copy_of(root, k8, "download_k8"))
    rc, out = search(w, "o")
    if engine_prunes():
        say("K=8 without --prune-pam: pruned by itself, same results as a dedicated index",
            rc == 0 and "searching with --prune-pam, as if it was given" in out
            and identical(table(w, "o"), b11), out[-600:])
    else:
        say("K=8 without --prune-pam on an engine without it: refused, saying what it needs",
            rc != 0 and "which the installed CRISPRitz cannot do" in out, out[-800:])
    rc, out = search(w, "o2", extra=("--prune-pam", "--index-path", "genome_library"))
    if engine_prunes():
        say("K=8 with --prune-pam: same results as a dedicated index",
            rc == 0 and identical(table(w, "o2"), b11), out[-600:])
    else:
        say("K=8 with --prune-pam on an engine without it: refused, saying so", rc != 0
            and "does not support it" in result_log(w, "o2"), out[-600:])
    rc, out = search(w, "o3", extra=("--prune-pam",))
    if engine_prunes():
        say("K=8 with --prune-pam and no --index-path: same results as a dedicated index",
            rc == 0 and identical(table(w, "o3"), b11), out[-600:])

    # build-index-only writes the name of 2.9.0, which carries the guide length and the side
    # of the PAM (20-NRG_2_gx...: a 20-nt guide, a 3' PAM), the INDELS index's too. Under the
    # older name, as a download or an earlier release has it, the same results; under the
    # name of a 5' PAM, not taken for this 3' search
    for prefix, cell in (("NRG", "the older name"), ("NRG-20", "a 5' PAM's name")):
        w = as_download(copy_of(root, nrg2, f"download_named_{prefix}"))
        lib = os.path.join(w, "genome_library")
        for name in sorted(os.listdir(lib)):
            os.rename(os.path.join(lib, name), os.path.join(lib, prefix + name[6:]))
        rc, out = search(w, "o")
        if prefix == "NRG":
            say(f"an index under {cell} ({prefix}_2_gx...): same results",
                rc == 0 and identical(table(w, "o"), b11), out[-600:])
        else:
            say(f"an index under {cell} ({prefix}_2_gx...): not taken for a 3' PAM",
                rc != 0 and "no variant index here serves this search" in out, out[-800:])

    # a PAM-less index under the 2.9.0 name with no .pam_build: the name records its side as
    # the file does (an older name with none is a 3' build at widths 3 and 8 only), so under
    # --index-path a 3'-named one of width 4 serves NGG, pruned by itself, and a 5'-named one a
    # 5' PAM of its width; both were refused
    ngg_ref = copy_of(root, ngg2, "truth_ngg_reference")
    search(ngg_ref, "o", pam="20bp-NGG-SpCas9.txt", variant=False)
    k4_3 = build(root, seed, "K4_3p", "24bp-NNNN-3p.txt", 1, variant=False)
    w = copy_of(root, k4_3, "named_3p_unrecorded")
    lib = os.path.join(w, "genome_library")
    os.remove(os.path.join(lib, "20-NNNN_2_gx", ".pam_build"))
    rc, out = search(w, "o", pam="20bp-NGG-SpCas9.txt", variant=False,
                     extra=("--index-path", "genome_library"))
    if engine_prunes():
        say("a 3'-named PAM-less index of width 4 with no .pam_build (20-NNNN_2_gx): NGG pruned "
            "by itself, same results as a dedicated index",
            rc == 0 and identical(table(w, "o"), table(ngg_ref, "o")), out[-600:])
    k4_5 = build(root, seed, "K4_5p", "24bp-NNNN-5p.txt", 1, variant=False)
    truth = copy_of(root, k4_5, "truth_vnnn")
    search(truth, "o", guide="g20_5p.txt", pam="24bp-VNNN-5p.txt", variant=False,
           extra=("--index-path", "genome_library"))
    w = copy_of(root, k4_5, "named_5p_unrecorded")
    lib = os.path.join(w, "genome_library")
    os.remove(os.path.join(lib, "NNNN-20_2_gx", ".pam_build"))
    rc, out = search(w, "o", guide="g20_5p.txt", pam="24bp-VNNN-5p.txt", variant=False,
                     extra=("--index-path", "genome_library"))
    say("a 5'-named PAM-less index with no .pam_build (NNNN-20_2_gx): a 5' PAM of its width, "
        "same results as the index with its record",
        rc == 0 and identical(table(w, "o"), table(truth, "o")), out[-600:])

    # E2: the indexes of one motif for the two sides of the guide, and for two guide lengths,
    # built in one library, each searched as it would be alone. Under the older name, which
    # says neither, the second build of each pair stopped on the first ("... is already there
    # ... The index for that PAM file has the same name: build it in another working directory")
    # and a search built its own over nothing it could use.
    w = copy_of(root, seed, "both_sides_one_library")
    rc3, out3 = cli(w, "build-index-only", "--genome", "Genomes/gx", "--pam", "PAMs/24bp-NNNN-3p.txt",
                    "--bDNA", "1", "--bRNA", "1", "--thread", THREADS)
    rc5, out5 = cli(w, "build-index-only", "--genome", "Genomes/gx", "--pam", "PAMs/24bp-NNNN-5p.txt",
                    "--bDNA", "1", "--bRNA", "1", "--thread", THREADS)
    built = sorted(os.listdir(os.path.join(w, "genome_library")))
    say("E2: build-index-only of a 3' and a 5' NNNN file in one library: both built, each under "
        "its name", rc3 == 0 and rc5 == 0 and built == ["20-NNNN_2_gx", "NNNN-20_2_gx"],
        f"{built} " + (out3 + out5)[-800:])
    alone3 = copy_of(root, k4_3, "both_sides_truth_3p")
    alone5 = copy_of(root, k4_5, "both_sides_truth_5p")
    for cell, guide, pam, alone, index in (
            ("3'", "g20_4.txt", "24bp-NNNN-3p.txt", alone3, "20-NNNN_2_gx"),
            ("5'", "g20_5p.txt", "24bp-NNNN-5p.txt", alone5, "NNNN-20_2_gx"),
            ("5' TTTV-like (VNNN)", "g20_5p.txt", "24bp-VNNN-5p.txt", alone5, "NNNN-20_2_gx")):
        out = f"o_{pam[5:9]}_{cell[0]}"
        search(alone, out, guide=guide, pam=pam, variant=False)
        rc, said = search(w, out, guide=guide, pam=pam, variant=False)
        say(f"E2: a {cell} search of that library: the results of its index alone, from {index}",
            rc == 0 and identical(table(w, out), table(alone, out))
            and f"Reference_index\t{index} (in " in result_log(w, out, ".indexes_searched.txt")
            and sorted(os.listdir(os.path.join(w, "genome_library"))) == built, said[-600:])
    # two guide lengths of NGG, with their variant indexes: 20 nt and 23 nt
    w = copy_of(root, seed, "two_lengths_one_library")
    for pam in ("20bp-NGG-SpCas9.txt", "23bp-NGG-X.txt"):
        rc, out = cli(w, "build-index-only", "--genome", "Genomes/gx", "--pam", f"PAMs/{pam}",
                      "--bDNA", "1", "--bRNA", "1", "--thread", THREADS, "--vcf", "list_vcf.txt",
                      "--samplesID", "list_samplesID.txt")
        say(f"E2: build-index-only of {pam} beside the other length: built", rc == 0, out[-600:])
    built = sorted(os.listdir(os.path.join(w, "genome_library")))
    say("E2: the two lengths' indexes, each under its name, the INDELS ones too",
        built == ["20-NGG_2_gx", "20-NGG_2_gx+gx_ds", "20-NGG_2_gx+gx_ds_INDELS", "23-NGG_2_gx",
                  "23-NGG_2_gx+gx_ds", "23-NGG_2_gx+gx_ds_INDELS"], str(built))
    alone23 = copy_of(root, build(root, seed, "NGG_2_23nt", "23bp-NGG-X.txt", 1), "truth_23nt")
    search(alone23, "b11", guide="g23.txt", pam="23bp-NGG-X.txt")
    for guide, pam, alone, index in (("g20.txt", "20bp-NGG-SpCas9.txt", ngg_truth, "20-NGG_2_gx"),
                                     ("g23.txt", "23bp-NGG-X.txt", alone23, "23-NGG_2_gx")):
        out = f"o{guide[1:3]}"
        rc, said = search(w, out, guide=guide, pam=pam)
        searched = result_log(w, out, ".indexes_searched.txt")
        say(f"E2: a {guide[1:3]}-nt variant search of that library: the results of its index "
            f"alone, from {index} and {index}+gx_ds",
            rc == 0 and identical(table(w, out), table(alone, "b11"))
            and f"Reference_index\t{index} (in " in searched
            and f"Variant_index:gx_ds\t{index}+gx_ds (in " in searched
            and f"its indels: {index}+gx_ds_INDELS" in searched
            and sorted(os.listdir(os.path.join(w, "genome_library"))) == built, said[-600:])
    # an index of the older name, as 2.8.1 built it, for the 3' side: a 5' build is put beside
    # it, and each side's search takes its own
    w = copy_of(root, k4_3, "older_name_beside")
    lib = os.path.join(w, "genome_library")
    os.rename(os.path.join(lib, "20-NNNN_2_gx"), os.path.join(lib, "NNNN_2_gx"))
    rc, out = cli(w, "build-index-only", "--genome", "Genomes/gx", "--pam", "PAMs/24bp-NNNN-5p.txt",
                  "--bDNA", "1", "--bRNA", "1", "--thread", THREADS)
    say("E2: a 5' build beside a 3' index of the older name: built beside it, saying so",
        rc == 0 and sorted(os.listdir(lib)) == ["NNNN-20_2_gx", "NNNN_2_gx"]
        and "NNNN_2_gx is an index for a PAM on the 3' side of the guide (as its .pam_build "
            "records): it is left as it is, and this index is built beside it, as NNNN-20_2_gx."
        in out, out[-600:])
    for cell, guide, pam, alone, index in (("3'", "g20_4.txt", "24bp-NNNN-3p.txt", alone3, "NNNN_2_gx"),
                                           ("5'", "g20_5p.txt", "24bp-NNNN-5p.txt", alone5, "NNNN-20_2_gx")):
        out = f"o_{cell[0]}"
        rc, said = search(w, out, guide=guide, pam=pam, variant=False)
        say(f"E2: beside the older name, a {cell} search: the results of its index alone, from {index}",
            rc == 0 and identical(table(w, out), table(alone, f"o_{pam[5:9]}_{cell[0]}"))
            and f"Reference_index\t{index} (in " in result_log(w, out, ".indexes_searched.txt"),
            said[-600:])

    # a search that enriches a dataset itself codes the fake-indel genome
    w = copy_of(root, ngg2, "on_demand")
    for d in ("Genomes/gx+gx_ds", "Genomes/gx+gx_ds_INDELS"):
        shutil.rmtree(os.path.join(w, d))
    rc, out = search(w, "o")
    say("on-demand enrichment: same results as a built index",
        rc == 0 and identical(table(w, "o"), b11), out[-600:])

    # the first search of a dataset that has only its VCFs, as the web's first 0-bulge one, which
    # enriches it: it finds the target that needs the indel and the SNP, each carrier's own
    # haplotype resolved (the indel post-analysis, 2.8.1), as on a built index
    w = copy_of(root, seed, "first_search")
    rc, out = search(w, "o", bulges=(0, 0), env={"CRISPRME_ALLOW_ONDEMAND_BUILD": "1"})
    say("first search of a dataset with only its VCFs: the cis SNP+indel site at CFD 0.353 with "
        "carrier S4 alone, and the genotype store",
        rc == 0 and the_cis_site(table(w, "o"))
        and os.path.isfile(os.path.join(w, "Dictionaries", "indel_genotypes_gx_ds",
                                        "gt_indel_chr3.tsv.gz")),
        f"{rows_at(table(w, 'o'), 'chr3', 20000)} " + out[-600:])

    # the guide length: served inside the rule, built outside it
    dedicated = copy_of(root, g21, "g21_truth")
    search(dedicated, "o", guide="g21.txt", pam="21bp-NGG-X.txt", variant=False)
    w = copy_of(root, ref3, "g21_on_N3")
    rc, out = search(w, "o", guide="g21.txt", pam="21bp-NGG-X.txt", variant=False)
    say("a 21-nt guide at 1/1 on a 20-nt N=3 index: same results as a dedicated one",
        rc == 0 and identical(table(w, "o"), table(dedicated, "o")), out[-600:])
    # (a 23-nt guide beside a 20-nt index, of the name a 2.8.1 build gave it, which the search
    # built its own under: it stopped on it, as the 23-nt index took that name too)
    w = copy_of(root, ngg2, "g23_on_N2")
    lib = os.path.join(w, "genome_library")
    for name in sorted(os.listdir(lib)):
        os.rename(os.path.join(lib, name), os.path.join(lib, name[3:]))
    record = os.path.join(lib, "NGG_2_gx", ".pam_build")
    with open(record) as fh:
        before = fh.read()
    dedicated = copy_of(root, seed, "g23_truth")
    search(dedicated, "o", guide="g23.txt", pam="23bp-NGG-X.txt", variant=False)
    rc, out = search(w, "o", guide="g23.txt", pam="23bp-NGG-X.txt", variant=False)
    say("a 23-nt guide beside a 20-nt index of the older name: its own built beside it, same "
        "results as a dedicated one", rc == 0 and identical(table(w, "o"), table(dedicated, "o"))
        and os.path.isdir(os.path.join(lib, "23-NGG_2_gx")), out[-800:])
    rc, out = cli(w, "build-index-only", "--genome", "Genomes/gx", "--pam", "PAMs/23bp-NGG-X.txt",
                  "--bDNA", "1", "--bRNA", "1", "--thread", THREADS)
    with open(record) as fh:
        say("build-index-only of the 23-nt PAM file there: the 23-nt index already present, the "
            "20-nt one's record kept", rc == 0 and "already present: " in out
            and "23-NGG_2_gx" in out and fh.read() == before, out[-600:])
    w = copy_of(root, seed, "g23_guide_on_20_pam")
    rc, out = search(w, "o", guide="g23.txt", pam="20bp-NGG-SpCas9.txt", variant=False)
    say("a guide longer than its PAM file's spacer: refused", rc != 0 and "silently drop" in out, out[-600:])

    # complete-test on a K=8 index, alone and beside the default; on K=8 it searches the
    # default PAM, NRR (NRG before 2.9.0)
    w = as_download(copy_of(root, k8, "complete_test_k8"))
    pams = sorted(os.listdir(os.path.join(w, "PAMs")))
    rc, out = cli(w, "complete-test", "--thread", THREADS)
    k8_table = table(w, "crisprme-test-out")
    n = as_download(copy_of(root, nrr2, "complete_test_nrr"))
    rc2, out2 = cli(n, "complete-test", "--thread", THREADS)
    if engine_prunes():
        say("complete-test on K=8: same results as on NRR, and no PAM file left",
            rc == 0 and rc2 == 0 and identical(k8_table, table(n, "crisprme-test-out"))
            and sorted(os.listdir(os.path.join(w, "PAMs"))) == pams, out[-800:])
    else:
        say("complete-test on K=8 on an engine without --prune-pam: refused, saying so, and no "
            "PAM file left", rc != 0 and "does not support it" in result_log(w, "crisprme-test-out")
            and sorted(os.listdir(os.path.join(w, "PAMs"))) == pams, out[-800:])
    both = as_download(copy_of(root, k8, "complete_test_k8_and_nrg"))
    for name in os.listdir(os.path.join(nrg2, "genome_library")):
        shutil.copytree(os.path.join(nrg2, "genome_library", name),
                        os.path.join(both, "genome_library", name), symlinks=True)
    rc, out = cli(both, "complete-test", "--chrom", "chr1", "--thread", THREADS)
    say("complete-test with K=8 beside NRG takes NRG", rc == 0 and "'20-NRG_2_gx+gx_ds'" in out, out[-600:])

    # a library the run cannot write to, as a user who may not write there
    ro = copy_of(root, ngg2, "read_only")
    for d, _, files in os.walk(ro):
        os.chmod(d, 0o755)
        for f in files:
            os.chmod(os.path.join(d, f), 0o644)
    uid = None
    if os.geteuid() == 0:  # root writes anywhere: the cells run as another user
        uid = 65534
        subprocess.run(["chown", "-R", f"{uid}:{uid}", ro], check=True)
        os.chmod(os.path.dirname(ro), 0o755)
        os.chmod(os.path.dirname(os.path.dirname(ro)), 0o755)
    for d, _, _files in os.walk(os.path.join(ro, "genome_library")):
        os.chmod(d, 0o555)
    rc, out = search(ro, "other_pam", variant=False, user=uid)
    say("read-only library lacking the index, with bulges: refused, not an empty result",
        rc != 0 and "is not installed, and this run cannot build it" in out, out[-800:])
    rc, out = search(ro, "same_pam", pam="20bp-NGG-SpCas9.txt", user=uid)
    w = copy_of(root, ngg2, "read_only_truth")
    search(w, "o", pam="20bp-NGG-SpCas9.txt")
    say("read-only library holding the index: same results",
        rc == 0 and identical(table(ro, "same_pam"), table(w, "o")), out[-600:])
    rc, out = search(ro, "no_bulges", bulges=(0, 0), variant=False, user=uid)
    w = copy_of(root, seed, "read_only_b0_truth")
    search(w, "o", bulges=(0, 0), variant=False)
    say("read-only library with no bulges: same results",
        rc == 0 and identical(table(ro, "no_bulges"), table(w, "o")), out[-600:])
    for d, _, _files in os.walk(os.path.join(ro, "genome_library")):
        os.chmod(d, 0o755)

    # a variant index in the --index-path library, the dataset's dictionaries here: read there
    w = as_download(copy_of(root, nrg2, "variant_library_elsewhere"))
    shutil.move(os.path.join(w, "genome_library"), os.path.join(w, "shared"))
    os.makedirs(os.path.join(w, "genome_library"))
    rc, out = search(w, "o", extra=("--index-path", "shared"))
    say("a variant index at --index-path: same results, nothing built here",
        rc == 0 and identical(table(w, "o"), b11)
        and os.listdir(os.path.join(w, "genome_library")) == [], out[-600:])

    # an index that lacks a chromosome's trees: refused, not searched as a whole one
    w = copy_of(root, nrg2, "index_lacking_a_chromosome")
    os.remove(os.path.join(w, "genome_library", "20-NRG_2_gx", "NRG_chr2_1.bin"))
    rc, out = search(w, "o", variant=False)
    say("a reference index lacking chr2's trees, with bulges: refused, naming it",
        rc != 0 and "lacks trees of FASTAs that hold sites of its PAM: chr2.fa" in out, out[-600:])

    # an install without the genotype store runs population-level, with no site lost
    w = as_download(copy_of(root, nrg2, "no_genotype_store"))
    for d in ("genotypes_gx_ds", "dictionaries_gx_ds"):
        shutil.rmtree(os.path.join(w, "Dictionaries", d), ignore_errors=True)
    rc, out = search(w, "o")
    rc2, out2 = search(w, "p", extra=("--population-level",))
    with open(os.path.join(w, "Results", "o", ".search_mode")) as fh:
        mode = fh.read()
    say("no genotype store: population-level, as --population-level gives it",
        rc == 0 and mode == "population-level" and identical(table(w, "o"), table(w, "p"))
        and "no genotype store" in out, out[-600:])
    # what it names, against a full install's --population-level run: the same rows but for their
    # carriers; a row of SNPs names none; the row of the cis site, which names the indel and the
    # SNP, lists every carrier of the indel, S4 and S5, as no genotype says which of them has the
    # SNP too (PUTATIVE), where a full install lists S4, whose haplotype holds both (CONFIRMED)
    search(full, "pop", extra=("--population-level",))
    mine, theirs = carriers(table(w, "o")), carriers(table(full, "pop"))
    cis = ("chr3", "20000")
    say("no genotype store: a full install's --population-level rows, no carrier on a row of SNPs, "
        "and every carrier of the indel at the cis SNP+indel site",
        without_carriers(table(w, "o")) is not None
        and without_carriers(table(w, "o")) == without_carriers(table(full, "pop"))
        and [r["Variant_rsID"] for r in rows_at(table(w, "o"), *cis)] == ["rs20011,rs20016"]
        and mine.get(cis) == {"S4,S5"} and theirs.get(cis) == {"S4"}
        and cooc(w, "o") == {("PUTATIVE", "NA")} and cooc(full, "pop") == {("CONFIRMED", "S4")}
        and len(mine) > 1
        and all(mine[k] == {""} and theirs.get(k, {""}) != {""} for k in mine if k != cis),
        f"no store {mine} {cooc(w, 'o')}, full {theirs} {cooc(full, 'pop')}")

    # a sites-only panel, as mega (its registry says so): --per-sample is the run without it,
    # and the WARNING says the panel has no genotypes, not to install a store it does not have
    w = as_download(copy_of(root, nrg2, "sites_only"))
    for d in ("genotypes_gx_ds", "dictionaries_gx_ds"):
        shutil.rmtree(os.path.join(w, "Dictionaries", d), ignore_errors=True)
    sites_only(w, "gx_ds")
    rc, out = search(w, "o", extra=("--per-sample",))
    rc2, out2 = search(w, "p")
    say("sites-only: --per-sample gives the results without it, saying the panel has no genotypes",
        rc == 0 and rc2 == 0 and identical(table(w, "o"), table(w, "p"))
        and "--per-sample has nothing to resolve here, as what it searches is sites-only, with no "
            "genotypes (gx_ds): the population-level analysis runs instead." in out
        and "install the store" not in out, out[-600:])

    # a sites-only copy of the panel beside it, with no store and no dictionaries, as mega beside
    # the default index: the phased dataset keeps its per-sample default, with the cis site's own
    # carrier, and the sites-only one's post-analysis is the population-level one
    w = as_download(copy_of(root, nrg2, "sites_only_beside_phased"))
    for t in ("registry", "log_indels"):
        shutil.copytree(os.path.join(w, "Dictionaries", f"{t}_gx_ds"),
                        os.path.join(w, "Dictionaries", f"{t}_gx_mega"))
    for suffix in ("", "_INDELS"):
        shutil.copytree(os.path.join(w, "genome_library", f"20-NRG_2_gx+gx_ds{suffix}"),
                        os.path.join(w, "genome_library", f"20-NRG_2_gx+gx_mega{suffix}"), symlinks=True)
    sites_only(w, "gx_mega")
    # its indel map names no carrier, as a sites-only VCF's (mega's) does
    logs = os.path.join(w, "Dictionaries", "log_indels_gx_mega")
    for name in os.listdir(logs):
        with gzip.open(os.path.join(logs, name), "rt") as fh:
            lines = fh.read().splitlines()
        with gzip.open(os.path.join(logs, name), "wt") as fh:
            for n, line in enumerate(lines):
                cells = line.split("\t")
                if n and len(cells) > 1:
                    cells[1] = ""
                fh.write("\t".join(cells) + "\n")
    with open(os.path.join(w, "list_vcf2.txt"), "w") as fh:
        fh.write("gx_ds\ngx_mega\n")
    rc, out = search(w, "o", vcf="list_vcf2.txt")
    mode = result_log(w, "o", ".search_mode")
    say("a sites-only panel beside the phased one: per-sample, the cis site's own carrier, and the "
        "sites-only panel's rows population-level",
        rc == 0 and mode == "per-sample" and result_log(w, "o", ".alt_mode") == "alt-alignments"
        and the_cis_site(table(w, "o"))
        and "Sites-only, with no genotypes to resolve: gx_mega." in out
        and "no genotype store" not in out
        and "Sites-only dataset on chr3: population-level post-analysis."
            in result_log(w, "o", "log_verbose.txt"),
        f"mode {mode}: {rows_at(table(w, 'o'), 'chr3', 20000)} " + out[-600:])

    # the NNNN_2 a 5' search of an earlier release leaves, under the older name with no record:
    # not taken for NGCG
    truth = copy_of(root, seed, "ngcg_truth")
    search(truth, "o", guide="g20_4.txt", pam="20bp-NGCG-SpCas9VRER.txt", variant=False)
    w = copy_of(root, seed, "unrecorded_5_prime_all_n")
    rc, out = search(w, "a", guide="g20_5p.txt", pam="24bp-NNNN-5p.txt", variant=False)
    os.rename(os.path.join(w, "genome_library", "NNNN-20_2_gx"), os.path.join(w, "genome_library", "NNNN_2_gx"))
    os.remove(os.path.join(w, "genome_library", "NNNN_2_gx", ".pam_build"))
    rc, out = search(w, "o", guide="g20_4.txt", pam="20bp-NGCG-SpCas9VRER.txt", variant=False)
    say("an unrecorded 5' all-N index of width 4: an NGCG search builds its own, same results",
        rc == 0 and identical(table(w, "o"), table(truth, "o"))
        and os.path.isdir(os.path.join(w, "genome_library", "20-NGCG_2_gx")), out[-600:])

    # build-index-only of a folder outside Genomes/ records it: another folder of its name is refused
    w = copy_of(root, seed, "outside_genomes")
    shutil.copytree(os.path.join(w, "Genomes", "gx"), os.path.join(w, "projA", "ref"))
    shutil.copytree(os.path.join(w, "Genomes", "gx"), os.path.join(w, "projB", "ref"))
    rc, out = cli(w, "build-index-only", "--genome", "projA/ref", "--pam", "PAMs/20bp-NRG-SpCas9.txt",
                  "--bDNA", "1", "--bRNA", "1", "--thread", THREADS)
    rc2, out2 = cli(w, "complete-search", "--genome", "projB/ref", "--guide", "g20.txt", "--pam",
                    "PAMs/20bp-NRG-SpCas9.txt", "--mm", "4", "--bDNA", "1", "--bRNA", "1",
                    "--output", "o", "--thread", THREADS)
    say("build-index-only of a folder outside Genomes/: recorded, and another folder of its name refused",
        rc == 0 and os.path.isfile(os.path.join(w, "genome_library", "20-NRG_2_ref", ".genome_source"))
        and rc2 != 0 and "holds the genome in" in out2, out2[-600:])

    print(f"{len(failed)} cell(s) failed" if failed else "INDEX FLEX PASS", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    work_root = os.path.abspath(sys.argv[1])
    os.makedirs(work_root)
    sys.exit(main(work_root))
