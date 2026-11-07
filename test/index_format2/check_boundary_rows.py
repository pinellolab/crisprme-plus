#!/usr/bin/env python3
"""Rows at contig ends and N runs, through the real CLI and the real engine.

    check_boundary_rows.py <work dir, which it creates>
    check_boundary_rows.py --compare <work dir of format 1> <work dir of format 2>

An index of format 2 (CRISPRitz 2.10.0 and later) holds the sites whose guide lies within a
few bases of a contig end or an N run, a PAM ending on a contig's last base, and on the minus
strand a PAM at a contig's first base (docs/INDEX_FORMAT.md, B1-B3); in a fake-indel genome,
whose windows are joined by one N, a site at a window's edge. Every step of CRISPRme that reads
the genome, a variant dictionary or the indel log around a row must hold at those rows, which
no index of format 1 ever returned. This builds a synthetic genome planted with such sites --
perfect, mismatched and bulged, on both strands, with SNPs and indels near them, phased over six
samples -- and runs complete-search on it (reference, and variant per-sample and
population-level, 4 mismatches and 1 + 1 bulges; NGG, and the 5' TTTV; and NGG under --prune-pam
off a PAM-less index of an 8-nt flank, in a working directory of its own), with the CRISPRitz on
PATH. It checks that every run ends well with nothing on stderr, and that every row of the
integrated results is the genome's own bases at its coordinates, in each of its three alignments
(highest CFD, fewest mismatches and bulges, highest CRISPR-Bulge score): a reference row its
reference bases, a variant row its reference bases but at its variants' positions (each a variant
of the dataset), an indel row a target within its chromosome. The --prune-pam search off the 8-nt
flank must give the rows of the dedicated NGG index (format 2; format 1 misses the PAMs within 5
nt of a contig end, IDX-4). Run inside an image with the tree under test in place, once with each
engine; --compare then checks that the format-2 results hold every format-1 row that is its
genome's bases (format 1 reports a minus-strand PAM at a contig's first base on the plus strand,
IDX-14), and that each row more lies at a boundary.

Exits non-zero on any failure. With CRISPRitz 2.9.1 (index format 1) four checks fail, as they
must: three searches report the minus-strand PAM at chrB's first base on the plus strand (IDX-14),
and so the dedicated index's rows differ from the 8-nt flank's by that row (measured: 12 passed,
4 failed).
"""

import glob
import gzip
import os
import random
import shutil
import subprocess
import sys

G = "CTAACAGTTGCTTTTATCAC"
TTTV_G = "TTTA" + "GACTTCAGATCGGTACCATG"  # a 5' TTTV site
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]
COMP = str.maketrans("ACGTN", "TGCAN")
# a row's three alignments, as the integrated results name them (resultIntegrator.py)
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
failed = []


def rc(seq):
    return seq.translate(COMP)[::-1]


def say(cell, ok, why=""):
    print(f"{'PASS' if ok else 'FAIL'}  {cell}{': ' + why if why and not ok else ''}", flush=True)
    if not ok:
        failed.append(cell)


def mutate(seq, at):
    """``seq`` with a mismatch at each position of ``at``."""
    s = list(seq)
    for i in at:
        s[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[i]]
    return "".join(s)


def genome():
    """{contig: sequence} and the VCF records [(contig, POS, REF, ALT, genotypes)]."""
    r = random.Random(20261106)
    rand = lambda n: [r.choice("ACGT") for _ in range(n)]  # noqa: E731
    seqs = {c: rand(3000) for c in ("chrA", "chrB", "chrC")}

    def plant(c, at, site):
        seqs[c][at:at + len(site)] = list(site)

    # chrA: perfect sites at the very ends, and on the minus strand a PAM at base 0
    plant("chrA", 0, G + "TGG")                       # + strand: guide from base 0, pad off
    plant("chrA", 3000 - 23, G + "AGG")               # + strand: PAM on the last base
    seqs["chrA"][1000:1100] = ["N"] * 100             # an N run
    plant("chrA", 1100, G + "CGG")                    # guide from the run's edge
    plant("chrA", 1000 - 23, rc(G + "TGG"))           # - strand: guide up to the run
    plant("chrA", 2000, G + "TGG")                    # interior control
    # chrB: a minus strand PAM at base 0 (B1), mismatched and bulged sites at the ends
    plant("chrB", 0, rc(mutate(G, [3]) + "TGG"))       # - strand, PAM at 0, 1 mismatch
    plant("chrB", 3000 - 24, G[:10] + "A" + G[10:] + "TGG")   # a DNA bulge, PAM on the last base
    seqs["chrB"][1500:1520] = ["N"] * 20
    plant("chrB", 1520, G[:7] + G[8:] + "AGG")         # an RNA bulge, guide from the run's edge
    plant("chrB", 1500 - 23, mutate(G, [0, 1]) + "CGG")  # 2 mismatches, PAM up to the run
    plant("chrB", 2500, G + "AGG")
    # chrC: sites 1 and 2 nt from N runs, a 5' TTTV site at each end, SNP/indel ground
    seqs["chrC"][500:520] = ["N"] * 20
    plant("chrC", 521, G + "TGG")                      # guide 1 nt from the run
    seqs["chrC"][1200:1210] = ["N"] * 10
    plant("chrC", 1212, rc(G + "AGG"))                 # - strand, 2 nt from the run
    plant("chrC", 0, TTTV_G)                           # 5' TTTV: PAM at base 0
    plant("chrC", 3000 - 24, TTTV_G)                   # and its guide on the last base
    plant("chrC", 2000, G + "TGG")
    seqs = {c: "".join(s) for c, s in seqs.items()}
    # SNPs in boundary sites (one makes a 1-mismatch site perfect, one makes a PAM) and
    # indels by them, phased
    recs = [
        ("chrB", 4 + 1, None, None, "1|0 0|0 1|1 0|0 0|1 0|0"),       # in the minus site at 0
        ("chrA", 3000 - 2, None, "C", "0|1 0|0 0|0 1|0 0|0 0|0"),     # in the last PAM
        ("chrA", 1100 + 5, None, None, "0|0 1|0 0|0 0|1 0|0 0|0"),    # by the N run
        ("chrB", 1500 - 23 + 1, None, None, "1|1 0|0 0|0 0|0 1|0 0|0"),
        ("chrC", 521 + 10, None, None, "0|0 0|1 1|0 0|0 0|0 0|0"),
        ("chrA", 2000 + 7, "indel+", None, "1|0 0|1 0|0 0|0 0|0 0|0"),   # an insertion
        ("chrA", 1100 + 30, "indel-", None, "0|0 0|0 1|0 0|1 0|0 0|0"),  # a deletion
        ("chrB", 2990, "indel-", None, "0|1 0|0 0|0 0|0 0|0 1|0"),       # one by the end
        ("chrC", 2000 + 15, "indel+", None, "0|0 0|0 0|0 1|1 0|0 0|0"),
    ]
    out = []
    for c, pos, kind, alt, gts in recs:
        ref = seqs[c][pos - 1]
        if kind == "indel+":
            ref, alt = ref, ref + "T"
        elif kind == "indel-":
            ref, alt = seqs[c][pos - 1:pos + 1], seqs[c][pos - 1]
        elif alt is None:
            alt = {"A": "G", "C": "T", "G": "A", "T": "C"}[ref]
        out.append((c, pos, ref, alt, gts))
    return seqs, out


def make(work):
    seqs, recs = genome()
    fold = os.path.join(work, "Genomes", "bx")
    os.makedirs(fold)
    for name, seq in seqs.items():
        with open(os.path.join(fold, f"{name}.fa"), "w") as fh:
            fh.write(f">{name}\n{seq}\n")
    vcfs = os.path.join(work, "VCFs", "bx_ds")
    os.makedirs(vcfs)
    for name in seqs:
        with gzip.open(os.path.join(vcfs, f"ds.{name}.vcf.gz"), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n"
                     f"##contig=<ID={name},length={len(seqs[name])}>\n"
                     '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(SAMPLES) + "\n")
            for c, pos, ref, alt, gts in sorted(r for r in recs if r[0] == name):
                fh.write(f"{c}\t{pos}\trs{c[-1]}{pos}\t{ref}\t{alt}\t.\tPASS\tAF=0.25\tGT\t"
                         + "\t".join(gts.split()) + "\n")
    os.makedirs(os.path.join(work, "samplesIDs"))
    with open(os.path.join(work, "samplesIDs", "bx_ds.samplesID.txt"), "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for s in SAMPLES:
            fh.write(f"{s}\tP1\tS1\tmale\n")
    os.makedirs(os.path.join(work, "PAMs"))
    for name, line in (("20bp-NGG-SpCas9.txt", "N" * 20 + "NGG 3\n"),
                       ("20bp-TTTV-Cas12a.txt", "TTTV" + "N" * 20 + " -4\n"),
                       ("20bp-NNNNNNNN-NO-PAM.txt", "N" * 28 + " 8\n")):
        with open(os.path.join(work, "PAMs", name), "w") as fh:
            fh.write(line)
    for name, text in (("list_vcf.txt", "bx_ds\n"), ("list_samplesID.txt", "bx_ds.samplesID.txt\n"),
                       ("g20.txt", G + "NNN\n"), ("g20_5p.txt", "NNNN" + TTTV_G[4:] + "\n")):
        with open(os.path.join(work, name), "w") as fh:
            fh.write(text)
    for d in ("Results", "Dictionaries", "Annotations", "genome_library"):
        os.makedirs(os.path.join(work, d), exist_ok=True)
    return seqs, recs


def cli(work, *args):
    done = subprocess.run(["crisprme.py", *args], cwd=work, capture_output=True, text=True,
                          env=dict(os.environ, HOME=work))
    return done.returncode, done.stdout, done.stderr


SEARCHES = {  # output: (PAM file, guide file, variant, mode)
    "ref_ngg": ("20bp-NGG-SpCas9.txt", "g20.txt", False, None),
    "var_ngg_ps": ("20bp-NGG-SpCas9.txt", "g20.txt", True, "--per-sample"),
    "var_ngg_pl": ("20bp-NGG-SpCas9.txt", "g20.txt", True, "--population-level"),
    "ref_tttv": ("20bp-TTTV-Cas12a.txt", "g20_5p.txt", False, None),
    "var_tttv_ps": ("20bp-TTTV-Cas12a.txt", "g20_5p.txt", True, "--per-sample"),
}


def table(work, output):
    found = glob.glob(os.path.join(work, "Results", output, "*integrated_results.tsv"))
    if not found:
        return None, []
    with open(found[0]) as fh:
        lines = fh.read().rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return head, [dict(zip(head, ln.split("\t"))) for ln in lines[1:] if ln]


def wrong_rows(seqs, recs, rows):
    """The rows of a run, as (block, chromosome, start, strand, target, why), that are not
    their chromosome's bases at their coordinates (REF), or those but at the positions of their
    variants (ALT), each a variant of the dataset; or, an indel row, not within its chromosome.
    Each of a row's blocks is checked, and a row that lacks one is wrong."""
    snps = {(c, p): (ref, alt) for c, p, ref, alt, _g in recs if len(ref) == len(alt) == 1}
    bad = []
    for row in rows:
        for block in BLOCKS:
            start = row.get(f"Start_coordinate_({block})")
            if start is None:
                bad.append((block, row.get("Chromosome"), None, None, None, "no such block"))
                continue
            if start in ("", "NA", "n"):
                continue
            chrom, strand = row["Chromosome"], row[f"Strand_({block})"]
            ref_al = row.get(f"Aligned_protospacer+PAM_REF_({block})", "")
            alt_al = row.get(f"Aligned_protospacer+PAM_ALT_({block})", "")
            origin = row.get(f"REF/ALT_origin_({block})", "")
            info = row.get(f"Variant_info_genome_({block})", "")
            start = int(start)
            target = (alt_al if origin == "alt" and alt_al not in ("", "NA") else ref_al)
            bases = target.replace("-", "").upper()
            if not bases:
                continue
            if not 0 <= start <= len(seqs[chrom]) - len(bases):
                bad.append((block, chrom, start, strand, target, "off the chromosome"))
                continue
            if any(len(t.split("_")[-2]) != len(t.split("_")[-1])
                   for t in info.split(",") if t.count("_") >= 3):
                continue  # an indel row: its bases are the haplotype's, not the reference's
            genomic = seqs[chrom][start:start + len(bases)]
            genomic = rc(genomic) if strand == "-" else genomic
            for i, (a, b) in enumerate(zip(bases, genomic)):
                if a == b:
                    continue
                pos = start + (len(bases) - 1 - i if strand == "-" else i) + 1
                if origin == "alt" and (chrom, pos) in snps:
                    want = snps[(chrom, pos)][1]
                    if (rc(want) if strand == "-" else want) == a:
                        continue
                bad.append((block, chrom, start, strand, target, f"base {i}: {a} vs {genomic[i]}"))
                break
    return bad


def check_rows(seqs, recs, output, rows):
    bad = wrong_rows(seqs, recs, rows)
    say(f"{output}: every row is its genome's bases ({len(rows)} rows)", not bad, str(bad[:4]))


def run(work):
    seqs, recs = make(work)
    for pam in ("20bp-NGG-SpCas9.txt", "20bp-TTTV-Cas12a.txt"):
        code, out, err = cli(work, "build-index-only", "--genome", "Genomes/bx", "--pam",
                             f"PAMs/{pam}", "--bDNA", "1", "--bRNA", "1", "--vcf", "list_vcf.txt",
                             "--samplesID", "list_samplesID.txt", "--thread", "2")
        say(f"build-index-only {pam}: ends well, nothing on stderr", code == 0 and not err,
            f"rc {code}: {err[-1500:]}")
    for output, (pam, guide, variant, mode) in SEARCHES.items():
        args = ["complete-search", "--genome", "Genomes/bx", "--guide", guide, "--pam",
                f"PAMs/{pam}", "--mm", "4", "--bDNA", "1", "--bRNA", "1", "--output", output,
                "--thread", "2"]
        if variant:
            args += ["--vcf", "list_vcf.txt", "--samplesID", "list_samplesID.txt", mode]
        code, out, err = cli(work, *args)
        log_error = ""
        if os.path.isfile(os.path.join(work, "Results", output, "log_error.txt")):
            with open(os.path.join(work, "Results", output, "log_error.txt")) as fh:
                log_error = fh.read()
        say(f"{output}: ends well, nothing on stderr", code == 0 and not log_error.strip(),
            f"rc {code}: {(log_error or err)[-1500:]}")
        _head, rows = table(work, output)
        check_rows(seqs, recs, output, rows)
    # NGG under --prune-pam off a PAM-less index of an 8-nt flank, the only index of a working
    # directory of its own: the dedicated index's rows (format 2)
    k8 = os.path.join(work, "k8")
    os.makedirs(k8)
    make(k8)
    code, out, err = cli(k8, "build-index-only", "--genome", "Genomes/bx", "--pam",
                         "PAMs/20bp-NNNNNNNN-NO-PAM.txt", "--bDNA", "1", "--bRNA", "1",
                         "--thread", "2")
    say("build-index-only of the 8-nt flank: ends well, nothing on stderr", code == 0 and not err,
        f"rc {code}: {err[-1500:]}")
    code, out, err = cli(k8, "complete-search", "--genome", "Genomes/bx", "--guide", "g20.txt",
                         "--pam", "PAMs/20bp-NGG-SpCas9.txt", "--mm", "4", "--bDNA", "1", "--bRNA",
                         "1", "--output", "ref_ngg_k8", "--thread", "2", "--prune-pam",
                         "--index-path", "genome_library")
    log_error = ""
    if os.path.isfile(os.path.join(k8, "Results", "ref_ngg_k8", "log_error.txt")):
        with open(os.path.join(k8, "Results", "ref_ngg_k8", "log_error.txt")) as fh:
            log_error = fh.read()
    say("ref_ngg_k8: ends well, nothing on stderr", code == 0 and not log_error.strip(),
        f"rc {code}: {(log_error or err)[-1500:]}")
    _head, rows = table(k8, "ref_ngg_k8")
    check_rows(seqs, recs, "ref_ngg_k8", rows)
    _head, dedicated = table(work, "ref_ngg")
    pruned, own = {key(r) for r in rows}, {key(r) for r in dedicated}
    say(f"ref_ngg_k8: the rows of the dedicated index ({len(pruned)} / {len(own)})", pruned == own,
        f"only pruned {sorted(pruned - own)[:3]}, only dedicated {sorted(own - pruned)[:3]}")


def key(row):
    return (row["Chromosome"], row["Start_coordinate_(highest_CFD)"],
            row["Strand_(highest_CFD)"], row.get("Aligned_protospacer+PAM_ALT_(highest_CFD)"),
            row.get("Aligned_protospacer+PAM_REF_(highest_CFD)"))


def boundary_distance(seq, start, length):
    """How far a target at start..start+length lies from a contig end or an N run."""
    left = seq.rfind("N", 0, start)
    right = seq.find("N", start + length)
    return min(start - (left + 1), (len(seq) if right < 0 else right) - (start + length))


def compare(one, two):
    seqs, recs = genome()
    for output in SEARCHES:
        _h1, rows1 = table(one, output)
        _h2, rows2 = table(two, output)
        # (format 1 reports a minus-strand PAM at a contig's first base on the plus strand,
        # at the max bulges: not its genome's bases, and not a row of format 2)
        wrong = {(c, str(st), sd) for _b, c, st, sd, _t, _w in wrong_rows(seqs, recs, rows1)}
        k1 = {key(r) for r in rows1 if (r["Chromosome"], r["Start_coordinate_(highest_CFD)"],
                                        r["Strand_(highest_CFD)"]) not in wrong}
        k2 = {key(r) for r in rows2}
        say(f"{output}: format 2 holds every format-1 row that is its genome's bases "
            f"({len(k1)} / {len(k2)}; {len(rows1) - len(k1)} not)", k1 <= k2,
            str(sorted(k1 - k2)[:4]))
        more = [r for r in rows2 if key(r) not in k1]
        far = []
        for r in more:
            start = int(r["Start_coordinate_(highest_CFD)"])
            length = len(r["Aligned_protospacer+PAM_REF_(highest_CFD)"].replace("-", ""))
            if boundary_distance(seqs[r["Chromosome"]], start, length) > 4:
                far.append(key(r))
        say(f"{output}: each of its {len(more)} rows more lies at a boundary", not far, str(far[:4]))


def main(argv):
    if argv[:1] == ["--compare"]:
        compare(os.path.abspath(argv[1]), os.path.abspath(argv[2]))
    else:
        work = os.path.abspath(argv[0])
        if os.path.exists(work):
            shutil.rmtree(work)
        os.makedirs(work)
        print(f"engine: {shutil.which('crispritz.py')}", flush=True)
        run(work)
    print(f"{len(failed)} failed" if failed else "all passed", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
