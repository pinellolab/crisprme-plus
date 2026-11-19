#!/usr/bin/env python3
"""crispritz_bins and index_fit against the real engine: the index files it writes, and the
sites it keeps.

    check_real_engine.py <work dir, which it creates> [seed] [q]

PostProcess/test_index_format2.py checks CRISPRme's reader of index format 2 on trees made as
buildTST makes them (a port of saveTST) and its site rule against buildTST's own loops (a port
of searchPAMonGenome and finishLeaf). This checks the ports' premises on the engine installed
on PATH (crispritz.py index-genome and index-info), which CI's light job does not have: run it
inside an image with the tree under test in place and the CRISPRitz under test installed.

For random contigs -- N runs, IUPAC codes, Q (the four-allele symbol, with the argument q, for
an engine that reads it: CRISPRitz 2.10.0), sites planted at contig ends and beside N runs --
and the PAMs of the boundary suite (NGG, NRG, a 6-nt
NNGRRT, the 5' TTTV, and PAM-less K = 3, 4 (5') and 8), each indexed at max bulges 0 to 3:

* every .bin is whole (bin_problem None), and of some, every cut of it is not;
* the number of targets of each contig (numLeaves summed over its parts) is
  crispritz_bins.count_targets by the rule of the index's format (index_format, from its own
  headers), and a contig has a .bin exactly where holds_target says it holds a target;
* index_fit.read_index_header returns, for every .bin, what crispritz.py index-info does
  (format 2 only: CRISPRitz 2.9.1 has no index-info), and for files damaged in one way each
  (cut in its header, a format version of 3, a flag bit it does not know, another header length,
  another guide length in its body, named for another PAM, no leaf, more leaves than it holds,
  no end magic);
* index_fit.write_index_meta records the format and build parameters the headers hold.

Exits non-zero on any failure.
"""

import json
import os
import random
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "PostProcess"))

import crispritz_bins  # noqa: E402
import index_fit  # noqa: E402

PAMS = ("N" * 20 + "NGG 3", "N" * 20 + "NRG 3", "N" * 21 + "NNGRRT 6", "TTTV" + "N" * 20 + " -4",
        "N" * 23 + " 3", "N" * 24 + " -4", "N" * 28 + " 8")
GUIDE = "CTAACAGTTGCTTTTATCAC"
failed = []


def say(cell, ok, why=""):
    print(f"{'PASS' if ok else 'FAIL'}  {cell}{': ' + why if why and not ok else ''}", flush=True)
    if not ok:
        failed.append(cell)


def contigs(r, n, q):
    """``n`` random contigs: plain, with N runs, with IUPAC codes (and Q, with ``q``), each
    with a site planted at its start, its end or beside an N run, now and then."""
    out = {}
    for i in range(n):
        size = r.choice((r.randint(20, 40), r.randint(40, 120), r.randint(120, 600)))
        mix = r.choice(("ACGT", "ACGT" * 8 + "N", "ACGTRYSWKMBDHV" + "ACGT" * 4)
                       + (("ACGT" * 6 + "Q",) if q else ()))
        seq = [r.choice(mix) for _ in range(size)]
        if size > 30 and r.random() < 0.5:
            a = r.randrange(size)
            run = len(seq[a:a + r.randint(1, 25)])
            seq[a:a + run] = ["N"] * run
        site = r.choice((GUIDE + "TGG", "CCA" + GUIDE, "TTTA" + GUIDE + "AAA", GUIDE + "AGAGT"))
        where = r.choice(("start", "end", "n", "none"))
        if where == "start":
            seq[:len(site)] = list(site)
        elif where == "end":
            seq[-len(site):] = list(site)
        elif where == "n" and size > 2 * len(site) + 4:
            a = r.randrange(size - len(site) - 2)
            seq[a:a + len(site) + 2] = list("N" + site + "N")
        out[f"c{i}"] = "".join(seq[:size])
    return out


def damaged_ones(work, index, name, cell):
    """read_index_header against index-info on copies of the format-2 file ``name`` of
    ``index``, each damaged in one way, each alone in a folder (index-info reads folders)."""
    import struct

    with open(os.path.join(index, name), "rb") as fh:
        data = fh.read()
    pam = name.split("_", 1)[0]
    header_len = 28 + len(pam)

    def put(at, value):
        return data[:at] + struct.pack("=i", value) + data[at + 4:]

    cases = {
        "cut in its header": (name, data[:20]),
        "format 3": (name, put(4, 3)),
        "flag bit 8": (name, put(12, struct.unpack_from("=i", data, 12)[0] | 8)),
        "header length": (name, put(8, header_len + 1)),
        "body guide length": (name, put(header_len + 4, struct.unpack_from("=i", data, 16)[0] + 1)),
        "another PAM's name": ("A" + name, data),
        "no leaf": (name, put(header_len, 0)),
        "more leaves than it holds": (name, put(header_len, 1 << 24)),
        "no end magic": (name, data[:-4]),
    }
    wrong = []
    for why, (file_name, damaged) in cases.items():
        folder = os.path.join(work, "damaged", why.replace(" ", "_"))
        shutil.rmtree(folder, ignore_errors=True)
        os.makedirs(folder)
        with open(os.path.join(folder, file_name), "wb") as fh:
            fh.write(damaged)
        info = subprocess.run(["crispritz.py", "index-info", folder], capture_output=True,
                              text=True)
        theirs = json.loads(info.stdout)["files"][file_name] if info.stdout.strip() else None
        ours = index_fit.read_index_header(os.path.join(folder, file_name))
        if ours != theirs or "error" not in ours:
            wrong.append((why, ours, theirs))
    say(cell + " read_index_header refuses a damaged file as index-info does", not wrong,
        str(wrong[:2]))


def main(argv):
    work = os.path.abspath(argv[0])
    seed = int(argv[1]) if len(argv) > 1 else 20261106
    q = argv[2:3] == ["q"]
    engine_format = index_fit.engine_index_format()
    print(f"engine: crispritz.py {index_fit.engine_version()} (index format {engine_format}) at "
          f"{shutil.which('crispritz.py')}", flush=True)
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(work)
    r = random.Random(seed)
    genome = os.path.join(work, "g")
    os.makedirs(genome)
    seqs = contigs(r, 60, q)
    for name, seq in seqs.items():
        with open(os.path.join(genome, f"{name}.fa"), "w") as fh:
            fh.write(f">{name}\n{seq}\n")
    cuts = 0
    for line in PAMS:
        motif, limit = line.split()
        limit = int(limit)
        pam = motif[:-limit] if limit < 0 else motif[-limit:]
        pam_file = os.path.join(work, f"{pam}_{limit}.txt")
        with open(pam_file, "w") as fh:
            fh.write(line + "\n")
        for mb in (0, 1, 2, 3):
            cell = f"{line} -bMax {mb}"
            lib = os.path.join(work, f"lib_{pam}_{limit}_{mb}")
            os.makedirs(lib)
            run = subprocess.run(["crispritz.py", "index-genome", "g", genome + "/", pam_file,
                                  "-bMax", str(mb), "-th", "2"], cwd=lib, capture_output=True,
                                 text=True)
            index = os.path.join(lib, "genome_library", f"{pam}_{mb}_g")
            if run.returncode or run.stderr.strip() or not os.path.isdir(index):
                say(cell + " builds", False, f"rc {run.returncode}: {run.stderr[-500:]}")
                continue
            fmt = crispritz_bins.index_format(index)
            say(cell + f" is in the engine's format ({fmt})", fmt == (engine_format or 1))
            bins = sorted(f for f in os.listdir(index) if f.endswith(".bin"))
            bad = [(b, crispritz_bins.bin_problem(os.path.join(index, b))) for b in bins]
            bad = [x for x in bad if x[1]]
            say(cell + f" every .bin is whole ({len(bins)})", not bad, str(bad[:3]))
            # every cut of the first two .bin is told
            for b in bins[:2]:
                path = os.path.join(index, b)
                with open(path, "rb") as fh:
                    data = fh.read()
                cut = os.path.join(work, "cut", b)
                os.makedirs(os.path.dirname(cut), exist_ok=True)
                missed = []
                for size in range(len(data)):
                    with open(cut, "wb") as fh:
                        fh.write(data[:size])
                    if crispritz_bins.bin_problem(cut) is None:
                        missed.append(size)
                cuts += len(data)
                say(cell + f" every cut of {b} is told ({len(data)})", not missed, str(missed[:5]))
            # the targets of each contig, and which contigs have a .bin
            leaves = {}
            for b in bins:
                contig = re.match(re.escape(pam) + r"_(.*)_[0-9]+\.bin\Z", b).group(1)
                leaves[contig] = leaves.get(contig, 0) + crispritz_bins._leaves(os.path.join(index, b))
            wrong = []
            for name, seq in seqs.items():
                want = crispritz_bins.count_targets(seq, pam, len(motif), limit < 0, mb, fmt)
                if want is not None and want != leaves.get(name, 0):
                    wrong.append((name, leaves.get(name, 0), want))
                if crispritz_bins.holds_target(seq, pam, len(motif), limit < 0, mb, fmt) \
                        != (name in leaves):
                    wrong.append((name, "holds", name in leaves))
            say(cell + " count_targets and holds_target are buildTST's", not wrong, str(wrong[:5]))
            missing, unread = crispritz_bins.missing_bins(index, genome, pam_file, mb)
            say(cell + " missing_bins finds nothing missing", not missing and not unread,
                str((missing, unread)))
            if engine_format == 2:
                info = subprocess.run(["crispritz.py", "index-info", index], capture_output=True,
                                      text=True)
                theirs = json.loads(info.stdout)["files"] if info.stdout.strip() else {}
                ours = {b: index_fit.read_index_header(os.path.join(index, b)) for b in bins}
                say(cell + " read_index_header is index-info's", ours == theirs,
                    str([(b, ours[b], theirs.get(b)) for b in bins if ours[b] != theirs.get(b)][:2]))
            if engine_format == 2 and bins:
                damaged_ones(work, index, bins[0], cell)
            why = index_fit.write_index_meta(index, pam_file, built=True)
            meta = index_fit.read_index_meta(index) or {}
            want = {"index_format": fmt, "build_pam": pam, "side": "5'" if limit < 0 else "3'",
                    "flank_len": len(pam), "guide_len": len(motif) - abs(limit), "max_bulges": mb}
            got = {k: meta.get(k) for k in want}
            say(cell + " index_meta records it", not why and got == want, f"{why} {got} {want}")
    print(f"{cuts} cuts told", flush=True)
    print(f"{len(failed)} failed" if failed else "all passed", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
