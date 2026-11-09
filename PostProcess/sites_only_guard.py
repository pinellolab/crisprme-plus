#!/usr/bin/env python
"""The sites-only density guard (CRISPRme 2.10.0): a sites-only dataset's candidates whose
window spells 2^bits sequences or more are withheld from its SNP and its indel
post-analysis, and the regions they lie in are listed.

A sites-only dataset (the mega: aggregate frequencies, no genotypes) has nothing to tell
which combinations of a window's alleles anyone carries: its post-analysis reports every
window by its worst case. With every allele coded (CRISPRitz 2.10.0) the immunoglobulin
loci and some tandem repeats make windows of 2^40 sequences and more, which match nearly
any guide: on chr14 the rebuilt mega gives 13x (SNP) and 23x (indel) the raw targets of
the shipped one, 92-93% of them in windows of 2^20 sequences or more, a search of five
guides at 4 mismatches + 1 + 1 bulges projected at 12-16 h against 46 min, and on chr22
the new rows at CFD 0.2 or more are combinations of 10 to 17 variants of AF product down
to 1e-33 (the judges of D2, DECISION_PANEL_2.10.0). The guard runs on the raw targets of
the search, before the post-analysis reads them, so the time they would cost is not
spent; the rows of the other windows, and the reference rows, are as before. It reads the
window as the post-analysis does: the letters of the raw target's aligned target (a SNP
target on the variant genome; an indel one on its fake-indel contig, whose flanks carry
the SNP codes), window_density.guard_withholds.

It applies to a sites-only dataset only -- its registry's manifest says so
(``aggregation`` ``info_af``, or ``data_type`` ``sites-only``), whatever the search's mode
(such a dataset's post-analysis is the population-level one in every search). A genotyped
dataset's rows are its carriers' own haplotypes (per sample) or are checked against its
phased genotypes by the carrier gate, which demotes and drops none (population-level):
the guard leaves them as they are.

``CRISPRME_SITES_ONLY_GUARD``: the bits (default window_density.SITES_ONLY_GUARD_BITS,
20); a negative value turns the guard off (complete-search --sites-only-density-guard N,
--no-sites-only-density-guard).

Output: the targets file rewritten without the withheld lines (atomically), and a BED of
the regions (``<out>``): one line per region of overlapping withheld windows of one guide
(0-based start, end, the guide, how many candidates, the most bits), its first lines
saying the bound, the dataset, the chromosome and the counts; for the indel side the
regions are on the reference, each window placed by its fake-indel contig's start (the
log's window), so an insertion's region can run a few bases past it. One line on stdout.

usage: sites_only_guard.py snp|indel <targets file> <Dictionaries folder> <dataset> <chrom>
                           <bed out>

Never fails the run: anything unexpected (no indel log to place the indel side's windows
by, among others) leaves the targets file as it was, unguarded, with a WARNING on stdout
(stderr is fatal to a post-analysis stage); a candidate it cannot place on the reference is
left to the post-analysis, so every one withheld is listed.
"""

import bisect
import gzip
import json
import os
import sys

import window_density as _wd

ENV = "CRISPRME_SITES_ONLY_GUARD"
BED_SUFFIX = ".sites_only_density_guard.bed"
BED_HEADER = "#chrom\tstart\tend\tguide\tcandidates_withheld\tmax_bits\n"


def guard_bits(environ=None):
    """The guard's bound in bits, from CRISPRME_SITES_ONLY_GUARD (default
    SITES_ONLY_GUARD_BITS); a negative one is off."""
    raw = (environ if environ is not None else os.environ).get(ENV, "")
    if raw.strip() == "":
        return _wd.SITES_ONLY_GUARD_BITS
    try:
        return int(raw)
    except ValueError:
        print(f"WARNING [sites-only guard]: {ENV}={raw!r} is no whole number: the guard "
              f"takes {_wd.SITES_ONLY_GUARD_BITS} bits", flush=True)
        return _wd.SITES_ONLY_GUARD_BITS


def registry_kind(dictionaries, dataset):
    """'sites-only' when the dataset's registry manifest says it has no genotypes
    (aggregation info_af, or data_type sites-only), else its data_type or None."""
    rdir = os.path.join(dictionaries, f"registry_{dataset}")
    try:
        names = sorted(os.listdir(rdir))
    except OSError:
        return None
    for name in names:
        if name.startswith("reg_") and name.endswith(".idx"):
            try:
                with open(os.path.join(rdir, name)) as fh:
                    man = json.load(fh)
            except (OSError, ValueError):
                continue
            if man.get("aggregation") == "info_af" or man.get("data_type") == "sites-only":
                return "sites-only"
            return man.get("data_type")
    return None


class FakeToReal(object):
    """Places a position of a chromosome's fake-indel genome on the reference: the start
    of the window of the contig it lies in (the indel log's CHR, <chrom>_<start>-<end>_<n>,
    and FAKEPOS, the contig's range) plus its offset into the contig."""

    def __init__(self, log_dir, chrom):
        self.starts, self.rows = [], []
        path = None
        for name in (f"log{chrom}.txt.gz", f"log{chrom}.txt"):
            if os.path.isfile(os.path.join(log_dir, name)):
                path = os.path.join(log_dir, name)
                break
        if path is None:
            return
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(header)}
            ci, fi = ix.get("CHR", 0), ix.get("FAKEPOS", 5)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    a = int(f[fi].split(",")[0])
                    win = int(f[ci].rsplit("_", 2)[-2].split("-")[0])
                except (IndexError, ValueError):
                    continue
                self.rows.append((a, win))
        self.rows.sort()
        self.starts = [a for a, _w in self.rows]

    def __call__(self, fake_pos):
        i = bisect.bisect_right(self.starts, fake_pos) - 1
        if i < 0:
            return None
        a, win = self.rows[i]
        return win + (fake_pos - a)


def run(kind, targets, dictionaries, dataset, chrom, bed_out, bits=None):
    """Withhold, in the raw targets file ``targets``, the candidates of a sites-only
    dataset whose window spells 2^bits sequences or more; write the regions to
    ``bed_out``. Returns (candidates, withheld), or None when the guard does not apply
    (another kind of dataset, or the guard off)."""
    if bits is None:
        bits = guard_bits()
    if bits < 0 or registry_kind(dictionaries, dataset) != "sites-only":
        return None
    place = None
    if kind == "indel":
        place = FakeToReal(os.path.join(dictionaries, f"log_indels_{dataset}"), chrom)
        if not place.rows:
            raise ValueError(f"no indel log of {chrom} in log_indels_{dataset}/ to place "
                             f"its windows on the reference by")
    regions = {}   # (guide, start) -> [end, max bits, candidates]
    total = withheld = 0
    tmp = targets + ".guard.tmp"
    with open(targets) as fh, open(tmp, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.split("\t", 7)
            if len(f) < 7:
                out.write(line)
                continue
            total += 1
            dna = f[2]
            if not _wd.guard_withholds(dna, bits):
                out.write(line)
                continue
            try:
                start = int(f[4])
            except ValueError:
                start = None
            if start is not None and place is not None:
                start = place(start)
            if start is None:
                # every candidate withheld is listed: one it cannot place on the
                # reference is left to the post-analysis
                out.write(line)
                continue
            withheld += 1
            end = start + len(dna) - dna.count("-")
            key = (f[1].replace("-", ""), start)
            have = regions.get(key)
            b = _wd.bits(_wd.sequences(dna))
            if have is None:
                regions[key] = [end, b, 1]
            else:
                if end > have[0]:
                    have[0] = end
                if b > have[1]:
                    have[1] = b
                have[2] += 1
    lines = []
    by_guide = {}
    for (guide, start), (end, b, n) in regions.items():
        by_guide.setdefault(guide, []).append((start, end, b, n))
    for guide in sorted(by_guide):
        cur = None
        for start, end, b, n in sorted(by_guide[guide]):
            if cur is not None and start <= cur[1]:
                cur[1] = max(cur[1], end)
                cur[2] = max(cur[2], b)
                cur[3] += n
            else:
                if cur is not None:
                    lines.append((cur[0], guide, cur))
                cur = [start, end, b, n]
        if cur is not None:
            lines.append((cur[0], guide, cur))
    lines.sort(key=lambda t: (t[0], t[1]))
    with open(bed_out + ".tmp", "w") as bed:
        bed.write(f"## CRISPRme sites-only density guard: candidates of the sites-only dataset "
                  f"{dataset} whose window spells 2^{bits} sequences or more, withheld from "
                  f"its {kind.upper() if kind == 'snp' else kind} post-analysis "
                  f"(window_density.guard_withholds)\n")
        bed.write(f"## bits={bits}\tdataset={dataset}\tchrom={chrom}\tkind={kind}\t"
                  f"candidates={total}\twithheld={withheld}\n")
        bed.write(BED_HEADER)
        for _s, guide, (start, end, b, n) in lines:
            bed.write(f"{chrom}\t{start}\t{end}\t{guide}\t{n}\t{b:.2f}\n")
    os.replace(bed_out + ".tmp", bed_out)
    os.replace(tmp, targets)
    what = "SNP" if kind == "snp" else "indel"
    print(f"Sites-only density guard ({dataset}, {chrom}, {what}): {withheld:,} of "
          f"{total:,} candidate(s) in windows of 2^{bits} sequences or more withheld, in "
          f"{len(lines):,} region(s); listed in {os.path.basename(bed_out)}", flush=True)
    return total, withheld


def read_bed_counts(path):
    """{bits, dataset, chrom, kind, candidates, withheld} of a guard BED's second line, or
    None (generate_report, check_indel_rows)."""
    try:
        with open(path) as fh:
            fh.readline()
            second = fh.readline()
    except OSError:
        return None
    if not second.startswith("## bits="):
        return None
    out = {}
    for kv in second[3:].rstrip("\n").split("\t"):
        k, _s, v = kv.partition("=")
        out[k] = int(v) if k in ("bits", "candidates", "withheld") else v
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 6 or argv[0] not in ("snp", "indel"):
        print("usage: sites_only_guard.py snp|indel <targets> <Dictionaries> <dataset> "
              "<chrom> <bed out>", flush=True)
        return 0
    kind, targets, dictionaries, dataset, chrom, bed_out = argv
    try:
        run(kind, targets, dictionaries, dataset, chrom, bed_out)
    except Exception as err:  # noqa: BLE001 - the guard never fails the run
        for p in (targets + ".guard.tmp", bed_out + ".tmp"):
            try:
                os.remove(p)
            except OSError:
                pass
        print(f"WARNING [sites-only guard]: {dataset} {chrom} {kind}: not applied ({err}); "
              f"its candidates are analysed unguarded", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
