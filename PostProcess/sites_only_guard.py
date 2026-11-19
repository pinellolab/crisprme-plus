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

The withheld windows go to the dense-region list as well, which resultIntegrator notes
every row by (High_variant_density_region) and the report bundles:
``<out stem>.high_variant_density_regions.bed``, the columns of the post-analysis's own
(new_simple_analysis), its ``greedy`` column ``withheld`` (no row of the window was
written: the guard withheld its candidates). The SNP side lists each window as the SNP
post-analysis lists a dense one, one line a window (start, end, guide; of the raw targets
of one window the line resultIntegrator reads first, the smallest variants, letters and
bits), so the list is the one the post-analysis writes without the guard, and a row the
guard keeps in such a window (a reference row, a genotyped dataset's) keeps its note; and
the four-allele positions of those windows go to ``<out stem>.four_allele_positions.bed``,
as the post-analysis lists those of every window it reads. The indel side, whose
post-analysis lists no dense window, lists one line a region of ``<out>`` (its densest
window's letters, on the fake-indel contig), its ``greedy`` column ``withheld-region``:
resultIntegrator notes a row by such a region only where no window of the list covers it,
so the note of a row in a window is the one the post-analysis gives it with the guard off
(a region, wider than the windows it holds, would come first: on the rebuilt mega's chr22,
155 rows got a region's note for their own window's).

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
# the dense-region list's columns (new_simple_analysis.py writes them; resultIntegrator and
# the report read them by place), and the greedy column's value for a window the guard withheld
HVDR_SUFFIX = ".high_variant_density_regions.bed"
HVDR_HEADER = ("#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\tiupac_protospacer\tbits"
               "\tfour_allele\tgreedy\n")
WITHHELD = "withheld"
# the indel side's lines are regions of windows on the fake-indel contigs, not windows:
# resultIntegrator notes a row by one only where no window covers it
WITHHELD_REGION = "withheld-region"
FOUR_ALLELE_SUFFIX = ".four_allele_positions.bed"
# the complement of a window's letters, IUPAC codes included (new_simple_analysis's table;
# Q, every base, is its own)
_COMPLEMENT = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv", "TGACYRSWKMDHVBtgacyrswkmdhvb")


def companion_paths(bed_out):
    """(dense-region list, four-allele positions) written beside the guard's BED ``bed_out``."""
    stem = bed_out[:-len(".bed")] if bed_out.endswith(".bed") else bed_out
    return stem + HVDR_SUFFIX, stem + FOUR_ALLELE_SUFFIX


def forward_letters(dna, strand):
    """A raw target's window as the SNP post-analysis lists it: its letters without the
    bulge gaps, on the forward strand (new_simple_analysis's replaceTarget)."""
    letters = dna.replace("-", "")
    return letters.translate(_COMPLEMENT)[::-1] if strand == "-" else letters


def hvdr_line(chrom, start, end, guide, letters, n_seq, flag=WITHHELD):
    """One line of the dense-region list for a window the guard withheld (or, ``flag``
    WITHHELD_REGION, a region of them): its span, guide, variant columns, no carriers (a
    sites-only dataset has none), letters, bits, four-allele positions, and ``flag`` where
    the post-analysis says whether a greedy worst case stands for the window."""
    return "%s\t%d\t%d\t%s\t%d\t.\t%s\t%.2f\t%d\t%s\n" % (
        chrom, start, end, guide, _wd.variant_columns(letters), letters, _wd.bits(n_seq),
        len(_wd.four_allele_offsets(letters)), flag)


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
    # (guide, start) -> [end, max bits, candidates, the densest window's letters and its
    # number of sequences]
    regions = {}
    windows = {}   # SNP side: (start, end, guide) -> (variants, letters, bits) to list it by
    quads = set()  # SNP side: the four-allele positions of the withheld windows
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
            guide = f[1].replace("-", "")
            key = (guide, start)
            have = regions.get(key)
            n_seq = _wd.sequences(dna)
            b = _wd.bits(n_seq)
            letters = forward_letters(dna, f[6] if len(f) > 6 else "+")
            if have is None:
                regions[key] = [end, b, 1, letters, n_seq]
            else:
                if end > have[0]:
                    have[0] = end
                if b > have[1]:
                    have[1], have[3], have[4] = b, letters, n_seq
                have[2] += 1
            if kind == "snp":
                # listed as the SNP post-analysis lists a dense window: one line a window,
                # the one resultIntegrator reads first (_record_dense_window)
                rep = ("%d" % _wd.variant_columns(letters), letters, "%.2f" % b)
                wkey = (start, start + len(letters), guide)
                if wkey not in windows or rep < windows[wkey][0]:
                    windows[wkey] = (rep, n_seq)
                # and its four-allele positions, as it lists those of every window it reads
                # (_note_four_allele: Q is its own complement, so a minus-strand window is
                # only read backwards)
                if "Q" in dna or "q" in dna:
                    window = dna.replace("-", "")
                    if f[6] == "-":
                        window = window[::-1]
                    quads.update(start + o for o in _wd.four_allele_offsets(window))
    lines = []
    by_guide = {}
    for (guide, start), (end, b, n, letters, n_seq) in regions.items():
        by_guide.setdefault(guide, []).append((start, end, b, n, letters, n_seq))
    for guide in sorted(by_guide):
        cur = None
        for start, end, b, n, letters, n_seq in sorted(by_guide[guide]):
            if cur is not None and start <= cur[1]:
                cur[1] = max(cur[1], end)
                if b > cur[2]:
                    cur[2], cur[4], cur[5] = b, letters, n_seq
                cur[3] += n
            else:
                if cur is not None:
                    lines.append((cur[0], guide, cur))
                cur = [start, end, b, n, letters, n_seq]
        if cur is not None:
            lines.append((cur[0], guide, cur))
    lines.sort(key=lambda t: (t[0], t[1]))
    hvdr_out, quads_out = companion_paths(bed_out)
    with open(hvdr_out + ".tmp", "w") as fh:
        fh.write(HVDR_HEADER)
        if kind == "snp":
            for wkey in sorted(windows):
                (_nvar, letters, _bits), n_seq = windows[wkey]
                fh.write(hvdr_line(chrom, wkey[0], wkey[1], wkey[2], letters, n_seq))
        else:
            for _s, guide, (start, end, _b, _n, letters, n_seq) in lines:
                fh.write(hvdr_line(chrom, start, end, guide, letters, n_seq, WITHHELD_REGION))
    if quads:
        with open(quads_out + ".tmp", "w") as fh:
            fh.write("#chrom\tstart\tend\n")
            for p in sorted(quads):
                fh.write("%s\t%d\t%d\n" % (chrom, p, p + 1))
    with open(bed_out + ".tmp", "w") as bed:
        bed.write(f"## CRISPRme sites-only density guard: candidates of the sites-only dataset "
                  f"{dataset} whose window spells 2^{bits} sequences or more, withheld from "
                  f"its {kind.upper() if kind == 'snp' else kind} post-analysis "
                  f"(window_density.guard_withholds)\n")
        bed.write(f"## bits={bits}\tdataset={dataset}\tchrom={chrom}\tkind={kind}\t"
                  f"candidates={total}\twithheld={withheld}\n")
        bed.write(BED_HEADER)
        for _s, guide, (start, end, b, n, _letters, _n_seq) in lines:
            bed.write(f"{chrom}\t{start}\t{end}\t{guide}\t{n}\t{b:.2f}\n")
    os.replace(bed_out + ".tmp", bed_out)
    os.replace(hvdr_out + ".tmp", hvdr_out)
    if quads:
        os.replace(quads_out + ".tmp", quads_out)
    os.replace(tmp, targets)
    what = "SNP" if kind == "snp" else "indel"
    listed = (f"{len(windows):,} window(s)" if kind == "snp" else f"the {len(lines):,} region(s)")
    print(f"Sites-only density guard ({dataset}, {chrom}, {what}): {withheld:,} of "
          f"{total:,} candidate(s) in windows of 2^{bits} sequences or more withheld, in "
          f"{len(lines):,} region(s); listed in {os.path.basename(bed_out)}, and {listed} "
          f"in the dense-region list {os.path.basename(hvdr_out)}", flush=True)
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
        for p in (targets + ".guard.tmp", bed_out + ".tmp",
                  *(c + ".tmp" for c in companion_paths(bed_out))):
            try:
                os.remove(p)
            except OSError:
                pass
        print(f"WARNING [sites-only guard]: {dataset} {chrom} {kind}: not applied ({err}); "
              f"its candidates are analysed unguarded", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
