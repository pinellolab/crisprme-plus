#!/usr/bin/env python
"""Where the indel post-analysis puts an indel row (analisi_indels_NNN.py): its
Start_coordinate, its Cluster_Position and its REF text, for every bulge type, on both
strands, with a 3' PAM (SpCas9, NRG) and a 5' PAM (Cas12a, TTTV).

A small reference chromosome holds, on its carriers' haplotypes only, a perfect target
with one or two RNA bulges (the guide less one or two bases), one RNA and one DNA bulge,
one DNA bulge, or none; the reference has three more bases inside it, which the carriers delete. The
fake-indel windows are cut as CRISPRitz add-variants cuts them (25 bp, the ALT, 27 bp),
and each planted alignment reaches the SHIPPED analyzer as the engine reports it:
Position is its first aligned base on the fake contig, the leftmost, on either strand;
Cluster_Position is the Position less the RNA bulges and plus the DNA bulges on the
strand CRISPRitz searches forward (+ with a 3' PAM, - with a 5' PAM), else the Position
(CRISPRitz 2.9.1 searchOnTST.cpp; a real search of these sites through crisprme.py gives
the same). Each row must start at the planted base, keep the engine's offset to its
cluster position, and show as REF the reference bases its target holds from that start,
gapped as its target, with the SNP+indel integration on and with CRISPRME_INDEL_SNP=0
(as on a fake-indel genome built without the SNP overlay) alike.

Before 2.8.1 an RNA-bulge row took the cluster position on the + strand with a 3' PAM
and on the - strand with a 5' PAM, so it started one base before its first aligned base,
and the REF text of an RNA,DNA row, and of a 5' PAM RNA-bulge row, was a base off.

Since 2.10.0 a target whose first base is one an insertion adds starts at the insertion's
anchor (INSIDE), and a PAM holding inserted bases, or a deletion's junction, is judged
wherever in its repeat the VCF writes the indel (indel_snp_cis.pam_reference_place): a
target whose PAM an insertion adds right after a reference protospacer does not move, and
its PAM is not created on a tandem duplication of the PAM's bases or where the reference
holds another PAM beside the protospacer (DUPS); the geometries where that rule reads
another PAM than 2.10.0's first fixes did (PLACED: a repeat insertion moved off the PAM, the
reference PAM beside the protospacer of a row that starts on its anchor, a deletion across
a - strand PAM, Cas12a on both strands).

Needs pandas and numpy (the analyzer imports them); intervaltree is stubbed when it is
not installed. Run: python -m unittest test_indel_row_coordinates -v
"""

import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

from test_indel_overlaid_rows import _STUB_INTERVALTREE, _SKIP_WC, _gnu_wc

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYZER = os.path.join(HERE, "analisi_indels_NNN.py")
CHROM = "chrT"
VCF = "synth"
G20 = "ACTGAAATCTGTAAGCAGGC"     # SpCas9 guide, NRG
G23 = "GATCCAGTCACGTTAGCATGCAA"  # Cas12a guide, TTTV
COMP = str.maketrans("ACGTacgt-", "TGCAtgca-")
JUNK = ("CAC", "GTG", "TCT", "AGA")

# (name, crRNA, DNA) of each planted alignment, PAM included, gaps '-': an RNA bulge is a
# gap in the DNA, a DNA bulge a gap in the crRNA. The guide's bulged bases sit in no run,
# so each alignment has one place for its gaps.
SITES = {
    "3p": [
        ("RNA", G20 + "NNN", G20[:8] + "-" + G20[9:] + "TGG"),
        ("RNA", G20 + "NNN", G20[:8] + "-" + G20[9:15] + "-" + G20[16:] + "AGG"),
        ("RNA,DNA", G20[:3] + "-" + G20[3:] + "NNN",
         G20[:3] + "C" + G20[3:15] + "-" + G20[16:] + "AGG"),
        ("DNA", G20[:11] + "-" + G20[11:] + "NNN", G20[:11] + "C" + G20[11:] + "TGG"),
        ("X", G20 + "NNN", G20 + "CGG"),
    ],
    "5p": [
        ("RNA", "NNNN" + G23, "TTTA" + G23[:9] + "-" + G23[10:]),
        ("RNA", "NNNN" + G23, "TTTC" + G23[:9] + "-" + G23[10:16] + "-" + G23[17:]),
        ("RNA,DNA", "NNNN" + G23[:5] + "-" + G23[5:],
         "TTTC" + G23[:5] + "T" + G23[5:15] + "-" + G23[16:]),
        ("DNA", "NNNN" + G23[:11] + "-" + G23[11:], "TTTG" + G23[:11] + "A" + G23[11:]),
        ("X", "NNNN" + G23, "TTTA" + G23),
    ],
}
PAMS = {"3p": ("NNNNNNNNNNNNNNNNNNNNNRG 3", G20 + "NNN"),
        "5p": ("TTTVNNNNNNNNNNNNNNNNNNNNNNN -4", "NNNN" + G23)}
MOTIF = {"3p": (slice(-3, None), "NRG"), "5p": (slice(0, 4), "TTTV")}
BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "V": "ACG", "N": "ACGT"}


def rc(s):
    return s.translate(COMP)[::-1]


class _Planted(object):
    """One PAM side's reference, indels, fake-indel windows and engine rows."""

    def __init__(self, side):
        rng = random.Random(4111)
        self.side = side
        seq = [rng.choice("ACGT") for _ in range(3000 + 600 * 2 * len(SITES[side]))]
        self.sites = []  # dict(bulge, strand, start0, crrna, dna, pos1, ref, alt, sample)
        start0 = 1000
        for bulge, crrna, dna in SITES[side]:
            for strand in "+-":
                text = dna.replace("-", "")
                fwd = text if strand == "+" else rc(text)
                q = len(fwd) // 2  # the deleted bases sit after fwd[:q]
                junk = next(j for j in JUNK if j[-1] != fwd[q - 1] and j[0] != fwd[q])
                for i, c in enumerate(fwd[:q] + junk + fwd[q:]):
                    seq[start0 + i] = c
                self.sites.append(dict(
                    bulge=bulge, strand=strand, start0=start0, crrna=crrna, dna=dna,
                    pos1=start0 + q, ref=fwd[q - 1] + junk, alt=fwd[q - 1],
                    sample="S%02d" % (len(self.sites) + 1)))
                start0 += 600
        self.ref = "".join(seq)

    def windows(self):
        """[(site, real window start, window)] as CRISPRitz add-variants cuts them."""
        out = []
        for s in self.sites:
            w0 = s["pos1"] - 26
            win = (self.ref[w0:w0 + 25] + s["alt"]
                   + self.ref[w0 + 25 + len(s["ref"]):w0 + 25 + len(s["ref"]) + 27])
            out.append((s, w0, win))
        return out

    def engine_rows(self, windows, starts):
        """Each planted alignment as CRISPRitz reports it on the fake contig."""
        rows = []
        for (s, w0, _win), f0 in zip(windows, starts):
            pos = f0 + s["start0"] - w0
            n_rna, n_dna = s["dna"].count("-"), s["crrna"].count("-")
            forward = (s["strand"] == "+") == (self.side == "3p")
            cluster = pos - n_rna + n_dna if forward else pos
            nb = n_rna + n_dna
            rows.append([s["bulge"], s["crrna"], s["dna"], "fake" + CHROM, str(pos),
                         str(cluster), s["strand"], "0", str(nb), str(nb)])
        return rows


# Targets whose first base (the leftmost, forward) is one an insertion adds (2.10.0):
# (strand, the target as shown, how many of its leftmost forward bases the insertion
# adds, the added bases before them, the reference base before the anchor). The row
# starts at the insertion's anchor; its PAM is reference bases, wholly inserted bases,
# or partly inserted ones, whose reference PAM is the one beside the protospacer: there,
# the base before the anchor decides whether the reference holds the motif.
INSIDE = {
    "3p": [
        ("+", G20 + "TGG", 6, "TTCAT", None),          # the PAM on the reference: TGG
        ("-", G20 + "TGG", 9, "GA", None),             # the PAM wholly inserted
        ("-", G20 + "AGG", 2, "GT", "C"),              # partly: CTT beside it, AAG on -
        ("-", G20 + "AGG", 2, "GT", "A"),              # partly: ATT, AAT on -: created
    ],
    "5p": [
        ("+", "TTTA" + G23, 6, "GA", None),            # the PAM wholly inserted
        ("+", "TTTC" + G23, 2, "G", "T"),              # partly: TTTC beside the spacer
        ("+", "TTTC" + G23, 2, "G", "G"),              # partly: GTTC: created
        ("-", "TTTG" + G23, 3, "CA", None),            # the PAM on the reference: TTTG
    ],
}


class _PlantedInsertions(object):
    """One PAM side's targets that start on a base an insertion adds (INSIDE): the
    reference, the insertions, their fake-indel windows and engine rows."""

    def __init__(self, side):
        rng = random.Random(2820)
        self.side = side
        seq = [rng.choice("ACGT") for _ in range(3000 + 600 * len(INSIDE[side]))]
        self.sites = []
        start0 = 1000
        for strand, shown, m, head, before in INSIDE[side]:
            fwd = shown if strand == "+" else rc(shown)
            anchor = start0 - 1  # 0-based
            if before is not None:
                seq[anchor - 1] = before
            seq[anchor] = "T"
            for i, c in enumerate(fwd[m:]):
                seq[start0 + i] = c
            self.sites.append(dict(
                strand=strand, shown=shown, fwd=fwd, m=m, anchor=anchor, pos1=start0,
                ref="T", alt="T" + head + fwd[:m], first_off=25 + 1 + len(head),
                sample="S%02d" % (len(self.sites) + 1)))
            start0 += 600
        self.ref = "".join(seq)

    windows = _Planted.windows

    def engine_rows(self, windows, starts):
        crrna = PAMS[self.side][1]
        return [["X", crrna, s["shown"], "fake" + CHROM, str(f0 + s["first_off"]),
                 str(f0 + s["first_off"]), s["strand"], "0", "0", "0"]
                for (s, _w0, _win), f0 in zip(windows, starts)]


# Targets whose PAM is wholly bases an insertion adds, on a row that starts on a reference
# base (2.10.0): the protospacer is reference bases up to the insertion's anchor, its last
# base, and the insertion adds the PAM's bases right after it. (strand, the target as
# shown, the reference's bases after the protospacer (forward), the expected PAM_gen.) The
# row does not move. Where the insertion copies the reference's own PAM there (a tandem
# duplication) the REF alignment is the ALT, and a placement past it leaves the PAM on the
# reference's own: not created; another valid PAM beside the protospacer: not created
# either (the PAM's place beside the protospacer, whose last base is a reference one); no
# PAM there, it is.
DUPS = {
    "3p": [
        ("+", G20 + "TGG", "TGG", "n"),        # a tandem duplication: REF = ALT
        ("+", G20 + "TGG", "AGG", "n"),        # the reference's AGG at the PAM's place
        ("+", G20 + "TGG", "TCA", "TGG"),      # no PAM there: created
    ],
    "5p": [
        ("-", "TTTA" + G23, "TAAA", "n"),      # a tandem duplication of TAAA (TTTA on -)
        ("-", "TTTA" + G23, "GAAA", "n"),      # TTTC on - at the PAM's place
        ("-", "TTTA" + G23, "TAAT", "TTTA"),   # ATTA on -: created
    ],
}


class _PlantedDuplications(object):
    """One PAM side's targets whose PAM an insertion adds after a reference protospacer
    (DUPS): the reference, the insertions, their fake-indel windows and engine rows."""

    def __init__(self, side):
        rng = random.Random(2822)
        self.side = side
        kp = len(MOTIF[side][1])
        seq = [rng.choice("ACGT") for _ in range(3000 + 600 * len(DUPS[side]))]
        self.sites = []
        start0 = 1000
        for strand, shown, after, pam_gen in DUPS[side]:
            fwd = shown if strand == "+" else rc(shown)
            prot, pam_fwd = fwd[:-kp], fwd[-kp:]
            for i, c in enumerate(prot + after):
                seq[start0 + i] = c
            anchor = start0 + len(prot) - 1  # 0-based: the protospacer's last base
            self.sites.append(dict(
                strand=strand, shown=shown, fwd=fwd, start0=start0, pos1=anchor + 1,
                ref=prot[-1], alt=prot[-1] + pam_fwd, pam_gen=pam_gen,
                first_off=25 - (len(prot) - 1), sample="S%02d" % (len(self.sites) + 1)))
            start0 += 600
        self.ref = "".join(seq)

    windows = _Planted.windows
    engine_rows = _PlantedInsertions.engine_rows


# PAMs judged wherever in its repeat the VCF writes the indel (2.10.0): (name, strand, the
# reference bases put at the site, the record (offset of its VCF POS's base, 0-based, from
# the site's first base; REF; ALT), the target as shown, the offset of its first base on the
# haplotype carrying the record, the offset the row starts at, the expected PAM_gen). Each
# one's PAM_gen differs from 2.10.0's first fixes' (in brackets), which judged a PAM holding
# inserted bases, or a deletion's junction, where the VCF writes the indel.
def _placed_sites():
    rc20 = rc(G20)
    f23 = rc(G23)
    return {
        "3p": [
            # (-) C>CCTGC repeats the 4 bases after its anchor; the alignment that starts on
            # the insertion's first base, its PAM CAG inserted bases, is moved off: n [CAG]
            ("A2", "-", "GA" + "C" + rc20[2:] + "GGAA", (2, "C", "CCTGC"),
             None, 3, 2, "n"),
            # (+) a tandem duplication of a 27-base unit holding a target: the copy's target
            # starts inside the insertion; moved one unit on, its PAM is the reference's: n
            # [TGG]
            ("A4", "+", "C" + "TT" + G20 + "TGG" + "AA" + "CGTACGTA",
             (0, "C", "C" + "TT" + G20 + "TGG" + "AA"), G20 + "TGG", 3, 0, "n"),
            # (-) C>CCA after a reference CTC (GAG on -), the row on the anchor: n [GCG
            # window: CGC]
            ("A1", "-", "CT" + "C" + rc20 + "TTAA", (2, "C", "CCA"), G20 + "TGG", 2, 2, "n"),
            # (+) G[16:] + TGG inserted after G[15], no repeat: created [the window AGG: n]
            ("A3", "+", G20[:16] + "TTTT" + "AGG" + "CATTAC",
             (15, G20[15], G20[15] + G20[16:] + "TGG"), G20 + "TGG", 0, 0, "TGG"),
            # (-) CA>C inside the PAM: the reference beside the protospacer CAG (CTG on -):
            # created [the window CCA, TGG on -: n]
            ("D1", "-", "CCAG" + rc20 + "TTGCA", (1, "CA", "C"), G20 + "CGG", 0, 0, "CGG"),
            # (-) TTCCA after the anchor G, a target from its 3rd base, the PAM wholly
            # inserted, the reference CCG (CGG on -) beside the protospacer: n [TGG]
            ("W2", "-", "CCG" + rc20 + "TTGCA", (2, "G", "GTTCCA"), G20 + "TGG", 5, 2, "n"),
        ],
        "5p": [
            # (+) T>TTT at the end of a T run (unnormalized): placed 3 bases left, the PAM
            # TTTA is reference bases: n [the window TAGA: TTTA]
            ("B1", "+", "CC" + "TT" + "T" + "A" + G23 + "GGCC", (4, "T", "TTT"),
             "TTTA" + G23, 4, 4, "n"),
            # (-) the spacer's last 4 forward bases and the PAM inserted, no repeat: created
            # [the window TAAA, TTTA on -: n]
            ("B3", "-", f23[:19] + "CCCC" + "TAAA" + "GTCAGT",
             (18, f23[18], f23[18] + f23[19:] + "TAAA"), "TTTA" + G23, 0, 0, "TTTA"),
        ],
    }


class _PlantedPlacements(object):
    """One PAM side's _placed_sites: the reference, the records, their fake-indel windows
    and engine rows."""

    def __init__(self, side):
        rng = random.Random(2100)
        self.side = side
        sites = _placed_sites()[side]
        seq = [rng.choice("ACGT") for _ in range(3000 + 600 * len(sites))]
        self.sites = []
        start0 = 1000
        for name, strand, put, (off, ref, alt), shown, first, at, pam_gen in sites:
            for i, c in enumerate(put):
                seq[start0 + i] = c
            hap = "".join(seq[:start0 + off]) + alt + "".join(seq[start0 + off + len(ref):])
            if shown is None:  # the alignment that starts on the first inserted base
                shown = rc(hap[start0 + first:start0 + first + 23])
            text = shown if strand == "+" else rc(shown)
            assert hap[start0 + first:start0 + first + len(text)] == text, name
            pos1 = start0 + off + 1
            self.sites.append(dict(
                name=name, strand=strand, shown=shown, pos1=pos1, ref=ref, alt=alt,
                start=start0 + at, pam_gen=pam_gen, first_off=25 + first - off,
                sample="S%02d" % (len(self.sites) + 1)))
            start0 += 600
        self.ref = "".join(seq)

    windows = _Planted.windows
    engine_rows = _PlantedInsertions.engine_rows


class RowCoordinates(unittest.TestCase):
    def analyze(self, p, side, gate="1"):
        """The SHIPPED analyzer on planted sites ``p`` (their reference, fake-indel
        windows and engine rows): {output suffix: {indel token: row columns}}."""
        tmp = tempfile.mkdtemp(prefix="indel_row_coordinates_")
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        windows = p.windows()
        fake, starts = [], []
        for _s, _w0, win in windows:
            if fake:
                fake.append("N")
            starts.append(sum(len(x) for x in fake))
            fake.append(win)
        rows = p.engine_rows(windows, starts)
        with open(os.path.join(tmp, CHROM + ".fa"), "w") as fh:
            fh.write(">%s\n%s\n" % (CHROM, p.ref))
        d = os.path.join(tmp, "Dictionaries")
        logd = os.path.join(d, "log_indels_" + VCF)
        os.makedirs(logd)
        with open(os.path.join(logd, "log%s.txt" % CHROM), "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            for k, (s, w0, win) in enumerate(windows):
                fh.write("%s_%d-%d_%d\t%s\trs%d\t0.05\t%s_%d_%s_%s\t%d,%d\t%s\n" % (
                    CHROM, w0, w0 + len(win), k + 1, s["sample"], k + 1, CHROM, s["pos1"],
                    s["ref"], s["alt"], starts[k], starts[k] + len(win),
                    p.ref[w0:w0 + len(win)]))
        igd = os.path.join(d, "indel_genotypes_" + VCF)
        os.makedirs(igd)
        with gzip.open(os.path.join(igd, "gt_indel_%s.tsv.gz" % CHROM), "wt") as fh:
            for s in sorted(p.sites, key=lambda x: x["pos1"]):
                fh.write("%d_%s_%s\t%s:1|0\n" % (s["pos1"], s["ref"], s["alt"], s["sample"]))
        pam_line, guide = PAMS[side]
        with open(os.path.join(tmp, "pam.txt"), "w") as fh:
            fh.write(pam_line + "\n")
        with open(os.path.join(tmp, "guides.txt"), "w") as fh:
            fh.write(guide + "\n")
        cluster = os.path.join(tmp, "targets.cluster.txt")
        with open(cluster, "w") as fh:
            fh.write("#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
                     "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\t"
                     "Annotation Type\tReal Guide\tID\tAF\tSNP\n")
            for r in rows:
                fh.write("\t".join(r + ["n", "n", "n", "n", guide, "n", "n", "n"]) + "\n")
        out = os.path.join(tmp, "out", "fake" + CHROM)
        os.makedirs(os.path.dirname(out))
        argv = [sys.executable, ANALYZER, "", cluster, out, logd,
                os.path.join(tmp, "pam.txt"), "4", os.path.join(tmp, CHROM + ".fa"),
                os.path.join(tmp, "guides.txt"), "1", "2"]
        env = {k: v for k, v in os.environ.items()
               if k not in ("MAMBA_EXE", "CONDA_EXE", "CRISPRME_CONDA_EXE")}
        env["PATH"] = os.path.dirname(sys.executable) + ":/usr/bin:/bin"
        env["CRISPRME_INDEL_SNP"] = gate
        try:
            import intervaltree  # noqa: F401
        except ImportError:
            stub = os.path.join(tmp, "stub")
            os.makedirs(stub)
            with open(os.path.join(stub, "intervaltree.py"), "w") as fh:
                fh.write(_STUB_INTERVALTREE)
            env["PYTHONPATH"] = stub + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(argv, cwd=HERE, env=env, capture_output=True, text=True,
                              timeout=600)
        self.assertEqual(proc.returncode, 0, proc.stdout[-2000:] + proc.stderr[-2000:])
        self.assertEqual(proc.stderr, "")
        got = {}
        for suffix in (".bestCFD_INDEL.txt", ".bestCRISPR_BULGE_INDEL.txt",
                       ".bestmmblg_INDEL.txt"):
            with open(out + suffix) as fh:
                next(fh)
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    got.setdefault(suffix, {})[c[17]] = c
        self.assertEqual(sorted(len(v) for v in got.values()), [len(p.sites)] * 3)
        return got

    def run_side(self, side, gate="1"):
        p = _Planted(side)
        got = self.analyze(p, side, gate)
        for s in p.sites:
            token = "%s_%d_%s_%s" % (CHROM, s["pos1"], s["ref"], s["alt"])
            n_rna, n_dna = s["dna"].count("-"), s["crrna"].count("-")
            forward = (s["strand"] == "+") == (side == "3p")
            ref = p.ref[s["start0"]:s["start0"] + len(s["dna"].replace("-", ""))]
            ref = list(ref if s["strand"] == "+" else rc(ref))
            for i, ch in enumerate(s["dna"]):
                if ch == "-":
                    ref.insert(i, "-")
            for suffix, by_token in got.items():
                c = by_token[token]
                where = "%s %s %s(%s) %s" % (side, suffix, s["bulge"], s["strand"], c[2])
                self.assertEqual(c[2], s["dna"], where)
                self.assertEqual(int(c[4]), s["start0"], where + ": Start_coordinate")
                self.assertEqual(int(c[4]) - int(c[5]), (n_rna - n_dna) if forward else 0,
                                 where + ": Cluster_Position")
                self.assertEqual(c[18].upper(), "".join(ref), where + ": REF text")
                # PAM_gen: the ALT's PAM when the reference bases it sits on, read through
                # the deletion (the junk lies inside the protospacer), lack the motif. Up to
                # 2.8.0 they were read from the row's start, three bases off past the junk.
                cut, motif = MOTIF[side]
                text = s["dna"].replace("-", "")
                shown = list(range(len(text)))[cut]  # the PAM's bases, as shown
                fwd_idx = sorted(j if s["strand"] == "+" else len(text) - 1 - j for j in shown)
                q, d = s["pos1"] - s["start0"], len(s["ref"]) - len(s["alt"])
                ref_pam = "".join(p.ref[s["start0"] + i + (d if i >= q else 0)] for i in fwd_idx)
                ref_pam = ref_pam if s["strand"] == "+" else rc(ref_pam)
                created = any(b not in BASES[m] for b, m in zip(ref_pam, motif))
                self.assertEqual(c[10], s["dna"][cut] if created else "n", where + ": PAM_gen")
                self.assertEqual(c[12], s["sample"], where)

    def run_inside(self, side, gate="1"):
        """A target whose first base is one an insertion adds starts at the insertion's
        anchor (its VCF POS), its REF text the reference from there, and its PAM_gen is
        the ALT's PAM when the PAM is wholly inserted bases, or, partly inserted, when the
        reference PAM beside the protospacer lacks the motif; a PAM on reference bases
        is read through the insertion. Up to 2.8.1 such a row started at the first base's
        place on the haplotype, past the anchor, and its PAM_gen read the REF text's
        window from there."""
        p = _PlantedInsertions(side)
        got = self.analyze(p, side, gate)
        cut, motif = MOTIF[side]
        kp = len(motif)
        for s in p.sites:
            token = "%s_%d_%s_%s" % (CHROM, s["pos1"], s["ref"], s["alt"])
            n, m, a = len(s["fwd"]), s["m"], s["anchor"]
            ref = p.ref[a:a + n]
            ref = ref if s["strand"] == "+" else rc(ref)
            # the PAM's forward indices on the target, and the reference PAM
            right = (s["strand"] == "+") == (side == "3p")
            idx = list(range(n - kp, n)) if right else list(range(kp))
            inserted = [i < m for i in idx]
            refpos = lambda i: a + 1 + i - m  # a forward reference base of the target
            if all(inserted):
                ref_pam = None
            elif any(inserted):
                sp = n - kp - 1 if right else kp  # the protospacer's base by the PAM
                if sp < m:
                    ref_pam = None
                else:
                    lo = refpos(sp) + 1 if right else refpos(sp) - kp
                    ref_pam = p.ref[lo:lo + kp]
            else:
                ref_pam = p.ref[refpos(idx[0]):refpos(idx[-1]) + 1]
            if ref_pam is not None and s["strand"] == "-":
                ref_pam = rc(ref_pam)
            created = ref_pam is None or any(b not in BASES[x] for b, x in zip(ref_pam, motif))
            for suffix, by_token in got.items():
                c = by_token[token]
                where = "%s inside %s(%s) m=%d %s" % (side, suffix, s["strand"], m, c[2])
                self.assertEqual(int(c[4]), a, where + ": Start_coordinate (the anchor)")
                self.assertEqual(int(c[5]), a, where + ": Cluster_Position")
                self.assertEqual(c[18].upper(), ref, where + ": REF text")
                self.assertEqual(c[10], s["shown"][cut] if created else "n", where + ": PAM_gen")

    def run_duplications(self, side, gate="1"):
        """A target whose PAM an insertion adds right after its protospacer, which is
        reference bases from the row's first base on: the row starts at that base (it does
        not move), its REF text is the reference from there, and its PAM_gen is NA where
        the insertion copies the reference's own PAM (REF = ALT) or the reference holds
        another valid PAM beside the protospacer, the ALT's PAM where it holds none, as up
        to 2.8.1, whose REF alignment read those bases. 2.10.0's first inserted-PAM rule
        took every such PAM for created (TRAC, HGDP + 1000 Genomes: chr22:20586277 (+),
        REF = ALT)."""
        p = _PlantedDuplications(side)
        got = self.analyze(p, side, gate)
        for s in p.sites:
            token = "%s_%d_%s_%s" % (CHROM, s["pos1"], s["ref"], s["alt"])
            n = len(s["fwd"])
            ref = p.ref[s["start0"]:s["start0"] + n]
            ref = ref if s["strand"] == "+" else rc(ref)
            for suffix, by_token in got.items():
                c = by_token[token]
                where = "%s duplication %s(%s) %s" % (side, suffix, s["strand"], c[2])
                self.assertEqual(c[2], s["shown"], where)
                self.assertEqual(int(c[4]), s["start0"], where + ": Start_coordinate")
                self.assertEqual(int(c[5]), s["start0"], where + ": Cluster_Position")
                self.assertEqual(c[18].upper(), ref, where + ": REF text")
                self.assertEqual(c[10], s["pam_gen"], where + ": PAM_gen")
                self.assertEqual(c[12], s["sample"], where)

    def run_placements(self, side, gate="1"):
        """2.10.0: a PAM holding inserted bases, or a deletion's junction, is judged wherever
        in its repeat the VCF writes the indel, through the shipped analyzer: each planted
        row starts where it should and gets the PAM_gen of _placed_sites, which 2.10.0's
        first fixes did not give it (a repeat insertion moved off the PAM, the reference PAM
        beside the protospacer on a row that starts on a reference base or inside an
        insertion, a deletion's junction on the - strand, Cas12a on both strands)."""
        p = _PlantedPlacements(side)
        got = self.analyze(p, side, gate)
        for s in p.sites:
            token = "%s_%d_%s_%s" % (CHROM, s["pos1"], s["ref"], s["alt"])
            for suffix, by_token in got.items():
                c = by_token[token]
                where = "%s %s %s(%s) %s" % (side, s["name"], suffix, s["strand"], c[2])
                self.assertEqual(c[2], s["shown"], where)
                self.assertEqual(int(c[4]), s["start"], where + ": Start_coordinate")
                self.assertEqual(c[10], s["pam_gen"], where + ": PAM_gen")
                self.assertEqual(c[12], s["sample"], where)

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_pams_judged_wherever_the_vcf_writes_the_indel_three_prime_pam(self):
        self.run_placements("3p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_pams_judged_wherever_the_vcf_writes_the_indel_five_prime_pam(self):
        """Cas12a (TTTV): the PAM at the target's left end on the + strand, forward, and at
        its right end on the - strand."""
        self.run_placements("5p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_pams_judged_wherever_the_vcf_writes_the_indel_gate_off(self):
        self.run_placements("3p", gate="0")
        self.run_placements("5p", gate="0")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_an_inserted_pam_on_a_row_that_does_not_move_three_prime_pam(self):
        self.run_duplications("3p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_an_inserted_pam_on_a_row_that_does_not_move_five_prime_pam(self):
        """Cas12a (TTTV): on the - strand the PAM is the target's right end, forward."""
        self.run_duplications("5p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_an_inserted_pam_on_a_row_that_does_not_move_gate_off(self):
        self.run_duplications("3p", gate="0")
        self.run_duplications("5p", gate="0")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_rows_starting_inside_an_insertion_three_prime_pam(self):
        self.run_inside("3p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_rows_starting_inside_an_insertion_five_prime_pam(self):
        """Cas12a (TTTV): on the + strand the PAM, at the target's 5' end, is the first
        base; on the - strand the protospacer's 3' end."""
        self.run_inside("5p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_rows_starting_inside_an_insertion_gate_off(self):
        self.run_inside("3p", gate="0")
        self.run_inside("5p", gate="0")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_three_prime_pam(self):
        self.run_side("3p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_five_prime_pam(self):
        self.run_side("5p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_three_prime_pam_gate_off(self):
        self.run_side("3p", gate="0")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_five_prime_pam_gate_off(self):
        self.run_side("5p", gate="0")


if __name__ == "__main__":
    unittest.main(verbosity=2)
