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
gapped as its target.

Before 2.8.1 an RNA-bulge row took the cluster position on the + strand with a 3' PAM
and on the - strand with a 5' PAM, so it started one base before its first aligned base,
and the REF text of an RNA,DNA row, and of a 5' PAM RNA-bulge row, was a base off.

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


class RowCoordinates(unittest.TestCase):
    def run_side(self, side):
        tmp = tempfile.mkdtemp(prefix="indel_row_coordinates_")
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        p = _Planted(side)
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
                # PAM_gen: the ALT's PAM when the reference there lacks the motif
                cut, motif = MOTIF[side]
                ref_pam = "".join(ref)[cut]
                created = any(b not in BASES[m] for b, m in zip(ref_pam, motif))
                self.assertEqual(c[10], s["dna"][cut] if created else "n", where + ": PAM_gen")
                self.assertEqual(c[12], s["sample"], where)

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_three_prime_pam(self):
        self.run_side("3p")

    @unittest.skipUnless(_gnu_wc(), _SKIP_WC)
    def test_five_prime_pam(self):
        self.run_side("5p")


if __name__ == "__main__":
    unittest.main(verbosity=2)
