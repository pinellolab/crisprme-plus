#!/usr/bin/env python
"""The real indel post-analysis (analisi_indels_NNN.py) on SNP-overlaid targets.

A hermetic, engine-free fixture: a small reference chromosome with planted SpCas9
(NRG) sites, a phased panel (indels and SNPs, cis / trans / unphased), the fake-indel
windows cut as CRISPRitz add-variants cuts them and SNP-overlaid by
overlay_indel_snps, the engine's targets on them (IUPAC codes in the DNA column, a
code matching when any allele fits), and the dictionaries a search reads: the
Tier-0 registry, the Tier-1 genotype store, the phased indel store and the indel log,
compiled with the repo's own writers. The SHIPPED analisi_indels_NNN.py is then run
on them and its rows checked:

  * every reported base comes from the listed carrier's own haplotype: the reported
    sequence (gaps removed) is on a phased slot of every CONFIRMED carrier, at the
    reported site, and every PUTATIVE carrier holds each allele it uses;
  * every PAM satisfies NRG, and no IUPAC code is left;
  * the planted sites give exactly the rows their carriers have (a PAM SNP in cis
    and in trans, a chimera only another sample's SNP would complete, a cis spacer
    SNP, an unphased carrier, a minus-strand PAM SNP, a second indel on the same
    slot: inside the target, inserted next to the indel's anchor, or overlapping
    the indel as part of one complex allele);
  * the SNP+indel co-occurrence companion names the SNPs a row requires, CONFIRMED
    or PUTATIVE, PAM-window SNPs included;
  * a classic dictionary install (no Tier-1 store: the per-sample SNP dictionary), a
    sites-only panel (no genotypes), CRISPRME_INDEL_SNP=0 and --max-total-edits.

On 2.8.0 / aa9c54b the chimera and the PAM-less rows (an R at the last PAM base read
as A) are reported, and a classic dictionary install leaves the codes in the rows.

Needs pandas and numpy (the analyzer imports them); intervaltree is stubbed when it
is not installed. Run: python -m unittest test_indel_overlaid_rows -v
"""

import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import overlay_indel_snps as ov
import tier0_registry as t0
import tier1_genotypes as t1
from bulge_cfd_penalty import cfd_score_bulge_aware

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYZER = os.path.join(HERE, "analisi_indels_NNN.py")
CHROM = "chrT"
VCF = "synth"
GUIDE = "ACTGAAATCTGTAAGCAGGC"
CODE = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
        frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M"}
ALLELES = {v: set(k) for k, v in CODE.items()}
COMP = str.maketrans("ACGTRYSWKMacgtryswkm", "TGCAYRSWMKtgcayrswmk")
SAMPLES = ["S%02d" % i for i in range(1, 13)]
SAMPLE_META = {s: ("synth", "POP", "female") for s in SAMPLES}


def _gnu_wc():
    try:
        return subprocess.run(["wc", "--version"], capture_output=True).returncode == 0
    except OSError:
        return False


_SKIP_WC = "analisi_indels_NNN.py reads the line count from GNU `wc -l`, as the pipeline image has"

_STUB_INTERVALTREE = '''
class Interval(tuple):
    """Hashed and compared by (begin, end), as intervaltree's Interval is."""
    def __new__(cls, begin, end, data=None):
        return tuple.__new__(cls, (begin, end, data))
    begin = property(lambda self: self[0])
    end = property(lambda self: self[1])
    data = property(lambda self: self[2])
    def __hash__(self):
        return hash((self[0], self[1]))
    def __eq__(self, other):
        return (self[0], self[1]) == (other[0], other[1])
    def __lt__(self, other):
        return (self[0], self[1]) < (other[0], other[1])

class IntervalTree(object):
    def __init__(self):
        self._iv = []
    def __setitem__(self, key, data):
        self._iv.append(Interval(key.start, key.stop, data))
    def __getitem__(self, point):
        return set(iv for iv in self._iv if iv.begin <= point < iv.end)
'''


def rc(s):
    return s.translate(COMP)[::-1]


class World(object):
    """The reference, the panel and the planted sites."""

    def __init__(self):
        rng = random.Random(28101)
        self.seq = [rng.choice("ACGT") for _ in range(1700)]
        self.snps = []    # (pos1, ref, alt, rsid, {sample: gt})
        self.indels = []  # (pos1, ref, alt, rsid, {sample: gt})
        g = GUIDE
        # A (+): a deletion of CAC restores a perfect target, PAM TGG; a G>A SNP at the
        # PAM's last base in cis (S02), in trans (S03) and with no indel (S04).
        a0 = 100
        self._put(a0, g[:10] + "CAC" + g[10:] + "TGG")
        self._indel("A", a0 + 10, g[9] + "CAC", g[9], {"S01": "1|0", "S02": "1|0", "S03": "1|0"})
        self._snp("A", a0 + 10 + 3 + 10 + 2 + 1, "G", "A", {"S02": "1|0", "S03": "0|1", "S04": "1|0"})
        # C (+): a chimera. With its deletion S05's haplotype is 4 mismatches + one
        # more at column 6 (C); a C>A SNP there would complete a 4-mismatch target,
        # but only S06 (no deletion) carries it.
        c0 = 300
        xc = "tgac" + g[4:6] + "C" + g[7:10]
        self._put(c0, xc.upper() + "T" + g[10:] + "AGG")
        self._indel("C", c0 + 10, xc[9].upper() + "T", xc[9].upper(), {"S05": "1|0"})
        self._snp("C", c0 + 6 + 1, "C", "A", {"S06": "0|1"})
        # D (+): a C>A spacer SNP the guide needs, in cis with the deletion in S07,
        # absent in S08.
        d0 = 500
        self._put(d0, g[:3] + "C" + g[4:10] + "G" + g[10:] + "CGG")
        self._indel("D", d0 + 10, g[9] + "G", g[9], {"S07": "1|0", "S08": "0|1"})
        self._snp("D", d0 + 3 + 1, "C", "G", {"S07": "1|0", "S09": "1|0"})
        # E (+): the same, unphased: S10 carries both, phase unknown -> PUTATIVE.
        e0 = 700
        self._put(e0, g[:3] + "C" + g[4:10] + "G" + g[10:] + "TGG")
        self._indel("E", e0 + 10, g[9] + "G", g[9], {"S10": "0/1", "S11": "0/1"})
        self._snp("E", e0 + 3 + 1, "C", "G", {"S10": "0/1"})
        # F (-): an insertion of TT restores the target on the minus strand; a C>T
        # SNP at the PAM's last base (forward C, first base of 'CCT') in cis (S11) and
        # in trans (S12).
        f0 = 900
        pre = rc(g[14:] + "AGG")
        post = rc(g[:12])
        self._put(f0, pre + post)
        self._indel("F", f0 + len(pre), pre[-1], pre[-1] + "TT", {"S11": "1|0", "S12": "1|0"})
        self._snp("F", f0 + 1, "C", "T", {"S11": "1|0", "S12": "0|1"})
        # G (+): a deletion plus an A>G SNP that CREATES the PAM (TGA -> TGG), in cis in
        # S03 and in trans in S04; the co-occurrence companion must name it.
        g0 = 1100
        self._put(g0, g[:10] + "C" + g[10:] + "TGA")
        self._indel("G", g0 + 10, g[9] + "C", g[9], {"S03": "1|0", "S04": "1|0"})
        self._snp("G", g0 + 10 + 1 + 10 + 2 + 1, "A", "G", {"S03": "1|0", "S04": "0|1"})
        # O (+), no SNP: a deletion restores a perfect target (S01, S02), but S01's
        # haplotype also carries a second deletion inside it (on the same slot), so it
        # does not have that target: the one-indel fake window cannot show it.
        o0 = 1250
        self._put(o0, g[:10] + "C" + g[10:] + "TGG")
        self._indel("O", o0 + 10, g[9] + "C", g[9], {"S01": "1|0", "S02": "0|1"})
        self._indel("O2", o0 + 3, g[2] + g[3], g[2], {"S01": "1|0"})
        # P (+), no SNP: inserting the PAM after the spacer's last base (the anchor) makes
        # a perfect target that ends inside the insertion; S09's slot also holds an A
        # inserted right before the anchor (P2), inside the target, S08's the other slot.
        p0 = 1400
        self._put(p0, g + "T" * 30)
        self._indel("P", p0 + 20, g[19], g[19] + "TGG", {"S08": "1|0", "S09": "1|0"})
        self._indel("P2", p0 + 19, g[18], g[18] + "A", {"S08": "0|1", "S09": "1|0"})
        # Q (+), no SNP: the same insertion; on S12's slot a deletion of the anchor and
        # the two bases before it (Q2) overlaps it, one complex allele split into two
        # records that the window does not hold. S10 carries both unphased: a
        # multi-allelic heterozygote, one allele per haplotype, so it stays listed.
        q0 = 1550
        self._put(q0, g + "T" * 30)
        self._indel("Q", q0 + 20, g[19], g[19] + "TGG", {"S07": "0|1", "S10": "0/1", "S12": "0|1"})
        self._indel("Q2", q0 + 17, g[16:20], g[16], {"S10": "0/1", "S12": "0|1"})
        self.ref = "".join(self.seq)

    def _put(self, start0, s):
        for i, c in enumerate(s):
            self.seq[start0 + i] = c

    def _indel(self, site, pos1, ref, alt, gts):
        assert "".join(self.seq[pos1 - 1:pos1 - 1 + len(ref)]) == ref, (site, ref)
        self.indels.append((pos1, ref, alt, "rs_%s_indel" % site, gts))

    def _snp(self, site, pos1, ref, alt, gts):
        assert self.seq[pos1 - 1] == ref, (site, pos1, ref, self.seq[pos1 - 1])
        self.snps.append((pos1, ref, alt, "rs_%s_snp" % site, gts))

    # -- the fake-indel genome, cut and overlaid as the build does --------------------
    def fake_windows(self):
        """[(indel, start_position, window)] -- 25 bp of plain reference, the ALT, then
        27 bp downstream (CRISPRitz add-variants), SNP-overlaid (overlay_indel_snps)."""
        enriched = list(self.ref)
        for pos1, ref, alt, _rs, _g in self.snps:
            enriched[pos1 - 1] = CODE[frozenset((ref, alt))]
        enriched = "".join(enriched)
        out = []
        for ind in self.indels:
            pos1, ref, alt = ind[:3]
            s = pos1 - 26
            sub = self.ref[s:s + 25] + alt + self.ref[s + 25 + len(ref):s + 25 + len(ref) + 27]
            out.append((ind, s, ov.overlay_sub_fasta(sub, s, ref, alt, enriched)))
        return out

    def haplotypes(self, sample):
        """The sample's two haplotypes as {slot: (sequence, real 0-based pos per base)};
        an unphased genotype's alt goes on slot 0 (the check below only reads phased
        carriers' slots)."""
        out = {}
        for h in (0, 1):
            variants = []
            for pos1, ref, alt, _rs, gts in self.snps + self.indels:
                gt = gts.get(sample)
                if gt is None:
                    continue
                toks = gt.replace("/", "|").split("|")
                if h < len(toks) and toks[h] == "1":
                    variants.append((pos1, ref, alt))
            seq, pos, i = [], [], 0
            for pos1, ref, alt in sorted(variants):
                p0 = pos1 - 1
                seq.extend(self.ref[i:p0])
                pos.extend(range(i, p0))
                seq.extend(alt)
                pos.extend([p0] + [None] * (len(alt) - 1) if len(alt) > len(ref)
                           else range(p0, p0 + len(alt)))
                i = p0 + len(ref)
            seq.extend(self.ref[i:])
            pos.extend(range(i, len(self.ref)))
            out[h] = ("".join(seq), pos)
        return out


def engine_targets(windows, max_mm=4):
    """The engine's X (no-bulge) NRG targets on the overlaid fake genome: a code
    matches the guide or the PAM when any of its alleles fits; mismatches lower
    case. Returns (fake chromosome, rows)."""
    fake, rows = [], []
    starts = []
    for _ind, _s, win in windows:
        if fake:
            fake.append("N")
        starts.append(sum(len(x) for x in fake))
        fake.append(win)
    fake = "".join(fake)

    def fits(base, want):
        return want in ALLELES.get(base.upper(), {base.upper()})

    for strand in "+-":
        for f0 in range(len(fake) - 22):
            seg = fake[f0:f0 + 23]
            if "N" in seg:
                continue
            shown = seg if strand == "+" else rc(seg)
            pam = shown[20:]
            if not (fits(pam[1], "A") or fits(pam[1], "G")) or not fits(pam[2], "G"):
                continue
            dna, mm = [], 0
            for g, b in zip(GUIDE, shown[:20]):
                if fits(b, g):
                    dna.append(b.upper())
                else:
                    dna.append(b.lower())
                    mm += 1
            if mm > max_mm:
                continue
            rows.append(["X", GUIDE + "NNN", "".join(dna) + pam.upper(), "fake" + CHROM,
                         str(f0), str(f0), strand, str(mm), "0", str(mm)])
    return fake, starts, rows


class _Fixture(object):
    """Write the world as a data folder and run the shipped analyzer on it."""

    def __init__(self, tmp, world, tiers=("registry", "tier1", "indelgt"), env=None,
                 max_total=None):
        self.tmp = tmp
        self.world = world
        d = os.path.join(tmp, "Dictionaries")
        os.makedirs(d)
        windows = world.fake_windows()
        fake, starts, rows = engine_targets(windows)
        self.rows_in = rows
        with open(os.path.join(tmp, CHROM + ".fa"), "w") as fh:
            fh.write(">%s\n%s\n" % (CHROM, world.ref))
        logd = os.path.join(d, "log_indels_" + VCF)
        os.makedirs(logd)
        with open(os.path.join(logd, "log%s.txt" % CHROM), "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            for k, ((pos1, ref, alt, rs, gts), s, win) in enumerate(windows):
                fh.write("%s_%d-%d_%d\t%s\t%s\t%s\t%s_%d_%s_%s\t%d,%d\t%s\n" % (
                    CHROM, s, s + len(win), k + 1, ",".join(sorted(gts)), rs,
                    "%.6g" % (sum(g.count("1") for g in gts.values()) / 24.0),
                    CHROM, pos1, ref, alt, starts[k], starts[k] + len(win),
                    world.ref[s:s + len(win)]))
        recs = [(p, r, a, rs, dict(g)) for p, r, a, rs, g in world.snps]
        if "registry" in tiers:
            rd = os.path.join(d, "registry_" + VCF)
            os.makedirs(rd)
            t0.compile_registry_panel(recs, SAMPLE_META, None, t0.autosomal_ploidy,
                                      os.path.join(rd, "reg_%s.bin" % CHROM),
                                      os.path.join(rd, "reg_%s.idx" % CHROM))
        if "tier1" in tiers:
            gd = os.path.join(d, "genotypes_" + VCF)
            os.makedirs(gd)
            t1.compile_genotypes([(p, a, dict(g), r) for p, r, a, _rs, g in world.snps],
                                 t1.build_sample_axis(SAMPLE_META),
                                 os.path.join(gd, "gt_%s.bin" % CHROM),
                                 os.path.join(gd, "gt_%s.idx" % CHROM))
        if "dict" in tiers:
            dd = os.path.join(d, "dictionaries_" + VCF)
            os.makedirs(dd)
            entries = {}
            for p, r, a, rs, g in world.snps:
                entries["%s,%d" % (CHROM, p)] = "%s;%s,%s;%s;%s" % (
                    ",".join("%s:%s" % kv for kv in sorted(g.items())), r, a, rs,
                    "%.6g" % (sum(x.count("1") for x in g.values()) / 24.0))
            with gzip.open(os.path.join(dd, "my_dict_%s.json.gz" % CHROM), "wt") as fh:
                json.dump(entries, fh)
        if "indelgt" in tiers:
            igd = os.path.join(d, "indel_genotypes_" + VCF)
            os.makedirs(igd)
            with gzip.open(os.path.join(igd, "gt_indel_%s.tsv.gz" % CHROM), "wt") as fh:
                for p, r, a, _rs, g in world.indels:
                    fh.write("%d_%s_%s\t%s\n" % (p, r, a, ",".join(
                        "%s:%s" % kv for kv in sorted(g.items()))))
        with open(os.path.join(tmp, "pam.txt"), "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        with open(os.path.join(tmp, "guides.txt"), "w") as fh:
            fh.write(GUIDE + "NNN\n")
        cluster = os.path.join(tmp, "targets.cluster.txt")
        with open(cluster, "w") as fh:
            fh.write("#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
                     "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\t"
                     "Annotation Type\tReal Guide\tID\tAF\tSNP\n")
            for r in rows:
                fh.write("\t".join(r + ["n", "n", "n", "n", GUIDE, "n", "n", "n"]) + "\n")
        self.out = os.path.join(tmp, "out", "fake" + CHROM)
        os.makedirs(os.path.dirname(self.out))
        argv = [sys.executable, ANALYZER, "", cluster, self.out, logd,
                os.path.join(tmp, "pam.txt"), "4", os.path.join(tmp, CHROM + ".fa"),
                os.path.join(tmp, "guides.txt"), "1", "1"]
        if max_total is not None:
            argv.append(str(max_total))
        run_env = {k: v for k, v in os.environ.items()
                   if k not in ("MAMBA_EXE", "CONDA_EXE", "CRISPRME_CONDA_EXE")}
        run_env["PATH"] = os.path.dirname(sys.executable) + ":/usr/bin:/bin"
        run_env.update(env or {})
        try:
            import intervaltree  # noqa: F401
        except ImportError:
            stub = os.path.join(tmp, "stub")
            os.makedirs(stub)
            with open(os.path.join(stub, "intervaltree.py"), "w") as fh:
                fh.write(_STUB_INTERVALTREE)
            run_env["PYTHONPATH"] = stub + os.pathsep + run_env.get("PYTHONPATH", "")
        self.proc = subprocess.run(argv, cwd=HERE, env=run_env, capture_output=True,
                                   text=True, timeout=600)

    def rows(self):
        """The CFD-ranked output rows: (crRNA, DNA, start, strand, mm, samples, rsID, AF,
        SNP, CFD)."""
        out = []
        with open(self.out + ".bestCFD_INDEL.txt") as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                out.append({"crRNA": c[1], "dna": c[2], "start": int(c[4]), "strand": c[6],
                            "mm": int(c[7]), "total": int(c[9]), "pam_gen": c[10],
                            "samples": c[12], "rsid": c[15], "af": c[16], "var": c[17],
                            "cfd": float(c[20])})
        return out

    def cooc(self):
        path = self.out + ".indel_snp_cooc.tsv"
        if not os.path.exists(path):
            return None
        with open(path) as fh:
            next(fh)
            return sorted(set(tuple(line.rstrip("\n").split("\t")) for line in fh))


class _Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_overlaid_rows_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.world = World()

    def run_analyzer(self, **kw):
        fx = _Fixture(self.tmp, self.world, **kw)
        self.assertEqual(fx.proc.returncode, 0, fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
        return fx

    def by_site(self, rows):
        """{site letter: [rows]} by the indel token each row starts with."""
        sites = {}
        for r in rows:
            pos1 = int(r["var"].split(",")[0].split("_")[1])
            site = next(rs for p, _r, _a, rs, _g in self.world.indels if p == pos1)
            sites.setdefault(site.split("_")[1], []).append(r)
        return sites

    def assert_on_carrier_haplotypes(self, rows, putative_ok=()):
        """Every listed carrier holds the reported sequence: phased carriers on one
        of their slots, at the reported strand; for the samples in putative_ok each
        allele the row names must be one they carry."""
        for r in rows:
            seq = r["dna"].replace("-", "").upper()
            self.assertTrue(set(seq) <= set("ACGT"), r)
            self.assertIn(seq[-2:], ("AG", "GG"), r)
            for s in r["samples"].split(","):
                if s in putative_ok:
                    continue
                haps = self.world.haplotypes(s)
                found = any(seq in (h if r["strand"] == "+" else rc(h))
                            for h, _pos in haps.values())
                self.assertTrue(found, "%s does not carry %s (%s)" % (s, r["dna"], r["var"]))


@unittest.skipUnless(_gnu_wc(), _SKIP_WC)
class FullTiers(_Base):
    """A genotyped, phased panel: registry + Tier-1 + phased indel store."""

    def test_rows_are_the_carriers_own_haplotypes(self):
        fx = self.run_analyzer()
        rows = fx.rows()
        self.assertTrue(rows)
        self.assert_on_carrier_haplotypes(rows, putative_ok=("S10", "S11"))
        site = self.by_site(rows)
        g = GUIDE
        best = lambda rs: max(rs, key=lambda r: r["cfd"])  # noqa: E731
        # A: the PAM SNP in cis destroys S02's PAM; S01 and S03 keep TGG
        a = best(site["A"])
        self.assertEqual((a["dna"], a["samples"], a["cfd"]), (g + "TGG", "S01,S03", 1.0))
        self.assertFalse(any(r["dna"].endswith("A") for r in site["A"]), site["A"])
        # C: the chimera (S05's deletion + S06's SNP) is not reported
        self.assertNotIn("C", site, site.get("C"))
        # D: the cis SNP gives S07 a perfect target, named; S08 has one mismatch
        d = {r["samples"]: r for r in site["D"] if r["mm"] == 0 or r["samples"] == "S08"}
        self.assertEqual(d["S07"]["dna"], g + "CGG")
        self.assertEqual(d["S07"]["var"].split(","),
                         ["%s_510_TG_T" % CHROM, "%s_504_C_G" % CHROM])
        self.assertEqual(d["S07"]["rsid"], "rs_D_indel,rs_D_snp")
        self.assertEqual(len(d["S07"]["af"].split(",")), 2)
        self.assertEqual((d["S08"]["mm"], d["S08"]["var"]), (1, "%s_510_TG_T" % CHROM))
        # E: unphased -- the worst case of S10's own alleles, S11 never gets the SNP
        e = {r["samples"]: r for r in site["E"]}
        self.assertEqual(e["S10"]["dna"], g + "TGG")
        self.assertIn("S11", e)
        self.assertEqual(e["S11"]["dna"][3], "c")
        # F: the minus-strand PAM SNP in cis (S11) destroys the PAM
        f = best(site["F"])
        self.assertEqual((f["dna"], f["samples"], f["cfd"]), (g + "AGG", "S12", 1.0))
        # G: the PAM-creating SNP in cis gives S03 a perfect target; S04 has no PAM
        gg = site["G"]
        self.assertEqual([(r["dna"], r["samples"], r["cfd"]) for r in gg],
                         [(g + "TGG", "S03", 1.0)])
        self.assertEqual(gg[0]["var"].split(",")[1], "%s_1124_A_G" % CHROM)
        self.assertEqual(gg[0]["pam_gen"], "TGG")  # the reference has no PAM there
        # O: no SNP code, but S01's slot also holds the second deletion: S02 only
        self.assertEqual([(r["dna"], r["samples"]) for r in site["O"]], [(g + "TGG", "S02")])
        self.assertNotIn("O2", site, site.get("O2"))
        # P: the target ends inside the insertion, and S09's slot has an A inserted
        # next to the anchor: S08 only; S09 is never listed for P2's window either
        self.assertEqual([(r["dna"], r["samples"]) for r in site["P"]], [(g + "TGG", "S08")])
        self.assertFalse([r for r in site.get("P2", ()) if "S09" in r["samples"].split(",")])
        # Q: S12's slot holds the overlapping deletion: listed on neither window
        self.assertEqual([(r["dna"], r["samples"]) for r in site["Q"]],
                         [(g + "TGG", "S07,S10")])
        self.assertFalse([r for r in site.get("Q2", ()) if "S12" in r["samples"].split(",")])
        # every recomputed CFD and mismatch count agrees with its alignment
        for r in rows:
            self.assertAlmostEqual(cfd_score_bulge_aware(r["crRNA"], r["dna"]), r["cfd"],
                                   places=3)
            self.assertEqual(r["mm"], sum(1 for c in r["dna"][:20] if c.islower()))

    def test_the_cooccurrence_companion(self):
        fx = self.run_analyzer()
        cooc = fx.cooc()
        got = sorted(set((c[1], c[6], c[7], c[8], c[10], c[11]) for c in cooc))
        self.assertEqual(got, [
            ("1110", "1124", "rs_G_snp", "CONFIRMED", "1", "S03"),  # a PAM-window SNP
            ("510", "504", "rs_D_snp", "CONFIRMED", "1", "S07"),
            ("710", "704", "rs_E_snp", "PUTATIVE", "1", "S10"),
        ])

    def test_the_total_edits_cap(self):
        """With --max-total-edits 0 only the perfect haplotypes stay."""
        fx = self.run_analyzer(max_total=0)
        rows = fx.rows()
        self.assertTrue(rows)
        self.assertTrue(all(r["total"] == 0 for r in rows), rows)


@unittest.skipUnless(_gnu_wc(), _SKIP_WC)
class OtherInstalls(_Base):
    def test_classic_dictionary_install(self):
        """No Tier-1 store and no registry: the per-sample SNP dictionary gives the
        genotypes; the rows match the full tiers' rows."""
        full = self.run_analyzer().rows()
        tmp2 = tempfile.mkdtemp(prefix="indel_overlaid_rows_dict_")
        self.addCleanup(shutil.rmtree, tmp2, ignore_errors=True)
        fx = _Fixture(tmp2, self.world, tiers=("dict", "indelgt"))
        self.assertEqual(fx.proc.returncode, 0, fx.proc.stdout[-2000:] + fx.proc.stderr[-2000:])
        rows = fx.rows()
        key = lambda r: (r["dna"], r["start"], r["strand"], r["samples"], r["var"], r["cfd"])  # noqa: E731
        self.assertEqual(sorted(map(key, rows)), sorted(map(key, full)))
        self.assertIsNone(fx.cooc())  # the companion needs the registry, as before

    def test_sites_only_panel(self):
        """No genotypes at all: one PUTATIVE worst-case row per target, the log's
        carriers, a valid PAM, no code left."""
        fx = self.run_analyzer(tiers=("registry",))
        rows = fx.rows()
        self.assertTrue(rows)
        for r in rows:
            seq = r["dna"].replace("-", "").upper()
            self.assertTrue(set(seq) <= set("ACGT"), r)
            self.assertIn(seq[-2:], ("AG", "GG"), r)
        a = max(self.by_site(rows)["A"], key=lambda r: r["cfd"])
        self.assertEqual((a["dna"], a["samples"]), (GUIDE + "TGG", "S01,S02,S03"))
        self.assertTrue(all(c[8] == "PUTATIVE" and c[10] == "NA" for c in fx.cooc() or []))

    def test_opting_out_reads_the_reference(self):
        """CRISPRME_INDEL_SNP=0 on an overlaid genome: each code is the reference
        base, as on a fake-indel genome without the overlay."""
        fx = self.run_analyzer(env={"CRISPRME_INDEL_SNP": "0"})
        rows = fx.rows()
        site = self.by_site(rows)
        a = max(site["A"], key=lambda r: r["cfd"])
        self.assertEqual((a["dna"], a["samples"], a["var"]),
                         (GUIDE + "TGG", "S01,S02,S03", "%s_110_TCAC_T" % CHROM))
        for r in rows:
            self.assertEqual(len(r["var"].split(",")), 1, r)  # no SNP is named
        self.assertIsNone(fx.cooc())


if __name__ == "__main__":
    unittest.main(verbosity=2)
