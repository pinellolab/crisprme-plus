"""PAM-creation on the INDEL path (analisi_indels_NNN.py).

The block under test is inline in the analyzer's main loop, not a function, and
the module's top level reads sys.argv and opens the fake genome -- so, following
test_indel_snp_overlay.py, the source is extracted and exec'd against a synthetic
namespace instead of imported.

Why this file exists: the first version of the fix wrote the created PAM to
final_result[11] instead of [10]. The predicate itself was right, so every
standalone check of "which rows should be flagged" passed, and only a full
end-to-end run exposed it. These tests pin the *destination index* as well as the
predicate, so that class of error is caught in milliseconds.

Written as a unittest.TestCase so the project's `python -m unittest` CI harness
collects it (a prior pytest-style module-level `def test_*` layout collected zero
tests under unittest).

The reference bases the predicate reads are those the row's PAM sits on, through the
indel (_reference_pam). Up to 2.8.0 they were the REF alignment's, read from the row's
start, so on the + strand with a 3' PAM (the - strand with a 5' one) a deletion or an
insertion in the protospacer read bases len(REF) - len(ALT) away, and a PAM the reference
already has was flagged as created: the WALK panel's rows chr2:10000 (rsW_J_del) and
chr3:20000 (rsW_N_del + rsW_N14) showed PAM_creation TGG over a reference TGG.

Since 2.8.2 a PAM holding bases an insertion adds follows the inserted-PAM rules (wholly
inserted: created; partly: the reference PAM beside the protospacer) on a row that starts
inside the insertion, at its anchor, and keeps the REF alignment's reading, as up to 2.8.1,
on a row that starts on a reference base: on a tandem duplication of a PAM's bases the REF
alignment is the ALT, and the PAM is not created. A row whose REF alignment is its ALT never
has a created PAM. The first version of those rules took a wholly inserted PAM for created
on every row (TRAC, HGDP + 1000 Genomes: chr22:20586277 (+), REF = ALT).
"""

import os
import unittest

import indel_snp_cis

HERE = os.path.dirname(__file__)
SRC = os.path.join(HERE, "analisi_indels_NNN.py")

# The real indel cluster row: 18 fields, NO Chromosome_fake column.
# 0 Bulge_type, 1 crRNA, 2 DNA(ALT target), 3 Chromosome, 4 Position,
# 5 Cluster Position, 6 Direction, 7 Mismatches, 8 Bulge_Size, 9 Total,
# 10 PAM_gen, 11 Var_uniq, 12 Samples, 13 Annotation Type, 14 Real Guide,
# 15 rsID, 16 AF, 17 indel info
PAM_GEN = 10
VAR_UNIQ = 11

IUPAC = {
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"},
    "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "A": {"A"}, "T": {"T"}, "C": {"C"}, "G": {"G"},
    "a": {"a"}, "t": {"t"}, "c": {"c"}, "g": {"g"},
    "N": {"A", "T", "G", "C"},
}


def _block():
    """Extract _reference_pam and the PAM-creation block of _queue_row (ref_pam
    assignment through the loop)."""
    with open(SRC) as fh:
        src = fh.read()
    a = src.index("def _reference_pam(")
    func = src[a:src.index("\ndef _queue_row(", a)]
    start = src.index("    ref_pam = _reference_pam(final_result, ref_aligned, hap_start)")
    end = src.index("    final_result.append(ref_aligned)", start)
    body = src[start:end]
    # dedent one level: the block lives inside _queue_row
    return func + "\n" + "\n".join(l[4:] if l.startswith("    ") else l
                                   for l in body.splitlines())


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq):
    return seq.translate(_COMP)[::-1]


def _row(alt_target, token="chr1_100_A_AT", start=100, strand="+"):
    r = ["RNA,DNA", "GUIDE", alt_target, "chr1", str(start), str(start), strand,
         "2", "0", "2", "n", "n", "SAMPLE1", "n", "GUIDE", "rs1", "0.01", token]
    assert len(r) == 18
    return r


def _ns(final_result, ref_aligned, pam="NRG", pam_begin=-3, pam_end=None, genome=None,
        five=False, hap_start=None):
    if genome is None:
        # the row's insertion at 100 puts the PAM of a 23-nt target one base past the
        # REF alignment's: give the reference the window's PAM there too
        genome = "A" * 119 + ref_aligned[pam_begin:pam_end].upper() + "A" * 40
    ns = {"iupac_code_set": IUPAC, "final_result": final_result,
          "ref_aligned": ref_aligned, "pam": pam, "pam_begin": pam_begin,
          "pam_end": pam_end, "pam_at_beginning": five, "genomeStr": genome,
          "reverse_complement_table": _revcomp, "hap_start": hap_start,
          "_isc": indel_snp_cis}
    exec(compile(_block(), "<pam_creation>", "exec"), ns)
    return ns


def _run(alt_target, ref_aligned, pam="NRG", pam_begin=-3, pam_end=None, five=False,
         genome=None, hap_start=None, **row):
    return _ns(_row(alt_target, **row), ref_aligned, pam, pam_begin, pam_end, genome,
               five, hap_start)["final_result"]


class IndelPamCreation(unittest.TestCase):

    def test_fires_into_pam_gen_and_leaves_var_uniq_alone(self):
        # reference window GAT is not NRG (position 2: R={A,G} vs T) -> creation
        r = _run("ACGTACGTACGTACGTACGTTGG", "acgtacgtacgtacgtacgtGAT")
        self.assertEqual(r[PAM_GEN], "TGG")
        # the original bug wrote here instead; Var_uniq must be untouched
        self.assertEqual(r[VAR_UNIQ], "n")

    def test_noop_when_reference_pam_already_valid(self):
        # CAG IS a valid NRG (N=C, R=A, G=G) -> not a creation event
        r = _run("ACGTACGTACGTACGTACGTAGG", "acgtacgtacgtacgtacgtCAG")
        self.assertEqual(r[PAM_GEN], "n")
        self.assertEqual(r[VAR_UNIQ], "n")

    def test_second_real_window_from_live_run(self):
        # GCA -> AGG, observed in the instrumented complete-test run
        r = _run("ACGTACGTACGTACGTACGTAGG", "acgtacgtacgtacgtacgtGCA")
        self.assertEqual(r[PAM_GEN], "AGG")

    def test_dash_in_pam_window_is_skipped_not_raised(self):
        # a bulge dash inside the PAM window: '-' is not in the IUPAC table
        r = _run("ACGTACGTACGTACGTACGTT-G", "acgtacgtacgtacgtacgtG-T")
        self.assertEqual(r[PAM_GEN], "n")

    def test_short_window_is_skipped(self):
        r = _run("ACGTTGG", "AC")  # ref_aligned shorter than the PAM
        self.assertEqual(r[PAM_GEN], "n")

    def test_lowercase_reference_window_is_not_a_false_positive(self):
        # lowercase IUPAC keys map to lowercase sets, so without .upper() this
        # would wrongly score as a creation event
        r = _run("ACGTACGTACGTACGTACGTAGG", "acgtacgtacgtacgtacgtcag")
        self.assertEqual(r[PAM_GEN], "n")

    def test_five_prime_pam_window(self):
        # Cas12a-style leading PAM: TTTN at [0:4], + strand, its first base the row's
        # first, the T that A>AT at 100 adds (at 100 on the haplotype): the row starts at
        # the insertion's anchor, 99 (hap_start 100). The reference PAM is the four
        # reference bases before the protospacer (its first base at 104 on the haplotype is
        # reference 103), 99..102: AAAA here, not TTTN, so created
        r = _run("TTTAACGTACGTACGTACGTACGT", "GGGAacgtacgtacgtacgtacgt",
                 pam="TTTN", pam_begin=0, pam_end=4, five=True, start=99, hap_start=100)
        self.assertEqual(r[PAM_GEN], "TTTA")
        # and where those four bases are TTTA the reference holds the PAM: not created
        genome = "A" * 99 + "TTTA" + "A" * 60
        r = _run("TTTAACGTACGTACGTACGTACGT", "GGGAacgtacgtacgtacgtacgt",
                 pam="TTTN", pam_begin=0, pam_end=4, five=True, genome=genome, start=99,
                 hap_start=100)
        self.assertEqual(r[PAM_GEN], "n")

    # -- the reference bases under the PAM are read through the indel ---------------

    GUIDE = "CTAACAGTTGCTTTTATCAC"
    # the WALK panel's genome (hg38_walk, Genomes/hg38): chr2 9990-10040, chr3 19990-20040
    WALK_CHR2 = "N" * 9990 + "TAACGCCCATCTAACAGTTGTCTTTTATCACTGGCGGTGATTGGGCTCAC"
    WALK_CHR3 = "N" * 19990 + "TGGCGTAGTACTACCAGTTGCATTTCATCACTGGAGTCGAAGTCTGATAT"

    def test_walk_j_a_protospacer_deletion_over_a_reference_pam(self):
        """chr2:10000 (+), rsW_J_del GT>G at 10010: the reference has TGG under the ALT
        PAM (10021-10023); the REF alignment, 23 bases from 10000, ends in CTG."""
        alt = self.GUIDE + "TGG"
        ref_aligned = "CTAACAGTTGtcTTTtatcaCTG"
        r = _run(alt, ref_aligned, genome=self.WALK_CHR2, token="chr2_10010_GT_G",
                 start=10000)
        self.assertEqual(r[PAM_GEN], "n")
        ns = _ns(_row(alt, "chr2_10010_GT_G", 10000), ref_aligned, genome=self.WALK_CHR2)
        self.assertEqual(ns["_reference_pam"](_row(alt, "chr2_10010_GT_G", 10000),
                                              ref_aligned), "TGG")

    def test_walk_n_a_deletion_and_a_snp_over_a_reference_pam(self):
        """chr3:20000 (+), rsW_N_del CA>C at 20011 with rsW_N14 C>T at 20016."""
        r = _run("CTAcCAGTTGCTTTTATCACTGG", "CTAcCAGTTGCaTTTcatcaCTG",
                 genome=self.WALK_CHR3, token="chr3_20011_CA_C,chr3_20016_C_T",
                 start=20000)
        self.assertEqual(r[PAM_GEN], "n")

    def test_a_protospacer_insertion_over_a_reference_pam(self):
        """+ strand: the guide lacks its 11th base on the reference, which an insertion
        puts back; the reference's TGG is one base before the REF alignment's PAM."""
        g = self.GUIDE
        genome = "N" * 50 + g[:10] + g[11:] + "TGG" + "TT" + "N" * 20
        token = "chr1_%d_%s_%s" % (50 + 10, g[9], g[9] + g[10])
        r = _run(g + "TGG", (g[:10] + g[11:] + "TGGT").lower(), genome=genome, token=token,
                 start=50)
        self.assertEqual(r[PAM_GEN], "n")
        # and where the reference has no PAM there, the insertion did create it
        genome = "N" * 50 + g[:10] + g[11:] + "TCA" + "TT" + "N" * 20
        r = _run(g + "TGG", (g[:10] + g[11:] + "TCAT").lower(), genome=genome, token=token,
                 start=50)
        self.assertEqual(r[PAM_GEN], "TGG")

    def test_a_protospacer_deletion_that_creates_the_pam(self):
        g = self.GUIDE
        genome = "N" * 50 + g[:10] + "C" + g[10:] + "TCA" + "N" * 20
        r = _run(g + "TGG", (g[:10] + "C" + g[10:] + "TC").lower(), genome=genome,
                 token="chr1_%d_%sC_%s" % (60, g[9], g[9]), start=50)
        self.assertEqual(r[PAM_GEN], "TGG")

    def test_minus_strand_rows_read_the_pam_at_their_start(self):
        """- strand, 3' PAM: the PAM is the forward bases at the row's start, before an
        indel inside the target: both readings agree."""
        g = self.GUIDE
        fwd = _revcomp(g[:10] + "C" + g[10:] + "TGG")  # CCA + the guide's complement
        genome = "N" * 50 + fwd + "N" * 20
        alt = g + "TGG"
        # the extra C sits at forward offset 3 + 9 = 12 from the start: GxC>G there
        x = 50 + 3 + 9
        token = "chr1_%d_%s_%s" % (x, genome[x - 1:x + 1], genome[x - 1])
        ref_aligned = _revcomp(genome[50:73])
        row = _row(alt, token, 50, "-")
        ns = _ns(row, ref_aligned, genome=genome)
        self.assertEqual(ns["_reference_pam"](row, ref_aligned), ref_aligned[-3:].upper())
        self.assertEqual(ns["final_result"][PAM_GEN], "n")

    def test_five_prime_pam_minus_strand_is_read_past_the_indel(self):
        """A 5' PAM (TTTV) on the - strand is the forward bases at the target's right
        end, past an indel in its protospacer."""
        spacer = "ACGTTGCAACGTTGCAACGTAC"
        disp_ref = "TTTA" + spacer[:8] + "G" + spacer[8:]   # the reference, shown on -
        genome = "N" * 50 + _revcomp(disp_ref) + "N" * 20
        alt = "TTTA" + spacer
        n = len(alt)
        # forward, the extra base is the complement of G at displayed offset 12
        x = 50 + len(disp_ref) - 1 - 12   # its 0-based forward position
        token = "chr1_%d_%s_%s" % (x, genome[x - 1:x + 1], genome[x - 1])
        ref_aligned = _revcomp(genome[50:50 + n])
        row = _row(alt, token, 50, "-")
        ns = _ns(row, ref_aligned, pam="TTTV", pam_begin=0, pam_end=4, genome=genome,
                 five=True)
        self.assertEqual(ns["_reference_pam"](row, ref_aligned), "TTTA")
        self.assertEqual(ns["final_result"][PAM_GEN], "n")
        self.assertNotEqual(ref_aligned[:4].upper(), "TTTA")  # the old reading's window

    def test_a_pam_holding_an_alt_base_or_a_junction_keeps_the_rows_window(self):
        g = self.GUIDE
        # an insertion of GG after the spacer makes the PAM: its bases are ALT bases
        genome = "N" * 50 + g + "TCA" + "N" * 20
        token = "chr1_%d_%s_%s" % (70, g[19], g[19] + "TG")
        r = _run(g + "TGG", g.lower() + "TCA", genome=genome, token=token, start=50)
        self.assertEqual(r[PAM_GEN], "TGG")
        # a deletion whose junction is inside the PAM: T|GG across it
        genome = "N" * 50 + g + "TAAGG" + "N" * 20
        token = "chr1_%d_TAA_T" % 71
        row = _row(g + "TGG", token, 50)
        ns = _ns(row, g.lower() + "TAA", genome=genome)
        self.assertEqual(ns["_reference_pam"](row, g.lower() + "TAA"), "TAA")
        self.assertEqual(ns["final_result"][PAM_GEN], "TGG")

    def test_a_row_starting_on_an_inserted_base_reads_its_pam_from_the_haplotype(self):
        """+ strand: the target's first six bases are the last six an insertion adds
        (TTCAT + the guide's first six), so the row starts at the insertion's anchor
        (indel_snp_cis.reference_start), not where the target's first base is on the
        haplotype. Its PAM, past the insertion, sits on the reference's TGG, found from
        that haplotype place (hap_start); read from the row's start it would be the
        guide's TAT, six bases before it, and the PAM would be flagged as created."""
        g = self.GUIDE
        x = "GATTACAGATTACA"
        genome = "N" * 50 + x + g[6:] + "TGG" + "N" * 20
        anchor = 50 + len(x) - 1  # 0-based
        token = "chr1_%d_%s_%s" % (anchor + 1, x[-1], x[-1] + "TTCAT" + g[:6])
        ref_aligned = genome[anchor:anchor + 23].lower()
        row = _row(g + "TGG", token, anchor)
        hap = anchor + 1 + 5  # the 6th added base, g[0]
        ns = _ns(row, ref_aligned, genome=genome, hap_start=hap)
        self.assertEqual(ns["_reference_pam"](row, ref_aligned, hap), "TGG")
        self.assertEqual(ns["final_result"][PAM_GEN], "n")
        self.assertEqual(ns["_reference_pam"](row, ref_aligned), g[14:17])  # the old reading
        ns = _ns(_row(g + "TGG", token, anchor), ref_aligned, genome=genome)
        self.assertEqual(ns["final_result"][PAM_GEN], "TGG")

    # -- a PAM made of bases an insertion adds (2.8.2) ----------------------------------

    def test_a_wholly_inserted_pam_was_created_by_the_insertion(self):
        """- strand, 3' PAM (the oracle's Z2): the target's leftmost bases, its PAM CCA
        (TGG on -), are the 3rd to 5th bases an insertion adds. The PAM has no reference
        bases: created, whatever the reference holds beside the anchor (up to 2.8.1 the
        REF alignment's window was read, from the row's start: on Z2 TGG on the reference
        hid it, and 2.8.2's first fix read it at the anchor)."""
        g = self.GUIDE
        fwd = _revcomp(g + "TGG")      # CCA + the guide's complement
        # the anchor and the reference after it: the REF alignment's window from the
        # anchor is no NRG (TGG: CCA on -), then one (CTA: TAG on -), which a reading
        # from the row's start would have taken for the reference's own PAM
        for anchor, beside in (("T", "GGT"), ("C", "TAG")):
            genome = "N" * 50 + anchor + beside + fwd[9:] + "N" * 20
            token = "chr1_51_%s_%sGA%s" % (anchor, anchor, fwd[:9])
            ref_aligned = _revcomp(genome[50:73])
            row = _row(g + "TGG", token, 50, "-")
            hap = 50 + 1 + 2           # the 3rd added base
            ns = _ns(row, ref_aligned, genome=genome, hap_start=hap)
            self.assertIsNone(ns["_reference_pam"](row, ref_aligned, hap))
            self.assertEqual(ns["final_result"][PAM_GEN], "TGG")
        # + strand: a target wholly inside a long insertion (the oracle's I)
        token = "chr1_51_T_TCAT" + g + "TGGCAT"
        genome = "N" * 50 + "T" + "TGGAAAACCCC" + "N" * 40
        row = _row(g + "TGG", token, 50, "+")
        ns = _ns(row, genome[50:73].lower(), genome=genome, hap_start=50 + 4)
        self.assertEqual(ns["final_result"][PAM_GEN], "TGG")

    def test_a_partly_inserted_pam_reads_the_reference_pam_beside_the_protospacer(self):
        """- strand, 3' PAM (the oracle's Z4, Z5; chr22:21052639 of the README search):
        the PAM's first two forward bases are the last two an insertion (T>TGTCC) adds,
        its third and the protospacer reference bases. The reference PAM is the one beside
        the protospacer: the three reference bases before its first forward base, here
        the base before the anchor, the anchor T and the PAM's A. CTA (TAG on -) holds
        NRG: not created; ATA (TAT) does not: created. Read from the row's start, the
        anchor (TAG, CTA on -), both would read as created."""
        g = self.GUIDE
        fwd = _revcomp(g + "TGG")
        for before, made in (("C", "n"), ("A", "TGG")):
            genome = "N" * 48 + before + "T" + fwd[2:] + "N" * 20
            token = "chr1_50_T_TGT" + fwd[:2]
            anchor = 49
            ref_aligned = _revcomp(genome[anchor:anchor + 23])
            row = _row(g + "TGG", token, anchor, "-")
            hap = anchor + 1 + 2       # the 3rd added base
            ns = _ns(row, ref_aligned, genome=genome, hap_start=hap)
            self.assertEqual(ns["_reference_pam"](row, ref_aligned, hap),
                             _revcomp(before + "TA"))
            self.assertEqual(ns["final_result"][PAM_GEN], made)
            self.assertEqual(ref_aligned[-3:].upper(), "CTA")  # the row start's window

    def test_a_pam_starting_inside_the_insertion_that_ends_the_protospacer(self):
        """+ strand, 3' PAM, a row starting inside an insertion that holds the protospacer
        (from its first base, the insertion's 2nd) and the PAM's first base, the PAM's
        other two being reference bases: the PAM has no place on the reference beside the
        protospacer, so the insertion created it, whatever the reference holds there (AGG
        here, which the REF alignment's window, read from the anchor, does not reach)."""
        g = self.GUIDE
        for after in ("TCA", "AGG"):
            genome = "N" * 50 + "A" + "GG" + after + "N" * 20
            token = "chr1_51_A_AC%sT" % g     # adds C, the guide and T
            ref_aligned = genome[50:73].lower()
            row = _row(g + "TGG", token, 50, "+")
            hap = 50 + 2                      # the 2nd added base, the guide's first
            ns = _ns(row, ref_aligned, genome=genome, hap_start=hap)
            self.assertIsNone(ns["_reference_pam"](row, ref_aligned, hap))
            self.assertEqual(ns["final_result"][PAM_GEN], "TGG")

    # -- a row that starts on a reference base keeps 2.8.1's reading (2.8.2) -------------

    def test_a_tandem_duplication_of_the_pam_on_a_row_that_does_not_move(self):
        """+ strand, 3' PAM: G>GTGG after the protospacer copies the three bases after it,
        the reference's own TGG, so on the haplotype the target's PAM is the inserted copy,
        bases the insertion adds, and the row starts on its first base, a reference one: it
        does not move. Its REF alignment, the reference from that start, is its ALT: the
        reference holds this very target there, with the same PAM, and the insertion
        created nothing. NA, as up to 2.8.1; 2.8.2's first version of the inserted-PAM rule
        took the PAM for created (TGG) wherever the row starts (TRAC, HGDP + 1000 Genomes:
        chr22:20586277 (+), alternative alignments with REF = ALT)."""
        g = self.GUIDE
        genome = "N" * 50 + g + "TGG" + "CATTAC" + "N" * 20
        token = "chr1_70_%s_%sTGG" % (g[19], g[19])
        ref_aligned = genome[50:73]
        self.assertEqual(ref_aligned, g + "TGG")
        for hap in (None, 50):  # hap_start: the row's start
            row = _row(g + "TGG", token, 50, "+")
            ns = _ns(row, ref_aligned, genome=genome, hap_start=hap)
            self.assertEqual(ns["_reference_pam"](row, ref_aligned, hap), "TGG")
            self.assertEqual(ns["final_result"][PAM_GEN], "n")

    def test_an_inserted_pam_on_a_row_that_does_not_move_reads_the_ref_alignment(self):
        """+ strand, 3' PAM: G>GTGG after the protospacer, the reference's next bases not
        TGG. The row starts on a reference base, and its PAM, the bases the insertion adds,
        is judged on its REF alignment, the reference from that start, as up to 2.8.1:
        AGG there is a valid PAM at the PAM's place beside the protospacer (not created,
        though REF and ALT differ); TCA is not (created)."""
        g = self.GUIDE
        token = "chr1_70_%s_%sTGG" % (g[19], g[19])
        for after, made in (("AGG", "n"), ("TCA", "TGG")):
            genome = "N" * 50 + g + after + "CATTAC" + "N" * 20
            ref_aligned = genome[50:73]
            row = _row(g + "TGG", token, 50, "+")
            ns = _ns(row, ref_aligned, genome=genome, hap_start=50)
            self.assertEqual(ns["_reference_pam"](row, ref_aligned, 50), after)
            self.assertEqual(ns["final_result"][PAM_GEN], made)

    def test_a_tandem_duplication_of_a_five_prime_pam_on_the_minus_strand(self):
        """- strand, 5' PAM (TTTV; the Cas12a oracle's C6): the PAM is the target's right
        end, forward, and an insertion after the protospacer's last forward base copies the
        reference's TAAA (TTTA on -) there. The row starts on a reference base and its REF
        alignment is its ALT: not created. Where the reference holds TAAT (ATTA) after the
        protospacer, the inserted TTTA was created (the oracle's C4)."""
        spacer = "ACGTTGCAACGTTGCAACGTAC"
        alt = "TTTA" + spacer
        fwd = _revcomp(alt)              # the protospacer's complement, then TAAA
        n = len(alt)
        for after, made in (("TAAA", "n"), ("TAAT", "TTTA")):
            genome = "N" * 50 + fwd[:-4] + after + "CCGT" + "N" * 20
            anchor = 50 + n - 4 - 1      # the protospacer's last forward base, 0-based
            token = "chr1_%d_%s_%sTAAA" % (anchor + 1, genome[anchor], genome[anchor])
            ref_aligned = _revcomp(genome[50:50 + n])
            row = _row(alt, token, 50, "-")
            ns = _ns(row, ref_aligned, pam="TTTV", pam_begin=0, pam_end=4, genome=genome,
                     five=True, hap_start=50)
            self.assertEqual(ns["_reference_pam"](row, ref_aligned, 50), _revcomp(after))
            self.assertEqual(ns["final_result"][PAM_GEN], made)

    def test_a_partly_inserted_pam_on_a_row_that_does_not_move_reads_the_ref_alignment(self):
        """The partly-inserted rule (the reference PAM beside the protospacer) is for rows
        that start inside an insertion; a row starting on a reference base keeps its REF
        alignment's window, as up to 2.8.1.
          * + strand, 3' PAM: the insertion holds the protospacer's last two bases and the
            PAM's first; the window holds TCA (created) or AGG (not created).
          * - strand, 3' PAM, the row starting on the insertion's anchor C, the PAM's far
            end (forward CCA: the anchor and the two bases C>CCA adds): the window is the
            REF alignment's, read from the anchor (ACG on the - strand: created), not the
            reference PAM beside the protospacer (forward CTC, GAG on the - strand)."""
        g = self.GUIDE
        token = "chr1_68_%s_%s" % (g[17], g[17] + g[18:] + "T")
        for after, made in (("TCA", "TGG"), ("AGG", "n")):
            genome = "N" * 50 + g[:18] + "GG" + after + "N" * 20
            ref_aligned = genome[50:73].lower()
            row = _row(g + "TGG", token, 50, "+")
            ns = _ns(row, ref_aligned, genome=genome, hap_start=50)
            self.assertEqual(ns["_reference_pam"](row, ref_aligned, 50), after)
            self.assertEqual(ns["final_result"][PAM_GEN], made)
        genome = "N" * 47 + "CT" + "C" + _revcomp(g) + "TTAA" + "N" * 20
        ref_aligned = _revcomp(genome[49:72])
        self.assertEqual(ref_aligned[-3:], "ACG")
        row = _row(g + "TGG", "chr1_50_C_CCA", 49, "-")
        ns = _ns(row, ref_aligned, genome=genome, hap_start=49)
        self.assertEqual(ns["_reference_pam"](row, ref_aligned, 49), "ACG")
        self.assertEqual(ns["final_result"][PAM_GEN], "TGG")

    def test_a_ref_alignment_equal_to_the_alt_is_not_created_on_a_row_that_moves(self):
        """+ strand, 3' PAM: an insertion after the guide's first base adds the guide and
        TGG, a copy of the 23 reference bases from that anchor on. A target starting on its
        first added base starts at the anchor, its PAM wholly inserted bases; its REF
        alignment, read from the anchor, is its ALT: the reference holds this very target
        at the row's start, and its PAM is not created (in every case NA when REF = ALT)."""
        g = self.GUIDE
        genome = "N" * 50 + g + "TGG" + "CATTAC" + "N" * 20
        token = "chr1_51_%s_%s%sTGG" % (g[0], g[0], g)
        ref_aligned = genome[50:73]
        row = _row(g + "TGG", token, 50, "+")
        ns = _ns(row, ref_aligned, genome=genome, hap_start=51)
        self.assertEqual(ns["_reference_pam"](row, ref_aligned, 51), "TGG")
        self.assertEqual(ns["final_result"][PAM_GEN], "n")
        # the same target where the reference from the anchor differs: created
        genome = "N" * 50 + g[0] + "A" + g[2:] + "TGG" + "CATTAC" + "N" * 20
        ref_aligned = genome[50:73]
        row = _row(g + "TGG", token, 50, "+")
        ns = _ns(row, ref_aligned, genome=genome, hap_start=51)
        self.assertIsNone(ns["_reference_pam"](row, ref_aligned, 51))
        self.assertEqual(ns["final_result"][PAM_GEN], "TGG")

    def test_a_partly_inserted_pam_after_a_reference_protospacer_plus_strand(self):
        """+ strand, 3' PAM: the protospacer ends at the insertion's anchor and the PAM is
        two added bases and a reference one. The row starts on a reference base and keeps
        its REF alignment's window, the three bases after the anchor, which here is the
        reference PAM beside the protospacer (no change in 2.8.2 for this geometry)."""
        g = self.GUIDE
        for after, made in (("AGG", "n"), ("GTA", "TGG")):
            genome = "N" * 50 + g + after + "N" * 20
            token = "chr1_70_%s_%s" % (g[19], g[19] + "TG")
            row = _row(g + "TGG", token, 50, "+")
            ns = _ns(row, genome[50:73].lower(), genome=genome)
            self.assertEqual(ns["_reference_pam"](row, genome[50:73].lower()), after)
            self.assertEqual(ns["final_result"][PAM_GEN], made)

    def test_source_targets_index_10_not_11(self):
        """Static guard: the assignment must name index 10.

        The sed header in analisi_indels_NNN.sh lists a Chromosome_fake column at 10,
        which invites re-deriving this as 11. The rows the analyzer reads have no such
        column (18 fields; verified by instrumenting a live run), and this function's
        own [12]/[15]/[16]/[17] assignments confirm the SNP layout.
        """
        with open(SRC) as fh:
            src = fh.read()
        self.assertIn("final_result[10] = final_result[2][pam_begin:pam_end].upper()", src)
        self.assertNotIn("final_result[11] = final_result[2][pam_begin:pam_end].upper()", src)


if __name__ == "__main__":
    unittest.main()
