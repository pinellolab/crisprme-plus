#!/usr/bin/env python
from math import trunc
from operator import index
import sys
import json
import gzip
import os
import pickle
import numpy as np
import pandas as pd
import time
from CRISTA_score import CRISTA_predict_list

# Tier-0 compact variant registry (CRISPRme+ dictless redesign). GUARDED /
# lazy: an old deploy may not ship tier0_registry, and most installs have no
# registry at all -- in either case ``t0_reg`` stays None and the legacy dict
# path is used unchanged. Import failure NEVER breaks the legacy path.
try:
    import tier0_registry as t0_reg
except Exception:  # module absent in an old deploy -> legacy path only
    t0_reg = None
# Tier-1 compact per-sample GENOTYPE tier (CRISPRme+ dictless redesign, Phase 2).
# GUARDED / lazy exactly like tier0_registry: an old deploy may not ship
# tier1_genotypes, and most installs have no genotype tier at all -- in either
# case ``t1_gt`` stays None and ``mygt`` stays None, so the Samples column is
# reconstructed from the per-sample dict (legacy/augment path) unchanged. Import
# failure NEVER breaks the legacy or registry paths.
try:
    import tier1_genotypes as t1_gt
except Exception:  # module absent in an old deploy -> no genotype tier
    t1_gt = None
# Pure, side-effect-free retrieveFromDict machinery (importable / unit-tested in
# isolation). This module has no hard dependency on tier0_registry, so its import
# cannot break the legacy path either.
import simple_analysis_registry as _sar
# ADDITIVE combination-aware population summary (Phase 3) + its companion-file writer
# (Phase 3c). GUARDED / lazy exactly like tier0_registry: an old deploy may not ship
# these -- in that case they stay None and the companion population-summary TSV is
# simply NOT written. Import failure NEVER breaks the legacy, registry, or genotype
# paths, and the companion file is ADDITIVE (a separate TSV next to the output stem,
# gated on ``myreg`` -- it never touches the bestMerge/altMerge/integrated_results
# columns or the Samples column).
try:
    import population_summary as _popsum
    import population_summary_companion as _popsum_companion
    import tier0_compile as _t0c
except Exception:  # any module absent in an old deploy -> no companion summary
    _popsum = None
    _popsum_companion = None
    _t0c = None
# ADDITIVE observed-haplotype enumerator (dict-less multi-variant path). GUARDED /
# lazy exactly like the tiers above: an old deploy may not ship it, and it is ONLY
# consulted on the ``mygt is not None`` dict-less branch. When absent, ``_obshap``
# stays None and iupac_decomposition takes the BYTE-IDENTICAL legacy dict path (dict
# decode + haplotype_check peel + greedy cap) unchanged. Its companion phase-
# confirmation TSV is ADDITIVE (a separate file next to the output stem) so it never
# touches the bestMerge/Samples columns or the downstream positional parsers.
try:
    import observed_haplotypes as _obshap
    import phase_confirmation_companion as _phase_companion
except Exception:  # module absent in an old deploy -> legacy dict path only
    _obshap = None
    _phase_companion = None
# Accumulator for the ADDITIVE phase-confirmation companion TSV (one row per emitted
# dict-less variant off-target: identity columns + CONFIRMED/PUTATIVE). Populated ONLY
# on the ``mygt is not None`` branch; dead/empty on every legacy install so the
# legacy path stays allocation-identical, not just byte-identical.
_phase_confirmation_rows = []
_phase_confirmation_keys = set()

# Module-level dictless state, DEFAULTED here in the import prologue so
# retrieveFromDict() / _collect_variant_off_target() / _write_population_summary_
# companion() always resolve these globals (they are read, not assigned, inside
# those functions). The real detection block further down (guarded on sys.argv,
# after inFasta is opened) RE-BINDS them to a live RegistryReader / GenotypeReader
# when a Tier-0 registry / Tier-1 genotype store is on disk; absent that, they stay
# None and every dictless hook is a byte-identical no-op on the legacy path. This
# prologue default also lets the pure-function unit tests (which AST-load only the
# code before ``inFasta = open(...)``) drive retrieveFromDict without the runtime
# module body -- e.g. test_phased_haplotype's iupac_decomposition round-trip.
myreg = None
mygt = None
# Registry-only install (Tier-0 registry present, NO per-sample dict, NO genotype
# tier -- e.g. a `--no-genotypes` download). Computed once after tier detection;
# when True the variant finalizer emits degraded off-target rows (empty Samples ->
# "NA") instead of dropping the site (#136). Default False keeps every other mode
# (legacy, registry+dict, dictless-with-genotypes) byte-identical.
registry_only_mode = False


# For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv", "TGACYRSWKMDHVBtgacyrswkmdhvb")

iupac_nucleotides = set("RYSWKMBDHVryswkmbdhv")
iupac_code_set = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "r": {"A", "G"},
    "y": {"C", "T"},
    "s": {"G", "C"},
    "w": {"A", "T"},
    "k": {"G", "T"},
    "m": {"A", "C"},
    "b": {"C", "G", "T"},
    "d": {"A", "G", "T"},
    "h": {"A", "C", "T"},
    "v": {"A", "C", "G"},
    "A": {"A"},
    "T": {"T"},
    "C": {"C"},
    "G": {"G"},
    "a": {"a"},
    "t": {"t"},
    "c": {"c"},
    "g": {"g"},
    "N": {"A", "T", "G", "C"},
}


iupac_code = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "r": ("A", "G"),
    "y": ("C", "T"),
    "s": ("G", "C"),
    "w": ("A", "T"),
    "k": ("G", "T"),
    "m": ("A", "C"),
    "b": ("C", "G", "T"),
    "d": ("A", "G", "T"),
    "h": ("A", "C", "T"),
    "v": ("A", "C", "G"),
    "N": ("A", "T", "C", "G"),
}


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


class reversor:
    """
    Nel caso debba ordinare più campi però con reverse diversi, eg uno True e l'altro False, posso usare questa classe nella chiave per
    simulare il contrario del reverse applicato
    """

    def __init__(self, obj):
        self.obj = obj

    def __eq__(self, other):
        return other.obj == self.obj

    def __lt__(self, other):
        return other.obj < self.obj


def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores, do_scores):
    # if do_scores == False:
    #     # print("NON CALCOLO")
    #     score = -1
    #     return score
    score = 1
    sg = sg.replace("T", "U")
    guide_seq = guide_seq.replace("T", "U")
    s_list = list(sg)
    guide_seq_list = list(guide_seq)

    for i, sl in enumerate(s_list):
        if guide_seq_list[i] == sl:
            score *= 1
        else:
            try:  # Catch exception if IUPAC character
                key = "r" + guide_seq_list[i] + ":d" + revcom(sl) + "," + str(i + 1)
            except Exception as e:
                score = 0
                break
            try:
                score *= mm_scores[key]
            except (
                Exception
            ) as e:  # If '-' is in first position, i do not have the score for that position
                pass

    # Guard against non-ATCG PAM bases (e.g. from real NIST/NCBI hg38):
    # a non-canonical PAM contributes 0 to the CFD product instead of raising
    # KeyError (issue #94). Canonical PAM bases (A/C/G/T) are unaffected.
    score *= pam_scores.get(pam, 0.0)
    return score


# ADDITIVE (Phase 3c) companion population-summary collection. We accumulate, per
# chromosome, one record per VARIANT off-target (identity columns + the creating SNP
# column) so a companion population-summary TSV can be written at the end WITHOUT
# perturbing the bestMerge output. Deduped by the bestMerge identity key so
# per-haplotype / per-score duplicate rows collapse to a single companion row. This
# collection is only ever CONSUMED when ``myreg`` is present (see the guarded write
# at the end); otherwise it is harmless dead data and no companion file is written.
_variant_off_targets = []
_variant_off_target_keys = set()

# Column indices in a finalized target line (see ``header`` below):
#   3=Chromosome  4=Position  6=Direction  1=crRNA  2=DNA  17=SNP
_OT_COL = {
    "Chromosome": 3,
    "Position": 4,
    "Direction": 6,
    "crRNA": 1,
    "DNA": 2,
    "SNP": 17,
}


def _collect_variant_off_target(final_line):
    """Record a VARIANT off-target's identity + creating-SNP columns for the
    companion population summary. PURE w.r.t. ``final_line`` (reads, never mutates).

    Guarded: any malformed line (short list, non-string SNP) is silently skipped so
    this ADDITIVE bookkeeping can never break the finalization loop. Deduped by the
    bestMerge identity key (Chromosome, Position, Direction, crRNA, DNA, SNP) so the
    per-haplotype / CFD-vs-CRISTA duplicate rows collapse to one companion row.

    Gated on the Tier-0 registry (and the companion module) being present: on a
    legacy / dict-only install the companion is never written, so accumulating this
    per-chromosome set+list would be dead RAM on the hot finalization loop -- which
    matters given this pipeline's OOM history on genome-wide searches. Returning
    early here keeps the legacy path ALLOCATION-identical, not just byte-identical.
    ``myreg`` (module scope) is bound before the target loop that drives this."""
    if myreg is None or _popsum_companion is None:
        return
    try:
        rec = {name: final_line[idx] for name, idx in _OT_COL.items()}
    except Exception:
        return
    key = (rec["Chromosome"], rec["Position"], rec["Direction"], rec["crRNA"],
           rec["DNA"], rec["SNP"])
    if key in _variant_off_target_keys:
        return
    _variant_off_target_keys.add(key)
    _variant_off_targets.append(rec)


def revcom(s):
    basecomp = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "-": "-"}
    letters = list(s[::-1])
    try:
        letters = [basecomp[base] for base in letters]
    except:
        return None  # If some IUPAC were not translated
    return "".join(letters)


def get_mm_pam_scores():
    # print(os.path.dirname(os.path.realpath(__file__)))
    try:
        mm_scores = pickle.load(
            open(
                os.path.dirname(os.path.realpath(__file__)) + "/mismatch_score.pkl",
                "rb",
            )
        )
        pam_scores = pickle.load(
            open(os.path.dirname(os.path.realpath(__file__)) + "/PAM_scores.pkl", "rb")
        )
        return (mm_scores, pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


def retrieveFromDict(chr_pos):
    """Return the 5-tuple (snp_list, sample_list, rsID_list, AF_list, snp_info_list)
    for a candidate off-target position.

    ADDITIVE Tier-0/Tier-1 wiring: the PURE selection-matrix logic lives in
    ``simple_analysis_registry.retrieve_5tuple`` (importable / unit-tested in
    isolation). Here we only supply the module-level registry reader (``myreg``,
    or None), the genotype-tier reader (``mygt``, or None), the raw ``mydict``
    entry for this position (or None if absent), and the current chromosome. When
    ``myreg`` is None this is BYTE-IDENTICAL to the legacy behavior; when a
    registry is present but the per-sample dict is ABSENT and ``mygt`` is present,
    the Samples column is reconstructed from the genotype tier (true dictless).
    The dict entry is detected with the SAME bare-subscript try/except as the
    original so "present" vs "absent" is decided identically."""
    # ``entry is None`` is the sentinel for "the dict had no entry for this
    # position" (legacy except-branch). A present entry (even a malformed one) is
    # passed through verbatim so the legacy decode -- including its IndexError on a
    # malformed value -- is unchanged.
    try:
        entry = mydict[current_chr + "," + str(chr_pos + 1)]
    except Exception:
        entry = None
    global_gid = t0_reg.GLOBAL_GROUP_ID if t0_reg is not None else "global"
    return _sar.retrieve_5tuple(
        myreg, entry, current_chr, chr_pos, global_gid, gtreader=mygt
    )


def _aligned_mm(seq_prerevert, realTarget, guide_no_pam, revert):
    """Mismatch count of an ungapped, pre-revert candidate sequence against the
    guide, mirroring the finalization below: reverse-complement on the minus
    strand, re-insert bulge gaps at realTarget's '-' positions, then compare the
    protospacer window to guide_no_pam. Used by the high-variant-density greedy
    representative (validated == brute-force argmin over all allele combinations)."""
    seq = (
        reverse_complement_table("".join(seq_prerevert))
        if revert
        else "".join(seq_prerevert)
    )
    t = list(seq)
    for pos, ch in enumerate(realTarget):
        if ch == "-":
            t.insert(pos, "-")
    window = t[pos_beg:pos_end]
    mm = 0
    for i, ch in enumerate(window):
        if i < len(guide_no_pam) and ch != "-" and ch.upper() != guide_no_pam[i]:
            mm += 1
    return mm


def _record_phase_confirmation(final_line, phase_state):
    """Record one dict-less variant off-target's identity + CONFIRMED/PUTATIVE flag
    for the ADDITIVE phase-confirmation companion TSV. PURE w.r.t. ``final_line``
    (reads, never mutates). Deduped by the SAME bestMerge identity key the population-
    summary companion uses, so per-haplotype / CFD-vs-CRISTA duplicates collapse to one
    row. GATED on the dict-less branch (mygt present); a no-op collection otherwise."""
    try:
        rec = {name: final_line[idx] for name, idx in _OT_COL.items()}
    except Exception:
        return
    key = (rec["Chromosome"], rec["Position"], rec["Direction"], rec["crRNA"],
           rec["DNA"], rec["SNP"])
    if key in _phase_confirmation_keys:
        return
    _phase_confirmation_keys.add(key)
    rec["Phase_Confirmation"] = phase_state
    _phase_confirmation_rows.append(rec)


def _finalize_observed_entry(split, realTarget, refSeq_prerevert,
                             refSeq_with_bulges, guide_no_pam, revert, seq_prerevert,
                             carriers, info, phase_state, cluster_to_save):
    """Materialize ONE observed haplotype into the finalized row + companion flag.

    Drives the EXACT same column writes (2/7/9/10/12/15/16/17), PAM check, PAM-
    creation, mm-threshold discard, Reference/33/tmp_pos appends and cluster append
    as the legacy finalization block (new_simple_analysis.py's totalDict path), so the
    scored output is byte-for-byte compatible with a dict-based variant row (same
    negative-index sentinels for CFD/CRISTA). The CONFIRMED/PUTATIVE flag is carried
    OUT OF BAND to the companion TSV -- it is NEVER appended to ``final_line`` (that
    would shift the target[-3]/target[-4] scoring sentinels).
    """
    if not carriers:
        return
    final_line = split.copy()

    # apply strand + re-insert bulge gaps, mirroring the legacy per-key block
    if revert:
        seq_joined = reverse_complement_table("".join(seq_prerevert))
    else:
        seq_joined = "".join(seq_prerevert)
    target_to_list = list(seq_joined)
    for pos, char in enumerate(realTarget):
        if char == "-":
            target_to_list.insert(pos, "-")

    mm_new_t = 0
    tmp_pos_mms = 0
    for position_t, char_t in enumerate(target_to_list[pos_beg:pos_end]):
        if char_t.upper() != guide_no_pam[position_t]:
            mm_new_t += 1
            tmp_pos_mms = position_t
            if guide_no_pam[position_t] != "-":
                target_to_list[pos_beg + position_t] = char_t.lower()

    pam_ok = True
    for pam_chr_pos, pam_chr in enumerate(target_to_list[pam_begin:pam_end]):
        if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
            pam_ok = False

    target_pam_ref = refSeq_with_bulges[pam_begin:pam_end]
    found_creation = False
    for pos_pam, pam_char in enumerate(target_pam_ref):
        if not iupac_code_set[pam[pos_pam]] & iupac_code_set[pam_char]:
            found_creation = True

    if mm_new_t - int(split[8]) > allowed_mms:
        return
    if not pam_ok:
        return

    final_line[2] = "".join(target_to_list)
    final_line[7] = str(mm_new_t - int(final_line[8]))
    final_line[9] = str(mm_new_t)
    if found_creation:
        final_line[10] = "".join(target_to_list[pam_begin:pam_end])
    final_line[12] = ",".join(carriers)
    tmp_matrix = np.array(info)
    if tmp_matrix.shape[0] > 1:
        final_line[15] = ",".join(tmp_matrix[:, 0])
        final_line[16] = ",".join(tmp_matrix[:, 1])
        final_line[17] = ",".join(tmp_matrix[:, 2])
    else:
        final_line[15] = str(tmp_matrix[0][0])
        final_line[16] = str(tmp_matrix[0][1])
        final_line[17] = str(tmp_matrix[0][2])
    final_line.append(refSeq_with_bulges)
    final_line.append(33)
    final_line.append(tmp_pos_mms)
    cluster_to_save.append(final_line)
    # ADDITIVE: companion phase-confirmation flag (out of band; never mutates
    # final_line so the bestMerge column count + scoring sentinels are unchanged) and
    # the ADDITIVE population-summary companion, both gated / deduped by identity.
    _record_phase_confirmation(final_line, phase_state)
    _collect_variant_off_target(final_line)


def _finalize_reference_entry(split, realTarget, refSeq_prerevert,
                              refSeq_with_bulges, guide_no_pam, revert,
                              cluster_to_save):
    """Materialize the REFERENCE off-target of an IUPAC candidate (dict-less path).

    LOCUS-COVERAGE FIX. Legacy ``iupac_decomposition`` ALWAYS represents a candidate's
    reference off-target: every alt row it emits embeds ``refSeq_with_bulges`` + the
    sentinel ``33``, and ``resultIntegrator`` resolves that embedded reference into the
    cluster's origin=ref alignment/score. The observed-haplotype path only emits real
    variant haplotypes (``enumerate_observed_haplotypes`` never yields the empty set),
    so a candidate with NO productive / in-budget observed haplotype lost its entire
    locus -- both its reference AND its alt off-targets. This helper restores the
    reference off-target for EVERY candidate, exactly as the SEPARATE reference-genome
    CRISPRitz search would (a pure ``n`` / sentinel-``55`` reference row), so the
    dict-less LOCUS coverage matches the stable reports.

    The row is a PURE reference alignment (no alt substituted). It carries the SAME
    column writes (2/7/9/10) and PAM check as ``_finalize_observed_entry``, but its
    Samples / rsID / AF / SNP markers are the LEGACY literal "n" (NOT "NA"), making the
    row byte-identical to the legacy non-IUPAC else-branch reference row that provably
    flows through the whole pipeline:
      * ``final_line[12]`` (Samples) = "n"  -> ``remove_contiguous_samples`` /
        ``merge_contiguous_targets`` classify it origin=ref via SNP(col18)=="n";
        ``remove_n_and_dots`` then rewrites the literal "n" into the literal string
        "NA", so ``resultIntegrator`` sees Samples == "NA" (origin=ref via
        target[13]=="NA") and SNP == "NA" (its target[18]!="NA" guard skips the
        variant parse). Writing "NA" here would break BOTH: pandas reads a written
        "NA" as NaN and re-serializes it as "" (empty) -> resultIntegrator enters the
        variant parse and `"".split("_")[2]` raises IndexError, and Samples ""!="NA"
        misclassifies the row origin=alt.
      * ``final_line[15/16/17]`` (rsID / AF / SNP) = "n"  -> a reference off-target
        has no creating variant; SNP[17] (col18 post-adjust) is the LOAD-BEARING
        marker the merge scripts and resultIntegrator key off.
      * the appended tail ``[ "n", 55, tmp_pos_mms ]`` -- the "n" Reference column and
        the ``55`` sentinel mark this DNA as REFERENCE (its CFD IS the ref CFD, no
        separate ref recompute), identical to the legacy non-IUPAC else-branch row.

    It applies the SAME ``mm_new_t - int(split[8]) > allowed_mms`` budget gate and
    ``pam_ok`` gate as the finalizer, so it keeps the SAME reference rows the legacy
    finalizer would keep (its level-0 seed rows share the reference allele at every
    non-variant column). It does NOT feed ``_record_phase_confirmation`` or
    ``_collect_variant_off_target`` -- the reference is not a variant off-target, so it
    must never populate the phase-confirmation or population-summary companions.
    """
    final_line = split.copy()

    if revert:
        seq_joined = reverse_complement_table("".join(refSeq_prerevert))
    else:
        seq_joined = "".join(refSeq_prerevert)
    target_to_list = list(seq_joined)
    for pos, char in enumerate(realTarget):
        if char == "-":
            target_to_list.insert(pos, "-")

    mm_new_t = 0
    tmp_pos_mms = 0
    for position_t, char_t in enumerate(target_to_list[pos_beg:pos_end]):
        if char_t.upper() != guide_no_pam[position_t]:
            mm_new_t += 1
            tmp_pos_mms = position_t
            if guide_no_pam[position_t] != "-":
                target_to_list[pos_beg + position_t] = char_t.lower()

    pam_ok = True
    for pam_chr_pos, pam_chr in enumerate(target_to_list[pam_begin:pam_end]):
        if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
            pam_ok = False

    # SAME budget + PAM gates as the finalizer (mirrors legacy :828/:830). An over-budget
    # or PAM-failing reference is dropped EXACTLY as the legacy seed rows would be, so we
    # keep the SAME reference rows stable keeps -- no more, no less.
    if mm_new_t - int(split[8]) > allowed_mms:
        return
    if not pam_ok:
        return

    final_line[2] = "".join(target_to_list)
    final_line[7] = str(mm_new_t - int(final_line[8]))
    final_line[9] = str(mm_new_t)
    # a pure reference row never creates the PAM (no alt substituted): leave [10] as-is.
    # The reference markers MUST be the legacy literal "n" (NOT "NA"): every downstream
    # consumer keys origin=ref off the LITERAL "n" in these columns. merge scripts bin a
    # row origin=ref only when SNP(col18 post-adjust) == "n"; remove_n_and_dots then
    # rewrites that literal "n" into the literal string "NA" (via chunk.replace({"n":
    # "NA"})), which resultIntegrator reads as a REAL string -> its `target[18]!="NA"`
    # and `target[13]=="NA"` guards both hold, so no variant-parse and origin=ref.
    # Writing "NA" here instead broke both: pandas reads the written "NA" as NaN and
    # re-serializes it as "" (empty), so resultIntegrator entered the variant parse and
    # `"".split("_")[2]` raised IndexError; and Samples "NA"->"" made target[13]!="NA"
    # so the reference row was misclassified origin=alt.
    final_line[12] = "n"   # Samples  -> origin=ref (literal "n" round-trips to "NA")
    final_line[15] = "n"   # rsID  (reference: no creating variant)
    final_line[16] = "n"   # AF
    final_line[17] = "n"   # SNP / snp_info  (LOAD-BEARING: col18 post-adjust)
    # Reference-row tail: "n" Reference column + sentinel 55 (DNA is REF, CFD is the ref
    # CFD, no separate ref recompute) + tmp_pos count. This is byte-identical to the
    # legacy non-IUPAC else-branch reference row, so every downstream consumer bins it
    # origin=ref.
    final_line.append("n")
    final_line.append(55)
    final_line.append(tmp_pos_mms)
    cluster_to_save.append(final_line)


def _iupac_decomposition_observed(split, guide_no_pam, cluster_to_save):
    """DICT-LESS observed-haplotype branch of iupac_decomposition.

    GATED by the caller on ``mygt is not None`` (a true dict-less install with a
    Tier-1 genotype store). Reads the window's genotypes ONCE, enumerates the DISTINCT
    haplotypes real individuals carry (via ``observed_haplotypes``), and feeds each
    into the shared finalization. Replaces the legacy per-SNP seed split, the 2^k
    growth lattice, the deferred phased peel AND the greedy IUPAC cap -- none of which
    run on this path. Cross-individual chimeras are excluded by construction; the
    IUPAC_CAP / high_variant_density BED are UNNECESSARY here (output is bounded by
    the number of distinct observed sets, not 2^k) and stay live only on the legacy
    path.
    """
    realTarget = split[2]
    replaceTarget = split[2].replace("-", "")
    refSeq = genomeStr[int(split[4]): int(split[4]) + len(replaceTarget)].upper()
    revert = False
    if split[6] == "-":
        revert = True
        replaceTarget = reverse_complement_table(replaceTarget)

    parse_haplotypes = _popsum.parse_haplotypes
    ploidy_of = _t0c.ploidy_of_for_chrom(current_chr)

    # 1) Collect the k ambiguity/variant columns, EXACTLY as the legacy loop gathers
    #    them (same boundary guard, same retrieveFromDict 5-tuple), but WITHOUT the
    #    per-SNP seed split. Each column carries per-alt carrier GT maps (from the
    #    genotype tier via retrieveFromDict's dict-less 5-tuple), the resolved per-alt
    #    alt-index (multiallelic-safe), and per-alt [rsID, AF, snp_info].
    positions = []
    for pos_c, c in enumerate(replaceTarget):
        if c not in iupac_code:
            continue
        if pos_c >= len(refSeq):
            # boundary guard: an IUPAC position past the (truncated) refSeq end would
            # index out of range -- skip it, exactly like the legacy path.
            continue
        snpToReplace, sampleSet, rsID, AF_var, snpInfo = retrieveFromDict(
            pos_c + int(split[4])
        )
        alts, carrier_gts, alt_index, info, ps = [], [], [], [], None
        for i, elem in enumerate(snpToReplace):
            # tokens are "sampleID:gt"; split each on its FIRST ':' (never on '|').
            gt_map = {}
            for tok in sampleSet[i]:
                if ":" in tok:
                    sid, gt = tok.split(":", 1)
                else:
                    sid, gt = tok, ""
                gt_map[sid.strip()] = gt.strip()
            alts.append(elem)
            carrier_gts.append(gt_map)
            # Resolve THIS (pos, alt) record's alt-index token (multiallelic-safe)
            # from its own carrier genotypes. Ploidy 2 is a safe upper bound for the
            # token scan (parse_haplotypes just returns more slots; the classification
            # only cares which non-ref token the record's carriers share).
            alt_index.append(
                _obshap._alt_index_for_record(gt_map.values(), parse_haplotypes, 2)
            )
            info.append([rsID[i], AF_var[i], snpInfo[i]])
        positions.append({
            "pos_c": pos_c,
            "alts": alts,
            "carrier_gts": carrier_gts,
            "alt_index": alt_index,
            "info": info,
            "ps": ps,  # PS not stored on-disk yet -> single-block assumption
        })

    if not positions:
        return

    # 2) Enumerate the distinct observed haplotypes (bounded by ~ploidy*n_carriers).
    #    NOTE: this may be EMPTY (a candidate whose registry variants are absent /
    #    unproductive in every individual). We must NOT early-return on an empty
    #    ``haplotypes`` -- the candidate's REFERENCE off-target is emitted below
    #    REGARDLESS, so its locus is never dropped (LOCUS-COVERAGE FIX; matches the
    #    legacy embedded-refSeq/sentinel-33 representation, which is present on every
    #    IUPAC candidate). The ``enumerate_observed_haplotypes`` contract never yields
    #    the empty/reference haplotype, so the reference row is our responsibility here.
    haplotypes = _obshap.enumerate_observed_haplotypes(
        positions, refSeq, parse_haplotypes, ploidy_of, current_chr,
        max_putative=IUPAC_CAP,
    )

    # 3) Build the reference-with-bulges row ONCE (mirrors the legacy block).
    refSeq_final = reverse_complement_table(refSeq) if revert else refSeq
    refSeq_with_bulges = list(refSeq_final)
    for pos, char in enumerate(realTarget):
        if char == "-":
            refSeq_with_bulges.insert(pos, "-")
    refSeq_with_bulges = "".join(refSeq_with_bulges)

    # 3a) Emit the candidate's REFERENCE off-target ONCE (LOCUS-COVERAGE FIX). This is
    #     emitted whether or not any observed variant haplotype survived, mirroring the
    #     legacy path (which always represents the reference via its alt rows' embedded
    #     refSeq). It is a pure reference alignment (no alt), origin=ref, in-budget/PAM
    #     gated exactly like the finalizer -- so the observed path's LOCUS coverage
    #     matches stable while still collapsing phantom multi-variant combinations.
    _finalize_reference_entry(
        split, realTarget, refSeq_final, refSeq_with_bulges, guide_no_pam,
        revert, cluster_to_save,
    )

    # SAFETY VALVE (visibility, not truncation): a dense/hypervariable window (many
    # ambiguity columns) still emits every REAL observed haplotype -- sensitivity-first,
    # never miss a putative region -- but we LOG it to the shared high-variant-density
    # BED (as the legacy cap does) so an MHC-scale window is surfaced. Unlike the legacy
    # 2^k cap this is O(#distinct observed sets), so no truncation is needed. It unions
    # hap.carriers, so it is a natural no-op when ``haplotypes`` is empty.
    if IUPAC_CAP >= 0 and len(positions) > IUPAC_CAP:
        _samples = set()
        for _h in haplotypes:
            _samples |= _h.carriers
        _start = int(split[4])
        try:
            # 7th field = the FULL IUPAC protospacer (every variant column as its
            # ref+alt ambiguity code) so a user can see all alleles in the window and
            # dig into the other possible alignments (resultIntegrator surfaces this
            # + a "dense region" flag as an integrated_results column).
            hvdr_bed.write(
                "%s\t%d\t%d\t%s\t%d\t%s\t%s\n"
                % (split[3], _start, _start + len(replaceTarget),
                   split[1].replace("-", ""), len(positions),
                   ",".join(sorted(_samples)) if _samples else ".",
                   replaceTarget)
            )
        except Exception:
            pass  # BED logging is best-effort; never break the run
        # SURFACE THE WORST CASE -- never leave a dense region with only a .bed entry
        # and no off-target ROW. The observed haplotypes above are per-individual
        # unions, which in a capped window may all land ABOVE threshold; the region
        # would then appear only in the high-variant-density BED, invisible to a reader
        # of the results table. Emit ONE greedy min-mismatch representative: at each
        # variant column pick the allele that most lowers the mismatch count (additive
        # per column, so this is the exact argmin over all 2^k combinations -- the SAME
        # best-case off-target the legacy dict cap reports at line ~775). This
        # guarantees every dense region ALWAYS surfaces >=1 off-target row, at parity
        # with the dict path. Tagged PUTATIVE (a synthetic worst-case combination, not a
        # confirmed per-individual haplotype); its carriers are the union of the chosen
        # alts' carriers (matching the legacy greedy's sample set). The shared finalizer
        # applies the identical mismatch/PAM budget gate, so this is reported only when
        # the dict path would report it too.
        greedy_seq = list(refSeq)  # pre-revert reference, like the enumerator's seqs
        greedy_samples, greedy_info = set(), []
        for _col in positions:
            _pc = _col["pos_c"]
            _best_mm = _aligned_mm(greedy_seq, realTarget, guide_no_pam, revert)
            _best_i = None
            for _i, _elem in enumerate(_col["alts"]):
                _trial = list(greedy_seq)
                _trial[_pc] = _elem
                _m = _aligned_mm(_trial, realTarget, guide_no_pam, revert)
                # strict improvement, or a tie preferring an alt (keeps PAM-creating
                # variants where the mismatch count is unaffected)
                if _m < _best_mm or (
                    _m == _best_mm
                    and (_best_i is None or _elem < _col["alts"][_best_i])
                ):
                    # deterministic tie-break: on an mm-neutral tie prefer the
                    # lexicographically-smaller alt, so the dict-less greedy rep
                    # matches the legacy one regardless of alt-enumeration order (#139)
                    _best_mm, _best_i = _m, _i
            if _best_i is not None:
                greedy_seq[_pc] = _col["alts"][_best_i]
                greedy_samples |= set(_col["carrier_gts"][_best_i].keys())
                greedy_info.append(_col["info"][_best_i])
        if greedy_info:  # at least one alt improved/held the alignment vs reference
            _finalize_observed_entry(
                split, realTarget, refSeq_final, refSeq_with_bulges,
                guide_no_pam, revert, "".join(greedy_seq),
                sorted(greedy_samples), greedy_info, _obshap.PUTATIVE,
                cluster_to_save,
            )

    # 3b) Finalize each observed VARIANT haplotype through the shared finalizer (may be
    #     an empty loop when no productive observed haplotype exists).
    for hap in haplotypes:
        _finalize_observed_entry(
            split, realTarget, refSeq_final, refSeq_with_bulges, guide_no_pam,
            revert, hap.seq, sorted(hap.carriers), hap.info, hap.phase_state,
            cluster_to_save,
        )


def iupac_decomposition(split, guide_no_bulge, guide_no_pam, cluster_to_save):
    # DICT-LESS observed-haplotype path: GATED on a Tier-1 genotype store being
    # present (mygt) AND the enumerator + its parsing helpers being importable. When
    # ANY is absent (every dict-based install, and any old deploy) fall through to the
    # BYTE-IDENTICAL legacy path below (per-SNP seed split + 2^k lattice + deferred
    # peel + greedy cap). This is the ONLY behavioral fork; nothing below changes.
    if mygt is not None and _obshap is not None and _popsum is not None \
            and _t0c is not None:
        _iupac_decomposition_observed(split, guide_no_pam, cluster_to_save)
        return

    realTarget = split[2]
    replaceTarget = split[2].replace("-", "")
    refSeq = genomeStr[int(split[4]) : int(split[4]) + len(replaceTarget)].upper()

    revert = False
    if split[6] == "-":
        revert = True
        replaceTarget = reverse_complement_table(replaceTarget)

    # dict with IUPAC scompositions
    totalDict = dict()
    totalDict[0] = dict()
    totalDict[0][0] = dict()
    if haplotype_check:
        totalDict[1] = dict()
        totalDict[1][0] = dict()
    countIUPAC = 0
    for pos_c, c in enumerate(replaceTarget):
        if c in iupac_code:
            if pos_c >= len(refSeq):
                # refSeq is genomeStr[pos : pos+len(replaceTarget)]; near a chromosome
                # boundary (or from an out-of-range position on the variant-enriched
                # combined index) it comes back truncated, shorter than replaceTarget.
                # An IUPAC position past its end would then hit listReplaceTarget[pos_c]
                # / refSeq[pos_c] / greedy_seq[pos_c] with an out-of-range index
                # (IndexError, aborting the whole post-analysis). Skip such positions.
                continue
            countIUPAC += 1
            snpToReplace, sampleSet, rsID, AF_var, snpInfo = retrieveFromDict(
                pos_c + int(split[4])
            )
            for i, elem in enumerate(snpToReplace):
                listReplaceTarget = list(refSeq)
                listReplaceTarget[pos_c] = elem
                listInfo = [[rsID[i], AF_var[i], snpInfo[i]]]
                if haplotype_check:
                    haploSamples = {0: [], 1: []}
                    for count, sample in enumerate(sampleSet[i]):
                        sampleInfo = sample.split(":")
                        gt_alleles = sampleInfo[1].split("|")
                        for haplo in totalDict:
                            allele = gt_alleles[haplo] if len(gt_alleles) > haplo else gt_alleles[0]  # haploid guard
                            if allele != "0":
                                haploSamples[haplo].append(sampleInfo[0])
                    totalDict[0][0][(pos_c, elem)] = [
                        listReplaceTarget,
                        set(haploSamples[0]),
                        listInfo,
                    ]
                    totalDict[1][0][(pos_c, elem)] = [
                        listReplaceTarget,
                        set(haploSamples[1]),
                        listInfo,
                    ]
                else:
                    sampleList = list()
                    for count, sample in enumerate(sampleSet[i]):
                        sampleList.append(sample.split(":")[0])
                    totalDict[0][0][(pos_c, elem)] = [
                        listReplaceTarget,
                        set(sampleList),
                        listInfo,
                    ]

    if countIUPAC > 0:  # if found valid alternative targets
        # High-variant-density cap: above IUPAC_CAP ambiguity codes the multi-SNP
        # combination below would enumerate up to 2^countIUPAC sequences. Instead
        # report a single GREEDY representative per haplotype -- at each ambiguity
        # position pick the min-mismatch allele (preferring an alt on ties, which
        # keeps PAM-region variant alleles so PAM-creation is preserved). This is the
        # exact min-mismatch / max-CFD haplotype (mismatch is additive per position,
        # so the greedy equals the argmin over all 2^k combinations) without the
        # blow-up. The region is logged to <out>.high_variant_density_regions.bed.
        capped = IUPAC_CAP >= 0 and countIUPAC > IUPAC_CAP
        if capped:
            _samples = set()
            for _cnt in totalDict:
                for _v in totalDict[_cnt][0].values():
                    _samples |= _v[1]
            _start = int(split[4])
            hvdr_bed.write(
                "%s\t%d\t%d\t%s\t%d\t%s\t%s\n"
                % (
                    split[3],
                    _start,
                    _start + len(replaceTarget),
                    split[1].replace("-", ""),
                    countIUPAC,
                    ",".join(sorted(_samples)) if _samples else ".",
                    replaceTarget,  # full IUPAC protospacer (dig-in aid)
                )
            )
            # Build the greedy representative per haplotype and REPLACE the per-SNP
            # level-0 entries with that single entry, so the finalization below scores
            # exactly one row (bulges/PAM/creation/CFD via the existing code path).
            for count in totalDict:
                # group candidate alt alleles by ambiguity position
                by_pos = {}
                for (pos_c, elem), v in totalDict[count][0].items():
                    by_pos.setdefault(pos_c, []).append((elem, v))
                greedy_seq = list(refSeq)  # pre-revert reference
                greedy_samples, greedy_info = set(), []
                for pos_c, cands in by_pos.items():
                    ref_allele = refSeq[pos_c]
                    ref_mm = _aligned_mm(greedy_seq, realTarget, guide_no_pam, revert)
                    best_elem, best_mm, best_v = ref_allele, ref_mm, None
                    for elem, v in cands:
                        trial = list(greedy_seq)
                        trial[pos_c] = elem
                        m = _aligned_mm(trial, realTarget, guide_no_pam, revert)
                        # strict improvement, or tie preferring an alt (keeps PAM-
                        # creating / present variants where mismatch is unaffected).
                        # Among alts, a mm-neutral tie prefers the lexicographically-
                        # smaller alt (deterministic), so this legacy greedy rep matches
                        # the dict-less one regardless of alt order (#139).
                        if m < best_mm or (
                            m == best_mm and (best_v is None or elem < best_elem)
                        ):
                            best_mm, best_elem, best_v = m, elem, v
                    if best_v is not None:
                        greedy_seq[pos_c] = best_elem
                        greedy_samples |= best_v[1]
                        greedy_info.extend(best_v[2])
                if not greedy_info:  # no allele changed anything: document the region anyway
                    any_v = next(iter(totalDict[count][0].values()), None)
                    if any_v is not None:
                        greedy_info = list(any_v[2])
                        greedy_samples = set(any_v[1])
                totalDict[count][0] = {
                    ("greedy", 0): [greedy_seq, greedy_samples, greedy_info]
                }
        if revert:
            refSeq = reverse_complement_table(refSeq)
        for count in totalDict:
            createdNewLayer = True
            # capped -> range(0): the single greedy entry is already in level 0
            for size in range(0 if capped else countIUPAC):  # else: full enumeration
                if createdNewLayer:
                    createdNewLayer = False
                else:
                    break
                totalDict[count][size + 1] = dict()
                # for each snp in target (fixpoint)
                for key in totalDict[count][size]:
                    # for each other snp in target (> fixpoint)
                    for newkey in totalDict[count][0]:
                        if newkey[-2] > key[-2]:
                            resultSet = totalDict[count][size][key][1].intersection(
                                totalDict[count][0][newkey][1]
                            )  # extract intersection of sample to generate possible multisnp target
                            if len(resultSet) > 0:  # if set is not null
                                createdNewLayer = True
                                # add new snp to preceding target seq with snp
                                replaceTarget1 = totalDict[count][0][newkey][0].copy()
                                replaceTarget2 = totalDict[count][size][key][0].copy()
                                replaceTarget2[newkey[0]] = replaceTarget1[newkey[0]]
                                listInfo2 = totalDict[count][size][key][2].copy()
                                listInfo2.extend(totalDict[count][0][newkey][2])
                                # add to next level the modified seq and set of samples and info of snp
                                combinedKey = key + newkey
                                totalDict[count][size + 1][combinedKey] = [
                                    replaceTarget2,
                                    resultSet,
                                    listInfo2,
                                ]
                                # NOTE (paired with the deferred peel below): the
                                # phased dedup subtraction USED to run HERE, emptying
                                # the parent combo and the level-0 seed as soon as a
                                # combo was formed. But deeper levels grow ONLY by
                                # crossing against the level-0 seeds (line above:
                                # `for newkey in totalDict[count][0]`); emptying those
                                # seeds mid-loop starved levels >= 2 (a sample carrying
                                # k>=4 cis alts got "used up" into disjoint 2-SNP pairs
                                # and its maximal k-combo could never form -> the true
                                # 0-mismatch haplotype was dropped). We keep the growth
                                # loop mutation-free so the FULL lattice forms, then do
                                # the dedup once, after the loop, in the deferred peel.
                                # DO NOT re-add any in-loop subtraction here or the
                                # under-report bug returns.

            # Deferred phased dedup ("peel"): now that the full combination lattice is
            # built (level-0 seeds were never emptied), attribute each sample to the
            # single MAXIMAL cis combo containing all of its co-occurring alt alleles,
            # and remove it from every shorter (ancestor) combo/seed. This reproduces
            # the exact maximal-combo-only attribution the old in-loop subtraction
            # intended, without starving the growth. Phased only: unphased cannot prove
            # cis, so it keeps every sample on every sub-combo (the conservative
            # superset) -- Phase B is skipped there, leaving that path byte-identical.
            # Capped (greedy) path has a single level-0 entry per haplotype (nothing to
            # peel), so it is skipped too.
            if haplotype_check and not capped:
                # position+allele set of every combo key: keys alternate
                # (pos, elem, pos, elem, ...), so pairs (k[0],k[1]),(k[2],k[3]),...
                # An ancestor is a combo whose (pos,elem) set is a STRICT subset -- the
                # (pos,elem) form (not pos-only) guarantees we peel a lower combo only
                # when its alleles match the higher one at the shared positions, so a
                # multiallelic / divergent-allele sibling that merely shares positions
                # is never wrongly peeled.
                posset = {}
                for lvl in totalDict[count]:
                    for k in totalDict[count][lvl]:
                        posset[(lvl, k)] = frozenset(zip(k[0::2], k[1::2]))
                levels_desc = sorted(totalDict[count].keys(), reverse=True)
                for hi in levels_desc:
                    for hk, hv in totalDict[count][hi].items():
                        if not hv[1]:  # no samples on this combo -> nothing to peel with
                            continue
                        hpos = posset[(hi, hk)]
                        for lo in totalDict[count]:
                            if lo >= hi:
                                continue
                            for lk in totalDict[count][lo]:
                                if posset[(lo, lk)] < hpos:  # strict allelic-prefix ancestor
                                    totalDict[count][lo][lk][1] = (
                                        totalDict[count][lo][lk][1] - hv[1]
                                    )

            refSeq_with_bulges = list(refSeq)
            for pos, char in enumerate(realTarget):
                if char == "-":
                    refSeq_with_bulges.insert(pos, "-")

            for position_t, char_t in enumerate(refSeq_with_bulges[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    tmp_pos_mms = position_t
                    if guide_no_pam[position_t] != "-":
                        refSeq_with_bulges[pos_beg + position_t] = char_t.lower()
            # ref sequence with bulges
            refSeq_with_bulges = "".join(refSeq_with_bulges)

            for level in totalDict[count]:
                for key in totalDict[count][level]:
                    # Normally a variant row needs >=1 carrier; on a registry-only
                    # install (#136) there is no carrier source, so emit the site
                    # anyway with a degraded ("NA") Samples column (set below).
                    if len(totalDict[count][level][key][1]) > 0 or registry_only_mode:
                        if revert:
                            totalDict[count][level][key][0] = reverse_complement_table(
                                "".join(totalDict[count][level][key][0])
                            )
                        else:
                            totalDict[count][level][key][0] = "".join(
                                totalDict[count][level][key][0]
                            )

                        final_line = split.copy()

                        target_to_list = list(totalDict[count][level][key][0])
                        for pos, char in enumerate(realTarget):
                            if char == "-":
                                target_to_list.insert(pos, "-")

                        mm_new_t = 0
                        tmp_pos_mms = 0
                        for position_t, char_t in enumerate(
                            target_to_list[pos_beg:pos_end]
                        ):
                            if char_t.upper() != guide_no_pam[position_t]:
                                mm_new_t += 1
                                tmp_pos_mms = position_t
                                if guide_no_pam[position_t] != "-":
                                    target_to_list[pos_beg + position_t] = (
                                        char_t.lower()
                                    )

                        # pam respect input PAM after IUPAC resolution
                        pam_ok = True
                        for pam_chr_pos, pam_chr in enumerate(
                            target_to_list[pam_begin:pam_end]
                        ):
                            if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                                pam_ok = False

                        target_pam_ref = refSeq_with_bulges[pam_begin:pam_end]
                        found_creation = False
                        for pos_pam, pam_char in enumerate(target_pam_ref):
                            # ref char not in set of general pam char
                            if (
                                not iupac_code_set[pam[pos_pam]]
                                & iupac_code_set[pam_char]
                            ):
                                found_creation = True
                        # value of mm and bulges is over allowed threshold, discard target
                        if mm_new_t - int(split[8]) > allowed_mms:
                            continue
                        elif pam_ok:
                            final_line[2] = "".join(target_to_list)
                            final_line[7] = str(mm_new_t - int(final_line[8]))
                            # total differences between targets and guide (mismatches + bulges)
                            final_line[9] = str(mm_new_t)
                            if found_creation:
                                final_line[10] = "".join(
                                    target_to_list[pam_begin:pam_end]
                                )
                            # degraded Samples sentinel on a registry-only install
                            # (empty carrier set) -- "NA" round-trips through the
                            # downstream integrator like the reference-row case,
                            # avoiding the ""->NaN->IndexError hazard (#136).
                            final_line[12] = (
                                ",".join(totalDict[count][level][key][1])
                                if totalDict[count][level][key][1]
                                else "NA"
                            )
                            tmp_matrix = np.array(totalDict[count][level][key][2])
                            if tmp_matrix.shape[0] > 1:
                                final_line[15] = ",".join(tmp_matrix[:, 0])
                                final_line[16] = ",".join(tmp_matrix[:, 1])
                                final_line[17] = ",".join(tmp_matrix[:, 2])
                            else:
                                final_line[15] = str(tmp_matrix[0][0])
                                final_line[16] = str(tmp_matrix[0][1])
                                final_line[17] = str(tmp_matrix[0][2])
                            # report ref DNA seq of target
                            final_line.append(refSeq_with_bulges)
                            # number to activate ref score calculation (active if target is alternative)
                            final_line.append(33)
                            # position of tmp_mms (removed later after processing)
                            final_line.append(tmp_pos_mms)
                            # append processed target to cluster to save
                            cluster_to_save.append(final_line)
                            # ADDITIVE (Phase 3c): record this VARIANT off-target's
                            # identity + creating-variant (SNP) columns for the
                            # companion population-summary TSV. GATED on ``myreg``
                            # (built later at module scope) so it is a no-op / empty
                            # collection on any legacy install; deduped by identity so
                            # per-haplotype / per-score duplicates collapse to one row.
                            # This reads ONLY final_line columns and never mutates it,
                            # so the bestMerge output is byte-identical.
                            _collect_variant_off_target(final_line)


def preprocess_CFD_score(target):
    # preprocess target then calculate CFD score
    if do_scores:
        if target[0] == "DNA":
            cfd_score = calc_cfd(
                target[1][int(target[bulge_pos]) :],
                target[2].upper()[int(target[bulge_pos]) : -3],
                target[2].upper()[-2:],
                mm_scores,
                pam_scores,
                do_scores,
            )
            # append to target the CFD score of the aligned sequence (alt or ref)
            target.append("{:.3f}".format(cfd_score))
            # -3 position is a placeholder for ref score
            if (
                target[-3] == 55
            ):  # if 55 sequence is ref so no score have to be calculated
                target[-3] = "{:.3f}".format(cfd_score)
            if (
                target[-3] == 33
            ):  # if 33 sequence is alt so ref score must be calculated
                cfd_ref_score = calc_cfd(
                    target[1][int(target[bulge_pos]) :],
                    target[-4].upper()[int(target[bulge_pos]) : -3],
                    target[-4].upper()[-2:],
                    mm_scores,
                    pam_scores,
                    do_scores,
                )
                target[-3] = "{:.3f}".format(cfd_ref_score)
        else:
            cfd_score = calc_cfd(
                target[1],
                target[2].upper()[:-3],
                target[2].upper()[-2:],
                mm_scores,
                pam_scores,
                do_scores,
            )
            target.append("{:.3f}".format(cfd_score))
            if target[-3] == 55:
                target[-3] = "{:.3f}".format(cfd_score)
            if target[-3] == 33:
                cfd_ref_score = calc_cfd(
                    target[1],
                    target[-4].upper()[:-3],
                    target[-4].upper()[-2:],
                    mm_scores,
                    pam_scores,
                    do_scores,
                )
                target[-3] = "{:.3f}".format(cfd_ref_score)
    else:
        # no score calculated, append -1 in CFD score and in position -3 insert -1 value (-1 means no score calculated)
        cfd_score = -1
        target.append("{:.3f}".format(cfd_score))
        target[-3] = "{:.3f}".format(cfd_score)
    # print(target)
    # print('cfd', cfd_score)
    return target


def preprocess_CRISTA_score(cluster_targets):
    # list with scored targets
    cluster_scored = list()
    index_to_null = list()

    # skip scoring for CRISTA, remove to activate scoring
    # do_scores = False

    if do_scores:
        pass
    else:
        for target in cluster_targets:
            target_CRISTA = target.copy()
            crista_score = -1  # null score
            target_CRISTA[-2] = "{:.3f}".format(crista_score)
            target_CRISTA.append("{:.3f}".format(crista_score))
            cluster_scored.append(target_CRISTA)
        return cluster_scored

    # preprocess target then calculate CRISTA score
    sgRNA_non_aligned_list = list()
    DNA_aligned_list = list()
    DNAseq_from_genome_list = list()
    # process all found targets
    for index, target in enumerate(cluster_targets):
        # list with non-aligned sgRNA
        sgRNA_non_aligned_list.append(str(target[1])[: len(str(target[1])) - 3] + "NGG")
        # list with aligned DNA
        DNA_aligned_list.append(str(target[2]))
        # first 5 nucleotide to add to protospacer
        pre_protospacer_DNA = genomeStr[int(target[4]) - 5 : int(target[4])].upper()
        # protospacer taken directly from the aligned target
        protospacerDNA = str(target[2]).replace("-", "")
        if target[6] == "-":
            protospacerDNA = reverse_complement_table(protospacerDNA)
        # last 5 nucleotides to add to protospacer
        post_protospacer_DNA = genomeStr[
            int(target[4]) + len(target[1]) : int(target[4]) + len(target[1]) + 5
        ].upper()

        # DNA seq extracted from genome and append to aligned DNA seq from CRISPRme
        complete_DNA_seq = (
            str(pre_protospacer_DNA) + protospacerDNA + str(post_protospacer_DNA)
        )

        for elem in iupac_nucleotides:
            if elem in complete_DNA_seq:
                complete_DNA_seq = complete_DNA_seq.replace(elem, "")

        # trim the 3' and 5' end to avoid sequences longer than 29
        len_DNA_seq = len(complete_DNA_seq)
        first_half = complete_DNA_seq[int(len_DNA_seq / 2) - 14 : int(len_DNA_seq / 2)]
        second_half = complete_DNA_seq[int(len_DNA_seq / 2) : int(len_DNA_seq / 2) + 15]
        complete_DNA_seq = first_half + second_half
        if target[6] == "-":
            complete_DNA_seq = reverse_complement_table(complete_DNA_seq)

        # if 'N' is present in the reference DNA seq, we must use a fake DNA seq to complete the aligned
        # that will be discarded after
        if (
            # A CRISTA window that isn't a full 29 nt of A/C/G/T cannot be scored:
            # near a chromosome boundary, or on the variant-enriched combined index a
            # window can be stripped down to empty by the IUPAC removal above, which
            # then hit `len(full_dna_seq)==0` -> ZeroDivisionError and aborted the whole
            # post-analysis. Null those targets (CRISTA score -1) like the N case.
            len(complete_DNA_seq) != 29
            or "N" in complete_DNA_seq
            or "n" in complete_DNA_seq
            or "N" in DNA_aligned_list[-1]
            or "n" in DNA_aligned_list[-1]
        ):
            complete_DNA_seq = "A" * 29
            DNA_aligned_list[-1] = "A" * len(str(target[2]))
            index_to_null.append(index)

        # append sequence to DNA list
        DNAseq_from_genome_list.append(complete_DNA_seq)

    # calculate scores for alt sequence
    crista_score_list_alt = list()
    if do_scores:
        crista_score_list_alt = CRISTA_predict_list(
            sgRNA_non_aligned_list, DNA_aligned_list, DNAseq_from_genome_list
        )

    # preprocess target then calculate CRISTA score
    sgRNA_non_aligned_list = list()
    DNA_aligned_list = list()
    DNAseq_from_genome_list = list()
    # process all ref sequences in targets
    for index, target in enumerate(cluster_targets):
        # list with non-aligned sgRNA
        sgRNA_non_aligned_list.append(str(target[1])[: len(str(target[1])) - 3] + "NGG")
        # list with aligned DNA
        if "n" not in target[-3]:
            DNA_aligned_list.append(str(target[-3]))
        else:
            DNA_aligned_list.append(str(target[2]))
        # first 5 nucleotide to add to protospacer
        pre_protospacer_DNA = genomeStr[int(target[4]) - 5 : int(target[4])]
        # protospacer taken directly from the ref genome
        protospacerDNA = genomeStr[int(target[4]) : int(target[4]) + len(target[1])]
        # last 5 nucleotides to add to protospacer
        post_protospacer_DNA = genomeStr[
            int(target[4]) + len(target[1]) : int(target[4]) + len(target[1]) + 5
        ]

        # DNA seq extracted from genome and append to aligned DNA seq from CRISPRme
        complete_DNA_seq = (
            str(pre_protospacer_DNA) + protospacerDNA + str(post_protospacer_DNA)
        )

        for elem in iupac_nucleotides:
            if elem in complete_DNA_seq:
                complete_DNA_seq = complete_DNA_seq.replace(elem, "")

        # trim the 3' and 5' end to avoid sequences longer than 29
        len_DNA_seq = len(complete_DNA_seq)
        first_half = complete_DNA_seq[int(len_DNA_seq / 2) - 14 : int(len_DNA_seq / 2)]
        second_half = complete_DNA_seq[int(len_DNA_seq / 2) : int(len_DNA_seq / 2) + 15]
        complete_DNA_seq = first_half + second_half
        if target[6] == "-":
            complete_DNA_seq = reverse_complement_table(complete_DNA_seq)

        # if 'N' is present in the reference DNA seq, we must use a fake DNA seq to complete the aligned
        # that will be discarded after
        if (
            # A CRISTA window that isn't a full 29 nt of A/C/G/T cannot be scored:
            # near a chromosome boundary, or on the variant-enriched combined index a
            # window can be stripped down to empty by the IUPAC removal above, which
            # then hit `len(full_dna_seq)==0` -> ZeroDivisionError and aborted the whole
            # post-analysis. Null those targets (CRISTA score -1) like the N case.
            len(complete_DNA_seq) != 29
            or "N" in complete_DNA_seq
            or "n" in complete_DNA_seq
            or "N" in DNA_aligned_list[-1]
            or "n" in DNA_aligned_list[-1]
        ):
            complete_DNA_seq = "A" * 29
            DNA_aligned_list[-1] = "A" * len(str(target[2]))
            index_to_null.append(index)

        # append sequence to DNA list
        DNAseq_from_genome_list.append(complete_DNA_seq)

    # calculate score
    crista_score_list_ref = list()
    if do_scores:
        crista_score_list_ref = CRISTA_predict_list(
            sgRNA_non_aligned_list, DNA_aligned_list, DNAseq_from_genome_list
        )

    for index, target in enumerate(cluster_targets):
        target_CRISTA = target.copy()
        # if any of the scored target is not valid, due to Ns in the sequence, return a -1 score
        if index in index_to_null:
            crista_score = -1  # null score
            target_CRISTA[-2] = "{:.3f}".format(crista_score)
            target_CRISTA.append("{:.3f}".format(crista_score))
        else:
            # else report the correct score
            if target_CRISTA[-2] == 55:  # reference target have duplicate score
                target_CRISTA[-2] = "{:.3f}".format(crista_score_list_alt[index])
                target_CRISTA.append("{:.3f}".format(crista_score_list_alt[index]))
            if target_CRISTA[-2] == 33:  # alternative target scoring
                target_CRISTA[-2] = "{:.3f}".format(crista_score_list_ref[index])
                target_CRISTA.append("{:.3f}".format(crista_score_list_alt[index]))
        # append to final score cluster
        cluster_scored.append(target_CRISTA)

    return cluster_scored


def calculate_scores(cluster_to_save):
    # function to calculate score for each input target
    # input is target line splitted in list format
    # list of functions to calculate specific score (to add a score, simply add your function to this call and update the clusters list in return)
    cluster_with_CFD_score = list()
    cluster_with_CRISTA_score = list()

    for target in cluster_to_save:  # calculate CFD score for each target
        target_CFD = target.copy()
        cluster_with_CFD_score.append(preprocess_CFD_score(target_CFD))

    # process score for each target in cluster, at the same time to improve execution time
    cluster_with_CRISTA_score = preprocess_CRISTA_score(cluster_to_save)

    # REMOVED TO CHECK IF FILE IS RETURN WITH IDENTICAL ROWS COUNT

    # analyze CFD scored targets, returning for each guide,chr,cluster_pos the highest scoring target (or multiple in case of equal)
    # df_CFD = pd.DataFrame(cluster_with_CFD_score, columns=['Bulge_type', 'crRNA', 'DNA', 'Chromosome',
    #                                                        'Position', 'Cluster_Position', 'Direction', 'Mismatches',
    #                                                        'Bulge_Size', 'Total', 'PAM_gen', 'Var_uniq', 'Samples', 'Annotation_Type',
    #                                                        'Real_Guide', 'rsID', 'AF', 'SNP', 'Reference_target', 'CFD',
    #                                                        'Seq_in_cluster', 'CFD_ref'])
    # # group by over real_guide,chr,cluster_pos to avoid mixing targets
    # # select lowest count of mm+bul
    # idx_fewest_mm_bul = df_CFD.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'Total'].transform(min) == df_CFD['Total']
    # df_CFD_fewest = df_CFD[idx_fewest_mm_bul]
    # df_CFD_fewest.drop_duplicates(
    #     ['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples', 'Mismatches', 'Bulge_Size'], inplace=True)
    # # select highest score
    # idx_max_score = df_CFD.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'CFD'].transform(max) == df_CFD['CFD']
    # df_CFD_best_score = df_CFD[idx_max_score]
    # # remove duplicate rows (possible due to haplotypes and variants)
    # df_CFD_best_score.drop_duplicates(['Real_Guide', 'Chromosome',
    #                                    'Cluster_Position', 'SNP', 'Samples', 'CFD'], inplace=True)
    # frames = [df_CFD_fewest, df_CFD_best_score]
    # df_CFD = pd.concat(frames)
    # cluster_with_CFD_score = df_CFD.values.tolist()

    # df_CRISTA = pd.DataFrame(cluster_with_CRISTA_score, columns=['Bulge_type', 'crRNA', 'DNA', 'Chromosome',
    #                                                              'Position', 'Cluster_Position', 'Direction', 'Mismatches',
    #                                                              'Bulge_Size', 'Total', 'PAM_gen', 'Var_uniq', 'Samples', 'Annotation_Type',
    #                                                              'Real_Guide', 'rsID', 'AF', 'SNP', 'Reference_target', 'CFD',
    #                                                              'Seq_in_cluster', 'CFD_ref'])
    # # select lowest count of mm+bul
    # idx_fewest_mm_bul = df_CRISTA.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'Total'].transform(min) == df_CRISTA['Total']
    # df_CRISTA_fewest = df_CRISTA[idx_fewest_mm_bul]
    # df_CRISTA_fewest.drop_duplicates(
    #     ['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples', 'Mismatches', 'Bulge_Size'], inplace=True)
    # # select highest score
    # idx_max_score = df_CRISTA.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'CFD'].transform(max) == df_CRISTA['CFD']
    # df_CRISTA_best_score = df_CRISTA[idx_max_score]
    # # remove duplicate rows (possible due to haplotypes and variants)
    # df_CRISTA_best_score.drop_duplicates(['Real_Guide', 'Chromosome',
    #                                       'Cluster_Position', 'SNP', 'Samples', 'CFD'], inplace=True)
    # frames = [df_CRISTA_fewest, df_CRISTA_best_score]
    # df_CRISTA = pd.concat(frames)
    # cluster_with_CRISTA_score = df_CRISTA.values.tolist()

    return [cluster_with_CFD_score, cluster_with_CRISTA_score]


def _collect_needed_dict_keys(cluster_path, chrom):
    """Return the SNP-dictionary keys ('<chrom>,<1-based pos>') this chromosome's
    targets will query. iupac_decomposition() calls retrieveFromDict() for each
    genomic position spanned by a target that contains an IUPAC code, so we take the
    full span of every such target -- a SUPERSET of the exact query set, so no
    queried position is ever missed. This lets us load only these entries instead of
    the whole per-chromosome dictionary (a 1000G+HGDP chr2 dict is ~13 GB on disk,
    ~26 GB via json.load; loading only the queried positions keeps peak RAM small)."""
    needed = set()
    try:
        with open(cluster_path) as fh:
            fh.readline()  # skip header
            for line in fh:
                split = line.rstrip("\n").split("\t")
                if len(split) <= 4:
                    continue
                target = split[2]
                if not any((c in iupac_nucleotides) for c in target):
                    continue
                try:
                    start = int(split[4])
                except (ValueError, IndexError):
                    continue
                span = len(target.replace("-", ""))
                for p in range(start, start + span):
                    needed.add(chrom + "," + str(p + 1))
    except OSError:
        pass
    return needed


def _load_dict_targeted(dict_path, needed_keys):
    """Load ONLY `needed_keys` from the per-chromosome SNP dictionary, streaming the
    JSON (ijson) so peak RAM is proportional to the queried off-target positions, not
    the whole file. Also derives the dataset phasing flag from the streamed genotype
    separators ('|' phased vs '/' unphased), exactly as the previous whole-dict scan
    did. Returns (mydict, haplotype_check). Falls back to a filtered json.load when
    ijson is unavailable (correct, but reads the whole file into RAM once).

    Prefers a gzip-compressed dict (``my_dict_<chrom>.json.gz``) when present, so a
    batteries install keeps the per-sample dicts compressed on disk (~40-50GB instead
    of ~152GB) and reads them on the fly; falls back to a plain ``.json`` for older /
    uncompressed installs. ijson and json.load both stream through a gzip file object
    transparently."""
    # resolve .gz-vs-plain (prefer the compressed file when it exists)
    if os.path.exists(dict_path + ".gz"):
        dict_path = dict_path + ".gz"
    _is_gz = dict_path.endswith(".gz")
    result = {}
    haplo = False
    decided = False
    try:
        import ijson
        with (gzip.open(dict_path, "rb") if _is_gz else open(dict_path, "rb")) as fh:
            for key, value in ijson.kvitems(fh, ""):
                if not decided and isinstance(value, str):
                    if "|" in value:
                        haplo, decided = True, True
                    elif "/" in value:
                        decided = True
                if key in needed_keys:
                    result[key] = value
        return result, haplo
    except ImportError:
        with (gzip.open(dict_path, "rt") if _is_gz else open(dict_path)) as fh:
            full = json.load(fh)
        for key, value in full.items():
            if not decided and isinstance(value, str):
                if "|" in value:
                    haplo, decided = True, True
                elif "/" in value:
                    decided = True
            if key in needed_keys:
                result[key] = value
        del full
        return result, haplo


def _resolve_registry_paths(dict_path, chrom):
    """Resolve a Tier-0 registry (bin, idx) for ``chrom`` from the SNP dict path,
    or return None if no registry is present.

    Layout mirrors the dicts: the search pipeline passes ``dict_path`` as
    ``<...>/Dictionaries/dictionaries_<vcf>/my_dict_<chrom>.json`` (see
    post_analisi_snp.sh). The registry lives alongside it as
    ``<...>/Dictionaries/registry_<vcf>/reg_<chrom>.bin`` + ``.idx``. We derive it
    by swapping the ``dictionaries_`` folder prefix for ``registry_`` and the
    ``my_dict_<chrom>.json[.gz]`` file for ``reg_<chrom>.bin``/``.idx``. Detection
    is by FILE EXISTENCE only -- absence => None => legacy path unchanged."""
    dict_dir = os.path.dirname(dict_path)
    parent = os.path.dirname(dict_dir)
    dict_folder_name = os.path.basename(dict_dir)
    if dict_folder_name.startswith("dictionaries_"):
        reg_folder_name = "registry_" + dict_folder_name[len("dictionaries_"):]
    else:
        # non-standard layout: try a sibling "registry_<same suffix>" then a
        # co-located registry dir; either way we only USE it if the files exist.
        reg_folder_name = "registry_" + dict_folder_name
    reg_dir = os.path.join(parent, reg_folder_name)
    bin_path = os.path.join(reg_dir, "reg_" + str(chrom) + ".bin")
    idx_path = os.path.join(reg_dir, "reg_" + str(chrom) + ".idx")
    if os.path.exists(bin_path) and os.path.exists(idx_path):
        return bin_path, idx_path
    return None


def _resolve_genotype_paths(dict_path, chrom):
    """Resolve a Tier-1 genotype store (bin, idx) for ``chrom`` from the SNP dict
    path, or return None if no genotype tier is present.

    MIRRORS ``_resolve_registry_paths`` exactly, only swapping the folder prefix
    (``dictionaries_`` -> ``genotypes_``) and the per-chromosome file naming
    (``my_dict_<chrom>.json[.gz]`` -> ``gt_<chrom>.bin`` + ``gt_<chrom>.idx``):
    the store lives alongside the dicts + registry as
    ``<...>/Dictionaries/genotypes_<vcf>/gt_<chrom>.bin`` + ``.idx``. Detection is
    by FILE EXISTENCE only -- absence => None => no genotype tier (Samples come
    from the dict, unchanged)."""
    dict_dir = os.path.dirname(dict_path)
    parent = os.path.dirname(dict_dir)
    dict_folder_name = os.path.basename(dict_dir)
    if dict_folder_name.startswith("dictionaries_"):
        gt_folder_name = "genotypes_" + dict_folder_name[len("dictionaries_"):]
    else:
        # non-standard layout: try a sibling "genotypes_<same suffix>"; either way
        # we only USE it if the files exist.
        gt_folder_name = "genotypes_" + dict_folder_name
    gt_dir = os.path.join(parent, gt_folder_name)
    bin_path = os.path.join(gt_dir, "gt_" + str(chrom) + ".bin")
    idx_path = os.path.join(gt_dir, "gt_" + str(chrom) + ".idx")
    if os.path.exists(bin_path) and os.path.exists(idx_path):
        return bin_path, idx_path
    return None


def _write_population_summary_companion():
    """ADDITIVE companion population-summary write. FULLY GUARDED + GATED on ``myreg``.

    Writes ``<outputFile>.population_summary.tsv`` -- a SEPARATE file next to the
    bestMerge / output stem, one row per VARIANT off-target keyed by the SAME identity
    columns the bestMerge uses (Chromosome, Position, Direction, crRNA, DNA) plus the
    creating SNP column, so it can be joined back. It NEVER touches the existing
    bestMerge / altMerge / integrated_results columns or the Samples column.

    Gate: when ``myreg`` is None (legacy / dict-only install, and any old deploy
    without the population-summary modules) this writes NOTHING (no file created) --
    BYTE-IDENTICAL legacy behavior. Any error is caught + logged + skipped here, so
    the companion write can NEVER break the run.

    The single dataset-wide ``phased`` flag is resolved by the companion writer
    (explicit -> gt manifest -> detect from GT strings -> conservative default), and
    we feed it the phasing already sniffed from the dict scan (``haplotype_check``) so
    a dict-backed run agrees with the dict decomposition's cis/trans handling.
    """
    if myreg is None:
        return  # GATE: no registry -> byte-identical legacy behavior (no file)
    if _popsum is None or _popsum_companion is None or _t0c is None:
        return  # modules absent (old deploy) -> no companion summary
    try:
        # The sample axis + ploidy model for THIS chromosome. Prefer the genotype
        # tier's axis (needed for exact multi-variant combinations); fall back to the
        # registry's axis if the reader exposes one, else None (single-variant rows
        # still resolve straight from the registry).
        axis = None
        if mygt is not None and hasattr(mygt, "_axis"):
            axis = mygt._axis
        elif hasattr(myreg, "axis"):
            try:
                axis = myreg.axis()
            except Exception:
                axis = None
        ploidy_of = _t0c.ploidy_of_for_chrom(current_chr)
        out_path = outputFile + ".population_summary.tsv"

        def _on_row_error(ot, err):
            print("population-summary companion: skipped one off-target -", err)

        wrote = _popsum_companion.write_companion(
            out_path,
            _variant_off_targets,
            myreg,
            mygt,
            axis,
            ploidy_of,
            _popsum,
            panel_cls=getattr(_popsum, "Panel", None),
            # dataset-wide phasing already sniffed from the dict scan; the writer
            # still resolves manifest/detect when this is False-by-default and a gt
            # tier is present.
            phased=(True if haplotype_check else None),
            global_group_id=(t0_reg.GLOBAL_GROUP_ID if t0_reg is not None else "global"),
            observed_gt_strings=None,
            on_error=_on_row_error,
        )
        if wrote:
            print(
                "Wrote population-summary companion (%d variant off-target row[s]) to %s"
                % (len(_variant_off_targets), out_path)
            )
    except Exception as _ps_err:  # ADDITIVE + guarded: never break the run
        print("population-summary companion skipped for", current_chr, "-", _ps_err)


def _write_phase_confirmation_companion():
    """ADDITIVE phase-confirmation companion write. FULLY GUARDED + GATED on the
    dict-less branch (``mygt``).

    Writes ``<outputFile>.phase_confirmation.tsv`` -- a SEPARATE joinable file next to
    the output stem, one row per dict-less variant off-target: identity columns + the
    CONFIRMED/PUTATIVE flag. It NEVER touches the bestMerge / Samples columns (the flag
    is out of band) so the byte-identical legacy path is unaffected.

    Gate: when ``mygt`` is None (legacy / dict-only install, or an old deploy without
    the enumerator) NOTHING is written -- byte-identical legacy behavior. Any error is
    caught + logged + skipped so this can NEVER break the run."""
    if mygt is None or _phase_companion is None:
        return  # GATE: no genotype tier -> no dict-less enumeration -> no companion
    if not _phase_confirmation_rows:
        return  # nothing enumerated (no variant off-targets on this chromosome)
    try:
        out_path = outputFile + ".phase_confirmation.tsv"
        n = _phase_companion.write_companion(out_path, _phase_confirmation_rows)
        print(
            "Wrote phase-confirmation companion (%d variant off-target row[s]) to %s"
            % (n, out_path)
        )
    except Exception as _pc_err:  # ADDITIVE + guarded: never break the run
        print("phase-confirmation companion skipped for", current_chr, "-", _pc_err)


# INPUT AND SETTINGS
# fasta of the reference chromosome
inFasta = open(sys.argv[1], "r")
current_chr = inFasta.readline().strip().replace(">", "")  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = "".join(genomeStr).upper()
# string of the whole chromosome on single line
genomeStr = genomeStr.replace("\n", "")
# targets clusterized by chr and ordered by position
inTarget = open(sys.argv[3], "r")
# text file with PAM sequence and length
inPAMfile = open(sys.argv[4], "r")
# outfile path
outputFile = sys.argv[5]
# max allowed mismatches in search (to validate ref targets in alternative case)
allowed_mms = int(sys.argv[6])
# column of bulges count
bulge_pos = 8
# header to insert into final file
header = "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tReference"
# cfd graphs pre-processing (deprecated)
cfd_for_graph = {"ref": [0] * 101, "var": [0] * 101}

# OUT BEST FILES FOR EACH SCORING SYSTEM

# file with best CFD targets
cfd_best = open(outputFile + ".bestCFD.txt", "w+")
cfd_best.write(header + "\tCFD\n")  # Write header

# file with best mm+bul targets
mmblg_best = open(outputFile + ".bestmmblg.txt", "w+")
mmblg_best.write(header + "\tCFD\n")  # Write header

# file with best CRISTA targets
crista_best = open(outputFile + ".bestCRISTA.txt", "w+")
crista_best.write(header + "\tCFD\n")  # Write header

# High-variant-density cap: a protospacer window overlapping many ambiguity codes
# (dense/hypervariable regions, e.g. MHC, and unbounded for unphased VCFs) would
# otherwise enumerate up to 2^k haplotype sequences. Above CRISPRME_IUPAC_CAP codes
# we report only the single-variant targets (k rows, not 2^k) and log the region,
# so a real off-target risk is still surfaced without the combinatorial blow-up.
IUPAC_CAP = int(os.environ.get("CRISPRME_IUPAC_CAP", "10"))
hvdr_bed = open(outputFile + ".high_variant_density_regions.bed", "w")
hvdr_bed.write(
    "#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\tiupac_protospacer\n"
)

# Load ONLY the SNP-dictionary entries this chromosome's targets actually query,
# streaming the (up to ~26 GB) per-chromosome dict so peak RAM stays proportional to
# the off-targets, not the dataset. Previously this json.load()ed the entire
# chromosome dictionary, which OOM-killed the run on genome-wide 1000G+HGDP searches
# (and cascaded into the cryptic "Killed ... EmptyDataError" downstream).
haplotype_check = False
mydict = {}
dict_tier_present = False  # True once a per-sample dict FILE loads (mode 1); stays
#                            False on a dict-less install (modes 2/3) -- the reliable
#                            "dict absent" signal for the registry-only guard (#136).
try:
    _needed_keys = _collect_needed_dict_keys(sys.argv[3], current_chr)
    mydict, haplotype_check = _load_dict_targeted(sys.argv[2], _needed_keys)
    dict_tier_present = True
    print(
        f"Loaded {len(mydict)} SNP dictionary entr(y/ies) for {current_chr} "
        f"(only the positions its targets span); haplotype processing {haplotype_check}"
    )
except Exception as _dict_err:  # keep going with no SNP annotation, never crash here
    print("No dict found (or unreadable) for", current_chr, "-", _dict_err)

# ADDITIVE Tier-0 registry detection (dictless redesign). If a registry exists for
# this chromosome (alongside the dict), open a module-level RegistryReader that
# retrieveFromDict() consults for corrected AF / rsID / metadata. When absent (the
# common case, and any old deploy without tier0_registry) ``myreg`` stays None and
# retrieveFromDict falls through to the byte-identical legacy dict path.
myreg = None
if t0_reg is not None:
    try:
        _reg_paths = _resolve_registry_paths(sys.argv[2], current_chr)
        if _reg_paths is not None:
            myreg = t0_reg.RegistryReader(_reg_paths[0], _reg_paths[1])
            print(
                f"Opened Tier-0 registry for {current_chr} "
                f"({len(myreg)} records) at {_reg_paths[0]}"
            )
    except Exception as _reg_err:  # never break the legacy path on a bad/old registry
        myreg = None
        print("No Tier-0 registry (or unreadable) for", current_chr, "-", _reg_err)

# ADDITIVE Tier-1 genotype-tier detection (dictless redesign, Phase 2). If a
# genotype store exists for this chromosome (alongside the dict / registry), open a
# module-level GenotypeReader that retrieveFromDict() consults to reconstruct the
# Samples column WITHOUT the per-sample dict (true dictless: Tier-0 registry +
# genotype tier, NO dict). When absent (the common case, and any old deploy without
# tier1_genotypes) ``mygt`` stays None and the Samples column comes from the dict
# (legacy/augment path) unchanged. This open is GUARDED exactly like the registry:
# import failure / missing files / open error -> mygt=None, never break the legacy
# or registry paths.
mygt = None
if t1_gt is not None:
    try:
        _gt_paths = _resolve_genotype_paths(sys.argv[2], current_chr)
        if _gt_paths is not None:
            mygt = t1_gt.GenotypeReader(_gt_paths[0], _gt_paths[1])
            print(
                f"Opened Tier-1 genotype store for {current_chr} "
                f"({len(mygt)} records) at {_gt_paths[0]}"
            )
    except Exception as _gt_err:  # never break legacy/registry on a bad/old store
        mygt = None
        print("No Tier-1 genotype store (or unreadable) for", current_chr, "-", _gt_err)

# Registry-only install (#136): Tier-0 registry present, but NEITHER a per-sample
# dict NOR a genotype tier (e.g. `download --no-genotypes`). Here every variant
# carrier list is empty, so the variant finalizer's `len(samples) > 0` guard would
# DROP the off-target site entirely. Detect the mode once so that finalizer can
# instead emit a DEGRADED row (Samples = "NA") -- the site + creating-variant
# rsID/AF are still surfaced, just without per-sample resolution. This is False on
# every other install (legacy: myreg None; registry+dict: dict_tier_present;
# dictless-with-genotypes: mygt set), so those paths stay byte-identical.
registry_only_mode = (myreg is not None) and (mygt is None) and (not dict_tier_present)
if registry_only_mode:
    print(
        f"Registry-only install for {current_chr}: emitting variant off-targets "
        f"with degraded (NA) Samples -- no genotype tier to resolve carriers."
    )

# check PAM position and relative coordinates on targets
pam_at_beginning = False
line = inPAMfile.read().strip()
pam = line.split(" ")[0]
len_pam = int(line.split(" ")[1])
guide_len = len(pam) - len_pam
pos_beg = 0
pos_end = None
pam_begin = 0
pam_end = len_pam * (-1)
if len_pam < 0:
    guide_len = len(pam) + len_pam
    pam = pam[: (len_pam * (-1))]
    len_pam = len_pam * (-1)
    pos_beg = len_pam
    pos_end = None
    pam_begin = 0
    pam_end = len_pam
    pam_at_beginning = True
else:
    pam = pam[(len_pam * (-1)) :]
    pos_beg = 0
    pos_end = len_pam * (-1)
    pam_begin = len_pam * (-1)
    pam_end = None

# start time counter
global_start = time.time()

# open mm and pam scores matrices for CFD
mm_scores, pam_scores = get_mm_pam_scores()

# if conditions, execute score (guidelen==20,pamlen==3,pam_at_beginning==FALSE)
do_scores = True
if len_pam != 3 or guide_len != 20 or pam_at_beginning:
    # sys.stderr.write('CFD SCORE IS NOT CALCULATED WITH GUIDES LENGTH != 20 OR PAM LENGTH !=3 OR UPSTREAM PAM')
    do_scores = False

# START TARGET PROCESSING

# keep track of current analyzed cluster (necessary to check if the cluster is terminated)
# current_guide_chr_pos_direction = ''
# count_cluster_dimension = 0

# skip header
inTarget.readline()

# list with clusterized targets in list format (contains ref seq and all other alternative targets)
cluster_to_save = list()
# read lines from target file
for line in inTarget:
    # split target into list
    split = line.strip().split("\t")
    # sgRNA sequence (with bulges and PAM)
    guide = split[1]
    # found target on DNA (with bulges, mismatches and PAM)
    target = split[2]
    guide_no_bulge = split[1].replace("-", "")
    guide_no_pam = guide[pos_beg:pos_end]

    # check if targets cointains IUPAC nucleotide
    if any((c in iupac_nucleotides) for c in target):
        iupac_decomposition(split, guide_no_bulge, guide_no_pam, cluster_to_save)
    else:
        # process_iupac = False
        # append to respect file format for post analysis
        # null ref sequence
        split.append("n")
        # specific value to represent a ref target to avoid recount score
        split.append(55)
        # count of mm_bul for ref sequence in case of alternative target
        split.append(0)
        cluster_to_save.append(split)

    if len(cluster_to_save) >= 100000:
        # after reading 100k lines from file and creating the cluster, start processing it
        clusters_with_scores = calculate_scores(cluster_to_save)

        for count, cluster in enumerate(clusters_with_scores):
            for target in cluster:
                if count == 0:  # CFD target
                    # remove count of tmp_mms
                    target.pop(-2)
                    # save CFD targets
                    cfd_best.write("\t".join(target) + "\t" + str(0) + "\n")
                    # save mm-bul targets
                    mmblg_best.write("\t".join(target) + "\t" + str(0) + "\n")
                if count == 1:  # CRISTA target
                    # remove count of tmp_mms
                    target.pop(-2)
                    # save CRISTA targets
                    crista_best.write("\t".join(target) + "\t" + str(0) + "\n")
        cluster_to_save = list()


if len(cluster_to_save):
    pass
else:
    # if cluster to save is empty, skip processing
    # close all files
    cfd_best.close()
    mmblg_best.close()
    crista_best.close()
    hvdr_bed.close()
    # rewrite header file
    os.system(
        "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
        + outputFile
        + ".bestCFD.txt"
    )
    os.system(
        "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
        + outputFile
        + ".bestmmblg.txt"
    )
    os.system(
        "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
        + outputFile
        + ".bestCRISTA.txt"
    )
    # cfd dataframe write
    cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
    cfd_dataframe.to_csv(outputFile + ".CFDGraph.txt", sep="\t", index=False)
    # ADDITIVE + guarded + gated-on-myreg: companion population-summary TSV.
    _write_population_summary_companion()
    # ADDITIVE + guarded + gated-on-mygt: dict-less phase-confirmation companion TSV.
    _write_phase_confirmation_companion()
    # print complete and exit with no error
    print("ANALYSIS COMPLETE IN", time.time() - global_start)
    exit(0)

# process cluster of targets if less then 1mln rows total
clusters_with_scores = calculate_scores(cluster_to_save)

for count, cluster in enumerate(clusters_with_scores):
    for target in cluster:
        # print(target)
        if count == 0:  # CFD target
            # remove count of tmp_mms
            target.pop(-2)
            # save CFD targets
            cfd_best.write("\t".join(target) + "\t" + str(0) + "\n")
            # save mm-bul targets
            mmblg_best.write("\t".join(target) + "\t" + str(0) + "\n")
        if count == 1:  # CRISTA target
            # remove count of tmp_mms
            target.pop(-2)
            # save CRISTA targets
            crista_best.write("\t".join(target) + "\t" + str(0) + "\n")

cfd_best.close()
mmblg_best.close()
crista_best.close()
hvdr_bed.close()

os.system(
    "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
    + outputFile
    + ".bestCFD.txt"
)
os.system(
    "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
    + outputFile
    + ".bestmmblg.txt"
)
os.system(
    "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
    + outputFile
    + ".bestCRISTA.txt"
)


cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
cfd_dataframe.to_csv(outputFile + ".CFDGraph.txt", sep="\t", index=False)

# ADDITIVE + guarded + gated-on-myreg: companion population-summary TSV.
_write_population_summary_companion()
# ADDITIVE + guarded + gated-on-mygt: dict-less phase-confirmation companion TSV.
_write_phase_confirmation_companion()

print("ANALYSIS COMPLETE IN", time.time() - global_start)
