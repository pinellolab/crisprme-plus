#!/usr/bin/env python
"""CRISPR-Bulge scores of one batch of finalized off-target rows (2.10.0).

The SNP and the indel post-analysis (new_simple_analysis.py, analisi_indels_NNN.py) score
their rows in batches: every 100,000 rows of a chromosome, and the rest at its end. The
batch is one request to the CRISPR-Bulge worker, and the model's floats depend on which
rows share it (the request is cut into mini-batches of the model's own size, and a row's
float can change in its last bits with the rows beside it), so a row's score can differ in
its last printed digit from one batching to another (measured by the 2.10.0 cost judges:
9 of 10,227 per-sample and 80 of 890,496 population-level chr22 rows, relative difference
up to 1.2e-5, when one chromosome's rows were scored in other batches; the judges' file has
890,497 lines with its header). That is why this
lives in one place: a post-analysis
cut into pieces (pool_post_analisi_snp.py, pool_post_analisi_indel.py) leaves the
CRISPR-Bulge scoring to the merge of its pieces (merge_pieces.py), which forms exactly the
batches one process forms (BATCH_ROWS, cut after the input line that reaches it) and
scores them with this same function, so the scores are those of one process, byte for byte.

``score_rows`` is the body both scripts' ``preprocess_CRISPR_BULGE_score`` had, unchanged
(their copies were identical); they now call it.
"""

# One process scores its rows every time this many have gathered, after the input line
# that reached it (the loops of new_simple_analysis.py and analisi_indels_NNN.py).
BATCH_ROWS = 100000

# The guide PAM the model is given. The ensemble scores the spacer; the guide's own PAM
# is normalized to NGG so a guide written with an IUPAC/N PAM (NRR, NNN, NRG...) reaches
# the model in the shape it was trained on.
PAM_MASK = "NGG"


def mask_guide_pam(crrna_aligned, pam_len=3, pam=PAM_MASK):
    """Normalize the guide's PAM to ``pam`` WITHOUT disturbing the alignment.

    The naive form of this -- ``crrna[:-pam_len] + pam`` -- is wrong on exactly one
    row shape, and wrong twice over: an in-place DNA-bulge alignment is ONE COLUMN
    LONGER than the guide, so when the bulge's gap column falls inside the PAM the
    guide's 3 PAM bases occupy 4 columns (``...NG-G``). Slicing a fixed 3 off the end
    then (a) DELETES the gap, so the caller's ``total_gaps`` test sees 0 gaps on a row
    that has one and scores a row it would otherwise have had to judge, and (b) hands
    the model ``...NNGG`` -- a gapless 24-mer with a phantom base that is not an
    alignment of anything. The score that comes back is not this site's.

    Overwriting the last ``pam_len`` NON-GAP columns instead leaves every gap column
    where the aligner put it, so the model sees the alignment it was given:

        ACGT...ACA NG-G   ->  ACGT...ACA NG-G   (PAM already NGG; bulge preserved)
        ACGT...ACA NN-N   ->  ACGT...ACA NG-G   (NRR/NNN guide PAM normalized)
        ACGT...ACA NNN    ->  ACGT...ACA NGG    (no gap: identical to the old form)

    With no gap in the PAM columns this is byte-identical to the slice form, which is
    every row of a CRISPRme search (the search anchors the PAM exactly, so a bulge
    never lands in it) -- verified on 52.2 M aligned spacer+PAM strings from the
    finished 2.10.0 runs, where the two forms differ on zero rows. They differ on
    externally supplied lists, which is where PAM-internal bulges do occur: 30,131 of
    the 2,940,588 rows (1.0%) of the CRISPR-Bulge external benchmark.

    ``pam_len`` is the caller's PAM width. When it does not match the mask (a
    non-SpCas9 ``--pam-length``) the old slice is kept verbatim: it changes the string's
    length, so the caller's ``len(off) != len(sg)`` rule nulls the row -- which is what
    a guide outside the 20+3 shape this ensemble scores should do, and is not something
    to start quietly scoring here.
    """
    if pam_len != len(pam):
        return crrna_aligned[: len(crrna_aligned) - pam_len] + pam
    chars = list(crrna_aligned)
    k = pam_len - 1
    for i in range(len(chars) - 1, -1, -1):
        if k < 0:
            break
        if chars[i] == "-":
            continue
        chars[i] = pam[k]
        k -= 1
    return "".join(chars)


def score_rows(cluster_targets, do_scores):
    """CRISPR-Bulge scoring (analogue of the retired legacy-scorer preprocessing path).

    Same output contract (target[-2] marker -> score, plus an appended score), same
    alt/ref two-pass + null semantics, so the write loop + report are unchanged. The
    ONLY differences vs CRISPR-Bulge: (1) the model needs just the aligned (sgRNA, off-target)
    pair -- NO 29-nt genomic window -- and (2) scoring runs in the dedicated `cbulge`
    conda env via the scorer-runner (persistent worker). The runner degrades gracefully
    to -1.0 if the env/model is absent, so a run never breaks.

    CRISPR-Bulge's encoder requires per-row len(sg)==len(off); for rows we cannot score
    (N in the aligned DNA, or a length mismatch) we substitute an equal-length dummy so
    the batch stays valid, then null them by index afterwards (mirrors the CRISPR-Bulge path,
    which substitutes 'A'*29).

    ``do_scores``: False where the guide or the PAM is not the 20-nt, 3' 3-nt shape the
    scores are defined for (each script decides it from its PAM file): every row then
    gets -1.
    """
    import scorer_runner  # lazy: only when CRISPR-Bulge is selected

    cluster_scored = list()
    index_to_null = list()

    if not do_scores:
        for target in cluster_targets:
            t = target.copy()
            t[-2] = "{:.6g}".format(-1)
            t.append("{:.6g}".format(-1))
            cluster_scored.append(t)
        return cluster_scored

    # This CRISPR-Bulge ensemble faithfully scores a protospacer+PAM alignment with AT MOST
    # ONE 1-bp bulge (model input seq_len=24). We deliberately score ONLY what the model can
    # represent and NULL (-1) everything else, so a value in the column is always trustworthy:
    #   * N in the aligned off-target, or sg/off length mismatch  -> null
    #   * >1 gap column (>=2 bulges, OR a single >=2-bp bulge)     -> null
    #   * over-length (> 24)                                        -> null
    # Nulled rows fall back to CFD (the PRIMARY score) + mismatches/bulges, which still flag
    # them. We deliberately do NOT collapse/reduce a multi-bulge alignment to a single bulge:
    # an adversarial review showed the naive reduction corrupts the alignment (deletes real
    # bases; mis-counts a 2-bp bulge as two events) and its "conservative" error direction is
    # unprovable for a non-monotonic GRU -- an unacceptable silent under-flag risk in a clinical
    # tool. Multi-bulge off-targets are the rare, heavily-disrupted (low-activity) tail; CFD +
    # edit distance carry them. A scoreable equal-length dummy is substituted for nulled rows so
    # one bad row can't fail the whole batch. CFD (the primary score) is unaffected.
    _MAXLEN = 24
    _SG_DUMMY = "A" * 20 + "NGG"
    _OFF_DUMMY = "A" * 20 + "AGG"

    def _pair(target, off_field, index):
        # gap-preserving PAM normalization: a DNA bulge inside the PAM must stay a
        # bulge here, or the row below is judged and scored as a gapless 24-mer that
        # is not this site (see mask_guide_pam).
        sg = mask_guide_pam(str(target[1]))
        off = str(off_field)
        total_gaps = sg.count("-") + off.count("-")
        if (("N" in off) or ("n" in off) or (len(off) != len(sg))
                or total_gaps > 1 or len(sg) > _MAXLEN):
            index_to_null.append(index)
            return _SG_DUMMY, _OFF_DUMMY
        return sg, off

    # ALT pass: aligned off-target = target[2]
    pairs_alt = [_pair(t, t[2], i) for i, t in enumerate(cluster_targets)]
    scores_alt = scorer_runner.CRISPR_BULGE_predict_list(
        [p[0] for p in pairs_alt], [p[1] for p in pairs_alt])

    # REF pass: aligned off-target = Reference_target (target[-3]) unless it carries an 'n'
    pairs_ref = []
    for i, target in enumerate(cluster_targets):
        ref = str(target[-3]) if "n" not in str(target[-3]) else str(target[2])
        pairs_ref.append(_pair(target, ref, i))
    scores_ref = scorer_runner.CRISPR_BULGE_predict_list(
        [p[0] for p in pairs_ref], [p[1] for p in pairs_ref])

    # CRISPR-Bulge scores are written at FULL precision ("{:.6g}"), NOT 3 decimals:
    # the sensitive operating point (cb_operating_point, 95%-recall default 1.9e-4) is
    # below the old 3-decimal floor (5e-4), so rounding would make it unreachable: on the
    # held-out benchmark a 5e-4 floor leaves the sensitive tier holding 267 of the 287
    # validated off-targets instead of 273 -- 6 sites, 2.1 points of recall, lost to the
    # rounding alone (docs/offtarget_risk_summary_METHODS.md section 3.5, which is the one
    # place that split is stated).
    # The rows to null as a set (2.10.0): a list's membership test read the whole list for
    # each row, and in a dense locus nearly every row has two bulges (a DNA and an RNA one)
    # and is nulled, so a batch of 100,000 rows took about a minute here alone.
    nulled = set(index_to_null)
    for index, target in enumerate(cluster_targets):
        t = target.copy()
        if index in nulled:
            s = -1
            t[-2] = "{:.6g}".format(s)
            t.append("{:.6g}".format(s))
        else:
            # mirror the CRISPR-Bulge marker convention exactly (55 = ref row, 33 = alt row)
            if t[-2] == 55:
                t[-2] = "{:.6g}".format(scores_alt[index])
                t.append("{:.6g}".format(scores_alt[index]))
            if t[-2] == 33:
                t[-2] = "{:.6g}".format(scores_ref[index])
                t.append("{:.6g}".format(scores_alt[index]))
        cluster_scored.append(t)

    return cluster_scored
