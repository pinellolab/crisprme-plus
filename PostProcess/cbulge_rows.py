#!/usr/bin/env python
"""CRISPR-Bulge scores of one batch of finalized off-target rows (2.10.0).

The SNP and the indel post-analysis (new_simple_analysis.py, analisi_indels_NNN.py) score
their rows in batches: every 100,000 rows of a chromosome, and the rest at its end. The
batch is one request to the CRISPR-Bulge worker, and the model's floats depend on which
rows share it (the request is cut into mini-batches of the model's own size, and a row's
float can change in its last bits with the rows beside it), so a row's score can differ in
its last printed digit from one batching to another (measured by the 2.10.0 cost judges:
9 of 10,227 per-sample and 80 of 890,497 population-level chr22 rows, relative difference
up to 1.2e-5, when one chromosome's rows were scored in other batches). That is why this
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
        sg = str(target[1])[: len(str(target[1])) - 3] + "NGG"
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
    # below the old 3-decimal floor (5e-4), so rounding would make it unreachable and
    # silently drop ~7% of validated off-targets from the sensitive tier.
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
