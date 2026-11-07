#!/usr/bin/env python

import sys
import time

import indel_snp_cis as _isc


def _rep_tiebreak(x):
    """Deterministic, intrinsic final tie-break for the per-cluster representative
    (issue #139).

    Appended as the LAST element of every sort key in ``get_best_targets``. It is only
    ever consulted when two rows tie EXACTLY on the pre-existing key elements (CFD,
    mm/bulges, and the cluster-position ``int(x[pos])`` which is constant within a
    cluster). In that exact-tie case the pre-fix winner fell to Python's stable-sort
    INPUT order -- i.e. the alt/enumeration order -- which differs between the dict
    (VCF order) and dict-less (sorted) paths, so they picked DIFFERENT representative
    alignments for the same cluster/carriers/CFD.

    This key is intrinsic to the row's own fields (aligned target DNA at col 2 -- unique
    per alignment within a cluster, but for the 2.8.1 note below -- plus the true genomic
    Position at col ``pos - 1`` as a secondary safety discriminator), never input order.
    Because it is consulted only when all earlier key elements are exactly equal, it
    CANNOT reorder any cluster whose best row is already unique: output stays byte-identical
    for those (~99.7%); for the exactly-tied clusters it makes both paths converge on the
    same representative.

    2.7.10: among rows that tie EXACTLY on CFD/severity/position, prefer the one carried
    by the MOST individuals ("highest CFD, ties -> most carriers"). This only ever
    reorders exactly-tied clusters (still a no-op for the unique-best ~99.7%), and makes
    the per-cluster headline the most prevalent of the equally-severe real haplotypes;
    the per-site ``observed_haplotypes`` companion lists every one regardless. The carrier
    count is read from the samples column (``true_guide - 2``); ref rows and empty/sentinel
    rosters count as 0, so this is inert for them and they fall through to the DNA/pos key.

    2.8.1: aligned DNA is not quite unique within a cluster. A row of SNPs and a row with
    an indel can give the same bases at the same position, and their aligned DNA is then
    the same string (rarely on the dictionary path: on chr21 and chr22 of the default
    index with three guides, no such tie). They fell through to
    input order, which is the whole-line sort of the merge step, so the REF alignment's
    letter case and then the sample IDs decided which row was reported. The last key
    element settles it: the row of SNPs comes first. A variant token is
    ``<contig>_<pos>_<ref>_<alt>``, split from the right because a contig name can hold
    '_'; a REF or ALT longer than one base comes from the indel path. It is consulted only
    after everything else ties, so no other choice of representative can change. It holds
    among rows with carriers only: the rows of a sites-only index (mega) have none and
    their REF alignment was never upper case, so their ties stay in the sort's order.

    2.8.2: two alignments of one indel target string can share their aligned DNA too, since
    the merge keys an indel row by its aligned spacer as well. Those with one carrier set are
    folded into one entry (_fold_same_carriers), and an exact tie between them keeps the
    input order, in which the merge kept the first of them before 2.8.2.
    """
    try:
        _samples = str(x[true_guide - 2])
        _ncarriers = (
            0 if _samples in ("", ".", "NA", "n", "NO_SAMPLES")
            else len([_c for _c in _samples.split(",") if _c and _c not in (".", "NA", "NO_SAMPLES")])
        )
    except Exception:  # never let the tie-break break the sort
        _ncarriers = 0
    _has_indel = _ncarriers > 0 and any(
        len(_var) == 4 and (len(_var[2]) > 1 or len(_var[3]) > 1)
        for _var in (_tok.rsplit("_", 3) for _tok in str(x[snp_info]).split(","))
    )
    # -_ncarriers => most carriers first; then the intrinsic deterministic discriminators;
    # False sorts first => a row of SNPs before a row with an indel.
    return (-_ncarriers, str(x[2]), str(x[pos - 1]), _has_indel)


def _merged_list(kept, other):
    """One list column (samples, rsID, AF) of a variant row that two datasets give. A
    sites-only dataset's row has no carriers and can have no allele frequency, which
    adjust_cols.py leaves empty: it adds no empty item to the other dataset's list."""
    return kept + "," + other if kept and other else kept or other


def _holds_indel_token(variants):
    """True iff a row's variant list (<contig>_<pos>_<ref>_<alt>,...) names a variant of
    the indel path: a REF or ALT longer than one base. That is an indel, or a multi-base
    record of equal length (an MNV, AC>TG), which CRISPRitz add-variants also puts in a
    fake-indel window, so its rows are resolved per carrier like indel rows. A token is
    split from the right because a contig name can hold '_'."""
    for tok in variants.split(","):
        p = tok.rsplit("_", 3)
        if len(p) == 4 and (len(p[2]) > 1 or len(p[3]) > 1):
            return True
    return False


def _merge_row(kept, ele, snp_info):
    """Add a merged row's samples, rsIDs and AFs to the row that stays (_merged_list: a
    sites-only dataset's empty lists add no empty item)."""
    kept[true_guide - 2] = _merged_list(
        kept[true_guide - 2], ele[true_guide - 2]
    )  # true_guide - 2 points to samples column
    kept[snp_info - 2] = _merged_list(
        kept[snp_info - 2], ele[snp_info - 2]
    )  # snp_info - 2 points to rsID column
    kept[snp_info - 1] = _merged_list(
        kept[snp_info - 1], ele[snp_info - 1]
    )  # snp_info - 1 points to AF column


def _first_base_offset(start, strand, variants, target, ref):
    """How far past its start an indel row's first aligned base is on the haplotype
    carrying its indel (indel_snp_cis.row_first_base), from the row's Position, strand,
    variants, aligned target and REF text; 0 for a row with no indel or one this cannot
    read."""
    try:
        toks = []
        for tok in variants.split(","):
            p = tok.rsplit("_", 3)
            toks.append((int(p[1]), p[2], p[3]))
        indel = next((t for t in toks if len(t[1]) != len(t[2])), None)
        if indel is None:
            return 0
        start = int(start)
        return _isc.row_first_base(strand, start, target, ref, indel,
                                   [t for t in toks if t is not indel]) - start
    except (IndexError, ValueError):
        return 0


def _haplotype_cluster_position(row, memo=None):
    """An indel row's cluster position on the haplotype carrying its indel: its
    Cluster_Position plus how far past its start its first aligned base is there.

    That is the Cluster_Position itself, except on a row whose first base is one an
    insertion adds. Since 2.8.2 such a row starts at the insertion's anchor
    (indel_snp_cis.reference_start), while its first base is the k-th added base, k bases
    past the anchor on the haplotype, which the row's bases tell
    (indel_snp_cis.row_first_base). So two alignments that start on different bases of
    one insertion have one Cluster_Position, while on the haplotype they are as many bases
    apart as 2.8.1 placed them; and two alignments before one PAM with different bulges,
    one place on the haplotype (on the + strand with a 3' PAM a DNA bulge starts one base
    earlier and the Cluster_Position adds it back), can have two Cluster_Positions when
    they start inside an insertion. The row's Position is the column before its
    Cluster_Position, its strand the one after, its aligned target and REF text columns
    2 and 3. A row with no indel, or one this cannot read: its Cluster_Position.

    ``memo`` (a dict, one per cluster) keeps the offset by what it is read from: the
    Position, strand and variants, and the bases of the target and REF text, gaps and case
    aside, as row_first_base reads them. The aligned spacer is not among them, and since
    2.8.2 the merge keys an indel row by it as well, so every alignment of a target comes
    here; without the memo each read the haplotype again, the whole ALT for a row at the
    anchor of a long insertion."""
    cluster = int(row[pos])
    try:
        what = (row[pos - 1], row[pos + 1], row[snp_info],
                row[2].replace("-", "").upper(), row[3].replace("-", "").upper())
    except IndexError:
        return cluster
    if memo is None:
        return cluster + _first_base_offset(*what)
    if what not in memo:
        memo[what] = _first_base_offset(*what)
    return cluster + memo[what]


def _fold_same_carriers(rows, fold_of, snp_info):
    """Alignments of one site naming the same variants with the same carriers are one
    entry, as for every other row: the first of them in ``rows``, which is sorted by the
    file's own criterion, so the entry is the best by it. It takes the rsIDs and AFs of
    the others. ``fold_of`` maps a row (by id) to its (cluster position on the haplotype
    carrying its indel, variants, carriers); other rows stay as they are. The carrier
    sets do not depend on the sort, so the three files keep one row count per cluster
    and their alternative alignments stay aligned line by line.

    The cluster position is the one on the haplotype (_haplotype_cluster_position): one
    site's alignments there, the bulged ones folding into the best, as up to 2.8.1. Since
    2.8.2 alignments that start on different bases of one insertion share its anchor and
    one site; on the Cluster_Position alone they would fold into one entry, and all but
    the best would leave both tables. They are other targets of the haplotype (2.8.1
    reported them at their own starts), and stay: the site's row and its alternatives."""
    kept = {}
    out = []
    for row in rows:
        ident = fold_of.get(id(row))
        if ident is None:
            out.append(row)
        elif ident in kept:
            _merge_row(kept[ident], row, snp_info)
        else:
            kept[ident] = row
            out.append(row)
    for row in kept.values():
        for col in (true_guide - 2, snp_info - 2, snp_info - 1):
            row[col] = ",".join(sorted(set(row[col].split(","))))
    return out


def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    # avoid crush when cluster is empty in the first call
    if not cluster:
        return
    list_ref = []
    seen_ref = set()
    dict_var = {}
    for ele in cluster:
        if ele[snp_info] == "n":
            # list (and count in #Seq_in_cluster) each reference alignment once. The
            # same alignment comes in more than once: every variant dataset's
            # post-analysis carries all the rows of the one reference search, and on
            # the dict-less tiers new_simple_analysis gives the reference alignment of
            # each variant window on the plus strand again (_finalize_reference_entry).
            # Copies are equal in every column but the count and the two scores: each
            # dataset scores the rows itself, and CRISPR-Bulge can give one alignment
            # two values. Keeping such a copy in that list alone would shift it against
            # the other two, which the pipeline pastes by position (pr -m). The merge
            # step's sort puts the higher score first, and the first copy is the one
            # kept. A row that differs in any other column is another alignment and is
            # kept.
            ref_key = tuple(ele[: cfd - 1] + ele[cfd + 2 :])
            if ref_key in seen_ref:
                continue
            seen_ref.add(ref_key)
            list_ref.append(ele)
        else:
            # merge samples of identical targets (coming from different VCF datasets).
            # Since 2.8.2 a row that starts on a base an insertion adds starts at its
            # anchor, so two copies of one target inside one insertion (a repeat it
            # holds), the same bases with the same carriers, are such identical rows too:
            # one entry, at the anchor (2.8.1 reported each at its own start).
            key = (ele[pos], ele[snp_info])
            # An indel row lists the carriers of its own sequence (each carrier haplotype
            # resolved on its own, analisi_indels_NNN.py), and two alignments at one site
            # can name the same variants with different carriers: an RNA-bulge alignment
            # does not reach the base a SNP sits on, so it holds carriers of that SNP
            # that the full alignment does not. So only rows with one sequence merge here.
            # And only rows with one alignment of it (2.8.2): a carrier resolved PUTATIVE
            # (phase unknown) gets the strongest target its own alleles allow on each
            # alignment (indel_snp_cis.worst_case_bases): the guide's base where the
            # alignment pairs a SNP's base with a guide base, the CFD's choice (the
            # reference base in a search that scores no CFD) where it leaves that base as a
            # DNA bulge. So two alignments of one target string can list it with different
            # bases, and with the target alone as the key it went onto the other alignment
            # with the first one's bases, a target weaker than its own there (the SNP+indel
            # oracle's S19). Rows of one alignment still merge (one target from two VCF
            # datasets); the alignments of a target with one carrier set are folded into
            # the best below (_fold_same_carriers), carrier sets compared whole.
            if _holds_indel_token(ele[snp_info]):
                key += (ele[2].upper(), ele[1].upper())
            if key in dict_var.keys():
                _merge_row(dict_var[key][0], ele, snp_info)
            else:
                dict_var[key] = [ele]
    # Alignments of one site naming the same variants with the same carriers are folded
    # into one entry once the variant rows are sorted (_fold_same_carriers), the site
    # being their cluster position on the haplotype carrying the indel, read once per
    # target, REF text, place and variants (the memo), not once per aligned spacer.
    fold_of = {}
    first_base_memo = {}
    for key in dict_var.keys():
        if len(key) == 4:  # an indel or MNV row (its target and aligned spacer in the key)
            row = dict_var[key][0]
            fold_of[id(row)] = (_haplotype_cluster_position(row, first_base_memo), key[1],
                                frozenset(row[true_guide - 2].split(",")))

    final_list_best_ref = list()
    var_only = False
    for target in list_ref:
        final_list_best_ref.append(target)
    if not final_list_best_ref:
        var_only = True

    final_list_best_var = list()
    # for each snp_info in dict, extract the targets
    for key in dict_var.keys():
        list_var = dict_var[key]
        # copy the targets in the variant list, adding unique if no ref target is found
        for target in list_var:
            if var_only:
                target[12] = "y"
            final_list_best_var.append(target)

    temp_final_list_best_var = list()
    # for target in final_list_best_var:
    for target in final_list_best_var:
        # remove duplicates; sorted() keeps the order deterministic across
        # processes (bare set() iteration is PYTHONHASHSEED-dependent).
        # remove duplicates into snp info col
        target[snp_info] = ",".join(sorted(set(target[snp_info].split(","))))
        # remove duplicate into rsID col
        target[snp_info - 2] = ",".join(sorted(set(target[snp_info - 2].split(","))))
        # remove duplicate into AF col
        target[snp_info - 1] = ",".join(sorted(set(target[snp_info - 1].split(","))))
        # remove duplicate into samples col
        target[true_guide - 2] = ",".join(sorted(set(target[true_guide - 2].split(","))))
        # append to temp list
        temp_final_list_best_var.append(target)

    # final list with polished targets (no duplicates in snp data)
    final_list_best_var = temp_final_list_best_var

    # check if lists are empty
    validity_check_ref = False
    validity_check_var = False
    if final_list_best_ref:
        validity_check_ref = True
    if final_list_best_var:
        validity_check_var = True

    # extract best target for each criteria
    if sort_order == "score":
        # sort per score (CFD or CRISPR-Bulge)
        if validity_check_ref:
            if sort_criteria == "mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (-float(x[cfd]), int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (-float(x[cfd]), int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (-float(x[cfd]), int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm+bulges,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,mm,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm+bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,mm+bulges,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,bulges,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )

        if validity_check_var:
            if sort_criteria == "mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (-float(x[cfd]), int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (-float(x[cfd]), int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (-float(x[cfd]), int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm+bulges,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,mm,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm+bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,mm+bulges,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,bulges,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        -float(x[cfd]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            final_list_best_var = _fold_same_carriers(final_list_best_var, fold_of, snp_info)

        if var_only:  # no ref found
            # count the residual targets in the list
            final_list_best_var[0][cfd - 1] = str(len(final_list_best_var) - 1)
            # append the best target to best_file
            fileOut.write("\t".join(final_list_best_var[0]))
            # pop the best target from the list
            bestTarget = final_list_best_var.pop(0)
        elif validity_check_ref and validity_check_var:  # ref and var targets found
            if float(final_list_best_ref[0][cfd]) >= float(final_list_best_var[0][cfd]):
                final_list_best_ref[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                fileOut.write("\t".join(final_list_best_ref[0]))
                bestTarget = final_list_best_ref.pop(0)
            else:
                final_list_best_var[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                fileOut.write("\t".join(final_list_best_var[0]))
                bestTarget = final_list_best_var.pop(0)
        else:  # only ref
            final_list_best_ref[0][cfd - 1] = str(len(final_list_best_ref) - 1)
            fileOut.write("\t".join(final_list_best_ref[0]))
            bestTarget = final_list_best_ref.pop(0)
        # write all the remaining targets in the alt file
        for count, elem in enumerate(final_list_best_ref):
            final_list_best_ref[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            fileOut_disc.write(("\t".join(elem)))
        for count, elem in enumerate(final_list_best_var):
            final_list_best_var[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            fileOut_disc.write(("\t".join(elem)))
    else:
        # sort for total (mm+bul) in target
        if validity_check_ref:
            if sort_criteria == "mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref, key=lambda x: (int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x))
                )
            elif sort_criteria == "bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref, key=lambda x: (int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x))
                )
            elif sort_criteria == "mm":
                final_list_best_ref = sorted(
                    final_list_best_ref, key=lambda x: (int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x))
                )
            elif sort_criteria == "mm+bulges,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (int(x[total]), int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm+bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (int(x[total]), int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "bulges,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (int(x[total - 1]), int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (int(x[total - 1]), int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (int(x[total - 2]), int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (int(x[total - 2]), int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm+bulges,bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,mm,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm+bulges,mm":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,mm+bulges,bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,bulges,mm+bulges":
                final_list_best_ref = sorted(
                    final_list_best_ref,
                    key=lambda x: (
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )

        if validity_check_var:
            if sort_criteria == "mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var, key=lambda x: (int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x))
                )
            elif sort_criteria == "bulges":
                final_list_best_var = sorted(
                    final_list_best_var, key=lambda x: (int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x))
                )
            elif sort_criteria == "mm":
                final_list_best_var = sorted(
                    final_list_best_var, key=lambda x: (int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x))
                )
            elif sort_criteria == "mm+bulges,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (int(x[total]), int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm+bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (int(x[total]), int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "bulges,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (int(x[total - 1]), int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (int(x[total - 1]), int(x[total - 2]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (int(x[total - 2]), int(x[total]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (int(x[total - 2]), int(x[total - 1]), int(x[pos]),
                        *_rep_tiebreak(x)),
                )
            elif sort_criteria == "mm+bulges,bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm+bulges,mm,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm+bulges,mm":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[total - 2]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "bulges,mm,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        int(x[total - 1]),
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,mm+bulges,bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        int(x[total - 2]),
                        int(x[total]),
                        int(x[total - 1]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            elif sort_criteria == "mm,bulges,mm+bulges":
                final_list_best_var = sorted(
                    final_list_best_var,
                    key=lambda x: (
                        int(x[total - 2]),
                        int(x[total - 1]),
                        int(x[total]),
                        int(x[pos]),
                        *_rep_tiebreak(x),
                    ),
                )
            final_list_best_var = _fold_same_carriers(final_list_best_var, fold_of, snp_info)

        if var_only:  # no ref found
            # count the residual targets in the list
            final_list_best_var[0][cfd - 1] = str(len(final_list_best_var) - 1)
            # append the best target to best_file
            fileOut.write("\t".join(final_list_best_var[0]))
            # pop the best target from the list
            bestTarget = final_list_best_var.pop(0)
        elif validity_check_ref and validity_check_var:  # ref and var targets found
            if int(final_list_best_ref[0][total]) <= int(final_list_best_var[0][total]):
                final_list_best_ref[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                fileOut.write("\t".join(final_list_best_ref[0]))
                bestTarget = final_list_best_ref.pop(0)
            else:
                final_list_best_var[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                fileOut.write("\t".join(final_list_best_var[0]))
                bestTarget = final_list_best_var.pop(0)
        else:  # only ref
            final_list_best_ref[0][cfd - 1] = str(len(final_list_best_ref) - 1)
            fileOut.write("\t".join(final_list_best_ref[0]))
            bestTarget = final_list_best_ref.pop(0)
        # write all the remaining targets in the alt file
        for count, elem in enumerate(final_list_best_ref):
            final_list_best_ref[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            fileOut_disc.write(("\t".join(elem)))
        for count, elem in enumerate(final_list_best_var):
            final_list_best_var[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            fileOut_disc.write(("\t".join(elem)))


tau = int(sys.argv[3])  # range in bp to merge targets
chrom = int(sys.argv[4]) - 1  # chromoso
pos = int(sys.argv[5]) - 1  # position of target
total = int(sys.argv[6]) - 1  # mm+bul value
true_guide = int(sys.argv[7]) - 1  # real guide used in the search
snp_info = int(sys.argv[8]) - 1  # snp_info (ref_alt_allele)
cfd = int(sys.argv[9]) - 1  # CFD score
sort_order = str(sys.argv[10])

# NEW
sort_criteria = str(sys.argv[11])

# -1 is to get the correct "python enumeration" from the bash script

start = time.time()
with open(sys.argv[1], "r") as fileIn:
    header = fileIn.readline()
    with open(sys.argv[2], "w") as fileOut:
        with open(sys.argv[2] + ".discarded_samples", "w") as fileOut_disc:
            fileOut.write(header)
            fileOut_disc.write(header)
            prev_pos = -(tau + 1)
            best_row = ""
            prev_guide = ""
            prev_chr = ""
            prev_snp = ""
            cluster = []
            for line in fileIn:
                splitted = line.split("\t")
                if (
                    prev_guide != splitted[true_guide]
                    or prev_chr != splitted[chrom]
                    or int(splitted[pos]) - prev_pos > tau
                ):
                    get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)
                    cluster = [splitted]
                else:
                    cluster.append(splitted)
                prev_guide = splitted[true_guide]
                prev_pos = int(splitted[pos])
                prev_chr = splitted[chrom]
                prev_snp = splitted[snp_info]

            get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)

print("Merging done in: " + str(time.time() - start))
