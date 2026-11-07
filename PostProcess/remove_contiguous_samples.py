#!/usr/bin/env python

import sys
import time


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


def _triples(row, snp_info):
    """(variant, rsID, AF) of each variant of a row as the post-analysis writes it, one
    to one, in its order: AF None where the row has no AF (a sites-only dataset's, which
    adjust_cols.py leaves empty). None where the lists are not one to one (no row the
    post-analysis writes)."""
    toks = row[snp_info].split(",")
    n = len(toks)
    rsids = row[snp_info - 2].split(",") if row[snp_info - 2] else ["."] * n
    afs = row[snp_info - 1].split(",") if row[snp_info - 1] else [None] * n
    if len(rsids) != n or len(afs) != n:
        return None
    return list(zip(toks, rsids, afs))


def _write_variant_lists(row, snp_info, triples):
    """A variant row's variant, rsID and AF columns from its (variant, rsID, AF)
    triples, those of the rows merged into it included: each variant once, sorted as text
    as before, and the rsID and AF lists made of the distinct triples, by variant, so
    the n-th rsID and AF are those of the n-th variant wherever the row holds one rsID
    and AF a variant (two datasets that give one variant two frequencies list both, one
    after the other). Sorting each list on its own, as text, and keeping each value
    once put an rsID or AF in another variant's place, and two variants of one frequency
    showed it once. A dataset with no AF (sites-only) adds no AF item where the variant
    has one from another dataset, nor an rsID item it repeats. A row whose lists are not
    one to one keeps the old rule. Returns the triples written (None on the old rule)."""
    if triples is None:
        row[snp_info] = ",".join(sorted(set(row[snp_info].split(","))))
        row[snp_info - 2] = ",".join(sorted(set(row[snp_info - 2].split(","))))
        row[snp_info - 1] = ",".join(sorted(set(row[snp_info - 1].split(","))))
        return None
    pairs = {}
    for tok, rsid, af in triples:
        have = pairs.setdefault(tok, [])
        if (rsid, af) not in have:
            have.append((rsid, af))
    any_af = any(af is not None for _tok, _rsid, af in triples)
    out = []
    for tok in sorted(pairs):
        have = pairs[tok]
        with_af = [p for p in have if p[1] is not None]
        named = {p[0] for p in with_af}
        kept = with_af + [p for p in have
                          if p[1] is None and p[0] != "." and p[0] not in named]
        for rsid, af in sorted(kept or have[:1], key=lambda p: (p[0], p[1] or "")):
            out.append((tok, rsid, af))
    row[snp_info] = ",".join(sorted(pairs))
    row[snp_info - 2] = ",".join(rsid for _tok, rsid, _af in out)
    row[snp_info - 1] = (",".join("NA" if af is None else af for _tok, _rsid, af in out)
                         if any_af else "")
    return out


def _merge_row(kept, ele, snp_info, triples):
    """Add a merged row's samples, rsIDs and AFs to the row that stays (_merged_list: a
    sites-only dataset's empty lists add no empty item). ``triples`` holds, by row (id),
    the (variant, rsID, AF) of the rows merged into it (_triples), which the merge
    writes back one to one (_write_variant_lists)."""
    kept[true_guide - 2] = _merged_list(
        kept[true_guide - 2], ele[true_guide - 2]
    )  # true_guide - 2 points to samples column
    mine = triples[id(kept)] if id(kept) in triples else _triples(kept, snp_info)
    theirs = triples[id(ele)] if id(ele) in triples else _triples(ele, snp_info)
    triples[id(kept)] = mine + theirs if mine is not None and theirs is not None else None
    kept[snp_info - 2] = _merged_list(
        kept[snp_info - 2], ele[snp_info - 2]
    )  # snp_info - 2 points to rsID column
    kept[snp_info - 1] = _merged_list(
        kept[snp_info - 1], ele[snp_info - 1]
    )  # snp_info - 1 points to AF column


def _fold_same_carriers(rows, fold_of, snp_info, triples):
    """Alignments of one site naming the same variants with the same carriers are one
    entry, as for every other row: the first of them in ``rows``, which is sorted by the
    file's own criterion, so the entry is the best by it. It takes the rsIDs and AFs of
    the others. ``fold_of`` maps a row (by id) to its (position, variants, carriers);
    other rows stay as they are. The carrier sets do not depend on the sort, so the
    three files keep one row count per cluster and their alternative alignments stay
    aligned line by line."""
    kept = {}
    out = []
    for row in rows:
        ident = fold_of.get(id(row))
        if ident is None:
            out.append(row)
        elif ident in kept:
            _merge_row(kept[ident], row, snp_info, triples)
        else:
            kept[ident] = row
            out.append(row)
    for row in kept.values():
        row[true_guide - 2] = ",".join(sorted(set(row[true_guide - 2].split(","))))
        triples[id(row)] = _write_variant_lists(
            row, snp_info, triples[id(row)] if id(row) in triples else _triples(row, snp_info))
    return out


def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    # avoid crush when cluster is empty in the first call
    if not cluster:
        return
    list_ref = []
    seen_ref = set()
    dict_var = {}
    triples = {}  # a variant row (id) -> (variant, rsID, AF) of it and the rows merged in
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
            # merge samples of identical targets (coming from different VCF datasets)
            key = (ele[pos], ele[snp_info])
            # An indel row lists the carriers of its own sequence (each carrier haplotype
            # resolved on its own, analisi_indels_NNN.py), and two alignments at one site
            # can name the same variants with different carriers: an RNA-bulge alignment
            # does not reach the base a SNP sits on, so it holds carriers of that SNP
            # that the full alignment does not. So only rows with one sequence merge here.
            if _holds_indel_token(ele[snp_info]):
                key += (ele[2].upper(),)
            if key in dict_var.keys():
                _merge_row(dict_var[key][0], ele, snp_info, triples)
            else:
                dict_var[key] = [ele]
    # Alignments of one site naming the same variants with the same carriers are folded
    # into one entry once the variant rows are sorted (_fold_same_carriers).
    fold_of = {}
    for key in dict_var.keys():
        if len(key) == 3:
            row = dict_var[key][0]
            fold_of[id(row)] = (key[0], key[1], frozenset(row[true_guide - 2].split(",")))

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
        # processes (bare set() iteration is PYTHONHASHSEED-dependent). The variant,
        # rsID and AF columns go together, each rsID and AF in its variant's place.
        triples[id(target)] = _write_variant_lists(
            target, snp_info,
            triples[id(target)] if id(target) in triples else _triples(target, snp_info))
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
            final_list_best_var = _fold_same_carriers(final_list_best_var, fold_of, snp_info,
                                                     triples)

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
            final_list_best_var = _fold_same_carriers(final_list_best_var, fold_of, snp_info,
                                                     triples)

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
