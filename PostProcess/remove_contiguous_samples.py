#!/usr/bin/env python

import heapq
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

    2.10.0: two alignments of one indel target string can share their aligned DNA too, since
    the merge keys an indel row by its aligned spacer as well. Those with one carrier set are
    folded into one entry (_fold_same_carriers), and an exact tie between them keeps the
    input order, in which the merge kept the first of them before 2.10.0.
    """
    try:
        _samples = str(x[true_guide - 2])
        _ncarriers = (
            0 if _samples in ("", ".", "NA", "n", "NO_SAMPLES")
            else len([_c for _c in _samples.split(",") if _c and _c not in (".", "NA", "NO_SAMPLES")])
        )
    except Exception:  # never let the tie-break break the sort
        _ncarriers = 0
    # (what an indel token is, as the merge's carrier resolution reads it: RV3-D)
    _has_indel = _ncarriers > 0 and _holds_indel_token(str(x[snp_info]))
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
    insertion adds. Since 2.10.0 such a row starts at the insertion's anchor
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
    2.10.0 the merge keys an indel row by it as well, so every alignment of a target comes
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


def _fold_same_carriers(rows, fold_of, snp_info, triples):
    """Alignments of one site naming the same variants with the same carriers are one
    entry, as for every other row: the first of them in ``rows``, which is sorted by the
    file's own criterion, so the entry is the best by it. It takes the rsIDs and AFs of
    the others. ``fold_of`` maps a row (by id) to its (cluster position on the haplotype
    carrying its indel, variants, carriers); other rows stay as they are. The carrier
    sets do not depend on the sort, so the three files keep one row count per cluster
    and their alternative alignments stay aligned line by line.

    The cluster position is the one on the haplotype (_haplotype_cluster_position): one
    site's alignments there, the bulged ones folding into the best, as up to 2.8.1. Since
    2.10.0 alignments that start on different bases of one insertion share its anchor and
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
            # merge samples of identical targets (coming from different VCF datasets).
            # Since 2.10.0 a row that starts on a base an insertion adds starts at its
            # anchor, so two copies of one target inside one insertion (a repeat it
            # holds), the same bases with the same carriers, are such identical rows too:
            # one entry, at the anchor (2.8.1 reported each at its own start).
            key = (ele[pos], ele[snp_info])
            # An indel row lists the carriers of its own sequence (each carrier haplotype
            # resolved on its own, analisi_indels_NNN.py), and two alignments at one site
            # can name the same variants with different carriers: an RNA-bulge alignment
            # does not reach the base a SNP sits on, so it holds carriers of that SNP
            # that the full alignment does not. So only rows with one sequence merge here.
            # And only rows with one alignment of it (2.10.0): a carrier resolved PUTATIVE
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
                kept = dict_var[key][0]
                _merge_row(kept, ele, snp_info, triples)
                # Two copies of one target inside one insertion are this one row, and the
                # indel post-analysis judges each copy's PAM on its own: one that a
                # placement of the insertion moves onto the reference's own PAM is not
                # created, the copy the insertion adds can be. The row is created when
                # either copy is (METHODS, "PAM_creation"), in whatever order the copies
                # come (2.10.0; the pre-merge sort's whole-line tie-break put a created
                # copy, its PAM in upper case, before an "n" one). PAM_gen is the column
                # after Total; one alignment's created PAM is the same string in every
                # copy, its ALT PAM.
                if len(key) == 4 and kept[total + 1] == "n" and ele[total + 1] != "n":
                    kept[total + 1] = ele[total + 1]
            else:
                dict_var[key] = [ele]
    # Alignments of one site naming the same variants with the same carriers are folded
    # into one entry once the variant rows are sorted (_fold_same_carriers), the site
    # being their cluster position on the haplotype carrying the indel, read once per
    # target, REF text, place and variants (the memo), not once per aligned spacer.
    fold_of = {}
    first_base_memo = _FIRST_BASE_MEMO
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


# The offsets of a row's first base past its start (_haplotype_cluster_position), read once
# per target and REF text: the merge's placement reads them, and the fold of each cluster
# again. Emptied at each guide and chromosome.
_FIRST_BASE_MEMO = {}


def _inside_place(row):
    """(anchor, place on the haplotype) of a row whose first base is one an insertion adds;
    None for any other row. Since 2.10.0 such a row starts at the insertion's anchor
    (indel_snp_cis.reference_start): its Position is the anchor whatever its bulges, and its
    place on the haplotype carrying the insertion (_haplotype_cluster_position) is past its
    Cluster_Position by the rank of its first base among the inserted ones, the same for
    one target's alignments before one PAM."""
    if not _holds_indel_token(row[snp_info]):
        return None
    place = _haplotype_cluster_position(row, _FIRST_BASE_MEMO)
    if place <= int(row[pos]):
        return None
    return int(row[pos - 1]), place


def _indel_tokens(variants):
    """The tokens of a row's variant list (<contig>_<pos>_<ref>_<alt>,...) whose REF and ALT
    differ in length, its indel's, in the list's order: the first is the one its place on
    the haplotype is read from (_first_base_offset). Split from the right, as
    _holds_indel_token reads them."""
    out = []
    for tok in variants.split(","):
        p = tok.rsplit("_", 3)
        if len(p) == 4 and len(p[2]) != len(p[3]):
            out.append(tok)
    return out


def _max_shift(path):
    """The largest difference between a row's Cluster_Position and its Position in the
    merge's input ``path``: the most an alignment's bulges move its Cluster_Position off its
    start, and so the most a row that starts inside an insertion is off the insertion's
    anchor. _clusters reads the rows by their place in the chain through a window this
    wide. One pass over the file, reading two columns of each line."""
    shift = 0
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.split("\t", pos + 1)
            try:
                d = abs(int(f[pos]) - int(f[pos - 1]))
            except (IndexError, ValueError):
                continue
            if d > shift:
                shift = d
    return shift


def _clusters(lines, shift=0):
    """The merge clusters of ``lines`` (sorted by guide, chromosome and Cluster_Position),
    each a list of rows in file order, in the order of their places.

    Rows within --merge bp (tau) of the previous one are one cluster, chained: one row of a
    population-level search. Since 2.10.0 a row that starts on a base an insertion adds
    starts at the insertion's anchor (indel_snp_cis.reference_start), up to the
    insertion's length before where 2.8.1 placed it, and the alignments that start inside
    one insertion are one site, at the anchor (decision (d), 2026-11-06: such a site can
    join a reference row's there).

    Chained like any other row, such a row bridged two clusters 2.8.1 kept apart, and a
    site of its own left the results for a neighbour's cluster: on hg38_mega chr22,
    population-level, chr22:36,838,374 + of VEGFA site 3 and chr22:47,379,217 + of EMX1
    (the insertion rows chr22_36838378_A_ATGTG... and chr22_47379212_A_AGGAGG... moved to
    their anchors, between them and the next cluster). Placed one by one, each by its own
    Cluster_Position, they split one target in two: its alignments before one PAM with
    different bulges have different Cluster_Positions (on the + strand with a 3' PAM a DNA
    bulge adds one, an RNA bulge takes one off), so with a site on either side the plain
    alignment joined one and its bulged alignment the other, and the fold
    (_fold_same_carriers), which works within a cluster, wrote both: chr22:45,890,432 (+)
    of VEGFA site 3 twice on HGDP + 1000 Genomes, CFD 0.179 and 0.172, one insertion and
    one carrier.

    So the chain is made of the other rows, and the rows that start inside an insertion
    at one anchor are placed together, one unit, which joins the cluster whose rows lie
    within tau of its rows (of their Cluster_Positions, the smallest to a cluster before
    the unit, the largest to one after it) -- of two, the one nearer its place on the
    haplotype (the nearest of its rows'), the next on a tie, where 2.8.1 had it -- or, with
    none, a cluster with the other such units within tau of it: a unit joins one site and
    never chains two. Its rows reach the merge in the order of their Cluster_Positions, up
    to a bulge's length (``shift``, _max_shift) either side of their anchor, so the rows are
    taken in the order of their place in the chain, the anchor for a unit's row and the
    Cluster_Position for any other, through a window of ``shift`` bp: every row of a unit
    is in when it is placed. A row's place, and so its cluster, does not depend on the
    order of the rows of one Cluster_Position, which the three merges sort differently
    (each must make the same clusters: their alternative alignments are pasted line by
    line). Without such rows the clusters are those of the plain chain, as up to 2.8.1.

    A unit is made of whole targets (review r2 of lane M). One target's alignments before
    one PAM start on different bases, so some can start inside the insertion and others on
    a reference base before it (on the + strand with a 3' PAM a DNA bulge starts one base
    earlier, an RNA bulge one base later): the anchor base itself, or a base before it.
    Such an alignment names the insertion and is a row of the chain like any other, at its
    Cluster_Position, which is the target's place on the haplotype, where the alignments
    inside are too (_haplotype_cluster_position). With those held in the unit, the unit
    could join the site before it while the target's alignment on a reference base opened
    the site after it, and the target was written twice: at 2 DNA + 2 RNA bulges, another
    target's 2-RNA-bulge alignment from the insertion's first base, one base after a site,
    took the unit there, and the target's 2-DNA-bulge alignment from the anchor base was a
    site of its own (the chain of every row, and 2.8.1, wrote it once). So a target that
    has an alignment on a reference base is placed by it: its alignments inside the
    insertion (the same strand, place on the haplotype and insertion) leave the unit and
    go with it, at its Cluster_Position, where 2.8.1 placed them all; the unit holds the
    targets that lie wholly inside the insertion. They reach the merge before it (their
    Cluster_Positions are smaller, by the rank of their first base among the inserted
    ones) and are still in the window then (it starts at the anchor or before, at most
    ``shift`` bp before its Cluster_Position), so every alignment of a target is in one
    cluster, whatever its bulges and --merge, and a row's cluster still does not depend on
    the order of the rows of one Cluster_Position."""
    group = None
    # (place in the chain, 1 for a unit's row, seq, row, cp, haplotype place, its target:
    # (strand, haplotype place, insertion) for a unit's row, None for a row of the chain)
    window = []
    current = None   # [last place of its chain, [(seq, row)]]
    pending = []     # [anchor, lowest cp, highest cp, nearest haplotype place, [(seq, row)]]
    inside = {}      # a target inside an insertion -> its rows in the window, in a unit
    moved = set()    # the seq of a unit's row that went with its target's row of the chain
    seq = 0

    def place(nxt_cp):
        """Clusters to emit before a chain row of Cluster_Position nxt_cp (None: the group's
        end) starts a cluster of its own, and the rows of the ``pending`` units that join
        it."""
        out, joins, alone = [], [], []
        for unit in pending:
            _at, lo, hi, hp, rows = unit
            to_cur = current is not None and lo - current[0] <= tau
            to_next = nxt_cp is not None and nxt_cp - hi <= tau
            if to_cur and to_next:
                if abs(nxt_cp - hp) <= abs(hp - current[0]):
                    joins.extend(rows)
                else:
                    current[1].extend(rows)
            elif to_cur:
                current[1].extend(rows)
            elif to_next:
                joins.extend(rows)
            else:
                alone.append(unit)
        if current is not None:
            out.append(current[1])
        chain, last = [], None
        for _at, lo, hi, _hp, rows in alone:
            if chain and lo - last > tau:
                out.append(chain)
                chain = []
            chain.extend(rows)
            last = hi if last is None else max(last, hi)
        if chain:
            out.append(chain)
        return out, joins

    def take(item):
        """Take one row in its place in the chain; return the clusters it completes."""
        nonlocal current, pending
        at, in_unit, s, row, cp, hp, _target = item
        if in_unit:
            if pending and pending[-1][0] == at:
                unit = pending[-1]
                unit[1], unit[2], unit[3] = min(unit[1], cp), max(unit[2], cp), min(unit[3], hp)
                unit[4].append((s, row))
            else:
                pending.append([at, cp, cp, hp, [(s, row)]])
            return []
        if current is not None and at - current[0] <= tau:
            # the chain goes on: the units held lie within it
            for unit in pending:
                current[1].extend(unit[4])
            pending = []
            current[0] = at
            current[1].append((s, row))
            return []
        done, joins = place(at)
        pending = []
        current = [at, joins + [(s, row)]]
        return done

    def pop():
        """Take the next row of the window; return the clusters it completes. A unit's row
        that went with its target's row of the chain is taken there, not in the unit; once
        a unit's row is taken, no row of the chain to come is its target's."""
        item = heapq.heappop(window)
        if item[1]:
            if item[2] in moved:
                moved.discard(item[2])
                return []
            inside.pop(item[6], None)
        return take(item)

    def drain():
        """The clusters of the rest of a guide and chromosome."""
        out = []
        while window:
            out.extend(pop())
        out.extend(place(None)[0])
        return out

    for line in lines:
        row = line.split("\t")
        key = (row[true_guide], row[chrom])
        if key != group:
            if group is not None:
                for c in drain():
                    yield [r for _s, r in sorted(c, key=lambda x: x[0])]
            current, pending, group = None, [], key
            inside.clear()
            moved.clear()
            _FIRST_BASE_MEMO.clear()
        cp = int(row[pos])
        at = _inside_place(row)
        if at is None:
            # a row of the chain; one naming an insertion can be an alignment of a target
            # whose other alignments start inside it: they go with it, at its place
            for tok in (_indel_tokens(row[snp_info]) if _holds_indel_token(row[snp_info])
                        else ()):
                for item in inside.pop((row[pos + 1], cp, tok), ()):
                    moved.add(item[2])
                    heapq.heappush(window, (cp, 0, item[2], item[3], item[4], item[5], None))
            heapq.heappush(window, (cp, 0, seq, row, cp, cp, None))
        else:
            toks = _indel_tokens(row[snp_info])
            target = (row[pos + 1], at[1], toks[0] if toks else None)
            item = (at[0], 1, seq, row, cp, at[1], target)
            inside.setdefault(target, []).append(item)
            heapq.heappush(window, item)
        seq += 1
        # every row still to come is at cp - shift or after in the chain
        while window and window[0][0] < cp - shift:
            for c in pop():
                yield [r for _s, r in sorted(c, key=lambda x: x[0])]
    if group is not None:
        for c in drain():
            yield [r for _s, r in sorted(c, key=lambda x: x[0])]


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
shift = _max_shift(sys.argv[1])
with open(sys.argv[1], "r") as fileIn:
    header = fileIn.readline()
    with open(sys.argv[2], "w") as fileOut:
        with open(sys.argv[2] + ".discarded_samples", "w") as fileOut_disc:
            fileOut.write(header)
            fileOut_disc.write(header)
            for cluster in _clusters(fileIn, shift):
                get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)

print("Merging done in: " + str(time.time() - start))
