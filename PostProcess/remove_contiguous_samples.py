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
# per target and REF text: the merge's chain reads them, and the fold of each cluster again.
# Emptied at each guide and chromosome.
_FIRST_BASE_MEMO = {}


def _inside_offset(row):
    """How far past its Cluster_Position a row is on the haplotype carrying its indel
    (_haplotype_cluster_position): more than 0 for a row whose first base is one an
    insertion adds, which since 2.10.0 starts at the insertion's anchor; 0 for any other."""
    if not _holds_indel_token(row[snp_info]):
        return 0
    return max(0, _haplotype_cluster_position(row, _FIRST_BASE_MEMO) - int(row[pos]))


def _clusters(lines):
    """The merge clusters of ``lines`` (sorted by guide, chromosome and Cluster_Position),
    each a list of rows in file order, in the order of their Cluster_Positions.

    Rows within --merge bp (tau) of the previous one are one cluster, chained: one row of a
    population-level search. Since 2.10.0 a row that starts on a base an insertion adds
    starts at the insertion's anchor (indel_snp_cis.reference_start), up to the
    insertion's length before where 2.8.1 placed it, and the alignments that start inside
    one insertion are one site, at the anchor (decision (d), 2026-11-06: an indel row can
    join a reference row's site there). Chained like any other row, such a row bridged two
    clusters 2.8.1 kept apart, and a site of its own left the results for a neighbour's
    cluster: on hg38_mega chr22, population-level, chr22:36,838,374 + of VEGFA site 3 and
    chr22:47,379,217 + of EMX1 (the insertion rows chr22_36838378_A_ATGTG... and
    chr22_47379212_A_AGGAGG... moved to their anchors, between them and the next cluster).
    So the chain is made of the other rows, and a row starting inside an insertion joins
    the cluster within tau of its Cluster_Position -- of two, the one nearer its place on
    the haplotype, the next on a tie (where 2.8.1 had it) -- or, with none, a cluster with
    the others of its kind within tau of it: it joins a site, and never chains two."""
    group = None
    current = None   # [last Cluster_Position of its chain, [(seq, row)]]
    pending = []     # (seq, row, Cluster_Position, place on the haplotype) to be placed
    seq = 0

    def place(nxt_cp):
        """Clusters to emit before a chain row of Cluster_Position nxt_cp (None: the group's
        end) starts a cluster of its own, and the rows of ``pending`` that join it."""
        out, joins, alone = [], [], []
        for item in pending:
            _s, _r, cp, hp = item
            to_cur = current is not None and cp - current[0] <= tau
            to_next = nxt_cp is not None and nxt_cp - cp <= tau
            if to_cur and to_next:
                if abs(nxt_cp - hp) <= abs(hp - current[0]):
                    joins.append(item)
                else:
                    current[1].append(item[:2])
            elif to_cur:
                current[1].append(item[:2])
            elif to_next:
                joins.append(item)
            else:
                alone.append(item)
        if current is not None:
            out.append(current[1])
        chain, last = [], None
        for item in alone:
            if chain and item[2] - last > tau:
                out.append(chain)
                chain = []
            chain.append(item[:2])
            last = item[2]
        if chain:
            out.append(chain)
        return out, [j[:2] for j in joins]

    for line in lines:
        row = line.split("\t")
        key = (row[true_guide], row[chrom])
        if key != group:
            if group is not None:
                done, _joins = place(None)
                for c in done:
                    yield [r for _s, r in sorted(c, key=lambda x: x[0])]
            current, pending, group = None, [], key
            _FIRST_BASE_MEMO.clear()
        cp = int(row[pos])
        if _inside_offset(row):
            pending.append((seq, row, cp, cp + _inside_offset(row)))
        elif current is not None and cp - current[0] <= tau:
            # the chain goes on: the inside rows held lie within it
            current[1].extend(item[:2] for item in pending)
            pending = []
            current[0] = cp
            current[1].append((seq, row))
        else:
            done, joins = place(cp)
            for c in done:
                yield [r for _s, r in sorted(c, key=lambda x: x[0])]
            pending = []
            current = [cp, joins + [(seq, row)]]
        seq += 1
    if group is not None:
        done, _joins = place(None)
        for c in done:
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
with open(sys.argv[1], "r") as fileIn:
    header = fileIn.readline()
    with open(sys.argv[2], "w") as fileOut:
        with open(sys.argv[2] + ".discarded_samples", "w") as fileOut_disc:
            fileOut.write(header)
            fileOut_disc.write(header)
            for cluster in _clusters(fileIn):
                get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)

print("Merging done in: " + str(time.time() - start))
