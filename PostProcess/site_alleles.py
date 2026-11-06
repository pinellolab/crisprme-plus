"""The alleles of one variant position, each read against its own genotype index.

CRISPRitz add-variants lists every SNP ALT of a position in the per-sample dictionary
(``my_dict_<chrom>.json``), one '$' entry per ALT, each with its own record's carriers
and their genotypes as the VCF writes them. A VCF split with ``bcftools norm -m -any``
(how every published panel is built) has one record per ALT, and its genotypes use 1
for that ALT. A record that lists several ALTs on one line (C>A,T) keeps its raw
genotype under every ALT its sample carries: S:2|1 is listed under A and under T.
CRISPRme reads "1" as "this slot holds this ALT" (the Tier-0 counts, the Tier-1
haplotypes, the dictionary path's haplotype seeds, the indel path's slots), so each
ALT's genotypes are first read against that ALT's own allele index: S:2|1 becomes 0|1
under A and 1|0 under T.

The allele index of an ALT is the token all its carriers share. When that leaves more
than one (every carrier holds two ALT indexes, e.g. 2|1), it is the one no other ALT of
the position needs. Where the genotypes cannot say (every carrier of two ALTs is 1|2),
no slot can be tied to the ALT: each carrier reads "1/0", one copy on an unknown slot,
and is resolved without phase (PUTATIVE). A position whose genotypes hold no index
other than 0 and 1 (split or biallelic records) comes back as it is.

Two records with the same REF>ALT at one position are one allele: their carriers are
merged, a slot carrying the ALT if either record puts it there (``merge_gt``). A slot
that ends up holding two different ALTs of the position is no haplotype; the readers
resolve such a sample without phase (PUTATIVE).

STDLIB ONLY and pure, like ``simple_analysis_registry``: imported by the dictionary and
tier readers (tier0_compile, simple_analysis_registry, indel_snp_cis, observed_haplotypes,
population_summary) and unit-tested on its own (test_site_alleles.py).
"""

# the tokens of a split or biallelic record's genotypes ("" an empty genotype)
_PLAIN = ("0", "1", ".", "")


def gt_tokens(gt):
    """(allele tokens, kind) of a genotype: kind is "phased" ('|'), "unphased" ('/')
    or "haploid" (one token, e.g. a male chrX call)."""
    gt = (gt or "").strip()
    if "|" in gt:
        return gt.split("|"), "phased"
    if "/" in gt:
        return gt.split("/"), "unphased"
    return [gt], "haploid"


def for_alt(gt, token):
    """A genotype as "this alt or not": ``token``, the allele index of the alt, becomes
    "1", any other allele "0"; a missing allele stays "."."""
    toks, kind = gt_tokens(gt)
    norm = ["1" if t == token else ("." if t == "." else "0") for t in toks]
    if kind == "haploid":
        return norm[0]
    return ("|" if kind == "phased" else "/").join(norm)


def _unplaced(gt):
    """One copy of the alt on a slot the genotypes do not tell: "1/0" (haploid "1")."""
    toks, _kind = gt_tokens(gt)
    return "/".join(["1"] + ["0"] * (len(toks) - 1))


def alt_tokens(entries):
    """The allele index of each entry of one position, or None where the genotypes
    cannot tell it.

    ``entries`` is a list of ``(alt, {sample: genotype})``, one per ALT record of the
    position (the same ALT may come twice, from two records). Returns a list of the
    same length: the token every carrier of the entry shares, narrowed by the tokens
    the position's other entries need when more than one is left; None when still
    more than one is left (every carrier holds both); "1" for an entry with no
    carriers."""
    shared = []
    for _alt, carriers in entries:
        common = None
        for g in carriers.values():
            toks = {t for t in gt_tokens(g)[0] if t not in ("0", ".", "")}
            common = toks if common is None else common & toks
        shared.append(common)
    settled = True
    while settled:
        settled = False
        taken = {next(iter(c)) for c in shared if c is not None and len(c) == 1}
        for i, c in enumerate(shared):
            if c is not None and len(c) > 1 and c & taken:
                shared[i] = c - taken
                settled = True
    out = []
    for c in shared:
        if c is None:
            out.append("1")
        elif len(c) == 1:
            out.append(next(iter(c)))
        else:
            out.append(None)
    return out


def plain(entries):
    """True when no genotype of the position holds an index other than 0 and 1: split
    records or biallelic ones, which every reader reads as they are."""
    return all(t in _PLAIN for _alt, carriers in entries for g in carriers.values()
               for t in gt_tokens(g)[0])


def per_alt_entries(entries):
    """Each entry's genotypes read against its own allele index (see the module
    docstring): a list of ``(alt, {sample: genotype})`` like ``entries``, "1" meaning
    "this slot holds this alt". ``entries`` itself when the position is plain."""
    if plain(entries):
        return entries
    out = []
    for (alt, carriers), token in zip(entries, alt_tokens(entries)):
        if token is None:
            out.append((alt, {s: _unplaced(g) for s, g in carriers.items()}))
        else:
            out.append((alt, {s: for_alt(g, token) for s, g in carriers.items()}))
    return out


def merge_gt(a, b):
    """One sample's genotype for one ALT that two records give (read against the ALT:
    "1" holds it): a slot holds the ALT if either record puts it there. Two phased
    genotypes of one ploidy (or two haploid ones) are merged slot by slot; otherwise
    which slot holds it is not known, and the sample reads unphased with the larger of
    the two copy counts."""
    ta, ka = gt_tokens(a)
    tb, kb = gt_tokens(b)
    if ka == kb and ka != "unphased" and len(ta) == len(tb):
        merged = []
        for x, y in zip(ta, tb):
            if x == "1" or y == "1":
                merged.append("1")
            elif x == "0" or y == "0":
                merged.append("0")
            else:
                merged.append(".")
        return merged[0] if ka == "haploid" else "|".join(merged)
    ploidy = max(len(ta), len(tb))
    copies = min(ploidy, max(ta.count("1"), tb.count("1")))
    if ploidy == 1:
        return "1" if copies else "0"
    return "/".join(["1"] * copies + ["0"] * (ploidy - copies))


def merge_carriers(a, b):
    """The carriers of one ALT that two records give, as one {sample: genotype}."""
    out = dict(a)
    for s, g in b.items():
        out[s] = merge_gt(out[s], g) if s in out else g
    return out


def per_alt(site):
    """{alt: {sample: genotype}} at one position, each genotype read against its own
    alt (``per_alt_entries`` on a mapping that holds each ALT once). ``site`` itself
    when the position is plain."""
    entries = list(site.items())
    got = per_alt_entries(entries)
    if got is entries:
        return site
    return dict(got)


def _af_number(af):
    """An allele-frequency string as a number, 0 when it does not read as one."""
    try:
        return float(af)
    except (TypeError, ValueError):
        return 0.0


def read_records(records):
    """The records of one position as the readers use them.

    ``records`` is a list of ``(ref, alt, rsid, af, {sample: genotype})``, one per
    dictionary entry in its order. Each entry's genotypes are read against its own
    allele index (``per_alt_entries``), then two records of one REF>ALT become one:
    their carriers merged (``merge_carriers``), with the rsID and AF of the more frequent
    record (on a tie, the earlier), as CRISPRitz 2.10.0's enricher merges them. Returns
    the same shape, one record per REF>ALT in first-seen order: ``records`` itself when
    the position is plain and repeats no REF>ALT."""
    keys = [(r[0], r[1]) for r in records]
    entries = [(r[1], r[4]) for r in records]
    got = per_alt_entries(entries)
    if got is entries and len(set(keys)) == len(keys):
        return records
    out = []
    where = {}
    for (ref, alt, rsid, af, _c), (_alt, carriers) in zip(records, got):
        i = where.get((ref, alt))
        if i is None:
            where[(ref, alt)] = len(out)
            out.append((ref, alt, rsid, af, carriers))
            continue
        pref, palt, prsid, paf, pcarriers = out[i]
        if _af_number(af) > _af_number(paf):
            prsid, paf = rsid, af
        out[i] = (pref, palt, prsid, paf, merge_carriers(pcarriers, carriers))
    return out


def two_alts_on_a_slot(carriers_by_alt):
    """The samples a phased slot of which holds two different ALTs of one position:
    {sample: set of such slots}. ``carriers_by_alt`` is {alt: {sample: genotype}}, the
    genotypes read against each alt (``per_alt``). Such a slot is no haplotype: the
    split records were phased one by one (672 sample-positions on chr22 of HGDP + 1000
    Genomes), or two records of one ALT were merged."""
    slots = {}
    out = {}
    for alt, carriers in carriers_by_alt.items():
        for s, g in carriers.items():
            toks, kind = gt_tokens(g)
            if kind != "phased":
                continue
            for h, t in enumerate(toks):
                if t != "1":
                    continue
                held = slots.setdefault((s, h), set())
                held.add(alt)
                if len(held) > 1:
                    out.setdefault(s, set()).add(h)
    return out
