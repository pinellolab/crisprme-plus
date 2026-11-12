"""The alleles of one variant position, each read against its own genotype index.

CRISPRitz add-variants lists every SNP ALT of a position in the per-sample dictionary
(``my_dict_<chrom>.json``), one '$' entry per ALT, each with its own record's carriers
and their genotypes. A VCF split with ``bcftools norm -m -any`` (how every published panel
is built) has one record per ALT, and its genotypes use 1 for that ALT. CRISPRitz 2.10.0
writes the entries of a record that lists several ALTs on one line (C>A,T) in that same
split form, 1 each entry's own ALT and 0 the reference or another ALT (S:2|1 is 0|1 under A
and 1|0 under T), so every position of its dictionaries comes back as it is (below). A
dictionary of an earlier CRISPRitz (2.9.1, or a 2.10.0 build up to c18c575) keeps such a
record's raw genotype under every ALT its sample carries: S:2|1 is listed under A and under
T, and the rest of this note is about reading it.
CRISPRme reads "1" as "this slot holds this ALT" (the Tier-0 counts, the Tier-1
haplotypes, the dictionary path's haplotype seeds, the indel path's slots), so each
ALT's genotypes are first read against that ALT's own allele index: S:2|1 becomes 0|1
under A and 1|0 under T.

The allele index of an ALT is the token all its carriers share. When that leaves more
than one (every carrier holds two ALT indexes, e.g. 2|1), it is the one no other ALT of
the position needs: the indexes the position's other entries of a line of several ALTs
settle are taken off it. An entry whose carriers hold no index other than 0 and 1 settles
none, for its "1" may be another line's numbering: CRISPRitz 2.10.0 lists every line of a
position (a split record, the listing it merges from two lines, written as a split
record's) beside a line of several ALTs, whose raw calls number its own ALTs, so C>G
(1|0) beside C>A,CT (1|2) would have tied A to 2, the insertion's number, and read A on
the other slot. Within one line of several ALTs such an entry never settles a number
another entry needs (a carrier of its ALT and of another one-base ALT is listed under
both, with both numbers), so the narrowing loses nothing there. Where the genotypes
cannot say (every carrier of two ALTs is 1|2, or of an SNV and an indel or a '*' of its
line), no slot can be tied to the ALT: each carrier reads "1/0", one copy on an unknown
slot, and is resolved without phase (PUTATIVE). A position whose genotypes hold no index
other than 0 and 1 (split or biallelic records) comes back as it is.

Two records with the same REF>ALT at one position are one allele: their carriers are
merged as CRISPRitz 2.10.0's enricher merges them into the allele's one dictionary entry
(its ``_gt_union``): two calls phased alike give a slot the ALT where either puts it there
(``merge_gt``), and two that are not (one unphased, or of another ploidy) give the call with
more copies of the ALT, the earlier on a tie, as it is. A slot that ends up holding two
different ALTs of the position is no haplotype; the readers resolve such a sample without
phase (PUTATIVE).

STDLIB ONLY and pure, like ``simple_analysis_registry``: imported by the dictionary and
tier readers (tier0_compile, simple_analysis_registry, indel_snp_cis, observed_haplotypes,
population_summary) and unit-tested on its own (test_site_alleles.py).
"""

import re

# the tokens of a split or biallelic record's genotypes ("" an empty genotype)
_PLAIN = ("0", "1", ".", "")
# the separators between the alleles of a genotype, either kind
_ALLELE_SEP = re.compile(r"[|/]")


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


def _numbers_another_line(carriers):
    """True when a carrier's genotype holds an allele index other than 0 and 1: the
    entry's genotypes are a line of several ALTs' raw calls (a split record's, and the
    entry CRISPRitz 2.10.0 merges from two records, hold only 0 and 1)."""
    return any(t not in _PLAIN for g in carriers.values() for t in gt_tokens(g)[0])


def alt_tokens(entries):
    """The allele index of each entry of one position, or None where the genotypes
    cannot tell it.

    ``entries`` is a list of ``(alt, {sample: genotype})``, one per ALT record of the
    position (the same ALT may come twice, from two records). Returns a list of the
    same length: the token every carrier of the entry shares, narrowed by the tokens
    the position's other entries of a line of several ALTs settle when more than one is
    left (an entry whose genotypes hold only 0 and 1 settles none: see the module
    docstring); None when still more than one is left (every carrier holds both), or
    none; "1" for an entry with no carriers."""
    shared = []
    numbered = []
    for _alt, carriers in entries:
        common = None
        for g in carriers.values():
            toks = {t for t in gt_tokens(g)[0] if t not in ("0", ".", "")}
            common = toks if common is None else common & toks
        shared.append(common)
        numbered.append(_numbers_another_line(carriers))
    settled = True
    while settled:
        settled = False
        taken = {next(iter(c)) for c, n in zip(shared, numbered)
                 if n and c is not None and len(c) == 1}
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
    """One sample's genotype for one ALT that two records give, ``a`` the earlier's (read
    against the ALT: "1" holds it), as CRISPRitz 2.10.0's enricher merges the two into
    the allele's one dictionary entry (its ``_gt_union``), so that a reader merging two
    listings reads them as the enricher would have. Two phased genotypes of one ploidy
    (or two haploid ones) are merged slot by slot: a slot holds the ALT if either record
    puts it there, else a called 0 if either has one. Otherwise the slots cannot be
    matched up, and the genotype with more copies of the ALT is kept as it is, the
    earlier on a tie (1|0 beside 0/1 stays 1|0; up to 2.10.0's integration this read
    "1/0", unphased, where the enricher kept 1|0)."""
    ta, ka = gt_tokens(a)
    tb, kb = gt_tokens(b)
    if ka == kb and ka != "unphased" and len(ta) == len(tb) and "/" not in a + b:
        merged = []
        for x, y in zip(ta, tb):
            if x == "1" or y == "1":
                merged.append("1")
            elif x == "0" or y == "0":
                merged.append("0")
            else:
                merged.append(x)
        return merged[0] if ka == "haploid" else "|".join(merged)
    copies_a = _ALLELE_SEP.split(a.strip()).count("1")
    copies_b = _ALLELE_SEP.split(b.strip()).count("1")
    return b if copies_b > copies_a else a


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
    their carriers merged (``merge_carriers``, a sample's two calls as ``merge_gt`` merges
    them), with the rsID and AF of the more frequent record (on a tie, the earlier), as
    CRISPRitz 2.10.0's enricher takes them. The REF is read in upper case, as VCF 4.x
    reads it: CRISPRitz codes C>G and c>G alike but lists them as two entries of G (its
    dictionary compares an SNV's REF,ALT as written), which stopped the Tier-1 compile
    ("duplicate (pos, alt) record"); the record is named as the first writes it. The ALT
    is read as written: CRISPRitz codes no ALT written in lower case. Returns the same
    shape, one record per REF>ALT in first-seen order: ``records`` itself when the
    position is plain and repeats no REF>ALT."""
    keys = [(r[0].upper(), r[1]) for r in records]
    entries = [(r[1], r[4]) for r in records]
    got = per_alt_entries(entries)
    if got is entries and len(set(keys)) == len(keys):
        return records
    out = []
    where = {}
    for key, (ref, alt, rsid, af, _c), (_alt, carriers) in zip(keys, records, got):
        i = where.get(key)
        if i is None:
            where[key] = len(out)
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
