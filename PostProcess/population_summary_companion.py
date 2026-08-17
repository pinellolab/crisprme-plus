"""ADDITIVE companion population-summary writer (CRISPRme+ dictless redesign,
Phase 3c wiring).

This module surfaces the combination-aware population summary
(``population_summary.summarize``) as a PIPELINE OUTPUT, the SAFE ADDITIVE way: it
writes a SEPARATE companion TSV next to the bestMerge / output stem, one row per
VARIANT off-target, keyed by the SAME identity columns the bestMerge uses so it can
be joined back. It NEVER touches the existing bestMerge/altMerge/integrated_results
columns or the Samples column, and when no Tier-0 registry is present it writes
NOTHING (byte-identical legacy behavior).

WHY a SEPARATE module (mirrors ``simple_analysis_registry``): the writer is a PURE,
importable function (no ``sys.argv`` / module-level file I/O) so it is unit-testable
in isolation. ``new_simple_analysis.py`` calls ``write_companion`` once per
chromosome, passing the module-level readers (``myreg`` / ``mygt``) and the list of
variant off-targets it already finalized; any error is caught by the caller (log +
skip) so the companion file can NEVER break a run.

HOW off-target creating (pos, alt) sets are gathered (see new_simple_analysis.py's
``iupac_decomposition`` finalization): every VARIANT off-target row already carries a
``SNP`` column (``final_line[17]``) that is the comma-joined ``snp_info`` tokens
"<chrom>_<pos1based>_<ref>_<alt>" for the single SNP or the multi-SNP COMBINATION
that CREATES the off-target. We parse those tokens back into (pos1based, alt) pairs
and hand them straight to ``population_summary.summarize`` -- the per-alt /
combination frequency correctness is DELEGATED to that module (we never re-derive a
frequency here). A reference / PAM-only row has NO SNP token (it went through the
non-IUPAC branch) and is SKIPPED.

PHASED CHOICE (documented): the ``phased`` flag passed to ``summarize`` controls the
k>=2 combination semantics (true-cis when phased vs the assume-cis upper bound + a
lower bound when unphased). We resolve it, in order:
  1. an explicit ``phased`` arg from the caller (e.g. a build-time manifest flag);
  2. else the gt-tier manifest's ``"phased"`` entry if the build recorded one;
  3. else DETECT from the genotype strings: any '|' seen across the sampled records
     => phased; only '/' => unphased. Default when NOTHING is observed: FALSE
     (the conservative path -- unphased reports bounds, never an over-confident cis
     allele frequency). 1000G is phased ('|'), so a 1000G store resolves to True.
Mixed-phasing (some '|' some '/') resolves to the CONSERVATIVE unphased path unless
the caller / manifest overrides, and we document that a single dataset-wide bool is
passed to ``summarize`` (per-record phase is out of scope for this companion).

STDLIB ONLY. The whole module is import-guarded at its single call site so a missing
``population_summary`` in an old deploy cannot break the legacy path.
"""

from __future__ import annotations


# The companion TSV header. The first five columns are the SAME identity columns the
# bestMerge uses (Chromosome, Position, Direction, crRNA/Spacer, DNA) so a downstream
# join is a straight key match; everything after is ADDITIVE population summary.
#
# JOIN DIRECTION: join FROM bestMerge (left) INTO the companion by (Chromosome,
# Position, Direction, crRNA, DNA[, SNP]). The companion captures EVERY finalized
# variant off-target (deduped by that identity key), whereas bestMerge/bestCFD/
# bestCRISTA keep only the best-scoring representative per cluster -- so the companion
# is a SUPERSET: every bestMerge variant row has a matching companion row, but the
# companion may contain extra rows for off-targets that did not win best-scoring.
# A naive inner/right join or a row-count equality check against bestMerge would be
# surprised; use bestMerge as the left side.
COMPANION_HEADER = [
    # ---- identity (join key; mirrors bestMerge columns) ----
    "Chromosome",
    "Position",
    "Direction",
    "crRNA",
    "DNA",
    # ---- creating variant(s) ----
    "SNP",            # the comma-joined "<chrom>_<pos>_<ref>_<alt>" (verbatim)
    "n_variants",     # k == number of creating (pos, alt) pairs (1 == single)
    # ---- record-level population summary (delegated to population_summary) ----
    "global_allele_freq",
    "global_carrier_freq",
    "global_carrier_n",
    "global_hom_freq",
    "max_subpop_af",
    "max_subpop_label",
    "observed",
    "allele_freq_defined",
    "phased",
    # ---- carrier lower bound (unphased combination only) ----
    # Clearly labeled: a LOWER BOUND on the assume-cis carrier fraction. Empty for
    # the exact (single-variant or phased) path where n_carrier is exact.
    "global_carrier_freq_lower_bound",
    # ---- compact per-(db x subpop) + per-db breakdowns (one column each) ----
    # Encoding "db::subpop=af|cf|n;db2::subpop2=...": af=allele_freq (NA if undefined,
    # i.e. an unphased combination), cf=carrier_freq, n=absolute carrier count.
    "allele_freq_by_group",
    "carrier_freq_by_group",
    "carrier_n_by_group",
]


def parse_snp_field(snp_field, chrom=None):
    """Parse a bestMerge ``SNP`` column value into creating (pos1based:int, alt:str).

    ``snp_field`` is the comma-joined ``snp_info`` tokens
    "<chrom>_<pos1based>_<ref>_<alt>" (a single token for a single-SNP off-target,
    several for a multi-SNP COMBINATION). We split on ',' and, per token, parse the
    LAST two underscore-delimited fields as ref, alt and the field before them as the
    1-based position -- so a chromosome name that itself contains '_' (rare, but e.g.
    "chr1_KI270711v1_random") does not corrupt the pos/ref/alt parse.

    A token that does not parse (empty, "n", a bare "." placeholder, or missing
    pos/alt) is skipped. ``chrom`` is accepted for symmetry/validation but not
    required. Returns a list of (pos, alt) pairs in the SNP field's order (which is
    the combination order the pipeline emitted); duplicates are preserved so the
    caller sees exactly what created the row.
    """
    if snp_field is None:
        return []
    s = str(snp_field).strip()
    if not s or s.lower() == "n" or s == ".":
        return []
    pairs = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        fields = tok.split("_")
        if len(fields) < 4:
            # need at least <chrom>_<pos>_<ref>_<alt>
            continue
        alt = fields[-1].strip()
        # ref = fields[-2] (unused here; alt is what creates the off-target)
        pos_str = fields[-3].strip()
        if not alt or not pos_str:
            continue
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        pairs.append((pos, alt))
    return pairs


def resolve_phased(explicit_phased, gt_reader, observed_gt_strings=None):
    """Resolve the single dataset-wide ``phased`` bool for ``summarize``.

    Order (see the module docstring):
      1. an explicit bool from the caller wins;
      2. else the gt-tier manifest's ``"phased"`` entry if present;
      3. else DETECT from the sampled genotype strings ('|' => phased, only '/' =>
         unphased); default FALSE (conservative) when nothing is observed.

    ``observed_gt_strings`` is an optional iterable of genotype strings already seen
    by the caller (e.g. from the dict-phasing scan); when None we sample a few
    records off ``gt_reader`` to sniff a separator. Never raises.
    """
    if explicit_phased is not None:
        return bool(explicit_phased)

    # (2) manifest flag, if the build recorded one.
    manifest = getattr(gt_reader, "manifest", None)
    if isinstance(manifest, dict) and "phased" in manifest:
        try:
            return bool(manifest["phased"])
        except Exception:
            pass

    # (3) detect from observed genotype strings.
    saw_pipe = False
    saw_slash = False

    def _scan(strings):
        nonlocal saw_pipe, saw_slash
        for gt in strings:
            if gt is None:
                continue
            if "|" in gt:
                saw_pipe = True
            elif "/" in gt:
                saw_slash = True

    if observed_gt_strings is not None:
        _scan(observed_gt_strings)

    if not saw_pipe and gt_reader is not None:
        # Sniff a handful of records straight off the tier's GT vocabulary (the
        # distinct genotype strings, stored once) -- cheap and representative.
        try:
            _scan(gt_reader.gt_vocab())
        except Exception:
            pass

    if saw_pipe:
        return True
    if saw_slash:
        return False
    # Nothing observed -> conservative default: unphased (reports bounds, never an
    # over-confident cis allele frequency).
    return False


def _fmt_af(value, defined):
    """Format an allele frequency: '%.6g', or 'NA' when undefined (unphased combo)."""
    if not defined or value is None:
        return "NA"
    try:
        import math
        if isinstance(value, float) and math.isnan(value):
            return "NA"
    except Exception:
        pass
    return "%.6g" % value


def _fmt_freq(value):
    if value is None:
        return "NA"
    return "%.6g" % value


def _encode_group_breakdowns(summary, global_group_id, sep):
    """Return (allele_freq_by_group, carrier_freq_by_group, carrier_n_by_group).

    Each is a compact ';'-joined "group=value" string over the (db x subpop) AND db
    groups the summary reports (SPARSE: only groups with >=1 carrier), in sorted
    group-id order (deterministic). GLOBAL is omitted here (it has dedicated
    columns). Allele freq is 'NA' for a group whose AF is undefined (unphased
    combination).
    """
    af_parts = []
    cf_parts = []
    n_parts = []
    for gid in sorted(summary.groups):
        if gid == global_group_id:
            continue
        gs = summary.groups[gid]
        af_parts.append("%s=%s" % (gid, _fmt_af(gs.allele_freq, gs.allele_freq_defined)))
        cf_parts.append("%s=%s" % (gid, _fmt_freq(gs.carrier_freq)))
        n_parts.append("%s=%d" % (gid, gs.n_carrier))
    return (";".join(af_parts), ";".join(cf_parts), ";".join(n_parts))


def summarize_off_target(off_target, tier0_reader, gt_reader, axis, ploidy_of,
                         phased, population_summary, panel=None,
                         global_group_id="global", sep="::"):
    """Compute the companion ROW dict for ONE variant off-target, or None to SKIP.

    ``off_target`` is a mapping/record with at least:
      Chromosome, Position, Direction, crRNA, DNA (identity columns), and SNP (the
      comma-joined creating snp_info tokens). It may be a dict OR any object exposing
      those via ``__getitem__``; we read with ``_get`` below.

    Returns an ordered dict aligned to ``COMPANION_HEADER``, or None if the row has
    NO creating variant (reference / PAM-only -> skipped, per spec). The frequency
    correctness is DELEGATED to ``population_summary.summarize``; we only format.
    """
    chrom = _get(off_target, "Chromosome")
    snp_field = _get(off_target, "SNP")
    pos_alts = parse_snp_field(snp_field, chrom)
    if not pos_alts:
        return None  # reference / PAM-only off-target -> skip (no companion row)

    summary = population_summary.summarize(
        pos_alts, tier0_reader, gt_reader, axis, ploidy_of, phased, panel=panel)

    af_by, cf_by, n_by = _encode_group_breakdowns(summary, global_group_id, sep)

    return {
        "Chromosome": _get(off_target, "Chromosome"),
        "Position": _get(off_target, "Position"),
        "Direction": _get(off_target, "Direction"),
        "crRNA": _get(off_target, "crRNA"),
        "DNA": _get(off_target, "DNA"),
        "SNP": "" if snp_field is None else str(snp_field),
        "n_variants": str(summary.k),
        "global_allele_freq": _fmt_af(summary.global_af, summary.allele_freq_defined),
        "global_carrier_freq": _fmt_freq(summary.global_carrier_freq),
        "global_carrier_n": str(summary.global_carrier_n),
        "global_hom_freq": _fmt_freq(summary.global_hom_freq),
        "max_subpop_af": _fmt_af(summary.max_subpop_af, True),
        "max_subpop_label": "" if summary.max_subpop_af_label is None
                            else str(summary.max_subpop_af_label),
        "observed": "1" if summary.observed else "0",
        "allele_freq_defined": "1" if summary.allele_freq_defined else "0",
        "phased": "1" if phased else "0",
        "global_carrier_freq_lower_bound": _global_carrier_lower(summary),
        "allele_freq_by_group": af_by,
        "carrier_freq_by_group": cf_by,
        "carrier_n_by_group": n_by,
    }


def _global_carrier_lower(summary):
    """The GLOBAL carrier-freq LOWER BOUND string (unphased combination only).

    Labeled in the header as 'lower bound on the assume-cis carrier fraction'. Empty
    for the exact path (single-variant or phased), where n_carrier is exact.
    """
    g = summary.global_summary
    if g is None:
        return ""
    lower = g.carrier_freq_lower
    if lower is None:
        return ""
    return "%.6g" % lower


def _get(record, key):
    """Read ``key`` from a dict-like or sequence-mapping record; '' if absent."""
    try:
        val = record[key]
    except (KeyError, IndexError, TypeError):
        return ""
    return "" if val is None else val


def build_rows(off_targets, tier0_reader, gt_reader, axis, ploidy_of, phased,
               population_summary, panel=None, global_group_id="global",
               sep="::", on_error=None):
    """Build companion rows for a list of off-targets (PURE; no file I/O).

    Skips reference/PAM-only rows (no creating variant). Per-row errors are isolated:
    if ``summarize`` raises on one off-target we skip THAT row (optionally reporting
    via ``on_error``) and keep going, so one bad row can never sink the file. Returns
    the list of row dicts (aligned to ``COMPANION_HEADER``).
    """
    rows = []
    for ot in off_targets:
        try:
            row = summarize_off_target(
                ot, tier0_reader, gt_reader, axis, ploidy_of, phased,
                population_summary, panel=panel, global_group_id=global_group_id,
                sep=sep)
        except Exception as err:  # isolate a single bad row; never sink the file
            if on_error is not None:
                on_error(ot, err)
            continue
        if row is not None:
            rows.append(row)
    return rows


def write_companion(out_path, off_targets, tier0_reader, gt_reader, axis,
                    ploidy_of, population_summary, panel_cls=None, panel=None,
                    phased=None, global_group_id="global", sep="::",
                    observed_gt_strings=None, on_error=None):
    """Write the companion population-summary TSV. ADDITIVE, gated on tier0_reader.

    Args:
      out_path: the companion file path (e.g. "<output>.population_summary.tsv").
      off_targets: an iterable of variant off-target records (dict-like with the
        identity + SNP columns). Reference/PAM-only rows are skipped internally.
      tier0_reader: the Tier-0 ``RegistryReader`` (``myreg``). GATE: if this is None
        we write NOTHING (no file created) and return False -- BYTE-IDENTICAL legacy
        behavior. This is the single hard gate the spec requires.
      gt_reader: the Tier-1 ``GenotypeReader`` (``mygt``), or None. When None only
        SINGLE-variant summaries are exact (they come straight from the registry);
        a multi-variant row cannot be computed without the genotype tier, so its
        companion row is marked 'requires genotype tier' (SNP kept, freqs NA) rather
        than dropped -- the caller still sees the row exists.
      axis: the sample axis (``gt_reader``'s axis, or a registry-derived axis when
        gt_reader is None -- caller supplies whatever it has).
      ploidy_of: callable(sample_id, sex) -> ploidy (from tier0_compile).
      population_summary: the ``population_summary`` module (injected so this module
        has no hard import dependency -- an old deploy without it simply never calls
        here).
      panel_cls / panel: build the Panel ONCE per chromosome and reuse it. Pass
        ``panel_cls=population_summary.Panel`` (or a prebuilt ``panel``); we build it
        from ``axis`` + ``ploidy_of`` when possible so every off-target shares one.
      phased: explicit dataset-wide phasing bool, or None to resolve (manifest ->
        detect -> conservative default False). See ``resolve_phased``.
      observed_gt_strings: optional genotype strings the caller already saw (dict
        phasing scan) to help resolve ``phased`` without re-sniffing.
      on_error: optional callable(off_target, exception) for per-row error logging.

    Returns True iff the companion file was written (registry present), else False.
    Any exception here propagates to the caller, which MUST guard the call (log +
    skip) so the companion write can never break the run.
    """
    # ---- the single hard GATE: no registry -> write NOTHING (legacy byte-identical).
    if tier0_reader is None:
        return False

    resolved_phased = resolve_phased(phased, gt_reader, observed_gt_strings)

    # Build the Panel ONCE per chromosome and reuse across every off-target.
    if panel is None and panel_cls is not None and axis is not None \
            and ploidy_of is not None:
        try:
            panel = panel_cls(axis, ploidy_of)
        except Exception:
            panel = None  # summarize() will build per-call as a fallback

    rows = build_rows(
        off_targets, tier0_reader, gt_reader, axis, ploidy_of, resolved_phased,
        population_summary, panel=panel, global_group_id=global_group_id, sep=sep,
        on_error=on_error)

    with open(out_path, "w") as fh:
        fh.write("#" + "\t".join(COMPANION_HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(col, "")) for col in COMPANION_HEADER)
                     + "\n")
    return True
