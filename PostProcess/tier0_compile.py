"""Compile a legacy SNP dict into a Tier-0 panel-aware registry (Phase 1 step 2).

This is the build-time bridge from the OLD ``my_dict_<chrom>.json[.gz]`` (carriers-
only per-ALT entries) to the compact, mmap-friendly registry in
``tier0_registry`` -- with FULL-PANEL allele numbers so AF = AC / AN uses the
whole panel as its denominator (the dict lists only carriers; every unlisted
panel sample is a hom-ref).

Ground truth for the dict format (PostProcess/new_simple_analysis.py 158-199 and
261-290):

  my_dict_<chrom>.json[.gz] :  "<chrom>,<pos1based>" -> '$'-joined per-ALT entries.
    each entry :  "<samples>;<ref,alt>;<rsID>;<AF>"
      <samples> = comma-joined "sampleID:genotype"  (genotype "a0|a1"/"a0/a1"/"a0")
                  ONLY alt-carriers are listed.
      <ref,alt> = comma-joined reference and this entry's alt allele.
      <rsID>    = variant id (".", "rsNNN", ...).
      <AF>      = precomputed allele frequency string (we IGNORE it; we recompute
                  a correct full-panel AF from the counts).

  samplesID file (per database) : tab-separated columns
      #SAMPLE_ID   POPULATION_ID   SUPERPOPULATION_ID   SEX      (SEX in {male,female})
    with a leading '#'-header line (skipped). We use SUPERPOPULATION_ID as the
    subpopulation grouping for Phase 1 (matches PopulationDistribution). The
    finer POPULATION_ID is a documented future option (see ``subpop_field``).

STDLIB ONLY (json, gzip; optional ijson for streaming) so it runs in the
lightweight unit-tests CI (no numpy / no pysam).
"""

from __future__ import annotations

import gzip
import json

import tier0_registry as t0

# Optional streaming JSON parser. If present we stream the (huge, genome-wide)
# dict key-by-key instead of json.load-ing the whole thing into RAM.
try:  # pragma: no cover - availability depends on the environment
    import ijson  # type: ignore
    _HAVE_IJSON = True
except Exception:  # pragma: no cover
    ijson = None
    _HAVE_IJSON = False


# --------------------------------------------------------------------------- #
# samplesID parsing -> sample_meta
# --------------------------------------------------------------------------- #
def read_samplesid(path, database, subpop_field="superpopulation"):
    """Read one per-database samplesID file into sample_meta rows.

    Returns dict sample_id -> (database, subpopulation, sex). ``subpop_field``
    selects the grouping column: "superpopulation" (col 2, the Phase-1 default,
    matching PopulationDistribution) or "population" (col 1, a finer future
    option). The leading '#'-header line is skipped.
    """
    if subpop_field == "superpopulation":
        col = 2
    elif subpop_field == "population":
        col = 1
    else:
        raise ValueError("subpop_field must be 'superpopulation' or 'population'")

    meta = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = True
        for line in fh:
            line = line.rstrip("\n")
            if first:
                first = False
                # Skip the header if it looks like one; otherwise treat as data.
                if line.startswith("#") or line.upper().startswith("SAMPLE_ID") \
                        or "SAMPLE_ID" in line.upper():
                    continue
            if not line.strip():
                continue
            fields = line.split("\t")
            if not fields or not fields[0].strip():
                continue
            sid = fields[0].strip()
            subpop = fields[col].strip() if len(fields) > col else ""
            sex = fields[3].strip() if len(fields) > 3 else ""
            meta[sid] = (database, subpop, sex)
    return meta


def build_sample_meta(db_to_samplesid, subpop_field="superpopulation",
                      genotyped_samples=None):
    """Read the ordered {database: samplesID_path} mapping into one sample_meta.

    Cross-db shared sample ids: for the current disjoint 1000G/HGDP panels this
    does not occur. If a sample_id DID appear in two databases, this dict is 1:1
    (sample_id -> single (db, subpop, sex)) so the LAST database read would win
    the per-db membership. To keep BOTH per-db memberships (and let GLOBAL dedup
    the canonical id once), the caller should instead disambiguate the sample ids
    per database upstream; the GLOBAL group in ``tier0_registry`` already dedups by
    the canonical (bare) sample_id. Phase-1 panels are disjoint so we keep the
    simple 1:1 mapping and note this as a documented limitation.
    """
    sample_meta = {}
    overlaps = []
    for database, path in db_to_samplesid.items():
        rows = read_samplesid(path, database, subpop_field=subpop_field)
        for sid, row in rows.items():
            if sid in sample_meta and sample_meta[sid][0] != database:
                overlaps.append(sid)
            sample_meta[sid] = row
    if genotyped_samples is not None:
        # #46 auto-fix: drop phantom samples listed in the samplesID file but
        # genotyped in NO VCF -- otherwise they are counted as hom-ref and inflate
        # the panel AN denominator, silently deflating every reported allele
        # frequency. Filtering here makes AN = ploidy x genotyped, not x listed.
        # No-op (byte-identical panel) when genotyped_samples is None or already a
        # superset -- the common matched-panel + batteries (VCF-absent) cases.
        gset = set(genotyped_samples)
        before = len(sample_meta)
        sample_meta = {s: r for s, r in sample_meta.items() if s in gset}
        dropped = before - len(sample_meta)
        if dropped:
            print("tier0: dropped %d phantom samplesID sample(s) not genotyped in "
                  "any VCF; panel AN now over the genotyped panel (#46)" % dropped,
                  flush=True)
    return sample_meta, overlaps


# --------------------------------------------------------------------------- #
# Chromosome -> ploidy_of selection
# --------------------------------------------------------------------------- #
def _norm_chrom(chrom):
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c.upper()


def ploidy_of_for_chrom(chrom):
    """Pick the ploidy_of callable for a chromosome.

    * autosomes (and PAR, TODO) -> ``autosomal_ploidy`` (everyone diploid).
    * chrX-nonPAR -> males haploid, females diploid.
    * chrY -> males haploid, females absent (caller omits chrY females; a female
      row simply never carries a chrY alt, and its baseline ploidy is set below).

    PAR (pseudoautosomal) handling is a documented TODO: within PAR, males are
    diploid on chrX/chrY. Phase 1 treats all of chrX-nonPAR / chrY as haploid for
    males; positions inside PAR would need a coordinate-aware override.
    """
    c = _norm_chrom(chrom)
    if c == "X":
        # chrX-nonPAR: haploid males, diploid females.
        return t0.make_chr_ploidy(haploid_male=True, haploid_female=False)
    if c == "Y":
        # chrY: haploid males; females carry NO chrY -> ABSENT (ploidy 0), so their
        # samplesID rows add no phantom alleles to the chrY AN denominator.
        return t0.make_chr_ploidy(haploid_male=True, absent_female=True)
    # Autosomes (and, TODO, PAR): diploid everyone.
    return t0.autosomal_ploidy


# --------------------------------------------------------------------------- #
# Dict entry parsing
# --------------------------------------------------------------------------- #
def parse_entry(entry):
    """Parse one '$'-delimited dict entry into (ref, alt, rsid, carrier_gts).

    entry = "<samples>;<ref,alt>;<rsID>;<AF>". ``carrier_gts`` maps sample_id ->
    genotype string (ONLY the listed alt-carriers). Returns None if the entry is
    malformed (too few ';' fields).
    """
    parts = entry.split(";")
    if len(parts) < 3:
        return None
    samples_field = parts[0].strip()
    refalt = parts[1].strip().split(",")
    if len(refalt) < 2:
        return None
    ref, alt = refalt[0], refalt[1]
    rsid = parts[2].strip() if len(parts) > 2 else "."

    carrier_gts = {}
    if samples_field:
        for tok in samples_field.split(","):
            tok = tok.strip()
            if not tok:
                continue
            # "sampleID:genotype"; genotype may itself be absent (defensive).
            if ":" in tok:
                sid, gt = tok.split(":", 1)
            else:
                sid, gt = tok, ""
            carrier_gts[sid.strip()] = gt.strip()
    return ref, alt, rsid, carrier_gts


def iter_dict_records(dict_path, chrom):
    """Stream (pos:int, ref, alt, rsid, carrier_gts) records from the SNP dict.

    Uses ijson if available (constant-memory streaming of the top-level object),
    else falls back to json.load. Handles gzip transparently by extension. Only
    keys whose "<chrom>," prefix matches ``chrom`` are yielded (the dict is
    per-chromosome, but we filter defensively so a mislabeled file cannot leak
    other chromosomes' positions into this chromosome's registry).
    """
    chrom = str(chrom)
    prefix = chrom + ","
    is_gz = str(dict_path).endswith(".gz")

    if _HAVE_IJSON:
        opener = gzip.open if is_gz else open
        with opener(dict_path, "rb") as fh:
            # kvitems over the root object yields (key, value) without building
            # the whole dict in memory.
            for key, value in ijson.kvitems(fh, ""):
                for rec in _emit_key(key, value, prefix):
                    yield rec
    else:
        opener = gzip.open if is_gz else open
        with opener(dict_path, "rt") as fh:
            data = json.load(fh)
        for key, value in data.items():
            for rec in _emit_key(key, value, prefix):
                yield rec


def _emit_key(key, value, prefix):
    """Yield parsed records for one dict key/value ('$'-joined entries)."""
    if not key.startswith(prefix):
        return
    try:
        pos = int(key.split(",", 1)[1])
    except (ValueError, IndexError):
        return
    if not value:
        return
    for entry in str(value).split("$"):
        parsed = parse_entry(entry)
        if parsed is None:
            continue
        ref, alt, rsid, carrier_gts = parsed
        yield pos, ref, alt, rsid, carrier_gts


# --------------------------------------------------------------------------- #
# Top-level compile
# --------------------------------------------------------------------------- #
def compile_from_dict(dict_path, db_to_samplesid, chrom, out_bin, out_idx,
                      *, alt_index="1", subpop_field="superpopulation",
                      compress=False, genotyped_samples=None):
    """Compile a legacy SNP dict into a panel-aware Tier-0 registry.

    Args:
      dict_path: path to ``my_dict_<chrom>.json`` or ``.json.gz``.
      db_to_samplesid: ORDERED mapping {database_name: samplesID_path}, e.g.
        {"1000G": ".../hg38_1000G.samplesID.txt",
         "HGDP":  ".../hg38_HGDP.samplesID.txt"}.
      chrom: the chromosome name as written in the dict keys (e.g. "chr1",
        "chrX"). Selects the ploidy model (autosomal vs chrX/Y haploid-male).
      out_bin, out_idx: output registry binary + json manifest paths.
      alt_index: the genotype token denoting each record's alt (default "1", the
        bcftools norm -m- biallelic convention). Passed through to the aggregator.
      subpop_field: "superpopulation" (Phase-1 default) or "population".

    SNP-first: Phase 1 supports single-character alts only. Non-SNP alts (multi-
    char, i.e. indels) are SKIPPED with a counter (never crash). Records whose
    carrier set is empty after parsing (no listed carriers) produce no group and
    are skipped too (they carry no off-target signal).

    Returns a dict of build stats:
      {"manifest": <manifest dict>, "n_written": int, "n_skipped_indel": int,
       "n_skipped_empty": int, "n_positions": int, "overlaps": [sample_id, ...]}.
    """
    sample_meta, overlaps = build_sample_meta(db_to_samplesid,
                                              subpop_field=subpop_field,
                                              genotyped_samples=genotyped_samples)
    ploidy_of = ploidy_of_for_chrom(chrom)

    # Build the PanelIndex ONCE (per-group hom-ref baselines) and reuse it.
    panel_index = t0.PanelIndex(sample_meta, ploidy_of)

    stats = {
        "n_written": 0,
        "n_skipped_indel": 0,
        "n_skipped_empty": 0,
        "n_positions": 0,
    }
    seen_positions = set()

    def record_stream():
        for (pos, ref, alt, rsid, carrier_gts) in iter_dict_records(dict_path, chrom):
            seen_positions.add(pos)
            # SNP-first: a real single-base substitution only. Skip indels
            # (multi-char alt or a multi-char ref) AND non-ACGT placeholders such
            # as the '*' spanning/overlapping-deletion marker, never crash.
            if not alt or len(alt) != 1 or (ref and len(ref) != 1) \
                    or alt.upper() not in ("A", "C", "G", "T"):
                stats["n_skipped_indel"] += 1
                continue
            if not carrier_gts:
                # No listed carriers -> no group -> nothing to write.
                stats["n_skipped_empty"] += 1
                continue
            stats["n_written"] += 1
            yield (pos, ref, alt, rsid, carrier_gts)

    manifest = t0.compile_registry_panel(
        record_stream(), sample_meta, None, ploidy_of, out_bin, out_idx,
        alt_index=alt_index, panel_index=panel_index, compress=compress,
    )

    stats["n_positions"] = len(seen_positions)
    stats["manifest"] = manifest
    stats["overlaps"] = overlaps
    return stats
