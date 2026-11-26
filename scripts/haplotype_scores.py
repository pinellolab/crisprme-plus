#!/usr/bin/env python
"""``haplotype_scores.py`` -- a row per observed haplotype, each with a score.

THE QUESTION THIS ANSWERS. "Is this sequence safe to test?" Every per-site product
CRISPRme+ writes collapses a SET of observed haplotypes to ONE representative, and each
picks it by a DIFFERENT key:

  * ``integrated_results`` / ``top1000``  -- by score, per axis (the displayed alignment
    is the highest-CFD one; the CRISPR-Bulge number beside it comes from the
    highest-CRISPR-Bulge one);
  * ``observed_haplotypes.tsv``           -- by variant set, at most two observations per
    source row, only at sites carrying >=2 distinct sets;
  * ``extract-flanks``                    -- by CARRIER COUNT (``(not found, -carriers,
    seq)``), a key that reads no score at all.

The one thing none of them carries is A SCORE ATTACHED TO A HAPLOTYPE. That is what this
writes: one row per site x haplotype x frame, so "everything above a threshold" and
"everything also in the top-1000" are filters the reader writes, not new columns we invent.

WHY A COMPANION FILE, NEVER A NEW COLUMN ON THE MERGED ROW. The finalized ``final_line``
is addressed by NEGATIVE index in the CFD and CRISPR-Bulge scorers
(``observed_haplotypes_companion.py``), so appending a trailing column silently corrupts
every variant row's score. Nothing here touches ``final_line``.

TWO TIERS.

  Tier 1, HARVEST (the default).  In per-sample mode the pipeline already finalizes and
  scores one row per observed haplotype (``new_simple_analysis.py``, ``for hap in
  haplotypes: _finalize_observed_entry(...)``). Reading those out costs ~0 and CANNOT
  disagree with the report, because they ARE the report's numbers.
  THE CAVEAT, SHIPPED AS A COLUMN: ``_finalize_observed_entry`` re-inserts the
  candidate's gap columns at fixed positions and counts mismatches positionally, so a
  harvested score is the haplotype's bases IN THE CANDIDATE'S ALIGNMENT FRAME, not the
  haplotype's own best alignment. Every harvested row says ``alignment_frame=candidate_frame``.

  Tier 2, ANCHORED RE-ALIGNMENT (opt-in, ``--score-source realign|auto``).  Only for rows
  the fixed frame cannot hold: an indel inside the scored footprint, ``target_found=no``,
  or no harvested score. ``scan_sequence`` (``scripts/check_indel_rows.py``) generates the
  in-budget candidates and we keep THE ONE MATCHING THE ROW'S PUBLISHED ALIGNMENT --
  never an argmax. A 32-base flank admits a median of ~190 in-budget alignments, and
  argmax CFD reproduces the published reference CFD on 4 of 186 rows while OVERSTATING it
  on 74. Containment does not license selection. Such rows say
  ``alignment_frame=haplotype_own``.

  Tier 2 is built but DEFAULT OFF. Measured on the delivering run: of 46 distinct indel
  variants among a 187-site shortlist, only 9 have their own independently scored row.
  That does not justify ``auto``.

TWO INVARIANTS.

  1. NEVER A BLANK. ``score_status`` is always set. Today a never-scored cell and the
     CRISPR-Bulge out-of-domain sentinel (-1.0) both render "NA", so a reader filtering
     ``CRISPR_Bulge >= 1.9e-4`` silently discards them and concludes the site is clean.
     That failure is what produced the question this file answers. -1.0 is never printed,
     never sorted on, never turned into a blank.
  2. COUNT BEFORE YOU CAP. Every per-site counter is computed on the full stream BEFORE
     any cap, so a cap may hide a sequence but never a number.

Stdlib only. Importable: the eventual home is ``extract-flanks --score``, so the logic
lives in functions and nothing runs at import.
"""

from __future__ import print_function

import argparse
import collections
import csv
import io
import os
import re
import sys
import zipfile

# --------------------------------------------------------------------------- #
# score_status -- the vocabulary. A cell is NEVER blank and NEVER -1.0.

ST_SCORED = "scored"
ST_OUT_OF_DOMAIN = "out-of-domain:%s"          # the model's own refusal rule
ST_NO_ALIGNMENT = "no-in-budget-alignment(%s)"  # budget quoted INSIDE the string
ST_NEEDS_REALIGN = "needs-realignment:%s"       # site_key; Tier 2 would resolve it
ST_UNRESOLVED = "unresolved-copy:%s"            # overlapping-records | unphased-genotype | ...
ST_WINDOW_DELETED = "window-deleted"
ST_NOT_CHECKED = "not-checked:no-vcf-on-contig"
ST_CAPPED = "capped-out"

FRAME_CANDIDATE = "candidate_frame"
FRAME_OWN = "haplotype_own"

# The CRISPR-Bulge null sentinel. It must never reach a column or a sort.
CB_NULL = -1.0

COLUMNS = [
    "site_key", "chrom", "site_start", "strand", "spacer",
    "haplotype_id", "variant_set",
    "aligned_spacer_pam", "aligned_protospacer_pam",   # make every number re-verifiable
    "alignment_frame",
    "cfd", "crispr_bulge", "score_status",
    "above_threshold", "in_top1000_cfd", "in_top1000_cb",
    "rank_in_site_by_score", "score_vs_site_reported",
    "n_carriers", "carriers",
    "oligo_seq", "oligo_len", "oligo_status",
    "oligo_group_id", "oligo_group_size", "oligo_group_carriers",
    "variants_in_window_n", "max_af",
    "scope_reason",
    # per-site roll-up, repeated on each of the site's rows so one row is self-contained
    "n_haplotypes_at_site", "n_above_threshold", "n_not_scored",
    "site_max_crispr_bulge", "row_reported_alignment_found",
]


class HaplotypeScoreError(Exception):
    """Refuse loudly. Never degrade into a complete-looking, unscored file."""


# --------------------------------------------------------------------------- #
# formatting


def fmt_score(value):
    """``%.6g``, never 3 decimal places.

    The sensitive operating point is 1.9e-4, BELOW the old 5e-4 display floor: at 3 dp it
    rounds to 0.000 and the threshold is destroyed in the file meant to carry it. A null
    or missing score is never a number here -- the caller writes MISSING and says why in
    ``score_status``.
    """
    if value is None:
        return "NA"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "NA"
    if f < 0:            # the -1.0 sentinel: NEVER printed as a number
        return "NA"
    return "%.6g" % f


# --------------------------------------------------------------------------- #
# the run's own budget -- read, never a flag


BUDGET_KEYS = (
    ("mismatches", "Mismatches"),
    ("dna", "DNA"),
    ("rna", "RNA"),
    ("max_total", "Max_total_edits"),
    ("pam", "Pam"),
    ("scorer", "Scorer"),
)


def read_budget(results_dir):
    """The search's OWN budget, from its ``.Params.txt``.

    Deliberately NOT a flag: a budget wider than the search's produces alignments, and so
    numbers, the report could never have produced. Unreadable means refuse.
    """
    path = os.path.join(results_dir, ".Params.txt")
    if not os.path.isfile(path):
        raise HaplotypeScoreError(
            "no .Params.txt in %s: the search's own budget (mismatches, bulges, PAM) is "
            "read from the run, never passed as a flag, because a different budget "
            "produces numbers the report could not have produced." % results_dir)
    raw = {}
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                raw[parts[0].strip()] = parts[1].strip()
    out = {}
    missing = []
    for key, name in BUDGET_KEYS:
        if name not in raw:
            missing.append(name)
        else:
            out[key] = raw[name]
    if missing:
        raise HaplotypeScoreError(
            "%s does not record %s: cannot establish the search's budget, so refusing "
            "rather than guessing one." % (path, ", ".join(missing)))
    for k in ("mismatches", "dna", "rna", "max_total"):
        try:
            out[k] = int(out[k])
        except ValueError:
            raise HaplotypeScoreError(
                "%s records %s as %r, which is not a number." % (path, k, out[k]))
    return out


def budget_text(budget):
    """``6mm/1DNA/1RNA/8total`` -- quoted INSIDE score_status, so a row that found no
    alignment says which budget failed to find one."""
    return "%dmm/%dDNA/%dRNA/%dtotal" % (
        budget["mismatches"], budget["dna"], budget["rna"], budget["max_total"])


# --------------------------------------------------------------------------- #
# the input source -- probe, record, or fail loudly


RAW_SIDECAR_RE = re.compile(r"_observed_haplotypes_raw\.(tsv|zip)$")
ALT_TABLE_RE = re.compile(r"_all_results_with_alternative_alignments\.(tsv|zip)$")


def probe_source(results_dir):
    """``(path, kind)`` for the richest haplotype source present.

    Order matters. ``*_observed_haplotypes_raw`` is written at search time and filters
    ONLY on origin -- an explicit UNSCOPED superset of the slim companion, which is
    additionally restricted to within 50 bp of a CFD-top-1000 site. Falling back to the
    full alt table is correct but larger and slower.

    Neither present is NOT "nothing found": it means this run cannot answer the question,
    and saying so loudly is the whole point.
    """
    raw, alt = [], []
    for name in sorted(os.listdir(results_dir)):
        full = os.path.join(results_dir, name)
        if not os.path.isfile(full):
            continue
        if RAW_SIDECAR_RE.search(name):
            raw.append(full)
        elif ALT_TABLE_RE.search(name):
            alt.append(full)
    if raw:
        return raw[0], "observed_haplotypes_raw"
    if alt:
        return alt[0], "all_results_with_alternative_alignments"
    raise HaplotypeScoreError(
        "neither *_observed_haplotypes_raw.(tsv|zip) nor "
        "*_all_results_with_alternative_alignments.(tsv|zip) is in %s. Per-haplotype "
        "enumeration is gated on a Tier-1 dict-less genotype store and is ABSENT "
        "ENTIRELY in population-level fast mode, so an empty companion here would read "
        "as 'no extra high-scoring haplotypes' when the truth is 'this run never "
        "enumerated any'. Refusing to write." % results_dir)


def open_table(path):
    """Stream a .tsv, or the single member of a .zip, as text lines."""
    if path.endswith(".zip"):
        zf = zipfile.ZipFile(path)
        names = [n for n in zf.namelist() if n.endswith(".tsv")]
        if not names:
            raise HaplotypeScoreError("%s holds no .tsv member" % path)
        handle = zf.open(names[0])
        return io.TextIOWrapper(handle, encoding="utf-8", errors="replace")
    return io.open(path, "r", encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------- #
# the three alignment blocks


BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def block_columns(header):
    """``{block: {logical: index}}`` for whichever of the three blocks the table carries.

    The raw per-haplotype tables carry all three; a curated table carries ONE alignment
    with both scores beside it. Both shapes are read.
    """
    idx = {name: i for i, name in enumerate(header)}
    out = {}
    for block in BLOCKS:
        cols = {}
        for logical, template in (
                ("start", "Start_coordinate_(%s)"),
                ("strand", "Strand_(%s)"),
                ("aligned_spacer", "Aligned_spacer+PAM_(%s)"),
                ("aligned_alt", "Aligned_protospacer+PAM_ALT_(%s)"),
                ("aligned_ref", "Aligned_protospacer+PAM_REF_(%s)"),
                ("variants", "Variant_info_genome_(%s)"),
                ("maf", "Variant_MAF_(%s)"),
                ("samples", "Variant_samples_(%s)"),
                ("origin", "REF/ALT_origin_(%s)"),
        ):
            name = template % block
            if name in idx:
                cols[logical] = idx[name]
        for logical, name in (("cfd", "CFD_score_(%s)" % block),
                              ("crispr_bulge", "CRISPR_BULGE_score_(%s)" % block)):
            if name in idx:
                cols[logical] = idx[name]
        if "start" in cols:
            out[block] = cols
    return out, idx


def site_key(chrom, start, strand, spacer):
    return "%s:%s:%s:%s" % (chrom, start, strand, spacer)


# --------------------------------------------------------------------------- #
# Tier 1 -- harvest


Observation = collections.namedtuple(
    "Observation",
    "site chrom start strand spacer block variants aligned_spacer aligned_alt "
    "cfd crispr_bulge samples maf origin")


def harvest_rows(stream, spacer_hint=None):
    """Yield one ``Observation`` per (row, block) the table actually names.

    Takes the score ONLY from the axis that named the haplotype. The integrated
    projections carry CFD on the CFD axis and CRISPR-Bulge on the CB axis with no
    cross-axis score, which is exactly why most rows of the delivered companion are NA on
    one of them. Promising both from this source would be inventing a number.
    """
    reader = csv.reader(stream, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration:
        return
    blocks, idx = block_columns(header)
    if not blocks:
        raise HaplotypeScoreError(
            "this table names none of the three alignment blocks (%s) -- it is not a "
            "per-haplotype results table." % ", ".join(BLOCKS))
    chrom_i = idx.get("Chromosome")
    spacer_i = idx.get("Spacer+PAM")

    def cell(row, cols, logical):
        i = cols.get(logical)
        if i is None or i >= len(row):
            return ""
        return row[i].strip()

    for row in reader:
        if not row or len(row) < 3:
            continue
        chrom = row[chrom_i].strip() if chrom_i is not None and chrom_i < len(row) else ""
        spacer = (row[spacer_i].strip() if spacer_i is not None and spacer_i < len(row)
                  else (spacer_hint or ""))
        for block, cols in blocks.items():
            start = cell(row, cols, "start")
            if not start:
                continue
            strand = cell(row, cols, "strand")
            variants = cell(row, cols, "variants")
            origin = cell(row, cols, "origin")
            yield Observation(
                site=site_key(chrom, start, strand, spacer),
                chrom=chrom, start=start, strand=strand, spacer=spacer, block=block,
                variants=variants,
                aligned_spacer=cell(row, cols, "aligned_spacer"),
                aligned_alt=cell(row, cols, "aligned_alt") or cell(row, cols, "aligned_ref"),
                cfd=cell(row, cols, "cfd"),
                crispr_bulge=cell(row, cols, "crispr_bulge"),
                samples=cell(row, cols, "samples"),
                maf=cell(row, cols, "maf"),
                origin=origin)


# --------------------------------------------------------------------------- #
# scoring status


def classify(obs, budget, axis_named):
    """``(cfd, crispr_bulge, score_status)`` for one harvested observation.

    ``axis_named`` says which axis named this haplotype. The other axis gets no number
    and says so -- never a blank, never a borrowed value.
    """
    def num(text):
        if text in ("", "NA", "n", "-"):
            return None
        try:
            v = float(text)
        except ValueError:
            return None
        if v <= CB_NULL:      # the model's null sentinel
            return None
        return v

    cfd = num(obs.cfd)
    cb = num(obs.crispr_bulge)
    if cfd is None and cb is None:
        return None, None, ST_NEEDS_REALIGN % obs.site
    if cb is None and axis_named == "highest_CRISPR_BULGE":
        # the scorer refused this alignment by its own rule
        return cfd, None, ST_OUT_OF_DOMAIN % "multi-gap-or-length"
    return cfd, cb, ST_SCORED


# --------------------------------------------------------------------------- #
# per-site accumulation -- counters BEFORE any cap


class Site(object):
    """One off-target site's haplotypes, with counters computed on the full stream.

    The cap is applied when rows are WRITTEN, never when they are COUNTED, so a capped
    file still reports the true totals. A cap may hide a sequence; it must never hide a
    number.
    """

    __slots__ = ("key", "chrom", "start", "strand", "spacer", "rows",
                 "haplotypes", "above", "scored_any", "max_cb",
                 "reported_alignments", "seen_alignments", "applicable", "has_report")

    def __init__(self, key, chrom, start, strand, spacer):
        self.key = key
        self.chrom, self.start, self.strand, self.spacer = chrom, start, strand, spacer
        self.rows = []
        # Counted per HAPLOTYPE, not per row. One haplotype appears once per alignment
        # block it was named on, so counting rows would multiply every total by the
        # number of frames -- a site of 3 haplotypes would report 9.
        self.haplotypes = set()
        self.above = set()
        self.scored_any = set()
        self.max_cb = None
        self.reported_alignments = set()
        self.seen_alignments = set()
        self.applicable = False
        self.has_report = False

    def add(self, row, cb_value, above, scored, aligned, hap_key):
        self.haplotypes.add(hap_key)                # COUNT first, and by haplotype
        if above:
            self.above.add(hap_key)
        if scored:
            self.scored_any.add(hap_key)
        if cb_value is not None and (self.max_cb is None or cb_value > self.max_cb):
            self.max_cb = cb_value
        if aligned:
            self.seen_alignments.add(aligned)
        self.rows.append(row)                       # cap applies to THIS, later

    @property
    def n_haplotypes(self):
        return len(self.haplotypes)

    @property
    def n_above(self):
        return len(self.above)

    @property
    def n_not_scored(self):
        """Haplotypes with no score on ANY axis -- not rows missing one axis."""
        return len(self.haplotypes - self.scored_any)

    @property
    def row_reported_alignment_found(self):
        """Whether the site's OWN reported alignment appears on a haplotype we wrote.

        The only mechanism that distinguishes "these are all the haplotypes" from "these
        are the ones we could reconstruct". It fires in real data.
        """
        # Three distinct answers, never collapsed:
        #   unknown -- nothing told us what this site reported (no --targets)
        #   n/a     -- it reported a REFERENCE haplotype, so there is no ALT haplotype
        #              for it to appear among; the check does not apply
        #   yes/no  -- it reported an ALT haplotype, and it is or is not among them
        if not self.has_report:
            return "unknown"
        if not self.applicable:
            return "n/a"
        return "yes" if (self.reported_alignments & self.seen_alignments) else "no"


# --------------------------------------------------------------------------- #
# the oligo column -- a second geometry, not a post-hoc slice


def oligo_from(sequence, offset, length, pad_up, pad_down):
    """The short orderable sequence around the scored footprint.

    NOT called a 32-mer: under the anchored policy the delivered "32-mer" file is 32 bases
    on 859 of 1,258 lines and ranges 22..34, because the target's reference span moves
    with bulges and indels. ``oligo_len`` is written beside it for that reason.
    """
    if sequence is None or offset is None or length is None:
        return "", 0, "no-geometry"
    try:
        offset, length = int(offset), int(length)
    except (TypeError, ValueError):
        return "", 0, "no-geometry"
    if not sequence:
        return "", 0, ST_WINDOW_DELETED
    lo = max(0, offset - pad_up)
    hi = min(len(sequence), offset + length + pad_down)
    if lo >= hi:
        return "", 0, "window-edge-unsafe"
    piece = sequence[lo:hi]
    status = "ok"
    if lo == 0 and offset - pad_up < 0:
        status = "window-edge-unsafe"
    elif hi == len(sequence) and offset + length + pad_down > len(sequence):
        status = "window-edge-unsafe"
    return piece, len(piece), status


# --------------------------------------------------------------------------- #
# the header -- preconditions travel WITH the data


def file_header(source_path, source_kind, budget, threshold, recall, cap, strict):
    lines = [
        "# CRISPRme+ per-haplotype scores. One row per SITE x HAPLOTYPE x FRAME.",
        "#",
        "# WHAT THIS IS. Every other per-site product collapses the observed haplotypes to",
        "# one representative, each by a different key (score per axis; variant set; carrier",
        "# count). This attaches a score to each haplotype so 'everything above a threshold'",
        "# is a filter you write, not a column we guessed.",
        "#",
        "# score_source=%s (%s)" % (source_kind, os.path.basename(source_path)),
        "# budget=%s (read from the run's .Params.txt, never a flag: a wider budget" % budget_text(budget),
        "#   produces numbers the report could not have produced)",
        "# pam=%s scorer=%s" % (budget.get("pam", "?"), budget.get("scorer", "?")),
        "#",
        "# THRESHOLD. above_threshold is the CRISPR-Bulge sensitive operating point,",
        "#   %s at %d%% recall, taken BY NAME from cb_operating_point.py." % (
            fmt_score(threshold), recall),
        "#   Its guide-cluster bootstrap 95%% CI is [5.8e-06, 3.1e-03]: it is an operating",
        "#   point, NOT a safe/unsafe line. Scores are written %.6g, never 3 dp -- at 3 dp",
        "#   the threshold itself rounds to 0.000.",
        "#",
        "# NEVER A BLANK. score_status is always set. A never-scored cell and the",
        "#   CRISPR-Bulge out-of-domain sentinel both used to render 'NA', so filtering",
        "#   on the score silently discarded them and the site looked clean. -1.0 is never",
        "#   printed, never sorted on.",
        "#",
        "# COUNT BEFORE CAP. n_haplotypes_at_site / n_above_threshold / n_not_scored are",
        "#   computed on the full stream before any cap (cap=%s). A cap may hide a" % (
            cap if cap else "none"),
        "#   sequence; it never hides a number.",
        "#",
        "# PRECONDITIONS. Per-haplotype enumeration is gated on a Tier-1 dict-less genotype",
        "#   store and is ABSENT ENTIRELY in population-level fast mode. An empty file from",
        "#   such a run would mean 'never enumerated', not 'nothing above threshold' -- this",
        "#   command refuses to write rather than let that ambiguity ship.",
        "#",
        "# alignment_frame=candidate_frame means the score is the haplotype's bases in the",
        "#   CANDIDATE's alignment shape (the pipeline re-inserts gap columns at fixed",
        "#   positions and counts mismatches positionally). haplotype_own means the",
        "#   haplotype was re-aligned on its own, anchored to the row's published alignment.",
        "#",
        "# RE-VERIFY ANY ROW: aligned_spacer_pam + aligned_protospacer_pam feed",
        "#   'crisprme.py targets-scoring' directly.",
        "#",
        "# row_reported_alignment_found=no means the site's OWN reported alignment is on",
        "#   none of the haplotypes written: the list is what we could reconstruct, not",
        "#   necessarily all of them.%s" % (
            " --strict makes this a non-zero exit." if strict else
            " (warned, not fatal; --strict makes it a non-zero exit.)"),
    ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# driver


_RSID = re.compile(r"\([^)]*\)")


def canonical_variants(text):
    """A comparable key for a variant set, whatever notation wrote it.

    THE RECONCILIATION CANNOT USE THE ALIGNMENT STRING. The curated table's
    ``Aligned_protospacer+PAM`` is the REPRESENTATIVE's own alignment; the per-haplotype
    rows carry the haplotype's bases re-inserted into the CANDIDATE's alignment shape.
    They are different objects by construction, and measurement agrees: over the
    delivering run's ALT-origin sites the two strings are equal on 6 of 52, and no
    normalisation (case, gaps, either column, any of the three blocks) lifts it above 9.
    Matching on them would report ~90% of sites as unreconciled, which is noise.

    The variant SET is comparable -- the same records, written in two notations:
        curated  ``chr11;6482966;T;G (rs10136)``
        per-row  ``chr11_6482966_T_G``
    so both are reduced to a sorted tuple of ``chrom:pos:REF:ALT``.
    """
    if not text:
        return ()
    text = _RSID.sub("", text).replace(" ", "")
    out = []
    for piece in text.split(","):
        piece = piece.strip()
        if not piece or piece.upper() == "NA":
            continue
        parts = re.split(r"[;_]", piece)
        if len(parts) >= 4:
            out.append("%s:%s:%s:%s" % (parts[0], parts[1],
                                        parts[2].upper(), parts[3].upper()))
    return tuple(sorted(out))


Reported = collections.namedtuple("Reported", "variants origin")


def read_reported(targets_path):
    """``{(chrom, start): Reported}`` -- each site's own reported haplotype.

    Keyed on the VARIANT SET (see ``canonical_variants``), and carrying the row's
    REF/ALT origin: a REF-origin site's representative is the reference sequence, so it
    has no ALT haplotype to appear among and reconciliation does not apply to it. Over
    the delivering run that is 115 of 187 sites -- checking them would manufacture
    failures.
    """
    out = {}
    with io.open(targets_path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return out
        idx = {n: i for i, n in enumerate(header)}
        c_i, p_i = idx.get("Chromosome"), idx.get("Position")
        v_i, o_i = idx.get("Variant"), idx.get("REF/ALT_origin")
        if c_i is None or p_i is None or v_i is None:
            raise HaplotypeScoreError(
                "%s is not a curated results table: it needs Chromosome, Position and "
                "Variant to say which haplotype each site reported." % targets_path)
        for row in reader:
            if max(c_i, p_i, v_i) >= len(row):
                continue
            origin = (row[o_i].strip().lower()
                      if o_i is not None and o_i < len(row) else "")
            out[(row[c_i].strip(), row[p_i].strip())] = Reported(
                variants=canonical_variants(row[v_i]), origin=origin)
    return out


def build(results_dir, out_path, cap=0, score_source="harvest", strict=False,
          anonymise=False, pad_up=3, pad_down=6, log=None, targets=None):
    """Write the companion. Returns ``(n_rows, n_sites, unreconciled)``."""
    say = log or (lambda m: print(m, file=sys.stderr))
    budget = read_budget(results_dir)
    source_path, source_kind = probe_source(results_dir)
    reported = read_reported(targets) if targets else {}
    if targets:
        say("haplotype-scores: reported alignments for %d site(s) from %s"
            % (len(reported), os.path.basename(targets)))
    else:
        say("haplotype-scores: no --targets given, so row_reported_alignment_found "
            "is 'unknown' on every row (never 'yes')")
    say("haplotype-scores: source %s (%s)" % (source_kind, os.path.basename(source_path)))
    say("haplotype-scores: budget %s from the run's .Params.txt" % budget_text(budget))

    threshold, recall = _threshold()
    sites = collections.OrderedDict()
    wanted = set(reported) if reported else None
    if wanted is not None:
        say("haplotype-scores: restricted to the %d site(s) --targets names" % len(wanted))

    with open_table(source_path) as stream:
        for obs in harvest_rows(stream):
            if obs.origin and obs.origin.lower() == "ref":
                continue                       # reference rows are not variant haplotypes
            if wanted is not None and (obs.chrom, obs.start) not in wanted:
                # Bounds memory to the sites actually asked for. Without --targets every
                # site of the run is held at once, which a genome-wide sidecar (millions
                # of rows x three blocks) will not fit in. Streaming with cluster close-out
                # is the next step; see the note in the PR.
                continue
            site = sites.get(obs.site)
            if site is None:
                site = Site(obs.site, obs.chrom, obs.start, obs.strand, obs.spacer)
                rep = reported.get((obs.chrom, obs.start))
                site.has_report = rep is not None
                if rep is not None:
                    # REF-origin: the representative IS the reference sequence, so it has
                    # no ALT haplotype to appear among. Not applicable, not a failure.
                    site.applicable = (rep.origin == "alt") and bool(rep.variants)
                    site.reported_alignments = {rep.variants} if site.applicable else set()
                sites[obs.site] = site
            cfd, cb, status = classify(obs, budget, obs.block)
            above = cb is not None and cb >= threshold
            row = {
                "site_key": obs.site, "chrom": obs.chrom, "site_start": obs.start,
                "strand": obs.strand, "spacer": obs.spacer,
                "haplotype_id": "%s|%s" % (obs.block, obs.variants or "-"),
                "variant_set": obs.variants or "-",
                "aligned_spacer_pam": obs.aligned_spacer,
                "aligned_protospacer_pam": obs.aligned_alt,
                "alignment_frame": FRAME_CANDIDATE,
                "cfd": fmt_score(cfd), "crispr_bulge": fmt_score(cb),
                "score_status": status,
                "above_threshold": "yes" if above else "no",
                "in_top1000_cfd": "unknown", "in_top1000_cb": "unknown",
                "rank_in_site_by_score": "",
                "score_vs_site_reported": "",
                "n_carriers": str(len([s for s in obs.samples.split(",") if s.strip()])),
                "carriers": "" if anonymise else obs.samples,
                "oligo_seq": "", "oligo_len": "0", "oligo_status": "no-geometry",
                "oligo_group_id": "", "oligo_group_size": "", "oligo_group_carriers": "",
                "variants_in_window_n": str(len([v for v in obs.variants.split(",")
                                                 if v.strip()])),
                "max_af": obs.maf or "NA",
                "scope_reason": source_kind,
            }
            site.add(row, cb, above, status == ST_SCORED,
                     canonical_variants(obs.variants),
                     hap_key=obs.variants or obs.aligned_alt or "-")

    # A mis-keyed comparison looks exactly like a catastrophic finding. Before reporting
    # per-site failures, ask whether the reported alignments resemble the haplotype
    # alignments AT ALL: if almost nothing matches, the counterpart column is wrong, not
    # the data. Reporting 106 of 114 sites as unreconciled would be noise presented as
    # signal, and a reader cannot tell the difference from the file.
    # The heuristic needs enough evidence to fire. On a handful of sites "none matched"
    # is indistinguishable from one genuine failure, and suppressing a real finding is
    # the worse error -- so below MIN_SITES_FOR_MISKEY we report what we measured.
    MIN_SITES_FOR_MISKEY = 10
    MATCH_RATE_FLOOR = 5          # fewer than 1 in 5 matching means the key, not the data
    checked = [s for s in sites.values() if s.reported_alignments]
    matched = [s for s in checked if s.reported_alignments & s.seen_alignments]
    mis_keyed = False
    if (len(checked) >= MIN_SITES_FOR_MISKEY
            and len(matched) * MATCH_RATE_FLOOR < len(checked)):
        mis_keyed = True
        say("haplotype-scores: WARNING the reported alignments match a written haplotype "
            "on only %d of %d site(s). That is a mis-keyed comparison, not %d real "
            "failures: the curated table's Aligned_protospacer+PAM is not the same "
            "representation as the per-haplotype alignment (measured: it equals the ALT "
            "column on 0 of 244 rows of the delivering run). Reporting "
            "row_reported_alignment_found=unknown rather than a wall of false 'no'."
            % (len(matched), len(checked), len(checked) - len(matched)))

    n_rows = 0
    unreconciled = []
    with io.open(out_path, "w", encoding="utf-8", newline="") as out:
        out.write(file_header(source_path, source_kind, budget, threshold, recall,
                              cap, strict))
        writer = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for site in sites.values():
            ranked = sorted(
                site.rows,
                key=lambda r: (-(float(r["crispr_bulge"])
                                 if r["crispr_bulge"] not in ("NA", "") else -1.0),
                               r["variant_set"]))
            kept = ranked if not cap else ranked[:cap]
            hidden = len(ranked) - len(kept)
            found = "unknown" if mis_keyed else site.row_reported_alignment_found
            if found == "no":
                unreconciled.append(site.key)
            for rank, row in enumerate(kept, 1):
                row["rank_in_site_by_score"] = str(rank)
                row["n_haplotypes_at_site"] = str(site.n_haplotypes)   # the TRUE total
                row["n_above_threshold"] = str(site.n_above)
                row["n_not_scored"] = str(site.n_not_scored)
                row["site_max_crispr_bulge"] = fmt_score(site.max_cb)
                row["row_reported_alignment_found"] = found
                writer.writerow(row)
                n_rows += 1
            if hidden:
                say("haplotype-scores: %s: %d haplotype(s) not written (cap=%d); the "
                    "counts above are the true totals" % (site.key, hidden, cap))

    for key in unreconciled:
        say("haplotype-scores: WARNING %s: the site's own reported alignment is on none "
            "of the haplotypes written -- this list is what could be reconstructed, not "
            "necessarily all of them" % key)
    if unreconciled:
        say("haplotype-scores: %d site(s) did not reconcile." % len(unreconciled))
    return n_rows, len(sites), unreconciled


def _postprocess_roots():
    """Where ``cb_operating_point`` may live: beside a source checkout, in an install, or
    wherever ``CRISPRME_POSTPROCESS`` says. Named, not guessed."""
    here = os.path.dirname(os.path.abspath(__file__))
    roots = []
    env = os.environ.get("CRISPRME_POSTPROCESS")
    if env:
        roots.append(env)
    roots.append(os.path.normpath(os.path.join(here, "..", "PostProcess")))
    roots.append(os.path.join(here, "PostProcess"))
    roots.append("/opt/conda/opt/crisprme/PostProcess")
    seen, out = set(), []
    for r in roots:
        if r and r not in seen:
            seen.add(r)
            out.append(r)
    return out


def _threshold():
    """The sensitive operating point, BY NAME. Never a literal.

    Resolved by name so a file written today still means "the sensitive threshold" after
    the operating point is re-fitted. If it cannot be found we REFUSE: a guessed number
    in this column is worse than no file.
    """
    roots = _postprocess_roots()
    for root in roots:
        if root not in sys.path:
            sys.path.append(root)
    try:
        import cb_operating_point as cbop
    except ImportError as exc:
        raise HaplotypeScoreError(
            "cb_operating_point.py is not importable, so the sensitive operating point "
            "cannot be resolved by name and this file will not invent one. Looked in: "
            "%s. Set CRISPRME_POSTPROCESS to the folder that holds it, or run from an "
            "install (%s)" % (", ".join(roots), exc))
    return cbop.sensitive_threshold(), cbop.recall_target()


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="haplotype_scores.py",
        description="A row per observed haplotype, each with a score and a never-blank "
                    "status.")
    p.add_argument("--results", required=True, help="a finished run's Results/<job> folder")
    p.add_argument("--out", required=True, help="the .tsv to write")
    p.add_argument("--max-haplotypes", type=int, default=0, dest="cap",
                   help="cap rows written per site (0 = every one). Counts are always "
                        "the true totals.")
    p.add_argument("--score-source", choices=("harvest", "realign", "auto"),
                   default="harvest",
                   help="harvest (default) reads the scores the pipeline already wrote. "
                        "realign/auto additionally re-align the rows the fixed frame "
                        "cannot hold; measured coverage does not yet justify auto.")
    p.add_argument("--targets", default=None,
                   help="the run's curated results table (Chromosome, Position, "
                        "Aligned_protospacer+PAM), so each site's own reported alignment "
                        "can be reconciled against the haplotypes written. Without it "
                        "row_reported_alignment_found is 'unknown', never 'yes'.")
    p.add_argument("--strict", action="store_true",
                   help="exit non-zero if any site's reported alignment is on none of "
                        "its written haplotypes (default: warn loudly, exit 0).")
    p.add_argument("--anonymise", "--anonymize", action="store_true", dest="anonymise",
                   help="write carrier counts but no sample identities.")
    args = p.parse_args(argv)

    try:
        n_rows, n_sites, unreconciled = build(
            args.results, args.out, cap=args.cap, score_source=args.score_source,
            strict=args.strict, anonymise=args.anonymise, targets=args.targets)
    except HaplotypeScoreError as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 2
    print("haplotype-scores: %d row(s) over %d site(s) -> %s"
          % (n_rows, n_sites, args.out), file=sys.stderr)
    if unreconciled and args.strict:
        print("haplotype-scores: --strict: %d site(s) did not reconcile."
              % len(unreconciled), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":      # a module-level guard: see backlog 57
    sys.exit(main())
