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

  NOT BUILT, AND NOT CLAIMED: anchored re-alignment. Rows the fixed frame cannot hold
  (an indel inside the scored footprint, the target not on this haplotype, no harvested
  score) would need re-aligning on the haplotype's own sequence, anchored to the row's
  published alignment -- never an argmax, which reproduces the published reference CFD on
  4 of 186 rows while OVERSTATING it on 74. That is designed and NOT implemented here.
  This file says so rather than shipping a flag that accepts the word and does nothing;
  such rows carry score_status=needs-realignment:<site> and no invented number.

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
import zlib

# --------------------------------------------------------------------------- #
# score_status -- the vocabulary. A cell is NEVER blank and NEVER -1.0.

ANONYMISED = "anonymised"       # --anonymise: counts kept, identities withheld
ST_SCORED = "scored"
# A row can be BOTH scored and have had an axis decline it. Two previous attempts let one
# of those facts erase the other: making the refusal win counted a haplotype scoring 20x
# the operating point as never scored, and then making the score win rendered a reached-
# and-declined axis identically to one that has no column at all. The status says both.
ST_SCORED_REFUSED = "scored(%s-model-refused)"
ST_OUT_OF_DOMAIN = "out-of-domain:%s"          # the model's own refusal rule
ST_NOT_SCORED_PAM = "not-scored:scores-not-defined-for-this-pam(%s)"
ST_NEEDS_REALIGN = "needs-realignment:%s"       # site_key; Tier 2 would resolve it

FRAME_CANDIDATE = "candidate_frame"

# The CRISPR-Bulge null sentinel. It must never reach a column or a sort.
CB_NULL = -1.0

COLUMNS = [
    "site_key", "chrom", "site_start", "strand", "spacer",
    "haplotype_id", "variant_set",
    "aligned_spacer_pam", "aligned_protospacer_pam",   # make every number re-verifiable
    "alignment_frame",
    "cfd", "crispr_bulge", "score_status",
    "above_threshold",
    "rank_in_site_by_score", "rank_score_axis",
    "n_carriers", "carriers_homozygous", "carriers_heterozygous",
    "carrier_panels", "carrier_populations", "carrier_superpopulations",
    "carriers",
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


# The budget is DESCRIPTIVE: it is rendered into the file header and quoted inside
# score_status, and nothing here computes with it. So a missing key costs its own value,
# not the file. Refusing outright made the companion impossible to produce for an
# assembly-search combined run, whose .Params.txt has no Max_total_edits by construction
# (pages/history_page.py and results_page.py both carry a fallback for that same gap),
# and the message blamed "the search's budget (mismatches, bulges, PAM)" -- sending the
# operator to look at search parameters that were all present.
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
    out["_missing"] = missing
    for k in ("mismatches", "dna", "rna", "max_total"):
        if k not in out:
            continue                      # absent: budget_text renders it "?"
        try:
            out[k] = int(out[k])
        except ValueError:
            # Present but not a number is a DIFFERENT fact from absent, and a value this
            # file cannot render honestly, so it still refuses.
            raise HaplotypeScoreError(
                "%s records %s as %r, which is not a number." % (path, k, out[k]))
    if missing:
        say_missing = ", ".join(missing)
        out["_note"] = (
            "%s does not record %s. Those are LABELS -- they are printed in this file's "
            "header and quoted in score_status, and nothing is computed from them -- so "
            "they render as '?' rather than costing you the file. An assembly-search "
            "combined run has no Max_total_edits by construction."
            % (path, say_missing))
    return out


def budget_text(budget):
    """``6mm/1DNA/1RNA/8total`` -- quoted INSIDE score_status, so a row that found no
    alignment says which budget failed to find one."""
    def v(k, suffix):
        # "?" where the run did not record it -- never a number this file invented.
        return ("%s%s" % (budget[k], suffix)) if budget.get(k) is not None else ("?" + suffix)
    return "".join((v("mismatches", "mm/"), v("dna", "DNA/"),
                    v("rna", "RNA/"), v("max_total", "total")))


# =========================================================================== #
# ONE VERDICT TABLE. Every value these three columns can hold, with the condition that
# produces it and what it means -- declared ONCE and rendered into the file header from
# this same table.
#
# Three of five review rounds found the header promising a value the code cannot emit, or
# defining one in terms that were no longer true, because the legend was prose written
# beside the code rather than generated from it. And four rounds found two genuinely
# different absences rendering alike -- which is hard to notice when the possible values
# are scattered across a dozen return statements and easy to see when they are a table.
#
# A value that is not in this table cannot reach the file: _verdict() asserts it.

VERDICTS = (
    # column, value, meaning (rendered into the header verbatim)
    ("score_status", ST_SCORED,
     "at least one axis produced a number"),
    ("score_status", "scored(<axis>-model-refused)",
     "scored, AND the named axis was reached and declined. UNREACHABLE with the schema "
     "CRISPRme writes today -- each block carries exactly one score column, so no "
     "observation can have a number on one axis and a refusal on the other. Kept "
     "because it is the correct answer if that ever changes, and listed here so the "
     "value is never a surprise; see test_the_composite_status_is_unreachable_today"),
    ("score_status", "out-of-domain:model-refused(n/length/multi-gap)",
     "reached, declined, and NO axis produced a number"),
    ("score_status", "not-scored:scores-not-defined-for-this-pam(<pam>)",
     "NO SCORE WAS COMPUTED for this run, on either axis: the pipeline scores only a "
     "20 nt guide with a 3 nt 3-prime PAM, and 11 of the 17 shipped PAM files are not "
     "that (every Cas12a and SaCas9 among them). The same -1 a model writes to decline "
     "a target also means this, and calling it a refusal asserts a judgement nothing "
     "made"),
    ("score_status", "needs-realignment:<site>",
     "the fixed alignment frame cannot hold it; no number, and none invented"),

    ("above_threshold", "yes",
     "measured, and at or above the CRISPR-Bulge operating point"),
    ("above_threshold", "no",
     "measured, and below it"),
    ("above_threshold", "unknown",
     "NOT measured: this frame carries no CRISPR-Bulge column at all"),
    ("above_threshold", "unknown(model-refused)",
     "measured in the sense that the model reached it, and DECLINED it. A different "
     "fact from never measured, and on a real run the commoner one"),
    ("above_threshold", "unknown(not-computed-for-this-pam)",
     "no score was computed for this RUN: the pipeline scores only a 20 nt guide with a "
     "3 nt 3-prime PAM, and nulls both axes otherwise. Neither measured nor declined"),

    ("row_reported_alignment_found", "yes",
     "the site's own reported haplotype is among those written"),
    ("row_reported_alignment_found", "no",
     "it is not -- this list is what could be reconstructed, not necessarily all"),
    ("row_reported_alignment_found", "n/a",
     "the curated table SAYS ref, so there is no ALT haplotype for it to appear "
     "among. A measured fact, never used for an absence"),
    ("row_reported_alignment_found", "unknown",
     "nothing said what this site reported: either no --targets, or --targets names "
     "the locus only under a DIFFERENT guide (stderr says which; this file cannot)"),
    ("row_reported_alignment_found", "unknown(origin-not-in-targets)",
     "the curated table has NO REF/ALT_origin column, so ref was never read"),
    ("row_reported_alignment_found", "unknown(unreadable-origin-in-targets)",
     "the curated table HAS that column and this row's value is neither ref nor alt. "
     "A different fact from the column being absent, and it used to share its value"),
    ("row_reported_alignment_found", "unknown(unparsable-variant-in-targets)",
     "an ALT row whose Variant cell this file cannot reduce to records"),
    ("row_reported_alignment_found", "unknown(ambiguous-in-targets)",
     "the curated table described that locus twice, incompatibly"),
    ("row_reported_alignment_found", "unknown(guide-ambiguous-in-targets)",
     "the source carries SEVERAL guides at this locus and the curated table names no "
     "guide column, so which of them its row describes cannot be known. Answering "
     "yes/no would be answering about a guide the table never described"),
    ("row_reported_alignment_found", "unknown(key-looks-wrong)",
     "the comparison looks mis-keyed across the run, so a 'no' is not reported as one"),
    ("row_reported_alignment_found", "yes(among-all-haplotypes-not-all-shown)",
     "found, but a cap hid the matching haplotype from this file"),
)

_VERDICT_VALUES = {}
for _col, _val, _why in VERDICTS:
    _VERDICT_VALUES.setdefault(_col, []).append(_val)


def _verdict(column, value):
    """Return ``value``, having checked the table declares it for ``column``.

    Values carrying a variable part are declared with a <placeholder>; the check is on
    the fixed prefix. The point is that adding a value to the code without adding it to
    the table -- and so to the header a reader is handed -- fails here, loudly, instead
    of shipping a cell nothing explains.
    """
    for declared in _VERDICT_VALUES.get(column, ()):
        if declared == value:
            return value
        head = declared.split("<")[0]
        if "<" in declared and value.startswith(head):
            return value
    raise HaplotypeScoreError(
        "internal: %r is not a declared value of %s. Add it to VERDICTS (with what it "
        "means) so the file header explains it, or do not write it." % (value, column))


def verdict_legend(column, width=86):
    """The header legend for one column, generated FROM the table.

    Every line starts with '#': the whole header is comment lines, and a continuation
    that did not would be read as data by anything skipping comments.
    """
    out = []
    for col, val, why in VERDICTS:
        if col != column:
            continue
        head = "#     %-40s -- " % val
        pad = "#" + " " * (len(head) - 1)
        cur = head
        for w in why.split():
            if len(cur) + 1 + len(w) > width and cur not in (head,):
                out.append(cur.rstrip())
                cur = pad + w
            else:
                cur = (cur + w) if cur in (head,) else (cur + " " + w)
        out.append(cur.rstrip())
    return out

# --------------------------------------------------------------------------- #
# carriers -- a reagent list, reported plainly


PANEL_1000G = re.compile(r"^(NA|HG)\d")
PANEL_HGDP = re.compile(r"^HGDP\d")


def carrier_breakdown(text, ancestry=None):
    """``(n, homozygous, heterozygous, panels)`` from a ``SAMPLE:slot`` carrier string.

    WHY THIS IS A COLUMN AND NOT A NOTE. 1000G and HGDP lines are commercially available, so
    a carrier id is something a reader can act on: test the haplotype in a genome that
    already has it instead of engineering it. Zygosity decides how they QUANTIFY, not
    whether they can test -- allele-specific quantification handles a heterozygote fine.
    So this is REPORTED, never ranked on and never used to filter or warn.

    A sample appearing on both slots is homozygous for that haplotype. ``panels`` says
    which collection an id belongs to: not every carrier id is a standard accession
    (1,177 of 39,446 tokens on the delivering run were neither 1000G- nor HGDP-shaped),
    so the reader can see what they are looking at rather than be told it is orderable.
    """
    per = collections.defaultdict(set)
    any_slot = False
    for token in (text or "").split(","):
        token = token.strip()
        # The pipeline writes a no-carrier cell as one of these, not as an empty string.
        # Counting them gave every reference-only or unlabelled row exactly 1 carrier --
        # a fabricated individual, in the column a reader would use to order a cell line.
        if not token or token.upper() in ("NA", "N", "-", ".", "NO_SAMPLES", "NONE"):
            continue
        if ":" in token:
            sample, slot = token.rsplit(":", 1)
            any_slot = True
        else:
            sample, slot = token, "?"
        per[sample.strip()].add(slot.strip())
    if not per:
        return 0, "0", "0", "-", "unknown", "unknown"
    if not any_slot:
        # THE SOURCE DOES NOT SAY. Only the SAMPLE:slot form (extract-flanks) records which
        # chromosome copy carries the haplotype; the per-haplotype results columns give bare
        # sample names. Printing 0 homozygous here would read as "none are", which is the
        # exact failure this file exists to prevent: a not-measured cell that looks like a
        # measured zero.
        return (len(per), "unknown", "unknown", _panel_text(per),
                _ancestry_text(per, ancestry, "population"),
                _ancestry_text(per, ancestry, "superpopulation"))
    hom = sum(1 for slots in per.values() if len(slots) >= 2)
    return (len(per), str(hom), str(len(per) - hom), _panel_text(per),
            _ancestry_text(per, ancestry, "population"),
            _ancestry_text(per, ancestry, "superpopulation"))


def read_ancestry(results_dir):
    """``{sample: {...}}`` from the run's own ``.sampleID.txt``, or ``{}``.

    The run writes the ancestry of every individual it searched, and no delivered file
    carries it: the report keeps only the superpopulation for one embedded chart, and the
    shared panel reader keeps only the sample id. Reading it here costs one small file and
    turns a carrier id into something a reader can reason about. An absent file means
    absent columns, never an invented label.
    """
    path = os.path.join(results_dir, ".sampleID.txt")
    if not os.path.isfile(path):
        return {}
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.append(here)
    try:
        from run_records import read_panel_meta
    except ImportError:
        return {}
    try:
        return read_panel_meta([path])
    except (OSError, ValueError):
        return {}


def _ancestry_text(per, ancestry, field):
    """``EUR:12,AFR:3`` over the carriers, or ``unknown`` when the panel file is absent.

    'unknown', never blank and never an empty count: an absent panel file must not read
    as an undiverse panel.
    """
    if not ancestry:
        return "unknown"
    counts = collections.Counter()
    for sample in per:
        label = (ancestry.get(sample) or {}).get(field) or "unlabelled"
        counts[label] += 1
    return ",".join("%s:%d" % (k, counts[k]) for k in sorted(counts)) or "-"


def _max_af(text):
    """The LARGEST allele frequency in the cell, not the cell.

    The column is named max_af and was the raw comma-joined per-variant list, so a
    haplotype of three variants printed "0.01,0.4,0.02" where a reader expected one
    number -- and any numeric filter on it silently matched nothing.
    """
    best = None
    for piece in str(text or "").split(","):
        piece = piece.strip()
        if not piece or piece.upper() in ("NA", "N", "-", "."):
            continue
        try:
            v = float(piece)
        except ValueError:
            continue
        if best is None or v > best:
            best = v
    return "NA" if best is None else "%.6g" % best


def _panel_text(per):
    panels = collections.Counter()
    for sample in per:
        if PANEL_HGDP.match(sample):
            panels["HGDP"] += 1
        elif PANEL_1000G.match(sample):
            panels["1000G"] += 1
        else:
            panels["other"] += 1
    return ",".join("%s:%d" % (k, panels[k]) for k in sorted(panels)) or "-"


# --------------------------------------------------------------------------- #
# the input source -- probe, record, or fail loudly


RAW_SIDECAR_RE = re.compile(r"_observed_haplotypes_raw\.(tsv|zip)$")
ALT_TABLE_RE = re.compile(r"_all_results_with_alternative_alignments\.(tsv|zip)$")


def probe_source(results_dir, say=None):
    """``(path, kind)`` for the richest haplotype source present.

    Order matters. ``*_observed_haplotypes_raw`` is written at search time and filters
    ONLY on origin -- an explicit UNSCOPED superset of the slim companion, which is
    additionally restricted to within 50 bp of a CFD-top-1000 site. Falling back to the
    full alt table is correct but larger and slower.

    A HEADER-ONLY candidate is skipped, not preferred. The raw sidecar is gated on a
    Tier-1 dict-less genotype store, so a run without that store writes the header and no
    rows -- and preferring it reported 0 haplotypes for a run holding 16,208 ALT rows with
    real carriers (measured on the chr1/7/22 NRG validation run). That is a FALSE CLEAN
    produced by the file written to prevent false cleans, so emptiness disqualifies a
    source rather than silencing the next one.

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
    # An unreadable candidate costs ITS OWN turn, not the fallback's. Raising here aborted
    # probe_source before the alt table was ever looked at, so a corrupt raw sidecar beside
    # a complete alt table refused the whole run.
    say = say or (lambda m: None)
    unreadable = []
    for paths, kind in ((raw, "observed_haplotypes_raw"),
                        (alt, "all_results_with_alternative_alignments")):
        for path in paths:
            try:
                if has_data_row(path):
                    # A SILENT FALLBACK IS STILL A DATA PROBLEM. An interrupted search
                    # leaves a truncated sidecar beside a complete alt table; using the
                    # second without a word would let a damaged file sit there unnoticed
                    # in a module whose whole doctrine is that an absence is stated.
                    for msg in unreadable:
                        say("haplotype-scores: WARNING %s" % msg)
                    if unreadable:
                        say("haplotype-scores: fell back to %s. The run SUCCEEDED, but "
                            "the preferred source above is damaged and should be "
                            "regenerated." % os.path.basename(path))
                    return path, kind
            except HaplotypeScoreError as exc:
                unreadable.append(str(exc))
    if unreadable and len(unreadable) == len(raw) + len(alt):
        raise HaplotypeScoreError(
            "no haplotype source in %s can be read. %s Refusing to write: an unreadable "
            "file is a different problem from an empty run."
            % (results_dir, " ".join(unreadable)))
    if unreadable:
        raise HaplotypeScoreError(
            "no haplotype source in %s could be used. %s%s Refusing to write: an "
            "unreadable file is a different problem from an empty run, and an empty "
            "companion would be read as the second."
            % (results_dir, " ".join(unreadable),
               (" The remaining candidates are header-only." if len(unreadable) < len(raw) + len(alt)
                else "")))
    if raw or alt:
        raise HaplotypeScoreError(
            "every haplotype source in %s is header-only (%s). A source with no data "
            "rows cannot distinguish 'this run found no variant haplotypes' from 'the "
            "writer never ran', and an empty companion would be read as the first. "
            "Refusing to write.\n"
            "The run's integrated_results.tsv is NOT a substitute: it carries one row "
            "per site holding each axis's single best representative, not every "
            "haplotype, so building this file from it would report three sequences per "
            "site and call that the enumeration. If the report shows variant rows but "
            "this refuses, the search ran without the Tier-1 dict-less genotype store "
            "(or in population-level fast mode) and the haplotypes were never written "
            "down; re-run the search with per-sample enumeration enabled."
            % (results_dir, ", ".join(os.path.basename(p) for p in (raw + alt))))
    raise HaplotypeScoreError(
        "neither *_observed_haplotypes_raw.(tsv|zip) nor "
        "*_all_results_with_alternative_alignments.(tsv|zip) is in %s. Per-haplotype "
        "enumeration is gated on a Tier-1 dict-less genotype store and is ABSENT "
        "ENTIRELY in population-level fast mode, so an empty companion here would read "
        "as 'no extra high-scoring haplotypes' when the truth is 'this run never "
        "enumerated any'. Refusing to write." % results_dir)


def has_data_row(path):
    """True if the table has at least one non-blank line after its header.

    An UNREADABLE source raises. Swallowing the error returned False, and probe_source
    then called a corrupt archive "header-only" and sent the operator off to re-run the
    search with per-sample enumeration -- for a bad zip.
    """
    try:
        with open_table(path) as stream:
            next(stream, None)                     # the header
            for line in stream:
                if line.strip():
                    return True
    except (OSError, EOFError, zipfile.BadZipFile, zlib.error) as exc:
        raise HaplotypeScoreError(
            "%s cannot be read (%s). This is not an empty run: the file is there and "
            "unreadable, which is a different problem from 'nothing was enumerated'. "
            "Refusing to guess which." % (path, exc))
    return False


def open_table(path):
    """Stream a .tsv, or the single member of a .zip, as text lines."""
    if path.endswith(".zip"):
        zf = zipfile.ZipFile(path)
        names = [n for n in zf.namelist() if n.endswith(".tsv")]
        if not names:
            zf.close()
            raise HaplotypeScoreError("%s holds no .tsv member" % path)
        handle = None
        try:
            handle = zf.open(names[0])
            wrapper = io.TextIOWrapper(handle, encoding="utf-8", errors="replace")
        except Exception:
            # ZipFile.close() does not close an already-open member, so both are needed.
            if handle is not None:
                handle.close()
            zf.close()
            raise
        # Closing the wrapper must close the archive too: probe_source now opens one per
        # candidate, so a leaked ZipFile per call is a real descriptor leak.
        _close = wrapper.close

        def close_both():
            try:
                _close()
            finally:
                zf.close()
        wrapper.close = close_both
        return wrapper
    return io.open(path, "r", encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------- #
# the three alignment blocks


# QUOTE_NONE EVERYWHERE, and that is the point. A TSV is not a CSV: a results table has
# no quoting convention, a double quote in a cell is a double quote, and csv's default
# dialect swallows it -- and where a cell opens a quote that closes on a later line, it
# JOINS those lines into one row, silently dropping every row between. One stray quote
# in a sample list or an alignment therefore truncated the companion with no error and
# no short-row to trip the truncation guard. PostProcess/extract_flanks.py:310 already
# documents this exact failure for these exact tables; this file simply had not followed
# it.
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
# The axis that defines a site's identity. The curated tables publish this block's
# coordinate as "Position", so a --targets list and these rows agree by construction.
ANCHOR_BLOCK = "highest_CFD"


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


def _is_coordinate(text):
    """A real genomic start, not one of this pipeline's several ways of saying 'unset'."""
    t = str(text or "").strip()
    return bool(t) and t.lstrip("+-").isdigit()


def site_key(chrom, start, strand, spacer):
    return "%s:%s:%s:%s" % (chrom, start, strand, spacer)


# --------------------------------------------------------------------------- #
# Tier 1 -- harvest


Observation = collections.namedtuple(
    "Observation",
    "site chrom start anchor_start strand spacer block variants aligned_spacer "
    "aligned_alt cfd crispr_bulge samples maf origin")


def harvest_rows(stream, spacer_hint=None):
    """Yield one ``Observation`` per (row, block) the table actually names.

    Takes the score ONLY from the axis that named the haplotype. The integrated
    projections carry CFD on the CFD axis and CRISPR-Bulge on the CB axis with no
    cross-axis score, which is exactly why most rows of the delivered companion are NA on
    one of them. Promising both from this source would be inventing a number.
    """
    reader = csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)
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

    width = len(header)
    for row in reader:
        if not row or (len(row) == 1 and not row[0].strip()):
            continue            # a genuinely blank line, the only one we skip
        # Anything else shorter than the header is a CUT, and the refusal below says so.
        # Skipping rows of fewer than three fields meant a file cut inside its first
        # cells was accepted silently and the companion written as if complete.
        if len(row) < width:
            # A SOURCE CUT MID-ROW IS NOT A ROW WITH NO SCORES. Short rows used to fall
            # through cell()'s bounds check and emit a complete-looking NA/unknown row,
            # so "the file was truncated here" rendered as "this haplotype was never
            # measured" -- two different absences, one of them invented.
            raise HaplotypeScoreError(
                "this table is truncated: a row has %d of %d fields (first cells: %s). "
                "Rows after the cut would be written as never-measured, which is a "
                "different fact from a damaged file. Refusing rather than reporting a "
                "cut as a result." % (len(row), width, ",".join(row[:3])))
        chrom = row[chrom_i].strip() if chrom_i is not None and chrom_i < len(row) else ""
        spacer = (row[spacer_i].strip() if spacer_i is not None and spacer_i < len(row)
                  else (spacer_hint or ""))
        # SITE IDENTITY COMES FROM ONE AXIS, FOR EVERY BLOCK OF THE ROW.
        #
        # The three blocks are independent projections and place the same physical
        # off-target at different starts: MEASURED at 2.1% of 73,994 rows on a genome-wide
        # per-sample run and 2.6% of 31,502 on a three-chromosome one, almost all of them
        # a 1-2 bp shift. Keying the site on the block being read split one site into
        # up to three, so the "per-site"
        # roll-ups were per-FRAME -- and worse, with --targets the filter keyed on the
        # per-block start while the curated table's Position is the highest-CFD start, so
        # the CRISPR-Bulge observation was dropped and a haplotype scoring 0.9, thousands
        # of times the operating point, rendered as having no CRISPR-Bulge score at all.
        # A false clean, inside the file written to prevent false cleans.
        #
        # The anchor is the highest-CFD block because that is what the curated tables
        # publish as Position, so --targets and these rows agree by construction. Each
        # observation still reports its OWN block's coordinate in site_start.
        anchor = blocks.get(ANCHOR_BLOCK) or {}
        anchor_start = cell(row, anchor, "start") if anchor else ""
        anchor_strand = cell(row, anchor, "strand") if anchor else ""
        # "NA" IS NOT A COORDINATE. cell() returns it as a truthy string, so
        # `anchor_start or start` kept it and every row with an unset anchor was keyed
        # chr1:NA:NA:<guide> -- merging physically distinct loci into one site, stamping
        # one locus's site_max_crispr_bulge onto another's rows, and letting a cap delete
        # whole loci from the file. This pipeline writes NA for "not set"; it is an
        # absence, and an absence cannot anchor anything.
        if not _is_coordinate(anchor_start):
            # ONE ROW IS ONE SITE, even when the preferred anchor block is unset.
            # Falling back to each block's OWN start re-split the row into up to three
            # sites and re-opened the --targets false clean the anchor was introduced to
            # close: the coarse filter compares the curated Position against the block's
            # start, so the CRISPR-Bulge observation was dropped and its score with it.
            # The next block that HAS a coordinate anchors the whole row instead.
            anchor_start, anchor_strand = "", ""
            for _b in BLOCKS:
                _c = blocks.get(_b)
                if not _c:
                    continue
                _s = cell(row, _c, "start")
                if _is_coordinate(_s):
                    anchor_start, anchor_strand = _s, cell(row, _c, "strand")
                    break
        for block, cols in blocks.items():
            start = cell(row, cols, "start")
            if not _is_coordinate(start):
                continue       # this block did not place the target anywhere
            strand = cell(row, cols, "strand")
            variants = cell(row, cols, "variants")
            origin = cell(row, cols, "origin")
            yield Observation(
                site=site_key(chrom, anchor_start or start,
                              anchor_strand or strand, spacer),
                chrom=chrom, start=start, anchor_start=anchor_start or start,
                strand=strand, spacer=spacer, block=block,
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


def scores_are_defined(budget):
    """Whether this run's shape lets the pipeline compute ANY score -- the ``do_scores``
    flag, reconstructed from the run's own PAM.

    BOTH AXES, not just CFD. ``new_simple_analysis.py`` sets ``do_scores = False`` when
    the PAM is not 3 nt, the guide is not 20 nt, or the PAM is 5-prime, and passes the
    same flag to both writers: ``preprocess_CFD_score`` writes -1.000, and
    ``cbulge_rows.py:103`` returns -1 for every row BEFORE the model is loaded --
    "every row then gets -1", in its own words. Eleven of the seventeen shipped PAM
    files are in that set.

    This was first written as cfd_is_defined and applied only where a CFD sentinel was
    seen. The CRISPR-Bulge block carries NO CFD column, so the check could never fire on
    the one frame holding the CRISPR-Bulge axis: on every Cas12a and SaCas9 run those
    rows still said the model reached the target and declined it -- naming three
    specific reasons (N / length / multi-gap) nothing had evaluated -- while the two
    frames beside them correctly said it was never computed, from the identical -1.
    """
    pam = str(budget.get("pam") or "").strip()
    return len(pam) == 3 if pam else True      # unknown PAM: do not invent a verdict


def classify(obs, budget, axis_named):
    """``(cfd, crispr_bulge, score_status, cb_refused)`` for one harvested observation.

    ``cb_refused`` is returned as its OWN value rather than left to be recovered from the
    status string. The caller used to test ``"cb-model-refused" in status``, which is
    false for the out-of-domain status -- so on the 139,656 real rows where the
    CRISPR-Bulge axis was reached and declined, above_threshold read a bare "unknown",
    which this file's header defines as NEVER MEASURED. The flagship invariant, broken by
    sniffing a string for a fact that was already known.

    ``axis_named`` says which axis named this haplotype. The other axis gets no number
    and says so -- never a blank, never a borrowed value.
    """
    def num(text):
        """``(value, refused)``. ``refused`` is True ONLY for the model's -1 sentinel.

        A cell that is blank or NA was never measured; a cell holding the sentinel WAS
        reached and the model declined it by its own rule (N in the target, a length
        mismatch, more than one gap column, over-length). Those are different facts and
        the old code collapsed both to None, which made the out-of-domain branch
        unreachable and typed a refusal as "needs-realignment".
        """
        if text in ("", "NA", "n", "-", "."):
            return None, False
        try:
            v = float(text)
        except ValueError:
            return None, False
        if v <= CB_NULL:
            return None, True
        return v, False

    cfd, cfd_refused = num(obs.cfd)
    cb, cb_refused = num(obs.crispr_bulge)
    # The CFD writer ALSO emits the -1 sentinel ("no score calculated, append -1" --
    # new_simple_analysis.py), and CFD is the one column present on all three blocks. On
    # the two blocks with no CB column a refused CFD left cfd=None, cb=None and typed the
    # row "needs-realignment", which says the frame could not hold it; the truth is the
    # scorer reached it and declined. Same conflation, on the commoner axis.
    # A REFUSAL ON ONE AXIS IS NOT A REFUSAL OF THE ROW. Marking the row out-of-domain
    # whenever either axis held the sentinel made a haplotype with a real CRISPR-Bulge
    # score of 0.004 -- twenty times the operating point -- count as NEVER SCORED, so the
    # same row was both above_threshold=yes and part of n_not_scored, and a reader keeping
    # only score_status=scored discarded the one row above the threshold. Out-of-domain is
    # for a row where nothing survived.
    if cfd is not None or cb is not None:
        # Scored -- and if the OTHER axis was reached and declined, the status says which.
        declined = [n for n, r in (("cb", cb_refused), ("cfd", cfd_refused)) if r]
        if declined:
            return cfd, cb, ST_SCORED_REFUSED % "+".join(declined), cb_refused
        return cfd, cb, ST_SCORED, cb_refused
    if (cfd_refused or cb_refused) and not scores_are_defined(budget):
        # The -1 is the "not calculated" one, on EITHER axis, so this is checked before
        # the refusal branch below and without asking which axis carried it.
        return (None, None,
                ST_NOT_SCORED_PAM % (str(budget.get("pam") or "?").strip() or "?"),
                False)          # NOT a refusal: nothing reached this row to refuse it
    if cb_refused or cfd_refused:
        return (None, None,
                ST_OUT_OF_DOMAIN % "model-refused(n/length/multi-gap)", cb_refused)
    return None, None, ST_NEEDS_REALIGN % obs.site, cb_refused


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
                 "reported_alignments", "seen_alignments", "applicable", "has_report",
                 "ambiguous", "unreadable_report", "anchor")

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
        self.ambiguous = False
        self.unreadable_report = ""
        # The site's ANCHOR coordinate, which is what every --targets set is keyed on.
        # site.start is the creating observation's own block start, and the two differ
        # on about 2% of real rows -- looking the site up by site.start is the exact
        # mistake that produced three earlier findings.
        self.anchor = ""

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
        #   unknown -- nothing told us what this site reported (no --targets), or the
        #              curated table described this locus twice, incompatibly
        #   n/a     -- it reported a REFERENCE haplotype, so there is no ALT haplotype
        #              for it to appear among; the check does not apply
        #   yes/no  -- it reported an ALT haplotype, and it is or is not among them
        if self.ambiguous:
            return "unknown(ambiguous-in-targets)"
        if self.unreadable_report:
            # NOT n/a. n/a asserts "it reported a REFERENCE haplotype, so there is no ALT
            # for it to appear among" -- a measured fact. Folding "the table has no origin
            # column" and "the Variant cell did not parse" into it asserted a reason
            # nobody read, left unreconciled empty, printed nothing, and made --strict --
            # the flag whose only job is to fail here -- exit 0 on a genuine failure.
            return "unknown(%s)" % self.unreadable_report
        if not self.has_report:
            return "unknown"
        if not self.applicable:
            return "n/a"
        return "yes" if (self.reported_alignments & self.seen_alignments) else "no"


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
        "#   It is an operating point, NOT a safe/unsafe line: at the shipped 95%% recall",
        "#   target its guide-cluster bootstrap CI spans [5.8e-06, 3.1e-03], three orders",
        "#   of magnitude. That interval belongs to the DEFAULT target; if",
        "#   CRISPRME_CB_RECALL_TARGET moved it, the threshold above is re-resolved and",
        "#   this interval is not. Scores are written %.6g, never 3 dp -- at 3 dp the",
        "#   threshold itself rounds to 0.000.",
        "#",
        "# THREE VALUES, NOT TWO -- and the third is the common one. MEASURED on a real",
        "#   NRR run: of 460,714 rows, 139,656 read unknown(model-refused) and 304,003 a",
        "#   bare unknown. Before they were one value, so 30% of the file said it had",
        "#   never measured rows the model had reached and declined.",
    ] + verdict_legend("above_threshold") + [
        "#",
        "#   row_reported_alignment_found:",
    ] + verdict_legend("row_reported_alignment_found") + [
        "#",
        "# RANK IS ONE ORDERING, NOT ONE PER AXIS. rank_in_site_by_score numbers the",
        "#   site's rows 1..n in the order a cap takes them: CRISPR-Bulge-bearing rows",
        "#   first (by CRISPR-Bulge), then CFD-only (by CFD), then unscored. So rank 1 is",
        "#   the site's top CRISPR-Bulge row where it has one, and the top CFD row ranks",
        "#   below it -- NOT also 1. rank_score_axis says which score fixed each row's",
        "#   place. The two scales are never compared (CFD runs 0..1; the CB operating",
        "#   point is 1.9e-4), which is why a cap can no longer discard a CRISPR-Bulge",
        "#   row in favour of a larger-looking CFD one.",
        "#",
        "# NEVER A BLANK. score_status is always set. A never-scored cell and the",
        "#   CRISPR-Bulge out-of-domain sentinel both used to render 'NA', so filtering",
        "#   on the score silently discarded them and the site looked clean. -1.0 is never",
        "#   printed, never sorted on. Its values (generated from VERDICTS, so this",
        "#   legend cannot drift from what the code can write):",
    ] + verdict_legend("score_status") + [
        "#   above_threshold reads unknown(model-refused) where the CRISPR-Bulge axis was",
        "#   reached and declined, as against a bare unknown for never measured.",
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
        "# alignment_frame is ALWAYS candidate_frame: the score is the haplotype's bases",
        "#   in the CANDIDATE's alignment shape (the pipeline re-inserts gap columns at",
        "#   fixed positions and counts mismatches positionally). Anchored re-alignment on",
        "#   the haplotype's own sequence is designed and NOT implemented, so no row here",
        "#   carries a haplotype_own frame; rows the fixed frame cannot hold say",
        "#   score_status=needs-realignment and carry no invented number.",
        "#",
        "# CARRIERS ARE A REAGENT LIST, REPORTED, NEVER RANKED. A heterozygous carrier is",
        "#   as orderable as a homozygous one -- allele-specific quantification handles it",
        "#   -- so zygosity is written beside each haplotype and never used to filter or",
        "#   to order. carriers_homozygous/heterozygous say 'unknown' when THE SOURCE DOES",
        "#   NOT SAY: only the SAMPLE:slot form records which allele, and the search",
        "#   sidecars write bare sample names (measured: 153,879 rows with real carriers",
        "#   on an NRR run, 0 of them slot-encoded). 'unknown' is that fact; 0 would read",
        "#   as 'none are homozygous'. For zygosity, use the extract-flanks output.",
        "#",
        "# RE-VERIFY ANY ROW: aligned_spacer_pam + aligned_protospacer_pam feed",
        "#   'crisprme.py targets-scoring' directly.",
        "#",
        "# --targets NAMES SITES THIS FILE MAY NOT HOLD. The curated Position is the",
        "#   merged representative's start; this source is written pre-merge, so some",
        "#   sites sit a few bp away and some were never enumerated. Sites that produced",
        "#   no rows are COUNTED AND NAMED on stderr, never fuzzy-matched -- a window",
        "#   would merge distinct sites. Absent here is 'not enumerated', not 'clean'.",
        "#",
        "# row_reported_alignment_found=no means the site's OWN reported alignment is on",
        "#   none of the haplotypes written: the list is what we could reconstruct, not",
        "#   necessarily all of them. Measured on an NRR run over 97 curated sites: 75",
        "#   yes, 22 no -- so this is a routine finding, not a rare alarm.%s" % (
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

    THE FIELDS ARE TAKEN FROM THE RIGHT. hg38 ships ~500 contigs whose NAME contains the
    same underscore that separates the fields -- chr1_KI270706v1_random, chrUn_KI270302v1
    -- and splitting left-to-right took the first four pieces, so
    ``chr1_KI270706v1_random_12345_A_T`` reduced to ``chr1:KI270706v1:RANDOM:12345`` with
    THE ALLELES DROPPED. Two different variants at that position then compared EQUAL, and
    a site reconciled 'yes' against a haplotype it does not carry. A false yes is the one
    answer this file must never give, because it says the enumeration is complete when it
    is not. Position/REF/ALT are the last three fields; everything before them is the
    contig, however many separators it contains.
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
            # from the RIGHT: pos, REF, ALT are the last three; the rest is the contig
            pos, ref, alt = parts[-3], parts[-2], parts[-1]
            chrom = "_".join(parts[:-3])
            if not pos.isdigit():
                # Not the shape we think it is. Keep the whole piece verbatim rather than
                # inventing a key that might collide with a real one.
                out.append(piece.upper())
                continue
            out.append("%s:%s:%s:%s" % (chrom, pos, ref.upper(), alt.upper()))
    return tuple(sorted(out))


Reported = collections.namedtuple("Reported", "variants origin origin_known")

# a locus the curated table describes twice, incompatibly: never reconciled against
_AMBIGUOUS = Reported(variants=(), origin="ambiguous", origin_known=True)
# the curated table had no REF/ALT_origin column at all, as against a cell we cannot read
_ORIGIN_NO_COLUMN = "<no-origin-column>"


# the curated table's guide column, under any of the names it has shipped with
_TARGET_SPACER_COLS = ("Spacer+PAM", "Spacer+PAM_(highest_CFD)", "crRNA", "Guide")


def read_reported(targets_path, say=None):
    """``{(chrom, start, spacer): Reported}`` -- each site's own reported haplotype.

    Keyed on the VARIANT SET (see ``canonical_variants``), and carrying the row's
    REF/ALT origin: a REF-origin site's representative is the reference sequence, so it
    has no ALT haplotype to appear among and reconciliation does not apply to it. Over
    the delivering run that is 115 of 187 sites -- checking them would manufacture
    failures.
    """
    say = say or (lambda m: None)
    out = {}
    ambiguous = set()
    with io.open(targets_path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        try:
            header = next(reader)
        except StopIteration:
            return out
        idx = {n: i for i, n in enumerate(header)}
        c_i, p_i = idx.get("Chromosome"), idx.get("Position")
        v_i, o_i = idx.get("Variant"), idx.get("REF/ALT_origin")
        s_i = next((idx[n] for n in _TARGET_SPACER_COLS if n in idx), None)
        if c_i is None or p_i is None or v_i is None:
            raise HaplotypeScoreError(
                "%s is not a curated results table: it needs Chromosome, Position and "
                "Variant to say which haplotype each site reported." % targets_path)
        if o_i is None:
            say("haplotype-scores: NOTE %s has no REF/ALT_origin column, so no site can "
                "be shown to have reported a REFERENCE haplotype. Reconciliation for "
                "every site is 'unknown(origin-not-in-targets)', NOT 'n/a' -- n/a would "
                "assert a reason that was never read."
                % os.path.basename(targets_path))
        if s_i is None:
            say("haplotype-scores: NOTE %s names no guide column (%s), so two guides "
                "reporting at the same locus cannot be told apart; any such locus is "
                "left 'unknown(ambiguous-in-targets)' rather than guessed"
                % (os.path.basename(targets_path), "/".join(_TARGET_SPACER_COLS)))
        short = 0
        for row in reader:
            if max(c_i, p_i, v_i) >= len(row):
                # COUNTED, NOT DROPPED. A short curated row silently shrank the
                # shortlist the run then reported on, so the per-site accounting was
                # complete about a list the user did not hand over.
                short += 1
                continue
            # The column being ABSENT and the cell being unreadable are different
            # facts, and folding them into one value made the legend true of only one.
            origin = (row[o_i].strip().lower()
                      if o_i is not None and o_i < len(row) else _ORIGIN_NO_COLUMN)
            spacer = (row[s_i].strip() if s_i is not None and s_i < len(row) else "")
            key = target_key(row[c_i], row[p_i], spacer)
            rep = Reported(variants=canonical_variants(row[v_i]), origin=origin,
                           origin_known=(o_i is not None and origin in ("ref", "alt")))
            prev = out.get(key)
            if prev is not None and prev != rep:
                # Two curated rows claim the same (locus, guide) with DIFFERENT reported
                # haplotypes. Last-one-wins would silently reconcile against whichever
                # row happened to be later in the file, so neither is used.
                ambiguous.add(key)
            out[key] = rep
    if short:
        say("haplotype-scores: WARNING %d row(s) of %s are too short to read (they do "
            "not carry Chromosome, Position and Variant) and are NOT among the %d "
            "site(s) reported on below. The shortlist this run answers about is "
            "smaller than the one you supplied."
            % (short, os.path.basename(targets_path), len(out)))
    for key in ambiguous:
        out[key] = _AMBIGUOUS
        say("haplotype-scores: WARNING %s:%s%s is named more than once in --targets with "
            "different reported haplotypes; reconciliation for it is "
            "'unknown(ambiguous-in-targets)'"
            % (key[0], key[1], (" guide %s" % key[2]) if key[2] else ""))
    return out


# =========================================================================== #
# ONE IDENTITY. Built in one place, used by every set, lookup and dict.
#
# This file previously held SIX parallel notions of "the same site": site_key, an ad-hoc
# (chrom, pos) tuple for `wanted`, another for `seen_anchors`, another for
# `ref_only_anchors`, a (chrom, pos, spacer) tuple for `reported`, and `_acct_key`.
# EVERY --targets defect found in five review rounds was two of them disagreeing:
#   - `reported` kept the guide while the accounting dropped it, so one guide's rows
#     marked a coordinate "seen" for every other guide curated at it;
#   - `wanted` used the anchor while the filter used the block start, so the
#     CRISPR-Bulge observation was discarded;
#   - the ALT-vs-REF check fell back to "any curated row here", so its verdict flipped
#     with the order of lines in the table.
# None of those is subtle once the key is ONE object. They were only possible because it
# was six.
#
# There are genuinely TWO keys, because the curated table has no strand column and the
# sites dict needs one. They are separate TYPES with one documented projection between
# them, not six tuples that happen to look alike.

TargetKey = collections.namedtuple("TargetKey", "chrom start spacer")


def target_key(chrom, start, spacer=""):
    """THE key for anything that talks to the curated table. Built nowhere else."""
    return TargetKey(str(chrom or "").strip(), str(start or "").strip(),
                     str(spacer or "").strip())


def coarse_key(k):
    """The guide-blind (chrom, start) projection, used ONLY for the read filter.

    Deliberately coarse: a row is read if its coordinate was asked for at all, and the
    guide then decides what is SAID about it. Narrowing this would silently DROP rows
    whose guide spelling differs from the table's, which is a warning (guide_mismatch),
    not a drop.
    """
    return (k.chrom, k.start)


def obs_target_key(reported, obs):
    """This observation's key, guide-qualified when the table named that (locus, guide).

    The one place an Observation becomes a TargetKey. A table with no guide column
    resolves to the guide-less key, so it behaves exactly as it always did.
    """
    exact = target_key(obs.chrom, obs.anchor_start, obs.spacer)
    return exact if exact in reported else target_key(obs.chrom, obs.anchor_start)


def lookup_reported(reported, chrom, start, spacer):
    """The curated row for this site, matching on the guide when the table named one."""
    if not reported:
        return None
    rep = reported.get(target_key(chrom, start, spacer))
    if rep is not None:
        return rep
    return reported.get(target_key(chrom, start))


def build(results_dir, out_path, cap=0, strict=False,
          anonymise=False, pad_up=3, pad_down=6, log=None, targets=None):
    """Write the companion. Returns ``(n_rows, n_sites, unreconciled)``."""
    say = log or (lambda m: print(m, file=sys.stderr))
    if cap is not None and cap < 0:
        # -1 is CRISPRme's own "off" (--max-total-edits, --max-degeneracy), so a user
        # asking for no cap got ranked[:-1] -- every site truncated by one row and a
        # single-row site emptied, silently, exit 0. Same meaning as 0 here, and said.
        say("haplotype-scores: --max-haplotypes %d means NO cap (as elsewhere in "
            "CRISPRme); writing every haplotype." % cap)
        cap = 0
    budget = read_budget(results_dir)
    source_path, source_kind = probe_source(results_dir, say)
    reported = read_reported(targets, say) if targets else {}
    if targets and not reported:
        # Branching on the truthiness of `reported` made a table that parsed to zero
        # usable rows indistinguishable from "no --targets given": the filter switched
        # off and the WHOLE run was written, with none of the per-site accounting the
        # flag exists for, and nothing said.
        raise HaplotypeScoreError(
            "%s yielded no usable rows: it has the required columns but every data row "
            "was unreadable or empty. Writing the whole run instead would silently give "
            "you the opposite of the shortlist you asked for. Refusing." % targets)
    ancestry = read_ancestry(results_dir)
    say("haplotype-scores: ancestry for %s individual(s) from the run's .sampleID.txt"
        % (len(ancestry) if ancestry else "no"))
    if targets:
        say("haplotype-scores: reported alignments for %d site(s) from %s"
            % (len(reported), os.path.basename(targets)))
    else:
        say("haplotype-scores: no --targets given, so row_reported_alignment_found "
            "is 'unknown' on every row (never 'yes')")
    say("haplotype-scores: source %s (%s)" % (source_kind, os.path.basename(source_path)))
    say("haplotype-scores: budget %s from the run's .Params.txt" % budget_text(budget))
    if budget.get("_note"):
        say("haplotype-scores: NOTE %s" % budget["_note"])

    threshold, recall = _threshold()
    sites = collections.OrderedDict()
    # THE ACCOUNTING IS PER (LOCUS, GUIDE), NOT PER LOCUS.
    #
    # read_reported and site_key have been per-guide since the multi-guide collision was
    # fixed, but these sets were left keyed on (chrom, start) -- so one guide's rows
    # marked the coordinate "seen" for EVERY other guide curated at it, and a requested
    # (locus, guide) that produced nothing was silently absent from a file whose header
    # promises such sites are counted and named. It also made the ALT-vs-REF check below
    # pick an arbitrary curated row per coordinate, so that verdict flipped with the
    # ORDER OF LINES in the curated table.
    #
    # The FILTER stays coarse on (chrom, start): a row is read if its coordinate was
    # asked for at all, and the guide decides what is SAID about it. Narrowing the filter
    # would silently drop rows whose guide spelling differs, which is a warning
    # (guide_mismatch), not a drop.
    # Keyed on whether --targets was GIVEN, not on whether it yielded anything.
    wanted = set(reported) if targets else None
    wanted_coarse = set(coarse_key(k) for k in wanted) if targets else None
    seen_anchors = set()        # anchors that produced an ALT haplotype row
    ref_only_anchors = set()    # anchors the source DOES hold, but only as reference
    guide_mismatch = set()      # the locus is in --targets, but under a different guide
    # Guides seen per coordinate. A curated table with NO guide column cannot say which
    # guide its row is about; if the source carries more than one there, every answer
    # but "unknown" is a claim about a guide the table never described.
    guides_at = collections.defaultdict(set)
    mismatched_coarse_all = set()
    # PER LOCUS, NOT PER TABLE. A run-wide any() meant ONE curated row with a blank
    # guide cell, in a table that names guides everywhere else, switched the
    # guide-ambiguity guard off for the WHOLE run -- so the file answered yes/no about
    # a guide the table never described, and --strict failed on it.
    guided_coarse = set(coarse_key(k) for k in (reported or {}) if k.spacer)
    if wanted is not None:
        say("haplotype-scores: restricted to the %d site(s) --targets names" % len(wanted))

    with open_table(source_path) as stream:
        for obs in harvest_rows(stream):
            is_ref = bool(obs.origin) and obs.origin.lower() == "ref"
            obs_key = obs_target_key(reported, obs) if wanted is not None else None
            if (wanted is not None and coarse_key(obs_key) in wanted_coarse
                    and obs_key not in wanted):
                # RECORDED FOR REF ROWS TOO. This lived inside the `site is None`
                # branch, which a reference observation never reaches (they are
                # filtered out just below), so a REF-only locus whose curated guide is
                # spelled differently was reported as "not in <source> at that
                # coordinate" -- false; it is in the source, under another spelling.
                guide_mismatch.add(target_key(obs.chrom, obs.anchor_start, obs.spacer))
                mismatched_coarse_all.add(coarse_key(obs_key))
            if wanted is not None and is_ref and coarse_key(obs_key) in wanted_coarse:
                # RECORDED BEFORE THE REF FILTER. A curated table is mostly REF-origin
                # (115 of 187 sites on the delivering run), and those sites are PRESENT in
                # the source -- they simply carry no ALT haplotype. Dropping them before
                # the anchor was noted made the gap warning call them "not enumerated
                # here", which is false and would have mislabelled most of a real list.
                ref_only_anchors.add(obs_key)
            if is_ref:
                continue                       # reference rows are not variant haplotypes
            if wanted is not None and coarse_key(obs_key) not in wanted_coarse:
                # Bounds memory to the sites actually asked for. Without --targets every
                # site of the run is held at once, which a genome-wide sidecar (millions
                # of rows x three blocks) will not fit in. Streaming with cluster close-out
                # is the next step; see the note in the PR.
                continue
            if wanted is not None:
                seen_anchors.add(obs_key)
                guides_at[coarse_key(obs_key)].add(obs.spacer)
            site = sites.get(obs.site)
            if site is None:
                site = Site(obs.site, obs.chrom, obs.start, obs.strand, obs.spacer)
                site.anchor = obs.anchor_start
                rep = lookup_reported(reported, obs.chrom, obs.anchor_start, obs.spacer)
                if rep is None and reported:
                    # The locus passed the coarse filter, so --targets DOES name it -- but
                    # no entry matched this guide. A curated table whose guide column is
                    # crRNA/Guide holds the 20-mer without the PAM while the source's
                    # Spacer+PAM includes it, so an innocent column-name choice silently
                    # turned reconciliation off for the whole run: every row 'unknown',
                    # no failures, nothing said. Say it.
                    guide_mismatch.add(target_key(obs.chrom, obs.anchor_start,
                                                  obs.spacer))
                site.has_report = rep is not None and rep is not _AMBIGUOUS
                if rep is _AMBIGUOUS:
                    site.ambiguous = True
                if site.has_report and not rep.origin_known:
                    site.unreadable_report = ("origin-not-in-targets"
                                              if rep.origin == _ORIGIN_NO_COLUMN
                                              else "unreadable-origin-in-targets")
                elif site.has_report and rep.origin == "alt" and not rep.variants:
                    # The table names an ALT site whose Variant cell this file cannot
                    # read (an rsID-only row, which generate_report._format_variant_id
                    # does emit as a documented fallback). The check APPLIES and could
                    # not be run; that is not "no ALT to check".
                    site.unreadable_report = "unparsable-variant-in-targets"
                if site.has_report:
                    # REF-origin: the representative IS the reference sequence, so it has
                    # no ALT haplotype to appear among. Not applicable, not a failure.
                    site.applicable = (rep.origin == "alt") and bool(rep.variants)
                    site.reported_alignments = {rep.variants} if site.applicable else set()
                sites[obs.site] = site
            cfd, cb, status, cb_refused = classify(obs, budget, obs.block)
            _cb = carrier_breakdown(obs.samples, ancestry)
            # NEVER ASSERT "no" ABOUT A ROW NOBODY SCORED. Two of the three frames never
            # carry a CRISPR-Bulge score at all, so a bare "no" here told a reader the
            # haplotype is below the operating point when the truth is that it was never
            # measured -- reviving, inside this file, the exact silent discard the file
            # exists to stop. Three values, like every other absence here.
            above = ("yes" if (cb is not None and cb >= threshold)
                     else ("no" if cb is not None
                           else ("unknown(model-refused)" if cb_refused
                                 else ("unknown(not-computed-for-this-pam)"
                                       if status.startswith("not-scored:") else
                                       "unknown"))))
            row = {
                "site_key": obs.site, "chrom": obs.chrom, "site_start": obs.start,
                "strand": obs.strand, "spacer": obs.spacer,
                "haplotype_id": "%s|%s" % (obs.block, obs.variants or "-"),
                "variant_set": obs.variants or "-",
                "aligned_spacer_pam": obs.aligned_spacer,
                "aligned_protospacer_pam": obs.aligned_alt,
                "alignment_frame": FRAME_CANDIDATE,
                "cfd": fmt_score(cfd), "crispr_bulge": fmt_score(cb),
                "score_status": _verdict("score_status", status),
                "above_threshold": _verdict("above_threshold", above),
                "rank_in_site_by_score": "", "rank_score_axis": "",
                "n_carriers": str(_cb[0]),
                "carriers_homozygous": _cb[1],
                "carriers_heterozygous": _cb[2],
                "carrier_panels": _cb[3],
                "carrier_populations": _cb[4],
                "carrier_superpopulations": _cb[5],
                # NEVER A BLANK -- including under --anonymise. A blank cell is an
                # absence with no stated reason, which is the one thing this file does
                # not do; it also failed the project's own gate
                # (scripts/check_haplotype_scores.py) on every row of an anonymised run.
                # The counts beside it are still real, so the cell says WHY it is empty.
                "carriers": ANONYMISED if anonymise else (obs.samples or "-"),
                "variants_in_window_n": str(len([
                    v for v in obs.variants.split(",")
                    if v.strip() and v.strip().upper() not in ("NA", "N", "-", ".")])),
                "max_af": _max_af(obs.maf),
                "scope_reason": source_kind,
            }
            hap_key = obs.variants or obs.aligned_alt or "-"
            # Carried, not reconstructed. _hap() used to rebuild this from the WRITTEN
            # variant_set, which is a lossy rendering of it ("-" means both "no variants"
            # and "the source wrote a literal dash"), so the cap's idea of a haplotype
            # disagreed with the counters' and a site reported 1 haplotype, wrote 1, and
            # said 1 was hidden. extrasaction="ignore" keeps this out of the file.
            row["_hap_key"] = hap_key
            site.add(row, cb, above.startswith("yes"), status.startswith(ST_SCORED),
                     canonical_variants(obs.variants), hap_key=hap_key)

    # A --targets SITE THAT PRODUCED NOTHING IS A GAP, AND GAPS ARE SAID OUT LOUD.
    #
    # The curated table's Position is the merged representative's start; the haplotype
    # source is written pre-merge, so the two agree on most sites but not all. Measured
    # on a real NRR run: of 200 curated top-CFD ALT sites, 97 appear at that exact start
    # in the source, 14 more appear only at another block's start, 47 only within 5 bp,
    # and 42 are absent entirely. Silently returning 97 lets a reader believe the other
    # 103 were examined and found clean.
    #
    # They are NOT fuzzy-matched to recover the near misses: a 5 bp window would merge
    # genuinely distinct sites and manufacture identity, which is worse than a named gap.
    if guide_mismatch:
        say("haplotype-scores: WARNING %d site(s) are named in --targets at this "
            "coordinate but under a DIFFERENT guide than the source carries, so their "
            "reconciliation is 'unknown' rather than yes/no. A curated 'crRNA'/'Guide' "
            "column usually holds the spacer WITHOUT the PAM while the source's "
            "Spacer+PAM includes it; if that is the case here, drop the guide column "
            "from --targets or write it as Spacer+PAM. First: %s"
            % (len(guide_mismatch),
               ", ".join("%s:%s (source guide %s)" % m
                         for m in sorted(guide_mismatch)[:3])))

    if wanted is not None:
        ref_coarse = set(coarse_key(k) for k in ref_only_anchors)
        ref_all = set(k for k in (wanted - seen_anchors)
                      if k in ref_only_anchors
                      or (coarse_key(k) in ref_coarse
                          and coarse_key(k) in mismatched_coarse_all))
        # A curated row that says ALT while the source holds only REF rows is NOT "no
        # variant haplotype here". The user's own table says there is one; this source
        # does not have it. Calling that a clean is the worst answer available.
        #
        # The lookup is EXACT on the accounting key. It used to fall back to "any curated
        # row at this coordinate", so at a locus curated under two guides the verdict was
        # whichever row the dict happened to yield -- the same data gave "GAP" or "no
        # variant haplotype here" depending on the ORDER OF LINES in the table.
        contradicted, ref_only, ambiguous_ref, origin_unknown = [], [], [], []
        for k in sorted(ref_all):
            rep = reported.get(k)
            if rep is _AMBIGUOUS:
                # Nothing is known about what this locus reported, so it belongs in
                # neither "gap" nor "clean". A bucket that asserts is worse than one
                # that says it does not know.
                ambiguous_ref.append(k)
            elif rep is not None and not rep.origin_known:
                # THE TABLE NEVER SAID. Saying "the curated table agrees they are
                # REF-origin" about a table with no origin column -- in the same run as
                # the NOTE saying that column is missing -- asserts the one thing that
                # was not read. Not a gap either; its own line.
                origin_unknown.append(k)
            elif rep is not None and rep.origin == "alt":
                contradicted.append(k)
            else:
                ref_only.append(k)
        if origin_unknown:
            say("haplotype-scores: %d of the %d site(s) --targets names are present in "
                "%s as reference rows only. Whether an ALT haplotype was EXPECTED there "
                "was never read -- the curated table does not say its origin -- so this "
                "is neither a gap nor 'no variant haplotype here'. First: %s"
                % (len(origin_unknown), len(wanted), os.path.basename(source_path),
                   ", ".join("%s:%s" % (m[0], m[1]) for m in origin_unknown[:5])))
        if ambiguous_ref:
            say("haplotype-scores: %d of the %d site(s) --targets names carry only "
                "reference rows in %s AND are described more than once, incompatibly, in "
                "the curated table -- so whether an ALT haplotype was expected there is "
                "unknown. Not counted as a gap and NOT as 'no variant haplotype here'. "
                "First: %s"
                % (len(ambiguous_ref), len(wanted), os.path.basename(source_path),
                   ", ".join("%s:%s" % (m[0], m[1]) for m in ambiguous_ref[:5])))
        if ref_only:
            say("haplotype-scores: %d of the %d site(s) --targets names are PRESENT in "
                "%s but carry only reference rows, and the curated table agrees they are "
                "REF-origin, so they have no ALT haplotype to list. That is 'no variant "
                "haplotype here', NOT 'not enumerated'. First: %s"
                % (len(ref_only), len(wanted), os.path.basename(source_path),
                   ", ".join("%s:%s" % (m[0], m[1]) for m in ref_only[:5])))
        if contradicted:
            say("haplotype-scores: WARNING %d of the %d site(s) --targets names are "
                "ALT-origin in the curated table but carry ONLY reference rows in %s. "
                "The ALT haplotype the curated row reports is NOT in this source -- this "
                "is a GAP, not 'no variant haplotype here'. First: %s"
                % (len(contradicted), len(wanted), os.path.basename(source_path),
                   ", ".join("%s:%s%s" % (m[0], m[1], (" guide %s" % m[2]) if m[2] else "")
                             for m in contradicted[:5])))
        # A LOCUS WHOSE COORDINATE WAS REACHED IS NOT ABSENT, whatever guide reached
        # it. When the curated table spells the guide differently -- a crRNA column
        # holding the 20-mer while the source's Spacer+PAM includes it, the documented
        # common case -- obs_target_key falls back to the guide-less key, so the exact
        # curated key never lands in seen_anchors. Every such site was then reported as
        # "produced NO rows and is absent from this file" IN THE SAME RUN as the warning
        # saying its rows were matched under a different guide, while its rows were in
        # the file. Two contradictory statements about one site; the second is false.
        # The suppression is narrow ON PURPOSE. Two cases look alike and are not:
        #   (a) the table names guide B at a locus where only guide A was enumerated.
        #       Guide B genuinely produced nothing: that IS a gap and must be said.
        #   (b) the table spells the guide differently from the source, so NO curated
        #       key matched and every observation fell back to the guide-less key. The
        #       rows ARE in the file; calling the locus absent contradicts them.
        # Only (b) is suppressed, and it is identified by guide_mismatch -- which is
        # populated exactly when a locus was reached but no curated guide matched it.
        # SUPPRESS ONLY WHERE NO CURATED GUIDE MATCHED ANYTHING AT ALL.
        #
        # The previous rule dropped every curated key sharing a coordinate with a
        # mismatch, so at a locus where guide A matched exactly and guide B genuinely
        # produced nothing, B was silently removed from the gap list -- the very case
        # the rule's own comment promised was still reported. The two situations differ
        # by whether ANY curated key at that coordinate was matched:
        #   matched there  -> the table does describe this locus, so a curated guide
        #                     that produced nothing really is a gap. Say it.
        #   nothing matched -> no curated key resolved; the rows are in the file under a
        #                     spelling the table does not use. That is the mismatch
        #                     warning's business, not a gap.
        mismatched_coarse = set(coarse_key(k) for k in guide_mismatch)
        matched_coarse = set(coarse_key(k) for k in seen_anchors if k in reported)
        mismatched_coarse -= matched_coarse
        missing = sorted(k for k in (wanted - seen_anchors - ref_all)
                         if coarse_key(k) not in mismatched_coarse)
        if missing:
            say("haplotype-scores: WARNING %d of the %d site(s) --targets names produced "
                "NO rows and are absent from this file: they are not in %s at that "
                "coordinate. This is not 'no haplotypes there' -- it is 'not enumerated "
                "here'. The curated Position is the merged representative's start while "
                "this source is written pre-merge, so some sites sit a few bp away, and "
                "they are NOT fuzzy-matched because a window would merge distinct sites. "
                "First missing: %s"
                % (len(missing), len(wanted), os.path.basename(source_path),
                   ", ".join("%s:%s%s" % (m[0], m[1], (" guide %s" % m[2]) if m[2] else "")
                             for m in missing[:5])))

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
        say("haplotype-scores: WARNING the curated table's Variant column and this "
            "source's Variant_info_genome reduce to the same records on only %d of %d "
            "site(s). That LOOKS like a mis-keyed comparison rather than %d real "
            "failures -- check the chromosome prefix, the coordinate base and the "
            "REF/ALT order in %s. The sites that did NOT match are reported as "
            "row_reported_alignment_found=unknown rather than a wall of 'no'; the ones "
            "that DID match keep their answer. If the key is right and these really are "
            "failures, that is what the per-site lines above say."
            % (len(matched), len(checked), len(checked) - len(matched),
               os.path.basename(targets) if targets else "the curated table"))

    n_rows = 0
    unreconciled = []
    with io.open(out_path, "w", encoding="utf-8", newline="") as out:
        out.write(file_header(source_path, source_kind, budget, threshold, recall,
                              cap, strict))
        # THE WRITER MATCHES THE READER. With the default dialect a cell holding a
        # double quote came back out as ""a""b,c"" -- a file written as CSV while every
        # reader in this project, including this file's own, parses it as a TSV. A TSV
        # has no quoting convention: quotechar=None says so, and no cell can hold a tab
        # because the source is tab-separated.
        writer = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore",
                                quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        for site in sites.values():
            def _rank_key(r):
                """Rank WITHIN an axis. The two scores are never compared.

                CFD and CRISPR-Bulge live on different scales -- CFD runs 0..1 while the
                CRISPR-Bulge operating point is 1.9e-4 -- and only the CRISPR-Bulge block
                carries a CB column. Ordering them in one sequence put every CFD-only row
                above every CB row, so cap=1 on a site with CFD 0.6 and CB 0.004 (twenty
                times the operating point) wrote the CFD row and DISCARDED the only
                haplotype that was above the threshold, while the same line still said
                n_above_threshold=1. A false clean produced by the cap.

                So: CB-bearing rows first, ordered by CB; then CFD-only rows, ordered by
                CFD; then the unscored. A cap takes from the bottom, and the rows this
                file exists to surface are never the ones it takes.
                """
                def val(col):
                    if r[col] in ("NA", ""):
                        return None
                    try:
                        return float(r[col])
                    except ValueError:
                        return None
                cb, cfd = val("crispr_bulge"), val("cfd")
                if cb is not None:
                    return (0, -cb, r["variant_set"])
                if cfd is not None:
                    return (1, -cfd, r["variant_set"])
                return (2, 0.0, r["variant_set"])

            def _rank_axis(r):
                if r["crispr_bulge"] not in ("NA", ""):
                    return "crispr_bulge"
                return "cfd" if r["cfd"] not in ("NA", "") else "none"

            ranked = sorted(site.rows, key=_rank_key)
            kept = ranked if not cap else ranked[:cap]
            # HAPLOTYPES, not rows. One haplotype appears once per frame that named it, so
            # counting dropped ROWS reported "2 haplotype(s) not written" for a site whose
            # own rows said n_haplotypes_at_site=1.
            def _hap(r):
                """The identity Site.add counted on, carried on the row itself."""
                return r["_hap_key"]
            hidden_haplotypes = len(set(_hap(r) for r in ranked)
                                    - set(_hap(r) for r in kept))
            hidden_rows = len(ranked) - len(kept)
            # The reconciliation speaks about the haplotypes THIS FILE SHOWS. If a cap
            # removed some, a "yes" earned by a hidden row would be a claim the reader
            # cannot check against anything in front of them.
            # Degrade ONLY the sites that did not match. Blanking the matched ones too
            # made "fully reconstructed" and "nobody could check this" read identically
            # -- the collapse this file exists to prevent, applied to its own column.
            found = site.row_reported_alignment_found
            # A GUIDE-LESS CURATED TABLE CANNOT SPEAK FOR TWO GUIDES. The shipped
            # curated shape has no guide column, so when the source carries several
            # guides at one coordinate, lookup_reported hands the SAME curated row to
            # all of them and each answers yes/no -- an answer about a guide the table
            # never described. One of them is necessarily wrong, and nothing said so.
            # n/a IS ALSO AN ANSWER. It asserts "the curated table says this site
            # reported a REFERENCE haplotype" -- about a guide the table may never have
            # named -- so it needs the same guard as yes/no.
            _coarse = (site.chrom, site.anchor)
            if (found in ("yes", "no", "n/a") and _coarse not in guided_coarse
                    and len(guides_at.get(_coarse, ())) > 1):
                found = "unknown(guide-ambiguous-in-targets)"
            if mis_keyed and found == "no":
                found = "unknown(key-looks-wrong)"
            # Only downgrade when the matching haplotype is actually out of sight. A cap
            # that hid other FRAMES of a shown haplotype hides nothing the reader needs.
            if found == "yes" and hidden_haplotypes:
                shown = set(canonical_variants(r["variant_set"]) for r in kept)
                if not (site.reported_alignments & shown):
                    found = "yes(among-all-haplotypes-not-all-shown)"
            if found == "no":
                unreconciled.append(site.key)
            for rank, row in enumerate(kept, 1):
                row["rank_in_site_by_score"] = str(rank)
                row["rank_score_axis"] = _rank_axis(row)
                row["n_haplotypes_at_site"] = str(site.n_haplotypes)   # the TRUE total
                row["n_above_threshold"] = str(site.n_above)
                row["n_not_scored"] = str(site.n_not_scored)
                row["site_max_crispr_bulge"] = fmt_score(site.max_cb)
                row["row_reported_alignment_found"] = _verdict(
                    "row_reported_alignment_found", found)
                writer.writerow(row)
                n_rows += 1
            if hidden_rows:
                say("haplotype-scores: %s: %d row(s) covering %d haplotype(s) not "
                    "written (cap=%d counts ROWS); the counts in the file are the true "
                    "totals" % (site.key, hidden_rows, hidden_haplotypes, cap))

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
                   help="cap rows written per site (0 or a negative number = every one, "
                        "matching CRISPRme's --max-total-edits/-1 'off' convention). "
                        "Counts are always the true totals.")
    p.add_argument("--targets", default=None,
                   help="the run's curated results table (Chromosome, Position, "
                        "Aligned_protospacer+PAM), so each site's own reported alignment "
                        "can be reconciled against the haplotypes written. Without it "
                        "row_reported_alignment_found is 'unknown', never 'yes'.")
    p.add_argument("--strict", action="store_true",
                   help="exit non-zero if any site's reported alignment is on none of "
                        "its written haplotypes (default: warn loudly, exit 0).")
    p.add_argument("--anonymise", "--anonymize", action="store_true", dest="anonymise",
                   help="write carrier counts but no sample identities: the carriers "
                        "column reads 'anonymised', never blank, so the absence has a "
                        "stated reason like every other absence in the file.")
    args = p.parse_args(argv)

    try:
        n_rows, n_sites, unreconciled = build(
            args.results, args.out, cap=args.cap,
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
