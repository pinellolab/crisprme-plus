#!/usr/bin/env python3
"""Assert the per-haplotype companion's own invariants against a REAL file.

WHY THIS EXISTS. Five review rounds of ``scripts/haplotype_scores.py`` found the same
class of defect repeatedly -- two genuinely different absences rendered as the same cell,
or a counter contradicting the rows beside it -- and every one was found by a person
reading code. The unit suite now enumerates the input space, but synthetic fixtures were
twice wrong about the real schema, and a fixture that is not the real shape is not
coverage. The strongest evidence available at any point came from running on a real run
and checking the output: that is what this does, mechanically, on any results folder.

It reads a companion TSV (and optionally the run it came from) and checks ONLY the things
the file promises about itself. It computes no biology and makes no judgement about
whether a score is right -- it asserts that the file cannot be MISREAD.

    python3 scripts/check_haplotype_scores.py --tsv <companion>.tsv
    python3 scripts/check_haplotype_scores.py --results Results/<job> --tsv out.tsv

Exit 0 when every invariant holds, 1 otherwise, with each violation named and counted and
a sample of offending rows. Intended for the release checklist and for anyone who wants to
trust a companion they did not generate.
"""

from __future__ import print_function

import argparse
import collections
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


class Violation(object):
    """One failed invariant, with examples kept for the report."""

    __slots__ = ("name", "why", "rows", "n")

    def __init__(self, name, why):
        self.name, self.why, self.rows, self.n = name, why, [], 0

    def hit(self, row):
        self.n += 1
        if len(self.rows) < 3:
            self.rows.append(row)


def _read(path):
    """Header comments and body rows of a companion TSV."""
    comments, header, body = [], None, []
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                comments.append(line)
                continue
            if header is None:
                header = line.split("\t")
                continue
            if line:
                body.append(dict(zip(header, line.split("\t"))))
    return comments, (header or []), body


def check(path, threshold=None):
    """Every invariant the file states about itself. Returns a list of Violations."""
    comments, header, body = _read(path)
    v = collections.OrderedDict()

    def mark(key, why):
        if key not in v:
            v[key] = Violation(key, why)
        return v[key]

    if not header:
        mark("no-header", "the file has no column header").hit({})
        return list(v.values())
    if not body:
        mark("no-rows", "the file has a header and no rows: that is 'never enumerated' "
                        "wearing the clothes of 'nothing found'").hit({})
        return list(v.values())

    # The threshold the file itself declares, so this script needs no module import.
    #
    # Read from the WHOLE comment block, not line by line. The real header writes
    # "...the CRISPR-Bulge sensitive operating point," and the number on the NEXT line,
    # so a per-line match found nothing on every real file -- and the check below was
    # silently skipped, which is how a gate passes 460,714 rows while inert. A check
    # that cannot run now says so instead of reporting success.
    if threshold is None:
        blob = " ".join(c.lstrip("# ").strip() for c in comments)
        anchor = blob.find("operating point")
        if anchor >= 0:
            for tok in blob[anchor:anchor + 200].replace(",", " ").split():
                try:
                    val = float(tok)
                except ValueError:
                    continue
                if 0.0 < val < 1.0:
                    threshold = val
                    break
    if threshold is None:
        mark("threshold-not-stated",
             "the file's header does not state its operating point, so "
             "above_threshold cannot be checked against the scores beside it. "
             "Pass --threshold to check it anyway; reporting success without this "
             "check would overstate what was verified").hit({})

    declared = set()
    for c in comments:
        s = c.lstrip("# ").strip()
        if s.startswith(("yes", "no", "n/a", "unknown", "scored", "out-of-domain",
                         "needs-realignment")) and "--" in s:
            declared.add(s.split("--")[0].strip())

    per_site = collections.defaultdict(list)
    for r in body:
        per_site[r.get("site_key", "")].append(r)

    for r in body:
        # 1. NEVER A BLANK. The file's first promise.
        blanks = [k for k, val in r.items() if val == ""]
        if blanks:
            mark("blank-cell",
                 "a blank cell: an absence with no stated reason").hit(r)

        # 2. The -1 sentinel is never printed and never sortable.
        #    Compared NUMERICALLY. This tested the literal "-1.0" while %.6g renders
        #    -1.0 as "-1" and the pipeline writes "-1.000", so it could not fire on any
        #    spelling anything actually produces -- the third check found here that
        #    reported success for work it did not do.
        for _col in ("cfd", "crispr_bulge"):
            _raw = r.get(_col, "NA")
            if _raw in ("NA", "", "-", None):
                continue
            try:
                _val = float(_raw)
            except ValueError:
                mark("score-not-a-number",
                     "a score column holds something that is neither a number nor a "
                     "stated absence, so no filter can treat it correctly").hit(r)
                continue
            if _val <= -1.0:
                mark("sentinel-printed",
                     "the model's refusal sentinel reached the file, where a numeric "
                     "filter reads it as a very low score").hit(r)

        # 3. A row is never simultaneously above the threshold and never scored.
        if r.get("above_threshold") == "yes" and not str(
                r.get("score_status", "")).startswith("scored"):
            mark("above-but-unscored",
                 "above_threshold=yes on a row whose status says it was not "
                 "scored -- a reader filtering on either one is misled").hit(r)

        # 4. above_threshold agrees with the number printed beside it.
        cb = r.get("crispr_bulge", "NA")
        at = r.get("above_threshold", "")
        if cb in ("NA", ""):
            if not at.startswith("unknown"):
                mark("verdict-without-number",
                     "above_threshold claims a measurement with no number beside "
                     "it").hit(r)
        elif threshold is not None:
            try:
                want = "yes" if float(cb) >= threshold else "no"
            except ValueError:
                want = None
            if want and at != want:
                mark("verdict-disagrees",
                     "above_threshold contradicts the score printed beside it "
                     "(threshold %g)" % threshold).hit(r)

        # 5. Scores are %.6g, never rounded to a width that destroys the threshold.
        for col in ("cfd", "crispr_bulge"):
            val = r.get(col, "NA")
            if val not in ("NA", "") and val.count(".") == 1:
                frac = val.split(".")[1]
                if len(frac) == 3 and "e" not in val and float(val) == 0.0:
                    mark("score-rounded-away",
                         "a score printed at 3 dp, where the 1.9e-4 operating point "
                         "itself rounds to 0.000").hit(r)

        # 6. Every verdict the file writes is one the file's own header explains.
        if declared:
            for col in ("score_status", "above_threshold",
                        "row_reported_alignment_found"):
                val = r.get(col)
                if val is None:
                    continue
                if val in declared:
                    continue
                if any(val.startswith(d.split("<")[0]) for d in declared if "<" in d):
                    continue
                if any(val.startswith(d.split(":")[0] + ":") for d in declared
                       if ":" in d):
                    continue
                mark("undeclared-verdict",
                     "a value the file's own header does not explain, so a reader "
                     "has nothing to interpret it with").hit(r)

    # 7b. The per-site CRISPR-Bulge roll-up must not contradict the rows under it.
    #     A cap only ever REMOVES rows, so site_max_crispr_bulge can never read NA while
    #     a row shows a number, and can never be below the largest shown. This was
    #     checked NOWHERE: setting it to NA on a site whose rows carry 0.004 -- twenty
    #     times the operating point -- passed with "OK: every invariant holds".
    for key, rows in per_site.items():
        nums = []
        for r in rows:
            try:
                nums.append(float(r.get("crispr_bulge", "NA")))
            except ValueError:
                pass
        if not nums:
            continue
        stated = rows[0].get("site_max_crispr_bulge", "NA")
        try:
            stated_v = float(stated)
        except ValueError:
            mark("site-max-missing",
                 "site_max_crispr_bulge says nothing was measured on a site whose own "
                 "rows carry a CRISPR-Bulge score").hit(rows[0])
            continue
        if stated_v < max(nums):
            mark("site-max-below-rows",
                 "site_max_crispr_bulge is below the largest score on the site's own "
                 "rows; a cap only removes rows, so the roll-up cannot be lower"
                 ).hit(rows[0])

    # 7. COUNT BEFORE CAP: the site's stated totals must not contradict its rows.
    for key, rows in per_site.items():
        shown = set(r.get("variant_set", "") for r in rows)
        try:
            total = int(rows[0].get("n_haplotypes_at_site", "0"))
            above = int(rows[0].get("n_above_threshold", "0"))
            unscored = int(rows[0].get("n_not_scored", "0"))
        except ValueError:
            mark("counter-not-a-number",
                 "a per-site counter is not a number").hit(rows[0])
            continue
        if total < len(shown):
            mark("total-below-shown",
                 "n_haplotypes_at_site is smaller than the number of distinct "
                 "haplotypes actually written for that site").hit(rows[0])
        if above > total or unscored > total:
            mark("counter-exceeds-total",
                 "a per-site counter exceeds the site's own total").hit(rows[0])
        # PER HAPLOTYPE, NOT PER ROW. One haplotype appears once per alignment frame
        # that named it, and can appear twice on ONE frame under different alignments --
        # measured on a real run: a site with n_above_threshold=2 legitimately shows 3
        # rows marked yes, because one variant set is there twice. Counting rows here
        # made this checker fail a correct file on its first run against real data.
        shown_above = len(set(r.get("variant_set", "") for r in rows
                              if r.get("above_threshold") == "yes"))
        if shown_above > above:
            mark("above-exceeds-counter",
                 "more distinct haplotypes are marked above_threshold than "
                 "n_above_threshold claims the site has").hit(rows[0])
        # every row of a site must agree about that site
        for col in ("n_haplotypes_at_site", "n_above_threshold", "n_not_scored",
                    "site_max_crispr_bulge", "row_reported_alignment_found"):
            if len(set(r.get(col) for r in rows)) > 1:
                mark("site-disagrees-with-itself",
                     "rows of one site carry different values in %s, so the "
                     "'one row is self-contained' promise is false" % col).hit(rows[0])

    return list(v.values())


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Assert the per-haplotype companion's own invariants on a real file.")
    ap.add_argument("--tsv", required=True, help="the companion .tsv to check")
    ap.add_argument("--threshold", type=float, default=None,
                    help="the operating point, if the file's header does not state it")
    ap.add_argument("--quiet", action="store_true", help="print only the verdict")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.tsv):
        print("no such file: %s" % args.tsv, file=sys.stderr)
        return 2

    violations = check(args.tsv, args.threshold)
    _c, header, body = _read(args.tsv)
    if not args.quiet:
        print("checked %s: %d row(s), %d column(s)" % (args.tsv, len(body), len(header)))
    if not violations:
        print("OK: every invariant holds.")
        return 0
    print("FAILED: %d invariant(s) violated." % len(violations))
    for v in violations:
        print("\n  [%s] %d row(s)\n    %s" % (v.name, v.n, v.why))
        for r in v.rows:
            bits = [(k, r.get(k)) for k in ("site_key", "variant_set", "cfd",
                                            "crispr_bulge", "score_status",
                                            "above_threshold") if k in r]
            print("      " + "  ".join("%s=%s" % kv for kv in bits))
    return 1


if __name__ == "__main__":
    sys.exit(main())
