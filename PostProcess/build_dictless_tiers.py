"""Build-time emission of the Tier-0 registry + Tier-1 genotype store (CRISPRme+
dictless redesign, Phase-3d BUILD wiring).

The SEARCH side already CONSUMES these tiers: ``new_simple_analysis.py`` resolves,
per chromosome, a Tier-0 registry (``registry_<vcf>/reg_<chrom>.bin`` + ``.idx``)
and a Tier-1 genotype store (``genotypes_<vcf>/gt_<chrom>.bin`` + ``.idx``) that
live as SIBLING directories of the per-sample SNP dicts
(``dictionaries_<vcf>/my_dict_<chrom>.json[.gz]``). What was MISSING is the BUILD
side that PRODUCES them. This module is that producer, made an IMPORTABLE, unit-
testable helper (mirroring ``simple_analysis_registry`` / ``population_summary_
companion``) so the build path in ``crisprme.py`` calls a single guarded function
per chromosome.

CRITICAL correctness property (the #1 risk): the directory + file names emitted
here MUST EXACTLY match what the search resolvers read, or a freshly-built install
silently ignores the tiers. To guarantee that, ``_sibling_tier_dir`` below derives
the sibling dir the SAME WAY as ``new_simple_analysis._resolve_registry_paths`` /
``_resolve_genotype_paths``: it swaps the ``dictionaries_`` folder-name prefix for
``registry_`` / ``genotypes_`` (falling back to prefixing a non-standard folder
name), and the per-chromosome files are named ``reg_<chrom>.bin``/``.idx`` and
``gt_<chrom>.bin``/``.idx``. The unit test asserts the emitted paths are exactly the
paths the resolvers would return for the same ``dict_path``.

ADDITIVE + GATED + GUARDED (see ``emit_dictless_tiers_guarded``):
  * ADDITIVE -- the per-sample dicts are untouched; the tiers are a NEW sibling
    output. Nothing about ``my_dict_<chrom>.json[.gz]`` changes.
  * GATED    -- absence of the tier modules (an old deploy) or a bad dict simply
    means the tiers are not produced; the dict build proceeds. The search side
    then falls back to the legacy dict path (its resolvers return None on absence).
  * GUARDED  -- ``emit_dictless_tiers_guarded`` wraps the emission in try/except and
    logs a clear WARNING to STDOUT (NEVER stderr -- stderr is fatal in this
    pipeline's ``[ -s $logerror ]`` stage checks) and NEVER raises.

gz DECISION: the emission reads whichever dict exists (plain or ``.json.gz``). Both
compile functions ultimately stream via ``tier0_compile.iter_dict_records``, which
is ALREADY gz-aware (it opens ``.gz`` with ``gzip`` transparently). So we do NOT
need to compile before the build gzips the dict; ``emit_dictless_tiers`` locates the
dict itself (``.json`` preferred, else ``.json.gz``) and hands the found path to the
compilers. This makes the caller order-independent: it can emit before OR after the
in-place gzip step and still get the same result.

STDLIB ONLY (os) plus the two tier modules (import-guarded at the call site).
"""

from __future__ import annotations

import os


def _sibling_tier_dir(dict_path, prefix):
    """Return the sibling tier directory for a dict path, mirroring the resolvers.

    ``dict_path`` is ``<...>/Dictionaries/dictionaries_<vcf>/my_dict_<chrom>.json``
    (or ``.json.gz``). This computes ``<...>/Dictionaries/<prefix><vcf>`` EXACTLY as
    ``new_simple_analysis._resolve_registry_paths`` / ``_resolve_genotype_paths`` do:

      dict_dir         = dirname(dict_path)        -> <...>/Dictionaries/dictionaries_<vcf>
      parent           = dirname(dict_dir)         -> <...>/Dictionaries
      dict_folder_name = basename(dict_dir)        -> dictionaries_<vcf>
      tier_folder_name = <prefix> + <vcf>          (swap the dictionaries_ prefix)
                         else <prefix> + dict_folder_name (non-standard layout)
      -> join(parent, tier_folder_name)

    ``prefix`` is "registry_" or "genotypes_".
    """
    dict_dir = os.path.dirname(dict_path)
    parent = os.path.dirname(dict_dir)
    dict_folder_name = os.path.basename(dict_dir)
    if dict_folder_name.startswith("dictionaries_"):
        tier_folder_name = prefix + dict_folder_name[len("dictionaries_"):]
    else:
        # non-standard layout: mirror the resolvers' fallback (prefix the whole
        # folder name). Only USED because the resolvers look here too.
        tier_folder_name = prefix + dict_folder_name
    return os.path.join(parent, tier_folder_name)


def registry_paths_for(dict_path, chrom):
    """(bin, idx) the Tier-0 registry for ``chrom`` MUST be emitted to.

    Byte-for-byte the paths ``_resolve_registry_paths`` will look for."""
    reg_dir = _sibling_tier_dir(dict_path, "registry_")
    return (os.path.join(reg_dir, "reg_" + str(chrom) + ".bin"),
            os.path.join(reg_dir, "reg_" + str(chrom) + ".idx"))


def genotype_paths_for(dict_path, chrom):
    """(bin, idx) the Tier-1 genotype store for ``chrom`` MUST be emitted to.

    Byte-for-byte the paths ``_resolve_genotype_paths`` will look for."""
    gt_dir = _sibling_tier_dir(dict_path, "genotypes_")
    return (os.path.join(gt_dir, "gt_" + str(chrom) + ".bin"),
            os.path.join(gt_dir, "gt_" + str(chrom) + ".idx"))


def _resolve_dict_path(dict_path):
    """Return an existing dict path (plain preferred, else .gz), or raise.

    ``dict_path`` may be given with or without the ``.gz`` suffix; the build gzips
    the dict in place, so at emission time only ONE of ``my_dict_<chrom>.json`` /
    ``my_dict_<chrom>.json.gz`` exists. We accept either, preferring the plain form
    when both are present (identical content; plain avoids a decompress). Raises
    FileNotFoundError with a clear message if neither exists, so the guarded wrapper
    logs it and the dict build continues.
    """
    if dict_path.endswith(".gz"):
        plain = dict_path[:-len(".gz")]
        candidates = [plain, dict_path]
    else:
        candidates = [dict_path, dict_path + ".gz"]
    for cand in candidates:
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(
        "build_dictless_tiers: no SNP dict found at %s (nor its .gz/plain "
        "variant); cannot emit dictless tiers" % (dict_path,))


def emit_dictless_tiers(dict_path, db_to_samplesid, chrom, dictionaries_dir=None,
                        *, subpop_field="superpopulation"):
    """Emit the Tier-0 registry + Tier-1 genotype store for ONE chromosome.

    Compiles both tiers FROM the per-sample SNP dict ``dict_path`` (the same dict
    the build already produced -- untouched), into the SIBLING ``registry_<vcf>`` /
    ``genotypes_<vcf>`` directories the search resolvers read. The dirs are created
    if absent. Returns a dict::

        {"registry_bin", "registry_idx", "genotype_bin", "genotype_idx",
         "registry_stats", "genotype_stats"}

    Args:
      dict_path: path to ``<...>/dictionaries_<vcf>/my_dict_<chrom>.json`` (or the
        ``.json.gz`` form -- either is accepted; the compilers are gz-aware).
      db_to_samplesid: ORDERED mapping {database_name: samplesID_path}. For the
        combined 1000G+HGDP panel this has BOTH entries (order preserved). Consumed
        by ``tier0_compile.build_sample_meta`` / ``tier1_genotypes`` to build the
        panel + sample axis.
      chrom: chromosome name AS WRITTEN in the dict keys (e.g. "chr1", "chrX"); it
        drives the ploidy model AND the emitted file names (reg_<chrom>/gt_<chrom>).
      dictionaries_dir: unused for path derivation (paths are derived from
        ``dict_path`` to guarantee the resolver match); accepted for call-site
        clarity / forward compatibility. If given, it is only used as a sanity check
        that ``dict_path`` lives under it.
      subpop_field: "superpopulation" (default) or "population".

    Raises on a missing dict / compile error; use ``emit_dictless_tiers_guarded`` at
    a build call site so a failure NEVER aborts the dict build.
    """
    # Import the tier modules lazily so a caller that only wants path helpers, or an
    # environment without them, is not forced to import numpy-free-but-present deps.
    import tier0_compile as t0c
    import tier1_genotypes as t1g

    resolved_dict = _resolve_dict_path(dict_path)

    if dictionaries_dir is not None:
        exp_dir = os.path.abspath(dictionaries_dir)
        got_dir = os.path.abspath(os.path.dirname(resolved_dict))
        if exp_dir != got_dir:
            # Not fatal -- the resolver derivation uses dict_path, not this arg -- but
            # surface a clear STDOUT note so a mis-wired caller is visible.
            print(
                "build_dictless_tiers: NOTE dict %s is not under the passed "
                "dictionaries_dir %s (paths derived from dict_path regardless)"
                % (resolved_dict, dictionaries_dir),
                flush=True,
            )

    # Derive the emit targets EXACTLY as the search resolvers will read them.
    reg_bin, reg_idx = registry_paths_for(dict_path, chrom)
    gt_bin, gt_idx = genotype_paths_for(dict_path, chrom)

    os.makedirs(os.path.dirname(reg_bin), exist_ok=True)
    os.makedirs(os.path.dirname(gt_bin), exist_ok=True)

    # Ship the v3 block-compressed registry (Issue #99): ~3.6x smaller on-disk,
    # logically identical, and the reader is backward-compatible (reads v2 + v3).
    # The build now emits v3 directly instead of requiring a manual post-build
    # transcode_registry pass.
    reg_stats = t0c.compile_from_dict(
        resolved_dict, db_to_samplesid, chrom, reg_bin, reg_idx,
        subpop_field=subpop_field, compress=True,
    )
    gt_stats = t1g.compile_genotypes_from_dict(
        resolved_dict, db_to_samplesid, chrom, gt_bin, gt_idx,
        subpop_field=subpop_field,
    )

    return {
        "registry_bin": reg_bin,
        "registry_idx": reg_idx,
        "genotype_bin": gt_bin,
        "genotype_idx": gt_idx,
        "registry_stats": reg_stats,
        "genotype_stats": gt_stats,
    }


def emit_dictless_tiers_guarded(dict_path, db_to_samplesid, chrom,
                                dictionaries_dir=None, *,
                                subpop_field="superpopulation"):
    """GUARDED ``emit_dictless_tiers``: never raises; logs to STDOUT on failure.

    Returns the paths dict on success, or None on ANY failure (missing dict, missing
    tier modules, compile error). A failure logs a clear WARNING to STDOUT (NEVER
    stderr -- stderr writes fail the pipeline's ``[ -s $logerror ]`` stage checks)
    and leaves the dict build untouched. This is the function build call sites use.
    """
    try:
        return emit_dictless_tiers(
            dict_path, db_to_samplesid, chrom, dictionaries_dir,
            subpop_field=subpop_field,
        )
    except Exception as exc:  # never let tier emission abort the dict build
        print(
            "WARNING [build_dictless_tiers]: could not emit dictless tiers for "
            "%s (chrom %s): %s -- the per-sample dict is unaffected; the search "
            "side will fall back to the legacy dict path for this chromosome."
            % (dict_path, chrom, exc),
            flush=True,
        )
        return None
