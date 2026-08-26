#!/usr/bin/env python
"""Regression: off-target dedup columns are emitted in a DETERMINISTIC order.

The SNP-info / rsID / AF / Samples / Annotation columns are deduped by joining a
Python ``set`` back into a comma string. A bare ``set`` iterates in
``PYTHONHASHSEED``-dependent order, so identical inputs produced byte-different
output run-to-run (issue #46). The fix wraps each dedup in ``sorted(...)``.

merge_contiguous_targets.remove_duplicate_targets is importable (its main is
__main__-guarded); it exercises the exact join pattern shared by
remove_contiguous_samples.py and annotate_final_results.py.

Run: python PostProcess/test_determinism.py   (or via unittest/pytest)
"""

import merge_contiguous_targets as mct


def test_remove_duplicate_targets_sorted_and_deduped():
    # cols: 0=snp_info, 1=rsID, 2=AF, 3=samples -- each has dupes + unordered
    target = ["6G>C,3C>A,6G>C", "rs2,rs1,rs2", "0.2,0.1,0.2", "HG3,HG1,HG2,HG1"]
    out = mct.remove_duplicate_targets(list(target), 0, 1, 2, 3)
    assert out[0] == "3C>A,6G>C", out[0]      # deduped + lexicographically sorted
    assert out[1] == "rs1,rs2", out[1]
    assert out[2] == "0.1,0.2", out[2]
    assert out[3] == "HG1,HG2,HG3", out[3]


def test_output_order_invariant_to_input_permutation():
    # the whole point: the emitted order must NOT depend on input/iteration order
    perms = ["HG2,HG1,HG3", "HG3,HG2,HG1", "HG1,HG3,HG2", "HG2,HG3,HG1"]
    outs = {
        mct.remove_duplicate_targets(["a", "b", "c", p], 0, 1, 2, 3)[3]
        for p in perms
    }
    assert outs == {"HG1,HG2,HG3"}, outs      # one deterministic result for all


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")
