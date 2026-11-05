#!/usr/bin/env python
"""Guard: --samplesID over-listing the VCF-genotyped panel is rejected with an ERROR (#46).

Phantom samplesID entries (listed but genotyped in no VCF) are counted as hom-ref
when the Tier-0 panel AN is built, inflating AN and silently DEFLATING every
reported allele frequency (Samples/AF/MAF) -- shared by the classic SNP path and
the indel+SNP joint-AF. validate_inputs.check_samplesid_overlisting surfaces it at
pre-flight, and is correct for multi-dataset runs (union of all VCF headers), a search
that names a dataset with no VCFs or no sample columns here among them.

Run: python PostProcess/test_samplesid_overlisting.py
"""

import gzip
import os
import tempfile

import validate_inputs as vi


def _mkvcf(d, name, samples):
    p = os.path.join(d, name)
    with gzip.open(p, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n")
    return p


def test_overlisting_errors_and_quantifies():
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2", "S3"])
    iss = vi.check_samplesid_overlisting([v], ["S1", "S2", "S3", "PH1", "PH2"])
    assert len(iss) == 1 and iss[0].severity == vi.ERROR, iss
    assert "2 sample" in iss[0].message and "registry_fix_an" in iss[0].message


def test_matched_panel_is_silent():
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2", "S3"])
    assert vi.check_samplesid_overlisting([v], ["S1", "S2", "S3"]) == []


def test_multidataset_union_no_false_positive():
    d = tempfile.mkdtemp()
    v1 = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2", "S3"])
    v2 = _mkvcf(d, "chr2.vcf.gz", ["S4", "S5"])
    # every sample is genotyped in SOME vcf -> no phantom, even though each is
    # absent from the other dataset's header
    assert vi.check_samplesid_overlisting([v1, v2], ["S1", "S2", "S3", "S4", "S5"]) == []


def test_multidataset_real_phantom_flagged():
    d = tempfile.mkdtemp()
    v1 = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2"])
    v2 = _mkvcf(d, "chr2.vcf.gz", ["S3", "S4"])
    iss = vi.check_samplesid_overlisting([v1, v2], ["S1", "S2", "S3", "S4", "GHOST"])
    assert len(iss) == 1 and "GHOST" in iss[0].message


def test_list_of_a_dataset_with_no_vcf_here_is_left_out():
    # one dataset's VCFs on disk, the other's index downloaded: its list names none of
    # the samples of the VCFs here, and goes with the dataset that has none here
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2"])
    lists = [["S1", "S2"], ["B4", "B5"]]
    flat = [s for one in lists for s in one]
    assert vi.check_samplesid_overlisting([v], flat, lists, 1) == []
    assert vi.check_samplesid_overlisting([v], flat, lists[::-1], 1) == []
    # the same two lists with every dataset's VCFs on disk: it belongs to none
    iss = vi.check_samplesid_overlisting([v], flat, lists, 0)
    assert len(iss) == 1 and "2 sample" in iss[0].message and "B4, B5" in iss[0].message, iss


def test_a_list_that_partly_overlaps_the_vcfs_still_fails():
    # (#46) it is the list of the VCFs that are here, with samples they lack
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2"])
    lists = [["S1", "S2", "GHOST"], ["B4", "B5"]]
    flat = [s for one in lists for s in one]
    iss = vi.check_samplesid_overlisting([v], flat, lists, 1)
    assert len(iss) == 1 and iss[0].severity == vi.ERROR, iss
    assert "1 sample" in iss[0].message and "(GHOST)" in iss[0].message, iss
    assert "listed 3 vs genotyped 2" in iss[0].message, iss


def test_no_more_lists_are_left_out_than_datasets_with_no_vcf_here():
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2"])
    lists = [["S1", "S2"], ["B4"], ["C5"]]
    flat = [s for one in lists for s in one]
    iss = vi.check_samplesid_overlisting([v], flat, lists, 1)
    assert len(iss) == 1 and "1 sample" in iss[0].message and "(C5)" in iss[0].message, iss
    assert vi.check_samplesid_overlisting([v], flat, lists, 2) == []


def test_a_sample_of_a_left_out_list_named_by_another_list_still_fails():
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2"])
    lists = [["B4"], ["S1", "S2", "B4"]]
    iss = vi.check_samplesid_overlisting([v], [s for one in lists for s in one], lists, 1)
    assert len(iss) == 1 and "(B4)" in iss[0].message, iss


def _datadir():
    """A data folder with genome, guide and PAM as a search reads them, two sample
    lists, and hg38_walk's VCF: returns (folder, --vcf dataset folders, the config)."""
    d = tempfile.mkdtemp()
    for sub in ("Genomes/hg38", "VCFs/hg38_walk", "samplesIDs", "Dictionaries"):
        os.makedirs(os.path.join(d, sub))
    with open(os.path.join(d, "Genomes/hg38/chr1.fa"), "w") as f:
        f.write(">chr1\nACGT\n")
    with open(os.path.join(d, "guide.txt"), "w") as f:
        f.write("CTAACAGTTGCTTTTATCACNNN\n")
    with open(os.path.join(d, "pam.txt"), "w") as f:
        f.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
    _mkvcf(os.path.join(d, "VCFs/hg38_walk"), "walk.chr1.vcf.gz", ["S1", "S2"])
    head = "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
    with open(os.path.join(d, "samplesIDs/hg38_walk.samplesID.txt"), "w") as f:
        f.write(head + "S1\tP\tS\tmale\nS2\tP\tS\tfemale\n")
    with open(os.path.join(d, "samplesIDs/hg38_walkB.samplesID.txt"), "w") as f:
        f.write(head + "B4\tP\tS\tmale\nB5\tP\tS\tfemale\n")
    config = os.path.join(d, "sid.txt")
    with open(config, "w") as f:
        f.write("hg38_walk.samplesID.txt\nhg38_walkB.samplesID.txt\n")
    return d, [os.path.join(d, "VCFs", n) for n in ("hg38_walk", "hg38_walkB")], config


def _lightweight(d, vcf_dirs, config):
    return vi.run_lightweight(
        os.path.join(d, "Genomes/hg38"), vcf_dirs, os.path.join(d, "guide.txt"),
        os.path.join(d, "pam.txt"), samplesfile=config, current_working_directory=d)


def _phantom_errors(report):
    return [line for line in report.render().splitlines()
            if "[ERROR]" in line and "genotyped in no VCF" in line]


def test_search_of_one_dataset_on_disk_beside_one_downloaded_as_an_index():
    # the documented steps: download --what vcf for one dataset, indexes of several in
    # the lists. The check the search runs refused it for the other's samples.
    d, vcf_dirs, config = _datadir()
    for sub in ("registry_hg38_walkB", "log_indels_hg38_walkB"):
        os.makedirs(os.path.join(d, "Dictionaries", sub))
    assert _phantom_errors(_lightweight(d, vcf_dirs, config)) == []


def test_search_of_one_dataset_on_disk_beside_a_sites_only_one():
    d, vcf_dirs, config = _datadir()
    os.makedirs(vcf_dirs[1])
    p = _mkvcf(vcf_dirs[1], "mega.chr1.vcf.gz", [])
    with gzip.open(p, "wt") as f:  # (8 columns: no FORMAT, no sample)
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    assert _phantom_errors(_lightweight(d, vcf_dirs, config)) == []


def test_search_with_every_dataset_on_disk_keeps_the_guard():
    d, vcf_dirs, config = _datadir()
    assert len(_phantom_errors(_lightweight(d, vcf_dirs[:1], config))) == 1  # a list of no dataset
    os.makedirs(vcf_dirs[1])
    _mkvcf(vcf_dirs[1], "walkB.chr1.vcf.gz", ["B4"])  # lists B4, B5; genotypes B4
    errors = _phantom_errors(_lightweight(d, vcf_dirs, config))
    assert len(errors) == 1 and "(B5)" in errors[0], errors


def test_a_partly_overlapping_list_beside_a_dataset_not_on_disk_still_fails():
    d, vcf_dirs, config = _datadir()
    for sub in ("registry_hg38_walkB", "log_indels_hg38_walkB"):
        os.makedirs(os.path.join(d, "Dictionaries", sub))
    with open(os.path.join(d, "samplesIDs/hg38_walk.samplesID.txt"), "a") as f:
        f.write("GHOST\tP\tS\tmale\n")
    errors = _phantom_errors(_lightweight(d, vcf_dirs, config))
    assert len(errors) == 1 and "(GHOST)" in errors[0], errors


def test_empty_inputs_are_silent():
    assert vi.check_samplesid_overlisting([], ["S1"]) == []
    assert vi.check_samplesid_overlisting(["x.vcf.gz"], []) == []


def test_panel_autofix_drops_phantoms():
    # #46 auto-fix: tier0_compile.build_sample_meta filters the panel to the
    # VCF-genotyped set (so AN = ploidy x genotyped), a strict no-op when None.
    import tier0_compile as t0c

    d = tempfile.mkdtemp()
    sid = os.path.join(d, "panel.samplesID.txt")
    with open(sid, "w") as f:
        f.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for s in ("S1", "S2", "PHANTOM"):
            f.write(f"{s}\tPOP\tSUPER\tM\n")
    db = {"DS": sid}
    all_meta, _ = t0c.build_sample_meta(db)
    noop_meta, _ = t0c.build_sample_meta(db, genotyped_samples=None)
    filt_meta, _ = t0c.build_sample_meta(db, genotyped_samples={"S1", "S2"})
    assert set(all_meta) == {"S1", "S2", "PHANTOM"}
    assert set(noop_meta) == set(all_meta)          # None -> byte-identical no-op
    assert set(filt_meta) == {"S1", "S2"}           # phantom dropped -> AN correct


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")
