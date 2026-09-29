"""Fast, hermetic unit tests for ``complete_test.run_crisprme_test`` dispatch.

REGRESSION GUARD (the exact situation this was written for): the v2.6.0 complete-test
rewrite made ``run_crisprme_test`` *always* require a prebuilt index and run a single
example-guide smoke, which silently dropped the ``validate-benchmarks`` CI's download +
per-benchmark on-demand-build flow -> ``validate-test`` had no predictions to compare
against the brute-force ground truth and the gate went red (unnoticed because the heavy
benchmark job only runs on PR/weekly-cron, not push-to-main).

These tests pin the two-mode contract WITHOUT any download / index build / subprocess, so
they run in the fast unit suite (unit-tests.yml, which DOES run on push-to-main):

  * CI / brute-force gate  (env CRISPRME_ALLOW_ONDEMAND_BUILD set): download the test data
    and run EVERY registered benchmark via complete-search with ON-DEMAND build (no
    --index-path), one per-benchmark output dir (crisprme-test-out_<name>).
  * default USER path (flag unset): require a prebuilt index, run a single smoke against it
    via --index-path, and NEVER download/build.
"""

import os
import unittest
from unittest import mock

import complete_test as ct


_REGISTRY = {
    "thresholds": {"mm": 4, "bDNA": 1, "bRNA": 1},
    "benchmarks": [
        {"name": "cas9_a", "nuclease": "SpCas9", "pam_name": "pA.txt",
         "pam_content": "NNNNNNNNNNNNNNNNNNNNNGG 3", "guide_file": "gA.txt",
         "guide_crisprme": "A" * 20 + "NGG"},
        {"name": "cas12a_b", "nuclease": "Cas12a", "pam_name": "pB.txt",
         "pam_content": "TTTVNNNNNNNNNNNNNNNNNNNNNNNN 3", "guide_file": "gB.txt",
         "guide_crisprme": "TTTV" + "C" * 23},
    ],
}


class DispatchContract(unittest.TestCase):
    def _run(self, allow_ondemand):
        """Run run_crisprme_test with every side-effecting dependency mocked; return the
        list of complete-search command strings + a flag of whether require_prebuilt_index
        was consulted."""
        calls = []
        env = dict(os.environ)
        env.pop("CRISPRME_ALLOW_ONDEMAND_BUILD", None)
        env.pop("CRISPRME_SKIP_HEAVY", None)
        if allow_ondemand:
            env["CRISPRME_ALLOW_ONDEMAND_BUILD"] = "1"
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(ct, "check_crisprme_directory_tree", lambda *_a, **_k: None), \
             mock.patch.object(ct, "check_output", lambda *_a, **_k: None), \
             mock.patch.object(ct, "subprocess") as msub, \
             mock.patch.object(ct, "require_prebuilt_index") as mreq, \
             mock.patch.object(ct, "download_genome_data", return_value="Genomes/hg38_chr22") as mdl, \
             mock.patch.object(ct, "download_vcf_data", lambda *_a, **_k: None), \
             mock.patch.object(ct, "write_vcf_config", return_value="vcf.config.txt"), \
             mock.patch.object(ct, "download_samples_ids_data", lambda *_a, **_k: None), \
             mock.patch.object(ct, "write_samplesids_config", return_value="samples.config.txt"), \
             mock.patch.object(ct, "download_annotation_data", return_value=("gencode.bed", "encode.bed")), \
             mock.patch.object(ct, "load_benchmarks", return_value=_REGISTRY), \
             mock.patch.object(ct, "write_pamfile", side_effect=lambda name, *_a, **_k: name), \
             mock.patch.object(ct, "write_guidefile", side_effect=lambda name, *_a, **_k: name), \
             mock.patch.object(ct, "write_sg1617_guidefile", return_value="sg1617.guide.txt"), \
             mock.patch.object(ct, "_parse_index_name", return_value=("NRG", "hg38_1000G2021", True)), \
             mock.patch.object(ct, "glob", return_value=["x.integrated_results.tsv"]), \
             mock.patch("os.path.isdir", return_value=True), \
             mock.patch("os.path.isfile", return_value=True):
            mreq.return_value = "NRG_3_hg38+hg38_1000G2021"
            msub.call.side_effect = lambda cmd, **_k: calls.append(cmd) or 0
            ct.run_crisprme_test("chr22", "1000G", 4, False)
        return calls, mreq.called, mdl.called

    def test_ci_mode_downloads_and_builds_per_benchmark_on_demand(self):
        """CRISPRME_ALLOW_ONDEMAND_BUILD set -> download + one on-demand-build complete-search
        per registered benchmark; NEVER require_prebuilt_index. This is the exact flow the
        brute-force benchmark gate needs; its loss is what this guard catches."""
        calls, req_called, dl_called = self._run(allow_ondemand=True)
        self.assertTrue(dl_called, "CI mode must download the test genome")
        self.assertFalse(req_called, "CI mode must NOT gate on a prebuilt index")
        self.assertEqual(len(calls), len(_REGISTRY["benchmarks"]),
                         "CI mode must run one complete-search per registered benchmark")
        for cmd in calls:
            self.assertIn("complete-search", cmd)
            self.assertNotIn("--index-path", cmd,
                             "CI mode must BUILD on demand (no --index-path), not reuse a prebuilt index")
        # each benchmark gets its own crisprme-test-out_<name> dir (validate-test reads these)
        for bench in _REGISTRY["benchmarks"]:
            self.assertTrue(any(f"{ct.COMPLETETESTRESDIR}_{bench['name']}" in c for c in calls),
                            f"missing per-benchmark output dir for {bench['name']}")

    def test_user_mode_requires_prebuilt_index_and_never_builds(self):
        """Flag unset -> require a prebuilt index + a single smoke via --index-path; NEVER
        download/build (Luca's rule: CRISPRme does not build an index for the user)."""
        calls, req_called, dl_called = self._run(allow_ondemand=False)
        self.assertTrue(req_called, "user mode must require a prebuilt index")
        self.assertFalse(dl_called, "user mode must NOT download/build test data")
        self.assertEqual(len(calls), 1, "user mode runs exactly one example smoke")
        self.assertIn("--index-path genome_library", calls[0],
                      "user-mode smoke must reuse the prebuilt index (--index-path)")


if __name__ == "__main__":
    unittest.main()
