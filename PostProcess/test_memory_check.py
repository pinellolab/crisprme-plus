#!/usr/bin/env python3
"""Tests for memory_check.memory_status (the web pre-Submit low-memory guard, B-3).

memory_status must be non-fatal, honor the env overrides, and only report low=True
when memory is reliably detected AND below the threshold — undetectable hosts stay
silent so the web form never scares users on machines we can't measure.
"""
import os
import unittest

import memory_check as mc


class TestMemoryStatus(unittest.TestCase):
    def setUp(self):
        # snapshot + clear the env overrides so each test controls them
        self._saved = {k: os.environ.get(k) for k in ("CRISPRME_SKIP_MEM_CHECK", "CRISPRME_MIN_GB")}
        for k in self._saved:
            os.environ.pop(k, None)

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def test_returns_well_formed_dict(self):
        s = mc.memory_status()
        self.assertEqual(set(s), {"total_gb", "required_gb", "low", "message"})
        self.assertIsInstance(s["low"], bool)
        self.assertEqual(s["required_gb"], float(mc.VARIANT_WGS_MIN_GB))

    def test_skip_env_forces_not_low(self):
        os.environ["CRISPRME_SKIP_MEM_CHECK"] = "1"
        s = mc.memory_status()
        self.assertFalse(s["low"])
        self.assertEqual(s["message"], "")

    def test_huge_threshold_triggers_low_with_message(self):
        # 10 PB: any real host is below this (if memory is detectable at all)
        os.environ["CRISPRME_MIN_GB"] = str(10 ** 7)
        s = mc.memory_status()
        if s["total_gb"] is None:
            self.skipTest("total memory undetectable on this host")
        self.assertTrue(s["low"])
        self.assertIn("variant-aware", s["message"])
        self.assertIn("reference-only", s["message"])

    def test_zero_threshold_never_low(self):
        os.environ["CRISPRME_MIN_GB"] = "0"
        s = mc.memory_status()
        self.assertFalse(s["low"])

    def test_detection_margin_tolerates_near_threshold(self):
        # A host provisioned to the recommended amount reports a few GB less inside
        # a container; the _DETECT_MARGIN_GB tolerance must not nag it. Set the
        # threshold just above the detected total (within the margin) -> not low;
        # well above (beyond the margin) -> low.
        base = mc.memory_status()
        if base["total_gb"] is None:
            self.skipTest("total memory undetectable on this host")
        total = base["total_gb"]
        os.environ["CRISPRME_MIN_GB"] = str(total + (mc._DETECT_MARGIN_GB - 1))
        self.assertFalse(mc.memory_status()["low"],
                         "a threshold within the detection margin must NOT warn")
        os.environ["CRISPRME_MIN_GB"] = str(total + mc._DETECT_MARGIN_GB + 6)
        self.assertTrue(mc.memory_status()["low"],
                        "a threshold well beyond the margin must warn")

    def test_undetectable_memory_stays_silent(self):
        orig = mc._total_memory_bytes
        mc._total_memory_bytes = lambda: None
        try:
            os.environ["CRISPRME_MIN_GB"] = str(10 ** 7)
            s = mc.memory_status()
            self.assertFalse(s["low"])
            self.assertIsNone(s["total_gb"])
            self.assertEqual(s["message"], "")
        finally:
            mc._total_memory_bytes = orig

    def test_never_raises_on_bad_detection(self):
        orig = mc._total_memory_bytes
        mc._total_memory_bytes = lambda: (_ for _ in ()).throw(RuntimeError("boom"))
        try:
            s = mc.memory_status()  # must not raise
            self.assertFalse(s["low"])
        finally:
            mc._total_memory_bytes = orig

    def test_required_gb_override_parsed(self):
        os.environ["CRISPRME_MIN_GB"] = "128"
        s = mc.memory_status()
        self.assertEqual(s["required_gb"], 128.0)


if __name__ == "__main__":
    unittest.main()
