"""Tests for the configurable concurrent-search-jobs setting (app_settings.py).

Default stays a conservative 2; env var overrides a persisted value; values are
clamped to [1, MAX_JOBS_CAP] so the pool is never oversubscribed."""

import os
import tempfile
import unittest

import app_settings as a


class TestMaxConcurrentJobs(unittest.TestCase):
    def setUp(self):
        os.environ.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)

    def tearDown(self):
        os.environ.pop("CRISPRME_MAX_CONCURRENT_JOBS", None)

    def test_default_is_two(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(a.get_max_concurrent_jobs(d), 2)

    def test_set_and_read_back(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(a.set_max_concurrent_jobs(4, d), 4)
            self.assertEqual(a.get_max_concurrent_jobs(d), 4)

    def test_clamped_to_cap_and_floor(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(a.set_max_concurrent_jobs(999, d), a.MAX_JOBS_CAP)
            self.assertEqual(a.set_max_concurrent_jobs(0, d), 1)

    def test_env_overrides_file(self):
        with tempfile.TemporaryDirectory() as d:
            a.set_max_concurrent_jobs(3, d)
            os.environ["CRISPRME_MAX_CONCURRENT_JOBS"] = "8"
            self.assertEqual(a.get_max_concurrent_jobs(d), 8)

    def test_corrupt_config_falls_back_to_default(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, a.CONFIG_FILE), "w") as f:
                f.write("{not valid json")
            self.assertEqual(a.get_max_concurrent_jobs(d), 2)


if __name__ == "__main__":
    unittest.main()
