"""The web app follows no link another user put in its data folder (the review of lane S,
round 5, S-R5-2; the live check of its refix, read_links_check.py, case by case).

A member of the data folder's group can write Results/, Settings/ and the job folders in them
under a umask of 002, and the web app listens on every network interface. It followed a link
put there: the Settings FAILED banner and Show full log, the search's log panel and the
results downloads (GET /Results/<path>) showed and served any file of the user's to anyone
who could reach it, and a Settings job wrote its log through one. The inventory of what
the web app reads, serves and writes found more of the kind, here too: the status page's
banner, the History table, the Settings index labels; the chunked and the dcc.Upload
uploads, the email settings (with the password), the tables' filter, the mode file and a
search's logs; and the tables' cache, pickles kept in Cache/ of the data folder, which
Flask-Caching loads.

Each case runs through the real Dash server (index.app's Flask test client: the requests a
browser sends, the callbacks' JSON) or the real function that writes, in a data folder whose
folders the group can write, with a secret (0600) and a victim file of the user's beside it.
A refused case answers without the secret, and leaves the victim as it was; legitimate files
are still read and served.

    python -m unittest test_web_no_links -v      (from PostProcess/; needs Dash)
"""

import base64
import json
import os
import secrets
import shutil
import stat
import sys
import tempfile
import unittest
import warnings
from unittest import mock

try:  # the pages need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
except ImportError:  # pragma: no cover
    dash = None
# Without Dash each class skips, not the module (a SkipTest raised at import is a traceback
# under `python -m unittest <module>`).
_NEEDS_DASH = unittest.skipUnless(dash, "dash not installed: web page tests skipped")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

if dash:
    sys.argv = ["crisprme.py", "--version"]  # (importing the pages starts the app object)
    import app as app_module  # noqa: E402
    import index  # noqa: E402
    import pages.history_page as hp  # noqa: E402
    import pages.load_page as lp  # noqa: E402
    import pages.main_page as mp  # noqa: E402
    import pages.pages_utils as pu  # noqa: E402
    import pages.results_page as rp  # noqa: E402
    import pages.settings_page as sp  # noqa: E402

_MODULES = ("app_module", "index", "hp", "lp", "mp", "pu", "rp", "sp")


class _DataFolder(unittest.TestCase):
    """The web app's data folder (index.check_directories), each folder writable by the
    group as under a umask of 002, beside the user's home: secret.txt (0600, a token a
    reader would show) and victim.txt (what a write through a link would change)."""

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)  # (files the page code leaves open)
        self.top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.top, ignore_errors=True)
        self.data = os.path.join(self.top, "data")
        self.home = os.path.join(self.top, "home")
        os.makedirs(self.home)
        os.makedirs(self.data)
        index.check_directories(self.data)
        for root, dirs, _files in os.walk(self.data):
            os.chmod(root, 0o775)
        self.token = "SECRET-" + secrets.token_hex(12)
        self.secret = os.path.join(self.home, "secret.txt")
        # each reader's leak, as the code before showed it: a line, an ERROR line of the
        # status page's, a .Params.txt field of the History table's
        with open(self.secret, "w") as fh:
            fh.write(f"{self.token}\n"
                     f"ERROR: {self.token} is unfinished: run build-index-only again to its "
                     f"end\nGenome_selected\t{self.token}\nGenome_idx\tNRG_2_hg38\n"
                     "Mismatches\t4\nDNA\t1\nRNA\t1\nPam\tNRG\n")
        os.chmod(self.secret, 0o600)
        self.victim = os.path.join(self.home, "victim.txt")
        with open(self.victim, "w") as fh:
            fh.write("VICTIM\n")
        for name in _MODULES:
            p = mock.patch.object(globals()[name], "current_working_directory",
                                  self.data + os.sep)
            p.start()
            self.addCleanup(p.stop)
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.data)  # (where the web app runs)
        self.client = index.app.server.test_client()
        self.deps = self.client.get("/_dash-dependencies").get_json()

    # --- the requests a browser sends -----------------------------------------------
    def _dep(self, input_id, input_prop, output=None):
        hits = [d for d in self.deps if any(i["id"] == input_id and i["property"] == input_prop
                                            for i in d["inputs"])
                and (output is None or output in d["output"])]
        self.assertEqual(len(hits), 1, f"callbacks on {input_id}.{input_prop}")
        return hits[0]

    @staticmethod
    def _outputs(dep):
        out = dep["output"]
        parts = out[2:-2].split("...") if out.startswith("..") else [out]
        res = [{"id": p.rsplit(".", 1)[0], "property": p.rsplit(".", 1)[1]} for p in parts]
        return res if out.startswith("..") else res[0]

    def call(self, input_id, input_prop, inputs, state=None, output=None):
        """POST /_dash-update-component for the callback on ``input_id.input_prop`` (the one
        whose output names ``output``, when several), with ``inputs`` and ``state``
        ({"id.prop": value}); (HTTP status, response text)."""
        dep = self._dep(input_id, input_prop, output)
        state = state or {}
        body = {
            "output": dep["output"],
            "outputs": self._outputs(dep),
            "inputs": [{"id": i["id"], "property": i["property"],
                        "value": inputs.get(f"{i['id']}.{i['property']}")}
                       for i in dep["inputs"]],
            "state": [{"id": s["id"], "property": s["property"],
                       "value": state.get(f"{s['id']}.{s['property']}")}
                      for s in dep["state"]],
            "changedPropIds": [f"{input_id}.{input_prop}"],
        }
        r = self.client.post("/_dash-update-component", json=body)
        return r.status_code, r.get_data(as_text=True)

    def page(self, pathname, search=""):
        """The page the address shows (index.change_page), as JSON text."""
        url = f"http://127.0.0.1:8080{pathname}{search}"
        return self.call("url", "pathname", {"url.href": url, "url.pathname": pathname,
                                              "url.search": search}, {"url.hash": ""})

    # --- what another user puts there -------------------------------------------------
    def folder(self, *parts):
        d = os.path.join(self.data, *parts)
        os.makedirs(d, exist_ok=True)
        os.chmod(d, 0o775)
        return d

    def write(self, path, text):
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def link(self, path, target=None):
        os.symlink(target or self.secret, path)
        return path

    def assertNoSecret(self, text):
        self.assertNotIn(self.token, text)

    def assertVictimUnchanged(self):
        with open(self.victim) as fh:
            self.assertEqual(fh.read(), "VICTIM\n")


# ===================================================================================
# The six cases of the live check (read_links_check.py)
# ===================================================================================
@_NEEDS_DASH
class TestTheLiveCheck(_DataFolder):
    FAILED = "Download genome\tQueued\nDownload genome\tStart\nDownload genome\tFAILED\n"

    def settings_job(self, job_id, logs_linked=True):
        job = self.folder("Settings", job_id)
        self.write(os.path.join(job, "log.txt"), self.FAILED)
        for name in ("log_error.txt", "log_verbose.txt"):
            if logs_linked:
                self.link(os.path.join(job, name))
            else:
                self.write(os.path.join(job, name), f"the job's own {name}\n")
        return job

    def test_1a_the_settings_failed_banner(self):
        self.settings_job("set_PLANTED000")
        code, text = self.call("settings-check", "n_intervals", {"settings-check.n_intervals": 1},
                               {"settings-active-job.data": "set_PLANTED000"})
        self.assertEqual(code, 200)
        self.assertIn("Download genome: failed.", text)
        self.assertNoSecret(text)
        self.assertIn("log_error.txt not shown", text)
        # a job's own log is shown, as before
        self.settings_job("set_OWNLOGS000", logs_linked=False)
        code, text = self.call("settings-check", "n_intervals", {"settings-check.n_intervals": 1},
                               {"settings-active-job.data": "set_OWNLOGS000"})
        self.assertIn("the job's own log_error.txt", text)

    def test_1b_the_settings_show_full_log(self):
        self.settings_job("set_PLANTED000")
        code, text = self.call("settings-show-log-btn", "n_clicks",
                               {"settings-show-log-btn.n_clicks": 1},
                               {"settings-full-log-container.hidden": True,
                                "settings-failed-job.data": "set_PLANTED000"})
        self.assertEqual(code, 200)
        self.assertNoSecret(text)
        self.assertIn("not shown", text)
        self.settings_job("set_OWNLOGS000", logs_linked=False)
        code, text = self.call("settings-show-log-btn", "n_clicks",
                               {"settings-show-log-btn.n_clicks": 1},
                               {"settings-full-log-container.hidden": True,
                                "settings-failed-job.data": "set_OWNLOGS000"})
        self.assertIn("the job's own log_verbose.txt", text)

    def test_2_the_search_log_panel(self):
        job = self.folder("Results", "PLANTEDJOB01")
        self.link(os.path.join(job, "log_error.txt"))
        self.write(os.path.join(job, "log_verbose.txt"), "the search's own step\n")
        code, text = self.call("button-show-log", "n_clicks", {"button-show-log.n_clicks": 1},
                               {"url.search": "?job=PLANTEDJOB01",
                                "full-log-container.hidden": True})
        self.assertEqual(code, 200)
        self.assertNoSecret(text)
        self.assertIn("the search's own step", text)  # (the job's own log, as before)

    def test_3a_a_download_of_a_link_in_a_job_folder(self):
        job = self.folder("Results", "PLANTEDJOB01")
        self.link(os.path.join(job, "planted.txt"))
        self.write(os.path.join(job, "results.tsv"), "the job's own results\n")
        r = self.client.get("/Results/PLANTEDJOB01/planted.txt")
        self.assertNoSecret(r.get_data(as_text=True))
        self.assertEqual(r.status_code, 403)
        r = self.client.get("/Results/PLANTEDJOB01/results.tsv")  # a file of the job's
        self.assertEqual((r.status_code, r.get_data(as_text=True)),
                         (200, "the job's own results\n"))
        self.assertIn("attachment", r.headers["Content-Disposition"])

    def test_3b_a_download_through_a_folder_that_is_a_link(self):
        self.link(os.path.join(self.data, "Results", "evil"), self.home)
        for path in ("/Results/evil/secret.txt", "/Results/evil/../evil/secret.txt"):
            with self.subTest(path=path):
                r = self.client.get(path)
                self.assertNoSecret(r.get_data(as_text=True))
                self.assertIn(r.status_code, (403, 404))
        # (an empty step: not the download route's, the app's page)
        r = self.client.get("/Results//evil/secret.txt")
        self.assertNoSecret(r.get_data(as_text=True))

    def test_4_a_settings_job_writes_no_log_through_a_link(self):
        job = self.folder("Settings", "set_WRITES0000")
        self.link(os.path.join(job, "log_verbose.txt"), self.victim)
        sp._run_settings_job(["echo", "WRITTEN-BY-THE-JOB"], job, "Download genome")
        self.assertVictimUnchanged()
        with open(os.path.join(job, "log.txt")) as fh:
            self.assertIn("Download genome\tFAILED", fh.read())  # (refused: not run)


# ===================================================================================
# What the inventory added
# ===================================================================================
@_NEEDS_DASH
class TestReadersTheInventoryAdded(_DataFolder):
    def test_the_status_page_banner(self):
        # the banner shows the ERROR lines of log_error.txt the pipeline writes
        job = self.folder("Results", "PLANTEDJOB02")
        self.write(os.path.join(job, ".Params.txt"), "Ref_comp\tFalse\n")
        self.write(os.path.join(job, "log.txt"), "Job\tStart\tFri\n")
        self.link(os.path.join(job, "log_error.txt"))
        code, text = self.call("load-page-check", "n_intervals",
                               {"load-page-check.n_intervals": 1},
                               {"url.search": "?job=PLANTEDJOB02"})
        self.assertEqual(code, 200)
        self.assertNoSecret(text)

    def test_the_history_table(self):
        # (a file of the user's whose lines read as the parameters the table shows)
        params = os.path.join(self.home, "params.txt")
        with open(params, "w") as fh:
            fh.write(f"Genome_selected\t{self.token}\nGenome_idx\tNRG_2_hg38\nMismatches\t4\n"
                     "DNA\t1\nRNA\t1\nPam\tNRG\n")
        os.chmod(params, 0o600)
        job = self.folder("Results", "PLANTEDJOB03")
        self.link(os.path.join(job, ".Params.txt"), params)
        own = self.folder("Results", "OWNJOB0004")
        self.write(os.path.join(own, ".Params.txt"),
                   "Genome_selected\tOWNGENOME\nGenome_idx\tNRG_2_hg38\nMismatches\t4\n"
                   "DNA\t1\nRNA\t1\nPam\tNRG\n")
        for d in (job, own):  # (the files the table reads besides: a job's guides and log)
            self.write(os.path.join(d, ".guides.txt"), "CTAACAGTTGCTTTTATCACNNN\n")
            self.write(os.path.join(d, "log.txt"), "Job\tStart\tTue Nov 3 09:49:01 UTC 2026\n")
        code, text = self.page("/history")
        self.assertEqual(code, 200)
        self.assertNoSecret(text)
        self.assertIn("OWNGENOME", text)  # (a job's own parameters, as before)

    def test_the_settings_index_labels(self):
        idx = self.folder("genome_library", "NRG_2_hg38")
        self.link(os.path.join(idx, ".display_label"))
        # (the page's catalogs of HuggingFace and HPRC, which it fetches, not fetched here)
        with mock.patch.object(sp, "_hf_catalog", lambda component: []), \
                mock.patch.object(sp, "_HPRC_SAMPLE_IDS", []):
            code, text = self.page("/settings")
        self.assertEqual(code, 200)
        self.assertNoSecret(text)
        self.assertIn("NRG_2_hg38", text)

    def test_a_job_folder_that_is_a_link_opens_no_page(self):
        # wherever it leads: out of Results/ (refused before too), or to another job folder
        # (followed before: no job folder is a link)
        outside = os.path.join(self.top, "outside_job")
        os.makedirs(outside)
        self.write(os.path.join(outside, "log.txt"), self.token + "\n")
        self.link(os.path.join(self.data, "Results", "LINKEDJOB05"), outside)
        real = self.folder("Results", "REALJOB0005")
        self.write(os.path.join(real, "log.txt"), "Job\tStart\tFri\n")
        self.link(os.path.join(self.data, "Results", "ALIASJOB005"), real)
        for job in ("LINKEDJOB05", "ALIASJOB005"):
            with self.subTest(job=job):
                self.assertIsNone(pu.results_job_dir(job))
                code, text = self.page("/load", f"?job={job}")
                self.assertNoSecret(text)
                self.assertIn("No job to display", text)
        self.assertEqual(pu.results_job_dir("REALJOB0005"), real)


@_NEEDS_DASH
class TestWritersTheInventoryAdded(_DataFolder):
    def test_a_chunked_upload_writes_no_part_through_a_link(self):
        uploads = self.folder("Settings", "uploads")
        self.link(os.path.join(uploads, "walk.bed.part"), self.victim)
        r = self.client.post("/settings/upload-chunk", data=b"chr1\t1\t2\tPLANTED\n",
                             headers={"X-Target": "annotation", "X-File-Name": "walk.bed",
                                      "X-Chunk-Index": "0", "X-Total-Chunks": "2"})
        self.assertVictimUnchanged()
        self.assertEqual(r.status_code, 400)
        # an upload of its own goes in, as before
        for i, chunk in enumerate((b"chr1\t1\t2\tone\n", b"chr1\t5\t9\ttwo\n")):
            r = self.client.post("/settings/upload-chunk", data=chunk,
                                 headers={"X-Target": "annotation", "X-File-Name": "own.bed",
                                          "X-Chunk-Index": str(i), "X-Total-Chunks": "2"})
        self.assertEqual(r.get_data(as_text=True), "complete:own.bed")
        with open(os.path.join(self.data, "Annotations", "own.bed")) as fh:
            self.assertEqual(fh.read(), "chr1\t1\t2\tone\nchr1\t5\t9\ttwo\n")

    def test_an_annotation_upload_writes_nothing_through_a_link(self):
        anns = os.path.join(self.data, "Annotations")
        self.link(os.path.join(anns, "walk.bed.uploading"), self.victim)
        self.link(os.path.join(anns, "walk.bed"), self.victim)
        contents = "data:text/plain;base64," + base64.b64encode(b"chr1\t1\t2\tx\n").decode()
        code, _text = self.call("annotation-upload", "contents",
                                {"annotation-upload.contents": contents},
                                {"annotation-upload.filename": "walk.bed"})
        self.assertEqual(code, 200)
        self.assertVictimUnchanged()
        with open(os.path.join(anns, "walk.bed")) as fh:  # (the link replaced)
            self.assertEqual(fh.read(), "chr1\t1\t2\tx\n")

    def test_the_email_settings_put_no_password_through_a_link(self):
        self.link(os.path.join(self.data, ".email.json.tmp"), self.victim)
        self.link(os.path.join(self.data, ".email.json"), self.victim)
        code, text = self.call("email-save-btn", "n_clicks", {"email-save-btn.n_clicks": 1},
                               {"email-smtp-host.value": "smtp.gmail.com",
                                "email-smtp-port.value": 465,
                                "email-sender.value": "me@example.org",
                                "email-password.value": "APP-PASSWORD",
                                "email-use-ssl.value": True})
        self.assertEqual(code, 200)
        self.assertIn("Saved.", text)
        self.assertVictimUnchanged()
        cfg = os.path.join(self.data, ".email.json")
        self.assertFalse(os.path.islink(cfg))
        self.assertEqual(stat.S_IMODE(os.stat(cfg).st_mode), 0o600)
        with open(cfg) as fh:
            self.assertEqual(json.load(fh)["password"], "APP-PASSWORD")

    def test_the_tables_filter_writes_nothing_through_a_link(self):
        job = self.folder("Results", "FILTERJOB01")
        self.link(os.path.join(job, ".dropdown.json"), self.victim)
        code, _text = self.call("target_filter_dropdown", "value",
                                {"target_filter_dropdown.value": "CFD"},
                                {"url.search": "?job=FILTERJOB01"}, output="store.data")
        self.assertEqual(code, 200)
        self.assertVictimUnchanged()
        self.assertEqual(pu.read_json("FILTERJOB01"), "CFD")

    def test_the_annotation_manifest_and_the_active_annotation(self):
        # the manifest went through .enabled.json.tmp, the active annotation through
        # .active.<genome>.bed.building and its sort, and its signature by name
        anns = os.path.join(self.data, "Annotations")
        for name, rows in (("a.hg38.bed", "chr1\t10\t500\tpromoter\n"),
                           ("b.hg38.bed", "chr1\t900\t1500\tenhancer\n")):
            self.write(os.path.join(anns, name), rows)
        for name in (".enabled.json.tmp", ".active.hg38.bed.building",
                     ".active.hg38.bed.tmp.sorted.bed", ".active.hg38.sig"):
            self.link(os.path.join(anns, name), self.victim)
        pu.write_enabled_annotations("hg38", ["a.hg38.bed", "b.hg38.bed"])
        self.assertEqual(sorted(pu.read_enabled_annotations("hg38")),
                         ["a.hg38.bed", "b.hg38.bed"])
        if shutil.which("sort-bed") and shutil.which("bgzip"):  # (the image's BEDOPS, htslib)
            active, _gencode = pu.build_active_annotation("hg38")
            self.assertEqual(active, ".active.hg38.bed.gz")
            self.assertTrue(os.path.isfile(os.path.join(anns, active)))
        self.assertVictimUnchanged()

    def test_the_mode_file_is_written_through_no_link(self):
        self.link(os.path.join(self.data, ".mode_type.txt"), self.victim)
        with mock.patch.object(index, "run_server"):  # (the server itself is not started)
            index.index()
        self.assertVictimUnchanged()
        with open(os.path.join(self.data, ".mode_type.txt")) as fh:
            self.assertEqual(fh.read(), "local")

    def test_a_search_writes_no_log_through_a_link(self):
        job = self.folder("Results", "SEARCHJOB01")
        self.link(os.path.join(job, "log_verbose.txt"), self.victim)
        mp._run_search_job(["echo", "WRITTEN-BY-THE-SEARCH"], {},
                           os.path.join(job, "log_verbose.txt"),
                           os.path.join(job, "log_error.txt"), self.data)
        self.assertVictimUnchanged()
        own = self.folder("Results", "SEARCHJOB02")  # a search's own logs, as before
        mp._run_search_job(["sh", "-c", "echo out; echo err >&2"], {"X": "1"},
                           os.path.join(own, "log_verbose.txt"),
                           os.path.join(own, "log_error.txt"), self.data)
        with open(os.path.join(own, "log_verbose.txt")) as fh:
            self.assertEqual(fh.read(), "out\n")
        with open(os.path.join(own, "log_error.txt")) as fh:
            self.assertEqual(fh.read(), "err\n")

    def test_the_tables_cache_is_a_folder_of_the_server_s_own(self):
        # Flask-Caching loads the pickles of its folder: Cache/ of the data folder, which the
        # group can write, let another user's pickle run as the web app's user
        with index.app.server.app_context():
            backend = app_module.cache.cache
        path = getattr(backend, "_path", None)
        if path is None:  # (a temporary folder others can write: in memory)
            self.assertEqual(type(backend).__name__, "SimpleCache")
        else:
            path = os.path.abspath(path)
            self.assertFalse(path.startswith(os.path.abspath(self.data) + os.sep), path)
            self.assertNotEqual(os.path.basename(path), "Cache")
            st = os.stat(path)
            self.assertEqual(stat.S_IMODE(st.st_mode), 0o700)
            self.assertEqual(st.st_uid, os.geteuid())


if __name__ == "__main__":
    unittest.main()
