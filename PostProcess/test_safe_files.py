"""safe_files: how the web app reads, writes and serves the data folder's files without
following a link someone else put there (the review of lane S, round 5, S-R5-2).

A member of the data folder's group can write its folders under a umask of 002: Results/,
Settings/ and the job folders in them, the email settings beside them. The web app, which
listens on every network interface, followed a link put there: it showed, and served for
download, any file of the user's, and wrote through it. Each case here is a link (or a FIFO,
or a hard link) put where the web app reads or writes, and what safe_files does with it;
the pages' own cases, through the real Dash server, are test_web_no_links.

Standard library only (pandas, sqlite3 and Flask where a case needs them), as safe_files
is: the hermetic unit-tests job runs it.

    python -m unittest test_safe_files -v      (from PostProcess/)
"""

import errno
import os
import shutil
import sqlite3
import stat
import sys
import tempfile
import threading
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import safe_files  # noqa: E402

SECRET = "SECRET-OF-THE-USER\n"


class _Data(unittest.TestCase):
    """A data folder (data/, group-writable, as under a umask of 002) with Results/<job>/,
    Settings/ and Genomes/, beside a folder of the user's elsewhere (home/) holding a secret
    and a file that must not be written."""

    def setUp(self):
        self.top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.top, ignore_errors=True)
        self.home = os.path.join(self.top, "home")
        self.data = os.path.join(self.top, "data")
        self.job = os.path.join(self.data, "Results", "JOB")
        for d in (self.home, self.job, os.path.join(self.data, "Settings"),
                  os.path.join(self.data, "Genomes")):
            os.makedirs(d)
        for d in (self.data, os.path.join(self.data, "Results"), self.job,
                  os.path.join(self.data, "Settings"), os.path.join(self.data, "Genomes")):
            os.chmod(d, 0o775)
        self.secret = os.path.join(self.home, "secret.txt")
        self.victim = os.path.join(self.home, "victim.txt")
        with open(self.secret, "w") as fh:
            fh.write(SECRET)
        with open(self.victim, "w") as fh:
            fh.write("VICTIM\n")
        with open(os.path.join(self.job, "log.txt"), "w") as fh:
            fh.write("Job\tDone\n")

    def victim_text(self):
        with open(self.victim) as fh:
            return fh.read()

    def plant(self, name, target):
        path = os.path.join(self.job, name)
        os.symlink(target, path)
        return path


class TestReadsFollowNoLink(_Data):
    def test_a_regular_file_is_read(self):
        self.assertEqual(safe_files.read_text(os.path.join(self.job, "log.txt"), self.data),
                         "Job\tDone\n")
        self.assertTrue(safe_files.is_file(os.path.join(self.job, "log.txt"), self.data))

    def test_a_link_in_a_job_folder_is_refused(self):
        link = self.plant("log_error.txt", self.secret)
        for read in (lambda: safe_files.read_text(link, self.data),
                     lambda: safe_files.read_bytes(link, self.data),
                     lambda: safe_files.open_fd(link, self.data),
                     lambda: safe_files.stat_file(link, self.data)):
            with self.assertRaises(safe_files.RefusedPath) as said:
                read()
            self.assertIn("log_error.txt is a link", str(said.exception))
        self.assertFalse(safe_files.is_file(link, self.data))

    def test_a_link_to_a_folder_at_any_step_is_refused(self):
        os.symlink(self.home, os.path.join(self.data, "Results", "evil"))
        evil = os.path.join(self.data, "Results", "evil")
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.read_text(os.path.join(evil, "secret.txt"), self.data)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.listdir(evil, self.data)
        self.assertFalse(safe_files.is_dir(evil, self.data))
        self.assertEqual(safe_files.glob(evil, "*", self.data), [])

    def test_a_job_folder_that_is_a_link_inside_results_is_refused(self):
        # (a link that stays inside Results/ is refused too: the job folders hold none)
        os.symlink(self.job, os.path.join(self.data, "Results", "ALIAS"))
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.read_text(os.path.join(self.data, "Results", "ALIAS", "log.txt"),
                                 self.data)

    def test_a_fifo_a_device_and_a_folder_are_not_regular_files(self):
        fifo = os.path.join(self.job, "fifo")
        os.mkfifo(fifo)
        for path in (fifo, os.path.join(self.data, "Results")):
            with self.assertRaises(safe_files.RefusedPath) as said:
                safe_files.read_text(path, self.data)  # (a FIFO would hold the request)
            self.assertIn("not a regular file", str(said.exception))

    def test_a_hard_link_is_refused_where_others_can_write(self):
        hard = os.path.join(self.job, "hard.txt")
        os.link(self.secret, hard)
        with self.assertRaises(safe_files.RefusedPath) as said:
            safe_files.read_text(hard, self.data)
        self.assertIn("hard link", str(said.exception))
        os.chmod(self.job, 0o755)  # where only the user writes, the user made it
        self.assertEqual(safe_files.read_text(hard, self.data), SECRET)

    def test_a_path_out_of_the_data_folder_or_with_dot_dot_is_refused(self):
        for path in (self.secret, self.data,
                     os.path.join(self.job, "..", "..", "..", "home", "secret.txt"),
                     # a step '..' after a name that is not a folder: never taken out
                     os.path.join(self.job, "nosuch..", "..", "log.txt")):
            with self.subTest(path=path):
                with self.assertRaises(safe_files.RefusedPath):
                    safe_files.read_text(path, self.data)

    def test_the_data_folder_s_own_folder_may_be_the_user_s_link_where_only_they_write(self):
        # Results/ linked onto a volume: followed when no one else can have put the link
        volume = os.path.join(self.top, "volume")
        shutil.move(os.path.join(self.data, "Results"), volume)
        os.symlink(volume, os.path.join(self.data, "Results"))
        log = os.path.join(self.data, "Results", "JOB", "log.txt")
        os.chmod(self.data, 0o755)
        self.assertEqual(safe_files.read_text(log, self.data), "Job\tDone\n")
        os.chmod(self.data, 0o775)  # the group can put its own there: refused
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.read_text(log, self.data)
        os.chmod(self.data, 0o1775)  # sticky: no one renames or removes the user's link
        self.assertEqual(safe_files.read_text(log, self.data), "Job\tDone\n")

    def test_own_links_are_followed_in_the_link_layout_folders_only_where_safe(self):
        shared = os.path.join(self.top, "shared_hg38")
        os.makedirs(shared)
        with open(os.path.join(shared, "chr1.fa"), "w") as fh:
            fh.write(">chr1\nACGT\n")
        os.symlink(shared, os.path.join(self.data, "Genomes", "hg38"))
        fasta = os.path.join(self.data, "Genomes", "hg38", "chr1.fa")
        os.chmod(os.path.join(self.data, "Genomes"), 0o755)
        self.assertEqual(safe_files.read_text(fasta, self.data), ">chr1\nACGT\n")
        with self.assertRaises(safe_files.RefusedPath):  # strict, asked so
            safe_files.read_text(fasta, self.data, links=safe_files.LINKS_NONE)
        os.chmod(os.path.join(self.data, "Genomes"), 0o775)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.read_text(fasta, self.data)

    def test_open_nofollow_refuses_a_link_at_the_last_step(self):
        link = self.plant("x_integrated_results.tsv", self.secret)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.open_nofollow(link)
        with safe_files.open_nofollow(os.path.join(self.job, "log.txt")) as fh:
            self.assertEqual(fh.read(), "Job\tDone\n")

    def test_opened_path_names_the_very_file_opened(self):
        log = os.path.join(self.job, "log.txt")
        with safe_files.opened_path(log, self.data) as name:
            os.rename(log, log + ".old")  # put aside, and a link put at its name
            os.symlink(self.secret, log)
            with open(name) as fh:
                self.assertEqual(fh.read(), "Job\tDone\n")


class TestWritesFollowNoLink(_Data):
    def test_replace_renames_over_a_link_never_writing_through_it(self):
        link = self.plant(".dropdown.json", self.victim)
        safe_files.write_text(link, self.data, "CFD")
        self.assertEqual(self.victim_text(), "VICTIM\n")
        self.assertFalse(os.path.islink(link))
        with open(link) as fh:
            self.assertEqual(fh.read(), "CFD")
        self.assertEqual(sorted(os.listdir(self.job)), [".dropdown.json", "log.txt"])

    def test_replace_takes_the_mode_asked_and_leaves_nothing_when_it_fails(self):
        path = os.path.join(self.data, ".email.json")
        safe_files.write_text(path, self.data, "{}", perm=0o600)
        self.assertEqual(stat.S_IMODE(os.stat(path).st_mode), 0o600)
        with self.assertRaises(RuntimeError):
            with safe_files.replace(path, self.data) as fh:
                fh.write("half")
                raise RuntimeError("stopped")
        with open(path) as fh:
            self.assertEqual(fh.read(), "{}")
        self.assertEqual([n for n in os.listdir(self.data) if n.endswith(".tmp")], [])

    def test_append_and_write_refuse_a_link_and_a_hard_link(self):
        link = self.plant("log_verbose.txt", self.victim)
        hard = os.path.join(self.job, "log_error.txt")
        os.link(self.victim, hard)
        for path in (link, hard):
            for opener in (safe_files.open_append, safe_files.open_write):
                with self.subTest(path=os.path.basename(path), opener=opener.__name__):
                    with self.assertRaises(safe_files.RefusedPath):
                        opener(path, self.data)
        self.assertEqual(self.victim_text(), "VICTIM\n")

    def test_append_and_write_of_a_regular_file(self):
        log = os.path.join(self.job, "log.txt")
        with safe_files.open_append(log, self.data) as fh:
            fh.write("more\n")
        self.assertEqual(safe_files.read_text(log, self.data), "Job\tDone\nmore\n")
        with safe_files.open_write(log, self.data) as fh:
            fh.write("new\n")
        self.assertEqual(safe_files.read_text(log, self.data), "new\n")
        safe_files.touch(os.path.join(self.job, "queue.txt"), self.data)
        self.assertTrue(safe_files.is_file(os.path.join(self.job, "queue.txt"), self.data))

    def test_create_refuses_anything_there(self):
        link = self.plant("guides.txt", self.victim)
        with self.assertRaises(FileExistsError):
            safe_files.create(link, self.data)
        self.assertEqual(self.victim_text(), "VICTIM\n")

    def test_a_write_through_a_folder_that_is_a_link_is_refused(self):
        os.symlink(self.home, os.path.join(self.data, "Results", "evil"))
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.write_text(os.path.join(self.data, "Results", "evil", "victim.txt"),
                                  self.data, "OVER")
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.mkdir(os.path.join(self.data, "Results", "evil", "new"), self.data)
        self.assertEqual(self.victim_text(), "VICTIM\n")
        self.assertEqual(sorted(os.listdir(self.home)), ["secret.txt", "victim.txt"])

    def test_rename_moves_a_regular_file_only(self):
        upload = os.path.join(self.data, "Settings", "x.bed.part")
        os.symlink(self.secret, upload)
        dest = os.path.join(self.job, "x.bed")
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.rename(upload, dest, self.data)
        self.assertFalse(os.path.lexists(dest))
        os.remove(upload)
        with open(upload, "w") as fh:
            fh.write("chr1\t1\t2\ta\n")
        os.symlink(self.victim, dest)  # a link where it goes: replaced
        safe_files.rename(upload, dest, self.data)
        self.assertFalse(os.path.islink(dest))
        self.assertEqual(self.victim_text(), "VICTIM\n")

    def test_rmtree_and_remove_act_in_the_folder_reached_and_never_through_a_link(self):
        os.symlink(self.home, os.path.join(self.data, "Results", "evil"))
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.rmtree(os.path.join(self.data, "Results", "evil"), self.data)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.rmtree(os.path.join(self.data, "Results", "evil", "x"), self.data)
        safe_files.remove(os.path.join(self.data, "Results", "evil"), self.data)  # the link
        self.assertEqual(sorted(os.listdir(self.home)), ["secret.txt", "victim.txt"])
        safe_files.rmtree(self.job, self.data)
        self.assertFalse(os.path.exists(self.job))

    def test_two_threads_replacing_one_file_each_leave_it_whole(self):
        path = os.path.join(self.job, "list.txt")
        texts = ["A" * 100000 + "\n", "B" * 100000 + "\n"]
        threads = [threading.Thread(target=safe_files.write_text, args=(path, self.data, t))
                   for t in texts]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertIn(safe_files.read_text(path, self.data), texts)


class TestWhatOtherCodeOpens(_Data):
    def test_sqlite_connects_read_only_to_a_regular_file_only(self):
        db = os.path.join(self.job, ".JOB.db")
        with sqlite3.connect(db) as conn:
            conn.execute("CREATE TABLE t (a)")
            conn.execute("INSERT INTO t VALUES (42)")
        conn = safe_files.sqlite_connect(db, self.data)
        self.assertEqual(conn.execute("SELECT a FROM t").fetchall(), [(42,)])
        with self.assertRaises(sqlite3.OperationalError):
            conn.execute("INSERT INTO t VALUES (1)")
        conn.close()
        link = self.plant(".LINK.db", db)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.sqlite_connect(link, self.data)
        missing = os.path.join(self.job, ".NONE.db")
        with self.assertRaises(FileNotFoundError):  # (sqlite3.connect made an empty one)
            safe_files.sqlite_connect(missing, self.data)
        self.assertFalse(os.path.exists(missing))

    def test_zip_file_zips_a_regular_file_and_never_a_link(self):
        with open(os.path.join(self.job, "t.tsv"), "w") as fh:
            fh.write("a\tb\n")
        safe_files.zip_file(os.path.join(self.job, "t.tsv"), os.path.join(self.job, "t.zip"),
                            self.data)
        with zipfile.ZipFile(os.path.join(self.job, "t.zip")) as zf:
            self.assertEqual(zf.namelist(), ["t.tsv"])
            self.assertEqual(zf.read("t.tsv"), b"a\tb\n")
        link = self.plant("x_integrated_results.tsv", self.secret)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.zip_file(link, os.path.join(self.job, "x.zip"), self.data)
        self.assertFalse(os.path.exists(os.path.join(self.job, "x.zip")))
        # a link where the zip goes is replaced, not read as an archive to add to
        os.symlink(self.victim, os.path.join(self.job, "u.zip"))
        safe_files.zip_file(os.path.join(self.job, "t.tsv"), os.path.join(self.job, "u.zip"),
                            self.data)
        self.assertEqual(self.victim_text(), "VICTIM\n")

    def test_read_csv_and_to_csv(self):
        try:
            import pandas as pd
        except ImportError:  # pragma: no cover
            self.skipTest("pandas not installed")
        link = self.plant("summary.txt", self.victim)
        safe_files.to_csv(pd.DataFrame({"a": [1, 2]}), link, self.data, index=False)
        self.assertEqual(self.victim_text(), "VICTIM\n")
        self.assertEqual(list(safe_files.read_csv(link, self.data)["a"]), [1, 2])
        os.remove(link)
        os.symlink(self.secret, link)
        with self.assertRaises(safe_files.RefusedPath):
            safe_files.read_csv(link, self.data)

    def test_send_file_serves_a_regular_file_and_refuses_a_link(self):
        try:
            import flask
        except ImportError:  # pragma: no cover
            self.skipTest("Flask not installed")
        app = flask.Flask(__name__)
        with app.test_request_context("/"):
            response = safe_files.send_file(os.path.join(self.job, "log.txt"), self.data)
            response.direct_passthrough = False
            self.assertEqual(response.get_data(), b"Job\tDone\n")
            self.assertIn("attachment", response.headers["Content-Disposition"])
            response.close()
            link = self.plant("planted.txt", self.secret)
            with self.assertRaises(safe_files.RefusedPath):
                safe_files.send_file(link, self.data)


if __name__ == "__main__":
    unittest.main()
