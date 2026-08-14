#!/usr/bin/env python3
"""Email-notification config tests: write/read round-trip (password not echoed,
kept on re-save), validation, 0600 perms, and send_mail resolving .email.json."""
import sys, tempfile, os, json
sys.path.insert(0, "/opt/conda/opt/crisprme")
from pages import pages_utils as pu  # noqa: E402

d = tempfile.mkdtemp(prefix="email_")
pu.current_working_directory = d + "/"
fails = []
def chk(c, m):
    print(("PASS " if c else "FAIL "), m)
    if not c: fails.append(m)

c = pu.read_email_config()
chk(c["smtp_host"] == "smtp.gmail.com" and c["password_set"] is False, "defaults when unset")
chk(pu.write_email_config("smtp.gmail.com", 465, "me@gmail.com", "app-pw-123", True) is None, "save valid")
c = pu.read_email_config()
chk(c["sender"] == "me@gmail.com" and c["password_set"] is True and "password" not in c, "password not echoed, flag set")
chk(pu.write_email_config("smtp.example.com", 587, "me@gmail.com", None, False) is None, "re-save w/ password=None")
raw = json.load(open(os.path.join(d, ".email.json")))
chk(raw["password"] == "app-pw-123" and raw["smtp_host"] == "smtp.example.com" and raw["use_ssl"] is False, "kept password, updated fields")
chk(pu.write_email_config("", 465, "x@y.com", "p") == "SMTP host is required.", "reject empty host")
chk(pu.write_email_config("h", "notaport", "x@y.com", "p") == "SMTP port must be a number.", "reject bad port")
chk(pu.write_email_config("h", 465, "bad", "p") == "Sender must be a valid email address.", "reject bad sender")
chk(oct(os.stat(os.path.join(d, ".email.json")).st_mode & 0o777) == "0o600", "0600 perms")

job = os.path.join(d, "Results", "J"); os.makedirs(job)
open(os.path.join(job, "email.txt"), "w").write("u@v.com\nhttp://link\n2026\n")
from pages import send_mail
cfg = send_mail._load_config(job)
chk(cfg["host"] == "smtp.example.com" and cfg["port"] == 587 and cfg["password"] == "app-pw-123" and cfg["use_ssl"] is False, "send_mail resolves .email.json")

print(f"\n{'FAILED' if fails else 'OK'}: {len(fails)} failure(s)")
os._exit(1 if fails else 0)
