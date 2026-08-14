"""Send email to notify the user that CRISPRme analysis was completed and
the corresponding results are available to be visualized and explored.

Structure of email.txt
    - Contact
    - link to CRISPRme job
    - Job submission date

Configuration (resolved in this order; first hit wins):
    1. <data_dir>/.email.json  -- written by Settings -> Email notifications
       (local mode). Keys: smtp_host, smtp_port, use_ssl, sender, password.
       <data_dir> is the CRISPRme working directory (two levels above the job
       results folder passed on argv[1], i.e. <data_dir>/Results/<job>).
    2. Environment: CRISPRME_SMTP_PSW (password/app password) and
       CRISPRME_SMTP_SENDER (sender address); host/port default to gmail SSL.

Email notifications are OPTIONAL and best-effort: with no password configured,
or if delivery fails, send_mail() no-ops gracefully (message to STDOUT only) so a
completed search is never lost or marked failed. IMPORTANT: this script must never
write to stderr or raise -- the pipeline treats any stderr output as a fatal job
error ([ -s log_error.txt ]).

TODO: avoid shell and call send_mail() in other python scripts
TODO: add run parameters to mail (job date + other params)
"""

import os
import sys
import ssl
import json
import smtplib
from email.message import EmailMessage

# default mail-server settings (used when neither the config file nor the
# environment override them). gmail over implicit SSL.
_DEFAULT_HOST = "smtp.gmail.com"
_DEFAULT_PORT = 465
_DEFAULT_SENDER = "crisprme.job@gmail.com"


def _load_config(output_folder: str) -> dict:
    """Resolve SMTP settings from <data_dir>/.email.json, then the environment.

    Returns a dict with keys host, port, use_ssl, sender, password. ``password``
    is None when nothing is configured (caller then no-ops).
    """
    cfg = {
        "host": _DEFAULT_HOST,
        "port": _DEFAULT_PORT,
        "use_ssl": True,
        "sender": os.environ.get("CRISPRME_SMTP_SENDER", _DEFAULT_SENDER),
        "password": os.environ.get("CRISPRME_SMTP_PSW"),
    }
    # <output_folder> = <data_dir>/Results/<job>  ->  <data_dir>
    data_dir = os.path.dirname(os.path.dirname(os.path.abspath(output_folder.rstrip("/"))))
    config_path = os.path.join(data_dir, ".email.json")
    try:
        if os.path.isfile(config_path):
            with open(config_path) as fh:
                filecfg = json.load(fh)
            if isinstance(filecfg, dict):
                if filecfg.get("smtp_host"):
                    cfg["host"] = str(filecfg["smtp_host"])
                if filecfg.get("smtp_port"):
                    cfg["port"] = int(filecfg["smtp_port"])
                if "use_ssl" in filecfg:
                    cfg["use_ssl"] = bool(filecfg["use_ssl"])
                if filecfg.get("sender"):
                    cfg["sender"] = str(filecfg["sender"])
                if filecfg.get("password"):
                    cfg["password"] = str(filecfg["password"])
    except Exception as exc:  # never let a bad config file fail the job
        print(f"Could not read email config ({config_path}): {exc}")
    return cfg


def _deliver(cfg: dict, msg: EmailMessage) -> None:
    """Open an SMTP connection per cfg and send one message."""
    context = ssl.create_default_context()
    if cfg["use_ssl"]:
        with smtplib.SMTP_SSL(cfg["host"], cfg["port"], context=context) as server:
            server.login(cfg["sender"], cfg["password"])
            server.send_message(msg)
    else:
        with smtplib.SMTP(cfg["host"], cfg["port"]) as server:
            server.starttls(context=context)
            server.login(cfg["sender"], cfg["password"])
            server.send_message(msg)


def send_mail() -> None:
    """Send the completion-notification email(s) for a finished CRISPRme job.

    Best-effort: prints to STDOUT and returns on any problem; never writes to
    stderr and never raises (a mail failure must not fail a completed search).
    """
    output_folder = sys.argv[1]
    cfg = _load_config(output_folder)
    if not cfg["password"]:
        print(
            "Email notifications not configured (Settings -> Email notifications, "
            "or set CRISPRME_SMTP_PSW) - skipping notification."
        )
        return
    try:
        with open(os.path.join(output_folder, "email.txt"), "r") as e:
            all_content = e.read().strip().split("--OTHEREMAIL--")
    except OSError as exc:
        print(f"No email recipients to notify ({exc}).")
        return

    for em in all_content:
        em = em.strip().split("\n")
        if not em or not em[0]:
            continue
        try:
            msg = EmailMessage()
            msg["To"] = em[0]
            job_link = em[1] if len(em) > 1 else ""
            msg["Subject"] = "CRISPRme - Job completed"
            msg["From"] = cfg["sender"]
            msg.set_content(
                "The requested job is completed, visit the following link "
                + job_link
                + " to view the report."
            )
            _deliver(cfg, msg)
            print(f"EMAIL SENT to {em[0]}")
        except Exception as exc:  # best-effort: log to stdout, keep going
            print(f"Email notification to {em[0]} failed: {exc}")


def main():
    """Call send_mail() to notify users of CRISPRme job completion."""
    try:
        send_mail()
    except Exception as exc:  # last-resort guard: never emit stderr / fail the job
        print(f"Email notification skipped ({exc}).")


# entry point
if __name__ == "__main__":
    main()
