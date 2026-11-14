#!/usr/bin/env python3
"""Playwright end-to-end check of the Settings / Data Manager uploads (Gate B walk
P4, P5, P11, P15) against a running web app, driven like a user in Chromium.

1. A genome uploaded as one FASTA of two sequences is offered on the cards (the
   upload splits it, which PostProcess/test_settings_progress.py checks), and a
   reference index of it builds (4).
2. The chunked VCF uploader sends the reference genome picked on the card (it read
   React-Select v2 classes that Dash 2's dcc.Dropdown does not have, and sent ''),
   so the dataset is offered on the search form for that genome even though its
   name does not hold the genome's, and labelled genotyped, not sites-only.
3. A BED whose file name names no installed genome can be enabled under 'Manage
   annotations', and a second upload whose bytes equal the first is still added.
4. 'Started' shows in the neutral colour, and goes once the job has ended; a
   reference build ends with a plain 'done.'.
5. A genome uploaded again under that name (its sequences swapped) is refused while
   the index built on it is installed, as searches of it would go on using that
   index with the old positions; once the index is removed under 'Remove installed
   data', it goes in (final review round 1, F2).

Standalone (no pytest); exits non-zero on failure. The web app must already be
serving at CRISPRME_WEB_URL (default http://localhost:8080), on a data folder this
may add a genome, a PAM, a VCF dataset and an annotation to:

    CRISPRME_WEB_URL=http://localhost:8080 python test/web/test_settings_upload_e2e.py
"""

import gzip
import os
import random
import re
import sys
import time

from playwright.sync_api import sync_playwright

BASE = os.environ.get("CRISPRME_WEB_URL", "http://localhost:8080").rstrip("/")
GENOME = "e2egenome"  # no other name on the page holds it
DATASET = "e2epanel"  # a name that does not hold the genome's, as the card's example

_fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL "), msg)
    sys.stdout.flush()
    if not cond:
        _fails.append(msg)


def goto(page, path):
    page.goto(BASE + path, wait_until="networkidle", timeout=60000)
    page.wait_for_timeout(1500)


def pick(page, dd_id, label_regex):
    """Select the dcc.Dropdown option whose label matches, as a user clicks it."""
    page.locator(f"#{dd_id} .Select-control").scroll_into_view_if_needed()
    page.locator(f"#{dd_id} .Select-control").click()
    page.wait_for_timeout(400)
    opts = page.locator(".Select-menu-outer .VirtualizedSelectOption, .Select-menu-outer .Select-option")
    labels = [opts.nth(i).inner_text() for i in range(opts.count())]
    for i, label in enumerate(labels):
        if re.search(label_regex, label):
            opts.nth(i).click()
            page.wait_for_timeout(500)
            return label
    page.keyboard.press("Escape")
    raise AssertionError(f"no option {label_regex!r} in #{dd_id}: {labels}")


def options(page, dd_id):
    page.locator(f"#{dd_id} .Select-control").click()
    page.wait_for_timeout(400)
    opts = page.locator(".Select-menu-outer .VirtualizedSelectOption, .Select-menu-outer .Select-option")
    labels = [opts.nth(i).inner_text() for i in range(opts.count())]
    page.keyboard.press("Escape")
    page.wait_for_timeout(200)
    return labels


def chunk_upload(page, target, name, data, timeout_s=60):
    widget = page.locator(f".crisprme-chunk-upload[data-target={target}]")
    widget.scroll_into_view_if_needed()
    widget.locator("input[type=file]").set_input_files(
        {"name": name, "mimeType": "application/octet-stream", "buffer": data})
    widget.locator("button").click()
    span = widget.locator("span")
    text = ""
    for _ in range(timeout_s * 4):
        text = span.inner_text()
        if text.startswith(("Upload complete", "Error", "Upload failed")):
            break
        page.wait_for_timeout(250)
    # a refusal in the cards' red, all else in grey (a refusal was grey too)
    colour = span.evaluate("el => getComputedStyle(el).color")
    refused = text.startswith(("Error", "Upload failed"))
    check(colour == ("rgb(187, 0, 0)" if refused else "rgb(85, 85, 85)"),
          f"the uploader's {'refusal' if refused else 'message'} is shown in "
          f"{'red' if refused else 'grey'} ({colour}: {text[:60]!r})")
    return text


def feedback(page, sel):
    """(text, colour) of a card's feedback, as shown."""
    return page.evaluate(
        """(sel) => { const el = document.querySelector(sel);
            const kid = el.firstElementChild || el;
            return [el.innerText, getComputedStyle(kid).color]; }""", sel)


def main():
    random.seed(281)
    seqs = {c: "".join(random.choice("ACGT") for _ in range(2000)) for c in ("chrA", "chrB")}
    fasta = "".join(f">{c}\n{s}\n" for c, s in seqs.items()).encode()
    ref = seqs["chrA"][99]
    alt = "ACGT"[("ACGT".index(ref) + 1) % 4]
    vcf = ("##fileformat=VCFv4.2\n##contig=<ID=chrA>\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\n"
           f"chrA\t100\trs1\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0|1\t1|1\n").encode()
    bed = b"chrA\t10\t50\tregion1\n"
    js_errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": 1400, "height": 1000})
        page.on("pageerror", lambda e: js_errors.append(str(e)))
        sent = []
        page.on("request", lambda r: sent.append(r.headers)
                if "/settings/upload-chunk" in r.url else None)

        # 1. a genome of two sequences in one FASTA file
        goto(page, "/settings")
        got = chunk_upload(page, "genome", f"{GENOME}.fa", fasta)
        check(got.startswith("Upload complete"), f"genome upload completes ({got!r})")
        goto(page, "/settings")
        check(any(o == GENOME for o in options(page, "vcf-ref-genome")),
              "the uploaded genome is offered on the VCF card")

        # 2. a VCF dataset whose name does not hold the genome's
        page.fill("#vcf-name", DATASET)
        pick(page, "vcf-ref-genome", rf"^{GENOME}$")
        got = chunk_upload(page, "vcf", f"{DATASET}.chrA.vcf.gz", gzip.compress(vcf))
        check(got.startswith("Upload complete"), f"VCF upload completes ({got!r})")
        genomes = [h.get("x-genome") for h in sent if h.get("x-target") == "vcf"]
        check(genomes == [GENOME], f"the uploader sends the genome picked (X-Genome {genomes})")
        goto(page, "/index")
        pick(page, "available-genome", rf"^{GENOME}$")
        page.wait_for_timeout(1500)
        variants = options(page, "variant-dataset")
        mine = [v for v in variants if v.startswith(DATASET)]
        check(bool(mine), f"the search form offers {DATASET} for {GENOME} ({variants})")
        check(bool(mine) and "genotyped, samples not indexed" in mine[0],
              f"labelled genotyped, not sites-only ({mine})")

        # 3. annotation BEDs
        goto(page, "/settings")
        page.locator("#annotation-upload input[type=file]").set_input_files(
            {"name": "e2e_regions.bed", "mimeType": "text/plain", "buffer": bed})
        page.wait_for_timeout(2500)
        first = page.locator("#annotation-feedback").inner_text()
        page.locator("#annotation-upload input[type=file]").set_input_files(
            {"name": "e2e_regions_copy.bed", "mimeType": "text/plain", "buffer": bed})
        page.wait_for_timeout(2500)
        second = page.locator("#annotation-feedback").inner_text()
        check("e2e_regions.bed" in first, f"first BED added ({first!r})")
        check("e2e_regions_copy.bed" in second,
              f"a second BED with the same bytes is added too ({second!r})")
        goto(page, "/settings")
        pick(page, "ann-manage-genome", rf"^{GENOME}$")
        page.wait_for_timeout(1200)
        boxes = page.locator("#ann-enabled-checklist label").all_inner_texts()
        check("e2e_regions.bed" in [b.strip() for b in boxes],
              f"a BED naming no genome can be enabled for {GENOME} ({boxes})")

        # 4. a job: a PAM, then a reference index build of the genome
        page.fill("#pam-name", "E2Ecas")
        page.fill("#pam-motif", "NRG")
        page.click("#pam-add-btn")
        page.wait_for_timeout(1500)
        goto(page, "/settings")
        pick(page, "index-build-genome", rf"^{GENOME}$")
        pick(page, "index-build-pam", r"NRG")
        page.click("#index-build-btn")
        page.wait_for_timeout(800)
        text, colour = feedback(page, "#index-feedback")
        check(text.startswith("Started") and colour == "rgb(85, 85, 85)",
              f"'Started' in the neutral colour ({text!r}, {colour})")
        done = ""
        deadline = time.time() + 600
        while time.time() < deadline:
            done = page.locator("#settings-progress").inner_text()
            if ": done." in done or "failed" in done.lower():
                break
            page.wait_for_timeout(1000)
        page.wait_for_timeout(2500)
        check(done.strip() == "Build index: done.", f"the build ends with a plain done ({done!r})")
        text, _colour = feedback(page, "#index-feedback")
        check(text == "", f"'Started' goes once the job has ended ({text!r})")

        # 5. the genome uploaded again, its sequences swapped, while its index is there
        swapped = f">chrA\n{seqs['chrB']}\n>chrB\n{seqs['chrA']}\n".encode()
        goto(page, "/settings")
        got = chunk_upload(page, "genome", f"{GENOME}.fa", swapped)
        index = [o[len("Index: "):] for o in options(page, "delete-item")
                 if o.startswith("Index: ") and o.endswith(f"_{GENOME}")]
        check(len(index) == 1 and got.startswith(
                  f"Error: Indexes or variant genomes built on a genome named '{GENOME}' are "
                  f"installed (genome_library/{index[0]})"),
              f"an upload is refused while an index of a genome of its name is there "
              f"({index}, {got!r})")
        page.on("dialog", lambda d: d.accept())
        pick(page, "delete-item", rf"^Index: .*_{GENOME}$")
        page.click("#delete-btn")
        page.wait_for_timeout(2000)
        gone = page.locator("#delete-feedback").inner_text()
        check(gone.startswith("Deleted index:"), f"the index is removed ({gone!r})")
        goto(page, "/settings")
        got = chunk_upload(page, "genome", f"{GENOME}.fa", swapped)
        check(got.startswith("Upload complete"), f"then the upload goes in ({got!r})")

        if js_errors:
            check(False, f"uncaught JS error(s): {js_errors}")
        browser.close()
    print(f"\n{'FAILED' if _fails else 'OK'}: {len(_fails)} failure(s)")
    sys.exit(1 if _fails else 0)


if __name__ == "__main__":
    main()
