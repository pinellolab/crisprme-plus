#!/usr/bin/env python3
"""The installers and the docs give the measured install and image sizes, and the same ones.

The one-line installers, the route the README and the website call the easiest, said the
engine image is ~800 MB: `docker pull pinellolab/crisprme:latest` downloads 4.84 GB on
amd64 (what Docker Desktop on Windows, Intel Macs and most Linux hosts pull) and 1.81 GB on
arm64, summed from the registry manifests of every layer, and it unpacks to ~7 GB or more.
They also asked for a data folder with ~45 GB free, which is the download, not the space it
takes. The README, the website, DATA_LAYOUT and the Quickstart said ~85 GB on disk and ~100 GB
free; TROUBLESHOOTING gave ~3 + 40 + 50 GB and ~150 GB free. Measured with `du -sb` on a clean
`download --what all` + `download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`:
97,924,848,937 bytes (genome 3.27 GB, indexes 41.49 GB, Dictionaries/ 52.44 GB, annotations
0.72 GB), and the install peaks at ~103 GB while the genotype store unpacks beside its archive.
So ~100 GB free was full before the first result was written.

Read as text. A re-measure changes the constants below and every surface together.
"""

import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

INSTALLED_GB = "98"   # the default variant-aware install, on disk
FREE_GB = "150"       # what every surface asks a user to keep free for it

# Every surface that sizes the default variant-aware install.
SURFACES = ("README.md", "docs/index.html", "docs/DOCKER_QUICKSTART.md", "docs/DATA_LAYOUT.md",
            "docs/TROUBLESHOOTING.md", "install/README.md", "docs/install.sh", "docs/install.ps1")
INSTALLERS = ("docs/install.sh", "docs/install.ps1")


def text(name):
    """The file as one line of prose: HTML entities and Markdown emphasis and quote markers
    dropped, whitespace runs folded."""
    s = (REPO / name).read_text(encoding="utf-8")
    for entity, char in (("&nbsp;", " "), ("&asymp;", "≈"), ("&rarr;", "→"), ("&mdash;", "—"),
                         ("&ndash;", "–")):
        s = s.replace(entity, char)
    s = re.sub(r"^\s*> ?", "", s, flags=re.M).replace("**", "")
    return " ".join(s.split())


def readme_image_download():
    """The amd64 and arm64 download sizes README §1.2 gives for the image."""
    m = re.search(r"The pull downloads about ([\d.]+) GB on amd64 and ([\d.]+) GB on Apple "
                  r"Silicon", text("README.md"))
    return m.group(1), m.group(2)


def image_line(name):
    """The line an installer prints before it pulls the engine image."""
    line = next(ln for ln in (REPO / name).read_text(encoding="utf-8").splitlines()
                if "Fetching the CRISPRme engine image" in ln)
    return line


class TheInstallersSayWhatThePullDownloads(unittest.TestCase):

    def test_the_readme_sentence_they_follow_is_found(self):
        self.assertEqual(readme_image_download(), ("5", "1.8"))

    def test_install_sh_gives_both_downloads(self):
        amd64, arm64 = readme_image_download()
        line = image_line("docs/install.sh")
        self.assertIn(f"~{amd64} GB download", line)
        self.assertIn(f"~{arm64} GB on Apple Silicon", line)
        self.assertNotIn("MB", line)

    def test_install_ps1_gives_the_amd64_download(self):
        # Docker Desktop on Windows pulls linux/amd64
        amd64, _ = readme_image_download()
        line = image_line("docs/install.ps1")
        self.assertIn(f"~{amd64} GB download", line)
        self.assertNotIn("MB", line)

    def test_no_installer_says_800_mb(self):
        for name in INSTALLERS:
            with self.subTest(installer=name):
                self.assertNotIn("800 MB", text(name))


class EverySurfaceGivesTheMeasuredInstallSize(unittest.TestCase):

    def test_each_names_the_installed_size_and_the_free_space(self):
        for name in SURFACES:
            with self.subTest(surface=name):
                body = text(name)
                self.assertIn(f"~{INSTALLED_GB} GB", body)
                self.assertRegex(body, rf"\b{FREE_GB} GB free")

    def test_every_free_space_figure_is_the_same(self):
        for name in SURFACES:
            with self.subTest(surface=name):
                said = set(re.findall(r"(\d+) GB free", text(name)))
                self.assertEqual(said, {FREE_GB})

    def test_every_installed_size_figure_is_the_same(self):
        for name in SURFACES:
            with self.subTest(surface=name):
                said = set(re.findall(r"~(\d+) GB (?:on disk|once installed)", text(name)))
                self.assertTrue(said <= {INSTALLED_GB}, said)

    def test_the_stale_figures_are_gone(self):
        for name in SURFACES:
            body = text(name)
            for stale in ("~85 GB", "≈ 100 GB", "~100 GB free", "~45 GB free",
                          "store ~45 GB of data"):
                with self.subTest(surface=name, stale=stale):
                    self.assertNotIn(stale, body)

    def test_the_breakdowns_list_every_dictionaries_folder(self):
        # the Quickstart left out the indel genotype store and the indel logs (14 GB)
        for name in ("docs/DOCKER_QUICKSTART.md", "docs/TROUBLESHOOTING.md"):
            body = text(name)
            for part in ("genotype store", "allele-frequency registry", "indel genotype store",
                         "indel logs"):
                with self.subTest(surface=name, part=part):
                    self.assertIn(part, body)


if __name__ == "__main__":
    unittest.main(verbosity=2)
