#!/usr/bin/env python3
"""kit_members.py - the member list of an enriched-genome kit, build-inputs/<dataset>.tar.gz: the
md5 of each file the kit holds, by its path in the kit, one line each in md5sum's format
("<md5>  <path>"), sorted in byte order (as LC_ALL=C sorts).

A kit is the panel's enriched genome and fake-indel genome, Genomes/<genome>+<dataset> and
Genomes/<genome>+<dataset>_INDELS of the working directory it was built in, packed by the rebuild
(seq_script/rebuild/index_set_2.10.0.tsv: tar, then pigz). Its tarball is not byte-reproducible:
tar records each file's time, owner and the order the file system lists them, so a kit packed
again from the same files is another tarball. Its SHA256SUMS (build-inputs/<dataset>.SHA256SUMS)
proves that a download is the published file; this list proves what the kit holds, the same for
any packing of the same files, and checks a kit once unpacked (md5sum -c, in the folder it was
unpacked in). It goes beside the kit, build-inputs/<dataset>.members.md5 of
lucapinello/crisprme-indexes-format2, and its own SHA-256 in docs/INDEX_PROVENANCE.md
(docs/RELEASE_CHECKLIST.md: the kits; seq_script/rebuild/README.md: The publish).

    python scripts/kit_members.py <working directory> <dataset>     # from the folders the kit packs:
                                                                     # <working directory>/<dataset>.members.md5
    python scripts/kit_members.py --kit <dataset>.tar.gz --out <list> # the same list, read from the tarball
    md5sum -c <dataset>.members.md5                                  # in the folder a kit was unpacked in

Each file is read once, in blocks; the list is written whole or not at all. A member that is not
a regular file or a folder (a link, a device) is said on stderr and the list is not written: a kit
holds none. Exit 0, or 1 with what is wrong; 2 on a usage error. Standard library only.
"""

import argparse
import hashlib
import os
import sys
import tarfile

BLOCK = 1 << 20


def _md5(stream):
    h = hashlib.md5()
    for block in iter(lambda: stream.read(BLOCK), b""):
        h.update(block)
    return h.hexdigest()


def _sorted(entries):
    """``entries`` [(path, md5)] in byte order of their paths, as LC_ALL=C sort orders lines."""
    return sorted(entries, key=lambda e: os.fsencode(e[0]))


def from_folders(workdir, dataset, genome="hg38"):
    """[(path in the kit, md5)] of the files the kit of ``dataset`` packs from ``workdir``:
    Genomes/<genome>+<dataset> and Genomes/<genome>+<dataset>_INDELS, as its row's tar names them
    (paths relative to the working directory). A ValueError for a folder that is not there or a
    member that is no regular file."""
    entries, odd = [], []
    for top in (f"Genomes/{genome}+{dataset}", f"Genomes/{genome}+{dataset}_INDELS"):
        root = os.path.join(workdir, top)
        if not os.path.isdir(root):
            raise ValueError(f"{root} is not there: the kit packs it")
        for folder, dirs, files in os.walk(root):
            dirs.sort()
            for name in dirs:
                if os.path.islink(os.path.join(folder, name)):
                    odd.append(os.path.relpath(os.path.join(folder, name), workdir))
            for name in files:
                path = os.path.join(folder, name)
                rel = os.path.relpath(path, workdir).replace(os.sep, "/")
                if os.path.islink(path) or not os.path.isfile(path):
                    odd.append(rel)
                    continue
                with open(path, "rb") as fh:
                    entries.append((rel, _md5(fh)))
    if odd:
        raise ValueError("%d member(s) that are no regular file (a kit holds none): %s"
                         % (len(odd), ", ".join(sorted(odd)[:5])))
    return _sorted(entries)


def from_kit(path):
    """[(path in the kit, md5)] of the regular files of the kit tarball ``path``, read as a stream
    (gzip, as pigz writes it, or plain tar): the same list as from_folders gives of the folders it
    was packed from. A ValueError for a member that is no regular file or folder."""
    entries, odd = [], []
    with tarfile.open(path, "r|*") as tf:
        for member in tf:
            name = member.name[2:] if member.name.startswith("./") else member.name
            if member.isdir():
                continue
            if not member.isfile():
                odd.append(name)
                continue
            entries.append((name, _md5(tf.extractfile(member))))
    if odd:
        raise ValueError("%d member(s) that are no regular file (a kit holds none): %s"
                         % (len(odd), ", ".join(sorted(odd)[:5])))
    return _sorted(entries)


def write_list(entries, out):
    """Writes ``entries`` to ``out`` in md5sum's format, whole (a temporary file, then a rename)."""
    tmp = out + ".part"
    with open(tmp, "w") as fh:
        for rel, digest in entries:
            fh.write(f"{digest}  {rel}\n")
    os.replace(tmp, out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("workdir", nargs="?", help="the working directory the kit was packed in")
    ap.add_argument("dataset", nargs="?", help="the kit's dataset, e.g. hg38_HGDP1kGP")
    ap.add_argument("--genome", default="hg38", help="the genome of its folders (hg38)")
    ap.add_argument("--kit", help="read the list from this kit tarball instead")
    ap.add_argument("--out", help="the list to write (<working directory>/<dataset>.members.md5)")
    args = ap.parse_args(argv)
    if bool(args.kit) == bool(args.workdir and args.dataset) or (args.kit and args.workdir):
        ap.error("give a working directory and a dataset, or --kit <tarball> with --out")
    try:
        if args.kit:
            if not args.out:
                ap.error("--kit needs --out")
            entries, out = from_kit(args.kit), args.out
        else:
            entries = from_folders(args.workdir, args.dataset, args.genome)
            out = args.out or os.path.join(args.workdir, f"{args.dataset}.members.md5")
    except (OSError, ValueError, tarfile.TarError) as exc:
        sys.stderr.write(f"kit_members: {exc}\n")
        return 1
    if not entries:
        sys.stderr.write("kit_members: the kit holds no file\n")
        return 1
    write_list(entries, out)
    print(f"{len(entries)} file(s) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
