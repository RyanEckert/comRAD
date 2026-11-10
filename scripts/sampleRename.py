#!/usr/bin/env python3
################################################################################
#Script Name    : sampleRename.py
#Description    : Rename files from a two-column csv of old name and sample number
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Rename files in the working directory from a two-column csv.

Column 1 is the existing name without its extension, column 2 the sample
number. The new name is <-n prefix><column 2>.<-f extension>.

Every rename is checked before any rename happens, because the failures are
silent otherwise:

  * two rows giving the same new name would have had the second overwrite the
    first
  * a new name that already belongs to a different file would be overwritten
  * a file present with the right extension but absent from the csv keeps its
    old name, and is then missed by whatever globs the new names afterwards

All three are reported and, for the first two, refuse to run. --dry-run shows
the whole plan and changes nothing.

A new name may be the old name of another file in the same batch (renumbering
SGM004 -> SGM003 while SGM003 -> SGM002, or swapping two names). That file is
being moved out of the way, so it is not a clash. Renames are done in two
passes, every old name to a temporary name first and then every temporary name
to its new name, so chains and swaps never overwrite anything.
"""
import argparse
import csv
import os
import sys

prog = os.path.basename(sys.argv[0])
tmp_suffix = ".sampleRename.tmp"


def die(msg):
    sys.exit(f"{prog}: {msg}")


def read_pairs(path, ext, prefix, header):
    """[(old, new)] from the csv, with malformed rows reported not skipped."""
    pairs = []
    with open(path, newline="") as fh:
        for lineno, row in enumerate(csv.reader(fh), start=1):
            if not row or not row[0].strip():
                continue
            if header and lineno == 1:
                continue
            if len(row) < 2 or not row[1].strip():
                die(f"{path} line {lineno} has no sample number: {row}")
            pairs.append((f"{row[0].strip()}.{ext}",
                          f"{prefix}{row[1].strip()}.{ext}"))
    if not pairs:
        die(f"{path} has no usable rows")
    return pairs


def main():
    p = argparse.ArgumentParser(
        prog=prog,
        description="Rename files based on columns of data from a specified "
                    ".csv file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # rename trimmed, demultiplexed reads to sample names
  {prog} -i sampleNames -f tr0.gz -n MGM_

  # see the whole plan first, change nothing
  {prog} -i sampleNames -f tr0.gz -n MGM_ --dry-run

  # csv has a header row
  {prog} -i sampleNames -f tr0.gz -n MGM_ --header

csv:
  Two columns, no quoting needed: the existing name without its extension, and
  the sample number.

    FG6-6_ACCA,001
    FG6-6_TGGT,002

  With -f tr0.gz -n MGM_ that renames FG6-6_ACCA.tr0.gz to MGM_001.tr0.gz.
  -f takes the whole extension without the leading dot, so two-part
  extensions like 'tr0.gz' or 'trim.gz' work.

  A new name may be another row's old name (renumbering, or swapping two
  samples); those are renamed through temporary names so nothing is
  overwritten.
""")
    p.add_argument("-f", type=str, default="fastq", metavar="EXT",
                   help="file extension, not including the leading '.'; "
                        "two-part extensions such as tr0.gz are fine "
                        "[default: %(default)s]")
    p.add_argument("-i", type=str, default="sampleList", metavar="CSV",
                   help="input .csv, with or without the .csv suffix, in the "
                        "working directory [default: %(default)s]")
    p.add_argument("-n", type=str, default="", metavar="PREFIX",
                   help="common name given to all samples, preceding the "
                        "sample number from the .csv")
    p.add_argument("--header", action="store_true",
                   help="skip the first row of the .csv")
    p.add_argument("--dry-run", action="store_true",
                   help="report the plan and rename nothing")
    a = p.parse_args()

    path = a.i if a.i.endswith(".csv") else a.i + ".csv"
    if not os.path.isfile(path):
        die(f"cannot open {path}")

    pairs = read_pairs(path, a.f, a.n, a.header)

    # every check runs before any rename
    seen = {}
    clash = []
    for old, new in pairs:
        if new in seen:
            clash.append((new, seen[new], old))
        seen[new] = old
    if clash:
        msg = "\n".join(f"  {new} <- both {a1} and {a2}"
                        for new, a1, a2 in clash)
        die(f"{len(clash)} new name(s) claimed by more than one row:\n{msg}")

    present = {f for f in os.listdir(".")
               if f.endswith("." + a.f) and os.path.isfile(f)}
    todo = [(old, new) for old, new in pairs if old in present]
    missing = [old for old, _ in pairs if old not in present]

    # a target is only taken if it exists and is not itself moving away
    moving = {old for old, _ in todo}
    blocked = [(old, new) for old, new in todo
               if os.path.exists(new) and new not in moving]
    if blocked:
        msg = "\n".join(f"  {new} already exists, wanted by {old}"
                        for old, new in blocked)
        die(f"{len(blocked)} target name(s) already taken:\n{msg}")

    todo = [(old, new) for old, new in todo if old != new]
    stale = [old + tmp_suffix for old, _ in todo
             if os.path.exists(old + tmp_suffix)]
    if stale:
        msg = "\n".join(f"  {t}" for t in stale[:10])
        die(f"{len(stale)} temporary name(s) already exist, left by an "
            f"interrupted run? Resolve these first:\n{msg}")

    renamed = {old for old, _ in pairs if old in present}
    extra = sorted(present - renamed)

    for old, new in todo:
        print(f"{old}\t{new}")

    if not a.dry_run:
        # pass 1: old -> temporary; on failure put them back, nothing changed
        parked = []
        try:
            for old, _ in todo:
                os.rename(old, old + tmp_suffix)
                parked.append(old)
        except OSError as e:
            for old in reversed(parked):
                os.rename(old + tmp_suffix, old)
            die(f"rename failed, all names restored: {e}")
        # pass 2: temporary -> new
        for i, (old, new) in enumerate(todo):
            try:
                os.rename(old + tmp_suffix, new)
            except OSError as e:
                left = "\n".join(f"  {o + tmp_suffix} -> {n}"
                                 for o, n in todo[i:])
                die(f"rename failed: {e}\n{len(todo) - i} file(s) still hold "
                    f"temporary names and should become:\n{left}")

    print(f"\n{prog}: {len(todo)} renamed"
          f"{' (dry run, nothing changed)' if a.dry_run else ''}",
          file=sys.stderr)
    if missing:
        print(f"{prog}: {len(missing)} row(s) in {path} had no matching file:",
              file=sys.stderr)
        for old in missing[:10]:
            print(f"  {old}", file=sys.stderr)
    if extra:
        print(f"{prog}: WARNING {len(extra)} *.{a.f} file(s) are not in "
              f"{path} and keep their old names -- anything globbing the new "
              f"names will miss them:", file=sys.stderr)
        for f in extra[:10]:
            print(f"  {f}", file=sys.stderr)


if __name__ == "__main__":
    main()
