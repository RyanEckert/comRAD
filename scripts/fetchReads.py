#!/usr/bin/env python3
################################################################################
#Script Name    : fetchReads.py
#Description    : Download pooled Illumina reads from BaseSpace, merge lane
#                 parts, and flatten to one file per pool
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Download reads from BaseSpace, optionally merge lane parts, flatten names.

Replaces downloadReads.sh and mergeReads.sh. Nothing is copied and nothing is
decompressed.

At this stage a file is a POOL of samples, told apart later by in-line barcode
in trim2bRAD.py. Sample identity is pool plus barcode, so two pools must never
end up in one file: barcode ACCA would then hold reads from two specimens with
no way to separate them. Everything below protects that.

  1. each --project is downloaded into its own subdirectory of --outdir. Only
     the read files are fetched -- not the run reports, logs, InterOp or
     sample sheets a project also holds -- and not Undetermined reads unless
     --undetermined is given
  2. every read file under --outdir is found, however deeply nested
  3. the Illumina fields are stripped from the END of each basename, so what
     is left is the pool name even if it contains underscores
  4. a pool name found in more than one project gets the project appended, so
     FG6-6 in two projects becomes FG6-6-JA23123 and FG6-6-JA23124. These are
     never merged -- identically named files in different projects are either
     the same file twice or two different pools, and guessing wrong silently
     mixes specimens
  5. lane parts -- same pool, same project, different _L00N_ -- are merged with
     --merge-reads. gzip streams concatenate, so this is a byte copy that zcat,
     python gzip, cutadapt and bowtie2 all read transparently
  6. every read lands flat in --outdir, then anything else bs left behind --
     the per-project json manifest -- is moved up too and the download
     subdirectories are removed, so nothing stays nested

  FG6-6_S22_L003_R1_001.fastq.gz  ->  FG6-6.fastq.gz   (--merge-reads)
  FG6-6_S22_L004_R1_001.fastq.gz  /
  FG6-6_S22_R1_001.fastq.gz       ->  FG6-6.fastq.gz   (no lane splitting)
  My_Pool_2_S3_L001_R1_001.fastq.gz -> My_Pool_2.fastq.gz

The download filter is bs's own --exclude/--include (BaseSpace CLI >= 1.4.0):
UNIX wildcard patterns matched against each file's path within the project,
applied in order starting from every file included, so
'--exclude=* --include=*.fastq.gz --exclude=*Undetermined*' keeps the reads
and then drops Undetermined.

Lane splitting is a demultiplexer setting, not something in the data: bcl2fastq
split per lane by default, while BCL Convert -- standard on NovaSeq X and
NextSeq 2000 -- is commonly run with it off. A pool with no _L00N_ is a single
file and needs no merge at all.

A read designator is kept only when a pool really has more than one read, so
single-end output is '<pool>.fastq.gz'. R1 and R2 are never merged.

A problem with one pool does not stop the run. That pool's files are moved to
--outdir under their original Illumina names -- long names mark them as
unprocessed -- the reason is printed, the remaining pools are processed, and
the exit status is non-zero so a launcher job reports the failure.
"""
import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys

READS_EXT = ".fastq.gz"
# Illumina fields, stripped from the end of the basename in this order
FIELDS = (re.compile(r"_001$"), re.compile(r"_(R[12])$"),
          re.compile(r"_L(\d{3})$"), re.compile(r"_S(\d+)$"))

prog = os.path.basename(sys.argv[0])
skipped = []


def die(msg):
    sys.exit(f"{prog}: {msg}")


def warn(msg):
    print(f"{prog}: {msg}", file=sys.stderr)


def skip(pool, reason):
    skipped.append((pool, reason))
    warn(f"skipping {pool}: {reason}")


def parse_name(basename, ext):
    """(pool, read, lane) from an Illumina filename, or (None, ..) if empty.

    Fields come off the end, so anything left is the pool name, underscores
    included. Absent fields are simply absent.
    """
    stem = basename[:-len(ext)]
    read = lane = None
    for pattern in FIELDS:
        m = pattern.search(stem)
        if not m:
            continue
        stem = stem[:m.start()]
        if pattern is FIELDS[1]:
            read = m.group(1)
        elif pattern is FIELDS[2]:
            lane = m.group(1)
    return (stem or None), read, lane


def download_filter(ext, undetermined):
    """bs --exclude/--include arguments: read files only, Undetermined opt-in."""
    args = ["--exclude=*", f"--include=*{ext}"]
    if not undetermined:
        args.append("--exclude=*Undetermined*")
    return args


def download(projects, outdir, dry, ext, undetermined):
    """One subdirectory per project, so a name clash can name its project."""
    if projects and not shutil.which("bs") and not dry:
        die("bs not found on PATH; install BaseSpace CLI, or omit -p to work "
            "on files already downloaded")
    for project in projects:
        dest = os.path.join(outdir, project)
        cmd = ["bs", "download", "project", "--concurrency=high", "-q",
               "-n", project, "-o", dest] + download_filter(ext, undetermined)
        warn(" ".join(shlex.quote(c) for c in cmd))
        if dry:
            continue
        os.makedirs(dest, exist_ok=True)
        if subprocess.run(cmd).returncode:
            die(f"bs download failed for project {project}")


def find_reads(outdir, ext):
    """[(path, project)] for every read file under outdir, however nested."""
    found = []
    for root, _, files in os.walk(outdir):
        for name in sorted(files):
            if not name.endswith(ext):
                continue
            path = os.path.join(root, name)
            rel = os.path.relpath(path, outdir).split(os.sep)
            found.append((path, rel[0] if len(rel) > 1 else ""))
    return sorted(found)


def build_groups(reads, ext):
    """{(pool, read): [(lane, path, project, basename)]}.

    A pool name seen in more than one project gets the project appended, and
    those are kept apart rather than merged.
    """
    parsed = []
    projects_per_pool = {}
    for path, project in reads:
        base = os.path.basename(path)
        pool, read, lane = parse_name(base, ext)
        if pool is None:
            skip(base, "nothing left of the name after stripping Illumina "
                       "fields")
            continue
        parsed.append((pool, read, lane, path, project, base))
        projects_per_pool.setdefault(pool, set()).add(project)

    shared = {p for p, projs in projects_per_pool.items() if len(projs) > 1}
    for pool in sorted(shared):
        where = ", ".join(pr or "top level"
                          for pr in sorted(projects_per_pool[pool]))
        warn(f"pool {pool} appears in more than one place ({where}); the "
             f"project is appended to keep them apart and they are not "
             f"merged -- check whether these are really the same pool")

    groups = {}
    for pool, read, lane, path, project, base in parsed:
        key_pool = f"{pool}-{project}" if pool in shared and project else pool
        groups.setdefault((key_pool, read), []).append(
            (lane or "", path, project, base))
    for key in groups:
        groups[key].sort()
    return groups


def free_name(outdir, base, pool, project):
    """A name at the top level that is not taken, keeping base recognisable."""
    dst = os.path.join(outdir, base)
    if not os.path.exists(dst):
        return dst
    suffix = f"-{project}" if project else "-2"
    cand = base.replace(pool, pool + suffix, 1)
    dst = os.path.join(outdir, cand)
    n = 2
    while os.path.exists(dst):
        cand = base.replace(pool, f"{pool}{suffix}-{n}", 1)
        dst = os.path.join(outdir, cand)
        n += 1
    return dst


def leave_flat(outdir, parts, dry):
    """Move an unprocessed pool's files to the top level, names unchanged."""
    for _, path, project, base in parts:
        pool = parse_name(base, READS_EXT)[0] or base
        dst = free_name(outdir, base, pool, project)
        if os.path.abspath(dst) == os.path.abspath(path):
            continue
        if os.path.basename(dst) != base:
            warn(f"  {base} already taken at the top level, moved as "
                 f"{os.path.basename(dst)}")
        if not dry:
            os.rename(path, dst)


def process(groups, outdir, ext, keep, dry, do_merge, handled):
    reads_per_pool = {}
    for pool, read in groups:
        reads_per_pool.setdefault(pool, set()).add(read)

    rows = []
    for (pool, read), parts in sorted(groups.items(),
                                      key=lambda kv: (kv[0][0],
                                                      kv[0][1] or "")):
        tag = f"_{read}" if read and len(reads_per_pool[pool]) > 1 else ""
        out = os.path.join(outdir, f"{pool}{tag}{ext}")
        paths = [p for _, p, _, _ in parts]
        handled.update(os.path.abspath(p) for p in paths)
        total = sum(os.path.getsize(p) for p in paths)
        label = f"{pool}{tag}"

        if None in reads_per_pool[pool] and len(reads_per_pool[pool]) > 1:
            skip(label, f"both {pool}{ext} and files still carrying a read "
                        f"designator are present")
            leave_flat(outdir, parts, dry)
            rows.append((label, len(paths), total, "skipped"))
            continue

        if len(paths) == 1 and os.path.abspath(paths[0]) == os.path.abspath(out):
            rows.append((label, 1, total, "already named"))
            continue

        if os.path.exists(out):
            skip(label, f"{os.path.basename(out)} already exists and is not "
                        f"one of its parts")
            leave_flat(outdir, parts, dry)
            rows.append((label, len(paths), total, "skipped"))
            continue

        if len(paths) > 1 and not do_merge:
            leave_flat(outdir, parts, dry)
            rows.append((label, len(paths), total,
                         "left as parts; pass --merge-reads to combine"))
            continue

        if dry:
            rows.append((label, len(paths), total,
                         "move" if len(paths) == 1 else "merge"))
            continue

        if len(paths) == 1:
            os.rename(paths[0], out)
            rows.append((label, 1, os.path.getsize(out), "moved"))
            continue

        with open(out, "wb") as fh:
            for p in paths:
                with open(p, "rb") as part:
                    shutil.copyfileobj(part, fh, 1024 * 1024)
        got = os.path.getsize(out)
        if got != total:
            skip(label, f"merged to {got} bytes but the parts total {total}; "
                        f"parts left in place")
            os.remove(out)
            leave_flat(outdir, parts, dry)
            rows.append((label, len(paths), total, "skipped"))
            continue
        action = f"merged {len(paths)}"
        if not keep:
            for p in paths:
                os.remove(p)
            action += ", parts removed"
        rows.append((label, len(paths), got, action))
    return rows


def flatten_rest(outdir, handled, dry):
    """Move everything else up too, so no download subdirectory survives.

    bs leaves a per-project json manifest alongside the reads; without this it
    would keep its project directory alive and leave files nested. `handled`
    is the read files already placed by process(), which in a dry run have not
    actually moved yet.
    """
    moved = 0
    for root, _, files in sorted(os.walk(outdir), reverse=True):
        if os.path.abspath(root) == os.path.abspath(outdir):
            continue
        for name in sorted(files):
            src = os.path.join(root, name)
            if os.path.abspath(src) in handled:
                continue
            dst = free_name(outdir, name, os.path.splitext(name)[0],
                            os.path.relpath(root, outdir).split(os.sep)[0])
            if os.path.basename(dst) != name:
                warn(f"  {name} already taken at the top level, moved as "
                     f"{os.path.basename(dst)}")
            if not dry:
                os.rename(src, dst)
            moved += 1
    return moved


def prune(outdir, dry):
    """Remove the now-empty download subdirectories."""
    if dry:
        return
    for root, _, _ in sorted(os.walk(outdir), reverse=True):
        if os.path.abspath(root) == os.path.abspath(outdir):
            continue
        try:
            os.rmdir(root)
        except OSError:
            pass


def main():
    p = argparse.ArgumentParser(
        prog=prog,
        description="Download pooled Illumina reads from BaseSpace, merge lane "
                    "parts without decompressing, and flatten to one file per "
                    "pool.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # download two projects and flatten to one file per pool
  {prog} -p JA23123 JA23124

  # same, and combine lane parts where a pool has them
  {prog} -p JA23123 JA23124 --merge-reads

  # also download the Undetermined reads
  {prog} -p JA26157 --undetermined

  # submit it as a single launcher job
  launcher_creator.py -b 'srun {prog} -p JA23123 JA23124' -n fetchReads \\
      -q shortq7 -t 06:00:00 -e you@fau.edu
  sbatch fetchReads.slurm

  # already downloaded; just merge and flatten what is there
  {prog} --merge-reads

  # show the plan and change nothing
  {prog} --dry-run

  # read counts afterwards, straight off the gz
  readCounts.py -f fastq.gz -o gom -j 16

what is downloaded:
  Only files ending in -e (default {READS_EXT}). A BaseSpace project also
  holds the demultiplexer's reports, logs, InterOp and sample sheet; those are
  not fetched. Undetermined reads are left out unless --undetermined is given.
  The filter is passed to bs as --exclude/--include, which needs BaseSpace CLI
  1.4.0 or later.

pools, not samples:
  A file here is a pool. Samples inside it are separated later by in-line
  barcode, so identity is pool plus barcode and two pools must never be
  combined. That is why identically named files in different projects get the
  project appended instead of being merged: they are either the same file
  twice or two different pools, and merging the second case would put two
  specimens under one barcode.

failures:
  A problem with one pool prints a reason, leaves that pool's files at the top
  level under their original Illumina names, and the run carries on. Long
  names mark them as unprocessed. The exit status is non-zero if anything was
  skipped. Only a bad --outdir, a failed download, or no read files at all
  stop the run.
""")
    p.add_argument("-p", "--project", nargs="+", metavar="NAME", default=[],
                   help="BaseSpace project name(s) to download, each into its "
                        "own subdirectory; omit to work on files already here")
    p.add_argument("-o", "--outdir", default=".", metavar="DIR",
                   help="directory to download into and flatten into "
                        "(default: cwd)")
    p.add_argument("-e", "--ext", default=READS_EXT, metavar="EXT",
                   help=f"read file extension; also the only files downloaded "
                        f"(default: {READS_EXT})")
    p.add_argument("--undetermined", action="store_true",
                   help="also download the Undetermined reads, which are left "
                        "out by default")
    p.add_argument("--merge-reads", dest="merge", action="store_true",
                   help="concatenate lane parts into one file per pool; "
                        "without it, pools with more than one part are left "
                        "unprocessed and reported")
    p.add_argument("--keep-parts", action="store_true",
                   help="keep the per-lane files after merging instead of "
                        "deleting them once the merged size checks out")
    p.add_argument("--dry-run", action="store_true",
                   help="report what would happen and change nothing")
    a = p.parse_args()

    if not os.path.isdir(a.outdir):
        die(f"{a.outdir} is not a directory")

    download(a.project, a.outdir, a.dry_run, a.ext, a.undetermined)

    reads = find_reads(a.outdir, a.ext)
    if not reads:
        die(f"no *{a.ext} in {a.outdir}")

    handled = set()
    rows = process(build_groups(reads, a.ext), a.outdir, a.ext,
                   a.keep_parts, a.dry_run, a.merge, handled)
    rest = flatten_rest(a.outdir, handled, a.dry_run)
    prune(a.outdir, a.dry_run)

    print("pool\tparts\tbytes\taction")
    for name, parts, size, action in rows:
        print(f"{name}\t{parts}\t{size}\t{action}")

    if rest:
        warn(f"moved {rest} other file(s) up out of the download "
             f"directories")
    warn(f"{len(rows)} pool(s)"
         f"{', ' + str(len(skipped)) + ' skipped' if skipped else ''}"
         f"{' (dry run, nothing changed)' if a.dry_run else ''}")
    if skipped:
        for pool, reason in skipped:
            warn(f"  {pool}: {reason}")
        sys.exit(1)


if __name__ == "__main__":
    main()
