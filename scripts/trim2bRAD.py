#!/usr/bin/env python3
################################################################################
#Script Name    : trim2bRAD.py
#Description    : Trim, deduplicate and demultiplex 2bRAD reads by in-line barcode
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
#
#Python port of 'trim2bRAD_2barcodes_dedup.pl' and its launcher
#'2bRAD_trim_launch_dedup.pl' by Misha Matz
#(https://github.com/z0on/2bRAD_denovo)
################################################################################
"""Trim/dedup/demultiplex 2bRAD reads, or write a job list of trim commands.

A read is kept only if it matches, from the 5' end and with no N anywhere:

    [degenerate leader][restriction site][in-line barcode][far-end adaptor]
     [ATGC][ATGC][AG][AT]   36 bp BcgI      [ATGC]{4}         AGATC?

Only the restriction fragment is written. Reads sharing a leader, the first 34
bases of the fragment, and a barcode are PCR duplicates; the first is kept.
Output is one fastq per barcode, '<sample>_<barcode>.demux.fastq' (--out-ext).

Differences from the Perl, all deliberate:

  * Quality strings line up with their bases. The Perl sliced quality from a
    hard-coded offset of 3 while taking sequence from the end of the matched
    leader -- 4 bases in the 2-barcode script, 6 in the 1-barcode one -- so
    every base carried its left neighbour's Q score, which propagated into
    bowtie2 and ANGSD genotype likelihoods. A literal that should have been
    derived from the match length.
  * --min-bc-count is inclusive. The Perl compared against a last-index rather
    than a count, so 'minBCcount=100000' actually required 100001 reads.
  * Reads are streamed to per-barcode files instead of being held in memory,
    and dedup keys are stored as 128-bit digests, so peak RAM does not scale
    with the number of reads kept.
  * gzipped input is read directly, and --gzip writes '<sample>_<bc>.demux.fastq.gz'
    at level 9, so nothing has to be compressed in a later pass.
  * The sample name comes from the basename, so a relative path like
    '../raw/x.fastq' no longer collapses the name to an empty string.
  * The file's last read is handled like every other one. The Perl's post-loop
    block dropped both the 'and $deduplicate' guard and the no-dedup branch.
  * In job-list mode the barcode pattern is passed as 'bc2', which is what the
    trimmer actually reads. The Perl launcher emitted 'bc=', which was ignored.
"""
import argparse
import glob
import gzip
import hashlib
import io
import os
import re
import shlex
import subprocess
import sys

DEGEN = "[ATGC][ATGC][AG][AT]CC"
MIN_YIELD = 0.01
SITE = ".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}"
BC2 = "[ATGC]{4}"
ADAPTOR = "AGATC?"
DUP_KEY_BASES = 34
OUT_EXT = ".demux.fastq"
LOG_EXT = ".demux.log"

prog = os.path.basename(sys.argv[0])


# --------------------------------------------------------------------- io ---
def open_reads(path):
    """Text handle for a fastq, gzipped or not. Returns (handle, proc)."""
    if not path.endswith(".gz"):
        return open(path, "r", encoding="ascii", errors="replace"), None
    for exe in ("pigz", "gzip"):
        which = subprocess.run(["which", exe], capture_output=True, text=True)
        if which.returncode == 0:
            p = subprocess.Popen([exe, "-cd", path], stdout=subprocess.PIPE)
            return io.TextIOWrapper(p.stdout, encoding="ascii",
                                    errors="replace"), p
    return gzip.open(path, "rt", encoding="ascii", errors="replace"), None


def fastq(handle, path):
    """Yield (header, seq, qual) with the 4-line structure enforced."""
    n = 0
    while True:
        head = handle.readline()
        if not head:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        n += 4
        if not qual:
            sys.exit(f"{prog}: {path} ends mid-record at line {n}")
        head, seq, qual = head.rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n")
        if not head.startswith("@") or not plus.startswith("+"):
            sys.exit(f"{prog}: {path} line {n - 3} is not a fastq header, or "
                     f"line {n - 1} is not '+'. Is this file wrapped or "
                     "truncated?")
        yield head, seq, qual


def sample_name(path, sample_id):
    """Perl's naming, but from the basename so paths with dots survive."""
    base = os.path.basename(path)
    if sample_id is None:
        return base.split(".", 1)[0]
    parts = base.split("_")
    if sample_id > len(parts):
        sys.exit(f"{prog}: --sample-id {sample_id} but '{base}' splits into "
                 f"{len(parts)} underscore-separated field(s)")
    return parts[sample_id - 1]


# ------------------------------------------------------------------ trim ---
def trim(a):
    read_re = re.compile(f"^({a.degen})({a.site})({a.bc2}){a.adaptor}")
    name = sample_name(a.input, a.sample_id)
    ext = a.out_ext + (".gz" if a.gzip else "")
    handle, proc = open_reads(a.input)

    seen = set()
    parts, counts, dups = {}, {}, {}
    total_reads = kept = 0

    def part_file(bc):
        if bc not in parts:
            path = f"{name}_{bc}{ext}.part"
            parts[bc] = (gzip.open(path, "wt", compresslevel=a.level)
                         if a.gzip else open(path, "w"))
            counts[bc] = 0
            dups[bc] = 0
        return parts[bc]

    def handle_read(head, seq, qual):
        nonlocal kept
        if "N" in seq:
            return
        m = read_re.match(seq)
        if not m:
            return
        leader, frag, bc = m.group(1), m.group(2), m.group(3)

        if a.dedup:
            read, start = frag, len(leader)
        else:
            read, start = leader + frag, 0

        if a.dedup:
            key = hashlib.blake2b(
                (leader + read[:DUP_KEY_BASES] + bc).encode("ascii"),
                digest_size=16).digest()
            if key in seen:
                part_file(bc)
                dups[bc] += 1
                return
            seen.add(key)

        stop = len(read) - a.clip
        if stop <= a.clip:
            sys.exit(f"{prog}: --clip {a.clip} removes the whole "
                     f"{len(read)} bp read")
        out_seq = read[a.clip:stop]
        out_qual = qual[start + a.clip:start + stop]

        fh = part_file(bc)
        fh.write(f"{head} bcd={bc}\n{out_seq}\n+\n{out_qual}\n")
        counts[bc] += 1
        kept += 1

    finished = False
    try:
        for head, seq, qual in fastq(handle, a.input):
            total_reads += 1
            handle_read(head, seq, qual)
        finished = True
    finally:
        for fh in parts.values():
            fh.close()
        handle.close()
        if proc is not None:
            proc.wait()
        if not finished:
            # a partial output is worse than none: it looks like a finished file
            for bc in parts:
                try:
                    os.remove(f"{name}_{bc}{ext}.part")
                except OSError:
                    pass

    keep_bc = []
    for bc in sorted(parts):
        if counts[bc] >= a.min_bc_count:
            keep_bc.append(bc)
        else:
            os.remove(f"{name}_{bc}{ext}.part")

    tot_good = tot_dup = 0
    for bc in keep_bc:
        os.replace(f"{name}_{bc}{ext}.part", f"{name}_{bc}{ext}")
        print(f"{name}_{bc}: goods: {counts[bc]} ; dups: {dups[bc]}")
        tot_good += counts[bc]
        tot_dup += dups[bc]

    dropped = [bc for bc in sorted(parts) if bc not in keep_bc]
    if dropped:
        print(f"{prog}: {len(dropped)} barcode(s) below --min-bc-count "
              f"{a.min_bc_count}, no file written:", file=sys.stderr)
        for bc in sorted(dropped, key=lambda b: -counts[b])[:10]:
            print(f"  {bc}\t{counts[bc]} reads", file=sys.stderr)

    print(f"\n{name}: total goods : {tot_good} ; total dups : {tot_dup}")
    print(f"{name}: {total_reads} reads in, {kept} kept "
          f"({100.0 * kept / total_reads:.1f}%)" if total_reads else
          f"{name}: empty input", file=sys.stderr)

    # A yield this low means the reads do not have the structure the regex
    # describes -- a wrong --degen length, the wrong enzyme, or the tag at an
    # offset other than len(degen). Random sequence matches at ~0.01%, so
    # silently reporting success here loses the whole run.
    if total_reads and kept / total_reads < MIN_YIELD:
        sys.exit(
            f"{prog}: {name}: only {kept} of {total_reads} reads matched "
            f"({100.0 * kept / total_reads:.3f}%) -- at or near the rate "
            f"random sequence matches by chance.\n"
            f"  The reads are not laid out as "
            f"[{a.degen}][{a.site}][{a.bc2}]{a.adaptor}.\n"
            f"  Check where the site motif actually falls:\n"
            f"    zcat {shlex.quote(a.input)} | awk 'NR%4==2' | head -200000 |"
            f" grep -boE 'CGA.{{6}}TGC|GCA.{{6}}TCG' | cut -d: -f1 |"
            f" sort -n | uniq -c | sort -rn | head\n"
            f"  The dominant offset must equal the leader length plus 12 "
            f"(18 for the default 6-base leader).")


# ------------------------------------------------------------- job list ---
def write_jobs(a):
    # outputs of an earlier run can match the input suffix (.fastq.gz); never trim those
    files = sorted(f for f in glob.glob(f"*{a.suffix}") if os.path.isfile(f)
                   and a.out_ext + "." not in f + "." and not f.endswith(a.out_ext))
    if not files:
        sys.exit(f"{prog}: no files matched *{a.suffix}")

    names = {}
    for f in files:
        names.setdefault(sample_name(f, a.sample_id), []).append(f)
    clash = {n: fs for n, fs in names.items() if len(fs) > 1}
    if clash:
        msg = "\n".join(f"  {n}: {', '.join(fs)}" for n, fs in
                        sorted(clash.items()))
        sys.exit(f"{prog}: these files resolve to the same sample name and "
                 f"would overwrite each other's output files:\n{msg}\n"
                 "Concatenate them first, or change --sample-id.")

    if a.logdir:
        os.makedirs(a.logdir, exist_ok=True)

    self_path = os.path.abspath(sys.argv[0]) if os.sep in sys.argv[0] else prog
    with open(a.name, "w") as fh:
        for f in files:
            cmd = [self_path, "-i", f,
                   "--degen", a.degen, "--site", a.site, "--bc2", a.bc2,
                   "--adaptor", a.adaptor, "--clip", str(a.clip),
                   "--min-bc-count", str(a.min_bc_count)]
            if a.sample_id is not None:
                cmd += ["--sample-id", str(a.sample_id)]
            if a.gzip:
                cmd += ["--gzip", "--level", str(a.level)]
            if a.out_ext != OUT_EXT:
                cmd += ["--out-ext", a.out_ext]
            if not a.dedup:
                cmd.append("--no-dedup")
            line = " ".join(shlex.quote(c) for c in cmd)
            if a.logdir:
                log = sample_name(f, a.sample_id) + a.log_ext
                line += f" > {shlex.quote(os.path.join(a.logdir, log))} 2>&1"
            fh.write(line + "\n")
    print(f"{a.name}: {len(files)} jobs", file=sys.stderr)

    if a.launcher:
        os.system("module load launcher")
        cmd = f"launcher_creator.py -j {a.name} -n {a.name}"
        if a.nodes:
            cmd += f" -N {a.nodes}"
        cmd += f" -q {a.queue} -t {a.time}"
        if a.wayness:
            cmd += f" -w {a.wayness}"
        if a.email:
            cmd += f" -e {a.email}"
        os.system(cmd)


# ------------------------------------------------------------------ main ---
def main():
    p = argparse.ArgumentParser(
        description="Trim, deduplicate and demultiplex 2bRAD reads by in-line "
                    "barcode, or write a list of trim commands to run in "
                    "parallel.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # trim one file; sample name is the first underscore field
  {prog} -i Sample_L001_R1.fastq --sample-id 1

  # gzipped in and out, keep any barcode with at least 1000 reads
  {prog} -i Sample_L001_R1.fastq.gz --sample-id 1 --min-bc-count 1000 --gzip

  # AlfI instead of BcgI
  {prog} -i Sample_L001_R1.fastq --site '.{{12}}GCA.{{6}}TGC.{{12}}'

  # clip 2 bases off each end, keep PCR duplicates
  {prog} -i Sample_L001_R1.fastq --clip 2 --no-dedup

  # write one trim command per .fastq here, then submit them
  {prog} -f fastq --sample-id 1 -n trims --launcher -e {'you@fau.edu'}

  # write the commands file only, submit it yourself
  {prog} -f .fastq.gz --sample-id 1 -n trims --gzip

  # collect the per-sample yields afterwards
  grep 'total goods' logs/*.demux.log
""")

    mode = p.add_argument_group("mode (give exactly one)")
    mode.add_argument("-i", "--input", metavar="FASTQ",
                      help="trim this one fastq (plain or .gz)")
    mode.add_argument("-f", "--suffix", metavar="EXT",
                      help="write one trim command per file ending in EXT in "
                           "the working directory, e.g. 'fastq' or '.fastq.gz'")

    trm = p.add_argument_group("trimming")
    trm.add_argument("--site", default=SITE, metavar="RE",
                     help="restriction site pattern (default: BcgI)")
    trm.add_argument("--bc2", default=BC2, metavar="RE",
                     help=f"in-line barcode following the fragment "
                          f"(default: {BC2})")
    trm.add_argument("--adaptor", default=ADAPTOR, metavar="RE",
                     help=f"far-end adaptor (default: {ADAPTOR})")
    trm.add_argument("--degen", default=DEGEN, metavar="RE",
                     help=f"degenerate 5' leader (default: {DEGEN})")
    trm.add_argument("--clip", type=int, default=0, metavar="N",
                     help="bases to remove from each end of the fragment "
                          "(default: 0)")
    trm.add_argument("--sample-id", type=int, metavar="N",
                     help="1-based position of the name-deriving field in the "
                          "underscore-separated filename; default is the "
                          "filename up to its first dot")
    trm.add_argument("--min-bc-count", type=int, default=100000, metavar="N",
                     help="write a barcode's file only if it has at least N "
                          "reads (default: 100000)")
    trm.add_argument("--gzip", action="store_true",
                     help="gzip the output (<out-ext>.gz); cutadapt reads it directly")
    trm.add_argument("--out-ext", default=OUT_EXT, metavar="EXT",
                     help=f"output name ending after <sample>_<barcode> (default: {OUT_EXT}; "
                          f"'.tr0' gives the 2bRAD_denovo names)")
    trm.add_argument("--level", type=int, default=5, choices=range(1, 10),
                     metavar="1-9",
                     help="gzip compression level for --gzip (default: 5). "
                          "9 costs about 6x the CPU of 5 for ~5%% less size")
    trm.add_argument("--no-dedup", dest="dedup", action="store_false",
                     help="keep PCR duplicates; output then retains the "
                          "degenerate leader, as in the Perl")

    lch = p.add_argument_group("job list")
    lch.add_argument("-n", "--name", default="trims", metavar="FILE",
                     help="commands file to write (default: trims)")
    lch.add_argument("--logdir", default="logs", metavar="DIR",
                     help="per-file logs, <sample><log-ext> (default: logs), '' to disable")
    lch.add_argument("--log-ext", default=LOG_EXT, metavar="EXT",
                     help=f"log name ending (default: {LOG_EXT})")
    lch.add_argument("--launcher", action="store_true",
                     help="also build the .slurm with launcher_creator.py")
    lch.add_argument("-N", "--nodes", metavar="N", help="nodes for launcher")
    lch.add_argument("-w", "--wayness", metavar="N",
                     help="commands per node for launcher")
    lch.add_argument("-q", "--queue", default="shortq7", help="SLURM partition")
    lch.add_argument("-t", "--time", default="6:00:00", help="walltime")
    lch.add_argument("-e", "--email", metavar="ADDR",
                     help="address for SLURM notifications")

    a = p.parse_args()

    if bool(a.input) == bool(a.suffix):
        p.error("give either -i (trim one file) or -f (write a job list)")
    if a.clip < 0:
        p.error("--clip must be >= 0")
    if a.sample_id is not None and a.sample_id < 1:
        p.error("--sample-id is 1-based")
    for name in ("site", "bc2", "adaptor", "degen"):
        try:
            re.compile(getattr(a, name))
        except re.error as err:
            p.error(f"--{name}: {err}")

    if a.input:
        if not os.path.isfile(a.input):
            sys.exit(f"{prog}: cannot open {a.input}")
        trim(a)
    else:
        write_jobs(a)


if __name__ == "__main__":
    main()
