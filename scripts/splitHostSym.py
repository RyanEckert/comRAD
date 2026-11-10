#!/usr/bin/env python3
################################################################################
#Script Name    : splitHostSym.py
#Description    : Split 2bRAD alignments into disjoint host and symbiont sets
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Separate host (coral) from symbiont 2bRAD reads.

Input is the two alignment files produced by mapping the SAME trimmed reads
independently against a host genome and a symbiont (meta)genome. Neither pass
may filter the other's input:

    bowtie2Launcher.py -g $SYMGENOME  -f .trim -n symMaps -s sym --bam \\
        -w 24 --launcher -e you@fau.edu
    bowtie2Launcher.py -g $HOSTGENOME -f .trim -n hostMaps -s host --bam \\
        --launcher -e you@fau.edu

A read that aligned to only one reference is written to that reference's
output. A read that aligned to both is written to neither, so the two outputs
are disjoint. The test is the read name: a read counts as aligned to both when
the same QNAME carries an aligned record in both files. No alignment score is
read -- a read mapping to both is ambiguous however lopsided the scores, and
which side would win a close race is partly decided by sequencing error. What
counts as "aligned" is therefore set entirely by bowtie2's --score-min at
mapping time.

Nothing is ever re-mapped; records are only dropped from files that exist.

Two subcommands, because masking is a reduction over ALL samples:

  scan    read every sample's pair and write the cross-mapping loci as one BED
          per sample per reference. No alignments are written.
  split   write the partitioned output. --bed additionally masks out the loci
          scan found, in every sample.

Run 'split' alone to just partition. Run 'scan' then 'split --bed' to also
mask. The two cannot be one command: split for one sample needs the BEDs from
every other sample, and launcher tasks are independent processes that cannot
wait for each other.

Why mask at all. --score-min is a hard cutoff and the verdict is per read, so
reads at one 2bRAD tag can differ. Sequencing error moves a marginal symbiont
alignment across the threshold, and -- more importantly -- so does genotype: if
one allele at a tag is a base closer to the symbiont, reads carrying it
cross-map and reads carrying the other allele do not. At a heterozygous
cross-mappable tag that discards one allele and leaves the other, so the
individual is called a homozygote. Not missing data, a wrong genotype, produced
deterministically within a single sample. Masking the locus replaces that with
no call. The cost is real: masking also drops reads that were unambiguous in
their own sample, and the BEDs only mark loci where some sample actually
produced a cross-mapping read, so it is a floor on the problem rather than a
complete mappability filter.

The host and symbiont BEDs are two independent sets in unrelated coordinate
systems -- there is no correspondence between symbiont contig positions and
host genome positions -- so a host BED is never applied to symbiont output.
"""
import argparse
import bisect
import glob
import os
import re
import shlex
import subprocess
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
SAM_EXTS = (".sam.gz", ".sam", ".bam", ".cram")

prog = os.path.basename(sys.argv[0])


# --------------------------------------------------------------------- io ---
def open_sam(path):
    """Text handle for a SAM/SAM.gz/BAM. Returns (handle, proc)."""
    if path.endswith((".bam", ".cram")):
        p = subprocess.Popen(["samtools", "view", "-h", path],
                             stdout=subprocess.PIPE, text=True)
        return p.stdout, p
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rt"), None
    return open(path, "rt"), None


def close_sam(handle, proc):
    handle.close()
    if proc is not None and proc.wait() != 0:
        sys.exit(f"{prog}: samtools failed reading input")


class SamWriter:
    """Write SAM text out as BAM, optionally coordinate sorted and indexed.

    The script only ever emits SAM lines; samtools does the compressing, and
    the sorting when asked, so no SAM is ever written to disk.
    """

    def __init__(self, path, sort=False, level=9):
        self.path = path
        self.sort = sort
        self.proc = None
        self.handle = None
        self.n = 0          # alignment records written, headers excluded
        if path is None:
            return
        cmd = (["samtools", "sort", "-l", str(level), "-o", path, "-"] if sort
               else ["samtools", "view", "-b", "--output-fmt-option",
                     f"level={level}", "-o", path, "-"])
        self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        self.handle = self.proc.stdin

    def write_header(self, line):
        if self.handle is not None:
            self.handle.write(line)

    def write(self, line):
        self.n += 1
        if self.handle is not None:
            self.handle.write(line)

    def close(self):
        if self.handle is not None:
            self.handle.close()
        if self.proc is not None and self.proc.wait() != 0:
            sys.exit(f"{prog}: samtools failed writing {self.path}")
        if self.path and self.sort:
            if subprocess.run(["samtools", "index", self.path]).returncode:
                sys.exit(f"{prog}: samtools index failed on {self.path}")


def is_aligned(fields):
    return not int(fields[1]) & 0x4


def ref_span(fields):
    """(chrom, bed_start, bed_end) covered by one aligned SAM record."""
    length = 0
    for n, op in CIGAR.findall(fields[5]):
        if op in "MDN=X":
            length += int(n)
    start = int(fields[3]) - 1
    return fields[2], start, start + length


def aligned_read_names(path):
    """Read names with at least one aligned record."""
    names = set()
    handle, proc = open_sam(path)
    for line in handle:
        if line.startswith("@"):
            continue
        fields = line.split("\t", 2)
        if not int(fields[1]) & 0x4:
            names.add(fields[0])
    close_sam(handle, proc)
    return names


# -------------------------------------------------------------------- bed ---
def write_bed(path, spans):
    with open(path, "wt") as fh:
        for chrom, start, end in sorted(spans):
            fh.write(f"{chrom}\t{start}\t{end}\n")


def load_bed(paths):
    """Read BED files into {chrom: (starts, ends)}, intervals merged.

    Inputs need not be sorted, merged or deduplicated, and any number of files
    can be pooled in one call, so no cat/sort/bedtools merge step is needed.
    """
    raw = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track",
                                                        "browser")):
                    continue
                f = line.split()
                raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    merged = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        keep = []
        for start, end in intervals:
            if keep and start <= keep[-1][1]:
                if end > keep[-1][1]:
                    keep[-1][1] = end
            else:
                keep.append([start, end])
        merged[chrom] = ([iv[0] for iv in keep], [iv[1] for iv in keep])
    return merged


def bed_size(bed):
    """(intervals, bases) covered by a loaded BED."""
    return (sum(len(v[0]) for v in bed.values()),
            sum(e - s for v in bed.values() for s, e in zip(*v)))


def overlaps(bed, fields):
    """True if this aligned record touches any merged interval."""
    entry = bed.get(fields[2])
    if entry is None:
        return False
    starts, ends = entry
    _, span_start, span_end = ref_span(fields)
    i = bisect.bisect_right(starts, span_start) - 1
    if i >= 0 and ends[i] > span_start:
        return True
    i += 1
    return i < len(starts) and starts[i] < span_end


def mask_beds(beddir, side, n_expected):
    """Every <beddir>/*.<side>.bed, with a warning if the count looks short.

    scan writes one BED per sample per side, so a count below the number of
    pairs being split means scan has not finished, or some of its tasks died.
    Masking against a partial set silently leaves cross-mapping loci in.
    """
    paths = sorted(glob.glob(os.path.join(beddir, f"*.{side}.bed")))
    if not paths:
        sys.exit(f"{prog}: no *.{side}.bed in {beddir} -- run scan first")
    if n_expected and len(paths) < n_expected:
        print(f"{prog}: warning: masking against {len(paths)} {side} BEDs but "
              f"splitting {n_expected} samples. If scan has not finished, or "
              f"some of its tasks failed, cross-mapping loci will be left "
              f"unmasked.", file=sys.stderr)
    return load_bed(paths)


# ------------------------------------------------------------------ paths ---
def resolve(args, what):
    """{sample name: file} for one quoted glob pattern.

    The sample name is the file name with the literal text after the pattern's
    last '*' removed: '*.trim.host.bt2.bam' names XXX13.1.trim.host.bt2.bam
    'XXX13.1'. Nothing about dots or sample naming is assumed.
    """
    example = f"{what} '*.trim.host.bt2.bam'"
    if len(args) != 1:
        sys.exit(f"{prog}: {what} got {len(args)} arguments. Quote the pattern "
                 f"so the shell does not expand it, e.g. {example}")
    pattern = args[0]
    fname = os.path.basename(pattern)
    if "*" not in fname:
        sys.exit(f"{prog}: {what} {pattern} has no '*' in its file name. The "
                 f"'*' stands for the sample name, e.g. {example}")
    suffix = fname.rsplit("*", 1)[1]
    if any(c in suffix for c in "?["):
        sys.exit(f"{prog}: {what} {pattern}: the text after the last '*' must "
                 f"be literal, since it is what gets removed to leave the "
                 f"sample name")
    hits = sorted(f for f in glob.glob(pattern) if os.path.isfile(f))
    if not hits:
        sys.exit(f"{prog}: {what} pattern matched nothing: {pattern}")
    files = {}
    for f in hits:
        name = os.path.basename(f)[:len(os.path.basename(f)) - len(suffix)]
        if not name:
            sys.exit(f"{prog}: {what} {f} leaves an empty sample name once "
                     f"'{suffix}' is removed")
        files.setdefault(name, []).append(f)
    dups = {n: fs for n, fs in files.items() if len(fs) > 1}
    if dups:
        msg = "\n".join(f"  {n}: {', '.join(fs)}"
                        for n, fs in sorted(dups.items()))
        sys.exit(f"{prog}: {what} {pattern} gives the same sample name to "
                 f"more than one file, and their outputs would overwrite each "
                 f"other:\n{msg}")
    return {n: fs[0] for n, fs in files.items()}


def in_dir(directory, name):
    """Join, treating '' and '.' as the working directory."""
    if not directory or directory == ".":
        return name
    return os.path.join(directory, name)


def make_dirs(*dirs):
    for d in dirs:
        if d and d != ".":
            os.makedirs(d, exist_ok=True)


def pair_up(hosts, syms):
    """[(name, host file, sym file)] by sample name, and the unpaired names.

    A name on only one side is a sample whose other mapping is missing. It is
    reported, not fatal, so the rest of the samples still run.
    """
    pairs = [(n, hosts[n], syms[n]) for n in sorted(hosts) if n in syms]
    unpaired = ([(n, "--host", hosts[n]) for n in sorted(hosts) if n not in syms]
                + [(n, "--sym", syms[n]) for n in sorted(syms) if n not in hosts])
    return pairs, unpaired


# ------------------------------------------------------------------- work ---
def run_pair(a, name, host_sam, sym_sam, n_pairs=0):
    """scan or split one pair. Both walk the same two files the same way."""
    scanning = a.cmd == "scan"

    if scanning:
        host_out = sym_out = None
        host_bed = in_dir(a.outdir, name + ".host.bed")
        sym_bed = in_dir(a.outdir, name + ".sym.bed")
    else:
        host_dir = a.outdir if a.outdir else os.path.dirname(host_sam)
        sym_dir = a.outdir if a.outdir else os.path.dirname(sym_sam)
        # the host set is the primary analysis product, so it keeps the bare
        # sample name; only the symbiont side is tagged. The BEDs stay tagged
        # on both sides -- they are internal plumbing, and --bed needs a clean
        # glob per side.
        host_out = in_dir(host_dir, name + ".bam")
        sym_out = in_dir(sym_dir, name + ".sym.bam")
        host_bed = sym_bed = None

    given = {os.path.abspath(host_sam), os.path.abspath(sym_sam)}
    for path in (host_out, sym_out, host_bed, sym_bed):
        if path and os.path.abspath(path) in given:
            sys.exit(f"{prog}: output {path} would overwrite an input; use "
                     f"--outdir")

    host_mask = mask_beds(a.bed, "host", n_pairs) if getattr(a, "bed", "") \
        else {}
    sym_mask = mask_beds(a.bed, "sym", n_pairs) if getattr(a, "bed", "") \
        else {}

    # Held in memory: read names aligning to the symbiont reference (the
    # smaller side of a coral 2bRAD library), plus those aligning to both.
    sym_names = aligned_read_names(sym_sam)

    both = set()
    host_spans = set()
    host_reads = host_masked = 0
    prev = None
    out = SamWriter(host_out, getattr(a, "sort", False), a.level)
    handle, proc = open_sam(host_sam)
    for line in handle:
        if line.startswith("@"):
            out.write_header(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if not is_aligned(fields):
            continue
        read = fields[0]
        if read != prev:
            # bowtie2 emits a read's records contiguously, so unique reads can
            # be counted by name change rather than held in a set.
            host_reads += 1
            prev = read
        if read in sym_names:
            both.add(read)
            if host_bed:
                host_spans.add(ref_span(fields))
            continue
        if host_mask and overlaps(host_mask, fields):
            host_masked += 1
            continue
        out.write(line)
    close_sam(handle, proc)
    out.close()
    host_written = out.n

    sym_spans = set()
    sym_masked = 0
    out = SamWriter(sym_out, getattr(a, "sort", False), a.level)
    handle, proc = open_sam(sym_sam)
    for line in handle:
        if line.startswith("@"):
            out.write_header(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if not is_aligned(fields):
            continue
        if fields[0] in both:
            if sym_bed:
                sym_spans.add(ref_span(fields))
            continue
        if sym_mask and overlaps(sym_mask, fields):
            sym_masked += 1
            continue
        out.write(line)
    close_sam(handle, proc)
    out.close()

    if host_bed:
        write_bed(host_bed, host_spans)
    if sym_bed:
        write_bed(sym_bed, sym_spans)

    if scanning:
        host_loci, host_bp = bed_size(load_bed([host_bed]))
        sym_loci, sym_bp = bed_size(load_bed([sym_bed]))
        if a.header:
            print("sample\thostAligned\tsymAligned\tbothRefs\thostLoci\t"
                  "hostBp\tsymLoci\tsymBp")
        print(f"{name}\t{host_reads}\t{len(sym_names)}\t{len(both)}\t"
              f"{host_loci}\t{host_bp}\t{sym_loci}\t{sym_bp}")
        return

    mask_loci, mask_bp = 0, 0
    for bed in (host_mask, sym_mask):
        loci, bp = bed_size(bed)
        mask_loci += loci
        mask_bp += bp
    if a.header:
        print("sample\thostAligned\tsymAligned\tbothRefs\thostWritten\t"
              "symWritten\thostMasked\tsymMasked\tmaskLoci\tmaskBp")
    print(f"{name}\t{host_reads}\t{len(sym_names)}\t{len(both)}\t"
          f"{host_written}\t{out.n}\t{host_masked}\t{sym_masked}\t"
          f"{mask_loci}\t{mask_bp}")


# --------------------------------------------------------------- job list ---
def self_path():
    return os.path.abspath(sys.argv[0]) if os.sep in sys.argv[0] else prog


def job_line(a, name):
    """One self-contained command: the same patterns, narrowed by --sample.

    The patterns are passed rather than the two files so each job derives the
    sample name exactly as the job list did.
    """
    cmd = [self_path(), a.cmd, "-H", a.host[0], "-S", a.sym[0],
           "--sample", name]
    if a.outdir:
        cmd += ["--outdir", a.outdir]
    if getattr(a, "bed", ""):
        cmd += ["--bed", a.bed]
    if getattr(a, "sort", False):
        cmd.append("--sort")
    if a.level != 9:
        cmd += ["--level", str(a.level)]
    line = " ".join(shlex.quote(c) for c in cmd)
    if a.logdir:
        log = os.path.join(a.logdir, name + f".{a.cmd}log")
        line += f" > {shlex.quote(log)} 2>&1"
    return line


def write_jobs(a, pairs):
    make_dirs(a.logdir, a.outdir)
    with open(a.name, "w") as fh:
        for name, _, _ in pairs:
            fh.write(job_line(a, name) + "\n")
    print(f"{a.name}: {len(pairs)} jobs", file=sys.stderr)

    if not a.launcher:
        return
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


# ------------------------------------------------------------------ args ---
def add_inputs(s):
    s.add_argument("-H", "--host", nargs="+", required=True,
                   metavar="PATTERN",
                   help="quoted glob of the host-genome alignments. The '*' "
                        "stands for the sample name: '*.trim.host.bt2.bam' "
                        "names XXX13.1.trim.host.bt2.bam 'XXX13.1'")
    s.add_argument("-S", "--sym", nargs="+", required=True, metavar="PATTERN",
                   help="quoted glob of the symbiont alignments of the same "
                        "reads, paired with -H by sample name")
    s.add_argument("--sample", default="", metavar="NAME",
                   help="run only this sample now instead of writing a job "
                        "list; the job list uses it for each line")
    s.add_argument("--level", type=int, default=9, metavar="N",
                   help="bam compression level, 0-9 (default: %(default)s). "
                        "9 keeps the output as small as htslib can make it, "
                        "which matters when these are the files you keep")
    s.add_argument("--header", action="store_true",
                   help="print the TSV column names before the summary")


def add_jobs(s, default_name):
    g = s.add_argument_group("job list")
    g.add_argument("-n", "--name", default=default_name, metavar="FILE",
                   help=f"commands file to write when more than one pair is "
                        f"given (default: {default_name})")
    g.add_argument("--logdir", default="logs", metavar="DIR",
                   help="per-job logs, which is where the TSV summaries land "
                        "(default: logs), '' to disable")
    g.add_argument("--launcher", action="store_true",
                   help="also build the .slurm with launcher_creator.py")
    g.add_argument("-N", "--nodes", metavar="N", help="nodes for launcher")
    g.add_argument("-w", "--wayness", metavar="N",
                   help="commands per node for launcher")
    g.add_argument("-q", "--queue", default="shortq7", help="SLURM partition")
    g.add_argument("-t", "--time", default="1:00:00", help="walltime")
    g.add_argument("-e", "--email", metavar="ADDR",
                   help="address for SLURM notifications")


def parse_args():
    p = argparse.ArgumentParser(
        prog=prog,
        description="Split host and symbiont 2bRAD alignments of the same "
                    "trimmed reads into two disjoint sets, dropping every read "
                    "that aligned to both references.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""the two subcommands:
  scan    find the cross-mapping loci, all samples, no alignments written
  split   write the partitioned bams, optionally masking what scan found

  They cannot be one command: masking is a reduction over ALL samples, so
  split for one sample needs every other sample's BEDs, and launcher tasks are
  independent processes that cannot wait for each other.

examples:
  # partition only, every sample, sorted and indexed
  {prog} split -H '*.trim.host.bt2.sam' -S 'symbiont/*.trim.sym.bt2.sam' \\
      --outdir final --sort --launcher -e you@fau.edu

  # the default path: find cross-mapping loci, then split with them masked
  {prog} scan  -H '*.trim.host.bt2.sam' -S 'symbiont/*.trim.sym.bt2.sam' \\
      --outdir bed --launcher -e you@fau.edu
  {prog} split -H '*.trim.host.bt2.sam' -S 'symbiont/*.trim.sym.bt2.sam' \\
      --bed bed --outdir final --sort --launcher -e you@fau.edu

  # one sample, right now
  {prog} split -H '*.trim.host.bt2.bam' -S '*.trim.sym.bt2.bam' \\
      --sample XXX13.1

  # collect the summaries
  cat logs/*.scanlog  | sort -k1,1
  cat logs/*.splitlog | sort -k1,1

sample names:
  The sample name is whatever the '*' in the pattern covers: the file name
  with the literal text after the last '*' removed. With -H
  '*.trim.host.bt2.bam',
    OGM001.trim.host.bt2.bam  -> OGM001
    XXX13.1.trim.host.bt2.bam -> XXX13.1
  The patterns must be quoted. If the shell expands them the script never sees
  the suffix, and it stops with an error rather than guess.

output names:
  OGM001 -> OGM001.bam (+ .bai with --sort) and OGM001.sym.bam
  scan writes OGM001.host.bed and OGM001.sym.bed. The host bams keep the bare
  sample name because they are the primary analysis product; only the
  symbiont side is tagged. List them with
    ls *.bam | grep -v '\\.sym\\.' > hostBams
    ls *.sym.bam > symBams

pairing:
  -H and -S are paired by sample name. A name matched by only one pattern is
  listed and skipped, the rest still run, and the exit status is non-zero. Two
  files giving the same name is an error, since their outputs would overwrite
  each other. A single pair, or --sample, runs immediately; more than one pair
  writes a job list.

before running either:
  Both mappings must have finished for every sample. A bowtie2 task that died
  mid-write leaves reads missing from its SAM, and a read missing from the
  symbiont SAM is indistinguishable from a host-only read. Check with
    ls *.trim | wc -l
    grep -lc 'overall alignment rate' logs/*.sym.bt2log  | wc -l
    grep -lc 'overall alignment rate' logs/*.host.bt2log | wc -l
""")
    subs = p.add_subparsers(dest="cmd", required=True, metavar="COMMAND")

    s = subs.add_parser(
        "scan", help="find the cross-mapping loci and write them as BEDs",
        description="For every pair, write the intervals covered by reads that "
                    "aligned to both references, as one BED per sample per "
                    "reference. No alignments are written. Feed the directory "
                    "to 'split --bed'.")
    add_inputs(s)
    s.add_argument("--outdir", default="bed", metavar="DIR",
                   help="directory for the BEDs (default: %(default)s)")
    add_jobs(s, "scanCmds")

    s = subs.add_parser(
        "split", help="write the two disjoint bams, optionally masked",
        description="For every pair, write the host-only and symbiont-only "
                    "alignments. Reads that aligned to both references are "
                    "dropped from both outputs.")
    add_inputs(s)
    s.add_argument("--outdir", default="", metavar="DIR",
                   help="directory for the output bams (default: beside each "
                        "input)")
    s.add_argument("--bed", default="", metavar="DIR",
                   help="directory of scan's BEDs; every sample's are pooled "
                        "and merged, and reads touching those loci are dropped "
                        "from every sample. Omit to partition without masking")
    s.add_argument("--sort", action="store_true",
                   help="coordinate sort and index the output instead of "
                        "leaving it in input order")
    add_jobs(s, "splitCmds")

    return p.parse_args()


def main():
    a = parse_args()
    pairs, unpaired = pair_up(resolve(a.host, "-H"), resolve(a.sym, "-S"))

    if a.sample:
        pairs = [p for p in pairs if p[0] == a.sample]
        if not pairs:
            sys.exit(f"{prog}: no host and symbiont pair named {a.sample}")
        unpaired = []
    elif not pairs:
        sys.exit(f"{prog}: no sample name is matched by both -H and -S")

    if unpaired:
        print(f"{prog}: {len(unpaired)} sample(s) matched by only one pattern "
              f"are skipped:", file=sys.stderr)
        for name, side, f in unpaired:
            print(f"  {name}\tonly {side}\t{f}", file=sys.stderr)

    if len(pairs) > 1:
        if getattr(a, "bed", ""):
            mask_beds(a.bed, "host", len(pairs))
            mask_beds(a.bed, "sym", len(pairs))
        write_jobs(a, pairs)
    else:
        make_dirs(a.outdir)
        run_pair(a, *pairs[0])

    if unpaired:
        sys.exit(1)


if __name__ == "__main__":
    main()
