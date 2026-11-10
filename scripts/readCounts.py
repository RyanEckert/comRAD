#!/usr/bin/env python3
################################################################################
#Script Name    : readCounts.py
#Description    : Generate a file with Illumina read counts from specified files
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, glob, gzip, math, os, sys
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Prints a tab-separated file of Illumina read counts for files with a specified extension.\n"
"Counts newlines in binary chunks and reads files in parallel, and handles gzipped input.\n"
"\n"
"Output is 'file<TAB>reads' under a header line, so it reads straight into R\n"
"or awk:\n"
"  read.delim('trimCounts')\n"
"  awk 'NR>1 {t+=$2} END {print t}' trimCounts\n"
"\n"
"usage examples:\n"
"  readCounts.py                                 # *fq in the working directory -> readCounts\n"
"  readCounts.py -f fastq -o raw                 # *fastq -> rawCounts\n"
"  readCounts.py -f fq.gz -j 8                   # gzipped reads, 8 files at a time\n"
"  readCounts.py -f trim -l 1 -o trim            # one line per read instead of four\n")
parser.add_argument("-f", type=str, help="File extension for files to add, don't need to include '.' [e.g. fastq]", action = "store", default='fq')
parser.add_argument("-o", type=str, help="Base name for Counts file [e.g. read -> readCounts]", action = "store", default='read')
parser.add_argument("-j", type=int, help="Number of files to count in parallel [default = 8]", action = "store", default=8)
parser.add_argument("-l", type=int, help="Lines per read [default = 4 for fastq; use 1 for one-line-per-read formats]", action = "store", default=4)

args = parser.parse_args()
f = args.f
o = args.o

CHUNK = 8 * 1024 * 1024


def lineCount(seqFile):
    """Count newlines in binary chunks. Roughly an order of magnitude faster
    than iterating decoded lines, and it reads .gz directly.

    Decompression stays on python's gzip module: it is zlib in C, and both
    'gzip -dc' and 'pigz -dc' through a pipe measured slower (1.7s and 1.2s
    against 0.8s on 55MB/6M lines) -- pigz cannot parallelise inflate of a
    single stream, it only overlaps io and crc.

    Returns (seqFile, count) or (seqFile, None) if the file could not be read,
    so one bad file does not take the whole run down with it."""
    opener = gzip.open if seqFile.endswith('.gz') else open
    try:
        n = 0
        last = b'\n'
        with opener(seqFile, 'rb') as fh:
            while True:
                block = fh.read(CHUNK)
                if not block:
                    break
                n += block.count(b'\n')
                last = block[-1:]
        if last != b'\n':       # file not newline-terminated
            n += 1
        return seqFile, n
    except (OSError, EOFError, gzip.BadGzipFile) as err:
        print(f"readCounts.py: {seqFile}: {err}", file=sys.stderr)
        return seqFile, None


seqFiles = sorted(glob.iglob(r'*' + f))
if not seqFiles:
    sys.exit("no files matched *" + f)

with ProcessPoolExecutor(max_workers=args.j) as pool:
    counts = list(pool.map(lineCount, seqFiles))

bad = [f for f, n in counts if n is None]

with open(o + "Counts", 'w') as outfile:
    outfile.write("file\treads\n")
    for seqFile, n in counts:
        if n is None:
            continue
        outfile.write(seqFile + '\t' + str(math.trunc(n / args.l)) + '\n')

print(f"readCounts.py: {len(counts) - len(bad)} of {len(counts)} file(s) counted"
      f" -> {o}Counts", file=sys.stderr)
if bad:
    print(f"readCounts.py: {len(bad)} file(s) failed and are absent from "
          f"{o}Counts:", file=sys.stderr)
    for f in bad:
        print('  ' + f, file=sys.stderr)
    sys.exit(1)
