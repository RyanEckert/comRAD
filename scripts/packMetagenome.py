#!/usr/bin/env python3
################################################################################
#Script Name : packMetagenome.py
#Description : Concatenate genomes into a metagenome of balanced, named chromosomes
#Author       : Ryan J Eckert
#Email         : ryan.j.eckert@gmail.com
################################################################################
"""
Concatenate multiple genomes into a metagenome where:
 - Each genome is split into a user-specified number of chromosomes.
 - Contigs are not split; each belongs entirely to one chromosome.
 - Chromosomes are balanced in size (~equal total length).
 - Chromosomes are named per genome, so alignment counts can be grouped by
   genome (or genus) from the chromosome name alone.
"""
import argparse, sys, os, gzip
from pathlib import Path

prog = os.path.basename(sys.argv[0])


def openMaybeGz(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path)


def parseFasta(fp):
    """Simple FASTA reader yielding (header, sequence)."""
    header = None
    seqLines = []
    for line in fp:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seqLines)
            header = line[1:].split()[0]
            seqLines = []
        else:
            seqLines.append(line)
    if header is not None:
        yield header, "".join(seqLines)


class FastaWriter:
    """Write a record incrementally, wrapped, without holding the whole sequence."""

    def __init__(self, fp, width=60):
        self.fp = fp
        self.width = width
        self.buf = ""

    def start(self, header):
        self.fp.write(">"+header+"\n")
        self.buf = ""

    def add(self, seq):
        self.buf += seq
        n = len(self.buf) - (len(self.buf) % self.width)
        if n:
            chunk, self.buf = self.buf[:n], self.buf[n:]
            for i in range(0, n, self.width):
                self.fp.write(chunk[i:i+self.width]+"\n")

    def end(self):
        if self.buf:
            self.fp.write(self.buf+"\n")
        self.buf = ""


def balanceContigs(contigs, k):
    """Distribute contigs into k bins of roughly equal total length.

    Greedy largest-contig-to-shortest-bin. Contigs are never split.
    Returns a list of bins, each a list of (contigId, seq).
    """
    bins = [[] for _ in range(k)]
    sizes = [0]*k
    contigs = sorted(contigs, key=lambda x: len(x[1]), reverse=True)
    for cid, seq in contigs:
        idx = min(range(k), key=lambda i: sizes[i])
        bins[idx].append((cid, seq))
        sizes[idx] += len(seq)
    return bins


parser = argparse.ArgumentParser(
    description="Pack genomes into a metagenome of balanced, per-genome named chromosomes.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='''chromosome naming:
  Chromosomes are named <prefix><number>, e.g. Sym001, Sym002, Brev001. The
  prefix comes from --prefixes, one per input genome in the same order as -i.
  Two inputs may share a prefix (two Cladocopium genomes, say); numbering then
  continues across them rather than restarting. Without --prefixes the input
  file stem is used, and with --plain-names the old chr1, chr2, ... naming is
  kept.

  The point of the naming is that samtools idxstats output can be summed per
  genome by name alone. The map TSV carries the same grouping in its genome and
  prefix columns, so either can be used downstream.

outputs:
  <output>   metagenome FASTA, one record per chromosome
  <map>      TSV: chr, prefix, genome, contig_id, contig_length,
             contig_start_on_chr, contig_end_on_chr

examples:
  # four Symbiodiniaceae genera, five chromosomes each
  {0} -i Symbiodinium.fasta Breviolum.fasta Cladocopium.fasta Durusdinium.fasta \\
      --prefixes Sym Brev Clad Dur --chrs-per-genome 5 \\
      -o symMetaGenome.fasta -m symMetaGenome.map.tsv

  # two Cladocopium assemblies share a prefix and number consecutively
  {0} -i CladA.fa.gz CladB.fa.gz --prefixes Clad Clad --chrs-per-genome 2 \\
      -o clad.fasta -m clad.map.tsv

  # then index it
  bowtie2-build symMetaGenome.fasta symMetaGenome.fasta
  samtools faidx symMetaGenome.fasta
'''.format(prog))
parser.add_argument("-i", "--inputs", nargs="+", required=True,
                    help="input FASTA files, one per genome, .gz ok (order = chromosome order)")
parser.add_argument("-o", "--output", required=True, help="output metagenome FASTA")
parser.add_argument("-m", "--map", required=True, help="output TSV mapping contigs to chromosomes")
parser.add_argument("--prefixes", nargs="+",
                    help="chromosome name prefix per input, same order as -i [e.g. Sym Brev Clad Dur]")
parser.add_argument("--chrs-per-genome", type=int, required=True, dest='chrsPerGenome',
                    help="number of chromosomes to split each genome into (>=1)")
parser.add_argument("--pad", type=int, default=3,
                    help="zero-padding width for the chromosome number")
parser.add_argument("--plain-names", action="store_true", dest='plainNames',
                    help="name chromosomes chr1, chr2, ... instead of by prefix")
parser.add_argument("--gap-seq", default="N"*50, dest='gapSeq',
                    help="gap sequence inserted between contigs")
args = parser.parse_args()

if args.chrsPerGenome < 1:
    sys.exit(prog+": --chrs-per-genome must be >= 1")
if args.prefixes and len(args.prefixes) != len(args.inputs):
    sys.exit(prog+": --prefixes has "+str(len(args.prefixes))+" values but -i has "
             + str(len(args.inputs))+" files; they must correspond")
if args.plainNames and args.prefixes:
    sys.exit(prog+": --plain-names and --prefixes are mutually exclusive")

counters = {}
chrIndex = 1
totals = {}

with open(args.output, "w") as outFasta, open(args.map, "w") as mapF:
    writer = FastaWriter(outFasta)
    mapF.write("\t".join(["chr", "prefix", "genome", "contig_id", "contig_length",
                          "contig_start_on_chr", "contig_end_on_chr"])+"\n")

    for gi, genomeFp in enumerate(args.inputs):
        genomeName = os.path.basename(genomeFp)
        while True:
            stem, ext = os.path.splitext(genomeName)
            if ext.lower() in ('.gz', '.bz2', '.fa', '.fas', '.fsa', '.fna', '.fasta'):
                genomeName = stem
            else:
                break
        prefix = args.prefixes[gi] if args.prefixes else genomeName

        with openMaybeGz(genomeFp) as fh:
            contigs = [(hid, seq) for hid, seq in parseFasta(fh)]
        if not contigs:
            sys.stderr.write("warning: no sequences in "+genomeFp+", skipped\n")
            continue
        if len(contigs) < args.chrsPerGenome:
            sys.stderr.write("warning: "+genomeFp+" has "+str(len(contigs))
                             + " contigs but --chrs-per-genome is "
                             + str(args.chrsPerGenome)+"; empty chromosomes are dropped\n")

        bins = balanceContigs(contigs, args.chrsPerGenome)

        for b in bins:
            if not b:
                continue
            if args.plainNames:
                chrName = "chr"+str(chrIndex)
            else:
                counters[prefix] = counters.get(prefix, 0)+1
                chrName = prefix+str(counters[prefix]).zfill(args.pad)
            chrIndex += 1

            writer.start(chrName)
            pos = 1
            for i, (cid, seq) in enumerate(b):
                if i > 0 and args.gapSeq:
                    writer.add(args.gapSeq)
                    pos += len(args.gapSeq)
                start = pos
                end = pos + len(seq) - 1
                mapF.write("\t".join(map(str, [
                    chrName, prefix, genomeName, cid, len(seq), start, end
                ]))+"\n")
                writer.add(seq)
                pos = end + 1
            writer.end()
            totals[prefix] = totals.get(prefix, 0) + (pos-1)

sys.stderr.write(prog+": wrote "+args.output+" and "+args.map+"\n")
for prefix in totals:
    sys.stderr.write("  "+prefix+"\t"+str(counters.get(prefix, 0) or chrIndex-1)
                     + " chromosomes\t"+format(totals[prefix], ",")+" bp\n")
