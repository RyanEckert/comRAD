#!/usr/bin/env python3
################################################################################
#Script Name	: bowtie2Launcher.py
#Description	: Write a list of bowtie2 commands, one per read file
#Note         	: default settings are tuned for 36 bp 2bRAD tags
#Author       	: Ryan J Eckert
#Email         	: ryan.j.eckert@gmail.com

#Modified from '2bRAD_bowtie2_launch.pl' by Misha Matz
#(https://github.com/z0on/2bRAD_denovo)
################################################################################
import os, argparse, sys, glob, shlex

prog = os.path.basename(sys.argv[0])
ALIGN = '--no-unal --score-min L,16,1 --local -L 16'

parser = argparse.ArgumentParser(description='Prints a list of bowtie2 commands, one per read file (for comrad-jobs or launcher_creator.py). Default alignment settings (' + ALIGN + ') are tuned for 36 bp 2bRAD tags.',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='''examples:
  # map every .trim in the working directory
  {0} -g $GENOME -f .trim -n maps --logdir logs

  # write BAM instead of SAM: same size as a gzip -9 SAM, but usable as-is
  {0} -g $GENOME -f .trim -n maps --bam --logdir logs

  # gzipped reads work too; bowtie2 decompresses them itself
  {0} -g $GENOME -f .trim.gz -n maps --bam --logdir logs

  # 4 bowtie2 threads per sample (ask the scheduler for 4 CPUs per task as well)
  {0} -g $GENOME -f .trim.gz -n maps --bam -p 4
  comrad-jobs maps -n maps --cpus 4 --mem 8GB

  # comRAD names: <sample>.host.bam and logs/<sample>.host.bowtie2.log, where the
  # sample is the file name without -f
  {0} -g $GENOME -f .filtered.fastq.gz -n hostMaps -s host --plain

  # different alignment settings (replaces the 2bRAD defaults entirely)
  {0} -g $GENOME -f .trim.gz -n maps --bam --args '--no-unal --local --very-sensitive-local'

  # two independent mappings of the SAME reads, for splitHostSym.py.
  # -s tags each pass so the two don't overwrite each other; --samdir is
  # optional on top of that.
  {0} -g $SYMGENOME  -f .trim -n symMaps  -s sym --bam
  {0} -g $HOSTGENOME -f .trim -n hostMaps -s host --bam

  # pull the alignment rates out of the logs afterwards
  grep -H 'overall alignment rate' logs/*.bt2log | sed 's|logs/||; s|\\\\.bt2log:| |; s|% overall alignment rate||'

  # rerun only the jobs whose log has no alignment rate (TAG is the -s value,
  # or '' if -s was not used)
  awk -v tag=.sym '{{if(!match($0, /-U [^ ]+/)) next;
        s=substr($0, RSTART+3, RLENGTH-3); sub(/\\\\.gz$/, "", s);
        f="logs/" s tag ".bt2log"; ok=0;
        while((getline l < f)>0) if(l ~ /overall alignment rate/) ok=1;
        close(f); if(!ok) print}}' symMaps > symMapsRerun

index:
  -g is the bowtie2 index basename, which for these genomes is the .fasta
  path itself: bowtie2-build writes <basename>.1.bt2 .. .rev.2.bt2 beside it
  (.bt2l for a large index). The job list is not written unless those exist,
  since bowtie2 would otherwise exit before emitting a single SAM line and
  every task would fail with samtools reporting 'fail to read the header'.

threads:
  -p sets bowtie2's own thread count. The scheduler's CPU request only reserves
  cores; bowtie2 uses one thread unless -p says otherwise.

exit status:
  With --bam the command is run through 'bash -o pipefail' so that a bowtie2
  that dies mid-run fails the task. Without it the pipeline reports samtools'
  status, and samtools exits 0 on a truncated stream -- leaving a valid but
  incomplete BAM, which splitHostSym.py cannot distinguish from reads that
  genuinely did not map.
'''.format(prog))
parser.add_argument("-f", type=str, help="File type to map, don't include [e.g. trim]", action = "store", default='trim')
parser.add_argument("-g", type=str, help="Path to reference genome", action = "store", default='')
parser.add_argument("-s", type=str, help="Optional tag inserted into the SAM and log names, e.g. -s sym gives SAMPLE.trim.sym.bt2.sam and SAMPLE.trim.sym.bt2log. Use a different -s and/or --samdir per reference so two mappings of the same reads don't overwrite each other", action = "store", default='')
parser.add_argument("-p", "--threads", type=int, help="bowtie2 threads per sample (bowtie2 -p) [default: 1]", action = "store", default=1)
parser.add_argument("--args", type=str, help="bowtie2 alignment options, replacing the default '" + ALIGN + "'", action = "store", default=ALIGN)
parser.add_argument("--plain", help="Name outputs <sample>.<tag>.bam (or .sam) and logs <sample>.<tag>.bowtie2.log, where <sample> is the file name without -f; needs -s", action = "store_true")
parser.add_argument("--samdir", type=str, help="Optional directory for alignment output [default: cwd]", action = "store", default='.')
parser.add_argument("--bam", help="If used, pipe bowtie2 through samtools and write .bam instead of .sam. Same size as a gzip -9 SAM but needs no decompressing, and splitHostSym.py reads it directly", action = "store_true")
parser.add_argument("-l", type=str, help="BAM compression level for --bam [default: 9]", action = "store", default='9')
parser.add_argument("--logdir", type=str, help="Directory for per-sample bowtie2 stderr logs [default: logs], '' to disable", action = "store", default='logs')
parser.add_argument("-n", type=str, help="Name of output file to write the commands to", action = "store", default='maps')
parser.add_argument("--launcher", help="If used, create .slurm file with launcher_creator.py", action = "store_true")
parser.add_argument("-N", type=str, help="With --launcher: the number of nodes to split job across in slurm", action = "store")
parser.add_argument("-w", type=str, help="With --launcher: wayness, the number of commands to run per node", action = "store")
parser.add_argument("-t", type=str, help="With --launcher: walltime [default: 6:00:00]", action = "store", default='6:00:00')
parser.add_argument("-e", help="With --launcher: email to send SLURM info to", action = "store")
parser.add_argument("-q", help="With --launcher: partition to submit job to", action = "store",default='shortq7')

args = parser.parse_args()
launcher = args.launcher
e = args.e
f = args.f
g = args.g
n = args.n
ld = args.logdir
sd = args.samdir
s = ('.'+args.s.lstrip('.')) if args.s else ''
bam = args.bam
lvl = args.l
N = args.N
w = args.w
t = args.t
q = args.q
threads = ' -p '+str(args.threads) if args.threads > 1 else ''
if args.plain and not args.s:
    sys.exit(f"{prog}: --plain needs -s (e.g. -s host), so the output is <sample>.host.bam")

# every task needs the index, so a missing one is checked once here rather
# than found by each job separately
if not g:
    sys.exit(f"{prog}: -g is required: the bowtie2 index basename")
built = os.path.expanduser(g)
if not (os.path.exists(built + '.1.bt2') or os.path.exists(built + '.1.bt2l')):
    sys.exit(f"{prog}: no bowtie2 index for {g} -- expected {g}.1.bt2 or "
             f"{g}.1.bt2l.\n  Build it with: bowtie2-build {g} {g}")

for d in [ld, sd]:
    if d and d != '.':
        os.makedirs(d, exist_ok=True)

outfile = open(n, 'w')

for file in sorted(glob.iglob(r'*'+f)):
    # a trailing .gz is dropped from the output and log names, so they are the
    # same whether or not the reads are gzipped
    base = file[:-3] if file.endswith('.gz') else file
    if args.plain:
        sample = file[:-len(f)] if file.endswith(f) else base
        out = (sd+'/' if sd and sd != '.' else '')+sample+s+'.'+('bam' if bam else 'sam')
        logname = sample+s+'.bowtie2.log'
    else:
        out = (sd+'/' if sd and sd != '.' else '')+base+s+'.bt2.'+('bam' if bam else 'sam')
        logname = base+s+'.bt2log'
    list = 'bowtie2 '+args.args.strip()+threads+' -x '+g+' -U '+file
    if not bam:
        list = list+' -S '+out
    if ld:
        list = list+' 2> '+ld+'/'+logname
    if bam:
        # samtools is last in the pipe and exits 0 on a truncated stream, so
        # without pipefail a dead bowtie2 leaves a valid-looking partial BAM
        list = ('bash -o pipefail -c '
                + shlex.quote(list+' | samtools view -b --output-fmt-option level='
                              +lvl+' -o '+out))
    outfile.write(list+'\n')

outfile.close()

if launcher:
    os.system("module load launcher");
    cmd = 'launcher_creator.py -j '+n+' -n '+n
    if N:
        cmd = cmd+' -N '+N
    cmd = cmd+' -q '+q+' -t '+t
    if w:
        cmd = cmd+' -w '+w
    if e:
        cmd = cmd+' -e '+e
    os.system(cmd)
