#!/usr/bin/env python3
################################################################################
#Script Name : zipper.py
#Description : Print list of files to gzip or gunzip with launcher_creator.py
#Author       : Ryan J Eckert
#Email         : ryan.j.eckert@gmail.com
################################################################################
import os, argparse, sys, glob

prog = os.path.basename(sys.argv[0])

parser = argparse.ArgumentParser(description="Prints a list of gzip/gunzip commands to be used with launcher_creator.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='''examples:
  # gzip every .fastq in the current directory, maximum compression
  {0} -f fastq -a -9 --launcher -e you@fau.edu

  # several file types at once
  {0} -f fastq tr0 trim sam -a -9 --launcher -e you@fau.edu

  # gunzip everything ending in .gz
  {0} --gunzip --launcher -e you@fau.edu

  # write the commands file only, then submit it yourself
  {0} -f sam -a -9 -n zipSams
'''.format(prog))

parser.add_argument("-f", type=str, nargs='+', help="File type(s) to gzip, don't include '.' [e.g. fastq tr0 sam]", action = "store")
parser.add_argument("--gunzip", help="If used, gunzip all .gz files", action = "store_true")
parser.add_argument("-a", type=str, help="Optional arguments to pass to g/unzip", action = "store", default='')
parser.add_argument("-n", type=str, help="Name of output file to write launcher commands to", action = "store", default='zip')
parser.add_argument("--launcher", help="If used, create .slurm file with launcher_creator.py", action = "store_true")
parser.add_argument("-e", help="If used, email to send SLURM info to", action = "store")
parser.add_argument("-w", help="Wayness for launcher", action = "store", default = "20")


args = parser.parse_args()
gunzip = args.gunzip
f = args.f
a = args.a
e = args.e
n = args.n
w = args.w
launcher = args.launcher

if not gunzip and not f:
    sys.exit(prog+": -f is required unless --gunzip is used")

outfile = open(n, 'w')

if gunzip:
    for file in sorted(glob.iglob(r'*'+'gz')):
        outfile.write('gunzip '+file+'\n')
else:
    seen = set()
    for suffix in f:
        for file in sorted(glob.iglob(r'*'+suffix)):
            if file not in seen:
                seen.add(file)
                outfile.write('gzip '+a+' '+file+'\n')


outfile.close()

if launcher:
    os.system("module load launcher");
    if e:
        os.system('launcher_creator.py -j '+n+' -n '+n+' -w '+w+' -q shortq7 -t 6:00:00 -e '+e)
    else:
        os.system('launcher_creator.py -j '+n+' -n '+n+' -w '+w+' -q shortq7 -t 6:00:00')
