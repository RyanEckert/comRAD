#!/usr/bin/env python3
################################################################################
#Script Name    : ngsAdmixLauncher.py
#Description    : Print list of ngsAdmix commands for specified K value range
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import os, argparse, sys, math

parser = argparse.ArgumentParser(description="Prints a list of NGSadmix commands for a range of K (for comrad-jobs or launcher_creator.py)",
formatter_class = argparse.RawDescriptionHelpFormatter,
epilog = "examples:\n"
"  ngsAdmixLauncher.py -f ../angsd/apal.beagle.gz --maxK 7 -r 20 -n apal\n"
"      beagle already MAF-filtered by ANGSD; NGSadmix keeps every site (-minMaf 0)\n"
"  ngsAdmixLauncher.py -f ../angsd/apal.beagle.gz --maxK 7 -r 20 -n apal --minMaf 0.05\n"
"      apply NGSadmix's own MAF filter (its built-in default is 0.05)\n"
"  ngsAdmixLauncher.py -f ../angsd/apal.beagle.gz --maxK 7 -r 20 -n apal -P 4\n"
"  comrad-jobs apalNgsAdmix -n apalNgsAdmix --cpus 4 --mem 8GB\n"
"      4 NGSadmix threads per run; ask the scheduler for the same number of CPUs\n"
"\n"
"threads:\n"
"  -P sets NGSadmix's own thread count. The scheduler's CPU request only reserves\n"
"  cores; NGSadmix uses what -P says.\n")

parser.add_argument("-f", type = str, help = "Input beagle file name", action = "store")
parser.add_argument("--minK", type = int, help = "Minium value of K", action = "store", default='1')
parser.add_argument("--maxK", type = int, help = "Maximum value of K", action = "store")
parser.add_argument("-n", type = str, help = "Common name of outputs", action = "store", default='ngsAdmix')
parser.add_argument("-P", "--threads", type = int, help = "NGSadmix threads per run (NGSadmix -P) [default: 10]", action = "store", default=10)
parser.add_argument("--minMaf", type = float, help = "NGSadmix -minMaf (default 0: keep all beagle sites so .fopt matches the beagle for evalAdmix/NGSremix)", action = "store", default=0)
parser.add_argument("--launcher", help = "If used, create .slurm file with launcher_creator.py", action = "store_true")
parser.add_argument("-e", help = "If used, email to send SLURM info to", action = "store")
parser.add_argument("-r", type = int, help = "Number of replicate simulations to run per K", action = "store", default='10')

args = parser.parse_args()
minK = args.minK
maxK = args.maxK
f = args.f
e = args.e
n = args.n
r = args.r
minMaf = args.minMaf
launcher = args.launcher
filename = n+'NgsAdmix'

outfile = open(filename, 'w')

for Run in range(1, (r+1)):
    for K in range(minK, (maxK+1)):
        outfile.write('NGSadmix -likes '+f+' -K '+str(K)+' -minMaf '+str(minMaf)+' -P '+str(args.threads)+' -o '+n+'_k'+str(K)+'_run'+str(Run)+'\n')

outfile.close()

if launcher:
    os.system("module load launcher/3.9");
    if e:
        os.system('launcher_creator.py -j '+filename+' -n '+filename+' -q shortq7 -t 6:00:00 -e '+e)
    else:
        os.system('launcher_creator.py -j '+filename+' -n '+filename+' -q shortq7 -t 6:00:00')
