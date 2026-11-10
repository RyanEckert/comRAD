#!/usr/bin/env python3
################################################################################
#Script Name    : cloneTestRun.py
#Description    : Submit the clone-step ANGSD run in several variants (one change
#                 each) so their clone separation can be compared
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, os, shlex, subprocess, sys

prog = os.path.basename(sys.argv[0])
BASE = "-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05"
VARIANTS = {
    "base": "the current [clones] options (control)",
    "hetbias": "+ -hetbias_pval 1e-5 (drop sites with skewed read ratios at heterozygous calls)",
    "maxdepth": "+ -setMaxDepth at --depth-mult x the median total depth (collapsed duplicates collect about twice the reads); the median from --maxdepth-from, else <sp>Depth.depthGlobal; or --maxdepth",
    "post2": "-doPost 2 (flat genotype prior) instead of 1 (allele-frequency prior)",
    "mind50": "-minInd at 50%% of samples instead of --min-ind",
    "maf01": "-minMaf 0.01 instead of 0.05",
    "dep3": "+ -setMinDepthInd 3 (a sample's data at a site counts only with 3 or more reads)",
    "dep5": "+ -setMinDepthInd 5",
}

ap = argparse.ArgumentParser(
    description="Submits one ANGSD job per variant of the clone step, each in its own folder under --outdir, "
                "with the IBS matrix and genotype probabilities that cloneTestScore.py compares. Run it from the "
                "species' angsd folder. Nothing in that folder is changed.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""variants:
""" + "\n".join(f"  {k:<9} {v}" for k, v in VARIANTS.items()) + f"""

examples:
  cd $koko/2bRAD/nwgom/past/angsd
  {prog} -s past -b clonesIBSnofilter/pastCloneBams --dry-run
      print the commands without submitting
  {prog} -s past -b clonesIBSnofilter/pastCloneBams
      submit every variant but maxdepth to shortq7 (8 CPUs, 100G each); outputs in cloneTests/<variant>/
  {prog} -s ofav -b clones/ofavCloneBams -v base post2 -p longq7
      two variants only, on another partition
  {prog} -s past -b clonesIBSnofilter/pastCloneBams -v maxdepth --maxdepth-from cloneTests/base/past.pos.gz
      when the depth pass's histogram is capped (deeply sequenced samples): the limit from base's site depths
  cloneTestScore.py -s past -r "PGM016.1 PGM016.2 PGM016.3 PGM076.1 PGM076.2 PGM076.3 PGM119.1 PGM119.2 PGM119.3"
      once they finish: score them
""")
ap.add_argument("-s", "--sp", required=True, help="species prefix")
ap.add_argument("-b", "--bams", required=True, help="bam list for the clone step (e.g. clones/<sp>CloneBams); paths as ANGSD reads them from here")
ap.add_argument("--depth-global", help="ANGSD .depthGlobal from the depth pass (default <sp>Depth.depthGlobal here)")
ap.add_argument("--depth-mult", type=float, default=2.0, help="maxdepth variant: limit = this x the median total depth (default 2)")
ap.add_argument("--maxdepth-from", metavar="POS.GZ", help="maxdepth variant: take the median total depth from this ANGSD .pos.gz (e.g. cloneTests/base/<sp>.pos.gz once base has run)")
ap.add_argument("--maxdepth", type=int, help="maxdepth variant: use this limit directly")
ap.add_argument("-v", "--variants", nargs="+", default=[v for v in VARIANTS if v != "maxdepth"], choices=list(VARIANTS),
                help="variants to run (default all but maxdepth, which needs a depth source: see --maxdepth-from)")
ap.add_argument("-o", "--outdir", default="cloneTests", help="folder for the runs (default cloneTests)")
ap.add_argument("--min-ind", type=float, default=0.75, help="-minInd as a proportion of samples, for every variant but mind50 (default 0.75)")
ap.add_argument("--args", default=BASE, help="the base ANGSD filters (default: comRAD's [clones] defaults)")
ap.add_argument("-p", "--partition", default="shortq7", help="SLURM partition (default shortq7)")
ap.add_argument("-c", "--cpus", type=int, default=8, help="ANGSD threads and SLURM CPUs per job (default 8)")
ap.add_argument("--mem", default="100G", help="memory per job (default 100G)")
ap.add_argument("--dry-run", action="store_true", help="print the commands; submit nothing")
a = ap.parse_args()

nb = sum(1 for l in open(a.bams) if l.strip())
mind = lambda p: max(1, int(nb * p))


def depth_quantiles(path):
    """median and 99th percentile of the depthGlobal histogram (counts per total depth 0..max; ANGSD puts every
    site at or above the depth pass's -maxDepth into the last bin), depth 0 excluded"""
    v = [float(x) for x in open(path).read().split()]
    tot = sum(v[1:])
    if not tot:
        sys.exit(f"{prog}: {path} has no sites with depth > 0")
    out, c = {}, 0
    for d, k in enumerate(v[1:], 1):
        c += k
        for q in (0.5, 0.99):
            if q not in out and c >= q * tot:
                out[q] = d
    return out[0.5], out[0.99], len(v) - 1, v[-1] / tot


dg = a.depth_global or f"{a.sp}Depth.depthGlobal"
maxd = None
if "maxdepth" in a.variants:
    if a.maxdepth:
        maxd = a.maxdepth
        print(f"maxdepth: limit {maxd} (--maxdepth)")
    elif a.maxdepth_from:
        import gzip, statistics
        dep = []
        with gzip.open(a.maxdepth_from, "rt") as fh:
            for l in fh:
                f = l.split()
                if len(f) >= 3 and f[2].isdigit():
                    dep.append(int(f[2]))
        if not dep:
            sys.exit(f"{prog}: no site depths in {a.maxdepth_from}")
        med = statistics.median(dep)
        maxd = int(round(a.depth_mult * med))
        print(f"maxdepth: median total depth {med:g} over {len(dep)} sites in {a.maxdepth_from}; limit {a.depth_mult:g} x median = {maxd}")
    else:
        if not os.path.exists(dg):
            sys.exit(f"{prog}: maxdepth needs {dg}, --maxdepth-from or --maxdepth; or leave maxdepth out")
        med, p99, cap, top = depth_quantiles(dg)
        if med >= cap:
            sys.exit(f"{prog}: {top:.0%} of sites in {dg} are at or above the depth pass ceiling ({cap}), so its median is unknown. "
                     f"Run the other variants first, then: {prog} -s {a.sp} -b {a.bams} -v maxdepth --maxdepth-from {a.outdir}/base/{a.sp}.pos.gz")
        maxd = int(round(a.depth_mult * med))
        print(f"maxdepth: median total depth {med}, limit {a.depth_mult:g} x median = {maxd}"
              f" (99th percentile {p99}{' = the depth pass ceiling, so not usable' if p99 >= cap else ''};"
              f" {top:.1%} of sites at or above the ceiling of {cap})")

out = ["-GL 1", "-doMajorMinor 1", "-doMaf 1", "-doCounts 1", "-dumpCounts 1", "-makeMatrix 1", "-doIBS 1", "-doGeno 8"]
for v in a.variants:
    filt, post, mi = a.args, 1, mind(a.min_ind)
    if v == "hetbias":
        filt += " -hetbias_pval 1e-5"
    elif v == "maxdepth":
        filt += f" -setMaxDepth {maxd}"
    elif v == "post2":
        post = 2
    elif v == "mind50":
        mi = mind(0.5)
    elif v in ("dep3", "dep5"):
        filt += f" -setMinDepthInd {v[3:]}"
    elif v == "maf01":
        filt = filt.replace("-minMaf 0.05", "-minMaf 0.01") if "-minMaf" in filt else filt + " -minMaf 0.01"
    d = os.path.join(a.outdir, v)
    os.makedirs(d, exist_ok=True)
    cmd = (f"angsd -b {a.bams} {' '.join(out)} -doPost {post} {filt} -minInd {mi} -P {a.cpus} "
           f"-out {d}/{a.sp}")
    open(os.path.join(d, "command.txt"), "w").write(cmd + "\n")
    sb = ["sbatch", "-p", a.partition, "-c", str(a.cpus), f"--mem={a.mem}", "-J", f"{a.sp}CloneTest_{v}",
          "-o", f"{d}/{a.sp}.%j.log", "--wrap", cmd]
    if a.dry_run:
        print(" ".join(shlex.quote(x) for x in sb))
    else:
        r = subprocess.run(sb, capture_output=True, text=True)
        print(f"{v:<9} {r.stdout.strip() or r.stderr.strip()}")
