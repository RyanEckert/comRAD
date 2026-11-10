#!/usr/bin/env python3
################################################################################
#Script Name    : localSaf.py
#Description    : Runs realSFS or thetaStat, failing on any bgzf read error;
#                 optionally on copies of their bgzipped inputs on node-local disk
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, glob, os, shutil, subprocess, sys, tempfile

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="Runs a realSFS or thetaStat command and exits 75 if it printed any bgzf read error, so the "
                "job is retried instead of passing a partial result on (realSFS can write partial output with "
                "exit 0 after one). With local copies on (--local-copy yes, or [cluster] local_copy = yes, "
                "which comrad exports as CFG_cluster_local_copy), every .saf.idx, .thetas.idx and .fst.idx it "
                "names is first copied, with the .gz files beside it, to node-local disk ($TMPDIR, else /tmp), "
                "and the command reads the copies: on some parallel file systems (BeeGFS) bgzf reads fail "
                "intermittently on intact files while plain copies read fine. Anything the command writes next "
                "to a copied input (thetaStat's .pestPG) is copied back beside the original.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  {prog} realSFS ../het/SGM184.saf.idx -fold 1 >../het/SGM184.ml
      one sample's SFS, read in place (local copies off, the default); fails on a bgzf error
  {prog} --local-copy yes realSFS ../het/SGM184.saf.idx -fold 1 >../het/SGM184.ml
      the same, read from a copy on node-local disk
  {prog} realSFS ofavCluster1SFS.saf.idx ofavCluster2SFS.saf.idx -maxIter 1000 -P 4 >ofavCluster1.Cluster2.sfs
      a 2D SFS
  {prog} thetaStat do_stat ../theta/sintCluster1.thetas.idx
      writes ../theta/sintCluster1.thetas.idx.pestPG as thetaStat would
  {prog} --local-copy yes --tmp /scratch/me realSFS fst stats ofavCluster1.Cluster2.fst.idx
      copies to /scratch/me instead of $TMPDIR

Output paths the command names (>, -outname, -fstout) are left alone and written where they point.
""")
ap.add_argument("--local-copy", choices=["yes", "no"],
                help="read copies on node-local disk (default: $CFG_cluster_local_copy, else no)")
ap.add_argument("--tmp", help="local folder for the copies (default $TMPDIR, else /tmp)")
ap.add_argument("cmd", nargs=argparse.REMAINDER, help="the realSFS or thetaStat command and its arguments")
a = ap.parse_args()
cmd = a.cmd[1:] if a.cmd[:1] == ["--"] else a.cmd
if not cmd:
    ap.error("give a command, e.g. realSFS x.saf.idx")

# index -> the bgzipped files read with it
local = (a.local_copy or os.environ.get("CFG_cluster_local_copy", "").strip() or "no") == "yes"
KINDS = {".saf.idx": (".saf.gz", ".saf.pos.gz"), ".thetas.idx": (".thetas.gz",), ".fst.idx": (".fst.gz",)}
tmp = tempfile.mkdtemp(prefix="localSaf.", dir=a.tmp or os.environ.get("TMPDIR") or None) if local else None
staged = {}                                   # original index -> its local copy
try:
    out = []
    for x in cmd:
        kind = next((k for k in KINDS if x.endswith(k)), None)
        if local and kind and os.path.exists(x) and x not in staged:
            sub = os.path.join(tmp, str(len(staged)))
            os.makedirs(sub)
            for f in [x] + [x[:-len(kind)] + e for e in KINDS[kind]]:
                if os.path.exists(f):
                    shutil.copy2(f, sub)
            staged[x] = os.path.join(sub, os.path.basename(x))
        out.append(staged.get(x, x))
    before = {os.path.dirname(loc): (os.path.dirname(x) or ".", set(os.listdir(os.path.dirname(loc))))
              for x, loc in staged.items()}
    p = subprocess.Popen(out, stderr=subprocess.PIPE, text=True, errors="replace")
    bad, missing = 0, ""
    for line in p.stderr:
        sys.stderr.write(line)
        if "bgzf" in line or line.startswith("[E::"):
            bad += 1
        elif "Problem opening file" in line:          # realSFS exits 0 after this
            missing = line.strip()
    rc = p.wait()
    for d, (orig, had) in before.items():           # e.g. thetaStat's .pestPG, written beside its input
        for f in set(os.listdir(d)) - had:
            shutil.copy2(os.path.join(d, f), os.path.join(orig, f))
    if missing:
        print(f"{prog}: {cmd[0]} could not open an input ({missing.lstrip('-> ')}); failing", file=sys.stderr)
        sys.exit(1)
    if bad:
        print(f"{prog}: {bad} bgzf read error(s) from {cmd[0]}" + (" even on local copies" if local else
              "; if the files pass gzip -t, set local_copy = yes in cluster.config") + "; failing so the job is retried",
              file=sys.stderr)
        sys.exit(75)
    sys.exit(rc)
finally:
    if tmp:
        shutil.rmtree(tmp, ignore_errors=True)
