#!/bin/bash
# clone-step regression on simulated pedigrees (simCloneKin.py -> cloneGroups.py, judged against truth.tsv).
# Expected: ~7x: genotype distance chosen, every clone grouped, no relative merged (also with a large clone of a
# highly heterozygous genet); ~2.5x: IBS chosen, every clone grouped. Usage: bash sim_clones.sh [workdir]
H=$(cd "$(dirname "$0")/.." && pwd)
W=${1:-$(mktemp -d)}; mkdir -p $W && cd $W
fail=0
judge(){ python3 - "$1" "$2" "$3" <<'PY'
import sys
d, want_metric, strict = sys.argv[1], sys.argv[2], sys.argv[3] == "strict"
g = [l.split() for l in open(f"{d}/simCloneGroups")]
grp = {s: i for i, x in enumerate(g) for s in x}
same = lambda a, b: a in grp and grp.get(a) == grp.get(b)
cl = tot = 0; bad = {}
for l in list(open(f"{d}/truth.tsv"))[1:]:
    a, b, r = l.split()
    if r in ("clone", "replicate"):
        tot += 1; cl += same(a, b)
    elif same(a, b):
        bad[r] = bad.get(r, 0) + 1
metric = dict(l.split("\t") for l in open(f"{d}/simCloneMetric.txt").read().splitlines())["metric"]
ok = cl == tot and metric == want_metric and (not strict or not bad)
print(f"{'ok  ' if ok else 'FAIL'} {d}: {metric}, clones {cl}/{tot}, relatives merged {bad or 0}")
sys.exit(0 if ok else 1)
PY
}
for spec in "simA::geno:strict" "simP:--big-clone 25 --fixed-het 0.3:geno:strict" "simL:--cov 2.5:ibs:loose"; do
  IFS=: read -r d opts metric strict <<<"$spec"
  python3 $H/scripts/simCloneKin.py -o $d $opts >/dev/null || { echo "FAIL $d: simulation"; fail=1; continue; }
  (cd $d && python3 $H/scripts/cloneGroups.py -i cloneTests/base/sim.ibsMat -g cloneTests/base/sim.geno.gz -b simBams \
     -r "$(cat techreps.txt)" -o sim >/dev/null) || { echo "FAIL $d: cloneGroups.py"; fail=1; continue; }
  judge $d $metric $strict || fail=1
done
exit $fail
