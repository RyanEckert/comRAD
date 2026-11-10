if inskip; then mkdir -p "$WD/relatedness" && : >"$WD/relatedness/${SP}KinDrop"; skip "in SKIP: no kin dropped (empty relatedness/${SP}KinDrop)"; fi
mkdir -p "$WD/relatedness" && cd "$WD/relatedness"
K=$(getK)
BEST=$(awk -F'\t' -v k=$K 'NR>1 && $1==k {print $NF}' ../ngsAdmix/${SP}_bestQ.txt)
need ../ngsAdmix/$BEST ../ngsAdmix/${BEST%.qopt}.fopt.gz
job "test -s ${SP}Remix" "NGSremix -beagle ../angsd/${SP}Unlinked.beagle.gz -qname ../ngsAdmix/$BEST -fname ../ngsAdmix/${BEST%.qopt}.fopt.gz -P $CFG_kin_cpus -o ${SP}Remix" >${SP}RemixJobs
comrad-jobs ${SP}RemixJobs -n ${SP}Remix --res kin || stop "NGSremix failed"
np=$(( $(wc -l <${SP}Remix) - 1 )); nb=$(wc -l <../angsd/${SP}Bams)
equal $np $(( nb * (nb - 1) / 2 )) "${SP}Remix pairs vs n(n-1)/2"

# pairs at theta >= KIN lose one member at a time: most remaining pairs first,
# then lower proportion of sites at >= QDEPTH x
python3 - ${SP}Remix ../angsd/${SP}Bams ../angsd/${SP}Quality${QDEPTH}X.txt $KIN ${SP} ../angsd/${SP}BamsClusters ${KIN_DROP:-yes} <<'PY' || stop "kin drop failed; see above"
import os, sys
remix, bams, qual, kin, sp, clusters, dodrop = sys.argv[1:8]
dodrop = dodrop != "no"
kin = float(kin)
key = lambda s: os.path.basename(s).removesuffix(".bam")
b = [l.strip() for l in open(bams) if l.strip()]
q = {key(f[0]): float(f[1]) for f in (l.split() for l in open(qual)) if len(f) >= 2}
lin = {key(f[0]): f[1] for f in (l.split() for l in open(clusters)) if len(f) >= 2}
miss = [x for x in b if key(x) not in q]
if miss:
    sys.exit(f"{len(miss)} bams not in {qual}, first {miss[0]}")
with open(remix) as fh:
    head = fh.readline().split()
    need = ["ind1", "ind2", "k0", "k1", "k2"]
    if not all(c in head for c in need):
        sys.exit(f"{remix} header is {head}; expected {need}")
    ix = {c: head.index(c) for c in need}
    rows = [l.split() for l in fh if l.strip()]
i1 = [int(r[ix["ind1"]]) for r in rows]
i2 = [int(r[ix["ind2"]]) for r in rows]
base = min(i1 + i2)
if base not in (0, 1) or max(i1 + i2) - base != len(b) - 1:
    sys.exit(f"{remix} indices run {base}..{max(i1 + i2)}; {bams} has {len(b)} lines")
print(f"NGSremix indices are {base}-based")
pairs = []
for r, x, y in zip(rows, i1, i2):
    th = float(r[ix["k1"]]) / 4 + float(r[ix["k2"]]) / 2
    if th >= kin:
        pairs.append((th, b[x - base], b[y - base], r[ix["k0"]], r[ix["k1"]], r[ix["k2"]]))
left, drop = list(pairs), []
while left:
    deg = {}
    for th, x, y, *_ in left:
        deg[x] = deg.get(x, 0) + 1
        deg[y] = deg.get(y, 0) + 1
    out = max(deg, key=lambda z: (deg[z], -q[key(z)]))
    drop.append(out)
    left = [r for r in left if out not in r[1:3]]
would = len(drop)
if not dodrop:
    drop = []
open(f"{sp}KinDrop", "w").write("".join(f"{x}\n" for x in drop))
if dodrop:
    if os.path.exists(f"{sp}KinKept"):
        os.remove(f"{sp}KinKept")
else:
    open(f"{sp}KinKept", "w").write(f"KIN_DROP = no: {len(pairs)} pairs at theta >= {kin} kept; {would} individuals would have been dropped\n")
with open(f"{sp}Kin", "w") as fh:
    fh.write("ind1\tind2\tlineage1\tlineage2\ttheta\tk0\tk1\tk2\tq5_1\tq5_2\tdropped\n")
    for th, x, y, k0, k1, k2 in sorted(pairs, reverse=True):
        d = ",".join(key(z) for z in (x, y) if z in drop)
        fh.write(f"{key(x)}\t{key(y)}\t{lin.get(key(x), '?')}\t{lin.get(key(y), '?')}\t{th:.4f}\t{k0}\t{k1}\t{k2}\t{q[key(x)]:.4f}\t{q[key(y)]:.4f}\t{d}\n")
cross = [p for p in pairs if lin.get(key(p[1])) != lin.get(key(p[2])) and "Admixed" not in (lin.get(key(p[1])), lin.get(key(p[2])))]
print(f"NOTE: {len(pairs)} pairs at theta >= {kin}; {len(drop)} of {len(b)} individuals dropped" if dodrop else
      f"NOTE: {len(pairs)} pairs at theta >= {kin}; KIN_DROP = no, so all {len(b)} kept ({would} would have been dropped)")
if cross:
    print(f"NOTE: {len(cross)} of those pairs are between two different lineages; check {sp}Kin before trusting them")
PY
cat ${SP}Kin

plots kin
