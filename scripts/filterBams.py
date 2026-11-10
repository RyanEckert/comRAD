#!/usr/bin/env python3
"""Filter a bam list by per-sample quality and/or clone groups."""
import argparse
import os

p = argparse.ArgumentParser(description="Drop low-quality bams and/or prune clone groups to their best member.")
p.add_argument("-b", "--bams",    required=True, help="bam list, one file/path per line")
p.add_argument("-q", "--quality", required=True, help="two columns: sample, value")
p.add_argument("-p", "--prop",    type=float, help="minimum value to keep; -p 0.25 -> >=0.25")
p.add_argument("-c", "--clones",  help="clone groups file, one group per line")
p.add_argument("-o", "--out",     help="output list (default: <bams>.q)")
a = p.parse_args()

if a.prop is None and not a.clones:
    p.error("nothing to do: give -p, -c, or both")

out = a.out or f"{a.bams}.q"

def key(s):
    """Sample name for a bam, whatever path or suffix it arrives with.

    The bamlist and the -q file carry whatever path ANGSD was given
    (../mappedReads/MGM001.bam); a clone groups file written from the
    dendrogram carries bare sample names. Stripping the directory is what
    lets the three key on the same thing."""
    s = os.path.basename(s)
    for suf in (".trim.bt2.bam", ".bam"):
        if s.endswith(suf):
            return s[:-len(suf)]
    return s

q = {}
for line in open(a.quality):
    f = line.split()
    if len(f) >= 2:
        q[key(f[0])] = float(f[1])

lines   = [l.strip() for l in open(a.bams) if l.strip()]
present = {key(b) for b in lines}
drops   = {}

if a.clones:
    for line in open(a.clones):
        m = [key(s) for s in line.split("#")[0].replace(",", " ").split()]
        m = [s for s in m if s in present]
        if len(m) < 2:
            continue
        best = max(m, key=lambda s: q.get(s, float("-inf")))
        for s in m:
            if s != best:
                drops[s] = f"clone of {best}"
        lost = ", ".join(f"{s} ({q.get(s, 0):.4f})" for s in m if s != best)
        print(f"clone: keep {best} ({q.get(best, 0):.4f})  drop {lost}")

keep, dropped, missing = [], [], []
for b in lines:
    k = key(b)
    if k in drops:                      dropped.append((b, drops[k]))
    elif a.prop is None:                keep.append(b)
    elif k not in q:                    missing.append(b); keep.append(b)
    elif q[k] < a.prop:                 dropped.append((b, f"{q[k]:.4f}"))
    else:                               keep.append(b)

with open(out, "w") as fh:
    fh.write("".join(b + "\n" for b in keep))

print(f"\n{out}: kept {len(keep)} of {len(lines)}, dropped {len(dropped)}")
for b, why in dropped:
    print(f"  {why:>18}  {b}")
if missing:
    print(f"WARNING: {len(missing)} bams absent from {a.quality} -- KEPT")
    for b in missing[:10]:
        print(f"  {b}")
