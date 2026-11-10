#!/usr/bin/env python3
################################################################################
#Script Name    : cloneGroups.py
#Description    : Clone groups from ANGSD output: genotype distance or IBS, the
#                 cut placed from technical replicates, a replicate check and a
#                 dendrogram; or (--check) a clone-free set checked against a cut
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, csv, gzip, os, re, sys
from collections import Counter

import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="Builds an average-linkage (UPGMA) tree of the samples from one of two distances, cuts it into "
                "clone groups, checks the technical replicates and draws the tree. The distances: IBS (1 - IBS from "
                "ANGSD -doIBS 1: one random read per sample and site) and genotype distance (the share of sites "
                "where two samples' genotype calls differ, among sites both called with probability above a "
                "confidence level, from ANGSD -doGeno 8 -doPost 1). With --metric auto the genotype distance is "
                "used when the data support it (enough shared confident sites for most pairs, and replicates "
                "tight compared with typical distances), at 0.99 confidence, else 0.95; otherwise IBS.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  {prog} -i ofavClones.ibsMat -g ofavClones.geno.gz -b ofavCloneBams -r "OGM024.1 OGM024.2 OGM066.1 OGM066.2" -o ofav
      automatic choice of distance and cut; writes ofavCloneGroups, ofavCloneReplicates.txt,
      ofavCloneMetric.txt, ofavCloneBorderline.txt and ofavCloneDend.pdf (plus the other distance's tree)
  {prog} -i ofavClones.ibsMat -b ofavCloneBams -r "OGM024.1 OGM024.2" --metric ibs --buffer 1.10 -o ofav
      IBS only, cut at the replicates' deepest join x 1.10 (treeTrimmer)
  {prog} -i pastClones.ibsMat -g pastClones.geno.gz -b pastCloneBams -r "$TECHREPS" --cut 0.075 -o past
      a cut you chose, in the units of the distance used (here the automatic one)
  {prog} -i ofavClones.ibsMat -g ofavClones.geno.gz -b ofavCloneBams -r "OGM024.1 OGM024.2" -o ofav \\
      --meta metadata.csv --meta-sample sample --meta-groups bank depthZone --quality ofavQuality5X.txt
      tips coloured by bank and shaped by depthZone; coverage reported for borderline joins
  {prog} -i ofav.ibsMat -g ofav.geno.gz -b ofavBams --metric geno --conf 0.99 --check 0.0676 -o ofav
      a clone-free set (the genotype step's output): tree with the clone cut, and any pair below it

distance and cut (--metric, --cut):
  auto    genotype distance at the first --conf level that passes both checks
          (--min-well of pairs with at least --min-shared shared confident sites; replicate
          tightness, largest replicate distance / median distance, at most --max-tight), else IBS
  geno    genotype distance at the first --conf level that passes, else the last level, with a warning
  ibs     IBS only
  both    as geno, and the IBS tree and grouping are written too, with the pairs they disagree on
  cut auto: for genotype distance, the widest gap between joins whose lower edge is within --window x
  the replicates' largest distance; for IBS, treeTrimmer (the replicates' deepest join x --buffer).
  gap and replicates force one of those; a number sets the cut directly.

replicates:
  -r names technical replicates as <sample>.<n>; a set is every name with the same <sample>.
  Names with no row are left out with a warning; a set needs two or more to count.

exit status:
  0 done; 3 no replicate set has two or more samples and the cut is not a number;
  4 the groups would drop more than --maxdrop of the samples (reason in <out>CloneStop.txt); 1 otherwise.
""")
ap.add_argument("-i", "--ibs", required=True, help="ANGSD .ibsMat (1 - IBS), one row per bam")
ap.add_argument("-g", "--geno", help="ANGSD .geno.gz with genotype probabilities (-doGeno 8 or 15, -doPost 1); without it, IBS only")
ap.add_argument("-b", "--bams", required=True, help="the bam list ANGSD read, in the same order")
ap.add_argument("-r", "--techreps", default="", help="technical replicates, space-separated <sample>.<n> names")
ap.add_argument("--metric", default="auto", choices=["auto", "geno", "ibs", "both"], help="distance (default auto)")
ap.add_argument("--cut", default="auto", help="auto, gap, replicates or a number (default auto)")
ap.add_argument("--conf", type=float, nargs="+", default=[0.99, 0.95], help="genotype call confidence levels, tried in order (default 0.99 0.95)")
ap.add_argument("--min-shared", type=int, default=1000, help="shared confident sites for a pair to count as well measured (default 1000)")
ap.add_argument("--min-well", type=float, default=0.90, help="share of pairs that must be well measured (default 0.90)")
ap.add_argument("--max-tight", type=float, default=0.15, help="largest replicate distance / median distance allowed (default 0.15)")
ap.add_argument("--window", type=float, default=3.0, help="gap cut: gaps must start within this x the replicates' largest distance (default 3)")
ap.add_argument("--buffer", type=float, default=1.10, help="treeTrimmer cut: the replicates' deepest join x this (default 1.10)")
ap.add_argument("--maxdrop", type=float, help="exit 4 if the groups would drop more than this proportion of the samples")
ap.add_argument("--check", type=float, metavar="CUT", help="check mode: tree of a clone-free set with this cut; warn for pairs below it")
ap.add_argument("--quality", help="<SP>Quality<N>X.txt (bam, proportion of sites at >= N x), for coverage in the borderline report")
ap.add_argument("--meta", help="sample metadata CSV, for the dendrogram's tip colours and shapes")
ap.add_argument("--meta-sample", default="sample", help="metadata column with the sample names (default sample)")
ap.add_argument("--meta-groups", nargs="*", default=[], help="metadata columns: the first fills the tips, the second shapes them")
ap.add_argument("--color", help="metadata column that colours the tips (default: first --meta-groups)")
ap.add_argument("--shape", help="metadata column that shapes the tips (default: second --meta-groups)")
ap.add_argument("--levels", default="", help="level order per column, as PLOT_LEVELS: 'bank: WFGB EFGB; depthZone: Shallow Mesophotic'")
ap.add_argument("--colors", default="", help="colours per level, as PLOT_COLORS: 'bank: WFGB=D53E4F EFGB=F88D51' (hex without #, or names)")
ap.add_argument("--shapes", default="", help="tip shapes per level, as PLOT_SHAPES: 'depthZone: Shallow=^ Mesophotic=v'")
ap.add_argument("-o", "--out", required=True, help="output prefix, e.g. the species code")
ap.add_argument("--dir", default=".", help="folder for the outputs (default here)")
a = ap.parse_args()


def note(*x):
    print("NOTE:", *x, flush=True)


def f(x):
    return os.path.join(a.dir, a.out + x)


base = lambda p: re.sub(r"\.bam$", "", os.path.basename(p))
names = [base(l.split()[0]) for l in open(a.bams) if l.strip()]
n = len(names)
pos = {s: i for i, s in enumerate(names)}
iu = np.triu_indices(n, 1)
M = np.loadtxt(a.ibs, ndmin=2)
if M.shape != (n, n):
    sys.exit(f"{prog}: {a.ibs} is {M.shape[0]} x {M.shape[1]}; {a.bams} has {n} bams")
techReps = a.techreps.split()
cov = {}
if a.quality and os.path.exists(a.quality):
    for l in open(a.quality):
        x = l.split()
        if len(x) >= 2:
            cov[base(x[0])] = float(x[1])


# ---- genotype distance
def read_geno(path):
    """called genotype (0/1/2) and its probability, per site and sample. Reads -doGeno 8 (3 probabilities per
    sample) and -doGeno 15 (major, minor, then per sample: call, call as letters, 3 probabilities)"""
    G, P = [], []
    layout = None
    with gzip.open(path, "rt") as fh:
        for l in fh:
            t = l.split()[2:]
            if layout is None:
                if len(t) == 3 * n:
                    layout = "8"
                elif len(t) == 2 + 5 * n:
                    layout = "15"
                else:
                    raise ValueError(f"{path}: {len(t)} fields after chr and pos; expected {3 * n} or {2 + 5 * n} for {n} samples")
            if layout == "8":
                q = np.array(t, float).reshape(n, 3)
            else:
                q = np.array(t[2:], dtype=object).reshape(n, 5)[:, 2:].astype(float)
            if not G and not np.allclose(q.sum(1), 1, atol=0.02):
                raise ValueError(f"{path}: the genotype probabilities of the first site do not sum to 1 (layout -doGeno {layout} assumed)")
            G.append(q.argmax(1).astype(np.int8)); P.append(q.max(1).astype(np.float32))
    if not G:
        raise ValueError(f"{path}: no sites")
    G, P = np.array(G), np.array(P)
    return G, P


def genodist(G, P, c):
    """share of sites with different calls among sites both samples called with probability > c, and the
    number of such shared sites, for every pair (by matrix products)"""
    C = P > c
    Cf = C.astype(np.float32)
    S = Cf.T @ Cf
    same = np.zeros_like(S)
    for g in (0, 1, 2):
        X = (C & (G == g)).astype(np.float32)
        same += X.T @ X
    S = np.rint(S).astype(np.int64); same = np.rint(same).astype(np.int64)
    D = (S - same) / np.maximum(S, 1)
    D = (D + D.T) / 2; np.fill_diagonal(D, 0)
    return D, S


# ---- replicate sets
sets = {}
gone = [t for t in techReps if t not in pos]
for t in techReps:
    if t in pos:
        sets.setdefault(t.split(".")[0], []).append(pos[t])
short = sorted(k for k, v in sets.items() if len(v) < 2)
sets = {k: v for k, v in sets.items() if len(v) >= 2}


def repmax(D):
    return max((D[i, j] for v in sets.values() for i in v for j in v if i < j), default=float("nan"))


def tight(D):
    m = float(np.median(D[iu]))
    return repmax(D) / m if sets and m > 0 else float("nan")


# ---- dendrogram
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

meta = {}
if a.meta and (a.meta_groups or a.color or a.shape) and os.path.exists(a.meta):
    with open(a.meta, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            meta[row.get(a.meta_sample, "")] = row


def spec(text, pairs):
    """'f: a b; g: c d' -> {f: [a, b]}; with pairs, 'f: a=x b=y' -> {f: {a: x}}"""
    out = {}
    for part in text.split(";"):
        if ":" in part:
            f, rest = part.split(":", 1)
            out[f.strip()] = dict(i.split("=", 1) for i in rest.split() if "=" in i) if pairs else rest.split()
    return out


LEVELS, COLORS, SHAPES = spec(a.levels, False), spec(a.colors, True), spec(a.shapes, True)


def ordered(col, vals):
    want = LEVELS.get(col, [])
    return [v for v in want if v in vals] + sorted(set(vals) - set(want))


def draw(Z, thresh, ylabel, path, title=""):
    fillCol = a.color or (a.meta_groups[0] if len(a.meta_groups) >= 1 else None)
    shapeCol = a.shape or (a.meta_groups[1] if len(a.meta_groups) >= 2 else None)
    val = lambda s, c: ((meta.get(s, {}).get(c) or "NA").strip() or "NA") if c else "all"
    fills = ordered(fillCol, {val(s, fillCol) for s in names})
    shapes = ordered(shapeCol, {val(s, shapeCol) for s in names})
    cmap = plt.get_cmap("tab10" if len(fills) <= 10 else "tab20")
    given = {k: ("#" + v if re.fullmatch(r"[0-9A-Fa-f]{6}", v) else v) for k, v in COLORS.get(fillCol, {}).items()}
    free = iter(c for i, c in enumerate(cmap(j % cmap.N) for j in range(40)))
    colour = {v: given.get(v) or next(free) for v in fills}
    gm = SHAPES.get(shapeCol, {})
    freem = iter(m for m in ["o", "^", "v", "s", "D", "P", "X", "*"] * 5 if m not in gm.values())
    marker = {v: gm.get(v) or next(freem) for v in shapes}
    top = float(Z[:, 2].max()) if len(Z) else 1.0
    drop = 0.02 * top                                    # leaves hang a little below their join
    fig, ax = plt.subplots(figsize=(max(12, n * .25), 9))
    dn = dendrogram(Z, labels=names, no_plot=True)
    leafY = {}
    for xs, ys in zip(dn["icoord"], dn["dcoord"]):
        ys = list(ys)
        if ys[0] == 0:
            ys[0] = ys[1] - drop; leafY[xs[0]] = ys[0]
        if ys[3] == 0:
            ys[3] = ys[2] - drop; leafY[xs[3]] = ys[3]
        ax.plot(xs, ys, color="black", lw=.5)
    ymin = 0
    for i, s in enumerate(dn["ivl"]):
        x = 5 + 10 * i
        y = leafY.get(x, 0)
        ymin = min(ymin, y)
        ax.scatter(x, y, s=40, marker=marker[val(s, shapeCol)], facecolor=colour[val(s, fillCol)], edgecolor="black", linewidth=.25, zorder=3)
        ax.text(x, y - (1.5 if s in techReps else 1.0) * drop, s, rotation=90, ha="center", va="top", fontsize=7,
                fontweight="bold" if s in techReps else "normal")
    if thresh == thresh:
        ax.axhline(thresh, color="red", ls=(0, (5, 3)), lw=.75)
    ax.set_xlim(0, 10 * n + 10)
    ax.set_ylim(min(ymin, 0) - 6 * drop, top * 1.03)
    ax.set_ylabel(ylabel)
    ax.set_xticks([])
    for sname in ("top", "right", "bottom"):
        ax.spines[sname].set_visible(False)
    handles = [Line2D([], [], marker="s", ls="", markerfacecolor=colour[v], markeredgecolor="black", label=v) for v in fills] if fillCol else []
    handles += [Line2D([], [], marker=marker[v], ls="", markerfacecolor="white", markeredgecolor="black", label=v) for v in shapes] if shapeCol else []
    if handles:
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.02), ncol=min(len(handles), 12), frameon=False,
                  title=" / ".join(c for c in (fillCol, shapeCol) if c))
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


LABEL = {"ibs": "Genetic distance (1 - IBS)", "geno": "Genotype distance (proportion of sites with different genotypes)"}

# ---- the genotype distance, if any
geno = {}
if a.geno and a.metric != "ibs":
    try:
        G, P = read_geno(a.geno)
        for c in a.conf:
            D, S = genodist(G, P, c)
            well = float((S[iu] >= a.min_shared).mean())
            few = S < 50                                   # too few shared sites to measure: set to the largest distance
            np.fill_diagonal(few, False)
            if few.any():
                D = D.copy(); D[few] = D[~few].max() if (~few).any() else 1
            geno[c] = dict(D=D, S=S, well=well, tight=tight(D), sites=len(G), shared=int(np.median(S[iu])), few=int(few[iu].sum()))
    except (OSError, ValueError) as e:
        note(f"WARNING: genotype distance not available ({e}); using IBS")

# ---- check mode: a clone-free set against a cut
if a.check is not None:
    if a.metric == "geno" and geno:
        c = a.conf[0]; D = geno[c]["D"]; lab = LABEL["geno"] + f", calls > {c:g}"
    else:
        D = M; lab = LABEL["ibs"]
    Z = linkage(squareform(D, checks=False), "average")
    close = sorted((D[i, j], names[i], names[j]) for i, j in zip(*iu))
    below = [x for x in close if x[0] < a.check]
    draw(Z, a.check, lab, f("CloneFreeDend.pdf"), title=f"{a.out}: clone-free samples; red: the clone cut ({a.check:.4f})")
    fmtc = lambda s: f" ({cov[s]:.2f})" if s in cov else ""
    note(f"clone-free check: closest pairs {', '.join(f'{x}{fmtc(x)}-{y}{fmtc(y)} {d:.4f}' for d, x, y in close[:5])}")
    if below:
        note(f"WARNING: {len(below)} pairs of clone-free samples are below the clone cut ({a.check:.4f}): "
             + ", ".join(f"{x}-{y} {d:.4f}" for d, x, y in below[:10])
             + "; the sites differ from the clone run, so look at these before removing anything (" + a.out + "CloneFreeDend.pdf)")
    else:
        note(f"clone-free check: no pair below the clone cut ({a.check:.4f})")
    sys.exit(0)

# ---- choose the distance
if gone:
    note("WARNING:", len(gone), "of", len(techReps), "TECHREPS are not in the clone step (no reads, or below CLONE_MINQ) and are left out:", *gone)
if short:
    note("WARNING: replicate sets with fewer than two bams are not used for the cut:", *short)
ibs_t = tight(M)
passing = [c for c in a.conf if c in geno and geno[c]["well"] >= a.min_well and not (geno[c]["tight"] > a.max_tight)]
why = []
for c in a.conf:
    if c in geno:
        g = geno[c]
        why.append(f"genotypes > {c:g}: {g['well']:.0%} of pairs with >= {a.min_shared} shared sites (median {g['shared']}), "
                   f"replicate tightness {g['tight']:.3f}")
why.append(f"IBS: replicate tightness {ibs_t:.3f}")
if a.metric == "ibs" or not geno:
    metric, conf = "ibs", None
elif a.metric == "auto":
    metric, conf = ("geno", passing[0]) if passing else ("ibs", None)
else:
    metric, conf = "geno", (passing[0] if passing else a.conf[-1])
    if not passing:
        note(f"WARNING: CLONE_METRIC = {a.metric}, but the genotype distance does not pass its checks at any confidence level; "
             f"using {conf:g} anyway")
D = geno[conf]["D"] if metric == "geno" else M
if metric == "geno" and geno[conf]["few"]:
    note(f"WARNING: {geno[conf]['few']} pairs share fewer than 50 confident sites; they are set to the largest distance (not grouped)")
note(f"clone distance: {'genotype distance, calls > ' + format(conf, 'g') if metric == 'geno' else 'IBS'} ({a.metric}); " + "; ".join(why))
if a.metric == "auto" and metric == "geno" and ibs_t == ibs_t:
    note(f"IBS was not used: genotype distance passes its checks (IBS replicate tightness {ibs_t:.3f}, "
         f"genotypes {geno[conf]['tight']:.3f}); both trees are drawn")
elif a.metric == "auto" and metric == "ibs" and geno:
    note("genotype distance did not pass its checks (too few shared confident sites, or replicates not tight): low coverage; using IBS")

Z = linkage(squareform(D, checks=False), "average")
heights = np.unique(Z[:, 2])
rmax = repmax(D) if sets else float("nan")


def cut_tree(h, Zx=None):
    return fcluster(Z if Zx is None else Zx, t=h, criterion="distance")


def treetrimmer(Zx=None):
    """treeTrimmer::tree_trimmer(method = "prefix", separator = "."): each set's lowest tree height at which
    its members share one cluster; the deepest across sets x --buffer"""
    raw = []
    for k, v in sets.items():
        idx = [i for i, s in enumerate(names) if re.match("^" + re.escape(k) + re.escape("."), s)]
        hs = heights if Zx is None else np.unique(Zx[:, 2])
        h = next((h for h in hs if len(set(cut_tree(h, Zx)[idx])) == 1), None)
        if h is None:
            note("WARNING: replicate set", k, "never joins in the tree")
        else:
            raw.append(h)
    return max(raw) * a.buffer if raw else None


def windowed_gap():
    up = [rmax] + sorted(h for h in heights if h > rmax)
    win = [(up[k + 1] - up[k], up[k], up[k + 1]) for k in range(len(up) - 1) if up[k] <= a.window * rmax]
    return max(win) if win else None


ctype = a.cut
try:
    thresh = float(ctype); how = "set by CLONE_CUT"; gap = None
except ValueError:
    if ctype not in ("auto", "gap", "replicates"):
        sys.exit(f"{prog}: --cut must be auto, gap, replicates or a number, not {ctype}")
    if not sets:
        sys.exit(3)
    if ctype == "auto":
        ctype = "gap" if metric == "geno" else "replicates"
    gap = None
    if ctype == "gap":
        gap = windowed_gap()
        if gap:
            thresh = gap[1] + 1e-9
            how = f"bottom of the widest gap near the replicates ({gap[1]:.4f} to {gap[2]:.4f})"
        else:
            note("WARNING: no gap between joins within the window above the replicates; using treeTrimmer")
            ctype = "replicates"
    if ctype == "replicates":
        t = treetrimmer()
        if t is None:
            sys.exit(3)
        thresh = t; how = f"treeTrimmer: replicates' deepest join x {a.buffer:g}"
print("clone threshold", f"{thresh:.7g}", how)
note(f"clone threshold {thresh:.4f} ({how})")

# ---- clone groups, numbered and listed in sample order (as R's cutree and split)
cl = cut_tree(thresh)
order = list(dict.fromkeys(cl))
groups = [[names[i] for i in range(n) if cl[i] == c] for c in order]
groups = [g for g in groups if len(g) > 1]
with open(f("CloneGroups"), "w") as fh:
    fh.write("".join(" ".join(g) + "\n" for g in groups))
nDrop = sum(len(g) - 1 for g in groups)
note(len(groups), "clone groups;", nDrop, "of", n, "samples would be dropped as clones")

# ---- borderline joins: grouped well above replicate level, or the nearest join above the cut
kids = {}
for k, (x, y, h, _) in enumerate(Z):
    kids[n + k] = (int(x), int(y), h)


def members(node):
    if node < n:
        return [node]
    x, y, _ = kids[node]
    return members(x) + members(y)


def side(node):
    m = [names[i] for i in members(node)]
    cv = [cov[s] for s in m if s in cov]
    return ("+".join(m[:4]) + (f"+{len(m) - 4} more" if len(m) > 4 else "")
            + (f" (lowest coverage {min(cv):.2f})" if cv else ""))


border = []
if rmax == rmax:
    above = sorted((h, k) for k, (x, y, h) in kids.items() if h > thresh)
    for k, (x, y, h) in sorted(kids.items(), key=lambda kv: kv[1][2]):
        if 1.5 * rmax < h <= thresh:
            border.append(("grouped", h, side(x), side(y)))
    if above:
        h, k = above[0]
        x, y, _ = kids[k]
        border.append(("not grouped (nearest join above the cut)", h, side(x), side(y)))
with open(f("CloneBorderline.txt"), "w") as fh:
    fh.write("call\tjoin\tside1\tside2\n" + "".join(f"{c}\t{h:.4f}\t{s1}\t{s2}\n" for c, h, s1, s2 in border))
for c, h, s1, s2 in border:
    if c == "grouped":
        note(f"borderline: {s1} joins {s2} at {h:.4f}, grouped as clones but more than 1.5 x the replicates' distance "
             f"({rmax:.4f}); a weak sample on one side points to a clone, good coverage on both to a close relative")

# ---- replicate check: nearest sample is one of its own replicates, or at least in its clone group
reps = [t for t in techReps if t in pos]
if reps:
    setOf = {r: r.split(".")[0] for r in reps}
    rows = []
    for r in reps:
        i = pos[r]
        others = [j for j in range(n) if j != i]
        d = D[i, others]
        k = int(np.argmin(d)); nn = others[k]
        own = [j for j in others if names[j] in setOf and setOf[names[j]] == setOf[r]]
        oth = [j for j in others if j not in own]
        ko = min(oth, key=lambda j: D[i, j]) if oth else None
        rows.append(dict(replicate=r, set=setOf[r], nearest=names[nn], nearestDist=round(float(D[i, nn]), 4),
                         ownMax=round(float(max(D[i, j] for j in own)), 4) if own else None,
                         nearestOther=names[ko] if ko is not None else "NA",
                         nearestOtherDist=round(float(D[i, ko]), 4) if ko is not None else None,
                         ok=nn in own, sameGroup=bool(cl[i] == cl[nn])))
    cols = ["replicate", "set", "nearest", "nearestDist", "ownMax", "nearestOther", "nearestOtherDist", "ok", "sameGroup"]
    fmt = lambda v: "NA" if v is None else ("TRUE" if v is True else "FALSE" if v is False else f"{v:g}" if isinstance(v, float) else str(v))
    with open(f("CloneReplicates.txt"), "w") as fh:
        fh.write("\t".join(cols) + "\n" + "".join("\t".join(fmt(r[c]) for c in cols) + "\n" for r in rows))
    print("\t".join(cols))
    for r in rows:
        print("\t".join(fmt(r[c]) for c in cols))
    gaps = [r["nearestOtherDist"] - r["ownMax"] for r in rows if r["ownMax"] is not None and r["nearestOtherDist"] is not None]
    note(f"replicates: {sum(r['ok'] for r in rows)} of {len(rows)} have one of their own replicates as nearest sample, "
         f"{sum(r['ok'] or r['sameGroup'] for r in rows)} a sample in their own clone group"
         + (f"; smallest gap to a sample outside their set {min(gaps):.4f}" if gaps else ""))
    bad = [f"{r['replicate']} ({r['nearest']})" for r in rows if not r["ok"] and not r["sameGroup"]]
    if bad:
        note("WARNING: replicates whose nearest sample is outside their clone group:", *bad, f"; see {a.out}CloneReplicates.txt")

# ---- the record of the choice
with open(f("CloneMetric.txt"), "w") as fh:
    fh.write(f"metric\t{metric}\nconfidence\t{conf if conf else 'NA'}\nrequested\t{a.metric}\ncut\t{thresh:.6f}\ncutRule\t{how}\n"
             f"replicateMax\t{rmax:.6f}\nibsTightness\t{ibs_t:.4f}\n"
             + "".join(f"geno{c:g}WellMeasured\t{g['well']:.3f}\ngeno{c:g}Tightness\t{g['tight']:.4f}\ngeno{c:g}MedianShared\t{g['shared']}\n"
                       for c, g in geno.items())
             + f"groups\t{len(groups)}\ndropped\t{nDrop}\n")

# ---- trees: the one used, and the other distance when both were built
title = (f"{a.out}: {'genotype distance, calls > ' + format(conf, 'g') if metric == 'geno' else '1 - IBS'}; "
         f"red: clone cut {thresh:.4f} ({len(groups)} groups, {nDrop} dropped)")
draw(Z, thresh, LABEL[metric] + (f", calls > {conf:g}" if metric == "geno" else ""), f("CloneDend.pdf"), title)
if geno and (metric == "geno" or a.metric == "both"):
    if metric == "geno":
        Zi = linkage(squareform(M, checks=False), "average")
        draw(Zi, float("nan"), LABEL["ibs"], f("CloneDendIBS.pdf"), f"{a.out}: 1 - IBS, for comparison (not used for the clone groups)")
    if a.metric == "both" and metric == "geno":
        Zi = linkage(squareform(M, checks=False), "average")
        ti = treetrimmer(Zi)
        if ti:
            ci = fcluster(Zi, t=ti, criterion="distance")
            diff = [(names[i], names[j]) for i, j in zip(*iu) if (cl[i] == cl[j]) != (ci[i] == ci[j])]
            note(f"CLONE_METRIC = both: IBS with treeTrimmer ({ti:.4f}) groups {len(diff)} pairs differently"
                 + (": " + ", ".join(f"{x}-{y}" for x, y in diff[:10]) if diff else ""))

if a.maxdrop is not None and nDrop > a.maxdrop * n:
    open(f("CloneStop.txt"), "w").write(
        f"clone grouping at {thresh:.4f} would drop {nDrop} of {n} samples (more than CLONE_MAXDROP = {a.maxdrop:g}); "
        f"look at {a.out}CloneDend.pdf, {a.out}CloneReplicates.txt and {a.out}CloneBorderline.txt, then set CLONE_CUT\n")
    sys.exit(4)
