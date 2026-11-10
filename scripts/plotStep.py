#!/usr/bin/env python3
################################################################################
#Script Name    : plotStep.py
#Description    : Quick-look plots for each comRAD step (reads, symbionts, PCA,
#                 admixture, kinship, outliers, neutral run, FST, heterozygosity,
#                 theta), drawn from the files the step wrote
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, csv, glob, gzip, os, re, sys

import numpy as np

prog = os.path.basename(sys.argv[0])
KINDS = {
    "reads": "filteredReads/<SP>Reads.pdf: reads per sample; host (and symbiont) alignment rate per sample",
    "pools": "RawCounts.pdf (reads pipeline, project folder): raw and kept reads per pool",
    "symbionts": "bams/<SP>Symbionts.pdf: share of symbiont reads per genus, per sample",
    "pca": "pcangsd/<SP>Pca.pdf: PC1-2 and PC3-4 of the PCAngsd covariance, % variance on the axes",
    "admix": "ngsAdmix/<SP>Admix.pdf: NGSadmix bars for K = 2..MAXK (best run, one sample order, colours matched "
             "across K), and PCAngsd's own admixture at its K; evalAdmix residual RMSD and log-likelihood against K",
    "kin": "relatedness/<SP>Kin.pdf: pairwise kinship (theta) from NGSremix, the KIN cut marked",
    "outliers": "angsd/lineages/<SP>Outliers.pdf: PCAngsd selection statistic along the genome per lineage "
                "(outliers in red); within-lineage PCAs",
    "neutral": "ngsAdmix/<SP>Neutral.pdf: neutral, kin-free bars next to the full run's at K; neutral PCA",
    "fst": "sfs/<SP>Fst.pdf: weighted FST between lineages (neutral and all sites) and between groups within lineages",
    "het": "het/<SP>Het.pdf: heterozygosity per sample by lineage and group; against coverage",
    "theta": "theta/<SP>Theta.pdf: pi, Watterson's theta and Tajima's D per lineage (and the THETA_N subsample)",
}
ap = argparse.ArgumentParser(
    description="Draws comRAD's quick-look plots for one step from the files that step wrote, as one PDF per step "
                "(several panels or pages where there is more than one thing to show), plus PNGs of each page in "
                "run/<pipeline>/plots/ for comrad report. Meant for a first look; final figures belong in R.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="kinds:\n" + "\n".join(f"  {k:<10} {v}" for k, v in KINDS.items()) + f"""

examples (from the species folder):
  {prog} pca -s sint --meta meta.csv --meta-groups bank depthZone
      PCA coloured by bank and shaped by depthZone
  {prog} pca -s mcav --meta meta.csv --meta-groups bank depthZone --color source
      coloured by another metadata column (PCA_COLOR in popgen.config), e.g. published vs new reads
  {prog} admix -s ofav --meta meta.csv --meta-groups bank depthZone -k 2
      bars for every K, samples ordered by group, then by their main cluster at K = 2
  {prog} kin -s sint --kin 0.125
  {prog} pools --dir /path/to/project
      the reads pipeline's per-pool yield

sample order and colours:
  Samples are grouped by the --order factors (default --meta-groups, in that order), with a divider and label
  per group, and ordered within each group by main cluster at -k (default 2), then by its share. Colours of
  clusters are matched across K (each K's clusters to the previous K's, by overlap), so a lineage keeps its colour.
  A plot that already exists is kept: a rerun writes <name>_v2.pdf, then _v3, and so on (the report shows the newest).
  --kcolor colours the admixture clusters from one colour toward white in MAXK steps (PLOT_KCOLOR). --levels, --colors and --shapes set level order, colours and shapes for named factors; levels not listed
  come after the listed ones, with default colours and shapes. In popgen.config these are PLOT_COLOR,
  PLOT_SHAPE, PLOT_ORDER, PLOT_LEVELS, PLOT_COLORS and PLOT_SHAPES.

  {prog} pca -s ofav --meta meta.csv --meta-groups bank depthZone \\
      --levels "bank: WFGB EFGB Bright McGrail" --colors "depthZone: Shallow=gold Meso=1f78b4" --color depthZone --shape bank
""")
ap.add_argument("kind", choices=list(KINDS), help="which step's plots")
ap.add_argument("-s", "--sp", default="", help="species prefix (SP)")
ap.add_argument("--dir", default=".", help="the pipeline's folder (default here)")
ap.add_argument("--meta", help="sample metadata CSV")
ap.add_argument("--meta-sample", default="sample", help="metadata column with the sample names (default sample)")
ap.add_argument("--meta-groups", nargs="*", default=[], help="metadata columns: sampling groups (colour, then shape)")
ap.add_argument("--color", help="factor (metadata column, or lineage) that colours points (default: first --meta-groups)")
ap.add_argument("--shape", help="factor that sets point shapes (default: second --meta-groups)")
ap.add_argument("--order", nargs="*", help="factors that group and order samples along bar plots (default: --meta-groups)")
ap.add_argument("--levels", default="", help="level order per factor: 'bank: WFGB EFGB Stetson; depthZone: Shallow Meso'")
ap.add_argument("--colors", default="", help="colours per level: 'bank: WFGB=1b9e77 EFGB=d95f02; lineage: Mc1=red' (hex without #, or names)")
ap.add_argument("--shapes", default="", help="point shapes per level: 'depthZone: Shallow=o Meso=^' (matplotlib markers)")
ap.add_argument("--kcolor", help="admixture bars: this colour for the first cluster, the rest lighter toward white in MAXK steps (hex without #, or a name)")
ap.add_argument("--prefix", default="", help="admix: draw the run whose files start <SP><prefix> (Neutral: the kin-free, neutral K run, "
                                               "written to ngsAdmix/<SP>NeutralAdmix.pdf)")
ap.add_argument("-k", type=int, help="admix: K that orders the samples (default 2); neutral: the chosen K")
ap.add_argument("--kin", type=float, default=0.125, help="kin: the KIN cut (default 0.125)")
ap.add_argument("--qdepth", default="5", help="het: QDEPTH, for the Quality<N>X.txt coverage file (default 5)")
ap.add_argument("--pipeline", default="popgen", help="pipeline name, for run/<pipeline>/plots (default popgen)")
a = ap.parse_args()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
W = os.path.abspath(a.dir)
SP = a.sp
P = lambda *x: os.path.join(W, *x)
base = lambda s: re.sub(r"\.bam$", "", os.path.basename(s.strip().split()[0]))


def lines(path):
    return [l.rstrip("\n") for l in open(path) if l.strip()]


def table(path):
    with open(path) as fh:
        r = list(csv.reader(fh, delimiter="\t"))
    return [dict(zip(r[0], x)) for x in r[1:] if x]


def names(path):
    return [base(l) for l in lines(path)]


PAL = plt.get_cmap("tab10").colors + plt.get_cmap("Dark2").colors + plt.get_cmap("Set1").colors
MARK = ["o", "^", "s", "D", "v", "P", "X", "*", "<", ">"]
meta = {}
if a.meta and os.path.exists(a.meta):
    with open(a.meta, newline="", encoding="utf-8-sig") as fh:
        meta = {r.get(a.meta_sample, ""): r for r in csv.DictReader(fh)}
lineage = {}
if SP and os.path.exists(P("angsd", f"{SP}BamsClusters")):
    lineage = {base(l): l.split()[1] for l in lines(P("angsd", f"{SP}BamsClusters")) if len(l.split()) > 1}


def mval(s, col):
    if not col:
        return "all"
    if col == "lineage":
        return lineage.get(s, "NA")
    return (meta.get(s, {}).get(col) or "NA").strip() or "NA"


def spec(text, pairs):
    """'f: a b; g: c d' -> {f: [a, b]}; with pairs, 'f: a=x b=y' -> {f: {a: x}}"""
    out = {}
    for part in text.split(";"):
        if ":" not in part:
            continue
        f, rest = part.split(":", 1)
        items = rest.split()
        out[f.strip()] = dict(i.split("=", 1) for i in items if "=" in i) if pairs else items
    return out


LEVELS = spec(a.levels, False)
COLORS = {f: {k: ("#" + v if re.fullmatch(r"[0-9A-Fa-f]{6}", v) else v) for k, v in d.items()}
          for f, d in spec(a.colors, True).items()}
SHAPES = spec(a.shapes, True)
CCOL = a.color or (a.meta_groups[0] if a.meta_groups else None)
SCOL = a.shape or (a.meta_groups[1] if len(a.meta_groups) > 1 else None)
ORDER = a.order if a.order else a.meta_groups


def lsort(f, vals):
    """levels of factor f in PLOT_LEVELS order, then the rest alphabetically"""
    want = LEVELS.get(f, [])
    vals = set(vals)
    return [v for v in want if v in vals] + sorted(vals - set(want))


def rank(f, v):
    """sort key of level v of factor f: PLOT_LEVELS order first, then alphabetical"""
    want = LEVELS.get(f, [])
    return (want.index(v), "") if v in want else (len(want), v)


def gkey(s):
    return tuple(rank(f, mval(s, f)) for f in ORDER)


def group(s):
    return " ".join(mval(s, c) for c in ORDER) if ORDER else ""


def colours(f, vals):
    lv = lsort(f, vals)
    d = COLORS.get(f, {})
    free = iter(c for c in PAL if c not in d.values())
    return {v: d.get(v) or next(free, "grey") for v in lv}


def markers(f, vals):
    lv = lsort(f, vals)
    d = SHAPES.get(f, {})
    free = iter(m for m in MARK if m not in d.values())
    return {v: d.get(v) or next(free, "o") for v in lv}


pages = []


def save(fig, path):
    """add a page to the PDF being written"""
    pages.append((fig, path))


def flush():
    if not pages:
        print(f"{prog}: nothing to plot for {a.kind}")
        return
    pngdir = os.path.join(W, "run", a.pipeline, "plots")
    os.makedirs(pngdir, exist_ok=True)
    tag = a.kind + a.prefix.lower()
    for old in glob.glob(os.path.join(pngdir, f"{tag}-*.png")):
        os.remove(old)
    out = {}
    for fig, path in pages:
        out.setdefault(path, []).append(fig)
    n = 0
    for path, figs in out.items():
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path):                      # earlier plots are kept: <name>_v2.pdf, _v3, ...
            stem, v = path[:-4], 2
            while os.path.exists(f"{stem}_v{v}.pdf"):
                v += 1
            path = f"{stem}_v{v}.pdf"
        with PdfPages(path) as pdf:
            for fig in figs:
                pdf.savefig(fig)
                n += 1
                fig.savefig(os.path.join(pngdir, f"{tag}-{n}.png"), dpi=70)
                plt.close(fig)
        print(f"wrote {os.path.relpath(path, W)} ({len(figs)} page{'s' if len(figs) > 1 else ''})")


def legend(ax, cols, shapes, ctitle, stitle):
    h = [Line2D([], [], marker="o", ls="", color=c, label=v) for v, c in cols.items()]
    h += [Line2D([], [], marker=m, ls="", color="grey", label=v) for v, m in shapes.items()] if len(shapes) > 1 else []
    if h:
        ax.legend(handles=h, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.01, .5),
                  title=" / ".join(t for t in (ctitle, stitle if len(shapes) > 1 else "") if t), title_fontsize=7)


# ---------------------------------------------------------------- pieces ---
def pca_axes(axs, cov, ids, title, ccol, scol):
    C = np.loadtxt(cov, ndmin=2)
    val, vec = np.linalg.eigh(C)
    o = np.argsort(val)[::-1]
    val, vec = val[o], vec[:, o]
    pct = 100 * val / val[val > 0].sum()
    cv = [mval(s, ccol) for s in ids]
    sv = [mval(s, scol) for s in ids]
    cols = colours(ccol, cv)
    shapes = markers(scol, sv)
    for ax, (x, y) in zip(axs, ((0, 1), (2, 3))):
        if y >= vec.shape[1]:
            ax.axis("off")
            continue
        for v in shapes:
            m = np.array([s == v for s in sv])
            ax.scatter(vec[m, x], vec[m, y], c=[cols[c] for c, k in zip(cv, m) if k], marker=shapes[v],
                       s=18, edgecolor="black", linewidth=.2, alpha=.85)
        ax.set_xlabel(f"PC{x + 1} ({pct[x]:.1f}%)")
        ax.set_ylabel(f"PC{y + 1} ({pct[y]:.1f}%)")
    axs[0].set_title(title, loc="left", fontsize=8)
    legend(axs[-1], cols, shapes, ccol or "", scol or "")


def align(Q, prev):
    """reorder Q's columns to best match prev's (Hungarian on overlap); unmatched columns go last"""
    from scipy.optimize import linear_sum_assignment
    if prev is None:
        return Q
    cost = -(prev.T @ Q)                                   # k-1 x k
    r, c = linear_sum_assignment(cost)
    order = [c[list(r).index(i)] for i in range(prev.shape[1])] + [j for j in range(Q.shape[1]) if j not in c]
    return Q[:, order]


def sample_order(ids, Q):
    """group, then main cluster at the ordering K, then that cluster's share (descending)"""
    main = Q.argmax(1)
    return sorted(range(len(ids)), key=lambda i: (gkey(ids[i]), main[i], -Q[i, main[i]]))


def kpal(n):
    """admixture cluster colours: --kcolor ramped toward white over MAXK steps (white itself never used), so a
    cluster keeps its colour in every K panel; else the palette"""
    if not a.kcolor:
        return [PAL[j % len(PAL)] for j in range(n)]
    from matplotlib.colors import to_rgb
    c = np.array(to_rgb("#" + a.kcolor if re.fullmatch(r"[0-9A-Fa-f]{6}", a.kcolor) else a.kcolor))
    m = max(n, KMAX())                                  # one ramp for every panel: MAXK steps before white
    return [tuple(c + (1 - c) * i / m) for i in range(n)]


def KMAX():
    f = P("ngsAdmix", f"{SP}_bestQ.txt")
    return max((int(r["K"]) for r in table(f)), default=1) if os.path.exists(f) else 1


def bars(ax, Q, order, ids, label, show_names, pal=None):
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    for j in range(Q.shape[1]):
        ax.bar(x, Q[order, j], bottom=bottom, width=1, color=(pal or kpal(Q.shape[1]))[j], linewidth=0)
        bottom += Q[order, j]
    ax.set_xlim(-.5, len(order) - .5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center")
    ax.spines["left"].set_visible(False)
    grp = [group(ids[i]) for i in order]
    for i in range(1, len(order)):
        if grp[i] != grp[i - 1]:
            ax.axvline(i - .5, color="white", lw=1.2)
    ax.set_xticks(x if show_names else [])
    if show_names:
        ax.set_xticklabels([ids[i] for i in order], rotation=90, fontsize=4)


def group_labels(ax, order, ids):
    grp = [group(ids[i]) for i in order]
    start = 0
    for i in range(1, len(grp) + 1):
        if i == len(grp) or grp[i] != grp[start]:
            ax.text((start + i - 1) / 2, 1.04, grp[start], ha="center", va="bottom", fontsize=6, rotation=0,
                    transform=ax.get_xaxis_transform())
            start = i


def qread(path):
    return np.loadtxt(path, ndmin=2)


# ----------------------------------------------------------------- kinds ---
def k_reads():
    rc = P("filteredReads", f"{SP}ReadCounts")
    rows = table(rc)
    ids = [re.sub(r"\.filtered\.fastq\.gz$", "", r["file"]) for r in rows]
    n = np.array([float(r["reads"]) for r in rows])
    rates = [(lab, f) for lab, f in (("host", f"{SP}HostAlignmentRates"), ("symbiont", f"{SP}SymbiontAlignmentRates"))
             if os.path.exists(P("filteredReads", f))]
    o = np.argsort(n)
    fig, axs = plt.subplots(1 + len(rates), 1, figsize=(max(8, len(ids) * .06), 3 + 2 * len(rates)), squeeze=False)
    ax = axs[0, 0]
    ax.bar(range(len(ids)), n[o] / 1e6, width=.9, color=PAL[0])
    ax.axhline(np.median(n) / 1e6, color="black", lw=.6, ls="--")
    ax.set_ylabel("reads (millions)")
    ax.set_title(f"{SP}: reads per sample after quality filtering (dashed: median {np.median(n) / 1e6:.2f} M)", loc="left")
    ax.set_xticks(range(len(ids)))
    ax.set_xticklabels([ids[i] for i in o], rotation=90, fontsize=4)
    for k, (lab, f) in enumerate(rates, 1):
        r = {x["sample"]: float(x["alignRate"]) for x in table(P("filteredReads", f))}
        ax = axs[k, 0]
        ax.bar(range(len(ids)), [r.get(ids[i], np.nan) for i in o], width=.9, color=PAL[k])
        ax.set_ylabel(f"{lab} alignment (%)")
        ax.set_ylim(0, 100)
        ax.set_xticks([])
    axs[-1, 0].set_xlabel("samples, in the order of the top panel")
    fig.tight_layout()
    save(fig, P("filteredReads", f"{SP}Reads.pdf"))


def k_pools():
    f = next((x for x in [P("RawCounts")] + sorted(glob.glob(P("*", "RawCounts"))) if os.path.exists(x)), None)
    if not f:
        sys.exit(f"{prog}: no RawCounts in {W} or a folder below it")
    rows = table(f)
    pools = [r["pool"] for r in rows]
    raw = np.array([float(r["rawReads"]) for r in rows]) / 1e6
    kept = np.array([float(r["kept"]) for r in rows]) / 1e6
    fig, ax = plt.subplots(figsize=(max(6, len(pools) * .3), 4))
    x = np.arange(len(pools))
    ax.bar(x, raw, width=.8, color="lightgrey", label="raw")
    ax.bar(x, kept, width=.8, color=PAL[0], label="kept")
    ax.set_xticks(x)
    ax.set_xticklabels(pools, rotation=90, fontsize=6)
    ax.set_ylabel("reads (millions)")
    ax.legend(frameon=False)
    ax.set_title("reads per pool: raw and kept after deduplication and demultiplexing", loc="left")
    fig.tight_layout()
    save(fig, os.path.join(os.path.dirname(f), "RawCounts.pdf"))


def k_symbionts():
    rows = table(P("bams", f"{SP}SymGenusReads"))
    ids = sorted({r["sample"] for r in rows}, key=lambda s: (gkey(s), s))
    gen = sorted({r["genus"] for r in rows})
    M = np.zeros((len(ids), len(gen)))
    for r in rows:
        M[ids.index(r["sample"]), gen.index(r["genus"])] += float(r["reads"])
    tot = M.sum(1)
    Q = M / np.where(tot > 0, tot, 1)[:, None]
    fig, axs = plt.subplots(2, 1, figsize=(max(8, len(ids) * .07), 5), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    order = list(range(len(ids)))
    bars(axs[0], Q, order, ids, "share", False, [PAL[j % len(PAL)] for j in range(Q.shape[1])])
    group_labels(axs[0], order, ids)
    axs[0].legend(handles=[Line2D([], [], marker="s", ls="", color=PAL[j % len(PAL)], label=g) for j, g in enumerate(gen)],
                  fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.01, .5))
    axs[0].set_title(f"{SP}: symbiont reads (MAPQ >= 20) by genus", loc="left")
    axs[1].bar(range(len(ids)), tot / 1e3, width=1, color="grey")
    axs[1].set_ylabel("reads (k)")
    axs[1].set_xticks(range(len(ids)))
    axs[1].set_xticklabels(ids, rotation=90, fontsize=4)
    fig.tight_layout()
    save(fig, P("bams", f"{SP}Symbionts.pdf"))


def k_pca():
    ids = lines(P("pcangsd", f"{SP}Pcangsd.keep"))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    pca_axes(axs, P("pcangsd", f"{SP}Pcangsd.cov"), ids, f"{SP}: PCAngsd, unlinked SNPs, all samples ({len(ids)})", CCOL, SCOL)
    fig.tight_layout()
    save(fig, P("pcangsd", f"{SP}Pca.pdf"))


def best_q(k, pre=""):
    f = P("ngsAdmix", f"{SP}{pre}_k{k}_best.qopt")
    if os.path.exists(f):
        return qread(f)
    rows = [r for r in table(P("ngsAdmix", f"{SP}{pre}_bestQ.txt")) if r["K"] == str(k)]
    return qread(P("ngsAdmix", rows[0]["best"])) if rows else None


def k_admix():
    pre = a.prefix
    ids = names(P("angsd", f"{SP}{pre}Bams"))
    ks = sorted(int(r["K"]) for r in table(P("ngsAdmix", f"{SP}{pre}_bestQ.txt")))
    Qs, prev = {}, None
    for k in ks:
        Q = best_q(k, pre)
        if Q is None or Q.shape[0] != len(ids):
            continue
        Q = align(Q, prev)
        Qs[k] = prev = Q
    show = [k for k in Qs if k >= 2] or list(Qs)
    if not show:
        sys.exit(f"{prog}: no best .qopt per K in ngsAdmix/")
    ko = a.k if a.k in Qs else (2 if 2 in Qs else show[0])
    order = sample_order(ids, Qs[ko])
    # PCAngsd's own admixture at its K (eigenvectors + 1), if the structure step wrote it
    pq = sorted(glob.glob(P("pcangsd", f"{SP}{pre}Pcangsd.admix.*.Q")))
    rows = [(f"K = {k}", Qs[k]) for k in show]
    if pq:
        kp = int(re.sub(r".*\.admix\.(\d+)\.Q$", r"\1", pq[-1]))
        Qp = qread(pq[-1])
        if Qp.shape[0] == len(ids):
            rows.append((f"PCAngsd\nK = {kp}", align(Qp, Qs[kp]) if kp in Qs and kp > 1 else Qp))
    fig, axs = plt.subplots(len(rows), 1, figsize=(max(8, len(ids) * .07), 1 + .9 * len(rows)), sharex=True, squeeze=False)
    for i, (lab, Q) in enumerate(rows):
        bars(axs[i, 0], Q, order, ids, lab, i == len(rows) - 1)
    group_labels(axs[0, 0], order, ids)
    axs[0, 0].set_title(f"{SP}{': kin-free, neutral SNPs' if pre else ''}: NGSadmix, best run per K; samples by {' x '.join(ORDER) or 'sample'}, "
                        f"then main cluster at K = {ko}", loc="left", pad=14)
    fig.tight_layout()
    save(fig, P("ngsAdmix", f"{SP}{pre}Admix.pdf"))
    fit = P("ngsAdmix", f"{SP}{pre}_fit.txt")
    bq = table(P("ngsAdmix", f"{SP}{pre}_bestQ.txt"))
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.4))
    if os.path.exists(fit):
        f = table(fit)
        x = [int(r["K"]) for r in f]
        axs[0].plot(x, [float(r["rmsd"]) for r in f], "o-", color=PAL[0])
        for r in f:
            if r.get("elbow") == "<-":
                axs[0].scatter(int(r["K"]), float(r["rmsd"]), s=80, facecolor="none", edgecolor="red")
        axs[0].set_xlabel("K")
        axs[0].set_ylabel("evalAdmix residual RMSD")
        axs[0].set_title("residual fit (red: elbow)", loc="left")
    lk = [(int(r["K"]), float(r["lik"])) for r in bq if r.get("lik") not in (None, "", "NA")]
    if lk:
        axs[1].plot(*zip(*lk), "o-", color=PAL[1])
        axs[1].set_xlabel("K")
        axs[1].set_ylabel("log-likelihood, best run")
        axs[1].set_title("NGSadmix likelihood", loc="left")
    bf = P("ngsAdmix", f"{SP}{pre}_bic.txt")
    if os.path.exists(bf):
        rows = [l.rstrip("\n").split("\t") for l in open(bf) if l.strip() and not l.startswith(("#", "K\t"))]
        axs[2].plot([int(r[0]) for r in rows], [float(r[1]) for r in rows], "o-", color=PAL[2])
        for r in rows:
            if r[-1] == "<-":
                axs[2].scatter(int(r[0]), float(r[1]), s=80, facecolor="none", edgecolor="red")
        axs[2].set_xlabel("K")
        axs[2].set_ylabel("BIC")
        axs[2].set_title("k-means on PCA scores (red: lowest BIC)", loc="left")
    else:
        axs[2].axis("off")
    fig.tight_layout()
    save(fig, P("ngsAdmix", f"{SP}{pre}Admix.pdf"))


def k_kin():
    rows = table(P("relatedness", f"{SP}Remix"))
    th = np.array([float(r["k1"]) / 4 + float(r["k2"]) / 2 for r in rows])
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(th, bins=60, color=PAL[0])
    ax.set_yscale("log")
    ax.axvline(a.kin, color="red", ls="--", lw=.8)
    ax.set_xlabel("kinship (theta = k1/4 + k2/2)")
    ax.set_ylabel("pairs (log scale)")
    ax.set_title(f"{SP}: NGSremix kinship, {len(th)} pairs; {int((th >= a.kin).sum())} at or above KIN = {a.kin:g} (red)",
                 loc="left")
    fig.tight_layout()
    save(fig, P("relatedness", f"{SP}Kin.pdf"))


def pops():
    return sorted({v for v in lineage.values() if v != "Admixed"})


def k_outliers():
    L = pops()
    fig, axs = plt.subplots(len(L), 1, figsize=(10, 1.8 * len(L) + .6), squeeze=False)
    for ax, pop in zip(axs[:, 0], L):
        f = P("angsd", "lineages", f"{SP}{pop}Selection.tsv")
        rows = table(f) if os.path.exists(f) else []
        out = set()
        of = P("angsd", "lineages", f"{SP}{pop}Outliers")
        if os.path.exists(of):
            out = {tuple(l.split()[:2]) for l in lines(of)}
        if not rows:
            ax.text(.5, .5, f"{pop}: no SNPs tested", transform=ax.transAxes, ha="center")
            ax.axis("off")
            continue
        chrs = list(dict.fromkeys(r["chr"] for r in rows))
        ci = {c: i for i, c in enumerate(chrs)}
        x = np.arange(len(rows))
        p = np.array([max(float(r["p1"]), 1e-300) for r in rows])
        y = -np.log10(p)
        col = [PAL[0] if ci[r["chr"]] % 2 else "grey" for r in rows]
        hit = np.array([(r["chr"], r["pos"]) in out for r in rows])
        ax.scatter(x[~hit], y[~hit], s=2, c=[c for c, h in zip(col, hit) if not h], linewidth=0)
        ax.scatter(x[hit], y[hit], s=6, c="red", linewidth=0)
        ax.set_ylabel("-log10 p (PC1)")
        ax.set_xticks([])
        ax.set_title(f"{pop}: {len(rows)} SNPs, {int(hit.sum())} outliers (red); "
                     f"{len(chrs)} chromosomes or contigs, in reference order", loc="left")
    axs[-1, 0].set_xlabel("SNPs along the reference")
    fig.tight_layout()
    save(fig, P("angsd", "lineages", f"{SP}Outliers.pdf"))
    have = [pop for pop in L if os.path.exists(P("pcangsd", f"{SP}{pop}Pcangsd.cov"))]
    if have:
        fig, axs = plt.subplots(len(have), 2, figsize=(10, 3.8 * len(have)), squeeze=False)
        for i, pop in enumerate(have):
            ids = names(P("angsd", f"{SP}{pop}Bams"))
            pca_axes(axs[i], P("pcangsd", f"{SP}{pop}Pcangsd.cov"), ids, f"{pop}: within-lineage PCA, outlier tags removed",
                     CCOL, SCOL)
        fig.tight_layout()
        save(fig, P("angsd", "lineages", f"{SP}Outliers.pdf"))


def k_neutral():
    nb = names(P("angsd", f"{SP}NeutralBams"))
    Qn = qread(P("ngsAdmix", f"{SP}Neutral_best.qopt"))
    ids = names(P("angsd", f"{SP}Bams"))
    k = a.k or Qn.shape[1]
    Qf = best_q(k)
    if Qf is None:
        sys.exit(f"{prog}: no best .qopt for K = {k}")
    pos = {s: i for i, s in enumerate(ids)}
    Qf = Qf[[pos[s] for s in nb]]
    Qn = align(Qn, Qf)
    order = sample_order(nb, Qf)
    fig, axs = plt.subplots(2, 1, figsize=(max(8, len(nb) * .07), 3), sharex=True)
    bars(axs[0], Qf, order, nb, "all SNPs", False)
    bars(axs[1], Qn, order, nb, "neutral,\nkin-free", True)
    group_labels(axs[0], order, nb)
    axs[0].set_title(f"{SP}: K = {k}, the same {len(nb)} samples (kin left out), colours matched", loc="left", pad=14)
    fig.tight_layout()
    save(fig, P("ngsAdmix", f"{SP}Neutral.pdf"))
    if os.path.exists(P("pcangsd", f"{SP}NeutralPcangsd.cov")):
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
        pca_axes(axs, P("pcangsd", f"{SP}NeutralPcangsd.cov"), nb, f"{SP}: PCA, neutral, kin-free; coloured by lineage",
                 "lineage", SCOL or CCOL)
        fig.tight_layout()
        save(fig, P("ngsAdmix", f"{SP}Neutral.pdf"))


def heat(ax, labs, val, title, vmax):
    n = len(labs)
    M = np.full((n, n), np.nan)
    for (x, y), v in val.items():
        i, j = labs.index(x), labs.index(y)
        M[i, j] = M[j, i] = v
    im = ax.imshow(M, cmap="YlOrRd", vmin=0, vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labs, rotation=90, fontsize=6)
    ax.set_yticklabels(labs, fontsize=6)
    for (i, j), v in np.ndenumerate(M):
        if v == v:
            ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=5)
    ax.set_title(title, loc="left", fontsize=7)
    return im


def k_fst():
    sets = []
    for lab, f in (("between lineages, neutral sites", P("sfs", f"{SP}LineageFst.out")),
                   ("between lineages, all sites", P("sfs", f"{SP}LineageFstAll.out"))):
        if os.path.exists(f):
            rows = table(f)
            v = {(r["pop1"], r["pop2"]): float(r["weightedFst"]) for r in rows}
            labs = lsort("lineage", {x for k in v for x in k})
            sets.append((lab, labs, v))
    g = P("angsd", "groups", f"{SP}GroupFst.out")
    if os.path.exists(g):
        rows = table(g)
        for L in sorted({r["lineage"] for r in rows}):
            v = {(r["group1"], r["group2"]): float(r["weightedFst"]) for r in rows if r["lineage"] == L}
            sets.append((f"{L}: between groups", sorted({x for k in v for x in k}), v))
    if not sets:
        return
    vmax = max(max(v.values()) for _, _, v in sets) or 1
    nc = min(3, len(sets))
    nr = -(-len(sets) // nc)
    fig, axs = plt.subplots(nr, nc, figsize=(3.6 * nc + 1, 3.4 * nr), squeeze=False)
    im = None
    for ax, (lab, labs, v) in zip(axs.flat, sets):
        im = heat(ax, labs, v, lab, vmax)
    for ax in list(axs.flat)[len(sets):]:
        ax.axis("off")
    fig.colorbar(im, ax=axs.ravel().tolist(), shrink=.6, label="weighted FST")
    fig.suptitle(f"{SP}: FST (realSFS, weighted), one colour scale", x=.01, ha="left", fontsize=8)
    save(fig, P("sfs", f"{SP}Fst.pdf"))


def k_het():
    rows = table(P("het", f"{SP}.het"))
    L = lsort("lineage", {r["lineage"] for r in rows})
    lc = colours("lineage", L)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), gridspec_kw={"width_ratios": [1, 2, 1.2]})
    rnd = np.random.default_rng(1)
    data = [[float(r["het"]) for r in rows if r["lineage"] == l] for l in L]
    axs[0].boxplot(data, showfliers=False)
    axs[0].set_xticks(range(1, len(L) + 1))
    axs[0].set_xticklabels(L)
    for i, d in enumerate(data, 1):
        axs[0].scatter(i + rnd.uniform(-.15, .15, len(d)), d, s=6, color=lc[L[i - 1]], alpha=.7)
    axs[0].set_ylabel("heterozygosity (per site)")
    axs[0].set_title(f"{SP}: per sample, by lineage", loc="left")
    first = {}
    for r in rows:
        first.setdefault((r["lineage"], group(r["sample"])), r["sample"])
    gl = sorted(first, key=lambda k: (rank("lineage", k[0]), gkey(first[k])))
    data = [[float(r["het"]) for r in rows if (r["lineage"], group(r["sample"])) == k] for k in gl]
    axs[1].boxplot(data, showfliers=False)
    for i, (k, d) in enumerate(zip(gl, data), 1):
        axs[1].scatter(i + rnd.uniform(-.15, .15, len(d)), d, s=6, color=lc[k[0]], alpha=.7)
    axs[1].set_xticks(range(1, len(gl) + 1))
    axs[1].set_xticklabels([f"{l} {g}".strip() for l, g in gl], rotation=90, fontsize=6)
    axs[1].set_title("by lineage and group", loc="left")
    q = P("angsd", f"{SP}Quality{a.qdepth}X.txt")
    if os.path.exists(q):
        cov = {base(l): float(l.split()[1]) for l in lines(q) if len(l.split()) > 1}
        for i, l in enumerate(L):
            pts = [(cov[r["sample"]], float(r["het"])) for r in rows if r["lineage"] == l and r["sample"] in cov]
            if pts:
                axs[2].scatter(*zip(*pts), s=8, color=lc[l], label=l)
        axs[2].set_xlabel(f"proportion of sites at >= {a.qdepth}x")
        axs[2].set_ylabel("heterozygosity")
        axs[2].legend(fontsize=6, frameon=False)
        axs[2].set_title("against coverage", loc="left")
    fig.tight_layout()
    save(fig, P("het", f"{SP}Het.pdf"))


def k_theta():
    rows = sorted(table(P("theta", f"{SP}Theta.out")), key=lambda r: rank("lineage", r["lineage"]))
    L = [r["lineage"] for r in rows]
    sub = "piSub" in rows[0]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    x = np.arange(len(L))
    w = .38 if sub else .7
    for ax, (col, lab) in zip(axs, (("pi", "pi (per site)"), ("thetaW", "Watterson's theta (per site)"), ("tajimaD", "Tajima's D"))):
        ax.bar(x - (w / 2 if sub else 0), [float(r[col]) for r in rows], width=w, color=PAL[0],
               label="all pure individuals")
        if sub:
            ax.bar(x + w / 2, [float(r[col + "Sub"]) for r in rows], width=w, color=PAL[1],
                   label=f"{rows[0]['nSub']} per lineage")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{l}\n(n = {r['n']})" for l, r in zip(L, rows)], fontsize=6)
        ax.set_title(lab, loc="left")
        if col == "tajimaD":
            ax.axhline(0, color="black", lw=.5)
    axs[-1].legend(fontsize=6, frameon=False)
    fig.suptitle(f"{SP}: diversity per lineage", x=.01, ha="left", fontsize=8)
    fig.tight_layout()
    save(fig, P("theta", f"{SP}Theta.pdf"))


globals()["k_" + a.kind]()
flush()
