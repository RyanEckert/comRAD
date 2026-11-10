#!/usr/bin/env python3
################################################################################
#Script Name    : momentsModels.py
#Description    : Fit two-population divergence models to a realSFS spectrum
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Fit two-population divergence models to a folded realSFS spectrum with moments.

Base histories
  SI    split, no gene flow                       nu1 nu2 T
  IM    split, constant symmetric gene flow        nu1 nu2 T m
  IMa   split, constant asymmetric gene flow       nu1 nu2 T m12 m21
  AM    gene flow for T1, then isolation for T2    nu1 nu2 T1 T2 m
  SC    isolation for T1, then gene flow for T2    nu1 nu2 T1 T2 m

Genomic heterogeneity (suffixes)
  2m    a fraction P of the genome exchanges at m, the rest at me
        (barrier loci / islands of differentiation). IM, AM, SC only.
  2N    a fraction Q of the genome has Ne scaled by hrf in both daughter
        populations (linked selection). All five bases.
  2m2N  both, the four components weighted P*(1-Q), P*Q, (1-P)*(1-Q), (1-P)*Q.

Ghost populations (only with --ghost, or by name with -m)
  IMghost  IM, plus an unsampled population G that split from the ancestor dT
           before the pop1/pop2 split and sends migrants into both at mG.
           nu1 nu2 nuG T dT m mG
  SCghost  SC with the same ghost.  nu1 nu2 nuG T1 T2 dT m mG

Units are moments' coalescent units: sizes relative to the ancestral Ne, times
in 2*Nref generations, migration as 2*Nref*m (m = fraction of the receiving
population replaced per generation). m12 is into pop1 from pop2. With --mu the
table also carries Nref = theta / (4 * mu * L), sizes in individuals, times in
generations and migration per generation; L is every site in the spectrum,
monomorphic included.

With --boot (realSFS -bootstrap output, one spectrum per line) the best run of
each model also gets, from the Godambe information matrix (Coffman et al. 2016,
Mol Biol Evol 33:591):
  CLAIC     composite-likelihood AIC, -2*ll + 2*tr(J H^-1) (Varin & Vidoni 2005,
            Biometrika 92:519). Plain AIC assumes independent sites and
            overstates support for complex models when SNPs are linked; CLAIC
            replaces the parameter count with the effective one.
  SEs       per-parameter standard errors, including theta, in <out>.godambe.tsv.
"""
import argparse
import functools
import os
import sys

import numpy as np

prog = os.path.basename(sys.argv[0])


# --------------------------------------------------------------- building ---
@functools.lru_cache(maxsize=None)
def _split_cached(n0, n1):
    import moments
    sts = moments.LinearSystem_1D.steady_state_1D(n0 + n1)
    return moments.Manips.split_1D_to_2D(moments.Spectrum(sts), n0, n1)


def _integrate(ns, nu1, nu2, epochs):
    """Equilibrium split, then one integration per (T, m12, m21) epoch."""
    fs = _split_cached(int(ns[0]), int(ns[1])).copy()
    for T, m12, m21 in epochs:
        fs.integrate([nu1, nu2], T, m=np.array([[0.0, m12], [m21, 0.0]]))
    return fs


# epochs as a function of the base parameters and the migration rate to use
BASE = {
    "SI":  (["nu1", "nu2", "T"],
            lambda d, m: [(d["T"], 0, 0)]),
    "IM":  (["nu1", "nu2", "T", "m"],
            lambda d, m: [(d["T"], m, m)]),
    "IMa": (["nu1", "nu2", "T", "m12", "m21"],
            lambda d, m: [(d["T"], d["m12"], d["m21"])]),
    "AM":  (["nu1", "nu2", "T1", "T2", "m"],
            lambda d, m: [(d["T1"], m, m), (d["T2"], 0, 0)]),
    "SC":  (["nu1", "nu2", "T1", "T2", "m"],
            lambda d, m: [(d["T1"], 0, 0), (d["T2"], m, m)]),
}

#          start   lower   upper
BOUNDS = {
    "nu1": (1,     1e-3,   100),
    "nu2": (1,     1e-3,   100),
    "nuG": (1,     1e-3,   100),
    "T":   (.5,    1e-3,   20),
    "T1":  (.5,    1e-3,   20),
    "T2":  (.5,    1e-3,   20),
    "dT":  (.5,    1e-3,   20),
    "m":   (1,     1e-4,   200),
    "m12": (1,     1e-4,   200),
    "m21": (1,     1e-4,   200),
    "me":  (.1,    1e-4,   200),
    "mG":  (.5,    1e-4,   200),
    "P":   (.5,    1e-3,   .999),
    "hrf": (.25,   1e-2,   1),
    "Q":   (.25,   1e-3,   .999),
}


def make_model(base, het):
    """Model function and parameter names for a base history + heterogeneity."""
    pnames, epochs = BASE[base]
    two_m, two_n = "2m" in het, "2N" in het
    names = list(pnames) + (["me", "P"] if two_m else []) + \
        (["hrf", "Q"] if two_n else [])

    def func(params, ns):
        d = dict(zip(names, params))
        m = d.get("m", 0)
        mig = [(1.0, m)] if not two_m else [(d["P"], m), (1 - d["P"], d["me"])]
        ne = [(1.0, 1.0)] if not two_n else [(1 - d["Q"], 1.0),
                                             (d["Q"], d["hrf"])]
        fs = None
        for wm, mm in mig:
            for wn, scale in ne:
                part = _integrate(ns, d["nu1"] * scale, d["nu2"] * scale,
                                  epochs(d, mm)) * (wm * wn)
                fs = part if fs is None else fs + part
        return fs

    return func, names


# Reference size for the demes graph of the ghost models. demes works in
# generations and per-generation migration fractions, which must be <= 1 (and
# sum to <= 1 into each deme). With an ancestral size of 1 a scaled rate M =
# 2*Nref*m becomes m = M/2, so any M above 2 was rejected ("invalid migration")
# and most starting points failed. A large Nref keeps every rate inside the
# bounds tiny; moments rescales to the same coalescent units, so the spectrum
# is unchanged (checked against the Nref = 1 graph and against direct IM
# integration with mG ~ 0).
GHOST_NREF = 1e4


def make_ghost(sc):
    """IM or SC between pop1 and pop2, plus an unsampled ghost population.

    Built with demes and integrated with moments.Spectrum.from_demes; the ghost
    is simulated and marginalised out. The ancestral size is GHOST_NREF, so
    sizes are nu * GHOST_NREF, times 2 * GHOST_NREF * T generations and
    per-generation rates M / (2 * GHOST_NREF).
    """
    import demes
    import moments
    N = GHOST_NREF
    names = (["nu1", "nu2", "nuG", "T1", "T2", "dT", "m", "mG"] if sc else
             ["nu1", "nu2", "nuG", "T", "dT", "m", "mG"])

    def func(params, ns):
        d = dict(zip(names, params))
        tsplit = d["T1"] + d["T2"] if sc else d["T"]
        b = demes.Builder(time_units="generations")
        b.add_deme("ANC", epochs=[dict(start_size=N, end_time=2 * N * tsplit)])
        b.add_deme("G", ancestors=["ANC"],
                   start_time=2 * N * (tsplit + d["dT"]),
                   epochs=[dict(start_size=d["nuG"] * N)])
        b.add_deme("pop1", ancestors=["ANC"],
                   epochs=[dict(start_size=d["nu1"] * N)])
        b.add_deme("pop2", ancestors=["ANC"],
                   epochs=[dict(start_size=d["nu2"] * N)])
        if sc:
            b.add_migration(demes=["pop1", "pop2"], rate=d["m"] / (2 * N),
                            start_time=2 * N * d["T2"], end_time=0)
        else:
            b.add_migration(demes=["pop1", "pop2"], rate=d["m"] / (2 * N))
        b.add_migration(source="G", dest="pop1", rate=d["mG"] / (2 * N),
                        start_time=2 * N * tsplit, end_time=0)
        b.add_migration(source="G", dest="pop2", rate=d["mG"] / (2 * N),
                        start_time=2 * N * tsplit, end_time=0)
        return moments.Spectrum.from_demes(
            b.resolve(), sampled_demes=["pop1", "pop2"],
            sample_sizes=[int(ns[0]), int(ns[1])])

    return func, names


MODELS = {}
for _b in BASE:
    MODELS[_b] = make_model(_b, "")
    MODELS[_b + "2N"] = make_model(_b, "2N")
    if _b in ("IM", "AM", "SC"):
        MODELS[_b + "2m"] = make_model(_b, "2m")
        MODELS[_b + "2m2N"] = make_model(_b, "2m2N")
GHOSTS = ["IMghost", "SCghost"]
ALL = list(MODELS) + GHOSTS


def get_model(name):
    if name == "IMghost":
        return make_ghost(False)
    if name == "SCghost":
        return make_ghost(True)
    return MODELS[name]


# ------------------------------------------------------------------- input ---
def read_spectra(path, sizes, many=False):
    shape = tuple(2 * n + 1 for n in sizes)
    want = int(np.prod(shape))
    out = []
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            if not line.split():
                continue
            vals = np.array([float(x) for x in line.split()])
            if vals.size != want:
                sys.exit(f"{prog}: {path} line {i} has {vals.size} values but "
                         f"-n {' '.join(map(str, sizes))} implies "
                         f"{' x '.join(map(str, shape))} = {want}. -n is "
                         f"individuals, not chromosomes.")
            out.append(vals.reshape(shape))
            if not many:
                break
    if not out:
        sys.exit(f"{prog}: {path} is empty")
    return out if many else out[0]


def prepare(raw, project):
    import moments
    fs = moments.Spectrum(raw)
    if project:
        fs = fs.project([2 * n for n in project])
    return fs.fold()


# ---------------------------------------------------------------- godambe ---
def godambe(func, popt, data, boots, eps):
    """CLAIC and Godambe SEs for a fitted model; theta is the last parameter."""
    import moments
    from moments.Godambe import _get_godambe
    ns = data.sample_sizes
    theta = moments.Inference.optimal_sfs_scaling(func(popt, ns), data)
    p0 = list(popt) + [theta]
    fex = lambda p, n: p[-1] * func(p[:-1], n)
    gim, H, J, _ = _get_godambe(fex, boots, p0, data, eps, log=False)
    if np.any(np.linalg.eigvalsh((H + H.T) / 2) <= 0):
        raise ValueError("Hessian is not positive definite, so the best run is "
                         "not at a local maximum (or sits on a bound); raise "
                         "-r / --maxiter or check the bounds")
    ll = moments.Inference.ll(fex(p0, ns), data)
    penalty = float(np.trace(J @ np.linalg.inv(H)))
    se = np.sqrt(np.abs(np.diag(np.linalg.inv(gim))))
    return -2 * ll + 2 * penalty, penalty, se, theta


# ------------------------------------------------------------------- main ---
COLS = ["nu1", "nu2", "nuG", "T", "T1", "T2", "dT",
        "m", "m12", "m21", "me", "mG", "P", "hrf", "Q"]
SCALED = ["Nref", "N1", "N2", "NG", "Tgen", "T1gen", "T2gen", "dTgen",
          "mGen", "m12Gen", "m21Gen", "meGen", "mGGen"]


def main():
    p = argparse.ArgumentParser(
        prog=prog,
        description="Fit two-population divergence models to a folded realSFS "
                    "spectrum with moments; optional Godambe CLAIC and SEs "
                    "from block-bootstrap spectra.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""models:
  {' '.join(MODELS)}
  ghost (with --ghost): {' '.join(GHOSTS)}

examples:
  # every non-ghost model, 20 restarts each, coalescent units only
  {prog} -i mcavC1.C2.sfs -n 51 15 -p C1 C2 -o C1.C2.all

  # one model per job, with Godambe CLAIC/SEs, scaled units and fit plots
  {prog} -i mcavC1.C2.sfs -n 51 15 -p C1 C2 -m SC2m -r 50 --mu 2e-8 \\
      --boot mcavC1.C2.boot.sfs --plot -o C1.C2.SC2m

  # ghost models only, down-projected to keep the 3-population integration fast
  {prog} -i mcavC1.C2.sfs -n 51 15 -p C1 C2 -m IMghost SCghost \\
      --project 20 10 --boot mcavC1.C2.boot.sfs -o C1.C2.ghost

output:
  <out>.tsv           one row per model per run
  <out>.best.tsv      best run of each model, sorted by CLAIC if --boot, else AIC
  <out>.godambe.tsv   with --boot: value and Godambe SE of every parameter
                      (coalescent units; theta included)
  <out>.<M>.pdf       with --plot: data / model / residual panels

notes:
  -n is INDIVIDUALS, in the order given to realSFS.
  deltaAIC / deltaCLAIC are within this invocation. When models are fitted in
  separate jobs, recompute them after combining the .best.tsv files.
  Convergence is the spread of ll across runs: if the top runs of a model
  disagree, raise -r before reading its parameters.
""")
    p.add_argument("-i", "--sfs", required=True, metavar="FILE",
                   help="realSFS 2d spectrum for the two populations")
    p.add_argument("-n", "--sizes", required=True, type=int, nargs=2,
                   metavar="N", help="individuals per population, in the "
                                     "order given to realSFS")
    p.add_argument("-p", "--pops", nargs=2, metavar="NAME",
                   default=["pop0", "pop1"], help="population names")
    p.add_argument("-o", "--out", required=True, metavar="PREFIX",
                   help="output prefix")
    p.add_argument("-m", "--models", nargs="+", choices=ALL, metavar="M",
                   help="models to fit (default: all non-ghost models)")
    p.add_argument("--ghost", action="store_true",
                   help="also fit the ghost-population models")
    p.add_argument("-r", "--runs", type=int, default=20, metavar="N",
                   help="perturbed starting points per model (default 20)")
    p.add_argument("--maxiter", type=int, default=100, metavar="N",
                   help="optimizer iterations per run (default 100)")
    p.add_argument("--fold", type=float, default=1.0, metavar="F",
                   help="how far to perturb each start, in factors of two "
                        "(default 1)")
    p.add_argument("--mu", type=float, metavar="RATE",
                   help="mutation rate per site per generation; without it "
                        "the scaled columns are left empty")
    p.add_argument("-L", "--sites", type=float, metavar="N",
                   help="sites the spectrum represents (default: its total, "
                        "monomorphic included)")
    p.add_argument("--project", type=int, nargs=2, metavar="N",
                   help="down-project to this many individuals per population "
                        "before folding; applied to --boot spectra too")
    p.add_argument("--boot", metavar="FILE",
                   help="realSFS -bootstrap output (one spectrum per line); "
                        "adds CLAIC and Godambe SEs for each model's best run")
    p.add_argument("--eps", type=float, default=.01, metavar="F",
                   help="fractional step for Godambe derivatives (default .01)")
    p.add_argument("--seed", type=int, default=1, help="rng seed (default 1)")
    p.add_argument("--plot", action="store_true",
                   help="write data/model/residual panels per model")
    a = p.parse_args()

    try:
        import moments
    except ImportError:
        sys.exit(f"{prog}: moments is not installed. The population-genetic "
                 f"package is 'moments-popgen' on PyPI; a plain "
                 f"'pip install moments' installs something else.\n"
                 f"  python3 -m pip install --user moments-popgen demes")

    models = a.models or (list(MODELS) + (GHOSTS if a.ghost else []))
    if a.project and any(a.project[i] > a.sizes[i] for i in (0, 1)):
        sys.exit(f"{prog}: --project {a.project[0]} {a.project[1]} is larger "
                 f"than -n {a.sizes[0]} {a.sizes[1]}")

    raw = read_spectra(a.sfs, a.sizes)
    L = a.sites if a.sites else raw.sum()
    data = prepare(raw, a.project)
    ns = data.sample_sizes
    print(f"{prog}: {a.sfs}  {a.pops[0]} x {a.pops[1]}  {ns[0]} x {ns[1]} "
          f"chromosomes  {L:.6g} sites, {data.S():.6g} segregating",
          file=sys.stderr)

    boots = None
    if a.boot:
        boots = [prepare(b, a.project)
                 for b in read_spectra(a.boot, a.sizes, many=True)]
        print(f"{prog}: {len(boots)} bootstrap spectra from {a.boot}",
              file=sys.stderr)

    np.random.seed(a.seed)
    rows = []
    for name in models:
        func, pnames = get_model(name)
        start = [BOUNDS[n][0] for n in pnames]
        lo = [BOUNDS[n][1] for n in pnames]
        hi = [BOUNDS[n][2] for n in pnames]

        def model_func(params, nsamp, _f=func):
            return _f(params, nsamp).fold()

        for run in range(1, a.runs + 1):
            pstart = moments.Misc.perturb_params(start, fold=a.fold,
                                                 lower_bound=lo,
                                                 upper_bound=hi)
            try:
                popt = moments.Inference.optimize_log(
                    pstart, data, model_func, lower_bound=lo, upper_bound=hi,
                    maxiter=a.maxiter, verbose=0)
                mod = model_func(popt, ns)
                ll = moments.Inference.ll_multinom(mod, data)
                theta = moments.Inference.optimal_sfs_scaling(mod, data)
            except Exception as err:                         # noqa: BLE001
                print(f"{prog}: {name} run {run} failed: {err}",
                      file=sys.stderr)
                continue
            if not np.isfinite(ll):
                print(f"{prog}: {name} run {run} gave a non-finite likelihood",
                      file=sys.stderr)
                continue
            k = len(pnames)
            rows.append(dict(model=name, run=run, ll=ll, k=k,
                             AIC=2 * k - 2 * ll, theta=theta,
                             popt=dict(zip(pnames, popt)), func=model_func,
                             pnames=pnames, CLAIC=None, pen=None, se=None))
        done = [r for r in rows if r["model"] == name]
        if done:
            best = max(done, key=lambda r: r["ll"])
            print(f"{prog}: {name}\t{len(done)}/{a.runs} runs\t"
                  f"best ll {best['ll']:.2f}\tAIC {best['AIC']:.2f}",
                  file=sys.stderr)
        else:
            print(f"{prog}: {name}\tno run converged", file=sys.stderr)

    if not rows:
        sys.exit(f"{prog}: nothing converged; check -n and the spectrum")

    bests = [max((r for r in rows if r["model"] == m), key=lambda r: r["ll"])
             for m in dict.fromkeys(r["model"] for r in rows)]

    for r in bests:
        at = [n for n in r["pnames"]
              if r["popt"][n] <= BOUNDS[n][1] * 1.1
              or r["popt"][n] >= BOUNDS[n][2] * .9]
        if at:
            print(f"{prog}: {r['model']} best run is at a bound for "
                  f"{' '.join(at)}", file=sys.stderr)

    failed = []
    if boots:
        for r in bests:
            try:
                r["CLAIC"], r["pen"], r["se"], _ = godambe(
                    r["func"], [r["popt"][n] for n in r["pnames"]], data,
                    boots, a.eps)
            except Exception as err:                         # noqa: BLE001
                print(f"{prog}: Godambe failed for {r['model']}: {err}",
                      file=sys.stderr)
                failed.append(r["model"])

    bestAIC = min(r["AIC"] for r in rows)
    cl = [r["CLAIC"] for r in bests if r["CLAIC"] is not None]
    bestCL = min(cl) if cl else None

    hdr = (["pop1", "pop2", "model", "run", "ll", "k", "AIC", "deltaAIC",
            "CLAIC", "deltaCLAIC", "effK", "theta"] + COLS + SCALED)

    def fmt(r):
        pv = r["popt"]
        out = [a.pops[0], a.pops[1], r["model"], r["run"], f"{r['ll']:.4f}",
               r["k"], f"{r['AIC']:.4f}", f"{r['AIC'] - bestAIC:.4f}"]
        if r["CLAIC"] is not None:
            out += [f"{r['CLAIC']:.4f}", f"{r['CLAIC'] - bestCL:.4f}",
                    f"{r['pen']:.4f}"]
        else:
            out += ["", "", ""]
        out += [f"{r['theta']:.6g}"]
        out += [f"{pv[c]:.6g}" if c in pv else "" for c in COLS]
        if a.mu:
            nref = r["theta"] / (4 * a.mu * L)
            def g(c, f):
                return f"{f(pv[c]):.6g}" if c in pv else ""
            size = lambda v: v * nref
            gens = lambda v: v * 2 * nref
            rate = lambda v: v / (2 * nref)
            out += [f"{nref:.6g}", g("nu1", size), g("nu2", size),
                    g("nuG", size), g("T", gens), g("T1", gens),
                    g("T2", gens), g("dT", gens), g("m", rate),
                    g("m12", rate), g("m21", rate), g("me", rate),
                    g("mG", rate)]
        else:
            out += [""] * len(SCALED)
        return "\t".join(str(x) for x in out)

    rows.sort(key=lambda r: (r["model"], -r["ll"]))
    with open(a.out + ".tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write(fmt(r) + "\n")

    bests.sort(key=lambda r: (r["CLAIC"] if r["CLAIC"] is not None
                              else float("inf"), r["AIC"]))
    with open(a.out + ".best.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in bests:
            fh.write(fmt(r) + "\n")

    if boots:
        with open(a.out + ".godambe.tsv", "w") as fh:
            fh.write("pop1\tpop2\tmodel\tparam\tvalue\tse\n")
            for r in bests:
                if r["se"] is None:
                    continue
                vals = [r["popt"][n] for n in r["pnames"]] + [r["theta"]]
                for n, v, s in zip(r["pnames"] + ["theta"], vals, r["se"]):
                    fh.write(f"{a.pops[0]}\t{a.pops[1]}\t{r['model']}\t{n}\t"
                             f"{v:.6g}\t{s:.6g}\n")

    print(f"\n{prog}: {a.out}.tsv ({len(rows)} rows), {a.out}.best.tsv"
          + (f", {a.out}.godambe.tsv" if boots else ""))
    print("\t".join(["model", "ll", "AIC", "CLAIC", "effK"]))
    for r in bests:
        print(f"{r['model']}\t{r['ll']:.2f}\t{r['AIC']:.2f}\t"
              + (f"{r['CLAIC']:.2f}\t{r['pen']:.2f}" if r["CLAIC"] is not None
                 else "\t"))

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for r in bests:
            mod = r["func"]([r["popt"][n] for n in r["pnames"]], ns)
            moments.Plotting.plot_2d_comp_multinom(
                mod, data, resid_range=3, pop_ids=a.pops, show=False)
            plt.gcf().savefig(f"{a.out}.{r['model']}.pdf")
            plt.close("all")
            print(f"{prog}: {a.out}.{r['model']}.pdf")

    if failed:
        sys.exit(f"{prog}: Godambe failed for {' '.join(failed)}; their CLAIC "
                 f"and SEs are empty")


if __name__ == "__main__":
    main()
