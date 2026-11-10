#!/usr/bin/env python3
"""Stand-ins for the HPC tools, for testing the pipelines' wiring.

Each tool is a symlink to this file; argv[0] picks the behaviour. Outputs have
the names and formats the steps read, filled with small synthetic data: five
contigs of 200 sites, every 10th site a SNP, lineages by sample number parity.
"""
import gzip, hashlib, os, random, re, shutil, subprocess, sys

tool = os.path.basename(sys.argv[0])
args = sys.argv[1:]
FAIL = os.environ.get("FAKE_FAIL", "")          # e.g. "angsd:once" to fail the first call


def opt(flag, default=None, n=1):
    if flag in args:
        i = args.index(flag)
        return args[i + 1] if n == 1 else args[i + 1:i + 1 + n]
    return default


def gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def lines(path):
    return [l.strip() for l in open(path) if l.strip()]


def name(b):
    return os.path.basename(b).removesuffix(".bam")


def lineage(s):
    """Cluster by the sample's number: odd -> 1, even -> 2; A01 is admixed."""
    m = re.search(r"(\d+)", s)
    k = int(m.group(1)) if m else 0
    return 0 if s.startswith("SA01") else (1 if k % 2 else 2)


SITES = [(f"chr{c}", p) for c in range(1, 6) for p in range(1, 201)]
SNPS = [s for s in SITES if s[1] % 10 == 0]


def maybe_fail():
    for f in FAIL.split(","):
        if f.startswith(tool + ":"):
            mark = f"/tmp/fakefail.{tool}.{os.getuid()}"
            if not os.path.exists(mark):
                open(mark, "w").close()
                if f.endswith(":timeout"):
                    print("slurmstepd: error: *** JOB 1 ON node CANCELLED AT now DUE TO TIME LIMIT ***", file=sys.stderr)
                    os.killpg(os.getpgid(0), 9)    # like SLURM: the whole job is killed, no footer gets written
                if f.endswith(":vanish"):              # node dies: no message, no footer
                    os.killpg(os.getpgid(0), 9)
                if f.endswith(":oom"):
                    print("slurmstepd: error: Detected 1 oom-kill event(s) in StepId=1.batch.", file=sys.stderr)
                    os.killpg(os.getpgid(0), 9)
                sys.exit(3)


def main():
    maybe_fail()
    if tool == "sbatch":
        return sbatch()
    if tool == "squeue":
        return
    if tool == "launcher_creator.py":
        j, n = opt("-j"), opt("-n")
        open(n + ".slurm", "w").write(
            f"#!/bin/bash\nid=${{SLURM_JOB_ID:-$RANDOM}}\nexec >>{n}.o$id 2>>{n}.e$id\nrc=0\n"
            f"while IFS= read -r l; do bash -c \"$l\" </dev/null || rc=1; done <{j}\nexit 0\n")
        return
    if tool == "mind":
        import math
        print(math.ceil(len(lines(args[0])) * float(args[1] if len(args) > 1 else .75) - 1e-9))
        return
    globals()["t_" + re.sub(r"\W", "_", tool.removesuffix(".py"))]()


def sbatch():
    wrap = opt("--wrap")
    o, e = opt("-o", "slurm-%j.out"), opt("-e", "slurm-%j.out")
    jid = str(random.randint(100000, 999999))
    o, e = o.replace("%j", jid).replace("%A", jid), e.replace("%j", jid).replace("%A", jid)
    print(f"Submitted batch job {jid}", flush=True)
    if os.environ.get("FAKE_SBATCH_LOG"):
        with open(os.environ["FAKE_SBATCH_LOG"], "a") as fh:
            fh.write(" ".join(x for x in args[:args.index("--wrap")] if x != "--wrap") + "\n" if "--wrap" in args else " ".join(args) + "\n")
    ex = [x for x in args if x.startswith("--exclude=")]
    if ex:
        print(f"(fake sbatch: {ex[0]})", flush=True)
    arr = [x for x in args if x.startswith("--array=")]
    if arr and wrap is not None:
        n = int(re.match(r"--array=1-(\d+)", arr[0]).group(1))
        rc = 0
        for i in range(1, n + 1):
            env = {**os.environ, "SLURM_ARRAY_TASK_ID": str(i), "SLURM_ARRAY_JOB_ID": jid, "SLURM_JOB_ID": f"{jid}{i}"}
            with open(o.replace("%a", str(i)), "a") as fo, open(e.replace("%a", str(i)), "a") as fe:
                rc = max(rc, subprocess.run(["bash", "-c", wrap], stdout=fo, stderr=fe, env=env, start_new_session=True).returncode)
        sys.exit(rc)
    if wrap is not None:
        with open(o, "a") as fo, open(e, "a") as fe:
            rc = subprocess.run(["bash", "-c", wrap], stdout=fo, stderr=fe, env={**os.environ, "SLURM_JOB_ID": jid}, start_new_session=True).returncode
    else:
        rc = subprocess.run(["bash", args[-1]], env={**os.environ, "SLURM_JOB_ID": jid}).returncode
    sys.exit(rc)


# ---------------------------------------------------------------- reads ---
def t_fetchReads():
    for p in ("pool1", "pool2", "pool3", "pool4"):
        gz(f"{p}.fastq.gz", "@r\nACGT\n+\nIIII\n")


def t_trim2bRAD():
    if opt("-f"):
        suf, n = opt("-f"), opt("-n", "trims")
        os.makedirs("logs", exist_ok=True)
        with open(n, "w") as fh:
            for f in sorted(x for x in os.listdir(".") if x.endswith(suf) and ".demux." not in x):
                fh.write(f"trim2bRAD.py -i {f} --gzip > logs/{f.split('.')[0]}.demux.log 2>&1\n")
        return
    f = opt("-i")
    pool = f.split(".")[0]
    for bc in ("AAAA", "CCCC", "GGGG", "TTTT"):
        gz(f"{pool}_{bc}.demux.fastq.gz", "@r\nACGT\n+\nIIII\n")
    print(f"{pool}: 1000 reads in, 950 kept (95.0%)")


def t_sampleRename():
    sheet = opt("-i")
    if not os.path.exists(sheet):
        sheet += ".csv"
    pairs = [l.split(",") for l in lines(sheet)]
    ext = opt("-f")
    for old, new in pairs:
        src = f"{old}.{ext}"
        if not os.path.exists(src):
            sys.exit(f"sampleRename: no {src}")
        if "--dry-run" in args:
            print(f"{src} -> {new}.{ext}")
        else:
            os.rename(src, f"{new}.{ext}")


# --------------------------------------------------------------- popgen ---
def t_cutadapt():
    out, inp = opt("-o"), args[-1]
    shutil.copy(inp, out)
    print("Total written (filtered):  1,000 bp (90%)")


def t_readCounts():
    o = opt("-o")
    with open(o + "Counts", "w") as fh:
        fh.write("file\treads\n")
        for f in sorted(x for x in os.listdir(".") if x.endswith("." + opt("-f"))):
            fh.write(f"{f}\t1000\n")


def t_bowtie2Launcher():
    g, n, tag, suf = opt("-g"), opt("-n"), opt("-s"), opt("-f")
    if "--plain" not in args or not tag:
        sys.exit("fake bowtie2Launcher: the pipeline should call --plain with -s")
    os.makedirs("logs", exist_ok=True)
    with open(n, "w") as fh:
        for f in sorted(x for x in os.listdir(".") if x.endswith(suf)):
            s = f[:-len(suf)]
            fh.write(f"bowtie2 -x {g} -U {f} --local > {s}.{tag}.bam 2> logs/{s}.{tag}.bowtie2.log\n")


def t_bowtie2():
    print("@HD\nr1\t0\tchr1\t10")
    print("95.00% overall alignment rate", file=sys.stderr)


def t_splitHostSym():
    cmd = args[0]
    if "--sample" not in args:
        n = opt("-n")
        os.makedirs("logs", exist_ok=True)
        pat = opt("-H")
        suf = pat.split("*")[-1]
        with open(n, "w") as fh:
            for f in sorted(x for x in os.listdir(".") if x.endswith(suf)):
                s = f[:-len(suf)]
                rest = " ".join(a for a in args[1:] if a not in ("-n", n))
                fh.write(f"splitHostSym.py {cmd} {rest} --sample {s} > logs/{s}.{cmd}log 2>&1\n")
        return
    s, out = opt("--sample"), opt("--outdir")
    os.makedirs(out, exist_ok=True)
    if cmd == "scan":
        print(f"{s}\t1000\t200\t10\t5\t100\t3\t60")
    else:
        for x in (f"{s}.bam", f"{s}.bam.bai", f"{s}.sym.bam", f"{s}.sym.bam.bai"):
            open(os.path.join(out, x), "w").write("bam\n")
        if s.endswith(("3", "7")):
            print("[bam_sort_core] merging from 1 files and 1 in-memory blocks...", file=sys.stderr)
        print(f"{s}\t1000\t200\t10\t980\t190\t5\t2\t8\t160")


def t_samtools():
    sub = args[0]
    if sub == "view":
        if "-c" in args:
            print(1000)
        else:
            for i in range(30):
                print(f"r{i}\t0\tSym{i % 3 + 1}\t10")
    elif sub == "sort":
        open(opt("-o"), "w").write("bam\n")
    elif sub == "index":
        open(args[-1] + ".bai", "w").write("bai\n")


def t_ngsParalog():
    open(opt("-outfile"), "w").write("lr\n")
    print("Finished!", file=sys.stderr)


def t_paralogFilter():
    d, sites, o = opt("--collect"), opt("-s"), opt("-o")
    for f in os.listdir(d):
        if f.endswith(".log") and "Finished!" not in open(os.path.join(d, f)).read():
            sys.exit(f"paralogFilter: {f} not finished")
    s = lines(sites)
    open(o + "ParalogSites", "w").write("\n".join(s[:5]) + "\n")
    open(o + "KeepSites", "w").write("\n".join(s[5:]) + "\n")
    open(o + "KeepRegions", "w").write("\n".join(sorted({x.split()[0] for x in s})) + "\n")


def samples_of(listfile):
    return [name(b) for b in lines(listfile)]


def t_angsd():
    if args[:2] == ["sites", "index"]:
        open(args[2] + ".idx", "w").write("idx\n")
        open(args[2] + ".bin", "w").write("bin\n")
        return
    out = opt("-out")
    bams = lines(opt("-b")) if opt("-b") else [opt("-i")]
    miss = [b for b in bams if not os.path.exists(b.split()[0])]
    if miss:                                  # as angsd: a bam that is not there (e.g. a path relative to another folder)
        sys.exit(f"angsd: cannot open {miss[0]} (from {os.getcwd()})")
    if opt("-b"):
        sm = samples_of(opt("-b"))
    else:
        sm = [name(opt("-i"))]
    n = len(sm)
    sites = SITES
    if opt("-sites"):
        keep = {tuple(l.split()[:2]) for l in lines(opt("-sites"))}
        sites = [s for s in SITES if (s[0], str(s[1])) in keep]
    snps = [s for s in sites if s[1] % 10 == 0]
    print(f"\t-> filtering: {len(snps)} sites retained", file=sys.stderr)
    if "-doDepth" in args:
        open(out + ".depthSample", "w").write("".join(" ".join(["0"] + ["5"] * 20) + "\n" for _ in sm))
        open(out + ".qs", "w").write("qscore\tcounts\n30\t100\n")
        gz(out + ".pos.gz", "chr\tpos\ttotDepth\n" + "".join(f"{c}\t{p}\t{5 * n}\n" for c, p in SITES))
        open(out + ".counts", "w").write("\t".join(f"ind{i}TotDepth" for i in range(n)) + "\t\n" +
                                         "".join("\t".join(["5"] * n) + "\t\n" for _ in SITES))
    if "-doMaf" in args:
        gz(out + ".mafs.gz", "chromo\tposition\tmajor\tminor\tref\tknownEM\tnInd\n" +
           "".join(f"{c}\t{p}\tA\tG\tA\t{.1 + (p % 7) / 20:.3f}\t{n}\n" for c, p in snps))
    if "-doGlf" in args:
        gz(out + ".beagle.gz", "marker\tallele1\tallele2\t" + "\t".join(f"Ind{i}\tInd{i}\tInd{i}" for i in range(n)) + "\n" +
           "".join(f"{c}_{p}\t0\t2\t" + "\t".join(".3\t.4\t.3" for _ in range(n)) + "\n" for c, p in snps))
    if "-doIBS" in args:
        rnd = random.Random(1)
        M = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                a, b = sm[i].split(".")[0], sm[j].split(".")[0]
                d = .04 + rnd.random() * .01 if a == b else .30 + rnd.random() * .05
                M[i][j] = M[j][i] = d
        open(out + ".ibsMat", "w").write("".join(" ".join(f"{v:.4f}" for v in r) + "\n" for r in M))
        open(out + ".covMat", "w").write("".join(" ".join("0.1" for _ in range(n)) + "\n" for _ in range(n)))
    if "-doGeno" in args:
        # genotypes shared within a sample prefix (clones and replicates), random otherwise; -doGeno 8: three
        # posteriors per sample; 15: major, minor, then per sample call, call as letters, three posteriors
        mode = int(opt("-doGeno"))
        rows = []
        for c, p in snps:
            f = [c, str(p)] + (["A", "G"] if mode & 1 else [])
            for s in sm:
                g = int(hashlib.md5(f"{s.split('.')[0]}{c}{p}".encode()).hexdigest(), 16) % 3
                post = ["0.980000" if k == g else "0.010000" for k in range(3)]
                f += ([str(g)] if mode & 2 else []) + ([["AA", "AG", "GG"][g]] if mode & 4 else []) + (post if mode & 8 else [])
            rows.append("\t".join(f) + "\n")
        gz(out + ".geno.gz", "".join(rows))
    if "-dumpCounts" in args:
        # depth 12 everywhere, but the first sample is below 10 at every 7th site
        gz(out + ".pos.gz", "chr\tpos\ttotDepth\n" + "".join(f"{c}\t{p}\t{12 * n}\n" for c, p in sites))
        gz(out + ".counts.gz", "".join(f"ind{i}TotDepth\t" for i in range(n)) + "\n" +
           "".join("".join(f"{8 if i == 0 and p % 7 == 0 else 12}\t" for i in range(n)) + "\n" for c, p in sites))
    if "-doSaf" in args:
        open(out + ".saf.idx", "w").write(f"n={n}\n")
        gz(out + ".saf.gz", "saf\n")
        gz(out + ".saf.pos.gz", "".join(f"{c}\t{p}\n" for c, p in sites))


def nsaf(idx):
    return int(open(idx).read().split("=")[1])


def em_log(n=1):
    for _ in range(n):
        print("\t-> maxiter:1000\nlik[57]=-123.4 diff=1.0e-11\nlikelihood: -123.4", file=sys.stderr)


def t_realSFS():
    if args[:2] == ["fst", "index"]:
        open(opt("-fstout") + ".fst.idx", "w").write("fst\n")
        return
    if args[:2] == ["fst", "stats"]:
        print("0.050000\t0.120000")
        return
    if args[0] == "saf2theta":
        open(opt("-outname") + ".thetas.idx", "w").write("thetas\n")
        return
    idx = [a for a in args if a.endswith(".saf.idx")]
    if args[0] == "print":
        for l in gzip.open(idx[0].removesuffix(".idx") + ".pos.gz", "rt"):
            print(l.rstrip("\n") + "\t0.000000\t-1.0\t-2.0")
        return
    if len(idx) == 1 and nsaf(idx[0]) == 1:          # one sample: heterozygosity
        ns = sum(1 for _ in gzip.open(idx[0].removesuffix(".idx") + ".pos.gz", "rt"))
        print(f"{ns * .8:.2f} {ns * .15:.2f} {ns * .05:.2f}")
        return
    if len(idx) == 2:
        na, nb = nsaf(idx[0]), nsaf(idx[1])
        reps = int(opt("-bootstrap", 1))
        em_log(reps)
        for _ in range(reps):
            print(" ".join(["1"] * ((2 * na + 1) * (2 * nb + 1))))
        return
    n = nsaf(idx[0])
    em_log()
    print(" ".join(["1"] * (2 * n + 1)))


def t_thetaStat():
    p = args[1]
    open(p + ".pestPG", "w").write("#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)\tChr\tWinCenter\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tnSites\n" +
                                   "".join(f"(0,1)(1,2)(0,2)\tchr{c}\t100\t2.0\t1.8\t1\t1\t1\t-0.5\t0\t0\t0\t0\t1000\n" for c in range(1, 6)))


def t_plotQC():
    pref, bams = opt("-p"), lines(opt("-b"))
    i = sys.argv.index("-q") + 1
    qs = []
    while i < len(sys.argv) and not sys.argv[i].startswith("-"):
        qs.append(sys.argv[i]); i += 1
    for q in qs:
        open(f"quality{q}X.txt", "w").write("".join(f"{b} {.3 + (int(hashlib.md5(b.encode()).hexdigest(), 16) % 60) / 100:.4f}\n" for b in bams))
    open(pref + ".pdf", "w").write("pdf\n")
    open(pref + ".info", "w").write("info\n")


def t_filterBams():
    bams, q, p, c, o = lines(opt("-b")), opt("-q"), float(opt("-p")), opt("-c"), opt("-o")
    qv = {name(l.split()[0]): float(l.split()[1]) for l in lines(q)}
    drop = set()
    for g in lines(c):
        m = g.split()
        best = max(m, key=lambda s: qv.get(s, 0))
        drop |= set(m) - {best}
    keep = [b for b in bams if name(b) not in drop and qv.get(name(b), 1) >= p]
    open(o, "w").write("".join(b + "\n" for b in keep))
    print(f"{o}: kept {len(keep)} of {len(bams)}")


def t_pruneSnps():
    m, b, x, o = opt("-m"), opt("-b"), opt("-x"), opt("-o")
    out = set()
    if x and os.path.exists(x):
        out = {tuple(l.split()[:2]) for l in lines(x)}
    rows = gzip.open(b, "rt").read().splitlines()
    keep = [rows[0]] + [r for r in rows[1:] if tuple(r.split("\t")[0].rsplit("_", 1)) not in out]
    gz(o + ".beagle.gz", "\n".join(keep) + "\n")
    open(o + ".sites", "w").write("".join(r.split("\t")[0].replace("_", "\t") + "\n" for r in keep[1:]))


def t_subsetBeagle():
    b, s, k, o = opt("-b"), lines(opt("-s")), lines(opt("-k")), opt("-o")
    idx = [s.index(x) for x in k]
    out = []
    for r in gzip.open(b, "rt").read().splitlines():
        f = r.split("\t")
        out.append("\t".join(f[:3] + [v for i in idx for v in f[3 + 3 * i:6 + 3 * i]]))
    gz(o, "\n".join(out) + "\n")
    open(o + ".order", "w").write("\n".join(k) + "\n")


def beagle_n(b):
    return (len(gzip.open(b, "rt").readline().split("\t")) - 3) // 3


def t_pcangsd():
    b, o = opt("-b"), opt("-o")
    n = beagle_n(b)
    open(o + ".cov", "w").write("".join(" ".join("0.1" for _ in range(n)) + "\n" for _ in range(n)))
    open(o + ".log", "w").write(f"Loaded {n} samples\nMAF filtering\nEM (MAF) converged in 10 iterations\nPCAngsd converged in 12 iterations using 1 eigenvectors.\n")
    if "--admix" in args:
        open(o + ".admix.2.Q", "w").write("".join(("0.9 0.1\n" if i % 2 else "0.1 0.9\n") for i in range(n)))
    if "--inbreed-samples" in args:
        open(o + ".inbreed.samples", "w").write("".join("0.01\n" for _ in range(n)))
    if "--selection" in args:
        open(o + ".selection", "w").write("1.0\n")
        open(o + ".sites", "w").write("1\n")


def t_pcangsdFisToText():
    open(opt("-o"), "w").write("bam\tlineage\tFis\n")


def t_ngsAdmixLauncher():
    f, K, r, n = opt("-f"), int(opt("--maxK")), int(opt("-r")), opt("-n")
    with open(n + "NgsAdmix", "w") as fh:
        for R in range(1, r + 1):
            for k in range(1, K + 1):
                fh.write(f"NGSadmix -likes {f} -K {k} -minMaf 0 -P 10 -o {n}_k{k}_run{R}\n")


def qmatrix(beagle, K, seed):
    """Rows follow the bamlist order behind the beagle: angsd/<sp>Bams, or the subset's .order."""
    order = beagle + ".order"
    if os.path.exists(order):
        sm = [name(x) for x in lines(order)]
    else:
        d = os.path.dirname(beagle) or "."
        sp = os.environ["SP"]
        sm = samples_of(os.path.join(d, f"{sp}Bams"))
    rows = []
    flip = seed % 2 == 0
    for s in sm:
        if K == 1:
            rows.append([1.0])
            continue
        l = lineage(s)
        q = [.02 / (K - 1)] * K
        if l == 0:
            q = [1 / K] * K
        else:
            c = (l - 1) if not flip else (2 - l)
            q[c] = .98
        rows.append(q)
    return rows


def t_NGSadmix():
    b, K, o = opt("-likes"), int(opt("-K")), opt("-o")
    seed = int(opt("-seed", re.search(r"run(\d+)", o).group(1) if re.search(r"run(\d+)", o) else 1))
    Q = qmatrix(b, K, seed)
    open(o + ".qopt", "w").write("".join(" ".join(f"{v:.4f}" for v in r) + "\n" for r in Q))
    gz(o + ".fopt.gz", "0.1\n")
    open(o + ".log", "w").write(f"best like=-{1000 - seed}.5 after 30 iterations\n")


def t_pickQ():
    g = opt("--glob")
    if "--fit" in args:
        ks = sorted(int(re.search(r"_k(\d+)_", f).group(1)) for f in os.listdir(".") if f.startswith(g + "_k") and f.endswith("_best.corres"))
        open(f"{g}_fit.txt", "w").write("K\trmsd\tmean\tmax\tn\tdrop\tinMode\telbow\n" +
                                        "".join(f"{k}\t0.1\t0.1\t0.2\t10\t\t20\t{'<-' if k == 2 else ''}\n" for k in ks))
        return
    lo, hi, ev = int(opt("--min")), int(opt("--max")), opt("--eval")
    rows = []
    for k in range(lo, hi + 1):
        best = f"{g}_k{k}_run1.qopt"
        shutil.copy(best, f"{g}_k{k}_best.qopt")
        rows.append(f"{k}\t20\t1\t20\t20\t1.0\t-1000\t{best}")
    open(f"{g}_bestQ.txt", "w").write("K\truns\tnModes\tmodeSizes\tinMode\tprop\tlik\tbest\n" + "\n".join(rows) + "\n")
    with open(f"{g}_evalAdmixJobs", "w") as fh:
        for k in range(lo, hi + 1):
            fh.write(f"evalAdmix -beagle {ev} -fname {g}_k{k}_run1.fopt.gz -qname {g}_k{k}_run1.qopt -o {g}_k{k}_best.corres -minMaf 0 -P 4\n")


def t_evalAdmix():
    """n x n residual correlations: two lineages lumped at K = 1 leave +/- blocks; K >= 2 fits"""
    b, q = opt("-beagle"), opt("-qname")
    K = int(re.search(r"_k(\d+)_", q).group(1))
    d = os.path.dirname(b) or "."
    sm = samples_of(os.path.join(d, f"{os.environ['SP']}{'Neutral' if 'Neutral' in os.path.basename(b) else ''}Bams"))
    lin = [lineage(s) for s in sm]
    rnd = random.Random(K)
    rows = []
    for i, a in enumerate(lin):
        r = []
        for j, c in enumerate(lin):
            if i == j:
                r.append("NA")
                continue
            v = rnd.uniform(-.002, .002)
            if K == 1 and a and c:
                v += .05 if a == c else -.05
            r.append(f"{v:.5f}")
        rows.append(" ".join(r))
    open(opt("-o"), "w").write("\n".join(rows) + "\n")


def t_NGSremix():
    o = opt("-o")
    sp = os.environ["SP"]
    sm = samples_of(f"../angsd/{sp}Bams")
    n = len(sm)
    out = ["ind1\tind2\tk0\tk1\tk2\tnIter"]
    for i in range(n):
        for j in range(i + 1, n):
            k = (1, 0, 0)
            if (sm[i], sm[j]) == ("S03", "S05"):
                k = (0, 1, 0)                     # theta .25, same lineage
            out.append(f"{i}\t{j}\t{k[0]}\t{k[1]}\t{k[2]}\t10")
    open(o, "w").write("\n".join(out) + "\n")


def t_pcangsdOutliers():
    o = opt("-o")
    open(o + "Selection.tsv", "w").write("chr\tpos\n")
    open(o + "Outliers", "w").write("chr1\t50\nchr2\t100\n" if o.endswith("Cluster1") else "")


def t_momentsModels():
    o, m, s = opt("-o"), opt("-m"), int(opt("--seed"))
    p = opt("-p", n=2)
    rnd = random.Random(hash((m, s)) % 1000)
    base = {"SI": -1900, "IM": -1790, "SC2N": -1639.5}.get(m, -1700)
    H = "pop1\tpop2\tmodel\trun\tll\tk\tAIC\tdeltaAIC\tCLAIC\tdeltaCLAIC\teffK\ttheta"
    rows = []
    for r in range(1, int(opt("-r")) + 1):
        ll = base - rnd.random() * .3
        rows.append([p[0], p[1], m, str(r), f"{ll:.4f}", "4", f"{-2 * ll + 8:.4f}", "0", f"{-2 * ll + 20:.4f}", "0", "10", "1"])
    rows.sort(key=lambda x: -float(x[4]))
    open(o + ".tsv", "w").write(H + "\n" + "\n".join("\t".join(x) for x in rows) + "\n")
    open(o + ".best.tsv", "w").write(H + "\n" + "\t".join(rows[0]) + "\n")
    open(o + ".godambe.tsv", "w").write(f"pop1\tpop2\tmodel\tparam\tvalue\tse\n{p[0]}\t{p[1]}\t{m}\tnu1\t1\t.1\n")
    open(o + ".pdf", "w").write("pdf\n")
    if m == "AM":
        print(f"momentsModels.py: AM best run is at a bound for T2", file=sys.stderr)


if __name__ == "__main__":
    main()
