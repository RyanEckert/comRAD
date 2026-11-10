"""What ran: tool versions, modules, comRAD itself, and the conf, recorded in the run diary.

record(pipeline) is called when a driver starts (and by `comrad step`). It writes
run/<pipeline>/versions.txt (tool, version, path; overwritten each time) and puts in the
diary either every line (first time) or only what changed since the last versions.txt.
The settings in effect (defaults, cluster.config and the config, merged) are saved as
run/<pipeline>/<pipeline>.config.<date> whenever they differ from
the last copy, and the diary gets one line per setting that changed.
"""
import glob
import hashlib
import importlib.metadata as md
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

SCRIPTS = ["trim2bRAD.py", "sampleRename.py", "fetchReads.py", "readCounts.py", "bowtie2Launcher.py",
           "splitHostSym.py", "packMetagenome.py", "plotQC.py", "filterBams.py", "paralogFilter.py",
           "pruneSnps.py", "pcangsdFisToText.py", "ngsAdmixLauncher.py", "pickQ.py", "pcangsdOutliers.py",
           "subsetBeagle.py", "momentsModels.py", "zipper.py", "mind", "cloneGroups.py", "admixResiduals.py",
           "cloneTestRun.py", "cloneTestScore.py", "simCloneKin.py", "plotStep.py", "kmeansBic.py"]
TOOLS = {
    "reads": ["cutadapt"],
    "popgen": ["cutadapt", "bowtie2", "samtools", "angsd", "realSFS", "thetaStat", "pcangsd", "NGSadmix",
               "evalAdmix", "NGSremix", "ngsParalog"],
    "demography": ["angsd", "realSFS"],
}
PY = {"reads": [], "popgen": ["numpy", "scipy", "matplotlib"],
      "demography": ["numpy", "scipy", "matplotlib", "moments-popgen", "demes"]}
RPKG = {}          # R is no longer used by the pipelines
# a line that names a version, or is nothing but one ("5.2", "v1.36.4", "samtools 1.19.2")
VLINE = re.compile(r"version", re.I)
VONLY = re.compile(r"^(\S+ )?v?\d+(\.\d+)+\S*$")


def short(path, n=12):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def script_sum(path):
    """checksum of a script without its first line (pip rewrites the #! line)"""
    with open(path, "rb") as fh:
        fh.readline()
        return hashlib.sha1(fh.read()).hexdigest()[:12]


def run(cmd, timeout=15):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, stdin=subprocess.DEVNULL)
        return (r.stdout or "") + (r.stderr or "")
    except (OSError, subprocess.TimeoutExpired):
        return ""


def tool_version(t, found):
    path = shutil.which(t)
    if not path:
        return "not on PATH", ""
    out = run([t, "--version"])
    for line in out.splitlines():
        s = line.strip().lstrip("->").strip()
        if s and (VLINE.search(s) or VONLY.match(s)) and "unknown" not in s.lower() \
                and "--version" not in s:          # tools that echo the flag back in an error
            s = re.sub(r"^\S*/", "", s)            # bowtie2 prints its own full path first
            return s[:120], path
    # realSFS and thetaStat have no version flag but come with angsd
    if t in ("realSFS", "thetaStat") and "angsd" in found and found["angsd"][1] \
            and os.path.dirname(found["angsd"][1]) == os.path.dirname(path):
        return f"from the same build as angsd ({found['angsd'][0]})", path
    st = os.stat(path)
    return f"no version flag; sha1 {short(path)}, modified {time.strftime('%Y-%m-%d', time.localtime(st.st_mtime))}", path


def collect(pipeline):
    rows = []
    try:
        rv = md.version("comrad")
    except md.PackageNotFoundError:
        rv = "?"
    h = hashlib.sha1()
    for f in sorted(glob.glob(os.path.join(HERE, "**", "*"), recursive=True)):
        if os.path.isfile(f) and "__pycache__" not in f:
            h.update(os.path.relpath(f, HERE).encode())
            h.update(open(f, "rb").read())
    rows.append(("comrad", f"{rv}, files {h.hexdigest()[:12]}", HERE))
    rows.append(("python", sys.version.split()[0], sys.executable))
    for p in PY.get(pipeline, []):
        try:
            rows.append((f"python:{p}", md.version(p), ""))
        except md.PackageNotFoundError:
            rows.append((f"python:{p}", "not installed", ""))
    found = {}
    for t in TOOLS.get(pipeline, []):
        found[t] = tool_version(t, found)
        rows.append((t, *found[t]))
    if RPKG.get(pipeline) and shutil.which("Rscript"):
        pk = ", ".join(f'"{p}"' for p in RPKG[pipeline])
        out = run(["Rscript", "-e", f'for (p in c({pk})) cat(p, tryCatch(as.character(packageVersion(p)), '
                                    f'error = function(e) "not installed"), "\\n")'], timeout=60)
        for line in out.splitlines():
            f = line.split()
            if len(f) >= 2 and f[0] in RPKG[pipeline]:
                rows.append((f"R:{f[0]}", " ".join(f[1:]), ""))
    for s in SCRIPTS:
        p = shutil.which(s)
        rows.append((s, f"sha1 {script_sum(p)}" if p else "not on PATH", p or ""))
    mods = [m for m in os.environ.get("LOADEDMODULES", "").split(":") if m]
    rows.append(("modules", " ".join(mods) if mods else "none recorded (LOADEDMODULES is empty)", ""))
    return rows


def read_versions(path):
    old = {}
    if os.path.exists(path):
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#"):
                old[f[0]] = f[1]
    return old


def render(settings, pl):
    """merged settings as a config file: top-level keys, then one block per section"""
    cl = pl.cluster_file()
    out = [f"# settings in effect: defaults, {cl or 'no cluster.config'}, {pl.conf_path}"]
    secs = {}
    for k, v, _ in settings:
        if "." in k:
            sec, key = k.split(".", 1)
            secs.setdefault(sec, []).append(f"{key} = {v}")
        else:
            out.append(f"{k} = {v}")
    for sec, lines in secs.items():
        out += ["", f"[{sec}]"] + lines
    return "\n".join(out) + "\n"


def record(pl, parse_conf):
    """pl: the comrad Pipeline (uses .name .dir .conf_path .log)"""
    vpath = os.path.join(pl.dir, "versions.txt")
    old = read_versions(vpath)
    rows = collect(pl.name)
    with open(vpath, "w") as fh:
        fh.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}  host {os.uname().nodename}  "
                 f"job {os.environ.get('SLURM_JOB_ID', '-')}\n# tool\tversion\tpath\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    if not old:
        for name, v, _ in rows:
            pl.log(f"version\t{name}: {v}")
    else:
        changed = [(n, v) for n, v, _ in rows if old.get(n) != v]
        for n, v in changed:
            pl.log(f"version\t{n}: {v} (was {old.get(n, 'not recorded')})")
        if not changed:
            pl.log(f"version\tunchanged since the last run ({os.path.relpath(vpath, pl.wd)})")

    # the settings in effect for this run: defaults, cluster.config and the config, merged
    if os.path.exists(pl.conf_path):
        copies = sorted(glob.glob(os.path.join(pl.dir, f"{pl.name}.conf.*")) +
                        glob.glob(os.path.join(pl.dir, f"{pl.name}.config.*")), key=os.path.getmtime)
        now = render(pl.settings(), pl)
        if not copies or open(copies[-1]).read() != now:
            dst = os.path.join(pl.dir, f"{pl.name}.config.{time.strftime('%Y%m%d-%H%M%S')}")
            open(dst, "w").write(now)
            if copies:
                a = {k: v for k, v, _ in parse_conf(copies[-1])}
                b = {k: v for k, v, _ in parse_conf(dst)}
                for k in sorted(set(a) | set(b)):
                    if a.get(k) != b.get(k):
                        pl.log(f"conf\t{k}: '{a.get(k, '')}' -> '{b.get(k, '')}'")
            pl.log(f"conf\tsaved as {os.path.relpath(dst, pl.wd)}")
