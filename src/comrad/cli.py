#!/usr/bin/env python3
################################################################################
#Script Name    : comrad
#Description    : Run a comRAD pipeline (reads, popgen, demography) from a project folder
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Run a comRAD pipeline from the folder it works in.

Pipelines:
  reads        project folder: download, trim and demultiplex, rename, split by species
  popgen       species folder: QC, mapping, clones, genotyping, structure, K, kin,
               lineages, outliers, SFS, FST, diversity, results tarball
  demography   species folder, after popgen: bootstrap SFS, moments, results tarball

Each pipeline reads <pipeline>.config in the current folder (`init` writes one to
fill in), after cluster.config (partitions, memory limit, email, extra sbatch
flags) from this folder or up to three above it (`init cluster` writes one).
Every job list has its own [section] in the config: cpus, mem, queue, slurm and
the tool's own options. State is kept in run/<pipeline>/: a .done marker and a
.log per step, and <pipeline>Summary.txt with one line per event. `run` submits
one small driver job (no time limit given) that runs the steps in order.
Steps submit their own jobs through comrad-jobs, which reruns failed jobs by
itself; a step stops the run only when something is wrong that a rerun would
not fix, and pauses it where a decision is needed (K). Submitting again picks
up after the last finished step.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import time

try:
    from comrad import provenance, setup_tools, report
    from comrad import jobs as rtjobs
except ImportError:                      # run from a source checkout
    import provenance
    import report
    import setup_tools
    import jobs as rtjobs

prog = os.path.basename(sys.argv[0])
HERE = os.path.dirname(os.path.abspath(__file__))
PIPES = os.path.join(HERE, "pipelines")
STOP, SKIP, PAUSE = 1, 98, 99


def die(msg):
    sys.exit(f"{prog}: {msg}")


# ------------------------------------------------------------------ config ---
def parse_conf(path, section=""):
    """KEY = value lines; '#' starts a comment; [name] starts a section, whose keys come
    back as 'name.key'. Returns [(key, value, comment)]."""
    out = []
    for n, line in enumerate(open(path), 1):
        s = line.rstrip("\n")
        if not s.strip() or s.lstrip().startswith("#"):
            continue
        m = re.match(r"\s*\[\s*([A-Za-z][A-Za-z0-9_]*)\s*\]\s*(#.*)?$", s)
        if m:
            section = m.group(1)
            continue
        m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$", s)
        if not m:
            die(f"{path} line {n}: expected KEY = value or [section], got: {s}")
        val, com = m.group(2), ""
        if "#" in val:
            val, com = val.split("#", 1)
        key = f"{section}.{m.group(1)}" if section else m.group(1)
        out.append((key, val.strip().strip('"').strip("'"), com.strip()))
    return out


def exports(conf):
    """section keys become CFG_<section>_<key>, top-level keys keep their name"""
    lines = []
    for k, v, _ in conf:
        v = v.replace("\\", "\\\\").replace('"', '\\"').replace("`", "\\`")
        k = "CFG_" + k.replace(".", "_") if "." in k else k
        lines.append(f'export {k}="{v}"')
    return "\n".join(lines)


def split_args(text):
    """'-a 1 --flag -b x' -> [('-a', '1'), ('--flag', None), ('-b', 'x')]. An option
    takes the next token as its value unless that token is itself an option."""
    opt = lambda t: t.startswith("-") and len(t) > 1 and not re.match(r"-[\d.]", t)
    toks, out, i = text.split(), [], 0
    while i < len(toks):
        t = toks[i]
        if not opt(t):
            out.append((t, None))               # a bare word: kept as it is
            i += 1
            continue
        if i + 1 < len(toks) and not opt(toks[i + 1]):
            out.append((t, toks[i + 1]))
            i += 2
        else:
            out.append((t, None))
            i += 1
    return out


def merge_args(default, changes):
    """the default options with the changes applied: an option given again replaces
    its value (or is added), '-option none' removes it"""
    merged = dict(split_args(default))
    for k, v in split_args(changes):
        if v == "none":
            merged.pop(k, None)
        else:
            merged[k] = v
    return " ".join(k if v is None else f"{k} {v}" for k, v in merged.items())


def _task_mem(v):
    return [(f"{p.split('=', 1)[0]}.mem", p.split("=", 1)[1]) for p in v.split() if "=" in p]


# settings from before the config had sections -> what they are now
RENAMED = {
    "MAX_RUNNING": lambda v: [("cluster.max_running", v)],
    "TRIM_ARGS": lambda v: [("demux.args", v)],
    "BOOT_QUEUE": lambda v: [("boot.queue", v)],
    "CONSTRAINT": lambda v: [("boot.slurm", f"constraint={v}"), ("moments.slurm", f"constraint={v}")] if v else [],
    "TASK_MEM": _task_mem,
    "MINQ5": lambda v: [("MINQ", v)],
    "CLONE_MINQ5": lambda v: [("CLONE_MINQ", v)],
    "CLONEH": lambda v: [("CLONE_CUT", v)] if v else [],
    "PLOT_KRAMP": lambda v: [("PLOT_KCOLOR", v)] if v else [],
    "HET_MINCOV": lambda v: [],              # shared-site heterozygosity, replaced by HET_DEPTH
}


def render_init(template, name):
    """the template as a config to fill in: each default 'args = ...' becomes a
    commented default line and an empty 'args =' for changes"""
    out = []
    for line in template.splitlines():
        m = re.match(r"args\s*=\s*([^#]*?)\s*(#.*)?$", line)
        if m and m.group(1):
            out.append(f"# default args: {m.group(1)}")
            out.append(f"args =                   # changes only: an option replaces the default's value, "
                       f"'-option none' removes it (comrad defaults {name})")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


# sections renamed since: old name -> new
SECTIONS_RENAMED = {"trim": "demux"}


def find_cluster(wd):
    """cluster.config in this folder or up to three folders above it"""
    d = wd
    for _ in range(4):
        f = os.path.join(d, "cluster.config")
        if os.path.exists(f):
            return f
        d = os.path.dirname(d)
    return None


class Pipeline:
    def __init__(self, name, wd):
        self.name = name
        self.wd = os.path.abspath(wd)
        self.src = os.path.join(PIPES, name)
        if not os.path.isdir(self.src):
            die(f"no pipeline '{name}'; there are: {', '.join(sorted(os.listdir(PIPES)))}")
        self.steps = [l.split("#")[0].strip() for l in open(os.path.join(self.src, "STEPS"))]
        self.steps = [s for s in self.steps if s]
        self.conf_path = os.path.join(self.wd, f"{name}.config")
        old = os.path.join(self.wd, f"{name}.conf")
        if not os.path.exists(self.conf_path) and os.path.exists(old):
            print(f"{prog}: reading {name}.conf; rename it {name}.config", file=sys.stderr)
            self.conf_path = old
        self.dir = os.path.join(self.wd, "run", name)
        self.summary = os.path.join(self.dir, f"{name}Summary.txt")
        self.job = f"{os.path.basename(self.wd)}.{name}"

    # ---- config
    def conf(self):
        if not os.path.exists(self.conf_path):
            die(f"no {self.conf_path}; write one with: {prog} init {self.name}")
        return parse_conf(self.conf_path)

    def template(self):
        """this pipeline's defaults, then the [cluster] defaults"""
        return (parse_conf(os.path.join(self.src, "template.config")) +
                parse_conf(os.path.join(PIPES, "cluster.config")))

    def cluster_file(self):
        return find_cluster(self.wd)

    def cluster(self):
        """cluster.config entries, as cluster.<key> (its keys may be top-level or under [cluster])"""
        f = self.cluster_file()
        if not f:
            return []
        return [(k if k.startswith("cluster.") else "cluster." + k.split(".")[-1], v, c)
                for k, v, c in parse_conf(f, section="cluster")]

    def own(self):
        """this config's lines, with settings from before sections renamed to theirs"""
        out = []
        for k, v, c in self.conf():
            sec = k.split(".", 1)[0] if "." in k else ""
            if sec in SECTIONS_RENAMED:
                k = SECTIONS_RENAMED[sec] + k[len(sec):]
            if k in RENAMED:
                for nk in RENAMED[k](v):
                    out.append((nk[0], nk[1], c))
            else:
                out.append((k, v, c))
        return out

    def settings(self):
        """template defaults, then cluster.config, then this config: later lines win.
        An empty [cluster] line in the config leaves cluster.config's value."""
        merged = {}
        for src, rows in (("t", self.template()), ("c", self.cluster()), ("o", self.own())):
            for k, v, c in rows:
                if src == "o" and k.startswith("cluster.") and not v and k in merged:
                    continue
                if src == "o" and k.endswith(".args") and k in merged:
                    v = merge_args(merged[k][1], v)      # args are changes to the defaults
                merged[k] = (k, v, c)
        return list(merged.values())

    def setting(self, key, default=""):
        return next((v for k, v, _ in self.settings() if k == key), "") or default

    def missing(self):
        """defaults in use because neither the config nor cluster.config sets them"""
        have = {k for k, _, _ in self.own()} | {k for k, _, _ in self.cluster()}
        return [(k, v, c) for k, v, c in self.template() if k not in have]

    def check(self, quiet=False):
        conf = dict((k, v) for k, v, _ in self.settings())
        bad = []
        for k, v, com in self.template():
            tags = com.split(":")[0].split()
            val = conf.get(k, "")
            if "required" in tags and not val:
                bad.append(f"{k} is required")
        if bad:
            for b in bad:
                print(f"  {b}")
            die(f"{self.conf_path} is not ready")
        known = {k for k, _, _ in self.template()}
        base = os.path.basename(self.conf_path)
        for k, v, _ in self.conf():
            sec = k.split(".", 1)[0] if "." in k else ""
            if sec in SECTIONS_RENAMED:
                print(f"  [{sec}] is now [{SECTIONS_RENAMED[sec]}]; read as {SECTIONS_RENAMED[sec]}{k[len(sec):]}")
                continue
            if k in RENAMED:
                new = ", ".join((f"[{a.split('.')[0]}] {a.split('.')[1]} = {b}" if "." in a else f"{a} = {b}") for a, b in RENAMED[k](v)) or "nothing (empty)"
                print(f"  {k} is from an older config: read as {new}; rename or move it")
            elif "." in k and k.split(".", 1)[1] in known:
                sec, key = k.split(".", 1)
                print(f"  {key} is under [{sec}] in {base}, so it is read as {k} and ignored; "
                      f"move it above the first [section] line")
            elif k not in known:
                print(f"  {k} in {base} is not a setting of this pipeline; ignored")
        for line in self.skip_notes():
            print(f"  {line}")
        mine = {k for k, v, _ in self.own() if k.endswith(".args") and v}
        args = [(k, v) for k, v, _ in self.settings() if k.endswith(".args")]
        if args:
            print("  options each step will run with (* = changed in your config):")
            for k, v in args:
                print(f"  {'*' if k in mine else ' '} [{k.split('.')[0]}] {v}")
        if not quiet:
            self.report_defaults(base)
        r = self.bash("check", open(os.path.join(self.src, "check.sh")).read(), log=None)
        if r != 0:
            die("check failed; see above")
        if not quiet:
            print(f"{self.conf_path}: ok")

    def skippable(self):
        """step -> (level, what changes), from the pipeline's SKIPPABLE file"""
        f = os.path.join(self.src, "SKIPPABLE")
        if not os.path.exists(f):
            return {}
        out = {}
        for l in open(f):
            if l.strip() and not l.startswith("#"):
                s, lev, msg = l.rstrip("\n").split("\t", 2)
                out[s] = (lev, msg)
        return out

    def skip_notes(self):
        """WARNING / note lines for the steps in SKIP; dies on a step that cannot be skipped"""
        want = self.setting("SKIP").split()
        if not want:
            return []
        ok = self.skippable()
        bad = [s for s in want if s not in self.steps]
        if bad:
            die(f"SKIP: no step {', '.join(bad)} in {self.name} (see: {prog} list {self.name})")
        need = [s for s in want if s not in ok]
        if need:
            die(f"SKIP: {', '.join(need)} cannot be skipped: later steps need {'its' if len(need) == 1 else 'their'} output. "
                f"Steps that can be: {' '.join(ok) or 'none'}")
        return [f"{'WARNING' if ok[s][0] == 'warn' else 'note'}: SKIP {s}: {ok[s][1]}" for s in want]

    def report_defaults(self, base):
        """which defaults apply because nothing sets them (grouped by section)"""
        top, secs = [], {}
        for k, v, _ in self.missing():
            if "." in k:
                sec, key = k.split(".", 1)
                if sec != "cluster":
                    secs.setdefault(sec, []).append(f"{key}={v or '(empty)'}")
            else:
                top.append(f"{k} = '{v}'")
        for t in top:
            print(f"  no {t.split(' =')[0]} line in {base}: using the default{t[len(t.split(' =')[0]) + 2:]}")
        for sec, kv in secs.items():
            print(f"  [{sec}] defaults: {' '.join(kv)}")
        cl = self.cluster_file()
        print(f"  cluster settings: {cl if cl else f'no cluster.config found (write one with: {prog} init cluster)'}")
        for k in ("queue", "long_queue", "max_mem", "max_running", "email", "slurm", "driver_queue", "driver_time", "local_copy"):
            print(f"    {k} = {os.path.expandvars(self.setting('cluster.' + k)) or '(empty)'}")

    # ---- state
    def path(self, step, ext):
        return os.path.join(self.dir, f"{step}.{ext}")

    def done(self, step):
        return os.path.exists(self.path(step, "done"))

    def log(self, msg):
        os.makedirs(self.dir, exist_ok=True)
        line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{msg}"
        with open(self.summary, "a") as fh:
            fh.write(line + "\n")
        print(line, flush=True)

    # ---- running
    def script(self, step, body):
        parts = ["#!/bin/bash",
                 f"export WD={self.wd}", f"export PIPELINE={self.name}", f"export COMRAD={HERE}",
                 f"export CONFIG={os.path.basename(self.conf_path)}",
                 f"cd {self.wd}",
                 "# ---- config (defaults, then cluster.config, then the config)",
                 exports(self.settings()),
                 "# ---- common", open(os.path.join(PIPES, "common.sh")).read()]
        pre = os.path.join(self.src, "prelude.sh")
        if os.path.exists(pre):
            parts += ["# ---- prelude", open(pre).read()]
        parts += [f"# ---- {step}", body, "exit 0"]
        return "\n".join(parts) + "\n"

    def bash(self, step, body, log):
        os.makedirs(self.dir, exist_ok=True)
        sh = self.path(step, "sh")
        open(sh, "w").write(self.script(step, body))
        if log is None:
            return subprocess.run(["bash", sh], cwd=self.wd).returncode
        with open(log, "a") as fh:
            fh.write(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')} {step}\n")
            fh.flush()
            p = subprocess.Popen(["bash", sh], cwd=self.wd, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True)
            for line in p.stdout:
                fh.write(line)
                fh.flush()
                if self.echo:
                    print(line, end="", flush=True)
            return p.wait()

    def since(self, step):
        """When this step's inputs last changed: the latest finish among the steps
        before it, or when it was rewound with `from`."""
        i = self.steps.index(step)
        times = [os.path.getmtime(self.path(s, "done")) for s in self.steps[:i] if self.done(s)]
        if os.path.exists(self.path(step, "fresh")):
            times.append(os.path.getmtime(self.path(step, "fresh")))
        return max(times) if times else 0.0

    def run_step(self, step):
        body = open(os.path.join(self.src, f"{step}.sh")).read()
        os.environ["RT_SINCE"] = f"{self.since(step):.6f}"
        os.environ["RT_STEP"] = step
        logf = self.path(step, "log")
        t0 = time.time()
        self.log(f"run\t{step}")
        rc = self.bash(step, body, logf)
        out = open(logf).read().rsplit("\n===== ", 1)[-1]
        for l in out.splitlines():
            if l.startswith("NOTE: "):
                self.log(f"note\t{step}: {l[6:]}")
        last = lambda tag: next((l[len(tag):] for l in reversed(out.splitlines()) if l.startswith(tag)), None)
        if rc in (0, SKIP) and os.path.exists(self.path(step, "fresh")):
            os.remove(self.path(step, "fresh"))
        if rc == 0:
            open(self.path(step, "done"), "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
            self.log(f"done\t{step} ({(time.time() - t0) / 3600:.2f} h)")
        elif rc == SKIP:
            open(self.path(step, "done"), "w").write("skipped\n")
            self.log(f"skip\t{step}: {last('SKIP: ') or ''}")
            rc = 0
        elif rc == PAUSE:
            self.log(f"PAUSE\t{step}: {last('PAUSE: ') or 'paused'}; then: {prog} run {self.name}")
        else:
            self.log(f"STOP\t{step}: {last('STOP: ') or f'exit {rc}'} (log: {logf})")
        return rc

    def drive(self):
        self.log(f"start\tdriver job {os.environ.get('SLURM_JOB_ID', '-')}")
        self.provenance()
        for line in self.skip_notes():
            self.log(f"note\t{line.removeprefix('note: ')}")
        for step in self.steps:
            if self.done(step):
                continue
            rc = self.run_step(step)
            if rc == PAUSE:
                return 0
            if rc != 0:
                return 1
        self.log("FINISHED")
        return 0

    def provenance(self):
        """tool versions and the conf into the diary; never stops a run"""
        try:
            provenance.record(self, parse_conf)
        except Exception as err:                                     # noqa: BLE001
            self.log(f"note\tversions not recorded: {err}")

    # ---- commands
    def listing(self):
        nxt = None
        for s in self.steps:
            st = "done" if self.done(s) else "pending"
            if st == "done" and open(self.path(s, "done")).read().startswith("skipped"):
                st = "skipped"
            mark = ""
            if nxt is None and st == "pending":
                nxt, mark = s, "  <- next"
            print(f"{st:<9} {s}{mark}")

    def status(self):
        if os.path.exists(self.summary):
            ev = [l for l in open(self.summary) if l.split("\t")[1:2] not in (["version"], ["conf"])]
            print("".join(ev[-10:]), end="")
        else:
            print(f"nothing run yet ({self.summary})")
        print(f"driver in the queue: {self.queued() or 'no'}")
        # array tasks of yours in the queue, by job name
        r = subprocess.run(["squeue", "-h", "-r", "-u", os.environ.get("USER", ""), "-o", "%j %T"],
                           capture_output=True, text=True)
        tally = {}
        for line in r.stdout.splitlines():
            f = line.split()
            if len(f) == 2 and f[0] != self.job:
                tally.setdefault(f[0], {}).setdefault(f[1].lower(), 0)
                tally[f[0]][f[1].lower()] += 1
        for n, c in sorted(tally.items()):
            print(f"  {n}: " + ", ".join(f"{v} {k}" for k, v in sorted(c.items())))

    def queued(self):
        return subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-n", self.job,
                               "-o", "%i %T %M"], capture_output=True, text=True).stdout.strip()

    def rewind(self, step, earlier=False):
        if step not in self.steps:
            die(f"no step '{step}' in {self.name}; steps: {' '.join(self.steps)}")
        i = self.steps.index(step)
        os.makedirs(self.dir, exist_ok=True)
        for s in self.steps[i:]:
            if os.path.exists(self.path(s, "done")):
                os.remove(self.path(s, "done"))
            open(self.path(s, "fresh"), "w").write("rewound: rerun every job\n")
        if earlier:
            for s in self.steps[:i]:
                if not self.done(s):
                    open(self.path(s, "done"), "w").write("marked done by --mark-earlier\n")
        self.log(f"rewind\tfrom {step}" + ("; earlier steps marked done" if earlier else ""))

    def submit(self, a):
        q = self.queued()
        if q:
            die(f"{self.job} is already queued or running: {q}")
        os.makedirs(self.dir, exist_ok=True)
        me = shutil.which("comrad") or os.path.abspath(sys.argv[0])
        st = lambda k: os.path.expandvars(self.setting("cluster." + k))
        part = a.partition or st("driver_queue") or st("long_queue") or st("queue")
        cmd = ["sbatch"] + (["-p", part] if part else []) + \
              ["-N", "1", "-n", "1", "-c", str(a.cpus or st("driver_cpus") or 8),
               f"--mem={a.mem or st('driver_mem') or '32GB'}"] + \
              (["-t", a.time or st("driver_time")] if (a.time or st("driver_time")) else []) + [
               "-J", self.job, "-o", os.path.join(self.dir, "driver.o%j"),
               "-e", os.path.join(self.dir, "driver.e%j")]
        if st("email"):
            cmd += ["--mail-type=END,FAIL", f"--mail-user={st('email')}"]
        cmd += rtjobs.slurm_flags(st("slurm"))
        cmd += ["--wrap", f"{me} _drive {self.name} --dir {self.wd}"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        print((r.stdout + r.stderr).strip())
        if r.returncode == 0:
            self.log(f"submit\t{r.stdout.strip()} ({part or 'default partition'})")
        return r.returncode


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  cd /path/to/myProject
  {prog} init cluster               # writes cluster.config (partitions, memory, email); once per project
  {prog} init reads                 # writes reads.config; fill it in
  {prog} check reads                # is everything it names there?
  {prog} run reads                  # submit; species folders appear when it finishes

  cd /path/to/myProject/sint
  {prog} init popgen && nano popgen.config
  {prog} run popgen                 # pauses at K
  {prog} status popgen              # where it is, or why it stopped
  nano popgen.config && {prog} run popgen   # K = 2 in the config
  {prog} from popgen clonegroups    # forget that step and everything after, resubmit
  {prog} from popgen kin --mark-earlier --no-submit   # earlier steps were run by hand
  {prog} step popgen fst            # run one step now, in this shell
  {prog} list demography            # the steps, done or pending
  {prog} report popgen              # run/popgen/popgenReport.html: warnings first, then the key tables
  {prog} summary popgen --dir mcav ofav sint past   # from the project folder: what each step said, and what is waiting on you
  {prog} setup --prefix /path/to/bin/comradDeps   # find or install the external tools
""")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="command")
    s = sub.add_parser("setup", help="find or install the tools the pipelines run",
                       description=setup_tools.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                       epilog=f"""examples:
  {prog} setup --prefix /path/to/bin/comradDeps      # everything missing
  {prog} setup --prefix /path/to/bin/comradDeps --only evalAdmix NGSremix
  {prog} setup --prefix /path/to/bin/comradDeps --no-conda
  echo 'source /path/to/bin/comradDeps/comrad.env' >> ~/.bashrc
""")
    s.add_argument("--prefix", help="folder to install into (bin/, src/, Rlib/, conda/ go there)")
    s.add_argument("--only", nargs="+", metavar="TOOL", help="just these tools (names as in docs/install.md)")
    s.add_argument("--no-modules", action="store_true", help="do not look for environment modules")
    s.add_argument("--no-conda", action="store_true", help="do not use conda/mamba, and do not install Miniforge")
    s.add_argument("--optional", action="store_true", help="also try optional tools (BaseSpace CLI)")
    s.add_argument("--write-docs", metavar="FILE", help=argparse.SUPPRESS)
    s = sub.add_parser("defaults", help="print every default setting of a pipeline",
                       epilog=f"examples:\n  {prog} defaults popgen\n  {prog} defaults popgen | grep -A4 '\\[clones\\]'\n  {prog} defaults cluster",
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    s.add_argument("pipeline", help="reads, popgen, demography or cluster")
    for c, h in [("init", "write <pipeline>.config (or cluster.config) to fill in"), ("check", "check <pipeline>.config and its inputs"),
                 ("run", "submit the driver job"), ("status", "last events, and whether the driver is queued"),
                 ("list", "every step, done or pending"), ("from", "clear the done markers from STEP on, then submit"),
                 ("step", "run one step now, in this shell"),
                 ("report", "write run/<pipeline>/<pipeline>Report.html: warnings, pauses and stops, then the key tables"),
                 ("summary", "print open pauses and stops, then each step's notes and warnings from its last run"),
                 ("_drive", argparse.SUPPRESS)]:
        s = sub.add_parser(c, help=h) if h != argparse.SUPPRESS else sub.add_parser(c)
        s.add_argument("pipeline", help="reads, popgen or demography" + (" (or cluster: write cluster.config)" if c == "init" else ""))
        if c == "summary":
            s.add_argument("--dir", nargs="+", default=["."], help="folder(s) the pipeline works in, e.g. one per species (default: here)")
        else:
            s.add_argument("--dir", default=".", help="folder the pipeline works in (default: here)")
        if c in ("from", "step"):
            s.add_argument("stepname", metavar="STEP", help="step name, as `list` prints it")
        if c == "from":
            s.add_argument("--mark-earlier", action="store_true",
                           help="mark every step before STEP done (they were run by hand)")
            s.add_argument("--no-submit", action="store_true", help="clear the markers only")
        if c == "init":
            s.add_argument("--force", action="store_true", help="overwrite an existing .config")
        if c in ("run", "from"):
            s.add_argument("-p", "--partition", help="driver partition (default: [cluster] driver_queue, else long_queue, else queue)")
            s.add_argument("-c", "--cpus", type=int, help="driver CPUs (default: [cluster] driver_cpus)")
            s.add_argument("--mem", help="driver memory (default: [cluster] driver_mem)")
            s.add_argument("-t", "--time", help="driver time limit, e.g. 2-00:00:00 (default: [cluster] driver_time, else none)")
            s.add_argument("--here", action="store_true", help="run the steps in this shell instead of submitting")
    a = p.parse_args()

    if a.cmd == "setup":
        sys.exit(setup_tools.run(a))
    if a.cmd == "defaults":
        src = os.path.join(PIPES, "cluster.config" if a.pipeline == "cluster" else os.path.join(a.pipeline, "template.config"))
        if not os.path.exists(src):
            die(f"no pipeline '{a.pipeline}' (reads, popgen, demography or cluster)")
        print(open(src).read(), end="")
        return
    if a.cmd == "init" and a.pipeline == "cluster":
        dst = os.path.join(os.path.abspath(a.dir), "cluster.config")
        if os.path.exists(dst) and not a.force:
            die(f"{dst} exists; --force to overwrite")
        shutil.copy(os.path.join(PIPES, "cluster.config"), dst)
        print(f"wrote {dst}; set your partitions, max_mem and email. Every pipeline run in this folder "
              f"or up to three folders below it uses it.")
        return
    if a.cmd == "summary":
        for i, d in enumerate(a.dir):
            try:
                print(("\n" if i else "") + report.summary(Pipeline(a.pipeline, d)), flush=True)
            except SystemExit as e:
                print(("\n" if i else "") + f"== {d}: {e.code}" if isinstance(e.code, str) else "", flush=True)
        return
    pl = Pipeline(a.pipeline, a.dir)
    pl.echo = a.cmd in ("step",) or getattr(a, "here", False)
    if a.cmd == "init":
        if os.path.exists(pl.conf_path) and not a.force:
            die(f"{pl.conf_path} exists; --force to overwrite")
        cl = pl.cluster_file()
        txt = render_init(open(os.path.join(pl.src, "template.config")).read(), a.pipeline)
        blk = [l for l in open(os.path.join(PIPES, "cluster.config")).read().splitlines() if l.startswith("[") or "=" in l.split("#")[0]]
        txt += ("\n# ---- cluster: these come from " + (cl if cl else f"cluster.config (none found yet: {prog} init cluster in the project folder)") +
                ".\n# To override one for this pipeline only, uncomment [cluster] and that line.\n" + "\n".join("# " + l for l in blk) + "\n")
        open(pl.conf_path, "w").write(txt)
        print(f"wrote {pl.conf_path}; fill in the required lines, then: {prog} check {a.pipeline}")
        if not cl:
            print(f"no cluster.config found here or above; write one in the project folder with: {prog} init cluster")
    elif a.cmd == "check":
        pl.check()
    elif a.cmd == "list":
        pl.listing()
    elif a.cmd == "status":
        pl.status()
    elif a.cmd == "report":
        print(report.write(pl))
    elif a.cmd == "_drive":
        sys.exit(pl.drive())
    elif a.cmd == "step":
        pl.check(quiet=True)
        pl.provenance()
        rc = pl.run_step(a.stepname) if a.stepname in pl.steps else die(f"no step '{a.stepname}'")
        sys.exit(0 if rc in (0, PAUSE) else 1)
    else:
        if a.cmd == "from":
            pl.rewind(a.stepname, a.mark_earlier)
            if a.no_submit:
                return
        pl.check(quiet=True)
        if a.here:
            sys.exit(pl.drive())
        sys.exit(pl.submit(a))


if __name__ == "__main__":
    main()
