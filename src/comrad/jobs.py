#!/usr/bin/env python3
################################################################################
#Script Name    : comrad-jobs
#Description    : Submit a job list, check every job's output, rerun only the failed ones
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Submit a list of commands to SLURM, wait, check each one, rerun the failures.

Each line of JOBS is one command. A line's check decides whether it is done:

  check<TAB>command      explicit: the check is the text before the first tab
  command                with -c TEMPLATE: the check is TEMPLATE, with {1}, {2}
                         ... filled from the -x regexes matched against the line

A check is a bash test (`test -s out.saf.idx`, `grep -q Finished! x.log`). A job
is done when its check passes and its command exited 0 (without a check, exit 0
alone). Lines already done are not submitted again, so a rerun only does what is
missing.

Under comrad, "already done" also needs the job to have succeeded after
RT_SINCE (when the steps before it last finished, or it was rewound): every
success leaves a stamp in .rtjobs/<name>/, named after the command, with its SLURM
job and log. .rtjobs/<name>/logs lists the log behind every job in the list
(`rtlogs NAME` in the pipelines).

How jobs are submitted:
  one line    a plain sbatch job, logs <name>.o<jobid> / <name>.e<jobid>
  more lines  a SLURM job array, one task per line, at most --max-running at once
              (default: cluster max_running, else 64), split into arrays of <= 1000. Each task
              writes its own header (node, start, end, exit code) into its log, in
              .rtjobs/<name>/tasks.<...>/; when the array ends they are combined
              into <name>.e<arrayid> / .o<arrayid> (a summary table, failed tasks
              first, then one section per task) and the task folder is removed.
              It is kept if combining fails.
--cpus and --mem are per job (per task for arrays). No time limit is asked for
unless -t is given: the partition's own default applies. No partition is asked
for unless -q is given (or the config sets one): the cluster's default applies.

Under comrad, --res SECTION takes cpus, mem, queue and slurm from that
section of the pipeline config, for any of them not given on the command line,
and [cluster] gives the defaults: queue, long_queue, max_mem, max_running,
email and slurm (extra sbatch flags for every job).

After each round, only the failed jobs are resubmitted, by what went wrong:
  time limit     -> on --long (cluster long_queue), if one is set; otherwise stop
  out of memory  -> that job's --mem doubled, up to --max-mem (cluster max_mem, default 250GB)
  node failure   -> again, avoiding the node it ran on (does not use up a retry)
  anything else  -> once more (--retries), then stop
Exits 1, listing the lines still failing, if any are left. Before submitting, it
waits for jobs of the same name still in the queue (from an interrupted run).

Under comrad, a diary note per job list gives the tasks' peak memory (from sacct, if it answers), longest run and
CPU-hours.
"""
import argparse
import concurrent.futures as cf
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time

prog = os.path.basename(sys.argv[0])
POLL = int(os.environ.get("RTJOBS_POLL", 120))
TIMEOUT = re.compile(r"DUE TO TIME LIMIT|TIMEOUT", re.I)
OOM = re.compile(r"oom[-_ ]kill|out of memory|OUT_OF_MEMORY|Killed process|exceeded memory limit", re.I)


def say(msg):
    print(f"{prog}: {msg}", flush=True)


def read_jobs(path, template, regexes):
    jobs = []
    for n, line in enumerate(open(path), 1):
        line = line.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if "\t" in line and not template:
            check, cmd = line.split("\t", 1)
        else:
            cmd, check = line, None
            if template:
                vals = []
                for rx in regexes:
                    m = re.search(rx, cmd)
                    if not m:
                        sys.exit(f"{prog}: line {n} of {path} does not match -x {rx!r}:\n  {cmd}")
                    vals += list(m.groups()) or [m.group(0)]
                check = template
                for i, v in enumerate(vals, 1):
                    check = check.replace("{%d}" % i, v)
        jobs.append({"cmd": cmd, "check": check, "done": False})
    return jobs


def passes(check):
    return subprocess.run(["bash", "-c", check], stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode == 0


# A job is only done if it also succeeded after its inputs last changed. comrad
# sets RT_SINCE to the time the steps before this one last finished (or this step
# was rewound with `from`); each success leaves a stamp named after the command, so
# old outputs, outputs from before an upstream rerun, and changed commands rerun.
SINCE = float(os.environ["RT_SINCE"]) if os.environ.get("RT_SINCE") else None


def stamp(name, cmd):
    return os.path.join(".rtjobs", re.sub(r"\W", "_", name), hashlib.sha1(cmd.encode()).hexdigest())


def fresh(name, cmd):
    if SINCE is None:
        return True
    s = stamp(name, cmd)
    return os.path.exists(s) and os.path.getmtime(s) > SINCE


def mark(name, cmd, job="", log=""):
    s = stamp(name, cmd)
    os.makedirs(os.path.dirname(s), exist_ok=True)
    open(s, "w").write(cmd + "\n" + (f"job\t{job}\nlog\t{os.path.abspath(log)}\n" if log else ""))


def stamp_log(name, cmd):
    """the .e log of the job that made this command's current output, if the stamp has one"""
    s = stamp(name, cmd)
    if os.path.exists(s):
        for line in open(s):
            if line.startswith("log\t"):
                return line.rstrip("\n").split("\t", 1)[1]
    return None


def done_now(j, name):
    if j["check"]:
        # a job that never got submitted cannot have made its output, even if an
        # old copy of that output passes the check
        ok = bool(j.get("jobid")) and passes(j["check"])
        # a plain sbatch job's exit code is the command's own (a launcher's is not),
        # so a failed command is not done even if an older output passes the check
        if j.get("plain") and j.get("rc") != 0:
            ok = False
    else:
        ok = j.get("rc") == 0
    if ok:
        mark(name, j["cmd"], j.get("jobid", ""), j.get("logfile", ""))
    return ok


def write_logs(name, jobs):
    """.rtjobs/<name>/logs: one .e log per line, for the jobs done now"""
    d = os.path.dirname(stamp(name, ""))
    os.makedirs(d, exist_ok=True)
    out, missing = [], 0
    for j in jobs:
        if not j["done"]:
            continue
        lg = stamp_log(name, j["cmd"])
        if lg and os.path.exists(lg):
            if lg not in out:
                out.append(lg)
        else:
            missing += 1
    with open(os.path.join(d, "logs"), "w") as fh:
        fh.write("".join(l + "\n" for l in out))
    if missing:
        say(f"{missing} done jobs have no recorded log (stamped before logs were recorded, "
            f"or the log is gone); steps fall back to {name}*.e* for them")
        open(os.path.join(d, "logs"), "a").write("FALLBACK\n")


def jobid(text):
    m = re.search(r"Submitted batch job (\d+)", text)
    return m.group(1) if m else ""


def queued(name):
    out = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%i %j"],
                         capture_output=True, text=True).stdout.split("\n")
    return [l.split()[0] for l in out if l.split()[1:] and
            (l.split()[1] == name or l.split()[1].startswith(name + "."))]


def wait_queue(name):
    said = 0
    while True:
        q = queued(name)
        if not q:
            return
        if time.time() - said > 3600:
            say(f"{len(q)} '{name}' jobs from an earlier run still queued; waiting")
            said = time.time()
        time.sleep(POLL)


def newest(pattern):
    f = glob.glob(pattern)
    return max(f, key=os.path.getmtime) if f else None


def text(path):
    try:
        return open(path, errors="replace").read()
    except OSError:
        return ""


NODEFAIL = re.compile(r"NODE_FAIL|PREEMPTED|BOOT_FAIL|launch failed|Unable to start", re.I)
GB = {"": 1e-9, "K": 1e-6, "M": 1e-3, "G": 1, "T": 1e3}
ARRAY_MAX = 1000


def gb(m):
    """'200GB' / '8G' / '512M' -> GB as a float; None when unset"""
    if not m:
        return None
    v = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([KMGT]?)B?\s*", str(m), re.I)
    return float(v.group(1)) * GB[v.group(2).upper()] if v else None


def mem_str(g):
    return f"{int(round(g))}GB"


def cfg(section, key):
    """a config setting exported by comrad as CFG_<section>_<key>; '' when unset"""
    return os.environ.get(f"CFG_{section}_{key}", "").strip()


def slurm_flags(text):
    """'constraint=amd qos=fast exclusive A=myacct' -> sbatch flags"""
    out = []
    for t in (text or "").split():
        if t.startswith("-"):
            out.append(t)
            continue
        k, eq, v = t.partition("=")
        if len(k) == 1:
            out += [f"-{k}", v] if eq else [f"-{k}"]
        else:
            out.append(f"--{k}={v}" if eq else f"--{k}")
    return out


def submit_wait(sb):
    r = subprocess.run(sb, capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    if out:
        print(out, flush=True)
    return r.returncode, jobid(r.stdout)


def sbatch_base(a, name, queue, cpus, mem, exclude):
    sb = ["sbatch", "--wait"] + (["-p", queue] if queue else []) + \
         ["-N", "1", "-n", "1", "-c", str(cpus), "-J", name]
    if a.time and queue == a.queue:
        sb += ["-t", a.time]
    elif a.long_time and queue == a.long:
        sb += ["-t", a.long_time]
    if mem:
        sb.append(f"--mem={mem}")
    if a.constraint:
        sb.append(f"--constraint={a.constraint}")
    if exclude:
        sb.append("--exclude=" + ",".join(sorted(exclude)))
    if a.email:
        sb += [f"--mail-user={a.email}", "--mail-type=FAIL"]
    return sb + a.slurm_flags


WRAP = ('echo "##rt start $(hostname) $(date +%s) ${{SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}}" >&2; '
        '{cmd}; rc=$?; echo "##rt end $rc $(date +%s)" >&2; exit $rc')


def parse_task(etext):
    """node, start, end, rc from the header and footer the task wrote itself"""
    node = start = end = rc = None
    m = re.search(r"^##rt start (\S+) (\d+)", etext, re.M)
    if m:
        node, start = m.group(1), int(m.group(2))
    m = re.search(r"^##rt end (-?\d+) (\d+)", etext, re.M)
    if m:
        rc, end = int(m.group(1)), int(m.group(2))
    return node, start, end, rc


def classify(j):
    """why a job failed: timeout, oom, node, other"""
    log = j.get("log", "")
    if TIMEOUT.search(log):
        return "timeout"
    if OOM.search(log):
        return "oom"
    if NODEFAIL.search(log) or NODEFAIL.search(j.get("state", "")) or (j.get("node") and not j.get("ended")) \
            or not j.get("node"):
        return "node"
    return "other"


def sacct_stats(ids):
    """{taskid: (state, maxrss GB)} from sacct; empty if it does not answer"""
    out = {}
    try:
        r = subprocess.run(["sacct", "-n", "-P", "-j", ",".join(ids), "-o", "JobID,State,MaxRSS"],
                           capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return out
    if r.returncode:
        return out
    for line in r.stdout.splitlines():
        f = line.split("|")
        if len(f) < 3:
            continue
        base = f[0].split(".")[0]
        st, rss = out.get(base, ("", 0.0))
        if "." not in f[0]:
            st = f[1]
        v = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)", f[2] or "")
        if v:
            rss = max(rss, float(v.group(1)) * GB[v.group(2)])
        out[base] = (st, rss)
    return out


def single_round(a, name, j, queue, mem, exclude):
    """one job: plain sbatch, so the exit code is the command's own"""
    sb = sbatch_base(a, name, queue, a.cpus, mem, exclude) + ["-o", name + ".o%j", "-e", name + ".e%j"]
    if a.exclusive:
        sb.append("--exclusive")
    sb += ["--wrap", WRAP.format(cmd=j["cmd"])]
    say(f"1 job as {name} ({queue or 'default partition'}, {a.cpus} CPU" + (f", --mem={mem}" if mem else "") +
        (", whole node" if a.exclusive else "") + ")")
    _, jid = submit_wait(sb)
    ef = f"{name}.e{jid}" if jid and os.path.exists(f"{name}.e{jid}") else (newest(f"{name}.e*") or "")
    j["log"] = text(ef) + text(newest(f"{name}.o*") or "")
    j["node"], st, en, j["rc"] = parse_task(j["log"])
    j["ended"] = j["rc"] is not None
    if j["rc"] is None:
        j["rc"] = -1
    j["jobid"], j["logfile"], j["plain"], j["mem"] = jid, ef, True, mem
    if jid:
        s = sacct_stats([jid]).get(jid)
        if s:
            j["state"], j["peak"] = s
    if st and en:
        j["secs"] = en - st
    j["cpus"] = a.cpus


def array_round(a, name, lines, queue, mem, exclude):
    """many jobs: one task each, in arrays of <= ARRAY_MAX; logs combined afterwards"""
    safe = re.sub(r"\W", "_", a.name)
    tmp = os.path.join(".rtjobs", safe, "tasks.%s.%d" % (re.sub(r"\W", "_", name), int(time.time())))
    chunks = [lines[i:i + ARRAY_MAX] for i in range(0, len(lines), ARRAY_MAX)]
    cap = max(1, a.max_running // len(chunks))
    say(f"{len(lines)} jobs as {name} (array{'s' if len(chunks) > 1 else ''}, {queue or 'default partition'}, {a.cpus} CPU" +
        (f", --mem={mem}" if mem else "") + f", at most {cap * len(chunks)} at once)")

    def run_chunk(c, part):
        d = os.path.join(tmp, str(c))
        os.makedirs(d, exist_ok=True)
        cmds = os.path.abspath(os.path.join(d, "cmds"))
        with open(cmds, "w") as fh:
            fh.write("".join(j["cmd"] + "\n" for j in part))
        sb = sbatch_base(a, name, queue, a.cpus, mem, exclude) + [
            f"--array=1-{len(part)}%{cap}", "-o", os.path.join(d, "%a.o"), "-e", os.path.join(d, "%a.e"),
            "--wrap", WRAP.format(cmd=f'bash -c "$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {cmds})"')]
        _, jid = submit_wait(sb)
        return d, jid

    with cf.ThreadPoolExecutor(max_workers=len(chunks)) as ex:
        res = list(ex.map(lambda p: run_chunk(*p), list(enumerate(chunks))))
    stats = sacct_stats([jid for _, jid in res if jid])
    for (d, jid), part in zip(res, chunks):
        for i, j in enumerate(part, 1):
            et, ot = text(os.path.join(d, f"{i}.e")), text(os.path.join(d, f"{i}.o"))
            j["log"], j["out"] = et, ot
            j["node"], st, en, rc = parse_task(et)
            j["ended"] = rc is not None
            j["rc"] = rc if rc is not None else -1
            j["secs"] = (en - st) if st and en else None
            j["jobid"], j["plain"], j["mem"], j["cpus"] = (f"{jid}_{i}" if jid else ""), True, mem, a.cpus
            j["state"], j["peak"] = stats.get(f"{jid}_{i}", ("", None))
            j["tasklog"] = os.path.abspath(os.path.join(d, f"{i}.e"))
    first = next((jid for _, jid in res if jid), str(int(time.time())))
    combine(name, first, lines, tmp)


def combine(name, jid, lines, tmp):
    """<name>.e<jid> / .o<jid>: summary table (failed first) then one section per task"""
    ef, of = f"{name}.e{jid}", f"{name}.o{jid}"
    order = sorted(range(len(lines)), key=lambda k: (lines[k]["rc"] == 0, k))
    def fmt_t(s):
        return "?" if s is None else f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"
    try:
        with open(ef + ".part", "w") as fe, open(of + ".part", "w") as fo:
            fe.write(f"# {name}: {len(lines)} tasks, {sum(1 for j in lines if j['rc'] != 0)} failed\n")
            fe.write("# task\texit\ttime\tpeakGB\tnode\tcommand\n")
            for k in order:
                j = lines[k]
                fe.write(f"# {k + 1}\t{j['rc'] if j['ended'] else 'killed'}\t{fmt_t(j['secs'])}\t"
                         f"{'?' if j.get('peak') is None else '%.1f' % j['peak']}\t"
                         f"{j.get('node') or '?'}\t{j['cmd'][:150]}\n")
            for k, j in enumerate(lines):
                head = (f"===== task {k + 1} of {len(lines)}   job {j['jobid']}   {j.get('node') or '?'}   "
                        f"exit {j['rc'] if j['ended'] else 'killed'}   {fmt_t(j['secs'])}\n{j['cmd']}\n")
                fe.write(head + j.get("log", "") + ("" if j.get("log", "").endswith("\n") else "\n"))
                fo.write(head + j.get("out", "") + ("" if j.get("out", "").endswith("\n") else "\n"))
        body = text(ef + ".part")
        if body.count("===== task ") != len(lines):
            raise OSError(f"{ef}: {body.count('===== task ')} sections for {len(lines)} tasks")
        os.replace(ef + ".part", ef)
        os.replace(of + ".part", of)
    except OSError as err:
        say(f"could not combine the task logs ({err}); they are still in {tmp}")
        for f in (ef + ".part", of + ".part"):
            try:
                os.remove(f)
            except OSError:
                pass
        for j in lines:
            j["logfile"] = j.get("tasklog", "")
        return
    shutil.rmtree(tmp, ignore_errors=True)
    for j in lines:
        j["logfile"] = ef
        j["log"] = j.get("log", "")


def usage_note(name, jobs):
    """NOTE for the diary: jobs, peak memory, longest run, CPU-hours"""
    ran = [j for j in jobs if j.get("jobid")]
    if not ran:
        return
    peaks = [(j["peak"], j.get("mem")) for j in ran if j.get("peak")]
    secs = [j["secs"] for j in ran if j.get("secs") is not None]
    cpuh = sum(j["secs"] * j.get("cpus", 1) for j in ran if j.get("secs") is not None) / 3600
    mem = (f"peak {max(p for p, _ in peaks):.1f} GB of {max(peaks)[1] or 'the default'}" if peaks
           else "peak memory not recorded (sacct did not answer)")
    t = f", longest {max(secs) / 3600:.2f} h, {cpuh:.1f} CPU-hours" if secs else ""
    print(f"NOTE: {name}: {len(ran)} job{'' if len(ran) == 1 else 's'}, {mem}{t}", flush=True)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # angsd per lineage; done when the .saf.idx exists
  {prog} sfsJobs -n sintSFS -c 'test -s {{1}}.saf.idx' -x='-out (\\S+)' --mem 32GB

  # NGSadmix job list from ngsAdmixLauncher.py (-P 10 each), at most 40 at once
  {prog} sintNgsAdmix -n sintNgsAdmix -c 'test -s {{1}}.qopt' -x='-o (\\S+)' --cpus 10 --mem 8GB --max-running 40

  # explicit checks: one 'check<TAB>command' per line
  printf 'test -s a.out\\tsleep 5; echo hi >a.out\\n' >jobs
  {prog} jobs -n demo

  # one NGSremix job: a plain sbatch job with 8 CPUs and 100 GB
  {prog} remixJob -n sintRemix --cpus 8 --mem 100GB

  # bootstrap spectra on AMD nodes, 20 CPUs each, on a named partition
  {prog} bootJobs -n sintBoot -c 'test -s {{1}}' -x='>(\\S+\\.sfs)' --cpus 20 --mem 8GB -q shortq7 --slurm 'constraint=amd'

  # under comrad: cpus, mem, queue and slurm from the config's [boot] section
  {prog} bootJobs -n sintBoot --res boot

  # a regex starting with '-' must be written -x='...' so it is not read as an option
""")
    p.add_argument("jobs", help="job list, one command per line (optionally check<TAB>command)")
    p.add_argument("-n", "--name", required=True, help="job name; logs are named after it")
    p.add_argument("-c", "--check", metavar="TEMPLATE", help="bash test per line, with {1} {2} ... from -x")
    p.add_argument("-x", "--extract", action="append", default=[], metavar="REGEX",
                   help="regex matched against each line; its groups fill {1} {2} ... in order (repeatable)")
    p.add_argument("--res", metavar="SECTION",
                   help="take cpus, mem, queue and slurm from this config section (under comrad)")
    p.add_argument("-q", "--queue", help="partition for the first attempt (default: cluster queue, else the cluster's default)")
    p.add_argument("-t", "--time", default="", help="time limit on --queue; default none (the partition's own)")
    p.add_argument("--long", help="partition after a time-limit failure (default: cluster long_queue; none = stop)")
    p.add_argument("--long-time", default="", help="time limit on --long; default none")
    p.add_argument("--cpus", type=int, help="CPUs per job (per array task) (default: the section's cpus, else 1)")
    p.add_argument("--mem", help="memory per job (per array task), e.g. 8GB (default: the section's mem)")
    p.add_argument("--max-mem", help="out-of-memory retries stop above this (default: cluster max_mem, else 250GB)")
    p.add_argument("--max-running", type=int,
                   help="array tasks running at once (default: cluster max_running, else 64)")
    p.add_argument("--slurm", default="", metavar="FLAGS",
                   help="extra sbatch flags as 'key=value key ...' (added to cluster and section slurm)")
    p.add_argument("--email", help="mail on failure to this address (default: cluster email)")
    p.add_argument("--constraint", help="sbatch --constraint, e.g. amd (same as --slurm constraint=amd)")
    p.add_argument("--exclusive", action="store_true",
                   help="a single job takes the whole node (only for jobs that need it)")
    p.add_argument("--retries", type=int, default=3, help="attempts after the first for anything but node failures (default %(default)s)")
    for old in ("-w", "--wayness", "-N", "--nodes", "--logdir"):
        p.add_argument(old, help=argparse.SUPPRESS)                 # launcher-era options, ignored
    p.add_argument("--each", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    r = a.res or "-"
    q = a.queue or cfg(r, "queue") or cfg("cluster", "queue")
    a.long = a.long or cfg("cluster", "long_queue")
    a.queue = a.long if q == "long" else q
    a.cpus = a.cpus or int(cfg(r, "cpus") or 1)
    a.mem = a.mem or cfg(r, "mem") or None
    a.max_mem = a.max_mem or cfg("cluster", "max_mem") or "250GB"
    a.max_running = a.max_running or int(cfg("cluster", "max_running") or 64)
    a.email = a.email if a.email is not None else cfg("cluster", "email")
    a.slurm_flags = slurm_flags(" ".join([cfg("cluster", "slurm"), cfg(r, "slurm"), a.slurm]))

    if a.extract and not a.check:
        p.error("-x needs -c")
    jobs = read_jobs(a.jobs, a.check, a.extract)
    if not jobs:
        say(f"{a.jobs} is empty; nothing to do")
        return
    wait_queue(a.name)

    base_mem = a.mem
    maxmem = gb(a.max_mem) or 250
    for j in jobs:
        j["done"] = fresh(a.name, j["cmd"]) and (passes(j["check"]) if j["check"] else SINCE is not None)
        j["queue"], j["mem"] = a.queue, base_mem
    stale = sum(1 for j in jobs if not j["done"] and j["check"] and passes(j["check"]))
    if stale:
        say(f"{stale} jobs have outputs from before their inputs last changed; rerunning them")

    exclude, other, nodefail, rnd = set(), 0, 0, 0
    while True:
        todo = [j for j in jobs if not j["done"]]
        if not todo:
            break
        if rnd:
            say(f"attempt {rnd + 1}: {len(todo)} of {len(jobs)} jobs to rerun")
        name = a.name if rnd == 0 else f"{a.name}.retry{rnd}"
        groups = {}
        for j in todo:
            groups.setdefault((j["queue"], j["mem"]), []).append(j)
        for (queue, mem), part in groups.items():
            gname = name if len(groups) == 1 else f"{name}.{len(part)}x{mem or 'mem'}"
            if len(part) == 1:
                single_round(a, gname, part[0], queue, mem, exclude)
            else:
                array_round(a, gname, part, queue, mem, exclude)
        for j in todo:
            j["done"] = done_now(j, a.name)
        failed = [j for j in todo if not j["done"]]
        usage_note(name, todo)
        rnd += 1
        if not failed:
            continue
        why = {}
        for j in failed:
            why.setdefault(classify(j), []).append(j)
        if "timeout" in why:
            if not a.long or all(j["queue"] == a.long for j in why["timeout"]):
                say(f"{len(why['timeout'])} hit the time limit" +
                    ("; set long_queue in [cluster] to move them to a longer partition" if not a.long else
                     f" on {a.long} as well"))
                break
            say(f"{len(why['timeout'])} hit the time limit; next attempt on {a.long}")
            for j in why["timeout"]:
                j["queue"] = a.long
        if "oom" in why:
            for j in why["oom"]:
                cur = gb(j["mem"]) or 4
                j["oom_cap"] = cur >= maxmem
                j["mem"] = mem_str(min(cur * 2, maxmem))
            capped = [j for j in why["oom"] if j.get("oom_cap")]
            say(f"{len(why['oom'])} ran out of memory; next attempt with " +
                ", ".join(sorted({j['mem'] for j in why['oom']})))
            if capped:
                say(f"{len(capped)} already had {a.max_mem} and still ran out; raise max_mem in [cluster] only if a node has more")
                break
        if "node" in why:
            nodes = {j["node"] for j in why["node"] if j.get("node")}
            exclude |= nodes
            nodefail += 1
            say(f"{len(why['node'])} failed on the cluster side" +
                (f" (nodes {', '.join(sorted(nodes))}, avoided from now on)" if nodes else ""))
            if nodefail > 5:
                break
        if "other" in why:
            other += 1
            say(f"{len(why['other'])} failed (exit " + ", ".join(sorted({str(j['rc']) for j in why['other']})) + ")" +
                ("; one more try" if other == 1 else ""))
            if other > 1 or rnd > a.retries + 6:
                break
        if rnd > a.retries + 12:                 # never loop forever
            break
    left = [j for j in jobs if not j["done"]]
    write_logs(a.name, jobs)
    if not left:
        say(f"all {len(jobs)} jobs done")
        return
    say(f"{len(left)} of {len(jobs)} jobs still failing:")
    for j in left[:10]:
        print(f"  {j['cmd'][:200]}")
        if j["check"]:
            print(f"    check: {j['check']}")
    sys.exit(1)


if __name__ == "__main__":
    main()
