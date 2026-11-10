"""comrad setup: find or install the external tools the pipelines run.

For each tool that is not on PATH, in order:
  1. an environment module (module avail): reported, and loaded by comrad.env;
  2. conda / mamba (bioconda, conda-forge), into <prefix>/conda. With no conda or
     mamba on PATH, Miniforge (conda-forge's installer, no admin rights needed) is
     installed into <prefix>/miniforge first; --no-conda skips both;
  3. a download or build from source (GitHub releases, git, make), in <prefix>/src;
  4. pip, for the Python tools.
Every program found or installed is linked into <prefix>/bin, and
<prefix>/comrad.env puts it on PATH.
When nothing works, the manual instructions for that tool are printed. The same
list writes docs/install.md (`comrad setup --write-docs docs/install.md`).
"""
import os
import platform
import shutil
import subprocess
import sys

prog = os.path.basename(sys.argv[0])

# name: what to call it; exes: programs it must provide; module: names to look for;
# conda: package(s); build: bash run in <prefix>/src, `bins` relative to it; pip: package;
# url: where it comes from; for: which pipelines need it
BOWTIE2 = "2.5.4"      # release versions for the download/build route
SAMTOOLS = "1.21"
MINIFORGE = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-{arch}.sh"

TOOLS = [
    dict(name="ANGSD", exes=["angsd", "realSFS", "thetaStat"], module=["angsd"], conda="bioconda::angsd",
         build="git clone --recurse-submodules https://github.com/ANGSD/angsd.git && cd angsd && make",
         bins=["angsd/angsd", "angsd/misc/realSFS", "angsd/misc/thetaStat"],
         libs="zlib, bzip2, xz (lzma) and libcurl development files",
         url="https://github.com/ANGSD/angsd", needed="popgen, demography"),
    dict(name="NGSadmix", exes=["NGSadmix"], module=["ngsadmix", "angsd"], conda="",
         build="[ -d angsd ] || git clone --recurse-submodules https://github.com/ANGSD/angsd.git; cd angsd && make",
         bins=["angsd/misc/NGSadmix"], libs="zlib, bzip2, xz (lzma) and libcurl development files",
         url="https://github.com/ANGSD/angsd (built in misc/)", needed="popgen"),
    dict(name="bowtie2", exes=["bowtie2", "bowtie2-build"], module=["bowtie2"], conda="bioconda::bowtie2",
         build=f"V={BOWTIE2}; D=bowtie2-$V-linux-$(uname -m); U=https://github.com/BenLangmead/bowtie2/releases/download/v$V/$D.zip; "
               "(curl -fsSL -o $D.zip $U || wget -q -O $D.zip $U) && rm -rf $D && (unzip -oq $D.zip || python3 -m zipfile -e $D.zip .) "
               "&& chmod +x $D/bowtie2* && ln -sfn $D bowtie2",
         bins=["bowtie2/bowtie2", "bowtie2/bowtie2-build"], libs="perl (bowtie2's wrapper script), and unzip or python3",
         url="https://github.com/BenLangmead/bowtie2 (prebuilt release binaries)", needed="popgen"),
    dict(name="samtools", exes=["samtools"], module=["samtools"], conda="bioconda::samtools",
         build=f"V={SAMTOOLS}; U=https://github.com/samtools/samtools/releases/download/$V/samtools-$V.tar.bz2; "
               "(curl -fsSL -o samtools-$V.tar.bz2 $U || wget -q -O samtools-$V.tar.bz2 $U) && tar -xjf samtools-$V.tar.bz2 "
               "&& cd samtools-$V && ./configure --without-curses --prefix=$PWD/inst && make -j4 && make install && cd .. "
               "&& ln -sfn samtools-$V samtools",
         bins=["samtools/inst/bin/samtools"], libs="zlib, bzip2 and xz (lzma) development files",
         url="https://www.htslib.org/download/", needed="popgen"),
    dict(name="cutadapt", exes=["cutadapt"], module=["cutadapt"], conda="bioconda::cutadapt", pip="cutadapt",
         build="", url="https://cutadapt.readthedocs.io", needed="popgen"),
    dict(name="PCAngsd", exes=["pcangsd"], module=["pcangsd"], conda="bioconda::pcangsd",
         pip="git+https://github.com/Rosemeis/pcangsd.git", build="",
         url="https://github.com/Rosemeis/pcangsd", needed="popgen"),
    dict(name="evalAdmix", exes=["evalAdmix"], module=["evaladmix"], conda="",
         build="git clone https://github.com/GenisGE/evalAdmix.git && cd evalAdmix && make",
         bins=["evalAdmix/evalAdmix"], url="https://github.com/GenisGE/evalAdmix", needed="popgen"),
    dict(name="NGSremix", exes=["NGSremix"], module=["ngsremix"], conda="",
         build="git clone https://github.com/KHanghoj/NGSremix.git && cd NGSremix/src && make -f CPP_Makefile",
         bins=["NGSremix/src/NGSremix"], url="https://github.com/KHanghoj/NGSremix", needed="popgen"),
    dict(name="ngsParalog", exes=["ngsParalog"], module=["ngsparalog"], conda="",
         build="git clone https://github.com/tplinderoth/ngsParalog.git && cd ngsParalog && make",
         bins=["ngsParalog/ngsParalog"], url="https://github.com/tplinderoth/ngsParalog", needed="popgen"),
    dict(name="BaseSpace CLI", exes=["bs"], module=["basespace", "bs"], conda="", build="", optional=True,
         url="https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview",
         needed="reads, only to download from BaseSpace (BS_PROJECTS)"),
]
PY_DEMOGRAPHY = ["moments-popgen", "demes"]
PY_IMPORTS = {"moments-popgen": "moments", "demes": "demes"}


def say(msg=""):
    print(msg, flush=True)


def sh(cmd, cwd=None, log=None):
    """run bash; return (rc, output)"""
    r = subprocess.run(["bash", "-c", cmd], cwd=cwd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {cmd}\n{out}\n")
    return r.returncode, out


def have(exe, path):
    return shutil.which(exe, path=path)


def modules(names):
    """module names the module system offers for any of these (Lmod or Environment Modules)"""
    rc, out = sh("type module >/dev/null 2>&1 || source /etc/profile >/dev/null 2>&1; "
                 "type module >/dev/null 2>&1 || exit 9; " +
                 "; ".join(f"module -t avail {n} 2>&1" for n in names))
    if rc == 9:
        return None
    found = []
    for l in out.splitlines():
        l = l.strip().removesuffix("(default)").strip()
        if l and not l.endswith(":") and any(l.lower().split("/")[0] == n.lower() for n in names):
            found.append(l)
    return sorted(set(found))


def conda_cmd():
    for c in ("mamba", "micromamba", "conda"):
        if shutil.which(c):
            return c
    return None


def miniforge(pre, log):
    """install Miniforge into <pre>/miniforge (or reuse it) and return its mamba; None if it can't"""
    mf = os.path.join(pre, "miniforge")
    mamba = os.path.join(mf, "bin", "mamba")
    if os.access(mamba, os.X_OK):
        return mamba
    if platform.system() != "Linux":
        say("   Miniforge: only set up here on Linux")
        return None
    url = MINIFORGE.format(arch=platform.machine())
    inst = os.path.join(pre, "src", "Miniforge3.sh")
    say(f"   no conda or mamba: installing Miniforge into {mf} ...")
    rc, out = sh(f"(curl -fsSL -o {inst} {url} || wget -q -O {inst} {url}) && bash {inst} -b -p {mf}", log=log)
    if rc == 0 and os.access(mamba, os.X_OK):
        return mamba
    say(f"   Miniforge did not install (see {log}); going on without conda")
    return None


def link(src, bindir):
    dst = os.path.join(bindir, os.path.basename(src))
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.abspath(src), dst)


def manual(t):
    lines = [f"{t['name']} ({', '.join(t['exes'])}): needed by {t['needed']}", f"  source: {t['url']}"]
    if t.get("module"):
        lines.append(f"  module: module avail {t['module'][0]}   then add  module load <that name>  to ~/.bashrc")
    if t.get("conda"):
        lines.append(f"  conda:  conda install -c conda-forge -c bioconda {t['conda'].split('::')[-1]}")
    if t.get("pip"):
        lines.append(f"  pip:    python3 -m pip install --user {t['pip']}")
    if t.get("build"):
        lines.append(f"  build:  cd <folder for programs> && {t['build']}")
        lines.append(f"          then put {', '.join(t['bins'])} on PATH")
        lines.append(f"          (needs git, make, gcc/g++{', ' + t['libs'] if t.get('libs') else ' and zlib development files'})")
    return "\n".join(lines)


def markdown():
    """docs/install.md, from the same list setup uses"""
    out = ["# Install", "",
           "comRAD needs a Linux cluster with SLURM and Python 3.9 or newer. Installing it is two steps:",
           "comRAD itself with pip, then the external programs it runs with `comrad setup`.", "",
           "## 1. comRAD", "",
           "```bash",
           "module load python                 # whichever Python you will run with (if your cluster uses modules)",
           "python3 -m pip install --user git+https://github.com/RyanEckert/comRAD.git",
           "comrad -h",
           "```", "",
           "This installs `comrad`, `comrad-jobs`, every script the pipelines call, and the Python packages",
           "they use (numpy, scipy, matplotlib, moments, demes, PCAngsd and cutadapt). R is not needed.", "",
           "!!! tip \"Small home folder?\"",
           "    `pip --user` installs into `~/.local`. To put it on bigger storage instead, set",
           "    `PYTHONUSERBASE` before installing, and put its `bin` on your `PATH`:",
           "",
           "    ```bash",
           "    echo 'export PYTHONUSERBASE=/path/to/storage/.local' >> ~/.bashrc",
           "    echo 'export PATH=$PYTHONUSERBASE/bin:$PATH' >> ~/.bashrc",
           "    source ~/.bashrc",
           "    ```", "",
           "Load the same Python before installing that you will run comRAD with: pip writes that",
           "interpreter into each script's first line. If pip itself is missing, `python3 -m ensurepip",
           "--upgrade --user` installs it. To update comRAD later, run the same `pip install` with",
           "`--upgrade`; to install a tagged version, add `@v0.1.0` after `.git`.", "",
           "## 2. The programs comRAD runs", "",
           "```bash",
           "comrad setup --prefix /path/to/storage/comradDeps",
           "echo 'source /path/to/storage/comradDeps/comrad.env' >> ~/.bashrc",
           "source ~/.bashrc",
           "```", "",
           "Below, `DIR` is that `--prefix` folder; put it on storage with space, not a small home folder.",
           "For each program that is not already on your `PATH`, `comrad setup` tries, in order:", "",
           "1. an **environment module** (`module avail`); `comrad.env` loads the ones it finds;",
           "2. **conda/mamba** from bioconda, into `DIR/conda`. With no conda or mamba on the cluster, it",
           "   first installs [Miniforge](https://github.com/conda-forge/miniforge) into `DIR/miniforge`",
           "   (no admin rights needed); `--no-conda` skips conda altogether;",
           "3. a **download or build from source** in `DIR/src` (needs git, make and a C/C++ compiler);",
           "4. **pip**, for the Python tools.", "",
           "Everything it installs is linked into `DIR/bin`, which `comrad.env` puts on your `PATH`.",
           "Anything it could not install, it prints the steps for; they are also below.",
           "`comrad check popgen` (in a species folder) lists anything still missing.", "",
           "This page is written by `comrad setup --write-docs docs/install.md`, from the same list setup uses.", "",
           "## The programs", "",
           "| tool | programs | needed by | source |", "|---|---|---|---|"]
    for t in TOOLS:
        out.append(f"| {t['name']} | {', '.join('`' + e + '`' for e in t['exes'])} | {t['needed']} | {t['url']} |")
    out += [""]
    for t in TOOLS:
        out += [f"### {t['name']}", ""]
        if t.get("module"):
            out += [f"Module (most clusters): `module avail {t['module'][0]}`, then add `module load <name>` to `~/.bashrc`.", ""]
        if t.get("conda"):
            out += ["Conda:", "", "```bash", f"conda install -c conda-forge -c bioconda {t['conda'].split('::')[-1]}", "```", ""]
        if t.get("pip"):
            out += ["pip:", "", "```bash", f"python3 -m pip install --user {t['pip']}", "```", ""]
        if t.get("build"):
            out += [f"From source (needs git, make, a C/C++ compiler and {t.get('libs', 'zlib development files')}):", "", "```bash",
                    "cd DIR/src", t["build"], "```", "",
                    "Then put " + ", ".join(f"`{b}`" for b in t["bins"]) + " on PATH (for example, link them into `DIR/bin`).", ""]
        if not (t.get("conda") or t.get("build") or t.get("pip")):
            out += [f"See {t['url']}.", ""]
    out += ["## Conda without conda", "",
            "With no conda or mamba on PATH, setup installs Miniforge (conda-forge's installer; no admin",
            "rights, conda-forge channels only); by hand that is:", "", "```bash",
            "curl -fsSL -o Miniforge3.sh " + MINIFORGE.format(arch="$(uname -m)"),
            "bash Miniforge3.sh -b -p DIR/miniforge",
            "DIR/miniforge/bin/mamba create -y -p DIR/conda -c conda-forge -c bioconda angsd bowtie2 samtools",
            "```", "",
            "## Python packages", "",
            "numpy, scipy, matplotlib, moments-popgen, demes, PCAngsd and cutadapt are installed with comRAD.",
            "R is not needed.", ""]
    return "\n".join(out)


def run(a):
    if a.write_docs:
        open(a.write_docs, "w").write(markdown() + "\n")
        say(f"wrote {a.write_docs}")
        return 0
    if not a.prefix:
        say(f"{prog} setup: --prefix DIR is required: where programs are installed (not a small home folder)")
        return 2
    pre = os.path.abspath(os.path.expanduser(a.prefix))
    bindir, src = (os.path.join(pre, d) for d in ("bin", "src"))
    for d in (bindir, src):
        os.makedirs(d, exist_ok=True)
    log = os.path.join(src, "setup.log")
    path = bindir + os.pathsep + os.environ.get("PATH", "")
    cenv = os.path.join(pre, "conda")
    conda = None if a.no_conda else conda_cmd()
    want = [t for t in TOOLS if not a.only or t["name"].lower() in [o.lower() for o in a.only]]
    result, mods_wanted = {}, []
    say(f"installing into {pre} (build log: {log})")
    if not conda and not a.no_conda and os.access(os.path.join(pre, "miniforge", "bin", "mamba"), os.X_OK):
        conda = os.path.join(pre, "miniforge", "bin", "mamba")          # from an earlier setup
    say(f"conda: {conda or ('not used (--no-conda)' if a.no_conda else 'not found; Miniforge is installed if a tool needs it')}"
        + (f" (env {cenv})" if conda else ""))
    tried_mf = False
    for t in want:
        miss = [e for e in t["exes"] if not have(e, path)]
        if not miss:
            result[t["name"]] = "found: " + ", ".join(have(e, path) for e in t["exes"])
            continue
        say(f"\n== {t['name']}: not on PATH ({', '.join(miss)})")
        if not a.no_modules and t.get("module"):
            m = modules(t["module"])
            if m:
                say(f"   module available: {', '.join(m)}")
                result[t["name"]] = f"module: add  module load {m[-1]}  to ~/.bashrc"
                mods_wanted.append(m[-1])
                continue
            if m is None:
                say("   no module system found")
        if t.get("optional") and not a.optional:
            result[t["name"]] = "not installed (optional; --optional to try)"
            say(manual(t))
            continue
        done = False
        if not conda and not a.no_conda and t.get("conda") and not tried_mf:
            tried_mf = True
            conda = miniforge(pre, log)
        if conda and t.get("conda"):
            pkg = t["conda"].split("::")[-1]
            verb = "install" if os.path.isdir(cenv) else "create"
            say(f"   conda {verb} {pkg} ...")
            rc, out = sh(f"{conda} {verb} -y -p {cenv} -c conda-forge -c bioconda {pkg}", log=log)
            if rc == 0 and all(os.path.exists(os.path.join(cenv, "bin", e)) for e in t["exes"]):
                for e in t["exes"]:
                    link(os.path.join(cenv, "bin", e), bindir)
                done = True
                result[t["name"]] = f"conda ({cenv})"
            else:
                say(f"   conda did not provide it (see {log})")
        if not done and t.get("build"):
            say("   building from GitHub ...")
            rc, out = sh(t["build"], cwd=src, log=log)
            got = [os.path.join(src, b) for b in t["bins"]]
            if rc == 0 and all(os.access(g, os.X_OK) for g in got):
                for g in got:
                    link(g, bindir)
                done = True
                result[t["name"]] = f"built in {src}"
            else:
                say("   build failed; last lines:")
                say("     " + "\n     ".join(out.strip().splitlines()[-8:]))
        if not done and t.get("pip"):
            say(f"   pip install --user {t['pip']} ...")
            rc, out = sh(f"{sys.executable} -m pip install --user {t['pip']}", log=log)
            if rc == 0:
                ub = subprocess.run([sys.executable, "-c", "import site; print(site.USER_BASE)"],
                                    capture_output=True, text=True).stdout.strip()
                for e in t["exes"]:
                    p = os.path.join(ub, "bin", e)
                    if os.path.exists(p):
                        link(p, bindir)
                done = all(have(e, path) for e in t["exes"])
                if done:
                    result[t["name"]] = "pip --user"
        if not done:
            result[t["name"]] = "MISSING: install by hand (below)"
            say(manual(t))

    # python packages for demography
    # moments-popgen and demes install with comRAD; this only reports a python that lacks them
    miss = [p for p in PY_DEMOGRAPHY
            if subprocess.run([sys.executable, "-c", f"import {PY_IMPORTS[p]}"], capture_output=True).returncode]
    if miss:
        say(f"\n== python packages: pip install --user {' '.join(miss)} ...")
        sh(f"{sys.executable} -m pip install --user {' '.join(miss)}", log=log)
        miss = [p for p in miss
                if subprocess.run([sys.executable, "-c", f"import {PY_IMPORTS[p]}"], capture_output=True).returncode]
    result["moments, demes"] = "found" if not miss else f"MISSING: python3 -m pip install --user {' '.join(miss)}"

    env = os.path.join(pre, "comrad.env")
    with open(env, "w") as fh:
        fh.write(f"# written by comrad setup\nexport PATH={bindir}:$PATH\n")
        for m in mods_wanted:
            fh.write(f"module load {m} 2>/dev/null\n")
    say("\n== summary")
    w = max(len(k) for k in result)
    for k, v in result.items():
        say(f"  {k:<{w}}  {v}")
    say(f"\nadd to ~/.bashrc:  source {env}")
    say("then check:        comrad check popgen   (in a species folder)")
    return 1 if any(v.startswith("MISSING") for v in result.values()) else 0
