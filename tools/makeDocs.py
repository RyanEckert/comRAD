#!/usr/bin/env python3
################################################################################
#Script Name    : makeDocs.py
#Description    : Write the documentation pages that come from the code: every
#                 setting, every script's help, and the tool install guide
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, os, re, subprocess, sys

prog = os.path.basename(sys.argv[0])
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(
    description="Write the documentation pages that are generated from the code, so they cannot drift from it: "
                "docs/reference/settings.md (every setting and job section, from the pipeline templates), "
                "docs/reference/scripts.md (every script's --help), docs/reference/commands.md (comrad and comrad-jobs --help) and docs/install.md (from comrad setup's tool list).",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  python3 tools/{prog}                  # from the repo folder; rewrites the three pages
  python3 tools/{prog} --check          # exit 1 if any page is out of date (for CI)
  python3 tools/{prog} && mkdocs serve  # then look at the site on http://127.0.0.1:8000

Run it with a python that has comRAD's dependencies (numpy, scipy, matplotlib,
moments-popgen, demes): each script's --help is taken by running the script.
""")
ap.add_argument("--check", action="store_true", help="write nothing; exit 1 if a page would change")
a = ap.parse_args()

PIPES = os.path.join(ROOT, "src", "comrad", "pipelines")
sys.path.insert(0, os.path.join(ROOT, "src"))


def parse_template(path):
    """(intro comments, [(heading, [(key, value, comment)])], [(section, comment, [(key, value, comment)])])"""
    intro, groups, sections = [], [["", []]], []
    cur = None
    started = False
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"^\[(\w+)\]\s*(?:#\s*(.*))?$", line)
        if m:
            cur = (m.group(1), m.group(2) or "", [])
            sections.append(cur)
            continue
        m = re.match(r"^([A-Za-z_]\w*)\s*=\s*(.*?)\s*(?:#\s*(.*))?$", line)
        if m:
            started = True
            row = (m.group(1), m.group(2), m.group(3) or "")
            (cur[2] if cur else groups[-1][1]).append(row)
            continue
        if line.startswith("#") and cur is None:
            text = line.lstrip("#").strip()
            if not started:
                intro.append(text)
            elif text and not text.startswith("----") and len(text.split()) <= 4 and not groups[-1][1] == []:
                groups.append([text, []])
            elif text and len(text.split()) <= 4:
                groups[-1][0] = text
    return intro, [g for g in groups if g[1]], sections


def cell(x):
    """table text: | escaped, and <...> kept as text rather than read as HTML"""
    return x.replace("|", "\\|").replace("<", "&lt;").replace(">", "&gt;") if x else ""


def code(x):
    return f"`{x.replace('|', '&#124;')}`" if x else "*(empty)*"


def settings_page():
    out = ["# Settings", "",
           "Every setting, its default and what it does, for each pipeline. This page is written from the",
           "pipeline templates by `tools/makeDocs.py`; `comrad defaults <pipeline>` prints the same thing in",
           "the terminal, and `comrad init <pipeline>` writes it as a config to fill in.", "",
           "A config is `KEY = value` lines, then one `[section]` per job list. In a section, `cpus` and `mem`",
           "are per job (per array task), and `cpus` is also the tool's own thread count. `args` lists",
           "**changes** to the tool's default options: an option you give replaces the default's value or is",
           "added, and `-option none` removes it. `queue` (`long` = the cluster's `long_queue`) and `slurm`",
           "(extra sbatch flags as `key=value`) can go in any section. Settings are read from the template's",
           "defaults, then `cluster.config`, then the pipeline's own config; later lines win.", ""]
    for name, title in (("reads", "reads"), ("popgen", "popgen"), ("demography", "demography")):
        intro, groups, sections = parse_template(os.path.join(PIPES, name, "template.config"))
        out += [f"## {title}", ""]
        out += [f"Run from {'the project folder' if name == 'reads' else 'the species folder' + (', after popgen' if name == 'demography' else '')}."
                " **Required** settings must be filled in; `comrad check " + name + "` says which are missing.", ""]
        for heading, rows in groups:
            if heading:
                out += [f"### {heading.capitalize()}", ""]
            out += ["| setting | default | what it does |", "|---|---|---|"]
            for k, v, c in rows:
                req = c.lower().startswith("required")
                c = re.sub(r"^required:?\s*", "", c, flags=re.I)
                out.append(f"| `{k}`{' **required**' if req else ''} | {code(v)} | {cell(c)} |")
            out.append("")
        out += [f"### Job sections", "",
                "| section | runs | cpus | mem | default args |", "|---|---|---|---|---|"]
        for sec, comment, rows in sections:
            d = {k: v for k, v, _ in rows}
            extra = [f"`{k} = {v}`" for k, v, _ in rows if k not in ("cpus", "mem", "args") and v]
            out.append(f"| `[{sec}]` | {cell(comment)}{'<br>' + ', '.join(extra) if extra else ''} | {d.get('cpus', '')} | "
                       f"{d.get('mem', '') or '*(see note)*'} | {code(d.get('args', '')) if 'args' in d else ''} |")
        out.append("")
        if name == "popgen":
            out += ["`[map]` with `mem` empty asks for the bowtie2 index size × 1.5 + 2 GB per task.", ""]
    intro, groups, sections = parse_template(os.path.join(PIPES, "cluster.config"))
    out += ["## cluster.config", "",
            "One per project, written by `comrad init cluster`. comrad reads it from the folder a pipeline runs in",
            "or up to three folders above, so one in the project folder covers every species folder.", "",
            "| setting | default | what it does |", "|---|---|---|"]
    for sec, comment, rows in sections:
        for k, v, c in rows:
            out.append(f"| `{k}` | {code(v)} | {cell(c)} |")
    out.append("")
    return "\n".join(out)


def script_help(s):
    path = os.path.join(ROOT, "scripts", s)
    first = open(path).readline()
    if "python" in first:
        r = subprocess.run([sys.executable, path, "-h"], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"{prog}: {s} -h failed:\n{r.stderr}")
        return r.stdout
    # a shell script: its usage comment lines
    return "\n".join(l.lstrip("#").strip() for l in open(path).read().splitlines()[1:6] if l.startswith("#")) + "\n"


def scripts_page():
    names = [l.strip() for l in open(os.path.join(ROOT, "tools", "SCRIPTS")) if l.strip() and not l.startswith("#")]
    helps = {s: script_help(s) for s in names}
    out = ["# Scripts", "",
           "comRAD installs these scripts on your `PATH`. The pipelines call them, and each one also works on its",
           "own. Below is every script's `--help`, as `tools/makeDocs.py` captured it.", "",
           "| script | what it does |", "|---|---|"]
    for s in names:
        h = helps[s]
        desc = ""
        m = re.search(r"\n\n(.+?)(\n\n|$)", h, re.S)        # argparse: usage, blank, description
        if m and "python" in open(os.path.join(ROOT, "scripts", s)).readline():
            desc = " ".join(m.group(1).split())
        else:
            desc = h.splitlines()[0] if h else ""
        desc = desc.split(". ")[0].rstrip(".") + "." if desc else ""
        out.append(f"| [`{s}`](#{re.sub(r'[^a-z0-9]+', '', s.lower())}) | {cell(desc)} |")
    out.append("")
    for s in names:
        out += [f"## {s}", "", "```text", helps[s].rstrip("\n"), "```", ""]
    return "\n".join(out)


def cli_help(mod, args):
    r = subprocess.run([sys.executable, "-c", f"import sys; sys.argv = {args!r}; from {mod} import main; main()"],
                       capture_output=True, text=True, env={**os.environ, "PYTHONPATH": os.path.join(ROOT, "src"), "COLUMNS": "100"})
    return r.stdout or r.stderr


def commands_page():
    out = ["# Commands", "",
           "Two commands: `comrad`, which runs the pipelines, and `comrad-jobs`, which submits and watches job",
           "lists (the pipelines call it; it also works on its own). Below is each one's `--help`, as",
           "`tools/makeDocs.py` captured it.", "",
           "## comrad", "", "```text", cli_help("comrad.cli", ["comrad", "-h"]).rstrip(), "```", ""]
    for c in ("init", "check", "run", "status", "summary", "list", "from", "step", "report", "defaults", "setup"):
        out += [f"### comrad {c}", "", "```text", cli_help("comrad.cli", ["comrad", c, "-h"]).rstrip(), "```", ""]
    out += ["## comrad-jobs", "", "```text", cli_help("comrad.jobs", ["comrad-jobs", "-h"]).rstrip(), "```", ""]
    return "\n".join(out)


def install_page():
    from comrad import setup_tools
    return setup_tools.markdown() + "\n"


pages = {os.path.join(ROOT, "docs", "reference", "settings.md"): settings_page(),
         os.path.join(ROOT, "docs", "reference", "scripts.md"): scripts_page(),
         os.path.join(ROOT, "docs", "reference", "commands.md"): commands_page(),
         os.path.join(ROOT, "docs", "install.md"): install_page()}
stale = []
for path, text in pages.items():
    old = open(path).read() if os.path.exists(path) else None
    if old != text:
        stale.append(os.path.relpath(path, ROOT))
        if not a.check:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").write(text)
if a.check:
    if stale:
        sys.exit(f"{prog}: out of date: {', '.join(stale)}; run python3 tools/{prog}")
    print(f"{prog}: all generated pages are up to date")
else:
    print(f"{prog}: wrote {', '.join(stale) if stale else 'nothing (all up to date)'}")
