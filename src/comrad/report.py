"""comrad report: one HTML page per pipeline run, with the diary's warnings, pauses
and stops first, then the key tables. Written to run/<pipeline>/<pipeline>Report.html;
the popgen and demography tarballs include it. It only reads files, so it can be run
at any point, including on a stopped or paused run."""
import glob, html, os, re, time

CSS = """
:root { --bg:#fff; --fg:#1d1d1f; --mute:#666; --line:#ddd; --head:#f3f3f3; --warn:#fff4d6; --stop:#fde2e2; --pause:#e3eefc; }
@media (prefers-color-scheme: dark) { :root { --bg:#161616; --fg:#e8e8e8; --mute:#9a9a9a; --line:#333; --head:#222;
  --warn:#3a3016; --stop:#3d1c1c; --pause:#18283d; } }
body { background:var(--bg); color:var(--fg); font:14px/1.45 -apple-system, system-ui, sans-serif; max-width:1100px; margin:0 auto; padding:16px; }
h1 { font-size:22px; margin:8px 0 2px } h2 { font-size:17px; margin:28px 0 8px; border-bottom:1px solid var(--line); padding-bottom:4px }
h3 { font-size:14px; margin:16px 0 6px } .mute { color:var(--mute) } .src { color:var(--mute); font-size:12px; font-weight:normal }
table { border-collapse:collapse; font-size:13px; margin:4px 0 8px; display:block; overflow-x:auto; max-width:100% }
td, th { border:1px solid var(--line); padding:3px 8px; text-align:left; white-space:nowrap; font-variant-numeric:tabular-nums }
th { background:var(--head) } tr.WARNING td { background:var(--warn) } tr.STOP td { background:var(--stop) } tr.PAUSE td { background:var(--pause) }
pre { background:var(--head); padding:8px; overflow-x:auto; font-size:12px } details summary { cursor:pointer; color:var(--mute) }
"""


def _rows(path):
    out = []
    for l in open(path, errors="replace"):
        l = l.rstrip("\n")
        if not l.strip() or l.startswith("#"):
            continue
        out.append(l.split("\t") if "\t" in l else l.split())
    return out


def _num(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def table(rows, header=None, cls=None):
    """header: True, False, or None to guess (first row has no numbers and the second does)"""
    if not rows:
        return '<p class="mute">(empty)</p>'
    if header is None:
        header = len(rows) > 1 and not any(_num(x) for x in rows[0]) and any(_num(x) for x in rows[1])
    h = ""
    if header:
        h = "<tr>" + "".join(f"<th>{html.escape(x)}</th>" for x in rows[0]) + "</tr>"
        rows = rows[1:]
    body = "".join(f'<tr class="{cls(r) if cls else ""}">' + "".join(f"<td>{html.escape(x)}</td>" for x in r) + "</tr>" for r in rows)
    return f"<table>{h}{body}</table>"


class Page:
    def __init__(self, wd, pipeline="popgen"):
        self.wd, self.parts, self.pipeline = wd, [], pipeline

    def find(self, *cands):
        for c in cands:
            hits = sorted(glob.glob(os.path.join(self.wd, c)))
            if hits:
                return hits[0]
        return None

    def h2(self, t):
        self.parts.append(f"<h2>{html.escape(t)}</h2>")

    def add(self, title, *cands, header=None, pre=False, fold=40, note="", last_block=False):
        """a file as a table (or preformatted text); folded when it is long; left out if absent"""
        f = self.find(*cands)
        if not f:
            return
        rel = os.path.relpath(f, self.wd)
        if pre:
            txt = open(f, errors="replace").read()
            if last_block and "\n=====" in "\n" + txt:
                txt = "=====" + ("\n" + txt).split("\n=====")[-1]
            body = f"<pre>{html.escape(txt)}</pre>" if txt.strip() else '<p class="mute">(empty)</p>'
            n = body.count("\n")
        else:
            rows = _rows(f)
            body, n = table(rows, header), len(rows)
        self.parts.append(f'<h3>{html.escape(title)} <span class="src">{html.escape(rel)}</span></h3>' +
                          (f'<p class="mute">{html.escape(note)}</p>' if note else ""))
        self.parts.append(f"<details><summary>{n} lines</summary>{body}</details>" if n > fold else body)

    def images(self, kinds):
        """the step plots' PNGs (plotStep.py writes them to run/<pipeline>/plots), inline"""
        import base64
        d = os.path.join(self.wd, "run", self.pipeline, "plots")
        for k in kinds:
            fs = sorted(glob.glob(os.path.join(d, f"{k}-*.png")), key=lambda f: int(re.sub(r"\D", "", f.rsplit("-", 1)[1]) or 0))
            if fs:
                self.parts.append(f"<h3>{html.escape(k)}</h3>" + "".join(
                    f'<img alt="{html.escape(k)}" style="max-width:100%" src="data:image/png;base64,'
                    f'{base64.b64encode(open(f, "rb").read()).decode()}">' for f in fs))

    def files(self, title, *pats):
        fs = sorted({os.path.relpath(f, self.wd) for p in pats for f in glob.glob(os.path.join(self.wd, p))})
        if fs:
            self.parts.append(f"<h3>{html.escape(title)}</h3><ul>" + "".join(f"<li>{html.escape(f)}</li>" for f in fs) + "</ul>")


def diary(path):
    """the diary as (time, kind, step, text) rows; kind STOP/PAUSE/WARNING are what to read first"""
    out = []
    for l in open(path, errors="replace"):
        f = l.rstrip("\n").split("\t", 2)
        if len(f) < 3 or f[1] == "version":
            continue
        t, kind, msg = f
        if kind in ("done", "run"):
            step, text = msg.split(" ", 1)[0], msg[len(msg.split(" ", 1)[0]):].strip()
        else:
            step, _, text = msg.partition(": ") if ": " in msg else ("", "", msg)
        if kind == "note" and "WARNING" in text:
            kind = "WARNING"
        elif kind in ("STOP", "PAUSE"):
            pass
        out.append([t, kind, step, text])
    return out


def build(pl, sp):
    p, wd, name = Page(pl.wd, pl.name), pl.wd, pl.name
    S = sp or ""
    qd = pl.setting("QDEPTH", "5")
    try:
        from importlib.metadata import version
        ver = version("comrad")
    except Exception:
        ver = "?"
    p.parts.append(f"<h1>{html.escape(S + ' ' if S else '')}{html.escape(name)}</h1>"
                   f'<p class="mute">{html.escape(wd)} · written {time.strftime("%Y-%m-%d %H:%M")} · comrad {html.escape(ver)}</p>')

    d = diary(pl.summary) if os.path.exists(pl.summary) else []
    # the last run: from the last driver start
    last = max([i for i, r in enumerate(d) if r[1] == "start"], default=0)
    # stops and pauses still open (the step has not finished since); each warning once, latest kept
    flags, seen = [], set()
    for i, r in enumerate(d):
        if r[1] in ("STOP", "PAUSE") and not any(x[1] == "done" and x[2] == r[2] for x in d[i + 1:]):
            flags.append(r)
        elif r[1] == "WARNING":
            later = [x for x in d[i + 1:] if x[1] == "WARNING" and x[2] == r[2] and x[3] == r[3]]
            if not later and (r[2], r[3]) not in seen:
                seen.add((r[2], r[3]))
                flags.append(r)
    p.h2("Read first")
    if d:
        end = d[-1]
        p.parts.append(f'<p>Last event: <b>{html.escape(end[1])}</b> {html.escape(end[2])} {html.escape(end[3])} <span class="mute">({html.escape(end[0])})</span></p>')
    p.parts.append(table([["time", "kind", "step", "message"]] + flags, True, lambda r: r[1]) if flags else
                   '<p class="mute">no warnings, and no open pauses or stops, in the diary</p>')
    sk = [r for r in d[last:] if r[1] == "skip"]
    if sk:
        p.parts.append("<h3>Skipped steps</h3>" + table([["time", "kind", "step", "reason"]] + sk, True))

    if name == "popgen" and S:
        p.h2("Samples")
        bc = p.find(f"angsd/{S}BamsClusters")
        if bc:
            from collections import Counter
            c = Counter(r[1] for r in _rows(bc) if len(r) > 1)
            p.parts.append("<h3>Lineages</h3>" + table([["lineage", "n"]] + [[k, str(v)] for k, v in sorted(c.items())], True))
        p.add(f"Samples left out of the clone step (proportion of sites at >= {qd}x)", f"angsd/clones/{S}CloneLowCoverage.txt", f"angsd/{S}CloneLowCoverage.txt", header=False)
        p.add("Technical replicates: nearest sample by the clone distance", f"angsd/clones/{S}CloneReplicates.txt", f"angsd/{S}CloneReplicates.txt")
        p.add("Clone distance and cut", f"angsd/clones/{S}CloneMetric.txt", f"angsd/{S}CloneMetric.txt", header=False)
        p.add("Clone groups", f"angsd/clones/{S}CloneGroups", f"angsd/{S}CloneGroups", header=False)
        p.add("Borderline joins (with each side's lowest coverage)", f"angsd/clones/{S}CloneBorderline.txt", f"angsd/{S}CloneBorderline.txt")
        p.add("Kin dropped", f"relatedness/{S}KinDrop", header=False)
        p.add(f"Proportion of sites at >= {qd}x, all samples", f"angsd/{S}Quality{qd}X.txt", header=False, fold=0)
        p.add("Coverage at the SNPs, per sample", f"angsd/{S}SnpDepth.txt", fold=0)

        p.h2("Structure")
        p.add("K (last run of the step)", "run/popgen/K.log", pre=True, fold=60, last_block=True)
        p.add("NGSadmix fit per K", f"ngsAdmix/{S}_fit.txt")
        p.add("Neutral, kin-free run against the lineages", f"ngsAdmix/{S}NeutralCompare.txt")
        p.add("PCAngsd inbreeding (Fis)", f"pcangsd/{S}Fis.txt", fold=0)

        p.h2("Differentiation and diversity")
        p.add("FST between lineages (neutral sites)", f"sfs/{S}LineageFst.out")
        p.add("FST between lineages (all sites)", f"sfs/{S}LineageFstAll.out")
        p.add("FST between sampling groups within lineages", f"angsd/groups/{S}GroupFst.out")
        p.add("Theta", f"theta/{S}Theta.out",
              note="pi compares well across lineages of different size; thetaW and Tajima's D depend on n, so compare the Sub columns when THETA_N is set")
        p.add("Observed heterozygosity per lineage", f"het/{S}Ho.txt")
        p.add("Heterozygosity against coverage", f"het/{S}HetDepth.txt")
        p.add("Heterozygosity per sample", f"het/{S}.het", fold=0)

        p.h2("Plots")
        p.images(["reads", "symbionts", "pca", "admix", "kin", "outliers", "neutral", "admixneutral", "fst", "het", "theta"])
        p.files("PDFs (open from the run folder or the tarball)", f"angsd/{S}Depth.pdf", f"angsd/clones/{S}CloneDend*.pdf",
                f"angsd/{S}CloneDend*.pdf", f"angsd/{S}CloneFreeDend.pdf", f"ngsAdmix/{S}*.pdf", f"pcangsd/{S}*.pdf",
                f"filteredReads/{S}*.pdf", f"bams/{S}*.pdf", f"relatedness/{S}*.pdf", f"angsd/lineages/{S}*.pdf",
                f"sfs/{S}*.pdf", f"het/{S}*.pdf", f"theta/{S}*.pdf")

    if name == "demography" and S:
        p.h2("Models")
        p.add("Models per pair, ranked", f"demo/{S}MomentsBest.tsv")
        p.add("Godambe uncertainty", f"demo/{S}MomentsGodambe.tsv", fold=0)
        p.add("Convergence", f"demo/{S}Convergence.txt")
        p.add("Nested-model tests", f"demo/{S}Nesting.txt")
        p.add("Projection sizes", f"demo/{S}PopN.txt")
        p.add("Warnings", f"demo/{S}Warnings.txt", pre=True)
        p.files("Plots", f"demo/*/{S}*.pdf", f"demo/*/{S}*.png")

    if name == "reads":
        p.h2("Reads")
        p.add("Reads per barcode", "RawCounts", "*/RawCounts")
        p.images(["pools"])
        p.files("Species folders", "*/demuxReads")

    p.h2("Settings")
    p.add("Config", os.path.relpath(pl.conf_path, wd), pre=True, fold=0)
    p.add("Tool versions", f"run/{name}/versions.txt", fold=0)
    p.h2("Diary")
    p.parts.append(f"<details><summary>{len(d)} events</summary>" + table([["time", "kind", "step", "message"]] + d, True, lambda r: r[1]) + "</details>")

    title = f"{S} {name} report".strip()
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
            f"<title>{html.escape(title)}</title><style>{CSS}</style></head><body>" + "\n".join(p.parts) + "</body></html>\n")


def write(pl):
    sp = pl.setting("SP")
    out = os.path.join(pl.dir, f"{pl.name}Report.html")
    os.makedirs(pl.dir, exist_ok=True)
    open(out, "w").write(build(pl, sp))
    return out


def summary(pl):
    """plain-text digest of the diary: open pauses and stops, then each step's notes from its last run,
    each warning once"""
    head = f"== {pl.setting('SP') or os.path.basename(os.path.abspath(pl.wd))} {pl.name} ({pl.wd})"
    if not os.path.exists(pl.summary):
        return head + "\nnothing run yet"
    d = diary(pl.summary)
    raw = [l.rstrip("\n").split("\t") for l in open(pl.summary, errors="replace") if l.strip()]
    if raw and raw[-1][1:2] == ["FINISHED"]:
        d.append([raw[-1][0], "FINISHED", "", ""])
    out = [head]
    opened = [r for i, r in enumerate(d) if r[1] in ("STOP", "PAUSE")
              and not any(x[1] in ("done", "rewind") and (x[2] == r[2] or x[1] == "rewind") for x in d[i + 1:])]
    if d and d[-1][1] == "FINISHED":
        out.append("finished " + d[-1][0])
    for r in opened:
        out.append(f"{r[1]} at {r[2]}: {r[3]}")
        if "then: comrad" not in r[3]:
            out.append(f"   when that is done: comrad run {pl.name}")
    last = {}
    for i, r in enumerate(d):
        if r[1] == "run":
            last[r[2]] = i
    rew = [(i, r[3].split()[-1]) for i, r in enumerate(d) if r[1] == "rewind" and r[3].split()]
    pos = {s: i for i, s in enumerate(pl.steps)}
    seen = set()
    for step in pl.steps:
        if step not in last:
            continue
        if any(i > last[step] and pos.get(s, 1e9) <= pos[step] for i, s in rew):
            continue                                   # cleared by a later `comrad from`; not run since
        rows = [r for r in d[last[step]:] if r[2] == step and r[1] in ("note", "WARNING", "skip")
                and not re.search(r"\d+ jobs?, peak memory", r[3])]   # job resource lines: in the report
        rows = [r for r in rows if not (r[3] in seen or seen.add(r[3]))]
        if rows:
            out.append(f"-- {step}")
            out += [("   ! " if r[1] == "WARNING" else "     ") + r[3] for r in rows]
    if not opened and not (d and d[-1][1] == "FINISHED"):
        out.append(f"last event: {' '.join(d[-1][1:]) if d else 'none'}")
    return "\n".join(out)
