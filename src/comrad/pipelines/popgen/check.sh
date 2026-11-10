tools comrad-jobs cutadapt bowtie2 samtools angsd realSFS thetaStat pcangsd NGSadmix evalAdmix NGSremix ngsParalog \
  readCounts.py bowtie2Launcher.py plotQC.py filterBams.py paralogFilter.py pruneSnps.py pcangsdFisToText.py ngsAdmixLauncher.py pickQ.py pcangsdOutliers.py subsetBeagle.py cloneGroups.py admixResiduals.py plotStep.py kmeansBic.py localSaf.py mind
hasSym && tools splitHostSym.py
idx(){ ls "$1".1.bt2 "$1".1.bt2l >/dev/null 2>&1 || ls "$1".1.bt2* >/dev/null 2>&1 || stop "$1 has no bowtie2 index"; [ -s "$1.fai" ] || stop "$1 has no .fai (samtools faidx)"; }
need "$GENOME_FASTA"; idx "$GENOME_FASTA"
if hasSym; then need "$SYMGENOME" "$SYMGENUS"; idx "$SYMGENOME"; fi
[ -z "$EXTRA_READS" ] || ls "$EXTRA_READS"/*.filtered.fastq.gz >/dev/null 2>&1 || stop "no .filtered.fastq.gz in EXTRA_READS $EXTRA_READS"
need "$META"
case $CLONE_METRIC in auto|geno|ibs|both) ;; *) stop "CLONE_METRIC = $CLONE_METRIC: give auto, geno, ibs or both" ;; esac
num(){ [[ $1 =~ ^[0-9]*\.?[0-9]+$ ]]; }
case $CLONE_CUT in auto|gap|replicates) [ -n "$TECHREPS" ] || stop "CLONE_CUT = $CLONE_CUT needs TECHREPS; set TECHREPS, or CLONE_CUT to a number" ;;
*) num "$CLONE_CUT" || stop "CLONE_CUT = $CLONE_CUT: give auto, gap, replicates or a number" ;; esac
num "$TREE_BUFFER" || stop "TREE_BUFFER = $TREE_BUFFER: give a number, e.g. 1.10"
case ${KIN_DROP:-yes} in yes|no) ;; *) stop "KIN_DROP = $KIN_DROP: give yes or no" ;; esac
[ "${KIN_DROP:-yes}" = no ] && echo "WARNING: KIN_DROP = no: related individuals stay in: SFS, theta, Tajima's D and FST will include close kin, which skews allele frequencies toward their families"
[ -z "$HET_DEPTH" ] || [[ $HET_DEPTH =~ ^[0-9]*\.?[0-9]+$ ]] || stop "HET_DEPTH = $HET_DEPTH: leave it empty or give a mean depth, e.g. 5"
case $THETA_N in ""|min) ;; *[!0-9]*) stop "THETA_N = $THETA_N: leave it empty, or give a whole number or min" ;; *) [ $THETA_N -ge 2 ] || stop "THETA_N must be at least 2" ;; esac
for v in QDEPTH:$QDEPTH MIN_GROUP:$MIN_GROUP; do case ${v#*:} in ""|*[!0-9]*) stop "${v%%:*} = ${v#*:}: give a whole number" ;; esac; done
python3 - "$META" "$META_SAMPLE" "$WD" $META_GROUPS <<'PY' || stop "metadata check failed"
import csv, glob, os, re, sys
f, col, wd, *groups = sys.argv[1:]
rows = list(csv.DictReader(open(f, newline="", encoding="utf-8-sig")))
miss = [c for c in [col] + groups if c not in (rows[0].keys() if rows else [])]
if miss:
    sys.exit(f"{f}: no column(s) {', '.join(miss)}; it has {', '.join(rows[0].keys()) if rows else 'no rows'}")
names = [r[col] for r in rows]
dup = sorted({n for n in names if names.count(n) > 1})
meta = {r[col]: r for r in rows}
reads = sorted({re.sub(r"\.(demux|filtered)\.fastq\.gz$", "", os.path.basename(p)) for p in
                glob.glob(f"{wd}/demuxReads/*.demux.fastq.gz") + glob.glob(f"{wd}/filteredReads/*.filtered.fastq.gz")})
row = meta.get
nomatch = [s for s in reads if not row(s)]
blank = [s for s in reads if row(s) and any(not row(s)[g].strip() for g in groups)]
print(f"{f}: {len(rows)} rows; {len(reads)} samples with reads, all with a row")
for s in nomatch[:10]: print(f"  no row for {s}")
for s in blank[:10]: print(f"  blank {'/'.join(groups)} for {s}")
for s in dup[:10]: print(f"  {s} appears more than once")
if nomatch or blank or dup:
    sys.exit(f"{len(nomatch)} samples without a row, {len(blank)} with a blank group column, {len(dup)} duplicated names")
PY
# plot settings: factors must be META columns or lineage; listed levels should exist
python3 - "$META" "$META_SAMPLE" <<'PY' || stop "plot settings: see above"
import csv, os, re, sys
f, col = sys.argv[1:3]
e = os.environ
rows = list(csv.DictReader(open(f, newline="", encoding="utf-8-sig")))
cols = set(rows[0].keys()) if rows else set()
bad, warn = [], []
for k in ("PLOT_COLOR", "PLOT_SHAPE", "PLOT_ORDER"):
    for x in e.get(k, "").split():
        if x != "lineage" and x not in cols:
            bad.append(f"{k} names {x}, which is not a column of {f}")
kr = e.get("PLOT_KCOLOR", "").strip()
if kr:
    from matplotlib.colors import is_color_like
    if not is_color_like("#" + kr if re.fullmatch(r"[0-9A-Fa-f]{6}", kr) else kr):
        bad.append(f"PLOT_KCOLOR = {kr} is not a colour (hex without the #, or a name)")
for k in ("PLOT_LEVELS", "PLOT_COLORS", "PLOT_SHAPES"):
    for part in e.get(k, "").split(";"):
        if not part.strip():
            continue
        if ":" not in part:
            bad.append(f"{k}: '{part.strip()}' needs the form  factor: level level  (or level=value)")
            continue
        fac, rest = (x.strip() for x in part.split(":", 1))
        if fac != "lineage" and fac not in cols:
            bad.append(f"{k} names {fac}, which is not a column of {f}")
            continue
        for it in rest.split():
            lv, _, val = it.partition("=")
            if k != "PLOT_LEVELS" and not val:
                bad.append(f"{k}: '{it}' needs level=value")
            if fac != "lineage" and lv not in {r[fac] for r in rows}:
                warn.append(f"{k}: {fac} has no level {lv} in {f}")
for w in warn:
    print("WARNING:", w)
if bad:
    sys.exit("\n".join(bad))
PY
python3 -c 'import numpy, scipy, matplotlib' 2>/dev/null || stop "numpy, scipy or matplotlib missing for $(command -v python3); reinstall comRAD with this python (python3 -m pip install --user ...)"
# samples: demultiplexed reads, or already quality-filtered ones (a run that starts later)
READS=$(ls "$WD"/demuxReads/*.demux.fastq.gz 2>/dev/null; ls "$WD"/filteredReads/*.filtered.fastq.gz 2>/dev/null)
READS=$(printf '%s\n' $READS | sed 's#.*/##; s#\.\(demux\|filtered\)\.fastq\.gz$##' | sort -u)
[ -n "$READS" ] || [ -e "$WD/run/popgen/readcounts.done" ] || stop "no demuxReads/*.demux.fastq.gz or filteredReads/*.filtered.fastq.gz in $WD"
m=$(for t in $TECHREPS; do grep -qx "$t" <<<"$READS" || printf '%s ' $t; done)
[ -z "$m" ] || echo "WARNING: TECHREPS with no reads (left out of the clone threshold): $m"
[ -z "$READS" ] || echo "reads: $(wc -l <<<"$READS") samples"
