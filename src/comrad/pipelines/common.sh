# Helpers every step can use. comrad puts this after the config exports.
#   stop  "why"   something is wrong that a rerun will not fix: the run stops
#   pause "what"  a decision is needed: the run stops, `comrad run` continues
#   skip  "why"   this step does not apply: marked done, the run goes on
#   note  "text"  copied into the run summary
stop(){ echo "STOP: $*" >&2; exit 1; }
pause(){ echo "PAUSE: $*" >&2; exit 99; }
skip(){ echo "SKIP: $*" >&2; exit 98; }
note(){ echo "NOTE: $*" >&2; }
# inskip [STEP]: true when STEP (default: this step) is in the config's SKIP list
inskip(){ [[ " $SKIP " == *" ${1:-$RT_STEP} "* ]]; }

# need FILE...: stop unless every file exists and is not empty
need(){ local f; for f in "$@"; do [ -s "$f" ] || stop "missing or empty: $f"; done; }

# equal A B what: stop unless two counts match
equal(){ [ "$1" -eq "$2" ] || stop "$3: $1 vs $2"; }

# yes VALUE: true for yes/true/1/on
yes(){ case "${1,,}" in yes|true|1|on) return 0;; *) return 1;; esac; }

# converged LOG...: every realSFS EM in these logs stopped before its -maxIter
# cap, or ended within 1e-3 of converging. Parallel tasks share a log, so each
# spectrum is read from its own line at the cap.
converged(){
  local e out
  [ $# -gt 0 ] || { echo "converged: no logs to check" >&2; return 0; }
  out=$(for e in "$@"; do [ -s "$e" ] || continue
    awk -v f="$e" 'match($0, /maxiter:[0-9]+/) {x = substr($0, RSTART+8, RLENGTH-8)+0; if (x > cap) cap = x}
      /^likelihood:/ {n++}
      /^lik\[[0-9]+\]/ {split($1, a, /[][]/); if (a[2]+0 == cap) {c++; d=$2; sub(/diff=/, "", d); if (d+0 > m) m = d+0}}
      END {s = (m > 1e-3) ? "CHECK" : "ok"; printf "%s\t%s\t%d spectra\t%d at the %d cap\tlargest diff at cap %g\n", s, f, n, c, cap, m}' "$e"
  done)
  echo "$out"
  if grep -q '^CHECK' <<<"$out"; then stop "a spectrum did not converge; raise -maxIter (see CHECK above)"; fi
}

# rtlogs NAME [GLOB]: the .e logs of the jobs that made NAME's current outputs (comrad-jobs
# writes .rtjobs/NAME/logs). Jobs stamped before logs were recorded fall back to GLOB
# (default NAME*.e*).
rtlogs(){
  local f=".rtjobs/${1//[^A-Za-z0-9_]/_}/logs" g=${2:-$1*.e*}
  if [ -s "$f" ]; then grep -v '^FALLBACK$' "$f"; if grep -qx FALLBACK "$f"; then ls $g 2>/dev/null; fi
  else ls $g 2>/dev/null; fi | sort -u
}

# job CHECK COMMAND: one line of an comrad-jobs job list (check<TAB>command)
job(){ printf '%s\t%s\n' "$1" "$2"; }

# abs PATH: PATH as absolute, relative ones taken from the working folder
abs(){ case "$1" in /*) echo "$1";; *) echo "$WD/$1";; esac; }

# tools NAME...: stop unless each is on PATH
tools(){ local t m=""; for t in "$@"; do command -v "$t" >/dev/null 2>&1 || m="$m $t"; done; [ -z "$m" ] || stop "not on PATH:$m (comrad setup --prefix DIR installs the external tools; see docs/install.md)"; }

# fststats FST.IDX: "unweighted<TAB>weighted" from realSFS fst stats; fails unless both are numbers
fststats(){
  local v; v=$(localSaf.py realSFS fst stats "$1" 2>/dev/null | tail -1 | tr -s ' \t' '\t' | sed 's/^\t//')
  [[ $v =~ ^-?[0-9.eE+-]+$'\t'-?[0-9.eE+-]+$ ]] || { echo "realSFS fst stats $1 gave '$v'" >&2; return 1; }
  echo "$v"
}
# sitesindex FILE: angsd sites index, with FILE dated two minutes back first. angsd refuses an index
# older than its sites file (and exits 0 without output); on network file systems a file just written
# can carry a later time than the index made right after it
sitesindex(){
  touch -d '2 minutes ago' "$1" && rm -f "$1.bin" "$1.idx" && angsd sites index "$1" 2>/dev/null &&
    [ -s "$1.bin" ] && [ -s "$1.idx" ] || stop "angsd sites index $1 failed"
}
