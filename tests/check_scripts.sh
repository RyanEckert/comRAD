#!/bin/bash
# every script in tools/SCRIPTS is in scripts/, executable, and starts with #!
cd "$(dirname "$0")/.." || exit 1
bad=0
while read -r s; do
case $s in ''|\#*) continue;; esac
f=scripts/$s
if [ ! -e "$f" ]; then echo "missing: $f"; bad=1; continue; fi
[ -x "$f" ] || { echo "not executable: $f"; bad=1; }
[ "$(head -c2 "$f")" = '#!' ] || { echo "no #! line: $f"; bad=1; }
done <tools/SCRIPTS
# and pip installs every one of them
for s in $(grep -v '^#' tools/SCRIPTS); do grep -q "\"scripts/$s\"" pyproject.toml || { echo "not in pyproject.toml script-files: $s"; bad=1; }; done
[ $bad -eq 0 ] && echo "all $(grep -cv '^#\|^$' tools/SCRIPTS) scripts present, executable, with #!"
exit $bad
