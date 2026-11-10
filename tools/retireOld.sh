#!/bin/bash
# After installing comRAD: move the old copies of its scripts out of a folder
# on PATH (your $koko/bin), so they cannot shadow the installed ones, then show
# which copy each name now runs.
#   bash tools/retireOld.sh $koko/bin
old=$1
[ -d "$old" ] || { echo "usage: bash tools/retireOld.sh <folder>" >&2; exit 1; }
cd "$(dirname "$0")/.." || exit 1
mkdir -p "$old/replacedByRadTools"
while read -r s; do
case $s in ''|\#*) continue;; esac
[ -e "$old/$s" ] && mv -n "$old/$s" "$old/replacedByRadTools/" && echo "moved   $old/$s -> $old/replacedByRadTools/"
done <tools/SCRIPTS
hash -r
bad=0
while read -r s; do
case $s in ''|\#*) continue;; esac
w=$(command -v $s) || w="nothing (not on PATH)"; case $w in "$HOME"/.local/bin/*) ;; *) echo "CHECK   $s runs $w"; bad=1;; esac
done <tools/SCRIPTS
[ $bad -eq 0 ] && echo "every script now runs from ~/.local/bin"
