_bc=$WD/angsd/${SP}BamsClusters
[ -s "$_bc" ] || stop "no $_bc: run popgen first"
export POPS=$(awk '$2 != "Admixed" {print $2}' "$_bc" | sort -u | tr '\n' ' ')
pairs(){ set -- $POPS; while [ $# -gt 1 ]; do a=$1; shift; for b in "$@"; do echo "$a $b"; done; done; }
# fields in one line of a pair's 2D spectrum: (2 n_a + 1) (2 n_b + 1)
nfields(){ echo $(( (2 * $(wc -l <$WD/angsd/${SP}${1}Bams) + 1) * (2 * $(wc -l <$WD/angsd/${SP}${2}Bams) + 1) )); }
