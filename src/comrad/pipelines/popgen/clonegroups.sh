cd "$WD/angsd"
export CLONEDIR=.; [ -e ${SP}CloneBams ] || CLONEDIR=clones     # filterbams moves them
rm -f $CLONEDIR/${SP}CloneStop.txt $CLONEDIR/${SP}CloneDendIBS.pdf
num(){ [[ $1 =~ ^[0-9]*\.?[0-9]+$ ]]; }
rc=0
cloneGroups.py -i $CLONEDIR/${SP}Clones.ibsMat -g $CLONEDIR/${SP}Clones.geno.gz -b $CLONEDIR/${SP}CloneBams -r "$TECHREPS" \
  --metric $CLONE_METRIC --cut $CLONE_CUT --buffer $TREE_BUFFER ${CLONE_MAXDROP:+--maxdrop $CLONE_MAXDROP} \
  --quality ${SP}Quality${QDEPTH}X.txt --meta "$META" --meta-sample "$META_SAMPLE" --meta-groups $META_GROUPS ${PLOT_COLOR:+--color "$PLOT_COLOR"} ${PLOT_SHAPE:+--shape "$PLOT_SHAPE"} --levels "$PLOT_LEVELS" --colors "$PLOT_COLORS" --shapes "$PLOT_SHAPES" -o $SP --dir $CLONEDIR || rc=$?
case $rc in
0) ;;
3) stop "no technical replicate set has two or more bams, so nothing sets the clone cut; set  CLONE_CUT = <number>  in $CONFIG" ;;
4) stop "$(cat $CLONEDIR/${SP}CloneStop.txt)" ;;
*) stop "clone groups failed; see above" ;;
esac
cat $CLONEDIR/${SP}CloneGroups
# a numeric CLONE_CUT = the cut was chosen, so no pause
if yes "$CLONE_PAUSE" && ! num "$CLONE_CUT"; then
[ ! -e "$WD/angsd/CLONES_OK" ] || stop "angsd/CLONES_OK is no longer read; set  CLONE_PAUSE = no  in $CONFIG (or CLONE_CUT to a number) and delete angsd/CLONES_OK"
pause "look at $(realpath --relative-to="$WD" "$CLONEDIR/${SP}CloneDend.pdf") and ${SP}CloneBorderline.txt; to accept these groups set  CLONE_PAUSE = no  in $CONFIG, or set  CLONE_CUT = <number>"
fi
