tools trim2bRAD.py sampleRename.py plotStep.py comrad-jobs
[ -s "$(abs "$SAMPLE_SHEET")" ] || [ -s "$(abs "$SAMPLE_SHEET").csv" ] || stop "SAMPLE_SHEET $SAMPLE_SHEET not found"
for p in $SPECIES; do [[ $p == ?*=?* ]] || stop "SPECIES entries must be species=prefix, got '$p'"; done
if [ -n "$BS_PROJECTS" ]; then
tools fetchReads.py bs
bs whoami >/dev/null 2>&1 || stop "bs is not authenticated; run: bs auth"
elif ! ls "$(abs "$RAW")"/*.fastq.gz >/dev/null 2>&1; then
stop "BS_PROJECTS is empty and there are no .fastq.gz in $(abs "$RAW")"
fi
