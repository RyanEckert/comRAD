[ -n "$BS_PROJECTS" ] || skip "BS_PROJECTS is empty; reads already in $RAW"
mkdir -p "$(abs "$RAW")" && cd "$(abs "$RAW")"
job 'ls *.fastq.gz >/dev/null 2>&1' "fetchReads.py -p $BS_PROJECTS" >FetchJobs
comrad-jobs FetchJobs -n fetchReads --res download || stop "download failed"
note "$(ls *.fastq.gz | wc -l) .fastq.gz"
