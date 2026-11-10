cd "$(abs "$RAW")" || stop "no $RAW"
n=$(ls *.fastq.gz 2>/dev/null | grep -vc '\.demux\.'); [ $n -gt 0 ] || stop "no .fastq.gz pools in $PWD"
trim2bRAD.py -f .fastq.gz -n DemuxJobs $CFG_demux_args || stop "trim2bRAD.py failed"
comrad-jobs DemuxJobs -n demux --res demux -c "grep -q 'reads in' {1}" -x="> '?([^ ']+\.demux\.log)" || stop "demultiplexing failed"
