cd "$(abs "$RAW")"
s=$(abs "$SAMPLE_SHEET")
sampleRename.py -i "$s" -f demux.fastq.gz $RENAME_ARGS --dry-run || stop "sampleRename.py --dry-run failed; nothing renamed"
sampleRename.py -i "$s" -f demux.fastq.gz $RENAME_ARGS || stop "sampleRename.py failed"
