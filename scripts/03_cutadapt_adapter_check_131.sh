#!/bin/bash
set -euo pipefail

cd "$HOME/PD_Gut_Microbiome"

RUN_LIST="data/metadata/PRJEB27564_Romano_matched_131_runs.txt"
INPUT_DIR="data/processed/PRJEB27564/cutadapt_trimmed_131"
LOG_DIR="results/qc/PRJEB27564/linked_adapter_diagnostic_131"

if [[ -e "$LOG_DIR" ]]; then
    echo "Diagnostic directory already exists. Stopping."
    exit 1
fi

if [[ $(wc -l < "$RUN_LIST" | tr -d ' ') -ne 131 ]]; then
    echo "Expected 131 sample IDs. Stopping."
    exit 1
fi

mkdir -p "$LOG_DIR"

count=0

while IFS= read -r run || [[ -n "$run" ]]; do
    [[ -z "$run" ]] && continue

    r1="$INPUT_DIR/${run}_R1_trimmed.fastq.gz"
    r2="$INPUT_DIR/${run}_R2_trimmed.fastq.gz"

    [[ -s "$r1" && -s "$r2" ]] || {
        echo "Missing input for $run. Stopping."
        exit 1
    }

    echo "Checking $run ($((count + 1))/131)"

    cutadapt \
        -j 4 -e 0.10 -O 45 \
        -a 'GGATTAGATACCCBDGTAGTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' \
        -a 'GGATTAGATACCCBDGTAGTCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' \
        -a 'GGATTAGATACCCBDGTAGTCAGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' \
        -a 'GGATTAGATACCCBDGTAGTCCACTCAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' \
        -A 'CTGCWGCCNCCCGTAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT' \
        -A 'CTGCWGCCNCCCGTAGGACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT' \
        -A 'CTGCWGCCNCCCGTAGGCTCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT' \
        -A 'CTGCWGCCNCCCGTAGGACACTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT' \
        --action=none \
        -o /dev/null -p /dev/null \
        "$r1" "$r2" \
        > "$LOG_DIR/${run}.log" 2>&1

    count=$((count + 1))

done < "$RUN_LIST"

echo "================================"
echo "Adapter diagnostic completed."
echo "Samples checked: $count"
echo "Logs: $LOG_DIR"
