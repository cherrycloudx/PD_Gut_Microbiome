#!/bin/bash
set -euo pipefail

cd "$HOME/PD_Gut_Microbiome"

RUN_LIST="data/metadata/PRJEB27564_Romano_matched_131_runs.txt"
RAW_DIR="data/raw/PRJEB27564"
OUT_DIR="data/processed/PRJEB27564/cutadapt_trimmed_131"
LOG_DIR="results/qc/PRJEB27564/cutadapt_trim_logs_131"

if [[ -e "$OUT_DIR" || -e "$LOG_DIR" ]]; then
    echo "Output directory already exists. Stopping to prevent overwriting."
    exit 1
fi

if [[ ! -f "$RUN_LIST" ]]; then
    echo "Sample list not found."
    exit 1
fi

sample_count=$(wc -l < "$RUN_LIST" | tr -d ' ')

if [[ "$sample_count" -ne 131 ]]; then
    echo "Expected 131 samples, found $sample_count."
    exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"

count=0

while IFS= read -r run || [[ -n "$run" ]]; do

    [[ -z "$run" ]] && continue

    r1="$RAW_DIR/${run}_1.fastq"
    r2="$RAW_DIR/${run}_2.fastq"

    out1="$OUT_DIR/${run}_R1_trimmed.fastq.gz"
    out2="$OUT_DIR/${run}_R2_trimmed.fastq.gz"

    log="$LOG_DIR/${run}.log"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "Missing input FASTQ for $run. Stopping."
        exit 1
    fi

    echo "Processing $run ($((count + 1))/131)"

    cutadapt \
        -j 4 \
        -e 0.15 \
        -O 12 \
        -g '^CCTACGGGNGGCWGCAG' \
        -g '^GTCCTACGGGNGGCWGCAG' \
        -g '^AGAGCCTACGGGNGGCWGCAG' \
        -g '^TAGTGTCCTACGGGNGGCWGCAG' \
        -G '^GACTACHVGGGTATCTAATCC' \
        -G '^AGACTACHVGGGTATCTAATCC' \
        -G '^TCTGACTACHVGGGTATCTAATCC' \
        -G '^CTGAGTGGACTACHVGGGTATCTAATCC' \
        -o "$out1" \
        -p "$out2" \
        "$r1" "$r2" \
        > "$log" 2>&1

    count=$((count + 1))

    echo "Completed $run"

done < "$RUN_LIST"

echo "================================"
echo "Primer trimming completed."
echo "Samples processed: $count"
echo "Output directory: $OUT_DIR"
echo "Log directory: $LOG_DIR"
