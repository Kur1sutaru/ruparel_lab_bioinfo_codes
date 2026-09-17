#!/usr/bin/env bash
# Download raw fastq files for GSE220448 from SRA.
# Requires: sra-tools (prefetch, fasterq-dump) on PATH.
# Reads sample_sheet.csv (columns: srr,sample_name,sex,condition,group) in this dir.

set -euo pipefail

SAMPLE_SHEET="sample_sheet.csv"
OUTDIR="fastq_raw"
THREADS=8

mkdir -p "$OUTDIR"

tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r srr sample_name sex condition group; do
    echo ">>> Downloading $srr ($sample_name / $group)"
    prefetch "$srr" -O "$OUTDIR"
    fasterq-dump "$OUTDIR/$srr/$srr.sra" \
        -O "$OUTDIR" \
        -e "$THREADS" \
        --split-3

    # Single-end data (this dataset is 50 bp SE) -> rename to sample_name for clarity
    if [[ -f "$OUTDIR/${srr}.fastq" ]]; then
        mv "$OUTDIR/${srr}.fastq" "$OUTDIR/${sample_name}.fastq"
        gzip -f "$OUTDIR/${sample_name}.fastq"
    fi

    rm -rf "$OUTDIR/$srr"
done

echo "Done. Raw fastq.gz files are in $OUTDIR/"
