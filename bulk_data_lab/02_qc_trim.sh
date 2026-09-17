#!/usr/bin/env bash
# QC (FastQC) + trimming (fastp) for single-end 50 bp reads.
# Requires: fastqc, fastp on PATH.

set -euo pipefail

SAMPLE_SHEET="sample_sheet.csv"
RAWDIR="fastq_raw"
TRIMDIR="fastq_trimmed"
QCDIR="qc/fastqc_raw"
QCDIR_TRIM="qc/fastqc_trimmed"
THREADS=8

mkdir -p "$TRIMDIR" "$QCDIR" "$QCDIR_TRIM"

# Raw QC (optional but recommended — check for adapter content / per-base quality
# before deciding trimming parameters, especially since SMART-seq2 can carry
# template-switching oligo adapter at the 5' end)
fastqc -t "$THREADS" -o "$QCDIR" "$RAWDIR"/*.fastq.gz

tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r srr sample_name sex condition group; do
    echo ">>> Trimming $sample_name"
    fastp \
        -i "$RAWDIR/${sample_name}.fastq.gz" \
        -o "$TRIMDIR/${sample_name}_trimmed.fastq.gz" \
        --thread "$THREADS" \
        --qualified_quality_phred 20 \
        --length_required 25 \
        --json "qc/${sample_name}_fastp.json" \
        --html "qc/${sample_name}_fastp.html"
done

fastqc -t "$THREADS" -o "$QCDIR_TRIM" "$TRIMDIR"/*.fastq.gz

echo "Done. Check qc/*_fastp.html and qc/fastqc_trimmed/ before proceeding —"
echo "if per-sample quality is uneven, adjust --qualified_quality_phred / --length_required accordingly."
