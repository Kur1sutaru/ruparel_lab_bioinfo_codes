#!/usr/bin/env bash
# Gene-level read counting with featureCounts (subread package).
# Produces raw counts + gene lengths in one table -> feeds directly into
# 05_deseq2_analysis.R for both DESeq2 (raw counts) and RPKM (counts + length).

set -euo pipefail

SAMPLE_SHEET="sample_sheet.csv"
ALIGNDIR="bam"
GTF="genome/gencode.vM_current.annotation.gtf"
OUTDIR="counts"
THREADS=16
STRANDEDNESS=0   # 0 = unstranded, 1 = stranded, 2 = reverse-stranded
                 # CONFIRM with RSeQC infer_experiment.py before running this —
                 # SMART-seq2 is usually 0, but check your actual libraries.

mkdir -p "$OUTDIR"

BAMS=()
NAMES=()
while IFS=',' read -r srr sample_name sex condition group; do
    BAMS+=("$ALIGNDIR/$sample_name/Aligned.sortedByCoord.out.bam")
    NAMES+=("$sample_name")
done < <(tail -n +2 "$SAMPLE_SHEET")

featureCounts \
    -a "$GTF" \
    -o "$OUTDIR/gene_counts.txt" \
    -T "$THREADS" \
    -s "$STRANDEDNESS" \
    -t exon \
    -g gene_id \
    --primary \
    "${BAMS[@]}"

echo "Done. counts/gene_counts.txt has Geneid, Chr, Start, End, Strand, Length, then one column per BAM."
echo "Rename the BAM-path column headers to your sample names before/inside R (handled in step 5)."
