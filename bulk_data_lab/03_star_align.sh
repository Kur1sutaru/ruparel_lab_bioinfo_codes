#!/usr/bin/env bash
# STAR index build (once) + alignment (per sample), single-end 50 bp reads.
# Requires: STAR, samtools on PATH.
#
# Reference: mouse GRCm39 (mm39) primary assembly + GENCODE vM current annotation.
# Swap for mm10 + matching GTF if your lab standardizes on that build instead —
# just keep genome FASTA and GTF from the same release.

set -euo pipefail

GENOME_DIR="genome/STAR_index"
GENOME_FASTA="genome/GRCm39.primary_assembly.genome.fa"
GTF="genome/gencode.vM_current.annotation.gtf"
READ_LENGTH=50           # this dataset: 50 bp SE reads
SJDB_OVERHANG=$((READ_LENGTH - 1))
THREADS=16

SAMPLE_SHEET="sample_sheet.csv"
TRIMDIR="fastq_trimmed"
ALIGNDIR="bam"

mkdir -p "$GENOME_DIR" "$ALIGNDIR"

# --- 1. Download reference (uncomment once; pick the current GENCODE mouse release) ---
# wget -P genome https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/GRCm39.primary_assembly.genome.fa.gz
# wget -P genome https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/gencode.vM_current.annotation.gtf.gz
# gunzip genome/*.gz

# --- 2. Build STAR index (run once) ---
if [[ ! -f "$GENOME_DIR/SAindex" ]]; then
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang "$SJDB_OVERHANG" \
        --runThreadN "$THREADS"
fi

# --- 3. Align each sample ---
tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r srr sample_name sex condition group; do
    echo ">>> Aligning $sample_name"
    mkdir -p "$ALIGNDIR/$sample_name"
    STAR \
        --runMode alignReads \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$TRIMDIR/${sample_name}_trimmed.fastq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$ALIGNDIR/$sample_name/" \
        --runThreadN "$THREADS" \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04

    samtools index "$ALIGNDIR/$sample_name/Aligned.sortedByCoord.out.bam"
done

echo "Done. Check bam/*/Log.final.out for alignment rate per sample before proceeding."
echo "Also worth running RSeQC infer_experiment.py on 1-2 BAMs to confirm library strandedness"
echo "(SMART-seq2 is typically unstranded, but verify before featureCounts -s)."
