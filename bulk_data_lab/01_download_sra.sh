#!/usr/bin/env bash
# Download raw fastq files for GSE220448 from SRA.
# Requires: sra-tools (prefetch, fasterq-dump) on PATH.
# Reads sample_sheet.csv (columns: srr,sample_name,sex,condition,group) in this dir.

GSE220448 accession cross-reference

Kept separate from sample_sheet.csv (which the pipeline scripts parse directly) so nothing here breaks the 5-column format the scripts expect. Use this table if you ever need to trace a sample back to its GEO/BioSample record.

srr	sample_name	GSM	BioSample	Library (SRX)
SRR22565240	MT1	GSM6803610	SAMN32093942	SRX18528542
SRR22565241	MT2	GSM6803609	SAMN32093943	SRX18528541
SRR22565242	MT3	GSM6803608	SAMN32093944	SRX18528540
SRR22565243	MN1	GSM6803607	SAMN32093945	SRX18528539
SRR22565244	MN2	GSM6803604	SAMN32093947	SRX18528537
SRR22565247	MN3	GSM6803606	SAMN32093946	SRX18528538
SRR22565245	FT1	GSM6803603	SAMN32093948	SRX18528536
SRR22565246	FT2	GSM6803602	SAMN32093949	SRX18528535
SRR22565248	FT3	GSM6803601	SAMN32093950	SRX18528534
SRR22565249	FN1	GSM6803600	SAMN32093951	SRX18528533
SRR22565250	FN2	GSM6803599	SAMN32093952	SRX18528532
SRR22565251	FN3	GSM6803598	SAMN32093953	SRX18528531

Assignment logic: GEO's condition field distinguishes tumor ("mouse tongue tumor induced by human oral squamous cell carcinoma cell lines HSC3") vs normal ("normal trigeminal ganglia (TG) tissues"); gender gives sex. SRR/GSM/BioSample fall into contiguous blocks of 3 per condition x sex combination, consistent with n=3/group in the published design.




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
