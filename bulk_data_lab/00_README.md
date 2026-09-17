# Bulk RNA-seq tutorial — GSE220448 (tongue-innervating TG neurons)

Dataset: Ibrahim et al. 2023, Sci Rep (Ruparel lab). FACS-sorted, WGA-retro-labeled
tongue-innervating trigeminal ganglion neurons from Balb/c nude mice.
4 groups, n = 3/group, SMART-seq2 libraries, 50 bp single-end reads.

| Group | Sex    | Tissue condition   | n |
|-------|--------|--------------------|---|
| MN    | Male   | Normal (naive)     | 3 |
| MT    | Male   | Tumor-bearing      | 3 |
| FN    | Female | Normal (naive)     | 3 |
| FT    | Female | Tumor-bearing      | 3 |

Note: the original paper aligned to mm9 with TopHat2/HTSeq/DESeq (v1). This tutorial
uses a modern equivalent (STAR + GENCODE vM current mouse assembly + DESeq2) — same
biology, current tools. Swap in mm10 instead of GRCm39/mm39 if your lab standardizes
on that build for other datasets; just keep the genome+GTF pair consistent.

## Pipeline order

1. `01_download_sra.sh` — get raw fastq from SRA for the samples in GSE220448
2. `02_qc_trim.sh` — FastQC + fastp adapter/quality trimming
3. `03_star_align.sh` — STAR genome index + alignment (single-end, 50 bp)
4. `04_featurecounts.sh` — gene-level raw counts + gene lengths (featureCounts)
5. `05_deseq2_analysis.R` — RPKM, DESeq2, annotation, DEG cutoffs, Excel export

## Before you start: build the sample sheet

Go to the GEO record and pull the SRA Run Selector table:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220448
(link near the bottom: "SRA Run Selector" → download `SraRunTable.csv`)

Match each SRR run to its GSM/sample title (MN1-3, MT1-3, FN1-3, FT1-3 or however
they're labeled) and save as `sample_sheet.csv` in this folder:

```csv
srr,sample_name,sex,condition,group
SRRxxxxxxx,MN1,Male,Normal,MN
SRRxxxxxxx,MN2,Male,Normal,MN
SRRxxxxxxx,MN3,Male,Normal,MN
SRRxxxxxxx,MT1,Male,Tumor,MT
...
SRRxxxxxxx,FT3,Female,Tumor,FT
```

All the scripts below read this file — fill it in once and everything downstream
is driven off it.

## DEG cutoff logic (tissue-dependent)

This is encoded as a configurable function in step 5, not hardcoded, because your
cutoffs change by tissue/cell type. For **tongue tissue / tongue-innervating
neurons** (this dataset), the rule is:

- RPKM > 1 in **at least one** of the two groups being compared
- FC > 1.5 **or** FC < 0.6
- padj < 0.05

If that yields too few DEGs (you decide the threshold — the script defaults to
flagging when n < 10), fall back to raw **p-value < 0.05** instead of padj, with
everything else the same. The script reports which cutoff was actually used per
comparison so it's never silent about it.

For a different tissue/cell type later, you'd just call the same function with a
different `rpkm_cutoff` (the original tongue-neuron paper itself used RPKM > 5 —
worth deciding as a lab which cutoff you standardize on, since it's stricter than
your RPKM > 1 rule).

## Output

`05_deseq2_analysis.R` produces one Excel workbook
(`GSE220448_DE_results.xlsx`) with:
- one sheet per comparison, all genes, columns: Gene symbol, Gene description,
  RPKM per group, normalized counts per group, FC, Log2FC, pvalue, padj, DEG flag
- one sheet per comparison filtered to DEGs only
- a summary sheet (comparison, cutoff used, n up, n down)
