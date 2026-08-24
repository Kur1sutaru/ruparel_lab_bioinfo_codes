# 00 — Bulk RNA-seq: Overview

## Goal
Understand the shape of the bulk RNA-seq pipeline before running any of it.

## The pipeline, end to end

```
FASTQs (from sequencer)
   │
   ▼
Read QC (FastQC/MultiQC) + trimming (fastp)
   │
   ▼
Alignment/quantification (Salmon, or STAR + featureCounts)
   │
   ▼
Gene-level count matrix
   │
   ▼
DESeq2: normalization + differential expression   ──► DEG table, volcano/MA plots
   │
   ▼
Enrichment analysis (GSEA / clusterProfiler)       ──► pathway-level interpretation
```

Each arrow is a tutorial in this folder (01 → 04).

## Key concepts to have in mind going in

- **Bulk RNA-seq measures the average expression across all cells in a sample.** You lose cell-type resolution but gain statistical power for detecting condition-level differences — this is the complementary strength to single-cell.
- **You need biological replicates, not technical ones, for statistics.** DESeq2's variance estimation assumes replicates capture real biological variability. Two conditions with n=1 each cannot be tested statistically — this needs to be planned before sequencing, not fixed afterward.
- **Pseudo-alignment (Salmon) vs. alignment (STAR)** — Salmon is much faster and sufficient for standard gene/transcript-level quantification; STAR gives you an actual BAM file, useful if you also need splice junction detection, variant calling, or visual inspection in a genome browser. For a straightforward DEG analysis, Salmon is the more common modern default.
- **GSEA vs. over-representation analysis (ORA)** are different questions. ORA (what `clusterProfiler::enrichGO`/`enrichKEGG` do) asks "are my significant DEGs enriched for this pathway?" using a hard significance cutoff. GSEA asks "is this pathway's gene set skewed toward one end of my whole ranked gene list?" without needing a cutoff — more sensitive to coordinated small shifts across many genes. Tutorial 04 covers both.

## Experimental design basics (read before you sequence, not after)

- Minimum 3 biological replicates per condition for DESeq2 to estimate dispersion reliably; more is better, especially for detecting subtle effects.
- If samples were processed/sequenced across multiple batches, record batch as a covariate — you'll include it in the DESeq2 design formula.
- Decide your comparisons of interest before sequencing so the design (paired vs. unpaired, how many conditions) actually supports the questions you want to ask.

## What you'll produce by the end of this track

- A gene-level count matrix
- A DESeq2 results table (DEGs with log2FC, p-values, adjusted p-values) per comparison
- Volcano and MA plots
- GSEA and ORA enrichment results (GO/KEGG) for your DEG set

Start with [`01_fastq_qc.md`](01_fastq_qc.md).
