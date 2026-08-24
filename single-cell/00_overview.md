# 00 — Single-cell RNA-seq: Overview

## Goal
Understand the shape of the whole single-cell pipeline before running any of it, so each later step makes sense in context.

## The pipeline, end to end

```
FASTQs (from sequencer)
   │
   ▼
Cell Ranger (10x Cloud)  ──► raw feature-barcode matrix, per-cell metrics
   │
   ▼
Seurat: load + QC          ──► remove empty droplets, doublets, dying cells
   │
   ▼
Normalization + clustering ──► PCA → neighbors → Louvain/Leiden clusters → UMAP
   │
   ▼
Annotation                 ──► assign biological cell type labels to clusters
   │
   ▼
Differential expression    ──► marker genes / DEGs per cluster (top 50 table)
```

Each arrow is a tutorial in this folder (01 → 05).

## Key concepts to have in mind going in

- **A "cell" in the raw matrix might not be a real cell.** Empty droplets, ambient RNA, and doublets (two cells captured together) all produce barcodes that look like cells until you QC them out. Most of tutorial 02 is about this.
- **Clustering resolution is a choice, not a fact.** There's no single "correct" number of clusters — it depends on the biological question. We use `clustree` to see how clusters split as resolution increases and pick a resolution where splits look biologically meaningful rather than arbitrary.
- **Annotation is inference, not measurement.** You're assigning cell type labels based on marker gene expression matching known biology. Always sanity-check against canonical markers, not just visual UMAP separation.
- **DEGs depend on what you're comparing.** "Marker genes per cluster" (one cluster vs. rest) is a different question from "DEGs between conditions within a cluster." Tutorial 05 covers the former; the same logic extends to the latter once you're comfortable.

## Experimental design basics (read before you sequence, not after)

- Include enough biological replicates to support statistics later — a single sample per condition cannot support formal statistical testing between conditions.
- Balance conditions across 10x lanes/chips where possible to avoid confounding condition with batch.
- Decide up front whether you need cell hashing / multiplexing (multiple samples per lane) — this changes the Cell Ranger run configuration in Tutorial 01.

## What you'll produce by the end of this track

- A clustered, annotated Seurat object (`.rds`)
- A UMAP figure colored by cell type
- A top-50-DEGs-per-cluster CSV table, ready to drop into a supplementary file

Start with [`01_cellranger_cloud.md`](01_cellranger_cloud.md).
