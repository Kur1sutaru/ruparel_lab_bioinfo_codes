# RNA-seq Analysis Tutorials

A step-by-step curriculum for learning single-cell and bulk RNA-seq analysis end to end — from raw FASTQs to publication-style figures. Written for a new lab member with some R experience but no prior single-cell/bulk RNA-seq background.

The single-cell track follows a **Seurat (R)** workflow, matching the analysis conventions used in this lab's cross-species work. Where useful, each tutorial links out to the equivalent chapter in [Single-cell best practices](https://www.sc-best-practices.org/preamble.html) (Python/Scanpy-oriented, but the underlying statistics and QC logic are framework-agnostic and worth reading alongside these tutorials).

## How to use this repo

Work through the two tracks in order. Each `.md` file is self-contained but assumes you've done the ones before it. Every tutorial has:
- **Goal** — what you'll be able to do by the end
- **Concepts** — the "why," briefly, before the "how"
- **Steps** — runnable code, in order
- **Checkpoint** — what your output should look like, so you know if something went wrong
- **Common pitfalls**

Start with [`environment/setup.md`](environment/setup.md) to get R, Bioconductor, and command-line tools installed before touching any data.

## Track 1 — Single-cell RNA-seq (Seurat)

| # | Tutorial | Covers |
|---|----------|--------|
| 0 | [`single-cell/00_overview.md`](single-cell/00_overview.md) | Big-picture pipeline map, experimental design basics |
| 1 | [`single-cell/01_cellranger_cloud.md`](single-cell/01_cellranger_cloud.md) | 10x Genomics Cloud Analysis: account setup, uploading FASTQs, running `cellranger count`, downloading outputs |
| 2 | [`single-cell/02_qc_filtering.md`](single-cell/02_qc_filtering.md) | Loading Cell Ranger output into Seurat, QC metrics, ambient RNA/doublet checks, filtering |
| 3 | [`single-cell/03_normalization_clustering.md`](single-cell/03_normalization_clustering.md) | Normalization, HVGs, PCA, neighbor graph, clustering, resolution selection, UMAP |
| 4 | [`single-cell/04_annotation.md`](single-cell/04_annotation.md) | Marker-based manual annotation, dot plots/feature plots, sanity checks |
| 5 | [`single-cell/05_degs_and_markers.md`](single-cell/05_degs_and_markers.md) | `FindAllMarkers`, exporting a top-50-DEGs-per-cluster table |

## Track 2 — Bulk RNA-seq

| # | Tutorial | Covers |
|---|----------|--------|
| 0 | [`bulk-rnaseq/00_overview.md`](bulk-rnaseq/00_overview.md) | Pipeline map, experimental design (replicates, batch) |
| 1 | [`bulk-rnaseq/01_fastq_qc.md`](bulk-rnaseq/01_fastq_qc.md) | FastQC, MultiQC, adapter/quality trimming |
| 2 | [`bulk-rnaseq/02_alignment_quantification.md`](bulk-rnaseq/02_alignment_quantification.md) | Salmon (pseudo-alignment) and STAR + featureCounts (alignment-based); building a count matrix |
| 3 | [`bulk-rnaseq/03_deseq2_degs.md`](bulk-rnaseq/03_deseq2_degs.md) | DESeq2 design formulas, DEG tables, MA/volcano plots |
| 4 | [`bulk-rnaseq/04_enrichment_gsea_clusterprofiler.md`](bulk-rnaseq/04_enrichment_gsea_clusterprofiler.md) | GSEA (fgsea) and over-representation analysis (clusterProfiler: GO, KEGG) |

## Further reading

- [Single-cell best practices](https://www.sc-best-practices.org/preamble.html) — the deeper conceptual reference for the single-cell track
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Harvard Chan Bioinformatics Core training materials](https://hbctraining.github.io/) — good companion exercises for both tracks

## License

See [`LICENSE`](LICENSE).
