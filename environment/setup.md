# Environment Setup

Do this once, before starting either track.

## 1. R and RStudio

- Install R (≥ 4.3) from [CRAN](https://cran.r-project.org/)
- Install [RStudio Desktop](https://posit.co/download/rstudio-desktop/)

## 2. Core R packages

```r
install.packages(c("Seurat", "tidyverse", "patchwork", "clustree", "BiocManager"))

BiocManager::install(c(
  "DESeq2", "clusterProfiler", "org.Hs.eg.db", "fgsea",
  "AnnotationDbi", "tximport", "DropletUtils", "celda"   # celda = decontX
))
```

- `Seurat` — single-cell analysis
- `clustree` — visualizing clustering resolution choices
- `DESeq2` — bulk differential expression
- `clusterProfiler` + `org.Hs.eg.db` — over-representation and GSEA enrichment (swap `org.Hs.eg.db` for the correct species annotation package if not working with human data, e.g. `org.Mm.eg.db` for mouse)
- `fgsea` — fast GSEA implementation used under the hood by `clusterProfiler::GSEA()`
- `tximport` — importing Salmon quantifications into R
- `celda` — provides `decontX` for ambient RNA estimation

## 3. Command-line tools (bulk RNA-seq track)

Install via conda (recommended) in a dedicated environment:

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
  fastqc multiqc fastp salmon star subread samtools
conda activate rnaseq
```

- `fastqc` / `multiqc` — read QC
- `fastp` — adapter/quality trimming
- `salmon` — pseudo-alignment + quantification
- `star` — splice-aware alignment (alternative/complement to Salmon)
- `subread` (provides `featureCounts`) — counting reads per gene from a BAM
- `samtools` — BAM sorting/indexing

## 4. 10x Genomics Cloud Analysis account

No local install needed — Cell Ranger runs in the browser. See [`single-cell/01_cellranger_cloud.md`](../single-cell/01_cellranger_cloud.md) for account setup.

## 5. A note on reproducibility

Once you're comfortable with the basic workflow, consider locking package versions with `renv::init()` in each project directory so analyses are reproducible months later. Not required to get started.
