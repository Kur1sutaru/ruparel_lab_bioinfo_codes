# 02 — Loading Data, QC, and Filtering

## Goal
Go from a Cell Ranger output folder to a clean, QC'd Seurat object.

## Concepts

Not every barcode in the filtered matrix is a healthy single cell. Three failure modes to catch:
- **Low-quality / dying cells** — few genes detected, high mitochondrial-gene percentage (leaking cytoplasmic RNA through a compromised membrane).
- **Doublets** — two cells captured in one droplet; look like a hybrid transcriptional profile, inflate gene/UMI counts.
- **Ambient RNA contamination** — free-floating RNA from lysed cells gets incorporated into other droplets, subtly contaminating otherwise-real cells with off-lineage transcripts.

## Step 1 — Load the data

```r
library(Seurat)
library(tidyverse)

counts <- Read10X_h5("path/to/filtered_feature_bc_matrix.h5")
# or, if you have the folder form instead of .h5:
# counts <- Read10X("path/to/filtered_feature_bc_matrix/")

obj <- CreateSeuratObject(counts = counts, project = "my_project", min.cells = 3, min.features = 200)
```

`min.cells = 3` drops genes detected in fewer than 3 cells; `min.features = 200` drops cells with fewer than 200 genes detected — a permissive first pass, not your final filter.

## Step 2 — Compute QC metrics

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   # "^mt-" for mouse
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

Look at the violin plots before picking any threshold — don't copy-paste generic cutoffs blindly. Typical starting points for a healthy tissue dissociation:
- `nFeature_RNA`: keep roughly 200–6,000 (upper bound flags likely doublets; adjust per tissue — some cell types are legitimately high-complexity)
- `percent.mt`: keep below ~10–20% (tissue-dependent; some cell types, e.g. metabolically active ones, run higher at baseline)

## Step 3 — Ambient RNA estimation with decontX

```r
library(celda)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(obj)
decon <- decontX(sce, z = obj$seurat_clusters)  # run after a first-pass clustering, or omit z for an unsupervised pass
obj$contamination <- decon$contamination
```

Flag clusters with high mean contamination for closer review rather than automatically discarding them — genuinely activated/secretory cell states can show elevated ambient-RNA-like signal without being contaminated. Cross-check against off-lineage marker module scores before deciding.

## Step 4 — Doublet checks (optional but recommended)

If you have expected multiplet rates from your loading concentration (10x provides these), consider `scDblFinder` or `DoubletFinder` to flag likely doublets explicitly, especially for high-loading experiments.

## Step 5 — Apply filters

```r
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
```

## Checkpoint

Re-run the violin plots after filtering — distributions should look substantially cleaner, with the long high-`percent.mt` and high-`nFeature_RNA` tails gone. Record how many cells you started with and how many survived filtering (a sentence for your methods section later, e.g. "12,450 of 14,200 barcodes passed QC").

## Common pitfalls

- Applying the same hard thresholds across every sample/species without looking at each one's own distributions first.
- Filtering on `percent.mt` before checking the mitochondrial gene prefix matches your species' annotation (`MT-` human, `mt-` mouse, but custom/non-model references may differ — check `grep("^MT", rownames(obj), ignore.case = TRUE)`).
- Treating decontX contamination score as a hard filter rather than a flag to investigate.

Next: [`03_normalization_clustering.md`](03_normalization_clustering.md)
