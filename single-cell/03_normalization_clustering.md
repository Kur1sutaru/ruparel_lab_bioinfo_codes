# 03 — Normalization, Dimensionality Reduction, and Clustering

## Goal
Take a filtered Seurat object to a UMAP with clusters you can defend the resolution choice for.

## Step 1 — Normalize and find variable features

```r
obj <- NormalizeData(obj)  # log-normalization, default: 10,000 scale factor
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(obj), 10)
VariableFeaturePlot(obj)
```

3,000 HVGs is a reasonable default; fewer (~2,000) if you're seeing noisy PCA structure downstream, more if you have very heterogeneous cell types and are worried about losing rare-population signal.

## Step 2 — Scale and run PCA

```r
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, npcs = 30, features = VariableFeatures(obj))
ElbowPlot(obj, ndims = 30)
```

Pick the number of PCs to carry forward where the elbow plot flattens out — commonly somewhere in the 15–30 range. Erring slightly high is safer than cutting off real variance.

## Step 3 — Neighbor graph and clustering

```r
obj <- FindNeighbors(obj, dims = 1:25)
```

### Picking a resolution with clustree

Rather than guessing one resolution, run a range and visualize how clusters split:

```r
library(clustree)

for (res in seq(0.1, 1.2, by = 0.1)) {
  obj <- FindClusters(obj, resolution = res)
}
clustree(obj, prefix = "RNA_snn_res.")
```

Read the tree bottom (low resolution, few clusters) to top (high resolution, many clusters). Look for the resolution just before clusters start splitting into groups with no clear marker-gene distinction (you'll check this in Tutorial 04) — that's usually your working resolution. It's fine to revisit this after annotation if two "clusters" turn out to be the same cell type split by a technical axis, or two true states got merged.

```r
obj <- FindClusters(obj, resolution = 0.5)  # replace with your chosen value
```

> **Housekeeping tip:** stale `RNA_snn_res.*` columns from earlier clustree runs can accumulate in `obj@meta.data` and cause `clustree()` to error on missing values if you ever subset/re-integrate the object later. If you hit an odd `clustree` NA error, check `colnames(obj@meta.data)` for leftover resolution columns from prior runs and remove them before re-running.

## Step 4 — UMAP

```r
obj <- RunUMAP(obj, dims = 1:25)
DimPlot(obj, reduction = "umap", label = TRUE) + NoLegend()
```

## A note on integration (multiple samples/species/batches)

If you're combining multiple samples, species, or batches, you generally do **not** want to cluster on raw merged PCA space — batch effects will dominate real biological structure. Two common approaches:

- **CCA integration** (`IntegrateLayers(..., method = CCAIntegration)`) — tends to work well when you expect the *same* cell types across batches/species and want them aligned even if compositions differ.
- **Harmony integration** (`IntegrateLayers(..., method = HarmonyIntegration)`) — often preserves finer activation-state substructure that CCA can collapse; a good choice when you care about within-cell-type states, not just coarse type identity.

A typical recipe (adapt `dims`/resolution to your data):

```r
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)  # or $species/$batch
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
obj <- IntegrateLayers(obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "cca")
obj <- FindNeighbors(obj, reduction = "cca", dims = 1:25)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, reduction = "cca", dims = 1:25)
```

Try both if you're unsure which suits your dataset, and compare whether known cell types separate cleanly without batch-driven sub-splitting.

## Checkpoint

A UMAP where clusters look visually separated and roughly match the number of expected major cell types for your tissue (you'll confirm identity in Tutorial 04). Save your progress:

```r
saveRDS(obj, "obj_clustered.rds")
```

## Common pitfalls

- Skipping the elbow plot and guessing a PC count.
- Picking clustering resolution purely by "how the UMAP looks" rather than checking (in Tutorial 04) whether each cluster has distinguishing marker genes.
- Clustering merged multi-sample/species data without integration and mistaking a batch effect for biology.
- Running `DotPlot()` directly on an integrated assay — it can hit a duplicate-rowname bug in some Seurat versions; if so, compute `AverageExpression()` manually and plot with `ggplot2` instead.

Next: [`04_annotation.md`](04_annotation.md)
