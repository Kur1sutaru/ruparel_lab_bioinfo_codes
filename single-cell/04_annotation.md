# 04 — Cluster Annotation

## Goal
Assign biologically meaningful cell type labels to each cluster, and be able to defend each call with markers.

## Concepts

Annotation is inference built on convergent evidence, not a single lookup. Combine: (1) canonical marker gene expression, (2) top DEGs per cluster (preview — full workflow in Tutorial 05), and (3) knowledge of what cell types are expected in your tissue.

## Step 1 — Build a canonical marker panel first, before looking at your data

Compile a short marker gene list per expected major cell type from literature/atlases for your tissue, e.g.:

```r
markers <- list(
  Epithelial  = c("EPCAM", "KRT5", "KRT14"),
  Immune      = c("PTPRC", "CD3D", "CD68"),
  Endothelial = c("PECAM1", "CDH5", "VWF"),
  Fibroblast  = c("COL1A1", "DCN", "LUM", "PDGFRA"),
  Muscle      = c("ACTA2", "MYH11", "DES"),
  Neural      = c("S100B", "SOX10", "MPZ")
)
```

Deciding markers up front (rather than post-hoc, after seeing which genes are differential) keeps annotation grounded in known biology instead of whatever happens to be a top DEG.

## Step 2 — Visualize markers across clusters

```r
DotPlot(obj, features = unlist(markers), group.by = "seurat_clusters") + RotatedAxis()
FeaturePlot(obj, features = c("EPCAM", "PTPRC", "PECAM1", "COL1A1"))
VlnPlot(obj, features = c("EPCAM", "PTPRC"), group.by = "seurat_clusters")
```

Look for dot plots where each expected cell type block "lights up" cleanly in a distinct set of clusters. Overlapping or diffuse signal across many clusters for one marker block usually means either the marker isn't specific enough for your tissue, or that cell type is more heterogeneous than one label captures (may need subclustering).

## Step 3 — Assign labels

```r
new_labels <- c(
  "0" = "Epithelial",
  "1" = "Immune",
  "2" = "Fibroblast",
  "3" = "Endothelial"
  # ... one entry per cluster number
)
obj$cell_type <- new_labels[as.character(obj$seurat_clusters)]
Idents(obj) <- obj$cell_type
DimPlot(obj, label = TRUE) + NoLegend()
```

## Step 4 — Sanity-check ambiguous or mixed clusters

If a cluster shows markers from two lineages at once, before assigning a label, check:
- Does it also show elevated ambient RNA contamination (Tutorial 02)? Likely a QC issue, not a real hybrid cell type.
- Does subclustering it (`FindSubCluster()` or re-running `FindClusters` on just that subset) resolve into two clean populations? Then it was under-clustered, not one ambiguous type.
- Is it a genuinely known transitional/doublet-prone population in your tissue (e.g. epithelial-mesenchymal transition zones)? Then the mixed signature may be real biology worth flagging explicitly rather than forcing a single label.

## Step 5 — Subclustering when one cluster hides multiple cell types

```r
sub <- subset(obj, idents = "Immune")
sub <- FindVariableFeatures(sub) |> ScaleData() |> RunPCA() |> FindNeighbors(dims = 1:15) |> FindClusters(resolution = 0.4)
DimPlot(sub, label = TRUE)
```

Re-run marker inspection (Steps 1–2) on the subcluster object with a finer, lineage-specific marker panel (e.g. T cell vs. B cell vs. myeloid markers within an "Immune" supercluster).

## Checkpoint

Every cluster has a label you can justify with at least 2–3 specific marker genes, and the DotPlot shows a clean, mostly non-overlapping block structure. Save:

```r
saveRDS(obj, "obj_annotated.rds")
```

## Common pitfalls

- Annotating from UMAP shape/proximity alone without checking markers.
- Using a marker gene that's actually broadly expressed (housekeeping-adjacent) as if it were cell-type-specific — check specificity, not just "is it high in this cluster."
- Forcing every cluster into a single confident label when the honest answer is "ambiguous, needs subclustering" — an unresolved or provisional label is more accurate than a wrong confident one.

Next: [`05_degs_and_markers.md`](05_degs_and_markers.md)
