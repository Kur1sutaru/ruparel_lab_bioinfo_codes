# 05 — Differential Expression and a Top-50-DEGs-per-Cluster Table

## Goal
Produce a clean, exportable table of the top 50 marker genes per cluster/cell type — the standard supplementary-table format for a paper.

## Concepts

"Marker genes" here means: genes differentially expressed in one cluster **vs. all other cells** (one-vs-rest), computed with a Wilcoxon rank-sum test by default in Seurat. This answers "what distinguishes this cluster" — a different question from "what changes between two conditions within a cluster" (same tools, different grouping variable — see the note at the end).

## Step 1 — Make sure you're on the right assay/join state

```r
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)   # required if the object went through split-layer integration (Tutorial 03)
```

## Step 2 — Run FindAllMarkers

```r
Idents(obj) <- obj$cell_type   # or seurat_clusters, if annotation isn't finalized yet

all_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,          # positive markers only — genes UP in this cluster
  min.pct = 0.25,           # expressed in at least 25% of cells in the cluster
  logfc.threshold = 0.25
)
```

## Step 3 — Extract top 50 per cluster

```r
library(dplyr)

top50 <- all_markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 50) |>
  ungroup()
```

If you want top-50-by-significance instead of by fold change:

```r
top50_by_pval <- all_markers |>
  group_by(cluster) |>
  arrange(p_val_adj, desc(avg_log2FC)) |>
  slice_head(n = 50) |>
  ungroup()
```

## Step 4 — Export

```r
dir.create("DEG_comparisons", showWarnings = FALSE)
write.csv(top50, "DEG_comparisons/top50_DEGs_per_cluster.csv", row.names = FALSE)
```

Keeping DEG outputs in a dedicated `DEG_comparisons/` folder (rather than scattered in the project root) makes it much easier to find things again once you have dozens of comparisons across cell types/conditions.

## Step 5 — Visual cross-check

```r
top5 <- top50 |> group_by(cluster) |> slice_head(n = 5) |> pull(gene) |> unique()
DotPlot(obj, features = top5, group.by = "cell_type") + RotatedAxis()
```

The top DEGs should visually correspond to the cell type identity you assigned in Tutorial 04. If a cluster's top DEGs look nothing like its assigned biology, revisit the annotation.

## Going further: DEGs between conditions (not just between clusters)

Once you have a clustered/annotated object with a condition variable (e.g. treated vs. control, tumor vs. normal), the same function works within a specific cell type by subsetting first:

```r
sub <- subset(obj, idents = "Fibroblast")
Idents(sub) <- sub$condition
condition_degs <- FindMarkers(sub, ident.1 = "treated", ident.2 = "control", min.pct = 0.1)
```

This is the more common analysis for hypothesis-driven questions ("what changes in fibroblasts with treatment?") versus the exploratory one-vs-rest table above.

## Checkpoint

`DEG_comparisons/top50_DEGs_per_cluster.csv` with columns including `cluster`, `gene`, `avg_log2FC`, `p_val_adj`, `pct.1`, `pct.2` — one block of ≤50 rows per cluster.

## Common pitfalls

- Running `FindAllMarkers` before `JoinLayers()` on an integrated object — will error or silently give wrong results.
- Reporting `avg_log2FC` without also checking `p_val_adj` — a large fold change in a handful of cells is not the same as a robust marker.
- Treating one-vs-rest markers as if they were condition DEGs, or vice versa — always be explicit in figure/table legends about which comparison was run.

This completes the single-cell track. See [`../bulk-rnaseq/00_overview.md`](../bulk-rnaseq/00_overview.md) for bulk RNA-seq.
