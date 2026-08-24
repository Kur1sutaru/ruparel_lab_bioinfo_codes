# 03 — Differential Expression with DESeq2

## Goal
Go from a count matrix to a DEG table with volcano and MA plots.

## Step 1 — Build the DESeqDataSet

```r
library(DESeq2)

# coldata: a data.frame, one row per sample, rownames matching counts colnames exactly,
# with a "condition" column (and "batch" if relevant)
stopifnot(all(colnames(counts) == rownames(coldata)))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition          # add batch: ~ batch + condition
)
```

Design formula order matters conceptually: list nuisance covariates (e.g. `batch`) before the variable of interest — DESeq2 tests the last term by default, controlling for the others.

## Step 2 — Prefilter low-count genes

```r
keep <- rowSums(counts(dds) >= 10) >= min(table(coldata$condition))
dds <- dds[keep, ]
```

Not strictly required (DESeq2 has independent filtering built in), but speeds up computation and avoids genes with essentially no signal cluttering the results.

## Step 3 — Run DESeq2

```r
dds <- DESeq(dds)
resultsNames(dds)   # confirm the coefficient name you want to extract
```

## Step 4 — Extract results for a comparison

```r
res <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
res <- lfcShrink(dds, contrast = c("condition", "treated", "control"), res = res, type = "ashr")
# lfcShrink moderates noisy log2FC estimates for low-count genes — use the shrunk
# values for ranking/plotting, not the raw results() log2FC

summary(res)
```

## Step 5 — Tidy and export

```r
library(tidyverse)

res_df <- as.data.frame(res) |>
  rownames_to_column("gene") |>
  arrange(padj)

dir.create("DEG_results", showWarnings = FALSE)
write.csv(res_df, "DEG_results/treated_vs_control_DEGs.csv", row.names = FALSE)

sig_degs <- res_df |> filter(padj < 0.05, abs(log2FoldChange) > 1)
```

## Step 6 — Diagnostic and results plots

```r
plotMA(res, ylim = c(-5, 5))

library(EnhancedVolcano)   # BiocManager::install("EnhancedVolcano")
EnhancedVolcano(res_df,
  lab = res_df$gene, x = "log2FoldChange", y = "padj",
  pCutoff = 0.05, FCcutoff = 1)

# PCA of samples — check that replicates cluster together and conditions separate
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")
```

Always check the PCA plot before trusting the DEG list — if replicates don't cluster together, or a suspected batch effect dominates PC1, address that (revisit the design formula, check for a sample swap/outlier) before interpreting results.

## Checkpoint

`DEG_results/treated_vs_control_DEGs.csv` with `gene`, `log2FoldChange`, `pvalue`, `padj` columns; a PCA plot showing sensible replicate clustering; a volcano plot.

## Common pitfalls

- Using raw (non-shrunk) `log2FoldChange` for ranking or plotting — genes with low counts can show huge, noisy fold changes without shrinkage.
- Forgetting to include a known batch/covariate in the design formula, inflating apparent condition effects.
- Filtering on `pvalue` instead of `padj` (multiple testing correction) when calling significance.
- Not checking that `coldata` rows are in the same order as `counts` columns — DESeq2 assumes this silently.

Next: [`04_enrichment_gsea_clusterprofiler.md`](04_enrichment_gsea_clusterprofiler.md)
