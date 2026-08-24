# 04 — Enrichment Analysis: GSEA and clusterProfiler (ORA)

## Goal
Move from a gene-level DEG table to pathway-level interpretation, using both GSEA (rank-based) and over-representation analysis (cutoff-based).

## Concepts recap

- **ORA (`enrichGO`/`enrichKEGG`)**: takes a fixed gene list (e.g. significant DEGs) and asks whether it overlaps a pathway's gene set more than expected by chance. Simple, but throws away magnitude/direction information and depends heavily on your significance cutoff.
- **GSEA (`clusterProfiler::GSEA`, powered by `fgsea`)**: takes your *entire* ranked gene list (ranked by, e.g., log2FC or a signed significance score) and asks whether a pathway's genes are skewed toward one end. No hard cutoff needed; more sensitive to coordinated small effects across many genes.

Run both — they answer related but distinct questions and often surface complementary pathways.

## Step 0 — Map gene IDs if needed

`clusterProfiler` functions generally expect Entrez IDs for KEGG, and can take SYMBOL or ENSEMBL for GO with the right `keyType`. Map once:

```r
library(clusterProfiler)
library(org.Hs.eg.db)   # swap for your species' annotation package

res_df$entrez <- mapIds(org.Hs.eg.db, keys = res_df$gene, keytype = "SYMBOL", column = "ENTREZID")
```

## Step 1 — Over-representation analysis (GO)

```r
sig_genes <- res_df |> filter(padj < 0.05, abs(log2FoldChange) > 1) |> pull(gene)

ego <- enrichGO(
  gene          = sig_genes,
  universe      = res_df$gene,        # background = all tested genes, not the whole genome
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",               # Biological Process; also "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

dotplot(ego, showCategory = 20)
```

Setting `universe` to the full tested gene set (not the default whole-genome background) matters — it avoids inflated significance from genes that were never expressed/tested in your experiment to begin with.

## Step 2 — Over-representation analysis (KEGG)

```r
ekegg <- enrichKEGG(
  gene         = na.omit(res_df$entrez[res_df$gene %in% sig_genes]),
  universe     = na.omit(res_df$entrez),
  organism     = "hsa",              # KEGG species code
  pvalueCutoff = 0.05
)

dotplot(ekegg, showCategory = 20)
```

## Step 3 — GSEA

```r
ranked <- res_df |>
  filter(!is.na(padj)) |>
  arrange(desc(log2FoldChange)) |>
  { \(df) setNames(df$log2FoldChange, df$gene) }()

gsea_go <- gseGO(
  geneList     = ranked,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05
)

gseaplot2(gsea_go, geneSetID = 1, title = gsea_go$Description[1])
dotplot(gsea_go, showCategory = 20)
```

Ranking by raw `log2FoldChange` is a common simple choice; ranking by a signed statistic that also accounts for significance (e.g. `sign(log2FoldChange) * -log10(pvalue)`) is an alternative worth trying if the two give meaningfully different pathway lists.

## Step 4 — Export results

```r
dir.create("enrichment_results", showWarnings = FALSE)
write.csv(as.data.frame(ego), "enrichment_results/GO_ORA_BP.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "enrichment_results/KEGG_ORA.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_go), "enrichment_results/GO_GSEA_BP.csv", row.names = FALSE)
```

## Checkpoint

Three enrichment result tables in `enrichment_results/`, plus dot plots for each, and at least one `gseaplot2` enrichment plot for a top pathway that you can explain in plain language (what the pathway is, and whether it's up or down in your comparison).

## Common pitfalls

- Using the whole-genome default background instead of your experiment's tested-gene universe for ORA — inflates significance.
- Running GSEA on a gene list with duplicate names or `NA`s in the ranking vector — clean these first (`ranked <- ranked[!is.na(ranked) & !duplicated(names(ranked))]`).
- Over-interpreting a single significant pathway without checking the leading-edge genes actually driving it (`gsea_go@result$core_enrichment`) — a pathway can be "significant" while being driven by only 2-3 genes out of a large gene set.
- Mixing up whether a positive `log2FoldChange` in your design means "up in treated" or "up in control" — always double check the `contrast` direction from Tutorial 03 before interpreting GSEA direction.

This completes the bulk RNA-seq track. See [`../README.md`](../README.md) for the full curriculum map.
