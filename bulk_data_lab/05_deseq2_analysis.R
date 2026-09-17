## ============================================================================
## GSE220448 — tongue-innervating TG neurons, bulk RNA-seq DE analysis
## Input : counts/gene_counts.txt (featureCounts output) + sample_sheet.csv
## Output: GSE220448_DE_results.xlsx
## ============================================================================

## ---- 0. Packages -----------------------------------------------------------
required_pkgs <- c("DESeq2", "org.Mm.eg.db", "AnnotationDbi", "openxlsx", "dplyr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}

library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(openxlsx)
library(dplyr)

## ---- 1. Load sample sheet and counts ---------------------------------------
sample_sheet <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)
rownames(sample_sheet) <- sample_sheet$sample_name
sample_sheet$group <- factor(sample_sheet$group, levels = c("MN", "MT", "FN", "FT"))

fc <- read.delim("counts/gene_counts.txt", comment.char = "#", check.names = FALSE)

# featureCounts BAM-path columns -> sample names (match order in sample_sheet)
bam_cols <- colnames(fc)[7:ncol(fc)]
# assumes bam paths look like bam/<sample_name>/Aligned.sortedByCoord.out.bam
matched_names <- sapply(bam_cols, function(x) {
  sample_sheet$sample_name[sapply(sample_sheet$sample_name, function(s) grepl(paste0("/", s, "/"), x))]
})
colnames(fc)[7:ncol(fc)] <- matched_names

gene_length <- setNames(fc$Length, fc$Geneid)      # bp, for RPKM
counts_mat  <- as.matrix(fc[, sample_sheet$sample_name])
rownames(counts_mat) <- fc$Geneid
storage.mode(counts_mat) <- "integer"

stopifnot(all(colnames(counts_mat) == rownames(sample_sheet)))

## ---- 2. DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sample_sheet,
  design    = ~ group
)
dds <- dds[rowSums(counts(dds)) > 0, ]   # drop genes with zero counts everywhere
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)

## ---- 3. RPKM ------------------------------------------------------------
## RPKM = 10^9 * count / (gene_length_bp * total_mapped_reads_for_sample)
lib_size <- colSums(counts_mat)
gene_len_matched <- gene_length[rownames(counts_mat)]
rpkm_mat <- sweep(counts_mat / (gene_len_matched / 1000), 2, lib_size / 1e6, "/")

## per-group means (RPKM and normalized counts), used in the output table
group_mean <- function(mat, group_vec) {
  sapply(levels(group_vec), function(g) rowMeans(mat[, group_vec == g, drop = FALSE]))
}
rpkm_group_mean  <- group_mean(rpkm_mat, sample_sheet$group)
norm_group_mean  <- group_mean(norm_counts, sample_sheet$group)
colnames(rpkm_group_mean) <- paste0("RPKM_", colnames(rpkm_group_mean))
colnames(norm_group_mean) <- paste0("NormCount_", colnames(norm_group_mean))

## ---- 4. Gene annotation (Ensembl ID -> symbol + description) ---------------
ensembl_ids <- sub("\\..*$", "", rownames(counts_mat))  # strip version suffix
gene_symbol <- mapIds(org.Mm.eg.db, keys = ensembl_ids, keytype = "ENSEMBL",
                       column = "SYMBOL", multiVals = "first")
gene_desc   <- mapIds(org.Mm.eg.db, keys = ensembl_ids, keytype = "ENSEMBL",
                       column = "GENENAME", multiVals = "first")
annot <- data.frame(
  Geneid          = rownames(counts_mat),
  `Gene symbol`   = ifelse(is.na(gene_symbol), rownames(counts_mat), gene_symbol),
  `Gene description` = ifelse(is.na(gene_desc), NA, gene_desc),
  check.names = FALSE
)

## ---- 5. DEG cutoff function (tissue-dependent, configurable) ---------------
## Applies: RPKM > rpkm_cutoff in >=1 of the two compared groups, FC outside
## [fc_down, fc_up], and padj < alpha — falling back to raw pvalue < alpha if
## too few genes pass padj (you decide the "too few" threshold).
flag_degs <- function(res_df, group1, group2,
                       rpkm_cutoff = 1, fc_up = 1.5, fc_down = 0.6,
                       alpha = 0.05, min_degs_for_padj = 10) {
  rpkm_cols  <- paste0("RPKM_", c(group1, group2))
  passes_rpkm <- apply(res_df[, rpkm_cols, drop = FALSE], 1, function(x) any(x > rpkm_cutoff, na.rm = TRUE))
  passes_fc   <- (res_df$FC > fc_up) | (res_df$FC < fc_down)

  passes_padj <- !is.na(res_df$padj) & res_df$padj < alpha
  n_padj <- sum(passes_rpkm & passes_fc & passes_padj, na.rm = TRUE)

  if (n_padj >= min_degs_for_padj) {
    res_df$DEG        <- passes_rpkm & passes_fc & passes_padj
    cutoff_used <- "padj"
  } else {
    message(sprintf("[%s vs %s] only %d DEGs at padj<%.2f — falling back to pvalue<%.2f",
                     group1, group2, n_padj, alpha, alpha))
    passes_pval <- !is.na(res_df$pvalue) & res_df$pvalue < alpha
    res_df$DEG  <- passes_rpkm & passes_fc & passes_pval
    cutoff_used <- "pvalue"
  }
  attr(res_df, "cutoff_used") <- cutoff_used
  res_df
}

## ---- 6. Build one results table per comparison ------------------------------
## Edit this list for the comparisons you actually want.
comparisons <- list(
  c("FN", "MN"),   # sex effect, naive
  c("MT", "MN"),   # tumor effect, males
  c("FT", "FN"),   # tumor effect, females
  c("FT", "MT")    # sex effect, tumor-bearing
)

build_comparison_table <- function(dds, group1, group2) {
  res <- results(dds, contrast = c("group", group1, group2))
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)

  out <- annot %>%
    left_join(res_df, by = "Geneid") %>%
    left_join(as.data.frame(rpkm_group_mean) %>%
                mutate(Geneid = rownames(rpkm_group_mean)),
              by = "Geneid") %>%
    left_join(as.data.frame(norm_group_mean) %>%
                mutate(Geneid = rownames(norm_group_mean)),
              by = "Geneid")

  out$FC     <- 2 ^ out$log2FoldChange           # linear fold change, group1/group2
  out$LogFC  <- out$log2FoldChange
  out <- out %>%
    rename(Pval = pvalue, padj = padj) %>%
    select(`Gene symbol`, `Gene description`,
           starts_with("RPKM_"), starts_with("NormCount_"),
           FC, LogFC, Pval, padj)

  out <- flag_degs(out, group1, group2)
  attr(out, "cutoff_used") <- attr(out, "cutoff_used")
  out
}

results_list <- lapply(comparisons, function(cp) build_comparison_table(dds, cp[1], cp[2]))
names(results_list) <- sapply(comparisons, paste, collapse = "_vs_")

## ---- 7. Write Excel workbook -------------------------------------------------
wb <- createWorkbook()

summary_df <- data.frame(
  Comparison  = character(), Cutoff_used = character(),
  n_up = integer(), n_down = integer(), stringsAsFactors = FALSE
)

for (cn in names(results_list)) {
  tab <- results_list[[cn]]
  cutoff_used <- attr(tab, "cutoff_used")

  addWorksheet(wb, cn)
  writeData(wb, cn, tab)

  deg_sheet <- paste0(cn, "_DEGs")
  addWorksheet(wb, deg_sheet)
  writeData(wb, deg_sheet, tab %>% filter(DEG))

  summary_df <- rbind(summary_df, data.frame(
    Comparison  = cn,
    Cutoff_used = cutoff_used,
    n_up   = sum(tab$DEG & tab$FC > 1, na.rm = TRUE),
    n_down = sum(tab$DEG & tab$FC < 1, na.rm = TRUE)
  ))
}

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
worksheetOrder(wb) <- c(which(sheets(wb) == "Summary"),
                         setdiff(seq_along(sheets(wb)), which(sheets(wb) == "Summary")))

saveWorkbook(wb, "GSE220448_DE_results.xlsx", overwrite = TRUE)
message("Done -> GSE220448_DE_results.xlsx")

## ---- Cutoffs used (tongue tissue, per your lab's convention) ----------------
## RPKM > 1 in at least one compared group
## FC > 1.5 or FC < 0.6
## padj < 0.05, falling back to pvalue < 0.05 if <10 DEGs pass padj
## Adjust rpkm_cutoff / fc_up / fc_down / alpha in flag_degs() for other tissues.
