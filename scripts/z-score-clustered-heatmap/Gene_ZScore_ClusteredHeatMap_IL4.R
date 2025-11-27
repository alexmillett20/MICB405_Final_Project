# MICB405 Final Project
# Created by: Harper Rapkin
# Last edited: 21/11/2025

# Load the required library
library(pheatmap)
library(DESeq2) # Assuming dds is a DESeqDataSet object

dds <- readRDS("./DEA_outputs/dds_final_ctrl_IL4.rds")

norm_counts <- counts(dds, normalized = TRUE)

res <- results(dds)

# Order by padj
res_ordered <- res[order(res$padj), ]

# Pick top 50 genes
top_genes <- rownames(res_ordered)[1:20]

mat <- norm_counts[top_genes, ]

# Z-score by gene (row)
mat_z <- t(scale(t(mat)))

sample_names <- colnames(mat_z)
groups <- gsub("_rep_.*", "", sample_names)

sort_df <- data.frame(sample = sample_names, group = groups)

desired_order_indices <- order(sort_df$group, sort_df$sample)

reordered_samples <- sample_names[desired_order_indices]


mat_z <- mat_z[, reordered_samples]

sample_groups <- data.frame(
  Condition = gsub("_rep_.*","", reordered_samples)
)
rownames(sample_groups) = reordered_samples

pheatmap(
  mat_z,
  annotation_col = sample_groups,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "Top DE Genes (Z-score standardized)",
  filename = "MICB405_IL4_Zscore_GeneHeatmap.png"
)


