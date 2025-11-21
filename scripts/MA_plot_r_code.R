library(DESeq2)
library(ggplot2)
library(ggpubr)

# load .rds to R
dds_control_13 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_no_ctrlrep8_IL13.rds")
dds_control_4 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_no_ctrlrep8_IL4.rds")
dds_13_4 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_IL13_IL4.rds")

# extract differential analysis results
res_control_13 <- results(dds_control_13)
res_df_control_13 <- as.data.frame(res_control_13)
res_df_control_13$gene_id <- rownames(res_df_control_13)

res_control_4 <- results(dds_control_4)
res_df_control_4 <- as.data.frame(res_control_4)
res_df_control_4$gene_id <- rownames(res_df_control_4)

res_13_4 <- results(dds_13_4)
res_df_13_4 <- as.data.frame(res_13_4)
res_df_13_4$gene_id <- rownames(res_df_13_4)

# Specify genes to label
up_genes_control_13 <- res_df_control_13 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

down_genes_control_13 <- res_df_control_13 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

up_genes_control_4 <- res_df_control_4 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

down_genes_control_4 <- res_df_control_4 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

up_genes_13_4 <- res_df_13_4 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

down_genes_13_4 <- res_df_13_4 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

# Plot
label_genes <- c(up_genes_control_13$gene_id, down_genes_control_13$gene_id)
ggmaplot(res_df_control_13, main = "MA Plot",
         fdr = 0.0001, fc.cut = 1, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df_control_13$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())

label_genes <- c(up_genes_control_4$gene_id, down_genes_control_4$gene_id)
ggmaplot(res_df_control_4, main = "MA Plot",
         fdr = 0.0001, fc.cut = 1, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df_control_4$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())

label_genes <- c(up_genes_13_4$gene_id, down_genes_13_4$gene_id)
ggmaplot(res_df_13_4, main = "MA Plot",
         fdr = 0.0001, fc.cut = 1, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df_13_4$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())
