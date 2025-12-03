# MICB 405 Final project: MA Plot
# Author: Jeongan Kim
# Date created: Oct 27, 2025
# Last updated: Dec 3, 2025
# ChatGPT was used to troubleshoot the errors in the code

# This script is to plot DESeq2 object in MA plot using ggpubr package

## load required packages ##
library(DESeq2)
library(tidyverse)
library(ggpubr)

## load .rds to R ##
# control vs. IL-13
dds_control_13 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_final_ctrl_IL13.rds")
# control vs. IL-4
dds_control_4 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_final_ctrl_IL4.rds")
# IL-13 vs. IL-4
dds_13_4 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_final_IL13_IL4.rds")

## extract differential analysis results ##
# control vs. IL-13
res_control_13 <- results(dds_control_13) # run statistical test
res_df_control_13 <- as.data.frame(res_control_13) # converts result to data frame
res_df_control_13$gene_id <- rownames(res_df_control_13) # adds gene_id column

# control vs. IL-4
res_control_4 <- results(dds_control_4)
res_df_control_4 <- as.data.frame(res_control_4)
res_df_control_4$gene_id <- rownames(res_df_control_4)

# IL-13 vs. IL-4
res_13_4 <- results(dds_13_4)
res_df_13_4 <- as.data.frame(res_13_4)
res_df_13_4$gene_id <- rownames(res_df_13_4)

## specify genes to label ##
# control vs. IL-13: filter upregulated gene by p value of 1e-5, save as variable
up_genes_control_13 <- res_df_control_13 %>% 
  dplyr::filter(log2FoldChange > 1, padj < 0.05) %>% 
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10) #limit number of labeled genes to top 10 up-regulated genes. without it, too many genes are labeled and not readable.

# control vs. IL-13: filter downregulated gene by p value of 1e-5, save as variable
down_genes_control_13 <- res_df_control_13 %>% 
  dplyr::filter(log2FoldChange < -1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10) #limit number of labeled genes to top 10 down-regulated genes. without it, too many genes are labeled and not readable.

# control vs. IL-4
up_genes_control_4 <- res_df_control_4 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

# control vs. IL-4
down_genes_control_4 <- res_df_control_4 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

# IL-13 vs. IL-4
up_genes_13_4 <- res_df_13_4 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

# IL-13 vs. IL-4
down_genes_13_4 <- res_df_13_4 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

## Plot ##
# control vs. IL-13
label_genes_1 <- c(up_genes_control_13$gene_id, down_genes_control_13$gene_id)
ggmaplot(res_df_control_13, main = "MA Plot: Control vs. IL13",
         fdr = 0.05, fc.cut = 1, size = 1.5,
         palette = c("red", "blue", "grey"),
         genenames = as.vector(res_df_control_13$gene_id),
         label.select = label_genes_1,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

# control vs. IL-4
label_genes_2 <- c(up_genes_control_4$gene_id, down_genes_control_4$gene_id)
ggmaplot(res_df_control_4, main = "MA Plot: Control vs. IL4",
         fdr = 0.05, fc.cut = 1, size = 1.5,
         palette = c("red", "blue", "grey"),
         genenames = as.vector(res_df_control_4$gene_id),
         label.select = label_genes_2,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

# IL-13 vs. IL-4
label_genes_3 <- c(up_genes_13_4$gene_id, down_genes_13_4$gene_id)
ggmaplot(res_df_13_4, main = "MA Plot: IL13 vs. IL4",
         fdr = 0.05, fc.cut = 1, size = 1.5,
         palette = c("red", "blue", "grey"),
         genenames = as.vector(res_df_13_4$gene_id),
         label.select = label_genes_3,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  guides(color = guide_legend(override.aes = list(size = 4)))
