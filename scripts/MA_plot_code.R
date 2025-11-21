# load required packages
library(DESeq2)
library(ggplot2)
library(ggpubr)

# load .rds to R
dds_control_13 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_final_ctrl_IL13.rds")
dds_control_4 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_final_ctrl_IL4.rds")
dds_13_4 <- readRDS("C:/Users/jeong/OneDrive/Desktop/UBC/25-26/MICB 405/Final assignemet/dds_final_IL13_IL4.rds")

# extract differential analysis results
res_control_13 <- results(dds_control_13) # run statistical test
res_df_control_13 <- as.data.frame(res_control_13) # converts result to data frame
res_df_control_13$gene_id <- rownames(res_df_control_13) # adds gene_id column

res_control_4 <- results(dds_control_4)
res_df_control_4 <- as.data.frame(res_control_4)
res_df_control_4$gene_id <- rownames(res_df_control_4)

res_13_4 <- results(dds_13_4)
res_df_13_4 <- as.data.frame(res_13_4)
res_df_13_4$gene_id <- rownames(res_df_13_4)

# Specify genes to label
# filter upregulated gene by p value of 1e-5, save as variable
up_genes_control_13 <- res_df_control_13 %>% 
  dplyr::filter(log2FoldChange > 1, padj < 0.05) %>% 
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10) #limit number of labeled genes to top 10 up-regulated genes. without it, too many genes are labeled and not readable.

# filter downregulated gene by p value of 1e-5, save as variable
down_genes_control_13 <- res_df_control_13 %>% 
  dplyr::filter(log2FoldChange < -1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10) #limit number of labeled genes to top 10 down-regulated genes. without it, too many genes are labeled and not readable.

up_genes_control_4 <- res_df_control_4 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

down_genes_control_4 <- res_df_control_4 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

up_genes_13_4 <- res_df_13_4 %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

down_genes_13_4 <- res_df_13_4 %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n=10)

# Plot
label_genes <- c(up_genes_control_13$gene_id, down_genes_control_13$gene_id)
ggmaplot(res_df_control_13, main = "MA Plot",
         fdr = 0.05, fc.cut = 1, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df_control_13$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

label_genes <- c(up_genes_control_4$gene_id, down_genes_control_4$gene_id)
ggmaplot(res_df_control_4, main = "MA Plot",
         fdr = 0.05, fc.cut = 1, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df_control_4$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

label_genes <- c(up_genes_13_4$gene_id, down_genes_13_4$gene_id)
ggmaplot(res_df_13_4, main = "MA Plot",
         fdr = 0.05, fc.cut = 1, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df_13_4$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  guides(color = guide_legend(override.aes = list(size = 4)))
