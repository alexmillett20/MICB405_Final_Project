# MA plot using dds_no_ctrolrep8_IL13.rds
# The MA plot can be used to display a difference in expression 
# patterns between two samples. The horizontal axis (A) shows 
# the average intensity while the vertical axis (M) shows 
# the intensity ratio between the two samples for the same 
# data point. In essence, an MA plot is a scatter plot tilted to 
# the side so that the differentially expressed probe(sets)/genes 
# are located above or below the 0 value of M. An MA plot is 
# also useful to visualize the results of normalization where 
# you would hope to see the median of the values follow a 
# horizontal line.

library(DESeq2)

# load .rds to R
dds <- readRDS("./DEA_outputs/dds_no_ctrlrep8_IL13.rds")

# plot
plotMA(dds,main="MA PLOT", alpha=0.0001, ylim=c(-5, 5))

############################################################
install.packages("ggpubr")
library(ggplot2)
library(ggpubr)
res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

#current
ggmaplot(res_df, main = "MA Plot",
         fdr = 0.0001, fc.cut = 2, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df$gene_id),
         top = 15,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())

#attempt 3:
up_genes <- res_df %>%
  dplyr::filter(log2FoldChange > 1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

down_genes <- res_df %>%
  dplyr::filter(log2FoldChange < -1, padj < 0.0001) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 15)

label_genes <- c(up_genes$gene_id, down_genes$gene_id)
ggmaplot(res_df, main = "MA Plot",
         fdr = 0.0001, fc.cut = 1, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_df$gene_id),
         label.select = label_genes,
         xlab = "Log2 mean expression",
         ylab = "Log2 fold change",
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic())
