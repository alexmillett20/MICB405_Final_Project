# Author: Alexandra Millett
# MICB405 Final Project DESeq2 Analysis
# Date created: October 31, 2025
# Last updated: October 31, 2025

# This script is to perform differential expression analysis for the MICB405 project.
# The structure of the script to perform DEA is based on tutorial 7. 

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


# Set path to a directory with files and load files into R, can do setwd("YourPath/")
path <- getwd()

##### assigning sample read counts to variables ----------------------------
# Control replicate 1
control_rep_1 <- read_tsv("reads_per_gene/control_rep_1_ReadsPerGene.out.tab", 
                           col_names = c("gene_id", "total", "antisense", "sense"), 
                           skip = 4) 
# Control replicate 3
control_rep_3 <- read_tsv("reads_per_gene/control_rep_3_ReadsPerGene.out.tab",
                           col_names = c("gene_id", "total","antisense", "sense"),
                           skip = 4)

# Control replicate 4
control_rep_4 <- read_tsv("reads_per_gene/control_rep_4_ReadsPerGene.out.tab",
                         col_names = c("gene_id", "total","antisense", "sense"),
                         skip = 4)

# Control replicate 6
control_rep_6 <- read_tsv("reads_per_gene/control_rep_6_ReadsPerGene.out.tab",
                         col_names = c("gene_id", "total","antisense", "sense"),
                         skip = 4)

# Control replicate 7
control_rep_7 <- read_tsv("reads_per_gene/control_rep_7_ReadsPerGene.out.tab",
                         col_names = c("gene_id", "total","antisense", "sense"),
                         skip = 4)

# Control replicate 8
control_rep_8 <- read_tsv("reads_per_gene/control_rep_8_ReadsPerGene.out.tab",
                         col_names = c("gene_id", "total","antisense", "sense"),
                         skip = 4)

# Control replicate 9
control_rep_9 <- read_tsv("reads_per_gene/control_rep_9_ReadsPerGene.out.tab",
                         col_names = c("gene_id", "total","antisense", "sense"),
                         skip = 4)

# IL4 replicate 2
IL4_rep_2 <- read_tsv("reads_per_gene/IL-4_rep_2_ReadsPerGene.out.tab",
                         col_names = c("gene_id", "total","antisense", "sense"),
                         skip = 4)

# IL4 replicate 3
IL4_rep_3 <- read_tsv("reads_per_gene/IL-4_rep_3_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)
# IL4 replicate 4
IL4_rep_4 <- read_tsv("reads_per_gene/IL-4_rep_4_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)

# IL4 replicate 5
IL4_rep_5 <- read_tsv("reads_per_gene/IL-4_rep_5_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)

# IL4 replicate 6
IL4_rep_6 <- read_tsv("reads_per_gene/IL-4_rep_6_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)

# IL4 replicate 8
IL4_rep_8 <- read_tsv("reads_per_gene/IL-4_rep_8_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)

# IL4 replicate 9
IL4_rep_9 <- read_tsv("reads_per_gene/IL-4_rep_9_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)

# IL13 replicate 1
IL13_rep_1 <- read_tsv("reads_per_gene/IL-13_rep_1_ReadsPerGene.out.tab",
                     col_names = c("gene_id", "total","antisense", "sense"),
                     skip = 4)

# IL13 replicate 2
IL13_rep_2 <- read_tsv("reads_per_gene/IL-13_rep_2_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 3
IL13_rep_3 <- read_tsv("reads_per_gene/IL-13_rep_3_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 4
IL13_rep_4 <- read_tsv("reads_per_gene/IL-13_rep_4_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 5
IL13_rep_5 <- read_tsv("reads_per_gene/IL-13_rep_5_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 6
IL13_rep_6 <- read_tsv("reads_per_gene/IL-13_rep_6_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 7
IL13_rep_7 <- read_tsv("reads_per_gene/IL-13_rep_7_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 8
IL13_rep_8 <- read_tsv("reads_per_gene/IL-13_rep_8_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)

# IL13 replicate 9
IL13_rep_9 <- read_tsv("reads_per_gene/IL-13_rep_9_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)


# redo dds without control rep 8 -----------------------------
combined_data_no_ctrlrep8 <- data.frame(row.names = control_rep_1$gene_id,
                            control_rep_1 = control_rep_1$total,
                            control_rep_3 = control_rep_3$total,
                            control_rep_4 = control_rep_4$total,
                            control_rep_6 = control_rep_6$total,
                            control_rep_7 = control_rep_7$total,
                            control_rep_9 = control_rep_9$total,
                            IL4_rep_2 = IL4_rep_2$total,
                            IL4_rep_3 = IL4_rep_3$total,
                            IL4_rep_4 = IL4_rep_4$total,
                            IL4_rep_5 = IL4_rep_5$total,
                            IL4_rep_6 = IL4_rep_6$total,
                            IL4_rep_8 = IL4_rep_8$total,
                            IL4_rep_9 = IL4_rep_9$total,
                            IL13_rep_1 = IL13_rep_1$total,
                            IL13_rep_2 = IL13_rep_2$total,
                            IL13_rep_3 = IL13_rep_3$total,
                            IL13_rep_4 = IL13_rep_4$total,
                            IL13_rep_5 = IL13_rep_5$total,
                            IL13_rep_6 = IL13_rep_6$total,
                            IL13_rep_7 = IL13_rep_7$total,
                            IL13_rep_8 = IL13_rep_8$total,
                            IL13_rep_9 = IL13_rep_9$total
)

# transform to matrix
combined_data_matrix_no_ctrlrep8 <- as.matrix(combined_data_no_ctrlrep8) 

# Metadata for samples
metadata_no_ctrlrep8 <- data.frame(row.names = colnames(combined_data_matrix_no_ctrlrep8), 
                       treatment = c("control", "control", "control", "control", "control", "control",
                                     "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4",
                                     "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13"
                       )
)


metadata_no_ctrlrep8$label <- rownames(metadata_no_ctrlrep8)

colnames(combined_data_matrix_no_ctrlrep8) == rownames(metadata_no_ctrlrep8)

# create dds_matrix
dds_matrix_no_ctrlrep8 <- DESeqDataSetFromMatrix(countData = combined_data_matrix_no_ctrlrep8,  
                                     colData = metadata_no_ctrlrep8, 
                                     design = ~treatment)


# Set control
dds_matrix_no_ctrlrep8$treatment <- relevel(dds_matrix_no_ctrlrep8$treatment, ref = "control")


dds_no_ctrlrep8 <- DESeq(dds_matrix_no_ctrlrep8)
saveRDS(dds_no_ctrlrep8, "./DEA_outputs/dds_no_ctrlrep8.rds")




##### generating plots -----------------------------------------------
##### all data

##### no control rep 8
# Perform log transformation on our count data
rld_no_ctrlrep8 <- rlog(dds_no_ctrlrep8)

# Generate a PCA plot with DESeq2's plotPCA function
PCA_plot_no_ctrlrep8 <- plotPCA(rld_no_ctrlrep8, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_no_ctrlrep8.png", PCA_plot_no_ctrlrep8, width = 10, height = 10)



# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)

# Names of the results that DESeq2 calculated
resultsNames(dds)

# Now we will extract the results for our comparison between the 12h timepoint and the 1h timepoint
res <- results(dds, name = "condition_12_hour_vs_1_hour") %>% as.data.frame() # we save it as a dataframe for easy manipulation with dplyr 
head(res)

results(dds, contrast = c("condition", "12h", "6h"))

glimpse(res)

view(res)

res_no_NA <- res %>% 
  drop_na()

# How many rows did we filter out!?
glimpse(res_no_NA) 

res_filtered <- res_no_NA %>% 
  filter(padj <= 0.05)

# How many rows did we filter out!?
glimpse(res_filtered)


res_filtered_final <- res_filtered %>% 
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>% # the '|' stand for OR here!
  rownames_to_column("gene_id") # Convert the rownames into a column so they can be saved in your CSV file

head(res_filtered_final)

# Top 10 upregulated genes (most positive log2FoldChange)
top10_genes <- res_filtered_final %>%
  arrange(desc(log2FoldChange)) %>% # NOTE that we use the desc() function to organize the column in descending order
  head(n = 10)
top10_genes

bot10_genes <- res_filtered_final %>%
  arrange(log2FoldChange) %>% # NOTE since we don't use desc(), the column is organized in ascending order
  head(n = 10)
bot10_genes

write_csv(res_filtered_final, "chlamy_results.csv")
