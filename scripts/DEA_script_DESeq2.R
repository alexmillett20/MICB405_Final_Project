# Author: Alexandra Millett
# MICB405 Final Project DESeq2 Analysis
# Date created: October 31, 2025
# Last updated: November 13, 2025

# This script is to perform differential expression analysis for the MICB405 project.
# The structure of the script to perform DEA is based on tutorial 7. 

library(tidyverse)
library(DESeq2)


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


#### data wrangling --------------------------------------------
combined_data <- data.frame(row.names = control_rep_1$gene_id,
                  control_rep_1 = control_rep_1$total,
                  control_rep_3 = control_rep_3$total,
                  control_rep_4 = control_rep_4$total,
                  control_rep_6 = control_rep_6$total,
                  control_rep_7 = control_rep_7$total,
                  control_rep_8 = control_rep_8$total,
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
combined_data_matrix<- as.matrix(combined_data) 

# Metadata for samples
metadata <- data.frame(row.names = colnames(combined_data_matrix), 
                       treatment = c("control", "control", "control", "control", "control", "control", "control",
                                     "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4",
                                     "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13"
                                     )
)

metadata$label <- rownames(metadata)

colnames(combined_data_matrix) == rownames(metadata)

# create dds_matrix
dds_matrix <- DESeqDataSetFromMatrix(countData = combined_data_matrix,  
                                     colData = metadata, 
                                     design = ~treatment)


# Set control
dds_matrix$treatment <- relevel(dds_matrix$treatment, ref = "control")

dds <- DESeq(dds_matrix)

saveRDS(dds, "./DEA_outputs/dds_all_samples.rds")


# redo dds without control rep 8 vs IL13 -----------------------------
combined_data_no_ctrlrep8 <- data.frame(row.names = control_rep_1$gene_id,
                            control_rep_1 = control_rep_1$total,
                            control_rep_3 = control_rep_3$total,
                            control_rep_4 = control_rep_4$total,
                            control_rep_6 = control_rep_6$total,
                            control_rep_7 = control_rep_7$total,
                            control_rep_9 = control_rep_9$total,
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

# dds without control rep 8 compared to IL13 -----------------------------
combined_data_no_ctrlrep8_IL13 <- data.frame(row.names = control_rep_1$gene_id,
                                        control_rep_1 = control_rep_1$total,
                                        control_rep_3 = control_rep_3$total,
                                        control_rep_4 = control_rep_4$total,
                                        control_rep_6 = control_rep_6$total,
                                        control_rep_7 = control_rep_7$total,
                                        control_rep_9 = control_rep_9$total,
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
combined_data_matrix_no_ctrlrep8_IL13 <- as.matrix(combined_data_no_ctrlrep8_IL13) 

# Metadata for samples
metadata_no_ctrlrep8_IL13 <- data.frame(row.names = colnames(combined_data_matrix_no_ctrlrep8_IL13), 
                                   treatment = c("control", "control", "control", "control", "control", "control",
                                                 "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13"
                                   ))


metadata_no_ctrlrep8_IL13$label <- rownames(metadata_no_ctrlrep8_IL13)

colnames(combined_data_matrix_no_ctrlrep8_IL13) == rownames(metadata_no_ctrlrep8_IL13)

# create dds_matrix
dds_matrix_no_ctrlrep8_IL13 <- DESeqDataSetFromMatrix(countData = combined_data_matrix_no_ctrlrep8_IL13,  
                                                 colData = metadata_no_ctrlrep8_IL13, 
                                                 design = ~treatment)


# Set control
dds_matrix_no_ctrlrep8_IL13$treatment <- relevel(dds_matrix_no_ctrlrep8_IL13$treatment, ref = "control")


dds_no_ctrlrep8_IL13 <- DESeq(dds_matrix_no_ctrlrep8_IL13)
saveRDS(dds_no_ctrlrep8_IL13, "./DEA_outputs/dds_no_ctrlrep8_IL13.rds")


## dds without control rep 8 compared to IL4 -----------------------------
combined_data_no_ctrlrep8_IL4 <- data.frame(row.names = control_rep_1$gene_id,
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
                                             IL4_rep_9 = IL4_rep_9$total)


# transform to matrix
combined_data_matrix_no_ctrlrep8_IL4 <- as.matrix(combined_data_no_ctrlrep8_IL4) 

# Metadata for samples
metadata_no_ctrlrep8_IL4 <- data.frame(row.names = colnames(combined_data_matrix_no_ctrlrep8_IL4), 
                                        treatment = c("control", "control", "control", "control", "control", "control",
                                                      "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4"
                                        ))


metadata_no_ctrlrep8_IL4$label <- rownames(metadata_no_ctrlrep8_IL4)

colnames(combined_data_matrix_no_ctrlrep8_IL4) == rownames(metadata_no_ctrlrep8_IL4)

# create dds_matrix
dds_matrix_no_ctrlrep8_IL4 <- DESeqDataSetFromMatrix(countData = combined_data_matrix_no_ctrlrep8_IL4,  
                                                      colData = metadata_no_ctrlrep8_IL4, 
                                                      design = ~treatment)


# Set control
dds_matrix_no_ctrlrep8_IL4$treatment <- relevel(dds_matrix_no_ctrlrep8_IL4$treatment, ref = "control")


dds_no_ctrlrep8_IL4 <- DESeq(dds_matrix_no_ctrlrep8_IL4)
saveRDS(dds_no_ctrlrep8_IL4, "./DEA_outputs/dds_no_ctrlrep8_IL4.rds")

## dds IL13 compared to IL4 -----------------------------
combined_data_IL13_IL4 <- data.frame(row.names = IL13_rep_1$gene_id,
                                            IL13_rep_1 = IL13_rep_1$total,
                                            IL13_rep_2 = IL13_rep_2$total,
                                            IL13_rep_3 = IL13_rep_3$total,
                                            IL13_rep_4 = IL13_rep_4$total,
                                            IL13_rep_5 = IL13_rep_5$total,
                                            IL13_rep_6 = IL13_rep_6$total,
                                            IL13_rep_7 = IL13_rep_7$total,
                                            IL13_rep_8 = IL13_rep_8$total,
                                            IL13_rep_9 = IL13_rep_9$total,
                                            IL4_rep_2 = IL4_rep_2$total,
                                            IL4_rep_3 = IL4_rep_3$total,
                                            IL4_rep_4 = IL4_rep_4$total,
                                            IL4_rep_5 = IL4_rep_5$total,
                                            IL4_rep_6 = IL4_rep_6$total,
                                            IL4_rep_8 = IL4_rep_8$total,
                                            IL4_rep_9 = IL4_rep_9$total)


# transform to matrix
combined_data_matrix_IL13_IL4 <- as.matrix(combined_data_IL13_IL4) 

# Metadata for samples
metadata_IL13_IL4 <- data.frame(row.names = colnames(combined_data_matrix_IL13_IL4), 
                                       treatment = c("IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13",
                                                     "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4"
                                       ))


metadata_IL13_IL4$label <- rownames(metadata_IL13_IL4)

colnames(combined_data_matrix_IL13_IL4) == rownames(metadata_IL13_IL4)

# create dds_matrix
dds_matrix_IL13_IL4 <- DESeqDataSetFromMatrix(countData = combined_data_matrix_IL13_IL4,  
                                                     colData = metadata_IL13_IL4, 
                                                     design = ~treatment)


# Set control
dds_matrix_IL13_IL4$treatment <- relevel(dds_matrix_IL13_IL4$treatment, ref = "IL13")


dds_IL13_IL4 <- DESeq(dds_matrix_IL13_IL4)
saveRDS(dds_IL13_IL4, "./DEA_outputs/dds_IL13_IL4.rds")

################ FINALIZED DATA -------------------------------------
### Based on the heatmap results from the dds objects generated above, we will remove the following replicates: control_rep3, control_rep8, IL13_rep6
# redo dds control (no rep 3 and 8) vs IL4 -----------------------------
combined_data_final_ctrl_IL4 <- data.frame(row.names = control_rep_1$gene_id,
                                        control_rep_1 = control_rep_1$total,
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
                                        IL4_rep_9 = IL4_rep_9$total)

# transform to matrix
combined_data_final_matrix_ctrl_IL4 <- as.matrix(combined_data_final_ctrl_IL4) 

# Metadata for samples
metadata_final_ctrl_IL4 <- data.frame(row.names = colnames(combined_data_final_matrix_ctrl_IL4), 
                                   treatment = c("control", "control", "control", "control", "control",
                                                 "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4"))


metadata_final_ctrl_IL4$label <- rownames(metadata_final_ctrl_IL4)

colnames(combined_data_final_matrix_ctrl_IL4) == rownames(metadata_final_ctrl_IL4)

# create dds_matrix
dds_matrix_final_ctrl_IL4 <- DESeqDataSetFromMatrix(countData = combined_data_final_matrix_ctrl_IL4,  
                                                 colData = metadata_final_ctrl_IL4, 
                                                 design = ~treatment)


# Set control
dds_matrix_final_ctrl_IL4$treatment <- relevel(dds_matrix_final_ctrl_IL4$treatment, ref = "control")


dds_final_ctrl_IL4 <- DESeq(dds_matrix_final_ctrl_IL4)
saveRDS(dds_final_ctrl_IL4, "./DEA_outputs/dds_final_ctrl_IL4.rds")

# redo dds control (no rep 3 and 8) vs IL13 (no rep 6) -----------------------------
combined_data_final_ctrl_IL13 <- data.frame(row.names = control_rep_1$gene_id,
                                            control_rep_1 = control_rep_1$total,
                                            control_rep_4 = control_rep_4$total,
                                            control_rep_6 = control_rep_6$total,
                                            control_rep_7 = control_rep_7$total,
                                            control_rep_9 = control_rep_9$total,
                                            IL13_rep_1 = IL13_rep_1$total,
                                            IL13_rep_2 = IL13_rep_2$total,
                                            IL13_rep_3 = IL13_rep_3$total,
                                            IL13_rep_4 = IL13_rep_4$total,
                                            IL13_rep_5 = IL13_rep_5$total,
                                            IL13_rep_7 = IL13_rep_7$total,
                                            IL13_rep_8 = IL13_rep_8$total,
                                            IL13_rep_9 = IL13_rep_9$total
)

# transform to matrix
combined_data_final_matrix_ctrl_IL13 <- as.matrix(combined_data_final_ctrl_IL13) 

# Metadata for samples
metadata_final_ctrl_IL13 <- data.frame(row.names = colnames(combined_data_final_matrix_ctrl_IL13), 
                                       treatment = c("control", "control", "control", "control", "control",
                                                     "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13"))


metadata_final_ctrl_IL13$label <- rownames(metadata_final_ctrl_IL13)

colnames(combined_data_final_matrix_ctrl_IL13) == rownames(metadata_final_ctrl_IL13)

# create dds_matrix
dds_matrix_final_ctrl_IL13 <- DESeqDataSetFromMatrix(countData = combined_data_final_matrix_ctrl_IL13,  
                                                     colData = metadata_final_ctrl_IL13, 
                                                     design = ~treatment)


# Set control
dds_matrix_final_ctrl_IL13$treatment <- relevel(dds_matrix_final_ctrl_IL13$treatment, ref = "control")


dds_final_ctrl_IL13 <- DESeq(dds_matrix_final_ctrl_IL13)
saveRDS(dds_final_ctrl_IL13, "./DEA_outputs/dds_final_ctrl_IL13.rds")



## dds IL13 (no rep 6) compared to IL4 -----------------------------
combined_data_final_IL13_IL4 <- data.frame(row.names = IL13_rep_1$gene_id,
                                     IL13_rep_1 = IL13_rep_1$total,
                                     IL13_rep_2 = IL13_rep_2$total,
                                     IL13_rep_3 = IL13_rep_3$total,
                                     IL13_rep_4 = IL13_rep_4$total,
                                     IL13_rep_5 = IL13_rep_5$total,
                                     IL13_rep_7 = IL13_rep_7$total,
                                     IL13_rep_8 = IL13_rep_8$total,
                                     IL13_rep_9 = IL13_rep_9$total,
                                     IL4_rep_2 = IL4_rep_2$total,
                                     IL4_rep_3 = IL4_rep_3$total,
                                     IL4_rep_4 = IL4_rep_4$total,
                                     IL4_rep_5 = IL4_rep_5$total,
                                     IL4_rep_6 = IL4_rep_6$total,
                                     IL4_rep_8 = IL4_rep_8$total,
                                     IL4_rep_9 = IL4_rep_9$total)


# transform to matrix
combined_data_final_matrix_IL13_IL4 <- as.matrix(combined_data_final_IL13_IL4) 

# Metadata for samples
metadata_final_IL13_IL4 <- data.frame(row.names = colnames(combined_data_final_matrix_IL13_IL4), 
                                treatment = c("IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13", "IL13",
                                              "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4"
                                ))


metadata_final_IL13_IL4$label <- rownames(metadata_final_IL13_IL4)

colnames(combined_data_final_matrix_IL13_IL4) == rownames(metadata_final_IL13_IL4)

# create dds_matrix
dds_matrix_final_IL13_IL4 <- DESeqDataSetFromMatrix(countData = combined_data_final_matrix_IL13_IL4,  
                                              colData = metadata_final_IL13_IL4, 
                                              design = ~treatment)


# Set control
dds_matrix_final_IL13_IL4$treatment <- relevel(dds_matrix_final_IL13_IL4$treatment, ref = "IL13")


dds_final_IL13_IL4 <- DESeq(dds_matrix_final_IL13_IL4)
saveRDS(dds_final_IL13_IL4, "./DEA_outputs/dds_final_IL13_IL4.rds")


##### generating plots -----------------------------------------------
##### all data
# Perform log transformation on our count data
rld <- rlog(dds)

# Generate a PCA plot with DESeq2's plotPCA function
PCA_plot <- plotPCA(rld, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_all_samples.png", PCA_plot, width = 10, height = 10)


##### no control rep 8
# Perform log transformation on our count data
rld_no_ctrlrep8 <- rlog(dds_no_ctrlrep8)

# Generate a PCA plot with DESeq2's plotPCA function
PCA_plot_no_ctrlrep8 <- plotPCA(rld_no_ctrlrep8, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_no_ctrlrep8.png", PCA_plot_no_ctrlrep8, width = 10, height = 10)

##### no control rep 8 and IL13 only
# Perform log transformation on our count data
rld_no_ctrlrep8_IL13 <- rlog(dds_no_ctrlrep8_IL13)

# Generate a PCA plot with DESeq2's plotPCA function
PCA_plot_no_ctrlrep8_IL13 <- plotPCA(rld_no_ctrlrep8_IL13, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_no_ctrlrep8_IL13.png", PCA_plot_no_ctrlrep8_IL13, width = 10, height = 10)


##### no control rep 8 and IL4 only
rld_no_ctrlrep8_IL4 <- rlog(dds_no_ctrlrep8_IL4)

# Generate a PCA plot with DESeq2's plotPCA function
PCA_plot_no_ctrlrep8_IL4 <- plotPCA(rld_no_ctrlrep8_IL4, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_no_ctrlrep8_IL4.png", PCA_plot_no_ctrlrep8_IL4, width = 10, height = 10)

#### IL13 and IL4 only
rld_IL13_IL4 <- rlog(dds_IL13_IL4)

# Generate a PCA plot with DESeq2's plotPCA function
PCA_plot_IL13_IL4 <- plotPCA(rld_IL13_IL4, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_IL13_IL4.png", PCA_plot_IL13_IL4, width = 10, height = 10)


#################### FINALIZED PLOTS ------------------------

### Control (no rep 8 or 3) and IL13 (no rep 6)
rld_final_ctrl_IL13 <- rlog(dds_final_ctrl_IL13)
PCA_plot_final_ctrl_IL13 <- plotPCA(rld_final_ctrl_IL13, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_final_ctrl_IL13.png", PCA_plot_final_ctrl_IL13, width = 10, height = 10)

### Control (no rep 8 or 3) and IL4
rld_final_ctrl_IL4 <- rlog(dds_final_ctrl_IL4)
PCA_plot_final_ctrl_IL4 <- plotPCA(rld_final_ctrl_IL4, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_final_ctrl_IL4.png", PCA_plot_final_ctrl_IL4, width = 10, height = 10)

### IL13 (no rep 6) and IL4
rld_final_IL13_IL4 <- rlog(dds_final_IL13_IL4)
PCA_plot_final_IL13_IL4 <- plotPCA(rld_final_IL13_IL4, intgroup = "treatment") +
  geom_text(aes(label = label))
ggsave("./DEA_outputs/PCA_plot_final_IL13_IL4.png", PCA_plot_final_IL13_IL4, width = 10, height = 10)


##################### Generating filtered up and down regulated gene lists ----------------------------
# Load finalized dds objects
dds_final_ctrl_IL4 <- readRDS("./DEA_outputs/dds_final_ctrl_IL4.rds")
dds_final_ctrl_IL13 <- readRDS("./DEA_outputs/dds_final_ctrl_IL13.rds")
dds_final_IL13_IL4 <- readRDS("./DEA_outputs/dds_final_IL13_IL4.rds")

# Get results for control vs IL4
res_final_ctrl_IL4 <- results(dds_final_ctrl_IL4, name = "treatment_IL4_vs_control") %>% as.data.frame()

res_final_ctrl_IL4_no_NA <- res_final_ctrl_IL4 %>%
  drop_na()

res_filtered_ctrl_IL_4 <- res_final_ctrl_IL4_no_NA %>%
  filter(padj < 0.0001)

res_filtered_final_ctrl_IL_4 <- res_filtered_ctrl_IL_4 %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>%
  rownames_to_column("gene_id")

write_csv(res_filtered_final_ctrl_IL_4, "./DEA_outputs/filtered_DE_genes_final_ctrl_vs_IL4.csv")

res_filtered_up_final_ctrl_IL_4 <- res_filtered_final_ctrl_IL_4 %>%
  filter(log2FoldChange > 1)

write_csv(res_filtered_up_final_ctrl_IL_4, "./DEA_outputs/filtered_upregulated_genes_final_ctrl_vs_IL4.csv")

res_filtered_down_final_ctrl_IL_4 <- res_filtered_final_ctrl_IL_4 %>%
  filter(log2FoldChange < -1)

write_csv(res_filtered_down_final_ctrl_IL_4, "./DEA_outputs/filtered_downregulated_genes_final_ctrl_vs_IL4.csv")

# Get results for control vs IL13
res_final_ctrl_IL13 <- results(dds_final_ctrl_IL13, name = "treatment_IL13_vs_control") %>% as.data.frame()
res_final_ctrl_IL13_no_NA <- res_final_ctrl_IL13 %>%
  drop_na()

res_filtered_ctrl_IL_13 <- res_final_ctrl_IL13_no_NA %>%
  filter(padj < 0.0001)

res_filtered_final_ctrl_IL_13 <- res_filtered_ctrl_IL_13 %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>%
  rownames_to_column("gene_id")

write_csv(res_filtered_final_ctrl_IL_13, "./DEA_outputs/filtered_DE_genes_final_ctrl_vs_IL13.csv")
res_filtered_up_final_ctrl_IL_13 <- res_filtered_final_ctrl_IL_13 %>%
  filter(log2FoldChange > 1)

write_csv(res_filtered_up_final_ctrl_IL_13, "./DEA_outputs/filtered_upregulated_genes_final_ctrl_vs_IL13.csv")

res_filtered_down_final_ctrl_IL_13 <- res_filtered_final_ctrl_IL_13 %>%
  filter(log2FoldChange < -1)

write_csv(res_filtered_down_final_ctrl_IL_13, "./DEA_outputs/filtered_downregulated_genes_final_ctrl_vs_IL13.csv")

# Get results for IL13 vs IL4
res_final_IL13_IL4 <- results(dds_final_IL13_IL4, name = "treatment_IL4_vs_IL13") %>% as.data.frame()

res_final_IL13_IL4_no_NA <- res_final_IL13_IL4 %>%
  drop_na()

res_filtered_IL13_IL4 <- res_final_IL13_IL4_no_NA %>%
  filter(padj < 0.0001)

res_filtered_final_IL13_IL4 <- res_filtered_IL13_IL4 %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>%
  rownames_to_column("gene_id")

write_csv(res_filtered_final_IL13_IL4, "./DEA_outputs/filtered_DE_genes_final_IL13_vs_IL4.csv")

res_filtered_up_final_IL13_IL4 <- res_filtered_final_IL13_IL4 %>%
  filter(log2FoldChange > 1)

write_csv(res_filtered_up_final_IL13_IL4, "./DEA_outputs/filtered_upregulated_genes_final_IL13_vs_IL4.csv")

res_filtered_down_final_IL13_IL4 <- res_filtered_final_IL13_IL4 %>%
  filter(log2FoldChange < -1)

write_csv(res_filtered_down_final_IL13_IL4, "./DEA_outputs/filtered_downregulated_genes_final_IL13_vs_IL4.csv") 

