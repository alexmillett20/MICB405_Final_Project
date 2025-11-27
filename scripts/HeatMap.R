# Harper Rapkin
# MICB405 Generating HeatMap
# 3/11/2025

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

#read in the dds file
dds <- readRDS("./DEA_outputs/dds_no_ctrlrep8_IL13.rds")

# Perform log transformation on our count data
rld <- rlog(dds)

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Generate a heatmap using the pheatmap package
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         filename = "MICB405_IL13_HeatMap.png")

