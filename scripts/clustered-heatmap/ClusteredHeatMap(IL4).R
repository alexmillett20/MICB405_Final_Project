# MICB405 Final Project
# Generate a clustered heat map for IL-4 compared to the control
# Created by: Harper Rapkin
# Last edited: Nov 28
# ChatGPT was used to cluster the groups on the final heatmap

# Load the required library
library(pheatmap)
library(DESeq2) 

# read in the dds file
dds <- readRDS("./DEA_outputs/dds_no_ctrlrep8_IL4.rds")

# Perform log transformation on our count data
rld <- rlog(dds)

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# 1. Restore the sample names from the rld object's columns
sample_names <- colnames(rld)
rownames(sample_dist_matrix) <- sample_names
colnames(sample_dist_matrix) <- sample_names 

# 2. Define the desired order
groups <- gsub("_rep_.*", "", sample_names)

# Create a data frame for sorting
sort_df <- data.frame(
  sample = sample_names,
  group = groups
)

# Sort first by group, then by sample name 
# The order function returns the indices to achieve the sorted order.
desired_order_indices <- order(sort_df$group, sort_df$sample)

# Get the list of reordered sample names
reordered_samples <- sample_names[desired_order_indices]

# 3. Apply the desired order to the distance matrix
reordered_dist_matrix <- sample_dist_matrix[reordered_samples, reordered_samples]

# Remove column names from the reordered matrix, 
colnames(reordered_dist_matrix) <- NULL

# Generate a heatmap using the pheatmap package.
# IMPORTANT: Set cluster_rows and cluster_cols to **FALSE** to force the defined order.
pheatmap(reordered_dist_matrix,
         # Set clustering to FALSE to force the order defined by the matrix input
         clustering_distance_rows = sample_dists, # You can keep the distance metrics, but clustering must be disabled
         clustering_distance_cols = sample_dists, # You can keep the distance metrics, but clustering must be disabled
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Heatmap with Forced IL4/Control Ordering",
         filename = "MICB405_IL4_HeatMap_Reordered.png")

