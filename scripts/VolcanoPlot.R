suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap)) # package for generating heatmaps

dds <- readRDS("./DEA_outputs/dds.rds")

rld <- rlog(dds)

# Get the PCA data as a dataframe with returnData = TRUE
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData = TRUE) 