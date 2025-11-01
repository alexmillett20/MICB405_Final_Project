suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap)) # package for generating heatmaps

dds <- readRDS("./DEA_outputs/dds.rds")

rld <- rlog(dds)

# Get the PCA data as a dataframe with returnData = TRUE
pcaData <- plotPCA(rld, intgroup=c("treatment"), returnData = TRUE) 

glimpse(pcaData)


# Using the attr() function, we will extract the percent variation explained of each axis of the pcaData object
percentVar <- round(100 * attr(pcaData, "percentVar"))
glimpse(percentVar)
