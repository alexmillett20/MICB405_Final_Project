# MICB405 Final Project Volcano Plot
# Date created: Parsa Nayyara, Nov 3, 2025
# Last updated: Parsa Nayyara, Nov 3, 2025

# This script is to plot DESeq2 object in volcano plots with EnhancedVolcano

# Install EnhancedVolcano
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")

# Import packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(EnhancedVolcano))

#Load DESeq2 object
dds <- readRDS("../DEA_outputs/dds_no_ctrlrep8_IL13.rds")

#extract results
res <- results(dds,
               contrast = c("treatment", "control", "IL13"))
#res <- lfcShrink(dds,
#                 contrast = c("treatment", "control", "IL13"),
#                 res=res, type = 'normal')

#plot most basic volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj')

