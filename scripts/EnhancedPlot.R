# MICB405 Final Project Volcano Plot
# Date created: Parsa Nayyara, Nov 5, 2025
# Last updated: Parsa Nayyara, Nov 5, 2025

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
               contrast = c("treatment", "IL13", "control"))
#res <- lfcShrink(dds,
#                 contrast = c("treatment", "IL13", "control"),
#                 res=res, type = 'normal')

#plot most basic volcano plot
#EnhancedVolcano(res,
#                lab = rownames(res),
#                x = 'log2FoldChange',
#                y = 'padj')

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change

keyvals.colour <- ifelse(
  res$log2FoldChange < -1.0, 'red',
  ifelse(res$log2FoldChange > 1.0, 'green',
         'grey'))
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == 'green'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low'

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = rownames(res)[which(names(keyvals.colour) %in% c('High', 'Low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'Control vs IL-13',
                pCutoff = 10e-5,
                FCcutoff = 1.0,
                pointSize = 1.5,
                labSize = 3.0,
                colCustom = keyvals.colour,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1.0,
                borderColour = 'black')

