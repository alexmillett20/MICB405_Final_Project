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
dds_c_13 <- readRDS("../DEA_outputs/dds_no_ctrlrep8_IL13.rds")
#dds_c_4 <- readRDS("../DEA_outputs/dds_no_ctrlrep8_IL13.rds")
#dds_4_13 <- readRDS("../DEA_outputs/dds_no_ctrlrep8_IL13.rds")

#extract results
res_c_13 <- results(dds_c_13,
               contrast = c("treatment", "IL13", "control"))
#res_c_13 <- lfcShrink(dds_c_13,
#                 contrast = c("treatment", "IL13", "control"),
#                 res=res, type = 'normal')

res_c_4 <- results(dds_c_4,
                    contrast = c("treatment", "IL4", "control"))
#res_c_4 <- lfcShrink(dds_c_4,
#                 contrast = c("treatment", "IL4", "control"),
#                 res=res, type = 'normal')

res_4_13 <- results(dds_4_13,
                    contrast = c("treatment", "IL13", "IL4"))
#res_4_13 <- lfcShrink(dds_4_13,
#                 contrast = c("treatment", "IL13", "IL4"),
#                 res=res, type = 'normal')

#plot most basic volcano plot
#EnhancedVolcano(res_c_13,
#                lab = rownames(res_c_13),
#                x = 'log2FoldChange',
#                y = 'padj')

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change

keyvals.colour <- ifelse(
  res_c_13$log2FoldChange < -1.0, 'red',
  ifelse(res_c_13$log2FoldChange > 1.0, 'green',
         'grey'))
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == 'green'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low'

EnhancedVolcano(res_c_13,
                lab = rownames(res_c_13),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = rownames(res_c_13)[which(names(keyvals.colour) %in% c('High', 'Low'))],
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

EnhancedVolcano(res_c_4,
                lab = rownames(res_c_4),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = rownames(res_c_4)[which(names(keyvals.colour) %in% c('High', 'Low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'Control vs IL-4',
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

EnhancedVolcano(res_4_13,
                lab = rownames(res_4_13),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = rownames(res_4_13)[which(names(keyvals.colour) %in% c('High', 'Low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'IL-4 vs IL-13',
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

