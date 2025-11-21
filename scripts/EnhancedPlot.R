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
dds_c_4 <- readRDS("../DEA_outputs/dds_no_ctrlrep8_IL4.rds")
dds_4_13 <- readRDS("../DEA_outputs/dds_IL13_IL4.rds")

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

# create custom key-value pairs for 'high', 'low', non-significant expression by fold-change

keyvals13.colour <- ifelse(
  res_c_13$log2FoldChange < -1.0, 'red',
  ifelse(res_c_13$log2FoldChange > 1.0, 'green',
         'grey'))
keyvals13.colour[is.na(keyvals13.colour)] <- 'grey'
names(keyvals13.colour)[keyvals13.colour == 'green'] <- 'upregulated'
names(keyvals13.colour)[keyvals13.colour == 'grey'] <- 'non-significant'
names(keyvals13.colour)[keyvals13.colour == 'red'] <- 'downregulated'

keyvals4.colour <- ifelse(
  res_c_4$log2FoldChange < -1.0, 'red',
  ifelse(res_c_4$log2FoldChange > 1.0, 'green',
         'grey'))

keyvals4.colour[is.na(keyvals4.colour)] <- 'grey'
names(keyvals4.colour)[keyvals4.colour == 'green'] <- 'upregulated'
names(keyvals4.colour)[keyvals4.colour == 'grey'] <- 'non-significant'
names(keyvals4.colour)[keyvals4.colour == 'red'] <- 'downregulated'

keyvals413.colour <- ifelse(
  res_4_13$log2FoldChange < -1.0, 'red',
  ifelse(res_4_13$log2FoldChange > 1.0, 'green',
         'grey'))
keyvals413.colour[is.na(keyvals413.colour)] <- 'grey'
names(keyvals413.colour)[keyvals413.colour == 'green'] <- 'upregulated'
names(keyvals413.colour)[keyvals413.colour == 'grey'] <- 'non-significant'
names(keyvals413.colour)[keyvals413.colour == 'red'] <- 'downregulated'

# only label selected genes
picked_genes = list('Zfp174', 'Ccl24', 'Prps1', 'Gm37510', 'Ccr7',
                    'Clec7a', 'mt-Cytb', 'mt-Nd5', 'mt-Nd4', 'mt-Atp6',
                    'Gm28437', 'mt-Co1', 'mt-Nd1', 'Enpep', 'Bnc2',
                    'Brnp5', 'GaInt16', 'Vit', 'Has1', 'Aox3')
  
EnhancedVolcano(res_c_13,
                lab = rownames(res_c_13),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = picked_genes,
                #selectLab = rownames(res_c_13)[which(names(keyvals.colour) %in% c('High', 'Low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'Control vs IL-13',
                pCutoff = 10e-5,
                FCcutoff = 1.0,
                pointSize = 1.5,
                labSize = 3.0,
                colCustom = keyvals13.colour,
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
                selectLab = picked_genes,
                #selectLab = rownames(res_c_4)[which(names(keyvals.colour) %in% c('High', 'Low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'Control vs IL-4',
                pCutoff = 10e-5,
                FCcutoff = 1.0,
                pointSize = 1.5,
                labSize = 3.0,
                colCustom = keyvals4.colour,
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
                selectLab = picked_genes,
                #selectLab = rownames(res_4_13)[which(names(keyvals.colour) %in% c('High', 'Low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'IL-4 vs IL-13',
                pCutoff = 10e-5,
                FCcutoff = 1.0,
                pointSize = 1.5,
                labSize = 3.0,
                colCustom = keyvals413.colour,
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

