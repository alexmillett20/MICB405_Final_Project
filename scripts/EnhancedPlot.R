# MICB405 Final Project Volcano Plot
# Date created: Parsa Nayyara, Nov 5, 2025
# Last updated: Parsa Nayyara, Dec 3, 2025

# This script is to plot DESeq2 object in volcano plots with EnhancedVolcano

# Import packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(EnhancedVolcano))

#Load DESeq2 object
dds_c_13 <- readRDS("./DEA_outputs/dds_final_ctrl_IL13.rds")
dds_c_4 <- readRDS("./DEA_outputs/dds_final_ctrl_IL4.rds")
dds_4_13 <- readRDS("./DEA_outputs/dds_final_IL13_IL4.rds")

#extract results
res_c_13 <- results(dds_c_13,
               contrast = c("treatment", "IL13", "control"))

res_c_4 <- results(dds_c_4,
                    contrast = c("treatment", "IL4", "control"))

res_4_13 <- results(dds_4_13,
                    contrast = c("treatment", "IL13", "IL4"))

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
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'Control vs IL-13',
                pCutoff = 0.05,
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
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'Control vs IL-4',
                pCutoff = 0.05,
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
                xlab = bquote(~Log[2]~ 'fold change'),
                axisLabSize = 10,
                title = 'IL-13 vs IL-4',
                pCutoff = 0.05,
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


# count upregulated and downregulated genes

print("Downregulated (control vs IL4)")
sum(res_c_4$log2FoldChange < -1 & res_c_4$padj < 0.05, na.rm=TRUE)

print("Upregulated (control vs IL4)")
sum(res_c_4$log2FoldChange > 1 & res_c_4$padj < 0.05, na.rm=TRUE)

print("Downregulated (control vs IL13)")
sum(res_c_13$log2FoldChange < -1 & res_c_13$padj < 0.05, na.rm=TRUE)

print("Upregulated (control vs IL13)")
sum(res_c_13$log2FoldChange > 1 & res_c_13$padj < 0.05, na.rm=TRUE)

print("Downregulated (IL13 vs IL4)")
sum(res_4_13$log2FoldChange < -1 & res_4_13$padj < 0.05, na.rm=TRUE)

print("Upregulated (IL13 vs IL4)")
sum(res_4_13$log2FoldChange > 1 & res_4_13$padj < 0.05, na.rm=TRUE)
