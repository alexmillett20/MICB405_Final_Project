# KEGG pathway analysis
# Author: Alexandra Millett
# Date created: Nov 27, 2025
# Last updated: Dec 2, 2025

# This script is to look at the pathways of interest (Glycolysis and EMP pathways) to see whether any up or downregulated DEGs
# map onto these pathways. 
# ChatGPT was used to help with coding of this script to visualise the KEGG pathways of interest

library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(pathview)

getwd() # must be top level of repository
output_path <- "./plots/kegg_pathway_plots" # relative to top level of repo
working_path <- "../.." # relative to output path (top level of repository)

##### IL13 vs IL4 -------------------------
# DEG results
DEGs_IL13_IL4 <- read.csv("./DEA_outputs/filtered_DE_genes_final_IL13_vs_IL4.csv")
# make symbol column from gene ID
DEGs_IL13_IL4$SYMBOL <- DEGs_IL13_IL4$gene_id
# gene names
gene_names <- DEGs_IL13_IL4$gene_id
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# add entrez IDs to merged df
merged_df <- merge(DEGs_IL13_IL4, gene_entrez_ids, by = "SYMBOL")
# Log2FC values
FC <- merged_df$log2FoldChange
# add names to FC data
names(FC) <- merged_df$ENTREZID

setwd(output_path)
# glycolysis
path <- pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00010", # pathway ID found in KEGG database
         species = "mmu",
         out.suffix = "IL13_IL4_glycolysis")
path$plot.data.gene
#Aldoc

# Oxidative phosphorylation
path <- pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00190", # pathway ID found in KEGG database
         species = "mmu",
         out.suffix = "IL13_IL4_OxPhos")



##### control vs IL13 -------------------------
setwd(working_path)
# DEG results
DEGs_ctrl_IL13 <- read.csv("./DEA_outputs/filtered_DE_genes_final_ctrl_vs_IL13.csv")
# make symbol column from gene ID
DEGs_ctrl_IL13$SYMBOL <- DEGs_ctrl_IL13$gene_id
# gene names
gene_names <- DEGs_ctrl_IL13$SYMBOL
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# 0.36% of input gene IDs failed to map. (2 genes)
# add entrez IDs to merged df
merged_df <- merge(DEGs_ctrl_IL13, gene_entrez_ids, by = "SYMBOL")

# check which gene did not map
setdiff(gene_names, merged_df$SYMBOL)
# Genes = LOC118567478, 0610040J01Rik

# Log2FC values
FC <- merged_df$log2FoldChange

# add names to FC data
names(FC) <- merged_df$ENTREZID

setwd(output_path)
# glycolysis
pathview(gene.data = FC,
          gene.idtype = "entrez",
          pathway.id = "mmu00010",
          species = "mmu",
          out.suffix = "ctrl_IL13_glycolysis")

# Oxidative phosphorylation
path <- pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00190",
         species = "mmu",
         out.suffix = "ctrl_IL13_OxPhos")

path$plot.data.gene
# Atp6v0d2


##### control vs IL4 -------------------------
setwd(working_path)
# DEG results
DEGs_ctrl_IL4 <- read.csv("./DEA_outputs/filtered_DE_genes_final_ctrl_vs_IL4.csv")
# make symbol column from gene ID
DEGs_ctrl_IL4$SYMBOL <- DEGs_ctrl_IL4$gene_id
# gene names
gene_names <- DEGs_ctrl_IL4$SYMBOL
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# add entrez IDs to merged df
merged_df <- merge(DEGs_ctrl_IL4, gene_entrez_ids, by = "SYMBOL")

# Log2FC values
FC <- merged_df$log2FoldChange

# add names to FC data
names(FC) <- merged_df$ENTREZID

setwd(output_path)
# glycolysis
path <- pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00010",
         species = "mmu",
         out.suffix = "ctrl_IL4_glycolysis")

# Oxidative phosphorylation
path <- pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00190",
         species = "mmu",
         out.suffix = "ctrl_IL4_OxPhos")

path$plot.data.gene
# Atp6v0d2

setwd(working_path)
