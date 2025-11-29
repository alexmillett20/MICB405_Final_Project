# KEGG pathway analysis
# Author: Alexandra Millett
# Date created: Nov 27, 2025
# Last updated: Nov 29, 2025

# This script is to look at the pathways of interest (Glycolysis and EMP pathways) to see whether any up or downregulated DEGs
# map onto these pathways. 

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

# gene names
gene_names <- DEGs_IL13_IL4$gene_id
# Log2FC values
FC <- DEGs_IL13_IL4$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# add names to FC data
names(FC) <- gene_entrez_ids$ENTREZID

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

# gene names
gene_names <- DEGs_ctrl_IL13$gene_id
# Log2FC values
FC <- DEGs_ctrl_IL13$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# 0.36% of input gene IDs failed to map. 

# check which gene did not map
setdiff(gene_names, gene_entrez_ids$SYMBOL)
# Gene = LOC118567478/0610040J01Rik

# add names to FC data
names(FC) <- gene_entrez_ids$ENTREZID

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

# gene names
gene_names <- DEGs_ctrl_IL4$gene_id
# Log2FC values
FC <- DEGs_ctrl_IL4$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# add names to FC data
names(FC) <- gene_entrez_ids$ENTREZID

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
