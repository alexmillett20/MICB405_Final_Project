# KEGG pathway analysis
# Author: Alexandra Millett
# Date created: Nov 27, 2025

# This script is to look at the pathways of interest (Glycolysis and EMP pathways) to see whether any up or downregulated DEGs
# map onto these pathways. 

library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(KEGGREST)
library(pathview)

output_path <- "./plots/kegg_pathway_plots" # relative to top level of repo
working_path <- "../.." # relative to output path

##### IL13 vs IL4 -------------------------
# upregulated gene results
DEG_results_up_IL13_IL4 <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_IL13_vs_IL4.csv")
# downregulated gene results
DEG_results_down_IL13_IL4 <- read.csv("./DEA_outputs/filtered_downregulated_genes_final_IL13_vs_IL4.csv")
# combine up and down regulated genes
DEGs_all <- bind_rows(DEG_results_up_IL13_IL4, DEG_results_down_IL13_IL4)

# gene names
gene_names <- DEGs_all$gene_id
# Log2FC values
FC <- DEGs_all$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# add names to FC data
names(FC) <- gene_entrez_ids$ENTREZID

setwd(output_path)
# glycolysis
pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00010",
         species = "mmu",
         out.suffix = "IL13_IL4_glycolysis")

# Oxidative phosphorylation
pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00190",
         species = "mmu",
         out.suffix = "IL13_IL4_OxPhos")


##### control vs IL13 -------------------------
setwd(working_path)
# upregulated gene results
DEG_results_up_ctrl_IL13 <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_ctrl_vs_IL13.csv")
# downregulated gene results
DEG_results_down_ctrl_IL13 <- read.csv("./DEA_outputs/filtered_downregulated_genes_final_ctrl_vs_IL13.csv")
# combine up and down regulated genes
DEGs_all <- bind_rows(DEG_results_up_ctrl_IL13, DEG_results_down_ctrl_IL13)

# gene names
gene_names <- DEGs_all$gene_id
# Log2FC values
FC <- DEGs_all$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# 0.36% of input gene IDs failed to map

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
pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00190",
         species = "mmu",
         out.suffix = "ctrl_IL13_OxPhos")


##### control vs IL4 -------------------------
setwd(working_path)
# upregulated gene results
DEG_results_up_ctrl_IL4 <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_ctrl_vs_IL4.csv")
# downregulated gene results
DEG_results_up_ctrl_IL4 <- read.csv("./DEA_outputs/filtered_downregulated_genes_final_ctrl_vs_IL4.csv")
# combine up and down regulated genes
DEGs_all <- bind_rows(DEG_results_up_ctrl_IL4, DEG_results_up_ctrl_IL4)

# gene names
gene_names <- DEGs_all$gene_id
# Log2FC values
FC <- DEGs_all$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# 0.36% of input gene IDs failed to map

# add names to FC data
names(FC) <- gene_entrez_ids$ENTREZID

setwd(output_path)
# glycolysis
pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00010",
         species = "mmu",
         out.suffix = "ctrl_IL4_glycolysis")

# Oxidative phosphorylation
pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00190",
         species = "mmu",
         out.suffix = "ctrl_IL4_OxPhos")

setwd(working_path)
