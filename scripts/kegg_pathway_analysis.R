library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(KEGGREST)
library(pathview)

# upregulated gene results
DEG_results_up_IL13_IL4 <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_IL13_vs_IL4.csv")
# gene names
gene_names <- DEG_results_up_IL13_IL4$gene_id
# Log2FC values
FC <- DEG_results_up_IL13_IL4$log2FoldChange
# convert to entrez ID from gene symbol
gene_entrez_ids <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# add names to FC data
names(FC) <- gene_entrez_ids$ENTREZID

pathview(gene.data = FC,
         gene.idtype = "entrez",
         pathway.id = "mmu00010",
         species = "mmu")
keggList("pathway", "mmu")[1:10]
