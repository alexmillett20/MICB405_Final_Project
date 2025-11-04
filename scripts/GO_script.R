# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("topGO")


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(ggplot2))

data <- read.csv("./DEA_outputs/ctrl_vs_il13_results.csv")

allGO2genes <- annFUN.org(
  whichOnto = 'BP',
  feasibleGenes = NULL,
  mapping = 'org.Mm.eg.db',
  ID = 'symbol')

allGO2genes

geneID2GO <- inverseList(allGO2genes)

geneUniverse <- names(geneID2GO)

up_gene <- read.csv("./DEA_outputs/ctrl_vs_il13_up.csv")
down_gene <- read.csv("./DEA_outputs/ctrl_vs_il13_down.csv")

upregulated_gene_names <- as.character(up_gene$gene_id)
downregulated_gene_names <- as.character(down_gene$gene_id)

up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_gene_names))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_gene_names))
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

up_GO_data <- new("topGOdata",
                  description = "MusMusculus_control_il13",
                  ontology = "BP",
                  allGenes = up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

down_GO_data <- new("topGOdata",
                    description = "MusMusculus_control_il13",
                    ontology = "BP",
                    allGenes = down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")

down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

up_summary <- GenTable(up_GO_data,
                       weight01 = up_result,
                       orderBy = "up_result",
                       ranksOf = "up_result",
                       topNodes = 30)

down_summary <- GenTable(down_GO_data,
                         weight01 = down_result,
                         orderBy = "down_result",
                         ranksOf = "down_result",
                         topNodes = 30)



up_plot <- ggplot(up_summary, aes(x = reorder(Term, -log10(as.numeric(weight01))),
                       y = -log10(as.numeric(weight01)))) +
  geom_col(fill = "#E64B35FF") +
  coord_flip() +
  labs(
    title = "Top 30 Enriched GO Terms\n (Upregulated Genes)",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("./plots/GO_up_ctrl_vs_il13.png", up_plot)

down_plot <- ggplot(down_summary, aes(x = reorder(Term, -log10(as.numeric(weight01))),
                         y = -log10(as.numeric(weight01)))) +
  geom_col(fill = "#4DBBD5FF") +
  coord_flip() +
  labs(
    title = "Top 30 Enriched GO Terms\n (Downregulated Genes)",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("./plots/GO_down_ctrl_vs_il13.png", down_plot)

