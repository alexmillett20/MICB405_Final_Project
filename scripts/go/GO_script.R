# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("topGO")


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(ggplot2))

# When running this script, you will want to generate the GO mapping object first

## Create geneID2GO mapping
allGO2genes <- annFUN.org(
  whichOnto = 'BP',
  feasibleGenes = NULL,
  mapping = 'org.Mm.eg.db',
  ID = 'symbol')

allGO2genes

geneID2GO <- inverseList(allGO2genes)

geneUniverse <- names(geneID2GO)

## Control vs IL4 ##################################
data <- read.csv("./DEA_outputs/filtered_DE_genes_final_ctrl_vs_IL4.csv")

up_gene <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_ctrl_vs_IL4.csv")
down_gene <- read.csv("./DEA_outputs/filtered_downregulated_genes_final_ctrl_vs_IL4.csv")

upregulated_gene_names <- as.character(up_gene$gene_id)
downregulated_gene_names <- as.character(down_gene$gene_id)

up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_gene_names))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_gene_names))
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

up_GO_data <- new("topGOdata", description = "MusMusculus_control_il4", ontology = "BP", allGenes = up_gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO)

down_GO_data <- new("topGOdata", description = "MusMusculus_control_il4", ontology = "BP", allGenes = down_gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO)

up_result <- runTest(up_GO_data, algorithm = "weight01", statistic = "fisher")

down_result <- runTest(down_GO_data, algorithm = "weight01", statistic = "fisher")

up_GO <- GenTable(up_GO_data, weight01 = up_result, orderBy = "up_result", ranksOf = "up_result", topNodes = 30, numChar = 200)

down_GO <- GenTable(down_GO_data, weight01 = down_result, orderBy = "down_result", ranksOf = "down_result", topNodes = 30, numChar = 200)


# ## Control vs IL13 ##############################
data <- read.csv("./DEA_outputs/filtered_DE_genes_final_ctrl_vs_IL13.csv")

up_gene <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_ctrl_vs_IL13.csv")
down_gene <- read.csv("./DEA_outputs/filtered_downregulated_genes_final_ctrl_vs_IL13.csv")

upregulated_gene_names <- as.character(up_gene$gene_id)
downregulated_gene_names <- as.character(down_gene$gene_id)

up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_gene_names))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_gene_names))
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

up_GO_data <- new("topGOdata", description = "MusMusculus_control_il13", ontology = "BP", allGenes = up_gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO)

down_GO_data <- new("topGOdata", description = "MusMusculus_control_il13", ontology = "BP", allGenes = down_gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO)

up_result <- runTest(up_GO_data, algorithm = "weight01", statistic = "fisher")

down_result <- runTest(down_GO_data, algorithm = "weight01", statistic = "fisher")

up_GO <- GenTable(up_GO_data, weight01 = up_result, orderBy = "up_result", ranksOf = "up_result", topNodes = 30, numChar = 200)

down_GO <- GenTable(down_GO_data, weight01 = down_result, orderBy = "down_result", ranksOf = "down_result", topNodes = 30, numChar = 200)

# ## IL13 vs IL4 ################################
data <- read.csv("./DEA_outputs/filtered_DE_genes_final_IL13_vs_IL4.csv")

up_gene <- read.csv("./DEA_outputs/filtered_upregulated_genes_final_IL13_vs_IL4.csv")
down_gene <- read.csv("./DEA_outputs/filtered_downregulated_genes_final_IL13_vs_IL4.csv")

upregulated_gene_names <- as.character(up_gene$gene_id)
downregulated_gene_names <- as.character(down_gene$gene_id)

up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_gene_names))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_gene_names))
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

up_GO_data <- new("topGOdata", description = "MusMusculus_il13_il4", ontology = "BP", allGenes = up_gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO)

down_GO_data <- new("topGOdata", description = "MusMusculus_il13_il4", ontology = "BP", allGenes = down_gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO)

up_result <- runTest(up_GO_data, algorithm = "weight01", statistic = "fisher")

down_result <- runTest(down_GO_data, algorithm = "weight01", statistic = "fisher")

up_GO <- GenTable(up_GO_data, weight01 = up_result, orderBy = "up_result", ranksOf = "up_result", topNodes = 30, numChar = 200)

down_GO <- GenTable(down_GO_data, weight01 = down_result, orderBy = "down_result", ranksOf = "down_result", topNodes = 30, numChar = 200)





# Visualize GO results ##############################

down_GO_filtered <- down_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

down_GO_filtered %>% 
  ggplot(aes(x = Term, y = GeneRatio)) +
  geom_col(width = 0.05) + 
  geom_point(size = 3) +
  coord_flip() # Flip the axes so the x-axis labels are readable

# First, let's arrange the data based on the enrichment ratio. 
down_GO_filtered_arranged <- down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Now let's extract the order of the term column
order_term <- down_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = order_term) + 
  scale_colour_gradient(low = "red", high = "blue")

down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  
  # Add these to make our plot prettier
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # this changes the scale of the axes


# Choose which plots to save
# ggsave("./plots/go-plots/ctrl_vs_il4/GO_down_ctrl_vs_il4.png", width = 8, height = 6)
# ggsave("./plots/go-plots/ctrl_vs_il13/GO_down_ctrl_vs_il13.png", width = 8, height = 6)
ggsave("./plots/go-plots/il13_vs_il4/GO_down_il13_vs_il4.png", width = 8, height = 6)

up_GO_filtered <- up_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

up_GO_filtered %>% 
  ggplot(aes(x = Term, y = GeneRatio)) +
  geom_col(width = 0.05) + 
  geom_point(size = 3) +
  coord_flip() # Flip the axes so the x-axis labels are readable

# First, let's arrange the data based on the enrichment ratio. 
up_GO_filtered_arranged <- up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term, levels = Term))

# Now let's extract the order of the term column
order_term_up <- up_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector
up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = order_term_up) + 
  scale_colour_gradient(low = "red", high = "blue") +
  # Add these to make our plot prettier
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # this changes the scale of the axes  

# Choose which plots to save
# ggsave("./plots/go-plots/ctrl_vs_il4/GO_up_ctrl_vs_il4.png", width = 8, height = 6)
# ggsave("./plots/go-plots/ctrl_vs_il13/GO_up_ctrl_vs_il13.png", width = 8, height = 6)
ggsave("./plots/go-plots/il13_vs_il4/GO_up_il13_vs_il4.png", width = 8, height = 6)

# Prepare data for combined plot
# Add labels to upregulated and downregulated dataframes
up_GO <- up_GO %>% 
  mutate(up_down = "UP")

down_GO <- down_GO %>% 
  mutate(up_down = "DOWN")

# Make a joined dataframe
joined_GO_filtered_arranged <- bind_rows(up_GO, down_GO) %>%
  filter(weight01 <= 0.05) %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term)) %>%
  head(n = 40)

# Extract the column order
order_term_joined <- joined_GO_filtered_arranged %>% 
  pull(Term)

joined_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_point(aes(size= Significant)) +
  coord_flip() +
  scale_x_discrete(limits = order_term_joined) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") +
  theme(panel.border = element_rect(color = "black"), 
        panel.grid = element_line(colour = "grey96"), 
        strip.background = element_rect(colour = "black")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  facet_grid(.~ up_down)

# Choose which plots to save 
# ggsave("./plots/go-plots/ctrl_vs_il4/GO_up_down_ctrl_vs_il4.png", width = 10, height = 6)
# ggsave("./plots/go-plots/ctrl_vs_il13/GO_up_down_ctrl_vs_il13.png", width = 10, height = 6)
ggsave("./plots/go-plots/il13_vs_il4/GO_up_down_il13_vs_il4.png", width = 10, height = 6)
