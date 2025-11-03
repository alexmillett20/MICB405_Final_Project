if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))

