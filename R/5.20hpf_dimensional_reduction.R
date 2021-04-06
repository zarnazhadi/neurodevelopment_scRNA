library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf20.cluster.rds")

sample <- RunTSNE(object = sample, dims.use = 1:70, do.fast = TRUE)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_tSNE_plot.pdf")
DimPlot(object = sample, reduction = "tsne")
dev.off()

# save as .rds file
saveRDS(sample, file = "hpf20_tSNE.cluster.rds")