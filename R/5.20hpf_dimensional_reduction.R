library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(patchwork)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf20.cluster.rds")

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_jack_straw_plot.pdf")
JackStrawPlot(object = sample, dims = 1:70, reduction = "pca", xmax = 0.0025) + NoLegend()
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_elbow_plot.pdf")
ElbowPlot(object = sample, ndims= 70)
dev.off()

sample <- RunTSNE(object = sample, dims.use = 1:70, do.fast = TRUE)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_tSNE_plot.pdf")
DimPlot(object = sample, reduction = "tsne") + plot_annotation(title = 'T-SNE Clusters for 70 PCs (20 hpf)')
dev.off()

# save as .rds file
saveRDS(sample, file = "hpf20_tSNE.cluster.rds")