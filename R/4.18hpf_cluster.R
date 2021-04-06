library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf18_combined_normalised.rds")

#------------------------------------———FIND VARIABLE GENES—————————————————————---------------------------- 
sample <- FindVariableFeatures(object = sample, mean.function = ExpMean, dispersion.function = LogVMR)
head(x = HVFInfo(object = sample))

#------------------------------------———SCALE DATA—————————————————————---------------------------- 
sample <- ScaleData(object = sample, vars.to.regress = c("nCount_RNA", "percent.mito")) 

#------------------------------------———DIMENSIONAL REDUCTION—————————————————————---------------------------- 
sample <- RunPCA(object = sample,  npcs = 120, verbose = FALSE)

#------------------------------------———PLOTS—————————————————————---------------------------- 
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_hmgn2_dimplot.pdf")
FeaturePlot(object = sample, features = "hmgn2")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_variable_feature_plot.pdf")
VariableFeaturePlot(object = sample)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_heat_map.pdf")
plot <- DimHeatmap(object = sample, reduction = "pca", cells = 200, balanced = TRUE)
dev.off()

#------------------------------------———DETERMINE STATISTICALLY SIGNIFICANT GENES—————————————————————---------------------------- 
sample <- JackStraw(object = sample, reduction = "pca", dims = 120, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)
sample <- ScoreJackStraw(object = sample, dims = 1:120, reduction = "pca")

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_jack_straw_plot.pdf")
JackStrawPlot(object = sample, dims = 1:120, reduction = "pca", xmax = 0.0025)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_elbow_plot.pdf")
ElbowPlot(object = sample, ndims= 120)
dev.off()

#------------------------------------———CLUSTER CELLS—————————————————————---------------------------- 
# find k-nearest neighbours
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:120)

# cluster cells
sample <- FindClusters(sample, resolution = 0.5, algorithm = 1)

# save as .rds file
saveRDS(sample, file = "hpf18.cluster.rds")
