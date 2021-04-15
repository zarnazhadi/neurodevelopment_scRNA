library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")

#------------------------------------———MERGING—————————————————————---------------------------- 
sample1 <- readRDS("S1a_hpf20_clean.rds")
sample2 <- readRDS("S1b_hpf20_clean.rds")
sample3 <- readRDS("S2_hpf20_clean.rds")

combined <- merge(sample1, y = c(sample2, sample3), add.cell.ids = c("sample1", "sample2", "sample3"), project = "combined")

# cell names have an added identifier 
head(colnames(combined))

# visualise
table(combined$orig.ident)

#------------------------------------———NORMALISATION—————————————————————---------------------------- 
combined.normalised <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssayData(combined.normalised[1:10, 1:15])

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-normalisation/post_normalisation_VLN_plot.pdf")
plot1 <- VlnPlot(combined.normalised, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
print(plot1)
dev.off()

# output file as .rds
saveRDS(combined.normalised, file = "hpf20_combined_normalised.rds")