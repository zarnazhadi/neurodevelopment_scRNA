library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")

#------------------------------------———MERGING—————————————————————---------------------------- 
sample1 <- readRDS("S3_hpf18_clean.rds")
sample2 <- readRDS("S4_hpf18_clean.rds")

combined <- merge(sample1, y = sample2, add.cell.ids = c("sample1", "sample2"), project = "combined")

# cell names have an added identifier 
head(colnames(combined))

# visualise
table(combined$orig.ident)

#------------------------------------———NORMALISATION—————————————————————---------------------------- 
combined.normalised <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssayData(combined.normalised[1:10, 1:15])

# output file as .rds
saveRDS(combined.normalised, file = "hpf18_combined_normalised.rds")