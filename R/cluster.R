library("Matrix")
library(tidyverse)
library("Seurat")
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

hpf18 <- readRDS("hpf18_combined_normalised.rds")
hpf20 <- readRDS("hpf20_combined_normalised.rds")
ex <- readRDS("GSE158142_zf10s_cc_filt.cluster.rds") # example

# find variable genes
hpf18 <- FindVariableFeatures(object = hpf18, mean.function = ExpMean, dispersion.function = LogVMR)
hpf20 <- FindVariableFeatures(object = hpf20, mean.function = ExpMean, dispersion.function = LogVMR)

head(x = HVFInfo(object = hpf18))
head(x = HVFInfo(object = hpf20))

# scaling data - error
hpf18 <- ScaleData(object = hpf18) # could not regress vars


