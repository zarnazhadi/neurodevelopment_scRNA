library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

# open clean data 
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")

sample.list <- list.files(pattern=glob2rx("*hpf18_clean.rds"))

for(sample in 1:length(sample.list)){
  sample.name <- substr(sample.list[sample], 1,3)
  if(grepl(sample.name, pattern="_")){
    sample.name <- substr(sample.list[sample], 1,2)
  }
  
  sample.obj <- readRDS(sample.list[sample])
# RNA scatter plot 
  pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/CLEAN_", sample.name,"_RNA_scatter_plot.pdf"))
  p <- FeatureScatter(object = sample.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p
  dev.off()
}