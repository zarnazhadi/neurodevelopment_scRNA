library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

# open clean data 
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")

s.list <- list.files(pattern=glob2rx("*hpf18_clean.rds"))

# S3
s <- 1
s.name <- substr(s.list[s], 1,3)
  if(grepl(s.name, pattern="_")){
  s.name <- substr(s.list[s], 1,2)
  }
  
s.obj <- readRDS(s.list[s])
  
# RNA scatter plot 
pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/CLEAN_", s.name, "_RNA_scatter_plot.pdf"))
FeatureScatter(object = s.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
  
# VLN plot
pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/CLEAN_", s.name, "_mito_vln_plot.pdf"))
VlnPlot(object = s.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# S4
s <- 2
s.name <- substr(s.list[s], 1,3)
if(grepl(s.name, pattern="_")){
  s.name <- substr(s.list[s], 1,2)
}

s.obj <- readRDS(s.list[s])

# RNA scatter plot 
pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/CLEAN_", s.name, "_RNA_scatter_plot.pdf"))
FeatureScatter(object = s.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# VLN plot
pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/CLEAN_", s.name, "_mito_vln_plot.pdf"))
VlnPlot(object = s.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()