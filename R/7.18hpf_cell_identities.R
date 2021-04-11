library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(data.table)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS(paste0("hpf18_tSNE.cluster.rds"))
sample.name <- "18hpf"
top2 <- read.csv(file = paste0(sample.name, "_top2_markers.csv"), header = TRUE)
top10 <- read.csv(file = paste0(sample.name, "_top10_markers.csv"), header = TRUE)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")
atlas <- as.data.frame(fread("ZFBrainAtlasMaster.tsv", header = TRUE))
atlas <- atlas %>% filter(STAGE == "18 hpf")

for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$`ENRICHED MARKERS` %like% top2[row,8], ]$`ASSIGNED CELL TYPE/STATE`)
  if(length(identity) == 0){
    identity <- "NULL"
  }
  print(row)
  top2$cell.identity[[row]] <- identity
}

identities.list <- list()
for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$`ENRICHED MARKERS` %like% top10[row,8], ]$`ASSIGNED CELL TYPE/STATE`)
  if(identity == ""){
    identity <- toString(top10[row, 7])
    print(identity)
  }
  top10$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

current.cluster.ids <- seq(from=0, to=33)
new.cluster.ids <- c("0", "Optic Vesicle", "Hindbrain", 
                     "Progenitor", "Rhombomere", "Telencephalon",
                     "Progenitors", "diencephalon", "comitted progenitors", 
                     "epidermis", "ventral midbrain", "neural crest", 
                     "diencephalon", "epidermis", "hindbrain",
                     "midbrain", "midbrain", "optic stalk", 
                     "mesoderm", "19", "hindbrain", 
                     "neural crest", "placode", "neurons", 
                     "lens", "RBI", "ventral diencephalon", 
                     "periderm", "heart primordium", "telencephalon",
                     "neurons", "blood precursors", "prechordal plate", 
                     "placode")

sample@active.ident <- plyr::mapvalues(x = sample@active.ident, from = current.cluster.ids, to = new.cluster.ids)

pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/NEW_tSNE_plot.pdf"))
pp <- DimPlot(object = sample, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 2) + theme_classic(base_size = 10) + guides(color = guide_legend(override.aes = list(size=1), ncol=2))
print(pp)
dev.off()
