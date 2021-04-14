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

identities.list <- list()
for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$`ENRICHED MARKERS` %like% top2[row,8], ]$`ASSIGNED CELL TYPE/STATE`)
  if(identity == ""){
    identity <- toString(top2[row, 7])
    print(identity)
  }
  top2$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

identities.list <- list()
for(row in seq(from= 1, to = dim(top10)[1])){
  identity <- toString(atlas[atlas$`ENRICHED MARKERS` %like% top10[row,8], ]$`ASSIGNED CELL TYPE/STATE`)
  if(identity == ""){
    identity <- toString(top10[row, 7])
    print(identity)
  }
  top10$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

current.cluster.ids <- seq(from=0, to=25)
new.cluster.ids <- c("Hindbrain (Floor plate/Boundary- R3)", "Hindbrain (Ant. Dorsal, Rhombic Lip)", "2", 
                     "Telencephalon*", "Optic Vesicle (Dorsal)", "Hindbrain- R7/Pharangeal Arch (posterior)*",
                     "Floor Plate", "Ventral Diencephalon", "Epidermis", 
                     "Placode (epibranchial)*", "Midbrain", "Neural Crest", 
                     "Hindbrain- R6", "Diencephalon", "Placode (Epibranchial)",
                     "Commited Progenitors", "Neural Crest*", "Heart Primordium", 
                     "Placode (Olfactory)", "Telencephalon (Pallium, Neuron)", "Lens (Differentiating)", 
                     "Rostral Blood Island (Myeloid)", "Placode (Octic)", "Periderm", 
                     "Prechordal Plate/Polster", "Blood Precursors")
sample@active.ident <- plyr::mapvalues(x = sample@active.ident, from = current.cluster.ids, to = new.cluster.ids)

pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/NEW_tSNE_plot.pdf"))
pp <- DimPlot(object = sample, reduction = "tsne", label = FALSE, pt.size = 0.5, label.size = 2) + theme_classic(base_size = 10) + 
  guides(color = guide_legend(override.aes = list(size=1), ncol=1))
print(pp)
dev.off()
