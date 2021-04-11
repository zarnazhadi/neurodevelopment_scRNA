library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(data.table)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample.name <- "20hpf"
sample <- readRDS(paste0("hpf20_tSNE.cluster.rds"))
top2 <- read.csv(file = paste0(sample.name, "_top2_markers.csv"), header = TRUE)
top10 <- read.csv(file = paste0(sample.name, "_top10_markers.csv"), header = TRUE)
atlas <- read.csv(file = paste0("atlas_", sample.name, ".csv"), header = TRUE)

# select rows from cells where the enriched gene markers match the gene markers in top2
for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$ENRICHED.MARKERS %like% top2[row,8], ]$ASSIGNED.CELL.TYPE.STATE)
  if(length(identity) == 0){
    identity <- "NULL"
  }
  print(row)
  top2$cell.identity[[row]] <- identity
}

identities.list <- list()
for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$ENRICHED.MARKERS %like% top10[row,8], ]$ASSIGNED.CELL.TYPE.STATE)
  if(identity == ""){
    identity <- toString(top10[row, 7])
    print(identity)
  }
  top10$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

current.cluster.ids <- seq(from=0, to=48)
new.cluster.ids <- c("Progenitors", "Hindbrain", "Pectoral Fin",
                     "Gabaergic", "Glutamatergic", "Neural Crest", 
                     "Tencephalon", "Optic Vesicle", "Placode", 
                     "Progenitors", "10", "Eye",
                     "Radial Glia", "Ventral Midbrain", "Midbrain Progenitors", 
                     "Diencephalon", "RBI", "17",
                     "Epidermis", "Retina", "Epidermis",
                     "21", "Glutamatergic", "Pharangeal Arch",
                     "24", "25", "Heart Primadorium", "Cerebellum",
                     "28", "Telencephalon", "Cartilage", 
                     "Retina", "32", "Ventral Diencephalon",
                     "Lens", "Glia", "Neurons", "Muscle", 
                     "Cornea", "Retina", "Progenitors", 
                     "41", "Neuromast", "Purkinje",
                     "Iridophore", "Dermal Bone", "Dorsal Habenula",
                     "47", "Polster")


sample@active.ident <- plyr::mapvalues(x = sample@active.ident, from = current.cluster.ids, to = new.cluster.ids)

pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/NEW_tSNE_plot.pdf"))
pp <- DimPlot(object = sample, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 3) + theme_classic(base_size = 10) + guides(color = guide_legend(override.aes = list(size=1), ncol=2) )
print(pp)
dev.off()
