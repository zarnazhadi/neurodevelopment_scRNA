library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(data.table)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample.name <- "20hpf"
sample <- readRDS("hpf20_tSNE.cluster.rds")
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data/markers")
top2 <- read.csv(file = paste0(sample.name, "_top2_markers.csv"), header = TRUE)
top10 <- read.csv(file = paste0(sample.name, "_top10_markers.csv"), header = TRUE)
atlas <- read.csv(file = paste0("atlas_", sample.name, ".csv"), header = TRUE)

# select rows from cells where the enriched gene markers match the gene markers in top2
for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$ENRICHED.MARKERS %like% top2[row,8], ]$ASSIGNED.CELL.TYPE.STATE)
  if(length(identity) == 0){
    identity <- "NULL"
  }
  top2$cell.identity[[row]] <- identity
}

top2 <- apply(top2,2,as.character)
write.csv(as.data.frame(top2), file = paste0(sample.name, "_top2_marker_names.csv")) 

identities.list <- list()
for(row in seq(from= 1, to = dim(top10)[1])){
  identity <- toString(atlas[atlas$ENRICHED.MARKERS %like% top10[row,8], ]$ASSIGNED.CELL.TYPE.STATE)
  if(identity == ""){
    identity <- toString(top10[row, 7])
  }
  top10$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

top10 <- top10 %>%
  filter(!(str_detect(cell.identity, "\\d")))
top10 <- apply(top10,2,as.character)
write.csv(top10, file = paste0(sample.name, "_top10_marker_names.csv"))

current.cluster.ids <- seq(from=0, to=40)
new.cluster.ids <- c("Progenitors", "Eye (Cornea)/Cartilage", "Glutamatergic Neurons*", #2
                     "GABAergic/\nGlutamatergic Neurons", "Neural Crest (Ectomesenchyme),\nMesoderm", "Telencephalon", 
                     "Hindbrain- R5/6", "Prechordal Plate/Polster", "Optic Vesicle (Ventral)",
                     "Radial Glia", "Progenitors (Eye),\nPigment Cell (Iridophore)", "Ventral Midbrain",
                     "Midbrain (Progenitors)", "Ventral Diencephalon", "Rostral Blood Island\n(Myeloid)", 
                     "Eye (Epidermis)", "Epidermis", "Epidermis",
                     "Telencephalon (Pallium, Neuron),\nGlutamatergic Neurons", "RPE", "Periderm", "Heart Primordium",
                     "Retina (Photoreceptor Precursor Cells)", "Macrophages", "Pharyngeal Arch (Posterior)", 
                     "Epidermal Ionocyte", "Cartilage", "Telencephalon (Pallium, Neuron),\nEpidermal Ionocyte", 
                     "Progenitors (Eye)", "Floor Plate", "Lens",
                     "Glia, Periderm", "Muscle", "Neural Crest*", 
                     "Retina (Muller Glia)", "Glutamatergic Neurons", "Oligodendrocytes",
                     "Neuromast (Otic)", "Pigment Cell (Iridophore)", "Dermal Bone", 
                     "Dorsal Habenula/Glutamatergic Neurons")


sample@active.ident <- plyr::mapvalues(x = sample@active.ident, from = current.cluster.ids, to = new.cluster.ids)

pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/NEW_tSNE_plot.pdf"))
pp <- DimPlot(object = sample, reduction = "tsne", label = FALSE, pt.size = 0.5, label.size = 0.5) + theme_classic(base_size = 10) + guides(color = guide_legend(override.aes = list(size=1), ncol=2) )
print(pp)
dev.off()
