library(dplyr)
library(data.table)

# subset atlas to 20 hours stage with cells only involved in neurodevelopment
setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

atlas <- as.data.frame(fread("ZFBrainAtlasMaster.tsv", header = TRUE))
atlas_20hpf <- atlas %>% filter(STAGE == "20 hpf")
cell_types <- atlas_20hpf$`ASSIGNED CELL TYPE/STATE`
as.data.frame(cell_types) # check that these are separated properly
 
groups <- c("brain",
            "brain",
            "unsure",
            "unsure",
            "unsure",
            "eye",
            "brain",
            "brain",
            "brain",
            "other",
            "other",
            "unsure",
            "eye",
            "unsure",
            "brain",
            "brain",
            "brain",
            "brain",
            "brain",
            "other",
            "brain",
            "brain",
            "other",
            "brain",
            "other",
            "brain",
            "eye",
            "other",
            "other",
            "brain",
            "brain",
            "other",
            "eye",
            "brain",
            "brain",
            "brain",
            "other",
            "other",
            "brain",
            "brain",
            "eye",
            "brain",
            "other",
            "brain",
            "eye",
            "brain",
            "eye",
            "unsure",
            "eye",
            "eye",
            "eye",
            "other",
            "brain",
            "brain",
            "other",
            "other",
            "other",
            "brain",
            "other",
            "other",
            "other",
            "brain",
            "unsure",
            "brain",
            "eye",
            "other",
            "other",
            "other",
            "brain",
            "eye")

cell_types <- cbind(cell_types, groups) # for reference
atlas_20hpf <- cbind(atlas_20hpf, groups)

# 20hpf and brain cells/neurons
atlas_20hpf_brain <- atlas_20hpf %>% filter(groups == "brain")

# write csv 
setwd("/rds/projects/v/vianaj-development-rna/appDir/R")
write.csv(atlas_20hpf, file = "atlas_20hpf.csv")


