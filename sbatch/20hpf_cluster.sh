#!/bin/bash
#SBATCH --ntasks 5
#SBATCH --time 60:00:00
#SBATCH --qos bbdefault
#SBATCH --mail-type ALL

set -e

module purge; module load bluebear
module load bear-apps/2020a
module load R/4.0.0-foss-2020a
export R_LIBS_USER=${HOME}/R/library/${EBVERSIONR}/${BB_CPU}


cd /rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/R # cd to R script

# run Rscripts
#Rscript 4.18hpf_cluster.R
Rscript 4.20hpf_cluster.R
