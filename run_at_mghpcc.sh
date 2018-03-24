#!/bin/bash

module load singularity/singularity-2.4.5

#BSUB -n 36
#BSUB -R rusage[mem=1024] # ask for 1GB per job slot, or 36GB total
#BSUB -W 2:00
#BSUB -q long # which queue we want to run in
#BSUB -R "span[hosts=1]" # All job slots on the same node (needed for threaded applications)
#BSUB -J staHB

singularity exec psadil-staHB-master-latest.simg Rscript -e "library(staHB); setup_job(jobs=6, n_chains=6, type_model = c('nmon','mon','full'), n_radius=3, n_radian=3)"
