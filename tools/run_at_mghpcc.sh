#!/bin/bash

module load singularity/singularity-current

#BSUB -n 48
#BSUB -R rusage[mem=1024] # ask for 1GB per job slot, or 36GB total
#BSUB -W 48:00
#BSUB -q long # which queue we want to run in
#BSUB -R "span[hosts=1]" # All job slots on the same node (needed for threaded applications)
#BSUB -J staHB

singularity exec psadil-staHB-master-latest.simg Rscript -e "library(staHB); cache1 <- drake::new_cache('.cmr'); setup_cart(jobs=6, n_chains=8, type_model =  c('nmon','mon','full'), n_condition_rho = 1, n_item = 20, n_subject = c(20,50), cache_dir = cache1)"
