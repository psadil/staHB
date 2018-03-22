#!/bin/sh

Rscript -e "library(staHB); prepare_long_run('bivariate_probit_mixed_onec_mo_nonc'); setup_job('future_lapply')"
