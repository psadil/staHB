#!/bin/sh

Rscript -e "rstan::stan_model(file.path(devtools::package_file(), 'inst', 'stan', 'bivariate_probit_mixed_onec_mo_nonc.stan'))"

Rscript -e "rmarkdown::render('mghpcc.Rmd', params = list(n_workers = 3))"
