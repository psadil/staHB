
devtools::load_all()
cache1 <- drake::new_cache(".cmr")

setup_cart(jobs=1, n_chains=1, type_model = c('full'), n_condition_rho = 1, n_item = 1, n_subject = 1, cache_dir = cache1)
