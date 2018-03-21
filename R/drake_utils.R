run_stan <- function(stan_data){
  pars <- c("theta_log", "theta_raw", "lps", "zeta_raw", "theta_raw","condition_mu_raw", "zeta",
            "Mu_unordered", "Mu_ordered",
            "subject_mu_raw", "item_mu_raw", "subject_mu",
            "subject_L", "item_L",
            "condition_omega")

  model <- stan_model(file.path(devtools::package_file(), "stan", paste0(params$model, ".stan")))

  post <- sampling(model,
                   data = stan_data,
                   iter = params$warmup + params$iter,
                   warmup = params$warmup,
                   chains = params$chains
                   , pars = pars
                   , include = FALSE
                   , init_r = 0.25
                   , control = list(adapt_delta = .99)
  )

  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }


  return(post)
}

apply_type <- function(d, model_type){
  d %<>% mutate(type = model_type)

  return(d)
}

collect_d <- function(d, stan_data, waic){
  d2 <- tibble(data = d,
               stan_data = stan_data,
               waic = waic)
  return(d2)
}



get_waic <- function(post){
  log_lik <- loo::extract_log_lik(post)
  w <- loo::waic(log_lik)

  return(w)
}
