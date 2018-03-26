#' @export
test_model <- function(cores = getOption("mc.cores", 1L)){

  data <- gen_dataset(n_item = 20,
                      n_subject = 20,
                      condition_rho = 0,
                      radian_mid = 0,
                      radius_mid = 0)
  data_typed <- apply_type(data, model_type = "mon")
  stan_data <- gen_stan_data(data_typed)

  post <- rstan::sampling(stanmodels$bivariate_probit_mixed_onec_mo_nonc,
                          data = stan_data,
                          iter = 20,
                          warmup = 10,
                          chains = cores
                          , cores = cores
                          , init_r = 0.25 # need tight initialization for moving around probits
                          # , control = list(adapt_delta = .9)
  )
}

#' @export
setup_job <- function(jobs = 1, parallelism = drake::default_parallelism(), n_chains = 1,
                      type_model = "full", reps = 1, n_condition_rho = 3,
                      n_subject = 20, n_item = 20, n_radius = 1, n_radian = 1,
                      iter = 500, warmup = 1000){


  wf_data <- drake::drake_plan(data = gen_dataset(n_item = N__ITEM,
                                                  n_subject = N__SUBJECT,
                                                  condition_rho = CONDITION__RHO,
                                                  radian_mid = RADIAN,
                                                  radius_mid = RADIUS)) %>%
    drake::evaluate_plan(., rules = list(N__ITEM = n_item,
                                         N__SUBJECT = n_subject,
                                         CONDITION__RHO = seq(from = -0.75, to = 0.75, length.out = n_condition_rho),
                                         RADIAN = seq(from = pi/6, to = pi/3, length.out = n_radian),
                                         RADIUS = seq(from = 0, to = 1, length.out = n_radius))) %>%
    drake::expand_plan(., values = stringr::str_c("rep", 1:reps))


  wf_apply_type <- drake::drake_plan(data_typed =
                                       apply_type(dataset__, model_type = "MODEL_TYPE"), strings_in_dots = "literals") %>%
    drake::evaluate_plan(plan = ., rules =
                           list(dataset__ = wf_data$target, MODEL_TYPE = type_model))

  wf_stan_data <- drake::drake_plan(stan_data = gen_stan_data(dataset__)) %>%
    drake::plan_analyses(plan = ., datasets = wf_apply_type)

  stan_data_reduce <- drake::gather_plan(wf_stan_data, target = "stan_data", gather = "list")
  typed_reduce <- drake::gather_plan(wf_apply_type, target = "typed_d", gather = "list")


  wf_stan <- drake::drake_plan(post = run_stan(dataset__
                                               , iter = iter
                                               , warmup = warmup
                                               , chains = n_chains)) %>%
    drake::plan_analyses(plan = ., datasets = wf_stan_data)

  wf_waic <- drake::drake_plan(waic1 = get_waic(post = POST)) %>%
    drake::evaluate_plan(plan = ., rules = list(POST = wf_stan$target))

  reduce_waic <- drake::gather_plan(wf_waic, target = "waic", gather = "list")


  assemble_d <- drake::drake_plan(d = collect_d(d = dataset__, stan_data=STAN_DATA, waic = WAIC)) %>%
    drake::evaluate_plan(plan = ., rules = list(dataset__ = typed_reduce$target,
                                                STAN_DATA = stan_data_reduce$target,
                                                WAIC = reduce_waic$target))


  wf_plan <- rbind(wf_data, wf_apply_type, wf_stan_data, stan_data_reduce,
                   typed_reduce, wf_stan, wf_waic, reduce_waic, assemble_d)

  con <- drake::drake_config(wf_plan)

  drake::make(wf_plan, jobs = jobs, parallelism = parallelism)

  return(con)
}

#' @export
run_stan <- function(stan_data, iter=1000, warmup=500, chains=1){

  pars <- c("theta_log", "theta_raw", "lps", "zeta_raw", "theta_raw","condition_mu_raw", "zeta",
            "Mu_unordered", "Mu_ordered",
            "subject_mu_raw", "item_mu_raw", "subject_mu",
            "subject_L", "item_L",
            "condition_omega")

  post <- rstan::sampling(stanmodels$bivariate_probit_mixed_onec_mo_nonc,
                          data = stan_data,
                          iter = warmup + iter,
                          warmup = warmup,
                          chains = chains
                          , cores = chains
                          , pars = pars
                          , include = FALSE
                          , init_r = 0.25 # need tight initialization for moving around probits
                          # , control = list(adapt_delta = .9)
  )

  return(post)
}



