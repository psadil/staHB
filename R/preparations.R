#' @export
test_model <- function(cores = getOption("mc.cores", 1L)){

  data <- gen_dataset(n_item = 20,
                      n_subject = 20,
                      condition_rho = 0,
                      radian_mid = 0,
                      radius_mid = 0)
  data_typed <- apply_type(data, model_type = "full")
  stan_data <- gen_stan_data(data_typed)

  post <- rstan::sampling(stanmodels$bivariate_probit_mixed_onec_mo_nonc,
                          data = stan_data,
                          iter = 20,
                          warmup = 10,
                          chains = cores
                          , cores = cores
                          , init_r = 0.25
                          , control = list(adapt_delta = .99)
  )
}


#' @export
setup_job <- function(jobs = 1, parallelism = drake::default_parallelism(), n_chains = 1){

  params <- list(data_dir = "data",
                 model = "bivariate_probit_mixed_onec_mo_nonc",
                 flag = "mghpcc-dev",
                 reps = 1,
                 n_condition = 3 ,
                 condition_rho = seq(from = -0.75, to = 0.75, length.out = 3),
                 n_subject = c(30),
                 n_item = 20,
                 type_model = c("full"),
                 radius = seq(from = 0, to = 1, length.out = 1),
                 radian = seq(from = pi/6, to = pi/3, length.out = 1),
                 chains = n_chains,
                 iter= 500,
                 warmup= 1000
  )


  run_stan <- function(stan_data){
    pars <- c("theta_log", "theta_raw", "lps", "zeta_raw", "theta_raw","condition_mu_raw", "zeta",
              "Mu_unordered", "Mu_ordered",
              "subject_mu_raw", "item_mu_raw", "subject_mu",
              "subject_L", "item_L",
              "condition_omega")


    post <- rstan::sampling(stanmodels$bivariate_probit_mixed_onec_mo_nonc,
                            data = stan_data,
                            iter = params$warmup + params$iter,
                            warmup = params$warmup,
                            chains = params$chains
                            , cores = params$chains
                            , pars = pars
                            , include = FALSE
                            , init_r = 0.25
                            , control = list(adapt_delta = .99)
    )


    return(post)
  }

  wf_data <- drake::drake_plan(data = gen_dataset(n_item = N__ITEM,
                                                  n_subject = N__SUBJECT,
                                                  condition_rho = CONDITION__RHO,
                                                  radian_mid = RADIAN,
                                                  radius_mid = RADIUS)) %>%
    drake::evaluate_plan(., rules = list(N__ITEM = params$n_item,
                                         N__SUBJECT = params$n_subject,
                                         CONDITION__RHO = params$condition_rho,
                                         RADIAN = params$radian,
                                         RADIUS = params$radius)) %>%
    drake::expand_plan(., values = stringr::str_c("rep", 1:params$reps))



  wf_apply_type <- drake::drake_plan(data_typed = apply_type(dataset__, model_type = "MODEL_TYPE"), strings_in_dots = "literals") %>%
    drake::evaluate_plan(plan = ., rules = list(dataset__ = wf_data$target, MODEL_TYPE = params$type_model))

  wf_stan_data <- drake::drake_plan(stan_data = gen_stan_data(dataset__)) %>%
    drake::plan_analyses(plan = ., datasets = wf_apply_type)

  stan_data_reduce <- drake::gather_plan(wf_stan_data, target = "stan_data", gather = "list")
  typed_reduce <- drake::gather_plan(wf_apply_type, target = "typed_d", gather = "list")


  wf_stan <- drake::drake_plan(post = run_stan(dataset__)) %>%
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

}





