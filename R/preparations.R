

#' @export
test_model <- function(cores = 1){

  options(mc.cores = cores)

  data <- gen_dataset(n_item = 20,
                      n_subject = 20,
                      condition_rho = 0,
                      radian_mid = 0,
                      radius_mid = 0)
  data_typed <- apply_type(data, model_type = "full")
  stan_data <- gen_stan_data(data_typed)

  post <- sampling(stanmodels$bivariate_probit_mixed_onec_mo_nonc,
                   data = stan_data,
                   iter = 20,
                   warmup = 10,
                   chains = cores
                   , init_r = 0.25
                   , control = list(adapt_delta = .99)
  )
}


#' @export
setup_job <- function(parallelism = NULL){


  params <- list(data_dir= "data",
                 model= "bivariate_probit_mixed_onec_mo_nonc",
                 flag= "mghpcc-dev",
                 reps= 1,
                 n_condition= 3,
                 condition_rho= seq(from = -0.75, to = 0.75, length.out = 3),
                 n_subject= c(30),
                 n_item= 20,
                 type_model= c("full"),
                 radius= seq(from = 0, to = 1, length.out = 1),
                 radian= seq(from = pi/6, to = pi/3, length.out = 1),
                 chains= 6,
                 iter= 500,
                 warmup= 1000,
                 compile= 0,
                 n_workers= 3
  )



  run_stan <- function(stan_data){
    pars <- c("theta_log", "theta_raw", "lps", "zeta_raw", "theta_raw","condition_mu_raw", "zeta",
              "Mu_unordered", "Mu_ordered",
              "subject_mu_raw", "item_mu_raw", "subject_mu",
              "subject_L", "item_L",
              "condition_omega")

    options(mc.cores = params$chains)
    # rstan::rstan_options("auto_write" = TRUE)



    post <- sampling(stanmodels$bivariate_probit_mixed_onec_mo_nonc,
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
    evaluate_plan(plan = ., rules = list(dataset__ = wf_data$target, MODEL_TYPE = params$type_model))

  wf_stan_data <- drake::drake_plan(stan_data = gen_stan_data(dataset__)) %>%
    plan_analyses(plan = ., datasets = wf_apply_type)

  stan_data_reduce <- drake::gather_plan(wf_stan_data, target = "stan_data", gather = "list")
  typed_reduce <- gather_plan(wf_apply_type, target = "typed_d", gather = "list")


  wf_stan <- drake_plan(post = run_stan(dataset__)) %>%
    plan_analyses(plan = ., datasets = wf_stan_data)

  wf_waic <- drake_plan(waic1 = get_waic(post = POST)) %>%
    evaluate_plan(plan = ., rules = list(POST = wf_stan$target))

  reduce_waic <- gather_plan(wf_waic, target = "waic", gather = "list")


  assemble_d <- drake_plan(d = collect_d(d = dataset__, stan_data=STAN_DATA, waic = WAIC)) %>%
    evaluate_plan(plan = ., rules = list(dataset__ = typed_reduce$target,
                                         STAN_DATA = stan_data_reduce$target,
                                         WAIC = reduce_waic$target))


  wf_plan <- rbind(wf_data, wf_apply_type, wf_stan_data, stan_data_reduce,
                   typed_reduce, wf_stan, wf_waic, reduce_waic, assemble_d)

  con <- drake::drake_config(wf_plan)

  if(is.null(parallelism)){
    drake::make(wf_plan)
  }else if(parallelism == "future_lapply"){
    future::plan(
      future.batchtools::batchtools_lsf,
      template = file.path(devtools::package_file(),"tests", "lsf.tmpl"),
      workers = params$n_workers
    )
    drake::make(wf_plan, parallelism = parallelism)
  }

}





