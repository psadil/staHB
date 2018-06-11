#' @export
setup_cart <- function(jobs = 1, parallelism = drake::default_parallelism(), n_chains = 1,
                      type_model = "full", reps = 1, n_condition_rho = 3,
                      n_subject = 1, n_item = 1, condition_x = c(-1, .7, -1/4, 1), condition_y = c(-1, -1/3, 1, 2/3 ),
                      subject_scale = 0, item_scale = 0, tau = 1,
                      iter = 500, warmup = 1000, cache_dir = NULL){


  wf_data <- drake::drake_plan(data = gen_dataset_cart(n_item = N__ITEM,
                                                  n_subject = N__SUBJECT,
                                                  condition_rho = CONDITION__RHO,
                                                  condition_x = condition_x,
                                                  condition_y = condition_y,
                                                  subject_scale = subject_scale,
                                                  item_scale = item_scale,
                                                  tau = tau)) %>%
    drake::evaluate_plan(., rules = list(N__ITEM = n_item,
                                         N__SUBJECT = n_subject,
                                         CONDITION__RHO = seq(from = -0.75, to = 0.75, length.out = n_condition_rho)
                                         )) %>%
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

  con <- drake::drake_config(wf_plan, cache = cache_dir)

  drake::make(wf_plan, jobs = jobs, parallelism = parallelism, cache = cache_dir)

  return(NULL)
}


to_cmr <- function(d){

  out <- d %>%
    dplyr::select(subject, condition, evidence_x, evidence_y, item) %>%
    dplyr::mutate(group = 1) %>%
    tidyr::gather(question, result, c("evidence_x","evidence_y")) %>%
    dplyr::mutate(question = plyr::mapvalues(question, from = unique(question), to = 1:2),
                  question = as.numeric(question)) %>%
    dplyr::group_by(subject, group, question, condition) %>%
    summarise(result = mean(result)) %>%
    spread(condition, result) %>%
    ungroup() %>%
    mutate(subject = as.numeric(as.character(subject)))

  return(out)
}
