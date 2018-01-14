
gen_grand_data <- function(params){

  d <- crossing(n_item = params$n_item
                , n_subject = params$n_subject
                , condition_rho = params$condition_rho
                , type = factor(params$type, levels = c("mon","nmon"))
                , subject_scale = params$subject_scale
                , item_scale = params$item_scale
                , subject_rho = params$subject_rho
                , item_rho = params$item_rho) %>%
    mutate(expt = factor(1:n())) %>%
    group_by(expt, condition_rho, type, subject_scale, subject_rho, item_scale, item_rho) %>% #expt specific params go here
    nest() %>%
    mutate(data = purrr::map(data, ~crossing(item = 1:.x$n_item
                                             , subject = 1:.x$n_subject
                                             , condition = 1:4))) %>%
    unnest(data) %>%
    mutate(condition_mu = map(condition, ~cbind(params$condition1_mu[.x], params$condition2_mu[.x]))) %>%
    mutate(subject_mu = map2(subject_scale,subject_rho, ~ mvtnorm::rmvnorm(1, mean = rep(0,2), sigma = gen_cor_mat(.x,.y)))) %>%
    mutate(item_mu = map2(item_scale, item_rho, ~ mvtnorm::rmvnorm(1, mean = rep(0,2), sigma = gen_cor_mat(.x,.y)))) %>%
    mutate(Mu = pmap(list(condition_mu, subject_mu, item_mu, condition_rho),
                     function(a,b,c,d) mvtnorm::rmvnorm(1, mean = a + b + c, sigma = matrix(c(1, d, d, 1), 2, 2)) )) %>%
    mutate(evidence_x = map_dbl(Mu, ~purrr::pluck(.x, 1)),
           evidence_y = map_dbl(Mu, ~purrr::pluck(.x, 2))
    ) %>%
    mutate(y_sim = map(Mu, ~if_else(.x > 0, 1, 0))) %>%
    select(-Mu) %>%
    mutate(y1 = map_dbl(y_sim, ~purrr::pluck(.x, 1)),
           y2 = map_dbl(y_sim, ~purrr::pluck(.x, 2))
    )

  return(d)
}


gen_stan_data <- function(d){

  stan_data <- d %>%
    select(condition, item, y1, y2, subject) %>%
    mutate(condition = factor(condition),
           item = factor(item),
           subject = factor(subject)) %>%
    tidybayes::compose_data() %>%
    c(.,
      D = 2,
      priors = list(c(1, 2, 1, 2, 1, 1.5, 1)),
      y  = list(cbind(d$y1, d$y2))
    )

  stan_data$X <- gen_X(d,unique(d$type))
  stan_data$n_orders <- dim(stan_data$X)[1]


  return(stan_data)
}

gen_cor_mat <- function(std, rho, D=2){

  stds <- rep(std, times=D)
  Omega <- matrix(c(1, rho, rho, 1), 2, 2) * (stds %*% t(stds))

  return(Omega)

}
