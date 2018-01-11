
gen_grand_data <- function(params){

  condition_stds <- c(1, 1)
  condition_Sigma <- matrix(c(1, params$condition_rho, params$condition_rho, 1), 2, 2) * (condition_stds %*% t(condition_stds))

  item_stds <- c(params$item_scale, params$item_scale)
  item_Sigma <- matrix(c(1, params$item_rho, params$item_rho, 1), 2, 2) * (item_stds %*% t(item_stds))

  subject_stds <- c(params$subject_scale, params$subject_scale)
  subject_Sigma <- matrix(c(1, params$subject_rho, params$subject_rho, 1), 2, 2) * (subject_stds %*% t(subject_stds))

  d <- expand.grid(n_item = seq(10, 50, length.out = 5)
                   , n_subject = seq(10, 50, length.out = 5)) %>%
    mutate(expt = 1:n()) %>%
    group_by(expt) %>%
    nest() %>%
    mutate(data2 = purrr::map(data, .x %>%
                               crossing(item = 1:n_item
                                        , subject = 1:n_subject
                                        , condition = 1:4))) %>%
    unnest(data2) %>%
    mutate(condition_mu = map(condition, ~ cbind(params$condition1_mu[.x], params$condition2_mu[.x]) + mvtnorm::rmvnorm(1, mean = rep(0,2), sigma = condition_Sigma))) %>%
    mutate(subject_mu = map(item, ~ mvtnorm::rmvnorm(1, mean = rep(0,2), sigma = subject_Sigma))) %>%
    mutate(item_mu = map(subject, ~ mvtnorm::rmvnorm(1, mean = rep(0,2), sigma = item_Sigma))) %>%
    mutate(Mu = pmap(list(condition_mu, subject_mu, item_mu), function(a,b,c) a + b + c)) %>%
    mutate(evidence_x = map_dbl(Mu, ~purrr::pluck(.x, 1)),
           evidence_y = map_dbl(Mu, ~purrr::pluck(.x, 2))
    ) %>%
    mutate(y_sim = map(Mu, ~if_else(.x > 0, 1, 0))) %>%
    select(-Mu) %>%
    mutate(y1 = map_dbl(y_sim, ~purrr::pluck(.x, 1)),
           y2 = map_dbl(y_sim, ~purrr::pluck(.x, 2))
    )%>%
    mutate(condition2 = plyr::mapvalues(condition, from = 1:4, to = c(1,3,2,4)))

  return(d)
}

