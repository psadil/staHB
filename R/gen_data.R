
#' @export
gen_dataset <- function(n_item, n_subject, n_condition = 3,
                        condition_rho, subject_scale = sqrt(.25), item_scale = sqrt(.25), subject_rho = 0, item_rho = 0,
                        radian_mid, radius_mid, theta2 = pi/4, base_length = 2,
                        tau = 1){

  # tau: scale of error. typically := 1, unless set to 0 to test exact condition effect

  # range between condition effects. used in determining std of item effects
  d <- tidyr::crossing(item = 1:n_item
                , subject = 1:n_subject
                , condition = 1:n_condition) %>%
    dplyr::mutate(item = factor(item)
           , subject = factor(subject)
           , condition = factor(condition)) %>%
    dplyr::mutate(radius = dplyr::case_when(condition == 1 ~ -base_length/2,
                              condition == 2 ~ radius_mid,
                              condition == 3 ~ base_length/2),
           radian = dplyr::case_when(condition == 1 ~ theta2,
                              condition == 2 ~ theta2 + radian_mid,
                              condition == 3 ~ theta2)) %>%
    dplyr::mutate(x = radius * cos(radian),
           y = radius * sin(radian)) %>%
    dplyr::mutate(condition_mu = purrr::map2(x, y, ~cbind(.x, .y))) %>%
    dplyr::mutate(subject_mu = purrr::map2(subject_scale, subject_rho, ~ r_mvn(n=1, tau = .x, rho = .y, mu=c(0,0))),
           item_mu = purrr::map2(item_scale, item_rho, ~ r_mvn(n=1, tau = .x, rho = .y, mu=c(0,0)))) %>%
    dplyr::mutate(Mu = purrr::pmap(list(condition_mu, subject_mu, item_mu, condition_rho, tau),
                     function(a,b,c,d, e) r_mvn(n=1, mu = a + b + c, tau = e, rho = d) )) %>%
    dplyr::mutate(evidence_x = purrr::map_dbl(Mu, ~purrr::pluck(.x, 1)),
           evidence_y = purrr::map_dbl(Mu, ~purrr::pluck(.x, 2))
    ) %>%
    dplyr::mutate(y_sim = purrr::map(Mu, ~dplyr::if_else(.x > 0, 1, 0))) %>%
    dplyr::select(-Mu) %>%
    dplyr::mutate(y1 = purrr::map_dbl(y_sim, ~purrr::pluck(.x, 1)),
           y2 = purrr::map_dbl(y_sim, ~purrr::pluck(.x, 2))
    )

  return(d)
}

#' @export
r_mvn <- function(n=1, mu = matrix(c(0,0),nrow=1), tau = 1, rho = 0){

  center <- matrix(mu,nrow=2)
  Omega <- matrix(c(1,rho,rho,1),nrow = 2)
  L_Omega <- t(chol(Omega))
  tau_diag <- diag(tau,2)
  error_raw <- rbind(rnorm(n),rnorm(n))
  shift <- do.call("rbind", rep(list(center), n))

  effect <- t(shift + ( (tau_diag %*% L_Omega) %*% error_raw ))

  return(effect)
}

#' @export
gen_stan_data <- function(d, degree="min", rho = NULL){

  stan_data <- d %>%
    dplyr::select(condition, item, y1, y2, subject) %>%
    dplyr::mutate(condition = factor(condition),
           item = factor(item),
           subject = factor(subject)) %>%
    tidybayes::compose_data() %>%
    c(.,
      D = 2,
      priors = list(c(1, 1, 1, 2, 1, 1, 1)),
      y  = list(cbind(d$y1, d$y2))
    )
  n_conditions <- dplyr::n_distinct(d$condition)

  stan_data$X <- gen_X(d, type = unique(d$type), n_conditions = n_conditions, degree = degree)
  stan_data$n_orders <- dim(stan_data$X)[1]

  if(!(is.null(rho))){
    stan_data$condition_omega <- matrix(c(1, rho, rho,1), nrow=2)
  }

  return(stan_data)
}

