

#' @export
apply_type <- function(d, model_type){
  d %<>% dplyr::mutate(type = model_type)

  return(d)
}

#' @export
collect_d <- function(d, stan_data, waic){
  d2 <- tibble::tibble(data = d,
               stan_data = stan_data,
               waic = waic)
  return(d2)
}


#' @export
get_waic <- function(post){
  log_lik <- loo::extract_log_lik(post)
  w <- loo::waic(log_lik)

  return(w)
}
