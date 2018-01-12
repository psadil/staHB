
#' @export
gen_orders <- function(n_conditions){
  library(magrittr)
  out <- gtools::permutations(n_conditions, n_conditions, v=1:n_conditions) %>%
    tibble::as_tibble()
  return(out)
}

#' @export
trim_orders <- function(order_X, base_group=1, end_group=4, degree = "none"){

  if(degree == "none"){
    out <- order_X
  }else if(degree == "min"){
    out <- order_X %>%
      filter(V1 == base_group)
  }else if(degree == "max"){
    out <- order_X %>%
      filter(V1 == base_group) %>%
      filter(V4 == end_group)
  }

  return(out)
}

#' @export
gen_X <- function(d, type){

  X1 <- d %>% modelr::model_matrix(., y1 ~ condition - 1)
  X2 <- d %>% modelr::model_matrix(., y1 ~ condition2 - 1)

  if(type=="mon"){
    X <- abind::abind(as.matrix(X1),as.matrix(X2), along = 3) %>%
      abind::abind(.,., along = 4) %>%
      aperm(.,perm = c(4,3,1,2))
  }else if(type == "nmon"){
    X <- abind::abind(as.matrix(X1),as.matrix(X2), along = 3) %>%
      abind::abind(.,abind::abind(as.matrix(X2),as.matrix(X1), along = 3), along = 4) %>%
      aperm(.,perm = c(4,3,1,2))
  }

  return(X)
}
