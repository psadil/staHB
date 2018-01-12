
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

  order_X <- gen_orders(4) %>%
    trim_orders(., degree = "max")

  X1 <- d %>% modelr::model_matrix(., ~ condition - 1)
  X2 <- d %>% modelr::model_matrix(., ~ condition2 - 1)

  if(type=="mon"){
    X <- gen_X_same_order(X1) %>% # dim 3 := question
      abind::abind(., gen_X_same_order(X2), along = 4) %>% # dim 4 := order
      aperm(.,perm = c(4,3,1,2))
  }else if(type == "nmon"){
    X <- gen_X_same_order(X1) %>% # dim 3 := question
      abind::abind(., gen_X_same_order(X2), along = 4) %>% # dim 4 := order
      abind::abind(., gen_X_mixed_order(X1,X2), along = 4) %>% # dim $ := order
      abind::abind(., gen_X_mixed_order(X2,X1), along = 4) %>% # dim 4 := order
      aperm(.,perm = c(4,3,1,2)) # final is (order, question, obsvervation, condition)
  }

  return(X)
}

gen_X_same_order <- function(X){
  out <- abind::abind(as.matrix(X),as.matrix(X), along = 3)
  return(out)
}

gen_X_mixed_order <- function(X1, X2){
  out <- abind::abind(as.matrix(X1),as.matrix(X2), along = 3)
  return(out)
}
