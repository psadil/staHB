
#' @export
gen_orders <- function(n_conditions){
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
      dplyr::filter(V1 == base_group)
  }else if(degree == "max"){
    out <- order_X %>%
      dplyr::filter(V1 == base_group) %>%
      dplyr::filter(V4 == end_group)
  }

  return(out)
}

#' @export
gen_X_one_order <- function(X1,X2){
  X <- abind::abind(as.matrix(X1),as.matrix(X2), along = 3)
  return(X)
}


#' @export
gen_X <- function(d, type, n_conditions=3, degree="min"){
  # this is the main function for generating order-dependent design matrices

  order_X <- gen_orders(n_conditions) %>%
    trim_orders(., degree = degree, base_group = 1, end_group = n_conditions)

  XX <- array(dim = c(nrow(order_X)^2, n_conditions, 2))
  perms <- gtools::permutations(n=6, r=2, v=1:6, repeats.allowed = TRUE)
  for(i in 1:nrow(XX)){
    XX[i,,1:2] <- abind::abind(order_X[perms[i,1],], order_X[perms[i,2],], along = 3)
  }

  if (type=="full"){
    order_X <- XX
  }else if(type=="mon"){
    order_X <- XX[rowSums( XX[,,1] == XX[,,2]) == n_conditions,,]
  }else if(type == "nmon"){
    order_X <- XX[rowSums( XX[,,1] != XX[,,2]) == n_conditions,,]
  }

  # careful! the default coding scheme for ordered factor is polynomial contrasts!
  # also, default coding scheme for unordered factor (contr.treatment) are not orthogonal to intercept
  for(order in 1:dim(order_X)[1]){
    X1 <- d %>%
      dplyr::mutate(condition_tmp = factor(condition, levels = order_X[order,,1], ordered = TRUE)) %>%
      modelr::model_matrix(~ condition_tmp, contrasts = list(condition_tmp = "contr.treatment"))
    X2 <- d %>%
      dplyr::mutate(condition_tmp = factor(condition, levels = order_X[order,,2], ordered = TRUE)) %>%
      modelr::model_matrix(~ condition_tmp, contrasts = list(condition_tmp = "contr.treatment"))
    if(order==1){
      X <- gen_X_one_order(X1,X2)
    }else{
      X <- X %>%
        abind::abind(.,gen_X_one_order(X1,X2), along=4)
    }
  }
  X <- aperm(X, perm=c(4,3,1,2))

  return(X)
}

