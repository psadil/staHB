

gen_orders <- function(n_conditions){
  library(magrittr)
  out <- gtools::permutations(n_conditions, n_conditions, v=1:n_conditions) %>%
    tibble::as_tibble()
  return(out)
}


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
