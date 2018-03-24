
#' @export
'%!in%' <- function(x,y)!('%in%'(x,y))


#' @export
findMatches <- function(d, cutoff=3){

  # distance of this response from each other potential target (all rows, all elements)

  distsOfAllTargets <- purrr::map(d$response_cue,
                                  .f=function(x) purrr::map(d$targets, .f=function(y)
                                    stringdist::stringdist(x, y, method="dl")))

  # index of minimum distance from each potential target
  indsOfCloseEnoughTargets <- purrr::at_depth(distsOfAllTargets,.depth=2, min) %>%
    purrr::map_int(., .f=function(x) ifelse(min(purrr::as_vector(x))<cutoff,
                                     which.min(purrr::as_vector(x)),
                                     NA_integer_)) %>%
    purrr:as_vector(.)
  return(indsOfCloseEnoughTargets)
}


#' @export
countSwaps <- function(d, condition){
  sapply(X=d$potentialMatches, FUN=function(x) d[x, ]$condition==condition) %>%
    unlist(.) %>%
    sum(.)
}

#' @export
listSwapConds <- function(d){

  sapply(X=d$potentialMatches, FUN=function(x) d[x, ]$condition) %>%
    unlist(.)
}

#' @export
nSwaps <- function(d,condition){
  purrr::map_int(d$condOfSwaps, .f=function(y)
    sum(equals(y,condition))
  ) %>%
    purrr::as_vector(.)
}



