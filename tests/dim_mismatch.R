
a <- 1
b <- 2
c <- 3
d <- 4

X <- array(0,dim=c(a,b,c,d))
y <- rnorm(c)

stan_data <- list(a=a,b=b,c=c,d=d,y=y,X=X)

# rstan::rstan_options("auto_write" = TRUE)
# rstan::stan_model("dim_mismatch.stan")

fit <- rstan::sampling(readr::read_rds("dim_mismatch.rds"), data = stan_data, verbose=TRUE)
sessionInfo()
