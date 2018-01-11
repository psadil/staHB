
library(rstan)

context("Probability function testing")

rstan::expose_stan_functions(file.path(devtools::package_file(), "stan", paste0("bivariate_probit_mixed", ".stan")))

z_0 <- matrix(c(0,0), ncol=2)
z_n20 <- matrix(c(-20,-20), ncol=2)
z_p20 <- matrix(c(20,20), ncol=2)
z_p1 <- matrix(c(1,1), ncol=2)
z_pn1 <- matrix(c(1,-1), ncol=2)

rho_0 <- 0
rho_1 <- 1
rho_5 <- 0.5
rho_n1 <- -1

testthat::test_that("binormal_cdf is close with z_0", {
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_0, rho_0)),
                         mvtnorm::pmvnorm(upper = as.vector(z_0), corr = matrix(c(1,rho_0,rho_0,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_0, rho_1)),
                         mvtnorm::pmvnorm(upper = as.vector(z_0), corr = matrix(c(1,rho_1,rho_1,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_0, rho_5)),
                         mvtnorm::pmvnorm(upper = as.vector(z_0), corr = matrix(c(1,rho_5,rho_5,1), ncol=2))[[1]] )
})

testthat::test_that("binormal_cdf is close with z_n20, z_p20", {
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_n20, rho_5)),
                         mvtnorm::pmvnorm(upper = as.vector(z_n20), corr = matrix(c(1,rho_5,rho_5,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_p20, rho_5)),
                         mvtnorm::pmvnorm(upper = as.vector(z_p20), corr = matrix(c(1,rho_5,rho_5,1), ncol=2))[[1]] )
})

testthat::test_that("binormal_cdf is close with z_p1, z_pn1", {
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_p1, rho_5)),
                         mvtnorm::pmvnorm(upper = as.vector(z_p1), corr = matrix(c(1,rho_5,rho_5,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_pn1, rho_5)),
                         mvtnorm::pmvnorm(upper = as.vector(z_pn1), corr = matrix(c(1,rho_5,rho_5,1), ncol=2))[[1]] )
})

testthat::test_that("binormal_cdf is close with z_p1, z_pn1", {
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_p1, rho_n1)),
                         mvtnorm::pmvnorm(upper = as.vector(z_p1), corr = matrix(c(1,rho_n1,rho_n1,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_pn1, rho_n1)),
                         mvtnorm::pmvnorm(upper = as.vector(z_pn1), corr = matrix(c(1,rho_n1,rho_n1,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_p1, rho_1)),
                         mvtnorm::pmvnorm(upper = as.vector(z_p1), corr = matrix(c(1,rho_1,rho_1,1), ncol=2))[[1]] )
  testthat::expect_equal(as.vector(binormal_cdf_vector(z_pn1, rho_1)),
                         mvtnorm::pmvnorm(upper = as.vector(z_pn1), corr = matrix(c(1,rho_1,rho_1,1), ncol=2))[[1]] )
})

