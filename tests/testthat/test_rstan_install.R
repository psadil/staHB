context("Rstan Installation")

fx <- inline::cxxfunction( signature(x = "integer", y = "numeric" ) ,
                           'return ScalarReal( INTEGER(x)[0] * REAL(y)[0] ) ;' )

testthat::test_that("# should be 10", {
  testthat::expect_equal(fx(2L,5), 10)
})
