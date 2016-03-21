library(evd)

n <- 10
randomGpdGumbel  <- rgpd(n, 0, 1, 0)
randomGpdFrechet <- rgpd(n, 0, 1, 0.3)
randomGpdWeibull <- rgpd(n, 0, 1, -0.3)

test_that("Likelihood calculations", {
  expect_equal(compute_nll_gpd(randomGpdGumbel,  rep(1, n), 0),
               -sum(log(dgpd(randomGpdGumbel, 0, 1, 0))))
  expect_equal(compute_nll_gpd(randomGpdFrechet, rep(1, n), 0.3),
               -sum(log(dgpd(randomGpdFrechet, 0, 1, 0.3))))
  expect_equal(compute_nll_gpd(randomGpdWeibull, rep(1, n), -0.3),
               -sum(log(dgpd(randomGpdWeibull, 0, 1, -0.3))))
})