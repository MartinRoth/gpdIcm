library(evd)
library(gpdIcm)

set.seed(1)
n <- 10
randomGpdGumbel  <- rgpd(n, 0, 1, 0)
randomGpdFrechet <- rgpd(n, 0, 1, 0.3)
randomGpdWeibull <- rgpd(n, 0, 1, -0.3)
randomGpdStepGumbel <- rgpd(n, 0, seq(1, 2, length.out = n), 0)

load("CETvalues.rda")



newScale <- make_gpd_admissible(scaleTest, yTest, -0.5) 
profileShape <- seq(-0.5, 0.3, by = 0.01)

context("Likelihood calculations")
test_that("Likelihood calculations", {
  expect_equal(compute_nll_gpd(randomGpdGumbel,  rep(1, n), 0),
               -sum(log(dgpd(randomGpdGumbel, 0, 1, 0))))
  expect_equal(compute_nll_gpd(randomGpdFrechet, rep(1, n), 0.3),
               -sum(log(dgpd(randomGpdFrechet, 0, 1, 0.3))))
  expect_equal(compute_nll_gpd(randomGpdWeibull, rep(1, n), -0.3),
               -sum(log(dgpd(randomGpdWeibull, 0, 1, -0.3))))
  expect_equal(compute_nll_gpd(yTest, scaleTest, -0.3),
               -sum(log(dgpd(yTest, 0, scaleTest, -0.3))))
})

context("Partial derivatives")
test_that("Partial derivative", {
  expect_equal(compute_pd1_scale_nll_gpd(randomGpdGumbel[1], 1, 0),
               1 - randomGpdGumbel[1] / 1)
  expect_equal(compute_pd1_scale_nll_gpd(randomGpdFrechet[1], 1, 0.3),
               1 - (1 + 0.3) * randomGpdFrechet[1] / (1 + 0.3 * randomGpdFrechet[1]))
  expect_equal(compute_pd1_scale_nll_gpd(randomGpdWeibull[1], 1, -0.3),
               1 - (1 - 0.3) * randomGpdWeibull[1] / (1 - 0.3 * randomGpdWeibull[1]))
})

context("Ensure admissibility")
test_that("Only admissable scale values", {
  expect_equal(make_gpd_admissible(-1, 1, 0), 1e-8)
  expect_equal(make_gpd_admissible(-1, 1, 0.1), 1e-8)
  expect_equal(make_gpd_admissible(-1, 1, -0.1), 0.1 + 1e-8)
  expect_equal(make_gpd_admissible(0.05, 1, -0.1), 0.1 + 1e-8)
  expect_equal(make_gpd_admissible(c(0.05, 0.05, 0.05), c(1, 0.9, 2), -0.1),
               c(0.1, 0.1, 0.2)+1e-8)
  expect_equal_to_reference(make_gpd_admissible(scaleTest, yTest, -0.5), "./outputTests/AdmissableScale.rds")
})

context("Isotonic fits")

test_that("GPD scale isotonic fit", {
  expect_equal_to_reference(gpd_scale_isotonic_fit(yTest, scaleTest,  0.1), "./outputTests/scaleFitFrechet.rds")
  expect_equal_to_reference(gpd_scale_isotonic_fit(yTest, scaleTest,  0.0), "./outputTests/scaleFitGumbel.rds")
  expect_equal_to_reference(gpd_scale_isotonic_fit(yTest, scaleTest, -0.1), "./outputTests/scaleFitWeibull.rds")
  expect_equal_to_reference(gpd_isotonic_scale_projected_gradient(yTest, scaleTest,  0.1), "./outputTests/scaleFitFrechet.rds")
  expect_equal_to_reference(gpd_isotonic_scale_projected_gradient(yTest, scaleTest,  0.0), "./outputTests/scaleFitGumbel.rds")
  expect_equal_to_reference(gpd_isotonic_scale_projected_gradient(yTest, scaleTest, -0.1), "./outputTests/scaleFitWeibull.rds")
})


test_that("Profile likelihood estimation", {
  expect_equal_to_reference(isotonic_scale_gpd_estimator(yTest, -0.5, 0.3), "./outputTests/ProfileLikelihoodMaximizer.rds")
  expect_error(isotonic_scale_gpd_estimator(yTest, 0.1, 0.3), "Zero should be included in the interval")
  
  yTestFrechet <- rgpd(500, 0, c(rep(1, 200), seq(1,1.1, length.out = 100), rep(1.1, 200)), 0.3)
  expect_equal_to_reference(isotonic_scale_gpd_estimator(yTestFrechet, -0.1, 0.4), "./outputTests/ProfileLikelihoodMaximizerFrechet.rds")
})

context("Failed Convergence")

load("badSimulation.rda")
test_that("Convergence fails", {
  startValue <-  isoreg(yBadTest)$yf
  tmp1 <- gpd_scale_isotonic_fit(yBadTest, startValue,  shapeBadTest)
  tmp2 <- gpd_isotonic_scale_projected_gradient(yBadTest, startValue, shapeBadTest)
  expect_false(tmp1$convergence)
  expect_false(tmp2$convergence)
})
