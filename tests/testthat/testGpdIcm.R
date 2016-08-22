library(evd)
library(gpdIcm)

set.seed(1)
n <- 10
randomGpdGumbel  <- rgpd(n, 0, 1, 0)
randomGpdFrechet <- rgpd(n, 0, 1, 0.3)
randomGpdWeibull <- rgpd(n, 0, 1, -0.3)
randomGpdStepGumbel <- rgpd(n, 0, seq(1, 2, length.out = n), 0)

load("CETvalues.rda")



newScale <- MakeScaleAdmissible(scaleTest, yTest, -0.5) 
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
  expect_equal(MakeScaleAdmissible(-1, 1, 0), 1e-8)
  expect_equal(MakeScaleAdmissible(-1, 1, 0.1), 1e-8)
  expect_equal(MakeScaleAdmissible(-1, 1, -0.1), 0.1 + 1e-8)
  expect_equal(MakeScaleAdmissible(0.05, 1, -0.1), 0.1 + 1e-8)
  expect_equal(MakeScaleAdmissible(c(0.05, 0.05, 0.05), c(1, 0.9, 2), -0.1),
               c(0.1, 0.1, 0.2)+1e-8)
  expect_equal_to_reference(MakeScaleAdmissible(scaleTest, yTest, -0.5), "./outputTests/AdmissableScale.rds")
})

context("Isotonic fits")

test_that("GPD scale isotonic fit", {
  scaleFitFrechet <- gpd_scale_isotonic_fit(yTest, scaleTest,  0.1)
  scaleFitGumbel  <- gpd_scale_isotonic_fit(yTest, scaleTest,  0.0)
  scaleFitWeibull <- gpd_scale_isotonic_fit(yTest, scaleTest, -0.1)
  expect_equal_to_reference(scaleFitFrechet, "./outputTests/scaleFitFrechet.rds")
  expect_equal_to_reference(scaleFitGumbel, "./outputTests/scaleFitGumbel.rds")
  expect_equal_to_reference(scaleFitWeibull, "./outputTests/scaleFitWeibull.rds")
  scaleFitFrechetPG <- gpd_isotonic_scale_projected_gradient(yTest, scaleTest,  0.1)
  scaleFitGumbelPG  <- gpd_isotonic_scale_projected_gradient(yTest, scaleTest,  0.0)
  scaleFitWeibullPG <- gpd_isotonic_scale_projected_gradient(yTest, scaleTest, -0.1)
  expect_equal(scaleFitFrechetPG$deviance, scaleFitFrechet$deviance)
  expect_equal(scaleFitGumbelPG$deviance,  scaleFitGumbel$deviance)
  expect_equal(scaleFitWeibullPG$deviance, scaleFitWeibull$deviance)
  expect_lt(max(abs(scaleFitFrechetPG$fitted.values - scaleFitFrechet$fitted.values)), 1e-4)
  expect_lt(max(abs(scaleFitGumbelPG$fitted.values  - scaleFitGumbel$fitted.values)),  1e-4)
  expect_lt(max(abs(scaleFitWeibullPG$fitted.values - scaleFitWeibull$fitted.values)), 1e-4)
})


test_that("Profile likelihood estimation", {
  expect_equal_to_reference(FitIsoScaleGPD(yTest, -0.5, 0.3), "./outputTests/ProfileLikelihoodMaximizer.rds")
  expect_error(FitIsoScaleGPD(yTest, 0.1, 0.3), "Zero should be included in the interval")
  
  yTestFrechet <- rgpd(500, 0, c(rep(1, 200), seq(1,1.1, length.out = 100), rep(1.1, 200)), 0.3)
  expect_equal_to_reference(FitIsoScaleGPD(yTestFrechet, -0.1, 0.4), "./outputTests/ProfileLikelihoodMaximizerFrechet.rds")
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
