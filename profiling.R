## profiling

library(microbenchmark)

load("tests/testthat/CETvalues.rda")

microbenchmark(
  gpd_scale_isotonic_fit(yTest, scaleTest,  0.1),
  gpd_isotonic_scale_projected_gradient(yTest, scaleTest, 0.1),
  times = 10)

load("tests/testthat/badSimulation.rda")

startValue <- isoreg(yBadTest)$yf

microbenchmark(
  gpd_scale_isotonic_fit(yBadTest, startValue,  shapeBadTest),
  gpd_isotonic_scale_projected_gradient(yBadTest, startValue,  shapeBadTest),
  times = 3)
