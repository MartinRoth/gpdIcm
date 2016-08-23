## profiling
library(profvis)
library(microbenchmark)

load("tests/testthat/CETvalues.rda")

microbenchmark(
  FitIsoScaleFixedICM(yTest, scaleTest,  0.1),
  FitIsoScaleFixedICM2(yTest, scaleTest,  0.1),
  FitIsoScaleFixedPG(yTest, scaleTest, 0.1),
  times = 10)

load("tests/testthat/badSimulation.rda")

startValue <- isoreg(yBadTest)$yf

microbenchmark(
  FitIsoScaleFixedICM(yBadTest, startValue,  shapeBadTest),
  FitIsoScaleFixedPG(yBadTest, startValue,  shapeBadTest),
  times = 3)

# profile FitIsoScaleFixedICM2

profvis({
  library(gpdIcm)
  FitIsoScaleFixedICM2(yTest, scaleTest,  0.1)
})
