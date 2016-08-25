## profiling
library(profvis)
library(microbenchmark)

load("tests/testthat/CETvalues.rda")

microbenchmark(
  FitIsoScaleFixedICM(yTest, scaleTest,  0.1),
  FitIsoScaleFixedPG(yTest, scaleTest, 0.1),
  times = 10)

load("tests/testthat/badSimulation.rda")

startValue <- isoreg(yBadTest)$yf

microbenchmark(
  FitIsoScaleFixedICM(yBadTest, startValue,  shapeBadTest),
  FitIsoScaleFixedPG(yBadTest, startValue,  shapeBadTest),
  times = 3)

profvis({
  library(gpdIcm)
  FitIsoScaleFixedICM(yTest, scaleTest,  0.1)
})

a <- c(rep(0, 100), seq(0, 2, length.out = 300), rep(2, 100), seq(2, 5, length.out = 200), seq(5, 15, length.out = 300))
b <- seq(1, 16, length.out = 1000)
y <- c(0, runif(1000, a, b))
x <- cumsum(c(0, sample(c(0.01, 0.1, 0.3, 0.5, 0.7, 1), 1000, replace = TRUE)))
plot(x, y)

library(fdrtool)

tmp1 <- GreatestConvexMinorant(x, y)
tmp2 <- gcmlcm(x, y)

tmp1$x.knots == tmp2$x.knots
max(abs(tmp1$y.knots - tmp2$y.knots)) < 1e-6
max(abs(tmp1$y.slopes - tmp2$slope.knots)) < 1e-6

microbenchmark(
  GreatestConvexMinorant(x, y),
  gcmlcm(x, y),
  times = 100
)

