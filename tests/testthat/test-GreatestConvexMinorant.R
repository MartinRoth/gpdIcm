context("GreatestConvexMinorant")

## TODO: Rename context
## TODO: Add more tests

library(fdrtool)

test_that("Left derivative of gcm", {
  x <- c(0, 0.5 , 1,   2   , 2.5, 3,   3.5)
  y <- c(0, 0.75, 0.5, 1.75, 1.5, 2.5, 2.5)
  expect_equal(LeftDerivativeGCM(x, y), c(0.5, 0.5, 2/3, 2/3, 1, 1), tolerance=1e-6)
})
