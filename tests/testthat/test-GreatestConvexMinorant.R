context("GreatestConvexMinorant")

## TODO: Rename context
## TODO: Add more tests

test_that("Left derivative of gcm", {
  x      <- c(0, 0.5 , 1,   2   , 2.5, 3,   3.5)
  y      <- c(0, 0.75, 0.5, 1.75, 1.5, 2.5, 2.5)
  result <- c(0.5, 0.5, 2/3, 2/3, 1, 1)
  expect_equal(GreatestConvexMinorant(x, y)$left.derivative, result)
  x      <- c(0, 0.5, 1,   1,   2)
  y      <- c(0, 1,   0.5, 1.5, 1.5)
  result <- c(0.5, 0.5, 0.5, 1)
  expect_equal(GreatestConvexMinorant(x, y)$left.derivative, result)
  x      <- c(0, 0,   1, 2)
  y      <- c(0, 0.5, 1, 1.5)
  result <- c(0, 0.75, 0.75)
  expect_equal(GreatestConvexMinorant(x, y)$left.derivative, result)
  x      <- c(0, 1, 2, 2)
  y      <- c(0, 1, 1, 2)
  result <- c(0.5, 0.5, 0.5)
  expect_equal(GreatestConvexMinorant(x, y)$left.derivative, result)
})
