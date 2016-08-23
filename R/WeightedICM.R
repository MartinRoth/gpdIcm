#' Left GreatestConvexMinorant
#' @param x Vector of x values
#' @param y Vector of y values
#' @export
#' @importFrom fdrtool gcmlcm
LeftDerivativeGCM <- function(x, y) {
  fit <- gcmlcm(x, y)
  LeftDerivative(x[-1], fit)
}

#' @importFrom purrr map_int
LeftDerivative <- function(x, fit) {
  knotsCount <- length(fit$x.knots)
  knotsIndex <- 1 : knotsCount
  fn <- function(x) {
    as.integer(max(knotsIndex[!fit$x.knots >= x]))
  }
  index <- map_int(x, fn) 
  return(fit$slope.knots[index])
}

# double beta = 0.5;
# double step = 1.0;
# double mu = 1e-4;
# int    e = 0;
# int    maxExponent = 25;
# double a = step;
# 
# NumericVector gradient = ComputeGradient(y, oldScale, shape); 
# 
# double nllOld = compute_nll_gpd(y, oldScale, shape);
# NumericVector direction = compute_next_icm_gpd(y, oldScale, shape) - oldScale;
# NumericVector newScale = MakeScaleAdmissible(oldScale + a * direction, y, shape);
# double nllNew = compute_nll_gpd(y, newScale, shape);
# double scalar = compute_scalar_product(gradient, oldScale - newScale);
# while ( (e < maxExponent) && (nllOld - nllNew < mu * scalar)) {
#   e += 1;
#   a *= beta;
#   newScale = MakeScaleAdmissible(oldScale + a * direction, y, shape);
#   nllNew = compute_nll_gpd(y, newScale, shape);
#   scalar = compute_scalar_product(gradient, oldScale - newScale);
# }
# return newScale;

#' @importFrom purrr map_dbl
LineSearchICM2 <- function(y, start, direction, gradient, shape) {
  beta        <- 0.5
  step        <- 1.0
  mu          <- 1e-4
  exponent    <- 0
  maxExponent <- 25
  alpha       <- step
  
  nll      <- compute_nll_gpd(y, start, shape)
  newScale <- MakeScaleAdmissible(start + alpha * direction, y, shape)
  newNll   <- compute_nll_gpd(y, newScale, shape);
  scalar   <- compute_scalar_product(gradient, start - newScale);
  
  while( (exponent < maxExponent) && (nll - newNll < mu * scalar)) {
    exponent <- exponent + 1
    alpha    <- alpha * beta
    newScale <- MakeScaleAdmissible(start + alpha * direction, y, shape)
    newNll   <- compute_nll_gpd(y, newScale, shape);
    scalar   <- compute_scalar_product(gradient, start - newScale);
  }
  
  return(list(scale = newScale, deviance = 2 * newNll))
  # alphas <- sort(c(seq(0.1, 1, by = 0.1), 5/10^(c(1:25)), 1/10^(c(2:25))),
  #                decreasing= TRUE)
  # 
  # fn <- function(alpha) {
  #   tmp_scale = MakeScaleAdmissible(start + alpha * direction, y, shape)
  #   compute_nll_gpd(y, tmp_scale, shape)
  # }
  # nll <- map_dbl(alphas, fn)
  # 
  # minIndex <- which.min(nll)
  # newScale <- MakeScaleAdmissible(start + alphas[minIndex] * direction,
  #                                          y, shape)
  # return(list(scale = newScale, deviance = 2 * nll[minIndex]))
}

#' @importFrom fdrtool gcmlcm
FitIsoScaleICMStep2 <- function(y, start, shape) {
  gradient <- ComputeGradient(y, start, shape)
  hesDiag  <- ComputeHessianDiagonal(y, start, shape)
  
  #posSubset           <- hesDiag > 0
  #negativeWeights     <- -mean(hesDiag[!posSubset])
  #hesDiag[!posSubset] <-  - hesDiag[!posSubset] # negativeWeights #
  #hesDiag <- rep(1, length(hesDiag))
  hesDiag[hesDiag < 0]    <- -hesDiag[hesDiag < 0]
  hesDiag[hesDiag < 1e-6] <- 1e-6
  
  points     <- cbind(c(0, cumsum(hesDiag)),
                      c(0, cumsum(start * hesDiag - gradient)))
  projection <- LeftDerivativeGCM(points[,1], points[, 2])
  
  direction <- projection - start
  
  LineSearchICM2(y, start, direction, gradient, shape)
}
  


#' Isotonic estimation (using an adapted version of the ICM algorithm)
#' @return isotonic scale parameter estimate and deviance
#' @inheritParams compute_nll_gpd
#' @param start Numeric vector of the initial scale parameters (will be forced to be admissible)
#' @param max_repetitions Maximal number of repetitions
#' @export
FitIsoScaleFixedICM2 <- function (y, start, shape, max_repetitions = 1e+5) {
  
  scale <- MakeScaleAdmissible(start, y, shape)
  value <- compute_nll_gpd(y, scale, shape) * 2
  
  nextIterate <- FitIsoScaleICMStep2(y, scale, shape)
  
  i <- 0;
  while ((value - nextIterate$deviance > 1e-6 || max(abs(scale - nextIterate$scale)) > 1e-6) && i < max_repetitions) {
    i = i + 1
    scale <- nextIterate$scale
    value <- nextIterate$deviance
    nextIterate <- FitIsoScaleICMStep2(y, scale, shape)
  }
  
  list(fitted.values = scale, 
       deviance = value,
       convergence = (i < max_repetitions),
       iterations = i);
}


