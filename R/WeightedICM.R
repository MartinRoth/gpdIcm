LeftDerivativeGCM <- function(x, fit) {
  knotsCount <- length(fit$x.knots)
  knotsIndex <- 1 : knotsCount
  index <- purrr:::map_int(x, function(x) {
    return(as.integer(min(c(max(knotsIndex[!fit$x.knots > x]), knotsCount-1))))})
  return(fit$slope.knots[index])
}

#' @import foreach
#' @import iterators
LineSearchICM2 <- function(y, start, direction, shape) {
  alphas <- sort(c(seq(0.1, 1, by = 0.1), 5/10^(c(1:25)), 1/10^(c(2:25))),
                 decreasing= TRUE)
  
  nll <- foreach(alpha = iter(alphas), .combine = "c") %do% {
    tmp_scale = MakeScaleAdmissible(start + alpha * direction, y, shape)
    compute_nll_gpd(y, tmp_scale, shape)
  }
  
  minIndex <- which.min(nll)
  newScale <- MakeScaleAdmissible(start + alphas[minIndex] * direction,
                                           y, shape)
  return(list(scale = newScale, deviance = 2 * nll[minIndex]))
}

#' @importFrom fdrtool gcmlcm
FitIsoScaleICMStep2 <- function(y, start, shape) {
  gradient <- ComputeGradient(y, start, shape)
  hesDiag  <- ComputeHessianDiagonal(y, start, shape)
  
  posSubset           <- hesDiag > 1e-6
  negativeWeights     <- -mean(hesDiag[!posSubset])
  hesDiag[!posSubset] <-  negativeWeights # - hesDiag[!posSubset] # 
  #hesDiag <- rep(1, length(hesDiag))
  
  points     <- cbind(c(0, cumsum(hesDiag)),
                      c(0, cumsum(start * hesDiag - gradient)))
  fit        <- fdrtool::gcmlcm(points[,1], points[, 2])
  projection <- LeftDerivativeGCM(points[-1, 1], fit)
  
  direction <- projection - start
  
  LineSearchICM2(y, start, direction, shape)
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
  while ((value - nextIterate$deviance > 1e-15 || max(abs(scale - nextIterate$scale)) > 1e-15) && i < max_repetitions) {
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
