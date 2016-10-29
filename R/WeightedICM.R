LineSearchICM <- function(y, start, direction, gradient, shape) {
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
}

FitIsoScaleICMStep <- function(y, start, shape) {
  gradient <- ComputeGradient(y, start, shape)
  hesDiag  <- ComputeHessianDiagonal(y, start, shape)
  
  hesDiag[hesDiag < 0] <- pmax(-hesDiag[hesDiag < 0], 1e-5)

  points     <- cbind(c(0, cumsum(hesDiag)),
                      c(0, cumsum(start * hesDiag - gradient)))
  projection <- GreatestConvexMinorant(points[,1], points[, 2])$left.derivative
  
  direction <- projection - start
  
  
  LineSearchICM(y, start, direction, gradient, shape)
}



#' Isotonic estimation (using an adapted version of the ICM algorithm)
#' @inheritParams compute_nll_gpd
#' @param start Numeric vector of the initial scale parameters (will be forced to be admissible)
#' @param max_repetitions Maximal number of repetitions
#' @return isotonic scale parameter estimate, deviance, convergence, iterations
#' and deviance trace
#' @export
FitIsoScaleFixedICM <- function (y, start, shape, max_repetitions = 1e+5) {
  
  i     <- 0
  scale <- MakeScaleAdmissible(start, y, shape)
  value <- compute_nll_gpd(y, scale, shape) * 2
  trace <- numeric(max_repetitions)
  
  
  nextIterate <- FitIsoScaleICMStep(y, scale, shape)
  while (((value - nextIterate$deviance)/2 > 1e-6 || max(abs(scale - nextIterate$scale)) > 1e-6) && i < max_repetitions) {
    i           <- i + 1
    trace[i]    <- value
    scale       <- nextIterate$scale
    value       <- nextIterate$deviance
    nextIterate <- FitIsoScaleICMStep(y, scale, shape)
  }
  
  list(fitted.values = scale, 
       deviance      = value,
       convergence   = (i < max_repetitions),
       iterations    = i,
       trace         = trace[1:i])
}

# FitIsoScaleGPD
#' Estimation of GPD parameters with constant shape parameter and non-decreasing
#' scale parameter 
#'
#' @param y input data
#' @param min_shape double minimum shape value
#' @param max_shape double maximum shape value
#' @param by double step size for the profile likelihood
#' @param max_repetitions Integer 
#' @return isotonic scale parameter estimate and deviance
#' @importFrom stats isoreg
#' @export
FitIsoScaleGPD <- function(y, min_shape, max_shape,by = 0.01,
                           max_repetitions = 1e+5) {
  
  xi         <- generate_shape_grid(min_shape, max_shape, by)
  ny         <- length(y)
  nxi        <- length(xi)
  deviance   <- numeric(nxi)
  isoReg     <- isoreg(y)$yf
  startValue <- isoReg;
  
  posZero = which.min(abs(xi))
  
  z_best            <- FitIsoScaleFixedICM(y, startValue, xi[posZero], max_repetitions)
  best_shape        <- xi[posZero]
  deviance[posZero] <- z_best$deviance
  min_deviance      <- deviance[posZero]
  
  for(i in (posZero + 1) : nxi) {
    z           <- FitIsoScaleFixedICM(y, startValue, xi[i], max_repetitions)
    deviance[i] <- z$deviance
    startValue  <- z$fitted.values;
    if (deviance[i] < min_deviance) {
      z_best       = z
      best_shape   = xi[i]
      min_deviance = deviance[i]
    }
    i <- i + 1
  }
  
  startValue = isoReg;
  
  for(i in (posZero - 1) : 1) {
    z           <- FitIsoScaleFixedICM(y, startValue, xi[i], max_repetitions)
    deviance[i] <- z$deviance
    startValue  <- z$fitted.values;
    if (deviance[i] < min_deviance) {
      z_best       = z
      best_shape   = xi[i]
      min_deviance = deviance[i]
    }
    i <- i + 1
  }
  
  list(scale = z_best$fitted.values, 
       shape = best_shape, 
       deviance = z_best$deviance,
       convergence = z_best$convergence)
}

