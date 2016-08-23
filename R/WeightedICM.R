#' Isotonic estimation (using an adapted version of the ICM algorithm)
#' @return isotonic scale parameter estimate and deviance
#' @inheritParams compute_nll_gpd
#' @param start Numeric vector of the initial scale parameters (will be forced to be admissible)
#' @param max_repetitions Maximal number of repetitions
#' @export
FitIsoScaleFixedICM2 <- function (y, start, shape, max_repetitions = 1e+5) {
  
  start = MakeScaleAdmissible(start, y, shape);
  
  value     = compute_nll_gpd(y, start, shape);
  scale     = start;
  new_scale = LineSearchICM(start, y, shape);
  new_value = compute_nll_gpd(y, new_scale, shape);
  
  
  i = 0;
  while ((value - new_value > 1e-6 || max(abs(scale - new_scale)) > 1e-6) && i < max_repetitions) {
    i = i + 1
    value     = new_value
    scale     = new_scale
    new_scale = LineSearchICM(new_scale, y, shape) 
    new_value = compute_nll_gpd(y, new_scale, shape)
  }
  
  nll = compute_nll_gpd(y, new_scale, shape)
  list(fitted.values = new_scale, 
       deviance = 2 * nll,
       convergence = (i < max_repetitions),
       iterations = i);
}
