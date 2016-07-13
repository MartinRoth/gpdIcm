#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "rcpp_convex_minorant.h"
#include "nll_gpd.h"

using namespace std;
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//// [[Rcpp::export]]
NumericVector compute_convex_minorant_of_cumsum(NumericVector x, NumericVector y) {
  vector<double> xx;
  vector<double> yy;
  
  xx.push_back(0);
  yy.push_back(0);
  partial_sum(x.begin(), x.end(), back_inserter(xx));
  partial_sum(y.begin(), y.end(), back_inserter(yy));
  
  NumericVector xxx = Rcpp::wrap(xx);
  NumericVector yyy = Rcpp::wrap(yy);
  
  NumericVector newBeta = diff(convexMinorant(xxx, yyy));
  return newBeta;
}



// ensure admissibility of new vectors
//' Ensure admissibility of the GPD scale parameter
//' up to now only for positive shape parameters
//'
//' @inheritParams compute_nll_gpd
//' @return isotonic scale parameter that does not violate the GPD constraint
//[[Rcpp::export]]
NumericVector make_gpd_admissible(NumericVector scale, NumericVector y, double shape) {
  NumericVector z(clone(scale));
  int           nz = z.length();
  
  if (shape >= 0) {
    if (z[0] <= 0 ) {
      for (int i = 0; i < nz; i++) {
        if (z[i] <= 0) z[i] = 1e-8;
      }
    }
  }
  else {
    double current_status = 0.00000001;
    for (int i = 0; i < nz; i++) {
      if (z[i] < current_status) z[i] = current_status;
      if (z[i] <= -shape * y[i]) z[i] = -shape *y[i]+1e-8;
      current_status = z[i];
    }
  }
  return z;
}

//// [[Rcpp::export]]
NumericVector compute_next_icm_gpd(NumericVector y, NumericVector scale, double shape) {
  int ny = y.length();
  NumericVector xx(ny, 1.0); 
  NumericVector yy(ny);
  for (int i = 0; i < ny; i++) {
    yy[i] = scale[i] - compute_pd1_scale_nll_gpd(y[i], scale[i], shape);
  }
  return compute_convex_minorant_of_cumsum(xx, yy);
}

//// search_line_icm_gpd
////[[Rcpp::export]]
//NumericVector search_line_icm_gpd(NumericVector y, NumericVector old_scale, NumericVector tmp_scale,
//    double shape, double value) {
//  
//  double        t        = 1.0;
//  NumericVector z        = old_scale + t * (tmp_scale - old_scale); 
//  z                      = make_gpd_admissible(z, y, shape);
//  double        newValue = compute_nll_gpd(y, z, shape);
//
//  int i = 0;
//  while (newValue > value && i < 30) {
//    i++;
//    t /= 2;
//    z = old_scale + t * (tmp_scale - old_scale);
//    z = make_gpd_admissible(z, y, shape);
//    newValue = compute_nll_gpd(y, z, shape);
//  }
//  return z;  
//}

double compute_scalar_product (NumericVector a, NumericVector b) {
  int n = a.length();
  double scalar = 0;
  for (int i=0; i<n; i++) {
      scalar = scalar + (a[i] * b[i]);
  }
  return scalar;
}

NumericVector ComputeGradient (NumericVector y, NumericVector scale, double shape) {
  int           ny = y.length();
  NumericVector gradient(ny);
  
  for (int i = 0; i < ny; i++) {
    gradient[i] = compute_pd1_scale_nll_gpd(y[i], scale[i], shape);
  }
   return gradient;
}

// Goldstein-Armijo search for ICM method
////[[Rcpp::export]]
NumericVector lineSearchICM (NumericVector oldScale, NumericVector y, double shape) {
  double b = 0.5;
  double s = 1.0;
  double m = 1e-4;
  int    e = 0;
  int    maxExponent = 25;
  double a = s;
  
  NumericVector gradient = ComputeGradient(y, oldScale, shape); 
  
  double nllOld = compute_nll_gpd(y, oldScale, shape);
  NumericVector direction = compute_next_icm_gpd(y, oldScale, shape) - oldScale;
  NumericVector newScale = make_gpd_admissible(oldScale + a * direction, y, shape);
  double nllNew = compute_nll_gpd(y, newScale, shape);
  double scalar = compute_scalar_product(gradient, oldScale - newScale);
  while ( (e < maxExponent) && (nllOld - nllNew < m * scalar)) {
    e += 1;
    a *= b;
    newScale = make_gpd_admissible(oldScale + a * direction, y, shape);
    nllNew = compute_nll_gpd(y, newScale, shape);
    scalar = compute_scalar_product(gradient, oldScale - newScale);
  }
  return newScale;
}

// gpd_Goldstein_Armijo_search
NumericVector gpd_Goldstein_Armijo_search (NumericVector y, NumericVector scale, NumericVector gradient, double shape, double ll) {
  
  int    ny = y.length();
  int    max_exponent = 25;
  int    exponent = 0;
  double beta = 0.5;
  double initial_step = 2;
  double ll_old = ll;
  
  NumericVector xx(ny, 1.0);
  NumericVector yy(ny, 0.0);
  NumericVector projection(ny, 0.0);
  
  do {
    exponent += 1;
    yy = scale - pow(beta, exponent) * initial_step * gradient;
    projection = compute_convex_minorant_of_cumsum(xx, yy);
  } while (min(projection) < 0); // use make gpd admissible instead
  
  double scalar = compute_scalar_product(gradient, scale - projection);
  
  ll = compute_nll_gpd(y, projection, shape);
  
  while ((ll_old - ll < 1e-4 * scalar) && exponent < max_exponent) {
    exponent  += 1;
    yy         = scale - pow(beta, exponent) * initial_step * gradient;
    projection = make_gpd_admissible(compute_convex_minorant_of_cumsum(xx, yy), y, shape);
    scalar     = compute_scalar_product(gradient, scale - projection);
    ll         = compute_nll_gpd(y, projection, shape);
  }
  
  return projection;
  //return List::create(Named("projection") = projection, Named("exponent") = exponent, Named("ll") = ll,
  //Named("ll_old") = ll_old, Named("scalar") = scalar);
}


// gpd_projected_gradient_next_step
NumericVector gpd_projected_gradient_next_step (NumericVector y, NumericVector scale, double shape) {
  
  int           ny = y.length();
  double        ll = compute_nll_gpd(y, scale, shape);
  NumericVector scale_old = scale;
  NumericVector gradient(ny);
  
  for (int i = 0; i < ny; i++) {
    gradient[i] = compute_pd1_scale_nll_gpd(y[i], scale[i], shape);
  }
  
  return gpd_Goldstein_Armijo_search(y, scale, gradient, shape, ll);
}

// gpd_scale_isotonic_fit
//' Isotonic estimation (using an adapted version of the ICM algorithm)
//'
//' @return isotonic scale parameter estimate and deviance
//' @inheritParams compute_nll_gpd
//' @param start Numeric vector of the initial scale parameters (will be forced to be admissible)
//' @param max_repetitions Maximal number of repetitions
//' @useDynLib gpdIcm
//' @importFrom Rcpp evalCpp
//' @export
//[[Rcpp::export]]
List gpd_scale_isotonic_fit (NumericVector y, NumericVector start, double shape, int max_repetitions = 1e+5) {

  start = make_gpd_admissible(start, y, shape);
  
  double        value     = compute_nll_gpd(y, start, shape);
  NumericVector scale     = start;
  NumericVector new_scale = lineSearchICM(start, y, shape);
  double        new_value = compute_nll_gpd(y, new_scale, shape);
  
  
  int i = 0;
  while ((value - new_value > 1e-6 || max(abs(scale - new_scale)) > 1e-6) && i < max_repetitions) {
    i++;
    value     = new_value;
    scale     = new_scale;
    new_scale = lineSearchICM(new_scale, y, shape); 
    new_value = compute_nll_gpd(y, new_scale, shape);
  }
  
  double nll = compute_nll_gpd(y, new_scale, shape);
  return List::create(Named("fitted.values") = new_scale,
                      Named("deviance") = 2 * nll,
                      Named("convergence") = (i < max_repetitions));
}


// gpd_isotonic_scale_projected_gradient
//' Isotonic estimation (using a projected grdient method)
//'
//' @note up to now only for positive shape parameters
//'
//' @inheritParams compute_nll_gpd
//' @inheritParams gpd_scale_isotonic_fit
//' @return isotonic scale parameter estimate and deviance
//' @export
//[[Rcpp::export]]
List gpd_isotonic_scale_projected_gradient (NumericVector y, NumericVector scale, double shape, int max_repetitions = 1e+5) {
  scale = make_gpd_admissible(scale, y, shape);
  
  NumericVector scale_new = scale;
  double        value     = compute_nll_gpd(y, scale, shape);
  double        new_value = compute_nll_gpd(y, scale, shape);
  int i = 0;
  do {
    value     = new_value;
    i        += 1;
    scale     = scale_new;
    scale_new = gpd_projected_gradient_next_step(y, scale, shape);
    new_value = compute_nll_gpd(y, scale_new, shape);
  } while ((value - new_value > 1e-6 || max(abs(scale - scale_new)) > 1e-6) && i < max_repetitions);
  
  double nll = compute_nll_gpd(y, scale_new, shape);
  return List::create(Named("fitted.values") = scale_new,
                      Named("deviance") = 2 * nll,
                      Named("convergence") = (i < max_repetitions));
}

struct add_multiple {
  int incr;
  int count;
  add_multiple(int incr)
    : incr(incr), count(0)
  {}
  inline int operator()(int d) {
    return d + incr * count++;
  }
};

// obtained from http://stackoverflow.com/questions/31336548/rcpp-version-of-base-r-seq-drops-values
NumericVector rcpp_seq(double from_, double to_, double by_ = 1.0) {
  int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
  int from = adjust * from_;
  int to = adjust * to_;
  int by = adjust * by_;
  
  std::size_t n = ((to - from) / by) + 1;
  Rcpp::IntegerVector res = Rcpp::rep(from, n);
  add_multiple ftor(by);
  
  std::transform(res.begin(), res.end(), res.begin(), ftor);
  return Rcpp::NumericVector(res) / adjust;
}

NumericVector generate_shape_grid(double from_, double to_, double by_ = 0.01) {
  if (from_ >= 0 || to_ <= 0) {
    stop("Zero should be included in the interval");
  }
  return rcpp_seq(from_, to_, by_);
}

// isotonic_scale_gpd_estimator
//' Estimation of GPD parameters with fixed shape parameter and non-decreasing
//' scale parameter 
//'
//' @inheritParams compute_nll_gpd
//' @inheritParams gpd_scale_isotonic_fit
//' @param min_shape double minimum shape value
//' @param max_shape double maximum shape value
//' @param by double step size for the profile likelihood 
//' @return isotonic scale parameter estimate and deviance
//' @export
//[[Rcpp::export]]
List isotonic_scale_gpd_estimator (NumericVector y, double min_shape,
                                   double max_shape, double by = 0.01,
                                   int max_repetitions = 1e+5) {
  
  // NumericVector xi = shape;
  NumericVector xi = generate_shape_grid(min_shape, max_shape, by);
  int           ny = y.length();
  int           nxi = xi.length();
  NumericVector xx(ny, 1.0);
  NumericVector log_likelihood(nxi);
  Rcpp::List    z;
  Rcpp::List    z_best;
  double        best_shape;
  double        max_log_likelihood;
  bool          convergence;
  NumericVector isoReg = compute_convex_minorant_of_cumsum(xx, y);
  NumericVector startValue = isoReg;
  
  int posZero = which_min(abs(xi));
  
  z_best             = gpd_scale_isotonic_fit(y, startValue, xi[posZero], max_repetitions); 
  best_shape         = xi[posZero];
  log_likelihood[posZero]  = -(float)z_best["deviance"]/2.0;
  max_log_likelihood = log_likelihood[posZero];
  
  for(int i = posZero + 1; i < nxi; i++) {
    z = gpd_scale_isotonic_fit(y, startValue, xi[i], max_repetitions);
    log_likelihood[i] = -(float)z["deviance"] / 2.0;
    startValue = z["fitted.values"];
    if (log_likelihood[i] > max_log_likelihood) {
      z_best = z;
      best_shape = xi[i];
      max_log_likelihood = log_likelihood[i];
    }
  }
  
  startValue = isoReg;
  
  for(int i = posZero - 1; i >= 0; i--) {
    z = gpd_scale_isotonic_fit(y, startValue, xi[i], max_repetitions);
    log_likelihood[i] = -(float)z["deviance"]/2.0;
    startValue = z["fitted.values"];
    if (log_likelihood[i] > max_log_likelihood) {
      z_best = z;
      best_shape = xi[i];
      max_log_likelihood = log_likelihood[i];
    }
  }
  
  return List::create(Named("scale") = z_best["fitted.values"], 
    Named("shape") = best_shape, 
    Named("deviance") = z_best["deviance"],
    Named("convergence") = z_best["convergence"]);
    // Named("posZero") = posZero);
}