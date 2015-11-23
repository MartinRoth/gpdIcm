#include <Rcpp.h>
#include <vector>
#include "rcpp_convex_minorant.h"
#include "nll_gpd.h"

using namespace std;
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
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

// [[Rcpp::export]]
NumericVector compute_next_icm_gpd(NumericVector y, NumericVector scale, double shape) {
  int ny = y.length();
  NumericVector xx(ny, 1.0); 
  NumericVector yy(ny);
  for (int i = 0; i < ny; i++) {
    yy[i] = scale[i] - compute_pd1_scale_nll_gpd(y[i], scale[i], shape);
  }
  return compute_convex_minorant_of_cumsum(xx, yy);
}

// ensure admissibility of new vectors
//' Ensure admissibility of the GPD scale parameter
//' up to now only for positive shape parameters
//'
//' @return isotonic scale parameter that does not violate the GPD constraint
//[[Rcpp::export]]
NumericVector make_gpd_admissible(NumericVector scale, NumericVector y, double shape) {
  NumericVector z(clone(scale));
  int           nz = z.length();
  
  if (shape >= 0) {
    if (z[0] <= 0 ) {
      for (int i = 0; i < nz; i++) {
        if (z[i] <= 0) z[i] = 0.0000001;
      }
    }
  }
  else {
    double current_status = 0.0000001;
    for (int i = 0; i < nz; i++) {
      if (z[i] < current_status) z[i] = current_status;
      if (z[i] < -shape * y[i]) z[i] = -shape *y[i];
      current_status = z[i];
    }
  }
  return z;
}

// search_line_icm_gpd
NumericVector search_line_icm_gpd(NumericVector y, NumericVector old_scale, NumericVector tmp_scale,
    double shape, double value) {
  
  double        t        = 1.0;
  NumericVector z        = old_scale + t * (tmp_scale - old_scale); 
  z                      = make_gpd_admissible(z, y, shape);
  double        newValue = compute_nll_gpd(y, z, shape);

  int i = 0;
  while (newValue > value && i < 30) {
    i++;
    t /= 2;
    z = old_scale + t * (tmp_scale - old_scale);
    z = make_gpd_admissible(z, y, shape);
    newValue = compute_nll_gpd(y, z, shape);
  }
  return z;  
}


// gpd_scale_isotonic_fit
//' Isotonic estimation (using an adapted version of the ICM algorithm)
//'
//' @return isotonic scale parameter estimate and deviance
//[[Rcpp::export]]
List gpd_scale_isotonic_fit (NumericVector y, NumericVector start, double shape) {

  int max_repetitions = 50000000;
  
  double        value     = compute_nll_gpd(y, start, shape);
  NumericVector tmp_scale = compute_next_icm_gpd(y, start, shape);
  NumericVector new_scale = search_line_icm_gpd(y, start, tmp_scale, shape, value);
  double        new_value = compute_nll_gpd(y, new_scale, shape);
  
  
  int i = 0;
  while (value - new_value > 0.0000000000000001 && i < max_repetitions) {
    i++;
    value     = new_value;
    tmp_scale = compute_next_icm_gpd(y, new_scale, shape);
    new_scale = search_line_icm_gpd(y, new_scale, tmp_scale, shape, value);
    new_value = compute_nll_gpd(y, new_scale, shape);
  }
  
  double nll = compute_nll_gpd(y, new_scale, shape);
  return List::create(Named("fitted.values") = new_scale, Named("deviance") = 2 * nll);
}

double compute_scalar_product (NumericVector a, NumericVector b) {
  int n = a.length();
  double scalar = 0;
  for (int i=0; i<n; i++) {
      scalar = scalar + (a[i] * b[i]);
  }
  return scalar;
}

// gpd_Goldstein_Armijo_search
NumericVector gpd_Goldstein_Armijo_search (NumericVector y, NumericVector scale, NumericVector gradient, double shape, double ll) {
  
  int    ny = y.length();
  int    max_exponent = 35;
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
  } while (min(projection) < 0);
  
  double scalar = compute_scalar_product(gradient, scale - projection);
  
  ll = compute_nll_gpd(y, projection, shape);
  
  while ((ll_old - ll < 0.0001 * scalar) && exponent < max_exponent) {
    exponent  += 1;
    yy         = scale - pow(beta, exponent) * initial_step * gradient;
    projection = compute_convex_minorant_of_cumsum(xx, yy);
    scalar     = compute_scalar_product(gradient, scale - projection);
    ll         = compute_nll_gpd(y, projection, shape);
  }
  
  return projection;
  //return List::create(Named("projection") = projection, Named("exponent") = exponent, Named("ll") = ll,
  //Named("ll_old") = ll_old, Named("scalar") = scalar);
}

// gpd_projected_gradient_next_step
//[[Rcpp::export]]
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

// gpd_isotonic_scale_projected_gradient
//' Isotonic estimation (using a projected grdient method)
//'
//' @note up to now only for positive shape parameters
//'
//' @return isotonic scale parameter estimate and deviance
//[[Rcpp::export]]
List gpd_isotonic_scale_projected_gradient (NumericVector y, NumericVector scale, double shape, int max_iterations) {
  NumericVector scale_old = scale;
  double        value     = compute_nll_gpd(y, scale, shape);
  double        new_value = compute_nll_gpd(y, scale, shape);
  int i = 0;
  do {
    value     = new_value;
    i        += 1;
    scale_old = scale;
    scale     = gpd_projected_gradient_next_step(y, scale_old, shape);
    new_value = compute_nll_gpd(y, scale, shape);
  //} while (is_true(any(scale != scale_old)) && i < max_iterations);
  } while (value - new_value > 0.0000000000000001 && i < max_iterations);
  
  double nll = compute_nll_gpd(y, scale, shape);
  return List::create(Named("fitted.values") = scale, Named("deviance") = 2 * nll);
}

// isotonic_scale_gpd_estimator
//' Estimation of GPD parameters with fixed shape parameter and non-decreasing scale parameter 
//'
//'
//' @return isotonic scale parameter estimate and deviance
//[[Rcpp::export]]
List isotonic_scale_gpd_estimator (NumericVector y, NumericVector xi) {
  
  int           ny = y.length();
  int           nxi = xi.length();
  NumericVector xx(ny, 1.0);
  NumericVector log_likelihood(nxi);
  Rcpp::List    z;
  Rcpp::List    z_best;
  double        best_shape;
  double        max_log_likelihood;
  NumericVector isoReg = compute_convex_minorant_of_cumsum(xx, y);
  
  z_best             = gpd_scale_isotonic_fit(y, isoReg, xi[0]);
  best_shape         = xi[0];
  log_likelihood[0]  = -(float)z_best["deviance"]/2.0;
  max_log_likelihood = log_likelihood[0];
  
  for(int i = 1; i < nxi; i++) {
    z = gpd_scale_isotonic_fit(y, isoReg, xi[i]);
    log_likelihood[i] = -(float)z["deviance"]/2.0;
    if (log_likelihood[i] > max_log_likelihood) {
      z_best = z;
      best_shape = xi[i];
      max_log_likelihood = log_likelihood[i];
    }
  }
  
  return List::create(Named("scale") = z_best["fitted.values"], 
    Named("shape") = best_shape, 
    Named("deviance") = z_best["deviance"]);
}

