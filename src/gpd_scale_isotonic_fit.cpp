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
//' @useDynLib gpdIcm
//[[Rcpp::export]]
NumericVector MakeScaleAdmissible(NumericVector scale, NumericVector y, double shape) {
  NumericVector z(clone(scale));
  int           nz = z.length();
  
  double step = 1e-8;
  double current = step;
  
  for (int i = 0; i < nz; i++) {
    if (z[i] < current) z[i] = current;
    else current = z[i];
    if (shape < 0) {
      if (z[i] <= -shape * y[i]) {
        current = -shape * y[i] + step;
        z[i] = current;
      }
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

//[[Rcpp::export]]
double compute_scalar_product (NumericVector a, NumericVector b) {
  int n = a.length();
  double scalar = 0;
  for (int i=0; i<n; i++) {
      scalar = scalar + (a[i] * b[i]);
  }
  return scalar;
}

//[[Rcpp::export]]
NumericVector ComputeGradient (NumericVector y, NumericVector scale, double shape) {
  int           ny = y.length();
  NumericVector gradient(ny);
  
  for (int i = 0; i < ny; i++) {
    gradient[i] = compute_pd1_scale_nll_gpd(y[i], scale[i], shape);
  }
   return gradient;
}

//[[Rcpp::export]]
NumericVector ComputeHessianDiagonal (NumericVector y, NumericVector scale, double shape) {
  int           ny = y.length();
  NumericVector hesDiag(ny);
  
  for (int i = 0; i < ny; i++) {
    hesDiag[i] = compute_pd2_scale_nll_gpd(y[i], scale[i], shape);
  }
  return hesDiag;
}


// LineSearchPG
NumericVector LineSearchPG (NumericVector y, NumericVector scale, NumericVector gradient, double shape, double ll) {
  // Goldstein-Armijo type choice of scaling factor
  int    ny = y.length();
  int    max_exponent = 31;
  int    exponent = 0;
  double beta = 0.5;
  double initial_step = 128;
  double ll_old = ll;
  double mu = 1e-4;
  
  NumericVector xx(ny, 1.0);
  NumericVector yy(ny, 0.0);
  NumericVector projection(ny, 0.0);
  
  exponent += 1;
  yy = scale - pow(beta, exponent) * initial_step * gradient;
  projection = MakeScaleAdmissible(compute_convex_minorant_of_cumsum(xx, yy), y, shape);
  
  double scalar = compute_scalar_product(gradient, scale - projection);
  
  ll = compute_nll_gpd(y, projection, shape);
  
  while ((ll_old - ll < mu * scalar) && exponent < max_exponent) {
    exponent  += 1;
    yy         = scale - pow(beta, exponent) * initial_step * gradient;
    projection = MakeScaleAdmissible(compute_convex_minorant_of_cumsum(xx, yy), y, shape);
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
  
  return LineSearchPG(y, scale, gradient, shape, ll);
}

// FitIsoScaleFixedPG
//' Isotonic estimation (using a projected gradient method)
//'
//' @note up to now only for positive shape parameters
//'
//' @inheritParams compute_nll_gpd
//' @inheritParams FitIsoScaleFixedICM
//' @return isotonic scale parameter estimate and deviance
//' @export
//[[Rcpp::export]]
List FitIsoScaleFixedPG (NumericVector y, NumericVector scale, double shape, int max_repetitions = 1e+5) {
  scale = MakeScaleAdmissible(scale, y, shape);
  
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
                      Named("convergence") = (i < max_repetitions),
                      Named("iterations") = i);
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

//[[Rcpp::export]]
NumericVector generate_shape_grid(double from_, double to_, double by_ = 0.01) {
  if (from_ >= 0 || to_ <= 0) {
    stop("Zero should be included in the interval");
  }
  return rcpp_seq(from_, to_, by_);
}
