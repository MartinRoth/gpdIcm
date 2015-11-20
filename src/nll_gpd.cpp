#include <Rcpp.h>
#include <vector>

using namespace std;
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//' Computes the negative log likelihood for the GPD
//' 
//' No threshold parameter
//' 
//' scale has same length as y and shape has length 1
// [[Rcpp::export]]
double compute_nll_gpd(NumericVector y, NumericVector scale, double shape) {
  int ny      = y.length();
  int n_scale = scale.length();
  
  if (ny != n_scale) stop("Error: Length y and length scale are different.");
  
  double nll = 0.0;
  
  //if (abs(shape) < 0.00000001) {
  if (shape == 0) {
    for (int i = 0; i < ny; i++) {
      nll += log(scale[i]) + y[i] / scale[i];
    }
  }
  else {
    double a = (1 + shape) / shape;
    for (int i = 0; i < ny; i++) {
      nll += log(scale[i]) + a * log(1 + shape * y[i] / scale[i]);
    }
  }
  return(nll);
}

//' Computes the first partial derivative (scale) of the negative log likelihood for the GPD 
//' 
//' No threshold parameter
//' 
//' scale has same length as y and shape has length 1
// [[Rcpp::export]]
double compute_pd1_scale_nll_gpd(double y, double scale, double shape) {
  return ( 1 / scale - (1 + shape) * y / (pow(scale, 2) * (1 + shape * y / scale)));
}


