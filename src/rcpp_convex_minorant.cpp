#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include "convexHull.h"


using namespace std;
using namespace Rcpp;

//Computes the slopes of the convex minorant given the points of the lower convex hull
//The arguments could be changed to double vectors
NumericVector compute_slopes(NumericVector x, NumericVector y) {
  double        startValue = y[0];
  NumericVector xDiff      = diff(x);
  NumericVector yDiff      = diff(y);
  int           nxd        = xDiff.length();

  vector<double> tmp;
  tmp.push_back(startValue);
  for (int i = 0; i < nxd; i++) {
    for (int j = 0; j < xDiff[i]; j++) {
      tmp.push_back( yDiff[i] / xDiff[i]);
    }
  }
  vector<double> tmp2;
  partial_sum(tmp.begin(), tmp.end(), back_inserter(tmp2));
  return Rcpp::wrap(tmp2);
}

//' Computes the convex minorant of a polygon.
//' @param x,y the coordinates of the polygon
//' @return vector of the y-coordinates of the convex minorant
NumericVector convexMinorant(NumericVector x, NumericVector y) {
  
  int ny = y.length();

  NumericVector XX = x; 
  NumericVector XY = y; 
  
  vector<Point> P(ny);
  for (int i = 0; i < ny; i++) {
    P[i].x = XX[i];
    P[i].y = XY[i];
  }
  
  vector<Point> convHull = convex_hull(P);
  
  int            nP = convHull.size();
  vector<int>    convHullX(nP); 
  vector<double> convHullY(nP); 
  for (int i = 0; i < nP; i++) {
    convHullX[i] = convHull.at(i).x;
    convHullY[i] = convHull.at(i).y;
  }
  
  NumericVector XXX     = Rcpp::wrap(convHullX);
  NumericVector XYY     = Rcpp::wrap(convHullY);
  NumericVector slopes  = compute_slopes(XXX, XYY);
  return slopes;
  //return List::create(Named("slopes") = newBeta, Named("x") = XXX, Named("y") = XYY);
}