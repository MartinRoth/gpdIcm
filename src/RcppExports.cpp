// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_convex_minorant_of_cumsum
NumericVector compute_convex_minorant_of_cumsum(NumericVector x, NumericVector y);
RcppExport SEXP gpdIcm_compute_convex_minorant_of_cumsum(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = Rcpp::wrap(compute_convex_minorant_of_cumsum(x, y));
    return __result;
END_RCPP
}
// compute_next_icm_gpd
NumericVector compute_next_icm_gpd(NumericVector y, NumericVector scale, double shape);
RcppExport SEXP gpdIcm_compute_next_icm_gpd(SEXP ySEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(compute_next_icm_gpd(y, scale, shape));
    return __result;
END_RCPP
}
// make_gpd_admissible
NumericVector make_gpd_admissible(NumericVector scale, NumericVector y, double shape);
RcppExport SEXP gpdIcm_make_gpd_admissible(SEXP scaleSEXP, SEXP ySEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(make_gpd_admissible(scale, y, shape));
    return __result;
END_RCPP
}
// search_line_icm_gpd
NumericVector search_line_icm_gpd(NumericVector y, NumericVector old_scale, NumericVector tmp_scale, double shape, double value);
RcppExport SEXP gpdIcm_search_line_icm_gpd(SEXP ySEXP, SEXP old_scaleSEXP, SEXP tmp_scaleSEXP, SEXP shapeSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type old_scale(old_scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmp_scale(tmp_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type value(valueSEXP);
    __result = Rcpp::wrap(search_line_icm_gpd(y, old_scale, tmp_scale, shape, value));
    return __result;
END_RCPP
}
// gpd_scale_isotonic_fit
List gpd_scale_isotonic_fit(NumericVector y, NumericVector start, double shape);
RcppExport SEXP gpdIcm_gpd_scale_isotonic_fit(SEXP ySEXP, SEXP startSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(gpd_scale_isotonic_fit(y, start, shape));
    return __result;
END_RCPP
}
// gpd_projected_gradient_next_step
NumericVector gpd_projected_gradient_next_step(NumericVector y, NumericVector scale, double shape);
RcppExport SEXP gpdIcm_gpd_projected_gradient_next_step(SEXP ySEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(gpd_projected_gradient_next_step(y, scale, shape));
    return __result;
END_RCPP
}
// gpd_isotonic_scale_projected_gradient
List gpd_isotonic_scale_projected_gradient(NumericVector y, NumericVector scale, double shape);
RcppExport SEXP gpdIcm_gpd_isotonic_scale_projected_gradient(SEXP ySEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(gpd_isotonic_scale_projected_gradient(y, scale, shape));
    return __result;
END_RCPP
}
// isotonic_scale_gpd_estimator
List isotonic_scale_gpd_estimator(NumericVector y, NumericVector shape);
RcppExport SEXP gpdIcm_isotonic_scale_gpd_estimator(SEXP ySEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    __result = Rcpp::wrap(isotonic_scale_gpd_estimator(y, shape));
    return __result;
END_RCPP
}
// compute_nll_gpd
double compute_nll_gpd(NumericVector y, NumericVector scale, double shape);
RcppExport SEXP gpdIcm_compute_nll_gpd(SEXP ySEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(compute_nll_gpd(y, scale, shape));
    return __result;
END_RCPP
}
// compute_pd1_scale_nll_gpd
double compute_pd1_scale_nll_gpd(double y, double scale, double shape);
RcppExport SEXP gpdIcm_compute_pd1_scale_nll_gpd(SEXP ySEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    __result = Rcpp::wrap(compute_pd1_scale_nll_gpd(y, scale, shape));
    return __result;
END_RCPP
}
// convexMinorant
NumericVector convexMinorant(NumericVector x, NumericVector y);
RcppExport SEXP gpdIcm_convexMinorant(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = Rcpp::wrap(convexMinorant(x, y));
    return __result;
END_RCPP
}
