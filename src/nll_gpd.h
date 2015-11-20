#ifndef NLL_GPD
#define NLL_GPD

double compute_nll_gpd(NumericVector y, NumericVector scale, double shape);
double compute_pd1_scale_nll_gpd(double y, double scale, double shape);

#endif