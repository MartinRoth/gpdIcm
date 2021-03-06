% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WeightedICM.R
\name{FitIsoScaleGPD}
\alias{FitIsoScaleGPD}
\title{Estimation of GPD parameters with constant shape parameter and non-decreasing
scale parameter}
\usage{
FitIsoScaleGPD(y, min_shape, max_shape, by = 0.01, max_repetitions = 1e+05)
}
\arguments{
\item{y}{input data}

\item{min_shape}{double minimum shape value}

\item{max_shape}{double maximum shape value}

\item{by}{double step size for the profile likelihood}

\item{max_repetitions}{Integer}
}
\value{
isotonic scale parameter estimate and deviance
}
\description{
Estimation of GPD parameters with constant shape parameter and non-decreasing
scale parameter
}

