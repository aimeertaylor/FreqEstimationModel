% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_prior_moi_ratio_nbinom.R
\name{log_prior_moi_ratio_nbinom}
\alias{log_prior_moi_ratio_nbinom}
\title{Log prior MOI ratio}
\usage{
log_prior_moi_ratio_nbinom(
  moi_proposed,
  moi_current,
  moi_hyperparameter,
  moi_size_hyperparameter
)
}
\description{
Function to calculate the log of the prior MOI ratio when the MOI prior distribution is negative binomial.
Note that the normalising constants due to truncation cancel.
}
