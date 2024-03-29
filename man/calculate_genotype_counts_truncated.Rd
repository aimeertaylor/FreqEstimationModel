% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_genotype_counts_truncated.R
\name{calculate_genotype_counts_truncated}
\alias{calculate_genotype_counts_truncated}
\title{Calculate legal genotype count moves}
\usage{
calculate_genotype_counts_truncated(
  genotype_counts,
  comp_genotypes,
  raw_data,
  moi_prior_min2
)
}
\description{
Function to calculate the genotype counts that can be removed without rendering the ensuing likelihood equal to zero
}
