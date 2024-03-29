% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_moi_initial.R
\name{draw_moi_initial}
\alias{draw_moi_initial}
\title{Draw initial MOIs}
\usage{
draw_moi_initial(moi_list, datasampleID, y_no_mxed, y_mxed)
}
\arguments{
\item{datasampleID}{is a vector of sample IDs; that is, a character vector such as \code{c('id1:1','id2:1','id3:1','id4:1')}.
It is required to ensure the order of the initial MOI estimates is consistent with the data.}

\item{y_no_mxed}{is a vector comprising a subset of \code{datasampleID}; for example \code{c('id1:1','id4:1')}.
The subset includes sample IDs for all samples for which there is no evidence of multiclonality; that is, samples for which an MOI of one is possible.}

\item{y_mxed}{is a vector comprising a subset of \code{datasampleID}; for example \code{c('id1:2','id4:3')}.
The subset includes sample IDs for all samples for which there is evidence of multiclonality; that is, samples for which an MOI of one is not possible.
Note that, \code{y_no_mxed} and \code{y_mxed} are disjoints sets that together union of is equal to the unordered equivalent of \code{datasampleID}.}

\item{moi_list$moi_prior}{is a character, either \code{'Poisson'}, \code{'nBinomial'}, \code{'Geometric'} or \code{'Uniform'}, that defines the prior distribution over the MOI.}

\item{moi_list$moi_max}{is an integer defining the maximum moi possible in the model}

\item{moi_list$moi_hyperparameter}{is a parameter of the prior distribution}

\item{moi_list$moi_size_hyperparameter}{is a parameter of the prior distribution}
}
\description{
Function for drawing initial MOIs.
}
