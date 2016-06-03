#'@title Draw initial MOIs
#'@name draw_moi_initial
#'@description Function for drawing initial MOIs.
#'@param moi_list$moi_prior is a character, either \code{'Poisson'}, \code{'nBinomial'}, \code{'Geometric'} or \code{'Uniform'}, that defines the prior distribution over the MOI.
#'@param datasampleID is a vector of sample IDs; that is, a character vector such as \code{c('id1:1','id2:1','id3:1','id4:1')}.
#'It is required to ensure the order of the initial MOI estimates is consistent with the data.
#'@param y_no_mxed is a vector comprising a subset of \code{datasampleID}; for example \code{c('id1:1','id4:1')}.
#'The subset includes sample IDs for all samples for which there is no evidence of multiclonality; that is, samples for which an MOI of one is possible.
#'@param y_mxed is a vector comprising a subset of \code{datasampleID}; for example \code{c('id1:2','id4:3')}.
#'The subset includes sample IDs for all samples for which there is evidence of multiclonality; that is, samples for which an MOI of one is not possible.
#'Note that, \code{y_no_mxed} and \code{y_mxed} are disjoints sets that together union of is equal to the unordered equivalent of \code{datasampleID}.
#'@param moi_list$moi_max is an integer defining the maximum moi possible in the model
#'@param moi_list$moi_hyperparameter is a parameter of the prior distribution
#'@param moi_list$moi_size_hyperparameter is a parameter of the prior distribution
#'@export
draw_moi_initial <- function(moi_list,
                             datasampleID,
                             y_no_mxed,
                             y_mxed){


  # Initialise MOIs by random draws
  m_initial <- vector('numeric', length = length(datasampleID))
  names(m_initial) <- datasampleID

  if(moi_list$moi_prior=='Uniform'){
    m_initial[y_no_mxed] <- sample(1:moi_list$moi_max, length(y_no_mxed), replace = TRUE)
    m_initial[y_mxed] <- sample(2:moi_list$moi_max, length(y_mxed), replace = TRUE)
  }

  if(moi_list$moi_prior=='Poisson'){
    m_initial[y_no_mxed]<- sample(1:moi_list$moi_max,length(y_no_mxed), replace = TRUE,
                                  prob = dpois(1:moi_list$moi_max, moi_list$moi_hyperparameter)/sum(dpois(1:moi_list$moi_max, moi_list$moi_hyperparameter)))
    m_initial[y_mxed] <- sample(2:moi_list$moi_max,length(y_mxed), replace=TRUE,
                                prob=dpois(2:moi_list$moi_max,moi_list$moi_hyperparameter)/sum(dpois(2:moi_list$moi_max,moi_list$moi_hyperparameter)))
  }

  if(moi_list$moi_prior=='Geometric'){
    m_initial[y_no_mxed]<-sample(1:moi_list$moi_max,length(y_no_mxed), replace=TRUE,
                                 prob=dgeom(1:moi_list$moi_max,1/moi_list$moi_hyperparameter)/sum(dgeom(1:moi_list$moi_max,1/moi_list$moi_hyperparameter)))
    m_initial[y_mxed]<-sample(2:moi_list$moi_max,length(y_mxed), replace=TRUE,
                              prob=dgeom(2:moi_list$moi_max,1/moi_list$moi_hyperparameter)/sum(dgeom(2:moi_list$moi_max,1/moi_list$moi_hyperparameter)))
  }

  if(moi_list$moi_prior=='nBinomial'){
    m_initial[y_no_mxed]<-sample(1:moi_list$moi_max,length(y_no_mxed), replace=TRUE,
                                 prob = dnbinom(1:moi_list$moi_max, mu = moi_list$moi_hyperparameter, size = moi_list$moi_size_hyperparameter)/sum(dnbinom(1:moi_list$moi_max,mu=moi_list$moi_hyperparameter,size=moi_list$moi_size_hyperparameter)))
    m_initial[y_mxed]<-sample(2:moi_list$moi_max,length(y_mxed),replace=TRUE,
                              prob = dnbinom(2:moi_list$moi_max,mu=moi_list$moi_hyperparameter, size = moi_list$moi_size_hyperparameter)/sum(dnbinom(2:moi_list$moi_max,mu=moi_list$moi_hyperparameter,size=moi_list$moi_size_hyperparameter)))
  }

  return(m_initial)
}



