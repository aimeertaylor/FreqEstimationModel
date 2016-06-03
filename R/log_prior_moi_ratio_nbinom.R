#'@title Log prior MOI ratio
#'@name log_prior_moi_ratio_nbinom
#'@description Function to calculate the log of the prior MOI ratio when the MOI prior distribution is negative binomial.
#' Note that the normalising constants due to truncation cancel.
#'@export
log_prior_moi_ratio_nbinom <- function(moi_proposed,
                                       moi_current,
                                       moi_hyperparameter,
                                       moi_size_hyperparameter){

  log_prior_moi_ratio <- dnbinom(moi_proposed, mu = moi_hyperparameter, size = moi_size_hyperparameter, log = TRUE) -
                         dnbinom(moi_current, mu = moi_hyperparameter, size = moi_size_hyperparameter, log = TRUE)

  return(log_prior_moi_ratio)
}
