#'@title Log prior MOI ratio
#'@name log_prior_moi_ratio_pois
#'@description Function to calculate the log of the prior MOI ratio when the MOI prior distribution is poisson.
#' Note that the normalising constants due to truncation cancel.
#'@export
log_prior_moi_ratio_pois <- function(moi_proposed,
                                     moi_current,
                                     moi_hyperparameter,
                                     moi_size_hyperparameter){

  log_prior_moi_ratio <- dpois(x = moi_proposed, lambda = moi_hyperparameter, log = TRUE) - dpois(x = moi_current, lambda = moi_hyperparameter, log = TRUE)

  return(log_prior_moi_ratio)
}
