#'@title lLg prior MOI ratio
#'@name log_prior_moi_ratio_geom
#'@description Function to calculate the log of the prior MOI ratio when the MOI prior distribution is geometric.
#' Note that the normalising constants due to truncation cancel.
#'@export
log_prior_moi_ratio_geom <- function(moi_proposed,
                                     moi_current,
                                     moi_hyperparameter,
                                     moi_size_hyperparameter){

  log_prior_moi_ratio <- dgeom(x = moi_proposed, prob = 1/moi_hyperparameter, log = TRUE) -
                         dgeom(x = moi_current, prob = 1/moi_hyperparameter, log = TRUE)

  return(log_prior_moi_ratio)
}


