#'@title Log moi prior
#'@name log_moi_prior_calculation_unif
#'@description Log moi prior calculation
#'@export
log_moi_prior_calculation_unif <- function(moi_store,
                                           moi_list,
                                           y_no_mxed,
                                           y_mxed){


  # Log pmf depending
  log_moi_prior <- rep(c(-log(moi_list$moi_max),-log(moi_list$moi_max-1)), c(length(y_no_mxed), length(y_mxed)))

  # Store in a matrix
  log_moi_prior_matrix <- matrix(log_moi_prior,
                                 nrow = nrow(moi_store),
                                 ncol = ncol(moi_store),
                                 byrow = TRUE)

  # Re-order
  y_order <- dimnames(moi_store)[[2]]
  colnames(log_moi_prior_matrix) <- c(y_no_mxed,y_mxed)
  log_moi_prior_store <- log_moi_prior_matrix[,y_order]

  # Return
  return(log_moi_prior_store)

}
