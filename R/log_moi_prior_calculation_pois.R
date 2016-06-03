#'@title Log moi prior
#'@name log_moi_prior_calculation_pois
#'@description Log moi prior calculation
#'@export
log_moi_prior_calculation_pois <- function(moi_store,
                                           moi_list,
                                           y_no_mxed,
                                           y_mxed){

    # Calculate the normalising constant dependent on whether demonstrably mixed or not
    normalising_pois_no_mxed <- sum(dpois(c(1:moi_list$moi_max), lambda = moi_list$moi_hyperparameter))
    normalising_pois_mxed <- sum(dpois(c(2:moi_list$moi_max), lambda = moi_list$moi_hyperparameter))

    # Concatenate and take logs
    log_normalising_pois_no_mxed_mxed <- log(rep(c(normalising_pois_no_mxed, normalising_pois_mxed), c(length(y_no_mxed),length(y_mxed))))

    # Log pmf
    log_moi_prior <- t(t(dpois(moi_store[,c(y_no_mxed,y_mxed)], moi_list$moi_hyperparameter, log = TRUE))
                       - log_normalising_pois_no_mxed_mxed)

    # Re-order
    y_order <- dimnames(moi_store)[[2]]
    colnames(log_moi_prior) <- c(y_no_mxed,y_mxed)
    log_moi_prior_store <- log_moi_prior[,y_order]

    # Return
    return(log_moi_prior_store)
}
