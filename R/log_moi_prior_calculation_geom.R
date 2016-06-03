#'@title Log moi prior
#'@name log_moi_prior_calculation_geom
#'@description Log moi prior calculation
#'@export
log_moi_prior_calculation_geom <- function(moi_store,
                                           moi_list,
                                           y_no_mxed,
                                           y_mxed){

    # Split data into y_mxed (MOI = 1 not possible) and y_non_mxed (MOI = 1 possible)
    normalising_geom_no_mxed <- sum(dgeom(1:moi_list$moi_max, 1/moi_list$moi_hyperparameter,log = FALSE))
    normalising_geom_mxed <- sum(dgeom(2:moi_list$moi_max, 1/moi_list$moi_hyperparameter,log = FALSE))

    log_normalising_geom_no_mxed_mxed <- log(rep(c(normalising_geom_no_mxed, normalising_geom_mxed), c(length(y_no_mxed), length(y_mxed))))

    # Log pmf
    log_moi_prior <- t(t(dgeom(moi_store[,c(y_no_mxed,y_mxed)], 1/moi_list$moi_hyperparameter, log = TRUE)) - log_normalising_geom_no_mxed_mxed)

    # Re-order
    y_order <- dimnames(moi_store)[[2]]
    colnames(log_moi_prior) <- c(y_no_mxed,y_mxed)
    log_moi_prior_store <- log_moi_prior[,y_order]

    # Return
    return(log_moi_prior_store)

}
