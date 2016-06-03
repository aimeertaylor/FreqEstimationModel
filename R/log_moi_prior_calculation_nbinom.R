#'@title Log moi prior
#'@name log_moi_prior_calculation_nbinom
#'@description Log moi prior calculation
#'@export
log_moi_prior_calculation_nbinom <- function(moi_store,
                                             moi_list,
                                             y_no_mxed,
                                             y_mxed){



    normalising_nbinom_no_mxed <- sum(dnbinom(1:moi_list$moi_max, mu = moi_list$moi_hyperparameter, size = moi_list$moi_size_hyperparameter, log = FALSE))
    normalising_nbinom_mxed <- sum(dnbinom(2:moi_list$moi_max, mu = moi_list$moi_hyperparameter, size = moi_list$moi_size_hyperparameter, log = FALSE))

    # Split data into y_mxed (MOI = 1 not possible) and y_non_mxed (MOI = 1 possible)
    log_normalising_nbinom_no_mxed_mxed <- log(rep(c(normalising_nbinom_no_mxed,normalising_nbinom_mxed),c(length(y_no_mxed),length(y_mxed))))

    # Log pmf
    log_moi_prior <- t(t(dnbinom(moi_store[,c(y_no_mxed,y_mxed)], mu = moi_list$moi_hyperparameter, size = moi_list$moi_size_hyperparameter, log = TRUE))
                       - log_normalising_nbinom_no_mxed_mxed)

    # Re-order
    y_order <- dimnames(moi_store)[[2]]
    colnames(log_moi_prior) <- c(y_no_mxed,y_mxed)
    log_moi_prior_store <- log_moi_prior[,y_order]

    # Return
    return(log_moi_prior_store)

}








