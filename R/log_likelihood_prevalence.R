#'@name log_likelihood_prevalence
#'@title Log likelihood, prevalence data
#'@description Calculate log likelihood of prevalence data
#'@export
log_likelihood_prevalence <- function(processed_data_list,
                                      genotype_counts){

    log_likelihood <- rep(0, processed_data_list$no_total)

    return(log_likelihood)
}
