# FUNCTION to calculate the likelihood of all the individuals in the data
# If the likelihood is impossible for a given individual this function
# should return -Inf for that individual
# This function could propbably be made faster
#
# Arguments:
# yij - data matrix
# zij - reference read matrix
# ai - the haplotype counts
# H_matrix - matrix of haplotypes
#
# Returns:
# log_likelihood - vector with the log likelihood for each indivdual
#------------------------------------------------------
#'@name log_likelihood_NGS
#'@title Log likelihood, NGS data
#'@description Calculate log likelihood of NGS data
#'@export
log_likelihood_NGS <- function(processed_data_list,
                               genotype_counts) {

  log_likelihood <- rep(NA, processed_data_list$no_total) # Allocate memory

  for(individual in 1:processed_data_list$no_total){

    # # Note that when pis <- (genotype_counts[individual,]/sum(genotype_counts[individual,])) %*% processed_data_list$comp_genotypes,
    # Errors produced because 1 not recognised as 1
    pis <- (genotype_counts[individual,] %*% processed_data_list$comp_genotypes)/sum(genotype_counts[individual,])

    log_likelihood[individual] <- sum(dbinom(processed_data_list$yij[individual,],
                                             size = processed_data_list$zij[individual,],
                                             prob = pis,
                                             log = TRUE))

  }

  return(log_likelihood)

}
