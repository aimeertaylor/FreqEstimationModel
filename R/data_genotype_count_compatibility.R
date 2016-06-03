#'@title Data genotype count compatibility
#'@description Data and genotype compatibility check function
#'@name data_genotype_count_compatibility
#'@export
data_genotype_count_compatibility <- function(comp_genotypes,
                                              genotype_counts,
                                              raw_data)
{
  no_markers <- ncol(raw_data)
  y_unnormalised <- genotype_counts %*% comp_genotypes
  y <- y_unnormalised/rowSums(genotype_counts) # Normalise by MOI
  y[!y %in% 0:1] <- 0.5 # Set mixed to 0.5
  y[raw_data == 99] <- 99 # Set missing to 99
  non_matching <- names(which((rowSums(y==raw_data)!=no_markers)==TRUE))
  return(non_matching)
}
