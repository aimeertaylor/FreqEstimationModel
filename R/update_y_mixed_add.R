#'@title Update mixed data samaple clone addition
#'@name update_y_mixed_add
#'@description Function for addition of clone to a mixed data sample
#'@export
update_y_mixed_add <- function(to_update,
                               genotype_counts_mixed_add,
                               frequency_truncated)
  {

    genotype_counts_proposed_mixed_add <- genotype_counts_mixed_add # Allocate space
    x2 <- unique(to_update) # Instead of looping over all one-by-one, group into classes

    for (i in 1:length(x2)) {

      x3 <- x2[i] # Extract unique data type
      x4 <- to_update[to_update == x3] # Find all bloodsamples that have the same observed data as the first unique data type
      x5 <- names(x4) # Extract the names of those bloodsamples
      x6 <- t(rmultinom(n = length(x5), size = 1, prob = 1 * frequency_truncated[x3,])) # For each bloodsample pick one proposed genotype from those compatible
      genotype_counts_proposed_mixed_add[x5,] <- genotype_counts_mixed_add[x5,] + x6 # Update genotype_counts_proposed_mixed_add

    }

    return(genotype_counts_proposed_mixed_add)
  }
