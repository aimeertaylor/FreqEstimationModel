#'@title Update mixed data samaple clone removal
#'@name update_y_mixed_sub
#'@description Function for removal of clone from mixed data sample
#'@export
######################################################################################
# Function for y_mxed m*= m-1
# We worried about removing a clone that would render the genotype count not compatible with the data
# y = 0 -> 0.5 not a problem because cannot happen when we're only allowing clone removal
# y = 1 -> 0.5 not a problem because cannot happen when we're only allowing clone removal
# y = 0.5 -> 0 is a problem if only one mutant clone contributes to the mixed nature of the SNP
# y =  0.5 -> 1 is a problem if MOI-1 mutant clones contribute to the mutant SNP
# We are, therefore, interested in finding the counts that either
# contribute a single mutant OR a single wild type
######################################################################################
update_y_mixed_sub <- function(to_update,
                             genotype_counts_mixed_sub,
                             genotype_counts_truncated)
{
  genotype_counts_proposed_mixed_sub <- genotype_counts_mixed_sub # Allocate space

  # Find patients where genotype count update is deterministic, since only one genotype to choose from (may be faster to incorporate into rmultinom)
  x1 <- names(which(rowSums(genotype_counts_truncated[to_update,,drop=FALSE]!=0) == 1))# Genotype counts of patients that have only one type of removable genotype

  # Update
  if(length(x1) > 0){
    genotype_counts_proposed_mixed_sub[x1,] <- genotype_counts_mixed_sub[x1,,drop=FALSE]-(genotype_counts_truncated[x1,,drop=FALSE]!=0)
  }

  # Update remaining patients for whom there is a choice as to which genotype to remove
  x3 <- names(which(rowSums(genotype_counts_truncated[to_update,,drop=FALSE]!=0) > 1)) # Names of non-deterministic update

  if(length(x3) > 0){
    x4 <- apply(genotype_counts_truncated[x3,,drop=FALSE],1, paste, sep='', collapse='')
    x5 <- unique(x4) # Instead of pick each one to remove one-by-one (test speed compared with one-by-one approach), group together classes that are the same
    for(i in 1:length(x5)){
      x6 <- names(which(x4==x5[i]))
      x7 <- t(rmultinom(n = length(x6), size=1, prob = genotype_counts_truncated[x6[1],]));

      genotype_counts_proposed_mixed_sub[x6,] <- genotype_counts_proposed_mixed_sub[x6,] - x7
    }
  }

  return(genotype_counts_proposed_mixed_sub)
}
