#'@title Calculate legal genotype count moves
#'@name calculate_genotype_counts_truncated
#'@description Function to calculate the genotype counts that can be removed without rendering the ensuing likelihood equal to zero
#'@export
calculate_genotype_counts_truncated <- function(genotype_counts,
                                                comp_genotypes,
                                                raw_data,
                                                moi_prior_min2)
{

  # If there are any samples for which MOI prior minimum is 2 due to flanking info,
  # find those that are on the boundary so that we can overwrite below
  if(!is.null(moi_prior_min2)){
    moi_prior_min2_count <- rowSums(genotype_counts[moi_prior_min2,])
    moi_prior_min2_count2 <- names(moi_prior_min2_count[moi_prior_min2_count == 2])
  }

  moi_problem <- rowSums(genotype_counts) - 1

  # To distingush zero counts from wild type alleles down stream, store a version of the counts where zero counts are labled missing (NA)
  genotype_counts_NA <- genotype_counts
  genotype_counts_NA[genotype_counts_NA == 0] <- NA

  # 2) Instead of doing the matrix calculation genotype_counts%*%comp_genotypes, which is what I would normally do to determine the observed data
  # Break it down, s.t. the colSums of genotype_counts1 are equally to the number of mutant markers contributing to the observed data at SNP1 etc.
  # That way, if only there is only one mutant, we don't need to unravel the matrix multiplication to work which genotype count it came from.
  genotype_counts_temp <- genotype_counts # Allocate memory
  genotype_counts_temp[] <- 1 # Template of zero and ones, to truncate any counts in genotype_counts that are non-removable

  for(i in 1:ncol(comp_genotypes)) # Repeat for each marker (maximum number 5)
  {
    not_missing <- raw_data[rownames(genotype_counts),i,drop=FALSE]!= 99 # logical vector, defining which bloodsamples are unflexible at marker 1
    genotype_counts1 <- t(comp_genotypes[,i]*t(genotype_counts_NA[not_missing,,drop=FALSE])) # Double transpose (compared with diag, sweep, apply) appears to be the fastest option
    genotype_counts1_lower <- genotype_counts1*(rowSums(genotype_counts1, na.rm = TRUE) == 1) # Find out which contribute only one mutant to SNP 1, and mutiply by itself, s.t only keep record of those that cannot be removed
    genotype_counts1_upper <- genotype_counts1[(rowSums(genotype_counts1,na.rm=TRUE)==(moi_problem[not_missing])),,drop=FALSE] # Find out which contribute MOI-1 mutants to SNP 1, and throw away the rest.
    # Fill in the NAs, so we can update the template.
    genotype_counts1_lower[is.na(genotype_counts1_lower)]<- 0 # In the case of a single mutant, we're looking for a solitary one, so can make NAs equal to zero
    genotype_counts1_upper[is.na(genotype_counts1_upper)]<- 9.5 # In the case of MOI-1 mutants, we're looking for the solitary zero, so we set the NAs to 9.5. We make it a decimal number to hightlight the fact it isn't a mutant count, all of which are discrete.
    genotype_counts_temp[not_missing,][genotype_counts1_lower == 1] <- 0  # Make any solitary wild or mutant types = zero in template s.t counts are removed in a_current
    x <- genotype_counts_temp[not_missing, , drop = FALSE][rownames(genotype_counts1_upper),,drop=FALSE]
    x[genotype_counts1_upper == 0] <- 0
    genotype_counts_temp[rownames(x),] <- x
  }

  genotype_counts_truncated <- genotype_counts * genotype_counts_temp

  # Overwrite the samples for which MOI prior minimum is 2 due to flanking info (some overlap since some already have zero if mixed)
  if(!is.null(moi_prior_min2)){
    genotype_counts_truncated[moi_prior_min2_count2,] <- 0 # zero to choose from
  }

  # Catch mois of 1 when the data are missing
  genotype_counts_truncated[rowSums(genotype_counts) == 1,] <- 0

  return(genotype_counts_truncated)
}





