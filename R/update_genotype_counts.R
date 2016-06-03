###################################################################################################################
# moi_diff                   - difference between new MOI and current MOI
# y_list <- list(
# processed_data_list$masked_data_character                - data, but in character form
# processed_data_list$y_mxed_andor_missing                    - names of discernibly multiclonal data (wild and mutant type markers have both been detected at at least one SNP)
# processed_data_list$y_pure_nomissing                     - names of non-discernible multiclonal data
# genotype_counts_current    - current genotype counts
# processed_data_list$frequency_truncated             - a list of genotypes compatible with each class of observed data (used to truncate frequency)
# Transpose fastest option for multiplying rows of a matrix by a vector according to http://stackoverflow.com/questions/3643555/multiply-rows-of-matrix-by-vector
###################################################################################################################
#'@title Update genotype counts
#'@name update_genotype_counts
#'@description Function to update genotype counts conditional on current genotype counts and new multiplicity of infection (MOI)
#'@export
update_genotype_counts <- function(moi_diff,
                                   genotype_counts_current,
                                   current_genotype_counts_truncated,
                                   processed_data_list)
{

  genotype_counts_proposed <- genotype_counts_current # allocate space for the new genotype counts
  deterministic_update <- moi_diff[processed_data_list$y_pure_nomissing]
  probablistic_update <- moi_diff[processed_data_list$y_mxed_andor_missing]


  #--------------------------------------------------------------------------------------------
  # 1) Addition processed_data_list$y_pure_nomissing
  x1 <- deterministic_update[deterministic_update == 1] # Extract MOI change + 1 for processed_data_list$y_pure_nomissing
  if(length(x1) > 0){
    x2 <- names(x1)
    x3 <- genotype_counts_current[x2,,drop=FALSE] # Extract current genotype counts for processed_data_list$y_pure_nomissing with m*=m+1
    genotype_counts_proposed[x2,]<- x3 + (x3 != 0) # Update genotype counts | m* (add one to non-zero count)
  }
  #--------------------------------------------------------------------------------------------


  #--------------------------------------------------------------------------------------------
  # 2) Addition processed_data_list$y_mxed_andor_missing
  x4 <- processed_data_list$masked_data_character[names(which(moi_diff[processed_data_list$y_mxed_andor_missing] == 1))] # Extract data, in character form, relating to processed_data_list$y_mxed_andor_missing m*=m+1
  if(length(x4) > 0){
    x5 <- names(x4)
    genotype_counts_proposed[x5,] <- update_y_mixed_add(to_update = x4,
                                                        genotype_counts_mixed_add = genotype_counts_current[x5,,drop=FALSE],
                                                        processed_data_list$frequency_truncated)
  }
  #--------------------------------------------------------------------------------------------


  #--------------------------------------------------------------------------------------------
  # 4) Subtraction processed_data_list$y_pure_nomissing
  x6 <- deterministic_update[deterministic_update == -1] # Extract MOI change - 1 for processed_data_list$y_pure_nomissing
  if(length(x6 > 0)){
    x7 <- names(x6) # Extract bloodsample id for processed_data_list$y_pure_nomissing m*=m-1
    x8 <- genotype_counts_current[x7,] # Extract current genotype counts for processed_data_list$y_pure_nomissing with m*=m+{-1,+1}
    genotype_counts_proposed[x7,] <- x8 - (x8 != 0) # Update genotype counts | m*
  }
  #--------------------------------------------------------------------------------------------


  #--------------------------------------------------------------------------------------------
  # 5) Subtraction processed_data_list$y_mxed_andor_missing
  x9 <- names(which(moi_diff[processed_data_list$y_mxed_andor_missing] == -1)) # Names of problem genotype counts
  if(length(x9 > 0)){
    genotype_counts_proposed[x9,] <- update_y_mixed_sub(to_update = x9,
                                                        genotype_counts_mixed_sub = genotype_counts_current[x9,,drop=FALSE],
                                                        current_genotype_counts_truncated)
  }
  #--------------------------------------------------------------------------------------------


  return(genotype_counts_proposed)

}

