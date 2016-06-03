#'@title Draw initial genotype counts
#'@name draw_genotype_counts_initial
#'@description Function for drawing initial genotype counts.
#'@export
draw_genotype_counts_initial <- function(moi_current,
                                         processed_data_list,
                                         frequency_current){

  genotype_counts_initial <- matrix(data = NA, # Allocate memory
                                    nrow = processed_data_list$no_total,
                                    ncol = processed_data_list$no_haplotypes,
                                    dimnames = list(processed_data_list$datasampleID,
                                                    processed_data_list$comp_genotypes_character)) # Matrix in which to store genotype counts

  # Draw intial genotype counts
  for(i in 1:processed_data_list$no_total){

    genotype_counts_initial[i,] <- as.vector(rmultinom(n = 1,
                                                       size = moi_current[i],
                                                       prob = frequency_current * processed_data_list$frequency_truncated[processed_data_list$masked_data_character[i],]))
  }

  non_matching <- data_genotype_count_compatibility(processed_data_list$comp_genotypes, genotype_counts_initial, processed_data_list$raw_data)# Check to see if inital A draw are compatible with the data

  # Assess data compatibility and re-draw if necessary
  count <- 1

  while(length(non_matching > 0)){

    if(count%%100 == 0){print(non_matching)}

    for(i in 1:length(non_matching)){

      genotype_counts_initial[non_matching[i],] <- as.vector(rmultinom(n = 1,
                                                                       size = moi_current[non_matching[i]],
                                                                       prob = frequency_current * processed_data_list$frequency_truncated[processed_data_list$masked_data_character[non_matching[i]],]))

    }

    non_matching <- data_genotype_count_compatibility(processed_data_list$comp_genotypes, genotype_counts_initial, processed_data_list$raw_data)

    count <- count + 1
  }

  return(genotype_counts_initial)
}
