#'@title Preprocess data
#'@name preprocess_data
#'@description Function to preprocess data
#'@export
preprocess_data <- function(data_summary,
                            log_like_zero,
                            NGS,
                            augment_missing_data,
                            moi_prior_min2){


  # Data preprocessing
  # Number of marker loci in the data
  nlocus <- length(grep('locus', colnames(data_summary$Data)))

  if(nlocus > 0){

    no_marker_loci <- nlocus

  } else {

    no_marker_loci <- ncol(data_summary$Data)

  }

  # Overwrite data if log_like_zero is TRUE
  if(log_like_zero){data_summary$Data[,1:no_marker_loci] <- 99}

  # Read counts if NGS
  if(NGS){

    yij <- data_summary$yij[, 1:no_marker_loci, drop=FALSE]
    zij <- data_summary$zij[, 1:no_marker_loci, drop=FALSE]

  } else {

    yij <- NA
    zij <- NA

  }

  raw_data_pre_augment_step <- data_summary$Data[, 1:no_marker_loci, drop = FALSE] # Extract data into matrix with row and col names.


  # Missing data (logical vector)
  ind_partial_initial <- apply(raw_data_pre_augment_step, 1, function(X){99 %in% X}) # Logical, with 'frequency' names. Locates datasamples with missing data.

  if(!augment_missing_data & sum(ind_partial_initial) > 0){

    raw_data <- raw_data_pre_augment_step[!ind_partial_initial,, drop = FALSE]# Throw away partially observed data only if there are any

  } else {

    raw_data <- raw_data_pre_augment_step

  }

  # IDs
  datasampleID <- rownames(raw_data)
  markerID <- colnames(raw_data)

  # Tag blood samples as:
  # 1) partial; 2) mixed; 3) pure, no mixed or missing (where partial and mixed are not mutually exclusive)
  ind_partial <- apply(raw_data, 1, function(X){99 %in% X}) # If datasamples with missing data are retained, separate datasamples with missing data
  ind_mixed <- apply(raw_data, 1, function(X){0.5 %in% X}) # Separate datasamples with discernibly multiclonal infections
  ind_non_mixed <- as.logical((!ind_partial)*(!ind_mixed)) # datasamples that are pure and have no missing

  # Using the tags above, separate into four types of blood sample:
  if(is.null(moi_prior_min2)){

    y_no_mxed <- datasampleID[!ind_mixed] # No mixed (inc. missing). Separate for draw_m_new (MOI=1 possible).
    y_mxed <- datasampleID[ind_mixed] # Least one mixed (inc. missing). Separate for draw_m_new (MOI=1 not possible).

    } else {

    y_mxed <- unique(c(datasampleID[ind_mixed],moi_prior_min2)) # mixed plus extra
    y_no_mxed <- datasampleID[match(datasampleID, y_mxed,  nomatch = 0L) == 0] # The remaining

  }
  y_pure_nomissing <- datasampleID[ind_non_mixed] # Neither mixed not missing. Separate for deterministic update_genotype_counts.
  y_mxed_andor_missing <- datasampleID[!ind_non_mixed] # Both mixed and or missing. Separate for probabilistic update_genotype_counts.



  # Numbers of data samples
  no_total <- nrow(raw_data) # Total number of blood samples to analyse

  # In case no_partial > 0, manipulate data before trying to ascertain the number of compatible genotypes and initial genotype counts
  masked_data <- raw_data
  masked_data[masked_data == 99] <- 0.5 # Because missing data (99s) are compatible with wild (0) or mutant (1) type markers, change all 99s to 0.5s before ascertaining compatible genotypes.
  masked_data_unique <- unique(masked_data) # To ascertain compatible genotypes
  masked_data_character <- apply(masked_data, 1, paste, sep='', collapse=',') # To access frequency_truncated within update_genotype_counts() and to initalise genotype counts
  masked_data_character_unique <- unique(masked_data_character) # Character vector of unique data types, to allocate frequency_truncated row names
  no_masked_data_character_unique <- length(masked_data_character_unique)

  # Ascertain compatible genotypes
  comp_genotypes <- c() # How many genotypes are compatible depends on the data, we cannot, therefore, preallocate matrix size.

  for(i in 1:no_masked_data_character_unique){# At most, length of for loop is 3^no.makers
    x <- compatible_genotypes(masked_data_unique[i,])
    comp_genotypes <- rbind(comp_genotypes,x)
  }

  comp_genotypes <- unique(comp_genotypes) # remove duplicates
  comp_genotypes_character <- apply(comp_genotypes, 1, paste, collapse='')
  dimnames(comp_genotypes) <- list(comp_genotypes_character, markerID)
  no_haplotypes <- nrow(comp_genotypes) # The number of genotypes compatible with the data, and hence number of genotype frequencies

  # Create frequency_truncated and neg_log_sum_frequency_truncated (only needs doing once since static conditional on data)
  # To truncate multinomial prob when adding a clone, work out which genotypes are compatible with each type of manipulated data (frequency_truncated)
  frequency_truncated <- matrix(0,
                                nrow = no_masked_data_character_unique,
                                ncol = no_haplotypes,
                                dimnames = list(masked_data_character_unique,
                                                comp_genotypes_character))

  for(i in 1:no_masked_data_character_unique){

    frequency_truncated[i, compatible_genotypes(unique(masked_data)[i,],as_numeric=FALSE)] <- 1

  }

  neg_log_sum_frequency_truncated <- -log(rowSums(frequency_truncated)) # This is used for calculating the proposal probability


  processed_data_list <- list(no_marker_loci = no_marker_loci,
                              yij = yij,
                              zij = zij,
                              raw_data = raw_data,
                              datasampleID = datasampleID,
                              markerID = markerID,
                              y_no_mxed = y_no_mxed,
                              y_mxed = y_mxed,
                              y_pure_nomissing = y_pure_nomissing,
                              y_mxed_andor_missing = y_mxed_andor_missing,
                              no_total = no_total,
                              masked_data_character = masked_data_character,
                              comp_genotypes = comp_genotypes,
                              comp_genotypes_character = comp_genotypes_character,
                              no_haplotypes = no_haplotypes,
                              frequency_truncated = frequency_truncated,
                              neg_log_sum_frequency_truncated = neg_log_sum_frequency_truncated)

  return(processed_data_list)
}

