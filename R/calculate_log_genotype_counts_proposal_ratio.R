#'@title Calculate the log genotype counts proposal ratio
#'@name calculate_log_genotype_counts_proposal_ratio
#'@description Function to calculate the log genotype counts proposal ratio
#'@export
calculate_log_genotype_counts_proposal_ratio <- function(moi_diff,
                                                         genotype_counts_diff,
                                                         neg_log_sum_frequency_truncated,
                                                         mh_genotype_counts_truncated,
                                                         y_character,
                                                         log_target_ratio,
                                                         datasampleID,
                                                         no_total)

{

  # Genotype count proposal
  log_genotype_counts_proposal_ratio <- array(0, dim = no_total, dimnames = list(datasampleID))
  # Calculate forward (defined in terms of m*=m+1)
  log_forwardadd_backwardremove_proposal <- neg_log_sum_frequency_truncated[y_character]
  # Vector constructed by con catonating col1, col2...
  How_many_to_choose_frm <- t(mh_genotype_counts_truncated)[abs(t(genotype_counts_diff)) == 1]
  # Calculate backward (defined in terms of m*=m+1)
  log_backwardadd_forwardremove_proposal <- log(How_many_to_choose_frm) - log(rowSums(mh_genotype_counts_truncated))
  # Calculate the log ratio defined in terms of m*=m+1)
  log_ratio <- log_backwardadd_forwardremove_proposal - log_forwardadd_backwardremove_proposal
  # If m*=m-1, forward and backward are reserved by moi_diff
  log_genotype_counts_proposal_ratio <- moi_diff * log_ratio
  return(log_genotype_counts_proposal_ratio)

}
