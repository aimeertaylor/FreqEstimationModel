#'@title Calculate log MOI proposal ratio
#'@name calculate_log_moi_proposal_ratio
#'@description Function to calculate the log MOI proposal ratio
#'@export
calculate_log_moi_proposal_ratio <- function(moi_max,
                                             moi_current,
                                             moi_proposed,
                                             moi_current_equal_moi_min,
                                             moi_proposed_equal_moi_min,
                                             datasampleID,
                                             no_total){

  # 0 for all log(p(m^t)/p(m*)) not on the boundry:
  log_moi_proposal_ratio <- array(0, dim = no_total, dimnames = list(datasampleID))
  # m^t = moi_max:
  log_moi_proposal_ratio[moi_current == moi_max] <- -log(2)
  # m* = moi_max:
  log_moi_proposal_ratio[moi_proposed == moi_max] <- log(2)
  # m^t = m_min, where m_min includes {m_i}_min^t, 2 and 1
  log_moi_proposal_ratio[moi_current_equal_moi_min] <- -log(2)
  # m* = m_min, where m_min includes {m_i}_min^t, 2 and 1
  log_moi_proposal_ratio[moi_proposed_equal_moi_min] <- log(2)
  return(log_moi_proposal_ratio)
}


