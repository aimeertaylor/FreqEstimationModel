#'@title Propose MOI
#'@name draw_moi_proposed
#'@description Function for drawing new MOI vector. Note that there is no option to stay still, and that Proposal density does not have equal probability over state space [m_min,..,m_max] and is not symmetric: Pr(m*|m)=Pr(m|m*)
#'@export
draw_moi_proposed <- function(moi_current,
                              moi_current_equal_moi_min,
                              moi_max){

    # Probabilistic update (do for all then overwrite at boundaries)
    moi_proposed <- moi_current + sample(c(-1,1), length(moi_current), replace = TRUE)

    # Overwrite by deterministic update (includes moi_min = {moi_i^t}_min, 2 and 1)
    moi_proposed[moi_current_equal_moi_min] <- moi_current[moi_current_equal_moi_min] + 1
    moi_proposed[moi_current == moi_max] <- moi_current[moi_current == moi_max] - 1

    return(moi_proposed)
  }
# #-----------------------------------------------------------------------------------------------------------------------
# # probability over state space [m_min,..,m_max] is not symmetric
# I<-1000; moi_min<-2; moi_max<-8
# moistore<-c()
# for(i in 1:I){moi<-draw_moi_proposed(moi,moi_min,moi_max); moistore<-c(moistore,moi)}
# barplot(table(moistore))
# #-----------------------------------------------------------------------------------------------------------------------
