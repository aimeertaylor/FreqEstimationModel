#'@title Draw frequency proposed child
#'@name draw_frequency_proposed_child
#'@description Function for drawing frequency proposed per child
#'@export
draw_frequency_proposed_child <- function(Alpha,
                                          genotype_counts_current){

  frequencies_proposed <- genotype_counts_current # Allocate memory
  frequencies_proposed[] <- NA # Frequencies_proposed is now a matrix

  children <- matrix(unlist(strsplit(rownames(frequencies_proposed), ':')), ncol = 2, byrow = TRUE)[,1]

  for(child in unique(children)){
    indexes <- which(child == children)
    Xtemp <- rdirichlet(n = 1, alpha = (colSums(genotype_counts_current[indexes,,drop = FALSE]) + Alpha))

    Xtemp[Xtemp == 0] <- .Machine$double.eps
    Xtemp[Xtemp == 1] <- 1-.Machine$double.eps
    TempVar2 <- Xtemp/sum(Xtemp)

    frequencies_proposed[indexes,]<-rep(TempVar2, each = length(indexes))
  }

  return(frequencies_proposed)
}
