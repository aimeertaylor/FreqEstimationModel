#'@title Compatible genotypes
#'@name compatible_genotypes
#'@description Function to ascertain the genotypes that are compatible with a given observation
#'@export
compatible_genotypes <- function(observation,
                                 as_numeric = TRUE)
{
  L <- length(observation) # Number of loci
  yl = vector('list',length = L) # yl is a list with entries of 0, 1, or c(0, 1)

  for (j in 1:L){    # convert data to 0s and 1s
    if (observation[j] == .5) {
      yl[[j]] = 0:1
    }   else    {
      yl[[j]] = observation[j]
    }
  }

  hapmat = as.matrix(expand.grid(yl)) # hapmat uses yl to construct matrix of all potential genotypes
  genotypes <- apply(hapmat,1,paste,sep = '',collapse = '')

  if(as_numeric == TRUE){

    return(hapmat)

  }else{

    return(genotypes)
  }
}
