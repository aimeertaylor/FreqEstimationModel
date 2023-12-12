#' @title Compatible genotypes
#' @name compatible_genotypes
#' @description Function to ascertain the genotypes that are compatible with a given observation
#'
#' @param observation A vector of allelic calls for a given individual, with one
#'   call per loci. Heteroallelic loci can either be encoded as 0.5.
#'   Alternatively, all alleles detected at heteroallelic loci can be listed.
#'   See examples below.
#'
#' @export
#'
#' @examples
#' # Example where heteroallelic loci are encoded as 0.5
#' y <- c(0, 1, 0.5, 0.5)
#' compatible_genotypes(y)
#' compatible_genotypes(y, FALSE)
#'
#' # Example where all alleles detected at heteroallelic loci are listed
#' y <- list(0, 1, c(0,1), c(0,1,2))
#' compatible_genotypes(y)
#' compatible_genotypes(y, FALSE)
compatible_genotypes <- function(observation,
                                 as_numeric = TRUE)
{

  L <- length(observation) # Number of loci

  if (!is.list(observation)) { # If yl is not a list
    if (!all(unique(observation) %in% c(0, 1, 0.5))) { # Check entries are in {0, 1, 0.5}
      stop("invalid allelic observations")
      } else { # If yl entries are in {0, 1, 0.5} convert to {0, 1, {0,1}}
        yl = vector('list',length = L)
        for (j in 1:L){    # convert data to 0s and 1s
          if (observation[j] == .5) {
            yl[[j]] = 0:1
          }   else    {
            yl[[j]] = observation[j]
          }
        }
      }
  } else {
    yl <- observation
  }

  hapmat = as.matrix(expand.grid(yl)) # hapmat uses yl to construct matrix of all potential genotypes
  genotypes <- apply(hapmat,1,paste,sep = '',collapse = '')

  if (as_numeric == TRUE){
    return(hapmat)
  } else {
    return(genotypes)
  }
}
