#'@title Log Multinomial Density
#'@name my_log_dmultinom
#'@description Function to calculate the density of multiple realisations of the multinomial distribution
#'simultaneously given either a single probability or a matrix of probability
#'@export
my_log_dmultinom <- function(x,
                             prob,
                             pop_fequency){

  if(pop_fequency == TRUE){

    ans <- drop(log(factorial(rowSums(x))) - rowSums(log(factorial(x))) + x %*% matrix(log(prob), ncol=1))

  } else {

    ans <- log(factorial(rowSums(x))) - rowSums(log(factorial(x))) + diag(x%*%log(t(prob)))

  }

  return(ans)
}



