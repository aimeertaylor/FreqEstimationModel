#'@title Random draw from the dirichlet
#'@name rdirichlet
#'@description Same as the rdirichlet function in MCMCpack but rewritten to minimise dependencies.
#'@details Need to add log option
#'@export
rdirichlet <- function(n, alpha)
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}
