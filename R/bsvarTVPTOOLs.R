

#' @title Extracts the posterior draws of the structural matrix \eqn{B} and 
#' returns the as an \code{array}
#' @description Extracts the posterior draws of the structural matrix \eqn{B} 
#' and returns the as an \code{array}
#' @param posterior Posterior draws of from a Markov-switching model
#' @return An array containing the posterior draws of the structural matrix \eqn{B}
structural_to_array <- function(posterior) {
  
  B_posterior = posterior$posterior$B_cpp
  S       = dim(B_posterior)[1]
  N       = dim(B_posterior[1,1][[1]])[1]
  M       = dim(B_posterior[1,1][[1]])[3]
  
  B_out   = array(NA, c(N, N, M, S))
  for (s in 1:S) {
    B_out[,,,s]  = B_posterior[s,1][[1]]
  }
  
  return(B_out)
}



#' @title Extracts the posterior draws of the autoregressive matrix \eqn{A} and 
#' returns the as an \code{array}
#' @description Extracts the posterior draws of the autoregressive matrix \eqn{A} 
#' and returns the as an \code{array}
#' @param posterior Posterior draws of from a Markov-switching model
#' @return An array containing the posterior draws of the autoregressive matrix \eqn{A}
autoregressive_to_array <- function(posterior) {
  
  A_posterior = posterior$posterior$A_cpp
  S       = dim(A_posterior)[1]
  N       = dim(A_posterior[1,1][[1]])[1]
  K       = dim(A_posterior[1,1][[1]])[2]
  M       = dim(A_posterior[1,1][[1]])[3]
  
  A_out   = array(NA, c(N, K, M, S))
  for (s in 1:S) {
    A_out[,,,s]  = A_posterior[s,1][[1]]
  }
  
  return(A_out)
}


