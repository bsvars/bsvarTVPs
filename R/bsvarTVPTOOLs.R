#' @title Computes regime-specific impulse response functions of the dependent variables 
#' to the structural shocks
#'
#' @description Computes regime-specific impulse response functions of the dependent variables 
#' to the structural shocks for SVAR models with Markov-switching structural matrix with \code{M} regimes
#' 
#' 
#' @param posterior a list with random draws from the posterior distribution of 
#' the Structural model with Markov-switching structural matrix
#' @param horizon a positive integer specifying the forecasting horizon for the impulse responses
#' 
#' @return An \code{N x N x horizon x M x S} array containing the posterior draws 
#' of the impulse responses 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvarTVPs_ir <- function(posterior, horizon) {
  
  stopifnot("Argument horizon must be a positive integer number." = horizon >= 1 & horizon %% 1 == 0)
  
  posterior_B = posterior$B
  posterior_A = posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  M           = dim(posterior_B[1,1][[1]])[3]
  S           = dim(posterior_A)[3]
  
  # compute IRFs
  irfs            = .Call(`_bsvarTVPs_bsvarTVPs_ir`, posterior_B, posterior_A, horizon, p)
  
  # transform the output to an array and return
  ir_posterior    = array(NA, c(N, N, horizon + 1, M, S))
  for (s in 1:S) {
    for (m in 1:M) {
      ir_posterior[,,,m,s]  = irfs[s,m][[1]]
    }
  }
  
  return(ir_posterior)
}
