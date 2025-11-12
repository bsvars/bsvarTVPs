

#' @title Computes posterior draws of impulse responses 
#'
#' @description Each of the draws from the posterior estimation of models from 
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the posterior distribution of the impulse responses. 
#' 
#' @param posterior posterior estimation outcome of class \code{PosteriorBSVARTVP} 
#' obtained by running the \code{estimate} function. 
#' The interpretation depends on the normalisation of the shocks
#' using function \code{normalise()}. Verify if the default settings are appropriate.
#' @param horizon a positive integer number denoting the forecast horizon for 
#' the impulse responses computations.
#' @param standardise a logical value. If \code{TRUE}, the impulse responses are 
#' standardised so that the variables' own shocks at horizon 0 are equal to 1. 
#' Otherwise, the parameter estimates determine this magnitude.
#' 
#' @return An object of class \code{PosteriorIR}, that is, an 
#' \code{NxNx(horizon+1)xS} array with attribute PosteriorIR 
#' containing \code{S} draws of the impulse responses.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @method compute_impulse_responses PosteriorBSVARTVP
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec   = specify_bsvarTVP$new(us_fiscal_lsuw)    # specify the model
#' burn   = estimate(spec, 5)                       # run the burn-in for convergence
#' post   = estimate(burn, 5)                       # estimate the model
#' irf    = compute_impulse_responses(post, horizon = 4) # compute impulse responses
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_fiscal_lsuw |>
#'   specify_bsvarTVP$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 5) -> post
#' post |> compute_impulse_responses(horizon = 4) -> irf
#' 
#' @export
compute_impulse_responses.PosteriorBSVARTVP <- function(posterior, horizon, standardise = FALSE) {
  
  posterior_B = posterior$posterior$B_cpp
  posterior_A = posterior$posterior$A_cpp
  N           = nrow(posterior_B[1,1][[1]])
  p           = posterior$last_draw$get_p()
  S           = dim(posterior_B)[1]
  M           = posterior$last_draw$get_M()
  
  # compute IRFs
  if ( posterior$last_draw$get_msa() ) {
    irfs      = .Call(`_bsvarTVPs_bsvarTVPs_ir_mssa`, posterior_B, posterior_A, horizon, p, standardise)
  } else {
    irfs      = .Call(`_bsvarTVPs_bsvarTVPs_ir_ms`, posterior_B, posterior_A, horizon, p, standardise)
  }
  
  # transform the output to an array and return
  ir_posterior    = array(NA, c(N, N, horizon + 1, M, S))
  for (s in 1:S) {
    for (m in 1:M) {
      ir_posterior[,,,m,s]  = irfs[s,m][[1]]
    }
  }
  class(ir_posterior) = "PosteriorIR"
  
  return(ir_posterior)
}


















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


