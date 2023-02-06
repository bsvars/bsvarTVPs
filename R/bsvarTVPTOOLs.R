#' @title Computes regime-specific impulse response functions of the dependent variables 
#' to the structural shocks
#'
#' @description Computes regime-specific impulse response functions of the dependent variables 
#' to the structural shocks for SVAR models with Markov-switching structural matrix with \code{M} regimes
#' 
#' 
#' @param posterior_output an object of classof class PosteriorBSVARSV-MSS5, 
#' PosteriorBSVARSV-MS, or PosteriorBSVARSV-S5 with random draws from the posterior 
#' distribution of the Structural model with Markov-switching structural matrix
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
compute_impulse_responses <- function(posterior_output, horizon) {
  
  stopifnot("Argument posterior must be of class PosteriorBSVARSV-MSS5, PosteriorBSVARSV-MS, or PosteriorBSVARSV-S5." = substr(class(posterior_output), 1, 17) == "PosteriorBSVARSV-")
  stopifnot("Argument horizon must be a positive integer number." = horizon >= 1 & horizon %% 1 == 0)
  
  is_MS       = substr(class(posterior_output), 1, 19) == "PosteriorBSVARSV-MS"
  
  posterior_B = posterior_output$posterior$B
  posterior_A = posterior_output$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = 1
  
  if ( is_MS ) M  = dim(posterior_B[1,1][[1]])[3]
  
  # compute IRFs
  if ( is_MS ) {
    irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir_ms`, posterior_B, posterior_A, horizon, p)
  } else {
    irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir`, posterior_B, posterior_A, horizon, p)
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


#' @title Computes regime and S5 component-specific impulse response functions of the dependent variables 
#' to the structural shocks
#'
#' @description Computes regime and S5 component-specific impulse response functions of the dependent variables 
#' to the structural shocks for SVAR models with Markov-switching structural matrix with \code{M} regimes 
#' and the S5 structure provided on \code{VB}.
#' 
#' 
#' @param ir_posterior an object of class \code{PosteriorIR} list with random draws from the posterior distribution of 
#' the impulse responses obtained applying function \code{compute_impulse_responses}
#' @param S5_posterior a \code{N x M x S}
#' @param VB a list containing matrices determining the S5 component specifications
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
compute_impulse_responses_by_components <- function(ir_posterior, S5_posterior, VB) {
  
  stopifnot("Argument ir_posterior must be an object of class PosteriorIR." = class(ir_posterior) == "PosteriorIR")
  
  M           = dim(ir_posterior)[4]
  S           = dim(ir_posterior)[5]
  N           = dim(ir_posterior)[1]
  comp        = VB[length(VB)][[1]]
  S5_ind      = which(comp > 1)
  comp        = comp[S5_ind]
  S5_length   = length(S5_ind)
  
  output      = vector("list", M)
  names(output)   = paste0(rep("regime", M), 1:M)
  S5_density  = vector("list", M)
  names(S5_density)   = paste0(rep("regime", M), 1:M)
  
  for (m in 1:M) {
    output[[m]]         = vector("list", S5_length)
    names(output[[m]])  = paste0(rep("equation", S5_length), S5_ind)
    S5_density[[m]]     = vector("list", S5_length)
    names(S5_density[[m]]) = paste0(rep("equation", S5_length), S5_ind)
    
    for (eq in 1:S5_length) {
      output[[m]][[eq]]         = vector("list", comp[eq])
      names(output[[m]][[eq]])  = paste0(rep("S5component", comp[eq]), 1:comp[eq])
      S5_density[[m]][[eq]]     = rep(NA, comp[eq])
      names(S5_density[[m]][[eq]])      = 1:comp[eq]
      
      for (component in 1:comp[eq]) {
        draws_count                     = S5_posterior[S5_ind[eq], m, ] == component
        output[[m]][[eq]][[component]]  = ir_posterior[,,,m,draws_count]
        S5_density[[m]][[eq]][component]= mean(draws_count)
      } # END component loop
    } # END eq loop
  } # END m loop
  
  output$density    = S5_density
  class(output)     = "PosteriorIRS5"
  
  return(output)
}
