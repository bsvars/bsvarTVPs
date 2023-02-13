

#' @title Computes regime-specific impulse response functions of the dependent variables 
#' to the structural shocks
#'
#' @description Computes regime-specific impulse response functions of the dependent variables 
#' to the structural shocks for SVAR models with Markov-switching structural matrix with \code{M} regimes
#' 
#' 
#' @param posterior_output an object of class of class PosteriorBSVARSVMSS5, 
#' PosteriorBSVARSVMS, or PosteriorBSVARSVS5 with random draws from the posterior 
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
  
  # check arguments
  stopifnot("Argument posterior must be of class PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5." = substr(class(posterior_output), 1, 16) == "PosteriorBSVARSV")
  stopifnot("Argument horizon must be a positive integer number." = horizon >= 1 & horizon %% 1 == 0)
  
  # call method
  UseMethod("compute_impulse_responses", posterior_output)
}


#' @inherit compute_impulse_responses
#' @method compute_impulse_responses PosteriorBSVARSVMSS5
#' @inheritParams compute_impulse_responses
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVMSS5 <- function(posterior_output, horizon) {
  
  posterior_B = posterior_output$posterior$B
  posterior_A = posterior_output$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = dim(posterior_B[1,1][[1]])[3]
  
  # compute IRFs
  irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir_ms`, posterior_B, posterior_A, horizon, p)
  
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


#' @inherit compute_impulse_responses
#' @method compute_impulse_responses PosteriorBSVARSVMS
#' @inheritParams compute_impulse_responses
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVMS <- function(posterior_output, horizon) {
  
  posterior_B = posterior_output$posterior$B
  posterior_A = posterior_output$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = dim(posterior_B[1,1][[1]])[3]
  
  # compute IRFs
  irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir_ms`, posterior_B, posterior_A, horizon, p)
  
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


#' @inherit compute_impulse_responses
#' @method compute_impulse_responses PosteriorBSVARSVS5
#' @inheritParams compute_impulse_responses
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVS5 <- function(posterior_output, horizon) {
  
  posterior_B = posterior_output$posterior$B
  posterior_A = posterior_output$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = 1
  
  # compute IRFs
  irfs        = .Call(`_bsvarTVPs_bsvarTVPs_ir`, posterior_B, posterior_A, horizon, p)
  
  # transform the output to an array and return
  ir_posterior    = array(NA, c(N, N, horizon + 1, M, S))
  for (s in 1:S) {
    ir_posterior[,,,1,s]  = irfs[s,1][[1]]
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
        if ( M == 1 ) {
          draws_count                   = S5_posterior[S5_ind[eq], ] == component
        } else {
          draws_count                   = S5_posterior[S5_ind[eq], m, ] == component
        }
        output[[m]][[eq]][[component]]  = ir_posterior[,,,m,draws_count]
        S5_density[[m]][[eq]][component] = mean(draws_count)
      } # END component loop
    } # END eq loop
  } # END m loop
  
  output$density    = S5_density
  class(output)     = "PosteriorIRS5"
  
  return(output)
}


#' @title Computes posterior draws of regime probabilities
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the regime probabilities. These represent either
#' the realisations of the regime indicators, when \code{type = "realized"}, filtered probabilities,
#' when \code{type = "filtered"}, forecasted regime probabilities, when \code{type = "forecasted"},
#' or the smoothed probabilities, when \code{type = "smoothed"}, .
#' 
#' @param posterior posterior estimation outcome of regime-dependent heteroskedastic models 
#' - an object of either of the classes: PosteriorBSVARSVMSS5, or PosteriorBSVARSVMS.
#' @param Y a \code{N x T} matrix of dependent variables
#' @param X a \code{K x T} matrix of lagged variables
#' @param type one of the values \code{"realized"}, \code{"filtered"}, \code{"forecasted"}, or \code{"smoothed"}
#' denoting the type of probabilities to be computed.
#' 
#' @return An object of class PosteriorRegimePr, that is, an \code{MxTxS} array with attribute PosteriorRegimePr 
#' containing \code{S} draws of the regime probabilities.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Song, Y., and Woźniak, T., (2021) Markov Switching. \emph{Oxford Research Encyclopedia of Economics and Finance}, Oxford University Press, \doi{https://doi.org/10.1093/acrefore/9780190625979.013.174}.
#'  
#' @export
compute_regime_probabilities <- function(posterior, Y, X, type = c("realized", "filtered", "forecasted", "smoothed")) {
  
  stopifnot("Argument posterior must be of class PosteriorBSVARSVMSS5 or PosteriorBSVARSVMS." = substr(class(posterior), 1, 18) == "PosteriorBSVARSVMS")
  
  type          = match.arg(type)
  
  posteriors    = posterior$posterior
  # Y             = posterior$last_draw$data_matrices$Y
  # X             = posterior$last_draw$data_matrices$X
  
  if (type == "realized") {
    probs       = posteriors$xi
  } else {
    if (type == "filtered") {
      forecasted  = FALSE
      smoothed    = FALSE
    } else if (type == "forecasted") {
      forecasted  = TRUE
      smoothed    = FALSE
    } else if (type == "smoothed") {
      forecasted  = FALSE
      smoothed    = TRUE
    }
    
    probs       = .Call(`_bsvarTVPs_bsvarTVPs_filter_forecast_smooth`, posteriors, Y, X, forecasted, smoothed)
  }
  
  class(probs)  = "PosteriorRegimePr"
  
  return(probs)
}


#' @title Computes posterior draws of dependent variables' fitted values
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the fitted values. 
#' 
#' @param posterior posterior estimation outcome - an object of either of the classes: 
#' PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5.
#' @param X a \code{K x T} matrix of vector autoregressive regressors.
#' 
#' @return An object of class PosteriorFitted, that is, an \code{NxTxS} array with attribute PosteriorFitted 
#' containing \code{S} draws of the fitted values.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @export
compute_fitted_values <- function(posterior, X) {
  
  stopifnot("Argument posterior must contain estimation output." = any(class(posterior)[1] == c("PosteriorBSVARSVMSS5", "PosteriorBSVARSVMS", "PosteriorBSVARSVS5")))
  
  posterior_A     = posterior$posterior$A
  # X               = posterior$last_draw$data_matrices$X
  
  fv              = .Call(`_bsvarTVPs_bsvarTVPs_fitted_values`, posterior_A, X)
  class(fv)       = "PosteriorFitted"
  
  return(fv)
}


#' @title Computes posterior draws of structural shocks
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the structural shocks. 
#' 
#' @param posterior posterior estimation outcome - an object of either of the classes: 
#' PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5.
#' @param Y a \code{N x T} matrix od dependent variables
#' @param X a \code{K x T} matrix of vector autoregressive regressors.
#' 
#' @return An object of class PosteriorShocks, that is, an \code{NxTxS} array with attribute PosteriorShocks 
#' containing \code{S} draws of the structural shocks.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @export
compute_structural_shocks <- function(posterior, Y, X) {
  
  stopifnot("Argument posterior must contain estimation output." = any(class(posterior)[1] == c("PosteriorBSVARSVMSS5", "PosteriorBSVARSVMS", "PosteriorBSVARSVS5")))
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  # Y               = posterior$last_draw$data_matrices$Y
  # X               = posterior$last_draw$data_matrices$X
  
  if ( substr(class(posterior), 1, 18) == "PosteriorBSVARSVMS" ) {
    posterior_xi    = posterior$posterior$xi
    ss              = .Call(`_bsvarTVPs_bsvarTVPs_structural_shocks`, posterior_B, posterior_A, posterior_xi, Y, X)
  } else {
    ss              = .Call(`_bsvarTVPs_bsvars_structural_shocks`, posterior_B, posterior_A, Y, X)
  }
  class(ss)       = "PosteriorShocks"
  
  return(ss)
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
        S5_density[[m]][[eq]][component] = mean(draws_count)
      } # END component loop
    } # END eq loop
  } # END m loop
  
  output$density    = S5_density
  class(output)     = "PosteriorIRS5"
  
  return(output)
}


#' @title Computes posterior draws of regime probabilities
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the regime probabilities. These represent either
#' the realisations of the regime indicators, when \code{type = "realized"}, filtered probabilities,
#' when \code{type = "filtered"}, forecasted regime probabilities, when \code{type = "forecasted"},
#' or the smoothed probabilities, when \code{type = "smoothed"}, .
#' 
#' @param posterior posterior estimation outcome of regime-dependent heteroskedastic models 
#' - an object of either of the classes: PosteriorBSVARSVMSS5, or PosteriorBSVARSVMS.
#' @param Y a \code{N x T} matrix of dependent variables
#' @param X a \code{K x T} matrix of lagged variables
#' @param type one of the values \code{"realized"}, \code{"filtered"}, \code{"forecasted"}, or \code{"smoothed"}
#' denoting the type of probabilities to be computed.
#' 
#' @return An object of class PosteriorRegimePr, that is, an \code{MxTxS} array with attribute PosteriorRegimePr 
#' containing \code{S} draws of the regime probabilities.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Song, Y., and Woźniak, T., (2021) Markov Switching. \emph{Oxford Research Encyclopedia of Economics and Finance}, Oxford University Press, \doi{https://doi.org/10.1093/acrefore/9780190625979.013.174}.
#'  
#' @export
compute_regime_probabilities <- function(posterior, Y, X, type = c("realized", "filtered", "forecasted", "smoothed")) {
  
  stopifnot("Argument posterior must be of class PosteriorBSVARSVMSS5 or PosteriorBSVARSVMS." = substr(class(posterior), 1, 18) == "PosteriorBSVARSVMS")
  
  type          = match.arg(type)
  
  posteriors    = posterior$posterior
  # Y             = posterior$last_draw$data_matrices$Y
  # X             = posterior$last_draw$data_matrices$X
  
  if (type == "realized") {
    probs       = posteriors$xi
  } else {
    if (type == "filtered") {
      forecasted  = FALSE
      smoothed    = FALSE
    } else if (type == "forecasted") {
      forecasted  = TRUE
      smoothed    = FALSE
    } else if (type == "smoothed") {
      forecasted  = FALSE
      smoothed    = TRUE
    }
    
    probs       = .Call(`_bsvarTVPs_bsvarTVPs_filter_forecast_smooth`, posteriors, Y, X, forecasted, smoothed)
  }
  
  class(probs)  = "PosteriorRegimePr"
  
  return(probs)
}


#' @title Computes posterior draws of dependent variables' fitted values
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the fitted values. 
#' 
#' @param posterior posterior estimation outcome - an object of either of the classes: 
#' PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5.
#' @param X a \code{K x T} matrix of vector autoregressive regressors.
#' 
#' @return An object of class PosteriorFitted, that is, an \code{NxTxS} array with attribute PosteriorFitted 
#' containing \code{S} draws of the fitted values.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @export
compute_fitted_values <- function(posterior, X) {
  
  stopifnot("Argument posterior must contain estimation output." = any(class(posterior)[1] == c("PosteriorBSVARSVMSS5", "PosteriorBSVARSVMS", "PosteriorBSVARSVS5")))
  
  posterior_A     = posterior$posterior$A
  # X               = posterior$last_draw$data_matrices$X
  
  fv              = .Call(`_bsvarTVPs_bsvarTVPs_fitted_values`, posterior_A, X)
  class(fv)       = "PosteriorFitted"
  
  return(fv)
}


#' @title Computes posterior draws of structural shocks
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the structural shocks. 
#' 
#' @param posterior posterior estimation outcome - an object of either of the classes: 
#' PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5.
#' @param Y a \code{N x T} matrix od dependent variables
#' @param X a \code{K x T} matrix of vector autoregressive regressors.
#' 
#' @return An object of class PosteriorShocks, that is, an \code{NxTxS} array with attribute PosteriorShocks 
#' containing \code{S} draws of the structural shocks.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @export
compute_structural_shocks <- function(posterior, Y, X) {
  
  stopifnot("Argument posterior must contain estimation output." = any(class(posterior)[1] == c("PosteriorBSVARSVMSS5", "PosteriorBSVARSVMS", "PosteriorBSVARSVS5")))
  
  # call method
  UseMethod("compute_structural_shocks", posterior)
}




#' @inherit compute_structural_shocks
#' @method compute_structural_shocks PosteriorBSVARSVMSS5
#' @inheritParams compute_structural_shocks
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSVMSS5 <- function(posterior, Y, X) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  
  posterior_xi    = posterior$posterior$xi
  ss              = .Call(`_bsvarTVPs_bsvarTVPs_structural_shocks`, posterior_B, posterior_A, posterior_xi, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
}


#' @inherit compute_structural_shocks
#' @method compute_structural_shocks PosteriorBSVARSVMS
#' @inheritParams compute_structural_shocks
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSVMS <- function(posterior, Y, X) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  
  posterior_xi    = posterior$posterior$xi
  ss              = .Call(`_bsvarTVPs_bsvarTVPs_structural_shocks`, posterior_B, posterior_A, posterior_xi, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
}


#' @inherit compute_structural_shocks
#' @method compute_structural_shocks PosteriorBSVARSVS5
#' @inheritParams compute_structural_shocks
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSVS5 <- function(posterior, Y, X) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  ss              = .Call(`_bsvarTVPs_bsvars_structural_shocks`, posterior_B, posterior_A, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
}




#' @title Extracts the posterior draws of the structural matrix \eqn{B} and returns the as an \code{array}
#'
#' @description Extracts the posterior draws of the structural matrix \eqn{B} and returns the as an \code{array}
#' 
#' @param posterior Posterior draws of from a Markov-switching model
#' 
#' @return An array containing the posterior draws of the structural matrix \eqn{B}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
structural_to_array <- function(posterior) {
  
  # call method
  UseMethod("structural_to_array", posterior)
}


#' @inherit structural_to_array
#' @method structural_to_array PosteriorBSVARSVMSS5
#' @inheritParams structural_to_array
#' 
#' @return An \code{N x N x M x S} array containing the posterior draws 
#' of the structural matrix \eqn{B}
#' 
#' @export
structural_to_array.PosteriorBSVARSVMSS5 <- function(posterior) {
  
  B_posterior = posterior$posterior$B
  S       = dim(B_posterior)[1]
  N       = dim(B_posterior[1,1][[1]])[1]
  M       = dim(B_posterior[1,1][[1]])[3]
  
  B_out   = array(NA, c(N, N, M, S))
  for (s in 1:S) {
    B_out[,,,s]  = B_posterior[1,1][[1]]
  }
  
  return(B_out)
}


#' @inherit structural_to_array
#' @method structural_to_array PosteriorBSVARSVMS
#' @inheritParams structural_to_array
#' 
#' @return An \code{N x N x M x S} array containing the posterior draws 
#' of the structural matrix \eqn{B}
#' 
#' @export
structural_to_array.PosteriorBSVARSVMS <- function(posterior) {
  
  B_posterior = posterior$posterior$B
  S       = dim(B_posterior)[1]
  N       = dim(B_posterior[1,1][[1]])[1]
  M       = dim(B_posterior[1,1][[1]])[3]
  
  B_out   = array(NA, c(N, N, M, S))
  for (s in 1:S) {
    B_out[,,,s]  = B_posterior[1,1][[1]]
  }
  
  return(B_out)
}