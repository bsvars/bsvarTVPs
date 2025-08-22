

#' @title Computes posterior draws of impulse responses 
#'
#' @description Each of the draws from the posterior estimation of models from 
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the posterior distribution of the impulse responses. 
#' 
#' @param posterior posterior estimation outcome obtained by running the \code{estimate} function. 
#' The interpretation depends on the normalisation of the shocks
#' using function \code{normalise_posterior()}. Verify if the default settings are appropriate.
#' @param horizon a positive integer number denoting the forecast horizon for the impulse responses computations.
#' @param standardise a logical value. If \code{TRUE}, the impulse responses are standardised 
#' so that the variables' own shocks at horizon 0 are equal to 1. Otherwise, the parameter estimates 
#' determine this magnitude.
#' 
#' @return An object of class PosteriorIR, that is, an \code{NxNx(horizon+1)xS} array with attribute PosteriorIR 
#' containing \code{S} draws of the impulse responses.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @method compute_impulse_responses PosteriorBSVARSVMSTVI
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVMSTVI <- function(posterior, horizon, standardise = FALSE) {
  
  posterior_B = posterior$posterior$B
  posterior_A = posterior$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = dim(posterior_B[1,1][[1]])[3]
  
  # compute IRFs
  irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir_ms`, posterior_B, posterior_A, horizon, p, standardise)
  
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




#' @inherit compute_impulse_responses.PosteriorBSVARSVMSTVI
#' @method compute_impulse_responses PosteriorBSVARSVMSATVI
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVMSATVI <- function(posterior, horizon, standardise = FALSE) {
  
  posterior_B = posterior$posterior$B
  posterior_A = posterior$posterior$A
  N           = dim(posterior_B[1,1][[1]])[1]
  p           = (dim(posterior_A[1,1][[1]])[2] - 1) / N
  S           = dim(posterior_B)[1]
  M           = dim(posterior_B[1,1][[1]])[3]
  
  # compute IRFs
  irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir_mssa`, posterior_B, posterior_A, horizon, p, standardise)
  
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







#' @inherit compute_impulse_responses.PosteriorBSVARSVMSTVI
#' @method compute_impulse_responses PosteriorBSVARSVMS
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVMS <- function(posterior, horizon, standardise = FALSE) {
  
  posterior_B = posterior$posterior$B
  posterior_A = posterior$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = dim(posterior_B[1,1][[1]])[3]
  
  # compute IRFs
  irfs          = .Call(`_bsvarTVPs_bsvarTVPs_ir_ms`, posterior_B, posterior_A, horizon, p, standardise)
  
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


#' @inherit compute_impulse_responses.PosteriorBSVARSVMSTVI
#' @method compute_impulse_responses PosteriorBSVARSVTVI
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSVTVI <- function(posterior, horizon, standardise = FALSE) {
  
  posterior_B = posterior$posterior$B
  posterior_A = posterior$posterior$A
  N           = nrow(posterior_A[,,1])
  p           = (ncol(posterior_A[,,1]) - 1) / N
  S           = dim(posterior_A)[3]
  M           = 1
  
  # compute IRFs
  irfs        = .Call(`_bsvarTVPs_bsvarTVPs_ir`, posterior_B, posterior_A, horizon, p, standardise)
  
  # transform the output to an array and return
  ir_posterior    = array(NA, c(N, N, horizon + 1, M, S))
  for (s in 1:S) {
    ir_posterior[,,,1,s]  = irfs[s,1][[1]]
  }
  class(ir_posterior) = "PosteriorIR"
  
  return(ir_posterior)
}




#' @title Computes regime and TVI component-specific impulse response functions of the dependent variables 
#' to the structural shocks
#'
#' @description Computes regime and TVI component-specific impulse response functions of the dependent variables 
#' to the structural shocks for SVAR models with Markov-switching structural matrix with \code{M} regimes 
#' and the TVI structure provided on \code{VB}.
#' 
#' 
#' @param ir_posterior an object of class \code{PosteriorIR} list with random draws from the posterior distribution of 
#' the impulse responses obtained applying function \code{compute_impulse_responses}
#' @param TVI_posterior a \code{N x M x S}
#' @param VB a list containing matrices determining the TVI component specifications
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
compute_impulse_responses_by_components <- function(ir_posterior, TVI_posterior, VB) {
  
  stopifnot("Argument ir_posterior must be an object of class PosteriorIR." = class(ir_posterior) == "PosteriorIR")
  
  M           = dim(ir_posterior)[4]
  S           = dim(ir_posterior)[5]
  N           = dim(ir_posterior)[1]
  comp        = VB[length(VB)][[1]]
  TVI_ind      = which(comp > 1)
  comp        = comp[TVI_ind]
  TVI_length   = length(TVI_ind)
  
  output      = vector("list", M)
  names(output)   = paste0(rep("regime", M), 1:M)
  
  for (m in 1:M) {
    output[[m]]         = vector("list", TVI_length)
    names(output[[m]])  = paste0(rep("equation", TVI_length), TVI_ind)
    
    for (eq in 1:TVI_length) {
      output[[m]][[eq]]         = vector("list", comp[eq])
      names(output[[m]][[eq]])  = paste0(rep("TVIcomponent", comp[eq]), 1:comp[eq])
      
      for (component in 1:comp[eq]) {
        if ( M == 1 ) {
          draws_count                   = TVI_posterior[TVI_ind[eq], ] == component
        } else {
          draws_count                   = TVI_posterior[TVI_ind[eq], m, ] == component
        }
        output[[m]][[eq]][[component]]  = ir_posterior[,,,m,draws_count]
      } # END component loop
    } # END eq loop
  } # END m loop
  
  class(output)     = "PosteriorIRTVI"
  return(output)
}




#' @title Computes posterior draws of structural shocks
#'
#' @description Each of the draws from the posterior estimation of models from
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the posterior distribution of the structural shocks. 
#' 
#' @param posterior posterior estimation outcome obtained by running the \code{estimate} function. 
#' The interpretation depends on the normalisation of the shocks
#' using function \code{normalise_posterior()}. Verify if the default settings are appropriate.
#' 
#' @return An object of class PosteriorShocks, that is, an \code{NxTxS} array with attribute PosteriorShocks 
#' containing \code{S} draws of the structural shocks.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @method compute_structural_shocks PosteriorBSVARSVMSTVI
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSVMSTVI <- function(posterior) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  Y               = posterior$last_draw$data_matrices$Y
  X               = posterior$last_draw$data_matrices$X
  
  posterior_xi    = posterior$posterior$xi
  ss              = .Call(`_bsvarTVPs_bsvarTVPs_structural_shocks`, posterior_B, posterior_A, posterior_xi, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
}


#' @inherit compute_structural_shocks.PosteriorBSVARSVMSTVI
#' @method compute_structural_shocks PosteriorBSVARSVMS
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSVMS <- function(posterior) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  Y               = posterior$last_draw$data_matrices$Y
  X               = posterior$last_draw$data_matrices$X
  
  posterior_xi    = posterior$posterior$xi
  ss              = .Call(`_bsvarTVPs_bsvarTVPs_structural_shocks`, posterior_B, posterior_A, posterior_xi, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
}


#' @inherit compute_structural_shocks.PosteriorBSVARSVMSTVI
#' @method compute_structural_shocks PosteriorBSVARSVTVI
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSVTVI <- function(posterior) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  Y               = posterior$last_draw$data_matrices$Y
  X               = posterior$last_draw$data_matrices$X
  
  ss              = .Call(`_bsvarTVPs_bsvars_structural_shocks`, posterior_B, posterior_A, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
}



#' @title Computes TVI components density
#'
#' @description Computes posterior density of (Markov state-specific) TVI components
#' based on the posterior draws of the TVI indicators
#' 
#' 
#' @param posterior estimation output of class PosteriorBSVARSVMSTVI or PosteriorBSVARSVTVI
#' @param VB a list containing matrices determining the TVI component specifications
#' 
#' @return A list containing TVI components density
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
compute_TVI_component_density <- function(posterior, VB) {
  
  stopifnot(
    "Argument posterior must contain estimation output of class PosteriorBSVARSVMSTVI, PosteriorBSVARSVTVI, or PosteriorBSVARSVMSATVI." 
    = any(class(posterior)[1] == c("PosteriorBSVARSVMSTVI","PosteriorBSVARSVTVI","PosteriorBSVARSVMSATVI"))
  )
  
  UseMethod("compute_TVI_component_density", posterior)
}



#' @inherit compute_TVI_component_density
#' @method compute_TVI_component_density PosteriorBSVARSVMSTVI
#' @inheritParams compute_TVI_component_density
#' 
#' @export
compute_TVI_component_density.PosteriorBSVARSVMSTVI <- function(posterior, VB) {
  
  TVI_posterior  = posterior$posterior$S4_indicator
  
  M           = dim(TVI_posterior)[2]
  S           = dim(TVI_posterior)[3]
  N           = dim(TVI_posterior)[1]
  comp        = VB[length(VB)][[1]]
  TVI_ind      = which(comp > 1)
  comp        = comp[TVI_ind]
  TVI_length   = length(TVI_ind)
  
  TVI_density        = vector("list", TVI_length)
  names(TVI_density) = paste0(rep("equation", TVI_length), TVI_ind)
  
  for (eq in 1:TVI_length) {
    TVI_density[[eq]]            = matrix(NA, M, comp[eq])
    rownames(TVI_density[[eq]])  = paste0(rep("regime", M), 1:M)
    colnames(TVI_density[[eq]])  = paste0(rep("component", comp[eq]), 1:comp[eq])
    
    for (m in 1:M) {
      for (component in 1:comp[eq]) {
        draws_count                   = TVI_posterior[TVI_ind[eq], m, ] == component
        TVI_density[[eq]][m,component] = mean(draws_count)
      } # END component loop
    } # END m loop
  } # END eq loop
  
  class(TVI_density)     = "TVIdensity"
  
  return(TVI_density)
}



#' @inherit compute_TVI_component_density
#' @method compute_TVI_component_density PosteriorBSVARSVMSATVI
#' @inheritParams compute_TVI_component_density
#' 
#' @export
compute_TVI_component_density.PosteriorBSVARSVMSATVI <- function(posterior, VB) {
  
  TVI_posterior  = posterior$posterior$S4_indicator
  
  M           = dim(TVI_posterior)[2]
  S           = dim(TVI_posterior)[3]
  N           = dim(TVI_posterior)[1]
  comp        = VB[length(VB)][[1]]
  TVI_ind      = which(comp > 1)
  comp        = comp[TVI_ind]
  TVI_length   = length(TVI_ind)
  
  TVI_density        = vector("list", TVI_length)
  names(TVI_density) = paste0(rep("equation", TVI_length), TVI_ind)
  
  for (eq in 1:TVI_length) {
    TVI_density[[eq]]            = matrix(NA, M, comp[eq])
    rownames(TVI_density[[eq]])  = paste0(rep("regime", M), 1:M)
    colnames(TVI_density[[eq]])  = paste0(rep("component", comp[eq]), 1:comp[eq])
    
    for (m in 1:M) {
      for (component in 1:comp[eq]) {
        draws_count                   = TVI_posterior[TVI_ind[eq], m, ] == component
        TVI_density[[eq]][m,component] = mean(draws_count)
      } # END component loop
    } # END m loop
  } # END eq loop
  
  class(TVI_density)     = "TVIdensity"
  
  return(TVI_density)
}





#' @inherit compute_TVI_component_density
#' @method compute_TVI_component_density PosteriorBSVARSVTVI
#' @inheritParams compute_TVI_component_density
#' 
#' @export
compute_TVI_component_density.PosteriorBSVARSVTVI <- function(posterior, VB) {
  
  TVI_posterior  = posterior$posterior$S4_indicator
  
  S           = dim(TVI_posterior)[2]
  N           = dim(TVI_posterior)[1]
  comp        = VB[length(VB)][[1]]
  TVI_ind      = which(comp > 1)
  comp        = comp[TVI_ind]
  TVI_length   = length(TVI_ind)
  
  TVI_density        = vector("list", TVI_length)
  names(TVI_density) = paste0(rep("equation", TVI_length), TVI_ind)
  
  for (eq in 1:TVI_length) {
    TVI_density[[eq]]            = matrix(NA, 1, comp[eq])
    colnames(TVI_density[[eq]])  = paste0(rep("component", comp[eq]), 1:comp[eq])
    
      for (component in 1:comp[eq]) {
        draws_count                   = TVI_posterior[TVI_ind[eq], ] == component
        TVI_density[[eq]][1,component] = mean(draws_count)
      } # END component loop
  } # END eq loop
  
  class(TVI_density)     = "TVIdensity"
  
  return(TVI_density)
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
#' @method structural_to_array PosteriorBSVARSVMSTVI
#' @inheritParams structural_to_array
#' 
#' @return An \code{N x N x M x S} array containing the posterior draws 
#' of the structural matrix \eqn{B}
#' 
#' @export
structural_to_array.PosteriorBSVARSVMSTVI <- function(posterior) {
  
  B_posterior = posterior$posterior$B
  S       = dim(B_posterior)[1]
  N       = dim(B_posterior[1,1][[1]])[1]
  M       = dim(B_posterior[1,1][[1]])[3]
  
  B_out   = array(NA, c(N, N, M, S))
  for (s in 1:S) {
    B_out[,,,s]  = B_posterior[s,1][[1]]
  }
  
  return(B_out)
}



#' @inherit structural_to_array
#' @method structural_to_array PosteriorBSVARSVMSATVI
#' @inheritParams structural_to_array
#' 
#' @return An \code{N x N x M x S} array containing the posterior draws 
#' of the structural matrix \eqn{B}
#' 
#' @export
structural_to_array.PosteriorBSVARSVMSATVI <- function(posterior) {
  
  B_posterior = posterior$posterior$B
  S       = dim(B_posterior)[1]
  N       = dim(B_posterior[1,1][[1]])[1]
  M       = dim(B_posterior[1,1][[1]])[3]
  
  B_out   = array(NA, c(N, N, M, S))
  for (s in 1:S) {
    B_out[,,,s]  = B_posterior[s,1][[1]]
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
    B_out[,,,s]  = B_posterior[s,1][[1]]
  }
  
  return(B_out)
}


#' @inherit structural_to_array
#' @method structural_to_array PosteriorBSVARSVTVI
#' @inheritParams structural_to_array
#' 
#' @return An \code{N x N x S} array containing the posterior draws 
#' of the structural matrix \eqn{B}
#' 
#' @export
structural_to_array.PosteriorBSVARSVTVI <- function(posterior) {
  
  return(posterior$posterior$B)
}



#' @title Extracts the posterior draws of the autoregressive matrix \eqn{A} and returns the as an \code{array}
#'
#' @description Extracts the posterior draws of the autoregressive matrix \eqn{A} and returns the as an \code{array}
#' 
#' @param posterior Posterior draws of from a Markov-switching model
#' 
#' @return An array containing the posterior draws of the autoregressive matrix \eqn{A}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
autoregressive_to_array <- function(posterior) {
  
  # call method
  UseMethod("autoregressive_to_array", posterior)
}


#' @inherit autoregressive_to_array
#' @method autoregressive_to_array PosteriorBSVARSVMSATVI
#' @inheritParams autoregressive_to_array
#' 
#' @return An \code{N x K x M x S} array containing the posterior draws 
#' of the autoregressive matrix \eqn{A}
#' 
#' @export
autoregressive_to_array.PosteriorBSVARSVMSATVI <- function(posterior) {
  
  A_posterior = posterior$posterior$A
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







#' @title Computes the posterior draws of data conditional covariances, correlations,
#' or standard deviations
#'
#' @description Computes the posterior draws of time-varying conditional 
#' covariances, correlations, or standard deviations given the past observations
#' 
#' @param posterior Estimation output of a Structural VAR model. An object of class 
#' PosteriorBSVARSVMSTVI, PosteriorBSVARSVTVI, PosteriorBSVARSVMS, PosteriorBSVARMSTVI, 
#' PosteriorBSVARTVI, or PosteriorBSVARMS
#' @param moment one of the values \code{"cov"}, \code{"cor"}, or \code{"sd"} denoting
#' the second-order moments to compute
#' @param thin a positive integer determining MCMC thinning used here to bypass the memory problems
#' 
#' @return An array of size \code{N x N x T x S} for \code{moment = "cov"} or \code{moment = "cor"},
#' or of size \code{N x T x S} for \code{moment = "sd"} with the draws of the selected
#' time-varying second-order moment
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
compute_conditional_cov <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  # call method
  UseMethod("compute_conditional_cov")
}


#' @inherit compute_conditional_cov
#' @method compute_conditional_cov PosteriorBSVARSVMSTVI
#' @inheritParams compute_conditional_cov
#' 
#' @export
compute_conditional_cov.PosteriorBSVARSVMSTVI <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  S_orig          = dim(posterior$posterior$xi)[3]
  keep            = seq(from = 1, to = S_orig, by = thin)
  
  posterior_B     = posterior$posterior$B[keep,]
  posterior_xi    = posterior$posterior$xi[,,keep]
  posterior_sigma = posterior$posterior$sigma[,,keep]
  N               = dim(posterior_sigma)[1]
  T               = dim(posterior_xi)[2]
  S               = dim(posterior_xi)[3]
  
  # compute conditional covariances
  covs_tmp        = .Call(`_bsvarTVPs_bsvarTVPs_covariances_rf_mssv`, posterior_B, posterior_xi, posterior_sigma)
  covs_out        = array(NA, c(N, N, T, S))
  
  if ( moment == "cor" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2cor`, covs_tmp)
  } else if ( moment == "sd" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2sd`, covs_tmp)
    covs_out      = covs_tmp
  }
  
  if ( moment == "cov" || moment == "cor" ) {
    for (s in 1:S) {
      for (t in 1:T) {
        covs_out[,,,s] = covs_tmp[s][[1]]
      } 
    }
  }
  class(covs_out) = "PosteriorMoments"
  return(covs_out)
}


#' @inherit compute_conditional_cov
#' @method compute_conditional_cov PosteriorBSVARSVMS
#' @inheritParams compute_conditional_cov
#' 
#' @export
compute_conditional_cov.PosteriorBSVARSVMS <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  S_orig          = dim(posterior$posterior$xi)[3]
  keep            = seq(from = 1, to = S_orig, by = thin)
  
  posterior_B     = posterior$posterior$B[keep,]
  posterior_xi    = posterior$posterior$xi[,,keep]
  posterior_sigma = posterior$posterior$sigma[,,keep]
  N               = dim(posterior_sigma)[1]
  T               = dim(posterior_xi)[2]
  S               = dim(posterior_xi)[3]
  
  # compute conditional covariances
  covs_tmp        = .Call(`_bsvarTVPs_bsvarTVPs_covariances_rf_mssv`, posterior_B, posterior_xi, posterior_sigma)
  covs_out        = array(NA, c(N, N, T, S))
  
  if ( moment == "cor" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2cor`, covs_tmp)
  } else if ( moment == "sd" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2sd`, covs_tmp)
    covs_out      = covs_tmp
  }
  
  if ( moment == "cov" || moment == "cor" ) {
    for (s in 1:S) {
      for (t in 1:T) {
        covs_out[,,,s] = covs_tmp[s][[1]]
      } 
    }
  }
  class(covs_out) = "PosteriorMoments"
  return(covs_out)
}


#' @inherit compute_conditional_cov
#' @method compute_conditional_cov PosteriorBSVARSVTVI
#' @inheritParams compute_conditional_cov
#' 
#' @export
compute_conditional_cov.PosteriorBSVARSVTVI <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  S_orig          = dim(posterior$posterior$sigma)[3]
  keep            = seq(from = 1, to = S_orig, by = thin)
  
  posterior_B     = posterior$posterior$B[,,keep]
  posterior_sigma = posterior$posterior$sigma[,,keep]
  N               = dim(posterior_sigma)[1]
  T               = dim(posterior_sigma)[2]
  S               = dim(posterior_sigma)[3]
  
  # compute conditional covariances
  covs_tmp        = .Call(`_bsvarTVPs_bsvarTVPs_covariances_rf_sv`, posterior_B, posterior_sigma)
  covs_out        = array(NA, c(N, N, T, S))
  
  if ( moment == "cor" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2cor`, covs_tmp)
  } else if ( moment == "sd" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2sd`, covs_tmp)
    covs_out      = covs_tmp
  }
  
  if ( moment == "cov" || moment == "cor" ) {
    for (s in 1:S) {
      for (t in 1:T) {
        covs_out[,,,s] = covs_tmp[s][[1]]
      } 
    }
  }
  class(covs_out) = "PosteriorMoments"
  return(covs_out)
}


#' @inherit compute_conditional_cov
#' @method compute_conditional_cov PosteriorBSVARMSTVI
#' @inheritParams compute_conditional_cov
#' 
#' @export
compute_conditional_cov.PosteriorBSVARMSTVI <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  S_orig          = dim(posterior$posterior$xi)[3]
  keep            = seq(from = 1, to = S_orig, by = thin)
  
  posterior_B     = posterior$posterior$B[keep,]
  posterior_xi    = posterior$posterior$xi[,,keep]
  N               = dim(posterior$posterior$B[1,1][[1]])[1]
  T               = dim(posterior_xi)[2]
  S               = dim(posterior_xi)[3]
  
  # compute conditional covariances
  covs_tmp        = .Call(`_bsvarTVPs_bsvarTVPs_covariances_rf_ms`, posterior_B, posterior_xi)
  covs_out        = array(NA, c(N, N, T, S))
  
  if ( moment == "cor" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2cor`, covs_tmp)
  } else if ( moment == "sd" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2sd`, covs_tmp)
    covs_out      = covs_tmp
  }
  
  if ( moment == "cov" || moment == "cor" ) {
    for (s in 1:S) {
      for (t in 1:T) {
        covs_out[,,,s] = covs_tmp[s][[1]]
      } 
    }
  }
  class(covs_out) = "PosteriorMoments"
  return(covs_out)
}


#' @inherit compute_conditional_cov
#' @method compute_conditional_cov PosteriorBSVARMS
#' @inheritParams compute_conditional_cov
#' 
#' @export
compute_conditional_cov.PosteriorBSVARMS <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  S_orig          = dim(posterior$posterior$xi)[3]
  keep            = seq(from = 1, to = S_orig, by = thin)
  
  posterior_B     = posterior$posterior$B[keep,]
  posterior_xi    = posterior$posterior$xi[,,keep]
  N               = dim(posterior$posterior$B[1,1][[1]])[1]
  T               = dim(posterior_xi)[2]
  S               = dim(posterior_xi)[3]
  
  # compute conditional covariances
  covs_tmp        = .Call(`_bsvarTVPs_bsvarTVPs_covariances_rf_ms`, posterior_B, posterior_xi)
  covs_out        = array(NA, c(N, N, T, S))
  
  if ( moment == "cor" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2cor`, covs_tmp)
  } else if ( moment == "sd" ) {
    covs_tmp      = .Call(`_bsvarTVPs_bsvarTVPs_cov2sd`, covs_tmp)
    covs_out      = covs_tmp
  }
  
  if ( moment == "cov" || moment == "cor" ) {
    for (s in 1:S) {
      for (t in 1:T) {
        covs_out[,,,s] = covs_tmp[s][[1]]
      } 
    }
  }
  class(covs_out) = "PosteriorMoments"
  return(covs_out)
}


#' @inherit compute_conditional_cov
#' @method compute_conditional_cov PosteriorBSVARTVI
#' @inheritParams compute_conditional_cov
#' 
#' @export
compute_conditional_cov.PosteriorBSVARTVI <- function(posterior, moment = c("cov", "cor", "sd"), thin = 1L) {
  
  moment <- match.arg(moment)
  
  S_orig          = dim(posterior$posterior$B)[3]
  keep            = seq(from = 1, to = S_orig, by = thin)
  
  posterior_B     = posterior$posterior$B[,,keep]
  N               = dim(posterior_B)[1]
  S               = dim(posterior_B)[3]
  
  # compute conditional covariances
  covs_tmp        = .Call(`_bsvarTVPs_bsvarTVPs_covariances_rf`, posterior_B)
  
  if ( moment == "cor" ) {
    covs_tmp      = array(apply(covs_tmp, 3, stats::cov2cor), c(N,N,S))
  } else if ( moment == "sd" ) {
    covs_tmp      = sapply(1:S, function(s){sqrt(diag(covs_tmp[,,s]))})
  }
  covs_out      = covs_tmp
  
  class(covs_out) = "PosteriorMoments"
  return(covs_out)
}


