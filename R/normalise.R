#' @title Waggoner & Zha (2003) row signs normalisation of the posterior draws for matrix \eqn{B}
#'
#' @description Normalises the sign of rows of matrix \eqn{B} MCMC draws 
#' Markov state-by-state and S5 component-by-component, 
#'  provided as the first argument \code{posterior_B}, relative to the matrices in
#'  \code{B_hat}, provided as the second argument of the function. If the second argument 
#'  is not provided, the function creates its own benchmark matrix. The implemented
#'  procedure proposed by Waggoner, Zha (2003) normalises the MCMC output in an
#'  optimal way leading to the unimodal posterior. Only normalised MCMC output is 
#'  suitable for the computations of the posterior characteristics of the \eqn{B}
#'  matrix elements and their functions such as the impulse response functions and other 
#'  economically interpretable values. 
#' 
#' @param posterior posterior estimation outcome - an object of either of classes: 
#' "PosteriorBSVARSVMSS5", "PosteriorBSVARSVMS", or "PosteriorBSVARSVS5"
#' containing, amongst other draws, the \code{S} draws from the posterior 
#' distribution of the \code{NxN} matrix of contemporaneous relationships \eqn{B}. 
#' These draws are to be normalised with respect to:
#' @param B_hat an \code{NxNxMxS5} array with \code{M} times \code{S5} \code{NxN} benchmark 
#' matrices specified by the user to be used to normalise the \eqn{B} matrices in the posterior output.
#' \code{S5} denotes the number of alternative S5 components.
#' If this argument is not provided the function creates its own proposal based on the 
#' last available state and component specific draws of \eqn{B} normalised to have 
#' positive elements on the main diagonal.
#' @param VB the list with matrices determining identification, including S5 identification
#' 
#' @return An object of class corresponding to the class of the first argument \code{posterior} 
#' with normalised draws of matrix \eqn{B}.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Waggoner, D.F., and Zha, T., (2003) Likelihood Preserving Normalization in Multiple Equation Models. 
#' \emph{Journal of Econometrics}, \bold{114}(2), 329--47, \doi{https://doi.org/10.1016/S0304-4076(03)00087-3}.
#'
#' @export
normalise <- function(posterior, B_hat, VB) {
  
  N = dim(B_hat)[1]
  
  # check the args
  stopifnot("Argument posterior must contain estimation output from the estimate function." = any(class(posterior)[1] == c("PosteriorBSVARSVMSS5", "PosteriorBSVARSVMS", "PosteriorBSVARSVS5")))
  # stopifnot("Argument B_hat must be a numeric array." = is.numeric(B_hat) & is.array(B_hat))
  stopifnot("Argument VB must be a list." = is.list(VB))
  
  # call method
  UseMethod("normalise", posterior)
}




#' @title Waggoner & Zha (2003) row signs normalisation of the posterior draws for matrix \eqn{B}
#'
#' @description Normalises the sign of rows of matrix \eqn{B} MCMC draws S5 component-by-component, 
#'  provided as the first argument \code{posterior_B}, relative to the matrices in
#'  \code{B_hat}, provided as the second argument of the function. If the second argument 
#'  is not provided, the function creates its own benchmark matrix. The implemented
#'  procedure proposed by Waggoner, Zha (2003) normalises the MCMC output in an
#'  optimal way leading to the unimodal posterior. Only normalised MCMC output is 
#'  suitable for the computations of the posterior characteristics of the \eqn{B}
#'  matrix elements and their functions such as the impulse response functions and other 
#'  economically interpretable values. 
#' 
#' @param posterior posterior estimation outcome - an object of class "PosteriorBSVARSVS5"
#' containing, amongst other draws, the \code{S} draws from the posterior 
#' distribution of the \code{NxN} matrix of contemporaneous relationships \eqn{B}. 
#' These draws are to be normalised with respect to:
#' @param B_hat an \code{NxNxS5} array with \code{S5} \code{NxN} benchmark 
#' matrices specified by the user to be used to normalise the \eqn{B} matrices in the posterior output.
#' \code{S5} denotes the number of alternative S5 components.
#' If this argument is not provided the function creates its own proposal based on the 
#' last available component specific draws of \eqn{B} normalised to have 
#' positive elements on the main diagonal.
#' @param VB the list with matrices determining identification, including S5 identification
#' 
#' @return An object of class corresponding to the class of the first argument \code{posterior} 
#' with normalised draws of matrix \eqn{B}.
#'
#' @method normalise PosteriorBSVARSVS5
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Waggoner, D.F., and Zha, T., (2003) Likelihood Preserving Normalization in Multiple Equation Models. 
#' \emph{Journal of Econometrics}, \bold{114}(2), 329--47, \doi{https://doi.org/10.1016/S0304-4076(03)00087-3}.
#'
#' @export
normalise.PosteriorBSVARSVS5 <- function(posterior, B_hat, VB) {
  
  S5_equation = 3
  
  M             = 1
  N             = dim(posterior$posterior$A)[1]
  comp          = VB[length(VB)][[1]]
  
  S5_indicator  = posterior$posterior$S4_indicator
  
  if ( missing(B_hat) ) {
    B_hat       = array(NA, c(N, N, comp[S5_equation]))
  }
  
  for (component in 1:comp[S5_equation]) {
    S5_indices    = which(S5_indicator[S5_equation, ] == component) 
    
    if ( length(S5_indices) == 0 ) next
    if ( missing(B_hat) ) {
      B_hat_tmp               = posterior$posterior$B[,,utils::tail(S5_indices,1)]
      B_hat[,,component]      = diag(sign(diag(B_hat_tmp))) %*% B_hat_tmp
    }
    B_to_normalise            = posterior$posterior$B[,,S5_indices]
    
    invisible(.Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, B_to_normalise, B_hat[,,component]))
    
    posterior$posterior$B[,,S5_indices] = B_to_normalise
    
    # last_draw
    if ( component == posterior$last_draw$S4_indicator[S5_equation,] ) {
      invisible(.Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, posterior$last_draw$B, B_hat[,,component]))
    }
  } # END component loop
  
return(posterior)
}



#' @title Waggoner & Zha (2003) row signs normalisation of the posterior draws for matrix \eqn{B}
#'
#' @description Normalises the sign of rows of matrix \eqn{B} MCMC draws S5 component-by-component, 
#'  provided as the first argument \code{posterior_B}, relative to the matrices in
#'  \code{B_hat}, provided as the second argument of the function. If the second argument 
#'  is not provided, the function creates its own benchmark matrix. The implemented
#'  procedure proposed by Waggoner, Zha (2003) normalises the MCMC output in an
#'  optimal way leading to the unimodal posterior. Only normalised MCMC output is 
#'  suitable for the computations of the posterior characteristics of the \eqn{B}
#'  matrix elements and their functions such as the impulse response functions and other 
#'  economically interpretable values. 
#' 
#' @param posterior posterior estimation outcome - an object of class "PosteriorBSVARSVMS"
#' containing, amongst other draws, the \code{S} draws from the posterior 
#' distribution of the \code{NxN} matrix of contemporaneous relationships \eqn{B}. 
#' These draws are to be normalised with respect to:
#' @param B_hat an \code{NxNxM} array with \code{M} \code{NxN} benchmark 
#' matrices specified by the user to be used to normalise the \eqn{B} matrices in the posterior output.
#' If this argument is not provided the function creates its own proposal based on the 
#' last available state-specific draws of \eqn{B} normalised to have 
#' positive elements on the main diagonal.
#' @param VB the list with matrices determining identification
#' 
#' @return An object of class corresponding to the class of the first argument \code{posterior} 
#' with normalised draws of matrix \eqn{B}.
#'
#' @method normalise PosteriorBSVARSVMS
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Waggoner, D.F., and Zha, T., (2003) Likelihood Preserving Normalization in Multiple Equation Models. 
#' \emph{Journal of Econometrics}, \bold{114}(2), 329--47, \doi{https://doi.org/10.1016/S0304-4076(03)00087-3}.
#'
#' @export
normalise.PosteriorBSVARSVMS <- function(posterior, B_hat, VB) {
  
  M             = dim(posterior$posterior$xi)[1]
  N             = dim(posterior$posterior$A)[1]
  S             = dim(posterior$posterior$A)[3]
  
  if ( missing(B_hat) ) {
    B_hat       = array(NA, c(N, N, M))
  }
  
  for (m in 1:M) {
    
    if ( missing(B_hat) ) {
      B_hat_tmp           = posterior$posterior$B[,,m]
      B_hat[,,m]          = diag(sign(diag(B_hat_tmp))) %*% B_hat_tmp
    }
    B_to_normalise        = array(NA, c(N, N, S))
    for (i in 1:S) {
      B_to_normalise[,,i]   = posterior$posterior$B[i,1][[1]][,,m]
    }
    
    invisible(.Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, B_to_normalise, B_hat[,,m]))
    
    for (i in 1:S) {
      posterior$posterior$B[i,1][[1]][,,m] = B_to_normalise[,,i]
    }
    
    # last_draw
    invisible(.Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, posterior$last_draw$B, B_hat[,,m]))
    
  } # END m loop
  
  return(posterior)
}



#' @title Waggoner & Zha (2003) row signs normalisation of the posterior draws for matrix \eqn{B}
#'
#' @description Normalises the sign of rows of matrix \eqn{B} MCMC draws Markov state-by-state
#' and S5 component-by-component relative to the matrices in
#'  \code{B_hat}, provided as the second argument of the function. If the second argument 
#'  is not provided, the function creates its own benchmark matrix. The implemented
#'  procedure proposed by Waggoner, Zha (2003) normalises the MCMC output in an
#'  optimal way leading to the unimodal posterior. Only normalised MCMC output is 
#'  suitable for the computations of the posterior characteristics of the \eqn{B}
#'  matrix elements and their functions such as the impulse response functions and other 
#'  economically interpretable values. 
#' 
#' @param posterior posterior estimation outcome - an object of class "PosteriorBSVARSVMSS5"
#' containing, amongst other draws, the \code{S} draws from the posterior 
#' distribution of the \code{NxN} matrix of contemporaneous relationships \eqn{B}. 
#' These draws are to be normalised with respect to:
#' @param B_hat an \code{NxNxMxS5} array with \code{M} times \code{S5} \code{NxN} benchmark 
#' matrices specified by the user to be used to normalise the \eqn{B} matrices in the posterior output.
#' \code{S5} denotes the number of alternative S5 components.
#' If this argument is not provided the function creates its own proposal based on the 
#' last available state and component specific draws of \eqn{B} normalised to have 
#' positive elements on the main diagonal.
#' @param VB the list with matrices determining identification, including S5 identification
#' 
#' @return An object of class corresponding to the class of the first argument \code{posterior} 
#' with normalised draws of matrix \eqn{B}.
#'
#' @method normalise PosteriorBSVARSVMSS5
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Waggoner, D.F., and Zha, T., (2003) Likelihood Preserving Normalization in Multiple Equation Models. 
#' \emph{Journal of Econometrics}, \bold{114}(2), 329--47, \doi{https://doi.org/10.1016/S0304-4076(03)00087-3}.
#'
#' @export
normalise.PosteriorBSVARSVMSS5 <- function(posterior, B_hat, VB) {
  
  S5_equation = 3
  
  M             = dim(posterior$posterior$xi)[1]
  N             = dim(posterior$posterior$A)[1]
  comp          = VB[length(VB)][[1]]
  S5_indicator  = posterior$posterior$S4_indicator
  
  if ( missing(B_hat) ) {
    B_hat       = array(NA, c(N, N, M, comp[S5_equation]))
  }
  
  for (m in 1:M) {
    for (component in 1:comp[S5_equation]) {
      
      S5_indices    = which(S5_indicator[S5_equation, m, ] == component) #
      if ( length(S5_indices) == 0 ) next
      if ( missing(B_hat) ) {
        B_hat_tmp               = posterior$posterior$B[utils::tail(S5_indices,1),1][[1]][,,m] #
        B_hat[,,m,component]    = diag(sign(diag(B_hat_tmp))) %*% B_hat_tmp
      }
      B_to_normalise = array(NA, c(N, N, length(S5_indices)))
      for (i in 1:length(S5_indices)) {
        B_to_normalise[,,i]   = posterior$posterior$B[S5_indices[i],1][[1]][,,m] # this is ridiculus!
      }
      
      invisible(.Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, B_to_normalise, B_hat[,,m,component]))
      
      for (i in 1:length(S5_indices)) {
        posterior$posterior$B[S5_indices[i],1][[1]][,,m] = B_to_normalise[,,i]
      }
      
      # last_draw
      if ( component == posterior$last_draw$S4_indicator[S5_equation,m] ) {
        invisible(.Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, posterior$last_draw$B, B_hat[,,m,component]))
      }
      
    } # END component loop
  } # END m loop
  
  return(posterior)
}

