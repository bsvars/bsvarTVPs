
#' @title Bayesian estimation of Structural Vector Autoregressions with TVI and 
#' Markov-Switching in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and Markov-switching and time-varying identification 
#' in the structural matrix
#' and a three-level local-global prior hierarchical prior structure for the structural and autoregressive matrices
#' 
#' @param SS a positive integer, the number of posterior draws to be generated
#' @param Y a \code{NxT} matrix with dependent variables
#' @param X a \code{KxT} matrix with regressors
#' @param prior a list containing the prior specification
#' @param VB a list providing the structural matrix specification
#' @param starting_values a list providing starting values to the estimated parameters
#' @param thin a positive integer determining MCMC thinning
#' @param centred_sv a logical value indicating whether the SV model should be 
#' in its centred form
#' @param finiteM a logical value, if \code{TRUE} a stationary Markov-switching 
#' model is estimated, if \code{FALSE} an over-fitting Markov-switching model is used
#' @param hyper_select an integer choosing th type of hyper-parameter hierarchy: 
#' \code{1} - horseshoe prior, \code{2} - 3-level hierarchy, \code{3} - fixed hyper-parameters.
#' @param studentt a logical value, if \code{TRUE} the model is estimated with 
#' t-distributed structural shocks, if \code{FALSE} the shocks are Gaussian
#' 
#' @return An object of class \code{PosteriorBSVARSVMSTVI} - a list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_mss_tvi_sv_boost <- function(SS, Y, X, prior, VB, starting_values, thin = 100L, centred_sv = FALSE, finiteM = TRUE, hyper_select = 1, studentt = TRUE) {
  output          = .Call(`_bsvarTVPs_bsvar_mss_s4_sv_boost_cpp`, SS, Y, X, prior, VB, starting_values, thin, centred_sv, finiteM, hyper_select, studentt)
  class(output)   = "PosteriorBSVARSVMSTVI"
  return(output)
}



#' @title Bayesian estimation of Structural Vector Autoregressions with
#' Markov-Switching in the Structural and Autoregressive Matrices and TVI
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and Markov-switching in the structural and autoregressive 
#' matrices and with time-varying identification, and a three-level local-global 
#' prior hierarchical prior structure for the structural and autoregressive matrices
#' 
#' @param SS a positive integer, the number of posterior draws to be generated
#' @param Y a \code{NxT} matrix with dependent variables
#' @param X a \code{KxT} matrix with regressors
#' @param prior a list containing the prior specification
#' @param VB a list providing the structural matrix specification
#' @param starting_values a list providing starting values to the estimated parameters
#' @param thin a positive integer determining MCMC thinning
#' @param centred_sv a logical value indicating whether the SV model should be 
#' in its centred form
#' @param finiteM a logical value, if \code{TRUE} a stationary Markov-switching 
#' model is estimated, if \code{FALSE} an over-fitting Markov-switching model is used
#' @param hyper_select an integer choosing th type of hyper-parameter hierarchy: 
#' \code{1} - horseshoe prior, \code{2} - 3-level hierarchy, \code{3} - fixed hyper-parameters.
#' @param studentt a logical value, if \code{TRUE} the model is estimated with 
#' t-distributed structural shocks, if \code{FALSE} the shocks are Gaussian
#' 
#' @return An object of class \code{PosteriorBSVARSVMSATVI} - a list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_mssa_tvi_sv_boost <- function(SS, Y, X, prior, VB, starting_values, thin = 100L, centred_sv = FALSE, finiteM = TRUE, hyper_select = 1, studentt = TRUE) {
  output          = .Call(`_bsvarTVPs_bsvar_mssa_s4_sv_boost_cpp`, SS, Y, X, prior, VB, starting_values, thin, centred_sv, finiteM, hyper_select, studentt)
  class(output)   = "PosteriorBSVARSVMSATVI"
  return(output)
}

