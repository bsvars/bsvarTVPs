

#' @title Bayesian estimation of Structural Vector Autoregressions with TVI and 
#' Markov-Switching in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and Markov-switching and time-varying identification 
#' in the structural matrix
#' 
#' 
#' @param SS a positive integer, the number of posterior draws to be generated
#' @param Y a \code{NxT} matrix with dependent variables
#' @param X a \code{KxT} matrix with regressors
#' @param prior a list containing the prior specification
#' @param VB a list providing the structural matrix specification
#' @param starting_values a list providing starting values to the estimated parameters
#' @param thin a positive integer determining MCMC thinning
#' 
#' @return An object of class \code{PosteriorBSVARSVMSTVI} - a list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_tvi_sv()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_mss_tvi_sv <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_mss_s4_sv_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSVMSTVI"
  return(output)
}





#' @title Bayesian estimation of Structural Vector Autoregressions with Markov-Switching 
#' in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and Markov-switching for the structural matrix
#' 
#' 
#' @param SS a positive integer, the number of posterior draws to be generated
#' @param Y a \code{NxT} matrix with dependent variables
#' @param X a \code{KxT} matrix with regressors
#' @param prior a list containing the prior specification
#' @param VB a list providing the structural matrix specification
#' @param starting_values a list providing starting values to the estimated parameters
#' @param thin a positive integer determining MCMC thinning
#' 
#' @return  An object of class \code{PosteriorBSVARSVMS} - a list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_tvi_sv()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_mss_sv <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_mss_sv_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSVMS"
  return(output)
}






#' @title Bayesian estimation of Structural Vector Autoregressions with TVI in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and time-varying identification for the structural matrix
#' 
#' 
#' @param SS a positive integer, the number of posterior draws to be generated
#' @param Y a \code{NxT} matrix with dependent variables
#' @param X a \code{KxT} matrix with regressors
#' @param prior a list containing the prior specification
#' @param VB a list providing the structural matrix specification
#' @param starting_values a list providing starting values to the estimated parameters
#' @param thin a positive integer determining MCMC thinning
#' 
#' @return  An object of class \code{PosteriorBSVARSVTVI} - a list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_tvi_sv()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_tvi_sv <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_s4_sv_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSVTVI"
  return(output)
}





#' @title Bayesian estimation of Structural Vector Autoregressions with TVI in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and time-varying identification for the structural matrix 
#' and a three-level local-global prior hierarchical prior structure for the structural and autoregressive matrices
#' 
#' 
#' @param SS a positive integer, the number of posterior draws to be generated
#' @param Y a \code{NxT} matrix with dependent variables
#' @param X a \code{KxT} matrix with regressors
#' @param prior a list containing the prior specification
#' @param VB a list providing the structural matrix specification
#' @param starting_values a list providing starting values to the estimated parameters
#' @param thin a positive integer determining MCMC thinning
#' 
#' @return  An object of class \code{PosteriorBSVARSVTVI} - a list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_tvi_sv_boost()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_tvi_sv_boost <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_s4_sv_boost_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSVTVIBOOST"
  return(output)
}
