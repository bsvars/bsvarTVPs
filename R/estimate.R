

#' @title Bayesian estimation of Structural Vector Autoregressions with S4 and 
#' Markov-Switching in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and Markov-switching and Stochastic Search Specification Selection 
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
#' @return A list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_s4_sv_cpp()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_mss_s4_sv_cpp <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_mss_s4_sv_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSV-MSS5"
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
#' @return A list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_s4_sv_cpp()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_mss_sv_cpp <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_mss_sv_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSV-MS"
  return(output)
}






#' @title Bayesian estimation of Structural Vector Autoregressions with S4 in the Structural Matrix
#'
#' @description Estimates a Structural Vector Autoregressions with Stochastic Volatility
#' heteroskedasticity and Stochastic Search Specification Selection for the structural matrix
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
#' @return A list containing the Bayesian estimation output in two elements:
#' 
#' \code{posterior} a list with a collection of \code{S} draws from the posterior distribution 
#' generated via Gibbs sampler containing many arrays and vectors whose selection depends on 
#' the model specification.
#' \code{last_draw} a list with the last draw of the current MCMC run as the starting value 
#' to be passed to the continuation of the MCMC estimation using \code{bsvar_s4_sv_cpp()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in 
#' Monetary Policy Shock Identification?
#' 
#' @export
bsvar_s4_sv_cpp <- function(SS, Y, X, prior, VB, starting_values, thin = 100L) {
  output          = .Call(`_bsvarTVPs_bsvar_s4_sv_cpp`, SS, Y, X, prior, VB, starting_values, thin)
  class(output)   = "PosteriorBSVARSV-S5"
  return(output)
}
