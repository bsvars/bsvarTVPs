
#' @title Bayesian forecasting of Structural Vector Autoregressions with TVI and 
#' Markov-Switching in the Structural Matrix (and Autoregressive matrix)
#'
#' @description Forecasts a Structural Vector Autoregressions with Stochastic 
#' Volatility heteroskedasticity and Markov-switching and time-varying 
#' identification in the structural matrix and a three-level local-global prior 
#' hierarchical prior structure for the structural and autoregressive matrices
#' 
#' @param posterior a list containing Bayesian estimation output
#' @param X a \code{KxT} matrix with regressors
#' @param horizon a positive integer specifying the forecast horizon
#' @param non_explosive a logical indicating whether only non explosive processes
#' should be used for forecasting
#' 
#' @return An object of class \code{Forecasts} - a list containing the Bayesian 
#' forecasting output with three elements:
#' \code{forecast} an \code{N x horizon x S} array with draws from predictive density
#' \code{forecast_mean} an \code{N x horizon x S} array with the mean of the
#' predictive density
#' \code{forecast_cov} an \code{N x N x horizon x S} array with the covariance 
#' of the predictive density
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2025) Time-Varying Identification of Structural Vector Autoregressions
#' 
#' @export
forecast_bsvar_mssa_sv <- function(posterior, X, horizon = 1L, non_explosive = FALSE) {
  
  T               = dim(posterior$posterior$h)[2]
  N               = dim(posterior$posterior$h)[1]
  S               = dim(posterior$posterior$h)[3]
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  posterior_PR_TR = posterior$posterior$PR_TR
  posterior_xi_T  = posterior$posterior$xi[,T,]
  posterior_h_T   = posterior$posterior$h[,T,]
  posterior_rho   = posterior$posterior$rho
  posterior_omega = posterior$posterior$omega
  X_T             = X[,T]
  
  output          = .Call(`_bsvarTVPs_forecast_mssa_sv`, posterior_B, posterior_A, posterior_PR_TR, posterior_xi_T, posterior_h_T, posterior_rho, posterior_omega, X_T, horizon, non_explosive)
  
  SS                  = dim(output$forecast)[3]
  forecast_covariance = array(NA, c(N, N, horizon, SS))
  for (s in 1:SS) forecast_covariance[,,,s] = output$forecast_cov[s,][[1]]
  output$forecast_covariance = forecast_covariance
  
  class(output)   = "Forecasts"
  return(output)
}


#' @title Bayesian forecasting of Structural Vector Autoregressions with TVI and 
#' Markov-Switching in the Structural Matrix
#'
#' @description Forecasts a Structural Vector Autoregressions with Stochastic 
#' Volatility heteroskedasticity and Markov-switching and time-varying 
#' identification in the structural matrix and a three-level local-global prior 
#' hierarchical prior structure for the structural and autoregressive matrices
#' 
#' @param posterior a list containing Bayesian estimation output
#' @param X a \code{KxT} matrix with regressors
#' @param horizon a positive integer specifying the forecast horizon
#' @param non_explosive a logical indicating whether only non explosive processes
#' should be used for forecasting
#' 
#' @return An object of class \code{Forecasts} - a list containing the Bayesian 
#' forecasting output with three elements:
#' \code{forecast} an \code{N x horizon x S} array with draws from predictive density
#' \code{forecast_mean} an \code{N x horizon x S} array with the mean of the
#' predictive density
#' \code{forecast_cov} an \code{N x N x horizon x S} array with the covariance 
#' of the predictive density
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2025) Time-Varying Identification of Structural Vector Autoregressions
#' 
#' @export
forecast_bsvar_mss_sv <- function(posterior, X, horizon = 1L, non_explosive = FALSE) {
  
  T               = dim(posterior$posterior$h)[2]
  N               = dim(posterior$posterior$h)[1]
  S               = dim(posterior$posterior$h)[3]
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  posterior_PR_TR = posterior$posterior$PR_TR
  posterior_xi_T  = posterior$posterior$xi[,T,]
  posterior_h_T   = posterior$posterior$h[,T,]
  posterior_rho   = posterior$posterior$rho
  posterior_omega = posterior$posterior$omega
  X_T             = X[,T]
  
  output          = .Call(`_bsvarTVPs_forecast_mss_sv`, posterior_B, posterior_A, posterior_PR_TR, posterior_xi_T, posterior_h_T, posterior_rho, posterior_omega, X_T, horizon, non_explosive)
  
  SS                  = dim(output$forecast)[3]
  forecast_covariance = array(NA, c(N, N, horizon, SS))
  for (s in 1:SS) forecast_covariance[,,,s] = output$forecast_cov[s,][[1]]
  output$forecast_covariance = forecast_covariance
  
  class(output)   = "Forecasts"
  return(output)
}




#' @title Bayesian forecasting of Structural Vector Autoregressions with TVI and 
#' Markov-Switching in the Structural Matrix, and homoskedastic structural shocks
#'
#' @description Forecasts a Structural Vector Autoregressions with 
#' Markov-switching and time-varying identification in the structural matrix,
#' homoskedastic structural shocks, and a three-level local-global prior 
#' hierarchical prior structure for the structural and autoregressive matrices
#' 
#' @param posterior a list containing Bayesian estimation output
#' @param X a \code{KxT} matrix with regressors
#' @param horizon a positive integer specifying the forecast horizon
#' @param non_explosive a logical indicating whether only non explosive processes
#' should be used for forecasting
#' 
#' @return An object of class \code{Forecasts} - a list containing the Bayesian 
#' forecasting output with three elements:
#' \code{forecast} an \code{N x horizon x S} array with draws from predictive density
#' \code{forecast_mean} an \code{N x horizon x S} array with the mean of the
#' predictive density
#' \code{forecast_cov} an \code{N x N x horizon x S} array with the covariance 
#' of the predictive density
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2025) Time-Varying Identification of Structural Vector Autoregressions
#' 
#' @export
forecast_bsvar_mss <- function(posterior, X, horizon = 1L, non_explosive = FALSE) {
  
  T               = dim(posterior$posterior$xi)[2]
  N               = dim(posterior$posterior$A)[1]
  S               = dim(posterior$posterior$A)[3]
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  posterior_PR_TR = posterior$posterior$PR_TR
  posterior_xi_T  = posterior$posterior$xi[,T,]
  X_T             = X[,T]
  
  output          = .Call(`_bsvarTVPs_forecast_mss`, posterior_B, posterior_A, posterior_PR_TR, posterior_xi_T, X_T, horizon, non_explosive)
  
  SS                  = dim(output$forecast)[3]
  forecast_covariance = array(NA, c(N, N, horizon, SS))
  for (s in 1:SS) forecast_covariance[,,,s] = output$forecast_cov[s,][[1]]
  output$forecast_covariance = forecast_covariance
  
  class(output)   = "Forecasts"
  return(output)
}



#' @title Bayesian forecasting of Structural Vector Autoregressions with TVI
#'
#' @description Forecasts a Structural Vector Autoregressions with Stochastic 
#' Volatility heteroskedasticity and time-varying identification in the 
#' structural matrix and a three-level local-global prior hierarchical prior 
#' structure for the structural and autoregressive matrices
#' 
#' @param posterior a list containing Bayesian estimation output
#' @param X a \code{KxT} matrix with regressors
#' @param horizon a positive integer specifying the forecast horizon
#' @param non_explosive a logical indicating whether only non explosive processes
#' should be used for forecasting
#' 
#' @return An object of class \code{Forecasts} - a list containing the Bayesian 
#' forecasting output with three elements:
#' \code{forecast} an \code{N x horizon x S} array with draws from predictive density
#' \code{forecast_mean} an \code{N x horizon x S} array with the mean of the
#' predictive density
#' \code{forecast_covariance} an \code{N x N x horizon x S} array with the covariance 
#' of the predictive density
#' \code{forecast_cov} an \code{S}-field of \code{N x N x horizon} arrays with 
#' the covariance of the predictive density
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2025) Time-Varying Identification of Structural Vector Autoregressions
#' 
#' @export
forecast_bsvar_sv <- function(posterior, X, horizon = 1L, non_explosive = FALSE) {
  
  T               = dim(posterior$posterior$h)[2]
  N               = dim(posterior$posterior$h)[1]
  S               = dim(posterior$posterior$h)[3]
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  posterior_h_T   = posterior$posterior$h[,T,]
  posterior_rho   = posterior$posterior$rho
  posterior_omega = posterior$posterior$omega
  X_T             = X[,T]
  
  output          = .Call(`_bsvarTVPs_forecast_sv`, posterior_B, posterior_A, posterior_h_T, posterior_rho, posterior_omega, X_T, horizon, non_explosive)
  
  SS                  = dim(output$forecast)[3]
  forecast_covariance = array(NA, c(N, N, horizon, SS))
  for (s in 1:SS) forecast_covariance[,,,s] = output$forecast_cov[s,][[1]]
  output$forecast_covariance = forecast_covariance
  
  class(output)   = "Forecasts"
  return(output)
}