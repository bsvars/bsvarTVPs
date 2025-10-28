
#' @export
generics::forecast

#' @title Bayesian forecasting of Structural Vector Autoregressions with TVI and 
#' Markov-Switching in the Structural Matrix (and Autoregressive matrix)
#'
#' @description Forecasts a Structural Vector Autoregressions with Stochastic 
#' Volatility heteroskedasticity and Markov-switching and time-varying 
#' identification in the structural matrix and a three-level local-global prior 
#' hierarchical prior structure for the structural and autoregressive matrices
#' 
#' @method forecast PosteriorBSVARTVP
#' 
#' @param object posterior estimation outcome of class \code{PosteriorBSVARTVP}
#' obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast forecasted values of the exogenous variables.
#' @param ... Not used.
#' 
#' @return A list of class \code{Forecasts} containing the
#' draws from the predictive density and for heteroskedastic models the draws 
#' from the predictive density of structural shocks conditional standard 
#' deviations and data. The output elements include:
#' \describe{
#'  \item{forecasts}{an \code{N x horizon x S} array with the draws from predictive density}
#'  \item{forecast_mean}{an \code{N x horizon x S} array with the mean of the predictive density}
#'  \item{forecast_cov}{an \code{N x N x horizon x S} array with the covariance of the predictive density}
#'  \item{Y}{an \eqn{NxT} matrix with the data on dependent variables}
#' }
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' Camehl, A. & Woźniak, T. (2025) Time-Varying Identification of Structural Vector Autoregressions, <doi:10.48550/arXiv.2502.19659>.
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec   = specify_bsvarTVP$new(us_fiscal_lsuw)    # specify the model
#' burn   = estimate(spec, 5)                       # run the burn-in for convergence
#' post   = estimate(burn, 10)                      # estimate the model
#' fore   = forecast(post, horizon = 2)             # forecast 2 periods ahead
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_fiscal_lsuw |>
#'   specify_bsvarTVP$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) -> post
#' post |> forecast(horizon = 2) -> fore
#' 
#' @export
forecast.PosteriorBSVARTVP <- function(
    object, 
    horizon = 1L, 
    exogenous_forecast = NULL,
    ...
) {
  
  T               = dim(object$posterior$h)[2]
  N               = dim(object$posterior$h)[1]
  S               = dim(object$posterior$h)[3]
  X_T             = object$last_draw$data_matrices$X[,T]
  K               = length(X_T)
  d               = K - N * object$last_draw$get_p() - 1
  posterior_B     = object$posterior$B_cpp
  posterior_A     = object$posterior$A_cpp
  posterior_PR_TR = object$posterior$PR_TR
  posterior_xi_T  = matrix(object$posterior$xi[,T,], ncol = S)
  posterior_h_T   = object$posterior$h[,T,]
  posterior_rho   = object$posterior$rho
  posterior_omega = object$posterior$omega
  posterior_df    = object$posterior$df
  
  # prepare forecasting with exogenous variables
  if (d == 0 ) {
    exogenous_forecast = matrix(NA, horizon, 1)
  } else {
    stopifnot("Forecasted values of exogenous variables are missing." = (d > 0) & !is.null(exogenous_forecast))
    stopifnot("The matrix of exogenous_forecast does not have a correct number of columns." = ncol(exogenous_forecast) == d)
    stopifnot("Provide exogenous_forecast for all forecast periods specified by argument horizon." = nrow(exogenous_forecast) == horizon)
    stopifnot("Argument exogenous has to be a matrix." = is.matrix(exogenous_forecast) & is.numeric(exogenous_forecast))
    stopifnot("Argument exogenous cannot include missing values." = sum(is.na(exogenous_forecast)) == 0 )
  }
  
  if (object$last_draw$get_sv()) {
    sv_select     = 1
  } else {
    sv_select     = 3
  }
  studentt        = !object$last_draw$get_normal()
  
  if (object$last_draw$get_msa()) {
    output        = .Call(`_bsvarTVPs_forecast_mssa_sv`, 
                          posterior_B, posterior_A, posterior_PR_TR, posterior_xi_T, posterior_h_T, posterior_rho, posterior_omega, posterior_df, 
                          X_T, exogenous_forecast, horizon, sv_select, studentt)
  } else {
    output        = .Call(`_bsvarTVPs_forecast_mss_sv`, 
                          posterior_B, posterior_A, posterior_PR_TR, posterior_xi_T, posterior_h_T, posterior_rho, posterior_omega, posterior_df, 
                          X_T, exogenous_forecast, horizon, sv_select, studentt)
  }  
  
  SS              = dim(output$forecast)[3]
  forecast_covariance         = array(NA, c(N, N, horizon, SS))
  for (s in 1:SS) forecast_covariance[,,,s] = output$forecast_cov[s,][[1]]
  output$forecast_covariance  = forecast_covariance
  
  output$Y        = object$last_draw$data_matrices$Y
  
  class(output)   = "Forecasts"
  return(output)
}
