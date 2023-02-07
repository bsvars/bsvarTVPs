
#' @title Computes posterior draws of regime probabilities
#'
#' @description Each of the draws from the posterior estimation of a model is transformed into
#' a draw from the posterior distribution of the regime probabilities. These represent either
#' the realisations of the regime indicators, when \code{type = "realized"}, filtered probabilities,
#' when \code{type = "filtered"}, forecasted regime probabilities, when \code{type = "forecasted"},
#' or the smoothed probabilities, when \code{type = "smoothed"}, .
#' 
#' @param posterior posterior estimation outcome of regime-dependent heteroskedastic models 
#' - an object of either of the classes: PosteriorBSVARSV-MSS5, or PosteriorBSVARSV-MS.
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
  
  stopifnot("Argument posterior must be of class PosteriorBSVARSV-MSS5 or PosteriorBSVARSV-MS." = substr(class(posterior), 1, 19) == "PosteriorBSVARSV-MS")
  
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

