

#' @title Drops Stochastic Volatility latent variables from the posterior output 
#' and returns draws on the remaining parameters
#'
#' @description Drops list elements \code{h} and \code{sigma} from the posterior 
#' output and returns draws on the remaining parameters. To be used for reducing 
#' the memory used to save the posterior output.
#' 
#' 
#' @param posterior an object of class of class PosteriorBSVARSVMSS5, 
#' PosteriorBSVARSVMS, or PosteriorBSVARSVS5 with random draws from the posterior 
#' distribution of the Structural model with Markov-switching structural matrix
#' 
#' @return An object of the same class as the input object class but with an extra 
#' attribute \code{"drop_sv"}
#' without the Stochastic Volatility latent processes draws.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @seealso \code{\link{drag_sv}}, \code{\link{djoin_sv}}
#' 
#' @export
drop_sv <- function(posterior) {
  
  # check arguments
  stopifnot("Argument posterior must be of class PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5." = substr(class(posterior), 1, 16) == "PosteriorBSVARSV")
  
  out = posterior$posterior
  out = within(out, rm("sigma", "h"))
  posterior$posterior = out
  class(posterior) = c(class(posterior), "drop_sv")
  
  return(posterior)
}


#' @title Returns only Stochastic Volatility latent variables from 
#' the posterior output
#' 
#' @description Returns a list with two elements, \code{h} and \code{sigma}, from the posterior 
#' output. To be used for managing the memory used to save the posterior output.
#' 
#' @param posterior an object of class of class PosteriorBSVARSVMSS5, 
#' PosteriorBSVARSVMS, or PosteriorBSVARSVS5 with random draws from the posterior 
#' distribution of the Structural model with Markov-switching structural matrix
#' 
#' @return An object of the same class as the input object class but with an extra 
#' attribute \code{"drag_sv"} - a list with two elements, \code{h} and \code{sigma}, 
#' from the posterior output.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @seealso \code{\link{drop_sv}}, \code{\link{djoin_sv}}
#' 
#' @export
drag_sv <- function(posterior) {
  
  # check arguments
  stopifnot("Argument posterior must be of class PosteriorBSVARSVMSS5, PosteriorBSVARSVMS, or PosteriorBSVARSVS5." = substr(class(posterior), 1, 16) == "PosteriorBSVARSV")
  
  out           = list(
    posterior   = list(
      h         = posterior$posterior$h,
      sigma     = posterior$posterior$sigma
    )
  )
  class(out) = c(class(posterior), "drag_sv")
  
  return(out)
}


#' @title Merges the posterior output of the Stochastic Volatility part with the rest of the parameters
#' 
#' @description Merges the posterior output of the Stochastic Volatility part created using 
#' function \code{drag_sv()} with the rest of the posterior draws created using function
#' \code{drop_sv()}.
#' 
#' @param drop_sv an object of class \code{drop_sv}
#' @param drag_sv an object of class \code{drag_sv}
#' 
#' @return Posterior draws from the estimation procedure.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @seealso \code{\link{drop_sv}}, \code{\link{drag_sv}}
#' 
#' @export
djoin_sv <- function(drop_sv, drag_sv) {
  
  # check arguments
  stopifnot("Argument drop_sv must be of class drop_sv." = any(class(drop_sv) == "drop_sv"))
  stopifnot("Argument drag_sv must be of class drag_sv." = any(class(drag_sv) == "drag_sv"))
  stopifnot("Both arguments must come from the estimation of the same model." = class(drag_sv)[1] == class(drop_sv)[1])
  
  drop_sv$posterior$h     = drag_sv$posterior$h
  drop_sv$posterior$sigma = drag_sv$posterior$sigma
  class(drop_sv) = class(drop_sv)[1]
  
  return(drop_sv)
}
