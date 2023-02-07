
#' @title Plots the median and highest posterior density region for a sequence of \code{K} random variables
#'
#' @description Plots the median and highest posterior density regions for a sequence 
#' of \code{K} random variables based on the \code{S} posterior draws provided for each of them. 
#' It is useful in plotting impulse responses or conditional standard deviations.
#' 
#' @param draws a \code{K x S} matrix with \code{S} posterior draws of 
#' \code{K} random variables, or a \code{K x S x N} array with \code{N} such matrices
#' @param probability a number from interval \code{(0,1)} denoting the probability content 
#' of the highest posterior density regions
#' @param col a colour of the plot
#' @param ylim the range of the \code{y} axis
#' @param ylab the label of the \code{y} axis
#' @param xlab the label of the \code{x} axis
#' @param start_at an integer to denote the beginning of the \code{x} axis range
#' @param add a logical value. If \code{TRUE} the current ribbon plot is added to an existing one
#' @param ... other graphical parameters to be passed to \code{base:plot}
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @export
ribbon_plot = function(
    draws,              # K x S x N
    probability = 0.9,  # integer
    col         = 1,
    ylim,
    ylab,
    xlab,
    start_at    = 0,
    add         = FALSE,
    ...
) {
  
  stopifnot("Argument draws must be a matrix or an array." = any(class(draws) == "matrix") || any(class(draws) == "array"))
  stopifnot("Argument probability must be a number from interval (0,1)." = is.numeric(probability) & length(probability) == 1 & probability > 0 & probability < 1)
  stopifnot("Argument start_at must be a matrix or an integer number." = start_at %% 1 == 0)
  stopifnot("Argument add must be a logical value." = is.logical(add) & length(add) == 1)
  
  mat_or_array_tmp  = length(dim(draws))
  
  K               = dim(draws)[1]
  S               = dim(draws)[2]
  if ( mat_or_array_tmp == 2 ) {
    N               = 1
    draws_tmp       = array(NA, c(K, S, 1)) 
    draws_tmp[,,1]  = draws
    draws           = draws_tmp
  } else if ( mat_or_array_tmp == 3 ) {
    N               = dim(draws)[3]
  }
  # compute characteristics of draws
  draws_median    = apply(draws, c(1,3), stats::median) # K x N
  draws_hdi       = apply(draws, c(1,3), HDInterval::hdi, credMass = probability) # 2 x K x N
  draws_range     = range(draws_hdi)
  
  # create col_ribbon
  col_ribbon_rgb  = grDevices::col2rgb(col)
  col_ribbon      = grDevices::rgb(col_ribbon_rgb[1], col_ribbon_rgb[2], col_ribbon_rgb[3], 100, maxColorValue = 255)
  
  # manage the arguments
  if ( missing(ylim) ) ylim = draws_range
  if ( missing(ylab) ) ylab = ""
  if ( missing(xlab) ) xlab = ""
  
  if ( !add ) {
    plot(
      x      = start_at:(K - 1 + start_at), 
      y      = draws_median[,1],
      type   = "n",
      ylim   = ylim,
      ylab   = ylab,
      xlab   = xlab,
      ...
    )
  }
  for (n in 1:N) {
    graphics::polygon(
      x      = c(start_at:(K - 1 + start_at), (K - 1 + start_at):start_at),
      y      = c(draws_hdi[1, 1:K, n], draws_hdi[2, K:1, n]),
      col    = col_ribbon,
      border = col_ribbon
    )
    graphics::lines(
      x      = start_at:(K - 1 + start_at), 
      y      = draws_median[,n],
      type   = "l",
      lwd    = 2,
      col    = col
    )
  }
}
