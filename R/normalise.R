
#' @title Waggoner & Zha (2003) row signs normalisation of the posterior draws 
#' for the structural matrix \eqn{B}
#'
#' @description Normalises the sign of rows of matrix \eqn{B} MCMC draws, 
#'  relative to matrix \code{B_benchmark}, provided as the second argument. The implemented
#'  procedure proposed by Waggoner, Zha (2003) normalises the MCMC output in an
#'  optimal way leading to the unimodal posterior. Only normalised MCMC output is 
#'  suitable for the computations of the posterior characteristics of the \eqn{B}
#'  matrix elements and their functions such as the impulse response functions and other 
#'  economically interpretable values. 
#' 
#' @param posterior posterior estimation outcome of class \code{PosteriorBSVARTVP} 
#' generated using the \code{estimate()} function, amongst other draws, 
#' the \code{S} draws from the posterior distribution of the \code{NxN} 
#' structural matrix of contemporaneous relationships \eqn{B}. These draws are 
#' to be normalised with respect to the matrix \code{B_benchmark}.
#' @param B_benchmark an \code{NxNxMxC} array containing the benchmark structural 
#' matrices specified by the user to have the desired row signs, where \code{C}
#' is the maximum number of identification for a row in the specified system. 
#' If not provided, it is set to the last draws of the structural matrix for a
#' particular regime and component with the row signs ensuring that the diagonal 
#' elements are positive.
#' 
#' @return An object of the same class as that provided as the input argument 
#' \code{posterior} containing the posterior draws including the draws of the 
#' normalised structural matrix.
#' 
#' @method normalise PosteriorBSVARTVP
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Waggoner, D.F., and Zha, T., (2003) Likelihood Preserving Normalization in Multiple Equation Models. 
#' \emph{Journal of Econometrics}, \bold{114}(2), 329--47, \doi{10.1016/S0304-4076(03)00087-3}.
#'
#' @examples
#' spec    = specify_bsvarTVP$new(us_fiscal_lsuw)  # specify the model
#' burn    = estimate(spec, 5)                     # run the burn-in
#' post    = estimate(burn, 5)                     # estimate the model
#' post    = normalise(post)                       # normalise the posterior
#' 
#' @export
normalise.PosteriorBSVARTVP <- function(posterior, B_benchmark = NULL) {
  
  VB            = posterior$last_draw$identification$VB
  TVI_equation  = which.max(VB[length(VB)][[1]])
  comp          = VB[length(VB)][[1]]
  
  M             = dim(posterior$posterior$xi)[1]
  N             = length(comp)
  
  set_B_benchmark = is.null(B_benchmark) 
  if (set_B_benchmark) {
    B_benchmark         = array(NA, c(N, N, M, comp[TVI_equation]))
  }
  
  TVI_indicator = posterior$posterior$S4_indicator
  
  for (m in 1:M) {
    for (component in 1:comp[TVI_equation]) {
      
      TVI_indices    = which(TVI_indicator[TVI_equation, m, ] == component) #
      if ( length(TVI_indices) == 0 ) next
      
      if (set_B_benchmark) {
        B_benchmark_tmp               = posterior$posterior$B_cpp[utils::tail(TVI_indices,1),1][[1]][,,m] #
        B_benchmark[,,m,component]    = diag(sign(diag(B_benchmark_tmp))) %*% B_benchmark_tmp
      }
      
      B_to_normalise = array(NA, c(N, N, length(TVI_indices)))
      for (i in 1:length(TVI_indices)) {
        B_to_normalise[,,i]   = posterior$posterior$B_cpp[TVI_indices[i],1][[1]][,,m] # this is ridiculous!
      }
      
      B_to_normalise          = .Call(`_bsvarTVPs_bsvars_normalisation_wz2003`, B_to_normalise, B_benchmark[,,m,component])
      
      for (i in 1:length(TVI_indices)) {
        posterior$posterior$B_cpp[TVI_indices[i],1][[1]][,,m] = B_to_normalise[,,i]
        posterior$posterior$B[,,m,TVI_indices[i]]             = B_to_normalise[,,i]
      }
      
      # last_draw
      if ( component == posterior$last_draw$starting_values$S4_indicator[TVI_equation,m] ) {
        posterior$last_draw$starting_values$B[,,m]    = .Call(`_bsvarTVPs_bsvars_normalisation_wz20031`, posterior$last_draw$starting_values$B[,,m], B_benchmark[,,m,component])
      }
      
    } # END component loop
  } # END m loop
  
  posterior$set_normalised()
  
  return(posterior)
} # normalise.PosteriorBSVARTVP

