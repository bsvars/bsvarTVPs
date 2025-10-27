
#' R6 Class Representing PriorBSVARTVP
#'
#' @description
#' The class PriorBSVARTVP presents a prior specification for the bsvarTVP model.
#' 
#' @examples 
#' # a prior for 3-variable example with one lag
#' prior = specify_prior_bsvarTVP$new(us_fiscal_lsuw[1:10,], N = 3, M = 2, p = 1)
#' prior$A # show autoregressive prior mean
#' 
#' @export
specify_prior_bsvarTVP = R6::R6Class(
  "PriorBSVARTVP",
  
  public = list(
    
    #' @field A an \code{NxK} matrix, the mean of the normal prior distribution for the parameter matrix \eqn{A}. 
    A          = matrix(),
    
    #' @field A_V_inv a \code{KxK} precision matrix of the normal prior distribution for each of 
    #' the row of the parameter matrix \eqn{A}. This precision matrix is equation invariant.
    A_V_inv    = matrix(),
    
    #' @field B_V_inv an \code{NxN} precision matrix of the generalised-normal prior distribution 
    #' for the structural matrix \eqn{B}. This precision matrix is equation invariant.
    B_V_inv    = matrix(),
    
    #' @field B_nu a positive integer greater of equal than \code{N}, a shape parameter of 
    #' the generalised-normal prior distribution for the structural matrix \eqn{B}.
    B_nu       = NA,
    
    #' @field hyper_nu_B a positive scalar, the shape parameter of the inverted-gamma 2 prior
    #' for the overall shrinkage parameter for matrix \eqn{B}.
    hyper_nu_B = NA,
    
    #' @field hyper_a_B a positive scalar, the shape parameter of the gamma prior
    #' for the second-level hierarchy for the overall shrinkage parameter for matrix \eqn{B}.
    hyper_a_B  = NA,
    
    #' @field hyper_s_BB a positive scalar, the scale parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{B}.
    hyper_s_BB  = NA,
    
    #' @field hyper_nu_BB a positive scalar, the shape parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{B}.
    hyper_nu_BB  = NA,
    
    #' @field hyper_nu_A a positive scalar, the shape parameter of the inverted-gamma 2 prior 
    #' for the overall shrinkage parameter for matrix \eqn{A}.
    hyper_nu_A  = NA,
    
    #' @field hyper_a_A a positive scalar, the shape parameter of the gamma prior
    #' for the second-level hierarchy for the overall shrinkage parameter for matrix \eqn{A}.
    hyper_a_A  = NA,
    
    #' @field hyper_s_AA a positive scalar, the scale parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{A}.
    hyper_s_AA  = NA,
    
    #' @field hyper_nu_AA a positive scalar, the shape parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{A}.
    hyper_nu_AA  = NA,
    
    #' @field sv_a_ a positive scalar, the shape parameter of the gamma prior in 
    #' the hierarchical prior for \eqn{\sigma^2_{\omega}}. 
    sv_a_      = numeric(),
    
    #' @field sv_s_ a positive scalar, the scale parameter of the gamma prior in 
    #' the hierarchical prior for \eqn{\sigma^2_{\omega}}.
    sv_s_      = numeric(),
    
    #' @field PR_TR an \code{MxM} matrix, the matrix of hyper-parameters of the 
    #' row-specific Dirichlet prior distribution for transition probabilities 
    #' matrix \eqn{P} of the Markov process \eqn{s_t}. 
    PR_TR      = matrix(),
    
    #' @field df_a a positive scalar, the rate parameter of the shifted exponential 
    #' prior distribution for the degrees of freedom parameters of the Student-t 
    #' distribution for structural shocks.
    df_a      = numeric(),
    
    #' @description
    #' Create a new prior specification PriorBSVARTVP
    #' @param train_data a \code{T_train x N} matrix of training sample for the 
    #' \eqn{B} matrix prior distribution.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param M a positive integer - the number of Markov regimes.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in 
    #' the model.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} 
    #' sets the prior mean for the autoregressive parameters of the \code{N}th 
    #' equation to the white noise process, otherwise to random walk.
    #' @return A new prior specification PriorBSVARTVP
    #' @examples 
    #' # a prior for 3-variable example with one lag
    #' prior = specify_prior_bsvarTVP$new(us_fiscal_lsuw[1:10,], N = 3, M = 2, p = 1)
    #' prior$A # show autoregressive prior mean
    #' 
    initialize = function(train_data = NULL, N, M, p, d = 0, stationary = rep(FALSE, N)){
      
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument M must be a positive integer number." = M > 0 & M %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))
      
      A1                = diag(as.numeric(!stationary))
      
      if (is.null(train_data)) {
        T_train         = 0
        S_inv           = diag(N)
      } else {
        stopifnot("Argument train_data has to be a matrix." = is.matrix(train_data) & is.numeric(train_data))
        stopifnot("Argument train_data has to contain at least 2 columns and 3 rows." = (ncol(train_data) >= 2 & nrow(train_data) >= 3))
        stopifnot("Argument train_data cannot include missing values." = sum(is.na(train_data)) == 0 )
        T_train         = nrow(train_data)
        resid           = train_data[2:T_train,] - train_data[1:(T_train - 1),] %*% t(A1)
        S_inv           = crossprod(resid)
      }
      
      K                 = N * p + 1 + d
      self$A            = cbind(A1, matrix(0, N, K - N))
      self$A_V_inv      = diag(c(kronecker((1:p)^2, rep(1, N) ), rep(1, d + 1)))
      self$B_V_inv      = S_inv
      self$B_nu         = T_train + N + 1
      self$hyper_nu_B   = 10
      self$hyper_a_B    = 10
      self$hyper_s_BB   = 100
      self$hyper_nu_BB  = 1
      self$hyper_nu_A   = 10
      self$hyper_a_A    = 10
      self$hyper_s_AA   = 10
      self$hyper_nu_AA  = 10
      self$sv_a_        = 1
      self$sv_s_        = 0.1
      self$PR_TR        = matrix(1, M, M) + 11 * diag(M)
      self$df_a         = 0.107
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBSVAR as a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 3-variable example
    #' prior = specify_prior_bsvarTVP$new(us_fiscal_lsuw[1:10,], N = 3, M = 2, p = 1)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        A        = self$A,
        A_V_inv  = self$A_V_inv,
        B_V_inv  = self$B_V_inv,
        B_nu     = self$B_nu,
        hyper_nu_B  = self$hyper_nu_B,
        hyper_a_B   = self$hyper_a_B,
        hyper_s_BB  = self$hyper_s_BB,
        hyper_nu_BB = self$hyper_nu_BB,
        hyper_nu_A  = self$hyper_nu_A,
        hyper_a_A   = self$hyper_a_A,
        hyper_s_AA  = self$hyper_s_AA,
        hyper_nu_AA = self$hyper_nu_AA,
        sv_a_       = self$sv_a_,
        sv_s_       = self$sv_s_,
        PR_TR       = self$PR_TR,
        df_a        = self$df_a
      )
    } # END get_prior
    
  ) # END public
) # END specify_prior_bsvarTVP





#' R6 Class Representing StartingValuesBSVARTVPMS
#'
#' @description
#' The class StartingValuesBSVARTVPMS presents starting values for the homoskedastic bsvar model.
#' 
#' @examples 
#' # starting values for a model for a 3-variable system
#' A = matrix(TRUE, 3, 4)
#' B = matrix(TRUE, 3, 3)
#' sv = specify_starting_values_bsvarTVPms$new(A = A, B = B, N = 3, M = 2, T = 120, p = 1)
#' 
#' @export
specify_starting_values_bsvarTVPms = R6::R6Class(
  "StartingValuesBSVARTVPMS",
  
  public = list(
    
    #' @field A an \code{NxK} matrix of starting values for the parameter \eqn{A}. 
    A             = matrix(),
    
    #' @field B an \code{NxNxM} array of starting values for the regime-dependent 
    #' parameter \eqn{B}. 
    B             = array(),
    
    #' @field hyper a list of starting values for the shrinkage hyper-parameters of the 
    #' hierarchical prior distribution. 
    hyper         = list(),
    
    #' @field h an \code{NxT} matrix with the starting values of the log-volatility processes.
    h             = matrix(),
    
    #' @field omega an \code{NxM} matrix with values of regime-dependent SV process 
    #' conditional standard deviations.
    omega         = matrix(),
    
    #' @field sigma2v an \code{NxM} matrix with values of regime-dependent SV process 
    #' conditional variances.
    sigma2v       = matrix(),
    
    #' @field rho an \code{Nx1} matrix with values of SV autoregressive parameters.
    rho           = matrix(),
    
    #' @field S an \code{NxT} integer matrix with the auxiliary mixture component indicators.
    S             = matrix(),
    
    #' @field sigma2_omega an \code{Nx1} matrix with variances of the zero-mean 
    #' normal prior for \eqn{\omega_n}.
    sigma2_omega  = numeric(),
    
    #' @field s_ an \code{Nx1} matrix of positive scalars with the scales of the 
    #' gamma prior of the hierarchical prior for \eqn{\sigma^2_{\omega}}.
    s_            = matrix(),
    
    #' @field S4_indicator a \code{NxM} matrix of initial identification allocations.
    S4_indicator  = matrix(),
    
    #' @field sigma a \code{NxT} matrix of starting values for conditional 
    #' standard deviations.
    sigma         = matrix(),
    
    #' @field PR_TR a \code{MxM} matrix of starting values for transition probabilities.
    PR_TR         = matrix(),
    
    #' @field pi_0 an \code{M}-vector of starting values for initial regime probabilities.
    pi_0          = numeric(),
    
    #' @field xi a \code{MxT} matrix of starting values for regime allocations.
    xi            = matrix(),
    
    #' @field lambda a \code{NxT} matrix of starting values for latent variables.
    lambda        = matrix(),
    
    #' @field df an \code{NxM} matrix of positive numbers with starting values 
    #' for the equation- and regime-specific degrees of freedom parameters of 
    #' the Student-t conditional distribution of structural shocks.
    df            = matrix(),
    
    #' @description
    #' Create new starting values StartingValuesBSVARTVPMS.
    #' @param A a logical \code{NxK} matrix containing value \code{TRUE} for the elements of 
    #' the autoregressive matrix \eqn{A} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the staructural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param M a positive integer - the number of Markov regimes.
    #' @param T a positive integer - the number of time periods in the data.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return Starting values StartingValuesBSVARTVPMS
    #' @examples 
    #' # starting values for a model with 4 lags for a 3-variable system
    #' A = matrix(TRUE, 3, 4)
    #' B = matrix(TRUE, 3, 3)
    #' sv = specify_starting_values_bsvarTVPms$new(A = A, B = B, N = 3, M = 2, T = 120, p = 1)
    #' 
    initialize = function(A, B, N, M, T, p, d = 0){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument M must be a positive integer number." = M > 0 & M %% 1 == 0)
      stopifnot("Argument T must be a positive integer number." = T > 0 & T %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)

      K                   = N * p + 1 + d
      self$B              = array(0, c(N, N, M))
      for (m in 1:M) {
        diag(self$B[,,m])[diag(B)] = runif(sum(diag(B)))
      }
      self$A              = matrix(0, N, K)
      diag(self$A)[diag(A[,1:N])] = runif(sum(diag(A[,1:N])))
      self$hyper          = matrix(10, 2 * N + 1, 2)
      
      self$S4_indicator   = matrix(1, N, M) 
      self$sigma          = matrix(1, N, T)
      
      self$lambda         = matrix(1, N, T)
      self$df             = matrix(3, N, M)
      
      self$h              = matrix(rnorm(N * T, sd = .01), N, T)
      self$rho            = rep(.5, N)
      self$omega          = rep(.1, N)
      self$sigma2v        = rep(.1^2, N)
      self$S              = matrix(1, N, T)
      self$sigma2_omega   = rep(1, N)
      self$s_             = rep(0.05, N)
      
      self$PR_TR          = diag(M)
      self$xi             = diag(M)[,sample(1:M, T, replace = TRUE)]
      self$pi_0           = rep(1/M, M)
    }, # END initialize
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBSVARTVPMS as a \code{list}.
    #' 
    #' @examples 
    #' # starting values for a model with 1 lag for a 3-variable system
    #' A = matrix(TRUE, 3, 4)
    #' B = matrix(TRUE, 3, 3)
    #' sv = specify_starting_values_bsvarTVPms$new(A = A, B = B, N = 3, M = 2, T = 120, p = 1)
    #' sv$get_starting_values()   # show starting values as list
    #' 
    get_starting_values   = function(){
      list(
        B                 = self$B,
        A                 = self$A,
        hyper             = self$hyper,
        lambda            = self$lambda,
        sigma             = self$sigma,
        df                = self$df,
        h                 = self$h,
        rho               = self$rho,
        omega             = self$omega,
        sigma2v           = self$sigma2v,
        S                 = self$S,
        sigma2_omega      = self$sigma2_omega,
        s_                = self$s_,
        S4_indicator      = self$S4_indicator,
        PR_TR             = self$PR_TR,
        xi                = self$xi,
        pi_0              = self$pi_0
      )
    }, # END get_starting_values
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBSVARTVPMS as a \code{list}.
    #' @param last_draw a list containing the last draw of elements \code{B} - an \code{NxN} matrix, 
    #' \code{A} - an \code{NxK} matrix, and \code{hyper} - a vector of 5 positive real numbers.
    #' @return An object of class StartingValuesBSVARTVPMS including the last draw of the current MCMC 
    #' as the starting value to be passed to the continuation of the MCMC estimation using \code{estimate()}.
    #' 
    #' @examples 
    #' # starting values for a model with 1 lag for a 3-variable system
    #' A = matrix(TRUE, 3, 4)
    #' B = matrix(TRUE, 3, 3)
    #' sv = specify_starting_values_bsvarTVPms$new(A = A, B = B, N = 3, M = 2, T = 120, p = 1)
    #' 
    #' # Modify the starting values by:
    #' sv_list = sv$get_starting_values()   # getting them as list
    #' sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    #' sv$set_starting_values(sv_list)      # providing to the class object
    #' 
    set_starting_values   = function(last_draw) {
        self$B            = last_draw$B
        self$A            = last_draw$A
        self$hyper        = last_draw$hyper
        self$lambda       = last_draw$lambda
        self$df           = last_draw$df
        self$sigma        = last_draw$sigma
        self$S4_indicator = last_draw$S4_indicator
        self$h            = last_draw$h
        self$rho          = last_draw$rho
        self$omega        = last_draw$omega
        self$sigma2v      = last_draw$sigma2v
        self$S            = last_draw$S
        self$sigma2_omega = last_draw$sigma2_omega
        self$s_           = last_draw$s_
        self$PR_TR        = last_draw$PR_TR
        self$xi           = last_draw$xi
        self$pi_0         = last_draw$pi_0
    } # END set_starting_values
  ) # END public
) # END specify_starting_values_bsvarTVPms





#' R6 Class Representing StartingValuesBSVARTVPMSA
#'
#' @description
#' The class StartingValuesBSVARTVPMSA presents starting values for the homoskedastic bsvar model.
#' 
#' @examples 
#' # starting values for a homoskedastic bsvar for a 3-variable system
#' A = matrix(TRUE, 3, 4)
#' B = matrix(TRUE, 3, 3)
#' sv = specify_starting_values_bsvarTVPmsa$new(A = A, B = B, N = 3, M = 2, T = 120, p = 1)
#' 
#' @export
specify_starting_values_bsvarTVPmsa = R6::R6Class(
  "StartingValuesBSVARTVPMSA",
  
  inherit = specify_starting_values_bsvarTVPms,
  
  public = list(
    
    #' @field A an \code{NxKxM} array of starting values for the regime-specific 
    #' parameter \eqn{A}. 
    A             = array(),
    
    #' @description
    #' Create new starting values StartingValuesBSVARTVPMSA.
    #' @param A a logical \code{NxK} matrix containing value \code{TRUE} for the elements of 
    #' the autoregressive matrix \eqn{A} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the staructural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param M a positive integer - the number of Markov regimes.
    #' @param T a positive integer - the number of time periods in the data.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return Starting values StartingValuesBSVARTVPMSA
    #' @examples 
    #' # starting values for a homoskedastic bsvar with 4 lags for a 3-variable system
    #' A = matrix(TRUE, 3, 4)
    #' B = matrix(TRUE, 3, 3)
    #' sv = specify_starting_values_bsvarTVPmsa$new(A = A, B = B, N = 3, M = 2, T = 120, p = 1)
    #' 
    initialize = function(A, B, N, M, T, p, d = 0){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument M must be a positive integer number." = M > 0 & M %% 1 == 0)
      stopifnot("Argument T must be a positive integer number." = T > 0 & T %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      
      K                   = N * p + 1 + d
      
      super$initialize(A, B, N, M, T, p, d)
      
      self$A              = array(0, c(N, K, M))
      for (m in 1:M) {
        diag(self$A[,,m])[diag(A[,1:N])] = runif(sum(diag(A[,1:N])))
      }
    }#, # END initialize
  ) # END public
) # END specify_starting_values_bsvarTVPmsa











#' R6 Class Representing IdentificationBSVARTVIs
#'
#' @description
#' The class IdentificationBSVARTVIs presents the identifying restrictions for the bsvar models.
#' 
#' @examples 
#' specify_identification_bsvarsTVI$new(N = 3, K = 4) # recursive specification for a 3-variable system
#' 
#' B = matrix(c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE), 3, 3); B
#' specify_identification_bsvarsTVI$new(B = B, N = 3, K = 4) # an alternative identification pattern
#' 
#' @export
specify_identification_bsvarsTVI = R6::R6Class(
  "IdentificationBSVARTVIs",
  
  public = list(
    
    #' @field VB a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{B}. 
    VB    = list(),
    
    #' @field VA a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{A}. 
    VA    = list(),
    
    #' @description
    #' Create new identifying restrictions IdentificationBSVARTVIs
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the structural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param A a logical \code{NxK} matrix containing value \code{TRUE} for the elements of 
    #' the autoregressive matrix \eqn{A} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param K a positive integer - the number of parameters in a row of autoregressive matrix.
    #' @return Identifying restrictions IdentificationBSVARTVIs
    initialize = function(B, A, N, K) {
      if (missing(B)) {
        B     = matrix(FALSE, N, N)
        B[lower.tri(B, diag = TRUE)] = TRUE
      }
      if (missing(A)) {
        A     = matrix(TRUE, N, K)
      }
      
      stopifnot("Argument B must be an NxN matrix with logical values." = is.logical(B) & is.matrix(B) & prod(dim(B) == N))
      stopifnot("Argument A must be an NxK matrix with logical values." = is.logical(A) & is.matrix(A) & prod(dim(A) == c(N, K)))
      
      self$VB          <- vector("list", N)
      self$VA          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
        self$VA[[n]]   <- matrix(diag(K)[A[n,],], ncol = K)
      }
      self$VB[[n + 1]] <- as.matrix(rep(1, N))
    }, # END initialize
    
    #' @description
    #' Returns the elements of the identification pattern IdentificationBSVARTVIs as a \code{list}.
    #' 
    #' @examples 
    #' B    = matrix(c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE), 3, 3); B
    #' spec = specify_identification_bsvarsTVI$new(B = B, N = 3, K = 4)
    #' spec$get_identification()
    #' 
    get_identification = function() {
      list(
        VB = self$VB,
        VA = self$VA
      )
    }, # END get_identification
    
    #' @description
    #' Set new starting values StartingValuesBSVAR.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the structural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param A a logical \code{NxK} matrix containing value \code{TRUE} for the elements of 
    #' the autoregressive matrix \eqn{A} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param K a positive integer - the number of parameters in a row of autoregressive matrix.
    #' 
    #' @examples 
    #' spec = specify_identification_bsvarsTVI$new(N = 3, K = 4) # specify a model with the default option
    #' B    = matrix(c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE), 3, 3); B
    #' spec$set_identification(B = B, N = 3, K = 4)  # modify an existing specification
    #' spec$get_identification()              # check the outcome
    set_identification = function(B, A, N, K) {
      if (missing(B)) {
        B     = matrix(FALSE, N, N)
        B[lower.tri(B, diag = TRUE)] = TRUE
      }
      if (missing(A)) {
        A     = matrix(TRUE, N, K)
      }
      
      stopifnot("Argument B must be an NxN matrix with logical values." = is.logical(B) & is.matrix(B) & prod(dim(B) == N))
      stopifnot("Argument A must be an NxK matrix with logical values." = is.logical(A) & is.matrix(A) & prod(dim(A) == c(N, K)))
      
      self$VB          <- vector("list", N)
      self$VA          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
        self$VA[[n]]   <- matrix(diag(K)[A[n,],], ncol = K)
      }
    }, # END set_identification
    
    
    #' @description
    #' Adds a new restriction to the restrictions identifying the shock as specified 
    #' in \code{shock} for the model to select from.
    #' @param restriction an \code{N}-logical vector with values \code{TRUE} for 
    #' the parameters to be estimated and \code{FALSE} for the parameters to be 
    #' restricted to zero.
    #' @param shock a positive integer specifying to which shock add the restriction.
    #' @examples 
    #' spec = specify_bsvarTVP$new(us_fiscal_lsuw)
    #' spec$add_restriction(c(TRUE, FALSE, TRUE), shock = 2)
    add_restriction = function(restriction, shock) {
      
      N     = length(self$VA)
      V     = length(self$VB)
      
      stopifnot(
        "The argument restriction must be a logical vector of length N." = 
        is.logical(restriction) & length(restriction) == N & all(is.na(restriction) == FALSE) 
      )
      stopifnot(
        "Argument shock has to be a positive integer in {1,..,N}." = 
        (shock %% 1) == 0 & shock > 0 & shock <= N
      )
      
      C  = tail(
        cumsum(self$VB[[V]][1:shock]),
        1
      )
      for (i in V:C) {
        self$VB[[i + 1]] = self$VB[[i]]
      }
      self$VB[[C]] = matrix( diag(N)[restriction,], ncol = N)
      self$VB[[V + 1]][shock] = self$VB[[V + 1]][shock] + 1
    } # END get_normal
    
  ) # END public
) # END specify_identification_bsvarsTVI





#' R6 Class representing the specification of the homoskedastic BSVAR model
#'
#' @description
#' The class BSVAR presents complete specification for the homoskedastic bsvar model.
#' 
#' @examples 
#' spec = specify_bsvarTVP$new(
#'    data = us_fiscal_lsuw,
#'    p = 4
#' )
#' 
#' @export
specify_bsvarTVP = R6::R6Class(
  "BSVARTVP",
  
  private = list(
    normal = TRUE,    # if TRUE - normal shocks, if FALSE - Student-t shocks
    msa    = FALSE    # if TRUE - MSSA, if FALSE - MSS
  ), # END private
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field M a positive integer specifying the number of Markov regimes in the model. 
    M                      = numeric(),
    
    #' @field identification an object IdentificationBSVARTVIs with the identifying restrictions. 
    identification         = list(),
    
    #' @field prior an object PriorBSVARTVP with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object DataMatricesBSVAR with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object StartingValuesBSVARTVPMS or 
    #' StartingValuesBSVARTVPMSA with the starting values.
    starting_values        = list(),
    
    #' @description
    #' Create a new specification of the homoskedastic bsvar model BSVAR.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param M a positive integer specifying the number of Markov regimes.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the structural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param A a logical \code{NxK} matrix containing value \code{TRUE} for the elements of 
    #' the autoregressive matrix \eqn{A} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param train_data a positive integer specifying the number of initial observations
    #' to be used as training sample to be used to train the prior distribution for \eqn{B}.
    #' @param distribution a character string specifying the conditional distribution 
    #' of structural shocks. Value \code{"norm"} sets it to the normal distribution, 
    #' while value \code{"t"} sets the Student-t distribution.
    #' @param ms4A a logical value - if \code{TRUE} the Markov switching is implemented
    #' in both matrices \eqn{A} and \eqn{B}, otherwise only in matrix \eqn{B}.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables. 
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets 
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process, 
    #' otherwise to random walk.
    #' @return A new complete specification for the homoskedastic bsvar model BSVAR.
    initialize = function(
    data,
    p = 1L,
    M = 2L,
    B,
    A,
    train_data = 0L,
    distribution = c("norm","t"),
    ms4A = FALSE,
    exogenous = NULL,
    stationary = rep(FALSE, ncol(data))
    ) {
      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      self$p     = p
      
      stopifnot("Argument M has to be a positive integer." = ((M %% 1) == 0 & M > 0))
      self$M     = M
      
      stopifnot("Argument train_data has to be a positive integer." = ((train_data %% 1) == 0 & train_data >= 0))
      stopifnot("Argument train_data has to much less than the data size." = train_data < nrow(data) / 2)
      stopifnot("Argument ms4A has to be a logical value." = is.logical(ms4A) & length(ms4A) == 1)
      
      distribution  = match.arg(distribution)
      
      TT            = nrow(data)
      T             = TT - self$p
      N             = ncol(data)
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous)
        exogenous   = exogenous[-(1:train_data),]
      }
      K             = N * p + 1 + d
      
      if (missing(B)) {
        message("The identification is set to the default option of lower-triangular structural matrix.")
        B     = matrix(FALSE, N, N)
        B[lower.tri(B, diag = TRUE)] = TRUE
      }
      stopifnot("Incorrectly specified argument B." = (is.matrix(B) & is.logical(B)) | (length(B) == 1 & is.na(B)))
      if (missing(A)) {
        A     = matrix(TRUE, N, K)
      }
      stopifnot("Incorrectly specified argument A." = (is.matrix(A) & is.logical(A)))
      
      if (distribution == "t") {
        private$normal = FALSE
      }
      
      if (train_data == 0) {
        data_train    = NULL
        data_estimate = data
      } else {
        data_train    = data[1:train_data,]
        data_estimate = data[-(1:train_data),]
      }
      
      self$data_matrices   = specify_data_matrices$new(data_estimate, p, exogenous)
      self$identification  = specify_identification_bsvarsTVI$new(B, A, N, K)
      self$prior           = specify_prior_bsvarTVP$new(data_train, N, M, p, d, stationary)
      if ( ms4A ) {
        self$starting_values = specify_starting_values_bsvarTVPmsa$new(A, B, N, M, T, p, d)
      } else {
        self$starting_values = specify_starting_values_bsvarTVPms$new(A, B, N, M, T, p, d)
      }
    }, # END initialize
    
    #' @description
    #' Adds a new restriction to the restrictions identifying the shock as specified 
    #' in \code{shock} for the model to select from.
    #' @param restriction an \code{N}-logical vector with values \code{TRUE} for 
    #' the parameters to be estimated and \code{FALSE} for the parameters to be 
    #' restricted to zero.
    #' @param shock a positive integer specifying to which shock add the restriction.
    #' @examples 
    #' spec = specify_bsvarTVP$new(us_fiscal_lsuw)
    #' spec$add_restriction(c(TRUE, FALSE, TRUE), shock = 2)
    add_restriction = function(restriction, shock) {
      self$identification$add_restriction(restriction, shock)
    }, # END get_normal
    
    #' @description
    #' Returns the logical value of whether the conditional shock distribution is normal.
    #' 
    #' @examples 
    #' spec = specify_bsvarTVP$new(us_fiscal_lsuw)
    #' spec$get_normal()
    get_normal = function() {
      private$normal
    }, # END get_normal
    
    #' @description
    #' Returns the data matrices as the DataMatricesBSVAR object.
    #' 
    #' @examples 
    #' spec = specify_bsvarTVP$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_data_matrices()
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the identifying restrictions as the IdentificationBSVARTVIs object.
    #' 
    #' @examples 
    #' spec = specify_bsvarTVP$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_identification()
    get_identification = function() {
      self$identification$clone()
    }, # END get_identification
    
    #' @description
    #' Returns the prior specification as the PriorBSVARTVP object.
    #' 
    #' @examples 
    #' spec = specify_bsvarTVP$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_prior()
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the StartingValuesBSVAR object.
    #' 
    #' @examples 
    #' spec = specify_bsvarTVP$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_starting_values()
    get_starting_values = function() {
      self$starting_values$clone()
    } # END get_starting_values
  ) # END public
) # END specify_bsvarTVP


















