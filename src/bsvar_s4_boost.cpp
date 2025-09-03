
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "sample_ABhyper.h"
#include "sample_sv_ms.h"
#include "sample_mst.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_s4_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100,
    const int                     hyper_select = 1, // {1 - horseshoe, 2 - boost, 3 - no hierarchy}
    const bool                    studentt = false 
) {
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvarTVPs: Bayesian Structural VARs with          |" << endl;
  Rcout << "  Markov-Switching, Time-Varying Identification   |" << endl;
  Rcout << "  and Stochastic Volatility                       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR model with            |" << endl;
  Rcout << "    time-varying identification                   |" << endl;
  Rcout << "    for the structural matrix                     |" << endl;
  if ( studentt ) {
    Rcout << "    and Student-t structural shocks               |" << endl;
  }
  Rcout << "**************************************************|" << endl;
  Rcout << " Progress of the MCMC simulation for " << SS << " draws" << endl;
  Rcout << "    Every " << thin << "th draw is saved via MCMC thinning" << endl;
  Rcout << " Press Esc to interrupt the computations" << endl;
  Rcout << "**************************************************|" << endl;
  Progress p(50, true);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  const int   K     = X.n_rows;
  const int   M     = 1;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  List  aux_hyper   = as<List>(starting_values["hyper"]);  // (2*N+1)x2 (gamma_0, gamma_+, s_0, s_+, s_)
  ivec  aux_SL      = as<ivec>(starting_values["S4_indicator"]) - 1;      // S4 indicator vector
  mat   aux_sigma(N, T, fill::ones);
  
  const int   S     = floor(SS / thin);
  
  cube  posterior_B(N, N, S);
  cube  posterior_A(N, K, S);
  List  posterior_hyper;
  imat  posterior_SL(N, S);
  
  cube  posterior_lambda(N, T, S);
  cube  posterior_df(N, M, S);
  
  // parameters for adaptive sampling of degrees of freedom
  NumericVector aag = {0.44, 0.6};
  const vec     adptive_alpha_gamma = as<vec>(aag);
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for df evaluated at df = 30
  double  adaptive_scale_init = pow(R::psigamma(15, 1) - 29 * pow(28, -2), -1) / (T / M);
  mat     adaptive_scale(N, M, fill::value(adaptive_scale_init));
  
  mat   aux_lambda(N, T, fill::ones);
  mat   aux_lambda_tmp(N, T, fill::ones);
  mat   aux_df = as<mat>(starting_values["df"]);
  mat   aux_xi(M, T, fill::ones);
  
  vec   acceptance_count(3);
  List  BSL;
  
  field<mat> precisionB;
  field<mat> precisionA;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_lambda and aux_df
    mat U;
    if ( studentt ) {
      
      U           = aux_B * (Y - aux_A * X);
      aux_lambda      = sample_lambda_ms(aux_df, aux_xi, U);
      aux_lambda_tmp  = pow(aux_lambda, 0.5);
      aux_sigma       = aux_lambda_tmp;
      
      List aux_df_tmp = sample_df_ms (aux_df, aux_lambda, aux_xi, U, prior, ss, adaptive_scale, adptive_alpha_gamma);
      aux_df          = as<mat>(aux_df_tmp["aux_df"]);
      adaptive_scale  = as<mat>(aux_df_tmp["adaptive_scale"]);
    } // END studentt
    
    
    // sample aux_hyper
    if ( hyper_select == 1 ) {
      
      aux_hyper           = sample_hyperparameter_horseshoe(aux_hyper, aux_B, aux_A, VB, aux_SL, prior);
      precisionB          = hyper2precisionB_horseshoe(aux_hyper);
      precisionA          = hyper2precisionA_horseshoe(aux_hyper);
      
    } else if ( hyper_select == 2 ) {
      
      aux_hyper           = sample_hyperparameter_boost_s4( aux_hyper, aux_B, aux_A, VB, aux_SL, prior, true);
      precisionB          = hyper2precisionB_boost(aux_hyper, prior);
      precisionA          = hyper2precisionA_boost(aux_hyper, prior);
      
    } else if ( hyper_select == 3 ) {
      
      aux_hyper           = sample_hyperparameter_boost_s4( aux_hyper, aux_B, aux_A, VB, aux_SL, prior, false);
      precisionB          = hyper2precisionB_boost(aux_hyper, prior);
      precisionA          = hyper2precisionA_boost(aux_hyper, prior);
      
    }
    
    // sample aux_B
    BSL     = List::create(
      _["aux_B"]      = aux_B,
      _["aux_SL"]     = aux_SL
    );
    BSL                 = sample_B_heterosk1_s4(aux_B, aux_SL, aux_A, precisionB, aux_sigma, Y, X, prior, VB);
    aux_B               = as<mat>(BSL["aux_B"]);
    aux_SL              = as<ivec>(BSL["aux_SL"]);
    
    // sample aux_A
    aux_A               = sample_A_heterosk1(aux_A, aux_B, precisionA, aux_sigma, Y, X, prior);
    
    if (ss % thin == 0) {
      posterior_B.slice(s)          = aux_B;
      posterior_A.slice(s)          = aux_A;
      // posterior_hyper(s)            = aux_hyper;
      posterior_SL.col(s)           = aux_SL;
      posterior_lambda.slice(s)     = aux_lambda;
      posterior_df.slice(s)         = aux_df;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["S4_indicator"] = aux_SL + 1,
      _["lambda"]   = aux_lambda,
      _["df"]       = aux_df
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["S4_indicator"] = posterior_SL + 1,
      _["lambda"]   = posterior_lambda,
      _["df"]       = posterior_df
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_s4_boost_cpp
