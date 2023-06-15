
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "bsvars.h"
#include "sample_ABhyper.h"

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
    const int                     thin = 100  // introduce thinning
) {
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvarTVPs: Bayesian Structural VARs with          |" << endl;
  Rcout << "  Markov-Switching and Time-Varying Identification|" << endl;
  Rcout << "  and Stochastic Volatility                       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR model with            |" << endl;
  Rcout << "    Stochastic Search Specification Selection     |" << endl;
  Rcout << "    for the structural matrix                     |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Progress of the MCMC simulation for " << SS << " draws" << endl;
  Rcout << "    Every " << thin << "th draw is saved via MCMC thinning" << endl;
  Rcout << " Press Esc to interrupt the computations" << endl;
  Rcout << "**************************************************|" << endl;
  Progress p(50, true);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  const int   K     = X.n_rows;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);  // (2*N+1)x2 (gamma_n, s_n, s) for both B and A
  ivec  aux_SL      = as<ivec>(starting_values["S4_indicator"]) - 1;      // S4 indicator vector
  mat   aux_sigma(N, T, fill::ones);
  
  const int   S     = floor(SS / thin);
  
  cube  posterior_B(N, N, S);
  cube  posterior_A(N, K, S);
  cube  posterior_hyper(2 * N + 1, 2, S);
  imat  posterior_SL(N, S);
  
  vec   acceptance_count(3);
  mat   aux_hyper_tmp       = aux_hyper;
  mat   aux_B_tmp           = aux_B;
  ivec  aux_SL_tmp          = aux_SL;
  List  BSL_tmp;
  mat   aux_A_tmp           = aux_A;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_hyper
    aux_hyper_tmp       = aux_hyper;
    try {
      aux_hyper_tmp     = sample_hyperparameter_boost_s4( aux_hyper, aux_B, aux_A, VB, aux_SL, prior);
    } catch (...) {
      acceptance_count(0)++;
    }
    aux_hyper           = aux_hyper_tmp;
    
    // sample aux_B
    aux_B_tmp           = aux_B;
    aux_SL_tmp          = aux_SL;
    try {
      BSL_tmp           = sample_B_heterosk1_s4_boost(aux_B, aux_SL, aux_A, aux_hyper, aux_sigma, Y, X, prior, VB);
    } catch (...) {
      acceptance_count(1)++;
    }
    aux_B               = as<mat>(BSL_tmp["aux_B"]);
    aux_SL              = as<ivec>(BSL_tmp["aux_SL"]);
    
    // sample aux_A
    aux_A_tmp           = aux_A;
    try {
      aux_A_tmp         = sample_A_heterosk1_boost(aux_A, aux_B, aux_hyper, aux_sigma, Y, X, prior);
    } catch (...) {
      acceptance_count(2)++;
    }
    aux_A               = aux_A_tmp;
    
    if (ss % thin == 0) {
      posterior_B.slice(s)          = aux_B;
      posterior_A.slice(s)          = aux_A;
      posterior_hyper.slice(s)      = aux_hyper;
      posterior_SL.col(s)           = aux_SL;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["S4_indicator"] = aux_SL + 1
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["S4_indicator"] = posterior_SL + 1
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_s4_boost_cpp
