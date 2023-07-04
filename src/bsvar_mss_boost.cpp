
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "bsvars.h"
#include "sample_ABhyper.h"
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_mss_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100  // introduce thinning
) {
  // // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvarTVPs: Bayesian Structural VARs with          |" << endl;
  Rcout << "  Markov-Switching, Time-Varying Identification   |" << endl;
  Rcout << "  and Stochastic Volatility                       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR model                 |" << endl;
  Rcout << "    with Markov-switching structural matrix       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Progress of the MCMC simulation for " << SS << " draws" << endl;
  Rcout << "    Every " << thin << "th draw is saved via MCMC thinning" << endl;
  Rcout << " Press Esc to interrupt the computations" << endl;
  Rcout << "**************************************************|" << endl;
  Progress p(50, true);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  const int   K     = X.n_rows;
  
  cube  aux_B       = as<cube>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);  // (2*N+1)x2 (gamma_0, gamma_+, s_0, s_+, s_)
  
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat  aux_xi       = as<mat>(starting_values["xi"]);
  
  mat   aux_sigma(N, T, fill::ones);
  
  const int M       = aux_PR_TR.n_cols;
  vec      Tm       = sum(aux_xi, 1);
  
  const int   S     = floor(SS / thin);
  
  field<cube> posterior_B(S);
  cube  posterior_A(N, K, S);
  cube  posterior_hyper(2 * N + 1, 2, S);
  
  cube  posterior_PR_TR(M, M, S);
  mat   posterior_pi_0(M,S);
  cube  posterior_xi(M, T, S);
  
  vec   acceptance_count(4);
  mat   aux_xi_tmp        = aux_xi;
  mat   aux_hyper_tmp     = aux_hyper;
  cube  aux_B_tmp         = aux_B;
  mat   aux_A_tmp         = aux_A;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_xi
    mat E             = Y - aux_A * X;
    aux_xi_tmp        = aux_xi;
    try {
      aux_xi_tmp      = sample_Markov_process_mss(aux_xi, E, aux_B, aux_sigma, aux_PR_TR, aux_pi_0);
    } catch (...) {
      acceptance_count(0)++;
    }
    aux_xi            = aux_xi_tmp;
    
    // sample aux_PR_TR
    bsvars::sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior);
    
    // sample aux_hyper
    aux_hyper_tmp     = aux_hyper;
    try {
      aux_hyper_tmp   = sample_hyperparameters_mss_boost( aux_hyper, aux_B, aux_A, VB, prior);
    } catch(...) {
      acceptance_count(1)++;
    }
    aux_hyper         = aux_hyper_tmp;
    
    // sample aux_B
    aux_B_tmp         = aux_B;
    try {
      aux_B_tmp       = sample_B_mss_boost(aux_B, aux_A, aux_hyper, aux_sigma, aux_xi, Y, X, prior, VB);
    } catch (...) {
      acceptance_count(2)++;
    }
    aux_B             = aux_B_tmp;
    
    // sample aux_A
    aux_A_tmp         = aux_A;
    try {
      aux_A_tmp       = sample_A_heterosk1_mss_boost(aux_A, aux_B, aux_xi, aux_hyper, aux_sigma, Y, X, prior);
    } catch (...) {
      acceptance_count(3)++;
    }
    aux_A             = aux_A_tmp;
    
    if (ss % thin == 0) {
      posterior_B(s)                = aux_B;
      posterior_A.slice(s)          = aux_A;
      posterior_hyper.slice(s)      = aux_hyper;
      posterior_PR_TR.slice(s)      = aux_PR_TR;
      posterior_pi_0.col(s)         = aux_pi_0;
      posterior_xi.slice(s)         = aux_xi;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["PR_TR"]    = aux_PR_TR,
      _["pi_0"]     = aux_pi_0,
      _["xi"]       = aux_xi
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["PR_TR"]    = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi"]       = posterior_xi
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_mss_boost_cpp
