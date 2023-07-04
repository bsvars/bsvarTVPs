
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "bsvars.h"
#include "sample_ABhyper.h"
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_mss_s4_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100  // introduce thinning
) {
  // // Progress bar setup
  
  Rcout << "check 1" << endl;
  
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvarTVPs: Bayesian Structural VARs with          |" << endl;
  Rcout << "  Markov-Switching, Time-Varying Identification   |" << endl;
  Rcout << "  and Stochastic Volatility                       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR model                 |" << endl;
  Rcout << "    with Markov-switching and regime-specific     |" << endl;
  Rcout << "    time-varying identification                   |" << endl;
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
  
  cube  aux_B       = as<cube>(starting_values["B"]);
  ivec  aux_SL      = as<ivec>(starting_values["S4_indicator"]) - 1;      // S4 indicator vector
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);  // (2*N+1)x2
  
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  
  const int M       = aux_PR_TR.n_cols;
  vec      Tm       = sum(aux_xi, 1);
  
  mat   aux_sigma(N, T, fill::ones);
  
  const int   S     = floor(SS / thin);
  
  field<cube> posterior_B(S);
  cube  posterior_A(N, K, S);
  cube  posterior_hyper(2 * N + 1, 2, S);
  
  cube  posterior_PR_TR(M, M, S);
  mat   posterior_pi_0(M,S);
  cube  posterior_xi(M, T, S);
  icube posterior_SL(N, M, S);
  
  vec   acceptance_count(4);
  mat   aux_xi_tmp        = aux_xi;
  mat   aux_hyper_tmp     = aux_hyper;
  mat   aux_A_tmp         = aux_A;
  List  BSL;
  
  int   s = 0;
  
  Rcout << "check 2" << endl;
  
  for (int ss=0; ss<SS; ss++) {
    
    Rcout << "   iteration: " << ss << endl;
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_xi
    
    Rcout << "   aux_xi" << endl;
    
    mat E = Y - aux_A * X;
    Rcout << "   aux_xi1" << endl;
    aux_xi_tmp        = aux_xi;
    Rcout << "   aux_xi2" << endl;
    try {
      aux_xi_tmp      = sample_Markov_process_mss(aux_xi, E, aux_B, aux_sigma, aux_PR_TR, aux_pi_0, true);
    } catch (...) {
      acceptance_count(0)++;
    }
    Rcout << "   aux_xi3" << endl;
    aux_xi            = aux_xi_tmp;
    
    
    
    // sample aux_PR_TR
    
    Rcout << "   aux_PR_TR" << endl;
    
    bsvars::sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior);
    
    // sample aux_hyper
    
    Rcout << "   aux_hyper" << endl;
    
    aux_hyper_tmp     = aux_hyper;
    try {
      aux_hyper_tmp   = sample_hyperparameters_mss_s4_boost( aux_hyper, aux_B, aux_A, VB, aux_SL, prior);
    } catch (...) {
      acceptance_count(1)++;
    }
    aux_hyper         = aux_hyper_tmp;
    
    // sample aux_B
    
    Rcout << "   aux_B" << endl;
    
    BSL     = List::create(
      _["aux_B"]      = aux_B,
      _["aux_SL"]     = aux_SL
    );
    try {
      BSL             = sample_B_mss_s4_boost(aux_B, aux_SL, aux_A, aux_hyper, aux_sigma, aux_xi, Y, X, prior, VB);
    } catch (...) {
      acceptance_count(2)++;
    }
    aux_B             = as<cube>(BSL["aux_B"]);
    aux_SL            = as<imat>(BSL["aux_SL"]);
    
    // sample aux_A
    
    Rcout << "   aux_A" << endl;
    
    aux_A_tmp         = aux_A;
    try {
      aux_A_tmp       = sample_A_heterosk1_mss_boost(aux_A, aux_B, aux_xi, aux_hyper, aux_sigma, Y, X, prior);
    } catch (...) {
      acceptance_count(3)++;
    }
    aux_A             = aux_A_tmp;
    
    
    Rcout << "   save" << endl;
    
    if (ss % thin == 0) {
      posterior_B(s)                = aux_B;
      posterior_A.slice(s)          = aux_A;
      posterior_hyper.slice(s)      = aux_hyper;
      posterior_PR_TR.slice(s)      = aux_PR_TR;
      posterior_pi_0.col(s)         = aux_pi_0;
      posterior_xi.slice(s)         = aux_xi;
      posterior_SL.slice(s)         = aux_SL;
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
      _["xi"]       = aux_xi,
      _["S4_indicator"] = aux_SL + 1
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["PR_TR"]    = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi"]       = posterior_xi,
      _["S4_indicator"] = posterior_SL + 1
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_mss_s4_boost_cpp
