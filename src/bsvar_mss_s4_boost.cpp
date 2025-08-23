
#include <RcppArmadillo.h>
#include "progress.hpp"

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
    const int                     thin = 100,  // introduce thinning
    const bool                    hyper_boost = true
) {
  // // Progress bar setup
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
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);  // (2*N+1)x2
  
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  
  const int M       = aux_PR_TR.n_cols;
  vec      Tm       = sum(aux_xi, 1);
  
  imat  aux_SL      = as<imat>(starting_values["S4_indicator"]) - 1;      // NxM S4 indicator matrix
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
  List  PR_TR_tmp;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_xi
    mat E = Y - aux_A * X;
    aux_xi            = sample_Markov_process_mss(aux_xi, E, aux_B, aux_sigma, aux_PR_TR, aux_pi_0);
    
    // sample aux_PR_TR and aux_pi_0
    PR_TR_tmp         = sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior);
    aux_PR_TR         = as<mat>(PR_TR_tmp["aux_PR_TR"]);
    aux_pi_0          = as<vec>(PR_TR_tmp["aux_pi_0"]);
    
    // sample aux_hyper
    aux_hyper         = sample_hyperparameters_mss_s4_boost( aux_hyper, aux_B, aux_A, VB, aux_SL, prior, true);
    
    field<mat> precisionB = hyper2precisionB_mss_boost(aux_hyper, prior);
    field<mat> precisionA = hyper2precisionA_boost(aux_hyper, prior);
    
    // sample aux_B
    BSL     = List::create(
      _["aux_B"]      = aux_B,
      _["aux_SL"]     = aux_SL
    );
    BSL               = sample_B_mss_s4(aux_B, aux_SL, aux_A, precisionB, aux_sigma, aux_xi, Y, X, prior, VB);
    aux_B             = as<cube>(BSL["aux_B"]);
    aux_SL            = as<imat>(BSL["aux_SL"]);
    
    // sample aux_A
    aux_A             = sample_A_heterosk1_mss(aux_A, aux_B, aux_xi, precisionA, aux_sigma, Y, X, prior);
    
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
