
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "sample_ABhyper.h"
#include "sample_sv_ms.h"

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
    const int                     hyper_select = 1 // {1 - horseshoe, 2 - boost, 3 - no hierarchy};
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
  List  aux_hyper   = as<List>(starting_values["hyper"]);  // (2*N+1)x2 (gamma_0, gamma_+, s_0, s_+, s_)
  ivec  aux_SL      = as<ivec>(starting_values["S4_indicator"]) - 1;      // S4 indicator vector
  mat   aux_sigma(N, T, fill::ones);
  
  const int   S     = floor(SS / thin);
  
  cube  posterior_B(N, N, S);
  cube  posterior_A(N, K, S);
  List  posterior_hyper;
  imat  posterior_SL(N, S);
  
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
