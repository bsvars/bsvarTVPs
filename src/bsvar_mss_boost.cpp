
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "sample_ABhyper.h"
#include "sample_sv_ms.h"
#include "sample_mst.h"

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
    const int                     thin = 100,  // introduce thinning
    const int                     hyper_select = 1,
    const bool                    studentt = false
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
  
  cube  aux_B       = as<cube>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  List  aux_hyper   = as<List>(starting_values["hyper"]);  // (2*N+1)x2 (gamma_0, gamma_+, s_0, s_+, s_)
  
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  
  const int M       = aux_PR_TR.n_cols;
  vec      Tm       = sum(aux_xi, 1);
  
  // parameters for adaptive sampling of degrees of freedom
  NumericVector aag = {0.44, 0.6};
  const vec     adptive_alpha_gamma = as<vec>(aag);
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for df evaluated at df = 30
  double  adaptive_scale_init = pow(R::psigamma(15, 1) - 29 * pow(28, -2), -1) / (T / M);
  mat     adaptive_scale(N, M, fill::value(adaptive_scale_init));
  
  mat   aux_lambda(N, T, fill::ones);
  mat   aux_sigma(N, T, fill::ones);
  mat   aux_df = as<mat>(starting_values["df"]);
  
  const int   S     = floor(SS / thin);
  
  field<cube> posterior_B(S);
  cube  posterior_A(N, K, S);
  cube  posterior_hyper(2 * N + 1, 2, S);
  
  cube  posterior_PR_TR(M, M, S);
  mat   posterior_pi_0(M,S);
  cube  posterior_xi(M, T, S);
  
  cube  posterior_lambda(N, T, S);
  cube  posterior_df(N, M, S);
  
  vec   acceptance_count(4);
  List  PR_TR_tmp;
  
  field<mat> precisionB;
  field<mat> precisionA;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_lambda and aux_df
    mat E             = (Y - aux_A * X);
    if ( studentt ) {
      
      mat U           = E;
      for (int t=0; t<T; t++) {
        int   m       = aux_xi.col(t).index_max();
        U.col(t)      = aux_B.slice(m) * (Y.col(t) - aux_A * X.col(t));
      }
      
      aux_lambda      = sample_lambda_ms(aux_df, aux_xi, U);
      aux_sigma       = pow(aux_lambda, 0.5);
      
      List aux_df_tmp = sample_df_ms (aux_df, aux_lambda, aux_xi, U, prior, ss, adaptive_scale, adptive_alpha_gamma);
      aux_df          = as<mat>(aux_df_tmp["aux_df"]);
      adaptive_scale  = as<mat>(aux_df_tmp["adaptive_scale"]);
    } // END studentt
    
    // sample aux_xi
    cube Z(N, T, M);
    for (int m=0; m<M; m++) {
      Z.slice(m)    = pow(aux_lambda, -0.5) % (aux_B.slice(m) * (Y - aux_A * X) );
    }
    aux_xi            = sample_Markov_process(Z, aux_xi, aux_PR_TR, aux_pi_0, true);
    
    // sample aux_PR_TR and aux_pi_0
    PR_TR_tmp         = sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior);
    aux_PR_TR         = as<mat>(PR_TR_tmp["aux_PR_TR"]);
    aux_pi_0          = as<vec>(PR_TR_tmp["aux_pi_0"]);
    
    // sample aux_hyper
    if ( hyper_select == 1 ) {
      
      aux_hyper           = sample_hyperparameter_mss_horseshoe(aux_hyper, aux_B, aux_A, VB, prior);
      precisionB          = hyper2precisionB_mss_horseshoe(aux_hyper);
      precisionA          = hyper2precisionA_horseshoe(aux_hyper);
      
    } else if ( hyper_select == 2 ) {
      
      aux_hyper         = sample_hyperparameters_mss_boost( aux_hyper, aux_B, aux_A, VB, prior, true);
      precisionB        = hyper2precisionB_mss_boost(aux_hyper, prior);
      precisionA        = hyper2precisionA_boost(aux_hyper, prior);
      
    } else if ( hyper_select == 3 ) {
      
      aux_hyper         = sample_hyperparameters_mss_boost( aux_hyper, aux_B, aux_A, VB, prior, false);
      precisionB        = hyper2precisionB_mss_boost(aux_hyper, prior);
      precisionA        = hyper2precisionA_boost(aux_hyper, prior);
      
    }
    
    // sample aux_B
    aux_B             = sample_B_mss(aux_B, aux_A, precisionB, aux_sigma, aux_xi, Y, X, prior, VB);
    
    // sample aux_A
    aux_A             = sample_A_heterosk1_mss(aux_A, aux_B, aux_xi, precisionA, aux_sigma, Y, X, prior);
    
    if (ss % thin == 0) {
      posterior_B(s)                = aux_B;
      posterior_A.slice(s)          = aux_A;
      // posterior_hyper(s)            = aux_hyper;
      posterior_PR_TR.slice(s)      = aux_PR_TR;
      posterior_pi_0.col(s)         = aux_pi_0;
      posterior_xi.slice(s)         = aux_xi;
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
      _["PR_TR"]    = aux_PR_TR,
      _["pi_0"]     = aux_pi_0,
      _["xi"]       = aux_xi,
      _["lambda"]   = aux_lambda,
      _["df"]       = aux_df
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["PR_TR"]    = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi"]       = posterior_xi,
      _["lambda"]   = posterior_lambda,
      _["df"]       = posterior_df
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_mss_boost_cpp
