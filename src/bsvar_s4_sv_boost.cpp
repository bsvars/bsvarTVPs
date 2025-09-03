
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "sample_ABhyper.h"
#include "sample_sv_ms.h"
#include "sample_mst.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_s4_sv_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100, // introduce thinning
    const bool                    centred_sv = false,
    const int                     hyper_select = 1,
    const bool                    studentt = false
) {
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvarTVPs: Bayesian Structural VARs with          |" << endl;
  Rcout << "  Markov-Switching, Time-Varying Identification   |" << endl;
  Rcout << "  and Stochastic Volatility                       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR-SV model with         |" << endl;
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
  mat   aux_h       = as<mat>(starting_values["h"]);
  vec   aux_rho     = as<vec>(starting_values["rho"]);
  vec   aux_omega   = as<vec>(starting_values["omega"]);
  vec   aux_sigma2v = as<vec>(starting_values["sigma2v"]);
  umat  aux_S       = as<umat>(starting_values["S"]);
  vec   aux_sigma2_omega = as<vec>(starting_values["sigma2_omega"]);
  vec   aux_s_      = as<vec>(starting_values["s_"]);
  ivec  aux_SL      = as<ivec>(starting_values["S4_indicator"]) - 1;      // S4 indicator vector
  mat   aux_sigma(N, T);
  
  // parameters for adaptive sampling of degrees of freedom
  NumericVector aag = {0.44, 0.6};
  const vec     adptive_alpha_gamma = as<vec>(aag);
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for df evaluated at df = 30
  double  adaptive_scale_init = pow(R::psigamma(15, 1) - 29 * pow(28, -2), -1) / (T / M);
  mat     adaptive_scale(N, M, fill::value(adaptive_scale_init));
  
  mat   aux_lambda  = as<mat>(starting_values["lambda"]);
  mat   aux_lambda_tmp = pow(aux_lambda, 0.5);
  mat   aux_df = as<mat>(starting_values["df"]);
  mat   aux_xi(M , T, fill::ones);
  
  if ( centred_sv ) {
    for (int n=0; n<N; n++) {
      aux_sigma.row(n) = exp(0.5 * aux_h.row(n));
    }
  } else {
    for (int n=0; n<N; n++) {
      aux_sigma.row(n) = exp(0.5 * aux_omega(n) * aux_h.row(n));
    }
  }
  mat   aux_sigma_tmp = aux_sigma;
  
  const int   S     = floor(SS / thin);
  
  cube  posterior_B(N, N, S);
  cube  posterior_A(N, K, S);
  List  posterior_hyper;
  
  cube  posterior_lambda(N, T, S);
  cube  posterior_df(N, M, S);
  
  cube  posterior_h(N, T, S);
  mat   posterior_rho(N, S);
  mat   posterior_omega(N, S);
  mat   posterior_sigma2v(N, S);
  ucube posterior_S(N, T, S);
  mat   posterior_sigma2_omega(N, S);
  mat   posterior_s_(N, S);
  imat  posterior_SL(N, S);
  cube  posterior_sigma(N, T, S);
  
  vec   acceptance_count(3+N);
  List  BSL;
  List  sv_n;
  
  field<mat> precisionB;
  field<mat> precisionA;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_lambda and aux_df
    mat U             = aux_B * (Y - aux_A * X);
    if ( studentt ) {
      
      U              /= aux_sigma_tmp;
      aux_lambda      = sample_lambda_ms(aux_df, aux_xi, U);
      aux_lambda_tmp  = pow(aux_lambda, 0.5);
      aux_sigma       = aux_sigma_tmp % aux_lambda_tmp;
      
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
    BSL                 = sample_B_heterosk1_s4(aux_B, aux_SL, aux_A, precisionB, aux_sigma, Y, X, prior, VB);
    aux_B               = as<mat>(BSL["aux_B"]);
    aux_SL              = as<ivec>(BSL["aux_SL"]);
    
    // sample aux_A
    aux_A               = sample_A_heterosk1(aux_A, aux_B, precisionA, aux_sigma, Y, X, prior);
    
    // sample aux_h, aux_omega and aux_S, aux_sigma2_omega
    U = aux_B * (Y - aux_A * X) / aux_lambda_tmp;
    
    for (int n=0; n<N; n++) {
      rowvec  h_tmp     = aux_h.row(n);
      double  rho_tmp   = aux_rho(n);
      double  omega_tmp = aux_omega(n);
      double  sigma2v_tmp = pow(aux_omega(n),2);
      urowvec S_tmp     = aux_S.row(n);
      rowvec  U_tmp     = U.row(n);
      double  s2o_tmp   = aux_sigma2_omega(n);
      double  s_n       = aux_s_(n);

      if ( centred_sv ) {
        sv_n            = svar_ce1( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, U_tmp, prior, true );
      } else {
        sv_n            = svar_nc1( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, U_tmp, prior, true );
      }
      
      aux_h.row(n)      = as<rowvec>(sv_n["aux_h_n"]);
      aux_rho(n)        = as<double>(sv_n["aux_rho_n"]);
      aux_omega(n)      = as<double>(sv_n["aux_omega_n"]);
      aux_sigma2v(n)    = as<double>(sv_n["aux_sigma2v_n"]);
      aux_S.row(n)      = as<urowvec>(sv_n["aux_S_n"]);
      aux_sigma2_omega(n)         = as<double>(sv_n["aux_sigma2_omega_n"]);
      aux_s_(n)         = as<double>(sv_n["aux_s_n"]);
      
      if ( centred_sv ) {
        aux_sigma_tmp.row(n)  = exp(0.5 * aux_h.row(n));
      } else {
        aux_sigma_tmp.row(n)  = exp(0.5 * aux_omega(n) * aux_h.row(n));
      }
    }
    aux_sigma       = aux_sigma_tmp % aux_lambda_tmp;
    
    if (ss % thin == 0) {
      posterior_B.slice(s)          = aux_B;
      posterior_A.slice(s)          = aux_A;
      // posterior_hyper(s)            = aux_hyper;
      posterior_lambda.slice(s)     = aux_lambda;
      posterior_df.slice(s)         = aux_df;
      posterior_h.slice(s)          = aux_h;
      posterior_rho.col(s)          = aux_rho;
      posterior_omega.col(s)        = aux_omega;
      posterior_sigma2v.col(ss)     = aux_sigma2v;
      posterior_S.slice(s)          = aux_S;
      posterior_sigma2_omega.col(s) = aux_sigma2_omega;
      posterior_s_.col(s)           = aux_s_;
      posterior_SL.col(s)           = aux_SL;
      posterior_sigma.slice(s)      = aux_sigma;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["h"]        = aux_h,
      _["rho"]      = aux_rho,
      _["omega"]    = aux_omega,
      _["sigma2v"]  = aux_sigma2v,
      _["S"]        = aux_S,
      _["sigma2_omega"] = aux_sigma2_omega,
      _["s_"]       = aux_s_,
      _["S4_indicator"] = aux_SL + 1,
      _["sigma"]    = aux_sigma,
      _["lambda"]   = aux_lambda,
      _["df"]       = aux_df
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["h"]        = posterior_h,
      _["rho"]      = posterior_rho,
      _["omega"]    = posterior_omega,
      _["sigma2v"]  = posterior_sigma2v,
      _["S"]        = posterior_S,
      _["sigma2_omega"] = posterior_sigma2_omega,
      _["s_"]        = posterior_s_,
      _["S4_indicator"] = posterior_SL + 1,
      _["sigma"]    = posterior_sigma,
      _["lambda"]   = posterior_lambda,
      _["df"]       = posterior_df
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_s4_sv_boost_cpp
