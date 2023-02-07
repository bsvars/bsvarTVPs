
#include <RcppArmadillo.h>
#include "bsvars.h"
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_ir_ms (
    arma::field<arma::cube>&  posterior_B,        // (S)(N, N, M)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p
) {
  
  const int       N = posterior_B(0).n_rows;
  const int       M = posterior_B(0).n_slices;
  const int       S = posterior_B.n_elem;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S, M);
  
  for (int s=0; s<S; s++) {
    for (int m=0; m<M; m++) {
      aux_irfs            = bsvars::bsvars_ir1( posterior_B(s).slice(m), posterior_A.slice(s), horizon, p );
      irfs(s, m)          = aux_irfs;
    } // END m loop
  } // END s loop
  
  return irfs;
} // END bsvarTVPs_ir_ms


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_ir (
    arma::cube&               posterior_B,        // (N, N, S)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p
) {
  
  const int       N = posterior_B.n_rows;
  const int       M = 1;
  const int       S = posterior_B.n_slices;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S, M);
  
  for (int s=0; s<S; s++) {
    aux_irfs            = bsvars::bsvars_ir1( posterior_B.slice(s), posterior_A.slice(s), horizon, p );
    irfs(s, 0)          = aux_irfs;
  } // END s loop
  
  return irfs;
} // END bsvarTVPs_ir


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarTVPs_filter_forecast_smooth (
    Rcpp::List&       posterior,
    const arma::mat&  Y,
    const arma::mat&  X,
    const bool        forecasted,
    const bool        smoothed
) {
  
  field<cube> posterior_B = as<field<cube>>(posterior["B"]);
  cube  posterior_A       = as<cube>(posterior["A"]);
  cube  posterior_sigma   = as<cube>(posterior["sigma"]);
  cube  posterior_PR_TR   = as<cube>(posterior["PR_TR"]);
  mat   posterior_pi_0    = as<mat>(posterior["pi_0"]);
  
  const int   N           = Y.n_rows;
  const int   M           = posterior_PR_TR.n_rows;
  const int   T           = Y.n_cols;
  const int   S           = posterior_A.n_slices;
  
  cube  filtered_probabilities(M, T, S);
  cube  for_smo_probabilities(M, T, S);
  cube  shocks(N, T, M);
  
  for (int s=0; s<S; s++) {
    for (int m=0; m<M; m++) {
      shocks.slice(m)     = pow(posterior_sigma.slice(s), -1) % (posterior_B(s).slice(m) * ( Y - posterior_A.slice(s) * X ));
    }
    
    filtered_probabilities.slice(s)   = filtering(
      shocks, 
      posterior_PR_TR.slice(s), 
      posterior_pi_0.col(s)
    );
    
    if (forecasted) {
      for_smo_probabilities.slice(s)  = posterior_PR_TR.slice(s) * filtered_probabilities.slice(s);
    } else if (smoothed) {
      for_smo_probabilities.slice(s)  = smoothing(filtered_probabilities.slice(s), posterior_PR_TR.slice(s));
    }
  } // END s loop
  
  cube  out     = filtered_probabilities;
  if (forecasted || smoothed) {
    out         = for_smo_probabilities;
  }
  
  return out;
} // END bsvarTVPs_filter_forecast_smooth


