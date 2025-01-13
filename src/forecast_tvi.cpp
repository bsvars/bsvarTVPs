#include <RcppArmadillo.h>
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;

Rcpp::List forecast_mssa_sv (
    const arma::field<arma::cube>&  posterior_B,          // (S)(N,N,M)
    const arma::field<arma::cube>&  posterior_A,          // (S)(N,K,M)
    const arma::cube&               posterior_PR_TR,      // (M,M,S)
    const arma::mat&                posterior_xi_T,       // (M,S)
    const arma::mat&                posterior_h_T,        // (N,S)
    const arma::mat&                posterior_rho,        // (N,S)
    const arma::cube&               posterior_omega,      // (N,M,S)
    const arma::vec&                X_T,                   // (K)
    const int&                      horizon
) {
  const int   N = posterior_h_T.n_rows;
  const int   K = X_T.n_elem;
  const int   M = posterior_PR_TR.n_rows;
  const int   S = posterior_h_T.n_cols;
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
    
  NumericVector zeroM = wrap(seq_len(M) - 1);
  
  for (int s=0; s<S; s++) {
    
    int       ST(M);
    vec       XT      = X_T;
    vec       PR_ST   = posterior_xi_T.col(s);
    vec       HT      = posterior_h_T.col(s);  
    vec       sigma2T(N);
    mat       BT_inv(N, N);
    cube      SigmaT(N, N, horizon);
    
    for (int h=0; h<horizon; h++) {
      
      PR_ST       = trans(posterior_PR_TR.slice(s)) * PR_ST;
      ST          = csample_num1(zeroM, wrap(PR_ST));
      HT          = posterior_omega.slice(s).col(ST) % HT + randn(N);
      sigma2T     = exp(posterior_omega.slice(s).col(ST) % HT); 
      BT_inv      = inv(posterior_B(s).slice(ST));
      SigmaT.slice(h) = BT_inv * diagmat(sigma2T) * BT_inv.t();
      out_forecast_mean.slice(s).col(h) = posterior_A(s).slice(ST) * XT;
      vec draw    = mvnrnd( out_forecast_mean.slice(s).col(h), SigmaT.slice(h) );
      out_forecast.slice(s).col(h) = draw;
      
      if ( h != horizon - 1 ) {
        XT        = join_cols( draw, XT.subvec(N, K - 1) );
      } // END if h
    } // END h loop
    
    out_forecast_cov(s) = SigmaT;
    
  } // END s loop
  
  return List::create(
    _["forecast"]       = out_forecast,
    _["forecast_mean"]  = out_forecast_mean,
    _["forecast_cov"]   = out_forecast_cov
  );
} // END forecast_mssa_sv