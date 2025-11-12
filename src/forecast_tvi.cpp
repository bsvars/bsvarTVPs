#include <RcppArmadillo.h>
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_mssa_sv (
    arma::field<arma::cube>&  posterior_B,          // (S)(N,N,M)
    arma::field<arma::cube>&  posterior_A,          // (S)(N,K,M)
    arma::cube&               posterior_PR_TR,      // (M,M,S)
    arma::mat&                posterior_xi_T,       // (M,S)
    arma::mat&                posterior_h_T,        // (N,S)
    arma::mat&                posterior_rho,        // (N,S)
    arma::cube&               posterior_omega,      // (N,M,S)
    arma::cube&               posterior_df,         // (N,M,S)
    arma::vec&                X_T,                  // (K)
    arma::mat&                exogenous_forecast, // (horizon, d)
    const int&                horizon,
    const int                 sv_select = 1,        // {1 - non-centred, 2 - centred, 3 - homoskedastic};
    const bool                studentt = false      // {true - normal, false - Student-t};
) {
  const int   N = posterior_h_T.n_rows;
  const int   K = X_T.n_elem;
  const int   M = posterior_PR_TR.n_rows;
  const int   S = posterior_h_T.n_cols;
  const int   d = exogenous_forecast.n_cols;
  
  bool        do_exog = exogenous_forecast.is_finite();
  vec         x_t;
  if ( do_exog ) {
    x_t       = X_T.rows(0, K - 1 - d);
  } else {
    x_t       = X_T.rows(0, K - 1);
  } // END if do_exog
  vec         Xt(K);
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
  
  vec         draw(N);  
  
  NumericVector zeroM = wrap(seq_len(M) - 1);
  
  for (int s=0; s<S; s++) {
    
    if ( do_exog ) {
      Xt          = join_cols(x_t, trans(exogenous_forecast.row(0)));
    } else {
      Xt          = x_t;
    } // END if do_exog
    
    int       ST(M);
    vec       PR_ST   = posterior_xi_T.col(s);
    vec       HT      = posterior_h_T.col(s);  
    mat       BT_inv_sigma(N, N);
    cube      SigmaT(N, N, horizon);
    vec       sigmaT(N, fill::ones);
    vec       sigma2_tmp(N, fill::ones);
    vec       lambda_tmp(N, fill::ones);
    
    for (int h=0; h<horizon; h++) {
      
      // predict MS regime
      PR_ST       = trans(posterior_PR_TR.slice(s)) * PR_ST;
      ST          = csample_num1(zeroM, wrap(PR_ST));
      
      // predict SV
      if (sv_select != 3) {
        HT          = posterior_omega.slice(s).col(ST) % HT + randn(N);
        sigma2_tmp  = exp(posterior_omega.slice(s).col(ST) % HT); 
      }
      
      // predict Student-t
      if (studentt) {
        vec df_s    = posterior_df.slice(s).col(ST);
        lambda_tmp  = df_s - 2;
        lambda_tmp  = lambda_tmp / chi2rnd(df_s);
      }
      
      sigmaT        = sqrt(sigma2_tmp % lambda_tmp);
      
      // predictive density
      BT_inv_sigma  = solve(posterior_B(s).slice(ST), diagmat(sigmaT));
      mat Sigma_tmp = BT_inv_sigma * BT_inv_sigma.t();
      SigmaT.slice(h) = 0.5 * (Sigma_tmp + Sigma_tmp.t());
      out_forecast_mean.slice(s).col(h) = posterior_A(s).slice(ST) * Xt;
      
      try {
        draw        = mvnrnd( out_forecast_mean.slice(s).col(h), SigmaT.slice(h) );
      } catch (std::logic_error &e) {
        continue;
      }
      out_forecast.slice(s).col(h) = draw;
      
      // create Xs
      if ( h != horizon - 1 ) {
        if ( do_exog ) {
          Xt          = join_cols( out_forecast.slice(s).col(h), Xt.subvec(N, K - 1 - d), trans(exogenous_forecast.row(h + 1)) );
        } else {
          Xt          = join_cols( out_forecast.slice(s).col(h), Xt.subvec(N, K - 1) );
        }
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


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_mss_sv (
    arma::field<arma::cube>&  posterior_B,          // (S)(N,N,M)
    arma::cube&               posterior_A,          // (N,K,S)
    arma::cube&               posterior_PR_TR,      // (M,M,S)
    arma::mat&                posterior_xi_T,       // (M,S)
    arma::mat&                posterior_h_T,        // (N,S)
    arma::mat&                posterior_rho,        // (N,S)
    arma::cube&               posterior_omega,      // (N,M,S)
    arma::cube&               posterior_df,         // (N,M,S)
    arma::vec&                X_T,                   // (K)
    arma::mat&                exogenous_forecast, // (horizon, d)
    const int&                horizon,
    const int                 sv_select = 1,        // {1 - non-centred, 2 - centred, 3 - homoskedastic};
    const bool                studentt = false      // {true - normal, false - Student-t};
) {
  const int   N = posterior_h_T.n_rows;
  const int   K = X_T.n_elem;
  const int   M = posterior_PR_TR.n_rows;
  const int   S = posterior_h_T.n_cols;
  const int   d = exogenous_forecast.n_cols;
  
  bool        do_exog = exogenous_forecast.is_finite();
  vec         x_t;
  if ( do_exog ) {
    x_t       = X_T.rows(0, K - 1 - d);
  } else {
    x_t       = X_T.rows(0, K - 1);
  } // END if do_exog
  vec         Xt(K);
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
  
  vec         draw(N);  
  
  NumericVector zeroM = wrap(seq_len(M) - 1);
  
  for (int s=0; s<S; s++) {
    
    if ( do_exog ) {
      Xt          = join_cols(x_t, trans(exogenous_forecast.row(0)));
    } else {
      Xt          = x_t;
    } // END if do_exog
    
    int       ST(M);
    vec       PR_ST   = posterior_xi_T.col(s);
    vec       HT      = posterior_h_T.col(s);  
    mat       BT_inv_sigma(N, N);
    cube      SigmaT(N, N, horizon);
    vec       sigmaT(N, fill::ones);
    vec       sigma2_tmp(N, fill::ones);
    vec       lambda_tmp(N, fill::ones);
    
    for (int h=0; h<horizon; h++) {
      
      // predict MS regime
      PR_ST       = trans(posterior_PR_TR.slice(s)) * PR_ST;
      ST          = csample_num1(zeroM, wrap(PR_ST));
      
      // predict SV
      if (sv_select != 3) {
        HT          = posterior_omega.slice(s).col(ST) % HT + randn(N);
        sigma2_tmp  = exp(posterior_omega.slice(s).col(ST) % HT); 
      }
      
      // predict Student-t
      if (studentt) {
        vec df_s    = posterior_df.slice(s).col(ST);
        lambda_tmp  = df_s - 2;
        lambda_tmp  = lambda_tmp / chi2rnd(df_s);
      }
      
      sigmaT        = sqrt(sigma2_tmp % lambda_tmp);
      
      // predictive density
      BT_inv_sigma  = solve(posterior_B(s).slice(ST), diagmat(sigmaT));
      mat Sigma_tmp = BT_inv_sigma * BT_inv_sigma.t();
      Sigma_tmp   = 0.5 * (Sigma_tmp + Sigma_tmp.t());
      SigmaT.slice(h) = Sigma_tmp;
      out_forecast_mean.slice(s).col(h) = posterior_A.slice(s) * Xt;
      
      try {
        draw        = mvnrnd( out_forecast_mean.slice(s).col(h), Sigma_tmp );
      } catch (std::logic_error &e) {
        continue;
      }
      out_forecast.slice(s).col(h) = draw;
      
      // create Xs
      if ( h != horizon - 1 ) {
        if ( do_exog ) {
          Xt          = join_cols( out_forecast.slice(s).col(h), Xt.subvec(N, K - 1 - d), trans(exogenous_forecast.row(h + 1)) );
        } else {
          Xt          = join_cols( out_forecast.slice(s).col(h), Xt.subvec(N, K - 1) );
        }
      } // END if h
      
    } // END h loop
    
    out_forecast_cov(s) = SigmaT;
  } // END s loop
  
  return List::create(
    _["forecast"]       = out_forecast,
    _["forecast_mean"]  = out_forecast_mean,
    _["forecast_cov"]   = out_forecast_cov
  );
} // END forecast_mss_sv


