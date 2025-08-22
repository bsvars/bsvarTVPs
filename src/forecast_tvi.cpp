#include <RcppArmadillo.h>
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


bool is_non_explosive_mat (
    mat& A
) {
  int N = A.n_rows;
  int p = (A.n_cols - 1) / N;
  
  mat AA = A.cols(0, N - 1);
  
  if (p > 1) {
    for (int i=1; i<p; i++) {
      AA += A.cols(i*N, (i+1)*N - 1);
    }
  }
  
  cx_vec eigval = eig_gen( AA );
  
  return max(abs(eigval)) <= 1;
} // END is_non_explosive_mat


bool is_non_explosive_cube (
    cube& A
) {
  int M = A.n_slices;
  
  bool non_explosive = true;
  for (int m=0; m<M; m++) {
    non_explosive *= is_non_explosive_mat(A.slice(m));
  }
  
  return non_explosive;
} // END is_non_explosive_cube


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
    arma::vec&                X_T,                   // (K)
    const int&                horizon, 
    const bool                non_explosive = false
) {
  const int   N = posterior_h_T.n_rows;
  const int   K = X_T.n_elem;
  const int   M = posterior_PR_TR.n_rows;
  const int   S = posterior_h_T.n_cols;
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
    
  NumericVector zeroM = wrap(seq_len(M) - 1);
  
  int ss = 0;
  
  for (int s=0; s<S; s++) {
    
    bool      ine  = is_non_explosive_cube(posterior_A(s));
    
    if (non_explosive && !ine) {
      continue;
    }
    
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
      mat Sigma_tmp = BT_inv * diagmat(sigma2T) * BT_inv.t();
      SigmaT.slice(h) = 0.5 * (Sigma_tmp + Sigma_tmp.t());
      out_forecast_mean.slice(ss).col(h) = posterior_A(s).slice(ST) * XT;
      vec draw    = mvnrnd( out_forecast_mean.slice(ss).col(h), SigmaT.slice(h) );
      out_forecast.slice(ss).col(h) = draw;
      
      if ( h != horizon - 1 ) {
        XT        = join_cols( draw, XT.subvec(N, K - 1) );
      } // END if h
    } // END h loop
    
    out_forecast_cov(ss) = SigmaT;

    ss++;
  } // END s loop
  
  Rcout << "Fraction of used MCMC draws: " << (double) 100 * ss / S << endl;
  
  cube _forecast = out_forecast.slices(0, ss - 1);
  cube _forecast_mean = out_forecast_mean.slices(0, ss - 1);
  field<cube> _forecast_cov = out_forecast_cov.rows(0, ss - 1);
  
  return List::create(
    _["forecast"]       = _forecast,
    _["forecast_mean"]  = _forecast_mean,
    _["forecast_cov"]   = _forecast_cov
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
    arma::vec&                X_T,                   // (K)
    const int&                horizon, 
    const bool                non_explosive = false
) {
  const int   N = posterior_h_T.n_rows;
  const int   K = X_T.n_elem;
  const int   M = posterior_PR_TR.n_rows;
  const int   S = posterior_h_T.n_cols;
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
  
  NumericVector zeroM = wrap(seq_len(M) - 1);
  
  int ss = 0;
  
  for (int s=0; s<S; s++) {
    
    bool      ine  = is_non_explosive_mat(posterior_A.slice(s));
    
    if (non_explosive && !ine) {
      continue;
    }
    
    int       ST(M);
    vec       XT      = X_T;
    vec       PR_ST   = posterior_xi_T.col(s);
    vec       HT      = posterior_h_T.col(s);  
    vec       sigmaT(N);
    mat       BT_inv_sigma(N, N);
    cube      SigmaT(N, N, horizon);
    
    for (int h=0; h<horizon; h++) {
      
      PR_ST       = trans(posterior_PR_TR.slice(s)) * PR_ST;
      ST          = csample_num1(zeroM, wrap(PR_ST));
      HT          = posterior_omega.slice(s).col(ST) % HT + randn(N);
      sigmaT      = exp(0.5 * posterior_omega.slice(s).col(ST) % HT); 
      BT_inv_sigma  = solve(posterior_B(s).slice(ST), diagmat(sigmaT));
      mat Sigma_tmp = BT_inv_sigma * BT_inv_sigma.t();
      Sigma_tmp   = 0.5 * (Sigma_tmp + Sigma_tmp.t());
      SigmaT.slice(h) = Sigma_tmp;
      out_forecast_mean.slice(ss).col(h) = posterior_A.slice(s) * XT;
      vec draw    = mvnrnd( out_forecast_mean.slice(ss).col(h), Sigma_tmp );
      out_forecast.slice(ss).col(h) = draw;
      
      if ( h != horizon - 1 ) {
        XT        = join_cols( draw, XT.subvec(N, K - 1) );
      } // END if h
    } // END h loop
    
    out_forecast_cov(ss) = SigmaT;
    
    ss++;
  } // END s loop
  
  Rcout << "Fraction of used MCMC draws: " << (double) 100 * ss / S << endl;
  
  cube _forecast = out_forecast.slices(0, ss - 1);
  cube _forecast_mean = out_forecast_mean.slices(0, ss - 1);
  field<cube> _forecast_cov = out_forecast_cov.rows(0, ss - 1);
  
  return List::create(
    _["forecast"]       = _forecast,
    _["forecast_mean"]  = _forecast_mean,
    _["forecast_cov"]   = _forecast_cov
  );
} // END forecast_mss_sv






// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_mss (
    arma::field<arma::cube>&  posterior_B,          // (S)(N,N,M)
    arma::cube&               posterior_A,          // (N,K,S)
    arma::cube&               posterior_PR_TR,      // (M,M,S)
    arma::mat&                posterior_xi_T,       // (M,S)
    arma::vec&                X_T,                   // (K)
    const int&                horizon, 
    const bool                non_explosive = false
) {
  const int   N = posterior_A.n_rows;
  const int   K = X_T.n_elem;
  const int   M = posterior_PR_TR.n_rows;
  const int   S = posterior_A.n_slices;
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
  
  NumericVector zeroM = wrap(seq_len(M) - 1);
  
  int ss = 0;
  
  for (int s=0; s<S; s++) {
    
    bool      ine  = is_non_explosive_mat(posterior_A.slice(s));
    
    if (non_explosive && !ine) {
      continue;
    }
    
    int       ST(M);
    vec       XT      = X_T;
    vec       PR_ST   = posterior_xi_T.col(s);
    vec       sigma2T(N, fill::ones);
    mat       BT_inv(N, N);
    cube      SigmaT(N, N, horizon);
    
    for (int h=0; h<horizon; h++) {
      
      PR_ST       = trans(posterior_PR_TR.slice(s)) * PR_ST;
      ST          = csample_num1(zeroM, wrap(PR_ST));
      BT_inv      = inv(posterior_B(s).slice(ST));
      mat Sigma_tmp = BT_inv * diagmat(sigma2T) * BT_inv.t();
      SigmaT.slice(h) = 0.5 * (Sigma_tmp + Sigma_tmp.t());
      out_forecast_mean.slice(ss).col(h) = posterior_A.slice(s) * XT;
      vec draw    = mvnrnd( out_forecast_mean.slice(ss).col(h), SigmaT.slice(h) );
      out_forecast.slice(ss).col(h) = draw;
      
      if ( h != horizon - 1 ) {
        XT        = join_cols( draw, XT.subvec(N, K - 1) );
      } // END if h
    } // END h loop
    
    out_forecast_cov(ss) = SigmaT;
    
    ss++;
  } // END s loop
  
  Rcout << "Fraction of used MCMC draws: " << (double) 100 * ss / S << endl;
  
  cube _forecast = out_forecast.slices(0, ss - 1);
  cube _forecast_mean = out_forecast_mean.slices(0, ss - 1);
  field<cube> _forecast_cov = out_forecast_cov.rows(0, ss - 1);
  
  return List::create(
    _["forecast"]       = _forecast,
    _["forecast_mean"]  = _forecast_mean,
    _["forecast_cov"]   = _forecast_cov
  );
} // END forecast_mss









// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_sv (
    arma::cube&               posterior_B,          // (N,N,S)
    arma::cube&               posterior_A,          // (N,K,S)
    arma::mat&                posterior_h_T,        // (N,S)
    arma::mat&                posterior_rho,        // (N,S)
    arma::mat&                posterior_omega,      // (N,S)
    arma::vec&                X_T,                   // (K)
    const int&                horizon, 
    const bool                non_explosive = false
) {
  const int   N = posterior_h_T.n_rows;
  const int   K = X_T.n_elem;
  const int   S = posterior_h_T.n_cols;
  
  cube        out_forecast(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  field<cube> out_forecast_cov(S);
  
  int ss = 0;
  
  for (int s=0; s<S; s++) {
    
    bool      ine  = is_non_explosive_mat(posterior_A.slice(s));
    
    if (non_explosive && !ine) {
      continue;
    }
    
    vec       XT      = X_T;
    vec       HT      = posterior_h_T.col(s);  
    vec       sigma2T(N);
    mat       BT_inv(N, N);
    cube      SigmaT(N, N, horizon);
    
    for (int h=0; h<horizon; h++) {
      
      HT          = posterior_omega.col(s) % HT + randn(N);
      sigma2T     = exp(posterior_omega.col(s) % HT); 
      BT_inv      = inv(posterior_B.slice(s));
      mat Sigma_tmp = BT_inv * diagmat(sigma2T) * BT_inv.t();
      SigmaT.slice(h) = 0.5 * (Sigma_tmp + Sigma_tmp.t());
      out_forecast_mean.slice(ss).col(h) = posterior_A.slice(s) * XT;
      vec draw    = mvnrnd( out_forecast_mean.slice(ss).col(h), SigmaT.slice(h) );
      out_forecast.slice(ss).col(h) = draw;
      
      if ( h != horizon - 1 ) {
        XT        = join_cols( draw, XT.subvec(N, K - 1) );
      } // END if h
    } // END h loop
    
    out_forecast_cov(ss) = SigmaT;
    
    ss++;
  } // END s loop
  
  Rcout << "Fraction of used MCMC draws: " << (double) 100 * ss / S << endl;
  
  cube _forecast = out_forecast.slices(0, ss - 1);
  cube _forecast_mean = out_forecast_mean.slices(0, ss - 1);
  field<cube> _forecast_cov = out_forecast_cov.rows(0, ss - 1);
  
  return List::create(
    _["forecast"]       = _forecast,
    _["forecast_mean"]  = _forecast_mean,
    _["forecast_cov"]   = _forecast_cov
  );
} // END forecast_sv

