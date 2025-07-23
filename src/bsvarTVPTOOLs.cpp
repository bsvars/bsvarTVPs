
#include <RcppArmadillo.h>
#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvars_ir1 (
    arma::mat&    aux_B,              // (N, N)
    arma::mat&    aux_A,              // (N, K)
    const int     horizon,
    const int     p,
    const bool    standardise = false
) {
  
  const int       N = aux_B.n_rows;
  cube            aux_irfs(N, N, horizon + 1);  // + 0 horizons
  mat             A_bold_tmp(N * (p - 1), N * p, fill::eye);
  
  mat   irf_0         = inv(aux_B);
  if ( standardise ) {
    irf_0             = irf_0 * diagmat(pow(diagvec(irf_0), -1));
  }
  mat   A_bold        = join_cols(aux_A.cols(0, N * p - 1), A_bold_tmp);
  mat   A_bold_power  = A_bold;
  
  aux_irfs.slice(0)   = irf_0;
  
  for (int h=1; h<horizon + 1; h++) {
    aux_irfs.slice(h) = A_bold_power.submat(0, 0, N-1, N-1) * irf_0;
    A_bold_power      = A_bold_power * A_bold;
  } // END h loop
  
  return aux_irfs;
} // END bsvars_ir1

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_ir_ms (
    arma::field<arma::cube>&  posterior_B,        // (S)(N, N, M)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p,
    const bool                standardise = false
) {
  
  const int       N = posterior_B(0).n_rows;
  const int       M = posterior_B(0).n_slices;
  const int       S = posterior_B.n_elem;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S, M);
  
  for (int s=0; s<S; s++) {
    for (int m=0; m<M; m++) {
      aux_irfs            = bsvars_ir1( posterior_B(s).slice(m), posterior_A.slice(s), horizon, p , standardise);
      irfs(s, m)          = aux_irfs;
    } // END m loop
  } // END s loop
  
  return irfs;
} // END bsvarTVPs_ir_ms



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_ir_mssa (
    arma::field<arma::cube>&  posterior_B,        // (S)(N, N, M)
    arma::field<arma::cube>&  posterior_A,        // (S)(N, K, S)
    const int                 horizon,
    const int                 p,
    const bool                standardise = false
) {
  
  const int       N = posterior_B(0).n_rows;
  const int       M = posterior_B(0).n_slices;
  const int       S = posterior_B.n_elem;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S, M);
  
  for (int s=0; s<S; s++) {
    for (int m=0; m<M; m++) {
      aux_irfs            = bsvars_ir1( posterior_B(s).slice(m), posterior_A(s).slice(m), horizon, p , standardise);
      irfs(s, m)          = aux_irfs;
    } // END m loop
  } // END s loop
  
  return irfs;
} // END bsvarTVPs_ir_mssa






// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_ir (
    arma::cube&               posterior_B,        // (N, N, S)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p,
    const bool                standardise = false
) {
  
  const int       N = posterior_B.n_rows;
  const int       M = 1;
  const int       S = posterior_B.n_slices;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S, M);
  
  for (int s=0; s<S; s++) {
    aux_irfs            = bsvars_ir1( posterior_B.slice(s), posterior_A.slice(s), horizon, p, standardise);
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




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarTVPs_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::mat&      X                   // KxT
) {
  
  const int   N = posterior_A.n_rows;
  const int   S = posterior_A.n_slices;
  const int   T = X.n_cols;
  
  cube    fitted_values(N, T, S);
  
  for (int s=0; s<S; s++) {
    fitted_values.slice(s) = posterior_A.slice(s) * X;
  } // END s loop
  
  return fitted_values;
} // END bsvarTVPs_fitted_values



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarTVPs_structural_shocks (
    const arma::field<arma::cube>&  posterior_B,    // (S)(N, N, M)
    const arma::cube&         posterior_A,    // (N, K, S)
    const arma::cube&         posterior_xi,   // (M, T, S)
    const arma::mat&          Y,              // NxT dependent variables
    const arma::mat&          X               // KxT dependent variables
) {
  
  const int       N = Y.n_rows;
  const int       T = Y.n_cols;
  const int       S = posterior_A.n_slices;
  
  cube            structural_shocks(N, T, S);
  
  for (int s=0; s<S; s++) {
    for (int t=0; t<T; t++) {
      int m  = index_max(posterior_xi.slice(s).col(t));
      structural_shocks.slice(s).col(t)    = posterior_B(s).slice(m) * (Y.col(t) - posterior_A.slice(s) * X.col(t));
    } // END t loop
  } // END s loop
  
  return structural_shocks;
} // END bsvarTVPs_structural_shocks



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_covariances_rf_mssv (
    const arma::field<arma::cube>&  posterior_B,  // (S)(N, N, M)
    const arma::cube&         posterior_xi,       // (M, T, S)
    const arma::cube&         posterior_sigma     // (N, T, S)
) {
  
  const int       N = posterior_sigma.n_rows;
  const int       T = posterior_sigma.n_cols;
  const int       S = posterior_sigma.n_slices;
  
  field<cube>     cov_rf(S);
  cube            cov_tmp(N, N, T);
  
  for (int s=0; s<S; s++) {
    for (int t=0; t<T; t++) {
      int m             = index_max(posterior_xi.slice(s).col(t));
      mat Bms_inv       = inv(posterior_B(s).slice(m));
      cov_tmp.slice(t)  = Bms_inv * diagmat(square(posterior_sigma.slice(s).col(t))) * Bms_inv.t();
    } // END t loop
    cov_rf(s)         = cov_tmp;
  } // END s loop
  
  return cov_rf;
} // END bsvarTVPs_covariances_rf_mssv



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_covariances_rf_sv (
    const arma::cube&         posterior_B,        // (N, N, S)
    const arma::cube&         posterior_sigma     // (N, T, S)
) {
  
  const int       N = posterior_sigma.n_rows;
  const int       T = posterior_sigma.n_cols;
  const int       S = posterior_sigma.n_slices;
  
  field<cube>     cov_rf(S);
  cube            cov_tmp(N, N, T);
  
  for (int s=0; s<S; s++) {
    for (int t=0; t<T; t++) {
      mat Bs_inv        = inv(posterior_B.slice(s));
      cov_tmp.slice(t)  = Bs_inv * diagmat(square(posterior_sigma.slice(s).col(t))) * Bs_inv.t();
    } // END t loop
    cov_rf(s)         = cov_tmp;
  } // END s loop
  
  return cov_rf;
} // END bsvarTVPs_covariances_rf_sv



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_covariances_rf_ms (
    const arma::field<arma::cube>&  posterior_B,  // (S)(N, N, M)
    const arma::cube&         posterior_xi        // (M, T, S)
) {
  
  const int       N = posterior_B(0).n_rows;
  const int       T = posterior_xi.n_cols;
  const int       S = posterior_xi.n_slices;
  
  field<cube>     cov_rf(S);
  cube            cov_tmp(N, N, T);
  
  for (int s=0; s<S; s++) {
    for (int t=0; t<T; t++) {
      int m             = index_max(posterior_xi.slice(s).col(t));
      mat Bms_inv       = inv(posterior_B(s).slice(m));
      cov_tmp.slice(t)  = Bms_inv * Bms_inv.t();
    } // END t loop
    cov_rf(s)         = cov_tmp;
  } // END s loop
  
  return cov_rf;
} // END bsvarTVPs_covariances_rf_ms


// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::cube bsvarTVPs_covariances_rf (
    const arma::cube&  posterior_B  // (N, N, S)
) {
  
  const int       S = posterior_B.n_slices;
  const int       N = posterior_B.n_rows;
  cube            cov_rf(N, N, S);
  
  for (int s=0; s<S; s++) {
      mat Bms_inv   = inv(posterior_B.slice(s));
      cov_rf.slice(s)  = Bms_inv * Bms_inv.t();
  } // END s loop
  
  return cov_rf;
} // END bsvarTVPs_covariances_rf



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_cov2cor (
    const arma::field<arma::cube>&  posterior_cov  // (S)(N, N, T)
) {
  
  const uword     N = posterior_cov(0).n_rows;
  const uword     T = posterior_cov(0).n_slices;
  const uword     S = posterior_cov.n_elem;
  
  field<cube>     cor_rf(S);
  cube            cov_tmp(N, N, T);
  
  for (int s=0; s<S; s++) {
    for (int t=0; t<T; t++) {
      mat pcov_tmp      = posterior_cov(s).slice(t);
      mat sd_inv        = diagmat( 1 / sqrt(pcov_tmp.diag()) );
      cov_tmp.slice(t)  = sd_inv * pcov_tmp * sd_inv;
      cor_rf(s)         = cov_tmp;
    } // END t loop
  } // END s loop
  
  return cor_rf;
} // END bsvarTVPs_cov2cor



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::cube bsvarTVPs_cov2sd (
    const arma::field<arma::cube>&  posterior_cov  // (S)(N, N, T)
) {
  
  const uword     N = posterior_cov(0).n_rows;
  const uword     T = posterior_cov(0).n_slices;
  const uword     S = posterior_cov.n_elem;
  
  cube            sd_rf(N, T, S);
  
  for (int s=0; s<S; s++) {
    for (int t=0; t<T; t++) {
      mat pcov_tmp      = posterior_cov(s).slice(t);
      sd_rf.slice(s).col(t)  = sqrt(pcov_tmp.diag());
    } // END t loop
  } // END s loop
  
  return sd_rf;
} // END bsvarTVPs_cov2sd




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvars_structural_shocks (
    const arma::cube&     posterior_B,    // (N, N, S)
    const arma::cube&     posterior_A,    // (N, K, S)
    const arma::mat&      Y,              // NxT dependent variables
    const arma::mat&      X               // KxT dependent variables
) {
  
  const int       N = Y.n_rows;
  const int       T = Y.n_cols;
  const int       S = posterior_A.n_slices;
  
  cube            structural_shocks(N, T, S);
  
  for (int s=0; s<S; s++) {
    structural_shocks.slice(s)    = posterior_B.slice(s) * (Y - posterior_A.slice(s) * X);
  } // END s loop
  
  return structural_shocks;
} // END bsvars_structural_shocks




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::rowvec normalisation_wz2003_s (
    const arma::mat& B,                   // NxN
    const arma::mat& B_hat_inv,           // NxN
    const arma::mat& Sigma_inv,           // NxN
    const arma::mat& diag_signs           // KxN
) {
  // returns a rowvec of signs
  const int N         = B.n_rows;
  const int K         = pow(2,N);
  vec       distance(K);
  
  for (int k=0; k<K; k++) {
    mat   B_tmp_inv   = inv(diagmat(diag_signs.row(k)) * B);
    mat   dist_tmp    = trans(B_tmp_inv - B_hat_inv);
    for (int n=0; n<N; n++) {
      distance(k)     += as_scalar(dist_tmp.row(n) * Sigma_inv * trans(dist_tmp.row(n)));
    } // END n loop
  } // END k loop
  
  rowvec out = diag_signs.row(distance.index_min());
  return out;
} // END normalisation_wz2003_s





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvars_normalisation_wz2003 (
    arma::cube        posterior_B,        // NxNxS
    const arma::mat&  B_hat               // NxN
) {
  // changes posterior_B by reference filling it with normalised values
  const int   N       = posterior_B.n_rows;
  const int   K       = pow(2, N);
  const int   S       = posterior_B.n_slices;
  
  mat B_hat_inv       = inv(B_hat);
  mat Sigma_inv       = B_hat.t() * B_hat;
  
  // create matrix diag_signs whose rows contain all sign combinations
  mat diag_signs(K, N);
  vec omo             = as<vec>(NumericVector::create(-1,1));
  for (int n=0; n<N; n++) {
    vec os(pow(2, n), fill::ones);
    vec oos(pow(2, N-1-n), fill::ones);
    diag_signs.col(n) = kron(os, kron(omo, oos));
  }
  
  // normalisation
  cube out(N, N, S);
  for (int s=0; s<S; s++) {
    rowvec sss            = normalisation_wz2003_s(posterior_B.slice(s), B_hat_inv, Sigma_inv, diag_signs);
    mat B_norm            = diagmat(sss) * posterior_B.slice(s);
    out.slice(s)  = B_norm;
  }
  
  return out;
} // END bsvars_normalisation_wz2003


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat bsvars_normalisation_wz20031 (
    arma::mat         aux_B,        // NxNxS
    const arma::mat&  B_hat               // NxN
) {
  // changes posterior_B by reference filling it with normalised values
  const int   N       = aux_B.n_rows;
  const int   K       = pow(2, N);
  
  mat B_hat_inv       = inv(B_hat);
  mat Sigma_inv       = B_hat.t() * B_hat;
  
  // create matrix diag_signs whose rows contain all sign combinations
  mat diag_signs(K, N);
  vec omo             = as<vec>(NumericVector::create(-1,1));
  for (int n=0; n<N; n++) {
    vec os(pow(2, n), fill::ones);
    vec oos(pow(2, N-1-n), fill::ones);
    diag_signs.col(n) = kron(os, kron(omo, oos));
  }
  
  // normalisation
  mat out(N, N);
  rowvec sss            = normalisation_wz2003_s(aux_B, B_hat_inv, Sigma_inv, diag_signs);
  mat B_norm            = diagmat(sss) * aux_B;
  out                   = B_norm;

  return out;
} // END bsvars_normalisation_wz20031



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec bsvars_normalisation_wz20031_diag (
    arma::mat         aux_B,        // NxNxS
    const arma::mat&  B_hat               // NxN
) {
  // changes posterior_B by reference filling it with normalised values
  const int   N       = aux_B.n_rows;
  const int   K       = pow(2, N);
  
  mat B_hat_inv       = inv(B_hat);
  mat Sigma_inv       = B_hat.t() * B_hat;
  
  // create matrix diag_signs whose rows contain all sign combinations
  mat diag_signs(K, N);
  vec omo             = as<vec>(NumericVector::create(-1,1));
  for (int n=0; n<N; n++) {
    vec os(pow(2, n), fill::ones);
    vec oos(pow(2, N-1-n), fill::ones);
    diag_signs.col(n) = kron(os, kron(omo, oos));
  }
  
  // normalisation
  rowvec sss            = normalisation_wz2003_s(aux_B, B_hat_inv, Sigma_inv, diag_signs);
  vec out               = sss.t();
  
  return out;
} // END bsvars_normalisation_wz20031
