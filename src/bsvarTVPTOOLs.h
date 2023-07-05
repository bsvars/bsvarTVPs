
#ifndef _BSVARTVPTOOLS_H_
#define _BSVARTVPTOOLS_H_

#include <RcppArmadillo.h>


arma::field<arma::cube> bsvarTVPs_ir_ms (
    arma::field<arma::cube>&  posterior_B,        // (S)(N, N, M)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p,
    const bool                standardise = false
);


arma::field<arma::cube> bsvarTVPs_ir (
    arma::cube&               posterior_B,        // (N, N, S)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p,
    const bool                standardise = false
);


arma::cube bsvarTVPs_filter_forecast_smooth (
    Rcpp::List&       posterior,
    const arma::mat&  Y,
    const arma::mat&  X,
    const bool        forecasted,
    const bool        smoothed
);


arma::cube bsvarTVPs_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::mat&      X                   // KxT
);


arma::field<arma::cube> bsvarTVPs_covariances_rf_mssv (
    const arma::field<arma::cube>&  posterior_B,  // (S)(N, N, M)
    const arma::cube&         posterior_xi,       // (M, T, S)
    const arma::cube&         posterior_sigma     // (N, T, S)
);


arma::field<arma::cube> bsvarTVPs_covariances_rf_sv (
    const arma::cube&         posterior_B,        // (N, N, S)
    const arma::cube&         posterior_sigma     // (N, T, S)
);


arma::field<arma::cube> bsvarTVPs_covariances_rf_ms (
    const arma::field<arma::cube>&  posterior_B,  // (S)(N, N, M)
    const arma::cube&         posterior_xi        // (M, T, S)
);


arma::cube bsvarTVPs_covariances_rf (
    const arma::cube&  posterior_B  // (N, N, S)
);


arma::field<arma::cube> bsvarTVPs_cov2cor (
    const arma::field<arma::cube>&  posterior_cov  // (S)(N, N, T)
);


arma::cube bsvarTVPs_cov2sd (
    const arma::field<arma::cube>&  posterior_cov  // (S)(N, N, T)
);


arma::cube bsvars_structural_shocks (
    const arma::field<arma::cube>&  posterior_B,    // (S)(N, N, M)
    const arma::cube&         posterior_A,    // (N, K, S)
    const arma::cube&         posterior_xi,   // (M, T, S)
    const arma::mat&          Y,              // NxT dependent variables
    const arma::mat&          X               // KxT dependent variables
);


arma::cube bsvars_structural_shocks (
    const arma::cube&     posterior_B,    // (N, N, S)
    const arma::cube&     posterior_A,    // (N, K, S)
    const arma::mat&      Y,              // NxT dependent variables
    const arma::mat&      X               // KxT dependent variables
);


arma::cube  bsvars_normalisation_wz2003 (
    arma::cube        posterior_B,        // NxNxS
    const arma::mat&  B_hat               // NxN
);


arma::mat bsvars_normalisation_wz20031 (
    arma::mat         aux_B,              // NxNxS
    const arma::mat&  B_hat               // NxN
);


#endif  // _BSVARTVPTOOLS_H_