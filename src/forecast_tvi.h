#ifndef _FORECAST_TVI_H_
#define _FORECAST_TVI_H_

#include <RcppArmadillo.h>

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
);


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
);


#endif  // _FORECAST_TVI_H_