#ifndef _FORECAST_TVI_H_
#define _FORECAST_TVI_H_


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
);


#endif  // _FORECAST_TVI_H_