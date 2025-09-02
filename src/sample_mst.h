
#ifndef _SAMPLE_MST_H_
#define _SAMPLE_MST_H_

#include <RcppArmadillo.h>


arma::mat sample_lambda_ms (
    const arma::mat&    aux_df,     // NxM
    const arma::mat&    aux_xi,     // MxT
    const arma::mat&    U           // NxT
);


double log_kernel_df_ms_nm (
    const double&         aux_df,
    const arma::rowvec&   aux_lambda,  // Tmx1
    const double&         prior_df_a
);


Rcpp::List sample_df_ms (
    arma::mat&        aux_df,             // NxM
    const arma::mat&  aux_lambda,         // NxT
    const arma::mat&  aux_xi,             // MxT
    const arma::mat&  U,                  // NxT
    const double&     prior_df_a,         // hyper-parameter for exponential prior for aux_df
    const int&        s,                  // MCMC iteration
    arma::mat&        adaptive_scale,     // NxM
    const arma::vec&  adptive_alpha_gamma // 2x1 vector with target acceptance rate and step size
);


#endif  // _SAMPLE_MST_H_