
#ifndef _BSVAR_MSS_SV_BOOST_H_
#define _BSVAR_MSS_SV_BOOST_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

Rcpp::List bsvar_mss_sv_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100, // introduce thinning
    const bool                    centred_sv = false,  // introduce thinning
    const int                     hyper_select = 1
);

#endif  // _BSVAR_MSS_SV_BOOST_H_